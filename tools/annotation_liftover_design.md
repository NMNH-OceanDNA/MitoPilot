# Annotation liftover: design

Branch: `annotation-liftover` (off `main`, bdeed7a4). Date: 2026-11-17.

## Intent

Annotators miss genes, most often tRNAs and rRNAs, varying by taxon. Today
users fix this by hand-editing the exported `.tbl` or by adding the gene in
Geneious. The reference genome shown in the annotate details modal usually
has the missing gene. Let the user copy it across: pick a gene on the
reference track, map it to the sample assembly, and add it as a sample
annotation when the mapping is good.

## Goals

- In the annotate details modal, pick a gene on the reference track of the
  BLAST Ref Synteny viewer and map it onto the sample assembly unit that is
  open (ID, path, scaffold).
- Support PCG, tRNA, and rRNA reference genes.
- Handle genes that wrap the origin of a circular assembly.
- Show every place the gene maps well (alternates), with quality numbers.
- Only let the user add a mapping that passes strict thresholds.
- Warn and require confirmation when the mapping substantially overlaps an
  existing sample annotation.
- New row gets `tool = "<accession>_liftover"` (e.g. `NC_012345.1_liftover`).
- Run only the new row through curation and validation, in the app's R
  process, not Nextflow.
- No new dependencies. Reuse existing code and patterns.

## Non-goals

- Control region liftover (not offered).
- Split reference gene models (intron or translational frameshift). Not
  selectable for liftover for now.
- Mapping to assembly units other than the one open in the modal.
- Batch liftover of many genes at once.
- Saving threshold settings per project.

## Anti-goals

- A poorly mapped gene ending up in a user's annotations. When in doubt,
  refuse.
- Mapping restricted or biased by where existing sample annotations sit.
  Existing annotations matter only at the overlap check, after mapping.

## Constraints

- `main` lacks unreleased `devel` work on this same modal (flag chips,
  Rotate button, annotation import). Expect merge conflicts in
  `R/app_annotate_details.R` later; keep the footprint there small and put
  new logic in a new file.
- Must work in both the container and the no-container install.
- ASCII only in R code. Minimal comments.

## Approaches

1. **blastn only.** Fast, finds all hits and both strands. Rejected: local
   alignment clips divergent gene ends, so coordinates come out short, worst
   on 70 bp tRNAs. Needs temp files and an external binary.
2. **R aligner only (chosen by user).** `pwalign::pairwiseAlignment` with
   `type = "global-local"`: the whole reference gene is forced end to end,
   free to land anywhere on the assembly. Exact start/stop, pure R, no
   files. Alternates by mask-and-repeat.
3. **blastn to seed, R aligner to refine.** Most robust, most code. Fallback
   only if option 2 proves too slow or misses real hits.

Option 2 was checked on simulated data (17 kb assembly; genes mutated ~12%
with indels): 70 bp tRNAs, 1.6 kb rRNA, 1.55 kb PCG, both strands, two
origin-wrapping cases, and a duplicated tRNA. All recovered at exact
coordinates. Time for 3 rounds x 2 strands: 0.1 s (tRNA) to 2.1 s (1.6 kb).
Absent genes scored 35-60% identity, far below the floor. Script:
scratchpad `liftover_check.R` (throwaway; logic reproduced below).

## Design

### Mapping (new pure function, no Shiny, no DB)

New file `R/annotation_liftover.R`. One function, roughly
`liftover_map(gene_seq, assembly_seq, circular, max_hits = 3)`, returning a
data frame of hits: `pos1`, `pos2`, `direction`, `pctid`, `gap_pct`, `score`.

- Reference gene sequence: slice `blast_ref_sequences.sequence` by the
  `blast_ref_annotations` row (`pos1`, `pos2`, `direction`); reverse
  complement for minus-strand genes so the query is in gene orientation.
  If the reference gene itself wraps the reference origin, splice the two
  pieces.
- Sample sequence: `get_assembly(ID, path, scaffold, con)` (the
  `assemblies` table). `annotations.pos1/pos2` are in that sequence's
  coordinates; the synteny alignment's rotation/strand play no part in
  mapping. Topology comes from the same call.
- Circular sample: subject is the assembly plus its first `nchar(gene)`
  bases appended. Fold hit coordinates back with modulo. A wrapping hit is
  stored `pos1 > pos2`, the existing convention. Length and sequence for
  such rows (including a wrapping PCG's codons and translation) come from
  the existing helpers (`circ_len`, `extract_circ_region`, `wrap_pos`), so
  they are computed on the spliced sequence.
- Scoring: same as `compute_blast_ref_alignment()` in
  `R/blast_ref_utils.R`: `nucleotideSubstitutionMatrix(match = 1,
  mismatch = -1, baseOnly = FALSE)`, `gapOpening = 10`, `gapExtension = 4`.
  Sanitise non-IUPAC characters to N the same way.
- Align the gene and its reverse complement; keep the higher score.
- Alternates: record the best hit, overwrite its bases (and the appended
  twin copy) with N, repeat. Stop at `max_hits` or when a hit fails the
  hard floor.
- `pctid` = matches / alignment length. `gap_pct` = total inserted plus
  deleted bases / reference gene length.
- `direction` of the hit is the strand of the gene on the sample.

### Thresholds

- Defaults: identity >= 85%, gaps <= 5%.
- User-adjustable in the liftover dialog, clamped to hard floors that
  cannot be passed: identity never below 80%, gaps never above 10%.
- Session-only; resets to defaults when the modal reopens (my call).
- Only hits passing the current thresholds can be added. Hits that fail
  the current thresholds but pass the hard floors are listed greyed, with
  their numbers, so the user can see why. Hits below the hard floors are
  not shown. If nothing passes the floors: "No acceptable mapping found."

### User flow

1. User clicks a gene on the reference track. Existing behavior (zoom to
   that gene) stays. The click now also tells the server which reference
   gene was picked (`synteny.js` currently posts to Shiny only for sample
   features).
2. A "Lift over <gene>" button appears near the synteny viewer controls.
   Hidden or disabled with a reason when: annotations are locked, the
   feature is a control region (`type == "ctrl"`), or the gene is a split
   model. Reference features come from NCBI GFF3, one row per CDS segment,
   with no segment/parent ID stored. Rule (my call): a reference gene is
   split when it is a PCG and `blast_ref_annotations` has more than one PCG
   row with that gene name on that accession. tRNA and rRNA rows sharing a
   name (e.g. two trnL or trnS) are separate copies and each is liftable.
   Known cost: a genuinely duplicated PCG with identical names is refused.
   A clicked reference feature is identified by (accession, gene, pos1),
   the table's primary key.
3. Click runs the mapping (spinner; up to a few seconds) and opens a dialog:

```
Lift over rrnS from NC_012345.1 (952 bp)

Min identity [85] %   Max gaps [5] %        (limits: 80%, 10%)

  #  Position        Strand  Identity  Gaps   Overlaps
( ) 1  1,204-2,153     +       93.1%    1.2%   none
( ) 2  9,877-10,820    +       86.0%    4.1%   trnV (12 bp)
    3  14,002-14,950   -       81.5%    7.9%   (below threshold)

                                   [Cancel]  [Add annotation]
```

4. "Add annotation" on a hit with a substantial overlap, or when a live
   sample annotation already has the same gene name, asks for confirmation
   first, naming the clashing gene(s).

### Overlap rule

Compare the chosen hit to every live (non-deleted) sample annotation on the
unit, strand-blind, wrap-aware. Substantial = shared bases > 25% of the
length of either feature. Smaller overlaps are listed in the table but need
no confirmation (tRNAs legitimately overlap neighbors by a few bases).

### The new row

- `ID`, `path`, `scaffold` of the open unit; `gene`, `type`, `product` from
  the reference annotation (builder: reference gene names are already
  normalised for the synteny plot; derive `product` the way existing rows
  of that gene/type get it); `pos1`, `pos2`, `direction`, `length` from the
  hit; `tool = "<accession>_liftover"`; `edited = 1`; `time_stamp` now.
- tRNA `anticodon`: left NA. The reference table stores none, and "NNN"
  would raise a false low-confidence warning.
- PCG `start_codon`, `stop_codon`, `translation`: computed from the sample
  sequence at the mapped coordinates with the sample's genetic code, using
  existing helpers.
- Key clash: the primary key is (ID, path, scaffold, gene, pos1). If a
  live row already has the same gene and pos1, refuse ("already annotated
  at this position"). Deleted rows carry `_DELETED_` in the gene name, so
  they do not clash.
- `notes`: "Liftover from <accession>: <pctid>% identity, <gap_pct>% gaps".

### Curate and validate the one row

Run curation and then validation with the unit's `curate_opts` params
(same rules WF2 uses) on a one-row annotations CSV in a temp dir, with the
unit's assembly FASTA and coverage file. Read the result back and only then
write the row to the DB and `rv$annotations` (nothing is inserted before
curation succeeds).

- Curation must never change the assembly here. `curate_mito_core()` can
  trim linear contigs to their annotations (`feature_trim`), and
  reverse-complement a contig when all its rRNAs are minus-strand
  (`flip_rRNA_minus_strand`, hardcoded TRUE in the clade wrappers) or by
  reference (`ref_based_rc`). On a one-row table these would wreck the
  unit. So (my call): call `curate_mito_core()` directly with all three
  set FALSE, ignoring the unit's own settings for them, and discard the
  assembly and coverage files it writes. As a guard, compare the returned
  assembly to the input and abort the liftover if it differs. Builder:
  check whether any clade wrapper does more than pass arguments through;
  if one does, keep that behavior.

- Curation needs the curation reference DBs and `blastp`/`blastn`; the
  modal already runs local BLAST against them (`get_top_hits_local`,
  `rv$local_db`). If they are not reachable, add the row uncurated with a
  warning "Liftover: curation skipped" rather than failing.
- Whole-genome checks inside curate/validate (missing genes, duplicates,
  gene counts) must not fire off a one-row table. Builder: find which
  checks are table-wide and keep their output from being written for this
  run; prefer filtering results over adding a flag to the core functions,
  but a small `single_feature` argument is acceptable if cleaner.
- Curation may move the row's ends or strand (PCG start/stop repair, rRNA
  trimming). That is expected. If curation deletes or empties the row,
  treat as failure: add nothing, tell the user.
- After the write, refresh the sample's summary counts the same way an
  existing edit/delete does (missing/extra genes, warnings).
- Undo is the existing delete-annotation control.

## Left to the builder

- Exact placement and styling of the button and dialog, following the
  modal's existing conventions.
- How the ref-gene click reaches the server (new input vs extending the
  existing one).
- Function names and internal structure of `R/annotation_liftover.R`.
- Known risk to measure, not solve up front: the 80% identity floor may
  refuse most tRNAs when the reference is a distant relative. The user set
  the floor deliberately (refuse when in doubt). The real-data check below
  should report how many reference tRNAs clear 85% and 80% on the fish and
  jellyfish test projects so the floor can be revisited with numbers.
- Tests: unit tests for the mapping function (forward, reverse, wrap,
  duplicate, absent, N-containing assembly) and the overlap rule; one
  real-data check on a shipped fish test sample (delete a tRNA and an rRNA,
  lift each back from the reference, compare coordinates).
- NEWS.md entry under the development heading.

## Addendum 2026-11-17: dialog views (user request)

The liftover dialog gains two views below the placement table. Dialog becomes
a large modal. One placement is selected at a time (radio; unchanged), and
one liftover adds exactly one placement.

### Alignment view

- Two rows only: the reference gene and the sample sequence at the selected
  placement. Changes when the user selects another placement.
- Same widget, options, and box styling as the Alignment view in the details
  modal (`msaR`, `.mp-msa-box`, same option set as `output$msa`). No edit
  controls.
- tRNA and rRNA: nucleotide alignment. Reuse the mapping alignment itself
  (the aligned gene and aligned sample strings from `pwalign`), so what is
  shown is exactly what was scored.
- PCG: protein alignment when the placement's first codon is a valid start
  codon for that gene (the gene's `start_codons` from the unit's curate
  params, `rules[[gene]]` over `default_rules$PCG`). Reference translated
  with the reference's genetic code (`blast_ref_sequences.genetic_code`),
  sample with the sample's. Otherwise show the nucleotide alignment with a
  visible warning line: "Placement does not start with a valid start codon
  (<codon>), showing nucleotide alignment."
- PCG extra line when the sample translation has internal stop codons:
  "<n> internal stop codon(s) in the mapped region." (my call; cheap and
  directly relevant to whether to add).
- Row labels: "<accession> <gene>" and "<sample ID> <pos1>-<pos2> (<strand>)".
- A greyed (below-threshold) placement cannot be selected, so it has no
  alignment view.

### Synteny view

- A second instance of the BLAST Ref Synteny canvas viewer, same payload as
  the one in the details modal (same alignment, same reference), plus:
  - every placement in the table (passing or greyed) drawn on the sample
    track as an extra feature named after the gene, with a bright highlight
    border;
  - the picked reference gene drawn with the same highlight border.
- Selecting a placement in the table zooms this viewer to that placement.
  The currently selected placement is visually distinct from the other
  highlighted placements (existing "selected" stroke on top of highlight).
- Clicking features in this viewer does nothing server-side (no row
  selection, no ref pick); pan/zoom controls work.
- Highlight colour: one constant, readable in the app's light theme on all
  feature type colours (user's call: bright green border, no amber and no
  blue; the selected placement gets a thicker green border).

## Addendum 2026-11-17 (2): curation of lifted PCGs and rRNAs through Nextflow

Problem: the curation reference database and `blastp`/`blastn` live inside the
container (or the native install), not on the machine running the app, so
in-app curation of a lifted PCG is skipped and a lifted rRNA gets no
reference hits.

Decisions (user):

- tRNA: curate and validate in the app's R process, as built. Instant.
- PCG and rRNA: curate and validate through a small dedicated Nextflow
  workflow, so the project's own configuration decides how tools are reached:
  Docker, Singularity, Apptainer, or a native (no-container) install.
- The workflow's processes always run on the machine where R is running
  (local executor), never submitted to a scheduler, whatever the project
  config's executor is.
- Still not a pipeline run: no Run Report, no entry in the project's run
  bookkeeping, no interaction with a pipeline run's work directory or
  Nextflow history.

Behavior:

- The one-row inputs are the same as the in-app path (pad row included);
  curation runs with `feature_trim`, `flip_rRNA_minus_strand`, and
  `ref_based_rc` all FALSE; the assembly-unchanged guard and all row
  post-processing stay in the app. The workflow only runs the two core
  functions (`curate_mito_core`, `validate_mito_core`) that already exist in
  released images, since the image's installed MitoPilot will not contain
  the liftover code.
- Curation reference database: the app makes it available as a local
  directory and reuses it across lifts. If the project's `ref_dir`/`ref_db`
  already point at a local directory, use it; if they point at a tarball or
  URL, download/unpack once into a cache folder inside the project and reuse.
- While the workflow runs the user sees a waiting state with the text
  "Curating and validating the new annotation, hold tight..." and a Cancel
  control. Cancel stops the workflow and adds nothing.
- If the workflow fails (engine unavailable, image pull fails, tool error,
  database download fails), tell the user why in plain words and offer to
  add the gene uncurated with the existing
  "Liftover: curation skipped (<reason>)" warning. Declining adds nothing.
- On success the row carries what a pipeline run would give that gene:
  curated coordinates/codons, reference hits, validation warnings.
- Scratch files for a lift are removed after success; kept after a failure
  and named in the error text.

## Addendum 2026-11-17 (3): validate the lifted gene against its neighbors

Decision (user): the lifted gene must get the neighbor-based validation
warnings a pipeline run would give it (overlaps, tRNA within a PCG or rRNA,
duplicates), instead of being validated alone.

- Validation of the new row runs against the unit's full live annotation
  table (existing non-deleted rows plus the new row), for both the in-app
  (tRNA) and the Nextflow (PCG, rRNA) paths. Curation still runs on the new
  row alone.
- The new row takes all warnings validation gives it.
- Existing rows are changed in one way only: a warning that validation
  raises on an existing row with the new gene present, and does not raise
  without it, is appended to that row. Nothing is removed from existing
  rows, and no other field of an existing row changes.
- Whole-table summaries (missing and extra genes, counts) keep coming from
  the modal's existing count refresh.
- Unchanged: the overlap confirm (more than a quarter of either gene), the
  same-position refusal, and the rule that uncurated rows do not get
  reference-comparison warnings.
