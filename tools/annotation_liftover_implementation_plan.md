# Annotation Liftover Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** In the annotate details modal, lift a gene from the reference track of the synteny viewer onto the sample assembly as a new, locally curated and validated annotation.

**Architecture:** Pure, Shiny-free logic in a new `R/annotation_liftover.R` (mapping, overlap, row build, one-row curate/validate). `synteny.js` gains one extra Shiny input for reference-gene clicks. `R/app_annotate_details.R` gets a thin wiring block (button, dialog, confirm, insert) that calls the pure functions.

**Tech Stack:** R, Shiny, pwalign/Biostrings (already imported), testthat, existing JS test harness under `tests/testthat/synteny`.

**Spec:** `tools/annotation_liftover_design.md` (read it first; it is authoritative).

## Global Constraints

- No new dependencies (DESCRIPTION Imports unchanged).
- ASCII only in R code. Minimal comments. No em dashes anywhere.
- Do NOT commit, push, or switch branches. Leave changes uncommitted on `annotation-liftover`.
- Keep the diff in `R/app_annotate_details.R` small; logic belongs in `R/annotation_liftover.R`.
- Thresholds: defaults identity >= 85, gaps <= 5; hard floors identity >= 80, gaps <= 10.
- `tool` value: `paste0(accession, "_liftover")`.
- Alignment scoring: `pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)`, `gapOpening = 10`, `gapExtension = 4`, `type = "global-local"`, gene as pattern.
- Curation must never alter the assembly: `feature_trim = FALSE`, `flip_rRNA_minus_strand = FALSE`, `ref_based_rc = FALSE`.
- Qualify `BiocGenerics::start/end` (bare `start`/`end` resolve to stats:: inside the package).
- Run tests with `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/<file>")'`.

## Review Focus

1. Assembly contains N or IUPAC codes (joined Path 0 assemblies): mapping must not error. Test in Task 1.
2. Reference gene longer than a linear sample scaffold: mapping returns zero hits, no error. Test in Task 1.
3. Reference gene itself wraps the reference origin (`pos1 > pos2`): sequence is spliced. Test in Task 1.
4. Curation reference databases or `blastp` not reachable from the app: row is still returned, uncurated, with warning "Liftover: curation skipped". Test in Task 2.
5. Hit touching an existing wrap-around annotation (`pos1 > pos2`): overlap computed around the circle. Test in Task 1.

---

### Task 1: Mapping and overlap (pure functions)

**Files:**
- Create: `R/annotation_liftover.R`
- Test: `tests/testthat/test-annotation-liftover.R`

**Interfaces (Produces):**

```r
LIFTOVER_DEFAULTS <- list(min_pctid = 85, max_gap_pct = 5)
LIFTOVER_FLOORS   <- list(min_pctid = 80, max_gap_pct = 10)

# Gene-orientation sequence of a reference feature. Splices when pos1 > pos2.
liftover_ref_gene_seq(ref_seq, pos1, pos2, direction)  # -> character(1)

# TRUE when the feature must not be lifted: ctrl, or PCG with >1 PCG row of that gene.
# blast_ref = rows of blast_ref_annotations for ONE accession.
liftover_blocked_reason(blast_ref, gene, pos1)  # -> NA_character_ or reason string

# Hits best-first. Only hits passing LIFTOVER_FLOORS are returned (0-row df if none).
# Columns: pos1, pos2 (pos1 > pos2 when wrapping), direction ("+"/"-"),
#          pctid, gap_pct, score.
liftover_map(gene_seq, assembly_seq, circular, max_hits = 3)  # -> data.frame

# Shared bases between a hit and each live annotation, strand-blind, wrap-aware.
# annotations: data.frame with gene, pos1, pos2 (rows with "_DELETED_" in gene ignored).
# Returns data.frame(gene, bp, substantial) for bp > 0; substantial = bp > 25% of
# either feature's length.
liftover_overlaps(pos1, pos2, annotations, asm_len)  # -> data.frame

# Clamp user thresholds to the floors.
liftover_clamp(min_pctid, max_gap_pct)  # -> list(min_pctid, max_gap_pct)
```

Algorithm for `liftover_map` (the tested throwaway version is in the spec's "Approaches" section and reproduced here in outline): sanitise non-IUPAC to N; if circular, subject = assembly + first `nchar(gene)` bases; if linear and gene longer than assembly, return empty; loop up to `max_hits`: align gene and its reverse complement, keep higher score; `pctid = 100 * nmatch / alignment length`; `gap_pct = 100 * (inserted + deleted bases) / nchar(gene)`; fold coordinates with `((x - 1) %% L) + 1`; stop when the hit fails the floors; otherwise record, overwrite the hit span and its twin copy (position +/- L within the extended subject) with N, repeat.

- [ ] **Step 1: Write failing tests.** Build a seeded random 5 kb assembly and a `mutate()` helper (substitutions + indels) inside the test file. Tests:
  - `liftover_map` recovers exact coordinates for: 70 bp forward, 70 bp reverse, 900 bp forward, a 70 bp gene wrapping the origin (expect `pos1 > pos2`), and both copies of a duplicated 70 bp gene (2 rows, best first).
  - Absent gene (fresh random 70 bp and 900 bp): 0 rows.
  - Linear assembly, gene at the very end: found; gene longer than assembly: 0 rows, no error.
  - Assembly with a 50 bp run of N and scattered R/Y codes: no error, true hit still found.
  - All returned rows satisfy `pctid >= 80` and `gap_pct <= 10`.
  - `liftover_ref_gene_seq`: plus strand, minus strand (reverse complement), and `pos1 > pos2` splice.
  - `liftover_blocked_reason`: ctrl blocked; PCG with two rows blocked; tRNA with two same-name rows NOT blocked; single PCG not blocked.
  - `liftover_overlaps`: no overlap -> 0 rows; 5 bp overlap of a 70 bp gene -> not substantial; 30 bp of 70 -> substantial; overlap with a wrap-around annotation (pos1 > pos2) counted around the circle; `_DELETED_` rows ignored.
  - `liftover_clamp(50, 40)` -> `list(min_pctid = 80, max_gap_pct = 10)`; `liftover_clamp(95, 2)` unchanged.
- [ ] **Step 2: Run, verify they fail** (functions not found).
- [ ] **Step 3: Implement** in `R/annotation_liftover.R`. Reuse `circ_len`, `wrap_pos`, `extract_circ_region` from `R/annotate_utils.R` where they fit.
- [ ] **Step 4: Run, verify pass.**

### Task 2: Build, curate, and validate one row (pure, no Shiny)

**Files:**
- Modify: `R/annotation_liftover.R`
- Test: `tests/testthat/test-annotation-liftover-curate.R`

**Interfaces:**
- Consumes: a one-row hit from `liftover_map()`.
- Produces:

```r
# Builds the new annotation row, runs curate then validate on it alone, returns
# list(row = <1-row data.frame with the annotations table columns>, curated = TRUE/FALSE)
# or stops with a user-readable message on failure.
liftover_annotate(
  hit,            # 1-row data.frame from liftover_map()
  ref_feature,    # 1-row: gene, type, pos1, pos2, direction (blast_ref_annotations)
  accession,
  unit,           # list(ID, path, scaffold)
  assembly,       # DNAStringSet from get_assembly(ID, path, scaffold, con), length 1
  coverage_fn,    # path or NULL
  curate_params,  # parsed curate_opts params (list with rules, default_rules, ref_dbs, ...)
  genetic_code,
  ref_dir,
  template        # 0+ row data.frame of the unit's existing annotations (column set / product lookup)
)
```

Requirements (spec sections "The new row" and "Curate and validate the one row"):
- Row fields per spec: `tool = paste0(accession, "_liftover")`, `edited = 1`, `anticodon = NA`, `notes` starts with `"Liftover from <accession>: <pctid>% identity, <gap_pct>% gaps"`, `length` via circular helpers; PCG `start_codon`/`stop_codon`/`translation` from the sample sequence with `genetic_code` (wrap-aware; find and reuse the helper curation/the modal already uses for this).
- `product`: copy from a `template` row with the same gene if present, else derive the way the annotation tools do for that type (find the existing lookup; do not invent a new table).
- Write the row to a temp dir in the CSV shape `curate_mito_core()` reads (study what WF2 passes: the `contig` column and FASTA header format; `get_assembly()` names are `"ID.path.scaffold topology"`). Call `curate_mito_core()` directly with the three flags FALSE, then `validate_mito_core()` on the curated row. First check whether the clade wrappers (`curate_<clade>_mito`) do anything beyond passing arguments; if any does, preserve that behavior and note it in your report.
- Guard: read back the assembly FASTA curation wrote; if it differs from the input, stop("Liftover aborted: curation changed the assembly").
- Curation returns zero rows or a row without coordinates: stop with a clear message.
- If curation errors because reference DBs or BLAST binaries are unreachable: return the uncurated row, `curated = FALSE`, with `warnings` containing `"Liftover: curation skipped"`, and still run validation.
- Keep only the row's own columns from validate's output; discard its summary and the table-wide `missing` list. Strip a "possible duplicate"-style table-wide artifact only if a one-row run can produce a false one (check; a one-row run should not).
- Returned row has exactly the columns of the `annotations` table (see `R/init_db.R:757`).

- [ ] **Step 1: Write failing tests.** Model fixtures on `tests/testthat/test-curate-circular-wraparound.R` (it already drives `curate_mito_core` from files). Tests:
  - tRNA hit on a circular assembly: row has the hit's coordinates, `tool == "NC_TEST.1_liftover"`, `is.na(anticodon)`, `edited == 1`, notes prefix as above.
  - Minus-strand rRNA hit: assembly passed in is unchanged afterward and returned `direction == "-"` (proves the flip is off).
  - Linear assembly, hit in the middle: returned coordinates equal the hit's (proves trimming is off).
  - PCG hit with `ref_dir` pointing at a nonexistent directory: returns a row, and either `curated == FALSE` with the skip warning or a curated row; never an error. `translation`, `start_codon`, `stop_codon` are filled.
  - Wrapping tRNA hit (`pos1 > pos2`): `length` equals the true gene length.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run both liftover test files plus `test-curate-circular-wraparound.R`; all pass.**

### Task 3: Reference-gene click reaches the server

**Files:**
- Modify: `inst/app/www/synteny.js` (click handler near line 335; constructor near 115)
- Test: extend `tests/testthat/test-synteny-js.R` / `tests/testthat/synteny/`

**Interfaces (Produces):** when payload field `refInput` is set, clicking a reference feature still calls `fitFeature` AND calls `Shiny.setInputValue(refInput, {gene, type, pos1, pos2, dir, nonce: Date.now()}, {priority: 'event'})`. `pos1/pos2/dir` are the feature's stored reference coordinates as received in `refFeatures` (not alignment columns). Sample-feature clicks unchanged.

- [ ] **Step 1: Failing JS test** following the existing harness: ref click posts the object above to `refInput`; sample click does not post to `refInput`; no `refInput` in payload -> no post, no error.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement** (store `p.refInput`; a few lines in `click`).
- [ ] **Step 4: Run `test-synteny-js.R` and `test-synteny-payload.R`; pass.**

### Task 4: Modal wiring

**Files:**
- Modify: `R/app_annotate_details.R` (synteny observer near line 1996; new block after it)
- Modify: `R/app_annotate_synteny.R` (`synteny_viewer_ui`: add the button slot and update the caption to mention clicking a reference gene)
- Modify: `NEWS.md` (entry under `# MitoPilot (development version)`)
- Test: `tests/testthat/test-annotation-liftover-ui.R` (only for any pure helper you extract, e.g. the hit-table builder)

**Interfaces (Consumes):** everything from Tasks 1-3; existing `get_assembly()`, `lock_blocked()`, `persist_annotations()`, `mp_confirm()`, `mp_toast()`, `mp_alert()`, `show_edit_waiter()`, `resolve_unit_coverage_file()`, `rv$blast_ref`, `rv$blast_ref_seq`, `active_ref_acc()`, `rv$updating`, `rv$annotations`.

Behavior (spec section "User flow"):
- Set `p$refInput <- ns("synteny_ref_pick")` beside `p$input`.
- On `input$synteny_ref_pick`: remember the picked feature; show an `actionButton` "Lift over <gene>" in the synteny controls row. Use `liftover_blocked_reason()`; when blocked (or `lock_blocked()`), show the button disabled with the reason as its title. Clear the pick when the reference accession changes or the modal opens another unit.
- On button click: waiter on; `liftover_map()` against `get_assembly(ID, path, scaffold, con)` (circular from its name/topology); waiter off; open a modal dialog per the spec sketch: two `numericInput`s (defaults 85 and 5, clamped through `liftover_clamp()`, labelled with the limits), radio list/table of hits with Position, Strand, Identity, Gaps, Overlaps (from `liftover_overlaps()`); hits failing the current thresholds are shown greyed and not selectable; zero hits -> text "No acceptable mapping found." and no Add button.
- "Add annotation": refuse if a live row already has the same gene and pos1 ("already annotated at this position"). If any overlap is substantial OR a live row has the same gene name, `mp_confirm()` naming the genes first. Then waiter on, `liftover_annotate()`, append the row to `rv$annotations` (it must carry the extra non-DB columns `rv$annotations` has, e.g. `faa`/`fas`; see how other code adds rows near lines 1771 and 3590), `persist_annotations()`, close dialog, select the new row, `mp_toast()` success (mention "curation skipped" when `curated` is FALSE). Errors -> `mp_alert()` with the message; nothing written.
- The existing `observeEvent(rv$annotations, ...)` already refreshes counts/missing/extra; verify it fires.
- Curate inputs: `curate_opts` params for the unit the same way `update_counts()` reads them (around line 3704); `genetic_code` and `ref_dir` from wherever WF2 gets them for this unit (trace `curate_opts` / `annotate_opts` columns); coverage via `resolve_unit_coverage_file()`.

- [ ] **Step 1:** Implement UI slot + JS payload field + pick observer + button.
- [ ] **Step 2:** Implement dialog and threshold clamping.
- [ ] **Step 3:** Implement add path with confirm, insert, persist.
- [ ] **Step 4:** `devtools::load_all()` clean; run full liftover + synteny test files; run `devtools::document()` only if roxygen exports changed (none expected; all new functions `@noRd`).
- [ ] **Step 5:** NEWS.md entry.

### Task 6: Dialog alignment and synteny views

Spec: "Addendum 2026-11-17: dialog views" in the design doc (binding; read it).

**Files:**
- Modify: `R/annotation_liftover.R` (hits carry alignment strings; new pure helper)
- Modify: `inst/app/www/synteny.js` (highlight flag)
- Modify: `R/app_annotate_synteny.R` (viewer UI reusable under a second id prefix; payload accepts extras)
- Modify: `R/app_annotate_details.R` (dialog only: the liftover block)
- Test: `tests/testthat/test-annotation-liftover.R`, `test-annotation-liftover-ui.R`, `test-synteny-js.R`, `test-synteny-payload.R`

**Interfaces (Produces):**

```r
# liftover_map() gains two character columns: aln_gene, aln_sample (gapped,
# equal length, gene orientation). Existing columns unchanged.

# Two-row alignment for the dialog. Pure.
# Returns list(aln = <XStringSet, 2 named rows>, kind = "nt" | "aa",
#              notes = character())  # warning lines to show above the view
liftover_alignment(hit, ref_feature, ref_seq, assembly_seq, type,
                   start_codons, sample_code, ref_code, labels)

# synteny_payload() gains optional args, default NULL, existing output
# unchanged when NULL:
#   extra_sample = data.frame(row, gene, type, pos1, pos2, direction)  # hl = TRUE
#   hl_ref = list(gene, pos1)   # marks that refFeature hl = TRUE
# synteny_viewer_ui(ns, prefix = "synteny")  # ids: <prefix>_canvas, <prefix>_tip
```

JS: features (sample or ref) with `hl: true` get a highlight border (thicker, constant colour) drawn under the existing selected stroke. Nothing else changes; payloads without `hl` render exactly as before.

Dialog wiring: `modalDialog(size = "l")`; below the table an `msaR::msaROutput` in a `.mp-msa-box` (same `msaR::msaR()` options as `output$msa` near line 1683, minus the editing bits) fed by `liftover_alignment()` for the selected placement, with `notes` rendered above it; below that `synteny_viewer_ui(ns, "liftover_syn")` fed by `synteny_payload(..., extra_sample = <all table placements, row = paste0("lift", i)>, hl_ref = <picked gene>)` with no `input`/`refInput` in the message; on radio change send `mpsyn_select` with `row = paste0("lift", i)` to zoom. Reference genetic code from `rv$blast_ref_seq$genetic_code`; start codons from the unit's curate params.

- [ ] **Step 1: Failing tests.**
  - `liftover_map`: `aln_gene`/`aln_sample` present, equal `nchar`, and removing "-" from `aln_gene` gives the gene sequence; all earlier tests still pass.
  - `liftover_alignment`: tRNA -> `kind == "nt"`, 2 rows, names equal `labels`; PCG placement starting with ATG and `start_codons = c("ATG")` -> `kind == "aa"`, no notes; PCG starting with a codon not in `start_codons` -> `kind == "nt"` and a note matching "does not start with a valid start codon \\(...\\)"; PCG with an internal stop -> note matching "1 internal stop codon"; wrapping placement (pos1 > pos2) works; minus-strand placement works.
  - `synteny_payload`: with `extra_sample` of 2 rows, `sampleFeatures` gains 2 features with `hl` TRUE and the given `row` values; `hl_ref` marks exactly one ref feature; with both NULL the payload equals the old payload (existing tests untouched).
  - JS: a feature with `hl: true` is drawn with the highlight stroke (assert via the harness's canvas inspection or an exposed state hook, following existing test style); `mpsyn.fit(id, "lift1")` zooms to that feature.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run all five liftover/synteny test files; `load_all` clean.**

### Task 7: Curate lifted PCGs and rRNAs through Nextflow

Spec: "Addendum 2026-11-17 (2)" in the design doc (binding; read it).

**Files:**
- Create: `inst/nextflow/liftover_curate.nf` (standalone entry script; do not touch `main.nf`)
- Modify: `R/annotation_liftover.R` (split `liftover_annotate()` into prepare / run / finish; Nextflow runner; DB cache)
- Modify: `R/app_annotate_details.R` (liftover add path only: async run, waiting state with Cancel, failure fallback)
- Modify: `NEWS.md`, `vignettes/Test-Project-Annotate.Rmd` (the curation caveat bullet), tests
- Test: `tests/testthat/test-annotation-liftover-curate.R` plus a new `test-annotation-liftover-nf.R`

**Interfaces (Produces):**

```r
# Existing signature and return value unchanged; used for tRNA and by tests.
liftover_annotate(...)

# Split used by both paths (names fixed; arguments are the builder's call):
liftover_prepare(...)   # builds the row, writes temp inputs (pad row), returns a state list
liftover_finish(state, curated_dir_or_NULL, skip_reason = NULL)  # guard, strip pad, filter warnings -> list(row, curated)

# TRUE when this gene type goes through Nextflow.
liftover_needs_nf(type)                 # "PCG", "rRNA" -> TRUE; "tRNA" -> FALSE

# Local directory holding the unpacked curation DB for the unit, cached under
# <project>/.liftover/ref_dbs/. Downloads/untars once. Errors with a plain message.
liftover_ref_db_dir(ref_dir, ref_db, project_dir)

# Command + env for the run; pure, unit-testable.
# Returns list(command, args, wd, env).
liftover_nf_command(state, db_dir, project_dir, config = file.path(project_dir, ".config"))
```

Fixed requirements:
- Launch directory is a per-lift scratch folder under `<project>/.liftover/` (so Nextflow's `.nextflow` history, logs, and work dir never mix with the project's); config passed with `-c <project>/.config`; executor forced local for every process (`-process.executor local` or an equivalent that is proven to win over the project config); no `-resume`.
- Nextflow launcher, pinned version, and native environment are resolved exactly as the app's existing run path does (`R/app_run_pipline.R` around lines 300-345: `native_env()`, `params.native_activate`, the NXF_VER pin). Reuse those helpers; do not re-implement.
- The process script calls only `MitoPilot::curate_mito_core()` and `MitoPilot::validate_mito_core()`. Confirm in the project's image that both accept the arguments used (`docker run --rm macguigand/mitopilot:1.5.10 Rscript -e 'args(MitoPilot::curate_mito_core)'`); if an argument is missing in the released image, stop and report rather than work around it.
- The project config sets process-wide closures (`cpus = { opts.cpus }`, `memory = { (opts?.memory?.GB ?: 1) * task.attempt }`) and loads the `nf-sqldb` plugin; the workflow must run under that config unmodified. Prove it by running, not by reading.
- No Run Report, nothing written under `.runs/`, no writes to the project `.sqlite` from the workflow (the app writes the row as today).
- Waiting state: text exactly "Curating and validating the new annotation, hold tight...", with a Cancel control that kills the process tree and adds nothing. The Shiny session must stay responsive (async process + polling), not block.
- Failure: plain-language reason, then a confirm offering to add the gene uncurated (existing "Liftover: curation skipped (<reason>)" warning); declining adds nothing.
- tRNA path unchanged and still instant.

- [ ] **Step 1: Failing unit tests.** `liftover_needs_nf`; `liftover_prepare` + `liftover_finish` reproduce `liftover_annotate()` results for a tRNA; `liftover_finish` with a skip reason gives the uncurated row and warning; `liftover_ref_db_dir` with a local directory returns it untouched, with a local tarball unpacks once and reuses (second call does not re-extract), with an unreachable URL errors plainly; `liftover_nf_command` contains the config path, the local-executor override, the scratch launch dir, no `-resume`, and native env vars when the config has `params.native_activate`.
- [ ] **Step 2: Run, verify fail. Step 3: Implement. Step 4: all liftover and synteny test files pass; `load_all` clean.**
- [ ] **Step 5: Live proof on this machine (Docker).** On a fresh scratch copy of `~/MitoPilot_scratch/docs_ann`: delete a PCG and an rRNA, lift each back through the real UI, and show: the waiting state and Cancel; the row added with NO "curation skipped" warning; refHits populated; coordinates and start/stop codons equal to the original pipeline annotation; wall time per lift (first lift including DB download, and a later lift); assembly sha256 unchanged; `.runs/` untouched; scratch folder removed. Also: Cancel mid-run adds nothing; a forced failure (e.g. nonexistent container image in a copy of the config) produces the fallback confirm and an uncurated add.
- [ ] **Step 6: Other environments.** Singularity is installed here (`/usr/local/bin/singularity`): run the workflow once from the command line with a Singularity variant of the config and report the result. Native: if a native install exists on this machine (see `dev/no_container/`), run once with its config; otherwise report NOT RUN and show, from `liftover_nf_command` unit tests, that the native env is applied.

### Task 5: Real-data check and app test

**Files:**
- Create: `dev/annotation_liftover/real_data_check.R` (gitignored area) and `dev/annotation_liftover/REPORT.md`

- [ ] **Step 1:** On a scratch copy of a completed fish test project and the Scyphozoa test project (never the originals): for every liftable gene on each unit's reference, run `liftover_map()` and record best-hit pctid/gap_pct, whether it clears 85/5 and 80/10, and whether its coordinates match the sample's existing annotation of the same gene (within a few bp). Report pass rates by type (tRNA, rRNA, PCG).
- [ ] **Step 2:** In the running app on the scratch project: delete one tRNA and one rRNA, lift each back from the reference, confirm coordinates, `tool` value, counts refresh, and that the assembly sequence in the DB is byte-identical before and after. Try an overlapping lift (confirm appears) and a split/ctrl gene (blocked).
- [ ] **Step 3:** Write findings to `REPORT.md`, including the tRNA pass-rate numbers the spec asks for.
