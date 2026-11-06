# Original user-annotation track: design

Date: 2026-11-06. Branch: `user-annot-import` (continues the unmerged
user-annotation import feature; see
`2026-11-02-user-annotation-import-design.md`).

## Intent

A user who imported their own annotation wants to see what MitoPilot's
curation did to it. Today the Annotation Details window shows only the
curated result. This feature adds a second gene track, the user's original
annotation exactly as imported, drawn above the curated track in both the
sequence viewer and the BLAST reference synteny viewer. Features that
curation changed stand out, so the user can spot boundary shifts and
dropped genes at a glance. The track uses the user's own raw gene names, so
mistakes in the gene-name table (a feature called `COI` that was mapped to
the wrong gene, say) are also visible.

## Goals

- Draw every imported feature of the sample as a muted track above the
  curated lanes, in both viewers, for user-annotation projects only.
- Label each original feature with the raw name from the user's file.
- Highlight original features that curation edited: boundary or strand
  changed, or feature absent from the curated table.
- Place original features in the CURRENT coordinate frame even after WF2
  rotation, app Linearize, or app Trim moved the sequence.
- A checkbox in each viewer hides and shows the track. Default on.
- Tooltip on hover explains the status of each original feature.

## Non-goals

- Clicking an original feature does not select a table row and does not
  edit anything. The track is read-only context.
- Name normalization is NOT a highlighted edit. `COI` becoming `COX1` is
  expected for every gene and would drown the real signal (user decision).
- Features curation ADDED (present only in the curated table) get no
  marker. They are visible by comparing tracks.
- No new track in the outlier review, MSA, circular/linear overview, or
  export views.
- No change to standard or userAsmb project flows; the track is absent there.
- Standalone diff table or report of edits. Viewer only.

## Anti-goals

- Drawing original features at wrong positions after the sequence frame
  moved. Wrong positions are worse than no track: when the frame cannot be
  recovered, hide the track and say why.
- Bookkeeping of rotation offsets scattered through WF2 and the asmb-edit
  code. Frame recovery is done by sequence search at view time, from stored
  sequences, not from recorded offsets.
- Visual clutter: the original track must read as background. Curated
  lanes keep their current look.

## User decisions (brainstorm 2026-11-06)

- Show ALL original features always (not only when something changed),
  highlight the edited ones.
- Use raw user names on the track.
- Boundary shifts are the main thing to see; dropped or added genes second.
- Ignore pure name changes for the highlight.
- Account for frame drift (WF2 rotation, Linearize, Trim) wherever the
  mapping can be recovered. Accepted limitation: a Trim that cuts through
  an original feature clips it.

## Data

### What exists

- `user_annotations` (one row per imported feature): `ID`, `feature_idx`,
  `raw_name`, `feature_type`, `pos1`, `pos2`, `direction`, `segment_group`,
  `contig_length`, more. Coordinates are in the frame of the user's FASTA.
- `gene_name_map`: `(raw_name, feature_type) -> gene`. NULL `gene` means
  the feature was dropped by the user's mapping (unmapped).
- Curated `annotations` rows for the unit. Rows derived from the import
  carry `tool = "user"`. There is no row-level link back to
  `user_annotations`; matching is by mapped gene + coordinates.
- `assemblies.sequence`: the CURRENT sequence. WF2 (`validate.nf`)
  overwrites it with the curated FASTA, so the user's original sequence is
  not in the DB after WF2.
- `assembly_backup`: written by the first app edit (Linearize or Trim),
  holds the pre-edit (post-WF2) sequence and `ops`. Removed by Restore.
  Helpers: `has_assembly_backup()`, `assembly_backup_ops()` in
  `R/app_annotate_asmb_edits.R`.

### New: store the original contig sequence

New table, created in `new_db_userAnnot()` (`R/init_db_userAnnot.R`) and
filled by `import_user_annotations()`, which already opens each FASTA for
its length:

```
CREATE TABLE user_contigs (
  ID TEXT PRIMARY KEY, sequence TEXT NOT NULL, length INTEGER NOT NULL
);
```

Uppercase sequence, one row per sample (this project type is one contig per
sample by construction). About 16 kb per row.

Projects created on this branch before this change lack the table. The
branch is unreleased; no migration. The app treats a missing table or row
as "frame unknown" and hides the track with the note below. Builder's call
whether `schema_gaps()` should mention it; not required.

## Frame recovery

Module: new file `R/user_annot_track.R`, pure functions, no Shiny.

### Probe search

`find_frame(from, to)` where both are uppercase DNA strings. Returns
`list(offset, rc)` or NULL.

- Take a probe of `k = 60` bases from `from` (start at position 1; if not
  found uniquely, try probes at a few more positions, e.g. every 1/8 of the
  length, so a Trim that removed the first bases still finds an anchor).
- Search the probe in `to` and in `reverseComplement(to)`. Require exactly
  one hit in exactly one orientation; otherwise return NULL (ambiguous:
  repeats, or sequence not derived from `from`).
- Also search across the wrap: when `nchar(from) == nchar(to)` append the
  first `k - 1` bases of `to` to itself before searching, so a probe that
  straddles the new origin still hits.
- `offset` is the signed shift such that position `p` in `from` lands at
  `p + offset` in `to` (before wrapping), `rc` TRUE when the match was on
  the reverse complement.

### Mapping one feature

`map_feature(pos1, pos2, dir, frame, L_from, L_to)`:

- Rotation (same length): `q = ((p - 1 + offset) mod L) + 1` for both ends.
  If `rc`: `q = L - q + 1`, swap ends, flip `dir`. A feature may legally end
  up with `pos1 > pos2` (spans the new origin); both viewers already draw
  those.
- Linear shift (lengths differ, i.e. a Trim happened): `q = p + offset`, no
  wrap. Clip to `[1, L_to]`; mark `clipped = TRUE` when either end moved.
  Drop the feature when it lies entirely outside. An origin-spanning
  feature (pos1 > pos2) in a trimmed, now-linear frame is split by the
  cut: keep the part(s) inside, flagged clipped (builder: simplest correct
  treatment, document the choice in a comment).

### Composition

`original_features(con, ID, path, scaffold)`:

1. `orig` = `user_contigs.sequence` for ID. NULL -> frame unknown.
2. `cur` = `assemblies.sequence` for the unit.
3. If `has_assembly_backup()`: `mid` = `assembly_backup.sequence`, else
   `mid = cur`.
4. `f1 = find_frame(orig, mid)` (WF2 rotation / RC). WF2 never changes the
   sequence length for this project type: `coverage_trim` is forced off for
   imported annotations (`R/annotate.R`), gap filling with `run_tools` adds
   annotation rows only, and curation never edits bases. So step 1 REQUIRES
   `nchar(orig) == nchar(mid)`; a length mismatch here means frame unknown.
5. If `mid` is not `cur`: `f2 = find_frame(mid, cur)` (Linearize rotation
   and/or Trim shift). Linearize-then-Trim is covered because `mid` is the
   pre-Linearize snapshot and `cur` is the result of both; the probe gives
   one offset that is applied as rotation when lengths match. When they do
   not, the offset holds rotation and trim combined, so each feature piece
   is tried shifted by that offset and by a full length either side, and
   the copy that lands inside the sequence is kept (at most one can, for a
   feature that does not straddle the Linearize cut). A feature cut by the
   trim, or straddling the Linearize cut, keeps its larger piece, clipped.
6. Apply `f1` then `f2` to every row of `user_annotations` for ID.
7. Any NULL frame -> return a sentinel (`list(ok = FALSE, reason = ...)`).

Unit keys: user-annotation projects are one contig per sample by
construction (import non-goal), so `user_contigs` is keyed by `ID` alone;
`assemblies` and `assembly_backup` are read for the unit the window is
showing (`rv$updating` ID, path, scaffold).

Probe fallback: try probes starting at 1, then at each 1/8 of the length,
then one ending at the last base. If none hits uniquely, frame unknown; no
exhaustive search.

Cheap: three SQL reads and a couple of `gregexpr`/`matchPattern` calls per
open of the details window. Recomputed on `asmb_edit_tick` and on
`annotations_modal` open, same triggers as the sequence viewer.

## Edit detection

`flag_original(orig_rows, curated)` where `orig_rows` are the mapped
original features joined to `gene_name_map`, `curated` is `rv$annotations`.

Which rows: every row of `user_annotations` for the sample. The parser has
already collapsed gene/CDS/exon hierarchies to one row per feature and
discarded unsupported types (`dropped_types` are never stored), so there is
no further type filter. `D-loop` and `misc_feature` control regions map to
`ctrl`, which curation does emit.

Matching is one-to-one. For each mapped gene, pair original and curated
rows greedily by descending circular overlap (bp); each curated row can be
claimed once. Duplicate genes (two trnL, two trnS, a doubled cox1) and
several raw names mapped to one gene therefore pair with their nearest
copy, not the first one found.

Boundary tolerance is zero. A 3 bp stop-codon completion IS a boundary
change the user wants to see (user's stated primary interest). Builder may
not add a tolerance.

Per original feature, status is the first that applies:

| status      | rule                                                               |
|-------------|--------------------------------------------------------------------|
| `unmapped`  | `gene_name_map.gene` is NULL (dropped by the user's own mapping)   |
| `clipped`   | frame mapping clipped it                                           |
| `unchanged` | a curated row with `pos1 > 0`, same mapped gene, same mapped pos1, pos2 and direction |
| `moved`     | a curated row with the same gene overlaps it (circular overlap)    |
| `dropped`   | nothing else                                                       |

Highlight (edited outline) is drawn for every status except `unchanged`.
`unmapped` is highlighted too (it is absent from the curated track) and the
tooltip says why, so the user can tell a mapping decision from a curation
decision.

Multi-segment originals: each `user_annotations` row is one segment and is
drawn as its own arrow with its own status; curated JOIN members are
separate rows too, so segments pair with segments through the same greedy
overlap rule. A group that curation merged into one row shows the surviving
segment as `moved` and the other as `dropped`; that is acceptable.

## Payload

Both payload builders gain one optional argument:

- `seqview_payload(..., original = NULL)`
- `synteny_payload(..., original = NULL)`

`original` is a data frame (or NULL) with columns `name` (raw name), `gene`
(mapped, may be NA), `type` (feature type of the mapped gene, or raw type
when unmapped, for the hue), `pos1`, `pos2`, `dir`, `status`, `orig1`,
`orig2` (coordinates as the user wrote them, for the tooltip). Emitted as
`p$original = list of features` only when there is at least one row;
absent otherwise, so existing payload tests and the non-user-annot path are
untouched. The synteny builder maps `pos1`/`pos2` to alignment columns in
the browser exactly as sample features are mapped today.

When the frame is unknown the server sends no `original` list and renders
the reason itself as a one-line caption under each viewer's controls
(`uiOutput`), so the JS never has to know about notes. Reasons: "Original sequence not stored for this
sample" and "Original annotation cannot be placed on the edited sequence".

## Drawing

Both `inst/app/www/seqviewer.js` and `inst/app/www/synteny.js`.

- Original features are lane-packed with the shared `window.mpseq.geom.lanes`
  into their own lane block, drawn ABOVE the curated lanes (sequence viewer:
  after the coverage and error tracks, before `drawLanes`; synteny: before
  `sLanes`, with the gutter label `user` next to it as `sample` and `ref`
  are labelled now). Height grows by `nOrigLanes * LANE_H + PAD` when the
  track is on and non-empty.
- Arrow shape identical to curated features. Muted style: fill is the type
  colour at about 20 percent alpha (curated uses `55`, i.e. 33 percent),
  outline `--mp-text-muted` 1 px. Label text in `--mp-text-muted`.
- Edited outline: `--mp-warning` (amber), 2 px, solid. Partial dashing is
  not drawn on the original track (no partial data is shown there).
- Hover tooltip (existing tip element): `Original <raw name> -> <gene or
  "unmapped"> <dir> <pos1>-<pos2> | <status text>`. When mapped coords
  differ from the user's, append `(file: <orig1>-<orig2>)`. Status text:
  "unchanged", "boundaries changed by curation", "dropped by curation",
  "unmapped name, dropped on import", "clipped by trim".
- Click on an original feature: no selection change, no Shiny input. The
  hit still exists for the tooltip.
- Toggle: `mp_checkbox(ns("show_orig"), "Original", value = TRUE)` added to
  both control rows, read by the same `show_*` binding loop the sequence
  viewer already has; synteny gains the same small loop. The checkbox is
  rendered only in user-annotation projects (`session$userData$user_annot`);
  with no checkbox the viewer treats the track as off.
- `window.mpseq.state(id)` gains `nOrigLanes` and `original` (the mapped
  list) so the chromote tests can assert on it, mirroring `features`.

## Server wiring

`R/app_annotate_details.R`:

- One new reactive `orig_track()` next to the seqview wiring: NULL unless
  `isTRUE(session$userData$user_annot)`; otherwise
  `original_features()` + `flag_original()` against `rv$annotations`,
  invalidated by `asmb_edit_tick()`, `gargoyle::watch("annotations_modal")`
  and `rv$annotations`.
- Passed as `original =` to both payload calls. `seqview_server()` gets it
  as one more argument (a reactive), like `selected`.
- `seqview_ui()` and `synteny_viewer_ui()` take a `user_annot` flag to
  render the checkbox.

## Tests

- `tests/testthat/test-user-annot-track.R` (new, pure R):
  - `find_frame`: identity, rotation, rotation + RC, probe over the new
    origin, trimmed start (first probe missing), ambiguous repeat -> NULL,
    unrelated sequence -> NULL.
  - `map_feature`: plus and minus strand through rotation and RC, feature
    that comes to span the origin, trim clipping both ends, feature dropped
    when outside.
  - `original_features` on an in-memory SQLite DB: no backup; backup after
    Linearize; backup after Linearize + Trim; missing `user_contigs` row.
  - `flag_original`: each status once; name-only change is `unchanged`;
    JOIN segments compared separately.
- `test-seqview-payload.R`, `test-synteny-payload.R`: `original` present
  only when given; fields carried; note when frame unknown.
- `test-seqview-js.R`, `test-synteny-js.R` (chromote): payload with
  `original` adds lanes above, toggling the checkbox removes them, clicking
  an original feature does not set the pick input, state exposes the list.
- `test-init-user-annot.R`: `user_contigs` is created and holds the FASTA
  sequence uppercased with the right length.
- Manual check on the demo project (`new_test_project_userAnnot()`):
  GB_DROP and GB_MITOS have a gene spanning the origin; Linearize one of
  them, then Trim, and confirm the original track follows.

## Docs

- `NEWS.md`: one bullet under the development heading.
- `vignettes/Your-Own-Annotations.Rmd`: a short paragraph and one
  screenshot of the two tracks. Screenshot can be taken after the feature
  is built; do not block the merge on it.

## Left to the builder

- Exact alpha values and the 2 px amber outline may be tuned on screen, as
  long as the original track reads as background and the edited outline is
  obvious at a glance.
- Probe length and probe positions, as long as the listed tests pass and
  the demo samples are placed correctly after Linearize + Trim.
- Whether `schema_gaps()` lists `user_contigs`.
- How to split an origin-spanning original across a Trim cut (keep the
  inside part(s), flag clipped).
