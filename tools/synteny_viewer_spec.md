# Interactive synteny viewer for the annotation details window

Design, 2026-10-29. Replaces the static BLAST Reference Synteny PNG and its
separate click-to-zoom panel with one canvas viewer that pans and zooms like the
Sequence viewer (`inst/app/www/seqviewer.js`, spec `tools/nt_viewer_spec.md`).

## 1. Intent

The synteny section today renders a gggenes/patchwork PNG (sample arrows,
rolling identity, reference arrows) at a fixed pixel-per-column scale, with HTML
label overlays stuck on top, a click handler that opens a second ggplot panel
showing a fixed-width base-pair window, and a numeric window-size box. The user
cannot see the whole picture and inspect a region in the same view, the zoom
panel eats vertical space, and its controls (click anchor, window size) differ
from the Sequence viewer's drag/scroll/pinch model one section above.

The new viewer is one canvas: whole alignment at one end of the zoom range,
individual bases of sample and reference at the other, continuous in between,
with the same gestures and button row as the Sequence viewer. It stays a
separate section from the Sequence viewer on purpose: the Sequence viewer's x
axis is the assembled sequence, this viewer's x axis is the pairwise alignment,
which has gaps on both sides.

## 2. Goals, non-goals, anti-goals

Goals

- See gene order of sample vs reference across the full alignment, then zoom
  into any region down to bases without leaving the view.
- Same interaction model as the Sequence viewer: drag or horizontal scroll to
  pan, wheel or pinch to zoom, Fit gene, Whole, zoom in, zoom out buttons,
  click a gene to select its table row, table row selection fits the view.
- Less vertical space than today's overview + zoom panel together.
- Gene labels behave as in the Sequence viewer: drawn inside the arrow when the
  name fits, otherwise shown in the hover tooltip.

Non-goals

- Editing from the viewer. All edits stay in the existing controls.
- Folding synteny into the Sequence viewer.
- Changing how the alignment is computed or stored (`compute_blast_ref_alignment`,
  `normalize_blast_ref_alignment`, `ensure_ref_alignment` are untouched).
- The no-alignment two-track fallback plot. When no alignment exists the section
  shows a one-line message instead (see 5.6).
- Image export, copy-as-FASTA, go-to-position input. Add later if asked.

Anti-goals (would make it a failure even if it works)

- A second set of controls that differ from the Sequence viewer's.
- Any blocking R-side rendering per pan or zoom. All drawing is in the browser
  from one payload per (unit, reference, annotations version).
- Coordinates that disagree with the tables: sample positions must be the
  stored sample coordinates, reference positions the GenBank record's own
  coordinates.

## 3. Constraints

- Shiny app, golem layout, plain JS (no new JS dependency; the Sequence viewer
  is hand-written canvas and this must match it in style and in code).
- Existing data: `rv$blast_ref_aln` (one row: `aligned_sample`, `aligned_ref`,
  `rotation`, `ref_length`, `ref_rc`, already normalised to the sample's stored
  orientation), `rv$blast_ref` (reference gene table: `gene`, `type`, `pos1`,
  `pos2`, `direction`, `ref_length`), `rv$annotations`, `rv$blast_ref_candidates`,
  `active_ref_acc()`.
- Reference picker, "Set ... as best reference", poor-reference flag toggle, the
  linear-topology warning and the "no usable annotations" alert stay exactly as
  they are, above the viewer.
- ASCII only in code. Minimal comments.
- Tests: pure R payload function under testthat; JS geometry under the existing
  chromote harness pattern (`tests/testthat/test-seqview-js.R`).

## 4. Data contract (R to browser)

A pure function `synteny_payload()` in a new file `R/app_annotate_synteny.R`
builds one list; the details module sends it with
`session$sendCustomMessage("mpsyn", payload)` whenever the alignment, the
reference annotations, the sample annotations, or the table selection change.

| field | type | source |
|---|---|---|
| `id` | string | canvas element id (namespaced) |
| `input` | string | Shiny input id for gene clicks |
| `sampleLabel` | string | sample ID |
| `refAcc` | string | active reference accession |
| `refLabel` | string | species (or accession when unknown) |
| `refRc` | bool | `aligned_ref` is reverse-complemented relative to the record |
| `s` | string | `aligned_sample`, upper case |
| `r` | string | `aligned_ref`, upper case |
| `refLength` | int | `ref_length` |
| `rotation` | int | 0-based rotation applied to the reference before aligning |
| `sampleFeatures` | array | one per row of `rv$annotations` with `pos1 > 0`: `row`, `gene`, `type`, `pos1`, `pos2`, `dir` (sample coordinates, as stored) |
| `refFeatures` | array | one per row of `rv$blast_ref`: `gene`, `type`, `pos1`, `pos2`, `dir` (GenBank record coordinates, as stored) |
| `selected` | int or null | `row` of the sample feature to highlight |
| `version` | int | bumps when `s` or `r` change, so the browser keeps its view across annotation-only refreshes and resets on a new alignment |

Coordinates are passed as stored. The browser maps them to alignment columns
with the same arithmetic the R plots use today (`s_to_pct`, `r_to_pct` in
`output$synteny_plot`), see 5.2. Both strings are the full padded alignment
(the reference side is always the whole reference), about 17 kb each for a
mitogenome, so the payload is roughly 40 KB.

## 5. Drawing and behaviour

### 5.1 Layout

One `<canvas>` in a `.mp-seqview` container, full container width, height from
the sum of its rows. Left gutter as in the Sequence viewer, with short row
labels; the long labels (sample ID, reference accession and species, the
"shown reverse-complemented" note in the existing red) go in one caption line
above the canvas.

Rows, top to bottom:

1. Sample ruler: tick labels are sample base positions.
2. Sample gene lanes: arrows, same lane packing, arrow shape, type colours and
   selection highlight as the Sequence viewer. Features that cross the sample
   origin split into two pieces at the alignment's ends.
3. Sample bases: one cell per column. At letter zoom the base letter in the
   Sequence viewer's base palette; a gap column draws "-" in grey. Below letter
   zoom the row collapses to a thin coloured bar, as the Sequence viewer does.
4. Identity band: zoomed out, a rolling-identity area plot binned per pixel
   (mean of match over the non-gap columns in that pixel), green on the light
   grey background used today. Zoomed in past letter threshold, one cell per
   column: green match, red mismatch, grey gap, the same three colours as the
   current zoom panel.
5. Reference bases: as row 3, for the reference string.
6. Reference gene lanes: as row 2. On a flipped reference the arrows mirror
   (a "+" record gene draws "-"), as today.
7. Reference ruler: tick labels are the record's own base positions, undoing
   rotation and flip so they match GenBank.

### 5.2 Coordinate model

The view is linear over alignment columns 1..L. Browser-side pure functions,
exported on `window.mpsyn.geom` for tests:

- `nongap(str)`: array of column indices where `str[i] != "-"`, so
  `nongap[k-1]` is the column of base k.
- Sample base p -> column `sNongap[clamp(p) - 1]`.
- Reference record base p -> rotated `pr = ((p - 1 - rotation) mod refLength) + 1`,
  flipped if `refRc` to `refLength - pr + 1`, then column `rNongap[pr - 1]`.
- A feature whose start column exceeds its end column (wraps the origin, or
  straddles the rotation point) is drawn as two pieces: `[c1, L]` and `[1, c2]`.
  On a flipped reference take min/max of the two mapped ends first.
- Column c -> sample base: count of non-gap sample columns up to c (null on a
  sample gap). Same for reference, then undo flip and rotation for display.
- Per-column class: `gap` if either side is "-", `match` if equal ignoring
  case, else `mismatch`.

Facts the model rests on (verified in `compute_blast_ref_alignment` and the
current `output$synteny_plot`):

- `aligned_sample` with its gaps removed is the whole stored unit sequence, base
  1 to the end, for every unit kind (the aligner is global on the sample side
  and pads the reference window out to the full record). No offset is needed.
- Origin-crossing features in both `rv$annotations` and `rv$blast_ref` arrive
  as one row with `pos1 > pos2`; a feature that straddles the reference
  rotation point maps to start column > end column after rotation. Both cases
  are the two-piece split above. Multi-exon joined genes are separate rows
  (one per piece, `JOIN:` in notes) and each row is its own arrow, as in the
  Sequence viewer.
- `rv$blast_ref` rows are record-forward (`pos1 <= pos2` unless origin-crossing)
  with `direction` separate. The existing `r_to_pct` maps both ends, and on a
  flipped reference takes min/max of the two mapped columns and inverts
  `direction`. Port that exactly; the tests pin it.

Pan and zoom state is `viewStart` (column) and `ppc` (pixels per column),
clamped to `[width / L, 14]` exactly as the Sequence viewer clamps `ppb`,
where L is the full padded alignment length. Whole therefore shows
reference-only stretches with empty sample rows for a short scaffold; that
matches today's reference-width canvas and is intended. Letter threshold 8 px
per column, same constant.

### 5.3 Gestures and buttons

Identical to the Sequence viewer: drag or horizontal wheel pans, vertical wheel
or pinch zooms around the pointer, buttons Fit gene, Whole, zoom in, zoom out
in the same `.mp-seqview-controls` row with the same `data-mpsyn` attribute
wiring, same caption text pattern under the controls.

Initial view on a new alignment version: fit the sample's aligned span (first to
last non-gap sample column) with a small margin, so a short scaffold aligned to
a long reference opens on the scaffold. Whole shows all L columns.

### 5.4 Hover and click

Tooltip (reusing the `.mp-maptoref-tip` element pattern):

- Over a sample arrow: gene, type, strand, sample pos1-pos2.
- Over a reference arrow: gene, type, strand (as drawn), record pos1-pos2.
- Over the base or identity rows: column, sample pos and base (or gap),
  reference record pos and base (or gap), match / mismatch / gap.

Click on a sample arrow: `Shiny.setInputValue(input, {row, nonce}, {priority:
"event"})`; the details module selects that row in the annotation table with
`reactable::updateReactable("table", selected = row)`, the same two-line
observer the Sequence viewer's `pick` uses (`row` is the index into
`rv$annotations`; reactable selection indices are data-row indices, unaffected
by sorting or filtering).
Click on a reference arrow: fit the view to it, no table change.

Table row selection: the module sends `mpsyn_select` with the row; the viewer
fits that sample feature and highlights it, as the Sequence viewer does.

### 5.5 Gene labels

As the Sequence viewer (`seqviewer.js` arrow drawing): gene name drawn centred
inside the arrow when `measureText(gene).width + 8` fits the arrow's on-screen
width, otherwise no label and the name is in the tooltip.

### 5.6 No alignment

When `rv$blast_ref_aln` is null or empty the section shows, in place of the
canvas and its controls, the `mp-coverage-caption` line "No alignment is
available for this reference." The picker and flag controls stay so another
reference can be chosen. The normalised two-track fallback PNG is removed.

## 6. What goes away

In `R/app_annotate_details.R`: `synteny_plot_w`, `synteny_proj`,
`synteny_overlay`, `output$synteny_plot` (renderImage), `output$synteny_labels`,
`output$synteny_zoom_ctrl`, `output$synteny_zoom_ui`, `zoom_overlay`,
`output$synteny_zoom_plot`, `output$synteny_zoom_labels`, `zoom_window_rv` and
its clamp observer, `zoom_click_col`, the `synteny_click` observer, the
`synteny_zoom` toggle observer, both `hScroll` auto-scroll observers, the
`synteny_zoom_window` keydown script and spin-button CSS in the modal footer,
and the `split_wrapped_genes` / `gene_label_overlay` helpers if nothing else
uses them. `@import patchwork` on that file goes; `gggenes`, `ggthemes`, `ragg`
leave DESCRIPTION if no other file uses them (patchwork stays, the coverage
details plot uses it). The `hScroll` message handler in `custom.js` goes if
unused elsewhere.

## 7. Files

- New: `inst/app/www/synteny.js` (viewer), `R/app_annotate_synteny.R`
  (`synteny_payload()`, `synteny_ui()` for the controls + canvas block).
- Shared helpers (`cssVar`, `typeColor`, `niceStep`, arrow drawing) are
  exported from `seqviewer.js` onto `window.mpseq` rather than copied.
- Edited: `R/app_annotate_details.R` (wire payload, remove old code),
  `inst/app/www/custom.css` (only if a new selector is needed), `DESCRIPTION`,
  `NEWS.md` (development version entry), `vignettes/Test-Project-Annotate.Rmd`
  (the synteny paragraph describes click-to-zoom; reword, and re-shoot
  `figures/get-started/annotate-synteny.png`).
- Tests: `tests/testthat/test-synteny-payload.R`,
  `tests/testthat/test-synteny-js.R` with a fixture page under
  `tests/testthat/synteny/`, covering column mapping (rotation, flip, wrap),
  per-column class, identity binning, and the ruler coordinate inversion.

## 8. Left to the builder

- Exact row heights and gutter width, within "less than today's 280 px overview
  plus 180 px zoom panel".
- Whether ruler ticks fall on round sample/reference positions or round
  columns (round positions preferred where cheap).
- Identity binning implementation detail (prefix sum over match counts is the
  obvious one).
- How to reuse the Sequence viewer's arrow and lane code (export vs small
  duplication) as long as the look is identical. The helpers are
  closure-private in `seqviewer.js` today; refactoring that file to export
  them on `window.mpseq` is allowed, provided `test-seqview-js.R` still passes.
