# Interactive Synteny Viewer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the static BLAST synteny PNG and its click-to-zoom panel with one canvas viewer that pans and zooms from whole alignment to bases, with the Sequence viewer's gestures.

**Architecture:** R builds one JSON payload per (unit, reference, annotations) from `rv$blast_ref_aln`, `rv$blast_ref`, `rv$annotations` and sends it with a custom message; a hand-written canvas viewer (`synteny.js`, mirroring `seqviewer.js`) owns all coordinate arithmetic over alignment columns and all drawing. The details module keeps the picker and flag controls and swaps the plot block for the viewer's controls + canvas.

**Tech Stack:** R/Shiny (golem), plain browser JS on canvas, testthat + chromote for JS tests.

**Spec:** `tools/synteny_viewer_spec.md`

## Global Constraints

- Branch `js-synteny-viewer` off `devel`; never push.
- ASCII only in code; minimal comments; no new JS or R dependencies.
- Same gestures, buttons, CSS classes (`.mp-seqview-controls`, `.mp-seqview`, `.mp-seqview-canvas`, `.mp-maptoref-tip`) and constants (letter threshold 8 px, max 14 px per column, bar threshold 3 px) as `seqviewer.js`.
- Sample coordinates are the stored sample coordinates; reference coordinates are GenBank record coordinates (undo rotation and flip for display).
- Reference picker, "Set ... as best reference", poor flag toggle, linear-topology warning and "no usable annotations" alert are untouched.
- `test-seqview-js.R` and `test-seqview-payload.R` keep passing.

## Review Focus

1. Short scaffold vs long reference (sample occupies 4% of columns): initial view must land on the scaffold, Whole must show all columns without error, sample ruler must print nothing over reference-only columns. Pinned in Task 3 (`sampleAt` null on gaps) and Task 4 (initial view test).
2. Flipped reference whose gene straddles the rotation point: must draw as two pieces at the ends, not one arrow across everything (today's R code gets this wrong). Pinned in Task 3 `refFeature` test.
3. Gap-only region under the pointer (both strings "-" cannot happen, but sample gap can): tooltip must say "gap" and never print `NaN`. Pinned in Task 3 `sampleAt`/`refAt` tests and Task 4 tooltip test.
4. Reference with no annotations but an alignment, or annotations but no alignment: section must show the right message and never send a payload with empty strings. Pinned in Task 2 payload test and Task 5 `has_aln` gating.
5. Re-render of the `uiOutput` (picker change) recreates the canvas: the new canvas must get a fresh viewer and the old payload must not be lost. Covered by the `withViewer` retry pattern and `get()` canvas-identity check; Task 4 test loads twice on a replaced canvas.

---

### Task 1: Export shared helpers from the Sequence viewer

**Files:**
- Modify: `inst/app/www/seqviewer.js` (second IIFE, near line 108-118 and the exports at the bottom)
- Test: `tests/testthat/test-seqview-js.R` (existing suite must pass)

**Interfaces:**
- Produces: `window.mpseq.ui = { cssVar(name, fallback), typeColor(type), niceStep(raw) }` and `window.mpseq.geom.lanes` (already exported).

- [ ] **Step 1: Add the export** after the existing `window.mpseq.state = ...` block: `window.mpseq.ui = { cssVar: cssVar, typeColor: typeColor, niceStep: niceStep };`
- [ ] **Step 2: Run** `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-seqview-js.R")'`. Expected: all pass (or skip if no Chrome; then run in a session with Chrome).
- [ ] **Step 3: Commit** `git add inst/app/www/seqviewer.js && git commit -m "Sequence viewer: export shared canvas helpers"`

### Task 2: R payload function

**Files:**
- Create: `R/app_annotate_synteny.R`
- Test: `tests/testthat/test-synteny-payload.R`

**Interfaces:**
- Produces: `synteny_payload(annotations, blast_ref, aln, sample_label, ref_acc, ref_label, selected = NULL, version = 0L)` returning a list with fields `sampleLabel, refAcc, refLabel, refRc, s, r, refLength, rotation, sampleFeatures, refFeatures, selected, version` per spec section 4. Returns `NULL` when `aln` is NULL, has zero rows, or either aligned string is empty.
- Feature objects: sample `list(row, gene, type, pos1, pos2, dir)`; reference `list(gene, type, pos1, pos2, dir)`. Integers as integer, strings as character, `dir` from `direction`.

- [ ] **Step 1: Write failing tests**

```r
# tests/testthat/test-synteny-payload.R
syn_ann <- function() data.frame(
  type = c("PCG", "tRNA", "PCG"), gene = c("nad1", "trnF", "cox1_DELETED_1"),
  pos1 = c(10L, 40L, 0L), pos2 = c(39L, 108L, 0L), direction = c("+", "-", "+"),
  stringsAsFactors = FALSE)
syn_ref <- function() data.frame(
  gene = c("nad1", "trnF"), type = c("PCG", "tRNA"), pos1 = c(5L, 90L), pos2 = c(34L, 20L),
  direction = c("+", "+"), ref_length = 100L, stringsAsFactors = FALSE)
syn_aln <- function(rc = FALSE) data.frame(
  aligned_sample = "acg-t", aligned_ref = "ACGTT", rotation = 4L, ref_length = 100L,
  ref_rc = rc, stringsAsFactors = FALSE)

test_that("payload carries alignment, labels, and live features", {
  p <- synteny_payload(syn_ann(), syn_ref(), syn_aln(TRUE), "S1", "NC_1", "Fish sp.", selected = 2L, version = 3L)
  expect_equal(p[c("sampleLabel", "refAcc", "refLabel")], list(sampleLabel = "S1", refAcc = "NC_1", refLabel = "Fish sp."))
  expect_true(p$refRc)
  expect_equal(p$s, "ACG-T"); expect_equal(p$r, "ACGTT")
  expect_equal(p$refLength, 100L); expect_equal(p$rotation, 4L)
  expect_equal(p$version, 3L); expect_equal(p$selected, 2L)
  expect_equal(vapply(p$sampleFeatures, `[[`, integer(1), "row"), c(1L, 2L))
  expect_equal(p$sampleFeatures[[2]][c("gene", "type", "pos1", "pos2", "dir")],
               list(gene = "trnF", type = "tRNA", pos1 = 40L, pos2 = 108L, dir = "-"))
  expect_equal(p$refFeatures[[2]][c("gene", "pos1", "pos2", "dir")], list(gene = "trnF", pos1 = 90L, pos2 = 20L, dir = "+"))
})

test_that("selection outside the live rows is dropped and missing alignment gives NULL", {
  expect_null(synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x", selected = 3L)$selected)
  expect_false(synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x")$refRc)
  expect_null(synteny_payload(syn_ann(), syn_ref(), NULL, "S1", "NC_1", "x"))
  expect_null(synteny_payload(syn_ann(), syn_ref(), syn_aln()[0, ], "S1", "NC_1", "x"))
  empty <- syn_aln(); empty$aligned_sample <- ""
  expect_null(synteny_payload(syn_ann(), syn_ref(), empty, "S1", "NC_1", "x"))
  expect_length(synteny_payload(syn_ann(), syn_ref()[0, ], syn_aln(), "S1", "NC_1", "x")$refFeatures, 0)
})
```

- [ ] **Step 2: Run** `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-synteny-payload.R")'`. Expected: FAIL, `synteny_payload` not found.
- [ ] **Step 3: Implement** `synteny_payload()` in `R/app_annotate_synteny.R` with a roxygen `@noRd` header pointing at `tools/synteny_viewer_spec.md` section 4. Upper-case both strings; `refRc <- isTRUE(aln$ref_rc[1])`; live rows are `!is.na(pos1) & pos1 > 0 & !is.na(pos2)`.
- [ ] **Step 4: Run** the test file. Expected: PASS.
- [ ] **Step 5: Commit** `git add R/app_annotate_synteny.R tests/testthat/test-synteny-payload.R && git commit -m "Synteny viewer: payload builder"`

### Task 3: JS geometry (pure functions)

**Files:**
- Create: `inst/app/www/synteny.js` (first IIFE: geometry only)
- Create: `tests/testthat/synteny/index.html` (fixture: Shiny shim as in `tests/testthat/seqview/index.html`, a `<details id="sy-section" open>` holding `.mp-seqview-controls` with buttons `data-mpsyn="zoom_in|zoom_out|whole|fit"`, a `.mp-seqview` div with `<canvas id="sy-canvas">` and `<div id="sy-tip" hidden>`, then `<script src=...seqviewer.js>` and `<script src=...synteny.js>`)
- Create: `tests/testthat/test-synteny-js.R` (copy `sv_page()` as `sy_page()`, rewriting both script `src` attributes to the installed files; check `typeof window.mpsyn === "object"`)

**Interfaces:**
- Produces `window.mpsyn.geom`:
  - `nongap(str) -> number[]` 1-based columns of non-gap characters.
  - `counts(str) -> Int32Array` length `str.length + 1`, `counts[c]` = non-gap characters in columns `1..c`.
  - `classify(s, r) -> Uint8Array` per column: 0 match (equal ignoring case), 1 mismatch, 2 gap (either side "-").
  - `sampleCol(p, sNongap) -> number` column of sample base `p` (clamped to `[1, sNongap.length]`).
  - `refCol(p, rNongap, refLength, rotation, rc) -> number` column of record base `p` after rotation and optional flip (spec 5.2).
  - `pieces(c1, c2, L) -> [[a,b],...]` one piece when `c1 <= c2`, else `[[c1, L], [1, c2]]`.
  - `refFeature(f, ctx) -> {c1, c2, dir}` with `ctx = {rNongap, refLength, rotation, rc}`: when `rc`, `c1 = refCol(f.pos2)`, `c2 = refCol(f.pos1)`, `dir` inverted; else `c1 = refCol(f.pos1)`, `c2 = refCol(f.pos2)`, `dir` as given.
  - `sampleAt(c, s, sCounts) -> number|null` sample base at column (null on sample gap).
  - `refAt(c, r, rCounts, refLength, rotation, rc) -> number|null` record position at column (null on reference gap): `k = rCounts[c]`, `pr = rc ? refLength - k + 1 : k`, `((pr - 1 + rotation) % refLength) + 1`.
  - `identity(c0, c1, matchCounts, bothCounts) -> number|null` fraction of matches over columns with both sides non-gap in `[c0, c1]`, null when that count is 0. `matchCounts`/`bothCounts` are prefix arrays built by `prefix(cls, pred)`.
  - `prefix(cls, value) -> Int32Array` prefix count of columns whose class equals `value`; `bothCounts = prefix of cls != 2` (expose as `prefixNot(cls, value)`).

- [ ] **Step 1: Write failing tests** in `test-synteny-js.R` (each via `js(b, ...)` returning JSON):

```r
test_that("column maps follow gaps, rotation and flip", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var s='AC-GT', rr='ACGT-';
    var sn=g.nongap(s), rn=g.nongap(rr);
    return JSON.stringify([sn, rn, g.sampleCol(3,sn), g.sampleCol(9,sn),
      g.refCol(1,rn,4,0,false), g.refCol(1,rn,4,1,false), g.refCol(1,rn,4,0,true),
      Array.from(g.classify(s,rr))]);})()")
  expect_equal(jsonlite::fromJSON(r, simplifyVector = FALSE),
    list(list(1,2,4,5), list(1,2,3,4), 4, 5, 1, 4, 4, list(0,0,2,2,2)))
})
test_that("features split at the ends and flip with the reference", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var rn=[1,2,3,4,5,6,7,8,9,10];
    var fwd=g.refFeature({pos1:2,pos2:4,dir:'+'},{rNongap:rn,refLength:10,rotation:0,rc:false});
    var rc=g.refFeature({pos1:2,pos2:4,dir:'+'},{rNongap:rn,refLength:10,rotation:0,rc:true});
    var straddle=g.refFeature({pos1:9,pos2:2,dir:'-'},{rNongap:rn,refLength:10,rotation:0,rc:true});
    return JSON.stringify([fwd, rc, straddle, g.pieces(straddle.c1,straddle.c2,10), g.pieces(3,5,10)]);})()")
  expect_equal(jsonlite::fromJSON(r, simplifyVector = FALSE),
    list(list(c1=2,c2=4,dir='+'), list(c1=7,c2=9,dir='-'), list(c1=9,c2=2,dir='+'),
         list(list(9,10), list(1,2)), list(list(3,5))))
})
test_that("positions under a column undo flip and rotation, null on a gap", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var s='AC-GT', rr='ACGT-';
    var sc=g.counts(s), rc=g.counts(rr);
    return JSON.stringify([g.sampleAt(3,s,sc), g.sampleAt(4,s,sc), g.refAt(5,rr,rc,4,0,false),
      g.refAt(1,rr,rc,4,2,false), g.refAt(1,rr,rc,4,0,true), g.refAt(2,rr,rc,4,1,true)]);})()")
  expect_equal(jsonlite::fromJSON(r, simplifyVector = FALSE), list(NULL, 3, NULL, 3, 4, 4))
})
test_that("identity is the match fraction over columns with both bases", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var cls=g.classify('ACGTA-','ACCTAA');
    var m=g.prefix(cls,0), both=g.prefixNot(cls,2);
    return JSON.stringify([g.identity(1,6,m,both), g.identity(1,2,m,both), g.identity(6,6,m,both)]);})()")
  expect_equal(jsonlite::fromJSON(r, simplifyVector = FALSE), list(0.8, 1, NULL))
})
```

- [ ] **Step 2: Run** `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-synteny-js.R")'`. Expected: FAIL (`window.mpsyn` undefined).
- [ ] **Step 3: Implement** the geometry IIFE at the top of `inst/app/www/synteny.js`, exporting `window.mpsyn.geom`. Clamp inside `sampleCol`/`refCol` as the R `pmin(pmax(...))` does.
- [ ] **Step 4: Run** the test file. Expected: PASS.
- [ ] **Step 5: Commit** `git add inst/app/www/synteny.js tests/testthat/synteny tests/testthat/test-synteny-js.R && git commit -m "Synteny viewer: column geometry"`

### Task 4: JS viewer (drawing, gestures, messages)

**Files:**
- Modify: `inst/app/www/synteny.js` (second IIFE)
- Modify: `tests/testthat/test-synteny-js.R`

**Interfaces:**
- Consumes: `window.mpsyn.geom` (Task 3), `window.mpseq.ui` and `window.mpseq.geom.lanes` (Task 1), payload shape (Task 2) plus `id`, `input`.
- Produces: Shiny message handlers `mpsyn` (load payload) and `mpsyn_select` (`{id, row}` fits + highlights; `row` null clears); `Shiny.setInputValue(input, {row, nonce}, {priority:'event'})` on sample-arrow click; `window.mpsyn.state(id) -> {L, version, viewStart, ppc, selected, height, sampleLanes, refLanes}`, `window.mpsyn.zoom(id, f)`, `whole(id)`, `fit(id, row)`, `click(id, px, py)`, `hitTest(id, px, py)`, `tip(id) -> string` (current tooltip text).

Constants: `MAX_PPC = 14`, `LETTER = 8`, `BAR = 3`, `RULER_H = 22`, `LANE_H = 22`, `NT_H = 20`, `IDENT_H = 40`, `GUTTER = 60`, `PAD = 4`, `GAP = 6`. Colours: bases from `BASE` palette as `seqviewer.js` (`A #4faf45, C #e0a53f, G #e0555a, T #4a90d9, N #666666`), match `#60BD68`, mismatch `#E55330`, gap `#CCCCCC`, identity background `#F0F0F0`.

Layout (y from top): sample ruler; sample lanes; sample nt row (if `ppc >= BAR`); identity band; ref nt row (if `ppc >= BAR`); ref lanes; ref ruler. Rulers: sample ticks at `niceStep(90/ppc)` round sample positions placed at `x(sampleCol(p))`; ref ticks at round record positions walked in drawn order from the first visible record position, direction `+1` (or `-1` when `rc`) with wrap modulo `refLength`, stopping when the drawn column leaves the view; labels drawn `textBaseline = 'bottom'` at the canvas bottom.

Behaviour:
- `load(p)`: `sameAln = version === p.version`; rebuild nongap/counts/class/prefix arrays only when not `sameAln`; map features to `{c1, c2, dir, gene, type, row?, pos1, pos2}` plus `pieces`; lane-pack each track with `mpseq.geom.lanes` on `{pos1: c1, pos2: c2}` with `len = L`; if not `sameAln` set the initial view to the sample span `[sNongap[0], sNongap[last]]` with 2% margin (via `fitCols(a, b)`); draw.
- `clamp`, `zoom`, `whole`, `fitCols`, `fit(row)` as in `seqviewer.js` with linear topology semantics.
- Identity: `ppc >= BAR` -> one cell per column coloured by class; else per pixel column `identity(c0, c1)` -> filled area height `IDENT_H * id`, skip null.
- Nt rows: `ppc >= LETTER` letters on base-coloured cells, gap columns draw "-" in grey on no fill; `BAR <= ppc < LETTER` coloured cells only.
- Arrows: same polygon, fill `col + '55'`, stroke `col` (selected: `--mp-primary`, width 2), label when `x1 - x0 > measureText(gene).width + 8`. Hits: `{x0, x1, y0, y1, f, track: 'sample'|'ref'}`.
- Gestures and buttons exactly as `seqviewer.js` `bindControls` minus copy/goto/checkboxes; attribute `data-mpsyn`.
- Hover text: arrow -> `gene (type) strand pos1-pos2` prefixed `Sample` or `Reference`; otherwise `Column c | sample P B | ref P B | match|mismatch|gap`, printing `gap` in place of `P B` on a null position.
- Click: sample hit -> select + setInputValue; ref hit -> `fitCols(piece)`; empty -> clear selection.

- [ ] **Step 1: Write failing tests** appended to `test-synteny-js.R`:

```r
sy_payload <- function(extra = "") sprintf("{id:'sy-canvas', input:'sy-pick', sampleLabel:'S1', refAcc:'NC_1', refLabel:'Fish', refRc:false,
  s:'----' + Array(200).join('ACGT').slice(0,400) + '----', r:Array(200).join('ACGT').slice(0,408), refLength:408, rotation:0, version:1, selected:null,
  sampleFeatures:[{row:1,gene:'nad1',type:'PCG',pos1:10,pos2:120,dir:'+'},{row:2,gene:'trnF',type:'tRNA',pos1:390,pos2:30,dir:'-'}],
  refFeatures:[{gene:'nad1',type:'PCG',pos1:14,pos2:124,dir:'+'}] %s}", extra)
test_that("the viewer opens on the sample span, fits, zooms and reports clicks", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); var s0 = window.mpsyn.state('sy-canvas');
    window.__handlers.mpsyn_select({id:'sy-canvas', row:1}); var s1 = window.mpsyn.state('sy-canvas');
    window.mpsyn.whole('sy-canvas'); var s2 = window.mpsyn.state('sy-canvas');
    window.mpsyn.zoom('sy-canvas', 2); var s3 = window.mpsyn.state('sy-canvas');
    window.mpsyn.zoom('sy-canvas', 100); var s4 = window.mpsyn.state('sy-canvas');
    var h = window.mpsyn.hitTest('sy-canvas', 100, 30);
    return JSON.stringify({s0:s0, s1:s1, s2:s2, s3:s3, s4:s4, h:h ? h.f.gene : null, inputs: window.__inputs});})()", sy_payload()))
  x <- jsonlite::fromJSON(r, simplifyVector = FALSE)
  expect_equal(x$s0$L, 408); expect_lt(x$s0$viewStart, 5); expect_gt(x$s0$ppc, x$s2$ppc)
  expect_equal(x$s1$selected, 1); expect_true(x$s1$viewStart <= 14 && x$s1$viewStart + 940 / x$s1$ppc >= 124)
  expect_equal(x$s2$viewStart, 1); expect_equal(x$s3$ppc, x$s2$ppc * 2); expect_equal(x$s4$ppc, 14)
  expect_equal(x$s0$sampleLanes, 1); expect_equal(x$s0$refLanes, 1)
})
test_that("clicking a sample arrow selects its row; a reference arrow only fits", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); window.mpsyn.whole('sy-canvas');
    var v = window.mpsyn.state('sy-canvas'); var xs = 60 + (60 - v.viewStart) * v.ppc;
    window.mpsyn.click('sy-canvas', xs, 22 + 4 + 8); var a = window.mpsyn.state('sy-canvas').selected;
    var n0 = window.__inputs.length;
    window.mpsyn.click('sy-canvas', xs, v.height - 22 - 4 - 8); var s = window.mpsyn.state('sy-canvas');
    return JSON.stringify({a:a, n0:n0, n1:window.__inputs.length, name:window.__inputs[0].name, row:window.__inputs[0].value.row, ppcAfter:s.ppc});})()", sy_payload()))
  x <- jsonlite::fromJSON(r)
  expect_equal(x$a, 1); expect_equal(x$n0, 1); expect_equal(x$n1, 1); expect_equal(x$name, "sy-pick"); expect_equal(x$row, 1)
  expect_gt(x$ppcAfter, 940 / 408)
})
test_that("the tooltip names gaps and never prints NaN", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); window.mpsyn.zoom('sy-canvas', 100); window.mpsyn.whole('sy-canvas'); window.mpsyn.zoom('sy-canvas', 100);
    var v = window.mpsyn.state('sy-canvas'); v.viewStart = 1; return JSON.stringify([window.mpsyn.tip('sy-canvas', 1), window.mpsyn.tip('sy-canvas', 5)]);})()", sy_payload()))
  x <- jsonlite::fromJSON(r)
  expect_match(x[1], "sample gap"); expect_match(x[1], "ref 1 A"); expect_match(x[2], "sample 1 A"); expect_false(any(grepl("NaN", x)))
})
test_that("a replaced canvas gets a fresh viewer", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); var old = document.getElementById('sy-canvas');
    var nw = document.createElement('canvas'); nw.id = 'sy-canvas'; old.parentNode.replaceChild(nw, old);
    window.__handlers.mpsyn(%s); return JSON.stringify(window.mpsyn.state('sy-canvas').L);})()", sy_payload(), sy_payload()))
  expect_equal(jsonlite::fromJSON(r), 408)
})
```

`window.mpsyn.tip(id, col)` returns the hover text for column `col` (a pure helper used by `hover`), so the test needs no mouse events.

- [ ] **Step 2: Run** the test file. Expected: new tests FAIL (`mpsyn` handler undefined).
- [ ] **Step 3: Implement** the viewer IIFE in `synteny.js` per the behaviour list above.
- [ ] **Step 4: Run** the test file and `test-seqview-js.R`. Expected: PASS.
- [ ] **Step 5: Commit** `git add inst/app/www/synteny.js tests/testthat/test-synteny-js.R && git commit -m "Synteny viewer: canvas viewer"`

### Task 5: Wire into the details module and remove the old plot

**Files:**
- Modify: `R/app_annotate_synteny.R` (add `synteny_viewer_ui(ns)`)
- Modify: `R/app_annotate_seqview.R` (return `sel = cur_sel` alongside `pick`)
- Modify: `R/app_annotate_details.R` (synteny block ~1083-2000, modal footer ~5210-5235, helpers 6-100, `@import patchwork` line 3)
- Modify: `inst/app/www/custom.js` (remove `hScroll` handler and the stale comment at ~154)
- Modify: `DESCRIPTION`, `NAMESPACE` (via `devtools::document()`)

**Interfaces:**
- Produces: `synteny_viewer_ui(ns)` returning a `tagList` of: `.mp-seqview-controls` div with buttons `data-mpsyn` = `fit` ("Fit gene"), `whole` ("Whole alignment"), `zoom_in`, `zoom_out` (icons + titles as `seqview_ui`); a `.mp-coverage-caption` with "Drag or scroll sideways to pan, scroll or pinch to zoom; click a sample gene to select its row. Bases appear when zoomed in."; a `.mp-seqview` div holding `tags$canvas(id = ns("synteny_canvas"), class = "mp-seqview-canvas")` and `tags$div(id = ns("synteny_tip"), class = "mp-maptoref-tip", hidden = NA)`.
- Consumes: `synteny_payload()` (Task 2), `sv$sel()` from `seqview_server`.

- [ ] **Step 1: `seqview_server`** returns `list(pick = reactive(input$pick), sel = cur_sel)`.
- [ ] **Step 2: In `output$synteny_ui`** replace the caption line (`lbl_col` column, `synteny_zoom_ctrl` row, `syntenyScrollDiv` image block and the `conditionalPanel` zoom block) with: a caption div (`sample_lbl`, " vs ", `ref_acc`, `ref_lbl`, and the existing red "shown reverse-complemented" note when `ref_rc`), then `if (has_aln) synteny_viewer_ui(ns) else div(class = "mp-coverage-caption", "No alignment is available for this reference.")`. Keep `has_ref` gating: the message also shows when `!has_ref`.
- [ ] **Step 3: Add the payload observer** next to the old render code:

```r
syn_version <- reactiveVal(0L); syn_last <- NULL
observeEvent(rv$blast_ref_aln, {
  a <- rv$blast_ref_aln
  key <- if (is.null(a) || nrow(a) == 0) "" else paste(a$aligned_sample[1], a$aligned_ref[1])
  if (!identical(key, syn_last)) { syn_last <<- key; syn_version(isolate(syn_version()) + 1L) }
}, ignoreNULL = FALSE)
observe({
  ctx <- req(fig_ctx()); req(rv$annotations, rv$blast_ref)
  acc <- active_ref_acc() %||% ctx$blast_accession
  cand <- rv$blast_ref_candidates
  lbl <- if (!is.null(cand)) cand$species[match(acc, cand$accession)] else NA
  p <- synteny_payload(rv$annotations, rv$blast_ref, rv$blast_ref_aln, ctx$ID, acc,
                       if (is.na(lbl) || !nzchar(lbl)) acc else lbl,
                       selected = sv$sel(), version = syn_version())
  req(p)
  p$id <- ns("synteny_canvas"); p$input <- ns("synteny_pick")
  session$sendCustomMessage("mpsyn", p)
})
observeEvent(sv$sel(), session$sendCustomMessage("mpsyn_select", list(id = ns("synteny_canvas"), row = sv$sel())), ignoreNULL = FALSE)
observeEvent(input$synteny_pick, reactable::updateReactable("table", selected = as.integer(input$synteny_pick$row)))
```

`sv` is defined at ~2854; move the synteny observers below it or reference `sv` lazily (observers run after module init, so placement after `sv <-` is enough).

- [ ] **Step 4: Delete** everything in spec section 6: `synteny_plot_w`, `synteny_proj`, `synteny_overlay`, `output$synteny_plot`, `output$synteny_labels`, `output$synteny_zoom_ctrl`, `output$synteny_zoom_ui`, `zoom_overlay`, `output$synteny_zoom_plot`, `output$synteny_zoom_labels`, `zoom_window_rv` + its observer, `zoom_click_col` (and its two uses in `selected`), the `synteny_click` observer, the `synteny_zoom` toggle observer, both `hScroll` observers, the modal footer `tags$script`/`tags$style` for `synteny_zoom_window` (keep the `annotation_table_details` toggle part of that script), `split_wrapped_genes`, `gene_label_overlay`, `gene_type_fill`, `gene_type_alpha`, and the `@import patchwork` line. Remove the `hScroll` handler and the "annotate-details views ... native scrolling" comment from `custom.js`. Remove `gggenes`, `ggthemes`, `ragg` from DESCRIPTION Imports.
- [ ] **Step 5: Verify nothing else used them:** `grep -rn "gggenes\|ggthemes\|ragg::\|gene_type_fill\|split_wrapped_genes\|gene_label_overlay\|hScroll\|zoom_click_col\|synteny_zoom" R/ inst/app/www/*.js` prints nothing.
- [ ] **Step 6: Run** `Rscript -e 'devtools::document(); devtools::load_all(); testthat::test_dir("tests/testthat", filter = "synteny|seqview")'` then `Rscript -e 'devtools::load_all()'` for parse errors and `R CMD build`-free check `Rscript -e 'devtools::check_man()'`. Expected: PASS, no NAMESPACE import of patchwork, no warnings about undefined globals from removed code (`Rscript -e 'devtools::load_all(); codetools::checkUsageEnv(asNamespace("MitoPilot"))' | grep -i synteny` prints nothing).
- [ ] **Step 7: Commit** `git add -A R inst DESCRIPTION NAMESPACE && git commit -m "Annotate details: interactive synteny viewer replaces static plot"`

### Task 6: Smoke test in the running app

**Files:** none (verification)

- [ ] **Step 1:** Open the app on the shipped fish test project (`tests/testthat` fixtures or `~/MitoPilot_scratch` per memory `scyphozoa-test-project` / `docs-screenshot-refresh`), open an annotation details window for a sample with a BLAST reference.
- [ ] **Step 2:** Check: viewer opens on the sample span; wheel zoom reaches letters with coloured mismatches; drag pans; Whole shows the full reference; selecting a table row fits the gene; clicking a sample arrow selects the row; switching the reference picker re-renders with the new reference; a reference with no alignment shows the message; no console errors.
- [ ] **Step 3:** Fix anything found, re-run the Task 4 and 5 test commands, commit.

### Task 7: NEWS and vignette

**Files:**
- Modify: `NEWS.md` (top of development version), `vignettes/Test-Project-Annotate.Rmd:229-240`
- Re-shoot: `vignettes/figures/get-started/annotate-synteny.png` using the `dev/docs_shots` harness if present and working; otherwise leave the PNG and say so in the final report.

- [ ] **Step 1: NEWS** entry: "**Annotate details.** The BLAST Reference Synteny plot is now an interactive viewer: drag or scroll to pan, scroll or pinch to zoom from the whole alignment down to individual bases, with mismatches and gaps coloured. The separate click-to-zoom panel and window-size box are gone, and the viewer takes less space."
- [ ] **Step 2: Vignette** paragraph: replace "Click anywhere in this plot to show a zoomed-in base-pair level alignment of your sample versus the reference." with "Drag or scroll sideways to pan and scroll or pinch to zoom; zoomed in far enough, the bases of your sample and the reference appear with mismatches and gaps coloured. Click a gene in the sample track to select it in the table."
- [ ] **Step 3: Screenshot** if the harness works; `git add` the PNG.
- [ ] **Step 4: Commit** `git add NEWS.md vignettes && git commit -m "Docs: interactive synteny viewer"`
