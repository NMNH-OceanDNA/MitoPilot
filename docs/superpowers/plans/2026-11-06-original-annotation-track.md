# Original User-Annotation Track Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Draw the user's original imported annotation as a muted, toggleable track above the curated lanes in the sequence and synteny viewers, placed in the current coordinate frame, with curation-edited features outlined in amber.

**Architecture:** Import stores each sample's original contig sequence. A pure R module recovers the frame by probe search (original -> WF2 output -> current after app edits), maps every `user_annotations` row, and flags it against `rv$annotations`. The two payload builders carry the result as `p$original`; the two canvas viewers draw it as an extra lane block.

**Tech Stack:** R (Shiny, DBI/RSQLite, Biostrings), vanilla JS canvas, testthat + chromote for the JS.

**Spec:** `docs/superpowers/specs/2026-11-06-original-annotation-track-design.md`

## Global Constraints

- Branch `user-annot-import`; commit per task; never push; no Claude attribution in commits.
- ASCII only in R and JS source. Minimal comments. No em dashes anywhere.
- Boundary tolerance for "unchanged" is exactly 0 bp. Name changes never count as edits.
- Track appears only when `isTRUE(session$userData$user_annot)`.
- Original features are never clickable for row selection and never send a Shiny input.
- When the frame cannot be recovered, send no `original` list and show a one-line note; never draw at guessed positions.
- Muted style: type colour fill at alpha `33`, outline `--mp-text-muted` 1 px, label `--mp-text-muted`. Edited outline `--mp-warning` 2 px solid.
- Run R tests with: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/<file>")'`

## Review Focus

1. An original feature that spans the origin (pos1 > pos2) after rotation: both viewers must draw it as two pieces, not a full-width bar. Pinned in Task 2 (`map_feature` wrap case) and Task 5 (chromote lane test with a wrapped original).
2. Reverse-complement rotation: a `+` original must come out `-` with ends swapped, and an unchanged gene must still flag `unchanged` against the curated row. Pinned in Task 2 and Task 3.
3. Two trnL (or two of any gene): each original must pair with its own curated copy; a false `moved` on one of them is a bug. Pinned in Task 3.
4. Trim that cuts through an original feature: clipped, flagged, still drawn inside the sequence; a feature wholly outside is absent, and the track still renders for the rest. Pinned in Task 2 and Task 3 (Linearize + Trim DB case).
5. Project older than this change (no `user_contigs` table) or a sample without a row: the window must open normally with the note, not error. Pinned in Task 3 and Task 7.

---

### Task 1: Store the original contig sequence at import

**Files:**
- Modify: `R/init_db_userAnnot.R` (table DDL next to `user_annotations`, ~line 48; fill in `import_user_annotations()`, ~line 87)
- Test: `tests/testthat/test-init-user-annot.R`

**Interfaces:**
- Produces: table `user_contigs (ID TEXT PRIMARY KEY, sequence TEXT NOT NULL, length INTEGER NOT NULL)`, one row per sample, sequence uppercased, written in the same transaction scope as `user_annotations`.

- [ ] **Step 1: Write the failing test** (append to the existing "builds tables" test or add a new one)

```r
test_that("new_project_userAnnot stores each original contig in user_contigs", {
  fx <- ua_fixture()   # whatever helper the file already uses to build the project
  con <- DBI::dbConnect(RSQLite::SQLite(), fx$db); withr::defer(DBI::dbDisconnect(con))
  uc <- DBI::dbReadTable(con, "user_contigs")
  expect_setequal(uc$ID, c("UA_LINEAR", "UA_CIRCULAR"))
  fa <- Biostrings::readDNAStringSet(file.path(fx$asm, "UA_CIRCULAR.fasta"))
  expect_identical(uc$sequence[uc$ID == "UA_CIRCULAR"], toupper(as.character(fa[[1]])))
  expect_identical(uc$length[uc$ID == "UA_CIRCULAR"], Biostrings::width(fa)[1])
})
```

- [ ] **Step 2: Run it, expect FAIL** ("no such table: user_contigs").
- [ ] **Step 3: Implement.** Add the `CREATE TABLE user_contigs` statement after the `gene_name_map` DDL. In `import_user_annotations()`, replace the `fasta.seqlengths` call with `Biostrings::readDNAStringSet()` of the same file, derive `L` from its width, and after `dbAppendTable(con, "user_annotations", ua)` append one row per sample to `user_contigs` with `toupper(as.character(seq[[1]]))`.
- [ ] **Step 4: Run `test-init-user-annot.R`, expect all PASS.**
- [ ] **Step 5: Commit** `feat: store original contig sequence for user annotation projects`

---

### Task 2: Frame recovery by probe search

**Files:**
- Create: `R/user_annot_track.R`
- Test: `tests/testthat/test-user-annot-track.R` (new)

**Interfaces:**
- Produces:
  - `find_frame(from, to, k = 60L)`: both uppercase DNA strings. Returns `list(offset = <integer>, rc = <logical>)` or `NULL`. Meaning: position `p` of `from` lands at `p + offset` in `to` (or in `reverseComplement(to)` when `rc`), wrapped modulo `nchar(to)` when the lengths are equal.
  - `map_feature(pos1, pos2, dir, frame, L_from, L_to)`: returns `list(pos1, pos2, dir, clipped = <logical>)` or `NULL` when the feature falls wholly outside after a linear shift.
  - `revcomp_chr(x)`: uppercase string -> reverse complement string (thin wrapper over Biostrings).

- [ ] **Step 1: Write failing tests**

```r
rot <- function(s, n) paste0(substr(s, n + 1, nchar(s)), substr(s, 1, n))
set.seed(1); S <- paste(sample(c("A","C","G","T"), 2000, TRUE), collapse = "")

test_that("find_frame: identity, rotation, rc, probe over the new origin", {
  expect_equal(find_frame(S, S), list(offset = 0L, rc = FALSE))
  expect_equal(find_frame(S, rot(S, 500)), list(offset = -500L, rc = FALSE))
  expect_equal(find_frame(S, revcomp_chr(S)), list(offset = 0L, rc = TRUE))
  expect_equal(find_frame(S, rot(S, 1990)), list(offset = -1990L, rc = FALSE))  # probe straddles origin
})
test_that("find_frame: trimmed start uses a later probe; repeats and unrelated give NULL", {
  expect_equal(find_frame(S, substr(S, 101, 1900)), list(offset = -100L, rc = FALSE))
  rep2 <- paste0(S, S)
  expect_null(find_frame(S, rep2))
  expect_null(find_frame(S, paste(rev(strsplit(S, "")[[1]]), collapse = "")))
})
test_that("map_feature: rotation keeps strand, rc flips, wrap allowed", {
  f <- list(offset = -500L, rc = FALSE)
  expect_equal(map_feature(600L, 700L, "+", f, 2000L, 2000L)[c("pos1","pos2","dir")], list(pos1 = 100L, pos2 = 200L, dir = "+"))
  expect_equal(map_feature(400L, 600L, "+", f, 2000L, 2000L)[c("pos1","pos2")], list(pos1 = 1900L, pos2 = 100L))
  r <- list(offset = 0L, rc = TRUE)
  expect_equal(map_feature(10L, 30L, "+", r, 2000L, 2000L)[c("pos1","pos2","dir")], list(pos1 = 1971L, pos2 = 1991L, dir = "-"))
})
test_that("map_feature: linear shift clips, flags, and drops", {
  f <- list(offset = -100L, rc = FALSE)
  expect_equal(map_feature(50L, 150L, "+", f, 2000L, 1800L), list(pos1 = 1L, pos2 = 50L, dir = "+", clipped = TRUE))
  expect_false(map_feature(500L, 600L, "+", f, 2000L, 1800L)$clipped)
  expect_null(map_feature(10L, 90L, "+", f, 2000L, 1800L))
  expect_equal(map_feature(1950L, 300L, "+", f, 2000L, 1800L)[c("pos1","pos2","clipped")], list(pos1 = 1L, pos2 = 200L, clipped = TRUE))  # keeps the larger surviving piece
})
```

- [ ] **Step 2: Run, expect FAIL** (functions not found).
- [ ] **Step 3: Implement in `R/user_annot_track.R`.** Probe positions: `1`, then `round(L * (1:7) / 8)`, then `L - k + 1`. Search with `gregexpr(probe, target, fixed = TRUE)` on `target` and on `revcomp_chr(target)`; when `nchar(from) == nchar(to)` search `paste0(target, substr(target, 1, k - 1))` instead so a probe over the new origin hits. Accept only exactly one hit in exactly one orientation. In `map_feature`, use the existing `wrap_pos()` from `R/annotate_utils.R` for the circular case; for the linear case split a `pos1 > pos2` input at the old origin and keep the larger surviving piece.
- [ ] **Step 4: Run, expect PASS.**
- [ ] **Step 5: Commit** `feat: recover original annotation frame by probe search`

---

### Task 3: Original features from the DB and edit flags

**Files:**
- Modify: `R/user_annot_track.R`
- Test: `tests/testthat/test-user-annot-track.R`

**Interfaces:**
- Consumes: Task 2 functions; `has_assembly_backup()`, `assembly_backup_ops()` (`R/app_annotate_asmb_edits.R`); `circ_overlap()`, `mito_gene_type()` (`R/annotate_utils.R`).
- Produces:
  - `original_features(con, ID, path, scaffold)`: returns `list(ok = TRUE, rows = <data.frame>)` or `list(ok = FALSE, reason = <string>)`. `rows` columns: `name` (raw_name), `gene` (mapped, NA when unmapped), `type` (`PCG`/`tRNA`/`rRNA`/`ctrl`, or from `feature_type` `CDS`->`PCG`, `tRNA`, `rRNA`, else `""`), `pos1`, `pos2`, `dir` (current frame), `orig1`, `orig2`, `odir` (as in the file), `clipped`, `segment_group`. Reasons (exact copy): `"Original sequence not stored for this sample"`, `"Original annotation cannot be placed on the edited sequence"`.
  - `flag_original(rows, curated)`: adds `status` in `c("unmapped","clipped","unchanged","moved","dropped")` using the spec's first-match order; `curated` is `rv$annotations` (uses `gene`, `pos1`, `pos2`, `direction`, rows with `pos1 > 0` only). Pairing is greedy one-to-one by descending circular overlap within a gene.

- [ ] **Step 1: Write failing tests** using an in-memory SQLite DB built by a small helper in the test file:

```r
track_db <- function(orig, cur, backup = NULL, ua, gm, with_contigs = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies (ID TEXT, path INTEGER, scaffold INTEGER, sequence TEXT, length INTEGER, topology TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES ('S1',1,1,?,?,'circular')", params = list(cur, nchar(cur)))
  if (with_contigs) { DBI::dbExecute(con, "CREATE TABLE user_contigs (ID TEXT PRIMARY KEY, sequence TEXT, length INTEGER)")
    DBI::dbExecute(con, "INSERT INTO user_contigs VALUES ('S1',?,?)", params = list(orig, nchar(orig))) }
  if (!is.null(backup)) { ensure_assembly_backup(con)
    DBI::dbExecute(con, "INSERT INTO assembly_backup (ID,path,scaffold,sequence,length,ops) VALUES ('S1',1,1,?,?,'linearize')", params = list(backup, nchar(backup))) }
  DBI::dbWriteTable(con, "user_annotations", ua); DBI::dbWriteTable(con, "gene_name_map", gm)
  con
}
ua_rows <- function(...) data.frame(ID = "S1", feature_idx = seq_along(list(...)[[1]]), raw_name = ..1, feature_type = "CDS", ...)[-1]  # builder: shape as needed
```

Tests (assert on the named outputs; use `S` and `rot()` from Task 2):

```r
test_that("original_features: no edits, WF2 rotation only", {
  # cur = rot(S, 500); one feature COI 600-700 + ; expect pos1 100 pos2 200, orig1 600, type PCG, gene cox1
})
test_that("original_features: backup after Linearize, then Trim", {
  # orig S; backup rot(S, 500); cur substr(rot(S,500), 51, 1800); feature 600-700 -> 50-150, clipped FALSE; feature 520-580 -> 1-30 clipped TRUE
})
test_that("original_features: missing table / row give the stored reasons", {
  # with_contigs = FALSE -> ok FALSE, reason "Original sequence not stored for this sample"
  # cur = unrelated sequence -> ok FALSE, reason "Original annotation cannot be placed on the edited sequence"
})
test_that("flag_original: one status each; name-only change is unchanged; duplicates pair one-to-one", {
  # unchanged: COI 100-200 + vs curated cox1 100-200 +
  # moved: 100-203 (3 bp stop completion) vs 100-200
  # dropped: no curated row of that gene overlapping
  # unmapped: gene NA
  # clipped: clipped TRUE wins over coordinates
  # two trnL originals 300-370 and 900-970; curated trnL 300-370 and 900-975 -> unchanged, moved (not moved, unchanged)
  # rc frame: original + becomes -, curated - at the mapped coords -> unchanged
})
```

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.** `original_features()`: read `user_contigs` row (missing table or row -> reason 1); read `assemblies.sequence` for the unit; `mid` = backup sequence when `has_assembly_backup()` else current; `f1 = find_frame(orig, mid)` requiring equal lengths; `f2 = find_frame(mid, cur)` only when a backup exists; any NULL -> reason 2; map every `user_annotations` row through `f1` then `f2`; join `gene_name_map` on `(raw_name, feature_type)`. `flag_original()`: per gene, compute circular overlap bp for every original x curated pair with `circ_overlap()` and lengths, sort descending, claim greedily.
- [ ] **Step 4: Run, expect PASS.**
- [ ] **Step 5: Commit** `feat: place and flag original user annotations in the current frame`

---

### Task 4: Payloads carry the original track

**Files:**
- Modify: `R/app_annotate_seqview.R` (`seqview_payload`, lines 9-55), `R/app_annotate_synteny.R` (`synteny_payload`, lines 8-37)
- Test: `tests/testthat/test-seqview-payload.R`, `tests/testthat/test-synteny-payload.R`

**Interfaces:**
- Consumes: `flag_original()` output shape (Task 3).
- Produces: both builders accept `original = NULL`. When a data frame with rows is given, `p$original` is an unnamed list of `list(name, gene, type, pos1, pos2, dir, status, orig1, orig2, odir)` with `gene` as `NULL` when NA. No `original` element otherwise.

- [ ] **Step 1: Write failing tests** (one per file)

```r
test_that("original features are carried only when given", {
  o <- data.frame(name = "COI", gene = "cox1", type = "PCG", pos1 = 10L, pos2 = 39L, dir = "+",
                  status = "unchanged", orig1 = 10L, orig2 = 39L, odir = "+", stringsAsFactors = FALSE)
  p <- seqview_payload(sv_ann(), "ACGT", "linear", "S1.1.1", original = o)
  expect_equal(p$original[[1]]$name, "COI"); expect_equal(p$original[[1]]$status, "unchanged")
  expect_null(seqview_payload(sv_ann(), "ACGT", "linear", "S1.1.1")$original)
  expect_null(seqview_payload(sv_ann(), "ACGT", "linear", "S1.1.1", original = o[0, ])$original)
})
```
Same shape for `synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x", original = o)`; also assert an unmapped row (`gene = NA`) emits `gene = NULL`.

- [ ] **Step 2: Run both files, expect the new tests FAIL.**
- [ ] **Step 3: Implement** a shared helper `original_payload(original)` in `R/user_annot_track.R` returning the list or NULL; both builders call it.
- [ ] **Step 4: Run both payload test files, expect PASS.**
- [ ] **Step 5: Commit** `feat: original track in viewer payloads`

---

### Task 5: Sequence viewer draws the original track

**Files:**
- Modify: `inst/app/www/seqviewer.js` (Viewer ctor ~120-138, `load` 140-156, `topH`/`laneY`/`height` 288-297, `draw` 298-328, new `drawOrigLanes` next to `drawLanes` 383-425, `hover` 598-612, `click`, `bindControls` checkbox loop 643-646, `state` 677-682)
- Modify: `tests/testthat/seqview/index.html` (add `<input type="checkbox" id="sv-show_orig" checked>`)
- Test: `tests/testthat/test-seqview-js.R`

**Interfaces:**
- Consumes: `p.original` (Task 4).
- Produces: `v.orig` (array), `v.nOrigLanes`, `v.showOrig`, `v.origOn()`; original lane block sits between the error track and the curated lanes; `laneY(lane)` of curated lanes shifts down by `nOrigLanes * LANE_H + PAD` when `origOn()`; `window.mpseq.state(id)` gains `nOrigLanes` and `original: [{name, pos1, pos2, lane, status}]`; hits for originals carry `orig: true` and `click` ignores them.

- [ ] **Step 1: Write failing chromote tests**

```r
test_that("original track packs above the curated lanes, toggles off, and is not clickable", {
  b <- sv_page()
  r <- js(b, "(function(){
    var seq = Array(2000).join('ACGT').slice(0, 2000);
    window.__handlers.mpseq({id:'sv-canvas', input:'sv-pick', unit:'S1.1.1', len:2000, topology:'circular', seq:seq, version:1, selected:null,
      features:[{row:1,type:'PCG',gene:'nad1',pos1:100,pos2:399,dir:'+',partial5:false,partial3:false,notes:''}],
      original:[{name:'ND1',gene:'nad1',type:'PCG',pos1:100,pos2:402,dir:'+',status:'moved',orig1:100,orig2:402,odir:'+'},
                {name:'CYTB',gene:'cob',type:'PCG',pos1:1950,pos2:60,dir:'+',status:'dropped',orig1:1950,orig2:60,odir:'+'}]});
    var s0 = window.mpseq.state('sv-canvas');
    var y = window.mpseq.laneY('sv-canvas', 0);
    var before = window.__inputs.length;
    var h = window.mpseq.hitTest('sv-canvas', 200, y - 22);   // inside the original lane above lane 0
    window.mpseq.click('sv-canvas', 200, y - 22);
    var box = document.getElementById('sv-show_orig'); box.checked = false; box.dispatchEvent(new Event('change'));
    var s1 = window.mpseq.state('sv-canvas');
    return JSON.stringify({n0:s0.nOrigLanes, n1:s1.nOrigLanes, orig:s0.original.length, hitOrig: !!(h && h.orig),
      sent: window.__inputs.length - before, sel: s1.selected, taller: s0.height > s1.height});})()")
  r <- jsonlite::fromJSON(r)
  expect_equal(r$n0, 1L); expect_equal(r$n1, 0L); expect_equal(r$orig, 2L)
  expect_true(r$hitOrig); expect_equal(r$sent, 0L); expect_null(r$sel); expect_true(r$taller)
})
```

- [ ] **Step 2: Run, expect FAIL** (`nOrigLanes` undefined).
- [ ] **Step 3: Implement.** Lane-pack `this.orig` with `G.lanes` (works on `pos1 > pos2`). `origOn = showOrig && orig.length > 0`. Draw with the same pentagon path as `drawLanes` but fill `col + '33'`, stroke `cssVar('--mp-text-muted')` 1 px, or `cssVar('--mp-warning', '#8a5a00')` 2 px when `status !== 'unchanged'`; no partial dashing; label in muted colour. Tooltip text for `h.orig`: `'Original ' + name + ' -> ' + (gene || 'unmapped') + ' ' + dir + ' ' + pos1 + '-' + pos2 + ' | ' + STATUS_TEXT[status]`, appending `' (file: ' + orig1 + '-' + orig2 + ')'` when they differ from the mapped coords. `STATUS_TEXT = {unchanged:'unchanged', moved:'boundaries changed by curation', dropped:'dropped by curation', unmapped:'unmapped name, dropped on import', clipped:'clipped by trim'}`. Expose `STATUS_TEXT` on `window.mpseq.ui` so synteny reuses it.
- [ ] **Step 4: Run `test-seqview-js.R`, expect PASS** (skips cleanly without chromote; if so, run once with chromote available before committing).
- [ ] **Step 5: Commit** `feat: original user-annotation track in the sequence viewer`

---

### Task 6: Synteny viewer draws the original track

**Files:**
- Modify: `inst/app/www/synteny.js` (ctor ~95-104, `load` 113-139, `ys` 184-193, `draw` 197-222, `drawLanes` 251-281 gains a style branch for `track === 'orig'`, `tipText` 321-333, `click` 334-343, `bindControls` gains the checkbox loop, `state` 418-423)
- Modify: `tests/testthat/synteny/index.html` (add `<input type="checkbox" id="sy-show_orig" checked>` to the controls)
- Test: `tests/testthat/test-synteny-js.R`

**Interfaces:**
- Consumes: `p.original`, `window.mpseq.ui.STATUS_TEXT` (Task 5).
- Produces: `v.oFeats` with `track: 'orig'`, `c1/c2` via `G.sampleCol` like sample features; `v.oLanes`; `ys()` gains `oLanes` before `sLanes` with gutter label `user`; `window.mpsyn.state(id)` gains `origLanes`; clicking an orig arrow does nothing.

- [ ] **Step 1: Write failing chromote test** in the style of the existing "clicking a sample arrow" test: load the payload used there plus `original:[{name:'ND1',gene:'nad1',type:'PCG',pos1:<sample nad1 pos1>,pos2:<pos2+3>,dir:'+',status:'moved',orig1:..,orig2:..,odir:'+'}]`; assert `state.origLanes == 1`, `state.sampleLanes` unchanged, clicking the orig arrow leaves `selected` null and sends nothing, unchecking `sy-show_orig` sets `origLanes` to 0 and lowers `height`, and `tip` over the orig arrow contains `'boundaries changed by curation'`.
- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement** mirroring Task 5. The checkbox id is `this.id.replace(/synteny_canvas$/, 'show_orig')` (the UI uses `ns("show_orig")`; mirror in the fixture as `sy-show_orig`).
- [ ] **Step 4: Run `test-synteny-js.R`, expect PASS.**
- [ ] **Step 5: Commit** `feat: original user-annotation track in the synteny viewer`

---

### Task 7: Wire the track into the Annotation Details window

**Files:**
- Modify: `R/app_annotate_seqview.R` (`seqview_ui(id)` -> `seqview_ui(id, user_annot = FALSE)`; `seqview_server(id, rv, tick, selected)` -> `seqview_server(id, rv, tick, selected, original = reactive(NULL))`)
- Modify: `R/app_annotate_synteny.R` (`synteny_viewer_ui(ns)` -> `synteny_viewer_ui(ns, user_annot = FALSE)`)
- Modify: `R/app_annotate_details.R` (seqview wiring ~1971-1976; synteny payload observer ~1989-2003; `synteny_viewer_ui(ns)` call ~1112; `seqview_ui(ns("seqview"))` call ~4375)
- Test: `tests/testthat/test-user-annot-track.R` (one `shiny::testServer` test on `seqview_server` is optional; the pure pieces are tested already). Manual check below is required.

**Interfaces:**
- Consumes: `original_features()`, `flag_original()` (Task 3), payload `original =` (Task 4).
- Produces: reactive `orig_track()` in the details module: `NULL` unless `isTRUE(session$userData$user_annot)`; otherwise `list(ok, reason, rows)` where `rows` already carry `status`. Both UIs render `mp_checkbox(ns("show_orig"), "Original", value = TRUE)` only when `user_annot`, plus `uiOutput(ns("orig_note"))` under the controls. `seqview_server` renders `output$orig_note` (the `reason` as a `mp-coverage-caption` div, or NULL); the details module renders `output$synteny_orig_note` the same way inside the synteny `renderUI`.

- [ ] **Step 1: Add the reactive** in `annotations_details_server` next to `seqview_sel`: depends on `asmb_edit_tick()`, `gargoyle::watch("annotations_modal")`, `rv$annotations`, `rv$updating`; wraps `original_features()` in `tryCatch` returning `list(ok = FALSE, reason = "Original annotation cannot be placed on the edited sequence")` on error.
- [ ] **Step 2: Pass it through.** `seqview_server("seqview", rv, asmb_edit_tick, seqview_sel, orig_track)`; inside, `original = if (isTRUE(o$ok)) o$rows` in the `seqview_payload()` call and `output$orig_note`. Same `original =` in the `synteny_payload()` call.
- [ ] **Step 3: UI flags.** `seqview_ui(ns("seqview"), user_annot = isTRUE(session$userData$user_annot))` in `annotate_details_modal()`; `synteny_viewer_ui(ns, isTRUE(session$userData$user_annot))` in the synteny `renderUI`.
- [ ] **Step 4: Manual check** on a fresh demo project (`new_test_project_userAnnot(no_raw_data = TRUE)` into a scratch dir, run WF2 for GB_DROP and GB_MITOS, open the app):
  - Both viewers show a muted `user` track with raw names above the curated lanes; `Original` checkbox hides it.
  - Hover an amber-outlined arrow: tooltip shows the status text and the file coordinates.
  - Linearize GB_MITOS, then Trim: the original track follows (cob still over its curated row); Restore: still aligned.
  - A standard (non user-annot) project shows no checkbox and no track.
  Record what was seen in the task's commit message body (one line per bullet).
- [ ] **Step 5: Run the full test suite:** `Rscript -e 'devtools::test()'`, expect no new failures.
- [ ] **Step 6: Commit** `feat: show original user annotations above curated tracks`

---

### Task 8: NEWS and vignette

**Files:**
- Modify: `NEWS.md` (under `# MitoPilot (development version)`)
- Modify: `vignettes/Your-Own-Annotations.Rmd` (short paragraph in the section that describes the Annotation Details window; screenshot optional, do not block)

- [ ] **Step 1: NEWS bullet.** One line: the Sequence and BLAST Reference Synteny viewers show the original imported annotation as a muted track above the curated one, with curation edits outlined; toggle with the `Original` checkbox.
- [ ] **Step 2: Vignette paragraph** saying what the track shows, what the amber outline means (coordinates or strand changed, feature dropped, or name unmapped; name normalisation alone is not marked), and that it follows Linearize and Trim except where Trim cuts through a feature.
- [ ] **Step 3: Commit** `docs: original annotation track`
