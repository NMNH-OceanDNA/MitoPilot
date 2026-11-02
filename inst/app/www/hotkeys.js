// Opt-in keyboard shortcuts. One document keydown listener resolves the
// context (open modal kind, else the current Step), looks the key up in that
// context's binding table, and acts by clicking the button the mouse would
// click, so every server-side guard stays in force. Off until the pill at the
// bottom right is switched on; remembered in localStorage.
(function () {
  var STORE = 'mp_hotkeys';
  var IGNORE = 'input:not([type=checkbox]):not([type=radio]):not([type=button]):not([type=submit]), ' +
    'textarea, select, [contenteditable], .selectize-input, .bootstrap-select';
  var TEXTY = 'input, textarea';
  var bindings = {};
  var helpEl = null;

  function enabled() { return localStorage.getItem(STORE) === '1'; }
  function setEnabled(on) {
    localStorage.setItem(STORE, on ? '1' : '0');
    var box = document.getElementById('mp-hotkeys-on');
    if (box) box.checked = !!on;
  }
  function keyOf(e) {
    if (e.key === ' ') return 'Space';
    if (e.key.length === 1) return e.key;
    return (e.shiftKey ? 'Shift+' : '') + e.key;
  }
  function has(root, sel) { return !!(root && root.querySelector(sel)); }
  // Module ids carry a runtime-nested prefix (e.g. annotate-annotate-annotations-),
  // so take the prefix from a marker element's id instead of hardcoding it.
  // mark.match identifies the element; mark.ctrl is the un-namespaced input id.
  function nsFrom(mark) {
    var el = document.querySelector('#shiny-modal [id$="' + mark.match + '"]');
    return el ? el.id.slice(0, el.id.length - mark.ctrl.length) : null;
  }
  var DETAILS_MARK = { match: 'annotations-close', ctrl: 'close' };
  var COVERAGE_MARK = { match: 'coverage_details-close_modal', ctrl: 'close_modal' };
  function context() {
    var m = document.getElementById('shiny-modal');
    if (nsFrom(DETAILS_MARK)) return 'details';
    if (nsFrom(COVERAGE_MARK)) return 'coverage';
    if (has(m, '[id$="-review_next"]')) return 'review';
    if (has(m, '[id$="-update_state"]')) return 'state';
    if (m) return 'modal';
    var r = document.querySelector('#mode input[type=radio]:checked');
    return 'page:' + (r ? r.value : 'Assemble');
  }
  // opts.show: text for the cheat sheet; opts.alias: hide this key from it.
  function bind(ctx, key, label, run, opts) {
    (bindings[ctx] || (bindings[ctx] = {}))[key] = { label: label, run: run, show: opts && opts.show, alias: !!(opts && opts.alias) };
  }
  function usable(el) { return !!el && el.offsetParent !== null && !el.disabled; }
  // Toolbar ids repeat across the Step panels; click the copy that is shown.
  function click(id) {
    var els = document.querySelectorAll('[id="' + id + '"]');
    for (var i = 0; i < els.length; i++) {
      if (usable(els[i])) { els[i].click(); return true; }
    }
    return false;
  }
  function clickIn(root, sel) {
    var el = root ? root.querySelector(sel) : null;
    if (!usable(el)) return false;
    el.click(); return true;
  }

  // ---- cheat sheet ----
  // A bare capital letter is typed with Shift; show it that way.
  function keyLabel(k) { return /^[A-Z]$/.test(k) ? 'Shift+' + k.toLowerCase() : k; }
  function help(show) {
    if (!helpEl) {
      helpEl = document.createElement('div');
      helpEl.id = 'mp-hotkeys-help'; helpEl.hidden = true;
      helpEl.addEventListener('click', function () { help(false); });
      document.body.appendChild(helpEl);
    }
    if (show === undefined) show = helpEl.hidden;
    if (!show) { helpEl.hidden = true; return; }
    var ctx = context(), b = bindings[ctx] || {};
    var rows = Object.keys(b).filter(function (k) { return !b[k].alias; }).map(function (k) {
      return '<tr><td><kbd>' + (b[k].show || keyLabel(k)) + '</kbd></td><td>' + b[k].label + '</td></tr>';
    }).join('');
    helpEl.innerHTML = '<div class="mp-hotkeys-help-box"><h4>Keyboard shortcuts: ' + ctx.replace('page:', '') +
      '</h4><table>' + rows + '</table><p class="text-muted">Press ? or Escape to close.</p></div>';
    helpEl.hidden = false;
  }

  // ---- listener ----
  function handle(e) {
    if (!enabled()) return;
    var k = keyOf(e);
    if (helpEl && !helpEl.hidden) {
      if (k === 'Escape' || k === '?') { help(false); e.preventDefault(); }
      return;
    }
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    var t = e.target && e.target.closest ? e.target : null;
    if (t && t.closest(IGNORE)) {
      if (k === 'Escape' && t.matches(TEXTY)) { t.blur(); e.preventDefault(); }
      return;
    }
    if (t && k === 'Space' && t.matches('input[type=checkbox], input[type=radio]')) return;
    if (document.querySelector('.swal2-container')) return;
    if (k === '?') { help(true); e.preventDefault(); return; }
    var b = bindings[context()], hit = b && b[k];
    if (hit && hit.run(e)) e.preventDefault();
  }
  document.addEventListener('keydown', handle, true);

  // ---- pill ----
  function pill() {
    if (document.getElementById('mp-hotkeys-pill')) return;
    var d = document.createElement('div');
    d.id = 'mp-hotkeys-pill';
    d.innerHTML = '<label><input type="checkbox" id="mp-hotkeys-on"> Hotkeys</label>' +
      '<button type="button" id="mp-hotkeys-help-btn" aria-label="Keyboard shortcuts">?</button>';
    document.body.appendChild(d);
    var box = d.querySelector('#mp-hotkeys-on');
    box.checked = enabled();
    box.addEventListener('change', function () { setEnabled(box.checked); });
    d.querySelector('#mp-hotkeys-help-btn').addEventListener('click', function () { help(true); });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', pill); else pill();


  // ---- reactable rows ----
  // Rows of a reactable body, minus pagination padding and rows hidden by CSS
  // (the lock/state pickers hide rows with display:none, still in the DOM).
  function rows(tableSel) {
    var t = document.querySelector(tableSel);
    if (!t) return [];
    return Array.prototype.filter.call(t.querySelectorAll('.rt-tbody .rt-tr:not(.rt-tr-pad)'), function (r) {
      return getComputedStyle(r).display !== 'none';
    });
  }
  // Single-select tables hide their selection column; a cell click selects.
  // scroll=false leaves the window where it is (the details window keeps the
  // user's place among the panels while the selected annotation changes).
  function selectRow(row, scroll) {
    if (!row) return false;
    var el = row.querySelector('input.rt-select-input') || row.querySelector('.rt-td');
    if (!el) return false;
    el.click();
    if (scroll !== false) row.scrollIntoView({ block: 'nearest' });
    return true;
  }
  function stepRow(tableSel, by, fromRow) {
    var rs = rows(tableSel);
    if (!rs.length) return null;
    var i = rs.indexOf(fromRow);
    if (i < 0) return by > 0 ? rs[0] : rs[rs.length - 1];
    var j = Math.max(0, Math.min(rs.length - 1, i + by));
    return j === i ? null : rs[j];
  }
  function post(id, value) {
    if (!window.Shiny) return false;
    window.Shiny.setInputValue(id, value, { priority: 'event' });
    return true;
  }

  // ---- shared helpers for modal bindings ----
  function summaryClick(el) {
    var sm = el && el.querySelector('summary');
    if (sm) sm.click();
  }
  // Toggle a <details> (or open + scroll to it with shift). Always via the
  // summary so an onclick on it (the Alignment panel's align trigger) fires.
  function panel(el, scroll) {
    if (typeof el === 'string') el = document.getElementById(el);
    if (!el) return false;
    if (!scroll) { summaryClick(el); return true; }
    if (!el.open) summaryClick(el);
    el.scrollIntoView({ block: 'start', behavior: 'smooth' });
    return true;
  }
  function openPanel(el) { if (el && !el.open) summaryClick(el); return el; }
  // Viewer keys act on the canvas the mouse last hovered; default sequence.
  var lastCanvas = null;
  document.addEventListener('mouseover', function (e) {
    var t = e.target;
    if (!t || !t.closest) return;
    if (t.classList.contains('mp-seqview-canvas')) lastCanvas = t.id.replace(/-reads$/, '-canvas');
    else if (t.closest('.msaR')) lastCanvas = 'msa';
  });
  function DET() { return nsFrom(DETAILS_MARK) || ''; }
  var msaApi = {
    pan: function (id, frac) {
      var el = document.querySelector('#shiny-modal .biojs_msa_rheader');
      if (!el) return;
      el.scrollLeft += frac * el.clientWidth;
      el.dispatchEvent(new Event('scroll'));
    }
  };
  function viewer() {
    if (lastCanvas === 'msa' && document.querySelector('#shiny-modal .msaR')) return { api: msaApi, id: 'msa', attr: null };
    var id = lastCanvas && document.getElementById(lastCanvas) ? lastCanvas : null;
    if (!id) {
      var first = document.querySelector('#shiny-modal .mp-seqview-canvas');
      id = first ? first.id : DET() + 'seqview-canvas';
    }
    var syn = /synteny_canvas$/.test(id);
    return { api: syn ? window.mpsyn : window.mpseq, id: id, attr: syn ? 'data-mpsyn' : 'data-mpseq' };
  }
  function viewerButton(action) {
    var v = viewer(), cv = document.getElementById(v.id);
    if (!v.attr) return false;
    var sec = cv && cv.closest('details');
    openPanel(sec);
    return clickIn(sec || document.getElementById('shiny-modal'), '[' + v.attr + '="' + action + '"]');
  }
  function zoom(f) { var v = viewer(); if (v.api && v.api.zoom) { v.api.zoom(v.id, f); return true; } return false; }
  // + or = zooms in, - or _ zooms out; the sheet lists each pair once.
  function zoomKeys(ctx, what) {
    bind(ctx, '+', 'Zoom in (' + what + ')', function () { return zoom(2); }, { show: '+ or =' });
    bind(ctx, '=', 'Zoom in (' + what + ')', function () { return zoom(2); }, { alias: true });
    bind(ctx, '-', 'Zoom out (' + what + ')', function () { return zoom(0.5); }, { show: '- or _' });
    bind(ctx, '_', 'Zoom out (' + what + ')', function () { return zoom(0.5); }, { alias: true });
  }
  function pan(f) { var v = viewer(); if (v.api && v.api.pan) { v.api.pan(v.id, f); return true; } return false; }
  function visible(id) { var el = document.getElementById(id); return !!el && el.offsetParent !== null; }
  function focusId(id) { var el = document.getElementById(id); if (!el) return false; el.focus(); return true; }

  // ---- Annotate details modal ----
  (function () {
    var D = DET, ctx = 'details';
    function b(key, label, run) { bind(ctx, key, label, run); }
    function c(key, label, id) { b(key, label, function () { return click(D() + id); }); }
    function notesDetails() { var n = document.getElementById(D() + 'notes'); return n && n.closest('details'); }
    function selectedRow() { return document.querySelector('#' + D() + 'table .rt-tr.rt-tr-selected'); }

    b('Escape', 'Close (cancels an open merge panel first)', function () {
      return visible(D() + 'merge_select_div') ? click(D() + 'cancel_merge') : click(D() + 'close');
    });
    c('q', 'Close', 'close');
    c('Q', 'Lock and close', 'lock');
    b('L', 'Lock / Unlock, stay open', function () { return click(D() + 'lock_only') || click(D() + 'unlock'); });
    c('Backspace', 'Back to Review (outlier review only)', 'back_to_review');
    c('Shift+ArrowRight', 'Next sample', 'nav_next');
    c('Shift+ArrowLeft', 'Previous sample', 'nav_prev');

    [['t', 'annotation_table_details', 'Annotation Table'], ['s', 'seqview-section', 'Sequence'],
     ['b', 'blast_synteny_details', 'BLAST Reference Synteny'], ['a', 'alignment_div', 'Alignment']
    ].forEach(function (p) {
      b(p[0], 'Toggle ' + p[2] + ' panel', function () { return panel(D() + p[1], false); });
      b(p[0].toUpperCase(), 'Open and scroll to ' + p[2], function () { return panel(D() + p[1], true); });
    });
    b('i', 'Toggle Notes panel', function () { return panel(notesDetails(), false); });
    b('I', 'Open Notes and start typing', function () {
      if (!panel(notesDetails(), true)) return false;
      return focusId(D() + 'notes');
    });

    c('v', 'Toggle ID verified', 'ID_verified');
    c('r', 'Toggle Reviewed', 'reviewed');
    c('p', 'Toggle Problematic', 'problematic');
    c('l', 'Toggle Partial', 'partial');

    c('d', 'Delete annotation', 'delete');
    c('D', 'Restore deleted annotation', 'restore');
    c('m', 'Merge', 'merge');
    b('Enter', 'Confirm merge (merge panel open)', function () {
      return visible(D() + 'merge_select_div') && click(D() + 'confirm_merge');
    });
    c('u', 'Un-join', 'unjoin');
    c('g', 'Assign gene name', 'assign_gene');
    c('o', 'Auto-assign ORFs', 'auto_assign_orfs');
    b('y', 'Copy nucleotide FASTA', function () { return clickIn(selectedRow(), 'button[onclick*="copy_fas"]'); });
    b('Y', 'Copy amino-acid FASTA', function () { return clickIn(selectedRow(), 'button[onclick*="copy_faa"]'); });

    zoomKeys(ctx, 'hovered viewer');
    b('ArrowLeft', 'Pan left (hovered viewer or alignment)', function () { return pan(-0.25); });
    b('ArrowRight', 'Pan right (hovered viewer or alignment)', function () { return pan(0.25); });
    b('F', 'Fit selected gene', function () { return viewerButton('fit'); });
    b('W', 'Whole genome / alignment', function () { return viewerButton('whole'); });
    b('C', 'Copy visible sequence', function () { return viewerButton('copy'); });

    // Annotation rows: single-select table, 50-row pages; walk DOM neighbours.
    function T() { return '#' + D() + 'table'; }
    function selRow() { return document.querySelector(T() + ' .rt-tr.rt-tr-selected'); }
    // Plain arrows scroll the window; Shift+arrows step through annotations.
    function scrollBox() {
      var cands = document.querySelectorAll('#shiny-modal .modal-body, #shiny-modal');
      for (var i = 0; i < cands.length; i++) if (cands[i].scrollHeight > cands[i].clientHeight) return cands[i];
      return null;
    }
    function scrollWin(dy) { var el = scrollBox(); if (!el) return false; el.scrollTop += dy; return true; }
    b('ArrowDown', 'Scroll down', function () { return scrollWin(80); });
    b('ArrowUp', 'Scroll up', function () { return scrollWin(-80); });
    b('Shift+ArrowDown', 'Next annotation', function () { return selectRow(stepRow(T(), 1, selRow()), false); });
    b('Shift+ArrowUp', 'Previous annotation', function () { return selectRow(stepRow(T(), -1, selRow()), false); });
    b('Home', 'First annotation', function () { var r = rows(T()); return selectRow(r[0], false); });
    b('End', 'Last annotation', function () { var r = rows(T()); return selectRow(r[r.length - 1], false); });

    // Edit mode. Boundary keys move by exactly one codon or base: they post the
    // value "one", which the R handlers read as step size 1. Codon buttons when
    // the segment owns that end (wrapper visible), else the nucleotide nudge.
    function editing() { return visible(D() + 'edit_mode_ctrls'); }
    function boundary(end, dir, step) {
      if (!editing()) return false;
      var codon = visible(D() + (end === '5' ? 'start_search_ctrl' : 'stop_search_ctrl'));
      var id = codon
        ? (end === '5' ? (dir === 'out' ? 'start-add' : 'start-minus') : (dir === 'out' ? 'stop-add' : 'stop-minus'))
        : 'rrna-' + end + '-' + dir;
      if (!codon && !document.getElementById(D() + id)) return false;
      return post(D() + id, step || 'one');
    }
    c('e', 'Edit mode', 'edit_mode');
    c('E', 'Save edits', 'save_edits');
    c('z', 'Reset edits', 'discard_edits');
    b(',', "5' boundary out by one", function () { return boundary('5', 'out'); });
    b('.', "5' boundary in by one", function () { return boundary('5', 'in'); });
    b('<', "3' boundary in by one", function () { return boundary('3', 'in'); });
    b('>', "3' boundary out by one", function () { return boundary('3', 'out'); });
    b('[', "5' boundary out by the step box", function () { return boundary('5', 'out', 'box'); });
    b(']', "5' boundary in by the step box", function () { return boundary('5', 'in', 'box'); });
    b('{', "3' boundary in by the step box", function () { return boundary('3', 'in', 'box'); });
    b('}', "3' boundary out by the step box", function () { return boundary('3', 'out', 'box'); });
    b('c', 'Toggle single-codon stepping', function () {
      return editing() && clickIn(document.getElementById(D() + 'edit_mode_ctrls'), 'input[onchange*="single_codon"]');
    });
    c('5', "Toggle partial 5'", 'toggle_partial_start');
    c('3', "Toggle partial 3'", 'toggle_partial_stop');
    c('P', 'poly-A stop', 'polyA_stop');
    c('f', 'Toggle Align fewer refs', 'reduce_align');

    b('G', 'Go to position', function () {
      openPanel(document.getElementById(D() + 'seqview-section'));
      return focusId(D() + 'seqview-goto');
    });
  })();


  // ---- keyboard cursor on reactable tables (keyed by table selector) ----
  var cursors = {}, lastToggled = {};
  function cursorRow(sel) {
    var r = cursors[sel];
    if (r && document.body.contains(r)) return r;
    delete cursors[sel]; return null;
  }
  function setCursor(sel, row) {
    var cur = cursorRow(sel);
    if (cur) cur.removeAttribute('data-mp-kb-cursor');
    if (!row) { delete cursors[sel]; return false; }
    row.setAttribute('data-mp-kb-cursor', '');
    row.scrollIntoView({ block: 'nearest' });
    cursors[sel] = row; return true;
  }
  function moveCursor(sel, by) {
    var next = stepRow(sel, by, cursorRow(sel));
    return next ? setCursor(sel, next) : !!cursorRow(sel);
  }
  function toggleCursor(sel) {
    var r = cursorRow(sel); if (!r) return false;
    lastToggled[sel] = r;
    return clickIn(r, 'input.rt-select-input');
  }
  // Check every row between the last Space-toggled row and the cursor.
  function extendCursor(sel) {
    var r = cursorRow(sel), rs = rows(sel), a = rs.indexOf(lastToggled[sel]), z = rs.indexOf(r);
    if (a < 0 || z < 0) return false;
    for (var i = Math.min(a, z); i <= Math.max(a, z); i++) {
      var box = rs[i].querySelector('input.rt-select-input');
      if (box && !box.checked) box.click();
    }
    return true;
  }
  function tableId(sel) { return sel.replace(/^#/, ''); }
  function reactable(fn) { if (!window.Reactable) return false; fn(window.Reactable); return true; }

  // ---- top-level pages: Assemble / Annotate / Export ----
  ['Assemble', 'Annotate', 'Export'].forEach(function (step) {
    var ctx = 'page:' + step, ns = step.toLowerCase(), T = '#' + ns + '-table';
    function b(key, label, run) { bind(ctx, key, label, run); }
    function c(key, label, id) { b(key, label, function () { return click(id); }); }
    ['Assemble', 'Annotate', 'Export'].forEach(function (st, i) {
      b(String(i + 1), 'Go to Step ' + st, function () {
        var r = document.querySelector('#mode input[value="' + st + '"]');
        if (!r) return false;
        (r.closest('.btn') || r).click(); return true;
      });
    });
    c('r', 'Refresh table', 'refresh');
    b('/', 'Search the table', function () { return focusIn(T + ' .rt-search'); });
    b('ArrowDown', 'Move row cursor down', function (e) { return moveCursor(T, e && e.repeat ? 3 : 1); });
    b('ArrowUp', 'Move row cursor up', function (e) { return moveCursor(T, e && e.repeat ? -3 : -1); });
    b('Space', 'Check / uncheck cursor row', function () { return toggleCursor(T); });
    b('x', 'Check / uncheck cursor row', function () { return toggleCursor(T); });
    b('X', 'Check rows from last toggled to cursor', function () { return extendCursor(T); });
    b('a', 'Check / uncheck all rows on the page', function () {
      return reactable(function (R) { R.toggleAllRowsSelected(tableId(T)); });
    });
    b('Escape', 'Clear row cursor', function () { return setCursor(T, null) || true; });
    b('Enter', 'Open details for cursor row', function () { return clickIn(cursorRow(T), 'button[onclick*="-details"]'); });
    b('o', 'Open details for cursor row', function () { return clickIn(cursorRow(T), 'button[onclick*="-details"]'); });
    b('[', 'Previous page', function () {
      return reactable(function (R) { R.gotoPage(tableId(T), R.getState(tableId(T)).page - 1); });
    });
    b(']', 'Next page', function () {
      return reactable(function (R) { R.gotoPage(tableId(T), R.getState(tableId(T)).page + 1); });
    });
    c('s', 'Set State (checked rows)', 'state');
    c('L', 'Lock / Unlock (checked rows)', 'lock');
    c('u', 'Update (run pipeline)', 'run_modal');
    c('m', 'Metadata fields', ns + '-meta_view_open');
    c('w', 'Work dirs', 'workdir_browser-open_browser');
    c('R', 'Run reports', 'run_reports-open');
    b('O', 'Open output folder of cursor row', function () { return clickIn(cursorRow(T), '[onclick*="-output"]'); });
    b('S', 'Specimen metadata of cursor row', function () { return clickIn(cursorRow(T), '[onclick*="specimen_open"]'); });
    if (step === 'Annotate') {
      c('v', 'Mark ID Verified (checked rows)', 'id_verified_top');
      c('p', 'Mark Problematic (checked rows)', 'problematic_top');
      c('l', 'Mark Partial (checked rows)', 'partial_top');
    }
    if (step === 'Export') {
      c('g', 'Assign Group (checked rows)', 'group');
      c('G', 'Clear Group (checked rows)', 'clear_group');
      c('e', 'Export Data', 'export');
    }
  });
  function focusIn(sel) { var el = document.querySelector(sel); if (!el) return false; el.focus(); return true; }

  // ---- Set State modal ----
  (function () {
    for (var n = 1; n <= 9; n++) (function (n) {
      bind('state', String(n), 'Pick state ' + n, function () {
        var rs = document.querySelectorAll('[id$="-new_state"] input[type=radio]');
        if (!rs[n - 1]) return false;
        rs[n - 1].click(); return true;
      });
    })(n);
    bind('state', 'Enter', 'Set state', function () {
      var el = document.querySelector('[id$="-update_state"]');
      if (!usable(el)) return false;
      el.click(); return true;
    });
  })();


  // ---- Assemble details modal ----
  (function () {
    function C() { return nsFrom(COVERAGE_MARK) || ''; }
    var ctx = 'coverage';
    function T() { return '#' + C() + 'table'; }
    function b(key, label, run) { bind(ctx, key, label, run); }
    function c(key, label, id) { b(key, label, function () { return click(C() + id); }); }
    function inModal(sel) { return clickIn(document.getElementById('shiny-modal'), sel); }
    c('Escape', 'Close', 'close_modal');
    c('q', 'Close', 'close_modal');
    b('L', 'Lock / Unlock', function () { return click(C() + 'lock') || click(C() + 'unlock'); });
    c('a', 'Align paths', 'align');
    c(']', 'Next conflict block', 'next_block');
    c('n', 'Next conflict block', 'next_block');
    c('[', 'Previous conflict block', 'prev_block');
    c('p', 'Previous conflict block', 'prev_block');
    [['path', "Use one path's bases"], ['majority', 'Coverage-majority base'],
     ['iupac', 'IUPAC ambiguity'], ['nmask', 'Mask with N']].forEach(function (mode, i) {
      b(String(i + 1), 'Resolve mode: ' + mode[1], function () {
        var r = document.querySelector('#' + C() + 'res_mode input[type=radio][value="' + mode[0] + '"]');
        if (!r) return false;
        r.click(); return true;
      });
    });
    c('B', 'Build resolved assembly', 'build_resolved');
    c('T', 'Trim to consensus', 'trim_consensus');
    b('i', 'Focus notes', function () { return focusId(C() + 'notes'); });
    b('I', 'Scroll to and focus notes', function () {
      var el = document.getElementById(C() + 'notes'); if (!el) return false;
      el.scrollIntoView({ block: 'start', behavior: 'smooth' }); el.focus(); return true;
    });
    b('ArrowDown', 'Move scaffold cursor down', function (e) { return moveCursor(T(), e && e.repeat ? 3 : 1); });
    b('ArrowUp', 'Move scaffold cursor up', function (e) { return moveCursor(T(), e && e.repeat ? -3 : -1); });
    b('Space', 'Check / uncheck cursor scaffold', function () { return toggleCursor(T()); });
    b('X', 'Include / ignore cursor scaffold', function () { return clickIn(cursorRow(T()), '[onclick*="-ignore"]'); });
    zoomKeys(ctx, 'MapToRef viewer');
    b('ArrowLeft', 'Pan left', function () { return pan(-0.25); });
    b('ArrowRight', 'Pan right', function () { return pan(0.25); });
    b('W', 'Whole genome', function () { return viewerButton('whole'); });
    b('G', 'Go to position', function () { return focusIn('#shiny-modal [id$="-goto"]'); });
    b('c', 'Toggle coverage layer', function () { return inModal('[id$="-show_cov"]'); });
    b('t', 'Toggle nucleotide layer', function () { return inModal('[id$="-show_nt"]'); });
    b('R', 'Toggle reads layer', function () { return inModal('[id$="-show_reads"]'); });
  })();

  // ---- Export outlier review modal ----
  (function () {
    var E = 'export-', ctx = 'review', T = '#' + E + 'review_table';
    function b(key, label, run) { bind(ctx, key, label, run); }
    function c(key, label, id) { b(key, label, function () { return click(E + id); }); }
    c('ArrowLeft', 'Previous gene', 'review_prev');
    c('[', 'Previous gene', 'review_prev');
    c('ArrowRight', 'Next gene', 'review_next');
    c(']', 'Next gene', 'review_next');
    b('ArrowDown', 'Move sample cursor down', function (e) { return moveCursor(T, e && e.repeat ? 3 : 1); });
    b('ArrowUp', 'Move sample cursor up', function (e) { return moveCursor(T, e && e.repeat ? -3 : -1); });
    b('h', 'Highlight cursor sample in the alignment', function () { return clickIn(cursorRow(T), '[onclick*="review_pick"]'); });
    b('x', 'Toggle resolved for cursor sample', function () { return clickIn(cursorRow(T), '[onclick*="toggle_resolved"]'); });
    b('Space', 'Toggle resolved for cursor sample', function () { return clickIn(cursorRow(T), '[onclick*="toggle_resolved"]'); });
    b('Enter', 'Edit cursor sample in Annotate', function () { return clickIn(cursorRow(T), '[onclick*="goto_annot"]'); });
    b('e', 'Edit cursor sample in Annotate', function () { return clickIn(cursorRow(T), '[onclick*="goto_annot"]'); });
    c('S', 'Mark whole gene resolved', 'skip_gene');
  })();

  window.mpHotkeys = {
    enabled: enabled, setEnabled: setEnabled, keyOf: keyOf, context: context,
    bindings: bindings, bind: bind, click: click, clickIn: clickIn, handle: handle, help: help,
    panel: panel, viewer: viewer, nsFrom: nsFrom,
    rows: rows, selectRow: selectRow, stepRow: stepRow,
    moveCursor: moveCursor, cursorRow: cursorRow, setCursor: setCursor
  };
})();
