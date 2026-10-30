// Opt-in keyboard shortcuts. One document keydown listener resolves the
// context (open modal kind, else the current Step), looks the key up in that
// context's binding table, and acts by clicking the button the mouse would
// click, so every server-side guard stays in force. Off until the pill at the
// bottom right is switched on; remembered in localStorage.
(function () {
  var STORE = 'mp_hotkeys';
  var IGNORE = 'input, textarea, select, [contenteditable], .selectize-input, .bootstrap-select';
  var MOUSE_ONLY = [
    'Continue export / Stop without exporting', 'Run pipeline start / stop',
    'Build joined assembly', 'Linearize, trim ends, restore assembly',
    'Circularize details and specimen tabs'
  ];
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
  function context() {
    var m = document.getElementById('shiny-modal');
    if (has(m, '#annotate-annotations-close')) return 'details';
    if (has(m, '#assemble-coverage_details-close_modal')) return 'coverage';
    if (has(m, '#export-review_next')) return 'review';
    if (has(m, '[id$="-update_state"]')) return 'state';
    if (has(m, '.modal')) return 'modal';
    var r = document.querySelector('#mode input[type=radio]:checked');
    return 'page:' + (r ? r.value : 'Assemble');
  }
  function bind(ctx, key, label, run) {
    (bindings[ctx] || (bindings[ctx] = {}))[key] = { label: label, run: run };
  }
  function usable(el) { return !!el && el.offsetParent !== null && !el.disabled; }
  function click(id) {
    var el = document.getElementById(id);
    if (!usable(el)) return false;
    el.click(); return true;
  }
  function clickIn(root, sel) {
    var el = root ? root.querySelector(sel) : null;
    if (!usable(el)) return false;
    el.click(); return true;
  }

  // ---- cheat sheet ----
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
    var rows = Object.keys(b).map(function (k) {
      return '<tr><td><kbd>' + k.replace('Shift+', 'Shift+') + '</kbd></td><td>' + b[k].label + '</td></tr>';
    }).join('');
    helpEl.innerHTML = '<div class="mp-hotkeys-help-box"><h4>Keyboard shortcuts: ' + ctx.replace('page:', '') +
      '</h4><table>' + rows + '</table><p class="text-muted">Mouse only: ' + MOUSE_ONLY.join('; ') +
      '.</p><p class="text-muted">Press ? or Escape to close.</p></div>';
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
    if (e.target && e.target.closest && e.target.closest(IGNORE)) return;
    if (document.querySelector('.swal2-container')) return;
    if (k === '?') { help(true); e.preventDefault(); return; }
    var b = bindings[context()], hit = b && b[k];
    if (hit && hit.run()) e.preventDefault();
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
    if (t && t.classList && t.classList.contains('mp-seqview-canvas')) lastCanvas = t.id;
  });
  var DET = 'annotate-annotations-';
  function viewer() {
    var syn = lastCanvas && /synteny_canvas$/.test(lastCanvas) && document.getElementById(lastCanvas);
    var id = syn ? lastCanvas : (document.getElementById(lastCanvas || '') ? lastCanvas : DET + 'seqview-canvas');
    return { api: syn ? window.mpsyn : window.mpseq, id: id, attr: syn ? 'data-mpsyn' : 'data-mpseq' };
  }
  function viewerButton(action) {
    var v = viewer(), cv = document.getElementById(v.id);
    var sec = cv && cv.closest('details');
    openPanel(sec);
    return clickIn(sec || document.getElementById('shiny-modal'), '[' + v.attr + '="' + action + '"]');
  }
  function zoom(f) { var v = viewer(); if (v.api && v.api.zoom) { v.api.zoom(v.id, f); return true; } return false; }
  function visible(id) { var el = document.getElementById(id); return !!el && el.offsetParent !== null; }
  function focusId(id) { var el = document.getElementById(id); if (!el) return false; el.focus(); return true; }

  // ---- Annotate details modal ----
  (function () {
    var D = DET, ctx = 'details';
    function b(key, label, run) { bind(ctx, key, label, run); }
    function c(key, label, id) { b(key, label, function () { return click(D + id); }); }
    function notesDetails() { var n = document.getElementById(D + 'notes'); return n && n.closest('details'); }
    function selectedRow() { return document.querySelector('#' + D + 'table .rt-tr.rt-tr-selected'); }

    b('Escape', 'Close (cancels an open merge panel first)', function () {
      return visible(D + 'merge_select_div') ? click(D + 'cancel_merge') : click(D + 'close');
    });
    c('q', 'Close', 'close');
    c('Q', 'Lock and close', 'lock');
    c('L', 'Lock, stay open', 'lock_only');
    c('Backspace', 'Back to Review (outlier review only)', 'back_to_review');
    c('Shift+ArrowRight', 'Next sample', 'nav_next');
    c('Shift+ArrowLeft', 'Previous sample', 'nav_prev');

    [['t', 'annotation_table_details', 'Annotation Table'], ['s', 'seqview-section', 'Sequence'],
     ['b', 'blast_synteny_details', 'BLAST Reference Synteny'], ['a', 'alignment_div', 'Alignment']
    ].forEach(function (p) {
      b(p[0], 'Toggle ' + p[2] + ' panel', function () { return panel(D + p[1], false); });
      b(p[0].toUpperCase(), 'Open and scroll to ' + p[2], function () { return panel(D + p[1], true); });
    });
    b('i', 'Toggle Notes panel', function () { return panel(notesDetails(), false); });
    b('I', 'Open Notes and start typing', function () {
      if (!panel(notesDetails(), true)) return false;
      return focusId(D + 'notes');
    });

    c('v', 'Toggle ID verified', 'ID_verified');
    c('r', 'Toggle Reviewed', 'reviewed');
    c('p', 'Toggle Problematic', 'problematic');
    c('l', 'Toggle Partial', 'partial');

    c('d', 'Delete annotation', 'delete');
    c('D', 'Restore deleted annotation', 'restore');
    c('m', 'Merge', 'merge');
    b('Enter', 'Confirm merge (merge panel open)', function () {
      return visible(D + 'merge_select_div') && click(D + 'confirm_merge');
    });
    c('u', 'Un-join', 'unjoin');
    c('g', 'Assign gene name', 'assign_gene');
    c('o', 'Auto-assign ORFs', 'auto_assign_orfs');
    b('y', 'Copy nucleotide FASTA', function () { return clickIn(selectedRow(), 'button[onclick*="copy_fas"]'); });
    b('Y', 'Copy amino-acid FASTA', function () { return clickIn(selectedRow(), 'button[onclick*="copy_faa"]'); });

    b('+', 'Zoom in (hovered viewer)', function () { return zoom(2); });
    b('=', 'Zoom in (hovered viewer)', function () { return zoom(2); });
    b('-', 'Zoom out (hovered viewer)', function () { return zoom(0.5); });
    b('F', 'Fit selected gene', function () { return viewerButton('fit'); });
    b('W', 'Whole genome / alignment', function () { return viewerButton('whole'); });
    b('C', 'Copy visible sequence', function () { return viewerButton('copy'); });
    b('G', 'Go to position', function () {
      openPanel(document.getElementById(D + 'seqview-section'));
      return focusId(D + 'seqview-goto');
    });
  })();

  window.mpHotkeys = {
    enabled: enabled, setEnabled: setEnabled, keyOf: keyOf, context: context,
    bindings: bindings, bind: bind, click: click, clickIn: clickIn, handle: handle, help: help,
    panel: panel, viewer: viewer, details: DET
  };
})();
