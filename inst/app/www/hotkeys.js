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

  window.mpHotkeys = {
    enabled: enabled, setEnabled: setEnabled, keyOf: keyOf, context: context,
    bindings: bindings, bind: bind, click: click, clickIn: clickIn, handle: handle, help: help
  };
})();
