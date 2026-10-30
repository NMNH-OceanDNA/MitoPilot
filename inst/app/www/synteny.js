// inst/app/www/synteny.js
// Interactive BLAST reference synteny viewer (tools/synteny_viewer_spec.md).
// The x axis is pairwise-alignment columns; seqviewer.js must load first.
(function () {
  'use strict';

  // ---- geometry: pure functions over alignment columns, testable without a canvas ----
  function nongap(str) {
    var out = [];
    for (var i = 0; i < str.length; i++) if (str.charAt(i) !== '-') out.push(i + 1);
    return out;
  }
  // counts[c] = non-gap characters in columns 1..c
  function counts(str) {
    var out = new Int32Array(str.length + 1);
    for (var i = 0; i < str.length; i++) out[i + 1] = out[i] + (str.charAt(i) === '-' ? 0 : 1);
    return out;
  }
  // 0 match, 1 mismatch, 2 gap on either side
  function classify(s, r) {
    var out = new Uint8Array(s.length);
    for (var i = 0; i < s.length; i++) {
      var a = s.charAt(i), b = r.charAt(i);
      out[i] = (a === '-' || b === '-') ? 2 : a.toUpperCase() === b.toUpperCase() ? 0 : 1;
    }
    return out;
  }
  function clampIdx(p, n) { return Math.min(Math.max(p, 1), n); }
  function sampleCol(p, sNongap) { return sNongap[clampIdx(p, sNongap.length) - 1]; }
  function rotFlip(p, refLength, rotation, rc) {
    var pr = ((p - 1 - rotation) % refLength + refLength) % refLength + 1;
    return rc ? refLength - pr + 1 : pr;
  }
  function refCol(p, rNongap, refLength, rotation, rc) {
    return rNongap[clampIdx(rotFlip(p, refLength, rotation, rc), rNongap.length) - 1];
  }
  function pieces(c1, c2, L) { return c1 <= c2 ? [[c1, c2]] : [[c1, L], [1, c2]]; }
  // Drawn column span of a reference gene. On a flipped reference the drawn
  // order reverses, so pos2 maps to the start and the strand inverts; a gene
  // straddling the rotation point then comes out as c1 > c2 like any wrap.
  function refFeature(f, ctx) {
    var a = refCol(f.pos1, ctx.rNongap, ctx.refLength, ctx.rotation, ctx.rc);
    var b = refCol(f.pos2, ctx.rNongap, ctx.refLength, ctx.rotation, ctx.rc);
    if (ctx.rc) return { c1: b, c2: a, dir: f.dir === '-' ? '+' : '-' };
    return { c1: a, c2: b, dir: f.dir };
  }
  function sampleAt(c, s, sCounts) { return s.charAt(c - 1) === '-' ? null : sCounts[c]; }
  function refAt(c, r, rCounts, refLength, rotation, rc) {
    if (r.charAt(c - 1) === '-') return null;
    var k = rCounts[c], pr = rc ? refLength - k + 1 : k;
    return ((pr - 1 + rotation) % refLength + refLength) % refLength + 1;
  }
  function prefixBy(cls, keep) {
    var out = new Int32Array(cls.length + 1);
    for (var i = 0; i < cls.length; i++) out[i + 1] = out[i] + (keep(cls[i]) ? 1 : 0);
    return out;
  }
  function prefix(cls, value) { return prefixBy(cls, function (v) { return v === value; }); }
  function prefixNot(cls, value) { return prefixBy(cls, function (v) { return v !== value; }); }
  function identity(c0, c1, matchCounts, bothCounts) {
    var both = bothCounts[c1] - bothCounts[c0 - 1];
    return both > 0 ? (matchCounts[c1] - matchCounts[c0 - 1]) / both : null;
  }

  window.mpsyn = window.mpsyn || {};
  window.mpsyn.geom = { nongap: nongap, counts: counts, classify: classify, sampleCol: sampleCol,
                        refCol: refCol, pieces: pieces, refFeature: refFeature, sampleAt: sampleAt,
                        refAt: refAt, prefix: prefix, prefixNot: prefixNot, identity: identity };
})();

(function () {
  'use strict';
  var G = window.mpsyn.geom;
  var MAX_PPC = 14, LETTER = 8, BAR = 3;
  var RULER_H = 22, LANE_H = 22, NT_H = 20, IDENT_H = 40, GUTTER = 60, PAD = 4, GAP = 6;
  // Same base shades as seqviewer.js; class colours as the old zoom panel.
  var BASE = { A: '#4faf45', C: '#e0a53f', G: '#e0555a', T: '#4a90d9', N: '#666666' };
  var CLS = ['#60BD68', '#E55330', '#CCCCCC'], CLS_NAME = ['match', 'mismatch', 'gap'];
  var viewers = {};
  function ui() { return window.mpseq.ui; }
  function cssVar(n, f) { return ui().cssVar(n, f); }

  function Viewer(id) {
    this.id = id;
    this.canvas = document.getElementById(id);
    this.wrap = this.canvas.parentElement;
    this.tip = document.getElementById(id.replace(/canvas$/, 'tip'));
    this.section = this.canvas.closest('details');
    this.L = 0; this.version = null; this.selected = null; this.viewStart = 1; this.ppc = 1;
    this.sFeats = []; this.rFeats = []; this.sLanes = 0; this.rLanes = 0; this.hits = [];
    this.bindControls();
    var self = this;
    if (window.ResizeObserver) new ResizeObserver(function () { self.draw(); }).observe(this.wrap);
    if (this.section) this.section.addEventListener('toggle', function () { self.draw(); });
  }

  function packLanes(feats, L) {
    var tmp = feats.map(function (f) { return { pos1: f.c1, pos2: f.c2 }; });
    var n = window.mpseq.geom.lanes(tmp, L);
    tmp.forEach(function (t, i) { feats[i].lane = t.lane; });
    return n;
  }

  Viewer.prototype.load = function (p) {
    var same = this.version === p.version && this.L === p.s.length;
    this.input = p.input; this.version = p.version;
    if (!same) {
      this.s = p.s; this.r = p.r; this.L = p.s.length;
      this.refLength = p.refLength; this.rotation = p.rotation; this.rc = !!p.refRc;
      this.sNongap = G.nongap(p.s); this.rNongap = G.nongap(p.r);
      this.sCounts = G.counts(p.s); this.rCounts = G.counts(p.r);
      this.cls = G.classify(p.s, p.r);
      this.matchCounts = G.prefix(this.cls, 0); this.bothCounts = G.prefixNot(this.cls, 2);
    }
    var self = this;
    this.sFeats = (p.sampleFeatures || []).map(function (f) {
      return { row: f.row, gene: f.gene, type: f.type, pos1: f.pos1, pos2: f.pos2, dir: f.dir, track: 'sample',
               c1: G.sampleCol(f.pos1, self.sNongap), c2: G.sampleCol(f.pos2, self.sNongap) };
    });
    var ctx = { rNongap: this.rNongap, refLength: this.refLength, rotation: this.rotation, rc: this.rc };
    this.rFeats = (p.refFeatures || []).map(function (f) {
      var m = G.refFeature(f, ctx);
      return { gene: f.gene, type: f.type, pos1: f.pos1, pos2: f.pos2, dir: m.dir, track: 'ref', c1: m.c1, c2: m.c2 };
    });
    this.sLanes = packLanes(this.sFeats, this.L); this.rLanes = packLanes(this.rFeats, this.L);
    this.selected = (p.selected === undefined || p.selected === null) ? null : p.selected;
    if (!same) this.initialView();
    this.draw();
  };

  // ---- view state (linear over columns 1..L) ----
  Viewer.prototype.width = function () { return Math.max(200, this.wrap.clientWidth - GUTTER); };
  Viewer.prototype.minPpc = function () { return this.width() / this.L; };
  Viewer.prototype.viewLen = function () { return this.width() / this.ppc; };
  Viewer.prototype.clamp = function () {
    this.ppc = Math.min(MAX_PPC, Math.max(this.minPpc(), this.ppc));
    this.viewStart = Math.max(1, Math.min(this.L - this.viewLen() + 1, this.viewStart));
  };
  Viewer.prototype.whole = function () { this.ppc = this.minPpc(); this.viewStart = 1; this.clamp(); this.draw(); };
  Viewer.prototype.zoom = function (factor, atCol) {
    var centre = atCol || (this.viewStart + this.viewLen() / 2);
    var np = Math.min(MAX_PPC, Math.max(this.minPpc(), this.ppc * factor));
    this.viewStart = centre - (centre - this.viewStart) * (this.ppc / np);
    this.ppc = np; this.clamp(); this.draw();
  };
  Viewer.prototype.fitCols = function (a, b, frac) {
    var span = b - a + 1, margin = Math.max(3, span * (frac || 0.05));
    this.ppc = Math.min(MAX_PPC, this.width() / (span + 2 * margin));
    this.viewStart = a - margin; this.clamp(); this.draw();
  };
  // A short scaffold aligned to a long reference opens on the scaffold.
  Viewer.prototype.initialView = function () {
    var n = this.sNongap.length;
    if (!n) { this.whole(); return; }
    this.fitCols(this.sNongap[0], this.sNongap[n - 1], 0.02);
  };
  Viewer.prototype.fitFeature = function (f) {
    if (f.c1 <= f.c2) this.fitCols(f.c1, f.c2); else this.whole();
  };
  Viewer.prototype.fit = function (row) {
    this.selected = (row === null || row === undefined) ? null : row;
    var f = this.sFeats.find(function (x) { return x.row === row; });
    if (f) this.fitFeature(f); else this.draw();
  };
  Viewer.prototype.x = function (col) { return GUTTER + (col - this.viewStart) * this.ppc; };
  Viewer.prototype.colF = function (px) { return this.viewStart + (px - GUTTER) / this.ppc; };
  Viewer.prototype.colAt = function (px) { return Math.floor(this.colF(px)); };
  Viewer.prototype.visible = function () {
    return [Math.max(1, Math.floor(this.viewStart)), Math.min(this.L, Math.ceil(this.viewStart + this.viewLen()))];
  };

  // ---- layout ----
  Viewer.prototype.ntOn = function () { return this.ppc >= BAR; };
  Viewer.prototype.ys = function () {
    var y = {}, cur = RULER_H + PAD, nt = this.ntOn();
    y.sLanes = cur; cur += this.sLanes * LANE_H + PAD;
    y.sNt = nt ? cur : null; if (nt) cur += NT_H;
    cur += GAP; y.ident = cur; cur += IDENT_H + GAP;
    y.rNt = nt ? cur : null; if (nt) cur += NT_H;
    cur += PAD; y.rLanes = cur; cur += this.rLanes * LANE_H + PAD;
    y.rRuler = cur; y.height = cur + RULER_H;
    return y;
  };
  Viewer.prototype.height = function () { return this.ys().height; };

  // ---- drawing ----
  Viewer.prototype.draw = function () {
    if (!this.L || !this.canvas.offsetParent) return;
    this.clamp();
    var dpr = window.devicePixelRatio || 1, W = this.wrap.clientWidth, y = this.ys(), H = y.height;
    this.canvas.width = W * dpr; this.canvas.height = H * dpr; this.canvas.style.height = H + 'px';
    var c = this.canvas.getContext('2d'); c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, W, H);
    c.font = '12px ' + cssVar('--mp-font-mono', 'monospace');
    this.hits = [];
    c.save(); c.beginPath(); c.rect(GUTTER, 0, W - GUTTER, H); c.clip();
    this.drawSampleRuler(c);
    this.drawLanes(c, W, this.sFeats, y.sLanes);
    if (y.sNt !== null) this.drawNt(c, y.sNt, this.s);
    this.drawIdentity(c, W, y.ident);
    if (y.rNt !== null) this.drawNt(c, y.rNt, this.r);
    this.drawLanes(c, W, this.rFeats, y.rLanes);
    this.drawRefRuler(c, y.rRuler);
    c.restore();
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.textAlign = 'right'; c.textBaseline = 'middle';
    if (this.sLanes) c.fillText('sample', GUTTER - 6, y.sLanes + (LANE_H - 6) / 2, GUTTER - 10);
    c.fillText('identity', GUTTER - 6, y.ident + IDENT_H / 2, GUTTER - 10);
    if (this.rLanes) c.fillText('ref', GUTTER - 6, y.rLanes + (LANE_H - 6) / 2, GUTTER - 10);
  };
  Viewer.prototype.tickStyle = function (c) {
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.strokeStyle = cssVar('--mp-border', '#ccc');
    c.textAlign = 'center'; c.textBaseline = 'top'; c.lineWidth = 1;
  };
  // Ticks at round sample positions, placed at the column that holds the base.
  Viewer.prototype.drawSampleRuler = function (c) {
    var step = ui().niceStep(90 / this.ppc), v = this.visible(), n = this.sNongap.length;
    var p0 = this.sCounts[v[0] - 1] + 1, p1 = Math.min(n, this.sCounts[v[1]]);
    this.tickStyle(c);
    for (var p = Math.max(1, Math.ceil(p0 / step) * step); p <= p1; p += step) {
      var x = this.x(this.sNongap[p - 1]) + this.ppc / 2;
      c.beginPath(); c.moveTo(x, RULER_H - 6); c.lineTo(x, RULER_H); c.stroke();
      c.fillText(p.toLocaleString(), x, 2);
    }
  };
  // Ticks at round record positions: walk the drawn reference bases in view
  // and label those whose record position is a multiple of the step.
  Viewer.prototype.drawRefRuler = function (c, y) {
    var step = ui().niceStep(90 / this.ppc), v = this.visible();
    var k0 = this.rCounts[v[0] - 1] + 1, k1 = this.rCounts[v[1]];
    this.tickStyle(c);
    for (var k = k0; k <= k1; k++) {
      var col = this.rNongap[k - 1], p = G.refAt(col, this.r, this.rCounts, this.refLength, this.rotation, this.rc);
      if (p % step !== 0) continue;
      var x = this.x(col) + this.ppc / 2;
      c.beginPath(); c.moveTo(x, y); c.lineTo(x, y + 6); c.stroke();
      c.fillText(p.toLocaleString(), x, y + 8);
    }
  };
  Viewer.prototype.drawLanes = function (c, W, feats, y0) {
    var self = this, vs = this.viewStart, ve = vs + this.viewLen();
    c.textAlign = 'center'; c.textBaseline = 'middle';
    feats.forEach(function (f) {
      var col = ui().typeColor(f.type), y = y0 + f.lane * LANE_H, h = LANE_H - 6;
      G.pieces(f.c1, f.c2, self.L).forEach(function (pc) {
        if (pc[1] + 1 < vs || pc[0] > ve) return;
        var x0 = Math.max(GUTTER, self.x(pc[0])), x1 = Math.min(W, self.x(pc[1] + 1));
        if (x1 - x0 < 1) return;
        var fwd = f.dir !== '-', head = Math.min(8, x1 - x0);
        var pts = fwd ? [[x0, y], [x1 - head, y], [x1, y + h / 2], [x1 - head, y + h], [x0, y + h]]
                      : [[x1, y], [x0 + head, y], [x0, y + h / 2], [x0 + head, y + h], [x1, y + h]];
        c.beginPath(); c.moveTo(pts[0][0], pts[0][1]);
        for (var i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]);
        c.closePath();
        c.fillStyle = col + '55'; c.fill();
        var sel = f.track === 'sample' && f.row === self.selected;
        c.lineWidth = sel ? 2 : 1; c.strokeStyle = sel ? cssVar('--mp-primary', '#337ab7') : col; c.stroke();
        if (x1 - x0 > c.measureText(f.gene).width + 8) {
          c.fillStyle = cssVar('--mp-text', '#333'); c.fillText(f.gene, (x0 + x1) / 2, y + h / 2);
        }
        self.hits.push({ x0: x0, x1: x1, y0: y, y1: y + h, f: f });
      });
    });
    c.lineWidth = 1;
  };
  Viewer.prototype.drawNt = function (c, y, str) {
    var v = this.visible(), letter = this.ppc >= LETTER, w = Math.max(1, this.ppc - (letter ? 1 : 0));
    c.textAlign = 'center'; c.textBaseline = 'middle';
    for (var col = v[0]; col <= v[1]; col++) {
      var ch = str.charAt(col - 1), x = this.x(col);
      if (ch === '-') {
        if (letter) { c.fillStyle = BASE.N; c.fillText('-', x + this.ppc / 2, y + NT_H / 2); }
        continue;
      }
      c.fillStyle = (BASE[ch] || BASE.N) + (letter ? '99' : 'cc'); c.fillRect(x, y + 2, w, NT_H - 4);
      if (letter) { c.fillStyle = '#000000'; c.fillText(ch, x + this.ppc / 2, y + NT_H / 2); }
    }
  };
  Viewer.prototype.drawIdentity = function (c, W, y) {
    c.fillStyle = '#F0F0F0'; c.fillRect(GUTTER, y, W - GUTTER, IDENT_H);
    var v = this.visible();
    if (this.ntOn()) {
      var w = Math.max(1, this.ppc - (this.ppc >= LETTER ? 1 : 0));
      for (var col = v[0]; col <= v[1]; col++) {
        c.fillStyle = CLS[this.cls[col - 1]]; c.fillRect(this.x(col), y + 2, w, IDENT_H - 4);
      }
      return;
    }
    // One bin per pixel column, collected into one path so fractional
    // devicePixelRatio does not stripe the fill.
    c.fillStyle = CLS[0]; c.beginPath();
    for (var px = GUTTER; px < W; px++) {
      var c0 = Math.max(1, Math.floor(this.colF(px))), c1 = Math.min(this.L, Math.ceil(this.colF(px + 1)) - 1);
      if (c0 > this.L) break;
      if (c1 < c0) c1 = c0;
      var id = G.identity(c0, c1, this.matchCounts, this.bothCounts);
      if (id === null) continue;
      c.rect(px, y + IDENT_H * (1 - id), 1, IDENT_H * id);
    }
    c.fill();
  };

  // ---- interaction ----
  Viewer.prototype.hitTest = function (px, py) {
    for (var i = this.hits.length - 1; i >= 0; i--) {
      var h = this.hits[i];
      if (px >= h.x0 && px <= h.x1 && py >= h.y0 && py <= h.y1) return h;
    }
    return null;
  };
  Viewer.prototype.tipText = function (col, h) {
    if (h) {
      var f = h.f;
      return (f.track === 'sample' ? 'Sample ' : 'Reference ') + f.gene + ' (' + f.type + ') ' + f.dir +
        ' ' + f.pos1.toLocaleString() + '-' + f.pos2.toLocaleString();
    }
    if (col < 1 || col > this.L) return '';
    var sp = G.sampleAt(col, this.s, this.sCounts);
    var rp = G.refAt(col, this.r, this.rCounts, this.refLength, this.rotation, this.rc);
    return 'Column ' + col.toLocaleString() +
      ' | sample ' + (sp === null ? 'gap' : sp.toLocaleString() + ' ' + this.s.charAt(col - 1)) +
      ' | ref ' + (rp === null ? 'gap' : rp.toLocaleString() + ' ' + this.r.charAt(col - 1)) +
      ' | ' + CLS_NAME[this.cls[col - 1]];
  };
  Viewer.prototype.click = function (px, py) {
    var h = this.hitTest(px, py);
    if (h && h.f.track === 'ref') { this.fitFeature(h.f); return; }
    this.selected = h ? h.f.row : null;
    if (h && this.input && window.Shiny) {
      window.Shiny.setInputValue(this.input, { row: h.f.row, nonce: Date.now() }, { priority: 'event' });
    }
    this.draw();
  };
  Viewer.prototype.bindControls = function () {
    var self = this, sec = this.section || document, cv = this.canvas, dragging = null;
    cv.addEventListener('wheel', function (e) {
      e.preventDefault();
      var rect = cv.getBoundingClientRect(), px = e.clientX - rect.left;
      var pan = e.shiftKey ? e.deltaY : Math.abs(e.deltaX) > Math.abs(e.deltaY) ? e.deltaX : 0;
      if (pan) { self.viewStart += pan / self.ppc; self.clamp(); self.draw(); return; }
      self.zoom(Math.pow(1.15, -e.deltaY / 100), self.colF(px));
    }, { passive: false });
    cv.addEventListener('mousedown', function (e) { dragging = { x: e.clientX, start: self.viewStart, moved: false }; });
    window.addEventListener('mousemove', function (e) {
      if (!dragging) return;
      var dx = e.clientX - dragging.x; if (Math.abs(dx) > 2) dragging.moved = true;
      self.viewStart = dragging.start - dx / self.ppc; self.clamp(); self.draw();
    });
    window.addEventListener('mouseup', function (e) {
      if (!dragging) return;
      var rect = cv.getBoundingClientRect();
      if (!dragging.moved) self.click(e.clientX - rect.left, e.clientY - rect.top);
      dragging = null;
    });
    cv.addEventListener('mousemove', function (e) {
      if (!self.tip) return;
      var rect = cv.getBoundingClientRect(), px = e.clientX - rect.left, py = e.clientY - rect.top;
      var t = px < GUTTER ? '' : self.tipText(self.colAt(px), self.hitTest(px, py));
      if (!t) { self.tip.hidden = true; return; }
      self.tip.textContent = t; self.tip.hidden = false;
      var tw = self.tip.offsetWidth, th = self.tip.offsetHeight;
      self.tip.style.left = Math.max(0, Math.min(px + 12, cv.clientWidth - tw)) + 'px';
      self.tip.style.top = (py + 12 + th > cv.clientHeight ? py - th - 8 : py + 12) + 'px';
    });
    cv.addEventListener('mouseleave', function () { if (self.tip) self.tip.hidden = true; });
    sec.querySelectorAll('[data-mpsyn]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var a = btn.getAttribute('data-mpsyn');
        if (a === 'zoom_in') self.zoom(2); else if (a === 'zoom_out') self.zoom(0.5);
        else if (a === 'whole') self.whole(); else if (a === 'fit' && self.selected !== null) self.fit(self.selected);
      });
    });
  };

  function get(id) {
    var el = document.getElementById(id);
    if (!el) return null;
    if (viewers[id] && viewers[id].canvas !== el) delete viewers[id];
    return viewers[id] || (viewers[id] = new Viewer(id));
  }
  // The canvas lives in a uiOutput and may reach the DOM after the message.
  function withViewer(id, fn, tries) {
    var v = get(id);
    if (v) { fn(v); return; }
    if ((tries || 0) >= 40) return;
    setTimeout(function () { withViewer(id, fn, (tries || 0) + 1); }, 250);
  }
  // The section's uiOutput re-renders (picker, flags, lock, reopen) after the
  // payload has already gone to the previous canvas, so keep the last payload
  // per id and reload it into whatever canvas now carries that id.
  var cache = {}, revivePending = null;
  function revive() {
    Object.keys(cache).forEach(function (id) {
      var el = document.getElementById(id);
      if (el && (!viewers[id] || viewers[id].canvas !== el)) get(id).load(cache[id]);
    });
  }
  if (window.MutationObserver) {
    new MutationObserver(function () {
      if (revivePending) return;
      revivePending = setTimeout(function () { revivePending = null; revive(); }, 50);
    }).observe(document.documentElement, { childList: true, subtree: true });
  }
  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler('mpsyn', function (p) { cache[p.id] = p; withViewer(p.id, function (v) { v.load(p); }); });
    window.Shiny.addCustomMessageHandler('mpsyn_select', function (p) { withViewer(p.id, function (v) { v.fit(p.row); }); });
  }
  window.mpsyn.state = function (id) {
    var v = viewers[id]; if (!v) return null;
    return { L: v.L, version: v.version, viewStart: v.viewStart, ppc: v.ppc, selected: v.selected,
             height: v.height(), sampleLanes: v.sLanes, refLanes: v.rLanes,
             live: v.canvas === document.getElementById(id) };
  };
  window.mpsyn.zoom = function (id, f) { var v = get(id); if (v) v.zoom(f); };
  window.mpsyn.whole = function (id) { var v = get(id); if (v) v.whole(); };
  window.mpsyn.fit = function (id, row) { var v = get(id); if (v) v.fit(row); };
  window.mpsyn.click = function (id, px, py) { var v = get(id); if (v) v.click(px, py); };
  window.mpsyn.hitTest = function (id, px, py) { var v = get(id); return v ? v.hitTest(px, py) : null; };
  window.mpsyn.tip = function (id, col) { var v = get(id); return v ? v.tipText(col, null) : ''; };
})();
