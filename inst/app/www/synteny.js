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
