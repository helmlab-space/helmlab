// Check a third-party JavaScript implementation against the Helmlab 2026 vectors.
// Usage: node run_js.mjs <colorjs|colorspace> <package root> <vectors dir> <css|astm>
// Prints one JSON object: per op, per level (U, R), the row count, the largest error,
// the worst row and, for level U, the rows over the level-U tolerance.
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const [kind, root, vecDir, frameArg] = process.argv.slice(2);
const frame = frameArg === 'astm' ? 'xyz-astm-d65' : 'xyz-css-d65';
const TOL_U = { fwd: 1e-10, inv: 1e-8, de: 1e-10 };
const imp = (p) => import(pathToFileURL(p).href);

let fwd, inv, dE, label;
if (kind === 'colorjs') {
  const pkgRoot = `${root}/node_modules/colorjs.io`;
  const Color = (await imp(`${pkgRoot}/src/index.js`)).default;
  label = `colorjs.io ${JSON.parse(readFileSync(`${pkgRoot}/package.json`)).version}`;
  const ids = { 'helmgen': 'helmgen', 'helmgenlch': 'helmgenlch', 'helmlab-metric': 'helmlab-metric' };
  fwd = (space, xyz) => ids[space] ? new Color('xyz-d65', xyz).to(ids[space]).coords.map(Number) : null;
  inv = (space, lab) => ids[space] ? new Color(ids[space], lab).to('xyz-d65').coords.map(Number) : null;
  dE = (x1, x2) => new Color('xyz-d65', x1).deltaE(new Color('xyz-d65', x2), 'Helmlab');
} else if (kind === 'colorspace') {
  const pkgRoot = `${root}/node_modules/color-space`;
  label = `color-space ${JSON.parse(readFileSync(`${pkgRoot}/package.json`)).version}`;
  const xyz = (await imp(`${pkgRoot}/spaces/xyz.js`)).default;
  const hg = (await imp(`${pkgRoot}/spaces/helmgen.js`)).default;
  const hgl = (await imp(`${pkgRoot}/spaces/helmgenlch.js`)).default;
  const hm = (await imp(`${pkgRoot}/spaces/helmlab-metric.js`)).default;
  const s100 = (v) => v.map((t) => t * 100), d100 = (v) => v.map((t) => t / 100);
  // color-space ships one MetricSpace form under the name helmlab-metric; it is checked
  // against both the base-form rows and the neutral-corrected rows.
  fwd = (space, x) => {
    if (space === 'helmgen') return xyz.helmgen(...s100(x));
    if (space === 'helmgenlch') return hg.helmgenlch(...xyz.helmgen(...s100(x)));
    if (space === 'helmlab-metric' || space === 'helmlab-metric-neutral') return xyz['helmlab-metric'](...s100(x));
    return null;
  };
  inv = (space, lab) => {
    if (space === 'helmgen') return d100(hg.xyz(...lab));
    if (space === 'helmgenlch') return d100(hg.xyz(...hgl.helmgen(...lab)));
    if (space === 'helmlab-metric' || space === 'helmlab-metric-neutral') return d100(hm.xyz(...lab));
    return null;
  };
  dE = null;
} else {
  throw new Error(`unknown kind ${kind}`);
}

const num = (v) => (typeof v === 'string' ? Number(v) : v);
function err(op, got, exp) {
  if (!got || got.length !== 3) return Infinity;
  if (op === 'helmgenlch.fwd') {
    const dL = Math.abs(got[0] - exp[0]), dC = Math.abs(got[1] - exp[1]);
    let chord = 0;
    if (Number.isFinite(got[2]) && Number.isFinite(exp[2])) {
      let dh = ((got[2] - exp[2]) * Math.PI) / 180;
      dh -= 2 * Math.PI * Math.round(dh / (2 * Math.PI));
      chord = Math.abs(exp[1]) * Math.abs(dh);
    } else if (Number.isNaN(got[2]) && Number.isFinite(exp[2])) {
      chord = Math.abs(exp[1]) * Math.PI; // an undefined hue is charged the largest chord, C * pi
    } else if (Math.abs(exp[1]) > 0) chord = Infinity;
    const e = Math.max(dL, dC, chord);
    return Number.isNaN(e) ? Infinity : e;
  }
  let e = 0;
  for (let i = 0; i < 3; i++) {
    const d = Math.abs(got[i] - exp[i]);
    if (!(d <= e)) e = Number.isNaN(d) ? Infinity : d;
  }
  return e;
}

const res = {};
const add = (op, level, e, id) => {
  const r = (res[op] ??= {});
  const s = (r[level] ??= { n: 0, max: 0, worst: null, over: 0 });
  s.n++;
  if (e > s.max || (s.worst === null)) { if (e >= s.max) { s.max = e; s.worst = id; } }
  const kindTol = op === 'deltaE' ? TOL_U.de : op.endsWith('.inv') ? TOL_U.inv : TOL_U.fwd;
  if (level === 'U' && !(e <= kindTol)) s.over++;
};

for (const name of ['core-1', 'conv-1', 'grey-1', 'edge-1', 'de-1', 'neutral-1']) {
  const rows = JSON.parse(readFileSync(`${vecDir}/${name}.json`)).rows;
  for (const r of rows) {
    if (r.level !== 'U' && r.level !== 'R') continue;
    const op = r.op;
    let e;
    try {
      if (op === 'deltaE') {
        if (!dE) continue;
        e = Math.abs(dE(r.in.xyz1[frame], r.in.xyz2[frame]) - num(r.out));
        if (Number.isNaN(e)) e = Infinity;
      } else {
        const space = op.slice(0, op.lastIndexOf('.'));
        if (op.endsWith('.fwd')) {
          const got = fwd(space, r.in[frame]);
          if (got === null) continue;
          e = err(op, got.map(Number), r.out.map(num));
        } else {
          const got = inv(space, r.in.map(num));
          if (got === null) continue;
          e = err(op, got.map(Number), r.out[frame].map(num));
        }
      }
    } catch (ex) {
      e = Infinity;
    }
    add(op, r.level, e, `${name} ${r.id}`);
  }
}
console.log(JSON.stringify({ implementation: label, frame, results: res }, (k, v) => (typeof v === "number" && !Number.isFinite(v) ? String(v) : v)));
