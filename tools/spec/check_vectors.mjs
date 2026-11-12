#!/usr/bin/env node
/*
 * Check a helmlab JavaScript build against the Helmlab 2026 test vectors at
 * conformance level R (spec/2026/spec.md, Section 8).
 *
 * Usage:
 *
 *     node tools/spec/check_vectors.mjs <path to dist/index.js> [vectors dir]
 *
 * The vectors dir defaults to spec/2026/vectors under the repository root
 * (two directories above this file). Every *.json file in it that has a
 * "rows" array is checked. A file without one (api-1.0-1.json, the regression
 * set of the package) is skipped. Rows with level "info" are counted and
 * skipped.
 *
 * Every row with level "U" or "R" is recomputed with the package, one colour
 * per call, and compared with the level R tolerances in the file header
 * (tolerances.R: forward and deltaE 1e-11, inverse 1e-8), absolute, per
 * channel. Forward rows take their xyz-astm-d65 input. Inverse rows are
 * compared in the xyz-astm-d65 frame. helmgenlch.fwd is compared as L, C and
 * the hue chord C*|dh|, with dh the hue difference in radians wrapped to
 * [-pi, pi] and C the expected chroma.
 *
 * Package calls, with gen = new Helmlab().gen, its GenSpace gen.space, and
 * p = compileParams(getDefaultParams()):
 *
 *   helmgen.fwd                   gen.space.fromXYZ(xyz)
 *   helmgen.inv                   gen.space.toXYZ(lab)
 *   helmgenlch.fwd                gen.toLch(gen.space.fromXYZ(xyz))
 *   helmgenlch.inv                gen.space.toXYZ(gen.fromLch(lch))
 *   helmlab-metric.fwd / .inv     new MetricSpace(p, { neutralCorrection: false }).fromXYZ / .toXYZ
 *   helmlab-metric-neutral.*      new MetricSpace(p, { neutralCorrection: true }).fromXYZ / .toXYZ
 *   deltaE                        new MetricSpace(p, { neutralCorrection: false, abRotateDeg: 0 })
 *                                   .distanceFromLab(lab1, lab2)
 *
 * deltaE: MetricSpace.distanceFromLab is an alias of MetricSpace.distance.
 * Before Steps D1 to D7 (spec Section 5.6) it undoes the rotation of its own
 * object and, if that object has the neutral correction on, adds the
 * correction back. On an object with neither (abRotateDeg 0, neutral
 * correction off) it evaluates Steps D1 to D7 on its two inputs as given.
 * That is the route counted here. The base-form object (rotation -28.2
 * degrees) gives the same value up to the rounding of the rotation; it is
 * printed as an info line and does not count towards pass or fail.
 *
 * Output: per file and op, the rows checked, the info rows skipped, the rows
 * whose difference is 0, the largest difference, the tolerance, the failures
 * and the worst rows; for an op with failures, the worst failing rows with
 * expected and actual values.
 *
 * Exit status: 0 when every checked row is within tolerance, 1 when a row
 * fails, 2 on a usage or input error.
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const REPO = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..');
const DEFAULT_VECTORS = path.join(REPO, 'spec', '2026', 'vectors');
const SPEC_R = { forward: 1e-11, inverse: 1e-8, deltaE: 1e-11 }; // spec Section 8; only used to flag a different header
const WORST = 3; // worst rows listed per file and op
const SHOW_FAIL = 5; // failing rows printed with their values per file and op

function fatal(message) {
  console.error(`check_vectors.mjs: ${message}`);
  process.exit(2);
}

// ── Package ─────────────────────────────────────────────────────────

async function loadPackage(file) {
  const url = pathToFileURL(path.resolve(file)).href;
  let H;
  try {
    H = await import(url);
  } catch (e) {
    fatal(`cannot load ${file}: ${e.message}`);
  }
  for (const name of ['Helmlab', 'MetricSpace', 'compileParams', 'getDefaultParams']) {
    if (typeof H[name] !== 'function') fatal(`${file} has no export ${name}`);
  }
  const gen = new H.Helmlab().gen;
  const genSpace = gen.space;
  if (!genSpace || typeof genSpace.fromXYZ !== 'function' || typeof genSpace.toXYZ !== 'function') {
    fatal('new Helmlab().gen.space has no fromXYZ and toXYZ');
  }
  if (typeof gen.toLch !== 'function' || typeof gen.fromLch !== 'function') {
    fatal('new Helmlab().gen has no toLch and fromLch');
  }
  const p = H.compileParams(H.getDefaultParams());
  const base = new H.MetricSpace(p, { neutralCorrection: false });
  const neutral = new H.MetricSpace(p, { neutralCorrection: true });
  const deltaE = new H.MetricSpace(p, { neutralCorrection: false, abRotateDeg: 0 });
  if (typeof deltaE.distanceFromLab !== 'function') fatal('MetricSpace has no distanceFromLab');

  // kind: which tolerance of the header applies; lch: compare the hue as a chord.
  const ops = new Map([
    ['helmgen.fwd', { kind: 'forward', run: (x) => genSpace.fromXYZ(x) }],
    ['helmgen.inv', { kind: 'inverse', run: (v) => genSpace.toXYZ(v) }],
    ['helmgenlch.fwd', { kind: 'forward', lch: true, run: (x) => gen.toLch(genSpace.fromXYZ(x)) }],
    ['helmgenlch.inv', { kind: 'inverse', run: (v) => genSpace.toXYZ(gen.fromLch(v)) }],
    ['helmlab-metric.fwd', { kind: 'forward', run: (x) => base.fromXYZ(x) }],
    ['helmlab-metric.inv', { kind: 'inverse', run: (v) => base.toXYZ(v) }],
    ['helmlab-metric-neutral.fwd', { kind: 'forward', run: (x) => neutral.fromXYZ(x) }],
    ['helmlab-metric-neutral.inv', { kind: 'inverse', run: (v) => neutral.toXYZ(v) }],
    ['deltaE', { kind: 'deltaE', run: (lab1, lab2) => deltaE.distanceFromLab(lab1, lab2) }],
  ]);
  // Info only: deltaE through the base-form object, which undoes its rotation first.
  const deltaEBaseObject = (lab1, lab2) => base.distanceFromLab(lab1, lab2);
  const genSpaceClass = typeof H.GenSpace === 'function' && genSpace instanceof H.GenSpace ? 'GenSpace' : 'not a GenSpace instance';
  return { H, ops, deltaEBaseObject, genSpaceClass };
}

// ── Values and differences ──────────────────────────────────────────

/** A JSON number, or one of the strings "NaN", "Infinity", "-Infinity" (allowed in info rows). */
function num(v, what) {
  if (typeof v === 'number') return v;
  if (v === 'NaN') return NaN;
  if (v === 'Infinity') return Infinity;
  if (v === '-Infinity') return -Infinity;
  throw new Error(`${what} is not a number: ${JSON.stringify(v)}`);
}

function triple(v, what) {
  if (!Array.isArray(v) || v.length !== 3) throw new Error(`${what} is not a list of 3 numbers`);
  return v.map((t) => num(t, what));
}

/** Largest absolute difference over the channels; NaN counts as Infinity. */
function maxAbsDiff(got, exp) {
  if (got.length !== exp.length) return Infinity;
  let worst = 0;
  for (let k = 0; k < exp.length; k++) {
    const d = Math.abs(got[k] - exp[k]);
    if (!(d <= worst)) worst = Number.isNaN(d) ? Infinity : d;
  }
  return worst;
}

/** L, C and the hue chord C*|dh|, dh in radians wrapped to [-pi, pi], C the expected chroma. */
function lchDiff(got, exp) {
  if (got.length !== 3) return Infinity;
  let dh = (got[2] - exp[2]) % 360; // degrees, in (-360, 360)
  if (dh > 180) dh -= 360;
  else if (dh < -180) dh += 360;
  const parts = [Math.abs(got[0] - exp[0]), Math.abs(got[1] - exp[1]), Math.abs(exp[1]) * Math.abs(dh * (Math.PI / 180))];
  if (parts.some(Number.isNaN)) return Infinity;
  return Math.max(...parts);
}

/** Recompute one U or R row. Returns { got, exp, diff }; throws on a malformed row or a package error. */
function evaluate(ops, row) {
  const op = ops.get(row.op);
  if (!op) throw new Error(`unknown op ${JSON.stringify(row.op)}`);
  let got, exp, diff;
  if (op.kind === 'forward') {
    const x = triple(row.in?.['xyz-astm-d65'], 'in.xyz-astm-d65');
    exp = triple(row.out, 'out');
    if (!exp.every(Number.isFinite)) throw new Error('expected output not finite in a U or R row');
    got = Array.from(op.run(x));
    diff = op.lch ? lchDiff(got, exp) : maxAbsDiff(got, exp);
  } else if (op.kind === 'inverse') {
    const v = triple(row.in, 'in');
    exp = triple(row.out?.['xyz-astm-d65'], 'out.xyz-astm-d65');
    if (!exp.every(Number.isFinite)) throw new Error('expected output not finite in a U or R row');
    got = Array.from(op.run(v));
    diff = maxAbsDiff(got, exp);
  } else {
    const lab1 = triple(row.in?.lab1, 'in.lab1');
    const lab2 = triple(row.in?.lab2, 'in.lab2');
    exp = [num(row.out, 'out')];
    if (!Number.isFinite(exp[0])) throw new Error('expected output not finite in a U or R row');
    got = [op.run(lab1, lab2)];
    diff = maxAbsDiff(got, exp);
  }
  return { got, exp, diff };
}

// ── Reporting helpers ───────────────────────────────────────────────

/** Keep the k items with the largest diff, largest first. */
function keepWorst(list, item, k) {
  if (list.length === k && !(item.diff > list[k - 1].diff)) return;
  list.push(item);
  list.sort((p, q) => q.diff - p.diff);
  if (list.length > k) list.pop();
}

const fmtDiff = (d) => (d === 0 ? '0' : d.toExponential(2));
const fmtTol = (t) => t.toExponential(0);
const fmtVals = (v) => `[${v.map((t) => String(t)).join(', ')}]`;

function newStats() {
  return { checked: 0, info: 0, equal: 0, largest: 0, failures: 0, worst: [], fails: [] };
}

// ── Main ────────────────────────────────────────────────────────────

async function main() {
  const args = process.argv.slice(2);
  const usage = 'usage: node tools/spec/check_vectors.mjs <path to dist/index.js> [vectors dir]';
  if (args.length === 1 && (args[0] === '-h' || args[0] === '--help')) {
    console.log(usage);
    process.exit(0);
  }
  if (args.length < 1 || args.length > 2) fatal(usage);
  const pkgFile = path.resolve(args[0]);
  const vectorDir = path.resolve(args[1] ?? DEFAULT_VECTORS);
  if (!fs.existsSync(pkgFile)) fatal(`no such file: ${pkgFile}`);
  if (!fs.existsSync(vectorDir) || !fs.statSync(vectorDir).isDirectory()) fatal(`no such directory: ${vectorDir}`);

  const { H, ops, deltaEBaseObject, genSpaceClass } = await loadPackage(pkgFile);
  console.log(`package: ${pkgFile}`);
  console.log(`version: ${H.VERSION ?? 'unknown'}, Node.js ${process.version}`);
  console.log(`vectors: ${vectorDir}`);
  console.log(`GenSpace: new Helmlab().gen.space (${genSpaceClass}); LCh: new Helmlab().gen.toLch / fromLch`);
  console.log('MetricSpace: base form { neutralCorrection: false }, neutral-corrected form { neutralCorrection: true }');
  console.log('deltaE: distanceFromLab of new MetricSpace(p, { neutralCorrection: false, abRotateDeg: 0 }) (Steps D1 to D7 on the input)');

  const files = fs.readdirSync(vectorDir).filter((f) => f.endsWith('.json')).sort();
  if (files.length === 0) fatal(`no .json files in ${vectorDir}`);

  const head = `${'file'.padEnd(16)} ${'op'.padEnd(27)} ${'checked'.padStart(7)} ${'info'.padStart(5)} ` +
    `${'equal'.padStart(6)} ${'largest'.padStart(9)} ${'tol'.padStart(6)} ${'fail'.padStart(5)}`;
  console.log('');
  console.log(head);

  const total = { files: 0, checked: 0, info: 0, failures: 0 };
  for (const name of files) {
    let data;
    try {
      data = JSON.parse(fs.readFileSync(path.join(vectorDir, name), 'utf8'));
    } catch (e) {
      fatal(`cannot read ${name}: ${e.message}`);
    }
    if (!data || !Array.isArray(data.rows)) {
      console.log(`${name.padEnd(16)} skipped: no "rows" (not in the vector format)`);
      continue;
    }
    const R = data.tolerances?.R;
    const tol = { forward: R?.forward, inverse: R?.inverse, deltaE: R?.deltaE };
    for (const [k, v] of Object.entries(tol)) {
      if (typeof v !== 'number' || !(v > 0)) fatal(`${name}: header has no tolerances.R.${k}`);
    }
    if (Object.keys(SPEC_R).some((k) => tol[k] !== SPEC_R[k])) {
      console.log(`${name.padEnd(16)} note: header tolerances ${JSON.stringify(tol)} differ from spec Section 8 ${JSON.stringify(SPEC_R)}`);
    }
    total.files++;

    const stats = new Map(); // op -> stats, in order of first appearance
    const baseObject = { checked: 0, largest: 0, worstId: '', sameAsCounted: 0 }; // info route for deltaE
    for (const row of data.rows) {
      const opName = String(row.op);
      if (!stats.has(opName)) stats.set(opName, newStats());
      const s = stats.get(opName);
      if (row.level === 'info') {
        s.info++;
        continue;
      }
      s.checked++;
      let res = null;
      let problem = null;
      if ((row.level !== 'U' && row.level !== 'R') || (row.level === 'U') !== (row.domain === 'in')) {
        problem = `level ${JSON.stringify(row.level)} with domain ${JSON.stringify(row.domain)}`;
      } else {
        try {
          res = evaluate(ops, row);
        } catch (e) {
          problem = `${e.name}: ${e.message}`;
        }
      }
      const diff = res ? res.diff : Infinity;
      const kind = ops.get(opName)?.kind;
      const limit = kind ? tol[kind] : 0;
      if (diff === 0) s.equal++;
      if (!(diff <= s.largest)) s.largest = diff;
      if (diff > 0) keepWorst(s.worst, { id: row.id, diff }, WORST);
      if (problem !== null || !(diff <= limit)) {
        s.failures++;
        keepWorst(s.fails, { id: row.id, diff, problem, res }, SHOW_FAIL);
      }
      if (opName === 'deltaE' && res) {
        let v;
        try {
          v = deltaEBaseObject(triple(row.in.lab1, 'in.lab1'), triple(row.in.lab2, 'in.lab2'));
        } catch {
          v = NaN; // counts as an infinite difference
        }
        const d = maxAbsDiff([v], res.exp);
        baseObject.checked++;
        if (!(d <= baseObject.largest)) {
          baseObject.largest = d;
          baseObject.worstId = row.id;
        }
        if (Object.is(v, res.got[0])) baseObject.sameAsCounted++;
      }
    }

    for (const [opName, s] of stats) {
      const kind = ops.get(opName)?.kind;
      const limit = kind ? tol[kind] : NaN;
      console.log(`${name.padEnd(16)} ${opName.padEnd(27)} ${String(s.checked).padStart(7)} ${String(s.info).padStart(5)} ` +
        `${String(s.equal).padStart(6)} ${fmtDiff(s.largest).padStart(9)} ${(kind ? fmtTol(limit) : '-').padStart(6)} ` +
        `${String(s.failures).padStart(5)}`);
      const worst = s.worst.map((w) => `${w.id} ${fmtDiff(w.diff)}`).join(', ');
      console.log(`${''.padEnd(17)}worst: ${worst || '-'}`);
      for (const f of s.fails) {
        if (f.problem !== null) {
          console.log(`${''.padEnd(17)}FAIL ${f.id}: ${f.problem}`);
        } else {
          console.log(`${''.padEnd(17)}FAIL ${f.id}: difference ${fmtDiff(f.diff)} > ${fmtTol(limit)}`);
          console.log(`${''.padEnd(22)}expected ${fmtVals(f.res.exp)}`);
          console.log(`${''.padEnd(22)}actual   ${fmtVals(f.res.got)}`);
        }
      }
      if (s.failures > s.fails.length) {
        console.log(`${''.padEnd(17)}(${s.failures - s.fails.length} more failing rows not shown)`);
      }
      total.checked += s.checked;
      total.info += s.info;
      total.failures += s.failures;
    }
    if (baseObject.checked > 0) {
      console.log(`${''.padEnd(17)}info, not counted: deltaE from the base-form object (rotation -28.2, undone inside ` +
        `distanceFromLab): largest difference ${fmtDiff(baseObject.largest)}` +
        `${baseObject.worstId ? ` at ${baseObject.worstId}` : ''}; equal to the counted route in ` +
        `${baseObject.sameAsCounted} of ${baseObject.checked} rows`);
    }
  }

  if (total.files === 0) fatal(`no vector files (with "rows") in ${vectorDir}`);
  console.log('');
  console.log(`total: ${total.files} files, ${total.checked} rows checked, ${total.info} info rows skipped, ` +
    `${total.failures} failures`);
  process.exit(total.failures > 0 ? 1 : 0);
}

await main();
