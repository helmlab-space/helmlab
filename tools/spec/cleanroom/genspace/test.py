"""Compare the clean-room GenSpace (impl.py) with the helmlab reference used as a black box.

Usage: python3 test.py [explicit|blas]   (matrix-vector mode of impl.py; default explicit)
"""
import json
import sys

import os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))   # reference code from the source tree
sys.path.insert(1, HERE)

import numpy as np
import helmlab
from helmlab import Helmlab
from helmlab.utils.srgb_convert import sRGB_to_XYZ    # the package's own sRGB -> XYZ

import impl

MODE = sys.argv[1] if len(sys.argv) > 1 else 'explicit'
impl.MATVEC_MODE = MODE

print('helmlab loaded from:', helmlab.__file__)
assert helmlab.__file__.startswith(os.path.join(ROOT, 'src') + os.sep), helmlab.__file__
print('impl MATVEC_MODE =', impl.MATVEC_MODE, '| INV_SOURCE =', impl.INV_SOURCE)

hl = Helmlab()
gen = hl.gen.space      # GenSpace with gen_params.json, neutral_correction=False  (= helmgen)

# ---------------- reference outputs (REFERENCE API, one batch per set) ----------------
k = np.arange(17) / 16
R, G, B = np.meshgrid(k, k, k, indexing='ij')
rgb = np.stack([R.ravel(), G.ravel(), B.ravel()], axis=-1)       # (4913, 3)
xyz_grid = sRGB_to_XYZ(rgb)
lab_grid = gen.from_XYZ(xyz_grid)                                  # helmgen forward
lch_grid = np.array([hl.gen.to_lch(lab) for lab in lab_grid])      # helmgenlch [L, C, h in degrees]

W = np.array([0.95047, 1.0, 1.08883])
Y = np.linspace(0.0, 10.0, 1001)
xyz_grey = np.outer(Y, W)
lab_grey = gen.from_XYZ(xyz_grey)

rng = np.random.default_rng(2026)
xyz_rand = rng.uniform(-0.05, 1.2, size=(20000, 3))
lab_rand = gen.from_XYZ(xyz_rand)

xyz_grid_inv = gen.to_XYZ(lab_grid)
xyz_rand_inv = gen.to_XYZ(lab_rand)

lab_grid_rt = np.array([np.asarray(hl.gen.from_lch(hl.gen.to_lch(lab))) for lab in lab_grid])

for nm, arr in [('xyz_grid', xyz_grid), ('lab_grid', lab_grid), ('lch_grid', lch_grid), ('lab_grey', lab_grey),
                ('lab_rand', lab_rand), ('xyz_grid_inv', xyz_grid_inv), ('xyz_rand_inv', xyz_rand_inv),
                ('lab_grid_rt', lab_grid_rt)]:
    arr = np.asarray(arr)
    print(f'  ref {nm:13s} shape={arr.shape} dtype={arr.dtype}')

# ---------------- my outputs ----------------
my_lab_grid = impl.from_XYZ(xyz_grid)
my_lch_grid = impl.to_lch(my_lab_grid)
my_lab_grey = impl.from_XYZ(xyz_grey)
my_lab_rand = impl.from_XYZ(xyz_rand)
my_xyz_grid_inv = impl.to_XYZ(lab_grid)
my_xyz_rand_inv = impl.to_XYZ(lab_rand)
my_lab_grid_rt = impl.from_lch(impl.to_lch(lab_grid))


def fmt(v):
    return '(' + ', '.join(repr(float(x)) for x in np.ravel(v)) + ')'


def compare(name, mine, ref, inputs, tol, input_label, comp_names, mask=None):
    """Max abs difference over all components; one-sided non-finite -> 1e300."""
    mine = np.asarray(mine, dtype=np.float64)
    ref = np.asarray(ref, dtype=np.float64)
    assert mine.shape == ref.shape, (name, mine.shape, ref.shape)
    n = mine.shape[0]
    fin_m = np.isfinite(mine)
    fin_r = np.isfinite(ref)
    both_nan = np.isnan(mine) & np.isnan(ref)
    same_inf = np.isinf(mine) & np.isinf(ref) & (np.sign(mine) == np.sign(ref))
    bad = ~(fin_m & fin_r) & ~both_nan & ~same_inf
    if mask is not None:
        bad &= mask
    with np.errstate(invalid='ignore'):
        diff = np.where(fin_m & fin_r, np.abs(mine - ref), 0.0)
    if mask is not None:
        diff = np.where(mask, diff, 0.0)
    nbad = int(bad.sum())
    if nbad:
        idx = np.argwhere(bad)[0]
        i, c = int(idx[0]), int(idx[1])
        worst = (f'{nbad} non-finite mismatches; first: row {i} {input_label}={fmt(inputs[i])} '
                 f'component {comp_names[c]}: mine={mine[i, c]!r} ref={ref[i, c]!r}')
        maxd = 1e300
    else:
        flat = int(np.argmax(diff))
        i, c = divmod(flat, diff.shape[1])
        maxd = float(diff[i, c])
        worst = (f'row {i} {input_label}={fmt(inputs[i])} component {comp_names[c]}: '
                 f'mine={mine[i, c]!r} ref={ref[i, c]!r}')
    ok = maxd <= tol
    print(f'{name:28s} n={n:6d} max_abs_diff={maxd:.3e} tol={tol:.0e} {"PASS" if ok else "FAIL"}  worst: {worst}')
    nexact = int(np.sum(np.all((mine == ref) | (np.isnan(mine) & np.isnan(ref)), axis=1)))
    print(f'{"":28s} rows bitwise equal: {nexact}/{n}')
    return dict(quantity=name, n=int(n), max_abs_diff=maxd, tolerance=tol, worst_input=worst), ok


results = []
allok = True

r, ok = compare('grid forward Lab (helmgen)', my_lab_grid, lab_grid, rgb, 1e-11, 'sRGB', ['L', 'a', 'b'])
results.append(r); allok &= ok

# LCh: L, C always; h when reference C >= 1e-4; a and b instead of h when C < 1e-4 (Section 8 rule quoted in 5.1.6)
lch_m = np.asarray(my_lch_grid, dtype=np.float64)
lch_r = np.asarray(lch_grid, dtype=np.float64)
small = lch_r[:, 1] < 1e-4
print(f'grid rows with reference C < 1e-4: {int(small.sum())}')
mine5 = np.column_stack([lch_m[:, 0], lch_m[:, 1], np.where(small, 0.0, lch_m[:, 2]),
                         np.where(small, my_lab_grid[:, 1], 0.0), np.where(small, my_lab_grid[:, 2], 0.0)])
ref5 = np.column_stack([lch_r[:, 0], lch_r[:, 1], np.where(small, 0.0, lch_r[:, 2]),
                        np.where(small, lab_grid[:, 1], 0.0), np.where(small, lab_grid[:, 2], 0.0)])
r, ok = compare('grid forward LCh (helmgenlch)', mine5, ref5, rgb, 1e-11, 'sRGB',
                ['L', 'C', 'h_deg', 'a (C<1e-4)', 'b (C<1e-4)'])
r['quantity'] = 'grid forward LCh (helmgenlch; h compared only where C >= 1e-4, a/b otherwise)'
results.append(r); allok &= ok
# raw h over all rows, for information
dh = np.abs(lch_m[:, 2] - lch_r[:, 2])
dh_circ = np.minimum(dh, 360.0 - dh)
i_raw = int(np.argmax(dh))
print(f'  info: raw |dh| over ALL grid rows max={dh.max():.3e} at row {i_raw} sRGB={fmt(rgb[i_raw])} '
      f'mine h={lch_m[i_raw, 2]!r} ref h={lch_r[i_raw, 2]!r} C_ref={lch_r[i_raw, 1]!r}; circular max={dh_circ.max():.3e}')
print(f'  info: raw |dh| over rows with C >= 1e-4: max={dh[~small].max():.3e}')
print(f'  info: h == 360 count mine={int(np.sum(lch_m[:, 2] == 360.0))} ref={int(np.sum(lch_r[:, 2] == 360.0))}; '
      f'-0 count mine={int(np.sum(np.signbit(lch_m[:, 2]) & (lch_m[:, 2] == 0)))} '
      f'ref={int(np.sum(np.signbit(lch_r[:, 2]) & (lch_r[:, 2] == 0)))}')
print(f'  info: LCh rows bitwise equal (all 3 cols incl. h of greys): '
      f'{int(np.sum(np.all(lch_m == lch_r, axis=1)))}/{len(lch_r)}')

r, ok = compare('grey forward Lab', my_lab_grey, lab_grey, xyz_grey, 1e-11, 'XYZ', ['L', 'a', 'b'])
results.append(r); allok &= ok

r, ok = compare('random forward Lab', my_lab_rand, lab_rand, xyz_rand, 1e-11, 'XYZ', ['L', 'a', 'b'])
results.append(r); allok &= ok

r, ok = compare('grid inverse XYZ (input = ref grid Lab)', my_xyz_grid_inv, xyz_grid_inv, lab_grid, 1e-8, 'Lab',
                ['X', 'Y', 'Z'])
results.append(r); allok &= ok

r, ok = compare('random inverse XYZ (input = ref random Lab)', my_xyz_rand_inv, xyz_rand_inv, lab_rand, 1e-8, 'Lab',
                ['X', 'Y', 'Z'])
results.append(r); allok &= ok

r, ok = compare('grid LCh round trip Lab->LCh->Lab', my_lab_grid_rt, lab_grid_rt, lab_grid, 1e-11, 'Lab',
                ['L', 'a', 'b'])
results.append(r); allok &= ok

print('\nPASS' if allok else '\nFAIL')
