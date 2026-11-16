"""Compare the clean-room implementation (impl.py) with the reference package,
which is used ONLY as a black box through the calls listed in REFERENCE API."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import impl  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
from helmlab.spaces.metric import MetricSpace, MetricParams  # noqa: E402
from helmlab.utils.srgb_convert import sRGB_to_XYZ  # noqa: E402

PATH = os.path.join(ROOT, 'spec', '2026', 'params', 'metric_params.json')
TOL = 1e-11
P = impl.load_params(PATH)

params = MetricParams.load(PATH)
base = MetricSpace(params)
disp = MetricSpace(params, neutral_correction=True)
W = np.array([0.95047, 1.0, 1.08883])

results = []
COMP = {3: ('L', 'a', 'b')}


def fmt(v):
    return '(' + ', '.join(repr(float(t)) for t in np.atleast_1d(v)) + ')'


def compare(quantity, mine, ref, inp, inp_name, comp_names, tol=TOL):
    mine = np.asarray(mine, dtype=np.float64)
    ref = np.asarray(ref, dtype=np.float64)
    assert mine.shape == ref.shape, (quantity, mine.shape, ref.shape)
    n = mine.shape[0]
    fm, fr = np.isfinite(mine), np.isfinite(ref)
    both_nan = np.isnan(mine) & np.isnan(ref)
    same_inf = np.isinf(mine) & np.isinf(ref) & (mine == ref)
    one_sided = (~fm | ~fr) & ~(both_nan | same_inf)
    if one_sided.any():
        pos = tuple(np.argwhere(one_sided)[0])
        row = pos[0]
        comp = comp_names[pos[1]] if len(pos) > 1 else comp_names[0]
        worst = (f'row {row}, {inp_name}={fmt(inp[row])}, {comp}: impl={float(mine[pos])!r}, '
                 f'ref={float(ref[pos])!r} (non-finite on one side; {int(one_sided.sum())} such values)')
        mad = 1e300
    else:
        diff = np.where(fm & fr, np.abs(mine - ref), 0.0)
        flat = int(np.argmax(diff))
        pos = np.unravel_index(flat, diff.shape)
        row = pos[0]
        comp = comp_names[pos[1]] if len(pos) > 1 else comp_names[0]
        mad = float(diff[pos])
        worst = (f'row {row}, {inp_name}={fmt(inp[row])}, {comp}: impl={float(mine[pos])!r}, '
                 f'ref={float(ref[pos])!r}')
        if mad == 0:
            worst = f'all {mine.size} values bit-identical (e.g. {worst})'
    rec = dict(quantity=quantity, n=int(n), max_abs_diff=mad, tolerance=tol, worst_input=worst)
    results.append(rec)
    status = 'PASS' if mad <= tol else 'FAIL'
    print(f'[{status}] {quantity}: n={n} max_abs_diff={mad:.3e} tol={tol:g}\n        worst: {worst}')
    return rec


# ------------------------------------------------------------- inputs
k = np.arange(17) / 16.0
R, G, B = np.meshgrid(k, k, k, indexing='ij')
rgb = np.stack([R.ravel(), G.ravel(), B.ravel()], axis=-1)
XYZ_grid = sRGB_to_XYZ(rgb)
XYZ_grey = np.outer(np.linspace(0.0, 1.0, 1001), W)
XYZ_rand = np.random.default_rng(2026).uniform(0.0, 1.0, (20000, 3))

# ----------------------------------------------------- reference (black box)
lab_grid = base.from_XYZ(XYZ_grid)
lab_grey = base.from_XYZ(XYZ_grey)
lab_rand = base.from_XYZ(XYZ_rand)
XYZ_inv_grid = base.to_XYZ(lab_grid)
XYZ_inv_rand = base.to_XYZ(lab_rand)
rng = np.random.default_rng(7)
idx = rng.integers(0, len(lab_grid), size=(5000, 2))
dE_ref = base.distance_from_lab(lab_grid[idx[:, 0]], lab_grid[idx[:, 1]])

# ------------------------------------------------------------ comparisons
LAB = ('L', 'a', 'b')
XYZn = ('X', 'Y', 'Z')
compare('forward base, grid (sRGB 17^3)', impl.forward(XYZ_grid, P), lab_grid,
        XYZ_grid, 'XYZ', LAB)
compare('forward base, grey (D65 * Y, 1001 in [0,1])', impl.forward(XYZ_grey, P), lab_grey,
        XYZ_grey, 'XYZ', LAB)
compare('forward base, random (20000, rng 2026)', impl.forward(XYZ_rand, P), lab_rand,
        XYZ_rand, 'XYZ', LAB)
compare('inverse base, grid (input = ref forward)', impl.inverse(lab_grid, P), XYZ_inv_grid,
        lab_grid, 'Lab', XYZn)
compare('inverse base, random (input = ref forward)', impl.inverse(lab_rand, P), XYZ_inv_rand,
        lab_rand, 'Lab', XYZn)
pairs = np.concatenate([lab_grid[idx[:, 0]], lab_grid[idx[:, 1]]], axis=1)
compare('delta-E helmlab, 5000 grid pairs (rng 7), on ref base coords',
        impl.delta_e(lab_grid[idx[:, 0]], lab_grid[idx[:, 1]], P), dE_ref,
        pairs, 'Lab1+Lab2', ('dE',))

# ------------------------------------------------------- neutral-corrected form
NT = impl.NeutralTable(P)
lab_grid_disp = disp.from_XYZ(XYZ_grid)
XYZ_grey3 = np.outer(np.linspace(0.0, 3.0, 601), W)
lab_grey3_disp = disp.from_XYZ(XYZ_grey3)
XYZ_inv_grid_disp = disp.to_XYZ(lab_grid_disp)
compare('forward neutral-corrected (neutral on), grid', impl.forward(XYZ_grid, P, neutral=NT),
        lab_grid_disp, XYZ_grid, 'XYZ', LAB)
compare('forward neutral-corrected (neutral on), grey Y in [0,3] (601)',
        impl.forward(XYZ_grey3, P, neutral=NT), lab_grey3_disp, XYZ_grey3, 'XYZ', LAB)
compare('inverse neutral-corrected (neutral on), grid (input = ref neutral-corrected forward) [extra]',
        impl.inverse(lab_grid_disp, P, neutral=NT), XYZ_inv_grid_disp,
        lab_grid_disp, 'Lab', XYZn)

# ------------------------------------------- extra end-to-end delta-E check
mine_grid = impl.forward(XYZ_grid, P)
compare('delta-E helmlab, same pairs, on impl base coords (end-to-end) [extra]',
        impl.delta_e(mine_grid[idx[:, 0]], mine_grid[idx[:, 1]], P), dE_ref,
        pairs, 'Lab1+Lab2', ('dE',))
xyz_pairs = np.concatenate([XYZ_grid[idx[:, 0]], XYZ_grid[idx[:, 1]]], axis=1)
compare('delta-E helmlab from XYZ: impl fwd+dE vs base.distance(XYZ_1, XYZ_2) [extra]',
        impl.delta_e(impl.forward(XYZ_grid[idx[:, 0]], P), impl.forward(XYZ_grid[idx[:, 1]], P), P),
        base.distance(XYZ_grid[idx[:, 0]], XYZ_grid[idx[:, 1]]), xyz_pairs, 'XYZ1+XYZ2', ('dE',))

# ----------------------------- extra: Section 5.6 sRGB -> XYZ from the text
M_SRGB = [[0.4124564390896922, 0.3575760776439511, 0.1804374832663989],
          [0.2126728514056226, 0.7151521552878178, 0.0721750036064596],
          [0.0193338955823293, 0.1191920258813418, 0.9503040785363679]]


def srgb_text(c):
    lin = lambda u: np.where(u <= 0.04045, u / 12.92, ((u + 0.055) / 1.055) ** 2.4)  # noqa: E731
    return np.stack(impl.matvec(M_SRGB, lin(c[:, 0]), lin(c[:, 1]), lin(c[:, 2])), axis=-1)


rgb_wide = np.concatenate([rgb, np.random.default_rng(3).uniform(-0.5, 1.5, (5000, 3))])
compare('Section 5.6 sRGB->XYZ from text vs sRGB_to_XYZ (grid + 5000 in [-0.5,1.5]^3) [extra]',
        srgb_text(rgb_wide), sRGB_to_XYZ(rgb_wide), rgb_wide, 'RGB', XYZn)

# ------------------- extra: Step N6 formulas vs scipy PchipInterpolator (named by the text)
from scipy.interpolate import PchipInterpolator  # noqa: E402
Lq = np.concatenate([np.linspace(0.0, NT.x[-1], 200001), NT.x])
mine_pc = np.stack([NT.Pa(Lq), NT.Pb(Lq)], axis=-1)
sp_pc = np.stack([PchipInterpolator(NT.x, NT.ya, extrapolate=False)(Lq),
                  PchipInterpolator(NT.x, NT.yb, extrapolate=False)(Lq)], axis=-1)
compare('Step N6 PCHIP from text vs scipy 1.17.1 PchipInterpolator on the same 380 knots [extra]',
        mine_pc, sp_pc, Lq[:, None], 'L', ('Pa', 'Pb'))

# ------------------------------------ checks of constants stated in the text
print('\n--- constants stated in the description ---')
M1i = np.linalg.inv(np.array(P['M1']))
M2i = np.linalg.inv(np.array(P['M2']))
print('M1_inv text == numpy.linalg.inv(M1):', np.array_equal(M1i, np.array(impl.M1_INV_DESC)),
      ' max diff', np.abs(M1i - np.array(impl.M1_INV_DESC)).max())
print('M2_inv text == numpy.linalg.inv(M2):', np.array_equal(M2i, np.array(impl.M2_INV_DESC)),
      ' max diff', np.abs(M2i - np.array(impl.M2_INV_DESC)).max())
phi, cr, sr = impl.rotation(P)
print('phi', repr(phi), phi.hex(), ' text -0.4921828490624009 -0x1.f7fec7e1801e7p-2',
      phi == -0.4921828490624009 and phi.hex() == '-0x1.f7fec7e1801e7p-2')
print('cr ', repr(cr), cr.hex(), ' text 0.8813034520649923 0x1.c33a34c0f12cep-1',
      cr == 0.8813034520649923 and cr.hex() == '0x1.c33a34c0f12cep-1')
print('sr ', repr(sr), sr.hex(), ' text -0.47255076486905395 -0x1.e3e459033fde9p-2',
      sr == -0.47255076486905395 and sr.hex() == '-0x1.e3e459033fde9p-2')
print('np.deg2rad(display_phi_deg) == phi:', np.deg2rad(P['display_phi_deg']) == phi)
bound = (1 / P['dist_compress']) ** P['dist_post_power']
print('dE bound', repr(bound), ' text 0.15003639078432002', bound == 0.15003639078432002)
print('max dE_ref', dE_ref.max())
print('neutral table: unique knots', NT.n_unique, '(text 255); appended', NT.n_ext,
      '(text 125); total', len(NT.x), '(text 380)')
print('  x_254', repr(NT.x254), '(text 1.3414218614893383)', NT.x254 == 1.3414218614893383)
print('  ya_254', repr(NT.ya254), '(text -0.0341747723961513)', NT.ya254 == -0.0341747723961513)
print('  yb_254', repr(NT.yb254), '(text 0.13436617089054925)', NT.yb254 == 0.13436617089054925)
print('  x_379', repr(NT.x[-1]), '(text 2.5914218614893394)', NT.x[-1] == 2.5914218614893394)
print('  first knot', NT.x[0], NT.ya[0], NT.yb[0])
Y_np = np.concatenate([[0.0], np.linspace(0.001, 0.01, 10), np.linspace(0.01, 0.1, 20),
                       np.linspace(0.1, 2.0, 226)])
print('  N1 linspace formula == numpy.linspace:', np.array_equal(Y_np, NT.Y))

# exact constants stated in the text, recorded as an extra quantity (tolerance 0)
computed = np.concatenate([M1i.ravel(), M2i.ravel(), [phi, cr, sr, bound, NT.x254, NT.ya254, NT.yb254,
                                                       NT.x[-1], NT.n_unique, NT.n_ext, len(NT.x)]])
stated = np.concatenate([np.ravel(impl.M1_INV_DESC), np.ravel(impl.M2_INV_DESC),
                         [-0.4921828490624009, 0.8813034520649923, -0.47255076486905395,
                          0.15003639078432002, 1.3414218614893383, -0.0341747723961513,
                          0.13436617089054925, 2.5914218614893394, 255, 125, 380]])
names = ([f'M1_inv[{i}][{j}]' for i in range(3) for j in range(3)]
         + [f'M2_inv[{i}][{j}]' for i in range(3) for j in range(3)]
         + ['phi', 'cr', 'sr', 'dE bound', 'x_254', 'ya_254', 'yb_254', 'x_379', 'n knots', 'n appended', 'n total'])
d = np.abs(computed - stated)
i = int(np.argmax(d))
results.append(dict(quantity='constants written in the text vs recomputed from params/description '
                    '(M1_inv, M2_inv via numpy.linalg.inv; phi, cr, sr; dE bound; N4/N5 table landmarks and counts) [extra]',
                    n=int(len(d)), max_abs_diff=float(d.max()), tolerance=0.0,
                    worst_input=('all bit-identical' if d.max() == 0 else
                                 f'{names[i]}: computed {computed[i]!r}, text {stated[i]!r}')))
print(f'[{"PASS" if d.max() == 0 else "FAIL"}] constants: n={len(d)} max_abs_diff={d.max():.3e}')

print('\nALL PASS' if all(r['max_abs_diff'] <= r['tolerance'] for r in results) else '\nSOME FAIL')
