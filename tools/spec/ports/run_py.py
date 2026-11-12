# Check ColorAide against the Helmlab 2026 vectors.
# Usage: python3 -I run_py.py <coloraide package dir> <vectors dir> <css|astm>
# Prints one JSON object in the same shape as run_js.mjs.
import json
import math
import sys

pkg, vec_dir, frame_arg = sys.argv[1:4]
sys.path.insert(0, pkg)
import coloraide  # noqa: E402
from coloraide import Color as Base  # noqa: E402
from coloraide.spaces.helmgen import Helmgen  # noqa: E402
from coloraide.spaces.helmgenlch import Helmgenlch  # noqa: E402
from coloraide.spaces.helmlab_metric import HelmlabMetric  # noqa: E402
from coloraide.distance.delta_e_helmlab import DEHelmlab  # noqa: E402


class Color(Base):
    pass


Color.register([Helmgen(), Helmgenlch(), HelmlabMetric(), DEHelmlab()])
frame = 'xyz-astm-d65' if frame_arg == 'astm' else 'xyz-css-d65'
TOL_U = {'fwd': 1e-10, 'inv': 1e-8, 'de': 1e-10}
SPACES = {'helmgen', 'helmgenlch', 'helmlab-metric'}


def num(v):
    return float(v)


def fwd(space, xyz):
    return Color('xyz-d65', xyz).convert(space).coords() if space in SPACES else None


def inv(space, lab):
    return Color(space, lab).convert('xyz-d65').coords() if space in SPACES else None


def err(op, got, exp):
    if got is None or len(got) != 3:
        return math.inf
    if op == 'helmgenlch.fwd':
        dL, dC = abs(got[0] - exp[0]), abs(got[1] - exp[1])
        if math.isfinite(got[2]) and math.isfinite(exp[2]):
            dh = math.remainder(math.radians(got[2] - exp[2]), 2 * math.pi)
            chord = abs(exp[1]) * abs(dh)
        else:
            # An undefined hue (NaN) is charged the largest possible chord, C * pi.
            chord = abs(exp[1]) * math.pi if math.isnan(got[2]) and math.isfinite(exp[2]) else math.inf
        e = max(dL, dC, chord)
        return math.inf if math.isnan(e) else e
    e = 0.0
    for g, x in zip(got, exp):
        d = abs(g - x)
        if math.isnan(d):
            return math.inf
        e = max(e, d)
    return e


res = {}


def add(op, level, e, rid):
    s = res.setdefault(op, {}).setdefault(level, {'n': 0, 'max': 0.0, 'worst': None, 'over': 0})
    s['n'] += 1
    if s['worst'] is None or e >= s['max']:
        s['max'], s['worst'] = e, rid
    tol = TOL_U['de'] if op == 'deltaE' else TOL_U['inv'] if op.endswith('.inv') else TOL_U['fwd']
    if level == 'U' and not e <= tol:
        s['over'] += 1


for name in ['core-1', 'conv-1', 'grey-1', 'edge-1', 'de-1', 'neutral-1']:
    rows = json.load(open(f'{vec_dir}/{name}.json'))['rows']
    for r in rows:
        if r['level'] not in ('U', 'R'):
            continue
        op = r['op']
        try:
            if op == 'deltaE':
                got = Color('xyz-d65', r['in']['xyz1'][frame]).delta_e(
                    Color('xyz-d65', r['in']['xyz2'][frame]), method='helmlab')
                e = abs(got - num(r['out']))
                e = math.inf if math.isnan(e) else e
            else:
                space = op.rsplit('.', 1)[0]
                if op.endswith('.fwd'):
                    got = fwd(space, r['in'][frame])
                    if got is None:
                        continue
                    e = err(op, [float(t) for t in got], [num(t) for t in r['out']])
                else:
                    got = inv(space, [num(t) for t in r['in']])
                    if got is None:
                        continue
                    e = err(op, [float(t) for t in got], [num(t) for t in r['out'][frame]])
        except Exception:
            e = math.inf
        add(op, r['level'], e, f"{name} {r['id']}")

out = {'implementation': f'ColorAide {coloraide.__version__}', 'frame': frame, 'results': res}
print(json.dumps(out, default=str).replace('Infinity', '"Infinity"'))
