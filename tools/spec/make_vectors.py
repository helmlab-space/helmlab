"""Generate the test vectors of Helmlab 2026 (spec/2026/spec.md, Section 9).

Usage, from the repository root:

    python3 -I tools/spec/make_vectors.py
        Writes spec/2026/vectors/*.json and rewrites spec/2026/SHA256SUMS: the
        two params lines first, then one line per vector file, sorted by path.
    python3 -I tools/spec/make_vectors.py --out DIR
        Writes the vector files and a SHA256SUMS of them ("<sha256>  <file>")
        into DIR. The repository is not touched.
    python3 -I tools/spec/make_vectors.py --core-table
        Prints core-1 as a Markdown table (id, XYZ, helmgen Lab, helmlab-metric
        Lab; every number as its Python repr) and writes nothing.

The expected outputs come from the Python reference package in src/helmlab,
one call with an input of shape (1, 3) per input (spec Appendix A.2). The
header of each file states the frames, the domain and level rules, how each
group of inputs was made, the excluded ids and the row counts. The report on
stdout gives, per file, the rows per op, level and domain, the excluded ids,
the size and the sha256, then the checks of the reference.

Before it computes anything the generator checks that src/helmlab equals tag
v1.0.0 apart from the version string in __init__.py and docstrings (a changed
Python file must have the same syntax tree as at v1.0.0 once docstrings are
removed; any other changed file stops the run), that src/helmlab has no
uncommitted change, and that both parameter files match
spec/2026/SHA256SUMS. It stops with status 1 if a check fails, and also if a
row inside the domain fails a level rule: an output that is not finite, a
one-row and a batch result that differ by more than the level R tolerance,
or, for an inverse row, a one-row forward of the expected XYZ that differs
from the row's input by more than 1e-8, or an output that moves by more than
1e-8 when one input channel changes by 1e-12. Nothing is written in either
case.
"""
import argparse
import ast
import decimal
import hashlib
import itertools
import json
import math
import platform
import re
import subprocess
import sys
import warnings
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True  # leave src/helmlab/__pycache__ as it is

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
SPEC = REPO / "spec" / "2026"
VECTOR_DIR = SPEC / "vectors"
SUMS = SPEC / "SHA256SUMS"
API_SOURCE = REPO / "packages" / "helmlab-js" / "tests" / "reference" / "reference-1.0.json"
GENERATOR = "tools/spec/make_vectors.py"

SPEC_REV = "Helmlab 2026, Spec Rev 1 (draft of 2026-11-12)"
VARIANT = "default: rotation -28.2 degrees, surround S = 0.5"

TOL = {"fwd": 1e-11, "inv": 1e-8, "de": 1e-11}  # level R, absolute, any channel
TRIP_TOL = 1e-8                                  # rule (c), inverse rows
NUDGE = 1e-12                                    # rule (d): change of one input channel
NUDGE_TOL = 1e-8                                 # rule (d): allowed change of the output
C0_MIN = 1e-4                                    # C0 rule, spec Section 5.5.2

# The ops in the order of the design table; rows of one input follow it.
OPS = ("helmgen.fwd", "helmgen.inv", "helmgenlch.fwd", "helmgenlch.inv",
       "helmlab-metric.fwd", "helmlab-metric.inv",
       "helmlab-metric-neutral.fwd", "helmlab-metric-neutral.inv", "deltaE")
OPS_SIX = OPS[:6]
OPS_FIVE = ("helmgen.fwd", "helmgen.inv", "helmgenlch.fwd",
            "helmlab-metric.fwd", "helmlab-metric.inv")
OPS_NEUTRAL = OPS[6:8]

CALLS = {
    "helmgen.fwd": "GenSpace().from_XYZ(xyz.reshape(1, 3))[0]",
    "helmgen.inv": "GenSpace().to_XYZ(lab.reshape(1, 3))[0]",
    "helmgenlch.fwd": "Helmlab().gen.to_lch(lab), lab the helmgen.fwd output (a plain ndarray)",
    "helmgenlch.inv": "GenSpace().to_XYZ(numpy.asarray(Helmlab().gen.from_lch(lch)).reshape(1, 3))[0]",
    "helmlab-metric.fwd": "MetricSpace().from_XYZ(xyz.reshape(1, 3))[0]",
    "helmlab-metric.inv": "MetricSpace().to_XYZ(lab.reshape(1, 3))[0]",
    "helmlab-metric-neutral.fwd": "MetricSpace(neutral_correction=True).from_XYZ(xyz.reshape(1, 3))[0]",
    "helmlab-metric-neutral.inv": "MetricSpace(neutral_correction=True).to_XYZ(lab.reshape(1, 3))[0]",
    "deltaE": "float(MetricSpace().distance_from_lab(lab1.reshape(1, 3), lab2.reshape(1, 3))[0])",
}

TOLERANCES = {
    "R": {"forward": 1e-11, "inverse": 1e-8, "deltaE": 1e-11},
    "U": {"forward": 1e-10, "inverse": 1e-8, "deltaE": 1e-10},
    "note": "absolute, per channel; forward: the coordinates, inverse: XYZ. Level R: the helmlab "
            "packages, XYZ in xyz-astm-d65. Level U: other implementations, XYZ (inputs of "
            "forward rows, outputs of inverse rows) in the frame of the implementation, "
            "xyz-astm-d65 or xyz-css-d65",
    "lch": "L and C as coordinates; hue as the chord C*|dh|, dh in radians wrapped to "
           "[-pi, pi], C the expected chroma, with the tolerance of C",
}

DOMAIN_RULES = (
    "helmgen and helmgenlch: in if 0 <= Y <= 10 and every component of M1*XYZ is >= 0 "
    "(M1 of params/gen_params.json, sign decided in exact rational arithmetic); "
    "helmlab-metric and helmlab-metric-neutral: in if the colour lies in the Display P3 or the "
    "Rec.2020 gamut, that is every component of (B*M)^-1 * XYZ is in "
    "[-1e-6, 1 + 1e-6] for M = lin_P3_to_XYZ or M = lin_2020_to_XYZ of CSS Color 4 (B of frames, "
    "exact rational arithmetic; the slack admits rounding, such as the sRGB white of spec "
    "Section 5.2, Y = 1.0000000102999; every sRGB colour of Section 5.2 lies in the Rec.2020 "
    "gamut, and some lie up to 1.1e-5 outside the Display P3 gamut); deltaE: in if both "
    "colours are in the MetricSpace domain; an inverse row has the flag of the forward row "
    "with the same id; XYZ and Y are the xyz-astm-d65 values; otherwise out."
)
LEVEL_RULES = (
    "A row in the domain has level U and is checked for every implementation (helmlab "
    "packages at level R, other implementations at level U); a row out of the domain has level R "
    "and is checked only for the helmlab packages, except level info (checked by no level) "
    "when (a) an output is not finite, (b) the one-row result and a batch result (one call "
    "with all inputs of the file for that op, shape (N, 3)) differ by more than the level R "
    "tolerance, (c) for an inverse row, the one-row forward of the expected XYZ differs "
    "from the row's input by more than 1e-8 (helmgenlch: L, C and the chord C*|dh|), or "
    "(d) for an inverse row, adding or subtracting 1e-12 to one channel of the input moves "
    "the output by more than 1e-8 (one batch call over the six changed inputs of each row)."
)

# Text that describes how XYZ was obtained (used in the "groups" headers).
CSS_SRC = ("CSS Color 4 conversions.js (Editor's Draft, https://drafts.csswg.org/css-color-4/conversions.js, "
           "downloaded 2026-11-09)")
TXT_GREY = "XYZ = (Y*0.95047, Y, Y*1.08883) in binary64 arithmetic"
TXT_SRGB = ("XYZ by spec Section 5.2: helmlab.utils.srgb_convert.sRGB_to_XYZ(numpy.array([[r, g, b]]))[0], "
            "one-row call")
TXT_P3 = ("XYZ by " + CSS_SRC + ": lin_sRGB of each channel with 50-digit decimal arithmetic, taken as "
          "an exact fraction, then the matrix lin_P3_to_XYZ and B (css to astm, see frames) in exact "
          "rational arithmetic, rounded once per component")
TXT_2020 = ("XYZ by " + CSS_SRC + ": lin_2020 (sign(v)*|v|**2.4, BT.1886 with b = 0) of each channel "
            "with 50-digit decimal arithmetic, taken as an exact fraction, then the matrix lin_2020_to_XYZ "
            "and B (css to astm, see frames) in exact rational arithmetic, rounded once per component")
TXT_2020_LIN = ("XYZ by " + CSS_SRC + ": the binary64 linear values taken as exact fractions, the "
                "matrix lin_2020_to_XYZ and B (css to astm, see frames) in exact rational arithmetic, "
                "rounded once per component")


def stop(message):
    print(f"make_vectors.py: {message}", file=sys.stderr)
    sys.exit(1)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    return sha256_bytes(Path(path).read_bytes())


# ---------------------------------------------------------------------------
# Exact 3x3 matrix arithmetic (Fractions)

def mat_vec(m, v):
    return [m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2] for i in range(3)]


def mat_mul(a, b):
    return [[a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] for j in range(3)]
            for i in range(3)]


def mat_inv(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    adj = [[e * i - f * h, c * h - b * i, b * f - c * e],
           [f * g - d * i, a * i - c * g, c * d - a * f],
           [d * h - e * g, b * g - a * h, a * e - b * d]]
    return [[x / det for x in row] for row in adj]


# ---------------------------------------------------------------------------
# Frames (design Section 3)

# Bradford matrix, Lam 1985 as tabulated by Lindbloom; decimal values, exact.
M_BFD = [[F("0.8951"), F("0.2664"), F("-0.1614")],
         [F("-0.7502"), F("1.7135"), F("0.0367")],
         [F("0.0389"), F("-0.0685"), F("1.0296")]]
W_ASTM = [F("0.95047"), F(1), F("1.08883")]
# CSS Color 4: [0.3127 / 0.3290, 1, (1 - 0.3127 - 0.3290) / 0.3290]
W_CSS = [F(3127, 3290), F(1), F(3583, 3290)]


def bradford(w_from, w_to):
    r_from = mat_vec(M_BFD, w_from)
    r_to = mat_vec(M_BFD, w_to)
    scale = [[r_to[i] / r_from[i] if i == j else F(0) for j in range(3)] for i in range(3)]
    return mat_mul(mat_inv(M_BFD), mat_mul(scale, M_BFD))


A = bradford(W_ASTM, W_CSS)  # xyz-astm-d65 -> xyz-css-d65
B = mat_inv(A)               # xyz-css-d65 -> xyz-astm-d65
assert mat_vec(A, W_ASTM) == W_CSS and mat_vec(B, W_CSS) == W_ASTM
A_FLOAT = [[float(x) for x in row] for row in A]


def css_of(xyz):
    """xyz-css-d65 value of an xyz-astm-d65 binary64 triple: A*xyz exactly, rounded once."""
    if all(math.isfinite(v) for v in xyz):
        return [float(c) for c in mat_vec(A, [F(v) for v in xyz])]
    # Only for outputs of info rows: binary64 arithmetic, NaN and infinities propagate.
    return [A_FLOAT[i][0] * xyz[0] + A_FLOAT[i][1] * xyz[1] + A_FLOAT[i][2] * xyz[2]
            for i in range(3)]


# ---------------------------------------------------------------------------
# CSS Color 4 input conversions (design Section 4)

DEC = decimal.Context(prec=50)
D = decimal.Decimal

M_P3 = [[F(608311, 1250200), F(189793, 714400), F(198249, 1000160)],
        [F(35783, 156275), F(247089, 357200), F(198249, 2500400)],
        [F(0), F(32229, 714400), F(5220557, 5000800)]]
M_2020 = [[F(63426534, 99577255), F(20160776, 139408157), F(47086771, 278816314)],
          [F(26158966, 99577255), F(472592308, 697040785), F(8267143, 139408157)],
          [F(0), F(19567812, 697040785), F(295819943, 278816314)]]
assert mat_vec(M_P3, [F(1)] * 3) == W_CSS and mat_vec(M_2020, [F(1)] * 3) == W_CSS
BM_P3 = mat_mul(B, M_P3)      # linear-light Display P3 -> xyz-astm-d65, exact
BM_2020 = mat_mul(B, M_2020)  # linear Rec.2020 -> xyz-astm-d65, exact
# MetricSpace domain: inside the Display P3 or the Rec.2020 gamut. The sRGB colours of spec
# Section 5.2 lie inside Rec.2020; towards blue and towards yellow some lie up to 1.1e-5
# outside Display P3, because Section 5.2 has its own matrix.
LIN_P3 = mat_inv(BM_P3)       # xyz-astm-d65 -> linear Display P3, exact
LIN_2020 = mat_inv(BM_2020)   # xyz-astm-d65 -> linear Rec.2020, exact
GAMUT_TOL = F(1, 10 ** 6)     # slack for rounding, e.g. the sRGB white of Section 5.2


def lin_srgb_css(v):
    """CSS Color 4 lin_sRGB (also lin_P3) of a binary64 value, 50 digits, as a Fraction."""
    d = D(v)                 # exact value of the binary64 input
    a = d.copy_abs()         # copy_abs and copy_negate do not round
    if a <= D("0.04045"):
        r = DEC.divide(d, D("12.92"))
    else:
        r = DEC.power(DEC.divide(DEC.add(a, D("0.055")), D("1.055")), D("2.4"))
        if d < 0:
            r = r.copy_negate()
    return F(r)


def lin_2020_css(v):
    """CSS Color 4 lin_2020: sign(v)*|v|**2.4, 50 digits, as a Fraction."""
    d = D(v)
    r = DEC.power(d.copy_abs(), D("2.4"))
    return F(r.copy_negate() if d < 0 else r)


def xyz_from_linear(bm, lin):
    """bm*lin in exact arithmetic, rounded once per component."""
    return [float(c) for c in mat_vec(bm, lin)]


def p3_xyz(rgb):
    return xyz_from_linear(BM_P3, [lin_srgb_css(v) for v in rgb])


def rec2020_xyz(rgb):
    return xyz_from_linear(BM_2020, [lin_2020_css(v) for v in rgb])


def rec2020_linear_xyz(rgb):
    return xyz_from_linear(BM_2020, [F(float(v)) for v in rgb])


def grey_xyz(y):
    return [y * 0.95047, y, y * 1.08883]


# ---------------------------------------------------------------------------
# Source checks (design Section 1)

def git(*args, check=True):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=check)


def syntax_without_docstrings(source):
    """ast.dump of a module with the docstrings of the module, its classes and
    its functions removed. Comments and layout are not part of the tree."""
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:] or [ast.Pass()]
    return ast.dump(tree)


def check_source():
    r = git("diff", "--name-only", "v1.0.0", "--", "src/helmlab", ":(exclude)src/helmlab/__init__.py")
    docstring_only = []
    for path in r.stdout.split():
        if not path.endswith(".py") or not (REPO / path).is_file():
            stop(f"{path} differs from tag v1.0.0")
        old = git("show", f"v1.0.0:{path}", check=False)
        if old.returncode != 0:
            stop(f"{path} is not in tag v1.0.0")
        if syntax_without_docstrings(old.stdout) != syntax_without_docstrings((REPO / path).read_text()):
            stop(f"{path} differs from tag v1.0.0 in more than docstrings")
        docstring_only.append(path)
    r = git("status", "--porcelain", "--untracked-files=all", "--", "src/helmlab")
    changed = [line for line in r.stdout.splitlines() if "__pycache__" not in line]
    if changed:
        stop("uncommitted changes under src/helmlab: " + "; ".join(changed))
    lines = SUMS.read_text().splitlines()
    listed = {}
    for line in lines:
        sha, _, path = line.partition("  ")
        listed[path] = sha
    params_lines = [line for line in lines if line.partition("  ")[2].startswith("params/")]
    sha = {}
    for name in ("gen_params.json", "metric_params.json"):
        want = listed.get("params/" + name)
        if want is None:
            stop(f"{SUMS} has no line for params/{name}")
        for path in (SRC / "helmlab" / "data" / name, SPEC / "params" / name):
            if sha256_file(path) != want:
                stop(f"sha256 of {path} differs from SHA256SUMS")
        sha[name] = want
    if len(params_lines) != 2:
        stop(f"{SUMS} must hold exactly two params lines")
    init = (SRC / "helmlab" / "__init__.py").read_text()
    m = re.search(r'__version__ = "([^"]+)"\s+# bump-version-fallback', init)
    if not m:
        stop("src/helmlab/__init__.py has no fallback version string")
    old = git("show", "v1.0.0:src/helmlab/__init__.py").stdout
    if old.replace('__version__ = "1.0.0"', f'__version__ = "{m.group(1)}"', 1) != init:
        stop("src/helmlab/__init__.py differs from tag v1.0.0 in more than the version string")
    return {
        "head": git("rev-parse", "HEAD").stdout.strip(),
        # The tree id names the content of src/helmlab (clean, checked above), so the
        # files do not change when unrelated commits move HEAD.
        "tree": git("rev-parse", "HEAD:src/helmlab").stdout.strip(),
        "tag": git("rev-parse", "v1.0.0^{commit}").stdout.strip(),
        "version": m.group(1),
        "docstring_only": docstring_only,
        "params_sha256": sha,
        "params_lines": params_lines,
    }


# ---------------------------------------------------------------------------
# Reference calls (design Section 2)

class Reference:
    """The reference package, its objects built once, and the calls of each op."""

    def __init__(self):
        sys.path.insert(0, str(SRC))
        import helmlab
        from helmlab.helmlab import Helmlab
        from helmlab.spaces.gen import GenSpace
        from helmlab.spaces.metric import MetricSpace
        from helmlab.utils.srgb_convert import sRGB_to_XYZ
        self.file = Path(helmlab.__file__).resolve()
        if SRC.resolve() not in self.file.parents:
            stop(f"helmlab was loaded from {self.file}, not from {SRC}")
        # helmlab.__version__ comes from package metadata; a stale src/helmlab.egg-info
        # can make it differ from the string in __init__.py. It is only reported.
        self.version_attr = helmlab.__version__
        self.gen = GenSpace()
        self.metric = MetricSpace()
        self.neutral = MetricSpace(neutral_correction=True)
        self.hl = Helmlab()
        self._srgb_to_xyz = sRGB_to_XYZ

    def srgb_xyz(self, rgb):
        return [float(v) for v in self._srgb_to_xyz(np.array([rgb], dtype=np.float64))[0]]

    def one(self, op, v):
        """One-row call of op on three numbers; returns three floats."""
        x = np.array(v, dtype=np.float64).reshape(1, 3)
        if op == "helmgen.fwd":
            r = self.gen.from_XYZ(x)[0]
        elif op == "helmgen.inv":
            r = self.gen.to_XYZ(x)[0]
        elif op == "helmgenlch.fwd":
            r = self.hl.gen.to_lch(self.gen.from_XYZ(x)[0])
        elif op == "helmgenlch.inv":
            lab = np.asarray(self.hl.gen.from_lch(x[0]))
            r = self.gen.to_XYZ(lab.reshape(1, 3))[0]
        elif op == "helmlab-metric.fwd":
            r = self.metric.from_XYZ(x)[0]
        elif op == "helmlab-metric.inv":
            r = self.metric.to_XYZ(x)[0]
        elif op == "helmlab-metric-neutral.fwd":
            r = self.neutral.from_XYZ(x)[0]
        elif op == "helmlab-metric-neutral.inv":
            r = self.neutral.to_XYZ(x)[0]
        else:
            raise ValueError(op)
        return [float(t) for t in r]

    def lch(self, lab):
        return [float(t) for t in self.hl.gen.to_lch(np.array(lab, dtype=np.float64))]

    def delta_e(self, lab1, lab2):
        a = np.array(lab1, dtype=np.float64).reshape(1, 3)
        b = np.array(lab2, dtype=np.float64).reshape(1, 3)
        return float(self.metric.distance_from_lab(a, b)[0])

    def distance(self, xyz1, xyz2):
        a = np.array(xyz1, dtype=np.float64).reshape(1, 3)
        b = np.array(xyz2, dtype=np.float64).reshape(1, 3)
        return float(self.metric.distance(a, b)[0])

    def batch(self, op, inputs):
        """One call with all inputs at once, shape (N, 3). helmgenlch applies to_lch and
        from_lch (which take one colour) row by row around the batch GenSpace call."""
        x = np.array(inputs, dtype=np.float64).reshape(-1, 3)
        if op == "helmgen.fwd":
            r = self.gen.from_XYZ(x)
        elif op == "helmgen.inv":
            r = self.gen.to_XYZ(x)
        elif op == "helmgenlch.fwd":
            r = [self.hl.gen.to_lch(lab) for lab in self.gen.from_XYZ(x)]
        elif op == "helmgenlch.inv":
            r = self.gen.to_XYZ(np.array([np.asarray(self.hl.gen.from_lch(t)) for t in x]))
        elif op == "helmlab-metric.fwd":
            r = self.metric.from_XYZ(x)
        elif op == "helmlab-metric.inv":
            r = self.metric.to_XYZ(x)
        elif op == "helmlab-metric-neutral.fwd":
            r = self.neutral.from_XYZ(x)
        elif op == "helmlab-metric-neutral.inv":
            r = self.neutral.to_XYZ(x)
        else:
            raise ValueError(op)
        return [[float(t) for t in row] for row in r]

    def batch_delta_e(self, labs1, labs2):
        a = np.array(labs1, dtype=np.float64).reshape(-1, 3)
        b = np.array(labs2, dtype=np.float64).reshape(-1, 3)
        return [float(t) for t in self.metric.distance_from_lab(a, b)]


# ---------------------------------------------------------------------------
# Comparisons

def abs_diff(got, exp):
    """Largest absolute channel difference; inf where only one side is finite or they differ."""
    worst = 0.0
    for g, e in zip(got, exp):
        if math.isfinite(g) and math.isfinite(e):
            worst = max(worst, abs(g - e))
        elif not (g == e or (math.isnan(g) and math.isnan(e))):
            return math.inf
    return worst


def lch_diff(got, exp):
    """L, C and the hue chord C*|dh| (spec Section 8); C is the expected chroma."""
    if not all(math.isfinite(v) for v in list(got) + list(exp)):
        return abs_diff(got, exp)
    dh = math.remainder(math.radians(got[2] - exp[2]), 2 * math.pi)
    return max(abs(got[0] - exp[0]), abs(got[1] - exp[1]), exp[1] * abs(dh))


def diff(op, got, exp):
    if op == "helmgenlch.fwd":
        return lch_diff(got, exp)
    if op == "deltaE":
        return abs_diff([got], [exp])
    return abs_diff(got, exp)


# ---------------------------------------------------------------------------
# Inputs and sets (design Sections 4, 5, 6 and 9)

class Builder:
    def __init__(self, ref, gen_params, metric_params):
        self.ref = ref
        self.m1_gen = [[F(v) for v in row] for row in gen_params["M1"]]
        self.mp = metric_params

    # Rules of one input --------------------------------------------------
    def gen_domain(self, xyz):
        if not 0 <= xyz[1] <= 10:
            return "out"
        if self.lms_negative(xyz):
            return "out"
        return "in"

    def lms_negative(self, xyz):
        x = [F(float(v)) for v in xyz]
        return any(row[0] * x[0] + row[1] * x[1] + row[2] * x[2] < 0 for row in self.m1_gen)

    def metric_domain(self, xyz):
        """In the Display P3 or the Rec.2020 gamut: every linear value of one of the two
        spaces within [-1e-6, 1 + 1e-6]."""
        x = [F(float(v)) for v in xyz]
        for m in (LIN_P3, LIN_2020):
            if all(-GAMUT_TOL <= v <= 1 + GAMUT_TOL for v in mat_vec(m, x)):
                return "in"
        return "out"

    def c0(self, xyz):
        """C0 = hypot(a0, b0) after MetricSpace Steps F1 to F3, binary64."""
        m1, m2, g = self.mp["M1"], self.mp["M2"], self.mp["gamma"]
        lms = [m1[i][0] * xyz[0] + m1[i][1] * xyz[1] + m1[i][2] * xyz[2] for i in range(3)]
        lc = [math.copysign(abs(v) ** p, v) if v != 0 else 0.0 for v, p in zip(lms, g)]
        a0 = m2[1][0] * lc[0] + m2[1][1] * lc[1] + m2[1][2] * lc[2]
        b0 = m2[2][0] * lc[0] + m2[2][1] * lc[1] + m2[2][2] * lc[2]
        return math.hypot(a0, b0)

    def make_input(self, id_, src, xyz):
        xyz = [float(v) for v in xyz]
        return {"id": id_, "src": src, "xyz": xyz, "css": css_of(xyz),
                "gen_domain": self.gen_domain(xyz),
                "metric_domain": self.metric_domain(xyz),
                "c0": self.c0(xyz),
                "black": xyz == [0.0, 0.0, 0.0]}

    # Groups ----------------------------------------------------------------
    def neg_lms(self):
        """Section 5.3.5 set: rows of default_rng(2026).uniform(-0.05, 1.2, (20000, 3))
        with a negative component of M1_gen*XYZ (exact sign), in row order."""
        x = np.random.default_rng(2026).uniform(-0.05, 1.2, (20000, 3))
        hits = [i for i in range(len(x)) if self.lms_negative(x[i])]
        return [[float(v) for v in x[i]] for i in hits], len(x)

    def core(self, neg_first):
        made = {}
        items = []

        def add(id_, src, xyz, text):
            items.append(self.make_input(id_, src, xyz))
            made[id_] = text

        for id_, y in (("core/black", 0.0), ("core/white", 1.0), ("core/grey-0.18", 0.18),
                       ("core/grey-0.5", 0.5)):
            add(id_, {"grey": y}, grey_xyz(y), f"grey, Y = {y!r}")
        primaries = (("red", [1.0, 0.0, 0.0]), ("green", [0.0, 1.0, 0.0]), ("blue", [0.0, 0.0, 1.0]))
        secondaries = (("cyan", [0.0, 1.0, 1.0]), ("magenta", [1.0, 0.0, 1.0]), ("yellow", [1.0, 1.0, 0.0]))
        for name, rgb in primaries + secondaries:
            add(f"core/srgb-{name}", {"srgb": rgb}, self.ref.srgb_xyz(rgb), f"encoded sRGB {rgb}")
        for name, rgb in primaries:
            add(f"core/p3-{name}", {"display-p3": rgb}, p3_xyz(rgb), f"encoded Display P3 {rgb}")
        for name, rgb in primaries:
            add(f"core/rec2020-{name}", {"rec2020": rgb}, rec2020_xyz(rgb), f"encoded Rec.2020 {rgb}")
        add("core/negative-xyz", None, [-0.05, -0.05, -0.05],
            "XYZ (-0.05, -0.05, -0.05), spec Section 5.3.5")
        add("core/negative-lms", None, neg_first, "the input of neg-lms/0 of edge-1")
        for id_, y in (("core/grey-2", 2.0), ("core/grey-10", 10.0)):
            add(id_, {"grey": y}, grey_xyz(y), f"grey, Y = {y!r}")
        red = next(t["xyz"] for t in items if t["id"] == "core/srgb-red")
        add("core/srgb-red-x2", None, [2.0 * v for v in red],
            "2 times the XYZ of core/srgb-red (exact doubling)")
        groups = {"core": {"count": len(items), "inputs": made,
                           "conversions": {"grey": TXT_GREY, "srgb": TXT_SRGB,
                                           "display-p3": TXT_P3, "rec2020": TXT_2020}}}
        return items, groups

    def conv(self):
        levels = [k / 8 for k in range(9)]
        items = []
        groups = {}
        for prefix, convert, text in (("srgb", self.ref.srgb_xyz, TXT_SRGB),
                                      ("display-p3", p3_xyz, TXT_P3),
                                      ("rec2020", rec2020_xyz, TXT_2020)):
            for i, j, k in itertools.product(range(9), repeat=3):
                rgb = [levels[i], levels[j], levels[k]]
                items.append(self.make_input(f"{prefix}/{i}-{j}-{k}", {prefix: rgb}, convert(rgb)))
            groups[prefix] = {"count": 729, "made": (
                f"9x9x9 grid of encoded values (i/8, j/8, k/8), i, j, k = 0..8, id {prefix}/i-j-k; "
                + text)}
        return items, groups

    def grey(self):
        items = []
        for k in range(1001):
            y = float(F(k * k, 62500))
            items.append(self.make_input(f"grey/{k}", {"grey": y}, grey_xyz(y)))
        groups = {"grey": {"count": 1001, "made": (
            "Y_k = float(Fraction(k*k, 62500)) for k = 0..1000 (Y from 0 to 16, Y = 1 at k = 250), "
            "id grey/k; " + TXT_GREY)}}
        return items, groups

    def edge(self, neg_rows, neg_draws):
        items = []
        groups = {}
        # neg-lms: Section 5.3.5 set
        for n, xyz in enumerate(neg_rows[:256]):
            items.append(self.make_input(f"neg-lms/{n}", None, xyz))
        groups["neg-lms"] = {"count": 256, "seed": 2026, "made": (
            "X = numpy.random.default_rng(2026).uniform(-0.05, 1.2, (20000, 3)) (spec Section 5.3.5); "
            "in row order, the first 256 rows with a negative component of M1*XYZ (M1 of "
            "params/gen_params.json, sign decided in exact rational arithmetic); XYZ as drawn"),
            "rows_with_negative_lms": len(neg_rows), "rows_drawn": neg_draws}
        # imaginary: negative X, Y or Z
        rng = np.random.default_rng(2027)
        draws = 0
        kept = 0
        while kept < 256:
            x = rng.uniform(-0.2, 1.2, 3)
            draws += 1
            if x[0] < 0 or x[1] < 0 or x[2] < 0:
                items.append(self.make_input(f"imaginary/{kept}", None, x))
                kept += 1
        groups["imaginary"] = {"count": 256, "seed": 2027, "made": (
            "rng = numpy.random.default_rng(2027); draws x = rng.uniform(-0.2, 1.2, 3), kept if "
            "X, Y or Z < 0, until 256 are kept; XYZ as drawn"), "draws": draws}
        # above-white: Rec.2020 linear s*u with Y > 1
        rng = np.random.default_rng(2028)
        draws = 0
        kept = 0
        while kept < 256:
            u = rng.uniform(0.0, 1.0, 3)
            s = rng.uniform(1.0, 16.0)
            draws += 1
            rgb = [float(v) for v in s * u]
            xyz = rec2020_linear_xyz(rgb)
            if xyz[1] > 1:
                items.append(self.make_input(f"above-white/{kept}", {"rec2020-linear": rgb}, xyz))
                kept += 1
        groups["above-white"] = {"count": 256, "seed": 2028, "made": (
            "rng = numpy.random.default_rng(2028); per draw u = rng.uniform(0, 1, 3), then "
            "s = rng.uniform(1, 16); Rec.2020 linear rgb = s*u in binary64; kept if Y > 1, until "
            "256 are kept; " + TXT_2020_LIN), "draws": draws}
        # rec2020-magenta: encoded magentas, scaled in linear light
        rng = np.random.default_rng(2029)
        for i in range(64):
            rgb = [float(v) for v in rng.uniform([0.75, 0.0, 0.75], [1.0, 0.25, 1.0])]
            lin = [lin_2020_css(v) for v in rgb]
            for s in (1, 2, 4, 8):
                xyz = xyz_from_linear(BM_2020, [s * t for t in lin])
                src = {"rec2020-scaled": rgb, "linear_scale": float(s)}
                items.append(self.make_input(f"rec2020-magenta/{i}-x{s}", src, xyz))
        groups["rec2020-magenta"] = {"count": 256, "seed": 2029, "made": (
            "rng = numpy.random.default_rng(2029); 64 encoded Rec.2020 colours, each "
            "(r, g, b) = rng.uniform([0.75, 0, 0.75], [1, 0.25, 1]); each linearized as in "
            "CSS Color 4 lin_2020 and multiplied by s = 1, 2, 4 and 8 in linear light (exact); "
            "id rec2020-magenta/i-xs; src gives the encoded colour (key rec2020-scaled) and "
            "s (key linear_scale): the XYZ is that of s times the linear colour; " + TXT_2020)}
        return items, groups

    def de(self):
        """de-1 pairs (design Section 9), all draws from one default_rng(2030)."""
        rng = np.random.default_rng(2030)
        redraws = {"de-near": 0, "de-far": 0}

        def colour(group):
            while True:
                u = [float(v) for v in rng.uniform(0.0, 1.0, 3)]
                xyz = rec2020_linear_xyz(u)
                if self.c0(xyz) >= C0_MIN:
                    return u, xyz
                redraws[group] += 1

        pairs = []
        for n in range(1000):
            u1, x1 = colour("de-near")
            while True:
                v = rng.standard_normal(3)
                d = v / math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
                e = rng.uniform(-4.0, -1.0)
                u2 = [float(t) for t in np.clip(np.array(u1) + 10.0 ** e * d, 0.0, 1.0)]
                x2 = rec2020_linear_xyz(u2)
                if self.c0(x2) >= C0_MIN:
                    break
                redraws["de-near"] += 1
            pairs.append((f"de-near/{n}", x1, x2))
        for n in range(1000):
            _, x1 = colour("de-far")
            _, x2 = colour("de-far")
            pairs.append((f"de-far/{n}", x1, x2))
        common = ("A colour is a Rec.2020 linear u = rng.uniform(0, 1, 3) (so Y <= 1), redrawn while "
                  "C0 < 1e-4; " + TXT_2020_LIN + ". All draws come from one "
                  "numpy.random.default_rng(2030): de-near/0..999 first, then de-far/0..999. "
                  "lab1 and lab2 are helmlab-metric.fwd outputs (one-row), deltaE is "
                  "distance_from_lab (one-row).")
        groups = {
            "de-near": {"count": 1000, "seed": 2030, "made": (
                common + " de-near/n: colour 1 is u; colour 2 is clip(u + 10**e*d, 0, 1) with "
                "d = v/|v|, v = rng.standard_normal(3), |v| = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) "
                "in binary64, then e = rng.uniform(-4, -1); v and e are redrawn while colour 2 has "
                "C0 < 1e-4."), "redraws": redraws["de-near"]},
            "de-far": {"count": 1000, "seed": 2030, "made": (
                common + " de-far/n: two colours drawn one after the other."),
                "redraws": redraws["de-far"]},
        }
        return pairs, groups


def build_rows(inputs, ops, ref):
    """Rows of a set: per input, the ops in table order (design Section 8)."""
    rows = []
    excluded = []
    black_ok = True
    uses_metric = any(op.startswith("helmlab-metric") for op in ops)
    for inp in inputs:
        metric_rows = inp["black"] or inp["c0"] >= C0_MIN
        if uses_metric and not metric_rows:
            excluded.append({"id": inp["id"], "C0": inp["c0"]})
        fwd = {}
        for op in ops:
            space, kind = op.split(".")
            if space.startswith("helmlab-metric") and not metric_rows:
                continue
            domain = inp["gen_domain"] if space.startswith("helmgen") else inp["metric_domain"]
            if kind == "fwd":
                if op == "helmgenlch.fwd":
                    out = ref.lch(fwd["helmgen"])
                else:
                    out = ref.one(op, inp["xyz"])
                fwd[space] = out
                if space.startswith("helmlab-metric") and inp["black"] and out != [0.0, 0.0, 0.0]:
                    black_ok = False
                rows.append({"id": inp["id"], "op": op, "kind": "fwd", "src": inp["src"],
                             "x": inp["xyz"], "x_css": inp["css"], "out": out, "domain": domain})
            else:
                x = fwd[space]
                rows.append({"id": inp["id"], "op": op, "kind": "inv", "x": x,
                             "out": ref.one(op, x), "domain": domain, "orig": inp["xyz"]})
    return rows, excluded, black_ok


def build_de_rows(pairs, ref, builder):
    rows = []
    same = 0
    for id_, x1, x2 in pairs:
        lab1 = ref.one("helmlab-metric.fwd", x1)
        lab2 = ref.one("helmlab-metric.fwd", x2)
        de = ref.delta_e(lab1, lab2)
        if ref.distance(x1, x2).hex() == de.hex():
            same += 1
        inside = all(builder.metric_domain(x) == "in" for x in (x1, x2))
        rows.append({"id": id_, "op": "deltaE", "kind": "de", "lab1": lab1, "lab2": lab2,
                     "xyz1": x1, "xyz2": x2, "out": de, "domain": "in" if inside else "out"})
    return rows, same


def assign_levels(rows, ref):
    """Level of each row (design Section 7). Returns the domain-in rows that fail a rule."""
    failures = []
    for op in OPS:
        rs = [r for r in rows if r["op"] == op]
        if not rs:
            continue
        if op == "deltaE":
            batch = ref.batch_delta_e([r["lab1"] for r in rs], [r["lab2"] for r in rs])
        else:
            batch = ref.batch(op, [r["x"] for r in rs])
        if op.endswith(".inv"):
            # Rule (d): the six inputs of each row with one channel changed by +-NUDGE.
            # Where the inverse is this sensitive, the expected XYZ depends on rounding
            # inside the implementation more than on the input.
            nudged = []
            for r in rs:
                for i in range(3):
                    for step in (NUDGE, -NUDGE):
                        x = list(r["x"])
                        x[i] = x[i] + step
                        nudged.append(x)
            outs = ref.batch(op, nudged)
            for k, r in enumerate(rs):
                r["d_nudge"] = max(abs_diff(o, r["out"]) for o in outs[6 * k:6 * k + 6])
        for r, b in zip(rs, batch):
            reasons = []
            values = r["out"] if isinstance(r["out"], list) else [r["out"]]
            if not all(math.isfinite(v) for v in values):
                reasons.append("(a) output not finite")
            r["d_batch"] = diff(op, b, r["out"])
            if r["d_batch"] > TOL[r["kind"]]:
                reasons.append(f"(b) one-row vs batch {r['d_batch']!r}")
            if r["kind"] == "inv":
                fop = op[:-4] + ".fwd"
                back = ref.one(fop, r["out"])
                r["d_trip"] = diff(fop, back, r["x"])
                if not r["d_trip"] <= TRIP_TOL:
                    reasons.append(f"(c) forward of the expected XYZ differs by {r['d_trip']!r}")
                if not r["d_nudge"] <= NUDGE_TOL:
                    reasons.append(f"(d) a change of 1e-12 in the input moves the output by "
                                   f"{r['d_nudge']!r}")
                r["d_xyz"] = abs_diff(r["out"], r["orig"])
            r["reasons"] = reasons
            if r["domain"] == "in":
                r["level"] = "U"
                if reasons:
                    failures.append(r)
            else:
                r["level"] = "info" if reasons else "R"
    return failures


# ---------------------------------------------------------------------------
# Files (design Section 8)

def num(v):
    v = float(v)
    if math.isfinite(v):
        return v
    return "NaN" if math.isnan(v) else ("Infinity" if v > 0 else "-Infinity")


def nums(values):
    return [num(v) for v in values]


def row_json(r):
    if r["kind"] == "fwd":
        d = {"id": r["id"], "op": r["op"]}
        if r["src"] is not None:
            d["src"] = r["src"]
        d["in"] = {"xyz-astm-d65": nums(r["x"]), "xyz-css-d65": nums(r["x_css"])}
        d["out"] = nums(r["out"])
        values = r["out"]
    elif r["kind"] == "inv":
        d = {"id": r["id"], "op": r["op"], "in": nums(r["x"]),
             "out": {"xyz-astm-d65": nums(r["out"]), "xyz-css-d65": nums(css_of(r["out"]))}}
        values = list(r["x"]) + list(r["out"])
    else:
        d = {"id": r["id"], "op": "deltaE",
             "in": {"lab1": nums(r["lab1"]), "lab2": nums(r["lab2"]),
                    "xyz1": {"xyz-astm-d65": nums(r["xyz1"]), "xyz-css-d65": nums(css_of(r["xyz1"]))},
                    "xyz2": {"xyz-astm-d65": nums(r["xyz2"]), "xyz-css-d65": nums(css_of(r["xyz2"]))}},
             "out": num(r["out"])}
        values = [r["out"]]
    if r["level"] != "info" and not all(math.isfinite(v) for v in values):
        raise AssertionError(f"non-finite number in a {r['level']} row: {r['id']} {r['op']}")
    d["domain"] = r["domain"]
    d["level"] = r["level"]
    return json.dumps(d, separators=(",", ":"), allow_nan=False)


def counts_of(rows):
    counts = {}
    for op in OPS:
        rs = [r for r in rows if r["op"] == op]
        if not rs:
            continue
        c = {}
        for level, domain in (("U", "in"), ("R", "out"), ("info", "out")):
            n = sum(1 for r in rs if r["level"] == level and r["domain"] == domain)
            if n:
                c.setdefault(level, {})[domain] = n
        counts[op] = c
    counts["total"] = len(rows)
    return counts


def file_text(header, rows):
    head = json.dumps(header, indent=1, allow_nan=False)
    assert head.endswith("\n}")
    body = ",\n".join(row_json(r) for r in rows)
    return head[:-2] + ',\n "rows": [\n' + body + "\n]}\n"


def environment():
    import scipy
    cfg = np.show_config(mode="dicts")
    deps = cfg.get("Build Dependencies", {})
    return {
        "python": sys.version,
        "numpy": np.__version__,
        "numpy_blas": deps.get("blas"),
        "numpy_lapack": deps.get("lapack"),
        "numpy_simd": cfg.get("SIMD Extensions"),
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "mac_ver": platform.mac_ver(),
    }


def frames_header():
    return {
        "xyz-astm-d65": {"white": [0.95047, 1.0, 1.08883],
                         "note": "reference frame; the reference calls take these values"},
        "xyz-css-d65": {"white": [float(v) for v in W_CSS],
                        "white_exact": "(3127/3290, 1, 3583/3290)"},
        "A": A_FLOAT,
        "A_definition": (
            "Bradford adaptation from the astm white to the css white, "
            "A = M_BFD^-1 * diag((M_BFD*W_css)_i / (M_BFD*W_astm)_i) * M_BFD with "
            "M_BFD = [[0.8951, 0.2664, -0.1614], [-0.7502, 1.7135, 0.0367], "
            "[0.0389, -0.0685, 1.0296]] (Lam 1985 as tabulated by Lindbloom), computed in exact "
            "rational arithmetic; A*W_astm = W_css exactly; the floats above are A rounded; "
            "B = A^-1 (exact) maps xyz-css-d65 to xyz-astm-d65"),
        "rule": (
            "Every xyz-css-d65 value in this file is A times the xyz-astm-d65 binary64 value, "
            "computed exactly and rounded once per component, and the expected outputs are "
            "computed from the xyz-astm-d65 values only."),
    }


# ---------------------------------------------------------------------------
# Report (design Section 10)

def fmt(x):
    return "-" if x is None else f"{x:.3g}"


def file_report(name, text, rows, excluded):
    data = text.encode("ascii")
    out = [f"{name}: {len(rows)} rows, {len(data)} bytes, sha256 {sha256_bytes(data)}"]
    out.append(f"  {'op':28s} {'U/in':>6s} {'R/out':>6s} {'info/out':>8s}   "
               f"{'max |1-row - batch| (in)':>24s} {'max round trip (c) (in)':>23s} "
               f"{'max |xyz rt| (in)':>17s} {'max change (d) (in)':>19s}")
    for op in OPS:
        rs = [r for r in rows if r["op"] == op]
        if not rs:
            continue
        n_u = sum(1 for r in rs if r["level"] == "U")
        n_r = sum(1 for r in rs if r["level"] == "R")
        n_i = sum(1 for r in rs if r["level"] == "info")
        inside = [r for r in rs if r["domain"] == "in"]
        b = max((r["d_batch"] for r in inside), default=None)
        c = max((r["d_trip"] for r in inside if "d_trip" in r), default=None)
        x = max((r["d_xyz"] for r in inside if "d_xyz" in r), default=None)
        n = max((r["d_nudge"] for r in inside if "d_nudge" in r), default=None)
        out.append(f"  {op:28s} {n_u:6d} {n_r:6d} {n_i:8d}   {fmt(b):>24s} {fmt(c):>23s} "
                   f"{fmt(x):>17s} {fmt(n):>19s}")
    info = [r for r in rows if r["level"] == "info"]
    if info:
        out.append("  info rows: " + "; ".join(f"{r['id']} {r['op']} {' '.join(r['reasons'])}"
                                               for r in info[:12])
                   + (f"; ... ({len(info)} in all)" if len(info) > 12 else ""))
    out.append("  excluded (C0 < 1e-4): " + (", ".join(f"{e['id']} (C0 = {e['C0']:.3g})"
                                                       for e in excluded) or "none"))
    return out


# ---------------------------------------------------------------------------

def core_table(items, rows):
    lines = ["| id | xyz-astm-d65 | helmgen (L, a, b) | helmlab-metric (L, a, b) |",
             "|---|---|---|---|"]
    notes = {}
    for inp in items:
        g = next(r for r in rows if r["id"] == inp["id"] and r["op"] == "helmgen.fwd")
        m = next((r for r in rows if r["id"] == inp["id"] and r["op"] == "helmlab-metric.fwd"), None)
        cells = [inp["id"], ", ".join(repr(v) for v in inp["xyz"]), ", ".join(repr(v) for v in g["out"]),
                 ", ".join(repr(v) for v in m["out"]) if m else f"no row (C0 = {inp['c0']!r})"]
        lines.append("| " + " | ".join(cells) + " |")
    for r in rows:
        if r["level"] != "U":
            notes.setdefault((r["op"], r["level"]), []).append(r["id"])
    lines.append("")
    lines.append("Rows not at level U (domain out):")
    for (op, level), ids in notes.items():
        lines.append(f"- {op}, level {level}: " + ", ".join(ids))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate the Helmlab 2026 test vectors.")
    parser.add_argument("--out", metavar="DIR", help="write the vector files and their SHA256SUMS "
                        "into DIR instead of the repository")
    parser.add_argument("--core-table", action="store_true",
                        help="print core-1 as a Markdown table and write nothing")
    args = parser.parse_args()

    src = check_source()
    gen_params = json.loads((SPEC / "params" / "gen_params.json").read_text())
    metric_params = json.loads((SPEC / "params" / "metric_params.json").read_text())
    ref = Reference()
    builder = Builder(ref, gen_params, metric_params)

    with np.errstate(all="ignore"), warnings.catch_warnings():
        # Out of the domain the reference warns (overflow, invalid power, the XYZ check of
        # MetricSpace.distance); the values are recorded as they are.
        warnings.simplefilter("ignore")

        neg_rows, neg_draws = builder.neg_lms()
        core_items, core_groups = builder.core(neg_rows[0])
        core_rows, core_excl, black_core = build_rows(core_items, OPS_SIX, ref)
        core_fail = assign_levels(core_rows, ref)
        if args.core_table:
            print(core_table(core_items, core_rows))
            for r in core_fail:
                print(f"FAIL (domain in): {r['id']} {r['op']}: {' '.join(r['reasons'])}")
            return 1 if core_fail else 0

        conv_items, conv_groups = builder.conv()
        grey_items, grey_groups = builder.grey()
        edge_items, edge_groups = builder.edge(neg_rows, neg_draws)
        pairs, de_groups = builder.de()

        conv_rows, conv_excl, black_conv = build_rows(conv_items, OPS_SIX, ref)
        grey_rows, grey_excl, black_grey = build_rows(grey_items, OPS_FIVE, ref)
        edge_rows, edge_excl, black_edge = build_rows(edge_items, OPS_FIVE, ref)
        de_rows, de_same = build_de_rows(pairs, ref, builder)
        neutral_items = grey_items + conv_items
        neutral_rows, neutral_excl, black_neutral = build_rows(neutral_items, OPS_NEUTRAL, ref)

        sets = [  # (set, file, rows, excluded, groups, ops)
            ("core", "core-1.json", core_rows, core_excl, core_groups, OPS_SIX),
            ("conv", "conv-1.json", conv_rows, conv_excl, conv_groups, OPS_SIX),
            ("grey", "grey-1.json", grey_rows, grey_excl, grey_groups, OPS_FIVE),
            ("edge", "edge-1.json", edge_rows, edge_excl, edge_groups, OPS_FIVE),
            ("de", "de-1.json", de_rows, [], de_groups, ("deltaE",)),
            ("neutral", "neutral-1.json", neutral_rows, neutral_excl,
             {**grey_groups, **conv_groups}, OPS_NEUTRAL),
        ]
        failures = [("core-1.json", r) for r in core_fail]
        for name, file_name, rows, *_ in sets[1:]:
            failures += [(file_name, r) for r in assign_levels(rows, ref)]

        # Checks of the reference reported at the end (design Sections 2 and 10).
        pick = [neutral_items[i * len(neutral_items) // 20] for i in range(20)]
        neutral_same = sum(
            1 for inp in pick
            if np.array(ref.one("helmlab-metric-neutral.fwd", inp["xyz"])).tobytes()
            == np.asarray(ref.hl.metric.from_xyz(np.array(inp["xyz"]).reshape(1, 3))[0]).tobytes())
        lch_same = sum(
            1 for r in core_rows if r["op"] == "helmgen.fwd"
            and np.asarray(ref.hl.gen.to_lch(np.array(r["out"]))).tobytes()
            == np.asarray(ref.hl.gen.to_lch(ref.hl.gen.lab(*r["out"]))).tobytes())
        n_lch = sum(1 for r in core_rows if r["op"] == "helmgen.fwd")

    common = {
        "params_sha256": src["params_sha256"],
        "reference": {
            "package": "helmlab (Python), loaded from src/helmlab of the repository",
            "version": src["version"],
            "source_tree": src["tree"],
            "note": "source_tree is the git tree id of src/helmlab. Its files equal tag v1.0.0 "
                    f"(commit {src['tag'][:7]}) except the version string in __init__.py "
                    f"(\"{src['version']}\")"
                    + (f" and docstrings in {', '.join(src['docstring_only'])} (same syntax tree "
                       "once docstrings are removed)" if src['docstring_only'] else "") + ".",
            "calling_convention": "One reference call with an input of shape (1, 3) per input "
                                  "(spec Appendix A.2); GenSpace(), MetricSpace(), "
                                  "MetricSpace(neutral_correction=True) and Helmlab() are built "
                                  "once and reused.",
        },
        "environment": environment(),
        "generator": {"path": GENERATOR, "sha256": sha256_file(Path(__file__))},
        "frames": frames_header(),
    }

    texts = {}
    for set_name, file_name, rows, excluded, groups, ops in sets:
        header = {
            "name": "Helmlab 2026 test vectors",
            "set": set_name,
            "file": file_name,
            "spec": SPEC_REV,
            "license": "CC0-1.0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
            "params_sha256": common["params_sha256"],
            "reference": dict(common["reference"], calls={op: CALLS[op] for op in ops}),
            "environment": common["environment"],
            "generator": common["generator"],
            "variant": VARIANT,
            "frames": common["frames"],
            "tolerances": TOLERANCES,
            "domain_rules": DOMAIN_RULES,
            "level_rules": LEVEL_RULES,
            "groups": groups,
            "excluded": excluded,
            "counts": counts_of(rows),
        }
        texts[file_name] = file_text(header, rows)
    api = API_SOURCE.read_bytes()

    # Report -------------------------------------------------------------------
    report = ["Helmlab 2026 test vectors: generator report", ""]
    report.append(f"reference: {ref.file} (version string in __init__.py: {src['version']}; "
                  f"helmlab.__version__ reports {ref.version_attr})")
    report.append(f"HEAD {src['head']}; src/helmlab tree {src['tree']} equals v1.0.0 "
                  f"({src['tag'][:7]}) apart from the version string in __init__.py"
                  + (f" and docstrings in {', '.join(src['docstring_only'])}" if src['docstring_only'] else "")
                  + "; no uncommitted change under src/helmlab")
    report.append("params sha256: " + ", ".join(f"{k} {v}" for k, v in src["params_sha256"].items()))
    report.append("frames: A*W_astm = W_css and B*W_css = W_astm in exact arithmetic; "
                  "CSS Color 4 lin_P3_to_XYZ and lin_2020_to_XYZ map (1, 1, 1) to W_css exactly")
    report.append("")
    for set_name, file_name, rows, excluded, groups, ops in sets:
        report += file_report(file_name, texts[file_name], rows, excluded)
        report.append("")
    report.append(f"api-1.0-1.json: byte copy of {API_SOURCE.relative_to(REPO)}, {len(api)} bytes, "
                  f"sha256 {sha256_bytes(api)}")
    report.append("")
    report.append("Over domain-in rows of all files, per op: max |one-row - batch|, "
                  "max round trip (c), max |xyz rt|, max change (d)")
    all_rows = [r for s in sets for r in s[2]]
    for op in OPS:
        inside = [r for r in all_rows if r["op"] == op and r["domain"] == "in"]
        if not inside:
            continue
        b = max(r["d_batch"] for r in inside)
        c = max((r["d_trip"] for r in inside if "d_trip" in r), default=None)
        x = max((r["d_xyz"] for r in inside if "d_xyz" in r), default=None)
        n = max((r["d_nudge"] for r in inside if "d_nudge" in r), default=None)
        report.append(f"  {op:28s} rows {len(inside):6d}   {fmt(b):>9s} {fmt(c):>9s} {fmt(x):>9s} "
                      f"{fmt(n):>9s}")
    report.append("")
    report.append(f"deltaE: MetricSpace().distance(xyz1, xyz2) equals distance_from_lab(lab1, lab2) "
                  f"bit for bit on {de_same} of {len(de_rows)} de-1 rows")
    report.append(f"neutral: MetricSpace(neutral_correction=True).from_XYZ equals "
                  f"Helmlab().metric.from_xyz bit for bit on {neutral_same} of {len(pick)} inputs "
                  f"({', '.join(p['id'] for p in pick)})")
    report.append(f"helmgenlch: Helmlab().gen.to_lch accepts a plain ndarray; equal to the call on "
                  f"hl.gen.lab(*lab) bit for bit on {lch_same} of {n_lch} core-1 inputs")
    black_ok = black_core and black_conv and black_grey and black_edge and black_neutral
    report.append(f"black (0, 0, 0): MetricSpace outputs exactly (0, 0, 0) in every file: "
                  f"{'yes' if black_ok else 'NO'}")
    report.append(f"neg-lms: {len(neg_rows)} of {neg_draws} rows have a negative component of "
                  f"M1*XYZ (spec Section 5.3.5 gives 566); edge-1 uses the first 256")
    print("\n".join(report))

    if failures or not black_ok:
        print("")
        for file_name, r in failures:
            print(f"FAIL (domain in) {file_name} {r['id']} {r['op']}: {' '.join(r['reasons'])}")
        if not black_ok:
            print("FAIL black input without output (0, 0, 0)")
        print("Nothing written.")
        return 1

    # Write ----------------------------------------------------------------------
    out_dir = Path(args.out).resolve() if args.out else VECTOR_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for file_name, text in texts.items():
        data = text.encode("ascii")
        (out_dir / file_name).write_bytes(data)
        written[file_name] = sha256_bytes(data)
    (out_dir / "api-1.0-1.json").write_bytes(api)
    written["api-1.0-1.json"] = sha256_file(out_dir / "api-1.0-1.json")
    if written["api-1.0-1.json"] != sha256_file(API_SOURCE):
        stop("api-1.0-1.json differs from its source")
    names = sorted(written)
    if args.out:
        sums = "".join(f"{written[n]}  {n}\n" for n in names)
        (out_dir / "SHA256SUMS").write_text(sums)
        print(f"\nwrote {len(names)} files and SHA256SUMS to {out_dir}")
    else:
        sums = "".join(line + "\n" for line in src["params_lines"])
        sums += "".join(f"{written[n]}  vectors/{n}\n" for n in names)
        SUMS.write_text(sums)
        print(f"\nwrote {len(names)} files to {out_dir} and updated {SUMS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
