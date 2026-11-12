"""Check the Python reference package against the Helmlab 2026 test vectors.

Usage, from the repository root:

    python3 -I tools/spec/check_vectors.py [--vectors DIR]

DIR defaults to spec/2026/vectors. Every row with level U or R is recomputed
with one-row calls of the package in src/helmlab (input of shape (1, 3), spec
Appendix A.2) and compared with the level R tolerances of spec Section 8,
absolute, per channel: forward and deltaE 1e-11, inverse 1e-8. helmgenlch
outputs are compared as L, C and the hue chord C*|dh|, dh in radians wrapped
to [-pi, pi], C the expected chroma. Inverse outputs are compared in the
xyz-astm-d65 frame. Rows with level info are counted and skipped. A JSON file
without "rows" (api-1.0-1.json, the regression set of the package) is
skipped.

Prints, per file and op, the rows checked, the info rows skipped, the largest
difference and the failures. Exit status 1 on any failure.
"""
import argparse
import json
import math
import sys
import warnings
from pathlib import Path

sys.dont_write_bytecode = True  # leave src/helmlab/__pycache__ as it is

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
TOL = {"fwd": 1e-11, "inv": 1e-8, "deltaE": 1e-11}  # level R, spec Section 8
SHOW = 5  # failing rows printed per file and op


def load_reference():
    sys.path.insert(0, str(SRC))
    import helmlab
    from helmlab.helmlab import Helmlab
    from helmlab.spaces.gen import GenSpace
    from helmlab.spaces.metric import MetricSpace
    path = Path(helmlab.__file__).resolve()
    if SRC.resolve() not in path.parents:
        print(f"check_vectors.py: helmlab was loaded from {path}, not from {SRC}", file=sys.stderr)
        sys.exit(1)
    return path, GenSpace(), MetricSpace(), MetricSpace(neutral_correction=True), Helmlab()


def number(v):
    """A JSON number, or one of the strings "NaN", "Infinity", "-Infinity" (float() reads both)."""
    return float(v)


def triple(values):
    return np.array([number(v) for v in values], dtype=np.float64).reshape(1, 3)


def max_abs(got, exp):
    worst = 0.0
    for g, e in zip(got, exp):
        d = abs(g - e)
        if not d <= worst:  # also catches NaN
            worst = d if not math.isnan(d) else math.inf
    return worst


def chord_diff(got, exp):
    """L, C and the hue chord C*|dh|; C is the expected chroma.

    An undefined hue (NaN) counts as dh = pi (spec.md, Section 8).
    """
    if not (math.isfinite(got[0]) and math.isfinite(got[1])) or math.isinf(got[2]):
        return math.inf  # max() below would drop a NaN that is not its first argument
    if math.isnan(got[2]):
        dh = math.pi
    else:
        dh = math.remainder(math.radians(got[2] - exp[2]), 2 * math.pi)
    d = max(abs(got[0] - exp[0]), abs(got[1] - exp[1]), exp[1] * abs(dh))
    return d if not math.isnan(d) else math.inf


def main():
    parser = argparse.ArgumentParser(description="Check helmlab against the Helmlab 2026 vectors.")
    parser.add_argument("--vectors", metavar="DIR", default=str(REPO / "spec" / "2026" / "vectors"))
    args = parser.parse_args()
    vector_dir = Path(args.vectors)

    path, gen, metric, neutral, hl = load_reference()

    def recompute(op, row):
        """Returns (got, expected, kind) for one row, with one-row reference calls."""
        if op == "deltaE":
            lab1, lab2 = triple(row["in"]["lab1"]), triple(row["in"]["lab2"])
            return [float(metric.distance_from_lab(lab1, lab2)[0])], [number(row["out"])], "deltaE"
        kind = op.rsplit(".", 1)[1]
        if kind == "fwd":
            x = triple(row["in"]["xyz-astm-d65"])
            exp = [number(v) for v in row["out"]]
            if op == "helmgen.fwd":
                got = gen.from_XYZ(x)[0]
            elif op == "helmgenlch.fwd":
                got = hl.gen.to_lch(gen.from_XYZ(x)[0])
            elif op == "helmlab-metric.fwd":
                got = metric.from_XYZ(x)[0]
            elif op == "helmlab-metric-neutral.fwd":
                got = neutral.from_XYZ(x)[0]
            else:
                raise KeyError(op)
        elif kind == "inv":
            v = triple(row["in"])
            exp = [number(t) for t in row["out"]["xyz-astm-d65"]]
            if op == "helmgen.inv":
                got = gen.to_XYZ(v)[0]
            elif op == "helmgenlch.inv":
                got = gen.to_XYZ(np.asarray(hl.gen.from_lch(v[0])).reshape(1, 3))[0]
            elif op == "helmlab-metric.inv":
                got = metric.to_XYZ(v)[0]
            elif op == "helmlab-metric-neutral.inv":
                got = neutral.to_XYZ(v)[0]
            else:
                raise KeyError(op)
        else:
            raise KeyError(op)
        return [float(t) for t in got], exp, kind

    print(f"reference: {path}")
    print(f"vectors: {vector_dir.resolve()}")
    print(f"{'file':16s} {'op':28s} {'checked':>8s} {'info':>6s} {'largest diff':>13s} {'tol':>7s} "
          f"{'failures':>8s}")
    total = {"checked": 0, "info": 0, "failures": 0}
    files = sorted(vector_dir.glob("*.json"))
    if not files:
        print(f"check_vectors.py: no .json files in {vector_dir}", file=sys.stderr)
        return 1
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the reference warns outside the domain
        for file in files:
            data = json.loads(file.read_text())
            if "rows" not in data:
                print(f"{file.name:16s} skipped: no rows (not in the vector format)")
                continue
            stats = {}  # op -> [checked, info, largest diff, failures, examples]
            for row in data["rows"]:
                op, level = row["op"], row["level"]
                s = stats.setdefault(op, [0, 0, 0.0, 0, []])
                if level == "info":
                    s[1] += 1
                    continue
                problem = None
                if level not in ("U", "R") or (level == "U") != (row["domain"] == "in"):
                    problem = f"level {level!r} with domain {row['domain']!r}"
                    d = math.inf
                    tol = 0.0
                else:
                    try:
                        got, exp, kind = recompute(op, row)
                    except (KeyError, ValueError) as e:  # unknown op, or a value check of the package
                        got, exp, kind = None, None, None
                        problem = f"{type(e).__name__}: {e}"
                    tol = TOL[kind] if kind else 0.0
                    if got is None:
                        d = math.inf
                    elif not all(math.isfinite(v) for v in exp):
                        d = math.inf
                        problem = "expected output not finite in a U or R row"
                    elif op == "helmgenlch.fwd":
                        d = chord_diff(got, exp)
                    else:
                        d = max_abs(got, exp)
                s[0] += 1
                s[2] = max(s[2], d)
                if problem is not None or not d <= tol:
                    s[3] += 1
                    if len(s[4]) < SHOW:
                        s[4].append(f"{row['id']}: " + (problem or f"difference {d!r} > {tol!r}"))
            for op, (checked, info, largest, failures, examples) in stats.items():
                kind = "deltaE" if op == "deltaE" else op.rsplit(".", 1)[1]
                print(f"{file.name:16s} {op:28s} {checked:8d} {info:6d} {largest:13.3g} "
                      f"{TOL.get(kind, 0.0):7.0e} {failures:8d}")
                for line in examples:
                    print(f"    FAIL {line}")
                total["checked"] += checked
                total["info"] += info
                total["failures"] += failures
    print(f"total: {total['checked']} rows checked, {total['info']} info rows skipped, "
          f"{total['failures']} failures")
    if total["checked"] + total["info"] == 0:
        print(f"check_vectors.py: no file in {vector_dir} has rows", file=sys.stderr)
        return 1
    return 1 if total["failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
