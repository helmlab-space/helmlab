"""Print the parameter tables of spec/2026/spec.md from the parameter files.

Usage: python3 tools/spec/param_tables.py spec/2026/params/gen_params.json spec/2026/params/metric_params.json

Each number is printed twice: the shortest decimal that round-trips to the same
IEEE 754 binary64 value (Python repr), and the exact value as a hexadecimal float.
"""
import hashlib
import json
import sys

SPACE = {"gen_params.json": "GenSpace", "metric_params.json": "MetricSpace"}


def rows(prefix, value):
    if isinstance(value, dict):
        for k, v in value.items():
            yield from rows(f"{prefix}.{k}" if prefix else k, v)
    elif isinstance(value, list):
        for i, v in enumerate(value):
            yield from rows(f"{prefix}[{i}]", v)
    else:
        yield prefix, value


def table(path):
    raw = open(path, "rb").read()
    params = json.loads(raw)
    sha = hashlib.sha256(raw).hexdigest()
    fname = path.split('/')[-1]
    out = [f"{SPACE[fname]}, `params/{fname}`, sha256 `{sha}`", "", "| Parameter | Decimal | Hexadecimal |", "|---|---|---|"]
    for name, v in rows("", params):
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            out.append(f"| `{name}` | `{json.dumps(v)}` | |")
        else:
            f = float(v)
            out.append(f"| `{name}` | `{f!r}` | `{f.hex()}` |")
    return "\n".join(out)


if __name__ == "__main__":
    print("\n\n".join(table(p) for p in sys.argv[1:]))
