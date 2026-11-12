# Checks of other implementations

These scripts produced the results under "Other implementations" in Section 8 of `spec/2026/spec.md`. Each one runs a third-party implementation on the rows at level U and R of every vector file and prints one JSON object: for each op and level, the number of rows, the largest difference, the row where it occurs and, for level U, the number of rows above the level U tolerance.

| Script | Implementations |
|---|---|
| `run_js.mjs` | Color.js (`colorjs.io`) and color-space |
| `run_py.py` | ColorAide |

The last argument selects the XYZ frame of the inputs and expected outputs: `css` for `xyz-css-d65`, `astm` for `xyz-astm-d65` (Section 9.3).

- color-space ships one MetricSpace form, under the name `helmlab-metric`. `run_js.mjs` checks it against both the base-form rows and the rows of `neutral-1.json`.
- An undefined hue in `helmgenlch` counts as a chord of C·π (Section 8).
- ΔE rows are computed from the XYZ of both colours (`xyz1`, `xyz2`). color-space has no ΔE `helmlab`.

The runs of 2026-11-12 used Node.js 26.8.1 and Python 3.14.7, with each package in a directory of its own. From this directory:

```sh
npm install --prefix colorjs-0.7.1 colorjs.io@0.7.1
node run_js.mjs colorjs colorjs-0.7.1 ../../../spec/2026/vectors css

npm install --prefix color-space-3.6.0 color-space@3.6.0
node run_js.mjs colorspace color-space-3.6.0 ../../../spec/2026/vectors css

python3 -m pip install --target coloraide-8.14.1 coloraide==8.14.1
python3 -I run_py.py coloraide-8.14.1 ../../../spec/2026/vectors css
```

The same commands with colorjs.io 0.7.0, color-space 3.5.0 and ColorAide 8.9 give the other rows of the table in Section 8.
