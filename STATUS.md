# Status

Last updated: 2026-11-16.

## Releases

| Package | Version | Notes |
|---|---|---|
| `helmlab` (npm, PyPI) | 1.0.2 | Same outputs as 1.0.0 (Helmlab 2026); documentation release |
| `postcss-helmlab` (npm) | 0.3.1 | README only |

`.github/workflows/deploy.yml` publishes to npm and PyPI (trusted publishing) when a package version changes. `.github/workflows/site.yml` deploys `site/` to helmlab.space.

## Definition

Helmlab 2026 is the output of `helmlab` 1.0.0. Its specification is Spec Rev 1 of 2026-11-16: `spec/2026/spec.md`, published at https://helmlab.space/spec/2026/ and tagged `spec-2026-rev1`. Test vectors are in `spec/2026/vectors/`; `tools/spec/check_vectors.py` and `tools/spec/check_vectors.mjs` check the packages against them, and `tools/spec/build_site.py` generates the site page from `spec.md`.

Parameter files:

- GenSpace: `src/helmlab/data/gen_params.json` and `packages/helmlab-js/src/data/gen-params.json`.
- MetricSpace: `src/helmlab/data/metric_params.json` and `packages/helmlab-js/src/data/params.json`. The JavaScript neutral-correction table is `packages/helmlab-js/src/data/neutral-lut.json`.

`checkpoints/` holds parameter files from development, among them `genspace_v0.11.1.json` and `metricspace_v21.json`.

## In preparation

- The rest of the documentation site at helmlab.space, with one evidence page. Performance figures will be published only there.

## Rules for texts in this repository

- No performance figures in READMEs, docstrings, type files or skills; link to the evidence page instead.
- Figures from ColorBench versions before v4 (for example GenSpace vs OKLab win counts) and the figures of the previous helmlab.space are withdrawn. Do not reuse them.
- MetricSpace was fitted to COMBVD, MacAdam 1974 and Helmlab's own Human Feedback data (weights 0.55, 0.15 and 0.30); He et al. 2022 was used in model selection. Results on these datasets are in-sample.
- Known limitations are listed in the README.

## Tests

```bash
PYTHONPATH=src python3 -m pytest tests/ -q
cd packages/helmlab-js && npx vitest run
cd packages/postcss-helmlab && node test.js
```

`tools/golden/` compares package outputs across versions.

## Other implementations

Color.js 0.7.x and ColorAide implement `helmlab-metric` in its base form. color-space 3.5.0 and 3.6.0 ship the neutral-corrected form under the name `helmlab-metric`. Conformance: Section 8 of the specification.
