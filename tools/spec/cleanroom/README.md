# Clean-room implementations of Section 5

On 2026-11-09 Section 5 of `spec/2026/spec.md` was checked by implementing it from the text alone. Each implementation was written without access to the reference code. Its test script compares it with the reference code, the Python package in `src/`, which it uses only through its public calls. The results are in Appendix A.4 of the spec.

| Folder | Implements | Written from |
|---|---|---|
| `genspace/` | GenSpace forward and inverse, LCh, sRGB input | Draft Section 5.1 |
| `metricspace/` | MetricSpace base form, ΔE, neutral-corrected form, sRGB input | Draft Sections 5.3 to 5.6 |

The implementations keep the numbering of the draft they were written from. The formulas are the same in the final text, but the numbering has changed:

| Draft | Final |
|---|---|
| 5.1.1 and 5.3.1 Notation | 5.1 |
| 5.1.4 and 5.6 sRGB and hex input | 5.2 |
| 5.1 GenSpace; steps F1 to F7, I1 to I7 | 5.3; steps GF1 to GF7, GI1 to GI7 |
| 5.1.6 LCh | 5.4 |
| 5.3 MetricSpace | 5.5 |
| 5.4 ΔE | 5.6 |
| 5.5 Neutral correction | Appendix B |

## Running

From the repository root (needs numpy; the MetricSpace test also needs scipy):

```sh
python3 tools/spec/cleanroom/genspace/test.py          # M · v as left-to-right sums, as written in Section 5.1
python3 tools/spec/cleanroom/genspace/test.py blas     # M · v as numpy computes it for the whole batch
python3 tools/spec/cleanroom/metricspace/test.py
```

Each script prints the largest absolute difference per test set and PASS or FAIL against the tolerances of Section 8. With `blas`, the GenSpace implementation matched the reference bit for bit in every row on 2026-11-09 (macOS arm64, Python 3.14.7, numpy 2.4.4 with Accelerate, scipy 1.17.1).
