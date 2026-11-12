# Helmlab 2026

## 1. Status

| | |
|---|---|
| Name | Helmlab 2026 |
| Identifiers | `helmgen`, `helmgenlch`, `helmlab-metric`, ΔE `helmlab`; reserved: `helmlab-metric-neutral` (Appendix B) |
| Status | Draft |
| Revision | Spec Rev 1 |
| Date | Set on publication |
| DOI | Assigned on publication |
| Parameter files | `params/gen_params.json` sha256 `310ad1105eeec0aeb72bb8bf24185cda15f62ea8f61822dceca9d61b8a69356f`; `params/metric_params.json` sha256 `a469b93c278e1708ff32aee3308375ea0fccedac75bdeba4d3b85441e279401b` |
| Reference code | helmlab 1.0.0 on npm and PyPI (2026-07-08), commit `96a0f16`, tag `v1.0.0` |
| Editor | Görkem Yıldız |

This is a single-author definition. It is not a CIE, ISO or W3C standard.

## 2. Scope

Helmlab 2026 defines two colour spaces and one colour difference:

| Identifier | What it names |
|---|---|
| `helmgen` | GenSpace in Lab coordinates, for generating colours |
| `helmgenlch` | GenSpace in LCh coordinates |
| `helmlab-metric` | MetricSpace in its base form, for measuring colour differences |
| ΔE `helmlab` | The colour difference computed from `helmlab-metric` coordinates |

By default `hl.metric` in the helmlab package (Python and JavaScript) returns the neutral-corrected form `helmlab-metric-neutral` (Appendix B), not `helmlab-metric`. ΔE `helmlab` is always computed from `helmlab-metric` coordinates.

The definition is the output of helmlab 1.0.0. The releases 1.0.1 and 1.0.2 give the same outputs, and so do later releases; a change within the level R tolerance of Section 8 is made only as class K3 (Section 11). Releases before 1.0.0 are out of scope.

Normative sources, in order of precedence:

1. The parameter files (Section 6) and the published test vectors (Section 9).
2. The reference code of helmlab 1.0.0.

If the reference code and a vector differ within the tolerance of Section 8, the vector governs. A larger difference is an erratum (Section 11).

The description of the mathematics in Section 5 is informative. Earlier descriptions, including the paper (arXiv:2602.23010), differ from the reference code in places (E-2026-0001).

Functions that the helmlab package builds on these spaces (gamut mapping, gradients, palettes, tokens, contrast helpers) are not part of the definition. Their outputs may change in minor releases; `CHANGELOG.md` lists such changes under the heading "Changed outputs".

## 3. Intended use

*This section is dedicated to the public domain under CC0 1.0, so that implementations can copy it.*

- GenSpace is for generating colours: gradients, palettes and gamut mapping.
- MetricSpace coordinates are a distance embedding. Use them to compute ΔE `helmlab`. They are not recommended for LCh or hue operations: in the base form a neutral grey has non-zero a and b (a ≈ 0.24, b ≈ 0.15 at Y = 0.18), and L reaches about 2.16 for Rec.2020 colours (Section 9.1, `core/grey-0.18`; `conv-1`, row `rec2020/8-0-8`).
- MetricSpace is for SDR content: colours inside the Display P3 or Rec.2020 gamut (Section 7). Above diffuse white, use GenSpace; its grey axis stays neutral to about 3e-15 up to at least Y = 10 (Section 5.3.3).

Links from these statements to the evidence page will be added when that page is published; adding them is an editorial change (class K1, Section 11).

## 4. Conventions

- The key words MUST, MUST NOT, SHOULD, SHOULD NOT and MAY are to be interpreted as described in BCP 14 ([RFC 2119](https://www.rfc-editor.org/rfc/rfc2119), [RFC 8174](https://www.rfc-editor.org/rfc/rfc8174)) when, and only when, they appear in all capitals.
- Arithmetic is IEEE 754 binary64.
- Input is CIE 1931 XYZ relative to D65, scaled so that diffuse white has Y = 1.
- The reference frame `xyz-astm-d65` is the one the reference code uses: D65 white (0.95047, 1, 1.08883). Many other libraries use `xyz-css-d65`, the D65 white and Bradford adaptation of CSS Color 4. The test vectors give every XYZ in both frames (Section 9.3).
- L, a and b are the coordinates of a space. C = √(a² + b²) and h = atan2(b, a) in degrees, in [0, 360]; the reference code can return exactly 360 (Section 5.4).

## 5. Definition

*This section is informative (Section 2).*

It describes the reference code of helmlab 1.0.0 step by step. Parameters appear under their names in the JSON files (Section 6). Where the Python and JavaScript packages compute differently, the text follows the Python package, which produces the test vectors (Appendix A.3). The draft of 2026-11-09 of this text was checked by implementing it without access to the reference code and comparing the results with the reference code (Appendix A.4).

### 5.1 Notation

- Arithmetic is IEEE 754 binary64. Each formula is evaluated in the order written: `**` first, then unary minus, then `*` and `/`, then `+` and `-`, left to right, inside parentheses first. The order follows the reference code. A different order can change the last bits; near the zero-chroma curve of MetricSpace (Section 5.5.2) it can change more. The tolerances of Section 8 allow for the last bits.
- `M · v` is the product of a 3×3 matrix and a column vector: `(M · v)[i] = M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2]`. A JSON matrix lists its rows: `M[i][j]` is row i, column j. The Python reference computes the product with numpy (BLAS), which does not fix the order of the sums or the use of fused multiply-add; the last bit can depend on the number of rows in the call (Appendix A.2).
- `x**2` is `x*x`, as numpy computes it on arrays. Every other power `x**y` is the C library function `pow(x, y)`. In particular `x**3` is `pow(x, 3)`, which can differ from `x*x*x` in the last bit. (A call with a single colour of shape (3,) evaluates some squares as `pow(x, 2)`; Appendix A.2.)
- `fmod(x, y)` is the remainder of x / y with the quotient truncated toward zero. It has the sign of x and is exact (C `fmod`, JavaScript `%`).
- `clip(x, lo, hi) = min(max(x, lo), hi)`. `min`, `max` and `clip` return NaN when an argument is NaN, as numpy does (IEEE 754-2019 `minimum` and `maximum`).
- `sign(x)` is −1, 0 or +1 for x < 0, x = 0 and x > 0.
- `pi` is the binary64 value `3.141592653589793`.
- `atan2(y, x)` is the C function, with results in [−pi, pi]. For signed zeros, `atan2(+0, +0) = +0`, `atan2(+0, −0) = +pi`, `atan2(−0, +0) = −0` and `atan2(−0, −0) = −pi`.
- `sinh` and `asinh` are the hyperbolic sine and its inverse. `hypot(x, y)` is `sqrt(x**2 + y**2)` computed as by the C function `hypot` (Python `math.hypot` uses another method and can differ in the last bit).
- Hue angles inside the pipelines are in radians and are not wrapped unless a step says so.
- Every comparison with NaN is false. NaN in the input gives NaN in the output, for each colour separately.
- Step labels: GenSpace uses GF1 to GF7 (forward) and GI1 to GI7 (inverse); MetricSpace uses F1 to F14 and I1 to I14; ΔE uses D1 to D7; the neutral correction of Appendix B uses N1 to N8.

### 5.2 sRGB and hex input

The spaces take XYZ. For hex strings and sRGB triples the reference code first converts to XYZ as below (Python `helmlab.utils.srgb_convert.sRGB_to_XYZ`, JavaScript `srgbToXyz`). `hl.gen.from_hex`, `hl.gen.from_srgb` and the colour-string parser of `hl.metric` use it. This conversion is not part of the definition. It is given here because test vectors that start from sRGB use it.

1. Hex: `#rrggbb` gives (R, G, B) = (rr/255, gg/255, bb/255), each pair two hexadecimal digits. `#rgb` is first expanded to `#rrggbb`. JavaScript checks the digits with a pattern. Python also accepts a sign or a space within a pair and several leading `#`; for example, `#-1-1-1` gives −1/255 in each channel.
2. Each channel v is linearized. There is no clip: a negative v uses the first branch and v > 1 the second.

```text
lin(v) = v / 12.92                        if v <= 0.04045
lin(v) = ((v + 0.055) / 1.055) ** 2.4     otherwise
```

3. The linear values are multiplied by a fixed matrix.

```text
XYZ = M_sRGB · (lin(R), lin(G), lin(B))

M_sRGB = [[0.4124564390896922, 0.3575760776439511, 0.1804374832663989],
          [0.2126728514056226, 0.7151521552878178, 0.0721750036064596],
          [0.0193338955823293, 0.1191920258813418, 0.9503040785363679]]
```

The rows of `M_sRGB` do not sum exactly to the reference white of Section 4: sRGB white (1, 1, 1) gives XYZ = (0.9504700000000421, 1.0000000102999, 1.088830000000039). The strings `color(display-p3 r g b)` and `color(rec2020 r g b)` use other matrices and transfer functions, which are not described here.

### 5.3 GenSpace (`helmgen`)

The steps follow `GenSpace.from_XYZ` and `GenSpace.to_XYZ` in `src/helmlab/spaces/gen.py`, with the parameters of `params/gen_params.json`.

`helmgen` is the GenSpace class of the reference code built with its default arguments: the parameters of `gen_params.json` and `neutral_correction = False`. `hl.gen` uses this object in both packages.

Each colour is computed on its own: for input of shape (3,) or (N, 3), no step depends on other colours in the same call. (An input of shape (N, 1, 3) gives an output of shape (N, N, 3); Appendix A.2.)

#### 5.3.1 Parameters

| Entry of `gen_params.json` | Used in | Role |
|---|---|---|
| `M1` (3×3) | GF1, GI7 | Matrix from XYZ to the responses LMS. The JSON array lists the rows: `M1[i][j]` is row i, column j. GI7 uses its inverse. |
| `M2` (3×3) | GF4, GI4 | Matrix from the compressed responses to (L, a, b). GI4 uses its inverse. |
| `transfer` = `"depcubic"` | GF2, GF3, GI5, GI6 | Selects the depressed-cubic transfer (GF2, GI6) and turns on the neutral blend (GF3, GI5). |
| `depcubic_alpha` | GF2, GI6 | The constant alpha of the cubic `y**3 + alpha*y = x`. |
| `L_corr_pw` (19 values) | GF5, GI3 | Shifts of the 19 inner knots of the piecewise-linear lightness correction. A non-empty list turns GF5 and GI3 on. |
| `L_corr_pw_step` | GF5, GI3 | Spacing of the knots on the input side. |
| `enrichment.type` = `"L_gated_hue"` | GF6, GI2 | Turns the lightness-gated hue rotation on, together with `abs(enrichment.amp) > 1e-10`. |
| `enrichment.amp` | GF6, GI2 | Largest rotation, in radians. |
| `enrichment.center_deg` | GF6, GI2 | Hue at the centre of the rotation window, in degrees. |
| `enrichment.sigma` | GF6, GI2 | Width of the Gaussian hue window, in radians. |
| `enrichment.L_lo`, `enrichment.L_hi` | GF6, GI2 | Lightness range of the gate. |
| `chroma_power` | GF7, GI1 | Exponent of the chroma compression. A value other than 1 turns GF7 and GI1 on. |
| `enrichment_gr`, `enrichment_rt` | none | Not read by the reference code (Python `GenParams.from_dict`, JavaScript `normalizeGenParams`). They have no effect on any output. They are top-level keys, not members of the `enrichment` object. |

The reference code has further stages. Their parameters are absent from `gen_params.json`, and the default values turn them off (Section 5.3.6).

#### 5.3.2 Derived constants

These values are computed once from the parameters. The decimal values are those of the Python reference.

```text
alpha = depcubic_alpha                     = 0.021
s     = sqrt(alpha / 3)                    = 0.08366600265340755
k     = 2 * s**3                           = 0.0011713240371477058
h_c   = enrichment.center_deg * (pi / 180) = 4.616395871525001
        (pi / 180 = 0.017453292519943295)
```

Knot tables of the lightness correction, with `n = 19` (the length of `L_corr_pw`) and `n + 2 = 21` entries each:

```text
P_in[i]    = i * L_corr_pw_step             for i = 0 .. n
P_in[n+1]  = 1.0
P_out[0]   = 0.0
P_out[i]   = P_in[i] + L_corr_pw[i-1]       for i = 1 .. n
P_out[n+1] = 1.0
```

| i | P_in[i] | P_out[i] |
|---|---|---|
| 0 | 0.0 | 0.0 |
| 1 | 0.05 | 0.009494013522189627 |
| 2 | 0.1 | 0.02564569838030986 |
| 3 | 0.15000000000000002 | 0.055259661658689105 |
| 4 | 0.2 | 0.10574901531227408 |
| 5 | 0.25 | 0.16055853320726027 |
| 6 | 0.30000000000000004 | 0.21405964892993756 |
| 7 | 0.35000000000000003 | 0.26786230508811226 |
| 8 | 0.4 | 0.3220435246104499 |
| 9 | 0.45 | 0.3739052098520243 |
| 10 | 0.5 | 0.43020997780918835 |
| 11 | 0.55 | 0.4835465162128873 |
| 12 | 0.6000000000000001 | 0.5399824670411353 |
| 13 | 0.65 | 0.5956710081330342 |
| 14 | 0.7000000000000001 | 0.6542161666450478 |
| 15 | 0.75 | 0.7115380216519989 |
| 16 | 0.8 | 0.7702762412711669 |
| 17 | 0.8500000000000001 | 0.8293313467712837 |
| 18 | 0.9 | 0.889406386197059 |
| 19 | 0.9500000000000001 | 0.9462829573474728 |
| 20 | 1.0 | 1.0 |

Both tables are strictly increasing. The slopes of the 20 segments lie between 0.18988 (first segment) and 1.2015.

`M1_inv` and `M2_inv` are the inverses of `M1` and `M2`. The Python reference computes them with `numpy.linalg.inv` (LU decomposition); another method can differ in the last bits. The values of the Python reference are:

```text
M1_inv = [[ 1.2326502723725545,    -0.5557414732179251,  0.2735612008453706  ],
          [-0.040766587876246804,   1.1122097325799587, -0.0714431447037118  ],
          [-0.07673215022426162,   -0.4216188546558442,  1.5871810048801054  ]]

M2_inv = [[ 0.9930001151336145,     0.32599327253052285,  0.12945085631713896 ],
          [ 0.9930001151336143,    -0.08708353111074632, -0.038613617430049534],
          [ 0.993000115133614,     -0.12386097008215027, -0.8351991365871065  ]]
```

#### 5.3.3 Forward: XYZ to Lab

Input: XYZ = (X, Y, Z). Output: (L, a, b). The steps run in the order GF1 to GF7:

```text
XYZ -> GF1 M1 and clip -> GF2 cubic transfer -> GF3 neutral blend -> GF4 M2
    -> GF5 lightness correction -> GF6 hue rotation -> GF7 chroma compression -> (L, a, b)
```

**GF1. LMS responses.** Multiplies XYZ by `M1` and sets negative components to 0.

```text
LMS    = M1 · XYZ
LMS[i] = max(LMS[i], 0)          for i = 0, 1, 2
```

Edge cases: the clip discards the negative part of each component, so the inverse cannot recover input with a negative component of `M1 · XYZ` (negative XYZ, imaginary colours). If no component is positive, the output is (0, 0, 0). A component equal to −inf is clipped to 0 like any other negative value. A component equal to +inf makes the output (NaN, NaN, NaN) (GF2, GF3).

**GF2. Transfer.** Replaces each component x = LMS[i] by the real root y of `y**3 + alpha*y = x`. The root is computed in closed form and refined by one Halley step.

```text
t   = x / k
y0  = 2 * s * sinh(asinh(t) / 3)
f   = y0**3 + alpha * y0 - x
fp  = 3 * y0**2 + alpha
fpp = 6 * y0
den = 2 * fp * fp - f * fpp
y   = y0 - 2 * f * fp / den      if abs(den) > 1e-30
y   = y0                         otherwise (also when den is NaN)
```

The three results are c = (c0, c1, c2).

Edge cases: x = 0 gives y = 0. There is no upper clip. For x above about 3e192 the products in the Halley step overflow, and GF2 returns +inf, −inf or NaN for some inputs; x = +inf gives y = +inf, because den is NaN. GF3 turns an infinite or NaN component into NaN for all three components, so the output is then (NaN, NaN, NaN), never ±inf (reference behaviour, no meaning guarantee).

**GF3. Neutral blend.** Moves the three components towards their mean when their relative spread is small.

```text
m  = (c0 + c1 + c2) / 3
r  = (max(c0, c1, c2) - min(c0, c1, c2)) / max(abs(m), 1e-30)
w  = exp(-((r / 1e-5)**2))
ci = ci + w * (m - ci)           for i = 0, 1, 2
```

Edge cases: for r below about 7.5e-14, w = 1 and the three components become equal to m: exactly when abs(m) >= 1e-30, and to rounding when the guard 1e-30 acts. For r above about 2.73e-4, w = 0 (exp underflows) and the components do not change. Between these values the spread is reduced but not removed. Grey input (Y times the reference white) has r near 1e-16, so its three components become equal. The guard 1e-30 acts whenever abs(m) < 1e-30: r is then smaller than the relative spread, so very dark colours are moved towards grey. For the sRGB primaries the blend starts below Y ≈ 1e-35, and below Y ≈ 1e-45 the output is grey to rounding. For black, m = 0, r = 0 and w = 1. An infinite or NaN component makes all three components NaN.

**GF4. Matrix M2.** Multiplies the blended responses by `M2`.

```text
(L, a, b) = M2 · (c0, c1, c2)
```

Edge cases: none. For grey input (c0 = c1 = c2), a and b are not exactly 0. They come from the sums of the rows M2[1] and M2[2], which are −2.220446049250313e-16 and −4.440892098500626e-16 in binary64, and from rounding in this product.

**GF5. Lightness correction.** Maps L through the piecewise-linear function with knots (P_in[i], P_out[i]) when 0 <= L <= 1, and leaves L unchanged otherwise.

```text
if 0 <= L <= 1:
    j = the largest index in 0 .. n with P_in[j] <= L
    t = clip((L - P_in[j]) / max(P_in[j+1] - P_in[j], 1e-30), 0, 1)
    L = P_out[j] + t * (P_out[j+1] - P_out[j])
```

Edge cases: L = 0 gives 0 and L = 1 gives 1 (j = 19, t = 1), so the function is continuous at both ends. L < 0, L > 1 (for example greys brighter than the reference white) and NaN pass through unchanged. The guard 1e-30 and the clip of t only absorb rounding.

**GF6. Lightness-gated hue rotation.** Rotates (a, b) by an angle that depends on the hue and on the corrected L. The chroma does not change.

```text
C = sqrt(a**2 + b**2)
if C < 1e-12: leave a and b unchanged and go to GF7
h  = atan2(b, a)
u  = clip((L - enrichment.L_lo) / (enrichment.L_hi - enrichment.L_lo), 0, 1)
g  = sin(pi * u)**2
d  = h - h_c
e  = fmod(d + pi, 2 * pi)
if e < 0: e = e + 2 * pi
d  = e - pi
G  = exp(-0.5 * (d / enrichment.sigma)**2)
h2 = h + enrichment.amp * g * G
a  = C * cos(h2)
b  = C * sin(h2)
```

L is the value after GF5. `enrichment.L_hi - enrichment.L_lo` = 0.63.

Edge cases:
- C < 1e-12 (greys and colours very near grey): a and b are kept as they are.
- L <= `enrichment.L_lo` gives g = 0 and h2 = h. a and b are still recomputed from C and h2, which can change their last bit.
- L >= `enrichment.L_hi` gives u = 1 and g = `sin(pi)**2` = 1.4997597826618576e-32, not 0.
- The wrap gives d in [-pi, pi]. e can round up to exactly 2 * pi, which gives d = pi.
- h2 is not wrapped and can lie outside [-pi, pi].

**GF7. Chroma compression.** Multiplies a and b by `C**(chroma_power - 1)`, so that the chroma becomes about `C**chroma_power`. The hue does not change.

```text
C = sqrt(a**2 + b**2 + 1e-30)
q = C ** (chroma_power - 1)
a = a * q
b = b * q
```

Edge cases: 1e-30 is added for every colour, not only near zero. `chroma_power - 1` = -0.02200000000000002. a = b = 0 stays 0. For `sqrt(a**2 + b**2)` below about 1e-14 the constant is not negligible, and the resulting a and b stay below 3e-14 in magnitude.

The result of GF7 is the output (L, a, b). For grey input with 0 <= Y <= 10, |a| and |b| are at most 2.9e-15 (reference code, 20,001 greys with Y evenly spaced from 0 to 10, one call of shape (1, 3) per grey; `core/grey-10` gives |a| = 2.8e-15).

#### 5.3.4 Inverse: Lab to XYZ

Input: (L, a, b). Output: XYZ. The steps GI1 to GI7 undo GF7 to GF1 in reverse order. GI2 is the only iterative step. GI5 repeats GF3 and does not undo it. LMS and XYZ are not clipped.

```text
(L, a, b) -> GI1 chroma -> GI2 hue rotation -> GI3 lightness -> GI4 M2_inv
          -> GI5 neutral blend -> GI6 cubic -> GI7 M1_inv -> XYZ
```

**GI1. Undo the chroma compression.** Multiplies a and b by `C**(1/chroma_power) / C`.

```text
C  = sqrt(a**2 + b**2 + 1e-30)
Cr = C ** (1 / chroma_power)
a  = a * (Cr / C)
b  = b * (Cr / C)
```

Edge cases: `1 / chroma_power` = 1.0224948875255624. a = b = 0 stays 0. Because of the 1e-30 term, GI1 is not the exact inverse of GF7 when `a**2 + b**2` is not much larger than 1e-30; the difference in a and b is then below 1e-17.

**GI2. Undo the hue rotation.** Finds the hue h that GF6 maps to the input hue, with 8 steps of Halley's method. L is the input L, which equals the L that GF6 used.

```text
C = sqrt(a**2 + b**2)
if C < 1e-12: leave a and b unchanged and go to GI3
ht = atan2(b, a)
u  = clip((L - enrichment.L_lo) / (enrichment.L_hi - enrichment.L_lo), 0, 1)
g  = sin(pi * u)**2
s2 = enrichment.sigma**2
A  = enrichment.amp * g
h  = ht
repeat exactly 8 times:
    d   = h - h_c
    e   = fmod(d + pi, 2 * pi)
    if e < 0: e = e + 2 * pi
    d   = e - pi
    G   = exp(-0.5 * (d / enrichment.sigma)**2)
    F   = h + A * G - ht
    Fp  = 1 + A * (G * (-d / s2))
    Fpp = A * (G * (-1 / s2 + d * d / (s2 * s2)))
    D   = 2 * Fp * Fp - F * Fpp
    if abs(D) < 1e-30: D = 1
    h   = h - 2 * F * Fp / D
a = C * cos(h)
b = C * sin(h)
```

`s2` = 0.48999999999999994. The start value is ht. The update is the Halley step for `F(h) = h + A * G(h) - ht = 0`, where Fp and Fpp are the first and second derivatives of F. There is no convergence test: the loop always runs 8 times. Stopping early when a step leaves h unchanged gives the same result.

Edge cases: for these parameters `abs(A * dG/dh) <= enrichment.amp / (enrichment.sigma * sqrt(exp(1)))` = 0.0503, so Fp lies between 0.949 and 1.051 and F has one root. The root can differ from the hue of GF6 by a multiple of 2 * pi, because atan2 wraps h2; cos and sin are not affected. C < 1e-12 and L <= `enrichment.L_lo` behave as in GF6.

**GI3. Undo the lightness correction.** Maps L through the inverse piecewise-linear function when 0 <= L <= 1, and leaves L unchanged otherwise.

```text
if 0 <= L <= 1:
    j = the largest index in 0 .. n with P_out[j] <= L
    t = clip((L - P_out[j]) / max(P_out[j+1] - P_out[j], 1e-30), 0, 1)
    L = P_in[j] + t * (P_in[j+1] - P_in[j])
```

Edge cases: the range test uses P_out[0] = 0 and P_out[n+1] = 1, so it is the same test as in GF5. L = 1 gives 1. Outside [0, 1] and for NaN, L passes through.

**GI4. Matrix M2_inv.** Multiplies (L, a, b) by `M2_inv`.

```text
(c0, c1, c2) = M2_inv · (L, a, b)
```

Edge cases: none. The components can be negative here.

**GI5. Neutral blend.** Applies the formula of GF3 to (c0, c1, c2) again.

```text
m  = (c0 + c1 + c2) / 3
r  = (max(c0, c1, c2) - min(c0, c1, c2)) / max(abs(m), 1e-30)
w  = exp(-((r / 1e-5)**2))
ci = ci + w * (m - ci)           for i = 0, 1, 2
```

Edge cases: this is not the inverse of GF3. It changes the components when r is below about 2.73e-4 (Section 5.3.5). m can be negative; r uses abs(m).

**GI6. Inverse transfer.** Evaluates the cubic.

```text
LMS[i] = ci**3 + alpha * ci      for i = 0, 1, 2
```

Edge cases: there is no clip. A negative ci gives a negative LMS[i].

**GI7. Matrix M1_inv.** Multiplies the responses by `M1_inv`.

```text
XYZ = M1_inv · LMS
```

#### 5.3.5 Out-of-range input and round trip

- Negative LMS: GF1 sets negative components of `M1 · XYZ` to 0. This loses information. Among 20,000 XYZ triples drawn uniformly from [-0.05, 1.2]^3 (`numpy.random.default_rng(2026)`), 566 give a negative component of `M1 · XYZ`; for them XYZ to Lab to XYZ changes XYZ by up to 0.227 (largest channel difference, reference code).
- Negative Y: there is no separate rule; GF1 decides. XYZ = (-0.05, -0.05, -0.05) gives (0, 0, 0).
- Small Y: the guard of GF3 moves colours with abs(m) below 1e-30 towards grey (Section 5.3.3). For the sRGB primaries this starts below Y ≈ 1e-35.
- Large Y: there is no upper clip. GF5 acts only for L in [0, 1]. For L >= 1 the gate g of GF6 is 1.5e-32, not 0. Grey L is 0.9999999999999998 at Y = 1 and 2.1663497793510853 at Y = 10 (a batch call with several rows gives 2.166349779351085; Appendix A.2). GF2 overflows for components above about 3e192.
- Infinite input: −inf components of `M1 · XYZ` are clipped to 0, so XYZ = (−inf, 0, 0) gives (0, 0, 0). A +inf component, and GF2 overflow, give (NaN, NaN, NaN). Because `M1[0][2]` < 0, this includes XYZ = (0, 0, −inf). The inverse returns (NaN, NaN, NaN) for any infinite input.
- Round trip away from grey: for input whose components of `M1 · XYZ` are all positive and whose spread r in GF3 is above about 2.73e-4, XYZ to Lab to XYZ returns the input up to rounding. On the 19,434 triples of the set above without a negative component, the largest channel difference is 3.1e-15.
- Round trip near grey: when r is below about 2.73e-4, GF3 reduces the spread and GI5 reduces it again. XYZ to Lab to XYZ then changes the colour. In a sweep of 3,630 near-greys (Y = 0.01, 0.1, 0.5, 1 and 3; each XYZ component in turn multiplied by 1 + q and by 1 − q, for 121 values of q spaced evenly on a log scale from 1e-9 to 1e-3), the largest channel difference is 2.5e-5 × Y. sRGB white, (0.9504700000000421, 1.0000000102999, 1.088830000000039), returns with a difference of 5.8e-9. Lab to XYZ to Lab changes such colours in the same way.

#### 5.3.6 Stages of the reference code that Helmlab 2026 does not use

The reference code contains further stages. Each one is skipped because its parameters are absent from `gen_params.json` and default to the values given here. An implementation of Helmlab 2026 does not need them.

- Power transfer `x**gamma[i]` (`gamma` = [1/3, 1/3, 1/3]) and softened cube root (`softcbrt_eps` = 0.001): replaced by GF2 because `transfer` is `"depcubic"`.
- Hue-dependent lightness correction (`hue_L_amp` = 0; `hue_L_center`, `hue_L_width`, `hue_L_knee` are then not used).
- Fourier hue correction (`hue_cos1`, `hue_sin1`, ..., `hue_cos4`, `hue_sin4` = 0).
- Cubic lightness correction (`L_corr_p1`, `L_corr_p2`, `L_corr_p3` = 0). A non-empty `L_corr_pw` also takes precedence over it.
- Dark lightness compression (`lp_dark`, `lp_dark_hcos`, `lp_dark_hsin` = 0).
- Lightness-dependent chroma scale (`lc1`, `lc2` = 0).
- Neutral correction: the constructor argument `neutral_correction` is False in the Python class, in `hl.gen` and in the JavaScript `GenSpace`.

### 5.4 GenSpace LCh (`helmgenlch`)

Lab to LCh, as `hl.gen.to_lch` (`to_lch` in `src/helmlab/helmlab.py`). L is unchanged.

```text
C = hypot(a, b)
h = atan2(b, a) * (180 / pi)
e = fmod(h, 360)
if e < 0: e = e + 360
h = e
```

`180 / pi` = 57.29577951308232. A result of -0 becomes +0.

Edge cases: h lies in [0, 360]. It equals 360 when atan2 returns a negative angle of magnitude below about 5e-16, because e + 360 then rounds to 360. When a = b = 0, h is 0 if a is +0 and 180 if a is -0. Grey colours have |a| and |b| near 1e-15 from rounding, so their h has no meaning. In general, h in degrees magnifies a difference in a and b by (180/pi)/C: between a one-row call and a batch call of the reference, h differs by up to 4e-10 degrees when C is between 1e-4 and 1e-3. Section 8 therefore compares the hue as the chord C·Δh.

LCh to Lab, as `hl.gen.from_lch` (`from_lch` in the same file). L is unchanged.

```text
a = C * cos(h * (pi / 180))
b = C * sin(h * (pi / 180))
```

Edge cases: h is not wrapped or checked, and any finite h works. C is not checked; a negative C gives the hue h + 180.

### 5.5 MetricSpace, base form (`helmlab-metric`)

The base form maps CIE XYZ (Section 4) to coordinates (L, a, b). It is the Python class `MetricSpace` (`src/helmlab/spaces/metric.py`) with these settings:

- neutral correction off (`neutral_correction=False`, the default of the class);
- the rotation by `display_phi_deg` = −28.2 applied (Step F14);
- surround S = 0.5, the default (Section 5.5.1).

All parameters come from `params/metric_params.json`. In the formulas each parameter appears under its JSON name. Section 5.5.5 gives the role of every entry. Every step gives its results new names. A name on the right-hand side of a formula is a value computed in an earlier step.

#### 5.5.1 Surround

The Python class takes a surround value S, with default `surround=0.5`. The Python class `Helmlab` clips S to [0, 1] and also defaults to 0.5. The code uses `dS = S − 0.5`. The JavaScript package has no surround input and does not read the 11 surround coefficients.

The 11 surround coefficients enter the forward and the inverse as follows:

```text
Steps F5, I10:  weight   = hk_weight + hk_weight_S*dS                    always
                exponent = clip(hk_power + hk_power_S*dS, 0.01, 10.0)    always
                F5       = F5 + hk_hue_S*dS*cos(h5)                      only if hk_hue_S ≠ 0
Steps F7, I8:   k        = k + lp_dark_S*dS + lp_dark_S2*dS**2           always
Steps F8, I7:   S8       = S8 * exp(cs_S_lin*dS + cs_S_quad*dS**2)       only if cs_S_lin ≠ 0 or cs_S_quad ≠ 0
Steps F10, I5:  T10      = T10 * exp(lc_S_lin*dS + lc_S_quad*dS**2)      only if lc_S_lin ≠ 0 or lc_S_quad ≠ 0
Steps F12, I3:  G12      = G12 * exp(hl_S_lin*dS)                        only if hl_S_lin ≠ 0
Steps F13, I2:  L3       = L3 + L_S_offset*dS    (inverse: minus)        only if L_S_offset ≠ 0
```

All 11 coefficients are 0. The terms marked "only if" are therefore skipped. The terms marked "always" add `0*dS` or `0*dS**2`. For S in [0, 1], and for any S with |S − 0.5| below about 1.3e154, this is +0 or −0, and adding it leaves the value unchanged. At the default S = 0.5, dS = 0 as well. The base form therefore does not depend on S, and the steps below omit the surround terms. (For larger |S − 0.5|, `dS**2` overflows and the Python class raises `OverflowError`.)

#### 5.5.2 Forward: XYZ to (L, a, b)

The input is XYZ = (X, Y, Z). The reference code skips a stage when all of its coefficients are 0. In `metric_params.json` every stage of Steps F4 to F12 has non-zero coefficients, so all of them run. The labels in parentheses are the stage numbers used in the comments of `metric.py`.

**Step F1. Cone-like response (stage 1).** A linear map gives three channel values.

```text
LMS = M1 · XYZ
```

Negative XYZ and negative LMS values are allowed. Nothing is clipped.

**Step F2. Compression (stage 2).** Each channel is raised to its own exponent, keeping its sign.

```text
LMS_c[i] = sign(LMS[i]) * abs(LMS[i]) ** gamma[i]        for i = 0, 1, 2
```

A channel equal to 0 gives 0.

**Step F3. Opponent coordinates (stage 3).** A second linear map gives the raw coordinates.

```text
(L0, a0, b0) = M2 · LMS_c
```

**Step F4. Hue shift (stage 3.5).** The point (a0, b0) is rotated by an angle that depends on its hue. Chroma is kept.

```text
h1    = atan2(b0, a0)
delta = hue_cos1*cos(h1) + hue_sin1*sin(h1) + hue_cos2*cos(2*h1) + hue_sin2*sin(2*h1)
        + hue_cos3*cos(3*h1) + hue_sin3*sin(3*h1) + hue_cos4*cos(4*h1) + hue_sin4*sin(4*h1)
C1    = sqrt(a0**2 + b0**2)
a1    = C1 * cos(h1 + delta)
b1    = C1 * sin(h1 + delta)
```

`h1 + delta` is not wrapped. When a0 = b0 = 0, C1 = 0 and a1 = b1 = 0. The function `h + delta(h)` is strictly increasing: its derivative lies between 0.41 and 1.72. Step I11 relies on this.

**Step F5. Chroma-dependent lightness (stage 3.7, "H-K" in the code).** A term that grows with chroma is added to L0. Chroma and hue are recomputed from a1 and b1; C5 can differ from C1 in the last bits. The hue h5 is used again in Steps F6, F7 and F8.

```text
C5   = sqrt(a1**2 + b1**2)
h5   = atan2(b1, a1)
P5   = clip(hk_power, 0.01, 10.0)
F5   = 1 + (hk_hue_mod*cos(h5) + hk_sin1*sin(h5) + hk_cos2*cos(2*h5) + hk_sin2*sin(2*h5))
L_hk = L0 + (hk_weight * C5 ** P5) * F5
```

`hk_power` lies inside [0.01, 10.0], so P5 = `hk_power`. When C5 = 0 the added term is 0. F5 is not clipped and can be negative; over all hues it ranges from about −0.41 to 3.08.

**Step F6. Lightness correction (stage 4).** A polynomial in t = L(1 − L), with one hue-dependent coefficient, corrects the lightness.

```text
t  = L_hk * (1 - L_hk)
L1 = L_hk + L_corr_p1*t + L_corr_p2*t*(0.5 - L_hk) + L_corr_p3*t*t + t*(Lh_cos1*cos(h5) + Lh_sin1*sin(h5))
```

Nothing is clipped. For L_hk outside [0, 1], t is negative.

**Step F7. Dark-lightness compression (stage 4.5).** A lightness below 1 is multiplied by a factor that depends on the lightness and the hue.

```text
k  = lp_dark + lp_dark_hcos*cos(h5) + lp_dark_hsin*sin(h5)
g  = k * L1 * max(0, 1 - L1) ** 2
L2 = L1 * exp(clip(g, -30, 30))
```

For finite L1 ≥ 1, g = 0 and L2 = L1. The clip to [−30, 30] has no effect for 0 ≤ L1 ≤ 1, where |g| < 0.21. If `t*t` overflows in Step F6 (|t| above about 1.3e154, for example for a grey with Y above about 1e150), L1 = +inf, g is NaN and the output is (NaN, NaN, NaN).

**Step F8. Hue-dependent chroma scale (stage 5).** a1 and b1 are multiplied by a positive factor that depends on h5.

```text
S8 = exp(cs_cos1*cos(h5) + cs_sin1*sin(h5) + cs_cos2*cos(2*h5) + cs_sin2*sin(2*h5)
         + cs_cos3*cos(3*h5) + cs_sin3*sin(3*h5) + cs_cos4*cos(4*h5) + cs_sin4*sin(4*h5))
a2 = a1 * S8
b2 = b1 * S8
```

Nothing is clipped.

**Step F9. Chroma power (stage 5.5).** Chroma is raised to a hue-dependent exponent. The hue is taken again, from a2 and b2.

```text
h9 = atan2(b2, a2)
C9 = sqrt(a2**2 + b2**2)
p9 = 1 + (cp_cos1*cos(h9) + cp_sin1*sin(h9) + cp_cos2*cos(2*h9) + cp_sin2*sin(2*h9))
D9 = C9 ** p9        if C9 > 0
D9 = 0               otherwise
a3 = D9 * cos(h9)
b3 = D9 * sin(h9)
```

Zero chroma stays 0. Over all hues the exponent p9 lies between about 0.66 and 1.26.

**Step F10. Lightness-dependent chroma scale (stage 6).** a3 and b3 are multiplied by a factor that depends on L2, the lightness after Step F7.

```text
u   = L2 - 0.5
T10 = exp(clip(lc1*u + lc2*u**2, -30, 30))
a4  = a3 * T10
b4  = b3 * T10
```

**Step F11. Hue-lightness chroma factor (stage 6.5).** a4 and b4 are multiplied by a factor that depends on the hue and on L2.

```text
h11 = atan2(b4, a4)
q11 = hlc_cos1*cos(h11) + hlc_sin1*sin(h11) + hlc_cos2*cos(2*h11) + hlc_sin2*sin(2*h11)
H11 = exp(clip((L2 - 0.5) * q11, -30, 30))
a5  = a4 * H11
b5  = b4 * H11
```

**Step F12. Hue-dependent lightness scale (stage 8).** L2 is multiplied by a positive factor that depends on the hue.

```text
h12 = atan2(b5, a5)
G12 = exp(hl_cos1*cos(h12) + hl_sin1*sin(h12) + hl_cos2*cos(2*h12) + hl_sin2*sin(2*h12))
L3  = L2 * G12
```

Nothing is clipped. The code has no stage 7.

**Step F13. Stages 9 and 10.** Neither changes the base form. Stage 9 adds `L_S_offset*dS` and is skipped because `L_S_offset` = 0 (Section 5.5.1). Stage 10 is the neutral correction, which is off in the base form (Appendix B).

**Step F14. Rotation (stage 11).** The (a, b) plane is rotated by the fixed angle `display_phi_deg`.

```text
phi = display_phi_deg * (pi / 180)
cr  = cos(phi)
sr  = sin(phi)
L   = L3
a   = a5*cr - b5*sr
b   = a5*sr + b5*cr
```

`pi / 180` is computed first. The values are phi = −0.4921828490624009 (`-0x1.f7fec7e1801e7p-2`), cr = 0.8813034520649923 (`0x1.c33a34c0f12cep-1`) and sr = −0.47255076486905395 (`-0x1.e3e459033fde9p-2`). The code applies this step only when phi ≠ 0. The rotation keeps C = sqrt(a² + b²).

The output is (L, a, b).

**Zero chroma at black.** XYZ = (0, 0, 0) gives a0 = b0 = +0. Every later a and b is then 0, every hue in Steps F4 to F12 is atan2(+0, +0) = 0, and the output is exactly (0, 0, 0).

**The zero-chroma curve.** A grey with Y > 0 does not have zero chroma in the base form (Section 3). a0 = b0 = 0 holds instead where LMS_c is a positive multiple of the cross product of the rows M2[1] and M2[2]:

```text
LMS_c = s * (2.5633498896562767, 3.0420833463457604, 4.462733874281098)        s > 0
```

These colours are greyish blues with chromaticity y ≈ 0.276 and x from about 0.18 (Y = 0.001) to 0.27 (Y = 1). For Y from about 0.01 to 0.62 they are inside the sRGB gamut. On the curve h1 is undefined, but L still depends on the hue, through Lh in Step F6, k in Step F7 and G12 in Step F12. The base form is therefore not continuous there: L tends to different values depending on the direction from which a point of the curve is approached. Over 360 directions the limits span 0.025 at Y = 0.007, 0.092 at Y = 0.077, 0.153 at Y = 0.31, 0.196 at Y = 0.49 and 0.244 at Y = 0.70.

On the curve itself a0 and b0 are rounding residue, and the reference returns different L for the same input: XYZ = (0.06564918031540655, 0.07736666851128433, 0.13699287588355347) gives L = 0.4453387746570464 in a one-row call and 0.4085666256969879 as one row of a 64-row call. Near the curve L depends on the rounding of a0 and b0. With C0 = sqrt(a0**2 + b0**2), a one-row call and a 64-row call differ in L by up to about 2e-12 at C0 = 1e-4, 1e-10 at C0 = 1e-6 and 1e-9 at C0 = 1e-7 (largest differences in samples of 20,000 colours of the domain at each C0, each colour at a random row of the 64-row call; these are not bounds). Among the 16,777,216 8-bit sRGB colours the smallest C0 is 7.17e-5 (`#456d85`). The test vectors keep C0 ≥ 1e-4 (Section 9).

#### 5.5.3 Inverse: (L, a, b) to XYZ

The inverse runs the stages in reverse order. Step In undoes Step F(15 − n). Steps I8, I9 and I11 use Newton iteration with a fixed number of iterations. The hues in Steps I3 to I8 are taken from the values at hand, which can differ in the last bits from the forward values.

For XYZ from the domain the inverse returns the input to rounding. The largest channel difference of XYZ → Lab → XYZ is 5.6e-15 over 400,000 random 8-bit sRGB colours, and 1.0e-13 over 200,000 random XYZ in [0, 1]³, where L reaches 10 (`numpy.random.default_rng(11)`, the sRGB colours drawn first; one batch call per set). Greys return to a relative 1e-12 up to Y ≈ 17.4. From Y ≈ 20.3 the relative error is above 1e-6, and Y = 50 returns Y = 196.9.

For a triple that is not an output of the forward, the inverse can return an XYZ whose forward is a different triple. At C = 0.05 and C = 0.2 (72 hues; L in steps of 0.01), forward(inverse(L, a, b)) returns the triple to 1e-12 for every hue when L is from −0.16 to 1.55 (C = 0.05) or to 1.61 (C = 0.2); it fails for some hues at L = −0.17, at 1.56 (C = 0.05) and at 1.62 (C = 0.2). For example, (−0.25, 0.05, 0) maps to XYZ (−1.818, −1.847, −3.036), whose forward is (5.346, 0, 0) up to rounding.

**Step I1. Undo the rotation (undoes F14).**

```text
a5 = a*cr + b*sr
b5 = -a*sr + b*cr
L3 = L
```

cr and sr are the values of Step F14. The step is skipped when phi = 0.

**Step I2. Stages 10 and 9 (undoes F13).** Nothing happens in the base form. In the neutral-corrected form the neutral correction is added back here (Appendix B).

**Step I3. Undo the hue-dependent lightness scale (undoes F12).**

```text
h12 = atan2(b5, a5)
G12 = exp(hl_cos1*cos(h12) + hl_sin1*sin(h12) + hl_cos2*cos(2*h12) + hl_sin2*sin(2*h12))
L2  = L3 / G12
```

**Step I4. Undo the hue-lightness chroma factor (undoes F11).** The hue is h12 of Step I3.

```text
q11 = hlc_cos1*cos(h12) + hlc_sin1*sin(h12) + hlc_cos2*cos(2*h12) + hlc_sin2*sin(2*h12)
H11 = exp(clip((L2 - 0.5) * q11, -30, 30))
a4  = a5 / H11
b4  = b5 / H11
```

**Step I5. Undo the lightness-dependent chroma scale (undoes F10).**

```text
u   = L2 - 0.5
T10 = exp(clip(lc1*u + lc2*u**2, -30, 30))
a3  = a4 / T10
b3  = b4 / T10
```

**Step I6. Undo the chroma power (undoes F9).**

```text
h9 = atan2(b3, a3)
C  = sqrt(a3**2 + b3**2)
p9 = 1 + (cp_cos1*cos(h9) + cp_sin1*sin(h9) + cp_cos2*cos(2*h9) + cp_sin2*sin(2*h9))
D  = C ** (1 / p9)   if C > 0
D  = 0               otherwise
a2 = D * cos(h9)
b2 = D * sin(h9)
```

**Step I7. Undo the hue-dependent chroma scale (undoes F8).**

```text
h8 = atan2(b2, a2)
S8 = exp(cs_cos1*cos(h8) + cs_sin1*sin(h8) + cs_cos2*cos(2*h8) + cs_sin2*sin(2*h8)
         + cs_cos3*cos(3*h8) + cs_sin3*sin(3*h8) + cs_cos4*cos(4*h8) + cs_sin4*sin(4*h8))
a1 = a2 / S8
b1 = b2 / S8
```

**Step I8. Undo the dark-lightness compression (undoes F7).** Newton iteration on `x * exp(g(x)) = L2`. The hue h5 is used again in Steps I9, I10 and I11.

```text
h5 = atan2(b1, a1)
k  = lp_dark + lp_dark_hcos*cos(h5) + lp_dark_hsin*sin(h5)
x  = L2                                        start value
repeat 12 times:
    m  = max(0, 1 - x)
    g  = k * x * m ** 2
    e  = exp(clip(g, -30, 30))
    f  = x*e - L2
    gd = k * m * (1 - 3*x)
    fd = e * (1 + x*gd)
    if abs(fd) < 1e-10: fd = 1
    x  = x - f / fd
L1 = x
```

The update runs exactly 12 times; there is no other stop rule. The derivative gd ignores the clip. JavaScript 1.0.1 and 1.0.2 stop as soon as an update leaves x unchanged, which gives the same result as running all 12 of their own updates; JavaScript 1.0.0 runs all 12. The JavaScript update expressions are associated differently from the ones above, so its iterates can differ in the last bits (Appendix A.3).

**Step I9. Undo the lightness correction (undoes F6).** Newton iteration on the polynomial of Step F6.

```text
Lh = Lh_cos1*cos(h5) + Lh_sin1*sin(h5)
x  = L1                                        start value
repeat 15 times:
    t  = x * (1 - x)
    dt = 1 - 2*x
    f  = x + L_corr_p1*t + L_corr_p2*t*(0.5 - x) + L_corr_p3*t*t + t*Lh - L1
    fd = 1 + (L_corr_p1 + Lh)*dt + L_corr_p2*(dt*(0.5 - x) - t) + L_corr_p3*2*t*dt
    if abs(fd) < 1e-10: fd = 1
    x  = x - f / fd
L_hk = x
```

The update runs exactly 15 times; there is no other stop rule. The JavaScript package stops early in the same way as in Step I8.

**Step I10. Undo the chroma-dependent lightness (undoes F5).**

```text
C5 = sqrt(a1**2 + b1**2)
P5 = clip(hk_power, 0.01, 10.0)
F5 = 1 + (hk_hue_mod*cos(h5) + hk_sin1*sin(h5) + hk_cos2*cos(2*h5) + hk_sin2*sin(2*h5))
L0 = L_hk - (hk_weight * C5 ** P5) * F5
```

**Step I11. Undo the hue shift (undoes F4).** Newton iteration on `x + delta(x) = h5`, where delta is the function of Step F4 and delta_d is its derivative.

```text
delta(x)   = hue_cos1*cos(x) + hue_sin1*sin(x) + hue_cos2*cos(2*x) + hue_sin2*sin(2*x)
             + hue_cos3*cos(3*x) + hue_sin3*sin(3*x) + hue_cos4*cos(4*x) + hue_sin4*sin(4*x)
delta_d(x) = -hue_cos1*sin(x) + hue_sin1*cos(x) - 2*hue_cos2*sin(2*x) + 2*hue_sin2*cos(2*x)
             - 3*hue_cos3*sin(3*x) + 3*hue_sin3*cos(3*x) - 4*hue_cos4*sin(4*x) + 4*hue_sin4*cos(4*x)
C1 = sqrt(a1**2 + b1**2)
x  = h5                                        start value
repeat 8 times:
    f  = x + delta(x) - h5
    fd = 1 + delta_d(x)
    if abs(fd) < 1e-10: fd = 1
    x  = x - f / fd
a0 = C1 * cos(x)
b0 = C1 * sin(x)
```

The update runs exactly 8 times; there is no other stop rule. The JavaScript package stops early in the same way as in Step I8. x is not wrapped. Because `1 + delta_d` is at least 0.41, the guard does not trigger for finite input.

**Step I12. Undo the opponent map (undoes F3).**

```text
LMS_c = M2_inv · (L0, a0, b0)
```

**Step I13. Undo the compression (undoes F2).**

```text
LMS[i] = sign(LMS_c[i]) * abs(LMS_c[i]) ** (1 / gamma[i])        for i = 0, 1, 2
```

`1 / gamma[i]` is rounded to binary64 before it is used as the exponent.

**Step I14. Undo the cone-like response (undoes F1).**

```text
XYZ = M1_inv · LMS
```

M1_inv and M2_inv are the inverses of M1 and M2 in binary64. The Python reference computes them once with `numpy.linalg.inv`; the JavaScript package uses Cramer's rule. Python gives:

```text
M1_inv = [[ 1.065107295808859,    -0.3150044075301121,  0.1803492381741039  ],
          [ 0.4721107713837796,    0.4271995765962455,  0.04150680489380985 ],
          [-0.09156391926309541,  -0.1736709363194979,  0.7719790382558294  ]]
M2_inv = [[ 0.9183897822815021,    0.5232051237088666,  0.7065804598090856  ],
          [ 1.0899090574433026,    0.07005324849041904, 0.4319776874787044  ],
          [ 1.598895729264209,     0.22061850068770233, 0.012506037355220951]]
```

Zero chroma: for an input with a = b = 0, Step I1 gives a5 and b5 equal to +0 or −0, with signs set by the arithmetic of Step I1. (a, b) = (+0, +0), (+0, −0) and (−0, −0) give hue +0 or −0 in Steps I3 to I8, and the same XYZ. (a, b) = (−0, +0) gives a5 = −0 and b5 = +0, so every atan2 in Steps I3 to I8 returns pi (Section 5.1), and the XYZ differs: at L = 0.9, (+0, +0) gives XYZ (0.4858091917200994, 0.5194949674376206, 0.8723817291333001) and (−0, +0) gives (1.7462612412095702, 1.760488782877648, 2.8468597479077937). In every case a0 = b0 = 0, so the result lies on the zero-chroma curve of Section 5.5.2, and its forward does not return the input: the XYZ from (0.5, +0, +0) has forward L = 0.5291.

#### 5.5.4 Reference entry points

- Python `MetricSpace().from_XYZ(XYZ)` and `MetricSpace().to_XYZ(lab)` give the base form. They take arrays of shape (..., 3). A result can depend on the number of rows in the call: in the last bits in general, and by more near the zero-chroma curve (Section 5.5.2). The vectors use one-row calls (Appendix A.2).
- The JavaScript class `AnalyticalSpace`, also exported as `MetricSpace`, gives the base form only with the option `{ neutralCorrection: false }`. Its default is the neutral-corrected form.
- The methods `base_from_XYZ` and `base_to_XYZ` (Python) and `baseFromXYZ`, `baseToXYZ` and `baseNeutralError` (JavaScript) compute only Steps F1 to F3 or a separate grey table. They are not the base form of this definition.

#### 5.5.5 Role of each entry of `metric_params.json`

The file has 76 entries: `M1` and `M2` (3×3 each), `gamma` (3 values) and 73 scalars.

| Entry | Steps | Role |
|---|---|---|
| `M1` | F1; I14 through M1_inv | XYZ to LMS matrix |
| `gamma` | F2, I13 | Exponent of each LMS channel |
| `M2` | F3; I12 through M2_inv | Compressed LMS to (L0, a0, b0) |
| `hue_cos1`, `hue_sin1`, `hue_cos2`, `hue_sin2`, `hue_cos3`, `hue_sin3`, `hue_cos4`, `hue_sin4` | F4, I11 | Fourier coefficients of the hue shift delta, harmonics 1 to 4 |
| `hk_weight` | F5, I10 | Weight of the chroma-dependent lightness term |
| `hk_power` | F5, I10 | Exponent of chroma in that term |
| `hk_hue_mod`, `hk_sin1`, `hk_cos2`, `hk_sin2` | F5, I10 | Coefficients of cos(h), sin(h), cos(2h), sin(2h) in the hue factor F5 |
| `L_corr_p1`, `L_corr_p2`, `L_corr_p3` | F6, I9 | Coefficients of t, t·(0.5 − L) and t² in the lightness correction |
| `Lh_cos1`, `Lh_sin1` | F6, I9 | Hue-dependent coefficient of t |
| `lp_dark`, `lp_dark_hcos`, `lp_dark_hsin` | F7, I8 | Coefficient k of the dark-lightness compression |
| `cs_cos1`, `cs_sin1`, `cs_cos2`, `cs_sin2`, `cs_cos3`, `cs_sin3`, `cs_cos4`, `cs_sin4` | F8, I7 | Logarithm of the hue-dependent chroma scale, harmonics 1 to 4 |
| `cp_cos1`, `cp_sin1`, `cp_cos2`, `cp_sin2` | F9, I6 | Hue-dependent part of the chroma exponent p9 |
| `lc1`, `lc2` | F10, I5 | Logarithm of the lightness-dependent chroma scale: coefficients of (L − 0.5) and (L − 0.5)² |
| `hlc_cos1`, `hlc_sin1`, `hlc_cos2`, `hlc_sin2` | F11, I4 | Hue factor q11, which multiplies (L − 0.5) in the logarithm of the chroma factor |
| `hl_cos1`, `hl_sin1`, `hl_cos2`, `hl_sin2` | F12, I3 | Logarithm of the hue-dependent lightness scale |
| `display_phi_deg` | F14, I1 | Rotation of the (a, b) plane, in degrees |
| `dist_sl` | D3 | Coefficient of the lightness weight SL |
| `dist_sc` | D4 | Coefficient of the chroma weight SC |
| `dist_wC` | D5 | Weight of the chromatic term |
| `dist_power` | D5 | The weighted sum is raised to `dist_power / 2` |
| `dist_compress` | D6 | Compression constant |
| `dist_post_power` | D7 | Final exponent |
| `dist_linear` | D6 | 0; drops out (Section 5.6) |
| `dist_nl`, `dist_sat` | — | 0; not used (Section 5.6) |
| `dist_sl_hcos1`, `dist_sl_hsin1`, `dist_sl_hcos2`, `dist_sl_hsin2` | D3 | 0; drop out (Section 5.6) |
| `dist_sc_hcos1`, `dist_sc_hsin1`, `dist_sc_hcos2`, `dist_sc_hsin2` | D4 | 0; drop out (Section 5.6) |
| `hk_weight_S`, `hk_power_S`, `hk_hue_S`, `lp_dark_S`, `lp_dark_S2`, `cs_S_lin`, `cs_S_quad`, `lc_S_lin`, `lc_S_quad`, `hl_S_lin`, `L_S_offset` | F5, F7, F8, F10, F12, F13; I10, I8, I7, I5, I3, I2 | 0; surround terms that drop out (Section 5.5.1) |

### 5.6 ΔE `helmlab`

ΔE `helmlab` is computed from two base-form triples (L_1, a_1, b_1) and (L_2, a_2, b_2). The Python reference is `MetricSpace.distance_from_lab`. The suffixes _1 and _2 refer to the two colours.

**Step D1. Differences.**

```text
dL2  = (L_1 - L_2) ** 2
dab2 = (a_1 - a_2) ** 2 + (b_1 - b_2) ** 2
```

**Step D2. Means.**

```text
Lm  = (L_1 + L_2) * 0.5
C_1 = sqrt(a_1**2 + b_1**2)
C_2 = sqrt(a_2**2 + b_2**2)
Cm  = (C_1 + C_2) * 0.5
```

**Step D3. Lightness weight.**

```text
SL   = 1 + dist_sl * (Lm - 0.5) ** 2
dL2w = dL2 / SL ** 2
```

For 0 ≤ Lm ≤ 1, SL lies between 0.7711 and 1. SL = 0 where (Lm − 0.5)² = −1/`dist_sl`, that is at Lm ≈ 1.5451 or −0.5451. The code has no guard: a computed SL of exactly 0 divides by 0, and ΔE is NaN.

**Step D4. Chroma weight.**

```text
SC    = 1 + dist_sc * Cm
dab2w = dab2 / SC ** 2
```

**Step D5. Weighted power sum.**

```text
E0 = (dL2w + dist_wC * dab2w) ** (dist_power / 2)
```

**Step D6. Compression.**

```text
E1 = E0 / (1 + dist_compress * E0)
```

**Step D7. Final power.**

```text
ΔE = E1 ** dist_post_power
```

ΔE ≥ 0. For equal triples ΔE = 0, unless SL computes to exactly 0, where ΔE is NaN. For finite E0, ΔE does not exceed a bound and tends to it as E0 grows:

```text
ΔE <= (1 / dist_compress) ** dist_post_power = 0.15003639078432002
```

In binary64 the bound is reached when E0 is large, for example for pairs with SL near 0 (Lm near 1.5451, outside the domain).

These ΔE parameters are 0 and drop out. The code would use them as follows; it skips each term when its parameters are 0.

```text
E1 = E0 * (1 + dist_linear*dist_compress*E0) / (1 + dist_compress*E0)        replaces Step D6 if dist_linear ≠ 0

hm = atan2((b_1 + b_2)*0.5, (a_1 + a_2)*0.5)                                   used if any of the 8 terms below ≠ 0
SL = 1 + (dist_sl + (dist_sl_hcos1*cos(hm) + dist_sl_hsin1*sin(hm) + dist_sl_hcos2*cos(2*hm) + dist_sl_hsin2*sin(2*hm))) * (Lm - 0.5)**2
SC = 1 + (dist_sc + (dist_sc_hcos1*cos(hm) + dist_sc_hsin1*sin(hm) + dist_sc_hcos2*cos(2*hm) + dist_sc_hsin2*sin(2*hm))) * Cm

E1 = E0 * exp(dist_nl * E0 / (1 + dist_sat*E0))     replaces Step D6 only if dist_compress = 0 and dist_nl ≠ 0
                                                    (dist_nl * E0 if dist_sat = 0)
```

For `dist_linear` = 0 the first formula equals Step D6. `dist_compress` is not 0, so `dist_nl` and `dist_sat` are never used.

ΔE uses only dL2, dab2, Lm and the two chromas. A rotation applied to both (a, b) planes changes none of them apart from rounding. ΔE from base-form coordinates therefore equals ΔE from the coordinates before Step F14, to within rounding.

Input checks of the reference code, not part of the definition: `distance_from_lab` raises `ValueError` if an input is not finite or if |L| > 3 for either input. `MetricSpace.distance(XYZ_1, XYZ_2)` computes both base-form triples with the forward (never with the neutral correction) and then ΔE; it warns, without failing, when the second or third component of an input is below −1e-9. `hl.metric.difference(color_a, color_b)` parses both colour strings to XYZ (Section 5.2) and calls `MetricSpace.distance`.

**Range of the fitting data.** The parameters of MetricSpace were fitted with CMA-ES to a weighted sum of STRESS on three datasets: COMBVD (3,813 pairs, weight 0.55), MacAdam 1974 (128 pairs, 0.15) and Helmlab's Human Feedback data (3,552 judgements of 48 colour pairs, made on the observers' own uncalibrated screens, 0.30). COMBVD pairs whose white is not (0.95047, 1, 1.08883) were adapted to that white with a Bradford transform. MacAdam 1974 is specified for D65 and the CIE 1964 10° observer; the fitting took its white as x = 0.31272, y = 0.32903, the D65 white of the 1931 observer, and did not adapt it. The Human Feedback colours are sRGB hex values, converted to XYZ (D65) with a 7-decimal sRGB matrix. He et al. 2022 (82 pairs) was used in model selection. The paper describes these data differently (E-2026-0002). The fitting script and the Human Feedback data are not published. The colour differences of the three fitting datasets:

| Dataset | Rows | ΔE `helmlab`: min, median, max | CIEDE2000: min, median, max |
|---|---|---|---|
| COMBVD | 3,813 | 0.00049, 0.023, 0.117 | 0.024, 1.26, 10.4 |
| MacAdam 1974 | 128 | 0.054, 0.110, 0.133 | 2.9, 8.8, 15.5 |
| Human Feedback | 3,552 | 0, 0.119, 0.148 | 0, 10.7, 85.0 |

ΔE `helmlab` is computed on the XYZ as the fitting used them, and CIEDE2000 on the original XYZ with the white the fitting took for each dataset; the Human Feedback rows include judgements of identical colours. Of the 3,989 colour pairs, 7 have a colour darker than L* 10 (4 in COMBVD, 3 in Human Feedback) and none has a colour lighter than L* 92.8 (Section 10).

The published datasets:

- COMBVD, the combined dataset from which CIEDE2000 was derived: M. R. Luo, G. Cui and B. Rigg, "The development of the CIE 2000 colour-difference formula: CIEDE2000", *Color Research & Application* 26(5), 340–350 (2001), [doi:10.1002/col.1049](https://doi.org/10.1002/col.1049).
- MacAdam 1974: D. L. MacAdam, "Uniform color scales", *Journal of the Optical Society of America* 64(12), 1691–1702 (1974), [doi:10.1364/JOSA.64.001691](https://doi.org/10.1364/JOSA.64.001691).
- He et al. 2022: R. He, K. Xiao, M. Pointer, M. Melgosa and Y. Bressler, "Optimizing parametric factors in CIELAB and CIEDE2000 color-difference formulas for 3D-printed spherical objects", *Materials* 15(12), 4055 (2022), [doi:10.3390/ma15124055](https://doi.org/10.3390/ma15124055); data: [doi:10.5281/zenodo.6502724](https://doi.org/10.5281/zenodo.6502724).

## 6. Parameters

The two JSON files in `params/` are normative. The tables below reproduce them, generated by `tools/spec/param_tables.py`. Each value is given as the shortest decimal that round-trips to the same binary64 value, and as a hexadecimal float.

In `metric_params.json`, the 11 surround coefficients (`hk_weight_S`, `hk_power_S`, `hk_hue_S`, `lp_dark_S`, `lp_dark_S2`, `cs_S_lin`, `cs_S_quad`, `lc_S_lin`, `lc_S_quad`, `hl_S_lin`, `L_S_offset`) are 0, so the surround parameter has no effect. The rotation `display_phi_deg` = −28.2 is part of the base form.

<!-- BEGIN GENERATED by tools/spec/param_tables.py; do not edit by hand -->

GenSpace, `params/gen_params.json`, sha256 `310ad1105eeec0aeb72bb8bf24185cda15f62ea8f61822dceca9d61b8a69356f`

| Parameter | Decimal | Hexadecimal |
|---|---|---|
| `M1[0][0]` | `0.8154374735648701` | `0x1.a1810541c9b7fp-1` |
| `M1[0][1]` | `0.3603221491264266` | `0x1.70f84a1a16f9dp-2` |
| `M1[0][2]` | `-0.12432703417946676` | `-0x1.fd3e581cf3a57p-4` |
| `M1[1][0]` | `0.03298391207546648` | `0x1.0e3446f510881p-5` |
| `M1[1][1]` | `0.9292940788255503` | `0x1.dbcc6ef9d8289p-1` |
| `M1[1][2]` | `0.03614494665290377` | `0x1.2819727946bfcp-5` |
| `M1[2][0]` | `0.048184113668356454` | `0x1.8ab9690c8e078p-5` |
| `M1[2][1]` | `0.26427748135788043` | `0x1.0e9ec18e0160ep-2` |
| `M1[2][2]` | `0.6336388271114471` | `0x1.446c4eefd6c3fp-1` |
| `M2[0][0]` | `0.21193779684470104` | `0x1.b20c7191dfcfcp-3` |
| `M2[0][1]` | `0.7992121834263127` | `0x1.993256dcc2f3dp-1` |
| `M2[0][2]` | `-0.00410075161564345` | `-0x1.0cbf3214056cdp-8` |
| `M2[1][0]` | `2.4672018828033475` | `0x1.3bcd4573a28f7p+1` |
| `M2[1][1]` | `-2.9877348024830788` | `-0x1.7e6e1810e4968p+1` |
| `M2[1][2]` | `0.520532919679731` | `0x1.0a834a75081c2p-1` |
| `M2[2][0]` | `-0.11390787868068575` | `-0x1.d291115b0b640p-4` |
| `M2[2][1]` | `1.3932982808117473` | `0x1.64af3235a8d3ap+0` |
| `M2[2][2]` | `-1.279390402131062` | `-0x1.4786211ff81d8p+0` |
| `transfer` | `"depcubic"` | |
| `depcubic_alpha` | `0.021` | `0x1.5810624dd2f1bp-6` |
| `L_corr_pw[0]` | `-0.040505986477810375` | `-0x1.4bd335e6dc7a5p-5` |
| `L_corr_pw[1]` | `-0.07435430161969014` | `-0x1.308e22dc60360p-4` |
| `L_corr_pw[2]` | `-0.09474033834131092` | `-0x1.840e71ec9b429p-4` |
| `L_corr_pw[3]` | `-0.09425098468772593` | `-0x1.820d520d97e0dp-4` |
| `L_corr_pw[4]` | `-0.08944146679273975` | `-0x1.6e5a2cec7f555p-4` |
| `L_corr_pw[5]` | `-0.0859403510700625` | `-0x1.6002fd540af71p-4` |
| `L_corr_pw[6]` | `-0.08213769491188778` | `-0x1.506f9d96a5695p-4` |
| `L_corr_pw[7]` | `-0.07795647538955015` | `-0x1.3f4f4a04f3e75p-4` |
| `L_corr_pw[8]` | `-0.07609479014797572` | `-0x1.37af2bb14dfb9p-4` |
| `L_corr_pw[9]` | `-0.06979002219081165` | `-0x1.1ddc246e58ae0p-4` |
| `L_corr_pw[10]` | `-0.06645348378711273` | `-0x1.1031873922783p-4` |
| `L_corr_pw[11]` | `-0.06001753295886477` | `-0x1.eba9e3a7d66fbp-5` |
| `L_corr_pw[12]` | `-0.054328991866965846` | `-0x1.bd10276962743p-5` |
| `L_corr_pw[13]` | `-0.045783833354952314` | `-0x1.770fa85e3db84p-5` |
| `L_corr_pw[14]` | `-0.03846197834800118` | `-0x1.3b149d649cac6p-5` |
| `L_corr_pw[15]` | `-0.029723758728833132` | `-0x1.e6fe7ae9e3ddbp-6` |
| `L_corr_pw[16]` | `-0.020668653228716388` | `-0x1.52a29d6adc648p-6` |
| `L_corr_pw[17]` | `-0.010593613802940943` | `-0x1.5b21ac6a403e3p-7` |
| `L_corr_pw[18]` | `-0.003717042652527224` | `-0x1.e7334142c93b1p-9` |
| `L_corr_pw_step` | `0.05` | `0x1.999999999999ap-5` |
| `enrichment.type` | `"L_gated_hue"` | |
| `enrichment.amp` | `0.058` | `0x1.db22d0e560419p-5` |
| `enrichment.center_deg` | `264.5` | `0x1.0880000000000p+8` |
| `enrichment.sigma` | `0.7` | `0x1.6666666666666p-1` |
| `enrichment.L_lo` | `0.37` | `0x1.7ae147ae147aep-2` |
| `enrichment.L_hi` | `1.0` | `0x1.0000000000000p+0` |
| `enrichment_gr` | `1.5240031403836274` | `0x1.862511def2f02p+0` |
| `enrichment_rt` | `8.326672684688674e-15` | `0x1.2c00000000000p-47` |
| `chroma_power` | `0.978` | `0x1.f4bc6a7ef9db2p-1` |

MetricSpace, `params/metric_params.json`, sha256 `a469b93c278e1708ff32aee3308375ea0fccedac75bdeba4d3b85441e279401b`

| Parameter | Decimal | Hexadecimal |
|---|---|---|
| `M1[0][0]` | `0.7212986433113499` | `0x1.714e0e4757c3fp-1` |
| `M1[0][1]` | `0.45344826541531813` | `0x1.d054bdf98c002p-2` |
| `M1[0][2]` | `-0.19288975751942616` | `-0x1.8b09c9023c048p-3` |
| `M1[1][0]` | `-0.788211869495579` | `-0x1.939081939ae57p-1` |
| `M1[1][1]` | `1.795241376757236` | `0x1.cb94f059992d0p+0` |
| `M1[1][2]` | `0.0876172451181785` | `0x1.66e1572592322p-4` |
| `M1[2][0]` | `-0.0917700599912156` | `-0x1.77e3e1e6a0aeap-4` |
| `M1[2][1]` | `0.45765588659459255` | `0x1.d4a3bea6fb8e4p-2` |
| `M1[2][2]` | `1.2922045513917677` | `0x1.4acdeadff8518p+0` |
| `gamma[0]` | `0.47229813098762524` | `0x1.e3a21f0a3710bp-2` |
| `gamma[1]` | `0.5149184096354483` | `0x1.07a362c2fc33bp-1` |
| `gamma[2]` | `0.5113233386366979` | `0x1.05cc2c324096dp-1` |
| `M2[0][0]` | `-0.26355622180094096` | `-0x1.0de1aea52b616p-2` |
| `M2[0][1]` | `0.4168322883703174` | `0x1.aad61559de922p-2` |
| `M2[0][2]` | `0.4926763141656403` | `0x1.f88023c36bc5cp-2` |
| `M2[1][0]` | `1.8897570508777322` | `0x1.e3c71e3ae7bfap+0` |
| `M2[1][1]` | `-3.1212232034205774` | `-0x1.8f843def1a971p+1` |
| `M2[1][2]` | `1.0421666921060384` | `0x1.0acb6fb393709p+0` |
| `M2[2][0]` | `0.3585108617962056` | `0x1.6f1d78aab3b0cp-2` |
| `M2[2][1]` | `1.7694028193790368` | `0x1.c4f7954aaee37p+0` |
| `M2[2][2]` | `-1.4120626067695372` | `-0x1.697cef5bfacbbp+0` |
| `hk_weight` | `0.2676231133101982` | `0x1.120bcb1d48b53p-2` |
| `hk_power` | `0.8934892185255707` | `0x1.c9776b39cab68p-1` |
| `hk_hue_mod` | `0.7173169828841472` | `0x1.6f442becb4a2cp-1` |
| `L_corr_p1` | `0.5385456675962418` | `0x1.13bc41fb74e9bp-1` |
| `L_corr_p2` | `0.12508858146241716` | `0x1.002e7133df76fp-3` |
| `L_corr_p3` | `0.6768950256217603` | `0x1.5a91fc1bbdb7cp-1` |
| `cs_cos1` | `-0.195370576218515` | `-0x1.901e72dbac6fap-3` |
| `cs_sin1` | `0.5330819227283227` | `0x1.10f01d2069f4bp-1` |
| `cs_cos2` | `0.08863325582067766` | `0x1.6b0ab47167adcp-4` |
| `cs_sin2` | `0.9365540137751136` | `0x1.df8401f8340b0p-1` |
| `cs_cos3` | `0.13789738139719568` | `0x1.1a69f13a70a0dp-3` |
| `cs_sin3` | `0.061650260197979936` | `0x1.f909f76ae3092p-5` |
| `lc1` | `-1.5239477450767043` | `-0x1.8621707deada1p+0` |
| `lc2` | `-1.751157310240011` | `-0x1.c04bd871a1d62p+0` |
| `hk_sin1` | `0.6915224124600773` | `0x1.620f39c3ef21dp-1` |
| `hk_cos2` | `0.48647127559605596` | `0x1.f22586ac838bep-2` |
| `hk_sin2` | `0.9853124591201782` | `0x1.f87adfe88666bp-1` |
| `hue_cos1` | `-0.02833024015436984` | `-0x1.d029a3bcdaa2bp-6` |
| `hue_sin1` | `-0.21131429516166544` | `-0x1.b0c58c972c53cp-3` |
| `hue_cos2` | `0.2189784817615645` | `0x1.c077ca4d8c997p-3` |
| `hue_sin2` | `-0.06871898981942523` | `-0x1.1979155e37a86p-4` |
| `hue_cos3` | `0.005506053349515315` | `0x1.68d83f10f273ep-8` |
| `hue_sin3` | `-0.0641329861299175` | `-0x1.06b04f605d75ap-4` |
| `hlc_cos1` | `-0.43576378069144767` | `-0x1.be38dc4b67829p-2` |
| `hlc_sin1` | `1.060094063845983` | `0x1.0f625316e6f8cp+0` |
| `hlc_cos2` | `0.47931193034584496` | `0x1.ead0bf25ac112p-2` |
| `hlc_sin2` | `-0.2622579649434462` | `-0x1.0c8d5a1a30c24p-2` |
| `hl_cos1` | `0.13610794232685908` | `0x1.16bfc2c8285fcp-3` |
| `hl_sin1` | `0.1168702235362288` | `0x1.deb34fbf6dd2ap-4` |
| `hl_cos2` | `-0.01617739641422492` | `-0x1.090ceb2225409p-6` |
| `hl_sin2` | `0.038145638815030566` | `0x1.387d33e6423a6p-5` |
| `cp_cos1` | `-0.09900209889026965` | `-0x1.9583398f814b6p-4` |
| `cp_sin1` | `0.059635520647228726` | `0x1.e888c05b8351ap-5` |
| `cp_cos2` | `-0.013586499967803128` | `-0x1.bd33d283b05b6p-7` |
| `cp_sin2` | `0.2253393118474472` | `0x1.cd7eb2771a7bap-3` |
| `lp_dark` | `-0.029053748937210654` | `-0x1.dc044160bca76p-6` |
| `lp_dark_hcos` | `1.3346761652952872` | `0x1.55ad5650b1a6ap+0` |
| `lp_dark_hsin` | `-0.1698908144723919` | `-0x1.5befb72065d72p-3` |
| `dist_power` | `1.9737081170404969` | `0x1.f944ef6689e69p+0` |
| `dist_wC` | `3.966003089807536` | `0x1.fba5fd3f4776fp+1` |
| `hue_cos4` | `-0.053592461436994296` | `-0x1.b70789a5e3966p-5` |
| `hue_sin4` | `-0.00954137464208059` | `-0x1.38a6da05f7a39p-7` |
| `Lh_cos1` | `-0.4963251525324449` | `-0x1.fc3ca9293c7aep-2` |
| `Lh_sin1` | `-0.09564696283240552` | `-0x1.87c51c153b053p-4` |
| `cs_cos4` | `0.0641970862504494` | `0x1.06f3861f1b538p-4` |
| `cs_sin4` | `-0.027401052793571013` | `-0x1.c0f05867f576dp-6` |
| `dist_nl` | `0.0` | `0x0.0p+0` |
| `dist_sat` | `0.0` | `0x0.0p+0` |
| `dist_compress` | `52.473130649294724` | `0x1.a3c8f8b8cba61p+5` |
| `dist_linear` | `0.0` | `0x0.0p+0` |
| `dist_post_power` | `0.47897301074925214` | `0x1.ea77e6a356ba2p-2` |
| `dist_sl` | `-0.9155125151657894` | `-0x1.d4be0e6f6e795p-1` |
| `dist_sc` | `2.9268353744941558` | `0x1.76a28aa31d395p+1` |
| `dist_sl_hcos1` | `0.0` | `0x0.0p+0` |
| `dist_sl_hsin1` | `0.0` | `0x0.0p+0` |
| `dist_sl_hcos2` | `0.0` | `0x0.0p+0` |
| `dist_sl_hsin2` | `0.0` | `0x0.0p+0` |
| `dist_sc_hcos1` | `0.0` | `0x0.0p+0` |
| `dist_sc_hsin1` | `0.0` | `0x0.0p+0` |
| `dist_sc_hcos2` | `0.0` | `0x0.0p+0` |
| `dist_sc_hsin2` | `0.0` | `0x0.0p+0` |
| `hk_weight_S` | `0.0` | `0x0.0p+0` |
| `hk_power_S` | `0.0` | `0x0.0p+0` |
| `hk_hue_S` | `0.0` | `0x0.0p+0` |
| `lp_dark_S` | `0.0` | `0x0.0p+0` |
| `lp_dark_S2` | `0.0` | `0x0.0p+0` |
| `cs_S_lin` | `0.0` | `0x0.0p+0` |
| `cs_S_quad` | `0.0` | `0x0.0p+0` |
| `lc_S_lin` | `0.0` | `0x0.0p+0` |
| `lc_S_quad` | `0.0` | `0x0.0p+0` |
| `hl_S_lin` | `0.0` | `0x0.0p+0` |
| `L_S_offset` | `0.0` | `0x0.0p+0` |
| `display_phi_deg` | `-28.2` | `-0x1.c333333333333p+4` |

<!-- END GENERATED -->

## 7. Domain

Outputs inside the domain are defined, except near the zero-chroma curve of MetricSpace (both forms), where L depends on the rounding of a0 and b0 (Section 5.5.2); the vectors have no rows there (C0 < 1e-4, Section 9.4). Outside the domain there is no guarantee of meaning. The test vectors record the reference code's behaviour there, at level R, or at level info where that behaviour is not stable (Section 9.4). Implementations MAY return NaN or clip outside the domain. The reference code may add warnings.

| Space | Domain |
|---|---|
| GenSpace (`helmgen`, `helmgenlch`) | 0 ≤ Y ≤ 10, and every component of M1·XYZ is ≥ 0 (M1 of `gen_params.json`) |
| MetricSpace (`helmlab-metric`; also `helmlab-metric-neutral`, Appendix B) | Inside the Display P3 or the Rec.2020 gamut: for P = `lin_P3_to_XYZ` or for P = `lin_2020_to_XYZ` of CSS Color 4, every component of (B·P)⁻¹·XYZ is in [−1e-6, 1 + 1e-6] |
| ΔE `helmlab` | Base coordinates of two colours that are in the MetricSpace domain |

XYZ and Y are `xyz-astm-d65` values, and B is the adaptation from `xyz-css-d65` to `xyz-astm-d65` (Section 9.3). The vectors decide each sign and each bound in exact rational arithmetic. The slack of 1e-6 admits rounding, such as the sRGB white of Section 5.2, whose Y is 1.0000000102999.

Reasons for these bounds (informative):

- GenSpace clips negative LMS values (Section 5.3.5), so XYZ → Lab → XYZ does not return such input. The bound Y ≤ 10 is the range of the statement on the grey axis in Section 3. Inside the domain, on the rows of the vectors, XYZ → Lab → XYZ changes XYZ by up to 5.8e-9, and a change of 1e-12 in one input coordinate moves the inverse by up to 2.1e-11.
- Every sRGB colour of Section 5.2 lies in the Rec.2020 gamut. Some lie up to 1.1e-5 outside the Display P3 gamut, towards blue and yellow, because Section 5.2 uses its own matrix. Display P3 colours near red and cyan lie outside Rec.2020 (22 of the 729 Display P3 inputs of `conv-1`). The union covers all three spaces.
- A bound on Y alone is not enough for MetricSpace. For 98 bright Rec.2020 magentas of `edge-1` with Y ≤ 1 (scaled by 2 and 4 in linear light, so outside both gamuts), a change of 1e-12 in one input coordinate moves the MetricSpace inverse by more than 1e-8, and by up to 1e25. Inside the domain, the largest such change in the vectors is 2.1e-9.

## 8. Conformance

Implementations declare conformance by running the published vectors themselves. There is no certification programme. Conformance is claimed for each identifier on its own. An implementation conforms for an identifier when it passes, at its level in the table below, the rows of that identifier's ops (Section 9.2; for ΔE `helmlab`, the op `deltaE`). Please describe an implementation as "conformant with Helmlab 2026" only together with the identifiers whose rows it passes, for example "conformant with Helmlab 2026 for `helmgen` and `helmgenlch`". `api-1.0-1.json` is not used for conformance.

| Level | Who, on which rows | Tolerance, absolute, per channel |
|---|---|---|
| R | The helmlab Python and JavaScript packages, on the rows at level U and R | Forward coordinates and ΔE ≤ 1e-11; inverses (XYZ) ≤ 1e-8 |
| U | Other implementations, on the rows at level U | Forward coordinates and ΔE ≤ 1e-10; inverses (XYZ) ≤ 1e-8 |
| G (optional) | Any implementation, for each space it provides, on all 16,777,216 colours `#rrggbb` | Each colour goes through the implementation's sRGB to XYZ conversion, the forward, the inverse and the inverse of that sRGB conversion, without gamut mapping. Each channel times 255, rounded to the nearest integer, gives `#rrggbb` again |

- An implementation reads the XYZ inputs, and compares the XYZ outputs, in its own frame: `xyz-astm-d65` or `xyz-css-d65` (Section 9.3). Values of the other frame change the forward outputs by up to 2.0e-4 on `conv-1` and 1.7e-3 on the other sets (measured for an implementation in `xyz-css-d65` given the `xyz-astm-d65` values), so they fail level U.
- LCh outputs are compared as L, C and the hue chord C·|Δh|, where Δh is the difference of the two hues in radians, wrapped to [−π, π], and C is the expected chroma. The chord has the tolerance of the coordinates. A hue in degrees magnifies a difference in a and b by (180/π)/C, so a fixed tolerance on h cannot be met near grey (Section 5.4). An undefined hue (NaN, or `none` in CSS) counts as Δh = π, a chord of C·π.
- Rows at level info are checked at no level (Section 9.4).
- Outputs of package functions outside the definition (Section 2), such as the cusp and the maximum chroma, are informative.
- The Python reference passes level G in `helmgen`, `helmgenlch`, `helmlab-metric` and `helmlab-metric-neutral`. The largest change is 2.6e-4 of an 8-bit step (`#0affff`), the same as that of its sRGB conversion (Section 5.2) alone.

**Other implementations, measured on 2026-11-12.** Each library was run on the rows at level U of every vector file, with the scripts in `tools/spec/ports/`. All three read XYZ in the frame `xyz-css-d65`; in `xyz-astm-d65` they miss the tolerances by up to 1.7e-3 (forward) and 3.4e-3 (inverse), as Section 9.3 leads one to expect.

| Implementation | Checked | Result at level U (largest difference) |
|---|---|---|
| Color.js 0.7.0 and 0.7.1 | `helmgen`, `helmgenlch`, `helmlab-metric`, ΔE `helmlab` | Passes every row: forward 1.1e-12, inverse 4.8e-14, ΔE 2.1e-13 |
| ColorAide 8.9 and 8.14.1 | `helmgen`, `helmgenlch`, `helmlab-metric`, ΔE `helmlab` | The three spaces pass every row: forward 2.8e-13, inverse 6.6e-13. ΔE `helmlab` fails all 2,000 rows, by up to 0.25: ColorAide uses the ΔE parameters of helmlab 0.2.0 to 0.9.0, which helmlab 0.9.1 replaced |
| color-space 3.5.0 and 3.6.0 | `helmgen`, `helmgenlch`, and its `helmlab-metric` on the rows of Appendix B | Passes every row: forward 9.5e-15, inverse 7.0e-13. Its `helmlab-metric` is the neutral-corrected form, so on the base-form rows it differs by up to 0.24 in the forward and by up to 42 in the inverse (XYZ, row `rec2020/8-0-8` of `conv-1`) |

ColorAide returns an undefined hue for greys in `helmgenlch`. On these rows the chord C·π stays below the tolerance.

## 9. Test vectors

The vectors are the files in `vectors/`. A published file never changes and never gains rows; new cases go into a new file, `<set>-<n>.json`. Every file is listed in `SHA256SUMS` and is dedicated to the public domain under CC0 1.0.

`tools/spec/make_vectors.py` produced them in the environment of Appendix A.1; run twice, it gives the same bytes. `tools/spec/check_vectors.py` and `tools/spec/check_vectors.mjs` check the Python and the JavaScript package against them at level R. On 2026-11-13 the Python package (1.0.2) and the npm packages 1.0.0, 1.0.1 and 1.0.2 (Node.js 26.8.1) each passed: 30,742 rows checked, the 1,007 rows at level info skipped, no failures.

### 9.1 Sets

| File | Inputs | Rows | Content |
|---|---|---|---|
| `core-1.json` | 21 | 126 | Black, white, greys at Y = 0.18 and 0.5, the sRGB primaries and secondaries, the Display P3 and Rec.2020 primaries, a negative XYZ, a negative LMS, greys at Y = 2 and 10, sRGB red at twice its luminance |
| `conv-1.json` | 2,187 | 13,122 | 9×9×9 grids of the encoded values k/8 in sRGB, Display P3 and Rec.2020 |
| `grey-1.json` | 1,001 | 5,005 | Greys with Y = 16·(k/1000)², k = 0 to 1000; Y = 1 at k = 250 |
| `edge-1.json` | 1,024 | 5,120 | 256 inputs each: negative LMS (the draws of Section 5.3.5), imaginary colours (a negative X, Y or Z), chromatic colours above white (Y from 1 to 16), bright Rec.2020 magentas (scaled by 1, 2, 4 and 8 in linear light) |
| `de-1.json` | 2,000 pairs | 2,000 | Synthetic pairs of Rec.2020 colours: 1,000 close pairs and 1,000 independent pairs |
| `neutral-1.json` | 3,188 | 6,376 | The neutral-corrected form (Appendix B) on the inputs of `grey-1` and `conv-1` |
| `api-1.0-1.json` | — | — | `reference-1.0.json` of the helmlab package, byte for byte |

Each input of `core-1` and `conv-1` has six rows: the forward and the inverse of `helmgen`, `helmgenlch` and `helmlab-metric`. `grey-1` and `edge-1` have no `helmgenlch.inv` rows. The ΔE pairs are synthetic for licensing reasons. `api-1.0-1.json` is the regression set of the helmlab package, not a definition vector: besides conversions it holds palette, gradient, cusp, contrast and token outputs. It has no header.

The inputs of `core-1` and their forward outputs. Each number is the shortest decimal that reads back as the same binary64 value.

| id | xyz-astm-d65 | helmgen (L, a, b) | helmlab-metric (L, a, b) |
|---|---|---|---|
| core/black | 0.0, 0.0, 0.0 | 0.0, 0.0, 0.0 | 0.0, 0.0, 0.0 |
| core/white | 0.95047, 1.0, 1.08883 | 0.9999999999999998, 2.37205375058874e-16, -4.74410750117748e-16 | 1.1211492319729706, 0.10163688399210657, 0.20041318295782515 |
| core/grey-0.18 | 0.1710846, 0.18, 0.19598939999999998 | 0.49045285137808425, -5.926275944761233e-16, -4.741020755808986e-16 | 0.7345175951675056, 0.24066071488702862, 0.15432534279871837 |
| core/grey-0.5 | 0.475235, 0.5, 0.544415 | 0.7590152487023349, 0.0, -9.47567870624053e-16 | 0.9748991198615568, 0.17398430780100327, 0.20596578731385257 |
| core/srgb-red | 0.4124564390896922, 0.2126728514056226, 0.0193338955823293 | 0.5632406520932904, 0.29471787551351664, 0.18550836212223776 | 0.9207283173164369, 0.940821340214519, -0.20629901228900596 |
| core/srgb-green | 0.3575760776439511, 0.7151521552878178, 0.1191920258813418 | 0.8463899596030668, -0.29349611360629385, 0.3038041098305225 | 0.7405747931790364, -0.1547720191626349, 0.8916585175791307 |
| core/srgb-blue | 0.1804374832663989, 0.0721750036064596, 0.9503040785363679 | 0.3648928088510713, -0.04191420507221431, -0.493523357523643 | 0.7623162469855196, -0.12541914692709508, -0.20263496227418212 |
| core/srgb-cyan | 0.53801356091035, 0.7873271588942774, 1.0694961044177096 | 0.8941173762514437, -0.19477396023688262, -0.05192276675649077 | 0.8824553686532568, -0.16707583683924032, 0.22456411141685242 |
| core/srgb-magenta | 0.5928939223560912, 0.2848478550120822, 0.9696379741186972 | 0.6504478096511678, 0.3527805572152495, -0.2825638949117671 | 1.3271843182860856, 0.094577409508457, -0.07810562793919387 |
| core/srgb-yellow | 0.7700325167336433, 0.9278250066934404, 0.1385259214636711 | 0.964838232885624, -0.08491388059686626, 0.32222641194591867 | 0.8997792183191262, 0.17785131197169998, 0.8759352384154622 |
| core/p3-red | 0.4865966042738641, 0.2289855262694398, -9.161322772031926e-07 | 0.5866938520827683, 0.34191892390080963, 0.21278811951194584 | 0.9724467568425391, 0.9573782826272872, -0.2681782928893235 |
| core/p3-green | 0.26569147550509964, 0.6917385161287015, 0.045110523202416826 | 0.8249345452959719, -0.38109903363728953, 0.3533615789478416 | 0.7051298217148839, -0.2387196582426005, 1.014335461008113 |
| core/p3-blue | 0.19818192022103628, 0.07927595760185867, 1.0437203929298604 | 0.3810098099292466, -0.04275154280044526, -0.5083473623628557 | 0.7804592198309552, -0.12013298407705753, -0.19447008743651217 |
| core/rec2020-red | 0.6369908659234614, 0.26271491004018815, -1.6762449594429903e-06 | 0.6323106732269393, 0.43067666146862993, 0.21295441502749626 | 1.3779512862778835, 0.3755779768463135, -0.1673351936222981 |
| core/rec2020-green | 0.14463594194092258, 0.677995186799333, 0.02807426584673224 | 0.802024269650607, -0.5212279501830006, 0.3725653050469163 | 0.7208279757639227, -0.3658127972188819, 1.0565232250305165 |
| core/rec2020-blue | 0.168843192135616, 0.05928990316047886, 1.0607574103982271 | 0.3353463190269114, -0.19794776307475045, -0.5384623286966809 | 0.8122848975865923, -0.1510209859849525, -0.2036955702133576 |
| core/negative-xyz | -0.05, -0.05, -0.05 | 0.0, 0.0, 0.0 | -0.3157585330879478, -0.2737183531818679, -0.2011960902300976 |
| core/negative-lms | -0.03721012230416233, 0.08556165376589255, 0.36374929771311004 | 0.26258301699364106, -0.9769899585398601, -0.18714918695925908 | 0.8571434638446472, -0.798980535518446, 0.0022195400546781707 |
| core/grey-2 | 1.90094, 2.0, 2.17766 | 1.2632074893724938, -7.086364231091119e-16, -1.4172728462182238e-15 | 1.3414218614893383, 0.03337649194059151, 0.13456668508162287 |
| core/grey-10 | 9.5047, 10.0, 10.8883 | 2.1663497793510853, -2.808797513891267e-15, -1.8725316759275113e-15 | 19.03829669317455, -0.00033810440984451823, 0.0019131891231165387 |
| core/srgb-red-x2 | 0.8249128781793844, 0.4253457028112452, 0.0386677911646586 | 0.7556349162876022, 0.3677325051600349, 0.23031571478524612 | 1.1147213954308688, 0.8484840459712474, -0.18314532189728533 |

Rows of `core-1` that are not at level U (all outside the domain):

- `helmgen` and `helmgenlch`, forward and inverse, level R: `core/negative-xyz`, `core/negative-lms`.
- `helmlab-metric.fwd`, level R: `core/negative-xyz`, `core/negative-lms`, `core/grey-2`, `core/grey-10`, `core/srgb-red-x2`.
- `helmlab-metric.inv`, level R: `core/negative-xyz`, `core/negative-lms`, `core/grey-2`, `core/srgb-red-x2`; level info: `core/grey-10`.

### 9.2 Operations

| `op` | Input | Output | Reference (Python, one call of shape (1, 3)) |
|---|---|---|---|
| `helmgen.fwd`, `helmgen.inv` | XYZ; (L, a, b) | (L, a, b); XYZ | `GenSpace().from_XYZ`, `GenSpace().to_XYZ` |
| `helmgenlch.fwd`, `helmgenlch.inv` | XYZ; (L, C, h) | (L, C, h); XYZ | `hl.gen.to_lch` of the GenSpace forward; `hl.gen.from_lch`, then `GenSpace().to_XYZ` |
| `helmlab-metric.fwd`, `helmlab-metric.inv` | XYZ; (L, a, b) | (L, a, b); XYZ | `MetricSpace().from_XYZ`, `MetricSpace().to_XYZ` |
| `helmlab-metric-neutral.fwd`, `.inv` | XYZ; (L, a, b) | (L, a, b); XYZ | `MetricSpace(neutral_correction=True)` |
| `deltaE` | two base-form triples | ΔE | `MetricSpace().distance_from_lab` |

The input of an inverse row is the expected output of the forward row with the same `id`.

Each file is one JSON object: the header, then `rows`, one row per line. A row has the keys `id`, `op`, `in`, `out`, `domain` and `level`. A forward row also has `src` when its XYZ was made from another value, such as an encoded colour or a grey level. A forward row gives the XYZ input in both frames and the coordinates as output. An inverse row gives the coordinates as input and the XYZ output in both frames. A ΔE row gives `lab1` and `lab2`, and the XYZ of both colours in both frames, as input, and ΔE as output. Numbers are the shortest decimals that read back as the same binary64 values. The format writes a non-finite value as the string `"NaN"`, `"Infinity"` or `"-Infinity"`, and only in rows at level info; these files contain none.

### 9.3 Inputs and frames

Each XYZ is given in both frames. The `xyz-astm-d65` value is what the reference code receives. The `xyz-css-d65` value is A·`xyz-astm-d65`, computed exactly and rounded once per component, where A is the Bradford adaptation from the white of `xyz-astm-d65`, W_astm = (0.95047, 1, 1.08883), to the CSS Color 4 white, W_css = (3127/3290, 1, 3583/3290):

```text
A = M_BFD^-1 · diag((M_BFD · W_css) / (M_BFD · W_astm)) · M_BFD

M_BFD = [[ 0.8951,  0.2664, -0.1614],
         [-0.7502,  1.7135,  0.0367],
         [ 0.0389, -0.0685,  1.0296]]
```

A·W_astm = W_css exactly, and B = A⁻¹ maps `xyz-css-d65` to `xyz-astm-d65`. Rounded to binary64:

```text
A = [[ 0.9999570254019493,     -2.0718742727340958e-05,  4.36173329246224e-05  ],
     [-2.694517763356012e-05,   1.000009385946498,       1.4900982235455368e-05],
     [ 7.94345929293135e-06,   -1.2878245967317602e-05,  1.000214063706999     ]]
```

With this rounded A and left-to-right sums in binary64, the `xyz-css-d65` values of the rows at level U and R are reproduced to within 3.6e-15.

Color.js 0.7 and color-space work in `xyz-css-d65` and apply the inverse of this step before M1. Given the `xyz-css-d65` values, they compute the function of Section 5 on the `xyz-astm-d65` values, to rounding. An emulation that applies B in binary64 before the reference code, and A after its inverse, differed from the expected outputs by up to 5.1e-15 on the rows at level U.

There is one set of expected outputs per row, computed by the Python reference from the `xyz-astm-d65` value with one call of shape (1, 3) (Appendix A.2). Inverse rows give the expected XYZ in both frames, the second as A times the first.

Where the inputs come from (informative; the `src` field of a row and the `groups` of the header):

- sRGB values: Section 5.2.
- Display P3 and Rec.2020 values: the conversions of CSS Color 4 (`conversions.js` of the Editor's Draft, 2026-11-09; Rec.2020 decodes as v^2.4), then B. The white (1, 1, 1) of each space then gives exactly W_astm. The reference code's own conversions for `color(display-p3 …)` and `color(rec2020 …)` are not used (Section 5.2).
- Greys: Y·W_astm in binary64.
- Drawn XYZ: the seeds and ranges are in each file header.

### 9.4 Labels

- `domain`: `in` or `out`, by the rules of Section 7. An inverse row has the flag of its forward row.
- MetricSpace inputs keep C0 ≥ 1e-4, away from the zero-chroma curve (Section 5.5.2). An input with C0 < 1e-4 gets no MetricSpace rows, except black, whose output is exactly (0, 0, 0). The header lists the excluded inputs. These files exclude none: the smallest C0 of an input other than black is 3.9e-3, at `grey/1`.
- `level`:
  - `U` for every row inside the domain. Every implementation is checked on it, at the tolerance of its level (Section 8).
  - `R` for rows outside the domain. Only the helmlab packages are checked; other implementations MAY return NaN or clip (Section 7). The values are reference behaviour, with no guarantee of meaning.
  - `info` for rows outside the domain where the reference behaviour is not stable: (a) an output is not finite; (b) one call with all inputs of the file for the op, of shape (N, 3), differs from the one-row call by more than the level R tolerance; (c) for an inverse row, the forward of the expected XYZ does not return the input to within 1e-8; or (d) for an inverse row, adding or subtracting 1e-12 to one input channel moves the output by more than 1e-8. No level checks these rows.
- The generator stops if a row inside the domain meets one of the conditions (a) to (d); none does. In these files the 1,007 rows at level info are all MetricSpace inverses, and each meets (c) or (d): the greys from Y ≈ 3.9 to 10.1 (`grey/494` to `grey/795`, in both forms), `core/grey-10`, and 402 rows of `edge-1`.
- The header gives the reference package and its git tree, the environment (Appendix A.1), the sha256 of the generator and of the parameter files, the frames and the matrix A, the tolerances of Section 8, the domain and level rules, how each group of inputs was made, the excluded inputs and the number of rows per op, level and domain.

### 9.5 Checking an implementation

1. Choose the frame of the implementation's XYZ. Read the inputs of forward rows, and the expected outputs of inverse rows, in that frame. ΔE rows give the base coordinates `lab1` and `lab2` as input.
2. Check every `U` row; the helmlab packages also check every `R` row. Skip `info` rows.
3. Compare each channel with the tolerance of the level (Section 8). For `helmgenlch`, compare L, C and the hue chord C·|Δh|; an undefined hue counts as Δh = π.
4. Level G (optional): run the round trip of Section 8 on all 16,777,216 colours, for each space the implementation provides.

## 10. Known limitations

| Topic | Limitation | Source |
|---|---|---|
| GenSpace lightness | A piecewise-linear lightness correction with 19 knots; its slope jumps by up to 0.42 (at L = 0.15), and the darkest 5% of L is compressed about five times (slope 0.19) | Section 5.3.2 (knot table) |
| GenSpace negative LMS | Negative LMS values are clipped. For such input, XYZ → Lab → XYZ changes XYZ by up to 0.227 (566 of 20,000 random XYZ in [−0.05, 1.2]³) | Section 5.3.5 |
| GenSpace inverse | Iterative (Halley's method). Near grey the round trip changes the colour by up to 2.5e-5 × Y; very dark colours (Y below about 1e-35) are moved towards grey | Sections 5.3.3 and 5.3.5 |
| GenSpace in sRGB | Small separate in-gamut regions near h ≈ 260°, from the roots of the cubic in the inverse (Step GI6) | helmlab 1.0.2 README |
| MetricSpace above white | L grows quickly: about 1.80 at Y = 3 and 19 at Y = 10. The inverse returns greys to a relative 1e-12 up to Y ≈ 17 and fails above Y ≈ 20; the neutral-corrected form only up to Y ≈ 4.3 | helmlab 1.0.2 README (L); Sections 5.5.3 and B.2 (round trip) |
| MetricSpace inverse outside the domain | Can be ill-conditioned. A change of 1e-12 in one input coordinate moved XYZ by up to 1e25 for bright Rec.2020 magentas with Y ≤ 1, and by up to 1.9 for greys from Y ≈ 3.9 to 10.1 (both forms). The vectors mark such rows `info` | Sections 7 and 9.4 |
| MetricSpace zero-chroma curve | L is not continuous on a curve of greyish blues (y ≈ 0.276) where the internal chroma is 0. Inside the sRGB gamut (Y from 0.01 to 0.62) the limits of L differ by up to about 0.2, and the reference returns different L for different call shapes | Section 5.5.2 |
| MetricSpace near black | The lightness axis is steep near black. ΔE `helmlab` is 0.0734 for `#010101` against `#000000` and 0.0290 for `#ffffff` against `#f7f7f7`; CIEDE2000 gives 0.157 and 1.60. Of the 3,989 colour pairs of the fitting data (Section 5.6), 7 have a colour darker than L* 10 and none has a colour lighter than L* 92.8, so below L* ≈ 10 ΔE `helmlab` rests on 7 pairs and above L* ≈ 93 on none | [helmlab#5](https://github.com/helmlab-space/helmlab/issues/5); Section 5.6 |
| Fitting data | GenSpace was selected with CMA-ES on ColorBench v1 win counts. MetricSpace was fitted to COMBVD, MacAdam 1974 and Helmlab's Human Feedback data, and He et al. 2022 was used in model selection (Section 5.6). Figures on these datasets are in-sample. | helmlab 1.0.2 README |
| Neutral-corrected form (Appendix B) | Greys are kept neutral only up to Y = 2 (grey chroma at most 8.0e-4). Above Y = 2 the grey chroma is 0.073 at Y = 2.59, 0.11 at Y = 3, about 0.14 for Y from 4 to 10, above 0.9 from Y ≈ 10.92 and 3.8 at Y = 100 | Appendix B; [helmlab#4](https://github.com/helmlab-space/helmlab/issues/4) |
| `jnd()` | Divides ΔE by 0.0356, from a fixed confidence model (47 pairs, one session). Not part of this definition. | Reference code |

## 11. Errata, changes and identifier history

### Errata

The errata are listed in `errata.yaml`.

| Identifier | Summary | Class | Status |
|---|---|---|---|
| E-2026-0001 | The paper describes the neutral correction (Stage 10) with 384 grey levels for Y from 0.001 to 2.59. The reference code uses 255 levels and holds the correction constant above L ≈ 1.341. The code governs. | K2 | Documented; outputs unchanged |
| E-2026-0002 | The paper says MetricSpace was trained on COMBVD and cross-validated on He et al. 2022 and MacAdam 1974, and calls MacAdam 1974 held-out. The parameters were fitted to COMBVD, MacAdam 1974 and Helmlab's Human Feedback data (Section 5.6), so results on all three are in-sample. | K2 | Documented; outputs unchanged |

### Change classes

| Class | Meaning | In a frozen definition |
|---|---|---|
| K1 Editorial | Wording, links, layout | Allowed |
| K2 Explanation | Documents existing behaviour | Allowed, with an erratum |
| K3 Clarification | Vector outputs unchanged; a conformance status becomes clear, or the reference code changes within the level R tolerance | Allowed after 14 days' notice; new Spec Rev and DOI |
| K4 Output change | Any vector output changes beyond tolerance | Not allowed; goes to the next definition |
| K5 Registration | New vector file, conformance row, reserved identifier, or a single ΔE replacing the previous one | Allowed; new Spec Rev and DOI |

Classes K1 and K2 need a commit only. A Spec Rev receives at most two rounds of K1 and K2 corrections; later corrections go into a new Spec Rev.

### Identifier history

| Implementation | Form | Conformance |
|---|---|---|
| Color.js 0.7.x ("experimental", registered by default) | Base | Level U, including ΔE (0.7.0 and 0.7.1, Section 8) |
| ColorAide 8.9 and later (not registered by default) | Base | Spaces at level U (8.9 and 8.14.1, Section 8); ΔE `helmlab` does not conform (parameters of helmlab 0.2.0 to 0.9.0) |
| helmlab: `hl.metric`; postcss-helmlab `helmlab()` and `helmlch()` | Neutral-corrected | Conforms to Appendix B |
| color-space 3.5.0 and 3.6.0, under the name `helmlab-metric` | Neutral-corrected | Appendix B at level U (Section 8) |

In drafts the reserved identifier of Appendix B was `helmlab-metric-display`. It was renamed `helmlab-metric-neutral` before publication, because in color-space the suffix `-display` marks display-encoded spaces ([helmlab#4](https://github.com/helmlab-space/helmlab/issues/4)).

Defaults in the reference code: `hl.metric` (Python and JavaScript) applies the neutral correction. The Python class `MetricSpace` does not (`src/helmlab/spaces/metric.py`, line 466); the JavaScript core class does (`packages/helmlab-js/src/core/analytical.ts`, line 32).

A rotation other than −28.2° (`abRotateDeg`) and the surround parameter give non-normative variants. The vectors cover only the default variant, which their header names (`variant`).

## 12. Notes for library maintainers

- Registering these spaces by default is not required, and leaving ΔE `helmlab` optional does not affect conformance.
- Use `helmlab-metric` only for the base form. The form with the neutral correction is `helmlab-metric-neutral` (Appendix B).
- A proposal for a later definition will include a way to remove the Helmlab 2026 identifiers, so that the footprint in a library does not grow.
- No candidate definition is proposed to any library before it passes a blind test.
- Bundle size matters to downstream users: in WordPress Gutenberg [PR #80272](https://github.com/WordPress/gutenberg/pull/80272), `build/scripts/theme/index.min.js` grew by 3.54 kB gzip (+15.67%), and a reviewer attributed the increase to Helmlab.

## 13. Licence, citation and evidence

- This document is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Section 3, the parameter files and the test vectors are dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). The reference code is under the MIT License.
- Patents: the editor has not applied for and will not apply for patents on the methods in this document. No employer or institution holds rights in this work.
- Citation: Görkem Yıldız, *Helmlab 2026*, Spec Rev 1. DOI: assigned on publication. Paper: [arXiv:2602.23010](https://arxiv.org/abs/2602.23010).
- This document contains no performance figures. They will be published on the evidence page of helmlab.space, each with its dataset, its size (n) and whether that data was used to fit, tune or design the model. The link will be added when the page is published (class K1).

## Appendix A. Reference environment

### A.1 Environment

The rows of the vectors of Section 9 were produced on 2026-11-09; the file headers were last regenerated on 2026-11-13, with the same rows. Environment:

| | |
|---|---|
| Reference code | `src/helmlab` of the helmlab repository, git tree `0bb705b9571ccdce897dabbec9199930ae0f5a3a`: the files of tag `v1.0.0` (commit `96a0f16`), with the version string "1.0.2" in `__init__.py` and changed docstrings in `helmlab.py` (the same syntax tree once docstrings are removed) |
| Python | 3.14.7 (Clang 21.0.0) |
| numpy | 2.4.4, with BLAS and LAPACK from Apple Accelerate |
| scipy | 1.17.1 |
| System | macOS 26.6.2 on arm64, with its C library |
| Generator | `tools/spec/make_vectors.py`; its sha256 is in each file header |

There is no lock file. Each vector file repeats these values in its header (`reference`, `environment`, `generator`). Another environment can give different last bits (A.2), and so files with different bytes; the vectors are the files listed in `SHA256SUMS`.

### A.2 Calling convention

The reference code can return different last bits for the same colour, depending on how it is called:

- Matrix products go through BLAS, whose summation order can depend on the number of rows in the call. On the editor's machine (macOS arm64, numpy 2.4.4 with Accelerate), a call with one row of shape (1, 3) gives the left-to-right sums of Section 5.1 bit for bit. An implementation of Section 5 with left-to-right sums matched one-row calls in every row of these sets: the sRGB 17-level grid, 1,001 greys from Y = 0 to 1 and 20,000 random XYZ in [0, 1]³, forward and inverse, in GenSpace and in the MetricSpace base form (103,656 rows). Batch calls differed from one-row calls in the last bits in most rows.
- On numpy arrays `x**2` is `x*x`. A call with one colour of shape (3,) works on numpy scalars, where `x**2` is `pow(x, 2)`; on the editor's machine `pow(x, 2)` differs from `x*x` in about 0.13% of arguments. Such calls differed from one-row calls in up to 0.16% of the rows of these sets, in the last bit.
- Near the zero-chroma curve of MetricSpace the differences exceed the last bits (Section 5.5.2).
- In GenSpace, an input of shape (N, 1, 3) gives an output of shape (N, N, 3): the mean of GF3 and GI5 is broadcast across the colours of the call.

The vectors are computed with one call of shape (1, 3) per input.

### A.3 Differences between the Python and JavaScript packages

Both packages are reference code at level R (Section 8). They compute in a different order in places, which changes the last bits:

- Matrix products: BLAS in Python, left-to-right sums in JavaScript. Matrix inverses: `numpy.linalg.inv` (LU decomposition) in Python, Cramer's rule in JavaScript.
- Elementary functions (`exp`, `pow`, `sin`, `cos`, `atan2`, `sinh`, `asinh`): the platform C library through numpy in Python, the engine's `Math` functions in JavaScript. ECMA-262 does not fix their accuracy.
- Several expressions are associated differently, and some steps are ordered differently where the order does not matter in exact arithmetic.
- JavaScript 1.0.1 and 1.0.2 stop their Newton and Halley loops when an update leaves the value unchanged. This gives the same result as the fixed counts, which JavaScript 1.0.0 runs.
- The JavaScript neutral-corrected form reads a correction table of 258 entries with values rounded to 12 decimal places (Appendix B.1).

Largest differences between the npm packages and the expected outputs of the vectors (Section 9), per channel, measured with Node.js 26.8.1. The packages 1.0.0, 1.0.1 and 1.0.2 give the same maxima.

| Operation | Rows at level U | Rows at level R |
|---|---|---|
| GenSpace, forward | 2.6e-15 | 2.6e-15 |
| GenSpace, inverse | 2.5e-14 | 3.2e-14 |
| GenSpace LCh, forward (L, C and the chord) | 3.2e-15 | 3.4e-15 |
| GenSpace LCh, inverse | 7.1e-15 | 5.0e-16 |
| MetricSpace base form, forward | 1.6e-15 | 7.8e-14 |
| MetricSpace base form, inverse | 4.9e-14 | 6.9e-14 |
| Neutral-corrected form, forward | 7.8e-13 | 7.2e-13 |
| Neutral-corrected form, inverse | 1.1e-9 | 4.0e-9 |
| ΔE | 0 | — |

- The neutral-corrected form differs more because its correction table is rounded to 12 decimal places (B.1). In the inverse the difference grows where the inverse is sensitive to its input: inside the domain, a change of 1e-12 in one input coordinate moves the output by up to 2.1e-9 (Section 7).
- ΔE: the JavaScript route for base coordinates is `distanceFromLab` on `new MetricSpace(p, { neutralCorrection: false, abRotateDeg: 0 })`, which evaluates Steps D1 to D7 on its input. The default base-form object first undoes the rotation of Step F14; its results differ by up to 1.7e-16.

### A.4 How Section 5 was checked

On 2026-11-09, one implementation of GenSpace and one of MetricSpace were written from that day's draft of Section 5 alone, which then also held what is now Appendix B (`tools/spec/cleanroom/README.md` maps the numbering), without access to the reference code, and compared with the reference code (Python, `src/helmlab` with the files of tag `v1.0.0` and the version string "1.0.1", batch calls). Largest absolute differences:

| Set | Rows | GenSpace | MetricSpace base form |
|---|---|---|---|
| sRGB 17-level grid, forward | 4,913 | 8.9e-16 | 1.6e-15 |
| Greys, forward | 1,001 | 1.6e-15 (Y from 0 to 10) | 1.3e-15 (Y from 0 to 1) |
| Random XYZ, forward | 20,000 | 3.6e-15 (XYZ in [−0.05, 1.2]³) | 2.1e-14 (XYZ in [0, 1]³) |
| Grid, inverse | 4,913 | 6.7e-16 | 6.7e-16 |
| Random XYZ, inverse | 20,000 | 1.6e-15 | 1.2e-15 |

- LCh (Section 5.4): Lab → LCh → Lab bit-identical on the grid; h within 1.1e-12 degrees where C ≥ 1e-4.
- ΔE (Section 5.6): bit-identical on 5,000 pairs of grid colours given the same coordinates; 7.1e-16 end to end.
- Neutral-corrected form (Appendix B): forward 1.9e-15 on the grid and 6.0e-15 on 601 greys from Y = 0 to 3; inverse 1.9e-14 on the grid. The interpolant of Step N6 was bit-identical to scipy 1.17.1 `PchipInterpolator` on 200,381 values.
- With matrix products computed as numpy computes them for the whole batch, the GenSpace implementation matched the reference bit for bit in every row.

The two implementations and their test scripts are in `tools/spec/cleanroom/`.

## Appendix B. Neutral-corrected form: `helmlab-metric-neutral`

The neutral-corrected form is MetricSpace with the neutral correction applied. Its identifier is reserved as `helmlab-metric-neutral`. Its outputs are those of helmlab 1.0.0 and change only as Section 2 allows; its vectors are the `neutral` set.

- `hl.metric` in the helmlab package, `helmlab()` and `helmlch()` in postcss-helmlab, and `helmlab-metric` in color-space 3.5.0 and 3.6.0 give this form.
- ΔE `helmlab` is computed from base coordinates, not from neutral-corrected coordinates (B.3).
- The correction table can be built with Steps N1 to N6 or read from the released file `neutral-lut.json` (B.1). The two differ by less than 1e-12 in each component of the correction.
- The correction keeps greys neutral only up to Y = 2, the grey of the last knot x_254. For Y ≤ 2 the grey chroma is at most 8.0e-4 (at Y ≈ 9e-5); it is at most 3.3e-5 for Y from 0.01 to 1 and 3.6e-7 for Y from 1 to 1.9916. In the last knot interval, Y from 1.9916 to 2, it reaches 1.5e-4 near Y = 1.997: Step N5 makes the derivative at x_254 zero while the grey values (a5, b5) still change with L there. Above Y = 2 the correction is constant while the grey values keep changing: the grey chroma is 0.073 at Y = 2.59, 0.11 at Y = 3, about 0.14 for Y from 4 to 10, above 0.9 from Y ≈ 10.92 and 3.8 at Y = 100 (Section 10).

### B.1 Neutral correction

The neutral-corrected form `helmlab-metric-neutral` is the base form with one more stage, stage 10 in `metric.py`. Between Step F12 and Step F14 (Section 5.5.2) it subtracts from (a5, b5) an offset that depends only on L3. The offset is interpolated from the values (a5, b5) of greys, before the rotation. The reference code applies it in Python `MetricSpace(neutral_correction=True)` and in `Helmlab().metric` (the default of the `Helmlab` constructor is `neutral_correction=True`). In JavaScript it is on by default, in `new AnalyticalSpace(params)` and in `new Helmlab().metric`. Both `Helmlab` classes pass the rotation angle −28.2 explicitly; it equals `display_phi_deg`.

#### Correction table (Python)

The Python package builds the table on first use (`_build_neutral_lut`, default `n_samples=256`).

**Step N1. Grey levels.** 257 values of Y, in this order:

```text
Y = [0] + linspace(0.001, 0.01, 10) + linspace(0.01, 0.1, 20) + linspace(0.1, 2.0, 226)

linspace(s, e, n)[j]     = j*((e - s)/(n - 1)) + s        for j = 0, ..., n - 2
linspace(s, e, n)[n - 1] = e
```

226 is `n_samples − 30`. The values 0.01 and 0.1 occur twice.

**Step N2. Grey colours.**

```text
XYZ_j = (Y_j*0.95047, Y_j*1.0, Y_j*1.08883)
```

**Step N3. Grey coordinates.** Each XYZ_j goes through Steps F1 to F12 of the base form. Neither the neutral correction nor the rotation of Step F14 is applied. The results are (L3_j, a5_j, b5_j). The Python reference evaluates the 257 greys in one call of shape (257, 3), so the triples carry the rounding of a batch call (Appendix A.2): they differ from the left-to-right sums of Section 5.1 by up to 2.2e-15 per component. An implementation that builds the table with those sums matches the reference only to within such differences; the tolerances of Section 8 allow for them.

**Step N4. Sorting and duplicates.** The 257 triples are sorted by L3 in ascending order. Then, going through the sorted list, the triple at position i ≥ 1 is kept only if

```text
L3_sorted[i] - L3_sorted[i - 1] > 1e-12
```

where i − 1 is the position just before it in the sorted list. This removes one of the two identical triples for Y = 0.01 and one for Y = 0.1 (which one does not matter: they are bit-identical). 255 knots remain, x_0 to x_254, with values ya_k (from a5) and yb_k (from b5). The first knot is (0, 0, 0), from Y = 0. The last is x_254 = 1.3414218614893383 with ya_254 = −0.0341747723961513 and yb_254 = 0.13436617089054925, from Y = 2.

**Step N5. Flat extension.** 125 knots are appended at `numpy.arange(x_254 + 0.01, 2.6, 0.01)`, that is from x_254 + 0.01 in steps of 0.01 while below 2.6. The last is x_379 = 2.5914218614893394. Each appended knot has the values ya_254 and yb_254. The table now has 380 knots. Because the appended values are all equal, the correction is constant for L3 ≥ x_254. The positions of the appended knots do not change the result. Their presence does, because it makes the derivative at x_254 equal to 0 (Step N6).

**Step N6. Interpolant.** `scipy.interpolate.PchipInterpolator` (with `extrapolate=False`) is built separately for ya and for yb over the knots x. For knots x_0 … x_(n−1) and values y_0 … y_(n−1) (n = 380), the derivatives d_k are:

```text
hk_k = x_(k+1) - x_k                                       for k = 0, ..., n - 2
m_k  = (y_(k+1) - y_k) / hk_k

interior knots, k = 1, ..., n - 2:
    w1  = 2*hk_k + hk_(k-1)
    w2  = hk_k + 2*hk_(k-1)
    d_k = 0                                        if m_(k-1) = 0 or m_k = 0 or sign(m_(k-1)) ≠ sign(m_k)
    d_k = 1 / ((w1/m_(k-1) + w2/m_k) / (w1 + w2))  otherwise

end knots:
    d_0     = edge(hk_0, hk_1, m_0, m_1)
    d_(n-1) = edge(hk_(n-2), hk_(n-3), m_(n-2), m_(n-3))

edge(h0, h1, m0, m1):
    e = ((2*h0 + h1)*m0 - h0*m1) / (h0 + h1)
    if sign(e) ≠ sign(m0):                           return 0
    if sign(m0) ≠ sign(m1) and abs(e) > 3*abs(m0):   return 3*m0
    return e
```

On the interval [x_k, x_(k+1)] the interpolant is a cubic in s = L − x_k:

```text
r_k  = (d_k + d_(k+1) - 2*m_k) / hk_k
c0_k = r_k / hk_k
c1_k = (m_k - d_k) / hk_k - r_k
c2_k = d_k
c3_k = y_k
P(L) = ((c3_k + c2_k*s) + c1_k*(s*s)) + c0_k*((s*s)*s)
```

The interval is the k with x_k ≤ L < x_(k+1). For L = x_(n−1), k = n − 2. Pa and Pb denote the interpolants of ya and yb.

**Step N7. Correction in the forward.** After Step F12 and before Step F14:

```text
Lc    = clip(L3, x_0, x_379)            x_0 = 0, x_379 = 2.5914218614893394
a5    = a5 - Pa(Lc)
b5    = b5 - Pb(Lc)
```

Step F14 then rotates the corrected (a5, b5). The clip keeps L inside the table. For L3 ≤ 0 the correction is (0, 0). For L3 ≥ x_254 it is the constant (ya_254, yb_254). A NaN L3 gives NaN.

**Step N8. Correction in the inverse.** In Step I2, after Step I1, the same correction is added back. It is evaluated at the input lightness L, which equals L3 because the correction does not change L:

```text
Lc = clip(L, x_0, x_379)
a5 = a5 + Pa(Lc)
b5 = b5 + Pb(Lc)
```

#### Correction table (JavaScript)

The JavaScript package stores the table in `packages/helmlab-js/src/data/neutral-lut.json`. The file holds three arrays, `L`, `a_err` and `b_err`, with 258 entries each:

- Entries 0 to 254 are the knots x_0 … x_254 and the values ya_k, yb_k of Step N4, each rounded to 12 decimal places.
- Entries 255 and 256 have L = 1.351421861489 and 1.361421861489, that is x_254 + 0.01 and x_254 + 0.02 rounded. Entry 257 has L = 2.59. These three entries have `a_err` = −0.034174772396 and `b_err` = 0.134366170891, the values of entry 254.

`core/neutral-lut.ts` computes the derivatives with the rules of Step N6 when the module loads. It writes the interior derivative as `(w1 + w2) / (w1/m_(k-1) + w2/m_k)`. It evaluates the cubic in Hermite form, with h = x_(k+1) − x_k, t = (L − x_k)/h, t2 = t*t and t3 = t2*t:

```text
P(L) = y_k*(2*t3 - 3*t2 + 1) + h*d_k*(t3 - 2*t2 + t) + y_(k+1)*(-2*t3 + 3*t2) + h*d_(k+1)*(t3 - t2)
```

For L ≤ 0 it returns entry 0. For L ≥ 2.59 it returns entry 257. Both table forms give (0, 0) for L ≤ 0 and a constant correction from the last grey knot upward. They differ only through the 12-decimal rounding of the stored values and the order of arithmetic, by less than 1e-12 in each component of the correction.

### B.2 Round trip

The neutral-corrected inverse returns greys to a relative 1e-12 only up to Y ≈ 4.3, and to 1e-6 up to Y ≈ 5.1; the base form does so to 1e-12 up to Y ≈ 17.4 and to 1e-6 up to Y ≈ 20.3 (Section 5.5.3). On a grid of step 0.01 (relative difference in each XYZ channel) the first failures are at Y = 4.31 and 5.15 in Python and at 4.28 and 5.17 in JavaScript 1.0.0, 1.0.1 and 1.0.2; in the base form they are at 17.40 and 20.31 in both packages.

On the vectors, the JavaScript neutral-corrected inverse differs from the expected outputs by up to 1.1e-9 on rows at level U (`rec2020/8-0-8`) and 4.0e-9 on rows at level R (`grey/796`). Level R allows 1e-8 for inverses (Section 8). For the greys from Y ≈ 3.9 to 10.1, a change of 1e-12 in the input moves the neutral-corrected inverse by more than 1e-8, and by up to 1.9. There the JavaScript and Python neutral-corrected inverses differ by up to 0.66 (`grey/573`, Y = 5.25), and the vectors mark these rows `info` (Section 9.4).

### B.3 ΔE from neutral-corrected coordinates

ΔE `helmlab` is defined on base-form coordinates. The reference code removes the neutral correction before measuring:

- Python `hl.metric.distance(lab_a, lab_b)` checks its inputs as `distance_from_lab` does (Section 5.6), converts both neutral-corrected triples to XYZ with the neutral-corrected inverse, and calls `MetricSpace.distance` (base-form forward, then ΔE). `hl.metric.difference` and `hl.metric.jnd` start from colour strings and compute base-form coordinates directly.
- JavaScript `distance(lab1, lab2)` undoes the rotation, adds the correction back with its own table, and evaluates Steps D1 to D7 on the result.

The two routes give the same ΔE where the neutral-corrected inverse returns its input, that is where forward(inverse(lab)) = lab. Outside that range they differ. For the neutral-corrected triples (−0.2, 0.05, 0) and (−0.25, 0.05, 0), Python `hl.metric.distance` gives 0.14965, while ΔE on the de-corrected coordinates (here equal to the input, because the correction is 0 for L ≤ 0) gives 0.08967.
