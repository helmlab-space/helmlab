"""Clean-room GenSpace (helmgen, helmgenlch) written only from spec Section 5.1.

Every formula is evaluated in the order written in the description:
'**' first, then unary minus, then '*' and '/', then '+' and '-', left to right.
x**2 is x*x; every other power is pow (np.power).

MATVEC_MODE selects how 'M . v' is evaluated:
  'explicit' : M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2], left to right (literal text, default)
  'blas'     : numpy matmul v @ M.T (the description says the reference uses numpy/BLAS)
"""
import json
import math
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
PARAMS_PATH = os.path.join(ROOT, 'spec', '2026', 'params', 'gen_params.json')

MATVEC_MODE = 'explicit'
# 'text' : use the M1_inv / M2_inv values printed in 5.1.3; 'numpy' : numpy.linalg.inv
INV_SOURCE = 'text'

with open(PARAMS_PATH) as _f:
    P = json.load(_f)

PI = math.pi                      # 3.141592653589793
TWO_PI = 2 * PI

# ---- 5.1.2 parameters ----
M1 = np.array(P['M1'], dtype=np.float64)
M2 = np.array(P['M2'], dtype=np.float64)
assert P['transfer'] == 'depcubic'
ALPHA = float(P['depcubic_alpha'])
L_CORR_PW = [float(v) for v in P['L_corr_pw']]
L_CORR_PW_STEP = float(P['L_corr_pw_step'])
ENR = P['enrichment']
assert ENR['type'] == 'L_gated_hue' and abs(ENR['amp']) > 1e-10
AMP = float(ENR['amp'])
CENTER_DEG = float(ENR['center_deg'])
SIGMA = float(ENR['sigma'])
L_LO = float(ENR['L_lo'])
L_HI = float(ENR['L_hi'])
CHROMA_POWER = float(P['chroma_power'])
assert CHROMA_POWER != 1

# ---- 5.1.3 derived constants ----
S = math.sqrt(ALPHA / 3)              # s = sqrt(alpha / 3)
K = 2 * S ** 3                        # k = 2 * s**3   (s**3 = pow(s, 3))
H_C = CENTER_DEG * (PI / 180)         # h_c = center_deg * (pi / 180)
CP_M1 = CHROMA_POWER - 1              # chroma_power - 1
CP_INV = 1 / CHROMA_POWER             # 1 / chroma_power
L_SPAN = L_HI - L_LO                  # enrichment.L_hi - enrichment.L_lo
S2 = SIGMA ** 2                       # s2 = sigma**2

N_PW = len(L_CORR_PW)                 # n = 19
P_IN = np.empty(N_PW + 2)
P_OUT = np.empty(N_PW + 2)
for _i in range(0, N_PW + 1):
    P_IN[_i] = _i * L_CORR_PW_STEP
P_IN[N_PW + 1] = 1.0
P_OUT[0] = 0.0
for _i in range(1, N_PW + 1):
    P_OUT[_i] = P_IN[_i] + L_CORR_PW[_i - 1]
P_OUT[N_PW + 1] = 1.0

M1_INV_TEXT = np.array([[1.2326502723725545, -0.5557414732179251, 0.2735612008453706],
                        [-0.040766587876246804, 1.1122097325799587, -0.0714431447037118],
                        [-0.07673215022426162, -0.4216188546558442, 1.5871810048801054]])
M2_INV_TEXT = np.array([[0.9930001151336145, 0.32599327253052285, 0.12945085631713896],
                        [0.9930001151336143, -0.08708353111074632, -0.038613617430049534],
                        [0.993000115133614, -0.12386097008215027, -0.8351991365871065]])
M1_INV_NP = np.linalg.inv(M1)
M2_INV_NP = np.linalg.inv(M2)


def _inv_mats():
    if INV_SOURCE == 'text':
        return M1_INV_TEXT, M2_INV_TEXT
    return M1_INV_NP, M2_INV_NP


# ---- 5.1.1 helpers ----
def matvec(M, v):
    """(M . v)[i] = M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2] for each row v of an (N, 3) array."""
    v = np.asarray(v, dtype=np.float64)
    if MATVEC_MODE == 'blas':
        return v @ M.T
    out = np.empty_like(v)
    for i in range(3):
        out[:, i] = M[i, 0] * v[:, 0] + M[i, 1] * v[:, 1] + M[i, 2] * v[:, 2]
    return out


def clip(x, lo, hi):
    """clip(x, lo, hi) = min(max(x, lo), hi), NaN-propagating."""
    return np.minimum(np.maximum(x, lo), hi)


def fmod_wrap(x, period):
    """e = fmod(x, period); if e < 0: e = e + period."""
    e = np.fmod(x, period)
    return np.where(e < 0, e + period, e)


# ---- forward steps ----
def F1(xyz):
    lms = matvec(M1, xyz)
    return np.maximum(lms, 0.0)


def F2(x):
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        t = x / K
        y0 = 2 * S * np.sinh(np.arcsinh(t) / 3)
        f = np.power(y0, 3.0) + ALPHA * y0 - x
        fp = 3 * (y0 * y0) + ALPHA
        fpp = 6 * y0
        den = 2 * fp * fp - f * fpp
        y_halley = y0 - 2 * f * fp / den
    return np.where(np.abs(den) > 1e-30, y_halley, y0)


def neutral_blend(c):
    """F3 / I5."""
    c0, c1, c2 = c[:, 0], c[:, 1], c[:, 2]
    m = (c0 + c1 + c2) / 3
    mx = np.maximum(np.maximum(c0, c1), c2)
    mn = np.minimum(np.minimum(c0, c1), c2)
    with np.errstate(over='ignore', invalid='ignore'):
        r = (mx - mn) / np.maximum(np.abs(m), 1e-30)
        q = r / 1e-5
        w = np.exp(-(q * q))
    out = np.empty_like(c)
    for i in range(3):
        ci = c[:, i]
        out[:, i] = ci + w * (m - ci)
    return out


def F4(c):
    return matvec(M2, c)


def _pw(L, xin, xout):
    """Piecewise-linear map used by F5 (xin=P_in, xout=P_out) and I3 (xin=P_out, xout=P_in)."""
    inr = (L >= 0.0) & (L <= 1.0)
    Ls = np.where(inr, L, 0.0)
    # j = the largest index in 0 .. n with xin[j] <= L
    j = np.searchsorted(xin, Ls, side='right') - 1
    j = np.clip(j, 0, N_PW)
    x0 = xin[j]
    x1 = xin[j + 1]
    y0 = xout[j]
    y1 = xout[j + 1]
    t = clip((Ls - x0) / np.maximum(x1 - x0, 1e-30), 0.0, 1.0)
    Lnew = y0 + t * (y1 - y0)
    return np.where(inr, Lnew, L)


def F5(L):
    return _pw(L, P_IN, P_OUT)


def _wrap_d(h):
    d = h - H_C
    e = fmod_wrap(d + PI, TWO_PI)
    return e - PI


def F6(L, a, b):
    C = np.sqrt(a * a + b * b)
    h = np.arctan2(b, a)
    u = clip((L - L_LO) / L_SPAN, 0.0, 1.0)
    sn = np.sin(PI * u)
    g = sn * sn
    d = _wrap_d(h)
    q = d / SIGMA
    G = np.exp(-0.5 * (q * q))
    h2 = h + AMP * g * G
    a2 = C * np.cos(h2)
    b2 = C * np.sin(h2)
    keep = C < 1e-12
    return np.where(keep, a, a2), np.where(keep, b, b2)


def F7(a, b):
    C = np.sqrt(a * a + b * b + 1e-30)
    q = np.power(C, CP_M1)
    return a * q, b * q


def from_XYZ(xyz):
    xyz = np.atleast_2d(np.asarray(xyz, dtype=np.float64))
    lms = F1(xyz)
    c = F2(lms)
    c = neutral_blend(c)
    lab = F4(c)
    L, a, b = lab[:, 0], lab[:, 1], lab[:, 2]
    L = F5(L)
    a, b = F6(L, a, b)
    a, b = F7(a, b)
    return np.stack([L, a, b], axis=-1)


# ---- inverse steps ----
def I1(a, b):
    C = np.sqrt(a * a + b * b + 1e-30)
    Cr = np.power(C, CP_INV)
    return a * (Cr / C), b * (Cr / C)


def I2(L, a, b, iters=8):
    C = np.sqrt(a * a + b * b)
    ht = np.arctan2(b, a)
    u = clip((L - L_LO) / L_SPAN, 0.0, 1.0)
    sn = np.sin(PI * u)
    g = sn * sn
    A = AMP * g
    h = ht.copy()
    for _ in range(iters):
        d = _wrap_d(h)
        q = d / SIGMA
        G = np.exp(-0.5 * (q * q))
        F = h + A * G - ht
        Fp = 1 + A * (G * (-d / S2))
        Fpp = A * (G * (-1 / S2 + d * d / (S2 * S2)))
        D = 2 * Fp * Fp - F * Fpp
        D = np.where(np.abs(D) < 1e-30, 1.0, D)
        h = h - 2 * F * Fp / D
    a2 = C * np.cos(h)
    b2 = C * np.sin(h)
    keep = C < 1e-12
    return np.where(keep, a, a2), np.where(keep, b, b2)


def I3(L):
    return _pw(L, P_OUT, P_IN)


def to_XYZ(lab):
    lab = np.atleast_2d(np.asarray(lab, dtype=np.float64))
    L, a, b = lab[:, 0].copy(), lab[:, 1].copy(), lab[:, 2].copy()
    a, b = I1(a, b)
    a, b = I2(L, a, b)
    L = I3(L)
    M1i, M2i = _inv_mats()
    c = matvec(M2i, np.stack([L, a, b], axis=-1))
    c = neutral_blend(c)
    lms = np.power(c, 3.0) + ALPHA * c
    return matvec(M1i, lms)


# ---- 5.1.6 LCh ----
RAD2DEG = 180 / PI
DEG2RAD = PI / 180


def to_lch(lab):
    lab = np.atleast_2d(np.asarray(lab, dtype=np.float64))
    L, a, b = lab[:, 0], lab[:, 1], lab[:, 2]
    C = np.hypot(a, b)
    h = np.arctan2(b, a) * RAD2DEG
    e = fmod_wrap(h, 360.0)
    h = e + 0.0          # -0 becomes +0
    return np.stack([L, C, h], axis=-1)


def from_lch(lch):
    lch = np.atleast_2d(np.asarray(lch, dtype=np.float64))
    L, C, h = lch[:, 0], lch[:, 1], lch[:, 2]
    a = C * np.cos(h * DEG2RAD)
    b = C * np.sin(h * DEG2RAD)
    return np.stack([L, a, b], axis=-1)


# ---- 5.1.4 sRGB (not part of GenSpace; used only for an extra check) ----
M_SRGB = np.array([[0.4124564390896922, 0.3575760776439511, 0.1804374832663989],
                   [0.2126728514056226, 0.7151521552878178, 0.0721750036064596],
                   [0.0193338955823293, 0.1191920258813418, 0.9503040785363679]])


def srgb_lin(v):
    v = np.asarray(v, dtype=np.float64)
    with np.errstate(invalid='ignore'):
        pw = np.power((v + 0.055) / 1.055, 2.4)
    return np.where(v <= 0.04045, v / 12.92, pw)


def srgb_to_xyz(rgb):
    rgb = np.atleast_2d(np.asarray(rgb, dtype=np.float64))
    return matvec(M_SRGB, srgb_lin(rgb))
