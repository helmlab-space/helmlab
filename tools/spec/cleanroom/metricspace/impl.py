"""Clean-room implementation of Helmlab 2026 MetricSpace (helmlab-metric) and
Delta-E helmlab, written ONLY from the description text (spec draft 5.3-5.5).

No code from src/helmlab or packages/helmlab-js was read.  Every step is named
as in the description (F1..F14, I1..I14, D1..D7, N1..N8).
"""
import json
import math
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
PARAMS_PATH = os.path.join(ROOT, 'spec', '2026', 'params', 'metric_params.json')

# Constants written in the description (not in the parameter file).
M1_INV_DESC = [[1.065107295808859, -0.3150044075301121, 0.1803492381741039],
               [0.4721107713837796, 0.4271995765962455, 0.04150680489380985],
               [-0.09156391926309541, -0.1736709363194979, 0.7719790382558294]]
M2_INV_DESC = [[0.9183897822815021, 0.5232051237088666, 0.7065804598090856],
               [1.0899090574433026, 0.07005324849041904, 0.4319776874787044],
               [1.598895729264209, 0.22061850068770233, 0.012506037355220951]]
PI = 3.141592653589793
EXP_CLIP = 30.0           # clip(g, -30, 30) in F7, F10, F11 (and I4, I5, I8)
HK_POW_LO, HK_POW_HI = 0.01, 10.0
NEWTON_GUARD = 1e-10
N_DARK, N_LCORR, N_HUE = 12, 15, 8
DUP_EPS = 1e-12


def load_params(path=PARAMS_PATH):
    with open(path) as f:
        return json.load(f)


# ----------------------------------------------------------------- helpers
def matvec(M, v0, v1, v2):
    """(M . v)[i] = M[i][0]*v[0] + M[i][1]*v[1] + M[i][2]*v[2] (5.3.1)."""
    return (M[0][0] * v0 + M[0][1] * v1 + M[0][2] * v2,
            M[1][0] * v0 + M[1][1] * v1 + M[1][2] * v2,
            M[2][0] * v0 + M[2][1] * v1 + M[2][2] * v2)


def signed_pow(x, e):
    """sign(x) * abs(x) ** e, sign(0) = 0."""
    return np.sign(x) * np.abs(x) ** e


def clip(x, lo, hi):
    """clip(x, lo, hi) = min(max(x, lo), hi)."""
    return np.minimum(np.maximum(x, lo), hi)


def hue_delta(P, x):
    return (P['hue_cos1'] * np.cos(x) + P['hue_sin1'] * np.sin(x)
            + P['hue_cos2'] * np.cos(2 * x) + P['hue_sin2'] * np.sin(2 * x)
            + P['hue_cos3'] * np.cos(3 * x) + P['hue_sin3'] * np.sin(3 * x)
            + P['hue_cos4'] * np.cos(4 * x) + P['hue_sin4'] * np.sin(4 * x))


def hue_delta_d(P, x):
    return (-P['hue_cos1'] * np.sin(x) + P['hue_sin1'] * np.cos(x)
            - 2 * P['hue_cos2'] * np.sin(2 * x) + 2 * P['hue_sin2'] * np.cos(2 * x)
            - 3 * P['hue_cos3'] * np.sin(3 * x) + 3 * P['hue_sin3'] * np.cos(3 * x)
            - 4 * P['hue_cos4'] * np.sin(4 * x) + 4 * P['hue_sin4'] * np.cos(4 * x))


def hk_factor(P, h5):
    return 1 + (P['hk_hue_mod'] * np.cos(h5) + P['hk_sin1'] * np.sin(h5)
                + P['hk_cos2'] * np.cos(2 * h5) + P['hk_sin2'] * np.sin(2 * h5))


def cs_scale(P, h):
    return np.exp(P['cs_cos1'] * np.cos(h) + P['cs_sin1'] * np.sin(h)
                  + P['cs_cos2'] * np.cos(2 * h) + P['cs_sin2'] * np.sin(2 * h)
                  + P['cs_cos3'] * np.cos(3 * h) + P['cs_sin3'] * np.sin(3 * h)
                  + P['cs_cos4'] * np.cos(4 * h) + P['cs_sin4'] * np.sin(4 * h))


def cp_exponent(P, h9):
    return 1 + (P['cp_cos1'] * np.cos(h9) + P['cp_sin1'] * np.sin(h9)
                + P['cp_cos2'] * np.cos(2 * h9) + P['cp_sin2'] * np.sin(2 * h9))


def hlc_q(P, h):
    return (P['hlc_cos1'] * np.cos(h) + P['hlc_sin1'] * np.sin(h)
            + P['hlc_cos2'] * np.cos(2 * h) + P['hlc_sin2'] * np.sin(2 * h))


def hl_scale(P, h):
    return np.exp(P['hl_cos1'] * np.cos(h) + P['hl_sin1'] * np.sin(h)
                  + P['hl_cos2'] * np.cos(2 * h) + P['hl_sin2'] * np.sin(2 * h))


def rotation(P):
    """Step F14 constants: phi = display_phi_deg * (pi / 180)."""
    phi = P['display_phi_deg'] * (PI / 180)
    return phi, math.cos(phi), math.sin(phi)


def hk_power_clipped(P):
    return min(max(P['hk_power'], HK_POW_LO), HK_POW_HI)


# ---------------------------------------------------------------- forward
def forward_to_F12(XYZ, P):
    """Steps F1..F12. Returns (L3, a5, b5)."""
    XYZ = np.asarray(XYZ, dtype=np.float64)
    X, Y, Z = XYZ[..., 0], XYZ[..., 1], XYZ[..., 2]
    # F1
    l, m, s = matvec(P['M1'], X, Y, Z)
    # F2
    g = P['gamma']
    lc = signed_pow(l, g[0])
    mc = signed_pow(m, g[1])
    sc = signed_pow(s, g[2])
    # F3
    L0, a0, b0 = matvec(P['M2'], lc, mc, sc)
    # F4
    h1 = np.arctan2(b0, a0)
    delta = hue_delta(P, h1)
    C1 = np.sqrt(a0 ** 2 + b0 ** 2)
    a1 = C1 * np.cos(h1 + delta)
    b1 = C1 * np.sin(h1 + delta)
    # F5
    C5 = np.sqrt(a1 ** 2 + b1 ** 2)
    h5 = np.arctan2(b1, a1)
    P5 = hk_power_clipped(P)
    F5 = hk_factor(P, h5)
    L_hk = L0 + (P['hk_weight'] * C5 ** P5) * F5
    # F6
    t = L_hk * (1 - L_hk)
    L1 = (L_hk + P['L_corr_p1'] * t + P['L_corr_p2'] * t * (0.5 - L_hk)
          + P['L_corr_p3'] * t * t
          + t * (P['Lh_cos1'] * np.cos(h5) + P['Lh_sin1'] * np.sin(h5)))
    # F7
    k = P['lp_dark'] + P['lp_dark_hcos'] * np.cos(h5) + P['lp_dark_hsin'] * np.sin(h5)
    gg = k * L1 * np.maximum(0.0, 1 - L1) ** 2
    L2 = L1 * np.exp(clip(gg, -EXP_CLIP, EXP_CLIP))
    # F8
    S8 = cs_scale(P, h5)
    a2 = a1 * S8
    b2 = b1 * S8
    # F9
    h9 = np.arctan2(b2, a2)
    C9 = np.sqrt(a2 ** 2 + b2 ** 2)
    p9 = cp_exponent(P, h9)
    with np.errstate(divide='ignore', invalid='ignore'):
        D9 = np.where(C9 > 0, C9 ** p9, 0.0)
    a3 = D9 * np.cos(h9)
    b3 = D9 * np.sin(h9)
    # F10
    u = L2 - 0.5
    T10 = np.exp(clip(P['lc1'] * u + P['lc2'] * u ** 2, -EXP_CLIP, EXP_CLIP))
    a4 = a3 * T10
    b4 = b3 * T10
    # F11
    h11 = np.arctan2(b4, a4)
    q11 = hlc_q(P, h11)
    H11 = np.exp(clip((L2 - 0.5) * q11, -EXP_CLIP, EXP_CLIP))
    a5 = a4 * H11
    b5 = b4 * H11
    # F12
    h12 = np.arctan2(b5, a5)
    G12 = hl_scale(P, h12)
    L3 = L2 * G12
    return L3, a5, b5


def forward(XYZ, P, neutral=None):
    """Base form (neutral=None) or neutral-corrected form (neutral = NeutralTable)."""
    L3, a5, b5 = forward_to_F12(XYZ, P)
    # F13: stage 9 skipped (L_S_offset = 0); stage 10 = neutral correction (N7)
    if neutral is not None:
        Lc = clip(L3, neutral.x[0], neutral.x[-1])
        a5 = a5 - neutral.Pa(Lc)
        b5 = b5 - neutral.Pb(Lc)
    # F14
    phi, cr, sr = rotation(P)
    L = L3
    if phi != 0:
        a = a5 * cr - b5 * sr
        b = a5 * sr + b5 * cr
    else:
        a, b = a5, b5
    return np.stack([L, a, b], axis=-1)


# ---------------------------------------------------------------- inverse
def inverse(Lab, P, neutral=None):
    Lab = np.asarray(Lab, dtype=np.float64)
    L, a, b = Lab[..., 0], Lab[..., 1], Lab[..., 2]
    # I1
    phi, cr, sr = rotation(P)
    if phi != 0:
        a5 = a * cr + b * sr
        b5 = -a * sr + b * cr
    else:
        a5, b5 = a, b
    L3 = L
    # I2 (neutral-corrected form only: N8)
    if neutral is not None:
        Lc = clip(L, neutral.x[0], neutral.x[-1])
        a5 = a5 + neutral.Pa(Lc)
        b5 = b5 + neutral.Pb(Lc)
    # I3
    h12 = np.arctan2(b5, a5)
    G12 = hl_scale(P, h12)
    L2 = L3 / G12
    # I4
    q11 = hlc_q(P, h12)
    H11 = np.exp(clip((L2 - 0.5) * q11, -EXP_CLIP, EXP_CLIP))
    a4 = a5 / H11
    b4 = b5 / H11
    # I5
    u = L2 - 0.5
    T10 = np.exp(clip(P['lc1'] * u + P['lc2'] * u ** 2, -EXP_CLIP, EXP_CLIP))
    a3 = a4 / T10
    b3 = b4 / T10
    # I6
    h9 = np.arctan2(b3, a3)
    C = np.sqrt(a3 ** 2 + b3 ** 2)
    p9 = cp_exponent(P, h9)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = np.where(C > 0, C ** (1 / p9), 0.0)
    a2 = D * np.cos(h9)
    b2 = D * np.sin(h9)
    # I7
    h8 = np.arctan2(b2, a2)
    S8 = cs_scale(P, h8)
    a1 = a2 / S8
    b1 = b2 / S8
    # I8
    h5 = np.arctan2(b1, a1)
    k = P['lp_dark'] + P['lp_dark_hcos'] * np.cos(h5) + P['lp_dark_hsin'] * np.sin(h5)
    x = L2
    for _ in range(N_DARK):
        m = np.maximum(0.0, 1 - x)
        g = k * x * m ** 2
        e = np.exp(clip(g, -EXP_CLIP, EXP_CLIP))
        f = x * e - L2
        gd = k * m * (1 - 3 * x)
        fd = e * (1 + x * gd)
        fd = np.where(np.abs(fd) < NEWTON_GUARD, 1.0, fd)
        x = x - f / fd
    L1 = x
    # I9
    Lh = P['Lh_cos1'] * np.cos(h5) + P['Lh_sin1'] * np.sin(h5)
    p1, p2, p3 = P['L_corr_p1'], P['L_corr_p2'], P['L_corr_p3']
    x = L1
    for _ in range(N_LCORR):
        t = x * (1 - x)
        dt = 1 - 2 * x
        f = x + p1 * t + p2 * t * (0.5 - x) + p3 * t * t + t * Lh - L1
        fd = 1 + (p1 + Lh) * dt + p2 * (dt * (0.5 - x) - t) + p3 * 2 * t * dt
        fd = np.where(np.abs(fd) < NEWTON_GUARD, 1.0, fd)
        x = x - f / fd
    L_hk = x
    # I10
    C5 = np.sqrt(a1 ** 2 + b1 ** 2)
    P5 = hk_power_clipped(P)
    F5 = hk_factor(P, h5)
    L0 = L_hk - (P['hk_weight'] * C5 ** P5) * F5
    # I11
    C1 = np.sqrt(a1 ** 2 + b1 ** 2)
    x = h5
    for _ in range(N_HUE):
        f = x + hue_delta(P, x) - h5
        fd = 1 + hue_delta_d(P, x)
        fd = np.where(np.abs(fd) < NEWTON_GUARD, 1.0, fd)
        x = x - f / fd
    a0 = C1 * np.cos(x)
    b0 = C1 * np.sin(x)
    # I12
    lc, mc, sc = matvec(M2_INV_DESC, L0, a0, b0)
    # I13
    g = P['gamma']
    l = signed_pow(lc, 1 / g[0])
    m = signed_pow(mc, 1 / g[1])
    s = signed_pow(sc, 1 / g[2])
    # I14
    X, Y, Z = matvec(M1_INV_DESC, l, m, s)
    return np.stack([X, Y, Z], axis=-1)


# ------------------------------------------------------------------ Delta-E
def delta_e(lab1, lab2, P):
    lab1 = np.asarray(lab1, dtype=np.float64)
    lab2 = np.asarray(lab2, dtype=np.float64)
    L_1, a_1, b_1 = lab1[..., 0], lab1[..., 1], lab1[..., 2]
    L_2, a_2, b_2 = lab2[..., 0], lab2[..., 1], lab2[..., 2]
    # D1
    dL2 = (L_1 - L_2) ** 2
    dab2 = (a_1 - a_2) ** 2 + (b_1 - b_2) ** 2
    # D2
    Lm = (L_1 + L_2) * 0.5
    C_1 = np.sqrt(a_1 ** 2 + b_1 ** 2)
    C_2 = np.sqrt(a_2 ** 2 + b_2 ** 2)
    Cm = (C_1 + C_2) * 0.5
    # D3
    SL = 1 + P['dist_sl'] * (Lm - 0.5) ** 2
    dL2w = dL2 / SL ** 2
    # D4
    SC = 1 + P['dist_sc'] * Cm
    dab2w = dab2 / SC ** 2
    # D5
    E0 = (dL2w + P['dist_wC'] * dab2w) ** (P['dist_power'] / 2)
    # D6 (dist_linear = 0, dist_compress != 0)
    E1 = E0 / (1 + P['dist_compress'] * E0)
    # D7
    return E1 ** P['dist_post_power']


# ------------------------------------------------------ neutral correction
def linspace_desc(s, e, n):
    """linspace(s, e, n)[j] = j*((e - s)/(n - 1)) + s, last = e (Step N1)."""
    step = (e - s) / (n - 1)
    out = [j * step + s for j in range(n - 1)] + [e]
    return np.array(out, dtype=np.float64)


def _sign(x):
    x = float(x)
    return int(x > 0) - int(x < 0)


def pchip_derivatives(x, y):
    """Step N6 derivative rules."""
    n = len(x)
    hk = [x[k + 1] - x[k] for k in range(n - 1)]
    mk = [(y[k + 1] - y[k]) / hk[k] for k in range(n - 1)]
    d = [0.0] * n
    for k in range(1, n - 1):
        w1 = 2 * hk[k] + hk[k - 1]
        w2 = hk[k] + 2 * hk[k - 1]
        if mk[k - 1] == 0 or mk[k] == 0 or _sign(mk[k - 1]) != _sign(mk[k]):
            d[k] = 0.0
        else:
            d[k] = 1 / ((w1 / mk[k - 1] + w2 / mk[k]) / (w1 + w2))

    def edge(h0, h1, m0, m1):
        e = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
        if _sign(e) != _sign(m0):
            return 0.0
        if _sign(m0) != _sign(m1) and abs(e) > 3 * abs(m0):
            return 3 * m0
        return e

    d[0] = edge(hk[0], hk[1], mk[0], mk[1])
    d[n - 1] = edge(hk[n - 2], hk[n - 3], mk[n - 2], mk[n - 3])
    return np.array(hk), np.array(mk), np.array(d)


class Pchip:
    def __init__(self, x, y):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.hk, self.mk, self.d = pchip_derivatives(list(self.x), list(self.y))
        n = len(self.x)
        # Per-interval coefficients (Step N6)
        hk, mk, d, y = self.hk, self.mk, self.d, self.y
        r = (d[:-1] + d[1:] - 2 * mk) / hk
        self.c0 = r / hk
        self.c1 = (mk - d[:-1]) / hk - r
        self.c2 = d[:-1].copy()
        self.c3 = y[:-1].copy()
        self.n = n

    def __call__(self, L):
        L = np.asarray(L, dtype=np.float64)
        k = np.searchsorted(self.x, L, side='right') - 1   # x_k <= L < x_{k+1}
        k = np.clip(k, 0, self.n - 2)                       # L = x_{n-1} -> n-2
        s = L - self.x[k]
        out = ((self.c3[k] + self.c2[k] * s) + self.c1[k] * (s * s)) + self.c0[k] * ((s * s) * s)
        return np.where(np.isnan(L), np.nan, out)


class NeutralTable:
    def __init__(self, P):
        # N1
        Y = np.concatenate([[0.0], linspace_desc(0.001, 0.01, 10),
                            linspace_desc(0.01, 0.1, 20),
                            linspace_desc(0.1, 2.0, 226)])
        assert len(Y) == 257
        self.Y = Y
        # N2
        XYZ = np.stack([Y * 0.95047, Y * 1.0, Y * 1.08883], axis=-1)
        # N3
        L3, a5, b5 = forward_to_F12(XYZ, P)
        # N4
        order = np.argsort(L3, kind='stable')
        Ls, As, Bs = L3[order], a5[order], b5[order]
        keep = [0] + [i for i in range(1, len(Ls)) if Ls[i] - Ls[i - 1] > DUP_EPS]
        x = Ls[keep]
        ya = As[keep]
        yb = Bs[keep]
        self.n_unique = len(x)
        self.x254 = x[-1]
        self.ya254 = ya[-1]
        self.yb254 = yb[-1]
        # N5
        ext = np.arange(x[-1] + 0.01, 2.6, 0.01)
        self.n_ext = len(ext)
        x = np.concatenate([x, ext])
        ya = np.concatenate([ya, np.full(len(ext), ya[-1])])
        yb = np.concatenate([yb, np.full(len(ext), yb[-1])])
        self.x = x
        self.ya = ya
        self.yb = yb
        # N6
        self.Pa = Pchip(x, ya)
        self.Pb = Pchip(x, yb)
