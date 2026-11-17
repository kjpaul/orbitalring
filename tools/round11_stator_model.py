#!/usr/bin/env python3
"""
Round 11, item 2: the model behind tools/round11_stator.py (a stator for the mass driver that carries a stated launch).
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 9.10, 9.11, 10.1 (Equations 10.1 to 10.3, Table 10.2), Section 14.5.
Run: imported by tools/round11_stator.py; `python tools/round11_stator_model.py` runs its self-checks only.

The winding is the full-pitch loop of Equation 9.3: N turns round a rectangle one pole pitch (tau) long and w wide, in a
two-layer three-phase winding with q coils per phase per pole pitch on each side of the sled.  What is modelled:

  thrust     F = sides x (1/2) K1 B0 sin(delta) w L_active, K1 = kd x 6 q N I / tau the fundamental of the current sheet (A/m).
             The same from the power balance: coils facing the sled x (1/2) E I sin(delta), E = 2 N kd v B0 w (Equation 9.3
             with the distribution factor kd of coil sides spread evenly over their 60 degree phase belts, kd = 3/pi).
  voltage    V = j omega (Psi_m + L_s I): Psi_m = N kd (2/pi) tau w B0 the sled's flux linkage, L_s the cyclic (synchronous)
             inductance of a coil: its self inductance plus the mutual inductance to every other coil of the side, each
             weighted by the cosine of its phase.  The inductances are Neumann's integral over straight filaments (closed
             form for parallel filaments, Grover, Inductance Calculations, 1946, ch. 6), with the geometric mean distance of
             the conductor band for the self terms.
  tape       lim_hysteresis.ic_tape (Robinson 70 K data of the AP tape), times a lift factor for 60, 50 and 30 K taken from the
             published Robinson measurement of a THEVA Pro-Line tape of 2019 (data/THEVA_ProLine_TPL4420c_Robinson_2019_...csv,
             Wimbush and Strickland, doi 10.6084/m9.figshare.3759327.v2): Ic(T, B, angle) / Ic(70 K, B, angle) of that sample.
             The rule is the bearing's: the worst tape at 50 percent of its limit in the field it sits in (the winding's own
             field at the peak current and the sled's field, amplitudes added).
  loss       Brandt and Indenbom's strip loss for the field through each tape's face plus Norris's transport loss
             (lim_hysteresis.q_strip_perp, q_norris_strip), tape by tape, in three regions: the active sides across the track,
             the end runs along the track (a flat band at each edge of the coil) and the legs that set the end runs back from
             the gap.  An eddy loss in the tape's metal is given beside it as a derived estimate (not in the totals).

The fields at the tapes are two-dimensional: the active sides as the phase belts of a periodic sheet (x-z plane), the end
runs as two flat bands carrying the net edge current kd K_belt tau / pi (y-z plane), each by lim_hysteresis.rect_field.
The sled's field is the book's: B0 sin(pi (x - v t) / tau) through the stator's face, uniform over the 2 m width.  Its
component along the track at the stator is taken as B0 kx / kappa, kappa = sqrt(kx^2 + ky^2), kx = pi / tau, ky = pi / w
(the lowest harmonic of a pole w wide) [UNVERIFIED: the sled's coils are not drawn anywhere in the repository].
"""
import functools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the repository root
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # round 13: this study was run on the old 130 m stator of lsm_config.py (its launch and its
                                                   # constants); the switch keeps it reproducing now that the 6 m stator is the default
import numpy as np

import lsm_config as cfg
import lim_hysteresis as lh

MU0 = 4e-7 * math.pi
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TPL19_CSV = os.path.join(HERE, "data", "THEVA_ProLine_TPL4420c_Robinson_2019_30K_50K_60K_70K.csv")
KD = math.sin(math.pi / 6) / (math.pi / 6)        # 3/pi, coil sides spread evenly over a 60 degree belt
SHARE = 0.5
W = cfg.W_COIL                                    # m, the width of the sled's field and of the coil
G = cfg.G_GAP
B0N = cfg.B_SLED_NOMINAL
VLIM = cfg.V_COIL_LIMIT
DELTA = cfg.DELTA_MAX
SITE = 500.0
ACTIVE = cfg.N_UNITS * cfg.N_POLES_PER_UNIT * cfg.TAU_P / cfg.L_SLED   # 0.975: the share of the sled's length that faces stator (39,000 m2 on the 10 km sled)
N_SITES = cfg.L_RING / SITE
T_HOT = 400.0
CARNOT_SHARE = 0.05                               # Chapter 14: the plant reaches 5 percent of Carnot between the cold side and 400 K
# the tape's metal for the eddy estimate [UNVERIFIED: layer build and resistivities, see the report]
T_HAST, RHO_HAST = 80e-6, 1.24e-6                 # m, ohm m (Hastelloy C-276 near 77 K, Lu et al., J. Appl. Phys. 103, 064908, 2008)
T_AG, RHO_AG = 2e-6, 3.0e-9                       # m, ohm m (silver cap near 77 K)
EDDY_G = T_HAST / RHO_HAST + T_AG / RHO_AG        # siemens, sheet conductance of the metal


def work_per_watt(T):
    """compressor work per watt lifted from T to the 400 K radiator at 5 percent of Carnot (Equation 14.1's rule)"""
    return (T_HOT - T) / T / CARNOT_SHARE


# ------------------------------------------------------------------------------------------------ the tape at other temperatures
@functools.lru_cache(maxsize=None)
def tpl19():
    import csv
    data = {}
    with open(TPL19_CSV, encoding="utf-8") as f:
        rd = csv.reader(f); next(rd)
        for T, B, a, i in rd:
            data.setdefault(float(T), {}).setdefault(float(B), []).append((float(a), float(i)))
    return {T: {B: (np.array([q[0] for q in sorted(v)]), np.array([q[1] for q in sorted(v)])) for B, v in d.items()} for T, d in data.items()}


def ic19(T, B, theta_deg):
    """Ic/w (A/cm) of the 2019 Pro-Line sample at T (30, 50, 60 or 70 K), folded and interpolated as lim_hysteresis.ic_tape does"""
    rob = tpl19()[float(T)]; Bs = np.array(sorted(rob))
    th = np.abs(np.asarray(theta_deg, float)); th = np.where(th > 90.0, 180.0 - th, th)
    B = np.asarray(B, float)
    per_B = np.array([np.minimum(np.interp(th, *rob[b]), np.interp(180.0 - th, *rob[b])) for b in Bs])
    Bc = np.clip(B, Bs[0], Bs[-1]); k = np.clip(np.searchsorted(Bs, Bc) - 1, 0, len(Bs) - 2)
    t = (Bc - Bs[k]) / (Bs[k + 1] - Bs[k]); rows = np.arange(len(np.atleast_1d(th)))
    return per_B[k, rows] * (1 - t) + per_B[k + 1, rows] * t


def lift(T, B, theta_deg):
    if abs(T - 70.0) < 1e-9:
        return np.ones_like(np.asarray(B, float))
    return ic19(T, B, theta_deg) / ic19(70.0, B, theta_deg)


def ic(B, theta_deg, w_t, T=70.0):
    """critical current (A) of the AP tape of width w_t: the measured 70 K value times the 2019 sample's lift factor at T"""
    B = np.atleast_1d(np.asarray(B, float)); th = np.atleast_1d(np.asarray(theta_deg, float))
    return lh.ic_tape(B, th, w_t) * lift(T, B, th)


# ------------------------------------------------------------------------------------------------ inductance (Neumann, straight filaments)
def _G(u, d):
    return u * np.arcsinh(u / d) - np.sqrt(u * u + d * d)


def m_par(a1, b1, a2, b2, d):
    """mutual inductance (H) of two parallel filaments, one from a1 to b1 and one from a2 to b2 along the same axis, a
    distance d apart: (mu0 / 4 pi) x double integral of dx1 dx2 / sqrt((x1 - x2)^2 + d^2), in closed form"""
    return MU0 / (4 * math.pi) * (_G(b2 - a1, d) - _G(b2 - b1, d) - _G(a2 - a1, d) + _G(a2 - b1, d))


def m_loops(X, tau, w, g_act, g_end):
    """mutual inductance of two one-turn rectangular loops tau x w in the same plane, the second shifted X along the track;
    X = 0 gives the self inductance.  g_act and g_end are the geometric mean distances of the active sides and of the end runs
    (used where two conductors lie on the same line)."""
    X = np.asarray(X, float)
    runs = 2 * m_par(0.0, tau, X, X + tau, g_end) - 2 * m_par(0.0, tau, X, X + tau, math.hypot(w, g_end))
    def side(dx):
        return m_par(-w / 2, w / 2, -w / 2, w / 2, np.maximum(np.abs(dx), g_act))
    sides = side(X) + side(X) - side(X + tau) - side(X - tau)       # left-left, right-right, left-right, right-left
    return runs + sides


def grover_rectangle(a, b, r):
    """self inductance of a one-turn rectangle a x b of a thin conductor with geometric mean radius r (Grover 1946, the
    rectangle of round wire without the internal term)"""
    dgl = math.hypot(a, b)
    return MU0 / math.pi * (a * math.log(2 * a / r) + b * math.log(2 * b / r) + 2 * dgl - a * math.asinh(a / b) - b * math.asinh(b / a) - 2 * (a + b))


def cyclic_inductance(tau, w, q, g_act, g_end, n_pp=12):
    """per-turn-squared cyclic inductance (H) of a coil in the two-layer winding with q coils per phase per pole pitch:
    sum over the coils of the side, shifted j tau / (3 q), of M(j) cos(pi j / (3 q)); returns (L_cyclic, L_self)"""
    nq = max(1, int(round(3 * q)))
    j = np.arange(-nq * n_pp, nq * n_pp + 1)
    X = j * tau / nq
    M = m_loops(X, tau, w, g_act, g_end)
    return float(np.sum(M * np.cos(math.pi * j / nq))), float(M[len(j) // 2])


# ------------------------------------------------------------------------------------------------ the fields at the tapes, per unit current
@functools.lru_cache(maxsize=None)
def active_unit_fields(tau, lay, w_t, n_pos=24, n_per=3):
    """amplitudes of the winding's own field at the tapes of an active side per unit belt current density (T per A/m):
    (through the face, along the face), each n_pos x 3 (both edges and the centre of each sample tape)"""
    beta = tau / 3.0
    xs = (np.arange(n_pos) + 0.5) / n_pos * beta
    if lay == "edge":                                  # tapes on edge, face normal along the track
        z1, z2 = -w_t / 2, w_t / 2
        px = np.repeat(xs, 3); pz = np.tile(np.array([-w_t / 2, 0.0, w_t / 2]), n_pos)
    else:                                              # tapes flat in the stator's face, normal toward the sled
        t = 1.0e-3; z1, z2 = -t / 2, t / 2
        px = (xs[:, None] + np.array([-w_t / 2, 0.0, w_t / 2])[None, :]).ravel(); pz = np.full(3 * n_pos, t / 2)
    Bx = np.zeros(len(px), complex); Bz = np.zeros(len(px), complex)
    for b in range(-6 * n_per, 6 * n_per + 6):
        ph = np.exp(-1j * math.pi * b / 3.0)
        bx, bz = lh.rect_field(px, pz, b * beta, (b + 1) * beta, z1, z2, beta)
        Bx += ph * bx; Bz += ph * bz
    ax = np.abs(Bx).reshape(n_pos, 3); az = np.abs(Bz).reshape(n_pos, 3)
    return (ax, az) if lay == "edge" else (az, ax)


@functools.lru_cache(maxsize=None)
def end_unit_fields(b_e, w_t, lay="flat", t_e=2.0e-3):
    """amplitudes of the two end bands' field at the tapes of one band per ampere of net band current (T per A):
    (through the face, along the face).  lay 'flat': the tapes lie flat side by side across the band (faces toward the sled),
    at the pitch of the flat lay.  lay 'edge': the tapes stand on edge, face to face across the band (faces across the track),
    the band as tall as the tape is wide, as in a racetrack pancake."""
    if lay == "flat":
        pitch = w_t + lh.KAPTON + lh.ADHESIVE
        n = max(4, int(b_e / pitch))
        ys = W / 2 - b_e / 2 + (np.arange(n) + 0.5) / n * b_e
        py = (ys[:, None] + np.array([-w_t / 2, 0.0, w_t / 2])[None, :]).ravel(); pz = np.full(3 * n, t_e / 2)
        z1, z2 = -t_e / 2, t_e / 2
    else:
        n = 64
        ys = W / 2 - b_e / 2 + (np.arange(n) + 0.5) / n * b_e
        py = np.repeat(ys, 3); pz = np.tile(np.array([-w_t / 2, 0.0, w_t / 2]), n)
        z1, z2 = -w_t / 2, w_t / 2
    by, bz = lh.rect_field(py, pz, W / 2 - b_e / 2, W / 2 + b_e / 2, z1, z2, 1.0)
    by2, bz2 = lh.rect_field(py, pz, -W / 2 - b_e / 2, -W / 2 + b_e / 2, z1, z2, -1.0)
    ay = np.abs(by + by2).reshape(n, 3); az = np.abs(bz + bz2).reshape(n, 3)
    return (az, ay) if lay == "flat" else (ay, az)


@functools.lru_cache(maxsize=None)
def flux_per_turn_unit(tau, b_e, q_m=4):
    """the stator's own flux through one turn per unit belt current density (Wb per A/m), from the cyclic inductance with
    q_m coils per phase per pole pitch: L_cyclic x (ampere-turns per coil) with N I = K_belt tau / (6 q_m)"""
    g_act = 0.2235 * tau / (3.0 * q_m); g_end = 0.2235 * b_e
    Lc, Ls = cyclic_inductance(tau, W, q_m, g_act, g_end)
    return Lc * tau / (6.0 * q_m), Lc, Ls


def region(BnK, BtK, K, Bn_s, Bt_s, i_t, w_t, T):
    """one region at sheet current K: the worst sample's limit, and the means over its tapes of the strip loss, the Norris
    loss (J/m per tape per cycle), the field through the face and its square"""
    Bn = BnK * K + Bn_s; Bt = BtK * K + Bt_s
    B = np.hypot(Bn, Bt); th = np.degrees(np.arctan2(Bt, Bn))
    ic_s = ic(B.ravel(), th.ravel(), w_t, T)
    Bm = B.mean(axis=1); Bnm = Bn.mean(axis=1)
    thm = np.degrees(np.arccos(np.clip(Bnm / np.maximum(Bm, 1e-30), 0, 1)))
    Ic = ic(Bm, thm, w_t, T)
    qp = lh.q_strip_perp(Ic, w_t, Bnm); qn = lh.q_norris_strip(i_t, Ic)
    return dict(ic_min=float(ic_s.min()), qp=float(qp.mean()), qn=float(qn.mean()), Bn=float(Bnm.mean()), Bn2=float((Bnm**2).mean()), Bmax=float(B.max()), ic_mean=float(Ic.mean()))


class Stator:
    """One stator: geometry and lay.  tau pole pitch (m); lay 'edge' or 'flat' for the active sides, end_lay the same for
    the end runs; w_t tape width (m);
    T temperature (K); b_e width of the flat end band at each edge (m); d setback of the end runs behind the face (m);
    sides 1 or 2; c_end the share of the sled's field that passes through the faces of the end-run tapes at the face."""

    def __init__(self, tau, lay="edge", w_t=0.012, T=70.0, b_e=0.5, d=0.0, sides=2, c_end=0.5, end_lay="flat"):
        self.tau, self.lay, self.w_t, self.T, self.b_e, self.d, self.sides, self.c_end = tau, lay, w_t, T, b_e, d, sides, c_end
        self.end_lay = end_lay
        self.kx = math.pi / tau; self.ky = math.pi / W; self.kappa = math.hypot(self.kx, self.ky)
        self.rx = self.kx / self.kappa                          # the sled's along-track field at the stator over B0
        self.rz = 1.0                                           # the sled's field through the stator's face over B0 (the book's sine wave)
        self.aBn, self.aBt = active_unit_fields(tau, lay, w_t)
        eBn, eBt = end_unit_fields(b_e, w_t, end_lay)
        self.i_net = KD * tau / math.pi                         # net band current per unit belt current density (A per A/m)
        self.eBn, self.eBt = eBn * self.i_net, eBt * self.i_net
        self.phi_unit, self.Lc, self.Ls = flux_per_turn_unit(tau, b_e)
        self.B_avg = self.phi_unit / (KD * (2 / math.pi) * tau * W)   # the stator's own mean field through a coil per unit K_belt
        kd_ = self.ky * d
        self.c_run = c_end * math.exp(-kd_)                     # sled field share at the end runs, set back d
        self.leg_decay = (1 - math.exp(-kd_)) / kd_ if kd_ > 0 else 1.0   # mean of the sled's field over the depth of a leg

    def state(self, K, B0, i_t):
        """the three regions at belt current density K (A/m), sled field B0 (T) and tape current i_t (A)"""
        own_n = np.maximum(self.aBt, self.B_avg) if self.lay == "edge" else np.maximum(self.aBn, self.B_avg)
        if self.lay == "edge":
            act = region(self.aBn, own_n, K, B0 * self.rx, B0 * self.rz, i_t, self.w_t, self.T)
        else:
            act = region(own_n, self.aBt, K, B0 * self.rz, B0 * self.rx, i_t, self.w_t, self.T)
        end = region(self.eBn, self.eBt, K, B0 * self.c_run, B0 * self.c_run, i_t, self.w_t, self.T)
        # the legs: tapes turned so that their faces look along the track (a quarter twist at each corner), so the sled's
        # field through the face is its along-track component; its other components lie along the face
        leg = region(self.eBn, self.eBt, K, B0 * self.rx * self.leg_decay, B0 * self.c_end * self.leg_decay, i_t, self.w_t, self.T) if self.d > 0 else end
        return act, end, leg

    def size(self, K_pk, B0=B0N):
        """the tape current at the 50 percent rule at the peak sheet current and the full sled field, and the tapes per meter"""
        act, end, leg = self.state(K_pk, B0, 1.0)
        ic_min = min(act["ic_min"], end["ic_min"], leg["ic_min"])
        i_t = SHARE * ic_min
        return i_t, K_pk / i_t, dict(act=act["ic_min"], end=end["ic_min"], leg=leg["ic_min"], B_act=act["Bmax"], B_end=end["Bmax"])

    def loss_per_m(self, K, B0, i_t, n_t):
        """loss per meter of track of one side per cycle (J/m), split: active, end runs, legs; and the sum of tape length x
        mean square field through the face (m T^2 per meter of track) for the eddy estimate"""
        act, end, leg = self.state(K, B0, i_t)
        qa = n_t * W * (act["qp"] + act["qn"]); qe = n_t * self.tau * (end["qp"] + end["qn"]); ql = n_t * 2 * self.d * (leg["qp"] + leg["qn"])
        b2 = n_t * (W * act["Bn2"] + self.tau * end["Bn2"] + 2 * self.d * leg["Bn2"])
        return qa, qe, ql, b2

    def tape_per_m(self, n_t):
        """meters of tape per meter of track of one side, and kg"""
        m = n_t * (W + self.tau + 2 * self.d)
        return m, m * self.w_t * lh.TAPE_T * lh.TAPE_RHO


def launch_profile(F_pk, P, m, v_f, n=3000):
    """the launch on the v axis: thrust F_pk until F v reaches P, then P / v.  Returns mid speeds, dt, thrust at each."""
    ve = np.linspace(0.0, v_f, n + 1); v = 0.5 * (ve[1:] + ve[:-1]); dv = ve[1] - ve[0]
    F = np.minimum(F_pk, P / v)
    dt = m * dv / F
    return v, dt, F


def kinematics(F_pk, P, m, v_f, marks=()):
    v, dt, F = launch_profile(F_pk, P, m, v_f)
    t = np.cumsum(dt); x = np.cumsum(v * dt)
    hit = {mk: (float(np.interp(mk, v, t)), float(np.interp(mk, v, x))) for mk in marks}
    return dict(t=float(t[-1]), x=float(x[-1]), marks=hit)


def evaluate(st, F_pk, P, L_sled, m_sc, v_f, alpha=0.0, n_s=8):
    """One stator on one launch.  F_pk peak thrust (N), P power cap (W), L_sled (m), m_sc spacecraft (kg), v_f release speed
    over the casing (m/s).  alpha splits the constant-power range between the sled's field and the stator's current:
    B0 = B0N s^alpha, K = K_pk s^(1 - alpha), s = min(1, v_c / v); alpha = 1 is the book's rule (the sled's field lowered,
    the coil current held), alpha = 0 holds the sled's field and lowers the coil current."""
    tau = st.tau
    L_act = L_sled * ACTIVE
    m = cfg.M_SLED_PER_M * L_sled + m_sc
    K1 = 2 * F_pk / (st.sides * B0N * math.sin(DELTA) * W * L_act)
    Kb = K1 / KD
    i_t, n_t, lim = st.size(Kb)
    v_c = P / F_pk
    s_min = min(1.0, v_c / v_f)
    sg = np.geomspace(s_min, 1.0, n_s) if s_min < 1 else np.array([1.0, 1.0])
    rows = np.array([st.loss_per_m(Kb * s**(1 - alpha), B0N * s**alpha, i_t * s**(1 - alpha), n_t) for s in sg])
    v, dt, F = launch_profile(F_pk, P, m, v_f)
    s = np.minimum(1.0, v_c / v)
    f = v / (2 * tau)
    def along(col):
        y = rows[:, col]
        if s_min >= 1:
            return np.full_like(v, y[-1])
        return np.exp(np.interp(np.log(s), np.log(sg), np.log(np.maximum(y, 1e-300))))
    qa, qe, ql, b2 = along(0), along(1), along(2), along(3)
    scale = st.sides * L_act
    E_a = scale * float(np.sum(qa * f * dt)); E_e = scale * float(np.sum(qe * f * dt)); E_l = scale * float(np.sum(ql * f * dt))
    E_h = E_a + E_e + E_l
    low = v < v_c
    E_low = scale * float(np.sum(((qa + qe + ql) * f * dt)[low])); x_c = float(np.sum((v * dt)[low]))
    p_eddy = (math.pi**2 / 6.0) * f**2 * st.w_t**3 * EDDY_G * b2          # W per meter of track of one side
    E_eddy = scale * float(np.sum(p_eddy * dt))
    # the voltage: flux per turn of the sled and of the stator along the launch
    B0v = B0N * s**alpha; Kv = Kb * s**(1 - alpha)
    phi_m = KD * (2 / math.pi) * tau * W * B0v; phi_s = st.phi_unit * Kv
    om = math.pi * v / tau
    v_turn_dem = om * np.sqrt(phi_m**2 + phi_s**2 - 2 * phi_m * phi_s * math.cos(DELTA))
    v_turn_mag = om * np.sqrt(phi_m**2 + phi_s**2 + 2 * phi_m * phi_s * math.cos(DELTA))
    e_turn = om * phi_m
    # the turns: the most that keep the terminal voltage under load, and the back-EMF alone (a coil that has lost its
    # current), inside the limit through the whole launch
    N = int(VLIM // float(max(v_turn_dem.max(), e_turn.max())))
    k_end = len(v) - 1; k_c = int(np.searchsorted(v, v_c)) if v_c < v_f else k_end
    vc_ = min(v_c, v_f); om_c = math.pi * vc_ / tau
    pm_c = KD * (2 / math.pi) * tau * W * B0N; ps_c = st.phi_unit * Kb                      # exactly at the corner speed
    Ec, Xc = om_c * pm_c, om_c * ps_c
    Vc_dem = math.sqrt(Ec**2 + Xc**2 - 2 * Ec * Xc * math.cos(DELTA)); Vc_mag = math.sqrt(Ec**2 + Xc**2 + 2 * Ec * Xc * math.cos(DELTA))
    kin = kinematics(F_pk, P, m, v_f)
    KE = 0.5 * m * v_f**2
    P_site = lambda k: st.sides * L_act * float((qa + qe + ql)[k] * f[k]) / (L_sled / SITE)
    tape_m, tape_kg = st.tape_per_m(n_t)
    wpw = work_per_watt(st.T)
    out = dict(tau=tau, lay=st.lay, end_lay=st.end_lay, w_t=st.w_t, T=st.T, b_e=st.b_e, d=st.d, sides=st.sides, alpha=alpha, c_end=st.c_end,
               F_pk=F_pk, P=P, L_sled=L_sled, m=m, v_f=v_f, v_c=v_c, K1=K1, K_belt=Kb, B_sheet=MU0 * K1 / 2, i_t=i_t, n_t=n_t,
               ic_act=lim["act"], ic_end=lim["end"], ic_leg=lim["leg"], B_act=lim["B_act"], B_end=lim["B_end"],
               tape_m_per_m=st.sides * tape_m, tape_kg_per_m=st.sides * tape_kg, tape_ring_t=st.sides * tape_kg * cfg.L_RING / 1e3,
               N=N, ampere_turns_pp=Kb * tau / 6.0, f_end=v_f / (2 * tau),
               E_turn_c=Ec, XI_turn_c=Xc, V_turn_c=Vc_dem, Vmag_turn_c=Vc_mag, E_turn_max=float(e_turn.max()),
               E_turn_f=float(om[k_end] * phi_m[k_end]), XI_turn_f=float(om[k_end] * phi_s[k_end]), V_turn_f=float(v_turn_dem[k_end]), Vmag_turn_f=float(v_turn_mag[k_end]),
               V_turn_max=float(v_turn_dem.max()), v_at_Vmax=float(v[int(v_turn_dem.argmax())]),
               t=kin["t"], x=kin["x"], laps=kin["x"] / cfg.L_RING, a_g=F_pk / m / cfg.G_0, KE=KE,
               E_hyst=E_h, E_act=E_a, E_end=E_e, E_leg=E_l, E_eddy=E_eddy, share=E_h / KE, E_site=E_h / N_SITES, E_low=E_low, x_c=x_c, E_site_start=E_low / max(x_c / SITE, 1.0),
               P_site_c=P_site(k_c), P_site_f=P_site(k_end), E_pass_f=P_site(k_end) * L_sled / v_f,
               work_per_watt=wpw, E_comp=E_h * wpw, comp_share=E_h * wpw / KE, Lc=st.Lc, Ls=st.Ls, B_avg_pk=st.B_avg * Kb,
               q_pk=float(rows[-1, :3].sum()), q_end_speed=float(rows[0, :3].sum()))
    return out


def coil_split(r, I_coil=cfg.I_PEAK):
    """a coil for the result r: layers that bring the coil current nearest I_coil, the coil current, the coils per phase per
    pole pitch per side, the coils facing the sled, and the voltages of a coil of N turns"""
    lay = max(1, int(round(I_coil / r["i_t"])))
    I = lay * r["i_t"]
    N = max(r["N"], 1)
    q = r["ampere_turns_pp"] / (N * I)
    n = r["sides"] * 3 * q * r["L_sled"] * ACTIVE / r["tau"]
    return dict(layers=lay, I=I, q=q, coils=n, E_c=N * r["E_turn_c"], XI_c=N * r["XI_turn_c"], V_c=N * r["V_turn_c"],
                E_f=N * r["E_turn_f"], XI_f=N * r["XI_turn_f"], V_f=N * r["V_turn_f"], V_max=N * r["V_turn_max"], E_max=N * r["E_turn_max"],
                P_balance_c=n * 0.5 * N * r["E_turn_c"] * I * math.sin(DELTA))


def selfcheck():
    print("SELF-CHECKS OF tools/round11_stator_model.py")
    # 1. the parallel-filament formula against a direct quadrature
    a1, b1, a2, b2, d = 0.0, 3.0, 1.0, 5.5, 0.4
    x1 = np.linspace(a1, b1, 4001); x2 = np.linspace(a2, b2, 4001)
    X1, X2 = np.meshgrid(0.5 * (x1[1:] + x1[:-1]), 0.5 * (x2[1:] + x2[:-1]))
    quad = MU0 / (4 * math.pi) * np.sum(1.0 / np.sqrt((X1 - X2)**2 + d * d)) * (x1[1] - x1[0]) * (x2[1] - x2[0])
    print(f"  parallel filaments: closed form {float(m_par(a1, b1, a2, b2, d)):.6e} H, quadrature {quad:.6e} H")
    # 2. the rectangle's self inductance against Grover's formula
    for a, b, r in ((130.0, 2.0, 0.146), (10.0, 2.0, 0.05), (2.0, 2.0, 0.01)):
        print(f"  rectangle {a:g} x {b:g} m, r {r:g} m: Neumann {float(m_loops(0.0, a, b, r, r))*1e6:.3f} uH, Grover {grover_rectangle(a, b, r)*1e6:.3f} uH")
    # 3. the cyclic inductance against the two-dimensional sheet for a pole pitch well under the width:
    #    flux per turn = kd (2/pi) tau w mu0 K1 / 2 with K1 = 6 q N I kd / tau, so L = (6/pi) mu0 q kd^2 w per turn squared
    for tau, w, q in ((0.25, 2.0, 4), (0.5, 2.0, 4)):
        kdq = math.sin(math.pi / 6) / (q * math.sin(math.pi / (6 * q)))
        Lc, Ls = cyclic_inductance(tau, w, q, 0.2235 * tau / (3 * q), 0.02)
        print(f"  cyclic inductance, tau {tau:g} m, w {w:g} m, q {q}: Neumann {Lc*1e6:.3f} uH, sheet (6/pi) mu0 q kd^2 w = {6/math.pi*MU0*q*kdq**2*w*1e6:.3f} uH")
    # 4. the field at an on-edge sheet: mu0 K / 2 at the tape edges
    ax, az = active_unit_fields(10.0, "edge", 0.012)
    print(f"  on-edge sheet, per A/m: through the face at the tape's edges {ax[12, 0]:.4e}, {ax[12, 2]:.4e}, centre {ax[12, 1]:.2e}; mu0 / 2 = {MU0/2:.4e}")
    # 5. the lift factor
    for T in (60.0, 50.0, 30.0):
        print(f"  lift factor at {T:.0f} K (2019 sample), field through the face: " + ", ".join(f"{B:g} T {float(lift(T, np.array([B]), np.array([0.0]))[0]):.2f}" for B in (0.0, 0.03, 0.1, 0.3, 1.0)))


if __name__ == "__main__":
    selfcheck()
