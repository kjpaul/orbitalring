#!/usr/bin/env python3
"""
The deployment LIM coil's winding as a set of tapes: the 2-D field each tape sees from the rest of the winding, the
measured critical current of each tape at its own field and angle (Robinson 70 K data), and the hysteresis loss of the
winding on the Brandt and Indenbom thin-strip model (the "brandt" option of lim_physics.calc_hysteresis_power).
Book: The Orbital Ring, Volume IV: Chapter 11 (Table 11.2, Sections 11.5 to 11.7), Chapter 18 (Table 18.2).
Run: imported by lim_config.py (when HYSTERESIS_MODEL = "brandt"), lim_physics.py and tools/lim_tape_limit.py; run alone for the checks:
     python lim_hysteresis.py
Environment: none (the Robinson data is data/THEVA_Robinson_70K_angle_dependence.csv)

The winding (tools/lim_tape_limit.py, 2026-10-26): lim_config.py budgets each turn at 1 mm of Kapton, 0.06 mm of adhesive
and 0.12 mm per layer of tape, the conductor being 0.08 mm per layer.  Two layouts: "flat", the tapes lying flat, edge to
edge, one turn beside the next (a sheet N x (w + 1.06 mm) wide); "pancake", the tapes on edge, face to face, one turn on the
next (a stack w wide and N x (1.06 + 0.12 x layers) mm high); and "compact", the first check's pancake with 0.05 mm of
insulation.  The field is the exact 2-D field of each conductor rectangle with uniform current density, summed over every
conductor but the tape's own (the measured limit was taken on a tape carrying its own current, so its self-field is inside
it).  The reaction plate is non-magnetic and its eddy currents are neglected (conservative); the other coil side 100 m away
and the neighbouring phase belts 33 m away add under 0.3 mT and are neglected.

Orientation (corrected 2026-10-27): in both layouts the bundle lies along the coil's surface (x) and is thin toward the plate (y),
and the field of a wide thin sheet of current runs along the sheet, that is along x.  The flat tapes have their faces normal
to y, so that field runs along their faces and only the end turns see much of it through the face; the pancake's tapes have
their faces normal to x, so the same field goes through every face.  The 2026-10-26 run of tools/lim_tape_limit.py had the
pancake's rectangles laid the other way (w along x, stacked along y) while taking x as the normal, so its pancake angles and
limits were those of a tape whose face is turned 90 degrees from where it is; the flat results were right.  The pancake
figures in that report (3 layers, 128.6 A worst tape, 184.6 A cap) are superseded by the ones in the 2026-10-27 report.

The loss (Brandt and Indenbom, Phys. Rev. B 48, 12893, 1993; Zeldov et al. 1994; reviewed by Grilli et al., IEEE Trans.
Appl. Supercond. 24, 8200433, 2014): a thin strip of half-width a, thickness d and critical current density Jc (Ic = 2 a d Jc)
in a perpendicular field of amplitude H_a loses, per metre and per cycle,

    Q_perp = 4 mu0 Jc d a^2 H_a [ (2/h) ln cosh h - tanh h ],   h = H_a / H_c,   H_c = Jc d / pi,

which is w Ic B_a for h >> 1 and goes as H_a^4 (h^3 times H_a) for h << 1.  For each tape H_a is the component of the
winding's field through the tape's face (B cos theta, theta from the normal), and Jc d = Ic / w with Ic the tape's measured
critical current at its own field and angle.  The component along the face gives a slab loss proportional to the film
thickness, under a thousandth of the strip loss at these fields, and is dropped.  The transport current's own loss is
Norris's strip formula per tape, added as an upper bound (the exact in-phase solution is lower).
"""
import math, os, functools
import numpy as np

MU0 = 4e-7 * math.pi
TAPE_T = 8e-05                     # m, the conductor thickness per layer (bearing_design.TAPE_T)
TAPE_RHO = 8890.0                  # kg/m3 (bearing_design.TAPE_RHO)
KAPTON, ADHESIVE, T_LAYER_BUDGET = 1.0e-3, 0.06e-3, 0.12e-3     # lim_config's per-turn budget
T_INS_COMPACT = 0.05e-3
HERE = os.path.dirname(os.path.abspath(__file__))
ROBINSON_CSV = os.path.join(HERE, "data", "THEVA_Robinson_70K_angle_dependence.csv")


# ---------------------------------------------------------------- the tape data
@functools.lru_cache(maxsize=None)
def robinson():
    """field (T) -> (angles deg, Ic/w A/cm) at 70 K from the Robinson dataset"""
    import csv
    data = {}
    with open(ROBINSON_CSV, encoding="utf-8-sig") as f:
        rd = csv.reader(f); next(rd)
        for B, a, i in rd:
            data.setdefault(float(B), []).append((float(a), float(i)))
    return {B: (np.array([q[0] for q in sorted(v)]), np.array([q[1] for q in sorted(v)])) for B, v in data.items()}


def ic_tape(B, theta_deg, w):
    """critical current (A) of a tape of width w (m) at field B (T) at theta degrees from its normal: the lower of theta and
    180 - theta in the Robinson data, linear between the measured fields (bearing_optimize.ic_angle_full)"""
    rob = robinson(); Bs = np.array(sorted(rob))
    th = np.abs(np.asarray(theta_deg, float)); th = np.where(th > 90.0, 180.0 - th, th)
    B = np.asarray(B, float)
    per_B = np.array([np.minimum(np.interp(th, *rob[b]), np.interp(180.0 - th, *rob[b])) for b in Bs])
    Bc = np.clip(B, Bs[0], Bs[-1]); k = np.clip(np.searchsorted(Bs, Bc) - 1, 0, len(Bs) - 2)
    t = (Bc - Bs[k]) / (Bs[k + 1] - Bs[k]); rows = np.arange(len(np.atleast_1d(th)))
    return (per_B[k, rows] * (1 - t) + per_B[k + 1, rows] * t) * w * 100.0


# ---------------------------------------------------------------- the field of the winding
def _G(u, v):
    r2 = u * u + v * v
    vv = np.where(v == 0, 1e-300, v)
    return 0.5 * u * np.log(np.maximum(r2, 1e-30)) - u + v * np.arctan(u / vv)


def rect_field(x, y, x1, x2, y1, y2, I):
    """Bx, By (T) at points x, y of a rectangle x1..x2, y1..y2 carrying current I (A) uniformly, straight along z"""
    J = I / ((x2 - x1) * (y2 - y1))
    def S(F):
        return F(x - x1, y - y1) - F(x - x2, y - y1) - F(x - x1, y - y2) + F(x - x2, y - y2)
    return -MU0 * J / (2 * math.pi) * S(_G), MU0 * J / (2 * math.pi) * S(lambda u, v: _G(v, u))


def layout(kind, N, L, w):
    """conductor rectangles (x1, x2, y1, y2) of the N turns (each L layers thick), the tape normal ('x' or 'y'), bundle width and height"""
    tc = L * TAPE_T
    if kind == "flat":
        pitch = w + KAPTON + ADHESIVE
        return [(k * pitch, k * pitch + w, 0.0, tc) for k in range(N)], "y", pitch * N, tc
    if kind == "pancake":                   # tapes on edge, face to face along x, one turn beside the next; width w along y (toward the plate)
        pitch = KAPTON + ADHESIVE + L * T_LAYER_BUDGET
        return [(k * pitch, k * pitch + tc, 0.0, w) for k in range(N)], "x", pitch * N, w
    pitch = L * TAPE_T + T_INS_COMPACT      # the compact pancake of the first check, the same orientation
    return [(k * pitch, k * pitch + tc, 0.0, w) for k in range(N)], "x", pitch * N, w


def winding_width(kind, N, L, w):
    """the width of the winding across the coil's surface, the cryostat band (m): the flat sheet's width or the pancake's stack width"""
    return layout(kind, N, L, w)[2]


@functools.lru_cache(maxsize=None)
def sample_fields(kind, N, L, w):
    """the field at every tape of the winding at a coil current of 1 A: three samples per tape (both edges and the centre of
    its width, on its mid-plane), from every conductor but the tape's own.  Returns tape index, normal component, tangential
    component (T per A), the tape count, bundle width and height."""
    rects, normal, bw, bh = layout(kind, N, L, w)
    tapes = []
    for (x1, x2, y1, y2) in rects:
        for j in range(L):
            tapes.append((x1, x2, y1 + j * TAPE_T, y1 + (j + 1) * TAPE_T) if normal == "y" else (x1 + j * TAPE_T, x1 + (j + 1) * TAPE_T, y1, y2))
    tapes = np.array(tapes); n = len(tapes)
    if normal == "y":
        px = np.concatenate([tapes[:, 0], 0.5 * (tapes[:, 0] + tapes[:, 1]), tapes[:, 1]]); py = np.tile(0.5 * (tapes[:, 2] + tapes[:, 3]), 3)
    else:
        px = np.tile(0.5 * (tapes[:, 0] + tapes[:, 1]), 3); py = np.concatenate([tapes[:, 2], 0.5 * (tapes[:, 2] + tapes[:, 3]), tapes[:, 3]])
    own = np.tile(np.arange(n), 3)
    Bx = np.zeros(len(px)); By = np.zeros(len(px))
    for k, (x1, x2, y1, y2) in enumerate(tapes):
        bx, by = rect_field(px, py, x1, x2, y1, y2, 1.0 / L)        # 1 A in the coil is 1/L A in each layer
        m = own != k
        Bx[m] += bx[m]; By[m] += by[m]
    Bn, Bt = (By, Bx) if normal == "y" else (Bx, By)
    return own, Bn, Bt, n, bw, bh


def tape_fields(kind, N, L, w, I_coil):
    """per tape at coil current I_coil: the field magnitude, the angle from the normal and the perpendicular component,
    each the mean of the tape's three samples (the strip model takes a uniform applied field), and the worst sample's values"""
    own, Bn1, Bt1, n, bw, bh = sample_fields(kind, N, L, w)
    Bn = np.abs(Bn1) * I_coil; Bt = np.abs(Bt1) * I_coil
    B = np.hypot(Bn, Bt); theta = np.degrees(np.arctan2(Bt, Bn))
    Bm = np.array([B[own == k].mean() for k in range(n)]); Bnm = np.array([Bn[own == k].mean() for k in range(n)])
    thm = np.degrees(np.arccos(np.clip(Bnm / np.maximum(Bm, 1e-30), 0, 1)))
    return dict(B=Bm, Bperp=Bnm, theta=thm, B_samples=B, theta_samples=theta, own=own, n=n, width=bw, height=bh)


def worst_tape(kind, N, L, w, I_coil):
    """the lowest measured critical current over every sample of every tape at coil current I_coil (the tape-limit rule)"""
    tf = tape_fields(kind, N, L, w, I_coil)
    ic = ic_tape(tf["B_samples"], tf["theta_samples"], w)
    k = int(ic.argmin())
    return dict(ic=float(ic[k]), B=float(tf["B_samples"][k]), theta=float(tf["theta_samples"][k]), Bmax=float(tf["B_samples"].max()),
                width=tf["width"], height=tf["height"], n_samples=len(ic), ic_perp=float(ic_tape(np.array([tf["B_samples"][k]]), np.array([0.0]), w)[0]))


def coil_limit(kind, N, L, w, share=0.5, I0=163.0):
    """the coil current at which the worst tape reaches `share` of its measured limit (the field grows with the current)"""
    I = I0
    for _ in range(60):
        I_new = share * worst_tape(kind, N, L, w, I)["ic"] * L
        if abs(I_new - I) < 1e-3:
            break
        I = 0.5 * (I + I_new)
    return I, worst_tape(kind, N, L, w, I)


# ---------------------------------------------------------------- the losses
def f_brandt(h):
    """the Brandt-Indenbom loss function (2/h) ln cosh h - tanh h, with its small-h series below 1e-3 to avoid cancellation"""
    h = np.asarray(h, float)
    small = h < 1e-3
    out = np.empty_like(h)
    hs = np.where(small, 1.0, h)
    out[~small] = (2.0 / hs[~small]) * (np.logaddexp(hs[~small], -hs[~small]) - math.log(2.0)) - np.tanh(hs[~small])
    out[small] = h[small]**3 / 6.0
    return out


def q_strip_perp(Ic, w, B_perp):
    """Brandt-Indenbom loss of a thin strip per metre per cycle (J/m): Ic (A) its critical current, w (m) its width, B_perp (T)
    the amplitude of the field through its face"""
    Jcd = Ic / w; Hc = Jcd / math.pi; Ha = B_perp / MU0
    h = Ha / np.maximum(Hc, 1e-30)
    return 4 * MU0 * Jcd * (w / 2)**2 * Ha * f_brandt(h)


def q_norris_strip(I, Ic):
    """Norris's self-field transport loss of a thin strip per metre per cycle (J/m)"""
    i = np.clip(np.asarray(I, float) / np.maximum(Ic, 1e-30), -0.999, 0.999)
    return (MU0 * Ic**2 / math.pi) * ((1 - i) * np.log(1 - i) + (1 + i) * np.log(1 + i) - i**2)


def winding_loss(kind, N, L, w, I_coil):
    """the loss of the whole winding at coil current I_coil, per metre of turn and per cycle (J/m), averaged over the N turns:
    (perpendicular-field strip loss summed over all tapes, Norris transport loss summed over all tapes, and per-tape arrays)"""
    tf = tape_fields(kind, N, L, w, I_coil)
    Ic = ic_tape(tf["B"], tf["theta"], w)
    qp = q_strip_perp(Ic, w, tf["Bperp"])
    qn = q_norris_strip(I_coil / L, Ic)
    return float(qp.sum() / N), float(qn.sum() / N), dict(Ic=Ic, qp=qp, qn=qn, Bperp=tf["Bperp"], B=tf["B"], theta=tf["theta"], own=tf["own"], n=tf["n"])


def loss_table(kind, N, L, w, I_max, n_pts=48):
    """Q(I) for the winding: coil currents from 0 to I_max and the loss per metre of turn per cycle, perpendicular plus Norris
    (the upper-bound sum the simulation uses) and the two parts; built once per configuration"""
    Is = np.linspace(0.0, I_max, n_pts)
    qp = np.zeros(n_pts); qn = np.zeros(n_pts)
    for k, I in enumerate(Is[1:], 1):
        qp[k], qn[k], _ = winding_loss(kind, N, L, w, I)
    return dict(I=Is, q=qp + qn, q_perp=qp, q_norris=qn, kind=kind, N=N, L=L, w=w)


def q_from_table(table, I):
    return float(np.interp(abs(I), table["I"], table["q"], right=table["q"][-1] * (abs(I) / table["I"][-1])**2))


def tape_mass_per_lim(N, L, w, w_coil=2.0, tau_p=100.0, phases=3, pitches=3):
    return 2 * (w_coil + tau_p) * N * L * w * TAPE_T * TAPE_RHO * phases * pitches


# ---------------------------------------------------------------- the checks (item 2 of the 2026-10-27 task)
def checks():
    print("CHECKS OF THE BRANDT-INDENBOM STRIP LOSS AS IMPLEMENTED")
    w, Ic = 3e-3, 223.0                                   # the 3 mm tape with no field, 70 K (Robinson)
    Jcd = Ic / w; Hc = Jcd / math.pi
    print(f"  3 mm tape, Ic {Ic:.0f} A: Jc d = Ic / w = {Jcd:.3e} A/m, H_c = Jc d / pi = {Hc:.3e} A/m, mu0 H_c = {MU0*Hc*1e3:.1f} mT (the same for every width)")
    for h in (0.01, 0.1, 1.0, 8.0, 30.0):
        B = MU0 * Hc * h
        q = float(q_strip_perp(Ic, w, B)); sat = w * Ic * B
        print(f"  h = {h:5.2f} (B_perp {B*1e3:7.2f} mT): Q = {q:.4e} J/m per cycle; f(h) = {float(f_brandt(h)):.5f}; Q / (w Ic B_a) = {q/sat:.5f}; Q / (w Ic B_a h^3/6) = {q/(sat*h**3/6):.4f}")
    print("  limits: f(h) -> 1 for h >> 1 (Q -> w Ic B_a, the loss proportional to width x Ic x field); f(h) -> h^3 / 6 for h << 1 (Q goes as B_a^4).")
    print(f"  at h = 1: f(1) = 2 ln cosh 1 - tanh 1 = {2*math.log(math.cosh(1))-math.tanh(1):.5f}; at h = 8 (the stack-end tapes of the pancake, 0.24 T): f = {float(f_brandt(8)):.4f}.")
    # the small-field series against the closed form just above the switch
    print(f"  continuity at the series switch: f(0.00099) = {float(f_brandt(0.00099)):.4e} (series), f(0.0011) = {float(f_brandt(0.0011)):.4e} (closed form), ratio to h^3/6: {float(f_brandt(0.0011))/(0.0011**3/6):.5f}")
    # Norris at zero applied field
    for i in (0.3, 0.5, 0.9):
        print(f"  Norris at zero applied field, I/Ic = {i}: Q_perp = {float(q_strip_perp(Ic, w, 0.0)):.1e} J/m, Q_norris = {float(q_norris_strip(i*Ic, Ic)):.4e} J/m: the combined loss is Norris's alone, by construction.")
    # the published form (Grilli et al. 2014, the preprint's page 16, Eq. 40): per unit volume Q_v = mu0 w Jc H_m [(2/x) ln cosh x - tanh x],
    # x = H_m / H_c, H_c = Jc d / pi, with Q_v -> mu0 w Jc H_m for H_m >> H_c; per unit length that is w d Q_v = 4 mu0 Jc d a^2 H_m f(x), the form
    # implemented.  Its Eq. 41: the peak of the loss over H_m^2 (the chi'' peak) is at H_peak = 0.7845 Jc d, that is x = 0.7845 pi = 2.465.
    x = np.linspace(0.5, 6.0, 55001); g = f_brandt(x) / x; k = int(g.argmax())
    print(f"  published form (Grilli et al. 2014 Eq. 40): Q_perp / (w d) -> mu0 w Jc H_m for h >> 1, here {float(q_strip_perp(Ic, w, MU0*Hc*30))/(w*TAPE_T):.4e} against {MU0*w*(Jcd/TAPE_T)*(Hc*30):.4e} J/m3 at h = 30 (f = {float(f_brandt(30)):.4f});")
    print(f"  published point (Eq. 41): the loss over H_m^2 peaks at H_peak = 0.7845 Jc d; here at x = {float(x[k]):.4f}, that is H_peak / (Jc d) = {float(x[k])/math.pi:.4f}.")


if __name__ == "__main__":
    checks()
