#!/usr/bin/env python3
"""
Round 10, item 3: the mass driver's stator (Chapters 9, 10 and 14).
Book: Orbital Ring Engineering, Volume III: Table 9.6, Section 10.1 (Equations 10.1 to 10.3, Table 10.2), Section 14.5.
Brief: Vol_III/Claude_Code_instructions_round10_2026-11-17.md, Section 3.
Run: python tools/round10_stator.py          (about two minutes; writes runs/round10_logs/round10_stator.txt and .json)

What it does:
  1. The winding.  States the coil that lsm_physics.calc_EMF assumes (a rectangle one pole pitch long and w_coil wide, N turns), and checks
     by a direct sum over its conductors (flux linkage for the back-EMF, the Lorentz force for the thrust, and their energy balance) what
     that coil gives.  Compares with the three things lsm_physics.py computes: calc_EMF, calc_thrust on calc_B_stator, and the tape length
     of calc_hysteresis_power.
  2. The tape limit.  lim_hysteresis.coil_limit on the flat winding (N turns of L layers of 12 mm tape, Robinson 70 K data) at the 50
     percent rule, in the winding's own field, then with the coil's return run and with the sled's field added, and the layers that
     keep the book's coil current.  The launch again for each case (the same time stepping as lsm_simulation.py).
  3. The hysteresis on the Brandt and Indenbom strip model (lim_hysteresis.q_strip_perp, q_norris_strip) for the full-pitch coil: per
     500 m site at 30 km/s, per pass, per launch, as a share of the launch's energy; the tape mass per meter of ring.
It changes no file but its own two logs.  It uses the round 10 additions at the end of lsm_physics.py (calc_turn_length_full_pitch,
calc_thrust_full_pitch, calc_power_limit_full_pitch), which change nothing in that module's existing functions.
"""
import contextlib
import io
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the repository root
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # round 13: this study was run on the old 130 m stator of lsm_config.py (its launch and its
                                                   # constants); the switch keeps it reproducing now that the 6 m stator is the default
import matplotlib
matplotlib.use("Agg")
import numpy as np

import lsm_config as cfg
import lsm_physics as phys
import lsm_simulation as sim
import lim_hysteresis as lh
import ring_altitude as ring

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round10_logs")
SITE = 500.0                                   # m, the spacing of the ring's sites
SHARE = 0.5                                    # the 50 percent rule of the bearing and the deployment LIMs
KIND = "flat"
SPARE_W = 11.9e3                               # W per site, the plant's spare capacity in a launch window (Chapter 14, Section 14.5; a book figure, not computed here)
COP_PENALTY = 94.0                             # W of compressor work per W lifted from 70 K (Chapter 14's figure, used only for one line)
_lines = []


def out(s=""):
    print(s)
    _lines.append(s)


# ------------------------------------------------------------------------------------------------ part 1: the winding, checked directly
def direct_check(N, I, B0, w, tau, v, delta, n_steps=2400, n_quad=4000):
    """One three-phase set of full-pitch coils (phase k's coil spans x_k .. x_k + tau, x_k = 2 k tau / 3) in the sled's field
    B(x, t) = B0 sin(pi (x - v t) / tau).  Returns the peak back-EMF of a coil (from the flux linkage, by numerical differences),
    the mean and the ripple of the force on the three coils (Lorentz force on the two sides of each), and the mean electrical power
    going into the back-EMFs (sum of e i)."""
    T = 2 * tau / v
    t = np.linspace(0.0, T, n_steps, endpoint=False); dt = t[1] - t[0]
    xq = (np.arange(n_quad) + 0.5) / n_quad * tau
    F = np.zeros(n_steps); P = np.zeros(n_steps); e_peak = 0.0
    for k in range(3):
        xk = 2.0 * k * tau / 3.0
        B = lambda x, tt: B0 * np.sin(math.pi * (x - v * tt) / tau)
        lam = np.array([N * w * np.sum(B(xk + xq, tt)) * (tau / n_quad) for tt in t])          # flux linkage, Wb-turns
        e = -(np.roll(lam, -1) - np.roll(lam, 1)) / (2 * dt)                                   # back-EMF of the coil
        e_peak = max(e_peak, float(np.abs(e).max()))
        # the current that gives the load angle delta: in phase with the back-EMF (which goes as sin ph) at delta = 90 degrees
        ph = math.pi * (xk - v * t) / tau
        i = I * np.cos(ph - delta)
        F += N * i * w * (B(xk, t) - B(xk + tau, t))                                           # the two sides, opposite currents
        P += -e * i                                                                            # the supply works against the back-EMF
    return dict(e_peak=e_peak, F_mean=float(F.mean()), F_ripple=float(F.max() - F.min()), P_mean=float(P.mean()), Fv=float(F.mean() * v))


# ------------------------------------------------------------------------------------------------ part 2: the tape of the flat winding
def flat_samples(N, L, w):
    """the sample points of lim_hysteresis.sample_fields for the flat layout, in the same order (x, y, tape index)"""
    rects, normal, bw, bh = lh.layout(KIND, N, L, w)
    tapes = np.array([(x1, x2, y1 + j * lh.TAPE_T, y1 + (j + 1) * lh.TAPE_T) for (x1, x2, y1, y2) in rects for j in range(L)])
    px = np.concatenate([tapes[:, 0], 0.5 * (tapes[:, 0] + tapes[:, 1]), tapes[:, 1]]); py = np.tile(0.5 * (tapes[:, 2] + tapes[:, 3]), 3)
    return px, py, rects, bw


def return_run_field(N, L, w, D):
    """the field per ampere of coil current, at the sample points of one run of the coil, of the coil's other run: the same flat sheet in
    the same plane, its center D away, carrying the opposite current.  Returns the signed normal and tangential components (T per A)."""
    px, py, rects, bw = flat_samples(N, L, w)
    Bx = np.zeros(len(px)); By = np.zeros(len(px))
    for (x1, x2, y1, y2) in rects:
        bx, by = lh.rect_field(px, py, x1 - D, x2 - D, y1, y2, -1.0)
        Bx += bx; By += by
    return By, Bx                                                                               # flat: the normal is y


def state(N, L, w, I, extra=None, B_sled=0.0):
    """The flat winding at coil current I: the worst tape's measured limit (every sample) and the strip-model loss per meter of turn per
    cycle, as lim_hysteresis.worst_tape and winding_loss do it, with two optional additions to the field at each tape:
    extra = (Bn, Bt) per ampere, signed (the coil's return run), and B_sled, an amplitude through the tape's face that is added to the
    amplitude of the winding's own normal field (the two are at the same frequency; adding the amplitudes is the upper bound)."""
    own, Bn1, Bt1, n, bw, bh = lh.sample_fields(KIND, N, L, w)
    Bn = Bn1 * I; Bt = Bt1 * I
    if extra is not None:
        Bn = Bn + extra[0] * I; Bt = Bt + extra[1] * I
    Bn = np.abs(Bn) + B_sled; Bt = np.abs(Bt)
    B = np.hypot(Bn, Bt); theta = np.degrees(np.arctan2(Bt, Bn))
    ic_s = lh.ic_tape(B, theta, w); k = int(ic_s.argmin())
    Bm = np.array([B[own == j].mean() for j in range(n)]); Bnm = np.array([Bn[own == j].mean() for j in range(n)])
    thm = np.degrees(np.arccos(np.clip(Bnm / np.maximum(Bm, 1e-30), 0, 1)))
    Ic = lh.ic_tape(Bm, thm, w)
    qp = lh.q_strip_perp(Ic, w, Bnm); qn = lh.q_norris_strip(I / L, Ic)
    return dict(ic=float(ic_s[k]), B=float(B[k]), theta=float(theta[k]), Bmax=float(B.max()), share=I / L / float(ic_s[k]),
                q_perp=float(qp.sum() / N), q_norris=float(qn.sum() / N), q=float((qp.sum() + qn.sum()) / N), width=bw)


def limit(N, L, w, share=SHARE, extra=None, B_sled=0.0):
    """the coil current at which the worst tape carries `share` of its measured limit (bisection; the limit falls as the current rises)"""
    lo, hi = 1.0, 2000.0 * L
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if share * state(N, L, w, mid, extra, B_sled)["ic"] * L > mid:
            lo = mid
        else:
            hi = mid
    I = 0.5 * (lo + hi)
    return I, state(N, L, w, I, extra, B_sled)


# ------------------------------------------------------------------------------------------------ the launch
def launch(F_peak, P_cap, m, v_target, dt=None, marks=()):
    """The launch as lsm_simulation.py steps it (x += v dt, then v += a dt): the thrust is F_peak until F_peak v reaches P_cap (the
    back-EMF limit) and P_cap / v after that; the ring's HVDC limit of lsm_config applies on top.  No cap on the acceleration is set
    here (lsm_config.A_MAX_G is 1 g; no case in this script passes it but the empty sled, which Section 10.1 gives at 1.6 g).
    Returns the time, the distance, the cycles of the supply (integral of f dt), and the time and distance at each speed in marks."""
    dt = cfg.DT if dt is None else dt
    if cfg.P_HVDC_MAX is not None:
        P_cap = min(P_cap, cfg.P_HVDC_MAX)
    t = 0.0; v = 0.0; x = 0.0; cyc = 0.0; hit = {}; todo = sorted(marks); vs = []
    while v < v_target:
        F = F_peak if F_peak * v <= P_cap else P_cap / v
        vs.append(v)
        cyc += v / (2 * cfg.TAU_P) * dt
        x += v * dt; v += F / m * dt; t += dt
        while todo and v >= todo[0]:
            hit[todo.pop(0)] = (t, x)
    return dict(t=t, x=x, cycles=cyc, marks=hit, v=np.array(vs), dt=dt)


def launch_rows(label, F_peak, P_cap, res_list, dt=None):
    m_full = cfg.M_TOTAL; m_empty = cfg.M_SLED_HARDWARE
    r = launch(F_peak, P_cap, m_full, cfg.V_LAUNCH, dt)
    v1 = cfg.V_ORBIT - ring.V_GROUND_SYNC; v2 = 9800.0 - ring.V_GROUND_SYNC      # speeds over the casing for 7,504 and 9,800 m/s in space (eastward launch)
    e = launch(F_peak, P_cap, m_empty, v2, dt, marks=(v1, v2))
    row = dict(label=label, F_peak_N=F_peak, P_cap_W=min(P_cap, cfg.P_HVDC_MAX), a_full_g=F_peak / m_full / cfg.G_0, t_min=r["t"] / 60, x_km=r["x"] / 1e3, laps=r["x"] / cfg.L_RING,
               cycles=r["cycles"], a_empty_g=F_peak / m_empty / cfg.G_0, t_orbit_min=e["marks"][v1][0] / 60, x_orbit_km=e["marks"][v1][1] / 1e3,
               t_98_min=e["marks"][v2][0] / 60, x_98_km=e["marks"][v2][1] / 1e3, v_cross=min(P_cap, cfg.P_HVDC_MAX) / F_peak)
    res_list.append(row)
    row["_v"] = r["v"]; row["_dt"] = r["dt"]
    return row


def fmt_launch(r):
    return (f"     {r['label']:<58} {r['F_peak_N']/1e6:8.2f} {r['P_cap_W']/1e9:8.2f} {r['a_full_g']:7.3f} {r['t_min']:9.1f} {r['t_min']/60:7.2f} {r['x_km']:11,.0f} {r['laps']:7.2f} |"
            f" {r['a_empty_g']:6.3f} {r['t_orbit_min']:8.1f} {r['x_orbit_km']:9,.0f} {r['t_98_min']:8.1f} {r['x_98_km']:9,.0f}")


def main():
    os.makedirs(LOGS, exist_ok=True)
    res = {}
    N = cfg.N_STATOR; Lyr = cfg.N_HTS_LAYERS; w = cfg.TAPE_WIDTH; tau = cfg.TAU_P; wc = cfg.W_COIL; B0 = cfg.B_SLED_NOMINAL
    I_book = cfg.I_PEAK; dmax = cfg.DELTA_MAX; Vlim = cfg.V_COIL_LIMIT
    v_cross = Vlim / (2 * N * B0 * wc)
    n_pp = cfg.N_UNITS * cfg.N_POLES_PER_UNIT * cfg.N_LIM_SIDES              # pole pitches of stator facing the sled, both sides
    n_code = cfg.LIM_PHASES * cfg.N_LIM_SIDES * cfg.N_UNITS                  # the coils calc_hysteresis_power counts
    n_full = cfg.LIM_PHASES * n_pp                                           # a coil of every phase starting in every pole pitch (two-layer winding)
    n_sites_sled = cfg.L_SLED / SITE
    t_pass = cfg.L_SLED / cfg.V_LAUNCH
    f_end = cfg.V_LAUNCH / (2 * tau)
    KE = 0.5 * cfg.M_TOTAL * cfg.V_LAUNCH**2
    out("=" * 150)
    out("  ROUND 10, ITEM 3: THE MASS DRIVER'S STATOR")
    out("=" * 150)
    out(f"  lsm_config.py: {N} turns, {Lyr} layers of {w*1e3:.0f} mm tape, pole pitch {tau:.0f} m, coil width {wc:.1f} m, gap {cfg.G_GAP*1e3:.0f} mm, {Vlim/1e3:.0f} kV, sled field {B0:.2f} T,")
    out(f"  I_PEAK {I_book:.1f} A, load angle limit {math.degrees(dmax):.0f} deg, sled {cfg.L_SLED/1e3:.0f} km = {cfg.N_UNITS} units of {cfg.N_POLES_PER_UNIT} pole pitches + {cfg.L_GAP:.0f} m, both sides: {n_pp} pole pitches,")
    out(f"  active area {cfg.A_ACTIVE_TOTAL:,.0f} m2, launch mass {cfg.M_TOTAL/1e3:,.0f} t, sled alone {cfg.M_SLED_HARDWARE/1e3:,.0f} t, ring {cfg.L_RING/1e3:,.1f} km.")

    # ================================================================================================ 1. the winding
    out()
    out("  1. THE WINDING")
    out("  1.1 What each function of lsm_physics.py assumes about the coil")
    turn_emf = phys.calc_turn_length_full_pitch(tau, wc)
    depth = N * 0.15e-3; turn_hyst = 2.0 * (wc + depth)
    B_code = phys.calc_B_stator(I_book, N, wc, cfg.G_GAP)
    F_code = phys.calc_thrust(B_code, B0, dmax, cfg.A_ACTIVE_TOTAL)
    E_code = phys.calc_EMF(N, v_cross, B0, wc)
    out(f"     calc_EMF (Equation 9.3): 'each stator coil spans one pole pitch with N turns and width w_coil'; EMF = 2 N v B w = {E_code/1e3:.1f} kV at {v_cross:,.0f} m/s.")
    out(f"        That coil is a rectangle {tau:.0f} m x {wc:.0f} m.  One turn is 2 x ({tau:.0f} + {wc:.0f}) = {turn_emf:.0f} m long.")
    out(f"     calc_hysteresis_power: turn = 2 x (w_coil + N x 0.15 mm) = {turn_hyst:.3f} m; coils under the sled = 3 phases x 2 sides x {cfg.N_UNITS} units = {n_code}.")
    out(f"        That is the two {wc:.0f} m sides of the same coil without the two {tau:.0f} m runs that join them: {turn_emf/turn_hyst:.1f} times too little tape.")
    out(f"     calc_B_stator (Equation 10.1): B = (2/pi) (mu0 N I / w_coil) atan(w_coil / 2g) = {B_code*1e3:.2f} mT at {I_book:.1f} A; calc_thrust (Equation 10.2) spreads")
    out(f"        that field as a sine wave over the whole active area: {F_code/1e6:.2f} MN at {math.degrees(dmax):.0f} deg.  mu0 N I / w_coil is the field beside a sheet that carries N I = {N*I_book/1e3:.1f} kA")
    out(f"        in every {wc:.0f} m of its length ({N*I_book/wc/1e3:.1f} kA/m).  The {tau:.0f} m coil has N I in every {tau/3:.1f} m (three phases), {3*N*I_book/tau/1e3:.2f} kA/m.")

    out()
    out("  1.2 The full-pitch coil checked directly (function direct_check: flux linkage and Lorentz force summed over the conductors, one three-phase set)")
    out(f"     {'speed m/s':>10} {'load angle':>10} {'EMF peak kV':>12} {'2NvBw kV':>9} {'force kN, 3 coils':>18} {'ripple kN':>10} {'3 N I B w sin d kN':>19} {'F v MW':>9} {'-sum e i MW':>11}")
    chk = []
    for v, d in ((v_cross, dmax), (v_cross, math.pi / 2), (1000.0, dmax)):
        c = direct_check(N, I_book, B0, wc, tau, v, d)
        ana = 3 * N * I_book * B0 * wc * math.sin(d)
        c.update(v=v, delta_deg=math.degrees(d), emf_formula=phys.calc_EMF(N, v, B0, wc), F_formula=ana); chk.append(c)
        out(f"     {v:10,.0f} {math.degrees(d):10.0f} {c['e_peak']/1e3:12.3f} {c['emf_formula']/1e3:9.3f} {c['F_mean']/1e3:18.3f} {c['F_ripple']/1e3:10.4f} {ana/1e3:19.3f} {c['Fv']/1e6:9.2f} {c['P_mean']/1e6:11.2f}")
    out(f"     The back-EMF of calc_EMF is confirmed.  The force on a coil is N I B w sin(delta) = {N*I_book*B0*wc*math.sin(dmax)/1e3:.2f} kN at {I_book:.1f} A and {math.degrees(dmax):.0f} deg, and force x speed")
    out(f"     equals the electrical power into the back-EMFs, (1/2) EMF I sin(delta) per coil.  At the {Vlim/1e3:.0f} kV limit a coil can pass (1/2) x {Vlim/1e3:.0f} kV x {I_book:.1f} A x sin {math.degrees(dmax):.0f} = {0.5*Vlim*I_book*math.sin(dmax)/1e6:.1f} MW.")
    res["check"] = chk

    out()
    out("  1.3 What the full-pitch coil gives against Equations 10.1 to 10.3 (lsm_physics.calc_thrust_full_pitch, calc_power_limit_full_pitch)")
    out(f"     {'coils facing the sled':<66} {'coils':>6} {'thrust MN':>10} {'power GW':>9} {'book thrust /':>13} {'equiv. B mT':>11}")
    wind = []
    n_need = F_code / (N * I_book * B0 * wc * math.sin(dmax))
    for lab, nc in ((f"as calc_hysteresis_power counts: one per phase per unit, per side", n_code),
                    (f"one per phase per two pole pitches (single layer, fully wound)", n_full // 2),
                    (f"one per phase per pole pitch (two layers, fully wound)", n_full),
                    (f"the count that Equation 10.2's {F_code/1e6:.1f} MN needs at {I_book:.0f} A", n_need)):
        Fw = phys.calc_thrust_full_pitch(N, I_book, B0, wc, dmax, nc); Pw = phys.calc_power_limit_full_pitch(Vlim, I_book, dmax, nc)
        B_eq = 2 * cfg.MU0 * Fw / (B0 * math.sin(dmax) * cfg.A_ACTIVE_TOTAL)
        wind.append(dict(label=lab, coils=nc, F_N=Fw, P_W=Pw, ratio=F_code / Fw, B_equiv_T=B_eq, coils_per_phase_per_pole_pitch_per_side=nc / 3 / n_pp))
        out(f"     {lab:<66} {nc:6.0f} {Fw/1e6:10.3f} {Pw/1e9:9.2f} {F_code/Fw:13.1f} {B_eq*1e3:11.2f}")
    out(f"     'equiv. B' is the stator field that Equation 10.2 would need to give the same thrust; Equation 10.1 gives {B_code*1e3:.1f} mT.")
    out(f"     The book's {F_code/1e6:.1f} MN and {F_code*v_cross/1e9:.0f} GW need {n_need:,.0f} coils of {N} turns at {I_book:.0f} A and {Vlim/1e3:.0f} kV facing the sled, {n_need/3/n_pp:.1f} per phase per pole pitch per side,")
    out(f"     or {n_need*N*I_book/1e6:,.0f} MA-turns in all against the {n_code*N*I_book/1e6:.1f} MA-turns of the {n_code} coils that the hysteresis function counts.")
    out(f"     Check on the book's power by its own voltage and current: {F_code*v_cross/1e9:.1f} GW / ((1/2) x {Vlim/1e3:.0f} kV x {I_book:.0f} A x sin 60) = {F_code*v_cross/(0.5*Vlim*I_book*math.sin(dmax)):,.0f} coils.")
    res["winding"] = dict(turn_full_pitch_m=turn_emf, turn_hysteresis_fn_m=turn_hyst, ratio=turn_emf / turn_hyst, B_stator_code_T=B_code, F_code_N=F_code, P_code_W=F_code * v_cross, v_cross=v_cross,
                          sheet_code_kA_per_m=N * I_book / wc / 1e3, sheet_winding_kA_per_m=3 * N * I_book / tau / 1e3, n_coils_code=n_code, n_coils_full=n_full, n_coils_needed=n_need, cases=wind)

    # rough inductive voltage of the coil (not in lsm_simulation.py)
    bw = lh.winding_width(KIND, N, Lyr, w)
    gmd = 0.2235 * (bw + Lyr * lh.TAPE_T)
    L_coil = cfg.MU0 / math.pi * N**2 * (tau * math.log(wc / gmd) + wc * math.log(tau / gmd))
    out()
    out(f"  1.4 For information [UNVERIFIED, rough]: the coil's own inductance, two-wire formula with the flat sheet's mean distance ({gmd*1e3:.0f} mm): {L_coil*1e3:.0f} mH;")
    out(f"     its reactive voltage 2 pi f L I at {I_book:.0f} A: {2*math.pi*v_cross/(2*tau)*L_coil*I_book/1e3:.0f} kV at {v_cross/1e3:.0f} km/s ({v_cross/(2*tau):.1f} Hz), {2*math.pi*f_end*L_coil*I_book/1e3:.0f} kV at 30 km/s ({f_end:.1f} Hz).")
    out(f"     lsm_simulation.py holds only the back-EMF to {Vlim/1e3:.0f} kV.  The mutual inductance of the other phases and of the sled's coils is not in this figure.")
    res["inductance_rough"] = dict(L_H=L_coil, V_at_cross=2 * math.pi * v_cross / (2 * tau) * L_coil * I_book, V_at_end=2 * math.pi * f_end * L_coil * I_book, flag="UNVERIFIED")

    # ================================================================================================ 2. the tape limit
    out()
    out("  2. THE TAPE LIMIT (70 K, Robinson data, flat winding, lim_hysteresis.py)")
    I_lh, wt_lh = lh.coil_limit(KIND, N, Lyr, w, share=SHARE, I0=1000.0)
    I50, st50 = limit(N, Lyr, w)
    s_book = state(N, Lyr, w, I_book)
    wl = lh.winding_loss(KIND, N, Lyr, w, I_book); wb = lh.worst_tape(KIND, N, Lyr, w, I_book)
    assert abs(s_book["q"] - (wl[0] + wl[1])) < 1e-9 * s_book["q"] and abs(s_book["ic"] - wb["ic"]) < 1e-9 * wb["ic"]
    out(f"  2.1 Check of this script's state() against lim_hysteresis at {I_book:.1f} A: worst tape {s_book['ic']:.1f} A (lim_hysteresis.worst_tape {wb['ic']:.1f}), loss {s_book['q']*1e3:.3f} mJ/m of turn")
    out(f"      (lim_hysteresis.winding_loss {(wl[0]+wl[1])*1e3:.3f}); lim_hysteresis.coil_limit at 50 percent: {I_lh:.1f} A; this script's bisection: {I50:.1f} A.")
    extra = return_run_field(N, Lyr, w, wc)
    B_sled_end = Vlim / (2 * N * cfg.V_LAUNCH * wc)
    cases = [("own field of one run of the coil (lim_hysteresis as it stands)", None, 0.0, ""),
             (f"plus the coil's return run {wc:.0f} m away in the same plane", extra, 0.0, ""),
             (f"plus the sled's {B0*1e3:.0f} mT through the face (below {v_cross/1e3:.0f} km/s)", extra, B0, "[UNVERIFIED: amplitudes added, an upper bound]")]
    out()
    out(f"  2.2 The flat winding of {N} turns, {Lyr} layers of {w*1e3:.0f} mm tape, sheet {bw*1e3:.0f} mm wide: the coil current at {SHARE*100:.0f} percent of the worst tape's limit")
    out(f"      {'field at the tapes':<66} {'coil A':>8} {'per tape A':>10} {'worst Ic A':>10} {'at T':>6} {'deg':>5} {'x book I':>8} | {'at the book 2,107 A: Ic A':>25} {'share %':>8}")
    lim_rows = []
    for lab, ex, Bs, flag in cases:
        I_c, st = limit(N, Lyr, w, extra=ex, B_sled=Bs)
        sb = state(N, Lyr, w, I_book, ex, Bs)
        lim_rows.append(dict(label=lab, I_50=I_c, I_tape=I_c / Lyr, ic_worst=st["ic"], B=st["B"], theta=st["theta"], ratio=I_c / I_book, ic_at_book=sb["ic"], share_at_book=sb["share"], flag=flag))
        out(f"      {lab:<66} {I_c:8.1f} {I_c/Lyr:10.1f} {st['ic']:10.1f} {st['B']:6.3f} {st['theta']:5.1f} {I_c/I_book:8.3f} | {sb['ic']:25.1f} {sb['share']*100:8.1f}  {flag}")
    res["limit"] = lim_rows
    I_A = lim_rows[0]["I_50"]; I_B = lim_rows[1]["I_50"]; I_C = lim_rows[2]["I_50"]

    out()
    out(f"  2.3 The layers of {w*1e3:.0f} mm tape, same {N} turns, and the coil current at {SHARE*100:.0f} percent (own field and return run; last columns with the sled's {B0*1e3:.0f} mT added)")
    out(f"      {'layers':>6} {'coil A':>9} {'per tape A':>10} {'worst Ic A':>10} {'x book I':>8} | {'coil A':>9} {'x book I':>8}")
    lay = []; L_keep = None; L_keep_s = None
    for L_ in range(5, 17):
        ex = return_run_field(N, L_, w, wc)
        Ia, sa = limit(N, L_, w, extra=ex); Ib, sb = limit(N, L_, w, extra=ex, B_sled=B0)
        lay.append(dict(layers=L_, I_50=Ia, ic=sa["ic"], I_50_sled=Ib, ic_sled=sb["ic"]))
        if L_keep is None and Ia >= I_book: L_keep = L_
        if L_keep_s is None and Ib >= I_book: L_keep_s = L_
        out(f"      {L_:6d} {Ia:9.1f} {Ia/L_:10.1f} {sa['ic']:10.1f} {Ia/I_book:8.3f} | {Ib:9.1f} {Ib/I_book:8.3f}")
    out(f"      Layers that keep the book's {I_book:.0f} A at {SHARE*100:.0f} percent: {L_keep} in the winding's own field, {L_keep_s} with the sled's field added (upper bound).")
    res["layers"] = lay; res["layers_keep"] = L_keep; res["layers_keep_sled"] = L_keep_s

    # ---- the launch again
    out()
    out("  2.4 The launch again (function launch: the time stepping of lsm_simulation.py, dt 0.1 s; thrust F until F v reaches the back-EMF limit, then that power)")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok = sim.run_lsm_simulation(quick_mode=False)
    t_sim = (sim.data_t[-1] + cfg.DT) / 60; x_sim = (sim.data_x[-1] + sim.data_v[-1] * cfg.DT) / 1e3
    F_sim = max(sim.data_F); P_sim = max(sim.data_P_mech)
    out(f"      lsm_simulation.run_lsm_simulation() as it stands: {t_sim:.1f} minutes, {x_sim:,.1f} km, {x_sim*1e3/cfg.L_RING:.2f} laps, peak thrust {F_sim/1e6:.2f} MN, peak power {P_sim/1e9:.2f} GW.")
    res["sim"] = dict(ok=bool(ok), t_min=t_sim, x_km=x_sim, laps=x_sim * 1e3 / cfg.L_RING, F_peak_N=F_sim, P_peak_W=P_sim)
    out(f"      {'case':<58} {'thrust':>8} {'power':>8} {'accel':>7} {'time':>9} {'time':>7} {'distance':>11} {'laps':>7} | {'empty sled':>6} {'to 7,504 m/s':>18} {'to 9,800 m/s':>18}")
    out(f"      (empty sled: the two speeds are in space, as Section 10.1 has them; the casing moves at {ring.V_GROUND_SYNC:.0f} m/s, so the sled is at {cfg.V_ORBIT-ring.V_GROUND_SYNC:,.0f} and {9800-ring.V_GROUND_SYNC:,.0f} m/s over it)")
    out(f"      {'':<58} {'MN':>8} {'GW':>8} {'g':>7} {'min':>9} {'h':>7} {'km':>11} {'':>7} | {'g':>6} {'min':>8} {'km':>9} {'min':>8} {'km':>9}")
    L_rows = []
    fp = lambda I_, nc: (phys.calc_thrust_full_pitch(N, I_, B0, wc, dmax, nc), phys.calc_power_limit_full_pitch(Vlim, I_, dmax, nc))
    code = lambda I_: (phys.calc_thrust(phys.calc_B_stator(I_, N, wc, cfg.G_GAP), B0, dmax, cfg.A_ACTIVE_TOTAL),) * 1
    def add(label, F, P, dt=None):
        r = launch_rows(label, F, P, L_rows, dt); out(" " + fmt_launch(r)); return r
    out("      A. Equations 10.1 and 10.2 as lsm_physics.py has them (thrust in proportion to the coil current)")
    rA0 = add(f"A0 book: {I_book:.0f} A", code(I_book)[0], code(I_book)[0] * v_cross)
    rA1 = add(f"A1 {I_A:.0f} A (50 percent, own field)", code(I_A)[0], code(I_A)[0] * v_cross)
    rA2 = add(f"A2 {I_B:.0f} A (50 percent, with the return run)", code(I_B)[0], code(I_B)[0] * v_cross)
    rA3 = add(f"A3 {I_C:.0f} A (50 percent, with the sled's field)", code(I_C)[0], code(I_C)[0] * v_cross)
    out("      B. The full-pitch coil (thrust N I B w sin 60 per coil, power (1/2) 100 kV I sin 60 per coil); time step 1 s")
    rB0 = add(f"B0 {n_code} coils (the hysteresis count), book {I_book:.0f} A", *fp(I_book, n_code), dt=1.0)
    rB1 = add(f"B1 {n_code} coils, {I_B:.0f} A", *fp(I_B, n_code), dt=1.0)
    rB2 = add(f"B2 {n_full} coils (one per phase per pole pitch), {I_B:.0f} A", *fp(I_B, n_full), dt=1.0)
    rB3 = add(f"B3 {n_full} coils, {I_C:.0f} A (with the sled's field)", *fp(I_C, n_full), dt=1.0)
    n_keep_B = F_code / fp(I_B, 1)[0]; n_keep_C = F_code / fp(I_C, 1)[0]
    rB4 = add(f"B4 {n_keep_B:,.0f} coils at {I_B:.0f} A (keeps the book's thrust)", *fp(I_B, n_keep_B))
    out(f"      Row A0 against lsm_simulation: time {rA0['t_min']:.1f} against {t_sim:.1f} minutes, distance {rA0['x_km']:,.1f} against {x_sim:,.1f} km.")
    out(f"      Coils that keep {F_code/1e6:.1f} MN on the full-pitch coil: {n_keep_B:,.0f} at {I_B:.0f} A ({n_keep_B/3/n_pp:.1f} per phase per pole pitch per side), {n_keep_C:,.0f} at {I_C:.0f} A ({n_keep_C/3/n_pp:.1f}).")
    out(f"      Rows B use a 1 s step for the long launches (row B4, the book's thrust, uses 0.1 s and matches A0).")
    res["launch"] = [{k: v for k, v in r.items() if not k.startswith("_")} for r in L_rows]
    res["coils_keep"] = dict(at_I_B=n_keep_B, at_I_C=n_keep_C, per_phase_per_pp_per_side_B=n_keep_B / 3 / n_pp, per_phase_per_pp_per_side_C=n_keep_C / 3 / n_pp)

    # ================================================================================================ 3. the hysteresis
    out()
    out("  3. THE HYSTERESIS ON THE STRIP MODEL, FULL-PITCH COIL, FLAT WINDING")
    out(f"     Loss per meter of turn per cycle q (Brandt and Indenbom for the field through each tape's face at its measured 70 K limit, plus Norris's transport loss,")
    out(f"     an upper bound): function state.  Tape: {turn_emf:.0f} m per turn x {N} turns = {turn_emf*N/1e3:.1f} km of turn per coil, {turn_emf*N*Lyr/1e3:.1f} km of {w*1e3:.0f} mm tape.")
    out(f"     Sites of {SITE:.0f} m; the sled covers {n_sites_sled:.0f}; pass time at 30 km/s {t_pass:.4f} s; supply {f_end:.2f} Hz; kinetic energy of the launch {KE/1e15:.3f} PJ.")
    m_coil = turn_emf * N * Lyr * w * lh.TAPE_T * lh.TAPE_RHO
    q_tab = []
    out(f"     {'current and field':<64} {'q perp mJ/m':>11} {'Norris mJ/m':>11} {'q mJ/m':>9}")
    qcases = [("book", I_book, None, 0.0, f"{I_book:.0f} A, own field (round 9's row)"),
              ("A", I_A, None, 0.0, f"{I_A:.0f} A, own field"),
              ("B", I_B, extra, 0.0, f"{I_B:.0f} A, own field and return run"),
              ("B_sled_end", I_B, extra, B_sled_end, f"{I_B:.0f} A, plus the sled's {B_sled_end*1e3:.1f} mT at 30 km/s [upper bound]"),
              ("C_sled_end", I_C, extra, B_sled_end, f"{I_C:.0f} A, plus the sled's {B_sled_end*1e3:.1f} mT at 30 km/s [upper bound]"),
              ("C_sled_full", I_C, extra, B0, f"{I_C:.0f} A, plus the sled's {B0*1e3:.0f} mT below {v_cross/1e3:.0f} km/s [upper bound]"),
              ("idle_sled_end", 1e-6, None, B_sled_end, f"no coil current, the sled's {B_sled_end*1e3:.1f} mT alone"),
              ("idle_sled_full", 1e-6, None, B0, f"no coil current, the sled's {B0*1e3:.0f} mT alone")]
    q = {}
    for key, I_, ex, Bs, lab in qcases:
        s = state(N, Lyr, w, I_, ex, Bs); q[key] = s["q"]
        q_tab.append(dict(key=key, label=lab, I=I_, B_sled=Bs, q_perp=s["q_perp"], q_norris=s["q_norris"], q=s["q"]))
        out(f"     {lab:<64} {s['q_perp']*1e3:11.3f} {s['q_norris']*1e3:11.3f} {s['q']*1e3:9.3f}")
    res["q"] = q_tab

    def q_launch(I_, ex, vs, dt, with_sled, Lq=None):
        Lq = Lyr if Lq is None else Lq
        """integral of q f dt over the launch, per meter of turn (J/m); with_sled adds the sled's field B_sled(v) at each step"""
        f = vs / (2 * tau)
        if not with_sled:
            return float(state(N, Lq, w, I_, ex, 0.0)["q"] * np.sum(f) * dt)
        Bg = np.linspace(B_sled_end, B0, 24); qg = np.array([state(N, Lq, w, I_, ex, b)["q"] for b in Bg])
        Bs = np.minimum(B0, Vlim / (2 * N * np.maximum(vs, 1.0) * wc))
        return float(np.sum(np.interp(Bs, Bg, qg) * f) * dt)

    out()
    out(f"     {'case':<50} {'coils':>6} {'coil A':>7} {'turn km':>8} {'tape kg/m':>9} {'ring tape t':>11} | {'kW per site':>11} {'kJ per pass':>11} {'launch J':>11} {'share':>9} {'kJ/site/launch':>14} | {'with sled field: kW':>19} {'launch J':>11} {'share':>9}")
    out(f"     {'':<50} {'':>6} {'':>7} {'under sled':>8} {'of ring':>9} {'':>11} | {'at 30 km/s':>11} {'':>11} {'':>11} {'of KE':>9} {'mean':>14} | {'per site':>19} {'':>11} {'of KE':>9}")
    hyst = []
    Lk = L_keep; extra_k = return_run_field(N, Lk, w, wc); s_k = state(N, Lk, w, I_book, extra_k)
    hcases = [("legacy count (4.015 m per turn), book current", n_code, I_book, None, rA0, turn_hyst, Lyr),
              ("B0: full pitch, 150 coils, book current", n_code, I_book, None, rB0, turn_emf, Lyr),
              ("B1: full pitch, 150 coils, 50 percent", n_code, I_B, extra, rB1, turn_emf, Lyr),
              ("B2: full pitch, 450 coils, 50 percent", n_full, I_B, extra, rB2, turn_emf, Lyr),
              ("B4: full pitch, coils for 83.3 MN, 50 percent", n_keep_B, I_B, extra, rB4, turn_emf, Lyr),
              (f"B5: full pitch, {Lk} layers, book current, 83.3 MN", n_need, I_book, extra_k, rA0, turn_emf, Lk)]
    for lab, nc, I_, ex, lr, tl, Lh in hcases:
        turn_m = tl * N * nc                                          # meters of turn facing the sled
        s0 = state(N, Lh, w, I_, ex, 0.0); s1 = state(N, Lh, w, I_, ex, B_sled_end)
        P_site = s0["q"] * turn_m * f_end / n_sites_sled; P_site_s = s1["q"] * turn_m * f_end / n_sites_sled
        E_l = q_launch(I_, ex, lr["_v"], lr["_dt"], False, Lh) * turn_m; E_ls = q_launch(I_, ex, lr["_v"], lr["_dt"], True, Lh) * turn_m
        kgm = tl * N * Lh * w * lh.TAPE_T * lh.TAPE_RHO * nc / cfg.L_SLED
        row = dict(label=lab, coils=nc, I=I_, layers=Lh, q_J_per_m=s0["q"], share_of_limit=s0["share"], hours_of_spare=E_l / (cfg.L_RING / SITE) / SPARE_W / 3600, turn_m_under_sled=turn_m, tape_kg_per_m=kgm, tape_ring_t=kgm * cfg.L_RING / 1e3, P_site_W=P_site, E_pass_J=P_site * t_pass, E_launch_J=E_l, share=E_l / KE,
                   E_site_mean_J=E_l / (cfg.L_RING / SITE), P_site_sled_W=P_site_s, E_launch_sled_J=E_ls, share_sled=E_ls / KE, t_launch_min=lr["t_min"], laps=lr["laps"],
                   compressor_share=E_l * COP_PENALTY / KE, compressor_share_sled=E_ls * COP_PENALTY / KE)
        hyst.append(row)
        out(f"     {lab:<50} {nc:6.0f} {I_:7.0f} {turn_m/1e3:8.1f} {kgm:9.2f} {kgm*cfg.L_RING/1e3:11,.0f} | {P_site/1e3:11.2f} {P_site*t_pass/1e3:11.2f} {E_l:11.3e} {E_l/KE:9.2e} {E_l/(cfg.L_RING/SITE)/1e3:14.1f} | {P_site_s/1e3:19.2f} {E_ls:11.3e} {E_ls/KE:9.2e}")
    out(f"     'per site' is the power of all the coils facing the sled divided by the {n_sites_sled:.0f} sites it covers, while the sled is over the site.  'with sled field' adds the sled's field")
    out(f"     through the tape faces to the winding's own (amplitudes added, {B0*1e3:.0f} mT below {v_cross/1e3:.0f} km/s falling to {B_sled_end*1e3:.1f} mT at 30 km/s): an upper bound [UNVERIFIED].")
    out(f"     Tape mass: conductor {lh.TAPE_T*1e6:.0f} micrometers thick at {lh.TAPE_RHO:.0f} kg/m3 (lim_hysteresis.py, the basis of Table 9.2's 1,293 kg per LIM); one coil holds {m_coil:.0f} kg.")
    out(f"     Compressor work at {COP_PENALTY:.0f} W per W, as a share of the launch's kinetic energy: " + "; ".join(f"{h['label'].split(':')[0]} {h['compressor_share']*100:.3f}%" for h in hyst[1:]) + ".")
    out(f"     Row B5: {Lk} layers at {I_book:.0f} A, worst tape at {s_k['share']*100:.1f} percent of its limit, q = {s_k['q']*1e3:.1f} mJ/m, {n_need:,.0f} coils.")
    out(f"     Time for the {SPARE_W/1e3:.1f} kW of spare plant capacity per site (Section 14.5's figure) to remove the mean heat of one launch: " + "; ".join(f"{h['label'].split(':')[0]} {h['hours_of_spare']:.2f} h" for h in hyst) + ".")
    res["hysteresis"] = hyst; res["coil_tape_kg"] = m_coil; res["KE_J"] = KE; res["t_pass_s"] = t_pass; res["f_end_Hz"] = f_end; res["B_sled_end_T"] = B_sled_end

    out()
    out("  Limits of this run: the field at each tape is two-dimensional (one straight run of the coil, with the return run where stated); the other two phases' runs, which lie")
    out("  along the same edges of the stator, are not added; the coil is taken 2 m wide between the centers of its 653 mm sheets; the sled's field is taken at its gap value over")
    out("  all of the tape; eddy loss in the tape's metal and coupling between the soldered layers are in neither model; the tape data is 70 K (lsm_config has the stator at 77 K).")

    with open(os.path.join(LOGS, "round10_stator.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, "round10_stator.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  written: {os.path.join(LOGS, 'round10_stator.txt')} and .json")


if __name__ == "__main__":
    main()
