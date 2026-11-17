#!/usr/bin/env python3
"""
Round 11, item 2: a stator for the mass driver that carries a stated launch (Chapters 9, 10 and 14).
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 9.10, 9.11, 10.1 (Equations 10.1 to 10.3, Table 10.2), Section 14.5.
Brief: Vol_III/Claude_Code_instructions_round11_2026-11-17.md, Section 2.
Run: python tools/round11_stator.py --jobs=5   (about 15 minutes with five worker processes, the default is four; writes
                                               runs/round11_logs/round11_stator_checks.txt, round11_stator_study.txt and their .json)
     python tools/round11_stator.py --quick    (parts 1 to 3 only, about a minute; writes round11_stator_checks.txt and .json)
     tools/round11_stator_crosscheck.py reads the study's json afterwards (a second loss model, braking, a racetrack sled).

Parts:
  1. The stator of Table 9.6 recomputed: thrust by the power balance, by the current sheet and by a direct sum over the
     conductors; the coil's self, mutual and cyclic inductance and its voltage; the tape limit and the strip loss of round 10.
  2. The tape below 70 K: what is measured, the lift factors used, the compressor work per watt.
  3. The strip loss against the brief's expectation (a colder tape saves tape but not heat): one tape in a field.
  4. The design study: four launches of the 5,000 tonne spacecraft, each on the stator that carries it with the least heat,
     over the pole pitch, the lay, the tape width, the end band, the setback and the split between the sled's field and the
     coil current above the corner speed; then the temperature, one side or two, the 20 km sled, the light launch.
  5. The heat of each and the plant.
The model is tools/round11_stator_model.py.  Nothing but the logs is written.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the repository root
import matplotlib
matplotlib.use("Agg")
import numpy as np

import lsm_config as cfg
import lsm_physics as phys
import lim_hysteresis as lh
import ring_altitude as ring
import round11_stator_model as m
import round10_stator as r10

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round11_logs")
SPARE_W = 11.9e3                # W per site at 70 K, the plant's spare in a launch window (Section 14.5; the book's figure)
BEARING_SITE_J = 6.8e6          # J per site per launch, the bearing's heat (Section 14.5; the book's figure)
MU0 = m.MU0
_lines = []
RES = {}


def out(s=""):
    print(s, flush=True)
    _lines.append(s)


def save(name):
    with open(os.path.join(LOGS, f"round11_stator_{name}.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, f"round11_stator_{name}.json"), "w", encoding="utf-8") as fh:
        json.dump(RES, fh, indent=1, default=float)


# ------------------------------------------------------------------------------------------------ a direct check of the spread winding
def direct_sheet_check(tau, N, I, q, B0, v, delta, w=m.W, n_t=360, n_u=24, n_x=240):
    """The two-layer winding with coil sides spread evenly over their 60 degree belts, by direct sums: the force on the
    belts in the sled's field (Lorentz force, per pole pair of one side), the back-EMF of one coil from its flux linkage
    (numerical differences), and the electrical power going into the back-EMFs of the coils of the pole pair."""
    beta = tau / 3.0; om = math.pi * v / tau; T = 2 * math.pi / om
    t = np.linspace(0.0, T, n_t, endpoint=False); dt = t[1] - t[0]
    Kb = 6 * q * N * I / tau
    B = lambda x, tt: B0 * np.sin(math.pi * x / tau - om * tt)
    F = np.zeros(n_t); P = np.zeros(n_t); e_pk = 0.0
    xq = (np.arange(n_x) + 0.5) / n_x
    uq = (np.arange(n_u) + 0.5) / n_u * beta
    for b in range(6):
        ph = math.pi * (b + 0.5) / 3.0
        i_b = I * np.cos(om * t - ph - delta)                 # the coils whose go side lies in belt b
        # force on belt b: go sides of these coils and return sides of the coils of belt b - 3 (equal current density, same sign)
        xb = (b + xq) * beta
        F += np.array([np.mean(B(xb, tt)) for tt in t]) * beta * (Kb * np.cos(om * t - ph - delta)) * w
        # flux linkage of one coil of belt b: go side at b beta + u, return side one pole pitch on, averaged over u
        lam = np.array([N * w * np.mean([np.mean(B(b * beta + u + xq * tau, tt)) * tau for u in uq]) for tt in t])
        e = -(np.roll(lam, -1) - np.roll(lam, 1)) / (2 * dt)
        e_pk = max(e_pk, float(np.abs(e).max()))
        P += q * (-e * i_b)
    return dict(F_mean=float(abs(F.mean())), F_ripple=float(F.max() - F.min()), e_peak=e_pk, P_mean=float(abs(P.mean())), Fv=float(abs(F.mean()) * v), K_belt=Kb)


# ------------------------------------------------------------------------------------------------ a coil with its inductive voltage counted
def launch_voltage_counted(n_coils, N, I_max, L_c, tau, w, B0_max, V, delta, mass, v_f, n=3000):
    """The launch of n_coils full-pitch coils of N turns and cyclic inductance L_c (H), each held to V at its terminals with
    the inductive voltage counted: at each speed the sled's field (up to B0_max) and the coil current (up to I_max) that give
    the most power with |V|^2 = E^2 + X^2 - 2 E X cos(delta) <= V^2, E = 2 N v B w, X = omega L_c I.  Returns the time,
    distance, peak thrust, the speed where the voltage first binds, and the power, current and sled field at v_f."""
    ve = np.linspace(0.0, v_f, n + 1); v = 0.5 * (ve[1:] + ve[:-1]); dv = ve[1] - ve[0]
    c = math.cos(delta)
    Pv = np.zeros(n); Iv = np.zeros(n); Bv = np.zeros(n); bind = None
    for k, vk in enumerate(v):
        E_max = 2 * N * vk * B0_max * w; X1 = math.pi * vk / tau * L_c
        if E_max**2 + (X1 * I_max)**2 - 2 * E_max * X1 * I_max * c <= V * V:
            E, I = E_max, I_max
        else:
            if bind is None:
                bind = vk
            Es = np.linspace(1e-3 * V, min(E_max, V / math.sin(delta) * 0.999999), 4000)
            X = Es * c + np.sqrt(np.maximum(V * V - Es * Es * (1 - c * c), 0.0))       # the larger root of X^2 - 2 E X c + E^2 - V^2 = 0
            Is = np.minimum(I_max, X / X1)
            j = int(np.argmax(Es * Is)); E, I = float(Es[j]), float(Is[j])
        Pv[k] = n_coils * 0.5 * E * I * math.sin(delta); Iv[k] = I; Bv[k] = E / (2 * N * vk * w)
    F = Pv / v
    dt = mass * dv / F
    return dict(t=float(dt.sum()), x=float((v * dt).sum()), F_pk=float(F.max()), v_bind=bind, P_end=float(Pv[-1]), P_max=float(Pv.max()), I_end=float(Iv[-1]), B_end=float(Bv[-1]))


# ------------------------------------------------------------------------------------------------ part 1
def part1():
    N = cfg.N_STATOR; tau = cfg.TAU_P; wc = cfg.W_COIL; B0 = cfg.B_SLED_NOMINAL; I = cfg.I_PEAK; d = cfg.DELTA_MAX; V = cfg.V_COIL_LIMIT
    L_act = cfg.N_UNITS * cfg.N_POLES_PER_UNIT * tau
    v_c = V / (2 * N * B0 * wc)
    B_code = phys.calc_B_stator(I, N, wc, cfg.G_GAP); F_code = phys.calc_thrust(B_code, B0, d, cfg.A_ACTIVE_TOTAL)
    out("=" * 150)
    out("  ROUND 11, ITEM 2, PART 1: THE STATOR OF TABLE 9.6 RECOMPUTED")
    out("=" * 150)
    out(f"  lsm_config.py: {N} turns, {cfg.N_HTS_LAYERS} layers of {cfg.TAPE_WIDTH_MM:.0f} mm tape, pole pitch {tau:.0f} m, coil {wc:.0f} m wide, gap {cfg.G_GAP*1e3:.0f} mm, {V/1e3:.0f} kV, sled field {B0:.2f} T, I_PEAK {I:.1f} A,")
    out(f"  load angle {math.degrees(d):.0f} deg, active length {L_act:,.0f} m on each of {cfg.N_LIM_SIDES} sides ({cfg.A_ACTIVE_TOTAL:,.0f} m2), launch mass {cfg.M_TOTAL/1e3:,.0f} t.")
    out(f"  Equations 10.1 and 10.2 as coded: B_stator {B_code*1e3:.2f} mT, thrust {F_code/1e6:.2f} MN, power at {v_c:,.0f} m/s {F_code*v_c/1e9:.2f} GW.")
    out()
    out("  1.1 The thrust three ways (concentrated coil sides, as round 10 took them: winding factor 1)")
    out(f"      {'coils per phase per pole pitch per side':<42} {'coils':>6} {'K1 A/m':>9} {'sheet field at the gap mT':>26} {'thrust, sheet MN':>17} {'thrust, power balance MN':>25} {'direct sum MN':>14} {'power GW':>9} {'book / this':>11}")
    rows = []
    c = r10.direct_check(N, I, B0, wc, tau, v_c, d)
    f_coil_direct = c["F_mean"] / 3.0
    for lab, q in (("1/3 (one coil per phase per 3-pitch unit)", 1.0 / 3.0), ("1/2 (single layer, fully wound)", 0.5), ("1 (two layers, fully wound)", 1.0)):
        n = cfg.N_LIM_SIDES * 3 * q * L_act / tau
        K1 = phys.calc_sheet_current_full_pitch(N, I, q, tau)
        B_gap = 0.5 * MU0 * K1 * math.exp(-math.pi * cfg.G_GAP / tau)
        F_sheet = phys.calc_thrust_current_sheet(K1, B0, d, cfg.A_ACTIVE_TOTAL)
        F_pb = n * 0.5 * phys.calc_EMF(N, v_c, B0, wc) * I * math.sin(d) / v_c
        F_dir = n * f_coil_direct
        rows.append(dict(q=q, coils=n, K1=K1, B_gap=B_gap, F_sheet=F_sheet, F_power_balance=F_pb, F_direct=F_dir, P=F_pb * v_c, ratio=F_code / F_sheet))
        out(f"      {lab:<42} {n:6.0f} {K1:9.0f} {B_gap*1e3:26.2f} {F_sheet/1e6:17.3f} {F_pb/1e6:25.3f} {F_dir/1e6:14.3f} {F_pb*v_c/1e9:9.2f} {F_code/F_sheet:11.1f}")
    K1_need = 2 * F_code / (B0 * math.sin(d) * cfg.A_ACTIVE_TOTAL)
    out(f"      K1 = 6 q N I / tau is the peak of the fundamental of the winding's current sheet (A per meter of track); the field of that sheet alone is mu0 K1 / 2 at its")
    out(f"      face.  Thrust = (1/2) K1 B_sled sin(delta) x area.  Equation 10.2's {F_code/1e6:.1f} MN needs K1 = {K1_need:,.0f} A/m, a sheet field of {0.5*MU0*K1_need*1e3:.1f} mT (mu0 K1 = {MU0*K1_need*1e3:.1f} mT, the")
    out(f"      figure Equation 10.1 returns); with {N} turns at {I:.0f} A that is q = {K1_need*tau/(6*N*I):.1f} coils per phase per pole pitch per side.")
    RES["table96_thrust"] = dict(rows=rows, K1_need=K1_need, F_code=F_code, B_code=B_code, direct_check=c)

    out()
    out("  1.2 The coil's inductance and its voltage (Neumann's integral over straight filaments; conductor band 653 mm wide, geometric mean distance 0.2235 x width)")
    bw = lh.winding_width("flat", N, cfg.N_HTS_LAYERS, cfg.TAPE_WIDTH)
    g = 0.2235 * (bw + cfg.N_HTS_LAYERS * lh.TAPE_T)
    L1 = float(m.m_loops(0.0, tau, wc, g, g)); Lg = m.grover_rectangle(tau, wc, g)
    L2w = MU0 / math.pi * (tau * math.log(wc / g) + wc * math.log(tau / g))
    out(f"      self inductance of one turn: Neumann {L1*1e6:.2f} uH, Grover's rectangle formula {Lg*1e6:.2f} uH, round 10's two-wire formula {L2w*1e6:.2f} uH;")
    out(f"      with {N} turns: {L1*N*N*1e3:.1f} mH (Neumann), {Lg*N*N*1e3:.1f} mH (Grover), {L2w*N*N*1e3:.1f} mH (round 10).")
    out(f"      mutual inductance to the other coils of the side, two-layer winding with one coil per phase per pole pitch (coils every {tau/3:.1f} m), {N} turns:")
    for j in (1, 2, 3, 4, 5, 6):
        Mj = float(m.m_loops(j * tau / 3.0, tau, wc, g, g))
        out(f"        coil {j} places on ({j*tau/3:6.1f} m, phase {60*j:3d} deg): M = {Mj*N*N*1e3:8.2f} mH, weighted by cos(phase): {Mj*N*N*1e3*math.cos(math.pi*j/3):8.2f} mH")
    Lc, Ls = m.cyclic_inductance(tau, wc, 1, g, g)
    out(f"      cyclic inductance (self plus every mutual times the cosine of its phase): {Lc*N*N*1e3:.1f} mH, {Lc/Ls:.2f} times the self inductance.")
    vb = []
    out(f"      {'speed':>8} {'f Hz':>7} {'sled field mT':>13} {'back-EMF kV':>12} {'own inductive kV':>17} {'cyclic inductive kV':>20} {'terminal kV, stator flux against the sled':>42} {'with it':>9}")
    for v in (v_c, cfg.V_LAUNCH):
        f = v / (2 * tau); om = 2 * math.pi * f
        Bs = min(B0, V / (2 * N * v * wc)); E = phys.calc_EMF(N, v, Bs, wc)
        Xs = om * Ls * N * N * I; Xc = om * Lc * N * N * I
        Vd = phys.calc_coil_voltage(E, Xc, d, demagnetizing=True); Vm = phys.calc_coil_voltage(E, Xc, d, demagnetizing=False)
        vb.append(dict(v=v, f=f, B_sled=Bs, E=E, X_self_I=Xs, X_cyclic_I=Xc, V_demag=Vd, V_mag=Vm))
        out(f"      {v:8,.0f} {f:7.2f} {Bs*1e3:13.1f} {E/1e3:12.1f} {Xs/1e3:17.1f} {Xc/1e3:20.1f} {Vd/1e3:42.1f} {Vm/1e3:9.1f}")
    out(f"      The terminal voltage is |E + j omega L I| with the current {90-math.degrees(d):.0f} degrees from the back-EMF (the {math.degrees(d):.0f} degree load angle), on either side of it.")
    RES["table96_inductance"] = dict(gmd=g, L_turn_neumann=L1, L_turn_grover=Lg, L_turn_round10=L2w, L_self_H=Ls * N * N, L_cyclic_H=Lc * N * N, voltage=vb)

    out()
    out("  1.3 The tape limit and the strip loss of round 10, run again (functions of tools/round10_stator.py and lim_hysteresis.py)")
    Lyr = cfg.N_HTS_LAYERS; w = cfg.TAPE_WIDTH
    extra = r10.return_run_field(N, Lyr, w, wc)
    I50, st50 = r10.limit(N, Lyr, w, extra=extra)
    I50s, st50s = r10.limit(N, Lyr, w, extra=extra, B_sled=B0)
    L9 = next(L_ for L_ in range(5, 14) if r10.limit(N, L_, w, extra=r10.return_run_field(N, L_, w, wc))[0] >= I)
    n150 = cfg.LIM_PHASES * cfg.N_LIM_SIDES * cfg.N_UNITS
    Fb = phys.calc_thrust_full_pitch(N, I50, B0, wc, d, n150); Pb = phys.calc_power_limit_full_pitch(V, I50, d, n150)
    la = r10.launch(Fb, Pb, cfg.M_TOTAL, cfg.V_LAUNCH, dt=1.0)
    turn = phys.calc_turn_length_full_pitch(tau, wc)
    q1 = r10.state(N, Lyr, w, I50, extra, 0.0)["q"]
    E_l = q1 * float(np.sum(la["v"] / (2 * tau)) * la["dt"]) * turn * N * n150
    KE = 0.5 * cfg.M_TOTAL * cfg.V_LAUNCH**2
    out(f"      five layers of 12 mm tape, flat, own field and return run: {I50:.1f} A at 50 percent ({I50/I:.3f} of {I:.0f} A); with the sled's {B0*1e3:.0f} mT added: {I50s:.1f} A; layers that keep {I:.0f} A: {L9}.")
    out(f"      150 coils at {I50:.0f} A: thrust {Fb/1e6:.2f} MN, power {Pb/1e9:.2f} GW, launch {la['t']/3600:.1f} h; strip loss {q1*1e3:.1f} mJ per meter of turn per cycle, {E_l:.3e} J per launch, {E_l/KE*100:.2f} percent of {KE/1e15:.3f} PJ.")
    RES["table96_tape"] = dict(I_50=I50, I_50_sled=I50s, layers_keep=L9, F=Fb, P=Pb, t_h=la["t"] / 3600, q_J_per_m=q1, E_launch=E_l, share=E_l / KE)
    out()
    out("  1.4 The stator of Table 9.6 with the coil's inductive voltage counted in the 100 kV (function launch_voltage_counted): at each speed the sled's field and the")
    out("      coil current that give the most power inside 100 kV at the terminals, the stator's flux against the sled's at the 60 degree load angle")
    out(f"      {'coils':>6} {'coil A at most':>14} {'inductance mH':>14} {'peak thrust MN':>15} {'voltage binds from km/s':>24} {'power at most GW':>17} {'power at 30 km/s GW':>20} {'coil A at 30 km/s':>18} {'sled field at 30 km/s mT':>25} {'launch h':>9} {'laps':>7}")
    g2 = g
    L_150 = (float(m.m_loops(0.0, tau, wc, g2, g2)) + 2 * float(m.m_loops(2 * tau / 3.0, tau, wc, g2, g2)) * math.cos(2 * math.pi / 3)) * N * N
    vc_rows = []
    for n_c, Imax, Lc_ in ((n150, I, L_150), (n150, I50, L_150), (3 * n150, I, Lc * N * N), (3 * n150, I50, Lc * N * N)):
        r = launch_voltage_counted(n_c, N, Imax, Lc_, tau, wc, B0, V, d, cfg.M_TOTAL, cfg.V_LAUNCH)
        vc_rows.append(dict(coils=n_c, I_max=Imax, L=Lc_, **r))
        out(f"      {n_c:6d} {Imax:14.0f} {Lc_*1e3:14.1f} {r['F_pk']/1e6:15.2f} {r['v_bind']/1e3:24.2f} {r['P_max']/1e9:17.2f} {r['P_end']/1e9:20.2f} {r['I_end']:18.0f} {r['B_end']*1e3:25.1f} {r['t']/3600:9.1f} {r['x']/cfg.L_RING:7.1f}")
    out(f"      The 150 coils are taken as round 10 placed them (three coils 2/3 of a pole pitch apart in each unit; cyclic inductance {L_150*1e3:.1f} mH); the 450 as the two-layer winding of 1.2.")
    RES["table96_voltage_counted"] = vc_rows
    return F_code, v_c


# ------------------------------------------------------------------------------------------------ part 2
def part2():
    out()
    out("=" * 150)
    out("  PART 2: THE TAPE BELOW 70 K")
    out("=" * 150)
    out("  Measured curves in data/: THEVA Pro-Line Advanced Pinning (the book's tape) at 70 K (field 0 to 8 T, every angle; Robinson) and at 77.5 K and 70 K (four angles; THEVA).")
    out("  There is no measured curve of that tape at 60, 50 or 30 K in the repository.")
    out("  Added this round: data/THEVA_ProLine_TPL4420c_Robinson_2019_30K_50K_60K_70K.csv, the Robinson Research Institute's published measurement of a THEVA Pro-Line tape")
    out("  of 2019 (4 mm wide, TPL4420c; Wimbush and Strickland, figshare, doi 10.6084/m9.figshare.3759327.v2, CC BY 4.0), 0 to 7 T, every angle, at 30, 50, 60 and 70 K.")
    out("  It is not the Advanced Pinning tape.  Used here only as a ratio: Ic(T, B, angle) / Ic(70 K, B, angle) of the 2019 sample, applied to the AP tape's 70 K curve [UNVERIFIED for the AP tape].")
    out()
    out("  2.1 The two samples at 70 K, A per cm of width (a check that the 2019 sample is a comparable tape)")
    out(f"      {'field T':>8} {'AP tape, through the face':>26} {'2019 sample':>12} {'AP tape, along the face':>24} {'2019 sample':>12}")
    cmp_rows = []
    for B in (0.0, 0.03, 0.1, 0.3, 1.0, 3.0):
        a0 = float(lh.ic_tape(np.array([B]), np.array([0.0]), 0.01)[0]); a90 = float(lh.ic_tape(np.array([B]), np.array([90.0]), 0.01)[0])
        s0 = float(m.ic19(70.0, np.array([B]), np.array([0.0]))[0]); s90 = float(m.ic19(70.0, np.array([B]), np.array([90.0]))[0])
        cmp_rows.append(dict(B=B, ap_0=a0, s19_0=s0, ap_90=a90, s19_90=s90))
        out(f"      {B:8.2f} {a0:26.0f} {s0:12.0f} {a90:24.0f} {s90:12.0f}")
    out()
    out("  2.2 The lift factor Ic(T) / Ic(70 K) of the 2019 sample, field through the face (0 deg) and along it (90 deg)")
    out(f"      {'field T':>8} | {'60 K':>6} {'50 K':>6} {'30 K':>6} | {'60 K':>6} {'50 K':>6} {'30 K':>6}")
    lf = []
    for B in (0.0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0):
        r0 = [float(m.lift(T, np.array([B]), np.array([0.0]))[0]) for T in (60.0, 50.0, 30.0)]
        r9 = [float(m.lift(T, np.array([B]), np.array([90.0]))[0]) for T in (60.0, 50.0, 30.0)]
        lf.append(dict(B=B, through=r0, along=r9))
        out(f"      {B:8.2f} | {r0[0]:6.2f} {r0[1]:6.2f} {r0[2]:6.2f} | {r9[0]:6.2f} {r9[1]:6.2f} {r9[2]:6.2f}")
    out("      Chapter 14 (Section 14.3) gives 1.7 at 60 K and 2.5 at 50 K against 70 K.  The 2019 sample gives those figures near 0.2 to 0.3 T; with no field it gives 1.52 and 2.04.")
    out()
    out("  2.3 Compressor work per watt lifted to the 400 K radiator at 5 percent of Carnot (Equation 14.1's rule), and the spare of Section 14.5 at each temperature")
    out(f"      {'T K':>5} {'Carnot COP':>11} {'W per W':>8} {'against 70 K':>12} {'the 11.9 kW spare as compressor power MW':>42} {'heat that power lifts kW':>25}")
    cw = []
    for T in (70.0, 60.0, 50.0, 30.0):
        wpw = m.work_per_watt(T); sp = SPARE_W * m.work_per_watt(70.0)
        cw.append(dict(T=T, carnot=T / (m.T_HOT - T), wpw=wpw, spare_heat_W=sp / wpw))
        out(f"      {T:5.0f} {T/(m.T_HOT-T):11.4f} {wpw:8.1f} {wpw/m.work_per_watt(70.0):12.2f} {sp/1e6:42.3f} {sp/wpw/1e3:25.2f}")
    RES["tape_T"] = dict(compare70=cmp_rows, lift=lf, compressor=cw)


# ------------------------------------------------------------------------------------------------ part 3
def part3():
    out()
    out("=" * 150)
    out("  PART 3: THE STRIP LOSS AGAINST THE BRIEF'S EXPECTATION (ONE TAPE IN A FIELD THROUGH ITS FACE)")
    out("=" * 150)
    out("  Brandt and Indenbom: Q = 4 mu0 (Ic/w) (w/2)^2 H f(h), h = H / H_c, mu0 H_c = mu0 Ic / (pi w) the penetration field; Q -> w Ic B for h >> 1, Q -> w Ic B h^3 / 6 for h << 1.")
    out("  'per kA rated' is the loss of the tapes that carry 1 kA at the 50 percent rule: Q / (0.5 Ic) x 1,000.  If a colder tape saved tape but not heat this column would not change with T.")
    out(f"      {'field mT':>9} {'tape':>6} {'T K':>4} {'Ic A':>7} {'penetration field mT':>21} {'h':>6} {'Q mJ/m per cycle':>17} {'Q / (w Ic B)':>13} {'tapes per kA':>13} {'mJ/m per cycle per kA rated':>28} {'against 70 K':>12}")
    rows = []
    for B in (0.010, 0.030, 0.050, 0.100, 0.300, 1.000):
        for wt in (0.012, 0.003):
            base = None
            for T in (70.0, 60.0, 50.0, 30.0):
                Ic = float(m.ic(np.array([B]), np.array([0.0]), wt, T)[0])
                Bp = MU0 * Ic / (math.pi * wt); h = B / Bp
                Q = float(lh.q_strip_perp(np.array([Ic]), wt, np.array([B]))[0])
                per_kA = Q / (0.5 * Ic) * 1000.0
                base = per_kA if base is None else base
                rows.append(dict(B=B, w_t=wt, T=T, Ic=Ic, B_pen=Bp, h=h, Q=Q, ratio_sat=Q / (wt * Ic * B), tapes_per_kA=1000.0 / (0.5 * Ic), Q_per_kA=per_kA, vs70=per_kA / base))
                out(f"      {B*1e3:9.0f} {wt*1e3:4.0f}mm {T:4.0f} {Ic:7.1f} {Bp*1e3:21.1f} {h:6.2f} {Q*1e3:17.3f} {Q/(wt*Ic*B):13.3f} {1000.0/(0.5*Ic):13.2f} {per_kA*1e3:28.2f} {per_kA/base:12.3f}")
    out("  Reading: at 0.3 and 1 T (h well above 1) the heat per kA rated changes little with the temperature, and the tape count falls: the brief's expectation holds there.")
    out("  At 10 to 50 mT (h near or below 1 at 70 K) the colder tape is penetrated less, and the heat per kA rated falls with the tape count: to 0.29 to 0.47 of the 70 K figure")
    out("  at 60 K, 0.12 to 0.24 at 50 K and 0.04 to 0.08 at 30 K.  At 100 mT it falls to 0.72, 0.50 and 0.23.  A 3 mm tape loses a quarter of what a 12 mm tape loses per kA rated")
    out("  at the same field and temperature (the loss per tape goes as the width squared, the current as the width).")
    RES["strip_check"] = rows


# ------------------------------------------------------------------------------------------------ part 4: the search
TAUS = (2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0, 130.0)
LAYS = (("edge", "edge", 0.0), ("edge", "flat", 0.0), ("edge", "flat", 1.0), ("flat", "flat", 0.0), ("flat", "flat", 1.0))
BES = (0.5, 1.0)
ALPHAS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
TAUS_COLD = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0)
LAYS_COLD = (("edge", "edge", 0.0), ("edge", "flat", 1.0), ("flat", "flat", 1.0))
if os.environ.get("ROUND11_STATOR_SMOKE"):          # a smoke test of the whole script on a tiny grid
    TAUS = (4.0, 6.0, 130.0); TAUS_COLD = (4.0, 6.0); BES = (1.0,); ALPHAS = (0.0, 0.4)
KEEP = ("tau", "lay", "end_lay", "w_t", "T", "b_e", "d", "sides", "alpha", "c_end", "F_pk", "P", "L_sled", "m", "v_f", "v_c", "K1", "K_belt", "B_sheet", "i_t", "n_t", "ic_act", "ic_end", "ic_leg",
        "B_act", "B_end", "tape_m_per_m", "tape_kg_per_m", "tape_ring_t", "N", "ampere_turns_pp", "f_end", "E_turn_c", "XI_turn_c", "V_turn_c", "Vmag_turn_c", "E_turn_f", "XI_turn_f",
        "V_turn_f", "Vmag_turn_f", "V_turn_max", "E_turn_max", "v_at_Vmax", "t", "x", "laps", "a_g", "KE", "E_hyst", "E_act", "E_end", "E_leg", "E_eddy", "share", "E_site", "E_low", "x_c", "E_site_start",
        "P_site_c", "P_site_f", "E_pass_f", "work_per_watt", "E_comp", "comp_share", "Lc", "Ls", "B_avg_pk", "q_pk", "q_end_speed")


def trim(r):
    return {k: r[k] for k in KEEP}


def search(F, P, w_t, T, L_sled=10000.0, m_sc=cfg.M_SPACECRAFT, v_f=cfg.V_LAUNCH, sides=2, taus=TAUS, lays=LAYS, bes=BES, alphas=ALPHAS, c_end=0.5, keep_all=False):
    """every stator of the grid on one launch; returns the best (least hysteresis heat) at each pole pitch, the best of all,
    and optionally every row.  A stator that cannot hold 100 kV with one turn is marked and never chosen."""
    by_tau = {}; allrows = []
    for tau in taus:
        best = None; n_bad = 0
        for lay, end_lay, d in lays:
            for be in bes:
                st = m.Stator(tau, lay, w_t, T, be, d, sides, c_end, end_lay)
                for al in alphas:
                    r = m.evaluate(st, F, P, L_sled, m_sc, v_f, al)
                    if r["N"] < 1:
                        n_bad += 1; continue
                    if keep_all:
                        allrows.append(trim(r))
                    if best is None or r["E_hyst"] < best["E_hyst"]:
                        best = r
        by_tau[tau] = (trim(best) if best else None, n_bad)
    ok = [b for b, _ in by_tau.values() if b is not None]
    top = min(ok, key=lambda r: r["E_hyst"]) if ok else None
    return by_tau, top, allrows


def _job(a):
    """one launch, tape width and temperature: the whole grid at 70 K, the shorter grid below it (for the worker pool)"""
    name, F, P, w_t, T = a
    if T == 70.0:
        by_tau, top, allrows = search(F, P, w_t, T, keep_all=(name == "book"))
    else:
        by_tau, top, allrows = search(F, P, w_t, T, taus=TAUS_COLD, lays=LAYS_COLD)
    return name, w_t, T, by_tau, top, allrows


def _best_job(kw):
    """the least-heat stator of the shorter pole-pitch grid for one case (3 mm tape, 70 K), for the worker pool"""
    by_tau, top, _ = search(kw["F"], kw["P"], 0.003, 70.0, L_sled=kw.get("L_sled", 10000.0), sides=kw.get("sides", 2), c_end=kw.get("c_end", 0.5), taus=TAUS_COLD)
    return top


def lay_name(r):
    a = "on edge" if r["lay"] == "edge" else "flat"; e = "on edge" if r["end_lay"] == "edge" else "flat"
    return f"{a} / {e}" + (f", set back {r['d']:.0f} m" if r["d"] > 0 else "")


def fmt_design(r):
    return (f"{r['tau']:5.0f} {lay_name(r):<26} {r['b_e']:4.1f} {r['alpha']:5.1f} | {r['i_t']:7.1f} {r['n_t']:7.0f} {r['tape_kg_per_m']:8.2f} {r['tape_ring_t']/1e3:8.0f} | {r['N']:4d} {r['f_end']:7.0f} |"
            f" {r['E_hyst']:10.3e} {r['share']*100:8.4f} {r['E_site']/1e6:9.2f} {r['comp_share']*100:8.3f} {r['E_eddy']/r['E_hyst']*100:7.1f}")


HDR = (f"{'tau m':>5} {'lay, active / end runs':<26} {'band':>4} {'alpha':>5} | {'A/tape':>7} {'tapes/m':>7} {'kg/m':>8} {'ring kt':>8} | {'N':>4} {'f Hz':>7} |"
       f" {'heat J':>10} {'% of KE':>8} {'MJ/site':>9} {'comp %KE':>8} {'eddy %':>7}")


def spare_launches(r):
    """launches a day on the 11.9 kW spare (taken as compressor power when the stator is colder than 70 K)"""
    cap = SPARE_W * m.work_per_watt(70.0) / r["work_per_watt"]
    return cap * 86400.0 / r["E_site"]


def part4(F_book, v_c):
    P_book = F_book * v_c
    targets = [("book", 1.0), ("half", 0.5), ("quarter", 0.25), ("tenth", 0.1)]
    out()
    out("=" * 150)
    out("  PART 4: THE DESIGN STUDY")
    out("=" * 150)
    out(f"  Launch: {cfg.M_SPACECRAFT/1e3:,.0f} t spacecraft on the {cfg.L_SLED/1e3:.0f} km sled ({cfg.M_TOTAL/1e3:,.0f} t) to {cfg.V_LAUNCH/1e3:.0f} km/s.  Targets: the book's {F_book/1e6:.2f} MN and {P_book/1e9:.2f} GW, and a half, a quarter and a tenth of both")
    out(f"  (the corner speed stays {v_c/1e3:.0f} km/s).  Winding: the full-pitch loop, two layers, coil sides spread over their belts (winding factor {m.KD:.4f}).")
    out(f"  Grid: pole pitch {', '.join(f'{t:g}' for t in TAUS)} m; lay of the active sides and of the end runs (on edge / on edge; on edge / flat; flat / flat; the last two with the end runs")
    out(f"  at the face or set back 1 m); end band {', '.join(f'{b:g}' for b in BES)} m wide; alpha {', '.join(f'{a:g}' for a in ALPHAS)}; tape 12 and 3 mm; 70, 60, 50, 30 K.")
    out(f"  Below 70 K the grid is shorter: pole pitch {', '.join(f'{t:g}' for t in TAUS_COLD)} m, and the lays on edge / on edge, on edge / flat set back, flat / flat set back.")
    out("  N is the most turns that keep both the terminal voltage under load and the back-EMF alone inside 100 kV through the whole launch.")
    out("  alpha: above the corner speed the sled's field is B0 s^alpha and the coil current I s^(1 - alpha), s = corner speed / speed.  alpha = 1 is the book's rule (field lowered, current held).")
    out("  Columns: A/tape at the 50 percent rule; tapes per meter of track in an active side; kg of tape per meter of ring (both sides); N the most turns that hold 100 kV;")
    out("  f the supply frequency at 30 km/s; heat the strip-model hysteresis per launch; comp %KE the compressor work for it as a share of the launch's kinetic energy; eddy % the")
    out("  eddy loss in the tape's metal (a derived estimate, not in the heat) as a percent of the hysteresis.")
    res = {name: {} for name, _ in targets}
    t0 = time.time()
    jobs = [(name, F_book * k, P_book * k, w_t, T) for name, k in targets for w_t in (0.012, 0.003) for T in (70.0, 60.0, 50.0, 30.0)]
    n_jobs = next((int(a.split("=")[1]) for a in sys.argv if a.startswith("--jobs=")), 4)
    import multiprocessing as mp
    with mp.Pool(n_jobs) as pool:
        for name, w_t, T, by_tau, top, allrows in pool.imap_unordered(_job, jobs):
            res[name][(w_t, T)] = dict(by_tau=by_tau, top=top, all=allrows)
            print(f"    [{time.time()-t0:6.0f} s] {name} {w_t*1e3:.0f} mm {T:.0f} K done", flush=True)
    # ---- 4.1 the pole pitch at 70 K
    for name, k in targets:
        out()
        out(f"  4.1 ({name}) {F_book*k/1e6:.2f} MN, {P_book*k/1e9:.2f} GW: the best stator at each pole pitch, 70 K")
        for w_t in (0.012, 0.003):
            out(f"      {w_t*1e3:.0f} mm tape")
            out("      " + HDR)
            for tau in TAUS:
                b, n_bad = res[name][(w_t, 70.0)]["by_tau"][tau]
                if b is None:
                    out(f"      {tau:5.0f} no stator of the grid holds 100 kV with one turn")
                else:
                    out("      " + fmt_design(b) + (f"   ({n_bad} of the grid's stators need less than one turn)" if n_bad else ""))
    # ---- 4.2 the lay and the tape at the book's launch
    out()
    out("  4.2 The lay of the winding and the tape's width at the book's launch, 70 K: the best over the end band and alpha at three pole pitches")
    out("      " + f"{'tape':>5} " + HDR)
    lay_tab = []
    for w_t in (0.012, 0.003):
        rows = res["book"][(w_t, 70.0)]["all"]
        for tau in (4.0, 6.0, 10.0):
            for lay, end_lay, d in LAYS:
                sel = [r for r in rows if r["tau"] == tau and r["lay"] == lay and r["end_lay"] == end_lay and r["d"] == d]
                if not sel:
                    continue
                b = min(sel, key=lambda r: r["E_hyst"]); lay_tab.append(b)
                out("      " + f"{w_t*1e3:3.0f}mm " + fmt_design(b))
    # ---- 4.3 the temperature
    out()
    out("  4.3 The temperature: the least-heat stator at each temperature and tape width, each launch")
    out("      " + f"{'launch':<8} {'tape':>5} {'T K':>4} " + HDR + f" {'launches/day on the spare':>25}")
    for name, k in targets:
        for w_t in (0.012, 0.003):
            for T in (70.0, 60.0, 50.0, 30.0):
                b = res[name][(w_t, T)]["top"]
                out("      " + f"{name:<8} {w_t*1e3:3.0f}mm {T:4.0f} " + fmt_design(b) + f" {spare_launches(b):25.1f}")
    # ---- 4.3b the temperature at a fixed geometry (what cooling alone does)
    out()
    out("  4.3b The same stator geometry cooled (the 70 K, 3 mm and 12 mm least-heat stators of the book's launch, rerated at each temperature)")
    out("      " + f"{'tape':>5} {'T K':>4} " + HDR + f" {'tape / 70 K':>11} {'heat / 70 K':>11} {'comp work / 70 K':>16}")
    fixed = []
    for w_t in (0.012, 0.003):
        b0 = res["book"][(w_t, 70.0)]["top"]
        for T in (70.0, 60.0, 50.0, 30.0):
            st = m.Stator(b0["tau"], b0["lay"], w_t, T, b0["b_e"], b0["d"], 2, 0.5, b0["end_lay"])
            r = trim(m.evaluate(st, F_book, P_book, 10000.0, cfg.M_SPACECRAFT, cfg.V_LAUNCH, b0["alpha"])); fixed.append(r)
            out("      " + f"{w_t*1e3:3.0f}mm {T:4.0f} " + fmt_design(r) + f" {r['tape_kg_per_m']/b0['tape_kg_per_m']:11.3f} {r['E_hyst']/b0['E_hyst']:11.3f} {r['E_comp']/b0['E_comp']:16.3f}")
    RES["search"] = {name: {f"{int(w*1e3)}mm_{int(T)}K": dict(by_tau={str(t): v[0] for t, v in d["by_tau"].items()}, top=d["top"]) for (w, T), d in res[name].items()} for name in res}
    RES["lay_table"] = lay_tab; RES["fixed_geometry_T"] = fixed
    return res, targets, P_book


def chosen(res, name):
    """the stator taken for a launch: 3 mm tape at 70 K (the plant as it is), least heat"""
    return res[name][(0.003, 70.0)]["top"]


def part5(res, targets, F_book, P_book):
    out()
    out("=" * 150)
    out("  PART 5: THE FOUR STATORS (3 mm TAPE, 70 K, LEAST HEAT), THEIR LAUNCHES AND THEIR HEAT")
    out("=" * 150)
    v1 = cfg.V_ORBIT - ring.V_GROUND_SYNC; v2 = 9800.0 - ring.V_GROUND_SYNC
    spec = {}
    out("  5.1 The stator and the coil (the coil's layers are chosen to bring its current nearest the 2,107 A of Table 9.6; q follows)")
    out(f"      {'launch':<8} {'tau m':>5} {'lay':<26} {'N':>4} {'layers':>6} {'coil A':>8} {'q per phase per pitch per side':>31} {'coils facing the sled':>22} {'K1 A/m':>8} {'sheet field mT':>14} {'f at 30 km/s Hz':>16} {'tape kg/m':>10} {'ring kt':>8}")
    for name, k in targets:
        r = chosen(res, name); c = m.coil_split(r); spec[name] = dict(r=r, c=c)
        out(f"      {name:<8} {r['tau']:5.0f} {lay_name(r):<26} {r['N']:4d} {c['layers']:6d} {c['I']:8.1f} {c['q']:31.2f} {c['coils']:22,.0f} {r['K1']:8.0f} {r['B_sheet']*1e3:14.2f} {r['f_end']:16.0f} {r['tape_kg_per_m']:10.2f} {r['tape_ring_t']/1e3:8.0f}")
    out()
    out("  5.2 The voltage budget of a coil of N turns (kV): back-EMF, the coil's inductive voltage (self and mutual, the cyclic inductance), and the terminal voltage")
    out(f"      {'launch':<8} {'at the corner speed: EMF':>25} {'inductive':>10} {'terminal':>9} {'at 30 km/s: EMF':>16} {'inductive':>10} {'terminal':>9} {'largest terminal':>17} {'at km/s':>8} {'largest EMF':>12} {'terminal if the stator flux adds, at 30 km/s':>44} {'sled field at 30 km/s mT':>25} {'coil A at 30 km/s':>18}")
    for name, k in targets:
        r, c = spec[name]["r"], spec[name]["c"]; N = r["N"]
        s = r["v_c"] / r["v_f"]
        out(f"      {name:<8} {c['E_c']/1e3:25.1f} {c['XI_c']/1e3:10.1f} {c['V_c']/1e3:9.1f} {c['E_f']/1e3:16.1f} {c['XI_f']/1e3:10.1f} {c['V_f']/1e3:9.1f} {c['V_max']/1e3:17.1f} {r['v_at_Vmax']/1e3:8.1f} {c['E_max']/1e3:12.1f} {N*r['Vmag_turn_f']/1e3:44.1f} {m.B0N*s**r['alpha']*1e3:25.1f} {c['I']*s**(1-r['alpha']):18.1f}")
    out()
    out("  5.3 The thrust two ways, and a direct sum over the conductors")
    out(f"      {'launch':<8} {'sheet: (1/2) K1 B0 sin60 x area MN':>34} {'power balance: coils x (1/2) E I sin60 at the corner GW':>56} {'target GW':>10} {'direct sum: force MN':>21} {'EMF per coil kV':>16} {'formula kV':>11} {'force x speed GW':>17} {'sum of e i GW':>14}")
    for name, k in targets:
        r, c = spec[name]["r"], spec[name]["c"]
        L_act = r["L_sled"] * m.ACTIVE
        F_sheet = r["sides"] * phys.calc_thrust_current_sheet(r["K1"], m.B0N, m.DELTA, m.W * L_act)
        dc = direct_sheet_check(r["tau"], r["N"], c["I"], c["q"], m.B0N, r["v_c"], m.DELTA)
        n_pp = r["sides"] * L_act / (2 * r["tau"])                       # pole pairs facing the sled
        spec[name]["direct"] = dc
        out(f"      {name:<8} {F_sheet/1e6:34.3f} {c['P_balance_c']/1e9:56.3f} {r['P']/1e9:10.3f} {dc['F_mean']*n_pp/1e6:21.3f} {dc['e_peak']/1e3:16.2f} {c['E_c']/1e3:11.2f} {dc['Fv']*n_pp/1e9:17.3f} {dc['P_mean']*n_pp/1e9:14.3f}")
    out()
    out("  5.4 The launches")
    out(f"      {'launch':<8} {'thrust MN':>10} {'power GW':>9} {'accel g':>8} {'time min':>9} {'h':>6} {'distance km':>12} {'laps':>6} | {'empty sled g':>12} {'to 7,504 m/s min':>17} {'km':>7} {'to 9,800 m/s min':>17} {'km':>7} | {'light: thrust MN':>17} {'GW':>6} {'g':>6} {'min':>6} {'km':>8} {'heat J':>10} {'% of KE':>8}")
    for name, k in targets:
        r = spec[name]["r"]
        e = m.kinematics(r["F_pk"], r["P"], cfg.M_SLED_HARDWARE, v2, marks=(v1, v2))
        st = m.Stator(r["tau"], r["lay"], r["w_t"], r["T"], r["b_e"], r["d"], r["sides"], r["c_end"], r["end_lay"])
        lt = m.evaluate(st, r["F_pk"] / 10.0, r["P"] / 10.0, 1000.0, 500e3, 12000.0, r["alpha"])
        spec[name]["empty"] = dict(a_g=r["F_pk"] / cfg.M_SLED_HARDWARE / cfg.G_0, t_orbit=e["marks"][v1][0], x_orbit=e["marks"][v1][1], t_98=e["marks"][v2][0], x_98=e["marks"][v2][1])
        spec[name]["light"] = trim(lt)
        out(f"      {name:<8} {r['F_pk']/1e6:10.2f} {r['P']/1e9:9.2f} {r['a_g']:8.3f} {r['t']/60:9.1f} {r['t']/3600:6.2f} {r['x']/1e3:12,.0f} {r['laps']:6.2f} | {r['F_pk']/cfg.M_SLED_HARDWARE/cfg.G_0:12.3f} {e['marks'][v1][0]/60:17.1f} {e['marks'][v1][1]/1e3:7,.0f} {e['marks'][v2][0]/60:17.1f} {e['marks'][v2][1]/1e3:7,.0f} |"
            f" {lt['F_pk']/1e6:17.3f} {lt['P']/1e9:6.2f} {lt['a_g']:6.3f} {lt['t']/60:6.1f} {lt['x']/1e3:8,.0f} {lt['E_hyst']:10.3e} {lt['share']*100:8.4f}")
    out(f"      The empty sled is the {cfg.M_SLED_HARDWARE/1e3:,.0f} t sled alone; the two speeds are in space ({v1:,.0f} and {v2:,.0f} m/s over the casing).  The light launch is 500 t to 12 km/s over the casing on a")
    out("      1 km sled of its own (525 t), on the same stator: a tenth of the thrust and of the power, the same current in each coil.")
    out()
    out("  5.5 The heat (strip model: Brandt and Indenbom for the field through each tape's face, Norris for the transport current; 70 K, 94.3 W per W)")
    out(f"      {'launch':<8} {'heat J':>10} {'% of KE':>8} {'active %':>8} {'end runs %':>10} {'legs %':>7} {'below corner %':>14} | {'mean per site MJ':>16} {'site kW, sled overhead at the corner':>37} {'at 30 km/s':>11} {'per pass at 30 km/s kJ':>23} {'start stretch MJ/site':>22} |"
        f" {'compressor J':>12} {'% of KE':>8} | {'launches/day on 11.9 kW':>24} {'with the bearing':>17} | {'14/day: kW/site':>15} {'comp MW/site':>12} {'ring GW':>8} | {'light launch MJ/site':>21} {'per day on spare':>17}")
    for name, k in targets:
        r = spec[name]["r"]; lt = spec[name]["light"]
        n_day = SPARE_W * 86400.0 / r["E_site"]; n_day_b = SPARE_W * 86400.0 / (r["E_site"] + BEARING_SITE_J)
        kw14 = 14 * r["E_site"] / 86400.0 / 1e3
        spec[name]["plant"] = dict(launches_per_day_spare=n_day, with_bearing=n_day_b, kW_site_14=kw14, comp_MW_site_14=kw14 * r["work_per_watt"] / 1e3, ring_GW_14=kw14 * r["work_per_watt"] * m.N_SITES / 1e6,
                                   light_per_day=SPARE_W * 86400.0 / lt["E_site"])
        out(f"      {name:<8} {r['E_hyst']:10.3e} {r['share']*100:8.4f} {r['E_act']/r['E_hyst']*100:8.1f} {r['E_end']/r['E_hyst']*100:10.1f} {r['E_leg']/r['E_hyst']*100:7.1f} {r['E_low']/r['E_hyst']*100:14.1f} | {r['E_site']/1e6:16.2f} {r['P_site_c']/1e3:37,.0f} {r['P_site_f']/1e3:11,.0f} {r['E_pass_f']/1e3:23,.0f} {r['E_site_start']/1e6:22.1f} |"
            f" {r['E_comp']:12.3e} {r['comp_share']*100:8.3f} | {n_day:24.1f} {n_day_b:17.1f} | {kw14:15.2f} {kw14*r['work_per_watt']/1e3:12.3f} {kw14*r['work_per_watt']*m.N_SITES/1e6:8.1f} | {lt['E_site']/1e6:21.3f} {SPARE_W*86400.0/lt['E_site']:17.0f}")
    out(f"      'mean per site' is the launch's heat over the ring's {m.N_SITES:,.0f} sites.  'start stretch' is the heat left below the corner speed divided by the sites of that stretch")
    out("      (the same sites every time if every launch starts at the same place).  'with the bearing' adds the bearing's 6.8 MJ per site per launch (Section 14.5).")
    out("      '14/day' is the mean heat per site that fourteen such launches a day leave, its compressor power per site, and that over the ring.")
    # ---- 5.6 the same four at the other tape and temperatures
    out()
    out("  5.6 The same table for the least-heat stator at each tape width and temperature (launches a day: the spare taken as 1.12 MW of compressor power per site)")
    out(f"      {'launch':<8} {'tape':>5} {'T K':>4} {'tau m':>5} {'kg/m':>8} {'ring kt':>8} {'heat J':>10} {'% of KE':>8} {'MJ/site':>9} {'W per W':>8} {'compressor J':>12} {'% of KE':>8} {'launches/day on the spare':>26} {'14/day: comp MW/site':>21} {'ring GW':>8}")
    alt = []
    for name, k in targets:
        for w_t in (0.012, 0.003):
            for T in (70.0, 60.0, 50.0, 30.0):
                b = res[name][(w_t, T)]["top"]
                mw = 14 * b["E_site"] / 86400.0 * b["work_per_watt"] / 1e6
                alt.append(dict(launch=name, w_t=w_t, T=T, tau=b["tau"], kg=b["tape_kg_per_m"], E=b["E_hyst"], share=b["share"], E_site=b["E_site"], E_comp=b["E_comp"], per_day=spare_launches(b), comp_MW_14=mw))
                out(f"      {name:<8} {w_t*1e3:3.0f}mm {T:4.0f} {b['tau']:5.0f} {b['tape_kg_per_m']:8.2f} {b['tape_ring_t']/1e3:8.0f} {b['E_hyst']:10.3e} {b['share']*100:8.4f} {b['E_site']/1e6:9.2f} {b['work_per_watt']:8.1f} {b['E_comp']:12.3e} {b['comp_share']*100:8.3f} {spare_launches(b):26.1f} {mw:21.3f} {mw*m.N_SITES/1e3:8.1f}")
    RES["chosen"] = {name: dict(design=spec[name]["r"], coil=spec[name]["c"], direct=spec[name]["direct"], empty=spec[name]["empty"], light=spec[name]["light"], plant=spec[name]["plant"]) for name in spec}
    RES["alternatives"] = alt
    return spec


def part6(res, spec, F_book, P_book):
    out()
    out("=" * 150)
    out("  PART 6: ONE SIDE OR TWO, THE 20 km SLED, THE THRUST AT A GIVEN POWER, AND WHAT THE RESULT RESTS ON")
    out("=" * 150)
    b = spec["book"]["r"]
    def one(label, **kw):
        a = dict(tau=b["tau"], lay=b["lay"], w_t=b["w_t"], T=b["T"], b_e=b["b_e"], d=b["d"], sides=2, c_end=0.5, end_lay=b["end_lay"]); a.update({k: v for k, v in kw.items() if k in a})
        st = m.Stator(a["tau"], a["lay"], a["w_t"], a["T"], a["b_e"], a["d"], a["sides"], a["c_end"], a["end_lay"])
        r = trim(m.evaluate(st, kw.get("F", F_book), kw.get("P", P_book), kw.get("L_sled", 10000.0), cfg.M_SPACECRAFT, cfg.V_LAUNCH, kw.get("alpha", b["alpha"])))
        out(f"      {label:<62} " + fmt_design(r) + f" {r['t']/60:8.1f} {r['laps']:6.2f} {r['a_g']:6.3f}")
        return dict(label=label, **r)
    hdr = f"      {'case':<62} " + HDR + f" {'time min':>8} {'laps':>6} {'g':>6}"
    rows = []
    out("  6.1 The book's launch on the chosen geometry, one thing changed at a time (not optimized again)")
    out(hdr)
    rows.append(one("as chosen (two sides, 10 km sled)"))
    rows.append(one("one side only (twice the sheet current on it)", sides=1))
    rows.append(one("20 km sled, same thrust and power (half the sheet current)", L_sled=20000.0))
    rows.append(one("20 km sled, same acceleration (thrust and power x 1.512)", L_sled=20000.0, F=F_book * 15.5 / 10.25, P=P_book * 15.5 / 10.25))
    rows.append(one("sled's field at the end runs: a quarter of B0 at the face", c_end=0.25))
    rows.append(one("sled's field at the end runs: all of B0 at the face", c_end=1.0))
    rows.append(one("end runs at the face (no setback)", d=0.0))
    rows.append(one("the book's control rule (alpha = 1: field lowered, current held)", alpha=1.0))
    rows.append(one("sled's field held (alpha = 0: coil current lowered)", alpha=0.0))
    out()
    out("  6.2 The best stator for each case, searched again (3 mm tape, 70 K)")
    out(hdr)
    cases2 = [("one side only", dict(sides=1)),
              ("20 km sled, same thrust and power", dict(L_sled=20000.0)),
              ("20 km sled, same acceleration", dict(L_sled=20000.0, F=F_book * 15.5 / 10.25, P=P_book * 15.5 / 10.25)),
              ("sled's field at the end runs: all of B0 at the face", dict(c_end=1.0))]
    cases3 = []
    for kP in (1.0, 0.5, 0.25, 0.1):
        for kF in (1.0, 0.5, 0.25, 0.1):
            if kF <= kP:
                cases3.append((f"{F_book*kF/1e6:6.2f} MN, {P_book*kP/1e9:7.2f} GW (corner at {P_book*kP/(F_book*kF)/1e3:5.1f} km/s)", dict(F=F_book * kF, P=P_book * kP)))
    n_jobs = next((int(a.split("=")[1]) for a in sys.argv if a.startswith("--jobs=")), 4)
    import multiprocessing as mp
    with mp.Pool(n_jobs) as pool:
        tops = pool.map(_best_job, [dict(F=F_book, P=P_book, **kw) if "F" not in kw else kw for _, kw in cases2 + cases3])
    rows2 = []; rows3 = []
    for (label, kw), top in zip(cases2, tops[:len(cases2)]):
        out(f"      {label:<62} " + fmt_design(top) + f" {top['t']/60:8.1f} {top['laps']:6.2f} {top['a_g']:6.3f}")
        rows2.append(dict(label=label, **top))
    out()
    out("  6.3 The thrust at a given power: the launch time is m v^2 / (2 P) + m P / (2 F^2), so the thrust adds little time and costs tape (3 mm, 70 K, searched again)")
    out(hdr)
    for (label, kw), top in zip(cases3, tops[len(cases2):]):
        out(f"      {label:<62} " + fmt_design(top) + f" {top['t']/60:8.1f} {top['laps']:6.2f} {top['a_g']:6.3f}")
        rows3.append(dict(label=label, **top))
    # ---- 6.4 the Table 9.6 pole pitch in the same model
    out()
    out("  6.4 The 130 m pole pitch of Table 9.6 in the same model (12 mm tape, 70 K, the best lay of the grid), each launch: what the rule and the voltage give")
    out(f"      {'launch':<8} {'K1 A/m':>8} {'net current in an end band kA':>30} {'field at the end band T':>24} {'worst tape Ic A':>16} {'A/tape':>7} {'kg/m':>9} {'ring kt':>9} {'N':>3} {'volts per turn at 30 km/s kV':>29} {'heat J':>10} {'% of KE':>8}")
    rows4 = []
    for name, k in (("book", 1.0), ("half", 0.5), ("quarter", 0.25), ("tenth", 0.1)):
        bt = res[name][(0.012, 70.0)]["by_tau"][130.0][0]
        if bt is None:
            out(f"      {name:<8} no stator of the grid holds 100 kV with one turn"); continue
        rows4.append(dict(launch=name, **bt))
        out(f"      {name:<8} {bt['K1']:8.0f} {m.KD*bt['K_belt']*130.0/math.pi/1e3:30.1f} {bt['B_end']:24.3f} {min(bt['ic_act'], bt['ic_end'], bt['ic_leg']):16.1f} {bt['i_t']:7.1f} {bt['tape_kg_per_m']:9.1f} {bt['tape_ring_t']/1e3:9.0f} {bt['N']:3d} {bt['V_turn_f']/1e3:29.2f} {bt['E_hyst']:10.3e} {bt['share']*100:8.3f}")
    RES["sensitivity"] = dict(one_at_a_time=rows, searched=rows2, thrust_at_power=rows3, tau130=rows4)


def main():
    os.makedirs(LOGS, exist_ok=True)
    quick = "--quick" in sys.argv
    t0 = time.time()
    F_book, v_c = part1()
    part2()
    part3()
    save("checks")
    if quick:
        print(f"\n  written: {os.path.join(LOGS, 'round11_stator_checks.txt')} and .json ({time.time()-t0:.0f} s)")
        return
    del _lines[:]; RES.clear()
    res, targets, P_book = part4(F_book, v_c)
    spec = part5(res, targets, F_book, P_book)
    part6(res, spec, F_book, P_book)
    out()
    out(f"  run time {time.time()-t0:.0f} s")
    save("study")
    print(f"\n  written: {os.path.join(LOGS, 'round11_stator_study.txt')} and .json")


if __name__ == "__main__":
    main()
