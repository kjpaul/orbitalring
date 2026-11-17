#!/usr/bin/env python3
"""
Round 12, items 1 to 3: the 6 m stator at 2,500 Hz (Chapters 9, 10 and 14).
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 9.3, 9.10, 9.11, 10.1, 14.5.
Brief: Vol_III/Claude_Code_instructions_round12_2026-11-17.md, Sections 1 to 3.
Run: python tools/round12_stator.py --jobs=5     (writes runs/round12_logs/round12_stator_voltage.txt, _heat.txt, _pitch.txt and .json)
     python tools/round12_stator.py --part=1     (one part only: 1 the voltage and the current, 2 the heat, 3 the pole pitch)

The model is tools/round12_stator_model.py on tools/round11_stator_model.py.  Nothing but the logs is written.
"""
import json
import math
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # round 13: this study was run on the old 130 m stator of lsm_config.py (its launch and its
                                                   # constants); the switch keeps it reproducing now that the 6 m stator is the default
import numpy as np

import lsm_config as cfg
import lsm_physics as phys
import round11_stator_model as m
import round12_stator_model as m12

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round12_logs")
MU0 = m.MU0
F_BOOK = phys.calc_thrust(phys.calc_B_stator(cfg.I_PEAK, cfg.N_STATOR, cfg.W_COIL, cfg.G_GAP), cfg.B_SLED_NOMINAL, cfg.DELTA_MAX, cfg.A_ACTIVE_TOTAL)
P_BOOK = F_BOOK * 5000.0
SPARE_W = 11.9e3; BEARING_SITE_J = 6.8e6
ALPHAS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
LAYS = (("edge", "flat", 1.0, 1.0), ("edge", "flat", 1.0, 0.5), ("edge", "edge", 0.0, 1.0), ("flat", "flat", 1.0, 1.0))
_lines = []; RES = {}


def out(s=""):
    print(s, flush=True); _lines.append(s)


def save(name):
    with open(os.path.join(LOGS, f"round12_stator_{name}.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, f"round12_stator_{name}.json"), "w", encoding="utf-8") as fh:
        json.dump(RES, fh, indent=1, default=float)
    del _lines[:]; RES.clear()


def stator11():
    """round 11's stator of the book's launch: 6 m, on edge, end runs flat in a 1 m band set back 1 m, 3 mm tape, 70 K"""
    return m.Stator(6.0, "edge", 0.003, 70.0, 1.0, 1.0)


def lay_name(lay, end_lay, d):
    return ("on edge" if lay == "edge" else "flat") + " / " + ("on edge" if end_lay == "edge" else "flat") + (f", set back {d:.0f} m" if d > 0 else "")


# ------------------------------------------------------------------------------------------------ worker
def _job(a):
    """one stator on the book's launch: (tau, lay, end_lay, d, b_e, alpha, sled spec or None, options)"""
    tau, lay, end_lay, d, b_e, alpha, sled, opt = a
    st = m.Stator(tau, lay, opt.get("w_t", 0.003), opt.get("T", 70.0), b_e, d, opt.get("sides", 2), 0.5, end_lay)
    if sled is not None:
        m12.apply_sled(st, m12.sled_ratios(tau, sled[0], sled[1]))
    r = m12.evaluate12(st, opt.get("F", F_BOOK), opt.get("P", P_BOOK), alpha=alpha, delta=opt.get("delta", m.DELTA), phi_scale=opt.get("phi_scale", 1.0), m_sc=opt.get("m_sc", cfg.M_SPACECRAFT))
    r["sled"] = sled
    return r


def best_of(rows, key):
    ok = [r for r in rows if r["N"] >= 1]
    return min(ok, key=key) if ok else None


def tot(r, build="no copper", soldered=True):
    return m12.total(r, build, soldered)[0]


# ------------------------------------------------------------------------------------------------ part 1
def flux_check(tau, w, q_m, g_act, g_end, n_pp=4, nx_per=40, ny=200):
    """the cyclic inductance a second way: the field of the winding's filaments by Biot-Savart, integrated over the area of
    one turn (the area starts g_act and g_end inside the turn's own conductors, as the geometric mean distance does)"""
    nq = 3 * q_m
    nx = nq * nx_per
    xs = (np.arange(nx) + 0.5) / nx * tau; ys = -w / 2 + g_end + (np.arange(ny) + 0.5) / ny * (w - 2 * g_end)
    keep = (xs > g_act) & (xs < tau - g_act)
    X, Y = np.meshgrid(xs[keep], ys); P = np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], axis=1)
    Bz = np.zeros(X.size)
    for j in range(-nq * n_pp, nq * n_pp + 1):
        x0 = j * tau / nq; amp = math.cos(math.pi * j / nq)
        c = [np.array([x0, -w / 2, 0.0]), np.array([x0 + tau, -w / 2, 0.0]), np.array([x0 + tau, w / 2, 0.0]), np.array([x0, w / 2, 0.0])]
        for k in range(4):
            Bz += amp * m12.seg_field(P, c[k], c[(k + 1) % 4])[:, 2]
    return float(Bz.sum() * (tau / nx) * ((w - 2 * g_end) / ny))


def part1(pool):
    st = stator11(); tau = 6.0; W = m.W
    r = m12.evaluate12(st, F_BOOK, P_BOOK, alpha=0.2)
    N = r["N"]; I = r["I_pk"]; q = r["q"]; Kb = r["K_belt"]
    out("=" * 150)
    out("  ROUND 12, ITEM 1: THE VOLTAGE AND THE CURRENT AT 2,500 Hz")
    out("=" * 150)
    out(f"  Round 11's stator of the book's launch: pole pitch {tau:g} m, {N} turns, {r['layers']} layers of 3 mm tape on edge, coil current {I:.1f} A peak, {q:.2f} coils per phase per pole pitch per side,")
    out(f"  {r['n_coils']:,.0f} coils facing the sled, belt current density {Kb:,.0f} A/m, load angle 60 deg, above 5 km/s the sled's field as s^0.2 and the coil current as s^0.8 (s = 5 km/s / v).")
    # ---- 1.1 inductance
    out()
    out("  1.1 The coil's inductance (Neumann's integral over straight filaments, round11_stator_model.m_par; a turn is 6 m x 2 m)")
    out(f"      {'conductor':<74} {'gmd, sides / runs m':>20} {'one turn uH':>12} {f'{N} turns mH':>12} {'cyclic, per turn^2 uH':>22} {'cyclic / self':>13}")
    ind = []
    for lab, ga, ge in (("a bundle 10 mm across (the book session's hand figure: 6 m x 1.9 m, 1.9 mH)", 0.0039, 0.0039),
                        ("a bundle 50 mm across", 0.0195, 0.0195),
                        ("a bundle 200 mm across", 0.078, 0.078),
                        ("the study's: sides spread over a 0.5 m slot, runs over the 1 m band", 0.2235 * tau / 12.0, 0.2235 * 1.0)):
        Ls = float(m.m_loops(0.0, tau, W, ga, ge)); Lc, _ = m.cyclic_inductance(tau, W, 4, ga, ge)
        ind.append(dict(label=lab, g_act=ga, g_end=ge, L_turn=Ls, L_coil=Ls * N * N, L_cyclic_turn=Lc))
        out(f"      {lab:<74} {ga:9.4f} /{ge:8.4f} {Ls*1e6:12.3f} {Ls*N*N*1e3:12.3f} {Lc*1e6:22.3f} {Lc/Ls:13.2f}")
    hand = m.grover_rectangle(6.0, 1.9, 0.0039)
    out(f"      Grover's rectangle formula for 6 m x 1.9 m with a 3.9 mm mean radius: {hand*1e6:.2f} uH a turn, {hand*N*N*1e3:.2f} mH for {N} turns; omega L I at 2,500 Hz and {I:.0f} A = {2*math.pi*2500*hand*N*N*I/1e3:.1f} kV.")
    ga, ge = 0.2235 * tau / 12.0, 0.2235
    Lc = st.Lc; L_eff = N * N * Lc * q / 4.0
    out(f"      The study's coil as the supply sees it: the flux through a turn from all the coils of its side is L_cyclic x (ampere-turns of a coil), {Lc*1e6:.3f} uH per turn^2 with")
    out(f"      four coils per phase per pole pitch; a coil of {N} turns carrying 1/{q:.2f} of a phase belt's current sees {L_eff*1e3:.3f} mH.  At {I:.0f} A and 2,500 Hz that would be {2*math.pi*2500*L_eff*I/1e3:.1f} kV;")
    out(f"      the launch never has both at once (table 1.2).")
    out(f"      mutual inductance of a turn with the turns of the other coils of the side (one coil per phase per pole pitch, {tau/3:.0f} m apart), uH per turn^2:")
    mut = []
    for j in range(1, 7):
        Mj = float(m.m_loops(j * tau / 3.0, tau, W, ga, ge)); mut.append(Mj)
        out(f"        {j} places on ({j*tau/3:4.0f} m, {60*j:3d} deg): {Mj*1e6:9.4f}; times cos(phase) {Mj*1e6*math.cos(math.pi*j/3):9.4f}")
    Lc2 = flux_check(tau, W, 4, ga, ge)
    out(f"      second way: the winding's field by Biot-Savart integrated over the turn's area (function flux_check): {Lc2*1e6:.3f} uH per turn^2 against Neumann's {Lc*1e6:.3f}.")
    # the other side and the sled
    out()
    out("      the other side of the stator, if the two sides face each other across the sled (its coils carry the same current in phase): cyclic mutual / own cyclic")
    oth = []
    for D in (0.3, 0.5, 1.0, 2.0):
        Mc = m12.cyclic(tau, W, 4, ga, ge, D); oth.append(dict(D=D, M=Mc, ratio=Mc / Lc))
        out(f"        faces {D:.1f} m apart: {Mc*1e6:8.3f} uH per turn^2, {Mc/Lc:6.3f} of the own cyclic inductance (a flat sheet gives exp(-pi D / tau) = {math.exp(-math.pi*D/tau):.3f})")
    M_sf, M_fs, L_f = m12.sled_coupling(tau, W, 4, ga, ge, m.G)
    rho = M_sf * M_fs / (L_f * Lc)
    At_coil = Kb * tau / 24.0
    q1 = m12.sled_ratios(6.0, 1, 0.2)
    out(f"      the sled's coils (one racetrack turn per pole, 6 m x 2 m, {m.G*1e3:.0f} mm from the stator, conductor 50 mm mean radius): flux through a stator turn per ampere in the sled's")
    out(f"      coils {M_sf*1e6:.3f} uH; flux through a sled coil per ampere-turn of each stator coil (wave aligned with the pole) {M_fs*1e6:.3f} uH; a sled coil's own cyclic inductance {L_f*1e6:.3f} uH.")
    out(f"      If the sled's coils hold their flux (persistent mode), the part of the stator's flux that lies along the sled's poles is cut by M_sf M_fs / (L_f L_cyclic) = {rho:.3f},")
    out(f"      and the sled's coils take up {M_fs*At_coil*math.cos(m.DELTA)/L_f/1e3:.2f} kA-turns per pole at the peak (the stator's ampere-turns along the poles x M_fs / L_f) against the {q1['At_for_0p1T']/1e3:.0f} kA-turns")
    out(f"      that make 0.10 T (round12_stator_model.sled_ratios).  The stator's field across the poles (the thrust component) links no flux with the sled's coils.")
    K1 = r["K1"]
    out(f"      the stator's own field at the sled's coils, {m.G*1e3:.0f} mm away: fundamental mu0 K1 / 2 x exp(-pi g / tau) = {0.5*MU0*K1*math.exp(-math.pi*m.G/tau)*1e3:.1f} mT at the peak current (steady in the sled's frame),")
    harm = []
    for nu in (5, 7, 11, 13):
        Bn = 0.5 * MU0 * K1 / nu * math.exp(-nu * math.pi * m.G / tau); k6 = nu + 1 if nu in (5, 11) else nu - 1
        harm.append(dict(nu=nu, B=Bn, mult=k6))
    out("      and the belt harmonics of a 60 degree winding (K1 / nu), which the sled sees as alternating fields: " + "; ".join(f"{h['nu']}th {h['B']*1e3:.2f} mT at {h['mult']} f" for h in harm) + f" (f = 417 Hz at 5 km/s).")
    RES["inductance"] = dict(rows=ind, grover_hand=hand, L_cyclic_turn=Lc, L_eff=L_eff, flux_check=Lc2, mutuals=mut, other_side=oth, sled=dict(M_sf=M_sf, M_fs=M_fs, L_f=L_f, rho=rho, At_induced=M_fs * At_coil * math.cos(m.DELTA) / L_f), harmonics=harm)
    # ---- 1.2 phasors
    speeds = (1000.0, 5000.0, 10000.0, 20000.0, 30000.0)
    def ph_table(rr, label):
        out()
        out(f"  {label}")
        out(f"      {'km/s':>5} {'f Hz':>7} {'sled field mT':>13} {'coil A':>8} {'back-EMF kV':>12} {'inductive kV':>13} {'terminal kV':>12} {'if flux adds kV':>15} {'coil MW':>8} {'coil MVA':>9} {'coil MVAr':>10} {'power factor':>12} {'V ahead of I deg':>16} | {'launch GW':>9} {'GVA':>8}")
        tab = m12.phasors(rr, speeds)
        for p in tab:
            out(f"      {p['v']/1e3:5.0f} {p['f']:7.1f} {p['B0']*1e3:13.1f} {p['I']:8.1f} {p['E']/1e3:12.2f} {p['X']/1e3:13.2f} {p['V']/1e3:12.2f} {p['V_mag']/1e3:15.2f} {p['P']/1e6:8.2f} {p['S']/1e6:9.2f} {p['Q']/1e6:10.2f} {p['pf']:12.3f} {p['angle_V_I']:16.1f} | {p['P']*rr['n_coils']/1e9:9.2f} {p['S']*rr['n_coils']/1e9:8.2f}")
        out(f"      {rr['N']} turns, {rr['n_coils']:,.0f} coils.  Largest terminal voltage {rr['V_max']/1e3:.1f} kV at {rr['v_at_Vmax']/1e3:.1f} km/s; largest back-EMF {rr['E_max']/1e3:.1f} kV; largest apparent power of the launch {rr['S_max']/1e9:.1f} GVA.")
        out(f"      An inverter that must give the largest voltage and the largest current, at different times: (1/2) x {rr['V_max']/1e3:.1f} kV x {rr['I_pk']:.0f} A = {rr['rating_coil']/1e6:.1f} MVA a coil, {rr['rating']/1e9:,.0f} GVA for the launch, {rr['rating']/rr['P']:.2f} times the {rr['P']/1e9:.0f} GW.")
        return tab
    RES["phasors_round11"] = ph_table(r, "1.2 The phasor budget of round 11's stator (current 30 degrees ahead of the back-EMF, so the stator's flux lies against the sled's at the 60 degree load angle)")
    out("      The terminal voltage is |E + j omega L I| = sqrt(E^2 + X^2 - 2 E X cos 60) with X = omega L I.  The real power is (1/2) E I sin 60; the reactive power the supply gives is")
    out("      (1/2) I (X - E cos 60): the coil's own (1/2) X I less what the sled's flux returns.  'V ahead of I' is the angle between the terminal voltage and the current.")
    # ---- 1.3 the control
    out()
    out("  1.3 What the control has to hold: the angle of the current.  The same coil and current amplitude with the load angle off by an error (thrust goes as sin, the voltage follows)")
    out(f"      {'load angle deg':>14} {'thrust / design':>15} {'terminal at 5 km/s kV':>22} {'terminal at 30 km/s kV':>23}")
    ctl = []
    for dd in (30.0, 45.0, 50.0, 60.0, 70.0, 75.0, 90.0, 120.0):
        tb = m12.phasors(r, (5000.0, 30000.0), delta=math.radians(dd))
        ctl.append(dict(delta=dd, thrust=math.sin(math.radians(dd)) / math.sin(m.DELTA), V5=tb[0]["V"], V30=tb[1]["V"]))
        out(f"      {dd:14.0f} {math.sin(math.radians(dd))/math.sin(m.DELTA):15.3f} {tb[0]['V']/1e3:22.2f} {tb[1]['V']/1e3:23.2f}")
    out("      120 degrees is the same thrust with the stator's flux adding to the sled's (the current 30 degrees behind the back-EMF).")
    RES["control"] = ctl
    # the split and the inverter
    out()
    out("      The split above 5 km/s (sled's field as s^alpha, coil current as s^(1 - alpha)), the load angle, and the inductance, on the same stator geometry (the book's launch; sine-wave sled as round 11;")
    out("      heat = hysteresis + eddy in a copper-free tape + soldered layers, items 2.1 and 2.2):")
    hdr = (f"      {'case':<66} {'N':>3} {'coil A':>7} {'coils':>7} {'largest terminal kV':>19} {'largest EMF kV':>14} {'at 30 km/s: A':>13} {'inductive kV':>12} {'pf':>6} {'largest GVA':>11} {'inverter GVA':>12} {'x power':>7} {'tape kg/m':>9} {'heat J':>10}")
    out(hdr)
    jobs = [("alpha = %.1f" % a + (" (round 11)" if a == 0.2 else " (the book's rule: field lowered, current held)" if a == 1.0 else " (field held, current lowered)" if a == 0.0 else ""), (6.0, "edge", "flat", 1.0, 1.0, a, None, {})) for a in ALPHAS]
    jobs += [("alpha = 0.2, load angle 90 deg (current in phase with the back-EMF)", (6.0, "edge", "flat", 1.0, 1.0, 0.2, None, dict(delta=math.pi / 2))),
             ("alpha = 0.6, load angle 90 deg", (6.0, "edge", "flat", 1.0, 1.0, 0.6, None, dict(delta=math.pi / 2)))]
    Lc_c, _ = m.cyclic_inductance(tau, W, 4, 0.0195, 0.0195)
    jobs += [(f"alpha = 0.2, conductors bundled 50 mm across (inductance x {Lc_c/Lc:.2f})", (6.0, "edge", "flat", 1.0, 1.0, 0.2, None, dict(phi_scale=Lc_c / Lc))),
             (f"alpha = 0.2, the two sides facing 0.3 m apart (inductance x {1+oth[0]['ratio']:.2f})", (6.0, "edge", "flat", 1.0, 1.0, 0.2, None, dict(phi_scale=1 + oth[0]["ratio"]))),
             (f"alpha = 0.6, the two sides facing 0.3 m apart", (6.0, "edge", "flat", 1.0, 1.0, 0.6, None, dict(phi_scale=1 + oth[0]["ratio"])))]
    res = pool.map(_job, [j[1] for j in jobs])
    split = []
    for (lab, _), rr in zip(jobs, res):
        p30 = m12.phasors(rr, (30000.0,))[0]
        split.append(dict(label=lab, N=rr["N"], I=rr["I_pk"], coils=rr["n_coils"], V_max=rr["V_max"], E_max=rr["E_max"], I30=p30["I"], X30=p30["X"], pf30=p30["pf"], S_max=rr["S_max"], rating=rr["rating"],
                          tape=rr["tape_kg_per_m"], heat=tot(rr), alpha=rr["alpha"], delta=rr["delta_deg"]))
        out(f"      {lab:<66} {rr['N']:3d} {rr['I_pk']:7.0f} {rr['n_coils']:7,.0f} {rr['V_max']/1e3:19.1f} {rr['E_max']/1e3:14.1f} {p30['I']:13.0f} {p30['X']/1e3:12.1f} {p30['pf']:6.3f} {rr['S_max']/1e9:11.1f} {rr['rating']/1e9:12.0f} {rr['rating']/rr['P']:7.2f} {rr['tape_kg_per_m']:9.2f} {tot(rr):10.3e}")
    RES["split"] = split
    r6 = res[3]                                  # alpha = 0.6
    RES["phasors_alpha06"] = ph_table(r6, "1.2b The same budget with the split at alpha = 0.6 (the smallest inverter of the table above)")
    # ---- 1.4 cancelling
    out()
    out("  1.4 What would cancel the inductive voltage: a series capacitor per coil, C = 1 / (omega^2 L), which resonates at one frequency only")
    out(f"      {'km/s':>5} {'f Hz':>7} {'coil L mH':>10} {'C uF':>8} {'capacitor volts kV':>19} {'capacitor MVAr (1/2 X I^2)':>27} {'energy in it J':>15} {'terminal with it kV':>20} {'without kV':>11}")
    cap = []
    for p in m12.phasors(r, speeds):
        om = 2 * math.pi * p["f"]; L = p["L_eff"]; C = 1 / (om * om * L)
        cap.append(dict(v=p["v"], f=p["f"], L=L, C=C, Vc=p["X"], Q=0.5 * p["X"] * p["I"], W=0.5 * C * p["X"]**2, V_with=p["E"], V_without=p["V"]))
        out(f"      {p['v']/1e3:5.0f} {p['f']:7.1f} {L*1e3:10.3f} {C*1e6:8.2f} {p['X']/1e3:19.2f} {0.5*p['X']*p['I']/1e6:27.2f} {0.5*C*p['X']**2:15.1f} {p['E']/1e3:20.2f} {p['V']/1e3:11.2f}")
    out("      With the inductive voltage cancelled the terminal voltage is the back-EMF, which is higher than the terminal voltage of table 1.2: at the 60 degree load angle the stator's flux")
    out("      already takes voltage off.  The capacitance needed runs from the first row's to the last row's in one launch, a range of " + f"{cap[0]['C']/cap[-1]['C']:,.0f} to 1.")
    RES["capacitor"] = cap
    out()
    out("      Fewer turns and more current, or a shorter coil: the ratio inductive / back-EMF does not depend on the turns (both go as N), nor on the frequency (both go as the speed).")
    out("      It is the stator's own flux through a turn over the sled's flux through it.  At the peak current and the full sled field, by pole pitch (same lay; 3 mm tape):")
    out(f"      {'tau m':>6} {'f at 30 km/s Hz':>16} {'stator flux / sled flux at the peak':>36} {'own mean field through a coil mT':>33}")
    rat = []
    for t_ in (3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0):
        s_ = m.Stator(t_, "edge", 0.003, 70.0, 1.0, 1.0)
        rat.append(dict(tau=t_, ratio=s_.B_avg * Kb / m.B0N, B_avg=s_.B_avg * Kb))
        out(f"      {t_:6.0f} {cfg.V_LAUNCH/(2*t_):16.0f} {s_.B_avg*Kb/m.B0N:36.3f} {s_.B_avg*Kb*1e3:33.1f}")
    RES["flux_ratio"] = rat
    # ---- 1.6 the time a coil has
    out()
    out("  1.6 The time a coil has, and the energy in its field")
    out(f"      {'km/s':>5} {'time under the 10 km sled s':>27} {'cycles in that time':>19} {'coil A':>8} {'coil L mH':>10} {'time to bring the current up at 100 kV us':>42} {'share of a period':>17} {'energy in a coil J':>19} {'in all coils MJ':>15} {'in and out at the sled ends MW':>30}")
    tm = []
    for p in m12.phasors(r, speeds):
        L = p["L_eff"]; t_pass = cfg.L_SLED / p["v"]; t_up = L * p["I"] / (m.VLIM); W1 = 0.5 * L * p["I"]**2
        Wall = 0.5 * W1 * r["n_coils"]                     # the mean of i^2 over the coils of a balanced set is half the peak
        tm.append(dict(v=p["v"], t_pass=t_pass, cycles=t_pass * p["f"], I=p["I"], L=L, t_up=t_up, W_coil=W1, W_all=Wall, P_ends=Wall / cfg.L_SLED * p["v"]))
        out(f"      {p['v']/1e3:5.0f} {t_pass:27.3f} {t_pass*p['f']:19.0f} {p['I']:8.1f} {L*1e3:10.3f} {t_up*1e6:42.1f} {t_up*p['f']:17.3f} {W1:19.1f} {Wall/1e6:15.2f} {Wall/cfg.L_SLED*p['v']/1e6:30.2f}")
    out("      'time to bring the current up' is L I / 100 kV, the coil alone with the whole 100 kV across it.  The energy in a coil is (1/2) L I^2 at its peak; over the coils of a")
    out("      balanced set the mean is half that.  That energy is put in where the sled's nose arrives and taken out behind its tail: stored energy per meter x speed.")
    RES["time"] = tm
    RES["round11"] = {k: v for k, v in r.items() if k != "eddy"}
    save("voltage")


# ------------------------------------------------------------------------------------------------ part 2
def heat_row(lab, rr, width=70):
    e = rr["eddy"]; KE = rr["KE"]
    return (f"      {lab:<{width}} {lay_name(rr['lay'], rr['end_lay'], rr['d']):<26} {rr['alpha']:5.1f} {rr['layers']:4d} | {rr['E_hyst']:10.3e} {e['no copper']['E']:10.3e} {e['20 um copper']['E']:10.3e} {e['40 um copper']['E']:10.3e} {rr['E_sold']:10.3e} {rr['E_film']:10.3e} |"
            f" {tot(rr):10.3e} {tot(rr)/KE*100:8.4f} {tot(rr, '20 um copper'):10.3e} {tot(rr, '40 um copper'):10.3e} {tot(rr, 'no copper', False):10.3e}")


HEAT_HDR = (f"{'lay, active / end runs':<26} {'alpha':>5} {'lyrs':>4} | {'hysteresis':>10} {'eddy no Cu':>10} {'eddy 20 um':>10} {'eddy 40 um':>10} {'soldered':>10} {'insulated':>10} |"
            f" {'total J':>10} {'% of KE':>8} {'with 20 um':>10} {'with 40 um':>10} {'insulated':>10}")


def part2(pool):
    out("=" * 150)
    out("  ROUND 12, ITEM 2: THE HEAT, WITH EVERYTHING IN IT (the book's launch on the 6 m stator, 3 mm tape, 70 K)")
    out("=" * 150)
    st = stator11(); r = m12.evaluate12(st, F_BOOK, P_BOOK, alpha=0.2)
    out("  2.1 and 2.2 Round 11's stator on round 11's sled field (a sine wave of 0.10 T), per launch, J")
    out("      hysteresis: Brandt and Indenbom for the field through each tape's face, plus Norris (round 11).  eddy: pi^2 f^2 B^2 w^3 (t / rho) / 6 per meter of tape for the field")
    out(f"      through the face, in the tape's metal: {m.T_HAST*1e6:.0f} um of Hastelloy at {m.RHO_HAST*1e6:.2f} uohm m and {m.T_AG*1e6:.0f} um of silver at {m.RHO_AG*1e9:.0f} nohm m ({m12.G_BARE:.0f} S), and with 20 or 40 um of copper at {m12.RHO_CU*1e9:.0f} nohm m")
    out(f"      ({m12.sheet_conductance(20e-6):,.0f} and {m12.sheet_conductance(40e-6):,.0f} S).  soldered: the {r['layers']} layers of a turn soldered into one stack {r['d_stack']*1e3:.2f} mm thick, in the field along the tapes' faces: one Bean slab of")
    out("      the stack's thickness and engineering current density (the bearing's model, bearing_round5.py Section 5).  insulated: the same field on layers not soldered (the film's own slab).")
    out(f"      {'case':<70} " + HEAT_HDR)
    out(heat_row("round 11's stator, sine-wave sled", r))
    RES["round11_sine"] = r
    out(f"      the eddy loss goes as the frequency squared: its share below 5 km/s is {r['eddy']['no copper']['E_low']/r['eddy']['no copper']['E']*100:.1f} percent, against {r['E_hyst_low']/r['E_hyst']*100:.1f} percent of the hysteresis.")
    # ---- 2.3 the sled's coils
    out()
    out("  2.3 The sled's field coil (Biot-Savart, round12_stator_model.sled_ratios): the field at the stator's tapes as a multiple of B0, peak over a cycle, averaged over the tape positions")
    out("      plain: one racetrack per pole, its pack laid flat.  stepped: n nested loops of equal current per pole, cross conductors where a sine-wave sheet has carried (j - 1/2) / n of its current.")
    out(f"      {'sled coil':<34} {'kA-turns per pole':>17} {'conductor kA m per m of sled':>29} | {'along the track: peak':>21} {'fundamental':>11} {'3rd':>6} {'5th':>6} {'7th':>6} | {'through the face: peak':>22} {'3rd':>6} {'5th':>6} {'7th':>6} | {'end band at the face':>20} {'1 m back':>9} {'down a leg':>10}")
    designs = [("plain, thin conductor", 1, 0.0), ("plain, pack 0.2 m", 1, 0.2), ("plain, pack 0.5 m", 1, 0.5), ("stepped, 3 loops", 3, 0.2), ("stepped, 6 loops", 6, 0.2), ("stepped, 13 loops", 13, 0.2), ("stepped, 26 loops", 26, 0.2), ("stepped, 13 loops, sides over 0.5 m", 13, 0.5)]
    sl = []
    for lab, n, pack in designs:
        q = m12.sled_ratios(6.0, n, pack); sl.append(dict(label=lab, **q))
        out(f"      {lab:<34} {q['At_for_0p1T']/1e3:17.1f} {q['Am_per_m']/1e3:29.1f} | {q['rx_peak']:21.2f} {q['rx_fund']:11.2f} {q['hx'][0]:6.3f} {q['hx'][1]:6.3f} {q['hx'][2]:6.3f} | {q['rz_peak']:22.2f} {q['hz'][0]:6.3f} {q['hz'][1]:6.3f} {q['hz'][2]:6.3f} | {q['cz_face']:20.2f} {q['cz_back']:9.2f} {q['leg_x']:10.2f}")
    q130 = m12.sled_ratios(130.0, 1, 0.2)
    out(f"      for comparison the 130 m racetrack of the present sled (pack 0.2 m): {q130['At_for_0p1T']/1e3:.1f} kA-turns per pole, {q130['Am_per_m']/1e3:.1f} kA m of conductor per meter of sled (one side).")
    out(f"      the study's sine wave has along the track {st.rx:.2f}, through the face 1.00, end band at the face 0.50, 1 m back {st.c_run:.2f}.")
    RES["sled_designs"] = sl; RES["sled_130"] = q130
    out()
    out("      the heat of the book's launch on each sled (the 6 m stator; the lay and the split chosen again for the least total with a copper-free tape and soldered layers), J")
    out(f"      {'sled coil':<70} " + HEAT_HDR)
    jobs = []
    for lab, n, pack in designs:
        for lay, end_lay, d, be in LAYS:
            for a in ALPHAS:
                jobs.append((6.0, lay, end_lay, d, be, a, (n, pack), {}))
    res = pool.map(_job, jobs)
    per = len(LAYS) * len(ALPHAS); heat = []
    for k, (lab, n, pack) in enumerate(designs):
        b = best_of(res[k * per:(k + 1) * per], tot); heat.append(dict(label=lab, n=n, pack=pack, **b))
        out(heat_row(lab, b))
    RES["sled_heat"] = heat
    # ---- 2.4 the total on the smooth sled
    chosen = next(h for h in heat if h["n"] == 13 and h["pack"] == 0.2)
    out()
    out("  2.4 The total on the stepped sled coil of 13 loops (copper-free tape, soldered layers; 94.3 W per W at 70 K)")
    brake = _job((6.0, chosen["lay"], chosen["end_lay"], chosen["d"], chosen["b_e"], chosen["alpha"], (13, 0.2), dict(m_sc=0.0)))
    E, E_low = m12.total(chosen); Eb, _ = m12.total(brake)
    n_sites = m.N_SITES; wpw = m.work_per_watt(70.0)
    site = E / n_sites; worst = E_low / chosen["n_start"]
    rows = [("hysteresis (field through the faces, and the transport current)", chosen["E_hyst"]), ("eddy loss in the tape's metal, no copper", chosen["eddy"]["no copper"]["E"]), ("soldered layers (field along the faces)", chosen["E_sold"])]
    for lab, val in rows:
        out(f"      {lab:<66} {val:10.3e} J  {val/chosen['KE']*100:8.4f} % of the launch's energy")
    out(f"      {'total per launch':<66} {E:10.3e} J  {E/chosen['KE']*100:8.4f} %      (round 11: 3.543e+12 J, 0.0768 %)")
    out(f"      {'braking the empty sled on the same stator':<66} {Eb:10.3e} J  ({Eb/E:.2f} of the launch)")
    out(f"      mean per site per launch {site/1e6:.2f} MJ; with braking {(E+Eb)/n_sites/1e6:.2f} MJ.  Worst sites: the {chosen['n_start']:,.0f} sites of the first {chosen['x_c']/1e3:,.0f} km (below 5 km/s) take {worst/1e6:.1f} MJ each from that pass,")
    out(f"      {chosen['site_kw_c']:,.0f} kW per site while the sled is overhead at 5 km/s (2 s) and {chosen['site_kw_f']:,.0f} kW at 30 km/s (1/3 s, {chosen['site_kw_f']*cfg.L_SLED/cfg.V_LAUNCH:,.0f} kJ per pass).")
    out(f"      compressor work per launch {E*wpw:.3e} J, {E*wpw/chosen['KE']*100:.2f} % of the launch's energy.")
    out(f"      heaviest launches a day on the 11.9 kW spare: stator alone {SPARE_W*86400/site:.1f}; with the bearing's 6.8 MJ {SPARE_W*86400/(site+BEARING_SITE_J):.1f}; with braking and the bearing {SPARE_W*86400/((E+Eb)/n_sites+BEARING_SITE_J):.1f}.")
    out(f"      fourteen a day: {14*site/86400/1e3:.2f} kW per site as a mean, {14*(E+Eb)/n_sites/86400/1e3:.2f} kW with braking; compressor {14*site/86400*wpw/1e6:.3f} MW per site.")
    for build in ("20 um copper", "40 um copper"):
        Ec = m12.total(chosen, build)[0]
        out(f"      with {build}: total {Ec:.3e} J ({Ec/chosen['KE']*100:.3f} %), {Ec/n_sites/1e6:.1f} MJ per site, {SPARE_W*86400/(Ec/n_sites):.1f} launches a day on the spare (stator alone).")
    RES["chosen_sled"] = dict(design=chosen, brake=brake, E=E, E_low=E_low, E_brake=Eb, site=site, worst=worst)
    save("heat")


# ------------------------------------------------------------------------------------------------ part 3
def n_steps_for(tau):
    """loops per pole that keep the cross conductors about 0.15 m apart where they are densest (the 6 m pole's 13)"""
    return max(3, int(round(tau / (math.pi * 0.15))))


def part3(pool):
    out("=" * 150)
    out("  ROUND 12, ITEM 3: THE POLE PITCH AGAIN, WITH ALL OF THAT COUNTED (the book's launch, 3 mm tape, 70 K, copper-free tape, soldered layers, stepped sled coil)")
    out("=" * 150)
    taus = (4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0)
    out("  The sled's stepped coil at each pole pitch (loops per pole chosen to keep the cross conductors about 0.15 m apart where they are densest):")
    out(f"      {'tau m':>6} {'loops':>6} {'kA-turns per pole':>18} {'conductor kA m per m of sled':>29} {'along the track: peak':>22} {'through the face: peak':>23} {'end band 1 m back':>18}")
    sleds = {}
    for t_ in taus:
        q = m12.sled_ratios(t_, n_steps_for(t_), 0.2); sleds[t_] = q
        out(f"      {t_:6.0f} {q['n_steps']:6d} {q['At_for_0p1T']/1e3:18.1f} {q['Am_per_m']/1e3:29.1f} {q['rx_peak']:22.2f} {q['rz_peak']:23.2f} {q['cz_back']:18.2f}")
    RES["sleds"] = {str(k): v for k, v in sleds.items()}
    jobs = []
    for t_ in taus:
        for lay, end_lay, d, be in LAYS:
            for a in ALPHAS:
                jobs.append((t_, lay, end_lay, d, be, a, (n_steps_for(t_), 0.2), {}))
    res = pool.map(_job, jobs)
    per = len(LAYS) * len(ALPHAS)
    hdr = (f"      {'tau m':>6} {'lay':<26} {'alpha':>5} {'N':>3} {'lyrs':>4} {'coils':>7} {'f Hz':>6} | {'terminal kV':>11} {'EMF kV':>7} {'largest GVA':>11} {'inverter GVA':>12} {'x power':>7} | {'tape kg/m':>9} {'ring kt':>8} |"
           f" {'hysteresis':>10} {'eddy':>10} {'soldered':>10} {'total J':>10} {'% of KE':>8} {'MJ/site':>8} {'per day':>8}")
    def row(b):
        E = tot(b)
        return (f"      {b['tau']:6.0f} {lay_name(b['lay'], b['end_lay'], b['d']):<26} {b['alpha']:5.1f} {b['N']:3d} {b['layers']:4d} {b['n_coils']:7,.0f} {b['f_end']:6.0f} | {b['V_max']/1e3:11.1f} {b['E_max']/1e3:7.1f} {b['S_max']/1e9:11.1f} {b['rating']/1e9:12.0f} {b['rating']/b['P']:7.2f} |"
                f" {b['tape_kg_per_m']:9.2f} {b['tape_ring_t']/1e3:8.0f} | {b['E_hyst']:10.3e} {b['eddy']['no copper']['E']:10.3e} {b['E_sold']:10.3e} {E:10.3e} {E/b['KE']*100:8.4f} {E/m.N_SITES/1e6:8.2f} {SPARE_W*86400/(E/m.N_SITES):8.1f}")
    sweep = {"heat": [], "inverter": [], "all": []}
    out()
    out("  3.1 The stator with the least total heat at each pole pitch ('per day' is the heaviest launches a day on the 11.9 kW spare, stator alone)")
    out(hdr)
    for k, t_ in enumerate(taus):
        b = best_of(res[k * per:(k + 1) * per], tot); sweep["heat"].append(b); out(row(b))
    out()
    out("  3.2 The stator with the smallest inverter rating at each pole pitch (largest terminal voltage x largest current)")
    out(hdr)
    for k, t_ in enumerate(taus):
        b = best_of(res[k * per:(k + 1) * per], lambda r: (round(r["rating"] / 1e9), tot(r))); sweep["inverter"].append(b); out(row(b))
    out()
    out("  3.3 The split at each pole pitch, lay on edge / flat set back 1 m, band 1 m (every alpha)")
    out(hdr)
    for k, t_ in enumerate(taus):
        for r in res[k * per:k * per + len(ALPHAS)]:
            sweep["all"].append(r)
            if r["N"] >= 1:
                out(row(r))
            else:
                out(f"      {t_:6.0f} alpha {r['alpha']:.1f}: one turn does not hold 100 kV")
    RES["sweep"] = sweep
    # the same sweep with 20 um of copper, to show what the frequency-squared loss does to the choice
    out()
    out("  3.4 The same with 20 um of copper on the tape (the least total heat at each pole pitch)")
    out(f"      {'tau m':>6} {'lay':<26} {'alpha':>5} {'f Hz':>6} {'tape kg/m':>9} {'hysteresis':>10} {'eddy':>10} {'soldered':>10} {'total J':>10} {'% of KE':>8}")
    cu = []
    for k, t_ in enumerate(taus):
        b = best_of(res[k * per:(k + 1) * per], lambda r: tot(r, "20 um copper")); cu.append(b); E = tot(b, "20 um copper")
        out(f"      {t_:6.0f} {lay_name(b['lay'], b['end_lay'], b['d']):<26} {b['alpha']:5.1f} {b['f_end']:6.0f} {b['tape_kg_per_m']:9.2f} {b['E_hyst']:10.3e} {b['eddy']['20 um copper']['E']:10.3e} {b['E_sold']:10.3e} {E:10.3e} {E/b['KE']*100:8.4f}")
    RES["sweep_cu20"] = cu
    save("pitch")


def main():
    os.makedirs(LOGS, exist_ok=True)
    part = next((a.split("=")[1] for a in sys.argv if a.startswith("--part=")), "all")
    n_jobs = next((int(a.split("=")[1]) for a in sys.argv if a.startswith("--jobs=")), 4)
    t0 = time.time()
    with mp.Pool(n_jobs) as pool:
        if part in ("1", "all"):
            part1(pool)
        if part in ("2", "all"):
            part2(pool)
        if part in ("3", "all"):
            part3(pool)
    print(f"  run time {time.time()-t0:.0f} s; logs in {LOGS}")


if __name__ == "__main__":
    main()
