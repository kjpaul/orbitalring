#!/usr/bin/env python3
"""
The adopted levitation bearing: air-core support and lift tracks, their ratings, the dynamics cases, the losses, the bundle and the stray field; and the ring renumbered with the bearing's hardware mass.
Book: The Orbital Ring, Volume IV: Chapter 8 (named; Sections 8.1 to 8.3, the 114 cases and the 50 mm smallest gap).
Run: python final_bearing.py [--parts=ABCD] (reads the caches in data/)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

The chosen bearing design run through, and the ring renumbered (Paul, 2026-10-12, third task of the day).

Design (requirement set 4 of bearing_mass_2026-10-12.md, air core): two casing support tracks, each two 0.5 m loops side by
side with opposite polarity, 12 mm tape in persistent mode, 200 mm plates 1.4 m wide, each able to carry 1.5 times the casing
weight alone or its weight share plus its load share; a driven 3 mm load compensation winding on each support track (25 m
coils, staggered); two lift compensation tracks of one 1 m loop each, 3 mm tape, driven (50 m coils, staggered), 200 mm plates;
no iron on the casing, no limiters, protection by the driven windings per face and per side.

Parts (Claude_Code_instructions_final_bearing_and_renumber_2026-10-12.md):
  A  the mass loop closed and every ring number of ring_altitude.py, old beside new, casing side inside or on top of 12,000 kg/m
  B  the 87 cases of bearing_design.part_F for this design, plus the double-winding losses and the one-lift-track speed limit
  C  tape, plate and ripple losses of the air-core windings; supply energy and VA per coil
  D  rounded and graded bundles (corner field, tape current, rating); stray field at the LIM stators and at 1, 3 and 10 m

Usage: python final_bearing.py [--parts=ABCD]
Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import levitation_config as lc
import lim_config
import lim_physics
import sled_passage_simulation as sp
import sled_passage_pulsed as pz
import compensation_tracks as ct
import compensation_fields as cf
import magnetostatic_2d as ms
import bearing_design as bd
import bearing_mass as bm

MU0 = lc.MU0
G_NET = ring.G_NET
M_LIM_PLATES = 1872.0
PLATE_T, PLATE_W = 0.20, 1.4
M_PLATE = PLATE_T * PLATE_W * bm.RHO_PLATE                       # 2,204 kg/m per plate
F_LIFT_TRACK = bm.F_LIFT_PAIR / 2                                 # 196.5 kN/m per lift track (the pair shares)
F_LOAD_TRACK = bm.F_PARKED / 2                                    # 198 kN/m per load winding (the pair shares)
CRYO_FACTOR = bm.CRYO_FACTOR
L_LIFT, L_LOAD = 50.0, 25.0
N_FAM = 20
RAMP = 20e-3
TAU_HW = bd.TAU_HW
DEPLOY_DAYS_OLD, DEPLOY_ENERGY_OLD_EJ = 567.3, 18.88743500        # lim_simulation.txt
DEPLOY_DEMAND_GW, OPS_DEMAND_GW = 711.57, 44.47                   # power_deployment.txt, power_ops.txt
# LIM stator position [UNVERIFIED]: reaction plates 1.2 m wide and 0.2 m thick on the cable's lateral faces, the stator 0.2 m
# outside them, 2 m tall, centered 1 m above the cable's horizontal centerline (Chapter 8 says "offset above the centerline")
STATOR_X = ring.CABLE_SIDE / 2 + 0.2 + 0.2 + 0.1
STATOR_Y = ring.CABLE_SIDE / 2 + 1.0
Y_TRACK = ring.CABLE_SIDE / 2 - PLATE_W / 2                       # 3.42 m: the plates at the cable edges


def hr(t=""):
    pz.hr(t)


def load_air():
    with open(bm.AIR_JSON, encoding="utf-8") as f:
        recs = json.load(f)
    sup = next(r for r in recs if r["meta"]["n_loops"] == 2 and r["meta"]["W"] == 0.5 and r["meta"]["plate_t"] == 0.20 and not r["meta"]["back_iron"] and not r["meta"]["plate_shape"])
    lift = next(r for r in recs if r["meta"]["n_loops"] == 1 and r["meta"]["W"] == 1.0 and r["meta"]["plate_t"] == 0.20 and not r["meta"]["back_iron"] and not r["meta"]["plate_shape"])
    return sup, lift


def curve_at(rec, F, key):
    rows = [r for r in rec["curve"] if r.get("nonlinear")]
    Fs = [r["F"] for r in rows]; NIs = [r["NI"] for r in rows]
    NI = float(np.interp(F, Fs, NIs))
    return NI, float(np.interp(NI, NIs, [r[key] for r in rows]))


# =============================================================================
# PART A
# =============================================================================

def design_masses(weight, theva, sup, lift):
    F_s = max(1.5 * weight, weight / 2 + F_LOAD_TRACK)
    t_sup = bm.air_tape(sup, F_s, theva, 0.012)
    t_load = bm.air_tape(sup, F_LOAD_TRACK, theva, 0.003)
    t_lift = bm.air_tape(lift, F_LIFT_TRACK, theva, 0.003)
    tape = 2 * (t_sup["mass"] + t_load["mass"]) + 2 * t_lift["mass"]
    return dict(F_s=F_s, t_sup=t_sup, t_load=t_load, t_lift=t_lift, tape=tape, casing_side=tape * (1 + CRYO_FACTOR), plates=4 * M_PLATE,
                hw=M_LIM_PLATES + 4 * M_PLATE)


def ring_numbers(m_load, m_hw):
    d = ring.design(m_load=m_load, m_hw=m_hw)
    d["weight"] = m_load * d["g_net"]
    d["v_wave"] = math.sqrt(d["tension"] / d["m_cable"])
    d["v_wave_struct"] = math.sqrt(d["tension"] / d["m_struct"])
    d["a_out"] = d["v_cable"]**2 / d["r"] - d["g"]
    d["sites"] = round(d["L_ring"] / 500.0)
    d["M_struct"] = d["m_struct"] * d["L_ring"]; d["M_cable"] = d["m_cable"] * d["L_ring"]; d["M_casing"] = m_load * d["L_ring"]
    d["deploy_days"] = DEPLOY_DAYS_OLD * d["E_deploy"] / (DEPLOY_ENERGY_OLD_EJ * 1e18)
    return d


def part_A(theva, sup, lift):
    hr("PART A: THE MASS LOOP CLOSED FOR THE CHOSEN DESIGN, AND THE RING NUMBERS")
    old = ring_numbers(ring.M_LOAD_M, ring.M_HARDWARE_M)
    # on top: the casing grows by the bearing's casing side; iterate
    m = ring.M_LOAD_M
    for it in range(30):
        dm = design_masses(m * G_NET, theva, sup, lift)
        m_new = ring.M_LOAD_M + dm["casing_side"]
        if abs(m_new - m) < 0.5:
            break
        m = m_new
    top = ring_numbers(m_new, dm["hw"]); dm_top = dm
    dm_in = design_masses(ring.M_LOAD_M * G_NET, theva, sup, lift)
    inside = ring_numbers(ring.M_LOAD_M, dm_in["hw"])
    print(f"  Design masses (on top of 12,000 kg/m, casing {m_new:,.0f} kg/m, weight {top['weight']/1e3:.1f} kN/m): support track requirement {dm_top['F_s']/1e3:.0f} kN/m per track")
    print(f"  (1.5 x weight = {1.5*top['weight']/1e3:.0f}, weight share + load share = {(top['weight']/2 + F_LOAD_TRACK)/1e3:.0f}); 12 mm support winding {dm_top['t_sup']['NI']/1e3:.0f} kA per leg, "
          f"{dm_top['t_sup']['I_op']:.0f} A, {dm_top['t_sup']['mass']:.1f} kg/m; 3 mm load winding {dm_top['t_load']['NI']/1e3:.0f} kA, {dm_top['t_load']['I_op']:.0f} A, {dm_top['t_load']['mass']:.1f} kg/m;")
    print(f"  3 mm lift winding {dm_top['t_lift']['NI']/1e3:.0f} kA, {dm_top['t_lift']['I_op']:.0f} A, {dm_top['t_lift']['mass']:.1f} kg/m per track.  Tape {dm_top['tape']:.0f} kg/m, with formers, cryostat and")
    print(f"  supports at one tape mass [UNVERIFIED] {dm_top['casing_side']:.0f} kg/m on the casing.  Plates 4 x {M_PLATE:,.0f} = {dm_top['plates']:,.0f} kg/m; hardware with the LIM plates {dm_top['hw']:,.0f} kg/m.")
    print()
    print(f"  Decision: the bearing's casing side ({dm_in['casing_side']:.0f} kg/m) is taken as INSIDE the 12,000 kg/m, because that figure was set in Chapter 3 as the whole")
    print(f"  casing load including its levitation tracks (Chapter 8 describes the coils as part of the casing) and {dm_in['casing_side']:.0f} kg/m is {dm_in['casing_side']/120:.1f} percent of it,")
    print(f"  inside the precision of a round estimate; the plates are the only new hardware.  Both versions follow; 'on top' is the conservative one.")
    print()
    rows = [("Casing mass", "kg/m", lambda d, ml, mh: ml, 0), ("Cable hardware", "kg/m", lambda d, ml, mh: mh, 0),
            ("Cable structural", "kg/m", lambda d, ml, mh: d["m_struct"], 0), ("Cable with hardware", "kg/m", lambda d, ml, mh: d["m_cable"], 0),
            ("Ring per meter", "kg/m", lambda d, ml, mh: d["m_ring"], 0), ("Cable area", "m^2", lambda d, ml, mh: d["area"], 2),
            ("Cable side", "m", lambda d, ml, mh: d["side"], 3), ("Cable tension", "GN", lambda d, ml, mh: d["tension"] / 1e9, 1),
            ("Cable speed", "m/s", lambda d, ml, mh: d["v_cable"], 1), ("Cable-casing relative speed", "m/s", lambda d, ml, mh: d["v_rel"], 1),
            ("Casing weight", "kN/m", lambda d, ml, mh: d["weight"] / 1e3, 2), ("Net g on the casing", "m/s^2", lambda d, ml, mh: d["g_net"], 4),
            ("Cable net outward accel.", "m/s^2", lambda d, ml, mh: d["a_out"], 4),
            ("Transverse wave speed, full mass", "m/s", lambda d, ml, mh: d["v_wave"], 0), ("Transverse wave speed, structural", "m/s", lambda d, ml, mh: d["v_wave_struct"], 0),
            ("Cable mass, structural", "Mt", lambda d, ml, mh: d["M_struct"] / 1e9, 2), ("Cable mass, with hardware", "Mt", lambda d, ml, mh: d["M_cable"] / 1e9, 2),
            ("Casing mass", "Mt", lambda d, ml, mh: d["M_casing"] / 1e9, 3), ("Ring mass", "Mt", lambda d, ml, mh: d["M_ring"] / 1e9, 2),
            ("Deployment energy", "EJ", lambda d, ml, mh: d["E_deploy"] / 1e18, 3), ("Deployment time at the same power", "days", lambda d, ml, mh: d["deploy_days"], 1),
            ("Deployment demand", "GW", lambda d, ml, mh: DEPLOY_DEMAND_GW, 2), ("Operating demand", "GW", lambda d, ml, mh: OPS_DEMAND_GW, 2),
            ("LIM sites at 500 m", "", lambda d, ml, mh: d["sites"], 0)]
    print(f"  {'quantity':>36} {'unit':>5} {'book':>14} {'inside (12,000)':>16} {'on top':>14}")
    for name, unit, fn, nd in rows:
        a, b, c = fn(old, ring.M_LOAD_M, ring.M_HARDWARE_M), fn(inside, ring.M_LOAD_M, dm_in["hw"]), fn(top, m_new, dm_top["hw"])
        print(f"  {name:>36} {unit:>5} {a:14,.{nd}f} {b:16,.{nd}f} {c:14,.{nd}f}")
    print(f"  Deployment time scales the 567.3 days of lim_simulation.txt with the deployment energy at the same LIM power (711.57 GW demand, power_deployment.txt);")
    print(f"  the operating demand (44.47 GW, power_ops.txt) is kept: it is LIM station-keeping, cryogenics and losses that do not depend on the cable mass [UNVERIFIED],")
    print(f"  plus the bearing's own losses of Part C.  Cable centrifugal and net outward accelerations, the cable speed and the relative speed move with the cable mass.")
    print()
    print(f"  What must change in ring_altitude.py (not changed; this run used overrides through ring.design()):")
    print(f"    M_HARDWARE_M = 4_076.0  ->  {dm_in['hw']:,.0f}  (LIM plates 1,872 + four 200 mm plates 1.4 m wide, {M_PLATE:,.0f} kg/m each)")
    print(f"    M_LOAD_M = 12_000.0     ->  unchanged (inside) or {m_new:,.0f} (on top)")
    print(f"  Downstream scripts that import ring_altitude and whose logs then need a rerun (run_700.sh): lsm_simulation, lim_simulation, levitation_simulation,")
    print(f"  md_simulation, anchor_line_analysis, j2_simulation, debris_simulation, tether_stress, climber_descent, power_deployment, power_ops, ring_altitude")
    print(f"  (table), sled_passage_simulation, sled_passage_limits, sled_passage_pulsed, sled_passage_spread, compensation_tracks, bearing_design, bearing_mass,")
    print(f"  final_bearing, and the Chapter 3 and 8 numbers (Table 8.1, Equations 3.x).  levitation_config.py also carries the plate thickness (0.10 m) and the")
    print(f"  40-coil array that the design no longer has.")
    return dict(old=old, inside=inside, top=top, dm_in=dm_in, dm_top=dm_top, m_top=m_new)


# =============================================================================
# PART B
# =============================================================================

def set_ring(m_casing, d):
    """Point the dynamics modules at the new casing and cable masses (module constants read at call time)."""
    pz.MU = m_casing; pz.W = m_casing * G_NET; pz.F_STEADY = pz.W; pz.MU_C = d["m_cable"]; pz.MATCH_RING = d["m_cable"] / (m_casing + d["m_cable"])
    bd.MU = m_casing; bd.W = pz.W; bd.MU_C = d["m_cable"]
    ct.MU = m_casing; ct.W = pz.W; ct.MU_C = d["m_cable"]
    sp.T_CABLE = d["tension"]


def part_B(A, theva, sup, lift, tabs):
    hr("PART B: THE CASES FOR THIS DESIGN")
    d = A["top"]; m_casing = A["m_top"]; dm = A["dm_top"]
    set_ring(m_casing, d)
    F_s = dm["F_s"]
    print(f"  Casing {m_casing:,.0f} kg/m (the on-top version, the conservative one), weight {pz.W/1e3:.1f} kN/m, cable {d['m_cable']:,.0f} kg/m, tension {d['tension']/1e9:.0f} GN.")
    print(f"  Caps: support one track {F_s/1e3:.0f} kN/m (12 mm persistent winding at half Ic); load winding {F_LOAD_TRACK/1e3:.0f} kN/m per track ({bm.F_PARKED/1e3:.0f} pair); lift {F_LIFT_TRACK/1e3:.0f} per track")
    print(f"  ({bm.F_LIFT_PAIR/1e3:.0f} pair).  Lift coils {L_LIFT:.0f} m, load coils {L_LOAD:.0f} m, {N_FAM} staggered families; planned pull ramped {RAMP*1e3:.0f} ms; protection 10 mm, 1 ms, 2 ms rise.")
    C = dict(F_one_cap=F_s, rec=dict(meta=dict(W=1.0)))
    D = dict(caps=dict(lift1=F_LIFT_TRACK, load1=F_LOAD_TRACK, lift2=2 * F_LIFT_TRACK, load2=2 * F_LOAD_TRACK),
             lift=dict(L_coil=L_LIFT, ramp=RAMP, W_track=1.0), load=dict(L_coil=L_LOAD, W_track=-0.4))
    gaps = bd.part_F(C, D, tabs)
    # ---- added cases
    print()
    print("  B.5  Added cases (point model, 30 km/s unless at rest).")
    u30 = 30000.0
    w_hd, w_hu = 2 * F_LIFT_TRACK / bd.P12, 2 * F_LOAD_TRACK / bd.P12

    def lane(u, **kw):
        L = L_LOAD if u < 7000 else L_LIFT
        dd = dict(L_sc=bd.L_SC, tau=TAU_HW, L_section=L / N_FAM, gap_dep=False, steady_loop_hz=1.0, T_spread=max(0.05, bd.L_SC / u), B_sat=1.2, n_coils=N_FAM,
                  L_coil=L, hd_width=w_hd, hu_width=w_hu, cmd_ramp=RAMP, steady_cap=2 * F_s)
        dd.update(kw)
        return pz.lane(**dd)

    def P(cap_up=2 * F_LOAD_TRACK, cap_down=2 * F_LIFT_TRACK):
        return ct.prot(0.010, 1e-3, cap_up=cap_up, cap_down=cap_down)

    cases = [("1. one support track loses both windings at rest", dict(no_sled=True, cancel=False, support_fail_t=0.0, support_fail_frac=0.5, steady_cap=F_s, prot=P(cap_up=F_LOAD_TRACK))),
             ("1. one support track loses both windings as the zone arrives", dict(support_fail_t=0.0, support_fail_frac=0.5, steady_cap=F_s, prot=P(cap_up=F_LOAD_TRACK))),
             ("1. one support track loses both windings 100 ms before", dict(support_fail_t=-0.1, support_fail_frac=0.5, steady_cap=F_s, prot=P(cap_up=F_LOAD_TRACK))),
             ("2. both support tracks lose both windings at rest", dict(no_sled=True, cancel=False, support_fail_t=0.0, support_fail_frac=1.0, prot=P(cap_up=0.0))),
             ("2. both support tracks lose both windings during the pass", dict(support_fail_t=0.0, support_fail_frac=1.0, prot=P(cap_up=0.0))),
             ("3. one lift track lost at 30 km/s (the other alone, capped)", dict(hd_width=F_LIFT_TRACK / bd.P12, prot=P(cap_down=F_LIFT_TRACK)))]
    res, hist = pz.point_pulsed([lane(u30, label=lab, **kw) for lab, kw in cases], u30, tabs, t_settle=3.0)
    print(f"  {'case':>66} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'contact':>7} {'at ms':>6}")
    extra = []
    for i, ((lab, kw), r) in enumerate(zip(cases, res)):
        t_c = float("nan")
        if r["contact"]:
            rel = hist["w"][i] - hist["wc"][i]
            j = int(np.argmax(np.abs(rel) >= bd.GAP0 - 1e-6))
            t_c = hist["t"][j] - (kw.get("support_fail_t") or 0.0)
        extra.append((lab, "point", r["min_gap_b"], r["min_gap_t"], r["contact"]))
        print(f"  {lab:>66} {r['max_exc_all']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {('yes' if r['contact'] else 'no'):>7} {t_c*1e3:6.0f}")
    print(f"  Case 2: with every up-pulling winding gone nothing holds the casing; it falls the 100 mm in about {math.sqrt(2*bd.GAP0/G_NET)*1e3:.0f} ms from rest (free fall at {G_NET:.2f} m/s^2).")
    print(f"  The four up-pulling windings are two persistent 12 mm support windings and two driven 3 mm load windings, one of each on each support track.")
    print(f"  For no single fault to take all four: each winding in its own cryostat with its own cooling loop (a cryostat breach or a cooling loss then")
    print(f"  takes one winding); its own flux pump or supply and its own energy store (the persistent windings ride through a ring power loss, the")
    print(f"  driven ones on their stores); its own gap sensors and control channel per face and per side; and no shared structure or feed that a")
    print(f"  single sled strike or quench can remove from both tracks.  A third independent winding (a second persistent 12 mm winding in a third")
    print(f"  cryostat on each track) would make the loss of any two windings survivable; the tape for it is another {dm['t_sup']['mass']:.0f} kg/m per track.")
    # the fastest launch with one lift track lost
    print()
    print(f"  Case 3: speed sweep, one lift track alone ({F_LIFT_TRACK/1e3:.0f} kN/m cap), 5,000 t; the fastest launch that keeps the casing within 20 mm:")
    us = [float(u) for u in range(10000, 30001, 1000)]
    lanes = [lane(u, hd_width=F_LIFT_TRACK / bd.P12, prot=P(cap_down=F_LIFT_TRACK), label=f"{u}") for u in us]
    excs = []
    for u, ln in zip(us, lanes):
        r, _ = pz.point_pulsed([ln], u, tabs)
        excs.append(r[0]["max_exc_all"])
        extra.append((f"3. one lift track alone at {u/1e3:.0f} km/s", "point", r[0]["min_gap_b"], r[0]["min_gap_t"], r[0]["contact"]))
    print("    " + "  ".join(f"{u/1e3:.0f}: {e*1e3:.1f}" for u, e in zip(us, excs)) + "  (km/s: mm)")
    ok = [u for u, e in zip(us, excs) if e <= 0.020]
    u_max = max(ok) if ok else float("nan")
    print(f"    fastest launch inside 20 mm with one lift track out: {u_max/1e3:.0f} km/s; at 30 km/s the casing moves {excs[-1]*1e3:.1f} mm (gap {min(extra[-1][2], extra[-1][3])*1e3:.1f} mm).")
    all_gaps = gaps + extra
    single = [g for g in all_gaps if "no per-side" not in g[0] and "no protection" not in g[0] and "support lost too" not in g[0] and not g[0].startswith("2.")]
    contacts = [g for g in all_gaps if g[4]]
    ws = min(single, key=lambda g: min(g[2], g[3]))
    print()
    print(f"  B.6  {len(all_gaps)} cases in all.  Contact in: " + ("; ".join(g[0] for g in contacts) if contacts else "none") + ".")
    print(f"  Smallest gap in any single failure: {min(ws[2], ws[3])*1e3:.1f} mm in '{ws[0]}'.  Case 2 (both tracks, both windings) is a double failure by design")
    print(f"  intent and the one case that touches; it is what the separation above is for.")
    return all_gaps


# =============================================================================
# PART C
# =============================================================================

def part_C(A, theva, sup, lift):
    hr("PART C: LOSSES FOR THE AIR-CORE WINDINGS")
    dm = A["dm_top"]; weight = A["top"]["weight"]
    cop = lim_physics.calc_cryo_cop(70.0, 400.0, lim_config.CRYO_EFF)
    t_p, v_p, P_p = pz.lsm_profile(); dt_p = np.diff(t_p)
    KE = 0.5 * (bd.M_SC0 + bd.M_SLED_M * bd.L_SLED) * 30000**2
    pull = bm.PULL_CURVE
    us = sorted(pull)
    w3 = 0.003; alpha = lim_config.ALPHA_TAPE
    # C.1 tape losses
    print("  C.1  Tape loss per pulse and per launch (3 mm tape, half Ic at the winding field, the pair sharing: lift tracks above 7 km/s on the 1 m")
    print("  loop, load windings below 7 km/s on the two 0.5 m loops; the winding field at each speed from the air-core curves at that speed's pull).")
    print(f"  Strip estimate 1/2 B_perp w Ic per meter of tape per pulse; book formula q_hts_loss_factor(I, N, w_coil, w_tape, alpha) with w_coil the loop width,")
    print(f"  N the turns of one coil (20 families), alpha {lim_config.ALPHA_PENETRATION_DEG:.0f} degrees; Norris transport loss for reference.  COP {cop:.4f}.")
    res = {}
    for name in ("strip", "book", "norris"):
        E = 0.0; q30 = 0.0; q1 = 0.0
        for i, dtt in enumerate(dt_p):
            vv = v_p[i + 1]
            Fu = float(np.interp(vv, us, [pull[u] for u in us])) / 2
            rec = lift if vv > 7000 else sup
            tp = bm.air_tape(rec, Fu, theva, w3) if Fu > 1e3 else None
            if tp is None:
                continue
            NI, By = curve_at(rec, Fu, "By_wm")
            W_loop = 1.0 if vv > 7000 else 0.5
            N_coil = (NI / tp["I_op"]) / N_FAM
            if name == "strip":
                q = 0.5 * By * w3 * tp["Ic"]
            elif name == "book":
                q = lim_physics.q_hts_loss_factor(tp["I_op"], N_coil, W_loop, w3, alpha)
            else:
                q = 0.5 * lim_physics.q_hysteresis_norris(tp["I_op"], tp["Ic"])
            q_m = q * tp["tape"] * 2
            if abs(vv - 30000) < 300 and q30 == 0.0:
                q30 = q_m
            if abs(vv - 1000) < 100 and q1 == 0.0:
                q1 = q_m
            E += q_m * vv * dtt
        res[name] = dict(E=E, q30=q30, q1=q1)
        print(f"    {name:>7}: {q1:7.1f} J per pulse per m of ring at 1 km/s (load pair), {q30:7.1f} at 30 km/s (lift pair); per launch {E/1e9:7.2f} GJ of heat, "
              f"{E/cop/1e9:8.1f} GJ electrical ({E/cop/KE*100:.3f} % of the launch energy, {E/cop/t_p[-1]/1e9:.2f} GW over the launch)")
    tl = dm["t_lift"]; ts = dm["t_load"]
    print(f"  At the design points: lift 1 m loop at {F_LIFT_TRACK/1e3:.0f} kN/m, {tl['NI']/1e3:.0f} kA per leg, winding {tl['Bw']:.2f} T, Ic {tl['Ic']:.0f} A, {tl['I_op']:.0f} A, {tl['tape']:,.0f} m of tape per m;")
    print(f"  load winding on the two 0.5 m loops at {F_LOAD_TRACK/1e3:.0f} kN/m, {ts['NI']/1e3:.0f} kA per leg, winding {ts['Bw']:.2f} T, Ic {ts['Ic']:.0f} A, {ts['I_op']:.0f} A, {ts['tape']:,.0f} m per m.")
    # C.2 plate losses
    print()
    print("  C.2  Plate losses per pass and per launch (200 mm plates, 1.4 m wide, hysteresis 50 J/m^3 per cycle at 1.2 T scaled as (B/1.2)^1.6 [UNVERIFIED],")
    print("  eddy by 1-D diffusion; plate field = the loop's peak field at that speed's pull; hold = energized length over u + V_rel).")
    mu_r, rho_p = lc.FERRO_MU_R, lc.FERRO_RESISTIVITY
    per_u = {}
    print(f"  {'u km/s':>7} {'pull per track kN/m':>19} {'plate peak T':>12} {'hold ms':>7} {'hyst J/m':>8} {'eddy J/m':>8} {'sum J/m':>8}")
    for u in us:
        Fu = pull[u] / 2
        rec = lift if u > 7000 else sup
        NI, Bp = curve_at(rec, Fu, "B_peak")
        L = L_LIFT if u > 7000 else L_LOAD
        L_on = u * max(0.05, bd.L_SC / u) + bd.L_SC + 2 * L + 2 * u * TAU_HW
        v_sw = u + bd.V_REL
        hyst = lc.FERRO_HYSTERESIS_COEFF * (Bp / 1.2) ** 1.6 * PLATE_T * PLATE_W * 2
        eddy = bd.eddy_cycle(Bp, bd.GAP0 / v_sw, L_on / v_sw, d_p=PLATE_T, mu_r=mu_r, rho_p=rho_p) * PLATE_W * 2
        per_u[u] = hyst + eddy
        print(f"  {u/1e3:7.0f} {Fu/1e3:19.0f} {Bp:12.2f} {L_on/v_sw*1e3:7.1f} {hyst:8.1f} {eddy:8.1f} {hyst+eddy:8.1f}")
    E_pl = 0.0
    for i, dtt in enumerate(dt_p):
        vv = v_p[i + 1]
        E_pl += float(np.interp(vv, us, [per_u[u] for u in us])) * (vv + bd.V_REL) * dtt
    print(f"  Per launch {E_pl/1e9:.2f} GJ in the plates ({E_pl/KE*100:.4f} % of the launch energy); the plate's peak field is at the legs, a third of the plate width.")
    # C.3 ripple
    print()
    print("  C.3  Ripple heating from the coil ends, support tracks (two 0.5 m loops, 12 mm winding at the weight share).")
    NI_w, Bp_w = curve_at(sup, weight / 2, "B_peak")
    print(f"  Model: a coil end in 3-D is beyond the 2-D solver.  Approximation: the field of a straight leg ends like that of a semi-infinite wire, which at")
    print(f"  the end plane is half the infinite-wire field and settles over about twice the distance from the leg to the plate [UNVERIFIED]; so each coil end")
    print(f"  (and each coil start) steps the plate field by 1/n of the loop's field over a length of {2*(bd.GAP0 + 0.05)*1e3:.0f} mm, with n coils per bundle.  The eddy")
    print(f"  energy of one step is half of a rise-and-fall cycle of that amplitude from the 1-D diffusion.  Steps per meter = 2 ends x 2 loops x n / L.")
    print(f"  Operating point: {weight/2e3:.0f} kN/m per track at {NI_w/1e3:.0f} kA per leg, plate peak field {Bp_w:.3f} T.")
    print(f"  {'coil length':>11} {'coils per bundle':>16} {'step T':>7} {'steps per m':>11} {'W/m, two tracks':>15} {'ring MW':>8} {'LIM power MW':>12}")
    ell = 2 * (bd.GAP0 + 0.05)
    for Lc in (1000.0, 2000.0, 5000.0, 10000.0, 20000.0):
        for n in (10, 20):
            dB = Bp_w / n
            E_step = 0.5 * bd.eddy_cycle(dB, ell / bd.V_REL, 0.0, d_p=PLATE_T, mu_r=mu_r, rho_p=rho_p)
            steps = 2 * 2 * n / Lc
            P = 2 * E_step * (0.5 + 0.4) * bd.V_REL * steps
            print(f"  {Lc/1e3:9.0f} km {n:16d} {dB:7.3f} {steps:11.3f} {P:15.2f} {P*bd.L_RING/1e6:8.1f} {P*bd.L_RING/1e6:12.1f}")
    print(f"  The drag this heating stands for is paid by the LIMs; the LIM power column equals the ring total.  Against the old limiter plates' 73.6 W/m,")
    print(f"  the 1 km coils with 10 per bundle are a tenth of that.")
    # C.4 supplies
    print()
    print("  C.4  Supply energy and VA per coil for the driven windings (stored energy from the solve at the operating ampere-turns; VA = 2 E / rise).")
    for tag, rec, F, L in (("lift", lift, F_LIFT_TRACK, L_LIFT), ("load", sup, F_LOAD_TRACK, L_LOAD)):
        NI, Wm = curve_at(rec, F, "W_mag")
        E = Wm * L / N_FAM
        tp = dm["t_lift"] if tag == "lift" else dm["t_load"]
        print(f"    {tag}: stored {Wm/1e3:.1f} kJ per m of track at {F/1e3:.0f} kN/m ({Wm/(F*bd.GAP0):.1f} x F gap), {E/1e6:.3f} MJ per {L:.0f} m coil; VA {2*E/TAU_HW/1e6:.0f} MVA at the 2 ms protection rise, "
              f"{2*E/(RAMP+TAU_HW)/1e6:.1f} MVA at the {RAMP*1e3:.0f} ms planned ramp; {tp['I_op']:.0f} A per tape: {2*E/TAU_HW/tp['I_op']/1e3:.0f} kV at 2 ms in one tape, "
              f"{2*E/TAU_HW/tp['I_op']/100/1e3:.1f} kV with 100 tapes in parallel; supplies per track {N_FAM*bd.L_RING/L/1e6:.1f} million")


# =============================================================================
# PART D
# =============================================================================

def build_shaped(shape, W, n_loops, NI, h=0.005, plate_t=PLATE_T, split=0.0):
    """One or two loops with a given bundle shape: 'square' 100 mm, 'round' (same area), 'graded' (round, current density
    falling as 1 - (r/R)^2, R = 80 mm), 'square 141' (double area)."""
    legs = []
    x0 = -n_loops * W / 2
    for k in range(n_loops):
        sgn = 1 if k % 2 == 0 else -1
        off = (k - (n_loops - 1) / 2) * split          # pull the loops apart so touching legs become two bundles
        legs.append((x0 + k * W + off, sgn * NI)); legs.append((x0 + (k + 1) * W + off, -sgn * NI))
    merged = {}
    for x, a in legs:
        merged[round(x, 6)] = merged.get(round(x, 6), 0.0) + a
    legs = sorted(merged.items())
    plate_w = n_loops * W + 0.4 + (n_loops - 1) * split
    g = ms.Grid(max(4.0, plate_w + 2.6), 3.0, h)
    g.add_iron(-plate_w / 2, plate_w / 2, bd.GAP0, bd.GAP0 + plate_t, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
    boxes = []
    for x, a in legs:
        if shape == "square":
            g.add_current(x - 0.05, x + 0.05, -0.10, 0.0, a); boxes.append((x - 0.05, x + 0.05, -0.10, 0.0))
        elif shape == "square 141":
            g.add_current(x - 0.0707, x + 0.0707, -0.1414, 0.0, a); boxes.append((x - 0.0707, x + 0.0707, -0.1414, 0.0))
        else:
            R = 0.0564 if shape == "round" else 0.080
            yc = -R
            rr = np.sqrt((g.X - x)**2 + (g.Y - yc)**2)
            m = rr <= R
            wgt = np.ones_like(g.X) if shape == "round" else np.clip(1 - (rr / R)**2, 0, None)
            wgt = wgt * m
            g.J[m] += a * wgt[m] / (wgt.sum() * g.h * g.h)
            boxes.append((x - R, x + R, -2 * R, 0.0))
    meta = dict(W=W, n_loops=n_loops, plate_t=plate_t, plate_shape=None, back_iron=0.0, bundle=0.10, gap=bd.GAP0, plate_w=plate_w, legs=legs)
    return g, boxes, meta


def part_D(A, theva, sup, lift):
    hr("PART D: BUNDLE SHAPE AND STRAY FIELD")
    dm = A["dm_top"]
    print("  D.1  Bundle shape, 5 mm cells, saturating plate, at the design ampere-turns.  Corner = largest cell field in the bundle; the tape current the")
    print("  corner allows is half of Ic at that field (THEVA 70 K); the rating change is the force at the same ampere-turns against the square bundle.")
    print(f"  {'track':>5} {'bundle':>11} {'NI kA':>5} {'force kN/m':>10} {'mean T':>6} {'corner T':>8} {'Ic at corner':>12} {'allowed A':>9} {'design A':>8} {'force vs square':>15}")
    shapes = ("square", "round", "graded", "square 141", "square, legs split 100 mm")
    print(f"  {'':>5} {'':>26} {'':>5} {'':>10} {'':>6} {'':>8} {'':>12} {'':>9} {'':>8} {'':>15} {'tape at fixed force':>19} {'force at fixed tape':>19}")
    for tag, rec, n_loops, W, F, w_tape, tp in (("support", sup, 2, 0.5, dm["F_s"], 0.012, dm["t_sup"]), ("lift", lift, 1, 1.0, F_LIFT_TRACK, 0.003, dm["t_lift"])):
        NI, _ = curve_at(rec, F, "F")
        F_sq = None
        for shape in shapes:
            if shape.startswith("square, legs") and n_loops == 1:
                continue
            g, boxes, meta = build_shaped(shape.split(",")[0], W, n_loops, NI, split=(0.1 if "split" in shape else 0.0))
            box = (-meta["plate_w"] / 2 - 0.1, meta["plate_w"] / 2 + 0.1, bd.GAP0 - 2 * g.h, bd.GAP0 + PLATE_T + 2 * g.h)
            g.solve_nonlinear(force_box=box)
            r = bm.evaluate_air(g, boxes, meta, NI)
            Bmax = 0.0; Bmean = 0.0; n = 0
            for b in boxes:
                m = g.box(*b) & (g.J != 0)
                Bm = np.sqrt(g.Bx[m]**2 + g.By[m]**2)
                Bmax = max(Bmax, float(Bm.max())); Bmean += float(Bm.sum()); n += int(m.sum())
            Bmean /= max(n, 1)
            Ic_c = bd.ic_tape(theva, Bmax, w_tape)
            if F_sq is None:
                F_sq = r["F"]
            I_allow = 0.5 * Ic_c
            print(f"  {tag:>5} {shape:>26} {NI/1e3:5.0f} {r['F']/1e3:10.1f} {Bmean:6.2f} {Bmax:8.2f} {Ic_c:12.0f} {I_allow:9.0f} {tp['I_op']:8.0f} {(r['F']/F_sq-1)*100:+14.1f} % "
                  f"{tp['I_op']/I_allow*100:17.0f} % {r['F']*min(1.0, I_allow/tp['I_op'])**2/1e3:16.0f} kN/m")
    print("  'mean' is over the conductor cells; the design current was set at half Ic at the mean field of the 10 mm square-bundle solve.  A series winding")
    print("  carries one current in every turn, so the worst turn (the corner) sets it: 'tape at fixed force' is the tape length relative to the design")
    print("  (the ampere-turns kept, the current cut to what the corner allows), 'force at fixed tape' the pull if the turns are kept and the current cut.")
    print("  The field at a bundle's surface is about mu0 x (its ampere-turns) / (its perimeter): the shared center leg of the two-loop support track carries")
    print("  2 NI and sets the worst field; splitting it into two bundles 100 mm apart halves that.  The 'graded' bundle (current density falling to zero")
    print("  at an 80 mm radius) and the larger square move the current's center away from the plate, which costs pull.  The tape's face orientation")
    print("  matters as much: in a bundle the field circles the bundle axis, so tapes wound as nested cylinders around that axis see it along their face,")
    print("  not through it [QUESTION FOR AUTHOR: the bundle's winding pattern].")
    # D.2 stray field: all four tracks in one coarse grid
    print()
    print(f"  D.2  Stray field with all four tracks in one 30 x 30 m box (50 mm cells, linear).  Cable {ring.CABLE_SIDE:.2f} m square from y = 0 (underside) to {ring.CABLE_SIDE:.2f} m;")
    print(f"  support loops under the plates at x = +-{Y_TRACK:.2f} m, lift loops above the top plates; LIM stators on the lateral faces at x = +-{STATOR_X:.2f} m,")
    print(f"  2 m tall, centered {STATOR_Y - ring.CABLE_SIDE/2:.1f} m above the cable's centerline [UNVERIFIED position].  Normal = support at the weight share, lift off;")
    print(f"  rating = support at {dm['F_s']/1e3:.0f} kN/m with the load winding at {F_LOAD_TRACK/1e3:.0f}, lift at {F_LIFT_TRACK/1e3:.0f} kN/m per track.")
    side = ring.CABLE_SIDE
    for mode in ("normal", "rating"):
        if mode == "normal":
            NI_s, _ = curve_at(sup, A["top"]["weight"] / 2, "F"); NI_l = 0.0
        else:
            NI_s, _ = curve_at(sup, dm["F_s"] + F_LOAD_TRACK, "F"); NI_l, _ = curve_at(lift, F_LIFT_TRACK, "F")
        h = 0.05
        g = ms.Grid(30.0, 30.0, h, x0=-15.0, y0=-10.0)
        for sx in (-1, 1):
            xc = sx * Y_TRACK
            g.add_iron(xc - 0.7, xc + 0.7, 0.0, PLATE_T, cf.PLATE_MU)                      # underside plate
            g.add_iron(xc - 0.7, xc + 0.7, side - PLATE_T, side, cf.PLATE_MU)              # top plate
            # support: two 0.5 m loops, legs at xc - 0.5, xc, xc + 0.5 carrying +NI, -2NI, +NI
            for dx, a in ((-0.5, NI_s), (0.0, -2 * NI_s), (0.5, NI_s)):
                g.add_current(xc + dx - 0.05, xc + dx + 0.05, -PLATE_T - 0.1, -PLATE_T, a)
            if NI_l:
                for dx, a in ((-0.5, NI_l), (0.5, -NI_l)):
                    g.add_current(xc + dx - 0.05, xc + dx + 0.05, side + 0.1, side + 0.2, a)
        g.solve()
        pts = {"LIM stator, center": (STATOR_X, STATOR_Y), "LIM stator, bottom": (STATOR_X, STATOR_Y - 1.0), "LIM stator, top": (STATOR_X, STATOR_Y + 1.0)}
        for d_ in (1.0, 3.0, 10.0):
            pts[f"{d_:.0f} m below the ring"] = (Y_TRACK, -PLATE_T - 0.1 - d_)
            pts[f"{d_:.0f} m beside the ring"] = (side / 2 + d_, side / 2)
            pts[f"{d_:.0f} m above the ring"] = (Y_TRACK, side + 0.2 + d_)
        print(f"  {mode}: support {NI_s/1e3:.0f} kA per leg, lift {NI_l/1e3:.0f} kA per leg")
        for k, (x, y) in pts.items():
            bx, by = g.field_at(x, y)
            print(f"    {k:>24}: {math.hypot(bx, by)*1e3:8.3f} mT")
        m_st = g.box(STATOR_X - 0.1, STATOR_X + 0.1, STATOR_Y - 1.0, STATOR_Y + 1.0)
        print(f"    {'largest over the stator':>24}: {float(np.sqrt(g.Bx[m_st]**2 + g.By[m_st]**2).max())*1e3:8.3f} mT")


# =============================================================================

def main():
    parts = "ABCD"
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
    t_all = time.time()
    hr("THE CHOSEN BEARING DESIGN, AND THE RING RENUMBERED (2026-10-12)")
    theva = bd.read_theva()
    sup, lift = load_air()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    A = part_A(theva, sup, lift)
    if "B" in parts:
        part_B(A, theva, sup, lift, tabs)
    if "C" in parts:
        part_C(A, theva, sup, lift)
    if "D" in parts:
        part_D(A, theva, sup, lift)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
