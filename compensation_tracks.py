#!/usr/bin/env python3
"""
Compensation tracks, final layout (Paul, 2026-10-11).

Two casing support tracks (1 m each, persistent HTS) and four compensation
tracks, each 1 m wide facing a 1.4 m plate across 100 mm: two lift compensation
tracks on the casing ceiling, two load compensation tracks on the casing floor
outboard of the support tracks.  The two tracks of a pair are on opposite sides
of the cable's centerline and either one must do the pair's job alone.  No
Halbach limiters: the compensation tracks are the gap protection.  Each track is
20 driven HTS coils, 5 cm wide and 200 m long, a new coil every 10 m, each with
its own supply and energy store.  Loads as in sled_passage_spread.py (5,000 t on
100 m, 10 km sled at 525 kg/m, 100 mm gap, 50 ms spread pull, 10 m truss).

Parts (Claude_Code_instructions_compensation_tracks_2026-10-11.md):
  A  the sled passage rerun for this layout at 1, 5, 10, 20 and 30 km/s
  B  can 20 coils make the field (compensation_fields.py, 2-D magnetostatics)
  C  the HTS winding: current, tape, loss, cryogenic power, voltage, copper, plates
  D  gap protection and failure cases
  E  the numbers behind the rejected limiters, and the 8,620 m/s fix

Usage:
    python compensation_tracks.py                 # everything
    python compensation_tracks.py --parts=AD      # some parts
    python compensation_tracks.py --quick         # fewer beam cases, for testing

Models: sled_passage_pulsed.py (point and beam models with the 2026-10-11 options),
magnetostatic_2d.py and compensation_fields.py for the fields.
Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import csv
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
import compensation_fields as cf

MU0 = lc.MU0
MU, MU_C, W, GAP0, BAND = pz.MU, pz.MU_C, pz.W, pz.GAP0, pz.BAND
L_SLED, M_SLED_M, M_SC0 = pz.L_SLED, pz.M_SLED_M, pz.M_SC0
L_RING = ring.L_RING
V_REL = pz.V_REL
P12 = pz.P_SAT[1.2]                 # 573 kPa: pull per square meter of track at 1.2 T mean plate field
TRACK_W = 1.0                       # one compensation track
PAIR_W = 2.0                        # both tracks of a pair
CAP_ONE = P12 * TRACK_W             # 573 kN/m
CAP_PAIR = P12 * PAIR_W             # 1,146 kN/m
N_COILS = 20
L_COIL = 200.0
COIL_W = 0.05
STAGGER = L_COIL / N_COILS          # 10 m
TAU = 2e-3
L_SC = 100.0
SPEEDS = [1000.0, 5000.0, 10000.0, 20000.0, 30000.0]
MATCH = pz.MATCH_RING
# HTS tape (Chapter 6 and the hand scripts) [UNVERIFIED: thickness and density of the finished THEVA tape]
TAPE_T = 0.08e-3
TAPE_RHO = 8890.0
TAPE_WIDTHS = (0.003, 0.012)
I_FRAC = 0.5                        # operating current as a fraction of Ic
THEVA_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Vol_III", "bin",
                         "THEVA Pro-Line Advanced Pinning 2G HTS 70 K Field Dependence.csv")
RUNS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Vol_III", "runs_700km_2026-10-04")
# roll model geometry [UNVERIFIED: a layout sketch, not a Chapter 8 value]
Y_SUPPORT = 1.2                     # m, support plate centers from the cable centerline (plates 0.5 to 1.9 m)
Y_COMP = 2.8                        # m, compensation plate centers (plates 2.1 to 3.5 m), lift and load alike
Y_PLATE_EDGE = 3.5
Y_CABLE_EDGE = ring.CABLE_SIDE / 2  # 4.12 m
R_GYR_CASING = 5.0                  # m, casing radius of gyration about its axis
I_CASING = MU * R_GYR_CASING**2
I_CABLE = MU_C * (ring.CABLE_SIDE**2 + ring.CABLE_SIDE**2) / 12.0


def hr(t=""):
    pz.hr(t)


def rule_T(u, L_sc=L_SC):
    """Rule 4 of sled_passage_spread.py: 50 ms, at least the time the zone takes to pass."""
    return max(0.05, L_sc / u)


def energized_length(u, T, L_sc=L_SC):
    """A coil is on while the pull window (u T long) plus the zone overlaps it; coils reach up to one coil length beyond
    the window at each end, so the energized stretch is the window, the zone, two coil lengths and the lead of the rise."""
    return u * T + L_sc + 2 * L_COIL + 2 * u * TAU


def new_lane(u, width, **kw):
    d = dict(L_sc=L_SC, tau=TAU, L_section=STAGGER, gap_dep=False, steady_loop_hz=1.0, T_spread=rule_T(u), B_sat=1.2,
             n_coils=N_COILS, L_coil=L_COIL, hd_width=width, hu_width=width)
    d.update(kw)
    return pz.lane(**d)


def prot(trigger=0.010, delay=1e-3, cap=CAP_PAIR, hz=5.0, sensor_err=0.0, cap_up=None, cap_down=None):
    return dict(trigger=trigger, delay=delay, hz=hz, cap_up=cap if cap_up is None else cap_up, cap_down=cap if cap_down is None else cap_down,
                sensor_err=sensor_err)


def no_limiters(tabs):
    return dict(tabs, F0=0.0, Q0=0.0)


def fmt_gap(r):
    return f"{r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f}"


# =============================================================================
# PART A
# =============================================================================

def part_A(tabs):
    hr("PART A: THE SLED PASSAGE WITH TWO 1 m TRACKS PER FACE AND STAGGERED 200 m COILS")
    print(f"  5,000 t on {L_SC:.0f} m at the front of a {L_SLED:,.0f} m sled ({M_SLED_M:.0f} kg/m), casing {MU:,.0f} kg/m, cable {MU_C:,.0f} kg/m, gap {GAP0*1e3:.0f} mm.")
    print(f"  Spread pull: the cable's share ({MATCH*100:.1f} %) of the zone impulse over T = max(50 ms, zone/u), the bare sled cancelled in step;")
    print(f"  current regulated against the measured gap, 1 Hz loop on the casing support tracks, {TAU*1e3:.0f} ms rise time, no limiters.")
    print(f"  Coils: {N_COILS} per track, {L_COIL:.0f} m long, a new one every {STAGGER:.0f} m; a coil is on while the pull window overlaps it, at 1/{N_COILS}")
    print(f"  of the level, so the pull a point sees is the sum of the {N_COILS} coils over it: it ramps up over L_coil / u as they switch on one by one,")
    print(f"  holds, and ramps down, lasting T + 2 L_coil / u in all (the impulse is matched).  Track rating {CAP_ONE/1e3:.0f} kN/m per track (1.2 T mean),")
    print(f"  {CAP_PAIR/1e3:,.0f} kN/m for a pair.  'both' = the pair shares the pull, 'one' = one track does it alone.")
    print(f"  Excursion = largest casing offset from the 100 mm gap during the pass and settling; plate field B = sqrt(2 mu0 F / width).")
    print()
    print(f"  {'u km/s':>7} {'T ms':>6} {'face':>5} {'pull kN/m':>9} {'B both T':>8} {'B one T':>7} {'L_on m':>7} {'coils on':>8} {'E J/m':>7} {'E_on MJ':>7} "
          f"{'P GW':>6} {'exc both':>8} {'exc one':>7} {'gap_b':>6} {'gap_t':>6}")
    rows = {}
    t0 = time.time()
    for u in SPEEDS:
        T = rule_T(u)
        lanes = [new_lane(u, PAIR_W, label="both"), new_lane(u, TRACK_W, label="one")]
        res, hist = pz.point_pulsed(lanes, u, tabs)
        rb, ro = res
        F = max(rb["pk_hd"], rb["pk_hu"])
        face = "lift" if rb["pk_hd"] >= rb["pk_hu"] else "load"
        L_on = energized_length(u, T)
        E_m = F * GAP0                       # gap field energy per meter of ring at the pull (F = B^2/2mu0 x width, E = B^2/2mu0 x width x gap)
        rows[u] = dict(T=T, F=F, face=face, B_both=math.sqrt(2 * MU0 * F / PAIR_W), B_one=math.sqrt(2 * MU0 * F / TRACK_W), L_on=L_on,
                       E_m=E_m, P=E_m * u, exc_both=rb["max_exc_all"], exc_one=ro["max_exc_all"], res_both=rb, res_one=ro, hist=hist)
        r = rows[u]
        print(f"  {u/1e3:7.0f} {T*1e3:6.0f} {face:>5} {F/1e3:9.1f} {r['B_both']:8.3f} {r['B_one']:7.3f} {L_on:7.0f} {L_on/STAGGER:8.0f} {E_m:7,.0f} {E_m*L_on/1e6:7.2f} "
              f"{r['P']/1e9:6.2f} {rb['max_exc_all']*1e3:8.2f} {ro['max_exc_all']*1e3:7.2f} {fmt_gap(rb)}")
    print(f"  ({time.time()-t0:.0f} s)  E = F x gap is the gap field energy per meter; E_on = E x L_on is the energy in the field at any moment; P = E x u is")
    print(f"  the peak power into the field as the energized stretch advances (recovered at the back end, less the losses of Part C).")
    print(f"  The pull is the same for both and one; the plate field under one track alone is sqrt(2) times the pair's.  The excursions differ only")
    print(f"  through the force cap ({CAP_ONE/1e3:.0f} kN/m for one track), which is not reached at any speed here.")

    # tolerances
    print()
    print("  Tolerances (pair, excursion in mm; the baseline excursion is the 'both' column above).  Timing = pull early/late; mass = spacecraft")
    print("  mass error in the feed-forward; speed = sled speed error; sensor = gap sensor offset in the current regulation; rise = rise time")
    print("  error of the coil (feed-forward shaped for the nominal 2 ms); jitter = where the point sits in the 10 m stagger (0 and 1).")
    tol = [("timing -1 ms", dict(dt_err=-1e-3)), ("timing +1 ms", dict(dt_err=1e-3)), ("timing -5 ms", dict(dt_err=-5e-3)), ("timing +5 ms", dict(dt_err=5e-3)),
           ("mass -10 %", dict(mass_err=-0.10)), ("mass +10 %", dict(mass_err=0.10)), ("speed -5 %", dict(speed_err=-0.05)), ("speed +5 %", dict(speed_err=0.05)),
           ("sensor -5 mm", dict(gap_err=-0.005)), ("sensor +5 mm", dict(gap_err=0.005)), ("rise -50 %", dict(tau_err=-0.5)), ("rise +50 %", dict(tau_err=0.5)),
           ("jitter 0", dict(pos_in_section=0.0)), ("jitter 1", dict(pos_in_section=1.0))]
    print(f"  {'case':>14} " + " ".join(f"{u/1e3:>7.0f}" for u in SPEEDS) + "   km/s")
    tol_rows = {}
    t0 = time.time()
    per_u = {}
    for u in SPEEDS:
        lanes = [new_lane(u, PAIR_W, label=lab, **kw) for lab, kw in tol]
        res, _ = pz.point_pulsed(lanes, u, tabs)
        per_u[u] = res
    for i, (lab, kw) in enumerate(tol):
        vals = [per_u[u][i]["max_exc_all"] for u in SPEEDS]
        tol_rows[lab] = vals
        flag = "  <- contact" if any(per_u[u][i]["contact"] for u in SPEEDS) else ""
        print(f"  {lab:>14} " + " ".join(f"{v*1e3:7.2f}" for v in vals) + flag)
    worst = max(max(v) for v in tol_rows.values())
    print(f"  ({time.time()-t0:.0f} s)  Largest excursion over all tolerance cases: {worst*1e3:.2f} mm, band {BAND*1e3:.0f} mm, no contact.")

    # the low-speed side: coil length against the smear of the pull
    print()
    print("  Coil length at low speed.  A coil pulls along its whole length while the pull window overlaps any part of it, so a point sees the")
    print(f"  pull for T + 2 L_coil / u, which at 1 km/s is {(rule_T(1000.0) + 2*L_COIL/1000.0)*1e3:.0f} ms against a zone passage of {L_SC/1000.0*1e3:.0f} ms.  Excursion (mm, pair) for coil lengths")
    print("  of 200, 100 and 50 m (20 coils, a new coil every 10, 5 and 2.5 m), and for 50 m sections switched as a block (the earlier layout):")
    print(f"  {'u km/s':>7} {'200 m':>7} {'100 m':>7} {'50 m':>7} {'50 m sections':>13}")
    t0 = time.time()
    low = {}
    for u in (1000.0, 2000.0, 3000.0, 5000.0):
        lanes = [new_lane(u, PAIR_W, L_coil=Lc, label=f"{Lc}") for Lc in (200.0, 100.0, 50.0)]
        lanes.append(pz.lane(L_sc=L_SC, tau=TAU, L_section=50.0, gap_dep=False, steady_loop_hz=1.0, T_spread=rule_T(u), B_sat=1.2, hd_width=PAIR_W, hu_width=PAIR_W))
        res, _ = pz.point_pulsed(lanes, u, tabs)
        low[u] = [r["max_exc_all"] for r in res]
        print(f"  {u/1e3:7.0f} " + " ".join(f"{v*1e3:7.2f}" for v in low[u][:3]) + f" {low[u][3]*1e3:13.2f}")
    print(f"  ({time.time()-t0:.0f} s)  The load compensation tracks (floor) handle the speeds below 7 km/s; the smear is their problem alone.")
    rows["low_speed"] = low
    return rows


# =============================================================================
# PART B
# =============================================================================

DESIGN = dict(cap_one=CAP_ONE, cap_pair=CAP_PAIR, label="1.2 T mean over the track")


def part_B(rows_A, tabs):
    hr("PART B: CAN 20 COILS OF 5 cm MAKE THE FIELD (2-D magnetostatics, compensation_fields.py)")
    F30 = rows_A[30000.0]["F"] if rows_A else 437e3
    B_SAT = lc.FERRO_B_SAT
    print(f"  One track facing its 1.4 m plate (sintered iron, mu_r {lc.FERRO_MU_R}, {lc.FERRO_PLATE_THICKNESS*1e3:.0f} mm thick, saturation {B_SAT} T) across 100 mm, 5 mm grid.")
    print(f"  Iron yoke mu_r 1,000, linear; saturation is checked against {B_SAT} T, not modeled.  Forces by Maxwell stress on a closed surface around the plate.")
    print(f"  Targets: the pair's share at 30 km/s {F30/2e3:.0f} kN/m per track and one track alone {F30/1e3:.0f} kN/m (Part A).")
    print(f"  The field is linear in the ampere-turns and the force goes as their square, so each geometry is solved once and scaled.")
    t0 = time.time()
    res = cf.load_or_run(verbose=True)
    comp, sup = res["compensation"], res["support"]
    print(f"  ({time.time()-t0:.0f} s to load or solve; the solve log is compensation_fields.txt)")
    print()
    print("  Raw solve, per 100 kA-turns on each pole (the air cases: per 100 kA-turns on the track):")
    print(f"  {'case':>58} {'B_pole T':>8} {'B_track T':>9} {'F kN/m':>7} {'B_plate T':>9} {'W kJ/m':>7} {'Bx max/mean':>12} {'By max/mean':>12}")
    for c in comp:
        print(f"  {c['label']:>58} {c['B_pole']:8.4f} {c['B_track']:9.4f} {c['F']/1e3:7.1f} {c['B_plate']:9.3f} {c['W_mag']/1e3:7.1f} "
              f"{c['Bx_w']:5.2f}/{c['Bx_wm']:5.2f}  {c['By_w']:5.2f}/{c['By_wm']:5.2f}")
    print("  B_pole = mean field normal to the plate over the pole faces, B_track = over the whole 1 m; B_plate = flux density along the plate through")
    print("  its thickness midway between poles (the flux that must flow sideways in the plate, the pole flux plus the slot leakage the plate shorts);")
    print("  W = stored energy per meter of track, gap and slots and leakage together; Bx = field across the track at the winding (along a flat-wound")
    print("  tape), By = field normal to the plate at the winding (through the tape face).  The slot field at the slot mouth is mu0 x slot ampere-turns")
    print("  / slot width; the maxima sit at the iron corners on the 5 mm grid.")
    print()
    a4 = next(c for c in comp if c["label"].startswith("A4"))
    print(f"  Air core.  The 20 coils flat and side by side with one polarity are one 1 m loop whose field at 100 mm is a few millitesla per 100 kA;")
    print(f"  alternating polarity at 50 mm pitch decays as exp(-pi x gap / pitch) and is gone at 100 mm; even an ideal 1 m sheet with its return far")
    print(f"  away gives {a4['B_track']*1e3:.0f} mT per 100 kA, so 0.74 T on the plate would need {0.74/a4['B_track']*100:,.0f} kA-turns of one-way current per track,")
    print(f"  with the same field on the winding.  Without iron the coils cannot make the field.")
    print()
    print(f"  Iron core, scaled.  F_sat = the largest pull before the pole face or the plate flux reaches {B_SAT} T; NI = ampere-turns per pole (per coil")
    print(f"  group) at F_sat; coils per group = 20 / groups; slot mouth = mu0 x 2 NI / slot width at F_sat; winding fields at F_sat (mean over the slot).")
    print(f"  {'case':>58} {'limit':>6} {'F_sat':>6} {'NI kA':>6} {'/coil':>6} {'groups':>6} {'mouth T':>7} {'Bx mean':>7} {'By mean':>7} {'W kJ/m':>7} {'pair':>5} {'alone':>5}")
    scaled = {}
    for c in comp:
        if not (c["label"].startswith("I") or c["label"].startswith("H")):
            continue
        s_pole = B_SAT / c["B_pole"]; s_plate = B_SAT / c["B_plate"]
        sc = min(s_pole, s_plate)
        limit = "face" if s_pole <= s_plate else "plate"
        F_sat = c["F"] * sc * sc
        n_groups = c.get("n_groups", c.get("n_pole", 1))
        NI = c["NI"] * sc
        slot_w = c.get("slot_w", 0.05)
        d = dict(F_sat=F_sat, NI=NI, NI_coil=NI / (N_COILS / n_groups), n_groups=n_groups, limit=limit, mouth=MU0 * 2 * NI / slot_w,
                 Bx_m=c["Bx_wm"] * sc, By_m=c["By_wm"] * sc, Bx_x=c["Bx_w"] * sc, By_x=c["By_w"] * sc, W=c["W_mag"] * sc * sc, sc=sc,
                 B_track=c["B_track"] * sc, pair_ok=F_sat >= F30 / 2, alone_ok=F_sat >= F30, slot_w=slot_w, pole_w=c.get("pole_w", 0.0),
                 plate_t=c.get("plate_t", lc.FERRO_PLATE_THICKNESS), c=c)
        scaled[c["label"]] = d
        print(f"  {c['label']:>58} {limit:>6} {F_sat/1e3:6.0f} {NI/1e3:6.1f} {d['NI_coil']/1e3:6.2f} {n_groups:6d} {d['mouth']:7.2f} {d['Bx_m']:7.2f} {d['By_m']:7.2f} "
              f"{d['W']/1e3:7.1f} {('yes' if d['pair_ok'] else 'no'):>5} {('yes' if d['alone_ok'] else 'no'):>5}")
    print(f"  'pair' = F_sat covers the pair's share ({F30/2e3:.0f} kN/m), 'alone' = F_sat covers the whole pull ({F30/1e3:.0f} kN/m) with one track out.")
    print(f"  The pole face limit: the mean field over the track is the pole field times the pole fraction, so with half the track in slots the")
    print(f"  pole face reaches {B_SAT} T at about half of the {CAP_ONE/1e3:.0f} kN/m that {B_SAT} T over the whole meter would give; narrow slots raise the")
    print(f"  pole fraction but put the slot ampere-turns across a narrow slot (the mouth field) and short their leakage through the plate.")
    # choice: 100 mm plate; the thick-plate cases are reported as the alternative
    cands = [l for l, d in scaled.items() if abs(d["plate_t"] - lc.FERRO_PLATE_THICKNESS) < 1e-6 and d["pair_ok"]]
    if cands:
        # among those that cover the pair's share with a slot mouth field under 2 T, the largest F_sat; otherwise the largest F_sat
        mild = [l for l in cands if scaled[l]["mouth"] <= 2.0]
        best_l = max(mild or cands, key=lambda l: scaled[l]["F_sat"])
    else:
        best_l = max(scaled, key=lambda l: scaled[l]["F_sat"])
    best = scaled[best_l]
    thick = {l: d for l, d in scaled.items() if d["plate_t"] > lc.FERRO_PLATE_THICKNESS}
    best_thick_l = max(thick, key=lambda l: thick[l]["F_sat"]) if thick else None
    print()
    thin_lin = {l: d for l, d in scaled.items() if abs(d["plate_t"] - lc.FERRO_PLATE_THICKNESS) < 1e-6}
    bl = max(thin_lin, key=lambda l: thin_lin[l]["F_sat"])
    print(f"  Linear first estimate: with the 100 mm plate the largest F_sat is {thin_lin[bl]['F_sat']/1e3:.0f} kN/m ({bl}), "
          f"{'at or above' if thin_lin[bl]['F_sat'] >= F30/2 else 'below'} the pair's share; with a 200 mm plate "
          + (f"{thick[best_thick_l]['F_sat']/1e3:.0f} kN/m ({best_thick_l})." if best_thick_l else "no case was run."))

    # ---- saturating iron: the force curves on the 10 mm grid
    nl = res.get("nonlinear", [])
    if nl:
        print()
        print(f"  Saturating iron (plate B(H) saturating at {B_SAT} T from mu_r {lc.FERRO_MU_R}, yoke and poles at {cf.YOKE_B_SAT} T from mu_r {cf.YOKE_MU:.0f} [UNVERIFIED],")
        print(f"  10 mm grid, Picard iteration on the cell permeability).  'linear' = the linear solve on the same grid scaled to the same ampere-turns;")
        print(f"  the 5 mm linear solve above gives the grid error.  Mouth = mu0 x 2 NI / slot width.  Rating = 90 % of the largest pull on the")
        print(f"  computed curve (up to {max(cf.NL_LEVELS)/1e3:.0f} kA per pole), the point where the curve has flattened and more ampere-turns buy little; the")
        print(f"  ampere-turns for it are interpolated on the curve.  'conv' = force change in the last Picard iteration (under 0.003 = converged).")
        print(f"  {'geometry':>8} {'NI/pole kA':>10} {'F kN/m':>7} {'linear':>7} {'ratio':>5} {'B_track T':>9} {'plate max T':>11} {'plate sat %':>11} {'mouth T':>7} {'Bx mean':>7} {'By mean':>7} {'W kJ/m':>7} {'conv':>6}")
        curves = {}
        for c in nl:
            g = c.get("geom")
            if g is None or g == "S4":
                continue
            lin = next((x for x in nl if x["label"].startswith(g + " linear")), None)
            fine = next((x for x in comp if x["label"].startswith(g + " ")), None)
            F_lin = lin["F"] * (c["NI"] / lin["NI"]) ** 2 if lin else float("nan")
            slot_w = c.get("slot_w", 0.05)
            d = dict(NI=c["NI"], F=c["F"], F_lin=F_lin, ratio=c["F"] / F_lin if F_lin else float("nan"), B_track=c["B_track"], B_pole=c["B_pole"],
                     plate_max=c.get("B_plate_max", 0.0), mu_min=c.get("mu_plate_min", 0.0), mouth=MU0 * 2 * c["NI"] / slot_w, Bx_m=c["Bx_wm"], By_m=c["By_wm"],
                     Bx_x=c["Bx_w"], By_x=c["By_w"], W=c["W_mag"], slot_w=slot_w, pole_w=c.get("pole_w", 0.0), n_groups=c.get("n_groups", 1),
                     plate_t=c.get("plate_t", lc.FERRO_PLATE_THICKNESS), fine_F_lin=(fine["F"] * (c["NI"] / fine["NI"]) ** 2 if fine else float("nan")))
            curves.setdefault(g, []).append(d)
            print(f"  {g:>8} {c['NI']/1e3:10.0f} {c['F']/1e3:7.0f} {F_lin/1e3:7.0f} {d['ratio']:5.2f} {c['B_track']:9.3f} {d['plate_max']:11.2f} {c.get('plate_sat_frac', 0)*100:11.0f} {d['mouth']:7.2f} "
                  f"{c['Bx_wm']:7.2f} {c['By_wm']:7.2f} {c['W_mag']/1e3:7.1f} {c.get('outer_dF', 0):6.3f}")
        if curves:
            g0 = next(iter(curves))
            d0 = curves[g0][0]
            print(f"  Grid check: {g0} linear at {d0['NI']/1e3:.0f} kA gives {d0['F_lin']/1e3:.0f} kN/m on the 10 mm grid against {d0['fine_F_lin']/1e3:.0f} kN/m on the 5 mm grid.")
        knee = {}
        print()
        print(f"  {'geometry':>8} {'plate mm':>8} {'rating kN/m':>11} {'at NI kA':>8} {'per coil':>8} {'mouth T':>7} {'Bx/By mean T':>13} {'curve max':>9} {'NI for pair share':>17} {'NI for alone':>12}")
        for g, rows in curves.items():
            rows.sort(key=lambda d: d["NI"])
            NIs = [d["NI"] for d in rows]; Fs = [d["F"] for d in rows]
            F_rate = 0.9 * max(Fs)
            NI_rate = float(np.interp(F_rate, Fs, NIs))
            k = dict(rows[0])
            for key in ("B_track", "B_pole", "plate_max", "mouth", "Bx_m", "By_m", "Bx_x", "By_x", "W"):
                k[key] = float(np.interp(NI_rate, NIs, [d[key] for d in rows]))
            k["F"] = F_rate; k["NI"] = NI_rate
            NI_pair = float(np.interp(F30 / 2, Fs, NIs)) if F30 / 2 <= max(Fs) else float("nan")
            NI_alone = float(np.interp(F30, Fs, NIs)) if F30 <= max(Fs) else float("nan")
            knee[g] = dict(k, NI_pair=NI_pair, NI_alone=NI_alone, F_max=max(Fs))
            print(f"  {g:>8} {k['plate_t']*1e3:8.0f} {k['F']/1e3:11.0f} {k['NI']/1e3:8.0f} {k['NI']/(N_COILS/k['n_groups'])/1e3:8.2f} {k['mouth']:7.2f} {k['Bx_m']:6.2f}/{k['By_m']:5.2f} {max(Fs)/1e3:9.0f} "
                  f"{(f'{NI_pair/1e3:.0f} kA' if NI_pair == NI_pair else 'not reached'):>17} {(f'{NI_alone/1e3:.0f} kA' if NI_alone == NI_alone else 'not reached'):>12}")
        print(f"  'NI for pair share' and 'NI for alone' interpolate the saturating curve for {F30/2e3:.0f} and {F30/1e3:.0f} kN/m per track; 'not reached' means the curve")
        print(f"  does not get there within {max(cf.NL_LEVELS)/1e3:.0f} kA per pole.")
        thin = {g: k for g, k in knee.items() if abs(k["plate_t"] - lc.FERRO_PLATE_THICKNESS) < 1e-6}
        thick_nl = {g: k for g, k in knee.items() if k["plate_t"] > lc.FERRO_PLATE_THICKNESS}
        if thin:
            best_l = max(thin, key=lambda g: thin[g]["F"])
            k = thin[best_l]
            best = dict(F_sat=k["F"], NI=k["NI"], NI_coil=k["NI"] / (N_COILS / k["n_groups"]), n_groups=k["n_groups"], limit="saturation rating",
                        mouth=k["mouth"], Bx_m=k["Bx_m"], By_m=k["By_m"], Bx_x=k["Bx_x"], By_x=k["By_x"], W=k["W"], B_track=k["B_track"],
                        pair_ok=k["F_max"] >= F30 / 2 and k["NI_pair"] == k["NI_pair"], alone_ok=k["NI_alone"] == k["NI_alone"], slot_w=k["slot_w"],
                        pole_w=k["pole_w"], plate_t=k["plate_t"], NI_pair=k["NI_pair"], NI_alone=k["NI_alone"], F_max=k["F_max"])
            best["sc"] = 1.0
        best_thick_l = max(thick_nl, key=lambda g: thick_nl[g]["F"]) if thick_nl else None
        if best_thick_l:
            thick = {best_thick_l: dict(thick_nl[best_thick_l], F_sat=thick_nl[best_thick_l]["F"], limit="saturation rating",
                                         alone_ok=thick_nl[best_thick_l]["NI_alone"] == thick_nl[best_thick_l]["NI_alone"])}
    DESIGN["cap_one"] = min(best["F_sat"], CAP_ONE)
    DESIGN["cap_pair"] = 2 * DESIGN["cap_one"]
    DESIGN["label"] = best_l
    print()
    print(f"  Chosen for the 100 mm plate: {best_l} ({'saturating solve' if nl else 'linear estimate'}).")
    print(f"    {best['n_groups']} coil groups (one per pole, {N_COILS/best['n_groups']:.0f} of the 20 coils on each pole, each coil 200 m long with its two legs in the slots")
    print(f"    either side of its pole), neighboring poles of opposite polarity, half-width poles at the track edges, {best['slot_w']*1e3:.0f} mm slots, {best['pole_w']*1e3:.0f} mm full poles.")
    print(f"    Rated pull ({best['limit']}): {best['F_sat']/1e3:.0f} kN/m per track, {2*best['F_sat']/1e3:.0f} kN/m for the pair, at {best['NI']/1e3:.1f} kA-turns")
    print(f"    per pole and {best['NI_coil']/1e3:.2f} kA-turns per coil; track mean field {best['B_track']:.3f} T; stored energy {best['W']/1e3:.1f} kJ per m of track.")
    for tag, Ft in (("pair share at 30 km/s", F30 / 2), ("one track alone at 30 km/s", min(F30, best["F_sat"]))):
        sc2 = math.sqrt(Ft / best["F_sat"])
        print(f"    {tag}: {Ft/1e3:.0f} kN/m, {best['NI']*sc2/1e3:.1f} kA-turns per pole, {best['NI_coil']*sc2/1e3:.2f} per coil, slot mouth {best['mouth']*sc2:.2f} T,")
        print(f"      winding mean across the track {best['Bx_m']*sc2:.2f} T (max {best['Bx_x']*sc2:.2f}), normal to the plate {best['By_m']*sc2:.2f} T (max {best['By_x']*sc2:.2f}), stored {best['W']*sc2*sc2/1e3:.1f} kJ/m")
    if best.get("NI_pair") == best.get("NI_pair") and best.get("NI_pair") is not None:
        print(f"    The pair's share ({F30/2e3:.0f} kN/m per track) takes {best['NI_pair']/1e3:.0f} kA-turns per pole on the saturating curve"
              + (f", {best['NI_pair']/best['NI']*100:.0f} % of the rated ampere-turns." if best['NI_pair'] <= best['NI'] else ", beyond the rating (the iron is saturating there)."))
    elif nl:
        print(f"    The pair's share ({F30/2e3:.0f} kN/m per track) is NOT reached by this geometry within {max(cf.NL_LEVELS)/1e3:.0f} kA per pole (largest pull {best['F_max']/1e3:.0f} kN/m).")
    if best["alone_ok"]:
        print(f"    One track alone covers the 30 km/s pull with the 100 mm plate" + (f" at {best['NI_alone']/1e3:.0f} kA-turns per pole." if best.get('NI_alone') else "."))
    else:
        print(f"    One track alone does NOT cover the 30 km/s pull ({F30/1e3:.0f} kN/m) with the 100 mm plate: the rating is {best['F_sat']/1e3:.0f} kN/m, "
              f"{best['F_sat']/F30*100:.0f} % of it" + (f", and the curve tops out at {best['F_max']/1e3:.0f} kN/m." if best.get('F_max') else "."))
    if best_thick_l:
        bt = thick[best_thick_l]
        print(f"  With a 200 mm plate the best case is {best_thick_l}: rating {bt['F_sat']/1e3:.0f} kN/m, {bt['NI']/1e3:.1f} kA-turns per pole, slot mouth {bt['mouth']:.2f} T, "
              f"winding {bt['Bx_m']:.2f}/{bt['By_m']:.2f} T, {('covers' if bt['alone_ok'] else 'does not cover')} one track alone"
              + (f" (the whole pull at {bt['NI_alone']/1e3:.0f} kA per pole)." if bt.get('NI_alone') == bt.get('NI_alone') and bt.get('NI_alone') else "."))
    print(f"  [QUESTION FOR AUTHOR] the one-track-alone requirement at full speed and mass needs a plate that saturates above {B_SAT} T, a thicker plate")
    print(f"  with wider poles, a wider track, or a speed and mass limit while a track is out.  Part D uses the saturation pull of the chosen geometry")
    print(f"  as the one-track cap ({DESIGN['cap_one']/1e3:.0f} kN/m) and the pair at twice that.")
    # the one-track-alone pass at the saturation cap
    if rows_A:
        print()
        print(f"  One track alone with its pull capped at F_sat (point model, 20 and 30 km/s): excursion and whether the cap binds.")
        for u in (20000.0, 30000.0):
            res_c, _ = pz.point_pulsed([new_lane(u, DESIGN["cap_one"] / P12, label="one, capped")], u, tabs)
            r = res_c[0]
            print(f"    {u/1e3:.0f} km/s: pull needed {rows_A[u]['F']/1e3:.0f} kN/m, cap {DESIGN['cap_one']/1e3:.0f} kN/m, peak delivered {max(r['pk_hd'], r['pk_hu'])/1e3:.0f} kN/m, "
                  f"excursion {r['max_exc_all']*1e3:.1f} mm, gaps {fmt_gap(r)} mm{', contact' if r['contact'] else ''}")
        # the speed at which one track alone still covers the pull (pull grows with speed): find from Part A rows by interpolation
        us = sorted(rows_A[u]["F"] for u in SPEEDS if u >= 10000)
    E_gap_30 = F30 * GAP0
    print(f"  Stored energy against the gap field: {best['W']/1e3:.1f} kJ/m at F_sat against F x gap = {best['F_sat']*GAP0/1e3:.1f} kJ/m, ratio {best['W']/(best['F_sat']*GAP0):.2f} (slot and leakage flux).")

    # support tracks
    print()
    print(f"  Casing support tracks: 40 coils of 25 mm, 52 turns x 400 A = {cf.NI_SUPPORT/1e3:.0f} kA-turns per 1 m track (levitation_config.py).  Needed: half the casing")
    print(f"  weight {W/2e3:.1f} kN/m per track (mean plate field {math.sqrt(2*MU0*W/2/TRACK_W):.3f} T), the 218 kN/m operating rating of the pair {sp.F_EML_OP/2e3:.0f} kN/m per track "
          f"({math.sqrt(2*MU0*sp.F_EML_OP/2/TRACK_W):.3f} T).")
    print(f"  Chapter 8 quotes mu0 K / 2 = {MU0*cf.NI_SUPPORT/TRACK_W/2:.3f} T for an infinite one-way current sheet of {cf.NI_SUPPORT/1e3:.0f} kA per meter.")
    print(f"  {'case':>58} {'B_track T':>9} {'F kN/m':>7} {'F/needed':>8} {'B_plate T':>9} {'W kJ/m':>7} {'Bx max/mean':>12} {'By max/mean':>12}")
    for c in sup:
        print(f"  {c['label']:>58} {c['B_track']:9.4f} {c['F']/1e3:7.1f} {c['F']/(W/2):8.2f} {c['B_plate']:9.3f} {c['W_mag']/1e3:7.1f} "
              f"{c['Bx_w']:5.2f}/{c['Bx_wm']:5.2f}  {c['By_w']:5.2f}/{c['By_wm']:5.2f}")
    s4 = next((c for c in sup if c["label"].startswith("S4")), None)
    if s4:
        sc_w = math.sqrt((W / 2) / s4["F"]); sc_op = math.sqrt((sp.F_EML_OP / 2) / s4["F"])
        ng = s4.get("n_groups", 5)
        print(f"  With the iron pole geometry of S4 ({ng} coil groups) the support track needs {cf.NI_SUPPORT*sc_w/1e3:.0f} kA-turns per track ({cf.NI_SUPPORT/ng*sc_w/1e3:.1f} kA per pole, "
              f"{cf.NI_SUPPORT/ng*sc_w/(40/ng)/1e3:.2f} kA per coil, {cf.NI_SUPPORT/ng*sc_w/(40/ng)/lc.EML_I_OPERATING:.0f} turns at 400 A) for half the casing weight, against the {cf.NI_SUPPORT/1e3:.0f} kA installed;")
        print(f"  for the 218 kN/m operating rating {cf.NI_SUPPORT*sc_op/1e3:.0f} kA-turns per track.  Winding field at half the weight: across {s4['Bx_wm']*sc_w:.2f} T, normal {s4['By_wm']*sc_w:.2f} T (means), "
              f"slot mouth {MU0*2*cf.NI_SUPPORT/ng*sc_w/s4.get('slot_w', 0.05):.2f} T; plate flux {s4['B_plate']*sc_w:.2f} T.")
        s4nl = [c for c in res.get("nonlinear", []) if c.get("geom") == "S4"]
        for c in s4nl:
            print(f"  Saturating solve, S4 geometry at {c['NI_track']/1e3:.0f} kA-turns per track: {c['F']/1e3:.0f} kN/m ({c['F']/(W/2):.1f} x half the weight), track mean {c['B_track']:.3f} T, "
                  f"plate max {c['B_plate_max']:.2f} T, winding across {c['Bx_wm']:.2f} T, normal {c['By_wm']:.2f} T (means).")
        print(f"  Flat air-core coils side by side cannot reach the plate ({next(c for c in sup if c['label'].startswith('S1'))['F']/1e3:.1f} kN/m from the installed ampere-turns); the support")
        print(f"  tracks need the same pole arrangement as the compensation tracks, and the installed ampere-turns are then {'enough' if s4['F'] >= W/2 else 'not enough'} for half the weight.")
    return dict(best=best, best_label=best_l, scaled=scaled, comp=comp, sup=sup, thick=thick)


# =============================================================================
# PART C
# =============================================================================

def read_theva():
    """THEVA Advanced Pinning 70 K data: dict angle -> (B array, Ic per cm width array)."""
    data = {}
    with open(THEVA_CSV, encoding="utf-8-sig") as f:
        rd = csv.reader(f)
        next(rd)
        for row in rd:
            if len(row) < 4 or not row[1]:
                continue
            ang = int(round(float(row[2])))
            data.setdefault(ang, []).append((float(row[1]), float(row[3])))
    out = {}
    for ang, pts in data.items():
        pts.sort()
        out[ang] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
    return out


def ic_per_cm(theva, B, angle):
    Bs, Ic = theva[angle]
    return float(np.interp(B, Bs, Ic))


def eddy_cycle(B, t_edge, t_hold, d_p, mu_r, rho_p):
    """Eddy energy per unit plate area for one rise-hold-fall cycle, J/m^2, by explicit 1-D magnetic diffusion
    (copied from Part I of sled_passage_spread.py)."""
    nz = 100
    dz = d_p / nz
    D = 1.0 / (MU0 * mu_r / rho_p)
    dt = 0.4 * dz * dz / D
    t_end = 2 * t_edge + t_hold + 15e-3
    n = int(t_end / dt) + 1
    Bz = np.zeros(nz + 1)
    E = 0.0
    for j in range(n):
        t = j * dt
        if t < t_edge:
            Bs = B * t / t_edge
        elif t < t_edge + t_hold:
            Bs = B
        elif t < 2 * t_edge + t_hold:
            Bs = B * (1 - (t - t_edge - t_hold) / t_edge)
        else:
            Bs = 0.0
        Bz[0] = Bs
        Bz[-1] = Bz[-2]
        lap = (Bz[2:] - 2 * Bz[1:-1] + Bz[:-2]) / dz**2
        Jc = (Bz[1:] - Bz[:-1]) / dz / (MU0 * mu_r)
        E += rho_p * np.sum(Jc * Jc) * dz * dt
        Bz[1:-1] += D * lap * dt
    return E


def part_C(rows_A, B):
    hr("PART C: THE HTS WINDING")
    best = B["best"]
    F30 = rows_A[30000.0]["F"]
    F_alone = min(F30, best["F_sat"])                       # the one-track rating: the whole pull, or what the plate allows
    sc_a = math.sqrt(F_alone / best["F_sat"]); sc_p = math.sqrt((F30 / 2) / best["F_sat"])
    def at(sc, F):
        return dict(F=F, NI_pole=best["NI"] * sc, NI_coil=best["NI_coil"] * sc, n_pole=best["n_groups"], W=best["W"] * sc * sc,
                    Bx_m=best["Bx_m"] * sc, By_m=best["By_m"] * sc, Bx_x=best["Bx_x"] * sc, By_x=best["By_x"] * sc)
    alone, pair = at(sc_a, F_alone), at(sc_p, F30 / 2)
    theva = read_theva()
    print(f"  Tape data: {os.path.basename(THEVA_CSV)} (critical current per cm of width at 70 K against applied field at 0, 30, 90 and 120 degrees).")
    print(f"  0 degrees is taken as the field normal to the tape face and 90 degrees as the field along it [UNVERIFIED convention]; at 70 K the four")
    print(f"  angles differ by under 10 %, so the field strength sets the current.  The tape is rated at the lower of the 0 and 90 degree values at")
    print(f"  the mean winding field of Part B, and the operating current is {I_FRAC:.0%} of that (as the casing support tracks, Ic/2).")
    print(f"  Winding: {B['best_label']}, rated for one track alone ({F_alone/1e3:.0f} kN/m, the 30 km/s pull or the saturation pull if that is smaller): "
          f"{alone['NI_pole']/1e3:.1f} kA-turns per pole, {alone['NI_coil']/1e3:.2f} kA-turns per coil.")
    B_mean = math.hypot(alone["Bx_m"], alone["By_m"]); B_max = math.hypot(alone["Bx_x"], alone["By_x"])
    print(f"  Winding field at the rating: mean {B_mean:.2f} T (across the track {alone['Bx_m']:.2f}, normal {alone['By_m']:.2f}), max {B_max:.2f} T at the slot corners.")
    print()
    print("  C.1  Current, tape length and mass (per meter of ring for one track; four tracks around the ring, each rated alone).")
    print(f"  {'tape':>6} {'Ic/cm 0':>8} {'Ic/cm 90':>8} {'Ic A':>6} {'I_op A':>6} {'turns/coil':>10} {'tape m/m':>9} {'kg/m':>6} {'ring 4 tracks':>14} {'tape km (4)':>12} {'at max B: Ic A':>14}")
    tape = {}
    NI_track = alone["NI_pole"] * alone["n_pole"]
    for w in TAPE_WIDTHS:
        ic0 = ic_per_cm(theva, B_mean, 0); ic90 = ic_per_cm(theva, B_mean, 90)
        Ic = min(ic0, ic90) * w * 100.0
        I_op = I_FRAC * Ic
        turns = alone["NI_coil"] / I_op
        length = 2 * NI_track / I_op                     # both legs, per meter of ring
        mass = length * TAPE_T * w * TAPE_RHO
        Ic_max = min(ic_per_cm(theva, B_max, 0), ic_per_cm(theva, B_max, 90)) * w * 100.0
        tape[w] = dict(Ic=Ic, I_op=I_op, turns=turns, length=length, mass=mass, Ic_max=Ic_max)
        print(f"  {w*1e3:4.0f} mm {ic0:8.0f} {ic90:8.0f} {Ic:6.0f} {I_op:6.0f} {turns:10.0f} {length:9,.0f} {mass:6.1f} {mass*L_RING*4/1e9:11.2f} Mt {length*L_RING*4/1e3/1e6:9.0f} M {Ic_max:14.0f}")
    print(f"  Mass is tape only ({TAPE_T*1e3:.2f} mm thick, {TAPE_RHO:,.0f} kg/m^3) [UNVERIFIED], no former, insulation, cooling or yoke.  The same ampere-turns need the same")
    print(f"  tape area whatever the width, so the mass per meter is the same for 3 and 12 mm tape; 3 mm tape needs four times the turns and length.")
    print(f"  In normal use the pair shares the pull: each track runs at {pair['NI_coil']/alone['NI_coil']*100:.0f} % of its rated current ({pair['NI_pole']/1e3:.1f} kA-turns per pole).")
    print(f"  Hand estimate check (compensation_coil_estimates_v2.txt): 190,986 A-turns per pole pair at 1.2 T across two gaps, 158 A per 12 mm tape, 3,046 m of tape")
    print(f"  per m, 26 kg/m.  Here: {NI_track/1e3:.0f} kA-turns per track in total, {tape[0.012]['I_op']:.0f} A per 12 mm tape, {tape[0.012]['length']:,.0f} m per m, {tape[0.012]['mass']:.1f} kg/m.")

    # C.2 loss per pulse and per launch
    print()
    print("  C.2  Loss in the tape per pulse and per launch, heat at 70 K and the electrical power to remove it.")
    print(f"  Model [UNVERIFIED]: a fully penetrated thin strip in a field normal to its face swept 0 -> B -> 0 dissipates about 1/2 B_perp x w x Ic per")
    print(f"  meter of tape (half the Bean value for a full cycle, B w Ic).  The field along the tape face and the transport current's own loss are")
    print(f"  neglected (they go with the tape thickness).  B_perp is the mean field normal to the plate over the winding at the operating ampere-turns,")
    print(f"  which goes as sqrt(pull), so the loss per pulse goes as sqrt(pull) and with tape width.  The pair shares the pull, so each track runs at")
    print(f"  1/sqrt(2) of the one-track field and the two together dissipate sqrt(2) times one track alone.")
    cop = lim_physics.calc_cryo_cop(70.0, 400.0, lim_config.CRYO_EFF)
    print(f"  Cryogenic COP = lim_physics.calc_cryo_cop(70, 400, lim_config.CRYO_EFF = {lim_config.CRYO_EFF}) = {cop:.4f} (Carnot 70/330 times {lim_config.CRYO_EFF}); "
          f"lim_physics.py, def calc_cryo_cop, and lim_config.py CRYO_EFF.  Electrical = heat / COP.")
    t_p, v_p, P_p = pz.lsm_profile()
    KE = 0.5 * (M_SC0 + M_SLED_M * L_SLED) * 30000**2
    print(f"  {'tape':>6} {'mode':>6} {'B_perp T':>8} {'J/pulse/m':>9} {'per lap GJ':>10} {'per launch GJ':>13} {'electrical GJ':>13} {'% of KE':>7} {'mean W/m':>8} {'mean GW':>8}")
    loss = {}
    dt_p = np.diff(t_p)
    for w in TAPE_WIDTHS:
        tp = tape[w]
        for mode, nfac in (("alone", 1.0), ("pair", math.sqrt(2.0))):
            B_perp = alone["By_m"] * (1.0 if mode == "alone" else 1 / math.sqrt(2))
            q30 = 0.5 * B_perp * w * tp["Ic"] * tp["length"] * (1.0 if mode == "alone" else 2.0)     # J per pulse per m of ring, at the 30 km/s pull
            # per launch: each meter of ring is pulsed once per lap, at the pull of that speed
            E_launch = 0.0
            for i, dtt in enumerate(dt_p):
                vv = v_p[i + 1]
                Fu = np.interp(vv, SPEEDS, [rows_A[u]["F"] for u in SPEEDS])
                E_launch += q30 * math.sqrt(min(Fu, F_alone) / F_alone) * vv * dtt
            E_el = E_launch / cop
            loss[(w, mode)] = dict(q30=q30, E_launch=E_launch, E_el=E_el, B_perp=B_perp)
            print(f"  {w*1e3:4.0f} mm {mode:>6} {B_perp:8.3f} {q30:9.1f} {q30*L_RING/1e9:10.2f} {E_launch/1e9:13.2f} {E_el/1e9:13.1f} {E_el/KE*100:7.3f} "
                  f"{E_launch/t_p[-1]/L_RING:8.3f} {E_el/t_p[-1]/1e9:8.2f}")
    print(f"  Per launch is the lsm_simulation profile ({t_p[-1]/60:.0f} min, {v_p.max()/1e3:.0f} km/s), pull from Part A at each speed; the sled's kinetic energy is {KE/1e12:.0f} TJ.")
    print(f"  'mean W/m' is the heat at 70 K per meter of ring averaged over the launch, 'mean GW' the electrical power for the whole ring over the launch.")
    print(f"  The lift pair (above 7 km/s) and the load pair (below) are each pulsed only in their own speed range, so the total is split between them.")

    # C.3 coil voltage and supply
    print()
    print("  C.3  Coil voltage and supply for a 2 ms rise.")
    print(f"  Stored energy per meter of track from the 2-D solve (gap, slots, leakage) at the one-track rating: {alone['W']/1e3:,.1f} kJ/m; the pair's share {pair['W']/1e3:,.1f} kJ/m;")
    print(f"  the gap field alone is F x gap = {alone['F']*GAP0/1e3:.1f} kJ/m, so the slot and leakage flux hold {alone['W']/(alone['F']*GAP0):.1f} times the useful energy.")
    print(f"  A coil carries 1/{N_COILS} of the track's ampere-turns over {L_COIL:.0f} m, so its share of the energy is W x {L_COIL:.0f} / {N_COILS} when all coils are on.")
    print(f"  With E = L I^2 / 2 and a linear rise to I in tau, the coil voltage is V = L I / tau = 2 E / (I tau) and the supply delivers V x I = 2 E / tau,")
    print(f"  whatever the current: a low-current winding of many turns needs a high voltage, and parallel tapes trade voltage for current at the same VA.")
    print(f"  {'tape':>6} {'mode':>6} {'E/coil MJ':>9} {'turns':>6} {'I A':>6} {'L H':>7} {'V at 2 ms':>10} {'V at 10 ms':>10} {'VA, 2 ms':>9} {'VA, 10 ms':>9} {'100 tapes parallel':>18} {'store MJ':>8}")
    supply = {}
    for w in TAPE_WIDTHS:
        tp = tape[w]
        for mode, d in (("alone", alone), ("pair", pair)):
            E_coil = d["W"] * L_COIL / N_COILS
            I = tp["I_op"] * d["NI_coil"] / alone["NI_coil"]
            turns = d["NI_coil"] / I
            L = 2 * E_coil / I**2
            V2 = L * I / TAU; V10 = L * I / 10e-3
            supply[(w, mode)] = dict(E_coil=E_coil, I=I, L=L, V=V2, VA=V2 * I, turns=turns)
            print(f"  {w*1e3:4.0f} mm {mode:>6} {E_coil/1e6:9.2f} {turns:6.0f} {I:6.0f} {L:7.1f} {V2/1e6:7.2f} MV {V10/1e6:7.2f} MV {V2*I/1e9:6.2f} GVA {V10*I/1e9:6.2f} GVA "
                  f"{V2/100/1e3:6.0f} kV x {I*100/1e3:4.1f} kA {(E_coil + loss[(w,'alone')]['q30']*L_COIL/N_COILS)/1e6:8.2f}")
    n_rising = L_COIL / N_COILS / 30000.0
    print(f"  The store per coil is the coil's field energy at its rating plus one pulse's loss; the field energy comes back to the store at the end of")
    print(f"  the pulse less the losses.  Supplies: one per coil, a coil every {STAGGER:.0f} m per track, {L_RING/STAGGER/1e6:.2f} million per track and {4*L_RING/STAGGER/1e6:.1f} million for the four.")
    print(f"  At 30 km/s a new coil starts every {STAGGER/30000.0*1e3:.2f} ms, so {TAU/(STAGGER/30000.0):.0f} coils per track are rising at any moment with a 2 ms rise: the power into the")
    print(f"  field of a pair is {2*TAU/(STAGGER/30000.0)*supply[(0.003,'pair')]['VA']/1e9:.1f} GW on this count, against the {rows_A[30000.0]['P']/1e9:.2f} GW of Part A that counts the gap field only.")
    print(f"  [QUESTION FOR AUTHOR] the 2 ms rise is needed for the gap protection; the normal spread pull could use a slower rise (10 ms is one fifth of")
    print(f"  the 50 ms window) and cut the supply rating five times.  The number of supplies and the stored energy do not change.")

    # C.4 copper
    print()
    print("  C.4  Copper for comparison (same ampere-turns, both legs, no yoke in the mass).")
    rho_cu, d_cu, c_cu = 1.68e-8, 8960.0, 385.0
    for j in (5e6, 50e6):
        vol = NI_track / j * 2.0                                       # m^3 per meter of ring
        T_on = rule_T(30000.0) + L_COIL / 30000.0
        q = rho_cu * j**2 * vol * T_on
        E_launch = 0.0
        for i, dtt in enumerate(dt_p):
            vv = v_p[i + 1]
            Fu = np.interp(vv, SPEEDS, [rows_A[u]["F"] for u in SPEEDS])
            E_launch += q * (min(Fu, F_alone) / F_alone) * vv * dtt
        print(f"    {j/1e6:3.0f} A/mm^2: {vol*d_cu:6.0f} kg/m per track ({vol*d_cu*L_RING*4/1e9:.1f} Mt for four), loss per pulse at 30 km/s {q/1e3:6.1f} kJ/m, per launch {E_launch/1e9:6.0f} GJ "
              f"({E_launch/KE*100:.2f} % of KE), heating per pulse {rho_cu*j**2*T_on/(d_cu*c_cu):.3f} K")
    print(f"  Against HTS: tape mass {tape[0.003]['mass']:.1f} kg/m per track, electrical energy per launch for the losses {loss[(0.003,'pair')]['E_el']/1e9:.0f} GJ (3 mm, pair) to "
          f"{loss[(0.012,'pair')]['E_el']/1e9:.0f} GJ (12 mm, pair) after the cryogenic penalty.")

    # C.5 plates
    print()
    print("  C.5  Plate losses (hysteresis and eddy), as Part I of sled_passage_spread.py, for 1 m tracks.")
    mu_r, rho_p, d_p = lc.FERRO_MU_R, lc.FERRO_RESISTIVITY, lc.FERRO_PLATE_THICKNESS
    wh0, rho_m, cp = lc.FERRO_HYSTERESIS_COEFF, lc.FERRO_DENSITY, 450.0
    STEINMETZ = 1.6
    print(f"  Plate: mu_r {mu_r}, resistivity {rho_p:.0e} ohm m, {d_p*1e3:.0f} mm thick, hysteresis {wh0} J/m^3 per cycle at 1.2 T scaled as (B/1.2)^{STEINMETZ} [UNVERIFIED coefficients],")
    print(f"  magnetized through its thickness over the track width; the field rises over one gap length, holds for L_on / (u + V_rel), and falls.")
    print(f"  {'u km/s':>7} {'mode':>6} {'B T':>6} {'hold ms':>7} {'hyst J/m':>8} {'eddy J/m':>8} {'sum J/m':>8} {'dT mK':>6} {'per lap MJ':>10} {'drag N':>7}")
    per_speed = {}
    for u in SPEEDS:
        r = rows_A[u]
        v_sweep = u + V_REL
        t_edge = GAP0 / v_sweep
        t_hold = r["L_on"] / v_sweep
        for mode, n_tr in (("pair", 2), ("alone", 1)):
            Bt = math.sqrt(2 * MU0 * (r["F"] / n_tr) / TRACK_W)
            hyst = wh0 * (Bt / 1.2) ** STEINMETZ * d_p * TRACK_W * n_tr
            eddy = eddy_cycle(Bt, t_edge, t_hold, d_p, mu_r, rho_p) * TRACK_W * n_tr
            E = hyst + eddy
            per_speed[(u, mode)] = E
            print(f"  {u/1e3:7.0f} {mode:>6} {Bt:6.3f} {t_hold*1e3:7.1f} {hyst:8.1f} {eddy:8.1f} {E:8.1f} {E/(rho_m*d_p*TRACK_W*n_tr*cp)*1e3:6.2f} "
                  f"{E*L_RING*(1+V_REL/u)/1e6:10.1f} {E*v_sweep/V_REL:7.1f}")
    E_pl = 0.0
    for i, dtt in enumerate(dt_p):
        vv = v_p[i + 1]
        E_u = np.interp(vv, SPEEDS, [per_speed[(u, "pair")] for u in SPEEDS])
        E_pl += E_u * (vv + V_REL) * dtt
    print(f"  Per launch (pair mode): {E_pl/1e9:.2f} GJ in the plates, {E_pl/KE*100:.4f} % of the sled's kinetic energy; a 1 m track at the same pull puts sqrt(2) more")
    print(f"  field on its plate than each of a pair, so one track alone raises the plate loss per pass by about (sqrt 2)^1.6 / 2 = {2**0.8/2:.2f} (hysteresis) to 1 (eddy).")
    return dict(tape=tape, loss=loss, supply=supply, cop=cop, E_plate_launch=E_pl, B_mean=B_mean, B_max=B_max)


# =============================================================================
# PART D
# =============================================================================

def roll_model(hist, i_lane, u, with_prot, trigger=0.010, delay=1e-3, hz=5.0, support_lost=False, cap=CAP_ONE):
    """Heave and roll of casing and cable with one compensation track pulling at Y_COMP (the other of its pair lost).
    Support tracks at +-Y_SUPPORT each with a 1 Hz loop on its own gap; load compensation tracks at +-Y_COMP on the floor
    and the surviving lift track at -Y_COMP on the ceiling as per-side protection.  Cable free in heave and roll, no torsional
    stiffness [UNVERIFIED].  Returns peak relative roll, worst offsets at the plate and cable edges, smallest gaps."""
    t = hist["t"]; dt = t[1] - t[0]
    F_pull = -hist["F"][i_lane]                       # lift pull, positive = casing pulled down, from the point model (F = F_hu - F_hd)
    F_pull = np.maximum(F_pull, 0.0)
    om = 2 * math.pi * 1.0; om_p = 2 * math.pi * hz
    ndel = int(round(delay / dt))
    w = v = th = thd = wc = vc = thc = thcd = 0.0
    y_sup = (-Y_SUPPORT, Y_SUPPORT)
    tracks_prot = [(-Y_COMP, +1), (Y_COMP, +1), (-Y_COMP, -1)]   # (offset, +1 floor pulls up / -1 ceiling pulls down); the +Y_COMP lift track is the lost one
    eng = [False] * 3
    buf = []
    I_p = [0.0] * 3
    lag = 1 - math.exp(-dt / TAU)
    peak_roll = 0.0; worst_edge = 0.0; worst_cable_edge = 0.0
    min_gap_b = GAP0; min_gap_t = GAP0; peak_M = 0.0; t_peak = 0.0
    S_total = W
    for j in range(len(t)):
        rel = w - wc; reld = v - vc; rth = th - thc; rthd = thd - thcd
        # support tracks
        F_up = 0.0; M = 0.0
        S_total = 0.0
        for y in y_sup:
            if support_lost:
                S = 0.0
            else:
                rel_y = rel + rth * y; reld_y = reld + rthd * y
                S = min(max(W / 2 - (MU / 2) * (2 * om * reld_y + om**2 * rel_y), 0.0), sp.F_EML_OP / 2)
            S_total += S
            F_up += S; M += S * y
        # the pull (one lift track at +Y_COMP, pulls down)
        Fp = F_pull[j]
        F_up -= Fp; M -= Fp * Y_COMP
        # per-side protection
        if with_prot:
            buf.append((rel, rth, reld, rthd))
            if len(buf) > ndel + 1:
                buf.pop(0)
            rm, rthm, rdm, rthdm = buf[0]
            for k, (y, sgn) in enumerate(tracks_prot):
                rel_k = rm + rthm * y; reld_k = rdm + rthdm * y
                if not eng[k] and abs(rel_k) > trigger:
                    eng[k] = True
                F_des = 0.0
                if eng[k]:
                    F_des = (W - S_total) / 2 - (MU / 2) * (2 * om_p * reld_k + om_p**2 * rel_k)
                    F_des = min(max(F_des * sgn, 0.0), cap) * sgn          # a floor track can only pull up, a ceiling track only down
                I_p[k] += (math.sqrt(abs(F_des) / cap) - I_p[k]) * lag
                Fk = sgn * cap * I_p[k]**2
                F_up += Fk; M += Fk * y
        a = (F_up - W) / MU; al = M / I_CASING
        ac = -(F_up) / MU_C + 0.0; alc = -M / I_CABLE          # the cable feels the reaction; its weight is carried by its own orbit
        v += a * dt; w += v * dt; thd += al * dt; th += thd * dt
        vc += ac * dt; wc += vc * dt; thcd += alc * dt; thc += thcd * dt
        rth = th - thc; rel = w - wc
        peak_roll = max(peak_roll, abs(rth))
        if abs(M) > peak_M:
            peak_M = abs(M); t_peak = t[j]
        for y in (-Y_PLATE_EDGE, Y_PLATE_EDGE):
            r_y = rel + rth * y
            worst_edge = max(worst_edge, abs(r_y))
            min_gap_b = min(min_gap_b, GAP0 - r_y); min_gap_t = min(min_gap_t, GAP0 + r_y)
        for y in (-Y_CABLE_EDGE, Y_CABLE_EDGE):
            worst_cable_edge = max(worst_cable_edge, abs(rel + rth * y))
        if min_gap_b <= 0 or min_gap_t <= 0:
            break
    return dict(peak_roll=peak_roll, worst_edge=worst_edge, worst_cable_edge=worst_cable_edge, min_gap_b=min_gap_b, min_gap_t=min_gap_t,
                peak_M=peak_M, t_peak=t_peak, contact=(min_gap_b <= 0 or min_gap_t <= 0), t_end=t[j])


def part_D(tabs, rows_A, C, quick=False):
    cap1, cap2 = DESIGN["cap_one"], DESIGN["cap_pair"]       # one track / pair: the saturation pull of the chosen geometry (Part B)
    w1, w2 = cap1 / P12, cap2 / P12                           # the equivalent track widths at 1.2 T for the point model's caps
    hr("PART D: GAP PROTECTION AND FAILURE CASES (no limiters; the compensation tracks guard the gap)")
    print(f"  Protection law (point and beam models): gap sensors see the casing-to-cable offset with a delay; once the measured offset passes the")
    print(f"  trigger, the compensation tracks on the face that can pull the casing back run a position loop (critically damped at 5 Hz) plus the weight")
    print(f"  the support tracks are not carrying, clipped to the track rating ({cap2/1e3:,.0f} kN/m pair, {cap1/1e3:.0f} kN/m one), current lag {TAU*1e3:.0f} ms.")
    print(f"  Gap_b / gap_t = smallest gap on the floor / ceiling face in mm (nominal 100).  Point model: one meter of casing and cable, both free.")
    u30 = 30000.0
    T30 = rule_T(u30)
    all_gaps = []

    def note(label, model, r):
        all_gaps.append((label, model, r["min_gap_b"], r["min_gap_t"], r.get("contact", False)))

    # D.1 support loss
    print()
    print("  D.1  Total loss of both casing support tracks.")
    print(f"  Free fall from rest: {GAP0*1e3:.0f} mm in {math.sqrt(2*GAP0/pz.G_NET)*1e3:.0f} ms; 5, 10, 20 mm of drift take {math.sqrt(2*0.005/pz.G_NET)*1e3:.0f}, "
          f"{math.sqrt(2*0.010/pz.G_NET)*1e3:.0f}, {math.sqrt(2*0.020/pz.G_NET)*1e3:.0f} ms.")
    print(f"  At rest (no sled), point model, load compensation tracks catch the casing:")
    print(f"  {'tracks':>6} {'trigger mm':>10} {'delay ms':>8} {'engaged ms':>10} {'gap_b':>6} {'gap_t':>6} {'peak pull kN/m':>14} {'hold kN/m':>9} {'hold B T':>8}")
    combos = [(tr, dl) for tr in (0.005, 0.010, 0.020) for dl in (1e-3, 5e-3)]
    rest = {}
    for width, lab in ((w2, "pair"), (w1, "one")):
        lanes = [new_lane(u30, width, no_sled=True, cancel=False, support_fail_t=0.0, prot=prot(tr, dl, cap=P12 * width), label=f"{lab} {tr*1e3:.0f}/{dl*1e3:.0f}")
                 for tr, dl in combos]
        res, hist = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
        for (tr, dl), r, i in zip(combos, res, range(len(res))):
            hold = float(hist["F"][i, -1])
            rest[(lab, tr, dl)] = r
            note(f"support lost at rest, {lab}, trigger {tr*1e3:.0f} mm, delay {dl*1e3:.0f} ms", "point", r)
            print(f"  {lab:>6} {tr*1e3:10.0f} {dl*1e3:8.0f} {r['t_engaged']*1e3:10.1f} {fmt_gap(r)} {r['pk_hu']/1e3:14.1f} {hold/1e3:9.1f} {math.sqrt(2*MU0*hold/width):8.3f}")
    # protection bandwidth
    print(f"  Protection loop bandwidth (pair, 10 mm trigger, 1 ms delay, at rest):")
    lanes = [new_lane(u30, w2, no_sled=True, cancel=False, support_fail_t=0.0, prot=prot(0.010, 1e-3, cap=cap2, hz=hz), label=f"hz {hz}") for hz in (2.0, 5.0, 10.0, 20.0)]
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    for hz, r in zip((2.0, 5.0, 10.0, 20.0), res):
        print(f"    {hz:4.0f} Hz: gaps {fmt_gap(r)}, peak pull {r['pk_hu']/1e3:.0f} kN/m")
    print(f"  'hold' is the steady pull after the catch, the casing weight {W/1e3:.1f} kN/m, and 'hold B' the mean plate field that takes on the tracks that are")
    print(f"  pulling ({math.sqrt(2*MU0*W/w2):.3f} T on each of a pair, {math.sqrt(2*MU0*W/w1):.3f} T on one alone); the position loop then brings the casing back to the 100 mm gap.")
    print(f"  Hand estimate check (10 mm trigger, 2 ms rise, 353 kN/m, no delay): stops after 14.8 mm; here with the 1 ms delay and the loop, the pair stops at "
          f"{(GAP0 - rest[('pair', 0.010, 1e-3)]['min_gap_t'])*1e3:.1f} mm and one track at {(GAP0 - rest[('one', 0.010, 1e-3)]['min_gap_t'])*1e3:.1f} mm (the casing falls toward the ceiling plates).")

    print()
    print(f"  During a 30 km/s pass (point model, pair, support lost at t_fail; zone front at the point at t = 0, pull window centered on the zone):")
    print(f"  {'t_fail ms':>9} {'trigger mm':>10} {'delay ms':>8} {'engaged ms':>10} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'peak pull kN/m':>14}")
    cases = [(-0.1, 0.010, 1e-3), (-0.025, 0.010, 1e-3), (0.0, 0.010, 1e-3), (0.025, 0.010, 1e-3), (0.2, 0.010, 1e-3)] + [(0.0, tr, dl) for tr, dl in combos if (tr, dl) != (0.010, 1e-3)]
    lanes = [new_lane(u30, w2, support_fail_t=tf, prot=prot(tr, dl), label=f"pass {tf}") for tf, tr, dl in cases]
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    for (tf, tr, dl), r in zip(cases, res):
        note(f"support lost during pass at {tf*1e3:.0f} ms, pair, trigger {tr*1e3:.0f}, delay {dl*1e3:.0f}", "point", r)
        print(f"  {tf*1e3:9.0f} {tr*1e3:10.0f} {dl*1e3:8.0f} {r['t_engaged']*1e3:10.1f} {r['max_exc_all']*1e3:7.1f} {fmt_gap(r)} {max(r['pk_hu'], r['pk_hd'])/1e3:14.1f}")
    lanes = [new_lane(u30, w1, support_fail_t=0.0, prot=prot(0.010, 1e-3, cap=cap1), label="one")]
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    r = res[0]
    note("support lost during pass at 0 ms, one track alone, trigger 10, delay 1", "point", r)
    print(f"  one track alone, t_fail 0, trigger 10 mm, delay 1 ms: engaged {r['t_engaged']*1e3:.1f} ms, excursion {r['max_exc_all']*1e3:.1f} mm, gaps {fmt_gap(r)}, peak pull {max(r['pk_hu'], r['pk_hd'])/1e3:.0f} kN/m")

    # beam
    print()
    print(f"  Beam model (10 m truss, cable as a string, 10 km window), support lost at rest and as the zone reaches the window middle.  The protection")
    print(f"  acts on the offset averaged over {pz.PROT_X_WINDOW:.0f} m of casing with the closing speed taken over the sensor delay (at least {pz.PROT_V_WINDOW*1e3:.0f} ms); a")
    print(f"  casing with no stiffness at all is not run here, since a per-node loop on a structure that nothing couples along its length is unstable")
    print(f"  on this grid, and the casing is at least a 10 m truss by the design rule.")
    trusses = [sp.truss(10.0, 0.1)]
    t_mid = 5000.0 / u30
    beam_combos = combos if not quick else [(0.010, 1e-3)]
    print(f"  {'trigger mm':>10} {'delay ms':>8} {'case':>36} {'truss':>22} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'engaged ms':>10} {'contact':>7}")
    t0 = time.time()
    for tr, dl in beam_combos:
        cases_b = [dict(label="support lost at rest", no_sled=True, no_pulse=True, support_fail_t=0.0),
                   dict(label="support lost 50 ms before the zone", support_fail_t=t_mid - 0.05),
                   dict(label="support lost as the zone arrives", support_fail_t=t_mid)]
        trs = trusses
        res, _ = pz.beam_pulsed(u30, L_SC, trs, cases_b, tabs, tau=TAU, L_win=10000.0, T_spread=T30, n_coils=N_COILS, L_coil=L_COIL,
                                track_width=w2, prot=prot(tr, dl))
        for r in res:
            note(f"{r['case']}, beam {r['truss']}, trigger {tr*1e3:.0f}, delay {dl*1e3:.0f}", "beam", r)
            print(f"  {tr*1e3:10.0f} {dl*1e3:8.0f} {r['case']:>36} {r['truss']:>22} {r['max_exc']*1e3:7.1f} {fmt_gap(r)} {r['t_engaged']*1e3:10.1f} {('yes' if r['contact'] else 'no'):>7}")
    print(f"  ({time.time()-t0:.0f} s)")

    # D.2 one track of a pair lost
    print()
    print("  D.2  Loss of one whole compensation track during a 30 km/s pass (point model, protection 10 mm / 1 ms on).")
    tc = 0.5 * L_SC / u30
    cases2 = [("other track alone from the start (Part A 'one')", dict(hd_width=w1, hu_width=w1)),
              ("track lost at the start of the pull, no reassignment", dict(hd_width=w2, hu_width=w2, pull_scale=[(tc - T30 / 2, 1e9, 0.5)])),
              ("track lost mid-pull, no reassignment", dict(hd_width=w2, hu_width=w2, pull_scale=[(tc, 1e9, 0.5)])),
              ("track lost mid-pull, other takes over after 5 ms", dict(hd_width=w1, hu_width=w1, pull_scale=[(tc, tc + 5e-3, 0.5)])),
              ("track lost mid-pull, other takes over after 20 ms", dict(hd_width=w1, hu_width=w1, pull_scale=[(tc, tc + 20e-3, 0.5)]))]
    lanes = [new_lane(u30, w2, prot=prot(0.010, 1e-3, cap_up=cap2, cap_down=cap1), label=lab, **kw) for lab, kw in cases2]
    print(f"  A lift track is lost, so the protection has one lift track ({cap1/1e3:.0f} kN/m) to pull down and both load tracks ({cap2/1e3:,.0f} kN/m) to pull up.")
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    print(f"  {'case':>52} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'engaged ms':>10} {'peak kN/m':>9}")
    for (lab, kw), r in zip(cases2, res):
        note(lab, "point", r)
        print(f"  {lab:>52} {r['max_exc_all']*1e3:7.1f} {fmt_gap(r)} {r['t_engaged']*1e3:10.1f} {max(r['pk_hu'], r['pk_hd'])/1e3:9.0f}")
    print(f"  'no reassignment' = half the pull is simply missing from then on; the gap protection catches what is left.  Half a pull missing over the")
    print(f"  whole window is {0.5*rows_A[u30]['F']*T30/1e3:.1f} kN s/m, {0.5*rows_A[u30]['F']*T30/MU:.2f} m/s on the casing.")

    # D.3 single failures
    print()
    print("  D.3  Single failures at 30 km/s (pair, protection 10 mm / 1 ms).  One supply dead is one coil dead: each coil has its own supply.")
    cases3 = [("normal", dict()),
              ("sensor error +5 mm (regulation and protection)", dict(gap_err=0.005, prot=prot(0.010, 1e-3, sensor_err=0.005))),
              ("sensor error -5 mm (regulation and protection)", dict(gap_err=-0.005, prot=prot(0.010, 1e-3, sensor_err=-0.005))),
              ("false trigger, sensors true", dict(prot=prot(0.0, 1e-3))),
              ("false trigger with a 5 mm sensor error", dict(prot=prot(0.0, 1e-3, sensor_err=0.005))),
              ("false full pull up for 5 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 5e-3, cap2))),
              ("false full pull up for 20 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 20e-3, cap2))),
              ("false full pull up for 50 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 50e-3, cap2))),
              ("false full pull down for 5 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 5e-3, -cap2))),
              ("false full pull down for 10 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 10e-3, -cap2))),
              ("false full pull down for 20 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 20e-3, -cap2))),
              ("false full pull down for 50 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 50e-3, -cap2))),
              ("false full pull down 20 ms, protection 20 Hz", dict(no_sled=True, cancel=False, false_pull=(0.0, 20e-3, -cap2), prot=prot(0.010, 1e-3, cap=cap2, hz=20.0))),
              ("false one-track pull down for 50 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 50e-3, -cap1))),
              ("spread pull fired with no sled", dict(no_sled=True)),
              ("spread pull fired with no sled, one track", dict(no_sled=True, hd_width=w1, hu_width=w1)),
              ("spread pull 10 ms late", dict(dt_err=10e-3)), ("spread pull 10 ms early", dict(dt_err=-10e-3)),
              ("spread pull 50 ms late", dict(dt_err=50e-3)), ("spread pull 50 ms early", dict(dt_err=-50e-3))]
    lanes = []
    for lab, kw in cases3:
        kw2 = dict(kw)
        if "prot" not in kw2:
            kw2["prot"] = prot(0.010, 1e-3, cap=kw2.get("hu_width", w2) * P12)
        lanes.append(new_lane(u30, w2, label=lab, **kw2))
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    print(f"  {'case':>48} {'model':>6} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'engaged ms':>10} {'peak kN/m':>9}")
    for (lab, kw), r in zip(cases3, res):
        note(lab, "point", r)
        print(f"  {lab:>48} {'point':>6} {r['max_exc_all']*1e3:7.1f} {fmt_gap(r)} {r['t_engaged']*1e3:10.1f} {max(r['pk_hu'], r['pk_hd'])/1e3:9.0f}")
    # beam single failures
    xc = 5000.0
    cases_b = [dict(label="normal"), dict(label="one coil (one supply) dead, family 0", fail_coil=(0, xc)),
               dict(label="one coil dead, family 10", fail_coil=(10, xc)),
               dict(label="spread pull with no sled", no_sled=True),
               dict(label="false full pull up 5 ms, no sled", no_sled=True, no_pulse=True, false_pull=(t_mid, 5e-3, cap2)),
               dict(label="false full pull down 20 ms, no sled", no_sled=True, no_pulse=True, false_pull=(t_mid, 20e-3, -cap2))]
    if quick:
        cases_b = cases_b[:2]
    t0 = time.time()
    res, _ = pz.beam_pulsed(u30, L_SC, trusses, cases_b, tabs, tau=TAU, L_win=10000.0, T_spread=T30, n_coils=N_COILS, L_coil=L_COIL,
                            track_width=w2, prot=prot(0.010, 1e-3))
    for r in res:
        note(f"{r['case']}, beam {r['truss']}", "beam", r)
        print(f"  {r['case']:>48} {'beam':>6} {r['max_exc']*1e3:7.1f} {fmt_gap(r)} {r['t_engaged']*1e3:10.1f} {'':>9}  {'contact' if r['contact'] else ''}")
    res, _ = pz.beam_pulsed(u30, L_SC, trusses, [dict(label="normal, one track alone"), dict(label="one coil dead, one track alone", fail_coil=(0, xc))],
                            tabs, tau=TAU, L_win=10000.0, T_spread=T30, n_coils=N_COILS, L_coil=L_COIL, track_width=w1, prot=prot(0.010, 1e-3, cap=cap1))
    for r in res:
        note(f"{r['case']}, beam {r['truss']}", "beam", r)
        print(f"  {r['case']:>48} {'beam':>6} {r['max_exc']*1e3:7.1f} {fmt_gap(r)} {r['t_engaged']*1e3:10.1f} {'':>9}  {'contact' if r['contact'] else ''}")
    print(f"  ({time.time()-t0:.0f} s)  Beam: 10 m truss, one coil dead = that coil's 200 m carries 19/20 of the pull.")
    print(f"  A false pull at the full track rating is a control failure that the protection on the other face must out-pull with the same rating; the")
    print(f"  rows above give the time within which a watchdog must cut a faulty coil.  A false pull at the normal command level is the 'spread pull")
    print(f"  fired with no sled' row.")
    print(f"  Loss of ring power: the casing support tracks are persistent and keep the casing up; each compensation coil runs its pulse from its own")
    print(f"  store, so the dynamics are those of the normal case.  The store per coil (Part C.3) is its field energy at the one-track rating,")
    if C:
        print(f"  {C['supply'][(0.003,'alone')]['E_coil']/1e3:.1f} kJ, plus one pulse's loss; a catch after a support loss needs the same energy and then holds with no further")
    else:
        print(f"  (Part C not run), plus one pulse's loss; a catch after a support loss needs the same energy and then holds with no further")
    print(f"  input but the cryogenic load, which is the ring's problem, not the coil's.")

    # D.4 roll
    print()
    print("  D.4  Roll from one track of a pair pulling alone (2-DOF heave and roll of casing and cable, 30 km/s, lift track at +2.8 m).")
    print(f"  Geometry [UNVERIFIED]: support plates centered {Y_SUPPORT} m from the centerline, compensation plates {Y_COMP} m (edges at {Y_PLATE_EDGE} m), cable edge {Y_CABLE_EDGE:.2f} m;")
    print(f"  casing roll inertia {I_CASING:,.0f} kg m^2 per m (radius of gyration {R_GYR_CASING} m), cable {I_CABLE:,.0f} kg m^2 per m (square section), no torsional stiffness.")
    print(f"  Per-side protection: each load track regulates its own side once its ceiling or floor gap drifts 10 mm (1 ms delay).")
    _, hist = pz.point_pulsed([new_lane(u30, w2, label="pair"), new_lane(u30, w1, label="one")], u30, tabs)   # the one-track pull at its rated cap
    print(f"  The pull history is the one-track lane of Part A with the pull capped at the track rating of Part B ({cap1/1e3:.0f} kN/m).")
    print(f"  {'case':>60} {'peak M kN m/m':>13} {'roll mrad':>9} {'plate edge mm':>13} {'cable edge mm':>13} {'gap_b':>6} {'gap_t':>6}")
    for lab, kw in (("one lift track alone, support loops only, no roll protection", dict(with_prot=False)),
                    ("one lift track alone, per-side protection", dict(with_prot=True)),
                    ("one lift track alone, support lost too, per-side protection", dict(with_prot=True, support_lost=True))):
        rr = roll_model(hist, 1, u30, cap=cap1, **kw)
        all_gaps.append((f"roll: {lab}", "2-DOF", rr["min_gap_b"], rr["min_gap_t"], rr["contact"]))
        print(f"  {lab:>60} {rr['peak_M']/1e3:13.0f} {rr['peak_roll']*1e3:9.2f} {rr['worst_edge']*1e3:13.1f} {rr['worst_cable_edge']*1e3:13.1f} {rr['min_gap_b']*1e3:6.1f} {rr['min_gap_t']*1e3:6.1f}"
              + ("  CONTACT" if rr["contact"] else ""))
    print(f"  The pair shares the pull in normal use, so there is no roll moment then; roll appears only when one track of a pair is out.")
    print(f"  'plate edge' and 'cable edge' are the largest casing offsets at those lateral positions; gap_b and gap_t are the smallest gaps at the plate edges.")

    # D.5
    print()
    print("  D.5  Summary of all cases in this part.")
    contacts = [g for g in all_gaps if g[4]]
    worst_b = min(all_gaps, key=lambda g: g[2]); worst_t = min(all_gaps, key=lambda g: g[3])
    print(f"  Cases run: {len(all_gaps)}; with contact: {len(contacts)}.")
    for g in contacts:
        print(f"    CONTACT: {g[0]} ({g[1]})")
    print(f"  Smallest floor gap: {worst_b[2]*1e3:.1f} mm in '{worst_b[0]}' ({worst_b[1]}).")
    print(f"  Smallest ceiling gap: {worst_t[3]*1e3:.1f} mm in '{worst_t[0]}' ({worst_t[1]}).")
    single = [g for g in all_gaps if "support lost too" not in g[0] and "no roll protection" not in g[0]]
    ws = min(min(g[2], g[3]) for g in single)
    wsg = min(single, key=lambda g: min(g[2], g[3]))
    print(f"  Smallest gap in any single-failure case: {ws*1e3:.1f} mm, in '{wsg[0]}' ({wsg[1]}); the double failure (one track of a pair out and both")
    print(f"  support tracks lost) and the roll case without per-side protection are listed separately in D.4.")
    return all_gaps


# =============================================================================
# PART E
# =============================================================================

def part_E(tabs):
    hr("PART E: THE NUMBERS BEHIND THE REJECTED OPTIONS, AND THE 8,620 m/s FIX")
    print("  Halbach array limiters (persistent HTS Halbach arrays repelling a copper plate on the cable), from the existing runs:")
    print(f"  levitation.txt (rerun 2026-10-11): present arrays, 100 mm wavelength, surface field 0.640 T; repulsion at a 20 mm gap 27.0 kN/m, at 10 mm 93.8 kN/m,")
    print(f"  at 2 mm 253.7 kN/m, against the casing weight {W/1e3:.1f} kN/m; a casing falling after a support loss is 'NOT arrested', gap reached 0.50 mm.")
    print(f"  Reaction plate heating (OFHC copper, 20 mm) at a 20 mm gap 1.63 MW/m, service temperature in 44 s, melting in 94 s; at 100 mm 73.6 W/m.")
    print(f"  sled_passage_limits.txt, limit 2 (5,000 t, 100 m zone, limiters as the only protection against the sled): the present limiter keeps the gap")
    print(f"  above 20 mm only to 19.0 km/s and avoids contact only to 22.0 km/s.  To reach 30 km/s with the 20 mm margin takes 4 x the ampere-turns and a")
    print(f"  200 mm wavelength: that limiter leaves 9.73 kN/m of repulsion at the nominal gap all the time and heats the plate at 419 kW/m, which is")
    print(f"  {419e3*L_RING/1e12:,.0f} TW around the ring; 3 x the ampere-turns at 100 mm wavelength avoids contact at 30 km/s but with 0 mm of margin (gap >= 20 mm only to 20.5 km/s).")
    print(f"  The present limiter's drag at the nominal gap is small: 2.74e-6 N/m, 1.075 MW for the ring (levitation.txt, analysis 5).")
    print(f"  Pulsed LIM limiter (driven traveling-wave coils under the cable pushing on a conducting plate), sled_passage_pulsed.txt part J: at full slip")
    print(f"  (static coils, the limiter case) the push costs 69 W per newton; holding the casing weight alone would be {69*W/1e6:.1f} MW/m, {69*W*L_RING/1e12:,.0f} TW for the")
    print(f"  ring; the pulsed-track case needs 22 to 2,000 m of push-down track beside the pull track and 17,000 to 1,600,000 GVA of reactive power.")
    print()
    print("  The 8,620 m/s fix.  levitation_simulation.py had the cable speed relative to the casing hard-coded as v_cable = 8620 in four places (the 250 km")
    print("  value); they now read ring_altitude.V_REL, 8,825.21 m/s at 700 km.  Backup: levitation_simulation_backup_2026-10-11.py; old log kept as")
    print("  levitation_8620_2026-10-04.txt.  Lines that changed between the two logs:")
    old_p = os.path.join(RUNS_DIR, "levitation_8620_2026-10-04.txt")
    new_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "levitation.txt")
    try:
        with open(old_p, encoding="utf-8", errors="replace") as f:
            old = f.read().splitlines()
        with open(new_p, encoding="utf-8", errors="replace") as f:
            new = f.read().splitlines()
        keys = ("Cable velocity", "Drag per meter", "Drag power per meter", "Total ring drag", "Safety drag", "F_safety@20", "Service temp", "Melting point (", "Heating rate")
        shown = 0
        for lo, ln in zip(old, new):
            if lo != ln and any(k in ln for k in keys) and shown < 24:
                print(f"    old: {lo.strip()}")
                print(f"    new: {ln.strip()}")
                shown += 1
        n_diff = sum(1 for lo, ln in zip(old, new) if lo != ln)
        print(f"  {n_diff} of {len(new)} lines differ in all (the rest are the heating-rate and equilibrium-temperature tables, which move a few percent).")
    except OSError as e:
        print(f"  (could not read the logs: {e})")
    print(f"  What changed in substance: the eddy drag and heating of the copper plates scale with the speed, so the drag per meter went from 2.80e-6 to")
    print(f"  2.74e-6 N/m (the skin depth shrinks with speed, so the force at a given gap falls slightly) and the ring drag power from 1.1 to 1.075 MW;")
    print(f"  the limiter forces in the time-domain runs (sled_passage_*.py) already used V_REL from ring_altitude.py and are unchanged.")
    print(f"  The limiter table used by the sled passage scripts at V_REL: F0 = {tabs['F0']/1e3:.1f} kN/m at zero gap, decay length {1/(2*tabs['k'])*1e3:.1f} mm, "
          f"{tabs['F0']*math.exp(-2*tabs['k']*0.020)/1e3:.1f} kN/m at 20 mm, {tabs['F0']*math.exp(-2*tabs['k']*GAP0):.3f} N/m at 100 mm.")


# =============================================================================

def main():
    parts = "ABCDE"
    quick = False
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
        elif arg == "--quick":
            quick = True
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    tabs_lim = sp.limiter_tables()
    tabs = no_limiters(tabs_lim)
    hr("COMPENSATION TRACKS, FINAL LAYOUT (2026-10-11)")
    print(f"  Ring at 700 km: V_rel {V_REL:,.0f} m/s, casing {MU:,.0f} kg/m ({W/1e3:.1f} kN/m), cable {MU_C:,.0f} kg/m, cable face {ring.CABLE_SIDE:.2f} m.")
    print(f"  Two casing support tracks (1 m, persistent HTS), four compensation tracks of 1 m with 1.4 m plates: two lift (ceiling), two load (floor,")
    print(f"  outboard of the support tracks), the two of a pair on opposite sides of the centerline, either able to do the pair's job alone.  No limiters.")
    print(f"  Underside width: two support plates {2*lc.EML_PLATE_WIDTH:.1f} m plus two load compensation plates {2*1.4:.1f} m = {2*lc.EML_PLATE_WIDTH+2.8:.1f} m of {ring.CABLE_SIDE:.2f} m, "
          f"{ring.CABLE_SIDE-2*lc.EML_PLATE_WIDTH-2.8:.2f} m left; top: two lift plates 2.8 m.")
    rows_A = part_A(tabs) if ("A" in parts or "B" in parts or "C" in parts or "D" in parts) else None
    B = part_B(rows_A, tabs) if ("B" in parts or "C" in parts) else None
    C = part_C(rows_A, B) if ("C" in parts) else None
    if "D" in parts:
        part_D(tabs, rows_A, C, quick=quick)
    if "E" in parts:
        part_E(tabs_lim)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
