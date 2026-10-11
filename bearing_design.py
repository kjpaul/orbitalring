#!/usr/bin/env python3
"""
Bearing field design, verification and sizing (Paul, 2026-10-12).

Parts (Claude_Code_instructions_bearing_design_2026-10-12.md):
  A  solver verification (the solver checks and the scikit-fem cross-check are in bearing_fem_check.py; here the plate
     B-H curve and the saturation sensitivity)
  B  pole geometry sweep: pole width, slot width, slot depth, number of poles, plate thickness and layout, track width
  C  casing support tracks sized for 95 kN/m alone with margin, 12 mm tape in persistent mode
  D  compensation tracks sized for one track alone (393 kN/m lift, the low-speed load case, the parked spacecraft),
     coil length and rise time
  E  tape losses with the book's own formulas (lim_physics.py) against the strip estimate
  F  the dynamics rerun with the chosen design

Usage:
    python bearing_design.py                 # everything (the sweep is cached in output/bearing_sweep.json)
    python bearing_design.py --parts=B       # some parts
    python bearing_design.py --resweep       # discard the cached sweep

Models: magnetostatic_2d.py (sparse direct solve, saturating iron), compensation_fields.py (track builder),
sled_passage_pulsed.py (point and beam dynamics), lim_physics.py (HTS loss formulas).
Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import csv
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
import compensation_fields as cf
import magnetostatic_2d as ms

MU0 = lc.MU0
MU, MU_C, W, GAP0, BAND = pz.MU, pz.MU_C, pz.W, pz.GAP0, pz.BAND
L_SLED, M_SLED_M, M_SC0 = pz.L_SLED, pz.M_SLED_M, pz.M_SC0
L_RING = ring.L_RING
V_REL = pz.V_REL
P12 = pz.P_SAT[1.2]
CABLE_FACE = ring.CABLE_SIDE                    # 8.24 m
L_SC = 100.0
TAU_HW = 2e-3                                   # s, the hardware rise the protection needs
MATCH = pz.MATCH_RING
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
SWEEP_JSON = os.path.join(OUT, "bearing_sweep.json")
RUNS_DIR = os.path.join(HERE, "..", "..", "Vol_III", "runs_700km_2026-10-04")
THEVA_CSV = os.path.join(HERE, "..", "..", "Vol_III", "bin", "THEVA Pro-Line Advanced Pinning 2G HTS 70 K Field Dependence.csv")
RHO_LAM = 7650.0                                # kg/m^3, non-oriented silicon steel laminations [UNVERIFIED: typical 7,650]
RHO_PLATE = lc.FERRO_DENSITY                    # 7,870 kg/m^3 (iron; the resin composite would be lighter) [UNVERIFIED]
TAPE_T = 0.08e-3; TAPE_RHO = 8890.0             # [UNVERIFIED]
LEVELS = (25e3, 50e3, 75e3, 100e3, 150e3, 200e3, 300e3)
H_SWEEP = 0.01
EFF_RATE = 0.5                                  # rating: the pull where the saturating force is half the linear force
F_LIFT_ALONE = 393e3                            # the 30 km/s pull (compensation_tracks.txt, Part A)
F_LOAD_ALONE = 353e3                            # the 1 km/s pull with short coils (sled_passage_spread.txt rule 4 at 1 km/s)
F_PARKED = M_SC0 * pz.G_NET / L_SC              # 396 kN/m, a parked 5,000 t on 100 m


def hr(t=""):
    pz.hr(t)


# =============================================================================
# PART B: THE GEOMETRY SWEEP
# =============================================================================

def rate_curve(curve, eff=EFF_RATE):
    """From a track_curve() output: the pull where the saturating force has fallen to `eff` of the linear force, with the
    ampere-turns and the other quantities interpolated there; plus the pull at 0.7 efficiency and the curve maximum."""
    rows = [r for r in curve if r.get("nonlinear")]
    NIs = np.array([r["NI"] for r in rows]); Fs = np.array([r["F"] for r in rows]); effs = Fs / np.array([r["F_lin"] for r in rows])
    out = dict(F_max=float(Fs.max()), NI_max=float(NIs[Fs.argmax()]))
    for tag, e in (("50", eff), ("70", 0.7)):
        if effs.min() > e:
            NI_e = float(NIs[-1])
        elif effs.max() < e:
            NI_e = float(NIs[0])
        else:
            # efficiency falls with NI; find the first crossing
            k = int(np.argmax(effs < e))
            k0 = max(k - 1, 0)
            NI_e = float(np.interp(e, [effs[k], effs[k0]], [NIs[k], NIs[k0]])) if k > 0 else float(NIs[0])
        d = {}
        for key in ("F", "B_track", "B_pole", "B_plate_max", "plate_sat_frac", "W_mag", "Bx_wm", "By_wm", "Bx_w", "By_w"):
            d[key] = float(np.interp(NI_e, NIs, [r[key] for r in rows]))
        d["NI"] = NI_e
        out[tag] = d
    return out


def geom_label(gm):
    return (f"W {gm['W']:.1f} m, {gm['n_full']} full + 2 half poles of {gm['pole_w']*1e3:.0f} mm, slots {gm['slot_w']*1e3:.0f} mm, depth {gm['slot_depth']*1e3:.0f} mm, "
            f"plate {gm['plate_t']*1e3:.0f} mm{' split' if gm['plate_split'] else ''}")


def sweep_geoms():
    """The list of geometries, in stages; later stages are filled from the results of the earlier ones at run time."""
    stage1 = []
    for n_full in (1, 2, 3):
        pitch = 1.0 / (n_full + 1)
        for f in (0.4, 0.5, 0.6, 0.7):
            stage1.append(dict(W=1.0, n_full=n_full, pole_w=round(f * pitch, 4), slot_w=round((1 - f) * pitch, 4), slot_depth=0.30, plate_t=0.20, plate_split=False))
    return stage1


def run_one(gm, verbose=True):
    t0 = time.time()
    curve = cf.track_curve(LEVELS, verbose=False, h=H_SWEEP, **gm)
    rt = rate_curve(curve)
    meta = {k: curve[0][k] for k in ("W", "n_full", "n_groups", "pole_w", "slot_w", "slot_depth", "plate_t", "plate_w", "plate_split", "iron_casing_m2", "iron_cable_m2", "winding_m2", "pole_frac")}
    rec = dict(gm=gm, meta=meta, rate=rt, F_lin100=curve[0]["F"], W_lin100=curve[0]["W_mag"],
               curve=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in r.items() if k not in ("strips",)} for r in curve],
               t=time.time() - t0)
    if verbose:
        r50 = rt["50"]
        print(f"    {geom_label(meta):90} rated {r50['F']/1e3:6.1f} kN/m at {r50['NI']/1e3:5.1f} kA/pole, max {rt['F_max']/1e3:6.1f}, "
              f"W/Fg {r50['W_mag']/(r50['F']*GAP0):4.1f}, winding {r50['Bx_wm']:.2f}/{r50['By_wm']:.2f} T, iron {meta['iron_casing_m2']*RHO_LAM:5.0f}+{meta['iron_cable_m2']*RHO_PLATE:5.0f} kg/m ({rec['t']:.0f} s)", flush=True)
    return rec


def run_sweep(verbose=True):
    if os.path.exists(SWEEP_JSON):
        with open(SWEEP_JSON, encoding="utf-8") as f:
            return json.load(f)
    os.makedirs(OUT, exist_ok=True)
    recs = []
    done = set()

    def key(gm):
        return (gm["W"], gm["n_full"], round(gm["pole_w"], 4), round(gm["slot_w"], 4), gm["slot_depth"], gm["plate_t"], gm["plate_split"])

    def run_list(gms, tag):
        if verbose:
            print(f"  Stage {tag}: {len(gms)} geometries", flush=True)
        for gm in gms:
            if key(gm) in done:
                continue
            recs.append(run_one(gm, verbose)); recs[-1]["stage"] = tag
            done.add(key(gm))

    t0 = time.time()
    run_list(sweep_geoms(), "1: pole count and pole fraction (1 m, 300 mm slots, 200 mm plate)")
    best = sorted(recs, key=lambda r: -r["rate"]["50"]["F"])[:3]
    stage2 = []
    for b in best:
        for pt in (0.10, 0.15, 0.20, 0.30):
            for split in (False, True):
                stage2.append(dict(b["gm"], plate_t=pt, plate_split=split))
    run_list(stage2, "2: plate thickness and layout for the best three")
    best = sorted([r for r in recs if r["gm"]["plate_t"] == 0.20 and not r["gm"]["plate_split"]], key=lambda r: -r["rate"]["50"]["F"])[:2]
    stage3 = [dict(b["gm"], slot_depth=d) for b in best for d in (0.20, 0.40)]
    run_list(stage3, "3: slot depth for the best two")
    b = best[0]["gm"]
    pitch = b["pole_w"] + b["slot_w"]
    stage4 = []
    for Wt in (1.5, 2.0):
        n_full = max(1, int(round(Wt / pitch)) - 1)
        for pt in (0.20, 0.30):
            stage4.append(dict(b, W=Wt, n_full=n_full, plate_t=pt))
    run_list(stage4, "4: wider tracks with the best pole pattern")
    if verbose:
        print(f"  sweep done, {len(recs)} geometries, {time.time()-t0:.0f} s", flush=True)
    with open(SWEEP_JSON, "w", encoding="utf-8") as f:
        json.dump(recs, f, indent=1, default=float)
    return recs


def part_B():
    hr("PART B: POLE GEOMETRY SWEEP (saturating iron, 10 mm cells, box at least 3.8 x 2.0 m)")
    print(f"  Track of width W facing its plate across 100 mm: n full poles and a half pole at each edge, slots between, windings filling the slots,")
    print(f"  laminated yoke 200 mm below; each pole wound with NI, inner slots carry 2 NI.  Plate: sintered iron in resin, mu_r {lc.FERRO_MU_R}, saturating at")
    print(f"  {lc.FERRO_B_SAT} T; yoke and poles mu_r 1,000 saturating at 1.8 T.  Current stepped 25 to 300 kA per pole with the permeability carried along.")
    print(f"  Rating = the pull where the saturating force has fallen to {EFF_RATE:.0%} of the linear force: beyond it half the ampere-turns go into saturated")
    print(f"  iron.  'W/Fg' = stored energy over F x gap (the leakage figure; 1 = gap field only).  Iron: yoke and poles at {RHO_LAM:,.0f} kg/m^3, plate at {RHO_PLATE:,.0f}.")
    t0 = time.time()
    recs = run_sweep(verbose=True)
    print(f"  ({time.time()-t0:.0f} s to run or load)")
    print()
    print(f"  {'geometry':>96} {'rated kN/m':>10} {'NI/pole kA':>10} {'B_track T':>9} {'plate sat %':>11} {'W/Fg':>5} {'Bx/By mean T':>13} {'Bx/By max T':>12} {'at 70 %':>8} {'curve max':>9} {'casing kg/m':>11} {'plate kg/m':>10}")
    for r in sorted(recs, key=lambda r: -r["rate"]["50"]["F"]):
        m = r["meta"]; r50 = r["rate"]["50"]; r70 = r["rate"]["70"]
        print(f"  {geom_label(m):>96} {r50['F']/1e3:10.1f} {r50['NI']/1e3:10.1f} {r50['B_track']:9.3f} {r50['plate_sat_frac']*100:11.0f} {r50['W_mag']/(r50['F']*GAP0):5.1f} "
              f"{r50['Bx_wm']:6.2f}/{r50['By_wm']:5.2f} {r50['Bx_w']:6.2f}/{r50['By_w']:5.2f} {r70['F']/1e3:8.1f} {r['rate']['F_max']/1e3:9.1f} "
              f"{m['iron_casing_m2']*RHO_LAM:11.0f} {m['iron_cable_m2']*RHO_PLATE:10.0f}")
    # per meter of width
    print()
    print("  Rated pull per meter of track width, and what the sweep says:")
    by_W = {}
    for r in recs:
        by_W.setdefault(round(r["meta"]["W"], 1), []).append(r)
    print("  (the 'split' rows repeat the unsplit ones: with three poles the plate is one strip per pole pair already, so the split does nothing;")
    print("  a split plate only differs for four poles and more, which the sweep did not reach because the three-pole tracks rate higher.)")
    for Wt, rs in sorted(by_W.items()):
        b = max(rs, key=lambda r: r["rate"]["50"]["F"])
        print(f"    W = {Wt:.1f} m: best {b['rate']['50']['F']/1e3:.0f} kN/m ({b['rate']['50']['F']/Wt/1e3:.0f} kN/m per meter of width): {geom_label(b['meta'])}")
    return recs


# =============================================================================
# PART A: THE PLATE B-H CURVE AND THE SATURATION SENSITIVITY
# =============================================================================

def part_A(recs):
    hr("PART A: SOLVER VERIFICATION (bearing_fem_check.txt) AND THE PLATE B-H CURVE")
    print("  The known-answer checks, the scikit-fem cross-check and the grid and box convergence are in bearing_fem_check.txt (Part A.1 to A.3).")
    print("  Box size: the 2026-10-11 runs used a 3.0 x 1.6 m box whose A = 0 walls sat 0.3 m under the yoke and 0.6 m over the plate; that")
    print("  read the H7 linear force 20 percent low (321 against 394 to 403 kN/m in boxes of 3.8 x 2.0 m and larger).  Every solve in this")
    print("  report uses a box of at least 3.8 x 2.0 m, so the ratings here supersede yesterday's.")
    Bs, mu0r = lc.FERRO_B_SAT, lc.FERRO_MU_R
    H_knee = Bs / (MU0 * (mu0r - 1))
    print(f"  Plate B-H curve: B(H) = mu0 H + B_sat tanh(mu0 (mu_r0 - 1) H / B_sat) with mu_r0 = {mu0r} and B_sat = {Bs} T (levitation_config.py); the")
    print(f"  initial permeability is {mu0r}, the knee (tanh argument 1) is at H = {H_knee:,.0f} A/m where B = {MU0*H_knee + Bs*math.tanh(1):.2f} T, and B is within 5 percent of")
    print(f"  B_sat at H = {2*H_knee:,.0f} A/m ({MU0*2*H_knee + Bs*math.tanh(2):.2f} T).  Above that the plate adds only mu0 H.  [UNVERIFIED] the curve shape; the 1.2 T is the chapter's")
    print(f"  number for sintered iron in resin, a composite whose saturation scales with its iron fraction.  Yoke and poles: the same form with")
    print(f"  mu_r0 1,000 and B_sat 1.8 T, standing for a non-oriented silicon steel lamination (EN 10106 M250-35A or AISI M-19 class, 0.35 mm,")
    print(f"  about 3 percent silicon): handbook values put B at 1.5 T near 1,000 A/m and 1.7 to 1.8 T at 5,000 to 10,000 A/m, with a saturation")
    print(f"  polarization near 2.0 T [UNVERIFIED: no datasheet on this machine; cite the EN 10106 datasheet in the chapter].")
    best = choose_best(recs, W_max=1.0)
    gm = dict(best["gm"])
    print(f"  Sensitivity of the rated pull to the plate saturation, geometry {geom_label(best['meta'])}:")
    print(f"  {'B_sat T':>8} {'rated kN/m':>10} {'at NI kA':>8} {'at 70 %':>8} {'curve max':>9} {'B_track at rating':>17}")
    sens = {}
    for bs in (1.0, 1.2, 1.5):
        mu_r = mu0r
        curve = cf.track_curve(LEVELS, h=H_SWEEP, plate_bsat=bs, **gm)
        rt = rate_curve(curve)
        sens[bs] = rt
        print(f"  {bs:8.1f} {rt['50']['F']/1e3:10.1f} {rt['50']['NI']/1e3:8.1f} {rt['70']['F']/1e3:8.1f} {rt['F_max']/1e3:9.1f} {rt['50']['B_track']:17.3f}")
    r10, r12, r15 = sens[1.0], sens[1.2], sens[1.5]
    print(f"  Rated pull against 1.2 T: {r10['50']['F']/r12['50']['F']*100:.0f} percent at 1.0 T, {r15['50']['F']/r12['50']['F']*100:.0f} percent at 1.5 T; curve maximum "
          f"{r10['F_max']/r12['F_max']*100:.0f} and {r15['F_max']/r12['F_max']*100:.0f} percent.  The 50 percent rule sits near the knee of the")
    print(f"  curve, where the yoke's own saturation and the slot leakage already count, so the rating moves less than the plate saturation; the")
    print(f"  curve maximum, which the plate sets, moves more.")
    return sens


def choose_best(recs, W_max=None, W_min=None, plate_max=None, F_need=None):
    """The geometry with the largest rated pull among those within the width and plate limits; with F_need, the narrowest
    (then thinnest plate) that reaches it."""
    cands = [r for r in recs if (W_max is None or r["meta"]["W"] <= W_max + 1e-6) and (W_min is None or r["meta"]["W"] >= W_min - 1e-6)
             and (plate_max is None or r["meta"]["plate_t"] <= plate_max + 1e-6)]
    if F_need is not None:
        ok = [r for r in cands if r["rate"]["50"]["F"] >= F_need]
        if ok:
            return min(ok, key=lambda r: (r["meta"]["W"], r["meta"]["plate_t"], -r["rate"]["50"]["F"]))
    return max(cands, key=lambda r: r["rate"]["50"]["F"])


# =============================================================================
# TAPE DATA
# =============================================================================

def read_theva():
    data = {}
    with open(THEVA_CSV, encoding="utf-8-sig") as f:
        rd = csv.reader(f); next(rd)
        for row in rd:
            if len(row) < 4 or not row[1]:
                continue
            data.setdefault(int(round(float(row[2]))), []).append((float(row[1]), float(row[3])))
    return {a: (np.array([q[0] for q in sorted(v)]), np.array([q[1] for q in sorted(v)])) for a, v in data.items()}


def ic_tape(theva, B, w):
    """Critical current of a tape of width w (m) at field B, the lower of the 0 and 90 degree curves (A)."""
    return min(float(np.interp(B, *theva[0])), float(np.interp(B, *theva[90]))) * w * 100.0


def eddy_cycle(B, t_edge, t_hold, d_p=lc.FERRO_PLATE_THICKNESS, mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY):
    """Eddy energy per unit plate area for one rise-hold-fall cycle, J/m^2 (1-D magnetic diffusion, as Part I of the spread-pull run)."""
    nz = 100; dz = d_p / nz
    D = 1.0 / (MU0 * mu_r / rho_p); dt = 0.4 * dz * dz / D
    n = int((2 * t_edge + t_hold + 15e-3) / dt) + 1
    Bz = np.zeros(nz + 1); E = 0.0
    for j in range(n):
        t = j * dt
        Bs = B * t / t_edge if t < t_edge else (B if t < t_edge + t_hold else (B * (1 - (t - t_edge - t_hold) / t_edge) if t < 2 * t_edge + t_hold else 0.0))
        Bz[0] = Bs; Bz[-1] = Bz[-2]
        lap = (Bz[2:] - 2 * Bz[1:-1] + Bz[:-2]) / dz**2
        Jc = (Bz[1:] - Bz[:-1]) / dz / (MU0 * mu_r)
        E += rho_p * np.sum(Jc * Jc) * dz * dt
        Bz[1:-1] += D * lap * dt
    return E


def curve_interp(rec, NI, key):
    rows = [r for r in rec["curve"] if r.get("nonlinear")]
    return float(np.interp(NI, [r["NI"] for r in rows], [r[key] for r in rows]))


def NI_for_force(rec, F):
    rows = [r for r in rec["curve"] if r.get("nonlinear")]
    Fs = [r["F"] for r in rows]; NIs = [r["NI"] for r in rows]
    if F > max(Fs):
        return float("nan")
    return float(np.interp(F, Fs, NIs))


# =============================================================================
# PART C: THE CASING SUPPORT TRACKS
# =============================================================================

def part_C(recs, theva):
    hr("PART C: THE CASING SUPPORT TRACKS (each carries 95.0 kN/m alone, 12 mm tape in persistent mode at half Ic)")
    F_need = 1.5 * W
    print(f"  Target: one track alone carries the casing weight {W/1e3:.1f} kN/m with a margin of at least 1.5, so {F_need/1e3:.0f} kN/m, at no more than half of the")
    print(f"  tape's critical current at the winding field; the pair then has a margin of 3.  Width 1 m first; the plate as thick as it needs to be.")
    cand = [r for r in recs if r["meta"]["W"] <= 1.0 + 1e-6]
    print(f"  {'geometry':>96} {'rated':>6} {'NI at 143':>9} {'B_wind':>6} {'Ic 12 mm':>8} {'turns/pole':>10} {'tape m/m':>8}")
    rows = []
    for r in sorted(cand, key=lambda r: -r["rate"]["50"]["F"]):
        NI = NI_for_force(r, F_need)
        if NI != NI:
            print(f"  {geom_label(r['meta']):>96} {r['rate']['50']['F']/1e3:6.0f} {'not reached':>9}")
            continue
        Bw = math.hypot(curve_interp(r, NI, "Bx_wm"), curve_interp(r, NI, "By_wm"))
        Ic = ic_tape(theva, Bw, 0.012); I_op = 0.5 * Ic
        turns = NI / I_op; tape = 2 * NI * r["meta"]["n_groups"] / I_op
        rows.append((r, NI, Bw, Ic, turns, tape))
        print(f"  {geom_label(r['meta']):>96} {r['rate']['50']['F']/1e3:6.0f} {NI/1e3:9.1f} {Bw:6.2f} {Ic:8.0f} {turns:10.0f} {tape:8,.0f}")
    # choose: the least tape among the geometries whose rating covers 143 kN/m, else the least tape overall; then size the
    # winding at the rated pull (50 percent efficiency), which gives more margin than the 1.5 minimum for a little more tape
    ok = [x for x in rows if x[0]["rate"]["50"]["F"] >= F_need]
    pick = min(ok or rows, key=lambda x: x[5])
    r, NI_min, Bw_min, Ic_min, turns_min, tape_min = pick
    m = r["meta"]
    F_rate = r["rate"]["50"]["F"]
    NI = r["rate"]["50"]["NI"]
    Bw = math.hypot(curve_interp(r, NI, "Bx_wm"), curve_interp(r, NI, "By_wm"))
    Ic = ic_tape(theva, Bw, 0.012); I_op = 0.5 * Ic
    turns = NI / I_op; tape = 2 * NI * m["n_groups"] / I_op
    F_ic = curve_interp(r, 2 * NI, "F") if 2 * NI <= max(LEVELS) else r["rate"]["F_max"]
    F_need_margin = F_need
    F_need = F_rate                                   # the design point: the rated pull at half Ic
    n_sub = 10                                   # coils per pole, with staggered ends [UNVERIFIED choice]
    print()
    print(f"  Chosen: {geom_label(m)}.")
    print(f"    Minimum winding for the 1.5 margin: {F_need_margin/1e3:.0f} kN/m at {NI_min/1e3:.1f} kA-turns per pole, {turns_min:.0f} turns per pole, {tape_min:,.0f} m of tape per m.")
    print(f"    Design point, the rated pull at half Ic: {F_need/1e3:.0f} kN/m at {NI/1e3:.1f} kA-turns per pole ({m['n_groups']} poles, {NI*m['n_groups']/1e3:.0f} kA-turns per track);")
    print(f"    winding field there {Bw:.2f} T mean, 12 mm tape Ic {Ic:.0f} A (THEVA 70 K), run at {0.5*Ic:.0f} A: {turns:.0f} turns per pole, {turns/n_sub:.0f} per coil with {n_sub} coils per pole.")
    print(f"    Normal operation, the pair sharing the weight: {W/2e3:.1f} kN/m per track at {NI_for_force(r, W/2)/1e3:.1f} kA-turns per pole, {NI_for_force(r, W/2)/turns:.0f} A in the tape.")
    print(f"    Rated pull (50 percent efficiency) {F_rate/1e3:.0f} kN/m per track, {2*F_rate/1e3:.0f} kN/m for the pair; curve maximum {r['rate']['F_max']/1e3:.0f} kN/m at {r['rate']['NI_max']/1e3:.0f} kA.")
    print(f"    At the tape's full Ic ({Ic:.0f} A, {2*NI/1e3:.0f} kA per pole) the pull would be {F_ic/1e3:.0f} kN/m, but persistent operation stays at half Ic.")
    print(f"    Margin on one track at half Ic: {F_need/W:.2f} times the weight, {(F_need - W)/1e3:.1f} kN/m = {(F_need - W)/pz.G_NET:,.0f} kg/m of extra local load; on the pair at half Ic: "
          f"{2*F_need/W:.2f} times, {(2*F_need - W)/1e3:.0f} kN/m = {(2*F_need - W)/pz.G_NET:,.0f} kg/m.")
    print(f"    Tape: {tape:,.0f} m per meter of ring per track, {tape*TAPE_T*0.012*TAPE_RHO:.1f} kg/m; two tracks around the ring {2*tape*L_RING/1e9:,.1f} million km, {2*tape*TAPE_T*0.012*TAPE_RHO*L_RING/1e9:.2f} Mt.")
    print(f"    Plate {m['plate_t']*1e3:.0f} mm thick, {m['plate_w']:.1f} m wide, {m['iron_cable_m2']*RHO_PLATE:,.0f} kg/m; yoke and poles {m['iron_casing_m2']*RHO_LAM:,.0f} kg/m of laminations per track.")
    print(f"    This replaces Equations 8.1 to 8.4: 218 kN/m (two tracks at 0.523 T), 871 kN/m (at Ic), 573 kPa (1.2 T) and 776 kN/m (margin) came from the")
    print(f"    one-way sheet formula, which the solves do not support (Part B of compensation_tracks_2026-10-11.md: a flat array gives 0.2 kN/m).")
    # ripple along the direction of travel
    print()
    print(f"  Field ripple along the ring from the coil ends.  Each pole carries {n_sub} coils of one tenth of its ampere-turns; a coil's end turn leaves a")
    print(f"  dip of 1/{n_sub} of the pole's ampere-turns over about one slot width ({m['slot_w']*1e3:.0f} mm) [UNVERIFIED end-turn length], so the plate field dips by")
    print(f"  about {100/n_sub:.0f} percent there and the pull by {100*(1-(1-1/n_sub)**2):.0f} percent, for {m['slot_w']/V_REL*1e6:.0f} us at {V_REL:,.0f} m/s.  The ends are staggered so no two")
    print(f"  coils end within a coil length divided by the number of coils of each other.  Eddy heating of the plate per dip crossing, from the 1-D diffusion,")
    print(f"  and the heating per meter of ring for several coil lengths (dips per meter = poles x coils per pole / coil length, two tracks):")
    B_op = curve_interp(r, NI_for_force(r, W / 2), "B_track")
    dB = B_op / n_sub
    E_dip = eddy_cycle(dB, 0.02 / V_REL, m["slot_w"] / V_REL, d_p=m["plate_t"])
    print(f"    operating plate field {B_op:.3f} T, dip {dB:.3f} T, eddy energy per dip {E_dip*1e3:.3f} mJ per m^2 of plate, over {m['pole_w']:.2f} m of pole width")
    print(f"    {'coil length m':>13} {'dips per m':>10} {'heating W/m':>11} {'flux pumps per track':>20} {'energy per coil kJ':>18} {'coil inductance H':>17}")
    W_op = curve_interp(r, NI_for_force(r, W / 2), "W_mag")
    for Lc in (200.0, 500.0, 1000.0, 2000.0, 5000.0):
        dips = m["n_groups"] * n_sub / Lc
        P = 2 * E_dip * m["pole_w"] * V_REL * dips
        E_coil = W_op * Lc / n_sub / m["n_groups"] * m["n_groups"]       # per coil: the pole's share over its length / coils per pole
        E_coil = W_op * Lc / (n_sub * m["n_groups"])
        I_tape = NI_for_force(r, W / 2) / turns
        print(f"    {Lc:13.0f} {dips:10.3f} {P:11.3f} {m['n_groups']*n_sub*L_RING/Lc/1e6:17.2f} M {E_coil/1e3:18.1f} {2*E_coil/I_tape**2:17.1f}")
    print(f"  Against the 73.6 W/m of the old limiter plates at the nominal gap, the ripple heating is small for any coil length over 200 m; the")
    print(f"  coil length for the support tracks is set by the flux pump count and the energy a coil holds (a quench dumps it), not by ripple.")
    print(f"  Recommendation [UNVERIFIED trade]: 1,000 m coils, {m['n_groups']*n_sub} per 1,000 m of track, {m['n_groups']*n_sub*L_RING/1000/1e6:.2f} million flux pumps per track.")
    return dict(rec=r, NI=NI, Bw=Bw, Ic=Ic, turns=turns, tape=tape, F_rate=F_rate, F_need=F_need, n_sub=n_sub, L_coil=1000.0, B_op=B_op,
                NI_op=NI_for_force(r, W / 2), F_one_cap=F_need)


# =============================================================================
# PART D: THE COMPENSATION TRACKS
# =============================================================================

def width_budget(W_s, W_load, W_lift):
    under = 2 * (W_s + 0.4) + 2 * (W_load + 0.4)
    top = 2 * (W_lift + 0.4)
    return under, top


def part_D(recs, theva, C, tabs):
    hr("PART D: THE COMPENSATION TRACKS (either track of a pair does the pair's job alone; 3 mm tape)")
    print(f"  Needs per track alone: lift {F_LIFT_ALONE/1e3:.0f} kN/m at 30 km/s; load {F_LOAD_ALONE/1e3:.0f} kN/m at 1 km/s (sled_passage_spread.txt rule 4, 50 m sections) and the parked")
    print(f"  5,000 t on 100 m {F_PARKED/1e3:.0f} kN/m; both must also hold the casing ({W/1e3:.0f} kN/m) with both support tracks lost.  Width budget: cable face {CABLE_FACE:.2f} m;")
    print(f"  underside two support plates and two load plates, top two lift plates, each plate 0.4 m wider than its track.")
    W_s = C["rec"]["meta"]["W"]
    lift = choose_best(recs, F_need=F_LIFT_ALONE)
    load = choose_best(recs, F_need=max(F_LOAD_ALONE, F_PARKED))
    if lift["rate"]["50"]["F"] < F_LIFT_ALONE:
        print(f"  No geometry in the sweep reaches {F_LIFT_ALONE/1e3:.0f} kN/m per track; the largest is {lift['rate']['50']['F']/1e3:.0f} kN/m ({geom_label(lift['meta'])}).")
    under, top = width_budget(W_s, load["meta"]["W"], lift["meta"]["W"])
    print(f"  Lift: {geom_label(lift['meta'])}: rated {lift['rate']['50']['F']/1e3:.0f} kN/m per track at {lift['rate']['50']['NI']/1e3:.0f} kA per pole.")
    print(f"  Load: {geom_label(load['meta'])}: rated {load['rate']['50']['F']/1e3:.0f} kN/m per track at {load['rate']['50']['NI']/1e3:.0f} kA per pole.")
    print(f"  Width used: underside {under:.2f} m of {CABLE_FACE:.2f} ({'fits' if under <= CABLE_FACE else 'DOES NOT FIT'}), top {top:.2f} m ({'fits' if top <= CABLE_FACE else 'DOES NOT FIT'}).")
    out = {}
    for tag, rec, F_need in (("lift", lift, F_LIFT_ALONE), ("load", load, max(F_LOAD_ALONE, F_PARKED))):
        m = rec["meta"]
        NI = NI_for_force(rec, F_need) if rec["rate"]["F_max"] >= F_need else rec["rate"]["50"]["NI"]
        F_one = min(F_need, rec["rate"]["50"]["F"]) if NI == NI else rec["rate"]["50"]["F"]
        Bw = math.hypot(curve_interp(rec, NI, "Bx_wm"), curve_interp(rec, NI, "By_wm"))
        Bx_x, By_x = curve_interp(rec, NI, "Bx_w"), curve_interp(rec, NI, "By_w")
        Ic = ic_tape(theva, Bw, 0.003); I_op = 0.5 * Ic
        turns = NI / I_op; tape = 2 * NI * m["n_groups"] / I_op
        W_mag = curve_interp(rec, NI, "W_mag")
        out[tag] = dict(rec=rec, NI=NI, F_one=F_one, Bw=Bw, Ic=Ic, I_op=I_op, turns=turns, tape=tape, W_mag=W_mag, By_m=curve_interp(rec, NI, "By_wm"),
                        Bx_m=curve_interp(rec, NI, "Bx_wm"), W_track=m["W"], plate_t=m["plate_t"], n_groups=m["n_groups"])
        print(f"  {tag} track at the one-track pull {F_one/1e3:.0f} kN/m: {NI/1e3:.1f} kA per pole, {NI*m['n_groups']/1e3:.0f} kA per track; winding {Bw:.2f} T mean (across {out[tag]['Bx_m']:.2f}, normal "
              f"{out[tag]['By_m']:.2f}; corner maxima {Bx_x:.2f}/{By_x:.2f}); 3 mm tape Ic {Ic:.0f} A, run at {I_op:.0f} A: {turns:.0f} turns per pole, {tape:,.0f} m of tape per m, "
              f"{tape*TAPE_T*0.003*TAPE_RHO:.1f} kg/m; stored {W_mag/1e3:.0f} kJ/m ({W_mag/(F_one*GAP0):.1f} x F gap); plate {m['plate_t']*1e3:.0f} mm, {m['iron_cable_m2']*RHO_PLATE:,.0f} kg/m; "
              f"yoke and poles {m['iron_casing_m2']*RHO_LAM:,.0f} kg/m.")
    # ---- coil length and switching: point model
    print()
    print("  Coil length.  Point model at 1, 3, 5, 10, 20 and 30 km/s (pair, caps at the pair rating, 2 ms hardware rise, no ramp): casing movement in mm")
    print("  for coils of 10, 25, 50, 100 and 200 m, staggered (20 coil families, a new coil every L/20) and switched as a block (sections of L).")
    import compensation_tracks as ct
    cap_lift2 = 2 * out["lift"]["F_one"]; cap_load2 = 2 * out["load"]["F_one"]
    w_hd, w_hu = cap_lift2 / P12, cap_load2 / P12
    speeds = (1000.0, 3000.0, 5000.0, 10000.0, 20000.0, 30000.0)
    lengths = (10.0, 25.0, 50.0, 100.0, 200.0)
    print(f"  {'u km/s':>7} " + " ".join(f"{'stag ' + str(int(L)):>9}" for L in lengths) + " " + " ".join(f"{'block ' + str(int(L)):>9}" for L in lengths))
    t0 = time.time()
    exc = {}; F_by_L = {}
    for u in speeds:
        T = max(0.05, L_SC / u)
        lanes = [pz.lane(L_sc=L_SC, tau=TAU_HW, L_section=L / 20, gap_dep=False, steady_loop_hz=1.0, T_spread=T, B_sat=1.2, n_coils=20, L_coil=L,
                         hd_width=w_hd, hu_width=w_hu, label=f"s{L}") for L in lengths]
        lanes += [pz.lane(L_sc=L_SC, tau=TAU_HW, L_section=0.0, gap_dep=False, steady_loop_hz=1.0, T_spread=T, B_sat=1.2, n_coils=1, L_coil=L,
                          hd_width=w_hd, hu_width=w_hu, label=f"b{L}") for L in lengths]
        res, _ = pz.point_pulsed(lanes, u, tabs)
        exc[u] = [r["max_exc_all"] for r in res]
        F_by_L[u] = {L: max(r["pk_hd"], r["pk_hu"]) for L, r in zip(lengths, res[:len(lengths)])}
        print(f"  {u/1e3:7.0f} " + " ".join(f"{v*1e3:9.2f}" for v in exc[u]))
    print(f"  ({time.time()-t0:.0f} s)")
    # supplies, energy and VA per coil
    print()
    print(f"  Supplies, stored energy and VA per coil (the energy is the track's stored energy at the one-track pull times the coil's share of length and")
    print(f"  ampere-turns; VA = 2 E / rise for a linear rise; the hardware rise is {TAU_HW*1e3:.0f} ms):")
    print(f"  {'track':>5} {'coil m':>6} {'layout':>7} {'supplies per track':>18} {'E per coil MJ':>13} {'VA at 2 ms GVA':>14} {'VA at 10 ms':>11} {'VA at 50 ms':>11}")
    for tag in ("lift", "load"):
        d = out[tag]
        for L in lengths:
            for layout, n_fam in (("stagger", 20), ("block", 1)):
                E = d["W_mag"] * L / n_fam
                n_sup = n_fam * L_RING / L
                print(f"  {tag:>5} {L:6.0f} {layout:>7} {n_sup/1e6:15.2f} M {E/1e6:13.2f} {2*E/TAU_HW/1e9:14.2f} {2*E/10e-3/1e9:11.2f} {2*E/50e-3/1e9:11.2f}")
    # single coil failure (beam) for the lengths
    print()
    print("  Single coil failure, beam model (10 m truss), lift at 30 km/s, one staggered family's coil dead at the window middle:")
    print(f"  {'coil m':>6} {'normal exc mm':>13} {'one coil dead exc mm':>20} {'floor gap':>9} {'ceiling gap':>11}")
    t0 = time.time()
    for L in (25.0, 50.0, 100.0, 200.0):
        res, _ = pz.beam_pulsed(30000.0, L_SC, [sp.truss(10.0, 0.1)], [dict(label="normal"), dict(label="coil dead", fail_coil=(0, 5000.0))], tabs,
                                tau=TAU_HW, L_win=10000.0, T_spread=0.05, n_coils=20, L_coil=L, track_width=w_hd, prot=ct.prot(0.010, 1e-3, cap_up=cap_load2, cap_down=cap_lift2))
        print(f"  {L:6.0f} {res[0]['max_exc']*1e3:13.1f} {res[1]['max_exc']*1e3:20.1f} {res[1]['min_gap_b']*1e3:9.1f} {res[1]['min_gap_t']*1e3:11.1f}")
    print(f"  ({time.time()-t0:.0f} s)")
    # choice of length
    L_lift, L_load = 50.0, 25.0
    print(f"  Choice: lift coils {L_lift:.0f} m staggered (movement within 11 mm from 3 km/s up, {20*L_RING/L_lift/1e6:.1f} million supplies per track), load coils {L_load:.0f} m")
    print(f"  staggered (the low-speed smear is the load pair's problem; {20*L_RING/L_load/1e6:.1f} million supplies per track).  [QUESTION FOR AUTHOR] the table shows the trade.")
    # rise time
    print()
    print(f"  Rise time for the planned pull: the command is ramped over 0 (hardware 2 ms only), 10, 20 and 50 ms; casing movement in mm (pair):")
    print(f"  {'u km/s':>7} {'2 ms':>7} {'10 ms':>7} {'20 ms':>7} {'50 ms':>7}")
    ramps = (0.0, 10e-3, 20e-3, 50e-3)
    t0 = time.time()
    rise_exc = {}
    for u in (1000.0, 3000.0, 10000.0, 30000.0):
        T = max(0.05, L_SC / u)
        L = L_load if u < 7000 else L_lift
        lanes = [pz.lane(L_sc=L_SC, tau=TAU_HW, L_section=L / 20, gap_dep=False, steady_loop_hz=1.0, T_spread=T, B_sat=1.2, n_coils=20, L_coil=L,
                         hd_width=w_hd, hu_width=w_hu, cmd_ramp=rp, label=f"r{rp}") for rp in ramps]
        res, _ = pz.point_pulsed(lanes, u, tabs)
        rise_exc[u] = [r["max_exc_all"] for r in res]
        print(f"  {u/1e3:7.0f} " + " ".join(f"{v*1e3:7.2f}" for v in rise_exc[u]))
    print(f"  ({time.time()-t0:.0f} s)")
    ramp_ok = max((rp for rp in ramps if all(rise_exc[u][ramps.index(rp)] <= 0.0105 for u in rise_exc)), default=0.0)
    print(f"  The longest ramp that keeps the planned pull inside the 10 mm strict band at every speed: {ramp_ok*1e3:.0f} ms (plus the 2 ms hardware rise); the")
    print(f"  tolerances then keep their room up to the 20 mm band.  The protection keeps the 2 ms hardware rise.")
    for tag in ("lift", "load"):
        d = out[tag]; L = L_lift if tag == "lift" else L_load
        E = d["W_mag"] * L / 20
        d.update(L_coil=L, E_coil=E, VA_hw=2 * E / TAU_HW, VA_ramp=2 * E / max(ramp_ok + TAU_HW, TAU_HW), ramp=ramp_ok)
        print(f"  {tag}: coil {L:.0f} m, {E/1e6:.2f} MJ per coil, supply {2*E/TAU_HW/1e9:.2f} GVA for the 2 ms protection rise, {2*E/(ramp_ok+TAU_HW)/1e9:.3f} GVA for the planned pull "
              f"({d['I_op']:.0f} A in one 3 mm tape: {2*E/TAU_HW/d['I_op']/1e6:.1f} MV at 2 ms; with 100 tapes in parallel {2*E/TAU_HW/d['I_op']/100/1e3:.0f} kV at {d['I_op']*100/1e3:.1f} kA).")
    # hold with both support tracks lost and the parked spacecraft
    print()
    print(f"  Holding the casing with both support tracks lost: {W/1e3:.0f} kN/m from the load pair, {W/2e3:.1f} kN/m per track at {NI_for_force(load, W/2)/1e3:.1f} kA per pole "
          f"({NI_for_force(load, W/2)/out['load']['NI']*100:.0f} percent of the one-track ampere-turns); one load track alone at {NI_for_force(load, W)/1e3:.1f} kA per pole.")
    print(f"  Parked 5,000 t on 100 m: {F_PARKED/1e3:.0f} kN/m from the load pair, {F_PARKED/2e3:.0f} per track ({NI_for_force(load, F_PARKED/2)/1e3:.1f} kA per pole); one track alone "
          f"{'reaches it' if NI_for_force(load, F_PARKED) == NI_for_force(load, F_PARKED) else 'does not reach it'} at {NI_for_force(load, F_PARKED)/1e3:.1f} kA per pole.")
    out["caps"] = dict(lift1=out["lift"]["F_one"], load1=out["load"]["F_one"], lift2=cap_lift2, load2=cap_load2)
    out["exc"] = exc; out["rise_exc"] = rise_exc
    out["F_u"] = {u: F_by_L[u][L_load if u < 7000 else L_lift] for u in speeds}
    return out


# =============================================================================
# PART E: TAPE LOSSES WITH THE BOOK'S FORMULAS
# =============================================================================

def part_E(D, theva):
    hr("PART E: TAPE LOSSES WITH THE BOOK'S OWN FORMULAS (lim_physics.py) AGAINST THE STRIP ESTIMATE")
    d = D["lift"]
    m = d["rec"]["meta"]
    I, Ic, N = d["I_op"], d["Ic"], d["turns"] / 20            # turns of one coil (20 coils per pole group of ... see below)
    N_coil = d["turns"] / (20 / m["n_groups"])                 # turns per coil when the 20 coils of a track are split over the poles
    w_tape = 0.003
    alpha = lim_config.ALPHA_TAPE
    w_coil = m["slot_w"]
    print(f"  Three estimates of the loss in one meter of tape for one pulse (0 to I to 0), lift track at its one-track pull:")
    print(f"  tape 3 mm, Ic {Ic:.0f} A at the winding field {d['Bw']:.2f} T, operating current {I:.0f} A, {N_coil:.0f} turns per coil in a {w_coil*1e3:.0f} mm slot.")
    q_strip = 0.5 * d["By_m"] * w_tape * Ic
    q_norris = lim_physics.q_hysteresis_norris(I, Ic)
    q_book = lim_physics.q_hts_loss_factor(I, N_coil, w_coil, w_tape, alpha)
    b_coil = MU0 * N_coil * I / w_coil
    print(f"    strip estimate (compensation_tracks.py): 1/2 B_perp w Ic with B_perp = {d['By_m']:.3f} T: {q_strip*1e3:.2f} mJ per m of tape per pulse")
    print(f"    q_hysteresis_norris(I, Ic) (Norris strip, transport current, per full cycle +-I): {q_norris*1e3:.3f} mJ per m per cycle, {q_norris/2*1e3:.3f} for a 0 to I to 0 pulse (half a cycle)")
    print(f"    q_hts_loss_factor(I, N, w_coil, w_tape, alpha): B_coil = mu0 N I / w_coil = {b_coil:.3f} T, times I w_tape sin(alpha) with alpha = {lim_config.ALPHA_PENETRATION_DEG:.0f} deg: "
          f"{q_book*1e3:.2f} mJ per m per cycle")
    print(f"  calc_hysteresis_power multiplies the per-cycle loss by the tape length and the supply frequency (lim_simulation.py line 312, switch NORRIS_HYSTERESIS")
    print(f"  in lim_config.py, {lim_config.NORRIS_HYSTERESIS} by default, so the LIM runs use q_hts_loss_factor).")
    print(f"  What they are: the Norris formula is the hysteresis loss of an isolated strip carrying its own transport current, which goes as Ic^2 times a")
    print(f"  function of I/Ic and does not depend on any external field; it is the loss of the self-field.  q_hts_loss_factor is a magnetization loss")
    print(f"  with the coil's own field across its width (mu0 N I / w_coil) standing in for the applied field and sin(alpha) for the fraction of the")
    print(f"  tape face that field crosses; it is of the form B x I x w.  The strip estimate is the magnetization loss of a fully penetrated strip in an")
    print(f"  applied field normal to its face, 1/2 B_perp w Ic for one 0 to B to 0 sweep: it uses Ic, not I, because the screening currents in a")
    print(f"  penetrated strip run at Ic whatever the transport current, and it uses the field normal to the tape, which the 2-D solve gives directly.")
    print(f"  The ratio strip / book = {q_strip/q_book:.1f}: the solve's normal field at the winding ({d['By_m']:.3f} T) against the coil-width field times sin(20 deg) ({b_coil*math.sin(alpha):.3f} T),")
    print(f"  and Ic against I.  In a pulsed coil the transport and magnetization losses add; the transport part (Norris) is {q_norris/2/q_strip*100:.2f} percent of the strip estimate.")
    # per pulse, per launch
    t_p, v_p, P_p = pz.lsm_profile()
    dt_p = np.diff(t_p)
    cop = lim_physics.calc_cryo_cop(70.0, 400.0, lim_config.CRYO_EFF)
    KE = 0.5 * (M_SC0 + M_SLED_M * L_SLED) * 30000**2
    print()
    print(f"  Per pulse and per launch, 3 mm tape, the lift and load tracks each in their speed range, each meter of ring pulsed once per lap; the")
    print(f"  loss per pulse scales with the normal winding field, which the saturating curve gives at each speed's ampere-turns (pair mode: each")
    print(f"  track at half the pull).  COP {cop:.4f} (lim_physics.calc_cryo_cop(70, 400, {lim_config.CRYO_EFF})).  Sled kinetic energy {KE/1e12:.0f} TJ.")
    rows_A = D["F_u"]
    print(f"  {'model':>22} {'J per pulse per m at 30 km/s':>28} {'heat per launch GJ':>18} {'electrical GJ':>13} {'% of KE':>7} {'mean GW':>8}")
    results = {}
    for name, per_m in (("strip, 1/2 B_perp w Ic", None), ("book q_hts_loss_factor", "book"), ("Norris (half cycle)", "norris")):
        E_launch = 0.0; q30 = 0.0
        for i, dtt in enumerate(dt_p):
            vv = v_p[i + 1]
            tag = "lift" if vv > 7000 else "load"
            dd = D[tag]; rec = dd["rec"]; mm = rec["meta"]
            Fu = float(np.interp(vv, sorted(rows_A), [rows_A[u] for u in sorted(rows_A)]))
            NI_u = NI_for_force(rec, Fu / 2)
            if NI_u != NI_u:
                NI_u = rec["rate"]["F_max"] and rec["rate"]["NI_max"]
            By = curve_interp(rec, NI_u, "By_wm")
            I_u = NI_u / dd["turns"]
            Nc = dd["turns"] / (20 / mm["n_groups"])
            if per_m is None:
                q = 0.5 * By * w_tape * dd["Ic"]
            elif per_m == "book":
                q = lim_physics.q_hts_loss_factor(I_u, Nc, mm["slot_w"], w_tape, alpha)
            else:
                q = 0.5 * lim_physics.q_hysteresis_norris(I_u, dd["Ic"])
            q_m = q * dd["tape"] * 2                      # two tracks of the pair per meter of ring
            if abs(vv - 30000) < 300 and q30 == 0.0:
                q30 = q_m
            E_launch += q_m * vv * dtt
        results[name] = dict(E=E_launch, E_el=E_launch / cop, q30=q30)
        print(f"  {name:>22} {q30:28.1f} {E_launch/1e9:18.2f} {E_launch/cop/1e9:13.1f} {E_launch/cop/KE*100:7.3f} {E_launch/cop/t_p[-1]/1e9:8.2f}")
    print(f"  THEVA angle convention: Chapter 6, Table 6.2 reproduces the 0 degree column of this data file as the field perpendicular to the tape; the")
    print(f"  file's 90 degree column is then the field along the tape face (the Advanced Pinning tape's in-plane peak).  Confirmed from the book and the")
    print(f"  file, not from the datasheet, which is not on this machine.")
    return results


def ct_rows_A():
    """The pull at each speed from the 2026-10-11 Part A run, read from its log (compensation_tracks.txt), kN/m -> N/m."""
    path = os.path.join(RUNS_DIR, "compensation_tracks.txt")
    rows = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    i0 = next(i for i, l in enumerate(lines) if "u km/s   T ms  face pull kN/m" in l)
    for l in lines[i0 + 1:i0 + 6]:
        parts = l.split()
        rows[float(parts[0]) * 1e3] = float(parts[3]) * 1e3
    return rows


# =============================================================================
# PART F: THE DYNAMICS WITH THE CHOSEN DESIGN
# =============================================================================

def part_F(C, D, tabs):
    hr("PART F: THE DYNAMICS WITH THE CHOSEN DESIGN")
    import compensation_tracks as ct
    caps = D["caps"]
    lift1, load1, lift2, load2 = caps["lift1"], caps["load1"], caps["lift2"], caps["load2"]
    L_lift, L_load = D["lift"]["L_coil"], D["load"]["L_coil"]
    ramp = D["lift"]["ramp"]
    S1 = C["F_one_cap"]
    print(f"  Caps: lift {lift1/1e3:.0f} kN/m per track ({lift2/1e3:.0f} pair), load {load1/1e3:.0f} ({load2/1e3:.0f}); support one track alone {S1/1e3:.0f} kN/m at half Ic, pair {2*S1/1e3:.0f}.")
    print(f"  Coils {L_lift:.0f} m (lift) and {L_load:.0f} m (load) staggered in 20 families; planned pull ramped over {ramp*1e3:.0f} ms plus the {TAU_HW*1e3:.0f} ms hardware rise; protection")
    print(f"  at 10 mm of drift, 1 ms sensor delay, 5 Hz, 2 ms rise, per face and per side.  Point model unless stated; beam = 10 m truss.")
    w_hd, w_hu = lift2 / P12, load2 / P12
    all_gaps = []

    def note(label, model, r):
        all_gaps.append((label, model, r["min_gap_b"], r["min_gap_t"], bool(r.get("contact", False))))

    def lane(u, **kw):
        L = L_load if u < 7000 else L_lift
        d = dict(L_sc=L_SC, tau=TAU_HW, L_section=L / 20, gap_dep=False, steady_loop_hz=1.0, T_spread=max(0.05, L_SC / u), B_sat=1.2, n_coils=20, L_coil=L,
                 hd_width=w_hd, hu_width=w_hu, cmd_ramp=ramp, steady_cap=2 * S1)
        d.update(kw)
        return pz.lane(**d)

    def P(u, cap_up=None, cap_down=None, **kw):
        return ct.prot(0.010, 1e-3, cap_up=load2 if cap_up is None else cap_up, cap_down=lift2 if cap_down is None else cap_down, **kw)

    # F.1 sled passage and tolerances
    print()
    print(f"  F.1  Sled passage, pair and one track alone (the one-track cap is the one-track rating):")
    print(f"  {'u km/s':>7} {'pull kN/m':>9} {'exc pair':>8} {'exc one':>7} {'gap_b':>6} {'gap_t':>6}")
    for u in (1000.0, 3000.0, 5000.0, 10000.0, 20000.0, 30000.0):
        res, _ = pz.point_pulsed([lane(u, label="pair"), lane(u, hd_width=lift1 / P12, hu_width=load1 / P12, label="one")], u, tabs)
        F = max(res[0]["pk_hd"], res[0]["pk_hu"])
        note(f"sled passage {u/1e3:.0f} km/s, pair", "point", res[0]); note(f"sled passage {u/1e3:.0f} km/s, one track", "point", res[1])
        print(f"  {u/1e3:7.0f} {F/1e3:9.1f} {res[0]['max_exc_all']*1e3:8.2f} {res[1]['max_exc_all']*1e3:7.2f} {min(res[0]['min_gap_b'], res[1]['min_gap_b'])*1e3:6.1f} {min(res[0]['min_gap_t'], res[1]['min_gap_t'])*1e3:6.1f}")
    tol = [("timing -5 ms", dict(dt_err=-5e-3)), ("timing +5 ms", dict(dt_err=5e-3)), ("mass -10 %", dict(mass_err=-0.10)), ("mass +10 %", dict(mass_err=0.10)),
           ("speed -5 %", dict(speed_err=-0.05)), ("speed +5 %", dict(speed_err=0.05)), ("sensor -5 mm", dict(gap_err=-0.005)), ("sensor +5 mm", dict(gap_err=0.005)),
           ("rise -50 %", dict(tau_err=-0.5)), ("rise +50 %", dict(tau_err=0.5))]
    print(f"  Tolerances (pair), excursion in mm at 1, 3, 10 and 30 km/s:")
    per = {}
    for u in (1000.0, 3000.0, 10000.0, 30000.0):
        res, _ = pz.point_pulsed([lane(u, label=lab, **kw) for lab, kw in tol], u, tabs)
        per[u] = res
        for lab, r in zip([t[0] for t in tol], res):
            note(f"tolerance {lab} at {u/1e3:.0f} km/s", "point", r)
    for i, (lab, kw) in enumerate(tol):
        print(f"    {lab:>14} " + " ".join(f"{per[u][i]['max_exc_all']*1e3:7.2f}" for u in per))
    # F.2 single failures at 30 km/s
    print()
    print(f"  F.2  Single failures at 30 km/s (pair, protection on):")
    u30 = 30000.0; tc = 0.5 * L_SC / u30
    cases = [("normal", dict()),
             ("sensor error +5 mm", dict(gap_err=0.005, prot=P(u30, sensor_err=0.005))), ("sensor error -5 mm", dict(gap_err=-0.005, prot=P(u30, sensor_err=-0.005))),
             ("false trigger with 5 mm sensor error", dict(prot=ct.prot(0.0, 1e-3, cap_up=load2, cap_down=lift2, sensor_err=0.005))),
             ("spread pull fired with no sled", dict(no_sled=True)), ("spread pull fired with no sled, one track", dict(no_sled=True, hd_width=lift1 / P12, hu_width=load1 / P12)),
             ("false full lift pull 20 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 20e-3, -lift2))),
             ("false full lift pull 50 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 50e-3, -lift2))),
             ("false full load pull 50 ms, no sled", dict(no_sled=True, cancel=False, false_pull=(0.0, 50e-3, load2))),
             ("pull 10 ms late", dict(dt_err=10e-3)), ("pull 50 ms late", dict(dt_err=50e-3)), ("pull 50 ms early", dict(dt_err=-50e-3)),
             ("one lift track lost mid-pull, nothing reassigned", dict(pull_scale=[(tc, 1e9, 0.5)], prot=P(u30, cap_down=lift1))),
             ("one lift track lost at the start, other takes over after 5 ms", dict(hd_width=lift1 / P12, pull_scale=[(tc - 0.025, tc - 0.020, 0.5)], prot=P(u30, cap_down=lift1))),
             ("both support tracks lost as the zone arrives", dict(support_fail_t=0.0)), ("both support tracks lost 100 ms before", dict(support_fail_t=-0.1)),
             ("one support track out for repair (cap one track), normal pass", dict(steady_cap=S1)),
             ("one support track out, support loop lost too", dict(steady_cap=S1, steady_loop_hz=0.0)),
             ("one support track lost as the zone arrives, other takes over", dict(steady_cap=S1, support_fail_t=0.0, support_fail_frac=0.5)),
             ("one support track lost 100 ms before, other takes over", dict(steady_cap=S1, support_fail_t=-0.1, support_fail_frac=0.5))]
    lanes = []
    for lab, kw in cases:
        kw2 = dict(kw)
        if "prot" not in kw2:
            kw2["prot"] = P(u30)
        lanes.append(lane(u30, label=lab, **kw2))
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    print(f"  {'case':>62} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'engaged ms':>10}")
    for (lab, kw), r in zip(cases, res):
        note(lab, "point", r)
        print(f"  {lab:>62} {r['max_exc_all']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['t_engaged']*1e3:10.1f}")
    # at rest: support lost, and one track out with the parked load
    print(f"  At rest (no sled):")
    cases_r = [("both support tracks lost, load pair catches", dict(no_sled=True, cancel=False, support_fail_t=0.0)),
               ("both support tracks lost, one load track catches", dict(no_sled=True, cancel=False, support_fail_t=0.0, prot=P(u30, cap_up=load1))),
               ("one support track lost at rest, other takes over (1 Hz loop)", dict(no_sled=True, cancel=False, steady_cap=S1, support_fail_t=0.0, support_fail_frac=0.5)),
               ("one support track lost at rest, other takes over, no protection", dict(no_sled=True, cancel=False, steady_cap=S1, support_fail_t=0.0, support_fail_frac=0.5, prot=None))]
    lanes = []
    for lab, kw in cases_r:
        kw2 = dict(kw)
        if "prot" not in kw2:
            kw2["prot"] = P(u30)
        lanes.append(lane(u30, label=lab, **kw2))
    res, _ = pz.point_pulsed(lanes, u30, tabs, t_settle=3.0)
    for (lab, kw), r in zip(cases_r, res):
        note(lab, "point", r)
        print(f"  {lab:>62} {r['max_exc_all']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['t_engaged']*1e3:10.1f}")
    # beam
    print(f"  Beam model, 10 m truss, 30 km/s, lift coils {L_lift:.0f} m:")
    cases_b = [dict(label="normal"), dict(label="one coil dead", fail_coil=(0, 5000.0)), dict(label="spread pull with no sled", no_sled=True),
               dict(label="both support tracks lost as the zone arrives", support_fail_t=5000.0 / u30),
               dict(label="both support tracks lost at rest", no_sled=True, no_pulse=True, support_fail_t=0.0),
               dict(label="false full lift pull 20 ms, no sled", no_sled=True, no_pulse=True, false_pull=(5000.0 / u30, 20e-3, -lift2))]
    t0 = time.time()
    res, _ = pz.beam_pulsed(u30, L_SC, [sp.truss(10.0, 0.1)], cases_b, tabs, tau=TAU_HW, L_win=10000.0, T_spread=0.05, n_coils=20, L_coil=L_lift,
                            track_width=w_hd, prot=P(u30))
    for r in res:
        note(f"{r['case']}, beam", "beam", r)
        print(f"  {r['case']:>62} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['t_engaged']*1e3:10.1f}{'  CONTACT' if r['contact'] else ''}")
    res, _ = pz.beam_pulsed(u30, L_SC, [sp.truss(10.0, 0.1)], [dict(label="normal, one lift track alone"), dict(label="one coil dead, one lift track alone", fail_coil=(0, 5000.0))],
                            tabs, tau=TAU_HW, L_win=10000.0, T_spread=0.05, n_coils=20, L_coil=L_lift, track_width=lift1 / P12, prot=P(u30, cap_down=lift1))
    for r in res:
        note(f"{r['case']}, beam", "beam", r)
        print(f"  {r['case']:>62} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['t_engaged']*1e3:10.1f}{'  CONTACT' if r['contact'] else ''}")
    print(f"  ({time.time()-t0:.0f} s)")
    # F.3 roll
    print()
    W_s = C["rec"]["meta"]["W"]; W_lo = D["load"]["W_track"]; W_li = D["lift"]["W_track"]
    y_load = CABLE_FACE / 2 - (W_lo + 0.4) / 2
    y_sup = CABLE_FACE / 2 - (W_lo + 0.4) - 0.2 - (W_s + 0.4) / 2
    y_lift = CABLE_FACE / 2 - (W_li + 0.4) / 2
    side = 10.0
    r_g = side / math.sqrt(3.0)
    ct.Y_SUPPORT = y_sup; ct.Y_COMP = y_lift; ct.Y_PLATE_EDGE = CABLE_FACE / 2; ct.Y_CABLE_EDGE = CABLE_FACE / 2
    ct.I_CASING = MU * r_g**2; ct.I_CABLE = MU_C * (CABLE_FACE**2 + CABLE_FACE**2) / 12.0
    print(f"  F.3  Roll with one lift track alone (2-DOF casing and cable).  Casing assumed a {side:.0f} m square truss box (the 10 m depth of the design rule")
    print(f"  taken as both height and width), its {MU:,.0f} kg/m on the perimeter, radius of gyration {r_g:.2f} m, roll inertia {MU*r_g**2:,.0f} kg m^2 per m [UNVERIFIED].")
    print(f"  Plates from the centerline: support at {y_sup:.2f} m, load at {y_load:.2f} m, lift at {y_lift:.2f} m (outer edges at the cable edge {CABLE_FACE/2:.2f} m); cable")
    print(f"  roll inertia {MU_C*(2*CABLE_FACE**2)/12:,.0f} kg m^2 per m (square section), no torsional stiffness.  Per-side regulation: each load track holds its own side.")
    _, hist = pz.point_pulsed([lane(u30, label="pair"), lane(u30, hd_width=lift1 / P12, hu_width=load1 / P12, label="one")], u30, tabs)
    print(f"  {'case':>62} {'M kN m/m':>8} {'roll mrad':>9} {'edge mm':>7} {'gap_b':>6} {'gap_t':>6}")
    for lab, kw in (("one lift track alone, no per-side regulation", dict(with_prot=False)), ("one lift track alone, per-side regulation", dict(with_prot=True)),
                    ("one lift track alone, support lost too, per-side regulation", dict(with_prot=True, support_lost=True))):
        rr = ct.roll_model(hist, 1, u30, cap=load1, **kw)
        all_gaps.append((f"roll: {lab}", "2-DOF", rr["min_gap_b"], rr["min_gap_t"], rr["contact"]))
        print(f"  {lab:>62} {rr['peak_M']/1e3:8.0f} {rr['peak_roll']*1e3:9.2f} {rr['worst_cable_edge']*1e3:7.1f} {rr['min_gap_b']*1e3:6.1f} {rr['min_gap_t']*1e3:6.1f}{'  CONTACT' if rr['contact'] else ''}")
    # F.4 summary
    print()
    contacts = [g for g in all_gaps if g[4]]
    single = [g for g in all_gaps if "support lost too" not in g[0] and "no per-side" not in g[0] and "no protection" not in g[0]]
    print(f"  (the rows 'no per-side regulation' and 'no protection' are what-ifs that show what the protection is for; they are not design cases)")
    ws = min(single, key=lambda g: min(g[2], g[3]))
    print(f"  F.4  {len(all_gaps)} cases; contact in {len(contacts)}" + (": " + "; ".join(g[0] for g in contacts) if contacts else "") + ".")
    print(f"  Smallest gap in any single failure: {min(ws[2], ws[3])*1e3:.1f} mm in '{ws[0]}' ({ws[1]}).")
    worst_all = min(all_gaps, key=lambda g: min(g[2], g[3]))
    print(f"  Smallest gap in any case at all: {min(worst_all[2], worst_all[3])*1e3:.1f} mm in '{worst_all[0]}' ({worst_all[1]}).")
    single_contact = [g for g in single if g[4]]
    print(f"  Does any single failure lead to contact: {'YES: ' + '; '.join(g[0] for g in single_contact) if single_contact else 'no'}.")
    return all_gaps


# =============================================================================

def main():
    parts = "ABCDEF"
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
        elif arg == "--resweep":
            if os.path.exists(SWEEP_JSON):
                os.remove(SWEEP_JSON)
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    hr("BEARING FIELD DESIGN, VERIFICATION AND SIZING (2026-10-12)")
    print(f"  Ring at 700 km: V_rel {V_REL:,.0f} m/s, casing {MU:,.0f} kg/m ({W/1e3:.1f} kN/m), cable {MU_C:,.0f} kg/m, cable face {CABLE_FACE:.2f} m.")
    print(f"  Field solver: magnetostatic_2d.py with the sparse direct solve (solver '{ms.SOLVER}'); verification in bearing_fem_check.txt.")
    recs = part_B()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    theva = read_theva()
    if "A" in parts:
        part_A(recs)
    C = part_C(recs, theva) if ("C" in parts or "D" in parts or "F" in parts) else None
    D = part_D(recs, theva, C, tabs) if ("D" in parts or "E" in parts or "F" in parts) else None
    if "E" in parts:
        part_E(D, theva)
    if "F" in parts:
        part_F(C, D, tabs)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
