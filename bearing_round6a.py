#!/usr/bin/env python3
"""
The bearing, round 6, part A (2026-11-14): the plate spacing and the load-support coupling (Section 3 of the brief), and which fast coils
keep their leads during the deployment (Section 4).
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 14 (the cold loop), Chapter 15 (the deployment).
Brief: Vol_III/Claude_Code_instructions_bearing_round6_2026-11-14.md.
Run (from Python/orbitalring):
  python bearing_round6a.py --part=spacing          the coupling against the spacing, round 5's plates, 0.2 m to the largest the underside allows
  python bearing_round6a.py --part=roll             the roll case (one lift track alone) with the support tracks at each spacing
  python bearing_round6a.py --part=cases --d=0.2,0.6   round 4's 10 ms cases rerun with the coupling carried at those spacings
  python bearing_round6a.py --part=deploy           the deployment's gap control, the coil count, the leads' heat
  python bearing_round6a.py --part=all
Reads runs/bearing_round5_final.json (the optimum plates).  Writes runs/round6_logs/*.json and the case logs runs/bearing_round6a_casing_*.txt.
Environment: ORBITAL_RING_OUT (default runs).
"""
import contextlib
import io
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_layout as bl
import bearing_round4 as r4
import bearing_round5 as r5
import compensation_tracks as ct
import final_bearing as fb
import sled_passage_pulsed as pz
import sled_passage_simulation as sp

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
LOGS = os.path.join(OUT, "round6_logs")
SIDE = ring.CABLE_SIDE                                     # 7.39 m (the book's hardware; the brief says 7.39)
W_HALF = be.W_HALF                                        # one support track's share of the casing weight, 47.5 kN/m
K_FLUX = r4.K_FLUX                                        # 0.33 kN/m per mm per track, the constant-flux stiffness
TAU = 10e-3
L_FAST, N_FAM = 25.0, 20
N_M = 3.2                                                 # fast coils per meter of ring (four tracks x 20 families / 25 m)
COILS_PER_SITE = 1600
D_LIST = (0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2)


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def load_final():
    fin = json.load(open(os.path.join(OUT, "bearing_round5_final.json"), encoding="utf-8"))
    return fin["results"]["support:opt135"], fin["results"]["fast:opt135"], fin


def geometry(d, rs, rl):
    """The plates on the underside: the load plates at the edges, the support plates inboard with d between plate edges."""
    x_load = SIDE / 2 - rl["Wp"] / 2
    x_sup = SIDE / 2 - rl["Wp"] - d - rs["Wp"] / 2
    return dict(d=d, x_load=x_load, x_sup=x_sup, x_lift=x_load, sep=2 * x_sup, free=2 * x_sup - rs["Wp"])


def d_max(rs, rl):
    return (SIDE - 2 * rl["Wp"] - 2 * rs["Wp"]) / 2


def d_list(rs, rl):
    dm = d_max(rs, rl)
    return [d for d in D_LIST if d < dm - 1e-9] + [round(dm, 4)]


# =============================================================================
# S3.1  the coupling against the spacing
# =============================================================================

def part_spacing():
    rs, rl, fin = load_final()
    hr("S3.1  THE LOAD-SUPPORT COUPLING AGAINST THE PLATE SPACING (round 5's plates: support %.1f m x %.0f mm with %.1f m loops, load %.1f m x %.0f mm with a %.1f m loop)"
       % (rs["Wp"], rs["T"] * 1e3, rs["Wl"], rl["Wp"], rl["T"] * 1e3, rl["Wl"]))
    dm = d_max(rs, rl)
    print(f"  The {SIDE:.2f} m underside: two load plates of {rl['Wp']:.1f} m at the edges, two support plates of {rs['Wp']:.1f} m inboard; {SIDE - 2*rl['Wp'] - 2*rs['Wp']:.2f} m free in all, so the")
    print(f"  spacing d between a load plate and the support plate beside it runs from 0.2 m (round 5, {geometry(0.2, rs, rl)['free']:.2f} m left between the support plates) to {dm:.3f} m")
    print(f"  (the support plates touching at the center).  Stage B's model: linear plates (mu {int(__import__('compensation_fields').PLATE_MU)}), 0.1 m bundles, 10 mm cells; the load track at its")
    print(f"  full {fb.F_LOAD_TRACK/1e3:.0f} kN/m ({rl['NI_req']*math.sqrt(fb.F_LOAD_TRACK/rl['F_req'])/1e3:.0f} kA) beside a support track at its normal {W_HALF/1e3:.1f} kN/m ({rs['NI_req']*math.sqrt(W_HALF/rs['F_req'])/1e3:.0f} kA); the support persistent, so")
    print(f"  its ampere-turns fall by M/L times the load's and its force by twice that fraction.  The force figures are of the support's normal force, not its cap.")
    print(f"  {'d m':>6} {'support centers m':>17} {'separation m':>12} {'free between m':>14} {'M/L %':>6} {'dF/F %':>7} {'kN/m per track':>14} {'kN/m pair':>9} {'as gap offset mm':>16} {'at half pull kN/m pair':>22} {'s':>5}")
    rows = []
    for d in d_list(rs, rl):
        t0 = time.time()
        g = geometry(d, rs, rl)
        cp = r5.coupling(rs, rl, d)
        dF = abs(cp["dF"]); per_track = dF * W_HALF; pair = 2 * per_track
        eq_mm = pair / (2 * K_FLUX) * 1e3
        rows.append(dict(d=d, x_sup=g["x_sup"], sep=g["sep"], free=g["free"], k=cp["k"], dF=dF, per_track=per_track, pair=pair, eq_mm=eq_mm, M=cp["M"], L_s=cp["L_s"],
                         NI_l_full=cp["NI_l_full"], NI_s_op=cp["NI_s_op"]))
        print(f"  {d:6.3f} {g['x_sup']:17.3f} {g['sep']:12.2f} {g['free']:14.2f} {cp['k']*100:6.2f} {dF*100:7.2f} {per_track/1e3:14.2f} {pair/1e3:9.2f} {eq_mm:16.1f} {pair*math.sqrt(0.5)/1e3:22.2f} {time.time()-t0:5.0f}", flush=True)
    print(f"  'as gap offset': the pair's force change over the support's constant-flux stiffness, 2 x {K_FLUX/1e6:.2f} kN/m per mm: the offset the support's regulation would read it as if it")
    print(f"  took the force change for a gap change.  Round 5 compared its 5.8 percent with 'Stage B's 0.9 percent at 0.5 m'; Stage B's 0.9 was the flux ratio M/L (0.87 percent),")
    print(f"  its force figure was 5.2 percent at 0.5 m with the 1.4 and 1.2 m plates (the round 5 report is corrected, Section 2 of this report).")
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, "round6a_spacing.json"), "w", encoding="utf-8") as f:
        json.dump(dict(rows=rows, d_max=dm, side=SIDE), f, indent=1)
    return rows


# =============================================================================
# S3.2  the roll case with the support tracks inboard
# =============================================================================

def casing_setup():
    """Round 4's casing dynamics setup: the 12,000 kg/m casing with the bearing inside it, 25 m fast coils in 20 families, the 10 ms rise."""
    A2, theva, sup, lift = r4.inside_A()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    fb.L_LIFT, fb.L_LOAD = L_FAST, L_FAST; fb.N_FAM = N_FAM; bd.N_FAM_F = N_FAM
    bd.TAU_HW = TAU; fb.TAU_HW = TAU
    fb.set_ring(A2["m_top"], A2["top"])
    return A2, theva, sup, lift, tabs


def roll_hist(A2, tabs, freeze=True):
    """The point-model histories bearing_design's F.3 feeds the roll model: the pair and one lift track alone at 30 km/s."""
    S1 = A2["dm_top"]["F_s"]; lift1 = fb.F_LIFT_TRACK; load1 = fb.F_LOAD_TRACK
    w_hd, w_hu = 2 * lift1 / bd.P12, 2 * load1 / bd.P12
    u30 = 30000.0

    def lane(**kw):
        d = dict(L_sc=bd.L_SC, tau=TAU, L_section=L_FAST / N_FAM, gap_dep=False, steady_loop_hz=0.0 if freeze else 1.0, T_spread=max(0.05, bd.L_SC / u30), B_sat=1.2,
                 n_coils=N_FAM, L_coil=L_FAST, hd_width=w_hd, hu_width=w_hu, cmd_ramp=fb.RAMP, steady_cap=2 * S1)
        d.update(kw)
        return pz.lane(**d)
    _, hist = pz.point_pulsed([lane(label="pair"), lane(hd_width=lift1 / bd.P12, hu_width=load1 / bd.P12, label="one")], u30, tabs)
    return hist, u30, load1


def part_roll():
    rs, rl, fin = load_final()
    hr("S3.2  THE ROLL CASE (one lift track alone at 30 km/s, per-side regulation) WITH THE SUPPORT TRACKS AT EACH SPACING")
    A2, theva, sup, lift, tabs = casing_setup()
    hist, u30, load1 = roll_hist(A2, tabs)
    side = 10.0; r_g = side / math.sqrt(3.0)
    ct.I_CASING = bd.MU * r_g**2; ct.I_CABLE = bd.MU_C * (2 * SIDE**2) / 12.0
    ct.Y_PLATE_EDGE = SIDE / 2; ct.Y_CABLE_EDGE = SIDE / 2
    print(f"  bearing_design F.3's 2-DOF heave-and-roll model (casing a {side:.0f} m square truss, radius of gyration {r_g:.2f} m, {bd.MU:,.0f} kg/m; cable {bd.MU_C:,.0f} kg/m, no torsional")
    print(f"  stiffness [UNVERIFIED]): one lift track pulls alone at {geometry(0.2, rs, rl)['x_lift']:.2f} m from the centerline, the support tracks at +-x_sup each with a 1 Hz loop on its own gap,")
    print(f"  the load tracks and the other lift track as per-side protection at 10 mm.  The support's roll authority goes as x_sup: the roll moment it can resist per unit")
    print(f"  of its force change falls as the tracks move inboard.  Round 4's figure at d = 0.2 m: 61.4 mm (per-side regulation).")
    print(f"  {'d m':>6} {'x_sup m':>8} {'separation m':>12} {'peak M kN m/m':>13} {'roll mrad':>9} {'cable edge mm':>13} {'gap_b mm':>8} {'gap_t mm':>8} {'no per-side regulation: gap mm':>30}")
    rows = []
    for d in d_list(rs, rl):
        g = geometry(d, rs, rl)
        ct.Y_SUPPORT = g["x_sup"]; ct.Y_COMP = g["x_lift"]
        rr = ct.roll_model(hist, 1, u30, with_prot=True, cap=load1)
        rn = ct.roll_model(hist, 1, u30, with_prot=False, cap=load1)
        rows.append(dict(d=d, x_sup=g["x_sup"], sep=g["sep"], peak_M=rr["peak_M"], roll=rr["peak_roll"], edge=rr["worst_cable_edge"], gap_b=rr["min_gap_b"], gap_t=rr["min_gap_t"],
                         contact=rr["contact"], gap_noprot=min(rn["min_gap_b"], rn["min_gap_t"]), contact_noprot=rn["contact"]))
        print(f"  {d:6.3f} {g['x_sup']:8.3f} {g['sep']:12.2f} {rr['peak_M']/1e3:13.0f} {rr['peak_roll']*1e3:9.2f} {rr['worst_cable_edge']*1e3:13.1f} {rr['min_gap_b']*1e3:8.1f} {rr['min_gap_t']*1e3:8.1f}"
              f"{'  CONTACT' if rr['contact'] else ''} {min(rn['min_gap_b'], rn['min_gap_t'])*1e3:30.1f}{'  CONTACT' if rn['contact'] else ''}", flush=True)
    print(f"  What else depends on the support tracks' separation: the roll stiffness of the support pair at constant flux (2 K_FLUX x_sup^2 per radian, {2*K_FLUX*geometry(0.2, rs, rl)['x_sup']**2/1e3:,.0f} kN m per rad per m")
    print(f"  at 0.2 m), the lateral (sideways) centering the support loops give the casing, which this model does not carry, and the 3.79 m in the round 5 cross-section drawing.")
    with open(os.path.join(LOGS, "round6a_roll.json"), "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))
    return rows


# =============================================================================
# S3.3  the cases with the coupling carried
# =============================================================================

def f1_block(txt):
    lines = txt.splitlines()
    i0 = next(i for i, l in enumerate(lines) if "F.1  Sled passage" in l)
    out = {}
    for l in lines[i0 + 2:i0 + 9]:
        p = l.split()
        if len(p) >= 6:
            try:
                out[float(p[0])] = dict(F=float(p[1]), pair=float(p[2]), one=float(p[3]), gap_b=float(p[4]), gap_t=float(p[5]))
            except ValueError:
                pass
    return out


def tol_block(txt):
    lines = txt.splitlines()
    i0 = next(i for i, l in enumerate(lines) if "Tolerances (pair)" in l)
    out = {}
    for l in lines[i0 + 1:i0 + 11]:
        p = l.strip().rsplit(None, 4)
        if len(p) == 5:
            try:
                out[p[0]] = [float(x) for x in p[1:]]
            except ValueError:
                pass
    return out


def case_gaps(gaps, sub):
    hits = [g for g in gaps if sub in g[0]]
    return min(min(g[2], g[3]) for g in hits) * 1e3 if hits else float("nan")


def part_cases(ds=(0.2,), spacing_rows=None):
    rs, rl, fin = load_final()
    hr("S3.3  ROUND 4'S 10 MS CASES RERUN WITH THE LOAD-SUPPORT COUPLING CARRIED")
    A2, theva, sup, lift, tabs = casing_setup()
    if spacing_rows is None:
        p = os.path.join(LOGS, "round6a_spacing.json")
        spacing_rows = json.load(open(p, encoding="utf-8"))["rows"] if os.path.exists(p) else part_spacing()
    cs = {round(r["d"], 4): r["dF"] for r in spacing_rows}
    print(f"  bearing_layout.run_cases (final_bearing's Part B, the 87 cases of round 4) with pz.lane's new 'couple': while the floor tracks pull F_hu the support force is multiplied")
    print(f"  by 1 - c sqrt(F_hu / F_full), c the coupling of one support track at one load track's full pull (S3.1), F_full that pull (one track) or the pair's; the support frozen,")
    print(f"  the 10 ms rise, the 10 mm protection in every lane, 25 m coils in 20 families, the 12,000 kg/m casing.  The lift pulls (above orbital speed) couple to nothing: the")
    print(f"  lift tracks face the top plates, where there is no support track; the protection's up-pull is a load-track pull and couples in every case.")
    base_path = os.path.join(OUT, "bearing_round4_casing_10ms_20fam.txt")
    if not os.path.exists(base_path):
        base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Vol_III", "runs_coil_optimization_2026-10-13", "round4_2026-11-13", "bearing_round4_casing_10ms_20fam.txt")
    base_txt = open(base_path, encoding="utf-8").read()
    runs = {"round 4, no coupling": dict(txt=base_txt, summary=bl.parse_summary(base_txt), gaps=None, records=None, d=None)}
    for d in ds:
        c = cs.get(round(d, 4))
        if c is None:
            c = abs(r5.coupling(rs, rl, d)["dF"])
        fname = f"bearing_round6a_casing_10ms_d{d*100:.0f}cm.txt"
        t0 = time.time()
        gaps, records, txt = bl.run_cases(TAU, True, A2, theva, sup, lift, tabs, fname, prot_all=True, head=f"[round 6a: spacing {d:.2f} m, coupling {c*100:.2f} percent]", couple=c)
        runs[f"d = {d:.2f} m, c = {c*100:.2f} %"] = dict(txt=txt, summary=bl.parse_summary(txt), gaps=gaps, records=records, d=d, c=c, seconds=time.time() - t0, log=fname)
        print(f"  ran d = {d:.2f} m (c = {c*100:.2f} percent): {time.time()-t0:.0f} s -> {fname}", flush=True)
    # F.1 and the tolerances at the load-pull speeds
    print()
    print(f"  F.1, the load-pull speeds (pair / one track alone, the excursion in mm; the smallest gap of the two):")
    print(f"  {'run':>28} " + " ".join(f"{'%.0f km/s' % u:>22}" for u in (1, 3, 5)) + f" {'10':>14} {'30':>14}")
    for name, r in runs.items():
        f1 = f1_block(r["txt"])
        cells = []
        for u in (1.0, 3.0, 5.0, 10.0, 30.0):
            x = f1.get(u)
            cells.append(f"{x['pair']:.2f} / {x['one']:.2f}, gap {min(x['gap_b'], x['gap_t']):.1f}" if x else "n/a")
        print(f"  {name:>28} " + " ".join(f"{c:>22}" if i < 3 else f"{c:>14}" for i, c in enumerate(cells)))
    print(f"  Tolerances (pair), the excursion in mm at 1 / 3 km/s (the load pulls):")
    tols = {name: tol_block(r["txt"]) for name, r in runs.items()}
    labs = list(next(iter(tols.values())).keys())
    print(f"  {'tolerance':>16} " + " ".join(f"{name:>28}" for name in runs))
    for lab in labs:
        print(f"  {lab:>16} " + " ".join(f"{tols[name][lab][0]:.2f} / {tols[name][lab][1]:.2f}".rjust(28) if lab in tols[name] else "n/a".rjust(28) for name in runs))
    # the failure cases a load pull touches, and the smallest single-failure gap
    print()
    items = [("sled passage 1 km/s, one track", "one load track lost: the one-track pass at 1 km/s"), ("false full load pull 50 ms", "a false full load pull of 50 ms"),
             ("both support tracks lost, load pair catches", "both support tracks lost at rest, the load pair catches"), ("both support tracks lost, one load track catches", "the same, one load track"),
             ("1. one support track loses both windings at rest", "one support track lost at rest"), ("sensor error +5 mm", "a gap sensor +5 mm (the protection's up-pull couples)"),
             ("pull 50 ms late", "a pull 50 ms late (the smallest single-failure gap of round 4)"), ("spread pull fired with no sled", "the spread pull with no sled")]
    print(f"  The cases a load-track pull touches (the smallest gap of the case, mm):")
    print(f"  {'case':>62} " + " ".join(f"{name:>28}" for name in runs))
    base_gaps = None
    for sub, lab in items:
        cells = []
        for name, r in runs.items():
            if r["gaps"] is None:
                cells.append("(see the round 4 log)")
            else:
                cells.append(f"{case_gaps(r['gaps'], sub):.1f}")
        print(f"  {lab:>62} " + " ".join(f"{c:>28}" for c in cells))
    print(f"  Smallest single-failure gap: " + "; ".join(f"{name}: {r['summary']['gap']} mm ({r['summary']['case'][:40]})" for name, r in runs.items()))
    # what the support's regulation and the gap sensors see: the support force history in the F.1 lanes
    print()
    print(f"  What the support's regulation sees on a pass (the support force history of the F.1 lanes, the dip below the casing weight and how long it lasts),")
    print(f"  and what the gap sensors see (the casing's excursion of the same lane):")
    print(f"  {'run':>28} {'u km/s':>6} {'lane':>5} {'dip kN/m pair':>13} {'as gap offset mm':>16} {'dip lasts ms':>12} {'excursion mm':>12}")
    out = dict(runs={})
    for name, r in runs.items():
        if r["records"] is None:
            continue
        seen = {}
        for lab, u, t, S, res in r["records"]:
            key = (u, lab)
            if lab not in ("pair", "one") or u not in (1000.0, 3000.0, 5000.0, 10000.0, 30000.0) or key in seen:
                continue
            seen[key] = True
            W = pz.W; dip = W - S.min(); dt = t[1] - t[0]
            above = (W - S) > dip / 2
            dur = above.sum() * dt
            print(f"  {name:>28} {u/1e3:6.0f} {lab:>5} {dip/1e3:13.2f} {dip/(2*K_FLUX)*1e3:16.1f} {dur*1e3:12.0f} {res['exc']*1e3:12.2f}")
            out["runs"].setdefault(name, []).append(dict(u=u, lane=lab, dip=dip, dur=dur, exc=res["exc"], gap_b=res["gap_b"], gap_t=res["gap_t"]))
        out["runs"][name + " summary"] = r["summary"]
    for name, r in runs.items():
        out["runs"][name + " f1"] = f1_block(r["txt"]); out["runs"][name + " tol"] = tol_block(r["txt"])
    with open(os.path.join(LOGS, "round6a_cases.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)
    return runs


# =============================================================================
# S4  the deployment's gap control and the leads
# =============================================================================

def part_deploy():
    rs, rl, fin = load_final()
    hr("S4  WHICH FAST COILS KEEP THEIR LEADS DURING THE DEPLOYMENT")
    rows, li = r5.lead_rows(rl)
    r10 = rows[10]; I_full = r10["I"]; A_m = r10["A_m"]; R_lead = be.RHO_CU_77 * 1.0 / A_m
    P_full = I_full**2 * R_lead
    q_idle = r10["q_idle_lead"]
    cold = json.load(open(os.path.join(OUT, "bearing_round5_coldload.json"), encoding="utf-8"))
    idle_site_all = cold["leads"]["10"]["idle_site"]
    W = ring.M_LOAD_M * ring.G_NET
    F_pair_per_mm = 2 * K_FLUX * 1e-3                         # N/m per mm of offset, the pair (0.66 kN/m per mm)
    print(f"  The support is attractive and unstable: at constant flux its stiffness is -{K_FLUX/1e6:.2f} kN/m per mm per track, -{F_pair_per_mm/1e3:.2f} kN/m per mm for the pair, e-folding")
    print(f"  {math.sqrt(ring.M_LOAD_M/(2*K_FLUX))*1e3:.0f} ms on the {ring.M_LOAD_M:,.0f} kg/m casing (round 4, 7.4).  The fast loops hold it: a floor (load) track pulls the casing toward the cable when it has sagged,")
    print(f"  a ceiling (lift) track pulls it back when it has risen; the pull a track needs for an offset z is {K_FLUX/1e6:.2f} kN/m per mm of z.  Full pull of a track {fb.F_LOAD_TRACK/1e3:.0f} kN/m at")
    print(f"  {I_full:,.0f} A per turn (10 turns), rating 1.46 x at 1.21 x the current; the force goes as the ampere-turns squared, so n of the 20 overlapping coils at a current")
    print(f"  fraction f give (n f / 20)^2 of the full pull, and a pull p (fraction of full) needs f = sqrt(p) x 20 / n per connected coil.")
    # S4.1 the disturbances
    print()
    print(f"  S4.1  The disturbances during the deployment (no launches), and the pull each takes:")
    GM = ring.GM; om_e = 7.2921159e-5
    r0 = 6905499.5; r1 = ring.R_ORBIT
    def g_eff(r):
        return GM / r**2 - (om_e * r)**2 / r
    W0 = ring.M_LOAD_M * g_eff(r0); W1 = ring.M_LOAD_M * g_eff(r1); days = 381.5
    dW_day = (W0 - W1) / days
    dist = []
    dist.append(("(a) the instability on an offset z: the support's negative stiffness", f"{F_pair_per_mm/1e3:.2f} kN/m per mm (pair)", f"{K_FLUX/1e6:.2f} kN/m per mm per track; at the 10 mm protection band {10*K_FLUX*1e-3/1e3:.1f} kN/m per track, {10*K_FLUX*1e-3/fb.F_LOAD_TRACK*100:.2f} % of full"))
    dist.append(("(b) the casing's weight falling as the ring rises (g_eff from %.3f to %.3f m/s2, %.0f to %.0f km)" % (g_eff(r0), g_eff(r1), (r0-6378.137e3)/1e3, (r1-6378.137e3)/1e3),
                 f"{(W0-W1)/1e3:.2f} kN/m over the {days:.0f} days, {dW_day:.1f} N/m per day", f"between trims of the support: a daily trim leaves {dW_day/2/1e3:.3f} kN/m per track, a weekly {7*dW_day/2/1e3:.2f}, a monthly {30*dW_day/2/1e3:.2f}"))
    dist.append(("(c) one support track lost (a quench, a cryostat) at rest", f"{W_HALF/1e3:.1f} kN/m, {W_HALF/2/1e3:.1f} per load track", f"until the other support track's supply has doubled its current: round 4's 1 percent per second is 100 s; its 1 kV dump-rated supply on the {81**2/55**2*21.3:.0f} H coil, {81**2/55**2*21.3*455/1000:.0f} s"))
    dist.append(("(d) the cable's shape passing at %.0f m/s: a bump of height b and length longer than the loop can follow" % ring.V_REL, "0.66 kN/m per mm of b (pair), as (a)",
                 f"the loop's 10 ms rise and 5 Hz follow features longer than about {ring.V_REL*0.2:,.0f} m; shorter ones must lie inside the band [the plates' straightness is a requirement, not a known value]"))
    dist.append(("(e) a connected coil failing", "its share: 1/n of the connected pull", "the one failure allowed for: n + 1 connected"))
    print(f"  {'disturbance':>100} | {'size':>45} | note")
    for a, b, c in dist:
        print(f"  {a:>100} | {b:>45} | {c}")
    # S4.2 the coil count, the current per coil, what the pulse-sized lead can carry
    print()
    print(f"  S4.2  The connected coils.  The pulse-sized lead ({A_m*1e6:.0f} mm2 of copper, 1 m, {R_lead*1e3:.3f} mohm at 77 K): Joule heat at full current {P_full:.0f} W; its conduction alone {q_idle:.2f} W into the cold.")
    c_cu = 400.0; m_lead = 8960.0 * A_m * 1.0
    print(f"  Adiabatic: a 20 K rise takes {20*c_cu*m_lead/P_full:.1f} s at full current ({c_cu:.0f} J/kg K as round 4's rule, {m_lead*1e3:.0f} g of copper), {20*c_cu*m_lead/P_full:.1f} / f^2 s at the fraction f.  Steady: the middle of a")
    print(f"  1 m lead cooled at both ends rises P / (8 k A / l) = {1/(8*400*A_m):.1f} K per W, so 20 K allows {20/(1/(8*400*A_m)):.2f} W of Joule heat, a steady current of {math.sqrt(20/(1/(8*400*A_m))/P_full)*100:.1f} % of full.")
    print(f"  {'pull per track':>30} {'fraction p':>10} {'n connected':>11} {'f per coil %':>12} {'A per coil':>10} {'Joule W per lead':>16} {'20 K in s':>9} {'steady?':>7}")
    pulls = [("(a) 1 mm offset, 0.33 kN/m", K_FLUX * 1e-3), ("(a) 5 mm, 1.65 kN/m", 5 * K_FLUX * 1e-3), ("(a) 10 mm band, 3.3 kN/m", 10 * K_FLUX * 1e-3),
             ("(b) weekly trim, %.2f kN/m" % (7 * dW_day / 2 / 1e3), 7 * dW_day / 2), ("(c) one support track lost, %.1f kN/m" % (W_HALF / 2 / 1e3), W_HALF / 2)]
    for lab, F in pulls:
        p = F / fb.F_LOAD_TRACK
        for n in (20, 10, 5, 4, 3):
            f = math.sqrt(p) * 20 / n
            PJ = f**2 * P_full; t20 = 20 * c_cu * m_lead / PJ if PJ > 0 else float("inf")
            steady = "yes" if PJ <= 20 / (1 / (8 * 400 * A_m)) else "no"
            print(f"  {lab:>30} {p*100:9.3f}% {n:11d} {f*100:12.1f} {f*I_full:10,.0f} {PJ:16.2f} {min(t20, 9e5):9.1f} {steady:>7}")
    # S4.3 the heat per site against n for a stated regulation duty, and the deployment
    print()
    print(f"  S4.3  The number that costs least heat.  With n of the 20 families connected (a fixed subset: a family tiles the track end to end, so n families put n coils over")
    print(f"  every point but the family's own coil ends, where n - 1; rotation buys nothing, the conduction is there while a lead is connected) the conduction is n/20 of all")
    print(f"  leads' {idle_site_all/1e3:.2f} kW per site, and the regulation's Joule heat goes as (20/n)^2 per coil, (20/n) per site: for a regulation that keeps the casing within a mean")
    print(f"  offset <|z|> the mean pull fraction is <p> = {K_FLUX*1e-3/fb.F_LOAD_TRACK*100:.3f} % per mm, the Joule heat per lead <p> (20/n)^2 x {P_full:.0f} W, plus the weekly trim's standing pull.")
    print(f"  {'<|z|> mm':>8} {'n':>3} {'conduction kW/site':>18} {'Joule kW/site':>13} {'total kW/site':>13} {'per lead: cond W':>16} {'Joule W':>8} {'steady 20 K?':>12}")
    best = {}
    table = []
    for z in (0.1, 0.3, 1.0):
        p_mean = K_FLUX * 1e-3 * z / fb.F_LOAD_TRACK + 7 * dW_day / 2 / fb.F_LOAD_TRACK
        for n in (20, 14, 10, 7, 5, 4, 3):
            cond = idle_site_all * n / 20
            PJ = p_mean * (20 / n)**2 * P_full
            joule = COILS_PER_SITE * n / 20 * 2 * PJ
            tot = cond + joule
            ok = PJ <= 20 / (1 / (8 * 400 * A_m))
            table.append(dict(z=z, n=n, cond=cond, joule=joule, tot=tot, PJ=PJ, ok=ok))
            if z not in best or tot < best[z][1]:
                best[z] = (n, tot)
            print(f"  {z:8.1f} {n:3d} {cond/1e3:18.2f} {joule/1e3:13.2f} {tot/1e3:13.2f} {q_idle:16.2f} {PJ:8.2f} {'yes' if ok else 'no':>12}")
    print(f"  Least heat: " + "; ".join(f"<|z|> = {z} mm: n = {n} ({tot/1e3:.2f} kW per site)" for z, (n, tot) in best.items()))
    # S4.4 leads sized for the deployment's currents
    print()
    print(f"  S4.4  Leads sized for the deployment's currents instead of the pulse.  For a steady current I the conduction-cooled copper lead has an optimum A/l at which its heat into the")
    print(f"  cold is least (be.lead_heat, Wiedemann-Franz, k constant); the pulse-sized lead is that optimum for about the current at which its own conduction and Joule heat are equal.")
    print(f"  {'steady A per coil':>17} {'as % of full':>12} {'pulse-sized lead W':>18} {'optimum lead W':>14} {'optimum mm2':>11} {'which':>30}")
    for frac in (0.01, 0.02, 0.05, 0.1, 0.2, 0.35):
        I = frac * I_full
        q_p = be.lead_heat(I, A_over_l=A_m)[0]
        q_o, aol = be.lead_heat(I)
        print(f"  {I:17,.0f} {frac*100:12.0f} {q_p:18.2f} {q_o:14.2f} {aol*1e6:11.1f} {('the pulse-sized lead is thinner' if aol > A_m else 'the pulse-sized lead is thicker'):>30}")
    print(f"  The deployment's regulation current per connected coil is f = sqrt(<p>) x 20 / n: at <|z|> = 0.1 mm and n = 5, {math.sqrt(K_FLUX*1e-3*0.1/fb.F_LOAD_TRACK)*20/5*100:.1f} % of full, {math.sqrt(K_FLUX*1e-3*0.1/fb.F_LOAD_TRACK)*20/5*I_full:,.0f} A; the pulse-sized lead's")
    print(f"  own heat there is set by its conduction, and a lead sized for that current would be thinner and leak less by the ratio of the sections, but could not carry the (c)")
    print(f"  hold or a launch; the connected coils' leads stay pulse-sized and the saving is in n, not in the lead.")
    # the deployment runs
    out = dict(I_full=I_full, A_m=A_m, R_lead=R_lead, P_full=P_full, q_idle=q_idle, idle_site_all=idle_site_all, W0=W0, W1=W1, dW_day=dW_day, table=table, best={str(k): v for k, v in best.items()})
    runs = {}
    for z, (n, tot) in best.items():
        runs[f"n{n}"] = dict(n=n, z=z, extra=tot)
    runs["n5"] = dict(n=5, z=0.1, extra=next(r["tot"] for r in table if r["z"] == 0.1 and r["n"] == 5))
    out["deploy_runs"] = runs
    print()
    print(f"  Deployment runs (tools/round5_deployment.py --custom=tag:extra_W, the optimum plates' hardware and 4.0 m cold width, the extra cold load per site):")
    for tag, r in runs.items():
        print(f"    {tag}: {r['n']} families connected, <|z|> = {r['z']} mm: {r['extra']:,.0f} W per site (conduction {idle_site_all*r['n']/20:,.0f} + Joule {r['extra']-idle_site_all*r['n']/20:,.0f})")
    with open(os.path.join(LOGS, "round6a_deploy.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    return out


def main():
    part = "all"; ds = (0.2,)
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--d="): ds = tuple(float(x) for x in a[4:].split(","))
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    print("=" * 150 + "\n  THE BEARING, ROUND 6 PART A (2026-11-14)\n" + "=" * 150)
    print(f"  Ring at 700 km, rising quadratic: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m with hardware {ring.M_HARDWARE_M:,.0f}, side {SIDE:.2f} m; the rise {TAU*1e3:.0f} ms.")
    t0 = time.time()
    if part in ("spacing", "all"):
        part_spacing()
    if part in ("roll", "all"):
        part_roll()
    if part in ("cases", "all"):
        part_cases(ds)
    if part in ("deploy", "all"):
        part_deploy()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
