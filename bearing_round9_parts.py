#!/usr/bin/env python3
"""
Round 9 (2026-11-16): the tables (run through bearing_round9.py --part=ontime|heat|table75|failures|lost).
Brief: Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md.  Reads runs/round7_logs/, runs/round8_logs/, runs/round9_logs/.
Writes runs/round9_logs/round9_*.json and prints the tables (redirect to runs/round9_logs/round9_*.txt).
"""
import json
import math
import os

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round5 as r5
import bearing_round6b as r6
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_round8 as r8
import bearing_round9 as r9
import final_bearing as fb

MU0 = 4e-7 * math.pi
COP = 0.0106; E_LAUNCH = 4.6e15
SPEEDS = (1, 3, 5, 10, 20, 30)


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


# =============================================================================
# correction 1: the fast coils' on-time
# =============================================================================

def part_ontime():
    hr("CORRECTION 1: THE FAST COILS' ON-TIME (the energized length over the sled's speed)")
    fast = r7.load("round7_fast_track.json"); ld = p7.lead(); A = ld["A"]
    out = {}
    print(f"  A fast coil is on the casing.  The pulled zone (the energized length L_on = u max(50 ms, L_sc / u) + L_sc + 2 L_coil + 2 u tau, with L_sc = {bd.L_SC:.0f} m, 25 m coils and the 10 ms")
    print(f"  rise) moves along the casing at the sled's speed u, so a coil is on for L_on / u.  Rounds 4 to 8 divided by u + V_rel ({bd.V_REL:,.0f} m/s), which is the time a point of the plate")
    print(f"  on the cable spends in the zone.  The plate's dwell is right for the plate's heat and stays there (bearing_round5.plate_heat_launch, bearing_round7_parts heat).")
    print(f"  {'u km/s':>7} {'pull per track kN/m':>19} {'of full':>8} {'zone m':>7} {'on-time, rounds 4-8 ms':>22} {'on-time, corrected ms':>21} {'ratio':>6} {'full-current s, old':>19} {'corrected':>10}")
    rows = []
    for uk in SPEEDS:
        u = uk * 1e3; F, kind = r4.pull_track_of_u(u); fr = F / fast["F_req"]
        t_new, L_on = r5.hold_time(u, frame="casing"); t_old = r5.hold_time(u, frame="plate")[0]
        rows.append(dict(u=u, F=F, kind=kind, frac=fr, L_on=L_on, t_old=t_old, t_new=t_new, s_old=fr * t_old, s_new=fr * t_new))
        print(f"  {uk:7d} {F/1e3:19.1f} {fr:8.2f} {L_on:7.0f} {t_old*1e3:22.1f} {t_new*1e3:21.1f} {t_new/t_old:6.2f} {fr*t_old:19.4f} {fr*t_new:10.4f}")
    out["rows"] = rows
    # the launch integrals, both frames
    li_o = r5.launch_integrals(fast["W"], fast["F_req"], frame="plate"); li_n = r5.launch_integrals(fast["W"], fast["F_req"], frame="casing")
    P_full = fast["I"]**2 * be.RHO_CU_77 * 1.0 / A
    E_o = li_o["lead_s"] * 2 * P_full; E_n = li_n["lead_s"] * 2 * P_full
    n_coils = r5.N_M * ring.L_RING
    print(f"\n  Over the heaviest launch (sled_passage_pulsed.lsm_profile, {li_n['x_run']/1e3:,.0f} km of run in {li_n['t_run']/60:.0f} minutes): the coils' full-current seconds summed around the ring")
    print(f"  {li_o['lead_s']:.4e} as rounds 4 to 8 had it and {li_n['lead_s']:.4e} corrected, {li_n['lead_s']/li_o['lead_s']:.3f} times.  Two leads per coil at {P_full:.0f} W each at full current")
    print(f"  ({fast['I']:,.0f} A in {A*1e6:.1f} mm2 by 1 m of copper at 77 K): the leads' Joule heat per launch {E_o/1e9:.1f} GJ before, {E_n/1e9:.1f} GJ corrected.")
    out.update(lead_s_old=li_o["lead_s"], lead_s_new=li_n["lead_s"], P_full=P_full, E_lead_old=E_o, E_lead_new=E_n, A_lead=A, x_run=li_n["x_run"], t_run=li_n["t_run"])
    # the lead against the corrected on-times
    print(f"\n  The 21 mm2 lead against the corrected on-times.  Its budget is the full-current time that warms the copper 20 K with no heat leaving (20 K x 400 J/kg K x 8,960 kg/m3 x A x 1 m")
    print(f"  = {ld['E20']:,.0f} J over {ld['P_full']:.0f} W): {ld['budget20']:.2f} s with all 20 families connected and {ld['budget10']:.2f} s with ten (each connected coil then carries twice the current).")
    print(f"  {'u km/s':>7} {'full-current s a pass uses':>26} {'of the 20-family budget':>23} {'copper warms K':>14} {'rounds 4-8: s':>13} {'of the budget':>13}")
    use = []
    for r in rows:
        if r["u"] in (1e3, 30e3):
            use.append(dict(u=r["u"], s=r["s_new"], share=r["s_new"] / ld["budget20"], dT=20 * r["s_new"] / ld["budget20"], s_old=r["s_old"], share_old=r["s_old"] / ld["budget20"]))
            print(f"  {r['u']/1e3:7.0f} {r['s_new']:26.3f} {r['s_new']/ld['budget20']*100:22.1f}% {20*r['s_new']/ld['budget20']:14.2f} {r['s_old']:13.4f} {r['s_old']/ld['budget20']*100:12.1f}%")
    out["budget"] = dict(b20=ld["budget20"], b10=ld["budget10"], use=use)
    print(f"  With ten families a launch is not flown (all 20 are connected for a launch); if it were, a pass would use four times these shares of the 20-family budget, which is the")
    print(f"  same seconds against {ld['budget10']:.2f} s: {rows[0]['s_new']/ld['budget10']*100:.0f} percent at 1 km/s and {rows[-1]['s_new']/ld['budget10']*100:.0f} percent at 30 km/s.")
    # the mean rise at 500 and at 14 launches a day
    print(f"\n  The lead's mean warming.  Round 5 thickened the pulse-sized lead until the day's mean Joule heat lifted its middle by 20 K at {r5.LAUNCHES_PER_DAY} launches a day (dT = P_mean l / (8 k A)")
    print(f"  written in the script as P_mean x 1 m / (2 x 400 x A), kept as it is).  With the corrected duty the same lead:")
    print(f"  {'launches a day':>14} {'duty, full-current ms per s per coil: old':>40} {'corrected':>10} {'mean rise K: old':>16} {'corrected':>10} {'lead by the same rule mm2':>25}")
    rl = r6.load_final()[1]; I5 = r5.fast_basics(rl, 10)["I"]; A_p = I5 / r4.J_PULSE_LEAD
    mr = []
    for n_day in (r5.LAUNCHES_PER_DAY, 14):
        d_o = li5(rl, "plate") * n_day / 86400 / n_coils; d_n = li5(rl, "casing") * n_day / 86400 / n_coils
        dT = lambda duty, Ax: I5**2 * (be.RHO_CU_77 / Ax) * duty / (2 * 400.0 * Ax)
        A_rule = A_p * math.sqrt(dT(d_n, A_p) / 20.0) if dT(d_n, A_p) > 20 else A_p
        mr.append(dict(n=n_day, duty_old=d_o, duty_new=d_n, dT_old=dT(d_o, A), dT_new=dT(d_n, A), A_rule=A_rule))
        print(f"  {n_day:14d} {d_o*1e3:40.3f} {d_n*1e3:10.3f} {dT(d_o, A):16.2f} {dT(d_n, A):10.2f} {A_rule*1e6:25.1f}")
    out["mean_rise"] = mr
    print(f"  (the rule was applied on round 5's optimum plate, {I5:,.0f} A; the lead is kept at {A*1e6:.1f} mm2 as ruled.  At {r5.LAUNCHES_PER_DAY} of the heaviest launches a day it would warm {mr[0]['dT_new']:.0f} K, not 20;")
    print(f"  at 14 a day {mr[1]['dT_new']:.2f} K.  The conduction with no current and the 20 K budget do not depend on the on-time.)")
    # every number of rounds 4 to 8 that the on-time moves
    print(f"\n  What the correction moves in rounds 4 to 8 (old to corrected; everything else in those rounds does not use a coil's on-time):")
    r_sup, r_lift = be.load_records()
    l4o = r4.launch_integrals(r_lift, frame="plate"); l4n = r4.launch_integrals(r_lift, frame="casing")
    n4 = 2 * r4.N_FAM / r4.L_FAST * ring.L_RING
    moved = [("round 4, S3.1: the coils' full-current on-time summed over the ring, million coil-seconds per launch", l4o["lead_s"] / 1e6, l4n["lead_s"] / 1e6, "{:.1f}"),
             ("round 4, S3.1: the same per coil, averaged over all coils, ms", l4o["lead_s"] / n4 * 1e3, l4n["lead_s"] / n4 * 1e3, "{:.0f}"),
             ("round 4, S3.1: a lift coil at 30 km/s is on for, ms", l4o["per_u"][30]["hold"] * 1e3, l4n["per_u"][30]["hold"] * 1e3, "{:.0f}"),
             ("round 4, S3.1 and 4.2; round 5, S6: the leads' pulse heat (every 'pulses' column) and the duty, as a ratio", 1.0, li5(rl, "casing") / li5(rl, "plate"), "{:.3f}"),
             ("round 5, S6: the duty at 500 launches a day, full-current ms per second per coil", mr[0]["duty_old"] * 1e3, mr[0]["duty_new"] * 1e3, "{:.2f}"),
             ("round 5, S6: the leads' pulse heat at 500 launches a day, W per meter of ring (10 turns)", pulse_m(rl, "plate"), pulse_m(rl, "casing"), "{:.2f}"),
             ("round 5, S6: the same per 500 m site, kW", pulse_m(rl, "plate") * 0.5, pulse_m(rl, "casing") * 0.5, "{:.2f}"),
             ("round 5, S4: a switch's hold at 30 km/s, ms (the thermal rating is still not reached)", r5.hold_time(30e3, frame="plate")[0] * 1e3, r5.hold_time(30e3, frame="casing")[0] * 1e3, "{:.0f}"),
             ("round 6b, Section 2: the leads' heat per launch in the table of every heat (each plate's row), as a ratio", 1.0, li_n["lead_s"] / li_o["lead_s"], "{:.3f}"),
             ("round 8: the leads' Joule heat of the heaviest launch, GJ", E_o / 1e9, E_n / 1e9, "{:.1f}"),
             ("round 8: the cold heat the heaviest launch leaves (tape and leads), GJ", (fast["Ep"] + fast["Es"] + E_o) / 1e9, (fast["Ep"] + fast["Es"] + E_n) / 1e9, "{:.0f}"),
             ("round 8: the same per 500 m site as a mean, MJ", (fast["Ep"] + fast["Es"] + E_o) / ring.L_RING * 500 / 1e6, (fast["Ep"] + fast["Es"] + E_n) / ring.L_RING * 500 / 1e6, "{:.2f}"),
             ("gap trade, item 5: the lowest speed at which a load coil under the spacecraft is inside the 20-family budget, m/s (already by the casing's frame)", 49.0, 49.0, "{:.0f}")]
    for lab, a, b, f in moved:
        print(f"    {lab:<150} {f.format(a):>9} -> {f.format(b):>9}")
    out["moved"] = [dict(label=l, old=a, new=b) for l, a, b, f in moved]
    print(f"  Not moved: the tape's pulse heat (a loss per cycle, {(fast['Ep']+fast['Es'])/1e9:.0f} GJ), the plate heat (the plate's dwell is the right time there, 247 GJ), the stores and the pulse power (the zone's")
    print(f"  front advances at u in both), the lead's section and its conduction (kept as ruled), the 20 K budgets, the change-over of a lost loop (it is the integral of the pull, not a pass),")
    print(f"  and the cases' movements (sled_passage_pulsed has its own timing and never used hold_time).")
    r9.save("round9_ontime.json", out)
    return out


def li5(rl, frame):
    return r5.launch_integrals(rl["W_mag"], rl["F_req"], frame=frame)["lead_s"]


def pulse_m(rl, frame):
    """Round 5's pulse heat per meter of ring at 500 launches a day, 10 turns, in the lead as ruled."""
    I = r5.fast_basics(rl, 10)["I"]; A = r5.lead_rows(rl)[0][10]["A_m"]
    duty = li5(rl, frame) * r5.LAUNCHES_PER_DAY / 86400 / (r5.N_M * ring.L_RING)
    return 2 * I**2 * be.RHO_CU_77 / A * duty * r5.N_M


# =============================================================================
# item 4: the heat of a launch
# =============================================================================

def part_heat():
    hr("ITEM 4: THE COLD HEAT OF THE HEAVIEST LAUNCH")
    fast = r7.load("round7_fast_track.json"); A = p7.lead()["A"]; P_full = fast["I"]**2 * be.RHO_CU_77 * 1.0 / A
    li_o = r5.launch_integrals(fast["W"], fast["F_req"], frame="plate"); li = r5.launch_integrals(fast["W"], fast["F_req"], frame="casing")
    E_tape = fast["Ep"] + fast["Es"]; E_lo = li_o["lead_s"] * 2 * P_full; E_l = li["lead_s"] * 2 * P_full
    laps = li["x_run"] / ring.L_RING
    print(f"  {'':>44} {'round 8':>12} {'round 9':>12}")
    rows = [("the tape's pulse heat, GJ", E_tape / 1e9, E_tape / 1e9), ("the leads' Joule heat, GJ", E_lo / 1e9, E_l / 1e9), ("their sum, GJ", (E_tape + E_lo) / 1e9, (E_tape + E_l) / 1e9),
            ("the compressor's work for it, TJ (COP 0.0106)", (E_tape + E_lo) / COP / 1e12, (E_tape + E_l) / COP / 1e12),
            ("its share of the launch's 4.6 PJ, percent", (E_tape + E_lo) / COP / E_LAUNCH * 100, (E_tape + E_l) / COP / E_LAUNCH * 100),
            ("per meter of ring as a mean, kJ/m", (E_tape + E_lo) / ring.L_RING / 1e3, (E_tape + E_l) / ring.L_RING / 1e3),
            ("per 500 m site as a mean, MJ", (E_tape + E_lo) / ring.L_RING * 500 / 1e6, (E_tape + E_l) / ring.L_RING * 500 / 1e6),
            ("the mean load at 14 launches a day, W per site", (E_tape + E_lo) * 14 / 86400 / ring.L_RING * 500, (E_tape + E_l) * 14 / 86400 / ring.L_RING * 500)]
    for lab, a, b in rows:
        print(f"  {lab:>44} {a:12.2f} {b:12.2f}")
    print(f"  The run is {li['x_run']/1e3:,.0f} km in {li['t_run']/60:.0f} minutes, {laps:.2f} laps of the ring.  The minutes the spare plant takes are in the sheet, on the run's own peak and launch-window loads.")
    out = dict(E_tape=E_tape, E_lead_old=E_lo, E_lead=E_l, x_run=li["x_run"], t_run=li["t_run"], laps=laps, per_m=(E_tape + E_l) / ring.L_RING, per_site=(E_tape + E_l) / ring.L_RING * 500,
               P_full=P_full, lead_s=li["lead_s"], lead_s_old=li_o["lead_s"], per_m_tape=E_tape / ring.L_RING, per_m_lead=E_l / ring.L_RING)
    r9.save("round9_heat.json", out)
    return out


# =============================================================================
# item 3: Table 7.5, the tolerances, the failures
# =============================================================================

def _setup():
    import bearing_round6a as a6
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    import compensation_tracks as ct
    A2, theva, sup, lift, tabs = a6.casing_setup()
    cp = r7.load("round7_coupling.json")["c_split"]
    S1 = A2["dm_top"]["F_s"]; lift1 = fb.F_LIFT_TRACK; load1 = fb.F_LOAD_TRACK

    def lane(u, prot=True, **kw):
        d = dict(L_sc=bd.L_SC, tau=a6.TAU, L_section=a6.L_FAST / a6.N_FAM, gap_dep=False, steady_loop_hz=0.0, T_spread=max(0.05, bd.L_SC / u), B_sat=1.2, n_coils=a6.N_FAM, L_coil=a6.L_FAST,
                 hd_width=2 * lift1 / bd.P12, hu_width=2 * load1 / bd.P12, cmd_ramp=fb.RAMP, steady_cap=2 * S1, couple=(cp, 2 * load1))
        if prot:
            d["prot"] = ct.prot(0.010, 1e-3, cap_up=2 * load1, cap_down=2 * lift1)
        d.update(kw)
        return pz.lane(**d)
    return pz, sp, ct, tabs, lane, dict(S1=S1, lift1=lift1, load1=load1, cp=cp)


def part_table75():
    hr("ITEM 3: TABLE 7.5, THE PLANNED PULL, THE FIELD AND THE CASING'S MOVEMENT")
    pz, sp, ct, tabs, lane, c = _setup()
    cases = r9.load("round9_cases_flux.json")
    print(f"  The pass of the 10 km sled with its {pz.M_SC0/1e3:,.0f} t spacecraft on {bd.L_SC:.0f} m, the pair of tracks; 25 m coils in 20 families, the 10 ms rise, the planned pull ramped over {fb.RAMP*1e3:.0f} ms,")
    print(f"  the support frozen at constant flux, the coupling {c['cp']*100:.1f} percent.  The sled's load on the casing per meter under the spacecraft is q = (sled + spacecraft per meter) x a_down(u):")
    print(f"  positive (down) below the speed at which the sled is in orbit, negative above it.  The planned pull is the feed-forward's: q times the share {pz.MATCH_RING:.3f} that leaves casing and")
    print(f"  cable moving together, spread over max(50 ms, L_sc / u), and thinned by the coil's own length (a 25 m coil is on for the window plus 25 m / u).  Its peak is read from a run with")
    print(f"  no gap protection, so that nothing but the planned pull is in it.  With the protection on, the 10 mm trigger adds its own pull when the casing has moved 10 mm; the 'pull' column")
    print(f"  of the case logs (bearing_design.part_F, F.1) is the largest pull of either kind on either face, max(pk_hd, pk_hu), which is why it is not the planned pull.")
    print(f"  {'u km/s':>7} {'q under the spacecraft kN/m':>27} {'window ms':>9} {'planned peak, pair kN/m':>23} {'track':>6} {'bo.PASS_PULLS':>13} {'largest pull, protection on':>27} {'case log F.1':>12} "
          f"{'even field T':>12} {'movement mm, protection on':>26} {'no protection':>13} {'smallest gap mm':>15}")
    rows = []
    old_log = {1: 395.9, 3: 164.9, 5: 126.8, 10: 120.1, 20: 246.0, 30: 393.0}
    for (uk, F_pp, kind_pp) in bo.PASS_PULLS:
        u = uk * 1e3
        res, hist = pz.point_pulsed([lane(u, True, label="pair"), lane(u, False, label="pair, no protection")], u, tabs)
        rp, rn = res
        planned = max(rn["pk_hd"], rn["pk_hu"]); kind = "lift" if rn["pk_hd"] > rn["pk_hu"] else "load"
        both = (rn["pk_hd"], rn["pk_hu"])
        w = 1.0; B = math.sqrt(2 * MU0 * (planned / 2) / w)
        qz = float(pz.q_zone(u, bd.L_SC)); T = max(0.05, bd.L_SC / u)
        rows.append(dict(u=u, q=qz, T=T, planned=planned, kind=kind, pk_hd_plan=both[0], pk_hu_plan=both[1], pass_pulls=F_pp * 1e3, peak_prot=max(rp["pk_hd"], rp["pk_hu"]), pk_hd=rp["pk_hd"], pk_hu=rp["pk_hu"],
                         B=B, exc=rp["max_exc_all"], exc_noprot=rn["max_exc_all"], gap=min(rp["min_gap_b"], rp["min_gap_t"]), log_old=old_log[uk] * 1e3, engaged=bool(rp["engaged"])))
        print(f"  {uk:7d} {qz/1e3:27.1f} {T*1e3:9.0f} {planned/1e3:23.1f} {kind:>6} {F_pp:13.0f} {max(rp['pk_hd'], rp['pk_hu'])/1e3:27.1f} {old_log[uk]:12.1f} {B:12.3f} {rp['max_exc_all']*1e3:26.2f} {rn['max_exc_all']*1e3:13.2f} "
              f"{min(rp['min_gap_b'], rp['min_gap_t'])*1e3:15.1f}")
    print(f"  (the even field: sqrt(2 mu0 F / w) with F half the pair's planned pull on a plate w = 1.0 m wide.  The other face's planned peak, the bare sled's part: "
          + "; ".join(f"{r['u']/1e3:.0f} km/s lift {r['pk_hd_plan']/1e3:.1f}, load {r['pk_hu_plan']/1e3:.1f}" for r in rows) + " kN/m.)")
    # the movements in the 114-case run for the same lanes
    cs = {r["label"]: r for r in cases["records"]}
    print(f"\n  The same passes in the 114-case run (`bearing_round9.py --part=cases`), the movement in mm: "
          + "; ".join(f"{uk} km/s {cs[f'pair']['exc']*1e3:.2f}" if False else "" for uk in ()) + "see the tolerances below.")
    r9.save("round9_table75.json", dict(rows=rows))
    return rows


def part_failures():
    hr("ITEM 3: THE TOLERANCES AND THE FAILURES THE CHAPTER QUOTES (the 114 cases, the support at constant flux in every model)")
    c9 = r9.load("round9_cases_flux.json"); rows = c9["rows"]
    try:
        c9r = r9.load("round9_cases_flux_roll1hz.json")
    except FileNotFoundError:
        c9r = None
    g_old = json.load(open(os.path.join(r7.OUT, "gap_trade_logs", "gap_trade_cases_100mm_flux.json"), encoding="utf-8"))
    g_no = json.load(open(os.path.join(r7.OUT, "gap_trade_logs", "gap_trade_cases_100mm.json"), encoding="utf-8")) if os.path.exists(os.path.join(r7.OUT, "gap_trade_logs", "gap_trade_cases_100mm.json")) else None
    def mg(r): return min(r["gap_b"], r["gap_t"]) * 1e3
    def mv(r): return 100.0 - mg(r)
    byl = lambda d: {r["label"]: r for r in d["rows"]}
    A = byl(c9); B = byl(g_old); C = byl(g_no) if g_no else {}; D = byl(c9r) if c9r else {}
    print(f"  {len(rows)} cases.  Columns: the smallest gap in mm in (1) rounds 4 to 7 (the support's pull independent of the gap), (2) the gap trade (constant flux in the point lanes only),")
    print(f"  (3) round 9 (constant flux in the point, beam and roll lanes, the roll model's support frozen), (4) round 9 with the roll model's own 1 Hz support loop left on.")
    print(f"  {'case':>76} {'model':>6} {'(1)':>7} {'(2)':>7} {'(3)':>7} {'(4)':>7} {'movement (3) mm':>15} {'single failure':>14}")
    for r in rows:
        l = r["label"]
        print(f"  {l[-76:]:>76} {r['model']:>6} {mg(C[l]) if l in C else float('nan'):7.1f} {mg(B[l]) if l in B else float('nan'):7.1f} {mg(r):7.1f} {mg(D[l]) if l in D else float('nan'):7.1f} {mv(r):15.1f} {'yes' if r['single'] else 'no':>14}"
              + ("   CONTACT" if r["contact"] else ""))
    singles = [r for r in rows if r["single"]]
    worst = sorted(singles, key=mg)[:6]
    print(f"\n  The smallest gaps in any single failure: " + "; ".join(f"{r['label']} {mg(r):.1f} mm" for r in worst))
    print(f"  Contacts among the single failures: {[r['label'] for r in singles if r['contact']] or 'none'}.  All contacts: {[r['label'] for r in rows if r['contact']]}")
    moved = sorted(((abs(mg(r) - mg(B[r['label']])), r["label"], mg(B[r['label']]), mg(r)) for r in rows if r["label"] in B), reverse=True)[:12]
    print(f"  The cases the beam and roll terms move most against the gap trade's run: " + "; ".join(f"{l} {a:.1f} -> {b:.1f}" for d, l, a, b in moved if d > 0.05))
    # the tolerances
    print(f"\n  The tolerances (the pair), the largest movement in mm over 1, 3, 10 and 30 km/s:")
    tol = {}
    for r in rows:
        if r["label"].startswith("tolerance "):
            key = r["label"][10:].rsplit(" at ", 1)[0]
            tol.setdefault(key, []).append((mv(r), r["label"].rsplit(" at ", 1)[1]))
    kinds = [("timing 5 ms", ("timing -5 ms", "timing +5 ms")), ("mass 10 percent", ("mass -10 %", "mass +10 %")), ("speed 5 percent", ("speed -5 %", "speed +5 %")),
             ("gap measurement 5 mm", ("sensor -5 mm", "sensor +5 mm")), ("the rise halved", ("rise -50 %",)), ("the rise half again as long", ("rise +50 %",))]
    tout = []
    for name, keys in kinds:
        best = max((x for k in keys for x in [(m, f"{k} at {at}") for m, at in tol.get(k, [])]), default=(float("nan"), ""))
        tout.append(dict(name=name, mm=best[0], case=best[1]))
        print(f"    {name:>28}: {best[0]:6.2f} mm ({best[1]})")
    r9.save("round9_failures.json", dict(rows=[dict(r, gap=mg(r), move=mv(r), gap_r47=mg(C[r['label']]) if r['label'] in C else None, gap_gt=mg(B[r['label']]) if r['label'] in B else None,
                                                      gap_roll1=mg(D[r['label']]) if r['label'] in D else None) for r in rows], tolerances=tout,
                                            worst=[dict(label=r["label"], gap=mg(r)) for r in worst]))
