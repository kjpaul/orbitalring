#!/usr/bin/env python3
"""
The report of the gap trade round: Vol_III/gap_trade_2026-11-16.md, written from the logs (every number from a run).
Brief: Vol_III/Claude_Code_instructions_gap_trade_2026-11-16.md.
Run (from Python/orbitalring): python tools/gap_trade_report.py [commit]
Reads runs/gap_trade_logs/ (bearing_gap_trade.py and its parts), runs/round7_logs/, runs/round8_logs/.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

import ring_altitude as ring
import bearing_electrical as be
import bearing_gap_trade as gt
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_round8 as r8

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_MD = os.path.join(HERE, "..", "..", "..", "Vol_III", "gap_trade_2026-11-16.md")
H = 0.005
GAPS = gt.GAPS
CASES4 = [("pull 50 ms late", "the pull 50 ms late"), ("spread pull fired with no sled", "the planned pull fired with no sled"), ("sled passage 1 km/s, one track", "one load track alone at 1 km/s"),
          ("false full load pull 50 ms, no sled", "a false load pull of 50 ms"), ("3. one lift track lost at 30 km/s (the other alone, capped)", "one lift track lost at 30 km/s")]


def pct(a, b):
    return "" if abs(a / b - 1) < 5e-4 else f"{(a/b-1)*100:+.0f} %"


def mass_rule(tape_sup, tape_fast, W, V, I, ld):
    """Round 8's rule (bearing_round8_parts, item 8): tape; formers, cryostats and supports at one tape mass; film capacitors at 1.5 kJ/kg for the four
    fast tracks' stored energy; IGBT modules of 4.3 kg, V / 2.8 kV in series, two switches and two diodes per coil, 3.2 coils per meter; the fast leads."""
    tape = 2 * tape_sup + 4 * tape_fast
    store = 4 * W / 1.5e3
    n_s = math.ceil(V / 2800.0); sw = 4 * n_s * 4.3 * 3.2
    leads = 2 * 8960.0 * ld["A"] * (I / ld["I"]) * 1.0 * 3.2
    return dict(tape=tape, formers=tape, store=store, n_s=n_s, switches=sw, leads=leads, total=2 * tape + store + sw + leads)


def main():
    commit = sys.argv[1] if len(sys.argv) > 1 else "COMMIT_PLACEHOLDER"
    A = gt.load("gap_trade_a.json"); B = gt.load("gap_trade_b.json") or {}; SH = gt.load("gap_trade_shape.json"); FACE = gt.load("gap_trade_face.json")
    ST = gt.load("gap_trade_stiff.json"); PU = gt.load("gap_trade_pump.json"); SA = gt.load("gap_trade_start.json")
    half = r7.load("round7_half_2p5.json"); fast7 = r7.load("round7_fast_track.json"); rows8 = r8.load("round8_rows.json")
    ld = dict(p7.lead(), I=fast7["I"])
    Fc = gt.F_cap()
    def sa(g): return A[gt.key("support", g, gt.SUP["Wl"], gt.SUP_BUNDLE[0], gt.SUP_BUNDLE[1], H)]
    def fa(g): return A[gt.key("fast", g, gt.FAST["Wl"], gt.FAST["bw"], gt.FAST["bh"], H)]
    RULED = "cells turned as ruled"; BEST = "each stack the better way"
    L = []
    w = L.append
    s100, f100 = sa(0.10), fa(0.10)
    m100 = mass_rule(s100["cap"]["tape"], f100["ways"][RULED]["tape"], f100["W"], f100["V"], f100["I"], ld)
    kl = {g: ST["local"][f"{g:.3f}"] for g in GAPS if f"{g:.3f}" in ST["local"]}
    kf = {g: dict(zip(("k_I", "k_phi"), gt.fit_stiffness(g, sa(g)["normal"]["F"]))) for g in GAPS}
    def efold(g): return math.sqrt(ring.M_LOAD_M / (2 * abs(kf[g]["k_phi"]))) * 1e3

    # ------------------------------------------------------------------ the cases
    def cases(g, flux):
        d = gt.load(f"gap_trade_cases_{g*1e3:.0f}mm{'_flux' if flux else ''}.json")
        if d is None:
            return None
        by = {r["label"]: r for r in d["rows"]}
        out = dict(n=len(d["rows"]), rows={})
        for lab, short in CASES4 + [("sled passage 1 km/s, pair", "pair 1"), ("sled passage 30 km/s, pair", "pair 30")]:
            r = by[lab]; left = min(r["gap_b"], r["gap_t"])
            out["rows"][lab] = dict(left=left, move=g - left, contact=r["contact"] or left <= 1e-6)
        sing = [r for r in d["rows"] if r["single"]]
        out["touch"] = [r["label"] for r in sing if r["contact"] or min(r["gap_b"], r["gap_t"]) <= 1e-6]
        ws = min(sing, key=lambda r: min(r["gap_b"], r["gap_t"])); out["worst"] = (ws["label"], min(ws["gap_b"], ws["gap_t"]))
        return out
    CF = {g: cases(g, True) for g in GAPS}; CN = {g: cases(g, False) for g in GAPS}
    base = CN.get(0.10)

    # ------------------------------------------------------------------ item 1 (b) picks
    pk = gt.picks(H)
    def full(v):
        k = gt.key(v["kind"], v["gap"], v["Wl"], v["bw"], v["bh"], v["h"])
        return B.get(k) if (B.get(k) or {}).get("full") else (A.get(k) if (A.get(k) or {}).get("full") else v)

    # ================================================================== header
    w("# The gap trade, the field at the plate, and the flux pump (2026-11-16)\n")
    w("Brief: `Vol_III/Claude_Code_instructions_gap_trade_2026-11-16.md`. Run from `localbrain`. No book chapter was edited and nothing was installed in `Vol_III/images/`. US spelling. "
      "Code in `Python/orbitalring`: `bearing_gap_trade.py` (the solves and the cases), `bearing_gap_trade_parts.py` (the flux pump and the start of a launch), `tools/gap_trade_figures.py` (the two figures), "
      "`tools/gap_trade_report.py` (this report, written from the logs). Logs in `Python/orbitalring/runs/gap_trade_logs/`, copied to `Vol_III/runs_coil_optimization_2026-10-13/gap_trade_2026-11-16/`. "
      "Figures in `Vol_III/figure_work_2026-11-16/`. Every number is from a run or a derivation shown; where no script gives a number, the text says so.\n")
    w("**The design** is the one of rounds 7 and 8, unchanged: two support tracks, each of two independent 0.4 m loops 50 mm apart at the center on a 1.2 m by 100 mm plate, bundles 150 by 67 mm, 12 mm tape; "
      "four fast tracks, each one 0.6 m loop with 204 by 50 mm bundles on a 1.0 m by 200 mm plate, 3 mm tape; 10 mm from the tape to the gap; the book's cable. "
      "The requirements are the same at every gap: a support track gives 47.5 kN/m in normal operation and 245.5 kN/m at its cap with the worst tape at 50 percent of its limit; a fast track gives 196.5 kN/m with the worst tape at 50 percent.\n")
    w(f"**The cells.** The gap trade, the grid of item 1 (b) and the bundle shapes of item 2 are solved at 5 mm cells, about half a minute a solve. The design's own row at 100 mm is in rounds 7 and 8 at 2.5 mm, and the two agree as follows: "
      f"the support's cap {s100['cap']['NI']/1e3:.1f} kA per leg at 5 mm against {half['cap']['NI']/1e3:.1f} at 2.5 mm, its tape {s100['cap']['tape']:.1f} against {half['cap']['tape']:.1f} kg/m per track; "
      f"the fast track {f100['NI']/1e3:.1f} against {fast7['NI_req']/1e3:.1f} kA and {f100['ways'][RULED]['tape']:.1f} against {fast7['tape']:.1f} kg/m. So the 5 mm tables read 1 to 4 percent low in tape. "
      "The ratios between gaps are what the book needs, and both rows of a ratio are at 5 mm. The field at the plate (item 3) is at 2.5 mm. "
      "Each solve is converged to within 1 percent of its required pull, so a tape figure carries about 1 percent of that as well.\n")

    # ================================================================== look at first
    def ratio(g, kind):
        return (sa(g)["cap"]["tape"] / s100["cap"]["tape"]) if kind == "s" else (fa(g)["ways"][RULED]["tape"] / f100["ways"][RULED]["tape"])
    def mass_g(g):
        f = fa(g); return mass_rule(sa(g)["cap"]["tape"], f["ways"][RULED]["tape"], f["W"], f["V"], f["I"], ld)
    w("## 1. What the book session should look at first\n")
    n = [0]
    def item(t):
        n[0] += 1; w(f"{n[0]}. {t}")
    m50, m75, m125, m150 = (mass_g(g) for g in (0.05, 0.075, 0.125, 0.15))
    item(f"**The figure the book has been missing.** With the tracks as designed and only the gap changed, the bearing's tape per meter of casing is {m50['tape']:.0f} kg at 50 mm, {m75['tape']:.0f} at 75 mm, {m100['tape']:.0f} at 100 mm, "
         f"{m125['tape']:.0f} at 125 mm and {m150['tape']:.0f} at 150 mm. Against 100 mm that is {(m50['tape']/m100['tape']-1)*100:+.0f}, {(m75['tape']/m100['tape']-1)*100:+.0f}, {(m125['tape']/m100['tape']-1)*100:+.0f} and {(m150['tape']/m100['tape']-1)*100:+.0f} percent. "
         f"Halving the gap takes {(1-sa(0.05)['cap']['NI']/s100['cap']['NI'])*100:.0f} percent off the support's ampere-turns and {(1-fa(0.05)['NI']/f100['NI'])*100:.0f} percent off the fast track's. It does not halve them, because the tape sits 10 mm behind the gap and the bundle is 50 to 67 mm deep: "
         f"at a 50 mm gap the support's tape is still {sa(0.05)['tape_near']*1e3:.0f} to {sa(0.05)['tape_far']*1e3:.0f} mm from the plate. "
         f"One caution for the book's sentence that a smaller gap would cut the current: it cuts the ampere-turns. The current in a support coil goes the other way if the coil keeps its ten layers, "
         f"{10*sa(0.05)['normal']['NI']/sa(0.05)['cap']['tapes_leg']:.0f} A in normal operation at 50 mm against {10*s100['normal']['NI']/s100['cap']['tapes_leg']:.0f} A at 100 mm, because each tape sits in a lower field and is allowed more. "
         f"The fast coil's current does fall, {fa(0.05)['I']:,.0f} A against {f100['I']:,.0f} A.")
    gains = []
    for g in GAPS:
        ps = pk.get(("support", g)); pf = pk.get(("fast", g))
        if ps and pf:
            gains.append((g, (2 * ps["cap"]["tape"] + 4 * pf["ways"][BEST]["tape"]) / (2 * sa(g)["cap"]["tape"] + 4 * fa(g)["ways"][BEST]["tape"]) - 1, ps, pf))
    if gains:
        item("**Choosing the loop width and the bundle again at each gap changes little.** By the optimization's rule (the least tape with a rating of 1.35 or more, the plates kept), the bearing's tape moves by "
             + ", ".join(f"{x*100:+.1f} percent at {g*1e3:.0f} mm" for g, x, ps, pf in gains) + " against the design's loops and bundles at the same gap. "
             + ("No point of the grid meets the rule at " + ", ".join(f"{g*1e3:.0f} mm" for g in GAPS if g not in [x[0] for x in gains]) + ". " if len(gains) < len(GAPS) else "")
             + "The gap is the lever; the loop and the bundle are fine adjustment. Table (a) is the one to quote.")
    item(f"**The field at the plate does not change with the gap.** At the cap the root mean square of the field across the support plate's face is {sa(0.05)['cap']['line']['rms_n']:.2f} T at 50 mm and {sa(0.15)['cap']['line']['rms_n']:.2f} T at 150 mm, "
         f"next to the even field of {sa(0.10)['cap']['line']['even']:.2f} T that the load alone sets. What the gap changes is the field at the tape: {sa(0.05)['cap']['B']:.2f} T at the worst tape at 50 mm, {sa(0.10)['cap']['B']:.2f} T at 100 mm, {sa(0.15)['cap']['B']:.2f} T at 150 mm. "
         f"The tape's limit falls with that field, so a wider gap costs twice: more ampere-turns, and less current in each tape.")
    if CF.get(0.10) and CF.get(0.075):
        t = [g for g in GAPS if CF.get(g) and CF[g]["touch"]]
        okg = [g for g in GAPS if CF.get(g) and not CF[g]["touch"]]
        item(f"**The other side of the trade.** The casing's movement in a failure does not get smaller when the gap does. With the support's own instability carried in the cases (the correction in the next item), a single failure makes contact at "
             + " and ".join(f"{g*1e3:.0f} mm" for g in t) + f" and does not at " + ", ".join(f"{g*1e3:.0f}" for g in okg) + f" mm. At 100 mm the worst single failure (the pull 50 ms late) leaves {CF[0.10]['rows']['pull 50 ms late']['left']*1e3:.1f} mm. "
             f"So the sentence for the book is: a 75 mm gap would save {abs(m75['tape']/m100['tape']-1)*100:.0f} percent of the tape, and the pull that arrives 50 ms late would then put the casing on the cable.")
        if base:
            item(f"**A correction to the figure the brief quotes.** The 114 cases of rounds 4 to 7 held the frozen support's pull at the casing's weight whatever the gap. A support coil that is frozen keeps its flux, and at constant flux its pull grows as the gap closes: "
                 f"{abs(kf[0.10]['k_phi'])/1e6:.2f} kN/m per millimeter per track at 100 mm, solved here. The cases never carried that. With it, the pull 50 ms late moves the casing {CF[0.10]['rows']['pull 50 ms late']['move']*1e3:.1f} mm, not {base['rows']['pull 50 ms late']['move']*1e3:.1f}, "
                 f"and leaves {CF[0.10]['rows']['pull 50 ms late']['left']*1e3:.1f} mm, not {base['rows']['pull 50 ms late']['left']*1e3:.1f}. The other large movers grow by a few millimeters each, and an ordinary pass hardly at all "
                 f"({CF[0.10]['rows']['sled passage 1 km/s, pair']['move']*1e3:.1f} mm at 1 km/s and {CF[0.10]['rows']['sled passage 30 km/s, pair']['move']*1e3:.1f} mm at 30 km/s, against {base['rows']['sled passage 1 km/s, pair']['move']*1e3:.1f} and {base['rows']['sled passage 30 km/s, pair']['move']*1e3:.1f}). "
                 f"No single failure touches at 100 mm in either model. The book's 24.3 mm should become {CF[0.10]['rows']['pull 50 ms late']['left']*1e3:.1f} mm. For Paul to rule on.")
    item("**The instability hardly depends on the gap.** At constant flux the support's stiffness is " + ", ".join(f"{abs(kf[g]['k_phi'])/1e6:.2f}" for g in GAPS) + " kN/m per millimeter per track at " + ", ".join(f"{g*1e3:.0f}" for g in GAPS)
         + f" mm, and the e-folding time with the support frozen is " + ", ".join(f"{efold(g):.0f}" for g in GAPS) + " ms. The 135 ms of the earlier rounds came from 0.33 kN/m per millimeter, the book's track on the book's plate; "
         f"the split support on its own plate gives {abs(kf[0.10]['k_phi'])/1e6:.2f} and {efold(0.10):.0f} ms at 100 mm. A smaller gap does not make the bearing more unstable for the same load. It only leaves less room.")
    shp = {(v["bw"], v["bh"]): v for v in SH.values()}
    sq = shp[(0.100, 0.100)]; de = shp[(0.150, 0.0667)]; fl = shp[(0.200, 0.050)]; vf = shp[(0.250, 0.040)]
    item(f"**The bundle's height is worth little.** At the 100 mm gap the support needs {sq['cap']['tape']:.1f} kg/m of tape per track with a 100 by 100 mm bundle, {de['cap']['tape']:.1f} with the design's 150 by 67, {fl['cap']['tape']:.1f} with 200 by 50 and {vf['cap']['tape']:.1f} with 250 by 40. "
         f"The flat bundle keeps its tape closer, but it also puts the two legs of a loop closer to each other, and it needs more ampere-turns ({fl['cap']['NI']/1e3:.0f} and {vf['cap']['NI']/1e3:.0f} kA against {de['cap']['NI']/1e3:.0f}). "
         f"The book can say that a flat bundle keeps its tape closer; it should not promise a saving. The design's shape is within {abs(fl['cap']['tape']/de['cap']['tape']-1)*100:.0f} percent of the best of the four.")
    fn, fc_, fr = FACE["support normal"], FACE["support cap"], FACE["fast 196.5"]
    item(f"**The three fields of the brief hold.** Solved at 2.5 mm, the root mean square across the face is {fn['line']['rms_n']:.3f} T in normal operation (even field {fn['line']['even']:.3f} T), {fc_['line']['rms_n']:.3f} T at the cap ({fc_['line']['even']:.3f} T) and "
         f"{fr['line']['rms_n']:.3f} T for a fast track at 196.5 kN/m ({fr['line']['even']:.3f} T). The peak on the face is {fn['line']['peak']:.2f}, {fc_['line']['peak']:.2f} and {fr['line']['peak']:.2f} T. "
         f"At the cap the root mean square is {(fc_['line']['rms_n']/fc_['line']['even']-1)*100:.0f} percent over the even field because the plate is {fc_['sat']*100:.0f} percent saturated and part of the field then runs along the face, where it does not pull. The figure is drawn.")
    item(f"**The flux pump: the session's arithmetic is right, its drift is not.** 1 mV moves {PU['check'][0]['value']:.2f} A a day in one loop ({PU['check'][1]['value']:.2f} with both loops of a track moving), and 371 A at 1 mV takes {PU['check'][4]['value']:.0f} days ({PU['check'][5]['value']:.0f}). "
         f"But the deployment's drift is not 0.05 A a day. The casing starts in orbit with no weight on the bearing, and the whole {max(r['w'] for r in PU['series'])/1e3:.1f} kN/m comes on while it is slowed. "
         f"The weight changes by up to {PU['dw_max']:.0f} N/m a day (day {PU['dw_max_day']:.0f}), not 12.7, and the coil current by {PU['dI_day1']:.0f} A a day on day 1, {PU['series'][3]['dI']:.1f} on day 10, {PU['series'][8]['dI']:.1f} on day 100 and {abs(PU['dI_end']):.2f} at the end. "
         f"Of the pumps published at hundreds of amperes, the 0.27 mV one keeps up from day {PU['pumps'][2]['d_rate']} and the 0.05 mV one from day {PU['pumps'][0]['d_rate']}. Before that the supply carries the current through the leads. This fits Paul's ruling, and the book should say the supply does the deployment's rise, not only the first charge.")
    item(f"**The regulation heat.** With the support trimmed continuously the fast coils' regulation heat is {PU['j_reg']:,.0f} W per site in place of {PU['joule']:,.0f} W, the extra cold load {PU['rows4'][1]['extra']:,.0f} W in place of {PU['extra']:,.0f} W, "
         f"and the 368.9 days would become about {PU['rows4'][1]['days']:.1f} ({PU['rows4'][1]['d']:+.1f}). The deployment was not rerun. Two things beside it: the {PU['joule']:,.0f} W was computed for weekly trims with 12.7 N/m a day, and with the run's real drift the same rule gives {PU['weekly_mean']/1e3:.0f} kW, "
         f"so weekly trims never worked; and while the supply holds the current its leads bring {PU['sup_leads']['site_normal']:,.0f} W per site into the cold, which round 8's run does not have.")
    item(f"**A thing to check in the earlier rounds.** The coil's on-time in a pass is computed as the energized length over the sled's speed plus the cable's speed, which is the dwell time of a point of the plate. The coils are on the casing. "
         f"By the casing's frame the leads' full-current coil-seconds of the heaviest launch are {SA['lead_s_casing']/SA['lead_s_old']:.2f} times the figure in use, so the leads' {rows8['cold_launch']['E_lead']/1e9:.1f} GJ of round 8 would be about {rows8['cold_launch']['E_lead']/1e9*SA['lead_s_casing']/SA['lead_s_old']:.0f} GJ. Nothing was changed. For Paul to rule on.")
    item(f"**At 150 mm the design stops fitting.** The support's tape fills {sa(0.15)['fill']*100:.0f} percent of its 150 by 67 mm bundle (the packing of Part 1b allows {be.F_CABLE*(1-be.F_COOL)*100:.0f}), and its rating is {sa(0.15)['rating']['margin']:.2f} times the requirement, under the 1.35 ruled in round 5.")
    w("")

    # ================================================================== item 1 (a)
    w("## 2. Item 1: the gap trade\n")
    w("### 2.1 Table (a): the loops and bundles as designed, only the gap changed\n")
    w("Script: `bearing_gap_trade.py --part=solve-a` (5 mm cells; `gap_trade_a.json`). The 2.5 mm row is rounds 7 and 8 (`bearing_round7.py --part=solve-half`, `--part=solve-fast`).\n")
    w("**One support track at its cap, 245.5 kN/m** (two loops, four legs):\n")
    w("| Gap | Ampere-turns per leg | Field and angle at the worst tape | Its limit at 70 K | Tapes per leg | Tape per meter of track | Against 100 mm | Coil current | Plate saturated | Peak field in the plate | Rating over requirement | Tape fills the bundle |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for g in GAPS:
        c = sa(g)["cap"]
        w(f"| {g*1e3:.0f} mm | {c['NI']/1e3:.1f} kA | {c['B']:.2f} T, {c['th']:.0f}° from its normal | {c['ic']:.1f} A | {c['tapes_leg']:,.0f} | {c['tape']:.1f} kg | {pct(c['tape'], s100['cap']['tape'])} | {10*0.5*c['ic']:,.0f} A | {c['sat']*100:.0f} % | {c['B_plate']:.2f} T | {sa(g)['rating']['margin']:.2f} | {sa(g)['fill']*100:.0f} % |")
    w(f"| 100 mm at 2.5 mm cells | {half['cap']['NI']/1e3:.1f} kA | {half['cap']['B']:.2f} T, {half['cap']['th']:.0f}° | {half['cap']['ic']:.1f} A | {half['cap']['tapes_leg']:,.0f} | {half['cap']['tape']:.1f} kg | | {10*0.5*half['cap']['ic']:,.0f} A | {half['cap']['sat']*100:.0f} % | {half['cap']['B_plate']:.2f} T | {rows8['rating_support']['margin']:.2f} | |")
    w("")
    w("The coil current is ten layers at half the worst tape's limit. The tape count is set at the cap. \"Plate saturated\" is the share of the plate's section whose permeability has fallen below a tenth of its unsaturated value.\n")
    w("**The same track in normal operation, 47.5 kN/m:**\n")
    w("| Gap | Ampere-turns per leg | Coil current | Worst tape, share of its limit | Plate saturated | Peak field in the plate | Pull under the square law | Nearest and farthest tape from the plate |")
    w("|---|---|---|---|---|---|---|---|")
    for g in GAPS:
        s = sa(g); nn = s["normal"]; T = s["cap"]["tapes_leg"]
        w(f"| {g*1e3:.0f} mm | {nn['NI']/1e3:.1f} kA | {10*nn['NI']/T:.0f} A | {nn['share']*100:.1f} % | {nn['sat']*100:.1f} % | {nn['B_plate']:.2f} T | {(1-nn['F']/(s['k_lin']*nn['NI']**2))*100:.0f} % | {s['tape_near']*1e3:.0f} and {s['tape_far']*1e3:.0f} mm |")
    w(f"| 100 mm at 2.5 mm cells | {half['normal']['NI']/1e3:.1f} kA | {10*half['normal']['NI']/half['cap']['tapes_leg']:.0f} A | {half['normal']['share']*100:.1f} % | {half['normal']['sat']*100:.1f} % | | {(1-rows8['support_points'][0]['dev'])*100:.0f} % | 110 and 177 mm |")
    w("")
    w("**One fast track at 196.5 kN/m** (one loop, two legs; the cells turned as ruled in round 6b):\n")
    w("| Gap | Ampere-turns per leg | Field and angle at the worst tape | Its limit | Tapes per leg | Tape per meter of track | Against 100 mm | Coil current | Plate saturated | Peak field in the plate | Energy of a 25 m coil | Voltage for the 10 ms rise | Tape's pulse heat, heaviest launch | Rating over requirement |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for g in GAPS:
        f = fa(g); d = f["ways"][RULED]
        w(f"| {g*1e3:.0f} mm | {f['NI']/1e3:.1f} kA | {d['B']:.2f} T, {d['th']:.0f}° | {d['ic']:.1f} A | {d['tapes_leg']:,.0f} | {d['tape']:.1f} kg | {pct(d['tape'], f100['ways'][RULED]['tape'])} | {f['I']:,.0f} A | {f['sat']*100:.1f} % | {f['B_plate']:.2f} T | {f['E_coil']/1e3:.0f} kJ | {f['V']/1e3:.1f} kV | {(d['Ep']+d['Es'])/1e9:.0f} GJ | {d['margin']:.2f} |")
    w(f"| 100 mm at 2.5 mm cells | {fast7['NI_req']/1e3:.1f} kA | {fast7['B']:.2f} T, {fast7['th']:.0f}° | {fast7['ic']:.1f} A | {fast7['tapes']:,.0f} | {fast7['tape']:.1f} kg | | {fast7['I']:,.0f} A | {fast7['sat']*100:.1f} % | {fast7['B_plate']:.2f} T | {fast7['E_coil']/1e3:.0f} kJ | {fast7['V']/1e3:.1f} kV | {(fast7['Ep']+fast7['Es'])/1e9:.0f} GJ | {fast7['margin']:.2f} |")
    w("")
    w("The coil current is the track's ampere-turns over 20 overlapping coils of 10 turns. The peak field in the fast plate is in a few cells above the bundles; the plate as a whole is not saturated at any gap. "
      "The pulse heat is round 6b's model (`bearing_round6b.tape_loss_launch`) on the 5 mm field, so it is to be read as a ratio between gaps.\n")
    w("**The bearing on the casing, by round 8's rule** (two support and four fast tracks; formers, cryostats and supports at one tape mass; film capacitors at 1.5 kJ/kg; 4.3 kg modules, one for each 2.8 kV):\n")
    w("| Gap | Tape per meter of casing | Against 100 mm | Formers and cryostats | Stores | Switches | Fast leads | Bearing mass on the casing | Against 100 mm |")
    w("|---|---|---|---|---|---|---|---|---|")
    for g in GAPS:
        m = mass_g(g)
        w(f"| {g*1e3:.0f} mm | {m['tape']:.0f} kg | {pct(m['tape'], m100['tape'])} | {m['formers']:.0f} kg | {m['store']:.0f} kg | {m['switches']:.0f} kg ({m['n_s']} modules in series) | {m['leads']:.1f} kg | {m['total']:,.0f} kg | {pct(m['total'], m100['total'])} |")
    m25 = rows8["mass"]
    w(f"| 100 mm at 2.5 mm cells (round 8) | {m25['tape']:.0f} kg | | {m25['formers']:.0f} kg | {m25['store']:.0f} kg | {m25['switches']:.0f} kg | {m25['leads_fast']:.1f} kg | {m25['total']:,.0f} kg | |")
    w("")
    w(f"The stores follow the fast track's stored energy, which falls with the gap as the tape does ({fa(0.05)['W']/1e3:.0f} kJ per meter of track at 50 mm, {fa(0.10)['W']/1e3:.0f} at 100 mm, {fa(0.15)['W']/1e3:.0f} at 150 mm). "
      "The plates are the same at every gap, so the cable's hardware does not change. The lost support loop, the stray field, the coupling to the load track and the coil ends were not solved again at the other gaps.\n")

    # ------------------------------------------------------------------ (b)
    w("### 2.2 Table (b): the loop width and the bundle chosen again at each gap\n")
    nb = len([v for v in B.values() if v.get("ok")])
    w(f"Script: `bearing_gap_trade.py --part=solve-b`, then `--part=solve-b2` for the chosen points (5 mm cells; `gap_trade_b.json`, {nb} points). **The rule.** Round 5 searched a grid of plate width, plate thickness and loop width and took the lightest point whose rating is at least 1.35 times the requirement; "
      "round 6b tried bundle shapes of one cross-section. Here the plates are the design's and stay fixed, so at each gap the grid is the loop width and the bundle shape, and the point taken is the one with the least tape among those that fit the plate and whose rating is at least 1.35 times the requirement. "
      "Support: loops of 0.25 to 0.45 m and bundles of 100 by 100, 150 by 67, 200 by 50 and 250 by 40 mm, the track no wider than its 1.2 m plate. Fast: loops of 0.4 to 0.8 m and bundles of 204 by 50, 150 by 68 and 300 by 34 mm, the track no wider than its 1.0 m plate, each stack laid the better of the two ways (round 6b's rule, since the ruled pattern of turned cells belongs to the 204 by 50 mm bundle). "
      "The bundle's cross-section is kept. No point of the grid changes the plate, the tape or the 10 mm.\n")
    if pk:
        w("| Gap | Support: loop, bundle | Cap, ampere-turns per leg | Tape per meter of track | Against (a) at that gap | Rating | Fast: loop, bundle | Ampere-turns per leg | Tape per meter of track | Against (a) at that gap | Rating | Tape per meter of casing | Against (a) |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for g in GAPS:
            ps = pk.get(("support", g)); pf = pk.get(("fast", g))
            if not ps or not pf:
                P = [v for v in list(B.values()) + list(A.values()) if v.get("ok") and v["kind"] == "support" and abs(v["gap"] - g) < 1e-9]
                lo = min(P, key=lambda v: v["cap"]["tape"]); hi = max(P, key=lambda v: v["rating"]["margin"])
                fcell = (f"{pf['Wl']:.1f} m, {pf['bw']*1e3:.0f} by {pf['bh']*1e3:.0f} mm | {pf['NI']/1e3:.1f} kA | {pf['ways'][BEST]['tape']:.1f} kg | {pct(pf['ways'][BEST]['tape'], fa(g)['ways'][BEST]['tape']) or 'the same point'} | {pf['ways'][BEST]['margin']:.2f}" if pf else " | | | | ")
                w(f"| {g*1e3:.0f} mm | no point reaches 1.35: the best rating is {hi['rating']['margin']:.2f}; the least tape is {lo['Wl']:.2f} m, {lo['bw']*1e3:.0f} by {lo['bh']*1e3:.0f} mm | {lo['cap']['NI']/1e3:.1f} kA | {lo['cap']['tape']:.1f} kg | {pct(lo['cap']['tape'], sa(g)['cap']['tape']) or 'the same point'} | {lo['rating']['margin']:.2f} | {fcell} | | |"); continue
            ta = sa(g)["cap"]["tape"]; tfa = fa(g)["ways"][BEST]["tape"]; tb = ps["cap"]["tape"]; tfb = pf["ways"][BEST]["tape"]
            w(f"| {g*1e3:.0f} mm | {ps['Wl']:.2f} m, {ps['bw']*1e3:.0f} by {ps['bh']*1e3:.0f} mm | {ps['cap']['NI']/1e3:.1f} kA | {tb:.1f} kg | {pct(tb, ta) or 'the same point'} | {ps['rating']['margin']:.2f} | {pf['Wl']:.1f} m, {pf['bw']*1e3:.0f} by {pf['bh']*1e3:.0f} mm | {pf['NI']/1e3:.1f} kA | {tfb:.1f} kg | {pct(tfb, tfa) or 'the same point'} | {pf['ways'][BEST]['margin']:.2f} | "
              f"{2*tb+4*tfb:.0f} kg | {pct(2*tb+4*tfb, 2*ta+4*tfa) or 'the same'} |")
        w("")
        # the grids
        for kind, shapes, wls, lab in (("support", gt.SUP_SHAPES, gt.SUP_WL, "Support track: tape per meter of track at the cap, kg (rating over requirement)"), ("fast", gt.FAST_SHAPES, gt.FAST_WL, "Fast track: tape per meter of track, kg, each stack the better way (rating over requirement)")):
            w(f"**{lab}:**\n")
            w("| Loop | Bundle | " + " | ".join(f"{g*1e3:.0f} mm" for g in GAPS) + " |")
            w("|---|---|" + "---|" * len(GAPS))
            for bw, bh in shapes:
                for Wl in wls:
                    cells = []
                    for g in GAPS:
                        v = B.get(gt.key(kind, g, Wl, bw, bh, H)) or A.get(gt.key(kind, g, Wl, bw, bh, H))
                        if not v or not v.get("ok"):
                            cells.append(""); continue
                        tp, mg = (v["cap"]["tape"], v["rating"]["margin"]) if kind == "support" else (v["ways"][BEST]["tape"], v["ways"][BEST]["margin"])
                        star = "**" if pk.get((kind, g)) is not None and pk[(kind, g)]["Wl"] == Wl and pk[(kind, g)]["bw"] == bw else ""
                        cells.append(f"{star}{tp:.1f} ({mg:.2f}){star}" if mg == mg else f"{tp:.1f} (rating above the solved range)")
                    if any(cells):
                        w(f"| {Wl:.2f} m | {bw*1e3:.0f} by {bh*1e3:.0f} mm | " + " | ".join(cells) + " |")
            w("")
        w("The chosen point of each gap is in bold. A support loop narrower than its bundle plus a clear window, or a track wider than its plate, is not in the grid. "
          "What the table says: choosing again moves the tape by a few percent at any gap; the gap itself moves it by tens of percent. The design's 0.4 m loops and 0.6 m loop were chosen at 100 mm on other grounds as well (the lost loop's trim, the coupling, the packing of the fast cells), and none of those was checked for the other points.\n")
    else:
        w("Not yet solved.\n")

    # ------------------------------------------------------------------ dynamics
    w("### 2.3 The other side of the trade: the casing's movement in a failure\n")
    w("**The stiffness.** Script: `bearing_gap_trade.py --part=stiff` (`gap_trade_stiff.json`). The support track as designed at each gap, at its solved normal ampere-turns, the plate saturating: four solves, the gap 10 mm either way and the ampere-turns 5 percent either way. "
      "From them the change of pull with the gap at constant current, and at constant flux (each loop keeps its flux, which is what a closed superconducting loop does). The pull grows as the gap closes in both, so both are unstable.\n")
    w("| Gap | At constant current, per track (four solves) | At constant flux, per track (four solves) | At constant current (the fitted sweep) | At constant flux (the fitted sweep) | At constant flux, the pair | Instability per millimeter, of the casing's weight | E-folding time, support frozen |")
    w("|---|---|---|---|---|---|---|---|")
    for g in GAPS:
        k = kl.get(g)
        w(f"| {g*1e3:.0f} mm | " + (f"{abs(k['k_I'])/1e6:.2f} kN/m per mm | {abs(k['k_phi'])/1e6:.2f} kN/m per mm" if k else " | ") + f" | {abs(kf[g]['k_I'])/1e6:.2f} kN/m per mm | {abs(kf[g]['k_phi'])/1e6:.2f} kN/m per mm | {2*abs(kf[g]['k_phi'])/1e6:.2f} kN/m per mm | {2*abs(kf[g]['k_phi'])*1e-3/(2*gt.W_HALF)*100:.2f} % | {efold(g):.0f} ms |")
    w("")
    w("Two estimates. The four solves are with the plate saturating, differenced over 20 mm. The fitted sweep is 27 solves of the linear plate at gaps from 10 to 270 mm, the logarithms of the pull and of the flux fitted with fourth-order curves, scaled to the solved normal pull. "
      "The four-solve values scatter by up to 30 percent: the raw sweep shows steps of about 3 percent in the pull where the plate's face crosses a cell edge of the nested grids, and a 20 mm difference catches them. The fit is the steadier figure and is the one used below. "
      "The e-folding time is the square root of the casing's 12,000 kg/m over the pair's stiffness. The earlier rounds' 135 ms used 0.33 kN/m per millimeter per track, solved in Part 1 for the book's track. "
      f"What can be said: at constant flux the stiffness is {min(abs(kf[g]['k_phi']) for g in GAPS)/1e6:.2f} to {max(abs(kf[g]['k_phi']) for g in GAPS)/1e6:.2f} kN/m per millimeter at every gap, because the pull at a given flux hardly depends on the gap. At constant current it falls a little as the gap opens.\n")
    w("**The cases.** Script: `bearing_gap_trade.py --part=cases --gap=... [--noflux]` (`gap_trade_cases_*.json` and `.txt`). These are round 7's 114 cases, run again with the nominal gap set to each value: the 12,000 kg/m casing, 25 m fast coils in 20 families, the 10 ms rise, the 10 mm protection in every lane, the support frozen, the split support's coupling. "
      "Two models. In the first the support's pull does not depend on the gap, as in every earlier round; the fast tracks' currents are set for the measured gap, so nothing in the model knows the nominal gap and the movement is the same at every gap. "
      "In the second the frozen support keeps its flux, and its pull follows the gap by the fitted sweep above (`sled_passage_pulsed.SUPPORT_OF_GAP`, new, off by default). The second is the truer one for a support that is frozen. "
      "The protection's band is 10 mm at every gap, and the coupling is the 100 mm one; neither was chosen again. The beam and roll cases do not carry the new term.\n")
    for flux, lab, C in ((False, "Model 1: the support's pull independent of the gap (rounds 4 to 7). Movement, and the gap left:", CN), (True, "Model 2: the frozen support at constant flux. Movement, and the gap left:", CF)):
        w(f"**{lab}**\n")
        w("| Gap | " + " | ".join(short.capitalize() for _, short in CASES4) + " | An ordinary pass, 1 and 30 km/s | Single failures in contact | Smallest gap in a single failure |")
        w("|---|" + "---|" * (len(CASES4) + 3))
        for g in GAPS:
            c = C.get(g)
            if c is None:
                if not flux and base:
                    cells = []
                    for lab_, _ in CASES4:
                        mv = base["rows"][lab_]["move"]; left = g - mv
                        cells.append(f"{mv*1e3:.1f} mm, {left*1e3:.1f} mm left" if left > 0 else f"contact ({mv*1e3:.1f} mm wanted)")
                    wm = max(base["rows"][l]["move"] for l, _ in CASES4)
                    w(f"| {g*1e3:.0f} mm (not run: the 100 mm movement) | " + " | ".join(cells) + f" | {base['rows']['sled passage 1 km/s, pair']['move']*1e3:.1f} and {base['rows']['sled passage 30 km/s, pair']['move']*1e3:.1f} mm | | {max(g - wm, 0)*1e3:.1f} mm |")
                else:
                    w(f"| {g*1e3:.0f} mm | not run | | | | | | |")
                continue
            cells = []
            for lab_, _ in CASES4:
                r = c["rows"][lab_]
                cells.append("**contact**" if r["contact"] else f"{r['move']*1e3:.1f} mm, {r['left']*1e3:.1f} mm left")
            w(f"| {g*1e3:.0f} mm | " + " | ".join(cells) + f" | {c['rows']['sled passage 1 km/s, pair']['move']*1e3:.1f} and {c['rows']['sled passage 30 km/s, pair']['move']*1e3:.1f} mm | {len(c['touch'])}" + (f": {'; '.join(c['touch'][:4])}" + (" and more" if len(c['touch']) > 4 else "") if c['touch'] else "") + f" | {max(c['worst'][1], 0)*1e3:.1f} mm ({c['worst'][0]}) |")
        w("")
    if CF.get(0.10):
        t = [g for g in GAPS if CF.get(g) and CF[g]["touch"]]
        w("**Plainly:** with the frozen support at constant flux, a single failure makes contact at " + (" and ".join(f"{g*1e3:.0f} mm" for g in t) if t else "none of the gaps run") + ". "
          + (f"By model 1, in which the 75.7 mm of the brief was computed, contact comes at every gap of 75.7 mm or less, so at 50 and 75 mm. " if base else "")
          + "At 100 mm and above no single failure touches in either model. In model 2 the movement at the larger gaps is not the same as at 100 mm, because the constant-flux pull is taken along the real curve of each gap. "
          "A smaller gap could be paired with a faster protection or a narrower band; that was not tried, and the 50 ms late pull is a fault of the planned pull's timing, which a narrower band does not shorten.\n")

    # ================================================================== item 2
    w("## 3. Item 2: the height of the bundle\n")
    w("Script: `bearing_gap_trade.py --part=solve-shape` (5 mm cells; `gap_trade_shape.json`). The support track at the 100 mm gap with 0.4 m loops, the center bundles 50 mm apart, the same 0.01 m2 of bundle in four shapes.\n")
    w("| Bundle | Nearest and farthest tape from the plate | Cap, ampere-turns per leg | Field and angle at the worst tape | Its limit | Tapes per leg | Tape per meter of track | Against the design | Coil current at the cap | Normal operation | Plate saturated at the cap | Peak in the plate | Rating | Track width |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for bw, bh in gt.SUP_SHAPES:
        v = shp[(bw, bh)]; c = v["cap"]
        w(f"| {bw*1e3:.0f} by {bh*1e3:.0f} mm{' (the design)' if abs(bw-0.15) < 1e-9 else ''} | {v['tape_near']*1e3:.0f} and {v['tape_far']*1e3:.0f} mm | {c['NI']/1e3:.1f} kA | {c['B']:.2f} T, {c['th']:.0f}° | {c['ic']:.1f} A | {c['tapes_leg']:,.0f} | {c['tape']:.1f} kg | {pct(c['tape'], de['cap']['tape'])} | "
          f"{10*0.5*c['ic']:,.0f} A | {v['normal']['NI']/1e3:.1f} kA, {10*v['normal']['NI']/c['tapes_leg']:.0f} A | {c['sat']*100:.0f} % | {c['B_plate']:.2f} T | {v['rating']['margin']:.2f} | {v['width']:.2f} m |")
    w("")
    w(f"The 200 by 50 and 250 by 40 mm bundles make the track {fl['width']:.2f} and {vf['width']:.2f} m wide, wider than the 1.2 m plate; they are solved as asked, on the same plate. "
      f"Going from the square bundle to the design's saves {abs(de['cap']['tape']/sq['cap']['tape']-1)*100:.0f} percent of the tape. Going flatter than the design saves {abs(fl['cap']['tape']/de['cap']['tape']-1)*100:.0f} percent at 200 by 50 and costs {(vf['cap']['tape']/de['cap']['tape']-1)*100:.0f} percent at 250 by 40. "
      "The reason is in the columns: the flat bundle's farthest tape is closer, so the worst tape sits in a lower field and each tape carries more, but the wide legs of a 0.4 m loop come closer to each other and to the next loop, so more ampere-turns are needed for the same pull. The two nearly cancel.\n")

    # ================================================================== item 3
    w("## 4. Item 3: the field at the plate, in tesla\n")
    w("Script: `bearing_gap_trade.py --part=solve-face` (2.5 mm cells; `gap_trade_face.json`). The design at 100 mm. The field is read in the gap on a line 5 mm under the plate's face, across the plate's width. "
      "The even field with the same pull is the square root of 2 mu0 F / w with w the plate's width.\n")
    w("| State | Ampere-turns per leg | Pull | Even field, sqrt(2 mu0 F / w) | Root mean square across the face, the part through the face | Root mean square, the whole field | Peak on the face | Peak inside the plate | Plate saturated |")
    w("|---|---|---|---|---|---|---|---|---|")
    names = {"support normal": "Support track, normal operation", "support cap": "Support track, the cap", "support rating": "Support track, the rating", "support one loop": "Support track, one loop alone carrying 47.5 kN/m",
             "fast 196.5": "Fast track, 196.5 kN/m", "fast rating": "Fast track, the rating"}
    for k in gt.FACE_STATES:
        s = FACE[k]; l = s["line"]
        w(f"| {names[k]} | {s['NI']/1e3:.1f} kA | {s['F']/1e3:.1f} kN/m | {l['even']:.3f} T | {l['rms_n']:.3f} T | {l['rms']:.3f} T | {l['peak']:.2f} T | {s['B_plate']:.2f} T | {s['sat']*100:.1f} % |")
    w("")
    w(f"**The check.** The brief's 0.32, 0.72 and 0.70 T are the even fields: {fn['line']['even']:.3f}, {fc_['line']['even']:.3f} and {fr['line']['even']:.3f} T. "
      f"The pull is the integral over the face of (Bn2 - Bt2) / 2 mu0, with Bn the field through the face and Bt the field along it. Where the plate is not saturated Bt is small and the root mean square of Bn is the even field: {fn['line']['rms_n']:.3f} against {fn['line']['even']:.3f} T in normal operation. "
      f"At the cap and the rating the saturated plate lets field run along its face, so the root mean square through the face is above the even field ({fc_['line']['rms_n']:.3f} against {fc_['line']['even']:.3f} T at the cap). "
      f"For the fast track it is a little under ({fr['line']['rms_n']:.3f} against {fr['line']['even']:.3f} T): the line under the plate carries {fr['line']['F_line']/1e3:.1f} of the {fr['F']/1e3:.1f} kN/m, and the rest comes on the plate's edges, where the field in the figure rises. "
      f"The peak inside the plate is in a few cells above the bundles and is not the field of the face.\n")
    w("**The figure**, a candidate for the book: `Vol_III/figure_work_2026-11-16/fig07_candidate_field_at_plate.svg` and `.png` (`tools/gap_trade_figures.py`). The field across the face of a support plate in normal operation and at the cap, and of a compensation plate at 196.5 kN/m, in tesla against position, "
      "the bundles drawn under the axis to the scale of the axis, the even field as a dashed horizontal line. 4.5 inches wide, the smallest text 8.5 points, the book's palette. I opened it.\n")

    # ================================================================== item 4
    w("## 5. Item 4: the flux pump for the trims\n")
    w("Script: `bearing_gap_trade.py --part=pump` (`gap_trade_pump.txt`, `gap_trade_pump.json`). The coil current against the pull is from round 8's 2.5 mm solves of the split support; the inductances from round 7.\n")
    w(f"### 5.1 The session's arithmetic\n")
    w("| Figure | The session | Checked | Note |")
    w("|---|---|---|---|")
    w(f"| One loop's inductance | 17.8 H | {PU['L_loop']:.2f} H | each of the ten coils of a loop, the ten moving together; {PU['L_loop']+PU['M_loop']:.2f} H when the loop beside it on the same plate moves too, which is the case in a trim |")
    w(f"| 1 mV moves | 4.9 A a day | {PU['check'][0]['value']:.2f} A a day | {PU['check'][1]['value']:.2f} with both loops of a track |")
    w(f"| 0.05 mV moves | 0.24 A a day | {PU['check'][2]['value']:.2f} A a day | {PU['check'][3]['value']:.2f} with both loops |")
    w(f"| The first charge of 371 A at 1 mV | 76 days | {PU['check'][4]['value']:.0f} days | {PU['check'][5]['value']:.0f} days with both loops |")
    w(f"| The deployment's drift | about 0.05 A a day | {PU['drift_127']:.3f} A a day for 12.7 N/m a day on the casing | but 12.7 N/m a day is not the deployment's drift: see 5.2 |")
    w("")
    w(f"At the normal point the coil current moves {PU['slope_A_per_N']*1e3:.2f} A for each kN/m of a track's pull, {PU['exponent']:.2f} percent for each percent. The square law would give 0.50.\n")
    w("### 5.2 The drift through the deployment\n")
    w(f"The 12.7 N/m a day of round 6a is the fall of gravity as the ring rises 173 km. It left out the larger thing: the casing starts in orbit, with no weight on the bearing, and its weight comes on as the LIMs slow it. "
      f"From the round 8 design run (`runs/round8_design/ch14/rising_radius_deployment_tension.csv`, the column `casing_net_weight_kN_m`), with each track carrying half:\n")
    w("| Day | Casing's weight on the bearing | Its change | Coil current | Its change | Voltage a coil's pump must give (21.5 H) |")
    w("|---|---|---|---|---|---|")
    for r in PU["series"]:
        w(f"| {r['day']} | {r['w']/1e3:.2f} kN/m | {r['dw']:+.0f} N/m a day | {r['I']:.0f} A | {r['dI']:+.3f} A a day | {abs(r['V2'])*1e6:,.0f} µV |")
    w("")
    w(f"The weight changes fastest on day {PU['dw_max_day']:.0f}, {PU['dw_max']:.0f} N/m a day. The current changes fastest at the very start, where it goes as the square root of the weight. Over the whole deployment the mean is {PU['I_end']/368.9:.2f} A a day. "
      f"At the end the weight is falling {abs(PU['dw_end']):.0f} N/m a day and the current {abs(PU['dI_end']):.3f} A a day: that last figure is the one the session's 0.05 A a day should be replaced by.\n")
    w("### 5.3 The pump\n")
    w("The voltages are the time-averaged outputs read in full for Stage A (`runs/bearing_layout_pumps.txt`, `bearing_layout.py --part=pumps`): 0.046 mV at 1,150 A (Geng et al. 2019, transformer-rectifier), 0.197 mV at 30 A (Bumby et al. 2016, dynamo), 0.267 mV at 320 A (Mallett et al. 2023, self-rectifier), "
      "1.2 mV at 58 A (Jiang et al., in Zhou et al. 2023, a dynamo at 600 Hz, its saturated voltage). One pump per coil, so one pump for each of the ten coils of each loop.\n")
    w("| Pump voltage | Moves | Keeps up with the deployment's rate from | Running alone from day 0, has caught up by | Largest shortfall on the way | Source |")
    w("|---|---|---|---|---|---|")
    for p in PU["pumps"]:
        w(f"| {p['V']*1e3:.3f} mV | {p['rate']:.2f} A a day | day {p['d_rate']} | {'day ' + str(p['d_catch']) if p['d_catch'] else 'never in the deployment'} | {max(p['short'], 0):.0f} A | {p['src']} |")
    w("")
    w(f"**What one pump has to give.** For the trims after the deployment: nothing standing. A trim of 1 percent ({0.01*PU['I_n']:.1f} A) takes {0.01*PU['I_n']*(PU['L_loop']+PU['M_loop'])/86400*1e6:.0f} µV if it is to be done in a day and {0.01*PU['I_n']*(PU['L_loop']+PU['M_loop'])/86400/7*1e6:.0f} µV in a week. "
      f"For the deployment from day 200 on: {abs(PU['series'][10]['V2'])*1e6:.0f} µV and less, which the 0.2 and 0.27 mV pumps give, and the 0.05 mV one from day {PU['pumps'][0]['d_rate']}. For its first hundred days: {abs(PU['series'][8]['V2'])*1e6:.0f} µV on day 100 rising to {abs(PU['series'][3]['V2'])*1e6:,.0f} µV on day 10, "
      f"which is above every pump read at the current it would have to carry. No pump read has given 1 mV at 371 A [the combination is UNVERIFIED]. So the supply through the leads carries the deployment's rise for about its first hundred to two hundred days, and the pumps take over when the rate has fallen to what they can give.\n")
    w(f"**The heat it puts into the cold.** {PU['pumps_per_site']:.0f} pumps per 500 m site (4 loops, 10 coils per kilometer). A pump's loss is a standing loss while it runs, not a share of its output (Stage A).\n")
    w("| Kind of pump | Loss per pump | Per site, every pump running | Around the ring | Source |")
    w("|---|---|---|---|---|")
    for hrow in PU["pump_heat"]:
        w(f"| {hrow['label'].split(':')[0].capitalize()} | {hrow['q']:.2f} W | {hrow['site']:.0f} W | {hrow['site']/500*ring.L_RING/1e6:.1f} MW | {hrow['label'].split(':', 1)[1].strip()} |")
    w("")
    w("Neither loss was measured on a pump carrying 371 A into 21.5 H [UNVERIFIED for this duty]. After the deployment the pumps run only for a trim, so their mean heat is far below the table.\n")
    w("### 5.4 The fast coils' regulation heat\n")
    w(f"Round 6a's rule: with ten families connected, the Joule heat in the fast leads goes as the mean pull a load track holds. That was 33 N/m for a 0.1 mm regulation plus 44 N/m for the weight a weekly trim leaves at 12.7 N/m a day, {PU['joule']:,.0f} W per site in all, of which the regulation is {PU['j_reg']:,.0f} W. "
      f"With the run's real drift the weekly trim's part is not 44 N/m: the same rule gives {PU['weekly_mean']/1e3:.1f} kW per site at the mean drift of {PU['dw_mean']:.0f} N/m a day and {PU['weekly_max']/1e3:.0f} kW at the largest. So weekly trims were never enough for the deployment, and round 8's {PU['joule']:,.0f} W is right only if the support follows the weight continuously, by the supply or by pumps.\n")
    w("| The support trimmed by | Regulation heat | Pumps | Fast leads' conduction | Extra cold load per site | Against round 8's 4,411 W | Days, at 2.0 per kW | Deployment |")
    w("|---|---|---|---|---|---|---|---|")
    for r in PU["rows4"]:
        w(f"| {r['label'].capitalize()} | {r['joule']:,.0f} W | {r['pumps']:.0f} W | {r['cond']:,.0f} W | {r['extra']:,.0f} W | {r['extra']-PU['extra']:+,.0f} W | {r['d']:+.1f} | {r['days']:.1f} days |")
    w("")
    w(f"The 2.0 days per kW of cold load per site is from the deployment runs of round 6a (8.8 days for 4.5 kW) and round 5 (12.3 days for 6.1 kW). **The deployment was not rerun.** The number for Paul: {PU['rows4'][1]['extra']:,.0f} W per site with dynamo pumps, about {abs(PU['rows4'][1]['d']):.1f} days less, {PU['rows4'][1]['days']:.1f} days. "
      f"One more term belongs beside it. While the supply carries the support's current its leads are attached: a lead sized for the cap brings {PU['sup_leads']['q0']:.1f} W into the cold with no current and {PU['sup_leads']['qn']:.1f} W at 371 A, "
      f"{PU['sup_leads']['site_normal']:,.0f} W per site for the {PU['sup_leads']['n']:.0f} leads, {PU['sup_leads']['site_normal']/1e3*2.0:.1f} days if they stayed attached for the whole deployment. Round 8's run has no support lead in its cold load. "
      f"If the leads are attached only for the first 100 to 200 days, the cost is in proportion, {PU['sup_leads']['site_normal']/1e3*2.0*100/368.9:.1f} to {PU['sup_leads']['site_normal']/1e3*2.0*200/368.9:.1f} days.\n")

    # ================================================================== item 5
    w("## 6. Item 5: for Volume III\n")
    w("### 6.1 The start of a launch\n")
    w(f"Script: `bearing_gap_trade.py --part=start` (`gap_trade_start.txt`). **This is arithmetic on round 8's numbers, not a run.** No script carries the stretch below 1 km/s. At 0.83 g the spacecraft reaches 1 km/s after {SA['t1']:.0f} s and {SA['x1']/1e3:.0f} km.\n")
    w("| Rule | Speed | Reached after | Distance |")
    w("|---|---|---|---|")
    w(f"| Below this the support can follow through its leads at 1 kV (normal to cap in {rows8['parked']['times']['1000']['T']:.1f} s, the load moving its own 100 m in that time) | {SA['u_sup']:.1f} m/s | {SA['u_sup']/SA['a']:.1f} s | {0.5*SA['u_sup']**2/SA['a']:.0f} m |")
    w(f"| The same at 210 V ({rows8['parked']['times']['210']['T']:.0f} s) | {SA['u_sup210']:.1f} m/s | {SA['u_sup210']/SA['a']:.1f} s | {0.5*SA['u_sup210']**2/SA['a']:.1f} m |")
    w(f"| Above this a 25 m load coil under the 100 m spacecraft is on for less than its leads' {SA['budget20']:.2f} s of full current | {SA['u_load']:.0f} m/s | {SA['u_load']/SA['a']:.1f} s | {0.5*SA['u_load']**2/SA['a']:.0f} m |")
    w(f"| The spacecraft leaves the kilometer that its raised 1 km support coils cover | {SA['u_leave']:.0f} m/s | {SA['t_leave']:.0f} s | 1,000 m |")
    w("")
    w(f"**There is a band between the two.** From {SA['u_sup']:.0f} to {SA['u_load']:.0f} m/s, the first {SA['u_load']/SA['a']:.0f} seconds and {0.5*SA['u_load']**2/SA['a']:.0f} m, the support cannot be moved fast enough through its leads and a load coil would be on for longer than its leads allow. "
      f"What round 7's scheme does there by itself: the support coils raised for the parked spacecraft are a kilometer long, so for the first {SA['t_leave']:.0f} seconds the spacecraft is still over support that is at its cap, and the cap does not have to move at all. "
      "That is the same open item as in round 8: how a cap is confined to 100 m with 1 km coils in ten staggered families, and what the support does for the kilometer either side where there is no load. "
      f"After the first kilometer the spacecraft is at {SA['u_leave']:.0f} m/s, inside the load tracks' budget, and the load tracks carry it as they do at 1 km/s. "
      "So by this arithmetic the band closes if the raised kilometer is accepted, and the raised coils come down behind the spacecraft at the supply's own pace. "
      "A run of this stretch needs the support's coils resolved along the track, which no script has. [The conclusion is UNVERIFIED until that run exists.]\n")
    w("### 6.2 The coil figure\n")
    w("`Vol_III/figure_work_2026-11-16/fig07_candidate_2_two_coils.svg` and `.png` (`tools/gap_trade_figures.py`). The support track is drawn as two separate loops, each of ten coils of one width, one behind the other and offset along the track; the fast track as one loop of twenty coils drawn the same way. "
      "The figure-eight and its crossover are gone, and no coil is drawn narrower than another. The smallest text is 8.5 points at 4.5 inches. I opened it. For Volume III, not for the book.\n")

    # ================================================================== not run
    w("## 7. What was not run or not computed\n")
    for t in ["**The gap rows at 2.5 mm cells.** The trade is at 5 mm; only the 100 mm row exists at 2.5 mm. The comparison of the two at 100 mm is given at the top.",
              "**At the other gaps:** the lost support loop and its trim, the stray field, the load track's coupling to the support, the coil ends, the plate heat, the dump. The cases at the other gaps use the 100 mm coupling.",
              "**The protection chosen again for a smaller gap.** The band stays 10 mm and the loop 5 Hz at every gap.",
              "**The constant-flux term in the beam and roll cases.** It is in the point cases only, and it uses the linear plate's sweep of the pull against the gap.",
              "**Item 1 (b) beyond the tape and the rating.** The chosen points were not checked for the lost loop, the coupling or the packing of the fast cells, and the bundle's cross-section was not changed with the tape count.",
              "**The deployment.** Not rerun, as the brief says.",
              "**A pump for this duty.** No source read gives a pump's voltage or loss at 371 A into 21.5 H. The pump is not sized.",
              "**The start of a launch as a run.** Arithmetic only.",
              "**The coil's on-time in a pass.** The figure in use was flagged, not changed; the leads' section and the leads' pulse heat of round 8 stand until Paul rules.",
              "**No new literature search.** The pump figures are the ones read in full for Stage A."]:
        w(f"- {t}")
    w("")
    w(f"Tests: {sys.argv[2] if len(sys.argv) > 2 else 'TESTS_PLACEHOLDER'}. Committed and pushed: {commit}.")
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    txt = "\n".join(L) + "\n"
    assert "—" not in txt and "–" not in txt, "a dash crept in"
    open(OUT_MD, "w", encoding="utf-8").write(txt)
    print(f"wrote {os.path.abspath(OUT_MD)}: {len(L)} lines")


if __name__ == "__main__":
    main()
