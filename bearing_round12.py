#!/usr/bin/env python3
"""
Round 12 (2026-11-17), item 5: the gap trade's failure column on the regulator of rounds 10 and 11.
Book: Orbital Ring Engineering, Volume III: Chapter 7, Table 7.4.  Brief: Vol_III/Claude_Code_instructions_round12_2026-11-17.md.

    python bearing_round12.py [--jobs=5]

At each gap of the trade (50, 75, 100, 125, 150 mm) the worst single failures again: the pull 50 ms early and late and a pull fired with no sled
(the point and the cross-section models, at 30 and 1 km/s), and one 25 m section of a lift track out at 30 km/s and of a load track at 1 km/s
(the along-track model, at the stated torsional stiffness and a tenth of it, lost as the pull starts and known before the pass).

What changes with the gap, as in the gap trade's own cases (bearing_gap_trade.part_cases): the nominal gap, and the frozen support's pull at
constant flux against its gap (bearing_gap_trade.flux_table; at 100 mm the split support's solved table, the default since round 9).  The tracks
are taken to give the same pulls and ratings at every gap (the gap trade solved them again at each gap to do so).  The load coils' coupling to
the support is the 100 mm figure at every gap (not solved at the other gaps).

Writes runs/round12_logs/round12_gaps.json; the table goes to stdout.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

import bearing_round7 as r7

LOGS = os.path.join(r7.OUT, "round12_logs")
GAPS = (0.05, 0.075, 0.10, 0.125, 0.15)


def set_gap(gap):
    import sled_passage_pulsed as pz
    import bearing_gap_trade as gtr
    pz.GAP0 = gap
    pz.SUPPORT_OF_GAP = pz._support_of_gap_default() if abs(gap - 0.10) < 1e-9 else gtr.flux_table(gap)


def one(c):
    import bearing_round11 as r11
    set_gap(c["gap"])
    if c["model"] == "track":
        o = r11.one(c)
        return dict(case=c, gap_left=o["gap"], contact=o["contact"], t_contact=o["t_contact"], roll=o["roll"], move=o["move_c"])
    import bearing_round10 as r10
    S = r11.S_(); pz = S["pz"]; u = c["u"]
    res, _ = pz.point_pulsed([S["lane"](u, label="x", **c["kw"])], u, S["tabs"], t_settle=3.0)
    r = res[0]
    o, _, _ = r10.section(S, u, label="x", **c["kw"])
    return dict(case=c, gap_left=min(r["min_gap_b"], r["min_gap_t"]), contact=bool(r["contact"]), move=r["max_exc_all"], sec_gap=o["gap"], sec_contact=o["contact"],
                sec_t=o["t_contact"], roll=o["roll"])


def main():
    jobs = 5
    for a in sys.argv[1:]:
        if a.startswith("--jobs="): jobs = int(a[7:])
        elif a in ("-h", "--help"): print(__doc__); return
    t0 = time.time()
    P = (("nothing wrong", {}), ("the pull 50 ms early", dict(dt_err=-50e-3)), ("the pull 50 ms late", dict(dt_err=50e-3)), ("the pull fired with no sled", dict(no_sled=True)))
    cases = []
    for gap in GAPS:
        for lab, kw in P:
            for uk in (30, 1):
                cases.append(dict(gap=gap, model="point", u=uk * 1e3, kw=kw, label=f"{lab}, {uk} km/s"))
        for scale in (1.0, 0.1):
            for uk in (30, 1):
                nm = "lift" if uk == 30 else "load"
                cases.append(dict(gap=gap, model="track", u=uk * 1e3, unit="track", Lz=25.0, t_rel="start", known=False, scale=scale, L=2000.0 if uk == 1 else None,
                                  label=f"one 25 m section of a {nm} track lost as the pull starts, {uk} km/s, torsion x {scale:g}"))
                cases.append(dict(gap=gap, model="track", u=uk * 1e3, unit="track", Lz=25.0, t_rel=None, known=True, plan="moved", scale=scale, L=2000.0 if uk == 1 else None,
                                  label=f"one 25 m section of a {nm} track out before the pass, {uk} km/s, torsion x {scale:g}"))
    with Pool(jobs) as p:
        out = p.map(one, cases, chunksize=1)
    labs = []
    for c in cases:
        if c["label"] not in labs: labs.append(c["label"])
    R = {(o["case"]["label"], o["case"]["gap"]): o for o in out}
    print("\n" + "=" * 150 + "\n  THE GAP LEFT IN THE WORST SINGLE FAILURES AT EACH GAP, WITH THE REGULATOR (mm; C = contact)\n" + "=" * 150)
    print(f"  {'case':<86} {'model':>14} " + " ".join(f"{g*1e3:>8.0f}" for g in GAPS) + "  mm gap")
    for lab in labs:
        if R[(lab, GAPS[0])]["case"]["model"] == "point":
            print(f"  {lab:<86} {'point':>14} " + " ".join((f"{'C':>8}" if R[(lab, g)]["contact"] else f"{R[(lab, g)]['gap_left']*1e3:8.1f}") for g in GAPS))
            print(f"  {'':<86} {'cross-section':>14} " + " ".join((f"{'C':>8}" if R[(lab, g)]["sec_contact"] else f"{R[(lab, g)]['sec_gap']*1e3:8.1f}") for g in GAPS))
        else:
            print(f"  {lab:<86} {'along track':>14} " + " ".join((f"{'C':>8}" if R[(lab, g)]["contact"] else f"{R[(lab, g)]['gap_left']*1e3:8.1f}") for g in GAPS))
    worst = {}
    for g in GAPS:
        rows = [(min(o["gap_left"], o.get("sec_gap", 1.0)) if not (o["contact"] or o.get("sec_contact")) else 0.0, o["case"]["label"]) for o in out
                if o["case"]["gap"] == g and "torsion x 0.1" not in o["case"]["label"]]
        rows_s = [(o["gap_left"] if not o["contact"] else 0.0, o["case"]["label"]) for o in out if o["case"]["gap"] == g and "torsion x 0.1" in o["case"]["label"]]
        worst[g] = dict(stated=min(rows), soft=min(rows_s))
        print(f"  gap {g*1e3:.0f} mm: the smallest gap left {worst[g]['stated'][0]*1e3:.1f} mm ({worst[g]['stated'][1]}); with a tenth of the torsion the sections leave {worst[g]['soft'][0]*1e3:.1f} mm ({worst[g]['soft'][1]})")
    os.makedirs(LOGS, exist_ok=True)
    json.dump(dict(gaps=GAPS, rows=out, worst={str(k): v for k, v in worst.items()}), open(os.path.join(LOGS, "round12_gaps.json"), "w", encoding="utf-8"), indent=1, default=float)
    print(f"\n  DONE ({len(cases)} runs, {time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
