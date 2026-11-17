#!/usr/bin/env python3
"""
Round 13 (2026-11-17), item 7: an unplanned release of the spacecraft, as the ring's bearing sees it.
Book: Orbital Ring Engineering, Volume III: Chapter 10 (10.1, A Failure of the Track) and Chapter 7.
Brief: Vol_III/Claude_Code_instructions_round13_2026-11-17.md.

    python bearing_round13.py

The sled passes at 30 km/s and at the crewed limit (15.8 km/s over the casing).  The planned pull fires as planned, for the 5,000 t spacecraft.
The spacecraft's load is absent (a) for the whole pass, (b) from a point partway through.  The sled itself still passes, so its own load is
there; round 10's "pull fired with no sled" had no load at all and is run beside it.  The control law is round 10's regulator, unchanged: it
follows the path the plan predicts, and the plan still expects the spacecraft.

Three models, as in round 10: the point model (sled_passage_pulsed.point_pulsed), the cross-section with the roll carried (gap_regulator.py),
and the beam model (sled_passage_pulsed.beam_pulsed, the 10 m truss, a 10 km window).  In the point and cross-section models "partway" is the
spacecraft released when it is a quarter, half and three quarters over the place.  In the beam model the spacecraft is released with its middle
at the middle of the window, so the one run holds every place: behind it the full load, under it part, ahead of it the planned pull and no load.

Writes runs/round13_logs/round13_release.json; the table goes to stdout.
"""
import json
import os
import time

import bearing_round7 as r7
import bearing_round10 as r10

LOGS = os.path.join(r7.OUT, "round13_logs")
SPEEDS = (30e3, 15.8e3)


def main():
    import bearing_design as bd
    import bearing_round6a as a6
    import sled_passage_simulation as sp
    S = r10.setup(); pz = S["pz"]; L = bd.L_SC
    t0 = time.time(); rows = []
    print("=" * 150 + "\n  AN UNPLANNED RELEASE: THE PLANNED PULL FIRES AND THE SPACECRAFT'S LOAD IS NOT THERE (the smallest gap in mm; round 10's regulator, unchanged)\n" + "=" * 150)
    print(f"  {'case':<74} {'km/s':>5} {'planned pull, pair kN/m':>23} {'point':>8} {'cross-section':>14} {'roll mrad':>9} {'beam':>8}")
    for u in SPEEDS:
        cases = [("nothing wrong", {}, {}),
                 ("the pull fired with no sled at all (round 10's case)", dict(no_sled=True), dict(no_sled=True)),
                 ("(a) the spacecraft absent for the whole pass, the sled passes", dict(release_t=-1e9), dict(release_t=-1e9)),
                 ("(b) released as its front reaches the place", dict(release_t=0.0), None),
                 ("(b) released a quarter of the way over the place", dict(release_t=0.25 * L / u), None),
                 ("(b) released half way over the place", dict(release_t=0.5 * L / u), None),
                 ("(b) released three quarters of the way over the place", dict(release_t=0.75 * L / u), None),
                 ("(b) released in the middle of the window (beam: every place in one run)", None, dict(release_t=(5000.0 + L / 2) / u))]
        bc = [dict(label=lab, **kb) for lab, kp, kb in cases if kb is not None]
        resb, _ = pz.beam_pulsed(u, L, [sp.truss(10.0, 0.1)], bc, S["tabs"], tau=a6.TAU, L_win=10000.0, T_spread=max(0.05, L / u), n_coils=a6.N_FAM, L_coil=a6.L_FAST,
                                 track_width=2 * S["c"]["lift1"] / bd.P12, prot=S["reg_pair"](), loop_hz=0.0)
        B = {r["case"]: r for r in resb}
        for lab, kp, kb in cases:
            row = dict(case=lab, u=u)
            if kp is not None:
                res, h = pz.point_pulsed([S["lane"](u, label="x", **kp)], u, S["tabs"], t_settle=3.0)
                r = res[0]
                o, _, _ = r10.section(S, u, label="x", **kp)
                row.update(point=min(r["min_gap_b"], r["min_gap_t"]), point_contact=bool(r["contact"]), planned=float(abs(h["F_ff"][0]).max()), section=o["gap"], section_contact=o["contact"],
                           roll=o["roll"], move=r["max_exc_all"])
            if kb is not None:
                b = B[lab]; row.update(beam=min(b["min_gap_b"], b["min_gap_t"]), beam_contact=bool(b["contact"]))
            rows.append(row)
            f = lambda k: (f"{'C':>8}" if row.get(k + "_contact") else f"{row[k]*1e3:8.1f}") if k in row else f"{'':>8}"
            print(f"  {lab:<74} {u/1e3:5.1f} {row.get('planned', float('nan'))/1e3:23.1f} {f('point')} {f('section'):>14} {row.get('roll', float('nan'))*1e3:9.2f} {f('beam')}", flush=True)
    print(f"\n  The loads at a place: the spacecraft and sled under the spacecraft {pz.q_zone(30e3, L)/1e3:,.0f} kN/m at 30 km/s and {pz.q_zone(15.8e3, L)/1e3:,.0f} kN/m at 15.8 km/s; the bare sled "
          f"{pz.q_bare(30e3)/1e3:.1f} and {pz.q_bare(15.8e3)/1e3:.1f} kN/m (negative: the sled lifts the casing).")
    os.makedirs(LOGS, exist_ok=True)
    json.dump(dict(rows=rows), open(os.path.join(LOGS, "round13_release.json"), "w", encoding="utf-8"), indent=1, default=float)
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
