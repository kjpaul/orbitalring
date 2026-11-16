#!/usr/bin/env python3
"""
Round 9 (2026-11-16): the reruns.  Brief: Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md.
Book: Orbital Ring Engineering, Volume III: Chapters 7, 8 and 14.

Three corrections to the model are in the modules themselves and are the defaults from this round on:
  1. bearing_round4.hold_time / bearing_round5.hold_time: a fast coil's on-time is the energized length over the sled's speed (the coil is on the casing);
  2. sled_passage_pulsed.SUPPORT_OF_GAP: the support's pull at constant flux follows its gap, in the point, beam and roll models;
  3. bearing_cold_load.py through lim_config.BEARING_COLD_LOAD: the bearing's cold load in the deployment from the run's own weight series.

Run (from Python/orbitalring):
  python bearing_round9.py --part=cases [--roll-hz=1] [--noflux]    the 114 cases (the support frozen; --roll-hz=1 keeps the roll model's own 1 Hz support loop)
  python bearing_round9.py --part=ontime | table75 | failures | heat | stray | stator      the tables (bearing_round9_parts.py)
Logs: runs/round9_logs/.
"""
import json
import os
import sys
import time

import ring_altitude as ring
import bearing_round7 as r7

LOGS = os.path.join(r7.OUT, "round9_logs")


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    json.dump(obj, open(os.path.join(LOGS, name), "w", encoding="utf-8"), indent=1, default=float)


def load(name):
    return json.load(open(os.path.join(LOGS, name), encoding="utf-8"))


def part_cases(roll_hz=0.0, flux=True, tag=None):
    """Round 7's 114 cases (final_bearing's Part B through bearing_layout.run_cases: the 12,000 kg/m casing, 25 m fast coils in 20 families, the 10 ms
    rise, the 10 mm protection in every lane, the support frozen, the split support's coupling) with the support's pull at constant flux in every
    model (the default from this round) and the roll model's support frozen as the other lanes' is."""
    import bearing_layout as bl
    import bearing_round6a as a6
    import sled_passage_pulsed as pz
    import compensation_tracks as ct
    cp = r7.load("round7_coupling.json")
    tag = tag or ("flux" if flux else "noflux") + ("" if roll_hz == 0 else f"_roll{roll_hz:g}hz")
    old = pz.SUPPORT_OF_GAP, ct.ROLL_SUPPORT_HZ
    if not flux:
        pz.SUPPORT_OF_GAP = None
    ct.ROLL_SUPPORT_HZ = roll_hz
    try:
        A2, theva, sup, lift, tabs = a6.casing_setup()
        t0 = time.time()
        gaps, records, txt = bl.run_cases(a6.TAU, True, A2, theva, sup, lift, tabs, os.path.join("round9_logs", f"round9_cases_{tag}.txt"), prot_all=True,
                                          head=f"[round 9: {'the support at constant flux in the point, beam and roll models' if flux else 'the support pull independent of the gap'}, "
                                               f"the roll model's support loop {roll_hz:g} Hz, coupling {cp['c_split']*100:.2f} percent]", couple=cp["c_split"])
    finally:
        pz.SUPPORT_OF_GAP, ct.ROLL_SUPPORT_HZ = old
    rows = [dict(label=g[0], model=str(g[1]), gap_b=float(g[2]), gap_t=float(g[3]), contact=bool(g[4]), single=bool(bl.is_single(g[0]))) for g in gaps]
    s = bl.parse_summary(txt)
    recs = [dict(label=r[0], u=float(r[1]), exc=float(r[4]["exc"]), gap_b=float(r[4]["gap_b"]), gap_t=float(r[4]["gap_t"]), contact=bool(r[4]["contact"]),
                 S_min=float(r[3].min()), S_max=float(r[3].max())) for r in records]
    save(f"round9_cases_{tag}.json", dict(flux=flux, roll_hz=roll_hz, rows=rows, records=recs, summary=s, seconds=time.time() - t0))
    print(f"  {tag}: {len(rows)} cases, smallest single-failure gap {s['gap']} mm in '{s['case']}', contacts: {s['contacts']} ({time.time()-t0:.0f} s)", flush=True)
    return rows


def main():
    part = "cases"; roll_hz = 0.0; flux = True
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--roll-hz="): roll_hz = float(a[10:])
        elif a == "--noflux": flux = False
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    t0 = time.time()
    if part == "cases":
        part_cases(roll_hz, flux)
    else:
        import bearing_round9_parts as p9
        getattr(p9, "part_" + part)()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
