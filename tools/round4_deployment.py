#!/usr/bin/env python3
"""
Round 4 (2026-11-13), S7: the design deployment run (150 turns, 12 MW per site, two LIMs, full steps, the rising quadratic) for each plate
candidate of Vol_III/Claude_Code_instructions_bearing_round4_2026-11-13.md, with the candidate's hardware mass in ring_altitude
(ORBITAL_RING_HW_KGM) and in lim_config (N_LEV_PLATES, D_LEV, W_LEV) and the cold width the cable radiates onto (LEV_COIL_WIDTH).
Book: Orbital Ring Engineering, Volume I, Chapters 11 and 14.

Run: python tools/round4_deployment.py --run=book_3.6,b,c,d [--jobs=4]   the runs, in parallel, into runs/round4_hw/<tag>/lim_design (about half an hour each)
     python tools/round4_deployment.py --copy                             copy the run of record (book_4.0) and Stage B's six-plate run (a) in, with their meta
     python tools/round4_deployment.py --list                             the candidates and their hardware
Environment: ORBITAL_RING_OUT (output folder, default runs).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json
import os
import shutil
import subprocess
import sys
import time

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
RHO_PLATE = 7870.0
M_LIM_PLATES = 1872.0
OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
BASE = os.path.join(OUT, "round4_hw")
T_D = (10686.0 - M_LIM_PLATES) / (6 * RHO_PLATE * 1.2)
# tag: (support plates (W, T, count), fast plates (W, T, count), cold width m, note)
CANDIDATES = {
    "book_4.0": ((1.4, 0.20, 2), (1.4, 0.20, 2), 4.0, "the run of record: four plates, the book's 1 m per track count"),
    "book_3.6": ((1.4, 0.20, 2), (1.4, 0.20, 2), 3.6, "the book's four plates with the loops' count of the cold width"),
    "a": ((1.2, 0.15, 2), (1.2, 0.15, 4), 5.2, "six plates of 1.2 m x 150 mm (Stage B's run)"),
    "b": ((1.2, 0.10, 2), (1.2, 0.20, 4), 5.2, "fast plates 1.2 m x 200 mm, support plates 1.2 m x 100 mm"),
    "c": ((1.2, 0.15, 2), (1.2, 0.20, 4), 5.2, "fast plates 1.2 m x 200 mm, support plates 1.2 m x 150 mm"),
    "d": ((1.2, T_D, 2), (1.2, T_D, 4), 5.2, "six equal plates at the book's 10,686 kg/m of hardware"),
}


def hw_of(tag):
    (ws, ts, ns), (wl, tl, nl), cold, note = CANDIDATES[tag]
    m_s = RHO_PLATE * ws * ts; m_l = RHO_PLATE * wl * tl
    hw = M_LIM_PLATES + ns * m_s + nl * m_l
    n = ns + nl; w_mean = (ns * ws + nl * wl) / n; t_eff = (ns * m_s + nl * m_l) / (n * RHO_PLATE * w_mean)
    return dict(tag=tag, hw=hw, n_plates=n, w_mean=w_mean, t_eff=t_eff, cold_width=cold, note=note, m_sup=m_s, m_lift=m_l)


def launch(tag):
    m = hw_of(tag)
    out = os.path.join(BASE, tag); os.makedirs(os.path.join(out, "lim_design"), exist_ok=True)
    with open(os.path.join(out, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(m, f, indent=1)
    over = {"MAX_SITE_POWER": 12e6, "N_TURNS": 150, "N_LEV_PLATES": m["n_plates"], "D_LEV": m["t_eff"], "W_LEV": m["w_mean"], "LEV_COIL_WIDTH": m["cold_width"]}
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; "
            f"runpy.run_path('lim_simulation.py', run_name='__main__')")
    env = dict(os.environ, ORBITAL_RING_HW_KGM=f"{m['hw']:.1f}", ORBITAL_RING_OUT=os.path.join(out, "lim_design"), MPLBACKEND="Agg", PYTHONUTF8="1")
    lf = open(os.path.join(out, "lim_design", "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace")
    p = subprocess.Popen([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    print(f"  launched {tag}: hardware {m['hw']:,.1f} kg/m, {m['n_plates']} plates (D_LEV {m['t_eff']*1e3:.1f} mm, W_LEV {m['w_mean']:.2f} m), cold width {m['cold_width']:.1f} m -> {out}", flush=True)
    return p, lf


def copy_in():
    for tag, src in (("book_4.0", os.path.join(OUT, "rising_radius_design", "L2_N150_P12MW")), ("a", os.path.join(OUT, "stageB_hw", "lim_design"))):
        s = os.path.join(src, "rising_radius_lim_simulation.txt")
        if not os.path.exists(s):
            print(f"  {tag}: {s} missing"); continue
        out = os.path.join(BASE, tag, "lim_design"); os.makedirs(out, exist_ok=True)
        shutil.copy(s, os.path.join(out, "rising_radius_lim_simulation.txt"))
        with open(os.path.join(BASE, tag, "meta.json"), "w", encoding="utf-8") as f:
            json.dump(dict(hw_of(tag), copied_from=s), f, indent=1)
        print(f"  {tag}: copied from {s}")


def main():
    jobs = 4; run = []; what = []
    for a in sys.argv[1:]:
        if a.startswith("--run="): run = a[6:].split(",")
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a in ("--help", "-h"): print(__doc__); return
        else: what.append(a)
    if "--list" in what or not (run or what):
        for tag in CANDIDATES:
            m = hw_of(tag)
            print(f"  {tag:>9}: hardware {m['hw']:,.1f} kg/m (support {m['m_sup']:,.0f} x {CANDIDATES[tag][0][2]}, fast {m['m_lift']:,.0f} x {CANDIDATES[tag][1][2]}, LIM {M_LIM_PLATES:,.0f}), cold width {m['cold_width']:.1f} m: {m['note']}")
    if "--copy" in what:
        copy_in()
    if run:
        t0 = time.time(); procs = []
        for tag in run:
            while len([p for p, lf in procs if p.poll() is None]) >= jobs:
                time.sleep(10)
            procs.append(launch(tag))
        for p, lf in procs:
            p.wait(); lf.close()
        print(f"  all done in {(time.time()-t0)/60:.0f} min")


if __name__ == "__main__":
    main()
