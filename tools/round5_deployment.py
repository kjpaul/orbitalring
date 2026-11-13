#!/usr/bin/env python3
"""
The design deployment run with round 5's optimum plates and the fast coils' idle lead heat (2026-11-13, Section 6 of the round 5 brief).
Book: Orbital Ring Engineering, Volume III: Chapter 14 (the cold loop), Chapter 15 (the deployment).
Run: python tools/round5_deployment.py --list                         the runs and their hardware, cold width and extra cold load
     python tools/round5_deployment.py --run=opt_0,opt_10,opt_14 [--jobs=3]
       opt_0   the optimum plates, no lead heat (the plates alone against round 4's runs)
       opt_10  the optimum plates and the idle heat of the fast coils' leads at 10 series turns
       opt_14  the same at 14 turns
     each run: lim_simulation --rising-radius with 150 turns, 12 MW per site, two LIMs, full steps, into runs/round5_hw/<tag>/lim_design (about 35 min)
Environment: ORBITAL_RING_OUT (default runs).  Reads runs/bearing_round5_final.json (the plates) and runs/bearing_round5_coldload.json (the leads).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json
import os
import subprocess
import sys
import time

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
BASE = os.path.join(OUT, "round5_hw")
M_LIM_PLATES = 1872.0
RHO_PLATE = 7870.0


def inputs():
    fin = json.load(open(os.path.join(OUT, "bearing_round5_final.json"), encoding="utf-8"))
    rs = fin["results"]["support:opt135"]; rl = fin["results"]["fast:opt135"]
    cold = json.load(open(os.path.join(OUT, "bearing_round5_coldload.json"), encoding="utf-8"))
    return rs, rl, cold


def runs():
    rs, rl, cold = inputs()
    m_s, m_l = rs["plate_m"], rl["plate_m"]
    hw = M_LIM_PLATES + 2 * m_s + 4 * m_l
    w_mean = (2 * rs["Wp"] + 4 * rl["Wp"]) / 6; t_eff = (2 * m_s + 4 * m_l) / (6 * RHO_PLATE * w_mean)
    cold_width = 2 * (2 * rs["Wl"]) + 4 * rl["Wl"]                        # the loops' count (round 4, S7): two loops per support track, one per fast track
    out = {}
    for tag, extra in (("opt_0", 0.0), ("opt_10", cold["leads"]["10"]["idle_site"]), ("opt_14", cold["leads"]["14"]["idle_site"])):
        out[tag] = dict(hw=hw, m_sup=m_s, m_lift=m_l, w_mean=w_mean, t_eff=t_eff, cold_width=cold_width, extra=extra,
                        sup=(rs["Wp"], rs["T"], rs["Wl"]), fast=(rl["Wp"], rl["T"], rl["Wl"]))
    return out


def launch(tag, cfg):
    out = os.path.join(BASE, tag); os.makedirs(os.path.join(out, "lim_design"), exist_ok=True)
    with open(os.path.join(out, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=1)
    over = {"MAX_SITE_POWER": 12e6, "N_TURNS": 150, "N_LEV_PLATES": 6, "D_LEV": cfg["t_eff"], "W_LEV": cfg["w_mean"],
            "LEV_COIL_WIDTH": cfg["cold_width"], "Q_COLD_EXTRA_PER_SITE": cfg["extra"]}
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; "
            f"runpy.run_path('lim_simulation.py', run_name='__main__')")
    env = dict(os.environ, ORBITAL_RING_HW_KGM=f"{cfg['hw']:.1f}", ORBITAL_RING_OUT=os.path.join(out, "lim_design"), MPLBACKEND="Agg", PYTHONUTF8="1")
    lf = open(os.path.join(out, "lim_design", "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace")
    p = subprocess.Popen([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    print(f"  launched {tag}: hardware {cfg['hw']:,.1f} kg/m, cold width {cfg['cold_width']:.2f} m, extra cold load {cfg['extra']:,.0f} W per site -> {out}", flush=True)
    return p


def main():
    what = None; jobs = 3
    for a in sys.argv[1:]:
        if a.startswith("--run="): what = a[6:].split(",")
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a in ("-h", "--help"): print(__doc__); return
    R = runs()
    if what is None:
        for tag, c in R.items():
            print(f"  {tag}: hardware {c['hw']:,.1f} kg/m (support {c['sup'][0]:.2f} x {c['sup'][1]*1e3:.0f} mm = {c['m_sup']:,.0f}, fast {c['fast'][0]:.2f} x {c['fast'][1]*1e3:.0f} mm = {c['m_lift']:,.0f}), "
                  f"D_LEV {c['t_eff']*1e3:.1f} mm, W_LEV {c['w_mean']:.2f} m, cold width {c['cold_width']:.2f} m, extra {c['extra']:,.0f} W per site")
        return
    procs = []; t0 = time.time()
    for tag in what:
        procs.append((tag, launch(tag, R[tag])))
        while sum(1 for _, p in procs if p.poll() is None) >= jobs:
            time.sleep(10)
    for tag, p in procs:
        p.wait(); print(f"  {tag} done (rc {p.returncode}) at {(time.time()-t0)/60:.0f} min", flush=True)


if __name__ == "__main__":
    main()
