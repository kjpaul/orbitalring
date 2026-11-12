#!/usr/bin/env python3
"""
The design point with the six-plate bearing of Stage B (2026-11-12): the cable hardware mass from the chosen plates, then the reruns the
numbers sheet needs, old (the present design point, 10,686 kg/m of hardware, runs/rising_quadratic/all) beside new.
Book: Orbital Ring Engineering, Volume I: the re-keying of the bearing's plates (Chapters 3, 7, 9, 14 and the numbers sheet).

Run: python tools/stageB_design_point.py --hw                      print the hardware mass of the chosen plates and stop
     python tools/stageB_design_point.py --ring                    ring_altitude's report under the new hardware (seconds)
     python tools/stageB_design_point.py --all-runs [--jobs=N]      run_all.py --no-lim under the new hardware into runs/stageB_hw (hours)
     python tools/stageB_design_point.py --lim                     the design deployment run (150 turns, 12 MW, two LIMs, full steps) with the
                                                                   six plates in lim_config, into runs/stageB_hw/lim_design (about half an hour)
     python tools/stageB_design_point.py --sheet                   the numbers sheet old beside new, runs/stageB_hw/numbers_stageB.md
     python tools/stageB_design_point.py --arith                   tools/design_point_numbers.py old and new, side by side
Environment: ORBITAL_RING_OUT (output folder, default runs).  The plates: --sup=W,T --lift=W,T in metres (defaults: the Stage B choice below).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import os
import subprocess
import sys
import time

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
RHO_PLATE = 7870.0
M_LIM_PLATES = 1872.0                       # the two LIM reaction plates (gamma-TiAl, 0.20 x 1.2 m)
# the Stage B choice (S6): support plates and the four pulsed plates (lift and load), width and thickness in metres
SUP = (1.2, 0.15)
LIFT = (1.2, 0.15)
OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
OLD = os.path.join(OUT, "rising_quadratic", "all")
NEW = os.path.join(OUT, "stageB_hw")


def hw(sup=SUP, lift=LIFT):
    m_sup = RHO_PLATE * sup[0] * sup[1]; m_lift = RHO_PLATE * lift[0] * lift[1]
    return M_LIM_PLATES + 2 * m_sup + 4 * m_lift, m_sup, m_lift


def main():
    sup, lift = SUP, LIFT; jobs = 4; what = []
    for a in sys.argv[1:]:
        if a.startswith("--sup="): sup = tuple(float(v) for v in a[6:].split(","))
        elif a.startswith("--lift="): lift = tuple(float(v) for v in a[7:].split(","))
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a in ("--help", "-h"): print(__doc__); return
        else: what.append(a)
    m_hw, m_sup, m_lift = hw(sup, lift)
    env = dict(os.environ, ORBITAL_RING_HW_KGM=f"{m_hw:.1f}", MPLBACKEND="Agg", PYTHONUTF8="1")
    print(f"Six plates: support {sup[0]:.2f} m x {sup[1]*1e3:.0f} mm = {m_sup:,.0f} kg/m (two), lift and load {lift[0]:.2f} m x {lift[1]*1e3:.0f} mm = {m_lift:,.0f} kg/m (four); LIM plates {M_LIM_PLATES:,.0f};")
    print(f"cable hardware {m_hw:,.1f} kg/m (was 10,686.4: four 1.4 m x 200 mm plates, 8,814, plus the LIM plates).")
    if "--hw" in what or not what:
        return
    if "--ring" in what:
        subprocess.run([sys.executable, "ring_altitude.py"], cwd=HERE, env=env)
    if "--all-runs" in what:
        os.makedirs(NEW, exist_ok=True)
        env2 = dict(env, ORBITAL_RING_OUT=NEW)
        t0 = time.time()
        with open(os.path.join(NEW, "_run_all_console.txt"), "w", encoding="utf-8", errors="replace") as lf:
            subprocess.run([sys.executable, "run_all.py", "--no-lim"], cwd=HERE, env=env2, stdout=lf, stderr=subprocess.STDOUT)
        print(f"run_all --no-lim done in {(time.time()-t0)/60:.0f} min -> {NEW}")
    if "--lim" in what:
        # lim_config's own plate model must carry the six plates too (its hardware and the levitation coils' cold area):
        # N_LEV_PLATES x DENSITY_IRON x D_LEV x W_LEV = the six plates' mass, with the lift plates' thickness and the mean width
        w_mean = (2 * sup[0] + 4 * lift[0]) / 6; t_eff = (2 * m_sup + 4 * m_lift) / (6 * RHO_PLATE * w_mean)
        over = {"MAX_SITE_POWER": 12e6, "N_TURNS": 150, "N_LEV_PLATES": 6, "D_LEV": t_eff, "W_LEV": w_mean,
                "LEV_COIL_WIDTH": 2 * 1.0 + 4 * 0.8}       # the cold area the cable radiates onto: two 1 m support tracks and four 0.8 m driven tracks
        out = os.path.join(NEW, "lim_design"); os.makedirs(out, exist_ok=True)
        code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; "
                f"runpy.run_path('lim_simulation.py', run_name='__main__')")
        env2 = dict(env, ORBITAL_RING_OUT=out)
        t0 = time.time()
        with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
            r = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env2, stdout=lf, stderr=subprocess.STDOUT)
        print(f"design run done (rc {r.returncode}) in {(time.time()-t0)/60:.0f} min -> {out}; lim_config overrides {over}")
    if "--sheet" in what:
        subprocess.run([sys.executable, "tools/number_sheet_hw10686.py", f"--old={OLD}", f"--new={NEW}", f"--out={os.path.join(NEW, 'numbers_stageB.md')}"], cwd=HERE, env=env)
    if "--arith" in what:
        print("---- the present hardware (10,686 kg/m):")
        subprocess.run([sys.executable, "tools/design_point_numbers.py"], cwd=HERE, env=dict(os.environ, PYTHONUTF8="1"))
        print(f"---- the six plates ({m_hw:,.1f} kg/m):")
        subprocess.run([sys.executable, "tools/design_point_numbers.py"], cwd=HERE, env=env)


if __name__ == "__main__":
    main()
