#!/usr/bin/env python3
"""
The deployment runs of round 8 (2026-11-15, Vol_III/Claude_Code_instructions_round8_book_cable_2026-11-15.md): round 7's set repeated on the
book's cable.  The cable keeps its hardware allowance of 10,686 kg/m (Paul's ruling); the bearing of round 7 uses 10,057 kg/m of it.  Everything
else is as in round 7: the plates as ruled for the cable's radiation, 4.0 m of cold width, ten fast families connected with the 0.1 mm
regulation's heat, the extra cold load now from the accepted fast track's own solve (2,258 A).
Book: Orbital Ring Engineering, Volume III: Chapters 3, 9, 11 and 14.
Run: python tools/round8_design_point.py --list
     python tools/round8_design_point.py --part=slip --run [--jobs=5]     the sweep: 120 to 280 turns at 12 MW, 150 to 250 at 8 MW
     python tools/round8_design_point.py --part=ch9 --run                 the Chapter 9 set (cases the sweep already has are copied)
     python tools/round8_design_point.py --part=ch14                      the Chapter 14 series of the design run
     python tools/round8_design_point.py --part=prograde                  the prograde run for Figure 3.2
The runs go to runs/round8_design/.  Environment: ORBITAL_RING_OUT (default runs).  Everything runs under SLIP_RULE optimal.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import json
import os
import subprocess
import sys

import ring_altitude as ring
import round7_design_point as r7

HW_BOOK = 10686.0                                                    # kg/m, the cable's hardware allowance (Table 7.1)
BASE = os.path.join(ring.out_dir(), "round8_design")
HERE = r7.HERE


def design_point():
    """Round 7's design point with the cable's hardware held at the book's allowance."""
    dp = r7._design_point_r7()
    dp["hw_bearing"] = dp["hw"]                                       # what the six plates and the LIM plates weigh
    dp["hw"] = HW_BOOK                                                # what the cable is sized for
    dp["margin"] = HW_BOOK - dp["hw_bearing"]
    return dp


def overrides(dp, N=165, mw=12):
    """lim_config's overrides.  The plates' dimensions are left at lim_config's defaults (four plates of 1.4 m x 200 mm), because the simulation takes the
    cable's hardware from them: 4 x 7,870 x 1.4 x 0.2 + 1,872 = 10,686 kg/m, the book's allowance.  They enter the run through that mass only."""
    return {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": N, "LEV_COIL_WIDTH": dp["cold_width"], "Q_COLD_EXTRA_PER_SITE": dp["extra"]}


def run_case(part, rule, name, spec, full=True):
    """round7_design_point.run_case with the book's hardware (see overrides)."""
    import shutil, time
    d = r7.d; DP = r7.DP
    lims, N, mw = spec
    out = os.path.join(d.case_dir(part, rule), name)
    os.makedirs(out, exist_ok=True)
    if part == "ch9":
        src = os.path.join(d.case_dir("slip", rule), f"N{N}_P{mw}MW")
        if lims == 2 and os.path.exists(os.path.join(src, "rising_radius_summary.json")):
            for f in d.FILES:
                shutil.copy2(os.path.join(src, f), os.path.join(out, f))
            return name, 0, 0.0
    over = overrides(DP, N, mw); hw = DP["hw"]
    if lims > 2:
        over["N_PLATES_PER_SIDE"] = lims // 2
        hw += r7.M_LIM_PLATES * (lims // 2 - 1)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_OUT=out, ORBITAL_RING_SLIP_RULE=rule, ORBITAL_RING_HW_KGM=f"{hw:.1f}")
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; runpy.run_path('lim_simulation.py', run_name='__main__')")
    t0 = time.time()
    with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
        r = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    return name, r.returncode, time.time() - t0


def main():
    part = next((a[7:] for a in sys.argv if a.startswith("--part=")), "slip")
    r7.BASE = BASE
    r7._design_point_r7 = r7.design_point
    r7.design_point = design_point
    r7.run_case = run_case
    if part in ("slip", "ch9") or "--list" in sys.argv:
        r7.main()
        return
    dp = json.load(open(os.path.join(BASE, "design_point.json"), encoding="utf-8"))
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_SLIP_RULE="optimal", ORBITAL_RING_HW_KGM=f"{dp['hw']:.1f}")
    if part == "ch14":
        out = os.path.join(BASE, "ch14"); os.makedirs(out, exist_ok=True)
        if os.path.exists(os.path.join(out, "started.txt")):
            print("ch14: already started or done"); return
        open(os.path.join(out, "started.txt"), "w").write("started")
        env["ORBITAL_RING_OUT"] = out
        args = [sys.executable, "tools/ch14_thermal_series.py", os.path.join(out, "ch14_thermal_series.csv"), "--rising"]
        for k, v in overrides(dp).items():
            args += ["--set", f"{k}={v!r}"]
        with open(os.path.join(out, "console.txt"), "w", encoding="utf-8", errors="replace") as lf:
            rc = subprocess.run(args, cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
        print("ch14 rc", rc)
    elif part == "prograde":
        out = os.path.join(BASE, "prograde_N165_P12MW"); os.makedirs(out, exist_ok=True)
        if os.path.exists(os.path.join(out, "started.txt")):
            print("prograde: already started or done"); return
        open(os.path.join(out, "started.txt"), "w").write("started")
        env["ORBITAL_RING_OUT"] = out; env["ORBITAL_RING_CABLE"] = "prograde"
        code = (f"import lim_config as c; c.configure(**{overrides(dp)!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; runpy.run_path('lim_simulation.py', run_name='__main__')")
        with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
            rc = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
        print("prograde rc", rc)


if __name__ == "__main__":
    main()
