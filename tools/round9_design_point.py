#!/usr/bin/env python3
"""
The deployment runs of round 9 (2026-11-16, Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md, item 2): round 8's set repeated with the
bearing's cold load as the design makes it (bearing_cold_load.py: the fast leads' conduction with ten families, the regulation's heat with the
support following the weight continuously, the support's leads until the pumps can keep up, then the pumps), in place of round 8's flat 4,411 W.
Everything else is round 8's: the book's cable (hardware 10,686 kg/m), 165 turns, 12 MW, two LIMs, SLIP_RULE optimal, 4.0 m of cold width.
Book: Orbital Ring Engineering, Volume III: Chapters 3, 9, 11 and 14.
Run: python tools/round9_design_point.py --part=slip --run [--jobs=6]    the sweep: 120 to 280 turns at 12 MW, 150 to 250 at 8 MW
     python tools/round9_design_point.py --part=ch9 --run                the Chapter 9 set
     python tools/round9_design_point.py --part=ch14 [--pump=dynamo]     the Chapter 14 series of the design run (the dynamo as a second case)
     python tools/round9_design_point.py --part=prograde                 the prograde run for Figure 3.2
     add --table / --figure as in round 7.  The runs go to runs/round9_design/.
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
import round8_design_point as r8

BASE = os.path.join(ring.out_dir(), "round9_design")
HERE = r7.HERE


def table_path(pump="mallett"):
    return os.path.abspath(os.path.join(ring.out_dir(), "round9_logs", f"bearing_cold_load_{pump}.json")).replace("\\", "/")


def overrides(dp, N=165, mw=12, pump="mallett"):
    return {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": N, "LEV_COIL_WIDTH": dp["cold_width"], "Q_COLD_EXTRA_PER_SITE": 0.0, "BEARING_COLD_LOAD": table_path(pump)}


def one(out, dp, pump="mallett", prograde=False, ch14=False):
    os.makedirs(out, exist_ok=True)
    if os.path.exists(os.path.join(out, "started.txt")):
        print(out, ": already started or done"); return
    open(os.path.join(out, "started.txt"), "w").write("started")
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_SLIP_RULE="optimal", ORBITAL_RING_HW_KGM=f"{dp['hw']:.1f}", ORBITAL_RING_OUT=out)
    if prograde:
        env["ORBITAL_RING_CABLE"] = "prograde"
    ov = overrides(dp, pump=pump)
    if ch14:
        code = (f"import sys; sys.path.insert(0, 'tools'); import lim_config as c; c.configure(**{ov!r}); import runpy; "
                f"sys.argv=['tools/ch14_thermal_series.py', {os.path.join(out, 'ch14_thermal_series.csv')!r}, '--rising']; runpy.run_path('tools/ch14_thermal_series.py', run_name='__main__')")
        log = "console.txt"
    else:
        code = (f"import lim_config as c; c.configure(**{ov!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; runpy.run_path('lim_simulation.py', run_name='__main__')")
        log = "rising_radius_lim_simulation.txt"
    with open(os.path.join(out, log), "w", encoding="utf-8", errors="replace") as lf:
        rc = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
    print(out, "rc", rc)


def main():
    part = next((a[7:] for a in sys.argv if a.startswith("--part=")), "slip")
    pump = next((a[7:] for a in sys.argv if a.startswith("--pump=")), "mallett")
    r8.BASE = BASE; r8.overrides = overrides
    r7.BASE = BASE
    r7._design_point_r7 = r7.design_point
    r7.design_point = r8.design_point
    r7.run_case = r8.run_case
    if part in ("slip", "ch9") or "--list" in sys.argv:
        r7.main()
        return
    dp = json.load(open(os.path.join(BASE, "design_point.json"), encoding="utf-8"))
    if part == "ch14":
        one(os.path.join(BASE, "ch14" if pump == "mallett" else f"ch14_{pump}"), dp, pump=pump, ch14=True)
    elif part == "prograde":
        one(os.path.join(BASE, "prograde_N165_P12MW"), dp, prograde=True)


if __name__ == "__main__":
    main()
