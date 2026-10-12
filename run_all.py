#!/usr/bin/env python3
"""
Run every simulation behind The Orbital Ring, Volume IV, at the design point of ring_altitude.py, and keep one
log per run.  This is the same list as run_700.sh and run_700_lim.sh, for any platform (no bash needed).

Usage:
    python run_all.py                 run everything (about 3 hours with the solver caches in data/)
    python run_all.py --list          print the runs and stop
    python run_all.py --only=a,b      run only the named logs (names as in --list, without .txt)
    python run_all.py --skip=a,b      run everything except the named logs
    python run_all.py --no-lim        skip the two LIM deployment runs (about 40 minutes)
    python run_all.py --lim-only      run only the two LIM deployment runs
    python run_all.py --dry-run       print the commands and stop

Logs go to <ORBITAL_RING_OUT>/<name>.txt (default runs/), figures and CSV files into subfolders of the same
place.  A timing table is written to <ORBITAL_RING_OUT>/_timing.txt and the file _done marks the end.

Environment: ORBITAL_RING_OUT (output folder, default runs); ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE,
ORBITAL_RING_HW_KGM, ORBITAL_RING_ARRAY_OFFSET_KM pass through to every script (see ring_altitude.py).
The runs that set their own environment (the prograde, 250 km and 800 km reference tables, the ring-mounted
power cases) do so on top of what is passed in.

Reference: "The Orbital Ring, Volume IV: Space Mass Transit Complete Design Review" by Paul G. de Jong
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# (log name, arguments after "python", environment on top of the caller's, stderr merged into the log)
RUNS = [
    ("ring_altitude_700", ["ring_altitude.py"], {}, False),
    ("ring_altitude_700_prograde", ["ring_altitude.py"], {"ORBITAL_RING_CABLE": "prograde"}, False),
    ("ring_altitude_250", ["ring_altitude.py"], {"ORBITAL_RING_ALT_KM": "250"}, False),
    ("ring_altitude_800", ["ring_altitude.py"], {"ORBITAL_RING_ALT_KM": "800"}, False),
]
for i in (0, 5, 10, 20, 30):
    RUNS.append((f"anchor_700_i{i}", ["anchor_line_analysis.py", f"--inclination={i}"], {}, True))
for i in (5, 10, 20, 30):
    RUNS.append((f"j2_retrograde_i{i}", ["j2_simulation.py", f"--inclination={i}"], {}, True))
    RUNS.append((f"j2_prograde_i{i}", ["j2_simulation.py", f"--inclination={i}", "--cable=prograde"], {}, True))
RUNS += [
    ("power_deployment", ["power_simulation.py"], {}, True),
    ("power_ops", ["power_simulation.py", "--demand=ops"], {}, True),
    ("power_deployment_ring_mounted_8MW", ["power_simulation.py"], {"ORBITAL_RING_ARRAY_OFFSET_KM": "0"}, True),
    ("power_deployment_ring_mounted_16MW", ["power_simulation.py", "--power=16"], {"ORBITAL_RING_ARRAY_OFFSET_KM": "0"}, True),
    ("levitation", ["levitation_simulation.py", "all"], {}, True),
    ("tether_stress", ["tether_stress.py"], {}, True),
    ("rendezvous", ["rendezvous_analysis.py"], {}, True),
    ("freefall", ["freefall_trajectory.py"], {}, True),
    ("atmospheric_loading", ["atmospheric_loading.py"], {}, True),
    ("climber_descent", ["climber_descent.py"], {}, True),
    ("lsm_chapter_metrics", ["lsm_chapter_metrics.py"], {}, True),
    ("lsm_n_optimization", ["lsm_n_optimization_demo.py"], {}, True),
    ("lsm_simulation", ["lsm_simulation.py"], {}, True),
    ("md_simulation", ["md_simulation.py", "--no-graphs"], {}, True),
    ("md_trade_study", ["md_trade_study.py"], {}, True),
    ("debris_simulation", ["debris_simulation.py"], {}, True),
    ("debris_orbits", ["debris_orbits.py"], {}, True),
    ("ch4_numbers", ["ch4_numbers.py"], {}, True),
    ("altitude_trade", ["altitude_trade.py"], {}, True),
    ("sled_passage", ["sled_passage_simulation.py"], {}, True),
    ("sled_passage_limits", ["sled_passage_limits.py"], {}, True),
    ("sled_passage_pulsed", ["sled_passage_pulsed.py"], {}, True),
    ("sled_passage_spread", ["sled_passage_spread.py"], {}, True),
    ("compensation_fields", ["compensation_fields.py"], {}, True),
    ("compensation_tracks", ["compensation_tracks.py"], {}, True),
    ("bearing_fem_check", ["bearing_fem_check.py"], {}, True),
    ("bearing_design", ["bearing_design.py"], {}, True),
    ("bearing_mass", ["bearing_mass.py"], {}, True),
    ("final_bearing", ["final_bearing.py"], {}, True),
]
# The LIM deployment at 8 MW per site, and at 16 MW with the same command run_700_lim.sh uses.
LIM_16MW = ("import lim_config; lim_config.MAX_SITE_POWER=16e6; import runpy,sys; sys.argv=['lim_simulation.py']; "
            "runpy.run_path('lim_simulation.py', run_name='__main__')")
LIM_RUNS = [
    ("lim_simulation", ["lim_simulation.py"], {}, True),
    ("lim_simulation_16MW", ["-c", LIM_16MW], {}, True),
]


def main():
    only = skip = None
    no_lim = lim_only = dry = False
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg == "--list":
            for name, args, env, _ in RUNS + LIM_RUNS:
                print(f"{name:36} {' '.join(f'{k}={v}' for k, v in env.items())} python {' '.join(args)}"[:160])
            return
        elif arg.startswith("--only="):
            only = set(arg[7:].split(","))
        elif arg.startswith("--skip="):
            skip = set(arg[7:].split(","))
        elif arg == "--no-lim":
            no_lim = True
        elif arg == "--lim-only":
            lim_only = True
        elif arg == "--dry-run":
            dry = True
        else:
            print(f"Unknown option {arg}"); return
    runs = ([] if lim_only else RUNS) + ([] if no_lim else LIM_RUNS)
    if only is not None:
        runs = [r for r in runs if r[0] in only]
    if skip is not None:
        runs = [r for r in runs if r[0] not in skip]

    out = os.environ.get("ORBITAL_RING_OUT", "runs")
    os.makedirs(out, exist_ok=True)
    base_env = dict(os.environ, MPLBACKEND="Agg", ORBITAL_RING_OUT=out)
    timing = os.path.join(out, "_timing.txt")
    if not dry:
        with open(timing, "w", encoding="utf-8") as tf:
            tf.write("script seconds\n")
    t_all = time.time()
    failed = []
    for name, args, env, merge in runs:
        log = os.path.join(out, name + ".txt")
        cmd = [sys.executable] + args
        print(f"{name:36} -> {log}", flush=True)
        if dry:
            continue
        t0 = time.time()
        with open(log, "w", encoding="utf-8", errors="replace") as lf:
            r = subprocess.run(cmd, cwd=HERE, env=dict(base_env, **env), stdout=lf,
                               stderr=subprocess.STDOUT if merge else None)
        dt = time.time() - t0
        with open(timing, "a", encoding="utf-8") as tf:
            tf.write(f"{name}.txt {dt:.0f}\n")
        if r.returncode != 0:
            failed.append(name)
            print(f"    exit code {r.returncode}", flush=True)
    if dry:
        return
    with open(timing, "a", encoding="utf-8") as tf:
        tf.write("ALLDONE\n")
    with open(os.path.join(out, "_done"), "w", encoding="utf-8") as f:
        f.write("ALLDONE\n")
    print(f"{len(runs)} runs in {(time.time() - t_all) / 60:.0f} min; failed: {', '.join(failed) if failed else 'none'}")


if __name__ == "__main__":
    main()
