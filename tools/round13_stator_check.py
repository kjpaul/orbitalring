#!/usr/bin/env python3
"""
Round 13, item 1: checks that the ruled stator as the default of lsm_config.py moves nothing it should not move.
Book: Orbital Ring Engineering, Volume III: Chapters 7, 9, 10 (the launch profile the bearing scripts take).
Run: python tools/round13_stator_check.py     (about two minutes; writes runs/round13_logs/round13_stator_check.txt and .json)

1. The launch profile that the bearing scripts import (sled_passage_pulsed.lsm_profile, lsm_simulation in quick mode,
   and the full run's data_t, data_v, data_x, data_P_mech that sled_passage_spread reads): each run in a fresh Python
   under both stators, and the arrays compared.
2. The constants other scripts read from lsm_config (L_SLED, M_SLED_PER_M, M_SPACECRAFT, P_HVDC_MAX, M_TOTAL, V_LAUNCH,
   A_ACTIVE_TOTAL) under both stators.
3. The five numbers of the default run against the launch the book prints.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round13_logs")
CHILD = r'''
import contextlib, io, json, sys
import numpy as np
import lsm_config as cfg
import lsm_simulation as ls
out = {"stator": cfg.STATOR}
for quick in (True, False):
    with contextlib.redirect_stdout(io.StringIO()):
        ls.run_lsm_simulation(quick_mode=quick)
    key = "quick" if quick else "full"
    out[key] = dict(t=list(ls.data_t), v=list(ls.data_v), x=list(ls.data_x), P=list(ls.data_P_mech), F=list(ls.data_F))
out["const"] = dict(L_SLED=cfg.L_SLED, M_SLED_PER_M=cfg.M_SLED_PER_M, M_SPACECRAFT=cfg.M_SPACECRAFT, P_HVDC_MAX=cfg.P_HVDC_MAX, M_TOTAL=cfg.M_TOTAL,
                    V_LAUNCH=cfg.V_LAUNCH, A_ACTIVE_TOTAL=cfg.A_ACTIVE_TOTAL, M_SLED_HARDWARE=cfg.M_SLED_HARDWARE, L_RING=cfg.L_RING)
json.dump(out, open(sys.argv[1], "w"))
'''


def run(stator):
    path = os.path.join(LOGS, f"_profile_{stator}.json")
    env = dict(os.environ, ORBITAL_RING_LSM_STATOR=stator, MPLBACKEND="Agg")
    subprocess.run([sys.executable, "-c", CHILD, path], cwd=HERE, env=env, check=True)
    with open(path) as fh:
        d = json.load(fh)
    os.remove(path)
    return d


def main():
    import numpy as np
    os.makedirs(LOGS, exist_ok=True)
    a = run("ruled"); b = run("table96")
    lines = ["ROUND 13: THE RULED STATOR AS THE DEFAULT, WHAT MOVES IN THE LAUNCH PROFILE OTHER SCRIPTS TAKE",
             f"  stators compared: {a['stator']} (the default) and {b['stator']} (the 130 m stator, what the code ran before round 13)"]
    res = {}
    for key, lab in (("quick", "lsm_simulation quick mode, 1 s steps (sled_passage_pulsed.lsm_profile: t, v, P_mech)"), ("full", "lsm_simulation full run, 0.1 s steps (sled_passage_spread: data_x, data_v)")):
        lines.append(f"  {lab}")
        res[key] = {}
        for q, unit in (("t", "s"), ("v", "m/s"), ("x", "m"), ("P", "W"), ("F", "N")):
            xa = np.array(a[key][q]); xb = np.array(b[key][q])
            same_len = len(xa) == len(xb)
            dmax = float(np.max(np.abs(xa - xb))) if same_len else float("nan")
            rel = float(np.max(np.abs(xa - xb) / np.maximum(np.abs(xb), 1e-300))) if same_len else float("nan")
            res[key][q] = dict(n=len(xa), n_old=len(xb), max_abs=dmax, max_rel=rel, last=float(xa[-1]), last_old=float(xb[-1]))
            lines.append(f"    {q:>2}: {len(xa):,} steps ({len(xb):,} before); largest difference {dmax:.3e} {unit}, {rel:.2e} of the value; last value {xa[-1]:.6f} against {xb[-1]:.6f}")
    lines.append("  constants read by other scripts (ruled, table96):")
    res["const"] = {}
    for k in a["const"]:
        res["const"][k] = (a["const"][k], b["const"][k])
        lines.append(f"    {k:<16} {a['const'][k]!r:>24} {b['const'][k]!r:>24} {'same' if a['const'][k] == b['const'][k] else 'DIFFERENT'}")
    f = a["full"]; dt = f["t"][1] - f["t"][0]
    t_end = (f["t"][-1] + dt) / 60.0; x_end = (f["x"][-1] + f["v"][-1] * dt) / 1e3
    lines.append("  the default run against the book's launch:")
    lines.append(f"    peak thrust {max(f['F'])/1e6:.2f} MN (83.29); peak power {max(f['P'])/1e9:.2f} GW (416.47); time {t_end:.1f} min (189.7); distance {x_end:,.1f} km (222,015.5); laps {x_end*1e3/a['const']['L_RING']:.2f} (4.99)")
    res["five"] = dict(F_MN=max(f["F"]) / 1e6, P_GW=max(f["P"]) / 1e9, t_min=t_end, x_km=x_end, laps=x_end * 1e3 / a["const"]["L_RING"])
    with open(os.path.join(LOGS, "round13_stator_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(LOGS, "round13_stator_check.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
