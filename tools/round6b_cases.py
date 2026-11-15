#!/usr/bin/env python3
"""Round 6b, Section 1.1: round 4 cases rerun with the coupling of the solve legs (run from Python/orbitalring after bearing_round6b.py --part=legs)."""
import json, os, sys
sys.path.insert(0, os.getcwd())
import bearing_round6a as a
rows = json.load(open("runs/round6b_logs/round6b_legs.json"))["rows"]
sr = [dict(d=r["d"], dF=r["dF"]) for r in rows if r["d"] in (0.2, 0.6)]
a.LOGS = "runs/round6b_logs"
import bearing_layout as bl
_rc = bl.run_cases
def rc(tau, prot, A2, theva, sup, lift, tabs, fname, **kw):
    return _rc(tau, prot, A2, theva, sup, lift, tabs, fname.replace("round6a", "round6b"), **kw)
bl.run_cases = rc
a.part_cases((0.2, 0.6), spacing_rows=sr)
