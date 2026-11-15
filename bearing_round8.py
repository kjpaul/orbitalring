#!/usr/bin/env python3
"""
The bearing, round 8 (2026-11-15): the rows the Chapter 7 rewrite needs (Section 3b of
Vol_III/Claude_Code_instructions_round8_book_cable_2026-11-15.md).
Book: Orbital Ring Engineering, Volume III: Chapter 7 (Tables 7.3, 7.4 and 7.6, Sections 7.5 to 7.9), Chapter 8, Chapter 14, Chapter 16.
Run (from Python/orbitalring):
  python bearing_round8.py --part=solve-sup --set=both|high|one   the split support track at 2.5 mm: the pull curve from 0.6 to 1.5 times the cap's
                                                                 ampere-turns with the worst tape at each level, and one loop alone at 1.0, 1.25 and 1.5
  python bearing_round8.py --part=rows                            the rows of Section 3b from the solves (bearing_round8_parts.py)
Writes runs/round8_logs/.  Environment: ORBITAL_RING_OUT (default runs).  US spelling in the printed text.
"""
import json
import math
import os
import sys
import time

import numpy as np

import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round6b as r6
import bearing_round7 as r7

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
LOGS = os.path.join(OUT, "round8_logs")
SUP = r7.SUP; GC = 0.05


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else (o.tolist() if isinstance(o, np.ndarray) else float(o)))


def load(name):
    p = os.path.join(LOGS, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


SETS = {"both": [("both", 0.6), ("both", 0.8), ("both", 1.0), ("lin", 0.0)],
        "high": [("both", 1.15), ("both", 1.3), ("both", 1.5)],
        "one": [("one", 1.0), ("one", 1.25), ("one", 1.5)],
        "norm": [("both", 0.358656), ("one", 0.586531)]}                 # the normal state (156.7 kA) and one loop alone carrying 47.5 kN/m (256.3 kA), for the plate's field


def job_sup(which, h=0.0025):
    """The split support track (two 0.4 m loops, 50 mm between the center bundles, the 1.2 m x 100 mm plate) at multiples of the cap's ampere-turns:
    the pull, the plate's saturation, the loops' flux and the worst 12 mm tape with its 70 K and 77.5 K limits; 'one' is loop 0 alone."""
    t0 = time.time()
    half = r7.load("round7_half_2p5.json"); NI_cap = half["cap"]["NI"]
    theva, theva77 = be.load_theva()
    geom = dict(bo.BEST_SUPPORT); geom["W"] = SUP["Wl"]; geom["gc"] = GC
    base = r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"], 0.005) or r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"])
    rows = []
    with r4.plate_patch(SUP["Wp"], SUP["T"]):
        tr = be.Track(geom, h, w_tape=0.012, I_ref=bo.I_REF_12, NI_ref=base["NI_req"])
        for kind, k in SETS[which]:
            if kind == "lin":
                g, s = tr.solve(tr.NI_design(1e3), linear=True)
                rows.append(dict(kind="lin", NI=1e3, F=s["F"], phi=s["phi"]))
                print(f"    linear: F / NI^2 = {s['F']/1e6:.6e} N/m per A^2 ({(time.time()-t0)/60:.0f} min)", flush=True)
                continue
            NI = k * NI_cap
            legs = tr.NI_design(NI) if kind == "both" else tr.NI_loop(0, NI)
            g, s = tr.solve(legs)
            wt, lim = be.worst_tape(g, tr, legs, 0.012, theva, theva77)
            ic77 = min((d["ic77_" + wt[3]] for d in lim if abs(d["x"] - wt[4]) < 1e-9 and abs(d["ic_" + wt[3]] - wt[0]) < 1e-9), default=None)
            rows.append(dict(kind=kind, k=k, NI=NI, F=s["F"], sat=s["sat"], B_plate=s["B_plate"], phi=s["phi"], W=s["W"], ic=wt[0], B=wt[1], th=wt[2], face=wt[3], x=wt[4], ic77=ic77))
            print(f"    {kind} {k}: {NI/1e3:.1f} kA, {s['F']/1e3:.1f} kN/m, sat {s['sat']*100:.1f} %, worst tape {wt[0]:.1f} A at {wt[1]:.2f} T, {wt[2]:.0f} deg (77.5 K: {ic77}) ({(time.time()-t0)/60:.0f} min)", flush=True)
            save(f"round8_sup_{which}.json", dict(h=h, gc=GC, NI_cap=NI_cap, rows=rows))
    save(f"round8_sup_{which}.json", dict(h=h, gc=GC, NI_cap=NI_cap, rows=rows, seconds=time.time() - t0))
    return rows


def job_face(h=0.005):
    """The support plate's field map in normal operation (both loops at the solved 156.7 kA per leg) and with one loop alone at that current, at 5 mm
    cells: the field in the cells at the gap's face across the plate's width, for the coil ends' heat."""
    t0 = time.time()
    half = r7.load("round7_half_2p5.json"); NI_n = half["normal"]["NI"]
    geom = dict(bo.BEST_SUPPORT); geom["W"] = SUP["Wl"]; geom["gc"] = GC
    base = r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"], 0.005) or r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"])
    out = dict(h=h, NI=NI_n)
    with r4.plate_patch(SUP["Wp"], SUP["T"]):
        tr = be.Track(geom, h, w_tape=0.012, I_ref=bo.I_REF_12, NI_ref=base["NI_req"])
        for tag, legs in (("both", tr.NI_design(NI_n)), ("one", tr.NI_loop(0, NI_n))):
            g, s = tr.solve(legs)
            m = g.box(-SUP["Wp"] / 2, SUP["Wp"] / 2, bd.GAP0, bd.GAP0 + SUP["T"])
            ii, jj = np.nonzero(m)
            B = np.hypot(g.Bx[m], g.By[m]); j0 = jj.min()
            f = jj == j0
            out[tag] = dict(F=s["F"], x=(g.x0 + ii[f] * g.h).tolist(), B_face=B[f].tolist(), B_max=float(B.max()), B_mean=float(B.mean()),
                            B_colmean=[float(B[ii == i].mean()) for i in sorted(set(ii.tolist()))])
            print(f"    {tag}: {s['F']/1e3:.1f} kN/m; plate field: peak {B.max():.2f} T, mean {B.mean():.2f} T; at the face: peak {B[f].max():.2f} T, rms {math.sqrt((B[f]**2).mean()):.2f} T ({(time.time()-t0)/60:.0f} min)", flush=True)
    save("round8_face.json", out)
    return out


def main():
    part = "rows"; which = "both"
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--set="): which = a[6:]
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    t0 = time.time()
    if part == "solve-sup":
        job_sup(which)
    elif part == "solve-face":
        job_face()
    else:
        import bearing_round8_parts as p8
        getattr(p8, "part_" + part)()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
