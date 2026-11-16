#!/usr/bin/env python3
"""
The lighter casing on the present design (round 9, 2026-11-16, item 5 of Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md):
Section 16.6 and Table 16.2 run again with the bearing of rounds 7 and 8 and the 165-turn coil.
Book: Orbital Ring Engineering, Volume III: Chapter 16 (Section 16.6, Tables 16.1 and 16.2); Chapters 3, 7, 9, 11 and 14 for the rules.

The cases:  design   12,000 kg/m, the book's cable (hardware allowance 10,686 kg/m, the bearing uses 10,057)
            same     6,000 kg/m, the plates unchanged (hardware allowance 10,686 kg/m)
            kept     6,000 kg/m, the support plates sized for the lighter casing, the standing load of 198 kN/m per support track kept
            halved   6,000 kg/m, the same with the standing load halved to 99 kN/m
The support track is the split support of rounds 7 and 8 (two separate 0.4 m loops 50 mm apart, 150 x 67 mm bundles of 12 mm tape, a 1.2 m plate,
the 100 mm gap); only the plate's thickness and the requirement change.  The rule is round 5's: the lightest plate whose rating (the pull with
the worst tape at 70 percent of its limit, the winding built at 50 percent at the requirement) is at least 1.35 times the requirement.
The fast tracks' requirement is the sled's and the spacecraft's, not the casing's, so their solve is the design's (bearing_round7 --part=solve-fast).

Run (from Python/orbitalring):
  python casing_6000_round9.py --part=solve [--h=0.005] [--jobs=2]      the support track against the plate's thickness, both requirements
  python casing_6000_round9.py --part=confirm                           the picked plates (and their neighbors) at 2.5 mm cells, in full
  python casing_6000_round9.py --part=tables                            the cold load tables for the deployment runs (after confirm)
  python casing_6000_round9.py --part=run [--jobs=6]                    the deployment runs (lim_simulation.py --rising-radius --m_load=6000)
  python casing_6000_round9.py --part=report                            Table 16.2 and the sentences of 16.6
Writes runs/round9_logs/casing6000_*.json and runs/round9_casing6000/.
"""
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round7 as r7
import bearing_gap_trade as gt
import bearing_round9 as r9

G = 7.918                                  # m/s2, the casing's net weight per kilogram at the ring (bearing_round5.W_HALF)
M_NEW = 6000.0
STAND = {"kept": 198e3, "halved": 99e3}
T_GRID = (0.100, 0.090, 0.080, 0.070, 0.060, 0.050, 0.040, 0.030)
MARGIN = 1.35
RHO_FE = 7870.0; M_LIM = 1872.0
BASE = os.path.join(r7.OUT, "round9_casing6000")
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "casing6000_support.json"


def key(case, T, h):
    return f"{case}|{T:.4f}|{h:.4f}"


def job(a):
    """The split support track on a 1.2 m plate of thickness T: the requirement Fc (half the casing's weight and the standing load), the worst tape
    there, the tape at 50 percent, the rating from higher levels; in full also normal operation (half the casing's weight), the linear plate's
    square law and the loops' inductances.  bearing_gap_trade.job_sup with the plate and the requirement free."""
    case, T, h, full, Fc, Wn = a
    t0 = time.time(); theva, theva77 = be.load_theva()
    Wp = gt.SUP["Wp"]; Wl = gt.SUP["Wl"]; bw, bh = gt.SUP_BUNDLE; gap = 0.10
    out = dict(case=case, T=T, h=h, full=full, F_req=Fc, W_half=Wn, plate_m=RHO_FE * Wp * T)
    try:
        with r4.plate_patch(Wp, T):
            tr = gt.sup_track(gap, Wl, bw, bh, h)

            def both(NI):
                g, s = tr.solve(tr.NI_design(NI))
                return s["F"], (g, s)
            NI_cap, F, (g, s) = r7.find_NI(both, Fc, 437e3 * math.sqrt(Fc / 245.5e3))
            wt, _ = be.worst_tape(g, tr, tr.NI_design(NI_cap), 0.012, theva, theva77)
            tapes = NI_cap / (0.5 * wt[0]); tape = len(tr.legs) * tapes * 0.012 * bd.TAPE_T * bd.TAPE_RHO
            out["cap"] = dict(NI=NI_cap, F=F, sat=s["sat"], B_plate=s["B_plate"], ic=wt[0], B=wt[1], th=wt[2], tapes_leg=tapes, tape=tape, line=gt.face_line(g, gap, Wp, F))
            out["fill"] = tapes * 0.012 * bd.TAPE_T / (bw * bh)
            NIs, Fs, ics = [NI_cap], [F], [wt[0]]
            for k in (1.25, 1.5, 2.0):
                g2, s2 = tr.solve(tr.NI_design(k * NI_cap)); w2, _ = be.worst_tape(g2, tr, tr.NI_design(k * NI_cap), 0.012, theva, theva77)
                NIs.append(k * NI_cap); Fs.append(s2["F"]); ics.append(w2[0])
                if (k * NI_cap / tapes) / w2[0] >= 0.7:
                    break
            NI70, F70, share = gt.rating(NIs, Fs, ics, tapes)
            out["rating"] = dict(NI=NI70, F=F70, margin=F70 / Fc, levels=[dict(NI=n, F=f, ic=i) for n, f, i in zip(NIs, Fs, ics)])
            if full:
                NI_n, Fn, (g, s) = r7.find_NI(both, Wn, NI_cap * math.sqrt(Wn / Fc))
                wn, _ = be.worst_tape(g, tr, tr.NI_design(NI_n), 0.012, theva, theva77)
                out["normal"] = dict(NI=NI_n, F=Fn, sat=s["sat"], B_plate=s["B_plate"], ic=wn[0], share=(NI_n / tapes) / wn[0], line=gt.face_line(g, gap, Wp, Fn))
                gl, sl = tr.solve(tr.NI_design(1e3), linear=True)
                out["k_lin"] = sl["F"] / 1e6
                gl, sl = tr.solve(tr.NI_loop(0, 1e3), linear=True)
                out["L11"] = sl["phi"][0] / 1e3; out["L21"] = sl["phi"][1] / 1e3
        out["ok"] = True
    except Exception as e:
        out["ok"] = False; out["error"] = repr(e)
    out["seconds"] = time.time() - t0
    return out


def cache():
    try:
        return r9.load(CACHE)
    except FileNotFoundError:
        return {}


def run_pool(tasks, jobs):
    import multiprocessing as mp
    done = cache()
    todo = [t for t in tasks if not (done.get(key(t[0], t[1], t[2])) or {}).get("ok") or (t[3] and not done[key(t[0], t[1], t[2])].get("full"))]
    print(f"  {len(tasks)} points, {len(todo)} to solve, {jobs} processes", flush=True)
    if todo:
        with mp.Pool(min(jobs, len(todo))) as pool:
            for res in pool.imap_unordered(job, todo):
                done = cache(); done[key(res["case"], res["T"], res["h"])] = res; r9.save(CACHE, done)
                v = (f"NI {res['cap']['NI']/1e3:.1f} kA, tape {res['cap']['tape']:.1f} kg/m, sat {res['cap']['sat']*100:.0f} %, margin {res['rating']['margin']:.2f}" if res["ok"] else "FAILED " + res["error"])
                print(f"    {res['case']:>7} plate {res['T']*1e3:.0f} mm at {res['h']*1e3:g} mm cells: {v} ({res['seconds']:.0f} s)", flush=True)
    return cache()


def requirement(case):
    return M_NEW * G / 2 + STAND[case], M_NEW * G / 2


def tasks(h, Ts=T_GRID, full=False):
    return [(c, T, h, full) + requirement(c) for c in ("kept", "halved") for T in Ts]


def pick(done, case, h):
    """The lightest plate of the grid whose rating is at least MARGIN times the requirement."""
    P = sorted([v for v in done.values() if v.get("ok") and v["case"] == case and abs(v["h"] - h) < 1e-9], key=lambda v: v["T"])
    ok = [v for v in P if v["rating"]["margin"] == v["rating"]["margin"] and v["rating"]["margin"] >= MARGIN]
    return (ok[0] if ok else None), P


def show(done, h):
    for case in ("kept", "halved"):
        Fc, Wn = requirement(case); best, P = pick(done, case, h)
        print(f"\n  6,000 kg/m, the standing load {case}: requirement {Fc/1e3:.1f} kN/m per support track ({h*1e3:g} mm cells)")
        print(f"  {'plate mm':>8} {'kg/m':>6} {'NI kA':>6} {'sat %':>6} {'B plate T':>9} {'worst tape A':>12} {'tapes/leg':>9} {'tape kg/m':>9} {'fill %':>6} {'rating kN/m':>11} {'x req':>6} {'even field T':>12}")
        for v in P:
            c = v["cap"]
            print(f"  {v['T']*1e3:8.0f} {v['plate_m']:6.0f} {c['NI']/1e3:6.0f} {c['sat']*100:6.1f} {c['B_plate']:9.2f} {c['ic']:12.0f} {c['tapes_leg']:9,.0f} {c['tape']:9.1f} {v['fill']*100:6.0f} {v['rating']['F']/1e3:11.0f} {v['rating']['margin']:6.2f} "
                  f"{c['line']['even']:12.3f}" + ("   <- the lightest plate that meets 1.35" if best is not None and v is best else ""))


def part_solve(h, jobs):
    show(run_pool(tasks(h), jobs), h)


def part_confirm(jobs):
    done = cache(); h = 0.0025; t = []
    for case in ("kept", "halved"):
        best, P = pick(done, case, 0.005)
        Ts = sorted(set([best["T"]] + [v["T"] for v in P if abs(v["T"] - best["T"]) <= 0.0101] + [0.100]))   # the pick, its neighbors and the design's 100 mm
        t += [(case, T, h, True) + requirement(case) for T in Ts]
    show(run_pool(t, jobs), h)


def main():
    part = "solve"; h = 0.005; jobs = 2
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--h="): h = float(a[4:])
        elif a.startswith("--jobs="): jobs = int(a[7:])
    os.makedirs(r9.LOGS, exist_ok=True)
    t0 = time.time()
    if part == "solve":
        part_solve(h, jobs)
    elif part == "confirm":
        part_confirm(jobs)
    elif part == "points":                                  # --points=kept:0.08,halved:0.05 at 2.5 mm cells, in full
        pts = next(a[9:] for a in sys.argv if a.startswith("--points=")).split(",")
        show(run_pool([(c, float(T), 0.0025, True) + requirement(c) for c, T in (x.split(":") for x in pts)], jobs), 0.0025)
    else:
        import casing_6000_round9_parts as cp
        getattr(cp, "part_" + part)(jobs)
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
