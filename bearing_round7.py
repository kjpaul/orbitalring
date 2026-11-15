#!/usr/bin/env python3
"""
The bearing, round 7 (2026-11-15): Vol_III/Claude_Code_instructions_bearing_round7_2026-11-15.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 14 (the cold loop).
Run (from Python/orbitalring):
  python bearing_round7.py --part=solve-half [--gc=0,0.05,0.1]   S1   a support track as two separate loops: the 5 mm solves (background, a few minutes each)
  python bearing_round7.py --part=solve-fast                     S2.2 the accepted fast track at 2.5 mm, its tape fields and its plate's field map saved (about an hour)
  python bearing_round7.py --part=support                        S1   the statics and the change-over of every arrangement, the LIM, the casing's stiffness
  python bearing_round7.py --part=dyn                            S1   the 2-DOF heave-and-roll run of a lost support part
  python bearing_round7.py --part=heat                           S2.1 the plate heat as a sum over the plate's section
  python bearing_round7.py --part=fast                           S2.2 the accepted fast track from the saved solve
  python bearing_round7.py --part=dump                           S2.3 the hot spot of a dump at the charging voltage
Writes runs/round7_logs/.  Environment: ORBITAL_RING_OUT (default runs).  US spelling in the printed text.
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
import bearing_layout as bl
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round5 as r5
import bearing_round6b as r6
import compensation_fields as cf
import final_bearing as fb
import levitation_config as lc
import lim_hysteresis as lh

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
LOGS = os.path.join(OUT, "round7_logs")
MU0 = 4e-7 * math.pi
W_HALF = be.W_HALF; W_CASING = 2 * W_HALF
F_FULL = fb.F_LOAD_TRACK
K_FLUX = r4.K_FLUX
SUP = dict(Wp=1.2, T=0.10, Wl=0.4)                  # the support plate and loops as ruled (round 5's)
FAST = dict(Wp=1.0, T=0.20, Wl=0.6, bw=0.204, bh=0.050)   # the fast plate, loop and bundle as ruled (round 7, Section 0)
E_LAUNCH = r6.E_LAUNCH; COP = r6.COP


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else (o.tolist() if isinstance(o, np.ndarray) else float(o)))


def load(name):
    p = os.path.join(LOGS, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def find_NI(solve_F, F_target, NI1, tol=r5.CONVERGE, n_max=5):
    """The ampere-turns at which solve_F(NI) -> (F, payload) gives F_target: a secant in the logs, as round 5's convergence step."""
    Fa, pa = solve_F(NI1); NIa = NI1
    NIb = NI1 * math.sqrt(F_target / Fa); Fb, pb = solve_F(NIb); n = 0
    while abs(Fb - F_target) / F_target > tol and n < n_max:
        NIn = math.exp(math.log(NIb) + (math.log(F_target) - math.log(Fb)) * (math.log(NIb) - math.log(NIa)) / (math.log(Fb) - math.log(Fa)))
        NIa, Fa = NIb, Fb; Fb, pb = solve_F(NIn); NIb = NIn; n += 1
    return NIb, Fb, pb


# =============================================================================
# S1  a support track as two separate loops (5 mm cells)
# =============================================================================

def job_half(gc, h=0.005):
    """The support track on its 1.2 m x 100 mm plate with its two 0.4 m loops gc apart at the center (0: the bundles touching, the present track):
    the cap with both loops, the normal state with both loops, one loop alone carrying the track's normal force, and the loops' linear inductances."""
    t0 = time.time()
    theva, theva77 = be.load_theva(); r_sup, r_lift = be.load_records()
    geom = dict(bo.BEST_SUPPORT); geom["W"] = SUP["Wl"]; geom["gc"] = gc
    F_cap = r_sup["F_req"]
    base = r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"], 0.005) or r6.cached("support", SUP["Wp"], SUP["T"], SUP["Wl"])
    with r4.plate_patch(SUP["Wp"], SUP["T"]):
        tr = be.Track(geom, h, w_tape=0.012, I_ref=bo.I_REF_12, NI_ref=base["NI_req"])
        def both(NI):
            g, s = tr.solve(tr.NI_design(NI))
            return s["F"], (g, s)
        def one(NI, k=0):
            g, s = tr.solve(tr.NI_loop(k, NI))
            return s["F"], (g, s)
        NI_cap, F1, (g, s) = find_NI(both, F_cap, base["NI_req"])
        wt, _ = be.worst_tape(g, tr, tr.NI_design(NI_cap), 0.012, theva, theva77)
        tapes_leg = NI_cap / (0.5 * wt[0]); tape = len(tr.legs) * tapes_leg * 0.012 * bd.TAPE_T * bd.TAPE_RHO
        cap = dict(NI=NI_cap, F=F1, sat=s["sat"], B_plate=s["B_plate"], ic=wt[0], B=wt[1], th=wt[2], tapes_leg=tapes_leg, tape=tape, W=s["W"])
        NI0, F0, (g, s) = find_NI(both, W_HALF, NI_cap * math.sqrt(W_HALF / F_cap))
        wt0, _ = be.worst_tape(g, tr, tr.NI_design(NI0), 0.012, theva, theva77)
        normal = dict(NI=NI0, F=F0, sat=s["sat"], ic=wt0[0], share=(NI0 / tapes_leg) / wt0[0], W=s["W"], phi=s["phi"])
        NI1, F1a, (g, s) = find_NI(one, W_HALF, NI0 * 1.6)
        wt1, _ = be.worst_tape(g, tr, tr.NI_loop(0, NI1), 0.012, theva, theva77)
        alone = dict(NI=NI1, F=F1a, sat=s["sat"], B_plate=s["B_plate"], ic=wt1[0], B=wt1[1], th=wt1[2], share=(NI1 / tapes_leg) / wt1[0], W=s["W"], phi=s["phi"])
        # one loop alone carrying the whole casing (the other three loops of the two tracks gone is not a single loss; this is the cap check of one loop at 2 x W_HALF, for information)
        gl, sl = tr.solve(tr.NI_loop(0, 1e3), linear=True)
        L11, L21 = sl["phi"][0] / 1e3, sl["phi"][1] / 1e3
        # the flux-conserving jump: the survivor's ampere-turns when its neighbor's go to zero, and the force then
        NI_jump = NI0 * (1 + L21 / L11)
        Fj, (g, s) = one(NI_jump)
        jump = dict(NI=NI_jump, F=Fj)
    return dict(gc=gc, h=h, legs=[(x, sg) for x, sg in tr.legs], cap=cap, normal=normal, alone=alone, L11=L11, L21=L21, jump=jump, seconds=time.time() - t0, solves=tr.solves)


def part_solve_half(gcs=(0.0, 0.05, 0.10), jobs=3):
    import multiprocessing as mp
    hr("S1  THE SUPPORT TRACK AS TWO SEPARATE LOOPS: THE 5 mm SOLVES")
    done = load("round7_half.json") or {}
    todo = [gc for gc in gcs if f"{gc:.3f}" not in done]
    if todo:
        with mp.Pool(min(jobs, len(todo))) as pool:
            for res in pool.imap_unordered(job_half, todo):
                done = load("round7_half.json") or {}; done[f"{res['gc']:.3f}"] = res; save("round7_half.json", done)
                print(f"  gc {res['gc']*1e3:.0f} mm solved in {res['seconds']:.0f} s ({res['solves']} solves)", flush=True)
    return done


# =============================================================================
# S2.2  the accepted fast track at 2.5 mm, with its fields saved
# =============================================================================

def job_fast(h=0.0025):
    """The 1.0 m x 200 mm plate, the 0.6 m loop, the 204 x 50 mm bundle: the requirement and 1.25, 1.5 and 2 times its ampere-turns; the field at every
    model tape at each level and the plate's field map at the requirement and at half of it, saved for the parts that follow."""
    t0 = time.time()
    r_sup, r_lift = be.load_records(); F_req = r_lift["F_req"]
    name = r6.shape_name(FAST["bw"], FAST["bh"])
    geom = dict(bo.BEST_LIFT); geom["W"] = FAST["Wl"]; geom["shape"] = name
    base = r6.cached("fast", FAST["Wp"], FAST["T"], FAST["Wl"])
    out = {}
    with r4.plate_patch(FAST["Wp"], FAST["T"]):
        tr = be.Track(geom, h, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=base["NI_req"])
        def plate_map(g):
            m = g.box(-FAST["Wp"] / 2, FAST["Wp"] / 2, bd.GAP0, bd.GAP0 + FAST["T"])
            ii, jj = np.nonzero(m)
            return dict(x=g.x0 + ii * g.h, y=g.y0 + jj * g.h, Bx=g.Bx[m], By=g.By[m], mu=g.mu_r[m])
        def solve(NI):
            g, s = tr.solve(tr.NI_design(NI))
            fields = [bo.field_at_points(g, xs, ys) for xs, ys in tr.lattices]
            return s["F"], (s, fields, plate_map(g))
        NI_req, F, (s, fields, pm) = find_NI(solve, F_req, base["NI_req"])
        levels = [(1.0, NI_req, s, fields, pm)]
        for k in (1.25, 1.5, 2.0, 0.5):
            Fk, (sk, fk, pk) = solve(k * NI_req)
            levels.append((k, k * NI_req, sk, fk, pk))
            print(f"    level {k}: {Fk/1e3:.1f} kN/m, sat {sk['sat']*100:.1f} % ({(time.time()-t0)/60:.0f} min)", flush=True)
    arrs = dict(leg_x=np.array([x for x, sg in tr.legs]), leg_s=np.array([sg for x, sg in tr.legs]))
    for i, (xs, ys) in enumerate(tr.lattices):
        arrs[f"xs{i}"] = xs; arrs[f"ys{i}"] = ys
    meta = dict(h=h, F_req=F_req, NI_req=NI_req, levels=[])
    for k, NI, sk, fk, pk in levels:
        tag = f"{k:.2f}"
        for i, (Bx, By) in enumerate(fk):
            arrs[f"Bx{i}_{tag}"] = Bx; arrs[f"By{i}_{tag}"] = By
        for key, v in pk.items():
            arrs[f"plate_{key}_{tag}"] = v
        meta["levels"].append(dict(k=k, NI=NI, F=sk["F"], sat=sk["sat"], B_plate=sk["B_plate"], W=sk["W"], phi=sk["phi"]))
    meta["seconds"] = time.time() - t0
    os.makedirs(LOGS, exist_ok=True)
    np.savez_compressed(os.path.join(LOGS, "round7_fast_fields.npz"), **arrs)
    save("round7_fast_solve.json", meta)
    print(f"  the accepted fast track solved in {(time.time()-t0)/60:.0f} min: {NI_req/1e3:.0f} kA for {F/1e3:.1f} kN/m", flush=True)
    return meta


def main():
    part = "support"; gcs = (0.0, 0.05, 0.10)
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--gc="): gcs = tuple(float(x) for x in a[5:].split(","))
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    print("=" * 150 + "\n  THE BEARING, ROUND 7\n" + "=" * 150)
    t0 = time.time()
    if part == "solve-half": part_solve_half(gcs)
    elif part == "solve-fast": job_fast()
    else:
        import bearing_round7_parts as p7
        getattr(p7, "part_" + part)()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
