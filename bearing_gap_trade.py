#!/usr/bin/env python3
"""
The gap trade, the field at the plate and the flux pump (2026-11-16): Vol_III/Claude_Code_instructions_gap_trade_2026-11-16.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing told through the field at the plate).
The design is the one of rounds 7 and 8: two support tracks, each of two independent 0.4 m loops 50 mm apart at the center on a 1.2 m x 100 mm
plate, bundles 150 x 67 mm; four fast tracks, each one 0.6 m loop with 204 x 50 mm bundles on a 1.0 m x 200 mm plate; 10 mm from the tape to the gap.

Run (from Python/orbitalring):
  python bearing_gap_trade.py --part=solve-a [--h=0.005] [--gaps=0.05,0.075,0.1,0.125,0.15]   item 1 (a): the tracks as designed, only the gap changed
  python bearing_gap_trade.py --part=solve-b [--h=0.005]                                       item 1 (b): the loop width and the bundle on a grid at each gap
  python bearing_gap_trade.py --part=solve-b2                                                  item 1 (b): the picked points solved in full
  python bearing_gap_trade.py --part=solve-shape [--h=0.005]                                   item 2: the support bundle's shape at the 100 mm gap
  python bearing_gap_trade.py --part=solve-face [--h=0.0025]                                   item 3: the field across the plate's face, six states
  python bearing_gap_trade.py --part=stiff                                                     item 1: the support's pull against the gap at constant flux
  python bearing_gap_trade.py --part=cases --gap=0.075 [--noflux]                              item 1: the 114 cases at a gap, the support's constant-flux pull carried
  python bearing_gap_trade.py --part=report|pump|start                                         the tables (bearing_gap_trade_parts.py)
Writes runs/gap_trade_logs/.  Environment: ORBITAL_RING_OUT (default runs).  US spelling in the printed text.
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
LOGS = os.path.join(OUT, "gap_trade_logs")
MU0 = 4e-7 * math.pi
SUP = r7.SUP; FAST = r7.FAST; GC = 0.05
GAPS = (0.05, 0.075, 0.10, 0.125, 0.15)
SUP_BUNDLE = (0.150, 0.0667)                                   # the design's support bundle (bearing_optimize's 'wide 150x67')
W_HALF = r7.W_HALF
F_FAST = 196.5e3
MARGIN = 1.35                                                  # round 5's ruling: the rating (worst tape at 70 percent) at least 1.35 times the requirement
SUP_SHAPES = ((0.100, 0.100), (0.150, 0.0667), (0.200, 0.050), (0.250, 0.040))
SUP_WL = (0.25, 0.30, 0.35, 0.40, 0.45)
FAST_SHAPES = ((0.204, 0.050), (0.150, 0.068), (0.300, 0.034))
FAST_WL = (0.4, 0.5, 0.6, 0.7, 0.8)


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else (o.tolist() if isinstance(o, np.ndarray) else float(o)))


def load(name):
    p = os.path.join(LOGS, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def key(kind, gap, Wl, bw, bh, h):
    return f"{kind}|{gap:.3f}|{Wl:.3f}|{bw:.3f}|{bh:.4f}|{h:.4f}"


def F_cap():
    return be.load_records()[0]["F_req"]


def face_line(g, gap, Wp, F=None, arrays=False):
    """The field in the gap two cells under the plate's face, across the plate's width: the root mean square of the component normal to the face,
    of the whole field, the peak, the even field with the same pull, and the pull the normal and tangential components give by the Maxwell stress
    on that line under the plate alone (the solve's pull is taken on a box 0.1 m wider on each side and around the plate)."""
    xs, By, Bx = g.line_By(gap - 2 * g.h, -Wp / 2, Wp / 2)
    B = np.hypot(By, Bx)
    out = dict(rms_n=float(np.sqrt((By**2).mean())), rms=float(np.sqrt((B**2).mean())), peak=float(B.max()), peak_n=float(np.abs(By).max()),
               F_line=float(((By**2 - Bx**2) / (2 * MU0)).sum() * g.h), y=gap - 2 * g.h)
    if F is not None:
        out["even"] = math.sqrt(2 * MU0 * F / Wp)
    if arrays:
        out.update(x=xs.tolist(), By=By.tolist(), Bx=Bx.tolist())
    return out


def sup_track(gap, Wl, bw, bh, h):
    geom = dict(bo.BEST_SUPPORT); geom["W"] = Wl; geom["gc"] = GC; geom["shape"] = r6.shape_name(bw, bh)
    return be.Track(geom, h, gap=gap, w_tape=0.012, I_ref=bo.I_REF_12, NI_ref=437e3)


def fast_track(gap, Wl, bw, bh, h):
    geom = dict(bo.BEST_LIFT); geom["W"] = Wl; geom["shape"] = r6.shape_name(bw, bh)
    return be.Track(geom, h, gap=gap, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=452e3)


def rating(NIs, Fs, ics, tapes):
    NIs, Fs, ics = np.array(NIs), np.array(Fs), np.array(ics)
    o = np.argsort(NIs); NIs, Fs, ics = NIs[o], Fs[o], ics[o]
    share = (NIs / tapes) / ics
    if share.max() < 0.7:
        return float("nan"), float("nan"), share.tolist()
    NI70 = float(np.interp(0.7, share, NIs))
    return NI70, math.exp(float(np.interp(math.log(NI70), np.log(NIs), np.log(Fs)))), share.tolist()


def job_sup(a):
    """One support track (two separate loops of Wl, GC apart at the center, four bundles bw x bh of 12 mm tape, the 1.2 m x 100 mm plate) at a gap:
    the cap (245.5 kN/m) with the worst tape, the tape at 50 percent, the rating from two higher levels; in full also normal operation (47.5 kN/m),
    the square law of the linear plate and the loops' inductances."""
    gap, Wl, bw, bh, h, full = a
    t0 = time.time(); theva, theva77 = be.load_theva(); Fc = F_cap()
    out = dict(kind="support", gap=gap, Wl=Wl, bw=bw, bh=bh, h=h, full=full, width=2 * Wl + 2 * bw + GC)
    try:
        with r4.plate_patch(SUP["Wp"], SUP["T"]):
            tr = sup_track(gap, Wl, bw, bh, h)
            def both(NI):
                g, s = tr.solve(tr.NI_design(NI))
                return s["F"], (g, s)
            NI_cap, F, (g, s) = r7.find_NI(both, Fc, 437e3 * (gap / 0.10)**0.8)
            wt, _ = be.worst_tape(g, tr, tr.NI_design(NI_cap), 0.012, theva, theva77)
            tapes = NI_cap / (0.5 * wt[0]); tape = len(tr.legs) * tapes * 0.012 * bd.TAPE_T * bd.TAPE_RHO
            ys = np.concatenate([y for x, y in tr.lattices])
            out["cap"] = dict(NI=NI_cap, F=F, sat=s["sat"], B_plate=s["B_plate"], ic=wt[0], B=wt[1], th=wt[2], face=wt[3], x=wt[4], tapes_leg=tapes, tape=tape, W=s["W"],
                              phi=s["phi"], line=face_line(g, gap, SUP["Wp"], F))
            out["tape_near"] = gap - float(ys.max()); out["tape_far"] = gap - float(ys.min())
            out["fill"] = tapes * 0.012 * bd.TAPE_T / (bw * bh)
            NIs, Fs, ics = [NI_cap], [F], [wt[0]]
            for k in (1.25, 1.5):
                g2, s2 = tr.solve(tr.NI_design(k * NI_cap)); w2, _ = be.worst_tape(g2, tr, tr.NI_design(k * NI_cap), 0.012, theva, theva77)
                NIs.append(k * NI_cap); Fs.append(s2["F"]); ics.append(w2[0])
                if (k * NI_cap / tapes) / w2[0] >= 0.7:
                    break
            NI70, F70, share = rating(NIs, Fs, ics, tapes)
            out["rating"] = dict(NI=NI70, F=F70, margin=F70 / Fc, levels=[dict(NI=n, F=f, ic=i) for n, f, i in zip(NIs, Fs, ics)])
            if full:
                NI_n, Fn, (g, s) = r7.find_NI(both, W_HALF, 0.36 * NI_cap)
                wn, _ = be.worst_tape(g, tr, tr.NI_design(NI_n), 0.012, theva, theva77)
                out["normal"] = dict(NI=NI_n, F=Fn, sat=s["sat"], B_plate=s["B_plate"], ic=wn[0], B=wn[1], th=wn[2], share=(NI_n / tapes) / wn[0], W=s["W"], phi=s["phi"],
                                     line=face_line(g, gap, SUP["Wp"], Fn))
                gl, sl = tr.solve(tr.NI_design(1e3), linear=True)
                out["k_lin"] = sl["F"] / 1e6
                gl, sl = tr.solve(tr.NI_loop(0, 1e3), linear=True)
                out["L11"] = sl["phi"][0] / 1e3; out["L21"] = sl["phi"][1] / 1e3
        out["ok"] = True
    except Exception as e:
        out["ok"] = False; out["error"] = repr(e)
    out["seconds"] = time.time() - t0
    return out


def job_fast(a):
    """One fast track (one loop of Wl, two bundles bw x bh of 3 mm tape, the 1.0 m x 200 mm plate) at a gap: the requirement (196.5 kN/m), the field at
    every model tape there and at 1.25 and 1.5 times the ampere-turns; for each way of laying the stacks the worst tape, the tape at 50 percent, the
    rating; the stored energy; in full also the tape's pulse loss of the heaviest launch."""
    import bearing_round7_parts as p7
    gap, Wl, bw, bh, h, full = a
    t0 = time.time(); theva, theva77 = be.load_theva()
    out = dict(kind="fast", gap=gap, Wl=Wl, bw=bw, bh=bh, h=h, full=full, width=Wl + bw)
    try:
        with r4.plate_patch(FAST["Wp"], FAST["T"]):
            tr = fast_track(gap, Wl, bw, bh, h)
            def solve(NI):
                g, s = tr.solve(tr.NI_design(NI))
                return s["F"], (g, s, [bo.field_at_points(g, xs, ys) for xs, ys in tr.lattices])
            NI_req, F, (g, s, f1) = r7.find_NI(solve, F_FAST, 452e3 * (gap / 0.10)**0.8)
            line = face_line(g, gap, FAST["Wp"], F)
            levels = [(NI_req, F, f1)]
            for k in (1.25, 1.5):
                Fk, (gk, sk, fk) = solve(k * NI_req); levels.append((k * NI_req, Fk, fk))
        z = dict(leg_x=np.array([x for x, sg in tr.legs]))
        for i, (xs, ys) in enumerate(tr.lattices):
            z[f"xs{i}"] = xs; z[f"ys{i}"] = ys
        design_bundle = abs(bw - FAST["bw"]) < 1e-9 and abs(bh - FAST["bh"]) < 1e-9
        ways = {"each stack the better way": lambda f: (f, "best"), "all faces horizontal": lambda f: (f, "h")}
        if design_bundle:
            ways["cells turned as ruled"] = lambda f: (p7.turned(z, f, (p7.PATTERN_NEG, p7.PATTERN_POS))[0], "h")
            ways["each cell the better way"] = lambda f: (p7.turned(z, f, None)[0], "h")
        W = s["W"]; Phi = 2 * W / NI_req
        out.update(NI=NI_req, F=F, sat=s["sat"], B_plate=s["B_plate"], W=W, E_coil=W * 25.0 / 20, V=10 * Phi * 25.0 / 10e-3, I=NI_req / 200.0, line=line, ways={})
        ys = np.concatenate([y for x, y in tr.lattices]); out["tape_near"] = gap - float(ys.max()); out["tape_far"] = gap - float(ys.min())
        li = dict(r6.launch_profile(W, F_FAST), F_req=F_FAST) if full else None
        for lab, fn in ways.items():
            wi = [r6.worst_ic(fn(fk)[0], fn(fk)[1], theva) for NI, Fk, fk in levels]
            ic, B, th = wi[0]
            tapes = NI_req / (0.5 * ic); tape = 2 * tapes * 0.003 * bd.TAPE_T * bd.TAPE_RHO
            NI70, F70, share = rating([l[0] for l in levels], [l[1] for l in levels], [w[0] for w in wi], tapes)
            d = dict(ic=ic, B=B, th=th, tapes_leg=tapes, tape=tape, NI70=NI70, F70=F70, margin=F70 / F_FAST, layers=(NI_req / 20) / (0.5 * ic) / 10, fill=tapes * 0.003 * bd.TAPE_T / (bw * bh))
            if full and lab in ("cells turned as ruled", "each stack the better way"):
                fm, face = fn(f1)
                d["Ep"], d["Es"] = r6.tape_loss_launch(fm, face, tapes, d["layers"], li)
            out["ways"][lab] = d
        out["ok"] = True
    except Exception as e:
        out["ok"] = False; out["error"] = repr(e)
    out["seconds"] = time.time() - t0
    return out


def _run(a):
    return job_sup(a[1:]) if a[0] == "support" else job_fast(a[1:])


def run_pool(tasks, fname, jobs=6):
    """Solve the tasks that the cache does not hold, in parallel; each result is written as it comes."""
    import multiprocessing as mp
    done = load(fname) or {}
    todo = [t for t in tasks if not (done.get(key(*t[:6])) or {}).get("ok") or (t[6] and not done[key(*t[:6])].get("full"))]
    print(f"  {len(tasks)} points, {len(todo)} to solve, {jobs} processes", flush=True)
    t0 = time.time()
    if todo:
        with mp.Pool(min(jobs, len(todo))) as pool:
            for res in pool.imap_unordered(_run, todo):
                done = load(fname) or {}; done[key(res["kind"], res["gap"], res["Wl"], res["bw"], res["bh"], res["h"])] = res; save(fname, done)
                if res["ok"]:
                    v = (f"cap {res['cap']['NI']/1e3:.1f} kA, tape {res['cap']['tape']:.1f} kg/m, margin {res['rating']['margin']:.2f}" if res["kind"] == "support" else
                         f"{res['NI']/1e3:.1f} kA, tape {res['ways']['each stack the better way']['tape']:.1f} kg/m, margin {res['ways']['each stack the better way']['margin']:.2f}")
                else:
                    v = "FAILED " + res["error"]
                print(f"    {res['kind']:>7} gap {res['gap']*1e3:.0f} mm, loop {res['Wl']:.2f} m, bundle {res['bw']*1e3:.0f} x {res['bh']*1e3:.0f}: {v} ({res['seconds']:.0f} s; {(time.time()-t0)/60:.1f} min)", flush=True)
    return load(fname) or {}


def tasks_a(h, gaps=GAPS):
    return [("support", g, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], h, True) for g in gaps] + [("fast", g, FAST["Wl"], FAST["bw"], FAST["bh"], h, True) for g in gaps]


def tasks_b(h, gaps=GAPS):
    t = []
    for g in gaps:
        for bw, bh in SUP_SHAPES:
            for Wl in SUP_WL:
                if 2 * Wl + 2 * bw + GC <= SUP["Wp"] + 1e-9 and Wl > bw + 1e-9:
                    t.append(("support", g, Wl, bw, bh, h, False))
        for bw, bh in FAST_SHAPES:
            for Wl in FAST_WL:
                if Wl + bw <= FAST["Wp"] + 1e-9 and Wl > bw + 1e-9:
                    t.append(("fast", g, Wl, bw, bh, h, False))
    return t


def picks(h, gaps=GAPS):
    """Item 1 (b)'s rule: at each gap, of the grid points whose rating is at least 1.35 times the requirement (round 5's ruling) and whose track fits
    its plate, the one with the least tape; the fast track with each stack laid the better way."""
    done = load("gap_trade_b.json") or {}; a = load("gap_trade_a.json") or {}
    out = {}
    for kind in ("support", "fast"):
        for g in gaps:
            P = [v for v in list(done.values()) + list(a.values()) if v.get("ok") and v["kind"] == kind and abs(v["gap"] - g) < 1e-9 and abs(v["h"] - h) < 1e-9]
            def tape(v): return v["cap"]["tape"] if kind == "support" else v["ways"]["each stack the better way"]["tape"]
            def margin(v): return v["rating"]["margin"] if kind == "support" else v["ways"]["each stack the better way"]["margin"]
            ok = [v for v in P if margin(v) == margin(v) and margin(v) >= MARGIN]
            if ok:
                out[(kind, g)] = min(ok, key=tape)
    return out


def part_solve_shape(h):
    tasks = [("support", 0.10, SUP["Wl"], bw, bh, h, True) for bw, bh in SUP_SHAPES]
    return run_pool(tasks, "gap_trade_shape.json", jobs=4)


FACE_STATES = ("support normal", "support cap", "support rating", "support one loop", "fast 196.5", "fast rating")


def job_face(a):
    """Item 3: one state of a track at the design's 100 mm gap; the field on the line two cells under the plate's face is saved."""
    state, h = a
    import bearing_round8 as r8
    t0 = time.time()
    half = r7.load("round7_half_2p5.json"); rows = r8.load("round8_rows.json"); fsol = r7.load("round7_fast_solve.json")
    gap = bd.GAP0
    if state.startswith("support"):
        NI = {"support normal": half["normal"]["NI"], "support cap": half["cap"]["NI"], "support rating": rows["rating_support"]["NI"], "support one loop": half["alone"]["NI"]}[state]
        with r4.plate_patch(SUP["Wp"], SUP["T"]):
            tr = sup_track(gap, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], h)
            legs = tr.NI_loop(0, NI) if state == "support one loop" else tr.NI_design(NI)
            g, s = tr.solve(legs)
        Wp = SUP["Wp"]
    else:
        NI = fsol["NI_req"] if state == "fast 196.5" else next(r["NI"] for r in rows["t74_fast"] if r["state"].startswith("the rating"))
        with r4.plate_patch(FAST["Wp"], FAST["T"]):
            tr = fast_track(gap, FAST["Wl"], FAST["bw"], FAST["bh"], h)
            g, s = tr.solve(tr.NI_design(NI))
        Wp = FAST["Wp"]
    out = dict(state=state, h=h, NI=NI, F=s["F"], sat=s["sat"], B_plate=s["B_plate"], Wp=Wp, legs=[(x, sg) for x, sg in tr.legs], bw=tr.wb, bh=tr.hb,
               line=face_line(g, gap, Wp, s["F"], arrays=True), seconds=time.time() - t0)
    # the same line 0.1 m wider on each side (the solve's force box), for the pull that falls outside the plate's width
    xs, By, Bx = g.line_By(gap - 2 * g.h, -Wp / 2 - 0.1, Wp / 2 + 0.1)
    out["F_line_box"] = float(((By**2 - Bx**2) / (2 * MU0)).sum() * g.h)
    return out


def part_solve_face(h):
    import multiprocessing as mp
    done = load("gap_trade_face.json") or {}
    todo = [(s, h) for s in FACE_STATES if s not in done]
    if todo:
        with mp.Pool(min(4, len(todo))) as pool:
            for res in pool.imap_unordered(job_face, todo):
                done = load("gap_trade_face.json") or {}; done[res["state"]] = res; save("gap_trade_face.json", done)
                l = res["line"]
                print(f"    {res['state']}: {res['NI']/1e3:.1f} kA, {res['F']/1e3:.1f} kN/m; face rms {l['rms_n']:.3f} T (normal), even {l['even']:.3f} T, peak {l['peak']:.2f} T; in the plate {res['B_plate']:.2f} T ({res['seconds']:.0f} s)", flush=True)
    return done


def part_stiff(h=0.005):
    """The support track as designed, linear plate: the loops' flux per ampere-turn and the pull per ampere-turn squared at gaps from 10 to 270 mm, from
    which the pull at constant flux against the gap; then at each gap of the trade, with the plate saturating and the solved normal ampere-turns, the
    stiffness at constant current and at constant flux from four solves (the gap 10 mm either way, the ampere-turns 5 percent either way)."""
    t0 = time.time()
    a = load("gap_trade_a.json") or {}
    out = load("gap_trade_stiff.json") or {}
    if "sweep" not in out:
        rows = []
        with r4.plate_patch(SUP["Wp"], SUP["T"]):
            for gp in np.arange(0.01, 0.2751, 0.01):
                tr = sup_track(float(gp), SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], h)
                g, s = tr.solve(tr.NI_design(1e3), linear=True)
                rows.append(dict(gap=float(gp), Lp=s["phi"][0] / 1e3, Fc=s["F"] / 1e6))
        out["sweep"] = rows; save("gap_trade_stiff.json", out)
        print(f"    linear sweep done ({(time.time()-t0)/60:.1f} min)", flush=True)
    out.setdefault("local", {})
    for g0 in GAPS:
        if f"{g0:.3f}" in out["local"]:
            continue
        ra = a.get(key("support", g0, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], h))
        if not ra or not ra.get("ok") or "normal" not in ra:
            print(f"    gap {g0*1e3:.0f} mm: no solve of item 1 (a) yet, skipped"); continue
        NI = ra["normal"]["NI"]; dg = 0.01; e = 0.05
        with r4.plate_patch(SUP["Wp"], SUP["T"]):
            def at(gp, n):
                tr = sup_track(gp, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], h)
                g, s = tr.solve(tr.NI_design(n))
                return s["F"], float(np.mean(np.abs(s["phi"])))
            Fm, pm = at(g0 - dg, NI); Fp, pp = at(g0 + dg, NI); Fl, pl = at(g0, NI * (1 - e)); Fh, ph = at(g0, NI * (1 + e))
        dFdg = (Fp - Fm) / (2 * dg); dpdg = (pp - pm) / (2 * dg); dFdN = (Fh - Fl) / (2 * e * NI); dpdN = (ph - pl) / (2 * e * NI)
        k_phi = dFdg - dFdN * dpdg / dpdN
        out["local"][f"{g0:.3f}"] = dict(gap=g0, NI=NI, k_I=dFdg, k_phi=k_phi, dFdN=dFdN, dpdg=dpdg, dpdN=dpdN)
        save("gap_trade_stiff.json", out)
        print(f"    gap {g0*1e3:.0f} mm: at constant current {dFdg/1e6:.3f}, at constant flux {k_phi/1e6:.3f} kN/m per mm per track ({(time.time()-t0)/60:.1f} min)", flush=True)
    return out


def sweep_fit():
    """The linear sweep's pull per ampere-turn squared and flux per ampere-turn against the gap, as fourth-order fits of their logarithms.  The raw
    sweep has steps of about 3 percent where the plate's face crosses a cell edge of the nested grids; the fit takes them out."""
    st = load("gap_trade_stiff.json")
    gs = np.array([r["gap"] for r in st["sweep"]])
    pf = np.polyfit(gs, np.log([r["Fc"] for r in st["sweep"]]), 4); pl = np.polyfit(gs, np.log([r["Lp"] for r in st["sweep"]]), 4)
    return gs, pf, pl


def fit_stiffness(g0, F):
    """From the fit: the stiffness of a track pulling F at the gap g0, at constant current and at constant flux (N/m per m, negative: unstable)."""
    gs, pf, pl = sweep_fit()
    dF = float(np.polyval(np.polyder(pf), g0)); dL = float(np.polyval(np.polyder(pl), g0))
    return dF * F, (dF - 2 * dL) * F


def flux_table(g0):
    """The support's pull at constant flux over its pull at the nominal gap g0, against the gap (the linear plate's sweep, fitted)."""
    gs, pf, pl = sweep_fit()
    g = np.arange(gs[0], gs[-1] + 1e-9, 0.0025)
    lnq = np.polyval(pf, g) - 2 * np.polyval(pl, g); q = np.exp(lnq - (np.polyval(pf, g0) - 2 * np.polyval(pl, g0)))
    return np.concatenate([[0.0], g]), np.concatenate([[q[0]], q])


def part_cases(gap, flux=True):
    """Round 7's 114 cases (final_bearing's Part B through bearing_layout.run_cases: the 12,000 kg/m casing, 25 m fast coils in 20 families, the 10 ms
    rise, the 10 mm protection in every lane, the support frozen, the split support's coupling) with the nominal gap set to gap.  flux: in the point
    lanes the frozen support's pull follows the gap at constant flux (sled_passage_pulsed.SUPPORT_OF_GAP); without it the support's pull does not
    depend on the gap, as in every earlier round."""
    import bearing_layout as bl
    import bearing_round6a as a6
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    import compensation_tracks as ct
    cp = r7.load("round7_coupling.json")
    old = pz.GAP0, bd.GAP0, sp.GAP0, be.GAP, bl.GAP, ct.GAP0
    pz.GAP0 = bd.GAP0 = sp.GAP0 = be.GAP = bl.GAP = ct.GAP0 = gap
    pz.SUPPORT_OF_GAP = flux_table(gap) if flux else None
    tag = f"{gap*1e3:.0f}mm" + ("_flux" if flux else "")
    out_dir = os.path.join("gap_trade_logs", f"gap_trade_cases_{tag}.txt")
    try:
        A2, theva, sup, lift, tabs = a6.casing_setup()
        t0 = time.time()
        gaps, records, txt = bl.run_cases(a6.TAU, True, A2, theva, sup, lift, tabs, out_dir, prot_all=True,
                                          head=f"[gap trade: nominal gap {gap*1e3:.0f} mm, {'the support at constant flux' if flux else 'the support pull independent of the gap'}, coupling {cp['c_split']*100:.2f} percent]",
                                          couple=cp["c_split"])
    finally:
        pz.GAP0, bd.GAP0, sp.GAP0, be.GAP, bl.GAP, ct.GAP0 = old; pz.SUPPORT_OF_GAP = None
    rows = [dict(label=g[0], model=str(g[1]), gap_b=float(g[2]), gap_t=float(g[3]), contact=bool(g[4]), single=bool(bl.is_single(g[0]))) for g in gaps]
    save(f"gap_trade_cases_{tag}.json", dict(gap=gap, flux=flux, table="fit", rows=rows, summary=bl.parse_summary(txt), seconds=time.time() - t0))
    s = bl.parse_summary(txt)
    print(f"  gap {gap*1e3:.0f} mm, flux {flux}: {len(rows)} cases, smallest single-failure gap {s['gap']} mm in '{s['case']}', contacts: {s['contacts']} ({time.time()-t0:.0f} s)", flush=True)
    return rows


def main():
    part = "report"; h = None; gaps = GAPS; gap = 0.10; flux = True
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--h="): h = float(a[4:])
        elif a.startswith("--gaps="): gaps = tuple(float(x) for x in a[7:].split(","))
        elif a.startswith("--gap="): gap = float(a[6:])
        elif a == "--noflux": flux = False
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    t0 = time.time()
    if part == "solve-a":
        run_pool(tasks_a(h or 0.005, gaps), "gap_trade_a.json", jobs=5)
    elif part == "solve-b":
        run_pool(tasks_b(h or 0.005, gaps), "gap_trade_b.json", jobs=6)
    elif part == "solve-b2":
        pk = picks(h or 0.005, gaps)
        run_pool([(k[0], v["gap"], v["Wl"], v["bw"], v["bh"], v["h"], True) for k, v in pk.items()], "gap_trade_b.json", jobs=5)
    elif part == "solve-shape":
        part_solve_shape(h or 0.005)
    elif part == "solve-face":
        part_solve_face(h or 0.0025)
    elif part == "stiff":
        part_stiff(h or 0.005)
    elif part == "cases":
        part_cases(gap, flux)
    else:
        import bearing_gap_trade_parts as pg
        getattr(pg, "part_" + part)()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
