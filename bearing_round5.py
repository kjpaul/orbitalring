#!/usr/bin/env python3
"""
The bearing, round 5 (2026-11-13): Vol_III/Claude_Code_instructions_bearing_round5_2026-11-13.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 8 (the stores), Chapter 14 (the cold load).
Run: python bearing_round5.py --part=<part> [--quick] [--jobs=N]
  sweep     S2  the plate search: every (plate width, thickness, loop width) point of a grid solved at a cell size, in parallel, cached
            in runs/bearing_round5_plates.json  (--h=0.01|0.005|0.0025, --kind=support|fast|both, --grid=coarse|refine|final)
  pick      S2  the optimum of the cached points for each margin (1.2, 1.35, 1.5) and for the 8.6 kg-of-cable objective
  final     S2  the chosen points at 2.5 mm with the full curve, the coupling at the resulting spacing, the plate heat of a launch
  stores    S3  the pulse power of the pulled zone, the three ways to feed it, capacitors against SMES, the local ride-through store
  switches  S4  what the switch does, the mass per meter at 10 ms for 10 and 14 turns from the datasheets read
  tapeloss  S5  the pulse loss of the soldered stacks by the field's angle on the optimum plates at 10 ms
  coldload  S6  the leads' idle heat, the cold load in operation
  cases     S7  the 10 ms cases on the book's casing with the optimum plates' caps
  charging  S8  the support coil's two closing orders
  all       everything but sweep (which is run first by hand, see the brief)
Environment: ORBITAL_RING_OUT (default runs).  The ring: ring_altitude's default (the rising quadratic, 12,000 kg/m casing).

Words as Paul ruled them: a tape; a layer (one tape in a soldered stack); a turn (a stack of soldered layers once around the loop);
a coil (its series turns); the coils that overlap in a bundle.  US spelling in the printed text (meter).
"""
import contextlib
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_optimize as bo
import bearing_electrical as be
import bearing_design as bd
import bearing_layout as bl
import bearing_stageB as sb
import bearing_round4 as r4
import final_bearing as fb
import sled_passage_pulsed as pz
import magnetostatic_2d as ms
import compensation_fields as cf
import lim_hysteresis as lh

OUT = ring.OUT_DIR
CACHE = os.path.join(OUT, "bearing_round5_plates.json")
L_SUP, L_FAST, N_FAM, N_SUP = 1000.0, 25.0, 20, 10
TAU = 10e-3                                      # Paul's ruling: the fast coils' rise is 10 ms
TURNS = (10, 14)                                 # the two series-turn counts carried
MARGINS = (1.2, 1.35, 1.5)                       # the rating margin: F70 / F_req; 1.35 is the ruling
CABLE_PER_CASING_KG = 8.6                        # a kilogram on the casing counted as 8.6 kg of cable, for information
LAUNCHES_PER_DAY = 500
Q_STANDING_W_PER_M = 29.0                        # Chapter 14's standing cold load, 14.5 kW per 500 m site

SPACING_MIN = 0.2                                # m between a support plate's edge and the load plate's beside it (Stage B S2.2)

# the grids of the search: plate width (m), plate thickness (m), loop width (m)
GRIDS = {
    "coarse": {
        "support": dict(Wp=(1.0, 1.2, 1.4, 1.6), T=(0.075, 0.10, 0.125, 0.15, 0.20), Wl=(0.4, 0.5, 0.6)),
        "fast": dict(Wp=(0.8, 1.0, 1.2, 1.4), T=(0.075, 0.10, 0.15, 0.20), Wl=(0.6, 0.8, 1.0)),
    },
}


def hr(t=""):
    print("\n" + "=" * 150)
    if t:
        print("  " + t)
        print("=" * 150)


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(c):
    os.makedirs(OUT, exist_ok=True)
    with open(CACHE, "w", encoding="utf-8") as f:
        json.dump(c, f, indent=1)


def key_of(kind, Wp, T, Wl, h):
    return f"{kind}|{Wp:.3f}|{T:.4f}|{Wl:.3f}|{h:.4f}"


# =============================================================================
# S2  one plate point
# =============================================================================

def solve_point(kind, Wp, T, Wl, h, ks=(1.5,), verbose=False):
    """One track (support: two loops of Wl in four 150 x 67 mm bundles of 12 mm tape; fast: one loop of Wl in two 200 x 50 mm bundles of
    3 mm tape) on a plate Wp x T, Part 3b's method as bearing_round4.solve_plate_case: the ampere-turns for the requirement, the worst tape,
    the turns at the 50 percent share, the tape; then the worst tape at ks x NI_req, from which the rating (70 percent) and the cap are
    interpolated.  Returns a plain dict (JSON-able)."""
    r_sup, r_lift = be.load_records(); theva, theva77 = be.load_theva()
    if kind == "support":
        r, F_req, w_tape, base = r_sup, r_sup["F_req"], 0.012, bo.BEST_SUPPORT
    else:
        r, F_req, w_tape, base = r_lift, r_lift["F_req"], 0.003, bo.BEST_LIFT
    geom = dict(base); geom["W"] = Wl
    I_ref = bo.I_REF_12 if w_tape == 0.012 else bo.I_REF_3
    t0 = time.time()
    with r4.plate_patch(Wp, T):
        tr = be.Track(geom, h, w_tape=w_tape, I_ref=I_ref, NI_ref=r["NI_req"])
        NI1 = r["NI_req"]
        g1, s1 = tr.solve(tr.NI_design(NI1))
        NI2 = NI1 * math.sqrt(F_req / s1["F"])
        g2, s2 = tr.solve(tr.NI_design(NI2))
        NI_req = math.exp(math.log(NI1) + (math.log(F_req) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI1)) / (math.log(s2["F"]) - math.log(s1["F"])))
        g3, s3 = tr.solve(tr.NI_design(NI_req))
        wt, lim = be.worst_tape(g3, tr, tr.NI_design(NI_req), w_tape, theva, theva77)
        ic_req = wt[0]; I50 = 0.5 * ic_req
        n_legs = len(tr.legs); turns_leg = NI_req / I50
        tape = n_legs * turns_leg * w_tape * bd.TAPE_T * bd.TAPE_RHO
        curve = [(NI_req, s3["F"], s3["sat"], s3["B_plate"], ic_req, s3["W"])]
        for k in ks:
            NIk = k * NI_req
            gk, sk = tr.solve(tr.NI_design(NIk))
            wtk, _ = be.worst_tape(gk, tr, tr.NI_design(NIk), w_tape, theva, theva77)
            curve.append((NIk, sk["F"], sk["sat"], sk["B_plate"], wtk[0], sk["W"]))
    NIs = np.array([c[0] for c in curve]); Fs = np.array([c[1] for c in curve]); ics = np.array([c[4] for c in curve])
    share = (NIs / turns_leg) / ics
    def at_share(s_target):
        if share.max() < s_target:
            return float("nan"), float("nan")
        NI_s = float(np.interp(s_target, share, NIs))
        F_s = math.exp(float(np.interp(math.log(NI_s), np.log(NIs), np.log(Fs))))
        return NI_s, F_s
    NI70, F70 = at_share(0.70); NI100, F100 = at_share(1.0)
    res = dict(kind=kind, Wp=Wp, T=T, Wl=Wl, h=h, NI_req=NI_req, F_req=F_req, F_solved=s3["F"], ic_req=ic_req, B=wt[1], th=wt[2], face=wt[3],
               turns_leg=turns_leg, n_legs=n_legs, tape=tape, sat=s3["sat"], B_plate=s3["B_plate"], W_mag=s3["W"], plate_m=bd.RHO_PLATE * Wp * T,
               curve=[list(c) for c in curve], share=share.tolist(), NI70=NI70, F70=F70, NI100=NI100, F100=F100, F_last=float(Fs[-1]),
               margin=F70 / F_req if F70 == F70 else float("nan"), seconds=time.time() - t0, ks=list(ks))
    if verbose:
        print(f"    {kind:>7} plate {Wp:.2f} m x {T*1e3:.0f} mm, loop {Wl:.2f} m, h {h*1e3:g} mm: {NI_req/1e3:.0f} kA for {F_req/1e3:.1f} kN/m (solved {s3['F']/1e3:.1f}); worst {w_tape*1e3:.0f} mm tape {wt[1]:.2f} T at {wt[2]:.0f} deg, {ic_req:.0f} A;")
        print(f"            {turns_leg:,.0f} turns per leg, {tape:.1f} kg/m; plate {s3['sat']*100:.1f} % saturated, {s3['B_plate']:.2f} T, {bd.RHO_PLATE*Wp*T:,.0f} kg/m; rating {F70/1e3:.0f} kN/m ({F70/F_req:.2f} x), cap {F100/1e3:.0f} ({time.time()-t0:.0f} s)", flush=True)
    return res


def _worker(args):
    kind, Wp, T, Wl, h, ks = args
    try:
        res = solve_point(kind, Wp, T, Wl, h, ks=tuple(ks))
        res["ok"] = True
    except Exception as e:                                               # a point that fails to converge is recorded, not fatal
        res = dict(kind=kind, Wp=Wp, T=T, Wl=Wl, h=h, ok=False, error=repr(e))
    return res


def part_sweep(h=0.01, kind="both", grid="coarse", jobs=6, points=None, ks=(1.5,)):
    """Solve every point of the grid (or the given points) at cell h that is not yet in the cache, in parallel."""
    import multiprocessing as mp
    hr(f"S2  THE PLATE SEARCH: {grid} grid at {h*1e3:g} mm cells, {kind}")
    cache = load_cache()
    todo = []
    kinds = ("support", "fast") if kind == "both" else (kind,)
    if points is None:
        points = []
        for kd in kinds:
            g = GRIDS[grid][kd]
            for Wp in g["Wp"]:
                for T in g["T"]:
                    for Wl in g["Wl"]:
                        points.append((kd, Wp, T, Wl))
    for kd, Wp, T, Wl in points:
        k = key_of(kd, Wp, T, Wl, h)
        if k in cache and cache[k].get("ok"):
            continue
        todo.append((kd, Wp, T, Wl, h, list(ks)))
    print(f"  {len(points)} points, {len(todo)} to solve, {jobs} processes", flush=True)
    t0 = time.time()
    if todo:
        with mp.Pool(jobs) as pool:
            for res in pool.imap_unordered(_worker, todo):
                k = key_of(res["kind"], res["Wp"], res["T"], res["Wl"], res["h"])
                cache = load_cache(); cache[k] = res; save_cache(cache)
                if res.get("ok"):
                    print(f"  {res['kind']:>7} {res['Wp']:.2f} x {res['T']*1e3:3.0f} mm, loop {res['Wl']:.2f}: NI {res['NI_req']/1e3:4.0f} kA, tape {res['tape']:6.1f} kg/m, sat {res['sat']*100:4.1f} %, "
                          f"rating {res['F70']/1e3:4.0f} kN/m = {res['margin']:.2f} x, plate {res['plate_m']:6,.0f} kg/m  ({res['seconds']:.0f} s; {(time.time()-t0)/60:.1f} min)", flush=True)
                else:
                    print(f"  {res['kind']:>7} {res['Wp']:.2f} x {res['T']*1e3:3.0f} mm, loop {res['Wl']:.2f}: FAILED {res['error']}", flush=True)
    print(f"  sweep done in {(time.time()-t0)/60:.1f} min", flush=True)


def part_pick(h=0.01, verbose=True):
    """From the cached points at cell h: for each kind and margin the lightest plate that meets it, and the point that is best when a kilogram
    on the casing is counted as 8.6 kg of cable."""
    hr(f"S2  THE OPTIMUM FROM THE CACHED POINTS AT {h*1e3:g} mm")
    cache = load_cache()
    pts = [v for v in cache.values() if v.get("ok") and abs(v["h"] - h) < 1e-6]
    picks = {}
    for kind in ("support", "fast"):
        P = [p for p in pts if p["kind"] == kind]
        if not P:
            continue
        print(f"\n  {kind}: {len(P)} points.  plate W x T, loop W | NI_req kA | worst tape A, T, deg | turns/leg | tape kg/m | sat % | rating kN/m, x | cap | plate kg/m | plate + 8.6 x 2 x tape")
        for p in sorted(P, key=lambda p: (p["Wl"], p["Wp"], p["T"])):
            alt = p["plate_m"] + CABLE_PER_CASING_KG * 2 * p["tape"]
            print(f"    {p['Wp']:.2f} x {p['T']*1e3:3.0f} mm, loop {p['Wl']:.2f} | {p['NI_req']/1e3:4.0f} | {p['ic_req']:4.0f} A {p['B']:.2f} T {p['th']:3.0f} | {p['turns_leg']:7,.0f} | {p['tape']:6.1f} | {p['sat']*100:4.1f} | "
                  f"{p['F70']/1e3:4.0f} {p['margin']:.2f} | {p['F100']/1e3:4.0f} | {p['plate_m']:6,.0f} | {alt:6,.0f}")
        picks[kind] = {}
        for m in MARGINS:
            ok = [p for p in P if p["margin"] == p["margin"] and p["margin"] >= m]
            if ok:
                best = min(ok, key=lambda p: p["plate_m"])
                alt = min(ok, key=lambda p: p["plate_m"] + CABLE_PER_CASING_KG * 2 * p["tape"])
                picks[kind][str(m)] = dict(best=best, alt=alt)
                print(f"  margin {m:.2f}: lightest plate {best['Wp']:.2f} x {best['T']*1e3:.0f} mm, loop {best['Wl']:.2f}: {best['plate_m']:,.0f} kg/m, tape {best['tape']:.1f} kg/m, rating {best['margin']:.2f} x, sat {best['sat']*100:.0f} %;"
                      f"  with the casing at 8.6: {alt['Wp']:.2f} x {alt['T']*1e3:.0f} mm, loop {alt['Wl']:.2f}: {alt['plate_m']:,.0f} + 8.6 x 2 x {alt['tape']:.1f} = {alt['plate_m'] + CABLE_PER_CASING_KG*2*alt['tape']:,.0f}")
            else:
                print(f"  margin {m:.2f}: no point of the grid meets it")
    return picks


# =============================================================================
# shared helpers
# =============================================================================
FINAL = os.path.join(OUT, "bearing_round5_final.json")
MU0 = 4e-7 * math.pi
TAPE_PITCH = bd.TAPE_T / be.F_CABLE                                   # 0.107 mm per soldered layer
N_M = 4 * N_FAM / L_FAST                                              # fast coils per meter of ring (four tracks): 3.2
W_HALF = ring.M_LOAD_M * 7.918 / 2                                    # the support track's normal share of the casing's weight, 47.5 kN/m [bearing_electrical.W_HALF]


def hold_time(u, L_coil=L_FAST, tau=TAU):
    """One coil's energized time in a pass at speed u and the energized length (final_bearing C.2 with the rise free)."""
    L_on = u * max(0.05, bd.L_SC / u) + bd.L_SC + 2 * L_coil + 2 * u * tau
    return L_on / (u + bd.V_REL), L_on


def launch_integrals(W_req, F_req, tau=TAU, L_coil=L_FAST, n_fam=N_FAM):
    """Per launch along the sled's run: the magnetic energy cycled through the pulsing pair (J), the full-current-equivalent coil-seconds per
    meter of ring, the run length; bearing_round4.launch_integrals with the rise and the requirement free."""
    t_p, v_p, P_p = pz.lsm_profile()
    E_cycle = 0.0; lead_s = 0.0; x_run = 0.0; t_run = float(t_p[-1])
    for i in range(1, len(t_p)):
        u = float(v_p[i]); dt = float(t_p[i] - t_p[i - 1]); dx = u * dt
        if u <= 0:
            continue
        F, kind = r4.pull_track_of_u(u)
        E_cycle += 2 * W_req * (F / F_req) * dx
        lead_s += (2 * n_fam / L_coil) * (F / F_req) * hold_time(u, L_coil, tau)[0] * dx
        x_run += dx
    return dict(E_cycle=E_cycle, lead_s=lead_s, x_run=x_run, t_run=t_run, t_p=t_p, v_p=v_p)


def load_final():
    if os.path.exists(FINAL):
        with open(FINAL, encoding="utf-8") as f:
            return json.load(f)
    return None


def fast_basics(res_fast, N):
    """The fast coil on a given plate result: ampere-turns per coil, current per turn, layers, voltage per turn and at the terminals at TAU."""
    Phi = 2 * res_fast["W_mag"] / res_fast["NI_req"]
    NI_coil = res_fast["NI_req"] / N_FAM; I = NI_coil / N
    tapes = NI_coil / (0.5 * res_fast["ic_req"]); layers = tapes / N
    V_turn = Phi * L_FAST / TAU; V = N * V_turn
    E_coil = res_fast["W_mag"] * L_FAST / N_FAM
    return dict(Phi=Phi, NI_coil=NI_coil, I=I, tapes=tapes, layers=layers, V_turn=V_turn, V=V, E_coil=E_coil, MVA=2 * E_coil / TAU)


# =============================================================================
# S2  the final points, the coupling and the plate heat
# =============================================================================

def two_track(d_edge, Wp_s, T_s, Wl_s, Wp_l, T_l, Wl_l, NI_s=0.0, NI_l=0.0, h=0.01):
    """Stage B's two_track_linear with the plates and loops free: the support track (two loops of Wl_s, plate Wp_s x T_s) at x = 0 and the load
    track (one loop of Wl_l, plate Wp_l x T_l) beside it, plate edge to plate edge d_edge, linear plates, uniform 0.1 m bundles."""
    g = ms.Grid(9.0, 3.0, h, x0=-3.0, y0=-1.0)
    mu = cf.PLATE_MU; B = sb.BUNDLE
    g.add_iron(-Wp_s / 2, Wp_s / 2, bd.GAP0, bd.GAP0 + T_s, mu)
    xl = Wp_s / 2 + d_edge + Wp_l / 2
    g.add_iron(xl - Wp_l / 2, xl + Wp_l / 2, bd.GAP0, bd.GAP0 + T_l, mu)
    legs_s = [(-Wl_s, NI_s), (0.0, -2 * NI_s), (Wl_s, NI_s)]
    legs_l = [(xl - Wl_l / 2, NI_l), (xl + Wl_l / 2, -NI_l)]
    for x, a in legs_s + legs_l:
        if a:
            g.add_current(x - B / 2, x + B / 2, -B, 0.0, a)
    g.solve()
    return g, legs_s, legs_l, xl


def coupling(rs, rl, d_edge=SPACING_MIN):
    """The load track's flux through the support's figure-eight per ampere-turn against the support's own, and what flux conservation does to
    the support force when the load track pulls its 198 kN/m beside a support at its normal current."""
    g_s, legs_s, _, xl = two_track(d_edge, rs["Wp"], rs["T"], rs["Wl"], rl["Wp"], rl["T"], rl["Wl"], NI_s=1e3)
    L_s = sb.fig8_flux(g_s, legs_s)[0] / 1e3
    g_l, legs_s, _, xl = two_track(d_edge, rs["Wp"], rs["T"], rs["Wl"], rl["Wp"], rl["T"], rl["Wl"], NI_l=1e3)
    M = sb.fig8_flux(g_l, legs_s)[0] / 1e3
    NI_l_full = rl["NI_req"] * math.sqrt(fb.F_LOAD_TRACK / rl["F_req"])
    NI_s_op = rs["NI_req"] * math.sqrt(W_HALF / rs["F_req"])
    dNI = -M / L_s * NI_l_full
    return dict(d_edge=d_edge, xl=xl, L_s=L_s, M=M, k=M / L_s, NI_l_full=NI_l_full, NI_s_op=NI_s_op, dNI=dNI, dF=2 * dNI / NI_s_op)


def B_plate_of(res, NI):
    """The plate's peak field at NI from the solved curve: linear below the requirement, log-log between the curve's points above it."""
    NIs = [c[0] for c in res["curve"]]; Bs = [c[3] for c in res["curve"]]
    if NI <= NIs[0]:
        return Bs[0] * NI / NIs[0]
    return math.exp(float(np.interp(math.log(NI), np.log(NIs), np.log(Bs))))


def plate_heat_launch(rl, tau=TAU, L_coil=L_FAST):
    """The hysteresis and eddy heat of a launch in the fast plates (Stage B P4's model with the optimum plate's dimensions and its solved peak
    field at each pass's ampere-turns), per plate per pass and per launch."""
    li = launch_integrals(rl["W_mag"], rl["F_req"], tau, L_coil)
    t_p, v_p = li["t_p"], li["v_p"]
    rows = []; per_u = {}
    for u_km in (1, 3, 5, 10, 20, 30):
        u = u_km * 1e3; F, kind = r4.pull_track_of_u(u)
        NI = rl["NI_req"] * math.sqrt(F / rl["F_req"]); Bp = B_plate_of(rl, NI)
        L_on = hold_time(u, L_coil, tau)[1]; v_sw = u + bd.V_REL
        import levitation_config as lc
        hyst = lc.FERRO_HYSTERESIS_COEFF * (Bp / 1.2) ** 1.6 * rl["T"] * rl["Wp"]
        eddy = bd.eddy_cycle(Bp, bd.GAP0 / v_sw, L_on / v_sw, d_p=rl["T"], mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY) * rl["Wp"]
        rows.append(dict(u=u, kind=kind, F=F, NI=NI, Bp=Bp, hyst=hyst, eddy=eddy, hold=L_on / v_sw))
        per_u[u] = 2 * (hyst + eddy)
    us = sorted(per_u); E = 0.0
    for i in range(1, len(t_p)):
        vv = float(v_p[i]); dtt = float(t_p[i] - t_p[i - 1])
        if vv < us[0]:
            continue
        E += float(np.interp(vv, us, [per_u[u] for u in us])) * (vv + bd.V_REL) * dtt
    return dict(rows=rows, E_launch=E, li={k: v for k, v in li.items() if k not in ("t_p", "v_p")})


def part_final(points, h=0.0025):
    """The chosen points at the design's resolution with the full curve; for the pair labelled 'opt135' the coupling, the plate heat, the hardware."""
    hr(f"S2  THE OPTIMUM PLATES AT {h*1e3:g} mm CELLS")
    cache = load_cache(); results = {}
    for kind, Wp, T, Wl, label in points:
        k = key_of(kind, Wp, T, Wl, h)
        if k in cache and cache[k].get("ok") and len(cache[k].get("ks", [])) >= 3:
            res = cache[k]
        else:
            res = solve_point(kind, Wp, T, Wl, h, ks=(1.25, 1.5, 2.0), verbose=True); res["ok"] = True
            cache = load_cache(); cache[k] = res; save_cache(cache)
        res["label"] = label; results[f"{kind}:{label}"] = res
    print(f"\n  {'track, plate, loop':>40} {'kA':>4} {'worst tape':>22} {'turns/leg':>9} {'tape kg/m':>9} {'sat %':>5} {'T in plate':>10} {'plate kg/m':>10} {'rating kN/m, x':>14} {'cap kN/m':>8} {'at 2x':>6} {'sat at 2x':>9} {'plate + 8.6 x 2 tape':>20}")
    for key, r in results.items():
        print(f"  {key + ': ' + f'{r['Wp']:.2f} x {r['T']*1e3:.0f} mm, loop {r['Wl']:.2f}':>40} {r['NI_req']/1e3:4.0f} {f'{r['B']:.2f} T at {r['th']:.0f} deg, {r['ic_req']:.0f} A':>22} {r['turns_leg']:9,.0f} {r['tape']:9.1f} {r['sat']*100:5.1f} {r['B_plate']:10.2f} {r['plate_m']:10,.0f} "
              f"{f'{r['F70']/1e3:.0f}, {r['margin']:.2f}':>14} {r['F100']/1e3:8.0f} {r['F_last']/1e3:6.0f} {r['curve'][-1][2]*100:9.0f} {r['plate_m'] + CABLE_PER_CASING_KG*2*r['tape']:20,.0f}")
    out = dict(results=results)
    rs = results.get("support:opt135"); rl = results.get("fast:opt135")
    if rs and rl:
        print()
        print(f"  The design pair (margin 1.35): support {rs['Wp']:.2f} m x {rs['T']*1e3:.0f} mm with {rs['Wl']:.2f} m loops; lift and load {rl['Wp']:.2f} m x {rl['T']*1e3:.0f} mm with a {rl['Wl']:.2f} m loop.")
        hw = fb.M_LIM_PLATES + 2 * rs["plate_m"] + 4 * rl["plate_m"]
        tape_load = rl["tape"] * math.sqrt(fb.F_LOAD_TRACK / rl["F_req"])
        tape_casing = 2 * rs["tape"] + 2 * tape_load + 2 * rl["tape"]
        side = ring.CABLE_SIDE; fit = 2 * rl["Wp"] + 2 * rs["Wp"] + 2 * SPACING_MIN
        print(f"  Hardware {hw:,.0f} kg/m (LIM plates {fb.M_LIM_PLATES:,.0f}, two support plates {rs['plate_m']:,.0f}, four fast plates {rl['plate_m']:,.0f}); the book's 10,686.  Tape per meter of casing")
        print(f"  (two support {rs['tape']:.1f}, two load {tape_load:.1f}, two lift {rl['tape']:.1f}): {tape_casing:.0f} kg/m; the bearing's casing side (tape x 2, Part 3b) {2*tape_casing:,.0f} kg/m; tape around the")
        print(f"  ring {tape_casing*ring.L_RING/1e9:.2f} Mt.  The plates on the {side:.2f} m underside: {fit:.2f} m of plate and spacing ({'fits' if fit <= side else 'DOES NOT FIT'}); loop widths {2*rs['Wl']:.2f} m inside a")
        print(f"  {rs['Wp']:.2f} m support plate ({(rs['Wp']-2*rs['Wl'])*1e3/2:.0f} mm of plate past each outer leg's centre) and {rl['Wl']:.2f} m inside a {rl['Wp']:.2f} m fast plate ({(rl['Wp']-rl['Wl'])*1e3/2:.0f} mm).")
        cp = coupling(rs, rl)
        print(f"  Coupling at {cp['d_edge']:.1f} m edge spacing (load plate centre {cp['xl']:.2f} m from the support's): the load loop's flux through the support's figure-eight {cp['M']*1e6:.4f} uWb/m per At")
        print(f"  against the support's own {cp['L_s']*1e6:.3f}: {cp['k']*100:.2f} percent; at the load track's {fb.F_LOAD_TRACK/1e3:.0f} kN/m ({cp['NI_l_full']/1e3:.0f} kA) flux conservation moves the support's {cp['NI_s_op']/1e3:.0f} kA by {cp['dNI']/1e3:.2f} kA,")
        print(f"  {cp['dF']*100:.2f} percent of its force (Stage B S2.1: 0.9 percent at 0.5 m with 1.4 m and 1.2 m plates).")
        ph = plate_heat_launch(rl)
        print(f"  Plate heat of a launch in the fast plates ({rl['Wp']:.2f} m x {rl['T']*1e3:.0f} mm, the solved peak field at each pass's ampere-turns, 25 m coils, the 10 ms rise):")
        print(f"  {'u km/s':>7} {'winding':>7} {'pull/track kN/m':>15} {'kA':>5} {'plate peak T':>12} {'hold ms':>7} {'hyst J/m':>8} {'eddy J/m':>8} {'per m of ring J/m':>17}")
        for r in ph["rows"]:
            print(f"  {r['u']/1e3:7.0f} {r['kind']:>7} {r['F']/1e3:15.0f} {r['NI']/1e3:5.0f} {r['Bp']:12.2f} {r['hold']*1e3:7.1f} {r['hyst']:8.1f} {r['eddy']:8.1f} {2*(r['hyst']+r['eddy']):17.1f}")
        print(f"  Per launch {ph['E_launch']/1e9:.1f} GJ in the plates (Stage B P4: 231 GJ on the book's plates, 185 on six 1.2 m x 150 mm); into the cable, which is warm, not the cold loop.")
        out.update(hw=hw, tape_casing=tape_casing, tape_load=tape_load, coupling=cp, plate_heat=dict(E_launch=ph["E_launch"], rows=ph["rows"]))
    with open(FINAL, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"  -> {FINAL}")
    return out


# =============================================================================
# S3  THE STORES
# =============================================================================

def part_stores():
    hr("S3  THE STORES: THE PULSE POWER OF THE PULLED ZONE, THE THREE WAYS TO FEED IT, CAPACITORS AGAINST SMES, THE RIDE-THROUGH STORE")
    fin = load_final(); r_sup, r_lift = be.load_records()
    if fin and "fast:opt135" in fin["results"]:
        rl = fin["results"]["fast:opt135"]; W_req, F_req, NI_req, src = rl["W_mag"], rl["F_req"], rl["NI_req"], "the optimum plate"
    else:
        W_req, F_req, NI_req, src = r_lift["W_req"], r_lift["F_req"], r_lift["NI_req"], "Part 3b's record (the final plates not yet run)"
    L_ring = ring.L_RING; n_sites = L_ring / 500.0
    print(f"  The fast track's stored energy at the requirement W = {W_req/1e3:.1f} kJ per meter of track ({src}); the energy is magnetic and comes back when the pull ends.")
    # S3.1 the pulse power of the zone
    print()
    print(f"  S3.1  The power at the front and the back of the pulled zone.  The energy per meter of track at a pass is W x F(u) / F_req (E goes as NI^2 and NI^2 as F);")
    print(f"  the front of the zone takes it in at the speed u and the back, the energized length behind, gives it back.  Per track and for the pair:")
    print(f"  {'u km/s':>7} {'winding':>7} {'pull/track kN/m':>15} {'kJ/m':>6} {'front GW/track':>14} {'pair GW':>8} {'zone m':>7} {'zone sites':>10} {'front past a site ms':>20} {'x a site s 12 MW':>16}")
    rows = []
    for uk in (1, 3, 5, 10, 20, 30):
        u = uk * 1e3; F, kind = r4.pull_track_of_u(u); Wu = W_req * F / F_req; P = Wu * u
        hold, L_on = hold_time(u)
        rows.append(dict(u=u, kind=kind, F=F, Wu=Wu, P=P, L_on=L_on))
        print(f"  {uk:7.0f} {kind:>7} {F/1e3:15.0f} {Wu/1e3:6.1f} {P/1e9:14.3f} {2*P/1e9:8.3f} {L_on:7.0f} {L_on/500:10.1f} {500/u*1e3:20.1f} {2*P/12e6:16.0f}")
    r30 = rows[-1]; r1 = rows[0]
    print(f"  Paul's estimate holds at 30 km/s: {r30['P']/1e9:.1f} GW per lift track, {2*r30['P']/1e9:.1f} for the pair, returned {r30['L_on']/1e3:.2f} km behind (the energized length with 25 m coils and the")
    print(f"  10 ms rise; 1.4 km with the 2 ms rise and no coil length), {500/30e3*1e3:.0f} ms past a 500 m site, {2*r30['P']/12e6:.0f} times its 12 MW.  At 1 km/s the load pair takes {2*r1['P']/1e6:.0f} MW at its front,")
    print(f"  {2*r1['P']/12e6:.0f} times a site (Part 3b's launch pull of 285 kN/m on the pair; the dynamics' spread pull at 1 km/s is the full 396 kN/m, {2*W_req*1e3/1e9:.2f} GW at the front).  Between the front and the back nothing is consumed but the losses (the converter's, the tape's, the plates'): the net draw of a launch")
    li = launch_integrals(W_req, F_req)
    print(f"  from the grid is the losses on the {li['E_cycle']/1e12:.1f} TJ cycled per launch, {li['E_cycle']/n_sites/1e6:.0f} MJ per site per launch, {li['E_cycle']*LAUNCHES_PER_DAY/86400/n_sites/1e3:.0f} kW per site averaged at 100 percent loss,")
    print(f"  {0.03*li['E_cycle']*LAUNCHES_PER_DAY/86400/n_sites/1e3:.0f} kW at 3 percent, {0.03*li['E_cycle']*LAUNCHES_PER_DAY/86400/n_sites/12e6*100:.2f} percent of the site's 12 MW.")
    # S3.2 three ways
    print()
    E_m_ring = 4 * W_req                                                  # four tracks, every coil its own full store
    E_m_shared = 2 * W_req                                                # one pair pulses at a time
    cap_kJ_kg = 1.5                                                       # Kerrigan CDE 2007, read in full in round 4
    I_bus = 2 * r30["P"] / 100e3
    n_bus_tapes = I_bus / 650.0
    m_bus = 2 * n_bus_tapes * 0.012 * bd.TAPE_T * bd.TAPE_RHO
    q_norris = float(lh.q_norris_strip(650.0, 892.0))                     # J per m per tape per cycle at the bus's 650 A of 892 A
    E_bus_launch = q_norris * 2 * n_bus_tapes * li["x_run"]               # every meter of bus sees one pulse per lap, both directions
    print(f"  S3.2  Three ways to feed the coils (the switches and the converter at every coil in all three, Section 4; the mass is per meter of casing):")
    print(f"  (a) a full local store at every coil, round 4: {E_m_ring/1e3:.0f} kJ per meter of ring ({W_req*L_FAST/N_FAM/1e6:.3f} MJ per coil), film capacitors at {cap_kJ_kg:.1f} kJ/kg: {E_m_ring/cap_kJ_kg/1e3:.0f} kg/m; warm, nothing into")
    print(f"      the cold but the leads (S6); the site's connection carries the losses only ({0.03*li['E_cycle']*LAUNCHES_PER_DAY/86400/n_sites/1e3:.0f} kW at 3 percent).")
    print(f"  (b) no store: the pulses drawn from the 100 kV section bus and returned to it.  The bus must carry the circulating {2*r30['P']/1e9:.1f} GW between the back and the front of")
    print(f"      the zone, {r30['L_on']/1e3:.2f} km apart inside one 10 km section: {I_bus/1e3:.0f} kA at 100 kV against the bus's 11,600 A at its 1,160 MW maximum (Chapter 8), {I_bus/11600:.1f} times it;")
    print(f"      in tape, {n_bus_tapes:.0f} of the 12 mm tapes at 650 A in each direction against 18, {m_bus:.1f} kg/m of tape; the bus's own pulse loss (Norris, 650 of 892 A) {q_norris*1e3:.1f} mJ per")
    print(f"      meter per tape per pulse, {E_bus_launch/1e9:.2f} GJ per launch, {E_bus_launch*LAUNCHES_PER_DAY/86400/L_ring:.2f} W per meter of ring at 500 launches a day, in the cold.  The converter at each coil then")
    print(f"      steps 100 kV down to the coil's 14 or 20 kV and must block the bus voltage (Section 4); the tap point from the HVDC line still sees only the losses.")
    print(f"      Mass: the bus's extra tape and the heavier converters; store none.  What is new in the cold: {E_bus_launch*LAUNCHES_PER_DAY/86400/L_ring:.2f} W/m against the standing {Q_STANDING_W_PER_M:.0f}.")
    print(f"  (c) a store shared along the casing: one pair pulses at a time (the load pair below orbital speed, the lift pair above), so a store shared by the")
    print(f"      four tracks holds {E_m_shared/1e3:.0f} kJ per meter, half of (a), {E_m_shared/cap_kJ_kg/1e3:.0f} kg/m in film capacitors; per 500 m site {E_m_shared*500/1e6:.0f} MJ.  Sharing further along the casing buys")
    print(f"      nothing: the zone is {r30['L_on']/500:.1f} sites long at 30 km/s and every meter of it is energized at once, so a store sized for the zone is a store for every meter")
    print(f"      of the zone.  Carrying the back's return to the front instead of storing it is (b): at the coil's 14 kV the circulating current is {2*r30['P']/14.3e3/1e3:.0f} kA, at 100 kV {I_bus/1e3:.0f}.")
    print(f"      Heat into the cold: the leads only, as (a).  The site's connection: the losses, as (a).")
    # S3.3 capacitors against SMES
    print()
    smes_kJ_kg = 12.5; smes_record = 11.0
    print(f"  S3.3  The store itself.  Film capacitors: {cap_kJ_kg:.1f} kJ/kg and 1.2 to 1.5 J/cc at 75 percent of the film's breakdown stress (Kerrigan, CDE, CARTS 2007, read in full");
    print(f"  in round 4), warm, no cold load, life set by the pulse count (500 a day is 1.8 million a year: the number Kerrigan's low-duty units are made for is hundreds to")
    print(f"  thousands, so the capacitor must be derated or run at a lower stress [UNVERIFIED for this pulse count]).  SMES: the energy per kilogram is set by the virial theorem,")
    print(f"  W = M_min sigma / d (Tixador 2008, read in full): {smes_kJ_kg:.1f} kJ/kg for a steel structure at 100 MPa, the CMS magnet {smes_record:.0f} kJ/kg, 'the order of ten kJ/kg'; the")
    print(f"  conductor, cryostat and converter come on top of the structure.  Per meter of ring: (a) {E_m_ring/smes_kJ_kg/1e3:.0f} kg/m of structure, (c) {E_m_shared/smes_kJ_kg/1e3:.0f}, against {E_m_ring/cap_kJ_kg/1e3:.0f} and")
    print(f"  {E_m_shared/cap_kJ_kg/1e3:.0f} in capacitors: {cap_kJ_kg/smes_kJ_kg*100:.0f} percent of the capacitors' mass.  But an SMES sits in the cold, and its AC loss on every charge and discharge, 'some low %' of the")
    print(f"  energy cycled (Tixador), lands in the 70 K loop: round 4's cold-converter rule, {0.01*li['E_cycle']*LAUNCHES_PER_DAY/86400/L_ring:.0f} W per meter of ring per percent lost, {0.02*li['E_cycle']*LAUNCHES_PER_DAY/86400/L_ring:.0f} W/m at 2 percent against")
    print(f"  the standing {Q_STANDING_W_PER_M:.0f}; an SMES store is excluded on the cold load unless its loss per cycle is well under a percent.  It is also a current source feeding an")
    print(f"  inductive coil: the converter between them is the same H-bridge as with a capacitor, so it saves nothing there.  [SPECULATIVE] A structure of the ring's")
    print(f"  own CNT at 12.5 GPa and 1,500 kg/m3 would put the virial bound at {12.5e9/1500/1e6:.1f} MJ/kg; the conductor and cryostat would then set the mass, not the structure.")
    print(f"  Chapter 8's SMES is for the grid's transients, the surge of a launch and the sled's passage at the section level, where its cold loop is the section's own;")
    print(f"  the coil-by-coil store of (a) and (c) is a different job and the capacitor does it warm.")
    # S3.4 the ride-through
    print()
    k_track = 0.33e3 * 1e3                                                # N/m per m of gap per track, the constant-flux negative stiffness
    for x_mm in (1.0, 5.0):
        F_corr = 2 * k_track * x_mm * 1e-3                                # both support tracks' pull-in at x, held by the fast pair
        frac = F_corr / (2 * F_req); E_corr = frac * (W_req * L_FAST / N_FAM)
        print(f"  {'S3.4  The ride-through store.  ' if x_mm <= 1 else '        '}A {x_mm:.0f} mm offset with the support frozen pulls {F_corr/1e3:.2f} kN/m (2 x 0.33 kN/m per mm); the fast pair holds it with {frac*100:.2f} percent of its full pull,")
        print(f"        {math.sqrt(frac)*100:.1f} percent of the full current, {E_corr:.0f} J borrowed per coil and returned.")
    print(f"  The energy of the control is not the issue; what drains a local store with the bus gone is the converter's own standing draw (gate drives, control, the")
    print(f"  current sensors) and the loss on each correction.  Neither is in a source read; the table sizes the store for a stated draw per coil and a stated time:")
    print(f"  {'draw W/coil':>11} | {'10 min kJ':>9} {'kg':>5} | {'1 h kJ':>7} {'kg':>5} | {'24 h kJ':>8} {'kg':>6} | {'the pulse store (a) lasts':>24}")
    E_pulse = W_req * L_FAST / N_FAM
    for P in (10.0, 30.0, 100.0):
        print(f"  {P:11.0f} | {P*600/1e3:9.1f} {P*600/cap_kJ_kg/1e3:5.1f} | {P*3600/1e3:7.0f} {P*3600/cap_kJ_kg/1e3:5.1f} | {P*86400/1e3:8.0f} {P*86400/cap_kJ_kg/1e3:6.0f} | {E_pulse/P/60:24.0f} min")
    print(f"  (kg at {cap_kJ_kg:.1f} kJ/kg per coil; {N_M:.1f} coils per meter of ring.)  The local pulse store of (a), {E_pulse/1e6:.3f} MJ, is also the ride-through store: at 30 W per coil it runs the gap")
    print(f"  control for {E_pulse/30/3600:.1f} hours with no bus, at 100 W {E_pulse/100/60:.0f} minutes.  With (c) the shared store does the same for its site.  With (b) there is no store and the")
    print(f"  ride-through must be added: 1 hour at 30 W is {30*3600/1e3:.0f} kJ per coil, {30*3600/cap_kJ_kg/1e3*N_M:.0f} kg per meter of ring in capacitors.  The stated time to size for is Paul's: the")
    print(f"  support coils hold the casing's weight for years (round 4, 2.6) and the gap needs only the fast loop's small corrections, so the store is sized for the")
    print(f"  time to restore the bus or to bring the casing onto its stops, not for the launches.")
    return dict(rows=rows, E_cycle=li["E_cycle"], E_m_ring=E_m_ring, E_m_shared=E_m_shared, I_bus=I_bus, m_bus=m_bus, E_bus_launch=E_bus_launch, E_pulse=E_pulse)


# =============================================================================
# S4  THE SWITCHES
# =============================================================================

# the devices read in full this round (datasheets): blocking V, DC or average current A, peak current A (and its time), mass kg
DEVICES = {
    "StakPak IGBT 5SNA 3000K452300 (ABB 5SYA 1450-03, 2017)": dict(V=4500.0, V_work=2800.0, I_dc=3000.0, I_pk=6000.0, t_pk=1e-3, I_diode_10ms=21000.0, m=4.3, turn_off=True,
                                                                  note="press-pack, 247 x 237 x 31.5 mm, diode inside (2:1); V_work the datasheet's test VCC"),
    "phase control thyristor 5STP 42U6500 (ABB 5SYA1043-03, 2007)": dict(V=6500.0, V_work=3300.0, I_dc=3460.0, I_pk=71.4e3, t_pk=10e-3, m=3.6, turn_off=False,
                                                                        note="VDWM 3300 V crest working; I_TSM 71.4 kA for a 10 ms half sine; cannot turn off against current"),
}
DEVICES_UNREAD = {
    "IGCT 5SHY 42L6500 [search summary, not read in full]": dict(V=6500.0, V_work=4000.0, I_dc=1290.0, I_pk=3800.0, t_pk=None, m=2.9, turn_off=True, note="VDC 4000 V, turn-off 3,800 A, 2.9 kg from a distributor's copy of the datasheet"),
}


def part_switches():
    hr("S4  THE SWITCHES: WHAT THEY DO, AND THEIR MASS PER METER AT 10 MS FOR 10 AND 14 TURNS")
    fin = load_final(); r_sup, r_lift = be.load_records()
    rl = fin["results"]["fast:opt135"] if fin and "fast:opt135" in fin["results"] else dict(W_mag=r_lift["W_req"], NI_req=r_lift["NI_req"], ic_req=r_lift["ic3"], F_req=r_lift["F_req"])
    print(f"  What the switch has to do, per coil: charge the coil from the store in {TAU*1e3:.0f} ms (+V across the coil), hold the current for the pass (zero voltage: the current")
    print(f"  free-wheels through one switch and one diode), and discharge it in {TAU*1e3:.0f} ms (-V: the energy back into the store); then carry the regulation's corrections.")
    print(f"  That is a two-quadrant converter (one current sign, both voltage signs): an asymmetric bridge of two switches and two diodes, or a full H-bridge.  A")
    print(f"  thyristor cannot do it: it turns off only when the current passes zero, so a thyristor and a capacitor make a half-sine ring, not a {TAU*1e3:.0f} ms ramp, a hold and")
    print(f"  a ramp down.  The 71 kA surge rating of a 6.5 kV thyristor (read in full) is of no use here; the switch is an IGBT or an IGCT.")
    print(f"  {'turns':>5} {'A per turn':>10} {'V terminals':>11} {'E coil MJ':>9} {'MVA':>5} {'hold at 30 km/s ms':>18} | device: series count x parallel x positions, kg per coil, kg per m of casing")
    out = {}
    for N in TURNS:
        fbk = fast_basics(rl, N)
        hold = hold_time(30e3)[0]
        line = f"  {N:5d} {fbk['I']:10,.0f} {fbk['V']/1e3:9.1f} kV {fbk['E_coil']/1e6:9.3f} {fbk['MVA']/1e6:5.0f} {hold*1e3:18.0f} |"
        print(line)
        out[N] = {}
        for name, d in list(DEVICES.items()) + list(DEVICES_UNREAD.items()):
            if not d["turn_off"]:
                continue
            I_rating = d["I_pk"] if d["t_pk"] is None else d["I_dc"]               # the IGCT by its turn-off current (its 1,290 A is a thermal average), the IGBT by its DC rating
            n_s = math.ceil(fbk["V"] / d["V_work"]); n_p = math.ceil(fbk["I"] / I_rating)
            for positions, lab in ((2, "asymmetric bridge, 2 switches + 2 diodes"), (4, "H-bridge, 4 switches")):
                m_coil = positions * n_s * n_p * d["m"] * (1.0 if d is not DEVICES.get("phase control thyristor 5STP 42U6500 (ABB 5SYA1043-03, 2007)") else 1.0)
                m_diodes = (2 * n_s * n_p * d["m"]) if positions == 2 else 0.0        # the diodes of the asymmetric bridge as press-packs of the same mass [assumption]
                m_tot = m_coil + m_diodes
                print(f"  {'':>62}  {name}: {n_s} x {n_p} x {positions} ({lab}): {m_tot:,.0f} kg per coil, {m_tot*N_M:,.0f} kg/m")
                out[N][f"{name}|{positions}"] = dict(n_s=n_s, n_p=n_p, m_coil=m_tot, m_per_m=m_tot * N_M)
    print(f"  The series count takes the device at its datasheet working voltage ({DEVICES['StakPak IGBT 5SNA 3000K452300 (ABB 5SYA 1450-03, 2017)']['V_work']:.0f} V of 4,500 for the StakPak, the test condition of its switching figures; the IGCT")
    print(f"  at its 4,000 V DC rating), one in parallel since the turn current is under the DC rating (3,000 A) and far under the 1 ms peak (6,000 A); the pulse is a")
    print(f"  {hold_time(30e3)[0]*1e3:.0f} ms hold at 30 km/s at a duty of 0.05 percent, so the thermal rating is not reached and the datasheet's DC current is the conservative choice.")
    print(f"  Clamps, gate units, bus bars and the cooling of a press-pack stack are not in the device mass [UNVERIFIED: a stack weighs of the order of the devices again].")
    print(f"  Round 4's 300 to 500 kg/m, from memory, lands inside this range; the figure now rests on the StakPak datasheet read in full.  With (b) of Section 3 the converter")
    print(f"  must also block the 100 kV bus: {math.ceil(100e3/2800)} StakPaks in series per position, {math.ceil(100e3/2800)*2*4.3*N_M:,.0f} kg/m for the two switch positions alone, four times the 14 kV bridge.")
    return out


# =============================================================================
# S5  THE TAPE'S PULSE LOSS
# =============================================================================

def q_bean_slab(Ba, Bp):
    """Bean critical-state slab in a parallel field cycling with amplitude Ba (T): loss per cycle per m^3 [textbook result, e.g. Wilson's
    Superconducting Magnets; not read in full this round]."""
    if Ba <= Bp:
        return (2.0 / 3.0) * Ba**3 / (MU0 * Bp)
    return (2.0 / MU0) * (Ba * Bp - 2.0 * Bp**2 / 3.0)


def part_tapeloss(h=0.005):
    hr(f"S5  THE PULSE LOSS OF THE SOLDERED STACKS ON THE OPTIMUM PLATES AT {TAU*1e3:.0f} ms ({h*1e3:g} mm cells)")
    fin = load_final(); r_sup, r_lift = be.load_records(); theva, theva77 = be.load_theva()
    if fin and "fast:opt135" in fin["results"]:
        rl = fin["results"]["fast:opt135"]
    else:
        print("  the final plates are not run yet; using Part 3b's 1.4 m x 200 mm plate and 0.8 m loop")
        rl = dict(Wp=1.4, T=0.2, Wl=0.8, NI_req=r_lift["NI_req"], F_req=r_lift["F_req"], ic_req=r_lift["ic3"], W_mag=r_lift["W_req"], turns_leg=r_lift["turns3"] / 2)
    geom = dict(bo.BEST_LIFT); geom["W"] = rl["Wl"]; w = 0.003
    with r4.plate_patch(rl["Wp"], rl["T"]):
        tr = be.Track(geom, h, w_tape=w, I_ref=bo.I_REF_3, NI_ref=rl["NI_req"])
        g, s = tr.solve(tr.NI_design(rl["NI_req"]))
        fields = [bo.field_at_points(g, xs, ys) for xs, ys in tr.lattices]
    print(f"  The lift track at its requirement on the {rl['Wp']:.2f} m x {rl['T']*1e3:.0f} mm plate: {rl['NI_req']/1e3:.0f} kA per leg (solved {s['F']/1e3:.1f} kN/m), {len(tr.legs)} legs of {len(tr.lattices[0][0])} model tapes each for")
    print(f"  {rl['turns_leg']:,.0f} real tapes per leg.  Each tape's field and angle from the solve; its critical current from the Robinson 70 K angle data (lim_hysteresis.ic_tape).")
    print(f"  Per pulse (0 to the pass's field and back, a cycle of half that amplitude): the perpendicular component through each layer's face on the Brandt-Indenbom")
    print(f"  strip model (per layer, the layers soldered or not); the component along the face on two models: (i) the layers insulated, the film's own slab loss, which")
    print(f"  Chapter 7 and lim_hysteresis drop as under a thousandth of the strip loss; (ii) the layers soldered as one slab of the stack's thickness with the stack's")
    print(f"  engineering current density (Bean slab), which is what point 4 of the layout brief asked.  Rate: both are hysteretic and do not depend on the rise;")
    print(f"  the solder's coupling time constant mu0 d^2 / (pi^2 rho) is of the order of 0.1 ms for a 7.5 mm stack at 1e-7 ohm m [estimate], far under the 10 ms rise,")
    print(f"  so the stack is fully coupled and (ii) applies to soldered layers.")
    tapes_per_leg = rl["turns_leg"]
    scale = tapes_per_leg / len(tr.lattices[0][0])
    out = {}
    for N, layers_lab in ((10, 70), (14, 50)):
        fbk = fast_basics(rl, N); layers = fbk["layers"]; d = layers * TAPE_PITCH
        res = {}
        for face in ("h", "v"):                                               # h: faces horizontal (normal along y); v: faces vertical (normal along x)
            q_perp = 0.0; q_par_film = 0.0; q_slab = 0.0; n = 0
            Bperp_list = []; Bpar_list = []
            for (Bx, By) in fields:
                Bm = np.hypot(Bx, By)
                Bn = np.abs(By) if face == "h" else np.abs(Bx)                 # through the face
                Bt = np.abs(Bx) if face == "h" else np.abs(By)                 # along the face
                th = np.degrees(np.arctan2(Bt, Bn))                           # from the normal
                Ic = lh.ic_tape(Bm, th, w)
                qp = lh.q_strip_perp(Ic, w, Bn / 2.0)                          # J/m per layer per pulse
                Jc_eff = Ic / (w * TAPE_PITCH); Bp = MU0 * Jc_eff * d / 2.0
                vol = w * TAPE_PITCH                                           # m^3 per m per layer
                qs = np.array([q_bean_slab(b / 2.0, bp) for b, bp in zip(Bt, Bp)]) * vol
                Jc_film = Ic / (w * 3.5e-6); Bp_f = MU0 * Jc_film * 3.5e-6 / 2.0
                qf = np.array([q_bean_slab(b / 2.0, bp) for b, bp in zip(Bt, Bp_f)]) * (w * 3.5e-6)
                q_perp += float(qp.sum()) * scale; q_slab += float(qs.sum()) * scale; q_par_film += float(qf.sum()) * scale; n += len(Bm)
                Bperp_list.append(Bn); Bpar_list.append(Bt)
            Bperp = np.concatenate(Bperp_list); Bpar = np.concatenate(Bpar_list)
            res[face] = dict(q_perp=q_perp, q_slab=q_slab, q_film=q_par_film, Bperp_mean=float(Bperp.mean()), Bpar_mean=float(Bpar.mean()), Bp=float(np.mean(Bp)))
        out[N] = res
        print(f"\n  {N} series turns of {layers:.0f} layers, the stack {d*1e3:.1f} mm thick, its penetration field about {res['h']['Bp']:.2f} T (at the mean tape limit):")
        print(f"  {'faces':>10} {'B through face, mean T':>22} {'B along face, mean T':>20} {'perp (Brandt) kJ/m':>19} {'along, film kJ/m':>16} {'along, soldered slab kJ/m':>25} {'soldered total kJ/m':>19}")
        for face, lab in (("h", "horizontal"), ("v", "vertical")):
            r = res[face]
            print(f"  {lab:>10} {r['Bperp_mean']:22.2f} {r['Bpar_mean']:20.2f} {r['q_perp']/1e3:19.2f} {r['q_film']/1e3:16.4f} {r['q_slab']/1e3:25.2f} {(r['q_perp']+r['q_slab'])/1e3:19.2f}")
    # per launch and per meter of ring: scale the fields with sqrt(F/F_req) along the run; the loss goes about as B^2 to B^3 below penetration, so recompute
    print()
    print(f"  Per launch along the run (the lift pair above orbital speed, the load pair below, both of this design; the pass's field scaled as sqrt(F / F_req) and the")
    print(f"  loss recomputed at each speed), and per meter of ring at {LAUNCHES_PER_DAY} launches a day against the standing {Q_STANDING_W_PER_M:.0f} W/m (Chapter 14):")
    li = launch_integrals(rl["W_mag"], rl["F_req"])
    t_p, v_p = li["t_p"], li["v_p"]
    speeds = [1e3, 2e3, 3e3, 5e3, 7e3, 10e3, 15e3, 20e3, 25e3, 30e3]
    launch = {}
    for N, layers_lab in ((10, 70), (14, 50)):
        fbk = fast_basics(rl, N); d = fbk["layers"] * TAPE_PITCH
        per_u = {face: {} for face in ("h", "v")}
        for u in speeds:
            F, kind = r4.pull_track_of_u(u); k = math.sqrt(F / rl["F_req"]) if F > 0 else 0.0
            for face in ("h", "v"):
                q_perp = 0.0; q_slab = 0.0
                if k > 0:
                    for (Bx, By) in fields:
                        Bm = np.hypot(Bx, By) * k
                        Bn = (np.abs(By) if face == "h" else np.abs(Bx)) * k; Bt = (np.abs(Bx) if face == "h" else np.abs(By)) * k
                        th = np.degrees(np.arctan2(Bt, Bn)); Ic = lh.ic_tape(Bm, th, w)
                        q_perp += float(lh.q_strip_perp(Ic, w, Bn / 2.0).sum()) * scale
                        Bp = MU0 * (Ic / (w * TAPE_PITCH)) * d / 2.0
                        q_slab += float(np.array([q_bean_slab(b / 2.0, bp) for b, bp in zip(Bt, Bp)]).sum()) * w * TAPE_PITCH * scale
                per_u[face][u] = (2 * q_perp, 2 * q_slab)                       # the pair
        E = {face: [0.0, 0.0] for face in ("h", "v")}
        for i in range(1, len(t_p)):
            vv = float(v_p[i]); dtt = float(t_p[i] - t_p[i - 1])
            if vv < speeds[0]:
                continue
            for face in ("h", "v"):
                us = sorted(per_u[face])
                E[face][0] += float(np.interp(vv, us, [per_u[face][u][0] for u in us])) * vv * dtt
                E[face][1] += float(np.interp(vv, us, [per_u[face][u][1] for u in us])) * vv * dtt
        launch[N] = E
        for face, lab in (("h", "horizontal"), ("v", "vertical")):
            Ep, Es = E[face]
            print(f"  {N} turns, faces {lab:>10}: perpendicular {Ep/1e9:6.1f} GJ per launch ({Ep*LAUNCHES_PER_DAY/86400/ring.L_RING:5.1f} W/m), soldered slab {Es/1e9:6.1f} GJ ({Es*LAUNCHES_PER_DAY/86400/ring.L_RING:5.1f} W/m), "
                  f"total {(Ep+Es)/1e9:6.1f} GJ, {(Ep+Es)*LAUNCHES_PER_DAY/86400/ring.L_RING:5.1f} W/m = {(Ep+Es)*LAUNCHES_PER_DAY/86400/ring.L_RING/Q_STANDING_W_PER_M:.2f} x the standing load")
    print(f"  (Chapter 7's bound of 605 GJ is all of the field through the face at the strip asymptote w Ic B; 53 GJ is its 5 degree case.  The slab term is new: with the")
    print(f"  layers soldered the field along the face is no longer free, and the bundle's winding has to be chosen on the sum, not on the perpendicular term alone.)")
    out["launch"] = launch
    with open(os.path.join(OUT, "bearing_round5_tapeloss.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)
    return out


# =============================================================================
# S6  THE COLD LOAD
# =============================================================================

def lead_rows(rl):
    """The fast coils' leads for 10 and 14 turns on the optimum plate: pulse-sized, resized for a 20 K mean rise at 500 launches a day (round 4's rule
    with the 10 ms rise's hold times); idle and pulse heat per meter of ring and per site."""
    li = launch_integrals(rl["W_mag"], rl["F_req"])
    duty = li["lead_s"] * LAUNCHES_PER_DAY / 86400 / (N_M * ring.L_RING)
    rows = {}
    for N in TURNS:
        fbk = fast_basics(rl, N); I = fbk["I"]
        A_p = I / r4.J_PULSE_LEAD; R_p = be.RHO_CU_77 * 1.0 / A_p
        P_mean = I**2 * R_p * duty; dT = P_mean * 1.0 / (2 * 400.0 * A_p)
        A_m = A_p * math.sqrt(dT / 20.0) if dT > 20.0 else A_p; R_m = be.RHO_CU_77 * 1.0 / A_m
        q_idle = be.lead_heat(0.0, A_over_l=A_m)[0]
        idle_m = 2 * q_idle * N_M; pulse_m = 2 * I**2 * R_m * duty * N_M
        rows[N] = dict(I=I, A_p=A_p, A_m=A_m, dT=dT, q_idle_lead=q_idle, idle_m=idle_m, pulse_m=pulse_m, idle_site=idle_m * 500, pulse_site=pulse_m * 500, duty=duty)
    return rows, li


def part_coldload():
    hr("S6  THE COLD LOAD: THE LEADS IDLE (THE DEPLOYMENT) AND IN OPERATION, THE PLATES, THE TAPE")
    fin = load_final(); r_sup, r_lift = be.load_records()
    rl = fin["results"]["fast:opt135"] if fin and "fast:opt135" in fin["results"] else dict(W_mag=r_lift["W_req"], NI_req=r_lift["NI_req"], ic_req=r_lift["ic3"], F_req=r_lift["F_req"], Wp=1.4, T=0.2)
    rows, li = lead_rows(rl)
    tl = None
    p = os.path.join(OUT, "bearing_round5_tapeloss.json")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            tl = json.load(f)
    E_pl = fin["plate_heat"]["E_launch"] if fin and "plate_heat" in fin else None
    print(f"  The fast coils' leads: a copper pair per coil, 1 m long, sized for the pulse (Part 1c's 20 K adiabatic rise in a 0.1 s pulse) and then thickened until its mean")
    print(f"  Joule heating at {LAUNCHES_PER_DAY} launches a day lifts its middle by 20 K at most (round 4, 4.2); the duty from the launch profile with the {TAU*1e3:.0f} ms rise: {rows[10]['duty']*1e3:.2f} ms of")
    print(f"  full-current-equivalent per second per coil.  {N_M:.1f} coils per meter of ring, 1,600 per 500 m site.")
    print(f"  {'turns':>5} {'A per turn':>10} {'lead mm2':>8} {'idle W per lead':>15} {'idle W/m of ring':>16} {'idle kW per site':>16} {'pulses W/m':>10} {'pulses kW/site':>14} {'idle + pulses, x standing':>24}")
    for N, r in rows.items():
        print(f"  {N:5d} {r['I']:10,.0f} {r['A_m']*1e6:8.0f} {r['q_idle_lead']:15.2f} {r['idle_m']:16.1f} {r['idle_site']/1e3:16.2f} {r['pulse_m']:10.1f} {r['pulse_site']/1e3:14.2f} {(r['idle_m']+r['pulse_m'])/Q_STANDING_W_PER_M:24.2f}")
    print(f"  Idle (the deployment, no launches): the standing cold load of a site goes from {Q_STANDING_W_PER_M*500/1e3:.1f} kW (Chapter 14, after the deployment) to {Q_STANDING_W_PER_M*500/1e3 + rows[10]['idle_site']/1e3:.1f} kW at 10 turns and")
    print(f"  {Q_STANDING_W_PER_M*500/1e3 + rows[14]['idle_site']/1e3:.1f} kW at 14; during the deployment it is added to the deployment's own cold load (26.5 to 30.8 kW per site, Table 14.2): tools/round5_deployment.py runs it.")
    print(f"  The leads could be withdrawn for the deployment as the support coils' are (round 4), which takes the idle heat to zero at the price of a contact per coil.")
    print()
    print(f"  In operation at {LAUNCHES_PER_DAY} launches a day, per meter of ring and per site (W/m, kW per site):")
    print(f"  {'item':>44} {'10 turns W/m':>12} {'kW/site':>8} {'14 turns W/m':>12} {'kW/site':>8} {'where it lands':>30}")
    def row(lab, a, b, where):
        print(f"  {lab:>44} {a:12.1f} {a*500/1e3:8.2f} {b:12.1f} {b*500/1e3:8.2f} {where:>30}")
    row("the standing load (Chapter 14)", Q_STANDING_W_PER_M, Q_STANDING_W_PER_M, "the 70 K loop")
    row("the leads, idle conduction", rows[10]["idle_m"], rows[14]["idle_m"], "the 70 K loop")
    row("the leads, the pulses' Joule heat", rows[10]["pulse_m"], rows[14]["pulse_m"], "the 70 K loop")
    tot10 = Q_STANDING_W_PER_M + rows[10]["idle_m"] + rows[10]["pulse_m"]; tot14 = Q_STANDING_W_PER_M + rows[14]["idle_m"] + rows[14]["pulse_m"]
    if tl:
        for face, lab in (("h", "horizontal"), ("v", "vertical")):
            Ep10, Es10 = tl["launch"]["10"][face]; Ep14, Es14 = tl["launch"]["14"][face]
            f = LAUNCHES_PER_DAY / 86400 / ring.L_RING
            row(f"the tape, faces {lab}: perpendicular (Brandt)", Ep10 * f, Ep14 * f, "the 70 K loop")
            row(f"the tape, faces {lab}: soldered slab", Es10 * f, Es14 * f, "the 70 K loop")
        best = min(("h", "v"), key=lambda fc: sum(tl["launch"]["10"][fc]))
        Et10 = sum(tl["launch"]["10"][best]) * LAUNCHES_PER_DAY / 86400 / ring.L_RING; Et14 = sum(tl["launch"]["14"][best]) * LAUNCHES_PER_DAY / 86400 / ring.L_RING
        tot10 += Et10; tot14 += Et14
        print(f"  (the total below takes the better orientation, faces {'horizontal' if best == 'h' else 'vertical'})")
    if E_pl is not None:
        f = LAUNCHES_PER_DAY / 86400 / ring.L_RING
        row("the fast plates, hysteresis and eddy", E_pl * f, E_pl * f, "the cable (warm), not the loop")
    row("cold load in operation, total", tot10, tot14, "the 70 K loop")
    print(f"  The plates are on the cable: their heat is radiated by the cable at its own temperature and reaches the coils only through the cable's radiation onto")
    print(f"  them (5 W/m2 at the reference, Chapter 14); it is not pumped at 70 K.  Compressor power: the cold load times 94 (Chapter 14's COP of 0.0106).")
    out = dict(leads={str(k): v for k, v in rows.items()}, tot10=tot10, tot14=tot14)
    with open(os.path.join(OUT, "bearing_round5_coldload.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=float)
    return out


# =============================================================================
# S7  THE CASES AT 10 MS
# =============================================================================

def part_cases():
    hr("S7  THE CASES AT 10 MS ON THE BOOK'S CASING WITH THE OPTIMUM PLATES' CAPS")
    fin = load_final()
    print(f"  The dynamics cap a lift track at its requirement ({fb.F_LIFT_TRACK/1e3:.1f} kN/m), a load track at its ({fb.F_LOAD_TRACK/1e3:.0f}) and a support track alone at {W_HALF*2*245.5/245.5/1e3*0+245.5:.1f} kN/m, the")
    print(f"  worst tape at half its limit (the caps line of every case log).  Those caps are requirements, set before the plate is chosen, and the optimum plates are")
    print(f"  solved to reach them with the worst tape at half its limit: the caps the cases use do not change with the plate, the plate's rating (1.35 x) is headroom")
    print(f"  the cases never draw on, and round 4's 10 ms run on the 12,000 kg/m casing stands.  Its results (runs/bearing_round4_casing_10ms_20fam.txt):")
    p = os.path.join(OUT, "bearing_round4_casing_10ms_20fam.txt")
    if not os.path.exists(p):
        print(f"  (log {p} not found)"); return {}
    lines = open(p, encoding="utf-8", errors="replace").read().splitlines()
    def block(start, n):
        i = next(k for k, l in enumerate(lines) if start in l)
        return lines[i:i + n]
    for l in block("F.1  Sled passage", 8): print("  " + l)
    for l in block("Tolerances (pair)", 11): print("  " + l)
    for l in block("F.2  Single failures", 23): print("  " + l)
    for l in block("At rest (no sled)", 5): print("  " + l)
    for l in block("Beam model", 10): print("  " + l)
    for l in block("F.3  Roll", 9): print("  " + l)
    for l in block("F.4", 4): print("  " + l)
    return dict(log=p)


# =============================================================================
# S8  THE SUPPORT COIL'S CHARGING ORDER
# =============================================================================

def part_charging():
    hr("S8  THE SUPPORT COIL'S TWO CLOSING ORDERS (55 turns of 10 layers, charged at about 110 V)")
    r_sup, r_lift = be.load_records(); c = bl.part1_constants(r_sup); sc = bl.switched_circuit(c, L_SUP)
    v = [x for x in r4.turn_variants(c, sc) if x["N"] == 55][0]
    I0, Lc, Ll = v["I0"], v["L_common"], v["L_leak"]
    V_ch = Lc * 0.01 * I0                                                 # 1 percent per second, ten together
    R_s = V_ch / (0.01 * I0)                                              # the strip takes 1 percent of the coil current while charging (Stage B's rule)
    V_d = 1.0
    t_close = (1.0, 10.0, 60.0)
    print(f"  The coil: {I0:.0f} A normal, {v['I_cap']:.0f} A at the cap; {Lc:.1f} H with the ten coils of the bundle charged together (the common mode), {Ll:.3f} H for one coil trimmed")
    print(f"  alone against nine persistent neighbours (its leakage).  Charging at 1 percent per second: {V_ch:.0f} V across the coil (ten together).  The strip's normal")
    print(f"  resistance by Stage B's rule, 1 percent of the coil current through it while charging: {R_s:.1f} ohm ({V_ch/R_s:.1f} A, {V_ch**2/R_s:.0f} W into the cold while the charge lasts).")
    print(f"  The free-wheel diode D2: {V_d:.1f} V at 70 K.  The strip's closing time: about 1 s with its own thermal link, 10 to 60 s for MRI-type switches (Ogawa 2014).")
    print()
    print(f"  Order (a), Paul's: the supply is removed, the coil free-wheels through D2, the strip is cooled.")
    print(f"  {'step':>44} {'V across strip':>14} {'I in strip':>12} {'I in D2':>10} {'coil current change, ten together':>34} {'one alone':>22}")
    print(f"  {'1. charging, strip warm (heater on), leads in':>44} {V_ch:11.0f} V {V_ch/R_s:10.1f} A {'0 (blocks)':>10} {'rising at 1 %/s':>34} {'rising, V / L_leak':>22}")
    print(f"  {'2. leads out: the coil free-wheels in D2':>44} {V_d:11.1f} V {V_d/R_s:10.3f} A {I0:8.0f} A {'':>34} {'':>22}")
    for t in t_close:
        print(f"  {f'   ... for a closing time of {t:.0f} s':>44} {'':>14} {'':>12} {'':>10} {f'-{V_d*t/Lc:.2f} A = {V_d*t/Lc/I0*100:.2f} %':>34} {(f'-{V_d*t/Ll:.0f} A = {V_d*t/Ll/I0*100:.0f} %' if V_d*t/Ll < I0 else 'all of it, before ' + f'{Ll*I0/V_d:.0f} s'):>22}")
    print(f"  {'3. the strip goes superconducting':>44} {'0':>14} {f'-> {I0:.0f} A':>12} {'-> 0':>10} {'none (the transfer is in the strip-D2 mesh, uH)':>34} {'none':>22}")
    print(f"  The allowance to charge to is the loss of step 2: negligible with the ten coils together ({V_d*60/Lc/I0*100:.2f} percent at a minute), but a lone trim loses")
    print(f"  {V_d/Ll:.0f} A per second, {V_d/Ll/I0*100:.1f} percent of its current per second, and needs the strip closed within a second or so, or order (b).  The supply's regulation")
    print(f"  does not enter (a): the supply is gone before the strip closes.")
    print()
    print(f"  Order (b): the supply holds the current while the strip cools, then ramps down.")
    print(f"  {'step':>44} {'V across strip':>14} {'I in strip':>12} {'I in D2':>10} {'coil current, ten together':>26} {'one alone':>22}")
    print(f"  {'1. charging, strip warm, leads in':>44} {V_ch:11.0f} V {V_ch/R_s:10.1f} A {'0':>10} {'rising':>26} {'rising':>22}")
    print(f"  {'2. hold at I0, current-regulated':>44} {'~0 (joint drops, mV)':>14} {'~0':>12} {'0':>10} {f'{I0:.0f} A, steady':>26} {'the same':>22}")
    print(f"  {'3. the strip goes superconducting':>44} {'0':>14} {'0':>12} {'0':>10} {'unchanged: no voltage, no change':>26} {'the same':>22}")
    print(f"  {'4. the supply ramps to zero over t_r':>44} {'0':>14} {f'0 -> {I0:.0f} A':>12} {'0':>10} {'unchanged (flux conserved)':>26} {'the same':>22}")
    print(f"  {'5. leads out':>44} {'0':>14} {f'{I0:.0f} A':>12} {'0':>10} {'persistent':>26} {'persistent':>22}")
    print(f"  No allowance at all, and the same for ten together and one alone, because the coil never sees a voltage after step 1: the current moves from the supply")
    print(f"  into the strip by flux conservation.  The cost: the leads carry {I0:.0f} A for the closing time, {2*be.lead_heat(I0)[0]:.0f} W for the pair, {2*be.lead_heat(I0)[0]*60/1e3:.1f} kJ at a minute; and the strip must")
    print(f"  carry the whole {I0:.0f} A superconducting (a 3 mm strip does not: it is a 12 mm strip or several in parallel, in the low field at the coil's end).")
    print()
    q = be.lead_heat(I0)[0]
    print(f"  Paul's concern, the supply regulating voltage when the strip closes.  In (b) the strip is a short across the supply's cold terminals the moment it goes")
    print(f"  superconducting.  A current-regulated supply sees its current unchanged (the coil holds it) and its voltage fall to the joint drops: nothing happens.  A")
    print(f"  voltage-regulated supply at V_set drives V_set / R_leads into the strip at once (the strip-lead mesh is microhenries): at {V_ch:.0f} V and 0.1 ohm of leads")
    print(f"  {V_ch/0.1:,.0f} A, limited only by the supply's current limit, past the strip's critical current in microseconds; the strip quenches back to {R_s:.1f} ohm, takes")
    print(f"  {V_ch**2/R_s:.0f} W until the supply is turned down, and the quench spot has the supply's energy in it.  The coil is not harmed (its current cannot change")
    print(f"  without a volt-second it does not get), the strip may be.  With ten coils together each strip does the same on its own supply; one alone the same.")
    print(f"  So (b) needs a current-regulated supply, or a voltage-regulated one set to the lead drop ({I0*0.1:.0f} V at 0.1 ohm, which is current regulation by another name),")
    print(f"  from before the strip closes to the end of the ramp; and the ramp-down in step 4 is a current ramp, dI/dt set by the supply, with no voltage on the coil.")
    print(f"  (a) is immune to the supply because the supply is gone; its price is the allowance and the lone-trim limit.  Recommendation: (a) for the bundle's ten")
    print(f"  coils charged together, (b) for a lone trim, with the supply in current regulation throughout.")
    return dict(I0=I0, L_common=Lc, L_leak=Ll, V_ch=V_ch, R_s=R_s)


def parse_points(s):
    pts = []
    for item in s.split(";"):
        kind, rest = item.split(":", 1)
        nums, label = rest.split(":") if ":" in rest else (rest, "opt135")
        Wp, T, Wl = (float(x) for x in nums.split(","))
        pts.append((kind, Wp, T, Wl, label))
    return pts


def main():
    part = "all"; jobs = 6; h = 0.01; kind = "both"; grid = "coarse"; points = None; ks = (1.5,)
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a.startswith("--h="): h = float(a[4:])
        elif a.startswith("--kind="): kind = a[7:]
        elif a.startswith("--grid="): grid = a[7:]
        elif a.startswith("--points="): points = a[9:]
        elif a.startswith("--ks="): ks = tuple(float(x) for x in a[5:].split(","))
        elif a in ("-h", "--help"): print(__doc__); return
    t0 = time.time()
    hr("THE BEARING, ROUND 5 (2026-11-13)")
    print(f"  Ring at 700 km, rising quadratic: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m with hardware {ring.M_HARDWARE_M:,.0f}, side {ring.CABLE_SIDE:.2f} m; the rise {TAU*1e3:.0f} ms, the margin 1.35.")
    if part == "sweep":
        pts = None
        if points:
            pts = [(k, Wp, T, Wl) for k, Wp, T, Wl, lab in parse_points(points)]
        part_sweep(h=h, kind=kind, grid=grid, jobs=jobs, points=pts, ks=ks)
    elif part == "pick":
        part_pick(h=h)
    elif part == "final":
        part_final(parse_points(points), h=h if h < 0.01 else 0.0025)
    else:
        parts = dict(stores=part_stores, switches=part_switches, tapeloss=part_tapeloss, coldload=part_coldload, cases=part_cases, charging=part_charging)
        for name in (parts if part == "all" else [part]):
            parts[name]()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
