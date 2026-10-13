#!/usr/bin/env python3
"""
Mass of the iron-core bearing against the cable it costs, and the air-core alternative (loop ratings, requirement sets, dynamics of the lightest design).
Book: The Orbital Ring, Volume IV: Chapter 8 (named; the field solver and the mass argument).
Run: python bearing_mass.py [--parts=ABCD] [--resolve] (reads data/bearing_sweep.json and data/bearing_air.json)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

Bearing mass, and an air-core alternative (Paul, 2026-10-12, second task of the day).

Parts (Claude_Code_instructions_bearing_mass_2026-10-12.md):
  A  mass as the objective with the ring-mass loop closed (ring_altitude.cable_mass_structural): the 2026-10-12 report design first
  B  air-core loops as wide as the track, with a saturating plate: pull against ampere-turns for loop widths and plate
     thicknesses, winding field, tape, variants (alternating loops, shaped plate, back iron), stray field, stored energy, ripple
  C  casing-side and cable-side mass for the best iron-core and the best air-core arrangement under four requirement sets
  D  the dynamics (bearing_design.part_F) for the lightest design of set 1 and of set 3

Usage:
    python bearing_mass.py                 # everything (the air-core solves are cached in output/bearing_air.json)
    python bearing_mass.py --parts=AB
    python bearing_mass.py --resolve       # discard the air-core cache

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import levitation_config as lc
import sled_passage_simulation as sp
import sled_passage_pulsed as pz
import compensation_fields as cf
import magnetostatic_2d as ms
import bearing_design as bd

MU0 = lc.MU0
AIR_JSON = os.path.join(ring.DATA_DIR, "bearing_air.json")   # Part B air-core cache (24 min to redo): --resolve recomputes it
G_NET = ring.G_NET
M_LIM_PLATES = 1872.0                 # kg/m, the LIM reaction plates (Chapter 3), kept
M_OLD_PLATES = 2204.0                 # kg/m, the two 100 mm levitation plates of Chapter 3, replaced
M_CASING_BASE = ring.M_LOAD_M         # 12,000 kg/m; no levitation allowance is identified in the book [UNVERIFIED]
RHO_LAM, RHO_PLATE = bd.RHO_LAM, bd.RHO_PLATE
TAPE_T, TAPE_RHO = bd.TAPE_T, bd.TAPE_RHO
CRYO_FACTOR = 1.0                     # casing-side mass for formers, cryostat and supports, as a multiple of the tape mass [UNVERIFIED]
F_LIFT_PAIR = 393e3                   # the 30 km/s pull (both tracks of the pair together)
F_PARKED = bd.F_PARKED                # 396 kN/m
MARGIN_DECIDED = None                 # the support margin 'as decided' is the rated pull of the chosen geometry, set in Part A
PULL_CURVE = {1000.0: 283.6e3, 3000.0: 165.6e3, 5000.0: 62.8e3, 10000.0: 67.8e3, 20000.0: 249.8e3, 30000.0: 424.8e3}   # bearing_design.txt F.1
AIR_LEVELS = (50e3, 100e3, 150e3, 200e3, 300e3, 400e3, 500e3, 700e3, 1000e3)
B_WIND_LIMIT = 1.0                    # T, mean field at the winding that caps an air-core rating (12 mm tape Ic 316 A there)
H_AIR = 0.01


def hr(t=""):
    pz.hr(t)


# =============================================================================
# THE RING MASS LOOP
# =============================================================================

def ring_masses(m_casing, m_hw):
    m_s = ring.cable_mass_structural(m_casing, m_hw)
    return dict(casing=m_casing, hw=m_hw, structural=m_s, cable=m_s + m_hw, tension=ring.SIGMA_CABLE * m_s / ring.RHO_CNT,
                weight=m_casing * G_NET, side=math.sqrt(m_s / ring.RHO_CNT))


def print_ring(label, rm):
    print(f"  {label:>58}: casing {rm['casing']:8,.0f} kg/m, hardware {rm['hw']:7,.0f}, cable {rm['cable']:9,.0f} kg/m ({rm['cable']/119592:.2f} x), "
          f"tension {rm['tension']/1e9:6.0f} GN, casing weight {rm['weight']/1e3:6.1f} kN/m, cable side {rm['side']:.2f} m")


# =============================================================================
# IRON-CORE CATALOG (from the bearing_design sweep)
# =============================================================================

def iron_catalog():
    with open(bd.SWEEP_JSON, encoding="utf-8") as f:
        recs = json.load(f)
    cat = []
    for r in recs:
        if r["gm"]["plate_split"]:
            continue
        m = r["meta"]
        cat.append(dict(kind="iron", rec=r, label=bd.geom_label(m), W=round(m["W"], 2), plate_t=m["plate_t"], rating=r["rate"]["50"]["F"],
                        lam=m["iron_casing_m2"] * RHO_LAM, plate=m["iron_cable_m2"] * RHO_PLATE, plate_w=m["plate_w"], n_groups=m["n_groups"]))
    return cat


def iron_tape(entry, F, theva, w_tape):
    """Tape per meter of ring and winding field for an iron geometry at pull F (half Ic at the winding field)."""
    rec = entry["rec"]
    NI = bd.NI_for_force(rec, F)
    if NI != NI:
        return None
    Bw = math.hypot(bd.curve_interp(rec, NI, "Bx_wm"), bd.curve_interp(rec, NI, "By_wm"))
    Ic = bd.ic_tape(theva, Bw, w_tape); I_op = 0.5 * Ic
    tape = 2 * NI * entry["n_groups"] / I_op
    return dict(NI=NI, Bw=Bw, Ic=Ic, I_op=I_op, tape=tape, mass=tape * TAPE_T * w_tape * TAPE_RHO, W_mag=bd.curve_interp(rec, NI, "W_mag"))


# =============================================================================
# PART B: AIR-CORE LOOPS
# =============================================================================

def build_air(W=1.0, plate_t=0.20, NI=1e5, bundle=0.10, h=H_AIR, n_loops=1, plate_shape=None, back_iron=0.0, gap=0.10,
              plate_over=0.20, domain=None):
    """n_loops loops of width W side by side with alternating polarity (so touching legs add); each leg a bundle x bundle
    conductor whose top face is at the gap; the plate above; optional shaped plate (plate_shape = thickness over the
    legs, plate_t between them) and a back-iron strip of thickness back_iron behind each bundle (lamination steel)."""
    legs = []
    x0 = -n_loops * W / 2
    for k in range(n_loops):
        sgn = 1 if k % 2 == 0 else -1
        legs.append((x0 + k * W, sgn * NI)); legs.append((x0 + (k + 1) * W, -sgn * NI))
    # merge touching legs
    merged = {}
    for x, a in legs:
        merged[round(x, 6)] = merged.get(round(x, 6), 0.0) + a
    legs = sorted(merged.items())
    plate_w = n_loops * W + 2 * plate_over
    Lx, Ly = domain or (max(4.0, plate_w + 2.6), max(3.0, 2 * (gap + max(plate_t, plate_shape or 0) + bundle + back_iron) + 2.0))
    g = ms.Grid(Lx, Ly, h)
    if plate_shape:
        # thicker between the legs: plate_t over the legs and the overhang, plate_shape between
        g.add_iron(-plate_w / 2, plate_w / 2, gap, gap + plate_t, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
        for (xa, _), (xb, _) in zip(legs[:-1], legs[1:]):
            g.add_iron(xa + bundle / 2, xb - bundle / 2, gap + plate_t, gap + plate_shape, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
    else:
        g.add_iron(-plate_w / 2, plate_w / 2, gap, gap + plate_t, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
    boxes = []
    for x, a in legs:
        g.add_current(x - bundle / 2, x + bundle / 2, -bundle, 0.0, a)
        boxes.append((x - bundle / 2, x + bundle / 2, -bundle, 0.0))
        if back_iron > 0:
            g.add_iron(x - bundle / 2 - 0.05, x + bundle / 2 + 0.05, -bundle - back_iron, -bundle, cf.YOKE_MU, B_sat=cf.YOKE_B_SAT)
    meta = dict(W=W, n_loops=n_loops, plate_t=plate_t, plate_shape=plate_shape, back_iron=back_iron, bundle=bundle, gap=gap, plate_w=plate_w, legs=legs,
                plate_m2=plate_w * plate_t + (sum((xb - xa - bundle) for (xa, _), (xb, _) in zip(legs[:-1], legs[1:])) * (plate_shape - plate_t) if plate_shape else 0.0),
                back_iron_m2=len(legs) * (bundle + 0.1) * back_iron, NI_legs=sum(abs(a) for _, a in legs) / 2 / NI)
    return g, boxes, meta


def evaluate_air(g, boxes, meta, NI):
    gap, plate_t, plate_w = meta["gap"], meta["plate_t"], meta["plate_w"]
    h = g.h
    top = gap + max(plate_t, meta["plate_shape"] or 0)
    F = g.force_on_box(-plate_w / 2 - 0.1, plate_w / 2 + 0.1, gap - 2 * h, top + 2 * h)
    xs, By, Bx = g.line_By(gap - 2 * h, -plate_w / 2, plate_w / 2)
    wf = [g.region_field(*b) for b in boxes]
    m_pl = g.box(-plate_w / 2, plate_w / 2, gap, top) & (g.mu_r != 1.0)
    Bm = np.sqrt(g.Bx[m_pl]**2 + g.By[m_pl]**2)
    # sideways flux: the largest mean |Bx| through the plate thickness over cuts between the legs and beyond them (at the
    # midline between two antiparallel legs the field is vertical by symmetry, so the flux of each leg flows out to both sides)
    j_lo = int(round((gap - g.y0) / h)) + 1; j_hi = int(round((top - g.y0) / h)) - 1
    cuts = []
    for x, _ in meta["legs"]:
        cuts += [x - meta["bundle"], x + meta["bundle"], x - 0.25 * meta["W"], x + 0.25 * meta["W"]]
    B_side = 0.0
    for xc in cuts:
        i = int(round((xc - g.x0) / h))
        if 0 <= i < g.nx:
            B_side = max(B_side, float(np.abs(g.Bx[i, j_lo:j_hi]).mean()))
    sat = float(np.mean(g.mu_r[m_pl] < 0.1 * cf.PLATE_MU)) if m_pl.any() else 0.0
    return dict(NI=NI, F=F, B_peak=float(np.abs(By).max()), B_mean=float(np.abs(By).mean()), B_plate_mid=B_side,
                B_track=float(np.abs(By).mean()), B_pole=float(np.abs(By).max()), plate_sat_frac=sat,
                B_plate_max=float(Bm.max()), sat_frac=sat, W_mag=float(0.5 * np.sum(g.A * g.J) * h * h),
                Bx_wm=float(np.mean([w["Bx_mean"] for w in wf])), By_wm=float(np.mean([w["By_mean"] for w in wf])),
                Bx_w=max(w["Bx_max"] for w in wf), By_w=max(w["By_max"] for w in wf))


def air_curve(levels=AIR_LEVELS, verbose=False, **geom):
    g, boxes, meta = build_air(NI=1e5, **geom)
    J_ref = g.J.copy()
    g.solve()
    lin = evaluate_air(g, boxes, meta, 1e5); lin["nonlinear"] = False
    out = [lin]
    for NI in levels:
        g.J = J_ref * (NI / 1e5)
        top = meta["gap"] + max(meta["plate_t"], meta["plate_shape"] or 0)
        g.solve_nonlinear(force_box=(-meta["plate_w"] / 2 - 0.1, meta["plate_w"] / 2 + 0.1, meta["gap"] - 2 * g.h, top + 2 * g.h))
        r = evaluate_air(g, boxes, meta, NI); r["nonlinear"] = True; r["F_lin"] = lin["F"] * (NI / 1e5) ** 2; r["outer_dF"] = g.outer_dF
        out.append(r)
        if verbose:
            print(f"      {NI/1e3:4.0f} kA: F {r['F']/1e3:6.1f} kN/m (linear {r['F_lin']/1e3:6.0f}), plate peak {r['B_peak']:.2f} T, sideways {r['B_plate_mid']:.2f} T, "
                  f"sat {r['sat_frac']*100:3.0f} %, winding {r['Bx_wm']:.2f}/{r['By_wm']:.2f} (max {r['Bx_w']:.2f}/{r['By_w']:.2f}) T, W {r['W_mag']/1e3:6.1f} kJ/m (dF {r['outer_dF']:.4f})", flush=True)
    return out, meta


def air_label(m):
    s = f"{m['n_loops']} loop{'s' if m['n_loops'] > 1 else ''} of {m['W']:.1f} m, plate {m['plate_t']*1e3:.0f} mm"
    if m["plate_shape"]:
        s += f" ({m['plate_shape']*1e3:.0f} mm between the legs)"
    if m["back_iron"]:
        s += f", back iron {m['back_iron']*1e3:.0f} mm"
    return s


def air_sweep(verbose=True):
    if os.path.exists(AIR_JSON):
        with open(AIR_JSON, encoding="utf-8") as f:
            return json.load(f)
    os.makedirs(ring.DATA_DIR, exist_ok=True)
    geoms = []
    for W in (0.5, 1.0, 1.5, 2.0):
        for pt in (0.10, 0.20, 0.30, 0.40):
            geoms.append(dict(W=W, plate_t=pt))
    geoms += [dict(W=1.0, plate_t=0.20, n_loops=2), dict(W=1.0, plate_t=0.20, n_loops=3), dict(W=0.5, plate_t=0.20, n_loops=2),
              dict(W=1.0, plate_t=0.10, plate_shape=0.30), dict(W=1.0, plate_t=0.15, plate_shape=0.40),
              dict(W=1.0, plate_t=0.20, back_iron=0.05), dict(W=1.0, plate_t=0.20, back_iron=0.10), dict(W=1.0, plate_t=0.30, back_iron=0.05)]
    recs = []
    t0 = time.time()
    for gm in geoms:
        t1 = time.time()
        curve, meta = air_curve(**gm)
        rt = bd.rate_curve(curve)
        # an air core has little iron to saturate: cap the rating where the mean winding field reaches B_WIND_LIMIT
        rows = [c for c in curve if c.get("nonlinear")]
        NIs = [c["NI"] for c in rows]; Bws = [math.hypot(c["Bx_wm"], c["By_wm"]) for c in rows]
        rt["limit"] = "knee"
        if max(Bws) >= B_WIND_LIMIT:
            NI_b = float(np.interp(B_WIND_LIMIT, Bws, NIs))
            if NI_b < rt["50"]["NI"]:
                d = {k: float(np.interp(NI_b, NIs, [c[k] for c in rows])) for k in ("F", "B_track", "B_pole", "B_plate_max", "plate_sat_frac", "W_mag", "Bx_wm", "By_wm", "Bx_w", "By_w")}
                d["NI"] = NI_b; rt["50"] = d; rt["limit"] = "winding 1.0 T"
        meta2 = {k: v for k, v in meta.items() if k != "legs"}
        recs.append(dict(gm=gm, meta=meta2, rate=rt, curve=curve, t=time.time() - t1))
        if verbose:
            r50 = rt["50"]
            print(f"    {air_label(meta):52} rated {r50['F']/1e3:6.1f} kN/m at {r50['NI']/1e3:5.0f} kA per leg ({rt['limit']}), max {rt['F_max']/1e3:6.1f}, W/Fg {r50['W_mag']/(r50['F']*0.1):4.1f}, "
                  f"winding {r50['Bx_wm']:.2f}/{r50['By_wm']:.2f} T, plate {meta['plate_m2']*RHO_PLATE:5.0f} kg/m, back iron {meta['back_iron_m2']*RHO_LAM:4.0f} kg/m ({time.time()-t1:.0f} s)", flush=True)
    if verbose:
        print(f"    air-core solves done, {len(recs)} geometries, {time.time()-t0:.0f} s")
    with open(AIR_JSON, "w", encoding="utf-8") as f:
        json.dump([{k: v for k, v in r.items() if k != "t"} for r in recs], f, indent=1, default=float)   # without the solve times, so the tracked cache is stable
    return recs


def air_tape(rec, F, theva, w_tape):
    rows = [r for r in rec["curve"] if r.get("nonlinear")]
    Fs = [r["F"] for r in rows]; NIs = [r["NI"] for r in rows]
    if F > max(Fs):
        return None
    NI = float(np.interp(F, Fs, NIs))
    Bw = math.hypot(float(np.interp(NI, NIs, [r["Bx_wm"] for r in rows])), float(np.interp(NI, NIs, [r["By_wm"] for r in rows])))
    Ic = bd.ic_tape(theva, Bw, w_tape); I_op = 0.5 * Ic
    n_legs = rec["meta"]["NI_legs"] * 2                    # leg-equivalents carrying NI each
    tape = NI * n_legs / I_op
    return dict(NI=NI, Bw=Bw, Ic=Ic, I_op=I_op, tape=tape, mass=tape * TAPE_T * w_tape * TAPE_RHO,
                W_mag=float(np.interp(NI, NIs, [r["W_mag"] for r in rows])), B_peak=float(np.interp(NI, NIs, [r["B_peak"] for r in rows])),
                Bmax=math.hypot(float(np.interp(NI, NIs, [r["Bx_w"] for r in rows])), float(np.interp(NI, NIs, [r["By_w"] for r in rows]))))


def stray_field(W, plate_t, NI, h=0.05):
    """Linear coarse solve in a 30 x 30 m box: the field at 1, 3 and 10 m below the track (casing side) and sideways."""
    g, boxes, meta = build_air(W=W, plate_t=plate_t, NI=NI, h=h, domain=(30.0, 30.0))
    g.solve()
    out = {}
    for d in (1.0, 3.0, 10.0):
        bx, by = g.field_at(0.0, -d); out[f"below {d:.0f} m"] = math.hypot(bx, by)
        bx, by = g.field_at(d + W / 2, 0.0); out[f"beside {d:.0f} m"] = math.hypot(bx, by)
        bx, by = g.field_at(0.0, meta["gap"] + plate_t + d); out[f"above the plate {d:.0f} m"] = math.hypot(bx, by)
    return out


def part_B(theva):
    hr("PART B: AIR-CORE LOOPS AS WIDE AS THE TRACK (saturating plate, 10 mm cells, no yoke)")
    print(f"  Each leg of the loop is a {100:.0f} x 100 mm bundle whose top face is 100 mm under the plate; the two legs carry +NI and -NI; the plate is 0.4 m wider than")
    print(f"  the loop, sintered iron saturating at {cf.PLATE_B_SAT} T from mu_r {cf.PLATE_MU}.  Current stepped 50 to 500 kA per leg with the permeability carried along.")
    print(f"  Rating = the pull where the saturating force is half the linear force (as in bearing_design.py) or, if that comes first, where the mean field at")
    print(f"  the winding reaches {B_WIND_LIMIT:.1f} T (an air core has no yoke to saturate; the tape's own field is its limit).  'sideways' = the largest mean flux")
    print(f"  density along the plate through its thickness at cuts beside the legs and a quarter width in from them.  Hand estimate to check (air_core_wide_loop_estimate.txt): 95 kN/m from a 1 m loop at 222 kA per leg,")
    print(f"  plate 241 mm thick at 1.2 T, about 0.9 T at the winding, 23 kg/m of 12 mm tape.")
    t0 = time.time()
    recs = air_sweep(verbose=True)
    print(f"  ({time.time()-t0:.0f} s to run or load)")
    print()
    print(f"  {'geometry':>52} {'rated kN/m':>10} {'NI/leg kA':>9} {'limit':>13} {'plate peak T':>12} {'sideways T':>10} {'sat %':>5} {'W/Fg':>5} {'winding across / normal T':>26} {'corner max T':>12} {'curve max':>9} {'plate kg/m':>10} {'back iron':>9}")
    for r in sorted(recs, key=lambda r: -r["rate"]["50"]["F"]):
        m = r["meta"]; r50 = r["rate"]["50"]
        rows = [c for c in r["curve"] if c.get("nonlinear")]
        NIs = [c["NI"] for c in rows]
        pk = float(np.interp(r50["NI"], NIs, [c["B_peak"] for c in rows])); mid = float(np.interp(r50["NI"], NIs, [c["B_plate_mid"] for c in rows]))
        sat = float(np.interp(r50["NI"], NIs, [c["sat_frac"] for c in rows]))
        print(f"  {air_label(m):>52} {r50['F']/1e3:10.1f} {r50['NI']/1e3:9.0f} {r['rate']['limit']:>13} {pk:12.2f} {mid:10.2f} {sat*100:5.0f} {r50['W_mag']/(r50['F']*0.1):5.1f} "
              f"{r50['Bx_wm']:12.2f} / {r50['By_wm']:4.2f} {r50['Bx_w']:5.2f}/{r50['By_w']:5.2f} {r['rate']['F_max']/1e3:9.0f} {m['plate_m2']*RHO_PLATE:10.0f} {m['back_iron_m2']*RHO_LAM:9.0f}")
    # pull against ampere-turns for the plain loops
    print()
    print("  Pull against ampere-turns (kN/m) for the plain loops:")
    print(f"  {'loop, plate':>22} " + " ".join(f"{NI/1e3:>7.0f}" for NI in AIR_LEVELS) + "   kA per leg")
    for r in recs:
        m = r["meta"]
        if m["n_loops"] != 1 or m["plate_shape"] or m["back_iron"]:
            continue
        rows = [c for c in r["curve"] if c.get("nonlinear")]
        print(f"  {m['W']:.1f} m, {m['plate_t']*1e3:3.0f} mm plate   " + " ".join(f"{c['F']/1e3:7.1f}" for c in rows))
    # tape at the hand estimate's point and at the ratings
    print()
    print("  Tape at half Ic at the winding field (THEVA 70 K): the 1 m loop with the 200 mm plate at 95 kN/m (the hand estimate's case), 12 mm tape;")
    print("  and each plain loop at its rating, 12 mm (support) and 3 mm (pulsed):")
    r1 = next(r for r in recs if r["meta"]["W"] == 1.0 and r["meta"]["plate_t"] == 0.20 and r["meta"]["n_loops"] == 1 and not r["meta"]["plate_shape"] and not r["meta"]["back_iron"])
    t = air_tape(r1, 95e3, theva, 0.012)
    if t:
        print(f"    1 m loop, 200 mm plate, 95 kN/m: {t['NI']/1e3:.0f} kA per leg (hand estimate 222), winding {t['Bw']:.2f} T (hand 0.9), Ic {t['Ic']:.0f} A, {t['I_op']:.0f} A per tape, "
              f"{t['tape']:,.0f} m of tape per m, {t['mass']:.1f} kg/m (hand 23), plate peak {t['B_peak']:.2f} T (hand 0.44)")
    print(f"  {'geometry':>52} {'rated kN/m':>10} {'12 mm: Ic A':>11} {'tape m/m':>8} {'kg/m':>6} {'3 mm: Ic A':>10} {'tape m/m':>8} {'kg/m':>6}")
    for r in sorted(recs, key=lambda r: -r["rate"]["50"]["F"]):
        F = r["rate"]["50"]["F"]
        t12 = air_tape(r, F, theva, 0.012); t3 = air_tape(r, F, theva, 0.003)
        if t12 and t3:
            print(f"  {air_label(r['meta']):>52} {F/1e3:10.1f} {t12['Ic']:11.0f} {t12['tape']:8,.0f} {t12['mass']:6.1f} {t3['Ic']:10.0f} {t3['tape']:8,.0f} {t3['mass']:6.1f}")
    # stray field
    print()
    print("  Stray field (linear, 50 mm cells, 30 x 30 m box) for the 1 m loop with the 200 mm plate at 200 kA per leg:")
    sf = stray_field(1.0, 0.20, 200e3)
    for k, v in sf.items():
        print(f"    {k:>22}: {v*1e3:8.3f} mT")
    print("  The loop's far field is that of a line dipole and falls as 1 over distance squared; the plate shorts the near field above it.")
    # ripple
    print()
    print("  Ripple along the direction of travel.  At a coil end the loop's whole ampere-turns cross from one leg to the other, so with n coils per")
    print("  bundle, each of NI/n and ends staggered, the plate field dips by 1/n over about half the loop width [UNVERIFIED].  Eddy heating per dip and")
    print("  per meter of ring (two tracks), 1 m loop, 200 mm plate at 95 kN/m, 10 coils per bundle:")
    if t:
        dB = t["B_peak"] / 10
        E_dip = bd.eddy_cycle(dB, 0.02 / pz.V_REL, 0.5 / pz.V_REL, d_p=0.20)
        print(f"    dip {dB:.3f} T, eddy energy per dip {E_dip*1e3:.3f} mJ per m^2 of plate over the 1.4 m plate width")
        for Lc in (200.0, 1000.0, 5000.0):
            dips = 2 * 10 / Lc
            print(f"    coil length {Lc:5.0f} m: {dips:.3f} dips per m, {2*E_dip*1.4*pz.V_REL*dips:.2f} W/m")
    return recs


# =============================================================================
# PART C: REQUIREMENTS AGAINST MASS
# =============================================================================

def pick_iron(cat, F, theva, w_tape, W_max=99.0):
    """The iron geometry with the least ring cost (10 x casing-side + cable-side mass) whose rating covers F; widths under 1 m are
    scaled from the 1 m geometries (rating and laminations proportional to width, plate to width + 0.4) [UNVERIFIED scaling]."""
    best = None
    for e in cat:
        if e["rating"] < F and e["W"] < 2.0:
            continue
        scale = max(F / e["rating"], 0.4) if e["W"] <= 1.0 else (F / e["rating"] if (e["W"] >= 2.0 and F > e["rating"]) else 1.0)
        if e["W"] * scale > W_max + 1e-6:
            continue
        tp = iron_tape(e, F / scale, theva, w_tape)
        if tp is None:
            continue
        lam = e["lam"] * scale; plate = e["plate_t"] * (e["W"] * scale + 0.4) * RHO_PLATE
        tape_mass = tp["mass"] * scale
        cost = 10 * (lam + tape_mass * (1 + CRYO_FACTOR)) + plate
        d = dict(entry=e, scale=scale, W=e["W"] * scale, lam=lam, plate=plate, tape=tape_mass, tape_len=tp["tape"] * scale, cost=cost, NI=tp["NI"], Bw=tp["Bw"],
                 I_op=tp["I_op"], plate_t=e["plate_t"], label=f"{e['label']}" + (f" scaled to {e['W']*scale:.2f} m" if abs(scale - 1) > 1e-6 else ""), rating=e["rating"] * scale,
                 Bmax=math.hypot(bd.curve_interp(e["rec"], tp["NI"], "Bx_w"), bd.curve_interp(e["rec"], tp["NI"], "By_w")))
        if best is None or cost < best["cost"]:
            best = d
    return best


def pick_air(recs, F, theva, w_tape, W_max=99.0):
    best = None
    for r in recs:
        if r["rate"]["50"]["F"] < F or r["meta"]["n_loops"] * r["meta"]["W"] > W_max + 1e-6:
            continue
        tp = air_tape(r, F, theva, w_tape)
        if tp is None:
            continue
        m = r["meta"]
        plate = m["plate_m2"] * RHO_PLATE; back = m["back_iron_m2"] * RHO_LAM
        cost = 10 * (back + tp["mass"] * (1 + CRYO_FACTOR)) + plate
        d = dict(rec=r, W=m["n_loops"] * m["W"], lam=back, plate=plate, tape=tp["mass"], tape_len=tp["tape"], cost=cost, NI=tp["NI"], Bw=tp["Bw"], I_op=tp["I_op"],
                 plate_t=m["plate_t"], label=air_label(m), rating=r["rate"]["50"]["F"], Bmax=tp["Bmax"], W_mag=tp["W_mag"])
        if best is None or cost < best["cost"]:
            best = d
    return best


REQ_SETS = {
    1: "as decided: each support track carries the casing alone with the decided margin; either compensation track does the pair's whole job alone",
    2: "support as decided; compensation pairs share their job (one track out limits the launch speed and the parked mass)",
    3: "as 2, the load compensation moved onto the support tracks as a second driven winding: four tracks",
    4: "as 3, support margin 1.5 on one track alone",
}


def requirements(set_no, weight, margin_decided):
    """Per-track force requirements for a requirement set at a casing weight: (support one-track, lift per track, load per track,
    n_support, n_lift, n_load).  The support requirement in set 3 and 4 includes the load share it must carry."""
    if set_no == 1:
        return dict(support=margin_decided * weight, lift=F_LIFT_PAIR, load=max(F_PARKED, F_LIFT_PAIR), n_lift=2, n_load=2, speed_limit=None)
    if set_no == 2:
        return dict(support=margin_decided * weight, lift=F_LIFT_PAIR / 2, load=max(F_PARKED, F_LIFT_PAIR) / 2, n_lift=2, n_load=2, speed_limit=True)
    if set_no == 3:
        # the support track alone must carry the weight with the decided margin, and in pair mode the weight share plus the load share
        return dict(support=max(margin_decided * weight, weight / 2 + F_PARKED / 2), lift=F_LIFT_PAIR / 2, load=0.0, n_lift=2, n_load=0, speed_limit=True)
    return dict(support=max(1.5 * weight, weight / 2 + F_PARKED / 2), lift=F_LIFT_PAIR / 2, load=0.0, n_lift=2, n_load=0, speed_limit=True)


def close_loop(set_no, picker, theva, margin_decided, verbose=False):
    """Iterate casing mass -> cable -> casing weight -> bearing forces -> bearing mass until the casing mass settles."""
    m_casing = M_CASING_BASE
    for it in range(30):
        weight = m_casing * G_NET
        req = requirements(set_no, weight, margin_decided)
        # width budget: the cable face is CABLE_SIDE wide; each plate is 0.4 m wider than its track; two support tracks and the load
        # tracks share the underside, the lift tracks the top.  Try every split of the underside and keep the cheapest that fits.
        face = ring.CABLE_SIDE
        lift = picker(req["lift"], theva, 0.003, face / 2 - 0.4) if req["n_lift"] else None
        best = None
        for W_s_cap in (0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
            W_l_cap = (face - 2 * (W_s_cap + 0.4)) / 2 - 0.4 if req["n_load"] else 0.0
            if req["n_load"] and W_l_cap < 0.5:
                continue
            if 2 * (W_s_cap + 0.4) > face + 1e-6:
                continue
            sup_c = picker(req["support"], theva, 0.012, W_s_cap)
            load_c = picker(req["load"], theva, 0.003, W_l_cap) if req["n_load"] else None
            if sup_c is None or (req["n_load"] and load_c is None):
                continue
            cost = 2 * sup_c["cost"] + (2 * load_c["cost"] if load_c else 0.0)
            if best is None or cost < best[0]:
                best = (cost, sup_c, load_c)
        if best is None or (req["n_lift"] and lift is None):
            return None
        _, sup, load = best
        # set 3 and 4: the support track carries a second, driven 3 mm winding for the load share (tape only, same poles)
        extra = 0.0
        if set_no in (3, 4):
            e2 = picker(F_PARKED / 2, theva, 0.003, 99.0)
            extra = e2["tape"] if e2 else 0.0
        casing_side = 2 * (sup["lam"] + (sup["tape"] + extra) * (1 + CRYO_FACTOR)) + req["n_lift"] * (lift["lam"] + lift["tape"] * (1 + CRYO_FACTOR)) \
            + (req["n_load"] * (load["lam"] + load["tape"] * (1 + CRYO_FACTOR)) if load else 0.0)
        cable_side = 2 * sup["plate"] + req["n_lift"] * lift["plate"] + (req["n_load"] * load["plate"] if load else 0.0)
        m_new = M_CASING_BASE + casing_side
        if abs(m_new - m_casing) < 1.0:
            break
        m_casing = m_new
    rm = ring_masses(m_new, M_LIM_PLATES + cable_side)
    return dict(req=req, sup=sup, lift=lift, load=load, casing_side=casing_side, cable_side=cable_side, ring=rm, iters=it + 1, extra_tape=extra)


def speed_limits(F_one):
    """With one compensation track out, the speeds the surviving track covers alone (pull curve of bearing_design.txt F.1) and the parked mass."""
    us = sorted(PULL_CURVE); Fs = [PULL_CURVE[u] for u in us]
    lo = None; hi = None
    for u, F in zip(us, Fs):
        if F <= F_one and lo is None:
            lo = u
    # upper: first speed above 10 km/s where the pull exceeds F_one, interpolated
    for u0, u1 in zip(us[:-1], us[1:]):
        if u0 >= 10000 and PULL_CURVE[u0] <= F_one < PULL_CURVE[u1]:
            hi = u0 + (F_one - PULL_CURVE[u0]) / (PULL_CURVE[u1] - PULL_CURVE[u0]) * (u1 - u0)
    # lower: between 1 and 3 km/s the pull falls with speed
    if PULL_CURVE[1000.0] > F_one >= PULL_CURVE[3000.0]:
        lo = 1000.0 + (PULL_CURVE[1000.0] - F_one) / (PULL_CURVE[1000.0] - PULL_CURVE[3000.0]) * 2000.0
    elif F_one >= PULL_CURVE[1000.0]:
        lo = 1000.0
    return lo, hi if hi else (30000.0 if F_one >= PULL_CURVE[30000.0] else None), F_one / F_PARKED * 5000.0


def part_C(cat, air, theva, margin_decided):
    hr("PART C: REQUIREMENTS AGAINST MASS, RING LOOP CLOSED")
    print(f"  Casing = {M_CASING_BASE:,.0f} kg/m plus the bearing's casing side (laminations or back iron, tape, and {CRYO_FACTOR:.0f} x the tape mass for formers, cryostat and")
    print(f"  supports [UNVERIFIED]); hardware = {M_LIM_PLATES:,.0f} kg/m of LIM plates plus the bearing plates.  Each design is sized at the casing weight it")
    print(f"  produces (iterated).  Support margin 'as decided' = {margin_decided:.2f} on one track alone (the rated pull of the 2026-10-12 geometry over the")
    print(f"  12,000 kg/m weight), kept as the casing grows; set 4 uses 1.5.  Compensation: lift pair {F_LIFT_PAIR/1e3:.0f} kN/m at 30 km/s, load pair the parked")
    print(f"  {F_PARKED/1e3:.0f} kN/m.  Geometries are picked from the catalogs for the least ring cost, 10 x casing side + cable side.")
    results = {}
    print(f"  Sets 1 to 3 are run twice: with the decided support margin {margin_decided:.2f} and with 1.5, so the cost of the support margin and the cost of")
    print(f"  the compensation requirements can be told apart (set 4 is set 3 at 1.5 by definition).")
    for set_no, desc in REQ_SETS.items():
        for margin in ((margin_decided, 1.5) if set_no < 4 else (1.5,)):
            print()
            print(f"  Set {set_no} (support margin {margin:.2f}): {desc}")
            for kind, picker in (("iron", lambda F, th, w, Wm=99.0: pick_iron(cat, F, th, w, Wm)), ("air", lambda F, th, w, Wm=99.0: pick_air(air, F, th, w, Wm))):
                res = close_loop(set_no, picker, theva, margin)
                results[(set_no, kind, round(margin, 2))] = res
                if res is None:
                    print(f"    {kind}: no geometry in the catalog covers the requirement within the width of the cable face")
                    continue
                print_pick(kind, res, rq_speed=True)
    print()
    print(f"  {'set':>3} {'margin':>6} {'kind':>4} {'casing kg/m':>11} {'hardware':>8} {'cable':>9} {'x book':>6} {'tension GN':>10} {'casing side':>11} {'cable side':>10}")
    for (set_no, kind, margin), res in results.items():
        if res:
            rm = res["ring"]
            print(f"  {set_no:3d} {margin:6.2f} {kind:>4} {rm['casing']:11,.0f} {rm['hw']:8,.0f} {rm['cable']:9,.0f} {rm['cable']/119592:6.2f} {rm['tension']/1e9:10.0f} {res['casing_side']:11,.0f} {res['cable_side']:10,.0f}")
    return results


def print_pick(kind, res, rq_speed=True):
    if True:
        if True:
            rq, rm = res["req"], res["ring"]
            rq, rm = res["req"], res["ring"]
            print(f"    {kind:>4}: support {res['sup']['label']} at {rq['support']/1e3:.0f} kN/m per track ({res['sup']['NI']/1e3:.0f} kA, {res['sup']['I_op']:.0f} A, winding {res['sup']['Bw']:.2f} T mean, "
                  f"{res['sup']['Bmax']:.2f} T corner, {res['sup']['tape']:.0f} kg/m tape, {res['sup']['lam']:,.0f} kg/m iron, {res['sup']['plate']:,.0f} kg/m plate)")
            under = 2 * (res['sup']['W'] + 0.4) + (rq['n_load'] * (res['load']['W'] + 0.4) if res['load'] else 0.0); top = rq['n_lift'] * (res['lift']['W'] + 0.4) if res['lift'] else 0.0
            print(f"          width: underside {under:.2f} m of {ring.CABLE_SIDE:.2f} ({'fits' if under <= ring.CABLE_SIDE else 'DOES NOT FIT'}), top {top:.2f} m ({'fits' if top <= ring.CABLE_SIDE else 'DOES NOT FIT'})")
            if res["lift"]:
                print(f"          lift x{rq['n_lift']} {res['lift']['label']} at {rq['lift']/1e3:.0f} kN/m per track ({res['lift']['tape']:.0f} kg/m tape, {res['lift']['lam']:,.0f} iron, {res['lift']['plate']:,.0f} plate)")
            if res["load"]:
                print(f"          load x{rq['n_load']} {res['load']['label']} at {rq['load']/1e3:.0f} kN/m per track ({res['load']['tape']:.0f} kg/m tape, {res['load']['lam']:,.0f} iron, {res['load']['plate']:,.0f} plate)")
            if res["extra_tape"]:
                print(f"          second driven 3 mm winding on each support track for the load share: {res['extra_tape']:.0f} kg/m of tape")
            print(f"          casing side {res['casing_side']:,.0f} kg/m, cable side {res['cable_side']:,.0f} kg/m ({res['iters']} iterations)")
            print_ring("ring", rm)
            if rq["speed_limit"]:
                F_one = rq["lift"]
                lo, hi, mass = speed_limits(F_one)
                print(f"          with one compensation track out ({F_one/1e3:.0f} kN/m alone): launches covered from {lo/1e3 if lo else float('nan'):.1f} to {hi/1e3 if hi else float('nan'):.1f} km/s at 5,000 t; "
                      f"parked mass limited to {mass:,.0f} t")


# =============================================================================

def part_A(cat, theva):
    hr("PART A: MASS AS THE OBJECTIVE, THE RING LOOP CLOSED")
    print(f"  ring_altitude.cable_mass_structural(m_casing, m_hardware): the book's casing {ring.M_LOAD_M:,.0f} kg/m and hardware {ring.M_HARDWARE_M:,.0f} kg/m (two 100 mm levitation plates")
    print(f"  {M_OLD_PLATES:,.0f} kg/m plus LIM plates {M_LIM_PLATES:,.0f} kg/m, Chapter 3).  No allowance for the old levitation coils inside the 12,000 kg/m is identified in the book, so")
    print(f"  the bearing's casing side is added on top [UNVERIFIED].")
    print_ring("book", ring_masses(ring.M_LOAD_M, ring.M_HARDWARE_M))
    # the report design
    lam = 2 * 3213 + 4 * 5317; tape = 2 * 31.9 + 4 * 31.0; plates = 2 * 2204 + 4 * 4486
    print_ring("2026-10-12 report design, tape and iron as reported (no cryostat)", ring_masses(ring.M_LOAD_M + lam + tape, M_LIM_PLATES + plates))
    print_ring("plates only, casing unchanged", ring_masses(ring.M_LOAD_M, M_LIM_PLATES + plates))
    print_ring("laminations only, cable hardware unchanged", ring_masses(ring.M_LOAD_M + lam, ring.M_HARDWARE_M))
    # closed loop for the report design: support width scales with the weight, compensation fixed
    sup_rating_per_m = 349e3; lam_per_m = 3213.0; tape_per_m = 31.9
    m = ring.M_LOAD_M + lam + tape
    for it in range(40):
        weight = m * G_NET
        W_s = 3.67 * weight / sup_rating_per_m
        m_new = ring.M_LOAD_M + 2 * (lam_per_m * W_s + tape_per_m * W_s) + 4 * 5317 + 4 * 31.0
        if abs(m_new - m) < 1:
            break
        m = m_new
    plates_cl = 2 * 0.2 * (W_s + 0.4) * RHO_PLATE + 4 * 4486
    print_ring(f"report design closed (support tracks {W_s:.2f} m to keep 3.67 x)", ring_masses(m, M_LIM_PLATES + plates_cl))
    print(f"  Every kilogram on the casing costs {(ring_masses(13000, ring.M_HARDWARE_M)['cable'] - 119592)/1000:.1f} kg of cable (the cable scales with the load it must carry), a kilogram on the")
    print(f"  cable costs {(ring_masses(ring.M_LOAD_M, ring.M_HARDWARE_M + 1000)['cable'] - 119592)/1000:.2f} kg: iron belongs on the cable side, the casing side should carry tape and little else.")
    margin = 349e3 / (ring.M_LOAD_M * G_NET)
    return margin


def part_D(results, theva):
    hr("PART D: DYNAMICS FOR THE LIGHTEST DESIGN OF SET 1 AND OF SET 3 (bearing_design.part_F)")
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    for set_no in (1, 3):
        cands = [(k, r) for k, r in results.items() if k[0] == set_no and r is not None]
        if not cands:
            print(f"  set {set_no}: nothing to run"); continue
        margins = sorted({k[2] for k, r in cands}, reverse=True)
        cands = [(k, r) for k, r in cands if k[2] == margins[0]]
        (sn, kind, margin), res = min(cands, key=lambda kr: kr[1]["ring"]["cable"])
        print(f"  (set {set_no} at support margin {margin:.2f}, the highest margin with a design that fits)")
        rq = res["req"]
        print()
        print(f"  Set {set_no}, lightest: {kind}.  Caps: support one track {rq['support']/1e3:.0f} kN/m; lift {rq['lift']/1e3:.0f} per track; load {rq['load']/1e3:.0f} per track"
              + ("" if rq["n_load"] else " (the load pull is a second winding on the support tracks, capped at the support rating less the weight share)") + ".")
        W_s = res["sup"]["W"]; W_li = res["lift"]["W"]; W_lo = res["load"]["W"] if res["load"] else W_s
        load_cap = rq["load"] if rq["n_load"] else max(rq["support"] - rq["support"] / 2, 1.0)
        C = dict(F_one_cap=rq["support"], rec=dict(meta=dict(W=W_s)))
        D = dict(caps=dict(lift1=rq["lift"], load1=load_cap, lift2=2 * rq["lift"], load2=2 * load_cap),
                 lift=dict(L_coil=50.0, ramp=20e-3, W_track=W_li), load=dict(L_coil=25.0, W_track=W_lo))
        if not rq["n_load"]:
            D["load"]["W_track"] = -0.4          # no load plates: the support plates sit at the cable edge
        bd.part_F(C, D, tabs)


def main():
    parts = "ABCD"
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
        elif arg == "--resolve" and os.path.exists(AIR_JSON):
            os.remove(AIR_JSON)
    t_all = time.time()
    hr("BEARING MASS, AND AN AIR-CORE ALTERNATIVE (2026-10-12)")
    theva = bd.read_theva()
    cat = iron_catalog()
    margin = part_A(cat, theva)
    air = part_B(theva) if ("B" in parts or "C" in parts or "D" in parts) else []
    results = part_C(cat, air, theva, margin) if ("C" in parts or "D" in parts) else {}
    if "D" in parts:
        part_D(results, theva)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
