#!/usr/bin/env python3
"""
Design history: the spread-pull step of the sled passage work that led to the adopted bearing in `final_bearing.py`.
The spread pull (the sled's load spread over 50 ms) as the design case for the bearing dynamics.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (the 50 ms spread pull).
Run: python sled_passage_spread.py [see --help]
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

The spread pull as the design case (Paul, 2026-10-10).

The lift compensation track does not match the spacecraft's force instant by instant.
It applies the same impulse (the cable's share, 90.9 percent) spread over a
time T centered on the pass, with the coil current regulated against the
measured gap and a slow (1 Hz) gap loop on the casing support tracks.  Below the
neutral speed the load compensation track does the same with the opposite sign.  The
bare sled is cancelled in step by the same track at its own low level.

Parts (section 3 of Claude_Code_instructions_spread_pull_2026-10-10.md; names of section 8):
  A  three rules for T over the launch, one recommended
  B  does the spacecraft load length and position drop out
  C  the bare sled and the low-speed side
  D  tolerances
  E  failure cases (beam model)
  F  stability margin and regulation bandwidth
  G  power, energy, plates and the width budget
  H  the beam model at 5, 10, 20 and 30 km/s
  I  hysteresis and eddy current losses in the plates and coils, added to the budget of G

Usage:
    python sled_passage_spread.py                  # everything, graphs on
    python sled_passage_spread.py --rule=1         # 1 = fixed 50 ms, 2 = 10 mm target, 3 = fixed 1.5 km (default: auto)
    python sled_passage_spread.py --parts=AB --no-graphs
    python sled_passage_spread.py --help

Models and conventions as in sled_passage_pulsed.py.
Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import math
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ring_altitude as ring
import levitation_config as lc
import sled_passage_simulation as sp
import sled_passage_pulsed as pz

MU, MU_C, W, GAP0, BAND = pz.MU, pz.MU_C, pz.W, pz.GAP0, pz.BAND
L_SLED, M_SLED_M, M_SC0 = pz.L_SLED, pz.M_SLED_M, pz.M_SC0
P12 = pz.P_SAT[1.2]
MATCH = pz.MATCH_RING
PALETTE = sp.PALETTE
TAU0 = 2e-3
L_SECTION = 50.0
L_MAIN = 100.0
RULE_NAMES = {1: "fixed T = 50 ms", 2: "T for a 10 mm excursion, at least the pass", 3: "fixed energized length 1.5 km",
              4: "fixed T = 50 ms, at least the pass"}
SPEEDS = [float(u) for u in range(1000, 30001, 1000)]


def hr(t=""):
    pz.hr(t)


def J_zone(u, L_sc, M_sc=M_SC0):
    """Zone impulse per meter the pull must cancel (cable's share), N s/m, magnitude."""
    return MATCH * abs(float(pz.q_zone(u, L_sc, M_sc))) * L_sc / u


def rule_T(rule, u, L_sc, M_sc=M_SC0):
    if rule == 1:
        return 0.05
    if rule == 2:
        kick = J_zone(u, L_sc, M_sc) / MU
        return max(8 * MU * 0.010 / (J_zone(u, L_sc, M_sc)) if kick > 0 else L_sc / u, L_sc / u)
    if rule == 3:
        return 1500.0 / u
    if rule == 4:
        return max(0.05, L_sc / u)
    raise ValueError(rule)


def energized_length(u, T, L_sc):
    return u * T + L_sc + L_SECTION + 2 * u * TAU0


def spread_lane(u, L_sc, T, **kw):
    d = dict(L_sc=L_sc, tau=TAU0, L_section=L_SECTION, gap_dep=False, steady_loop_hz=1.0, T_spread=T, B_sat=1.2)
    d.update(kw)
    return pz.lane(**d)


def main():
    rule_arg = "auto"
    parts = "ABCDEFGHI"
    graphs = True
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--rule="):
            rule_arg = arg.split("=")[1]
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
        elif arg == "--no-graphs":
            graphs = False
        else:
            print(f"Unknown option {arg}"); return
    out_dir = ring.out_dir("graphs_sled_passage")   # under ORBITAL_RING_OUT (default runs)
    t_all = time.time()
    tabs = sp.limiter_tables()

    hr("THE SPREAD PULL AS THE DESIGN CASE")
    print(f"  5,000 t, sled {L_SLED:,.0f} m at {M_SLED_M:.0f} kg/m, casing {MU:,.0f} kg/m, cable {MU_C:,.0f} kg/m, plate field 1.2 T "
          f"({P12/1e3:.0f} kPa), {L_SECTION:.0f} m sections, {TAU0*1e3:.0f} ms rise time, pulse regulated against the gap, 1 Hz casing support loop,")
    print(f"  pulse cancels the cable's share {MATCH*100:.1f} % of the impulse, spacecraft at the front of the sled, zone {L_MAIN:g} m unless stated.")
    print(f"  Point model (one meter of casing and of cable, both free) for A to D, F, G; beam model (10 m truss, cable as a string) for E and H.")
    print(f"  Excursion = largest casing offset from the 100 mm gap during the pass; the band is 20 mm; the limiter must stay negligible.")

    # ------------------------------------------------------------------ A
    rows_A = {}
    if "A" in parts or True:
        hr("PART A: ONE RULE FOR T OVER THE LAUNCH (zone 100 m)")
        print(f"  Rules: 1 {RULE_NAMES[1]}; 2 {RULE_NAMES[2]}, T = 8 mu x 0.010 / kick speed; 3 {RULE_NAMES[3]}, T = 1,500 / u.")
        print(f"  force = peak pulsed force (zone pull plus the bare sled's level where they overlap); width at 1.2 T; L_on = u T + zone")
        print(f"  + section + 2 u tau; power = 573 kPa x 0.1 m x width x u; F_lim = peak limiter force.")
        t0 = time.time()
        rows_A = {1: [], 2: [], 3: [], 4: []}
        for u in SPEEDS:
            Ts = {k: rule_T(k, u, L_MAIN) for k in (1, 2, 3, 4)}
            res, _ = pz.point_pulsed([spread_lane(u, L_MAIN, Ts[k], label=f"rule{k}") for k in (1, 2, 3, 4)], u, tabs)
            for k, r in zip((1, 2, 3, 4), res):
                F = max(r["pk_hd"], r["pk_hu"])
                r.update(T=Ts[k], F=F, width=F / P12, L_on=energized_length(u, Ts[k], L_MAIN), P=P12 * GAP0 * (F / P12) * u, J=J_zone(u, L_MAIN))
                rows_A[k].append(r)
        print(f"  ({time.time()-t0:.0f} s)")
        for rule in (1, 2, 3, 4):
            print(f"\n  Rule {rule}: {RULE_NAMES[rule]}")
            print(f"  {'u km/s':>7} {'J kN s/m':>9} {'T ms':>7} {'F kN/m':>8} {'width m':>8} {'L_on m':>7} {'P GW':>7} {'exc mm':>7} {'F_lim kN/m':>10} {'band':>5}")
            for r in rows_A[rule]:
                print(f"  {r['u']/1e3:7.0f} {r['J']/1e3:9.2f} {r['T']*1e3:7.1f} {r['F']/1e3:8.0f} {r['width']:8.2f} {r['L_on']:7.0f} {r['P']/1e9:7.2f} "
                      f"{r['max_exc']*1e3:7.1f} {r['pk_lim']/1e3:10.2f} {('ok' if r['ok_band'] else 'FAIL'):>5}")
        summary = {}
        for rule in (1, 2, 3, 4):
            rr = rows_A[rule]
            ok = all(r["ok_band"] for r in rr)
            summary[rule] = dict(ok=ok, max_exc=max(r["max_exc"] for r in rr), max_w=max(r["width"] for r in rr),
                                 max_P=max(r["P"] for r in rr), max_Lon=max(r["L_on"] for r in rr),
                                 w_hd=max(r["width"] for r in rr if r["u"] > 7000), w_hu=max(r["width"] for r in rr if r["u"] < 7000))
            s = summary[rule]
            print(f"  Rule {rule}: {'all speeds inside the band' if ok else 'FAILS the band'}; largest excursion {s['max_exc']*1e3:.1f} mm; "
                  f"lift compensation width {s['w_hd']:.2f} m, load compensation width {s['w_hu']:.2f} m; peak power {s['max_P']/1e9:.2f} GW; longest energized {s['max_Lon']:,.0f} m.")
        if rule_arg == "auto":
            ok_rules = [k for k in (1, 2, 3, 4) if summary[k]["ok"]]
            rule = (min(ok_rules, key=lambda k: (round(summary[k]["max_Lon"] / 100.0), summary[k]["max_exc"])) if ok_rules
                    else min((1, 2, 3, 4), key=lambda k: summary[k]["max_exc"]))
        else:
            rule = int(rule_arg)
        print(f"\n  Chosen rule: {rule} ({RULE_NAMES[rule]}){' (set with --rule)' if rule_arg != 'auto' else ', the rule inside the band with the shortest energized track, then the smallest excursion'}.")
        print(f"  Rule 2 has no upper bound: as the load goes to zero at the neutral speed its T grows without limit, so it needs a cap to")
        print(f"  be a design rule; rule 3 spreads a low-speed pull too thin.  A fixed T with the pass as its floor (rule 4) is the")
        print(f"  one-sentence rule: the pull lasts 50 ms, or as long as the spacecraft zone takes to pass if that is longer.")
        print(f"  Single lift compensation width for the whole launch: {summary[rule]['w_hd']:.2f} m at 1.2 T; load compensation width {summary[rule]['w_hu']:.2f} m.")
        print(f"  Hand values at 30 km/s (first run): J 18.9 kN s/m; T = 50 ms: 378 kN/m, 0.66 m, 1,720 m, 1.1 GW, 9.0 mm; T = 100 ms: 189 kN/m, 0.33 m, 3,220 m, 0.57 GW, 18.8 mm.")
        for T in (0.05, 0.10):
            res, _ = pz.point_pulsed([spread_lane(30000.0, L_MAIN, T)], 30000.0, tabs)
            r = res[0]
            Jz = J_zone(30000.0, L_MAIN)
            print(f"    this run, T = {T*1e3:.0f} ms: J {Jz/1e3:.2f} kN s/m, zone pull {Jz/T/1e3:.0f} kN/m ({Jz/T/P12:.2f} m), peak with the bare level "
                  f"{r['pk_hd']/1e3:.0f} kN/m ({r['pk_hd']/P12:.2f} m), L_on {energized_length(30000.0, T, L_MAIN):,.0f} m, "
                  f"power {P12*GAP0*(Jz/T/P12)*30000/1e9:.2f} GW, excursion {r['max_exc']*1e3:.1f} mm, closed form J T / 8 mu = {Jz*T/(8*MU)*1e3:.1f} mm.")
        print(f"  The width in this run counts the bare sled's level added to the zone pull where they overlap; the first run's 0.66 m is the zone pull alone.")

    # ------------------------------------------------------------------ B
    rows_B = []
    if "B" in parts:
        hr("PART B: DOES THE LOAD LENGTH DROP OUT (rule %d, zones 25 to 400 m, front, middle, rear)" % rule)
        t0 = time.time()
        zones = (25.0, 50.0, 100.0, 200.0, 400.0)
        for u in SPEEDS:
            lanes = []
            for L in zones:
                T = rule_T(rule, u, L)
                for pos in pz.POSITIONS:
                    lanes.append(spread_lane(u, L, T, pos=pos, label=f"L{L:g} {pos}"))
            res, _ = pz.point_pulsed(lanes, u, tabs)
            for r in res:
                r["width"] = max(r["pk_hd"], r["pk_hu"]) / P12
            rows_B.extend(res)
        print(f"  ({time.time()-t0:.0f} s)")
        print(f"  Largest excursion over the launch (mm) and lift compensation width (m, largest over speeds above 7 km/s):")
        print(f"  {'zone m':>7} " + " ".join(f"{pos+' exc':>10} {pos+' w':>8}" for pos in pz.POSITIONS))
        worst = {}
        for L in zones:
            line = f"  {L:7.0f} "
            for pos in pz.POSITIONS:
                sel = [r for r in rows_B if r["L_sc"] == L and r["pos"] == pos]
                e = max(r["max_exc"] for r in sel); wd = max(r["width"] for r in sel if r["u"] > 7000)
                worst[(L, pos)] = (e, wd)
                line += f"{e*1e3:10.1f} {wd:8.2f} "
            print(line)
        es = [v[0] for v in worst.values()]; ws = [v[1] for v in worst.values()]
        print(f"  Largest difference between any two cases: excursion {(max(es)-min(es))*1e3:.1f} mm, width {max(ws)-min(ws):.2f} m "
              f"(excursion {min(es)*1e3:.1f} to {max(es)*1e3:.1f} mm, width {min(ws):.2f} to {max(ws):.2f} m).")
        print(f"  At 30 km/s by zone (front): " + ", ".join(f"{L:g} m: {next(r for r in rows_B if r['L_sc']==L and r['pos']=='front' and r['u']==30000)['max_exc']*1e3:.1f} mm" for L in zones))
        print(f"  At 1 km/s by zone (front): " + ", ".join(f"{L:g} m: {next(r for r in rows_B if r['L_sc']==L and r['pos']=='front' and r['u']==1000)['max_exc']*1e3:.1f} mm" for L in zones))

    # ------------------------------------------------------------------ C
    if "C" in parts:
        hr("PART C: THE BARE SLED AND THE LOW-SPEED SIDE (rule %d)" % rule)
        print(f"  The bare sled's load lasts L_sled / u ({L_SLED/1000:.0f} s at 1 km/s, {L_SLED/30000:.2f} s at 30 km/s), far longer than T, so it is")
        print(f"  not spread: the same track cancels it in step at {MATCH*100:.1f} % of its level while the bare sled is over the section,")
        print(f"  with the coil current regulated against the gap (open loop it runs away, growth time sqrt(g mu / 2F)).")
        print(f"  {'u km/s':>7} {'q_bare kN/m':>11} {'cancel kN/m':>11} {'width m':>8} {'duration s':>10} {'growth ms':>9}")
        for u in (1000.0, 3000.0, 7000.0, 10000.0, 20000.0, 30000.0):
            qb = float(pz.q_bare(u)); F = MATCH * abs(qb)
            tg = math.sqrt(GAP0 * MU / (2 * F)) if F > 0 else float('inf')
            print(f"  {u/1e3:7.0f} {qb/1e3:11.1f} {F/1e3:11.1f} {F/P12:8.3f} {L_SLED/u:10.2f} {tg*1e3:9.0f}")
        print(f"\n  Load compensation track below the neutral speed, rule {rule} (from Part A):")
        print(f"  {'u km/s':>7} {'T ms':>7} {'F kN/m':>8} {'width m':>8} {'exc mm':>7}")
        for r in rows_A[rule]:
            if r["u"] < 7000:
                print(f"  {r['u']/1e3:7.0f} {r['T']*1e3:7.1f} {r['F']/1e3:8.0f} {r['width']:8.2f} {r['max_exc']*1e3:7.1f}")
        print(f"\n  The 12.7 mm at 1 km/s of the first run came from the 50 m section steps of a matched 100 ms pulse.  At 1 km/s, zone 100 m:")
        lanes = [pz.lane(L_sc=L_MAIN, tau=TAU0, L_section=Ls, gap_dep=False, steady_loop_hz=1.0, label=f"matched, {Ls:g} m sections") for Ls in (10.0, 25.0, 50.0, 100.0)]
        lanes += [spread_lane(1000.0, L_MAIN, rule_T(rule, 1000.0, L_MAIN), label=f"spread rule {rule}")]
        lanes += [spread_lane(1000.0, L_MAIN, rule_T(rule, 1000.0, L_MAIN), L_section=10.0, label=f"spread rule {rule}, 10 m sections")]
        res, _ = pz.point_pulsed(lanes, 1000.0, tabs)
        for r in res:
            print(f"    {r['label']:34} excursion {r['max_exc']*1e3:6.1f} mm, force {max(r['pk_hd'], r['pk_hu'])/1e3:5.0f} kN/m")

    # ------------------------------------------------------------------ D
    rows_D = {}
    if "D" in parts:
        hr("PART D: TOLERANCES WITH THE SPREAD PULL (rule %d, zone 100 m, at 5, 10, 20 and 30 km/s)" % rule)
        t0 = time.time()
        tests = []
        for de in (1e-3, 2e-3, 5e-3, 10e-3, 20e-3):
            tests += [(f"timing {de*1e3:+.0f} ms", dict(dt_err=de)), (f"timing {-de*1e3:+.0f} ms", dict(dt_err=-de))]
        for e in (0.01, 0.05):
            tests += [(f"mass {e*100:+.0f} %", dict(mass_err=e)), (f"mass {-e*100:+.0f} %", dict(mass_err=-e)),
                      (f"speed {e*100:+.0f} %", dict(speed_err=e)), (f"speed {-e*100:+.0f} %", dict(speed_err=-e))]
        for tau in (1e-3, 2e-3, 5e-3, 10e-3):
            tests.append((f"rise time {tau*1e3:.0f} ms", dict(tau=tau)))
        tests += [("rise time +20 %", dict(tau_err=0.2)), ("rise time -20 %", dict(tau_err=-0.2))]
        for Ls in (10.0, 25.0, 50.0, 100.0):
            tests.append((f"sections {Ls:.0f} m, point at section end", dict(L_section=Ls, pos_in_section=1.0)))
        for ge in (1e-3, 5e-3):
            tests += [(f"gap sensor {ge*1e3:+.0f} mm", dict(gap_err=ge)), (f"gap sensor {-ge*1e3:+.0f} mm", dict(gap_err=-ge))]
        us = (5000.0, 10000.0, 20000.0, 30000.0)
        for u in us:
            lanes = [spread_lane(u, L_MAIN, rule_T(rule, u, L_MAIN), label=name, **kw) for name, kw in tests]
            res, _ = pz.point_pulsed(lanes, u, tabs)
            rows_D[u] = {r["label"]: r for r in res}
        print(f"  ({time.time()-t0:.0f} s)")
        print(f"  Excursion (mm) at each speed; '*' outside the 20 mm band or limiter above 1 kN/m:")
        print(f"  {'error':>36} " + " ".join(f"{u/1e3:>6.0f} km/s" for u in us) + f" {'worst':>7}")
        for name, _ in tests:
            vals = [rows_D[u][name] for u in us]
            line = f"  {name:>36} " + " ".join(f"{r['max_exc']*1e3:10.1f}{'*' if not r['ok_band'] else ' '}" for r in vals)
            print(line + f" {max(r['max_exc'] for r in vals)*1e3:7.1f}")
        def band(prefix, items):
            ok = [v for v, name in items if all(rows_D[u][name]["ok_band"] for u in us)]
            return f"{prefix}: " + (f"up to {max(ok):g}" if ok else "none")
        print("  Tolerance band (largest error of each kind inside 20 mm at all four speeds):")
        print("    " + band("timing error (ms)", [(de, f"timing {s*de:+.0f} ms") for de in (1, 2, 5, 10, 20) for s in (1, -1)]))
        print("    " + band("mass error (%)", [(e, f"mass {s*e:+.0f} %") for e in (1, 5) for s in (1, -1)]))
        print("    " + band("speed error (%)", [(e, f"speed {s*e:+.0f} %") for e in (1, 5) for s in (1, -1)]))
        print("    " + band("rise time (ms)", [(t, f"rise time {t:.0f} ms") for t in (1, 2, 5, 10)]))
        print("    " + band("section length (m)", [(L, f"sections {L:.0f} m, point at section end") for L in (10, 25, 50, 100)]))
        print("    " + band("gap sensor error (mm)", [(g, f"gap sensor {s*g:+.0f} mm") for g in (1, 5) for s in (1, -1)]))
        print("  Matched pulse for comparison at 30 km/s (regulated, 1 Hz loop), timing error:")
        lanes = [pz.lane(L_sc=L_MAIN, tau=TAU0, L_section=L_SECTION, gap_dep=False, steady_loop_hz=1.0, dt_err=de, label=f"{de*1e3:+.0f}")
                 for de in (0, 1e-3, -1e-3, 2e-3, -2e-3, 5e-3, -5e-3, 10e-3, -10e-3, 20e-3, -20e-3)]
        resM, _ = pz.point_pulsed(lanes, 30000.0, tabs)
        rows_D["matched30"] = {r["label"]: r for r in resM}
        print("    " + ", ".join(f"{r['label']} ms: {r['max_exc']*1e3:.1f}" for r in resM) + " mm")

    # ------------------------------------------------------------------ F (before E so the bandwidth is known)
    if "F" in parts:
        hr("PART F: STABILITY MARGIN (rule %d)" % rule)
        print(f"  Open-loop growth time of the attractive pull sqrt(g mu / 2F) against T, and the regulation bandwidth, taken as ten")
        print(f"  times the growth rate (the loop must be well above 1 / growth time).  The bare sled's level is the long pull.")
        print(f"  {'u km/s':>7} {'T ms':>7} {'F kN/m':>8} {'growth ms':>9} {'T/growth':>8} {'bare kN/m':>9} {'bare growth ms':>14} {'bare T s':>8}")
        bw = 0.0
        for r in rows_A[rule]:
            u = r["u"]; F = r["F"]
            tg = math.sqrt(GAP0 * MU / (2 * F)) if F > 0 else float('inf')
            Fb = MATCH * abs(float(pz.q_bare(u))); tgb = math.sqrt(GAP0 * MU / (2 * Fb)) if Fb > 0 else float('inf')
            bw = max(bw, 10.0 / tg, 10.0 / tgb)
            if u in (1000, 2000, 3000, 5000, 7000, 8000, 10000, 15000, 20000, 25000, 30000):
                print(f"  {u/1e3:7.0f} {r['T']*1e3:7.1f} {F/1e3:8.0f} {tg*1e3:9.1f} {r['T']/tg:8.2f} {Fb/1e3:9.1f} {tgb*1e3:14.0f} {L_SLED/u:8.2f}")
        print(f"  Regulation bandwidth for the chapter: {bw:.0f} Hz (ten times the fastest growth rate over the launch).")
        rows_A["bw"] = bw

    # ------------------------------------------------------------------ I (losses, before E and G so the budget can use them)
    losses = None
    t_p, v_p, P_p = pz.lsm_profile()
    import lsm_simulation as ls
    if "I" in parts:
        hr("PART I: HYSTERESIS AND EDDY CURRENT LOSSES (rule %d, both compensation tracks)" % rule)
        mu_r, rho_p, d_p = lc.FERRO_MU_R, lc.FERRO_RESISTIVITY, lc.FERRO_PLATE_THICKNESS
        wh0, rho_m, cp = lc.FERRO_HYSTERESIS_COEFF, lc.FERRO_DENSITY, 450.0
        STEINMETZ = 1.6
        w_lift = max(r["width"] for r in rows_A[rule] if r["u"] > 7000)
        w_load = max(r["width"] for r in rows_A[rule] if r["u"] < 7000)
        print(f"  Plate on the cable (levitation_config.py, sintered iron in resin): mu_r {mu_r}, resistivity {rho_p:.0e} ohm m, {d_p*1e3:.0f} mm thick,")
        print(f"  hysteresis {wh0} J/m^3 per cycle at 1.2 T scaled as (B/1.2)^{STEINMETZ} (the coefficient and exponent are design guesses, as is")
        print(f"  the resistivity), density {rho_m} kg/m^3, specific heat {cp:.0f} J/kg K assumed.  The plate is taken as the track width wide and")
        print(f"  magnetized through its full thickness.  Track widths fixed over the launch: lift compensation {w_lift:.2f} m, load compensation {w_load:.2f} m,")
        print(f"  so the plate field at each speed is B = sqrt(2 mu0 F / width), 1.2 T only at the sizing speed.")
        print(f"  Each piece of plate crosses the energized stretch once per pass at u + {pz.V_REL:,.0f} m/s: the field rises over one gap length")
        print(f"  (a few microseconds), holds for L_on / (u + V_rel), and falls.  Eddy loss from a 1-D magnetic diffusion of that cycle into the")
        print(f"  plate (field applied on the pole face, no flux leaving the far face; diffusivity 1/(mu0 mu_r sigma) = {1/(lc.MU0*mu_r/rho_p):.3f} m^2/s,")
        print(f"  diffusion time through the plate {lc.MU0*mu_r/rho_p*d_p**2/math.pi**2*1e3:.1f} ms).  The gap regulation holds the force, so the plate field is")
        print(f"  steady during the hold and the regulation ripple adds nothing to the plate; it rides on the coil current (below).")

        def eddy_cycle(B, t_edge, t_hold):
            """Eddy energy per unit plate area for one rise-hold-fall cycle, J/m^2, by explicit 1-D diffusion."""
            nz = 100
            dz = d_p / nz
            D = 1.0 / (lc.MU0 * mu_r / rho_p)
            dt = 0.4 * dz * dz / D
            t_end = 2 * t_edge + t_hold + 15e-3
            n = int(t_end / dt) + 1
            Bz = np.zeros(nz + 1)
            E = 0.0
            for j in range(n):
                t = j * dt
                if t < t_edge:
                    Bs = B * t / t_edge
                elif t < t_edge + t_hold:
                    Bs = B
                elif t < 2 * t_edge + t_hold:
                    Bs = B * (1 - (t - t_edge - t_hold) / t_edge)
                else:
                    Bs = 0.0
                Bz[0] = Bs
                Bz[-1] = Bz[-2]
                lap = (Bz[2:] - 2 * Bz[1:-1] + Bz[:-2]) / dz**2
                Jc = (Bz[1:] - Bz[:-1]) / dz / (lc.MU0 * mu_r)      # current density between cells
                E += rho_p * np.sum(Jc * Jc) * dz * dt
                Bz[1:-1] += D * lap * dt
            return E

        print(f"\n  {'u km/s':>7} {'track':>5} {'B T':>5} {'edge us':>7} {'hold ms':>7} {'hyst J/m':>8} {'eddy J/m':>8} {'sum J/m':>8} {'dT mK':>6} "
              f"{'per lap MJ':>10} {'drag N':>7} {'LIM kW':>7}")
        per_speed = {}
        for r in rows_A[rule]:
            u = r["u"]
            lift = u > 7000
            w = w_lift if lift else w_load
            F = r["F"]
            if F <= 0:
                per_speed[u] = dict(E=0.0, hyst=0.0, eddy=0.0, B=0.0, w=w)
                continue
            B = math.sqrt(2 * lc.MU0 * F / w)
            v_sweep = u + pz.V_REL
            t_edge = GAP0 / v_sweep
            t_hold = r["L_on"] / v_sweep
            hyst = wh0 * (B / 1.2) ** STEINMETZ * d_p * w
            eddy = eddy_cycle(B, t_edge, t_hold) * w
            E = hyst + eddy
            per_lap = E * ring.L_RING * (1 + pz.V_REL / u)
            drag = E * v_sweep / pz.V_REL
            per_speed[u] = dict(E=E, hyst=hyst, eddy=eddy, B=B, w=w, per_lap=per_lap, drag=drag)
            if u in (1000, 2000, 3000, 5000, 8000, 10000, 15000, 20000, 25000, 30000):
                print(f"  {u/1e3:7.0f} {('lift' if lift else 'load'):>5} {B:5.2f} {t_edge*1e6:7.2f} {t_hold*1e3:7.1f} {hyst:8.1f} {eddy:8.1f} {E:8.1f} "
                      f"{E/(rho_m*d_p*w*cp)*1e3:6.2f} {per_lap/1e6:10.1f} {drag:7.1f} {drag*pz.V_REL/1e3:7.1f}")
        print(f"  dT is the bulk plate temperature rise per pass; drag is the force on the cable that the plate losses amount to, energy")
        print(f"  lost per meter of cable travel = E (u + V_rel) / V_rel; LIM is the power to replace it, E (u + V_rel).")
        # per launch from the profile
        dt_p = np.diff(t_p)
        E_launch_plate = 0.0; H_launch = 0.0; EY_launch = 0.0
        split = {"load": 0.0, "lift": 0.0}
        for i, dtt in enumerate(dt_p):
            vv = v_p[i + 1]
            ps_ = min(per_speed.items(), key=lambda kv: abs(kv[0] - vv))[1]
            E_launch_plate += ps_["E"] * (vv + pz.V_REL) * dtt
            H_launch += ps_["hyst"] * (vv + pz.V_REL) * dtt
            EY_launch += ps_["eddy"] * (vv + pz.V_REL) * dtt
            split["lift" if vv > 7000 else "load"] += ps_["E"] * (vv + pz.V_REL) * dtt
        KE = 0.5 * (M_SC0 + M_SLED_M * L_SLED) * 30000**2
        print(f"  Per launch (profile of lsm_simulation, {t_p[-1]/60:.0f} min): hysteresis {H_launch/1e9:.2f} GJ, eddy {EY_launch/1e9:.2f} GJ, "
              f"total {E_launch_plate/1e9:.2f} GJ in the plates, {E_launch_plate/KE*100:.4f} % of the sled's kinetic energy;")
        print(f"  of that, the load compensation plate (below 7 km/s) {split['load']/1e9:.2f} GJ ({split['load']/KE*100:.4f} %), "
              f"the lift compensation plate {split['lift']/1e9:.2f} GJ ({split['lift']/KE*100:.4f} %).")
        # casing side: coil conductor and core
        j_cu = 50e6; rho_cu = 1.68e-8; pitch = 0.5
        NI_peak = 1.2 * 2 * GAP0 / lc.MU0
        print(f"\n  Casing side.  Conductor: copper at {j_cu/1e6:.0f} A/mm^2 (one 50 ms pulse per lap at any point, duty under 1e-4, so no case for a")
        print(f"  cryogenic driven coil; adiabatic heating per pulse {rho_cu*j_cu**2*0.05/(8960*385):.2f} K), ampere-turns per pole pair {NI_peak:,.0f} at 1.2 T, turn length")
        print(f"  2 x (track width + {pitch} m pole pitch) [UNVERIFIED winding layout].  Core: laminated steel yoke, 20 J/m^3 per cycle at 1.2 T scaled as")
        print(f"  (B/1.2)^{STEINMETZ}, 0.3 m deep [UNVERIFIED], eddy loss in laminations neglected.")
        E_cu_launch = 0.0; E_core_launch = 0.0
        print(f"  {'u km/s':>7} {'track':>5} {'copper kg/m':>11} {'coil loss J/pulse/m':>19} {'core J/pulse/m':>14}")
        coil = {}
        for r in rows_A[rule]:
            u = r["u"]; lift = u > 7000; w = w_lift if lift else w_load
            B = per_speed[u]["B"]
            V_cu = (NI_peak / j_cu) * 2 * (w + pitch)
            T_on = r["T"] + L_SECTION / u
            loss = rho_cu * (j_cu * B / 1.2) ** 2 * V_cu * T_on
            core = 20.0 * (B / 1.2) ** STEINMETZ * w * 0.3
            coil[u] = dict(loss=loss, core=core)
            if u in (1000, 3000, 5000, 10000, 20000, 30000):
                print(f"  {u/1e3:7.0f} {('lift' if lift else 'load'):>5} {V_cu*8960:11.0f} {loss:19.1f} {core:14.1f}")
        cu_split = {"load": 0.0, "lift": 0.0}
        for i, dtt in enumerate(dt_p):
            vv = v_p[i + 1]
            c = min(coil.items(), key=lambda kv: abs(kv[0] - vv))[1]
            E_cu_launch += (c["loss"] + c["core"]) * vv * dtt      # every meter of ring passed is pulsed once
            cu_split["lift" if vv > 7000 else "load"] += (c["loss"] + c["core"]) * vv * dtt
        print(f"  Per launch: coil and core {E_cu_launch/1e9:.2f} GJ, of which the load compensation track {cu_split['load']/1e9:.2f} GJ "
              f"and the lift compensation track {cu_split['lift']/1e9:.2f} GJ.")
        print(f"  The copper loss per pulse is rho j^2 V T_on with V = NI / j x turn length, so it scales with j; halving the current density halves")
        print(f"  the loss and doubles the copper.  Same winding at other current densities (lift compensation track, per meter of ring):")
        for jj in (50e6, 20e6, 10e6, 5e6):
            print(f"    {jj/1e6:3.0f} A/mm^2: copper {(NI_peak/jj)*2*(w_lift+pitch)*8960:5.0f} kg/m, coil loss per launch {E_cu_launch*jj/j_cu/1e9:6.0f} GJ, "
                  f"heating per pulse {rho_cu*jj**2*0.05/(8960*385):.3f} K")
        losses = dict(plate=E_launch_plate, hyst=H_launch, eddy=EY_launch, coil=E_cu_launch, per_speed=per_speed, w_lift=w_lift, w_load=w_load)

    # ------------------------------------------------------------------ E
    rows_E = []
    if "E" in parts:
        hr("PART E: FAILURE CASES WITH THE SPREAD PULL (beam model, rule %d, zone 100 m, 10 m truss and no stiffness)" % rule)
        trusses = [sp.truss(10.0, 0.1), sp.ZERO_STIFFNESS]
        xc = 5000.0
        for u in (10000.0, 20000.0, 30000.0):
            T = rule_T(rule, u, L_MAIN)
            L_on = energized_length(u, T, L_MAIN)
            cases = [dict(label="normal"),
                     dict(label="one 50 m section out", fail_x=(xc, xc + 50)),
                     dict(label="200 m out", fail_x=(xc - 100, xc + 100)),
                     dict(label="1 km supply out", fail_x=(xc - 500, xc + 500)),
                     dict(label="whole pull missing at one place", fail_x=(xc - L_on / 2, xc + L_on / 2)),
                     dict(label="pull 10 ms late", late=10e-3), dict(label="pull 10 ms early", late=-10e-3),
                     dict(label="pull 50 ms late", late=50e-3), dict(label="pull 50 ms early", late=-50e-3),
                     dict(label="pull 100 ms late", late=100e-3), dict(label="pull 100 ms early", late=-100e-3),
                     dict(label="spurious pull, no sled", spurious=True),
                     dict(label="one section open loop", open_x=(xc, xc + 50)),
                     dict(label="all sections open loop", regulated=False),
                     dict(label="casing support loop lost", no_loop=True),
                     dict(label="no pull at all", no_pulse=True)]
            t0 = time.time()
            res, _ = pz.beam_pulsed(u, L_MAIN, trusses, cases, tabs, L_section=L_SECTION, tau=TAU0, L_win=10000.0, T_spread=T)
            rows_E.extend([dict(r, u=u) for r in res])
            print(f"\n  u = {u/1e3:.0f} km/s, T = {T*1e3:.0f} ms, energized {L_on:.0f} m, track {res[0]['width']:.2f} m ({time.time()-t0:.0f} s):")
            print(f"  {'case':>32} {'truss':>22} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'v_pk':>6} {'F_lim kN/m':>10} {'contact':>7}")
            for r in res:
                print(f"  {r['case']:>32} {r['truss']:>22} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['v_pk']:6.3f} "
                      f"{r['pk_lim']/1e3:10.1f} {('yes' if r['contact'] else 'no'):>7}")
        print(f"\n  Survived without contact (10 m truss) at every speed: " + ", ".join(sorted({r['case'].split(',')[0] for r in rows_E if r['truss'].startswith('h=10') and not r['contact']} -
              {r['case'].split(',')[0] for r in rows_E if r['truss'].startswith('h=10') and r['contact']})))
        print(f"  Contact (10 m truss) in at least one speed: " + ", ".join(sorted({r['case'].split(',')[0] for r in rows_E if r['truss'].startswith('h=10') and r['contact']})))
        Jsp = J_zone(30000.0, L_MAIN) + MATCH * abs(float(pz.q_bare(30000.0))) * (L_SLED - L_MAIN) / 30000.0
        print(f"  A spurious spread pull at 30 km/s delivers {Jsp/1e3:.1f} kN s/m (zone share plus the bare level over 10 km), {Jsp/MU:.2f} m/s on the casing;")
        print(f"  the limiters catch {math.sqrt(2*tabs['F0']/(2*tabs['k'])/MU):.2f} m/s, the casing support track's spare 123 kN/m over 100 mm another "
              f"{math.sqrt(2*123e3*GAP0/MU):.2f} m/s, so it is {'still fatal' if Jsp/MU > math.sqrt(2*(tabs['F0']/(2*tabs['k'])+123e3*GAP0)/MU) else 'survivable'} on the point model; the beam row above is the run.")

    # ------------------------------------------------------------------ G
    if "G" in parts:
        hr("PART G: POWER, ENERGY AND PLATES (rule %d, 1.2 T, fringe factor %.1f, core energy neglected, recovery %.0f %%, coil loss %.0f %%)"
           % (rule, pz.FRINGE_FACTOR, pz.RECOVERY_EFF * 100, pz.COIL_LOSS_FRACTION * 100))
        NI = 1.2 * 2 * GAP0 / lc.MU0
        laps = ls.data_x[-1] / ring.L_RING
        print(f"  Ampere-turns per pole pair for 1.2 T across two 100 mm gaps: {NI:,.0f}.  Launch: {t_p[-1]/60:.0f} min, {laps:.1f} laps, LSM peak {P_p.max()/1e9:.0f} GW.")
        print(f"  {'u km/s':>7} {'width m':>8} {'E J/m':>8} {'L_on m':>7} {'E_on MJ':>8} {'P GW':>6} {'P/LSM':>7} {'E_sec MJ':>8} {'E_net/lap GJ':>12}")
        E_net_launch = 0.0
        for r in rows_A[rule]:
            u = r["u"]; wd = r["width"]
            E_m = P12 * GAP0 * wd * pz.FRINGE_FACTOR
            P = E_m * u
            E_net = E_m * ring.L_RING * ((1 - pz.RECOVERY_EFF) + pz.COIL_LOSS_FRACTION)
            if u in (1000, 3000, 5000, 8000, 10000, 15000, 20000, 25000, 30000):
                print(f"  {u/1e3:7.0f} {wd:8.2f} {E_m:8,.0f} {r['L_on']:7.0f} {E_m*r['L_on']/1e6:8.2f} {P/1e9:6.2f} {P/P_p.max():7.3f} {E_m*L_SECTION/1e6:8.2f} {E_net/1e9:12.1f}")
        # energy per launch: integrate E_net per lap over the laps, with the width at each speed
        dx = np.diff(ls.data_x); v_mid = np.array(ls.data_v)[1:]
        E_launch = 0.0
        for dxx, vv in zip(dx, v_mid):
            rr = min(rows_A[rule], key=lambda r: abs(r["u"] - vv))
            E_launch += P12 * GAP0 * rr["width"] * ((1 - pz.RECOVERY_EFF) + pz.COIL_LOSS_FRACTION) * dxx
        print(f"  Net electrical energy per launch (width at each speed, every meter passed energized once): {E_launch/1e9:,.0f} GJ, "
              f"{E_launch/(0.5*(M_SC0 + M_SLED_M*L_SLED)*30000**2)*100:.2f} % of the sled's kinetic energy.")
        if losses:
            E_field_launch = 0.0
            for dxx, vv in zip(dx, v_mid):
                rr = min(rows_A[rule], key=lambda r: abs(r["u"] - vv))
                E_field_launch += P12 * GAP0 * rr["width"] * dxx
            tot = losses["plate"] + losses["coil"]
            print(f"  Losses (Part I) as their own line: plates on the cable {losses['plate']/1e9:.2f} GJ per launch (hysteresis {losses['hyst']/1e9:.2f}, "
                  f"eddy {losses['eddy']/1e9:.2f}), coils and cores on the casing {losses['coil']/1e9:.2f} GJ; together {tot/1e9:.2f} GJ,")
            print(f"  {tot/E_field_launch*100:.2f} % of the gap field energy cycled per launch ({E_field_launch/1e9:,.0f} GJ) and {tot/E_launch*100:.1f} % of the net electrical energy")
            print(f"  computed with the 2 % coil-loss assumption above.  With the Part I losses in place of that assumption the net electrical energy")
            print(f"  per launch is the unrecovered field energy {(1-pz.RECOVERY_EFF)*E_field_launch/1e9:,.0f} GJ plus the losses {tot/1e9:,.0f} GJ = "
                  f"{((1-pz.RECOVERY_EFF)*E_field_launch + tot)/1e9:,.0f} GJ, {((1-pz.RECOVERY_EFF)*E_field_launch + tot)/(0.5*(M_SC0 + M_SLED_M*L_SLED)*30000**2)*100:.2f} % "
                  f"of the sled's kinetic energy; at 10 A/mm^2 the copper loss falls to {losses['coil']/5/1e9:,.0f} GJ.")
        w_hd, w_hu = rows_A["summary"]["w_hd"] if "summary" in rows_A else (0, 0), 0
        s = max(r["width"] for r in rows_A[rule] if r["u"] > 7000); s_hu = max(r["width"] for r in rows_A[rule] if r["u"] < 7000)
        print(f"\n  Width budget.  Cable top face {ring.CABLE_SIDE:.2f} m: lift compensation track {s:.2f} m plus a plate 0.4 m wider "
              f"({s+0.4:.2f} m), as the casing support plates are.")
        used = 2 * lc.EML_PLATE_WIDTH + 2 * lc.SAFETY_PLATE_WIDTH
        print(f"  Cable underside {ring.CABLE_SIDE:.2f} m: two casing support plates {2*lc.EML_PLATE_WIDTH:.1f} m, two limiter plates {2*lc.SAFETY_PLATE_WIDTH:.1f} m, "
              f"load compensation plate {s_hu+0.4:.2f} m ({s_hu:.2f} m of track): {used + s_hu + 0.4:.2f} m used, {ring.CABLE_SIDE - used - s_hu - 0.4:.2f} m left.")

    # ------------------------------------------------------------------ H
    hists = {}
    profiles = {}
    if "H" in parts:
        hr("PART H: THE FULL BEAM MODEL (rule %d, zone 100 m, 10 m truss and no stiffness, cable as a string)" % rule)
        trusses = [sp.truss(10.0, 0.1), sp.ZERO_STIFFNESS]
        print(f"  {'u km/s':>7} {'T ms':>7} {'truss':>22} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'v_pk':>6} {'F_lim':>6} {'point model exc':>15}")
        for u in (5000.0, 10000.0, 20000.0, 30000.0):
            T = rule_T(rule, u, L_MAIN)
            res, prof = pz.beam_pulsed(u, L_MAIN, trusses, [dict(label="normal")], tabs, L_section=L_SECTION, tau=TAU0, L_win=10000.0, T_spread=T, snap=True)
            profiles[u] = prof
            rp = next(r for r in rows_A[rule] if r["u"] == u)
            for r in res:
                print(f"  {u/1e3:7.0f} {T*1e3:7.1f} {r['truss']:>22} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['v_pk']:6.3f} "
                      f"{r['pk_lim']/1e3:6.2f} {rp['max_exc']*1e3:15.1f}")
            if prof is not None:
                rel = prof["rel"][0]
                print(f"    along the window when the zone is at its middle: casing offset from {rel.min()*1e3:+.1f} to {rel.max()*1e3:+.1f} mm, "
                      f"pulled region {u*T:.0f} m long")
            res2, h = pz.point_pulsed([spread_lane(u, L_MAIN, T), spread_lane(u, L_MAIN, T, cancel=False)], u, tabs)
            hists[u] = h

    # ------------------------------------------------------------------ graphs
    if graphs:
        hr("GRAPHS")
        fig, axes = plt.subplots(1, 2, figsize=(sp.GRAPH_W * 1.3, sp.GRAPH_H))
        for k, col in zip((1, 2, 3, 4), PALETTE):
            rr = rows_A[k]
            axes[0].plot([r["u"] / 1e3 for r in rr], [r["F"] / 1e3 for r in rr], color=col, marker="o", ms=3, lw=1.4, label=f"rule {k}: {RULE_NAMES[k]}")
            axes[1].plot([r["u"] / 1e3 for r in rr], [r["width"] for r in rr], color=col, marker="o", ms=3, lw=1.4, label=f"rule {k}")
        for ax in axes:
            ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.grid(alpha=0.3, which="both"); ax.set_yscale("log")
        axes[0].set_ylabel("Peak pulsed track force (kN/m)"); axes[1].set_ylabel("Track width at 1.2 T (m)")
        axes[0].legend(fontsize=7); axes[1].legend(fontsize=7)
        sp._save(fig, "fig_sled_19_spread_rules.png", out_dir)
        if rows_B:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            for L, col in zip((25.0, 50.0, 100.0, 200.0, 400.0), PALETTE):
                for pos, ls in (("front", "-"), ("middle", "--"), ("rear", ":")):
                    sel = sorted([r for r in rows_B if r["L_sc"] == L and r["pos"] == pos], key=lambda r: r["u"])
                    ax.plot([r["u"] / 1e3 for r in sel], [r["max_exc"] * 1e3 for r in sel], color=col, ls=ls, lw=1.2, label=f"{L:g} m, {pos}")
            ax.axhline(BAND * 1e3, color="gray", ls="--", lw=1); ax.axhline(10, color="gray", ls=":", lw=1)
            ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel(f"Largest casing excursion from the 100 mm gap, rule {rule} (mm)")
            ax.set_ylim(0, 25); ax.grid(alpha=0.3); ax.legend(fontsize=6.5, ncol=3)
            sp._save(fig, "fig_sled_20_spread_excursion.png", out_dir)
        if rows_D:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            des = [-20, -10, -5, -2, -1, 1, 2, 5, 10, 20]
            for u, col in zip((5000.0, 10000.0, 20000.0, 30000.0), PALETTE):
                ax.plot(des, [rows_D[u][f"timing {d:+.0f} ms"]["max_exc"] * 1e3 for d in des], color=col, marker="o", ms=3, lw=1.3, label=f"spread pull, {u/1e3:g} km/s")
            dm = [-20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20]
            ax.plot(dm, [rows_D["matched30"][f"{d:+.0f}"]["max_exc"] * 1e3 for d in dm], color="black", marker="s", ms=3, lw=1.3, ls="--", label="matched pulse, 30 km/s")
            ax.axhline(BAND * 1e3, color="gray", ls="--", lw=1)
            ax.set_xlabel("Timing error of the pull (ms, positive = late)"); ax.set_ylabel("Largest casing excursion (mm)")
            ax.set_ylim(0, 105); ax.grid(alpha=0.3); ax.legend(fontsize=7)
            sp._save(fig, "fig_sled_21_spread_timing.png", out_dir)
        if hists:
            fig, axes = plt.subplots(2, 1, figsize=(sp.GRAPH_W, sp.GRAPH_H * 1.4))
            for (u, h), col in zip(sorted(hists.items()), PALETTE):
                tm = h["t"] * 1e3
                axes[0].plot(tm, (h["w"][0] - h["wc"][0]) * 1e3, color=col, lw=1.3, label=f"{u/1e3:g} km/s, spread pull")
                axes[0].plot(tm, (h["w"][1] - h["wc"][1]) * 1e3, color=col, lw=1.0, ls=":", label=f"{u/1e3:g} km/s, no pull")
                axes[1].plot(tm, h["F"][0] / 1e3, color=col, lw=1.3)
            axes[0].axhspan(-20, 20, color="gray", alpha=0.12); axes[0].set_ylim(-105, 105)
            axes[0].set_ylabel("Casing offset from the cable (mm, up positive)"); axes[0].legend(fontsize=7, ncol=2)
            axes[1].set_ylabel("Pulsed track force on the casing (kN/m, up positive)")
            axes[1].set_xlabel(f"Time after the sled front passes the point (ms), rule {rule}, zone 100 m")
            for a_ in axes:
                a_.set_xlim(-400, 1200); a_.grid(alpha=0.3)
            sp._save(fig, "fig_sled_22_spread_time_history.png", out_dir)
        if rows_E:
            cases = []
            for r in rows_E:
                c = r["case"]
                if c not in cases: cases.append(c)
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W * 1.2, sp.GRAPH_H))
            n = len(cases); width = 0.25
            for k, (u, col) in enumerate(zip((10000.0, 20000.0, 30000.0), PALETTE)):
                vals = [min(next(r for r in rows_E if r["u"] == u and r["case"] == c and r["truss"].startswith("h=10"))[g] for g in ("min_gap_b", "min_gap_t")) * 1e3 for c in cases]
                ax.bar(np.arange(n) + (k - 1) * width, vals, width, color=col, label=f"{u/1e3:g} km/s, 10 m truss")
            ax.axhline(BAND * 1e3, color="gray", ls="--", lw=1); ax.axhline(80, color="gray", ls=":", lw=1)
            ax.set_xticks(np.arange(n)); ax.set_xticklabels([c.replace(", ", "\n") for c in cases], fontsize=6.5, rotation=60, ha="right")
            ax.set_ylabel("Smallest gap, either side (mm)"); ax.set_ylim(0, 105); ax.grid(alpha=0.3, axis="y"); ax.legend(fontsize=7)
            sp._save(fig, "fig_sled_23_spread_failures.png", out_dir)
        if profiles:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            for (u, prof), col in zip(sorted(profiles.items()), PALETTE):
                if prof is None: continue
                ax.plot(prof["x"], prof["rel"][0] * 1e3, color=col, lw=1.3, label=f"{u/1e3:g} km/s, 10 m truss, pulled region {u*rule_T(rule, u, L_MAIN):.0f} m")
            ax.axhspan(-20, 20, color="gray", alpha=0.12)
            ax.set_xlabel("Distance along the ring from the spacecraft zone (m)"); ax.set_ylabel("Casing offset from the cable (mm, up positive)")
            ax.grid(alpha=0.3); ax.legend(fontsize=7)
            sp._save(fig, "fig_sled_24_spread_profile.png", out_dir)

    hr(f"DONE in {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
