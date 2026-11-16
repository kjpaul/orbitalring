#!/usr/bin/env python3
"""
A local standing load on the staggered support coils (round 9 addendum, A1, 2026-11-16: Vol_III/Claude_Code_instructions_round9_addendum_2026-11-16.md).
Book: Orbital Ring Engineering, Volume III: Chapter 7 (Section 7.7), Table 19.2.

A parked 5,000 tonne spacecraft on 100 m of sled presses 396 kN/m on the casing there, 198 kN/m on each support track on top of the track's
47.5 kN/m.  The support coils of a loop are L_c long, ten to a bundle, each starting a tenth of L_c after the one before, so the ampere-turns at a
place are the sum of the ten coils that cover it.  This script asks what coil currents best carry the load and the casing's weight everywhere else.

The model (a table and a paragraph are wanted, not a design):
  - one support track, both its loops driven alike; cells a tenth of L_c long; the load on whole cells; a window of 6 L_c (at least 3 km for the long coils, 1.2 km for coils under 500 m), periodic;
  - the pull of a cell from its ampere-turns per leg by the solved pull curve of the split support at 2.5 mm cells (round 8, bearing_gap_trade_parts
    .support_current's points), so the saturation of the plate is in it; no end effects along the track [the field of a coil's end is not modeled];
  - the coil currents by bounded least squares on the leftover (pull less load), each coil between 0 and the cap's coil current;
  - the gap's movement under the leftover: the casing a Timoshenko beam (sled_passage_simulation.truss(10 m, 0.1), the beam model's truss), the cable
    a string at its tension, the two support tracks' pull following the gap at constant flux (2 x bearing_round4.K_FLUX, scaled with the local pull),
    and the fast tracks holding the gap with the protection's law (5 Hz on the casing's mass, compensation_tracks.prot) [UNVERIFIED as the law of
    the standing regulation; it is the only gap law in the code].  Solved mode by mode on the window.
Run (from Python/orbitalring): python bearing_round9_local.py [--pass]      writes runs/round9_logs/round9_local.json and prints the tables
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import least_squares

import ring_altitude as ring
import bearing_round4 as r4
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_round8 as r8
import bearing_round9 as r9
import sled_passage_simulation as sp

W_T = r7.W_HALF                      # one support track's share of the casing's weight, N/m
F_LOAD = 198e3                       # the standing load on one support track, N/m
L_LOAD = 100.0                       # m
HZ = 5.0


def curve():
    rows = r8.load("round8_rows.json"); pts = rows["support_points"]; k = rows["k_lin"]
    NIs = np.array([p["NI"] for p in pts]); Fs = np.array([p["F"] for p in pts]); d_n = pts[0]["dev"]; NI_n = NIs[0]

    def F_of(NI):
        NI = np.maximum(NI, 0.0)
        low = k * NI**2 * (1 - (1 - d_n) * NI / NI_n)
        hi = np.exp(np.interp(np.log(np.maximum(NI, 1.0)), np.log(NIs), np.log(Fs)))
        return np.where(NI < NI_n, low, hi)
    return F_of, rows["n_t"], pts[0]["I_coil"], rows["parked"]["I_cap"], NIs, Fs


def solve(Lc, F_of, n_t, I_n, I_cap, verbose=False):
    s = Lc / 10.0; N = int(round(max(6 * Lc, 3000.0 if Lc >= 500 else 1200.0) / s)); m = max(1, int(round(L_LOAD / s))); c0 = N // 2 - m // 2
    w = np.full(N, W_T); w[c0:c0 + m] += F_LOAD * min(1.0, L_LOAD / (m * s))
    A = np.zeros((N, N))
    for j in range(N):
        for q in range(10):
            A[(j + q) % N, j] = n_t                         # coil j covers cells j .. j+9
    def jac(I):
        NI = A @ I; d = 50.0
        return ((F_of(NI + d) - F_of(NI - d)) / (2 * d) / 1e3)[:, None] * A
    res = least_squares(lambda I: (F_of(A @ I) - w) / 1e3, np.full(N, I_n), jac=jac, bounds=(0.0, I_cap), xtol=1e-10, ftol=1e-10, max_nfev=300)
    I = res.x; F = F_of(A @ I); r = F - w
    x = (np.arange(N) - c0 - m / 2 + 0.5) * s             # cell centers from the load's center
    return dict(Lc=Lc, s=s, N=N, m=m, c0=c0, x=x, I=I, F=F, w=w, r=r)


def movement(sol):
    """The gap's movement under the leftover on both support tracks, mode by mode: rel = f C / (1 + Ks C), C = 1/K_beam + 1/(T k^2)."""
    tr = sp.truss(10.0, 0.1); EI, GA = tr["EI"], tr["GA"]; T = ring.CABLE_TENSION
    N, s = sol["N"], sol["s"]; f = 2 * sol["r"]                                  # both tracks
    K = 2 * r4.K_FLUX * float(np.mean(sol["F"])) / W_T                           # the support's constant-flux pull against the gap, N/m per m (its mean over the window)
    k_reg = ring.M_LOAD_M * (2 * math.pi * HZ)**2
    fk = np.fft.rfft(f); kk = 2 * math.pi * np.fft.rfftfreq(N, d=s)
    C = np.full(len(kk), np.inf)
    nz = kk > 0
    Kb = 1.0 / (1.0 / (EI * kk[nz]**4) + 1.0 / (GA * kk[nz]**2))
    C[nz] = 1.0 / Kb + 1.0 / (T * kk[nz]**2)
    out = {}
    for name, Ks in (("fast tracks holding", k_reg - K), ("structure alone", -K)):
        den = 1 + Ks * C
        relk = np.where(np.isinf(C), fk / Ks if Ks != 0 else 0.0, fk * C / np.where(np.abs(den) < 1e-12, 1e-12, den))
        rel = np.fft.irfft(relk, n=N)
        out[name] = dict(rel=rel, unstable=bool((den[nz] <= 0).any()) if Ks < 0 else False, hold=k_reg * rel if Ks > 0 else None)
    return out, dict(K=K, k_reg=k_reg, EI=EI, GA=GA, T=T)


def part_pass():
    """Item 4: the heaviest sled at 30 km/s over a place that carries the standing load (the point model; the casing's mass and weight per meter there
    include the parked spacecraft, the support is at its cap and follows the gap at constant flux)."""
    import sled_passage_pulsed as pz
    import bearing_round9_parts as p9
    pz_, sp_, ct, tabs, lane, c = p9._setup()
    old = pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING
    m_load = 5e6 / L_LOAD
    rows = []
    try:
        for u in (30e3, 1e3):
            base, _ = pz.point_pulsed([lane(u, True, label="pair")], u, tabs)
            pz.MU = old[0] + m_load; pz.W = old[1] + 2 * F_LOAD; pz.F_STEADY = pz.W; pz.MATCH_RING = pz.MU_C / (pz.MU + pz.MU_C)
            res, _ = pz.point_pulsed([lane(u, True, label="pair, the standing load on the casing"), lane(u, True, match=old[3], label="the same, the planned pull as for the bare casing")], u, tabs)
            pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING = old
            rows.append(dict(u=u, bare=base[0]["max_exc_all"], bare_gap=min(base[0]["min_gap_b"], base[0]["min_gap_t"]), loaded=res[0]["max_exc_all"], loaded_gap=min(res[0]["min_gap_b"], res[0]["min_gap_t"]),
                             loaded_contact=bool(res[0]["contact"]), loaded_bareplan=res[1]["max_exc_all"], loaded_bareplan_gap=min(res[1]["min_gap_b"], res[1]["min_gap_t"]), contact2=bool(res[1]["contact"])))
    finally:
        pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING = old
    print(f"\n  ITEM 4.  A sled passing a place that carries the standing load (point model, the pair, 10 mm protection; casing {old[0]+m_load:,.0f} kg/m and {(old[1]+2*F_LOAD)/1e3:.0f} kN/m there, the support at its cap and at constant flux):")
    print(f"  {'u km/s':>7} {'bare casing: movement mm':>24} {'gap mm':>7} | {'loaded place, the pull planned for its mass: mm':>47} {'gap':>6} | {'loaded place, the pull planned as for the bare casing: mm':>56} {'gap':>6}")
    for r in rows:
        print(f"  {r['u']/1e3:7.0f} {r['bare']*1e3:24.2f} {r['bare_gap']*1e3:7.1f} | {r['loaded']*1e3:47.2f} {r['loaded_gap']*1e3:6.1f}{' CONTACT' if r['loaded_contact'] else ''} | {r['loaded_bareplan']*1e3:56.2f} {r['loaded_bareplan_gap']*1e3:6.1f}{' CONTACT' if r['contact2'] else ''}")
    return rows


def main():
    F_of, n_t, I_n, I_cap, NIs, Fs = curve()
    ld = p7.lead()
    out = dict(I_n=I_n, I_cap=I_cap, n_t=n_t)
    print("=" * 150 + "\n  A1: A LOCAL STANDING LOAD ON THE STAGGERED SUPPORT COILS\n" + "=" * 150)
    print(f"  One support track: {W_T/1e3:.1f} kN/m everywhere and {F_LOAD/1e3:.0f} kN/m more on {L_LOAD:.0f} m.  A coil has {n_t:.1f} turns of ten layers; {I_n:.0f} A normal, {I_cap:,.0f} A at the cap.")
    rows = []
    for Lc in (1000.0, 500.0, 200.0, 100.0, 50.0, 25.0):
        sol = solve(Lc, F_of, n_t, I_n, I_cap); mv, par = movement(sol)
        r = sol["r"]; i_hi = int(np.argmax(r)); i_lo = int(np.argmin(r)); sl = slice(sol["c0"], sol["c0"] + sol["m"])
        fh = mv["fast tracks holding"]; rel = fh["rel"]; hold = fh["hold"] / 2                       # per load track
        row = dict(Lc=Lc, s=sol["s"], I_max=float(sol["I"].max()), I_min=float(sol["I"].min()), F_load=float(sol["F"][sl].mean()), short=float(-r[sl].mean()), over=float(r[i_hi]), x_over=float(sol["x"][i_hi]),
                   under=float(r[i_lo]), x_under=float(sol["x"][i_lo]), net=float(r.sum() * sol["s"]), reach=float(np.abs(sol["x"][np.abs(r) > 1e3]).max()) if (np.abs(r) > 1e3).any() else 0.0,
                   rel_max=float(np.abs(rel).max()), hold_max=float(np.abs(hold).max()), alone_unstable=mv["structure alone"]["unstable"], coils_km=4 * 10 * 1000.0 / Lc)
        rows.append(row)
        if Lc == 1000.0:
            out["profile_1km"] = dict(x=sol["x"].tolist(), pull=sol["F"].tolist(), load=sol["w"].tolist(), leftover=r.tolist(), I=sol["I"].tolist(), rel=rel.tolist(), hold=hold.tolist())
            print(f"\n  ITEM 1.  The 1 km coils as designed (cells of 100 m, the load on one cell; window {sol['N']*sol['s']/1e3:.0f} km).  The pull along the track and what is left over, per support track:")
            print(f"  {'from the load center m':>22} {'load kN/m':>9} {'pull kN/m':>9} {'leftover kN/m':>13} {'gap moves mm (fast tracks holding)':>34} {'each load track holds kN/m':>26}")
            for i in range(sol["N"]):
                if abs(sol["x"][i]) <= 1250:
                    print(f"  {sol['x'][i]:22.0f} {sol['w'][i]/1e3:9.1f} {sol['F'][i]/1e3:9.1f} {r[i]/1e3:13.1f} {rel[i]*1e3:34.2f} {hold[i]/1e3:26.1f}")
            cov = [j % sol["N"] for j in range(sol["c0"] - 9, sol["c0"] + 1)]
            print(f"  The ten coils that cover the load: " + ", ".join(f"{sol['I'][j]:.0f}" for j in cov) + f" A.  The others between {np.delete(sol['I'], cov).min():.0f} and {np.delete(sol['I'], cov).max():.0f} A.")
            print(f"  (stiffnesses per meter of ring: the support's constant-flux pull {par['K']/1e6:.2f} MN/m2, destabilizing; the fast tracks' law {par['k_reg']/1e6:.1f} MN/m2; casing EI {par['EI']:.2e} N m2, GA {par['GA']:.2e} N; cable {par['T']/1e9:.0f} GN)")
    p_steady = (0.048)**2 * p7.F_FULL
    print(f"\n  ITEMS 1 TO 3.  Against the coil's length (ten coils to a bundle, each a tenth of its length after the one before; the load on 100 m):")
    print(f"  {'coil m':>6} {'largest coil A':>14} {'of the cap':>10} {'pull under the load kN/m':>24} {'short there kN/m':>16} {'largest surplus kN/m':>20} {'at m':>6} {'leftover reaches m':>18} {'gap moves mm':>12} "
          f"{'a load track holds kN/m':>23} {'coils per km (4 loops)':>22}")
    for r in rows:
        print(f"  {r['Lc']:6.0f} {r['I_max']:14.0f} {r['I_max']/I_cap:10.2f} {r['F_load']/1e3:24.1f} {r['short']/1e3:16.1f} {r['over']/1e3:20.1f} {r['x_over']:6.0f} {r['reach']:18.0f} {r['rel_max']*1e3:12.2f} {r['hold_max']/1e3:23.1f} {r['coils_km']:22.0f}")
    print(f"  (the gap's movement with the fast tracks holding it by the 5 Hz law; what a load track then holds steadily is set against what its leads carry steadily, {p_steady/1e3:.2f} kN/m:")
    print(f"  round 6a's 4.8 percent of full current, squared, of {p7.F_FULL/1e3:.0f} kN/m.  With the structure alone the long waves are unstable in every row, as they are with no load: {all(r['alone_unstable'] for r in rows)}.)")
    out["rows"] = rows; out["p_steady"] = p_steady
    if "--pass" in sys.argv:
        out["pass"] = part_pass()
    r9.save("round9_local.json", out)


if __name__ == "__main__":
    main()
