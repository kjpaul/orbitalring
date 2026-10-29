#!/usr/bin/env python3
"""
Design history: the limits study of the sled passage work that led to the adopted bearing in `final_bearing.py`.
What the three levers (coil length, truss depth, control) buy against the sled passage, and the rejected limiters.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (the limits that led to the compensation tracks).
Run: python sled_passage_limits.py [see --help]
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

Sled passage limits: what the three levers buy.

Companion to sled_passage_simulation.py.  That run showed that above about
5 km/s the casing stiffness does not change the smallest gap, because the
sled outruns every wave in the casing and every meter gets the same impulse.
So this script uses a single point of casing (12,000 kg/m, no stiffness) with
the same forces, which is exact in that regime and runs in seconds, and
checks itself against the full run at a few points.

Three limits are mapped:

  1. EML track width.  The EML pulls only, so it matters at low speed, where
     the spacecraft's weight must be carried.  Required force is closed form:
     F = casing weight + (525 + M/L_sc) * a_down(u).  Width scales the force
     at fixed current and field (Chapter 7: two 1 m tracks give 218 kN/m at
     400 A, 871 kN/m at 800 A, 1,146 kN/m at plate saturation).

  2. Bottom limiter strength and reach.  Strength is the ampere-turns per
     meter K (F0 goes as K^2); reach is the Halbach wavelength (decay length
     lambda / 2 pi; energy stored F0 lambda / 4 pi).  Both are scanned, and
     for each pair the top sled speed for a 5,000 t spacecraft is found, with
     the costs at the nominal gap: steady force, eddy heating and drag.

  3. Spacecraft mass against sled speed (and radial g) for the present
     limiter and for upgraded ones, both for the 20 mm forbidden zone and for
     no contact.

Usage:
    python sled_passage_limits.py              # all three, graphs on
    python sled_passage_limits.py --no-graphs
    python sled_passage_limits.py --l_sc=200
    python sled_passage_limits.py --help

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
import levitation_physics as lp
import sled_passage_simulation as sp

MU = sp.MU
W = sp.W_CASING
G0 = sp.G0_STD
GAP0 = sp.GAP0
GAPF = sp.GAP_FORBIDDEN
L_SLED = sp.L_SLED
M_SLED_M = sp.M_SLED_M
M_SC0 = sp.M_SC
PALETTE = sp.PALETTE

# EML scaling per meter of total track width (two tracks of 1 m give the Chapter 7 values)
EML_WIDTH0 = 2.0 * lc.EML_COIL_WIDTH            # m of track, present design
EML_PER_M_400 = sp.F_EML_OP / EML_WIDTH0        # N/m per meter of track at 400 A
EML_PER_M_800 = sp.F_EML_MAX / EML_WIDTH0       # at 800 A (stated maximum)
P_SAT = lc.FERRO_B_SAT**2 / (2 * lc.MU0)        # plate saturation pressure, Pa
EML_PER_M_SAT = P_SAT                           # N/m per meter of track at saturation

# =============================================================================
# LIMITER DESIGN
# =============================================================================

def limiter_design(k_mult=1.0, wavelength=lc.SAFETY_WAVELENGTH, material=lc.SAFETY_PLATE_MATERIAL):
    """Limiter law for an ampere-turn multiplier k_mult (sheet current K = k_mult * K0)
    and a Halbach wavelength.  Returns F0, Q0 at contact, k, the values at the nominal
    gap, the drag at the nominal gap, and the energy stored over the whole travel."""
    coils_per_m = 4.0 / wavelength
    K0 = (4.0 / lc.SAFETY_WAVELENGTH) * lc.SAFETY_TURNS_PER_COIL * lc.SAFETY_I_PERSISTENT   # A/m, present
    K = k_mult * K0
    # halbach_surface_field takes turns, current and coils per meter; give it K through turns * current
    B = lp.halbach_surface_field(1, K / coils_per_m, lc.SAFETY_COIL_WIDTH, wavelength, coils_per_meter=coils_per_m)
    mat = lp.get_material(material)
    sig = 1.0 / mat["resistivity"]
    k = 2 * math.pi / wavelength
    F0 = lp.safety_repulsive_force(B, 0.0, wavelength, sig, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, sp.V_REL, n_arrays=2)
    Q0 = lp.eddy_current_heating_halbach(B, 0.0, wavelength, sp.V_REL, mat, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, n_arrays=2)
    F_nom = lp.safety_repulsive_force(B, GAP0, wavelength, sig, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, sp.V_REL, n_arrays=2)
    Q_nom = lp.eddy_current_heating_halbach(B, GAP0, wavelength, sp.V_REL, mat, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, n_arrays=2)
    D_nom = lp.safety_drag_force(B, GAP0, wavelength, sig, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, sp.V_REL, n_arrays=2)
    U_all = F0 / (2 * k) * (1 - math.exp(-2 * k * GAP0))
    U_to_20 = F0 / (2 * k) * (math.exp(-2 * k * GAPF) - math.exp(-2 * k * GAP0))
    return dict(k_mult=k_mult, wavelength=wavelength, K=K, K0=K0, B=B, F0=F0, Q0=Q0, k=k, F_nom=F_nom, Q_nom=Q_nom,
                D_nom=D_nom, U_all=U_all, U_to_20=U_to_20, mat=mat)


# =============================================================================
# SINGLE POINT MODEL, VECTORIZED OVER (mass, speed)
# =============================================================================

def point_response(M_sc, u, L_sc, lim, F_max=sp.F_EML_MAX, pd_freq_hz=sp.PD_FREQ_HZ, zeta=1.0,
                   t_settle=1.5, steps_per_pass=12):
    """Smallest bottom gap and contact flag for one meter of casing (no stiffness)
    under the ideal controller, for every pair of spacecraft mass M_sc (array) and
    sled speed u (array).  Same forces and controller as sled_passage_simulation,
    without the beam and without the one-node smearing of the zone edges."""
    M = np.asarray(M_sc, dtype=float)[:, None]
    U = np.asarray(u, dtype=float)[None, :]
    ad = sp.a_down(U)                                 # down positive
    q_zone = (M_SLED_M + M / L_sc) * ad               # (n_m, n_u)
    q_bare = M_SLED_M * ad                            # (1, n_u)
    t_zone = L_sc / U
    t_sled = L_SLED / U
    dt = min(float(L_sc / np.max(U) / steps_per_pass), 2e-3)
    t_end = float(np.max(t_sled)) + t_settle
    n = int(math.ceil(t_end / dt)) + 1
    om = 2 * math.pi * pd_freq_hz
    k2, F0 = 2 * lim["k"], lim["F0"]
    w = np.zeros(np.broadcast(M, U).shape)
    v = np.zeros_like(w)
    gmin = np.full_like(w, GAP0)
    hit = np.zeros(w.shape, dtype=bool)
    for i in range(n):
        t = i * dt
        in_zone = (t < t_zone)
        in_bare = (t >= t_zone) & (t < t_sled)
        q = np.where(in_zone, q_zone, 0.0) + np.where(in_bare, q_bare, 0.0)
        f_known = -W - q
        F_eml = np.clip(-f_known - MU * (2 * zeta * om * v + om * om * w), 0.0, F_max)
        F_lb = F0 * np.exp(-k2 * np.maximum(GAP0 - w, 0.0))
        F_lt = F0 * np.exp(-k2 * np.maximum(GAP0 + w, 0.0))
        a = (f_known + F_eml + F_lt - F_lb) / MU
        v += a * dt
        w += v * dt
        hb = w >= GAP0
        if hb.any():
            w[hb] = GAP0
            v[hb] = np.minimum(v[hb], 0.0)
            hit |= hb
        ht = w <= -GAP0
        if ht.any():
            w[ht] = -GAP0
            v[ht] = np.maximum(v[ht], 0.0)
        np.minimum(gmin, GAP0 - w, out=gmin)
    return gmin, hit


def top_speed(lim, M_sc, L_sc, target, F_max=sp.F_EML_MAX, u_lo=8000.0, u_hi=sp.U_MAX, n=12):
    """Largest u with the smallest bottom gap >= target (bisection; the gap falls with u above the neutral speed)."""
    g, _ = point_response([M_sc], [u_hi], L_sc, lim, F_max)
    if g[0, 0] >= target:
        return u_hi
    g, _ = point_response([M_sc], [u_lo], L_sc, lim, F_max)
    if g[0, 0] < target:
        return u_lo
    for _ in range(n):
        mid = 0.5 * (u_lo + u_hi)
        g, _ = point_response([M_sc], [mid], L_sc, lim, F_max)
        if g[0, 0] >= target:
            u_lo = mid
        else:
            u_hi = mid
    return 0.5 * (u_lo + u_hi)


def mass_envelope(lim, u_arr, L_sc, masses, F_max=sp.F_EML_MAX):
    """Largest spacecraft mass at each speed with gap >= 20 mm, and with no contact."""
    gmin, hit = point_response(masses, u_arr, L_sc, lim, F_max)
    M20 = np.full(len(u_arr), np.nan)
    Mnc = np.full(len(u_arr), np.nan)
    for j in range(len(u_arr)):
        ok20 = np.where(gmin[:, j] >= GAPF)[0]
        oknc = np.where(~hit[:, j])[0]
        if ok20.size:
            M20[j] = masses[ok20.max()]
        if oknc.size:
            Mnc[j] = masses[oknc.max()]
    return M20, Mnc, gmin, hit


def eml_required(M_sc, u, L_sc):
    """EML force per meter needed to hold the gap exactly under the spacecraft zone, N/m."""
    return W + (M_SLED_M + M_sc / L_sc) * sp.a_down(u)


def eml_mass_limit(F_max, u, L_sc):
    """Largest spacecraft mass the EML can carry on L_sc at speed u (below the neutral speed)."""
    ad = sp.a_down(np.asarray(u, dtype=float))
    with np.errstate(divide='ignore', invalid='ignore'):
        m = L_sc * ((F_max - W) / ad - M_SLED_M)
    return np.where(ad > 0, m, np.inf)


# =============================================================================
# MAIN
# =============================================================================

def hr(t=""):
    print("\n" + "=" * 78)
    if t:
        print(t)
        print("=" * 78)


def main():
    L_sc = 100.0
    graphs = True
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__)
            return
        elif arg.startswith("--l_sc="):
            L_sc = float(arg.split("=")[1])
        elif arg == "--no-graphs":
            graphs = False
        else:
            print(f"Unknown option {arg}")
            return
    out_dir = ring.out_dir("graphs_sled_passage")   # under ORBITAL_RING_OUT (default runs)
    t_all = time.time()

    hr("SLED PASSAGE LIMITS")
    print(f"  Design point {ring.ALTITUDE_KM:.0f} km; casing {MU:,.0f} kg/m, weight {W/1e3:.1f} kN/m; spacecraft zone {L_sc:g} m; "
          f"sled {L_SLED:,.0f} m at {M_SLED_M:.0f} kg/m.")
    print(f"  Point model: one meter of casing, no stiffness, ideal controller (feed-forward plus {sp.PD_FREQ_HZ:.0f} Hz return),")
    print(f"  pull-only EML, limiters top and bottom, contact a dead stop.  Valid above about 5 km/s (full run, 2026-10-09).")
    lim0 = limiter_design()

    # ------------------------------------------------------------ validation
    hr("CHECK AGAINST THE FULL RUN (zero stiffness lane, 5,000 t, L_sc = 100 m, ideal 871)")
    g, h = point_response([M_SC0], [10000.0, 15000.0, 20000.0, 22000.0, 23000.0, 30000.0], 100.0, lim0)
    full = {10000: 95.0, 15000: 63.8, 20000: 10.2, 22000: None, 23000: 0.0, 30000: 0.0}
    for j, u in enumerate([10000.0, 15000.0, 20000.0, 22000.0, 23000.0, 30000.0]):
        f = full[int(u)]
        print(f"  u = {u/1e3:4.0f} km/s: point model {g[0,j]*1e3:6.1f} mm{' contact' if h[0,j] else '        '}"
              f"   full run {('%.1f mm' % f) if f is not None else 'n/a'}")
    print(f"  Top speed, present design, 5,000 t: gap >= 20 mm {top_speed(lim0, M_SC0, 100.0, GAPF)/1e3:.1f} km/s "
          f"(full run bisection 19.0), no contact {top_speed(lim0, M_SC0, 100.0, 0.0005)/1e3:.1f} km/s (full run 22.0)")

    # ------------------------------------------------------------ 1. EML width
    hr("LIMIT 1: EML TRACK WIDTH (low speed, the EML must carry the spacecraft)")
    print(f"  Per meter of total track width: {EML_PER_M_400/1e3:.1f} kN/m at 400 A, {EML_PER_M_800/1e3:.1f} kN/m at 800 A, "
          f"{EML_PER_M_SAT/1e3:.1f} kN/m at plate saturation ({lc.FERRO_B_SAT} T).  Present: {EML_WIDTH0:.0f} m of track.")
    print(f"  Chapter 6 note: at 70 K the 12 mm tape carries about 400 A at 0.5 T and 316 A at 1 T, so the 400 A column is the")
    print(f"  practical one and 800 A is the stated maximum only.")
    print(f"\n  Required EML force (kN/m) and track width at 400 A (m) for a 5,000 t spacecraft, by zone length and speed:")
    print(f"  {'u km/s':>7} " + " ".join(f"{'L='+str(int(L))+' kN/m':>12} {'width m':>8}" for L in (50.0, 100.0, 200.0)))
    for u in (0.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0):
        line = f"  {u/1e3:7.1f} "
        for L in (50.0, 100.0, 200.0):
            F = float(eml_required(M_SC0, u, L))
            line += f"{F/1e3:12.1f} {F/EML_PER_M_400:8.2f} "
        print(line)
    print(f"\n  Largest spacecraft mass (t) the EML can carry on L_sc = {L_sc:g} m, by total track width at 400 A, below the neutral speed:")
    widths = [2.0, 4.0, 6.0, 8.0, 12.0]
    print(f"  {'u km/s':>7} " + " ".join(f"{'%g m' % wd:>10}" for wd in widths) + f"   {'8 m at 800 A':>13}")
    for u in (0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0):
        line = f"  {u/1e3:7.1f} "
        for wd in widths:
            m = float(eml_mass_limit(wd * EML_PER_M_400, u, L_sc))
            line += f"{m/1e3:10,.0f} "
        line += f"  {float(eml_mass_limit(8.0*EML_PER_M_800, u, L_sc))/1e3:13,.0f}"
        print(line)
    print(f"  The worst case is standing still (u = 0, the sled at rest on the casing).  The width that carries 5,000 t at rest:")
    for L in (50.0, 100.0, 200.0):
        F = float(eml_required(M_SC0, 0.0, L))
        print(f"    L_sc = {L:4g} m: {F/1e3:7.1f} kN/m -> {F/EML_PER_M_400:5.2f} m of track at 400 A, {F/EML_PER_M_800:5.2f} m at 800 A, "
              f"{F/EML_PER_M_SAT:5.2f} m at plate saturation")

    # ------------------------------------------------------------ 2. limiter strength and reach
    hr("LIMIT 2: BOTTOM LIMITER STRENGTH (ampere-turns) AND REACH (wavelength), 5,000 t, L_sc = %g m" % L_sc)
    print(f"  Present: K0 = {lim0['K0']/1e3:,.0f} kA/m (40 coils/m x 24 turns x 750 A), wavelength {lim0['wavelength']*1e3:.0f} mm, "
          f"F0 = {lim0['F0']/1e3:.1f} kN/m, stored {lim0['U_all']/1e3:.2f} kJ/m.")
    print(f"  F0 goes as K^2.  Reach: decay length wavelength / 2 pi; the force left at the nominal gap and its heating and drag are")
    print(f"  the price of reach.  Top speed is the largest u with the bottom gap >= 20 mm (first number) and with no contact (second).")
    k_mults = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
    wls = [0.10, 0.15, 0.20, 0.30, 0.40]
    print(f"\n  {'K mult':>6} {'K kA/m':>7} {'lambda mm':>9} {'F0 kN/m':>8} {'U kJ/m':>7} {'U to 20mm':>9} {'F nom kN/m':>10} "
          f"{'Q nom kW/m':>10} {'drag nom N/m':>12} {'u 20mm':>7} {'u touch':>8}")
    lim_map = {}
    t0 = time.time()
    for km in k_mults:
        for wl in wls:
            lim = limiter_design(km, wl)
            u20 = top_speed(lim, M_SC0, L_sc, GAPF)
            unc = top_speed(lim, M_SC0, L_sc, 0.0005)
            lim_map[(km, wl)] = (lim, u20, unc)
            print(f"  {km:6.1f} {lim['K']/1e3:7,.0f} {wl*1e3:9.0f} {lim['F0']/1e3:8.0f} {lim['U_all']/1e3:7.2f} {lim['U_to_20']/1e3:9.2f} "
                  f"{lim['F_nom']/1e3:10.2f} {lim['Q_nom']/1e3:10.2f} {lim['D_nom']:12.2f} {u20/1e3:7.1f} {unc/1e3:8.1f}")
    print(f"  ({time.time()-t0:.0f} s)")
    print(f"  Q nom is the steady eddy heating of the reaction plate at the nominal gap, paid all the time, not per pass.")
    print(f"  U to 20 mm is the energy the limiter can take before the gap passes 20 mm; the 20 mm criterion is what makes long")
    print(f"  reach pay, because with a 100 mm wavelength only {lim0['U_to_20']/lim0['U_all']*100:.0f} % of the stored energy is available above 20 mm.")

    # ------------------------------------------------------------ 3. mass against speed
    hr("LIMIT 3: SPACECRAFT MASS AGAINST SLED SPEED, L_sc = %g m, ideal control, 871 kN/m ceiling" % L_sc)
    u_env = np.arange(8000.0, sp.U_MAX + 1, 1000.0)
    masses = np.concatenate([np.arange(100e3, 2e6, 50e3), np.arange(2e6, 20.001e6, 100e3)])
    designs = [("present (K0, 100 mm)", lim0),
               ("K x2, 100 mm", limiter_design(2.0, 0.10)),
               ("K x2, 200 mm", limiter_design(2.0, 0.20)),
               ("K x3, 300 mm", limiter_design(3.0, 0.30))]
    env = {}
    t0 = time.time()
    for name, lim in designs:
        env[name] = mass_envelope(lim, u_env, L_sc, masses)
    print(f"  ({time.time()-t0:.0f} s)")
    print(f"  Largest spacecraft mass (t) for the bottom gap >= 20 mm / for no contact.  Radial load a/g0 is set by u alone.")
    print(f"  {'u km/s':>7} {'v_i km/s':>8} {'g':>6} " + " ".join(f"{name:>26}" for name, _ in designs))
    for j, u in enumerate(u_env):
        line = f"  {u/1e3:7.0f} {(u+sp.V_GS)/1e3:8.2f} {-float(sp.a_down(u))/G0:6.2f} "
        for name, _ in designs:
            M20, Mnc = env[name][0][j], env[name][1][j]
            def lab(m):
                if not np.isfinite(m):
                    return "-"
                if m >= masses[-1]:
                    return ">%.0f" % (m / 1e3)
                return "%.0f" % (m / 1e3)
            line += f"{('%7s / %7s' % (lab(M20), lab(Mnc))):>26}"
        print(line)
    print(f"  '-' means not even {masses[0]/1e3:.0f} t passes; '>' means the top of the mass grid, {masses[-1]/1e3:,.0f} t, passes.")
    print(f"  The low-speed side of the envelope is the EML (limit 1): below the neutral speed the mass limit is")
    print(f"  L_sc * ((F_max - weight)/a_down - 525), which with {EML_WIDTH0:.0f} m of track at 400 A and L_sc = {L_sc:g} m is "
          f"{float(eml_mass_limit(sp.F_EML_OP, 0.0, L_sc))/1e3:,.0f} t at rest and {float(eml_mass_limit(sp.F_EML_OP, 3000.0, L_sc))/1e3:,.0f} t at 3 km/s.")

    # ------------------------------------------------------------ graphs
    if graphs:
        hr("GRAPHS")
        # 8: EML width
        fig, axes = plt.subplots(1, 2, figsize=(sp.GRAPH_W * 1.3, sp.GRAPH_H))
        u_lo = np.linspace(0, 7000.0, 300)
        for L, col in zip((50.0, 100.0, 200.0), PALETTE):
            axes[0].plot(u_lo / 1e3, eml_required(M_SC0, u_lo, L) / EML_PER_M_400, color=col, lw=1.5, label=f"L_sc = {L:g} m, 400 A")
            axes[0].plot(u_lo / 1e3, eml_required(M_SC0, u_lo, L) / EML_PER_M_800, color=col, lw=1.2, ls="--", label=f"L_sc = {L:g} m, 800 A")
        axes[0].axhline(EML_WIDTH0, color="gray", ls=":", lw=1, label="present, 2 m of track")
        axes[0].set_xlabel("Sled speed relative to the casing (km/s)")
        axes[0].set_ylabel("Total EML track width needed for a 5,000 t spacecraft (m)")
        axes[0].set_ylim(0, 10)
        axes[0].grid(alpha=0.3)
        axes[0].legend(fontsize=7)
        for wd, col in zip(widths, PALETTE):
            m = eml_mass_limit(wd * EML_PER_M_400, u_lo, L_sc)
            axes[1].plot(u_lo / 1e3, m / 1e3, color=col, lw=1.5, label=f"{wd:g} m of track at 400 A")
        axes[1].axhline(M_SC0 / 1e3, color="gray", ls=":", lw=1, label="5,000 t")
        axes[1].set_xlabel("Sled speed relative to the casing (km/s)")
        axes[1].set_ylabel(f"Largest spacecraft mass the EML can carry, L_sc = {L_sc:g} m (t)")
        axes[1].set_ylim(0, 20000)
        axes[1].grid(alpha=0.3)
        axes[1].legend(fontsize=7)
        sp._save(fig, "fig_sled_08_eml_width.png", out_dir)

        # 9: limiter map
        fig, axes = plt.subplots(1, 2, figsize=(sp.GRAPH_W * 1.3, sp.GRAPH_H))
        for wl, col in zip(wls, PALETTE):
            axes[0].plot(k_mults, [lim_map[(km, wl)][1] / 1e3 for km in k_mults], color=col, marker="o", ms=4, lw=1.4,
                         label=f"wavelength {wl*1e3:.0f} mm, gap >= 20 mm")
            axes[0].plot(k_mults, [lim_map[(km, wl)][2] / 1e3 for km in k_mults], color=col, marker="s", ms=4, lw=1.2, ls="--",
                         label=f"wavelength {wl*1e3:.0f} mm, no contact")
            axes[1].plot(k_mults, [lim_map[(km, wl)][0]["Q_nom"] / 1e3 for km in k_mults], color=col, marker="o", ms=4, lw=1.4,
                         label=f"wavelength {wl*1e3:.0f} mm")
        axes[0].axhline(30, color="gray", ls=":", lw=1)
        axes[0].set_xlabel("Limiter ampere-turns, multiple of the present 720 kA/m")
        axes[0].set_ylabel("Top sled speed for a 5,000 t spacecraft (km/s)")
        axes[0].set_ylim(8, 31)
        axes[0].grid(alpha=0.3)
        axes[0].legend(fontsize=6.5, ncol=2)
        axes[1].set_xlabel("Limiter ampere-turns, multiple of the present 720 kA/m")
        axes[1].set_ylabel("Steady reaction plate heating at the nominal gap (kW/m)")
        axes[1].set_yscale("log")
        axes[1].grid(alpha=0.3, which="both")
        axes[1].legend(fontsize=7)
        sp._save(fig, "fig_sled_09_limiter_map.png", out_dir)

        # 10: mass envelope
        fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
        for (name, _), col in zip(designs, PALETTE):
            M20, Mnc = env[name][0], env[name][1]
            ax.plot(u_env / 1e3, M20 / 1e3, color=col, lw=1.5, marker="o", ms=3, label=f"{name}, gap >= 20 mm")
            ax.plot(u_env / 1e3, Mnc / 1e3, color=col, lw=1.2, ls="--", marker="s", ms=3, label=f"{name}, no contact")
        ax.axhline(M_SC0 / 1e3, color="gray", ls=":", lw=1, label="5,000 t")
        ax.set_xlabel("Sled speed relative to the casing (km/s)")
        ax.set_ylabel(f"Largest spacecraft mass (t), L_sc = {L_sc:g} m")
        ax.set_yscale("log")
        ax.set_ylim(100, 20000)
        ax.grid(alpha=0.3, which="both")
        ax2 = ax.secondary_xaxis("top", functions=(lambda u: -sp.a_down(u * 1e3) / G0,
                                                   lambda g: (np.sqrt(np.maximum((sp.G_LOCAL + g * G0) * sp.R_ORBIT, 0.0)) - sp.V_GS) / 1e3))
        ax2.set_xlabel("Radial load on the sled (g)")
        ax.legend(fontsize=7, ncol=2)
        sp._save(fig, "fig_sled_10_mass_envelope.png", out_dir)

    hr(f"DONE in {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
