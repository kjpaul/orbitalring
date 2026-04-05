#!/usr/bin/env python3
"""
Minimum Safe Altitude — Orbital Ring Catastrophic Failure Analysis

Computes the minimum orbital ring altitude such that a catastrophic
structural failure does not deliver more than 50 million tonnes of
carbon debris into the atmosphere (the threshold for a nuclear-winter-
class event).

Physics:
  1. Ring at altitude H orbits at v_orb = sqrt(GM/r).
  2. Catastrophic cable-casing collision at v_rel ~ 9 km/s shatters
     the structure into fragments.
  3. Fragments receiving retrograde delta_v lower their orbital perigee.
  4. If perigee < 500 km, atmospheric drag deorbits the fragment.
  5. Ejecta follow a power-law mass-velocity distribution:
       M(> delta_v) = M_retro * (delta_v / v_min)^(-beta)
     where M_retro = f_retrograde * M_total.

Adjustable parameters (see EJECTA DISTRIBUTION section below):
  BETA         Power-law exponent. Range 1.2-2.0.
                 1.2 = violent fragmentation (more mass at high velocity)
                 2.0 = mild fragmentation (mass concentrated at low velocity)
  V_MIN        Minimum fragment delta_v (m/s). Default 10 m/s.
  V_MAX        Maximum fragment delta_v (m/s). Default 4500 m/s.
  F_RETROGRADE Fraction of total mass scattered retrograde. Default 0.5.

Usage:
    python debris_simulation.py [keywords] [options]

    Keywords:
      sweep        Mass vs. altitude for multiple beta values
      safe         Minimum safe altitude with threshold crossing
      delta_v      Critical delta_v vs. altitude
      parametric   Safe altitude vs. beta exponent
      all          Generate all plots

    Options:
      --save       Save plots to graphs_debris/
      --show       Display plots interactively
      --beta=N     Override power-law exponent (e.g. --beta=1.8)
      --help, -h   Show this help

References:
  Housen & Holsapple (2011). Ejecta from impact craters.
    Icarus 211(1), 856-875.
  Melosh (1989). Impact Cratering: A Geologic Process.
    Oxford Univ. Press. Ch. 7 (ejecta scaling laws).
"""

import sys
import os
import math
import csv

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt


# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

GM = 3.986004418e14          # Earth gravitational parameter (m^3/s^2)
R_EARTH = 6_378_137          # Earth equatorial radius (m)


# =============================================================================
# RING PARAMETERS
# =============================================================================

M_RING = 5.7e12              # Total ring mass (kg)
V_REL = 9000.0               # Cable-casing relative collision velocity (m/s)
L_RING = 41_646_000          # Ring circumference (m)
E_PER_M = 446e9              # Collision energy per meter (J/m)
H_DRAG = 500e3               # Atmospheric drag threshold altitude (m)
M_THRESHOLD = 5.0e10         # Nuclear winter threshold (kg) = 50 Mt


# =============================================================================
# EJECTA DISTRIBUTION PARAMETERS
#
# The power-law M(>v) = M_retro * (v/v_min)^(-beta) is the cumulative
# mass of fragments with retrograde delta_v greater than v.
#
# To make the fragmentation MORE violent (more mass at high delta_v):
#   - Decrease BETA toward 1.2
#   - Increase V_MAX toward V_REL (9000 m/s)
#   - Increase F_RETROGRADE toward 1.0
#
# To make it LESS violent:
#   - Increase BETA toward 2.0 or higher
#   - Decrease V_MAX
#   - Decrease F_RETROGRADE
# =============================================================================

BETA = 1.5                   # Power-law exponent (default: moderate)
V_MIN = 10.0                 # Minimum fragment delta_v (m/s)
V_MAX = 4500.0               # Maximum fragment delta_v ~ V_REL / 2 (m/s)
F_RETROGRADE = 0.5           # Fraction of mass scattered retrograde


# =============================================================================
# GRAPH OUTPUT
# =============================================================================

GRAPH_DPI = 200
GRAPH_DIR = "graphs_debris"
SAVE_GRAPHS = True


# =============================================================================
# SECTION 1: ORBITAL MECHANICS
# =============================================================================

def v_orbital(h):
    """Keplerian circular orbital velocity at altitude h (m) above surface."""
    return math.sqrt(GM / (R_EARTH + h))


def perigee_altitude(h_ring, delta_v):
    """Perigee altitude (m) for a fragment at h_ring given retrograde delta_v.

    The fragment's velocity becomes v_frag = v_orbit - delta_v (tangential).
    The collision point is now the apogee of the new elliptical orbit.
    Vis-viva gives the semi-major axis; perigee = 2a - r_apogee.
    """
    r = R_EARTH + h_ring
    v_orb = math.sqrt(GM / r)
    v_frag = v_orb - delta_v

    if v_frag <= 0:
        return -R_EARTH

    # Semi-major axis from vis-viva: v^2 = GM * (2/r - 1/a)
    a = 1.0 / (2.0 / r - v_frag**2 / GM)

    if a <= 0:
        return -R_EARTH

    # Perigee: for an ellipse, r_perigee + r_apogee = 2a
    r_p = 2 * a - r
    return r_p - R_EARTH


def delta_v_critical(h_ring, h_drag=H_DRAG):
    """Minimum retrograde delta_v to drop a fragment's perigee to h_drag.

    For an orbit with apogee at r_ring and perigee at r_drag:
      a = (r_ring + r_drag) / 2
      v_apogee = sqrt(GM * (2/r_ring - 1/a))
      delta_v = v_orbit - v_apogee
    """
    if h_ring <= h_drag:
        return 0.0

    r_ring = R_EARTH + h_ring
    r_drag = R_EARTH + h_drag

    a = (r_ring + r_drag) / 2.0
    v_orb = math.sqrt(GM / r_ring)
    v_apogee = math.sqrt(GM * (2.0 / r_ring - 1.0 / a))

    return v_orb - v_apogee


# =============================================================================
# SECTION 2: EJECTA MASS-VELOCITY DISTRIBUTION
# =============================================================================

def mass_below_drag_zone(h_ring, beta=BETA, v_min=V_MIN, v_max=V_MAX,
                         f_retro=F_RETROGRADE, m_total=M_RING):
    """Total debris mass whose perigee drops below H_DRAG.

    The cumulative mass of retrograde fragments with delta_v > v is:
      M(>v) = M_retro * (v / v_min)^(-beta)    for v_min <= v <= v_max

    At v = v_min, M = M_retro (all retrograde mass).
    At v > v_max, M = 0 (no fragments that fast).
    """
    if h_ring <= H_DRAG:
        return f_retro * m_total

    dv = delta_v_critical(h_ring)

    if dv <= v_min:
        # Even the slowest fragments drop perigee below drag zone
        return f_retro * m_total
    if dv >= v_max:
        return 0.0

    m_retro = f_retro * m_total
    return m_retro * (dv / v_min) ** (-beta)


def ejecta_energy_check(beta=BETA, v_min=V_MIN, v_max=V_MAX,
                        f_retro=F_RETROGRADE, m_total=M_RING):
    """Compute mean square velocity and total KE of the ejecta distribution.

    Returns (v_rms, ke_total, ke_fraction) where ke_fraction is the
    fragment KE as a fraction of the total collision energy.
    """
    # The PDF of the power law M(>v) = C * v^(-beta) is:
    #   p(v) = beta * v_min^beta * v^(-beta-1)   for v_min <= v <= v_max
    #
    # <v^2> = integral of v^2 * p(v) dv from v_min to v_max

    if abs(beta - 2.0) < 1e-10:
        # Special case: beta = 2 gives a log integral
        mean_v2 = 2.0 * v_min**2 * math.log(v_max / v_min)
    else:
        k = 2.0 - beta
        mean_v2 = (beta / k) * v_min**beta * (v_max**k - v_min**k)

    v_rms = math.sqrt(mean_v2)

    # Total KE of retrograde fragments (same for prograde by symmetry)
    m_retro = f_retro * m_total
    ke_retro = 0.5 * m_retro * mean_v2
    ke_total = 2 * ke_retro  # both directions

    e_collision = E_PER_M * L_RING
    ke_fraction = ke_total / e_collision

    return v_rms, ke_total, ke_fraction


# =============================================================================
# SECTION 3: ROOT FINDING — MINIMUM SAFE ALTITUDE
# =============================================================================

def find_safe_altitude(target_mass=M_THRESHOLD, beta=BETA,
                       h_low=501e3, h_high=20_000e3):
    """Find altitude where debris mass equals target_mass.

    Uses Brent's method on:  mass_below_drag_zone(h) - target_mass = 0
    """
    def objective(h):
        return mass_below_drag_zone(h, beta=beta) - target_mass

    f_low = objective(h_low)
    f_high = objective(h_high)

    if f_low <= 0:
        return h_low
    if f_high >= 0:
        return h_high

    return brentq(objective, h_low, h_high, rtol=1e-12)


# =============================================================================
# SECTION 4: CONSOLE OUTPUT
# =============================================================================

def print_parameters(beta=BETA):
    """Print simulation parameters and key results."""
    sys.stdout.reconfigure(encoding='utf-8')

    print("=" * 70)
    print("MINIMUM SAFE ALTITUDE — CATASTROPHIC FAILURE ANALYSIS")
    print("=" * 70)

    print(f"\n{'Ring Parameters':}")
    print(f"  M_RING          {M_RING:.3e} kg ({M_RING/1e9:.1f} billion kg)")
    print(f"  V_REL           {V_REL:.0f} m/s")
    print(f"  L_RING          {L_RING:,.0f} m ({L_RING/1e3:.0f} km)")
    print(f"  E_PER_M         {E_PER_M/1e9:.0f} GJ/m")
    print(f"  E_TOTAL         {E_PER_M * L_RING:.3e} J")
    print(f"  H_DRAG          {H_DRAG/1e3:.0f} km")
    print(f"  M_THRESHOLD     {M_THRESHOLD:.2e} kg ({M_THRESHOLD/1e9:.0f} Mt)")

    print(f"\n{'Ejecta Distribution':}")
    print(f"  BETA            {beta:.2f}")
    print(f"  V_MIN           {V_MIN:.0f} m/s")
    print(f"  V_MAX           {V_MAX:.0f} m/s")
    print(f"  F_RETROGRADE    {F_RETROGRADE:.2f}")

    # Energy cross-check
    v_rms, ke_total, ke_frac = ejecta_energy_check(beta=beta)
    print(f"\n{'Energy Cross-Check':}")
    print(f"  v_rms (fragments)    {v_rms:.1f} m/s")
    print(f"  KE (all fragments)   {ke_total:.3e} J")
    print(f"  KE / E_collision     {ke_frac:.4f} ({ke_frac*100:.2f}%)")
    if ke_frac > 1.0:
        print(f"  WARNING: fragment KE exceeds collision energy!")
        print(f"  Reduce V_MAX or increase BETA for energy consistency.")

    # Delta_v at key altitudes
    print(f"\n{'Delta_v to reach 500 km perigee':}")
    print(f"  {'Altitude (km)':>15}  {'delta_v (m/s)':>13}  {'Mass below (Mt)':>15}")
    print(f"  {'-'*15}  {'-'*13}  {'-'*15}")
    for h_km in [600, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000]:
        h = h_km * 1e3
        dv = delta_v_critical(h)
        m = mass_below_drag_zone(h, beta=beta)
        flag = " <-- THRESHOLD" if abs(m - M_THRESHOLD) / M_THRESHOLD < 0.5 else ""
        print(f"  {h_km:>15,}  {dv:>13.1f}  {m/1e9:>15.1f}{flag}")

    # Find and report safe altitude
    h_safe = find_safe_altitude(beta=beta)
    dv_safe = delta_v_critical(h_safe)
    print(f"\n{'RESULT':}")
    print(f"  Minimum safe altitude:  {h_safe/1e3:.1f} km")
    print(f"  Critical delta_v:       {dv_safe:.1f} m/s")
    print(f"  v_orbital at that alt:  {v_orbital(h_safe):.1f} m/s")

    # Parametric: safe altitude vs beta
    print(f"\n{'Sensitivity to beta (power-law exponent)':}")
    print(f"  {'beta':>6}  {'Safe altitude (km)':>18}  {'dv_crit (m/s)':>13}")
    print(f"  {'-'*6}  {'-'*18}  {'-'*13}")
    for b in [1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0]:
        try:
            hs = find_safe_altitude(beta=b)
            dvc = delta_v_critical(hs)
            print(f"  {b:>6.1f}  {hs/1e3:>18.1f}  {dvc:>13.1f}")
        except ValueError:
            print(f"  {b:>6.1f}  {'> 20,000':>18}  {'N/A':>13}")

    print()
    return h_safe


# =============================================================================
# SECTION 5: PLOTTING
# =============================================================================

def save_or_show(fig, filename, save=True, show=False):
    """Save figure to GRAPH_DIR or show interactively."""
    if save:
        os.makedirs(GRAPH_DIR, exist_ok=True)
        path = os.path.join(GRAPH_DIR, filename)
        fig.savefig(path, dpi=GRAPH_DPI, bbox_inches='tight')
        print(f"  Saved: {path}")
    if show:
        plt.show()
    plt.close(fig)


def plot_mass_vs_altitude(save=True, show=False, beta=BETA):
    """Plot 1: Mass entering drag zone vs. ring altitude for multiple beta."""
    print("Generating: mass vs. altitude sweep...")

    h_km = np.linspace(550, 10000, 2000)
    betas = [1.2, 1.5, 1.8, 2.0]
    if beta not in betas:
        betas.append(beta)
        betas.sort()

    fig, ax = plt.subplots(figsize=(10, 6))

    for b in betas:
        mass = np.array([mass_below_drag_zone(h * 1e3, beta=b) for h in h_km])
        mass_mt = mass / 1e9  # megatonnes
        lw = 2.5 if abs(b - beta) < 0.01 else 1.2
        ax.plot(h_km, mass_mt, linewidth=lw, label=f"beta = {b:.1f}")

    # Threshold line
    ax.axhline(M_THRESHOLD / 1e9, color='red', linestyle='--', linewidth=1.5,
               label=f"Extinction threshold ({M_THRESHOLD/1e9:.0f} Mt)")

    # Mark safe altitude for the default beta
    h_safe = find_safe_altitude(beta=beta)
    ax.axvline(h_safe / 1e3, color='gray', linestyle=':', linewidth=1.0, alpha=0.7)
    ax.plot(h_safe / 1e3, M_THRESHOLD / 1e9, 'ko', markersize=8, zorder=5)
    ax.annotate(f"  {h_safe/1e3:.0f} km (beta={beta})",
                xy=(h_safe / 1e3, M_THRESHOLD / 1e9),
                fontsize=10, va='bottom', ha='left')

    ax.set_yscale('log')
    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ax.set_ylabel("Mass Dropping Below 500 km (Mt)", fontsize=12)
    ax.set_title("Debris Mass Entering Drag Zone vs. Ring Altitude", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(500, 10000)
    ax.set_ylim(1e-1, 5e3)

    save_or_show(fig, "01-mass_vs_altitude.png", save, show)


def plot_safe_altitude(save=True, show=False, beta=BETA):
    """Plot 2: Zoomed view of threshold crossing with clear labeling."""
    print("Generating: safe altitude detail...")

    h_safe = find_safe_altitude(beta=beta)
    h_center = h_safe / 1e3
    h_lo = max(550, h_center - 500)
    h_hi = h_center + 500

    h_km = np.linspace(h_lo, h_hi, 1000)
    mass = np.array([mass_below_drag_zone(h * 1e3, beta=beta) for h in h_km])
    mass_mt = mass / 1e9

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(h_km, mass_mt, 'b-', linewidth=2.5, label=f"Debris mass (beta={beta})")

    ax.axhline(M_THRESHOLD / 1e9, color='red', linestyle='--', linewidth=1.5,
               label=f"Extinction threshold ({M_THRESHOLD/1e9:.0f} Mt)")

    ax.axvline(h_safe / 1e3, color='green', linestyle='-', linewidth=2.0, alpha=0.7,
               label=f"Min safe altitude = {h_safe/1e3:.0f} km")
    ax.plot(h_safe / 1e3, M_THRESHOLD / 1e9, 'r*', markersize=15, zorder=5)

    # Shade danger zone
    ax.axvspan(h_lo, h_safe / 1e3, alpha=0.08, color='red')
    ax.axvspan(h_safe / 1e3, h_hi, alpha=0.08, color='green')
    ax.text(h_lo + 20, ax.get_ylim()[0] * 1.5 if ax.get_yscale() == 'log' else 0,
            "DANGER", fontsize=14, color='red', alpha=0.4, va='bottom')
    ax.text(h_hi - 20, ax.get_ylim()[0] * 1.5 if ax.get_yscale() == 'log' else 0,
            "SAFE", fontsize=14, color='green', alpha=0.4, va='bottom', ha='right')

    ax.set_yscale('log')
    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ax.set_ylabel("Mass Dropping Below 500 km (Mt)", fontsize=12)
    ax.set_title("Minimum Safe Altitude — Threshold Crossing Detail", fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    save_or_show(fig, "02-safe_altitude.png", save, show)


def plot_delta_v(save=True, show=False):
    """Plot 3: Critical delta_v vs. altitude."""
    print("Generating: delta_v vs. altitude...")

    h_km = np.linspace(550, 10000, 1000)
    dv = np.array([delta_v_critical(h * 1e3) for h in h_km])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(h_km, dv, 'b-', linewidth=2.0)

    # Mark V_MAX
    ax.axhline(V_MAX, color='red', linestyle='--', linewidth=1.0,
               label=f"V_MAX = {V_MAX:.0f} m/s")
    ax.axhline(V_MIN, color='gray', linestyle=':', linewidth=1.0,
               label=f"V_MIN = {V_MIN:.0f} m/s")

    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ax.set_ylabel("Retrograde delta_v to Reach 500 km Perigee (m/s)", fontsize=12)
    ax.set_title("Critical Delta-v for Atmospheric Re-entry", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    save_or_show(fig, "03-delta_v_vs_altitude.png", save, show)


def plot_parametric(save=True, show=False):
    """Plot 4: Minimum safe altitude vs. beta exponent."""
    print("Generating: parametric beta sweep...")

    betas = np.linspace(1.0, 3.0, 200)
    h_safe_km = []

    for b in betas:
        try:
            hs = find_safe_altitude(beta=b)
            h_safe_km.append(hs / 1e3)
        except (ValueError, RuntimeError):
            h_safe_km.append(np.nan)

    h_safe_km = np.array(h_safe_km)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(betas, h_safe_km, 'b-', linewidth=2.5)

    # Mark default
    h_def = find_safe_altitude(beta=BETA) / 1e3
    ax.plot(BETA, h_def, 'ro', markersize=10, zorder=5,
            label=f"Default: beta={BETA}, H={h_def:.0f} km")

    # Annotate typical ranges
    ax.axvspan(1.2, 2.0, alpha=0.08, color='blue',
               label="Typical hypervelocity range (1.2-2.0)")

    ax.set_xlabel("Power-Law Exponent (beta)", fontsize=12)
    ax.set_ylabel("Minimum Safe Altitude (km)", fontsize=12)
    ax.set_title("Sensitivity of Minimum Safe Altitude to Fragmentation Violence", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    save_or_show(fig, "04-parametric_beta.png", save, show)


# =============================================================================
# SECTION 6: CSV EXPORT
# =============================================================================

def write_csv(beta=BETA):
    """Write results to CSV files."""
    os.makedirs(GRAPH_DIR, exist_ok=True)

    # Altitude sweep
    path = os.path.join(GRAPH_DIR, "debris_altitude_sweep.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["altitude_km", "delta_v_crit_m_s", "mass_below_500km_kg",
                     "mass_below_500km_Mt"])
        for h_km in range(550, 10001, 10):
            h = h_km * 1e3
            dv = delta_v_critical(h)
            m = mass_below_drag_zone(h, beta=beta)
            w.writerow([h_km, f"{dv:.2f}", f"{m:.4e}", f"{m/1e9:.4f}"])
    print(f"  Saved: {path}")

    # Beta sensitivity
    path = os.path.join(GRAPH_DIR, "debris_beta_sensitivity.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["beta", "safe_altitude_km", "delta_v_crit_m_s"])
        for b100 in range(100, 301, 5):
            b = b100 / 100.0
            try:
                hs = find_safe_altitude(beta=b)
                dv = delta_v_critical(hs)
                w.writerow([f"{b:.2f}", f"{hs/1e3:.1f}", f"{dv:.1f}"])
            except (ValueError, RuntimeError):
                w.writerow([f"{b:.2f}", "N/A", "N/A"])
    print(f"  Saved: {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    sys.stdout.reconfigure(encoding='utf-8')

    show_graphs = []
    save = SAVE_GRAPHS
    show = False
    beta = BETA

    # Parse arguments
    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--save":
                save = True
            elif arg == "--show":
                show = True
                save = False
            elif arg.startswith("--beta="):
                try:
                    beta = float(arg.split("=")[1])
                except ValueError:
                    print(f"Invalid beta value: {arg}")
            else:
                show_graphs.append(arg.lower())

    # Always print parameters and results
    h_safe = print_parameters(beta=beta)

    # Generate requested plots
    if not show_graphs:
        return

    keywords = set(show_graphs)
    do_all = "all" in keywords

    print("Generating graphs...")

    if do_all or "sweep" in keywords:
        plot_mass_vs_altitude(save, show, beta)

    if do_all or "safe" in keywords:
        plot_safe_altitude(save, show, beta)

    if do_all or "delta_v" in keywords:
        plot_delta_v(save, show)

    if do_all or "parametric" in keywords:
        plot_parametric(save, show)

    # Write CSV
    if save:
        write_csv(beta)

    print("Done.")


if __name__ == "__main__":
    main()
