#!/usr/bin/env python3
"""
Tether Stress and Taper Profile
================================

Computes the stress distribution, cross-section taper profile, and mass
of a space elevator tether under uniform-stress (optimal taper) design.
Analyses both GEO-to-ground and GEO-to-275 km configurations, including
the effect of orbital inclination on taper ratio (Gassend correction).

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python tether_stress.py                    Print results only
    python tether_stress.py all                Generate all graphs
    python tether_stress.py stress_profile     Generate stress profile graph
    python tether_stress.py taper_profile      Generate taper profile graph
    python tether_stress.py taper_vs_strength  Generate taper vs strength graph
    python tether_stress.py --show             Display graphs interactively
    python tether_stress.py --save             Save graphs to files (default)
    python tether_stress.py --help, -h         Show this help

Graph keywords:
    stress_profile     Stress vs altitude for both configurations
    taper_profile      Cross-section ratio vs altitude (log scale)
    taper_vs_strength  Taper ratio vs material breaking strength
    all                Generate all graphs
"""

import sys
import os
import math
import csv

import matplotlib.pyplot as plt

import tether_config as cfg


# =============================================================================
# NET RADIAL ACCELERATION
# =============================================================================

def a_net(r):
    """Net downward acceleration at radius r (positive toward Earth).

    a_net(r) = GM/r^2 - omega^2 * r

    Below GEO, gravity dominates so a_net > 0 (net force toward Earth).
    Above GEO, centrifugal dominates so a_net < 0 (net force outward).
    At GEO, a_net = 0 by definition.

    Args:
        r: Orbital radius (m)

    Returns:
        Net downward acceleration (m/s^2), positive toward Earth.
    """
    return cfg.GM / r**2 - cfg.OMEGA_SIDEREAL**2 * r


# =============================================================================
# ANALYTICAL INTEGRAL OF a_net
# =============================================================================

def phi_integral(r_bottom, r_top):
    """Analytical integral of a_net(r) from r_bottom to r_top.

    Integral of (GM/r^2 - omega^2 * r) dr from r_bottom to r_top:

        [-GM/r - omega^2 * r^2 / 2] evaluated from r_bottom to r_top
        = GM * (1/r_bottom - 1/r_top) - omega^2 * (r_top^2 - r_bottom^2) / 2

    This integral appears in the exponent of the uniform-stress taper
    profile.  Its sign determines whether the cross-section increases
    or decreases from bottom to top.

    Args:
        r_bottom: Lower radius (m)
        r_top:    Upper radius (m)

    Returns:
        Value of the integral (m^2/s^2), same units as specific energy.
    """
    omega = cfg.OMEGA_SIDEREAL
    return (cfg.GM * (1.0 / r_bottom - 1.0 / r_top)
            - omega**2 * (r_top**2 - r_bottom**2) / 2.0)


# =============================================================================
# TAPER RATIO
# =============================================================================

def taper_ratio(r_bottom, sigma=None, rho=None):
    """Taper ratio from r_bottom to GEO for uniform-stress tether.

    The taper ratio is A(GEO) / A(r_bottom) for a tether designed so
    that every cross-section operates at the same stress sigma.

        TR = exp(rho / sigma * Phi)

    where Phi = integral of a_net from r_bottom to R_GEO.

    For the below-GEO segment, Phi > 0 and the cross-section increases
    from the bottom to GEO (the tether must support its own weight plus
    the payload below).

    Args:
        r_bottom: Bottom attachment radius (m)
        sigma:    Operating stress (Pa), default cfg.SIGMA_OPERATING
        rho:      Material density (kg/m^3), default cfg.RHO_CNT

    Returns:
        Taper ratio A(GEO) / A(bottom), dimensionless (>= 1).
    """
    if sigma is None:
        sigma = cfg.SIGMA_OPERATING
    if rho is None:
        rho = cfg.RHO_CNT

    Phi = phi_integral(r_bottom, cfg.R_GEO)
    return math.exp(rho / sigma * Phi)


# =============================================================================
# TAPER PROFILE
# =============================================================================

def taper_profile(r_bottom, r_top, A_bottom, n_points=None):
    """Compute cross-section area vs radius for uniform-stress tether.

    The uniform-stress design sets:
        A(r) = A_bottom * exp(rho/sigma * integral(a_net, r_bottom, r))

    For the below-GEO segment (r_bottom < r < R_GEO), the cross-section
    increases with radius because the tether must support all material
    below.

    For the above-GEO segment, call with r_bottom = R_GEO and r_top > R_GEO.
    The outward centrifugal force dominates, and the cross-section increases
    with radius above GEO as well (tether must support material above).

    Args:
        r_bottom:  Starting radius (m)
        r_top:     Ending radius (m)
        A_bottom:  Cross-section area at r_bottom (m^2)
        n_points:  Number of points, default cfg.N_POINTS_TETHER

    Returns:
        (r_array, A_array): Lists of radius (m) and area (m^2).
    """
    if n_points is None:
        n_points = cfg.N_POINTS_TETHER

    rho = cfg.RHO_CNT
    sigma = cfg.SIGMA_OPERATING

    r_array = []
    A_array = []
    dr = (r_top - r_bottom) / (n_points - 1)

    for i in range(n_points):
        r = r_bottom + i * dr
        Phi = phi_integral(r_bottom, r)
        A = A_bottom * math.exp(rho / sigma * Phi)
        r_array.append(r)
        A_array.append(A)

    return r_array, A_array


# =============================================================================
# TETHER MASS
# =============================================================================

def tether_mass(r_bottom, A_bottom, n_points=None):
    """Total tether mass from r_bottom to GEO (uniform-stress design).

    M = integral(rho * A(r) dr) from r_bottom to R_GEO

    Uses trapezoidal numerical integration over the taper profile.

    Args:
        r_bottom: Bottom radius (m)
        A_bottom: Cross-section area at r_bottom (m^2)
        n_points: Number of integration points, default cfg.N_POINTS_TETHER

    Returns:
        Total tether mass (kg).
    """
    if n_points is None:
        n_points = cfg.N_POINTS_TETHER

    r_arr, A_arr = taper_profile(r_bottom, cfg.R_GEO, A_bottom, n_points)
    rho = cfg.RHO_CNT

    # Trapezoidal integration
    mass = 0.0
    for i in range(len(r_arr) - 1):
        dr = r_arr[i + 1] - r_arr[i]
        mass += 0.5 * rho * (A_arr[i] + A_arr[i + 1]) * dr

    return mass


# =============================================================================
# MAXIMUM PAYLOAD
# =============================================================================

def max_payload(r_bottom, A_bottom):
    """Maximum payload supportable at r_bottom for given bottom cross-section.

    The bottom cross-section must support the payload weight against the
    net downward acceleration at r_bottom:

        F_payload = m_payload * a_net(r_bottom)
        A_bottom  = F_payload / sigma_operating
        m_payload = A_bottom * sigma_operating / a_net(r_bottom)

    Args:
        r_bottom: Bottom attachment radius (m)
        A_bottom: Cross-section area at bottom (m^2)

    Returns:
        Maximum payload mass (kg).
    """
    g_eff = a_net(r_bottom)
    if g_eff <= 0:
        return float('inf')  # No net downward force (at or above GEO)
    return A_bottom * cfg.SIGMA_OPERATING / g_eff


# =============================================================================
# INCLINED ORBIT TAPER RATIO (GASSEND CORRECTION)
# =============================================================================

def taper_ratio_inclined(inclination_deg, r_bottom, sigma=None, rho=None):
    """Taper ratio for inclined geosynchronous orbit (Gassend approximation).

    For non-equatorial geosynchronous orbits, the tether is not purely
    radial — it traces a figure-eight as seen from the ground.  This
    increases the effective path length and adds a lateral stress component.

    The Gassend correction factor is approximately:

        TR_inclined = TR_equatorial * exp(rho * omega^2 * R_GEO^2 * sin^2(i) / (2 * sigma))

    This is an approximation valid for small to moderate inclinations.
    At large inclinations the actual geometry becomes more complex, and
    the tether may need active steering or multiple anchor points.

    Note: This approximation comes from Gassend's analysis of inclined
    space elevator tethers.  The correction captures the additional stress
    from the tether's lateral displacement but does not account for all
    higher-order effects.

    Args:
        inclination_deg: Orbital inclination (degrees)
        r_bottom:        Bottom attachment radius (m)
        sigma:           Operating stress (Pa), default cfg.SIGMA_OPERATING
        rho:             Material density (kg/m^3), default cfg.RHO_CNT

    Returns:
        Taper ratio for inclined orbit, dimensionless.
    """
    if sigma is None:
        sigma = cfg.SIGMA_OPERATING
    if rho is None:
        rho = cfg.RHO_CNT

    TR_eq = taper_ratio(r_bottom, sigma, rho)

    i_rad = math.radians(inclination_deg)
    omega = cfg.OMEGA_SIDEREAL
    correction = math.exp(rho * omega**2 * cfg.R_GEO**2
                          * math.sin(i_rad)**2 / (2.0 * sigma))

    return TR_eq * correction


# =============================================================================
# STRESS PROFILE (for graphing)
# =============================================================================

def stress_profile(r_bottom, r_top, n_points=None):
    """Compute stress vs radius for a uniform-stress tether.

    In the uniform-stress design, every cross-section operates at sigma_operating.
    However, it is instructive to compute what the stress WOULD be in a
    uniform cross-section tether (no tapering).  This shows why tapering is
    essential.

    For a uniform cross-section tether hanging from GEO, the stress at
    radius r is:

        sigma(r) = rho * |integral(a_net, r, R_GEO)|

    This is the stress from supporting the weight of tether material
    between r and GEO.

    Args:
        r_bottom: Lower radius (m)
        r_top:    Upper radius (m)
        n_points: Number of points, default cfg.N_POINTS_TETHER

    Returns:
        (r_array, sigma_array): Lists of radius (m) and stress (Pa).
    """
    if n_points is None:
        n_points = cfg.N_POINTS_TETHER

    rho = cfg.RHO_CNT

    r_array = []
    sigma_array = []
    dr = (r_top - r_bottom) / (n_points - 1)

    for i in range(n_points):
        r = r_bottom + i * dr
        # Stress at r is from supporting tether between r and GEO
        Phi = phi_integral(r, cfg.R_GEO)
        sigma = rho * abs(Phi)
        r_array.append(r)
        sigma_array.append(sigma)

    return r_array, sigma_array


# =============================================================================
# GRAPH: STRESS PROFILE
# =============================================================================

def plot_stress_profile(show_graph):
    """Plot stress vs altitude for uniform cross-section tether.

    Shows the stress that would develop in a constant-area tether
    for both GEO-to-ground and GEO-to-275km configurations, with
    the operating stress limit marked.
    """
    # Configuration 1: GEO to ground
    r_arr_g, sigma_arr_g = stress_profile(cfg.R_BOTTOM_GROUND, cfg.R_GEO)
    alt_g = [(r - cfg.R_E) / 1000.0 for r in r_arr_g]
    sigma_g_gpa = [s / 1e9 for s in sigma_arr_g]

    # Configuration 2: GEO to 275 km
    r_arr_s, sigma_arr_s = stress_profile(cfg.R_BOTTOM_SHORT, cfg.R_GEO)
    alt_s = [(r - cfg.R_E) / 1000.0 for r in r_arr_s]
    sigma_s_gpa = [s / 1e9 for s in sigma_arr_s]

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.plot(alt_g, sigma_g_gpa, 'b-', linewidth=1.5,
            label='GEO to ground')
    ax.plot(alt_s, sigma_s_gpa, 'r-', linewidth=1.5,
            label=f'GEO to {cfg.ALT_BOTTOM_SHORT/1000:.0f} km')
    ax.axhline(y=cfg.SIGMA_OPERATING / 1e9, color='green', linestyle='--',
               linewidth=1.0,
               label=f'Operating stress ({cfg.SIGMA_OPERATING/1e9:.1f} GPa)')
    ax.axhline(y=cfg.SIGMA_BREAK_CNT / 1e9, color='red', linestyle=':',
               linewidth=1.0,
               label=f'Breaking stress ({cfg.SIGMA_BREAK_CNT/1e9:.0f} GPa)')

    ax.set_xlabel('Altitude (km)', fontsize=14)
    ax.set_ylabel('Stress in uniform cross-section tether (GPa)', fontsize=14)
    ax.set_title('Tether Stress vs Altitude (uniform cross-section, no tapering)',
                 fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, cfg.ALT_GEO / 1000.0)

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"01-stress_profile.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: TAPER PROFILE
# =============================================================================

def plot_taper_profile(show_graph):
    """Plot cross-section taper ratio vs altitude (log scale).

    Shows A(r)/A_min for both configurations, where A_min is the minimum
    cross-section (at the bottom of the tether).
    """
    A_ref = 1.0  # Reference area = 1 m^2

    # Configuration 1: GEO to ground
    r_arr_g, A_arr_g = taper_profile(cfg.R_BOTTOM_GROUND, cfg.R_GEO, A_ref)
    alt_g = [(r - cfg.R_E) / 1000.0 for r in r_arr_g]
    ratio_g = [A / A_ref for A in A_arr_g]

    # Configuration 2: GEO to 275 km
    r_arr_s, A_arr_s = taper_profile(cfg.R_BOTTOM_SHORT, cfg.R_GEO, A_ref)
    alt_s = [(r - cfg.R_E) / 1000.0 for r in r_arr_s]
    ratio_s = [A / A_ref for A in A_arr_s]

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.semilogy(alt_g, ratio_g, 'b-', linewidth=1.5,
                label='GEO to ground')
    ax.semilogy(alt_s, ratio_s, 'r-', linewidth=1.5,
                label=f'GEO to {cfg.ALT_BOTTOM_SHORT/1000:.0f} km')

    ax.set_xlabel('Altitude (km)', fontsize=14)
    ax.set_ylabel('Cross-section ratio A(r) / A_bottom', fontsize=14)
    ax.set_title('Uniform-Stress Taper Profile (CNT, safety factor = '
                 f'{cfg.SAFETY_FACTOR:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0, cfg.ALT_GEO / 1000.0)

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"02-taper_profile.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: TAPER RATIO VS MATERIAL STRENGTH
# =============================================================================

def plot_taper_vs_strength(show_graph):
    """Plot taper ratio vs material breaking strength.

    Sweeps breaking strength from 5 to 130 GPa (with safety factor 2)
    and computes the resulting taper ratio for GEO-to-ground.
    """
    sigma_break_gpa = cfg.SIGMA_SWEEP_GPA
    tr_ground = []
    tr_short = []

    for s_gpa in sigma_break_gpa:
        sigma_op = s_gpa * 1e9 / cfg.SAFETY_FACTOR
        tr_g = taper_ratio(cfg.R_BOTTOM_GROUND, sigma=sigma_op)
        tr_s = taper_ratio(cfg.R_BOTTOM_SHORT, sigma=sigma_op)
        tr_ground.append(tr_g)
        tr_short.append(tr_s)

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.semilogy(sigma_break_gpa, tr_ground, 'b-o', linewidth=1.5,
                markersize=5, label='GEO to ground')
    ax.semilogy(sigma_break_gpa, tr_short, 'r-s', linewidth=1.5,
                markersize=5, label=f'GEO to {cfg.ALT_BOTTOM_SHORT/1000:.0f} km')

    # Mark CNT operating point
    ax.axvline(x=cfg.SIGMA_BREAK_CNT / 1e9, color='green', linestyle='--',
               linewidth=1.0,
               label=f'CNT ({cfg.SIGMA_BREAK_CNT/1e9:.0f} GPa)')

    # Horizontal line at taper ratio = 1 (ideal)
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.5)

    ax.set_xlabel('Material breaking strength (GPa)', fontsize=14)
    ax.set_ylabel('Taper ratio A(GEO) / A(bottom)', fontsize=14)
    ax.set_title(f'Taper Ratio vs Material Strength '
                 f'(safety factor = {cfg.SAFETY_FACTOR:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"03-taper_vs_strength.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# CONSOLE OUTPUT
# =============================================================================

def print_results():
    """Print summary of tether stress and taper analysis."""
    sep = "=" * 70

    # --- Tether parameters ---
    print(f"\n{sep}")
    print("TETHER STRESS AND TAPER PROFILE")
    print(sep)
    print(f"  {'Material':28} CNT")
    print(f"  {'Density':28} {cfg.RHO_CNT:>12,.0f} kg/m^3")
    print(f"  {'Breaking stress':28} {cfg.SIGMA_BREAK_CNT/1e9:>12.1f} GPa")
    print(f"  {'Safety factor':28} {cfg.SAFETY_FACTOR:>12.1f}")
    print(f"  {'Operating stress':28} {cfg.SIGMA_OPERATING/1e9:>12.1f} GPa")
    print(f"  {'Characteristic velocity':28} {cfg.V_CHAR_CNT:>12.0f} m/s")
    print(f"  {'Earth sidereal rate':28} {cfg.OMEGA_SIDEREAL:>12.4e} rad/s")
    print(f"  {'GEO radius':28} {cfg.R_GEO/1e6:>12.3f} Mm")
    print(f"  {'GEO altitude':28} {cfg.ALT_GEO/1000:>12.0f} km")

    # --- Configuration comparison ---
    configs = [
        ("GEO to ground", cfg.R_BOTTOM_GROUND, cfg.ALT_BOTTOM_GROUND),
        (f"GEO to {cfg.ALT_BOTTOM_SHORT/1000:.0f} km",
         cfg.R_BOTTOM_SHORT, cfg.ALT_BOTTOM_SHORT),
    ]

    A_ref = 1.0  # 1 m^2 reference bottom area

    print(f"\n{sep}")
    print("TAPER ANALYSIS — BELOW-GEO SEGMENT")
    print(sep)
    print(f"\n  Reference bottom cross-section: {A_ref:.0f} m^2")
    print()

    print(f"  {'Quantity':36} ", end="")
    for label, _, _ in configs:
        print(f"{label:>16}", end="")
    print()
    print(f"  {'-'*36} ", end="")
    for _ in configs:
        print(f"{'-'*16}", end="")
    print()

    results = []
    for label, r_bot, alt_bot in configs:
        tr = taper_ratio(r_bot)
        mass = tether_mass(r_bot, A_ref)
        m_pay = max_payload(r_bot, A_ref)
        g_eff = a_net(r_bot)
        length_km = (cfg.R_GEO - r_bot) / 1000.0
        results.append({
            'label': label,
            'r_bottom': r_bot,
            'alt_bottom': alt_bot,
            'taper_ratio': tr,
            'mass': mass,
            'max_payload': m_pay,
            'g_eff': g_eff,
            'length_km': length_km,
        })

    # Print rows
    row_data = [
        ('Bottom altitude (km)',
         [f"{r['alt_bottom']/1000:>16.0f}" for r in results]),
        ('Bottom radius (Mm)',
         [f"{r['r_bottom']/1e6:>16.3f}" for r in results]),
        ('Tether length (km)',
         [f"{r['length_km']:>16.0f}" for r in results]),
        ('g_eff at bottom (m/s^2)',
         [f"{r['g_eff']:>16.4f}" for r in results]),
        ('Taper ratio A(GEO)/A(bot)',
         [f"{r['taper_ratio']:>16.4f}" for r in results]),
        ('log10(taper ratio)',
         [f"{math.log10(r['taper_ratio']):>16.4f}" for r in results]),
        ('Mass per m^2 bottom (kg)',
         [f"{r['mass']:>16,.0f}" for r in results]),
        ('Mass per m^2 bottom (tonnes)',
         [f"{r['mass']/1000:>16,.1f}" for r in results]),
        ('Max payload per m^2 (kg)',
         [f"{r['max_payload']:>16,.0f}" for r in results]),
        ('Max payload per m^2 (tonnes)',
         [f"{r['max_payload']/1000:>16,.1f}" for r in results]),
    ]

    for row_label, row_vals in row_data:
        print(f"  {row_label:36} ", end="")
        for v in row_vals:
            print(v, end="")
        print()

    # --- Stress check at bottom ---
    print(f"\n{sep}")
    print("STRESS CHECK — UNIFORM CROSS-SECTION (NO TAPERING)")
    print(sep)
    print(f"\n  Stress at bottom if tether had constant cross-section:")

    for r in results:
        Phi = phi_integral(r['r_bottom'], cfg.R_GEO)
        sigma_bottom = cfg.RHO_CNT * abs(Phi)
        ratio_break = sigma_bottom / cfg.SIGMA_BREAK_CNT
        print(f"    {r['label']:28} {sigma_bottom/1e9:8.2f} GPa"
              f"  ({ratio_break:.1f}x breaking stress)")

    # --- Inclined orbit ---
    print(f"\n{sep}")
    print("INCLINED ORBIT — GASSEND CORRECTION")
    print(sep)
    print(f"\n  Taper ratio increase for inclined geosynchronous orbit.")
    print(f"  Approximation: TR_incl = TR_eq * exp(rho*omega^2*R_GEO^2*sin^2(i) / (2*sigma))")
    print(f"  Note: Valid for small-to-moderate inclinations; higher-order")
    print(f"  effects not captured.\n")

    print(f"  {'Inclination':>12} ", end="")
    for label, _, _ in configs:
        print(f"  {'TR ' + label:>20}", end="")
    print(f"  {'Correction':>12}")
    print(f"  {'-'*12} ", end="")
    for _ in configs:
        print(f"  {'-'*20}", end="")
    print(f"  {'-'*12}")

    for i_deg in cfg.INCLINATION_ANGLES_DEG:
        print(f"  {i_deg:>10.0f} deg ", end="")
        for label, r_bot, _ in configs:
            tr_incl = taper_ratio_inclined(i_deg, r_bot)
            print(f"  {tr_incl:>20.4f}", end="")
        # Correction factor (same for both configs — depends only on inclination)
        i_rad = math.radians(i_deg)
        omega = cfg.OMEGA_SIDEREAL
        corr = math.exp(cfg.RHO_CNT * omega**2 * cfg.R_GEO**2
                        * math.sin(i_rad)**2
                        / (2.0 * cfg.SIGMA_OPERATING))
        print(f"  {corr:>12.6f}")

    # --- Material strength sweep ---
    print(f"\n{sep}")
    print("TAPER RATIO VS MATERIAL STRENGTH")
    print(sep)
    print(f"\n  Safety factor = {cfg.SAFETY_FACTOR:.0f}")
    print(f"  {'sigma_break':>14} {'sigma_op':>12} {'TR (ground)':>14} {'TR (275km)':>14}")
    print(f"  {'(GPa)':>14} {'(GPa)':>12} {'':>14} {'':>14}")
    print(f"  {'-'*56}")

    for s_gpa in cfg.SIGMA_SWEEP_GPA:
        sigma_op = s_gpa * 1e9 / cfg.SAFETY_FACTOR
        tr_g = taper_ratio(cfg.R_BOTTOM_GROUND, sigma=sigma_op)
        tr_s = taper_ratio(cfg.R_BOTTOM_SHORT, sigma=sigma_op)
        # Format large numbers in scientific notation
        if tr_g > 1e6:
            tr_g_str = f"{tr_g:>14.3e}"
        else:
            tr_g_str = f"{tr_g:>14.4f}"
        if tr_s > 1e6:
            tr_s_str = f"{tr_s:>14.3e}"
        else:
            tr_s_str = f"{tr_s:>14.4f}"
        marker = " <-- CNT" if s_gpa == cfg.SIGMA_BREAK_CNT / 1e9 else ""
        print(f"  {s_gpa:>14.0f} {s_gpa/cfg.SAFETY_FACTOR:>12.1f}"
              f" {tr_g_str} {tr_s_str}{marker}")

    print(f"\n{sep}\n")


# =============================================================================
# CSV EXPORT
# =============================================================================

def export_csv():
    """Export taper profile data as CSV for both configurations."""
    configs = [
        ("ground", cfg.R_BOTTOM_GROUND),
        (f"{cfg.ALT_BOTTOM_SHORT/1000:.0f}km", cfg.R_BOTTOM_SHORT),
    ]

    A_ref = 1.0  # 1 m^2 reference

    for label, r_bot in configs:
        r_arr, A_arr = taper_profile(r_bot, cfg.R_GEO, A_ref)
        r_sig, sigma_arr = stress_profile(r_bot, cfg.R_GEO)

        filename = f"taper_profile_{label}.csv"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'altitude_km',
                'radius_m',
                'A_ratio (A/A_bottom)',
                'stress_MPa (uniform cross-section)',
            ])
            for i in range(len(r_arr)):
                alt_km = (r_arr[i] - cfg.R_E) / 1000.0
                A_ratio = A_arr[i] / A_ref
                sigma_mpa = sigma_arr[i] / 1e6
                writer.writerow([
                    f"{alt_km:.3f}",
                    f"{r_arr[i]:.1f}",
                    f"{A_ratio:.6e}",
                    f"{sigma_mpa:.3f}",
                ])

        print(f"  Saved CSV: {filename}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main entry point."""
    sys.stdout.reconfigure(encoding='utf-8')

    # Parse CLI arguments
    show_graphs = []

    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--save":
                cfg.SAVE_GRAPHS = True
            elif arg == "--show":
                cfg.SAVE_GRAPHS = False
            else:
                show_graphs.append(arg)

    # Create output directory if saving
    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    # Print results to console
    print_results()

    # Export CSV
    export_csv()

    # Generate graphs
    if show_graphs:
        print("Generating graphs...")
        if "stress_profile" in show_graphs or "all" in show_graphs:
            plot_stress_profile(show_graphs)
        if "taper_profile" in show_graphs or "all" in show_graphs:
            plot_taper_profile(show_graphs)
        if "taper_vs_strength" in show_graphs or "all" in show_graphs:
            plot_taper_vs_strength(show_graphs)
        print("Done.")


if __name__ == "__main__":
    main()
