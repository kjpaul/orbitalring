#!/usr/bin/env python3
"""
Rendezvous Velocity Analysis (Section 5)
==========================================

Computes rendezvous delta-v requirements for spacecraft transferring from
the orbital ring (250 km altitude) to various points along the space elevator
tether.  The tether co-rotates with Earth at omega_sidereal, so its velocity
at radius r is simply omega * r.  A Hohmann transfer from the ring arrives
at radius r with a different velocity, and the difference is the rendezvous
delta-v that must be supplied to dock with the tether.

Key results: the rendezvous delta-v is minimised where the Hohmann arrival
velocity crosses the tether velocity curve.  GEO is the unique altitude where
tether velocity equals circular orbital velocity (both = omega * R_GEO), but
the Hohmann arrival from 250 km is sub-circular at GEO, requiring a burn.
The total mission delta-v (departure + rendezvous) varies with target altitude.

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python rendezvous_analysis.py                    Print results only
    python rendezvous_analysis.py all                Generate all graphs
    python rendezvous_analysis.py rendezvous_dv      Generate rendezvous dv graph
    python rendezvous_analysis.py velocity_comparison Generate velocity comparison graph
    python rendezvous_analysis.py direct_match       Generate direct match comparison graph
    python rendezvous_analysis.py --show             Display interactively
    python rendezvous_analysis.py --save             Save graphs to files (default)
    python rendezvous_analysis.py --help, -h         Show this help

Graph keywords:
    rendezvous_dv        Rendezvous and departure dv vs altitude
    velocity_comparison  Tether, Hohmann arrival, and circular velocity vs altitude
    direct_match         Hohmann total dv vs direct velocity-matching dv
    all                  Generate all graphs
"""

import sys
import os
import math
import csv

import matplotlib.pyplot as plt

import tether_config as cfg


# =============================================================================
# TETHER VELOCITY
# =============================================================================

def v_tether(r):
    """Tether velocity at radius r.

    The tether co-rotates with Earth at the sidereal rate, so every point
    on it moves at omega * r in the prograde direction.

    Args:
        r: Orbital radius (m)

    Returns:
        Tether tangential velocity (m/s)
    """
    return cfg.OMEGA_SIDEREAL * r


# =============================================================================
# CIRCULAR ORBITAL VELOCITY
# =============================================================================

def v_circular(r):
    """Circular orbital velocity at radius r.

    From vis-viva with a = r (circular orbit): v = sqrt(GM/r).

    Args:
        r: Orbital radius (m)

    Returns:
        Circular orbital velocity (m/s)
    """
    return math.sqrt(cfg.GM / r)


# =============================================================================
# HOHMANN TRANSFER
# =============================================================================

def hohmann_transfer(r_from, r_to):
    """Compute Hohmann transfer parameters between two circular orbits.

    Uses the vis-viva equation: v = sqrt(GM * (2/r - 1/a)).

    The departure burn occurs at r_from, the arrival at r_to.  If r_to > r_from
    the transfer is outward (r_from is periapsis); if r_to < r_from the transfer
    is inward (r_from is apoapsis).

    Args:
        r_from: Departure orbit radius (m)
        r_to:   Target orbit radius (m)

    Returns:
        dict with keys:
            v_circ_from     - circular velocity at departure orbit (m/s)
            v_depart        - departure velocity on transfer orbit (m/s)
            v_arrive        - arrival velocity on transfer orbit (m/s)
            v_circ_to       - circular velocity at target orbit (m/s)
            delta_v_depart  - magnitude of departure burn (m/s)
            delta_v_arrive  - magnitude of circularisation burn at target (m/s)
            transfer_time   - half-period of transfer ellipse (s)
            semi_major_axis - semi-major axis of transfer ellipse (m)
    """
    a = (r_from + r_to) / 2.0

    v_circ_from = math.sqrt(cfg.GM / r_from)
    v_circ_to = math.sqrt(cfg.GM / r_to)

    # Vis-viva at departure point
    v_depart = math.sqrt(cfg.GM * (2.0 / r_from - 1.0 / a))

    # Vis-viva at arrival point
    v_arrive = math.sqrt(cfg.GM * (2.0 / r_to - 1.0 / a))

    # Delta-v magnitudes
    delta_v_depart = abs(v_depart - v_circ_from)
    delta_v_arrive = abs(v_circ_to - v_arrive)

    # Transfer time = half the orbital period of the transfer ellipse
    transfer_time = math.pi * math.sqrt(a**3 / cfg.GM)

    return {
        'v_circ_from': v_circ_from,
        'v_depart': v_depart,
        'v_arrive': v_arrive,
        'v_circ_to': v_circ_to,
        'delta_v_depart': delta_v_depart,
        'delta_v_arrive': delta_v_arrive,
        'transfer_time': transfer_time,
        'semi_major_axis': a,
    }


# =============================================================================
# RENDEZVOUS DELTA-V (HOHMANN FROM RING)
# =============================================================================

def rendezvous_dv(r_target):
    """Delta-v to rendezvous with tether at r_target via Hohmann from ring.

    The spacecraft departs the orbital ring on a Hohmann transfer to r_target,
    then must match the tether velocity (omega * r_target) at arrival.  The
    rendezvous delta-v is the difference between the Hohmann arrival velocity
    and the tether velocity at that radius.

    Args:
        r_target: Target radius on tether (m)

    Returns:
        dict with keys:
            delta_v_depart  - departure burn from ring circular orbit (m/s)
            delta_v_arrive  - rendezvous burn = |v_tether - v_arrive| (m/s)
            delta_v_total   - sum of departure and rendezvous burns (m/s)
            v_tether        - tether velocity at r_target (m/s)
            v_arrive        - Hohmann arrival velocity at r_target (m/s)
            transfer_time   - transfer time (s)
    """
    ht = hohmann_transfer(cfg.R_ORBIT_RING, r_target)

    v_teth = v_tether(r_target)
    dv_rendezvous = abs(v_teth - ht['v_arrive'])

    return {
        'delta_v_depart': ht['delta_v_depart'],
        'delta_v_arrive': dv_rendezvous,
        'delta_v_total': ht['delta_v_depart'] + dv_rendezvous,
        'v_tether': v_teth,
        'v_arrive': ht['v_arrive'],
        'transfer_time': ht['transfer_time'],
    }


# =============================================================================
# DIRECT VELOCITY-MATCHING
# =============================================================================

def direct_match_dv(r_target):
    """Delta-v for direct velocity-matching launch from ring.

    Instead of a Hohmann transfer followed by a rendezvous burn, the
    spacecraft departs the ring on a trajectory where its velocity at
    r_target exactly equals the tether velocity omega * r.  From energy
    conservation:

        v_depart^2 = v_arrive^2 + 2*GM*(1/r_target - 1/r_ring)
        v_arrive = omega * r_target   (by construction)

    so: v_depart = sqrt(omega^2 * r_target^2 + 2*GM*(1/r_ring - 1/r_target))

    If the expression under the square root is negative, this trajectory
    is not physically realisable (the spacecraft cannot reach r_target
    with enough energy to match tether velocity).

    Args:
        r_target: Target radius on tether (m)

    Returns:
        dict with keys:
            v_depart       - required departure velocity from ring (m/s)
            delta_v_depart - |v_depart - v_circ_ring| (m/s)
            feasible       - True if the trajectory is physically realisable
    """
    r_ring = cfg.R_ORBIT_RING
    omega = cfg.OMEGA_SIDEREAL
    v_circ_ring = math.sqrt(cfg.GM / r_ring)

    v_arrive_sq = (omega * r_target) ** 2
    energy_term = 2.0 * cfg.GM * (1.0 / r_ring - 1.0 / r_target)
    v_depart_sq = v_arrive_sq + energy_term

    if v_depart_sq < 0:
        return {
            'v_depart': float('nan'),
            'delta_v_depart': float('nan'),
            'feasible': False,
        }

    v_depart = math.sqrt(v_depart_sq)
    delta_v = abs(v_depart - v_circ_ring)

    return {
        'v_depart': v_depart,
        'delta_v_depart': delta_v,
        'feasible': True,
    }


# =============================================================================
# ALTITUDE SWEEP
# =============================================================================

def sweep_altitudes(alt_min_km=250, alt_max_km=100000, n_points=500):
    """Compute rendezvous delta-v vs altitude.

    Sweeps from alt_min_km to alt_max_km (above Earth's surface) and computes
    all relevant velocities and delta-v values at each point.

    Args:
        alt_min_km: Minimum altitude (km), default 250 (ring altitude)
        alt_max_km: Maximum altitude (km), default 100,000
        n_points:   Number of sample points, default 500

    Returns:
        dict with arrays (lists):
            alt_km          - altitude (km)
            r               - radius (m)
            v_tether        - tether velocity (m/s)
            v_circ          - circular orbital velocity (m/s)
            v_arrive        - Hohmann arrival velocity (m/s)
            dv_depart       - Hohmann departure burn (m/s)
            dv_rendezvous   - rendezvous burn (m/s)
            dv_total        - total Hohmann delta-v (m/s)
            dv_direct       - direct velocity-matching delta-v (m/s)
            transfer_time_h - transfer time (hours)
    """
    alts = []
    radii = []
    vt_arr = []
    vc_arr = []
    va_arr = []
    dvd_arr = []
    dvr_arr = []
    dvt_arr = []
    dvm_arr = []
    tt_arr = []

    for i in range(n_points):
        alt_km = alt_min_km + (alt_max_km - alt_min_km) * i / (n_points - 1)
        r = cfg.R_E + alt_km * 1000.0

        # Skip if target equals ring radius (Hohmann with zero transfer)
        if abs(r - cfg.R_ORBIT_RING) < 1.0:
            alts.append(alt_km)
            radii.append(r)
            vt_arr.append(v_tether(r))
            vc_arr.append(v_circular(r))
            va_arr.append(v_circular(r))  # at ring, arrival = circular
            dvd_arr.append(0.0)
            dvr_arr.append(abs(v_tether(r) - v_circular(r)))
            dvt_arr.append(abs(v_tether(r) - v_circular(r)))
            dm = direct_match_dv(r)
            dvm_arr.append(dm['delta_v_depart'] if dm['feasible'] else float('nan'))
            tt_arr.append(0.0)
            continue

        rdv = rendezvous_dv(r)
        dm = direct_match_dv(r)

        alts.append(alt_km)
        radii.append(r)
        vt_arr.append(rdv['v_tether'])
        vc_arr.append(v_circular(r))
        va_arr.append(rdv['v_arrive'])
        dvd_arr.append(rdv['delta_v_depart'])
        dvr_arr.append(rdv['delta_v_arrive'])
        dvt_arr.append(rdv['delta_v_total'])
        dvm_arr.append(dm['delta_v_depart'] if dm['feasible'] else float('nan'))
        tt_arr.append(rdv['transfer_time'] / 3600.0)

    return {
        'alt_km': alts,
        'r': radii,
        'v_tether': vt_arr,
        'v_circ': vc_arr,
        'v_arrive': va_arr,
        'dv_depart': dvd_arr,
        'dv_rendezvous': dvr_arr,
        'dv_total': dvt_arr,
        'dv_direct': dvm_arr,
        'transfer_time_h': tt_arr,
    }


# =============================================================================
# GRAPH: RENDEZVOUS DELTA-V
# =============================================================================

def plot_rendezvous_dv(show_graph):
    """Plot rendezvous and departure delta-v vs altitude.

    Top panel: rendezvous delta-v (burn to match tether velocity at arrival),
    showing the minimum near GEO.
    Bottom panel: departure delta-v from the ring.
    """
    data = sweep_altitudes()
    alt_geo_km = (cfg.R_GEO - cfg.R_E) / 1000.0

    fig, (ax1, ax2) = plt.subplots(2, 1,
                                    figsize=(cfg.GRAPH_WIDTH_INCHES,
                                             cfg.GRAPH_HEIGHT_INCHES * 1.4))

    # Top panel: rendezvous delta-v
    ax1.plot(data['alt_km'], data['dv_rendezvous'], 'b-', linewidth=1.5,
             label='Rendezvous $\\Delta v$ (match tether)')
    ax1.axvline(x=alt_geo_km, color='red', linestyle='--', linewidth=1.0,
                alpha=0.7, label=f'GEO ({alt_geo_km:.0f} km)')

    # Find and mark minimum
    min_idx = min(range(len(data['dv_rendezvous'])),
                  key=lambda i: data['dv_rendezvous'][i])
    min_alt = data['alt_km'][min_idx]
    min_dv = data['dv_rendezvous'][min_idx]
    ax1.plot(min_alt, min_dv, 'ro', markersize=8, zorder=5)
    ax1.annotate(f'Min: {min_dv:.1f} m/s\nat {min_alt:.0f} km',
                 xy=(min_alt, min_dv),
                 xytext=(min_alt + 5000, min_dv + 200),
                 arrowprops=dict(arrowstyle='->', color='red'),
                 fontsize=10, color='red')

    ax1.set_ylabel('Rendezvous $\\Delta v$ (m/s)', fontsize=13)
    ax1.set_title('Rendezvous $\\Delta v$ vs Altitude (Hohmann from ring)',
                  fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(data['alt_km'][0], data['alt_km'][-1])

    # Bottom panel: departure delta-v
    ax2.plot(data['alt_km'], data['dv_depart'], 'g-', linewidth=1.5,
             label='Departure $\\Delta v$ from ring')
    ax2.axvline(x=alt_geo_km, color='red', linestyle='--', linewidth=1.0,
                alpha=0.7, label=f'GEO ({alt_geo_km:.0f} km)')

    ax2.set_xlabel('Target Altitude (km)', fontsize=13)
    ax2.set_ylabel('Departure $\\Delta v$ (m/s)', fontsize=13)
    ax2.set_title('Departure $\\Delta v$ from Ring vs Target Altitude',
                  fontsize=14)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(data['alt_km'][0], data['alt_km'][-1])

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"01-rendezvous_dv.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: VELOCITY COMPARISON
# =============================================================================

def plot_velocity_comparison(show_graph):
    """Plot tether, Hohmann arrival, and circular orbital velocity vs altitude.

    Shows how the tether velocity (linear in r) crosses the circular orbital
    velocity (proportional to 1/sqrt(r)) at GEO.  The Hohmann arrival
    velocity from the ring is also shown.
    """
    data = sweep_altitudes()
    alt_geo_km = (cfg.R_GEO - cfg.R_E) / 1000.0

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.plot(data['alt_km'], data['v_tether'], 'r-', linewidth=1.5,
            label='Tether velocity ($\\omega \\times r$)')
    ax.plot(data['alt_km'], data['v_arrive'], 'b-', linewidth=1.5,
            label='Hohmann arrival velocity')
    ax.plot(data['alt_km'], data['v_circ'], 'g--', linewidth=1.5,
            label='Circular orbital velocity')
    ax.axvline(x=alt_geo_km, color='orange', linestyle=':', linewidth=1.0,
               alpha=0.8, label=f'GEO ({alt_geo_km:.0f} km)')

    # Mark GEO intersection
    v_geo = v_tether(cfg.R_GEO)
    ax.plot(alt_geo_km, v_geo, 'ko', markersize=8, zorder=5)
    ax.annotate(f'GEO: {v_geo:.0f} m/s',
                xy=(alt_geo_km, v_geo),
                xytext=(alt_geo_km + 5000, v_geo + 300),
                arrowprops=dict(arrowstyle='->', color='black'),
                fontsize=10)

    ax.set_xlabel('Altitude (km)', fontsize=13)
    ax.set_ylabel('Velocity (m/s)', fontsize=13)
    ax.set_title('Velocity Comparison: Tether vs Orbital vs Hohmann Arrival',
                 fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(data['alt_km'][0], data['alt_km'][-1])

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"02-velocity_comparison.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: DIRECT MATCH COMPARISON
# =============================================================================

def plot_direct_match(show_graph):
    """Compare Hohmann total delta-v vs direct velocity-matching delta-v.

    Shows which approach is more fuel-efficient at each altitude.
    """
    data = sweep_altitudes()
    alt_geo_km = (cfg.R_GEO - cfg.R_E) / 1000.0

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.plot(data['alt_km'], data['dv_total'], 'b-', linewidth=1.5,
            label='Hohmann total $\\Delta v$ (depart + rendezvous)')
    ax.plot(data['alt_km'], data['dv_direct'], 'r--', linewidth=1.5,
            label='Direct velocity-matching $\\Delta v$')
    ax.axvline(x=alt_geo_km, color='green', linestyle=':', linewidth=1.0,
               alpha=0.7, label=f'GEO ({alt_geo_km:.0f} km)')

    ax.set_xlabel('Target Altitude (km)', fontsize=13)
    ax.set_ylabel('Total $\\Delta v$ (m/s)', fontsize=13)
    ax.set_title('Hohmann Transfer vs Direct Velocity-Matching',
                 fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(data['alt_km'][0], data['alt_km'][-1])

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"03-direct_match.{cfg.GRAPH_FORMAT}"
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
    """Print rendezvous analysis results to console."""
    sep = "=" * 120

    print(f"\n{sep}")
    print("  SECTION 5: RENDEZVOUS VELOCITY ANALYSIS")
    print(f"  Reference: \"Orbital Ring Engineering\" by Paul G de Jong")
    print(f"{sep}\n")

    # Ring parameters
    alt_ring_km = (cfg.R_ORBIT_RING - cfg.R_E) / 1000.0
    v_circ_ring = v_circular(cfg.R_ORBIT_RING)
    v_teth_ring = v_tether(cfg.R_ORBIT_RING)
    v_teth_ground = v_tether(cfg.R_E)
    v_teth_geo = v_tether(cfg.R_GEO)
    v_circ_geo = v_circular(cfg.R_GEO)

    print(f"  Orbital ring altitude:          {alt_ring_km:.0f} km")
    print(f"  Orbital ring radius:            {cfg.R_ORBIT_RING:,.0f} m")
    print(f"  Ring circular orbit velocity:   {v_circ_ring:,.1f} m/s")
    print(f"  Tether velocity at ring:        {v_teth_ring:,.1f} m/s")
    print(f"  Tether velocity at ground:      {v_teth_ground:,.1f} m/s")
    print(f"  Tether velocity at GEO:         {v_teth_geo:,.1f} m/s")
    print(f"  Circular orbit velocity at GEO: {v_circ_geo:,.1f} m/s")
    print(f"  Sidereal angular rate:          {cfg.OMEGA_SIDEREAL:.7e} rad/s")
    print()

    # Key altitudes table
    key_alts_km = [250, 500, 1000, 2000, 5000, 10000, 20000,
                   (cfg.R_GEO - cfg.R_E) / 1000.0,
                   50000, 75000, 100000]

    print(f"  {'Altitude':>12s}  {'V_tether':>12s}  {'V_arrive':>12s}"
          f"  {'dv_rendez':>12s}  {'dv_depart':>12s}  {'dv_total':>12s}"
          f"  {'t_transfer':>12s}")
    print(f"  {'(km)':>12s}  {'(m/s)':>12s}  {'(m/s)':>12s}"
          f"  {'(m/s)':>12s}  {'(m/s)':>12s}  {'(m/s)':>12s}"
          f"  {'(hours)':>12s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}"
          f"  {'-'*12}  {'-'*12}  {'-'*12}"
          f"  {'-'*12}")

    # Track minimum rendezvous dv
    min_dv_rendez = float('inf')
    min_dv_alt = 0.0

    for alt_km in key_alts_km:
        r = cfg.R_E + alt_km * 1000.0

        if abs(r - cfg.R_ORBIT_RING) < 1.0:
            # At the ring itself
            vt = v_tether(r)
            vc = v_circular(r)
            dv_r = abs(vt - vc)
            label = " <-- RING"
            print(f"  {alt_km:>12.0f}  {vt:>12.1f}  {vc:>12.1f}"
                  f"  {dv_r:>12.1f}  {0.0:>12.1f}  {dv_r:>12.1f}"
                  f"  {0.0:>12.2f}{label}")
            if dv_r < min_dv_rendez:
                min_dv_rendez = dv_r
                min_dv_alt = alt_km
            continue

        rdv = rendezvous_dv(r)
        is_geo = abs(alt_km - (cfg.R_GEO - cfg.R_E) / 1000.0) < 1.0
        label = " <-- GEO" if is_geo else ""

        print(f"  {alt_km:>12.0f}  {rdv['v_tether']:>12.1f}"
              f"  {rdv['v_arrive']:>12.1f}  {rdv['delta_v_arrive']:>12.1f}"
              f"  {rdv['delta_v_depart']:>12.1f}"
              f"  {rdv['delta_v_total']:>12.1f}"
              f"  {rdv['transfer_time']/3600:>12.2f}{label}")

        if rdv['delta_v_arrive'] < min_dv_rendez:
            min_dv_rendez = rdv['delta_v_arrive']
            min_dv_alt = alt_km

    print()
    print(f"  ** Minimum rendezvous dv: {min_dv_rendez:.1f} m/s"
          f" at altitude {min_dv_alt:.0f} km **")

    # Find precise minimum by fine sweep around the coarse minimum
    fine_lo = max(300, min_dv_alt - 5000)
    fine_hi = min(100000, min_dv_alt + 5000)
    fine_data = sweep_altitudes(fine_lo, fine_hi, 2000)
    fine_min_idx = min(range(len(fine_data['dv_rendezvous'])),
                       key=lambda i: fine_data['dv_rendezvous'][i])
    fine_min_alt = fine_data['alt_km'][fine_min_idx]
    fine_min_dv = fine_data['dv_rendezvous'][fine_min_idx]
    print(f"  ** Fine search minimum:   {fine_min_dv:.2f} m/s"
          f" at altitude {fine_min_alt:.1f} km **")

    # GEO-specific analysis
    alt_geo_km = (cfg.R_GEO - cfg.R_E) / 1000.0
    rdv_geo = rendezvous_dv(cfg.R_GEO)
    print()
    print(f"  At GEO ({alt_geo_km:.0f} km):")
    print(f"    Tether velocity = circular orbital velocity = {v_teth_geo:.1f} m/s")
    print(f"    Hohmann arrival velocity from ring         = {rdv_geo['v_arrive']:.1f} m/s")
    print(f"    Rendezvous dv (circularisation burn)       = {rdv_geo['delta_v_arrive']:.1f} m/s")
    print(f"    Departure dv from ring                     = {rdv_geo['delta_v_depart']:.1f} m/s")
    print(f"    Note: GEO is where tether velocity = circular orbit velocity,")
    print(f"    but the Hohmann arrival from 250 km is slower, requiring a burn.")
    print()

    # Direct velocity-matching comparison
    print(f"  {'':->96}")
    print(f"  DIRECT VELOCITY-MATCHING COMPARISON")
    print(f"  {'':->96}")
    print()
    print(f"  {'Altitude':>12s}  {'dv_Hohmann':>12s}  {'dv_direct':>12s}"
          f"  {'Difference':>12s}  {'Better':>12s}")
    print(f"  {'(km)':>12s}  {'(m/s)':>12s}  {'(m/s)':>12s}"
          f"  {'(m/s)':>12s}  {'':>12s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}"
          f"  {'-'*12}  {'-'*12}")

    for alt_km in key_alts_km:
        r = cfg.R_E + alt_km * 1000.0
        dm = direct_match_dv(r)

        if abs(r - cfg.R_ORBIT_RING) < 1.0:
            vt = v_tether(r)
            vc = v_circular(r)
            dv_h = abs(vt - vc)
        else:
            rdv = rendezvous_dv(r)
            dv_h = rdv['delta_v_total']

        if dm['feasible']:
            dv_d = dm['delta_v_depart']
            diff = dv_h - dv_d
            better = "Direct" if dv_d < dv_h else "Hohmann"
            print(f"  {alt_km:>12.0f}  {dv_h:>12.1f}  {dv_d:>12.1f}"
                  f"  {diff:>+12.1f}  {better:>12s}")
        else:
            print(f"  {alt_km:>12.0f}  {dv_h:>12.1f}  {'N/A':>12s}"
                  f"  {'---':>12s}  {'---':>12s}")

    print(f"\n{sep}\n")


# =============================================================================
# CSV EXPORT
# =============================================================================

def export_csv():
    """Export rendezvous analysis data as CSV."""
    data = sweep_altitudes(alt_min_km=250, alt_max_km=100000, n_points=500)

    filename = "rendezvous_analysis.csv"
    filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'altitude_km',
            'radius_m',
            'v_tether_m_s',
            'v_circular_m_s',
            'v_hohmann_arrive_m_s',
            'dv_depart_m_s',
            'dv_rendezvous_m_s',
            'dv_total_hohmann_m_s',
            'dv_direct_match_m_s',
            'transfer_time_hours',
        ])
        for i in range(len(data['alt_km'])):
            writer.writerow([
                f"{data['alt_km'][i]:.3f}",
                f"{data['r'][i]:.1f}",
                f"{data['v_tether'][i]:.3f}",
                f"{data['v_circ'][i]:.3f}",
                f"{data['v_arrive'][i]:.3f}",
                f"{data['dv_depart'][i]:.3f}",
                f"{data['dv_rendezvous'][i]:.3f}",
                f"{data['dv_total'][i]:.3f}",
                f"{data['dv_direct'][i]:.3f}",
                f"{data['transfer_time_h'][i]:.4f}",
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
        if "rendezvous_dv" in show_graphs or "all" in show_graphs:
            plot_rendezvous_dv(show_graphs)
        if "velocity_comparison" in show_graphs or "all" in show_graphs:
            plot_velocity_comparison(show_graphs)
        if "direct_match" in show_graphs or "all" in show_graphs:
            plot_direct_match(show_graphs)
        print("Done.")


if __name__ == "__main__":
    main()
