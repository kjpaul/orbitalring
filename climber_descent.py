#!/usr/bin/env python3
"""
Section 3: Climber Descent from GEO
====================================

Computes forces, braking power, energy dissipation, and descent timelines
for a climber pod descending from geosynchronous orbit to the ground along
a space elevator tether.

The climber must brake throughout the descent: above GEO it brakes against
centrifugal pull outward, below GEO it brakes against gravitational pull
downward.  All braking energy can be regenerated as electricity.

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python climber_descent.py                    Print results only
    python climber_descent.py all                Generate all graphs
    python climber_descent.py braking_power      Generate braking power graph
    python climber_descent.py --show             Display interactively
    python climber_descent.py --speed=100        Override descent speed (m/s)
    python climber_descent.py --mass=50000       Override pod mass (kg)
"""

import sys
import os
import csv
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import tether_config as cfg

# === PHYSICAL CONSTANTS (from config) ===

GM = cfg.GM
R_E = cfg.R_E
R_GEO = cfg.R_GEO
OMEGA = cfg.OMEGA_SIDEREAL
M_POD_DEFAULT = cfg.M_POD_DEFAULT
V_DESCENT_DEFAULT = cfg.V_DESCENT_DEFAULT
N_POINTS_DEFAULT = cfg.N_POINTS_TETHER


# === FORCE AND POWER FUNCTIONS ===

def f_net(r, m_pod=None):
    """Net force on pod at radius r (positive = toward Earth).

    F_net(r) = m_pod * (GM/r^2 - omega^2 * r)

    Below GEO: gravity > centrifugal, F_net > 0 (downward)
    Above GEO: centrifugal > gravity, F_net < 0 (outward)
    At GEO: F_net = 0 (weightless)
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    return m_pod * (GM / r**2 - OMEGA**2 * r)


def braking_power(r, v_descent, m_pod=None):
    """Braking power required at radius r for given descent speed.

    P_brake = |F_net(r)| * v_descent

    The climber always brakes regardless of direction of net force.
    """
    return abs(f_net(r, m_pod)) * v_descent


def effective_potential(r):
    """Effective potential per unit mass in the rotating frame.

    Phi_eff(r) = -GM/r - omega^2 * r^2 / 2
    """
    return -GM / r - OMEGA**2 * r**2 / 2.0


def total_energy(r_top, r_bottom, m_pod=None):
    """Total braking energy for descent (analytical).

    The climber must dissipate energy equal to the change in effective
    potential.  Above GEO the pod moves inward against centrifugal pull;
    below GEO it moves inward with gravity.  The total braking energy
    is the integral of |F_net| dr over the full path.

    For a descent from r_top to r_bottom (both below or above GEO, or
    spanning GEO), we split at GEO and sum absolute work in each segment.
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT

    phi_top = effective_potential(r_top)
    phi_bottom = effective_potential(r_bottom)

    # If the path crosses GEO, split into two segments
    if r_bottom <= R_GEO <= r_top:
        phi_geo = effective_potential(R_GEO)
        # Above GEO segment: pod descends from r_top to R_GEO
        # Centrifugal exceeds gravity, pod loses centrifugal PE
        e_above = m_pod * abs(phi_geo - phi_top)
        # Below GEO segment: pod descends from R_GEO to r_bottom
        # Gravity exceeds centrifugal, pod loses gravitational PE
        e_below = m_pod * abs(phi_bottom - phi_geo)
        return e_above + e_below
    else:
        return m_pod * abs(phi_bottom - phi_top)


# === DESCENT PROFILE FUNCTIONS ===

def descent_constant_speed(r_top, r_bottom, v_descent, m_pod=None, n_points=None):
    """Compute descent profile at constant speed.

    Returns dict with arrays:
        r              - radius (m)
        alt            - altitude above surface (m)
        time           - elapsed time (s)
        f_net          - net force (N), positive = toward Earth
        p_brake        - braking power (W)
        energy_cumulative - cumulative energy dissipated (J)
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    if n_points is None:
        n_points = N_POINTS_DEFAULT

    r_arr = np.linspace(r_top, r_bottom, n_points)
    alt_arr = r_arr - R_E
    dr = abs(r_arr[1] - r_arr[0])  # positive step size

    # Time: constant speed, so t = distance / v
    dist_from_start = r_top - r_arr  # cumulative distance descended
    time_arr = dist_from_start / v_descent

    # Force and power at each point
    f_net_arr = np.array([f_net(r, m_pod) for r in r_arr])
    p_brake_arr = np.abs(f_net_arr) * v_descent

    # Cumulative energy: integrate |F_net| dr using trapezoidal rule
    energy_arr = np.zeros(n_points)
    abs_f = np.abs(f_net_arr)
    for i in range(1, n_points):
        energy_arr[i] = energy_arr[i - 1] + 0.5 * (abs_f[i - 1] + abs_f[i]) * dr

    return {
        'r': r_arr,
        'alt': alt_arr,
        'time': time_arr,
        'f_net': f_net_arr,
        'p_brake': p_brake_arr,
        'energy_cumulative': energy_arr,
        'v_descent': v_descent,
        'm_pod': m_pod,
    }


def descent_constant_power(r_top, r_bottom, p_max, m_pod=None, n_points=None):
    """Compute descent profile at constant power.

    The descent speed varies: v(r) = P / |F_net(r)|.
    Near GEO where F_net -> 0, speed is capped at 500 m/s to avoid
    singularity (physically, near-zero force means negligible braking).

    Returns dict with arrays:
        r              - radius (m)
        alt            - altitude above surface (m)
        time           - elapsed time (s)
        v_descent      - descent speed at each point (m/s)
        f_net          - net force (N)
        energy_cumulative - cumulative energy dissipated (J)
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    if n_points is None:
        n_points = N_POINTS_DEFAULT

    v_max_cap = 500.0  # speed cap near GEO (m/s)

    r_arr = np.linspace(r_top, r_bottom, n_points)
    alt_arr = r_arr - R_E
    dr = abs(r_arr[1] - r_arr[0])

    f_net_arr = np.array([f_net(r, m_pod) for r in r_arr])
    abs_f = np.abs(f_net_arr)

    # Speed at each point: v = P / |F|, capped
    v_arr = np.where(abs_f > 1e-3, p_max / abs_f, v_max_cap)
    v_arr = np.minimum(v_arr, v_max_cap)

    # Time: integrate dr / v(r)
    time_arr = np.zeros(n_points)
    for i in range(1, n_points):
        v_avg = 0.5 * (v_arr[i - 1] + v_arr[i])
        time_arr[i] = time_arr[i - 1] + dr / v_avg

    # Cumulative energy
    energy_arr = np.zeros(n_points)
    for i in range(1, n_points):
        energy_arr[i] = energy_arr[i - 1] + 0.5 * (abs_f[i - 1] + abs_f[i]) * dr

    return {
        'r': r_arr,
        'alt': alt_arr,
        'time': time_arr,
        'v_descent': v_arr,
        'f_net': f_net_arr,
        'energy_cumulative': energy_arr,
        'p_max': p_max,
        'm_pod': m_pod,
    }


# === CONSOLE OUTPUT ===

def print_descent_table(results):
    """Print descent profile at key altitudes."""
    key_alts_km = [0, 100, 1000, 5000, 10000, 20000,
                   cfg.ALT_GEO / 1000, 50000, 60000]
    key_alts_m = [a * 1000 for a in key_alts_km]

    r_arr = results['r']
    alt_arr = results['alt']
    f_arr = results['f_net']

    is_const_speed = 'v_descent' in results and not isinstance(results['v_descent'], np.ndarray)

    print(f"\n  {'Alt (km)':>12}  {'Radius (km)':>12}  {'F_net (kN)':>12}"
          f"  {'F/m (N/kg)':>12}  {'P_brake (MW)':>12}")
    print(f"  {'-' * 64}")

    m_pod = results['m_pod']
    for alt_target in key_alts_m:
        # Find closest point
        idx = np.argmin(np.abs(alt_arr - alt_target))
        r_val = r_arr[idx]
        alt_val = alt_arr[idx]
        f_val = f_arr[idx]
        f_per_m = f_val / m_pod

        if is_const_speed:
            v = results['v_descent']
        else:
            v = results['v_descent'][idx]
        p_val = abs(f_val) * v

        label = ""
        if abs(alt_val - cfg.ALT_GEO) < 100e3:
            label = "  <-- GEO"

        print(f"  {alt_val/1000:12.1f}  {r_val/1000:12.1f}  {f_val/1000:12.3f}"
              f"  {f_per_m:12.4f}  {p_val/1e6:12.3f}{label}")


def print_force_profile():
    """Print force profile at key altitudes."""
    key_alts_km = [0, 100, 1000, 5000, 10000, 20000,
                   cfg.ALT_GEO / 1000, 50000, 60000]

    sep = "=" * 72
    print(f"\n{sep}")
    print("FORCE PROFILE: Net force on climber vs altitude")
    print(sep)
    print(f"  Pod mass: {M_POD_DEFAULT:,.0f} kg")
    print(f"  GM = {GM:.6e} m^3/s^2")
    print(f"  omega = {OMEGA:.7e} rad/s")
    print(f"  R_GEO = {R_GEO:,.0f} m  ({cfg.ALT_GEO/1000:,.1f} km altitude)")
    print()

    print(f"  {'Alt (km)':>12}  {'Radius (km)':>14}  {'F_net (kN)':>12}"
          f"  {'F/m (m/s^2)':>12}  {'Direction':>12}")
    print(f"  {'-' * 66}")

    for alt_km in key_alts_km:
        r = R_E + alt_km * 1000
        f = f_net(r)
        f_per_m = f / M_POD_DEFAULT
        if abs(f) < 1.0:
            direction = "~ zero"
        elif f > 0:
            direction = "downward"
        else:
            direction = "outward"

        label = ""
        if abs(alt_km - cfg.ALT_GEO / 1000) < 1:
            label = "  <-- GEO"

        print(f"  {alt_km:12.1f}  {r/1000:14.1f}  {f/1000:12.3f}"
              f"  {f_per_m:12.4f}  {direction:>12}{label}")

    # Verify GEO zero crossing
    f_geo = f_net(R_GEO)
    print(f"\n  Verification: F_net at GEO = {f_geo:.6f} N"
          f" (should be ~0, residual from R_GEO rounding)")


def print_descent_summary(speeds, m_pod):
    """Print descent summary for multiple speeds."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("DESCENT SUMMARY: Constant-speed descents from GEO to ground")
    print(sep)
    print(f"  Pod mass: {m_pod:,.0f} kg")
    print(f"  Descent range: GEO ({cfg.ALT_GEO/1000:,.1f} km) to surface")
    print()

    print(f"  {'Speed':>8}  {'Speed':>10}  {'Time':>10}  {'Time':>8}"
          f"  {'Peak P':>10}  {'Energy':>12}  {'E/kg':>10}")
    print(f"  {'(m/s)':>8}  {'(km/h)':>10}  {'(hours)':>10}  {'(days)':>8}"
          f"  {'(MW)':>10}  {'(GJ)':>12}  {'(MJ/kg)':>10}")
    print(f"  {'-' * 72}")

    for v in speeds:
        res = descent_constant_speed(R_GEO, R_E, v, m_pod)
        t_total = res['time'][-1]
        p_peak = np.max(res['p_brake'])
        e_total = res['energy_cumulative'][-1]

        print(f"  {v:8.1f}  {v*3.6:10.1f}  {t_total/3600:10.1f}  {t_total/86400:8.2f}"
              f"  {p_peak/1e6:10.3f}  {e_total/1e9:12.3f}  {e_total/m_pod/1e6:10.3f}")


def print_one_week_check(m_pod):
    """Print one-week descent feasibility check."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("ONE-WEEK DESCENT CHECK")
    print(sep)

    distance = R_GEO - R_E
    t_week = 7 * 24 * 3600
    v_avg = distance / t_week

    print(f"  Distance GEO to ground:  {distance/1000:,.1f} km")
    print(f"  Time:                    7 days = {t_week:,} s")
    print(f"  Required avg speed:      {v_avg:.2f} m/s = {v_avg*3.6:.1f} km/h")
    print(f"  Config V_DESCENT_WEEK:   {cfg.V_DESCENT_WEEK:.2f} m/s")

    # Constant-speed descent at v_avg
    res = descent_constant_speed(R_GEO, R_E, v_avg, m_pod)
    p_peak = np.max(res['p_brake'])
    e_total = res['energy_cumulative'][-1]

    print(f"\n  At constant {v_avg:.2f} m/s:")
    print(f"    Peak braking power:    {p_peak/1e6:.3f} MW")
    print(f"    Total energy:          {e_total/1e9:.3f} GJ")
    print(f"    Energy per kg:         {e_total/m_pod/1e6:.3f} MJ/kg")


def print_energy_summary(m_pod):
    """Print energy recovery summary."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("ENERGY RECOVERY SUMMARY")
    print(sep)

    # Total descent energy (analytical)
    e_descent = total_energy(R_GEO, R_E, m_pod)

    # Effective potential difference (lift energy)
    phi_geo = effective_potential(R_GEO)
    phi_surface = effective_potential(R_E)
    delta_phi = phi_geo - phi_surface
    e_lift = m_pod * abs(delta_phi)

    # Numerical check
    res = descent_constant_speed(R_GEO, R_E, 50.0, m_pod)
    e_numerical = res['energy_cumulative'][-1]

    print(f"  Pod mass:                 {m_pod:,.0f} kg")
    print(f"  Phi_eff(surface):         {phi_surface/1e6:.4f} MJ/kg")
    print(f"  Phi_eff(GEO):             {phi_geo/1e6:.4f} MJ/kg")
    print(f"  Delta Phi_eff:            {delta_phi/1e6:.4f} MJ/kg")
    print()
    print(f"  Descent energy (analyt):  {e_descent/1e9:.4f} GJ")
    print(f"  Descent energy (numer):   {e_numerical/1e9:.4f} GJ")
    print(f"  Lift energy (GEO->sfc):   {e_lift/1e9:.4f} GJ")
    print(f"  Ratio descent/lift:       {e_descent/e_lift:.6f}")
    print(f"  Energy per kg (descent):  {e_descent/m_pod/1e6:.4f} MJ/kg")
    print(f"  Energy per kg (lift):     {e_lift/m_pod/1e6:.4f} MJ/kg")
    print()

    # Note about energy conservation
    # Above GEO: centrifugal > gravity. Moving inward (descending), the pod
    # loses centrifugal PE and gains gravitational PE. The NET effective PE
    # change from GEO to surface equals the integral of |F_net| dr only if
    # we account for the sign change at GEO.
    if abs(e_descent - e_lift) / e_lift < 0.01:
        print(f"  Descent and lift energies match (within 1%): energy is conserved.")
    else:
        print(f"  Note: descent energy ({e_descent/e_lift:.4f}x lift) differs because")
        print(f"  the climber brakes in both segments (above and below GEO).")
        print(f"  The sum of |work| exceeds |Delta Phi_eff| when the path")
        print(f"  crosses the GEO equilibrium point.")


# === GRAPH FUNCTIONS ===

def _save_or_show(fig, filename, save_mode):
    """Save figure or show interactively."""
    if save_mode:
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR,
                                f"{filename}.{cfg.GRAPH_FORMAT}")
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        print(f"  Saved: {filepath}")
    else:
        plt.show()


def plot_force_profile(m_pod, save_mode=True):
    """01-force_profile: F_net/m vs altitude, showing zero crossing at GEO."""
    n = 2000
    # Range from surface to 60,000 km altitude (beyond GEO)
    r_arr = np.linspace(R_E, R_E + 60000e3, n)
    alt_arr = (r_arr - R_E) / 1000  # km

    f_per_m = np.array([(GM / r**2 - OMEGA**2 * r) for r in r_arr])
    f_total = np.array([f_net(r, m_pod) for r in r_arr])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(cfg.GRAPH_WIDTH_INCHES,
                                                    cfg.GRAPH_HEIGHT_INCHES * 1.5),
                                   sharex=True)

    # Panel 1: Force per unit mass
    ax1.plot(alt_arr, f_per_m, 'b-', linewidth=1.5)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=cfg.ALT_GEO / 1000, color='r', linestyle='--', alpha=0.7,
                label=f'GEO ({cfg.ALT_GEO/1000:,.0f} km)')
    ax1.set_ylabel('Net acceleration (m/s$^2$)', fontsize=12)
    ax1.set_title('Net Force on Climber vs Altitude', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Annotate regions
    ax1.annotate('Gravity dominates\n(brakes against fall)',
                 xy=(10000, f_per_m[n // 6]), fontsize=9,
                 ha='center', color='#0072B2')
    ax1.annotate('Centrifugal dominates\n(brakes against outward pull)',
                 xy=(50000, f_per_m[-n // 6]), fontsize=9,
                 ha='center', color='#CC0000')

    # Panel 2: Total force for pod mass
    ax2.plot(alt_arr, f_total / 1000, 'b-', linewidth=1.5)
    ax2.axhline(y=0, color='k', linewidth=0.5)
    ax2.axvline(x=cfg.ALT_GEO / 1000, color='r', linestyle='--', alpha=0.7,
                label=f'GEO ({cfg.ALT_GEO/1000:,.0f} km)')
    ax2.set_xlabel('Altitude (km)', fontsize=12)
    ax2.set_ylabel(f'Net force on {m_pod/1000:.0f}t pod (kN)', fontsize=12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    _save_or_show(fig, "01-force_profile", save_mode)


def plot_braking_power(m_pod, save_mode=True):
    """02-braking_power: P_brake vs altitude for multiple descent speeds."""
    n = 2000
    r_arr = np.linspace(R_E, R_GEO, n)
    alt_arr = (r_arr - R_E) / 1000

    speeds = [20, 40, 60, 80, 100]
    colors = ['#0072B2', '#009E73', '#E69F00', '#CC0000', '#882288']

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    for v, color in zip(speeds, colors):
        p_arr = np.array([braking_power(r, v, m_pod) for r in r_arr])
        ax.plot(alt_arr, p_arr / 1e6, color=color, linewidth=1.5,
                label=f'{v} m/s ({v*3.6:.0f} km/h)')

    ax.set_xlabel('Altitude (km)', fontsize=12)
    ax.set_ylabel('Braking power (MW)', fontsize=12)
    ax.set_title(f'Braking Power vs Altitude (pod mass = {m_pod/1000:.0f} t)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    _save_or_show(fig, "02-braking_power", save_mode)


def plot_descent_timeline(m_pod, save_mode=True):
    """03-descent_timeline: Altitude vs time for various descent profiles."""
    speeds = [20, 40, 60, 80, 100]
    colors = ['#0072B2', '#009E73', '#E69F00', '#CC0000', '#882288']

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    # Constant-speed descents
    for v, color in zip(speeds, colors):
        res = descent_constant_speed(R_GEO, R_E, v, m_pod, n_points=2000)
        ax.plot(res['time'] / 3600, res['alt'] / 1000, color=color,
                linewidth=1.5, label=f'Const {v} m/s')

    # Constant-power descent: use power that gives ~60 m/s at surface
    f_surface = abs(f_net(R_E, m_pod))
    p_const = f_surface * 60.0  # power matched to 60 m/s at surface
    res_cp = descent_constant_power(R_GEO, R_E, p_const, m_pod, n_points=2000)
    ax.plot(res_cp['time'] / 3600, res_cp['alt'] / 1000, 'k--', linewidth=2,
            label=f'Const power {p_const/1e6:.1f} MW')

    ax.axhline(y=cfg.ALT_GEO / 1000, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Descent Timeline: GEO to Ground', fontsize=14)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    _save_or_show(fig, "03-descent_timeline", save_mode)


def plot_energy_profile(m_pod, save_mode=True):
    """04-energy_profile: Cumulative energy dissipated vs altitude."""
    v_ref = 60.0  # reference speed
    res = descent_constant_speed(R_GEO, R_E, v_ref, m_pod, n_points=5000)

    # Analytical lift energy from surface upward
    alt_arr = res['alt']
    e_lift_arr = np.array([
        m_pod * abs(effective_potential(R_E + a) - effective_potential(R_E))
        for a in alt_arr
    ])

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    ax.plot(alt_arr / 1000, res['energy_cumulative'] / 1e9, 'b-', linewidth=2,
            label='Cumulative braking energy (descent)')
    ax.plot(alt_arr / 1000, e_lift_arr / 1e9, 'r--', linewidth=1.5,
            label='Lift energy from surface (reference)')

    # Mark GEO
    ax.axvline(x=cfg.ALT_GEO / 1000, color='gray', linestyle=':', alpha=0.5,
               label=f'GEO ({cfg.ALT_GEO/1000:,.0f} km)')

    # Mark total
    e_total_val = res['energy_cumulative'][-1]
    ax.annotate(f'Total: {e_total_val/1e9:.2f} GJ',
                xy=(alt_arr[0] / 1000, e_total_val / 1e9),
                fontsize=10, ha='left', va='bottom',
                xytext=(5000, e_total_val / 1e9 * 0.85),
                arrowprops=dict(arrowstyle='->', color='blue'))

    ax.set_xlabel('Altitude (km)', fontsize=12)
    ax.set_ylabel('Cumulative energy (GJ)', fontsize=12)
    ax.set_title(f'Energy Profile: Descent from GEO (pod = {m_pod/1000:.0f} t)',
                 fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.invert_xaxis()  # GEO on left, ground on right

    fig.tight_layout()
    _save_or_show(fig, "04-energy_profile", save_mode)


# === CSV EXPORT ===

def export_csv(results, filename="climber_descent_profile.csv"):
    """Save descent profile data to CSV."""
    filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    is_const_power = isinstance(results.get('v_descent'), np.ndarray)

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        if is_const_power:
            writer.writerow([
                'radius (m)', 'altitude (m)', 'time (s)',
                'v_descent (m/s)', 'f_net (N)', 'energy_cumulative (J)'
            ])
            for i in range(len(results['r'])):
                writer.writerow([
                    f"{results['r'][i]:.2f}",
                    f"{results['alt'][i]:.2f}",
                    f"{results['time'][i]:.2f}",
                    f"{results['v_descent'][i]:.4f}",
                    f"{results['f_net'][i]:.4f}",
                    f"{results['energy_cumulative'][i]:.2f}",
                ])
        else:
            writer.writerow([
                'radius (m)', 'altitude (m)', 'time (s)',
                'f_net (N)', 'p_brake (W)', 'energy_cumulative (J)'
            ])
            for i in range(len(results['r'])):
                writer.writerow([
                    f"{results['r'][i]:.2f}",
                    f"{results['alt'][i]:.2f}",
                    f"{results['time'][i]:.2f}",
                    f"{results['f_net'][i]:.4f}",
                    f"{results['p_brake'][i]:.4f}",
                    f"{results['energy_cumulative'][i]:.2f}",
                ])

    print(f"  Saved CSV: {filepath}")


# === MAIN ===

def main():
    """Main entry point."""
    sys.stdout.reconfigure(encoding='utf-8')

    # Defaults
    m_pod = M_POD_DEFAULT
    v_descent = V_DESCENT_DEFAULT
    save_mode = True
    graph_keywords = []

    # === ARGUMENT PARSING ===
    for arg in sys.argv[1:]:
        if arg in ('--help', '-h'):
            print(__doc__)
            return
        elif arg == '--save':
            save_mode = True
        elif arg == '--show':
            save_mode = False
            matplotlib.use('TkAgg')
        elif arg.startswith('--speed='):
            v_descent = float(arg.split('=')[1])
        elif arg.startswith('--mass='):
            m_pod = float(arg.split('=')[1])
        else:
            graph_keywords.append(arg.lower())

    do_all = 'all' in graph_keywords

    # === CONSOLE OUTPUT ===
    sep = "=" * 72
    print(f"\n{sep}")
    print("SECTION 3: CLIMBER DESCENT FROM GEO")
    print(f"Reference: Orbital Ring Engineering by Paul G de Jong")
    print(sep)
    print(f"  Pod mass:          {m_pod:>12,.0f} kg ({m_pod/1000:.0f} t)")
    print(f"  Descent speed:     {v_descent:>12.2f} m/s ({v_descent*3.6:.1f} km/h)")
    print(f"  Earth radius:      {R_E:>12,.0f} m")
    print(f"  GEO radius:        {R_GEO:>12,.0f} m")
    print(f"  GEO altitude:      {cfg.ALT_GEO:>12,.0f} m ({cfg.ALT_GEO/1000:,.1f} km)")
    print(f"  omega (sidereal):  {OMEGA:>16.7e} rad/s")

    # === FORCE PROFILE ===
    print_force_profile()

    # === DESCENT TABLE (default speed) ===
    print(f"\n{sep}")
    print(f"DESCENT PROFILE: Constant speed = {v_descent:.2f} m/s ({v_descent*3.6:.1f} km/h)")
    print(sep)
    res_default = descent_constant_speed(R_GEO, R_E, v_descent, m_pod)
    print_descent_table(res_default)

    # === MULTI-SPEED SUMMARY ===
    speeds = cfg.DESCENT_SPEED_SWEEP
    print_descent_summary(speeds, m_pod)

    # === ONE-WEEK CHECK ===
    print_one_week_check(m_pod)

    # === ENERGY SUMMARY ===
    print_energy_summary(m_pod)

    # === CONSTANT-POWER DESCENT ===
    f_surface = abs(f_net(R_E, m_pod))
    p_const = f_surface * 60.0
    res_cp = descent_constant_power(R_GEO, R_E, p_const, m_pod)

    print(f"\n{sep}")
    print(f"CONSTANT-POWER DESCENT: P = {p_const/1e6:.2f} MW")
    print(sep)
    print(f"  Total time:        {res_cp['time'][-1]/3600:.1f} hours"
          f" ({res_cp['time'][-1]/86400:.2f} days)")
    print(f"  Speed at surface:  {res_cp['v_descent'][-1]:.2f} m/s")
    print(f"  Speed at GEO-1km:  {res_cp['v_descent'][1]:.2f} m/s (capped near GEO)")
    print(f"  Total energy:      {res_cp['energy_cumulative'][-1]/1e9:.3f} GJ")

    # === CSV EXPORT ===
    export_csv(res_default, "climber_descent_const_speed.csv")
    export_csv(res_cp, "climber_descent_const_power.csv")

    # === GRAPHS ===
    if graph_keywords:
        print(f"\n{sep}")
        print("GENERATING GRAPHS")
        print(sep)

        if do_all or 'force_profile' in graph_keywords:
            plot_force_profile(m_pod, save_mode)

        if do_all or 'braking_power' in graph_keywords:
            plot_braking_power(m_pod, save_mode)

        if do_all or 'descent_timeline' in graph_keywords:
            plot_descent_timeline(m_pod, save_mode)

        if do_all or 'energy_profile' in graph_keywords:
            plot_energy_profile(m_pod, save_mode)

    print(f"\n{'=' * 72}")
    print("Done.")
    print(f"{'=' * 72}")


if __name__ == '__main__':
    main()
