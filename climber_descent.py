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
    python climber_descent.py regen              Regenerative braking trade-off
    python climber_descent.py eddy_force         Eddy current F(v) curves
    python climber_descent.py eddy_temperature   Eddy current plate temperature
    python climber_descent.py --show             Display interactively
    python climber_descent.py --speed=100        Override descent speed (m/s)
    python climber_descent.py --mass=50000       Override pod mass (kg)
"""

import sys
import os
import csv
import math

import numpy as np
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


# === REGENERATIVE BRAKING — CONDUCTIVE CORE ===

def regen_braking_at_altitude(h, v_descent, m_pod=None, A_cond=None, V_bus=None):
    """Regenerative braking parameters when pod is at altitude h.

    The pod transmits braking power as DC current through the tether's
    conductive CNT core to the ground station.

    Returns dict with: P_brake, I, R_cond, P_loss, efficiency,
    q_per_meter (I²R heating per meter), T_eq (conductor temperature).
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    if A_cond is None:
        A_cond = cfg.A_CONDUCTOR_DEFAULT
    if V_bus is None:
        V_bus = cfg.V_BUS

    r = R_E + h
    P = abs(f_net(r, m_pod)) * v_descent
    I = P / V_bus if V_bus > 0 else 0.0

    # Conductor resistance from pod to ground station
    rho_e = cfg.RHO_ELEC_CNT
    R_cond = rho_e * h / A_cond if h > 0 else 0.0

    # I²R total loss
    P_loss = I**2 * R_cond

    # Efficiency
    eta = 1.0 - P_loss / P if P > 0 else 1.0

    # I²R heating per meter of conductor (uniform along tether below pod)
    q_per_m = I**2 * rho_e / A_cond

    # Conductor diameter (circular cross-section)
    d_cond = 2.0 * math.sqrt(A_cond / math.pi)

    # Equilibrium temperature: radiative cooling in vacuum
    # q = ε × σ_SB × T⁴ × π × d
    eps = cfg.EPSILON_CNT_IR
    sigma_sb = cfg.SIGMA_SB
    denom = eps * sigma_sb * math.pi * d_cond
    T_eq = (q_per_m / denom) ** 0.25 if q_per_m > 0 and denom > 0 else 0.0

    return {
        'P_brake': P,
        'I': I,
        'R_cond': R_cond,
        'P_loss': P_loss,
        'efficiency': eta,
        'q_per_meter': q_per_m,
        'T_eq': T_eq,
        'd_cond': d_cond,
    }


def regen_profile(v_descent, m_pod=None, A_cond=None, V_bus=None, n_points=500):
    """Compute regenerative braking profile vs altitude.

    Returns dict with arrays for altitude, power, current, loss, efficiency, temperature.
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT

    alts = np.linspace(1000, cfg.ALT_GEO, n_points)  # skip h=0 (R_cond=0)
    results = {
        'alt': alts,
        'P_brake': np.zeros(n_points),
        'I': np.zeros(n_points),
        'P_loss': np.zeros(n_points),
        'efficiency': np.zeros(n_points),
        'q_per_meter': np.zeros(n_points),
        'T_eq': np.zeros(n_points),
    }

    for i, h in enumerate(alts):
        r = regen_braking_at_altitude(h, v_descent, m_pod, A_cond, V_bus)
        for key in ['P_brake', 'I', 'P_loss', 'efficiency', 'q_per_meter', 'T_eq']:
            results[key][i] = r[key]

    return results


def conductor_tradeoff(v_descent, h_ref, m_pod=None, V_bus=None):
    """Conductor cross-section trade-off at a reference altitude.

    Returns list of dicts, one per A_cond value.
    """
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    results = []
    for a_cm2 in cfg.A_CONDUCTOR_SWEEP_CM2:
        A_cond = a_cm2 * 1e-4  # cm² to m²
        r = regen_braking_at_altitude(h_ref, v_descent, m_pod, A_cond, V_bus)
        r['A_cond_cm2'] = a_cm2
        r['A_cond'] = A_cond
        # Mass per meter of conductor
        r['mass_per_m'] = A_cond * cfg.RHO_CNT
        results.append(r)
    return results


# === EDDY CURRENT BRAKING — REACTION PLATE ===

def eddy_critical_speed(d_plate, sigma_plate=None, L_pole=None):
    """Critical speed (maximum braking force) for eddy current brake.

    v_c where skin depth δ = plate thickness d_plate.
    δ = sqrt(2 L_pole / (π μ₀ σ v))
    Setting δ = d: v_c = 2 L_pole / (π μ₀ σ d²)
    """
    if sigma_plate is None:
        sigma_plate = cfg.SIGMA_PLATE_AL
    if L_pole is None:
        L_pole = cfg.L_POLE
    return 2.0 * L_pole / (math.pi * cfg.MU_0 * sigma_plate * d_plate**2)


def eddy_skin_depth(v, sigma_plate=None, L_pole=None):
    """Electromagnetic skin depth at speed v.

    f_eff = v / (2 L_pole)
    δ = sqrt(2 / (μ₀ σ 2π f_eff))
    """
    if sigma_plate is None:
        sigma_plate = cfg.SIGMA_PLATE_AL
    if L_pole is None:
        L_pole = cfg.L_POLE
    f_eff = v / (2.0 * L_pole)
    if f_eff <= 0:
        return float('inf')
    return math.sqrt(2.0 / (cfg.MU_0 * sigma_plate * 2.0 * math.pi * f_eff))


def eddy_force_per_meter(v, d_plate, B_gap=None, w_plate=None,
                         sigma_plate=None, L_pole=None):
    """Eddy current braking force per meter of reaction plate track.

    F(v) = f₀ × (v/v_c) / (1 + (v/v_c)²)

    where f₀ = B₀² × w / (2μ₀) and peak force = f₀/2 at v = v_c.

    Returns force in N per meter of track.
    """
    if B_gap is None:
        B_gap = cfg.B_GAP
    if w_plate is None:
        w_plate = cfg.PLATE_WIDTH
    if sigma_plate is None:
        sigma_plate = cfg.SIGMA_PLATE_AL
    if L_pole is None:
        L_pole = cfg.L_POLE

    v_c = eddy_critical_speed(d_plate, sigma_plate, L_pole)
    f_0 = B_gap**2 * w_plate / (2.0 * cfg.MU_0)

    u = v / v_c if v_c > 0 else 0.0
    return f_0 * u / (1.0 + u**2)


def eddy_total_force(v, d_plate, L_magnet=None, **kwargs):
    """Total eddy current braking force on pod (N).

    Force per meter × magnet array length.
    """
    if L_magnet is None:
        L_magnet = cfg.L_MAGNET_POD
    return eddy_force_per_meter(v, d_plate, **kwargs) * L_magnet


def eddy_plate_temperature(v, d_plate, B_gap=None, w_plate=None,
                           sigma_plate=None, L_pole=None,
                           epsilon=None):
    """Equilibrium temperature of the reaction plate at speed v.

    Power dissipated per meter of track = f(v) × v
    Radiative cooling per meter = ε × σ_SB × T⁴ × (w + 2d)

    Returns T_eq in Kelvin.
    """
    if w_plate is None:
        w_plate = cfg.PLATE_WIDTH
    if epsilon is None:
        epsilon = cfg.EPSILON_CNT_IR

    f = eddy_force_per_meter(v, d_plate, B_gap, w_plate,
                             sigma_plate, L_pole)
    P_dissipated = f * v  # W per meter of track

    # Radiative cooling surface: one face (w) + two edges (2d)
    perimeter = w_plate + 2.0 * d_plate
    denom = epsilon * cfg.SIGMA_SB * perimeter
    if denom <= 0 or P_dissipated <= 0:
        return 0.0
    return (P_dissipated / denom) ** 0.25


def eddy_max_continuous_speed(d_plate, T_max=None, **kwargs):
    """Maximum continuous descent speed for eddy current braking.

    Finds speed where plate equilibrium temperature = T_max.
    Uses bisection search.
    """
    if T_max is None:
        T_max = cfg.T_PLATE_MAX

    # Search from 0.1 to 10000 m/s
    v_lo, v_hi = 0.1, 10000.0
    T_hi = eddy_plate_temperature(v_hi, d_plate, **kwargs)

    # If even max speed is below T_max, return inf
    if T_hi < T_max:
        return float('inf')

    # Check low end
    T_lo = eddy_plate_temperature(v_lo, d_plate, **kwargs)
    if T_lo >= T_max:
        return 0.0

    # Bisection
    for _ in range(100):
        v_mid = (v_lo + v_hi) / 2.0
        T_mid = eddy_plate_temperature(v_mid, d_plate, **kwargs)
        if T_mid < T_max:
            v_lo = v_mid
        else:
            v_hi = v_mid
        if abs(v_hi - v_lo) < 0.01:
            break

    return (v_lo + v_hi) / 2.0


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


def plot_regen_efficiency(m_pod, save_mode=True):
    """05-regen_efficiency: Efficiency and temperature vs conductor area."""
    v_ref = 200.0  # m/s reference speed
    h_ref = 10000e3  # 10,000 km altitude (mid-tether)

    areas_cm2 = cfg.A_CONDUCTOR_SWEEP_CM2
    tradeoff = conductor_tradeoff(v_ref, h_ref, m_pod)

    fig, (ax1, ax2) = plt.subplots(1, 2,
                                    figsize=(cfg.GRAPH_WIDTH_INCHES,
                                             cfg.GRAPH_HEIGHT_INCHES))

    etas = [r['efficiency'] * 100 for r in tradeoff]
    temps = [r['T_eq'] for r in tradeoff]

    ax1.bar(range(len(areas_cm2)), etas, color='#0072B2', alpha=0.8)
    ax1.set_xticks(range(len(areas_cm2)))
    ax1.set_xticklabels([f"{a}" for a in areas_cm2])
    ax1.set_xlabel('Conductor area (cm²)', fontsize=11)
    ax1.set_ylabel('Transmission efficiency (%)', fontsize=11)
    ax1.set_title(f'Efficiency at 10,000 km, {v_ref:.0f} m/s', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_ylim(0, 105)

    ax2.bar(range(len(areas_cm2)), temps, color='#CC0000', alpha=0.8)
    ax2.set_xticks(range(len(areas_cm2)))
    ax2.set_xticklabels([f"{a}" for a in areas_cm2])
    ax2.set_xlabel('Conductor area (cm²)', fontsize=11)
    ax2.set_ylabel('Conductor temperature (K)', fontsize=11)
    ax2.set_title(f'Equilibrium T at surface, {v_ref:.0f} m/s', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')

    fig.suptitle(f'Regenerative Braking Conductor Trade-off ({m_pod/1000:.0f}t pod)',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    _save_or_show(fig, "05-regen_efficiency", save_mode)


def plot_eddy_force(save_mode=True):
    """06-eddy_force: Braking force vs speed for various plate thicknesses."""
    v_arr = np.logspace(-1, 3, 500)  # 0.1 to 1000 m/s

    fig, (ax1, ax2) = plt.subplots(2, 1,
                                    figsize=(cfg.GRAPH_WIDTH_INCHES,
                                             cfg.GRAPH_HEIGHT_INCHES * 1.4))
    colors = ['#0072B2', '#009E73', '#E69F00', '#CC0000', '#882288']

    for d_mm, color in zip([1, 2, 5, 10, 20], colors):
        d = d_mm / 1000.0
        v_c = eddy_critical_speed(d)
        f_arr = np.array([eddy_force_per_meter(v, d) for v in v_arr])
        F_total = np.array([eddy_total_force(v, d) for v in v_arr])

        ax1.loglog(v_arr, f_arr / 1000, color=color, linewidth=1.5,
                   label=f'd={d_mm} mm (v_c={v_c:.1f} m/s)')
        ax2.loglog(v_arr, F_total / 1000, color=color, linewidth=1.5,
                   label=f'd={d_mm} mm')

        # Mark critical speed
        f_peak = eddy_force_per_meter(v_c, d)
        ax1.plot(v_c, f_peak / 1000, 'o', color=color, markersize=6)

    ax1.set_ylabel('Force per meter of track (kN/m)', fontsize=11)
    ax1.set_title('Eddy Current Braking Force vs Speed', fontsize=14)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, which='both')

    ax2.set_xlabel('Pod speed (m/s)', fontsize=11)
    ax2.set_ylabel(f'Total pod force (kN), L_mag={cfg.L_MAGNET_POD:.0f}m',
                   fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, which='both')

    # Add reference line for gravity at surface
    f_gravity = abs(f_net(R_E)) / 1000
    ax2.axhline(y=f_gravity, color='gray', linestyle='--', alpha=0.5,
                label=f'Surface gravity ({f_gravity:.0f} kN)')
    ax2.legend(fontsize=9)

    fig.tight_layout()
    _save_or_show(fig, "06-eddy_force", save_mode)


def plot_eddy_temperature(save_mode=True):
    """07-eddy_temperature: Plate equilibrium temperature vs speed."""
    v_arr = np.logspace(-1, 3, 500)

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))
    colors = ['#0072B2', '#009E73', '#E69F00', '#CC0000', '#882288']

    for d_mm, color in zip([1, 2, 5, 10, 20], colors):
        d = d_mm / 1000.0
        T_arr = np.array([eddy_plate_temperature(v, d) for v in v_arr])
        ax.semilogx(v_arr, T_arr, color=color, linewidth=1.5,
                    label=f'd={d_mm} mm')

    # Temperature limits
    ax.axhline(y=cfg.T_PLATE_MAX, color='red', linestyle='--', alpha=0.7,
               label=f'T_max = {cfg.T_PLATE_MAX:.0f} K')
    ax.axhline(y=933, color='darkred', linestyle=':', alpha=0.5,
               label='Al melting (933 K)')

    ax.set_xlabel('Pod speed (m/s)', fontsize=12)
    ax.set_ylabel('Plate equilibrium temperature (K)', fontsize=12)
    ax.set_title('Eddy Current Brake — Plate Temperature vs Speed', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_ylim(0, 1200)

    fig.tight_layout()
    _save_or_show(fig, "07-eddy_temperature", save_mode)


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

    # === REGENERATIVE BRAKING ===
    print(f"\n{sep}")
    print("FAST DESCENT WITH REGENERATIVE BRAKING")
    print(sep)
    print(f"  Conductor material:        CNT")
    print(f"  Electrical conductivity:   {cfg.SIGMA_ELEC_CNT/1e6:.0f} MS/m")
    print(f"  Resistivity:               {cfg.RHO_ELEC_CNT:.3e} Ω·m")
    print(f"  Bus voltage:               {cfg.V_BUS/1e3:.0f} kV DC")
    print(f"  Default conductor area:    {cfg.A_CONDUCTOR_DEFAULT*1e4:.2f} cm²")
    print(f"  Default conductor dia:     {cfg.D_CONDUCTOR_DEFAULT*1000:.1f} mm")
    print(f"  IR emissivity:             {cfg.EPSILON_CNT_IR}")
    print()

    print(f"  {'Speed':>8}  {'P_brake':>10}  {'I_peak':>10}  {'I²R at':>12}"
          f"  {'Efficiency':>12}  {'q/m':>10}  {'T_cond':>8}")
    print(f"  {'(m/s)':>8}  {'(MW)':>10}  {'(A)':>10}  {'10kkm (MW)':>12}"
          f"  {'at 10kkm':>12}  {'(W/m)':>10}  {'(K)':>8}")
    print(f"  {'-' * 76}")

    for v in cfg.FAST_DESCENT_SPEEDS:
        # Peak values at surface
        r_surface = regen_braking_at_altitude(1000, v, m_pod)
        # Efficiency at 10,000 km (representative mid-tether)
        r_mid = regen_braking_at_altitude(10000e3, v, m_pod)

        print(f"  {v:8.0f}  {r_surface['P_brake']/1e6:10.2f}"
              f"  {r_surface['I']:10.1f}"
              f"  {r_mid['P_loss']/1e6:12.3f}"
              f"  {r_mid['efficiency']*100:11.1f}%"
              f"  {r_surface['q_per_meter']:10.2f}"
              f"  {r_surface['T_eq']:8.0f}")

    # Conductor cross-section trade-off
    print(f"\n{sep}")
    print("CONDUCTOR CROSS-SECTION TRADE-OFF")
    print(sep)
    v_ref = 200.0
    h_ref = 10000e3
    print(f"  Reference: v = {v_ref:.0f} m/s, pod at 10,000 km altitude")
    print()

    tradeoff = conductor_tradeoff(v_ref, h_ref, m_pod)
    print(f"  {'A_cond':>8}  {'Diameter':>10}  {'Mass/m':>10}"
          f"  {'Efficiency':>12}  {'I²R loss':>10}  {'T_cond':>8}")
    print(f"  {'(cm²)':>8}  {'(mm)':>10}  {'(kg/m)':>10}"
          f"  {'(%)':>12}  {'(MW)':>10}  {'(K)':>8}")
    print(f"  {'-' * 64}")

    for r in tradeoff:
        print(f"  {r['A_cond_cm2']:8.0f}  {r['d_cond']*1000:10.1f}"
              f"  {r['mass_per_m']:10.2f}"
              f"  {r['efficiency']*100:11.1f}%"
              f"  {r['P_loss']/1e6:10.3f}"
              f"  {r['T_eq']:8.0f}")

    # === EDDY CURRENT BRAKING ===
    print(f"\n{sep}")
    print("EDDY CURRENT BRAKING ANALYSIS")
    print(sep)
    print(f"  Plate material:            Aluminum")
    print(f"  Plate conductivity:        {cfg.SIGMA_PLATE_AL/1e6:.1f} MS/m")
    print(f"  Magnetic field B_gap:      {cfg.B_GAP:.1f} T")
    print(f"  Pole pitch L_pole:         {cfg.L_POLE*1000:.0f} mm")
    print(f"  Plate width:               {cfg.PLATE_WIDTH*1000:.0f} mm")
    print(f"  Pod magnet length:         {cfg.L_MAGNET_POD:.0f} m")
    print(f"  Max plate temperature:     {cfg.T_PLATE_MAX:.0f} K")
    print()

    f_0 = cfg.B_GAP**2 * cfg.PLATE_WIDTH / (2.0 * cfg.MU_0)
    print(f"  Characteristic force f₀:   {f_0:.0f} N/m")
    print(f"  Peak force (f₀/2):         {f_0/2:.0f} N/m")
    print()

    print(f"  {'d_plate':>8}  {'v_crit':>10}  {'δ at 100':>10}"
          f"  {'F_peak/m':>10}  {'F_pod':>10}  {'T at v_c':>8}"
          f"  {'v_max':>10}")
    print(f"  {'(mm)':>8}  {'(m/s)':>10}  {'m/s (mm)':>10}"
          f"  {'(N/m)':>10}  {'(kN)':>10}  {'(K)':>8}"
          f"  {'(m/s)':>10}")
    print(f"  {'-' * 72}")

    for d_mm in [1, 2, 5, 10, 20]:
        d = d_mm / 1000.0
        v_c = eddy_critical_speed(d)
        delta_100 = eddy_skin_depth(100.0)
        f_peak = eddy_force_per_meter(v_c, d)
        F_pod = eddy_total_force(v_c, d)
        T_vc = eddy_plate_temperature(v_c, d)
        v_max = eddy_max_continuous_speed(d)
        v_max_str = f"{v_max:10.1f}" if v_max < 9999 else "      >10k"

        print(f"  {d_mm:8.0f}  {v_c:10.2f}  {delta_100*1000:10.2f}"
              f"  {f_peak:10.0f}  {F_pod/1000:10.1f}  {T_vc:8.0f}"
              f"  {v_max_str}")

    # Radiative cooling capacity
    print(f"\n  Radiative cooling capacity at various temperatures:")
    print(f"  (plate width {cfg.PLATE_WIDTH*1000:.0f} mm, thickness 10 mm)")
    print()
    d_ref = 0.010
    perim = cfg.PLATE_WIDTH + 2 * d_ref
    for T in [300, 500, 700, 800, 933, 1200]:
        P_rad = cfg.EPSILON_CNT_IR * cfg.SIGMA_SB * T**4 * perim
        print(f"    T = {T:6.0f} K:  P_rad = {P_rad:8.1f} W/m")

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

        if do_all or 'regen' in graph_keywords:
            plot_regen_efficiency(m_pod, save_mode)

        if do_all or 'eddy_force' in graph_keywords:
            plot_eddy_force(save_mode)

        if do_all or 'eddy_temperature' in graph_keywords:
            plot_eddy_temperature(save_mode)

    print(f"\n{'=' * 72}")
    print("Done.")
    print(f"{'=' * 72}")


if __name__ == '__main__':
    main()
