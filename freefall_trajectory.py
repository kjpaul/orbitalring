#!/usr/bin/env python3
"""
Free-fall Trajectory from 275 km Release
=========================================

Simulates the descent of a pod released from the orbital ring at 275 km
altitude.  The pod starts at rest relative to Earth's surface (co-rotating
at sidereal rate), so in the inertial frame it has tangential velocity
v_0 = omega * r_0.  After release it follows a ballistic trajectory,
initially in vacuum, then through the atmosphere where aerodynamic drag
and heating become significant.

The 2D polar-coordinate treatment preserves angular-momentum conservation
during the vacuum phase and correctly captures the small but non-zero
ground-track drift.

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python freefall_trajectory.py                Print results only
    python freefall_trajectory.py all            Generate all graphs
    python freefall_trajectory.py velocity       Generate velocity profile
    python freefall_trajectory.py --show         Display interactively
    python freefall_trajectory.py --beta=1000    Override ballistic coefficient
    python freefall_trajectory.py --mass=50000   Override pod mass (kg)
"""

import sys
import os
import math
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import tether_config as cfg

# === CONSTANTS ===

GM = cfg.GM
R_E = cfg.R_E
OMEGA = cfg.OMEGA_SIDEREAL
G0 = cfg.G_SURFACE

ALT_RELEASE = 275_000.0                   # Release altitude (m)
R_RELEASE = R_E + ALT_RELEASE              # Release radius (m)
V_TANGENTIAL_0 = OMEGA * R_RELEASE         # Inertial tangential velocity (m/s)

# Pod defaults
BETA_DEFAULT = cfg.BALLISTIC_COEFF         # kg/m^2
M_POD_DEFAULT = cfg.M_POD_DEFAULT          # kg

# Sweep values
BETA_SWEEP = cfg.BETA_SWEEP                # [100, 200, 500, 1000, 2000]

# Graph numbering prefix
GRAPH_PREFIX = ""

# Key altitudes for the summary table (m)
TABLE_ALTITUDES = [275e3, 200e3, 150e3, 100e3, 80e3, 60e3, 50e3, 40e3, 30e3, 20e3, 10e3, 0.0]


# === PHYSICS HELPERS ===

def vacuum_velocity(r):
    """Vacuum velocity at radius r from energy conservation.

    v(r) = sqrt(v0^2 + 2*GM*(1/r - 1/r0))
    """
    return math.sqrt(V_TANGENTIAL_0**2 + 2.0 * GM * (1.0 / r - 1.0 / R_RELEASE))


def atmosphere_density(alt):
    """Return atmospheric density at altitude alt (m)."""
    _, _, rho = cfg.atmosphere(alt)
    return rho


def local_sound_speed(alt):
    """Return speed of sound at altitude alt (m)."""
    T, _, _ = cfg.atmosphere(alt)
    return cfg.speed_of_sound(T)


# === RK4 INTEGRATOR ===

def rk4_step(state, dt, beta):
    """Advance state = [r, theta, v_r, v_theta] by one RK4 step.

    Equations of motion in polar coordinates (inertial frame):
        dr/dt     = v_r
        dtheta/dt = v_theta / r
        dv_r/dt   = -GM/r^2 + v_theta^2/r - drag_r/m
        dv_theta/dt = -v_r*v_theta/r - drag_theta/m

    Drag is computed using velocity relative to the co-rotating atmosphere.
    """
    def derivs(s):
        r, theta, vr, vt = s
        alt = r - R_E
        if alt < 0:
            alt = 0.0

        # Atmosphere-relative velocity
        vr_rel = vr
        vt_rel = vt - OMEGA * r
        v_rel = math.sqrt(vr_rel**2 + vt_rel**2)

        # Drag acceleration magnitude: F_drag/m = 0.5 * rho * v_rel^2 / beta
        rho = atmosphere_density(alt)
        drag_accel = 0.5 * rho * v_rel / beta if v_rel > 1e-12 else 0.0
        # drag_accel * v_rel_component gives component deceleration
        # Direction: oppose v_rel
        drag_r = drag_accel * vr_rel
        drag_t = drag_accel * vt_rel

        dr = vr
        dtheta = vt / r
        dvr = -GM / (r * r) + vt * vt / r - drag_r
        dvt = -vr * vt / r - drag_t

        return np.array([dr, dtheta, dvr, dvt])

    s = np.array(state, dtype=float)
    k1 = derivs(s)
    k2 = derivs(s + 0.5 * dt * k1)
    k3 = derivs(s + 0.5 * dt * k2)
    k4 = derivs(s + dt * k3)
    return s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


# === SIMULATION ===

def simulate_descent(beta=None, m_pod=None, dt=None):
    """Simulate full descent from 275 km to ground.

    Parameters
    ----------
    beta : float
        Ballistic coefficient m/(C_D * A) in kg/m^2.
    m_pod : float
        Pod mass in kg (used only for reporting; dynamics depend on beta).
    dt : float
        Integration time step in seconds.

    Returns
    -------
    dict with time-series numpy arrays:
        t, r, theta, v_r, v_theta, v_total, alt, mach,
        q_dot, q_total, rho_local, g_load, v_rel
    """
    if beta is None:
        beta = BETA_DEFAULT
    if m_pod is None:
        m_pod = M_POD_DEFAULT
    if dt is None:
        dt = cfg.DT_TRAJECTORY

    # Initial conditions (inertial frame)
    r0 = R_RELEASE
    theta0 = 0.0
    vr0 = 0.0           # No radial velocity at release
    vt0 = V_TANGENTIAL_0  # Tangential = omega * r0

    state = [r0, theta0, vr0, vt0]

    # Pre-allocate storage
    max_steps = cfg.N_POINTS_TRAJECTORY
    t_arr = np.zeros(max_steps)
    r_arr = np.zeros(max_steps)
    theta_arr = np.zeros(max_steps)
    vr_arr = np.zeros(max_steps)
    vt_arr = np.zeros(max_steps)
    vtot_arr = np.zeros(max_steps)
    alt_arr = np.zeros(max_steps)
    mach_arr = np.zeros(max_steps)
    qdot_arr = np.zeros(max_steps)
    qtot_arr = np.zeros(max_steps)
    rho_arr = np.zeros(max_steps)
    gload_arr = np.zeros(max_steps)
    vrel_arr = np.zeros(max_steps)

    # Store initial values
    alt0 = r0 - R_E
    rho0 = atmosphere_density(alt0)
    a0 = local_sound_speed(alt0)
    vrel0 = math.sqrt(vr0**2 + (vt0 - OMEGA * r0)**2)

    t_arr[0] = 0.0
    r_arr[0] = r0
    theta_arr[0] = theta0
    vr_arr[0] = vr0
    vt_arr[0] = vt0
    vtot_arr[0] = math.sqrt(vr0**2 + vt0**2)
    alt_arr[0] = alt0
    mach_arr[0] = vrel0 / a0 if a0 > 0 else 0.0
    qdot_arr[0] = 0.0
    qtot_arr[0] = 0.0
    rho_arr[0] = rho0
    gload_arr[0] = 0.0
    vrel_arr[0] = vrel0

    t = 0.0
    q_total = 0.0
    n_steps = 1

    for i in range(1, max_steps):
        state_new = rk4_step(state, dt, beta)
        t += dt

        r, theta, vr, vt = state_new
        alt = r - R_E
        if alt < 0:
            alt = 0.0

        # Atmosphere-relative velocity
        vr_rel = vr
        vt_rel = vt - OMEGA * r
        v_rel = math.sqrt(vr_rel**2 + vt_rel**2)
        v_total = math.sqrt(vr**2 + vt**2)

        # Atmospheric properties
        rho = atmosphere_density(alt)
        T_atm, _, _ = cfg.atmosphere(alt)
        a_sound = cfg.speed_of_sound(T_atm)

        # Mach number (relative to atmosphere)
        mach = v_rel / a_sound if a_sound > 0 else 0.0

        # Stagnation-point heat flux: q = k * sqrt(rho / r_nose) * v_rel^3
        q_dot = cfg.K_SUTTON_GRAVES * math.sqrt(rho / cfg.R_NOSE) * v_rel**3
        q_total += q_dot * dt

        # g-load: drag deceleration / g0
        drag_accel = 0.5 * rho * v_rel**2 / beta
        g_load = drag_accel / G0

        # Store
        t_arr[i] = t
        r_arr[i] = r
        theta_arr[i] = theta
        vr_arr[i] = vr
        vt_arr[i] = vt
        vtot_arr[i] = v_total
        alt_arr[i] = alt
        mach_arr[i] = mach
        qdot_arr[i] = q_dot
        qtot_arr[i] = q_total
        rho_arr[i] = rho
        gload_arr[i] = g_load
        vrel_arr[i] = v_rel

        n_steps = i + 1
        state = state_new

        # Stop when altitude reaches zero
        if alt <= 0:
            break

    # Trim arrays
    sl = slice(0, n_steps)
    return {
        'beta': beta,
        'm_pod': m_pod,
        't': t_arr[sl],
        'r': r_arr[sl],
        'theta': theta_arr[sl],
        'v_r': vr_arr[sl],
        'v_theta': vt_arr[sl],
        'v_total': vtot_arr[sl],
        'alt': alt_arr[sl],
        'mach': mach_arr[sl],
        'q_dot': qdot_arr[sl],
        'q_total': qtot_arr[sl],
        'rho': rho_arr[sl],
        'g_load': gload_arr[sl],
        'v_rel': vrel_arr[sl],
    }


# === VACUUM REFERENCE CURVE ===

def vacuum_reference_curve():
    """Return (alt_array, v_array) for vacuum free-fall from release."""
    alts = np.linspace(0, ALT_RELEASE, 2000)
    vs = np.array([vacuum_velocity(R_E + a) for a in alts])
    return alts, vs


# === CONSOLE OUTPUT ===

def print_validation():
    """Print physics validation checks."""
    print("=" * 72)
    print("VALIDATION CHECKS")
    print("=" * 72)

    v_surface_vac = vacuum_velocity(R_E)
    rho_0 = atmosphere_density(0)
    rho_50 = atmosphere_density(50e3)
    rho_100 = atmosphere_density(100e3)

    print(f"  Vacuum velocity at surface:  {v_surface_vac:10.1f} m/s   (expect ~2324)")
    print(f"  Atmosphere rho at sea level:  {rho_0:12.4f} kg/m^3 (expect ~1.225)")
    print(f"  Atmosphere rho at 50 km:      {rho_50:12.2e} kg/m^3 (expect ~1.0e-3)")
    print(f"  Atmosphere rho at 100 km:     {rho_100:12.2e} kg/m^3 (expect ~5.6e-7)")
    print(f"  Release radius r_0:          {R_RELEASE:12.0f} m")
    print(f"  Tangential velocity v_0:     {V_TANGENTIAL_0:12.1f} m/s")
    print()


def print_trajectory_table(results):
    """Print key trajectory points at selected altitudes."""
    beta = results['beta']
    alt = results['alt']
    t = results['t']
    v_rel = results['v_rel']
    mach = results['mach']
    qdot = results['q_dot']
    gload = results['g_load']

    print(f"  Trajectory summary  (beta = {beta:.0f} kg/m^2)")
    print("-" * 90)
    print(f"  {'Alt (km)':>10s}  {'Time (s)':>10s}  {'v_rel (m/s)':>12s}  "
          f"{'Mach':>8s}  {'q_dot (kW/m^2)':>16s}  {'g-load':>8s}")
    print("-" * 90)

    for target_alt in TABLE_ALTITUDES:
        # Find the index closest to this altitude (descending)
        idx = np.argmin(np.abs(alt - target_alt))
        # Only print if we actually passed through this altitude
        if abs(alt[idx] - target_alt) < 5000:  # within 5 km tolerance
            print(f"  {alt[idx]/1e3:10.1f}  {t[idx]:10.1f}  {v_rel[idx]:12.1f}  "
                  f"{mach[idx]:8.2f}  {qdot[idx]/1e3:16.2f}  {gload[idx]:8.2f}")

    print("-" * 90)


def print_peak_values(results):
    """Print peak values from a trajectory."""
    beta = results['beta']
    alt = results['alt']
    t = results['t']
    v_rel = results['v_rel']
    qdot = results['q_dot']
    gload = results['g_load']

    idx_v = np.argmax(v_rel)
    idx_q = np.argmax(qdot)
    idx_g = np.argmax(gload)

    print(f"\n  Peak values  (beta = {beta:.0f} kg/m^2)")
    print("-" * 60)
    print(f"  Peak velocity:    {v_rel[idx_v]:10.1f} m/s   at {alt[idx_v]/1e3:7.1f} km")
    print(f"  Peak heat flux:   {qdot[idx_q]/1e3:10.2f} kW/m^2 at {alt[idx_q]/1e3:7.1f} km")
    print(f"  Peak g-load:      {gload[idx_g]:10.2f} g     at {alt[idx_g]/1e3:7.1f} km")
    print(f"  Total descent:    {t[-1]:10.1f} s   ({t[-1]/60:.1f} min)")
    print(f"  Total heat load:  {results['q_total'][-1]/1e6:10.2f} MJ/m^2")
    print("-" * 60)


# === CSV EXPORT ===

def save_csv(results):
    """Export full trajectory data to CSV."""
    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
    beta = results['beta']
    filename = f"freefall_beta{int(beta)}.csv"
    filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            't (s)', 'alt (m)', 'r (m)', 'theta (rad)',
            'v_r (m/s)', 'v_theta (m/s)', 'v_total (m/s)', 'v_rel (m/s)',
            'Mach', 'rho (kg/m^3)', 'q_dot (W/m^2)', 'q_total (J/m^2)',
            'g_load (g)'
        ])

        n = len(results['t'])
        for i in range(n):
            writer.writerow([
                f"{results['t'][i]:.3f}",
                f"{results['alt'][i]:.1f}",
                f"{results['r'][i]:.1f}",
                f"{results['theta'][i]:.8f}",
                f"{results['v_r'][i]:.3f}",
                f"{results['v_theta'][i]:.3f}",
                f"{results['v_total'][i]:.3f}",
                f"{results['v_rel'][i]:.3f}",
                f"{results['mach'][i]:.4f}",
                f"{results['rho'][i]:.6e}",
                f"{results['q_dot'][i]:.3f}",
                f"{results['q_total'][i]:.3f}",
                f"{results['g_load'][i]:.4f}",
            ])

    print(f"  Saved CSV: {filepath}")


# === GRAPH HELPERS ===

def _save_or_show(fig, filename, show_graphs):
    """Save figure to file or display interactively."""
    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        print(f"  Saved graph: {filepath}")
    else:
        plt.show()


# === GRAPHS ===

def plot_velocity_profile(results_dict, show_graphs):
    """Plot velocity vs altitude for multiple ballistic coefficients."""
    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    # Vacuum reference
    alt_vac, v_vac = vacuum_reference_curve()
    ax.plot(v_vac, alt_vac / 1e3, 'k--', linewidth=1.5, label='Vacuum', alpha=0.7)

    # Each beta
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results_dict)))
    for (beta, res), color in zip(sorted(results_dict.items()), colors):
        ax.plot(res['v_rel'], res['alt'] / 1e3, color=color, linewidth=1.5,
                label=f'$\\beta$ = {beta:.0f} kg/m$^2$')

    ax.set_xlabel('Velocity relative to atmosphere (m/s)', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Free-fall Velocity Profile from 275 km Release', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 280)

    _save_or_show(fig, f"01-velocity_profile.{cfg.GRAPH_FORMAT}", show_graphs)


def plot_mach_profile(results_dict, show_graphs):
    """Plot Mach number vs altitude for default beta."""
    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results_dict)))
    for (beta, res), color in zip(sorted(results_dict.items()), colors):
        ax.plot(res['mach'], res['alt'] / 1e3, color=color, linewidth=1.5,
                label=f'$\\beta$ = {beta:.0f} kg/m$^2$')

    ax.set_xlabel('Mach Number', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Mach Number vs Altitude during Free-fall Descent', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 280)

    _save_or_show(fig, f"02-mach_profile.{cfg.GRAPH_FORMAT}", show_graphs)


def plot_heat_flux(results_dict, show_graphs):
    """Plot stagnation heat flux vs altitude for multiple beta."""
    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results_dict)))
    for (beta, res), color in zip(sorted(results_dict.items()), colors):
        ax.plot(res['q_dot'] / 1e3, res['alt'] / 1e3, color=color, linewidth=1.5,
                label=f'$\\beta$ = {beta:.0f} kg/m$^2$')

    ax.set_xlabel('Stagnation Heat Flux (kW/m$^2$)', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Stagnation-Point Heat Flux vs Altitude', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 150)

    _save_or_show(fig, f"03-heat_flux_profile.{cfg.GRAPH_FORMAT}", show_graphs)


def plot_altitude_time(results_dict, show_graphs):
    """Plot altitude vs time for multiple beta."""
    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(results_dict)))
    for (beta, res), color in zip(sorted(results_dict.items()), colors):
        ax.plot(res['t'], res['alt'] / 1e3, color=color, linewidth=1.5,
                label=f'$\\beta$ = {beta:.0f} kg/m$^2$')

    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Altitude vs Time during Free-fall Descent', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 280)

    _save_or_show(fig, f"04-altitude_vs_time.{cfg.GRAPH_FORMAT}", show_graphs)


def plot_parametric(show_graphs):
    """Parametric sweep: peak velocity, peak heat flux, total heat load vs beta."""
    betas = np.array([50, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000])
    peak_v = np.zeros(len(betas))
    peak_q = np.zeros(len(betas))
    total_q = np.zeros(len(betas))
    peak_g = np.zeros(len(betas))

    for i, b in enumerate(betas):
        res = simulate_descent(beta=b)
        peak_v[i] = np.max(res['v_rel'])
        peak_q[i] = np.max(res['q_dot'])
        total_q[i] = res['q_total'][-1]
        peak_g[i] = np.max(res['g_load'])

    fig, axes = plt.subplots(1, 3, figsize=(cfg.GRAPH_WIDTH_INCHES * 1.4,
                                             cfg.GRAPH_HEIGHT_INCHES))

    # Peak velocity
    axes[0].semilogx(betas, peak_v, 'b-o', markersize=4)
    axes[0].set_xlabel('Ballistic Coefficient (kg/m$^2$)', fontsize=10)
    axes[0].set_ylabel('Peak Velocity (m/s)', fontsize=10)
    axes[0].set_title('Peak Velocity vs $\\beta$', fontsize=11)
    axes[0].grid(True, alpha=0.3)

    # Peak heat flux
    axes[1].semilogx(betas, peak_q / 1e3, 'r-o', markersize=4)
    axes[1].set_xlabel('Ballistic Coefficient (kg/m$^2$)', fontsize=10)
    axes[1].set_ylabel('Peak Heat Flux (kW/m$^2$)', fontsize=10)
    axes[1].set_title('Peak Heat Flux vs $\\beta$', fontsize=11)
    axes[1].grid(True, alpha=0.3)

    # Total heat load
    axes[2].semilogx(betas, total_q / 1e6, 'g-o', markersize=4)
    axes[2].set_xlabel('Ballistic Coefficient (kg/m$^2$)', fontsize=10)
    axes[2].set_ylabel('Total Heat Load (MJ/m$^2$)', fontsize=10)
    axes[2].set_title('Total Heat Load vs $\\beta$', fontsize=11)
    axes[2].grid(True, alpha=0.3)

    fig.suptitle('Parametric Sweep: Free-fall from 275 km', fontsize=13, y=1.02)
    fig.tight_layout()

    _save_or_show(fig, f"05-parametric_sweep.{cfg.GRAPH_FORMAT}", show_graphs)


# === COMMAND-LINE INTERFACE ===

def main():
    """Main entry point."""
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')

    # Defaults
    beta_override = None
    mass_override = None
    show_graphs = []

    # Parse CLI arguments
    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg.startswith("--beta="):
                beta_override = float(arg.split("=")[1])
            elif arg.startswith("--mass="):
                mass_override = float(arg.split("=")[1])
            elif arg == "--save":
                cfg.SAVE_GRAPHS = True
            elif arg == "--show":
                cfg.SAVE_GRAPHS = False
            elif arg.startswith("--outdir="):
                cfg.GRAPH_OUTPUT_DIR = arg.split("=")[1]
            elif arg.startswith("--dpi="):
                cfg.GRAPH_DPI = int(arg.split("=")[1])
            else:
                show_graphs.append(arg)

    # === VALIDATION ===

    print_validation()

    # === DEFAULT TRAJECTORY ===

    beta_default = beta_override if beta_override else BETA_DEFAULT
    mass_default = mass_override if mass_override else M_POD_DEFAULT

    print("=" * 72)
    print(f"FREE-FALL TRAJECTORY FROM {ALT_RELEASE/1e3:.0f} km")
    print("=" * 72)
    print(f"  Release altitude:           {ALT_RELEASE/1e3:.0f} km")
    print(f"  Release radius:             {R_RELEASE:.0f} m")
    print(f"  Initial tangential v:       {V_TANGENTIAL_0:.1f} m/s")
    print(f"  Orbital velocity at r_0:    {math.sqrt(GM/R_RELEASE):.1f} m/s")
    print(f"  v_0 / v_orbital:            {V_TANGENTIAL_0/math.sqrt(GM/R_RELEASE):.4f}")
    print(f"  Default ballistic coeff:    {beta_default:.0f} kg/m^2")
    print(f"  Default pod mass:           {mass_default:.0f} kg")
    print(f"  Time step:                  {cfg.DT_TRAJECTORY} s")
    print()

    # Run default trajectory
    res_default = simulate_descent(beta=beta_default, m_pod=mass_default)
    print_trajectory_table(res_default)
    print_peak_values(res_default)
    print()

    # === MULTI-BETA SWEEP ===

    print("=" * 72)
    print("BALLISTIC COEFFICIENT SWEEP")
    print("=" * 72)

    # Always include default sweep betas
    sweep_betas = BETA_SWEEP
    if beta_override and beta_override not in sweep_betas:
        sweep_betas = sorted(set(sweep_betas + [beta_override]))

    results_dict = {}
    for b in sweep_betas:
        res = simulate_descent(beta=b, m_pod=mass_default)
        results_dict[b] = res

        idx_v = np.argmax(res['v_rel'])
        idx_q = np.argmax(res['q_dot'])
        idx_g = np.argmax(res['g_load'])

        print(f"  beta={b:5.0f}  |  v_max={res['v_rel'][idx_v]:7.1f} m/s "
              f"@ {res['alt'][idx_v]/1e3:5.1f} km  |  "
              f"q_max={res['q_dot'][idx_q]/1e3:8.2f} kW/m^2 "
              f"@ {res['alt'][idx_q]/1e3:5.1f} km  |  "
              f"g_max={res['g_load'][idx_g]:6.2f} g  |  "
              f"t_total={res['t'][-1]:7.1f} s  |  "
              f"Q_total={res['q_total'][-1]/1e6:7.2f} MJ/m^2")

    print()

    # === CSV EXPORT ===

    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
    save_csv(res_default)

    # === GRAPHS ===

    if show_graphs:
        print()
        print("Generating graphs...")
        if "velocity" in show_graphs or "all" in show_graphs:
            plot_velocity_profile(results_dict, show_graphs)
        if "mach" in show_graphs or "all" in show_graphs:
            plot_mach_profile(results_dict, show_graphs)
        if "heat_flux" in show_graphs or "all" in show_graphs:
            plot_heat_flux(results_dict, show_graphs)
        if "altitude" in show_graphs or "all" in show_graphs:
            plot_altitude_time(results_dict, show_graphs)
        if "parametric" in show_graphs or "all" in show_graphs:
            plot_parametric(show_graphs)
        print("Done.")


if __name__ == "__main__":
    main()
