#!/usr/bin/env python3
"""
Lunar Mass Driver Plot Helpers — Optional plotting utilities

Standalone plotting functions for generating publication-quality figures
from lunar mass driver data. These can be used independently or imported
by lunar_simulation.py.

Usage:
    python lunar_plots.py  (generates example plots from a quick simulation)

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
import os


def plot_centrifugal_profile(v_max=15000, R_moon=1_737_400, g_moon=1.623,
                              m_per_m=900, save_path=None):
    """Plot centrifugal acceleration and rail force vs velocity.

    Shows the transition at v_orbital where centrifugal exceeds gravity.

    Args:
        v_max: Maximum velocity for plot (m/s)
        R_moon: Lunar radius (m)
        g_moon: Lunar surface gravity (m/s²)
        m_per_m: Total mass per metre of sled (kg/m)
        save_path: If set, save figure to this path
    """
    import matplotlib.pyplot as plt

    v_orbital = math.sqrt(g_moon * R_moon)

    velocities = [v_max * i / 500 for i in range(501)]
    a_cent = [v**2 / R_moon for v in velocities]
    g_eff = [g_moon - ac for ac in a_cent]
    F_rail = [m_per_m * max(0, ac - g_moon) for ac in a_cent]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Acceleration panel
    ax1.plot([v/1000 for v in velocities], a_cent, 'b-', linewidth=1.5,
             label='Centrifugal (v²/R)')
    ax1.axhline(g_moon, color='gray', ls='--', lw=1, label=f'g_moon = {g_moon} m/s²')
    ax1.axvline(v_orbital/1000, color='red', ls=':', lw=1, alpha=0.7,
                label=f'v_orbital = {v_orbital/1000:.2f} km/s')
    ax1.fill_between([v/1000 for v in velocities],
                     [g_moon]*len(velocities), a_cent,
                     where=[ac > g_moon for ac in a_cent],
                     alpha=0.15, color='red', label='Outward force zone')
    ax1.set_ylabel('Acceleration (m/s²)', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Centrifugal vs Gravity on Lunar Surface', fontsize=13,
                  fontweight='bold')

    # Rail force panel
    ax2.plot([v/1000 for v in velocities], [F/1000 for F in F_rail],
             '#D55E00', linewidth=1.5)
    ax2.axvline(v_orbital/1000, color='red', ls=':', lw=1, alpha=0.7)
    ax2.set_xlabel('Velocity (km/s)', fontsize=12)
    ax2.set_ylabel('Rail Force (kN/m)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Outward Rail Force per Metre', fontsize=13, fontweight='bold')

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close(fig)
    else:
        plt.show()


def plot_moon_phase_sweep(save_path=None):
    """Plot launch velocity vs Moon phase for Mars and Jupiter missions.

    Args:
        save_path: If set, save figure to this path
    """
    import matplotlib.pyplot as plt

    MU_SUN = 1.32712440018e20
    MU_EARTH = 3.986004418e14
    AU = 1.495978707e11
    A_EARTH = AU
    A_MARS = 1.524 * AU
    A_JUPITER = 5.204 * AU
    A_MOON_ORBIT = 384_400_000
    V_MOON = math.sqrt(MU_EARTH / A_MOON_ORBIT)
    V_ESC_E = math.sqrt(2.0 * MU_EARTH / A_MOON_ORBIT)
    V_ESC_M = math.sqrt(2.0 * 1.623 * 1_737_400)

    def v_inf(a_dep, a_arr):
        a_t = (a_dep + a_arr) / 2.0
        return math.sqrt(MU_SUN * (2.0/a_dep - 1.0/a_t)) - math.sqrt(MU_SUN/a_dep)

    def launch_v(vi_earth, phase):
        v_need = math.sqrt(vi_earth**2 + V_ESC_E**2)
        phi = math.radians(phase)
        b = 2.0 * V_MOON * math.cos(phi)
        c = V_MOON**2 - v_need**2
        vi_moon = max(0, (-b + math.sqrt(b**2 - 4*c)) / 2.0)
        return math.sqrt(vi_moon**2 + V_ESC_M**2)

    phases = list(range(0, 181, 5))
    vi_mars = v_inf(A_EARTH, A_MARS)
    vi_jup = v_inf(A_EARTH, A_JUPITER)

    v_mars = [launch_v(vi_mars, p)/1000 for p in phases]
    v_jup = [launch_v(vi_jup, p)/1000 for p in phases]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(phases, v_mars, 'b-', linewidth=2, label='Mars Hohmann')
    ax.plot(phases, v_jup, 'r-', linewidth=2, label='Jupiter Hohmann')
    ax.axhline(5.0, color='gray', ls='--', lw=1, alpha=0.7,
               label='v_cross = 5 km/s')
    ax.set_xlabel('Moon Orbital Phase (degrees)', fontsize=12)
    ax.set_ylabel('Required Launch Velocity (km/s)', fontsize=12)
    ax.set_title('Launch Velocity vs Moon Phase\n'
                 '(0° = optimal/prograde, 180° = worst/retrograde)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close(fig)
    else:
        plt.show()


def plot_three_regime_diagram(save_path=None):
    """Plot the three-regime thrust/power model.

    Shows F and P vs velocity with constant-thrust, voltage-limited,
    and grid-limited regimes annotated.

    Args:
        save_path: If set, save figure to this path
    """
    import matplotlib.pyplot as plt

    F_per_m = 7500.0  # N/m at N=50
    L_sled = 10000.0  # m
    F_total = F_per_m * L_sled  # 75 MN
    v_cross = 5000.0  # m/s
    P_total = F_total * v_cross  # 375 GW
    P_hvdc = 200e9  # 200 GW
    v_power = P_hvdc / F_total if P_hvdc < P_total else None

    v_max = 15000
    velocities = [v_max * i / 1000 for i in range(1, 1001)]
    thrust = []
    power = []

    for v in velocities:
        if v <= v_cross:
            F = F_total
            P = F * v
        else:
            P = P_total
            F = P_total / v
        if P_hvdc is not None and P > P_hvdc:
            P = P_hvdc
            F = P_hvdc / v
        thrust.append(F)
        power.append(P)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Thrust
    ax1.plot([v/1000 for v in velocities], [F/1e6 for F in thrust],
             'purple', linewidth=2)
    ax1.axvline(v_cross/1000, color='orange', ls='--', lw=1, alpha=0.7,
                label=f'v_cross = {v_cross/1000:.0f} km/s')
    if v_power:
        ax1.axvline(v_power/1000, color='red', ls='--', lw=1, alpha=0.7,
                    label=f'v_power = {v_power/1000:.1f} km/s')
    ax1.set_ylabel('Thrust (MN)', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Three-Regime Motor Model (N=50, 200 GW HVDC)',
                  fontsize=13, fontweight='bold')
    ax1.annotate('Constant\nThrust', xy=(2.5, F_total/1e6*0.8),
                fontsize=10, ha='center', color='purple')

    # Power
    ax2.plot([v/1000 for v in velocities], [P/1e9 for P in power],
             'green', linewidth=2)
    ax2.axhline(P_hvdc/1e9, color='red', ls='--', lw=1, alpha=0.7,
                label=f'P_HVDC = {P_hvdc/1e9:.0f} GW')
    ax2.axhline(P_total/1e9, color='orange', ls=':', lw=1, alpha=0.5,
                label=f'P_motor = {P_total/1e9:.0f} GW')
    ax2.axvline(v_cross/1000, color='orange', ls='--', lw=1, alpha=0.7)
    ax2.set_xlabel('Velocity (km/s)', fontsize=12)
    ax2.set_ylabel('Power (GW)', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    print("Generating lunar mass driver reference plots...")
    out_dir = "./graphs_lunar"

    plot_centrifugal_profile(save_path=os.path.join(out_dir, "centrifugal_profile.png"))
    plot_moon_phase_sweep(save_path=os.path.join(out_dir, "moon_phase_sweep.png"))
    plot_three_regime_diagram(save_path=os.path.join(out_dir, "three_regime_model.png"))

    print("Done.")
