#!/usr/bin/env python3
"""
Lunar Mass Driver Simulation — Main simulation loop and plotting

Simulates an LSM mass driver launching a 10 km sled with a 5,000-tonne
spacecraft from rest to target velocity on the Moon's equatorial surface.

The motor uses a three-regime model:
  1. Constant thrust (v < v_cross): F = F_per_m × L_sled
  2. Constant power (v > v_cross): P = P_per_m × L_sled  (voltage-limited)
  3. Grid-limited (P > P_HVDC_MAX): P = P_HVDC_MAX

Key lunar-specific features:
  - Centrifugal force tracked (exceeds gravity above v_orbital = 1,679 m/s)
  - Rail force per metre tracked for structural assessment
  - Lap counting around 10,917 km lunar circumference
  - No g-load throttling by default (unmanned cargo)
  - Optional rail force structural limit

Usage:
    python lunar_simulation.py [options] [graph keywords]

Options:
  --v_launch=N        Target velocity in m/s (default: 5000)
  --mass=N            Spacecraft mass in tonnes (default: 5000)
  --sled-mass=N       Sled mass per metre in kg/m (default: 400)
  --sled-length=N     Sled length in m (default: 10000)
  --p_hvdc_max=N      HVDC power limit in watts (default: 200e9)
  --no_power_limit    Disable HVDC power limit
  --g-limit=N         Set g-load limit (default: no limit)
  --rail-limit=N      Rail force limit in N/m (default: no limit)
  --quick             Use 1-second timesteps
  --no-graphs         Suppress graph output
  --save-csv          Export time-series data

Graph keywords:
  combined            12-panel combined plot
  all                 All individual plots
  velocity, accel, thrust, power, frequency
  gload, rail_force, centrifugal, laps
  hysteresis, radiator_width, E_hyst

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import math
import csv
import os

# Ensure UTF-8 output on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import lunar_config as cfg
import lunar_physics as phys

# =============================================================================
# DATA COLLECTION
# =============================================================================

data_t = []
data_v = []
data_x = []
data_a = []
data_F = []
data_P_mech = []
data_f_supply = []
data_a_centrifugal = []
data_a_net_outward = []
data_a_occupant_g = []
data_F_rail_per_m = []
data_KE = []
data_laps = []
data_P_hysteresis = []
data_E_hyst_cumulative = []
data_radiator_width = []


def clear_data():
    """Clear all data arrays for a fresh run."""
    global data_t, data_v, data_x, data_a, data_F, data_P_mech
    global data_f_supply, data_a_centrifugal, data_a_net_outward
    global data_a_occupant_g, data_F_rail_per_m, data_KE, data_laps
    global data_P_hysteresis, data_E_hyst_cumulative, data_radiator_width

    data_t = []
    data_v = []
    data_x = []
    data_a = []
    data_F = []
    data_P_mech = []
    data_f_supply = []
    data_a_centrifugal = []
    data_a_net_outward = []
    data_a_occupant_g = []
    data_F_rail_per_m = []
    data_KE = []
    data_laps = []
    data_P_hysteresis = []
    data_E_hyst_cumulative = []
    data_radiator_width = []


# =============================================================================
# MAIN SIMULATION LOOP
# =============================================================================

def run_simulation(quick_mode=False):
    """Run the lunar mass driver simulation.

    Uses the three-regime motor model:
      1. v < v_cross: F = F_total (constant thrust)
      2. v > v_cross: P = P_total (constant power, voltage-limited)
      3. P > P_HVDC_MAX: P = P_HVDC_MAX (grid-limited)

    Args:
        quick_mode: If True, use larger timesteps

    Returns:
        True if target velocity reached
    """
    clear_data()

    dt = cfg.DT_QUICK if quick_mode else cfg.DT

    # Motor parameters
    F_total = cfg.F_TOTAL
    P_total = cfg.P_TOTAL
    v_cross = cfg.V_CROSSOVER
    m_total = cfg.M_TOTAL
    m_per_m_total = cfg.M_PER_M_TOTAL

    # Physics parameters for hysteresis
    params = cfg.get_physics_params()
    n_units = cfg.N_UNITS
    active_length = cfg.N_UNITS * cfg.N_POLES_PER_UNIT * cfg.TAU_P

    print("\n" + "=" * 80)
    print("LUNAR MASS DRIVER SIMULATION")
    print("=" * 80)
    print_parameters()

    # Initialize state
    t = 0.0
    v = 0.0
    x = 0.0
    KE = 0.0
    E_hyst_cum = 0.0
    current_lap = 0

    # Max tracking
    max_values = {
        'thrust': 0.0,
        'power': 0.0,
        'g_load': 0.0,
        'rail_force': 0.0,
        'frequency': 0.0,
        'hysteresis': 0.0,
        'radiator_width': 0.0,
    }

    # Flags
    voltage_limited = False
    power_limited = False
    power_limit_engage_v = None
    rail_limit_hit = False
    g_load_limited = False
    g_load_engage_v = None

    print("\nStarting simulation...")
    print(f"Target velocity:  {cfg.V_LAUNCH:,.0f} m/s ({cfg.V_LAUNCH/1000:.1f} km/s)")
    print(f"v_crossover:      {v_cross:,.0f} m/s ({v_cross/1000:.1f} km/s)")
    print(f"F_total:          {F_total/1e6:.2f} MN")
    print(f"P_total:          {P_total/1e9:.2f} GW")
    print(f"Total mass:       {m_total:,.0f} kg ({m_total/1000:,.0f} t)")
    print(f"Peak accel:       {F_total/m_total/cfg.G_0:.3f} g")
    if cfg.P_HVDC_MAX is not None:
        print(f"HVDC power limit: {cfg.P_HVDC_MAX/1e9:.1f} GW")
    else:
        print("HVDC power limit: None (disabled)")
    if cfg.G_LOAD_MAX is not None:
        print(f"G-load limit:     {cfg.G_LOAD_MAX:.1f} g")
    else:
        print("G-load limit:     None (unmanned cargo)")
    if cfg.F_RAIL_MAX_PER_M is not None:
        print(f"Rail force limit: {cfg.F_RAIL_MAX_PER_M/1000:.1f} kN/m")
    else:
        print("Rail force limit: None (reporting only)")
    print(f"Timestep:         {dt} s")
    print(f"Track length:     {cfg.C_MOON/1000:,.0f} km (lunar circumference)")
    print(f"v_orbital:        {cfg.V_ORBITAL:,.1f} m/s (centrifugal = gravity)")
    print(f"v_escape:         {cfg.V_ESCAPE:,.1f} m/s")
    print()

    step_count = 0
    next_report_v = 500.0  # Report every 500 m/s

    # Main simulation loop
    while v < cfg.V_LAUNCH:
        # Determine thrust from three-regime model
        if v <= v_cross or v <= 0:
            F = F_total
            P_mech = F * max(v, 0)
        else:
            voltage_limited = True
            P_mech = P_total
            F = P_total / v

        # HVDC power limit
        if cfg.P_HVDC_MAX is not None and P_mech > cfg.P_HVDC_MAX and v > 0:
            power_limited = True
            if power_limit_engage_v is None:
                power_limit_engage_v = v
            P_mech = cfg.P_HVDC_MAX
            F = P_mech / v

        # Tangential acceleration
        a = F / m_total

        # G-load limiting (optional, for crewed missions)
        a_net_out = phys.calc_net_outward_acceleration(v, cfg.R_MOON, cfg.G_MOON)
        a_total_g = phys.calc_g_load(a, a_net_out, cfg.G_0)

        if cfg.G_LOAD_MAX is not None and a_total_g > cfg.G_LOAD_MAX:
            g_load_limited = True
            if g_load_engage_v is None:
                g_load_engage_v = v
            a_radial = abs(a_net_out)
            a_max_tang_sq = (cfg.G_LOAD_MAX * cfg.G_0)**2 - a_radial**2
            if a_max_tang_sq > 0:
                a = min(a, math.sqrt(a_max_tang_sq))
                F = a * m_total
                P_mech = F * v
            else:
                print(f"\nRadial g-load alone exceeds {cfg.G_LOAD_MAX:.1f}g at "
                      f"v={v/1000:.2f} km/s")
                break
            a_total_g = phys.calc_g_load(a, a_net_out, cfg.G_0)

        # Rail force check
        F_rail = phys.calc_rail_force_per_m(v, cfg.R_MOON, cfg.G_MOON, m_per_m_total)
        if cfg.F_RAIL_MAX_PER_M is not None and F_rail > cfg.F_RAIL_MAX_PER_M:
            if not rail_limit_hit:
                rail_limit_hit = True
                print(f"\n  Rail force limit ({cfg.F_RAIL_MAX_PER_M/1000:.1f} kN/m) "
                      f"reached at v = {v/1000:.2f} km/s")
                print(f"  Maximum velocity for this rail rating reached.")
            break

        # Centrifugal and supply frequency
        a_cent = phys.calc_centrifugal_acceleration(v, cfg.R_MOON)
        f_supply = phys.calc_supply_frequency(v, cfg.TAU_P)

        # HTS hysteresis
        P_hyst = phys.calc_hysteresis_power(f_supply, cfg.I_TARGET, params, n_units)
        E_hyst_cum += P_hyst * dt

        # Radiator width for HTS cooling
        epsilon = 0.9
        sigma = 5.67e-8
        T_rad = 100.0
        q_rad = epsilon * sigma * T_rad**4
        radiator_w = P_hyst / (2.0 * q_rad * active_length) if active_length > 0 else 0.0

        # Update velocity and position
        v_new = v + a * dt
        x_new = x + v * dt
        KE = phys.calc_kinetic_energy(m_total, v_new)

        # Lap counting
        new_lap = int(x_new / cfg.C_MOON)
        if new_lap > current_lap:
            current_lap = new_lap

        # Track maximums
        max_values['thrust'] = max(max_values['thrust'], F)
        max_values['power'] = max(max_values['power'], P_mech)
        max_values['g_load'] = max(max_values['g_load'], a_total_g)
        max_values['rail_force'] = max(max_values['rail_force'], F_rail)
        max_values['frequency'] = max(max_values['frequency'], f_supply)
        max_values['hysteresis'] = max(max_values['hysteresis'], P_hyst)
        max_values['radiator_width'] = max(max_values['radiator_width'], radiator_w)

        # Data collection
        data_t.append(t)
        data_v.append(v)
        data_x.append(x)
        data_a.append(a)
        data_F.append(F)
        data_P_mech.append(P_mech)
        data_f_supply.append(f_supply)
        data_a_centrifugal.append(a_cent)
        data_a_net_outward.append(a_net_out)
        data_a_occupant_g.append(a_total_g)
        data_F_rail_per_m.append(F_rail)
        data_KE.append(KE)
        data_laps.append(x / cfg.C_MOON)
        data_P_hysteresis.append(P_hyst)
        data_E_hyst_cumulative.append(E_hyst_cum)
        data_radiator_width.append(radiator_w)

        # Update state
        v = v_new
        x = x_new
        t += dt
        step_count += 1

        # Periodic reporting (by velocity milestones)
        if v >= next_report_v:
            laps = x / cfg.C_MOON
            regime = "thrust" if v <= v_cross else "V-lim"
            if cfg.P_HVDC_MAX is not None and P_mech >= cfg.P_HVDC_MAX * 0.99:
                regime = "P-lim"
            rail_str = f" rail={F_rail/1000:.1f}kN/m" if F_rail > 0 else ""
            print(f"  v = {v/1000:7.2f} km/s | t = {t/60:6.1f} min | "
                  f"F = {F/1e6:6.2f} MN | P = {P_mech/1e9:6.1f} GW | "
                  f"g = {a_total_g:.2f} | laps = {laps:.2f} [{regime}]{rail_str}")
            next_report_v += 500.0

        # Safety check
        if t > 72000:
            print("\nSimulation time exceeded 20 hours — stopping")
            break

    # Final report
    total_time_min = t / 60.0
    total_time_hr = t / 3600.0
    total_distance_km = x / 1000.0
    total_laps = x / cfg.C_MOON

    print("\n" + "=" * 80)
    print("SIMULATION COMPLETE")
    print("=" * 80)
    print(f"Final velocity:       {v:,.1f} m/s ({v/1000:.2f} km/s)")
    print(f"Target velocity:      {cfg.V_LAUNCH:,.1f} m/s ({cfg.V_LAUNCH/1000:.1f} km/s)")
    print(f"Achievement:          {100 * v / cfg.V_LAUNCH:.1f}%")
    print(f"Total time:           {total_time_min:.1f} min ({total_time_hr:.3f} hr)")
    print(f"Distance traveled:    {total_distance_km:,.1f} km")
    print(f"Laps of Moon:         {total_laps:.2f}")
    print(f"Final KE:             {KE/1e12:.3f} TJ ({KE/3.6e9:.1f} GWh)")
    print()
    print(f"Peak thrust:          {max_values['thrust']/1e6:.2f} MN")
    print(f"Peak power:           {max_values['power']/1e9:.2f} GW")
    print(f"Peak acceleration:    {max_values['thrust']/m_total/cfg.G_0:.3f} g")
    print(f"Peak g-load:          {max_values['g_load']:.3f} g")
    print(f"Peak frequency:       {max_values['frequency']:.1f} Hz")
    print(f"Peak rail force:      {max_values['rail_force']/1000:.1f} kN/m")
    print()
    print(f"Total HTS hysteresis: {E_hyst_cum:.3e} J ({E_hyst_cum/3.6e6:.1f} kWh)")
    print(f"Peak HTS hysteresis:  {max_values['hysteresis']/1e3:.2f} kW")
    print(f"Peak radiator width:  {max_values['radiator_width']:.4f} m")

    # Per-coil hysteresis breakdown
    n_phases = cfg.N_PHASES
    n_sides = cfg.N_STATOR_SIDES
    n_coils_per_unit = n_phases * n_sides
    total_track_coils = int(cfg.C_MOON / (cfg.N_POLES_PER_UNIT * cfg.TAU_P + cfg.L_GAP)) * n_coils_per_unit * cfg.N_POLES_PER_UNIT
    if total_laps > 0 and total_track_coils > 0:
        # Each fixed stator coil only sees the sled for L_sled/v seconds per pass
        n_coils_under_sled = n_units * n_coils_per_unit * cfg.N_POLES_PER_UNIT
        n_passes = max(1.0, total_laps)
        coils_touched = total_track_coils * min(1.0, n_passes)
        E_per_coil = E_hyst_cum / coils_touched if coils_touched > 0 else 0.0
        exposure_at_final = cfg.L_SLED / v if v > 0 else 0.0
        print()
        print(f"  Per-coil hysteresis breakdown:")
        print(f"    Track stator coils: ~{total_track_coils:,.0f}")
        print(f"    Coils under sled:   ~{n_coils_under_sled:,.0f}")
        print(f"    Laps completed:     {total_laps:.2f}")
        print(f"    Avg energy/coil:    {E_per_coil:.1f} J (over {n_passes:.1f} passes)")
        print(f"    Exposure per pass:  {exposure_at_final:.2f} s (at final v={v/1000:.1f} km/s)")
    print()
    print(f"Voltage limited:      {voltage_limited}")
    if power_limited:
        print(f"HVDC power limit:     {cfg.P_HVDC_MAX/1e9:.1f} GW "
              f"(engaged at {power_limit_engage_v/1000:.2f} km/s)")
    elif cfg.P_HVDC_MAX is not None:
        print(f"HVDC power limit:     {cfg.P_HVDC_MAX/1e9:.1f} GW (not reached)")
    else:
        print("HVDC power limit:     None (disabled)")
    if g_load_limited:
        print(f"G-load limited:       True (engaged at {g_load_engage_v/1000:.2f} km/s)")
    elif cfg.G_LOAD_MAX is not None:
        print(f"G-load limited:       False (limit: {cfg.G_LOAD_MAX:.1f}g)")
    else:
        print("G-load limited:       N/A (no limit set)")
    if rail_limit_hit:
        print(f"Rail force limit:     Hit at v = {v/1000:.2f} km/s")
    elif cfg.F_RAIL_MAX_PER_M is not None:
        print(f"Rail force limit:     {cfg.F_RAIL_MAX_PER_M/1000:.1f} kN/m (not reached)")
    else:
        print("Rail force limit:     None (reporting only)")

    # Centrifugal summary
    print()
    if v > cfg.V_ORBITAL:
        print(f"v > v_orbital:        Yes — sled above orbital velocity")
        print(f"  v_orbital:          {cfg.V_ORBITAL:,.1f} m/s")
        print(f"  Peak outward rail:  {max_values['rail_force']/1000:.1f} kN/m")
    else:
        print(f"v < v_orbital:        Sled below orbital velocity (gravity holds sled)")

    if v > cfg.V_ESCAPE:
        print(f"v > v_escape:         Yes — sled exceeds lunar escape velocity")

    print("=" * 80)

    return v >= cfg.V_LAUNCH * 0.99


# =============================================================================
# PARAMETER DISPLAY
# =============================================================================

def print_parameters():
    """Display configuration parameters."""
    derived = cfg.calc_derived()

    print(f"\nConfiguration:")
    print(f"  Track:                  Lunar equator ({cfg.C_MOON/1000:,.0f} km)")
    print(f"  Lunar gravity:          {cfg.G_MOON:.3f} m/s²")
    print(f"  v_orbital:              {cfg.V_ORBITAL:,.1f} m/s")
    print(f"  v_escape:               {cfg.V_ESCAPE:,.1f} m/s")
    print()
    print(f"  Pole pitch:             {cfg.TAU_P:.1f} m")
    print(f"  Stator turns:           {cfg.N_STATOR}")
    print(f"  Coil height:            {cfg.W_COIL:.2f} m")
    print(f"  Voltage limit:          {cfg.V_COIL_LIMIT / 1000:.0f} kV")
    print(f"  B_sled (nominal):       {cfg.B_SLED_NOMINAL:.3f} T")
    print()
    print(f"  F/L (at I_TARGET):      {derived['F_per_m']:,.0f} N/m")
    print(f"  v_crossover:            {derived['V_crossover']:,.0f} m/s")
    print(f"  P/L (invariant):        {derived['P_per_m']/1e6:.1f} MW/m")
    print(f"  F_total (sled):         {derived['F_total']/1e6:.2f} MN")
    print(f"  P_total (sled):         {derived['P_total']/1e9:.2f} GW")
    print()
    print(f"  Sled length:            {cfg.L_SLED:,.0f} m")
    print(f"  Sled mass/m:            {cfg.M_SLED_PER_M} kg/m")
    print(f"  Sled hardware mass:     {derived['M_sled_hardware']:,.0f} kg "
          f"({derived['M_sled_hardware']/1000:,.0f} t)")
    print(f"  Spacecraft mass:        {cfg.M_SPACECRAFT:,.0f} kg "
          f"({cfg.M_SPACECRAFT/1000:,.0f} t)")
    print(f"  Total mass:             {derived['M_total']:,.0f} kg "
          f"({derived['M_total']/1000:,.0f} t)")
    print(f"  Total mass/m:           {derived['M_per_m_total']:.0f} kg/m")
    print()
    print(f"  Target velocity:        {cfg.V_LAUNCH:,.0f} m/s ({cfg.V_LAUNCH/1000:.1f} km/s)")
    print(f"  Target KE:              {derived['KE_target'] / 1e12:.3f} TJ")
    if cfg.P_HVDC_MAX is not None:
        print(f"  HVDC power limit:       {cfg.P_HVDC_MAX/1e9:.1f} GW")
    else:
        print("  HVDC power limit:       None (disabled)")
    print()


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def _param_subtitle():
    """Build parameter subtitle string."""
    return (f"(Lunar, N={cfg.N_STATOR}, "
            f"\u03c4\u209a={cfg.TAU_P:.0f} m, "
            f"m={cfg.M_SPACECRAFT/1000:,.0f} t, "
            f"v={cfg.V_LAUNCH/1000:.1f} km/s)")


def _time_axis(data_t):
    """Choose time unit based on duration."""
    if not data_t:
        return data_t, "Time (min)"
    t_total_min = data_t[-1] / 60.0
    if t_total_min > 120:
        return [t / 3600 for t in data_t], "Time (hr)"
    else:
        return [t / 60 for t in data_t], "Time (min)"


def plot_combined():
    """Create a 12-panel combined plot."""
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(16, 16))
    gs = gridspec.GridSpec(4, 3, hspace=0.3, wspace=0.3)
    tx, t_label = _time_axis(data_t)

    # 1. Velocity
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(tx, [v / 1000 for v in data_v], 'b-', linewidth=1)
    ax.set_ylabel('Velocity (km/s)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Velocity')

    # 2. Acceleration
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(tx, [a / cfg.G_0 for a in data_a], 'r-', linewidth=1)
    ax.set_ylabel('Acceleration (g)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Acceleration')

    # 3. Thrust
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(tx, [F / 1e6 for F in data_F], 'purple', linewidth=1)
    ax.set_ylabel('Thrust (MN)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Thrust')

    # 4. Power
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(tx, [P / 1e9 for P in data_P_mech], 'green', linewidth=1)
    ax.set_ylabel('Power (GW)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Mechanical Power')

    # 5. G-load
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(tx, data_a_occupant_g, 'darkred', linewidth=1)
    if cfg.G_LOAD_MAX is not None:
        ax.axhline(cfg.G_LOAD_MAX, color='red', ls='--', alpha=0.5, lw=1,
                    label=f'{cfg.G_LOAD_MAX:.0f}g limit')
        ax.legend()
    ax.set_ylabel('G-load (g)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('G-Load')

    # 6. Supply frequency
    ax = fig.add_subplot(gs[1, 2])
    ax.plot(tx, data_f_supply, 'brown', linewidth=1)
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Supply Frequency')

    # 7. Rail force per metre
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(tx, [F / 1000 for F in data_F_rail_per_m], '#D55E00', linewidth=1)
    ax.set_ylabel('Rail Force (kN/m)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Outward Rail Force')

    # 8. Centrifugal vs gravity
    ax = fig.add_subplot(gs[2, 1])
    ax.plot(tx, data_a_centrifugal, 'b-', linewidth=1, label='Centrifugal')
    ax.axhline(cfg.G_MOON, color='gray', ls='--', alpha=0.7, lw=1, label='g_moon')
    ax.plot(tx, [max(0, a) for a in data_a_net_outward], 'r-', linewidth=1,
            label='Net outward')
    ax.set_ylabel('Acceleration (m/s²)')
    ax.set_xlabel(t_label)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_title('Centrifugal vs Gravity')

    # 9. Laps
    ax = fig.add_subplot(gs[2, 2])
    ax.plot(tx, data_laps, '#009E73', linewidth=1)
    ax.set_ylabel('Laps')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Laps of Moon')

    # 10. HTS Hysteresis
    ax = fig.add_subplot(gs[3, 0])
    ax.plot(tx, [P / 1000 for P in data_P_hysteresis], '#CC79A7', linewidth=1)
    ax.set_ylabel('HTS Hyst (kW)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('HTS Hysteresis Power')

    # 11. Radiator width
    ax = fig.add_subplot(gs[3, 1])
    ax.plot(tx, data_radiator_width, '#0072B2', linewidth=1)
    ax.set_ylabel('Radiator Width (m)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Cryo Radiator Width')

    # 12. Cumulative hysteresis energy
    ax = fig.add_subplot(gs[3, 2])
    ax.plot(tx, [E / 3.6e6 for E in data_E_hyst_cumulative], '#009E73', linewidth=1)
    ax.set_ylabel('Energy (kWh)')
    ax.set_xlabel(t_label)
    ax.grid(True, alpha=0.3)
    ax.set_title('Cumulative Hysteresis Energy')

    fig.suptitle(f"Lunar Mass Driver Launch Profile\n{_param_subtitle()}",
                 fontsize=14, fontweight='bold')

    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR,
                                f"lunar_combined.{cfg.GRAPH_FORMAT}")
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight')
        print(f"Saved combined plot: {filepath}")
        plt.close(fig)
    else:
        plt.show()


def plot_individual(keyword):
    """Create individual plot based on keyword."""
    import matplotlib.pyplot as plt

    tx, t_label = _time_axis(data_t)

    plots = {
        'velocity': (data_v, 'Velocity (km/s)', 'Velocity', 'blue',
                      lambda y: y / 1000),
        'accel': (data_a, 'Acceleration (g)', 'Acceleration', 'red',
                   lambda y: y / cfg.G_0),
        'thrust': (data_F, 'Thrust (MN)', 'Thrust', 'purple',
                    lambda y: y / 1e6),
        'power': (data_P_mech, 'Power (GW)', 'Mechanical Power', 'green',
                   lambda y: y / 1e9),
        'frequency': (data_f_supply, 'Frequency (Hz)', 'Supply Frequency',
                       'brown', lambda y: y),
        'gload': (data_a_occupant_g, 'G-load (g)', 'G-Load', 'darkred',
                   lambda y: y),
        'rail_force': (data_F_rail_per_m, 'Rail Force (kN/m)',
                        'Outward Rail Force', '#D55E00', lambda y: y / 1000),
        'centrifugal': (data_a_centrifugal, 'Acceleration (m/s²)',
                         'Centrifugal Acceleration', 'blue', lambda y: y),
        'laps': (data_laps, 'Laps', 'Laps of Moon', '#009E73', lambda y: y),
        'hysteresis': (data_P_hysteresis, 'HTS Hyst (kW)',
                        'HTS Hysteresis Power', '#CC79A7', lambda y: y / 1000),
        'radiator_width': (data_radiator_width, 'Radiator Width (m)',
                            'Cryo Radiator Width', '#0072B2', lambda y: y),
        'E_hyst': (data_E_hyst_cumulative, 'Hyst Energy (kWh)',
                    'Cumulative Hysteresis Energy', '#009E73',
                    lambda y: y / 3.6e6),
    }

    if keyword not in plots:
        print(f"Unknown plot keyword: {keyword}")
        return

    data, ylabel, title, color, transform = plots[keyword]
    y_plot = [transform(y) for y in data]

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))
    ax.plot(tx, y_plot, color=color, linewidth=1.5)
    ax.set_xlabel(t_label, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"{title}\n{_param_subtitle()}", fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR,
                                f"lunar_{keyword}.{cfg.GRAPH_FORMAT}")
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
        plt.close(fig)
    else:
        plt.show()


def save_csv():
    """Export time-series data to CSV."""
    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
    filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, "lunar_data.csv")

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            't (s)', 'v (m/s)', 'x (m)', 'a (m/s^2)', 'F (N)', 'P_mech (W)',
            'f_supply (Hz)', 'a_centrifugal (m/s^2)', 'a_net_outward (m/s^2)',
            'g_load (g)', 'F_rail_per_m (N/m)', 'KE (J)', 'laps',
            'P_hysteresis (W)', 'E_hyst_cumulative (J)', 'radiator_width (m)'
        ])

        for i in range(len(data_t)):
            writer.writerow([
                data_t[i], data_v[i], data_x[i], data_a[i], data_F[i],
                data_P_mech[i], data_f_supply[i], data_a_centrifugal[i],
                data_a_net_outward[i], data_a_occupant_g[i], data_F_rail_per_m[i],
                data_KE[i], data_laps[i], data_P_hysteresis[i],
                data_E_hyst_cumulative[i], data_radiator_width[i]
            ])

    print(f"Saved CSV data: {filepath}")


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def main():
    """Main entry point with command-line argument parsing."""
    quick_mode = False
    save_csv_flag = False
    no_graphs = False
    plot_keywords = []

    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--quick":
                quick_mode = True
            elif arg == "--save-csv":
                save_csv_flag = True
            elif arg == "--no-graphs":
                no_graphs = True
            elif arg.startswith("--v_launch="):
                cfg.V_LAUNCH = float(arg.split("=")[1])
            elif arg.startswith("--mass="):
                cfg.M_SPACECRAFT = float(arg.split("=")[1]) * 1000
            elif arg.startswith("--sled-mass="):
                cfg.M_SLED_PER_M = float(arg.split("=")[1])
            elif arg.startswith("--sled-length="):
                cfg.L_SLED = float(arg.split("=")[1])
            elif arg.startswith("--p_hvdc_max="):
                cfg.P_HVDC_MAX = float(arg.split("=")[1])
            elif arg == "--no_power_limit":
                cfg.P_HVDC_MAX = None
            elif arg.startswith("--g-limit="):
                cfg.G_LOAD_MAX = float(arg.split("=")[1])
            elif arg.startswith("--rail-limit="):
                cfg.F_RAIL_MAX_PER_M = float(arg.split("=")[1])
            else:
                plot_keywords.append(arg)

    # Recalculate derived parameters
    cfg.calc_derived()

    # Run simulation
    success = run_simulation(quick_mode)

    # Save CSV if requested
    if save_csv_flag:
        save_csv()

    # Generate plots
    if not no_graphs and plot_keywords:
        if "combined" in plot_keywords:
            plot_combined()
        if "all" in plot_keywords:
            for kw in ['velocity', 'accel', 'thrust', 'power', 'frequency',
                        'gload', 'rail_force', 'centrifugal', 'laps',
                        'hysteresis', 'radiator_width', 'E_hyst']:
                plot_individual(kw)
        else:
            for kw in plot_keywords:
                if kw != "combined":
                    plot_individual(kw)


if __name__ == "__main__":
    main()
