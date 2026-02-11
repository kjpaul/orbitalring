#!/usr/bin/env python3
"""
LSM Mass Driver Simulation - Main simulation loop and plotting

This simulates an LSM (Linear Synchronous Motor) mass driver launching a
5,000 m sled with a 5,000-tonne spacecraft from rest to 30 km/s on an orbital
ring at 250 km altitude.

Key differences from LIM:
  - No slip (synchronous operation)
  - No eddy current losses (no sled heating)
  - Load angle δ is the control variable (not slip velocity)
  - Voltage limit is the primary constraint at high speed
  - Single pole pitch works across entire velocity range
  - No thermal limit on the sled

Usage:
    python lsm_simulation.py [options] [graph keywords]

Options:
  --v_launch=N      Target velocity in m/s (default: 15000)
  --accel=N         Maximum acceleration in g (default: 0.5)
  --mass=N          Spacecraft mass in tonnes (default: 500)
  --B_sled=N        Sled DC field in T (default: 0.10)
  --tau_p=N         Pole pitch in m (default: 20)
  --N_stator=N      Stator turns (default: 10)
  --adjustable-B    Allow controller to reduce B_sled
  --sled-length=N   Sled length in m (default: 5000)
  --quick           Use 1-second timesteps
  --save-csv        Export time-series data

Graph keywords:
  combined          12-panel combined plot
  all               All individual plots
  velocity, accel, thrust, power, delta, frequency, voltage, bfield, gload
  hysteresis, radiator_width, E_hyst

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import math
import csv
import os

# Ensure UTF-8 output on Windows (τ_p, etc.)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Import configuration and physics modules
import lsm_config as cfg
import lsm_physics as phys

# =============================================================================
# DATA COLLECTION
# =============================================================================

# Time-series data arrays
data_t = []
data_v = []
data_x = []
data_a = []
data_F = []
data_P_mech = []
data_I_stator = []
data_B_stator = []
data_B_sled = []
data_delta = []
data_f_supply = []
data_V_coil = []
data_a_centrifugal = []
data_a_radial_net = []
data_a_occupant_g = []
data_KE = []
data_x_fraction = []
data_P_hysteresis = []
data_E_hyst_cumulative = []
data_radiator_width = []


def clear_data():
    """Clear all data arrays for a fresh run."""
    global data_t, data_v, data_x, data_a, data_F, data_P_mech
    global data_I_stator, data_B_stator, data_B_sled, data_delta
    global data_f_supply, data_V_coil, data_a_centrifugal
    global data_a_radial_net, data_a_occupant_g, data_KE, data_x_fraction
    global data_P_hysteresis, data_E_hyst_cumulative, data_radiator_width

    data_t = []
    data_v = []
    data_x = []
    data_a = []
    data_F = []
    data_P_mech = []
    data_I_stator = []
    data_B_stator = []
    data_B_sled = []
    data_delta = []
    data_f_supply = []
    data_V_coil = []
    data_a_centrifugal = []
    data_a_radial_net = []
    data_a_occupant_g = []
    data_KE = []
    data_x_fraction = []
    data_P_hysteresis = []
    data_E_hyst_cumulative = []
    data_radiator_width = []


# =============================================================================
# MAIN SIMULATION LOOP
# =============================================================================

def run_lsm_simulation(quick_mode=False):
    """Run the LSM mass driver simulation.

    Args:
        quick_mode: If True, use larger timesteps (1 s instead of 0.1 s)

    Returns:
        success: True if simulation completed successfully
    """
    clear_data()

    # Get physics parameters
    params = cfg.get_physics_params()

    # Time step
    dt = cfg.DT_QUICK if quick_mode else cfg.DT

    print("\n" + "=" * 80)
    print("LSM MASS DRIVER SIMULATION")
    print("=" * 80)
    print_parameters()

    # Initialize state
    t = 0.0
    v = 0.0  # Start from rest
    x = 0.0  # Distance traveled
    KE = 0.0
    E_hyst_cum = 0.0

    # Controller state
    I_stator = cfg.I_MIN
    B_sled = cfg.B_SLED_NOMINAL
    delta = 0.0
    F = 0.0
    a = 0.0

    # Target acceleration (in m/s²)
    a_target = cfg.A_MAX_G * cfg.G_0

    # Max/min tracking
    max_values = {
        'thrust': 0.0,
        'power': 0.0,
        'current': 0.0,
        'voltage': 0.0,
        'frequency': 0.0,
        'delta': 0.0,
        'g_load': 0.0,
        'hysteresis': 0.0,
        'radiator_width': 0.0,
    }

    # Flags
    voltage_limited = False
    thrust_limited = False

    print(f"\nStarting simulation...")
    print(f"Target velocity: {cfg.V_LAUNCH:,.0f} m/s")
    print(f"Target acceleration: {cfg.A_MAX_G:.2f} g = {a_target:.3f} m/s²")
    print(f"Total mass: {cfg.M_TOTAL:,.0f} kg")
    print("Stator: Air-core")
    print(f"B_sled adjustable: {cfg.B_SLED_ADJUSTABLE}")
    print(f"Timestep: {dt} s")
    print()

    step_count = 0
    next_report_time = 600.0  # Report every 10 minutes

    # Main simulation loop
    while v < cfg.V_LAUNCH:

        # Controller: adjust current to achieve desired thrust
        # This iterates to find the right current and load angle
        for _ in range(5):  # Iterate to converge
            # Calculate supply frequency (synchronous)
            f_supply = phys.calc_frequency(v, cfg.TAU_P)

            # Calculate back-EMF
            V_coil = phys.calc_EMF(cfg.N_STATOR, v, B_sled, cfg.W_COIL)

            # Check voltage limit
            if V_coil > cfg.V_COIL_LIMIT:
                voltage_limited = True
                if cfg.B_SLED_ADJUSTABLE:
                    # Reduce B_sled to meet voltage limit
                    B_sled = phys.calc_max_B_sled_for_voltage(
                        cfg.V_COIL_LIMIT, cfg.N_STATOR, v, cfg.W_COIL
                    )
                    V_coil = cfg.V_COIL_LIMIT
                else:
                    # Can't reduce B_sled - voltage limit prevents further acceleration
                    print(f"\nVoltage limit reached at v = {v:.1f} m/s")
                    print(f"Cannot achieve target velocity without adjustable B_sled")
                    break

            # Calculate stator field
            B_stator = phys.calc_B_stator(
                I_stator, cfg.N_STATOR, cfg.W_COIL, cfg.G_GAP
            )

            # Desired thrust for target acceleration
            F_desired = cfg.M_TOTAL * a_target

            # Calculate required load angle for desired thrust
            delta, limited = phys.calc_delta_for_thrust(
                F_desired, B_stator, B_sled, cfg.A_ACTIVE_TOTAL, cfg.DELTA_MAX
            )

            if limited:
                thrust_limited = True
                # Need more current to achieve target thrust
                if I_stator < cfg.I_PEAK:
                    I_stator = min(cfg.I_PEAK, I_stator * 1.1)
                    continue  # Recalculate with new current
                else:
                    # Already at max current, accept limited thrust
                    F = phys.calc_thrust(B_stator, B_sled, cfg.DELTA_MAX, cfg.A_ACTIVE_TOTAL)
                    a = F / cfg.M_TOTAL
                    break
            else:
                # Can achieve desired thrust
                F = F_desired
                a = a_target

                # Check if we can reduce current (operating below target load angle)
                if delta < cfg.DELTA_TARGET and I_stator > cfg.I_MIN:
                    I_stator = max(cfg.I_MIN, I_stator * 0.99)
                break

        # Exit if voltage limit hit
        if V_coil > cfg.V_COIL_LIMIT and not cfg.B_SLED_ADJUSTABLE:
            break

        # Calculate power
        P_mech = phys.calc_power_mechanical(F, v)

        # HTS hysteresis losses and radiator width
        P_hyst = phys.calc_hysteresis_power(f_supply, I_stator, params,
                                            cfg.N_UNITS, cfg.NORRIS_HYSTERESIS)
        E_hyst_cum += P_hyst * dt
        active_length = cfg.N_UNITS * cfg.N_POLES_PER_UNIT * cfg.TAU_P
        radiator_w = phys.calc_radiator_width(
            P_hyst, cfg.T_STATOR, cfg.T_RADIATOR_HOT,
            cfg.CRYO_EFF, cfg.EM_HEATSINK, active_length, cfg.T_SPACE)

        # Calculate accelerations and g-loads
        a_centrifugal = phys.calc_centrifugal_acceleration(v, cfg.R_ORBIT)
        a_radial_net = phys.calc_net_radial_acceleration(v, cfg.R_ORBIT, cfg.G_250)
        a_occupant_g = phys.calc_occupant_g_load(a, a_radial_net, cfg.G_0)

        # Update velocity and position
        v_new = v + a * dt
        x_new = x + v * dt  # Use old velocity for distance (midpoint would be more accurate)

        # Update kinetic energy
        KE = phys.calc_kinetic_energy(cfg.M_TOTAL, v_new)

        # Track maximums
        max_values['thrust'] = max(max_values['thrust'], F)
        max_values['power'] = max(max_values['power'], P_mech)
        max_values['current'] = max(max_values['current'], I_stator)
        max_values['voltage'] = max(max_values['voltage'], V_coil)
        max_values['frequency'] = max(max_values['frequency'], f_supply)
        max_values['delta'] = max(max_values['delta'], delta)
        max_values['g_load'] = max(max_values['g_load'], a_occupant_g)
        max_values['hysteresis'] = max(max_values['hysteresis'], P_hyst)
        max_values['radiator_width'] = max(max_values['radiator_width'], radiator_w)

        # Data collection
        data_t.append(t)
        data_v.append(v)
        data_x.append(x)
        data_a.append(a)
        data_F.append(F)
        data_P_mech.append(P_mech)
        data_I_stator.append(I_stator)
        data_B_stator.append(B_stator)
        data_B_sled.append(B_sled)
        data_delta.append(delta)
        data_f_supply.append(f_supply)
        data_V_coil.append(V_coil)
        data_a_centrifugal.append(a_centrifugal)
        data_a_radial_net.append(a_radial_net)
        data_a_occupant_g.append(a_occupant_g)
        data_KE.append(KE)
        data_x_fraction.append(x / cfg.L_RING)
        data_P_hysteresis.append(P_hyst)
        data_E_hyst_cumulative.append(E_hyst_cum)
        data_radiator_width.append(radiator_w)

        # Update state
        v = v_new
        x = x_new
        t += dt
        step_count += 1

        # Periodic reporting
        if t >= next_report_time:
            minutes = t / 60.0
            km = x / 1000.0
            print(f"t = {minutes:7.1f} min | v = {v:8.1f} m/s | x = {km:10.1f} km | "
                  f"F = {F/1e6:6.2f} MN | P = {P_mech/1e9:6.2f} GW | "
                  f"I = {I_stator:6.1f} A | delta = {math.degrees(delta):5.1f} deg | g = {a_occupant_g:.3f}")
            next_report_time += 600.0  # Next report in 10 minutes

        # Safety check: stop if time exceeds reasonable limit (10 hours)
        if t > 36000:
            print("\nSimulation time exceeded 10 hours - stopping")
            break

    # Final report
    total_time_min = t / 60.0
    total_time_hr = t / 3600.0
    total_distance_km = x / 1000.0
    fraction_of_ring = x / cfg.L_RING

    print("\n" + "=" * 80)
    print("SIMULATION COMPLETE")
    print("=" * 80)
    print(f"Final velocity:       {v:,.1f} m/s")
    print(f"Target velocity:      {cfg.V_LAUNCH:,.1f} m/s")
    print(f"Achievement:          {100 * v / cfg.V_LAUNCH:.1f}%")
    print(f"Total time:           {total_time_min:.1f} minutes ({total_time_hr:.2f} hours)")
    print(f"Distance traveled:    {total_distance_km:,.1f} km")
    print(f"Fraction of ring:     {fraction_of_ring:.1%}")
    print(f"Final KE:             {KE/1e12:.3f} TJ")
    print()
    print(f"Peak thrust:          {max_values['thrust']/1e6:.2f} MN")
    print(f"Peak power:           {max_values['power']/1e9:.2f} GW")
    print(f"Peak frequency:       {max_values['frequency']:.1f} Hz")
    print(f"Peak voltage:         {max_values['voltage']/1e3:.1f} kV")
    print(f"Peak current:         {max_values['current']:.1f} A")
    print(f"Peak load angle:      {math.degrees(max_values['delta']):.1f}°")
    print(f"Peak g-load:          {max_values['g_load']:.3f} g")
    print()
    print(f"Total HTS hysteresis: {E_hyst_cum:.3e} J  ({E_hyst_cum/3.6e6:.1f} kWh)")
    print(f"Peak HTS hysteresis:  {max_values['hysteresis']/1e3:.2f} kW")
    print(f"Peak radiator width:  {max_values['radiator_width']:.4f} m")
    print()
    print(f"Voltage limited:      {voltage_limited}")
    print(f"Thrust limited:       {thrust_limited}")
    print("=" * 80)

    success = (v >= cfg.V_LAUNCH * 0.99)  # Within 1% of target

    return success


# =============================================================================
# PARAMETER DISPLAY
# =============================================================================

def print_parameters():
    """Display configuration parameters."""
    derived = cfg.calc_derived()

    print(f"\nConfiguration:")
    print(f"  Pole pitch:             {cfg.TAU_P:.1f} m")
    print(f"  Stator turns:           {cfg.N_STATOR}")
    print(f"  Air gap:                {cfg.G_GAP * 1000:.1f} mm")
    print(f"  Coil height:            {cfg.W_COIL:.2f} m")
    print("  Stator type:            Air-core")
    print(f"  Voltage limit:          {cfg.V_COIL_LIMIT / 1000:.0f} kV")
    print()
    print(f"  Sled length:            {cfg.L_SLED:,.0f} m")
    print(f"  Repeating unit length:  {derived['L_unit']:.1f} m")
    print(f"  Active fraction:        {derived['f_active']:.1%}")
    print(f"  Number of units:        {derived['n_units']}")
    print(f"  Total active area:      {derived['A_active_total']:,.1f} m²")
    print()
    print(f"  B_sled (nominal):       {cfg.B_SLED_NOMINAL:.3f} T")
    print(f"  B_sled adjustable:      {cfg.B_SLED_ADJUSTABLE}")
    print()
    print(f"  Sled hardware mass:     {derived['M_sled_hardware']:,.0f} kg")
    print(f"  Spacecraft mass:        {cfg.M_SPACECRAFT:,.0f} kg")
    print(f"  Total mass:             {derived['M_total']:,.0f} kg")
    print()
    print(f"  Target velocity:        {cfg.V_LAUNCH:,.0f} m/s")
    print(f"  Max acceleration:       {cfg.A_MAX_G:.2f} g")
    print(f"  Target KE:              {derived['KE_target'] / 1e12:.3f} TJ")
    print(f"  Max supply frequency:   {derived['f_max_supply']:.1f} Hz")
    print()


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def _param_subtitle():
    """Build parameter subtitle string from current config."""
    return (f"(\u03c4\u209a: {cfg.TAU_P:.0f} m, N: {cfg.N_STATOR}, "
            f"{cfg.N_LAYERS}\u00d7{cfg.TAPE_WIDTH*1000:.0f} mm HTS, "
            f"m: {cfg.M_SPACECRAFT/1000:,.0f} t, "
            f"v: {cfg.V_LAUNCH/1000:.1f} km/s)")

def plot_combined():
    """Create a 9-panel combined plot."""
    fig = plt.figure(figsize=(16, 16))
    gs = gridspec.GridSpec(4, 3, hspace=0.3, wspace=0.3)

    # Convert time to minutes
    t_min = [t / 60.0 for t in data_t]

    # 1. Velocity
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t_min, [v / 1000 for v in data_v], 'b-', linewidth=1)
    ax1.set_ylabel('Velocity (km/s)')
    ax1.set_xlabel('Time (min)')
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Velocity')

    # 2. Acceleration
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(t_min, [a / cfg.G_0 for a in data_a], 'r-', linewidth=1)
    ax2.set_ylabel('Acceleration (g)')
    ax2.set_xlabel('Time (min)')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Acceleration')

    # 3. Thrust
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(t_min, [F / 1e6 for F in data_F], 'purple', linewidth=1)
    ax3.set_ylabel('Thrust (MN)')
    ax3.set_xlabel('Time (min)')
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Thrust')

    # 4. Power
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(t_min, [P / 1e9 for P in data_P_mech], 'green', linewidth=1)
    ax4.set_ylabel('Power (GW)')
    ax4.set_xlabel('Time (min)')
    ax4.grid(True, alpha=0.3)
    ax4.set_title('Mechanical Power')

    # 5. Load angle
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(t_min, [math.degrees(d) for d in data_delta], 'orange', linewidth=1)
    ax5.set_ylabel('Load Angle (°)')
    ax5.set_xlabel('Time (min)')
    ax5.grid(True, alpha=0.3)
    ax5.set_title('Load Angle δ')

    # 6. Supply frequency
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.plot(t_min, data_f_supply, 'brown', linewidth=1)
    ax6.set_ylabel('Frequency (Hz)')
    ax6.set_xlabel('Time (min)')
    ax6.grid(True, alpha=0.3)
    ax6.set_title('Supply Frequency')

    # 7. Current and Voltage (dual axis)
    ax7 = fig.add_subplot(gs[2, 0])
    ax7_twin = ax7.twinx()
    ax7.plot(t_min, data_I_stator, 'b-', linewidth=1, label='Current')
    ax7_twin.plot(t_min, [V / 1000 for V in data_V_coil], 'r-', linewidth=1, label='Voltage')
    ax7.set_ylabel('Current (A)', color='b')
    ax7_twin.set_ylabel('Voltage (kV)', color='r')
    ax7.set_xlabel('Time (min)')
    ax7.grid(True, alpha=0.3)
    ax7.set_title('Stator Current & Coil Voltage')

    # 8. Magnetic fields
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.plot(t_min, data_B_stator, 'b-', linewidth=1, label='B_stator')
    ax8.plot(t_min, data_B_sled, 'r--', linewidth=1, label='B_sled')
    ax8.set_ylabel('Magnetic Field (T)')
    ax8.set_xlabel('Time (min)')
    ax8.legend()
    ax8.grid(True, alpha=0.3)
    ax8.set_title('Magnetic Fields')

    # 9. Occupant g-load
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.plot(t_min, data_a_occupant_g, 'darkred', linewidth=1)
    ax9.set_ylabel('G-load (g)')
    ax9.set_xlabel('Time (min)')
    ax9.grid(True, alpha=0.3)
    ax9.set_title('Occupant G-Load')

    # 10. HTS Hysteresis Power
    ax10 = fig.add_subplot(gs[3, 0])
    ax10.plot(t_min, [P / 1000 for P in data_P_hysteresis], '#CC79A7', linewidth=1)
    ax10.set_ylabel('HTS Hyst (kW)')
    ax10.set_xlabel('Time (min)')
    ax10.grid(True, alpha=0.3)
    ax10.set_title('HTS Hysteresis Power')

    # 11. Radiator Width
    ax11 = fig.add_subplot(gs[3, 1])
    ax11.plot(t_min, data_radiator_width, '#0072B2', linewidth=1)
    ax11.set_ylabel('Radiator Width (m)')
    ax11.set_xlabel('Time (min)')
    ax11.grid(True, alpha=0.3)
    ax11.set_title('Cryo Radiator Width')

    # 12. Cumulative Hysteresis Energy
    ax12 = fig.add_subplot(gs[3, 2])
    ax12.plot(t_min, [E / 3.6e6 for E in data_E_hyst_cumulative], '#009E73', linewidth=1)
    ax12.set_ylabel('Energy (kWh)')
    ax12.set_xlabel('Time (min)')
    ax12.grid(True, alpha=0.3)
    ax12.set_title('Cumulative Hysteresis Energy')

    # Overall title
    fig.suptitle(f"LSM Mass Driver Launch Profile\n{_param_subtitle()}",
                 fontsize=14, fontweight='bold')

    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, f"lsm_combined.{cfg.GRAPH_FORMAT}")
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight')
        print(f"Saved combined plot: {filepath}")
        plt.close(fig)
    else:
        plt.show()


def plot_individual(keyword):
    """Create individual plot based on keyword."""
    t_min = [t / 60.0 for t in data_t]

    plots = {
        'velocity': (data_v, 'Velocity (m/s)', 'Velocity', 'blue', lambda y: y),
        'accel': (data_a, 'Acceleration (g)', 'Acceleration', 'red', lambda y: y / cfg.G_0),
        'thrust': (data_F, 'Thrust (MN)', 'Thrust', 'purple', lambda y: y / 1e6),
        'power': (data_P_mech, 'Power (GW)', 'Mechanical Power', 'green', lambda y: y / 1e9),
        'delta': (data_delta, 'Load Angle (°)', 'Load Angle δ', 'orange', lambda y: math.degrees(y)),
        'frequency': (data_f_supply, 'Frequency (Hz)', 'Supply Frequency', 'brown', lambda y: y),
        'voltage': (data_V_coil, 'Voltage (kV)', 'Coil Voltage', 'red', lambda y: y / 1000),
        'current': (data_I_stator, 'Current (A)', 'Stator Current', 'blue', lambda y: y),
        'bfield': (data_B_stator, 'Magnetic Field (T)', 'Stator Field', 'navy', lambda y: y),
        'gload': (data_a_occupant_g, 'G-load (g)', 'Occupant G-Load', 'darkred', lambda y: y),
        'hysteresis': (data_P_hysteresis, 'HTS Hyst (kW)', 'HTS Hysteresis Power', '#CC79A7', lambda y: y / 1000),
        'radiator_width': (data_radiator_width, 'Radiator Width (m)', 'Cryo Radiator Width', '#0072B2', lambda y: y),
        'E_hyst': (data_E_hyst_cumulative, 'Hyst Energy (kWh)', 'Cumulative Hysteresis Energy', '#009E73', lambda y: y / 3.6e6),
    }

    if keyword not in plots:
        print(f"Unknown plot keyword: {keyword}")
        return

    data, ylabel, title, color, transform = plots[keyword]
    y_plot = [transform(y) for y in data]

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES, cfg.GRAPH_HEIGHT_INCHES))
    ax.plot(t_min, y_plot, color=color, linewidth=1.5)
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"{title}\n{_param_subtitle()}", fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, f"lsm_{keyword}.{cfg.GRAPH_FORMAT}")
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight')
        print(f"Saved plot: {filepath}")
        plt.close(fig)
    else:
        plt.show()


def save_csv():
    """Export time-series data to CSV."""
    os.makedirs("./output", exist_ok=True)
    filepath = "./output/lsm_data.csv"

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            't (s)', 'v (m/s)', 'x (m)', 'a (m/s²)', 'F (N)', 'P_mech (W)',
            'I_stator (A)', 'B_stator (T)', 'B_sled (T)', 'delta (rad)',
            'f_supply (Hz)', 'V_coil (V)', 'a_centrifugal (m/s²)',
            'a_radial_net (m/s²)', 'a_occupant_g (g)', 'KE (J)', 'x_fraction',
            'P_hysteresis (W)', 'E_hyst_cumulative (J)', 'radiator_width (m)'
        ])

        for i in range(len(data_t)):
            writer.writerow([
                data_t[i], data_v[i], data_x[i], data_a[i], data_F[i], data_P_mech[i],
                data_I_stator[i], data_B_stator[i], data_B_sled[i], data_delta[i],
                data_f_supply[i], data_V_coil[i], data_a_centrifugal[i],
                data_a_radial_net[i], data_a_occupant_g[i], data_KE[i], data_x_fraction[i],
                data_P_hysteresis[i], data_E_hyst_cumulative[i], data_radiator_width[i]
            ])

    print(f"Saved CSV data: {filepath}")


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def main():
    """Main entry point with command-line argument parsing."""
    quick_mode = False
    save_csv_flag = False
    plot_keywords = []

    # Parse arguments
    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--quick":
                quick_mode = True
            elif arg == "--save-csv":
                save_csv_flag = True
            elif arg.startswith("--v_launch="):
                cfg.V_LAUNCH = float(arg.split("=")[1])
            elif arg.startswith("--accel="):
                cfg.A_MAX_G = float(arg.split("=")[1])
            elif arg.startswith("--mass="):
                cfg.M_SPACECRAFT = float(arg.split("=")[1]) * 1000  # Convert tonnes to kg
            elif arg.startswith("--B_sled="):
                cfg.B_SLED_NOMINAL = float(arg.split("=")[1])
            elif arg.startswith("--tau_p="):
                cfg.TAU_P = float(arg.split("=")[1])
            elif arg.startswith("--N_stator="):
                cfg.N_STATOR = int(arg.split("=")[1])
            elif arg == "--adjustable-B":
                cfg.B_SLED_ADJUSTABLE = True
            elif arg.startswith("--sled-length="):
                cfg.L_SLED = float(arg.split("=")[1])
            else:
                # Assume it's a plot keyword
                plot_keywords.append(arg)

    # Recalculate derived parameters if anything changed
    cfg.calc_derived()

    # Run simulation
    success = run_lsm_simulation(quick_mode)

    # Save CSV if requested
    if save_csv_flag:
        save_csv()

    # Generate plots
    if plot_keywords:
        if "combined" in plot_keywords:
            plot_combined()
        if "all" in plot_keywords:
            for kw in ['velocity', 'accel', 'thrust', 'power', 'delta',
                      'frequency', 'voltage', 'current', 'bfield', 'gload',
                      'hysteresis', 'radiator_width', 'E_hyst']:
                plot_individual(kw)
        else:
            for kw in plot_keywords:
                if kw != "combined":
                    plot_individual(kw)


if __name__ == "__main__":
    main()
