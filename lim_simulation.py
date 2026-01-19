#!/usr/bin/env python3
"""
LIM Deployment Simulation - Main simulation loop and plotting

This is the main entry point for the orbital ring deployment simulation.

Usage:
    python lim_simulation.py [--model=N] [--thermal=cable|cryo] [graph_options]
    
    --model=1  Narrow plate eddy current (default)
    --model=2  Goodness factor model (Laithwaite)
    --model=3  Slip × pressure (theoretical maximum)
    
    --thermal=cable  Eddy heat goes to cable
    --thermal=cryo   Eddy heat goes to cryogenic system (default)

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import datetime
import math
import matplotlib.pyplot as plt

try:
    import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False
    print("Warning: tabulate not installed. Tables will be simplified.")

# Import configuration and physics modules
import lim_config as cfg
import lim_physics as phys


# =============================================================================
# DATA COLLECTION
# =============================================================================

data_current = []
data_voltage = []
data_v_slip = []
data_slip_ratio = []
data_f_slip = []
data_thrust = []
data_thrust_power = []
data_v_rel = []
data_b_field = []
data_p_eddy = []
data_skin_depth_eff = []
data_skin_depth_calc = []
data_p_cryo = []
data_p_hyst = []
data_p_heat = []
data_p_lim = []
data_p_site = []
data_temp_plate = []
data_E_site = []
data_E_total = []
data_cable_temp = []
data_radiator_width = []
data_q_cryo_cold = []


def clear_data():
    """Clear all data lists for a fresh run."""
    global data_current, data_voltage, data_v_slip, data_slip_ratio, data_f_slip
    global data_thrust, data_thrust_power, data_v_rel, data_b_field, data_p_eddy
    global data_skin_depth_eff, data_skin_depth_calc, data_p_cryo, data_p_hyst
    global data_p_heat, data_p_lim, data_p_site, data_temp_plate
    global data_E_site, data_E_total, data_cable_temp, data_radiator_width, data_q_cryo_cold
    
    data_current = []
    data_voltage = []
    data_v_slip = []
    data_slip_ratio = []
    data_f_slip = []
    data_thrust = []
    data_thrust_power = []
    data_v_rel = []
    data_b_field = []
    data_p_eddy = []
    data_skin_depth_eff = []
    data_skin_depth_calc = []
    data_p_cryo = []
    data_p_hyst = []
    data_p_heat = []
    data_p_lim = []
    data_p_site = []
    data_temp_plate = []
    data_E_site = []
    data_E_total = []
    data_cable_temp = []
    data_radiator_width = []
    data_q_cryo_cold = []


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def make_month_ticks(data_list, total_time):
    """Generate x-axis tick positions and labels in months."""
    months = total_time / cfg.MONTH
    num_ticks = int(months) + 1
    n_points = len(data_list)
    
    tick_positions = [n_points * i / max(months, 1) for i in range(num_ticks)]
    tick_labels = [str(i) for i in range(num_ticks)]
    return tick_positions, tick_labels


def annotate_final_value(data_list, unit=""):
    """Add annotation showing final value on plot."""
    fmt = ".3f"
    if not data_list:
        return
    x = len(data_list)
    y = data_list[-1]
    
    if y > 1e18:
        unit = "E" + unit
        y = y / 1e18
    elif y > 1e15:
        unit = "P" + unit
        y = y / 1e15
    elif y > 1e12:
        unit = "T" + unit
        y = y / 1e12
    elif y > 1e9:
        unit = "G" + unit
        y = y / 1e9
    elif y > 1e6:
        unit = "M" + unit
        y = y / 1e6
    elif y > 1e3:
        unit = "k" + unit
        y = y / 1e3
    
    if y > 100:
        fmt = ".1f"
    elif y > 10:
        fmt = ".2f"
    
    print(f"{y:{fmt}} {unit}".strip())
    label = f"{y:{fmt}} {unit}".strip()
    plt.annotate(
        label, xy=(x, y), xytext=(-90, 30),
        textcoords="offset points", fontsize=25, ha="left",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0)
    )


# =============================================================================
# MAIN SIMULATION LOOP
# =============================================================================

def run_deployment_simulation(v_slip_init, i_peak_init, thrust_model=1, eddy_to_cable=False):
    """Run the deployment simulation.
    
    Args:
        v_slip_init: Initial slip velocity (m/s)
        i_peak_init: Initial peak current (A)
        thrust_model: 1, 2, or 3
        eddy_to_cable: If True, eddy heat stays in cable
    
    Returns:
        (param_string, success, total_time)
    """
    clear_data()
    
    # Get physics parameters
    params = cfg.get_physics_params()
    
    # Print configuration
    print(f"VOLTS_MAX: {cfg.VOLTS_MAX/1000:.0f} kV")
    print(f"THRUST_MODEL: {thrust_model}", end=" ")
    model_names = {
        1: "(Narrow plate eddy current - geometry corrected)",
        2: "(Goodness factor - Laithwaite)",
        3: "(Slip × pressure - theoretical maximum)"
    }
    print(model_names.get(thrust_model, "(Unknown)"))
    print(f"T_PLATE: {cfg.T_PLATE * 1000:.1f} mm")
    print(f"EDDY_HEAT_TO_CABLE: {eddy_to_cable}")
    if eddy_to_cable:
        print("  → Eddy current heat stays in cable (equilibrium temperature mode)")
    else:
        print("  → Eddy current heat goes to cryogenic system")

    # Initialize max/min tracking
    max_min = {
        "volts_max": [0, 0, "", 0],
        "current_max": [0, 0, "", 0],
        "v_slip_max": [0, 0, "", 0],
        "p_cryo_max": [0, 0, "", 0],
        "p_eddy_max": [0, 0, "", 0],
        "p_hyst_max": [0, 0, "", 0],
        "p_heat_max": [0, 0, "", 0],
        "f_slip_max": [0, 0, "", 0],
        "f_supply_max": [0, 0, "", 0],
        "b_field_max": [0, 0, "", 0],
        "thrust_max": [0, 0, "", 0],
        "p_thrust_max": [0, 0, "", 0],
        "volts_lim_max": [0, 0, "", 0],
        "delta_t_max": [0, 0, "", 0],
        "plate_t_max": [0, 0, "", 0],
        "skin_d_max": [0, 0, "", 0],
        "skin_d_min": [100, 0, "", 0],
        "min_heatsink_area": [0, 0, "", 0],
        "min_cryo_area": [0, 0, "", 0],
        "p_heat_load": [0, 0, "", 0],
        "plate_temp_out": [0, 0, "", 0],
        "lim_power_max": [0, 0, "", 0],
        "site_power_max": [0, 0, "", 0],
        "cable_temp_max": [0, 0, "", 0],
        "radiator_width_max": [0, 0, "", 0],
        "q_cryo_cold_max": [0, 0, "", 0],
    }
    exit_msg = "PASSED"

    # Initialize state
    dt = cfg.DT1
    v_cable = cfg.V_ORBIT
    v_casing = cfg.V_ORBIT
    time = 0
    v_slip = v_slip_init
    i_peak = i_peak_init
    plate_temp = 70.0
    cable_temp = 70.0
    heatsink_area = 0

    # Controller state
    v_slip_cmd = v_slip
    i_peak_prev = i_peak
    site_power = 0.0

    # Controller tuning
    CTRL_TAU = 15 * cfg.DAY
    DV_SLIP_MAX = 0.0001
    DI_MAX = 0.1

    # Energy accumulators
    E_site = 0.0
    E_total = 0.0

    # Tracking
    sample_time = 0
    power_track = []
    month_count = 0

    param_str = f"τp={cfg.TAU_P}m, N={cfg.N_TURNS}, Gap={cfg.GAP*1000:.0f}mm, Spacing={cfg.LIM_SPACING:.0f}m, HTS width={cfg.HTS_TAPE_WIDTH_MM:.0f}mm, I_max={cfg.I_TARGET:.0f}A"

    # Main loop
    while v_casing > cfg.V_GROUND_STATIONARY:
        
        # Controller iterations for self-consistency
        for _ in range(10):
            v_rel = phys.calc_relative_velocity(v_cable, v_casing)

            # Predictive power control
            power_margin = 1.0 - (site_power / cfg.MAX_SITE_POWER)
            MARGIN_COMFORTABLE = 0.15
            MARGIN_TIGHT = 0.05

            if power_margin >= MARGIN_COMFORTABLE:
                s_target = cfg.SLIP_RATIO_NORMAL
            elif power_margin <= MARGIN_TIGHT:
                s_target = cfg.SLIP_RATIO_REDUCED
            else:
                t = (power_margin - MARGIN_TIGHT) / (MARGIN_COMFORTABLE - MARGIN_TIGHT)
                s_target = cfg.SLIP_RATIO_REDUCED + t * (cfg.SLIP_RATIO_NORMAL - cfg.SLIP_RATIO_REDUCED)

            # Current headroom factor
            current_ratio = i_peak / cfg.I_TARGET
            if current_ratio < 0.7:
                blend = min(1.0, max(0.0, (0.7 - current_ratio) / 0.2))
                s_target = s_target + blend * (cfg.SLIP_RATIO_REDUCED - s_target)

            s_target = max(1e-6, min(0.999, s_target))

            # Calculate target slip velocity
            v_slip_target = (s_target / (1.0 - s_target)) * max(v_rel, 0.0)
            v_slip_target = max(cfg.V_SLIP_MIN, min(cfg.V_SLIP_MAX, v_slip_target))

            # Smooth actuator response
            alpha = 1.0 - math.exp(-dt / max(CTRL_TAU, 1e-9))
            v_slip_cmd = v_slip_cmd + alpha * (v_slip_target - v_slip_cmd)
            max_step = DV_SLIP_MAX * dt
            v_slip_cmd = max(v_slip - max_step, min(v_slip + max_step, v_slip_cmd))
            v_slip = v_slip_cmd

            # Calculate operating point
            v_wave = v_rel + v_slip
            slip = phys.calc_slip_ratio(v_slip, v_rel)
            f_slip = phys.calc_slip_frequency(v_slip, cfg.TAU_P)
            f_supply = phys.calc_supply_frequency(v_wave, cfg.TAU_P)

            b_plate = phys.calc_b_field_at_plate(i_peak, cfg.N_TURNS, cfg.W_COIL, cfg.GAP)
            b_coil = phys.calc_b_field_in_coil(i_peak, cfg.N_TURNS, cfg.W_COIL)

            # Use cable temp for thrust calculation when eddy heat goes to cable
            calc_temp = cable_temp if eddy_to_cable else plate_temp
            thrust = phys.calc_thrust(f_slip, f_supply, i_peak, calc_temp, params, thrust_model)
            p_thrust = phys.calc_thrust_power(thrust, v_rel)

            volts = phys.calc_coil_voltage(i_peak, f_supply, b_coil, cfg.A_COIL, cfg.N_TURNS)

            p_eddy = phys.calc_eddy_loss_from_thrust(v_slip, thrust)
            p_hyst = phys.calc_hysteresis_power(f_supply, i_peak, params, cfg.NORRIS_HYSTERESIS)
            p_heat = p_eddy + p_hyst

            # Temperature calculations depend on thermal mode
            if eddy_to_cable:
                cable_temp = phys.calc_cable_equilibrium_temperature(
                    p_eddy, cfg.LIMS_PER_SITE, cfg.LIM_SITES,
                    cfg.CABLE_EMISSIVITY, cfg.CABLE_SURFACE_AREA_PER_M
                )
                temp_plate_avg = cable_temp
                delta_temp = phys.calc_cable_local_temperature(
                    p_eddy, v_rel, cfg.L_ACTIVE, cfg.LIM_SPACING,
                    cfg.PLATE_DENSITY, cfg.T_PLATE, cfg.W_PLATE, cfg.PLATE_CP
                )
                plate_temp = temp_plate_avg + delta_temp
            else:
                temp_plate_avg = phys.calc_plate_temperature(
                    v_rel, p_eddy, cfg.L_ACTIVE, cfg.LIM_SPACING,
                    cfg.A_LIM, cfg.PLATE_EM, cfg.EM_HEATSINK, cfg.T_LN2_BOIL
                )
                volume = phys.calc_plate_eddy_volume(
                    f_slip, plate_temp, cfg.W_COIL, cfg.L_ACTIVE, cfg.T_PLATE,
                    cfg.PLATE_RHO_293K, cfg.PLATE_ALPHA, cfg.PLATE_MATERIAL, cfg.V_SLIP_MIN
                )
                delta_temp = phys.calc_plate_temp_rise(
                    v_rel, p_eddy, cfg.L_ACTIVE, volume, cfg.PLATE_DENSITY, cfg.PLATE_CP
                )
                plate_temp = temp_plate_avg + delta_temp
                cable_temp = plate_temp

            # Calculate cryogenic heat load based on thermal mode
            q_cryo_cold = phys.calc_cryo_heat_load(
                p_eddy, p_hyst, cfg.LIMS_PER_SITE, cfg.Q_ABSORBED_PER_SITE,
                eddy_to_cable, cfg.Q_COIL_ENVIRONMENT
            )

            # Calculate radiator width
            radiator_width = phys.calc_radiator_width(
                q_cryo_cold, cfg.T_LN2_SUPPLY, cfg.T_RADIATOR_HOT,
                cfg.CRYO_EFF, cfg.EM_HEATSINK, cfg.LIM_SPACING
            )

            heat_load = phys.calc_heat_load(
                p_eddy, cfg.L_ACTIVE, cfg.HEATSINK_LENGTH, cfg.Q_ABSORBED_PER_SITE
            )

            if v_rel > 0:
                p_cryo = phys.calc_cryo_power(
                    q_cryo_cold, cfg.T_LN2_BOIL, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF
                )
            else:
                p_cryo = 0.0

            p_lim = (p_thrust + p_eddy + p_hyst) / max(cfg.LIM_EFF, 1e-6)
            site_power = (cfg.LIMS_PER_SITE * p_lim + p_cryo) / max(cfg.INV_EFF, 1e-6)

            # Enforce limits
            changed = False

            if volts > cfg.VOLTS_MAX:
                scale = cfg.VOLTS_MAX / max(volts, 1.0)
                i_peak = max(cfg.I_MIN, i_peak * max(0.0, min(1.0, scale)))
                changed = True

            if temp_plate_avg > cfg.T_MAX_PLATE or (temp_plate_avg + delta_temp) > cfg.T_MAX_PLATE:
                i_peak = max(cfg.I_MIN, i_peak * 0.95)
                v_slip = max(cfg.V_SLIP_MIN, v_slip * 0.95)
                changed = True

            if site_power > cfg.MAX_SITE_POWER:
                scale = math.sqrt(cfg.MAX_SITE_POWER / max(site_power, 1.0))
                i_new = max(cfg.I_MIN, i_peak * max(0.0, min(1.0, scale)))
                if i_new < i_peak:
                    i_peak = i_new
                    changed = True
                if i_peak <= cfg.I_MIN * 1.0001 and v_slip > cfg.V_SLIP_MIN:
                    v_slip = max(cfg.V_SLIP_MIN, v_slip * 0.95)
                    changed = True

            # Current slew limit
            i_step = DI_MAX * dt
            i_peak = max(i_peak_prev - i_step, min(i_peak_prev + i_step, i_peak))
            i_peak_prev = i_peak

            heatsink_area = phys.calc_heatsink_area_required(
                p_heat, plate_temp, cfg.T_LN2_BOIL, cfg.PLATE_EM, cfg.EM_HEATSINK
            )
            min_cryo_area = phys.calc_radiative_heat_transfer(
                plate_temp, cfg.LIM_SPACING * cfg.W_COIL, cfg.T_LN2_BOIL,
                cfg.PLATE_EM, cfg.EM_HEATSINK
            )

            # Ramp up if under limits
            if (not changed and i_peak < cfg.I_TARGET
                and site_power < cfg.POWER_HEADROOM * cfg.MAX_SITE_POWER
                and volts < cfg.POWER_HEADROOM * cfg.VOLTS_MAX
                and temp_plate_avg < cfg.POWER_HEADROOM * cfg.T_MAX_PLATE):
                i_peak = min(cfg.I_TARGET, i_peak * cfg.CURRENT_UPRATE)
                changed = True

            if not changed:
                break

        # Apply acceleration
        v_cable = phys.calc_cable_velocity(v_cable, thrust, dt, cfg.LIM_SITES, cfg.M_CABLE_TOTAL)
        v_casing = phys.calc_casing_velocity(v_casing, thrust, dt, cfg.LIM_SITES, cfg.M_LOAD_TOTAL)

        # Accumulate energy
        E_site += p_thrust * cfg.LIMS_PER_SITE * dt
        E_total += p_thrust * cfg.LIM_SITES * cfg.LIMS_PER_SITE * dt

        # Startup ramp
        if time < cfg.HR:
            if (temp_plate_avg < cfg.T_MAX_PLATE * 0.8
                and volts < cfg.VOLTS_MAX * 0.8
                and site_power < cfg.MAX_SITE_POWER * 0.8):
                i_peak += (cfg.I_TARGET - i_peak) * 0.01

        # Data collection
        skin_eff = phys.calc_effective_plate_depth(
            f_slip, plate_temp, cfg.T_PLATE, cfg.PLATE_RHO_293K,
            cfg.PLATE_ALPHA, cfg.PLATE_MATERIAL, cfg.V_SLIP_MIN
        )
        skin_calc = phys.calc_skin_depth(
            f_slip, plate_temp, cfg.PLATE_RHO_293K,
            cfg.PLATE_ALPHA, cfg.PLATE_MATERIAL, cfg.V_SLIP_MIN
        )

        if time > sample_time and time < cfg.SAMPLE_TIME_MAX:
            sample_time += cfg.SKIP
            data_current.append(i_peak)
            data_voltage.append(volts)
            data_v_slip.append(v_slip)
            data_slip_ratio.append(slip * 100)
            data_f_slip.append(f_slip)
            data_thrust.append(thrust * 2)
            data_thrust_power.append(p_thrust * 2)
            data_v_rel.append(v_rel)
            data_b_field.append(b_plate)
            data_p_eddy.append(p_eddy * 2)
            data_skin_depth_eff.append(skin_eff * 1000)
            data_skin_depth_calc.append(skin_calc * 1000)
            data_p_cryo.append(p_cryo)
            data_p_hyst.append(p_hyst * 2)
            data_p_heat.append(p_heat * 2)
            data_p_lim.append(p_lim)
            data_p_site.append(site_power)
            data_temp_plate.append(temp_plate_avg)
            data_E_site.append(E_site)
            data_E_total.append(E_total)
            data_cable_temp.append(cable_temp)
            data_radiator_width.append(radiator_width)
            data_q_cryo_cold.append(q_cryo_cold)

        # Monthly progress display
        if time > month_count * cfg.MONTH:
            progress = phys.calc_deployment_progress(v_casing)
            if month_count == 0:
                print("\nMonth | Progress | Voltage | Current | Thrust | Site Power | Radiator W")
                print("-" * 75)
            print(f"{month_count:5} | {progress:6.1f}% | {volts:7.0f} V | {i_peak:5.1f} A | {thrust:6.0f} N | {site_power/1e6:5.2f} MW | {radiator_width:6.2f} m")

            power_track.append([
                f"{month_count} mth", f"{round(v_casing, 1)} m/s", f"{round(i_peak, 1)} A",
                f"{round(volts/1e3, 2)} kV", f"{round(v_slip, 1)} m/s", f"{round(slip*100, 1)}%",
                f"{round(f_supply, 1)} Hz", f"{round(p_eddy/1e3, 2)} kW", f"{round(p_hyst/1e3, 2)} kW",
                f"{round(thrust, 1)} N", f"{round(p_thrust/1e6, 3)} MW",
                f"{round(p_cryo/1e6, 3)} MW", f"{round(site_power/1e6, 3)} MW"
            ])
            month_count += 1

        # Check for failures
        if volts > cfg.VOLTS_MAX * 1.1:
            exit_msg = f"FAIL over voltage limit: {volts:.0f} V"
            print(f"\nFAILURE: Voltage exceeded limit ({volts:.0f} V)")
            return (exit_msg, False, time)

        if site_power > cfg.MAX_SITE_POWER * 1.1:
            exit_msg = f"FAIL max site power exceeded: {site_power:.0f} W"
            print(f"\nFAILURE: Site power exceeded limit ({site_power/1e6:.1f} MW)")
            return (exit_msg, False, time)

        if cfg.TAU_P * cfg.PITCH_COUNT > cfg.LIM_SPACING:
            exit_msg = f"FAIL LIM longer than spacing"
            print(f"\nFAILURE: LIM length ({cfg.TAU_P * cfg.PITCH_COUNT}m) > spacing ({cfg.LIM_SPACING}m)")
            return (exit_msg, False, time)

        # Collect min-max data
        if v_rel != 0 and time > cfg.DAY:
            if max_min["volts_max"][3] < volts:
                max_min["volts_max"] = [round(volts, 2), time, "V", volts]
            if max_min["current_max"][3] < i_peak:
                max_min["current_max"] = [round(i_peak, 2), time, "A", i_peak]
            if max_min["v_slip_max"][3] < v_slip:
                max_min["v_slip_max"] = [round(v_slip, 6), time, "m/s", v_slip]
            if max_min["f_slip_max"][3] < f_slip:
                max_min["f_slip_max"] = [round(f_slip, 6), time, "Hz", f_slip]
            if max_min["f_supply_max"][3] < f_supply:
                max_min["f_supply_max"] = [round(f_supply, 6), time, "Hz", f_supply]
            if max_min["b_field_max"][3] < b_plate:
                max_min["b_field_max"] = [round(b_plate*1e3, 6), time, "mT", b_plate]
            if max_min["thrust_max"][3] < thrust:
                max_min["thrust_max"] = [round(thrust, 2), time, "N", thrust]
            if max_min["p_thrust_max"][3] < p_thrust:
                max_min["p_thrust_max"] = [round(p_thrust/1e6, 3), time, "MW", p_thrust]
            if max_min["p_eddy_max"][3] < p_eddy:
                max_min["p_eddy_max"] = [round(p_eddy/1e3, 2), time, "kW", p_eddy]
            if max_min["p_hyst_max"][3] < p_hyst:
                max_min["p_hyst_max"] = [round(p_hyst/1e3, 2), time, "kW", p_hyst]
            if max_min["p_heat_max"][3] < p_heat:
                max_min["p_heat_max"] = [round(p_heat/1e3, 2), time, "kW", p_heat]
            if max_min["plate_t_max"][3] < temp_plate_avg:
                max_min["plate_t_max"] = [round(temp_plate_avg, 2), time, "K", temp_plate_avg]
            if max_min["delta_t_max"][3] < delta_temp:
                max_min["delta_t_max"] = [round(delta_temp, 6), time, "K", delta_temp]
            if max_min["lim_power_max"][3] < p_lim:
                max_min["lim_power_max"] = [round(p_lim/1e6, 3), time, "MW", p_lim]
            if max_min["site_power_max"][3] < site_power:
                max_min["site_power_max"] = [round(site_power/1e6, 3), time, "MW", site_power]
            if max_min["p_cryo_max"][3] < p_cryo:
                max_min["p_cryo_max"] = [round(p_cryo/1e6, 2), time, "MW", p_cryo]
            if max_min["cable_temp_max"][3] < cable_temp:
                max_min["cable_temp_max"] = [round(cable_temp, 2), time, "K", cable_temp]
            if max_min["radiator_width_max"][3] < radiator_width:
                max_min["radiator_width_max"] = [round(radiator_width, 2), time, "m", radiator_width]
            if max_min["q_cryo_cold_max"][3] < q_cryo_cold:
                max_min["q_cryo_cold_max"] = [round(q_cryo_cold/1e3, 2), time, "kW", q_cryo_cold]

        # Time stepping
        if time < cfg.DAY:
            dt = cfg.DT1
        elif time < cfg.SAMPLE_TIME_MAX:
            dt = cfg.DT2
        else:
            dt = cfg.DT3

        time += dt

        if time > 100 * cfg.YR:
            print("\nWARNING: Deployment time exceeded 100 years")
            break

    # Final entry in power track
    power_track.append([
        "final", f"{round(v_casing, 1)} m/s", f"{round(i_peak, 1)} A",
        f"{round(volts/1e3, 2)} kV", f"{round(v_slip, 1)} m/s", f"{round(slip*100, 1)}%",
        f"{round(f_supply, 1)} Hz", f"{round(p_eddy/1e3, 2)} kW", f"{round(p_hyst/1e3, 2)} kW",
        f"{round(thrust, 1)} N", f"{round(p_thrust/1e6, 3)} MW",
        f"{round(p_cryo/1e6, 3)} MW", f"{round(site_power/1e6, 3)} MW"
    ])

    # Sanity check for KE
    ke_cable = phys.get_ke(cfg.M_CABLE_M, cfg.V_ORBIT, v_cable)
    ke_load = phys.get_ke(cfg.M_LOAD_M, cfg.V_ORBIT, v_casing)
    ke_added = ke_cable + ke_load
    ke_diff = ke_added - E_total

    # Final results
    print("\n" + "=" * 70)
    print("DEPLOYMENT COMPLETE")
    print("=" * 70)
    print(f"Deployment time: {time/cfg.DAY:.1f} days ({time/cfg.YR:.2f} years)")
    print(f"Final cable velocity: {v_cable:.2f} m/s")
    print(f"Final casing velocity: {v_casing:.2f} m/s")
    print(f"Energy per site: {E_site/1e12:.6f} TJ")
    print(f"Total energy delivered: {E_total/1e18:.8f} EJ")
    print(f"Total energy from 1/2 m v^2: {ke_added/1e18:.8f} EJ")
    print(f"These values differ by: {ke_diff/1e18:.10f} EJ")
    print(f"\nThermal Mode: {'Cable absorbs eddy heat' if eddy_to_cable else 'Cryo handles all heat'}")
    print(f"Max cable temperature: {max_min['cable_temp_max'][0]:.1f} K")
    print(f"Max radiator width: {max_min['radiator_width_max'][0]:.2f} m")
    print(f"Max cryo cold-side load: {max_min['q_cryo_cold_max'][0]:.1f} kW")
    print("=" * 70)

    # Print max/min data
    print("\n" + "-" * 70)
    for key, value in max_min.items():
        if value[1] < 60:
            t = f"{value[1]:8.2f} seconds"
        elif value[1] < cfg.HR:
            t = f"{value[1]/60:8.2f} minutes"
        elif value[1] < cfg.DAY:
            t = f"{value[1]/cfg.HR:8.2f} hours"
        elif value[1] < cfg.MONTH:
            t = f"{value[1]/cfg.DAY:8.2f} days"
        elif value[1] < cfg.YR:
            t = f"{value[1]/cfg.MONTH:8.2f} months"
        else:
            t = f"{value[1]/cfg.YR:8.2f} years"
        print(f"\t{key:18} {value[0]:18.2f} {value[2]:8} {t}")
    print("-" * 70)

    # Print power track table
    if HAS_TABULATE:
        print("\n")
        print(tabulate.tabulate(power_track, headers=[
            "Time", "V_Shell", "Amps", "Volts", "V_Slip", "Slip",
            "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "Cryo Power", "Site Power"
        ]))
    
    print("-" * 70)

    # Write to file
    if cfg.WRITE_FILE:
        try:
            import os
            os.makedirs("./output", exist_ok=True)
            with open(f"./output/_orbital_ring_model{thrust_model}.txt", "a") as file:
                timestamp = datetime.datetime.now()
                file.write(f"\n{'+'*70}\nTimestamp: {timestamp}\n")
                file.write(f"Deployment time: {time/cfg.DAY:.1f} days ({time/cfg.YR:.2f} years)\n")
                file.write(f"Thermal Mode: {'Cable' if eddy_to_cable else 'Cryo'}\n")
                file.write(f"Max cable temp: {max_min['cable_temp_max'][0]:.1f} K\n")
                file.write(f"Max radiator width: {max_min['radiator_width_max'][0]:.2f} m\n")
        except Exception as e:
            print(f"Warning: Could not write to file: {e}")

    return (param_str, True, time)


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_results(show_graphs, param_str, total_time):
    """Generate requested plots."""
    
    plots = {
        "current": (data_current, "Current", "A", "blue"),
        "volts": (data_voltage, "Voltage", "V", "red"),
        "v_slip": (data_v_slip, "Slip Velocity", "m/s", "green"),
        "thrust": (data_thrust, "Thrust per site", "N", "purple"),
        "p_eddy": (data_p_eddy, "Eddy Current Losses per site", "W", "darkblue"),
        "v_rel": (data_v_rel, "Relative Velocity", "m/s", "olive"),
        "f_slip": (data_f_slip, "Slip Frequency", "Hz", "orange"),
        "slip": (data_slip_ratio, "Slip Ratio", "%", "cyan"),
        "p_thrust": (data_thrust_power, "Thrust Power per site", "W", "navy"),
        "b_peak": (data_b_field, "Magnetic Field", "T", "brown"),
        "hyst": (data_p_hyst, "Hysteresis Losses per site", "W", "brown"),
        "p_heat": (data_p_heat, "Heat Losses per site", "W", "brown"),
        "cryo": (data_p_cryo, "Cryogenic Power per site", "W", "teal"),
        "power": (data_p_site, "Site Power", "W", "darkslategray"),
        "lim_power": (data_p_lim, "LIM Power", "W", "darkslategray"),
        "plate_temp": (data_temp_plate, "Plate Temperature", "K", "lime"),
        "ke_site": (data_E_site, "Site Kinetic Energy", "J", "green"),
        "ke_all": (data_E_total, "Total Kinetic Energy", "J", "darkgreen"),
        "skin": (data_skin_depth_eff, "Effective Skin Depth", "mm", "magenta"),
        "skin_calc": (data_skin_depth_calc, "Calculated Skin Depth", "mm", "darkmagenta"),
        "cable_temp": (data_cable_temp, "Cable Equilibrium Temperature", "K", "red"),
        "radiator_width": (data_radiator_width, "Required Radiator Width", "m", "blue"),
        "q_cryo": (data_q_cryo_cold, "Cryo Cold-Side Heat Load", "W", "cyan"),
    }

    time_str = f"{total_time/cfg.DAY:.1f} days ({total_time/cfg.YR:.2f} years)"
    month_str = f"{total_time/cfg.MONTH:.1f}"

    for name, (data, label, unit, color) in plots.items():
        if name in show_graphs or "all" in show_graphs:
            if not data:
                continue

            tick_pos, tick_labels = make_month_ticks(data, total_time)

            plt.figure(figsize=(14, 7))
            plt.scatter(range(len(data)), data, c=color, s=1)
            plt.xlabel(f"{month_str} Months", fontsize=24)
            plt.ylabel(f"{label} ({unit})", fontsize=24)
            if "fulldata" in show_graphs:
                plt.title(f"{label} - {param_str}", fontsize=22)
            elif "timeonly" in show_graphs:
                plt.title(f"{label} - {time_str}", fontsize=24)
            else:
                plt.title(f"{label}", fontsize=24)
            plt.xticks(tick_pos, tick_labels)
            annotate_final_value(data, unit=unit)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    """Main entry point."""
    thrust_model = cfg.THRUST_MODEL
    eddy_to_cable = cfg.EDDY_HEAT_TO_CABLE
    show_graphs = []

    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            print("""
Graph Options:
  all          Show all graphs
  current      Coil current (A)
  volts        Induced voltage (V)
  v_slip       Slip velocity (m/s)
  thrust       Thrust force (N)
  p_thrust     Thrust power (W)
  p_eddy       Eddy current losses (W)
  power        Total site power (W)
  plate_temp   Reaction plate temperature (K)
  cable_temp   Cable equilibrium temperature (K)
  radiator_width  Required radiator width (m)
  q_cryo       Cryo cold-side heat load (W)
  skin         Skin depth (mm)
  slip         Slip ratio (%)
  f_slip       Slip frequency (Hz)
  v_rel        Relative velocity (m/s)
  hyst         Hysteresis losses (W)
  p_heat       Heat losses (W)
  cryo         Cryogenic power (W)
  ke_site      Site kinetic energy (J)
  ke_all       Total kinetic energy (J)

Examples:
  python lim_simulation.py --model=1 thrust power
  python lim_simulation.py --model=2 --thermal=cable all
  python lim_simulation.py --thermal=cryo cable_temp radiator_width
            """)
            return

        for arg in sys.argv[1:]:
            if arg.startswith("--model="):
                try:
                    model = int(arg.split("=")[1])
                    if model in [1, 2, 3]:
                        thrust_model = model
                    else:
                        print(f"Invalid model {model}, using default (1)")
                except ValueError:
                    print(f"Invalid model argument: {arg}")
            elif arg.startswith("--thermal="):
                mode = arg.split("=")[1].lower()
                if mode == "cable":
                    eddy_to_cable = True
                    print(f"Thermal mode set to {mode}: EDDY_HEAT_TO_CABLE = True")
                elif mode == "cryo":
                    eddy_to_cable = False
                    print(f"Thermal mode set to {mode}: EDDY_HEAT_TO_CABLE = False")
                else:
                    print(f"Invalid thermal mode {mode}, using default (cryo)")
            else:
                show_graphs.append(arg)

    # Print configuration
    cfg.print_parameters()
    print(f"Using Thrust Model {thrust_model}")
    print(f"Thermal Mode: {'cable' if eddy_to_cable else 'cryo'}")

    # Run simulation
    v_slip = cfg.V_SLIP_MAX
    i_peak = cfg.I_MIN

    param_str, success, total_time = run_deployment_simulation(
        v_slip, i_peak, thrust_model, eddy_to_cable
    )

    # Plot results if requested
    if success and show_graphs:
        plot_results(show_graphs, param_str, total_time)


if __name__ == "__main__":
    main()
