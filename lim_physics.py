#!/usr/bin/env python3
"""
LIM Physics Module - Calculation functions for orbital ring deployment simulation

This module contains all the physics calculations for:
  - Orbital dynamics
  - Electromagnetic fields and thrust
  - Thermal calculations
  - Cryogenic system sizing
  - Radiator calculations

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7        # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.670374e-8  # Stefan-Boltzmann constant (W/m²K⁴)

# Orbital parameters at 250 km altitude
V_ORBIT = 7754.866              # Orbital velocity (m/s)
V_GROUND_STATIONARY = 483.331   # Ground-stationary velocity at 250 km (m/s)
L_RING = 41_645_813.012         # Ring circumference (m)

# Deep space
T_SPACE = 2.7                   # Deep space temperature (K)


# =============================================================================
# ORBITAL DYNAMICS
# =============================================================================

def calc_cable_velocity(v0, thrust_per_lim, dt, lim_sites, m_cable_total):
    """Update cable velocity based on thrust."""
    thrust_total = thrust_per_lim * lim_sites * 2
    return v0 + (thrust_total / m_cable_total) * dt


def calc_casing_velocity(v0, thrust_per_lim, dt, lim_sites, m_load_total):
    """Update casing velocity based on reaction force."""
    thrust_total = thrust_per_lim * lim_sites * 2
    return v0 - (thrust_total / m_load_total) * dt


def calc_relative_velocity(v_cable, v_casing):
    """Relative velocity between cable (rotor) and casing (stator)."""
    return v_cable - v_casing


def calc_deployment_progress(v_casing):
    """Percentage completion of deployment."""
    total_delta_v = V_ORBIT - V_GROUND_STATIONARY
    current_delta_v = V_ORBIT - v_casing
    return 100.0 * current_delta_v / total_delta_v


# =============================================================================
# SLIP AND FREQUENCY
# =============================================================================

def calc_slip_ratio(v_slip, v_rel):
    """Calculate slip ratio: s = v_slip / v_wave."""
    v_wave = v_slip + v_rel
    if v_wave < 1:
        v_wave = 1.0
    return v_slip / v_wave


def calc_slip_frequency(v_slip, tau_p):
    """Slip frequency in Hz."""
    return v_slip / (2 * tau_p)


def calc_supply_frequency(v_wave, tau_p):
    """Electrical supply frequency for a given wave speed."""
    return v_wave / (2 * tau_p)


# =============================================================================
# MAGNETIC FIELDS
# =============================================================================

def calc_b_field_at_plate(i_peak, n_turns, w_coil, gap):
    """Traveling wave magnetic field amplitude at the reaction plate."""
    b_single = (2 * MU0 * n_turns * i_peak / (math.pi * w_coil)) * math.atan(w_coil / (2 * gap))
    return 1.5 * b_single


def calc_b_field_in_coil(i_peak, n_turns, w_coil):
    """Peak magnetic field inside the coil."""
    return MU0 * n_turns * i_peak / w_coil


# =============================================================================
# MATERIAL PROPERTIES
# =============================================================================

def calc_resistivity(temp_K, rho_293K, alpha, material="titanium", rho_70K_alu=4.853e-9):
    """Temperature-dependent resistivity of reaction plate material."""
    rho = rho_293K * (1 + alpha * (temp_K - 293))
    if material == "aluminum":
        return max(rho, rho_70K_alu)
    return rho


def calc_skin_depth(f_slip, temp_K, rho_293K, alpha, material, v_slip_min):
    """Electromagnetic skin depth in the reaction plate."""
    rho = calc_resistivity(temp_K, rho_293K, alpha, material)
    if f_slip <= 0:
        f_slip = v_slip_min / 200  # Approximate for tau_p=100
    return math.sqrt(rho / (math.pi * MU0 * f_slip))


def calc_effective_plate_depth(f_slip, temp_K, t_plate, rho_293K, alpha, material, v_slip_min):
    """Effective conducting depth: minimum of plate thickness and skin depth."""
    delta = calc_skin_depth(f_slip, temp_K, rho_293K, alpha, material, v_slip_min)
    return min(t_plate, delta)


# =============================================================================
# EDDY CURRENT LOSSES
# =============================================================================

def calc_eddy_loss_from_thrust(v_slip, thrust):
    """Reaction plate ohmic loss consistent with thrust (per LIM)."""
    if v_slip <= 0 or thrust <= 0:
        return 0.0
    return thrust * v_slip


def calc_plate_eddy_volume(f_slip, temp_K, w_coil, l_active, t_plate, rho_293K, alpha, material, v_slip_min):
    """Volume of reaction plate participating in eddy currents."""
    depth = calc_effective_plate_depth(f_slip, temp_K, t_plate, rho_293K, alpha, material, v_slip_min)
    return w_coil * l_active * depth


# =============================================================================
# GOODNESS FACTOR (for Model 2)
# =============================================================================

def calc_goodness_factor(f_slip, temp_K, tau_p, t_plate, rho_293K, alpha, material, v_slip_min):
    """Laithwaite's goodness factor G."""
    if f_slip <= 0:
        return 0.0
    
    rho = calc_resistivity(temp_K, rho_293K, alpha, material)
    sigma = 1.0 / rho
    omega = 2 * math.pi * f_slip
    delta_eff = calc_effective_plate_depth(f_slip, temp_K, t_plate, rho_293K, alpha, material, v_slip_min)
    
    return (omega * MU0 * sigma * delta_eff * tau_p) / math.pi


def calc_slip_efficiency(slip_ratio, G):
    """Slip efficiency from goodness factor model."""
    if slip_ratio <= 0 or G <= 0:
        return 0.0
    sG = slip_ratio * G
    return (2 * sG) / (1 + sG ** 2)


# =============================================================================
# THRUST MODELS
# =============================================================================

def calc_thrust_force_max(b_field, a_lim):
    """Upper bound on thrust from magnetic pressure."""
    return (b_field ** 2) * a_lim / (2 * MU0)


def calc_loop_inductance(tau_p, w_coil):
    """Inductance of eddy current loop for narrow plate geometry."""
    return (MU0 * tau_p / math.pi) * math.log(tau_p / w_coil)


def calc_thrust_model1(f_slip, f_supply, i_peak, temp_K, params):
    """Thrust from narrow plate eddy current model."""
    if i_peak <= 0 or f_supply <= 0 or f_slip <= 0:
        return 0.0
    
    B = calc_b_field_at_plate(i_peak, params['n_turns'], params['w_coil'], params['gap'])
    v_slip = f_slip * 2 * params['tau_p']
    if v_slip <= 0:
        return 0.0
    
    rho = calc_resistivity(temp_K, params['rho_293K'], params['alpha'], params['material'])
    delta = calc_skin_depth(f_slip, temp_K, params['rho_293K'], params['alpha'], 
                           params['material'], params['v_slip_min'])
    d_eff = min(params['t_plate'], delta)
    
    EMF = (4.0 / math.pi) * v_slip * B * params['w_coil']
    
    omega = 2 * math.pi * f_slip
    R = rho * params['tau_p'] / (d_eff * params['w_coil'])
    L = calc_loop_inductance(params['tau_p'], params['w_coil'])
    X = omega * L
    Z = math.sqrt(R**2 + X**2)
    
    I_eddy = EMF / Z
    N_loops = params['l_active'] / params['tau_p']
    P = N_loops * I_eddy**2 * R
    
    F = P / v_slip
    return params['thrust_efficiency'] * F


def calc_thrust_model2(f_slip, f_supply, i_peak, temp_K, params):
    """Thrust from goodness factor model."""
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    
    b_plate = calc_b_field_at_plate(i_peak, params['n_turns'], params['w_coil'], params['gap'])
    slip = f_slip / f_supply if f_supply > 0 else 0
    slip = max(0.001, min(1.0, slip))
    
    G = calc_goodness_factor(f_slip, temp_K, params['tau_p'], params['t_plate'],
                            params['rho_293K'], params['alpha'], params['material'],
                            params['v_slip_min'])
    eta_slip = calc_slip_efficiency(slip, G)
    F_max = calc_thrust_force_max(b_plate, params['a_lim'])
    
    return params['thrust_efficiency'] * F_max * eta_slip


def calc_thrust_model3(f_slip, f_supply, i_peak, temp_K, params):
    """Thrust from slip × pressure model."""
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    
    b_plate = calc_b_field_at_plate(i_peak, params['n_turns'], params['w_coil'], params['gap'])
    slip = f_slip / f_supply if f_supply > 0 else 0
    slip = max(0.0, min(1.0, slip))
    
    F_max = calc_thrust_force_max(b_plate, params['a_lim'])
    return params['thrust_efficiency'] * slip * F_max


def calc_thrust(f_slip, f_supply, i_peak, temp_K, params, model=1):
    """Calculate thrust using the selected model."""
    if model == 1:
        return calc_thrust_model1(f_slip, f_supply, i_peak, temp_K, params)
    elif model == 2:
        return calc_thrust_model2(f_slip, f_supply, i_peak, temp_K, params)
    elif model == 3:
        return calc_thrust_model3(f_slip, f_supply, i_peak, temp_K, params)
    else:
        return calc_thrust_model1(f_slip, f_supply, i_peak, temp_K, params)


# =============================================================================
# ELECTRICAL CALCULATIONS
# =============================================================================

def calc_coil_voltage(i_peak, f_supply, b_field, a_coil, n_turns):
    """Induced voltage in LIM coil from changing flux linkage."""
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    phi_peak = b_field * a_coil
    omega = 2 * math.pi * f_supply
    return omega * n_turns * phi_peak


def calc_thrust_power(thrust, v_rel):
    """Mechanical power delivered to the cable."""
    return thrust * v_rel


# =============================================================================
# HTS HYSTERESIS LOSSES
# =============================================================================

def q_hts_loss_factor(i_peak, n_turns, w_coil, w_tape, alpha_tape):
    """Hysteresis loss factor for HTS tape."""
    b_coil = MU0 * n_turns * i_peak / w_coil
    return b_coil * i_peak * w_tape * math.sin(alpha_tape)


def q_hysteresis_norris(i_now, i_c):
    """Norris formula for hysteresis loss."""
    ii = i_now / i_c
    if ii >= 1:
        ii = 0.999
    if ii <= -1:
        ii = -0.999
    return (MU0 * i_c**2 / math.pi) * ((1 - ii) * math.log(1 - ii) + (1 + ii) * math.log(1 + ii) - ii**2)


def calc_hysteresis_power(f_supply, i_peak, params, norris=False):
    """Total hysteresis power loss per LIM."""
    if norris:
        q = q_hysteresis_norris(i_peak, params['i_c'])
    else:
        q = q_hts_loss_factor(i_peak, params['n_turns'], params['w_coil'], 
                              params['w_tape'], params['alpha_tape'])
    p_coil = q * params['l_hts_coil'] * f_supply
    return p_coil * params['lim_phases'] * params['pitch_count']


# =============================================================================
# THERMAL CALCULATIONS
# =============================================================================

def calc_effective_emissivity(em1, em2):
    """Effective emissivity between two parallel surfaces."""
    return 1 / (1/em1 + 1/em2 - 1)


def calc_radiative_heat_transfer(T_hot, area, T_cold, em1, em2):
    """Radiative heat transfer between two surfaces."""
    em_eff = calc_effective_emissivity(em1, em2)
    return em_eff * STEFAN_BOLTZMANN * area * (T_hot**4 - T_cold**4)


def calc_heatsink_area_required(p_heat, T_hot, T_cold, em1, em2):
    """Required heatsink area to radiate given power."""
    em_eff = calc_effective_emissivity(em1, em2)
    if T_hot <= T_cold:
        T_hot = T_cold + 1
    return p_heat / (em_eff * STEFAN_BOLTZMANN * (T_hot**4 - T_cold**4))


def calc_plate_temperature(v_rel, p_eddy, l_active, lim_spacing, a_lim, plate_em, em_heatsink, t_ln2_boil, v_rel_min=10):
    """Equilibrium temperature of reaction plate."""
    v_rel_adj = max(v_rel_min, v_rel)
    time_under_lim = l_active / v_rel_adj
    time_between_lims = lim_spacing / v_rel_adj
    time_cycle = time_under_lim + time_between_lims
    
    em_eff = calc_effective_emissivity(plate_em, em_heatsink)
    term = (p_eddy * time_under_lim) / (em_eff * STEFAN_BOLTZMANN * (2 * a_lim) * time_cycle)
    return (term + t_ln2_boil**4) ** 0.25


def calc_plate_temp_rise(v_rel, p_eddy, l_active, volume, density, cp, v_rel_min=10):
    """Temperature rise of plate during one LIM passage."""
    v_rel_adj = max(v_rel_min, v_rel)
    time_under_lim = l_active / v_rel_adj
    mass = volume * density
    if mass <= 0:
        return 0.0
    return (p_eddy * time_under_lim) / (mass * cp)


def calc_heat_load(p_eddy, l_active, heatsink_length, q_absorbed_per_site):
    """Total heat load including environmental absorption."""
    duty_cycle = l_active / heatsink_length
    p_average = p_eddy * duty_cycle
    return p_average + q_absorbed_per_site


# =============================================================================
# CABLE EQUILIBRIUM TEMPERATURE (NEW)
# =============================================================================

def calc_cable_equilibrium_temperature(p_eddy_per_site, lims_per_site, lim_sites, 
                                       cable_emissivity, cable_surface_area_per_m):
    """Calculate cable equilibrium temperature when eddy heat stays in cable."""
    p_eddy_total = p_eddy_per_site * lims_per_site * lim_sites
    q_in_per_m = p_eddy_total / L_RING
    
    radiative_coeff = cable_emissivity * STEFAN_BOLTZMANN * cable_surface_area_per_m
    
    if radiative_coeff <= 0:
        return 500.0  # Fallback
    
    T4 = q_in_per_m / radiative_coeff + T_SPACE**4
    T_eq = T4 ** 0.25
    
    return T_eq


def calc_cable_local_temperature(p_eddy, v_rel, l_active, lim_spacing, 
                                 plate_density, t_plate, w_plate, plate_cp, v_rel_min=10):
    """Calculate local cable temperature rise as it passes through LIM."""
    if v_rel <= 0:
        return 0.0
    
    v_rel_adj = max(v_rel_min, v_rel)
    time_under_lim = l_active / v_rel_adj
    
    mass_per_m = plate_density * t_plate * w_plate
    mass_heated = mass_per_m * l_active
    
    if mass_heated <= 0:
        return 0.0
    
    delta_T_heat = (p_eddy * time_under_lim) / (mass_heated * plate_cp)
    return delta_T_heat


# =============================================================================
# CRYOGENIC SYSTEM
# =============================================================================

def calc_cryo_cop(T_cold, T_hot, efficiency):
    """Coefficient of performance for cryogenic system."""
    if T_hot <= T_cold:
        T_hot = T_cold + 1
    cop_carnot = T_cold / (T_hot - T_cold)
    return cop_carnot * efficiency


def calc_cryo_power(heat_load, T_cold, T_hot, efficiency):
    """Electrical power required to pump heat from cold to hot side."""
    cop = calc_cryo_cop(T_cold, T_hot, efficiency)
    if cop <= 0:
        return heat_load * 100  # Fallback
    return heat_load / cop


def calc_cryo_heat_load(p_eddy, p_hyst, lims_per_site, q_absorbed_per_site,
                        eddy_to_cable, q_coil_environment):
    """Calculate heat load that must go through cryogenic system."""
    if eddy_to_cable:
        # Only hysteresis in coils + heat leak through MLI to coils
        q_cold = lims_per_site * p_hyst + q_coil_environment
    else:
        # Everything goes through cryo
        total_heat = lims_per_site * (p_eddy + p_hyst) + q_absorbed_per_site
        q_cold = total_heat
    
    return q_cold


# =============================================================================
# RADIATOR CALCULATIONS (NEW)
# =============================================================================

def calc_radiator_width(q_cold, T_cold, T_hot, efficiency, em_heatsink, radiator_length):
    """Calculate required radiator width to reject cryogenic heat."""
    if q_cold <= 0:
        return 0.0
    
    if T_hot <= T_cold:
        T_hot = T_cold + 1
    cop_carnot = T_cold / (T_hot - T_cold)
    cop_real = cop_carnot * efficiency
    
    if cop_real <= 0:
        cop_real = 0.01
    
    # Heat rejected at hot side
    q_reject = q_cold * (1 + 1/cop_real)
    
    # Radiator power density (one-sided)
    p_per_m2 = em_heatsink * STEFAN_BOLTZMANN * (T_hot**4 - T_SPACE**4)
    
    if p_per_m2 <= 0:
        return 1000.0  # Fallback
    
    area_required = q_reject / p_per_m2
    width = area_required / radiator_length
    
    return width


def calc_radiator_width_direct(p_heat, T_radiator, em_heatsink, radiator_length):
    """Calculate radiator width for direct (non-cryogenic) heat rejection."""
    if p_heat <= 0:
        return 0.0
    
    p_per_m2 = em_heatsink * STEFAN_BOLTZMANN * (T_radiator**4 - T_SPACE**4)
    
    if p_per_m2 <= 0:
        return 1000.0
    
    area_required = p_heat / p_per_m2
    width = area_required / radiator_length
    
    return width


# =============================================================================
# KINETIC ENERGY
# =============================================================================

def get_ke(mass_m, v_init, v_final):
    """Calculate kinetic energy change."""
    return 0.5 * mass_m * (v_final**2 - v_init**2) * L_RING
