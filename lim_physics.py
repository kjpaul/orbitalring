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
R_ORBIT = 6_628_137               # radius of 250 km orbit (m)
V_GROUND_STATIONARY = 483.331   # Ground-stationary velocity at 250 km (m/s)
L_RING = 41_645_813.012         # Ring circumference (m)
A_250_KM = 9.038                # net acceleration at 250 km geostationary orbit (m/s²)

# Deep space
T_SPACE = 2.7                   # Deep space temperature (K)

# Material properties
RHO_CNT = 1700                  # CNT density (kg/m³)


# =============================================================================
# CABLE MASS
# =============================================================================

def calc_cable_mass(load_mass=12000, sigma_target=12.633E9):
    """Calculate cable mass based on --m_load=NUMBER from argv. m_load must be > 999."""
    F_load_m = A_250_KM * load_mass 
    return F_load_m * R_ORBIT * RHO_CNT / sigma_target    

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


def calc_skin_depth(f_slip, temp_K, rho_293K, alpha, material, v_slip_min, tau_p=100):
    """Electromagnetic skin depth in the reaction plate."""
    rho = calc_resistivity(temp_K, rho_293K, alpha, material)
    if f_slip <= 0:
        f_slip = v_slip_min / (2 * tau_p)
    return math.sqrt(rho / (math.pi * MU0 * f_slip))


def calc_effective_plate_depth(f_slip, temp_K, t_plate, rho_293K, alpha, material, v_slip_min, tau_p=200):
    """Effective conducting depth: minimum of plate thickness and skin depth."""
    delta = calc_skin_depth(f_slip, temp_K, rho_293K, alpha, material, v_slip_min, tau_p)
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


def calc_thrust_power(thrust, v_rel):
    """Mechanical power delivered to the cable."""
    return thrust * v_rel

# =============================================================================
# ELECTRICAL CALCULATIONS
# =============================================================================

def calc_coil_voltage(i_peak, f_supply, b_field, a_coil, n_turns):
    """Induced voltage in LIM coil from changing flux linkage."""
    # this check is needed since i_peak is set after volts. Neither should be negative.
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    phi_peak = b_field * a_coil
    omega = 2 * math.pi * f_supply
    return omega * n_turns * phi_peak


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

def calc_cable_equilibrium_temperature(p_eddy_per_lim, lims_per_site, lim_sites, 
                                       cable_emissivity, cable_surface_area_per_m,
                                       q_env_per_m=0.0):
    """Calculate cable equilibrium temperature when eddy heat stays in cable.
    
    The cable must radiate away both:
    1. Eddy current heat from LIM operation (distributed around the ring)
    2. Environmental heat absorbed through MLI shielding (sun + earth)
    
    Args:
        p_eddy_per_lim: Eddy current losses per LIM (W)
        lims_per_site: Number of LIMs per site
        lim_sites: Total number of LIM sites around the ring
        cable_emissivity: Emissivity of cable surface
        cable_surface_area_per_m: Radiating surface area per meter of cable (m²/m)
        q_env_per_m: Environmental heat absorbed per meter (W/m)
    
    Returns:
        Equilibrium temperature (K)
    """
    # Eddy heat distributed around the ring
    p_eddy_total = p_eddy_per_lim * lims_per_site * lim_sites
    q_eddy_per_m = p_eddy_total / L_RING
    
    # Total heat input per meter = eddy + environmental
    q_in_per_m = q_eddy_per_m + q_env_per_m
    
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
# TWO-NODE WARM THERMAL MODEL (thermosyphon)
# =============================================================================

def calc_two_node_thermal(p_eddy_per_lim, lims_per_site, lim_sites,
                          cable_emissivity, cable_radiating_area_per_m,
                          wall_emissivity, warm_radiator_area_per_m,
                          warm_radiator_emissivity, q_env_per_m=0.0):
    """Two-node thermal model for cable and casing wall temperatures.
    
    The warm thermal management system uses a gravity-driven thermosyphon:
    a working fluid evaporates at warm surfaces inside the casing, vapor
    rises to external radiator panels where it condenses, and liquid 
    returns by gravity. This moves heat with negligible electrical power.
    
    Node 1 (T_cable): The cable interior.
        Heat in:  eddy current losses (ring-averaged) + environmental MLI leak
        Heat out: radiation to casing inner wall
        
    Node 2 (T_wall): Casing wall / thermosyphon evaporator.
        Heat in:  radiation from cable
        Heat out: thermosyphon carries heat to external radiators → space
        
    The thermosyphon is efficient enough that T_wall ≈ T_radiator.
    At equilibrium, all heat flows through both nodes to space.
    
    Note: inverter waste heat and cryo hot-side rejection are handled
    separately — inverters are local to each LIM site on the casing,
    and the cryo system has its own dedicated radiators.
    
    Args:
        p_eddy_per_lim: Eddy current losses per LIM (W)
        lims_per_site: Number of LIMs per site
        lim_sites: Total number of LIM sites around the ring
        cable_emissivity: Emissivity of cable surface
        cable_radiating_area_per_m: Cable surface facing casing wall (m²/m)
        wall_emissivity: Emissivity of casing inner wall
        warm_radiator_area_per_m: External radiator area per meter (m²/m)
        warm_radiator_emissivity: Emissivity of external warm radiators
        q_env_per_m: Environmental heat absorbed per meter through MLI (W/m)
    
    Returns:
        (T_cable, T_wall): Equilibrium temperatures (K)
    """
    # Total heat input per meter (ring-averaged)
    p_eddy_total = p_eddy_per_lim * lims_per_site * lim_sites
    q_eddy_per_m = p_eddy_total / L_RING
    q_total_per_m = q_eddy_per_m + q_env_per_m
    
    # At equilibrium, all heat exits through the external radiators.
    # External radiator: q_total = ε_rad * σ * A_rad * (T_wall⁴ - T_space⁴)
    # Solve for T_wall:
    rad_coeff_ext = warm_radiator_emissivity * STEFAN_BOLTZMANN * warm_radiator_area_per_m
    
    if rad_coeff_ext <= 0:
        return (500.0, 500.0)  # Fallback: no radiator
    
    T_wall_4 = q_total_per_m / rad_coeff_ext + T_SPACE**4
    T_wall = T_wall_4 ** 0.25
    
    # Cable-to-wall radiation link:
    # q_total = ε_eff * σ * A_cable * (T_cable⁴ - T_wall⁴)
    # where ε_eff = 1/(1/ε_cable + 1/ε_wall - 1) for parallel surfaces
    em_eff = 1.0 / (1.0/cable_emissivity + 1.0/wall_emissivity - 1.0)
    rad_coeff_cw = em_eff * STEFAN_BOLTZMANN * cable_radiating_area_per_m
    
    if rad_coeff_cw <= 0:
        return (500.0, T_wall)  # Fallback
    
    T_cable_4 = q_total_per_m / rad_coeff_cw + T_wall**4
    T_cable = T_cable_4 ** 0.25
    
    return (T_cable, T_wall)


def calc_warm_radiator_width(q_total_per_m, warm_radiator_emissivity, 
                              warm_radiator_length, T_wall):
    """Calculate required warm radiator width for given wall temperature.
    
    This is the width of external radiator panels needed to reject the
    warm-loop heat at the thermosyphon operating temperature.
    
    Args:
        q_total_per_m: Total warm-loop heat per meter of ring (W/m)
        warm_radiator_emissivity: Emissivity of external radiators
        warm_radiator_length: Length of radiator per LIM spacing (m)
        T_wall: Casing wall / radiator temperature (K)
    
    Returns:
        Required radiator width (m)
    """
    if q_total_per_m <= 0:
        return 0.0
    
    q_per_m2 = warm_radiator_emissivity * STEFAN_BOLTZMANN * (T_wall**4 - T_SPACE**4)
    
    if q_per_m2 <= 0:
        return 1000.0  # Fallback
    
    # Total heat per site spacing
    q_per_site_length = q_total_per_m * warm_radiator_length
    area_required = q_per_site_length / q_per_m2
    return area_required / warm_radiator_length


# =============================================================================
# LEVITATION COIL HEAT LOAD
# =============================================================================

def calc_levitation_coil_heat(T_cable, lev_coil_area, q_ref, T_ref, T_cryo=77.4):
    """Calculate heat load on levitation coils from warm cable above.
    
    The DC levitation coils are wrapped in MLI but still absorb heat from
    the warm cable through radiative transfer. The heat flux scales with
    the T^4 difference.
    
    Args:
        T_cable: Cable temperature (K)
        lev_coil_area: Levitation coil area per site (m²)
        q_ref: Reference MLI heat flux at T_ref (W/m²)
        T_ref: Reference hot-side temperature (K)
        T_cryo: Cryogenic temperature (K), default 77.4 K (LN2)
    
    Returns:
        Heat load on levitation coils (W)
    """
    # MLI heat flux scales approximately with T_hot^4 - T_cold^4
    ref_term = T_ref**4 - T_cryo**4
    actual_term = T_cable**4 - T_cryo**4
    
    if ref_term <= 0:
        return 0.0
    
    q_actual = q_ref * (actual_term / ref_term)
    return q_actual * lev_coil_area


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
                        eddy_to_cable, q_coil_environment, q_lev_coil=0.0):
    """Calculate heat load that must go through cryogenic system.
    
    Args:
        p_eddy: Eddy current losses per LIM (W)
        p_hyst: Hysteresis losses per LIM (W)
        lims_per_site: Number of LIMs per site
        q_absorbed_per_site: Environmental heat absorbed per site (W)
        eddy_to_cable: True if eddy heat goes to cable, False if to cryo
        q_coil_environment: Heat leak through LIM coil MLI (W)
        q_lev_coil: Heat load on levitation coils from warm cable (W)
    
    Returns:
        Total heat load for cryogenic system (W)
    """
    if eddy_to_cable:
        # Eddy heat stays in cable, cryo handles:
        # - HTS hysteresis losses
        # - Heat leak through LIM coil MLI
        # - Heat from warm cable to levitation coils
        q_cold = lims_per_site * p_hyst + q_coil_environment + q_lev_coil
    else:
        # Everything goes through cryo (legacy mode)
        total_heat = lims_per_site * (p_eddy + p_hyst) + q_absorbed_per_site + q_lev_coil
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
