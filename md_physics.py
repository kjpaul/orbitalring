import math
import numpy as np
import md_config as cfg

MU0 = 4 * math.pi * 1e-7

def calc_resistivity(material_props, temp_K):
    rho_0 = material_props['rho_293']
    alpha = material_props['alpha']
    return rho_0 * (1 + alpha * (temp_K - 293.0))

def calc_norris_loss(i_peak, tape_cfg):
    """Calculates hysteresis loss per meter of tape."""
    ic = tape_cfg.max_safe_current / 0.8 # Back-calculate true Ic from safe limit
    if ic <= 0: return 0.0
    
    i_ratio = i_peak / ic
    i_ratio = max(-0.999, min(0.999, i_ratio))
    
    # Loss J/cycle/m
    loss = (MU0 * ic**2 / math.pi) * ((1 - i_ratio) * math.log(1 - i_ratio) + 
                                      (1 + i_ratio) * math.log(1 + i_ratio) - 
                                      i_ratio**2)
    return loss * tape_cfg.layers # Loss scales with layers

def solve_stage_performance(stage, v_sled, v_slip, temp_K, sled_cfg):
    """
    Calculates Thrust, Voltage, Power for ONE instance of a LIM stage.
    """
    mat = sled_cfg.mat_props
    
    # 1. Frequencies
    v_sync = v_sled + v_slip
    f_supply = v_sync / (2 * stage.tau_p)
    f_slip = v_slip / (2 * stage.tau_p)
    omega_supply = 2 * math.pi * f_supply
    omega_slip = 2 * math.pi * f_slip
    
    # Check Handoff / Cutoff
    if f_supply > stage.handoff_freq:
        return 0, 0, 0, 0, 0, 0  # Stage disabled (too fast)
    
    # 2. Impedance (Model 1)
    rho = calc_resistivity(mat, temp_K)
    
    # Skin depth
    if f_slip > 0:
        delta = math.sqrt(rho / (math.pi * MU0 * f_slip))
    else:
        delta = 999.0
        
    d_eff = min(sled_cfg.thickness, delta)
    
    # Loop Resistance & Inductance
    # Resistance of the path in the plate
    R_plate = rho * stage.tau_p / (d_eff * sled_cfg.height)
    
    # Loop Inductance (Rectangular loop approx)
    # L = mu0/pi * tau * ln(2*tau/w)
    # Clamp log term to avoid negatives if geometry is weird
    geom_factor = max(1.1, 2 * stage.tau_p / sled_cfg.height)
    L_plate = (MU0 * stage.tau_p / math.pi) * math.log(geom_factor)
    
    X_plate = omega_slip * L_plate
    Z_plate = math.sqrt(R_plate**2 + X_plate**2)
    
    # 3. Optimize Current (Binary Search)
    # Maximize current such that V < 100kV and I < I_limit
    
    i_limit = stage.tape.max_safe_current
    i_min, i_max = 0.0, i_limit
    
    best_I = 0.0
    best_F = 0.0
    best_V = 0.0
    
    # Stator Inductance (for voltage calc)
    # L_stator approx N^2 * mu0 * Area / Gap
    area_coil = stage.tau_p * stage.w_coil
    L_stator = (MU0 * stage.turns**2 * area_coil) / stage.gap
    
    for _ in range(10):
        i_try = (i_min + i_max) / 2
        
        # B-field
        B_gap = (MU0 * stage.turns * i_try) / stage.gap
        
        # Induced EMF in Plate
        # EMF = v_slip * B * Height * (coupling factor 2/pi)
        EMF = (2/math.pi) * v_slip * B_gap * sled_cfg.height
        
        # Eddy Current
        I_eddy = EMF / Z_plate
        
        # Thrust (Power / v_slip)
        # Active loops per stage = Length_Active / Tau
        n_loops = stage.length_active / stage.tau_p
        P_mech_stage = n_loops * I_eddy**2 * R_plate
        F_stage = P_mech_stage / v_slip if v_slip > 0 else 0
        
        # Voltage Check (Inductive drop)
        V_ind = omega_supply * L_stator * i_try
        
        if V_ind > cfg.VOLTS_MAX:
            i_max = i_try
        else:
            i_min = i_try
            best_I = i_try
            best_F = F_stage
            best_V = V_ind

    # 4. Power Calcs
    # Cryo / Hysteresis
    loss_per_m = calc_norris_loss(best_I, stage.tape)
    # Total tape length in this stage
    # (Turns * 2 * (Tau + W_coil)) * Phases(3)
    len_tape = stage.turns * 2 * (stage.tau_p + stage.w_coil) * 3
    P_hyst = loss_per_m * len_tape * f_supply
    P_cryo_wall = P_hyst * cfg.CRYO_PENALTY_RATIO
    
    P_mech_total = best_F * v_sled
    P_total_wall = P_mech_total + P_cryo_wall
    
    return best_F, best_I, best_V, P_total_wall, P_mech_total, P_cryo_wall