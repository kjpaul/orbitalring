"""
Mass Driver Physics Module - Model 1 & Cryo Constraints
"""
import math
import md_config as cfg

MU0 = 4 * math.pi * 1e-7

def get_resistivity(temp_K):
    # Gamma-TiAl resistivity model
    rho_293 = 7.5e-7 
    alpha = 0.0012
    return rho_293 * (1 + alpha * (temp_K - 293))

def q_hysteresis_norris(i_peak, i_c):
    if i_c <= 0: return 0.0
    ii = i_peak / i_c
    ii = max(-0.999, min(0.999, ii))
    loss = (MU0 * i_c**2 / math.pi) * ((1 - ii) * math.log(1 - ii) + (1 + ii) * math.log(1 + ii) - ii**2)
    return loss

def calc_thrust_model1(v_sled, v_slip, i_peak, temp_K):
    v_sync = v_sled + v_slip
    f_supply = v_sync / (2 * cfg.TAU_P)
    f_slip = v_slip / (2 * cfg.TAU_P)
    
    if f_supply <= 0.001: return 0.0, 0.0, 0.0, 0.0

    b_coil = (MU0 * cfg.N_TURNS * i_peak) / cfg.W_COIL
    b_plate = b_coil * math.exp(-math.pi * cfg.GAP / cfg.TAU_P)

    rho = get_resistivity(temp_K)
    
    # Model 1 Impedance
    delta = math.sqrt(rho / (math.pi * MU0 * f_slip)) if f_slip > 0 else 999.0
    d_eff = min(cfg.T_PLATE, delta)
    
    R_loop = rho * cfg.TAU_P / (d_eff * cfg.W_COIL)
    L_loop = (MU0 * cfg.TAU_P / math.pi) * math.log(2 * cfg.TAU_P / cfg.W_COIL)
    X_loop = 2 * math.pi * f_slip * L_loop
    Z_loop = math.sqrt(R_loop**2 + X_loop**2)
    
    # Thrust
    EMF = 4.0 * v_slip * b_plate * cfg.W_COIL / math.pi
    I_eddy = EMF / Z_loop
    n_loops = cfg.SLED_LENGTH / cfg.TAU_P
    P_mech = n_loops * I_eddy**2 * R_loop
    
    F_thrust = P_mech / v_slip if v_slip > 0 else 0.0
    return F_thrust, f_supply, I_eddy, b_plate

def calc_power_balance(v_sled, v_slip, i_peak, thrust, f_supply):
    p_mech = thrust * v_sled
    loss_J_m = q_hysteresis_norris(i_peak, cfg.I_C_TOTAL)
    len_coil = 2 * (cfg.TAU_P + cfg.W_COIL) * cfg.N_TURNS
    n_active_coils = (cfg.SLED_LENGTH / cfg.TAU_P) * 3 
    p_hysteresis = loss_J_m * len_coil * n_active_coils * f_supply
    p_cryo_wall = p_hysteresis * cfg.CRYO_COP_PENALTY
    return p_mech + p_cryo_wall, p_cryo_wall, p_hysteresis

def calc_voltage(i_peak, f_supply):
    omega = 2 * math.pi * f_supply
    A_coil = cfg.TAU_P * cfg.W_COIL
    L_stator = (MU0 * cfg.N_TURNS**2 * A_coil) / cfg.GAP 
    return omega * L_stator * i_peak