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
    """
    Norris elliptical model for AC loss in superconducting strips.
    """
    if i_c <= 0: return 0.0
    ii = i_peak / i_c
    # Avoid log domain errors
    ii = max(-0.999, min(0.999, ii))
    
    # Loss per cycle per meter length
    loss = (MU0 * i_c**2 / math.pi) * ((1 - ii) * math.log(1 - ii) + (1 + ii) * math.log(1 + ii) - ii**2)
    return loss

def calc_thrust_model1(v_sled, v_slip, i_peak, temp_K):
    """
    Model 1: Narrow Plate Eddy Current Model.
    Focuses on loop inductance and resistance.
    """
    # 1. Frequency setup
    v_sync = v_sled + v_slip
    f_supply = v_sync / (2 * cfg.TAU_P)
    f_slip = v_slip / (2 * cfg.TAU_P)
    omega_slip = 2 * math.pi * f_slip
    
    if f_supply <= 0.001: return 0.0, 0.0, 0.0, 0.0

    # 2. Field Calculation
    # B at plate (approximate with spacing decay)
    b_coil = (MU0 * cfg.N_TURNS * i_peak) / cfg.W_COIL
    # Decay over gap (simple exponential approx for validation)
    b_plate = b_coil * math.exp(-math.pi * cfg.GAP / cfg.TAU_P)

    # 3. Model 1 Parameters
    rho = get_resistivity(temp_K)
    
    # Skin depth check (though Model 1 assumes current fills the "effective" depth)
    delta = math.sqrt(rho / (math.pi * MU0 * f_slip)) if f_slip > 0 else 999.0
    d_eff = min(cfg.T_PLATE, delta)
    
    # R (Resistance of eddy loop) & L (Inductance)
    # R = rho * Length / Area
    R_loop = rho * cfg.TAU_P / (d_eff * cfg.W_COIL)
    
    # L = Inductance of rectangular loop
    L_loop = (MU0 * cfg.TAU_P / math.pi) * math.log(2 * cfg.TAU_P / cfg.W_COIL)
    
    X_loop = omega_slip * L_loop
    Z_loop = math.sqrt(R_loop**2 + X_loop**2)
    
    # 4. Induced Current & Thrust
    # EMF induced in the plate loop
    # EMF = v_slip * B * Width (approx)
    EMF = 4.0 * v_slip * b_plate * cfg.W_COIL / math.pi
    
    I_eddy = EMF / Z_loop
    
    # Power dissipated in Plate (Thrust Power)
    # Number of loops active under the sled
    n_loops = cfg.SLED_LENGTH / cfg.TAU_P
    
    P_thrust_mech = n_loops * I_eddy**2 * R_loop
    
    if v_slip > 0:
        F_thrust = P_thrust_mech / v_slip
    else:
        F_thrust = 0.0
        
    return F_thrust, f_supply, I_eddy, b_plate

def calc_power_balance(v_sled, v_slip, i_peak, thrust, f_supply):
    """
    Calculates total site power including the Massive Cryo Penalty.
    """
    # 1. Mechanical Power (delivered to sled)
    p_mech = thrust * v_sled
    
    # 2. Hysteresis Loss (The Cryo Killer)
    # Loss per meter of tape * total tape length * frequency
    loss_J_m = q_hysteresis_norris(i_peak, cfg.I_C_TOTAL)
    
    # Total HTS length involved per site (active coils)
    # Approx: Turns * Circumference of coil * Phases
    len_coil = 2 * (cfg.TAU_P + cfg.W_COIL) * cfg.N_TURNS
    # Only active coils under sled? Or whole sector? 
    # Standard: Losses occur in active coils. Sled covers SLED_LENGTH.
    # Coils active = Sled_Length / Tau_P * Phases
    n_active_coils = (cfg.SLED_LENGTH / cfg.TAU_P) * 3 # 3-phase
    
    p_hysteresis = loss_J_m * len_coil * n_active_coils * f_supply
    
    # 3. Cryogenic Penalty
    # The heat (p_hysteresis) is at 70K. It must be pumped to 300K.
    # COP penalty applies to the electrical wall power needed.
    p_cryo_wall = p_hysteresis * cfg.CRYO_COP_PENALTY
    
    # 4. Total Electrical Load
    p_total = p_mech + p_cryo_wall # Neglecting resistive copper losses for now
    
    return p_total, p_cryo_wall, p_hysteresis

def calc_voltage(i_peak, f_supply):
    """
    Inductive voltage drop check.
    V ~ L * di/dt ~ omega * L * I
    """
    omega = 2 * math.pi * f_supply
    # Inductance of the stator coil itself
    # Approx solenoid: mu0 * N^2 * A / l
    A_coil = cfg.TAU_P * cfg.W_COIL
    L_stator = (MU0 * cfg.N_TURNS**2 * A_coil) / cfg.GAP # Rough air-core approx
    
    V_peak = omega * L_stator * i_peak
    return V_peak