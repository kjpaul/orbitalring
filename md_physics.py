"""
Mass Driver Physics Module
Calculates thrust, drag, and thermal loads for the launch sled.
"""
import math
import md_config as cfg

# Physical Constants
MU0 = 4 * math.pi * 1e-7
STEFAN_BOLTZMANN = 5.670374e-8

def get_resistivity(material, temp_K):
    """
    Returns resistivity in Ohm-m.
    Gamma-TiAl is significantly more resistive than Aluminum.
    """
    if material == "gamma_titanium":
        # Approx 30 micro-ohm-cm at room temp -> 3.0e-7 Ohm-m
        # Rises with temperature
        base_rho = 3.0e-7
        temp_coeff = 0.001 # Approximation per degree K
        return base_rho * (1 + temp_coeff * (temp_K - 293))
    elif material == "aluminum":
        return 2.65e-8 * (1 + 0.004 * (temp_K - 293))
    return 1e-6 # Default fallback

def calc_goodness_factor(freq, tau_p, gap, rho_plate, plate_thickness):
    """
    Laithwaite's Goodness Factor (G).
    G < 1 means the motor is effectively just a heater.
    We need G >> 1 for efficient thrust.
    """
    omega = 2 * math.pi * freq
    # Equivalent conductivity calculation for the sheet rotor
    sigma_surf = plate_thickness / rho_plate
    
    # G = (2 * mu0 * f * tau^2) / (pi * rho/d * g_eff)
    # Simplified proportional metric:
    g_factor = (2 * MU0 * freq * tau_p**2 * sigma_surf) / (math.pi * gap)
    return g_factor

def calc_thrust_force(v_sled, v_slip, temp_K):
    """
    Calculates total thrust on the 1000m sled.
    """
    # 1. Electromagnetic Setup
    v_sync = v_sled + v_slip
    freq = v_sync / (2 * cfg.TAU_P)
    omega = 2 * math.pi * freq
    k = math.pi / cfg.TAU_P
    
    rho = get_resistivity(cfg.PLATE_MATERIAL, temp_K)
    
    # 2. Magnetic Field (Stator Surface)
    # B_peak approx mu0 * N * I / Gap (simplified for air-cored/iron-less)
    # We apply a factor for winding distribution
    b_peak = (MU0 * cfg.N_TURNS * cfg.I_PEAK) / (cfg.GAP * 1.1) 
    
    # 3. Goodness Factor
    G = calc_goodness_factor(freq, cfg.TAU_P, cfg.GAP, rho, cfg.PLATE_THICKNESS)
    
    # 4. Thrust Calculation (Rotor Efficiency Model)
    # F_x = (Length * Width) * (B^2 / 2*mu0) * (Factor depending on G and Slip)
    # Ideally: Thrust peaks when s * G ~ 1 (for some machine types) or increases with s
    
    # Calculate Slip S
    s = v_slip / v_sync
    
    # Classical induction motor thrust curve approximation per m^2
    # F_density = 0.5 * B^2/mu0 * (sG / (1 + (sG)^2)) * GeometricFactors
    # Note: For linear motors with large gaps, we use a simplified Lorentz model
    
    magnetic_pressure = (b_peak**2) / (2 * MU0)
    
    # The "Induction Factor" determines how much B penetrates and pushes
    # For high resistivity plates, we need high slip to get current.
    induction_efficiency = (s * G) / (1 + (s * G)**2)
    
    force_density = magnetic_pressure * induction_efficiency * 2 # *2 for two sides (top/bottom or left/right)
    
    active_area = cfg.SLED_LENGTH * cfg.PLATE_WIDTH
    total_thrust = force_density * active_area
    
    # 5. Drag (Atmospheric is negligible at 250km, but magnetic drag exists)
    # We assume drag is handled by the efficiency loss in net force
    
    return total_thrust, G, b_peak

def calc_heating(thrust, v_slip):
    """
    Heat generated in the sled plates = Thrust * Slip Velocity
    P_heat = F * v_slip
    """
    return thrust * v_slip