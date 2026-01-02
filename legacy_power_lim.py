# flake8: noqa
"""
    Calculations for "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import math
import tabulate as tabulate
import matplotlib.pyplot as plt


# -------------------------------------------------
# VARIABLE PARAMETERS, TO BE SET EXTERNALLY
# -------------------------------------------------

# -------------------------------------------------
# HTS TAPE AND FREQUENTLY CHANGED PARAMETERS
# -------------------------------------------------
# HTS tape comes in standard widths: 12mm, 6mm, 4mm, 3mm
# Critical current (Ic) scales linearly with tape width
# Multiple LIMs can be placed on each side of the cable
N_TURNS = 40                               # turns per phase coil (int)
LIMS_PER_SIDE = 1                           # How many LIMs on each side of cable (1, 2, 3, ...)
HTS_TAPE_WIDTH_MM = 12                       # tape width in mm (12, 6, 4, or 3)
HTS_TAPE_LAYERS = 2                         # number of tape layers (1 or 2)
V_SLIP_MAX = 200.0                          # maximum slip velocity (m/s)
V_SLIP_MIN = 5.0                            # minimum slip velocity (m/s)
SLIP_RATIO_NORMAL = 0.02                    # target slip ratio at full current (2%)
SLIP_RATIO_REDUCED = 0.005                  # target slip ratio when current-limited (1%)
TAU_P = 100.0                               # pole-pitch (m)
W_COIL = 0.5                                # LIM width (m)
GAP = 0.05                                  # coil-to-plate gap (m)
LIM_SPACING = 500.0                         # distance at which LIMs are place (m)
T_PLATE = 0.05                              # aluminium thickness (m)
IC_PER_MM_PER_LAYER = 66.7                  # Ic per mm of tape width per layer, 66 max (A/mm)
MAX_SITE_POWER = 16.0e6                     # power limit per LIM site (W)

# graph and run-time variables 
WRITE_FILE = True
MAKE_GRAPHS = True

# -------------------------------------------------
# TIME CONTROL
# -------------------------------------------------
SKIP = 200   # if 1 = no skip
DT1 = 1   # loop dt for first day
DT2 = 10    # loop dt until sample time max, which also ends data collection for graphs
DT3 = 50    # loop dt until end

HR = 60 * 60
DAY = HR * 24
WEEK = 7 * DAY
YR = round(DAY * 365.33)
MONTH = 30 * DAY
SAMPLE_TIME_MAX = 5 * YR      # length of sample (must be longer than deployment)
                                               
# -------------------------------------------------
# CONTROL / COUPLING PARAMETERS (tunable)
# -------------------------------------------------
# Thrust efficiency factor: multiplier on the idealized model F = s * F_max
# The idealized model (THRUST_EFFICIENCY = 1.0) assumes resistive regime and is
# likely conservative for this high magnetic Reynolds number LIM (R_m >> 1).
# Values of 5-20 may be realistic based on shielding effects, but require
# detailed analysis or FEM simulation to validate for this specific geometry.
THRUST_EFFICIENCY = 1.0     # 1.0 = conservative idealized model

# -------------------------------------------------
# THRUST MODEL SELECTION
# -------------------------------------------------
# Select which thrust model to use:
#   1 = Narrow plate eddy current model (geometry-corrected for W << tau_p)
#   2 = Goodness factor model (Laithwaite, assumes W > tau_p)
#   3 = Slip x pressure model (F = slip * F_max, theoretical maximum)
# 
# Use command line: --model=1, --model=2, or --model=3
THRUST_MODEL = 1  # Default to narrow plate model

# Slip control notes:
# - V_SLIP_MIN ensures minimum thrust at low v_rel (startup regime)
# - V_SLIP_MAX limits losses at high v_rel
# - SLIP_RATIO_NORMAL (2%) used when current is at target
# - SLIP_RATIO_REDUCED (1%) used when current is power-limited
# - The slip target tapers between these based on current headroom

SLIP_MIN = 0.005
CURRENT_UPRATE = 1.01       # multiplicative current ramp per controller iteration when under limits
POWER_HEADROOM = 0.98       # aim to use this fraction of MAX_SITE_POWER
MAX_HEATSINK_AREA = (LIM_SPACING) * 2 * W_COIL # the heatsink extends under the coils, since they are 99.99 % open space. 
A_COIL = TAU_P * W_COIL     # coil area

# Derived HTS parameters
W_TAPE = HTS_TAPE_WIDTH_MM / 1000               # HTS tape width (m)
I_C = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * HTS_TAPE_LAYERS  # Critical current (A)
I_PEAK = 0.875 * I_C                            # I_peak, typically ~87.5% of Ic (A)
I_TARGET = 0.8125 * I_C                         # I_target, typically ~81.25% of Ic (A)

# -------------------------------------------------
# OTHER LIM PARAMETERS
# -------------------------------------------------
# I_C, I_PEAK, I_TARGET, W_TAPE are now computed from HTS_TAPE_WIDTH_MM and HTS_TAPE_LAYERS above
I_MIN = 10.0                                # lower limit on current. Reduce slip instead. (A)
P_HEAT_MAX = 100000
HTS_D = 80.0                                # HTS thickness in micrometers (kept for reference)
VOLTS_MAX = 100e3                           # absolute peak coil voltage limit, per your assumption (V)
ALPHA_ANGLE_DEG = 20.0                      # magnetic penetration angle of HTS in coils (deg)
HEAT_SINK_L = LIM_SPACING                   # shorter heat sink mean higher ave reaction plate temp.
LIMS_PER_SITE = 2 * LIMS_PER_SIDE               # Total LIMs per site

# Variables that could change, but probably won't
PITCH_COUNT = 3                             # pitches per LIM (int)
CASING_OUTER_W = 10.0                       # external width of orbital ring casing (m)
LIM_PHASES = 3                              # power supply phases into LIMs (int)
INV_EFF = 0.90                              # DC to AC inverter efficiency factor (int)
LIM_EFF = 0.95                              # unaccounted for LIM losses efficiency factor (int)

# -------------------------------------------------
# CABLE and LOAD masses are dependant on each other 
# according to Newtown's 3rd Law and the cable tension.
# -------------------------------------------------
M_CABLE_M = 99_198                         # orbital ring cable mass per meter (kg/m)
M_LOAD_M = 12_000                            # orbital ring casing + load mass per meter (kg/m)

# lable data for plots
SAMPLE_TIME_STR = "-- 100 tonne cable, 12 tonne load" # graph title
PARAM_STR1 = f"τp={TAU_P} m, N={N_TURNS}, V_Slip_max={V_SLIP_MAX} m/s, d={LIM_SPACING} m"

# -------------------------------------------------
# BASIC CONSTANTS & USER TUNABLE DESIGN PARAMETERS
# -------------------------------------------------
V_ORBIT = 7754.866                          # 250 km orbital velocity (m/s)
GEO_250_ORBIT = 483.331                     # 250 km ground stationary velocity (m/s)
MU0 = 4 * math.pi * 1e-7                    # permeability of a vacuum (kg m s^-2 A^-2)
SIGMA = 5.670374e-8                         # Stephan-Boltzmann constant (kg s^-3 K^-4)
RHO_ALU_E_70K = 4.853e-9                    # Al temp starts at 70 K(Ωm)
RHO_ALU_M3 = 2700                           # mass per m^3 of aluminium (kg/m^3)
L_ACTIVE = TAU_P * PITCH_COUNT              # length of LIM (m)
A_LIM = L_ACTIVE * W_COIL                   # area under LIM (m^2)
L_HTS_COIL = 2 * (W_COIL + TAU_P) * N_TURNS # HTS length in one coil (m)
L_HTS_LIM = L_HTS_COIL * LIM_PHASES * PITCH_COUNT # HTS in one LIM (m)
L_RING_250 = 41645813.012                   # length of ring in meters (m)
LIM_SITES = round(L_RING_250 / LIM_SPACING) # total number of LIMs (int)
ALPHA_TAPE = ALPHA_ANGLE_DEG * math.pi / 180 # field penetration angle of HTS (rad)
M_CABLE_T = M_CABLE_M * L_RING_250          # orbital ring cable mass total (m)
M_LOAD_T = M_LOAD_M * L_RING_250            # orbital ring casing + load mass total (m)

# -------------------------------------------------
# HEAT RELATED CONSTANTS
# -------------------------------------------------
Q_SUN_1AU = 1361        # solar energy at 1 AU (W/m^2)
Q_EARTH_DAY = 650       # average daytime energy reflected up from earth's surface (W/m^2)
Q_SHIELDING = 0.005     # percentage of heat radiation that passes through multi-layer alu shielding (%/100)
Q_ABS_M = (Q_SUN_1AU + Q_EARTH_DAY) * CASING_OUTER_W * Q_SHIELDING # external heat absorbed by ring (W/m)
Q_ABS_LIM = Q_ABS_M * LIM_SPACING # external heat absorbed by ring (W)
T_FREE_SPACE = 2.7      # temperature of free space (K)
CRYO_EFF = 0.18         # Cryo efficiency: 0.031 -> 32 W are needed for Cryo per W of heat (1/W)
C_P_ALU = 900           # heat capacity or aluminium (J/kg K)
H_CONV_LN2 = 90         # Convective Heat Transfer Coefficient LN2 (W m^-2 K^-1)
T2_LN2_BOIL = 77.4      # boiling point of liquid nitrogen at 1 atm (K)
T1_LN2_CRYO = 70        # temperature of LN2 as it leaves cryogenic system (K)
T_N2_HOT = 300          # N2 temperature as it leave the compressor and enters the radiator (K)
T_MAX_REACTION_PLATE = 500  # maximum temp for aluminium reaction plate (K)
K_ALU = 205             # thermal conductivity of aluminium (W/m*K)
C_P_LN2 = 2040          # Specific heat capacity of LN2 (J/kg*K)
L_V_LN2 = 199000        # Latent heat of vaporization LN2 (J/kg)
EM_ALU = 0.85           # emissivity  or reaction plate (black anodized alu) (num)
EM_HEAT_SINK = 0.9      # emissivity  of heat sink surface (num)
PA_LN2_70K = 0.4        # pressure that turns boiling point of LN2 to 70 K (atm)
V_REL_MIN = 10          # fudge factor compensate of misleading initial heat at low v_rel (m/s)


PARAM_LIST = {
    # HTS Configuration
    "N_TURNS             #": N_TURNS,
    "TAU_P               m": TAU_P,
    "W_COIL              m": W_COIL,
    "SLIP_RATIO_NORMAL   %": SLIP_RATIO_NORMAL * 100,
    "SLIP_RATIO_REDUCED  %": SLIP_RATIO_REDUCED * 100,
    "GAP                mm": GAP * 1000,
    "HTS_TAPE_WIDTH_MM  mm": HTS_TAPE_WIDTH_MM,
    "HTS_TAPE_LAYERS     #": HTS_TAPE_LAYERS,
    "W_TAPE              m": W_TAPE,
    "I_C                 A": round(I_C),
    "I_PEAK              A": round(I_PEAK),
    "I_TARGET            A": round(I_TARGET),
    "V_SLIP_MAX        m/s": V_SLIP_MAX,
    "V_SLIP_MIN        m/s": V_SLIP_MIN,
    "SLIP_MIN            %": SLIP_MIN * 100,
    # LIM Configuration
    "LIMS_PER_SIDE       #": LIMS_PER_SIDE,
    "LIMS_PER_SITE       #": LIMS_PER_SITE,
    "PITCH_COUNT         #": PITCH_COUNT,
    "LIM_SPACING         m": LIM_SPACING,
    # Operating limits
    "THRUST_EFFICIENCY   %": THRUST_EFFICIENCY * 100,
    "I_MIN               A": I_MIN,
    "VOLTS_MAX          kV": VOLTS_MAX/1000,
    "P_HEAT_MAX         kW": P_HEAT_MAX/1000,
    "MAX_SITE_POWER     MW": MAX_SITE_POWER/1e6,
    # Thermal
    "T_PLATE            mm": T_PLATE * 1000,
    "ALPHA_ANGLE_DEG   deg": ALPHA_ANGLE_DEG,
    "CRYO_EFF            %": CRYO_EFF,
    "EM_ALU              %": EM_ALU * 100,
    "EM_HEAT_SINK        %": EM_HEAT_SINK * 100,
    "HEAT_SINK_L         m": HEAT_SINK_L,
    "MAX_HEATSINK_AREA  m²": MAX_HEATSINK_AREA,
    "V_REL_MIN         m/s": V_REL_MIN,
    "Q_ABS_LIM           W": Q_ABS_LIM,
    # Mass
    "CASING_OUTER_WIDTH  m": CASING_OUTER_W,
    "M_CABLE_M       tonne": M_CABLE_M/1000,
    "M_LOAD_M        tonne": M_LOAD_M/1000,
}
for key, value in PARAM_LIST.items():
    str1 = (f"\t{key:20}\t{value:8}")
    print(str1)
print("--------------------------------------------------------------------")


# chart lists data collectors
list_i_peak = []
list_volts = []
list_v_slip = []
list_slip = []
list_f_slip = []
list_thrust = []
list_thrust_power = []
list_v_rel = []
list_b_peak = []
list_p_eddy = []
list_skin_depth_eff = []
list_skin_depth_calc = []
list_p_cryo = []
list_p_hyst = []
list_p_lim = []
list_p_lim_site = []
list_temp_plate_ave = []
list_E_site_ke = []      # Track site KE over time (TJ)
list_E_total_ke = []     # Track total KE over time (EJ)
list_goodness_G = []     # Track goodness factor
list_eta_slip = []       # Track slip efficiency

"""
    -------------------------------------------------
    CALCULATIONS
    -------------------------------------------------
"""
"""
    orbital dynamics and ring velocities
"""
def get_v_cable(v0, thr, dt):
    thrust_total = thr * (LIM_SITES * 2)
    return v0 + (thrust_total / M_CABLE_T) * dt


def get_v_casing(v0, thr, dt):
    thrust_total = thr * (LIM_SITES * 2)
    return v0 - (thrust_total / M_LOAD_T) * dt


def get_v_rel(cable_velocity, casing_velocity):
    vrel = cable_velocity - casing_velocity
    if vrel < 0:
        sys.stdout.write("ALERT! VREL = ")
        sys.stdout.write(str(vrel))
        sys.stdout.flush()
    return vrel


"""
    velocity, slip and frequency equations
"""
def get_v_slip(slip, v_wave):
    v_slip = slip * v_wave
    if slip < 0 or slip > 1:
        sys.stdout.write("ALERT! slip = ")
        sys.stdout.write(str(slip))
        sys.stdout.flush()
    return v_slip


def get_slip(v_slip, v_rel):
    v_wave = v_slip + v_rel
    if v_wave < 1:
        v_wave = 1.0
    return (v_slip/v_wave)


def get_slip_f(f_slip, f_supply):
    return f_slip / f_supply


def get_slip_frequency(v_slip):
    return v_slip / (2 * TAU_P)


def get_supply_frequency(v_wave):
    # Electrical supply frequency (Hz) for a given travelling-wave speed.
    return v_wave / (2 * TAU_P)


"""
    magnetic fields
"""
def get_b_plate_peak(i_peak):
    """Traveling wave amplitude at reaction plate for 3-phase LIM.
    
    For a single-phase coil modeled as a rectangular current sheet:
        B_single = (2μ₀NI / πw) * arctan(w / 2g)
    
    For a balanced 3-phase system, three sinusoidal fields displaced by 120°
    in both time and space combine to form a constant-amplitude traveling wave.
    The backward-traveling components cancel; the forward components add:
        B_traveling = (3/2) * B_single * cos(ωt - kx)
    
    Since the traveling wave has constant amplitude (no pulsation), we use
    B_traveling directly in force calculations with no RMS conversion needed.
    """
    b_single_phase = (2 * MU0 * N_TURNS * i_peak / (math.pi * W_COIL)) * math.atan(W_COIL / (2 * GAP))
    return 1.5 * b_single_phase  # 3-phase traveling wave amplitude


def get_b_coil_peak(i_peak):
    """Peak field inside individual coil (for voltage/hysteresis calculations).
    
    This remains the single-phase value since each coil experiences its own
    field for purposes of flux linkage and hysteresis loss.
    """
    return MU0 * N_TURNS * i_peak / W_COIL


"""
    eddy currents & skin depth
"""
def get_plate_eddy_loss_from_thrust(v_slip, thrust):
    """Secondary (reaction plate) ohmic loss for one LIM, consistent with thrust.

    In a linear induction motor, the power dissipated in the secondary corresponds to the slip power:
        P_secondary = F * v_slip

    This ties losses to the force model and prevents the classic runaway that occurs when the
    secondary back-reaction is not solved self-consistently.

    Returns: Watts per LIM
    """
    if v_slip <= 0 or thrust <= 0:
        return 0.0
    return thrust * v_slip


def get_plate_eddy_loss(v_slip, i_peak, tempK):
    """LEGACY: retained for reference only (not used by default).

    This thin-plate / sinusoidal-field scaling can greatly over-predict at low f_slip and cryogenic rho
    unless it is solved self-consistently with field loading. Use get_plate_eddy_loss_from_thrust().
    """
    b_plate = get_b_plate_peak(i_peak)
    f_slip = get_slip_frequency(v_slip)
    eff_t_plate = get_eff_plate_depth(f_slip, tempK)
    rho_tempK = get_rho_alu(tempK)
    ce = math.pi**2 / (6 * rho_tempK)
    v_plate_eddy = get_plate_eddy_volume(f_slip, tempK)
    # return ((math.pi**2 * b_plate**2 * t_plate**2 * f_slip**2) / (6 * RHO_ALU_E_70K)) * v_plate_eddy
    return ce * b_plate**2 * eff_t_plate**2 * f_slip**2 * v_plate_eddy


def get_skin_depth_eddy(f_slip, tempK):
    # Skin depth penetration of eddy currents in reaction plate
    rho_tempK = get_rho_alu(tempK)
    if f_slip <= 0:
        f_slip = get_slip_frequency(V_SLIP_MIN)
    return math.sqrt((rho_tempK) / (math.pi * MU0 * f_slip))


def get_eff_plate_depth(f_slip, tempK):
    delta_depth = get_skin_depth_eddy(f_slip, tempK)
    return min(T_PLATE, delta_depth)


def get_plate_eddy_volume(f_slip, tempK):
    depth = get_eff_plate_depth(f_slip, tempK)
    return W_COIL * L_ACTIVE * depth


def get_goodness_factor(f_slip, tempK):
    """Calculate Laithwaite's goodness factor for the LIM.
    
    The goodness factor G determines the operating regime:
    - G << 1: Resistive regime (thrust ∝ slip, low efficiency)
    - G ~ 1: Transitional (peak thrust at optimal slip)
    - G >> 1: Inductive regime (thrust limited by reaction field)
    
    G = (ω_slip * μ₀ * σ * δ_eff * τ_p) / π
    
    For this orbital ring LIM with τ_p = 100m, G is typically very high (500-1000),
    meaning we operate in the inductive regime. This has critical implications:
    - Optimal slip ratio s_opt = 1/G ≈ 0.1%
    - At startup (s = 100%), thrust is only ~0.2% of F_max
    - Thrust improves dramatically as slip decreases during deployment
    
    Reference: Laithwaite (1965), Boldea & Nasar (1976)
    """
    if f_slip <= 0:
        return 0.0
    
    rho = get_rho_alu(tempK)
    sigma = 1.0 / rho
    omega_slip = 2 * math.pi * f_slip
    delta_eff = get_eff_plate_depth(f_slip, tempK)
    
    return (omega_slip * MU0 * sigma * delta_eff * TAU_P) / math.pi


def get_slip_efficiency(slip, G):
    """Calculate slip efficiency factor from goodness factor model.
    
    η_slip = (2sG) / (1 + s²G²)
    
    This unified formula correctly captures both regimes:
    - Low G (resistive): η_slip ≈ 2sG, thrust ∝ slip
    - High G (inductive): η_slip ≈ 2/(sG), thrust drops at high slip
    - Peak at s = 1/G: η_slip = 1 (maximum thrust = F_max)
    
    For this orbital ring LIM at 70K with v_slip = 50 m/s:
    - G ≈ 908
    - s_opt = 0.11%
    - At s = 1% (operational): η_slip ≈ 22%
    - At s = 100% (startup): η_slip ≈ 0.2%
    """
    if slip <= 0 or G <= 0:
        return 0.0
    
    sG = slip * G
    return (2 * sG) / (1 + sG ** 2)


"""
    thrust and power
"""
def get_thrust_force_max(b_plate):
    """Hard upper bound on thrust from magnetic energy density.

    F_max = (B^2 / (2*mu0)) * A  [N]
    where A is the active coupling area under the LIM.
    """
    return (b_plate ** 2) * A_LIM / (2 * MU0)


def get_thrust_power_max(b_plate, v_slip):
    """Hard upper bound on transferable power from magnetic energy flux.

    P_max = (B^2 / (2*mu0)) * A * v_slip  [W]
    """
    return (b_plate ** 2) * A_LIM * v_slip / (2 * MU0)


# =============================================================================
# THRUST MODEL 1: Narrow Plate Eddy Current Model
# =============================================================================
def get_loop_inductance():
    """Calculate inductance of eddy current loop for narrow plate (W << tau_p).
    
    For a rectangular loop of length tau_p and width W:
        L = (mu_0 * tau_p / pi) * ln(tau_p / W)
    """
    return (MU0 * TAU_P / math.pi) * math.log(TAU_P / W_COIL)


def get_loop_resistance(tempK):
    """Calculate resistance of eddy current loop for narrow plate.
    
    Current must return along the length tau_p, through cross-section d x W.
    This return path dominates the loop resistance:
        R = rho * tau_p / (d_eff * W)
    
    where d_eff = min(T_PLATE, skin_depth)
    """
    rho = get_rho_alu(tempK)
    f_slip = get_slip_frequency(V_SLIP_MIN)  # Use minimum for skin depth estimate
    delta = get_skin_depth_eddy(f_slip, tempK)
    d_eff = min(T_PLATE, delta)
    return rho * TAU_P / (d_eff * W_COIL)


def get_f_thrust_model1(f_slip, f_supply, i_peak, tempK=77.0):
    """Thrust Model 1: Narrow Plate Eddy Current Model.
    
    For the orbital ring geometry (W << tau_p), eddy currents form elongated
    loops with high return-path resistance. This model accounts for the
    actual current path geometry.
    
    EDDY CURRENT LOOPS (not related to coil turns N_TURNS):
    --------------------------------------------------------
    The traveling magnetic field B(x) = B0*cos(pi*x/tau_p) creates alternating
    regions of positive and negative field. Eddy currents flow in opposite
    directions under opposite poles, forming loops that span one pole pitch:
    
        ←←←←←← current under negative pole ←←←←←←
        ↓                                        ↑
        ↓ return path (length τ_p)               ↑
        ↓                                        ↑
        →→→→→→ current under positive pole →→→→→→
        
        |←—————————— τ_p = 100m ——————————→|
    
    Number of loops = L_ACTIVE / tau_p = 300m / 100m = 3
    (This has NOTHING to do with N_TURNS = 195)
    
    Physics:
    1. EMF per loop ≈ v_slip * B * W * 2  (factor of 2: both legs contribute)
    2. Loop impedance Z = sqrt(R² + X²) where R dominates for thin plates
    3. Loop current I = EMF / Z
    4. Power dissipated P = N_loops * I² * R
    5. Thrust F = P / v_slip (from power balance)
    
    Key results:
    - Thin plates (d < 2mm): F increases with thickness (R limits current)
    - Thick plates (d > 10mm): F decreases with thickness (inductance limits)
    - Optimal around d = 2mm where R ~ X
    """
    if i_peak <= 0 or f_supply <= 0 or f_slip <= 0:
        return 0.0
    
    # Get magnetic field (this DOES depend on N_TURNS via i_peak)
    B = get_b_plate_peak(i_peak)
    
    # Slip velocity
    v_slip = f_slip * 2 * TAU_P
    if v_slip <= 0:
        return 0.0
    
    # Get material properties at temperature
    rho = get_rho_alu(tempK)
    delta = get_skin_depth_eddy(f_slip, tempK)
    d_eff = min(T_PLATE, delta)
    
    # EMF per current loop
    # Each loop spans one pole pitch. The field varies sinusoidally:
    #   B(x) = B₀ cos(πx/τ_p)
    # At x=0: B = +B₀, at x=τ_p: B = -B₀
    # Both legs across W contribute to EMF (they see opposite fields):
    #   EMF = 2 * v_slip * B_avg * W
    # For sinusoidal variation, B_avg = (2/π) * B₀
    # Combined: EMF = (4/π) * v_slip * B₀ * W ≈ 1.27 * v_slip * B * W
    EMF = (4.0 / math.pi) * v_slip * B * W_COIL
    
    # Loop impedance
    omega = 2 * math.pi * f_slip
    R = rho * TAU_P / (d_eff * W_COIL)
    L = get_loop_inductance()
    X = omega * L
    Z = math.sqrt(R**2 + X**2)
    
    # Eddy current Loop
    I = EMF / Z
    
    # Number of eddy current loops in active region
    # Each loop spans one pole pitch (tau_p), NOT related to coil turns!
    # With PITCH_COUNT = 3, we have L_ACTIVE = 3 * tau_p = 300m
    # This gives N_loops = 3
    N_loops = L_ACTIVE / TAU_P
    
    # Power dissipated in loops
    P = N_loops * I**2 * R
    
    # Thrust from power balance: F * v_slip = P
    F = P / v_slip
    
    # Apply efficiency factor
    return THRUST_EFFICIENCY * F


# =============================================================================
# THRUST MODEL 2: Goodness Factor Model (Laithwaite)
# =============================================================================
def get_f_thrust_model2(f_slip, f_supply, i_peak, tempK=77.0):
    """Thrust Model 2: Goodness Factor Model.
    
    Uses Laithwaite's goodness factor to determine operating regime:
    
    F = F_max * eta_slip
    
    where:
        F_max = B^2 * A / (2*mu_0) = magnetic pressure limit
        eta_slip = (2*s*G) / (1 + s^2*G^2) = slip efficiency
        s = slip ratio
        G = goodness factor = (omega * mu_0 * sigma * d_eff * tau_p) / pi
    
    WARNING: This model assumes W > tau_p (wide plate). For the orbital ring
    geometry where W << tau_p, this model may not be accurate.
    
    At high G (inductive regime):
    - Optimal slip s_opt = 1/G
    - At startup (s=100%): thrust is very low
    - As slip decreases: thrust improves
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    
    b_plate = get_b_plate_peak(i_peak)
    slip = get_slip_f(f_slip, f_supply)
    if slip < 0:
        slip = 0.0
    if slip > 1:
        slip = 1.0
    
    # Ensure minimum slip to avoid division issues
    slip = max(slip, 0.001)
    
    # Calculate goodness factor
    G = get_goodness_factor(f_slip, tempK)
    
    # Calculate slip efficiency
    eta_slip = get_slip_efficiency(slip, G)
    
    # Maximum thrust from magnetic pressure
    F_max = get_thrust_force_max(b_plate)
    
    return THRUST_EFFICIENCY * F_max * eta_slip


# =============================================================================
# THRUST MODEL 3: Slip x Pressure Model (Theoretical Maximum)
# =============================================================================
def get_f_thrust_model3(f_slip, f_supply, i_peak, tempK=77.0):
    """Thrust Model 3: Slip x Pressure Model.
    
    Simple theoretical model:
        F = slip * F_max
    
    where:
        slip = f_slip / f_supply
        F_max = B^2 * A / (2*mu_0) = magnetic pressure limit
    
    This assumes:
    - Resistive regime (perfect phase alignment)
    - No inductance effects
    - Linear relationship between slip and thrust
    
    This is the theoretical maximum thrust for a given slip ratio.
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0

    b_plate = get_b_plate_peak(i_peak)
    slip = get_slip_f(f_slip, f_supply)
    if slip < 0:
        slip = 0.0
    if slip > 1:
        slip = 1.0

    F_max = get_thrust_force_max(b_plate)
    return THRUST_EFFICIENCY * slip * F_max


# =============================================================================
# UNIFIED THRUST FUNCTION (dispatches to selected model)
# =============================================================================
def get_f_thrust(f_slip, f_supply, i_peak, tempK=77.0):
    """Calculate thrust using the selected model.
    
    Model selection via global THRUST_MODEL:
        1 = Narrow plate eddy current model (geometry-corrected)
        2 = Goodness factor model (Laithwaite)
        3 = Slip x pressure model (theoretical maximum)
    
    Use --model=N on command line to select.
    """
    if THRUST_MODEL == 1:
        return get_f_thrust_model1(f_slip, f_supply, i_peak, tempK)
    elif THRUST_MODEL == 2:
        return get_f_thrust_model2(f_slip, f_supply, i_peak, tempK)
    elif THRUST_MODEL == 3:
        return get_f_thrust_model3(f_slip, f_supply, i_peak, tempK)
    else:
        # Default to model 1
        return get_f_thrust_model1(f_slip, f_supply, i_peak, tempK)


def get_thrust_power(thr, vrel):
    return thr * vrel


def get_coil_volts_peak(i_peak_now, f_supply, b_plate):
    if i_peak_now <= 0 or f_supply <= 0:
        return 0.0
    phi_peak = b_plate * A_COIL
    omega = 2 * math.pi * f_supply
    return omega * N_TURNS * phi_peak


# Note: set_r_r removed - replaced by get_rotor_resistance() used in thrust models


"""
    hysteresis
"""
def set_q_hts(i_peak):
    b_coil = get_b_coil_peak(i_peak)
    return b_coil * i_peak * W_TAPE * math.sin(ALPHA_TAPE)


def get_p_hyst_lim(f_supply, i_peak):
    q = set_q_hts(i_peak)
    p_coil_hyst = q * L_HTS_COIL * f_supply        # per coil
    return p_coil_hyst * LIM_PHASES * PITCH_COUNT  # per LIM


"""
    Heat
"""
def get_conduction(d_path, A_cross, T_hot, T_cold, k=K_ALU):
    # gives rate of heat transfer due to thermal conduction (W)
    # T_cold is the cold temperature (K)
    # T_hot is the hot temperature (K)
    # d_path is the path length of conduction (m)
    # A_cross is the cross-sectional area (m^2)
    # k is the thermal conductivity of the material (default aluminium) (W/m*K)
    # - is for direction, i.e. heat is being removed
    return -k * A_cross * (T_cold - T_hot)/d_path


def get_rad(A_heatsink, T_hot=T_N2_HOT, T_cold=T2_LN2_BOIL, em1=EM_ALU, em2=EM_HEAT_SINK):
    # gives rate of radiative heat transfer (W)
    # T_hot & T_cold cold & hot resp. T_hot = reaction plate (K), T_cold = casing 
    # A_heatsink is the effective area of transmission (m^2)
    # em1 & em2 are emissivity of the resp radiating surfaces where 1 = black body (num)
    # em_eff is effective radiation absorption coefficient between em1 and em2
    em_eff = get_em_eff(em1, em2)
    return em_eff * SIGMA * A_heatsink * (T_hot**4 - T_cold**4)


def get_em_eff(em1=EM_ALU, em2=EM_HEAT_SINK):
    return 1/(1/em1 +1/em2 - 1)


def get_convection(A, T_surface=T2_LN2_BOIL, T_fluid=T1_LN2_CRYO, h=H_CONV_LN2):
    # gives rate of heat transfer to cryo-fluid (W)
    # Tsurface is Temperature of the surface (K)
    # Tfl is Temperature of the fluid far away from surface (K)
    # A is surface area (m^2)
    # h is Convective heat transfer coefficient (W/m²·K)
    return h * A * (T_surface - T_fluid)


def get_LN2_mass_flow_rate(Q_heat, T1=T1_LN2_CRYO):
    # gives minimum mass flow rate of LN2 liquid based on W of heat dissipated (kg/s)
    # Q_heat is the heat generated by the coil (W)
    # T1 is starting temp of LN2 (K)
    return Q_heat / (C_P_LN2 * (T2_LN2_BOIL - T1) + L_V_LN2)


def get_cop_actual(T_cold, T_hot=T_N2_HOT, eff=CRYO_EFF):
    # effective coefficient of performance (COP) of the cryogenic system (num)
    # eff is the overall efficiency of the cryogenic system. 0.05 = 1/20 ratio of performance (num)
    if T_hot == T_cold:
        T_hot = T_cold + 1
    return (T_cold/(T_hot - T_cold))*eff


def get_p_cryo(p_heat, T_hot=T_N2_HOT, T_cold=T2_LN2_BOIL, eff=CRYO_EFF):
    # gives power needed to opperate cryogenic system (W)
    # Q is the amount of heat (W)
    # T1 & T2 cold & hot resp. (K)
    # eff is cryo efficiency factor. 1=100% efficient. 0.05 -> need 20 W for 1 W of heat (num)
    cop = get_cop_actual(T_cold, T_hot, eff)
    return p_heat/cop


def get_v_rel_min(v_rel):
    """
        The windings or the coils are really only a few mmm wide and
        center of the coil act as both a heat sink and a source of the
        eddy currents. The code does not really account for this, so the
        early temperature sikes are not accurate. We fudge the initial 
        value of v_rel to v_rel_min only 
    """
    return  max(V_REL_MIN, v_rel)


def get_p_load(p_heat):
    duty = L_ACTIVE / HEAT_SINK_L    # shorter heat sink mean higher ave reaction plate temp.
    # the coils are 99.99% open space, so heatsink extends under them.
    p_ave = p_heat * duty
    p_load = p_ave + Q_ABS_LIM       # Q_ABS_LIM is the daytime heat from the environment, daytime is what matters
    return p_load


def get_heatsink_area(p_heat, T_hot=T_N2_HOT, T_cold=T2_LN2_BOIL, em1=EM_ALU, em2=EM_HEAT_SINK):
    # returns area of heatsink (m^2) for a single LIM
    # T_cold = casing, T_hot = reaction plate (K)
    # A is cross-sectional area (m^2)
    # em1 & em2 are emissivity of the resp radiating surfaces where 1 = black body (num)
    # em_eff is effective radiation absorption coefficient between em1 and em2
    #print(p_heat)
    em_eff = get_em_eff(em1, em2)
    if T_hot == T_cold:
        T_hot = T_hot + 1
    return p_heat / (em_eff * SIGMA * (T_hot**4 - T_cold**4))


def get_T_min_ambient(A, p_heat, em1=EM_ALU, em2=EM_HEAT_SINK):
    # returns the minimum ambient temperature inside ring casing
    em_eff = get_em_eff(em1, em2)
    return (p_heat / (A * em_eff * SIGMA) + T2_LN2_BOIL**4)**0.25


def get_plate_temp(v_rel, p_eddy):
    # equilibrium temperature of the reaction plate
    v_rel_adjusted = get_v_rel_min(v_rel) 
    time_lim = L_ACTIVE / v_rel_adjusted
    #heatsink_length = min(HEAT_SINK_L, LIM_SPACING - L_ACTIVE)    # the heatsink should never be longer than the distance between the LIMs
    heatsink_length = LIM_SPACING               # the coils are 99.99% open space
    time_gap = heatsink_length / v_rel_adjusted
    time_cycle = time_gap + time_lim
    em_eff = get_em_eff()
    # this is get_rad solved for T_hot with the ratio of lim_time/cycle_time added
    return ((p_eddy * time_lim) / (em_eff * SIGMA * (2 * A_LIM) * time_cycle) + T2_LN2_BOIL**4)**(0.25) # (K)
    

def get_delta_t_lim(v_rel, f_slip, p_eddy, tempK):
    # returns the increase in reaction plate temp as it passes thru LIM.
    v_rel_adjusted = get_v_rel_min(v_rel) 
    time_lim = L_ACTIVE / v_rel_adjusted
    v_plate_eddy = get_plate_eddy_volume(f_slip, tempK)
    m_plate_eddy = v_plate_eddy * RHO_ALU_M3
    return (p_eddy * time_lim) / (m_plate_eddy * C_P_ALU) # (K)

def get_rho_alu(tempK):
    # calculate the current resistivity of aluminium based on tempreature in kelvin
    rho_alu_293K = 2.65E-8
    alpha_alu_e = 3.663E-3
    rho_alu = rho_alu_293K * (1 + alpha_alu_e * (tempK - 293))
    if rho_alu < RHO_ALU_E_70K:
        rho_alu = RHO_ALU_E_70K
    return rho_alu


"""
    helper equations
"""
def get_total_for_ring(x):
    # units are passed thru unchanged
    return x * LIM_SITES


def get_progress(vcase):
    return round((1 - (vcase - GEO_250_ORBIT) / (V_ORBIT - GEO_250_ORBIT)) * 100, 2)


def make_month_ticks(data_list, total_time):
    """Generate tick positions and labels based on deployment time in months."""
    months = total_time / MONTH
    num_ticks = int(months) + 1  # No cap - show all months
    list_len = len(data_list)
    
    tick_positions = [list_len * i / max(months, 1) for i in range(num_ticks)]
    tick_labels = [str(i) for i in range(num_ticks)]
    
    return tick_positions, tick_labels


def annotate_final(data_list, unit="", fmt=".1f"):
    """Add annotation showing the final value on the plot."""
    if not data_list:
        return
    x = len(data_list) - 1
    y = data_list[-1]
    label = f"{y:{fmt}} {unit}".strip()
    plt.annotate(
        label,
        xy=(x, y),
        xytext=(-40, 10),
        textcoords="offset points",
        fontsize=10,
        ha="left",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0),
    )



"""
    # -------------------------------------------------
    # MAIN LOOP. CALCULATE DEPLOYMENT TIME
    # -------------------------------------------------
"""
def get_deployment_time(v_slip, i_peak_now):
    # loop variables
    dt = DT1        # initial number of seconds per loop
    i_target = i_peak_now # set initial value of i_min

    print("VOLTS_MAX: ", VOLTS_MAX)
    print(f"THRUST_MODEL: {THRUST_MODEL} ", end="")
    if THRUST_MODEL == 1:
        print("(Narrow plate eddy current - geometry corrected)")
    elif THRUST_MODEL == 2:
        print("(Goodness factor - Laithwaite)")
    elif THRUST_MODEL == 3:
        print("(Slip x pressure - theoretical maximum)")
    print(f"T_PLATE: {T_PLATE * 1000:.1f} mm")

    # variable with starting conditions
    vcable = V_ORBIT
    vcasing = V_ORBIT
    time = 0
    yrs = 0
    mts = 0
    second = 0
    minute = 60
    hours = 1
    days = 1
    months = 1
    max_min = {         # [value, time]
        "volts_max":          [0, 0, "",0],
        "current_max":        [0, 0, "",0],
        "v_slip_max":         [0, 0, "",0],
        "p_cryo_max":         [0, 0, "",0],    
        "p_eddy_max":         [0, 0, "",0],
        "p_hyst_max":         [0, 0, "",0],
        "p_heat_max":         [0, 0, "",0],
        "f_slip_max":         [0, 0, "",0],
        "f_supply_max":       [0, 0, "",0],
        "b_field_max":        [0, 0, "",0],
        "thrust_max":         [0, 0, "",0],
        "p_thrust_max":       [0, 0, "",0],
        "volts_lim_max":      [0, 0, "",0],
        "delta_t_max":        [0, 0, "",0],
        "plate_t_max":        [0, 0, "",0],
        "skin_d_max":         [0, 0, "",0],
        "skin_d_min":         [100, 0, "",0],
        "min_heatsink_area":  [0, 0, "",0],
        "min_ambient_T":      [0, 0, "",0],
        "min_cryo_area":      [0, 0, "",0],
        "p_heat_load":        [0, 0, "", 0],
        "alu_temp_out":       [0, 0, "", 0],
        "lim_power_max":      [0, 0, "",0],
        "site_power_max":     [0, 0, "",0],
    }
    power_track = []
    p_cryo = 0
    # Energy accumulators (Joules)
    E_site_ke = 0.0      # Kinetic energy per site (thrust power only)
    E_total_ke = 0.0     # Total kinetic energy all sites
    count = 0
    sample_time = 0
    skip = 1             # record first round, then start skip
    sample_time_max = SAMPLE_TIME_MAX

    param_str1 = f"{SAMPLE_TIME_STR}\n{PARAM_STR1}"
    exit_msg = "PASSED"
    make_graphs = MAKE_GRAPHS
    alu_temp_out = 70.0  # the entire orbital ring starts at cryogenic temperature

    # --- Controller smoothing (does not change physics, only actuator behavior) ---
    CTRL_TAU_SLIP = 15 * DAY          # time constant for v_slip response (tune: 3–12 hr)
    DVSLIP_MAX_PER_S = 0.0001         # max v_slip slew rate (m/s^2). tune: 0.005–0.05
                                     # 0.01 = max 36 m/s change per hour

    v_slip_cmd = v_slip             # controller state (commanded v_slip)

    DI_MAX_PER_S = 0.1   # max current slew rate (A/s). tune: 0.5–10 A/s

    i_peak_prev = i_peak_now
    
    # Initialize power tracking for predictive control
    lim_site_power = 0.0  # Will be calculated properly in first iteration

    """
        The two control variables are i_peak_now and v_slip. i_peak_now 
        is the current flowing into the LIM coils, and v_slip controls 
        the power supply frequency (f_supply). If v_slip is 0, then f_supply 
        exactly matches v_rel (the relative velocity between the cable and casing).

        Decreasing f_supply below v_slip = 0 (a negative value for v_slip) 
        would function as magnetic braking. Excessive v_slip and/or i_peak_now 
        will spike the induced voltage in the superconducting coils, causing 
        a quenching event or worse.
    """
    while vcasing > GEO_250_ORBIT:
        # time varying parameters that need to be set per loop
        #
        # IMPORTANT: We enforce voltage and power limits *before* applying this timestep's acceleration.
        # This prevents the early-time power spikes from "spending" illegal power for one step.

        # --- Controller inner loop (a few quick iterations for self-consistency) ---

        for _ctrl in range(10):
            # Relative velocity between cable and casing (your v_rel definition)
            v_rel = get_v_rel(vcable, vcasing)

            # ============================================================
            # PREDICTIVE POWER-MARGIN CONTROL
            # ============================================================
            # Instead of reacting to current reduction, we proactively adjust
            # v_slip based on how close we are to the power limit.
            #
            # Key insight: p_heat drives p_cryo which dominates site_power.
            # By monitoring power margin directly, we can anticipate and
            # smoothly adjust before hitting hard limits.
            # ============================================================
            
            # Calculate current power margin (1.0 = plenty of headroom, 0.0 = at limit)
            power_margin = 1.0 - (lim_site_power / MAX_SITE_POWER)
            
            # Define margin thresholds for proactive control
            MARGIN_COMFORTABLE = 0.15  # Above this: use optimal slip ratio
            MARGIN_TIGHT = 0.05        # Below this: aggressively increase v_slip
            
            if power_margin >= MARGIN_COMFORTABLE:
                # Plenty of headroom - use optimal slip ratio for thrust
                s_tgt = SLIP_RATIO_NORMAL
            elif power_margin <= MARGIN_TIGHT:
                # Very tight - maximize v_slip to reduce heat
                s_tgt = SLIP_RATIO_REDUCED
            else:
                # Linear interpolation based on power margin
                t = (power_margin - MARGIN_TIGHT) / (MARGIN_COMFORTABLE - MARGIN_TIGHT)
                s_tgt = SLIP_RATIO_REDUCED + t * (SLIP_RATIO_NORMAL - SLIP_RATIO_REDUCED)
            
            # Also factor in current headroom (secondary consideration)
            current_ratio = i_peak_now / I_TARGET
            if current_ratio < 0.7:
                # If current is significantly reduced, blend toward higher slip
                current_blend = (0.7 - current_ratio) / 0.2  # 0 at 70%, 1 at 50%
                current_blend = max(0.0, min(1.0, current_blend))
                s_tgt = s_tgt + current_blend * (SLIP_RATIO_REDUCED - s_tgt)

            # Clamp slip ratio to valid range
            s_tgt = max(1e-6, min(0.999, s_tgt))

            # Compute v_slip from slip ratio: s = v_slip/(v_rel + v_slip) => v_slip = s/(1-s) * v_rel
            v_slip_target = (s_tgt / (1.0 - s_tgt)) * max(v_rel, 0.0)
            
            # Apply floor (ensures thrust at low v_rel) and ceiling (limits losses)
            #v_slip = max(V_SLIP_MIN, min(V_SLIP_MAX, v_slip_target))
            # Compute v_slip from slip ratio target
            v_slip_target = (s_tgt / (1.0 - s_tgt)) * max(v_rel, 0.0)

            # First apply bounds to the target (still the same physics constraints)
            v_slip_target = max(V_SLIP_MIN, min(V_SLIP_MAX, v_slip_target))

            # --- Smooth actuator: first-order lag + slew-rate limit ---
            # Low-pass coefficient that remains stable across big dt changes
            alpha = 1.0 - math.exp(-dt / max(CTRL_TAU_SLIP, 1e-9))

            # Low-pass toward target
            v_slip_cmd = v_slip_cmd + alpha * (v_slip_target - v_slip_cmd)

            # Slew-rate limit (prevents one-step jumps when dt is large)
            max_step = DVSLIP_MAX_PER_S * dt
            v_slip_cmd = max(v_slip - max_step, min(v_slip + max_step, v_slip_cmd))

            v_slip = v_slip_cmd


            v_wave = v_rel + v_slip
            slip = get_slip(v_slip, v_rel)

            # Frequencies
            f_slip = get_slip_frequency(v_slip)
            f_supply = get_supply_frequency(v_wave)

            # Magnetic field at plate
            b_plate = get_b_plate_peak(i_peak_now)
            b_coil = get_b_coil_peak(i_peak_now)

            # Thrust and thrust power (per LIM)
            # Pass alu_temp_out for temperature-dependent calculations
            thrust = get_f_thrust(f_slip, f_supply, i_peak_now, alu_temp_out)
            p_thrust = get_thrust_power(thrust, v_rel)

            # Induced coil voltage (RMS) from flux linkage
            volts_lim = get_coil_volts_peak(i_peak_now, f_supply, b_coil)

            # Secondary and hysteresis losses (per LIM)
            p_eddy = get_plate_eddy_loss_from_thrust(v_slip, thrust)
            # p_eddy = get_plate_eddy_loss(v_slip, i_peak_now, alu_temp_out)


            # Hard clamp: secondary power cannot exceed magnetic energy flux for same B
            p_eddy_cap = get_thrust_power_max(b_plate, v_slip)
            #if p_eddy > p_eddy_cap:
            #    p_eddy = p_eddy_cap

            p_hyst = get_p_hyst_lim(f_supply, i_peak_now)

            # Heating (per LIM) used for thermal and cryo sizing
            # p_heat includes both eddy currents (in reaction plate) and hysteresis (in HTS)
            p_heat = p_eddy + p_hyst

            # Reaction plate temperature estimate
            temp_plate_ave = get_plate_temp(v_rel, p_eddy)
            delta_temp_lim = get_delta_t_lim(v_rel, f_slip, p_eddy, alu_temp_out)
            alu_temp_out = temp_plate_ave + delta_temp_lim

            # Heat dissipation surfaces (per LIM)
            # Note: p_heat_load is for REACTION PLATE heatsink only (p_eddy)
            # HTS hysteresis (p_hyst) goes directly to cryo system, not this heatsink
            p_heat_load = get_p_load(p_eddy)
            min_heatsink_area = get_heatsink_area(p_heat_load, alu_temp_out)
            min_T_ambient = get_T_min_ambient((LIMS_PER_SITE * MAX_HEATSINK_AREA), p_heat_load)  # one heatsink per LIM
            
            # Cryo radiator calculation:
            # The cryo system removes heat_t from cold side (77K) using p_cryo electrical power
            # The radiator must reject Q_hot = heat_t + p_cryo to space
            # This is Q_cold * (1 + 1/COP) = Q_cold * (COP + 1) / COP
            if v_rel != 0:
                heat_t = LIMS_PER_SITE * p_heat + Q_ABS_LIM  # Cold-side heat from all LIMs at site
                p_cryo_temp = get_p_cryo(heat_t)
                Q_hot_cryo = heat_t + p_cryo_temp  # Total heat to radiate to space
                min_cryo_area = get_heatsink_area(Q_hot_cryo, T_N2_HOT, T_FREE_SPACE, EM_HEAT_SINK, 1.0)
            else:
                min_cryo_area = 0.0

            # Total LIM-side power bookkeeping (per LIM, then per site)
            p_total_lim = (p_thrust + p_heat_load + p_hyst) / max(LIM_EFF, 1e-6)

            # Cryo power is per LIM site; include external absorbed heat via p_heat_load already.
            if v_rel != 0:
                heat_t = LIMS_PER_SITE * p_heat + Q_ABS_LIM
                p_cryo = get_p_cryo(heat_t)
            else:
                p_cryo = 0.0

            lim_site_power = (LIMS_PER_SITE * p_total_lim + p_cryo) / max(INV_EFF, 1e-6)

            # --- Enforce limits by scaling current/slip ---
            changed = False

            # Voltage limit: volts_lim scales ~ I (because B ~ I and Phi ~ B)
            if volts_lim > VOLTS_MAX:
                scale_v = VOLTS_MAX / max(volts_lim, 1.0)
                # keep scaling sane
                scale_v = max(0.0, min(1.0, scale_v))
                i_peak_now = max(I_MIN, i_peak_now * scale_v)
                changed = True

            # Thermal limit: conservative proportional throttling
            if (temp_plate_ave > T_MAX_REACTION_PLATE) or ((temp_plate_ave + delta_temp_lim) > T_MAX_REACTION_PLATE):
                i_peak_now = max(I_MIN, i_peak_now * 0.95)
                v_slip = max(V_SLIP_MIN, v_slip * 0.95)
                changed = True

            # Site power limit: lim_site_power scales ~ I^2 (dominant terms), so use sqrt scaling
            if lim_site_power > MAX_SITE_POWER:
                scale_p = math.sqrt(MAX_SITE_POWER / max(lim_site_power, 1.0))
                scale_p = max(0.0, min(1.0, scale_p))
                i_new = max(I_MIN, i_peak_now * scale_p)
                if i_new < i_peak_now:
                    i_peak_now = i_new
                    changed = True
                # If current is already at minimum and still too much, reduce slip
                if i_peak_now <= I_MIN * 1.0001 and v_slip > V_SLIP_MIN:
                    v_slip = max(V_SLIP_MIN, v_slip * 0.95)
                    changed = True

            # Current slew limit (actuator realism)
            i_step = DI_MAX_PER_S * dt
            i_peak_now = max(i_peak_prev - i_step, min(i_peak_prev + i_step, i_peak_now))
            i_peak_prev = i_peak_now


            # If we are comfortably below all limits, ramp current upward to use available power headroom.
            # Without this, the controller can get "stuck" at a low current set early in the run.
            if (not changed
                    and i_peak_now < I_TARGET
                    and lim_site_power < POWER_HEADROOM * MAX_SITE_POWER
                    and volts_lim < POWER_HEADROOM * VOLTS_MAX
                    and temp_plate_ave < POWER_HEADROOM * T_MAX_REACTION_PLATE):
                i_peak_now = min(I_TARGET, i_peak_now * CURRENT_UPRATE)
                changed = True

            if not changed:
                break

        # --- Apply this timestep's acceleration using the LIMIT-COMPLIANT thrust ---
        vcable = get_v_cable(vcable, thrust, dt)
        vcasing = get_v_casing(vcasing, thrust, dt)

        # Accumulate kinetic energy (thrust power only, no losses)
        # Each site has LIMS_PER_SITE LIMs
        E_site_ke += p_thrust * LIMS_PER_SITE * dt
        E_total_ke += p_thrust * LIM_SITES * LIMS_PER_SITE * dt

        # Startup ramp (only if comfortably below limits)
        if time < HR:  # First hour
            if (temp_plate_ave < T_MAX_REACTION_PLATE * 0.8 and volts_lim < VOLTS_MAX * 0.8 and lim_site_power < MAX_SITE_POWER * 0.8):
                i_peak_now += (I_TARGET - i_peak_now) * 0.01
                # v_slip is controlled by slip ratio; do not force-ramp it here
        # Starting values for table, ignore 0. Only run once.
        if time == 1: 
            power_track.append(["START", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])

        
        skin_depth_eff = get_eff_plate_depth(f_slip, alu_temp_out)
        assert skin_depth_eff > 0.0
        skin_depth_calc = get_skin_depth_eddy(f_slip, alu_temp_out)
        assert skin_depth_calc > 0.0
        
        if time > sample_time and time < sample_time_max:
            sample_time += SKIP
            count += 1
            
            list_i_peak.append(i_peak_now)
            list_volts.append(volts_lim )
            list_v_slip.append(v_slip)
            list_slip.append(slip*100)
            list_f_slip.append(f_slip)
            list_thrust.append(thrust)
            list_thrust_power.append(p_thrust)
            list_v_rel.append(v_rel)
            list_b_peak.append(b_plate)
            list_p_eddy.append(p_eddy)
            list_skin_depth_eff.append(skin_depth_eff*1000)
            list_skin_depth_calc.append(skin_depth_calc*1000)
            list_p_cryo.append(p_cryo)
            list_p_hyst.append(p_hyst)
            list_p_lim.append(p_total_lim)
            list_p_lim_site.append(lim_site_power)
            list_temp_plate_ave.append(temp_plate_ave)
            list_E_site_ke.append(E_site_ke)
            list_E_total_ke.append(E_total_ke)

        # collect data at various early intervals
        if time > second and time < 60:
            power_track.append([f"{round(time)} sec", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            second += 5
        if time > minute and time < HR:
            power_track.append([f"{round(time/60)} min", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            if minute == 60:
                minute += 60*4
            else:
                minute += 60*5
        if time > hours * HR and time < 24 * HR:
            power_track.append([f"{round(time/HR)} hrs", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            hours += 1
        if time > days * DAY and time < 31 * DAY:
            power_track.append([f"{round(time/DAY)} day", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            days += 1

        # show progress to stdout and collect data monthly
        if time > months * MONTH:
            power_track.append([f"{round(time/MONTH)} mth", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            months += 1
            # show progress to screen once a month
            #sys.stdout.write(".")
            sys.stdout.flush()

        # show yearly progress to stdout
        if time > mts * MONTH:
            if mts < 10:
                mthstr = str(mts)+"  "
            else:
                mthstr = str(mts)+" "
            sys.stdout.write(mthstr)
            prog = get_progress(vcasing)
            if mts == 0:
                print(f"dt =  {dt}")
                print(f"Months| Prog | Volts  | I_Peak  | V_Slip  | Slip |F_Supply|  Eddy |  Hyst | Cryo Power|Site Power| Thrust")
                print(f"{mthstr}mts, {prog:.1f}%, {volts_lim:.2f} V,  {i_peak_now:.2f} A,  {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:.2f} Hz, {p_eddy:.2f} W, {p_hyst:.2f} W, {p_cryo/1e6:.3f} MW, {lim_site_power/1e6:.3f} MW, {thrust:.3f} N")
            else:
                print(f"mts, {prog:.0f}%, {volts_lim:.0f} V, {i_peak_now:.2f} A, {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:2.2f} Hz, {p_eddy:.0f} W, {p_hyst:.0f} W, {p_cryo/1e6:.3f} MW, {lim_site_power/1e6:.1f} MW, {thrust:.0f} N")
            mts += 1
 
        # collect min-max data
        if v_rel != 0 and time > DAY: # things are out of wack at the beginning and a real deployment would account for it.
            if max_min["volts_max"][3] < volts_lim:
                max_min["volts_max"] = [round(volts_lim,2), time, "V", volts_lim]

            if max_min["current_max"][3] < i_peak_now:
                max_min["current_max"] = [round(i_peak_now,2), time, "A", i_peak_now]

            if max_min["v_slip_max"][3] < v_slip:
                max_min["v_slip_max"] = [round(v_slip,6), time, "m/s", v_slip]

            if max_min["f_slip_max"][3] < f_slip:
                max_min["f_slip_max"] = [round(f_slip,6), time, "Hz", f_slip]

            if max_min["f_supply_max"][3] < f_supply:
                max_min["f_supply_max"] = [round(f_supply,6), time, "Hz", f_slip]

            if max_min["b_field_max"][3] < b_plate:
                max_min["b_field_max"] = [round(b_plate*1e3,6), time, "mT", b_plate]

            if max_min["thrust_max"][3] < thrust:
                max_min["thrust_max"] = [round(thrust,2), time, "N", thrust]

            if max_min["p_thrust_max"][3] < p_thrust:
                max_min["p_thrust_max"] = [round(p_thrust), time, "W", p_thrust]

            if max_min["volts_lim_max"][3] < volts_lim:
                max_min["volts_lim_max"] = [round(volts_lim,2), time, "V", volts_lim]

            if max_min["p_eddy_max"][3] < p_eddy:
                max_min["p_eddy_max"] = [round(p_eddy,2), time, "W", p_eddy]

            if max_min["p_hyst_max"][3] < p_hyst:
                max_min["p_hyst_max"] =  [round(p_hyst,2), time, "W", p_hyst]

            if max_min["plate_t_max"][3] < temp_plate_ave:
                max_min["plate_t_max"] = [round(temp_plate_ave,2), time, "K", temp_plate_ave]

            if max_min["delta_t_max"][3] < delta_temp_lim:
                max_min["delta_t_max"] = [round(delta_temp_lim,6), time, "K", delta_temp_lim]

            if max_min["p_heat_max"][3] < p_heat:
                max_min["p_heat_max"] = [round(p_heat,2), time, "W", p_heat]

            if max_min["lim_power_max"][3] < p_total_lim:
                max_min["lim_power_max"] = [round(p_total_lim/1e6,3), time, "MW", p_total_lim]

            if max_min["site_power_max"][3] < lim_site_power:
                max_min["site_power_max"] = [round(lim_site_power/1e6,3), time, "MW", lim_site_power]

            if max_min["min_heatsink_area"][3] < min_heatsink_area:
                max_min["min_heatsink_area"] = [round(min_heatsink_area), time, "m²/LIM", min_heatsink_area]

            if max_min["min_ambient_T"][3] < min_T_ambient:
                max_min["min_ambient_T"] = [round(min_T_ambient), time, "K", min_T_ambient]

            if max_min["min_cryo_area"][3] < min_cryo_area:
                max_min["min_cryo_area"] = [round(min_cryo_area), time, "m²", min_cryo_area]

            if max_min["p_heat_load"][3] < p_heat_load:
                max_min["p_heat_load"] = [round(p_heat_load), time, "W", p_heat_load]

            if max_min["alu_temp_out"][3] < alu_temp_out:
                max_min["alu_temp_out"] = [round(alu_temp_out), time, "W", alu_temp_out]

            if max_min["p_cryo_max"][3] < p_cryo:                  
                max_min["p_cryo_max"] = [round(p_cryo/1e6,2), time, "MW", p_cryo]

            if max_min["skin_d_max"][3] < skin_depth_eff:
                max_min["skin_d_max"] = [round(skin_depth_eff*1000,2), time, "mm", skin_depth_eff]

            if max_min["skin_d_min"][3] > skin_depth_eff:
                max_min["skin_d_min"] = [round(skin_depth_eff*1000,2), time, "mm", skin_depth_eff]


        if volts_lim > VOLTS_MAX * 1.1:
            exit_msg = f"FAIL over voltage limit: {volts_lim:.0f} V"
            make_graphs = False
            break
        if lim_site_power > MAX_SITE_POWER * 1.1:
            exit_msg = f"FAIL max site power exceeded: {lim_site_power:.0f} W"
            make_graphs = False
            break
        if TAU_P * PITCH_COUNT > LIM_SPACING:
            exit_msg = f"FAIL The LIM ia longer than the LIM_SPACING: LIM={TAU_P * PITCH_COUNT} m, SPACING={LIM_SPACING} m\n"
            make_graphs = False
            break
        # increment loop variable for next loop and break it time exceeded
        if time < DAY:
            dt = DT1
        elif time < sample_time_max + 60:
            dt = DT2
        else:
            dt = DT3
        time += dt
        if time > 100 * YR:
            print(">>>>>>>>> DEPLOYMENT TIME EXCEEDED.<<<<<<<<<<<")
            break

    # Cryo radiator must reject Q_hot = Q_cold + W = Q_cold * (1 + 1/COP)
    # where Q_cold = heat removed from cold side, W = electrical power
    cop = get_cop_actual(T2_LN2_BOIL, T_N2_HOT)
    if cop > 0:
        Q_hot_factor = 1 + 1/cop  # Multiply cold-side heat by this to get Q_hot
    else:
        Q_hot_factor = 60  # Fallback
    
    # p_cryo is electrical power, so Q_cold = p_cryo * COP, and Q_hot = p_cryo * (COP + 1)
    p_cryo_max = max_min["p_cryo_max"][3]
    Q_hot_max = p_cryo_max * (cop + 1) if cop > 0 else p_cryo_max * 60
    cryo_radiator_size_min = get_heatsink_area(Q_hot_max, T_N2_HOT, T_FREE_SPACE, EM_HEAT_SINK, 1.0)
    cryo_radiator_width = cryo_radiator_size_min / LIM_SPACING
    
    print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(f"{exit_msg}\t\t{time/YR:.2f} yrs")
    # post loop data collection
    power_track.append([f"{round(time/YR,2)} yrs", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,3), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_cryo/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
    power_track.append(["Time", "V_Shell", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "Cryo Power", "Site Power"])

    str2 = f"Deployment time: {time / DAY:.2f} days, {time / YR:.2f} years"
    str3 = f"Cable velocity: {vcable:.2f} m/s"
    str4 = f"Casing velocity: {vcasing:.2f} m/s"

    # post loop wrap-up. show data.
    lines = ["\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", exit_msg+"\n", str2+"\n", str3+"\n", str4+"\n"]
    lines.append("\n--------------------------------------------------------------------\n")
    for key, value in PARAM_LIST.items():
        str1 = (f"\t{key:20}\t{value:8}\n")
        lines.append(str1)

    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")
    # print out max_min data
    for key, value in max_min.items():
        if value[1] < 60:
            t = f"{value[1]:8.2f} seconds"
        elif value[1] < HR:
            t = f"{value[1]/60:8.2f} minutes"
        elif value[1] < DAY:
            t = f"{value[1]/HR:8.2f} hours"
        elif value[1] < MONTH:
            t = f"{value[1]/DAY:8.2f} days"
        elif value[1] < YR:
            t = f"{value[1]/MONTH:8.2f} months"
        else:
            t = f"{value[1]/YR:8.2f} years"
        str1 = (f"\t{key:18} {value[0]:18.2f} {value[2]:8} {t}")
        print(str1)
        lines.append(str1+"\n")
    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")
   
    print("")
    print(tabulate.tabulate(power_track, headers=["V_Shell", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "Cryo Power", "Site Power"]))
    lines.append(tabulate.tabulate(power_track, headers=["V_Shell", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "Cryo Power", "Site Power"]))
    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")

 
    print("")

    print(str2)
    print(str3)
    print(str4)
    print("Count = ", count)
    print(f"Site KE ({LIMS_PER_SITE} LIMs): {E_site_ke/1e12:.2f} TJ")
    print(f"Total KE (all sites): {E_total_ke/1e18:.4f} EJ")
    print(f"Minimum cryo radiator width: {cryo_radiator_width:.1f} m (area: {cryo_radiator_size_min:.0f} m²)")

    if exit_msg != "PASSED":
        print(f"{exit_msg}. time in seconds: {time:.3f} s, years: {time/YR:.2f}")

    param_str1 = param_str1 + f", Model: {THRUST_MODEL}, F_Supply: {round(f_supply,1)} Hz, Deployment time: {round(time / YR, 2)} years."

    tick_count = round(time / YR)

    if WRITE_FILE:
        with open(f"./output/_orbital_ring_model{THRUST_MODEL}.txt", "a") as file:
            file.writelines(lines)

    return [tick_count, f"{param_str1} ", make_graphs, time]

def main() -> None:
    global THRUST_MODEL
    
    show = []
    if len(sys.argv) > 1:
        show = sys.argv
        
        # Check for help
        if "--help" in sys.argv or "-h" in sys.argv:
            print("""
Orbital Ring Deployment Simulation
===================================

Usage: python legacy_power_lim.py [options] [graphs]

Thrust Model Selection:
  --model=1    Narrow plate eddy current model (geometry-corrected for W << tau_p)
  --model=2    Goodness factor model (Laithwaite, assumes W > tau_p)
  --model=3    Slip x pressure model (theoretical maximum, F = slip * F_max)

Graph Options:
  all          Show all graphs
  current      Current (Amps)
  volts        Voltage
  v_slip       Slip velocity
  thrust       Thrust
  p_thrust     Thrust power
  p_eddy       Eddy current losses
  power        Site power used
  plate_temp   Reaction plate temperature
  skin         Skin depth
  slip         Slip ratio
  f_slip       Slip frequency
  v_rel        Relative velocity
  hyst         Hysteresis losses
  cryo         Cryogenic power
  ke_site      Site kinetic energy
  ke_all       Total kinetic energy

Examples:
  python legacy_power_lim.py --model=1 thrust power
  python legacy_power_lim.py --model=2 all
  python legacy_power_lim.py --model=3
            """)
            return
        
        # Parse --model=N argument
        for arg in sys.argv:
            if arg.startswith("--model="):
                try:
                    model_num = int(arg.split("=")[1])
                    if model_num in [1, 2, 3]:
                        THRUST_MODEL = model_num
                        print(f"Using thrust model {THRUST_MODEL}:")
                        if THRUST_MODEL == 1:
                            print("  Model 1: Narrow plate eddy current (geometry-corrected)")
                        elif THRUST_MODEL == 2:
                            print("  Model 2: Goodness factor (Laithwaite)")
                        elif THRUST_MODEL == 3:
                            print("  Model 3: Slip x pressure (theoretical maximum)")
                    else:
                        print(f"Invalid model number: {model_num}. Using default model 1.")
                        THRUST_MODEL = 1
                except ValueError:
                    print(f"Invalid model argument: {arg}. Using default model 1.")
                    THRUST_MODEL = 1
    
    print(f"\n*** THRUST MODEL: {THRUST_MODEL} ***\n")
    
    v_slip = V_SLIP_MAX
    i_peaks = I_MIN
    param = get_deployment_time(v_slip, i_peaks)



    if param[2]:
        param_str = param[1]
        total_time = param[3]

        x_lable = SAMPLE_TIME_STR
        
        #"""

        # 1 Current in Amps
        if "current" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_i_peak, total_time)
            plt.scatter(range(len(list_i_peak)), list_i_peak, c="blue")
            plt.xlabel("Months")
            plt.ylabel("Current in Amps")
            plt.legend(["Amps"])
            plt.title(f"Current {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_i_peak, unit="A", fmt=".0f")
            plt.show()



        # 2 Volts
        if "volts" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_volts, total_time)
            plt.scatter(range(len(list_volts)), list_volts, c="red")
            plt.xlabel("Months")
            plt.xlabel("Voltage")
            plt.ylabel("Voltage")
            plt.legend(["Volts"])
            plt.title(f"Volts {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_volts, unit="V", fmt=".0f")
            plt.show()


        # 3 V_Slip Velocity
        if "v_slip" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_v_slip, total_time)
            plt.scatter(range(len(list_v_slip)), list_v_slip, c="green")
            plt.xlabel("Slip Velocity")
            plt.ylabel("Slip Velocity in m/s")
            plt.legend(["Slip Velocity"])
            plt.title(f"V_Slip {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_v_slip, unit="m/s", fmt=".0f")
            plt.show()
        

        # 4 Skin Calc & Eff
        if "skin" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_skin_depth_calc, total_time)
            plt.scatter(range(len(list_skin_depth_calc)), list_skin_depth_calc, c="darkgreen")
            plt.scatter(range(len(list_skin_depth_eff)), list_skin_depth_eff, c="magenta")
            plt.xlabel("Skin Depth (Calc & Eff)")
            plt.ylabel("Skin Depth/Plate Thickness in mm")
            #    plt.yticks([0,1,2,3,4,5,6])
            plt.legend(["Skin Depth Calc", "Skin Depth Eff"])
            plt.title(f"Effective Reaction Plate Plate Thickness in mm {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_skin_depth_eff, unit="mm", fmt=".1f")
            annotate_final(list_skin_depth_calc, unit="mm", fmt=".1f")
            plt.show()


        # 5 Thrust
        if "thrust" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_thrust, total_time)
            plt.scatter(range(len(list_thrust)), list_thrust, c="purple")
            plt.xlabel("Thrust")
            plt.ylabel("Thrust in Newtons")
            plt.legend(["Thrust"])
            plt.title(f"Thrust {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_thrust, unit="N", fmt=".0f")
            plt.show()

        
        # 6 P_Eddy
        if "p_eddy" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_p_eddy, total_time)
            plt.scatter(range(len(list_p_eddy)), list_p_eddy, c="darkblue")
            plt.xlabel("Eddy Currents")
            plt.ylabel("Eddy Currents in Watts")
            plt.legend(["Eddy Currents"])
            plt.title(f"Eddy Currents in Reaction Plate {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_p_eddy, unit="W", fmt=".0f")
            plt.show()


        # 7 V_Rel
        if "v_rel" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_v_rel, total_time)
            plt.scatter(range(len(list_v_rel)), list_v_rel, c="olive")
            plt.xlabel("Relative Velocity")
            plt.ylabel("Relative Velocity Cable/Casing")
            plt.legend(["V_Rel"])
            plt.title(f"V_Rel {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_v_rel, unit="m/s", fmt=".0f")
            plt.show()


        # 8 F_Slip Frequency
        if "f_slip" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_f_slip, total_time)
            plt.scatter(range(len(list_f_slip)), list_f_slip, c="orange")
            plt.xlabel("Slip Frequency")
            plt.ylabel("F Slip Frequency Hz")
            #S    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
            plt.legend(["F Slip Frequency Hz"])
            plt.title(f"F Slip Frequency Hz {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_f_slip, unit="Hz", fmt=".1f")
            plt.show()


        # 9 Slip
        if "slip" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_slip, total_time)
            plt.scatter(range(len(list_slip)), list_slip, c="cyan")
            plt.xlabel("Slip")
            plt.ylabel("Slip %")
            #S    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
            plt.legend(["Slip %"])
            plt.title(f"Slip {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_slip, unit="%", fmt=".1f")
            plt.show()


        # 10 Thrust Power
        if "p_thrust" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_thrust_power, total_time)
            plt.scatter(range(len(list_thrust_power)), list_thrust_power, c="navy")
            plt.xlabel("Thrust Power")
            plt.ylabel("Thrust Power in Watts")
            plt.legend(["Thrust Power"])
            plt.title(f"Thrust Power {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final([p/1e6 for p in list_thrust_power], unit="MW", fmt=".2f")
            plt.show()


        # 11 B Peak
        if "b_peak" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_b_peak, total_time)
            plt.scatter(range(len(list_b_peak)), list_b_peak, c="brown")
            plt.xlabel("Magnetic Field")
            plt.ylabel("Magnetic Field in Tesla")
            plt.legend(["B Peak Reaction Plate"])
            plt.title(f"Magnetic Field at Reaction Plate {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_b_peak, unit="T", fmt=".4f")
            plt.show()


        # 13 Hysteresis
        if "hyst" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_p_hyst, total_time)
            plt.scatter(range(len(list_p_hyst)), list_p_hyst, c="brown")
            plt.xlabel("Hysteresis Losses")
            plt.ylabel("Hysteresis Losses in Watts")
            plt.legend(["Hysteresis Losses"])
            plt.title(f"Hysteresis Losses {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_p_hyst, unit="W", fmt=".0f")
            plt.show()


        # 14 Cryo Power
        if "cryo" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_p_cryo, total_time)
            plt.scatter(range(len(list_p_cryo)), list_p_cryo, c="teal")
            plt.xlabel("Cryogenic Power")
            plt.ylabel("Cryogenic Power in Watts")
            plt.legend(["Cryogenic Power"])
            plt.title(f"Cryogenic Power {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_p_cryo, unit="W", fmt=".0f")
            plt.show()


        # 15 Power Used
        if "power" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_p_lim_site, total_time)
            plt.scatter(range(len(list_p_lim_site)), list_p_lim_site, c="#2F4F4F")
            plt.xlabel("Power Used")
            plt.ylabel("Power Used in Watts")
            plt.legend(["Power Used"])
            plt.title(f"Power Used {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_p_lim_site, unit="W", fmt=".0f")
            plt.show()

        # 15 LIM Power
        if "lim_power" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_p_lim, total_time)
            plt.scatter(range(len(list_p_lim)), list_p_lim, c="#2F4F4F")
            plt.xlabel("LIM Power" )
            plt.ylabel("LIM Power in Watts")
            plt.legend(["LIM Power"])
            plt.title(f"LIM Power {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_p_lim, unit="W", fmt=".0f")
            plt.show()

        # 1 Reaction Plate Temperature in Kelvin
        if "plate_temp" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_temp_plate_ave, total_time)
            plt.scatter(range(len(list_temp_plate_ave)), list_temp_plate_ave, c="lime")
            plt.xlabel("Reaction Plate Temperature")
            plt.ylabel("Reaction Plate Temperature in Kelvin")
            #    plt.yticks([0,1,2,3,4,5,6])
            plt.legend(["Reaction Plate Temperature"])
            plt.title(f"Reaction Plate Temperature in Kelvin {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final(list_temp_plate_ave, unit="K", fmt=".0f")
            plt.show()

        #"""

        # Site Kinetic Energy (TJ)
        if "ke_site" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_E_site_ke, total_time)
            plt.scatter(range(len(list_E_site_ke)), [e/1e12 for e in list_E_site_ke], c="green")
            plt.xlabel("Site Kinetic Energy")
            plt.ylabel("Site Kinetic Energy in TJ")
            plt.legend(["Site KE (thrust only)"])
            plt.title(f"Cumulative Kinetic Energy per Site {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final([e/1e12 for e in list_E_site_ke], unit="TJ", fmt=".1f")
            plt.show()

        # Total Kinetic Energy (EJ)
        if "ke_all" in show or "all" in show:
            tick_pos, tick_labels = make_month_ticks(list_E_total_ke, total_time)
            plt.scatter(range(len(list_E_total_ke)), [e/1e18 for e in list_E_total_ke], c="darkgreen")
            plt.xlabel("Total Kinetic Energy")
            plt.ylabel("Total Kinetic Energy in EJ")
            plt.legend(["Total KE (thrust only)"])
            plt.title(f"Cumulative Kinetic Energy All Sites {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final([e/1e18 for e in list_E_total_ke], unit="EJ", fmt=".2f")
            plt.show()


if __name__ == "__main__":
    main()

