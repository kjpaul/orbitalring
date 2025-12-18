# flake8: noqa
"""
    Calculations for "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import math
import tabulate as tabulate
import matplotlib.pyplot as plt


"""
    Loop timing variables for graphs and data collection
"""
HR = 60 * 60
DAY = HR * 24
WEEK = 7 * DAY
YR = round(DAY * 365.33)
MONTH = 30 * DAY

# graph variables
SAMPLE_PERIOD = 1             # sample step time must me >= DT
SAMPLE_TIME_MAX = 5 * YR      # length of sample (must be longer than deployment)
WRITE_FILE = True
MAKE_GRAPHS = True
SAMPLE_TIME_STR = "-- 100 tonne cable, 12 tonne load" # graph title
SKIP = 40                      # if 1 = no skip
V_REL_LOCK = 5
DT1 = 10  # First day
DT2 = 10    # Data collection
DT3 = 10    #

# -------------------------------------------------
# VARIABLE PARAMETERS, TO BE SET EXTERNALLY
# -------------------------------------------------
N_TURNS = 100                               # turns per phase coil (int)
V_SLIP_MAX = 200                            # 
V_SLIP_MIN = 10
TAU_P = 50.0                                # pole-pitch (m)
W_COIL = 2.0                                # LIM width (m)
GAP = 0.20                                   # coil-to-plate gap (m)
LIM_SPACING = 500                           # distance at which LIMs are place (m)
HTS_D = 80                                # HTS thickness in micrometers (kept for reference)
# NOTE: We are intentionally NOT using a k_fill / Kapton packing model here.
# The prior K_FILL + Kapton thickness logic was mixing insulation/geometry assumptions into inductance/limits
# in a way that was not physically well-founded for this stage of the model.
VOLTS_MAX = 100e3                         # absolute coil voltage limit (RMS), per your assumption (V)
I_C = 800                                   # HTS, field & temperature dependant Ic (A)
I_PEAK = 700                                # I_peak, typically ~75% of Ic (A)
I_TARGET = 650
I_MIN = 10                                  # lower limit on current. Reduce slip instead. (A)
SLIP_MIN = 0.01
P_HEAT_MAX = 100000
T_PLATE = 0.200                             # aluminium thickness (m)
W_TAPE = 0.012                              # HTS tape width is 3mm (m)
ALPHA_ANGLE_DEG = 20                        # magnetic penettration angle of HTS in coils (deg)
HEAT_SINK_L = 30                            # shorter heat sink mean higher ave reaction plate temp.

# Variables that could change, but probably won't
PITCH_COUNT = 3                             # pitches per LIM (int)
CASING_OUTER_W = 10.0                       # external width of orbital ring casing (m)
INV_EFF = 0.90                              # DC to AC inverter efficiency factor (int)
LIM_EFF = 0.95                              # unaccounted for LIM losses efficiency factor (int)

# -------------------------------------------------
# CABLE and LOAD masses are dependant on each other 
# according to Newtown's 3rd Law and the cable tension.
# -------------------------------------------------
M_CABLE_M = 99_198                         # orbital ring cable mass per meter (m)
M_LOAD_M = 12_000                            # orbital ring casing + load mass per meter (m)

# lable data for plots
PARAM_STR1 = f"τp={TAU_P} m, N={N_TURNS}, V_Slip_max={V_SLIP_MAX} m/s, d={LIM_SPACING} m"

# -------------------------------------------------
# BASIC CONSTANTS & USER TUNABLE DESIGN PARAMETERS
# -------------------------------------------------
V_ORBIT = 7754.866                          # 250 km orbital velocity (m/s)
GEO_250_ORBIT = 483.331                     # 250 km ground stationary velocity (m/s)
MU0 = 4 * math.pi * 1e-7                    # permeability of a vacuum (kg m s^-2 A^-2)
SIGMA = 5.670374e-8                         # Stephan-Boltzmann constant (kg s^-3 K^-4)
RHO_ALU_E_106K = 8.35e-9                    # resistivity of Al at 106K (Ωm) (FIX THIS LATER--SHOULD BE CALCULATED)
RHO_ALU_M3 = 2700                           # mass per m^3 of aluminium (kg/m^3)
L_ACTIVE = TAU_P * PITCH_COUNT              # length of LIM (m)
A_LIM = L_ACTIVE * W_COIL                   # area under LIM (m^2)
LIM_PHASES = 3                              # power supply phases into LIMs (int)
L_HTS_COIL = 2 * (W_COIL + TAU_P) * N_TURNS  # HTS length in one coil (m)
L_HTS_LIM = L_HTS_COIL * LIM_PHASES * PITCH_COUNT # HTS in one LIM (m)
L_RING = 41645813.012                       # length of ring in meters (m)
LIM_SITES = round(L_RING / LIM_SPACING)     # total number of LIMs (int)
ALPHA_TAPE = ALPHA_ANGLE_DEG * math.pi / 180 # field penetration angle of HTS (rad)
M_CABLE_T = M_CABLE_M * L_RING              # orbital ring cable mass total (m)
M_LOAD_T = M_LOAD_M * L_RING                # orbital ring casing + load mass total (m)

# -------------------------------------------------
# HEAT RELATED CONSTANTS
# -------------------------------------------------
Q_SUN_1AU = 1361        # solar energy at 1 AU (W/m^2)
Q_EARTH_DAY = 650       # average daytime energy reflected up from earth's surface (W/m^2)
Q_SHEILDING = 0.005     # percentage of heat radiation that passes through multi-layer alu shielding (%/100)
Q_ABS_M = (Q_SUN_1AU + Q_EARTH_DAY) * CASING_OUTER_W * Q_SHEILDING # external heat absorbed by ring (W/m)
Q_ABS_LIM = Q_ABS_M * LIM_SPACING # external heat absorbed by ring (W)
T_FREE_SPACE = 2.7      # temperature of free space (K)
CRYO_EFF = 0.031         # Cryo efficiency: 0.031 -> 32 W are needed for Cryo per W of heat (1/W)
C_P_ALU = 900           # heat capacity or aluminium (J/kg K)
H_CONV_LN2 = 90         # Convective Heat Transfer Coefficient LN2 (W m^-2 K^-1)
T2_LN2_BOIL = 77.4      # boiling point of liquid nitroger at 1 atm (K)
T1_LN2_CRYO = 70        # temperature of LN2 as it leaves cryogenic system (K)
T_N2_HOT = 300          # N2 temperature as it leave the compressor and enters the radiator (K)
T_MAX_Reaction_Plate = 500  # maximum temp for aluminium reaction plate (K)
K_ALU = 205             # thermal conductivity of aluminium (W/m*K)
C_P_LN2 = 2040          # Specific heat capacity of LN2 (J/kg*K)
L_V_LN2 = 199000        # Latent heat of vaporization LN2 (J/kg)
EM_ALU = 0.85           # emisivity or reaction plate (black anodized alu) (num)
EM_HEAT_SINK = 0.9      # emisivity of heat sink surface (num)
PA_LN2_70K = 0.4        # pressure that turns boiling point of LN2 to 70 K (atm)
V_REL_MIN = 10          # fudge factor compensate of missleading initial heat at low v_rel (m/s)
MAX_SITE_POWER = 8.0e6  # power limit per LIM site (W)
MAX_HEATSINK_AREA = (LIM_SPACING) * 2 * W_COIL # the heatsink extends under the coils, since they are 99.99 % open space. 

TEST = True

# Derived helper: mathematically derived constant appearing in equations below
K_B = math.sqrt((24 * RHO_ALU_E_106K * TAU_P**2) / (math.pi**2 * T_PLATE**2))
C_E = (math.pi**2 * T_PLATE**2) / (6 * RHO_ALU_E_106K)
K_F = W_COIL * L_ACTIVE / (2 * MU0)
L_P = math.pi / 2 * TAU_P

#inductance = MU0 * N_TURNS**2 * A_COIL / GAP    # coil inductance
#L_TURN = MU0 * (TAU_P + W_COIL) * (math.log(2 * (TAU_P + W_COIL) / W_TAPE) - 1.5)  # inductance one turn (H) 
#L_COIL = K_FILL * N_TURNS * L_TURN                                            # inductance for coil (H)
A_COIL = TAU_P * W_COIL
L_COIL = MU0 * N_TURNS**2 * A_COIL / W_COIL * (2 / math.pi) * math.atan(W_COIL/(2 * GAP))


PARAM_LIST = {
    "N_TURNS": N_TURNS,
    "V_SLIP_MAX": V_SLIP_MAX,
    "V_SLIP_MIN": V_SLIP_MIN,
    "TAU_P": TAU_P,
    "I_TARGET": I_TARGET,
    "I_MIN": I_MIN,
    "W_COIL": W_COIL,
    "LIM_SPACING": LIM_SPACING,
    "VOLTS_MAX": VOLTS_MAX,
    "I_C": I_C,
    "I_PEAK": I_PEAK,
    "SLIP_MIN": SLIP_MIN,
    "P_HEAT_MAX": P_HEAT_MAX,
    "Q_ABS_LIM": Q_ABS_LIM,
    "GAP": GAP,
    "T_PLATE": T_PLATE,
    "ALPHA_ANGLE_DEG": ALPHA_ANGLE_DEG,
    "PITCH_COUNT": PITCH_COUNT,
    "CASING_OUTER_W": CASING_OUTER_W,
    "M_CABLE_M": M_CABLE_M,
    "M_LOAD_M": M_LOAD_M,
    "CRYO_EFF": CRYO_EFF,
    "EM_ALU": EM_ALU,
    "EM_HEAT_SINK": EM_HEAT_SINK,
    "HEAT_SINK_L": HEAT_SINK_L,
    "MAX_HEATSINK_AREA": MAX_HEATSINK_AREA,
    "V_REL_MIN": V_REL_MIN,
}

# chart lists
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
    """Magnetic field at reaction plate surface for rectangular current sheet."""
    return (2 * MU0 * N_TURNS * i_peak / (math.pi * W_COIL)) * math.atan(W_COIL / (2 * GAP))
#def get_b_plate_peak(i_peak):
#    return MU0 * N_TURNS * i_peak / GAP


def get_b_coil_peak(i_peak):
    #return MU0 * N_TURNS * i_peak / (2 * math.pi * W_COIL / 2)
    return MU0 * N_TURNS * i_peak / (math.pi * (W_COIL + TAU_P))


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


def get_plate_eddy_loss(v_slip, i_peak):
    """LEGACY: retained for reference only (not used by default).

    This thin-plate / sinusoidal-field scaling can greatly over-predict at low f_slip and cryogenic rho
    unless it is solved self-consistently with field loading. Use get_plate_eddy_loss_from_thrust().
    """
    b_plate = get_b_plate_peak(i_peak)
    f_slip = get_slip_frequency(v_slip)
    eff_t_plate = get_eff_plate_depth(f_slip)
    ce = math.pi**2 / (6 * RHO_ALU_E_106K)
    v_plate_eddy = get_plate_eddy_volume(f_slip)
    # return ((math.pi**2 * b_plate**2 * t_plate**2 * f_slip**2) / (6 * RHO_ALU_E_106K)) * v_plate_eddy
    return ce * b_plate**2 * eff_t_plate**2 * f_slip**2 * v_plate_eddy


def get_skin_depth_eddy(f_slip):
    # Skin depth penetration of eddy currents in reaction plate
    if f_slip <= 0:
        f_slip = get_slip_frequency(V_SLIP_MIN)
    return math.sqrt((RHO_ALU_E_106K) / (math.pi * MU0 * f_slip))


def get_eff_plate_depth(f_slip):
    delta_depth = get_skin_depth_eddy(f_slip)
    return min(T_PLATE, delta_depth)


def get_plate_eddy_volume(f_slip):
    depth = get_eff_plate_depth(f_slip)
    return W_COIL * L_ACTIVE * depth


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


def get_slip_efficiency(slip):
    """Dimensionless coupling/efficiency factor vs slip.

    Shape: eta = 2*s / (1 + s^2), peaks at eta=1 when s=1.
    Enforces eta -> 0 as s -> 0 and as s -> inf, and caps to [0,1].
    """
    if slip <= 0:
        return 0.0
    eta = (2.0 * slip) / (1.0 + slip ** 2)
    if eta < 0.0:
        return 0.0
    if eta > 1.0:
        return 1.0
    return eta


def get_f_thrust(f_slip, f_supply, i_peak):
    """Thrust model (dimensionally consistent, physically bounded).

        s = f_slip / f_supply
        F_max = (B^2 / (2*mu0)) * A
        F = eta(s) * F_max

    This deliberately avoids runaway results when secondary back-reaction is not
    solved self-consistently.
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0

    b_plate = get_b_plate_peak(i_peak)

    slip = get_slip_f(f_slip, f_supply)
    if slip < 0:
        slip = 0.0

    return get_slip_efficiency(slip) * get_thrust_force_max(b_plate)
def get_thrust_power(thr, vrel):
    return thr * vrel


def get_coil_volts(i_peak_now, f_supply, b_plate):
    """Induced coil voltage (RMS) from flux linkage (engineering estimate).

    We model the per-phase induced EMF as:
        V_rms = (omega * N * Phi_peak) / sqrt(2)  ≈ 4.44 * N * f * Phi_peak

    with a simple linked flux estimate:
        Phi_peak ≈ B_plate_peak * (W_COIL * TAU_P)

    This is intentionally a *flux-linkage* voltage model (insulation/turn-to-turn stress),
    not the reactive drop I*omega*L.
    """
    if i_peak_now <= 0 or f_supply <= 0:
        return 0.0
    phi_peak = b_plate * W_COIL * TAU_P
    return 4.44 * N_TURNS * f_supply * phi_peak


def set_r_r(f_slip):
    t_plate = get_eff_plate_depth(f_slip)
    return (RHO_ALU_E_106K * L_P) / (t_plate * TAU_P)


"""
    hysteresis
"""
def set_q_hts(i_peak):
    b_coil = get_b_coil_peak(i_peak)
    return b_coil * I_C * W_TAPE * math.sin(ALPHA_TAPE)


def get_p_hyst_lim(f_supply, i_peak):
    q = set_q_hts(i_peak)
    p_coil_hyst = q * L_HTS_COIL * f_supply
    return p_coil_hyst * LIM_PHASES * PITCH_COUNT


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
#    sink_l = min(LIM_SPACING - L_ACTIVE, HEAT_SINK_L)
    sink_l = min(LIM_SPACING, HEAT_SINK_L) # the coils are 99.99% open space, so heatsink extends under them.
    duty = L_ACTIVE / (sink_l + L_ACTIVE)
    p_ave = p_heat * duty
    p_load = p_ave + 0.5 * Q_ABS_LIM # Q_ABS_LIM is the daytime heat from the environment
    return p_load


def get_heatsink_area(p_heat, T_hot=T_N2_HOT, T_cold=T2_LN2_BOIL, em1=EM_ALU, em2=EM_HEAT_SINK):
    # returns area of heatsink (m^2) for a single LIM
    # T_cold = casing, T_hot = reaction plate (K)
    # A is cross-sectional area (m^2)
    # em1 & em2 are emissivity of the resp radiating surfaces where 1 = black body (num)
    # em_eff is effective radiation absorption coefficient between em1 and em2
    #print(p_heat)
    em_eff = get_em_eff(em1, em2)
    if T_hot < T_N2_HOT:
        T_hot = T_N2_HOT
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
    heatsink_length = min(HEAT_SINK_L, LIM_SPACING)    # the coils are 99.99% open space
    time_gap = heatsink_length / v_rel_adjusted
    time_cycle = time_gap + time_lim
    em_eff = get_em_eff()
    # this is get_rad solved for T_hot with the ratio of lim_time/cycle_time added
    return ((p_eddy * time_lim) / (em_eff * SIGMA * (2 * A_LIM) * time_cycle) + T2_LN2_BOIL**4)**(0.25) # (K)
    

def get_delta_t_lim(v_rel, f_slip, p_eddy):
    # returns the increase in reaction plate temp as it passes thru LIM.
    v_rel_adjusted = get_v_rel_min(v_rel) 
    time_lim = L_ACTIVE / v_rel_adjusted
    v_plate_eddy = get_plate_eddy_volume(f_slip)
    m_plate_eddy = v_plate_eddy * RHO_ALU_M3
    return (p_eddy * time_lim) / (m_plate_eddy * C_P_ALU) # (K)



"""
    space elevator
"""
"""
mu      = 3.986e14        # m^3 s^-2
omega   = 7.2921159e-5    # rad s^-1
rho     = 1700            # kg m^-3
sigma   = 7e9             # Pa (14 GPa / 2)

def F(r):                 # antiderivative
    return mu/r + 0.5*omega**2*r**2

k = math.exp((rho/sigma)*(F(r_bottom) - F(r_geo)))   # taper factor
A_geo = A_bottom / k
d_geo = math.sqrt(4*A_geo/math.pi)
"""


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
def get_deployment_time(v_slip, i_peak_now, param_str1):
    # loop variables
    dt = DT1        # initial number of seconds per loop
    i_target = i_peak_now # set initial value of i_min

    print("VOLTS_MAX: ", VOLTS_MAX)

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
    run_once = True
    run_once2 = True
    count = 0
    sample_period = SAMPLE_PERIOD
    sample_time = 0
    skip = 1    # record first round, then start skip
    sample_time_max = SAMPLE_TIME_MAX

    param_str1 = f"{SAMPLE_TIME_STR}\n{param_str1}"
    exit_msg = "PASSED"
    make_graphs = MAKE_GRAPHS

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
        for _ctrl in range(6):
            # Relative velocity between cable and casing (your v_rel definition)
            v_rel = get_v_rel(vcable, vcasing)

            # Travelling wave speed and slip kinematics
            v_wave = v_rel + v_slip
            slip = get_slip(v_slip, v_rel)

            # Frequencies
            f_slip = get_slip_frequency(v_slip)
            f_supply = get_supply_frequency(v_wave)

            # Magnetic field at plate
            b_plate = get_b_plate_peak(i_peak_now)

            # Thrust and thrust power (per LIM)
            thrust = get_f_thrust(f_slip, f_supply, i_peak_now)
            p_thrust = get_thrust_power(thrust, v_rel)

            # Induced coil voltage (RMS) from flux linkage
            volts_lim = get_coil_volts(i_peak_now, f_supply, b_plate)

            # Secondary and hysteresis losses (per LIM)
            p_eddy = get_plate_eddy_loss_from_thrust(v_slip, thrust)

            # Hard clamp: secondary power cannot exceed magnetic energy flux for same B
            p_eddy_cap = get_thrust_power_max(b_plate, v_slip)
            if p_eddy > p_eddy_cap:
                p_eddy = p_eddy_cap

            p_hyst = get_p_hyst_lim(f_supply, i_peak_now)

            # Heating (per LIM) used for thermal and cryo sizing
            p_heat = p_eddy + p_hyst

            # Reaction plate temperature estimate
            temp_plate_ave = get_plate_temp(v_rel, p_eddy)
            delta_temp_lim = get_delta_t_lim(v_rel, f_slip, p_eddy)
            alu_temp_out = temp_plate_ave + delta_temp_lim / 2

            # Heat dissipation surfaces (per LIM)
            p_heat_load = get_p_load(p_eddy)
            min_heatsink_area = get_heatsink_area(p_heat_load, alu_temp_out)
            min_T_ambient = get_T_min_ambient((2 * MAX_HEATSINK_AREA), p_heat_load)  # 2 heatsinks per LIM site

            # Total LIM-side power bookkeeping (per LIM, then per site)
            p_total_lim = (p_thrust + p_heat_load + p_hyst) * (2 - LIM_EFF)

            # Cryo power is per LIM site (2 LIMs); include external absorbed heat via p_heat_load already.
            if v_rel != 0:
                heat_t = 2 * p_heat_load
                p_cryo = get_p_cryo(heat_t, alu_temp_out)
            else:
                p_cryo = 0.0

            lim_site_power = (2 * p_total_lim + p_cryo) * (2 - INV_EFF)

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
            if (temp_plate_ave > T_MAX_Reaction_Plate) or ((temp_plate_ave + delta_temp_lim) > T_MAX_Reaction_Plate):
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

            if not changed:
                break

        # --- Apply this timestep's acceleration using the LIMIT-COMPLIANT thrust ---
        vcable = get_v_cable(vcable, thrust, dt)
        vcasing = get_v_casing(vcasing, thrust, dt)

        # Accumulate kinetic energy (thrust power only, no losses)
        # Each site has 2 LIMs
        E_site_ke += p_thrust * 2.0 * dt
        E_total_ke += p_thrust * LIM_SITES * 2.0 * dt

        # Startup ramp (only if comfortably below limits)
        if time < HR:  # First hour
            if (temp_plate_ave < T_MAX_Reaction_Plate * 0.8 and volts_lim < VOLTS_MAX * 0.8 and lim_site_power < MAX_SITE_POWER * 0.8):
                i_peak_now += (I_TARGET - i_peak_now) * 0.01
                v_slip += (V_SLIP_MAX - v_slip) * 0.01
        # Starting values for table, ignore 0. Only run once.
        if time == 1: 
            power_track.append(["START", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])

        if time > sample_time and time < sample_time_max:
            sample_time += sample_period
            skip -= 1
            skin_depth_eff = get_eff_plate_depth(f_slip)
            assert skin_depth_eff > 0.0
            skin_depth_calc = get_skin_depth_eddy(f_slip)
            assert skin_depth_calc > 0.0
            # collect data from graphs
            if skip == 0:
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


                count += 1
                skip = SKIP

        # collect data at various early intervals
        if time > second and time < 60:
            power_track.append([f"{round(time)} sec", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            second += 5
        if time > minute and time < HR:
            power_track.append([f"{round(time/60)} min", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            if minute == 60:
                minute += 60*4
            else:
                minute += 60*5
        if time > hours * HR and time < 24 * HR:
            power_track.append([f"{round(time/HR)} hrs", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            hours += 1
        if time > days * DAY and time < 31 * DAY:
            power_track.append([f"{round(time/DAY)} day", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            days += 1

        # show progress to stdout and collect data monthly
        if time > months * MONTH:
            power_track.append([f"{round(time/MONTH)} mth", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            months += 1
            # show progress to screen once a month
            #sys.stdout.write(".")
            sys.stdout.flush()

        # show yearly progress to stdout
        if time > mts * MONTH:
            mthstr = str(mts)+" "
            sys.stdout.write(mthstr)
            prog = get_progress(vcasing)
            if mts == 0:
                print(f"dt =  {dt}")
                print(f"Months | Prog |   Volts   | I_Peak  |  V_Slip  |  Slip  |F_Supply|  Eddy  |  Hyst  |LIM Power  | Site Power| Thrust")
                print(f"{mthstr} mts, {prog:.2f}%, {volts_lim:.2f} V, {i_peak_now:.2f} A, {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:.2f} Hz, {p_eddy:.2f} W, {p_hyst:.2f} W, {p_total_lim/1e6:.3f} MW, {lim_site_power/1e6:.3f} MW, {thrust:.3f} N")
            else:
                print(f"mts, {prog:.2f}%, {volts_lim:.2f} V, {i_peak_now:.2f} A, {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:.2f} Hz, {p_eddy:.2f} W, {p_hyst:.2f} W, {p_total_lim/1e6:.3f} MW, {lim_site_power/1e6:.3f} MW, {thrust:.3f} N")
            mts += 1
 
        # collect min-max data
        if v_rel != 0:
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
                max_min["min_heatsink_area"] = [round(min_heatsink_area), time, "m^2/LIM", min_heatsink_area]

            if max_min["min_ambient_T"][3] < min_T_ambient:
                max_min["min_ambient_T"] = [round(min_T_ambient), time, "K", min_T_ambient]

            if max_min["p_heat_load"][3] < p_heat_load:
                max_min["p_heat_load"] = [round(p_heat_load), time, "W", p_heat_load]

            if max_min["alu_temp_out"][3] < alu_temp_out:
                max_min["alu_temp_out"] = [round(alu_temp_out), time, "W", alu_temp_out]

            if max_min["p_cryo_max"][3] < p_cryo:
                max_min["p_cryo_max"] = [round(p_cryo,2), time, "W", p_cryo]

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

    print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(f"{exit_msg}\t\t{time/YR:.2f} yrs")
    # post loop data collection
    power_track.append([f"{round(time/YR,2)} yrs", round(vcasing,1), round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
    power_track.append(["Time", "V_Shell", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "LIM Power", "Site Power"])

    str2 = f"Deployment time: {time / DAY:.2f} days, {time / YR:.2f} years"
    str3 = f"Cable velocity: {vcable:.2f} m/s"
    str4 = f"Casing velocity: {vcasing:.2f} m/s"

    # post loop wrap-up. show data.
    lines = ["\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", exit_msg+"\n", str2+"\n", str3+"\n", str4+"\n"]
    lines.append("\n--------------------------------------------------------------------\n")
    print("--------------------------------------------------------------------")
    for key, value in PARAM_LIST.items():
        str1 = (f"\t{key:15}\t{value:8}")
        print(str1)
        lines.append(str1+"\n")
    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")
    for key, value in max_min.items():
        if value[1] < 60:
            t = f"{value[1]:7.2f} seconds"
        elif value[1] < HR:
            t = f"{value[1]/60:7.2f} minutes"
        elif value[1] < DAY:
            t = f"{value[1]/HR:7.2f} hours"
        elif value[1] < MONTH:
            t = f"{value[1]/DAY:7.2f} days"
        elif value[1] < YR:
            t = f"{value[1]/MONTH:7.2f} months"
        else:
            t = f"{value[1]/YR:7.2f} years"
        str1 = (f"\t{key:15} {value[0]:18.2f} {value[2]:4} {t}")
        print(str1)
        lines.append(str1+"\n")
    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")
   
    print("")
    print(tabulate.tabulate(power_track, headers=["Time", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "LIM Power", "Site Power"]))
    lines.append(tabulate.tabulate(power_track, headers=["Time", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "LIM Power", "Site Power"]))
    print("--------------------------------------------------------------------")
    lines.append("\n--------------------------------------------------------------------\n")

    # print out max_min data
 
    print("")

    print(str2)
    print(str3)
    print(str4)
    print("Count = ", count)
    print(f"Site KE (2 LIMs): {E_site_ke/1e12:.2f} TJ")
    print(f"Total KE (all sites): {E_total_ke/1e18:.4f} EJ")

    if exit_msg != "PASSED":
        print(f"{exit_msg}. time in seconds: {time:.3f} s, years: {time/YR:.2f}")

    param_str1 = param_str1 + f", F_Supply: {round(f_supply,1)} Hz, Deployment time: {round(time / YR, 2)} years."

    tick_count = round(time / YR)

    if WRITE_FILE:
        with open("./output/_orbital_ring_options_03.txt", "a") as file:
            file.writelines(lines)

    return [tick_count, f"{param_str1} ", make_graphs, time]

def main() -> None:
    show = []
    if len(sys.argv) > 1:
        show = sys.argv
    v_slip = V_SLIP_MAX
    i_peaks = I_MIN
    param = get_deployment_time(v_slip, i_peaks, PARAM_STR1)



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
            annotate_final(list_skin_depth_calc, unit="mm", fmt=".1f")
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

