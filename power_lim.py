# flake8: noqa
"""
    Calculations for "How to Build an Orbital Ring" by Paul de Jong
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
SAMPLE_TIME_MAX = YR          # length of sample
TICKS = 6                      # ticks on x-axis of graphs
WRITE_FILE = True
MAKE_GRAPHS = False
SAMPLE_TIME_STR = "Full Deployment in 6 Months" # graph title
SKIP = 40                      # if 1 = no skip
V_REL_LOCK = 5
DT1 = 1  # First day
DT2 = 1    # Data collection
DT3 = 1    #

# -------------------------------------------------
# VARIABLE PARAMETERS, TO BE SET EXTERNALLY
# -------------------------------------------------
N_TURNS = 8                                 # turns per phase coil (int)
V_SLIP_MAX = 460                            # 
V_SLIP_MIN = 10
TAU_P = 50.0                                # pole-pitch (m)
W_COIL = 1.0                                # LIM width (m)
LIM_SPACING = 500                           # distance at which LIMs are place (m)
D_KAPTON = 0.01                             # Kapton tape thickness in mm (mm)
K_FILL = 0.002 / (0.002 + D_KAPTON)         # 0.002 is HTS thickness (int)
KAPTON_SAFE_V = 1e5 * D_KAPTON              # Kapton breakdown 2e5 V/mm safe = 1e5 V/mm (V)
KAPTON_V = (N_TURNS - 1) * KAPTON_SAFE_V    # Max voltage based on D_KAPTON (V)
VOLTS_MAX = min(2000, KAPTON_V)            # set a reasonable limit on voltage (V)
THRUST_TARGET = 1000                        # 
I_C = 800                                   # HTS, field & temperature dependant Ic (A)
I_PEAK = 700                                # I_peak, typically ~75% of Ic (A)
I_TARGET = 650
I_MIN = 10                                  # lower limit on current. Reduce slip instead. (A)
SLIP_MIN = 0.01
P_HEAT_MAX = 100000
GAP = 0.2                                   # coil-to-plate gap (m)
T_PLATE = 0.040                             # aluminium thickness (m)
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
M_CABLE_M = 26300                         # orbital ring cable mass per meter (m)
M_LOAD_M = 5000                            # orbital ring casing + load mass per meter (m)

# lable data for plots
PARAM_SRT1 = f"Pitch={TAU_P} m, N={N_TURNS}, Tape={W_TAPE*1000} mm, Spacing={LIM_SPACING} m"
PARAM_STR2 = f"Ic={I_C} A, I_peak={I_PEAK} A, V_max={VOLTS_MAX} V, V_Slip_max={V_SLIP_MAX} m/s"

# -------------------------------------------------
# BASIC CONSTANTS & USER TUNABLE DESIGN PARAMETERS
# -------------------------------------------------
V_ORBIT = 7754.866                          # 250 km orbital velocity (m/s)
GEO_250_ORBIT = 483.331                     # 250 km ground stationary velocity (m/s)
MU0 = 4 * math.pi * 1e-7                    # permeability of a vacuum (kg m s^-2 A^-2)
SIGMA = 5.670374e-8                         # Stephan-Boltzmann constant (kg s^-3 K^-4)
RHO_ALU_E = 2.86e-8                         # resistivity of Al at 300K (Ωm)
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
Q_SHEILDING = 0.005      # percentage of heat radiation that passes through multi-layer alu shielding (%/100)
Q_ABS_M = (Q_SUN_1AU + Q_EARTH_DAY) * CASING_OUTER_W * Q_SHEILDING # external heat absorbed by ring (W/m)
Q_ABS_LIM = Q_ABS_M * LIM_SPACING # external heat absorbed by ring (W)
T_FREE_SPACE = 2.7      # temperature of free space (K)
CRYO_EFF = 0.05          # Cryo efficiency: 0.05 -> 20 W are needed for Cryo per W of heat (1/W)
C_P_ALU = 900           # heat capacity or aluminium (J/kg K)
H_CONV_LN2 = 90         # Convective Heat Transfer Coefficient LN2 (W m^-2 K^-1)
T2_LN2_BOIL = 77.4      # boiling point of liquid nitroger at 1 atm (K)
T1_LN2_CRYO = 70        # temperature of LN2 as it leaves cryogenic system (K)
T_N2_HOT = 200          # N2 temperature as it leave the compressor and enters the radiator (K)
T_MAX_Reaction_Plate = 500      # maximum temp for aluminium reaction plate (K)
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
K_B = math.sqrt((24 * RHO_ALU_E * TAU_P**2) / (math.pi**2 * T_PLATE**2))
C_E = (math.pi**2 * T_PLATE**2) / (6 * RHO_ALU_E)
C_H = (MU0 / (2 * math.pi * I_C**2)) * K_FILL
K_F = W_COIL * L_ACTIVE / (2 * MU0)
L_P = math.pi / 2 * TAU_P

#inductance = MU0 * N_TURNS**2 * A_COIL / GAP    # coil inductance
#L_TURN = MU0 * (TAU_P + W_COIL) * (math.log(2 * (TAU_P + W_COIL) / W_TAPE) - 1.5)  # inductance one turn (H) 
#L_COIL = K_FILL * N_TURNS * L_TURN                                            # inductance for coil (H)
A_COIL = TAU_P * W_COIL
L_COIL = N_TURNS**2 * MU0 * (A_COIL / GAP) * K_FILL


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
    "D_KAPTON": D_KAPTON,
    "K_FILL": round(K_FILL,3),
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
    return MU0 * N_TURNS * i_peak / GAP


def get_b_coil_peak(i_peak):
    #return MU0 * N_TURNS * i_peak / (2 * math.pi * W_COIL / 2)
    return MU0 * N_TURNS * i_peak / (math.pi * (W_COIL + TAU_P))


"""
    eddy currents & skin depth
"""
def get_plate_eddy_loss(v_slip, i_peak):
    # Eddy current losses in the reaction plate of one LIM
    b_plate = get_b_plate_peak(i_peak)
    f_slip = get_slip_frequency(v_slip)
    eff_t_plate = get_eff_plate_depth(f_slip)
    ce = math.pi**2 / (6 * RHO_ALU_E)
    v_plate_eddy = get_plate_eddy_volume(f_slip)
    # return ((math.pi**2 * b_plate**2 * t_plate**2 * f_slip**2) / (6 * RHO_ALU_E)) * v_plate_eddy
    return ce * b_plate**2 * eff_t_plate**2 * f_slip**2 * v_plate_eddy


def get_skin_depth_eddy(f_slip):
    # Skin depth penetration of eddy currents in reaction plate
    if f_slip <= 0:
        f_slip = get_slip_frequency(V_SLIP_MIN)
    return math.sqrt((2 * RHO_ALU_E) / (2 * math.pi * MU0 * f_slip))


def get_eff_plate_depth(f_slip):
    delta_depth = get_skin_depth_eddy(f_slip)
    return min(T_PLATE, delta_depth)


def get_plate_eddy_volume(f_slip):
    depth = get_eff_plate_depth(f_slip)
    return W_COIL * L_ACTIVE * depth


"""
    thrust and power
"""
def get_f_thrust(f_slip, f_supply, i_peak):
    b_plate = get_b_plate_peak(i_peak)
    rr = set_r_r(f_slip)
    slip = get_slip_f(f_slip, f_supply)
    v_plate_eddy = get_plate_eddy_volume(f_slip)
    # return ((b_plate**2 * slip) / rr) * v_plate_eddy * f_slip
    return ((b_plate**2 * slip/(1+slip**2)) / rr) * v_plate_eddy


def get_thrust_power(thr, vrel):
    return thr * vrel


def get_coil_volts(i_peak_now, f_supply, b_plate):
    if i_peak_now == 0:
        return 0
    return 2 * math.pi * f_supply * L_COIL * i_peak_now / math.sqrt(3)
    #phi_max = b_plate * W_COIL * TAU_P
    #return 4.44 * N_TURNS * f_supply * phi_max # 2*pi/sqrt(2) ~ 4.44


def set_r_r(f_slip):
    t_plate = get_eff_plate_depth(f_slip)
    return (RHO_ALU_E * L_P) / (t_plate * TAU_P)


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


def get_cop_actual(T_cold, T_hot=T_N2_HOT, eff=0.05):
    # effective coefficient of performance (COP) of the cryogenic system (num)
    # eff is the overall efficiency of the cryogenic system. 0.05 = 1/20 ratio of performance (num)
    if T_hot == T_cold:
        T_hot = T_cold + 1
    return (T_cold/(T_hot - T_cold))*eff


def get_p_cryo(p_heat, T_hot=T_N2_HOT, T_cold=T2_LN2_BOIL, eff=0.05):
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
    p_load = p_ave + 0.5 * Q_ABS_LIM
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
    return round((vcase - GEO_250_ORBIT) / (V_ORBIT - GEO_250_ORBIT) * 100, 2)


"""
    # -------------------------------------------------
    # MAIN LOOP. CALCULATE DEPLOYMENT TIME
    # -------------------------------------------------
"""
def get_deployment_time(v_slip, i_peak_now, param_str1, param_str2):
    # loop variables
    dt = DT1        # initial number of seconds per loop
    i_target = i_peak_now # set initial value of i_min

    print("VOLTS_MAX: ", VOLTS_MAX)
    print("K_FILL: ", K_FILL)

    # variable with starting conditions
    vcable = V_ORBIT
    vcasing = V_ORBIT
    time = 0
    yrs = 0
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
    run_once = True
    run_once2 = True
    count = 0
    sample_period = SAMPLE_PERIOD
    sample_time = 0
    skip = 1    # record first round, then start skip
    sample_time_max = SAMPLE_TIME_MAX

    param_str1 = f"{SAMPLE_TIME_STR}\n{param_str1}"
    param_str2 = param_str2 + f", dt1={DT1} sec, dt2={DT2} sec"
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

        # electric phase velocity & slip
        v_rel = get_v_rel(vcable, vcasing)

        v_wave = v_rel + v_slip
        slip = get_slip(v_slip, v_rel)

        # frequency
        f_slip = get_slip_frequency(v_slip)
        f_supply = get_supply_frequency(v_wave)

        # magnetic field
        #b_coil = get_b_coil_peak(i_peak_now)
        b_plate = get_b_plate_peak(i_peak_now)

        # eddy current penetration depth, volume and mass
        skin_depth_eff = get_eff_plate_depth(f_slip)
        skin_depth_calc = get_skin_depth_eddy(f_slip)

        # thrust one LIM
        thrust = get_f_thrust(f_slip, f_supply, i_peak_now)
        p_thrust = get_thrust_power(thrust, v_rel)

        # induced coil voltage is a critical limiting factor
        #volts_lim = get_coil_volts(i_peak_now, p_thrust)
        volts_lim = get_coil_volts(i_peak_now, f_supply, b_plate)

        # get the new cable and casing velocities based on thrust * delta time
        vcable = get_v_cable(vcable, thrust, dt)        # uses total mass and 
        vcasing = get_v_casing(vcasing, thrust, dt)     # thrust for all LIMs

        # eddy current and hysteresis losses for one LIM
        p_eddy = get_plate_eddy_loss(v_slip, i_peak_now)
        p_hyst = get_p_hyst_lim(f_supply, i_peak_now)
        p_heat = p_eddy + p_hyst

        # equilibrium temperature of reaction plate
        temp_plate_ave = get_plate_temp(v_rel, p_eddy)
        # delta tempurture of reaction plate as it exits LIM
        delta_temp_lim = get_delta_t_lim(v_rel, f_slip, p_eddy) 
        alu_temp_out = temp_plate_ave + delta_temp_lim / 2
        if temp_plate_ave > T_MAX_Reaction_Plate or (temp_plate_ave + delta_temp_lim) > T_MAX_Reaction_Plate:
            i_peak_now *= 0.95  # Reduce current by 5%
            v_slip *= 0.95
        if time < HR:  # First hour
            if temp_plate_ave < T_MAX_Reaction_Plate * 0.8 and volts_lim < VOLTS_MAX * 0.8:
                i_peak_now += (I_TARGET - i_peak_now) * 0.01  # Gradual increase
                v_slip += (V_SLIP_MAX - v_slip) * 0.01

        # calculate heat dissipation surfaces
        p_heat_load = get_p_load(p_eddy)
        min_heatsink_area = get_heatsink_area(p_heat_load, alu_temp_out)
        min_T_ambient = get_T_min_ambient((2 * MAX_HEATSINK_AREA), p_heat_load) # 2 heatsinks per LIM site


        # sum up all the losses for one LIM (cryo to be added)
        p_total_lim = (p_thrust + p_heat_load + p_hyst) * (2 - LIM_EFF)

        if v_rel != 0:
            heat_t = 2 * p_heat_load         # heat from 2 LIMs + what is absobed from outside
            p_cryo = get_p_cryo(heat_t, alu_temp_out)

        lim_site_power = (2 * p_total_lim + p_cryo) * (2 - INV_EFF)
        #if thrust < THRUST_TARGET:
        if lim_site_power < MAX_SITE_POWER * 0.9:
            if i_peak_now < I_TARGET:
                i_peak_now += (I_TARGET - i_peak_now) * 0.001
            if v_slip < V_SLIP_MAX:
                v_slip += (V_SLIP_MAX - v_slip) * 0.001
                f_slip = get_slip_frequency(v_slip)
        elif lim_site_power > MAX_SITE_POWER:
            if i_peak_now > I_MIN:
                i_peak_now *= 0.95
            if v_slip > V_SLIP_MIN:
                v_slip *= 0.95
                f_slip = get_slip_frequency(v_slip)

        # Starting values for table, ignore 0. Only run once.
        if time == 1: 
            power_track.append(["START", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])

        if time > sample_time and time < sample_time_max:
            sample_time += sample_period
            skip -= 1
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


                count += 1
                skip = SKIP

        # collect data at various early intervals
        if time > second and time < 60:
            power_track.append([f"{round(time)} sec", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            second += 5
        if time > minute and time < HR:
            power_track.append([f"{round(time/60)} min", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            if minute == 60:
                minute += 60*4
            else:
                minute += 60*5
        if time > hours * HR and time < 24 * HR:
            power_track.append([f"{round(time/HR)} hrs", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            hours += 1
        if time > days * DAY and time < 31 * DAY:
            power_track.append([f"{round(time/DAY)} day", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            days += 1

        # show progress to stdout and collect data monthly
        if time > months * MONTH:
            power_track.append([f"{round(time/MONTH)} mth", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
            months += 1
            # show progress to screen once a month
            sys.stdout.write(".")
            sys.stdout.flush()

        # show yearly progress to stdout
        if time > yrs * YR:
            years = str(yrs)+" "
            sys.stdout.write(years)
            prog = get_progress(vcasing)
            if yrs == 0:
                print(f"dt =  {dt}")
                print(f"............Years  | Prog  | Volts | I_Peak |  V_Slip  |  Slip  |F_Supply|  Eddy | Hyst  |LIM Power  | Site Power| Thrust")
                print(f"............{years} yrs, {prog:.2f}%, {volts_lim:.2f} V, {i_peak_now:.2f} A, {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:.2f} Hz, {p_eddy:.2f} W, {p_hyst:.2f} W, {p_total_lim/1e6:.3f} MW, {lim_site_power/1e6:.3f} MW, {thrust:.3f} N")
            else:
                print(f"yrs, {prog:.2f}%, {volts_lim:.2f} V, {i_peak_now:.2f} A, {v_slip:.2f} m/s, {slip*100:.2f}%, {f_supply:.2f} Hz, {p_eddy:.2f} W, {p_hyst:.2f} W, {p_total_lim/1e6:.3f} MW, {lim_site_power/1e6:.3f} MW, {thrust:.3f} N")
            yrs += 1
 
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
    power_track.append([f"{round(time/YR,2)} yrs", round(i_peak_now,1), round(volts_lim,1), round(v_slip,1), f"{round(slip*100,1)}%", round(f_supply,1), round(p_eddy,1), round(p_hyst,1), round(thrust,1), f"{round(p_thrust/1e6,3)} MW", f"{round(p_total_lim/1e6,3)} MW", f"{round(lim_site_power/1e6,3)} MW"])
    power_track.append(["Time", "Amps","Volts", "V_Slip", "Slip", "F_Supply", "Eddy", "Hyst", "Thrust", "P_Thrust", "LIM Power", "Site Power"])

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
    if exit_msg != "PASSED":
        print(f"{exit_msg}. time in seconds: {time:.3f} s, years: {time/YR:.2f}")

    param_str1 = param_str1 + f", F_Supply: {round(f_supply,1)} Hz, Deployment time: {round(time / YR, 2)} years."

    tick_count = round(time / YR)

    if WRITE_FILE:
        with open("_orbital_ring_options_03.txt", "a") as file:
            file.writelines(lines)

    return [tick_count, f"{param_str1}\n{param_str2}, steps: {count}", make_graphs]

def main() -> None:
    v_slip = 50
    i_peaks = 50
    param = get_deployment_time(v_slip, i_peaks, PARAM_SRT1, PARAM_STR2)

    if param[2]:
        param_str = param[1]

        x_lable = SAMPLE_TIME_STR
        #ticks = param[0]
        ticks = TICKS

        print("ticks: ",ticks)
        
        # 1 Current in Amps
        list_len = len(list_i_peak)
        list_len_tick = list_len/ticks
        plt.scatter(range(list_len), list_i_peak, c="blue")
        plt.xlabel("Current: " + x_lable)
        plt.ylabel("Current in Amps")
        plt.legend(["Amps"])
        plt.title(f"Current {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 2 Volts
        list_len = len(list_volts)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_volts)), list_volts, c="red")
        plt.xlabel("Voltage: " + x_lable)
        plt.ylabel("Voltage")
        plt.legend(["Volts"])
        plt.title(f"Volts {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 3 V_Slip Velocity
        list_len = len(list_v_slip)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_v_slip)), list_v_slip, c="green")
        plt.xlabel("Slip Velocity: " + x_lable)
        plt.ylabel("Slip Velocity in m/s")
        plt.legend(["Slip Velocity"])
        plt.title(f"V_Slip {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()
        

        # 4 Skin Calc & Eff
        list_len = len(list_skin_depth_eff)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_skin_depth_calc)), list_skin_depth_calc, c="darkgreen")
        plt.scatter(range(len(list_skin_depth_eff)), list_skin_depth_eff, c="magenta")
        plt.xlabel("Skin Depth (Calc & Eff): " + x_lable)
        plt.ylabel("Skin Depth/Plate Thickness in mm")
    #    plt.yticks([0,1,2,3,4,5,6])
        plt.legend(["Skin Depth Calc", "Skin Depth Eff"])
        plt.title(f"Effective Reaction Plate Plate Thickness in mm {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 5 Thrust
        list_len = len(list_thrust)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_thrust)), list_thrust, c="purple")
        plt.xlabel("Thrust: " + x_lable)
        plt.ylabel("Thrust in Newtons")
        plt.legend(["Thrust"])
        plt.title(f"Thrust {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        #"""
        # 6 P_Eddy
        list_len = len(list_p_eddy)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_p_eddy)), list_p_eddy, c="darkblue")
        plt.xlabel("Eddy Currents: " + x_lable)
        plt.ylabel("Eddy Currents in Watts")
        plt.legend(["Eddy Currents"])
        plt.title(f"Eddy Currents in Reaction Plate {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        #"""

        # 7 V_Rel
        list_len = len(list_v_rel)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_v_rel)), list_v_rel, c="olive")
        plt.xlabel("Relative Velocity: " + x_lable)
        plt.ylabel("Relative Velocity Cable/Casing")
        plt.legend(["V_Rel"])
        plt.title(f"V_Rel {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 8 F_Slip Frequency
        list_len = len(list_f_slip)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_f_slip)), list_f_slip, c="orange")
        plt.xlabel("Slip Frequency: " + x_lable)
        plt.ylabel("F Slip Frequency Hz")
    #S    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
        plt.legend(["F Slip Frequency Hz"])
        plt.title(f"F Slip Frequency Hz {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 9 Slip
        list_len = len(list_slip)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_slip)), list_slip, c="cyan")
        plt.xlabel("Slip: " + x_lable)
        plt.ylabel("Slip %")
    #S    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
        plt.legend(["Slip %"])
        plt.title(f"Slip {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 10 Thrust Power
        list_len = len(list_thrust_power)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_thrust_power)), list_thrust_power, c="navy")
        plt.xlabel("Thrust Power: " + x_lable)
        plt.ylabel("Thrust Power in Watts")
        plt.legend(["Thrust Power"])
        plt.title(f"Thrust Power {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 11 B
        list_len = len(list_b_peak)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_b_peak)), list_b_peak, c="brown")
        plt.xlabel("Magnetic Field: " + x_lable)
        plt.ylabel("Magnetic Field in Tesla")
        plt.legend(["B Peak Reaction Plate"])
        plt.title(f"Magnetic Field at Reaction Plate {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 13 Hysteresis
        list_len = len(list_p_hyst)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_p_hyst)), list_p_hyst, c="brown")
        plt.xlabel("Hysteresis Losses: " + x_lable)
        plt.ylabel("Hysteresis Losses in Watts")
        plt.legend(["Hysteresis Losses"])
        plt.title(f"Hysteresis Losses {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 14 Cryo Power
        list_len = len(list_p_cryo)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_p_cryo)), list_p_cryo, c="teal")
        plt.xlabel("Cryogenic Power: " + x_lable)
        plt.ylabel("Cryogenic Power in Watts")
        plt.legend(["Cryogenic Power"])
        plt.title(f"Cryogenic Power {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        # 15 Hysteresis
        list_len = len(list_p_lim_site)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_p_lim_site)), list_p_lim_site, c="#2F4F4F")
        plt.xlabel("Power Used: " + x_lable)
        plt.ylabel("Power Used in Watts")
        plt.legend(["Power Used"])
        plt.title(f"Power Used {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()

        # 15 LIM Power
        list_len = len(list_p_lim)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_p_lim)), list_p_lim, c="#2F4F4F")
        plt.xlabel("LIM Power: " + x_lable)
        plt.ylabel("LIM Power in Watts")
        plt.legend(["LIM Power"])
        plt.title(f"LIM Power {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()

        # 1 Reaction Plate Temperature in Kelvin
        list_len = len(list_temp_plate_ave)
        list_len_tick = list_len/ticks
        plt.scatter(range(len(list_temp_plate_ave)), list_temp_plate_ave, c="lime")
        plt.xlabel("Reaction Plate Temperature: " + x_lable)
        plt.ylabel("Reaction Plate Temperature in Kelvin")
    #    plt.yticks([0,1,2,3,4,5,6])
        plt.legend(["Reaction Plate Temperature"])
        plt.title(f"Reaction Plate Temperature in Kelvin {param_str}")
        plt.xticks([
                0,list_len_tick,list_len_tick*2,list_len_tick*3,
                list_len_tick*4,list_len_tick*5,list_len_tick*6,
                #list_len_tick*7,list_len_tick*8, list_len_tick*9,
                #list_len_tick*10,list_len_tick*11,list_len_tick*12,
                #list_len_tick*13,list_len_tick*14,list_len_tick*15,
                #list_len_tick*16,list_len_tick*17,list_len_tick*18,
                #list_len_tick*19,list_len_tick*20,list_len_tick*21,
                #list_len_tick*22,list_len_tick*23,list_len_tick*24,
                #list_len_tick*25,list_len_tick*26,list_len_tick*27,
                #list_len_tick*28,list_len_tick*29,list_len_tick*30,
            ],[
                "0", "1", "2", "3", "4", "5", "6",
                #"7", "8", "9", "10","11","12",
                #"13", "14", "15", "16", "17",
                #"18", "19", "20", "21", "22", "23", "24",
                #"25", "26", "27", "28", "29", "30",
            ])
        plt.show()


        #"""

if __name__ == "__main__":
    main()

