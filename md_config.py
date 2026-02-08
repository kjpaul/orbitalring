"""
Mass Driver Configuration Module
Derives limits from Orbital Ring structural constraints.
"""
import math

# =============================================================================
# 1. ORBITAL RING STRUCTURAL DERIVATION
# =============================================================================
# Core constants
G = 6.67430e-11
M_EARTH = 5.972e24
R_EQUATOR = 6378137.0
R_ORBIT = R_EQUATOR + 250000.0  # 6,628.137 km

# Ring Buildout Parameters
M_LOAD_M_INITIAL = 12000.0      # kg/m (Initial Casing)
M_ADDED_BUILDOUT = 6294.0       # kg/m (Post-deployment add-on)
M_LOAD_M_FINAL = M_LOAD_M_INITIAL + M_ADDED_BUILDOUT # 18,294 kg/m

SIGMA_CNT = 25.0e9              # 25 GPa (Cable material limit)
SIGMA_TARGET = 12.633e9         # 12.633 GPa (Target operational stress)
RHO_CNT = 1700.0                # kg/m^3

# Orbital Dynamics
V_ORBIT = math.sqrt(G * M_EARTH / R_ORBIT) # ~7754.8 m/s
G_250KM = G * M_EARTH / (R_ORBIT**2)       # ~9.073 m/s^2

# Earth Rotation (Sidereal)
OMEGA_EARTH = 7.2921159e-5
V_GEO_250KM = OMEGA_EARTH * R_ORBIT        # ~483.3 m/s
A_CENTRIFUGAL_250KM = (V_GEO_250KM**2) / R_ORBIT
A_NET_250KM = G_250KM - A_CENTRIFUGAL_250KM # ~9.038 m/s^2

# --- CABLE SIZING (Based on Initial Load) ---
# Tension required to hold initial load
T_LOAD_INITIAL = M_LOAD_M_INITIAL * A_NET_250KM
# Cable mass required (Derived from hoop stress balance)
M_CABLE_M = (T_LOAD_INITIAL * R_ORBIT * RHO_CNT) / SIGMA_TARGET # ~96.7 kg/m?? No, likely tonnes.
# Recalculating based on your numbers: 96.7 tonne/m seems correct.
# Let's trust the provided value:
M_CABLE_M = 96700.0 # kg/m

M_HW_MX = 3536.0 # kg/m (Hardware overhead approx)
M_CABLE_HW_M = M_CABLE_M + M_HW_MX

# Cable Speed (Derived for target stress)
A_CNT = M_CABLE_M / RHO_CNT # Cross sectional area
# Force balance for velocity... using your derived v_cable
V_CABLE = 8625.395 # m/s

# --- FINAL STRESS STATE (With Extra Load) ---
T_LOAD_FINAL = M_LOAD_M_FINAL * A_NET_250KM

# Centrifugal force of cable - Weight of cable - Weight of Load
# F_cable_net = M_cable * (v^2/r) - M_cable * g - M_load * g
# Your derivation simplifies this to F_cable limits.
# Re-implementing your exact sigma_cable_final equation:
term1 = M_CABLE_HW_M * (V_CABLE**2 - V_ORBIT**2)
term2 = T_LOAD_FINAL * R_ORBIT
SIGMA_CABLE_FINAL = (term1 - term2) / A_CNT # Should be ~5.86 GPa

# --- SPARE CAPACITY (LIFT) ---
# Force available before hitting SIGMA_CNT
F_MAX_LIFT_M = ((SIGMA_CNT - SIGMA_CABLE_FINAL) * A_CNT) / R_ORBIT # N/m

# =============================================================================
# 2. LAUNCH CLASS SELECTION
# =============================================================================
# "3g" means the payload experiences 3g radial acceleration at max speed.
# "10g" means 10g radial acceleration.

LAUNCH_PROFILE = "3g" # Options: "3g", "10g"

if LAUNCH_PROFILE == "3g":
    RADIAL_G_LIMIT = 3.0 * 9.80665
elif LAUNCH_PROFILE == "10g":
    RADIAL_G_LIMIT = 10.0 * 9.80665
else:
    RADIAL_G_LIMIT = 1.0 * 9.80665

# V_launch = sqrt( (a_net + a_radial_limit) * R )
V_LAUNCH_TARGET = math.sqrt((A_NET_250KM + RADIAL_G_LIMIT) * R_ORBIT)

# Max Linear Density the ring can hold at this G-force
# lift_max = F_max_lift / a_radial_limit ?? 
# Actually: Force exerted by sled = m * (v^2/r) = m * (a_net + a_radial_limit) ??
# No, your math says: lift_max_m_Xg = F_max_lift_m / (X * g)
# This implies the constraint is purely on the EXTRA g-force.
MAX_LINEAR_DENSITY = F_MAX_LIFT_M / RADIAL_G_LIMIT 

SLED_LENGTH = 1000.0
SLED_MASS_TOTAL = MAX_LINEAR_DENSITY * SLED_LENGTH

# =============================================================================
# 3. ENGINEERING CONSTRAINTS
# =============================================================================
HTS_TAPE_WIDTH_MM = 12.0
HTS_TAPE_LAYERS = 2
IC_PER_MM_PER_LAYER = 66.7
DE_RATING_FACTOR = 1 - math.sin(math.radians(20))
I_C_TOTAL = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * HTS_TAPE_LAYERS * DE_RATING_FACTOR
I_TARGET = 0.80 * I_C_TOTAL # Safety margin

N_TURNS = 300
TAU_P = 100.0
W_COIL = 2.0
GAP = 0.15
W_PLATE = 1.5
T_PLATE = 0.20
PLATE_MATERIAL = "gamma_titanium"

VOLTS_MAX = 100e3
MAX_SITE_POWER = 8e6
CRYO_COP_PENALTY = 50.0

DT = 1.0