#!/usr/bin/env python3
"""
LIM Configuration Module - Parameters and material properties

This module contains all configurable parameters for the orbital ring
deployment simulation.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: USER-CONFIGURABLE PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 HTS (High-Temperature Superconductor) Tape Configuration
# -----------------------------------------------------------------------------
HTS_TAPE_WIDTH_MM = 3       # Standard widths: 12, 6, 4, or 3 mm
HTS_TAPE_LAYERS = 3         # Layers of tape in parallel: 1, 2, or 3
IC_PER_MM_PER_LAYER = 66.7  # Critical current density: ~66.7 A/mm per layer
NORRIS_HYSTERESIS = False   # Use Norris Formula for hysteresis loss

# -----------------------------------------------------------------------------
# 1.2 LIM Geometry
# -----------------------------------------------------------------------------
N_TURNS = 80            # Turns per phase coil (typically 50-200)
TAU_P = 100.0           # Pole pitch in meters
W_COIL = 0.5            # LIM coil width in meters
GAP = 0.20              # Air gap between coil and reaction plate (m)
T_PLATE = 0.2           # Reaction plate thickness (m)
PITCH_COUNT = 3         # Number of pole pitches per LIM

# -----------------------------------------------------------------------------
# 1.2b Plate Geometry
# -----------------------------------------------------------------------------
W_PLATE = 1.2               # Reaction plate width (m)
N_PLATES_PER_SIDE = 1       # LIM reaction plates per side of cable
N_LEV_PLATES = 1            # Levitation plates per side
D_LEV = 0.10                # Levitation plate thickness (m)
W_LEV = 1.4                 # Levitation plate width (m)

# Reaction plate material: "aluminum", "cuni7030", "titanium", "alpha_titanium", "gamma_titanium"
PLATE_MATERIAL = "titanium"

# -----------------------------------------------------------------------------
# 1.3 LIM Spacing and Site Configuration
# -----------------------------------------------------------------------------
LIM_SPACING = 500.0     # Distance between LIM sites (m)

# -----------------------------------------------------------------------------
# 1.4 Operating Limits
# -----------------------------------------------------------------------------
VOLTS_MAX = 100e3       # Maximum induced coil voltage (V)
MAX_SITE_POWER = 16e6   # Maximum power per LIM site (W)

# -----------------------------------------------------------------------------
# 1.5 Slip Control Parameters
# -----------------------------------------------------------------------------
V_SLIP_MIN = 5.0        # Minimum slip velocity (m/s)
V_SLIP_MAX = 200.0      # Maximum slip velocity (m/s)
SLIP_RATIO_NORMAL = 0.02    # Target slip ratio at full current (2%)
SLIP_RATIO_REDUCED = 0.005  # Reduced slip when power-limited (0.5%)

# -----------------------------------------------------------------------------
# 1.6 Thrust Model Selection
# -----------------------------------------------------------------------------
THRUST_MODEL = 1        # Select model: 1, 2, or 3
THRUST_EFFICIENCY = 1.0 # Multiplier on calculated thrust

# -----------------------------------------------------------------------------
# 1.7 Mass Configuration
# -----------------------------------------------------------------------------
M_CABLE_STRUCTURAL = 96_700     # Structural cable mass per meter (kg/m)
M_LOAD_M = 12_000               # Casing + payload mass per meter (kg/m)

# -----------------------------------------------------------------------------
# 1.8 Simulation Control
# -----------------------------------------------------------------------------
SAMPLE_TIME_MAX = 5 * 365.33 * 24 * 3600  # Maximum simulation time (5 years)
DT1 = 1                 # Time step for first day (s)
DT2 = 10                # Time step until sample_time_max (s)
DT3 = 50                # Time step after sample_time_max (s)
SKIP = 200              # Data collection interval

# -----------------------------------------------------------------------------
# 1.9 Output Control
# -----------------------------------------------------------------------------
WRITE_FILE = True
MAKE_GRAPHS = True

# -----------------------------------------------------------------------------
# 1.10 Thermal Mode Configuration (NEW)
# -----------------------------------------------------------------------------
EDDY_HEAT_TO_CABLE = False      # True = cable heats up, False = cryo handles it

# Cable thermal properties
CABLE_EMISSIVITY = 0.85
CABLE_SURFACE_AREA_PER_M = 0.5  # Radiating surface area per meter (m²/m)

# Coil thermal isolation
COIL_MLI_EFFECTIVENESS = 0.001
COIL_SURFACE_AREA_PER_SITE = 10  # Approximate coil cryostat surface area (m²)


# =============================================================================
# SECTION 2: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7
STEFAN_BOLTZMANN = 5.670374e-8

# Orbital parameters at 250 km altitude
V_ORBIT = 7754.866
V_GROUND_STATIONARY = 483.331
L_RING = 41_645_813.012

# =============================================================================
# SECTION 3: MATERIAL PROPERTIES
# =============================================================================

# Aluminum
RHO_ALU_293K = 2.65e-8
RHO_ALU_70K = 4.853e-9
DENSITY_ALU = 2700
C_P_ALU = 900
K_ALU = 205
EM_ALU = 0.85
ALPHA_ALU = 3.663e-3

# CuNi 70/30 (Cupronickel)
RHO_CUNI_293K = 38e-8
DENSITY_CUNI = 8900
C_P_CUNI = 377
K_CUNI = 29
EM_CUNI = 0.65
ALPHA_CUNI = 0.0004

# Pure Titanium
RHO_TI_293K = 42e-8
DENSITY_TI = 4500
C_P_TI = 520
K_TI = 22
EM_TI = 0.60
ALPHA_TI = 0.0035

# Alpha-2 titanium aluminide, α₂-Ti₃Al
RHO_aTI_293K = 50e-8
DENSITY_aTI = 4200
C_P_aTI = 550
K_aTI = 17
EM_aTI = 0.60
ALPHA_aTI = 0.0015

# Gamma titanium aluminide, Ti-48Al-2Cr-2Nb
RHO_gTI_293K = 75e-8
DENSITY_gTI = 3900
C_P_gTI = 570
K_gTI = 15
EM_gTI = 0.60
ALPHA_gTI = 0.0012

# Iron (for levitation plates)
DENSITY_IRON = 7870

# Liquid nitrogen
T_LN2_BOIL = 77.4
T_LN2_SUPPLY = 70
C_P_LN2 = 2040
L_V_LN2 = 199000
H_CONV_LN2 = 90

# Thermal environment
Q_SUN = 1361
Q_EARTH_ALBEDO = 650
Q_SHIELDING = 0.005
T_SPACE = 2.7
T_RADIATOR_HOT = 300
T_MAX_PLATE = 500

# Cryogenic system
CRYO_EFF = 0.18
EM_HEATSINK = 0.9

# HTS tape
HTS_THICKNESS_UM = 80
ALPHA_PENETRATION_DEG = 20.0


# =============================================================================
# SECTION 4: DERIVED PARAMETERS
# =============================================================================

# Time constants
HR = 3600
DAY = 24 * HR
WEEK = 7 * DAY
MONTH = 30 * DAY
YR = round(365.33 * DAY)

# HTS current ratings
W_TAPE = HTS_TAPE_WIDTH_MM / 1000
I_C = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * HTS_TAPE_LAYERS
I_PEAK = 0.875 * I_C
I_TARGET = 0.8125 * I_C
I_MIN = 10.0

# LIM geometry
LIMS_PER_SITE = 2 * N_PLATES_PER_SIDE
L_ACTIVE = TAU_P * PITCH_COUNT
A_LIM = L_ACTIVE * W_COIL
A_COIL = TAU_P * W_COIL
L_HTS_COIL = 2 * (W_COIL + TAU_P) * N_TURNS
L_HTS_LIM = L_HTS_COIL * 3 * PITCH_COUNT
LIM_PHASES = 3
LIM_SITES = round(L_RING / LIM_SPACING)

# Angular conversion
ALPHA_TAPE = ALPHA_PENETRATION_DEG * math.pi / 180

# Material selection
def get_material_properties(material):
    """Return material properties dict for given material name."""
    if material == "cuni7030":
        return {
            'density': DENSITY_CUNI, 'rho_293K': RHO_CUNI_293K,
            'alpha': ALPHA_CUNI, 'cp': C_P_CUNI, 'k': K_CUNI, 'em': EM_CUNI
        }
    elif material == "titanium":
        return {
            'density': DENSITY_TI, 'rho_293K': RHO_TI_293K,
            'alpha': ALPHA_TI, 'cp': C_P_TI, 'k': K_TI, 'em': EM_TI
        }
    elif material == "alpha_titanium":
        return {
            'density': DENSITY_aTI, 'rho_293K': RHO_aTI_293K,
            'alpha': ALPHA_aTI, 'cp': C_P_aTI, 'k': K_aTI, 'em': EM_aTI
        }
    elif material == "gamma_titanium":
        return {
            'density': DENSITY_gTI, 'rho_293K': RHO_gTI_293K,
            'alpha': ALPHA_gTI, 'cp': C_P_gTI, 'k': K_gTI, 'em': EM_gTI
        }
    else:  # aluminum
        return {
            'density': DENSITY_ALU, 'rho_293K': RHO_ALU_293K,
            'alpha': ALPHA_ALU, 'cp': C_P_ALU, 'k': K_ALU, 'em': EM_ALU
        }

# Get current material properties
MATERIAL_PROPS = get_material_properties(PLATE_MATERIAL)
PLATE_DENSITY = MATERIAL_PROPS['density']
PLATE_RHO_293K = MATERIAL_PROPS['rho_293K']
PLATE_ALPHA = MATERIAL_PROPS['alpha']
PLATE_CP = MATERIAL_PROPS['cp']
PLATE_K = MATERIAL_PROPS['k']
PLATE_EM = MATERIAL_PROPS['em']

# Mass calculations
M_LIM_PLATE = PLATE_DENSITY * T_PLATE * W_PLATE
M_LEV_PLATE = DENSITY_IRON * D_LEV * W_LEV
M_HARDWARE = (2 * N_LEV_PLATES * M_LEV_PLATE + 2 * N_PLATES_PER_SIDE * M_LIM_PLATE)
M_CABLE_M = M_CABLE_STRUCTURAL + M_HARDWARE

# Total masses
M_CABLE_TOTAL = M_CABLE_M * L_RING
M_LOAD_TOTAL = M_LOAD_M * L_RING

# Thermal
CASING_WIDTH = 10.0
Q_ABSORBED_PER_M = (Q_SUN + Q_EARTH_ALBEDO) * CASING_WIDTH * Q_SHIELDING
Q_ABSORBED_PER_SITE = Q_ABSORBED_PER_M * LIM_SPACING
MAX_HEATSINK_AREA = LIM_SPACING * 2 * W_COIL
HEATSINK_LENGTH = LIM_SPACING
V_REL_MIN_FUDGE = 10

# Coil environmental heat leak
Q_COIL_ENVIRONMENT = (Q_SUN + Q_EARTH_ALBEDO) * COIL_SURFACE_AREA_PER_SITE * COIL_MLI_EFFECTIVENESS

# Efficiency factors
INV_EFF = 0.90
LIM_EFF = 0.95

# Controller parameters
SLIP_MIN = 0.005
CURRENT_UPRATE = 1.01
POWER_HEADROOM = 0.98


# =============================================================================
# SECTION 5: PARAMETER DICTIONARY FOR PHYSICS FUNCTIONS
# =============================================================================

def get_physics_params():
    """Return parameter dictionary for physics functions."""
    return {
        'n_turns': N_TURNS,
        'tau_p': TAU_P,
        'w_coil': W_COIL,
        'gap': GAP,
        't_plate': T_PLATE,
        'pitch_count': PITCH_COUNT,
        'w_plate': W_PLATE,
        'l_active': L_ACTIVE,
        'a_lim': A_LIM,
        'a_coil': A_COIL,
        'l_hts_coil': L_HTS_COIL,
        'lim_phases': LIM_PHASES,
        'rho_293K': PLATE_RHO_293K,
        'alpha': PLATE_ALPHA,
        'material': PLATE_MATERIAL,
        'density': PLATE_DENSITY,
        'cp': PLATE_CP,
        'em': PLATE_EM,
        'v_slip_min': V_SLIP_MIN,
        'thrust_efficiency': THRUST_EFFICIENCY,
        'w_tape': W_TAPE,
        'alpha_tape': ALPHA_TAPE,
        'i_c': I_C,
    }


# =============================================================================
# SECTION 6: DISPLAY DICTIONARY
# =============================================================================

PARAM_DISPLAY = {
    "N_TURNS             #": N_TURNS,
    "TAU_P               m": TAU_P,
    "W_COIL              m": W_COIL,
    "SLIP_RATIO_NORMAL   %": SLIP_RATIO_NORMAL * 100,
    "SLIP_RATIO_REDUCED  %": SLIP_RATIO_REDUCED * 100,
    "GAP                mm": GAP * 1000,
    "HTS_TAPE_WIDTH     mm": HTS_TAPE_WIDTH_MM,
    "HTS_TAPE_LAYERS     #": HTS_TAPE_LAYERS,
    "I_C                 A": round(I_C),
    "I_PEAK              A": round(I_PEAK),
    "I_TARGET            A": round(I_TARGET),
    "V_SLIP_MAX        m/s": V_SLIP_MAX,
    "V_SLIP_MIN        m/s": V_SLIP_MIN,
    "N_PLATES_PER_SIDE   #": N_PLATES_PER_SIDE,
    "LIMS_PER_SITE       #": LIMS_PER_SITE,
    "PITCH_COUNT         #": PITCH_COUNT,
    "LIM_SPACING         m": LIM_SPACING,
    "LIM_SITES           #": LIM_SITES,
    "THRUST_EFFICIENCY   %": THRUST_EFFICIENCY * 100,
    "VOLTS_MAX          kV": VOLTS_MAX / 1000,
    "MAX_SITE_POWER     MW": MAX_SITE_POWER / 1e6,
    "T_PLATE            mm": T_PLATE * 1000,
    "NORRIS_HYSTERESIS    ": NORRIS_HYSTERESIS,
    "EDDY_HEAT_TO_CABLE   ": EDDY_HEAT_TO_CABLE,
    "M_CABLE_STRUCT   kg/m": M_CABLE_STRUCTURAL,
    "M_HARDWARE       kg/m": M_HARDWARE,
    "M_CABLE_TOTAL    kg/m": M_CABLE_M,
    "M_LOAD           kg/m": M_LOAD_M,
}


def print_parameters():
    """Display current parameter configuration."""
    streq = "="*70
    print(f"\n{streq}")
    print("SIMULATION PARAMETERS")
    print(streq)
    print(f"  {'PLATE_MATERIAL':24} {PLATE_MATERIAL:>12}")
    print(f"  {'EDDY_HEAT_TO_CABLE':24} {str(EDDY_HEAT_TO_CABLE):>12}")
    for key, value in PARAM_DISPLAY.items():
        if isinstance(value, bool):
            print(f"  {key:24} {str(value):>12}")
        else:
            print(f"  {key:24} {value:>12.2f}")
    print(f"{streq}\n")
