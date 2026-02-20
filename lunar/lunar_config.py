#!/usr/bin/env python3
"""
Lunar Mass Driver Configuration — Parameters for lunar surface LSM simulation

This module contains all configurable parameters for a lunar surface mass driver
using LSM (Linear Synchronous Motor) technology identical to the orbital ring system.

Key differences from orbital ring:
  - Moon gravity (1.623 m/s²) instead of microgravity + centrifugal
  - Track follows lunar circumference (10,917 km loop)
  - Sled guided by rails (no magnetic bearings needed)
  - Above v_orbital (1,679 m/s), centrifugal exceeds gravity → rails hold sled down
  - Solar-powered HVDC grid (200 GW default)
  - Launch targets: lunar escape + interplanetary transfer velocities
  - No g-load throttling by default (unmanned cargo)

Motor physics are identical to the orbital ring LSM:
  - Same HTS stator coils (N=50, 12 mm × 5 layers)
  - Same sled SC field coils (B_sled = 0.10 T)
  - Same voltage limit (100 kV)
  - F/L = 7,500 N/m at I_TARGET, v_cross = 5,000 m/s, P/L = 37.5 MW/m

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: LUNAR PARAMETERS
# =============================================================================

R_MOON = 1_737_400          # m, lunar mean radius
G_MOON = 1.623              # m/s², surface gravitational acceleration
G_0 = 9.807                 # m/s², standard gravity at Earth sea level
C_MOON = 2 * math.pi * R_MOON  # m, lunar circumference (~10,917,000 m)
V_ORBITAL = math.sqrt(G_MOON * R_MOON)  # m/s, lunar orbital velocity (~1,679 m/s)
V_ESCAPE = math.sqrt(2.0 * G_MOON * R_MOON)  # m/s, lunar escape velocity (~2,375 m/s)

# =============================================================================
# SECTION 2: LSM STATOR PARAMETERS (identical to orbital ring)
# =============================================================================

# Pole pitch and layout
TAU_P = 130.0               # m, pole pitch — large pitch minimizes HTS hysteresis
N_POLES_PER_UNIT = 3        # pole pitches per repeating unit
L_GAP = 2.0                 # m, structural gap between repeating units
W_COIL = 2.0                # m, coil height (tangential extent)
G_GAP = 0.100               # m, air gap (stator face to sled coil face)

# Stator coil configuration
N_STATOR = 50               # turns per stator coil (fixed)

# HTS tape specification (Gömöry model) — 12 mm × 5 layers
ALPHA_TAPE_FIELD = 20.0                      # degrees — perpendicular field angle
SIN_ALPHA = math.sin(math.radians(ALPHA_TAPE_FIELD))
DE_RATING_FACTOR = 1 - SIN_ALPHA
TAPE_WIDTH_MM = 12.0                         # mm — REBCO tape width
TAPE_WIDTH = TAPE_WIDTH_MM / 1000.0          # m
I_C_PER_MM_LAYER = 66.7                      # A/mm-width/layer
N_HTS_LAYERS = 5
I_C_PER_LAYER_NOMINAL = I_C_PER_MM_LAYER * TAPE_WIDTH_MM   # 800 A
I_C_PER_LAYER = I_C_PER_LAYER_NOMINAL * DE_RATING_FACTOR   # 526.4 A
I_CRITICAL = I_C_PER_LAYER * N_HTS_LAYERS    # 2633 A
I_TARGET = 0.70 * I_CRITICAL                 # 1843 A normal operating current
I_PEAK = 0.80 * I_CRITICAL                   # 2107 A peak current
I_MIN = 10.0                # A, minimum current
N_STATOR_SIDES = 2           # Stator on both sides of sled
N_PHASES = 3                 # Three-phase system

# Voltage limit
V_COIL_LIMIT = 100_000      # V, coil insulation limit

# Sled SC field coils
B_SLED_NOMINAL = 0.10       # T, sled DC field at air gap
B_SLED_ADJUSTABLE = True    # Controller reduces B_sled at high speed for voltage limit

# =============================================================================
# SECTION 3: MOTOR CALIBRATION (simplified three-regime model)
#
# At I_TARGET (1,843 A), the thrust per metre of sled at the reference point
# (N=10, B_sled=0.10 T, delta_max=60°) is 1,500 N/m. Thrust scales linearly
# with both N and B_sled, giving:
#
#   F/L = 1500 × (N/10) × (B_sled/0.10)  N/m
#   v_cross = V_max / (2 × N × B_sled × w_coil)  m/s
#   P/L = F/L × v_cross  =  37.5 MW/m  (INVARIANT with N and B_sled)
#
# The three regimes:
#   1. v < v_cross: constant thrust F = F/L × L_sled
#   2. v_cross < v: constant power P = F × v_cross × L_sled  (voltage-limited)
#   3. P > P_HVDC_MAX: constant power P = P_HVDC_MAX  (grid-limited)
# =============================================================================

F_PER_M_CAL = 1500.0        # N/m baseline thrust per metre at N=10, B=0.10 (I_TARGET)
N_CAL = 10                   # reference stator turns for calibration
B_CAL = 0.10                 # reference sled field for calibration

# =============================================================================
# SECTION 4: SLED PARAMETERS
# =============================================================================

L_SLED = 10_000             # m, sled length (10 km)
M_SLED_PER_M = 400          # kg/m, sled hardware linear density
# Breakdown (lighter than orbital ring — rails replace magnetic bearings):
#   SC field coils (76 LSM modules)          30 kg/m
#   Rail guide system (replaces mag bearings) 40 kg/m
#   Cryostats (LSM coils)                    60 kg/m
#   LN2 reservoir                              5 kg/m
#   Structural frame                         200 kg/m
#   Spacecraft rails + locking mechanisms     35 kg/m
#   SC coil structural mounts + alignment     15 kg/m
#   Misc (wiring, sensors, thermal shields)   15 kg/m
M_SPACECRAFT = 5_000_000    # kg, payload mass (5,000 tonnes)

# Rail structural limit — maximum outward force per metre of track
# Set to None to disable (just report maximum rail force)
F_RAIL_MAX_PER_M = None     # N/m

# =============================================================================
# SECTION 5: HVDC POWER GRID (Solar)
# =============================================================================

# The lunar mass driver is powered by solar panels along the equatorial track.
# A 100 m wide strip of panels running the full 10,917 km circumference
# generates ~200 GW at lunar noon. Power varies with solar angle and is
# stored/buffered for launch operations.
P_HVDC_MAX = 200e9          # W (200 GW), set to None to disable

# =============================================================================
# SECTION 6: MISSION PARAMETERS
# =============================================================================

V_LAUNCH = 5_000            # m/s, default target velocity
A_MAX_G = 10.0              # maximum tangential acceleration in g (thrust-limited)
G_LOAD_MAX = None           # g-load limit (None = no limit, for unmanned cargo)
DT = 0.1                    # s, timestep
DT_QUICK = 1.0              # s, timestep for --quick mode

# Load angle limits
DELTA_MAX = math.pi / 3     # rad, maximum load angle (60°)
DELTA_TARGET = math.pi / 6  # rad, target load angle (30°)

# =============================================================================
# SECTION 7: CELESTIAL MECHANICS CONSTANTS
# =============================================================================

MU_SUN = 1.32712440018e20   # m³/s², Sun gravitational parameter
MU_EARTH = 3.986004418e14   # m³/s², Earth gravitational parameter
MU_MOON = 4.9048695e12      # m³/s², Moon gravitational parameter
AU = 1.495978707e11         # m, astronomical unit

# Planetary semi-major axes
A_EARTH = 1.000 * AU        # m, Earth
A_MARS = 1.524 * AU         # m, Mars
A_JUPITER = 5.204 * AU      # m, Jupiter
A_SATURN = 9.583 * AU       # m, Saturn

V_EARTH_ORBITAL = math.sqrt(MU_SUN / A_EARTH)  # ~29,783 m/s

# Moon orbital parameters (around Earth)
A_MOON_ORBIT = 384_400_000  # m, Moon semi-major axis around Earth
V_MOON_ORBITAL = math.sqrt(MU_EARTH / A_MOON_ORBIT)  # ~1,018 m/s
V_ESC_EARTH_AT_MOON = math.sqrt(2.0 * MU_EARTH / A_MOON_ORBIT)  # ~1,440 m/s

# =============================================================================
# SECTION 8: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7    # H/m, permeability of free space
STEFAN_BOLTZMANN = 5.670374e-8  # W/m²K⁴

# =============================================================================
# SECTION 9: DERIVED PARAMETERS
# =============================================================================

def calc_derived():
    """Calculate all derived quantities from base parameters."""
    global L_UNIT, F_ACTIVE, N_UNITS, A_ACTIVE_TOTAL, W_ACTIVE_PER_UNIT
    global M_SLED_HARDWARE, M_TOTAL, M_PER_M_TOTAL, KE_TARGET
    global F_MAX_SUPPLY, A_STATOR_COIL, L_TRACK
    global F_PER_M, V_CROSSOVER, P_PER_M, F_TOTAL, P_TOTAL

    # Repeating unit geometry
    L_UNIT = N_POLES_PER_UNIT * TAU_P + L_GAP
    W_ACTIVE_PER_UNIT = N_POLES_PER_UNIT * TAU_P
    F_ACTIVE = W_ACTIVE_PER_UNIT / L_UNIT

    # Number of units on sled
    N_UNITS = int(L_SLED / L_UNIT)

    # Total active area (both sides)
    A_ACTIVE_TOTAL = 2 * N_UNITS * W_ACTIVE_PER_UNIT * W_COIL
    A_STATOR_COIL = TAU_P * W_COIL

    # Sled mass
    M_SLED_HARDWARE = L_SLED * M_SLED_PER_M
    M_TOTAL = M_SLED_HARDWARE + M_SPACECRAFT
    M_PER_M_TOTAL = M_SLED_PER_M + M_SPACECRAFT / L_SLED  # total kg/m for rail force

    # Motor calibration (three-regime model)
    F_PER_M = F_PER_M_CAL * (N_STATOR / N_CAL) * (B_SLED_NOMINAL / B_CAL)
    V_CROSSOVER = V_COIL_LIMIT / (2.0 * N_STATOR * B_SLED_NOMINAL * W_COIL)
    P_PER_M = F_PER_M * V_CROSSOVER  # invariant: 37.5 MW/m

    F_TOTAL = F_PER_M * L_SLED
    P_TOTAL = P_PER_M * L_SLED

    # Target KE
    KE_TARGET = 0.5 * M_TOTAL * V_LAUNCH**2

    # Max supply frequency
    F_MAX_SUPPLY = V_LAUNCH / (2 * TAU_P)

    # Track length
    L_TRACK = C_MOON

    return {
        'L_unit': L_UNIT,
        'f_active': F_ACTIVE,
        'n_units': N_UNITS,
        'A_active_total': A_ACTIVE_TOTAL,
        'M_sled_hardware': M_SLED_HARDWARE,
        'M_total': M_TOTAL,
        'M_per_m_total': M_PER_M_TOTAL,
        'KE_target': KE_TARGET,
        'f_max_supply': F_MAX_SUPPLY,
        'A_stator_coil': A_STATOR_COIL,
        'L_track': L_TRACK,
        'F_per_m': F_PER_M,
        'V_crossover': V_CROSSOVER,
        'P_per_m': P_PER_M,
        'F_total': F_TOTAL,
        'P_total': P_TOTAL,
    }

# Initialize derived parameters
calc_derived()


def get_physics_params():
    """Return parameter dictionary for physics functions."""
    return {
        'n_stator': N_STATOR,
        'tau_p': TAU_P,
        'w_coil': W_COIL,
        'g_gap': G_GAP,
        'b_sled_nominal': B_SLED_NOMINAL,
        'v_coil_limit': V_COIL_LIMIT,
        'a_active_total': A_ACTIVE_TOTAL,
        'a_stator_coil': A_STATOR_COIL,
        'delta_max': DELTA_MAX,
        'm_total': M_TOTAL,
        'r_moon': R_MOON,
        'g_moon': G_MOON,
        # HTS Gömöry hysteresis parameters
        'n_turns': N_STATOR,
        'tape_width_mm': TAPE_WIDTH_MM,
        'n_hts_layers': N_HTS_LAYERS,
        'i_c_per_layer': I_C_PER_LAYER,
        'sin_alpha': SIN_ALPHA,
        'n_phases': N_PHASES,
        'n_stator_sides': N_STATOR_SIDES,
    }

# =============================================================================
# SECTION 10: DISPLAY CONFIGURATION
# =============================================================================

SAVE_GRAPHS = True
GRAPH_OUTPUT_DIR = "./graphs_lunar"
GRAPH_DPI = 300
GRAPH_WIDTH_INCHES = 10
GRAPH_HEIGHT_INCHES = 6
GRAPH_FORMAT = "png"
