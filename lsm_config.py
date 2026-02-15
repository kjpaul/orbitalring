#!/usr/bin/env python3
"""
LSM Configuration Module - Parameters for LSM mass driver simulation

This module contains all configurable parameters for the orbital ring
LSM (Linear Synchronous Motor) mass driver launch system.

Key difference from LIM: The LSM uses superconducting DC field coils on the sled
that interact directly with the track's AC stator. No eddy currents, no slip,
no thermal limit on the sled. The voltage limit is the primary constraint.

Design choices:
  - τ_p = 130 m: Large pole pitch minimizes HTS hysteresis losses, which scale
    as ~1/τ_p³ (B_coil_peak ~ N×I/l_coil where l_coil = τ_p/3).
  - The LSM requires periodic alternating N/S SC poles on the sled (not uniform
    DC field). The stator AC field locks to the sled pole pattern synchronously.

Launch modes:
  - "crew": v_target = 15,000 m/s (Mars transfer), 3g human physiological limit
  - "cargo": v_target = 30,000 m/s (deep solar system), 13g structural limit

Sled architecture:
  - "fixed" (Option A): Single pole pitch, sled SC coils at τ_p. Simple, default.
  - "reconfig" (Option B): Reconfigurable coils at base pitch τ_base, grouped to
    create effective pitch = n × τ_base. Multi-stage operation like the LIM.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: ORBITAL RING PARAMETERS
# =============================================================================

# Orbital parameters at 250 km altitude
R_ORBIT = 6_628_137         # m, orbital radius (Earth radius + 250 km)
G_250 = 9.073               # m/s², gravitational acceleration at 250 km
G_0 = 9.807                 # m/s², standard gravity at sea level
V_ORBIT = 7755              # m/s, orbital velocity at 250 km
L_RING = 41_646_000         # m, ring circumference

# =============================================================================
# SECTION 2: LSM STATOR PARAMETERS
# =============================================================================

# Pole pitch and layout
TAU_P = 130.0               # m, pole pitch — large pitch minimizes hysteresis (P ~ 1/τ_p³)
N_POLES_PER_UNIT = 3        # pole pitches per repeating unit
L_GAP = 2.0                 # m, structural gap between repeating units
W_COIL = 2.0                # m, coil height (tangential extent)
G_GAP = 0.100               # m, air gap (stator face to sled coil face)

# Stator coil configuration
N_STATOR = 10               # turns per stator coil

# HTS tape specification (Gömöry model)
# 12 mm × 4 layers: simpler winding than LIM's 3 mm × 16 layers.
# Hysteresis is 4× higher but still negligible for LSM (kW not GW).
ALPHA_TAPE_FIELD = 20.0                      # degrees — Gömöry perpendicular field angle
SIN_ALPHA = math.sin(math.radians(ALPHA_TAPE_FIELD))  # 0.34202
DE_RATING_FACTOR = 1 - SIN_ALPHA             # 0.658 — Ic reduction per layer
TAPE_WIDTH_MM = 12.0                         # mm — 12 mm tape (simpler winding, 4 layers)
TAPE_WIDTH = TAPE_WIDTH_MM / 1000.0          # m
I_C_PER_MM_LAYER = 66.7                      # A/mm-width/layer (critical current density)
N_HTS_LAYERS = 4                             # layers of 12 mm tape
I_C_PER_LAYER_NOMINAL = I_C_PER_MM_LAYER * TAPE_WIDTH_MM   # 800 A nominal
I_C_PER_LAYER = I_C_PER_LAYER_NOMINAL * DE_RATING_FACTOR   # 526.4 A effective
I_CRITICAL = I_C_PER_LAYER * N_HTS_LAYERS    # 2106 A total
I_TARGET = 0.80 * I_CRITICAL                 # 1685 A operating current
I_PEAK = 0.90 * I_CRITICAL                   # 1895 A peak current
I_MIN = 10.0                # A, minimum current
N_LIM_SIDES = 2                # Stator on both sides of sled
LIM_PHASES = 3                 # Three-phase system

# Voltage limit
V_COIL_LIMIT = 100_000      # V, coil insulation limit

# =============================================================================
# SECTION 3: SLED SC FIELD COILS
# =============================================================================

# Sled superconducting DC field
B_SLED_NOMINAL = 0.10       # T, sled DC field at air gap (base value)
B_SLED_ADJUSTABLE = True    # Controller reduces B_sled at high speed to meet voltage limit

# Sled architecture: "fixed" (Option A) or "reconfig" (Option B)
SLED_ARCHITECTURE = "fixed"

# Option B reconfigurable coil parameters
TAU_BASE = 43.0              # m, base SC coil spacing on sled
# Effective pitches available: 43, 86, 129 m (1×, 2×, 3× base)

RECONFIG_STAGES = [
    {"name": "RC1", "effective_pitch": TAU_BASE,       "f_handoff_hz": 30.0},
    {"name": "RC2", "effective_pitch": 2 * TAU_BASE,   "f_handoff_hz": 60.0},
    {"name": "RC3", "effective_pitch": 3 * TAU_BASE,   "f_handoff_hz": None},  # final stage
]

# =============================================================================
# SECTION 4: SLED MASS
# =============================================================================

L_SLED = 10_000             # m, sled length (10 km for adequate thrust)
M_SLED_PER_M = 525          # kg/m, total sled hardware linear density
# Breakdown:
#   SC field coils (76 LSM modules)          30 kg/m
#   Bearing system (EML+EDL, top+bottom)    110 kg/m
#   Cryostats (LSM coils + EML bearings)     60 kg/m
#   LN2 reservoir                             5 kg/m
#   Structural frame                        250 kg/m
#   Spacecraft rails + locking mechanisms    40 kg/m
#   SC coil structural mounts + alignment    15 kg/m
#   Misc (wiring, sensors, thermal shields)  15 kg/m
M_SPACECRAFT = 5_000_000    # kg, payload mass (5,000 tonnes)

# =============================================================================
# SECTION 5: MISSION PARAMETERS
# =============================================================================

# Launch mode: "crew" or "cargo"
LAUNCH_MODE = "cargo"

# Mode-dependent parameters
if LAUNCH_MODE == "crew":
    V_LAUNCH = 15_000       # m/s — Mars transfer velocity
    G_LOAD_MAX = 3.0        # g — human physiological limit
elif LAUNCH_MODE == "cargo":
    V_LAUNCH = 30_000       # m/s — deep solar system resupply
    G_LOAD_MAX = 13.0       # g — structural limit only

A_MAX_G = 0.5               # maximum tangential acceleration in g
DT = 0.1                    # s, timestep (default)
DT_QUICK = 1.0              # s, timestep for --quick mode

# Load angle limits
DELTA_MAX = math.pi / 3     # rad, maximum load angle (60°) for stability margin
DELTA_TARGET = math.pi / 6  # rad, target load angle (30°) for normal operation

# =============================================================================
# SECTION 6: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7    # H/m, permeability of free space
STEFAN_BOLTZMANN = 5.670374e-8  # W/m²K⁴, Stefan-Boltzmann constant
T_SPACE = 2.7               # K, deep space temperature
T_STATOR = 77.0             # K, stator coil operating temperature (LN2)

# =============================================================================
# SECTION 7: DERIVED PARAMETERS
# =============================================================================

def calc_derived():
    """Calculate all derived quantities from base parameters."""
    global L_UNIT, F_ACTIVE, N_UNITS, A_ACTIVE_TOTAL, W_ACTIVE_PER_UNIT
    global M_SLED_HARDWARE, M_TOTAL, KE_TARGET, F_MAX_SUPPLY
    global A_STATOR_COIL

    # Repeating unit geometry
    L_UNIT = N_POLES_PER_UNIT * TAU_P + L_GAP  # m, length of one repeating unit
    W_ACTIVE_PER_UNIT = N_POLES_PER_UNIT * TAU_P  # m, active length per unit
    F_ACTIVE = W_ACTIVE_PER_UNIT / L_UNIT  # active fraction

    # Number of units on sled
    N_UNITS = int(L_SLED / L_UNIT)

    # Total active area (both sides of sled)
    A_ACTIVE_TOTAL = 2 * N_UNITS * W_ACTIVE_PER_UNIT * W_COIL  # m²

    # Stator coil area
    A_STATOR_COIL = TAU_P * W_COIL  # m²

    # Sled mass
    M_SLED_HARDWARE = L_SLED * M_SLED_PER_M  # kg
    M_TOTAL = M_SLED_HARDWARE + M_SPACECRAFT  # kg

    # Target kinetic energy
    KE_TARGET = 0.5 * M_TOTAL * V_LAUNCH**2  # J

    # Maximum supply frequency (at target velocity)
    F_MAX_SUPPLY = V_LAUNCH / (2 * TAU_P)  # Hz

    return {
        'L_unit': L_UNIT,
        'f_active': F_ACTIVE,
        'n_units': N_UNITS,
        'A_active_total': A_ACTIVE_TOTAL,
        'M_sled_hardware': M_SLED_HARDWARE,
        'M_total': M_TOTAL,
        'KE_target': KE_TARGET,
        'f_max_supply': F_MAX_SUPPLY,
        'A_stator_coil': A_STATOR_COIL,
        'w_active_per_unit': W_ACTIVE_PER_UNIT,
    }

# Initialize derived parameters
calc_derived()

# =============================================================================
# SECTION 8: PARAMETER DICTIONARY FOR PHYSICS FUNCTIONS
# =============================================================================

def get_physics_params():
    """Return parameter dictionary for physics functions."""
    return {
        'n_stator': N_STATOR,
        'tau_p': TAU_P,
        'w_coil': W_COIL,
        'g_gap': G_GAP,
        'b_sled_nominal': B_SLED_NOMINAL,
        'b_sled_adjustable': B_SLED_ADJUSTABLE,
        'v_coil_limit': V_COIL_LIMIT,
        'a_active_total': A_ACTIVE_TOTAL,
        'a_stator_coil': A_STATOR_COIL,
        'delta_max': DELTA_MAX,
        'm_total': M_TOTAL,
        'r_orbit': R_ORBIT,
        'g_250': G_250,
        # HTS Gömöry hysteresis parameters
        'n_turns': N_STATOR,
        'tape_width_mm': TAPE_WIDTH_MM,
        'n_hts_layers': N_HTS_LAYERS,
        'i_c_per_layer': I_C_PER_LAYER,
        'sin_alpha': SIN_ALPHA,
        'lim_phases': LIM_PHASES,
        'n_lim_sides': N_LIM_SIDES,
    }

def apply_launch_mode():
    """Apply launch mode settings (call after CLI overrides)."""
    global V_LAUNCH, G_LOAD_MAX
    if LAUNCH_MODE == "crew":
        V_LAUNCH = 15_000
        G_LOAD_MAX = 3.0
    elif LAUNCH_MODE == "cargo":
        V_LAUNCH = 30_000
        G_LOAD_MAX = 13.0

# =============================================================================
# SECTION 9: DISPLAY CONFIGURATION
# =============================================================================

# Output control
SAVE_GRAPHS = True
GRAPH_OUTPUT_DIR = "./graphs_lsm"
GRAPH_DPI = 300
GRAPH_WIDTH_INCHES = 10
GRAPH_HEIGHT_INCHES = 6
GRAPH_FORMAT = "png"
