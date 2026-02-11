#!/usr/bin/env python3
"""
LSM Configuration Module - Parameters for LSM mass driver simulation

This module contains all configurable parameters for the orbital ring
LSM (Linear Synchronous Motor) mass driver launch system.

Key difference from LIM: The LSM uses superconducting DC field coils on the sled
that interact directly with the track's AC stator. No eddy currents, no slip,
no thermal limit on the sled. The voltage limit is the primary constraint.

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
TAU_P = 20.0                # m, pole pitch (single value for entire launch)
N_POLES_PER_UNIT = 3        # pole pitches per repeating unit
L_GAP = 2.0                 # m, structural gap between repeating units
W_COIL = 2.0                # m, coil height (tangential extent)
G_GAP = 0.100               # m, air gap (stator face to sled coil face)

# Stator coil configuration
N_STATOR = 10               # turns per stator coil

# HTS tape specification (same as LIM)
TAPE_WIDTH = 0.012          # m (12 mm)
N_LAYERS = 3                # layers of tape in parallel
IC_PER_MM_LAYER = 66.7      # A/mm-width/layer (critical current density)
I_CRITICAL = IC_PER_MM_LAYER * (TAPE_WIDTH * 1000) * N_LAYERS  # = 2400 A
I_TARGET = 0.80 * I_CRITICAL   # 1920 A operating current
I_PEAK = 0.90 * I_CRITICAL     # 2160 A peak current
I_MIN = 10.0                # A, minimum current

# Voltage limit
V_COIL_LIMIT = 100_000      # V, coil insulation limit

# =============================================================================
# SECTION 3: SLED SC FIELD COILS
# =============================================================================

# Sled superconducting DC field
B_SLED_NOMINAL = 0.10       # T, sled DC field at air gap (base value)
B_SLED_ADJUSTABLE = False   # If True, controller can reduce B_sled to meet voltage limit

# =============================================================================
# SECTION 4: SLED MASS
# =============================================================================

L_SLED = 5000               # m, sled length
M_SC_COILS_PER_M = 30       # kg/m, SC field coil linear density
M_STRUCTURE_PER_M = 100     # kg/m, structural frame + bearings + cryostat
M_SPACECRAFT = 500_000      # kg, payload mass

# =============================================================================
# SECTION 5: MISSION PARAMETERS
# =============================================================================

V_LAUNCH = 15_000           # m/s, target launch velocity
A_MAX_G = 0.5               # maximum acceleration in g
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

# =============================================================================
# SECTION 7: DERIVED PARAMETERS
# =============================================================================

def calc_derived():
    """Calculate all derived quantities from base parameters."""
    global L_UNIT, F_ACTIVE, N_UNITS, A_ACTIVE_TOTAL, W_ACTIVE_PER_UNIT
    global M_SLED_HARDWARE, M_TOTAL, KE_TARGET, F_MAX_SUPPLY
    global A_STATOR_COIL, G_EFF

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

    # Effective gap (including coil thickness effects)
    G_EFF = G_GAP  # For now, just the physical gap

    # Sled mass
    M_SLED_HARDWARE = L_SLED * (M_SC_COILS_PER_M + M_STRUCTURE_PER_M)  # kg
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
        'g_eff': G_EFF,
        'b_sled_nominal': B_SLED_NOMINAL,
        'b_sled_adjustable': B_SLED_ADJUSTABLE,
        'v_coil_limit': V_COIL_LIMIT,
        'a_active_total': A_ACTIVE_TOTAL,
        'a_stator_coil': A_STATOR_COIL,
        'delta_max': DELTA_MAX,
        'm_total': M_TOTAL,
        'r_orbit': R_ORBIT,
        'g_250': G_250,
    }

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
