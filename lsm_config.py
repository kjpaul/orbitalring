#!/usr/bin/env python3
"""
Parameters of the LSM mass driver.
Book: Orbital Ring Engineering, Volume I: Chapter 10 (imported by lsm_simulation.py).
Run: imported; not run by itself
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

LSM Configuration Module - Parameters for LSM mass driver simulation

This module contains all configurable parameters for the orbital ring
LSM (Linear Synchronous Motor) mass driver launch system.

Key difference from LIM: The LSM uses superconducting DC field coils on the sled
that interact directly with the track's AC stator. No eddy currents, no slip,
no thermal limit on the sled. The voltage limit is the primary constraint.

Design choices:
  - N_stator = 50: Fixed stator turn count. Optimization analysis shows the
    product F×v_cross = P/L is invariant with N and B_sled.
    N=50 gives only 2.3% time overhead vs the theoretical minimum (N→∞).
    Switchable windings provide zero benefit: cargo always wants max N,
    crewed missions are g-limited so N is irrelevant.
  - τ_p = 130 m: Large pole pitch minimizes HTS hysteresis losses, which scale
    as ~1/τ_p³ (B_coil_peak ~ N×I/l_coil where l_coil = τ_p/3).
  - The LSM requires periodic alternating N/S SC poles on the sled (not uniform
    DC field). The stator AC field locks to the sled pole pattern synchronously.

Launch modes:
  - "crew": v_target = 15,000 m/s (Mars transfer), 3g human physiological limit
  - "cargo": v_target = 30,000 m/s (deep solar system), 13g structural limit

The stator (round 13, 2026-11-17).  Two stators are held here, chosen by one switch:

  ORBITAL_RING_LSM_STATOR=ruled     (the default) the stator ruled on 17 November 2026 and written into Chapters 6, 9, 10,
                                    14 and 19: 6 m pole pitch, 17 turns of 30 soldered layers of 3 mm tape standing on
                                    edge, 2,088 A, 14,190 coils facing the 10 km sled, two sides; above 5 km/s the power is
                                    held by lowering the sled's field as (5 km/s / v)^0.4 and the coil current as
                                    (5 km/s / v)^0.6; 60 degree load angle.  Section 6b below.
  ORBITAL_RING_LSM_STATOR=table96   the 130 m stator of the earlier Table 9.6 (50 turns, five layers of 12 mm tape,
                                    2,107 A, Equations 10.1 and 10.2 as they were coded, the back-EMF held at 100 kV above
                                    5 km/s by lowering the sled's field alone).  Kept for the record: with it
                                    lsm_simulation.py prints what it printed before round 13, line for line.

Set the variable before Python starts (bash: ORBITAL_RING_LSM_STATOR=table96 python lsm_simulation.py; PowerShell:
$env:ORBITAL_RING_LSM_STATOR = 'table96'; python lsm_simulation.py), or pass --stator=table96 to lsm_simulation.py.
Both stators give the same launch: 83.29 MN, 416.47 GW, 189.7 minutes, 222,015.5 km, 4.99 laps.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import json
import math
import ring_altitude as ring
import os

STATOR = os.environ.get("ORBITAL_RING_LSM_STATOR", "ruled").strip().lower()
if STATOR not in ("ruled", "table96"):
    raise ValueError(f"ORBITAL_RING_LSM_STATOR must be 'ruled' or 'table96', not {STATOR!r}")

# =============================================================================
# SECTION 1: ORBITAL RING PARAMETERS
# =============================================================================

# Orbital parameters at the design altitude (ring_altitude.py)
R_ORBIT = ring.R_ORBIT      # m, orbital radius
G_250 = ring.G_LOCAL        # m/s², gravitational acceleration at the ring (name kept for compatibility)
G_0 = 9.807                 # m/s², standard gravity at sea level
V_ORBIT = ring.V_ORBIT      # m/s, orbital velocity at the ring
L_RING = ring.L_RING        # m, ring circumference

# HVDC power grid capacity
# The ring's total orbit-averaged generation is 8 MW at each LIM site (666 GW at 250 km, 722 GW at 800 km).
# The actual power available for the mass driver depends on how much of the HVDC grid
# has been upgraded for mass driver service during post-deployment buildout.
# This is the maximum instantaneous power that can be delivered to the active stator zone.
# Set to None to disable the power limit (for parameter studies).
P_HVDC_MAX = round(L_RING / 500.0) * 8e6  # W, maximum power delivery to mass driver (8 MW at each LIM site, 500 m spacing)

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
N_STATOR = 50               # turns per stator coil (fixed; see lsm_n_optimization_demo.py)

# HTS tape specification (Gömöry model)
# 12 mm x 5 layers: simpler winding than LIM's 3 mm x 16 layers.
# 5 layers gives better operating margin than 4.
ALPHA_TAPE_FIELD = 20.0                      # degrees — Gömöry perpendicular field angle
SIN_ALPHA = math.sin(math.radians(ALPHA_TAPE_FIELD))  # 0.34202
DE_RATING_FACTOR = 1 - SIN_ALPHA             # 0.658 — Ic reduction per layer
TAPE_WIDTH_MM = 12.0                         # mm — 12 mm REBCO tape
TAPE_WIDTH = TAPE_WIDTH_MM / 1000.0          # m
I_C_PER_MM_LAYER = 66.7                      # A/mm-width/layer (critical current density)
N_HTS_LAYERS = 5                             # layers of 12 mm tape
I_C_PER_LAYER_NOMINAL = I_C_PER_MM_LAYER * TAPE_WIDTH_MM   # 800 A nominal
I_C_PER_LAYER = I_C_PER_LAYER_NOMINAL * DE_RATING_FACTOR   # 526.4 A effective
I_CRITICAL = I_C_PER_LAYER * N_HTS_LAYERS    # 2633 A total
I_TARGET = 0.70 * I_CRITICAL                 # 1843 A normal operating current
I_PEAK = 0.80 * I_CRITICAL                   # 2107 A peak current (thrust-limited cases)
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

A_MAX_G = 1.0               # maximum tangential acceleration in g (system is thrust-limited at ~0.75g)
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
# SECTION 6b: THE RULED STATOR (the default since round 13; STATOR == "ruled")
#
# Sections 2 and 3 above hold the 130 m stator of the earlier Table 9.6.  When
# the switch is "ruled" the values below replace them.  The tape current, the
# coil's inductance and the loss table come from the tape-by-tape model of
# rounds 11 and 12, run once by tools/round13_stator_table.py and kept in
# data/lsm_stator_ruled.json.
# =============================================================================

K_W = 1.0                   # winding factor in the back-EMF and the current sheet (1 for the old stator's formulas)
ALPHA_SLED = 1.0            # above the corner speed the sled's field goes as s^ALPHA_SLED, the coil current as s^(1 - ALPHA_SLED)
V_CORNER = None             # m/s, where the power limit begins (the old stator: where the back-EMF reaches the voltage limit)
ACTIVE_FRACTION = None      # share of the sled's length that faces stator (the old stator: from its repeating units)
COILS_Q = None              # coils per phase per pole pitch on each side (the old stator's formulas do not use it)
L_COIL = None               # H, the coil's cyclic inductance (self and mutual)
F_PEAK_DESIGN = None        # N
STATOR_TABLE = None         # the loss table of the ruled stator

if STATOR == "ruled":
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "lsm_stator_ruled.json"), encoding="utf-8") as _f:
        STATOR_TABLE = json.load(_f)
    TAU_P = 6.0                                  # m, pole pitch
    N_POLES_PER_UNIT = 1                         # the sled's poles are counted one by one
    L_GAP = 0.0
    N_STATOR = 17                                # turns per stator coil
    TAPE_WIDTH_MM = 3.0                          # mm, tape without copper, standing on edge
    TAPE_WIDTH = TAPE_WIDTH_MM / 1000.0
    N_HTS_LAYERS = 30                            # soldered layers in a turn
    K_W = 3.0 / math.pi                          # coil sides spread evenly over their 60 degree belts
    ALPHA_SLED = 0.4
    V_CORNER = 5000.0                            # m/s
    ACTIVE_FRACTION = 0.975                      # 9,750 m of the 10 km sled, 39,000 m2 on two sides
    F_PEAK_DESIGN = 83_294_324.34775715          # N, the design launch's peak thrust (the book's 83.29 MN, to the last digit the old equations gave)
    T_STATOR = 70.0                              # K
    assert (STATOR_TABLE["tau_p"], STATOR_TABLE["n_turns"], STATOR_TABLE["alpha"], STATOR_TABLE["layers"]) == (TAU_P, N_STATOR, ALPHA_SLED, N_HTS_LAYERS)
    I_C_PER_LAYER = 2.0 * STATOR_TABLE["tape_current"]       # A, the worst tape's measured limit in the field it sits in (the tape runs at half of it)
    I_CRITICAL = I_C_PER_LAYER * N_HTS_LAYERS
    I_PEAK = STATOR_TABLE["tape_current"] * N_HTS_LAYERS     # A, 2,087.6: the 50 percent rule
    I_TARGET = I_PEAK
    L_COIL = STATOR_TABLE["L_coil_cyclic"]                   # H, 3.365 mH
    # the coils that give the design thrust: F = sides x 3 q (L_active / tau) x N k_w I B w sin(delta)
    COILS_Q = F_PEAK_DESIGN / (N_LIM_SIDES * 3.0 * (ACTIVE_FRACTION * 10_000.0 / TAU_P) * N_STATOR * K_W * I_PEAK * B_SLED_NOMINAL * W_COIL * math.sin(DELTA_MAX))
    P_STATOR_MAX = F_PEAK_DESIGN * V_CORNER                  # W, the rating of the stator and its supply (416.47 GW)

# =============================================================================
# SECTION 7: DERIVED PARAMETERS
# =============================================================================

def calc_derived():
    """Calculate all derived quantities from base parameters."""
    global L_UNIT, F_ACTIVE, N_UNITS, A_ACTIVE_TOTAL, W_ACTIVE_PER_UNIT
    global M_SLED_HARDWARE, M_TOTAL, KE_TARGET, F_MAX_SUPPLY
    global A_STATOR_COIL

    if STATOR == "ruled":
        # the sled's poles, one pole pitch each, over the active share of its length
        N_UNITS = int(round(L_SLED * ACTIVE_FRACTION / TAU_P))
        W_ACTIVE_PER_UNIT = N_POLES_PER_UNIT * TAU_P
        L_UNIT = L_SLED / N_UNITS
        F_ACTIVE = W_ACTIVE_PER_UNIT / L_UNIT
    else:
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
        # round 13: the stator switch and the ruled stator's values (None or 1 for the old stator)
        'stator': STATOR,
        'winding_factor': K_W,
        'coils_per_phase_per_pole_pitch': COILS_Q,
        'alpha_sled': ALPHA_SLED,
        'v_corner': V_CORNER,
        'l_coil': L_COIL,
        'stator_table': STATOR_TABLE,
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
GRAPH_OUTPUT_DIR = os.path.join(ring.OUT_DIR, "graphs_lsm")   # under ORBITAL_RING_OUT (default runs)
GRAPH_DPI = 300
GRAPH_WIDTH_INCHES = 10
GRAPH_HEIGHT_INCHES = 6
GRAPH_FORMAT = "png"
