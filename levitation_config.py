#!/usr/bin/env python3
"""
Levitation System Configuration Module

Parameters for the orbital ring levitation bearing simulation,
including EML attractive bearings, Halbach array vs solenoid comparison,
and passive superconducting Halbach safety backstop.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: ORBITAL RING GEOMETRY
# =============================================================================

R_ORBIT = 6_628_137              # Orbital radius at 250 km (m)
L_RING = 41_645_813.012          # Ring circumference (m)
G_NET = 9.038                    # Net downward acceleration at 250 km (m/s^2)
                                 # = gravitational - centrifugal at ground-sync speed
M_LOAD_M = 12_000               # Casing mass per meter (kg/m)
M_CABLE_STRUCTURAL = 96_700      # Structural CNT cable mass per meter (kg/m)

# =============================================================================
# SECTION 2: EML LEVITATION BEARING (Attractive)
# =============================================================================

# --- Halbach array configuration (baseline) ---
EML_N_TRACKS = 2                 # Number of levitation tracks (one each side of cable underside)
EML_COIL_WIDTH = 1.0             # Width of each HTS Halbach array (m)
EML_PLATE_WIDTH = 1.4            # Ferromagnetic plate width (m), 0.4 m wider than coil
EML_NOMINAL_GAP = 0.100          # Nominal operating gap (m)
EML_GAP_MIN = 0.020              # Minimum gap before safety engages (m)
EML_GAP_MAX = 0.180              # Maximum gap before safety engages (m)

# HTS tape for levitation coils
EML_HTS_TAPE_WIDTH_MM = 12      # Theva TPL5000 tape width (mm)
EML_HTS_IC = 800                 # Critical current per tape at 70 K (A)
EML_I_OPERATING = 400            # Normal operating current (A) = Ic/2
EML_COILS_PER_METER = 40         # Coils per meter of array length
EML_TURNS_PER_COIL = 52          # Turns per double pancake coil
EML_COIL_SIZE = 0.025            # Coil cross-section dimension (m), 25 mm x 25 mm

# Ferromagnetic plate properties (sintered iron in resin)
FERRO_B_SAT = 1.2                # Saturation field of sintered iron (T)
FERRO_DENSITY = 7870             # Iron density (kg/m^3)
FERRO_PLATE_THICKNESS = 0.10     # Plate thickness (m)
FERRO_MU_R = 200                 # Relative permeability (effective, sintered in resin)
                                 # Lower than bulk iron (~5000) due to resin gaps
FERRO_HYSTERESIS_COEFF = 50      # Hysteresis loss coefficient (J/m^3 per cycle)
FERRO_RESISTIVITY = 1e-4         # Effective resistivity (ohm-m), high due to resin isolation

# --- Solenoid configuration (comparison) ---
SOL_N_TRACKS = 2
SOL_COIL_WIDTH = 1.0
SOL_TURNS_PER_COIL = 52
SOL_COILS_PER_METER = 40
SOL_I_OPERATING = 400
SOL_IC = 800
SOL_NOMINAL_GAP = 0.100
# Solenoid uses same HTS tape but arranged as simple solenoid coils
# rather than Halbach phased arrays

# =============================================================================
# SECTION 3: SAFETY BACKSTOP (Passive Repulsive Halbach Arrays)
# =============================================================================
# Superconducting Halbach arrays on the casing, top and bottom.
# These create repulsive eddy currents in a reaction plate on the cable.
# They engage only if the cable enters the "forbidden zone" beyond the
# nominal gap limits.

SAFETY_N_ARRAYS = 4              # 2 above cable, 2 below cable
SAFETY_COIL_WIDTH = 1.0          # Width of each safety Halbach array (m)
                                 # Wider arrays = more force area
SAFETY_HTS_IC = 800              # Critical current of safety HTS tape (A)
SAFETY_I_PERSISTENT = 750        # Persistent current in safety coils (A)
                                 # Set to ~94% of Ic. These coils are in persistent
                                 # mode (superconducting loop, no power supply needed).
                                 # Higher than EML operating current because the safety
                                 # MUST work when everything else has failed.
SAFETY_WAVELENGTH = 0.10         # Halbach array wavelength (m)
                                 # Short wavelength concentrates field on front face
                                 # and attenuates rapidly with distance.
                                 # decay length = lambda/(2*pi) = 15.9 mm
SAFETY_TURNS_PER_COIL = 24      # Turns per safety coil
                                 # 4 coils per wavelength, 25 mm coil pitch at 100mm lambda
                                 # Each coil is 25 mm wide x 50 mm deep,
                                 # wound with 3 mm HTS tape: 8 turns/layer x 3 layers
SAFETY_COIL_THICKNESS = 0.050    # Thickness of Halbach array assembly (m)
SAFETY_STANDOFF = 0.020          # Distance from casing wall to safety array face (m)

# The safety reaction plate is a separate plate on the cable, positioned
# at the top and bottom (not the same as the EML ferromagnetic plate).
# This plate must be electrically conductive for eddy current repulsion.

# Forbidden zone boundaries (measured from cable center to casing wall)
FORBIDDEN_ZONE_INNER = 0.020     # If gap < 20 mm, safety pushes cable away
FORBIDDEN_ZONE_OUTER = 0.180     # If gap > 180 mm, upper safety pushes cable back

# =============================================================================
# SECTION 4: SAFETY REACTION PLATE MATERIAL CANDIDATES
# =============================================================================
# These are the materials being evaluated for the safety backstop
# reaction plates. The key properties are:
#   - Electrical conductivity (higher = stronger eddy currents = more repulsion)
#   - Thermal conductivity (higher = better heat spreading)
#   - Melting/service temperature (higher = more margin before failure)
#   - Density (lower = less cable mass penalty)

REACTION_PLATE_MATERIALS = {
    "aluminum_6061": {
        "name": "Aluminum 6061-T6",
        "resistivity": 3.99e-8,       # ohm-m at 293 K
        "temp_coeff_rho": 0.0039,     # /K resistivity temperature coefficient
        "density": 2700,              # kg/m^3
        "specific_heat": 896,         # J/(kg K)
        "thermal_conductivity": 167,  # W/(m K)
        "melting_point": 855,         # K (solidus)
        "max_service_temp": 573,      # K (above this, strength degrades severely)
        "emissivity": 0.15,           # Radiation emissivity (oxidized Al)
        "young_modulus": 69e9,        # Pa
        "notes": "Excellent conductivity, low melting point. Standard aerospace alloy."
    },
    "copper_OFHC": {
        "name": "OFHC Copper",
        "resistivity": 1.68e-8,       # ohm-m at 293 K
        "temp_coeff_rho": 0.00393,    # /K
        "density": 8960,              # kg/m^3
        "specific_heat": 385,         # J/(kg K)
        "thermal_conductivity": 401,  # W/(m K)
        "melting_point": 1358,        # K
        "max_service_temp": 773,      # K (softening begins)
        "emissivity": 0.78,           # Radiation emissivity (oxidized Cu)
        "young_modulus": 117e9,       # Pa
        "notes": "Best conductivity, highest eddy current force. Heavy."
    },
    "gamma_TiAl": {
        "name": "gamma-TiAl (Ti-48Al-2Cr-2Nb)",
        "resistivity": 75e-8,         # ohm-m at 293 K
        "temp_coeff_rho": 0.001,      # /K (estimated)
        "density": 3900,              # kg/m^3
        "specific_heat": 600,         # J/(kg K)
        "thermal_conductivity": 22,   # W/(m K)
        "melting_point": 1733,        # K
        "max_service_temp": 1073,     # K (750 C service temperature)
        "emissivity": 0.45,           # Radiation emissivity (oxidized)
        "young_modulus": 176e9,       # Pa
        "notes": "Already used for LIM reaction plates. High temp, moderate conductivity."
    },
    "molybdenum": {
        "name": "Molybdenum",
        "resistivity": 5.34e-8,       # ohm-m at 293 K
        "temp_coeff_rho": 0.0047,     # /K
        "density": 10220,             # kg/m^3
        "specific_heat": 251,         # J/(kg K)
        "thermal_conductivity": 138,  # W/(m K)
        "melting_point": 2896,        # K
        "max_service_temp": 1923,     # K
        "emissivity": 0.30,           # Radiation emissivity (oxidized)
        "young_modulus": 329e9,       # Pa
        "notes": "Extremely high melting point. Good conductor. Very heavy."
    },
    "graphene_CNT_composite": {
        "name": "Graphene/CNT Composite",
        "resistivity": 3.33e-8,       # ohm-m (30 MS/m target, same as HVDC conductor)
        "temp_coeff_rho": -0.0005,    # /K (negative: metallic CNT bundles)
        "density": 1700,              # kg/m^3
        "specific_heat": 700,         # J/(kg K)
        "thermal_conductivity": 200,  # W/(m K) (axial much higher, transverse lower)
        "melting_point": 3800,        # K (graphitization, not true melting)
        "max_service_temp": 2500,     # K (in vacuum, no oxidation)
        "emissivity": 0.85,           # Radiation emissivity (CNT forest)
        "young_modulus": 200e9,       # Pa (composite, not single tube)
        "notes": "Lightest option. Highest temp tolerance in vacuum. Speculative conductivity."
    },
}

# Default reaction plate for safety backstop
SAFETY_PLATE_MATERIAL = "copper_OFHC"
SAFETY_PLATE_THICKNESS = 0.020    # m (thinner than LIM plate, only needs to survive transient)
SAFETY_PLATE_WIDTH = 1.0          # m (matches safety array width)

# =============================================================================
# SECTION 5: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7         # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.670374e-8    # Stefan-Boltzmann constant (W/(m^2 K^4))
T_SPACE = 2.7                     # Cosmic microwave background (K)
T_CASING_WALL = 247               # Casing wall temperature in steady state (K)
T_CABLE_STEADY = 254              # Cable steady-state temperature (K)

# =============================================================================
# SECTION 6: SIMULATION PARAMETERS
# =============================================================================

GAP_SWEEP_MIN = 0.005             # Minimum gap for force-vs-gap sweep (m)
GAP_SWEEP_MAX = 0.300             # Maximum gap for force-vs-gap sweep (m)
GAP_SWEEP_POINTS = 200            # Number of points in gap sweep

TEMP_SWEEP_MIN = 250              # Minimum temperature for thermal sweep (K)
TEMP_SWEEP_MAX = 3000             # Maximum temperature for thermal sweep (K)
TEMP_SWEEP_POINTS = 200           # Number of points in thermal sweep

# Transient simulation for safety engagement
TRANSIENT_DT = 1e-4               # Time step for transient sim (s)
TRANSIENT_DURATION = 5.0          # Total transient simulation time (s)
