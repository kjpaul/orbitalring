#!/usr/bin/env python3
"""
Power Configuration Module - Solar generation and HVDC distribution parameters

This module contains all configurable parameters for the orbital ring
power generation and circumferential HVDC distribution simulation.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: USER-CONFIGURABLE PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Solar Panel Configuration
# -----------------------------------------------------------------------------
PANEL_WIDTH = 53.0              # Total panel width across casing (m)
CELL_EFFICIENCY = 0.45          # Multi-junction cell efficiency (45%)
PANEL_PACKING = 0.95            # Panel area packing factor (95% coverage)

# -----------------------------------------------------------------------------
# 1.2 Solar Environment
# -----------------------------------------------------------------------------
SOLAR_CONSTANT = 1361.0         # W/m^2 at 1 AU
EARTH_ALBEDO = 0.30             # Average Earth albedo
ALBEDO_FRACTION = 0.30          # Fraction of albedo reaching panels (view factor)

# -----------------------------------------------------------------------------
# 1.3 HVDC Transmission Configuration
# -----------------------------------------------------------------------------
V_HVDC = 10e6                   # Pole-to-pole voltage (V) — 10 MV bipolar
LOSS_BUDGET = 0.05              # Target fractional loss (5%)
CABLE_DIAMETER = 0.10           # Individual cable diameter (m)

# -----------------------------------------------------------------------------
# 1.4 Conductor Material — Graphene-Enhanced CNT Fiber
# -----------------------------------------------------------------------------
SIGMA_CONDUCTOR = 30e6          # Electrical conductivity (S/m)
RHO_CONDUCTOR = 1700.0          # Density (kg/m^3)

# -----------------------------------------------------------------------------
# 1.5 Power Demand Configuration
# -----------------------------------------------------------------------------
LIM_POWER_PER_SITE = 8e6        # LIM power per site during deployment (W)
OPS_POWER_PER_SITE = 500e3      # Post-deployment operations power per site (W)
LIM_SPACING = 500.0             # Distance between LIM sites (m)

# -----------------------------------------------------------------------------
# 1.6 Simulation Resolution
# -----------------------------------------------------------------------------
N_POINTS = 3600                 # Angular resolution (points around ring)

# -----------------------------------------------------------------------------
# 1.7 Parametric Sweep Ranges
# -----------------------------------------------------------------------------
PANEL_WIDTH_MIN = 53.0          # Minimum panel width for sweep (m)
PANEL_WIDTH_MAX = 400.0         # Maximum panel width for sweep (m)
PANEL_WIDTH_STEPS = 50          # Number of steps in panel width sweep

V_HVDC_MIN = 2e6                # Minimum HVDC voltage for sweep (V)
V_HVDC_MAX = 20e6               # Maximum HVDC voltage for sweep (V)
V_HVDC_STEPS = 50               # Number of steps in voltage sweep

# -----------------------------------------------------------------------------
# 1.8 Output Control
# -----------------------------------------------------------------------------
SAVE_GRAPHS = True
GRAPH_OUTPUT_DIR = "./graphs_power"
GRAPH_DPI = 300
GRAPH_WIDTH_INCHES = 10
GRAPH_HEIGHT_INCHES = 5
GRAPH_FORMAT = "png"


# =============================================================================
# SECTION 2: PHYSICAL CONSTANTS
# =============================================================================

R_EARTH = 6_371_000.0           # Earth mean radius (m) — WGS-84
ALTITUDE = 250_000.0            # Orbital altitude (m)
R_ORBIT = R_EARTH + ALTITUDE    # Orbital radius (m)
L_RING = 41_645_813.012         # Ring circumference (m) — canonical value from Chapter 3


# =============================================================================
# SECTION 3: DERIVED PARAMETERS
# =============================================================================

# LIM site count
LIM_SITES = round(L_RING / LIM_SPACING)

# HVDC conductor resistivity
RHO_ELEC = 1.0 / SIGMA_CONDUCTOR  # Electrical resistivity (ohm-m)

# Per-pole voltage (bipolar system)
V_POLE = V_HVDC / 2.0

# Specific conductivity (figure of merit)
SPECIFIC_CONDUCTIVITY = SIGMA_CONDUCTOR / RHO_CONDUCTOR  # S*m^2/kg

# Shadow geometry
# At 250 km altitude, shadow starts at phi = 180 - arcsin(R_E / r_orbit)
SHADOW_HALF_ANGLE = math.pi - math.asin(R_EARTH / R_ORBIT)  # rad from subsolar

# Total power demand
P_DEMAND_DEPLOYMENT = LIM_SITES * LIM_POWER_PER_SITE       # W
P_DEMAND_OPS = LIM_SITES * OPS_POWER_PER_SITE              # W

# Reference average power output (Chapter 7 cross-check)
P_AVG_REFERENCE = 300.0        # W/m^2, orbit-averaged bifacial output

# Peak power output at local noon (cross-check)
P_PEAK_REFERENCE = 844.0       # W/m^2, at subsolar point
