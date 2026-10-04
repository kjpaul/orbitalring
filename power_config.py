#!/usr/bin/env python3
"""
Power Configuration Module - Solar generation and HVDC distribution parameters

This module contains all configurable parameters for the orbital ring
power generation and circumferential HVDC distribution simulation.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
import ring_altitude as ring

# =============================================================================
# SECTION 1: USER-CONFIGURABLE PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Solar Panel Configuration
# -----------------------------------------------------------------------------
PANEL_WIDTH = 53.0              # Total panel width across casing (m)
CELL_EFFICIENCY = 0.45          # Multi-junction cell efficiency (45%)
PANEL_PACKING = 1.0             # Panel area packing factor — Ch7 works in
                                # "per square meter of solar array" which
                                # already accounts for cell layout. Leave at
                                # 1.0 to match the published 300/844 numbers.
AR_BOOST = 1.0785               # Anti-reflection nanostructure boost factor.
                                # Standard AR capture ~89% of geometric limit;
                                # moth-eye / inverted-pyramid nanostructures
                                # push this to ~96%. Net ~7.85% multiplier.
                                # Fit to match Ch7 peak of 844 W/m^2 exactly.

# -----------------------------------------------------------------------------
# 1.2 Solar Environment
# -----------------------------------------------------------------------------
SOLAR_CONSTANT = 1361.0         # W/m^2 at 1 AU
EARTH_ALBEDO = 0.30             # Average Earth albedo

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
ALTITUDE = ring.ARRAY_ALTITUDE  # Altitude of the solar array (m): ring altitude minus the array offset (ring_altitude.py)
R_ORBIT = R_EARTH + ALTITUDE    # Radius of the array (m); sets shadow and albedo geometry
L_ARRAY = ring.L_ARRAY          # Circumference of the array (m)
L_RING = ring.L_RING            # Ring circumference (m)


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
# SHADOW_HALF_ANGLE is the half-angle of the SUNLIT arc measured from
# solar noon, so at 250 km it is about 105.8 degrees. The albedo model
# uses the (1+cos(phi))/2 phase formula, which extends smoothly through
# the shadow boundary (as the physics actually does), so direct flux is
# the only term that uses this cutoff.
SHADOW_HALF_ANGLE = math.pi - math.asin(R_EARTH / R_ORBIT)  # rad from subsolar

# Earth-view factor from a flat plate above Earth at the orbital radius.
# This is the fraction of the downward hemisphere filled by Earth, used
# to scale the peak Earth-albedo flux on the back face of a panel.
VIEW_FACTOR_EARTH = (R_EARTH / R_ORBIT) ** 2

# Peak Earth-albedo flux on the back face of a panel at solar noon.
# This is the "378 W/m^2" value in Chapter 7.
ALBEDO_PEAK = EARTH_ALBEDO * SOLAR_CONSTANT * VIEW_FACTOR_EARTH  # W/m^2

# Total power demand
P_DEMAND_DEPLOYMENT = LIM_SITES * LIM_POWER_PER_SITE       # W
P_DEMAND_OPS = LIM_SITES * OPS_POWER_PER_SITE              # W

# Reference average power output (Chapter 7 cross-check)
# 1361 * 0.45 * 0.4572 * 1.0785 = 301.99 W/m^2, text rounds to 300.
P_AVG_REFERENCE = 302.0        # W/m^2, orbit-averaged bifacial output

# Peak power output at local noon (cross-check)
# (1361 + 378) * 0.45 * 1.0785 = 844.0 W/m^2
P_PEAK_REFERENCE = 844.0       # W/m^2, at subsolar point
