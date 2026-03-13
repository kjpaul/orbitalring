#!/usr/bin/env python3
"""
Space Elevator Tether Configuration
====================================

All physical constants and configurable parameters for the space elevator
tether dynamics simulation.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# 1. PHYSICAL CONSTANTS
# =============================================================================

GM = 3.986004418e14             # Earth gravitational parameter (m^3/s^2)
R_E_EQUATOR = 6_378_137.0      # Earth equatorial radius (m)
R_E_POLES = 6_356_752.0         # Earth polar radius (m)
R_E = R_E_EQUATOR               # Default Earth radius (m)
R_ORBIT_RING = 6_628_137.0     # Orbital ring radius (m), 250 km altitude
R_GEO = 42_164_000.0           # Geosynchronous orbit radius (m)
OMEGA_SIDEREAL = 7.2921159e-5  # Earth sidereal rotation rate (rad/s)
G_SURFACE = 9.80665            # Standard gravity (m/s^2)

# Verify GEO consistency: GM/R_GEO^2 should equal OMEGA^2 * R_GEO
_a_grav_geo = GM / R_GEO**2
_a_cent_geo = OMEGA_SIDEREAL**2 * R_GEO
# These differ by ~0.1% due to rounding in R_GEO; acceptable.

# =============================================================================
# 2. MATERIAL PROPERTIES — CNT TETHER
# =============================================================================

SIGMA_BREAK_CNT = 25.0e9       # CNT breaking stress (Pa)
RHO_CNT = 1_700.0              # CNT density (kg/m^3)
SAFETY_FACTOR = 2.0            # Safety factor (dimensionless)
SIGMA_OPERATING = SIGMA_BREAK_CNT / SAFETY_FACTOR  # Operating stress (Pa)

# Characteristic velocity for CNT
V_CHAR_CNT = math.sqrt(SIGMA_OPERATING / RHO_CNT)  # ~2711 m/s

# =============================================================================
# 3. TETHER CONFIGURATIONS
# =============================================================================

# Configuration 1: GEO to ground
R_BOTTOM_GROUND = R_E_EQUATOR   # Bottom of ground-anchored tether (m)
ALT_BOTTOM_GROUND = 0.0         # Altitude (m)

# Configuration 2: GEO to 275 km
ALT_BOTTOM_SHORT = 275_000.0    # Short tether termination altitude (m)
R_BOTTOM_SHORT = R_E + ALT_BOTTOM_SHORT  # 6,653,137 m

# GEO altitude
ALT_GEO = R_GEO - R_E          # ~35,786 km

# =============================================================================
# 4. POD / PAYLOAD PARAMETERS
# =============================================================================

M_POD_DEFAULT = 10_000.0       # Default pod mass (kg) = 10 tonnes
M_PAYLOAD_DEFAULT = 5_000.0    # Default payload mass (kg) = 5 tonnes

# Freefall pod parameters
C_D_POD = 1.2                  # Drag coefficient (blunt body)
BALLISTIC_COEFF = 500.0        # Default ballistic coefficient (kg/m^2)
R_NOSE = 1.0                   # Nose radius for heat flux calc (m)

# Sutton-Graves constant for air
K_SUTTON_GRAVES = 1.7415e-4    # W/(m^1.5 * (m/s)^3 * (kg/m^3)^0.5)

# =============================================================================
# 5. CLIMBER PARAMETERS
# =============================================================================

V_DESCENT_DEFAULT = 55.56      # Default descent speed (m/s) = 200 km/h
V_DESCENT_WEEK = 59.15         # ~35786 km / 168 hr ≈ 59.15 m/s for 1 week

# =============================================================================
# 6. US STANDARD ATMOSPHERE 1976 — LAYER DEFINITIONS
# =============================================================================

# Layer boundaries: (h_base [m], T_base [K], lapse_rate [K/m], p_base [Pa])
# Pressure at each layer base is computed from the layer below.
_M_AIR = 0.0289644             # Molar mass of air (kg/mol)
_R_GAS = 8.31447               # Universal gas constant (J/(mol·K))

# Layer table: (h_base_m, T_base_K, lapse_rate_K_per_m)
ATMO_LAYERS = [
    (0,      288.150,  -0.0065),    # Troposphere
    (11000,  216.650,   0.0),       # Tropopause
    (20000,  216.650,   0.001),     # Stratosphere
    (32000,  228.650,   0.0028),    # Stratosphere upper
    (47000,  270.650,   0.0),       # Stratopause
    (51000,  270.650,  -0.0028),    # Mesosphere
    (71000,  214.650,  -0.002),     # Mesosphere upper
]

# Pre-compute pressure at each layer base
_P_SEA_LEVEL = 101325.0  # Pa
_ATMO_P_BASE = [_P_SEA_LEVEL]
for i in range(len(ATMO_LAYERS) - 1):
    h_b, T_b, L = ATMO_LAYERS[i]
    h_next = ATMO_LAYERS[i + 1][0]
    dh = h_next - h_b
    if abs(L) < 1e-10:
        # Isothermal
        p_next = _ATMO_P_BASE[i] * math.exp(
            -G_SURFACE * _M_AIR * dh / (_R_GAS * T_b)
        )
    else:
        # Gradient layer
        T_next = T_b + L * dh
        exponent = G_SURFACE * _M_AIR / (_R_GAS * L)
        p_next = _ATMO_P_BASE[i] * (T_next / T_b) ** (-exponent)
    _ATMO_P_BASE.append(p_next)


def atmosphere(h):
    """US Standard Atmosphere 1976 — temperature, pressure, density.

    Args:
        h: Geometric altitude (m), 0 to 86 km.
            Above 86 km returns exponential extrapolation.

    Returns:
        (T [K], p [Pa], rho [kg/m^3])
    """
    if h < 0:
        h = 0.0

    # Above 86 km: exponential extrapolation
    if h > 86000:
        # Use scale height model above 86 km
        T_86 = 186.87  # K at 86 km (USSA76)
        # Density at 86 km from the layered model
        T_71, p_71, _ = _atmo_layer(71000, 6)
        L = ATMO_LAYERS[6][2]
        dh = 86000 - 71000
        T_top = T_71 + L * dh
        exponent = G_SURFACE * _M_AIR / (_R_GAS * L)
        p_86 = _ATMO_P_BASE[6] * (T_top / T_71) ** (-exponent)
        rho_86 = p_86 * _M_AIR / (_R_GAS * T_top)

        # Exponential decay with scale height
        H_scale = _R_GAS * T_86 / (_M_AIR * G_SURFACE)  # ~5.5 km
        rho = rho_86 * math.exp(-(h - 86000) / H_scale)
        T = T_86
        p = rho * _R_GAS * T / _M_AIR
        return T, p, rho

    # Find the layer
    layer_idx = 0
    for i in range(len(ATMO_LAYERS) - 1, -1, -1):
        if h >= ATMO_LAYERS[i][0]:
            layer_idx = i
            break

    return _atmo_layer(h, layer_idx)


def _atmo_layer(h, layer_idx):
    """Compute atmosphere within a specific layer."""
    h_b, T_b, L = ATMO_LAYERS[layer_idx]
    p_b = _ATMO_P_BASE[layer_idx]
    dh = h - h_b

    if abs(L) < 1e-10:
        T = T_b
        p = p_b * math.exp(-G_SURFACE * _M_AIR * dh / (_R_GAS * T_b))
    else:
        T = T_b + L * dh
        exponent = G_SURFACE * _M_AIR / (_R_GAS * L)
        p = p_b * (T / T_b) ** (-exponent)

    rho = p * _M_AIR / (_R_GAS * T)
    return T, p, rho


def speed_of_sound(T):
    """Speed of sound in air at temperature T (K).

    Uses gamma = 1.4 for diatomic ideal gas.

    Returns:
        Speed of sound (m/s)
    """
    gamma = 1.4
    return math.sqrt(gamma * _R_GAS * T / _M_AIR)


# =============================================================================
# 7. WIND PROFILE MODEL
# =============================================================================

def wind_speed(h):
    """Simple wind speed profile vs altitude (m/s).

    Troposphere: increases to jet stream peak ~60 m/s at 10-12 km.
    Stratosphere: decreases to ~10 m/s.
    Mesosphere: variable, ~20 m/s.
    Above 100 km: negligible (tether moves through very thin air).

    Args:
        h: Altitude (m)

    Returns:
        Wind speed (m/s)
    """
    h_km = h / 1000.0
    if h_km < 1:
        return 10.0 * h_km          # Surface boundary layer
    elif h_km < 10:
        return 10.0 + 50.0 * (h_km - 1) / 9.0  # Rising to jet stream
    elif h_km < 13:
        return 60.0                  # Jet stream core
    elif h_km < 20:
        return 60.0 - 40.0 * (h_km - 13) / 7.0  # Declining
    elif h_km < 50:
        return 20.0 - 10.0 * (h_km - 20) / 30.0  # Stratosphere
    elif h_km < 80:
        return 10.0 + 10.0 * (h_km - 50) / 30.0  # Mesosphere
    elif h_km < 100:
        return 20.0 * (100 - h_km) / 20.0  # Declining to zero
    else:
        return 0.0


# =============================================================================
# 8. SIMULATION RESOLUTION
# =============================================================================

N_POINTS_TETHER = 10000        # Number of integration points along tether
N_POINTS_TRAJECTORY = 100000   # Max integration steps for trajectory
DT_TRAJECTORY = 0.1            # Time step for trajectory integration (s)

# =============================================================================
# 9. OUTPUT CONTROL
# =============================================================================

SAVE_GRAPHS = True
GRAPH_OUTPUT_DIR = "./graphs_tether"
GRAPH_DPI = 300
GRAPH_WIDTH_INCHES = 10
GRAPH_HEIGHT_INCHES = 6
GRAPH_FORMAT = "png"

# =============================================================================
# 10. INCLINED ORBIT PARAMETERS
# =============================================================================

INCLINATION_ANGLES_DEG = [0, 5, 10, 20, 30, 45]  # Degrees

# =============================================================================
# 11. PARAMETRIC SWEEP RANGES
# =============================================================================

BETA_SWEEP = [100, 200, 500, 1000, 2000]          # Ballistic coefficients (kg/m^2)
POD_MASS_SWEEP = [1000, 5000, 10000, 50000]        # Pod masses (kg)
DESCENT_SPEED_SWEEP = [20, 40, 60, 80, 100]        # Descent speeds (m/s)
SIGMA_SWEEP_GPA = [5, 10, 15, 20, 25, 30, 40, 50,
                   60, 80, 100, 130]                # Material strengths (GPa)
