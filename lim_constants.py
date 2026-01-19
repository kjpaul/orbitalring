import math

# =============================================================================
# SECTION 2: PHYSICAL CONSTANTS
# =============================================================================
# These values should not be changed unless you have good reason.

# Fundamental constants
MU0 = 4 * math.pi * 1e-7        # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.670374e-8  # Stefan-Boltzmann constant (W/m²K⁴)

# Orbital parameters at 250 km altitude
V_ORBIT = 7754.866              # Orbital velocity (m/s)
V_GROUND_STATIONARY = 483.331   # Ground-stationary velocity at 250 km (m/s)
L_RING = 41_645_813.012         # Ring circumference (m)

# Liquid nitrogen properties
T_LN2_BOIL = 77.4               # Boiling point at 1 atm (K)
T_LN2_SUPPLY = 70               # Supply temperature from cryo system (K)
C_P_LN2 = 2040                  # Specific heat capacity (J/kg·K)
L_V_LN2 = 199000                # Latent heat of vaporization (J/kg)
H_CONV_LN2 = 90                 # Convective heat transfer coefficient (W/m²·K)

# Thermal environment
Q_SUN = 1361                    # Solar flux at 1 AU (W/m²)
Q_EARTH_ALBEDO = 650            # Earth albedo contribution (W/m²)
Q_SHIELDING = 0.005             # MLI shielding effectiveness (fraction transmitted)
T_SPACE = 2.7                   # Deep space temperature (K)
T_RADIATOR_HOT = 300            # Cryo radiator hot side temperature (K)
T_MAX_PLATE = 500               # Maximum reaction plate temperature (K)

# Cryogenic system
CRYO_EFF = 0.18                 # Cryo system efficiency (fraction of Carnot)
EM_HEATSINK = 0.9               # Heatsink emissivity

# HTS tape properties
HTS_THICKNESS_UM = 80           # Tape thickness in micrometers
ALPHA_PENETRATION_DEG = 20.0    # Magnetic field penetration angle (degrees)


# =============================================================================
# SECTION 3: DERIVED PARAMETERS
# =============================================================================
# Calculated from user parameters and constants. Do not edit directly.

# Time constants
HR = 3600
DAY = 24 * HR
WEEK = 7 * DAY
MONTH = 30 * DAY
YR = round(365.33 * DAY)