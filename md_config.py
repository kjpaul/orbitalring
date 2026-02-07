"""
Mass Driver Configuration Module
Adapted for high-mass, high-velocity sled launch.
"""
import math

# =============================================================================
# 1. PAYLOAD & SLED CONFIGURATION
# =============================================================================
# From mass_driver_section_v2.md: 1g operational max is 8,372 tonnes
SLED_MASS_TOTAL = 8_372_000  # kg (Sled + Payload)
SLED_LENGTH = 1000.0         # meters (Active length of reaction plates)
SLED_WIDTH = 2.5             # meters (Physical width)

# Reaction plate (Secondary)
# Gamma-TiAl has higher resistivity than Al, but high thermal limit.
# We compensate with higher Stator B-fields.
PLATE_MATERIAL = "gamma_titanium" 
PLATE_WIDTH = 1.5            # Active magnetic width (m)
PLATE_THICKNESS = 0.10       # 50mm thick plates for heat sinking

# =============================================================================
# 2. STATOR (LIM) CONFIGURATION
# =============================================================================
# To push 8000+ tonnes, we need massive B-fields (0.5 - 1.0 T in the gap).
# We assume HTS Roebel cables capability.

N_TURNS = 100            # Turns per coil (High current, lower turns for low inductance)
I_PEAK = 5000.0          # Peak Amps per turn (Requires actively cooled HTS cable)
TAU_P = 50.0             # Pole Pitch (m). Shorter than ring to maintain frequency at low V.
W_COIL = 2.5             # Coil width (m)
GAP = 0.15               # Air gap (m) - Tighter control than the ring cable

# Derived Stator Properties
LIM_EFFICIENCY = 0.90    # Assumes regenerative capture efficiency
POWER_MAX_GW = 20.0      # Max power draw limit per sector (Gigawatts)

# =============================================================================
# 3. LAUNCH PARAMETERS
# =============================================================================
# Launch targets from markdown
G_FORCE_LIMIT = 29.43    # 3g limit for human-rated (or 98.1 for 10g cargo)
V_LAUNCH_TARGET = 11000  # m/s (Interplanetary injection velocity)
DT = 0.1                 # Simulation time step (s)

# Control Logic
# We need high slip at low speeds to penetrate the resistive TiAl plate
V_SLIP_MIN = 20.0        # Minimum slip velocity (m/s)
V_SLIP_MAX = 300.0       # Max slip to prevent excessive frequency loss
TARGET_SLIP_RATIO = 0.15 # Higher slip ratio target for resistive plates