"""
Mass Driver Configuration Module - Validated Engineering Constraints
"""
import math

# =============================================================================
# 1. HTS TAPE CONFIGURATION (THEVA TLP AP @ 70K)
# =============================================================================
HTS_TAPE_WIDTH_MM = 12.0     # Using 12mm for max current capacity
HTS_TAPE_LAYERS = 2          # 2 layers to boost amp-turns without insane width
IC_PER_MM_PER_LAYER = 66.7   # Theva TLP AP spec at 70K
DE_RATING_FACTOR = 1 - math.sin(math.radians(20)) # Anisotropy/Self-field penalty
NORRIS_HYSTERESIS = True     # Critical: 12mm tape has high AC losses

# Current Calculations
# I_c total = Width * Layers * Ic/mm * De-rating
# We apply the 80% safety margin HERE.
_eff_layers = float(HTS_TAPE_LAYERS)
if HTS_TAPE_LAYERS > 1:
    _eff_layers *= DE_RATING_FACTOR

I_C_TOTAL = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * _eff_layers
I_TARGET = 0.80 * I_C_TOTAL  # <--- The 80% Hard Limit
I_PEAK_MAX = I_TARGET

# =============================================================================
# 2. STATOR (LIM) GEOMETRY
# =============================================================================
# To get thrust with lower current, we need more turns and optimization
N_TURNS = 300            # Increased turns to compensate for lower Amps
TAU_P = 100.0            # Pole pitch (m)
W_COIL = 2.0             # Coil width (m)
GAP = 0.15               # Air gap (m)
L_RING = 41_645_813.012  # Ring Circumference

# =============================================================================
# 3. REACTION PLATE (Gamma-TiAl)
# =============================================================================
PLATE_MATERIAL = "gamma_titanium"
T_PLATE = 0.20           # 200mm thickness (thermal mass)
W_PLATE = 1.5            # Plate width (m)
SLED_MASS = 8372000      # 8,372 Tonnes (1g Limit Payload)
SLED_LENGTH = 1000.0     # Active length of sled

# =============================================================================
# 4. SYSTEM LIMITS
# =============================================================================
VOLTS_MAX = 100e3        # 100 kV Insulation limit
MAX_SITE_POWER = 8e6     # 8 MW per site limit
CRYO_COP_PENALTY = 50.0  # W_electric / W_thermal (20:1 to 95:1 range)

# =============================================================================
# 5. LAUNCH PROFILE
# =============================================================================
V_LAUNCH_TARGET = 11000  # m/s
DT = 1.0                 # Time step (s) - slower dynamics, larger step