#!/usr/bin/env python3
"""
Parameters and material properties of the LIM deployment simulation.
Book: The Orbital Ring, Volume IV: Chapters 11 and 12 (imported by lim_simulation.py).
Run: imported; not run by itself
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

LIM Configuration Module - Parameters and material properties

This module contains all configurable parameters for the orbital ring
deployment simulation.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
import ring_altitude as ring
import os

# =============================================================================
# SECTION 1: USER-CONFIGURABLE PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 HTS (High-Temperature Superconductor) Tape Configuration
# -----------------------------------------------------------------------------
HTS_TAPE_WIDTH_MM = 3       # Standard widths: 12, 6, 4, or 3 mm
HTS_TAPE_LAYERS_LEGACY = 1  # Legacy rule: layers of tape in parallel: 1, 2, or 3 (the measured basis below sets its own count; HTS_TAPE_LAYERS is derived)
IC_PER_MM_PER_LAYER = 66.7  # Legacy rule: critical current density ~66.7 A/mm per layer, with no field around the tape
DE_RATING_FACTOR = 1 - math.sin(math.radians(20)) # Legacy rule: compensates for self-inductance in multilayer HTS configurations.
NORRIS_HYSTERESIS = False   # Use Norris Formula for hysteresis loss

# The coil current limit on the measured tape limit (2026-10-26, tools/lim_tape_limit.py, report
# Vol_III/runs_700km_hw10686_2026-10-12/lim_tape_limit_report_2026-10-26.md).  The current a tape carries falls with the field
# around it and with the field's angle to the tape, and the winding makes its own field; the limit here is the measured critical
# current of the worst-placed tape at its own field and angle at 70 K (Robinson dataset, data/THEVA_Robinson_70K_angle_dependence.csv,
# Pantoja and Strickland 2024, doi 10.6084/m9.figshare.25893079.v2), and the winding runs that tape at 50 percent, the bearing rule.
# Winding: the pancake of this file's own turn budget (COIL_WINDING_WIDTH: 1 mm Kapton, 0.06 mm adhesive, 0.12 mm per layer per turn),
# 208 turns on edge, 3 layers of 3 mm tape soldered per turn, bundle 3 x 295 mm.  Worst tape: at a stack end, where the field of
# the other conductors is 0.269 T at 25.7 degrees from the tape's normal at the limit and its measured critical current is 123.0 A;
# at 50 percent, 61.5 A per tape and 184.6 A in the coil.  The flat winding (tapes flat, edge to edge, 2 layers) gives 183.1 A with
# less tape but a cryostat three times wider; the compact pancake of the first check (0.05 mm insulation, 4 layers) 184.9 A.
# TAPE_LIMIT_MODE = "legacy" restores the 66.7 A/mm x layers rule (environment ORBITAL_RING_TAPE_LIMIT=legacy).
TAPE_LIMIT_MODE = os.environ.get("ORBITAL_RING_TAPE_LIMIT", "measured")
WINDING_LAYOUT_MEASURED = "pancake"
HTS_TAPE_LAYERS_MEASURED = 3
I_TAPE_LIMIT_MEASURED = 123.0        # A per 3 mm tape, the worst tape at the coil limit's field and angle
B_TAPE_WORST_MEASURED = 0.269        # T at that tape, from the other conductors of the winding
THETA_TAPE_WORST_MEASURED = 25.7     # degrees between that field and the tape's normal
I_COIL_LIMIT_MEASURED = 184.6        # A, the coil current at which the worst tape reaches 50 percent of its limit

# -----------------------------------------------------------------------------
# 1.2 LIM Geometry
# -----------------------------------------------------------------------------
N_TURNS = 208           # Turns per phase coil (typically 50-200)
TAU_P = 100.0           # Pole pitch in meters
W_COIL = 2.0            # LIM coil width in meters
GAP = 0.20              # Air gap between coil and reaction plate (m)
T_PLATE = 0.2           # Reaction plate thickness (m)
PITCH_COUNT = 3         # Number of pole pitches per LIM

# -----------------------------------------------------------------------------
# 1.2b Plate Geometry
# -----------------------------------------------------------------------------
W_PLATE = 1.2               # Reaction plate width (m)
N_PLATES_PER_SIDE = 1       # LIM reaction plates per side of cable
N_LEV_PLATES = 4            # Ferromagnetic bearing plates on the cable: two support plates below, two lift plates on top (2026-10-12)
D_LEV = 0.20                # Bearing plate thickness (m), 200 mm (was 0.10 with two plates until 2026-10-12)
D_LEV_COIL = 1.0            # Levitation coild witdth (m)
W_LEV = D_LEV_COIL + 0.4    # Levitation plate width (m)

# Reaction plate material: "aluminum", "cuni7030", "titanium", "alpha_titanium", "gamma_titanium"
PLATE_MATERIAL = "gamma_titanium"

# -----------------------------------------------------------------------------
# 1.3 LIM Spacing and Site Configuration
# -----------------------------------------------------------------------------
LIM_SPACING = 500.0     # Distance between LIM sites (m)

# -----------------------------------------------------------------------------
# 1.4 Operating Limits
# -----------------------------------------------------------------------------
VOLTS_MAX = 100e3       # Maximum induced coil voltage (V)
MAX_SITE_POWER = 8e6   # Maximum power per LIM site (W)

# -----------------------------------------------------------------------------
# 1.5 Slip Control Parameters
# -----------------------------------------------------------------------------
V_SLIP_MIN = 5.0        # Minimum slip velocity (m/s)
V_SLIP_MAX = 200.0      # Maximum slip velocity (m/s)
SLIP_RATIO_NORMAL = 0.02    # Target slip ratio at full current (2%)
SLIP_RATIO_REDUCED = 0.005  # Reduced slip when power-limited (0.5%)

# -----------------------------------------------------------------------------
# 1.6 Thrust Model Selection
# -----------------------------------------------------------------------------
THRUST_MODEL = 1        # Select model: 1, 2, or 3
THRUST_EFFICIENCY = 1.0 # Multiplier on calculated thrust

# -----------------------------------------------------------------------------
# 1.7 Mass Configuration
# -----------------------------------------------------------------------------
M_CABLE_STRUCTURAL_DEFAULT = ring.M_CABLE_STRUCTURAL  # quadratic sizing at m_load=12000, design altitude (ring_altitude.py)
M_LOAD_M_DEFAULT = 12_000               # Casing + payload mass per meter (kg/m)
M_CABLE_STRUCTURAL = M_CABLE_STRUCTURAL_DEFAULT # Session proofing
M_LOAD_M = M_LOAD_M_DEFAULT
SIGMA_CABLE = 12.5e9                    # Operating stress (Pa), safety factor 2.0 on 25 GPa

# -----------------------------------------------------------------------------
# 1.8 Simulation Control
# -----------------------------------------------------------------------------
SAMPLE_TIME_MAX = 5 * 365.33 * 24 * 3600  # Maximum simulation time (5 years)
DT1 = 1                 # Time step for first day (s)
DT2 = 10                # Time step until sample_time_max (s)
DT3 = 50                # Time step after sample_time_max (s)
SKIP = 200              # Data collection interval

# -----------------------------------------------------------------------------
# 1.9 Output Control
# -----------------------------------------------------------------------------
WRITE_FILE = True
MAKE_GRAPHS = True

# -----------------------------------------------------------------------------
# 1.9b Graph Output Configuration
# -----------------------------------------------------------------------------
SAVE_GRAPHS = True                      # Save graphs to files instead of displaying
GRAPH_OUTPUT_DIR = os.path.join(ring.OUT_DIR, "graphs_lim")   # figures go under ORBITAL_RING_OUT (default runs)
GRAPH_DPI = 300                         # Resolution for print (300 DPI is standard for print)
GRAPH_WIDTH_INCHES = 10                 # Width in inches (10" at 300 DPI = 3000 pixels)
GRAPH_HEIGHT_INCHES = 5                 # Height in inches (6" at 300 DPI = 1800 pixels)
GRAPH_FORMAT = "png"                    # Output format: "png", "pdf", "svg", "jpg"

# -----------------------------------------------------------------------------
# 1.10 Thermal Mode Configuration
# -----------------------------------------------------------------------------
EDDY_HEAT_TO_CABLE = True      # True = cable heats up, False = cryo handles it
HEATSINK_LENGTH_GAP = 50       # reduces the length of the radiators to allow for other stuctures (purely a design issue) (m)

# Thermal model selection (when EDDY_HEAT_TO_CABLE = True):
#   "legacy"   = Original single-node model (cable radiates to space through fudge factor)
#   "two_node" = Two-node model (cable → casing wall → external radiators via thermosyphon)
WARM_THERMAL_MODEL = "two_node"

# Cable thermal properties
CABLE_EMISSIVITY = 0.85
CABLE_SURFACE_AREA_PER_M = 0.5  # Legacy model: effective radiating area (fudge factor) (m²/m)

# Two-node model: cable-to-wall radiation
# The cable is ~7.5 m wide with a CNT structural cross-section of ~57 m².
# Exposed cable perimeter facing the casing inner wall, excluding MLI-wrapped
# levitation coils (which are at 77 K and handled by the cryo system).
# Conservative estimate: top and bottom faces of the structural cable.
CABLE_RADIATING_AREA_PER_M = 7.5  # Cable surface area facing casing wall (7.5 m ~ top surface width) (m²/m)

# Two-node model: casing wall and external radiators
CASING_WALL_EMISSIVITY = 0.85      # Emissivity of casing inner wall
WARM_RADIATOR_AREA_PER_M = 0.75    # External warm-loop radiator area (m²/m)
WARM_RADIATOR_EMISSIVITY = 0.90    # Emissivity of external warm radiators

# LIM coil thermal isolation
COIL_MLI_EFFECTIVENESS = 0.001   # MLI transmission for LIM coils (0.1%)
# COIL_WINDING_WIDTH, COIL_SURFACE_AREA_PER_COIL and COIL_SURFACE_AREA_PER_SITE are derived in Section 4 (they depend on the layer count):
# Kapton tape + adhesive + HTS thickness per turn, times the turns (m); the cryostat area per coil and per site from it.

# Levitation coil thermal properties
# The DC levitation coils run continuously along the ring (no gaps).
# They absorb heat from the warm cable above through their MLI insulation.
LEV_COIL_WIDTH = D_LEV_COIL * N_LEV_PLATES # Combined width of levitation coil rings (m)
LEV_COIL_MLI_HEAT_FLUX_REF = 5.0  # MLI heat flux at 300K reference (W/m²)

# Note: Environmental heat (Q_ABSORBED_PER_M, calculated in Section 4) is added to the
# cable's heat budget when EDDY_HEAT_TO_CABLE = True. This represents solar and earth
# radiation that leaks through the MLI shielding (~100 W/m with current parameters).
# This sets a minimum cable temperature of ~254 K even with zero eddy losses.


# =============================================================================
# SECTION 2: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7
STEFAN_BOLTZMANN = 5.670374e-8

# Orbital parameters at the design altitude (ring_altitude.py)
V_ORBIT = ring.V_ORBIT
# Casing end velocity in the ring's launch direction (negative for the
# retrograde ring: the casing passes through zero and ends eastward).
V_GROUND_STATIONARY = ring.V_CASING_FINAL_LAUNCH_FRAME
L_RING = ring.L_RING

# =============================================================================
# SECTION 3: MATERIAL PROPERTIES
# =============================================================================

# Aluminum
RHO_ALU_293K = 2.65e-8
RHO_ALU_70K = 4.853e-9
DENSITY_ALU = 2700
C_P_ALU = 900
K_ALU = 205
EM_ALU = 0.85
ALPHA_ALU = 3.663e-3

# CuNi 70/30 (Cupronickel)
RHO_CUNI_293K = 38e-8
DENSITY_CUNI = 8900
C_P_CUNI = 377
K_CUNI = 29
EM_CUNI = 0.65
ALPHA_CUNI = 0.0004

# Pure Titanium
RHO_TI_293K = 42e-8
DENSITY_TI = 4500
C_P_TI = 520
K_TI = 22
EM_TI = 0.60
ALPHA_TI = 0.0035

# Alpha-2 titanium aluminide, α₂-Ti₃Al
RHO_aTI_293K = 50e-8
DENSITY_aTI = 4200
C_P_aTI = 550
K_aTI = 17
EM_aTI = 0.60
ALPHA_aTI = 0.0015

# Gamma titanium aluminide, Ti-48Al-2Cr-2Nb
RHO_gTI_293K = 75e-8
DENSITY_gTI = 3900
C_P_gTI = 570
K_gTI = 15
EM_gTI = 0.60
ALPHA_gTI = 0.0012

# Iron (for levitation plates)
DENSITY_IRON = 7870

# Liquid nitrogen
T_LN2_BOIL = 77.4
T_LN2_SUPPLY = 70
C_P_LN2 = 2040
L_V_LN2 = 199000
H_CONV_LN2 = 90

# Thermal environment
Q_SUN = 1361
Q_EARTH_ALBEDO = 650
Q_SHIELDING = 0.005
T_SPACE = 2.7
T_RADIATOR_HOT = 400
T_MAX_PLATE = 500

# Cryogenic system
CRYO_EFF = 0.05
EM_HEATSINK = 0.9

# HTS tape
HTS_THICKNESS_UM = 80
ALPHA_PENETRATION_DEG = 20.0


# =============================================================================
# SECTION 4: DERIVED PARAMETERS
# =============================================================================

# Time constants
HR = 3600
DAY = 24 * HR
WEEK = 7 * DAY
MONTH = 30 * DAY
YR = round(365.33 * DAY)

# Material selection
def get_material_properties(material):
    """Return material properties dict for given material name."""
    if material == "cuni7030":
        return {
            'density': DENSITY_CUNI, 'rho_293K': RHO_CUNI_293K,
            'alpha': ALPHA_CUNI, 'cp': C_P_CUNI, 'k': K_CUNI, 'em': EM_CUNI
        }
    elif material == "titanium":
        return {
            'density': DENSITY_TI, 'rho_293K': RHO_TI_293K,
            'alpha': ALPHA_TI, 'cp': C_P_TI, 'k': K_TI, 'em': EM_TI
        }
    elif material == "alpha_titanium":
        return {
            'density': DENSITY_aTI, 'rho_293K': RHO_aTI_293K,
            'alpha': ALPHA_aTI, 'cp': C_P_aTI, 'k': K_aTI, 'em': EM_aTI
        }
    elif material == "gamma_titanium":
        return {
            'density': DENSITY_gTI, 'rho_293K': RHO_gTI_293K,
            'alpha': ALPHA_gTI, 'cp': C_P_gTI, 'k': K_gTI, 'em': EM_gTI
        }
    else:  # aluminum
        return {
            'density': DENSITY_ALU, 'rho_293K': RHO_ALU_293K,
            'alpha': ALPHA_ALU, 'cp': C_P_ALU, 'k': K_ALU, 'em': EM_ALU
        }


def _derive():
    """The derived parameters of Section 4, from the user parameters of Section 1; configure() reruns it."""
    global W_TAPE, eff_hts_tape_layers, I_C, I_PEAK, I_TARGET, I_MIN, LIMS_PER_SITE, L_ACTIVE, A_LIM, A_COIL, L_HTS_COIL, L_HTS_LIM, LIM_PHASES, LIM_SITES, ALPHA_TAPE, MATERIAL_PROPS, PLATE_DENSITY, PLATE_RHO_293K, PLATE_ALPHA, PLATE_CP, PLATE_K, PLATE_EM, M_LIM_PLATE, M_LEV_PLATE, M_HARDWARE, M_CABLE_M, M_CABLE_TOTAL, M_LOAD_TOTAL, CASING_WIDTH, Q_ABSORBED_PER_M, Q_ABSORBED_PER_SITE, MAX_HEATSINK_AREA, HEATSINK_LENGTH, V_REL_MIN_FUDGE, Q_COIL_ENVIRONMENT, LEV_COIL_AREA_PER_SITE, INV_EFF, LIM_EFF, SLIP_MIN, CURRENT_UPRATE, POWER_HEADROOM, COIL_WINDING_WIDTH, COIL_SURFACE_AREA_PER_COIL, COIL_SURFACE_AREA_PER_SITE, HTS_TAPE_LAYERS, I_TAPE_LIMIT, WINDING_LAYOUT
    # HTS current ratings
    W_TAPE = HTS_TAPE_WIDTH_MM / 1000

    if TAPE_LIMIT_MODE == "measured":
        # the measured basis (Section 1.1): the layer count the limit computation chose, the coil cap at the worst tape's 50 percent
        HTS_TAPE_LAYERS = HTS_TAPE_LAYERS_MEASURED
        WINDING_LAYOUT = WINDING_LAYOUT_MEASURED
        I_TAPE_LIMIT = I_TAPE_LIMIT_MEASURED
        eff_hts_tape_layers = float(HTS_TAPE_LAYERS)
        I_C = I_TAPE_LIMIT * HTS_TAPE_LAYERS          # the coil current at which the worst tape reaches its measured limit
        I_TARGET = I_COIL_LIMIT_MEASURED              # the cap the controller holds: the worst tape at 50 percent
        I_PEAK = I_TARGET                             # no transient above the cap
    else:
        # the legacy rule: a tape rating that ignores the winding's field, 66.7 A/mm per layer, de-rated for more than one layer
        WINDING_LAYOUT = "legacy"
        HTS_TAPE_LAYERS = HTS_TAPE_LAYERS_LEGACY
        eff_hts_tape_layers = float(HTS_TAPE_LAYERS)
        if HTS_TAPE_LAYERS > 1:
            eff_hts_tape_layers = DE_RATING_FACTOR * eff_hts_tape_layers
        I_C = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * eff_hts_tape_layers
        I_TAPE_LIMIT = I_C / max(eff_hts_tape_layers, 1.0)
        I_PEAK = 0.875 * I_C
        I_TARGET = 0.8125 * I_C
    I_MIN = 10.0

    # LIM coil cryostat (thermal): Kapton tape + adhesive + HTS thickness per turn, times the turns (m); the pancake's stack height
    COIL_WINDING_WIDTH = (1 + 0.06 + 0.12*HTS_TAPE_LAYERS)*N_TURNS/1000
    COIL_SURFACE_AREA_PER_COIL = 2 * (TAU_P + W_COIL) * COIL_WINDING_WIDTH   # Approximate LIM coil cryostat surface area (m²)
    COIL_SURFACE_AREA_PER_SITE = COIL_SURFACE_AREA_PER_COIL * 3 * PITCH_COUNT * 2 * N_PLATES_PER_SIDE  # area per coil * phases * pitches * sides * LIMs per side (m²)

    # LIM geometry
    LIMS_PER_SITE = 2 * N_PLATES_PER_SIDE
    L_ACTIVE = TAU_P * PITCH_COUNT
    A_LIM = L_ACTIVE * W_COIL
    A_COIL = TAU_P * W_COIL
    L_HTS_COIL = 2 * (W_COIL + TAU_P) * N_TURNS
    L_HTS_LIM = L_HTS_COIL * 3 * PITCH_COUNT
    LIM_PHASES = 3
    LIM_SITES = round(L_RING / LIM_SPACING)

    # Angular conversion
    ALPHA_TAPE = ALPHA_PENETRATION_DEG * math.pi / 180

    # Get current material properties
    MATERIAL_PROPS = get_material_properties(PLATE_MATERIAL)
    PLATE_DENSITY = MATERIAL_PROPS['density']
    PLATE_RHO_293K = MATERIAL_PROPS['rho_293K']
    PLATE_ALPHA = MATERIAL_PROPS['alpha']
    PLATE_CP = MATERIAL_PROPS['cp']
    PLATE_K = MATERIAL_PROPS['k']
    PLATE_EM = MATERIAL_PROPS['em']

    # Mass calculations
    M_LIM_PLATE = PLATE_DENSITY * T_PLATE * W_PLATE
    M_LEV_PLATE = DENSITY_IRON * D_LEV * W_LEV
    M_HARDWARE = (N_LEV_PLATES * M_LEV_PLATE + 2 * N_PLATES_PER_SIDE * M_LIM_PLATE)
    M_CABLE_M = M_CABLE_STRUCTURAL + M_HARDWARE

    # Total masses
    M_CABLE_TOTAL = M_CABLE_M * L_RING
    M_LOAD_TOTAL = M_LOAD_M * L_RING

    # Thermal
    CASING_WIDTH = 10.0
    Q_ABSORBED_PER_M = (Q_SUN + Q_EARTH_ALBEDO) * CASING_WIDTH * Q_SHIELDING
    Q_ABSORBED_PER_SITE = Q_ABSORBED_PER_M * LIM_SPACING
    MAX_HEATSINK_AREA = LIM_SPACING * 2 * W_COIL
    HEATSINK_LENGTH = LIM_SPACING - HEATSINK_LENGTH_GAP
    V_REL_MIN_FUDGE = 10

    # Coil environmental heat leak (LIM coils only)
    Q_COIL_ENVIRONMENT = (Q_SUN + Q_EARTH_ALBEDO) * COIL_SURFACE_AREA_PER_SITE * COIL_MLI_EFFECTIVENESS

    # Levitation coil area per site (for heat load calculation)
    LEV_COIL_AREA_PER_SITE = LEV_COIL_WIDTH * LIM_SPACING  # m²

    # Efficiency factors
    INV_EFF = 0.90
    LIM_EFF = 0.95

    # Controller parameters
    SLIP_MIN = 0.005
    CURRENT_UPRATE = 1.01
    POWER_HEADROOM = 0.98




_derive()


def configure(**overrides):
    """Set user parameters (N_TURNS, MAX_SITE_POWER, TAPE_LIMIT_MODE, HTS_TAPE_LAYERS, ...) and recompute the derived ones,
    for runs that vary a parameter: import lim_config; lim_config.configure(N_TURNS=150, MAX_SITE_POWER=16e6)."""
    g = globals()
    for k, v in overrides.items():
        if k not in g:
            raise KeyError(f"lim_config has no parameter {k}")
        g[k] = v
    _derive()


# =============================================================================
# SECTION 5: PARAMETER DICTIONARY FOR PHYSICS FUNCTIONS
# =============================================================================

def get_physics_params():
    """Return parameter dictionary for physics functions."""
    return {
        'n_turns': N_TURNS,
        'tau_p': TAU_P,
        'w_coil': W_COIL,
        'gap': GAP,
        't_plate': T_PLATE,
        'pitch_count': PITCH_COUNT,
        'w_plate': W_PLATE,
        'l_active': L_ACTIVE,
        'a_lim': A_LIM,
        'a_coil': A_COIL,
        'l_hts_coil': L_HTS_COIL,
        'lim_phases': LIM_PHASES,
        'rho_293K': PLATE_RHO_293K,
        'alpha': PLATE_ALPHA,
        'material': PLATE_MATERIAL,
        'density': PLATE_DENSITY,
        'cp': PLATE_CP,
        'em': PLATE_EM,
        'v_slip_min': V_SLIP_MIN,
        'thrust_efficiency': THRUST_EFFICIENCY,
        'w_tape': W_TAPE,
        'alpha_tape': ALPHA_TAPE,
        'i_c': I_TAPE_LIMIT if TAPE_LIMIT_MODE == "measured" else I_C,   # per tape (measured) or the legacy coil-level figure
        'layers': HTS_TAPE_LAYERS if TAPE_LIMIT_MODE == "measured" else 1,
    }


