#!/usr/bin/env python3
"""
Linear Induction Motor Deployment Simulation for Orbital Ring Systems

This simulation models the deployment phase of an orbital ring, where Linear
Induction Motors (LIMs) accelerate a magnetically levitated cable from orbital
velocity (7.75 km/s) down to ground-stationary velocity (483 m/s at 250 km).

The simulation tracks:
  - Electromagnetic thrust from three different physics models
  - Thermal dynamics of the aluminium reaction plate  
  - Cryogenic power requirements for HTS coil cooling
  - Power consumption at each LIM site
  - Deployment timeline

Reference: "Orbital Ring Engineering" by Paul G de Jong
Companion code for: "Ion Propulsion Engineering" by Paul G de Jong

Usage:
    python lim_deployment_sim.py [--model=N] [graph_options]
    
    --model=1  Narrow plate eddy current (default, geometry-corrected)
    --model=2  Goodness factor model (Laithwaite)
    --model=3  Slip × pressure (theoretical maximum)
    
    Run with --help for full options.
"""

import sys
import math
# import tabulate
import matplotlib.pyplot as plt


# =============================================================================
# SECTION 1: USER-CONFIGURABLE PARAMETERS
# =============================================================================
# These are the parameters you'll most likely want to modify when exploring
# different LIM designs. Each parameter includes valid ranges and trade-offs.

# -----------------------------------------------------------------------------
# 1.1 HTS (High-Temperature Superconductor) Tape Configuration
# -----------------------------------------------------------------------------
# Commercial HTS tape (e.g., SuperPower, AMSC) comes in standard widths.
# Critical current scales linearly with width and number of layers.
# 
# Trade-offs:
#   - Wider tape / more layers → higher current → more thrust, but higher cost
#   - Narrower tape → lower current → longer deployment, but cheaper

HTS_TAPE_WIDTH_MM = 3       # Standard widths: 12, 6, 4, or 3 mm
HTS_TAPE_LAYERS = 3         # Layers of tape in parallel: 1, 2, or 3
IC_PER_MM_PER_LAYER = 66.7  # Critical current density: ~66.7 A/mm per layer

# -----------------------------------------------------------------------------
# 1.2 LIM Geometry
# -----------------------------------------------------------------------------
# These parameters define the physical dimensions of each Linear Induction Motor.
#
# Key relationships:
#   - Pole pitch (TAU_P) sets the wavelength of the traveling magnetic field
#   - Longer pole pitch → lower frequency for same slip velocity
#   - Coil width (W_COIL) must be << TAU_P for narrow-plate model validity
#   - Gap affects field strength: B ∝ arctan(W/2g), smaller gap = stronger field

N_TURNS = 80            # Turns per phase coil (typically 50-200)
TAU_P = 100.0           # Pole pitch in meters (wavelength/2 of traveling field)
W_COIL = 0.5            # LIM coil width in meters (must be << TAU_P)
GAP = 0.05              # Air gap between coil and reaction plate (m)
T_PLATE = 0.1           # Aluminium reaction plate thickness (m)
PITCH_COUNT = 3         # Number of pole pitches per LIM (typically 2-4)

# -----------------------------------------------------------------------------
# 1.3 LIM Spacing and Site Configuration
# -----------------------------------------------------------------------------
# LIMs are distributed around the ring. Each "site" has LIMs on both sides
# of the cable for balanced forces.
#
# Trade-offs:
#   - Closer spacing → more LIMs → faster deployment, higher capital cost
#   - Wider spacing → fewer LIMs → slower deployment, lower capital cost

LIM_SPACING = 500.0     # Distance between LIM sites (m)
LIMS_PER_SIDE = 1       # LIMs on each side of cable at each site (1, 2, or 3)

# -----------------------------------------------------------------------------
# 1.4 Operating Limits
# -----------------------------------------------------------------------------
# These limits protect the HTS tape from quenching and prevent thermal runaway.
#
# Critical limits:
#   - VOLTS_MAX: Prevents dielectric breakdown and excessive dV/dt stress
#   - MAX_SITE_POWER: Total electrical power available at each LIM site

VOLTS_MAX = 100e3       # Maximum induced coil voltage (V) - prevents quench
MAX_SITE_POWER = 16e6   # Maximum power per LIM site (W) - 16 MW default

# -----------------------------------------------------------------------------
# 1.5 Slip Control Parameters
# -----------------------------------------------------------------------------
# Slip is the difference between the traveling wave speed and the reaction 
# plate speed. It determines thrust and losses.
#
# Key insight: Higher slip → more eddy current losses but potentially more thrust
#              Lower slip → less loss but may reduce thrust in resistive regime

V_SLIP_MIN = 5.0        # Minimum slip velocity (m/s) - ensures startup thrust
V_SLIP_MAX = 200.0      # Maximum slip velocity (m/s) - limits thermal load
SLIP_RATIO_NORMAL = 0.02    # Target slip ratio at full current (2%)
SLIP_RATIO_REDUCED = 0.005  # Reduced slip when power-limited (0.5%)

# -----------------------------------------------------------------------------
# 1.6 Thrust Model Selection
# -----------------------------------------------------------------------------
# Three physics models are available, each with different assumptions:
#
#   Model 1 (default): Narrow plate eddy current model
#       - Accounts for W << TAU_P geometry of orbital ring
#       - Models elongated eddy current loops with high return-path resistance
#       - Most physically accurate for this application
#
#   Model 2: Goodness factor model (Laithwaite)
#       - Classic LIM theory assuming W > TAU_P
#       - May overestimate thrust for narrow-plate geometry
#       - Useful for comparison with literature
#
#   Model 3: Slip × pressure model
#       - Theoretical maximum: F = slip × (B²A/2μ₀)
#       - Upper bound assuming perfect phase alignment
#       - Useful for understanding limits

THRUST_MODEL = 1        # Select model: 1, 2, or 3
THRUST_EFFICIENCY = 1.0 # Multiplier on calculated thrust (1.0 = no adjustment)

# -----------------------------------------------------------------------------
# 1.7 Mass Configuration
# -----------------------------------------------------------------------------
# These masses come from the orbital ring design in "Orbital Ring Engineering".
# They determine the acceleration achievable for a given thrust.

M_CABLE_M = 99_198      # Cable mass per meter (kg/m) - the rotor
M_LOAD_M = 12_000       # Casing + payload mass per meter (kg/m) - the stator

# -----------------------------------------------------------------------------
# 1.8 Simulation Control
# -----------------------------------------------------------------------------
# Time stepping and data collection parameters.

SAMPLE_TIME_MAX = 5 * 365.33 * 24 * 3600  # Maximum simulation time (5 years)
DT1 = 1                 # Time step for first day (s)
DT2 = 10                # Time step until sample_time_max (s)  
DT3 = 50                # Time step after sample_time_max (s)
SKIP = 200              # Data collection interval (1 = every step)

# -----------------------------------------------------------------------------
# 1.9 Output Control
# -----------------------------------------------------------------------------

WRITE_FILE = True       # Write results to file
MAKE_GRAPHS = True      # Generate matplotlib graphs


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

# Aluminium properties
RHO_ALU_70K = 4.853e-9          # Resistivity at 70 K (Ω·m)
RHO_ALU_MASS = 2700             # Density (kg/m³)
C_P_ALU = 900                   # Specific heat capacity (J/kg·K)
K_ALU = 205                     # Thermal conductivity (W/m·K)
EM_ALU = 0.85                   # Emissivity (black anodized)

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

# HTS current ratings
# The 87.5% and 81.25% factors provide safety margin below critical current
# to prevent thermal runaway near Ic where losses increase sharply.
W_TAPE = HTS_TAPE_WIDTH_MM / 1000
I_C = IC_PER_MM_PER_LAYER * HTS_TAPE_WIDTH_MM * HTS_TAPE_LAYERS
I_PEAK = 0.875 * I_C            # Peak operating current (87.5% of Ic)
I_TARGET = 0.8125 * I_C         # Target steady-state current (81.25% of Ic)
I_MIN = 10.0                    # Minimum current (A)

# LIM geometry
LIMS_PER_SITE = 2 * LIMS_PER_SIDE
L_ACTIVE = TAU_P * PITCH_COUNT
A_LIM = L_ACTIVE * W_COIL
A_COIL = TAU_P * W_COIL
L_HTS_COIL = 2 * (W_COIL + TAU_P) * N_TURNS
L_HTS_LIM = L_HTS_COIL * 3 * PITCH_COUNT  # 3 phases
LIM_PHASES = 3
LIM_SITES = round(L_RING / LIM_SPACING)

# Angular conversion
ALPHA_TAPE = ALPHA_PENETRATION_DEG * math.pi / 180

# Total masses
M_CABLE_TOTAL = M_CABLE_M * L_RING
M_LOAD_TOTAL = M_LOAD_M * L_RING

# Thermal
CASING_WIDTH = 10.0             # External casing width (m)
Q_ABSORBED_PER_M = (Q_SUN + Q_EARTH_ALBEDO) * CASING_WIDTH * Q_SHIELDING
Q_ABSORBED_PER_SITE = Q_ABSORBED_PER_M * LIM_SPACING
MAX_HEATSINK_AREA = LIM_SPACING * 2 * W_COIL
HEATSINK_LENGTH = LIM_SPACING
V_REL_MIN_FUDGE = 10            # Minimum v_rel for thermal calculations (m/s)

# Efficiency factors
INV_EFF = 0.90                  # DC-AC inverter efficiency
LIM_EFF = 0.95                  # Miscellaneous LIM efficiency

# Controller parameters
SLIP_MIN = 0.005
CURRENT_UPRATE = 1.01
POWER_HEADROOM = 0.98


# =============================================================================
# SECTION 4: DATA COLLECTION
# =============================================================================
# Lists to store time-series data for plotting.

data_current = []
data_voltage = []
data_v_slip = []
data_slip_ratio = []
data_f_slip = []
data_thrust = []
data_thrust_power = []
data_v_rel = []
data_b_field = []
data_p_eddy = []
data_skin_depth_eff = []
data_skin_depth_calc = []
data_p_cryo = []
data_p_hyst = []
data_p_lim = []
data_p_site = []
data_temp_plate = []
data_E_site = []
data_E_total = []


# =============================================================================
# SECTION 5: PHYSICS FUNCTIONS - Orbital Dynamics
# =============================================================================

def calc_cable_velocity(v0, thrust_per_lim, dt):
    """Update cable velocity based on thrust.
    
    Newton's second law: F = ma, so dv = (F_total/m) * dt
    Total thrust = thrust_per_lim × LIM_SITES × 2 (both sides of ring)
    """
    thrust_total = thrust_per_lim * LIM_SITES * 2
    return v0 + (thrust_total / M_CABLE_TOTAL) * dt


def calc_casing_velocity(v0, thrust_per_lim, dt):
    """Update casing velocity based on reaction force.
    
    Newton's third law: The thrust accelerating the cable creates an equal
    and opposite reaction on the casing (which carries the LIM stators).
    """
    thrust_total = thrust_per_lim * LIM_SITES * 2
    return v0 - (thrust_total / M_LOAD_TOTAL) * dt


def calc_relative_velocity(v_cable, v_casing):
    """Relative velocity between cable (rotor) and casing (stator)."""
    v_rel = v_cable - v_casing
    if v_rel < 0:
        print(f"WARNING: Negative v_rel = {v_rel:.2f} m/s")
    return v_rel


# =============================================================================
# SECTION 6: PHYSICS FUNCTIONS - Slip and Frequency
# =============================================================================

def calc_slip_velocity(slip_ratio, v_wave):
    """Calculate slip velocity from slip ratio and wave speed."""
    if slip_ratio < 0 or slip_ratio > 1:
        print(f"WARNING: Invalid slip ratio = {slip_ratio:.4f}")
    return slip_ratio * v_wave


def calc_slip_ratio(v_slip, v_rel):
    """Calculate slip ratio: s = v_slip / v_wave where v_wave = v_slip + v_rel."""
    v_wave = v_slip + v_rel
    if v_wave < 1:
        v_wave = 1.0
    return v_slip / v_wave


def calc_slip_frequency(v_slip):
    """Slip frequency in Hz.
    
    The traveling field moves at v_wave = 2 * TAU_P * f_supply.
    The slip frequency is f_slip = v_slip / (2 * TAU_P).
    """
    return v_slip / (2 * TAU_P)


def calc_supply_frequency(v_wave):
    """Electrical supply frequency for a given wave speed."""
    return v_wave / (2 * TAU_P)


# =============================================================================
# SECTION 7: PHYSICS FUNCTIONS - Magnetic Fields
# =============================================================================

def calc_b_field_at_plate(i_peak):
    """Traveling wave magnetic field amplitude at the reaction plate.
    
    For a rectangular current sheet (the coil) at distance g (gap) from the plate:
        B_single = (2μ₀NI/πW) × arctan(W/2g)
    
    For balanced 3-phase excitation, the traveling wave amplitude is:
        B_traveling = 1.5 × B_single
    
    This is a constant-amplitude wave (no pulsation), so no RMS conversion needed.
    """
    b_single = (2 * MU0 * N_TURNS * i_peak / (math.pi * W_COIL)) * math.atan(W_COIL / (2 * GAP))
    return 1.5 * b_single


def calc_b_field_in_coil(i_peak):
    """Peak magnetic field inside the coil (for voltage and hysteresis calculations).
    
    Each coil experiences its own field for flux linkage purposes.
    """
    return MU0 * N_TURNS * i_peak / W_COIL


# =============================================================================
# SECTION 8: PHYSICS FUNCTIONS - Material Properties
# =============================================================================

def calc_resistivity_alu(temp_K):
    """Temperature-dependent resistivity of aluminium.
    
    Uses linear temperature coefficient:
        ρ(T) = ρ₀ × [1 + α × (T - T₀)]
    
    where ρ₀ = 2.65×10⁻⁸ Ω·m at 293 K, α = 3.663×10⁻³ K⁻¹
    
    Clamped to minimum value at 70 K to prevent unphysical extrapolation.
    """
    rho_293K = 2.65e-8
    alpha = 3.663e-3
    rho = rho_293K * (1 + alpha * (temp_K - 293))
    return max(rho, RHO_ALU_70K)


def calc_skin_depth(f_slip, temp_K):
    """Electromagnetic skin depth in the reaction plate.
    
    δ = √(ρ / πμ₀f)
    
    Eddy currents are confined to approximately this depth at frequency f.
    """
    rho = calc_resistivity_alu(temp_K)
    if f_slip <= 0:
        f_slip = calc_slip_frequency(V_SLIP_MIN)
    return math.sqrt(rho / (math.pi * MU0 * f_slip))


def calc_effective_plate_depth(f_slip, temp_K):
    """Effective conducting depth: minimum of plate thickness and skin depth."""
    delta = calc_skin_depth(f_slip, temp_K)
    return min(T_PLATE, delta)


# =============================================================================
# SECTION 9: PHYSICS FUNCTIONS - Eddy Current Losses
# =============================================================================

def calc_eddy_loss_from_thrust(v_slip, thrust):
    """Reaction plate ohmic loss consistent with thrust (per LIM).
    
    In any induction machine, slip power goes to the secondary (reaction plate):
        P_secondary = F × v_slip
    
    This self-consistent approach prevents the runaway that occurs when
    eddy losses are calculated independently of the field loading.
    """
    if v_slip <= 0 or thrust <= 0:
        return 0.0
    return thrust * v_slip


def calc_plate_eddy_volume(f_slip, temp_K):
    """Volume of reaction plate participating in eddy currents."""
    depth = calc_effective_plate_depth(f_slip, temp_K)
    return W_COIL * L_ACTIVE * depth


# =============================================================================
# SECTION 10: PHYSICS FUNCTIONS - Goodness Factor (for Model 2)
# =============================================================================

def calc_goodness_factor(f_slip, temp_K):
    """Laithwaite's goodness factor G.
    
    G determines the LIM operating regime:
        G << 1: Resistive regime (thrust ∝ slip)
        G ~ 1:  Transitional (peak thrust at optimal slip)
        G >> 1: Inductive regime (thrust limited by back-EMF)
    
    G = (ω_slip × μ₀ × σ × δ_eff × τ_p) / π
    
    For orbital ring LIMs with τ_p = 100m, G is typically 500-1000,
    placing operation firmly in the inductive regime.
    
    Reference: Laithwaite (1965), Boldea & Nasar (1976)
    """
    if f_slip <= 0:
        return 0.0
    
    rho = calc_resistivity_alu(temp_K)
    sigma = 1.0 / rho
    omega = 2 * math.pi * f_slip
    delta_eff = calc_effective_plate_depth(f_slip, temp_K)
    
    return (omega * MU0 * sigma * delta_eff * TAU_P) / math.pi


def calc_slip_efficiency(slip_ratio, G):
    """Slip efficiency from goodness factor model.
    
    η_slip = 2sG / (1 + s²G²)
    
    This captures both regimes:
        Low G: η_slip ≈ 2sG (resistive, thrust ∝ slip)
        High G: η_slip ≈ 2/(sG) (inductive, thrust drops at high slip)
        Peak at s = 1/G: η_slip = 1
    """
    if slip_ratio <= 0 or G <= 0:
        return 0.0
    sG = slip_ratio * G
    return (2 * sG) / (1 + sG ** 2)


# =============================================================================
# SECTION 11: THRUST MODELS
# =============================================================================

def calc_thrust_force_max(b_field):
    """Upper bound on thrust from magnetic pressure.
    
    F_max = (B²/2μ₀) × A
    
    This is the magnetic pressure times the coupling area.
    """
    return (b_field ** 2) * A_LIM / (2 * MU0)


def calc_thrust_power_max(b_field, v_slip):
    """Upper bound on power transfer from magnetic energy flux."""
    return (b_field ** 2) * A_LIM * v_slip / (2 * MU0)


def calc_loop_inductance():
    """Inductance of eddy current loop for narrow plate geometry.
    
    For a rectangular loop of length τ_p and width W:
        L = (μ₀ × τ_p / π) × ln(τ_p / W)
    
    This is significant when W << τ_p (the orbital ring case).
    """
    return (MU0 * TAU_P / math.pi) * math.log(TAU_P / W_COIL)


# -----------------------------------------------------------------------------
# Model 1: Narrow Plate Eddy Current Model
# -----------------------------------------------------------------------------
def calc_thrust_model1(f_slip, f_supply, i_peak, temp_K=77.0):
    """Thrust from narrow plate eddy current model.
    
    For the orbital ring geometry (W << τ_p), eddy currents form elongated
    loops spanning one pole pitch. The return path along length τ_p dominates
    the loop resistance.
    
    Eddy current loop geometry:
    
        ←←←←←← current under negative pole ←←←←←←
        ↓                                        ↑
        ↓ return path (length τ_p)               ↑
        ↓                                        ↑
        →→→→→→ current under positive pole →→→→→→
        
        |←—————————— τ_p = 100m ——————————→|
    
    Number of loops = L_ACTIVE / τ_p (NOT related to N_TURNS)
    
    Physics:
        1. EMF per loop ≈ (4/π) × v_slip × B × W (both legs contribute)
        2. Loop impedance Z = √(R² + X²)
        3. Loop current I = EMF / Z
        4. Power dissipated P = N_loops × I² × R
        5. Thrust F = P / v_slip (from power balance)
    """
    if i_peak <= 0 or f_supply <= 0 or f_slip <= 0:
        return 0.0
    
    B = calc_b_field_at_plate(i_peak)
    v_slip = f_slip * 2 * TAU_P
    if v_slip <= 0:
        return 0.0
    
    # Material properties
    rho = calc_resistivity_alu(temp_K)
    delta = calc_skin_depth(f_slip, temp_K)
    d_eff = min(T_PLATE, delta)
    
    # EMF per loop (sinusoidal field average gives 4/π factor)
    EMF = (4.0 / math.pi) * v_slip * B * W_COIL
    
    # Loop impedance
    omega = 2 * math.pi * f_slip
    R = rho * TAU_P / (d_eff * W_COIL)
    L = calc_loop_inductance()
    X = omega * L
    Z = math.sqrt(R**2 + X**2)
    
    # Eddy current and power
    I = EMF / Z
    N_loops = L_ACTIVE / TAU_P
    P = N_loops * I**2 * R
    
    # Thrust from power balance
    F = P / v_slip
    return THRUST_EFFICIENCY * F


# -----------------------------------------------------------------------------
# Model 2: Goodness Factor Model (Laithwaite)
# -----------------------------------------------------------------------------
def calc_thrust_model2(f_slip, f_supply, i_peak, temp_K=77.0):
    """Thrust from goodness factor model.
    
    F = F_max × η_slip
    
    where η_slip = 2sG / (1 + s²G²)
    
    WARNING: This model assumes W > τ_p (wide plate). May overestimate
    thrust for the narrow-plate orbital ring geometry.
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    
    b_plate = calc_b_field_at_plate(i_peak)
    slip = f_slip / f_supply if f_supply > 0 else 0
    slip = max(0.001, min(1.0, slip))
    
    G = calc_goodness_factor(f_slip, temp_K)
    eta_slip = calc_slip_efficiency(slip, G)
    F_max = calc_thrust_force_max(b_plate)
    
    return THRUST_EFFICIENCY * F_max * eta_slip


# -----------------------------------------------------------------------------
# Model 3: Slip × Pressure Model (Theoretical Maximum)
# -----------------------------------------------------------------------------
def calc_thrust_model3(f_slip, f_supply, i_peak, temp_K=77.0):
    """Thrust from slip × pressure model.
    
    F = slip × F_max
    
    Assumes resistive regime with perfect phase alignment.
    This is a theoretical upper bound.
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    
    b_plate = calc_b_field_at_plate(i_peak)
    slip = f_slip / f_supply if f_supply > 0 else 0
    slip = max(0.0, min(1.0, slip))
    
    F_max = calc_thrust_force_max(b_plate)
    return THRUST_EFFICIENCY * slip * F_max


# -----------------------------------------------------------------------------
# Unified Thrust Function
# -----------------------------------------------------------------------------
def calc_thrust(f_slip, f_supply, i_peak, temp_K=77.0):
    """Calculate thrust using the selected model."""
    if THRUST_MODEL == 1:
        return calc_thrust_model1(f_slip, f_supply, i_peak, temp_K)
    elif THRUST_MODEL == 2:
        return calc_thrust_model2(f_slip, f_supply, i_peak, temp_K)
    elif THRUST_MODEL == 3:
        return calc_thrust_model3(f_slip, f_supply, i_peak, temp_K)
    else:
        return calc_thrust_model1(f_slip, f_supply, i_peak, temp_K)


# =============================================================================
# SECTION 12: ELECTRICAL CALCULATIONS
# =============================================================================

def calc_coil_voltage(i_peak, f_supply, b_field):
    """Induced voltage in LIM coil from changing flux linkage.
    
    V = ω × N × Φ_peak where Φ = B × A_coil
    """
    if i_peak <= 0 or f_supply <= 0:
        return 0.0
    phi_peak = b_field * A_COIL
    omega = 2 * math.pi * f_supply
    return omega * N_TURNS * phi_peak


def calc_thrust_power(thrust, v_rel):
    """Mechanical power delivered to the cable."""
    return thrust * v_rel


# =============================================================================
# SECTION 13: HTS HYSTERESIS LOSSES
# =============================================================================

def calc_hts_loss_factor(i_peak):
    """Hysteresis loss factor for HTS tape.
    
    Q = B × I × W_tape × sin(α)
    
    where α is the magnetic field penetration angle.
    """
    b_coil = calc_b_field_in_coil(i_peak)
    return b_coil * i_peak * W_TAPE * math.sin(ALPHA_TAPE)


def calc_hysteresis_power(f_supply, i_peak):
    """Total hysteresis power loss per LIM."""
    q = calc_hts_loss_factor(i_peak)
    p_coil = q * L_HTS_COIL * f_supply
    return p_coil * LIM_PHASES * PITCH_COUNT


# =============================================================================
# SECTION 14: THERMAL CALCULATIONS
# =============================================================================

def calc_effective_emissivity(em1=EM_ALU, em2=EM_HEATSINK):
    """Effective emissivity between two parallel surfaces."""
    return 1 / (1/em1 + 1/em2 - 1)


def calc_radiative_heat_transfer(area, T_hot, T_cold, em1=EM_ALU, em2=EM_HEATSINK):
    """Radiative heat transfer between two surfaces."""
    em_eff = calc_effective_emissivity(em1, em2)
    return em_eff * STEFAN_BOLTZMANN * area * (T_hot**4 - T_cold**4)


def calc_heatsink_area_required(p_heat, T_hot, T_cold=T_LN2_BOIL, em1=EM_ALU, em2=EM_HEATSINK):
    """Required heatsink area to radiate given power."""
    em_eff = calc_effective_emissivity(em1, em2)
    if T_hot <= T_cold:
        T_hot = T_cold + 1
    return p_heat / (em_eff * STEFAN_BOLTZMANN * (T_hot**4 - T_cold**4))


def calc_plate_temperature(v_rel, p_eddy):
    """Equilibrium temperature of reaction plate.
    
    The plate heats as it passes through the LIM and cools as it radiates
    in the gap between LIMs.
    """
    v_rel_adj = max(V_REL_MIN_FUDGE, v_rel)
    time_under_lim = L_ACTIVE / v_rel_adj
    time_between_lims = LIM_SPACING / v_rel_adj
    time_cycle = time_under_lim + time_between_lims
    
    em_eff = calc_effective_emissivity()
    # Solve Stefan-Boltzmann for T_hot
    term = (p_eddy * time_under_lim) / (em_eff * STEFAN_BOLTZMANN * (2 * A_LIM) * time_cycle)
    return (term + T_LN2_BOIL**4) ** 0.25


def calc_plate_temp_rise(v_rel, f_slip, p_eddy, temp_K):
    """Temperature rise of plate during one LIM passage."""
    v_rel_adj = max(V_REL_MIN_FUDGE, v_rel)
    time_under_lim = L_ACTIVE / v_rel_adj
    volume = calc_plate_eddy_volume(f_slip, temp_K)
    mass = volume * RHO_ALU_MASS
    return (p_eddy * time_under_lim) / (mass * C_P_ALU)


def calc_heat_load(p_eddy):
    """Total heat load including environmental absorption."""
    duty_cycle = L_ACTIVE / HEATSINK_LENGTH
    p_average = p_eddy * duty_cycle
    return p_average + Q_ABSORBED_PER_SITE


# =============================================================================
# SECTION 15: CRYOGENIC SYSTEM
# =============================================================================

def calc_cryo_cop(T_cold, T_hot=T_RADIATOR_HOT, efficiency=CRYO_EFF):
    """Coefficient of performance for cryogenic system.
    
    COP_actual = COP_carnot × efficiency
    COP_carnot = T_cold / (T_hot - T_cold)
    """
    if T_hot <= T_cold:
        T_hot = T_cold + 1
    cop_carnot = T_cold / (T_hot - T_cold)
    return cop_carnot * efficiency


def calc_cryo_power(heat_load, T_cold=T_LN2_BOIL, T_hot=T_RADIATOR_HOT, efficiency=CRYO_EFF):
    """Electrical power required to pump heat from cold to hot side."""
    cop = calc_cryo_cop(T_cold, T_hot, efficiency)
    return heat_load / cop


# =============================================================================
# SECTION 16: UTILITY FUNCTIONS
# =============================================================================

def calc_deployment_progress(v_casing):
    """Percentage completion of deployment (0% at start, 100% at ground-stationary)."""
    total_delta_v = V_ORBIT - V_GROUND_STATIONARY
    current_delta_v = V_ORBIT - v_casing
    return 100.0 * current_delta_v / total_delta_v


def make_month_ticks(data_list, total_time):
    """Generate x-axis tick positions and labels in months."""
    months = total_time / MONTH
    num_ticks = int(months) + 1
    n_points = len(data_list)
    
    tick_positions = [n_points * i / max(months, 1) for i in range(num_ticks)]
    tick_labels = [str(i) for i in range(num_ticks)]
    return tick_positions, tick_labels


def annotate_final_value(data_list, unit="", fmt=".1f"):
    """Add annotation showing final value on plot."""
    if not data_list:
        return
    x = len(data_list) - 1
    y = data_list[-1]
    label = f"{y:{fmt}} {unit}".strip()
    plt.annotate(
        label, xy=(x, y), xytext=(-40, 10),
        textcoords="offset points", fontsize=10, ha="left",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0)
    )


# =============================================================================
# SECTION 17: PARAMETER DISPLAY
# =============================================================================

PARAM_DISPLAY = {
    "N_TURNS             #": N_TURNS,
    "TAU_P               m": TAU_P,
    "W_COIL              m": W_COIL,
    "SLIP_RATIO_NORMAL   %": SLIP_RATIO_NORMAL * 100,
    "SLIP_RATIO_REDUCED  %": SLIP_RATIO_REDUCED * 100,
    "GAP                mm": GAP * 1000,
    "HTS_TAPE_WIDTH     mm": HTS_TAPE_WIDTH_MM,
    "HTS_TAPE_LAYERS     #": HTS_TAPE_LAYERS,
    "I_C                 A": round(I_C),
    "I_PEAK              A": round(I_PEAK),
    "I_TARGET            A": round(I_TARGET),
    "V_SLIP_MAX        m/s": V_SLIP_MAX,
    "V_SLIP_MIN        m/s": V_SLIP_MIN,
    "LIMS_PER_SIDE       #": LIMS_PER_SIDE,
    "LIMS_PER_SITE       #": LIMS_PER_SITE,
    "PITCH_COUNT         #": PITCH_COUNT,
    "LIM_SPACING         m": LIM_SPACING,
    "LIM_SITES           #": LIM_SITES,
    "THRUST_EFFICIENCY   %": THRUST_EFFICIENCY * 100,
    "VOLTS_MAX          kV": VOLTS_MAX / 1000,
    "MAX_SITE_POWER     MW": MAX_SITE_POWER / 1e6,
    "T_PLATE            mm": T_PLATE * 1000,
    "M_CABLE         kg/m": M_CABLE_M,
    "M_LOAD          kg/m": M_LOAD_M,
}


def print_parameters():
    """Display current parameter configuration."""
    print("\n" + "="*70)
    print("SIMULATION PARAMETERS")
    print("="*70)
    for key, value in PARAM_DISPLAY.items():
        print(f"  {key:24} {value:>12.2f}")
    print("="*70 + "\n")


# =============================================================================
# SECTION 18: MAIN SIMULATION LOOP
# =============================================================================

def run_deployment_simulation(v_slip_init, i_peak_init):
    """Run the deployment simulation.
    
    Returns: (param_string, make_graphs, total_time)
    """
    global THRUST_MODEL
    
    # Print configuration
    print(f"VOLTS_MAX: {VOLTS_MAX/1000:.0f} kV")
    print(f"THRUST_MODEL: {THRUST_MODEL}", end=" ")
    model_names = {
        1: "(Narrow plate eddy current - geometry corrected)",
        2: "(Goodness factor - Laithwaite)", 
        3: "(Slip × pressure - theoretical maximum)"
    }
    print(model_names.get(THRUST_MODEL, "(Unknown)"))
    print(f"T_PLATE: {T_PLATE * 1000:.1f} mm")
    
    # Initialize state
    dt = DT1
    v_cable = V_ORBIT
    v_casing = V_ORBIT
    time = 0
    v_slip = v_slip_init
    i_peak = i_peak_init
    plate_temp = 70.0  # Start at cryogenic temperature
    
    # Controller state
    v_slip_cmd = v_slip
    i_peak_prev = i_peak
    site_power = 0.0
    
    # Controller tuning
    CTRL_TAU = 15 * DAY
    DV_SLIP_MAX = 0.0001
    DI_MAX = 0.1
    
    # Energy accumulators
    E_site = 0.0
    E_total = 0.0
    
    # Tracking
    sample_time = 0
    power_track = []
    # max_values = {}
    
    # Time display
    month_count = 0
    
    param_str = f"τp={TAU_P}m, N={N_TURNS}, Gap={GAP*1000:.0f}mm"
    
    # Main loop
    while v_casing > V_GROUND_STATIONARY:
        
        # Controller iterations for self-consistency
        for _ in range(10):
            v_rel = calc_relative_velocity(v_cable, v_casing)
            
            # Predictive power control
            power_margin = 1.0 - (site_power / MAX_SITE_POWER)
            MARGIN_COMFORTABLE = 0.15
            MARGIN_TIGHT = 0.05
            
            if power_margin >= MARGIN_COMFORTABLE:
                s_target = SLIP_RATIO_NORMAL
            elif power_margin <= MARGIN_TIGHT:
                s_target = SLIP_RATIO_REDUCED
            else:
                t = (power_margin - MARGIN_TIGHT) / (MARGIN_COMFORTABLE - MARGIN_TIGHT)
                s_target = SLIP_RATIO_REDUCED + t * (SLIP_RATIO_NORMAL - SLIP_RATIO_REDUCED)
            
            # Current headroom factor
            current_ratio = i_peak / I_TARGET
            if current_ratio < 0.7:
                blend = min(1.0, max(0.0, (0.7 - current_ratio) / 0.2))
                s_target = s_target + blend * (SLIP_RATIO_REDUCED - s_target)
            
            s_target = max(1e-6, min(0.999, s_target))
            
            # Calculate target slip velocity
            v_slip_target = (s_target / (1.0 - s_target)) * max(v_rel, 0.0)
            v_slip_target = max(V_SLIP_MIN, min(V_SLIP_MAX, v_slip_target))
            
            # Smooth actuator response
            alpha = 1.0 - math.exp(-dt / max(CTRL_TAU, 1e-9))
            v_slip_cmd = v_slip_cmd + alpha * (v_slip_target - v_slip_cmd)
            max_step = DV_SLIP_MAX * dt
            v_slip_cmd = max(v_slip - max_step, min(v_slip + max_step, v_slip_cmd))
            v_slip = v_slip_cmd
            
            # Calculate operating point
            v_wave = v_rel + v_slip
            slip = calc_slip_ratio(v_slip, v_rel)
            f_slip = calc_slip_frequency(v_slip)
            f_supply = calc_supply_frequency(v_wave)
            
            b_plate = calc_b_field_at_plate(i_peak)
            b_coil = calc_b_field_in_coil(i_peak)
            
            thrust = calc_thrust(f_slip, f_supply, i_peak, plate_temp)
            p_thrust = calc_thrust_power(thrust, v_rel)
            
            volts = calc_coil_voltage(i_peak, f_supply, b_coil)
            
            p_eddy = calc_eddy_loss_from_thrust(v_slip, thrust)
            p_hyst = calc_hysteresis_power(f_supply, i_peak)
            p_heat = p_eddy + p_hyst
            
            temp_plate_avg = calc_plate_temperature(v_rel, p_eddy)
            delta_temp = calc_plate_temp_rise(v_rel, f_slip, p_eddy, plate_temp)
            plate_temp = temp_plate_avg + delta_temp
            
            heat_load = calc_heat_load(p_eddy)
            
            if v_rel > 0:
                total_heat = LIMS_PER_SITE * p_heat + Q_ABSORBED_PER_SITE
                p_cryo = calc_cryo_power(total_heat)
            else:
                p_cryo = 0.0
            
            p_lim = (p_thrust + heat_load + p_hyst) / max(LIM_EFF, 1e-6)
            site_power = (LIMS_PER_SITE * p_lim + p_cryo) / max(INV_EFF, 1e-6)
            
            # Enforce limits
            changed = False
            
            if volts > VOLTS_MAX:
                scale = VOLTS_MAX / max(volts, 1.0)
                i_peak = max(I_MIN, i_peak * max(0.0, min(1.0, scale)))
                changed = True
            
            if temp_plate_avg > T_MAX_PLATE or (temp_plate_avg + delta_temp) > T_MAX_PLATE:
                i_peak = max(I_MIN, i_peak * 0.95)
                v_slip = max(V_SLIP_MIN, v_slip * 0.95)
                changed = True
            
            if site_power > MAX_SITE_POWER:
                scale = math.sqrt(MAX_SITE_POWER / max(site_power, 1.0))
                i_new = max(I_MIN, i_peak * max(0.0, min(1.0, scale)))
                if i_new < i_peak:
                    i_peak = i_new
                    changed = True
                if i_peak <= I_MIN * 1.0001 and v_slip > V_SLIP_MIN:
                    v_slip = max(V_SLIP_MIN, v_slip * 0.95)
                    changed = True
            
            # Current slew limit
            i_step = DI_MAX * dt
            i_peak = max(i_peak_prev - i_step, min(i_peak_prev + i_step, i_peak))
            i_peak_prev = i_peak
            
            # Ramp up if under limits
            if (not changed and i_peak < I_TARGET 
                and site_power < POWER_HEADROOM * MAX_SITE_POWER
                and volts < POWER_HEADROOM * VOLTS_MAX
                and temp_plate_avg < POWER_HEADROOM * T_MAX_PLATE):
                i_peak = min(I_TARGET, i_peak * CURRENT_UPRATE)
                changed = True
            
            if not changed:
                break
        
        # Apply acceleration
        v_cable = calc_cable_velocity(v_cable, thrust, dt)
        v_casing = calc_casing_velocity(v_casing, thrust, dt)
        
        # Accumulate energy
        E_site += p_thrust * LIMS_PER_SITE * dt
        E_total += p_thrust * LIM_SITES * LIMS_PER_SITE * dt
        
        # Startup ramp
        if time < HR:
            if (temp_plate_avg < T_MAX_PLATE * 0.8 
                and volts < VOLTS_MAX * 0.8 
                and site_power < MAX_SITE_POWER * 0.8):
                i_peak += (I_TARGET - i_peak) * 0.01
        
        # Data collection
        skin_eff = calc_effective_plate_depth(f_slip, plate_temp)
        skin_calc = calc_skin_depth(f_slip, plate_temp)
        
        if time > sample_time and time < SAMPLE_TIME_MAX:
            sample_time += SKIP
            data_current.append(i_peak)
            data_voltage.append(volts)
            data_v_slip.append(v_slip)
            data_slip_ratio.append(slip * 100)
            data_f_slip.append(f_slip)
            data_thrust.append(thrust)
            data_thrust_power.append(p_thrust)
            data_v_rel.append(v_rel)
            data_b_field.append(b_plate)
            data_p_eddy.append(p_eddy)
            data_skin_depth_eff.append(skin_eff * 1000)
            data_skin_depth_calc.append(skin_calc * 1000)
            data_p_cryo.append(p_cryo)
            data_p_hyst.append(p_hyst)
            data_p_lim.append(p_lim)
            data_p_site.append(site_power)
            data_temp_plate.append(temp_plate_avg)
            data_E_site.append(E_site)
            data_E_total.append(E_total)
        
        # Monthly progress display
        if time > month_count * MONTH:
            progress = calc_deployment_progress(v_casing)
            if month_count == 0:
                print(f"\nMonth | Progress | Voltage | Current | Thrust | Site Power")
                print("-" * 65)
            print(f"{month_count:5} | {progress:6.1f}% | {volts:7.0f} V | {i_peak:5.1f} A | {thrust:6.0f} N | {site_power/1e6:5.2f} MW")
            
            power_track.append([
                f"{month_count} mth", round(v_casing, 1), round(i_peak, 1),
                round(volts, 1), round(v_slip, 1), f"{round(slip*100, 1)}%",
                round(f_supply, 1), round(p_eddy, 1), round(p_hyst, 1),
                round(thrust, 1), f"{round(p_thrust/1e6, 3)} MW",
                f"{round(p_cryo/1e6, 3)} MW", f"{round(site_power/1e6, 3)} MW"
            ])
            month_count += 1
        
        # Check for failures
        if volts > VOLTS_MAX * 1.1:
            print(f"\nFAILURE: Voltage exceeded limit ({volts:.0f} V)")
            return (param_str, False, time)
        
        if site_power > MAX_SITE_POWER * 1.1:
            print(f"\nFAILURE: Site power exceeded limit ({site_power/1e6:.1f} MW)")
            return (param_str, False, time)
        
        if TAU_P * PITCH_COUNT > LIM_SPACING:
            print(f"\nFAILURE: LIM length ({TAU_P * PITCH_COUNT}m) > spacing ({LIM_SPACING}m)")
            return (param_str, False, time)
        
        # Time stepping
        if time < DAY:
            dt = DT1
        elif time < SAMPLE_TIME_MAX:
            dt = DT2
        else:
            dt = DT3
        
        time += dt
        
        if time > 100 * YR:
            print("\nWARNING: Deployment time exceeded 100 years")
            break
    
    # Final results
    print("\n" + "=" * 70)
    print("DEPLOYMENT COMPLETE")
    print("=" * 70)
    print(f"Deployment time: {time/DAY:.1f} days ({time/YR:.2f} years)")
    print(f"Final cable velocity: {v_cable:.2f} m/s")
    print(f"Final casing velocity: {v_casing:.2f} m/s")
    print(f"Total energy delivered: {E_total/1e18:.3f} EJ")
    print(f"Energy per site: {E_site/1e12:.2f} TJ")
    print("=" * 70)
    
    return (param_str, True, time)


# =============================================================================
# SECTION 19: PLOTTING FUNCTIONS
# =============================================================================

def plot_results(show_graphs, param_str, total_time):
    """Generate requested plots."""
    
    plots = {
        "current": (data_current, "Current", "A", "blue"),
        "volts": (data_voltage, "Voltage", "V", "red"),
        "v_slip": (data_v_slip, "Slip Velocity", "m/s", "green"),
        "thrust": (data_thrust, "Thrust", "N", "purple"),
        "p_eddy": (data_p_eddy, "Eddy Current Losses", "W", "darkblue"),
        "v_rel": (data_v_rel, "Relative Velocity", "m/s", "olive"),
        "f_slip": (data_f_slip, "Slip Frequency", "Hz", "orange"),
        "slip": (data_slip_ratio, "Slip Ratio", "%", "cyan"),
        "p_thrust": (data_thrust_power, "Thrust Power", "W", "navy"),
        "b_peak": (data_b_field, "Magnetic Field", "T", "brown"),
        "hyst": (data_p_hyst, "Hysteresis Losses", "W", "brown"),
        "cryo": (data_p_cryo, "Cryogenic Power", "W", "teal"),
        "power": (data_p_site, "Site Power", "W", "darkslategray"),
        "lim_power": (data_p_lim, "LIM Power", "W", "darkslategray"),
        "plate_temp": (data_temp_plate, "Plate Temperature", "K", "lime"),
        "ke_site": (data_E_site, "Site Kinetic Energy", "J", "green"),
        "ke_all": (data_E_total, "Total Kinetic Energy", "J", "darkgreen"),
        "skin": (data_skin_depth_eff, "Skin Depth", "mm", "magenta"),
    }
    
    for name, (data, label, unit, color) in plots.items():
        if name in show_graphs or "all" in show_graphs:
            if not data:
                continue
            
            tick_pos, tick_labels = make_month_ticks(data, total_time)
            
            plt.figure(figsize=(10, 6))
            plt.scatter(range(len(data)), data, c=color, s=1)
            plt.xlabel("Months")
            plt.ylabel(f"{label} ({unit})")
            plt.title(f"{label} - {param_str}")
            plt.xticks(tick_pos, tick_labels)
            annotate_final_value(data, unit=unit)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()


# =============================================================================
# SECTION 20: MAIN ENTRY POINT
# =============================================================================

def main():
    """Main entry point."""
    global THRUST_MODEL
    
    # Parse command line
    show_graphs = []
    
    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            print("""
Graph Options:
  all          Show all graphs
  current      Coil current (A)
  volts        Induced voltage (V)
  v_slip       Slip velocity (m/s)
  thrust       Thrust force (N)
  p_thrust     Thrust power (W)
  p_eddy       Eddy current losses (W)
  power        Total site power (W)
  plate_temp   Reaction plate temperature (K)
  skin         Skin depth (mm)
  slip         Slip ratio (%)
  f_slip       Slip frequency (Hz)
  v_rel        Relative velocity (m/s)
  hyst         Hysteresis losses (W)
  cryo         Cryogenic power (W)
  ke_site      Site kinetic energy (J)
  ke_all       Total kinetic energy (J)

Examples:
  python lim_deployment_sim.py --model=1 thrust power
  python lim_deployment_sim.py --model=2 all
  python lim_deployment_sim.py --model=3
            """)
            return
        
        # Parse model selection
        for arg in sys.argv[1:]:
            if arg.startswith("--model="):
                try:
                    model = int(arg.split("=")[1])
                    if model in [1, 2, 3]:
                        THRUST_MODEL = model
                    else:
                        print(f"Invalid model {model}, using default (1)")
                except ValueError:
                    print(f"Invalid model argument: {arg}")
            else:
                show_graphs.append(arg)
    
    # Print configuration
    print_parameters()
    print(f"Using Thrust Model {THRUST_MODEL}")
    
    # Run simulation
    v_slip = V_SLIP_MAX
    i_peak = I_MIN
    
    param_str, success, total_time = run_deployment_simulation(v_slip, i_peak)
    
    # Plot results if requested
    if success and show_graphs:
        plot_results(show_graphs, param_str, total_time)


if __name__ == "__main__":
    main()
