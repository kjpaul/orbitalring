#!/usr/bin/env python3
"""
LSM Physics Module - Calculation functions for LSM mass driver simulation

This module contains all the physics calculations for the LSM (Linear Synchronous Motor):
  - Magnetic field calculations (air-core)
  - Thrust from Maxwell stress
  - Back-EMF and voltage constraints
  - Power calculations
  - Orbital dynamics and occupant g-loads

Key differences from LIM physics:
  - No slip (synchronous operation)
  - No eddy current losses
  - No thermal calculations for sled
  - Direct magnetic interaction (not induced currents)
  - Voltage limit is the primary constraint at high speed

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7    # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.670374e-8  # Stefan-Boltzmann constant (W/m²K⁴)
PI = math.pi

# =============================================================================
# STATOR MAGNETIC FIELD
# =============================================================================

def calc_B_stator_air_core(I, N, w_coil, g_gap):
    """Calculate stator magnetic field at the air gap (air-core stator).

    For an air-core stator, the field is calculated using the arctangent formula
    for a finite-width coil:

    B = (2/π) × (μ₀ × N × I / w_coil) × arctan(w_coil / (2 × g_gap))

    Args:
        I: Stator current (A)
        N: Turns per stator coil
        w_coil: Coil height/width (m)
        g_gap: Effective air gap (m)

    Returns:
        Stator field at air gap (T)
    """
    if I <= 0:
        return 0.0

    B = (2.0 / PI) * (MU0 * N * I / w_coil) * math.atan(w_coil / (2.0 * g_gap))
    return B


def calc_B_stator(I, N, w_coil, g_gap):
    """Calculate stator magnetic field at the air gap (air-core).

    Dispatches to air-core formula.

    Args:
        I: Stator current (A)
        N: Turns per stator coil
        w_coil: Coil height (m)
        g_gap: Effective air gap (m)

    Returns:
        Stator field at air gap (T)
    """
    return calc_B_stator_air_core(I, N, w_coil, g_gap)


# =============================================================================
# BACK-EMF AND VOLTAGE
# =============================================================================

def calc_EMF(N, v, B_sled, w_coil):
    """Calculate back-EMF per stator coil from moving sled field.

    Derivation from Faraday's law:
      Each stator coil spans one pole pitch τ_p with N turns and width w_coil.
      The sled's sinusoidal field B₀sin(πx/τ_p) gives peak flux linkage:
        Λ_peak = N × w_coil × (2τ_p/π) × B₀
      With ω = πv/τ_p, taking dΛ/dt:
        EMF_peak = ω × Λ_peak = (πv/τ_p) × N × w_coil × (2τ_p/π) × B₀
                 = 2 × N × v × B₀ × w_coil
      The τ_p cancels. The coefficient is 2, not 2π.

    This is INDEPENDENT of pole pitch τ_p and grows linearly with velocity,
    making it the primary voltage constraint at high speed.

    Args:
        N: Turns per stator coil
        v: Sled velocity (m/s)
        B_sled: Sled DC field at air gap (T)
        w_coil: Coil height (m)

    Returns:
        Back-EMF per coil (V)
    """
    return 2.0 * N * v * B_sled * w_coil


def calc_max_B_sled_for_voltage(V_limit, N, v, w_coil):
    """Calculate maximum B_sled allowed by voltage limit at given velocity.

    Inverts the EMF formula: B_sled_max = V_limit / (2 × N × v × w_coil)

    Args:
        V_limit: Coil voltage limit (V)
        N: Turns per stator coil
        v: Sled velocity (m/s)
        w_coil: Coil height (m)

    Returns:
        Maximum B_sled (T), or inf if v is zero
    """
    if v <= 0:
        return float('inf')

    return V_limit / (2.0 * N * v * w_coil)


# =============================================================================
# THRUST CALCULATIONS
# =============================================================================

def calc_thrust(B_stator, B_sled, delta, A_active):
    """Calculate LSM thrust from Maxwell stress with spatial averaging.

    When both B_stator and B_sled are peak fundamental amplitudes of
    sinusoidal field distributions, the time-and-space-averaged shear
    stress includes a factor of 1/2 from averaging sin²(kx) over a
    pole pitch:

    σ_shear = (B_stator × B_sled × sin(δ)) / (2 × μ₀)
    F = σ_shear × A_active

    where:
      - B_stator = peak stator AC field amplitude (T)
      - B_sled = peak sled DC field amplitude (T)
      - δ = load angle (electrical angle between stator and sled fields)
      - A_active = total active magnetic area (m²)

    Maximum thrust occurs at δ = π/2 (90°), but this is the pull-out limit.
    Safe operation is typically δ < π/3 (60°).

    Args:
        B_stator: Stator field (T)
        B_sled: Sled field (T)
        delta: Load angle (radians)
        A_active: Total active area, both sides (m²)

    Returns:
        Thrust force (N)
    """
    if B_stator <= 0 or B_sled <= 0 or A_active <= 0:
        return 0.0

    sigma_shear = (B_stator * B_sled * math.sin(delta)) / (2.0 * MU0)
    F = sigma_shear * A_active
    return F


def calc_delta_for_thrust(F_desired, B_stator, B_sled, A_active, delta_max):
    """Calculate required load angle to achieve desired thrust.

    Inverts the spatially-averaged thrust equation:

    sin(δ) = (F_desired × 2 × μ₀) / (B_stator × B_sled × A_active)

    The factor of 2 in the numerator matches the 1/2 in the thrust
    denominator from spatial averaging of sinusoidal fields.

    Args:
        F_desired: Desired thrust (N)
        B_stator: Stator field (T)
        B_sled: Sled field (T)
        A_active: Total active area (m²)
        delta_max: Maximum allowed load angle (rad)

    Returns:
        (delta, limited) where:
          - delta = required load angle (rad)
          - limited = True if thrust is limited by delta_max
    """
    if B_stator <= 0 or B_sled <= 0 or A_active <= 0:
        return (0.0, True)

    sin_delta_required = (F_desired * 2.0 * MU0) / (B_stator * B_sled * A_active)

    sin_delta_max = math.sin(delta_max)

    if sin_delta_required > sin_delta_max:
        # Thrust limited
        return (delta_max, True)
    else:
        # Can achieve desired thrust
        delta = math.asin(max(0.0, min(1.0, sin_delta_required)))
        return (delta, False)


# =============================================================================
# SUPPLY FREQUENCY
# =============================================================================

def calc_frequency(v, tau_p):
    """Calculate synchronous supply frequency.

    For synchronous operation (no slip), the electrical frequency must match
    the mechanical frequency:

    f = v / (2 × τ_p)

    Args:
        v: Sled velocity (m/s)
        tau_p: Pole pitch (m)

    Returns:
        Supply frequency (Hz)
    """
    if tau_p <= 0:
        return 0.0

    return v / (2.0 * tau_p)


# =============================================================================
# POWER CALCULATIONS
# =============================================================================

def calc_power_mechanical(F, v):
    """Calculate mechanical power delivered to sled.

    For synchronous operation, there is no slip, so:

    P = F × v

    This is the power transferred to the sled's kinetic energy. Unlike the LIM,
    there are no slip losses (eddy current heating in the reaction plate).

    Args:
        F: Thrust force (N)
        v: Sled velocity (m/s)

    Returns:
        Mechanical power (W)
    """
    return F * v


# =============================================================================
# ORBITAL DYNAMICS AND G-LOADS
# =============================================================================

def calc_centrifugal_acceleration(v, R_orbit):
    """Calculate centrifugal acceleration.

    a_c = v² / R_orbit

    Args:
        v: Velocity (m/s)
        R_orbit: Orbital radius (m)

    Returns:
        Centrifugal acceleration (m/s²)
    """
    return v**2 / R_orbit


def calc_net_radial_acceleration(v, R_orbit, g_250):
    """Calculate net radial acceleration (centrifugal - gravitational).

    a_radial_net = v²/R - g_250

    Positive means net outward (centrifugal exceeds gravity).

    Args:
        v: Velocity (m/s)
        R_orbit: Orbital radius (m)
        g_250: Gravitational acceleration at 250 km (m/s²)

    Returns:
        Net radial acceleration (m/s²), positive = outward
    """
    a_c = calc_centrifugal_acceleration(v, R_orbit)
    return a_c - g_250


def calc_tangential_acceleration(F, m_total):
    """Calculate tangential acceleration from thrust.

    a_t = F / m_total

    Args:
        F: Thrust force (N)
        m_total: Total mass (kg)

    Returns:
        Tangential acceleration (m/s²)
    """
    if m_total <= 0:
        return 0.0

    return F / m_total


def calc_occupant_g_load(a_tangential, a_radial_net, g_0):
    """Calculate total occupant g-load (magnitude).

    The occupant experiences both tangential and radial accelerations.
    The total magnitude is:

    a_total = sqrt(a_tangential² + a_radial_net²)

    Expressed in g:

    g_load = a_total / g_0

    Args:
        a_tangential: Tangential acceleration (m/s²)
        a_radial_net: Net radial acceleration (m/s²)
        g_0: Standard gravity (m/s²)

    Returns:
        Occupant g-load (dimensionless, in units of g₀)
    """
    a_total = math.sqrt(a_tangential**2 + a_radial_net**2)
    return a_total / g_0


# =============================================================================
# KINETIC ENERGY
# =============================================================================

def calc_kinetic_energy(m, v):
    """Calculate kinetic energy.

    KE = ½ m v²

    Args:
        m: Mass (kg)
        v: Velocity (m/s)

    Returns:
        Kinetic energy (J)
    """
    return 0.5 * m * v**2


# =============================================================================
# HTS HYSTERESIS LOSSES — GÖMÖRY PERPENDICULAR FIELD MODEL
#
# The LSM stator uses HTS coils at 77 K. As the sled passes, the AC
# magnetization cycle generates hysteresis losses in the superconducting
# tape. The Gömöry model calculates:
#
#   Q_hyst = B_coil_peak × Ic_per_layer × w_tape × sin(α_tape_field)
#
# where B_coil_peak = μ₀ × N × I / l_coil is the field inside the winding.
# =============================================================================

def calc_hysteresis_power(f_supply, I_stator, params, n_units):
    """Total HTS hysteresis power using Gömöry perpendicular field model.

    Args:
        f_supply: Electrical supply frequency (Hz)
        I_stator: Stator operating current (A)
        params: Physics parameter dict from lsm_config.get_physics_params()
        n_units: Number of repeating units on the sled

    Returns:
        Total hysteresis power (W) for all active coils
    """
    n_turns = params['n_turns']
    tau_p = params['tau_p']
    w_coil = params['w_coil']

    # Coil self-field (field inside the winding)
    l_coil = tau_p / 3.0
    B_coil_peak = MU0 * n_turns * I_stator / l_coil

    # Tape geometry per layer
    coil_depth = n_turns * 0.15e-3          # winding depth (m)
    turn_perimeter = 2.0 * (w_coil + coil_depth)
    L_tape_per_layer = n_turns * turn_perimeter

    # Gömöry loss per meter of tape per cycle (J/m/cycle)
    w_tape = params['tape_width_mm'] * 1e-3
    Q_hyst = B_coil_peak * params['i_c_per_layer'] * w_tape * params['sin_alpha']

    # Power per phase coil = layers × tape_length × Q_hyst × frequency
    P_hts_per_coil = params['n_hts_layers'] * L_tape_per_layer * Q_hyst * f_supply

    # Total: phases × sides × n_units
    n_coils_active = params['lim_phases'] * params['n_lim_sides'] * n_units
    return P_hts_per_coil * n_coils_active


# =============================================================================
# CRYOGENIC RADIATOR WIDTH
#
# LN2 absorbs HTS hysteresis heat at 77 K; radiators at ~100 K reject
# it to space (~3 K background). Two-sided radiator panels run along
# the stator section.
# =============================================================================

def calc_radiator_width(P_hts_total, L_stator_active):
    """Required radiator width for LN2-cooled HTS heat rejection.

    Args:
        P_hts_total: Total HTS hysteresis power to reject (W)
        L_stator_active: Length of active stator section (m)

    Returns:
        Required radiator width (m) perpendicular to ring
    """
    if P_hts_total <= 0 or L_stator_active <= 0:
        return 0.0

    epsilon = 0.9       # radiator emissivity
    sigma = 5.67e-8     # Stefan-Boltzmann
    T_rad = 100.0       # K, radiator temperature (slightly above LN2)
    q_rad = epsilon * sigma * T_rad**4   # W/m² per side
    # Two-sided radiator:
    return P_hts_total / (2.0 * q_rad * L_stator_active)
