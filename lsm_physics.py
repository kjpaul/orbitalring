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

    As the sled moves at velocity v, its DC field B_sled sweeps past the
    stationary stator coils, inducing a back-EMF:

    EMF = N × 2π × v × B_sled × w_coil

    Key insight: This is INDEPENDENT of pole pitch τ_p. The pole pitch cancels
    when you substitute f = v/(2τ_p) and Φ = B_sled × 2τ_p × w_coil.

    This EMF grows linearly with velocity and is the primary voltage constraint
    at high speed.

    Args:
        N: Turns per stator coil
        v: Sled velocity (m/s)
        B_sled: Sled DC field at air gap (T)
        w_coil: Coil height (m)

    Returns:
        Back-EMF per coil (V)
    """
    return N * 2.0 * PI * v * B_sled * w_coil


def calc_max_B_sled_for_voltage(V_limit, N, v, w_coil):
    """Calculate maximum B_sled allowed by voltage limit at given velocity.

    At high velocities, the back-EMF can exceed the coil insulation limit.
    This function calculates the maximum sled field that keeps the voltage
    within limits:

    B_sled_max = V_limit / (N × 2π × v × w_coil)

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

    return V_limit / (N * 2.0 * PI * v * w_coil)


# =============================================================================
# THRUST CALCULATIONS
# =============================================================================

def calc_thrust(B_stator, B_sled, delta, A_active):
    """Calculate LSM thrust from Maxwell stress.

    The thrust comes from the magnetic shear stress at the air gap:

    σ_shear = (B_stator × B_sled × sin(δ)) / μ₀
    F = σ_shear × A_active

    where:
      - B_stator = stator AC field amplitude (T)
      - B_sled = sled DC field amplitude (T)
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

    sigma_shear = (B_stator * B_sled * math.sin(delta)) / MU0
    F = sigma_shear * A_active
    return F


def calc_delta_for_thrust(F_desired, B_stator, B_sled, A_active, delta_max):
    """Calculate required load angle to achieve desired thrust.

    Inverts the thrust equation to find the load angle:

    sin(δ) = F_desired × μ₀ / (B_stator × B_sled × A_active)

    If the required sin(δ) exceeds sin(delta_max), the thrust is limited.

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

    sin_delta_required = (F_desired * MU0) / (B_stator * B_sled * A_active)

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
# HTS HYSTERESIS LOSSES
#
# The LSM stator uses HTS coils at 77 K. As the sled passes, the AC
# magnetization cycle generates hysteresis losses in the superconducting
# tape. These must be pumped out by the cryogenic system.
# =============================================================================

def q_hts_loss_factor(i_peak, n_turns, w_coil, w_tape, alpha_tape):
    """Hysteresis loss factor for HTS tape (Bean model approximation).

    Returns energy loss per meter of tape per cycle (J/m/cycle).
    """
    b_coil = MU0 * n_turns * i_peak / w_coil
    return b_coil * i_peak * w_tape * math.sin(alpha_tape)


def q_hysteresis_norris(i_now, i_c):
    """Norris strip formula for hysteresis loss per meter per cycle.

    More accurate than loss-factor model at high I/Ic ratios.
    Returns energy loss per meter of tape per cycle (J/m/cycle).
    """
    ii = i_now / i_c
    if ii >= 1:
        ii = 0.999
    if ii <= -1:
        ii = -0.999
    return (MU0 * i_c**2 / math.pi) * (
        (1 - ii) * math.log(1 - ii) + (1 + ii) * math.log(1 + ii) - ii**2
    )


def calc_hysteresis_power(f_supply, I_stator, params, n_units, norris=False):
    """Total HTS hysteresis power for all active stator coils.

    Args:
        f_supply: Electrical supply frequency (Hz)
        I_stator: Stator operating current (A)
        params: Physics parameter dict from lsm_config.get_physics_params()
        n_units: Number of repeating units on the sled
        norris: Use Norris formula if True, else loss-factor model

    Returns:
        Total hysteresis power (W) for all active coils
    """
    if norris:
        q = q_hysteresis_norris(I_stator, params['i_c'])
    else:
        q = q_hts_loss_factor(I_stator, params['n_turns'], params['w_coil'],
                              params['w_tape'], params['alpha_tape'])

    # Power per phase coil = q (J/m/cycle) × tape_length (m) × frequency (Hz)
    p_coil = q * params['l_hts_coil'] * f_supply

    # Total: phases × sides × n_units
    n_coils = params['lim_phases'] * params['n_lim_sides'] * n_units
    return p_coil * n_coils


# =============================================================================
# CRYOGENIC RADIATOR WIDTH
# =============================================================================

def calc_radiator_width(q_cold, T_cold, T_hot, efficiency, em_heatsink,
                        radiator_length, T_space):
    """Required cryogenic radiator width to reject hysteresis heat.

    The cryo system pumps heat from T_cold (77 K stator) to T_hot
    (400 K radiator). The radiator must reject both the cold-side heat
    AND the compressor work.

    Args:
        q_cold: Cold-side heat load (W)
        T_cold: Cold temperature (K), typically 77 K
        T_hot: Hot-side radiator temperature (K)
        efficiency: Fraction of Carnot COP achieved
        em_heatsink: Radiator emissivity
        radiator_length: Length of radiator along ring (m)
        T_space: Background space temperature (K)

    Returns:
        Required radiator width (m) perpendicular to ring
    """
    if q_cold <= 0:
        return 0.0
    if T_hot <= T_cold:
        T_hot = T_cold + 1

    cop_carnot = T_cold / (T_hot - T_cold)
    cop_real = cop_carnot * efficiency
    if cop_real <= 0:
        cop_real = 0.01

    q_reject = q_cold * (1 + 1 / cop_real)
    p_per_m2 = em_heatsink * STEFAN_BOLTZMANN * (T_hot**4 - T_space**4)

    if p_per_m2 <= 0:
        return 1000.0

    area_required = q_reject / p_per_m2
    return area_required / radiator_length
