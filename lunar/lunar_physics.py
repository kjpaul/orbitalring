#!/usr/bin/env python3
"""
Lunar Mass Driver Physics — Calculation functions for lunar surface LSM

This module contains all physics calculations for the lunar mass driver:
  - Motor physics (simplified three-regime calibration model)
  - HTS hysteresis losses (Gömöry model)
  - Lunar centrifugal and rail force calculations
  - G-load calculations
  - Celestial mechanics (Hohmann transfers, Moon phase effects, Kepler)

Motor physics use the calibration model: F/L scales linearly with N and B_sled
from a reference point (F_PER_M_CAL = 1,500 N/m at N=10, B=0.10, I_TARGET).
The power per metre P/L = 37.5 MW/m is invariant with N and B_sled.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7    # Permeability of free space (H/m)
PI = math.pi

# =============================================================================
# MOTOR PHYSICS — CALIBRATION MODEL
# =============================================================================

def thrust_per_m(N, B_sled, f_per_m_cal=1500.0, n_cal=10, b_cal=0.10):
    """Thrust per metre of sled (linear scaling from calibration point).

    F/L = F_cal × (N / N_cal) × (B_sled / B_cal)

    At the reference point (N=10, B=0.10, I_TARGET):
      F/L = 1,500 N/m → F_total = 15 MN for 10 km sled

    At N=50 (default):
      F/L = 7,500 N/m → F_total = 75 MN for 10 km sled

    Args:
        N: Stator turns
        B_sled: Sled DC field at air gap (T)
        f_per_m_cal: Calibration thrust per metre (N/m)
        n_cal: Calibration turn count
        b_cal: Calibration sled field (T)

    Returns:
        Thrust per metre of sled (N/m)
    """
    return f_per_m_cal * (N / n_cal) * (B_sled / b_cal)


def v_crossover(N, B_sled, v_max=100_000, w_coil=2.0):
    """Voltage crossover velocity where back-EMF = V_max.

    v_cross = V_max / (2 × N × B_sled × w_coil)

    Above this velocity, the motor transitions from constant-thrust to
    constant-power operation (voltage-limited regime).

    Args:
        N: Stator turns
        B_sled: Sled DC field (T)
        v_max: Coil voltage limit (V)
        w_coil: Coil width (m)

    Returns:
        Crossover velocity (m/s)
    """
    return v_max / (2.0 * N * B_sled * w_coil)


def power_per_m(N, B_sled, f_per_m_cal=1500.0, n_cal=10, b_cal=0.10,
                v_max=100_000, w_coil=2.0):
    """Power per metre at crossover (F/L × v_cross).

    This is INVARIANT with N and B_sled:
      P/L = F_cal/N_cal/B_cal × V_max / (2 × w_coil)
          = 1500/10/0.10 × 100,000 / (2 × 2.0)
          = 150,000 × 25,000 = 37,500,000 = 37.5 MW/m

    Args:
        N: Stator turns
        B_sled: Sled DC field (T)
        (other args: calibration and voltage parameters)

    Returns:
        Power per metre of sled (W/m)
    """
    return thrust_per_m(N, B_sled, f_per_m_cal, n_cal, b_cal) * \
           v_crossover(N, B_sled, v_max, w_coil)


# =============================================================================
# HTS HYSTERESIS LOSSES — GÖMÖRY MODEL
# =============================================================================

def calc_hysteresis_power(f_supply, I_stator, params, n_units):
    """Total HTS hysteresis power using Gömöry perpendicular field model.

    Identical to the orbital ring calculation. The stator HTS coils
    experience AC magnetization as the sled passes, generating hysteresis
    losses proportional to frequency.

    Args:
        f_supply: Electrical supply frequency (Hz)
        I_stator: Stator operating current (A)
        params: Physics parameter dict from lunar_config.get_physics_params()
        n_units: Number of repeating units on the sled

    Returns:
        Total hysteresis power (W) for all active coils
    """
    n_turns = params['n_turns']
    tau_p = params['tau_p']
    w_coil = params['w_coil']

    # Coil self-field
    l_coil = tau_p / 3.0
    B_coil_peak = MU0 * n_turns * I_stator / l_coil

    # Tape geometry per layer
    coil_depth = n_turns * 0.15e-3
    turn_perimeter = 2.0 * (w_coil + coil_depth)
    L_tape_per_layer = n_turns * turn_perimeter

    # Gömöry loss per meter of tape per cycle (J/m/cycle)
    w_tape = params['tape_width_mm'] * 1e-3
    Q_hyst = B_coil_peak * params['i_c_per_layer'] * w_tape * params['sin_alpha']

    # Power per phase coil
    P_hts_per_coil = params['n_hts_layers'] * L_tape_per_layer * Q_hyst * f_supply

    # Total: phases × sides × n_units
    n_coils_active = params['n_phases'] * params['n_stator_sides'] * n_units
    return P_hts_per_coil * n_coils_active


def calc_supply_frequency(v, tau_p):
    """Synchronous supply frequency: f = v / (2 × τ_p).

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
# LUNAR CENTRIFUGAL AND RAIL FORCES
# =============================================================================

def calc_centrifugal_acceleration(v, R_moon):
    """Centrifugal acceleration at velocity v on lunar surface.

    a_cent = v² / R_moon

    At v_orbital (~1,679 m/s), a_cent = g_moon.

    Args:
        v: Sled velocity (m/s)
        R_moon: Lunar radius (m)

    Returns:
        Centrifugal acceleration (m/s²)
    """
    return v**2 / R_moon


def calc_net_outward_acceleration(v, R_moon, g_moon):
    """Net outward (centrifugal - gravitational) acceleration.

    a_out = v²/R - g_moon

    Positive means centrifugal exceeds gravity (above v_orbital).
    Negative means gravity exceeds centrifugal (below v_orbital).

    Args:
        v: Sled velocity (m/s)
        R_moon: Lunar radius (m)
        g_moon: Lunar surface gravity (m/s²)

    Returns:
        Net outward acceleration (m/s²), positive = outward
    """
    return v**2 / R_moon - g_moon


def calc_rail_force_per_m(v, R_moon, g_moon, m_per_m_total):
    """Outward rail force per metre of sled (for structural loading).

    Below v_orbital: force is inward (gravity holds sled on track, returns 0).
    Above v_orbital: force is outward (rails must hold sled down).

    F_rail/m = m_per_m_total × max(0, v²/R - g)

    where m_per_m_total includes both sled hardware and distributed payload mass.

    Args:
        v: Sled velocity (m/s)
        R_moon: Lunar radius (m)
        g_moon: Lunar surface gravity (m/s²)
        m_per_m_total: Total mass per metre of sled (kg/m)

    Returns:
        Outward rail force per metre (N/m), zero if below v_orbital
    """
    a_net = v**2 / R_moon - g_moon
    if a_net <= 0:
        return 0.0
    return m_per_m_total * a_net


def calc_effective_gravity(v, R_moon, g_moon):
    """Effective gravity felt by sled (gravity minus centrifugal).

    g_eff = g_moon - v²/R

    Positive means net downward (sled stays on track by gravity).
    Negative means net outward (rails needed).

    Args:
        v: Sled velocity (m/s)
        R_moon: Lunar radius (m)
        g_moon: Lunar surface gravity (m/s²)

    Returns:
        Effective gravity (m/s²), positive = downward
    """
    return g_moon - v**2 / R_moon


# =============================================================================
# G-LOAD CALCULATIONS
# =============================================================================

def calc_g_load(a_tangential, a_net_outward, g_0):
    """Total g-load magnitude from tangential and radial accelerations.

    g_load = sqrt(a_tangential² + a_net_outward²) / g_0

    On the Moon, the net outward acceleration is small at typical launch
    velocities (a few km/s), so g-load is dominated by tangential acceleration.

    Args:
        a_tangential: Tangential (thrust) acceleration (m/s²)
        a_net_outward: Net outward acceleration (m/s²)
        g_0: Standard gravity (m/s²)

    Returns:
        Total g-load (dimensionless, in units of g₀)
    """
    a_total = math.sqrt(a_tangential**2 + a_net_outward**2)
    return a_total / g_0


# =============================================================================
# CELESTIAL MECHANICS — HOHMANN TRANSFERS
# =============================================================================

def hohmann_v_inf(a_departure, a_arrival, mu_sun=1.32712440018e20):
    """Calculate v_inf at departure for a Hohmann transfer orbit.

    For a transfer from inner to outer planet:
      a_transfer = (a_departure + a_arrival) / 2
      v_dep = sqrt(μ_sun × (2/a_dep - 1/a_t))
      v_circ = sqrt(μ_sun / a_dep)
      v_inf = v_dep - v_circ

    Args:
        a_departure: Departure planet semi-major axis (m)
        a_arrival: Arrival planet semi-major axis (m)
        mu_sun: Sun gravitational parameter (m³/s²)

    Returns:
        v_inf at departure (m/s) — excess velocity over circular orbit
    """
    a_t = (a_departure + a_arrival) / 2.0
    v_dep = math.sqrt(mu_sun * (2.0 / a_departure - 1.0 / a_t))
    v_circ = math.sqrt(mu_sun / a_departure)
    return v_dep - v_circ


def hohmann_transfer_time(a_departure, a_arrival, mu_sun=1.32712440018e20):
    """Transfer time for a Hohmann orbit (half the ellipse period).

    T = π × sqrt(a_transfer³ / μ_sun)

    Args:
        a_departure: Departure planet semi-major axis (m)
        a_arrival: Arrival planet semi-major axis (m)
        mu_sun: Sun gravitational parameter (m³/s²)

    Returns:
        Transfer time (s)
    """
    a_t = (a_departure + a_arrival) / 2.0
    return PI * math.sqrt(a_t**3 / mu_sun)


# =============================================================================
# CELESTIAL MECHANICS — MOON PHASE AND LAUNCH VELOCITY
# =============================================================================

def v_at_moon_orbit(v_inf_earth, v_esc_earth_at_moon):
    """Speed at Moon's orbital distance needed for given v_inf at Earth.

    v_at_moon = sqrt(v_inf_earth² + v_esc_earth²)

    This is the vis-viva equation at Moon's distance from Earth.

    Args:
        v_inf_earth: Required v_inf at Earth (m/s)
        v_esc_earth_at_moon: Escape velocity from Earth at Moon's distance (m/s)

    Returns:
        Required speed relative to Earth at Moon's orbit (m/s)
    """
    return math.sqrt(v_inf_earth**2 + v_esc_earth_at_moon**2)


def v_inf_moon_for_phase(v_inf_earth, phase_deg,
                         v_moon_orbital=1018.0,
                         v_esc_earth_at_moon=1440.0):
    """Required v_inf at Moon's SOI for a given Earth departure v_inf and Moon phase.

    The Moon orbits Earth at v_moon. Depending on the Moon's orbital phase
    relative to the desired departure direction:
      - Phase 0° (prograde): Moon's velocity adds → minimum v_inf_moon needed
      - Phase 180° (retrograde): Moon's velocity opposes → maximum v_inf_moon

    Solved from the vector addition:
      v_at_moon² = v_moon² + v_inf_moon² + 2 × v_moon × v_inf_moon × cos(φ)

    Args:
        v_inf_earth: Required v_inf at Earth (m/s)
        phase_deg: Moon orbital phase (0° = optimal/prograde, 180° = worst)
        v_moon_orbital: Moon's orbital velocity around Earth (m/s)
        v_esc_earth_at_moon: Escape velocity from Earth at Moon's distance (m/s)

    Returns:
        Required v_inf at Moon's SOI (m/s)
    """
    v_needed = v_at_moon_orbit(v_inf_earth, v_esc_earth_at_moon)
    phase_rad = math.radians(phase_deg)
    cos_phi = math.cos(phase_rad)

    # Quadratic: v_inf² + 2*v_moon*cos(φ)*v_inf + (v_moon² - v_needed²) = 0
    a = 1.0
    b = 2.0 * v_moon_orbital * cos_phi
    c = v_moon_orbital**2 - v_needed**2

    discriminant = b**2 - 4.0 * a * c
    if discriminant < 0:
        return float('inf')

    v_inf_moon = (-b + math.sqrt(discriminant)) / (2.0 * a)
    return max(0.0, v_inf_moon)


def launch_velocity_from_moon(v_inf_moon, v_escape_moon):
    """Launch velocity from Moon's surface for given v_inf at Moon's SOI.

    v_launch = sqrt(v_inf_moon² + v_escape_moon²)

    This is the vis-viva equation at Moon's surface (r = R_moon).

    Args:
        v_inf_moon: Required v_inf at Moon's SOI (m/s)
        v_escape_moon: Lunar escape velocity (m/s)

    Returns:
        Required launch velocity from Moon's surface (m/s)
    """
    return math.sqrt(v_inf_moon**2 + v_escape_moon**2)


def mission_launch_velocity(a_departure, a_arrival, phase_deg,
                            v_escape_moon, v_moon_orbital, v_esc_earth_at_moon,
                            mu_sun=1.32712440018e20):
    """Full chain: planet pair → required launch velocity from Moon surface.

    Combines: Hohmann v_inf → Moon phase → launch velocity

    Args:
        a_departure: Departure planet SMA (m) — usually A_EARTH
        a_arrival: Arrival planet SMA (m)
        phase_deg: Moon orbital phase (0° = optimal)
        v_escape_moon: Lunar escape velocity (m/s)
        v_moon_orbital: Moon's orbital velocity (m/s)
        v_esc_earth_at_moon: Escape velocity from Earth at Moon's distance (m/s)
        mu_sun: Sun gravitational parameter (m³/s²)

    Returns:
        Required launch velocity from Moon's surface (m/s)
    """
    v_inf_earth = hohmann_v_inf(a_departure, a_arrival, mu_sun)
    v_inf_moon = v_inf_moon_for_phase(v_inf_earth, phase_deg,
                                       v_moon_orbital, v_esc_earth_at_moon)
    return launch_velocity_from_moon(v_inf_moon, v_escape_moon)


# =============================================================================
# KEPLER'S EQUATION
# =============================================================================

def solve_kepler(M, e, tol=1e-12, max_iter=100):
    """Solve Kepler's equation M = E - e × sin(E) for eccentric anomaly E.

    Uses Newton-Raphson iteration.

    Args:
        M: Mean anomaly (radians)
        e: Eccentricity (0 ≤ e < 1)
        tol: Convergence tolerance
        max_iter: Maximum iterations

    Returns:
        Eccentric anomaly E (radians)
    """
    # Initial guess
    E = M + e * math.sin(M) if e < 0.8 else PI

    for _ in range(max_iter):
        dE = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break

    return E


def transfer_orbit_elements(a_departure, a_arrival, mu_sun=1.32712440018e20):
    """Calculate Hohmann transfer orbit elements.

    Args:
        a_departure: Departure orbit SMA (m)
        a_arrival: Arrival orbit SMA (m)
        mu_sun: Central body gravitational parameter (m³/s²)

    Returns:
        dict with 'a' (SMA), 'e' (eccentricity), 'T' (period)
    """
    a_t = (a_departure + a_arrival) / 2.0
    e = abs(a_arrival - a_departure) / (a_arrival + a_departure)
    T = 2.0 * PI * math.sqrt(a_t**3 / mu_sun)
    return {'a': a_t, 'e': e, 'T': T}


# =============================================================================
# KINETIC ENERGY
# =============================================================================

def transfer_time_to_orbit(v_inf_earth, a_departure, a_arrival,
                           mu_sun=1.32712440018e20):
    """Calculate heliocentric transfer time for arbitrary departure v_inf.

    For a Hohmann transfer, v_inf equals the Hohmann excess velocity and the
    transfer orbit is tangent to both orbits.  For v_inf > Hohmann, the
    transfer orbit is a faster ellipse or hyperbola.

    Method:
      1. Compute departure heliocentric speed: v_dep = v_circ + v_inf
      2. Determine transfer orbit energy  →  semi-major axis  →  eccentricity
      3. For elliptic (a > 0): find eccentric anomaly at arrival, use Kepler
      4. For hyperbolic (a < 0): find hyperbolic anomaly, use Kepler analog

    Args:
        v_inf_earth: Excess velocity at Earth departure (m/s)
        a_departure: Departure heliocentric SMA (m), e.g. A_EARTH
        a_arrival: Target heliocentric SMA (m), e.g. A_MARS
        mu_sun: Sun gravitational parameter (m³/s²)

    Returns:
        dict with:
          't_transfer': transfer time (s)
          'a_transfer': transfer SMA (m), negative if hyperbolic
          'e_transfer': eccentricity
          'v_arrival': heliocentric speed at arrival orbit (m/s)
          'v_inf_arrival': excess speed over circular at arrival (m/s)
    """
    r1 = a_departure
    r2 = a_arrival

    # Departure heliocentric speed (prograde injection)
    v_circ_dep = math.sqrt(mu_sun / r1)
    v_dep = v_circ_dep + v_inf_earth

    # Specific orbital energy and angular momentum
    eps = 0.5 * v_dep**2 - mu_sun / r1
    h = r1 * v_dep  # purely tangential at departure

    # Semi-major axis (negative if hyperbolic)
    if abs(eps) < 1e-6:
        # Parabolic (edge case) — treat as very large ellipse
        a_t = 1e20
    else:
        a_t = -mu_sun / (2.0 * eps)

    # Semi-latus rectum and eccentricity
    p = h**2 / mu_sun
    e = math.sqrt(max(0.0, 1.0 - p / a_t)) if a_t > 0 else math.sqrt(1.0 + p / abs(a_t))

    # Check if orbit reaches r2
    if e < 1.0:
        r_apo = a_t * (1.0 + e)
        if r2 > r_apo * 1.001:
            return {
                't_transfer': float('inf'),
                'a_transfer': a_t,
                'e_transfer': e,
                'v_arrival': 0.0,
                'v_inf_arrival': 0.0,
            }

    # Speed at r2 (vis-viva)
    v_at_r2 = math.sqrt(max(0.0, 2.0 * (eps + mu_sun / r2)))
    v_circ_arr = math.sqrt(mu_sun / r2)
    v_inf_arr = abs(v_at_r2 - v_circ_arr)

    # True anomaly at departure (r1 = perihelion for prograde injection)
    # At perihelion, true anomaly = 0
    nu1 = 0.0

    # True anomaly at arrival
    cos_nu2 = (p / r2 - 1.0) / e if e > 1e-10 else 0.0
    cos_nu2 = max(-1.0, min(1.0, cos_nu2))
    nu2 = math.acos(cos_nu2)

    if e < 1.0:
        # Elliptic transfer
        # Eccentric anomaly from true anomaly
        E1 = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(nu1 / 2.0),
                                math.sqrt(1.0 + e) * math.cos(nu1 / 2.0))
        E2 = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(nu2 / 2.0),
                                math.sqrt(1.0 + e) * math.cos(nu2 / 2.0))

        # Mean anomaly
        M1 = E1 - e * math.sin(E1)
        M2 = E2 - e * math.sin(E2)

        # Transfer time
        n = math.sqrt(mu_sun / a_t**3)  # mean motion
        dt = (M2 - M1) / n
        if dt < 0:
            dt += 2.0 * PI / n  # wrap around
    else:
        # Hyperbolic transfer
        a_h = abs(a_t)
        # Hyperbolic anomaly from true anomaly
        # tanh(H/2) = sqrt((e-1)/(e+1)) × tan(nu/2)
        tan_half_nu2 = math.tan(nu2 / 2.0)
        tanh_H2 = math.sqrt((e - 1.0) / (e + 1.0)) * tan_half_nu2
        H2 = 2.0 * math.atanh(max(-0.999999, min(0.999999, tanh_H2)))

        # Hyperbolic mean anomaly: M_h = e × sinh(H) - H
        M_h2 = e * math.sinh(H2) - H2

        # Departure at perihelion: H1 = 0, M_h1 = 0
        n_h = math.sqrt(mu_sun / a_h**3)
        dt = M_h2 / n_h

    return {
        't_transfer': abs(dt),
        'a_transfer': a_t,
        'e_transfer': e,
        'v_arrival': v_at_r2,
        'v_inf_arrival': v_inf_arr,
    }


def calc_kinetic_energy(m, v):
    """KE = ½mv²"""
    return 0.5 * m * v**2
