#!/usr/bin/env python3
"""
Mass Driver Physics Module - Calculation functions for mass driver launch

This module contains physics calculations for:
  - Electromagnetic field calculations (B-field, skin depth, impedance)
  - Three thrust models (eddy current, goodness factor, slip×pressure)
  - Thermal calculations (transient heating of reaction plate)
  - Stage handoff logic (frequency-based)
  - Centrifugal loading on the orbital ring
  - HTS tape voltage and current limits

The mass driver LIM differs from the orbital ring LIM in important ways:
  - Sled moves past FIXED stators (vs moving cable past fixed casing)
  - Three stages with different pole pitches operate SEQUENTIALLY
  - Thermal is TRANSIENT (plate heats up over minutes, not steady-state)
  - Power is UNLIMITED (draws from entire ring grid)
  - Multiple repeating units interact with sled simultaneously

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
import md_config as cfg

MU0 = cfg.MU0
G_ACCEL = cfg.G_ACCEL


# =============================================================================
# MATERIAL PROPERTIES (temperature-dependent)
# =============================================================================

def resistivity_at_T(T, mat=None):
    """Electrical resistivity at temperature T (K).

    Uses linear temperature coefficient: ρ(T) = ρ_293 × [1 + α × (T - 293)]

    Args:
        T: Temperature (K)
        mat: Material dict (default: configured plate material)

    Returns:
        Resistivity (Ω·m)
    """
    if mat is None:
        mat = cfg.get_material()
    return mat["rho_e"] * (1.0 + mat["alpha_rho"] * (T - cfg.T_AMBIENT))


def skin_depth(rho_e, f_slip):
    """Classical skin depth.

    δ = sqrt(ρ / (π × μ₀ × f))

    Args:
        rho_e: Resistivity (Ω·m)
        f_slip: Slip frequency (Hz) — must be > 0

    Returns:
        Skin depth (m)
    """
    if f_slip <= 0:
        return 1e6  # Effectively infinite
    return math.sqrt(rho_e / (math.pi * MU0 * f_slip))


# =============================================================================
# ELECTROMAGNETIC FIELD CALCULATIONS
# =============================================================================

def peak_B_field(stage, I_peak):
    """Peak magnetic field at the reaction plate surface.

    B_peak = μ₀ × N × I_peak / (2 × g_eff)

    where g_eff is the effective magnetic gap (air gap + plate thickness/2
    for field penetration).

    Args:
        stage: LIMStageConfig object
        I_peak: Peak phase current (A)

    Returns:
        B_peak (T)
    """
    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0
    return MU0 * stage.n_turns * I_peak / (2.0 * g_eff)


def slip_velocity(v_sled, slip_ratio):
    """Slip velocity of magnetic wave relative to reaction plate.

    v_slip = v_wave - v_sled, where v_wave = v_sled / (1 - s)
    For small s: v_slip ≈ s × v_sled

    But at v_sled = 0, we need a minimum slip velocity to have any
    thrust at all, so we use the exact relation.

    Args:
        v_sled: Sled velocity (m/s)
        slip_ratio: Slip ratio s (0 to 1)

    Returns:
        v_slip (m/s)
    """
    if slip_ratio >= 1.0:
        return v_sled  # Stall condition
    v_wave = v_sled / (1.0 - slip_ratio) if v_sled > 0 else slip_ratio * 100.0
    return v_wave - v_sled


def supply_frequency(v_wave, tau_p):
    """Electrical supply frequency for the traveling magnetic wave.

    f = v_wave / (2 × τ_p)

    Args:
        v_wave: Wave velocity (m/s)
        tau_p: Pole pitch (m)

    Returns:
        Supply frequency (Hz)
    """
    return v_wave / (2.0 * tau_p)


def slip_frequency(v_slip, tau_p):
    """Slip frequency — the frequency of the field as seen by the plate.

    f_slip = v_slip / (2 × τ_p)
    """
    return v_slip / (2.0 * tau_p)


def induced_voltage(stage, B_peak, v_slip):
    """Induced EMF in the stator coils.

    V_ind ≈ 2π × f_supply × N × B_peak × τ_p × w_coil

    This is the back-EMF that the power supply must overcome.

    Args:
        stage: LIMStageConfig
        B_peak: Peak B field (T)
        v_slip: Slip velocity (m/s)

    Returns:
        Induced voltage (V)
    """
    # The wave velocity = v_sled + v_slip, but for voltage we care about
    # the rate of flux change through the coils.
    # V = N × dΦ/dt = N × B × A × ω where A = τ_p × w_coil, ω = 2π × f
    # Simplified: V ≈ N × B_peak × w_coil × v_wave (for traveling field)
    # But we need to be more careful. The coils see the wave at f_supply.
    # V_coil = 2π × f_supply × N × B_peak × τ_p × w_coil / π
    #        = 2 × f_supply × N × B_peak × τ_p × w_coil
    # This is an upper estimate; the actual depends on coil geometry.
    # For now, use V = N × B × w_coil × v_wave (back-EMF scaling)
    # v_wave is not directly available here, so we use the standard formula:
    # V_induced ≈ 4.44 × f × N × Φ_max where Φ_max = B × τ_p × w_coil / π
    # This gives V = 4.44 × f × N × B × τ_p × w_coil / π
    # For our purposes, we track this to enforce the 100 kV limit.
    # Use a simpler direct scaling that matches the orbital ring code:
    f_supply_local = v_slip / (2.0 * stage.tau_p)  # This is f_slip, not f_supply
    # Actually, the back-EMF depends on the FULL wave speed, not just slip.
    # But we don't have v_sled here. We'll compute this in the simulation
    # where we have access to v_sled. Return 0 as placeholder.
    return 0.0  # Computed in simulation loop


def coil_back_emf(stage, B_peak, v_wave):
    """Total back-EMF in stator coils (without power factor correction).

    V_back ≈ N × B_peak × w_coil × v_wave × (2/π)

    This is the FULL terminal voltage including reactive magnetizing
    component. With PFC capacitors, the supply only needs to provide
    the real power component (see coil_voltage_with_pfc).

    Args:
        stage: LIMStageConfig
        B_peak: Peak B at plate (T)
        v_wave: Wave velocity (m/s) = v_sled + v_slip

    Returns:
        Full back-EMF voltage (V)
    """
    return stage.n_turns * B_peak * stage.w_coil * v_wave * (2.0 / math.pi)


def coil_voltage_with_pfc(stage, B_peak, v_slip):
    """Supply voltage with power factor correction capacitors.

    With PFC, the reactive magnetizing current is supplied by capacitor
    banks. The power supply only provides real power. The voltage scales
    with v_SLIP (tens of m/s), not v_WAVE (km/s).

    V_supply ≈ N × B_peak × w_coil × v_slip × (2/π)

    This is why high-speed LIMs are practical despite enormous wave
    velocities — the power electronics only see the slip frequency.

    Args:
        stage: LIMStageConfig
        B_peak: Peak B at plate (T)
        v_slip: Slip velocity (m/s)

    Returns:
        Supply voltage with PFC (V)
    """
    return stage.n_turns * B_peak * stage.w_coil * v_slip * (2.0 / math.pi)


# =============================================================================
# THRUST MODELS
#
# Three models matching the orbital ring code, adapted for mass driver geometry.
# =============================================================================

def thrust_model_1_eddy(stage, B_peak, v_slip, T_plate, n_active):
    """Model 1: Narrow plate eddy current model.

    The reaction plate is narrow relative to the pole pitch (w_plate < τ_p).
    Eddy currents form elongated loops spanning one pole pitch. The return
    path resistance dominates.

    F = (B² × w_plate × t_eff × v_slip) / (ρ × k_geometry)

    where t_eff = min(t_plate, δ) and k_geometry accounts for the eddy
    current path length (~2 × τ_p for the return path).

    Args:
        stage: LIMStageConfig
        B_peak: Peak B field at plate (T)
        v_slip: Slip velocity (m/s)
        T_plate: Plate temperature (K)
        n_active: Number of simultaneously active LIM segments of this type

    Returns:
        Tuple of (thrust_N, P_eddy_W) — thrust force and eddy current losses
    """
    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)

    f_slip = v_slip / (2.0 * stage.tau_p)
    delta = skin_depth(rho, f_slip)
    t_eff = min(cfg.PLATE_THICKNESS, delta)

    # EMF induced in one eddy current loop spanning one pole pitch:
    # ε = B_peak × w_plate × v_slip (for a strip of width w_plate)
    emf = B_peak * cfg.PLATE_HEIGHT * v_slip

    # Resistance of one eddy current loop:
    # The current path is: across the plate (w_plate), down τ_p,
    # back across the plate, and up τ_p.
    # R_loop = ρ × path_length / (t_eff × cross_section_width)
    # path_length ≈ 2 × τ_p + 2 × w_plate ≈ 2 × τ_p (since τ_p >> w_plate)
    # cross_section_width: the eddy current spreads over ~w_plate
    path_length = 2.0 * stage.tau_p + 2.0 * cfg.PLATE_HEIGHT
    R_loop = rho * path_length / (t_eff * cfg.PLATE_HEIGHT)

    # Inductance of the loop (simplified)
    L_loop = MU0 * stage.tau_p * t_eff / cfg.PLATE_HEIGHT

    # Impedance
    omega_slip = 2.0 * math.pi * f_slip
    Z_loop = math.sqrt(R_loop**2 + (omega_slip * L_loop)**2)

    # Current in the eddy loop
    I_eddy = emf / Z_loop if Z_loop > 0 else 0.0

    # Power dissipated in one loop
    P_one_loop = I_eddy**2 * R_loop

    # Number of eddy current loops per active LIM length
    # Each pole pitch hosts one full loop
    n_loops_per_lim = stage.pitches

    # Total number of active loops (across all simultaneously active segments)
    n_loops_total = n_loops_per_lim * n_active

    # Total eddy current power
    P_eddy = P_one_loop * n_loops_total

    # Thrust from power balance: F = P_eddy / v_slip
    thrust = P_eddy / v_slip if v_slip > 0 else 0.0

    return thrust, P_eddy


def thrust_model_2_goodness(stage, B_peak, v_slip, T_plate, n_active):
    """Model 2: Goodness factor model (Laithwaite).

    Uses dimensionless goodness factor G = μ₀ × σ × v_slip × τ_p² / (π × g_eff)
    Thrust per unit area = B² × s × G / (μ₀ × (1 + (s×G)²))

    This model works best when plate width > τ_p.

    Args/Returns: Same as model 1.
    """
    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)
    sigma = 1.0 / rho

    f_slip = v_slip / (2.0 * stage.tau_p)
    delta = skin_depth(rho, f_slip)
    t_eff = min(cfg.PLATE_THICKNESS, delta)

    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0

    # Goodness factor with effective conducting depth
    # G = μ₀ × σ_eff × v_slip × τ_p / (π × g_eff)
    # σ_eff = σ × t_eff / g_eff (accounts for thin plate)
    sigma_eff = sigma * t_eff / g_eff
    G = MU0 * sigma_eff * v_slip * stage.tau_p / (math.pi * g_eff)

    # Slip ratio for the goodness factor model
    # s = v_slip / v_wave; for the model, we use s × G as the key parameter
    # At low sG: F ≈ (B²/μ₀) × A × sG
    # At high sG: F ≈ (B²/μ₀) × A / (sG)
    sG = G  # Since we're using v_slip directly in G, sG ≡ G here

    # Thrust per unit area of air gap
    # f_a = B² / (2μ₀) × 2sG / (1 + sG²)
    f_area = B_peak**2 / (2 * MU0) * 2 * sG / (1 + sG**2)

    # Active area per LIM segment
    A_active = stage.L_active * cfg.PLATE_HEIGHT

    # Total thrust
    thrust = f_area * A_active * n_active

    # Eddy current power = thrust × v_slip (power deposited in plate)
    P_eddy = thrust * v_slip

    return thrust, P_eddy


def thrust_model_3_slip_pressure(stage, B_peak, v_slip, T_plate, n_active):
    """Model 3: Slip × magnetic pressure (theoretical maximum).

    F = s × F_max where F_max = B²/(2μ₀) × A
    s = v_slip / v_wave

    This is the upper bound — actual thrust is always less due to
    impedance mismatch, end effects, and geometry.

    Args/Returns: Same as model 1.
    """
    # Magnetic pressure
    P_mag = B_peak**2 / (2 * MU0)

    # Active area
    A_active = stage.L_active * cfg.PLATE_HEIGHT * n_active

    # F_max
    F_max = P_mag * A_active

    # For this model, slip ratio s is approximate
    # v_wave = v_sled + v_slip, but we don't have v_sled here
    # Use s ≈ v_slip / (v_sled + v_slip); at low slip, s ≈ v_slip/v_sled
    # As an upper bound model, just use s = v_slip / max(v_slip, 1)
    # Actually, the caller should provide v_sled. For now, use sG approach.
    # We'll compute this properly in the simulation. As placeholder:
    # Assume s is small, so F ≈ F_max × s where s is treated as a fraction
    # This model is for comparison only.

    # Use a normalized slip: assume v_wave ≈ v_sled ≈ 1000 m/s as reference
    # The simulation will replace this with actual values.
    # For standalone: s = v_slip / (v_slip + 1000) as rough estimate
    s_approx = v_slip / (v_slip + 5000.0)  # Will be overridden in sim

    thrust = F_max * s_approx
    P_eddy = thrust * v_slip

    return thrust, P_eddy


def calc_thrust(stage, I_peak, v_slip, v_sled, T_plate, n_active, model=1):
    """Calculate thrust and eddy losses using the selected model.

    This is the main entry point for thrust calculations.

    Args:
        stage: LIMStageConfig
        I_peak: Peak phase current (A)
        v_slip: Slip velocity (m/s)
        v_sled: Sled velocity (m/s)
        T_plate: Plate temperature (K)
        n_active: Number of simultaneously active LIM segments of this type
        model: Thrust model (1, 2, or 3)

    Returns:
        dict with keys: thrust, P_eddy, B_peak, f_supply, f_slip, delta,
                        v_wave, V_back, goodness
    """
    B = peak_B_field(stage, I_peak)
    v_wave = v_sled + v_slip
    f_sup = supply_frequency(v_wave, stage.tau_p)
    f_sl = slip_frequency(v_slip, stage.tau_p)

    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)
    delta = skin_depth(rho, f_sl)

    # Back-EMF (full, without PFC)
    V_back = coil_back_emf(stage, B, v_wave)
    # Supply voltage with PFC (scales with v_slip, not v_wave)
    V_pfc = coil_voltage_with_pfc(stage, B, v_slip)

    # Goodness factor (for reporting, even if not using model 2)
    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0
    t_eff = min(cfg.PLATE_THICKNESS, delta)
    sigma = 1.0 / rho
    sigma_eff = sigma * t_eff / g_eff
    G = MU0 * sigma_eff * v_slip * stage.tau_p / (math.pi * g_eff)

    if model == 1:
        thrust, P_eddy = thrust_model_1_eddy(stage, B, v_slip, T_plate, n_active)
    elif model == 2:
        thrust, P_eddy = thrust_model_2_goodness(stage, B, v_slip, T_plate, n_active)
    elif model == 3:
        # For model 3, fix the slip ratio using actual v_sled
        P_mag = B**2 / (2 * MU0)
        A_active = stage.L_active * cfg.PLATE_HEIGHT * n_active
        F_max = P_mag * A_active
        s = v_slip / v_wave if v_wave > 0 else 1.0
        thrust = F_max * s
        P_eddy = thrust * v_slip
    else:
        thrust, P_eddy = thrust_model_1_eddy(stage, B, v_slip, T_plate, n_active)

    return {
        "thrust": thrust,
        "P_eddy": P_eddy,
        "B_peak": B,
        "f_supply": f_sup,
        "f_slip": f_sl,
        "delta": delta,
        "v_wave": v_wave,
        "V_back": V_back,
        "V_pfc": V_pfc,
        "goodness": G,
        "t_eff": t_eff,
    }


# =============================================================================
# THERMAL CALCULATIONS
# =============================================================================

def plate_temp_rise(P_eddy, dt, T_plate):
    """Temperature rise of the reaction plate from eddy current heating.

    Transient model: dT = P × dt / (m × Cp)
    where m is the total plate mass and Cp is specific heat.

    This assumes uniform heating — valid when the sled passes many
    LIM segments and the heat diffuses quickly through the thin plate.

    Args:
        P_eddy: Total eddy current power deposited in plate (W)
        dt: Time step (s)
        T_plate: Current plate temperature (K)

    Returns:
        New plate temperature (K)
    """
    d = cfg.calc_derived()
    mat = d["mat"]
    m_plate = d["plate_mass"]

    # Energy deposited this timestep
    dE = P_eddy * dt

    # Temperature rise
    dT = dE / (m_plate * mat["Cp"])

    return T_plate + dT


def radiation_cooling(T_plate, dt):
    """Radiative cooling of the plate to space.

    P_rad = ε × σ × A_surface × (T⁴ - T_space⁴)

    The plate radiates from both sides (2 × L × h).

    Args:
        T_plate: Plate temperature (K)
        dt: Time step (s)

    Returns:
        New plate temperature after cooling (K)
    """
    d = cfg.calc_derived()
    mat = d["mat"]

    # Radiating surface area (both sides of the plate)
    A_surface = 2.0 * cfg.SLED_LENGTH * cfg.PLATE_HEIGHT

    # Radiative power
    P_rad = mat["emissivity"] * cfg.STEFAN_BOLTZMANN * A_surface * (
        T_plate**4 - cfg.T_SPACE**4)

    # Energy removed
    dE = P_rad * dt

    # Temperature drop
    dT = dE / (d["plate_mass"] * mat["Cp"])

    return T_plate - dT


# =============================================================================
# STAGE SELECTION AND HANDOFF
# =============================================================================

def select_active_stage(v_sled, v_slip, stages=None):
    """Determine which LIM stage should be active based on supply frequency.

    Each stage has a maximum supply frequency (f_handoff) at which it hands
    off to the next stage. The supply frequency depends on the wave velocity
    and the pole pitch: f = v_wave / (2 × τ_p).

    The sled starts on Stage 1 (LV). When f_supply reaches f_handoff for
    that stage, control transfers to Stage 2 (MV), and so on.

    Args:
        v_sled: Current sled velocity (m/s)
        v_slip: Current slip velocity (m/s)
        stages: List of LIMStageConfig (default: cfg.LIM_STAGES)

    Returns:
        Index into stages list (0, 1, or 2)
    """
    if stages is None:
        stages = cfg.LIM_STAGES

    v_wave = v_sled + v_slip

    # Check stages in order — use the lowest-index stage whose
    # supply frequency hasn't exceeded its handoff limit
    for i, stage in enumerate(stages):
        f = supply_frequency(v_wave, stage.tau_p)
        if stage.f_handoff_hz is None or f <= stage.f_handoff_hz:
            return i

    # If all stages are over their handoff, use the last one
    return len(stages) - 1


def count_active_segments(stage_index, stages=None):
    """Count how many LIM segments of the active stage type interact
    with the sled simultaneously.

    The sled is 5 km long. The repeating unit is ~990 m. Within each
    repeating unit, each stage occupies its L_active length.

    The number of simultaneous interactions = sled_length / L_repeating_unit
    (This is the number of repeating units the sled spans.)

    Each repeating unit has exactly one of each stage.

    Args:
        stage_index: Index of the active stage (0, 1, 2)
        stages: List of LIMStageConfig

    Returns:
        Number of simultaneous active segments (float, can be fractional)
    """
    return cfg.SLED_LENGTH / cfg.L_REPEATING_UNIT_WITH_GAPS


# =============================================================================
# CENTRIFUGAL LOADING
# =============================================================================

def centrifugal_force(v_sled, m_total):
    """Centrifugal force on the orbital ring from the accelerating sled.

    At velocity v on a ring of radius R, centrifugal acceleration = v²/R.
    The net outward force depends on whether v > V_orbit (centrifugal
    exceeds gravity) or v < V_orbit.

    F_centrifugal = m × v² / R
    F_gravity = m × g_local
    F_net = F_centrifugal - F_gravity (positive = outward)

    Args:
        v_sled: Sled velocity (m/s)
        m_total: Total sled mass (kg)

    Returns:
        dict with: F_centrifugal, F_gravity, F_net, net_g
    """
    F_cent = m_total * v_sled**2 / cfg.R_ORBIT
    F_grav = m_total * cfg.G_LOCAL
    F_net = F_cent - F_grav
    net_g = F_net / (m_total * G_ACCEL) if m_total > 0 else 0.0

    return {
        "F_centrifugal": F_cent,
        "F_gravity": F_grav,
        "F_net": F_net,
        "net_g": net_g,
    }


def occupant_g_force(thrust, m_total):
    """G-force experienced by occupants from the LIM thrust.

    This is the tangential acceleration only. Centripetal effects are
    handled separately.

    Args:
        thrust: Total thrust (N)
        m_total: Total mass (kg)

    Returns:
        Acceleration in g
    """
    if m_total <= 0:
        return 0.0
    return thrust / (m_total * G_ACCEL)


# =============================================================================
# CONTROLLER HELPERS
# =============================================================================

def throttle_for_voltage(stage, B_peak, v_sled, v_slip, I_current):
    """Reduce current to keep back-EMF below V_max.

    If V_back > V_max, we need to reduce I (which reduces B, which
    reduces V_back). V_back ∝ B ∝ I, so I_new = I × (V_max / V_back).

    Args:
        stage: LIMStageConfig
        B_peak: Current B field (T)
        v_sled: Sled velocity (m/s)
        v_slip: Slip velocity (m/s)
        I_current: Current operating current (A)

    Returns:
        Throttled current (A)
    """
    v_wave = v_sled + v_slip
    V_back = coil_back_emf(stage, B_peak, v_wave)

    if V_back > cfg.VOLTS_MAX and V_back > 0:
        ratio = cfg.VOLTS_MAX / V_back
        return I_current * ratio * 0.95  # 5% margin
    return I_current


def throttle_for_accel(thrust, m_total, max_g, I_current, stage, T_plate,
                       v_slip, v_sled, n_active, model):
    """Reduce current to keep acceleration below max_g.

    If thrust/m > max_g × g, reduce current. Since thrust ∝ B² ∝ I²,
    I_new = I × sqrt(F_target / F_current).

    Returns:
        Throttled current (A)
    """
    a = thrust / m_total if m_total > 0 else 0
    a_max = max_g * G_ACCEL

    if a > a_max and thrust > 0:
        ratio = math.sqrt(a_max * m_total / thrust)
        return I_current * ratio * 0.98  # 2% margin
    return I_current
