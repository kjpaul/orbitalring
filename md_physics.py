#!/usr/bin/env python3
"""
Mass Driver Physics Module — Model 1 (eddy current) only

Design:
  - Stator coils are on the ring at 77 K (HTS, cryo-cooled)
  - Single reaction plate on the sled, LIMs on BOTH sides
  - Plate starts at ~293 K, heated by eddy currents, cooled by
    radiation to cryo stator surfaces
  - γ-TiAl is paramagnetic — no hysteresis losses, only eddy currents

Thermal:
  - Plate radiates to stator faces at 77 K (view factor ≈ 1.0)
  - Cryo system absorbs radiated heat via LN2 reserves
  - Each ring section has a full orbit (1437 min) to recover
  - HTS coil AC losses are ~10 kW total — negligible

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""
import math
import md_config as cfg

MU0 = cfg.MU0
G_ACCEL = cfg.G_ACCEL

# ═════════════════════════════════════════════════════════════════════
# MATERIAL PROPERTIES
# ═════════════════════════════════════════════════════════════════════

def resistivity_at_T(T, mat=None):
    """ρ(T) = ρ_293 × [1 + α × (T - 293)]"""
    if mat is None:
        mat = cfg.get_material()
    return mat["rho_e"] * (1.0 + mat["alpha_rho"] * (T - cfg.T_AMBIENT))

def skin_depth(rho_e, f_slip):
    """δ = √(ρ / (π μ₀ f))"""
    if f_slip <= 0:
        return 1e6
    return math.sqrt(rho_e / (math.pi * MU0 * f_slip))

# ═════════════════════════════════════════════════════════════════════
# ELECTROMAGNETIC
# ═════════════════════════════════════════════════════════════════════

def peak_B_field(stage, I_peak):
    """B field at plate midplane from ONE side of the LIM.

    Air-core model using 2D rectangular current sheet:
      B = (2/π) × (μ₀ N I / w_coil) × arctan(w_coil / (2 g_eff))

    Accounts for finite coil height and air return path.
    No saturation limit. Typically B ≈ 0.3–0.7 T.

    Returns B at the plate midplane from ONE side. The double-sided
    geometry is handled by N_LIM_SIDES in calc_thrust.
    """
    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0

    return (2.0 / math.pi) * cfg.MU0 * stage.n_turns * I_peak / (
        stage.w_coil) * math.atan(stage.w_coil / (2.0 * g_eff))

def supply_frequency(v_wave, tau_p):
    """f = v_wave / (2 τ_p)"""
    return v_wave / (2.0 * tau_p)

def slip_frequency(v_slip, tau_p):
    """f_slip = v_slip / (2 τ_p)"""
    return v_slip / (2.0 * tau_p)


# ═════════════════════════════════════════════════════════════════════
# THRUST MODEL 1 — NARROW PLATE EDDY CURRENT
#
# Eddy current loops span one pole pitch in the plate. The induced
# EMF drives current around the loop; the return path along τ_p
# dominates the loop resistance for large τ_p.
#
# γ-TiAl is paramagnetic (magnetic susceptibility ~ 10⁻⁵), so there
# are NO hysteresis losses. All EM losses are ohmic (eddy current).
#
# Each loop: ε = B × h × v_slip
#            R = ρ × (2τ_p + 2h) / (t_eff × h)
#            I_eddy = ε / Z_loop
#            P = I² × R
#            F = ΣP / v_slip
# ═════════════════════════════════════════════════════════════════════

def calc_thrust(stage, I_peak, v_slip, v_sled, T_plate, n_active):
    """Calculate thrust and eddy losses using Model 1 (eddy current).

    The plate is paramagnetic γ-TiAl — no hysteresis, only eddy
    current losses contribute to thrust.

    Args:
        stage: LIMStageConfig
        I_peak: Peak phase current (A)
        v_slip: Slip velocity (m/s)
        v_sled: Sled velocity (m/s)
        T_plate: Plate temperature (K)
        n_active: Number of simultaneously active repeating units

    Returns:
        dict with thrust, P_eddy, B_peak, f_supply, f_slip, delta,
             v_wave, t_eff, and eddy loop diagnostics
    """
    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)

    B = peak_B_field(stage, I_peak)
    v_wave = v_sled + v_slip
    f_sup = supply_frequency(v_wave, stage.tau_p)
    f_sl = slip_frequency(v_slip, stage.tau_p)
    delta = skin_depth(rho, f_sl)
    t_eff = min(cfg.PLATE_THICKNESS, delta)

    # ── Stator coil inductance and voltage ───────────────────────
    # L = μ₀ N² A_field / l_coil
    coil_length = stage.tau_p / 3.0  # each phase coil spans τ_p/3
    g_eff_field = stage.gap + cfg.PLATE_THICKNESS / 2.0
    A_field = stage.w_coil * 2 * g_eff_field  # both gaps
    L_coil = MU0 * stage.n_turns**2 * A_field / coil_length

    # V_coil = inductive voltage across coil insulation (must be < 100 kV)
    V_coil = 2.0 * math.pi * f_sup * L_coil * I_peak

    # HTS hysteresis losses (Gömöry perpendicular field model)
    P_hts_total = calc_hysteresis_power_stage(stage, I_peak, f_sup, n_active)

    # ── Eddy current loop model (in reaction plate) ──────────────
    # EMF per loop: ε = B × h_plate × v_slip
    emf = B * cfg.PLATE_HEIGHT * v_slip

    # Loop resistance: R = ρ × path_length / (cross_section_area)
    # Path: across plate height + along τ_p + across + back along τ_p
    # Cross section of current path: t_eff × h_plate
    path_length = 2.0 * stage.tau_p + 2.0 * cfg.PLATE_HEIGHT
    R_loop = rho * path_length / (t_eff * cfg.PLATE_HEIGHT)

    # Loop inductance (approximate: rectangular loop)
    L_loop = MU0 * stage.tau_p * t_eff / cfg.PLATE_HEIGHT

    # Impedance at slip frequency
    omega_slip = 2.0 * math.pi * f_sl
    Z_loop = math.sqrt(R_loop**2 + (omega_slip * L_loop)**2)

    I_eddy = emf / Z_loop if Z_loop > 0 else 0.0
    P_one_loop = I_eddy**2 * R_loop

    # Total: loops per stage × active units × both sides
    n_loops_per_lim = stage.pitches
    n_loops_total = n_loops_per_lim * n_active * cfg.N_LIM_SIDES

    P_eddy = P_one_loop * n_loops_total
    thrust = P_eddy / v_slip if v_slip > 0 else 0.0

    return {
        "thrust": thrust, "P_eddy": P_eddy,
        "B_peak": B, "f_supply": f_sup, "f_slip": f_sl,
        "delta": delta, "v_wave": v_wave,
        "V_coil": V_coil, "L_coil": L_coil, "P_hts_coil": P_hts_total,
        "t_eff": t_eff, "R_loop": R_loop, "Z_loop": Z_loop,
        "I_eddy": I_eddy, "emf_loop": emf,
    }

# ═════════════════════════════════════════════════════════════════════
# HTS HYSTERESIS LOSS — GÖMÖRY PERPENDICULAR FIELD MODEL
#
# HTS tape at 77 K experiences AC magnetization losses as the stator
# field cycles. The Gömöry model calculates the energy loss per meter
# of tape per cycle from the coil self-field and critical current:
#
#   Q_hyst = B_coil_peak × Ic_per_layer × w_tape × sin(α_tape_field)
#
# where B_coil_peak = μ₀ × N × I / l_coil is the field inside the
# winding, and α_tape_field = 20° is the perpendicular field angle.
# ═════════════════════════════════════════════════════════════════════

def calc_hysteresis_power_stage(stage, I_operating, f_supply, n_active):
    """Total HTS hysteresis power using Gömöry perpendicular field model.

    Args:
        stage: LIMStageConfig for the active stage
        I_operating: Operating current (A)
        f_supply: Electrical supply frequency (Hz)
        n_active: Number of repeating units the sled spans

    Returns:
        Total hysteresis power (W) for all active coils in this stage
    """
    # Coil self-field (field inside the winding)
    l_coil = stage.tau_p / 3.0                    # coil slot length
    B_coil_peak = MU0 * stage.n_turns * I_operating / l_coil

    # Tape geometry per layer
    coil_depth = stage.n_turns * 0.15e-3          # winding depth (m)
    turn_perimeter = 2.0 * (stage.w_coil + coil_depth)
    L_tape_per_layer = stage.n_turns * turn_perimeter

    # Gömöry loss per meter of tape per cycle (J/m/cycle)
    w_tape = stage.hts_width_mm * 1e-3
    Q_hyst = B_coil_peak * stage.I_c_per_layer * w_tape * cfg.SIN_ALPHA

    # Power per phase coil = layers × tape_length × Q_hyst × frequency
    P_hts_per_coil = stage.hts_layers * L_tape_per_layer * Q_hyst * f_supply

    # Total: 3 phases × 2 sides × n_active repeating units
    n_coils_active = 3 * cfg.N_LIM_SIDES * n_active
    return P_hts_per_coil * n_coils_active


# ═════════════════════════════════════════════════════════════════════
# CRYOGENIC RADIATOR WIDTH
#
# LN2 absorbs HTS hysteresis heat at 77 K; radiators at ~100 K reject
# it to space (~3 K background). Two-sided radiator panels run along
# the stator section.
# ═════════════════════════════════════════════════════════════════════

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

# ═════════════════════════════════════════════════════════════════════
# THERMAL MODEL
#
# Plate heating: eddy current ohmic losses (P = I²R in plate loops)
# Plate cooling: radiation to cryo stator faces at 77 K
#
# γ-TiAl is paramagnetic → NO hysteresis losses in the plate.
# All heating is from induced eddy currents.
# ═════════════════════════════════════════════════════════════════════

def plate_temp_rise(P_eddy, dt, T_plate):
    """Transient plate heating from eddy currents: dT = P×dt / (m×Cp)"""
    d = cfg.calc_derived()
    mat = d["mat"]
    dE = P_eddy * dt
    dT = dE / (d["plate_mass"] * mat["Cp"])
    return T_plate + dT

def plate_sink_temperature():
    """Temperature of the surface the plate radiates to.

    With water cooling: thermal shield at T_SHIELD_WATER (~350 K).
    With LN2: plate faces cryo stator surfaces at T_STATOR (77 K).
    """
    if cfg.PLATE_COOLANT == "water":
        return cfg.T_SHIELD_WATER
    return cfg.T_STATOR

def calc_stator_heat_flux(T_plate):
    """Radiation heat flux from hot plate to thermal shield / stator surfaces.

    Both sides of the plate radiate to the nearest surface (thermal shield
    or cryo stator face). Returns total power and per-metre heat flux.

    Returns:
        (P_rad, q_per_m): total radiated power (W) and heat flux
        per metre of stator currently under the sled (W/m)
    """
    mat = cfg.get_material()
    T_sink = plate_sink_temperature()
    A_surface = 2.0 * cfg.SLED_LENGTH * cfg.PLATE_HEIGHT
    P_rad = mat["emissivity"] * cfg.STEFAN_BOLTZMANN * A_surface * (
        T_plate**4 - T_sink**4)
    q_per_m = P_rad / cfg.SLED_LENGTH
    return P_rad, q_per_m

def radiation_cooling(T_plate, dt):
    """Radiative cooling from plate to thermal shield / stator surfaces.

    View factor ≈ 1.0 (parallel plates, 100mm gap, extending km).
    """
    d = cfg.calc_derived()
    P_rad, _ = calc_stator_heat_flux(T_plate)
    dE = P_rad * dt
    dT = dE / (d["plate_mass"] * d["mat"]["Cp"])
    return T_plate - dT

def calc_coolant_absorption():
    """Absorption capacity of the plate radiation coolant.

    Returns:
        (capacity, name, density): J/kg, display name, kg/m³
    """
    if cfg.PLATE_COOLANT == "water":
        sensible = cfg.C_P_WATER * (cfg.T_WATER_BOIL - cfg.T_WATER_SUPPLY)
        capacity = sensible + cfg.L_V_WATER
        return capacity, "Water", cfg.RHO_WATER
    sensible = cfg.C_P_LN2 * (cfg.T_LN2_BOIL - cfg.T_LN2_SUPPLY)
    capacity = sensible + cfg.L_V_LN2
    return capacity, "LN2", cfg.RHO_LN2

def calc_energy_per_m_per_pass(q_per_m, v_sled):
    """Energy deposited per metre of stator per sled pass.

    A fixed stator section sees the sled for t_exposure = SLED_LENGTH / v_sled.
    During that time it absorbs q_per_m watts per metre.

    Returns:
        Energy per metre per pass (J/m)
    """
    if v_sled <= 0:
        return 0.0
    return q_per_m * cfg.SLED_LENGTH / v_sled

# ═════════════════════════════════════════════════════════════════════
# STAGE SELECTION
# ═════════════════════════════════════════════════════════════════════

def select_active_stage(v_sled, v_slip):
    """Select stage by supply frequency: use lowest-index stage whose
    f_supply hasn't exceeded its handoff limit."""
    v_wave = v_sled + v_slip
    for i, stage in enumerate(cfg.LIM_STAGES):
        f = supply_frequency(v_wave, stage.tau_p)
        if stage.f_handoff_hz is None or f <= stage.f_handoff_hz:
            return i
    return len(cfg.LIM_STAGES) - 1

def count_active_segments(stage=None):
    """Number of repeating units the sled spans for the active stage.

    Each stage's coils tile the ring independently. When S1 is active,
    only S1 coils are energized — the sled passes over S2/S3/S4 coils
    without activating them. So n_active depends on the active stage's
    own L_active, not the combined repeating unit.

    If stage is None, uses the combined repeating unit (legacy behavior).
    """
    if stage is not None:
        return cfg.SLED_LENGTH / (stage.L_active + cfg.L_STAGE_GAP)
    return cfg.SLED_LENGTH / cfg.L_REPEATING_UNIT_WITH_GAPS

# ═════════════════════════════════════════════════════════════════════
# OCCUPANT G-FORCE
#
# The occupant experiences two acceleration components:
# 1. Thrust (tangential): a_thrust = F / m (along track direction)
# 2. Net radial: a_radial = v²/R - g_local
#    - At rest: -g_local ≈ -9.073 m/s² (toward Earth)
#    - At v_orbit ≈ 7755 m/s: zero (free fall)
#    - At 15 km/s: +24.9 m/s² (away from Earth) = +2.54g
#
# Total g-load = √(a_thrust² + a_radial²)
# ═════════════════════════════════════════════════════════════════════

def occupant_forces(thrust, v_sled, m_total):
    """Calculate all g-force components felt by occupants.

    Returns dict with:
        a_thrust: tangential acceleration (m/s²)
        a_radial: net radial acceleration (m/s²) = centrifugal - gravity
        g_thrust: tangential in g
        g_radial: radial in g (negative = toward Earth)
        g_total:  magnitude of total g-vector
    """
    a_thrust = thrust / m_total if m_total > 0 else 0.0
    a_centrifugal = v_sled**2 / cfg.R_ORBIT
    a_radial = a_centrifugal - cfg.G_LOCAL  # positive = away from Earth

    g_thrust = a_thrust / G_ACCEL
    g_radial = a_radial / G_ACCEL
    g_total = math.sqrt(g_thrust**2 + g_radial**2)

    return {
        "a_thrust": a_thrust, "a_radial": a_radial,
        "g_thrust": g_thrust, "g_radial": g_radial,
        "g_total": g_total,
        "a_centrifugal": a_centrifugal,
    }

def centrifugal_force(v_sled, m_total):
    """Net radial force on ring from sled (for ring structural loading)."""
    F_cent = m_total * v_sled**2 / cfg.R_ORBIT
    F_grav = m_total * cfg.G_LOCAL
    F_net = F_cent - F_grav
    net_g = F_net / (m_total * G_ACCEL) if m_total > 0 else 0.0
    return {"F_centrifugal": F_cent, "F_gravity": F_grav,
            "F_net": F_net, "net_g": net_g}
