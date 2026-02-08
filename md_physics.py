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
    """B_peak = μ₀ N I / (2 g_eff)
    g_eff = air gap + half plate thickness (field measured at plate center)
    """
    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0
    return MU0 * stage.n_turns * I_peak / (2.0 * g_eff)

def supply_frequency(v_wave, tau_p):
    """f = v_wave / (2 τ_p)"""
    return v_wave / (2.0 * tau_p)

def slip_frequency(v_slip, tau_p):
    """f_slip = v_slip / (2 τ_p)"""
    return v_slip / (2.0 * tau_p)

def coil_supply_voltage(stage, B_peak, v_slip):
    """Supply voltage with power factor correction.

    With PFC capacitor banks on the stator, the reactive magnetizing
    current is handled locally. The power supply only provides real
    power, so the voltage scales with v_slip (tens of m/s) rather
    than v_wave (km/s).

    V_supply ≈ N × B × w_coil × v_slip × (2/π)
    """
    return stage.n_turns * B_peak * stage.w_coil * v_slip * (2.0 / math.pi)

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
             v_wave, V_supply, t_eff, and eddy loop diagnostics
    """
    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)

    B = peak_B_field(stage, I_peak)
    v_wave = v_sled + v_slip
    f_sup = supply_frequency(v_wave, stage.tau_p)
    f_sl = slip_frequency(v_slip, stage.tau_p)
    delta = skin_depth(rho, f_sl)
    t_eff = min(cfg.PLATE_THICKNESS, delta)
    V_supply = coil_supply_voltage(stage, B, v_slip)

    # ── Stator coil inductance and voltage ───────────────────────
    # L = μ₀ N² A_field / l_coil
    coil_length = stage.tau_p / 3.0  # each phase coil spans τ_p/3
    g_eff_field = stage.gap + cfg.PLATE_THICKNESS / 2.0
    A_field = stage.w_coil * 2 * g_eff_field  # both gaps
    L_coil = MU0 * stage.n_turns**2 * A_field / coil_length

    # V_coil = inductive voltage across coil insulation (must be < 100 kV)
    V_coil = 2.0 * math.pi * f_sup * L_coil * I_peak

    # HTS AC losses per coil (hysteresis + coupling + stabilizer eddy)
    # ~0.5 J/m per cycle at I/Ic ≈ 0.8
    tape_per_coil = stage.n_turns * 2.0 * (stage.w_coil + 0.1)
    P_hts_per_coil = 0.5 * tape_per_coil * f_sup
    n_coils = 3 * cfg.N_LIM_SIDES * n_active
    P_hts_total = P_hts_per_coil * n_coils

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
        "delta": delta, "v_wave": v_wave, "V_supply": V_supply,
        "V_coil": V_coil, "L_coil": L_coil, "P_hts_coil": P_hts_total,
        "t_eff": t_eff, "R_loop": R_loop, "Z_loop": Z_loop,
        "I_eddy": I_eddy, "emf_loop": emf,
    }

# ═════════════════════════════════════════════════════════════════════
# HTS COIL THERMAL
#
# HTS tape at 77 K has AC losses (coupling + magnetization) but they
# are very small: ~10 kW total for all active coils at peak frequency.
# The cryo system easily handles this with LN2 reserves.
#
# No hysteresis in the stator core — the core is non-magnetic (the
# field is created by HTS coils alone, no iron).
# ═════════════════════════════════════════════════════════════════════

def hts_ac_loss_estimate(stage, I_peak, f_supply, n_active):
    """Estimate HTS tape AC losses in active stator coils.

    Model: P_ac ≈ p₀ × (I/I_ref) × (f/f_ref) per meter of tape
    p₀ ≈ 0.5 mW/m at I_ref=100 A, f_ref=50 Hz (conservative YBCO)

    Returns:
        Total AC loss power (W) for all active coils
    """
    p_0 = 0.5e-3        # W/m at reference conditions
    I_ref = 100.0        # A
    f_ref = 50.0         # Hz

    P_per_m = p_0 * (I_peak / I_ref) * (max(f_supply, 0.1) / f_ref)

    # Tape length: turns × perimeter per turn × phases × sides × units
    perimeter_per_turn = 2.0 * (stage.w_coil + 0.1)  # approx coil depth 100mm
    tape_per_coil = stage.n_turns * perimeter_per_turn
    n_coils = 3 * cfg.N_LIM_SIDES * n_active  # 3 phases × 2 sides × units
    total_tape = n_coils * tape_per_coil

    return P_per_m * total_tape

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

def radiation_cooling(T_plate, dt):
    """Radiative cooling from plate to cryo-cooled stator surfaces.
    
    Both sides of the plate face stator surfaces at T_STATOR (77 K).
    View factor ≈ 1.0 (parallel plates, 100mm gap, extending km).
    The cryo system absorbs this heat — it has the full ring grid's
    power and LN2 reserves, and only 5 km of 41,646 km ring is active.
    """
    d = cfg.calc_derived()
    mat = d["mat"]
    A_surface = 2.0 * cfg.SLED_LENGTH * cfg.PLATE_HEIGHT
    P_rad = mat["emissivity"] * cfg.STEFAN_BOLTZMANN * A_surface * (
        T_plate**4 - cfg.T_STATOR**4)
    dE = P_rad * dt
    dT = dE / (d["plate_mass"] * mat["Cp"])
    return T_plate - dT

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

def count_active_segments():
    """Number of repeating units the sled spans."""
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
