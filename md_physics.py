#!/usr/bin/env python3
"""
Mass Driver Physics Module — Model 1 (eddy current) only

Key corrections from previous version:
  - Stator coils are on the ring at 77 K, plate is on the sled at ~293 K
  - Thermal model: each plate section passes each stator briefly
    (exposure time = L_stator / v_sled), accumulating heat across
    all stator encounters along the full launch track
  - Cryo system cools stators between passes (full orbit to recover)
  - Model 1 only — goodness factor does not work for narrow coils

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
    """B_peak = μ₀ N I / (2 g_eff)"""
    g_eff = stage.gap + cfg.PLATE_THICKNESS / 2.0
    return MU0 * stage.n_turns * I_peak / (2.0 * g_eff)

def supply_frequency(v_wave, tau_p):
    """f = v_wave / (2 τ_p)"""
    return v_wave / (2.0 * tau_p)

def slip_frequency(v_slip, tau_p):
    """f_slip = v_slip / (2 τ_p)"""
    return v_slip / (2.0 * tau_p)

def coil_voltage_pfc(stage, B_peak, v_slip):
    """Supply voltage with PFC — scales with v_slip, not v_wave.
    V ≈ N × B × w × v_slip × (2/π)"""
    return stage.n_turns * B_peak * stage.w_coil * v_slip * (2.0 / math.pi)

# ═════════════════════════════════════════════════════════════════════
# THRUST MODEL 1 — NARROW PLATE EDDY CURRENT
#
# Eddy current loops span one pole pitch. The return path (2×τ_p)
# dominates the loop resistance. Each loop carries current I_eddy
# determined by the induced EMF and loop impedance.
#
# F = Σ(P_eddy_per_loop) / v_slip
# ═════════════════════════════════════════════════════════════════════

def calc_thrust(stage, I_peak, v_slip, v_sled, T_plate, n_active):
    """Calculate thrust and eddy losses using Model 1 (eddy current).

    Args:
        stage: LIMStageConfig
        I_peak: Peak phase current (A)
        v_slip: Slip velocity (m/s)
        v_sled: Sled velocity (m/s)
        T_plate: Plate temperature (K)
        n_active: Number of simultaneously active repeating units

    Returns:
        dict with thrust, P_eddy, B_peak, f_supply, f_slip, delta,
             v_wave, V_pfc, t_eff
    """
    mat = cfg.get_material()
    rho = resistivity_at_T(T_plate, mat)

    B = peak_B_field(stage, I_peak)
    v_wave = v_sled + v_slip
    f_sup = supply_frequency(v_wave, stage.tau_p)
    f_sl = slip_frequency(v_slip, stage.tau_p)
    delta = skin_depth(rho, f_sl)
    t_eff = min(cfg.PLATE_THICKNESS, delta)
    V_pfc = coil_voltage_pfc(stage, B, v_slip)

    # ── Eddy current loop model ──────────────────────────────────
    # EMF per loop: ε = B × w_plate × v_slip
    emf = B * cfg.PLATE_HEIGHT * v_slip

    # Loop resistance: R = ρ × path_length / (t_eff × w_path)
    # Path: across plate (w) + down τ_p + across + up τ_p
    path_length = 2.0 * stage.tau_p + 2.0 * cfg.PLATE_HEIGHT
    R_loop = rho * path_length / (t_eff * cfg.PLATE_HEIGHT)

    # Loop inductance
    L_loop = MU0 * stage.tau_p * t_eff / cfg.PLATE_HEIGHT

    # Impedance
    omega_slip = 2.0 * math.pi * f_sl
    Z_loop = math.sqrt(R_loop**2 + (omega_slip * L_loop)**2)

    I_eddy = emf / Z_loop if Z_loop > 0 else 0.0
    P_one_loop = I_eddy**2 * R_loop

    # Loops per active LIM: one per pole pitch
    n_loops_per_lim = stage.pitches
    n_loops_total = n_loops_per_lim * n_active * cfg.N_LIM_SIDES  # both sides

    P_eddy = P_one_loop * n_loops_total
    thrust = P_eddy / v_slip if v_slip > 0 else 0.0

    return {
        "thrust": thrust, "P_eddy": P_eddy,
        "B_peak": B, "f_supply": f_sup, "f_slip": f_sl,
        "delta": delta, "v_wave": v_wave, "V_pfc": V_pfc,
        "t_eff": t_eff, "R_loop": R_loop, "Z_loop": Z_loop,
        "I_eddy": I_eddy, "emf_loop": emf,
    }

# ═════════════════════════════════════════════════════════════════════
# THERMAL MODEL
#
# Each section of the reaction plate passes over many stator segments
# as the sled travels along the ring. The exposure time per stator =
# L_stator / v_sled. The total heat accumulated = Σ(P_eddy × dt) over
# all stator encounters.
#
# The plate radiates to space from both sides.
# The cryo system cools the stators — NOT the plate directly.
# ═════════════════════════════════════════════════════════════════════

def plate_temp_rise(P_eddy, dt, T_plate):
    """Transient plate heating: dT = P × dt / (m_plate × Cp)"""
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
# CENTRIFUGAL / G-FORCE
# ═════════════════════════════════════════════════════════════════════

def centrifugal_force(v_sled, m_total):
    """Net force on ring from sled: F_cent - F_grav"""
    F_cent = m_total * v_sled**2 / cfg.R_ORBIT
    F_grav = m_total * cfg.G_LOCAL
    F_net = F_cent - F_grav
    net_g = F_net / (m_total * G_ACCEL) if m_total > 0 else 0.0
    return {"F_centrifugal": F_cent, "F_gravity": F_grav,
            "F_net": F_net, "net_g": net_g}

def occupant_g(thrust, m_total):
    if m_total <= 0:
        return 0.0
    return thrust / (m_total * G_ACCEL)
