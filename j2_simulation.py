#!/usr/bin/env python3
"""
J2 Perturbation Simulation for an Inclined Orbital Ring

Models the effect of Earth's J2 oblateness on an orbital ring built at a
non-equatorial inclination.  The ring orbits at 250 km altitude and consists
of a fast-spinning internal cable held in place by a ground-synchronous
casing through magnetic levitation.

Four phases are modelled:
  Phase 1 — Free-flying cable loop (pre-deployment)
  Phase 2 — Deployment (LIM spin-up of cable, casing deceleration)
  Phase 3 — Stabilization before anchor drop (LIMs resist precession)
  Phase 4 — Anchor lines deployed (forces transmitted to ground)

Usage:
    python j2_simulation.py [options]

Options:
    --inclination=N   Ring inclination in degrees (default 5.0)
    --cable=DIR       Cable direction: prograde or retrograde (default prograde)
    --anchors=N       Number of anchor stations (default 800)
    --plots           Generate all plots
    --all             Show all output and generate all plots

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
import sys
import os

# Ensure UTF-8 output on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

PI = math.pi

# =============================================================================
# CONSTANTS
# =============================================================================

# --- Earth ---
G_CONST = 6.674e-11            # gravitational constant (m^3 kg^-1 s^-2)
M_EARTH = 5.972e24             # kg
GM = 3.986004418e14            # m^3/s^2
R_E = 6_378_137                # equatorial radius (m)
J2 = 1.08263e-3                # Earth's J2 coefficient
OMEGA_EARTH = 7.2921159e-5     # Earth's rotation rate (rad/s)
SIDEREAL_DAY = 86164.1         # seconds
G_SURFACE = 9.807              # m/s^2

# --- Orbital parameters at 250 km altitude ---
H_ORBIT = 250_000              # altitude (m)
R_ORBIT = R_E + H_ORBIT        # orbital radius (m)
V_ORBIT = math.sqrt(GM / R_ORBIT)              # circular velocity (m/s)
N_ORBIT = math.sqrt(GM / R_ORBIT**3)           # mean motion (rad/s)
T_ORBIT = 2 * PI / N_ORBIT                     # orbital period (s)
C_RING = 2 * PI * R_ORBIT                      # ring circumference (m)
G_ALT = GM / R_ORBIT**2                        # gravity at altitude (m/s^2)

# --- Ring mass properties ---
M_CABLE_PER_M = 100_811        # cable + attached hardware (kg/m)
M_CASING_PER_M = 12_000        # casing + payload track (kg/m)
M_TOTAL_PER_M = M_CABLE_PER_M + M_CASING_PER_M

M_CABLE_TOTAL = M_CABLE_PER_M * C_RING
M_CASING_TOTAL = M_CASING_PER_M * C_RING
M_RING_TOTAL = M_TOTAL_PER_M * C_RING

# --- Post-deployment velocities ---
V_CABLE_PROGRADE = 8620.0      # m/s
V_CABLE_RETROGRADE = 8735.5    # m/s (magnitude; direction is retrograde)
V_CASING = 483.3               # m/s (ground-synchronous)

# --- LIM parameters ---
LIM_SPACING = 500.0            # m
MAX_SITE_POWER = 8e6           # W (8 MW per LIM site)
N_LIM_SITES = round(C_RING / LIM_SPACING)
P_RING_TOTAL = 666e9           # total generation capacity (W)

# --- Anchor line parameters ---
N_ANCHORS_DEFAULT = 800
ANCHOR_LENGTH = 250_000        # m (approximate)
SIGMA_CNT = 14e9               # Pa (CNT tensile strength)
SAFETY_FACTOR = 2.0
SIGMA_OPERATING = SIGMA_CNT / SAFETY_FACTOR  # 7 GPa
RHO_CNT = 1700                 # kg/m^3

# --- Default ---
I_DEFAULT = 5.0                # degrees


# =============================================================================
# PHYSICS HELPER FUNCTIONS
# =============================================================================

def j2_precession_rate(i_rad, r=R_ORBIT):
    """Nodal precession rate for a point mass in a circular orbit (rad/s).

    Omega_dot = -(3/2) n J2 (R_E/a)^2 cos(i)
    Negative = westward regression for prograde orbits (i < 90 deg).
    """
    n = math.sqrt(GM / r**3)
    return -1.5 * n * J2 * (R_E / r)**2 * math.cos(i_rad)


def j2_torque_on_ring(i_rad, M_ring, r=R_ORBIT):
    """J2 gravitational torque on a circular ring.

    tau = 3 GM J2 R_E^2 M sin(i) cos(i) / (2 r^3)

    This is the torque about the line of nodes, causing nodal precession.
    """
    return 3 * GM * J2 * R_E**2 * M_ring * math.sin(i_rad) * math.cos(i_rad) / (2 * r**3)


def j2_out_of_plane_accel(i_rad, u, r=R_ORBIT):
    """Out-of-plane J2 acceleration on a ring element at argument of latitude u.

    a_W = -(3 GM J2 R_E^2) / (2 r^4) * sin(2i) * sin(u)
    """
    return -(3 * GM * J2 * R_E**2) / (2 * r**4) * math.sin(2 * i_rad) * math.sin(u)


def ring_angular_momentum(m_cable_per_m, v_cable, m_casing_per_m, v_casing, r=R_ORBIT):
    """Total angular momentum of the ring system.

    L = (m_cable * v_cable + m_casing * v_casing) * r * C
    """
    C = 2 * PI * r
    return (m_cable_per_m * v_cable + m_casing_per_m * v_casing) * r * C


# =============================================================================
# PHASE 1: FREE-FLYING CABLE (PRE-DEPLOYMENT)
# =============================================================================

def run_phase1(i_deg):
    """Phase 1: Free-flying cable loop at orbital velocity."""
    i_rad = math.radians(i_deg)

    # Ring properties (cable only, no casing yet)
    M_ring = M_CABLE_TOTAL
    L_ring = M_ring * V_ORBIT * R_ORBIT

    # Point-mass precession rate
    omega_prec_pt = j2_precession_rate(i_rad)
    omega_prec_pt_deg_day = omega_prec_pt * 180 / PI * 86400

    # J2 torque on the ring
    tau_j2 = j2_torque_on_ring(i_rad, M_ring)

    # Precession from torque: Omega_dot = -tau / (L sin i)
    omega_prec_ring = -tau_j2 / (L_ring * math.sin(i_rad))
    omega_prec_ring_deg_day = omega_prec_ring * 180 / PI * 86400

    # Peak lateral force per meter
    f_max_per_m = j2_out_of_plane_accel(i_rad, PI / 2) * M_CABLE_PER_M
    # (negative sign means toward equator; report magnitude)
    f_max_per_m = abs(f_max_per_m)

    # Total lateral force (absolute, integrated)
    F_total_lateral = 4 * tau_j2 / (PI * R_ORBIT)

    print("\n" + "=" * 78)
    print("PHASE 1: FREE-FLYING CABLE (PRE-DEPLOYMENT)")
    print("=" * 78)

    print(f"\n  Inclination:          {i_deg:.1f}°")
    print(f"  Altitude:             {H_ORBIT / 1000:.0f} km")
    print(f"  Orbital radius:       {R_ORBIT / 1e6:.4f} Mm ({R_ORBIT:,.0f} m)")
    print(f"  Orbital velocity:     {V_ORBIT:,.1f} m/s")
    print(f"  Orbital period:       {T_ORBIT / 60:.1f} min")
    print(f"  Ring circumference:   {C_RING / 1000:,.0f} km")
    print(f"  Cable mass/m:         {M_CABLE_PER_M:,} kg/m")
    print(f"  Total cable mass:     {M_CABLE_TOTAL:.3e} kg ({M_CABLE_TOTAL / 1e12:.2f} Tt)")
    print(f"  Angular momentum:     {L_ring:.3e} kg·m²/s")

    # --- Derivation ---
    print(f"\n  --- DERIVATION: Ring vs Point-Mass Precession ---")
    print(f"\n  For a ring element at argument of latitude u on an inclined ring:")
    print(f"    sin(φ) = sin(i) × sin(u)")
    print(f"\n  The out-of-plane J2 acceleration on this element is:")
    print(f"    a_W(u) = -(3 GM J2 R_E²) / (2 r⁴) × sin(2i) × sin(u)")
    print(f"\n  [Derived by projecting the J2 perturbation acceleration vector")
    print(f"   in ECI coordinates onto the orbit-normal (W) direction.]")
    print(f"\n  The torque about the line of nodes is:")
    print(f"    dτ = a_W(u) × dm × r × sin(u)")
    print(f"         where dm = λ × r × du  (λ = linear mass density)")
    print(f"\n  Integrating over the full ring (u = 0 to 2π):")
    print(f"    τ = -(3μJ₂R_E²λ)/(2r²) × sin(2i) × ∫ sin²(u) du")
    print(f"      = -(3μJ₂R_E²λ)/(2r²) × sin(2i) × π")
    print(f"\n  Using λ = M/(2πr) and sin(2i) = 2 sin(i) cos(i):")
    print(f"    τ = -(3 GM J₂ R_E² M sin(i) cos(i)) / (2 r³)")
    print(f"\n  The precession rate is:")
    print(f"    Ω̇ = -τ / (L sin i) = -(3/2) n J₂ (R_E/r)² cos(i)")
    print(f"\n  This is IDENTICAL to the point-mass formula.")
    print(f"\n  Physical reason: a continuous ring simultaneously occupies all")
    print(f"  arguments of latitude. Integrating J2 over the ring is equivalent")
    print(f"  to orbit-averaging for a point mass, because ∫ sin²(u) du / (2π)")
    print(f"  gives the same factor for both calculations.")

    # --- Numerical results ---
    print(f"\n  --- NUMERICAL RESULTS ---")
    print(f"\n  Point-mass precession rate:")
    print(f"    Ω̇ = -(3/2) × {N_ORBIT:.4e} × {J2} × ({R_E/R_ORBIT:.6f})² × cos({i_deg}°)")
    print(f"    Ω̇ = {omega_prec_pt:.6e} rad/s")
    print(f"    Ω̇ = {omega_prec_pt_deg_day:.4f} °/day")

    print(f"\n  Ring-integrated precession rate:")
    print(f"    τ_J2 = {tau_j2:.4e} N·m")
    print(f"    L   = {L_ring:.4e} kg·m²/s")
    print(f"    Ω̇  = -τ/(L sin i) = {omega_prec_ring:.6e} rad/s")
    print(f"    Ω̇  = {omega_prec_ring_deg_day:.4f} °/day")

    print(f"\n  Comparison:")
    diff_pct = abs(omega_prec_ring - omega_prec_pt) / abs(omega_prec_pt) * 100
    print(f"    Point-mass:  {omega_prec_pt_deg_day:+.4f} °/day")
    print(f"    Ring:        {omega_prec_ring_deg_day:+.4f} °/day")
    print(f"    Difference:  {diff_pct:.6f}% (identical to numerical precision)")

    print(f"\n  J2 out-of-plane force on ring:")
    print(f"    Peak lateral acceleration:  {abs(j2_out_of_plane_accel(i_rad, PI/2)):.4e} m/s²")
    print(f"    Peak lateral force/m:       {f_max_per_m:.4f} N/m")
    print(f"    Total lateral force:        {F_total_lateral/1e9:.3f} GN (absolute, both halves)")

    # --- Ground track drift table ---
    print(f"\n  --- GROUND TRACK DRIFT (ascending node longitude) ---")
    print(f"\n  In the Earth-fixed frame, the node drifts at:")
    print(f"    Ω̇_ground = Ω̇_J2 - ω_Earth = {omega_prec_pt:.4e} - {OMEGA_EARTH:.4e}")
    omega_ground = omega_prec_pt - OMEGA_EARTH
    omega_ground_deg_day = omega_ground * 180 / PI * 86400
    print(f"              = {omega_ground:.4e} rad/s = {omega_ground_deg_day:.2f} °/day")
    print(f"\n  The ground track regresses ~{abs(omega_ground_deg_day):.0f}°/day, or one full revolution")
    print(f"  every {360 / abs(omega_ground_deg_day):.2f} days. This is dominated by Earth's rotation,")
    print(f"  not J2 (Earth rotation = {OMEGA_EARTH * 180/PI * 86400:.1f} °/day, J2 = {abs(omega_prec_pt_deg_day):.1f} °/day).")

    print(f"\n  {'Day':>6}  {'Ω_inertial':>14}  {'Ω_ground':>14}  {'Drift':>10}")
    print(f"  {'':>6}  {'(°)':>14}  {'(°)':>14}  {'(°)':>10}")
    print(f"  " + "-" * 50)
    for day in range(0, 31, 3):
        omega_inertial_deg = omega_prec_pt_deg_day * day
        omega_ground_deg = omega_ground_deg_day * day
        drift = omega_ground_deg % 360
        if drift > 180:
            drift -= 360
        print(f"  {day:6d}  {omega_inertial_deg:14.2f}  {omega_ground_deg:14.1f}  {drift:10.1f}")

    print(f"\n  Angular momentum vector:")
    print(f"    |L| = {L_ring:.4e} kg·m²/s")
    print(f"    Direction: {i_deg:.1f}° from polar axis")
    print(f"    L precesses around the polar axis at {abs(omega_prec_pt_deg_day):.2f} °/day")
    print(f"    Full precession cycle: {360 / abs(omega_prec_pt_deg_day):.1f} days")

    return {
        'omega_prec': omega_prec_pt,
        'omega_prec_deg_day': omega_prec_pt_deg_day,
        'tau_j2': tau_j2,
        'L': L_ring,
        'f_max_per_m': f_max_per_m,
    }


# =============================================================================
# PHASE 2: DEPLOYMENT (LIM SPIN-UP)
# =============================================================================

def run_phase2(i_deg, cable_dir="prograde"):
    """Phase 2: Cable accelerated to super-orbital, casing to ground-sync."""
    i_rad = math.radians(i_deg)

    v_cable_final = V_CABLE_PROGRADE if cable_dir == "prograde" else -V_CABLE_RETROGRADE
    v_casing_final = V_CASING

    # J2 torque is constant (depends on total mass, not velocity)
    tau_j2 = j2_torque_on_ring(i_rad, M_RING_TOTAL)

    print("\n" + "=" * 78)
    print(f"PHASE 2: DEPLOYMENT ({cable_dir.upper()} CABLE)")
    print("=" * 78)

    print(f"\n  Cable: {V_ORBIT:.1f} → {v_cable_final:+.1f} m/s ({cable_dir})")
    print(f"  Casing: {V_ORBIT:.1f} → {v_casing_final:.1f} m/s (ground-synchronous)")
    print(f"  Total ring mass: {M_RING_TOTAL:.3e} kg")

    # Angular momentum at start
    L_start = ring_angular_momentum(M_CABLE_PER_M, V_ORBIT, M_CASING_PER_M, V_ORBIT)

    # Angular momentum at end
    L_end = ring_angular_momentum(M_CABLE_PER_M, v_cable_final, M_CASING_PER_M, v_casing_final)

    print(f"\n  Angular momentum (L = Σ m_i v_i r × C):")
    print(f"    Start (all at v_orbit):  L = {L_start:.4e} kg·m²/s")
    print(f"    End (deployed):          L = {L_end:.4e} kg·m²/s")

    # Check momentum conservation
    l_start_per_m = M_TOTAL_PER_M * V_ORBIT
    l_end_per_m = M_CABLE_PER_M * v_cable_final + M_CASING_PER_M * v_casing_final
    delta_l = l_end_per_m - l_start_per_m
    print(f"\n  Momentum per metre:")
    print(f"    Start: {l_start_per_m:,.0f} kg·m/s per m")
    print(f"    End:   {l_end_per_m:,.0f} kg·m/s per m")
    print(f"    Δ:     {delta_l:+,.0f} kg·m/s per m ({delta_l/l_start_per_m*100:+.4f}%)")

    if cable_dir == "prograde":
        print(f"\n  LIMs are internal forces → angular momentum conserved (Δ ≈ 0).")
        print(f"  Cable gains +{M_CABLE_PER_M * (V_CABLE_PROGRADE - V_ORBIT)/1e6:.1f} MN·s/m,")
        print(f"  casing loses {M_CASING_PER_M * (V_ORBIT - V_CASING)/1e6:.1f} MN·s/m. These balance.")
    else:
        delta_cable = M_CABLE_PER_M * (-V_CABLE_RETROGRADE - V_ORBIT)
        delta_casing = M_CASING_PER_M * (V_CASING - V_ORBIT)
        print(f"\n  RETROGRADE CABLE: angular momentum is NOT conserved.")
        print(f"  Cable Δp/m: {delta_cable/1e6:+,.0f} MN·s/m (cable reverses direction)")
        print(f"  Casing Δp/m: {delta_casing/1e6:+,.0f} MN·s/m")
        print(f"  Net Δp/m:    {(delta_cable + delta_casing)/1e6:+,.0f} MN·s/m")
        print(f"  External angular momentum source required (anchor lines to Earth,")
        print(f"  or phased deployment with early anchor drops).")

    # Sweep deployment progress
    print(f"\n  --- PRECESSION RATE DURING DEPLOYMENT ---")
    print(f"  J2 torque is constant: τ = {tau_j2:.4e} N·m")
    print(f"  (Depends on total mass and inclination, not on velocity)")
    print()

    fractions = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    print(f"  {'Prog':>6}  {'v_cable':>10}  {'v_casing':>10}  {'L_per_m':>14}  "
          f"{'L_total':>12}  {'Ω̇':>10}")
    print(f"  {'(%)':>6}  {'(m/s)':>10}  {'(m/s)':>10}  {'(kg·m/s/m)':>14}  "
          f"{'(kg·m²/s)':>12}  {'(°/day)':>10}")
    print(f"  " + "-" * 74)

    for f in fractions:
        v_c = V_ORBIT + f * (v_cable_final - V_ORBIT)
        v_cas = V_ORBIT + f * (v_casing_final - V_ORBIT)
        l_per_m = M_CABLE_PER_M * v_c + M_CASING_PER_M * v_cas
        L = l_per_m * R_ORBIT * C_RING
        if abs(L * math.sin(i_rad)) > 1e10:
            omega_prec = -tau_j2 / (L * math.sin(i_rad))
            omega_deg_day = omega_prec * 180 / PI * 86400
        else:
            omega_deg_day = float('inf') if L >= 0 else float('-inf')

        # Cap display for divergent values
        if abs(omega_deg_day) > 1e6:
            omega_str = f"{'DIVERGES':>10}"
        else:
            omega_str = f"{omega_deg_day:10.2f}"

        print(f"  {f*100:6.0f}  {v_c:10.1f}  {v_cas:10.1f}  {l_per_m:14.0f}  "
              f"{L:12.3e}  {omega_str}")

    if cable_dir == "prograde":
        print(f"\n  For prograde cable: L is conserved, Ω̇ stays constant at "
              f"{j2_precession_rate(i_rad) * 180/PI * 86400:.2f} °/day.")
        print(f"  Motor parameters (N, v_cable, etc.) do NOT affect the precession rate")
        print(f"  because the total angular momentum doesn't change.")
    else:
        print(f"\n  For retrograde cable: L passes through zero at ~50% deployment!")
        print(f"  At the zero-crossing, precession rate diverges (τ/L → ∞).")
        print(f"  After crossing, L becomes retrograde and precession reverses")
        print(f"  (prograde advance instead of regression).")
        print(f"  This zero-crossing is a critical stability point requiring careful")
        print(f"  management during retrograde deployment.")


# =============================================================================
# PHASE 3: STABILIZATION BEFORE ANCHOR DROP
# =============================================================================

def run_phase3(i_deg, cable_dir="prograde"):
    """Phase 3: Earth-tracking forces (computed here, applied by anchors in Phase 4).

    IMPORTANT: The LIM is an INTERNAL force between cable and casing. It cannot
    change the system's precession rate. Earth tracking requires an EXTERNAL
    torque, provided by anchor lines in Phase 4.

    During Phase 3 (free-floating after deployment, before anchor drop), the ring
    precesses freely at the J2 rate. The anchor drop must be timed so the ring's
    ascending node aligns with the ground anchor stations.

    This phase computes the forces that will be needed once anchored.
    """
    i_rad = math.radians(i_deg)

    v_cable = V_CABLE_PROGRADE if cable_dir == "prograde" else -V_CABLE_RETROGRADE

    # System angular momentum
    L_system = ring_angular_momentum(M_CABLE_PER_M, v_cable, M_CASING_PER_M, V_CASING)

    # J2 precession rate (for the coupled system)
    tau_j2 = j2_torque_on_ring(i_rad, M_RING_TOTAL)
    omega_j2 = -tau_j2 / (L_system * math.sin(i_rad))
    omega_j2_deg_day = omega_j2 * 180 / PI * 86400

    # Required: ascending node must track Earth's rotation
    # In inertial frame: Omega_dot = omega_earth
    # J2 gives: Omega_dot = omega_j2 (negative for prograde cable)
    # Correction needed: omega_correction = omega_earth - omega_j2
    omega_correction = OMEGA_EARTH - omega_j2
    omega_correction_deg_day = omega_correction * 180 / PI * 86400

    # Correction torque
    tau_correction = abs(omega_correction) * abs(L_system) * math.sin(i_rad)

    # Peak lateral ACCELERATION (sinusoidal distribution a(u) = a_max * sin(u))
    # Derivation: tau = a_max * M * r / 2 (from integrating a_max * sin^2(u) dm * r)
    # Therefore: a_max = 2 * tau / (M * r)
    a_lat_max = 2 * tau_correction / (M_RING_TOTAL * R_ORBIT)

    # Peak lateral FORCE per meter of ring
    F_lat_max_per_m = a_lat_max * M_TOTAL_PER_M

    # Total absolute lateral force (integrated around ring)
    # F_total = integral of F_peak * |sin(u)| * ds = F_peak_per_m * r * 4
    # Equivalently: F_total = 4 * tau / (pi * r)
    F_total = 4 * tau_correction / (PI * R_ORBIT)

    # Coriolis cross-check: a_lat = 2 * Omega_Earth * v_cable * sin(i)
    # This is the Coriolis acceleration on the cable in the rotating frame
    a_coriolis_cable = 2 * OMEGA_EARTH * abs(v_cable) * math.sin(i_rad)
    F_coriolis_cable = a_coriolis_cable * M_CABLE_PER_M
    a_coriolis_casing = 2 * OMEGA_EARTH * V_CASING * math.sin(i_rad)
    F_coriolis_casing = a_coriolis_casing * M_CASING_PER_M

    # Bearing loads (the bearing transmits the lateral force from anchor to cable)
    F_vert_bearing = M_CASING_PER_M * (G_ALT - V_CASING**2 / R_ORBIT)  # vertical (N/m)
    F_bearing_total = math.sqrt(F_lat_max_per_m**2 + F_vert_bearing**2)
    bearing_angle = math.degrees(math.atan2(F_lat_max_per_m, F_vert_bearing))

    print("\n" + "=" * 78)
    print(f"PHASE 3: EARTH-TRACKING FORCE ANALYSIS ({cable_dir.upper()})")
    print("=" * 78)

    print(f"\n  --- FREE-FLOATING PERIOD ---")
    print(f"\n  Before anchor drop, the ring is free-floating. The LIM is an")
    print(f"  INTERNAL force and cannot change the system's precession rate.")
    print(f"  The ring precesses freely at the J2 rate of {abs(omega_j2_deg_day):.2f} deg/day")
    print(f"  with a period of {360/abs(omega_j2_deg_day):.1f} days.")
    print(f"\n  The anchor drop must be timed so that the ascending node aligns")
    print(f"  with the pre-built ground anchor stations (temporal alignment).")

    print(f"\n  --- WHAT EARTH TRACKING REQUIRES ---")
    print(f"\n  The ring's ground track is fixed when the ascending node is")
    print(f"  fixed in the Earth-fixed (rotating) frame. In the inertial frame,")
    print(f"  this requires dΩ/dt = ω_Earth = {OMEGA_EARTH:.4e} rad/s")
    print(f"  ({OMEGA_EARTH * 180/PI * 86400:.2f} deg/day = 360 deg/sidereal day).")
    print(f"\n  J2 causes regression: dΩ/dt_J2 = {omega_j2_deg_day:+.4f} deg/day")
    print(f"  Required correction:  dΩ/dt_corr = {omega_correction_deg_day:+.4f} deg/day")
    print(f"\n  The correction has two components:")
    print(f"    1. Cancel J2 regression:  {abs(omega_j2_deg_day):.4f} deg/day")
    print(f"    2. Track Earth rotation:  {OMEGA_EARTH * 180/PI * 86400:.2f} deg/day")
    print(f"  Earth rotation dominates by {OMEGA_EARTH / abs(omega_j2) if omega_j2 != 0 else 0:.0f}x.")
    print(f"\n  This force must come from the ANCHOR LINES (external torque),")
    print(f"  not from the LIMs (internal force). The bearing transmits the")
    print(f"  lateral force from the casing/anchors to the cable.")

    print(f"\n  --- TORQUE AND FORCE ---")
    print(f"\n  System angular momentum: L = {L_system:.4e} kg*m^2/s")
    print(f"  J2 torque (on total mass): tau_J2 = {tau_j2:.4e} N*m")
    print(f"\n  Correction torque:")
    print(f"    tau_corr = dΩ/dt_corr x |L| x sin(i)")
    print(f"             = {abs(omega_correction):.4e} x {abs(L_system):.4e} x sin({i_deg} deg)")
    print(f"             = {tau_correction:.4e} N*m")

    print(f"\n  Lateral force distribution: F(u) = F_max x sin(u)")
    print(f"    Peak lateral acceleration: {a_lat_max:.4f} m/s^2 ({a_lat_max/9.81:.4f} g)")
    print(f"    Peak lateral force/m:      {F_lat_max_per_m/1000:.2f} kN/m")
    print(f"    Total lateral force:       {F_total/1e9:.3f} GN (absolute, both halves)")

    # --- Coriolis cross-check ---
    print(f"\n  --- CORIOLIS CROSS-CHECK ---")
    print(f"\n  In the Earth-fixed frame, the Coriolis acceleration on a mass")
    print(f"  moving at velocity v is a = 2 x Omega_Earth x v x sin(i).")
    print(f"\n  Cable (v = {abs(v_cable):.0f} m/s):")
    print(f"    a_Coriolis = 2 x {OMEGA_EARTH:.4e} x {abs(v_cable):.0f} x sin({i_deg} deg)")
    print(f"              = {a_coriolis_cable:.4f} m/s^2 ({a_coriolis_cable/9.81:.4f} g)")
    print(f"    F_Coriolis = {F_coriolis_cable/1000:.2f} kN/m")
    print(f"\n  Casing (v = {V_CASING:.0f} m/s):")
    print(f"    a_Coriolis = {a_coriolis_casing:.6f} m/s^2")
    print(f"    F_Coriolis = {F_coriolis_casing:.1f} N/m (negligible)")
    print(f"\n  Combined Coriolis: {(F_coriolis_cable + F_coriolis_casing)/1000:.2f} kN/m")
    print(f"  Torque-derived:    {F_lat_max_per_m/1000:.2f} kN/m")
    print(f"  (These differ slightly because the torque formula uses the full")
    print(f"   precession integral, while the Coriolis formula is a peak value.)")

    # --- Bearing load analysis ---
    print(f"\n  --- BEARING LOAD ANALYSIS ---")
    print(f"\n  The levitation bearing must transmit both vertical and lateral loads.")
    print(f"  Vertical (gravity support):   {F_vert_bearing/1000:.1f} kN/m")
    print(f"  Lateral (Earth tracking peak): {F_lat_max_per_m/1000:.1f} kN/m")
    print(f"  Resultant:                    {F_bearing_total/1000:.1f} kN/m")
    print(f"  Bearing angle from vertical:  {bearing_angle:.1f} deg")
    print(f"  Lateral/vertical ratio:       {F_lat_max_per_m/F_vert_bearing:.2f}")
    print(f"\n  The Chapter 5 bearing was designed for vertical loads only (109 kN/m).")
    if F_lat_max_per_m / F_vert_bearing > 0.1:
        print(f"  At {i_deg} deg inclination, the lateral load is {F_lat_max_per_m/F_vert_bearing*100:.0f}%")
        print(f"  of the vertical load. The bearing MUST be redesigned for 2D loading.")
    else:
        print(f"  At {i_deg} deg inclination, the lateral load is {F_lat_max_per_m/F_vert_bearing*100:.1f}%")
        print(f"  of the vertical load. This is within the margin of the existing design.")

    # --- Power analysis ---
    print(f"\n  --- POWER ANALYSIS ---")
    print(f"\n  The mechanical power for ideal gyroscopic precession is ZERO.")
    print(f"\n  Proof: the lateral force f(u) is proportional to sin(u) and the lateral")
    print(f"  velocity from precession v_lat(u) is proportional to cos(u). The power")
    print(f"  integral:")
    print(f"    P = integral of f(u) x v_lat(u) dm = f_0 x lambda x r^2 x integral sin(u)cos(u) du = 0")
    print(f"  (sin x cos integrates to zero over a full circle)")
    print(f"\n  The Coriolis force is always perpendicular to velocity, confirming")
    print(f"  P = F dot v = 0 at every point. No energy is needed for Earth tracking.")
    print(f"\n  The bearing transmits a large lateral force ({F_lat_max_per_m/1000:.1f} kN/m peak)")
    print(f"  but does zero net work, just like the vertical bearing force.")

    # --- Inclination sweep ---
    print(f"\n  --- EARTH-TRACKING FORCE VS INCLINATION ---")
    print(f"\n  {'i':>6}  {'a_lat':>10}  {'F_lat/m':>10}  {'F_vert/m':>10}  "
          f"{'Lat/Vert':>10}  {'Bearing':>10}")
    print(f"  {'(deg)':>6}  {'(m/s^2)':>10}  {'(kN/m)':>10}  {'(kN/m)':>10}  "
          f"{'ratio':>10}  {'angle(deg)':>10}")
    print(f"  " + "-" * 66)

    for i_d in range(0, 46, 2):
        i_r = math.radians(i_d)
        if i_d == 0:
            print(f"  {i_d:6d}  {0:10.4f}  {0:10.2f}  {F_vert_bearing/1000:10.1f}  "
                  f"{0:10.4f}  {0:10.1f}")
            continue
        tau_j = j2_torque_on_ring(i_r, M_RING_TOTAL)
        om_j = -tau_j / (L_system * math.sin(i_r))
        om_c = OMEGA_EARTH - om_j
        tau_c = abs(om_c) * abs(L_system) * math.sin(i_r)
        am = 2 * tau_c / (M_RING_TOTAL * R_ORBIT)
        flm = am * M_TOTAL_PER_M
        ratio = flm / F_vert_bearing
        angle = math.degrees(math.atan2(flm, F_vert_bearing))
        print(f"  {i_d:6d}  {am:10.4f}  {flm/1000:10.2f}  {F_vert_bearing/1000:10.1f}  "
              f"{ratio:10.4f}  {angle:10.1f}")

    return {
        'tau_correction': tau_correction,
        'a_lat_max': a_lat_max,
        'F_lat_max_per_m': F_lat_max_per_m,
        'F_total': F_total,
        'omega_correction': omega_correction,
        'L_system': L_system,
    }


# =============================================================================
# PHASE 4: ANCHOR LINES DEPLOYED
# =============================================================================

def run_phase4(i_deg, cable_dir="prograde", n_anchors=N_ANCHORS_DEFAULT):
    """Phase 4: J2 forces transmitted through anchor lines."""
    i_rad = math.radians(i_deg)

    v_cable = V_CABLE_PROGRADE if cable_dir == "prograde" else -V_CABLE_RETROGRADE

    # System angular momentum
    L_system = ring_angular_momentum(M_CABLE_PER_M, v_cable, M_CASING_PER_M, V_CASING)

    # Correction torque (same as Phase 3)
    tau_j2 = j2_torque_on_ring(i_rad, M_RING_TOTAL)
    omega_j2 = -tau_j2 / (L_system * math.sin(i_rad))
    omega_correction = OMEGA_EARTH - omega_j2
    tau_correction = abs(omega_correction) * abs(L_system) * math.sin(i_rad)
    a_lat_max = 2 * tau_correction / (M_RING_TOTAL * R_ORBIT)  # peak lateral acceleration (m/s^2)
    F_lat_max_per_m = a_lat_max * M_TOTAL_PER_M  # peak lateral force per meter (N/m)

    # Net vertical force per metre (supports payload + structural load)
    f_cable_up = M_CABLE_PER_M * (v_cable**2 / R_ORBIT - G_ALT)
    f_casing_down = M_CASING_PER_M * (G_ALT - V_CASING**2 / R_ORBIT)
    f_net_up = f_cable_up - f_casing_down

    # Per anchor station
    segment_length = C_RING / n_anchors
    F_vert_per_anchor = f_net_up * segment_length
    F_lat_per_anchor_peak = F_lat_max_per_m * segment_length  # CORRECTED: uses force/m, not accel

    # Anchor line cross-section (sized for vertical load)
    A_anchor = abs(F_vert_per_anchor) / SIGMA_OPERATING if F_vert_per_anchor > 0 else 1.0

    # Deflection angle from vertical
    theta_rad = math.atan2(F_lat_per_anchor_peak, abs(F_vert_per_anchor))
    theta_deg = math.degrees(theta_rad)
    theta_arcsec = theta_deg * 3600

    # Total tension with lateral component
    T_total = math.sqrt(F_vert_per_anchor**2 + F_lat_per_anchor_peak**2)
    sigma_total = T_total / A_anchor
    sigma_j2_only = F_lat_per_anchor_peak / A_anchor

    print("\n" + "=" * 78)
    print(f"PHASE 4: ANCHOR LINES DEPLOYED ({n_anchors} stations)")
    print("=" * 78)

    print(f"\n  Anchor stations:        {n_anchors}")
    print(f"  Segment per anchor:     {segment_length/1000:.1f} km")
    print(f"  Anchor line length:     {ANCHOR_LENGTH/1000:.0f} km")

    print(f"\n  --- VERTICAL FORCE BALANCE ---")
    print(f"\n  Cable centrifugal excess: {f_cable_up:,.0f} N/m")
    print(f"    (v_cable²/r = {v_cable**2/R_ORBIT:.2f} m/s², g = {G_ALT:.2f} m/s², "
          f"excess = {v_cable**2/R_ORBIT - G_ALT:.2f} m/s²)")
    print(f"  Casing gravity excess:    {f_casing_down:,.0f} N/m")
    print(f"    (g = {G_ALT:.2f} m/s², v_cas²/r = {V_CASING**2/R_ORBIT:.4f} m/s², "
          f"excess = {G_ALT - V_CASING**2/R_ORBIT:.2f} m/s²)")
    print(f"  Net upward force:         {f_net_up:,.0f} N/m ({f_net_up/1000:.1f} kN/m)")
    print(f"  Per anchor (vertical):    {F_vert_per_anchor/1e9:.3f} GN")

    print(f"\n  --- EARTH-TRACKING LATERAL FORCE ON ANCHORS ---")
    print(f"\n  Peak lateral accel:       {a_lat_max:.4f} m/s^2 ({a_lat_max/9.81:.4f} g)")
    print(f"  Peak lateral force/m:     {F_lat_max_per_m/1000:.2f} kN/m")
    print(f"  Peak lateral per anchor:  {F_lat_per_anchor_peak/1e6:.2f} MN ({F_lat_per_anchor_peak/1e9:.4f} GN)")

    print(f"\n  --- ANCHOR LINE STRESS ---")
    print(f"\n  Cross-section (sized for vertical): {A_anchor:.4f} m^2 "
          f"(d = {math.sqrt(4*A_anchor/PI)*1000:.0f} mm)")
    print(f"  Vertical stress:      {SIGMA_OPERATING/1e9:.1f} GPa (= operating limit by design)")
    print(f"  Lateral stress:       {sigma_j2_only/1e6:.2f} MPa ({sigma_j2_only/1e9:.4f} GPa)")
    print(f"  Total stress:         {sigma_total/1e9:.4f} GPa")
    print(f"  Lateral / operating:  {sigma_j2_only/SIGMA_OPERATING*100:.2f}%")
    if sigma_j2_only / SIGMA_OPERATING > 0.01:
        print(f"\n  The lateral stress is significant ({sigma_j2_only/SIGMA_OPERATING*100:.1f}% of limit).")
        print(f"  Anchor lines must be oversized to handle the combined load.")
    else:
        print(f"\n  The lateral stress is negligible for anchor line design.")

    print(f"\n  Deflection angle:     {theta_deg:.6f}° = {theta_arcsec:.3f} arcsec")
    print(f"  The anchor lines remain essentially vertical.")

    # --- Inclination sweep ---
    print(f"\n  --- ANCHOR FORCES VS INCLINATION ---")
    print(f"\n  {'i':>6}  {'F_lat/m':>10}  {'F_lat/anchor':>14}  {'theta':>10}  "
          f"{'sigma_lat':>12}  {'lat/op':>10}")
    print(f"  {'(deg)':>6}  {'(kN/m)':>10}  {'(MN)':>14}  {'(deg)':>10}  "
          f"{'(MPa)':>12}  {'(%)':>10}")
    print(f"  " + "-" * 68)

    for i_d in range(0, 46, 2):
        i_r = math.radians(i_d)
        if i_d == 0:
            print(f"  {i_d:6d}  {0:10.2f}  {0:14.2f}  {0:10.4f}  {0:12.2f}  {0:10.4f}")
            continue
        tau_j = j2_torque_on_ring(i_r, M_RING_TOTAL)
        om_j = -tau_j / (L_system * math.sin(i_r))
        om_c = OMEGA_EARTH - om_j
        tau_c = abs(om_c) * abs(L_system) * math.sin(i_r)
        am = 2 * tau_c / (M_RING_TOTAL * R_ORBIT)
        flm = am * M_TOTAL_PER_M  # force per meter (N/m)
        fl_a = flm * segment_length  # force per anchor (N)
        th_d = math.degrees(math.atan2(fl_a, abs(F_vert_per_anchor)))
        sig_j = fl_a / A_anchor
        print(f"  {i_d:6d}  {flm/1000:10.2f}  {fl_a/1e6:14.2f}  {th_d:10.4f}  "
              f"{sig_j/1e6:12.2f}  {sig_j/SIGMA_OPERATING*100:10.4f}")

    return {
        'F_vert_per_anchor': F_vert_per_anchor,
        'F_lat_per_anchor': F_lat_per_anchor_peak,
        'theta_deg': theta_deg,
        'sigma_j2': sigma_j2_only,
        'A_anchor': A_anchor,
    }


# =============================================================================
# PLOT FUNCTIONS
# =============================================================================

def generate_all_plots(i_deg, cable_dir, n_anchors):
    """Generate all 6 plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available — skipping plots")
        return

    graph_dir = os.path.dirname(os.path.abspath(__file__))
    i_rad = math.radians(i_deg)

    v_cable = V_CABLE_PROGRADE if cable_dir == "prograde" else -V_CABLE_RETROGRADE
    L_system = ring_angular_momentum(M_CABLE_PER_M, v_cable, M_CASING_PER_M, V_CASING)

    # ---- Plot 1: Precession rate vs inclination ----
    print("\n  Generating plot 1: precession_rate_vs_inclination.png ...")
    fig, ax = plt.subplots(figsize=(10, 6))
    inc = [x * 0.5 for x in range(1, 61)]  # 0.5 to 30 degrees
    omega_pt = [j2_precession_rate(math.radians(i)) * 180 / PI * 86400 for i in inc]
    omega_ring = []
    for i in inc:
        i_r = math.radians(i)
        tau = j2_torque_on_ring(i_r, M_CABLE_TOTAL)
        L = M_CABLE_TOTAL * V_ORBIT * R_ORBIT
        om = -tau / (L * math.sin(i_r))
        omega_ring.append(om * 180 / PI * 86400)

    ax.plot(inc, [abs(o) for o in omega_pt], 'b-', linewidth=2, label='Point mass')
    ax.plot(inc, [abs(o) for o in omega_ring], 'r--', linewidth=2, label='Ring (integrated)')
    ax.set_xlabel('Inclination (°)', fontsize=12)
    ax.set_ylabel('|Precession rate| (°/day)', fontsize=12)
    ax.set_title('J2 Nodal Precession Rate vs Inclination\n(250 km altitude, circular orbit)',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 30)
    fp = os.path.join(graph_dir, "precession_rate_vs_inclination.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    # ---- Plot 2: Ground track drift ----
    print("  Generating plot 2: ground_track_drift.png ...")
    fig, ax = plt.subplots(figsize=(10, 6))
    omega_prec = j2_precession_rate(i_rad)
    days = [d * 0.5 for d in range(61)]  # 0 to 30 days, 0.5 day steps
    omega_inertial = [omega_prec * 180 / PI * 86400 * d for d in days]
    ax.plot(days, omega_inertial, 'b-', linewidth=2)
    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel('Ascending node longitude, inertial (°)', fontsize=12)
    ax.set_title(f'Ground Track Drift — Free-Flying Ring at i = {i_deg}°\n'
                 f'(Ω̇ = {omega_prec * 180/PI * 86400:.2f} °/day)',
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 30)
    fp = os.path.join(graph_dir, "ground_track_drift.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    # ---- Plot 3: Deployment precession ----
    print("  Generating plot 3: deployment_precession.png ...")
    fig, ax = plt.subplots(figsize=(10, 6))

    tau_j2 = j2_torque_on_ring(i_rad, M_RING_TOTAL)
    fracs = [f / 200 for f in range(201)]

    for direction, ls, color, label in [
        ("prograde", '-', 'blue', 'Prograde cable'),
        ("retrograde", '--', 'red', 'Retrograde cable')
    ]:
        v_final = V_CABLE_PROGRADE if direction == "prograde" else -V_CABLE_RETROGRADE
        rates = []
        for f in fracs:
            v_c = V_ORBIT + f * (v_final - V_ORBIT)
            v_cas = V_ORBIT + f * (V_CASING - V_ORBIT)
            l_pm = M_CABLE_PER_M * v_c + M_CASING_PER_M * v_cas
            L = l_pm * R_ORBIT * C_RING
            if abs(L * math.sin(i_rad)) > 1e10:
                om = -tau_j2 / (L * math.sin(i_rad))
                om_dd = om * 180 / PI * 86400
            else:
                om_dd = float('nan')
            # Clamp for plotting
            if abs(om_dd) > 100:
                om_dd = float('nan')
            rates.append(om_dd)

        ax.plot([f * 100 for f in fracs], rates, ls, color=color,
                linewidth=2, label=label)

    ax.set_xlabel('Deployment Progress (%)', fontsize=12)
    ax.set_ylabel('Precession rate (°/day)', fontsize=12)
    ax.set_title(f'J2 Precession During Deployment (i = {i_deg}°)\n'
                 f'Retrograde: diverges at L → 0 (~50% deployment)',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 100)
    fp = os.path.join(graph_dir, "deployment_precession.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    # ---- Plot 4: Stabilization force vs inclination ----
    print("  Generating plot 4: stabilization_force_vs_inclination.png ...")
    fig, ax1 = plt.subplots(figsize=(10, 6))

    inc_stab = [x * 0.5 for x in range(1, 61)]
    f_maxes_kN = []
    taus = []
    for i_d in inc_stab:
        i_r = math.radians(i_d)
        tau_j = j2_torque_on_ring(i_r, M_RING_TOTAL)
        om_j = -tau_j / (L_system * math.sin(i_r))
        om_c = OMEGA_EARTH - om_j
        tau_c = abs(om_c) * abs(L_system) * math.sin(i_r)
        am = 2 * tau_c / (M_RING_TOTAL * R_ORBIT)
        flm = am * M_TOTAL_PER_M  # N/m
        f_maxes_kN.append(flm / 1000)  # kN/m
        taus.append(tau_c)

    ax1.plot(inc_stab, f_maxes_kN, 'b-', linewidth=2, label='Peak force/m')
    ax1.set_xlabel('Inclination (°)', fontsize=12)
    ax1.set_ylabel('Peak lateral force (kN/m)', fontsize=12, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_xlim(0, 30)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(inc_stab, [t / 1e18 for t in taus], 'r--', linewidth=2, label='Correction torque')
    ax2.set_ylabel('Correction torque (10¹⁸ N·m)', fontsize=12, color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=11, loc='upper left')

    ax1.set_title('Stabilization Force and Torque vs Inclination\n'
                  '(P_mechanical = 0 for ideal precession)', fontsize=13)
    fp = os.path.join(graph_dir, "stabilization_force_vs_inclination.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    # ---- Plot 5: Anchor force vs inclination ----
    print("  Generating plot 5: anchor_force_vs_inclination.png ...")
    fig, ax1 = plt.subplots(figsize=(10, 6))

    segment = C_RING / n_anchors
    f_cable_up = M_CABLE_PER_M * (v_cable**2 / R_ORBIT - G_ALT)
    f_casing_down = M_CASING_PER_M * (G_ALT - V_CASING**2 / R_ORBIT)
    f_net_up = f_cable_up - f_casing_down
    F_vert = f_net_up * segment

    fl_anchors = []
    thetas = []
    for i_d in inc_stab:
        i_r = math.radians(i_d)
        tau_j = j2_torque_on_ring(i_r, M_RING_TOTAL)
        om_j = -tau_j / (L_system * math.sin(i_r))
        om_c = OMEGA_EARTH - om_j
        tau_c = abs(om_c) * abs(L_system) * math.sin(i_r)
        am = 2 * tau_c / (M_RING_TOTAL * R_ORBIT)
        flm = am * M_TOTAL_PER_M  # N/m
        fl = flm * segment  # N per anchor
        fl_anchors.append(fl / 1e6)  # MN
        th = math.degrees(math.atan2(fl, abs(F_vert)))  # degrees
        thetas.append(th)

    ax1.plot(inc_stab, fl_anchors, 'b-', linewidth=2, label='Force per anchor (peak)')
    ax1.set_xlabel('Inclination (°)', fontsize=12)
    ax1.set_ylabel('Lateral force per anchor (MN)', fontsize=12, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_xlim(0, 30)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(inc_stab, thetas, 'r--', linewidth=2, label='Deflection angle')
    ax2.set_ylabel('Anchor deflection (deg)', fontsize=12, color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=11, loc='upper left')

    ax1.set_title(f'Earth-Tracking Lateral Force Per Anchor vs Inclination\n'
                  f'({n_anchors} stations, {segment/1000:.0f} km spacing)',
                  fontsize=13)
    fp = os.path.join(graph_dir, "anchor_force_vs_inclination.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    # ---- Plot 6: Anchor stress margin ----
    print("  Generating plot 6: anchor_stress_margin.png ...")
    fig, ax = plt.subplots(figsize=(10, 6))

    A_anchor = abs(F_vert) / SIGMA_OPERATING
    stresses = []
    for i_d in inc_stab:
        i_r = math.radians(i_d)
        tau_j = j2_torque_on_ring(i_r, M_RING_TOTAL)
        om_j = -tau_j / (L_system * math.sin(i_r))
        om_c = OMEGA_EARTH - om_j
        tau_c = abs(om_c) * abs(L_system) * math.sin(i_r)
        am = 2 * tau_c / (M_RING_TOTAL * R_ORBIT)
        flm = am * M_TOTAL_PER_M
        fl = flm * segment
        sig = fl / A_anchor
        stresses.append(sig / 1e6)  # MPa

    ax.semilogy(inc_stab, stresses, 'b-', linewidth=2, label='J2 lateral stress')
    ax.axhline(SIGMA_OPERATING / 1e6, color='red', ls='--', linewidth=2,
               label=f'Operating limit ({SIGMA_OPERATING/1e9:.0f} GPa)')
    ax.set_xlabel('Inclination (°)', fontsize=12)
    ax.set_ylabel('Stress (MPa)', fontsize=12)
    ax.set_title('Anchor Line J2 Lateral Stress vs Operating Limit\n'
                 '(log scale — J2 stress is negligible)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0, 30)
    ax.set_ylim(1e-4, 1e5)
    fp = os.path.join(graph_dir, "anchor_stress_margin.png")
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {fp}")

    print(f"\n  All 6 plots saved to: {graph_dir}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    i_deg = I_DEFAULT
    cable_dir = "prograde"
    n_anchors = N_ANCHORS_DEFAULT
    make_plots = False

    for arg in sys.argv[1:]:
        if arg.startswith("--inclination="):
            i_deg = float(arg.split("=")[1])
        elif arg.startswith("--cable="):
            cable_dir = arg.split("=")[1].lower()
        elif arg.startswith("--anchors="):
            n_anchors = int(arg.split("=")[1])
        elif arg == "--plots":
            make_plots = True
        elif arg == "--all":
            make_plots = True
        elif arg in ("--help", "-h"):
            print(__doc__)
            return

    if cable_dir not in ("prograde", "retrograde"):
        print(f"Error: --cable must be 'prograde' or 'retrograde', got '{cable_dir}'")
        return

    print("=" * 78)
    print("J2 PERTURBATION SIMULATION FOR AN INCLINED ORBITAL RING")
    print("=" * 78)
    print(f"\n  Inclination:    {i_deg}°")
    print(f"  Cable:          {cable_dir}")
    print(f"  Anchor stations: {n_anchors}")

    # Computed constants
    print(f"\n  --- ORBITAL PARAMETERS ---")
    print(f"  R_orbit:        {R_ORBIT:,} m ({R_ORBIT/1e6:.4f} Mm)")
    print(f"  v_orbit:        {V_ORBIT:,.1f} m/s")
    print(f"  Orbital period: {T_ORBIT:.1f} s ({T_ORBIT/60:.1f} min)")
    print(f"  Mean motion:    {N_ORBIT:.6e} rad/s")
    print(f"  Circumference:  {C_RING/1000:,.0f} km")
    print(f"  g at altitude:  {G_ALT:.3f} m/s²")
    print(f"  (R_E/r)²:       {(R_E/R_ORBIT)**2:.6f}")
    print(f"  J2:             {J2}")
    print(f"  ω_Earth:        {OMEGA_EARTH:.4e} rad/s "
          f"({OMEGA_EARTH * 180/PI * 86400:.2f} °/day)")

    # Run all four phases
    run_phase1(i_deg)
    run_phase2(i_deg, cable_dir)
    run_phase3(i_deg, cable_dir)
    run_phase4(i_deg, cable_dir, n_anchors)

    # Summary
    print("\n" + "=" * 78)
    print("SUMMARY OF FINDINGS")
    print("=" * 78)

    i_rad = math.radians(i_deg)
    omega_prec = j2_precession_rate(i_rad)
    v_cable = V_CABLE_PROGRADE if cable_dir == "prograde" else -V_CABLE_RETROGRADE
    L_system = ring_angular_momentum(M_CABLE_PER_M, v_cable, M_CASING_PER_M, V_CASING)
    tau_j2 = j2_torque_on_ring(i_rad, M_RING_TOTAL)
    omega_correction = OMEGA_EARTH - (-tau_j2 / (L_system * math.sin(i_rad)))
    tau_correction = abs(omega_correction) * abs(L_system) * math.sin(i_rad)
    a_lat_max = 2 * tau_correction / (M_RING_TOTAL * R_ORBIT)
    F_lat_max_per_m = a_lat_max * M_TOTAL_PER_M
    segment = C_RING / n_anchors
    f_cable_up = M_CABLE_PER_M * (v_cable**2 / R_ORBIT - G_ALT)
    f_casing_down = M_CASING_PER_M * (G_ALT - V_CASING**2 / R_ORBIT)
    f_net_up = f_cable_up - f_casing_down
    F_vert = f_net_up * segment
    F_lat = F_lat_max_per_m * segment
    F_vert_bearing_pm = M_CASING_PER_M * (G_ALT - V_CASING**2 / R_ORBIT)
    lat_vert_ratio = F_lat_max_per_m / F_vert_bearing_pm

    print(f"""
  At inclination i = {i_deg} deg:

  1. PRECESSION: The ring precesses at {abs(omega_prec * 180/PI * 86400):.2f} deg/day (westward).
     This is IDENTICAL to a point-mass satellite at the same altitude.
     Ring integration confirms: the result is exact, not approximate.

  2. DEPLOYMENT: For prograde cable, angular momentum is conserved by the
     LIMs. The precession rate stays constant throughout deployment.
     For retrograde cable, L passes through zero at ~50% deployment,
     causing a precession singularity requiring external stabilization.

  3. EARTH TRACKING: To hold the ground track fixed, the correction torque is
     {tau_correction:.3e} N*m, dominated by Earth-tracking ({OMEGA_EARTH * 180/PI * 86400:.0f} deg/day)
     rather than J2 cancellation ({abs(omega_prec * 180/PI * 86400):.1f} deg/day).
     Peak lateral acceleration: {a_lat_max:.4f} m/s^2 ({a_lat_max/9.81:.4f} g)
     Peak lateral force: {F_lat_max_per_m/1000:.1f} kN/m
     Mechanical power: ZERO (Coriolis force perpendicular to velocity).

  4. BEARING LOAD: The bearing must handle lateral AND vertical forces.
     Vertical: {F_vert_bearing_pm/1000:.1f} kN/m (gravity support)
     Lateral:  {F_lat_max_per_m/1000:.1f} kN/m (Earth tracking, peak)
     Ratio:    {lat_vert_ratio:.2f} (lateral/vertical)
     {"BEARING REDESIGN REQUIRED for 2D loading." if lat_vert_ratio > 0.1 else "Within existing design margin."}

  5. ANCHORS: Peak lateral force per anchor = {F_lat/1e6:.2f} MN.
     Vertical load per anchor = {F_vert/1e9:.2f} GN.
     Lateral / vertical = {F_lat/abs(F_vert)*100:.2f}%
     Anchor deflection: {math.degrees(math.atan2(F_lat, abs(F_vert))):.2f} deg from vertical.
     Lateral stress: {F_lat / (abs(F_vert)/SIGMA_OPERATING) / 1e6:.1f} MPa vs 7 GPa limit.

  CONCLUSION: The J2 perturbation force itself is small, but the Earth-tracking
  force (Coriolis) is significant. At {i_deg} deg inclination, the bearing must
  handle a lateral load of {lat_vert_ratio*100:.0f}% of vertical. The force requires
  no energy (P = 0) but demands a bearing designed for 2D loading and anchor
  lines that can transmit lateral tension.
""")

    print("=" * 78)

    # Generate plots if requested
    if make_plots:
        generate_all_plots(i_deg, cable_dir, n_anchors)


if __name__ == "__main__":
    main()
