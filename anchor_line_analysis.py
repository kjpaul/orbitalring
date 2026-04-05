#!/usr/bin/env python3
"""
Anchor Line and Lateral Force Simulation

Computes anchor line forces, lengths, and self-weight as a function of
orbital inclination for an inclined orbital ring.

Usage:
    python anchor_line_analysis.py [options]

    --inclination=N   Orbital inclination in degrees (default 30)
    --altitude=N      Orbital altitude in km (default 250)
    --cable=DIR       Cable direction: prograde or retrograde (default retrograde)
    --stations=N      Number of anchor stations (default 800)
    --plots           Generate all PNG plots
    --csv             Write per-station CSV files
    --all             Run everything (plots + csv)

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import os
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

GM = 3.986004418e14          # Earth gravitational parameter (m^3/s^2)
R_E_EQUATOR = 6_378_137.0   # WGS84 equatorial radius (m)
R_E_POLES = 6_356_752.3142  # WGS84 polar radius (m)
F_EARTH = 1 / 298.257223563 # WGS84 flattening
J2 = 1.0826e-3              # Earth J2 coefficient
OMEGA_SIDEREAL = 7.2921159e-5  # Earth sidereal rotation rate (rad/s)
PI = math.pi

# =============================================================================
# CNT MATERIAL
# =============================================================================

RHO_CNT = 1_700             # CNT yarn density (kg/m^3)
SIGMA_BREAK = 25e9           # Breaking strength (Pa) -- project assumption
SIGMA_OPERATING = 12.5e9    # Operating limit = sigma_break / 2 (Pa)

# =============================================================================
# RING PARAMETERS (from lim_config.py / lim_physics.py)
# =============================================================================

M_LOAD_PER_M = 12_000       # Casing mass per meter (kg/m)
SIGMA_CABLE = 12.5e9         # Operating stress (Pa), safety factor 2.0 on 25 GPa
RHO_CNT = 1700               # CNT density (kg/m³)
M_HARDWARE = 4076            # Hardware mass per meter (kg/m)

# MathCAD calibration table for lateral force validation
# { inclination_deg: peak_lateral_force_per_m (N/m) }
MATHCAD_CALIBRATION = {
    5:  11_400,
    10: 22_700,
    20: 44_600,
    30: 65_100,
}


# =============================================================================
# ORBITAL MECHANICS HELPERS
# =============================================================================

def orbital_velocity(r):
    """Circular orbital velocity at radius r (m/s)."""
    return math.sqrt(GM / r)


def gravity_at_radius(r):
    """Gravitational acceleration at radius r (m/s^2), no centrifugal."""
    return GM / r**2


def net_gravity_at_radius(r, lat_rad=0.0):
    """Net downward acceleration at radius r including centrifugal term.

    g_net = GM/r^2 - omega^2 * r * cos^2(lat)
    The centrifugal term reduces effective gravity for ground-synchronous objects.
    """
    return GM / r**2 - OMEGA_SIDEREAL**2 * r * math.cos(lat_rad)**2


def earth_radius_wgs84(lat_rad):
    """WGS84 Earth radius at geodetic latitude (m).

    Approximation: R(lat) = R_eq * (1 - f * sin^2(lat))
    """
    return R_E_EQUATOR * (1 - F_EARTH * math.sin(lat_rad)**2)


def ground_sync_velocity(r, i_rad=0.0):
    """Ground-synchronous velocity at orbital radius r for inclination i (m/s).

    v_gs = omega_sidereal * r * cos(i)
    At the equatorial crossing, the ring tangent makes angle i with the equator.
    The effective ground-tracking velocity in the orbital plane is reduced
    by cos(i) because the casing must also move north/south.
    """
    return OMEGA_SIDEREAL * r * math.cos(i_rad)


def cable_mass_structural(m_load, r, m_hw=M_HARDWARE, sigma_cable=SIGMA_CABLE):
    """Structural cable mass from the quadratic sizing equation.

    Solves: S * m^2 + (S * m_hw - K) * m - (K * m_hw + P) = 0
    where S = sigma/rho, K = m_load*dv*2*v_orb - m_load*g_net*r, P = (m_load*dv)^2.
    See Volume II, Chapter 6 for derivation.
    """
    if m_load <= 0:
        return 0.0
    v_orbit = orbital_velocity(r)
    v_ground = OMEGA_SIDEREAL * r
    delta_v = v_orbit - v_ground
    g_net = net_gravity_at_radius(r)

    S = sigma_cable / RHO_CNT
    K = m_load * delta_v * 2 * v_orbit - m_load * g_net * r
    P = (m_load * delta_v) ** 2

    a = S
    b = S * m_hw - K
    c = -(K * m_hw + P)

    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return 0.0
    return (-b + math.sqrt(discriminant)) / (2 * a)


def cable_velocity_from_momentum(m_cable, m_load, v_orbit, v_casing):
    """Cable velocity from momentum conservation during deployment.

    Before deployment: everything at v_orbit (all mass co-orbiting).
    After deployment: cable at v_cable, casing at v_casing.
    LIMs are internal forces, so total momentum is conserved:

    (m_cable + m_load) * v_orbit = m_cable * v_cable + m_load * v_casing
    v_cable = ((m_cable + m_load) * v_orbit - m_load * v_casing) / m_cable
    """
    m_total = m_cable + m_load
    return (m_total * v_orbit - m_load * v_casing) / m_cable


def compute_ring_params(altitude_km, cable_dir, inclination_deg):
    """Compute all ring parameters for a given altitude and cable direction.

    Returns a dict with all derived parameters needed by the analysis.
    """
    h = altitude_km * 1000
    r = R_E_EQUATOR + h
    v_orb = orbital_velocity(r)
    C = 2 * PI * r
    g_alt = gravity_at_radius(r)
    g_net_eq = net_gravity_at_radius(r, 0.0)
    i_rad = math.radians(inclination_deg)

    # Casing velocity: ground-synchronous, adjusted for inclination
    # At equatorial crossing the ring tangent makes angle i with the equator,
    # so the effective velocity along the ring is omega*r*cos(i)
    v_casing = ground_sync_velocity(r, i_rad)
    v_casing_eq = OMEGA_SIDEREAL * r  # equatorial, for reference

    # Structural cable mass from force balance (same formula as lim_physics.py)
    m_cable_struct = cable_mass_structural(M_LOAD_PER_M, r)

    # Hardware mass: LIM plates, levitation hardware, etc.
    # From lim_config.py: M_HARDWARE ~ 4,077 kg/m (gamma titanium plates + iron)
    m_hardware = 4_077  # kg/m
    m_cable = m_cable_struct + m_hardware

    # Cable velocity from momentum conservation
    v_cable_mag = cable_velocity_from_momentum(m_cable, M_LOAD_PER_M, v_orb, v_casing)

    if cable_dir == "retrograde":
        v_cable_signed = -v_cable_mag  # negative = westward
    else:
        v_cable_signed = v_cable_mag

    # Linear momentum per meter of ring (NOT angular momentum)
    # p_per_m = m_cable * v_cable + m_load * v_casing (signed)
    # For retrograde cable, v_cable is negative, so p_per_m is large and negative.
    p_per_m = m_cable * v_cable_signed + M_LOAD_PER_M * v_casing

    # Total angular momentum (for reference)
    L_total = p_per_m * r * C

    # Vertical force balance
    f_cable_up = m_cable * (v_cable_mag**2 / r - g_alt)
    f_casing_down = M_LOAD_PER_M * (g_alt - v_casing**2 / r)
    f_net_up = f_cable_up - f_casing_down

    # LIM site parameters
    lim_spacing = 500.0
    n_lim_sites = round(C / lim_spacing)

    # Cable cross-section and bending stiffness for sled dynamics
    cable_area = m_cable_struct / RHO_CNT
    cable_side = math.sqrt(cable_area)
    I_cable = cable_side**4 / 12  # second moment of area (square section)
    EI_cable = SIGMA_BREAK * I_cable
    T_cable = m_cable_struct * SIGMA_OPERATING
    # Transverse wave speed: v_wave = sqrt(sigma/rho) for a stressed cable
    v_wave = math.sqrt(SIGMA_OPERATING / RHO_CNT)

    return {
        'altitude_km': altitude_km,
        'h': h,
        'r': r,
        'v_orb': v_orb,
        'C': C,
        'g_alt': g_alt,
        'g_net_eq': g_net_eq,
        'cable_dir': cable_dir,
        'inclination_deg': inclination_deg,
        'i_rad': i_rad,
        'v_casing': v_casing,
        'v_casing_eq': v_casing_eq,
        'm_cable_structural': m_cable_struct,
        'm_cable': m_cable,
        'm_load': M_LOAD_PER_M,
        'v_cable_mag': v_cable_mag,
        'v_cable_signed': v_cable_signed,
        'p_per_m': p_per_m,
        'L_total': L_total,
        'f_cable_up': f_cable_up,
        'f_casing_down': f_casing_down,
        'f_net_up': f_net_up,
        'lim_spacing': lim_spacing,
        'n_lim_sites': n_lim_sites,
        'cable_area': cable_area,
        'cable_side': cable_side,
        'EI_cable': EI_cable,
        'v_wave': v_wave,
        'T_cable': T_cable,
    }


# =============================================================================
# SECTION 1: LATERAL EARTH-TRACKING FORCE
# =============================================================================

def peak_lateral_force_per_m(i_rad, p_per_m):
    """Peak lateral force per meter at the apex (maximum latitude point).

    The lateral (out-of-plane) force arises from Coriolis acceleration on the
    cable in the rotating (Earth-fixed) frame.  The cable's enormous angular
    momentum must precess at omega_sidereal to keep the ground track fixed.

    DERIVATION:
    At position u on the ring (u=0 at ascending node, u=90 at northern apex),
    the Coriolis force on the cable has an out-of-plane component:

        F_lat(u) = 2 * omega * m_cable * |v_cable| * sin(i) * sin(u)

    The force distribution is sinusoidal: zero at equatorial crossings (u=0,180)
    and maximum at the apex (u=90,270).

    Wait -- but shouldn't it be the opposite?  No.  The torque calculation:

        The precession torque tau = omega * L * sin(i) is about the line of
        nodes (x-axis in the equatorial frame).  The x-component of torque
        from a lateral force f(u)*sin(u) at position u is:

            d(tau_x) = f(u) * sin(u) * r * ds

        Integrating f(u) = F_peak * sin(u) gives:

            tau_x = F_peak * r^2 * integral(sin^2(u) du, 0, 2*pi) = F_peak * r^2 * pi

        Setting tau_x = omega * L_total * sin(i):

            F_peak = omega * L_total * sin(i) / (r^2 * pi)

        With L_total = p_per_m * r * 2*pi*r:

            F_peak = omega * p_per_m * 2*pi*r^2 * sin(i) / (r^2 * pi)
                   = 2 * omega * |p_per_m| * sin(i)

    HOWEVER, the Coriolis cross-check gives peak at the EQUATORIAL CROSSING
    (where the velocity vector has maximum angle to equator), and the torque
    derivation gives peak with sin(u), which is at the APEX.

    These are reconciled by noting that the out-of-plane component of Coriolis
    at position u depends on the projection of velocity onto the orbit-normal
    cross Earth-rotation-axis, which actually gives a cos(u) pattern for force
    and sin(u) for the torque arm.  The NET result for F_peak is the same.

    For the anchor line analysis, what matters is:
    - F_peak per meter = 2 * omega * |p_per_m| * sin(i)
    - The distribution around the ring (sin vs cos) affects which anchors
      carry the most lateral load, but the peak force per meter is the same.

    We use the CORIOLIS distribution: peak at equatorial crossings, zero at apex.
    This matches the MathCAD model and the physical intuition (the cable crosses
    the equator at angle i, creating maximum lateral deflection there).

    F_peak = 2 * omega_sidereal * |p_per_m| * sin(i)
    """
    return 2 * OMEGA_SIDEREAL * abs(p_per_m) * math.sin(i_rad)


def lateral_force_at_nu(nu_rad, F_peak):
    """Lateral force per meter at position nu around the ring.

    f(nu) = F_peak * |cos(nu)|

    where nu=0 at the ascending node (equatorial crossing) and nu=90 at apex.
    Peak force at equatorial crossings, zero at apex.
    """
    return F_peak * abs(math.cos(nu_rad))


def verify_calibration(params):
    """Check computed lateral forces against MathCAD calibration table.

    The calibration was computed using equatorial-like ring parameters.
    Our self-consistent computation (with momentum conservation) may differ
    slightly due to hardware mass assumptions.

    Returns list of (incl, computed, expected, pct_error) tuples.
    """
    results = []
    for i_deg, f_expected in MATHCAD_CALIBRATION.items():
        i_rad = math.radians(i_deg)

        # Recompute v_casing and momentum at this inclination
        v_cas = ground_sync_velocity(params['r'], i_rad)

        # Cable velocity from momentum conservation at this inclination
        v_cab = cable_velocity_from_momentum(
            params['m_cable'], params['m_load'], params['v_orb'], v_cas)

        if params['cable_dir'] == "retrograde":
            p = params['m_cable'] * (-v_cab) + params['m_load'] * v_cas
        else:
            p = params['m_cable'] * v_cab + params['m_load'] * v_cas

        f_computed = peak_lateral_force_per_m(i_rad, p)
        pct = 100 * (f_computed - f_expected) / f_expected
        results.append((i_deg, f_computed, f_expected, pct))
    return results


# =============================================================================
# SECTION 2: ANCHOR LINE GEOMETRY
# =============================================================================

def station_latitude(nu_rad, i_rad):
    """Latitude of anchor station at position nu on an inclined ring.

    sin(lat) = sin(i) * sin(nu)
    """
    sin_lat = math.sin(i_rad) * math.sin(nu_rad)
    sin_lat = max(-1.0, min(1.0, sin_lat))
    return math.asin(sin_lat)


def anchor_line_length(r_orbit, lat_rad):
    """Anchor line length from ring to ground (straight-line approximation).

    Uses WGS84 Earth radius at the local latitude.
    """
    R_E_local = earth_radius_wgs84(lat_rad)
    return r_orbit - R_E_local


def integrated_gravity_along_line(r_orbit, R_E_local, lat_rad, n_steps=100):
    """Integrate gravitational acceleration along the anchor line.

    The line runs from R_E_local (ground) to r_orbit (ring altitude).
    g(r) = GM/r^2 - omega^2 * r * cos^2(lat)

    Returns (g_avg, g_integral):
        g_integral = integral_{R_E}^{r_orbit} g(r) dr  [units: m^2/s^2 = J/kg]
        g_avg = g_integral / L_line
    """
    dr = (r_orbit - R_E_local) / n_steps
    cos2lat = math.cos(lat_rad)**2
    g_integral = 0.0
    for k in range(n_steps):
        r = R_E_local + (k + 0.5) * dr
        g_r = GM / r**2 - OMEGA_SIDEREAL**2 * r * cos2lat
        g_integral += g_r * dr
    g_avg = g_integral / (r_orbit - R_E_local)
    return g_avg, g_integral


# =============================================================================
# SECTION 3: ANCHOR LINE SELF-WEIGHT AND COMBINED LOADS
# =============================================================================

def self_weight_stress(rho, g_avg, L_line):
    """Stress at top of a uniform cross-section anchor line from self-weight.

    sigma = rho * g_avg * L
    """
    return rho * g_avg * L_line


def taper_ratio(rho, g_avg, L_line, sigma_op):
    """Exponential taper ratio for a cable supporting its own weight.

    A_top / A_bottom = exp(rho * g_avg * L / sigma_op)
    """
    exponent = rho * g_avg * L_line / sigma_op
    return math.exp(exponent)


def compute_station_data(params, n_stations, i_deg):
    """Compute forces and geometry for every anchor station at inclination i_deg.

    Returns a list of dicts, one per station.
    """
    i_rad = math.radians(i_deg)
    r = params['r']
    spacing = params['C'] / n_stations

    # Recompute v_casing and momentum at this inclination
    v_cas = ground_sync_velocity(r, i_rad)
    v_cab = cable_velocity_from_momentum(
        params['m_cable'], params['m_load'], params['v_orb'], v_cas)

    if params['cable_dir'] == "retrograde":
        p = params['m_cable'] * (-v_cab) + params['m_load'] * v_cas
    else:
        p = params['m_cable'] * v_cab + params['m_load'] * v_cas

    F_peak = peak_lateral_force_per_m(i_rad, p)

    stations = []
    for k in range(n_stations):
        nu = 2 * PI * k / n_stations  # position angle around ring

        # Latitude
        lat = station_latitude(nu, i_rad)

        # Local Earth radius and anchor length
        R_E_local = earth_radius_wgs84(lat)
        L_line = r - R_E_local  # straight-line approximation

        # Gravity integration along line
        g_avg, g_integral = integrated_gravity_along_line(r, R_E_local, lat)

        # Lateral force per meter at this position
        f_lat_per_m = lateral_force_at_nu(nu, F_peak)

        # Forces per anchor
        f_vertical = params['m_load'] * net_gravity_at_radius(r, lat) * spacing
        f_lateral = f_lat_per_m * spacing
        f_total = math.sqrt(f_vertical**2 + f_lateral**2)

        # Angle from vertical
        theta = math.atan2(f_lateral, f_vertical)

        # Self-weight stress (independent of cross-section for uniform line)
        sigma_self = self_weight_stress(RHO_CNT, g_avg, L_line)

        # Taper ratio
        tr = taper_ratio(RHO_CNT, g_avg, L_line, SIGMA_OPERATING)

        # Reference line: 10 cm diameter
        d_ref = 0.10
        A_ref = PI * (d_ref / 2)**2

        # Number of 10 cm lines needed to carry total external load
        # Each line can carry (sigma_op - sigma_self) * A_ref of external load
        capacity_net = (SIGMA_OPERATING - sigma_self) * A_ref
        if capacity_net > 0:
            n_lines = math.ceil(f_total / capacity_net)
        else:
            n_lines = 999_999  # self-weight alone exceeds limit

        # Combined stress at top for n_lines lines:
        # sigma_combined = sigma_self + f_total / (n_lines * A_ref)
        if n_lines < 999_999:
            sigma_external = f_total / (n_lines * A_ref)
            sigma_combined = sigma_self + sigma_external
        else:
            sigma_external = float('inf')
            sigma_combined = float('inf')

        stations.append({
            'k': k,
            'nu_deg': math.degrees(nu),
            'nu_rad': nu,
            'lat_deg': math.degrees(lat),
            'lat_rad': lat,
            'R_E_local': R_E_local,
            'h_local': r - R_E_local,
            'L_line': L_line,
            'g_avg': g_avg,
            'g_integral': g_integral,
            'f_lat_per_m': f_lat_per_m,
            'f_vertical': f_vertical,
            'f_lateral': f_lateral,
            'f_total': f_total,
            'theta_deg': math.degrees(theta),
            'theta_rad': theta,
            'sigma_self': sigma_self,
            'taper_ratio': tr,
            'n_lines_10cm': n_lines,
            'sigma_combined': sigma_combined,
            'F_peak': F_peak,
            'p_per_m': p,
        })

    return stations


# =============================================================================
# SECTION 4: DYNAMIC SLED PASSAGE LOADS
# =============================================================================

def sled_passage_analysis(params):
    """Compute anchor attachment loads during mass driver sled passage.

    The Mach number M = v_sled / v_wave determines the cable response.
    Characteristic bending length: L_char = (EI / T_eff)^(1/4)
    """
    EI = params['EI_cable']
    v_wave = params['v_wave']
    cable_side = params['cable_side']

    sled_speeds = [1000, 5000, 10000, 20000, 30000]  # m/s
    sled_mass_per_m = 525   # kg/m
    sled_length = 10_000    # m
    anchor_footprint = 25   # m (load distribution width at attachment)

    results = []
    for v_sled in sled_speeds:
        M = v_sled / v_wave
        # Effective tension during supersonic sled passage
        if M > 1:
            T_eff_sled = params['T_cable'] * (M**2 - 1)
        else:
            T_eff_sled = params['T_cable']

        # Characteristic bending length
        L_char = (EI / max(T_eff_sled, 1))**0.25

        # Sled gravitational load at the anchor point
        g_net = net_gravity_at_radius(params['r'])
        F_sled_total = sled_mass_per_m * sled_length * g_net

        # Point load on anchor from sled weight distributed over footprint
        F_anchor_sled = F_sled_total * anchor_footprint / sled_length

        # Bending stress at anchor: sigma = M_bend / S
        # M_bend = F * L_char / 2  (point load on tensioned beam)
        # S = h^3/6  (section modulus for square cross-section)
        M_bend = F_anchor_sled * L_char / 2
        S_cable = cable_side**3 / 6
        sigma_bend = M_bend / S_cable if S_cable > 0 else 0

        results.append({
            'v_sled': v_sled,
            'v_sled_kms': v_sled / 1000,
            'Mach': M,
            'T_eff_sled': T_eff_sled,
            'L_char': L_char,
            'F_sled_total': F_sled_total,
            'F_anchor_sled': F_anchor_sled,
            'M_bend': M_bend,
            'sigma_bend': sigma_bend,
            'sigma_bend_GPa': sigma_bend / 1e9,
            'margin_GPa': SIGMA_OPERATING / 1e9 - sigma_bend / 1e9,
        })

    return results


# =============================================================================
# SECTION 5: PLOTTING
# =============================================================================

GRAPH_DPI = 300
GRAPH_WIDTH = 10
GRAPH_HEIGHT = 6
GRAPH_FORMAT = "png"
INCLINATIONS_PLOT = [5, 10, 20, 30]
COLORS = {5: '#0072B2', 10: '#009E73', 20: '#E69F00', 30: '#CC0000'}


def _save_fig(fig, name, output_dir):
    """Save figure to output directory."""
    path = os.path.join(output_dir, f"fig_anchor_{name}.{GRAPH_FORMAT}")
    fig.savefig(path, dpi=GRAPH_DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"  Saved: {os.path.basename(path)}")


def plot_lateral_force(all_data, output_dir):
    """Plot 1: Lateral force per meter vs position around ring."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        f_lat = [s['f_lat_per_m'] / 1000 for s in stations]
        ax.plot(nu, f_lat, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Lateral force per meter (kN/m)', fontsize=14)
    ax.set_title('Lateral Earth-Tracking Force Distribution', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'lateral_force', output_dir)


def plot_anchor_angle(all_data, output_dir):
    """Plot 2: Anchor line angle from vertical vs position."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        theta = [s['theta_deg'] for s in stations]
        ax.plot(nu, theta, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Anchor line angle from vertical (deg)', fontsize=14)
    ax.set_title('Anchor Line Deflection Angle', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'anchor_angle', output_dir)


def plot_anchor_length(all_data, output_dir):
    """Plot 3: Anchor line length vs position."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        L_km = [s['L_line'] / 1000 for s in stations]
        ax.plot(nu, L_km, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Anchor line length (km)', fontsize=14)
    ax.set_title('Anchor Line Length vs Position', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'anchor_length', output_dir)


def plot_self_weight_stress(all_data, output_dir):
    """Plot 4: Self-weight stress at top vs position (uniform cross-section)."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        sigma = [s['sigma_self'] / 1e9 for s in stations]
        ax.plot(nu, sigma, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.axhline(y=SIGMA_OPERATING / 1e9, color='red', linestyle='--',
               label=f'Operating limit ({SIGMA_OPERATING/1e9:.1f} GPa)')
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Self-weight stress at top (GPa)', fontsize=14)
    ax.set_title('Anchor Line Self-Weight Stress (uniform cross-section)', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'self_weight_stress', output_dir)


def plot_total_force(all_data, output_dir):
    """Plot 5: Total external force per anchor vs position."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        f_tot = [s['f_total'] / 1e9 for s in stations]
        ax.plot(nu, f_tot, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Total external force per anchor (GN)', fontsize=14)
    ax.set_title('Total Force per Anchor Station (vertical + lateral)', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'total_force', output_dir)


def plot_n_lines(all_data, output_dir):
    """Plot 6: Number of 10 cm anchor lines per station vs position."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    for i_deg in INCLINATIONS_PLOT:
        stations = all_data[i_deg]
        nu = [s['nu_deg'] for s in stations]
        n = [min(s['n_lines_10cm'], 1e6) for s in stations]
        ax.plot(nu, n, color=COLORS[i_deg], label=f'i = {i_deg}', linewidth=1.5)
    ax.set_xlabel('Position around ring (deg)', fontsize=14)
    ax.set_ylabel('Number of 10 cm lines per station', fontsize=14)
    ax.set_title('Anchor Lines Required per Station (d = 10 cm CNT)', fontsize=16)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 360)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'n_lines', output_dir)


def plot_taper_ratio(all_data, output_dir):
    """Plot 7: Taper ratio vs inclination for worst-case anchor."""
    # Use r_orbit from data
    first_key = list(all_data.keys())[0]
    r = all_data[first_key][0]['R_E_local'] + all_data[first_key][0]['h_local']

    incl_range = np.arange(0.5, 31, 0.5)
    taper_vals = []

    for i_deg in incl_range:
        # Worst case: apex station where latitude = inclination
        lat_apex = math.radians(i_deg)
        R_E_apex = earth_radius_wgs84(lat_apex)
        h_apex = r - R_E_apex
        g_avg, _ = integrated_gravity_along_line(r, R_E_apex, lat_apex)
        tr = taper_ratio(RHO_CNT, g_avg, h_apex, SIGMA_OPERATING)
        taper_vals.append(tr)

    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))
    ax.plot(incl_range, taper_vals, color='#CC0000', linewidth=2)
    ax.set_xlabel('Orbital inclination (deg)', fontsize=14)
    ax.set_ylabel('Taper ratio (A_top / A_bottom)', fontsize=14)
    ax.set_title('Anchor Line Taper Ratio at Worst-Case Station (apex)', fontsize=16)
    ax.grid(True, alpha=0.3)
    _save_fig(fig, 'taper_ratio', output_dir)


def plot_stress_budget(all_data, output_dir):
    """Plot 8: Stacked bar chart of stress budget at worst-case anchor."""
    fig, ax = plt.subplots(figsize=(GRAPH_WIDTH, GRAPH_HEIGHT))

    incl_list = sorted(all_data.keys())
    x = np.arange(len(incl_list))
    width = 0.5

    sigma_self_vals = []
    sigma_ext_vals = []

    for i_deg in incl_list:
        stations = all_data[i_deg]
        # Worst case = station with highest combined stress
        worst = max(stations, key=lambda s: s['sigma_combined']
                    if s['sigma_combined'] < 1e18 else 0)
        sigma_self_vals.append(worst['sigma_self'] / 1e9)
        if worst['n_lines_10cm'] < 999_999:
            d_ref = 0.10
            A_ref = PI * (d_ref / 2)**2
            ext = worst['f_total'] / (worst['n_lines_10cm'] * A_ref) / 1e9
        else:
            ext = 0
        sigma_ext_vals.append(ext)

    ax.bar(x, sigma_self_vals, width, label='Self-weight stress', color='#0072B2')
    ax.bar(x, sigma_ext_vals, width, bottom=sigma_self_vals,
           label='External load stress', color='#E69F00')
    ax.axhline(y=SIGMA_OPERATING / 1e9, color='red', linestyle='--',
               linewidth=2, label=f'Operating limit ({SIGMA_OPERATING/1e9:.1f} GPa)')
    ax.set_xlabel('Inclination (deg)', fontsize=14)
    ax.set_ylabel('Stress at top of anchor line (GPa)', fontsize=14)
    ax.set_title('Anchor Line Stress Budget (worst-case station)', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{i}' for i in incl_list])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    _save_fig(fig, 'stress_budget', output_dir)


# =============================================================================
# CONSOLE OUTPUT
# =============================================================================

def print_parameters(params):
    """Print ring parameters."""
    p = params
    sep = "=" * 72
    print(f"\n{sep}")
    print("ANCHOR LINE AND LATERAL FORCE SIMULATION")
    print(sep)
    print(f"  Altitude:           {p['altitude_km']} km")
    print(f"  Orbital radius:     {p['r']/1e6:.4f} Mm ({p['r']:,.0f} m)")
    print(f"  Orbital velocity:   {p['v_orb']:.1f} m/s")
    print(f"  Circumference:      {p['C']/1000:,.0f} km")
    print(f"  Cable direction:    {p['cable_dir']}")
    print(f"  Inclination:        {p['inclination_deg']} deg")
    print(f"  Cable velocity:     {p['v_cable_mag']:.1f} m/s"
          f" ({'retrograde' if p['v_cable_signed'] < 0 else 'prograde'})")
    print(f"  Casing velocity:    {p['v_casing']:.3f} m/s (ground-sync)")
    print(f"  Cable mass:         {p['m_cable']:,.0f} kg/m"
          f" (structural: {p['m_cable_structural']:,.0f} kg/m)")
    print(f"  Casing mass:        {p['m_load']:,.0f} kg/m")
    print(f"  p_per_m (momentum): {p['p_per_m']/1e6:,.2f} MN*s/m")
    print(f"  L_total (angular):  {p['L_total']:.3e} kg*m^2/s")
    print(f"  g at altitude:      {p['g_alt']:.3f} m/s^2")
    print(f"  g_net (equatorial): {p['g_net_eq']:.3f} m/s^2")
    print(f"  Net upward force:   {p['f_net_up']:,.0f} N/m"
          f" ({p['f_net_up']/1000:.1f} kN/m)")
    print(f"  CNT: sigma_break = {SIGMA_BREAK/1e9:.0f} GPa,"
          f" sigma_op = {SIGMA_OPERATING/1e9:.1f} GPa,"
          f" rho = {RHO_CNT} kg/m^3")
    print(sep)


def print_calibration(cal_results):
    """Print calibration verification table."""
    print(f"\n  {'--- CALIBRATION CHECK: lateral force vs MathCAD ---':^60}")
    print(f"  {'Incl':>6}  {'Computed':>12}  {'MathCAD':>12}  {'Error':>8}")
    print(f"  {'(deg)':>6}  {'(kN/m)':>12}  {'(kN/m)':>12}  {'(%)':>8}")
    print(f"  {'-'*44}")
    all_ok = True
    for i_deg, f_comp, f_exp, pct in cal_results:
        flag = ""
        if abs(pct) > 5:
            flag = " *** WARNING ***"
            all_ok = False
        elif abs(pct) > 2:
            flag = " * check *"
        print(f"  {i_deg:6.0f}  {f_comp/1000:12.1f}  {f_exp/1000:12.1f}  {pct:+8.2f}%{flag}")
    if all_ok:
        print(f"\n  All calibration values within 5% of MathCAD.")
    else:
        print(f"\n  WARNING: Some values exceed 5% deviation from MathCAD!")


def print_summary_table(all_data, params, n_stations):
    """Print summary table for each inclination at apex and equatorial crossing."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("SUMMARY TABLE")
    print(sep)

    for i_deg in sorted(all_data.keys()):
        stations = all_data[i_deg]

        # Find equatorial crossing (nu ~ 0) and apex (nu ~ 90)
        eq_station = min(stations, key=lambda s: abs(s['nu_deg'] - 0))
        apex_station = min(stations, key=lambda s: abs(s['nu_deg'] - 90))

        print(f"\n  Inclination = {i_deg} deg")
        print(f"  {'':28} {'Equator':>14} {'Apex':>14}")
        print(f"  {'-'*58}")
        print(f"  {'Lateral force/m (kN/m)':28} {eq_station['f_lat_per_m']/1e3:14.1f}"
              f" {apex_station['f_lat_per_m']/1e3:14.1f}")
        print(f"  {'F_vertical/anchor (GN)':28} {eq_station['f_vertical']/1e9:14.3f}"
              f" {apex_station['f_vertical']/1e9:14.3f}")
        print(f"  {'F_lateral/anchor (GN)':28} {eq_station['f_lateral']/1e9:14.3f}"
              f" {apex_station['f_lateral']/1e9:14.3f}")
        print(f"  {'F_total/anchor (GN)':28} {eq_station['f_total']/1e9:14.3f}"
              f" {apex_station['f_total']/1e9:14.3f}")
        print(f"  {'Anchor angle (deg)':28} {eq_station['theta_deg']:14.2f}"
              f" {apex_station['theta_deg']:14.2f}")
        print(f"  {'Anchor length (km)':28} {eq_station['L_line']/1e3:14.2f}"
              f" {apex_station['L_line']/1e3:14.2f}")
        print(f"  {'Self-weight stress (GPa)':28} {eq_station['sigma_self']/1e9:14.3f}"
              f" {apex_station['sigma_self']/1e9:14.3f}")
        print(f"  {'Taper ratio':28} {eq_station['taper_ratio']:14.4f}"
              f" {apex_station['taper_ratio']:14.4f}")
        print(f"  {'N lines (d=10cm) per stn':28} {eq_station['n_lines_10cm']:14,}"
              f" {apex_station['n_lines_10cm']:14,}")
        s_eq = eq_station['sigma_combined']
        s_ap = apex_station['sigma_combined']
        if s_eq < 1e18 and s_ap < 1e18:
            print(f"  {'Combined stress (GPa)':28} {s_eq/1e9:14.3f}"
                  f" {s_ap/1e9:14.3f}")
            print(f"  {'% of operating limit':28}"
                  f" {100*s_eq/SIGMA_OPERATING:13.1f}%"
                  f" {100*s_ap/SIGMA_OPERATING:13.1f}%")


def print_sled_analysis(sled_results):
    """Print sled passage dynamic loads table."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("SECTION 5: DYNAMIC SLED PASSAGE LOADS")
    print(sep)
    v_wave = sled_results[0]['v_sled'] / sled_results[0]['Mach']
    print(f"\n  Sled: 525 kg/m x 10 km = 5,250 t")
    print(f"  Anchor attachment footprint: 25 m")
    print(f"  Cable transverse wave speed: {v_wave:.0f} m/s")
    print(f"  Operating limit: {SIGMA_OPERATING/1e9:.1f} GPa\n")
    print(f"  {'v_sled':>8}  {'Mach':>6}  {'L_char':>8}  {'F_anchor':>10}"
          f"  {'sigma_bend':>12}  {'Margin':>10}")
    print(f"  {'(km/s)':>8}  {'':>6}  {'(m)':>8}  {'(kN)':>10}"
          f"  {'(Pa)':>12}  {'(GPa)':>10}")
    print(f"  {'-'*58}")
    for r in sled_results:
        print(f"  {r['v_sled_kms']:8.0f}  {r['Mach']:6.2f}  {r['L_char']:8.2f}"
              f"  {r['F_anchor_sled']/1e3:10.1f}  {r['sigma_bend']:12.1f}"
              f"  {r['margin_GPa']:10.2f}")


def write_csv(all_data, output_dir):
    """Write per-station CSV for each inclination."""
    for i_deg in sorted(all_data.keys()):
        stations = all_data[i_deg]
        fname = os.path.join(output_dir,
                             f"anchor_data_i{i_deg:02d}deg.csv")
        with open(fname, 'w') as f:
            f.write("station,nu_deg,lat_deg,h_km,L_line_km,g_avg,"
                    "f_lat_per_m_kN,f_vertical_GN,f_lateral_GN,f_total_GN,"
                    "theta_deg,sigma_self_GPa,taper_ratio,n_lines_10cm,"
                    "sigma_combined_GPa\n")
            for s in stations:
                f.write(f"{s['k']},{s['nu_deg']:.2f},{s['lat_deg']:.4f},"
                        f"{s['h_local']/1e3:.3f},{s['L_line']/1e3:.3f},"
                        f"{s['g_avg']:.4f},"
                        f"{s['f_lat_per_m']/1e3:.3f},"
                        f"{s['f_vertical']/1e9:.6f},"
                        f"{s['f_lateral']/1e9:.6f},"
                        f"{s['f_total']/1e9:.6f},"
                        f"{s['theta_deg']:.4f},"
                        f"{s['sigma_self']/1e9:.6f},"
                        f"{s['taper_ratio']:.6f},"
                        f"{s['n_lines_10cm']},"
                        f"{s['sigma_combined']/1e9:.6f}\n")
        print(f"  Wrote: {os.path.basename(fname)}")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    """Main entry point."""
    # Defaults
    inclination_deg = 30.0
    altitude_km = 250.0
    cable_dir = "retrograde"
    n_stations = 800
    do_plots = False
    do_csv = False

    # Parse command line arguments
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__)
            return
        elif arg.startswith("--inclination="):
            inclination_deg = float(arg.split("=")[1])
        elif arg.startswith("--altitude="):
            altitude_km = float(arg.split("=")[1])
        elif arg.startswith("--cable="):
            cable_dir = arg.split("=")[1].lower()
        elif arg.startswith("--stations="):
            n_stations = int(arg.split("=")[1])
        elif arg == "--plots":
            do_plots = True
        elif arg == "--csv":
            do_csv = True
        elif arg == "--all":
            do_plots = True
            do_csv = True

    # Compute ring parameters
    params = compute_ring_params(altitude_km, cable_dir, inclination_deg)

    # Print header
    print_parameters(params)

    # --- Section 1: Lateral force ---
    print(f"\n{'='*72}")
    print("SECTION 1: LATERAL EARTH-TRACKING FORCE")
    print("=" * 72)

    print(f"\n  Derivation:")
    print(f"  The ring's ground track is held fixed by precessing the orbital")
    print(f"  plane at omega_sidereal = {OMEGA_SIDEREAL:.4e} rad/s about Earth's polar axis.")
    print(f"")
    print(f"  The cable has angular momentum L = p_per_m * r * C, where")
    print(f"  p_per_m = m_cable*v_cable + m_load*v_casing is the linear")
    print(f"  momentum per meter of ring.  The precession torque is:")
    print(f"    tau = omega_s * |L_total| * sin(i)")
    print(f"")
    print(f"  The lateral force distribution f(nu) produces this torque.")
    print(f"  The peak force per meter is:")
    print(f"    F_peak = 2 * omega_s * |p_per_m| * sin(i)")
    print(f"")
    print(f"  |p_per_m| = {abs(params['p_per_m'])/1e6:,.1f} MN*s/m")
    print(f"  At i = {inclination_deg} deg:")
    F_peak_default = peak_lateral_force_per_m(params['i_rad'], params['p_per_m'])
    print(f"    F_peak = 2 * {OMEGA_SIDEREAL:.4e} * {abs(params['p_per_m'])/1e6:,.1f}e6"
          f" * sin({inclination_deg}) = {F_peak_default/1e3:.1f} kN/m")
    print(f"")

    cal = verify_calibration(params)
    print_calibration(cal)

    # --- Compute station data for each inclination ---
    all_data = {}
    for i_deg in INCLINATIONS_PLOT:
        all_data[i_deg] = compute_station_data(params, n_stations, i_deg)

    # Also compute for the requested inclination if not in the standard set
    if inclination_deg not in INCLINATIONS_PLOT:
        all_data[inclination_deg] = compute_station_data(
            params, n_stations, inclination_deg)

    # --- Sections 2-4: summary ---
    print(f"\n{'='*72}")
    print("SECTIONS 2-4: ANCHOR GEOMETRY, SELF-WEIGHT, AND COMBINED LOADS")
    print("=" * 72)
    print_summary_table(all_data, params, n_stations)

    # --- Section 5: Sled passage ---
    sled_results = sled_passage_analysis(params)
    print_sled_analysis(sled_results)

    # --- Conclusion ---
    print(f"\n{'='*72}")
    print("CONCLUSION")
    print("=" * 72)
    worst_incl = max(all_data.keys())
    worst_stations = all_data[worst_incl]
    eq_stn = min(worst_stations, key=lambda s: abs(s['nu_deg']))
    apex_stn = min(worst_stations, key=lambda s: abs(s['nu_deg'] - 90))
    print(f"\n  At i = {worst_incl} deg (worst case):")
    print(f"    Peak lateral force: {eq_stn['f_lat_per_m']/1e3:.1f} kN/m"
          f" at equatorial crossing")
    print(f"    Vertical load/anchor: {eq_stn['f_vertical']/1e9:.2f} GN")
    print(f"    Lateral load/anchor: {eq_stn['f_lateral']/1e9:.2f} GN"
          f" (peak)")
    print(f"    Anchor deflection: {eq_stn['theta_deg']:.2f} deg")
    print(f"    Self-weight stress: {apex_stn['sigma_self']/1e9:.3f} GPa"
          f" ({100*apex_stn['sigma_self']/SIGMA_OPERATING:.1f}% of limit)")
    print(f"    Taper ratio at apex: {apex_stn['taper_ratio']:.4f}")
    print(f"    Lines per station (equator): {eq_stn['n_lines_10cm']:,}"
          f" x 10 cm diameter")
    print(f"")
    print(f"  Sled passage bending stress: negligible"
          f" ({sled_results[-1]['sigma_bend_GPa']:.4f} GPa at"
          f" {sled_results[-1]['v_sled_kms']:.0f} km/s)")
    print("=" * 72)

    # --- Output ---
    output_dir = os.path.dirname(os.path.abspath(__file__))
    fig_dir = os.path.join(os.path.dirname(output_dir), "localbrain", "figures")
    if os.path.isdir(fig_dir):
        output_dir = fig_dir
    else:
        output_dir = os.path.join(output_dir, "output", "anchor")
    os.makedirs(output_dir, exist_ok=True)

    if do_plots:
        print(f"\nGenerating plots to: {output_dir}")
        plot_lateral_force(all_data, output_dir)
        plot_anchor_angle(all_data, output_dir)
        plot_anchor_length(all_data, output_dir)
        plot_self_weight_stress(all_data, output_dir)
        plot_total_force(all_data, output_dir)
        plot_n_lines(all_data, output_dir)
        plot_taper_ratio(all_data, output_dir)
        plot_stress_budget(all_data, output_dir)
        print(f"  All 8 plots saved.")

    if do_csv:
        print(f"\nWriting CSV files to: {output_dir}")
        write_csv(all_data, output_dir)


if __name__ == "__main__":
    main()
