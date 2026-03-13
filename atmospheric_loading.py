#!/usr/bin/env python3
"""
Atmospheric Loading on Ground-Anchored Tethers
================================================

Section 4: Computes wind drag loading, lateral deflection, and thermal
environment for a ground-anchored space elevator tether passing through
Earth's atmosphere (0-100 km).

The ground-anchored tether is STATIONARY relative to Earth's surface.
Both the tether and the lower atmosphere co-rotate with Earth, so
there is no orbital-velocity aerodynamic heating. Wind loading from
natural atmospheric winds is the dominant lateral force.

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python atmospheric_loading.py                 Print results only
    python atmospheric_loading.py all             Generate all graphs
    python atmospheric_loading.py drag_loading    Generate drag loading graph
    python atmospheric_loading.py deflection      Generate lateral deflection graph
    python atmospheric_loading.py temperature     Generate temperature profile graph
    python atmospheric_loading.py --show          Display interactively
    python atmospheric_loading.py --save          Save graphs to files (default)
    python atmospheric_loading.py --area=1e-3     Override bottom cross-section area (m^2)
    python atmospheric_loading.py --help, -h      Show this help

Graph keywords:
    drag_loading    Wind drag per meter vs altitude with atmospheric context
    deflection      Lateral deflection and deflection angle vs altitude
    temperature     Equilibrium tether temperature vs altitude
    all             Generate all graphs
"""

import sys
import os
import math
import csv

import matplotlib.pyplot as plt

import tether_config as cfg


# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

C_D_CYLINDER = 1.2          # Drag coefficient for cylinder in crossflow
C_P_AIR = 1005.0            # Specific heat capacity of air (J/(kg·K))
SIGMA_SB = 5.670374419e-8   # Stefan-Boltzmann constant (W/(m^2·K^4))
SOLAR_FLUX = 1361.0          # Solar constant (W/m^2)
ALPHA_SOLAR = 0.9            # Solar absorptivity (dark CNT)
EPSILON_IR = 0.9             # Infrared emissivity (CNT)
H_ATM_TOP = 100_000.0       # Top of significant atmosphere (m)


# =============================================================================
# NET ACCELERATION AND TAPER INTEGRAL
# =============================================================================

def a_net(r):
    """Net downward acceleration at radius r (positive toward Earth).

    a_net(r) = GM/r^2 - omega^2 * r
    """
    return cfg.GM / r**2 - cfg.OMEGA_SIDEREAL**2 * r


def phi_integral(r_bottom, r_top):
    """Analytical integral of a_net(r) from r_bottom to r_top.

    Integral of (GM/r^2 - omega^2 * r) dr:
        GM * (1/r_bottom - 1/r_top) - omega^2 * (r_top^2 - r_bottom^2) / 2
    """
    omega = cfg.OMEGA_SIDEREAL
    return (cfg.GM * (1.0 / r_bottom - 1.0 / r_top)
            - omega**2 * (r_top**2 - r_bottom**2) / 2.0)


# =============================================================================
# TETHER CROSS-SECTION AND DIAMETER
# =============================================================================

def tether_area(h, A_bottom=1e-4):
    """Tether cross-section at altitude h using uniform-stress taper profile.

    A(h) = A_bottom * exp(rho_CNT / sigma_operating * Phi(r_bottom, r))

    where Phi is the integral of a_net from the ground to radius R_E + h.

    Args:
        h: Altitude above Earth surface (m)
        A_bottom: Cross-section area at ground level (m^2), default 1 cm^2

    Returns:
        Cross-section area (m^2)
    """
    r_bottom = cfg.R_E
    r = cfg.R_E + h
    if r <= r_bottom:
        return A_bottom
    Phi = phi_integral(r_bottom, r)
    return A_bottom * math.exp(cfg.RHO_CNT / cfg.SIGMA_OPERATING * Phi)


def tether_diameter(h, A_bottom=1e-4):
    """Tether diameter at altitude h for given bottom cross-section area.

    d(h) = sqrt(4 * A(h) / pi)

    Args:
        h: Altitude above Earth surface (m)
        A_bottom: Cross-section area at ground level (m^2), default 1 cm^2

    Returns:
        Tether diameter (m)
    """
    A = tether_area(h, A_bottom)
    return math.sqrt(4.0 * A / math.pi)


# =============================================================================
# WIND DRAG
# =============================================================================

def wind_drag_per_meter(h, A_bottom=1e-4, C_D=C_D_CYLINDER):
    """Lateral wind drag force per meter of tether at altitude h.

    f_drag = 0.5 * rho_atm * v_wind^2 * C_D * d(h)

    The tether is stationary relative to Earth's surface. The only
    aerodynamic force comes from natural atmospheric winds.

    Args:
        h: Altitude (m)
        A_bottom: Bottom cross-section area (m^2)
        C_D: Drag coefficient for cylinder in crossflow

    Returns:
        Drag force per unit length (N/m)
    """
    T_atm, p_atm, rho_atm = cfg.atmosphere(h)
    v_wind = cfg.wind_speed(h)
    d = tether_diameter(h, A_bottom)
    return 0.5 * rho_atm * v_wind**2 * C_D * d


# =============================================================================
# TETHER TENSION
# =============================================================================

def tether_tension(h, A_bottom=1e-4):
    """Tether tension at altitude h.

    In the uniform-stress design, the stress is sigma_operating everywhere.
    The tension at altitude h is:
        T(h) = sigma_operating * A(h)

    Args:
        h: Altitude (m)
        A_bottom: Bottom cross-section area (m^2)

    Returns:
        Tension (N)
    """
    A = tether_area(h, A_bottom)
    return cfg.SIGMA_OPERATING * A


# =============================================================================
# LATERAL DEFLECTION
# =============================================================================

def lateral_deflection_profile(A_bottom=1e-4, n_points=1000):
    """Compute lateral deflection vs altitude due to wind loading.

    At each altitude, the local deflection angle is:
        theta(h) = f_drag(h) / (T(h) per unit length)
                 = f_drag(h) * 1m / T(h)

    The cumulative lateral deflection is:
        x(h) = integral(0 to h) tan(theta(h')) dh'

    For small angles, tan(theta) ~ theta.

    Args:
        A_bottom: Bottom cross-section area (m^2)
        n_points: Number of integration points

    Returns:
        (h_array, deflection_array, angle_array)
        h_array: altitudes (m)
        deflection_array: cumulative lateral deflection (m)
        angle_array: local deflection angle (rad)
    """
    h_array = []
    deflection_array = []
    angle_array = []

    dh = H_ATM_TOP / (n_points - 1)
    cumulative_x = 0.0

    for i in range(n_points):
        h = i * dh
        f_drag = wind_drag_per_meter(h, A_bottom)
        T = tether_tension(h, A_bottom)

        # Deflection angle: ratio of lateral force per meter to tension
        if T > 0:
            theta = f_drag / T  # radians (small angle)
        else:
            theta = 0.0

        # Accumulate lateral displacement
        if i > 0:
            cumulative_x += math.tan(theta) * dh

        h_array.append(h)
        deflection_array.append(cumulative_x)
        angle_array.append(theta)

    return h_array, deflection_array, angle_array


# =============================================================================
# THERMAL ENVIRONMENT
# =============================================================================

def equilibrium_temperature(h):
    """Equilibrium temperature of tether at altitude h.

    The ground-anchored tether is stationary relative to the atmosphere,
    so there is NO aerodynamic heating from the tether's own velocity.
    Wind-induced stagnation temperature rise is negligible (< 2 K at
    jet stream speeds).

    Heating sources:
        - Solar radiation (absorbed on projected area = d * 1m)
        - Atmospheric convection (significant in lower atmosphere)

    Cooling:
        - Infrared radiation (emitted from full circumference = pi * d * 1m)
        - Atmospheric convection (can also cool in lower atmosphere)

    Above ~30 km, convection is negligible and temperature is set by
    radiative equilibrium:
        T_eq = (alpha_s * S / (pi * epsilon * sigma_SB))^(1/4)

    In the lower atmosphere, convective cooling brings the tether closer
    to the local air temperature.

    Args:
        h: Altitude (m)

    Returns:
        Equilibrium temperature (K)
    """
    T_atm, p_atm, rho_atm = cfg.atmosphere(h)
    v_wind = cfg.wind_speed(h)

    # Radiative equilibrium temperature (no convection)
    T_rad = (ALPHA_SOLAR * SOLAR_FLUX / (math.pi * EPSILON_IR * SIGMA_SB)) ** 0.25

    # Convective heat transfer coefficient estimate (forced convection on cylinder)
    # Using simplified Hilpert correlation: h_conv ~ k/d * C * Re^m * Pr^(1/3)
    # For simplicity, use a representative convective coefficient
    # At sea level in moderate wind: h_conv ~ 50-200 W/(m^2·K)
    # Scale with density and wind speed: h_conv ~ h_conv_ref * (rho/rho_ref)^0.5 * (v/v_ref)^0.6

    rho_ref = 1.225   # sea-level density (kg/m^3)
    v_ref = 10.0      # reference wind speed (m/s)
    h_conv_ref = 80.0  # reference convective coefficient (W/(m^2·K))

    if rho_atm > 1e-6 and v_wind > 0.1:
        h_conv = h_conv_ref * (rho_atm / rho_ref) ** 0.5 * (v_wind / v_ref) ** 0.6
    else:
        h_conv = 0.0

    # Stagnation temperature rise from wind (negligible at wind speeds)
    delta_T_stag = v_wind**2 / (2.0 * C_P_AIR)

    # Energy balance per unit projected area (width d, length 1m):
    # Solar absorbed: alpha_s * S (per unit projected area, on one side)
    # Convective: h_conv * pi (full circumference/projected width) * (T_atm - T)
    #   (cooling when T > T_atm, heating when T < T_atm)
    # Radiative emission: epsilon * sigma_SB * T^4 * pi (full circumference)
    #
    # Balance: alpha_s * S + h_conv * pi * T_atm = epsilon * sigma_SB * pi * T^4 + h_conv * pi * T
    # Rearranging: epsilon * sigma_SB * pi * T^4 + h_conv * pi * T = alpha_s * S + h_conv * pi * T_atm

    # Solve iteratively
    T_eq = T_rad  # initial guess
    for _ in range(50):
        q_solar = ALPHA_SOLAR * SOLAR_FLUX
        q_conv_in = h_conv * math.pi * (T_atm + delta_T_stag)
        q_rad_out = EPSILON_IR * SIGMA_SB * math.pi * T_eq**4
        q_conv_out = h_conv * math.pi * T_eq

        # f(T) = q_rad_out + q_conv_out - q_solar - q_conv_in = 0
        f = q_rad_out + q_conv_out - q_solar - q_conv_in
        # f'(T) = 4 * epsilon * sigma_SB * pi * T^3 + h_conv * pi
        f_prime = 4.0 * EPSILON_IR * SIGMA_SB * math.pi * T_eq**3 + h_conv * math.pi

        if abs(f_prime) < 1e-30:
            break
        T_eq = T_eq - f / f_prime
        if T_eq < 50.0:
            T_eq = 50.0

    return T_eq


# =============================================================================
# FULL ATMOSPHERIC PROFILE
# =============================================================================

def atmospheric_profile(A_bottom=1e-4, n_points=1000):
    """Compute full atmospheric loading profile from ground to 100 km.

    Returns:
        Dictionary with arrays:
            h:           altitude (m)
            rho:         air density (kg/m^3)
            wind:        wind speed (m/s)
            drag:        drag force per meter (N/m)
            tension:     tether tension (N)
            diameter:    tether diameter (m)
            area:        tether cross-section (m^2)
            deflection:  cumulative lateral deflection (m)
            angle:       local deflection angle (rad)
            temperature: equilibrium temperature (K)
    """
    h_arr, defl_arr, angle_arr = lateral_deflection_profile(A_bottom, n_points)

    rho_arr = []
    wind_arr = []
    drag_arr = []
    tension_arr = []
    diameter_arr = []
    area_arr = []
    temp_arr = []

    for h in h_arr:
        T_atm, p_atm, rho_atm = cfg.atmosphere(h)
        v_wind = cfg.wind_speed(h)
        f_drag = wind_drag_per_meter(h, A_bottom)
        T_ten = tether_tension(h, A_bottom)
        d = tether_diameter(h, A_bottom)
        A = tether_area(h, A_bottom)
        T_eq = equilibrium_temperature(h)

        rho_arr.append(rho_atm)
        wind_arr.append(v_wind)
        drag_arr.append(f_drag)
        tension_arr.append(T_ten)
        diameter_arr.append(d)
        area_arr.append(A)
        temp_arr.append(T_eq)

    return {
        'h': h_arr,
        'rho': rho_arr,
        'wind': wind_arr,
        'drag': drag_arr,
        'tension': tension_arr,
        'diameter': diameter_arr,
        'area': area_arr,
        'deflection': defl_arr,
        'angle': angle_arr,
        'temperature': temp_arr,
    }


# =============================================================================
# GRAPH: DRAG LOADING
# =============================================================================

def plot_drag_loading(profile):
    """Plot wind drag per meter vs altitude with atmospheric context.

    Three subplots: wind speed, air density (log), and drag force per meter.
    """
    h_km = [h / 1000.0 for h in profile['h']]

    fig, axes = plt.subplots(1, 3, figsize=(cfg.GRAPH_WIDTH_INCHES * 1.5,
                                             cfg.GRAPH_HEIGHT_INCHES))

    # Subplot 1: Wind speed
    ax1 = axes[0]
    ax1.plot(profile['wind'], h_km, 'b-', linewidth=1.5)
    ax1.set_xlabel('Wind speed (m/s)', fontsize=12)
    ax1.set_ylabel('Altitude (km)', fontsize=12)
    ax1.set_title('Wind Profile', fontsize=13)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 100)

    # Subplot 2: Air density (log scale)
    ax2 = axes[1]
    ax2.semilogx(profile['rho'], h_km, 'r-', linewidth=1.5)
    ax2.set_xlabel('Air density (kg/m³)', fontsize=12)
    ax2.set_title('Atmospheric Density', fontsize=13)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_ylim(0, 100)

    # Subplot 3: Drag force per meter
    ax3 = axes[2]
    ax3.plot([d * 1000 for d in profile['drag']], h_km, 'g-', linewidth=1.5)
    ax3.set_xlabel('Drag per meter (mN/m)', fontsize=12)
    ax3.set_title('Wind Drag Loading', fontsize=13)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 100)

    fig.suptitle('Atmospheric Drag Loading on Ground-Anchored Tether',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"01-drag_loading.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: LATERAL DEFLECTION
# =============================================================================

def plot_deflection(profile):
    """Plot lateral deflection and deflection angle vs altitude."""
    h_km = [h / 1000.0 for h in profile['h']]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(cfg.GRAPH_WIDTH_INCHES * 1.2,
                                                    cfg.GRAPH_HEIGHT_INCHES))

    # Subplot 1: Cumulative lateral deflection
    ax1.plot(profile['deflection'], h_km, 'b-', linewidth=1.5)
    ax1.set_xlabel('Lateral deflection (m)', fontsize=12)
    ax1.set_ylabel('Altitude (km)', fontsize=12)
    ax1.set_title('Cumulative Lateral Deflection', fontsize=13)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 100)

    # Subplot 2: Local deflection angle
    angle_urad = [a * 1e6 for a in profile['angle']]
    ax2.plot(angle_urad, h_km, 'r-', linewidth=1.5)
    ax2.set_xlabel('Deflection angle (urad)', fontsize=12)
    ax2.set_title('Local Deflection Angle', fontsize=13)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 100)

    fig.suptitle('Wind-Induced Lateral Deflection of Ground-Anchored Tether',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"02-lateral_deflection.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# GRAPH: TEMPERATURE PROFILE
# =============================================================================

def plot_temperature(profile):
    """Plot equilibrium tether temperature vs altitude."""
    h_km = [h / 1000.0 for h in profile['h']]
    T_celsius = [t - 273.15 for t in profile['temperature']]

    # Also get atmospheric temperature for comparison
    T_atm_celsius = []
    for h in profile['h']:
        T_atm, _, _ = cfg.atmosphere(h)
        T_atm_celsius.append(T_atm - 273.15)

    fig, ax = plt.subplots(figsize=(cfg.GRAPH_WIDTH_INCHES,
                                     cfg.GRAPH_HEIGHT_INCHES))

    ax.plot(T_celsius, h_km, 'r-', linewidth=1.5, label='Tether equilibrium')
    ax.plot(T_atm_celsius, h_km, 'b--', linewidth=1.0, alpha=0.7,
            label='Atmospheric temperature')
    ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.5)

    ax.set_xlabel('Temperature (deg C)', fontsize=12)
    ax.set_ylabel('Altitude (km)', fontsize=12)
    ax.set_title('Tether Equilibrium Temperature vs Altitude\n'
                 '(solar heated, radiative + convective cooling)',
                 fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 100)

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"03-temperature_profile.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


# =============================================================================
# CONSOLE OUTPUT
# =============================================================================

def print_results(A_bottom):
    """Print summary of atmospheric loading analysis."""
    sep = "=" * 78

    print(f"\n{sep}")
    print("ATMOSPHERIC LOADING ON GROUND-ANCHORED TETHER")
    print(sep)
    print(f"  {'Reference bottom area':30} {A_bottom*1e4:>12.2f} cm^2")
    print(f"  {'Bottom diameter':30} {math.sqrt(4*A_bottom/math.pi)*1000:>12.3f} mm")
    print(f"  {'Drag coefficient (cylinder)':30} {C_D_CYLINDER:>12.1f}")
    print(f"  {'CNT operating stress':30} {cfg.SIGMA_OPERATING/1e9:>12.1f} GPa")
    print(f"  {'CNT density':30} {cfg.RHO_CNT:>12.0f} kg/m^3")

    # Key point about ground-anchored tether
    print(f"\n  NOTE: Ground-anchored tether is STATIONARY relative to Earth.")
    print(f"  Both tether and lower atmosphere co-rotate with Earth.")
    print(f"  No orbital-velocity aerodynamic heating. Wind is the only")
    print(f"  source of lateral aerodynamic loading.")

    # Table at key altitudes
    key_alts_km = [0, 1, 5, 10, 20, 50, 80, 100]

    print(f"\n{sep}")
    print("ATMOSPHERIC LOADING AT KEY ALTITUDES")
    print(sep)

    header = (f"  {'Alt':>6}  {'rho':>11}  {'Wind':>7}  {'Drag/m':>10}  "
              f"{'Diameter':>10}  {'Tension':>12}  {'Defl angle':>11}")
    units =  (f"  {'(km)':>6}  {'(kg/m^3)':>11}  {'(m/s)':>7}  {'(N/m)':>10}  "
              f"{'(mm)':>10}  {'(kN)':>12}  {'(urad)':>11}")
    print(header)
    print(units)
    print(f"  {'-'*6}  {'-'*11}  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*11}")

    for alt_km in key_alts_km:
        h = alt_km * 1000.0
        T_atm, p_atm, rho_atm = cfg.atmosphere(h)
        v_wind = cfg.wind_speed(h)
        f_drag = wind_drag_per_meter(h, A_bottom)
        d = tether_diameter(h, A_bottom)
        T_ten = tether_tension(h, A_bottom)
        theta = f_drag / T_ten if T_ten > 0 else 0.0

        print(f"  {alt_km:>6.0f}  {rho_atm:>11.4e}  {v_wind:>7.1f}  "
              f"{f_drag:>10.4e}  {d*1000:>10.4f}  {T_ten/1000:>12.3f}  "
              f"{theta*1e6:>11.4f}")

    # Find peak drag
    print(f"\n{sep}")
    print("PEAK DRAG AND MAXIMUM DEFLECTION")
    print(sep)

    n_search = 2000
    max_drag = 0.0
    max_drag_alt = 0.0
    for i in range(n_search):
        h = i * H_ATM_TOP / (n_search - 1)
        f_drag = wind_drag_per_meter(h, A_bottom)
        if f_drag > max_drag:
            max_drag = f_drag
            max_drag_alt = h

    print(f"  {'Peak drag altitude':30} {max_drag_alt/1000:>12.1f} km")
    print(f"  {'Peak drag per meter':30} {max_drag:>12.4e} N/m")
    print(f"  {'Peak drag per meter':30} {max_drag*1000:>12.4f} mN/m")

    # Get full deflection profile for max deflection
    h_arr, defl_arr, angle_arr = lateral_deflection_profile(A_bottom, 2000)
    max_defl = max(defl_arr)
    max_defl_alt = h_arr[defl_arr.index(max_defl)]

    print(f"  {'Maximum lateral deflection':30} {max_defl:>12.6f} m")
    print(f"  {'Maximum deflection':30} {max_defl*1000:>12.4f} mm")
    print(f"  {'At altitude':30} {max_defl_alt/1000:>12.1f} km")

    # Wind stagnation temperature rise
    v_jet = 60.0  # jet stream speed
    delta_T_jet = v_jet**2 / (2.0 * C_P_AIR)
    print(f"\n  {'Jet stream stagnation heating':30} {delta_T_jet:>12.2f} K (negligible)")

    # Temperature range
    print(f"\n{sep}")
    print("THERMAL ENVIRONMENT")
    print(sep)
    print(f"  {'Solar absorptivity (CNT)':30} {ALPHA_SOLAR:>12.2f}")
    print(f"  {'IR emissivity (CNT)':30} {EPSILON_IR:>12.2f}")
    print(f"  {'Solar flux':30} {SOLAR_FLUX:>12.0f} W/m^2")

    T_min = float('inf')
    T_max = float('-inf')
    T_min_alt = 0.0
    T_max_alt = 0.0

    temp_alts = []
    for i in range(201):
        h = i * H_ATM_TOP / 200
        T_eq = equilibrium_temperature(h)
        if T_eq < T_min:
            T_min = T_eq
            T_min_alt = h
        if T_eq > T_max:
            T_max = T_eq
            T_max_alt = h
        temp_alts.append((h, T_eq))

    print(f"\n  {'Minimum tether temperature':30} {T_min:>8.1f} K  ({T_min-273.15:>7.1f} deg C)"
          f"  at {T_min_alt/1000:.0f} km")
    print(f"  {'Maximum tether temperature':30} {T_max:>8.1f} K  ({T_max-273.15:>7.1f} deg C)"
          f"  at {T_max_alt/1000:.0f} km")

    # Radiative equilibrium (no convection, for reference)
    T_rad = (ALPHA_SOLAR * SOLAR_FLUX / (math.pi * EPSILON_IR * SIGMA_SB)) ** 0.25
    print(f"\n  {'Radiative equilibrium (space)':30} {T_rad:>8.1f} K  ({T_rad-273.15:>7.1f} deg C)")
    print(f"  This is the temperature in vacuum with direct solar illumination.")
    print(f"  In lower atmosphere, convective cooling brings temperature closer")
    print(f"  to ambient air temperature.")

    # Temperature at key altitudes
    print(f"\n  {'Alt (km)':>10}  {'T_atm (K)':>10}  {'T_tether (K)':>12}  {'T_tether (C)':>12}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}")
    for alt_km in key_alts_km:
        h = alt_km * 1000.0
        T_atm, _, _ = cfg.atmosphere(h)
        T_eq = equilibrium_temperature(h)
        print(f"  {alt_km:>10.0f}  {T_atm:>10.1f}  {T_eq:>12.1f}  {T_eq-273.15:>12.1f}")

    print(f"\n{sep}\n")


# =============================================================================
# CSV EXPORT
# =============================================================================

def export_csv(profile, A_bottom):
    """Export atmospheric loading profile data as CSV."""
    filename = f"atmospheric_loading_A{A_bottom:.0e}.csv"
    filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
    os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'altitude_km',
            'air_density_kg_per_m3',
            'wind_speed_m_per_s',
            'drag_per_meter_N_per_m',
            'tether_diameter_mm',
            'tether_area_m2',
            'tether_tension_kN',
            'deflection_angle_urad',
            'cumulative_deflection_m',
            'equilibrium_temperature_K',
        ])
        for i in range(len(profile['h'])):
            writer.writerow([
                f"{profile['h'][i]/1000:.4f}",
                f"{profile['rho'][i]:.6e}",
                f"{profile['wind'][i]:.2f}",
                f"{profile['drag'][i]:.6e}",
                f"{profile['diameter'][i]*1000:.6f}",
                f"{profile['area'][i]:.6e}",
                f"{profile['tension'][i]/1000:.6f}",
                f"{profile['angle'][i]*1e6:.6f}",
                f"{profile['deflection'][i]:.8f}",
                f"{profile['temperature'][i]:.2f}",
            ])

    print(f"  Saved CSV: {filename}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main entry point."""
    sys.stdout.reconfigure(encoding='utf-8')

    # Parse CLI arguments
    show_graphs = []
    A_bottom = 1e-4  # default: 1 cm^2

    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--save":
                cfg.SAVE_GRAPHS = True
            elif arg == "--show":
                cfg.SAVE_GRAPHS = False
            elif arg.startswith("--area="):
                try:
                    A_bottom = float(arg.split("=", 1)[1])
                except ValueError:
                    print(f"Error: invalid area value: {arg}")
                    return
            else:
                show_graphs.append(arg)

    # Create output directory if saving
    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    # Print results to console
    print_results(A_bottom)

    # Compute full profile for graphs and CSV
    profile = atmospheric_profile(A_bottom, n_points=1000)

    # Export CSV
    export_csv(profile, A_bottom)

    # Generate graphs
    if show_graphs:
        print("Generating graphs...")
        if "drag_loading" in show_graphs or "all" in show_graphs:
            plot_drag_loading(profile)
        if "deflection" in show_graphs or "all" in show_graphs:
            plot_deflection(profile)
        if "temperature" in show_graphs or "all" in show_graphs:
            plot_temperature(profile)
        print("Done.")


if __name__ == "__main__":
    main()
