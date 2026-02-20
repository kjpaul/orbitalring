#!/usr/bin/env python3
"""
Lunar Mass Driver Mission Analysis — 11-Part Demonstration

Self-contained analysis script demonstrating the physics and capabilities
of a lunar surface mass driver using LSM technology.

Parts:
  1. Motor Summary — P/L invariance, F/L and v_cross at different N
  2. Centrifugal Table — acceleration vs velocity, v_orbital, rail forces
  3. Launch Profiles — simulated launches for escape, Mars, Jupiter targets
  4. Moon Phase Sensitivity — how orbital phase affects Mars mission velocity
  5. 3g Crew Analysis — launch time with human g-load limit
  6. 30 km/s Equivalent — what extreme velocity looks like on the Moon
  7. HVDC Sensitivity — power limit sweep from 50 to 500 GW
  8. Sled Mass Sensitivity — effect of sled mass on acceleration and time
  9. HVDC Dominance — why N, B_sled, τ_p are irrelevant for lunar driver
 10. G-Force Table — Earth-equivalent g-forces experienced by cargo
 11. Fast Transfers — trip time vs launch velocity, spacecraft Δv needs

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# Ensure UTF-8 output on Windows
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# =============================================================================
# CONSTANTS (self-contained — no imports from other lunar_*.py modules)
# =============================================================================

# Moon
R_MOON = 1_737_400          # m
G_MOON = 1.623              # m/s²
G_0 = 9.807                 # m/s²
C_MOON = 2 * math.pi * R_MOON  # ~10,917,000 m
V_ORBITAL = math.sqrt(G_MOON * R_MOON)  # ~1,679 m/s
V_ESCAPE = math.sqrt(2.0 * G_MOON * R_MOON)  # ~2,375 m/s

# Motor calibration (at I_TARGET = 1,843 A)
V_MAX = 100_000             # V, coil voltage limit
W_COIL = 2.0                # m
F_PER_M_CAL = 1500.0        # N/m at N=10, B=0.10 (I_TARGET)
N_CAL = 10
B_CAL = 0.10

# Sled
L_SLED = 10_000             # m
M_SLED_PER_M = 400          # kg/m
M_SPACECRAFT = 5_000_000    # kg
M_SLED_HW = M_SLED_PER_M * L_SLED  # 4,000,000 kg
M_TOTAL = M_SLED_HW + M_SPACECRAFT  # 9,000,000 kg
M_PER_M_TOTAL = M_SLED_PER_M + M_SPACECRAFT / L_SLED  # 900 kg/m

# Default motor (N=50, B=0.10)
N_DEFAULT = 50
B_DEFAULT = 0.10

# Celestial mechanics
MU_SUN = 1.32712440018e20   # m³/s²
MU_EARTH = 3.986004418e14   # m³/s²
AU = 1.495978707e11         # m
A_EARTH = 1.000 * AU
A_MARS = 1.524 * AU
A_JUPITER = 5.204 * AU
A_SATURN = 9.583 * AU
A_MOON_ORBIT = 384_400_000  # m
V_MOON_ORBITAL = math.sqrt(MU_EARTH / A_MOON_ORBIT)  # ~1,018 m/s
V_ESC_EARTH_AT_MOON = math.sqrt(2.0 * MU_EARTH / A_MOON_ORBIT)  # ~1,440 m/s

DT = 0.5                    # s, simulation timestep


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def thrust_per_m(N, B_sled):
    """Thrust per metre of sled."""
    return F_PER_M_CAL * (N / N_CAL) * (B_sled / B_CAL)

def v_crossover(N, B_sled):
    """Voltage crossover velocity."""
    return V_MAX / (2.0 * N * B_sled * W_COIL)

def power_per_m(N, B_sled):
    """Power per metre (invariant)."""
    return thrust_per_m(N, B_sled) * v_crossover(N, B_sled)


def simulate_launch(v_target, m_total=M_TOTAL, l_sled=L_SLED,
                    N=N_DEFAULT, B_sled=B_DEFAULT,
                    p_hvdc_max=200e9, g_limit=None, dt=DT):
    """Simulate a lunar mass driver launch.

    Returns:
        dict with keys: t, v_final, reached, distance, laps, KE,
        peak_thrust, peak_power, peak_g, peak_rail_force,
        power_limit_v, regime_at_end
    """
    F_per_m = thrust_per_m(N, B_sled)
    F_total = F_per_m * l_sled
    v_cross = v_crossover(N, B_sled)
    P_total = F_total * v_cross
    m_per_m = m_total / l_sled

    v = 0.0
    t = 0.0
    x = 0.0
    max_time = 72000.0
    power_limit_v = None
    peak_thrust = 0.0
    peak_power = 0.0
    peak_g = 0.0
    peak_rail = 0.0
    regime = "thrust"

    while v < v_target and t < max_time:
        # Three-regime model
        if v <= v_cross or v <= 0:
            F = F_total
            P = F * max(v, 0)
            regime = "thrust"
        else:
            P = P_total
            F = P_total / v
            regime = "V-lim"

        # HVDC power limit
        if p_hvdc_max is not None and P > p_hvdc_max and v > 0:
            if power_limit_v is None:
                power_limit_v = v
            P = p_hvdc_max
            F = P / v
            regime = "P-lim"

        a = F / m_total

        # G-load limiting
        a_net_out = v**2 / R_MOON - G_MOON
        a_total_g = math.sqrt(a**2 + a_net_out**2) / G_0

        if g_limit is not None and a_total_g > g_limit:
            a_radial = abs(a_net_out)
            a_max_sq = (g_limit * G_0)**2 - a_radial**2
            if a_max_sq > 0:
                a = min(a, math.sqrt(a_max_sq))
                F = a * m_total
                P = F * v
            else:
                break
            a_total_g = math.sqrt(a**2 + a_net_out**2) / G_0

        # Rail force
        F_rail = m_per_m * max(0, v**2 / R_MOON - G_MOON)

        peak_thrust = max(peak_thrust, F)
        peak_power = max(peak_power, P)
        peak_g = max(peak_g, a_total_g)
        peak_rail = max(peak_rail, F_rail)

        v += a * dt
        x += v * dt
        t += dt

    KE = 0.5 * m_total * v**2
    return {
        't': t,
        'v_final': v,
        'reached': v >= v_target * 0.99,
        'distance': x,
        'laps': x / C_MOON,
        'KE': KE,
        'peak_thrust': peak_thrust,
        'peak_power': peak_power,
        'peak_g': peak_g,
        'peak_rail': peak_rail,
        'power_limit_v': power_limit_v,
        'regime': regime,
    }


def hohmann_v_inf(a_dep, a_arr, mu=MU_SUN):
    """v_inf at departure for Hohmann transfer."""
    a_t = (a_dep + a_arr) / 2.0
    v_dep = math.sqrt(mu * (2.0 / a_dep - 1.0 / a_t))
    v_circ = math.sqrt(mu / a_dep)
    return v_dep - v_circ

def hohmann_time(a_dep, a_arr, mu=MU_SUN):
    """Transfer time for Hohmann orbit (half period)."""
    a_t = (a_dep + a_arr) / 2.0
    return math.pi * math.sqrt(a_t**3 / mu)

def v_inf_moon_for_phase(v_inf_earth, phase_deg):
    """Required v_inf at Moon for given Earth v_inf and phase."""
    v_needed = math.sqrt(v_inf_earth**2 + V_ESC_EARTH_AT_MOON**2)
    phi = math.radians(phase_deg)
    cos_phi = math.cos(phi)
    b = 2.0 * V_MOON_ORBITAL * cos_phi
    c = V_MOON_ORBITAL**2 - v_needed**2
    disc = b**2 - 4.0 * c
    if disc < 0:
        return float('inf')
    return max(0.0, (-b + math.sqrt(disc)) / 2.0)

def launch_v_from_moon(v_inf_moon):
    """Launch velocity from Moon surface."""
    return math.sqrt(v_inf_moon**2 + V_ESCAPE**2)

def mission_launch_v(a_arr, phase_deg):
    """Full chain: target planet → launch velocity from Moon."""
    vi_earth = hohmann_v_inf(A_EARTH, a_arr)
    vi_moon = v_inf_moon_for_phase(vi_earth, phase_deg)
    return launch_v_from_moon(vi_moon)


# =============================================================================
# PART 1: MOTOR SUMMARY — P/L INVARIANCE
# =============================================================================

print("=" * 76)
print("LUNAR MASS DRIVER MISSION ANALYSIS")
print("=" * 76)

print("\n" + "-" * 76)
print("PART 1: MOTOR SUMMARY — P/L IS INVARIANT WITH N AND B_sled")
print("-" * 76)
print()
print("  Same LSM motor as orbital ring: air-core stator + SC sled coils")
print("  Operating at I_TARGET (1,843 A, 70% of I_critical)")
print()
print("  F/L = 1500 x (N/10) x (B_sled/0.10)  N/m")
print("  v_cross = 100,000 / (2 x N x B_sled x 2.0)  m/s")
print("  P/L = F/L x v_cross = 37.5 MW/m  (INVARIANT)")
print()

N_values = [10, 20, 50, 100, 200]
print(f"{'N':>6}  {'F/L (N/m)':>12}  {'v_cross (m/s)':>14}  {'P/L (MW/m)':>12}  "
      f"{'F_total (MN)':>12}  {'a_peak (g)':>10}")
print("-" * 76)

for N in N_values:
    fl = thrust_per_m(N, 0.10)
    vc = v_crossover(N, 0.10)
    pl = power_per_m(N, 0.10)
    ft = fl * L_SLED
    a_g = ft / M_TOTAL / G_0
    print(f"{N:6d}  {fl:12,.0f}  {vc:14,.0f}  {pl/1e6:12.1f}  "
          f"{ft/1e6:12.2f}  {a_g:10.3f}")

print()
print(f"  Default: N={N_DEFAULT}, F/L = {thrust_per_m(N_DEFAULT, B_DEFAULT):,.0f} N/m, "
      f"v_cross = {v_crossover(N_DEFAULT, B_DEFAULT):,.0f} m/s")
print(f"  Sled: {L_SLED/1000:.0f} km, {M_SLED_HW/1000:,.0f} t hw + "
      f"{M_SPACECRAFT/1000:,.0f} t payload = {M_TOTAL/1000:,.0f} t total")


# =============================================================================
# PART 2: CENTRIFUGAL TABLE
# =============================================================================

print("\n" + "-" * 76)
print("PART 2: CENTRIFUGAL ACCELERATION AND RAIL FORCES")
print("-" * 76)
print()
print(f"  v_orbital = sqrt(g x R) = {V_ORBITAL:,.1f} m/s  ({V_ORBITAL/1000:.2f} km/s)")
print(f"  v_escape  = sqrt(2) x v_orbital = {V_ESCAPE:,.1f} m/s  ({V_ESCAPE/1000:.2f} km/s)")
print()
print(f"  Below v_orbital: gravity > centrifugal (sled stays on track)")
print(f"  Above v_orbital: centrifugal > gravity (rails hold sled down)")
print()

v_table = [500, 1000, V_ORBITAL, V_ESCAPE, 3000, 5000, 7500, 10000, 15000, 20000]

print(f"{'v (m/s)':>10}  {'v (km/s)':>10}  {'a_cent':>10}  {'g_eff':>10}  "
      f"{'F_rail':>10}  {'F_rail':>10}")
print(f"{'':>10}  {'':>10}  {'(m/s²)':>10}  {'(m/s²)':>10}  "
      f"{'(N/m)':>10}  {'(kN/m)':>10}")
print("-" * 68)

for v in v_table:
    a_cent = v**2 / R_MOON
    g_eff = G_MOON - a_cent
    F_rail = M_PER_M_TOTAL * max(0, a_cent - G_MOON)
    marker = " <-- v_orbital" if abs(v - V_ORBITAL) < 1 else \
             " <-- v_escape" if abs(v - V_ESCAPE) < 1 else ""
    print(f"{v:10,.0f}  {v/1000:10.2f}  {a_cent:10.3f}  {g_eff:10.3f}  "
          f"{F_rail:10,.0f}  {F_rail/1000:10.1f}{marker}")

print()
print(f"  At v = 5 km/s:  rail force = "
      f"{M_PER_M_TOTAL * (5000**2/R_MOON - G_MOON)/1000:.1f} kN/m")
print(f"  At v = 10 km/s: rail force = "
      f"{M_PER_M_TOTAL * (10000**2/R_MOON - G_MOON)/1000:.1f} kN/m")


# =============================================================================
# PART 3: LAUNCH PROFILES
# =============================================================================

print("\n" + "-" * 76)
print("PART 3: LAUNCH PROFILES (N=50, 200 GW HVDC, no g-limit)")
print("-" * 76)
print()

# Mission targets
missions = [
    ("Lunar escape", V_ESCAPE),
    ("Mars (optimal phase)", mission_launch_v(A_MARS, 0)),
    ("Mars (90° phase)", mission_launch_v(A_MARS, 90)),
    ("Mars (worst phase)", mission_launch_v(A_MARS, 180)),
    ("Jupiter (optimal)", mission_launch_v(A_JUPITER, 0)),
    ("Jupiter (worst)", mission_launch_v(A_JUPITER, 180)),
    ("Saturn (optimal)", mission_launch_v(A_SATURN, 0)),
    ("10 km/s cargo", 10_000),
    ("15 km/s fast cargo", 15_000),
]

print(f"{'Mission':>24}  {'v_target':>10}  {'t_launch':>10}  {'distance':>10}  "
      f"{'laps':>6}  {'peak P':>8}  {'peak g':>7}  {'rail kN/m':>10}")
print(f"{'':>24}  {'(km/s)':>10}  {'(min)':>10}  {'(km)':>10}  "
      f"{'':>6}  {'(GW)':>8}  {'':>7}  {'':>10}")
print("-" * 100)

for name, v_target in missions:
    r = simulate_launch(v_target, p_hvdc_max=200e9)
    print(f"{name:>24}  {v_target/1000:10.2f}  {r['t']/60:10.1f}  "
          f"{r['distance']/1000:10,.0f}  {r['laps']:6.2f}  "
          f"{r['peak_power']/1e9:8.1f}  {r['peak_g']:7.2f}  "
          f"{r['peak_rail']/1000:10.1f}")


# =============================================================================
# PART 4: MOON PHASE SENSITIVITY — MARS HOHMANN
# =============================================================================

print("\n" + "-" * 76)
print("PART 4: MOON PHASE SENSITIVITY — MARS HOHMANN TRANSFER")
print("-" * 76)
print()

v_inf_earth_mars = hohmann_v_inf(A_EARTH, A_MARS)
t_transfer_mars = hohmann_time(A_EARTH, A_MARS)

print(f"  Earth-Mars Hohmann:")
print(f"    v_inf at Earth:     {v_inf_earth_mars:,.0f} m/s ({v_inf_earth_mars/1000:.2f} km/s)")
print(f"    Transfer time:      {t_transfer_mars/86400:.0f} days ({t_transfer_mars/86400/30.44:.1f} months)")
print(f"    Moon orbital vel:   {V_MOON_ORBITAL:,.0f} m/s")
print(f"    Earth escape at Moon: {V_ESC_EARTH_AT_MOON:,.0f} m/s")
print()

phases = list(range(0, 181, 15))
print(f"{'Phase':>8}  {'v_inf_moon':>12}  {'v_launch':>10}  {'t_launch':>10}  "
      f"{'peak P':>8}  {'regime':>8}")
print(f"{'(deg)':>8}  {'(m/s)':>12}  {'(km/s)':>10}  {'(min)':>10}  "
      f"{'(GW)':>8}  {'':>8}")
print("-" * 64)

for phase in phases:
    vi_moon = v_inf_moon_for_phase(v_inf_earth_mars, phase)
    v_launch = launch_v_from_moon(vi_moon)
    r = simulate_launch(v_launch, p_hvdc_max=200e9)
    print(f"{phase:8d}  {vi_moon:12,.0f}  {v_launch/1000:10.2f}  "
          f"{r['t']/60:10.1f}  {r['peak_power']/1e9:8.1f}  {r['regime']:>8}")

print()
v_opt = mission_launch_v(A_MARS, 0)
v_worst = mission_launch_v(A_MARS, 180)
print(f"  Optimal phase (0°):   {v_opt/1000:.2f} km/s launch → below v_cross")
print(f"  Worst phase (180°):   {v_worst/1000:.2f} km/s launch → "
      f"{'above' if v_worst > v_crossover(N_DEFAULT, B_DEFAULT) else 'below'} v_cross")
print(f"  v_cross at N=50:      {v_crossover(N_DEFAULT, B_DEFAULT)/1000:.1f} km/s")
print(f"  Phase effect:         {(v_worst - v_opt)/v_opt*100:.0f}% increase in "
      f"launch velocity")


# =============================================================================
# PART 5: 3g CREW ANALYSIS
# =============================================================================

print("\n" + "-" * 76)
print("PART 5: 3g CREW ANALYSIS — HUMAN-RATED LAUNCHES")
print("-" * 76)
print()
print("  Crewed missions impose a 3g physiological limit.")
print(f"  Motor peak thrust: {thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/1e6:.2f} MN")
print(f"  Unconstrained peak: "
      f"{thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.3f} g")
print()

crew_targets = [
    ("Lunar escape", V_ESCAPE),
    ("Mars (optimal)", mission_launch_v(A_MARS, 0)),
    ("Mars (worst)", mission_launch_v(A_MARS, 180)),
]

print(f"{'Mission':>20}  {'v_target':>10}  {'t (no lim)':>10}  {'t (3g)':>10}  "
      f"{'overhead':>10}")
print(f"{'':>20}  {'(km/s)':>10}  {'(min)':>10}  {'(min)':>10}  {'(%)':>10}")
print("-" * 68)

for name, v_target in crew_targets:
    r_unlim = simulate_launch(v_target, p_hvdc_max=200e9, g_limit=None)
    r_crew = simulate_launch(v_target, p_hvdc_max=200e9, g_limit=3.0)

    overhead = (r_crew['t'] - r_unlim['t']) / r_unlim['t'] * 100
    print(f"{name:>20}  {v_target/1000:10.2f}  {r_unlim['t']/60:10.1f}  "
          f"{r_crew['t']/60:10.1f}  {overhead:10.1f}")

print()
print("  Motor peak accel is only "
      f"{thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.2f}g,")
print("  well below the 3g limit. Crew missions are NOT g-limited —")
print("  the motor simply cannot push hard enough to exceed 3g.")


# =============================================================================
# PART 6: 30 km/s EQUIVALENT — EXTREME VELOCITY
# =============================================================================

print("\n" + "-" * 76)
print("PART 6: 30 km/s ON THE MOON — EXTREME VELOCITY ANALYSIS")
print("-" * 76)
print()
print("  What would it take to reach orbital ring velocities (30 km/s)")
print("  using a lunar mass driver?")
print()

r_30 = simulate_launch(30_000, p_hvdc_max=200e9)
r_30_nolim = simulate_launch(30_000, p_hvdc_max=None)

print(f"  With 200 GW HVDC limit:")
print(f"    Time:               {r_30['t']/3600:.2f} hr ({r_30['t']/60:.0f} min)")
print(f"    Distance:           {r_30['distance']/1000:,.0f} km")
print(f"    Laps of Moon:       {r_30['laps']:.1f}")
print(f"    Peak rail force:    {r_30['peak_rail']/1000:.1f} kN/m")
print(f"    Final KE:           {r_30['KE']/1e12:.1f} TJ")
print()
print(f"  Without power limit:")
print(f"    Time:               {r_30_nolim['t']/3600:.2f} hr ({r_30_nolim['t']/60:.0f} min)")
print(f"    Distance:           {r_30_nolim['distance']/1000:,.0f} km")
print(f"    Laps of Moon:       {r_30_nolim['laps']:.1f}")
print()

F_rail_30 = M_PER_M_TOTAL * (30000**2 / R_MOON - G_MOON)
print(f"  Rail force at 30 km/s: {F_rail_30/1000:.0f} kN/m ({F_rail_30/1e6:.1f} MN/m)")
print(f"  This is extreme — {F_rail_30/1000/50:.0f}x a 50 kN/m rail rating.")
print(f"  Conclusion: 30 km/s is impractical on the Moon. The lunar driver")
print(f"  is best suited for interplanetary injection (2-10 km/s).")


# =============================================================================
# PART 7: HVDC POWER SENSITIVITY
# =============================================================================

print("\n" + "-" * 76)
print("PART 7: HVDC POWER LIMIT SENSITIVITY (Mars optimal, 3.3 km/s)")
print("-" * 76)
print()

v_mars_opt = mission_launch_v(A_MARS, 0)
P_total = thrust_per_m(N_DEFAULT, B_DEFAULT) * L_SLED * v_crossover(N_DEFAULT, B_DEFAULT)

print(f"  Target: Mars optimal phase = {v_mars_opt/1000:.2f} km/s")
print(f"  Motor P_total: {P_total/1e9:.1f} GW")
print(f"  v_cross: {v_crossover(N_DEFAULT, B_DEFAULT)/1000:.1f} km/s")
print(f"  Since {v_mars_opt/1000:.2f} km/s < v_cross ({v_crossover(N_DEFAULT, B_DEFAULT)/1000:.1f} km/s),")
print(f"  the motor never enters voltage-limited mode for this mission.")
print(f"  Power at v_target: F_total x v = {thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED*v_mars_opt/1e9:.0f} GW")
print()

# Also do Jupiter for more interesting power effects
v_jup_opt = mission_launch_v(A_JUPITER, 0)
print(f"  Also showing Jupiter optimal = {v_jup_opt/1000:.2f} km/s (above v_cross)")
print()

P_limits = [50e9, 100e9, 150e9, 200e9, 300e9, 500e9, None]

# Mars sweep
print(f"  Mars (v = {v_mars_opt/1000:.2f} km/s):")
print(f"  {'P_HVDC':>12}  {'t_launch':>10}  {'limit @':>12}")
print(f"  {'(GW)':>12}  {'(min)':>10}  {'(km/s)':>12}")
print(f"  " + "-" * 40)

for P_lim in P_limits:
    r = simulate_launch(v_mars_opt, p_hvdc_max=P_lim)
    p_str = f"{P_lim/1e9:12.0f}" if P_lim is not None else f"{'unlimited':>12}"
    pv_str = f"{r['power_limit_v']/1000:12.1f}" if r['power_limit_v'] else f"{'N/A':>12}"
    print(f"  {p_str}  {r['t']/60:10.1f}  {pv_str}")

# Jupiter sweep
print()
print(f"  Jupiter (v = {v_jup_opt/1000:.2f} km/s):")
print(f"  {'P_HVDC':>12}  {'t_launch':>10}  {'limit @':>12}")
print(f"  {'(GW)':>12}  {'(min)':>10}  {'(km/s)':>12}")
print(f"  " + "-" * 40)

for P_lim in P_limits:
    r = simulate_launch(v_jup_opt, p_hvdc_max=P_lim)
    p_str = f"{P_lim/1e9:12.0f}" if P_lim is not None else f"{'unlimited':>12}"
    pv_str = f"{r['power_limit_v']/1000:12.1f}" if r['power_limit_v'] else f"{'N/A':>12}"
    print(f"  {p_str}  {r['t']/60:10.1f}  {pv_str}")


# =============================================================================
# PART 8: SLED MASS SENSITIVITY
# =============================================================================

print("\n" + "-" * 76)
print("PART 8: SLED MASS SENSITIVITY (Mars optimal, 200 GW)")
print("-" * 76)
print()
print(f"  F_total is fixed ({thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/1e6:.2f} MN). "
      f"Heavier sled = lower acceleration.")
print()

mass_per_m_sweep = [200, 300, 400, 500, 600, 800]
v_target = mission_launch_v(A_MARS, 0)

print(f"  Target: Mars optimal = {v_target/1000:.2f} km/s")
print()
print(f"{'m_sled/m':>10}  {'m_total':>10}  {'a_peak':>8}  {'t_launch':>10}  "
      f"{'peak rail':>10}")
print(f"{'(kg/m)':>10}  {'(t)':>10}  {'(g)':>8}  {'(min)':>10}  {'(kN/m)':>10}")
print("-" * 56)

for mpm in mass_per_m_sweep:
    m_hw = mpm * L_SLED
    m_tot = m_hw + M_SPACECRAFT
    r = simulate_launch(v_target, m_total=m_tot, p_hvdc_max=200e9)
    a_peak = thrust_per_m(N_DEFAULT, B_DEFAULT) * L_SLED / m_tot / G_0
    print(f"{mpm:10d}  {m_tot/1000:10,.0f}  {a_peak:8.3f}  "
          f"{r['t']/60:10.1f}  {r['peak_rail']/1000:10.1f}")

print()
print(f"  Lighter sled (200 kg/m): faster acceleration but less structural margin")
print(f"  Heavier sled (800 kg/m): more mass to accelerate, longer launch")
print(f"  Default (400 kg/m): good balance for rail-guided lunar operations")


# =============================================================================
# PART 9: HVDC DOMINANCE — WHY N, B_sled, τ_p DON'T MATTER
# =============================================================================

print("\n" + "-" * 76)
print("PART 9: HVDC DOMINANCE — WHY N, B_sled, \u03c4_p ARE IRRELEVANT")
print("-" * 76)
print()
print("  The three-regime motor model has a key invariant:")
print("    P/L = F/L × v_cross = 37.5 MW/m  (independent of N and B_sled)")
print()
print("  With P_HVDC = 200 GW and F_total = 75 MN:")
print(f"    v_HVDC_engage = P_HVDC / F_total = "
      f"{200e9 / (thrust_per_m(N_DEFAULT, B_DEFAULT) * L_SLED) / 1000:.2f} km/s")
print(f"    v_cross = {v_crossover(N_DEFAULT, B_DEFAULT) / 1000:.1f} km/s")
print()
print("  HVDC engages at 2.67 km/s — well BELOW v_cross (5.0 km/s).")
print("  The voltage limit is NEVER reached. The motor always operates in either:")
print("    - Constant thrust (v < 2.67 km/s)")
print("    - HVDC-limited constant power (v > 2.67 km/s)")
print()
print("  Proof: any N gives the same launch time for HVDC-limited missions.")
print()

N_sweep = [10, 20, 50, 100, 200]
v_target_9 = mission_launch_v(A_MARS, 0)

print(f"  Mars optimal phase (v = {v_target_9/1000:.2f} km/s), 200 GW HVDC:")
print(f"  {'N':>6}  {'v_cross':>10}  {'F_total':>10}  {'HVDC @':>10}  {'t_launch':>10}  {'P/L':>10}")
print(f"  {'':>6}  {'(km/s)':>10}  {'(MN)':>10}  {'(km/s)':>10}  {'(min)':>10}  {'(MW/m)':>10}")
print("  " + "-" * 62)

for N in N_sweep:
    fl = thrust_per_m(N, 0.10)
    vc = v_crossover(N, 0.10)
    ft = fl * L_SLED
    pl = power_per_m(N, 0.10)
    v_hvdc_engage = 200e9 / ft if ft > 0 else 0
    r = simulate_launch(v_target_9, N=N, B_sled=0.10, p_hvdc_max=200e9)
    print(f"  {N:6d}  {vc/1000:10.1f}  {ft/1e6:10.1f}  {v_hvdc_engage/1000:10.2f}  "
          f"{r['t']/60:10.1f}  {pl/1e6:10.1f}")

print()
print("  At low N (e.g. 10): HVDC engages far above v_target → constant thrust")
print("  only, and low F_total means slow acceleration. At high N (e.g. 200):")
print("  HVDC engages early → most of launch is P-limited → fast.")
print("  Increasing N beyond ~50 gives diminishing returns — the thrust phase")
print("  shortens but HVDC phase dominates. Once HVDC-limited, more thrust")
print("  can't help; you need more power.")
print()
print("  N=50 is chosen to match the orbital ring stator coils —")
print("  reuse the same HTS coil tooling for both systems.")


# =============================================================================
# PART 10: G-FORCE TABLE — EARTH-EQUIVALENT FORCES ON CARGO
# =============================================================================

print("\n" + "-" * 76)
print("PART 10: G-FORCE TABLE — EARTH-EQUIVALENT FORCES ON CARGO")
print("-" * 76)
print()
print("  Cargo experiences three accelerations simultaneously:")
print("    1. Tangential (motor thrust): forward along track")
print("    2. Centrifugal (v\u00b2/R):       outward from lunar surface")
print("    3. Gravity (g_moon):           inward toward Moon center")
print()
print("  Total g-load = sqrt(a_tangential\u00b2 + a_net_radial\u00b2) / g_Earth")
print(f"  where a_net_radial = v\u00b2/R - g_moon")
print()
print(f"  At rest:       a_tang = {thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL:.2f} m/s\u00b2 = "
      f"{thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.3f} g")
print(f"  Motor peak:    {thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.3f} g "
      f"(well below human 3g limit)")
print()

# Build the table
g_table_velocities = [0, 500, 1000, V_ORBITAL, V_ESCAPE, 3000, 4000, 5000,
                      7500, 10000, 15000, 20000]

F_per_m_default = thrust_per_m(N_DEFAULT, B_DEFAULT)
F_total_default = F_per_m_default * L_SLED
P_hvdc = 200e9

print(f"{'v':>10}  {'a_tang':>8}  {'a_cent':>8}  {'g_eff':>8}  {'a_net_rad':>9}  "
      f"{'g_total':>8}  {'regime':>8}")
print(f"{'(km/s)':>10}  {'(g)':>8}  {'(m/s\u00b2)':>8}  {'(m/s\u00b2)':>8}  {'(m/s\u00b2)':>9}  "
      f"{'(g)':>8}  {'':>8}")
print("-" * 72)

for v in g_table_velocities:
    # Determine thrust at this velocity (three-regime + HVDC)
    v_cross_default = v_crossover(N_DEFAULT, B_DEFAULT)
    P_total_default = F_total_default * v_cross_default

    if v <= 0:
        F = F_total_default
        P = 0
        regime = "thrust"
    elif v <= v_cross_default:
        F = F_total_default
        P = F * v
        regime = "thrust"
    else:
        P = P_total_default
        F = P / v
        regime = "V-lim"

    if P_hvdc is not None and P > P_hvdc and v > 0:
        P = P_hvdc
        F = P / v
        regime = "P-lim"

    a_tang = F / M_TOTAL
    a_cent = v**2 / R_MOON
    g_eff = G_MOON - a_cent
    a_net_rad = a_cent - G_MOON  # positive = outward
    g_total = math.sqrt(a_tang**2 + a_net_rad**2) / G_0

    marker = ""
    if abs(v - V_ORBITAL) < 1:
        marker = " <-- v_orb"
    elif abs(v - V_ESCAPE) < 1:
        marker = " <-- v_esc"

    print(f"{v/1000:10.2f}  {a_tang/G_0:8.3f}  {a_cent:8.3f}  {g_eff:8.3f}  "
          f"{a_net_rad:9.3f}  {g_total:8.3f}  {regime:>8}{marker}")

print()
print("  Key observations:")
print(f"  - Motor thrust gives only {F_total_default/M_TOTAL/G_0:.3f}g at full power")
print(f"  - At 5 km/s: centrifugal = {5000**2/R_MOON:.1f} m/s\u00b2 "
      f"({5000**2/R_MOON/G_0:.2f}g) but net radial = "
      f"{5000**2/R_MOON - G_MOON:.1f} m/s\u00b2 ({(5000**2/R_MOON - G_MOON)/G_0:.2f}g)")
print(f"  - At 10 km/s: net radial dominates at "
      f"{(10000**2/R_MOON - G_MOON)/G_0:.1f}g — purely structural load")
print(f"  - Cargo never experiences more than ~{F_total_default/M_TOTAL/G_0:.1f}g from motor")
print(f"  - Rail structural loads at high v are the real engineering challenge")


# =============================================================================
# PART 11: FAST TRANSFERS — BEYOND HOHMANN
# =============================================================================

print("\n" + "-" * 76)
print("PART 11: FAST TRANSFERS — TRIP TIME VS LAUNCH VELOCITY")
print("-" * 76)
print()
print("  Hohmann transfers minimize \u0394v but take months/years.")
print("  Faster launches inject into steeper transfer orbits, trading")
print("  mass driver energy for shorter trip times.")
print()
print("  The spacecraft gets a velocity boost from the mass driver, then")
print("  uses its own propulsion for course corrections and arrival braking.")
print("  v_inf_arrival is the speed the spacecraft must shed (or deflect)")
print("  at the destination.")
print()


def fast_transfer_time(v_inf_earth, a_dep, a_arr, mu=MU_SUN):
    """Calculate transfer time for arbitrary v_inf (elliptic or hyperbolic)."""
    r1 = a_dep
    r2 = a_arr

    v_circ = math.sqrt(mu / r1)
    v_dep = v_circ + v_inf_earth

    eps = 0.5 * v_dep**2 - mu / r1
    h = r1 * v_dep

    if abs(eps) < 1e-6:
        a_t = 1e20
    else:
        a_t = -mu / (2.0 * eps)

    p = h**2 / mu
    if a_t > 0:
        e = math.sqrt(max(0.0, 1.0 - p / a_t))
        r_apo = a_t * (1.0 + e)
        if r2 > r_apo * 1.001:
            return float('inf'), 0.0, 0.0
    else:
        e = math.sqrt(1.0 + p / abs(a_t))

    # Speed at r2
    v_at_r2 = math.sqrt(max(0.0, 2.0 * (eps + mu / r2)))
    v_circ_arr = math.sqrt(mu / r2)
    v_inf_arr = abs(v_at_r2 - v_circ_arr)

    # True anomaly at r2 (departure at perihelion, nu=0)
    cos_nu2 = (p / r2 - 1.0) / e if e > 1e-10 else 0.0
    cos_nu2 = max(-1.0, min(1.0, cos_nu2))
    nu2 = math.acos(cos_nu2)

    if e < 1.0:
        E2 = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(nu2 / 2.0),
                                math.sqrt(1.0 + e) * math.cos(nu2 / 2.0))
        M2 = E2 - e * math.sin(E2)
        n = math.sqrt(mu / a_t**3)
        dt = M2 / n
        if dt < 0:
            dt += 2.0 * math.pi / n
    else:
        a_h = abs(a_t)
        tan_half = math.tan(nu2 / 2.0)
        tanh_H = math.sqrt((e - 1.0) / (e + 1.0)) * tan_half
        tanh_H = max(-0.999999, min(0.999999, tanh_H))
        H2 = 2.0 * math.atanh(tanh_H)
        M_h = e * math.sinh(H2) - H2
        n_h = math.sqrt(mu / a_h**3)
        dt = M_h / n_h

    return abs(dt), v_inf_arr, e


# --- Mars fast transfers ---
print("  MARS TRANSFERS (Earth → Mars, a = 1.524 AU)")
print()

v_inf_hohmann_mars = hohmann_v_inf(A_EARTH, A_MARS)
t_hohmann_mars = hohmann_time(A_EARTH, A_MARS)

print(f"  Hohmann baseline: v_inf = {v_inf_hohmann_mars/1000:.2f} km/s, "
      f"trip = {t_hohmann_mars/86400:.0f} days ({t_hohmann_mars/86400/30.44:.1f} months)")
print()

# Sweep v_inf from Hohmann up to 20 km/s
v_inf_sweep = [v_inf_hohmann_mars]
for vi in [3000, 5000, 7500, 10000, 15000, 20000, 30000]:
    if vi > v_inf_hohmann_mars:
        v_inf_sweep.append(float(vi))

print(f"{'v_inf_Earth':>12}  {'trip':>8}  {'trip':>8}  {'v_inf_arr':>10}  "
      f"{'v_launch':>10}  {'t_launch':>10}  {'e':>6}")
print(f"{'(km/s)':>12}  {'(days)':>8}  {'(months)':>8}  {'(km/s)':>10}  "
      f"{'(km/s)':>10}  {'(min)':>10}  {'':>6}")
print("-" * 78)

for vi in v_inf_sweep:
    dt, v_arr, e_t = fast_transfer_time(vi, A_EARTH, A_MARS)
    # Launch velocity from Moon (optimal phase)
    vi_moon = v_inf_moon_for_phase(vi, 0)
    v_launch = launch_v_from_moon(vi_moon)
    r = simulate_launch(v_launch, p_hvdc_max=200e9)

    if dt < float('inf'):
        print(f"{vi/1000:12.2f}  {dt/86400:8.0f}  {dt/86400/30.44:8.1f}  "
              f"{v_arr/1000:10.2f}  {v_launch/1000:10.2f}  "
              f"{r['t']/60:10.1f}  {e_t:6.3f}")
    else:
        print(f"{vi/1000:12.2f}  {'orbit too small':>20}")

print()
print(f"  At v_inf = 5 km/s:  trip drops to ~4 months, spacecraft needs "
      f"~{fast_transfer_time(5000, A_EARTH, A_MARS)[1]/1000:.1f} km/s braking")
print(f"  At v_inf = 10 km/s: trip drops to ~3 months, hyperbolic approach")
print(f"  At v_inf = 20 km/s: trip drops to ~2 months — fast cargo express")

# --- Jupiter fast transfers ---
print()
print("  JUPITER TRANSFERS (Earth → Jupiter, a = 5.204 AU)")
print()

v_inf_hohmann_jup = hohmann_v_inf(A_EARTH, A_JUPITER)
t_hohmann_jup = hohmann_time(A_EARTH, A_JUPITER)

print(f"  Hohmann baseline: v_inf = {v_inf_hohmann_jup/1000:.2f} km/s, "
      f"trip = {t_hohmann_jup/86400:.0f} days ({t_hohmann_jup/86400/365.25:.1f} years)")
print()

v_inf_sweep_j = [v_inf_hohmann_jup]
for vi in [10000, 15000, 20000, 30000, 50000]:
    if vi > v_inf_hohmann_jup:
        v_inf_sweep_j.append(float(vi))

print(f"{'v_inf_Earth':>12}  {'trip':>8}  {'trip':>8}  {'v_inf_arr':>10}  "
      f"{'v_launch':>10}  {'t_launch':>10}  {'e':>6}")
print(f"{'(km/s)':>12}  {'(days)':>8}  {'(years)':>8}  {'(km/s)':>10}  "
      f"{'(km/s)':>10}  {'(min)':>10}  {'':>6}")
print("-" * 78)

for vi in v_inf_sweep_j:
    dt, v_arr, e_t = fast_transfer_time(vi, A_EARTH, A_JUPITER)
    vi_moon = v_inf_moon_for_phase(vi, 0)
    v_launch = launch_v_from_moon(vi_moon)
    r = simulate_launch(v_launch, p_hvdc_max=200e9)

    if dt < float('inf'):
        print(f"{vi/1000:12.2f}  {dt/86400:8.0f}  {dt/86400/365.25:8.1f}  "
              f"{v_arr/1000:10.2f}  {v_launch/1000:10.2f}  "
              f"{r['t']/60:10.1f}  {e_t:6.3f}")
    else:
        print(f"{vi/1000:12.2f}  {'orbit too small':>20}")

print()
print("  Hohmann to Jupiter takes ~2.7 years. With mass driver boost:")
print("  - v_inf = 15 km/s:  ~1.5 years, spacecraft brakes ~10 km/s at arrival")
print("  - v_inf = 30 km/s:  ~1 year, but arrival v_inf is extreme")
print()
print("  Key insight: the mass driver provides FREE \u0394v from solar power.")
print("  Spacecraft propellant budget is only for course correction + braking,")
print("  not for the initial heliocentric injection.")


# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 76)
print("SUMMARY")
print("=" * 76)
print()
print(f"  Lunar mass driver: {L_SLED/1000:.0f} km sled, {M_SPACECRAFT/1000:,.0f} t payload, "
      f"N={N_DEFAULT}")
print(f"  P/L = 37.5 MW/m (invariant), F/L = 7,500 N/m, v_cross = 5,000 m/s")
print(f"  Track: {C_MOON/1000:,.0f} km (full lunar circumference)")
print()
print(f"  Key velocities:")
print(f"    v_orbital:    {V_ORBITAL/1000:.2f} km/s (centrifugal = gravity)")
print(f"    v_escape:     {V_ESCAPE/1000:.2f} km/s")
print(f"    Mars (opt):   {mission_launch_v(A_MARS, 0)/1000:.2f} km/s")
print(f"    Mars (worst): {mission_launch_v(A_MARS, 180)/1000:.2f} km/s")
print(f"    Jupiter:      {mission_launch_v(A_JUPITER, 0)/1000:.2f} — "
      f"{mission_launch_v(A_JUPITER, 180)/1000:.2f} km/s")
print()
print(f"  Moon phase matters: Mars launch velocity varies by "
      f"{(mission_launch_v(A_MARS, 180) - mission_launch_v(A_MARS, 0))/1000:.1f} km/s "
      f"depending on Moon's orbital phase.")
print()
print(f"  The 3g crew limit is never reached — motor peak is only "
      f"{thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.2f}g.")
print()
print(f"  HVDC power limit (200 GW) dominates performance — motor parameters")
print(f"  (N, B_sled, \u03c4_p) are irrelevant once HVDC-limited.")
print()
print(f"  Cargo g-forces: motor peak is {thrust_per_m(N_DEFAULT, B_DEFAULT)*L_SLED/M_TOTAL/G_0:.2f}g "
      f"— no g-load concern for any cargo.")
print(f"  Rail structural loads dominate at high velocity.")
print()
print(f"  Fast transfers: mass driver provides free \u0394v from solar power.")
v_inf_5k_dt, v_inf_5k_arr, _ = fast_transfer_time(5000, A_EARTH, A_MARS)
print(f"  Mars at v_inf=5 km/s: {v_inf_5k_dt/86400:.0f} day trip, "
      f"spacecraft brakes {v_inf_5k_arr/1000:.1f} km/s")
print(f"  Best suited for 2-10 km/s launches (interplanetary injection).")
print(f"  Above ~15 km/s, rail forces become extreme. 30 km/s is impractical.")
print()
print("=" * 76)
