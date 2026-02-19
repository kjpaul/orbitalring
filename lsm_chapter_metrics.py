#!/usr/bin/env python3
"""Compute chapter metrics from baseline LSM simulation results.

Baseline: N=50 stator turns, tau_p=130 m, 525 kg/m sled, 10 km, 5000t payload.
"""
import math

# =============================================================================
# BASELINE SIMULATION RESULTS (N=50, tau_p=130, 525 kg/m, 10 km sled)
# =============================================================================

V_FINAL = 30_000.0          # m/s
T_LAUNCH = 210.8 * 60       # s (210.8 minutes = 3.51 hr)
DIST = 246_682.2e3           # m
L_RING = 41_646_000          # m
N_LAPS = DIST / L_RING

M_SLED_PER_M = 525          # kg/m
L_SLED = 10_000              # m
M_SLED_HW = M_SLED_PER_M * L_SLED  # 5,250,000 kg
M_SPACECRAFT = 5_000_000    # kg
M_TOTAL = M_SLED_HW + M_SPACECRAFT  # 10,250,000 kg

F_PEAK = 74.96e6            # N — thrust-limited at delta=60 deg, I_PEAK=1896 A
P_PEAK = 374.82e9           # W (from sim)
G_0 = 9.807

# Peak acceleration (thrust-limited at delta_max=60 deg, I_peak)
A_NOM = F_PEAK / M_TOTAL    # 7.313 m/s^2 = 0.746 g
G_NOM = A_NOM / G_0

# Voltage crossover: EMF = V_max at v = V_max / (2*N*B_sled*w_coil)
# = 100,000 / (2*50*0.10*2.0) = 5,000 m/s
V_CROSSOVER = 5_000         # m/s

# HTS hysteresis (from simulation at N=50)
P_HTS_PEAK = 82.58e3        # W (82.58 kW)
E_HTS_TOTAL = 6.791e8       # J (188.6 kWh)

print("=" * 76)
print("LSM CHAPTER METRICS (N=50, tau_p=130, 525 kg/m sled, 10 km, 5000t)")
print("=" * 76)

# =============================================================================
# 1. LAUNCH PROFILE
# =============================================================================

print("\n1. LAUNCH PROFILE")
print(f"   Final velocity:          {V_FINAL/1000:.1f} km/s")
print(f"   Total launch time:       {T_LAUNCH/60:.1f} min  ({T_LAUNCH/3600:.2f} hr)")
print(f"   Distance traveled:       {DIST/1000:.1f} km")
print(f"   Ring circumferences:     {N_LAPS:.1f}")
print(f"   Peak thrust:             {F_PEAK/1e6:.2f} MN  (delta_max=60 deg, I_peak=1896 A)")
print(f"   Tangential acceleration: {A_NOM:.3f} m/s^2  ({G_NOM:.3f} g)")
print(f"   Peak power:              {P_PEAK/1e9:.2f} GW")
print(f"   Voltage crossover:       {V_CROSSOVER/1000:.0f} km/s")
print(f"   Total KE:                {0.5*M_TOTAL*V_FINAL**2/1e12:.1f} TJ")
print(f"   Payload KE:              {0.5*M_SPACECRAFT*V_FINAL**2/1e12:.1f} TJ")
print(f"   Energy efficiency:       {M_SPACECRAFT/M_TOTAL*100:.1f}%")
print(f"   KE wasted on sled:       {0.5*M_SLED_HW*V_FINAL**2/1e12:.1f} TJ")

# =============================================================================
# 2. PER 500 m STATOR SITE
# =============================================================================

SITE_LENGTH = 500            # m
N_SITES_ACTIVE = L_SLED / SITE_LENGTH  # 20 sites active simultaneously
F_PER_SITE = F_PEAK / N_SITES_ACTIVE
P_PER_SITE = P_PEAK / N_SITES_ACTIVE
P_HTS_PER_SITE = P_HTS_PEAK / N_SITES_ACTIVE
T_RECOVERY = (L_RING - L_SLED) / V_FINAL

print("\n2. PER 500 m STATOR SITE")
print(f"   Active sites (simultaneous): {N_SITES_ACTIVE:.0f}")
print(f"   Peak thrust per site:        {F_PER_SITE/1e6:.3f} MN  ({F_PER_SITE/1e3:.0f} kN)")
print(f"   Peak power per site:         {P_PER_SITE/1e9:.2f} GW")
print(f"   Peak HTS hysteresis/site:    {P_HTS_PER_SITE/1000:.2f} kW")
print(f"   Recovery time (final v):     {T_RECOVERY:.0f} s  ({T_RECOVERY/60:.1f} min)")
print(f"   Total passes per site:       {N_LAPS:.1f}")
print(f"   Exposure per pass (final v): {L_SLED/V_FINAL:.2f} s")

# =============================================================================
# 3. HVDC POWER GRID
# =============================================================================

HVDC_CAPACITY = 159e9        # W (159 GW from Chapter 6)

# Power = F * v. Grid capacity exceeded when F * v > 159 GW
# In constant-thrust phase (v < v_cross): P = F_peak * v
V_GRID_LIMIT = HVDC_CAPACITY / F_PEAK

print("\n3. HVDC POWER GRID")
print(f"   Peak power to active zone:   {P_PEAK/1e9:.2f} GW")
print(f"   Peak power per 500m site:    {P_PER_SITE/1e9:.2f} GW")
print(f"   HVDC grid capacity (Ch.6):   {HVDC_CAPACITY/1e9:.0f} GW")
print(f"   Oversubscription ratio:      {P_PEAK/HVDC_CAPACITY:.2f}x")
print(f"   Power deficit at peak:       {(P_PEAK - HVDC_CAPACITY)/1e9:.1f} GW")
print(f"   Grid capacity exceeded at:   v = {V_GRID_LIMIT/1000:.1f} km/s")
print(f"   Time at grid limit:          ~{V_GRID_LIMIT/A_NOM/60:.0f} min into launch")

# =============================================================================
# 4. SLED CONSUMABLES (LN2 for sled cryostat)
# =============================================================================

# Sled SC coils are DC — no hysteresis. Heat load is parasitic thermal leaks.
Q_MLI_PER_M2 = 1.0          # W/m^2
A_CRYO_PER_M = 2.0          # m^2/m
Q_PARASITIC_PER_M = Q_MLI_PER_M2 * A_CRYO_PER_M  # W/m
Q_PARASITIC_TOTAL = Q_PARASITIC_PER_M * L_SLED
Q_SUPPORTS_PER_M = 0.5      # W/m
Q_SUPPORTS_TOTAL = Q_SUPPORTS_PER_M * L_SLED
Q_SLED_CRYO_TOTAL = Q_PARASITIC_TOTAL + Q_SUPPORTS_TOTAL

# LN2 absorption: sensible (70->77.4 K) + latent = 214 kJ/kg
LN2_CAPACITY = 2040 * 7.4 + 199_000  # J/kg

E_CRYO_SLED = Q_SLED_CRYO_TOTAL * T_LAUNCH
M_LN2_BOILED = E_CRYO_SLED / LN2_CAPACITY

# LN2 reservoir (from mass budget: 5 kg/m)
M_LN2_RESERVOIR = 5.0 * L_SLED  # kg
V_LN2_RESERVOIR = M_LN2_RESERVOIR / 807  # m^3

print("\n4. SLED CRYOSTAT CONSUMABLES")
print(f"   MLI heat leak:               {Q_PARASITIC_PER_M:.1f} W/m  ({Q_PARASITIC_TOTAL/1000:.1f} kW total)")
print(f"   Support conduction:          {Q_SUPPORTS_PER_M:.1f} W/m  ({Q_SUPPORTS_TOTAL/1000:.1f} kW total)")
print(f"   Total parasitic load:        {Q_SLED_CRYO_TOTAL/1000:.1f} kW")
print(f"   Launch duration:             {T_LAUNCH/3600:.2f} hr")
print(f"   Total heat over launch:      {E_CRYO_SLED/1e6:.1f} MJ  ({E_CRYO_SLED/3.6e6:.1f} kWh)")
print(f"   LN2 boiled off:             {M_LN2_BOILED:.0f} kg  ({M_LN2_BOILED/1000:.1f} t)")
print(f"   LN2 reservoir (5 kg/m):      {M_LN2_RESERVOIR/1000:.0f} t  ({V_LN2_RESERVOIR:.0f} m^3)")
print(f"   Reservoir margin:            {M_LN2_RESERVOIR/M_LN2_BOILED:.0f}x boil-off")

# =============================================================================
# 5. PER-METER THRUST AND STRUCTURAL LOADING
# =============================================================================

# At peak thrust (74.96 MN, before voltage crossover)
f_per_m = F_PEAK / L_SLED
a_sled = F_PEAK / M_TOTAL
f_inertia_per_m = M_SLED_PER_M * a_sled
f_surplus_per_m = f_per_m - f_inertia_per_m

T_frame_peak = f_surplus_per_m * L_SLED
F_spacecraft_check = M_SPACECRAFT * a_sled

print("\n5. PER-METER THRUST AND STRUCTURAL LOADING (at peak thrust)")
print(f"   Peak thrust (total):         {F_PEAK/1e6:.2f} MN")
print(f"   Sled acceleration:           {a_sled:.3f} m/s^2  ({a_sled/G_0:.3f} g)")
print(f"   Thrust per metre of sled:    {f_per_m:.1f} N/m  ({f_per_m/1000:.3f} kN/m)")
print(f"   Inertia per metre (sled hw): {f_inertia_per_m:.1f} N/m")
print(f"   Net surplus per metre:       {f_surplus_per_m:.1f} N/m  ({f_surplus_per_m/1000:.3f} kN/m)")
print(f"   Surplus fraction:            {f_surplus_per_m/f_per_m*100:.1f}%  (= payload mass fraction)")
print(f"   Peak frame tension:          {T_frame_peak/1e6:.2f} MN  (at spacecraft attachment)")
print(f"   (check: M_sc x a =           {F_spacecraft_check/1e6:.2f} MN)")

# After voltage crossover
F_FINAL = 12.52e6  # from sim at v~29.9 km/s (last report line)
a_final = F_FINAL / M_TOTAL
f_per_m_final = F_FINAL / L_SLED
f_inertia_final = M_SLED_PER_M * a_final
f_surplus_final = f_per_m_final - f_inertia_final
T_frame_final = f_surplus_final * L_SLED

print(f"\n   At end of launch (~30 km/s, voltage-limited):")
print(f"   Thrust (total):              {F_FINAL/1e6:.2f} MN")
print(f"   Sled acceleration:           {a_final:.3f} m/s^2  ({a_final/G_0:.3f} g)")
print(f"   Thrust per metre:            {f_per_m_final:.1f} N/m")
print(f"   Net surplus per metre:       {f_surplus_final:.1f} N/m")
print(f"   Peak frame tension:          {T_frame_final/1e6:.2f} MN")

print(f"\n{'='*76}")
