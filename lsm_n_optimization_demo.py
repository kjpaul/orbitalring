#!/usr/bin/env python3
"""
LSM Stator Turn Count Optimization — Why N=50 Fixed is Optimal

This script demonstrates the analysis that led to choosing N=50 fixed stator
turns for the LSM mass driver, and why switchable/reconfigurable windings were
rejected. It is self-contained (no imports from other lsm_*.py modules) so
readers can run it standalone.

Key findings:
  1. P/L = 41.6 MW/m is invariant with both N and B_sled (at I_PEAK)
  2. N=50 gives only ~2.3% time overhead vs theoretical minimum
  3. Cargo missions always want max N -> switching never activates
  4. Crewed missions are g-limited -> N is irrelevant
  5. At 666 GW HVDC limit, heavy cargo (416 GW) is not power-limited
  6. Conclusion: Fixed N=50, simple and serviceable

Reference: "Orbital Ring Engineering" by Paul G de Jong, Chapter 7
"""

import math
import numpy as np

# Ensure UTF-8 output on Windows
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# =============================================================================
# CONSTANTS
# =============================================================================

G_LOCAL = 9.073             # m/s^2, gravitational acceleration at 250 km
G_0 = 9.807                # m/s^2, standard gravity at sea level
R_ORBIT = 6_628_000        # m, orbital radius at 250 km
V_MAX = 100_000            # V, coil insulation limit
W_COIL = 2.0               # m, coil width (tangential extent)
V_START = 483.0             # m/s, casing velocity (sled starts here)
MU0 = 4 * math.pi * 1e-7   # H/m, permeability of free space

# Thrust calibration at I_PEAK = 2,107 A (5-layer HTS, 80% of Ic):
# F/L = 1666 N/m at N=10, B_sled=0.10
F_PER_M_BASE = 1666.0      # N/m baseline thrust per metre of sled
N_BASE = 10                 # reference stator turns for calibration
B_BASE = 0.10               # reference sled field for calibration

# Note: at I_TARGET = 1,843 A (70% of Ic), the baseline would be 1,458 N/m,
# giving P/L = 36.4 MW/m. The optimization conclusions are the same either way
# since they depend on the ratio, not the absolute value. We use I_PEAK here
# because the heavy cargo sled (the reference case) is always thrust-limited.

# Cargo reference case
V_TARGET_CARGO = 30_000.0   # m/s
L_SLED = 10_000.0           # m
M_TOTAL_CARGO = 10_250_000  # kg (5,250 t sled + 5,000 t payload)
G_LIMIT_CARGO = 13.0        # g structural limit

# Crewed reference case
V_TARGET_CREW = 15_000.0    # m/s
M_TOTAL_CREW = 10_250_000   # kg (same sled)
G_LIMIT_CREW = 3.0          # g human physiological limit

DT = 0.5                    # s, simulation timestep


def thrust_per_m(N, B_sled):
    """Thrust per metre of sled (linear scaling from calibration point)."""
    return F_PER_M_BASE * (N / N_BASE) * (B_sled / B_BASE)


def v_crossover(N, B_sled):
    """Voltage crossover velocity: EMF = V_max."""
    return V_MAX / (2.0 * N * B_sled * W_COIL)


def power_per_m(N, B_sled):
    """Power per metre at crossover (F/L * v_cross)."""
    return thrust_per_m(N, B_sled) * v_crossover(N, B_sled)


def simulate_launch(N, B_sled_nom, m_total, v_target, g_limit,
                    p_hvdc_max=None, dt=DT):
    """Simulate a launch with given parameters.

    Args:
        N: Stator turns
        B_sled_nom: Nominal sled field (T)
        m_total: Total mass (kg)
        v_target: Target velocity (m/s)
        g_limit: G-load limit (g)
        p_hvdc_max: HVDC power limit (W), or None for unlimited
        dt: Timestep (s)

    Returns:
        (t_total, v_final, reached_target, power_limit_v)
        where power_limit_v is the velocity where HVDC limit engaged, or None.
    """
    F_per_m = thrust_per_m(N, B_sled_nom)
    F_total = F_per_m * L_SLED
    v_cross = v_crossover(N, B_sled_nom)
    P_total = F_total * v_cross  # constant power above v_cross

    v = V_START
    t = 0.0
    max_time = 72000.0  # 20 hours safety
    power_limit_v = None

    while v < v_target and t < max_time:
        # Determine thrust
        if v <= v_cross:
            F = F_total  # constant thrust phase
        else:
            F = P_total / v  # constant power phase (voltage-limited)

        # Tangential acceleration
        a_tang = F / m_total

        # G-load limiting
        a_cent = v**2 / R_ORBIT - G_LOCAL  # net radial
        a_total_g = math.sqrt(a_tang**2 + a_cent**2) / G_0

        if a_total_g > g_limit:
            a_radial = abs(a_cent)
            a_max_tang_sq = (g_limit * G_0)**2 - a_radial**2
            if a_max_tang_sq > 0:
                a_tang = min(a_tang, math.sqrt(a_max_tang_sq))
            else:
                break  # radial alone exceeds g-limit

        # HVDC power limit
        P_mech = F * v if F == F_total else P_total
        if a_tang < F / m_total:
            P_mech = a_tang * m_total * v  # g-limited power
        if p_hvdc_max is not None and P_mech > p_hvdc_max and v > 0:
            if power_limit_v is None:
                power_limit_v = v
            F_capped = p_hvdc_max / v
            a_tang = min(a_tang, F_capped / m_total)

        v += a_tang * dt
        t += dt

    return t, v, (v >= v_target * 0.99), power_limit_v


# =============================================================================
# PART 1: POWER INVARIANCE DEMONSTRATION
# =============================================================================

print("=" * 76)
print("LSM STATOR TURN COUNT OPTIMIZATION")
print("Why N=50 fixed turns is optimal; why switchable windings are useless")
print("=" * 76)

print("\n" + "-" * 76)
print("PART 1: POWER PER METRE IS INVARIANT WITH N AND B_sled")
print("-" * 76)
print()
print("  EMF_peak = 2 * N * v * B_sled * w_coil")
print("  v_cross  = V_max / (2 * N * B_sled * w_coil)")
print("  F/L      = 1666 * (N/10) * (B_sled/0.10)  N/m  [at I_PEAK]")
print("  P/L      = F/L * v_cross  ==>  the N and B_sled terms CANCEL")
print()

N_values = [10, 20, 50, 100, 200, 500]
B_values = [0.10, 0.20]

header = f"{'N':>6}  {'B_sled':>8}"
for B in B_values:
    header += f"  {'F/L (N/m)':>12}  {'v_cross':>10}  {'P/L (MW/m)':>12}"
print(header)
print("-" * len(header))

for N in N_values:
    row = f"{N:6d}  "
    for i, B in enumerate(B_values):
        fl = thrust_per_m(N, B)
        vc = v_crossover(N, B)
        pl = power_per_m(N, B)
        if i == 0:
            row += f"{B:8.2f}  {fl:12,.0f}  {vc:10,.0f}  {pl/1e6:12.1f}"
        else:
            row += f"  {fl:12,.0f}  {vc:10,.0f}  {pl/1e6:12.1f}"
    print(row)

print()
print("  Result: P/L = 41.6 MW/m in EVERY case (at I_PEAK).")
print("  N and B_sled only trade thrust for crossover velocity;")
print("  the power envelope is fixed.")

# =============================================================================
# PART 2: LAUNCH TIME VS N (CARGO, NO G-LIMIT IN PRACTICE)
# =============================================================================

print("\n" + "-" * 76)
print("PART 2: LAUNCH TIME VS N (Cargo: 30 km/s, 10 km sled, 10,250 t)")
print("-" * 76)
print()

# Theoretical minimum time: all energy delivered at constant power
P_total_fixed = power_per_m(10, 0.10) * L_SLED  # = 41.6 MW/m * 10 km = 416 GW
# Time at constant power from v_start to v_target:
# P = m*v*dv/dt => dt = m*v*dv/P => t = m/(2P) * (v_f^2 - v_s^2)
t_min = M_TOTAL_CARGO / (2 * P_total_fixed) * (V_TARGET_CARGO**2 - V_START**2)

print(f"  Constant power:     {P_total_fixed/1e9:.1f} GW (invariant)")
print(f"  Theoretical t_min:  {t_min/60:.1f} min ({t_min/3600:.2f} hr)")
print(f"  (limit as v_cross -> 0, i.e. N -> infinity)")
print()

N_sweep = [10, 20, 30, 50, 75, 100, 150, 200, 250]

print(f"{'N':>6}  {'v_cross':>10}  {'F_total':>10}  {'a_peak':>8}  "
      f"{'t_launch':>10}  {'t_min':>8}  {'overhead':>10}")
print(f"{'':>6}  {'(km/s)':>10}  {'(MN)':>10}  {'(g)':>8}  "
      f"{'(hr)':>10}  {'(hr)':>8}  {'(%)':>10}")
print("-" * 76)

for N in N_sweep:
    B = 0.10
    fl = thrust_per_m(N, B)
    F_total = fl * L_SLED
    vc = v_crossover(N, B)
    a_peak = F_total / M_TOTAL_CARGO / G_0

    t_sim, v_sim, ok, _ = simulate_launch(N, B, M_TOTAL_CARGO,
                                           V_TARGET_CARGO, G_LIMIT_CARGO)

    overhead = (t_sim - t_min) / t_min * 100

    print(f"{N:6d}  {vc/1000:10.1f}  {F_total/1e6:10.2f}  {a_peak:8.3f}  "
          f"{t_sim/3600:10.2f}  {t_min/3600:8.2f}  {overhead:10.1f}")

print()
print("  N=50: v_cross = 5.0 km/s, overhead ~ 2.3%")
print("  Beyond N=100: diminishing returns (< 0.6% improvement per doubling)")
print("  N=50 is the sweet spot: simple coils, low overhead, easy to service.")

# =============================================================================
# PART 3: WHY SWITCHABLE WINDINGS DON'T HELP — CARGO CASE
# =============================================================================

print("\n" + "-" * 76)
print("PART 3: WHY SWITCHABLE WINDINGS DON'T HELP -- CARGO")
print("-" * 76)
print()
print("  For cargo launches, the objective is minimum launch time.")
print("  More turns = more thrust = shorter constant-thrust phase.")
print("  A switchable system with stages [N=10, N=25, N=50] would:")
print()

stages_sw = [10, 25, 50]
for N in stages_sw:
    fl = thrust_per_m(N, 0.10)
    vc = v_crossover(N, 0.10)
    print(f"    N={N:3d}: F/L = {fl:,.0f} N/m, v_cross = {vc/1000:.1f} km/s")

print()
print("  But the controller ALWAYS wants maximum thrust to minimize time.")
print("  It would immediately select the highest N (50) and never switch.")
print("  The lower stages are never used. Switching adds complexity for")
print("  zero benefit.")
print()

# Prove it: simulate fixed N=50 vs "switchable" starting at N=10
t_fixed, _, _, _ = simulate_launch(50, 0.10, M_TOTAL_CARGO,
                                    V_TARGET_CARGO, G_LIMIT_CARGO)

# Simulate a switching scenario: N=10 until 25 km/s, then N=50
# (this is strictly worse than N=50 throughout)
v = V_START
t = 0.0
N_lo, N_hi = 10, 50
vc_lo = v_crossover(N_lo, 0.10)
F_lo = thrust_per_m(N_lo, 0.10) * L_SLED
P_lo = F_lo * vc_lo
vc_hi = v_crossover(N_hi, 0.10)
F_hi = thrust_per_m(N_hi, 0.10) * L_SLED
P_hi = F_hi * vc_hi

while v < V_TARGET_CARGO and t < 72000:
    if v < 25000:
        F = F_lo if v <= vc_lo else P_lo / v
    else:
        F = F_hi if v <= vc_hi else P_hi / v
    a = F / M_TOTAL_CARGO
    v += a * DT
    t += DT

t_switch = t

print(f"  Fixed N=50 launch time:      {t_fixed/3600:.2f} hr")
print(f"  Switched N=10->50 time:      {t_switch/3600:.2f} hr  (SLOWER)")
print(f"  Switching penalty:           +{(t_switch - t_fixed)/60:.0f} min")
print()
print("  Conclusion: Cargo ALWAYS benefits from max N. Switching is useless.")

# =============================================================================
# PART 4: WHY SWITCHABLE WINDINGS DON'T HELP — CREWED CASE
# =============================================================================

print("\n" + "-" * 76)
print("PART 4: WHY SWITCHABLE WINDINGS DON'T HELP -- CREWED (3g limit)")
print("-" * 76)
print()
print("  Crewed launches are g-limited at 3g. The sled adjusts B_sled to")
print("  stay within the g-limit. Since P/L is invariant with B_sled,")
print("  the launch time doesn't change regardless of N.")
print()

N_crew_sweep = [50, 100, 200, 500]

print(f"{'N':>6}  {'v_cross':>10}  {'F_peak':>10}  {'a_peak':>8}  "
      f"{'t_launch':>10}  {'v_final':>10}")
print(f"{'':>6}  {'(km/s)':>10}  {'(MN)':>10}  {'(g)':>8}  "
      f"{'(hr)':>10}  {'(km/s)':>10}")
print("-" * 70)

for N in N_crew_sweep:
    B = 0.10
    fl = thrust_per_m(N, B)
    F_total = fl * L_SLED
    a_peak_g = F_total / M_TOTAL_CREW / G_0

    t_sim, v_sim, ok, _ = simulate_launch(N, B, M_TOTAL_CREW,
                                           V_TARGET_CREW, G_LIMIT_CREW)

    vc = v_crossover(N, B)
    print(f"{N:6d}  {vc/1000:10.1f}  {F_total/1e6:10.2f}  {a_peak_g:8.3f}  "
          f"{t_sim/3600:10.2f}  {v_sim/1000:10.1f}")

print()
print("  All cases reach 15 km/s in nearly identical time.")
print("  Higher N gives more thrust, but the g-limit caps it anyway.")
print("  The sled simply reduces B_sled to stay at 3g, and since P is")
print("  invariant with B_sled, the power envelope (and thus the launch")
print("  time) is the same regardless of N.")
print()
print("  Conclusion: Crewed missions don't care about N. Switching useless.")

# =============================================================================
# PART 5: POWER LIMIT SENSITIVITY
# =============================================================================

print("\n" + "-" * 76)
print("PART 5: HVDC POWER LIMIT SENSITIVITY (Cargo: 30 km/s, 10 km sled)")
print("-" * 76)
print()
print("  The ring's HVDC grid has finite capacity. How does a power cap")
print("  affect launch time for the heavy cargo reference case?")
print()

# Unlimited baseline
t_unlimited, _, _, _ = simulate_launch(50, 0.10, M_TOTAL_CARGO,
                                        V_TARGET_CARGO, G_LIMIT_CARGO,
                                        p_hvdc_max=None)

P_limits = [100e9, 200e9, 300e9, 416e9, 500e9, 666e9, None]

print(f"{'P_HVDC_MAX':>12}  {'t_launch':>10}  {'vs unlimited':>14}  {'limit engages':>16}")
print(f"{'(GW)':>12}  {'(hr)':>10}  {'(%)':>14}  {'(km/s)':>16}")
print("-" * 60)

for P_lim in P_limits:
    t_sim, v_sim, ok, pv = simulate_launch(50, 0.10, M_TOTAL_CARGO,
                                            V_TARGET_CARGO, G_LIMIT_CARGO,
                                            p_hvdc_max=P_lim)

    if P_lim is not None:
        p_str = f"{P_lim/1e9:12.0f}"
    else:
        p_str = f"{'unlimited':>12}"

    overhead = (t_sim - t_unlimited) / t_unlimited * 100

    if pv is not None:
        pv_str = f"{pv/1000:16.1f}"
    else:
        pv_str = f"{'N/A':>16}"

    print(f"{p_str}  {t_sim/3600:10.2f}  {overhead:13.1f}%  {pv_str}")

print()
P_voltage_limited = thrust_per_m(50, 0.10) * L_SLED * v_crossover(50, 0.10)
print(f"  Voltage-limited power (F_peak x v_cross): {P_voltage_limited/1e9:.1f} GW")
print(f"  At 666 GW (full ring capacity): heavy cargo is NOT power-limited")
print(f"  At 416 GW: matches the voltage-limited power exactly (no effect)")
print(f"  Below ~416 GW: launch time increases as grid constrains power")

# =============================================================================
# PART 6: CONCLUSION SUMMARY
# =============================================================================

print("\n" + "-" * 76)
print("PART 6: CONCLUSION")
print("-" * 76)
print()
print("  1. P/L = 41.6 MW/m is INVARIANT with N and B_sled (at I_PEAK)")
print("     (N and B_sled cancel in F * v_cross)")
print()
print("  2. N=50 gives only ~2.3% time overhead vs theoretical minimum")
print("     (v_cross = 5 km/s vs 30 km/s target)")
print()
print("  3. CARGO: Always wants maximum N to minimize launch time.")
print("     A switchable system would never switch down. Useless.")
print()
print("  4. CREWED: G-limited at 3g. The sled adjusts B_sled, not N.")
print("     Since P is invariant with B_sled, N doesn't matter. Useless.")
print()
print(f"  5. HVDC POWER LIMIT: At 666 GW (full ring capacity), the heavy"
      f"\n     cargo case ({P_voltage_limited/1e9:.0f} GW) is NOT power-limited."
      f"\n     Grid investment below ~{P_voltage_limited/1e9:.0f} GW starts to impact launch time.")
print()
print("  6. CONCLUSION: Fixed N=50 stator turns is optimal.")
print("     Simple construction. Easy to service. No switching needed.")
print("     The only trade-off (2.3% longer launch) is negligible.")
print()
print("=" * 76)
