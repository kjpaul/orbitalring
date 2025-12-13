#!/usr/bin/env python3
"""
LIM Parameter Optimization Tool (v4 - Adaptive Slip Strategy)
==============================================================
For: "Ion Propulsion Engineering" by Paul de Jong

This version explores the adaptive slip velocity strategy:
- When THRUST LIMITED: Use higher slip for stability
- When POWER LIMITED: Reduce slip to increase efficiency

Key insight: At low slip, efficiency is much higher, so more of your
limited electrical power becomes useful thrust instead of heat.

================================================================================
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================
MU0 = 4 * math.pi * 1e-7
RHO_AL_106K = 8.35e-9
SIGMA_AL_106K = 1.0 / RHO_AL_106K


# ==============================================================================
# DESIGN PARAMETERS
# ==============================================================================
@dataclass
class LIMParams:
    """LIM Design Parameters"""
    # Geometry
    N: int = 100
    W: float = 2.0
    tau_p: float = 50.0
    pitch_count: int = 3
    gap: float = 0.05
    t_plate: float = 0.08
    
    # Electrical
    I_target: float = 650.0
    v_slip: float = 200.0       # Default/max slip velocity
    v_slip_min: float = 50.0    # Minimum slip velocity for control stability
    
    # Tape construction
    d_HTS: float = 80e-6
    d_kapton: float = 1e-3
    
    # Limits
    I_c: float = 800.0
    E_kapton_safe: float = 100e6
    P_limit_site: float = 8e6
    
    @property
    def P_limit_lim(self) -> float:
        return self.P_limit_site / 2
    
    @property
    def k_fill(self) -> float:
        return self.d_HTS / (self.d_HTS + self.d_kapton)
    
    @property
    def A_coil(self) -> float:
        return self.tau_p * self.W
    
    @property
    def A_active(self) -> float:
        return self.tau_p * self.pitch_count * self.W
    
    @property
    def V_kapton_limit(self) -> float:
        return self.E_kapton_safe * self.d_kapton


# ==============================================================================
# CORE PHYSICS
# ==============================================================================

def calc_B_field(p: LIMParams, I: float) -> float:
    """Magnetic field at reaction plate (T)"""
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N * I / (math.pi * p.W)) * arctan_factor


def calc_inductance(p: LIMParams) -> float:
    """Coil inductance (H)"""
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N**2 * p.A_coil * p.k_fill / (math.pi * p.W)) * arctan_factor


def calc_frequencies(p: LIMParams, v_rel: float, v_slip: float) -> dict:
    """Calculate frequencies for given slip velocity"""
    v_wave = v_rel + v_slip
    f_supply = v_wave / (2 * p.tau_p)
    f_slip = v_slip / (2 * p.tau_p)
    slip = f_slip / f_supply if f_supply > 0 else 0
    return {
        'v_wave': v_wave,
        'f_supply': f_supply,
        'f_slip': f_slip,
        'slip': slip
    }


def calc_voltage(p: LIMParams, v_rel: float, v_slip: float, I: float) -> float:
    """RMS coil voltage (V)"""
    L = calc_inductance(p)
    freq = calc_frequencies(p, v_rel, v_slip)
    omega = 2 * math.pi * freq['f_supply']
    return omega * L * I / math.sqrt(2)


def calc_skin_depth(f_slip: float) -> float:
    """Skin depth in aluminum (m)"""
    if f_slip <= 0:
        return 0.1
    return math.sqrt(RHO_AL_106K / (math.pi * MU0 * f_slip))


def calc_goodness_factor(p: LIMParams, f_slip: float) -> float:
    """Goodness factor G"""
    if f_slip <= 0:
        return 0
    omega_slip = 2 * math.pi * f_slip
    delta = calc_skin_depth(f_slip)
    delta_eff = min(delta, p.t_plate)
    return (omega_slip * MU0 * SIGMA_AL_106K * delta_eff * p.tau_p) / math.pi


def calc_slip_efficiency(slip: float, G: float) -> float:
    """Slip efficiency η = 2sG / (1 + s²G²)"""
    if G <= 0:
        return 0
    sG = slip * G
    return 2 * sG / (1 + sG**2)


def calc_thrust(p: LIMParams, v_rel: float, v_slip: float, I: float) -> float:
    """Thrust force (N)"""
    B = calc_B_field(p, I)
    freq = calc_frequencies(p, v_rel, v_slip)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    magnetic_pressure = B**2 / (2 * MU0)
    return magnetic_pressure * p.A_active * eta


def calc_eddy_losses(p: LIMParams, v_slip: float, I: float) -> float:
    """
    Approximate eddy current losses in reaction plate (W)
    
    P_eddy ∝ B² × f² × volume
    """
    B = calc_B_field(p, I)
    f_slip = v_slip / (2 * p.tau_p)
    delta = calc_skin_depth(f_slip)
    delta_eff = min(delta, p.t_plate)
    
    # P = (π²/6ρ) × B² × δ² × f² × V
    ce = math.pi**2 / (6 * RHO_AL_106K)
    V_eddy = p.W * p.A_active * delta_eff
    
    return ce * B**2 * delta_eff**2 * f_slip**2 * V_eddy


# ==============================================================================
# ADAPTIVE CONTROL STRATEGIES
# ==============================================================================

def evaluate_fixed_slip(p: LIMParams, v_rel: float) -> Dict:
    """
    STRATEGY 1: Fixed slip velocity
    
    Controller keeps v_slip constant, reduces current when power-limited.
    This is the baseline/simple approach.
    """
    v_slip = p.v_slip
    
    # Calculate unconstrained thrust
    F_unconstrained = calc_thrust(p, v_rel, v_slip, p.I_target)
    P_unconstrained = F_unconstrained * v_rel
    
    # Determine if power-limited
    if P_unconstrained <= p.P_limit_lim:
        I_actual = p.I_target
        power_limited = False
    else:
        # Reduce current: F ∝ I², so I_new = I_target × sqrt(P_limit/P_unconstrained)
        reduction = math.sqrt(p.P_limit_lim / P_unconstrained)
        I_actual = p.I_target * reduction
        power_limited = True
    
    # Calculate actual performance
    F = calc_thrust(p, v_rel, v_slip, I_actual)
    P_mech = F * v_rel
    V = calc_voltage(p, v_rel, v_slip, I_actual)
    P_eddy = calc_eddy_losses(p, v_slip, I_actual)
    
    freq = calc_frequencies(p, v_rel, v_slip)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    
    return {
        'strategy': 'Fixed Slip',
        'v_slip': v_slip,
        'I_actual': I_actual,
        'F_kN': F / 1000,
        'P_mech_MW': P_mech / 1e6,
        'P_eddy_kW': P_eddy / 1000,
        'V_kV': V / 1000,
        'eta_%': eta * 100,
        'slip_%': freq['slip'] * 100,
        'G': G,
        'power_limited': power_limited,
    }


def evaluate_adaptive_slip(p: LIMParams, v_rel: float) -> Dict:
    """
    STRATEGY 2: Adaptive slip velocity
    
    When power-limited, REDUCE slip velocity to increase efficiency.
    This extracts more thrust from the same power budget.
    
    Controller logic:
    1. Start with max slip (p.v_slip)
    2. If power-limited, reduce slip to improve efficiency
    3. Stop at minimum slip (p.v_slip_min) for control stability
    """
    # First check if we're power-limited at max slip
    F_max_slip = calc_thrust(p, v_rel, p.v_slip, p.I_target)
    P_max_slip = F_max_slip * v_rel
    
    if P_max_slip <= p.P_limit_lim:
        # Not power-limited - use max slip for stability
        v_slip = p.v_slip
        I_actual = p.I_target
        power_limited = False
    else:
        # Power-limited - find optimal slip
        # Try reducing slip to get more thrust from power budget
        
        best_thrust = 0
        best_v_slip = p.v_slip_min
        
        # Search from min slip to max slip
        for v_slip_test in range(int(p.v_slip_min), int(p.v_slip) + 1, 5):
            # At this slip, what thrust can we get within power limit?
            F_test = calc_thrust(p, v_rel, v_slip_test, p.I_target)
            P_test = F_test * v_rel
            
            if P_test <= p.P_limit_lim:
                # Can run at full current
                I_test = p.I_target
                F_actual = F_test
            else:
                # Must reduce current
                reduction = math.sqrt(p.P_limit_lim / P_test)
                I_test = p.I_target * reduction
                F_actual = calc_thrust(p, v_rel, v_slip_test, I_test)
            
            if F_actual > best_thrust:
                best_thrust = F_actual
                best_v_slip = v_slip_test
                best_I = I_test
        
        v_slip = best_v_slip
        I_actual = best_I
        power_limited = True
    
    # Calculate actual performance at chosen operating point
    F = calc_thrust(p, v_rel, v_slip, I_actual)
    P_mech = F * v_rel
    V = calc_voltage(p, v_rel, v_slip, I_actual)
    P_eddy = calc_eddy_losses(p, v_slip, I_actual)
    
    freq = calc_frequencies(p, v_rel, v_slip)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    
    return {
        'strategy': 'Adaptive Slip',
        'v_slip': v_slip,
        'I_actual': I_actual,
        'F_kN': F / 1000,
        'P_mech_MW': P_mech / 1e6,
        'P_eddy_kW': P_eddy / 1000,
        'V_kV': V / 1000,
        'eta_%': eta * 100,
        'slip_%': freq['slip'] * 100,
        'G': G,
        'power_limited': power_limited,
    }


def find_crossover_velocity(p: LIMParams) -> float:
    """Find v_rel where system transitions from thrust to power limited"""
    v_low, v_high = 10, 10000
    while v_high - v_low > 1:
        v_mid = (v_low + v_high) / 2
        F = calc_thrust(p, v_mid, p.v_slip, p.I_target)
        P = F * v_mid
        if P < p.P_limit_lim:
            v_low = v_mid
        else:
            v_high = v_mid
    return v_mid


# ==============================================================================
# OUTPUT FORMATTING
# ==============================================================================

def print_comparison_table(results_fixed: List[Dict], results_adaptive: List[Dict], 
                           v_rel_values: List[float]):
    """Print side-by-side comparison of strategies"""
    
    print("\n" + "─" * 100)
    print(f"{'v_rel':>8} │ {'FIXED SLIP':^40} │ {'ADAPTIVE SLIP':^40} │ {'Thrust':>8}")
    print(f"{'(m/s)':>8} │ {'v_slip':>8} {'I':>6} {'F':>8} {'η':>6} {'P_eddy':>8} │ {'v_slip':>8} {'I':>6} {'F':>8} {'η':>6} {'P_eddy':>8} │ {'Gain':>8}")
    print("─" * 100)
    
    for i, v in enumerate(v_rel_values):
        rf = results_fixed[i]
        ra = results_adaptive[i]
        
        gain = (ra['F_kN'] / rf['F_kN'] - 1) * 100 if rf['F_kN'] > 0 else 0
        
        mode_f = "*" if rf['power_limited'] else " "
        mode_a = "*" if ra['power_limited'] else " "
        
        print(f"{v:>8.0f} │ {rf['v_slip']:>7.0f}{mode_f} {rf['I_actual']:>6.0f} {rf['F_kN']:>8.2f} {rf['eta_%']:>5.1f}% {rf['P_eddy_kW']:>7.0f} │ " +
              f"{ra['v_slip']:>7.0f}{mode_a} {ra['I_actual']:>6.0f} {ra['F_kN']:>8.2f} {ra['eta_%']:>5.1f}% {ra['P_eddy_kW']:>7.0f} │ {gain:>+7.1f}%")
    
    print("─" * 100)
    print("* = power limited")


def print_strategy_report(p: LIMParams, v_rel: float):
    """Print detailed comparison at a single operating point"""
    
    rf = evaluate_fixed_slip(p, v_rel)
    ra = evaluate_adaptive_slip(p, v_rel)
    
    print()
    print("=" * 70)
    print(f"STRATEGY COMPARISON at v_rel = {v_rel} m/s")
    print("=" * 70)
    
    print("\n┌─ FIXED SLIP STRATEGY ────────────────────────────────────────────────")
    print(f"│ Slip velocity:                {rf['v_slip']:>10.0f} m/s")
    print(f"│ Current:                      {rf['I_actual']:>10.0f} A")
    print(f"│ Slip ratio:                   {rf['slip_%']:>10.2f} %")
    print(f"│ Efficiency:                   {rf['eta_%']:>10.1f} %")
    print(f"│ Thrust:                       {rf['F_kN']:>10.2f} kN")
    print(f"│ Mechanical power:             {rf['P_mech_MW']:>10.2f} MW")
    print(f"│ Eddy losses:                  {rf['P_eddy_kW']:>10.0f} kW")
    print(f"│ Coil voltage:                 {rf['V_kV']:>10.2f} kV")
    mode = "POWER LIMITED" if rf['power_limited'] else "THRUST LIMITED"
    print(f"│ Mode:                         {mode:>10}")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ ADAPTIVE SLIP STRATEGY ─────────────────────────────────────────────")
    print(f"│ Slip velocity:                {ra['v_slip']:>10.0f} m/s")
    print(f"│ Current:                      {ra['I_actual']:>10.0f} A")
    print(f"│ Slip ratio:                   {ra['slip_%']:>10.2f} %")
    print(f"│ Efficiency:                   {ra['eta_%']:>10.1f} %")
    print(f"│ Thrust:                       {ra['F_kN']:>10.2f} kN")
    print(f"│ Mechanical power:             {ra['P_mech_MW']:>10.2f} MW")
    print(f"│ Eddy losses:                  {ra['P_eddy_kW']:>10.0f} kW")
    print(f"│ Coil voltage:                 {ra['V_kV']:>10.2f} kV")
    mode = "POWER LIMITED" if ra['power_limited'] else "THRUST LIMITED"
    print(f"│ Mode:                         {mode:>10}")
    print("└──────────────────────────────────────────────────────────────────────")
    
    if ra['F_kN'] > rf['F_kN']:
        gain = (ra['F_kN'] / rf['F_kN'] - 1) * 100
        print(f"\n  → Adaptive slip gives {gain:.0f}% MORE THRUST for same power!")
        print(f"  → Eddy losses reduced by {(1 - ra['P_eddy_kW']/rf['P_eddy_kW'])*100:.0f}%")
    else:
        print(f"\n  → Both strategies give same result (thrust limited)")
    
    print("=" * 70)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "=" * 70)
    print("ADAPTIVE SLIP VELOCITY STRATEGY")
    print("For: Ion Propulsion Engineering")
    print("=" * 70)
    
    p = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        gap=0.05,
        I_target=650,
        v_slip=200,          # Max slip (for stability when thrust-limited)
        v_slip_min=50,       # Min slip (control limit when power-limited)
        d_kapton=1e-3,
        P_limit_site=8e6,
    )
    
    v_crossover = find_crossover_velocity(p)
    
    print(f"""
THE PROBLEM:
When power-limited (v_rel > {v_crossover:.0f} m/s), the controller reduces current
to stay within the {p.P_limit_lim/1e6:.0f} MW limit. But at high slip velocity, 
efficiency is low (~8%), so most power goes to heating the reaction plate.

THE SOLUTION:
Reduce slip velocity when power-limited. Lower slip → higher efficiency.
More of the limited power budget becomes useful thrust.

CONSTRAINTS:
- Max slip: {p.v_slip} m/s (used when thrust-limited, for stability)
- Min slip: {p.v_slip_min} m/s (control stability limit)
""")
    
    # ==== Deployment progression comparison ====
    print("\n" + "=" * 70)
    print("1. DEPLOYMENT PROGRESSION: Fixed vs Adaptive Slip")
    print("=" * 70)
    
    v_rel_values = [100, 250, 500, 1000, 1500, 2000, 3000, 4000, 6000, 8000]
    
    results_fixed = [evaluate_fixed_slip(p, v) for v in v_rel_values]
    results_adaptive = [evaluate_adaptive_slip(p, v) for v in v_rel_values]
    
    print_comparison_table(results_fixed, results_adaptive, v_rel_values)
    
    # ==== Detailed comparison at key points ====
    print("\n" + "=" * 70)
    print("2. DETAILED COMPARISON AT KEY OPERATING POINTS")
    print("=" * 70)
    
    print("\n--- Early Deployment (v_rel = 500 m/s) - THRUST LIMITED ---")
    print_strategy_report(p, 500)
    
    print("\n--- Mid Deployment (v_rel = 2000 m/s) - POWER LIMITED ---")
    print_strategy_report(p, 2000)
    
    print("\n--- Late Deployment (v_rel = 8000 m/s) - POWER LIMITED ---")
    print_strategy_report(p, 8000)
    
    # ==== Total deployment benefit ====
    print("\n" + "=" * 70)
    print("3. CUMULATIVE BENEFIT OVER DEPLOYMENT")
    print("=" * 70)
    
    # Integrate thrust over velocity to estimate momentum transfer
    total_impulse_fixed = 0
    total_impulse_adaptive = 0
    total_eddy_fixed = 0
    total_eddy_adaptive = 0
    
    v_step = 100
    for v in range(100, 8001, v_step):
        rf = evaluate_fixed_slip(p, v)
        ra = evaluate_adaptive_slip(p, v)
        
        # Impulse ∝ F × time, and time ∝ 1/F for same delta-v
        # So compare F directly at each velocity
        total_impulse_fixed += rf['F_kN']
        total_impulse_adaptive += ra['F_kN']
        total_eddy_fixed += rf['P_eddy_kW']
        total_eddy_adaptive += ra['P_eddy_kW']
    
    benefit = (total_impulse_adaptive / total_impulse_fixed - 1) * 100
    eddy_reduction = (1 - total_eddy_adaptive / total_eddy_fixed) * 100
    
    print(f"""
Summing thrust contributions across deployment (v_rel = 100 to 8000 m/s):

  Fixed slip total:      {total_impulse_fixed:>10.1f} kN (relative units)
  Adaptive slip total:   {total_impulse_adaptive:>10.1f} kN (relative units)
  
  THRUST IMPROVEMENT:    {benefit:>+10.1f} %
  
  Fixed slip eddy heat:  {total_eddy_fixed:>10.0f} kW (relative units)
  Adaptive slip heat:    {total_eddy_adaptive:>10.0f} kW (relative units)
  
  HEAT REDUCTION:        {eddy_reduction:>+10.1f} %
""")
    
    # ==== Implementation notes ====
    print("\n" + "=" * 70)
    print("4. IMPLEMENTATION NOTES")
    print("=" * 70)
    print(f"""
CONTROL LOGIC:

1. Measure v_rel (cable-casing relative velocity)
2. Calculate unconstrained thrust at max slip ({p.v_slip} m/s)
3. If P_thrust < P_limit ({p.P_limit_lim/1e6:.0f} MW):
   - Run at max slip ({p.v_slip} m/s) for stability
   - Use full target current ({p.I_target} A)
4. Else (power-limited):
   - Reduce slip velocity toward {p.v_slip_min} m/s
   - Maintain current at {p.I_target} A if possible
   - Only reduce current if still over power limit at min slip

BENEFITS:
- More thrust from same power budget
- Less heating of reaction plate
- Faster deployment

RISKS:
- Low slip is harder to control
- Disturbances could push into regenerative braking
- Need robust slip feedback control

RECOMMENDED MINIMUM SLIP: {p.v_slip_min} m/s
- Provides {100/calc_goodness_factor(p, p.v_slip_min/(2*p.tau_p)):.1f}× safety margin above optimal slip
- Still practical for control systems
""")


if __name__ == "__main__":
    main()