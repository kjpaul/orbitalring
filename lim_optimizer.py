#!/usr/bin/env python3
"""
LIM Parameter Optimization Tool (with Power Limit)
===================================================
For: "Ion Propulsion Engineering" by Paul de Jong

This version correctly includes the POWER LIMIT constraint, which
dominates at high v_rel (end of deployment).

Key insight: The system transitions from thrust-limited (early) to 
power-limited (late) during deployment. The crossover happens when
the unconstrained thrust power exceeds the supply limit.

================================================================================
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================
MU0 = 4 * math.pi * 1e-7          # Permeability of free space: 1.257e-6 H/m
RHO_AL_106K = 8.35e-9             # Resistivity of aluminum at 106K: 8.35e-9 Ω·m
SIGMA_AL_106K = 1.0 / RHO_AL_106K # Conductivity: 1.20e8 S/m


# ==============================================================================
# DESIGN PARAMETERS
# ==============================================================================
@dataclass
class LIMParams:
    """
    LIM Design Parameters
    """
    # Coil geometry
    N: int = 100              # Turns per coil
    W: float = 2.0            # Coil width perpendicular to motion (m)
    tau_p: float = 50.0       # Pole pitch (m)
    pitch_count: int = 3      # Number of pole pitches per LIM
    
    # Gap and reaction plate
    gap: float = 0.05         # Air gap (m)
    t_plate: float = 0.08     # Reaction plate thickness (m)
    
    # Electrical
    I_target: float = 650.0   # Target peak current (A) - controller tries to reach this
    v_slip: float = 200.0     # Slip velocity (m/s)
    
    # HTS tape construction
    d_HTS: float = 80e-6      # HTS tape thickness: 80 μm
    d_kapton: float = 1e-3    # Kapton insulation: 1 mm
    
    # Safety limits
    I_c: float = 800.0        # HTS critical current (A)
    E_kapton_safe: float = 100e6  # Safe Kapton field: 100 kV/mm
    
    # POWER LIMIT - this is the key constraint!
    P_limit_site: float = 8e6     # Power limit per LIM site = 8 MW
    # Site has 2 LIMs, so per-LIM limit is P_limit_site / 2
    
    @property
    def P_limit_lim(self) -> float:
        """Power limit per single LIM (W)"""
        return self.P_limit_site / 2
    
    @property
    def k_fill(self) -> float:
        """Fill factor"""
        return self.d_HTS / (self.d_HTS + self.d_kapton)
    
    @property
    def A_coil(self) -> float:
        """Single coil area (m²)"""
        return self.tau_p * self.W
    
    @property
    def A_active(self) -> float:
        """Total active LIM area (m²)"""
        return self.tau_p * self.pitch_count * self.W
    
    @property
    def V_kapton_limit(self) -> float:
        """Maximum coil voltage from Kapton (V) - total coil, not per layer"""
        return self.E_kapton_safe * self.d_kapton


# ==============================================================================
# CORE PHYSICS FUNCTIONS
# ==============================================================================

def calc_B_field(p: LIMParams, I: float) -> float:
    """
    Magnetic field at reaction plate (T)
    
    B = (2μ₀NI)/(πW) × arctan(W/2g)
    """
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N * I / (math.pi * p.W)) * arctan_factor


def calc_inductance(p: LIMParams) -> float:
    """
    Coil inductance (H)
    
    L = (2μ₀N²×A_coil×k_fill)/(πW) × arctan(W/2g)
    """
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N**2 * p.A_coil * p.k_fill / (math.pi * p.W)) * arctan_factor


def calc_frequencies(p: LIMParams, v_rel: float) -> dict:
    """Calculate slip and supply frequencies"""
    v_wave = v_rel + p.v_slip
    f_supply = v_wave / (2 * p.tau_p)
    f_slip = p.v_slip / (2 * p.tau_p)
    slip = f_slip / f_supply if f_supply > 0 else 0
    return {
        'v_wave': v_wave,
        'f_supply': f_supply,
        'f_slip': f_slip,
        'slip': slip
    }


def calc_voltage(p: LIMParams, v_rel: float, I: float) -> float:
    """
    RMS coil voltage (V)
    
    V_rms = ωLI / √2
    """
    L = calc_inductance(p)
    freq = calc_frequencies(p, v_rel)
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
    """Slip efficiency η_slip = 2sG / (1 + s²G²)"""
    if G <= 0:
        return 0
    sG = slip * G
    return 2 * sG / (1 + sG**2)


def calc_thrust(p: LIMParams, v_rel: float, I: float) -> float:
    """
    Thrust force (N)
    
    F = (B²/2μ₀) × A_active × η_slip
    """
    B = calc_B_field(p, I)
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    
    magnetic_pressure = B**2 / (2 * MU0)
    return magnetic_pressure * p.A_active * eta


# ==============================================================================
# POWER-LIMITED OPERATION
# ==============================================================================

def calc_power_limited_current(p: LIMParams, v_rel: float) -> float:
    """
    Calculate the current needed to stay within power limit.
    
    At high v_rel, unconstrained thrust would exceed power limit.
    Controller must reduce current to maintain P = F × v_rel ≤ P_limit.
    
    Since F ∝ I², we have P ∝ I² × v_rel
    So I_limited = I_target × sqrt(P_limit / P_unconstrained)
    """
    # First, calculate unconstrained thrust at target current
    F_unconstrained = calc_thrust(p, v_rel, p.I_target)
    P_unconstrained = F_unconstrained * v_rel
    
    if P_unconstrained <= p.P_limit_lim:
        # Not power limited - use target current
        return p.I_target
    else:
        # Power limited - reduce current
        # F ∝ I², so to reduce power by factor k, reduce I by sqrt(k)
        reduction_factor = math.sqrt(p.P_limit_lim / P_unconstrained)
        I_limited = p.I_target * reduction_factor
        return I_limited


def find_crossover_velocity(p: LIMParams) -> float:
    """
    Find v_rel where system transitions from thrust-limited to power-limited.
    
    This is where P_unconstrained = P_limit
    """
    # Binary search
    v_low, v_high = 10, 10000
    
    while v_high - v_low > 1:
        v_mid = (v_low + v_high) / 2
        F = calc_thrust(p, v_mid, p.I_target)
        P = F * v_mid
        
        if P < p.P_limit_lim:
            v_low = v_mid
        else:
            v_high = v_mid
    
    return v_mid


# ==============================================================================
# COMPREHENSIVE EVALUATION
# ==============================================================================

def evaluate_design(p: LIMParams, v_rel: float, enforce_power_limit: bool = True) -> Dict:
    """
    Evaluate LIM design at a given operating point.
    
    If enforce_power_limit=True (default), current is reduced to stay within
    power budget at high v_rel. This matches real controller behavior.
    """
    # Determine operating current
    if enforce_power_limit:
        I_actual = calc_power_limited_current(p, v_rel)
    else:
        I_actual = p.I_target
    
    # Calculate all values at actual current
    B = calc_B_field(p, I_actual)
    L = calc_inductance(p)
    V = calc_voltage(p, v_rel, I_actual)
    F = calc_thrust(p, v_rel, I_actual)
    P_mech = F * v_rel
    
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    delta = calc_skin_depth(freq['f_slip'])
    
    # Also calculate unconstrained values for comparison
    F_unconstrained = calc_thrust(p, v_rel, p.I_target)
    P_unconstrained = F_unconstrained * v_rel
    
    # Constraint checks
    V_margin = (p.V_kapton_limit - V) / p.V_kapton_limit * 100
    I_margin = (p.I_c - I_actual) / p.I_c * 100
    P_margin = (p.P_limit_lim - P_mech) / p.P_limit_lim * 100
    
    is_power_limited = P_unconstrained > p.P_limit_lim
    
    return {
        # Operating point
        'I_actual_A': I_actual,
        'I_target_A': p.I_target,
        'I_reduction_%': (1 - I_actual / p.I_target) * 100,
        
        # Magnetic
        'B_mT': B * 1000,
        
        # Electrical
        'L_mH': L * 1000,
        'V_kV': V / 1000,
        'V_limit_kV': p.V_kapton_limit / 1000,
        'V_margin_%': V_margin,
        
        # Frequency
        'f_slip_Hz': freq['f_slip'],
        'f_supply_Hz': freq['f_supply'],
        'slip_%': freq['slip'] * 100,
        
        # Efficiency
        'G': G,
        'eta_%': eta * 100,
        'delta_mm': delta * 1000,
        
        # Performance (actual, power-limited)
        'F_kN': F / 1000,
        'P_mech_MW': P_mech / 1e6,
        'P_limit_MW': p.P_limit_lim / 1e6,
        'P_margin_%': P_margin,
        
        # Unconstrained (what physics allows)
        'F_unconstrained_kN': F_unconstrained / 1000,
        'P_unconstrained_MW': P_unconstrained / 1e6,
        
        # Status
        'power_limited': is_power_limited,
        'voltage_ok': V_margin > 0,
        'current_ok': I_margin > 0,
        'feasible': V_margin > 0 and I_margin > 0,
    }


# ==============================================================================
# OPTIMIZATION
# ==============================================================================

def find_optimal_N(p: LIMParams, v_rel: float, N_range: range = range(50, 501, 10)) -> Tuple[int, Dict]:
    """
    Find optimal turn count.
    
    With power limit enforced, higher N doesn't always help because
    current gets reduced anyway. Find the N that maximizes thrust
    while staying within ALL constraints.
    """
    best_N = N_range[0]
    best_thrust = 0
    best_result = None
    
    original_N = p.N
    
    for N in N_range:
        p.N = N
        result = evaluate_design(p, v_rel, enforce_power_limit=True)
        
        if result['feasible'] and result['F_kN'] > best_thrust:
            best_thrust = result['F_kN']
            best_N = N
            best_result = result
    
    p.N = original_N
    return best_N, best_result


# ==============================================================================
# OUTPUT FORMATTING
# ==============================================================================

def print_design_report(p: LIMParams, v_rel: float):
    """Print comprehensive design report"""
    result = evaluate_design(p, v_rel, enforce_power_limit=True)
    
    print()
    print("=" * 70)
    print("LIM DESIGN REPORT")
    print("=" * 70)
    
    print("\n┌─ DESIGN PARAMETERS ──────────────────────────────────────────────────")
    print(f"│ Turns per coil (N):           {p.N:>10}")
    print(f"│ Coil width (W):               {p.W:>10.2f} m")
    print(f"│ Pole pitch (τp):              {p.tau_p:>10.1f} m")
    print(f"│ Air gap (g):                  {p.gap*100:>10.1f} cm")
    print(f"│ Target current:               {p.I_target:>10.0f} A")
    print(f"│ Slip velocity:                {p.v_slip:>10.0f} m/s")
    print(f"│ Relative velocity:            {v_rel:>10.0f} m/s")
    print(f"│ Fill factor:                  {p.k_fill*100:>10.1f} %")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ LIMITS ─────────────────────────────────────────────────────────────")
    print(f"│ Power limit (per LIM):        {p.P_limit_lim/1e6:>10.1f} MW")
    print(f"│ Voltage limit (Kapton):       {p.V_kapton_limit/1000:>10.0f} kV")
    print(f"│ Current limit (Ic):           {p.I_c:>10.0f} A")
    print("└──────────────────────────────────────────────────────────────────────")
    
    # Highlight if power-limited
    if result['power_limited']:
        print("\n┌─ ⚠️  POWER LIMITED OPERATION ─────────────────────────────────────────")
        print(f"│ Actual current:               {result['I_actual_A']:>10.0f} A")
        print(f"│ Current reduction:            {result['I_reduction_%']:>10.1f} %")
        print(f"│ Unconstrained thrust:         {result['F_unconstrained_kN']:>10.1f} kN")
        print(f"│ Unconstrained power:          {result['P_unconstrained_MW']:>10.1f} MW")
        print("└──────────────────────────────────────────────────────────────────────")
    else:
        print("\n┌─ THRUST LIMITED OPERATION ────────────────────────────────────────────")
        print(f"│ Operating at target current:  {result['I_actual_A']:>10.0f} A")
        print(f"│ Power headroom:               {result['P_margin_%']:>10.1f} %")
        print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ ACTUAL PERFORMANCE ─────────────────────────────────────────────────")
    print(f"│ Magnetic field (B):           {result['B_mT']:>10.1f} mT")
    print(f"│ Inductance (L):               {result['L_mH']:>10.1f} mH")
    print(f"│ Coil voltage:                 {result['V_kV']:>10.2f} kV")
    print(f"│ Voltage margin:               {result['V_margin_%']:>10.1f} %")
    print(f"│ Supply frequency:             {result['f_supply_Hz']:>10.1f} Hz")
    print(f"│ Slip efficiency:              {result['eta_%']:>10.1f} %")
    print(f"│ Thrust:                       {result['F_kN']:>10.2f} kN")
    print(f"│ Mechanical power:             {result['P_mech_MW']:>10.2f} MW")
    print("└──────────────────────────────────────────────────────────────────────")
    
    status = "✓ FEASIBLE" if result['feasible'] else "✗ CONSTRAINT VIOLATED"
    mode = "POWER LIMITED" if result['power_limited'] else "THRUST LIMITED"
    print(f"\n  Status: {status} ({mode})")
    print("=" * 70)


def print_sweep_table(results: List[Dict], columns: List[str]):
    """Print formatted results table"""
    # Header
    header = " │ ".join(f"{c:>12}" for c in columns)
    separator = "─" * len(header)
    print(separator)
    print(header)
    print(separator)
    
    # Rows
    for r in results:
        row = []
        for c in columns:
            val = r.get(c, 0)
            if isinstance(val, bool):
                row.append(f"{'Yes':>12}" if val else f"{'No':>12}")
            elif isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    row.append(f"{val:>12.2e}")
                elif abs(val) > 10000:
                    row.append(f"{val:>12.0f}")
                else:
                    row.append(f"{val:>12.2f}")
            else:
                row.append(f"{val:>12}")
        print(" │ ".join(row))
    
    print(separator)


# ==============================================================================
# MAIN: OPTIMIZATION GUIDE
# ==============================================================================

def main():
    """Generate optimization guide with power limit"""
    
    print("\n" + "=" * 70)
    print("LIM PARAMETER OPTIMIZATION GUIDE (WITH POWER LIMIT)")
    print("For: Ion Propulsion Engineering")
    print("=" * 70)
    
    # Base parameters matching your code
    p = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        gap=0.05,
        I_target=650,
        v_slip=200,
        d_kapton=1e-3,
        E_kapton_safe=100e6,
        P_limit_site=8e6,  # 8 MW per site = 4 MW per LIM
    )
    
    # Find crossover velocity
    v_crossover = find_crossover_velocity(p)
    
    print(f"\n┌─ SYSTEM LIMITS ──────────────────────────────────────────────────────")
    print(f"│ Power limit per site:         {p.P_limit_site/1e6:>10.1f} MW")
    print(f"│ Power limit per LIM:          {p.P_limit_lim/1e6:>10.1f} MW")
    print(f"│ Voltage limit (Kapton):       {p.V_kapton_limit/1000:>10.0f} kV")
    print(f"│ Current limit (Ic):           {p.I_c:>10.0f} A")
    print(f"│")
    print(f"│ Crossover velocity:           {v_crossover:>10.0f} m/s")
    print(f"│   Below this: THRUST LIMITED (can use full current)")
    print(f"│   Above this: POWER LIMITED (must reduce current)")
    print("└──────────────────────────────────────────────────────────────────────")
    
    # ==== Show base design at different deployment stages ====
    print("\n" + "=" * 70)
    print("1. DEPLOYMENT PROGRESSION (with power limit enforced)")
    print("=" * 70)
    print(f"""
As v_rel increases during deployment:
- Frequency increases → voltage increases
- Thrust capacity increases (higher efficiency at lower slip ratio)
- BUT power limit caps the actual thrust at high v_rel

The controller reduces current to stay within the {p.P_limit_lim/1e6:.0f} MW limit.
""")
    
    v_rel_values = [100, 250, 500, 1000, 2000, 4000, 6000, 8000]
    results = []
    for v in v_rel_values:
        result = evaluate_design(p, v, enforce_power_limit=True)
        result['v_rel'] = v
        results.append(result)
    
    print(f"N = {p.N}, gap = {p.gap*100:.0f} cm, I_target = {p.I_target} A:")
    print_sweep_table(results, ['v_rel', 'I_actual_A', 'f_supply_Hz', 'V_kV', 
                                 'eta_%', 'F_kN', 'P_mech_MW', 'power_limited'])
    
    # ==== Show effect of N with power limit ====
    print("\n" + "=" * 70)
    print("2. EFFECT OF TURN COUNT (N) - with power limit")
    print("=" * 70)
    print("""
With power limit enforced, increasing N has DIMINISHING RETURNS:
- Higher N means higher unconstrained thrust
- But controller reduces current to stay within power limit
- So actual thrust is capped regardless of N

This is different from voltage-limited case!
""")
    
    print(f"\nAt v_rel = 8000 m/s (end of deployment, POWER LIMITED):")
    N_values = [50, 100, 150, 200, 250, 300]
    results = []
    original_N = p.N
    for N in N_values:
        p.N = N
        result = evaluate_design(p, 8000, enforce_power_limit=True)
        result['N'] = N
        results.append(result)
    p.N = original_N
    
    print_sweep_table(results, ['N', 'I_actual_A', 'B_mT', 'V_kV', 
                                 'F_kN', 'P_mech_MW', 'power_limited'])
    
    print(f"\nAt v_rel = 200 m/s (early deployment, THRUST LIMITED):")
    results = []
    for N in N_values:
        p.N = N
        result = evaluate_design(p, 200, enforce_power_limit=True)
        result['N'] = N
        results.append(result)
    p.N = original_N
    
    print_sweep_table(results, ['N', 'I_actual_A', 'B_mT', 'V_kV', 
                                 'F_kN', 'P_mech_MW', 'power_limited'])
    
    # ==== Slip velocity effect ====
    print("\n" + "=" * 70)
    print("3. EFFECT OF SLIP VELOCITY - with power limit")
    print("=" * 70)
    print("""
Lower slip velocity → higher efficiency → more thrust per amp.
But at high v_rel, you're power-limited anyway, so the benefit is capped.

The real benefit of lower slip is at EARLY deployment where you're
thrust-limited and every bit of efficiency helps.
""")
    
    print(f"\nAt v_rel = 200 m/s (THRUST LIMITED):")
    v_slip_values = [50, 100, 150, 200, 300, 400]
    results = []
    original_v_slip = p.v_slip
    for vs in v_slip_values:
        p.v_slip = vs
        result = evaluate_design(p, 200, enforce_power_limit=True)
        result['v_slip'] = vs
        results.append(result)
    p.v_slip = original_v_slip
    
    print_sweep_table(results, ['v_slip', 'slip_%', 'eta_%', 'F_kN', 
                                 'P_mech_MW', 'power_limited'])
    
    print(f"\nAt v_rel = 8000 m/s (POWER LIMITED):")
    results = []
    for vs in v_slip_values:
        p.v_slip = vs
        result = evaluate_design(p, 8000, enforce_power_limit=True)
        result['v_slip'] = vs
        results.append(result)
    p.v_slip = original_v_slip
    
    print_sweep_table(results, ['v_slip', 'slip_%', 'eta_%', 'F_kN', 
                                 'P_mech_MW', 'power_limited'])
    
    # ==== Full design reports ====
    print("\n" + "=" * 70)
    print("4. DETAILED DESIGN REPORTS")
    print("=" * 70)
    
    print("\n--- Early Deployment (v_rel = 200 m/s) ---")
    print_design_report(p, 200)
    
    print("\n--- Mid Deployment (v_rel = 2000 m/s) ---")
    print_design_report(p, 2000)
    
    print("\n--- End Deployment (v_rel = 8000 m/s) ---")
    print_design_report(p, 8000)
    
    # ==== Summary ====
    print("\n" + "=" * 70)
    print("OPTIMIZATION SUMMARY")
    print("=" * 70)
    print(f"""
KEY INSIGHT: The system has TWO operating regimes:

1. THRUST LIMITED (v_rel < {v_crossover:.0f} m/s)
   - Controller runs at target current ({p.I_target} A)
   - More turns (N) → more thrust
   - Voltage is the binding constraint
   
2. POWER LIMITED (v_rel > {v_crossover:.0f} m/s)  
   - Controller reduces current to stay within {p.P_limit_lim/1e6:.0f} MW
   - More turns (N) → controller reduces current more → same thrust
   - Power supply is the binding constraint

PRACTICAL IMPLICATIONS:

- Early deployment: Maximize N to get more thrust (voltage permitting)
- Late deployment: N doesn't matter much - you're power limited anyway
- The crossover at {v_crossover:.0f} m/s is when most of the momentum transfer
  shifts from "thrust-limited acceleration" to "power-limited cruising"

RECOMMENDED APPROACH:
1. Choose N based on voltage limit at end of deployment
2. Accept that late deployment will be power-limited
3. Focus optimization on early deployment performance
4. Consider whether power limit can be increased (more solar panels?)
""")


if __name__ == "__main__":
    main()