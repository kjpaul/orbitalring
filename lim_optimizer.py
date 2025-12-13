#!/usr/bin/env python3
"""
LIM Parameter Optimization Tool
================================
Uses corrected arctan formulas for B-field and inductance.

Key insight: With the corrected formulas, increasing N no longer causes
voltage to blow up as fast, because both B and L saturate at small gaps.

Author: Paul de Jong / Claude
For: Ion Propulsion Engineering
"""

import math
import numpy as np
from dataclasses import dataclass
from typing import Tuple, List, Dict

# Physical constants
MU0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)
SIGMA_AL_106K = 1.0 / 8.35e-9  # Conductivity of Al at 106K (S/m)

@dataclass
class LIMParams:
    """LIM design parameters"""
    N: int = 100            # Turns per coil
    W: float = 2.0          # Coil width (m)
    tau_p: float = 50.0     # Pole pitch (m)
    gap: float = 0.05       # Air gap (m)
    I_peak: float = 650.0   # Peak current (A)
    v_slip: float = 200.0   # Slip velocity (m/s)
    t_plate: float = 0.08   # Reaction plate thickness (m)
    d_HTS: float = 80e-6    # HTS tape thickness (m)
    d_kapton: float = 1e-3  # Kapton thickness (m)
    pitch_count: int = 3    # Number of pole pitches per LIM
    
    # Limits
    I_c: float = 800.0      # Critical current (A)
    V_max: float = 100e3    # Maximum voltage (V)
    E_kapton: float = 50e6  # Safe Kapton field strength (V/m) - 10 kV/mm
    
    @property
    def k_fill(self) -> float:
        """Fill factor"""
        return self.d_HTS / (self.d_HTS + self.d_kapton)
    
    @property
    def A_coil(self) -> float:
        """Coil area (m²)"""
        return self.tau_p * self.W
    
    @property
    def A_active(self) -> float:
        """Active LIM area (m²)"""
        return self.tau_p * self.pitch_count * self.W
    
    @property
    def V_breakdown(self) -> float:
        """Voltage limit from Kapton breakdown (V)"""
        return (self.N - 1) * self.E_kapton * self.d_kapton


def calc_B_field(p: LIMParams) -> float:
    """
    Corrected magnetic field at reaction plate (T)
    
    B = (2μ₀NI)/(πW) × arctan(W/2g)
    
    Limiting behavior:
    - gap → 0: B → μ₀NI/W (finite, saturates)
    - gap → ∞: B → μ₀NI/(πg) (recovers 1/g scaling)
    """
    return (2 * MU0 * p.N * p.I_peak / (math.pi * p.W)) * math.atan(p.W / (2 * p.gap))


def calc_inductance(p: LIMParams) -> float:
    """
    Corrected coil inductance (H)
    
    L = (2μ₀N²A_coil×k_fill)/(πW) × arctan(W/2g)
    
    Limiting behavior:
    - gap → 0: L → μ₀N²A_coil×k_fill/W (finite)
    - gap → ∞: L → μ₀N²A_coil×k_fill/(πg)
    """
    return (2 * MU0 * p.N**2 * p.A_coil * p.k_fill / (math.pi * p.W)) * math.atan(p.W / (2 * p.gap))


def calc_frequency(p: LIMParams, v_rel: float) -> Tuple[float, float]:
    """
    Calculate slip and supply frequencies (Hz)
    
    Returns: (f_slip, f_supply)
    """
    f_slip = p.v_slip / (2 * p.tau_p)
    v_wave = v_rel + p.v_slip
    f_supply = v_wave / (2 * p.tau_p)
    return f_slip, f_supply


def calc_voltage(p: LIMParams, v_rel: float) -> float:
    """
    RMS coil voltage (V)
    
    V_rms = ωLI/√2 = 2πf × L × I / √2
    """
    L = calc_inductance(p)
    _, f_supply = calc_frequency(p, v_rel)
    omega = 2 * math.pi * f_supply
    return omega * L * p.I_peak / math.sqrt(2)


def calc_skin_depth(f_slip: float, rho: float = 8.35e-9) -> float:
    """Skin depth in aluminum (m)"""
    if f_slip <= 0:
        return 0.1  # Return large value for DC
    return math.sqrt(rho / (math.pi * MU0 * f_slip))


def calc_goodness_factor(p: LIMParams, f_slip: float) -> float:
    """
    Goodness factor G
    
    G = (ω_slip × μ₀ × σ × δ × τp) / π
    
    Large G means efficient energy transfer to rotor.
    """
    if f_slip <= 0:
        return 0
    omega_slip = 2 * math.pi * f_slip
    delta = calc_skin_depth(f_slip)
    delta_eff = min(delta, p.t_plate)
    return (omega_slip * MU0 * SIGMA_AL_106K * delta_eff * p.tau_p) / math.pi


def calc_slip_efficiency(slip: float, G: float) -> float:
    """
    Slip efficiency η_slip
    
    η_slip = 2sG / (1 + s²G²)
    
    Maximum η = 1 occurs at s = 1/G when G >> 1
    """
    if G <= 0:
        return 0
    denominator = 1 + (slip * G)**2
    return 2 * slip * G / denominator


def calc_thrust(p: LIMParams, v_rel: float) -> float:
    """
    Thrust force (N)
    
    F = (B²/2μ₀) × A_active × η_slip
    """
    B = calc_B_field(p)
    f_slip, f_supply = calc_frequency(p, v_rel)
    
    if f_supply <= 0:
        return 0
    
    slip = f_slip / f_supply
    G = calc_goodness_factor(p, f_slip)
    eta = calc_slip_efficiency(slip, G)
    
    magnetic_pressure = B**2 / (2 * MU0)
    return magnetic_pressure * p.A_active * eta


def calc_thrust_power(p: LIMParams, v_rel: float) -> float:
    """Mechanical power delivered (W)"""
    return calc_thrust(p, v_rel) * v_rel


def evaluate_design(p: LIMParams, v_rel: float) -> Dict:
    """
    Comprehensive evaluation of a LIM design point.
    
    Returns dict with all key metrics and constraint status.
    """
    B = calc_B_field(p)
    L = calc_inductance(p)
    V = calc_voltage(p, v_rel)
    F = calc_thrust(p, v_rel)
    P = calc_thrust_power(p, v_rel)
    
    f_slip, f_supply = calc_frequency(p, v_rel)
    slip = f_slip / f_supply if f_supply > 0 else 0
    G = calc_goodness_factor(p, f_slip)
    eta = calc_slip_efficiency(slip, G)
    delta = calc_skin_depth(f_slip)
    
    # Constraint margins
    V_limit = min(p.V_max, p.V_breakdown)
    V_margin = (V_limit - V) / V_limit * 100
    I_margin = (p.I_c - p.I_peak) / p.I_c * 100
    
    return {
        'B_mT': B * 1000,
        'L_mH': L * 1000,
        'V_kV': V / 1000,
        'V_limit_kV': V_limit / 1000,
        'V_margin_%': V_margin,
        'I_margin_%': I_margin,
        'F_kN': F / 1000,
        'P_MW': P / 1e6,
        'f_slip_Hz': f_slip,
        'f_supply_Hz': f_supply,
        'slip_%': slip * 100,
        'G': G,
        'eta_%': eta * 100,
        'delta_mm': delta * 1000,
        'k_fill': p.k_fill,
        'feasible': V_margin > 0 and I_margin > 0
    }


def find_optimal_N(p: LIMParams, v_rel: float, N_range: range = range(50, 501, 10)) -> Tuple[int, Dict]:
    """
    Find optimal turn count that maximizes thrust while staying within voltage limits.
    
    Returns: (optimal_N, results_dict)
    """
    best_N = N_range[0]
    best_thrust = 0
    best_result = None
    
    for N in N_range:
        p.N = N
        result = evaluate_design(p, v_rel)
        
        if result['feasible'] and result['F_kN'] > best_thrust:
            best_thrust = result['F_kN']
            best_N = N
            best_result = result
    
    p.N = best_N
    return best_N, best_result


def parameter_sweep(base_params: LIMParams, v_rel: float, 
                    param_name: str, values: List[float]) -> List[Dict]:
    """
    Sweep a single parameter and return results for each value.
    """
    results = []
    p = LIMParams(**vars(base_params))  # Copy
    
    for val in values:
        setattr(p, param_name, val)
        result = evaluate_design(p, v_rel)
        result[param_name] = val
        results.append(result)
    
    return results


def print_sweep_table(results: List[Dict], param_name: str, 
                      columns: List[str] = None):
    """Print a formatted table of sweep results."""
    if columns is None:
        columns = [param_name, 'B_mT', 'L_mH', 'V_kV', 'V_limit_kV', 
                   'V_margin_%', 'F_kN', 'P_MW', 'eta_%']
    
    # Header
    header = " | ".join(f"{c:>12}" for c in columns)
    print(header)
    print("-" * len(header))
    
    # Rows
    for r in results:
        row = []
        for c in columns:
            val = r.get(c, 0)
            if isinstance(val, float):
                if abs(val) < 0.01 or abs(val) > 1000:
                    row.append(f"{val:>12.2e}")
                else:
                    row.append(f"{val:>12.2f}")
            else:
                row.append(f"{val:>12}")
        print(" | ".join(row))


def print_design_summary(p: LIMParams, v_rel: float):
    """Print a complete design summary."""
    result = evaluate_design(p, v_rel)
    
    print("=" * 60)
    print("LIM DESIGN SUMMARY")
    print("=" * 60)
    print(f"\nDesign Parameters:")
    print(f"  Turns (N):          {p.N}")
    print(f"  Coil width (W):     {p.W} m")
    print(f"  Pole pitch (τp):    {p.tau_p} m")
    print(f"  Air gap (g):        {p.gap*100:.1f} cm")
    print(f"  Peak current:       {p.I_peak} A")
    print(f"  Slip velocity:      {p.v_slip} m/s")
    print(f"  Relative velocity:  {v_rel} m/s")
    print(f"  Fill factor:        {p.k_fill:.4f}")
    print(f"  Kapton thickness:   {p.d_kapton*1000:.2f} mm")
    
    print(f"\nCalculated Values:")
    print(f"  Magnetic field:     {result['B_mT']:.2f} mT")
    print(f"  Inductance:         {result['L_mH']:.2f} mH")
    print(f"  Coil voltage:       {result['V_kV']:.2f} kV")
    print(f"  Voltage limit:      {result['V_limit_kV']:.0f} kV")
    print(f"  Voltage margin:     {result['V_margin_%']:.1f}%")
    
    print(f"\nPerformance:")
    print(f"  Thrust:             {result['F_kN']:.2f} kN")
    print(f"  Thrust power:       {result['P_MW']:.3f} MW")
    print(f"  Slip frequency:     {result['f_slip_Hz']:.1f} Hz")
    print(f"  Supply frequency:   {result['f_supply_Hz']:.1f} Hz")
    print(f"  Goodness factor:    {result['G']:.0f}")
    print(f"  Slip efficiency:    {result['eta_%']:.1f}%")
    print(f"  Skin depth:         {result['delta_mm']:.1f} mm")
    
    status = "✓ FEASIBLE" if result['feasible'] else "✗ INFEASIBLE"
    print(f"\nStatus: {status}")
    print("=" * 60)


def generate_optimization_guide(v_rel_start: float = 100, v_rel_end: float = 8000):
    """
    Generate comprehensive optimization guidance for the book.
    Shows how to select parameters at different deployment stages.
    """
    print("\n" + "=" * 70)
    print("LIM PARAMETER OPTIMIZATION GUIDE")
    print("=" * 70)
    
    # Base parameters
    p = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        gap=0.05,
        I_peak=650,
        v_slip=200,
        d_kapton=1e-3,  # 1mm Kapton
    )
    
    print("\n1. EFFECT OF TURN COUNT (N)")
    print("-" * 40)
    print(f"   At v_rel = {v_rel_end} m/s (end of deployment)")
    print()
    
    N_values = [50, 100, 150, 200, 250, 300, 400, 500]
    results = parameter_sweep(p, v_rel_end, 'N', N_values)
    print_sweep_table(results, 'N', ['N', 'B_mT', 'V_kV', 'V_limit_kV', 'V_margin_%', 'F_kN'])
    
    # Find optimal
    p_opt = LIMParams(**vars(p))
    opt_N, opt_result = find_optimal_N(p_opt, v_rel_end, range(50, 501, 10))
    print(f"\n   → Optimal N = {opt_N} (max thrust within voltage limit)")
    
    print("\n2. EFFECT OF GAP")
    print("-" * 40)
    print(f"   At N = {p.N}, v_rel = {v_rel_end} m/s")
    print()
    
    gap_values = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
    results = parameter_sweep(p, v_rel_end, 'gap', gap_values)
    # Convert gap to cm for display
    for r in results:
        r['gap_cm'] = r['gap'] * 100
    print_sweep_table(results, 'gap_cm', ['gap_cm', 'B_mT', 'L_mH', 'V_kV', 'F_kN', 'eta_%'])
    
    print("\n   Key insight: Smaller gap → higher B and thrust, but also higher L and V")
    print("   With corrected formulas, B and L saturate as gap→0 (no blowup)")
    
    print("\n3. EFFECT OF SLIP VELOCITY")
    print("-" * 40)
    print(f"   At N = {p.N}, gap = {p.gap*100:.0f} cm, v_rel = {v_rel_end} m/s")
    print()
    
    v_slip_values = [50, 100, 150, 200, 300, 400, 500]
    results = parameter_sweep(p, v_rel_end, 'v_slip', v_slip_values)
    print_sweep_table(results, 'v_slip', ['v_slip', 'f_slip_Hz', 'f_supply_Hz', 'slip_%', 'G', 'eta_%', 'F_kN', 'V_kV'])
    
    print("\n   Key insight: Higher slip → higher frequency → higher voltage")
    print("   But also higher slip efficiency (up to a point)")
    
    print("\n4. DEPLOYMENT PROGRESSION")
    print("-" * 40)
    print(f"   N = {p.N}, gap = {p.gap*100:.0f} cm, v_slip = {p.v_slip} m/s")
    print()
    
    v_rel_values = [100, 500, 1000, 2000, 4000, 6000, 8000]
    results = []
    for v in v_rel_values:
        result = evaluate_design(p, v)
        result['v_rel'] = v
        results.append(result)
    print_sweep_table(results, 'v_rel', ['v_rel', 'f_supply_Hz', 'V_kV', 'slip_%', 'eta_%', 'F_kN', 'P_MW'])
    
    print("\n   Key insight: Voltage increases with v_rel (frequency increases)")
    print("   Design must work at END of deployment (highest v_rel)")
    
    print("\n5. KAPTON THICKNESS TRADE-OFF")
    print("-" * 40)
    print()
    
    print("   Thicker Kapton → higher voltage limit, but lower k_fill")
    print()
    d_kapton_values = [0.1e-3, 0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3]
    results = []
    for d in d_kapton_values:
        p.d_kapton = d
        result = evaluate_design(p, v_rel_end)
        result['d_kapton_mm'] = d * 1000
        result['k_fill'] = p.k_fill
        results.append(result)
    print_sweep_table(results, 'd_kapton_mm', ['d_kapton_mm', 'k_fill', 'V_limit_kV', 'L_mH', 'V_kV', 'V_margin_%'])
    
    # Reset
    p.d_kapton = 1e-3
    
    print("\n" + "=" * 70)
    print("OPTIMIZATION STRATEGY")
    print("=" * 70)
    print("""
1. START with voltage constraint at END of deployment (v_rel = 8000 m/s)
   - This is when frequency and voltage are highest

2. CHOOSE Kapton thickness first:
   - Thicker = higher voltage headroom, but lower fill factor
   - 1 mm is a good starting point (k_fill ≈ 0.07, V_limit = 100 kV)

3. SELECT turn count (N) to maximize thrust within voltage limit:
   - Use find_optimal_N() function
   - More turns = more thrust, but also more voltage

4. VERIFY at START of deployment (v_rel = 100 m/s):
   - Thrust will be lower (lower efficiency at low slip)
   - But voltage will also be lower (more margin)

5. ADJUST slip velocity for thermal management:
   - Higher slip = more eddy losses = more heating
   - May need to reduce slip at start when thermal margin is tight
""")


def main():
    """Run example optimization analysis."""
    
    # Example: Find optimal design for end of deployment
    print("\n" + "=" * 70)
    print("EXAMPLE: Optimizing for End of Deployment")
    print("=" * 70)
    
    p = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        gap=0.05,
        I_peak=650,
        v_slip=200,
        d_kapton=1e-3,
    )
    
    v_rel = 8000  # End of deployment
    
    print(f"\nStarting parameters: N={p.N}, gap={p.gap*100:.0f} cm")
    print_design_summary(p, v_rel)
    
    print("\nSearching for optimal N...")
    opt_N, opt_result = find_optimal_N(p, v_rel)
    p.N = opt_N
    print_design_summary(p, v_rel)
    
    # Generate full guide
    generate_optimization_guide()


if __name__ == "__main__":
    main()