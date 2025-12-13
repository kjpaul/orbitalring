#!/usr/bin/env python3
"""
LIM Design Parameter Optimizer
==============================
For: "Ion Propulsion Engineering" by Paul de Jong

PURPOSE: Find optimal values for key design parameters:
  1. Number of turns (N)
  2. Coil width (W)
  3. Pole pitch (τp)
  4. Kapton thickness (d_kapton)

CONSTRAINTS:
  - Voltage limit (Kapton breakdown)
  - Current limit (HTS critical current)
  - Power limit (site electrical capacity)

================================================================================
"""

import math
from dataclasses import dataclass, field
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
    """LIM Design Parameters - the knobs we want to optimize"""
    
    # === PRIMARY DESIGN VARIABLES (what we're optimizing) ===
    N: int = 100                  # Turns per coil
    W: float = 2.0                # Coil width (m)
    tau_p: float = 50.0           # Pole pitch (m)
    d_kapton: float = 1e-3        # Kapton thickness (m)
    
    # === SECONDARY PARAMETERS (fixed or constrained) ===
    pitch_count: int = 3          # Pole pitches per LIM
    gap: float = 0.05             # Air gap (m)
    t_plate: float = 0.08         # Reaction plate thickness (m)
    d_HTS: float = 80e-6          # HTS tape thickness (m)
    
    # === OPERATING PARAMETERS ===
    I_target: float = 650.0       # Target current (A)
    v_slip: float = 200.0         # Slip velocity (m/s)
    
    # === LIMITS ===
    I_c: float = 800.0            # HTS critical current (A)
    E_kapton_safe: float = 100e6  # Safe Kapton field: 100 kV/mm = 100e6 V/m
    P_limit_site: float = 8e6     # Site power limit (W)
    
    # === DERIVED PROPERTIES ===
    @property
    def P_limit_lim(self) -> float:
        """Power limit per LIM (W)"""
        return self.P_limit_site / 2
    
    @property
    def k_fill(self) -> float:
        """Fill factor: fraction of coil that is superconductor"""
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
    def V_limit(self) -> float:
        """Voltage limit from Kapton (V) - total coil voltage"""
        return self.E_kapton_safe * self.d_kapton
    
    @property
    def L_HTS_coil(self) -> float:
        """HTS tape length per coil (m)"""
        return 2 * (self.W + self.tau_p) * self.N
    
    @property
    def L_HTS_lim(self) -> float:
        """Total HTS tape per LIM (m)"""
        return self.L_HTS_coil * 3 * self.pitch_count  # 3 phases


# ==============================================================================
# PHYSICS CALCULATIONS
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
    """RMS coil voltage (V): V = ωLI/√2"""
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
    """Slip efficiency η = 2sG / (1 + s²G²)"""
    if G <= 0:
        return 0
    sG = slip * G
    return 2 * sG / (1 + sG**2)


def calc_thrust(p: LIMParams, v_rel: float, I: float) -> float:
    """Thrust force (N): F = (B²/2μ₀) × A_active × η"""
    B = calc_B_field(p, I)
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    magnetic_pressure = B**2 / (2 * MU0)
    return magnetic_pressure * p.A_active * eta


def calc_power_limited_current(p: LIMParams, v_rel: float) -> Tuple[float, bool]:
    """
    Calculate actual operating current considering power limit.
    Returns (I_actual, is_power_limited)
    """
    F_unconstrained = calc_thrust(p, v_rel, p.I_target)
    P_unconstrained = F_unconstrained * v_rel
    
    if P_unconstrained <= p.P_limit_lim:
        return p.I_target, False
    else:
        reduction = math.sqrt(p.P_limit_lim / P_unconstrained)
        return p.I_target * reduction, True


# ==============================================================================
# DESIGN EVALUATION
# ==============================================================================

def evaluate_design(p: LIMParams, v_rel: float) -> Dict:
    """
    Comprehensive evaluation of a LIM design at given operating point.
    """
    # Get power-limited current
    I_actual, power_limited = calc_power_limited_current(p, v_rel)
    
    # Core calculations
    B = calc_B_field(p, I_actual)
    L = calc_inductance(p)
    V = calc_voltage(p, v_rel, I_actual)
    F = calc_thrust(p, v_rel, I_actual)
    P_mech = F * v_rel
    
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    delta = calc_skin_depth(freq['f_slip'])
    
    # Unconstrained values
    F_unconstrained = calc_thrust(p, v_rel, p.I_target)
    P_unconstrained = F_unconstrained * v_rel
    
    # Margins
    V_margin = (p.V_limit - V) / p.V_limit * 100 if p.V_limit > 0 else -100
    I_margin = (p.I_c - I_actual) / p.I_c * 100
    
    return {
        # Design parameters
        'N': p.N,
        'W': p.W,
        'tau_p': p.tau_p,
        'd_kapton_mm': p.d_kapton * 1000,
        'k_fill_%': p.k_fill * 100,
        'A_active': p.A_active,
        'L_HTS_km': p.L_HTS_lim / 1000,
        
        # Operating point
        'v_rel': v_rel,
        'v_slip': p.v_slip,
        'I_actual': I_actual,
        'I_target': p.I_target,
        
        # Electrical
        'B_mT': B * 1000,
        'L_mH': L * 1000,
        'V_kV': V / 1000,
        'V_limit_kV': p.V_limit / 1000,
        'V_margin_%': V_margin,
        
        # Frequency
        'f_slip_Hz': freq['f_slip'],
        'f_supply_Hz': freq['f_supply'],
        'slip_%': freq['slip'] * 100,
        
        # Efficiency
        'G': G,
        'eta_%': eta * 100,
        'delta_mm': delta * 1000,
        
        # Performance
        'F_kN': F / 1000,
        'F_unconstrained_kN': F_unconstrained / 1000,
        'P_mech_MW': P_mech / 1e6,
        'P_unconstrained_MW': P_unconstrained / 1e6,
        'P_limit_MW': p.P_limit_lim / 1e6,
        
        # Status
        'power_limited': power_limited,
        'voltage_ok': V_margin > 0,
        'feasible': V_margin > 0 and I_margin > 0,
    }


# ==============================================================================
# PARAMETER SWEEPS
# ==============================================================================

def sweep_parameter(base: LIMParams, param_name: str, values: List, 
                    v_rel: float) -> List[Dict]:
    """Sweep a single parameter and return results"""
    results = []
    for val in values:
        p = LIMParams(**{k: v for k, v in vars(base).items()})
        setattr(p, param_name, val)
        result = evaluate_design(p, v_rel)
        results.append(result)
    return results


def find_optimal_N(p: LIMParams, v_rel: float, N_range: range = range(50, 501, 10)) -> int:
    """Find N that maximizes thrust within voltage limit"""
    best_N = N_range[0]
    best_thrust = 0
    
    original_N = p.N
    for N in N_range:
        p.N = N
        result = evaluate_design(p, v_rel)
        if result['feasible'] and result['F_kN'] > best_thrust:
            best_thrust = result['F_kN']
            best_N = N
    p.N = original_N
    return best_N


# ==============================================================================
# OUTPUT FORMATTING
# ==============================================================================

def print_table(results: List[Dict], columns: List[Tuple[str, str, str]]):
    """
    Print formatted table.
    columns: list of (key, header, format_spec)
    """
    # Build header
    headers = [h for _, h, _ in columns]
    widths = [max(len(h), 10) for h in headers]
    
    header_row = " │ ".join(f"{h:>{w}}" for h, w in zip(headers, widths))
    separator = "─" * len(header_row)
    
    print(separator)
    print(header_row)
    print(separator)
    
    for r in results:
        row_parts = []
        for (key, _, fmt), w in zip(columns, widths):
            val = r.get(key, 0)
            if isinstance(val, bool):
                s = "Yes" if val else "No"
            elif isinstance(val, float):
                s = f"{val:{fmt}}"
            else:
                s = str(val)
            row_parts.append(f"{s:>{w}}")
        print(" │ ".join(row_parts))
    
    print(separator)


def print_design_summary(p: LIMParams, v_rel: float, title: str = ""):
    """Print a design summary box"""
    r = evaluate_design(p, v_rel)
    
    print()
    print("=" * 70)
    if title:
        print(title)
        print("=" * 70)
    
    print(f"""
┌─ DESIGN PARAMETERS ──────────────────────────────────────────────────
│ Turns (N):              {p.N:>10}        Pole pitch (τp):    {p.tau_p:>8.1f} m
│ Coil width (W):         {p.W:>10.2f} m     Kapton thickness:   {p.d_kapton*1000:>8.2f} mm
│ Air gap:                {p.gap*100:>10.1f} cm    Fill factor:        {p.k_fill*100:>8.1f} %
│ Active area:            {p.A_active:>10.0f} m²    HTS per LIM:        {p.L_HTS_lim/1000:>8.1f} km
└──────────────────────────────────────────────────────────────────────

┌─ OPERATING POINT ────────────────────────────────────────────────────
│ Relative velocity:      {v_rel:>10.0f} m/s   Slip velocity:      {p.v_slip:>8.0f} m/s
│ Current (actual):       {r['I_actual']:>10.0f} A     Current (target):   {p.I_target:>8.0f} A
│ Supply frequency:       {r['f_supply_Hz']:>10.1f} Hz    Slip frequency:     {r['f_slip_Hz']:>8.2f} Hz
└──────────────────────────────────────────────────────────────────────

┌─ PERFORMANCE ────────────────────────────────────────────────────────
│ Magnetic field:         {r['B_mT']:>10.1f} mT    Goodness factor:    {r['G']:>8.0f}
│ Inductance:             {r['L_mH']:>10.1f} mH    Slip efficiency:    {r['eta_%']:>8.1f} %
│ Coil voltage:           {r['V_kV']:>10.2f} kV    Voltage limit:      {r['V_limit_kV']:>8.0f} kV
│ Voltage margin:         {r['V_margin_%']:>10.1f} %     
│ Thrust:                 {r['F_kN']:>10.2f} kN    Power:              {r['P_mech_MW']:>8.2f} MW
│ Power limit:            {r['P_limit_MW']:>10.1f} MW    Power limited:      {'Yes' if r['power_limited'] else 'No':>8}
└──────────────────────────────────────────────────────────────────────
""")
    status = "✓ FEASIBLE" if r['feasible'] else "✗ INFEASIBLE"
    print(f"  Status: {status}")
    print("=" * 70)


# ==============================================================================
# MAIN OPTIMIZATION ANALYSIS
# ==============================================================================

def main():
    print("\n" + "=" * 70)
    print("LIM DESIGN PARAMETER OPTIMIZATION")
    print("For: Ion Propulsion Engineering")
    print("=" * 70)
    
    # Baseline design
    baseline = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        d_kapton=1e-3,
        gap=0.05,
        I_target=650,
        v_slip=200,
        P_limit_site=8e6,
    )
    
    # We'll analyze at two key points: early and late deployment
    v_early = 500    # Thrust-limited regime
    v_late = 8000    # Power-limited regime
    
    print(f"""
ANALYSIS APPROACH:
- Sweep each design parameter independently
- Evaluate at v_rel = {v_early} m/s (early deployment, thrust-limited)
- Evaluate at v_rel = {v_late} m/s (late deployment, power-limited)
- Find optimal values considering both regimes
""")
    
    # =========================================================================
    # 1. NUMBER OF TURNS (N)
    # =========================================================================
    print("\n" + "=" * 70)
    print("1. NUMBER OF TURNS (N)")
    print("=" * 70)
    print("""
Effect of N:
  - B-field ∝ N (more turns = stronger field)
  - Inductance ∝ N² (increases faster than B)
  - Voltage ∝ N² × f (inductance × frequency)
  - Thrust ∝ B² ∝ N² (when not power-limited)
  - HTS tape ∝ N (more turns = more tape)

Trade-off: More turns gives more thrust, but voltage increases faster.
The voltage limit caps how many turns you can use.
""")
    
    N_values = [50, 100, 150, 200, 250, 300, 350, 400]
    
    print(f"\nAt v_rel = {v_early} m/s (EARLY deployment, thrust-limited):")
    results = sweep_parameter(baseline, 'N', N_values, v_early)
    print_table(results, [
        ('N', 'N', '.0f'),
        ('B_mT', 'B (mT)', '.1f'),
        ('L_mH', 'L (mH)', '.1f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('V_limit_kV', 'V_lim', '.0f'),
        ('V_margin_%', 'Margin%', '.0f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('P_mech_MW', 'P (MW)', '.2f'),
        ('L_HTS_km', 'HTS(km)', '.1f'),
        ('feasible', 'OK?', ''),
    ])
    
    print(f"\nAt v_rel = {v_late} m/s (LATE deployment, power-limited):")
    results = sweep_parameter(baseline, 'N', N_values, v_late)
    print_table(results, [
        ('N', 'N', '.0f'),
        ('I_actual', 'I (A)', '.0f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('V_margin_%', 'Margin%', '.0f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('P_mech_MW', 'P (MW)', '.2f'),
        ('power_limited', 'P_lim?', ''),
        ('feasible', 'OK?', ''),
    ])
    
    opt_N_early = find_optimal_N(baseline, v_early)
    opt_N_late = find_optimal_N(baseline, v_late)
    print(f"\n→ Optimal N at early deployment: {opt_N_early}")
    print(f"→ Optimal N at late deployment: {opt_N_late}")
    print(f"→ Recommendation: Use N = {min(opt_N_early, opt_N_late)} (limited by late deployment voltage)")
    
    # =========================================================================
    # 2. COIL WIDTH (W)
    # =========================================================================
    print("\n" + "=" * 70)
    print("2. COIL WIDTH (W)")
    print("=" * 70)
    print("""
Effect of W:
  - B-field: complex (arctan term), roughly ∝ 1/W for W >> gap
  - Active area ∝ W (wider = more area)
  - Inductance ∝ W (through A_coil)
  - Thrust: B² × A, net effect depends on geometry

Trade-off: Wider coils have more area but weaker field per amp.
There's an optimal width that balances these effects.
""")
    
    W_values = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    
    print(f"\nAt v_rel = {v_early} m/s (thrust-limited):")
    results = sweep_parameter(baseline, 'W', W_values, v_early)
    print_table(results, [
        ('W', 'W (m)', '.1f'),
        ('A_active', 'Area', '.0f'),
        ('B_mT', 'B (mT)', '.1f'),
        ('L_mH', 'L (mH)', '.1f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('L_HTS_km', 'HTS(km)', '.1f'),
        ('feasible', 'OK?', ''),
    ])
    
    print(f"\nAt v_rel = {v_late} m/s (power-limited):")
    results = sweep_parameter(baseline, 'W', W_values, v_late)
    print_table(results, [
        ('W', 'W (m)', '.1f'),
        ('I_actual', 'I (A)', '.0f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('P_mech_MW', 'P (MW)', '.2f'),
        ('power_limited', 'P_lim?', ''),
    ])
    
    # =========================================================================
    # 3. POLE PITCH (τp)
    # =========================================================================
    print("\n" + "=" * 70)
    print("3. POLE PITCH (τp)")
    print("=" * 70)
    print("""
Effect of τp:
  - Wave velocity: v_sync = 2 × τp × f_supply
  - Frequency: f = v_wave / (2 × τp), so larger τp = lower frequency
  - Lower frequency = lower voltage (V ∝ f)
  - Goodness factor ∝ τp (larger pitch = better efficiency)
  - Active area ∝ τp (more area per LIM)
  - But: longer LIM = more material, more complex

Trade-off: Longer pole pitch gives lower frequency (less voltage),
better efficiency, but requires longer LIM structures.
""")
    
    tau_p_values = [20, 30, 40, 50, 60, 75, 100]
    
    print(f"\nAt v_rel = {v_early} m/s (thrust-limited):")
    results = sweep_parameter(baseline, 'tau_p', tau_p_values, v_early)
    print_table(results, [
        ('tau_p', 'τp (m)', '.0f'),
        ('f_supply_Hz', 'f (Hz)', '.1f'),
        ('G', 'G', '.0f'),
        ('eta_%', 'η (%)', '.1f'),
        ('L_mH', 'L (mH)', '.0f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('A_active', 'Area', '.0f'),
    ])
    
    print(f"\nAt v_rel = {v_late} m/s (power-limited):")
    results = sweep_parameter(baseline, 'tau_p', tau_p_values, v_late)
    print_table(results, [
        ('tau_p', 'τp (m)', '.0f'),
        ('f_supply_Hz', 'f (Hz)', '.1f'),
        ('G', 'G', '.0f'),
        ('eta_%', 'η (%)', '.1f'),
        ('I_actual', 'I (A)', '.0f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('F_kN', 'F (kN)', '.2f'),
    ])
    
    # =========================================================================
    # 4. KAPTON THICKNESS (d_kapton)
    # =========================================================================
    print("\n" + "=" * 70)
    print("4. KAPTON THICKNESS (d_kapton)")
    print("=" * 70)
    print("""
Effect of d_kapton:
  - Voltage limit = E_safe × d_kapton (thicker = higher limit)
  - Fill factor k_fill = d_HTS / (d_HTS + d_kapton) (thicker = lower k_fill)
  - Inductance ∝ k_fill (lower k_fill = lower inductance = lower voltage)
  - Net effect: thicker Kapton allows higher voltage limit AND
    reduces actual voltage (double benefit for voltage margin)

Trade-off: Thicker Kapton is better for voltage, but means less
superconductor per unit volume (lower k_fill).
""")
    
    d_kapton_values = [0.25e-3, 0.5e-3, 0.75e-3, 1.0e-3, 1.5e-3, 2.0e-3]
    
    print(f"\nAt v_rel = {v_late} m/s (where voltage matters most):")
    results = sweep_parameter(baseline, 'd_kapton', d_kapton_values, v_late)
    print_table(results, [
        ('d_kapton_mm', 'd (mm)', '.2f'),
        ('k_fill_%', 'k_fill%', '.1f'),
        ('V_limit_kV', 'V_lim', '.0f'),
        ('L_mH', 'L (mH)', '.1f'),
        ('V_kV', 'V (kV)', '.2f'),
        ('V_margin_%', 'Margin%', '.0f'),
        ('F_kN', 'F (kN)', '.2f'),
        ('feasible', 'OK?', ''),
    ])
    
    # =========================================================================
    # BASELINE DESIGN SUMMARY
    # =========================================================================
    print("\n" + "=" * 70)
    print("BASELINE DESIGN SUMMARY")
    print("=" * 70)
    
    print_design_summary(baseline, v_early, f"Early Deployment (v_rel = {v_early} m/s)")
    print_design_summary(baseline, v_late, f"Late Deployment (v_rel = {v_late} m/s)")
    
    # =========================================================================
    # OPTIMIZATION SUMMARY
    # =========================================================================
    print("\n" + "=" * 70)
    print("OPTIMIZATION SUMMARY")
    print("=" * 70)
    print(f"""
BASELINE PARAMETERS:
  N = {baseline.N} turns
  W = {baseline.W} m
  τp = {baseline.tau_p} m
  d_kapton = {baseline.d_kapton*1000} mm

KEY FINDINGS:

1. NUMBER OF TURNS (N):
   - More turns → more thrust (when not power-limited)
   - But voltage ∝ N², so limited by Kapton breakdown
   - At v_rel = {v_late} m/s: max N ≈ {opt_N_late} before voltage limit
   - Increasing N doesn't help when power-limited (controller reduces current)

2. COIL WIDTH (W):
   - Wider → more active area, but weaker field
   - Sweet spot around W = 2-3 m for these parameters
   - Wider coils need more HTS tape

3. POLE PITCH (τp):
   - Longer τp → lower frequency → lower voltage
   - Also better efficiency (higher G)
   - Trade-off: longer structures, more material
   - τp = 50 m is a good compromise

4. KAPTON THICKNESS (d_kapton):
   - Thicker → higher voltage limit AND lower actual voltage
   - 1 mm gives good margin with reasonable k_fill (7.4%)
   - Could go thinner if voltage margin is comfortable

POWER LIMIT IS THE DOMINANT CONSTRAINT AT HIGH v_rel:
   At v_rel = {v_late} m/s, thrust = P_limit / v_rel = {baseline.P_limit_lim/1e6:.1f} MW / {v_late} m/s = {baseline.P_limit_lim/v_late/1000:.2f} kN
   This is independent of N, W, τp - only power supply size matters!

RECOMMENDATIONS:
   - Keep N = {baseline.N}-300 (voltage-limited at high v_rel)
   - Keep W = 2 m (good balance of field and area)
   - Keep τp = 50 m (practical size, good efficiency)
   - Keep d_kapton = 1 mm (safe voltage margin)
   - To get more thrust at high v_rel: increase power limit (more solar panels)
""")


if __name__ == "__main__":
    main()