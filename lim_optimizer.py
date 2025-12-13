#!/usr/bin/env python3
"""
LIM Parameter Optimization Tool
================================
For: "Ion Propulsion Engineering" by Paul de Jong

This tool helps find optimal Linear Induction Motor parameters for
orbital ring deployment. All formulas use the corrected arctan model
that prevents division-by-zero at small gaps.

Key Design Constraints:
1. Voltage must stay below Kapton insulation limit
2. Current must stay below HTS critical current
3. Thrust must be sufficient for deployment timeline
4. Losses must be manageable for thermal system

================================================================================
PHYSICAL MODEL OVERVIEW
================================================================================

The LIM accelerates the orbital ring cable through electromagnetic induction.
A traveling magnetic wave in the stator (coils) induces eddy currents in the
rotor (aluminum reaction plate), creating thrust.

Key relationships:
- More turns (N) → stronger B-field → more thrust, but also higher voltage
- Smaller gap → stronger B-field → more thrust, but harder mechanically
- Higher slip → more eddy losses → more heating
- Higher frequency → higher voltage (V ∝ ωL)

================================================================================
"""

import math
from dataclasses import dataclass
from typing import Tuple, List, Dict

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================
MU0 = 4 * math.pi * 1e-7      # Permeability of free space: 1.257e-6 H/m
RHO_AL_106K = 8.35e-9         # Resistivity of aluminum at 106K: 8.35e-9 Ω·m
SIGMA_AL_106K = 1.0 / RHO_AL_106K  # Conductivity: 1.20e8 S/m


# ==============================================================================
# DESIGN PARAMETERS
# ==============================================================================
@dataclass
class LIMParams:
    """
    LIM Design Parameters
    
    These are the knobs you can turn to optimize your design.
    """
    # Coil geometry
    N: int = 100              # Turns per coil (more turns = stronger field)
    W: float = 2.0            # Coil width perpendicular to motion (m)
    tau_p: float = 50.0       # Pole pitch - distance for one magnetic cycle (m)
    pitch_count: int = 3      # Number of pole pitches per LIM
    
    # Gap and reaction plate
    gap: float = 0.05         # Air gap between coil and reaction plate (m)
    t_plate: float = 0.08     # Reaction plate thickness (m)
    
    # Electrical
    I_peak: float = 650.0     # Peak current in coils (A)
    v_slip: float = 200.0     # Slip velocity: how much slower plate moves than wave (m/s)
    
    # HTS tape construction
    d_HTS: float = 80e-6      # HTS tape thickness: 80 μm (m)
    d_kapton: float = 1e-3    # Kapton insulation thickness: 1 mm (m)
    
    # Safety limits
    I_c: float = 800.0        # HTS critical current - must not exceed (A)
    E_kapton_safe: float = 100e6  # Safe Kapton field: 100 kV/mm = 100e6 V/m
                                   # (Breakdown is ~200 kV/mm, this is 2× safety factor)
    
    # ==== Derived Properties ====
    
    @property
    def k_fill(self) -> float:
        """
        Fill Factor: fraction of coil cross-section that is superconductor.
        
        k_fill = d_HTS / (d_HTS + d_kapton)
        
        Thicker Kapton → lower k_fill → lower inductance → lower voltage
        But also less superconductor per unit volume.
        
        Typical values:
        - 1mm Kapton: k_fill = 0.074 (7.4%)
        - 0.1mm Kapton: k_fill = 0.44 (44%)
        """
        return self.d_HTS / (self.d_HTS + self.d_kapton)
    
    @property
    def A_coil(self) -> float:
        """Single coil area (m²)"""
        return self.tau_p * self.W
    
    @property
    def A_active(self) -> float:
        """Total active LIM area (m²) = tau_p × pitch_count × W"""
        return self.tau_p * self.pitch_count * self.W
    
    @property
    def V_kapton_limit(self) -> float:
        """
        Maximum Coil Voltage from Kapton Insulation (V)
        
        This is a TOTAL COIL voltage limit, not per-layer.
        
        Rationale: In a fault or transient, the full coil voltage could
        appear across a single insulation layer. Design conservatively.
        
        V_limit = E_safe × d_kapton
        
        Example: 100 kV/mm × 1 mm = 100 kV
        """
        return self.E_kapton_safe * self.d_kapton


# ==============================================================================
# MAGNETIC FIELD
# ==============================================================================
def calc_B_field(p: LIMParams) -> float:
    """
    Magnetic Field at Reaction Plate Surface (T)
    
    CORRECTED FORMULA (prevents infinity at small gaps):
    
        B = (2μ₀NI)/(πW) × arctan(W/2g)
    
    Physical interpretation:
    - The arctan factor accounts for the finite width of the current sheet
    - As gap → 0: arctan(W/2g) → π/2, so B → μ₀NI/W (finite!)
    - As gap → ∞: arctan(W/2g) → W/2g, so B → μ₀NI/(πg) (1/g falloff)
    
    OLD (WRONG) FORMULA for reference:
        B = μ₀NI/g  ← This blows up as g → 0!
    
    Typical values for orbital ring LIM:
    - N=100, I=650A, W=2m, g=5cm → B ≈ 40 mT
    - This is weaker than a refrigerator magnet (~5 mT at surface)
    - But it acts over a huge area (300 m²), giving substantial force
    """
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N * p.I_peak / (math.pi * p.W)) * arctan_factor


# ==============================================================================
# INDUCTANCE
# ==============================================================================
def calc_inductance(p: LIMParams) -> float:
    """
    Coil Inductance (H)
    
    CORRECTED FORMULA (prevents infinity at small gaps):
    
        L = (2μ₀N²×A_coil×k_fill)/(πW) × arctan(W/2g)
    
    Physical interpretation:
    - Inductance is flux linkage per ampere: L = NΦ/I
    - The arctan factor comes from integrating the B-field
    - Fill factor k_fill accounts for only part of coil being conductor
    
    As gap → 0: L → μ₀N²×A_coil×k_fill/W (finite!)
    As gap → ∞: L → μ₀N²×A_coil×k_fill/(πg) (1/g falloff)
    
    Why inductance matters:
    - Voltage is V = ωLI = 2πfLI
    - Higher L → higher voltage at same frequency
    - Voltage is limited by Kapton insulation
    """
    arctan_factor = math.atan(p.W / (2 * p.gap))
    return (2 * MU0 * p.N**2 * p.A_coil * p.k_fill / (math.pi * p.W)) * arctan_factor


# ==============================================================================
# FREQUENCY AND SLIP
# ==============================================================================
def calc_frequencies(p: LIMParams, v_rel: float) -> dict:
    """
    Calculate Slip and Supply Frequencies
    
    The traveling magnetic wave moves at v_wave = 2 × tau_p × f_supply
    
    Definitions:
    - v_rel: Relative velocity between cable and casing (m/s)
    - v_slip: How much slower the cable moves than the wave (m/s)
    - v_wave: Speed of the magnetic wave = v_rel + v_slip
    
    Frequencies:
    - f_supply = v_wave / (2 × tau_p)  [Hz] - what the power supply produces
    - f_slip = v_slip / (2 × tau_p)    [Hz] - frequency of induced currents
    
    Slip ratio:
    - s = v_slip / v_wave = f_slip / f_supply
    - s = 0: No thrust (wave moves with cable)
    - s = 1: Stalled (wave moves, cable stationary)
    - Typical: s = 0.02 to 0.05 (2-5%)
    """
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


# ==============================================================================
# VOLTAGE
# ==============================================================================
def calc_voltage(p: LIMParams, v_rel: float) -> float:
    """
    RMS Coil Voltage (V)
    
    Formula:
        V_rms = ωLI / √2 = 2πf × L × I / √2
    
    This is the induced voltage from changing flux linkage.
    
    Why it increases with v_rel:
    - Higher v_rel → higher f_supply → higher ω
    - V ∝ ω, so voltage increases throughout deployment
    - Design must work at END of deployment (highest v_rel = ~8000 m/s)
    
    Constraint:
    - Must stay below V_kapton_limit to avoid insulation breakdown
    """
    L = calc_inductance(p)
    freq = calc_frequencies(p, v_rel)
    omega = 2 * math.pi * freq['f_supply']
    return omega * L * p.I_peak / math.sqrt(2)


# ==============================================================================
# SKIN DEPTH
# ==============================================================================
def calc_skin_depth(f_slip: float, rho: float = RHO_AL_106K) -> float:
    """
    Skin Depth in Reaction Plate (m)
    
    Formula:
        δ = √(ρ / (π × μ₀ × f))
    
    Physical meaning:
    - Eddy currents concentrate near the surface
    - δ is the depth where current density falls to 1/e (~37%)
    - Lower frequency → deeper penetration
    
    Typical values:
    - f = 2 Hz (end of deployment): δ ≈ 33 mm
    - f = 0.5 Hz (early deployment): δ ≈ 65 mm
    
    If δ > t_plate, the reaction plate is "magnetically thin"
    and eddy currents use the full thickness.
    """
    if f_slip <= 0:
        return 0.1  # Large value for DC (full penetration)
    return math.sqrt(rho / (math.pi * MU0 * f_slip))


# ==============================================================================
# GOODNESS FACTOR
# ==============================================================================
def calc_goodness_factor(p: LIMParams, f_slip: float) -> float:
    """
    Goodness Factor G (dimensionless)
    
    Formula:
        G = (ω_slip × μ₀ × σ × δ_eff × τp) / π
    
    Physical meaning:
    - G measures how well the LIM converts electrical to mechanical power
    - Large G (>100) means efficient energy transfer
    - G depends on frequency, conductivity, and geometry
    
    For orbital ring LIM:
    - G ≈ 500-1000 (very large due to long pole pitch)
    - This means optimal slip is very small: s_opt ≈ 1/G ≈ 0.1-0.2%
    """
    if f_slip <= 0:
        return 0
    
    omega_slip = 2 * math.pi * f_slip
    delta = calc_skin_depth(f_slip)
    delta_eff = min(delta, p.t_plate)
    
    return (omega_slip * MU0 * SIGMA_AL_106K * delta_eff * p.tau_p) / math.pi


# ==============================================================================
# SLIP EFFICIENCY
# ==============================================================================
def calc_slip_efficiency(slip: float, G: float) -> float:
    """
    Slip Efficiency η_slip (dimensionless, 0 to 1)
    
    Formula:
        η_slip = 2sG / (1 + s²G²)
    
    Physical meaning:
    - Fraction of air-gap power converted to mechanical thrust
    - Maximum η = 1 occurs at s = 1/G
    - Our design typically operates below optimal slip
    
    Example with G = 1000:
    - Optimal slip: s = 1/1000 = 0.1%
    - At s = 2.5% (typical): η ≈ 8%
    
    Why we don't run at optimal slip:
    - Optimal slip (0.1%) is hard to control precisely
    - Higher slip gives more stable operation
    - We accept lower efficiency for controllability
    """
    if G <= 0:
        return 0
    sG = slip * G
    return 2 * sG / (1 + sG**2)


# ==============================================================================
# THRUST
# ==============================================================================
def calc_thrust(p: LIMParams, v_rel: float) -> float:
    """
    Thrust Force (N)
    
    Formula:
        F = (B²/2μ₀) × A_active × η_slip
    
    Components:
    - B²/2μ₀: Magnetic pressure (Pa) - force per unit area from B-field
    - A_active: Area where field acts on reaction plate (m²)
    - η_slip: Efficiency of converting magnetic pressure to thrust
    
    Typical values:
    - B = 40 mT, A = 300 m², η = 8%
    - Magnetic pressure = (0.04)²/(2×1.26e-6) = 637 Pa
    - Raw force = 637 × 300 = 191 kN
    - Actual thrust = 191 × 0.08 = 15 kN per LIM
    
    With ~83,000 LIM sites × 2 LIMs each:
    - Total thrust ≈ 2.5 MN (250,000 tonnes-force equivalent)
    """
    B = calc_B_field(p)
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    
    magnetic_pressure = B**2 / (2 * MU0)
    return magnetic_pressure * p.A_active * eta


# ==============================================================================
# COMPREHENSIVE EVALUATION
# ==============================================================================
def evaluate_design(p: LIMParams, v_rel: float) -> Dict:
    """
    Evaluate all aspects of a LIM design at a given operating point.
    
    Returns a dictionary with:
    - All calculated values
    - Constraint margins (positive = OK, negative = violated)
    - Feasibility flag
    """
    # Calculate all values
    B = calc_B_field(p)
    L = calc_inductance(p)
    V = calc_voltage(p, v_rel)
    F = calc_thrust(p, v_rel)
    P = F * v_rel  # Mechanical power
    
    freq = calc_frequencies(p, v_rel)
    G = calc_goodness_factor(p, freq['f_slip'])
    eta = calc_slip_efficiency(freq['slip'], G)
    delta = calc_skin_depth(freq['f_slip'])
    
    # Check constraints
    V_margin = (p.V_kapton_limit - V) / p.V_kapton_limit * 100
    I_margin = (p.I_c - p.I_peak) / p.I_c * 100
    
    return {
        # Magnetic
        'B_mT': B * 1000,
        'magnetic_pressure_Pa': B**2 / (2 * MU0),
        
        # Electrical
        'L_mH': L * 1000,
        'V_kV': V / 1000,
        'V_limit_kV': p.V_kapton_limit / 1000,
        'V_margin_%': V_margin,
        'I_margin_%': I_margin,
        
        # Frequency
        'f_slip_Hz': freq['f_slip'],
        'f_supply_Hz': freq['f_supply'],
        'slip_%': freq['slip'] * 100,
        
        # Efficiency
        'G': G,
        'eta_%': eta * 100,
        'delta_mm': delta * 1000,
        'delta_eff_mm': min(delta, p.t_plate) * 1000,
        
        # Performance
        'F_kN': F / 1000,
        'P_MW': P / 1e6,
        
        # Design info
        'k_fill_%': p.k_fill * 100,
        
        # Status
        'feasible': V_margin > 0 and I_margin > 0
    }


# ==============================================================================
# OPTIMIZATION FUNCTIONS
# ==============================================================================
def find_optimal_N(p: LIMParams, v_rel: float, 
                   N_range: range = range(50, 501, 10)) -> Tuple[int, Dict]:
    """
    Find the turn count that maximizes thrust while staying within voltage limit.
    
    Strategy: Sweep N from low to high, keep the highest N that's still feasible.
    More turns = more thrust, but also more voltage.
    """
    best_N = N_range[0]
    best_thrust = 0
    best_result = None
    
    original_N = p.N
    
    for N in N_range:
        p.N = N
        result = evaluate_design(p, v_rel)
        
        if result['feasible'] and result['F_kN'] > best_thrust:
            best_thrust = result['F_kN']
            best_N = N
            best_result = result
    
    p.N = original_N  # Restore
    return best_N, best_result


def parameter_sweep(base_params: LIMParams, v_rel: float,
                    param_name: str, values: List[float]) -> List[Dict]:
    """
    Sweep a single parameter and return results for each value.
    Useful for generating tables for the book.
    """
    results = []
    
    for val in values:
        # Create a copy with modified parameter
        p = LIMParams(**{k: v for k, v in vars(base_params).items()})
        setattr(p, param_name, val)
        
        result = evaluate_design(p, v_rel)
        result[param_name] = val
        results.append(result)
    
    return results


# ==============================================================================
# OUTPUT FORMATTING
# ==============================================================================
def print_design_report(p: LIMParams, v_rel: float):
    """
    Print a comprehensive, book-ready design report.
    """
    result = evaluate_design(p, v_rel)
    
    print()
    print("=" * 70)
    print("LINEAR INDUCTION MOTOR DESIGN REPORT")
    print("=" * 70)
    
    print("\n┌─ DESIGN PARAMETERS ─────────────────────────────────────────────────")
    print(f"│ Turns per coil (N):           {p.N:>10}")
    print(f"│ Coil width (W):               {p.W:>10.2f} m")
    print(f"│ Pole pitch (τp):              {p.tau_p:>10.1f} m")
    print(f"│ Pole pitches per LIM:         {p.pitch_count:>10}")
    print(f"│ Air gap (g):                  {p.gap*100:>10.1f} cm")
    print(f"│ Peak current (I):             {p.I_peak:>10.0f} A")
    print(f"│ Slip velocity (v_slip):       {p.v_slip:>10.0f} m/s")
    print(f"│ Relative velocity (v_rel):    {v_rel:>10.0f} m/s")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ TAPE CONSTRUCTION ──────────────────────────────────────────────────")
    print(f"│ HTS tape thickness:           {p.d_HTS*1e6:>10.0f} μm")
    print(f"│ Kapton thickness:             {p.d_kapton*1000:>10.2f} mm")
    print(f"│ Fill factor (k_fill):         {result['k_fill_%']:>10.1f} %")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ MAGNETIC FIELD ─────────────────────────────────────────────────────")
    print(f"│ B at reaction plate:          {result['B_mT']:>10.1f} mT")
    print(f"│ Magnetic pressure:            {result['magnetic_pressure_Pa']:>10.0f} Pa")
    print(f"│ Active area:                  {p.A_active:>10.0f} m²")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ ELECTRICAL ─────────────────────────────────────────────────────────")
    print(f"│ Inductance (L):               {result['L_mH']:>10.1f} mH")
    print(f"│ Supply frequency:             {result['f_supply_Hz']:>10.1f} Hz")
    print(f"│ Coil voltage:                 {result['V_kV']:>10.1f} kV")
    print(f"│ Voltage limit (Kapton):       {result['V_limit_kV']:>10.0f} kV")
    print(f"│ Voltage margin:               {result['V_margin_%']:>10.1f} %")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ SLIP AND EFFICIENCY ────────────────────────────────────────────────")
    print(f"│ Slip frequency:               {result['f_slip_Hz']:>10.2f} Hz")
    print(f"│ Slip ratio:                   {result['slip_%']:>10.2f} %")
    print(f"│ Goodness factor (G):          {result['G']:>10.0f}")
    print(f"│ Optimal slip (1/G):           {100/result['G'] if result['G'] > 0 else 0:>10.2f} %")
    print(f"│ Slip efficiency (η):          {result['eta_%']:>10.1f} %")
    print(f"│ Skin depth:                   {result['delta_mm']:>10.1f} mm")
    print("└──────────────────────────────────────────────────────────────────────")
    
    print("\n┌─ PERFORMANCE ────────────────────────────────────────────────────────")
    print(f"│ Thrust per LIM:               {result['F_kN']:>10.1f} kN")
    print(f"│ Mechanical power per LIM:     {result['P_MW']:>10.2f} MW")
    print("└──────────────────────────────────────────────────────────────────────")
    
    status = "✓ FEASIBLE" if result['feasible'] else "✗ VOLTAGE LIMIT EXCEEDED"
    print(f"\n  Status: {status}")
    print("=" * 70)


def print_sweep_table(results: List[Dict], sweep_param: str, columns: List[str]):
    """
    Print a formatted table of parameter sweep results.
    """
    # Header
    print()
    header = " │ ".join(f"{c:>12}" for c in columns)
    print("─" * len(header))
    print(header)
    print("─" * len(header))
    
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
    
    print("─" * len(header))


# ==============================================================================
# MAIN: GENERATE OPTIMIZATION GUIDE
# ==============================================================================
def main():
    """
    Generate a comprehensive optimization guide for the book.
    """
    print("\n" + "=" * 70)
    print("LIM PARAMETER OPTIMIZATION GUIDE")
    print("For: Ion Propulsion Engineering")
    print("=" * 70)
    
    # Base parameters matching legacy_power_lim.py
    p = LIMParams(
        N=100,
        W=2.0,
        tau_p=50.0,
        gap=0.05,          # 5 cm
        I_peak=650,
        v_slip=200,
        d_kapton=1e-3,     # 1 mm
        E_kapton_safe=100e6,  # 100 kV/mm
    )
    
    v_rel_end = 8000  # End of deployment
    
    # Show base design
    print("\n" + "=" * 70)
    print("BASE DESIGN (N=100, gap=5cm, v_rel=8000 m/s)")
    print("=" * 70)
    print_design_report(p, v_rel_end)
    
    # ==== Turn Count Sweep ====
    print("\n" + "=" * 70)
    print("1. EFFECT OF TURN COUNT (N)")
    print("=" * 70)
    print("""
More turns creates a stronger magnetic field (B ∝ N), which increases
thrust (F ∝ B²). However, inductance also increases (L ∝ N²), and
voltage increases with it (V ∝ L). The Kapton insulation limits voltage.

Key insight: Thrust scales as N², voltage scales as N². So the thrust
you can achieve before hitting the voltage limit is set by geometry
and operating frequency, not by how many turns you wind.
""")
    
    N_values = [50, 100, 150, 200, 250, 300, 350, 400]
    results = parameter_sweep(p, v_rel_end, 'N', N_values)
    print(f"At v_rel = {v_rel_end} m/s (end of deployment):")
    print_sweep_table(results, 'N', ['N', 'B_mT', 'L_mH', 'V_kV', 'V_limit_kV', 'V_margin_%', 'F_kN', 'feasible'])
    
    opt_N, opt_result = find_optimal_N(p, v_rel_end)
    print(f"\n→ Optimal N = {opt_N} (maximum thrust within voltage limit)")
    print(f"  Thrust = {opt_result['F_kN']:.1f} kN, Voltage = {opt_result['V_kV']:.1f} kV")
    
    # ==== Gap Sweep ====
    print("\n" + "=" * 70)
    print("2. EFFECT OF AIR GAP")
    print("=" * 70)
    print("""
Smaller gaps increase the magnetic field at the reaction plate.
With the CORRECTED formula, B saturates as gap → 0 instead of blowing up.

The practical limit on gap is mechanical: maintaining 2 cm clearance
between a moving cable and stationary coils over 40,000 km is challenging.
""")
    
    gap_values = [0.02, 0.03, 0.05, 0.10, 0.15, 0.20]
    results = parameter_sweep(p, v_rel_end, 'gap', gap_values)
    for r in results:
        r['gap_cm'] = r['gap'] * 100
    print(f"At N = {p.N}, v_rel = {v_rel_end} m/s:")
    print_sweep_table(results, 'gap_cm', ['gap_cm', 'B_mT', 'L_mH', 'V_kV', 'F_kN'])
    
    # ==== Slip Velocity Sweep ====
    print("\n" + "=" * 70)
    print("3. EFFECT OF SLIP VELOCITY")
    print("=" * 70)
    print("""
Slip velocity determines how much slower the cable moves than the
magnetic wave. Higher slip → higher slip frequency → more eddy losses.

The goodness factor G is large (~1000), so optimal slip is tiny (~0.1%).
We run at higher slip (2-3%) for controllability, accepting lower efficiency.

SURPRISE: Lower slip velocity gives MORE thrust! This is because
we're closer to optimal slip where efficiency peaks.
""")
    
    v_slip_values = [50, 100, 150, 200, 250, 300, 400]
    results = parameter_sweep(p, v_rel_end, 'v_slip', v_slip_values)
    print(f"At N = {p.N}, gap = {p.gap*100:.0f} cm, v_rel = {v_rel_end} m/s:")
    print_sweep_table(results, 'v_slip', ['v_slip', 'slip_%', 'G', 'eta_%', 'F_kN', 'V_kV'])
    
    # ==== Deployment Progression ====
    print("\n" + "=" * 70)
    print("4. DEPLOYMENT PROGRESSION")
    print("=" * 70)
    print("""
During deployment, v_rel increases from ~100 m/s to ~8000 m/s.
This means supply frequency increases, and with it, voltage.

The design must work at the END of deployment when voltage is highest.
At the START, voltage is low but so is thrust (due to low efficiency).
""")
    
    v_rel_values = [100, 500, 1000, 2000, 4000, 6000, 8000]
    results = []
    for v in v_rel_values:
        result = evaluate_design(p, v)
        result['v_rel'] = v
        results.append(result)
    print(f"At N = {p.N}, gap = {p.gap*100:.0f} cm, v_slip = {p.v_slip} m/s:")
    print_sweep_table(results, 'v_rel', ['v_rel', 'f_supply_Hz', 'V_kV', 'slip_%', 'eta_%', 'F_kN', 'P_MW'])
    
    # ==== Kapton Thickness ====
    print("\n" + "=" * 70)
    print("5. KAPTON THICKNESS TRADE-OFF")
    print("=" * 70)
    print("""
Kapton insulation thickness has two competing effects:

1. Thicker Kapton → Higher voltage limit (V_limit = E_safe × d_kapton)
2. Thicker Kapton → Lower fill factor → Lower inductance → Lower voltage

At 100 kV/mm safe field strength:
- 0.5 mm Kapton → 50 kV limit
- 1.0 mm Kapton → 100 kV limit
- 2.0 mm Kapton → 200 kV limit

But thicker Kapton also means less superconductor per unit volume.
""")
    
    d_kapton_values = [0.25e-3, 0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3]
    results = []
    for d in d_kapton_values:
        p_test = LIMParams(**{k: v for k, v in vars(p).items()})
        p_test.d_kapton = d
        result = evaluate_design(p_test, v_rel_end)
        result['d_kapton_mm'] = d * 1000
        results.append(result)
    
    print(f"At N = {p.N}, v_rel = {v_rel_end} m/s:")
    print_sweep_table(results, 'd_kapton_mm', ['d_kapton_mm', 'k_fill_%', 'V_limit_kV', 'L_mH', 'V_kV', 'V_margin_%'])
    
    # ==== Summary ====
    print("\n" + "=" * 70)
    print("OPTIMIZATION SUMMARY")
    print("=" * 70)
    print("""
VOLTAGE LIMIT (most important constraint):
- Total coil voltage must stay below Kapton breakdown
- V_limit = E_safe × d_kapton (e.g., 100 kV/mm × 1mm = 100 kV)
- This is the TOTAL COIL voltage, not per-layer
- Design for worst case: end of deployment (v_rel = 8000 m/s)

TO MAXIMIZE THRUST:
1. Increase turns (N) until voltage limit is reached
2. Use smallest practical gap (2-5 cm)
3. Use lower slip velocity (more efficient, but harder to control)

RECOMMENDED DESIGN PROCESS:
1. Fix Kapton thickness based on voltage limit needed
2. Fix gap based on mechanical constraints
3. Choose slip velocity for controllability (~200 m/s)
4. Find maximum N that keeps V < V_limit at end of deployment
5. Verify thermal margins at start of deployment

With base parameters (gap=5cm, d_kapton=1mm, v_slip=200 m/s):
- Voltage limit: 100 kV
- Optimal turns: N ≈ 300
- Maximum thrust: ~140 kN per LIM
""")


if __name__ == "__main__":
    main()
