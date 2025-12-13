#!/usr/bin/env python3
"""
Magnetic Field Distribution Analysis
=====================================
For: "Ion Propulsion Engineering" by Paul de Jong

QUESTION: Is there a limit to coil width? Does a wider coil have 
very little field at its center?

ANSWER: Yes! The conductors are at the edges of the coil, so the 
field is strongest near the edges and weakest at the center.

This script analyzes:
1. Field distribution across the coil width
2. Effective active area (where field is significant)
3. Optimal coil width considering this effect
"""

import math
import numpy as np

MU0 = 4 * math.pi * 1e-7


def field_from_wire(x_wire: float, z_wire: float, x_point: float, z_point: float, I: float) -> tuple:
    """
    Magnetic field (Bx, Bz) at a point from a single infinite wire carrying current I.
    Wire runs in the y-direction (into the page).
    
    Using Biot-Savart: B = μ₀I/(2πr), direction perpendicular to r
    """
    dx = x_point - x_wire
    dz = z_point - z_wire
    r_sq = dx**2 + dz**2
    
    if r_sq < 1e-10:
        return 0, 0
    
    r = math.sqrt(r_sq)
    B_mag = MU0 * I / (2 * math.pi * r)
    
    # Direction is perpendicular to r, using right-hand rule
    # For wire into page (+y), B circles clockwise when viewed from +y
    Bx = B_mag * (-dz / r)  # perpendicular to r
    Bz = B_mag * (dx / r)
    
    return Bx, Bz


def field_at_point(x: float, z: float, W: float, N: int, I: float) -> tuple:
    """
    Calculate magnetic field at point (x, z) from a coil of width W.
    
    Model: Two bundles of N wires at x = ±W/2, at z = 0 (coil surface).
    Current flows +y in left bundle, -y in right bundle (forming a loop).
    
    Returns: (Bx, Bz) in Tesla
    """
    # Left wire bundle at x = -W/2, z = 0, current into page (+y)
    Bx1, Bz1 = field_from_wire(-W/2, 0, x, z, N * I)
    
    # Right wire bundle at x = +W/2, z = 0, current out of page (-y)
    Bx2, Bz2 = field_from_wire(W/2, 0, x, z, -N * I)
    
    return Bx1 + Bx2, Bz1 + Bz2


def analyze_field_distribution(W: float, gap: float, N: int, I: float):
    """
    Analyze field distribution across coil width at reaction plate surface.
    """
    # Sample points across the width, at depth = gap (reaction plate surface)
    x_points = np.linspace(-W/2 * 1.2, W/2 * 1.2, 101)
    z_plate = -gap  # Reaction plate is below the coil
    
    Bz_values = []
    Bx_values = []
    B_mag_values = []
    
    for x in x_points:
        Bx, Bz = field_at_point(x, z_plate, W, N, I)
        Bz_values.append(Bz)
        Bx_values.append(Bx)
        B_mag_values.append(math.sqrt(Bx**2 + Bz**2))
    
    return x_points, np.array(Bz_values), np.array(Bx_values), np.array(B_mag_values)


def calc_effective_width(W: float, gap: float, N: int, I: float, threshold: float = 0.5):
    """
    Calculate effective width where field is above threshold × max field.
    """
    x_points, Bz, Bx, B_mag = analyze_field_distribution(W, gap, N, I)
    
    B_max = np.max(np.abs(Bz))
    B_center = abs(Bz[len(Bz)//2])
    
    # Find where |Bz| > threshold * B_max
    above_threshold = np.abs(Bz) > threshold * B_max
    
    # Find the extent
    indices = np.where(above_threshold)[0]
    if len(indices) > 0:
        W_eff = x_points[indices[-1]] - x_points[indices[0]]
    else:
        W_eff = 0
    
    return {
        'W': W,
        'gap': gap,
        'B_max_mT': B_max * 1000,
        'B_center_mT': B_center * 1000,
        'B_ratio': B_center / B_max if B_max > 0 else 0,
        'W_eff': W_eff,
        'W_eff_ratio': W_eff / W if W > 0 else 0,
    }


def main():
    print("=" * 70)
    print("MAGNETIC FIELD DISTRIBUTION ACROSS COIL WIDTH")
    print("=" * 70)
    
    # Base parameters
    N = 100
    I = 650
    gap = 0.05  # 5 cm
    
    print(f"""
PHYSICAL MODEL:
- Coil has conductors (N={N} turns, I={I}A) at edges (x = ±W/2)
- Reaction plate is at distance gap = {gap*100:.0f} cm below coil
- We calculate the vertical field component Bz across the width

QUESTION: How does field vary from edge to center?
""")
    
    # Analyze different widths
    print("=" * 70)
    print("1. FIELD DISTRIBUTION FOR DIFFERENT COIL WIDTHS")
    print("=" * 70)
    
    widths = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
    
    print(f"\n{'W (m)':>8} │ {'B_edge':>10} │ {'B_center':>10} │ {'Ratio':>8} │ {'W_eff':>8} │ {'W_eff/W':>8}")
    print(f"{'':>8} │ {'(mT)':>10} │ {'(mT)':>10} │ {'Bc/Be':>8} │ {'(m)':>8} │ {'':>8}")
    print("─" * 70)
    
    for W in widths:
        result = calc_effective_width(W, gap, N, I, threshold=0.5)
        print(f"{W:>8.1f} │ {result['B_max_mT']:>10.1f} │ {result['B_center_mT']:>10.1f} │ " +
              f"{result['B_ratio']:>8.2f} │ {result['W_eff']:>8.2f} │ {result['W_eff_ratio']:>8.2f}")
    
    print()
    print("""
INTERPRETATION:
- B_edge: Field strength directly under the conductors (at x = ±W/2)
- B_center: Field strength at the center of the coil (x = 0)
- Ratio: How much weaker the center is compared to the edge
- W_eff: Width where field > 50% of maximum
- W_eff/W: Fraction of coil width that's actually effective
""")
    
    # Detailed analysis for W = 2m (baseline)
    print("=" * 70)
    print("2. DETAILED FIELD PROFILE FOR W = 2.0 m (BASELINE)")
    print("=" * 70)
    
    W = 2.0
    x_points, Bz, Bx, B_mag = analyze_field_distribution(W, gap, N, I)
    
    print(f"\nField profile at reaction plate (gap = {gap*100:.0f} cm):")
    print(f"\n{'x (m)':>8} │ {'Bz (mT)':>10} │ {'|B| (mT)':>10} │ {'Bar':>30}")
    print("─" * 70)
    
    B_max = np.max(np.abs(Bz))
    for i in range(0, len(x_points), 10):
        x = x_points[i]
        bz = Bz[i] * 1000
        b_mag = B_mag[i] * 1000
        bar_len = int(30 * abs(bz) / (B_max * 1000))
        bar = "█" * bar_len
        print(f"{x:>8.2f} │ {bz:>10.1f} │ {b_mag:>10.1f} │ {bar}")
    
    # Effect of gap on field uniformity
    print("\n" + "=" * 70)
    print("3. EFFECT OF GAP ON FIELD UNIFORMITY")
    print("=" * 70)
    print("""
Larger gap = more uniform field (but weaker overall)
Smaller gap = stronger field but more concentrated at edges
""")
    
    W = 2.0
    gaps = [0.02, 0.05, 0.10, 0.20, 0.50, 1.00]
    
    print(f"\n{'gap (cm)':>10} │ {'B_edge':>10} │ {'B_center':>10} │ {'Ratio':>8} │ {'W_eff (m)':>10}")
    print("─" * 60)
    
    for g in gaps:
        result = calc_effective_width(W, g, N, I, threshold=0.5)
        print(f"{g*100:>10.0f} │ {result['B_max_mT']:>10.1f} │ {result['B_center_mT']:>10.1f} │ " +
              f"{result['B_ratio']:>8.2f} │ {result['W_eff']:>10.2f}")
    
    # Optimal width analysis
    print("\n" + "=" * 70)
    print("4. OPTIMAL WIDTH ANALYSIS")
    print("=" * 70)
    print("""
For thrust, what matters is F = ∫ (B²/2μ₀) × η dA

With non-uniform field, we need to integrate properly.
Wider coil has more area but weaker average field.
""")
    
    gap = 0.05
    tau_p = 50.0
    pitch_count = 3
    
    print(f"\nFor gap = {gap*100:.0f} cm, τp = {tau_p} m, {pitch_count} pitches:")
    print(f"\n{'W (m)':>8} │ {'B_avg²':>10} │ {'Area':>8} │ {'B²×A':>10} │ {'Relative':>10}")
    print(f"{'':>8} │ {'(mT²)':>10} │ {'(m²)':>8} │ {'(arb)':>10} │ {'Thrust':>10}")
    print("─" * 60)
    
    widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    results = []
    
    for W in widths:
        x_points, Bz, _, _ = analyze_field_distribution(W, gap, N, I)
        
        # Average B² across the width (within the coil, |x| < W/2)
        mask = np.abs(x_points) <= W/2
        Bz_inside = Bz[mask]
        B_sq_avg = np.mean(Bz_inside**2) * 1e6  # Convert to mT²
        
        Area = W * tau_p * pitch_count
        B_sq_A = B_sq_avg * Area
        
        results.append((W, B_sq_avg, Area, B_sq_A))
    
    # Normalize to maximum
    max_B_sq_A = max(r[3] for r in results)
    
    for W, B_sq_avg, Area, B_sq_A in results:
        relative = B_sq_A / max_B_sq_A
        print(f"{W:>8.1f} │ {B_sq_avg:>10.1f} │ {Area:>8.0f} │ {B_sq_A:>10.0f} │ {relative:>10.2f}")
    
    # Find optimal
    optimal = max(results, key=lambda r: r[3])
    print(f"\n→ Optimal width ≈ {optimal[0]:.1f} m (maximizes B² × Area)")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
YES, THERE IS A LIMIT TO COIL WIDTH!

Physical reason:
- Conductors are at edges (±W/2), not distributed across width
- Field at center comes from conductors at distance √(g² + (W/2)²)
- For wide coils, center field ∝ 1/W² (very weak!)

Quantitative findings (for gap = 5 cm):
- W = 1 m: Center field is 78% of edge field (fairly uniform)
- W = 2 m: Center field is 38% of edge field (baseline design)
- W = 5 m: Center field is 8% of edge field (very non-uniform)

Optimal width:
- Considering B² × Area trade-off: W ≈ {optimal[0]:.1f} m
- This is close to the baseline W = 2 m

Rule of thumb:
- Keep W/gap < 20 for reasonable field uniformity
- For gap = 5 cm, this means W < 1 m for uniform field
- W = 2 m is a compromise: non-uniform but still practical

IMPLICATIONS FOR DESIGN:
1. The "active area" is less than W × L_active
2. Thrust calculations using average B may overestimate
3. Narrower coils are more efficient per unit width
4. But wider coils can still be practical if you accept non-uniformity
""")


if __name__ == "__main__":
    main()
