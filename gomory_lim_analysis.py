#!/usr/bin/env python3
"""
Analysis: Connecting Gömöry's HTS Research to Orbital Ring LIM Coil Design
===========================================================================

This script works backwards from experimental HTS coil parameters found in
Prof. Dr. Fedor Gömöry's research to derive appropriate values for:
- Fill factor (k_fill)
- Voltage limits
- Inductance scaling

For the Ion Propulsion Engineering textbook.

Updated: January 2026 to match current lim_config.py defaults
"""

import numpy as np

# =============================================================================
# SECTION 1: TYPICAL HTS TAPE PARAMETERS (from Gömöry's work and suppliers)
# =============================================================================

print("=" * 70)
print("SECTION 1: HTS TAPE LAYER STRUCTURE")
print("=" * 70)

# Typical 2G REBCO coated conductor structure (from SuperPower, AMSC, etc.)
# All dimensions in micrometers (μm)

tape_layers = {
    "Substrate (Hastelloy C-276)": 50.0,      # structural support
    "Buffer layers (MgO, etc.)": 0.5,          # texture template
    "REBCO superconductor": 2.0,               # the actual superconductor
    "Silver cap layer": 2.0,                   # current sharing
    "Copper stabilizer (each side)": 20.0,     # thermal/electrical stabilizer
}

total_tape_thickness = sum(tape_layers.values()) + 20.0  # +20 for other Cu side
print(f"\nHTS Tape Layer Stack:")
for layer, thickness in tape_layers.items():
    print(f"  {layer:35s}: {thickness:6.1f} μm")
print(f"  {'TOTAL TAPE THICKNESS':35s}: {total_tape_thickness:6.1f} μm")

# What fraction is superconductor?
sc_fraction = tape_layers["REBCO superconductor"] / total_tape_thickness
print(f"\nSuperconductor fraction of tape: {sc_fraction:.1%}")

# Typical tape dimensions
tape_width_mm = 12.0  # mm (common width for power applications)
tape_thickness_mm = total_tape_thickness / 1000  # convert to mm

print(f"\nTypical tape: {tape_width_mm} mm × {tape_thickness_mm:.3f} mm")

# =============================================================================
# SECTION 2: INSULATION OPTIONS AND FILL FACTOR
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 2: INSULATION AND FILL FACTOR CALCULATION")
print("=" * 70)

# Kapton co-wind insulation (used by BNL, SuperPower, and in Gömöry's experiments)
kapton_options = {
    "Thin Kapton (25 μm)": 25.0,
    "Standard Kapton (38 μm)": 38.0,
    "Thick Kapton (50 μm)": 50.0,
}

# Additional adhesive layer (often required)
adhesive_thickness = 60.0  # μm (can range 50-100 μm)

print("\nFill Factor Calculation for Different Insulation Options:")
print("-" * 70)
print(f"{'Insulation Type':30s} {'Kapton':>10s} {'Adhesive':>10s} {'Total Turn':>12s} {'k_fill':>10s}")
print("-" * 70)

for name, kapton_um in kapton_options.items():
    total_insulation = kapton_um + adhesive_thickness
    total_turn_thickness = total_tape_thickness + total_insulation
    k_fill = total_tape_thickness / total_turn_thickness
    print(f"{name:30s} {kapton_um:>10.0f} μm {adhesive_thickness:>10.0f} μm {total_turn_thickness:>12.0f} μm {k_fill:>10.3f}")

# The "65% conductor packing factor" mentioned in IEEE papers
print("\nNote: Literature reports '65% conductor packing factor' for optimized 2G HTS coils")
print("This aligns with k_fill ≈ 0.5 to 0.65 depending on insulation thickness.")

# RECOMMENDED VALUE
k_fill_recommended = 0.50
print(f"\n>>> RECOMMENDED k_fill = {k_fill_recommended} (conservative, uses standard Kapton + adhesive)")

# =============================================================================
# SECTION 3: INDUCTANCE ANALYSIS - GÖMÖRY'S EXPERIMENTAL COILS
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 3: INDUCTANCE FROM GÖMÖRY'S EXPERIMENTAL COILS")
print("=" * 70)

# Examples from literature search
gomory_coils = [
    {"name": "10-turn pancake (SuperPower 12mm)", "N": 10, "ID_mm": 60, "OD_mm": 80, "L_mH": None},
    {"name": "62-turn split coil", "N": 62, "ID_mm": 50, "OD_mm": 80, "L_mH": None},
    {"name": "Double pancake coil (BiSCCO)", "N": 364, "ID_mm": 108, "OD_mm": 220, "L_mH": 28},
    {"name": "DPC with iron core", "N": None, "ID_mm": None, "OD_mm": None, "L_mH": 45},
]

# For a pancake coil, approximate inductance:
# L ≈ μ₀ × N² × A_mean / (2π × r_mean) × geometry_factor

mu_0 = 4 * np.pi * 1e-7  # H/m

print("\nAnalysis of experimental coil inductances:")
print("-" * 70)

# Let's estimate what the 10-turn coil would have
# For a single-layer pancake coil with N turns:
# L ≈ μ₀ × N² × r_mean × [ln(8×r_mean/a) - 2]
# where a is the wire radius

for coil in gomory_coils:
    if coil["ID_mm"] is not None and coil["OD_mm"] is not None and coil["N"] is not None:
        r_inner = coil["ID_mm"] / 2000  # m
        r_outer = coil["OD_mm"] / 2000  # m
        r_mean = (r_inner + r_outer) / 2
        N = coil["N"]
        
        # Simple approximation for pancake coil inductance
        # Using Neumann formula approximation
        radial_width = r_outer - r_inner
        height = tape_thickness_mm / 1000 * N  # approximate stack height
        
        # Wheeler's formula approximation for multilayer coil
        # L ≈ 31.6 × N² × r_mean² / (6×r_mean + 9×height + 10×radial_width) [μH, dimensions in inches]
        # Convert to SI
        r_mean_inch = r_mean * 39.37
        height_inch = height * 39.37
        width_inch = radial_width * 39.37
        
        L_uH = 31.6 * N**2 * r_mean_inch**2 / (6*r_mean_inch + 9*height_inch + 10*width_inch)
        L_mH = L_uH / 1000
        
        print(f"\n{coil['name']}:")
        print(f"  N = {N}, r_mean = {r_mean*1000:.0f} mm")
        if coil["L_mH"]:
            print(f"  Measured L = {coil['L_mH']} mH")
        print(f"  Estimated L ≈ {L_mH:.2f} mH")

# =============================================================================
# SECTION 4: VOLTAGE LIMITS FOR HTS APPLICATIONS
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 4: VOLTAGE LIMITS FOR HTS COILS")
print("=" * 70)

voltage_applications = """
Application                          Typical Voltage      Notes
---------------------------------------------------------------------------
HTS Research pancake coils           10-100 V            Low voltage, lab testing
HTS Motor/Generator field coils      100-1000 V          Mid-range, insulated
HTS Power cables (Japan M-PACC)      66 kV, 275 kV       Specialized PPLP insulation
Conventional motor hi-pot test       2 × U_rated + 1 kV  IEC/IEEE standards
Kapton breakdown (per mil)           ~7 kV/mil           ~275 kV/mm
"""
print(voltage_applications)

# Kapton dielectric strength
kapton_breakdown_kV_per_mil = 7.0  # kV per mil (25.4 μm)
kapton_breakdown_kV_per_mm = kapton_breakdown_kV_per_mil / 0.0254

# For 38 μm Kapton (1.5 mil)
kapton_38um_breakdown = 1.5 * kapton_breakdown_kV_per_mil

print(f"Kapton dielectric strength: {kapton_breakdown_kV_per_mm:.0f} kV/mm")
print(f"38 μm Kapton can withstand: {kapton_38um_breakdown:.1f} kV (before breakdown)")
print(f"With safety factor of 3: {kapton_38um_breakdown/3:.1f} kV working voltage per turn")

# For our LIM application
print("\n>>> RECOMMENDED VOLTAGE LIMIT:")
print("    For Kapton-insulated HTS LIM coils: 2 kV is conservative and appropriate")
print("    This allows for:")
print("    - Safety factor of ~5 on Kapton breakdown")
print("    - Margin for high-frequency effects (skin effect in Al)")
print("    - Conservative approach for novel space application")

# =============================================================================
# SECTION 5: APPLYING TO ORBITAL RING LIM COILS
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 5: ORBITAL RING LIM COIL PARAMETERS")
print("=" * 70)

# LIM parameters from lim_config.py (January 2026)
tau_p = 100.0       # pole pitch [m]
w_coil = 2.0        # coil width [m]
gap = 0.20          # magnetic gap [m]
I_peak = 163        # peak current [A] - 81.25% of I_c for 3mm × 1 layer
N_turns = 200       # turns per coil

# Calculate supply frequency at typical operating point
v_rel_typical = 3000  # m/s (mid-deployment)
v_slip = 60           # m/s (2% slip)
v_wave = v_rel_typical + v_slip
f_supply = v_wave / (2 * tau_p)  # ~15 Hz at mid-deployment

print(f"\nLIM Operating Parameters (from lim_config.py):")
print(f"  Pole pitch (τ_p):     {tau_p} m")
print(f"  Coil width:           {w_coil} m")
print(f"  Magnetic gap:         {gap} m")
print(f"  Peak current:         {I_peak} A (3mm × 1 layer HTS)")
print(f"  Number of turns:      {N_turns}")
print(f"  Supply frequency:     {f_supply:.1f} Hz (at v_rel = {v_rel_typical} m/s)")

# Coil inductance formula (from legacy code analysis)
# L_coil = N² × μ₀ × (A_coil/gap) × k_fill
# where A_coil is the coil cross-sectional area

print("\n" + "-" * 70)
print("VOLTAGE ANALYSIS FOR DIFFERENT TURN COUNTS AND FILL FACTORS")
print("-" * 70)

# The phase voltage equation
# V_phase = (2π × f × L_coil × I_peak) / √3

def calculate_voltage(N, k_fill, tau_p, w_coil, gap, f, I_peak):
    """Calculate phase voltage for given parameters."""
    # Coil area - assuming single-layer flat coil spanning one pole
    # A_coil ≈ tau_p × w_coil (face area)
    # But for inductance we need the flux linkage area
    # This is a simplified model
    A_coil = tau_p * w_coil  # m²
    
    L_coil = N**2 * mu_0 * (A_coil / gap) * k_fill
    
    # Phase voltage (assuming star connection)
    V_phase = (2 * np.pi * f * L_coil * I_peak) / np.sqrt(3)
    
    return L_coil, V_phase

print(f"\nAt f_supply = {f_supply:.1f} Hz (mid-deployment):")
print(f"\n{'k_fill':>8s} {'N':>6s} {'L_coil (mH)':>12s} {'V_phase (kV)':>12s} {'Status':>15s}")
print("-" * 55)

V_limit = 100000  # V (100 kV max from lim_config.py)

for k_fill in [0.05, 0.10, 0.20, 0.50]:
    for N in [100, 150, 200, 250, 300]:
        L_coil, V_phase = calculate_voltage(N, k_fill, tau_p, w_coil, gap, f_supply, I_peak)
        
        if V_phase < 0.5 * V_limit:
            status = "✓ OK"
        elif V_phase < V_limit:
            status = "⚠ Near limit"
        else:
            status = "✗ OVER"
            
        if N == 200 or (k_fill == 0.10 and N in [100, 150, 200, 250]):  # Show key cases
            print(f"{k_fill:>8.3f} {N:>6d} {L_coil*1000:>12.1f} {V_phase/1000:>12.2f} {status:>15s}")
    if k_fill != 0.50:
        print()

# =============================================================================
# SECTION 6: RECOMMENDATIONS
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 6: RECOMMENDATIONS FOR ION PROPULSION ENGINEERING TEXTBOOK")
print("=" * 70)

recommendations = f"""
Based on analysis of Prof. Dr. Fedor Gömöry's experimental work and 
HTS coil design literature, updated for current baseline (January 2026):

CURRENT BASELINE (from lim_config.py):
  τp = {tau_p} m, W = {w_coil} m, gap = {gap} m
  N = {N_turns} turns
  I_peak = {I_peak} A (3mm × 1 layer HTS tape)
  V_limit = 100 kV

1. FILL FACTOR (k_fill)
   ─────────────────────
   Legacy value:     k_fill = 0.167 (too low - assumed only REBCO layer)
   Realistic range:  k_fill = 0.05-0.20 for space-rated insulation
   
   Justification:
   • HTS tape: ~100 μm total thickness
   • Space-rated insulation needs higher voltage standoff
   • At orbital ring scale (100 kV class), conservative k_fill is appropriate

2. VOLTAGE LIMIT
   ─────────────
   Current baseline: V_limit = 100 kV
   
   Justification:
   • Large-scale HTS systems can achieve high voltages with proper insulation
   • Japanese HTS cables achieve 275 kV with specialized PPLP insulation
   • 100 kV is achievable with proper vacuum insulation design
   • The large gap (200 mm) provides additional safety margin

3. OPTIMAL TURN COUNT
   ──────────────────
   Current baseline: N = 200 turns
   
   This is validated by the deployment simulation which shows:
   • Voltage stays well within limits throughout deployment
   • Thrust is adequate for ~14 month deployment
   • Higher N would increase thrust but risk voltage limits at high v_rel

4. KEY INSIGHT FROM GÖMÖRY'S WORK
   ───────────────────────────────
   His AC loss studies show that hysteresis loss dominates in HTS coils
   under AC conditions. For the LIM operating at ~15-40 Hz (varying with v_rel):
   
   • AC losses scale with frequency and (I/I_c)³ approximately
   • Striation of wide tapes (dividing into filaments) reduces AC loss
   • The Norris thin-strip model provides a theoretical basis for loss calculation
   
   The simulation includes both Norris-based and simplified hysteresis models.
"""

print(recommendations)

# Final summary table
print("\n" + "=" * 70)
print("SUMMARY: CURRENT LIM COIL PARAMETERS")
print("=" * 70)

k_fill_current = 0.10  # Typical for space-rated insulation
N_current = 200

L_coil_current, V_phase_current = calculate_voltage(N_current, k_fill_current, tau_p, w_coil, gap, f_supply, I_peak)

print(f"""
Current Baseline Parameters (lim_config.py January 2026):
─────────────────────────────────────────────────────────────────────
  Pole pitch (τp):        {tau_p} m
  Coil width (W):         {w_coil} m  
  Magnetic gap:           {gap} m
  Plate thickness:        0.20 m (γ-TiAl)
  Number of turns (N):    {N_current}
  HTS tape:               3mm × 1 layer → I_c = 200 A
  Peak current:           {I_peak} A (81.25% of I_c)
  Fill factor (k_fill):   ~{k_fill_current} (estimated)
  
At mid-deployment (v_rel = 3000 m/s, f = {f_supply:.1f} Hz):
  Coil inductance:        {L_coil_current*1000:.1f} mH
  Phase voltage:          {V_phase_current/1000:.2f} kV
  Voltage margin:         {(1 - V_phase_current/V_limit)*100:.0f}% below 100 kV limit
─────────────────────────────────────────────────────────────────────
""")

print("=" * 70)
