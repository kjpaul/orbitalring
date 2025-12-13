# LIM Physics Formulas - Verified Summary for Orbital Ring Engineering

This document summarizes all physics formulas used in the orbital ring LIM simulation,
with verification that each formula has correct limiting behavior (no division by zero
or infinity issues).

## 1. Magnetic Field at Reaction Plate

**Corrected Formula:**
```
B_plate = (2μ₀NI)/(πW) × arctan(W/2g)
```

**Why the correction?** The naive formula B = μ₀NI/g assumes an infinite current sheet
and blows up as g→0. The corrected formula accounts for finite coil width W.

**Limiting behavior:**
- As g → 0: arctan(W/2g) → π/2, so B → μ₀NI/W (finite)
- As g → ∞: arctan(W/2g) ≈ W/2g, so B → μ₀NI/(πg) (correct far-field)

**Numerical comparison (N=100, I=650A, W=2m):**
| Gap (cm) | Naive B (T) | Corrected B (mT) | Ratio |
|----------|-------------|------------------|-------|
| 20       | 0.41        | 35.7             | 11×   |
| 10       | 0.82        | 38.2             | 21×   |
| 5        | 1.63        | 39.5             | 41×   |
| 2        | 4.08        | 40.3             | 101×  |
| → 0      | → ∞         | 40.8             | ∞     |

**Simulation code: CORRECT** (uses arctan formula in b_plate_peak function)

---

## 2. Coil Inductance

**Corrected Formula:**
```
L_coil = (2μ₀N²A_coil×k_fill)/(πW) × arctan(W/2g)
```

where A_coil = τp × W (pole pitch × coil width)

**Derivation:** L = NΦ/I where Φ = B × A_coil. Substituting the corrected B formula
and simplifying gives the corrected L formula.

**Limiting behavior:**
- As g → 0: L → μ₀N²A_coil×k_fill/W (finite)
- As g → ∞: L → μ₀N²A_coil×k_fill/(πg) (correct far-field)

**Numerical comparison (N=100, τp=50m, W=2m, k_fill=0.074):**
| Gap (cm) | Naive L (H) | Corrected L (mH) | Ratio |
|----------|-------------|------------------|-------|
| 20       | 0.46        | 40.7             | 11×   |
| 10       | 0.93        | 43.6             | 21×   |
| 5        | 1.86        | 45.1             | 41×   |
| 2        | 4.66        | 46.0             | 101×  |

**Simulation code:** Not used directly (voltage uses flux linkage, which is equivalent)

---

## 3. Induced Voltage

**Formula (two equivalent methods):**
```
Method 1: V_rms = ωLI/√2
Method 2: V_rms = ωNBA_coil×k_fill/√2
```

Both methods give identical results when using the corrected B and L formulas.

**Frequency relationship:**
```
f_supply = v_supply/(2τp)    where v_supply = v_rel + v_slip
```

**Simulation code: CORRECT** (uses Method 2 with corrected B)

---

## 4. Fill Factor

**Formula:**
```
k_fill = d_HTS / (d_HTS + d_Kapton)
```

**Physical meaning:** Fraction of winding cross-section that is conductor vs insulation.

**Typical values:**
| Insulation | d_HTS | d_Kapton | k_fill |
|------------|-------|----------|--------|
| Thin (38μm film + 60μm adhesive) | 80 μm | 98 μm | 0.449 |
| Thick (1mm) | 80 μm | 1000 μm | 0.074 |

**⚠️ CONFIG NOTE:** The YAML uses k_fill=0.449 with d_kapton=1mm. These are 
inconsistent (k_fill=0.449 is for thin Kapton). This is conservative (overestimates 
voltage) but should be made consistent.

---

## 5. Voltage Limits

**Breakdown voltage:**
```
V_breakdown = (N-1) × E_breakdown × d_Kapton
```

where E_breakdown ≈ 200 kV/mm for Kapton.

**Design voltage (with safety factor):**
```
V_design = (N-1) × E_design × d_Kapton
```

where E_design = 10-50 kV/mm typically.

**Example (N=100, d_Kapton=1mm):**
- V_breakdown = 99 × 200 kV/mm × 1mm = 19.8 MV
- V_design = 99 × 10 kV/mm × 1mm = 990 kV

**Simulation code: CORRECT** (volts_max_allowed function)

---

## 6. Skin Depth

**Formula:**
```
δ = √(ρ/(πμ₀f_slip))
```

**Effective depth:** δ_eff = min(δ, d_plate)

**Limiting behavior:**
- f_slip → 0: δ → ∞ (full penetration at DC)
- f_slip → ∞: δ → 0 (surface currents only)

**Typical values (aluminum at 150K, ρ = 1.26×10⁻⁸ Ω·m):**
| f_slip (Hz) | δ (mm) |
|-------------|--------|
| 1           | 56     |
| 2           | 40     |
| 5           | 25     |
| 10          | 18     |

**Simulation code: CORRECT** (skin_depth function handles f≤0)

---

## 7. Aluminum Resistivity (Temperature Dependent)

**Formula:**
```
ρ(T) = ρ_20°C × [1 + α × (T - 293)]
```

where α ≈ 0.00366 K⁻¹ (≈ 1/273)

**Typical values:**
| T (K) | ρ (×10⁻⁸ Ω·m) | Relative |
|-------|---------------|----------|
| 293   | 2.65          | 1.00     |
| 200   | 1.75          | 0.66     |
| 150   | 1.26          | 0.48     |
| 100   | 0.77          | 0.29     |

**Simulation code: CORRECT** (rho_alu_at_temp function with minimum floor)

---

## 8. Goodness Factor

**Formula:**
```
G = (ω_slip × μ₀ × σ × δ_eff × τp) / π
```

where σ = 1/ρ is conductivity.

**Physical meaning:** Ratio of induced EMF to resistive voltage in reaction plate.
- G >> 1: Good motor behavior
- G << 1: Resistive losses dominate

**Optimal slip ratio:** s_opt = 1/G

**Typical values (v_slip=200 m/s, τp=50m, T=150K):**
G ≈ 800 (very large due to large τp)

**Simulation code: CORRECT** (goodness_factor function)

---

## 9. Slip Ratio and Efficiency

**Slip ratio:**
```
s = v_slip / (v_rel + v_slip) = f_slip / f_supply
```

**Slip efficiency:**
```
η_slip = 2sG / (1 + s²G²)
```

**Properties:**
- Maximum η = 1 when s = 1/G (and G >> 1)
- η → 0 as s → 0 (synchronous, no thrust)
- η → 0 as s → ∞ (excessive slip, all losses)

**Limiting behavior:** Denominator (1 + s²G²) ≥ 1, so no division issues.

**Simulation code: CORRECT** (slip_ratio and slip_efficiency in compute_lim_outputs)

---

## 10. Thrust Force

**Formula:**
```
F = (B²/2μ₀) × A_active × η_slip
```

where A_active = W × L_active = W × (τp × pitch_count)

**Physical interpretation:** Magnetic pressure times area times efficiency factor.

**Maximum magnetic pressure:**
```
P_mag = B²/(2μ₀)
```

For B = 39.5 mT: P_mag = 620 Pa ≈ 0.006 atm

**Simulation code: CORRECT** (uses corrected B, proper η_slip)

---

## 11. Eddy Current Losses

**Formula:**
```
P_eddy = (π²/6ρ) × B² × δ_eff² × f_slip² × V_eddy
```

where V_eddy = W × L_active × δ_eff

**Limiting behavior:**
- ρ → 0: δ → 0, so losses remain bounded (self-regulating)
- f → 0: P_eddy → 0 (no losses at DC)

**Simulation code: CORRECT** (resistivity has minimum floor)

---

## 12. Hysteresis Losses (HTS AC Losses)

**Formula:**
```
P_hyst = q × L_HTS × f_supply × 3 × pitch_count
```

where:
- q = B_coil × Ic × w_tape × sin(α) (loss per unit length per cycle)
- L_HTS = 2(W + τp) × N (total HTS length per phase)
- Factor of 3 for three phases

**Simulation code: CORRECT**

---

## 13. Equilibrium Plate Temperature

**Formula:**
```
T = (P/(εσA))^(1/4)
```

where σ = 5.67×10⁻⁸ W/(m²·K⁴) is Stefan-Boltzmann constant.

**Simulation code: CORRECT** (minimum temperature floor enforced)

---

## 14. Total Site Power

**Formula:**
```
P_site = (P_thrust + P_eddy + P_hyst + P_cryo) / (η_inv × η_lim)
```

**Physical interpretation:** Total electrical power required to:
- Deliver mechanical thrust power P_thrust
- Supply eddy current losses (radiated from plate)
- Supply hysteresis losses (removed by cryocoolers)
- Run cryogenic system
- Account for inverter and wiring losses

**Power limiting:** The controller reduces current to keep P_site ≤ P_max (8 MW).

**Simulation code: CORRECT**

---

## Summary of Issues Found

1. **k_fill / d_kapton inconsistency in YAML:**
   - YAML has k_fill=0.449 (thin Kapton) but d_kapton=1mm (thick)
   - Should use k_fill=0.074 for 1mm Kapton, or d_kapton=0.098mm for thin
   - Current config is conservative (overestimates voltage)

2. **All physics formulas verified correct:**
   - B-field uses arctan correction ✓
   - Voltage uses corrected B via flux linkage ✓
   - No division-by-zero issues ✓
   - All minimum floors in place ✓

---

## Recommended YAML Configuration (Consistent)

For 1mm Kapton insulation:
```yaml
lim:
  d_kapton_mm: 1.0
  k_fill: 0.074        # = 80/(80+1000)
  e_allow_kv_per_mm: 10.0
```

For thin Kapton insulation (38μm + 60μm adhesive):
```yaml
lim:
  d_kapton_mm: 0.098   # 98 μm total
  k_fill: 0.449        # = 80/(80+98)
  e_allow_kv_per_mm: 10.0
```
