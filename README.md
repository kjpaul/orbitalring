# Orbital Ring LIM Deployment Simulator

Physics-based simulator for deploying an orbital ring using linear induction motors (LIMs) with high-temperature superconducting (HTS) coils. The model includes electromagnetic thrust, temperature-dependent aluminum resistivity, HTS hysteresis losses, cryogenic system sizing, and thermal radiation.

This code accompanies *Orbital Ring Engineering* (Volume I) and *Ion Propulsion Engineering* (Volume II) from the *Astronomy's Shocking Twist* technical series.

---

## Quick Start

### Install Dependencies

```bash
pip install tabulate matplotlib
```

### Run the Simulator

```bash
python legacy_power_lim.py
```

The simulator runs with the parameters defined at the top of the file and produces:
- Console output with deployment time and key metrics
- PNG graphs (if `MAKE_GRAPHS = True`)
- Data file (if `WRITE_FILE = True`)

---

## Configuration Overview

All parameters are defined as constants at the top of `legacy_power_lim.py`. The key configuration sections are:

### 1. HTS Tape Configuration

```python
HTS_TAPE_WIDTH_MM = 12      # Tape width: 12, 6, 4, or 3 mm
HTS_TAPE_LAYERS = 2         # Number of tape layers: 1 or 2
IC_PER_MM_PER_LAYER = 33.33 # Critical current per mm width per layer (A/mm)
```

These automatically compute:
- `W_TAPE` — tape width in meters
- `I_C` — critical current (A)
- `I_PEAK` — peak operating current, 87.5% of I_C (A)
- `I_TARGET` — target current, 81.25% of I_C (A)

**Standard HTS Tape Options:**

| Tape Width | Layers | I_C | I_PEAK | I_TARGET |
|------------|--------|-----|--------|----------|
| 12 mm | 2 | 800 A | 700 A | 650 A |
| 12 mm | 1 | 400 A | 350 A | 325 A |
| 6 mm | 2 | 400 A | 350 A | 325 A |
| 4 mm | 2 | 267 A | 233 A | 217 A |
| 3 mm | 2 | 200 A | 175 A | 162 A |
| 3 mm | 1 | 100 A | 87 A | 81 A |

### 2. LIM Site Configuration

```python
LIMS_PER_SIDE = 1           # LIMs on each side of cable (1, 2, 3, ...)
LIMS_PER_SITE = 2 * LIMS_PER_SIDE  # Total LIMs per site (computed)
```

LIMs are placed symmetrically on both sides of the cable. Setting `LIMS_PER_SIDE = 2` gives 4 total LIMs per site.

**Multi-LIM Trade-offs:**

More LIMs per site with narrower tape can be more efficient:
- Narrower tape → lower hysteresis per unit thrust
- More LIMs → higher total thrust for same cryo load
- Trade-off: more HTS tape length required

### 3. LIM Geometry

```python
N_TURNS = 25                # Turns per phase coil
TAU_P = 100.0               # Pole pitch (m)
W_COIL = 0.5                # Coil width (m)
GAP = 0.20                  # Air gap to reaction plate (m)
PITCH_COUNT = 3             # Pole pitches per LIM
LIM_SPACING = 500.0         # Distance between LIM sites (m)
```

Derived quantities:
- `L_ACTIVE = TAU_P × PITCH_COUNT` — active LIM length (m)
- `A_LIM = L_ACTIVE × W_COIL` — active area under LIM (m²)
- `LIM_SITES = L_RING / LIM_SPACING` — total number of sites around the ring

### 4. Power and Electrical Limits

```python
MAX_SITE_POWER = 16.0e6     # Power limit per site (W)
VOLTS_MAX = 100e3           # Peak coil voltage limit (V)
I_MIN = 10.0                # Minimum current before reducing slip (A)
```

### 5. Reaction Plate (Aluminum)

```python
T_PLATE = 0.150             # Aluminum plate thickness (m)
T_MAX_Reaction_Plate = 500  # Maximum plate temperature (K)
EM_ALU = 0.85               # Plate emissivity (black anodized)
```

### 6. Cryogenic System

```python
CRYO_EFF = 0.05             # Cryo efficiency (fraction of Carnot)
T2_LN2_BOIL = 77.4          # LN2 boiling point (K)
T_N2_HOT = 300              # Radiator temperature (K)
EM_HEAT_SINK = 0.9          # Radiator emissivity
```

The cryogenic coefficient of performance (COP) is:
```
COP = (T_cold / (T_hot - T_cold)) × CRYO_EFF
```

With default values: COP ≈ 0.017, meaning ~58 W electrical power per 1 W of heat removed at 77 K.

### 7. Thrust Model

```python
THRUST_EFFICIENCY = 1.0     # Multiplier on idealized thrust model
```

The thrust model uses:
```
F = THRUST_EFFICIENCY × slip × F_max
```

where `F_max = B² × A / (2μ₀)` is the magnetic pressure limit.

- `THRUST_EFFICIENCY = 1.0` is conservative (idealized resistive regime)
- Higher values (5-20) may be realistic for high magnetic Reynolds number operation
- Requires FEM simulation to validate for specific geometry

### 8. Controller Parameters

```python
SLIP_TARGET = 0.01          # Target slip ratio (v_slip / v_wave)
SLIP_TARGET_START = 0.10    # Initial slip ratio
CURRENT_UPRATE = 1.01       # Current ramp multiplier per iteration
POWER_HEADROOM = 0.98       # Fraction of MAX_SITE_POWER to target
```

### 9. Mass Model

```python
M_CABLE_M = 99_198          # Cable mass per meter (kg/m)
M_LOAD_M = 12_000           # Casing + load mass per meter (kg/m)
```

---

## Example Configurations

### Configuration A: Baseline (2 LIMs, 12mm tape)

```python
HTS_TAPE_WIDTH_MM = 12
HTS_TAPE_LAYERS = 2
LIMS_PER_SIDE = 1
N_TURNS = 25
```

Results: I_C = 800 A, 2 LIMs per site

### Configuration B: High-Efficiency (4 LIMs, 3mm tape)

```python
HTS_TAPE_WIDTH_MM = 3
HTS_TAPE_LAYERS = 2
LIMS_PER_SIDE = 2
N_TURNS = 100  # Increase turns to maintain B field
```

Results: I_C = 200 A, 4 LIMs per site, ~4× better thrust/hysteresis efficiency

### Configuration C: Low-Power (4 LIMs, 3mm tape, same turns)

```python
HTS_TAPE_WIDTH_MM = 3
HTS_TAPE_LAYERS = 2
LIMS_PER_SIDE = 2
N_TURNS = 25   # Same turns → lower B field
```

Results: Much lower thrust but dramatically reduced cryo requirements

---

## Physics Summary

### Magnetic Field at Reaction Plate

Uses finite-width current sheet model:

```
B_plate = (2μ₀ × N × I_peak) / (π × W) × arctan(W / 2g)
```

This correctly saturates as gap → 0, unlike the naive μ₀NI/g formula.

### Temperature-Dependent Resistivity

Aluminum resistivity varies with temperature:

```
ρ(T) = ρ_ref × [1 + α × (T - T_ref)]
```

where α ≈ 0.00366 K⁻¹. At 77 K, resistivity is ~30% of room temperature value.

### Skin Depth

```
δ = √(ρ / (π × μ₀ × f_slip))
```

Effective depth is min(δ, t_plate). Lower temperature → lower resistivity → shallower skin depth.

### Thrust Model

```
F = THRUST_EFFICIENCY × s × F_max
F_max = B² × A_active / (2μ₀)
```

where s = slip ratio = f_slip / f_supply.

### Eddy Current Losses

```
P_eddy = F × v_slip
```

This is the power dissipated in the reaction plate, which must be radiated away.

### HTS Hysteresis Losses

```
q = B_coil × I × W_tape × sin(α)
P_hyst = q × L_HTS × f_supply
```

Key insight: **Hysteresis per unit thrust scales with tape width.** Narrower tape is inherently more efficient.

### Cryogenic Power

```
P_cryo = P_heat / COP
COP = (T_cold / (T_hot - T_cold)) × CRYO_EFF
```

The cryo radiator must reject:
```
Q_hot = P_heat + P_cryo = P_heat × (1 + 1/COP)
```

### Radiator Sizing

```
A_radiator = Q_hot / (ε × σ × (T_hot⁴ - T_space⁴))
```

At 300 K radiating to 3 K space with ε = 0.9:
- Radiates ~410 W/m²
- For 10 MW rejection: ~24,000 m² required

---

## Output Metrics

The simulator prints:

| Metric | Description |
|--------|-------------|
| Deployment time | Time to reach orbital velocity (years) |
| Cable velocity | Final cable velocity (m/s) |
| Casing velocity | Final casing velocity (m/s) |
| Site KE | Kinetic energy per site (TJ) |
| Total KE | Total kinetic energy all sites (EJ) |
| Cryo radiator width | Minimum radiator width per site (m) |

### Generated Plots

When `MAKE_GRAPHS = True`, produces:
- Current and voltage vs time
- Thrust and thrust power vs time
- Eddy and hysteresis losses vs time
- Plate temperature vs time
- Cryo power vs time
- Relative velocity vs time

---

## Key Scaling Relationships

### Thrust
```
F ∝ B² ∝ (N × I)²
```

To maintain thrust with lower current, increase turns proportionally.

### Hysteresis
```
P_hyst ∝ B × I × W_tape × L_HTS × f
      ∝ (N × I) × I × W_tape × N × f
      ∝ N × I² × W_tape × f
```

For constant B (N × I = const):
```
P_hyst ∝ W_tape
```

**Narrower tape gives proportionally lower hysteresis for same thrust.**

### Efficiency Improvement

Switching from 12mm to 3mm tape (4× narrower):
- Same N × I product maintains B field
- Hysteresis drops by 4×
- Thrust/hysteresis efficiency improves 4×
- Requires 4× more turns (4× more tape length)
- Net tape cost similar if $/m scales with width

---

## File Structure

```
legacy_power_lim.py     # Main simulator
README.md               # This file
```

Output files (generated):
```
graphs/                 # Plot images (PNG)
data/                   # Simulation data files
```

---

## Dependencies

- Python 3.7+
- `tabulate` — table formatting
- `matplotlib` — plotting

---

## Related Resources

This simulation supports the engineering analysis in:
- *Orbital Ring Engineering* (Volume I) — Cable mechanics, material science, LIM fundamentals
- *Ion Propulsion Engineering* (Volume II) — Propulsion systems, mission analysis

Part of the *Astronomy's Shocking Twist* series combining technical textbooks with hard science fiction.

---

## License

MIT License
