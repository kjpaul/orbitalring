# LIM Deployment Simulator for Orbital Ring Systems

A physics-based simulator for modeling the deployment of an orbital ring using Linear Induction Motors (LIMs). This code accompanies **"Ion Propulsion Engineering"** (Technical Book II of the *Astronomy's Shocking Twist* series) by Paul G de Jong.

The simulator models the months-long process of decelerating a magnetically levitated cable from orbital velocity (7.75 km/s) to ground-stationary velocity (483 m/s at 250 km altitude) using distributed LIM sites around the ring.

---

## Quick Start

### Installation

```bash
pip install matplotlib tabulate
```

### Run with Defaults

```bash
python lim_simulation.py
```

### Run with Graphs

```bash
python lim_simulation.py thrust power current
```

### Select Thrust Model

```bash
python lim_simulation.py --model=1    # Narrow plate eddy current (default)
python lim_simulation.py --model=2    # Goodness factor (Laithwaite)
python lim_simulation.py --model=3    # Slip × pressure (theoretical max)
```

### Select Thermal Mode

```bash
python lim_simulation.py --thermal=cryo    # Cryo handles all heat (default)
python lim_simulation.py --thermal=cable   # Cable absorbs eddy heat
```

### Quick Test Run

```bash
python lim_simulation.py --quick           # Fast run with larger time steps
```

### Full Help

```bash
python lim_simulation.py --help
```

---

## Code Structure

The simulator is organized into three modules:

| File | Description |
|------|-------------|
| `lim_simulation.py` | Main entry point, simulation loop, plotting |
| `lim_config.py` | All configurable parameters and material properties |
| `lim_physics.py` | Physics calculation functions |

This modular structure allows you to:
- Modify parameters without touching physics code
- Import physics functions for standalone calculations
- Test individual components independently

### Section Organization

**lim_config.py** (Configuration):
- Section 1: User-configurable parameters
- Section 2: Physical constants
- Section 3: Material properties
- Section 4: Derived parameters
- Section 5: Physics parameter dictionary
- Section 6: Display utilities

**lim_physics.py** (Physics):
- Orbital dynamics
- Slip and frequency calculations
- Magnetic field calculations
- Material properties (temperature-dependent)
- Eddy current losses
- Thrust models (3 options)
- Thermal calculations
- Cryogenic system sizing
- Radiator calculations

**lim_simulation.py** (Simulation):
- Data collection arrays
- Main simulation loop
- Controller logic
- Plotting functions
- Command-line interface

---

## Configuration

All parameters are configured by editing `lim_config.py`. The most frequently-modified parameters appear at the top with detailed comments explaining valid ranges and trade-offs.

### HTS Tape Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HTS_TAPE_WIDTH_MM` | 3 | Tape width (3, 4, 6, or 12 mm) |
| `HTS_TAPE_LAYERS` | 2 | Parallel tape layers (1–3) |
| `IC_PER_MM_PER_LAYER` | 66.7 | Critical current density (A/mm/layer) |

**Trade-off:** Wider tape → higher current capacity → more thrust, but more hysteresis losses.

### LIM Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_TURNS` | 100 | Turns per phase coil |
| `TAU_P` | 100 m | Pole pitch (traveling wave wavelength / 2) |
| `W_COIL` | 0.5 m | Coil width |
| `GAP` | 0.20 m | Air gap to reaction plate |
| `T_PLATE` | 0.2 m | Reaction plate thickness |
| `PITCH_COUNT` | 3 | Pole pitches per LIM |

**Key constraint:** `W_COIL << TAU_P` for narrow-plate model validity.

### Reaction Plate Material

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PLATE_MATERIAL` | "titanium" | Material selection |

**Options:** `"aluminum"`, `"cuni7030"`, `"titanium"`, `"alpha_titanium"`, `"gamma_titanium"`

The code includes full material properties for each option, including temperature-dependent resistivity.

### Site Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LIM_SPACING` | 500 m | Distance between LIM sites |
| `N_PLATES_PER_SIDE` | 1 | Reaction plates per side of cable |

**Trade-off:** Closer spacing → faster deployment, higher plate temperature, limits TAU_P length, less solar power per site.

### Operating Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `VOLTS_MAX` | 100 kV | Maximum induced coil voltage |
| `MAX_SITE_POWER` | 16 MW | Power limit per LIM site |
| `V_SLIP_MIN` | 5 m/s | Minimum slip velocity |
| `V_SLIP_MAX` | 200 m/s | Maximum slip velocity |

### Thermal Mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EDDY_HEAT_TO_CABLE` | False | Thermal management strategy |

**Options:**
- `False` (default): All heat goes through cryogenic system at 70 K
- `True`: Eddy current heat stays in cable, which warms to equilibrium temperature

See the **Thermal Modes** section below for details.

### Mass Model

| Parameter | Default | Description |
|-----------|---------|-------------|
| `M_CABLE_STRUCTURAL` | 96,700 kg/m | Structural cable mass per meter |
| `M_LOAD_M` | 12,000 kg/m | Casing + payload per meter |

---

## Quick Mode

For rapid testing, use `--quick` to run with larger time steps:

```bash
python lim_simulation.py --quick thrust power
python lim_simulation.py --quick --thermal=cable radiator_width
```

Quick mode uses:
- 10× larger time steps (10s → 100s → 500s instead of 1s → 10s → 50s)
- More frequent data collection for smooth plots

This reduces run time from several minutes to under a minute, with some loss of accuracy in the transient behavior. The final steady-state values are typically within a few percent of the full simulation.

---

## Thermal Modes

The simulator supports two thermal management strategies, selectable via `--thermal=` or by setting `EDDY_HEAT_TO_CABLE` in the config.

### Cryogenic Mode (`--thermal=cryo`)

All heat—eddy current losses, hysteresis losses, and environmental absorption—is removed through the cryogenic system at 70 K. This requires large radiators due to the thermodynamic penalty of lifting heat from 70 K to 300 K (approximately 22.9 W rejected per W of cooling).

**Advantages:**
- Cable stays cold, minimizing thermal expansion
- Well-understood thermal path

**Disadvantages:**
- Requires massive radiators (~50 m² per watt of 70 K cooling)
- High cryocooler power consumption

### Cable Heat Sink Mode (`--thermal=cable`)

Eddy current heat remains in the cable, which warms until radiative equilibrium is reached (typically 300–350 K). Only hysteresis losses and coil environmental heat leak go through the cryogenic system.

**Advantages:**
- Dramatically smaller radiators (roughly 1/3 the size)
- Lower cryocooler power consumption
- Simpler heat path for eddy current losses

**Disadvantages:**
- Thermal expansion of cable (~0.3% strain at 350 K vs 70 K)
- Must manage radiative coupling between warm cable and cold coils

### Comparing Modes

Run both modes and compare radiator requirements:

```bash
python lim_simulation.py --thermal=cryo radiator_width cable_temp
python lim_simulation.py --thermal=cable radiator_width cable_temp
```

---

## Thrust Models

The simulator includes three physics models for LIM thrust, selectable via `--model=N`:

### Model 1: Narrow Plate Eddy Current (Default)

Accounts for the orbital ring geometry where coil width W is much smaller than pole pitch τ_p. Eddy currents form elongated loops with high return-path resistance and significant loop inductance.

```
←←←←←← current under negative pole ←←←←←←
↓                                        ↑
↓ return path (length τ_p = 100m)        ↑
↓                                        ↑
→→→→→→ current under positive pole →→→→→→
```

The loop impedance Z = √(R² + X²) where:
- R = ρ × τ_p / (δ_eff × W) — dominated by the long return path
- X = ω × L_loop — significant due to the elongated loop geometry

**Best for:** Orbital ring design where W << τ_p. Most physically accurate for this application.

### Model 2: Goodness Factor (Laithwaite)

Classic LIM theory using Laithwaite's goodness factor G:

$$G = \frac{\omega_{slip} \cdot \mu_0 \cdot \sigma \cdot \delta_{eff} \cdot \tau_p}{\pi}$$

Thrust efficiency: $\eta = \frac{2sG}{1 + s^2G^2}$

**Best for:** Comparison with literature; assumes W > τ_p (wide plate).

### Model 3: Slip × Pressure (Theoretical Maximum)

Simple upper bound assuming perfect phase alignment:

$$F = s \times \frac{B^2 A}{2\mu_0}$$

**Best for:** Understanding theoretical limits.

---

## Available Plots

Specify plots as command-line arguments:

| Name | Description |
|------|-------------|
| `all` | Show all available plots |
| `current` | Coil current (A) |
| `volts` | Induced voltage (V) |
| `v_slip` | Slip velocity (m/s) |
| `thrust` | Thrust force (N) |
| `p_thrust` | Thrust power (W) |
| `p_eddy` | Eddy current losses (W) |
| `power` | Total site power (W) |
| `plate_temp` | Reaction plate temperature (K) |
| `cable_temp` | Cable equilibrium temperature (K) |
| `radiator_width` | Required radiator width (m) |
| `q_cryo` | Cryo cold-side heat load (W) |
| `skin` | Electromagnetic skin depth (mm) |
| `slip` | Slip ratio (%) |
| `f_slip` | Slip frequency (Hz) |
| `v_rel` | Relative velocity cable–casing (m/s) |
| `hyst` | HTS hysteresis losses (W) |
| `cryo` | Cryogenic power requirement (W) |
| `ke_site` | Cumulative kinetic energy per site (J) |
| `ke_all` | Cumulative kinetic energy all sites (J) |

**Example:**
```bash
python lim_simulation.py --model=1 thrust power plate_temp cable_temp radiator_width
```

---

## Output

The simulator prints monthly progress during deployment:

```
Month | Progress | Voltage | Current | Thrust | Site Power | Radiator W
---------------------------------------------------------------------------
    0 |    0.0% |    1234 V |  48.7 A |   134 N |  2.45 MW |   3.21 m
    1 |    2.3% |    2456 V | 195.2 A |   538 N |  8.12 MW |   5.67 m
   ...
```

Final summary includes:
- Deployment time (days, years)
- Final cable and casing velocities
- Total energy delivered (EJ)
- Energy per LIM site (TJ)
- Thermal mode and key thermal metrics
- Maximum cable temperature
- Maximum radiator width
- Maximum cryogenic cold-side heat load

---

## Physics Summary

### Magnetic Field at Reaction Plate

For a rectangular current sheet at gap distance g, with 3-phase traveling wave:

$$B = \frac{3\mu_0 N I_{peak}}{\pi W} \cdot \arctan\left(\frac{W}{2g}\right)$$

The factor of 1.5 (giving 3μ₀/π) accounts for the 3-phase traveling wave amplitude.

### Temperature-Dependent Resistivity

Resistivity varies with temperature:

$$\rho(T) = \rho_{293K} \times [1 + \alpha(T - 293)]$$

| Material | ρ₂₉₃ₖ (nΩ·m) | α (1/K) |
|----------|-------------|---------|
| Aluminum | 26.5 | 0.00366 |
| Pure Ti | 420 | 0.0035 |
| γ-TiAl | 750 | 0.0012 |

### Skin Depth

$$\delta = \sqrt{\frac{\rho}{\pi \mu_0 f_{slip}}}$$

Effective depth = min(δ, plate_thickness).

### Slip Power Balance

Secondary (reaction plate) losses equal slip power:

$$P_{eddy} = F \times v_{slip}$$

This self-consistent approach prevents unphysical runaway in the thermal model.

### Cable Equilibrium Temperature

When `EDDY_HEAT_TO_CABLE = True`, the cable temperature is found from radiative equilibrium:

$$T_{eq} = \left( \frac{P_{eddy,total} / L_{ring}}{\varepsilon \sigma A_{surface/m}} + T_{space}^4 \right)^{0.25}$$

### Cryogenic Power and Radiator Sizing

COP with practical efficiency:

$$COP = \frac{T_{cold}}{T_{hot} - T_{cold}} \times \eta_{cryo}$$

Heat rejected at radiator (must include both heat lifted and compressor work):

$$Q_{reject} = Q_{cold} \times \left(1 + \frac{1}{COP}\right)$$

Radiator width for given length L:

$$W_{radiator} = \frac{Q_{reject}}{\varepsilon \sigma (T_{hot}^4 - T_{space}^4) \times L}$$

Default: η_cryo = 0.18 (18% of Carnot), T_cold = 70 K, T_hot = 300 K.

---

## Typical Results

With default parameters (3mm × 3-layer HTS tape, 80 turns, 16 MW site power, titanium plate):

| Metric | Cryo Mode | Cable Mode |
|--------|-----------|------------|
| Deployment time | ~18 months | ~17.5 months |
| Peak thrust per LIM | ~500 N | ~500 N |
| Peak site power | 16 MW | 16 MW |
| Max cable temperature | ~77 K | ~350 K |
| Max radiator width | ~15 m | ~5 m |
| Total energy | ~15 EJ | ~15 EJ |

Results vary significantly with thrust model selection and configuration.

---

## Legacy Code

The original single-file version of the simulator is available as `lim_deployment_sim.py`. To use the legacy code:

```bash
python lim_deployment_sim.py --model=1 thrust power
```

The legacy code has all parameters and physics in a single file. It does not include the thermal mode switch or radiator width calculations. Use it if you prefer a single-file solution or need to compare with earlier results.

**Note:** The legacy code uses `PLATE_MATERIAL = "titanium"` by default but references aluminum in some comments. The modular code has been updated for consistency.

---

## Related Resources

- **Orbital Ring Engineering** (Volume I) — Cable mechanics, anchor systems, material science
- **Ion Propulsion Engineering** (Volume II) — LIM physics, deployment analysis
- Orbital Ring Engineering website: https://www.orbitalring.space/coding-hub/

---

## License

MIT License — see LICENSE file for details.

---

## Citation

If you use this code in academic work, please cite:

> de Jong, P.G. (2025). *Ion Propulsion Engineering: Linear Induction Motors for Orbital Ring Deployment, Volume II*. Astronomy's Shocking Twist Series.

> de Jong, P.G. (2025). *Orbital Ring Engineering: Mechanics and Material Science for Space Launch Mass Transit Systems, Volume I*. Astronomy's Shocking Twist Series.