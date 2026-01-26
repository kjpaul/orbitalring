# LIM Deployment Simulator for Orbital Ring Systems

A physics-based simulator for modeling the deployment of an orbital ring using Linear Induction Motors (LIMs). This code accompanies **"Ion Propulsion Engineering"** (Technical Book II of the *Astronomy's Shocking Twist* series) by Paul G de Jong.

The simulator models the months-long process of decelerating a magnetically levitated cable from orbital velocity (7,755 m/s) to ground-stationary velocity (483 m/s at 250 km altitude) using distributed LIM sites around the ring.

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

### Run with Graphs (Interactive Display)

```bash
python lim_simulation.py --show thrust power current
```

### Save Graphs to Files

```bash
python lim_simulation.py --save all clean
```

### Select Thrust Model

```bash
python lim_simulation.py --model=1    # Narrow plate eddy current (default)
python lim_simulation.py --model=2    # Goodness factor (Laithwaite)
python lim_simulation.py --model=3    # Slip × pressure (theoretical max)
```

### Select Thermal Mode

```bash
python lim_simulation.py --thermal=cable   # Cable absorbs eddy heat (default)
python lim_simulation.py --thermal=cryo    # Cryo handles all heat
```

### Quick Test Run

```bash
python lim_simulation.py --quick           # Fast run with larger time steps
```

### Custom Mass Configuration

```bash
python lim_simulation.py --m_load=50000    # Set cable mass based on new load mass of 50,000 kg/m
```

### Full Help

```bash
python lim_simulation.py --help
```

---

## Graph Output

The simulator can save publication-quality graphs to files or display them interactively.

### Saving Graphs

```bash
python lim_simulation.py --save all clean              # Save all graphs to ./graphs/
python lim_simulation.py --outdir=./figures all clean  # Custom output directory
python lim_simulation.py --dpi=600 all clean           # Higher resolution for print
```

Graphs are saved with numbered filenames in logical order:
- `01-current.png`
- `02-volts.png`
- `03-v_slip.png`
- ... and so on

### Graph Configuration

Edit `lim_config.py` to set defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SAVE_GRAPHS` | True | Save to files (True) or display (False) |
| `GRAPH_OUTPUT_DIR` | "./graphs" | Output directory for saved graphs |
| `GRAPH_DPI` | 300 | Resolution (300 DPI is standard for print) |
| `GRAPH_WIDTH_INCHES` | 10 | Figure width in inches |
| `GRAPH_HEIGHT_INCHES` | 5 | Figure height in inches |
| `GRAPH_FORMAT` | "png" | Output format: "png", "pdf", "svg", "jpg" |

### Command-Line Overrides

| Option | Description |
|--------|-------------|
| `--save` | Force save to files |
| `--show` | Force interactive display |
| `--outdir=PATH` | Set output directory |
| `--dpi=N` | Set resolution |

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

### Module Organization

**lim_config.py** (Configuration):
- Section 1: User-configurable parameters
- Section 2: Physical constants
- Section 3: Material properties (5 materials)
- Section 4: Derived parameters
- Section 5: Physics parameter dictionary
- Section 6: Display utilities

**lim_physics.py** (Physics):
- Orbital dynamics (cable/casing velocity, deployment progress)
- Slip and frequency calculations
- Magnetic field calculations (B at plate, B in coil)
- Material properties (temperature-dependent resistivity)
- Eddy current losses (from thrust power balance)
- Thrust models (3 options)
- Thermal calculations (plate temperature, cable equilibrium)
- Cryogenic system sizing (COP, heat loads)
- Radiator calculations (width for heat rejection)
- Kinetic energy tracking

**lim_simulation.py** (Simulation):
- Data collection arrays (23 tracked quantities)
- Main simulation loop with 10-iteration controller
- Predictive power control with slip ratio adjustment
- Limit enforcement (voltage, power, temperature)
- Monthly progress display
- Plotting functions with peak/final annotations
- High-resolution graph export
- Command-line interface

---

## Configuration

All parameters are configured by editing `lim_config.py`. The most frequently-modified parameters appear at the top with detailed comments explaining valid ranges and trade-offs.

### HTS Tape Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HTS_TAPE_WIDTH_MM` | 3 | Tape width (3, 4, 6, or 12 mm) |
| `HTS_TAPE_LAYERS` | 1 | Parallel tape layers (1–3) |
| `IC_PER_MM_PER_LAYER` | 66.7 | Critical current density (A/mm/layer) |
| `DE_RATING_FACTOR` | 0.9 | Multilayer self-inductance compensation |
| `NORRIS_HYSTERESIS` | False | Use Norris formula for hysteresis loss |

With 3 mm tape and 1 layer, I_c = 200 A, giving I_target ≈ 163 A.

**Trade-off:** More layers → higher current capacity → more thrust, but more hysteresis losses and diminishing returns due to de-rating.

### LIM Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_TURNS` | 200 | Turns per phase coil |
| `TAU_P` | 100 m | Pole pitch (traveling wave wavelength / 2) |
| `W_COIL` | 2.0 m | Coil width |
| `GAP` | 0.20 m | Air gap to reaction plate |
| `T_PLATE` | 0.2 m | Reaction plate thickness |
| `PITCH_COUNT` | 3 | Pole pitches per LIM |

**Key constraint:** The narrow-plate model (Model 1) assumes `W_COIL << TAU_P`.

### Plate Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `W_PLATE` | 1.2 m | Reaction plate width |
| `N_PLATES_PER_SIDE` | 1 | LIM reaction plates per side of cable |
| `N_LEV_PLATES` | 1 | Levitation plates per side |
| `D_LEV` | 0.10 m | Levitation plate thickness |
| `W_LEV` | 1.4 m | Levitation plate width |

### Reaction Plate Material

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PLATE_MATERIAL` | `"gamma_titanium"` | Material selection |

**Options:** `"aluminum"`, `"cuni7030"`, `"titanium"`, `"alpha_titanium"`, `"gamma_titanium"`

The code includes full material properties for each option, including temperature-dependent resistivity.

| Material | ρ₂₉₃ₖ (nΩ·m) | α (1/K) | Density (kg/m³) |
|----------|-------------|---------|-----------------|
| Aluminum | 26.5 | 0.00366 | 2700 |
| CuNi 70/30 | 380 | 0.0004 | 8900 |
| Pure Ti | 420 | 0.0035 | 4500 |
| α₂-Ti₃Al | 500 | 0.0015 | 4200 |
| γ-TiAl | 750 | 0.0012 | 3900 |

### Site Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LIM_SPACING` | 500 m | Distance between LIM sites |

**Trade-off:** Closer spacing → faster deployment, higher plate temperature, limits τ_p length, less solar power per site.

### Operating Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `VOLTS_MAX` | 100 kV | Maximum induced coil voltage |
| `MAX_SITE_POWER` | 8 MW | Power limit per LIM site |
| `V_SLIP_MIN` | 5 m/s | Minimum slip velocity |
| `V_SLIP_MAX` | 200 m/s | Maximum slip velocity |

### Slip Control Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SLIP_RATIO_NORMAL` | 0.02 | Target slip ratio at full current (2%) |
| `SLIP_RATIO_REDUCED` | 0.005 | Reduced slip when power-limited (0.5%) |

### Thermal Mode Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EDDY_HEAT_TO_CABLE` | True | Thermal management strategy |
| `CABLE_EMISSIVITY` | 0.85 | Cable surface emissivity |
| `CABLE_SURFACE_AREA_PER_M` | 0.5 m²/m | Radiating surface area per meter |
| `COIL_MLI_EFFECTIVENESS` | 0.001 | MLI insulation effectiveness |
| `COIL_SURFACE_AREA_PER_SITE` | 10 m² | Approximate coil cryostat surface area |

### Mass Model

| Parameter | Default | Description |
|-----------|---------|-------------|
| `M_CABLE_STRUCTURAL` | 96,700 kg/m | Structural cable mass per meter |
| `M_LOAD_M` | 12,000 kg/m | Casing + payload per meter |
| `SIGMA_TARGET` | 12.633 GPa | Target post-deployment cable tension |

The `--m_load=` command-line option recalculates cable mass using the physics function `calc_cable_mass()` which determines the CNT cable cross-section needed to support the specified load at 250 km orbital altitude.

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

### Cable Heat Sink Mode (`--thermal=cable`) — Default

Eddy current heat remains in the cable, which warms until radiative equilibrium is reached (typically 300–450 K). Only hysteresis losses and coil environmental heat leak go through the cryogenic system.

**Advantages:**
- Dramatically smaller radiators (roughly 1/3 the size)
- Lower cryocooler power consumption
- Simpler heat path for eddy current losses

**Disadvantages:**
- Thermal expansion of cable (~0.3% strain at 350 K vs 70 K)
- Must manage radiative coupling between warm cable and cold coils

### Cryogenic Mode (`--thermal=cryo`)

All heat—eddy current losses, hysteresis losses, and environmental absorption—is removed through the cryogenic system at 70 K. This requires large radiators due to the thermodynamic penalty of lifting heat from 70 K to 300 K.

**Advantages:**
- Cable stays cold, minimizing thermal expansion
- Well-understood thermal path

**Disadvantages:**
- Requires massive radiators
- High cryocooler power consumption

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

Slip efficiency: $\eta = \frac{2sG}{1 + s^2G^2}$

**Best for:** Comparison with literature; assumes W > τ_p (wide plate).

### Model 3: Slip × Pressure (Theoretical Maximum)

Simple upper bound assuming perfect phase alignment:

$$F = s \times \frac{B^2 A}{2\mu_0}$$

**Best for:** Understanding theoretical limits.

---

## Available Plots

Specify plots as command-line arguments. When saving, plots are numbered in this order:

| # | Name | Description |
|---|------|-------------|
| 01 | `current` | Coil current (A) |
| 02 | `volts` | Induced voltage (V) |
| 03 | `v_slip` | Slip velocity (m/s) |
| 04 | `slip` | Slip ratio (%) |
| 05 | `f_slip` | Slip frequency (Hz) |
| 06 | `v_rel` | Relative velocity cable–casing (m/s) |
| 07 | `b_peak` | Magnetic field at plate (T) |
| 08 | `thrust` | Thrust force per site (N) |
| 09 | `p_thrust` | Thrust power per site (W) |
| 10 | `p_eddy` | Eddy current losses per site (W) |
| 11 | `hyst` | HTS hysteresis losses per site (W) |
| 12 | `p_heat` | Total heat losses per site (W) |
| 13 | `lim_power` | LIM power (W) |
| 14 | `power` | Total site power (W) |
| 15 | `cryo` | Cryogenic power requirement (W) |
| 16 | `q_cryo` | Cryo cold-side heat load (W) |
| 17 | `plate_temp` | Reaction plate temperature (K) |
| 18 | `cable_temp` | Cable equilibrium temperature (K) |
| 19 | `radiator_width` | Required radiator width (m) |
| 20 | `skin` | Effective skin depth (mm) |
| 21 | `skin_calc` | Calculated skin depth (mm) |
| 22 | `ke_site` | Cumulative kinetic energy per site (J) |
| 23 | `ke_all` | Cumulative kinetic energy all sites (J) |

Use `all` to generate all plots.

### Plot Display Options

| Option | Description |
|--------|-------------|
| `fulldata` | Show parameter list in title bar |
| `timeonly` | Show deployment days and years in title bar |
| `clean` | Only show title and axis labels |

**Examples:**
```bash
python lim_simulation.py --save all clean
python lim_simulation.py --show thrust power plate_temp cable_temp radiator_width
python lim_simulation.py --outdir=./book_figures --dpi=600 all timeonly
```

---

## Output

### Monthly Progress Display

The simulator prints monthly progress during deployment:

```
Month | Progress | Voltage | Current | Thrust | Site Power | Radiator W
---------------------------------------------------------------------------
    0 |    0.0% |    1234 V |  48.7 A |   134 N |  2.45 MW |   3.21 m
    1 |    2.3% |    2456 V | 163.0 A |   538 N |  5.12 MW |   4.67 m
   ...
```

### Final Summary

Final summary includes:
- Deployment time (days, years)
- Final cable and casing velocities
- Energy per site (TJ)
- Total energy delivered (EJ)
- Energy verification via ½mv² calculation
- Thermal mode and key thermal metrics
- Maximum cable temperature
- Maximum radiator width
- Maximum cryogenic cold-side heat load

### Max/Min Tracking

The simulation tracks peak values for 19 quantities including:
- Voltage, current, slip velocity
- Thrust, thrust power
- Eddy current and hysteresis losses
- Plate and cable temperatures
- Site power, cryogenic power
- Radiator width requirements

### Output File

Results are appended to `./output/_orbital_ring_modelN.txt` where N is the thrust model number.

---

## Physics Summary

### Orbital Parameters at 250 km

| Parameter | Value |
|-----------|-------|
| Orbital velocity | 7,754.866 m/s |
| Ground-stationary velocity | 483.331 m/s |
| Ring circumference | 41,645,813 m |
| Net acceleration | 9.038 m/s² |
| Orbital radius | 6,628,137 m |

### Magnetic Field at Reaction Plate

For a rectangular current sheet at gap distance g, with 3-phase traveling wave:

$$B = \frac{3\mu_0 N I_{peak}}{\pi W} \cdot \arctan\left(\frac{W}{2g}\right)$$

The factor of 1.5 (giving 3μ₀/π) accounts for the 3-phase traveling wave amplitude.

### Temperature-Dependent Resistivity

Resistivity varies with temperature:

$$\rho(T) = \rho_{293K} \times [1 + \alpha(T - 293)]$$

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

Default: η_cryo = 0.05 (5% of Carnot), T_cold = 70 K, T_hot = 400 K.

### Kinetic Energy Verification

The simulation tracks cumulative thrust power and verifies against analytical kinetic energy:

$$KE = \frac{1}{2} m_{per\_m} \cdot (v_{final}^2 - v_{initial}^2) \cdot L_{ring}$$

This provides a sanity check that energy is conserved to high precision.

---

## Typical Results

With default parameters (3mm × 1-layer HTS tape, 200 turns, 8 MW site power, γ-TiAl plate):

| Metric | Cable Mode |
|--------|------------|
| Deployment time | ~14 months |
| Peak thrust per site | ~2,800 N |
| Peak site power | 8 MW |
| Max cable temperature | ~370 K |
| Max radiator width | ~5 m |
| Total energy | ~15 EJ |

Results vary significantly with thrust model selection and configuration.

---

## Controller Logic

The simulation uses a proportional controller with predictive power management:

1. **Power margin assessment**: Calculate headroom relative to MAX_SITE_POWER
2. **Slip ratio targeting**: Interpolate between SLIP_RATIO_NORMAL (2%) and SLIP_RATIO_REDUCED (0.5%) based on power margin
3. **Current headroom blending**: Further reduce slip target when current is well below target
4. **Smooth actuator response**: Exponential smoothing with configurable time constant
5. **Limit enforcement**: Scale back current and slip if voltage, power, or temperature limits exceeded
6. **Startup ramp**: Gradual current increase during first hour

The controller iterates up to 10 times per time step to achieve self-consistency between thrust, losses, temperature, and power.

---

## Legacy Code

The original single-file version of the simulator is available as `lim_deployment_sim.py`. To use the legacy code:

```bash
python lim_deployment_sim.py --model=1 thrust power
```

The legacy code has all parameters and physics in a single file. It does not include the thermal mode switch, radiator width calculations, graph saving, or the `--m_load` option. Use it if you prefer a single-file solution or need to compare with earlier results.

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

> de Jong, P.G. (2025). *Orbital Ring Engineering: Mechanics and Material Science for Space Launch Mass Transit Systems, Volume I*. Astronomy's Shocking Twist Series.
