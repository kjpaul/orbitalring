# LIM Deployment Simulator for Orbital Ring Systems

A physics-based simulator for modeling the deployment of an orbital ring using Linear Induction Motors (LIMs). This code accompanies **"Ion Propulsion Engineering"** (Technical Book I of the *Astronomy's Shocking Twist* series) by Paul G de Jong.

The simulator models the months-long process of decelerating a magnetically levitated cable from orbital velocity (7.75 km/s) to ground-stationary velocity (483 m/s at 250 km altitude) using distributed LIM sites around the ring.

---

## Quick Start

### Installation

```bash
pip install matplotlib tabulate
```

### Run with Defaults

```bash
python lim_deployment_sim.py
```

### Run with Graphs

```bash
python lim_deployment_sim.py thrust power current
```

### Select Thrust Model

```bash
python lim_deployment_sim.py --model=1    # Narrow plate eddy current (default)
python lim_deployment_sim.py --model=2    # Goodness factor (Laithwaite)
python lim_deployment_sim.py --model=3    # Slip × pressure (theoretical max)
```

### Full Help

```bash
python lim_deployment_sim.py --help
```

---

## Configuration

All parameters are configured by editing constants in Section 1 of `lim_deployment_sim.py`. The code is organized so that frequently-modified parameters appear at the top with detailed comments explaining valid ranges and trade-offs.

### HTS Tape Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HTS_TAPE_WIDTH_MM` | 3 | Tape width (3, 4, 6, or 12 mm) |
| `HTS_TAPE_LAYERS` | 3 | Parallel tape layers (1–3) |
| `IC_PER_MM_PER_LAYER` | 66.7 | Critical current density (A/mm/layer) |

**Trade-off:** Wider tape → higher current capacity → more thrust, more hysteresis loses.
**Trade-off:** More layers → higher current capacity → more thrust, model may be less accurate.

### LIM Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_TURNS` | 80 | Turns per phase coil |
| `TAU_P` | 100 m | Pole pitch (traveling wave wavelength / 2) |
| `W_COIL` | 0.5 m | Coil width |
| `GAP` | 0.05 m | Air gap to reaction plate |
| `T_PLATE` | 0.05 m | Aluminium plate thickness |
| `PITCH_COUNT` | 3 | Pole pitches per LIM |

**Key constraint:** `W_COIL << TAU_P` for narrow-plate model validity.

### Site Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LIM_SPACING` | 500 m | Distance between LIM sites |
| `LIMS_PER_SIDE` | 1 | LIMs per side of cable (×2 for site total) |

**Trade-off:** Closer spacing → faster deployment, higher plate temperature, limits length of TAU_P, less solar power per LIM site.

### Operating Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `VOLTS_MAX` | 100 kV | Maximum induced coil voltage |
| `MAX_SITE_POWER` | 16 MW | Power limit per LIM site |
| `V_SLIP_MIN` | 5 m/s | Minimum slip velocity |
| `V_SLIP_MAX` | 200 m/s | Maximum slip velocity |

### Mass Model

| Parameter | Default | Description |
|-----------|---------|-------------|
| `M_CABLE_M` | 99,198 kg/m | Cable mass per meter |
| `M_LOAD_M` | 12,000 kg/m | Casing + payload per meter |

---

## Thrust Models

The simulator includes three physics models for LIM thrust, selectable via `--model=N`:

### Model 1: Narrow Plate Eddy Current (Default)

Accounts for the orbital ring geometry where coil width W is much smaller than pole pitch τ_p. Eddy currents form elongated loops with high return-path resistance:

```
←←←←←← current under negative pole ←←←←←←
↓                                        ↑
↓ return path (length τ_p = 100m)        ↑
↓                                        ↑
→→→→→→ current under positive pole →→→→→→
```

**Best for:** Orbital ring design where W << τ_p. This is probably the most realistic model.

### Model 2: Goodness Factor (Laithwaite)

Classic LIM theory using Laithwaite's goodness factor G:

$$G = \frac{\omega_{slip} \cdot \mu_0 \cdot \sigma \cdot \delta_{eff} \cdot \tau_p}{\pi}$$

Thrust efficiency: $\eta = \frac{2sG}{1 + s^2G^2}$

**Best for:** Comparison with literature; assumes W > τ_p. Needs FEA to test validity. Fastest deployment times.

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
python lim_deployment_sim.py --model=1 thrust power plate_temp cryo
```

---

## Output

The simulator prints monthly progress during deployment:

```
Month | Progress | Voltage | Current | Thrust | Site Power
-----------------------------------------------------------------
    0 |    0.0% |    1234 V |  48.7 A |   134 N |  2.45 MW
    1 |    2.3% |    2456 V | 195.2 A |   538 N |  8.12 MW
   ...
```

Final summary includes:
- Deployment time (days, years)
- Final cable and casing velocities
- Total energy delivered (EJ)
- Energy per LIM site (TJ)

---

## Physics Summary

### Magnetic Field at Reaction Plate

For a rectangular current sheet at gap distance g:

$$B = \frac{3\mu_0 N I_{peak}}{\pi W} \cdot \arctan\left(\frac{W}{2g}\right)$$

The factor of 1.5 (giving 3/2) accounts for the 3-phase traveling wave.

### Temperature-Dependent Resistivity

Aluminium resistivity varies with temperature:

$$\rho(T) = \rho_{293K} \times [1 + \alpha(T - 293)]$$

where ρ₂₉₃ₖ = 2.65×10⁻⁸ Ω·m and α = 3.663×10⁻³ K⁻¹.

At cryogenic temperatures (~70 K), resistivity drops to ~18% of room temperature.

### Skin Depth

$$\delta = \sqrt{\frac{\rho}{\pi \mu_0 f_{slip}}}$$

Effective depth = min(δ, plate_thickness).

### Slip Power Balance

Secondary (reaction plate) losses equal slip power:

$$P_{eddy} = F \times v_{slip}$$

This self-consistent approach prevents unphysical runaway in the thermal model.

### Cryogenic Power

Using Carnot efficiency with a practical efficiency factor:

$$COP = \frac{T_{cold}}{T_{hot} - T_{cold}} \times \eta_{cryo}$$

$$P_{cryo} = \frac{Q_{heat}}{COP}$$

Default: η_cryo = 0.18 (18% of Carnot), T_cold = 77 K, T_hot = 300 K.

---

## Code Structure

The code is organized into 20 clearly-labeled sections:

| Sections | Content |
|----------|---------|
| 1–3 | Configuration (user parameters, constants, derived values) |
| 4–6 | Orbital dynamics, slip/frequency calculations |
| 7–8 | Magnetic fields, material properties |
| 9–11 | Eddy losses, goodness factor, thrust models |
| 12–13 | Electrical calculations, HTS hysteresis |
| 14–15 | Thermal model, cryogenic system |
| 16–17 | Utilities, parameter display |
| 18–20 | Main simulation loop, plotting, entry point |

All physics functions include docstrings with equations and references.

---

## Typical Results

With default parameters (3mm × 3-layer HTS tape, 80 turns, 16 MW site power):

| Metric | Value |
|--------|-------|
| Deployment time | ~18 months |
| Peak thrust per LIM | ~500 N |
| Peak site power | 16 MW (limit) |
| Total energy (all sites) | ~15 EJ |

Results vary significantly with thrust model selection and HTS tape configuration.

---

## Related Resources

- **Orbital Ring Engineering** (Volume I) — Cable mechanics, anchor systems, material science
- **Ion Propulsion Engineering** (Volume II) — LIM physics, deployment analysis
- Legacy simulation: `legacy_power_lim.py` with YAML configuration (see `README_LEGACY.md`)

---

## License

MIT License — see LICENSE file for details.

---

## Citation

If you use this code in academic work, please cite:

> de Jong, P.G. (2025). *Orbital Ring Engineering: Mechanics and Material Science for Space Lauch Mass Transit Systems, Volume I*. Astronomy's Shocking Twist Series.
