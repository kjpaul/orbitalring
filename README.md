# Orbital Ring Deployment Simulator

This repository contains a physics-based simulator for deploying an orbital ring using linear induction motors (LIMs). The model includes coil geometry, voltage/current/power limits, proportional controller, temperature-dependent aluminum resistivity, thermal radiation, and full cumulative energy accounting.

It accompanies the *Astronomy's Shocking Twist* technical volumes on large-scale space infrastructure.

---

## How to Run the Simulator

The simulator requires a YAML configuration file. **You must pass it using `--config`**. Positional arguments (e.g., `python orbital_ring.py example.yaml`) are **not supported**.

### 1. Install Dependencies

```bash
pip install pyyaml matplotlib
```

### 2. Basic Run

```bash
python orbital_ring.py --config example_lim_config.yaml
```

This runs the full simulation and prints a JSON result summary to stdout.

### 3. Plotting

To display specific plots:
```bash
python orbital_ring.py --config example_lim_config.yaml --plot velocities,stress
```

To save PNGs instead of showing them:
```bash
python orbital_ring.py --config example_lim_config.yaml \
    --plot velocities,stress --save-plots
```

You may include any set of comma-separated plot names.

### 4. Full Example Run

```bash
python orbital_ring.py \
    --config example_lim_config.yaml \
    --plot velocities,stress,power,plate_temp,losses \
    --save-plots \
    --outdir results
```

### 5. Helpful Developer Flags

Verbose output including governing equations:
```bash
python orbital_ring.py --config example_lim_config.yaml --verbose --explain
```

Skip the first N seconds of simulation when plotting:
```bash
python orbital_ring.py --config example_lim_config.yaml --plot-skip-s 10000
```

---

## Configuration File Overview

All simulation parameters are stored in a YAML file such as `example_lim_config.yaml`.

### Key Parameter Groups

**Mass Model**
| Parameter | Description | Units |
|-----------|-------------|-------|
| `m_cable_per_m` | Cable linear mass | kg/m |
| `m_casing_per_m` | Casing linear mass | kg/m |
| `cable_area_m2` | Cable load-bearing cross-section | m² |
| `sigma_break` | Breaking stress | Pa |

**LIM Geometry**
| Parameter | Description | Units |
|-----------|-------------|-------|
| `n_turns` | Turns per phase coil | — |
| `tau_p` | Pole pitch | m |
| `w_coil` | Coil width | m |
| `gap` | Air gap to reaction plate | m |
| `pitch_count` | Number of pole pitches per LIM | — |
| `spacing` | Distance between LIM sites | m |

**Electrical Limits**
| Parameter | Description | Units |
|-----------|-------------|-------|
| `i_peak_target` | Target peak current | A |
| `ic` | HTS critical current | A |
| `volts_max_user` | User-set voltage limit | V |
| `max_site_power` | Maximum power per site | W |
| `d_kapton_mm` | Kapton insulation thickness | mm |
| `e_allow_kv_per_mm` | Allowable electric field | kV/mm |

**Reaction Plate & Thermal**
| Parameter | Description | Units |
|-----------|-------------|-------|
| `t_plate` | Aluminum plate thickness | m |
| `t_plate_min` | Minimum plate temperature (floor) | K |
| `rho_alu_e_20C` | Aluminum resistivity at 20°C | Ω·m |
| `alpha_alu_e` | Temperature coefficient of resistivity | 1/K |
| `em_alu` | Plate emissivity | — |
| `heat_sink_l` | Heat sink extension length | m |

**Controller Ramps**
| Parameter | Description | Units |
|-----------|-------------|-------|
| `di_max_per_s` | Maximum current ramp rate | A/s |
| `dvslip_max_per_s` | Maximum slip velocity ramp rate | m/s² |

**Run Controls**
| Parameter | Description |
|-----------|-------------|
| `dt` | Timestep (s) |
| `max_time` | Maximum simulated time (s) |
| `plot` | List of enabled plots |
| `plot_skip_s` | Skip early transient data |
| `save_plots` | Save plots to files |
| `outdir` | Output directory |
| `graph_dir` | Subdirectory for graphs |

---

## Available Plots

Enable any of these with `--plot name1,name2,...` or define them in the YAML file.

| Plot Name | Description |
|-----------|-------------|
| `velocities` | Casing & cable speeds (m/s) with end labels |
| `stress` | Cable stress T/A in GPa |
| `max_load_per_m` | max(0, w′)/g in t/m (maximum supported load/m) |
| `accel_excess` | Mass-weighted (a_out − g) in m/s² |
| `net_weight_sign` | Indicator −1/0/+1 (line turns red if inward force occurs) |
| `power` | Per-site input electrical power (MW) with limit line |
| `voltage` | Phase voltage RMS (kV) with limit line |
| `current_slip` | Peak coil current (A) and slip velocity (m/s), dual axis |
| `relative_speed` | v_cable − v_casing (m/s) |
| `frequency` | Supply and slip frequencies (Hz) |
| `losses` | Eddy current (kW) and hysteresis (W) losses, dual axis |
| `plate_temp` | Reaction plate equilibrium temperature (K) with max annotation |
| `cryo` | Cryogenic system heat load (kW) |
| `site_energy` | Cumulative site electrical energy (TJ), incl. losses |
| `total_energy` | Cumulative electrical energy all sites (EJ) |
| `site_ke` | Cumulative kinetic energy per site (TJ), thrust-only |
| `total_ke` | Cumulative kinetic energy all sites (EJ), thrust-only |

---

## Output Metrics

The simulator prints a JSON object containing key results:

| Metric | Meaning |
|--------|---------|
| `months_elapsed` | Deployment time in 30-day months |
| `v_casing_end` | Final casing velocity (m/s) |
| `v_cable_end` | Final cable velocity (m/s) |
| `stress_max_GPa` | Peak cable stress (GPa) |
| `site_energy_TJ_2LIMs` | Total electrical energy per site (TJ), including losses |
| `total_energy_EJ_all` | Electrical energy all sites (EJ), including losses |
| `site_ke_TJ_2LIMs` | Mechanical kinetic energy per site (TJ), thrust only |
| `total_ke_EJ_all` | Mechanical kinetic energy all sites (EJ), thrust only |

**Energy separation:**
- Electrical input power → thrust + eddy + hysteresis + cryo + inefficiencies
- `ke` metrics integrate only **F × v_rel**, the mechanical work done

---

## Physics Summary

### Magnetic Field at Reaction Plate

The simulator uses a finite-width current sheet model that correctly saturates as gap → 0:

$$B_{plate} = \frac{2\mu_0 N I_{peak}}{\pi W} \cdot \arctan\left(\frac{W}{2g}\right)$$

This avoids the unphysical B → ∞ behavior of the naive μ₀NI/g formula.

### Temperature-Dependent Resistivity

Aluminum resistivity varies significantly with temperature:

$$\rho(T) = \rho_{20°C} \times \left(1 + \alpha (T - 293)\right)$$

where α ≈ 1/273 K⁻¹. At cryogenic temperatures, resistivity drops to ~30% of room temperature values, affecting skin depth and eddy current losses.

### Skin Depth

$$\delta = \sqrt{\frac{\rho_{alu}(T)}{\pi \mu_0 f_{slip}}}$$

The effective depth is min(δ, t_plate). Lower resistivity at cold temperatures means shallower skin depth.

### LIM Goodness Factor

$$G = \frac{\omega_{slip} \mu_0 \sigma_{alu} \delta_{eff} \tau_p}{\pi}$$

This dimensionless ratio determines the optimal slip point (s_opt = 1/G) and the slip-dependent efficiency.

### Thrust Model

$$F = \frac{B^2}{2\mu_0} \times A_{active} \times \eta_{slip}$$

where the slip efficiency factor is:

$$\eta_{slip} = \frac{2sG}{1 + s^2 G^2}$$

### Reaction Plate Temperature

The plate radiates eddy current heat to space via Stefan-Boltzmann:

$$T_{plate} = \left(\frac{P_{eddy}}{\varepsilon \sigma A_{radiating}}\right)^{0.25}$$

A minimum temperature floor (default 100 K) is enforced to reflect thermal contact with the cryogenic structure.

### Net Outward Weight per Meter

$$w' = m'_{cable}\left(\frac{v_{cable}^2}{r} - g\right) + m'_{casing}\left(\frac{v_{casing}^2}{r} - g\right)$$

### Hoop Tension (when w′ > 0)

$$T \approx w' \cdot r, \qquad \sigma = T/A_{cable}$$

Breaking stress assumed at **25 GPa**.

### Mechanical Work from LIM Thrust

$$P_{thrust} = F \cdot v_{rel}$$

The LIM transfers momentum between casing and cable:
- Cable gains ~29.8 EJ kinetic energy
- Casing loses ~15.0 EJ kinetic energy  
- LIMs supply the ~14.8 EJ difference (not the sum)

This matches `total_ke_EJ_all` from the simulator.

---

## Code Architecture

The refactored simulator uses a clean separation of concerns:

### State and Output Dataclasses

```python
@dataclass
class LIMState:
    """State variables that evolve each timestep."""
    i_peak: float       # Current peak current (A)
    v_slip: float       # Current slip velocity (m/s)
    T_plate: float      # Reaction plate temperature (K)

@dataclass
class LIMOutputs:
    """Computed outputs for current state (derived quantities)."""
    F_site: float       # Thrust per LIM (N)
    P_thrust: float     # Thrust power (W)
    P_eddy: float       # Eddy current losses (W)
    P_hyst: float       # Hysteresis losses (W)
    # ... and more
```

### Single Source of Truth

All LIM physics calculations happen in one function:

```python
def compute_lim_outputs(state: LIMState, v_rel: float, p: Params) -> LIMOutputs:
    """Compute all LIM outputs from current state."""
    # Temperature-dependent resistivity
    # Skin depth, goodness factor, slip efficiency
    # Thrust, eddy losses, hysteresis losses
    # Voltage estimate, total power
    # ... all in one place, no duplication
```

### Proportional Controller

The controller uses proportional feedback with a dead band to avoid oscillations:

- Over limit (>105%): Scale back proportionally
- Near limit (95-105%): Hold steady
- Under limit (80-95%): Ramp up gradually
- Well under (<80%): Ramp to target

---

## Repository Layout

```
orbitalring/
├── orbital_ring.py            # main simulator (refactored)
├── legacy_power_lim.py        # reference version
├── example_lim_config.yaml    # sample configuration
├── .github/workflows/ci.yml   # smoke test
└── README.md
```

---

## Related

This simulation toolkit supports the engineering material presented in the *Astronomy's Shocking Twist* series on orbital ring megastructures.

---

## License

MIT License - see LICENSE file for details.
