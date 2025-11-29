# Orbital Ring Deployment Simulator

This repository contains a physics-based simulator for deploying an orbital ring using linear induction motors (LIMs). The model includes coil geometry, voltage/current/power limits, controller ramps, thermal behavior, net outward weight, tension, and full cumulative energy accounting.

It accompanies the *Astronomy's Shocking Twist* technical volumes on large-scale space infrastructure.

---

# **How to Run the Simulator**

The simulator requires a YAML configuration file. **You must pass it using `--config`**. Positional arguments (e.g., `python orbital_ring.py example.yaml`) are **not supported**.

## **1. Install Dependencies**

```bash
pip install pyyaml matplotlib
```

## **2. Basic Run**

```bash
python orbital_ring.py --config example_lim_config.yaml
```

This runs the full simulation and prints a JSON result summary to stdout.

## **3. Plotting**

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

## **4. Full Example Run**

```bash
python orbital_ring.py \
    --config example_lim_config.yaml \
    --plot velocities,stress,max_load_per_m,total_ke \
    --save-plots \
    --outdir results
```

## **5. Helpful Developer Flags**

Verbose output including governing equations:
```bash
python orbital_ring.py --config example_lim_config.yaml --verbose --explain
```

Skip the first N seconds of simulation when plotting:
```bash
python orbital_ring.py --config example_lim_config.yaml --plot-skip-s 10000
```

---

# **Configuration File Overview**

All simulation parameters are stored in a YAML file such as `example_lim_config.yaml`.

### **Key Parameter Groups**

**Mass Model**
- `m_cable_per_m` — cable linear mass [kg/m]
- `m_casing_per_m` — casing linear mass [kg/m]
- `cable_area_m2` — cable load-bearing cross-section [m²]
- `sigma_break` — breaking stress [Pa]

**LIM Geometry & Limits**
- `n_turns`, `tau_p`, `w_coil`, `gap`, `pitch_count`, `spacing`
- `volts_max_user`, `max_site_power`
- Coil tape, fill factor, thermal emissivities

**Controller Ramps**
- `di_max_per_s` — max dI/dt
- `dvslip_max_per_s` — max change in slip-speed per second

**Run Controls**
- `dt` — timestep
- `max_time` — total simulated time (s)
- `plot` — list of enabled plots
- `plot_skip_s` — skip early transient data
- `save_plots`, `outdir`

---

# **Available Plots**

Enable any of these with `--plot name1,name2,...` or define them in the YAML file.

| Plot Name | Description |
|-----------|-------------|
| `velocities` | Casing & cable speeds (m/s) with end labels |
| `stress` | Cable stress T/A in GPa |
| `max_load_per_m` | max(0, w′)/g in t/m (maximum supported load/m) |
| `accel_excess` | Mass-weighted (a_out − g) in m/s² |
| `net_weight_sign` | Indicator −1/0/+1 (line turns red if inward force occurs) |
| `power` | Per-site input electrical power (MW) |
| `voltage` | Coil inductive voltage (V) |
| `current_slip` | Peak coil current (A) and slip velocity (m/s) |
| `relative_speed` | v_cable − v_casing (m/s) |
| `site_energy` | Cumulative site electrical energy (TJ), incl. losses |
| `total_energy` | Cumulative electrical energy all sites (EJ) |
| `site_ke` | Cumulative kinetic energy per site (TJ), thrust-only |
| `total_ke` | Cumulative kinetic energy all sites (EJ), thrust-only |

---

# **Output Metrics**

The simulator prints a JSON object containing key results:

| Metric | Meaning |
|--------|---------|
| `site_energy_TJ_2LIMs` | Total electrical energy per site (TJ), including losses |
| `total_energy_EJ_all` | Electrical energy all sites (EJ), including losses |
| `site_ke_TJ_2LIMs` | Mechanical kinetic energy per site (TJ), thrust only |
| `total_ke_EJ_all` | Mechanical kinetic energy all sites (EJ), thrust only |

**Energy separation:**
- Electrical input power → thrust + eddy + hysteresis + cryo + inefficiencies
- `ke` metrics integrate only **F × v_rel**, the mechanical work done

---

# **Physics Summary**

### **Net Outward Weight per Meter**

$$w' = m'_{cable}\left(\frac{v_{cable}^2}{r} - g\right) + m'_{casing}\left(\frac{v_{casing}^2}{r} - g\right)$$

### **Hoop Tension (when w′ > 0)**

$$T \approx w' r, \qquad \sigma = T/A_{cable}$$

Breaking stress assumed at **25 GPa**.

### **Mechanical Work from LIM Thrust**

$$P_{thrust} = F v_{rel}$$

The LIM transfers momentum between casing and cable:
- Cable gains ~29.8 EJ
- Casing loses ~15.0 EJ
- LIMs supply the ~14.8 EJ difference (not the sum)

This matches `total_ke_EJ_all` from the simulator.

---

# **Repository Layout**

```
orbitalring/
├── orbital_ring.py            # main simulator
├── legacy_power_lim.py        # reference version
├── example_lim_config.yaml    # sample configuration
├── .github/workflows/ci.yml   # smoke test
└── README.md
```

---

# **Related**

This simulation toolkit supports the engineering material presented in the *Astronomy's Shocking Twist* series on orbital ring megastructures.
