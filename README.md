# Orbital Ring Deployment Simulator

This repository contains:
- **`orbital_ring.py`** — _current simulator_ (recommended): LIM-based deployment with coil geometry, voltage/current/power caps, controller ramps, and improved plotting/units.
- **`legacy_power_lim.py`** — _original reference code_ (unchanged except for the rename).

A minimal CI "smoke" run verifies that `orbital_ring.py` executes and emits JSON.

[![CI (smoke)](https://github.com/kjpaul/orbitalring/actions/workflows/ci.yml/badge.svg)](../../actions)

---

## Quick Start
```bash
pip install pyyaml matplotlib
python orbital_ring.py --config example_lim_config.yaml \
  --plot velocities,stress,max_load_per_m,total_ke
```

- Add `--save-plots` to write PNGs (no GUI display).
- `--verbose --explain` prints governing equations and periodic status (1 s → 1 h → 1 d → 1 wk → monthly).

> **Altitude setpoint:** 250 km (equator).  
> `v_orbit = 7754.866 m/s`, `v_casing_final = 483.331 m/s`, `g = 9.073 m/s²`.

---

## Configuration

Edit `example_lim_config.yaml`. Key groups:

- **Mass model:** `m_cable_per_m`, `m_casing_per_m`, `cable_area_m2`, `sigma_break`.
- **LIM geometry & limits:** `n_turns`, `tau_p`, `w_coil`, `gap`, `pitch_count`, `spacing`, `volts_max_user`, `max_site_power`, etc.
- **Controller ramps:** `di_max_per_s`, `dvslip_max_per_s`.
- **Run controls:** `dt`, `max_time`, `plot`, `save_plots`, `plot_skip_s`.

---

## Plots (opt-in via `--plot`)

| Plot name | Description |
|-----------|-------------|
| `velocities` | Casing/cable speeds (m/s) with end labels |
| `stress` | Cable stress T/A in **GPa** |
| `max_load_per_m` | max(0, w′)/g in **t/m** (maximum supported load per meter) |
| `accel_excess` | Mass-weighted (a_out − g) in **m/s²** (positive = net outward) |
| `net_weight_sign` | Indicator −1, 0, +1 (line is **red** if any inward intervals occur) |
| `power` | Per-site input power in **MW** |
| `voltage` | Coil inductive voltage in **V** |
| `current_slip` | Peak current (A) and slip velocity (m/s) |
| `relative_speed` | v_cable − v_casing in **m/s** |
| `site_energy` | Cumulative total energy per LIM site in **TJ** (includes all losses) |
| `total_energy` | Cumulative total energy all sites in **EJ** (includes all losses) |
| `site_ke` | Cumulative kinetic energy per LIM site in **TJ** (thrust power only) |
| `total_ke` | Cumulative kinetic energy all sites in **EJ** (thrust power only) |

---

## Output Metrics

The simulator outputs JSON with key results:

| Metric | Description |
|--------|-------------|
| `site_energy_TJ_2LIMs` | Total electrical energy per site (TJ), including losses |
| `total_energy_EJ_all` | Total electrical energy all sites (EJ), including losses |
| `site_ke_TJ_2LIMs` | Kinetic energy per site (TJ), thrust power only |
| `total_ke_EJ_all` | Kinetic energy all sites (EJ), thrust power only |

The "ke" (kinetic energy) metrics track only the mechanical work done accelerating the cable and decelerating the casing (F × v_rel), excluding eddy current losses, hysteresis, cryogenic overhead, and inverter inefficiencies.

---

## Energy Physics

The simulation confirms a key result from Newton's third law applied to momentum and energy:

- **Cable gains:** ~29.8 EJ of kinetic energy (accelerating from v_orbit to v_cable)
- **Casing loses:** ~15.0 EJ of kinetic energy (decelerating from v_orbit to v_casing_final)
- **LIMs supply:** ~14.8 EJ (the difference)

Newton's third law guarantees equal and opposite *impulses*, not equal energy changes. The LIM acts at the interface between cable and casing, transferring momentum (and the energy carried by that momentum) from the casing to the cable. The LIMs only supply the 14.8 EJ shortfall—not the sum of both energy changes.

The `total_ke_EJ_all` output confirms this: integrating thrust power (F × v_rel) over the full deployment yields **14.8 EJ**.

---

## Physics Summary

Net outward "weight" per meter (N/m):

$$w' = m'_\text{cable}\left(\frac{v_\text{cable}^2}{r} - g\right) + m'_\text{casing}\left(\frac{v_\text{casing}^2}{r} - g\right)$$

Hoop tension proxy (only when w′ > 0):

$$T \approx w' \cdot r, \qquad \sigma = T/A_\text{cable} \quad (\text{break at } 25\ \text{GPa})$$

LIM thrust power (mechanical work into cable/casing):

$$P_\text{thrust} = F \times v_\text{rel}$$

Total site power includes thrust + eddy losses + hysteresis + cryo, divided by inverter and LIM efficiencies. See comments in `orbital_ring.py` for the exact formulae.

---

## Repository Layout
```
orbitalring/
├── orbital_ring.py            # ← current simulator you should run
├── legacy_power_lim.py        # ← original code for reference
├── example_lim_config.yaml
├── .github/workflows/ci.yml   # smoke test (PyYAML + short run)
└── README.md                  # you are here
```

---

## Related

This simulator accompanies the book *[Astronomy's Shocking Twist](https://orbitalring.space)* series exploring orbital ring megastructure engineering.