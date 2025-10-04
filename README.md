# Orbital Ring Deployment Simulator

This repository contains:
- **`orbital_ring.py`** — _current simulator_ (recommended): LIM-based deployment with coil geometry, voltage/current/power caps, controller ramps, and improved plotting/units.
- **`legacy_power_lim.py`** — _original reference code_ (unchanged except for the rename).

A minimal CI “smoke” run verifies that `orbital_ring.py` executes and emits JSON.

[![CI (smoke)](https://github.com/kjpaul/orbitalring/actions/workflows/ci.yml/badge.svg)](../../actions)

---

## Quick Start

```bash
pip install pyyaml
python orbital_ring.py --config example_lim_config.yaml   --plot velocities,stress,max_load_per_m,accel_excess
```

- Add `--save-plots` to write PNGs (no plots by default).
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

- `velocities` — casing/cable speeds (m/s) with end labels.
- `stress` — cable stress = \(T/A\) in **GPa**.
- `max_load_per_m` — \(\max(0,w')/g\) in **t/m** (_maximum supported load per meter_).
- `accel_excess` — mass-weighted \((a_\text{out} - g)\) in **m/s²** (positive = net outward).
- `net_weight_sign` — indicator \(-1,0,+1\) (line is **red** if any inward intervals occur).
- Diagnostics: `power` (MW), `voltage` (V), `current_slip` (A & m/s), `relative_speed` (m/s).

---

## Physics summary (v2)

Net outward “weight” per meter (N/m):
\[
w' = m'_\text{cable}\!\left(\frac{v_\text{cable}^2}{r} - g\right) +
     m'_\text{casing}\!\left(\frac{v_\text{casing}^2}{r} - g\right)
\]

Hoop tension proxy (only when \(w'>0\)):
\[
T \approx w' \, r, \qquad \sigma = T/A_\text{cable} \quad (\text{break at } 25\ \text{GPa})
\]

LIM thrust and power flow depend on coil geometry, slip/supply frequency, eddy coupling, and caps (voltage/current/inverter power). See comments in `orbital_ring.py` for the exact formulae.

---

## Repository layout

```
orbitalring/
├── orbital_ring.py            # ← current simulator you should run
├── legacy_power_lim.py        # ← original code for reference
├── example_lim_config.yaml
├── .github/workflows/ci.yml   # smoke test (PyYAML + short run)
└── README.md                  # you are here
```
