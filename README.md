# Orbital Ring Deployment Simulator (CLI)

**`orbital_ring_sim.py`** is a clean, commented, command‑line tool for simulating the **tangential momentum exchange** between an orbital ring’s **outer casing** and **inner cable** during deployment. It prints what it’s doing, explains the equations (with `--explain`), and can optionally plot results.

> **Fixed kinematics by altitude (default 250 km circular orbit):**
> - Initial orbital speed: **7,754.866 m/s**
> - Final casing (ground‑sync) speed: **483.331 m/s**
> - Local gravity at altitude: **9.073 m/s²**
> - Orbital radius: **~6,628,137 m** (Earth radius + 250 km)
> - Circumference: **~41,645,813 m**

## Quick Start

```bash
# 1) Grab the script and optional example config
python orbital_ring_sim.py --help

# 2) Run with simple constant force (prints per‑second status and explanations)
python orbital_ring_sim.py --model constant_force --force 1e6 --dt 0.5 --verbose --explain

# 3) Or run with constant power and plot velocities
python orbital_ring_sim.py --model constant_power --power 5e9 --dt 0.25 --plot velocities
```

If you prefer editing a file instead of many CLI flags:

```bash
python orbital_ring_sim.py --config example_config.yaml --plot velocities
```

If `PyYAML` is not installed, either install it (`pip install pyyaml`) or pass parameters via CLI flags.

## What the Code Simulates

- A time‑dependent loop integrates
  \\[
    \\frac{dv_\\text{casing}}{dt} = \\frac{F_\\text{casing}}{M_\\text{casing}}, \\qquad
    \\frac{dv_\\text{cable}}{dt}  = \\frac{F_\\text{cable}}{M_\\text{cable}}
  \\]
  with **Newton’s third law** (internal LIM force): \\(F_\\text{cable} = -F_\\text{casing}\\).
- Two thrust models:
  - **`constant_force`**: apply a fixed |F| until the casing hits target.
  - **`constant_power`**: use \\(F = P / v\\) (capped by `--fmax` if set).
- The loop **stops** when the casing slows to **483.331 m/s** (within tolerance).

As a cross‑check, the script also prints the **final cable speed** implied by **angular‑momentum conservation** in a closed system:
\\[
v_{\\text{cable,final}}=\\frac{(M_\\text{total} v_\\text{orbit}) - (M_\\text{casing} v_\\text{casing,final})}{M_\\text{cable}}.
\\]

## Plots (optional)

- `--plot velocities` — time histories of casing and cable speeds
- `--plot tension` — **first‑order** total cable tension proxy (placeholder; see TODO)
- `--plot weight` — sign of net radial load (+outward, −inward), **schematic**

Save plots instead of showing them interactively:
```bash
python orbital_ring_sim.py --plot velocities --save-plots --outdir results
```

## Parameters You’ll Likely Edit

- `--m-cable-per-m` and `--m-casing-per-m` (kg/m): mass budgets
- `--model`, `--force` or `--power`, `--fmax` (N): actuation model
- `--dt`, `--max-time`, `--tol-v` (s): integration controls
- `--plot`, `--save-plots`, `--outdir`: outputs

Or keep these in a YAML file (see `example_config.yaml`).

## Notes & Roadmap

- **Tension & Net Weight**: The current implementation includes *placeholders* for total tension and a simple radial sign check. You can extend:
  - Use distributed ring mechanics to compute section tensions (Fourier‑mode approach).
  - Include gravity vs centrifugal terms per subsystem; incorporate coupler loads.
- **Losses & EM Detail**: Add LIM slip, electrical/cryogenic losses, and power budgets.
- **Controllers**: Add higher‑level controllers (e.g., accelerate within a power envelope to reach the target in a specified time).

## License

MIT — feel free to fork, extend, and PR improvements.
