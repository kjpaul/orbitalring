# Space Infrastructure Lab (open-source)

**Numerical tools for designing and understanding orbital rings, space elevators, and mass drivers — built for accuracy, clarity, and learning.**

---

## Why this exists

While writing *How to Build an Orbital Ring, Vol. 1*, I began with a set of assumptions that turned out to be wrong. Some issues yielded to clean, analytic solutions, but the time-dependent processes (e.g., orbital-ring deployment) stubbornly demanded **iterative, numerical** approaches. Running simple Python scripts over and over led to a series of surprising, testable results (see the book’s Chapter 6). The code was crude and magical in equal measure.

This repository is the next step: an open-source project to **model, optimize, and explain** space-infrastructure systems, starting with **orbital rings**, **space elevators**, and **mass drivers**, and expanding to **trajectory design** (e.g., 4 g launch profiles from a specific ring) and related problems.

The core idea: a **Python interface** that anyone can run locally. It shows the equations, explains what they mean, and lets contributors improve models and algorithms. If the community likes the improvement, it will be added to the repository. The tool is both a **design environment** and a **learning environment** for curious minds who want to see how these systems actually work.

---

## Scope (initial + planned)

* **Orbital ring deployment**: 1-D/2-D time-dependent models, linear induction motor parameters, casing/cable interactions, energy budgets.
* **Space elevator tether**: tapered-cable sizing, stress margins, counterweight strategies, J2 and drag perturbation options.
* **Mass drivers & coilguns**: staging, eddy-current coupling, thermal limits, power electronics envelopes.
* **Trajectory design (planned)**: ring-to-LEO transfers, 4 g launch profiles, interplanetary windows, plane change economics.
* **Education mode**: inline derivations, parameter sensitivity plots, and “explain this equation” toggles.

---

## Design principles

1. **Transparency over black boxes**: Equations are shown next to code paths.
2. **Determinism**: Simulations are reproducible; seeds and integrator tolerances are explicit.
3. **Engineering-grade defaults**: Units, material properties, and environmental constants are sane, cited, and overrideable.
4. **Performance where it matters**: Pure Python first; vectorize/numba only when fidelity or speed demands it.

---

## Package layout (proposed)

```
space_infra_lab/
  __init__.py
  constants.py          # GM, ?, atmosphere models, material DB hooks
  materials/            # Cable & conductor property tables; user-extensible
  dynamics/
    gravity.py          # Point mass + J2
    rotation.py         # Earth rotation, frames
    drag.py             # Rarefied atmo models for low alt segments
  systems/
    orbital_ring.py     # Ring kinematics & deployment integrators
    space_elevator.py   # Taper sizing, stress profiles, counterweight models
    mass_driver.py      # Electromechanical models, staging, losses
  optimize/
    optimizers.py       # Grid search, CMA-ES, gradient-ish hooks
    objectives.py       # Mass, energy, throughput, cost functions
  io/
    config.py           # YAML/JSON schema & validation
    results.py          # HDF5/Parquet writers; plotting helpers
  viz/
    plots.py            # Matplotlib renderers; explain-mode charts
cli/
  sil.py                # `sil` command-line entry point
examples/
  ring_deploy.yaml
  elevator_taper.yaml
  mass_driver.yaml
```

---

## Getting started

### Requirements

* Python **3.11+**
* `numpy`, `scipy`, `matplotlib`, `pydantic`, `pyyaml`, `pandas` (installed automatically)

### Install (development mode)

```bash
git clone https://github.com/<your-org>/space-infra-lab.git
cd space-infra-lab
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e ".[dev]"    # installs package + dev tools (black, ruff, pytest)
```

### Quickstart: simulate a simple ring deployment

```python
from space_infra_lab.systems.orbital_ring import RingDeployment
from space_infra_lab.constants import EARTH

ring = RingDeployment(
    planet=EARTH,
    initial_speed=7755.0,          # m/s at ~250 km
    target_speed=18000.0,          # m/s cable target
    lim_spacing=1000.0,            # m
    lim_efficiency=0.85,
    casing_mass_per_m=1200.0,      # kg/m
    cable_mass_per_m=2700.0        # kg/m (example)
)

result = ring.run(dt=0.2, t_max=30_000)  # seconds
result.plot()  # shows velocity, tension, power vs time
```

### CLI (after install)

```bash
sil sim ring --config examples/ring_deploy.yaml --out runs/ring_001
sil sim elevator --config examples/elevator_taper.yaml --out runs/elev_001
```

---

## How the math shows up

Every model surfaces the governing equations in docstrings and `explain()` calls. For example, the elevator taper is implemented as:

$$
A(r) = A_0 \exp\!\left[\frac{\rho}{\eta\,\sigma}
\left(\tfrac{1}{2}\omega^2(r^2-r_0^2) - GM\Bigl(\frac{1}{r}-\frac{1}{r_0}\Bigr)\right)\right]
$$

with a companion `explain_taper()` that walks through assumptions, units, and typical parameter ranges.

---

## Contributing

Contributions are very welcome — from bug fixes to full subsystems.

1. **Fork** and create a feature branch.
2. **Code style**: `ruff` + `black`; **typing**: `mypy` (strict where practical).
3. **Tests**: add `pytest` coverage for new behavior; include a minimal example in `examples/`.
4. **Docs**: expand docstrings and `explain()` text for any new equations.
5. **PR**: open against `main` with a brief design note (what changed, why, references).

We follow the **Contributor Covenant v2.1** for community conduct.

---

## Roadmap

* **v0.1**: Orbital ring 1-D deployment integrator + power/energy accounting
* **v0.2**: Space elevator taper & counterweight explorer (GEO ± ?r)
* **v0.3**: Mass-driver stage model with eddy-current/heating limits
* **v0.4**: Trajectory module (ring-assisted 4 g launch; patched-conic seeds)
* **v0.5**: Optimization API (multi-objective: mass × cost × throughput)

---

## Educational mode

Toggle `explain=True` or run:

```bash
sil explain elevator --param sigma=1e11 rho=1300 eta=0.5
```

and the tool prints equations, assumptions, and sensitivity tables alongside plots.

---

## License

**Apache-2.0** (permissive, patent-grant).
If you prefer MIT or a dual-license, update `LICENSE` and `pyproject.toml` accordingly before your first public release.

---

## Disclaimer

This software is for **research and education**. It is **not** a substitute for professional engineering review, certification, or safety analysis. Real-world systems demand margins, fault trees, and environmental testing that exceed the scope of this codebase.

---

## Citation

If you use this project in academic or public work, please cite:

> Dutchy, *How to Build an Orbital Ring, Vol. 1* (forthcoming).
> Space Infrastructure Lab (GitHub repository), version X.Y.Z.

A BibTeX entry will be added after the first tagged release.

---

## Acknowledgments

Thanks to the early Python notebooks that proved what the equations wouldn’t — and to the community that will make them rigorous, reusable, and teachable.

---

### Contact

Open an issue or discussion on GitHub. For collaboration inquiries, please file a “Proposal” discussion with a short abstract and references.



