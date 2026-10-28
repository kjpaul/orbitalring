# The Orbital Ring, Volume IV: Space Mass Transit Complete Design Review

Simulation code by Paul G. de Jong, part of *Astronomy's Shocking Twist*, a ten-part series of five technical works and five novels.

This repository holds the calculations behind *The Orbital Ring*, the first of the five technical works, published in four volumes: Volume I, *The Economic Case for Space Mass Transit*; Volume II, *First Principles: Orbital Mechanics, Materials Science, and Electrical Engineering*; Volume III, *Mass Transit Launch Technologies: Levitation, Linear Motors, Mass Drivers, and Space Elevators*; Volume IV, *Space Mass Transit Complete Design Review*. The scripts here produce the numbers, tables and figures of Volume IV: an orbital ring at 700 km, a cable moving retrograde inside a stationary casing, held off the casing by a superconducting levitation bearing, deployed by linear induction motors, and used as a launch track by a linear synchronous motor mass driver. The first edition of the book, which the manuals in `docs/` were written for, is *Orbital Ring Engineering: Mechanics and Material Science for Space Launch Mass Transit Systems* (October 16, 2025).

Anyone can clone this repository, install the listed packages, run one command per subject, and get the numbers printed in the book.

## The design point

Set in `ring_altitude.py`, which every script reads (values from its log, `ring_altitude_700.txt`):

| Quantity | Value |
|---|---|
| Ring altitude | 700 km, equatorial |
| Cable direction | retrograde (westward) |
| Ring circumference | 44,473 km |
| Casing mass per meter | 12,000 kg/m |
| Cable hardware per meter (motor plates and bearing plates) | 10,686 kg/m |
| Cable mass per meter, structural | 115,005 kg/m |
| Cable mass per meter, with hardware | 125,691 kg/m |
| Cable speed | 8,270 m/s |
| Cable speed relative to the casing | 8,786 m/s |
| Cable tension | 845.6 GN |
| Cable cross-section | 67.65 m² (8.22 m square) |
| Ring mass | 6,124 million tonnes |
| Deployment energy | 18.8 EJ |
| Suspended solar array | 600 km, 100 km below the ring |

## Install

Python 3.10 or later. Then, from the repository root:

```bash
pip install -r requirements.txt
```

This installs numpy, scipy, matplotlib, pymsis (the NRLMSIS atmosphere model, for `altitude_trade.py`) and scikit-fem (the independent check of the field solver in `bearing_fem_check.py`). The logs behind the book were made with Python 3.14, numpy 2.4, scipy 1.18, matplotlib 3.11, pymsis 0.13 and scikit-fem 12.0, and a clean installation of the same list reproduced them.

Optional packages, only for the scripts named: `tabulate` for nicer tables in `lim_simulation.py` (the book's logs were made without it, so they carry the line "tabulate not installed"); `aep8`, `astropy`, `ppigrf`, `netCDF4` and `geant4_pybind` for the radiation chain in `radiation/` (the Geant4 and AP-8 results are tracked, so the Chapter 5 tables can be rebuilt without them); `cairosvg` to turn the bearing cross-section figure into a PNG; a LaTeX installation with `dvisvgm` for the equation images; `pytest` for the reference check.

## Reproduce the book

```bash
python run_all.py
```

This runs the 48 cases behind Chapters 1 to 9 and writes one log per case into `runs/` (`runs/ring_altitude_700.txt`, `runs/final_bearing.txt` and so on), with the figures and CSV files in subfolders. It took 3 hours 18 minutes on a laptop, with the solver caches in `data/` in place; the LIM deployment and the levitation bearing dynamics take most of it. `python run_all.py --list` prints the cases, `--only=name,name` runs a few, `--no-lim` skips the two deployment runs (about 20 minutes). The same list is in `run_700.sh` and `run_700_lim.sh` for bash.

## Subjects

Run each script from the repository root. Chapters are chapters of Volume IV; "log" is the file `run_all.py` writes.

### The ring

| Script | What it prints | Chapter |
|---|---|---|
| `ring_altitude.py` | the design point: radius, speeds, cable mass and tension, ring mass, deployment energy (log `ring_altitude_700.txt`; the prograde, 250 km and 800 km tables through the environment variables) | 1, 3, 8 (Table 8.1), and every chapter's ring numbers |
| `j2_simulation.py --inclination=20` | J2 precession of an inclined ring and the torque the anchor lines hold (logs `j2_retrograde_i*.txt`, `j2_prograde_i*.txt`) | 2 |
| `anchor_line_analysis.py --inclination=20` | anchor lines of an inclined ring: lateral force, line angle, length, stress, taper, number of lines per station (logs `anchor_700_i*.txt`) | 2 (Table 2.1), 15 |
| `lim_simulation.py` | deployment of the ring by linear induction motors: months of decelerating the cable, power, current and thermal limits (logs `lim_simulation.txt`, `lim_simulation_16MW.txt`, on the measured tape limit since 2026-10-26 and the flat winding with the Brandt-Indenbom hysteresis since 2026-10-27; `run_700_lim_tapelimit.sh` and `run_700_lim_brandt.sh` for the trades; manual in `docs/LIM_README.md`) | 3 (deployment time and energy), 10, 12, 13 |
| `debris_simulation.py` | where the released cable and the falling casing go after a break; drag lifetimes; minimum safe altitude (log `debris_simulation.txt`) | 4, 19 |
| `debris_orbits.py` | orbits of the released cable: apogee for each recoil case, the two halves (log `debris_orbits.txt`) | 4 |
| `altitude_trade.py` | released-cable orbits and drag lifetimes at several ring altitudes (log `altitude_trade.txt`) | 4, 19 |
| `radiation/` | the trapped-proton and cosmic-ray chain of Chapter 5: AP-8 spectra, the anomaly sector, Geant4 depth-dose, dose tables, Figure 5.1 (see `radiation/README_radiation.md`) | 5 |

### The levitation bearing

The bearing has four tracks of air-core superconducting coils: two support tracks under the cable that carry the casing, each with a driven compensation winding, and two lift tracks above the cable that pull against a passing spacecraft. There are no Halbach arrays and no iron on the casing.

| Script | What it prints | Chapter |
|---|---|---|
| `final_bearing.py` | the adopted bearing: track ratings, 114 dynamics cases (no single failure leads to contact, smallest gap 50 mm), losses, bundle and stray field, and the ring renumbered with the bearing's plates (log `final_bearing.txt`) | 8 (Sections 8.1 to 8.3, Tables 8.2 to 8.4) |
| `bearing_mass.py` | the mass of an iron-core bearing against the cable it costs, and the air-core alternative (log `bearing_mass.txt`) | 8 |
| `bearing_design.py` | iron-core field design: pole geometries solved with the 2-D solver, ratings, tape, coil lengths, protection (log `bearing_design.txt`) | 8 (the route that mass ruled out) |
| `bearing_fem_check.py` | the field solver checked against known answers and against scikit-fem (log `bearing_fem_check.txt`) | 8 |
| `compensation_fields.py` | the field cases of the compensation and support tracks; rewrites the cache `data/compensation_fields.json` (log `compensation_fields.txt`) | 8 |
| `sled_passage_simulation.py` | the load a passing spacecraft puts on the bearing: 12.6 g radial, the neutral speed, the casing's movement (log `sled_passage.txt`) | 8 (Section 8.3) |
| `magnetostatic_2d.py` | the two-dimensional field solver the bearing scripts import | 8 |
| `lsm_simulation.py` | the mass driver launch the bearing scripts take their speed profile from (see below) | 8, 14 |

| `bearing_optimize.py --part=1`, `--part=2 --stage=screen` (then `confirm25`, `confirm125`, `report`) | Part 1: the field solution checked for the present bearing (grid convergence, the tape limit from every tape of a modelled winding, an image-method check; two hours). Part 2: the support track search over loop width, center spacing, bundle cross-section, cryostat wall and loop count at 5 mm, confirmed at 2.5 and 1.25 mm (a day of computing; results cached in the output folder) | 8 |

### Launch and power

| Script | What it prints | Chapter |
|---|---|---|
| `lsm_simulation.py` | linear synchronous motor mass driver: the casing carries the stator coils and the sled carries DC superconducting magnets; laps, time, power and acceleration to the launch speed (log `lsm_simulation.txt`; manual in `docs/LSM_README.md`) | 14 |
| `md_simulation.py --no-graphs` | four-stage linear induction mass driver, eddy current model (log `md_simulation.txt`) | 14 |
| `md_trade_study.py` | the induction driver's payload against length, plate and thermal margin, with figures (log `md_trade_study.txt`) | 14 |
| `power_simulation.py` | solar generation on the ring and the suspended array, HVDC distribution, deployment and operations demand (logs `power_*.txt`) | 10, 16 |
| `tether_stress.py`, `freefall_trajectory.py`, `climber_descent.py`, `atmospheric_loading.py`, `rendezvous_analysis.py` | the ground-anchored tether: taper, free fall from release, descent from GEO, wind loading, rendezvous speeds (logs `tether_stress.txt`, `freefall.txt`, `climber_descent.txt`, `atmospheric_loading.txt`, `rendezvous.txt`; `tether_run_all.py` runs the five) | 17 |
| `elevator_system.py`, `radiation/elevator_dose_rates.py` | the elevator as one system: pods of 25 to 400 seats, their water shield from the dose along the tether, the descent with a 1.5 g stop, the graded CNT brake track, the Halbach magnets, the heat into the water, the tether with the pods and the track as loads, free fall after a loss of the tether, cargo drop pods (logs `elevator_system.txt`, `elevator_dose_rates.txt`; the dose script needs `aep8` and writes `elevator_dose_rates.npz` to the output folder, which the system script reads) | 17 |
| `elevator_system2.py` | the second pass on the elevator system: the track's thickness from the heat of a pass and its conductivity profile, the ribbon to 2,000 K with a shielded carriage, where the cruise energy goes (passive ribbon, regenerating stator, free fall then one stop), 4,000 m/s and the stop acceleration, the facing-pair and iron-core magnets, throughput against pod spacing, and the safety case (Allen-Eggers, a break below or above the pod, the impulse to a rescue orbit); imports `elevator_system.py` (log `elevator_system2.txt`) | 17 |
| `lunar/` | a mass driver on the Moon (its own configuration, physics and simulation; unchanged since the first edition). The Moon is covered in Volumes I and IV of *The Orbital Ring*. | 22 |

### Tools (reproduce a table or figure)

| Script | What it does | Chapter |
|---|---|---|
| `tools/ch4_numbers.py` | the numbers of the Chapter 4 tables, written to `runs/ch4_numbers.json` (log `ch4_numbers.txt`) | 4 (Tables 4.2 to 4.5) |
| `tools/ch4_extra.py` | Chapter 4 numbers the tables do not carry | 4 (Table 4.1) |
| `tools/ch4_figures.py` | Figures 4.1 to 4.6, as SVG and PNG in `runs/figures_ch04/` | 4 |
| `tools/fig03_02_deployment_tension.py` | Figure 3.2, the cable tension and the bearing load through the deployment, from the time series `lim_simulation.py` writes | 3 (Figure 3.2) |
| `tools/generate_equations_latex.py` | the worked equations with every number from `ring_altitude.py`, as LaTeX in `runs/equations/` (SVG when LaTeX is installed) | 1 to 4, 10 |
| `tools/generate_fig_bearing_cross_section.py` | Figure 8.1, the bearing cross section, in `runs/figures_ch08/` | 8 |
| `tools/ch10_numbers.py` (with `ORBITAL_RING_ARRAY_OFFSET_KM=0`), `tools/ch13_numbers.py` | the Chapter 10 power numbers and the Chapter 13 cable-direction rows | 10 (Table 10.1), 13 (Table 13.1) |
| `tools/lim_tape_current_check.py` | the first check of the tape current of the deployment LIM coil against the THEVA limit, by winding layout and layers (superseded by the next row) | 11 (Table 11.2), 12 |
| `lim_hysteresis.py` | the deployment LIM coil's winding as a set of tapes: the 2-D field each tape sees from the rest of the winding, its measured critical current at that field and angle, and the hysteresis loss on the Brandt-Indenbom strip model (the "brandt" option, the default since 2026-10-27); run alone for its checks | 11 (Sections 11.5 to 11.7), 18 (Table 18.2) |
| `tools/lim_hysteresis_model.py` | the three hysteresis models compared at the design point, the width, winding and turns trade on the Brandt model, and `--tables` for the Table 11.3 and 18.2 forms of the candidate runs (log `lim_hysteresis_model.txt`) | 11, 18 |
| `tools/ch18_thermal_series.py`, `tools/fig18_figures.py` | the thermal time series of the default deployment (cable and wall temperature, cold-side load, radiator width, the losses of a site) and Figures 18.1 and 18.2 drawn from it with the values Chapter 18 quotes | 18 (Figures 18.1, 18.2, Table 18.1) |
| `tools/lim_tape_limit.py` | the deployment LIM coil on the measured tape limit: the field at every tape of its own winding, the worst tape's critical current at its own field and angle (Robinson 70 K data), the layers, the coil current limit, the hysteresis by layer and width, the turns trade; `--tables` reads the deployment logs into the forms of Tables 11.3 and 18.2 (log `lim_tape_limit.txt`) | 11 (Table 11.2, Sections 11.5 to 11.7), 18 (Table 18.2) |
| `tools/fig11_01_traveling_wave.py`, `tools/fig11_01_load_angle.py` | Figures 11.1 and 11.2, in `runs/figures_ch11/` | 11 |
| `tools/numbers_700_hand.py` | hand calculations that check quoted numbers (Table 2.1, horizon, array shading) | 2, 3, 4, 15, 16 |
| `tools/number_sheet_hw10686.py` | the number sheet of the cable hardware change, from two folders of logs | the renumbering pass |
| `tools/lsm_chapter_metrics.py`, `tools/lsm_n_optimization_demo.py` | energies and cryogenic loads of the launch; why a fixed stator turn count is best (logs `lsm_chapter_metrics.txt`, `lsm_n_optimization.txt`) | 14 |
| `tools/make_bifacial_24hr_graph.py`, `tools/make_solar_hammock_plots.py` | the solar panel and suspended array figures, in `runs/graphs_power/` | 10, 16 |
| `tools/sled_passage_config_table.py` | the configuration table of the sled passage problem | review only |
| `tools/graph_plate_thickness_vs_thrust.py`, `tools/hts_coil_inductance_analysis.py` | thrust against plate thickness; HTS coil inductance with eddy currents | 12 (not named by the chapter) |

### Design history

Kept and runnable, labelled in their first line. The bearing they describe is not the adopted one.

| Script | What it is | Chapter |
|---|---|---|
| `levitation_simulation.py all`, `levitation_physics.py`, `levitation_config.py` | the Halbach limiter bearing that Chapter 8 considered and did not adopt (log `levitation.txt`) | 7, 8 |
| `compensation_tracks.py` | the six-track layout with Halbach and LIM limiters, also not adopted (log `compensation_tracks.txt`) | 8 |
| `sled_passage_limits.py`, `sled_passage_pulsed.py`, `sled_passage_spread.py` | the steps of the sled passage work that led to the adopted bearing (logs `sled_passage_limits.txt`, `sled_passage_pulsed.txt`, `sled_passage_spread.txt`) | 8 |

### Not yet mapped to a revised chapter

`gomory_lim_analysis.py`, `field_distribution.py`, `reentry.py`, `hts_coil_inductance_analysis.py` (in `tools/`), `gyro_simulation_v02.py`, `tools/lim_optimizer_alu.py`, `tools/lim_optimizer_gTiAl.py`. Chapters 14 onward have not yet been revised for the 700 km design; `docs/CHAPTER_MAP.md` maps Chapters 1 to 13.

## Environment variables

| Variable | Meaning | Default |
|---|---|---|
| `ORBITAL_RING_ALT_KM` | ring altitude in km | 700 |
| `ORBITAL_RING_CABLE` | `retrograde` or `prograde` | retrograde |
| `ORBITAL_RING_HW_KGM` | cable hardware mass per meter, kg/m | 10686 |
| `ORBITAL_RING_ARRAY_OFFSET_KM` | how far below the ring the solar array hangs, km (0 puts the panels on the ring) | 100 |
| `ORBITAL_RING_OUT` | the folder every script writes to | runs |

Every script that imports `ring_altitude.py` reads the first four; each script's docstring says which it uses. For example, `ORBITAL_RING_ALT_KM=800 python ring_altitude.py` prints the 800 km reference table, and `ORBITAL_RING_CABLE=prograde python j2_simulation.py --inclination=20` runs the prograde ring.

## Tests

```bash
pip install pytest
python -m pytest tests
```

About a minute. The test asserts the design point of `ring_altitude.py` to the precision the book quotes, runs the 22 fast simulations and compares their headline lines with `tests/reference/` (copied from the reference run), and checks the 13 slow simulations against the logs of your own `python run_all.py` run (skipped, with the command to make them, until that run exists). Any difference is a changed number: the reference lines carry no run-time stamps and no paths.

## Data and documentation

- `data/`: input files the scripts read, all tracked: the THEVA tape data (used with permission of THEVA), the three solver caches with their recompute commands, and three design-history logs. See `data/README.md`.
- `docs/LIM_README.md` and `docs/LSM_README.md`: the manuals of the deployment simulator and the mass driver, written for the first edition at 250 km.
- `docs/CHAPTER_MAP.md`: for Chapters 1 to 13, the script and log behind each table and figure.
- `radiation/README_radiation.md`: the Chapter 5 radiation chain.

## License

MIT License, see `LICENSE`.

## Citation

Paul G. de Jong, *The Orbital Ring, Volume IV: Space Mass Transit Complete Design Review*. https://orbitalring.space and https://github.com/kjpaul/orbitalring
