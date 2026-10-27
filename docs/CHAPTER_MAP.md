# Chapter map: which script and log is behind each table and figure

*The Orbital Ring, Volume IV: Space Mass Transit Complete Design Review*, Chapters 1 to 13. Logs are the files `python run_all.py` writes into the output folder (`runs/` by default); the equation images come from `tools/generate_equations_latex.py`, which takes every number from `ring_altitude.py`. Figures marked "drawn" are diagrams made by hand, not by a script. Chapters 14 onward have not yet been revised for the 700 km design and are not mapped here.

## Chapter 1: The Ring and Its Orbit

| Item | Script | Log or file |
|---|---|---|
| Table 1.1, orbital parameters at 700 km with 250 km for comparison | `ring_altitude.py` (`ORBITAL_RING_ALT_KM=250` for the second column) | `ring_altitude_700.txt`, `ring_altitude_250.txt` |
| Equations 1.3 to 1.8 (orbital velocity, gravity, ground-synchronous velocity, net acceleration, circumference at 700 km) | `tools/generate_equations_latex.py` | `equations/eq03_vorbit_700.tex` to `eq08_circumference_700.tex` |
| Equations 1.1, 1.2, 1.6 (sidereal day, angular velocity, equatorial surface speed) | constants, no script | |

## Chapter 2: Equatorial versus Inclined Orbits

| Item | Script | Log or file |
|---|---|---|
| Equation 2.1, J2 precession at 700 km and 30 degrees | `j2_simulation.py --inclination=30` | `j2_retrograde_i30.txt` |
| Table 2.1, lateral Earth-tracking forces by inclination | `anchor_line_analysis.py --inclination=5, 10, 20, 30`; hand check in `tools/numbers_700_hand.py` | `anchor_700_i5.txt` to `anchor_700_i30.txt` |
| Equation 2.2, lateral kinetic energy of the ring | `tools/generate_equations_latex.py` | `equations/eq20_lateral_KE.tex` |
| Figures 2.1 to 2.5 and 2.7 (great circle against latitude line, rotation radius, Earth spinning inside the ring, apex and node velocities, staged atmospheric braking) | drawn | `images/figures/fig02_01` to `fig02_07` in the book |
| Figure 2.6, relative ground velocity around the ring at 30 degrees | no script in this repository found for this plot [UNVERIFIED] | |
| Casing velocity at the apex and the nodes (447 m/s at 30 degrees) | `j2_simulation.py` | `j2_retrograde_i30.txt` |

## Chapter 3: The Cable

| Item | Script | Log or file |
|---|---|---|
| Table 3.1, prograde against retrograde cable at 700 km | `ring_altitude.py` (`ORBITAL_RING_CABLE=prograde` for the second column) | `ring_altitude_700.txt`, `ring_altitude_700_prograde.txt` |
| Equations 3.6 to 3.10 (delta-v, cable velocity, tension, second moment of area, bending length) | `tools/generate_equations_latex.py` | `equations/eq09_dv_retro_eq.tex`, `eq15_vcable_retro_eq.tex`, `eq16_tension_retro_eq.tex`, `eq17_second_moment.tex`, `eq18_Lchar.tex` |
| Equations 3.1 to 3.5 (momentum conservation, tension, stress, the sizing quadratic) | symbolic, no script; the quadratic is `ring_altitude.cable_mass_structural` | |
| Section 3.4, deployment time and energy (499.7 days at 8 MW per site, 427.3 at 16 MW, 18.8 EJ, on the flat winding at the measured tape limit with the Brandt-Indenbom hysteresis, the default since 2026-10-27; 597.8 and 424.6 on the 2026-10-26 pancake with the book's loss factor; 599.7 and 461.0 on the old 163 A cap) | `lim_simulation.py`; `ring_altitude.py` for the energy | `lim_simulation.txt`, `lim_simulation_16MW.txt` (the same runs as `lim_simulation_brandt_flat3mm_N208_8MW.txt`, `..._16MW.txt`), `ring_altitude_700.txt` |
| Figure 3.1, hoop tension force cancellation | drawn | |
| Figure 3.2, cable tension, stress, casing net weight and cable excess force through the deployment; the deployment tension table of Section 3.4 (0, 13.1, 118.3, 328.6, 740.3, 845.6 GN at casing slowdowns of 0 to 8,020 m/s) | `lim_simulation.py` (the "DEPLOYMENT TENSION" section and `lim_deployment_tension.csv`), drawn by `tools/fig03_02_deployment_tension.py` from the default run's series (since 2026-10-27 the flat winding on the Brandt model, `deployment_tension_700_brandt_flat.csv`, 499.7 days, the casing at rest on day 433); `ORBITAL_RING_CABLE=prograde` for the prograde curve | `lim_simulation.txt`, `figures_ch03/fig03_02_deployment_tension_700.svg` |

## Chapter 4: Why 700 Kilometers? Catastrophic Failure and the Design Altitude

| Item | Script | Log or file |
|---|---|---|
| Table 4.1, collision energy against casing mass participation | `tools/ch4_extra.py` | printed to the console |
| Equations 4.3 to 4.5, 4.8, 4.9, 4.11 (centrifugal and net acceleration, gap crossing time, apogee at 700 km, lifetime estimate, critical delta-v at 700 km) | `tools/generate_equations_latex.py` | `equations/ch2_eq03_a_centrifugal.tex` to `ch2_eq05_t_gap.tex`, `eq04_08_apogee_700.tex`, `eq04_09_lifetime.tex`, `ch2_eq09_dv_crit_700.tex` |
| Tables 4.2 to 4.5 (released cable's orbit by altitude, fragment lifetimes, the 250 km cases, prompt soot by altitude and beta) | `tools/ch4_numbers.py`, built on `debris_orbits.py` and `debris_simulation.py` | `ch4_numbers.txt`, `ch4_numbers.json` |
| Table 4.6, 250 km against 700 km | `ring_altitude.py` at both altitudes, `debris_simulation.py`, `tools/ch4_numbers.py` | `ring_altitude_250.txt`, `ring_altitude_700.txt`, `debris_simulation.txt`, `ch4_numbers.txt` |
| Apogee of the released cable (4,565 km; 2,131 and 8,113 km for the two halves) | `debris_orbits.py` | `debris_orbits.txt` |
| Drag lifetimes by altitude (Section 4.6) | `altitude_trade.py` (needs pymsis) | `altitude_trade.txt` |
| Minimum safe altitude and the ejecta populations | `debris_simulation.py` | `debris_simulation.txt` |
| Figure 4.1 (orbit of the released cable, `fig04_kepler_orbit`), Figures 4.2 to 4.6 (`fig04_01_released_cable_orbits` to `fig04_05_ring_mass_vs_altitude`, with Figure 4.5 from `fig04_04_proton_free_share`) and Figure 4.7 (`fig04_06_design_space`) | `tools/ch4_figures.py` (Figure 4.5 reads `radiation/igrf_altscan_clean.txt`) | `figures_ch04/*.svg` and `.png` |

## Chapter 5: The Radiation Environment at 700 km

All scripts run from the `radiation/` folder and write there; see `radiation/README_radiation.md` for the order and the packages (the Geant4 and AP-8 results are tracked, so the tables can be rebuilt without them).

| Item | Script | Log or file |
|---|---|---|
| Table 5.1, the anomaly sector on the equator as of 2026 | `radiation/igrf_sector.py`, `radiation/igrf_altscan.py` | `igrf_sector.txt`, `igrf_altscan_clean.txt` |
| Table 5.2, trapped proton flux on the equator by altitude; Figure 5.1, flux around the ring | `radiation/ch5_figure.py` | `fig05_01_proton_flux_around_ring.svg` and the rows it prints |
| Table 5.3, range of protons in water, aluminium and carbon | NIST PSTAR tables `radiation/pstar_*.txt`, read by `radiation/response.py` | |
| Table 5.4, dose from trapped protons behind a water shield | `radiation/tables_trapped.py` (`fold_trapped.py`, Geant4 runs in `out_trap/`) | `tables_trapped_igrf2026.txt`, `tables_trapped_epoch.txt` |
| Table 5.5, solar cell degradation; Table 5.6, dose in silicon behind aluminium | `radiation/solar_cal.py`, `radiation/ch5_tables.py` | `solar_cal.txt`, `ch5_tables.txt` |
| Table 5.7, cosmic ray dose behind a water shield | `radiation/gcr.py` | `gcr_per_sr.txt` |
| Tables 5.8 and 5.9, dose under the cable (trapped protons, cosmic rays) | `radiation/sweetspot_trapped.py`, `radiation/sweetspot_run.py`, `radiation/sweetspot_scan.py` | `sweetspot_trapped.txt`, `sweetspot_gcr.txt` |
| Table 5.10, dose limits and reference values | from the sources cited in the chapter, no script | |
| Equation 5.1, the geomagnetic cutoff rigidity at 700 km | not in `tools/generate_equations_latex.py`; a hand-made image [UNVERIFIED] | |
| Low-energy protons and atomic oxygen checks (Section 5.10) | `radiation/ch5_checks.py`, `radiation/lowE_protons.py` | printed to the console |
| Figure 5.2, the sky seen from the room under the cable (shares of the sky from the sweet-spot logs) | `radiation/fig05_02_sky_under_cable.py` (reads `sweetspot_gcr.txt`, `sweetspot_trapped.txt`) | `fig05_02_sky_under_the_cable.svg` |
| Figure 5.3, dose against shield thickness, trapped protons and cosmic rays | `radiation/fig05_03_dose_vs_shield.py` (reads `tables_trapped_igrf2026.txt`, `sweetspot_gcr.txt`) | `fig05_03_dose_vs_shield.svg` |
| Figure 5.4, the anomaly sector on the equator as a strip map | `radiation/fig05_04_anomaly_sector_map.py` (reads `ap8_igrf2026_spectra.npz`; land from the Wikipedia "Equator" crossing table) | `fig05_04_anomaly_sector_map.svg` |
| Solar particle event dose (Section 5.9): reference events, quiet and storm-time cutoff, dose per event behind 5, 20, 100 cm of water and under the cable | `radiation/spe_700.py` (King 1974 and Koldobskiy et al. 2021 spectra, Stormer cutoff of `sweetspot.py`, Geant4 proton responses) | `spe_700.txt` |

## Chapter 6: Core Superconducting Technologies

| Item | Script | Log or file |
|---|---|---|
| Table 6.1, layers of the Theva TPL4000 tape | from the manufacturer's data, no script | |
| Table 6.2, critical current of a 12 mm Theva AP tape in a magnetic field | the 0 degree column of `data/THEVA_Pro-Line_Advanced_Pinning_2G_HTS_70K_Field_Dependence.csv` (the file `bearing_design.read_theva` reads), no script | |
| Table 6.3, solenoid and Halbach geometries on the ring | no script found [UNVERIFIED]; the chapter is a technology survey | |
| Figures 6.1 to 6.4 | drawings and photographs | |

## Chapter 7: Electromagnetic Levitation, EML vs EDL

| Item | Script | Log or file |
|---|---|---|
| Table 7.1, the three families of magnetic levitation | survey, no script | |
| Figures 7.1 and 7.2 (`fig05_01_eml_vs_edl_mechanism`, `fig05_02_edl_lift_drag`) | no script in this repository found for these two figures [UNVERIFIED]; `levitation_simulation.py` (design history) draws the related `fig_levitation_*` figures | `graphs_levitation/` |

## Chapter 8: The Levitation System

The adopted bearing: four tracks of air-core superconducting coils (two support tracks under the cable with driven compensation windings, two lift tracks above it), four plates on the cable, no iron on the casing.

| Item | Script | Log or file |
|---|---|---|
| Table 8.1, reference case for the levitation system | `ring_altitude.py`; `final_bearing.py` Part A (the ring renumbered with the bearing's plates) | `ring_altitude_700.txt`, `final_bearing.txt` |
| Table 8.2, tape current in a support track against the corner limit | `final_bearing.py` Part D.1 (bundle shape, the corner field) | `final_bearing.txt` |
| Table 8.3, pull of a 1 m air-core loop on a 200 mm plate and the turns per bundle | `final_bearing.py` Part A, from the air-core cache `bearing_mass.py` makes (`data/bearing_air.json`) | `final_bearing.txt`, `bearing_mass.txt` |
| Table 8.4, compensation pull and casing movement during a launch | `final_bearing.py` Parts B and F (the dynamics cases) | `final_bearing.txt` |
| The 114 cases, no single failure leads to contact, smallest gap 50 mm | `final_bearing.py` Part B.6 | `final_bearing.txt` |
| The field solver and its check (named at Section 8.1) | `magnetostatic_2d.py`, `bearing_fem_check.py` | `bearing_fem_check.txt` |
| The mass argument against iron (named at Section 8.1) | `bearing_mass.py` | `bearing_mass.txt` |
| Section 8.3, the sled passage: 12.6 g radial load, the neutral speed, the casing's momentum | `sled_passage_simulation.py` | `sled_passage.txt` |
| The rejected options (Halbach limiters, LIM limiters, 200 m coils) | `compensation_tracks.py`, `sled_passage_limits.py` (design history) | `compensation_tracks.txt`, `sled_passage_limits.txt` |
| Figure 8.1, the bearing cross section | `tools/generate_fig_bearing_cross_section.py` | `figures_ch08/fig10_bearing_cross_section.svg` |

## Chapter 9: AC, DC, and Three-Phase Currents

| Item | Script | Log or file |
|---|---|---|
| Table 9.1, where DC lives on the ring; Figures 9.1 and 9.2 | survey and drawings, no script | |

## Chapter 10: Powering the System

| Item | Script | Log or file |
|---|---|---|
| Table 10.1, power budget summary (ring-mounted panels for deployment, demand during operations) | `power_simulation.py` with `ORBITAL_RING_ARRAY_OFFSET_KM=0`, `--power=16` for the 16 MW case, `--width=` for the panel width (the chapter uses 55 m); `--demand=ops` for operations. The chapter's numbers come from `ORBITAL_RING_ARRAY_OFFSET_KM=0 python tools/ch10_numbers.py` (built on `power_config.py` and `power_physics.py`). | `power_deployment_ring_mounted_8MW.txt`, `power_deployment_ring_mounted_16MW.txt`, `power_ops.txt` |
| Figure 10.1, power around the ring during deployment (`fig10_01_ring_power_profile_700`) | `power_simulation.py --save ring_profile` with `ORBITAL_RING_ARRAY_OFFSET_KM=0` | `graphs_power/01-ring_profile.png` |
| Equations 10.1 and 10.2, deployment power and expanded power | `tools/generate_equations_latex.py` | `equations/eq10_01_deploy_power_700.tex`, `eq10_02_expanded_power_700.tex` |

## Chapter 11: Linear Induction and Linear Synchronous Motors

| Item | Script | Log or file |
|---|---|---|
| Table 11.1, choosing between a LIM and an LSM | the comparison draws on the three motor runs below; no single script | `lim_simulation.txt`, `lsm_simulation.txt`, `md_simulation.txt` |
| Table 11.2, LIM parameters for deployment | `lim_simulation.py` (parameters in `lim_config.py`, physics in `lim_physics.py`); the 16 MW column from the second run | `lim_simulation.txt`, `lim_simulation_16MW.txt` |
| Table 11.2, the tape current rows (layers of 3 mm tape, the limit in the winding's field, the coil current limit) | `tools/lim_tape_limit.py` (the field at every tape of the coil's own winding, the Robinson 70 K angle data of `data/`; the first check `tools/lim_tape_current_check.py` is superseded) | `lim_tape_limit.txt` |
| Table 11.3, deployment results on the measured tape limit (8 and 16 MW per site); Section 11.7, the turns trade (100, 150, 208, 300 turns) | `lim_simulation.py` with `lim_config.TAPE_LIMIT_MODE = "measured"` (the default since 2026-10-26); `run_700_lim_tapelimit.sh`; `tools/lim_tape_limit.py --tables` for the table forms | `lim_simulation_tapelimit.txt`, `lim_simulation_tapelimit_16MW.txt`, `lim_simulation_tapelimit_N100_8MW.txt` to `..._N300_16MW.txt`; report `lim_tape_limit_report_2026-10-26.md` (its pancake figures superseded by the next row) |
| The stator hysteresis on the Brandt-Indenbom strip model, the width, winding and turns trade on it, and the deployment of the candidates (Sections 11.5 to 11.7, Table 11.3; the flat winding and the Brandt model are the defaults since 2026-10-27) | `lim_hysteresis.py` (the field at every tape, the loss, its checks), `tools/lim_hysteresis_model.py` (comparison, trade; `--tables` for the table forms), `run_700_lim_brandt.sh` | `lim_hysteresis_model.txt`, `lim_hysteresis_model_tables.txt`, `lim_simulation_brandt_flat3mm_N208_8MW.txt` and the other seven candidates; report `lim_hysteresis_model_report_2026-10-27.md` |
| Table 11.3, mass driver LSM parameters | `lsm_simulation.py` (parameters in `lsm_config.py`); `tools/lsm_chapter_metrics.py` for the energies and cryogenic loads | `lsm_simulation.txt`, `lsm_chapter_metrics.txt` |
| Figure 11.1, how three phases make a traveling wave; Figure 11.2, the load angle of a linear synchronous motor | `tools/fig11_01_traveling_wave.py`, `tools/fig11_01_load_angle.py` (no ring numbers in them) | `figures_ch11/fig11_01_traveling_wave.png`, `figures_ch11/fig11_02_load_angle.png` (and .svg) |
| Equations 11.1 to 11.3, slip frequency, LSM thrust, back EMF (symbolic) | `tools/generate_equations_latex.py` | `equations/eq11_01_slip_frequency.tex` to `eq11_03_back_emf.tex` |

## Chapter 12: The Linear Induction Motors

| Item | Script | Log or file |
|---|---|---|
| The deployment at 8 MW and 16 MW per site: time, current, slip, plate heating, the voltage limit (Sections 12.1 to 12.5) | `lim_simulation.py`, the scripts the chapter names (`lim_config.py`, `lim_physics.py`); manual in `docs/LIM_README.md` | `lim_simulation.txt`, `lim_simulation_16MW.txt` |
| The tape current | `tools/lim_tape_limit.py`, as Table 11.2 above (the 2-D field of the winding and the measured angle data) | `lim_tape_limit.txt` |
| Equations `eq26_V_lim` and `eq27_V_N2` (the LIM voltage) | images kept from the first edition, not in the generator | |

## Chapter 13: The Retrograde Cable Question

| Item | Script | Log or file |
|---|---|---|
| Table 13.1, prograde and retrograde cables at 700 km as the subsystems see them | `ring_altitude.py` for both directions (`ORBITAL_RING_CABLE=prograde` for the second); the chapter's table rows come from `python tools/ch13_numbers.py`, run once for each direction | `ring_altitude_700.txt`, `ring_altitude_700_prograde.txt` |
| The inclined-ring consequences of the cable direction | `j2_simulation.py --cable=prograde` | `j2_prograde_i5.txt` to `j2_prograde_i30.txt` |
| The deployment with a prograde cable | `ORBITAL_RING_CABLE=prograde python lim_simulation.py` (not in the run set; run on demand) | |

## Chapter 18: Cryogenics

| Item | Script | Log or file |
|---|---|---|
| Figures 18.1 (cable and wall temperature) and 18.2 (cold-side load and radiator width) through the deployment; Table 18.1's eddy loss per site | `tools/ch18_thermal_series.py` (reruns the default deployment and writes the thermal series), `tools/fig18_figures.py` (the two figures and the values the chapter quotes) | `ch18_thermal_series.csv`, `ch18_thermal_series_run.txt`, `fig18_values.txt` |
| Table 18.2, the cold loop at one LIM site (hysteresis in the two LIMs, the cryostat leak, the levitation coils, the compressor power and the radiator) | from the deployment run's peak and its months 12 to 15 with the constants of `lim_config.py` and `lim_physics.py`; `tools/lim_hysteresis_model.py --tables` prints the table's form from the Brandt-model logs (`tools/lim_tape_limit.py --tables` from the 2026-10-26 ones) | `lim_simulation_brandt_flat3mm_N208_8MW.txt` (the default run), `lim_hysteresis_model_tables.txt` |

## Chapters 14 onward

Not yet revised for the 700 km design, so not mapped (Chapter 18's Table 18.2 above is the exception). The scripts that will serve them are listed in the README under "Launch and power" (Chapters 14, 16, 17), "Tools" and "Not yet mapped to a revised chapter". The lunar mass driver in `lunar/` belongs to Chapter 22, Rings at Other Worlds and the Portable Ring.
