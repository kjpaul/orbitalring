# Chapter map: which script and log is behind each table and figure

*The Orbital Ring, Volume IV: Space Mass Transit Complete Design Review*, Chapters 1 to 9. Logs are the files `python run_all.py` writes into the output folder (`runs/` by default); the equation images come from `tools/generate_equations_latex.py`, which takes every number from `ring_altitude.py`. Figures marked "drawn" are diagrams made by hand, not by a script. Chapters 10 onward have not yet been revised for the 700 km design and are not mapped here.

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
| Section 3.4, deployment time and energy (599.7 days at 8 MW per site, 461.0 at 16 MW, 18.8 EJ) | `lim_simulation.py`; `ring_altitude.py` for the energy | `lim_simulation.txt`, `lim_simulation_16MW.txt`, `ring_altitude_700.txt` |
| Figure 3.1, hoop tension force cancellation | drawn | |

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

## Chapters 10 onward

Not yet revised for the 700 km design, so not mapped. The scripts that will serve them are listed in the README under "Launch and power" (Chapters 10, 12, 14, 16, 17), "Tools" and "Not yet mapped to a revised chapter". The lunar mass driver in `lunar/` belongs to Chapter 22, Rings at Other Worlds and the Portable Ring.
