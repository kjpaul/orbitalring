# The numbers sheet for the rising-quadratic design point (2026-11-07)

Asked for in `Vol_III/Claude_Code_instructions_numbers_sheet_2026-11-06.md` and part 4 of `Claude_Code_instructions_slip_rule_and_sheet_2026-11-07.md`. The design point: `ring_altitude.py` with the rising quadratic (round 4, commit 31e364d): structural cable 92,883 kg/m, 103,569 with hardware, 8,330.8 m/s (8,847 past the casing), 683.0 GN, 54.64 m^2, 7.39 m side, 4,606 Mt of cable, 5,140 Mt ring; prograde 71,386 / 82,072 kg/m, 8,420.5 m/s (7,904), 524.9 GN, 41.99 m^2, 6.48 m. The deployment design since round 5 and 6a: two LIMs per site at 12 MW with 150 turns, 381.1 days (Sections 14 and 15 of `deployment_rising_radius_2026-11.md`). Every number below comes from a script run under the two quadratics (`ORBITAL_RING_QUADRATIC=fixed` for the October value, the default `rising` for the new one), or from `tools/design_point_numbers.py`, which does the arithmetic the book did by hand from `ring_altitude.py`'s values ("arithmetic" in the script column). "Book now" is what the chapter file says on 2026-11-07, since Fable has re-keyed some of it already: "new" means the chapter already carries the new value, "October" that it still carries the old one.

The logs: October in `Vol_III/runs_700km_hw10686_2026-10-12/` (and `runs/design_point_sheet/fixed/` for the tools run under the fixed flag today), new in `runs/rising_quadratic/all/` and `runs/design_point_sheet/rising/`; the deployment runs in `runs/rising_radius_design/`.

**Hard-coded October values found and fixed (all read `ring_altitude.py` now):** `debris_simulation.ring_state` computed the cable with its own copy of the fixed-radius quadratic and momentum sum (through `lim_physics.calc_cable_mass`), so `debris_orbits.py`, `altitude_trade.py`, `debris_simulation.py`, `tools/ch4_numbers.py`, `tools/ch4_extra.py` and the Chapter 4 figures did not follow the design point; it now calls `ring.design()` (the October log reproduces under the fixed flag to the line). `radiation/sweetspot.py` had the cable's side as 8.20 m and `radiation/fig05_04_sky_under_cable.py` as 8.22 m; both take `ring.CABLE_SIDE` (the October 8.20 is kept under the fixed flag so the October dose files reproduce). `tools/ch03_direction_numbers.py` told the cable's direction from `DELTA_V_LOAD`, which fails under the rising quadratic (it printed the prograde cable with the retrograde relative speed); it reads `ring.RETROGRADE`. `tools/ch08_numbers.py` had 8,786.2 m/s typed into its LIM frequency line and no 12 MW case; both fixed. `casing_6000.py` carries the October check values in a label only.

## Chapter 1

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| Table 1.1, cable velocity | 8,270 m/s westward, 8,786 m/s relative | 8,331 m/s, 8,847 m/s | new | `ring_altitude.py` |
| Table 1.1, cable mass | 115,005 structural, 125,691 with hardware | 92,883; 103,569 kg/m | new | `ring_altitude.py` |
| Table 1.1, cable cross-section | 67.7 m^2, 8.22 m | 54.6 m^2, 7.39 m | new | `ring_altitude.py` |
| Table 1.1, cable tension | 846 GN | 683 GN | new | `ring_altitude.py` |
| Table 1.1, deployment | 647 days at 8 MW, 712 GW, 23.9 EJ (October cable, rising radius) | 381 days at 12 MW, 1,067 GW, 19.2 EJ | new | `lim_simulation.py --rising-radius` (Section 15) |
| Figure 1.1, the cross section | "8.22 m square", "27.5 m a side at 8 MW" | "7.39 m square", "41.25 m a side at 12 MW" | redrawn | the SVG's text (the figure is drawn, not scripted) |

## Chapter 3

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| Table 3.1, structural cable, retro / pro | 115,005 / 89,266 kg/m | 92,883 / 71,386 | new | `ring_altitude.py` |
| Table 3.1, cable with hardware | 125,691 / 99,952 | 103,569 / 82,072 | new | `ring_altitude.py` |
| Table 3.1, cable velocity | 8,270 / 8,343 m/s | 8,331 / 8,420 | new | `ring_altitude.py` |
| Table 3.1, relative velocity | 8,786 / 7,827 | 8,847 / 7,904 | new | `ring_altitude.py` |
| Table 3.1, cross-section, side | 67.65 m^2, 8.22 m / 52.51, 7.25 | 54.6, 7.39 / 42.0, 6.48 | new | `ring_altitude.py` |
| Table 3.1, tension | 846 / 656 GN | 683 / 525 | new | `ring_altitude.py` |
| Table 3.1, ring mass per meter and total | 137,691 kg/m, 6,124 Mt / 111,952, 4,979 | 115,569, 5,140 / 94,072, 4,184 | new | `ring_altitude.py` |
| Table 3.1, mass ratio structural / casing | 9.58 / 7.44 | 7.74 / 5.95 | new | arithmetic |
| Table 3.1, wave speed sqrt(T/m), Mach of the cable past the casing | 2,594 / 2,563 m/s, 3.39 / 3.05 | 2,568 / 2,529 m/s, 3.45 / 3.53 | new | `tools/ch03_direction_numbers.py` |
| Table 3.1, LIM supply frequency | 43.9 / 39.1 Hz | 44.2 / 39.5 Hz | new (44 / 40) | `tools/ch03_direction_numbers.py` |
| Equations 3.6, 3.7 (the momentum constraint) | 8,020.4 m/s; 7,504.3 + 12,000 x 8,020.4 / 125,691 = 8,270.0 | 7,928.4 m/s; 7,412.2 + 12,000 x 7,928.4 / 103,569 = 8,330.8 (and Eq 3.5b, v_mom) | new images | `tools/generate_equations_latex.py` |
| Equation 3.8, tension | 12.5e9 x 67.65 = 845.6 GN | 12.5e9 x 54.64 = 683.0 GN | new image | the same |
| Equation 3.9, second moment; the bullet "I is ... 381 m^4" | 8.22^4 / 12 = 381 m^4 | 7.39^4 / 12 = 249 m^4 | image new; bullet October (381) | the same; arithmetic |
| Equation 3.10 and its bullets, L_char; "T is the cable tension, 846 GN" | sqrt(500e9 x 381 / 845.6e9) = 15.0 m; 846 GN | sqrt(500e9 x 249 / 683.0e9) = 13.5 m; 683 GN | image and the text's 13.5 m new; the bullet October (846 GN) | the same; `sled_passage.py` prints the static L_char 13.5 m |
| 3.2, the hoop stress with no casing load, retro / pro | 22.4 / 25.2 GPa | 24.8 / 28.5 GPa | check | arithmetic, (m_c/A)(v_c^2 - v_orb^2) |
| 3.2, the added load at which the hoop stress is zero, retro / pro | 15,000 / 12,000 kg/m (15,088 / 11,711) | 12,200 / 9,400 kg/m (12,185 / 9,365) | check (Chapter 16 L72 says "about 15,000 kg/m") | arithmetic |
| 3.2 and 9, the deployment energy | 18.8 EJ (fixed radius) | 19.18 EJ (the rising closed form: kinetic 11.57, lift 7.24, elastic 0.37) | new | `ring_altitude.py` |
| 3.2, the build orbit | none | 527.36 km, 7,597.5 m/s, v_mom 7,412.2 m/s, rise 172.64 km, strain 2.5 % | new | `ring_altitude.py` |
| Westward launch band (Chapter 3 L19, Chapter 10): subsonic between | 6,200 and 11,400 m/s (6,192 to 11,380) | 6,279 to 11,415 m/s | check | `tools/ch03_direction_numbers.py` |
| Prograde cable, the same band | 5,300 to 10,400 m/s (5,265 to 10,390) | 5,375 to 10,433 m/s | check | the same, `ORBITAL_RING_CABLE=prograde` |

## Chapter 4

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 4.1, v_rel and the collision | 8,786 m/s | 8,847 m/s | new | `tools/ch4_extra.py` |
| Table 4.2, reduced mass and energy, whole casing | 10,954 kg/m, 423 GJ/m, 35.2 MJ/kg | 10,754 kg/m, 421 GJ/m, 35.1 MJ/kg | check | arithmetic, `tools/ch4_extra.py` (Table 4.1 block) |
| Table 4.2, half casing | 5,727 kg/m, 221 GJ/m, 36.8 MJ/kg (book: 1.68 MJ/kg of casing share) | 5,671, 222, 37.0 (2.03 MJ/kg) | October (5,727, 221, 1.68) | the same |
| Table 4.2, a tenth | 1,189 kg/m, 46 GJ/m | 1,186, 46 | check | the same |
| Table 4.2, ring total | 1.880e19 J = 4.49 Gt TNT | 1.872e19 J = 4.47 Gt | check | `tools/ch4_extra.py` |
| 4.2, the gap time and tangential run | 0.1439 s, 1,264 m, 1.39 m/s radial | 0.1428 s, 1,264 m, 1.40 m/s | check | `tools/ch4_extra.py` |
| 4.2, the outward force on the roof | 214.5 kN/m = 95.0 + 119.5; 26.1 kPa | 191.5 kN/m = 95.0 + 96.5; 25.9 kPa | check | `tools/ch4_extra.py` |
| 4.2, Mach wedge, fragment share of the side | 27.1 deg; 0.28 % | 26.9 deg; 0.31 % | check | `tools/ch4_extra.py` |
| 4.2, ejecta: 0.3 % and 5 % of the cable | 377 and 6,285 kg/m; +1,200 = 7,485 (5.4 % of the ring); +3,000 = 9,285 (6.7 %) | 311 and 5,178 kg/m; 6,378 (5.5 %); 8,178 (7.1 %) | check | `tools/ch4_extra.py` |
| 4.3, "N m/s faster than a free satellite" (v_cable - v_orbit) | 766 m/s | 827 m/s | check | `tools/ch4_extra.py` ("excess 826.5") |
| Figure 4.1, Equation 4.8: apogee of the released cable | 4,565 km altitude (radius 10,943 km) at 8,270 m/s | 4,986 km (radius 11,364 km) at 8,331 m/s | new | `debris_orbits.py`, `tools/ch4_numbers.py` |
| Table 4.4, no recoil | 8,270 m/s, 700 x 4,565 km | 8,331, 700 x 4,986 | October | `tools/ch4_extra.py` |
| Table 4.4, fast half +429 | 8,699 m/s, 700 x 8,113 km | 8,760, 700 x 8,751 | October | the same |
| Table 4.4, slow half -429 | 7,841 m/s, 700 x 2,131 km | 7,902, 700 x 2,429 | October | the same |
| Table 4.4, slow half with the radial kick (the floor) | 7,841 m/s, 429 m/s, 587 x 2,301 km | 7,902, 429, 602 x 2,585 km | October | the same |
| Table 4.4, slow half + 1,200 kg/m roof; with the kick | 7,762 m/s, 700 x 1,760; 560 x 1,960 | 7,806, 700 x 1,963; 576 x 2,142 | October | the same |
| Table 4.4, slow half + 3,000 kg/m roof; with the kick | 7,646 m/s, 700 x 1,260; 492 x 1,520 | 7,665, 700 x 1,341; 506 x 1,585 | October | the same |
| 4.4, "N m/s to lose" (slowdown to a 500 km perigee) | 820 m/s | 880 m/s | check | `tools/ch4_extra.py` |
| 4.4, the casing the slow half must sweep up to reach 500 km | 6,170 kg/m, 51 % | 5,873 kg/m, 49 % | check | `tools/ch4_extra.py` |
| 4.4, "between about 590 km and 8,300 km", "8,100 km" | 587 km to 8,113 km | 602 km to 8,751 km | October | `tools/ch4_extra.py` |
| Table 4.5 (and Figure 4.2), the released cable against ring altitude: 700 km row | 8,270; 7,504; 766; 4,565; 2,131 to 8,113; 587; 820 | 8,331; 7,504; 827; 4,986; 2,429 to 8,751; 602; 880 | check | `tools/ch4_extra.py` (every altitude row in the log) |
| Table 4.5, 250 km row | 8,506; 751; 3,631; 1,481 to 6,707; 144 | 8,571; 816; 4,023; 1,762 to 7,289; 159 | check | the same |
| Table 4.6, the full-section ballistic coefficient | 6,400 kg/m^2 (6,356) with the 8.22 m section | 5,700 kg/m^2 (5,712) with 7.39 m | October | `tools/ch4_numbers.py` ("Brod") |
| Table 4.6 (and Figure 4.3), full-section lifetimes: 150 / 200 / 250 / 300 / 400 / 500 / 600 / 700 / 800 km | 10 d / 133 d / 1.9 y / 6.9 y / 55 y / 330 y / 1,700 y / 7,400 y / 27,000 y | 9 d / 120 d / 1.7 y / 6.2 y / 49 y / 300 y / 1,500 y / 6,700 y / 24,000 y (0.0246, 0.328, 1.70, 6.21, 49.3, 296, 1,510, 6,685, 24,272 years) | October | `tools/ch4_numbers.py` ("circ_rod"; the 0.8 m and 8 cm columns do not change) |
| Table 4.7, 700 km slow half + kick: orbit, lifetimes | 587 x 2,300 km; 39,400 y / 3,700 / 370 (log: 65,950 / 6,265 / 627 cycle-averaged by the tool's other rule) | 602 x 2,585 km; 88,194 / 9,306 / 931 (the tool's rule) | October; the book's 39,400 does not match either log, see the note | `tools/ch4_numbers.py` ("life_w_*") |
| Table 4.7, + 1,200 kg/m roof; + 3,000 kg/m roof | 560 x 1,960; 492 x 1,520 km | 576 x 2,142; 506 x 1,585 km | October | `tools/ch4_extra.py`; lifetimes to be rerun with `debris_simulation.py` for those orbits (the tool prints the slow-half-plus-kick case only) |
| Table 4.8 (Figure 4.4), prompt soot, 700 km row | 7.8 / 5.6 / 3.4 / 1.5 / 0.3 Mt | 7.1 / 5.1 / 3.1 / 1.3 / 0.25 Mt | October | `tools/ch4_numbers.py` (the ejecta scales with the cable) |
| Table 4.8, 600 km row; 550 km | 15.4 / 12.6 / 9.3 / 5.6 / 2.1; (550 not in the log) | 14.0 / 11.5 / 8.5 / 5.1 / 1.9 | October | the same |
| Table 4.8, 250 km (5 % ejecta) | 82 Mt | 67 Mt | October | `tools/ch4_numbers.py` |
| 4.5, Figure 4.6, the ring mass against altitude | 5,978 Mt at 250, 6,124 at 700, 6,222 at 1,000 km; 5,590 Mt of cable at 700 | 4,965 at 250, 5,140 at 700, 5,220 at 1,000 (from `ring.design`); 4,606 Mt of cable | 4.5 October (6,124); Figure redrawn | `ring_altitude.design`, `tools/ch4_numbers.py` |
| Table 4.10, 250 km orbits: no recoil / slow / fast / radial / slow+radial / fast+radial | 250 x 3,630 / 1,480 / 6,710; 202 x 3,740; 144 x 1,640; 220 x 6,830 km | 250 x 4,023 / 1,762 / 7,289; 205 x 4,135; 159 x 1,904; 221 x 7,413 km | October | `tools/ch4_numbers.py` |
| Table 4.10, lifetimes, full section: nominal / slow / fast / radial / slow+radial / fast+radial | 470 / 136 / 920 / 143 / 5.6 / 438 years | 477 / 160 / 897 / 160 / 12.8 / 444 years | October | the same |
| Table 4.11, the 250 km cable: section, mass | 8.43 m section, 5,033 Mt structural (131,549 kg/m with hardware, 5,478 Mt) | 7.94 m, 4,019 Mt structural (107,229 kg/m, 4,466 Mt); soot at 5 % 67 Mt (0.3 %: 4.0 Mt) | October | `tools/ch4_numbers.py` (`250 km:` line); the soot rows scale with the cable |
| Table 4.12, 250 vs 700: perigee, lifetimes, soot, cable mass, ring mass | 144 to 250 / 587 to 700 km; 82 Mt / 0.3 to 7.8; 120,863 / 115,005 kg/m; 5,978 / 6,124 Mt | 159 to 250 / 602 to 700 km; 67 Mt / 0.25 to 7.1; 96,543 / 92,883 kg/m; 4,965 / 5,140 Mt | mixed (cable row new, the rest October) | `tools/ch4_numbers.py`, `ring_altitude.design(250)` |
| 4.6, "5,590 million tonnes of cable fragments" | 5,590 Mt | 4,606 Mt | check | `ring_altitude.py` |
| 4.6, the 3 degree tilt (429 m/s across the cable's speed) | 429 / 8,270 = 2.97 deg | 429 / 8,331 = 2.95 deg | unchanged at "about 3 degrees" | arithmetic |
| Figure 4.7, the design space | the floor at 600 km from the 587 km perigee | 602 km perigee (the floor stays at 600 km) | redrawn | `tools/ch4_figures.py` |

Note on Table 4.7: the book's lifetimes (39,400 / 3,700 / 370 years for the slow half with the kick) are not the numbers the current `tools/ch4_numbers.py` prints for that orbit under the fixed quadratic (65,950 / 6,265 / 627, cycle-averaged in its `life_w_*` rule) nor the quick-circular ones; they were made with an earlier version of the lifetime rule or by hand from Figure 4.3, so this row needs a decision on which rule the book uses before it is re-keyed. The new orbit (602 x 2,585 km) is higher, so every lifetime in Tables 4.7 and 4.10 rises; the ratios new to October in the tool's rule are 1.34 (full section), 1.49 (0.8 m) and 1.48 (8 cm) at 700 km.

## Chapter 5

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 5.7, "It is 8.22 m of carbon through its center" | 8.22 m | 7.39 m | October | `ring_altitude.py` |
| Figure 5.4, the sky fractions at 1 m below the cable: cable / Earth / open | 42.4 % / 28.3 % / 29.3 % | 41.6 % / 28.3 % / 30.1 % | caption new (7.39 m); figure redrawn with the new shares | `radiation/sweetspot_run.py` (`sweetspot_gcr.txt`) |
| Figure 5.4 and Table 5.12, the cosmic ray dose under the cable, thin walls, 1 m; with no cable | 23.2 mSv/yr (book: 23); 39.6 (book: 40) | 24.9 mSv/yr; 39.6 | October (23) | the same |
| Table 5.12, all rows: cable / Earth / open and dose at 0.5 / 1 / 1.5 / 2 / 3 / 5 m | 46/28/26, 20; 42/28/29, 23; 39/28/33, 25; 36/28/36, 27; 30/28/42, 29; 22/28/50, 30 | 46/28/26, 21.5; 42/28/30, 24.9; 38/28/34, 26.9; 34/28/38, 28.2; 28/28/43, 29.6; 20/28/51, 31.1 | October | the same |
| Table 5.11, trapped-proton dose under the cable, peak longitude 56 W: open 5 g/cm2 / under the cable thin / 1 m water / 3 m water / no cable 1 m / cable + 3 m floor | 65 / 26 / 0.39 / 0.039 / 0.72 / 0.047 Sv/yr | 64.5 / 25.9 / 0.43 / 0.076 / 0.72 / 0.085 Sv/yr | October | `radiation/sweetspot_trapped.py` (`sweetspot_trapped.txt`); the narrower cable lets more of the flux reach the thick-wall cases |
| Table 5.11, the other longitudes | 100 W 0.056 / 0.0009 / 0.0000; 90 W 2.14 / 0.015 / 0.0007; 80 W 4.18 / 0.051 / 0.0023; 70 W 8.96 / 0.11 / 0.0077; 60 W 20.1 / 0.25 / 0.026 | 0.056 / 0.0009 / 0.0000; 2.14 / 0.015 / 0.0006; 4.18 / 0.053 / 0.0033; 8.97 / 0.12 / 0.015; 20.2 / 0.28 / 0.052 | October | the same |
| 5.7, the share of the trapped flux arriving through the cable at the flux peak; from below the horizontal | 48 %; 50 % | 48 %; 50 % (unchanged to the printed digit) | same | the same |

## Chapter 7

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 7.1, Table 7.1, 7.2, 7.4: the relative velocity | 8,786 m/s | 8,847 m/s | new | `ring_altitude.py` |
| 7.1, the cable's speed westward | 8,270 m/s | 8,331 m/s | new | `ring_altitude.py` |
| 7.2, Equation 7.1, the drag to lift ratio at v | 0.0011 at 8,786 m/s | 0.0011 at 8,847 m/s (the ratio to two digits does not move) | new | `levitation_simulation.py` |
| 7.5 to 7.8, the bearing: the plates on the cable face | 4.33 m of the 8.22 m underside, 1.90 m on top | 4.33 m of the 7.39 m underside, 1.90 m on top (fits) | not quoted | `bearing_mass.py` |
| 7.5 to 7.8, the plates' edges from the centerline | 4.11 m cable edge | 3.70 m cable edge | not quoted | `final_bearing.py` |
| 7.x, the eddy heating gap: service temperature, melting | 51.3 mm, 33.3 mm | 51.4 mm, 33.4 mm | not quoted | `levitation_simulation.py` |
| Tables 7.3 to 7.5, the track ratings, tape, the 114 dynamics cases, smallest gap | unchanged (the casing's weight 95.0 kN/m sets them) | unchanged; smallest gap 50.2 mm | same | `final_bearing.py` |
| Figure 7.4, the bearing cross section | the 8.22 m cable face | the 7.39 m face (the plates and tracks the same) | redrawn | `tools/generate_fig_bearing_cross_section.py` |

## Chapter 8

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 8.x, Equation 8.1, the deployment demand | 88,946 x 8 MW = 712 GW | 88,946 x 12 MW = 1,067 GW | to check | `tools/ch08_numbers.py` (12 MW case) |
| 8.x, the deployment panel per site at 291 W/m^2 | 27,500 m^2, 55 m (27.5 m a side) | 41,282 m^2, 82.6 m (41.3 m a side) | to check | `tools/ch08_numbers.py`: "matched width 82.56 m (41.28 per side), 41,282 m^2" |
| 8.x, noon output per site, crossover | 22.6 MW, +/- 73.5 deg (8 MW) | 33.9 MW, +/- 73.5 deg | to check | the same |
| 8.x, peak flow, conductor | 171.5 GW, 86.5 kg/m (8 MW) | 257.3 GW, 25.7 kA, 129.7 kg/m, 5.77 Mt, 12.9 GW loss; 1.08 % of the casing | to check | the same |
| 8.x, after deployment: generation, operations, export | 711.6 - 44.5 - 8.6 = 658.5 GW | 1,067.4 - 44.5 - 12.9 = 1,010 GW | to check | the same |
| 8.x, per 10 km tap section | 160 MW, 1.6 kA | 240 MW, 2.40 kA | to check | the same |
| 8.x, export per anchor station | 741 MW (8 MW case) | 1,136 MW (1,010 GW / 889) | to check | arithmetic from the same |
| Table 8.1, the cryogenic plant: 88,946 sites at 1.53 MW = 136 GW | 1.53 MW, 136 GW | 1.37 MW, 122 GW (the standing load of the 150-turn site, Chapter 14) | October | `tools/ch14_thermal_series.py` (Chapter 14 below) |
| 8.x, the LIM supply frequency v_rel / 200 | 43.9 Hz | 44.2 Hz | to check | `tools/ch08_numbers.py` |

## Chapter 9 (the deployment; the design is 150 turns, 12 MW, two LIMs)

See `deployment_rising_radius_2026-11.md`: Section 13.2 for the 8 MW run of the new cable, Section 14 for the turns, Section 15 for the corrected slip rule and the power curve of the design. The rows the chapter quotes are listed at the end of this sheet under "Chapter 9 and 14, from the design run" once the reruns of Section 15 are in.

## Chapter 10

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| The static wave speed and the circuit time | 2,568 m/s, 4.8 hours (done, Paul's note) | 2,568 m/s (the same to the meter: sqrt(T/m) with T and m both scaled), 4.8 h | new | `tools/ch03_direction_numbers.py` |
| The prograde cable past the casing (catch-up speed) | 7,827 m/s | 7,904 m/s | new (Chapter 3) | `ring_altitude.py` |
| The subsonic band, prograde cable | 5,300 to 10,400 m/s | 5,375 to 10,433 m/s | check | `tools/ch03_direction_numbers.py` |
| The subsonic band, westward launches on the retrograde cable | 6,200 to 11,400 m/s | 6,279 to 11,415 m/s | check | the same |
| The sled passage: peak radial load, the hand check | 12.2 g, 12.6 g, 12.15 g, 3.7 % | the same (the passage loads do not depend on the cable's section) | same | `sled_passage.py` |
| The sled passage: the cable's design line, static L_char | 125,691 kg/m, 845.6 GN, 8.22 m, 8,270 m/s; L_char 15.0 m | 103,569 kg/m, 683.0 GN, 7.39 m, 8,330.8 m/s; L_char 13.5 m | the L_char table of `sled_passage.txt` changes with it (every u row) | `sled_passage.py` |

## Chapter 11

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| Table 2.1 and 11.x, the lateral force per anchor at 20 degrees | 5.331 GN (retrograde column) | 5.317 GN | not quoted as such (the chapter gives the force per meter) | `anchor_line_analysis.py --inclination=20` |
| The anchor lines sized from the ring's tension or mass | the line angles, lengths, stresses and tapers of `anchor_700_i*.txt` | 42 of 135 lines change in each log; the lines are sized from the casing's weight and the lateral force, which moves 0.3 % | to check line by line against `anchor_700_i20.txt` | `anchor_line_analysis.py` |
| 11.x, L_char in the anchor loads | 15.0 m static | 13.5 m | new (13.5) | `sled_passage.py` |

## Chapter 12

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 12.14, the cable per kilogram of casing | 10.5 kg (10.47) | 8.6 kg (8.63) | new | arithmetic |
| 12.14, the array's 1,200 kg/m claims | 12,600 kg/m of cable, 560 Mt | 10,400 kg/m (10,357), 460 Mt (461) | new | arithmetic |
| 12.14, a tenth of the deployment energy | 1.9 EJ of 18.8 | 1.9 EJ of 19.2 | new | arithmetic |
| 12.x, the deployment panel strip | 55 m or 110 m | 83 m (12 MW) or 110 m (16 MW) | October (55 m) | `tools/ch08_numbers.py` |

## Chapter 15

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| 15.x, "The cable moves past the casing interior at" | 8,786 m/s | 8,847 m/s | October | `ring_altitude.py` |
| The drag heating and erosion around the ring (30 MW, 41 MW, half a micrometer a year) | no script: `atmospheric_loading.py` is the tether's wind loading, not the ring's drag; nothing in the repository computes these | the quantities scale with the cable's wetted perimeter (4 x 7.39 = 29.6 m against 32.9 m, 0.90) and with the speed past the gas (8,847 against 8,786 m/s, 1.007; the heating as the cube, 1.021) | October | no script; arithmetic: 30 MW x 0.90 x 1.021 = 27.5 MW, 41 x 0.90 x 1.014 = 37.4 MW, the erosion per unit area the same (0.5 um/yr) if it is per unit area |

## Chapter 16

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| Table 16.1, deployment solar panels "55 m each side", 55 to 110 kg/m; HVDC lines for the 55 m panels, 43 kg/m | 55 m, 43 kg/m (8 MW) | 83 m (41 m each side), 83 to 165 kg/m at 0.5 to 1.0 kg/m^2; the HVDC conductor 129.7 kg/m at 12 MW | October | `tools/ch08_numbers.py` |
| Table 16.1, LIM stator tape, two LIMs per site | 6.5 kg/m (3,259 kg per site, 208 turns) | 4.7 kg/m (2,350 kg per site, 150 turns) | October | `lim_hysteresis.tape_mass_per_lim` (Section 14.2) |
| 16.x, "T is the design tension of the cable, 845.6 GN, ... A (67.65 m^2)" | 845.6 GN, 67.65 m^2 | 683.0 GN, 54.64 m^2 | October | `ring_altitude.py` |
| 16.x, "The result is about 15,000 kg/m, a casing and a quarter" (the added load at zero hoop stress) | 15,000 kg/m | 12,200 kg/m, about a casing | October | arithmetic |
| 16.x, "the number that sets the cable at 10.5 kg per kilogram of casing" | 10.5 | 8.6 | October | arithmetic |
| Table 16.2, the design case | 115,005 / 125,691 kg/m, 845.6 GN, 18.8 EJ, 6.124 Tt | 92,883 / 103,570 kg/m, 683.0 GN, 19.18 EJ, 5.140 Tt | October | `casing_6000.py --table` |
| Table 16.2, (a) 6,000 kg/m, hardware unchanged | 57,135 / 67,821 kg/m, 420.1 GN | 45,231 / 55,917 kg/m, 332.6 GN, 9.51 EJ, 2.754 Tt | October | the same |
| Table 16.2, (b) hardware halved | 57,170 / 67,306 kg/m, 420.4 GN | 46,442 / 51,785 kg/m, 341.5 GN, 9.59 EJ, 2.570 Tt | October | the same |
| Table 16.2, (c) plates sized, standing load halved | 57,280 / 65,763 kg/m, 421.2 GN | 45,724 / 54,207 kg/m, 336.2 GN, 9.54 EJ, 2.678 Tt (standing load kept: 45,353 / 55,489, 333.5 GN) | October | the same |
| 16.6, "halving the casing halves the cable" | (a) 49.7 %, (b) 49.7 %, (c) 49.8 % of the design case's structural cable | (a) 48.7 %, (b) 50.0 %, (c) 48.8 % kept, 49.2 % halved | holds | the same |
| Table 16.2, the deployment rows (days, LIM work for the 6,000 kg/m casing) | the October fixed-radius runs of `run_700_casing_6000.sh` | not rerun (they need rising-radius runs of the lighter casings with the 150-turn, 12 MW design) | October | `lim_simulation.py --rising-radius` with `configure(M_LOAD_M=6000, ...)`, not done |

## Chapter 17

| Where | October | New | Book now | Script |
|---|---|---|---|---|
| Table 17.1 | the orbital velocity, ring length and stress ratio by body; no cable mass or tension column | unchanged (the Earth row's 7,504 m/s, 44,473 km, 12.5 GPa do not depend on the quadratic) | same | no script (hand) |
| Table 17.2, the portable ring | no Earth numbers | unchanged | same | no script |
| The velocity-squared scaling sentences | none in the chapter quotes 846 GN or 125,691 kg/m | | same | |

## Figures that carry the cable

Redrawn from the new runs on 2026-11-07, originals in `Vol_III/images/figures_backup_2026-11-06b/`, the same layouts as the current files (the v14 scripts of `Vol_III/figure_scripts_2026-11-03/` for the figures proof v14 re-laid, copied with their paths fixed to `Vol_III/figure_scripts_2026-11-07/`; the repository tools for the rest):

| Figure | What changed | Drawn by |
|---|---|---|
| 1.1, the ring in cross section | the labels "8.22 m square" to "7.39 m square" and "27.5 m a side at 8 MW" to "41.25 m a side at 12 MW" (two labels) | the SVG's text edited (the figure is hand-drawn) |
| 3.2, the deployment tension | the design run: 150 turns, 12 MW, the new cable, 381 days, no anchored phase; the prograde dashed curve from the prograde design run | `Vol_III/figure_scripts_2026-11-07/gen_fig03_02_rising.py` (the v14 layout), pending the prograde run |
| 4.1, the released cable's orbit | 8,331 m/s, apogee 11,364 km radius | `figure_scripts_2026-11-07/ch4_figures_v14.py` |
| 4.2, 4.3, 4.4, 4.7 | the orbits, the full-section lifetime (B 5,712 kg/m^2), the ejecta soot, the design space (the 602 km perigee) | `tools/ch4_figures.py` (the October layout these four keep) |
| 4.5, the proton-free share | does not depend on the cable; not replaced | |
| 4.6, the ring mass against altitude | 104 t/m and 5,140 Mt at 700 km, 4,965 Mt at 250 km | `ch4_figures_v14.py` |
| 5.4, the sky under the cable | the 7.39 m cable, 41.6 / 28.3 / 30.1 %, 24.9 against 39.6 mSv/yr, 48 % | `figure_scripts_2026-11-07/f54_v14.py` on the new `sweetspot_*.txt` |
| 7.4, the bearing cross section | the 7.39 m cable face | `tools/generate_fig_bearing_cross_section.py` (PNG through the headless browser) |
| 14.1, 14.2, the thermal history | the design run's series (Chapter 14 below) | `tools/fig14_figures.py` |
| 16.1 to 16.4, 17.1 | the array geometry and the stress scaling: no cable numbers; no script in the repository; not redrawn | |
| 3.1, 7.1 to 7.3 | drawn figures without cable numbers | |

## Chapter 14, from the design run (150 turns, 12 MW, two LIMs; `tools/ch14_thermal_series.py --rising --set N_TURNS=150 --set MAX_SITE_POWER=12e6`)

The book's Table 14.2 carries the 208-turn, 8 MW, October-cable run (26.4 kW cold side and 4.28 m at seven weeks, 18.0 kW steady, 16.2 kW and 1.53 MW after the deployment). The design's history is different in shape: the 150-turn coil runs at its 183 A cap for 133 days and then at the 12 MW cap to the end, so its hysteresis and eddy loss peak at five months, not seven weeks, and stay high to the end. The compressor power here is the book's rule, the cold-side load over the COP of 0.0106; the simulation's own cryogenic power (`p_cryo`, its CRYO_EFF model) is given beside it.

| Table 14.2 row | October (book) | New: peak, at 5.0 months | New: steady state, months 6 to 10 | New: after deployment (standing load) |
|---|---|---|---|---|
| Hysteresis in the two LIMs | 10.2 / 1.8 / 0 kW | 16.3 kW | 12.0 kW | 0 |
| Leak through the cryostat insulation | 6.2 kW | 4.5 kW (the 150-turn coil's cryostat is 2,236 m^2 against 3,101) | 4.5 | 4.5 |
| Heat onto the levitation coils | 10.0 kW | 10.0 | 10.0 | 10.0 |
| Cold-side load | 26.4 / 18.0 / 16.2 kW | 30.8 kW | 26.5 kW | 14.5 kW |
| Compressor power (COP 0.0106) | 2.49 / 1.70 / 1.53 MW | 2.90 MW (simulation 2.57) | 2.50 MW | 1.37 MW (simulation 1.21) |
| Heat rejected at 400 K | 2.52 / 1.72 / 1.55 MW | 2.93 MW | 2.53 MW | 1.38 MW |
| Cryogenic radiator area | 1,925 / 1,314 / 1,183 m^2 | 2,247 m^2 | 1,935 m^2 | 1,058 m^2 |
| Radiator width over 450 m | 4.28 / 2.92 / 2.63 m | 4.99 m | 4.30 m | 2.35 m |

Every row changes. The peak is 17 % higher and comes at five months (the hysteresis peak, 16.3 kW per two LIMs, is at the same time); the steady state is 47 % higher because the coil is at full power to the end; the standing load after the deployment is 10 % lower because the narrower coil leaks less, so the ring's largest standing load (Chapter 8) is 1.37 MW per site, 122 GW. The cable's temperature (Figure 14.1, Table 14.3): 331 K on day 1 (not 382), 278 K at 1.5 months, a second rise to 313 K at four months while the coil is at its current cap, then 270 K from month 5 to the end (not 257 K), the wall 263 K; the end-of-run eddy loss is 20.6 kW per LIM at 468 N and 0.5 % slip against 11 kW before. The cold-side load starts at 14.5 kW, is 17.5 kW at one month, 25.0 at three, 30.8 at the peak (5.0 months), 26.3 at eight, 24.4 kW and 3.95 m at the end (12.7 months). Figures 14.1 and 14.2 are redrawn from this series.

## Chapter 9, from the design runs (Sections 14 and 15 of `deployment_rising_radius_2026-11.md`; the corrected slip rule)

| Where | October (the book's Table 9.4 and 9.8 text: the October cable, 208 turns, the October slip rule) | New | Book now | Script |
|---|---|---|---|---|
| Table 9.2, turns, peak current, LIMs | 208 turns, 183 A | 150 turns, 183 A, two LIMs | new | `tools/rising_radius_turns_sweep.py` |
| Table 9.3, the design column: days, limits, voltage at the end, thrust, energy | 647 days at 8 MW (the October cable) | 381 days (12.5 months) at 12 MW; current cap to day 133, power cap from day 130 to the end, 100 kV never; 99 kV and 125 A at the end; 1,013 N peak, 468 N at the end; 19.2 EJ | new | `tools/rising_radius_design_power.py` (L2_N150_P12MW) |
| Table 9.3, the 8 MW column | | 535 days (17.6 months); current cap to day 114, 8 MW from day 105; 79 kV, 100 A; 300 N at the end | new | the same (L2_N150_P8MW) |
| Table 9.4, 8 / 12 / 16 MW with 208 turns and two LIMs | 647 / 573 / 574 days, 100 kV from day 353 / 91 / 60, the last (anchored) phase 146 days, 5,590 Mt | with the design's 150 turns: 535 / 381 / 355 days (and 353 at 20 and 24 MW), 100 kV never / never / day 250, no anchored phase, 4,606 Mt; with 208 turns at 12 and 16 MW: 436 / 428 days, 100 kV from day 89 / 59 | October | the same |
| Table 9.4, four LIMs: 8 MW 208; 16 MW 208; 16 MW 150; 24 MW 208; 24 MW 150 | 643; 316; 310 days (and the six-LIM rows 210, 206) | 16 MW: 247 (208 turns), 255 (150); 24 MW: 220 (208), 188 (150) days; no anchored phase; 4,670 Mt; the 8 MW four-LIM and the six-LIM cases not rerun with the new cable | October | the same |
| Table 9.4, panel width and tape | 55 / 83 / 110 / 165 m; 3.3 / 6.5 t (208 turns), 4.7 / 7.1 t (150) | the same widths at the same powers; 2.35 t (150 turns, two LIMs), 3.26 (208, two), 4.70 (150, four), 6.52 (208, four) | mixed | `lim_hysteresis.tape_mass_per_lim`, `tools/ch08_numbers.py` |
| 9.8, "power alone saturates ... about 74 days ... 12 MW or 16" | the October cable with 208 turns | with 150 turns: 8 to 12 MW saves 153 days, 12 to 16 another 26, nothing beyond 20 MW; the knee of the 2 % rule is 16 MW (355 days); the 208-turn coil saturates at 12 MW (436, 428 days) | October | `tools/rising_radius_design_power.py` |
| 9.7, the turns sweep at 12 MW | the book's 100 / 150 / 208 / 300 at 8 MW | 650, 448, 381, 402, 401, 437, 541, 719 days at 120 to 280 turns; minimum 150 | new | `tools/rising_radius_turns_sweep.py` |
| Table 9.3 and 14.x, the cryogenic plant of the design | 2.49 MW peak at seven weeks, 1.53 MW after | 2.90 MW (the book's COP) at five months, 1.37 MW after; the simulation's own 2.57 and 1.21 MW | October | `tools/ch14_thermal_series.py --rising --set N_TURNS=150 --set MAX_SITE_POWER=12e6` |
