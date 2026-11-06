# The deployment with the rising radius

Report of the run asked for in `Vol_III/Claude_Code_instructions_deployment_radius_2026-11-06.md`, with the physics of `Vol_III/Vol_III_Strain_Radius_Plan_2026-11-06.md`. Date: 2026-11-06. Code: `rising_radius.py` and the `--rising-radius` flag of `lim_simulation.py`; the tools `tools/rising_radius_hanging.py`, `tools/fig_rising_radius.py` and `tools/fig03_02_deployment_tension.py --rising`. Every number below is from the run's log and files, `runs/rising_radius_lim_simulation.txt`, `rising_radius_daily.csv`, `rising_radius_deployment_tension.csv`, `rising_radius_summary.json` and `rising_radius_hanging.csv`, and from `ring_altitude.py`. Nothing in the book, in the design point or in the default fixed-radius run was changed; the default run reproduces its log (Section 8).

## 1. The problem

The cable's operating strain is sigma / E = 12.5 GPa / 500 GPa = 2.5 %. A closed loop is longer by its strain, C = C0 (1 + eps), and nothing but the cable's own speed can stretch a 43,000 km loop, so the ring is assembled unstressed at r0 = R / (1 + 0.025) and rises to R = 7,078,137 m (700 km) as the tension builds. The fixed-radius run of Chapters 3 and 9 (`lim_simulation.py` as it stands) keeps the ring at 700 km from the first day and conserves linear momentum per meter. With the radius rising, angular momentum is the conserved quantity, the circular velocity at the start is that of the build orbit, the LIMs must also lift the ring, and the design point carries more angular momentum than the ring has at the build orbit, which only the Earth can supply through the anchor lines. The run answers how much the deployment's time, energy, power and the final tensioning change.

## 2. The method

The rising-radius mode changes only the ring's kinematics. The LIM model is exactly the fixed-radius one: the same thrust model, slip control, current and voltage limits, plate and cable temperatures, cryogenic load and hysteresis, with the same time steps (1 s on the first day, then 10 s) and the same number of sites, 88,946. The state of the ring is the cable's angular momentum L_c and the casing's L_l about the Earth's axis (totals for the whole ring), the radius r, and the unstressed circumference C0 = 2 pi r0, which is fixed. The total masses are fixed: M_r = m_rotor x 2 pi R = 5,589.9 Mt for the cable with its hardware (125,691 kg/m at the design point) and M_l = 533.7 Mt for the casing (12,000 kg/m). The per-meter masses at any radius are M / (2 pi r).

Each step:

1. The total LIM thrust F (the thrust of one LIM from the site model, times two LIMs per site, times 88,946 sites) moves angular momentum from the casing to the cable: dL_c = +F r dt, dL_l = -F r dt. Westward is positive, the cable is retrograde, the casing slows from westward orbital speed through zero to eastward ground speed, as in the fixed-radius run. No external torque in phase 1.
2. The speeds from the momenta: v_c = L_c / (M_r r), v_l = L_l / (M_l r).
3. Radial equilibrium at the current r, the same force balance as the fixed-radius `deployment_tension()` with the casing's own centrifugal term and g evaluated at r:

   T = m_r (v_c^2 - GM / r) - m_l (GM / r - v_l^2), with m_r = M_r / (2 pi r) and m_l = M_l / (2 pi r).

   With the momenta given this is T(r) = a / r^3 - b / r^2, a = (L_c^2 / M_r + L_l^2 / M_l) / (2 pi), b = GM (M_r + M_l) / (2 pi).
4. Strain and radius: eps = T / (E A) with E = 500 GPa and A = 67.65 m^2, the structural CNT area at the design point (the area the stress of the book is figured on); r = r0 (1 + eps). Steps 2 to 4 are a fixed point in r. The code solves r - r0 (1 + T(r) / E A) = 0 by Newton's method from the previous step's radius, to 1e-6 m; it never needed more than 2 iterations. The rise is treated as quasi-static (the radial period is of the order of the orbit, the deployment is hundreds of days; no radial dynamics).
5. Phase 2 begins when the casing reaches the ground-synchronous speed, v_l = -omega_E r (omega_E = 7.2921159e-5 rad/s). From then on the anchor lines hold the casing at that speed: L_l = -M_l omega_E r^2 each step, the LIMs keep pushing the cable, and the difference between what the LIM reaction takes from the casing and what its momentum requires, dL_E = L_l(new) - (L_l(old) - F r dt), is what the Earth supplies. The external torque is dL_E / dt and the anchor reaction is that over r, as a total around the ring and over round(L_RING / 50 km) = 889 stations. Phase 2 ends when the cable reaches the design speed (8,270.0 m/s, `ring_altitude.V_CABLE`) or the tension reaches the design tension (845.6 GN, `ring_altitude.CABLE_TENSION`), whichever comes first.
6. Energy, both ways. The LIM work is the sum of thrust times relative speed times dt, as in the fixed-radius run. The mechanical energy is E = 1/2 M_r v_c^2 + 1/2 M_l v_l^2 - GM (M_r + M_l) / r + 1/2 E A eps^2 C0, the last term the elastic energy of a linear loop stretched by eps C0, which equals 1/2 T eps C0. In phase 2 the Earth's work through the anchors is the torque times the casing's angular rate, -omega_E (negative: the casing turns eastward with the ground while the Earth pulls it westward, so energy passes from the ring to the Earth). Because the radius is always at equilibrium, dE/dr = 0 at fixed momenta, and the balance "LIM work = delta E - Earth's work" holds exactly in continuous time; in the run it holds to the discretization of the 10 s step, as the fixed-radius run's two energies do.

Three choices were made differently from, or beyond, the brief:

- The elastic energy uses the unstressed circumference C0 (1/2 T eps C0), not the stretched one. For a linear loop the stretch is eps C0 and the stored energy is the integral of T over it; with C instead the balance would miss by 2.5 % of the elastic term, about 0.01 EJ. (The brief wrote 1/2 T eps C.)
- The spreader arm in the hanging calculation is the arm as Section 12.11 hangs it during the deployment, vertical with its lower end held at 150 km: 550 km below a 700 km casing (the same 550 km `tools/anchor_lines_deployment.py` pays out). The arm's full length, 100 km to the array and 893 km below the kink (Section 12.3), is longer than the ring is high and is never hung straight down, so it was not used.
- The LIM site spacing is 487.8 m at the build and 500 m at the design point, but the LIM model keeps `LIM_SPACING` = 500 m throughout, as the brief asks (the LIM physics untouched); the spacing enters the model only through the plate's local heating duty, a 2.5 % effect on a term that is itself small.

The brief's estimate of the final tensioning, 650 MN over 10 days or 220 MN over 30, assumed the Earth's share could be pushed in quickly. The run shows it cannot: the LIMs are what push the cable in phase 2, and their thrust at the end of the deployment is about 45 MN in total (250 N per LIM at the 100 kV voltage limit), so phase 2 takes as long as that thrust needs. See Section 4.

## 3. The build orbit

| Quantity | Value |
|---|---|
| Design radius and altitude | 7,078,137 m, 700.0 km |
| Design strain, sigma / E | 12.5 GPa / 500 GPa = 2.500 % |
| Build radius r0 = R / (1 + eps) | 6,905,499.5 m |
| Build altitude | 527.36 km |
| Rise to the design point | 172.64 km |
| Circular velocity at the build orbit | 7,597.5 m/s (7,504.3 m/s at 700 km) |
| Unstressed circumference C0 | 43,388.5 km (44,473.0 km at the design point) |
| Stretch at the design point, eps C0 | 1,084.7 km |
| LIM site spacing at the build, C0 / 88,946 | 487.8 m (500.0 m at the design point) |
| Anchor stations, round(L_RING / 50 km) | 889 |
| Angular momentum of the ring at the build | 3.2127e23 kg m^2/s |
| Angular momentum at the design point | 3.2526e23 kg m^2/s, 1.242 % more |

## 4. The run

The deployment takes 647.1 days (1.77 years, 21.3 months) against 499.7 days at the fixed radius: 501.1 days of phase 1 until the casing is at the ground-synchronous speed, then 146.1 days of phase 2 against the Earth. The design point is reached exactly: 8,270.0 m/s, 516.1 m/s eastward, 8,786.2 m/s relative, 845.63 GN, 12.500 GPa, 2.5000 % strain, 700.00 km, casing net weight 95.0 kN/m, cable excess force 214.5 kN/m. Phase 2 ended on the design tension, with the cable 0.002 m/s short of the design speed; the two criteria are the same state and the run reached them in the same step to round-off, so there is no meaningful "which first". The Newton solve of the radius took at most 2 iterations in any step.

The ring through the run (the daily table is `rising_radius_daily.csv`; the marks interpolated on the sampled series):

| Casing slowdown | Casing speed | Cable speed | Tension | Stress | Altitude | Casing weight | Cable excess | Day |
|---|---|---|---|---|---|---|---|---|
| 0 | 7,598 W | 7,598 | 0.0 GN | 0.00 GPa | 527.4 km | 0.0 kN/m | 0.0 kN/m | 0.0 |
| 1,000 | 6,597 W | 7,690 | 10.8 GN | 0.16 GPa | 529.6 km | 25.2 kN/m | 26.8 kN/m | 25.1 |
| 3,000 | 4,597 W | 7,860 | 96.2 GN | 1.42 GPa | 547.0 km | 64.5 kN/m | 78.4 kN/m | 87.3 |
| 5,000 | 2,597 W | 8,010 | 264.2 GN | 3.91 GPa | 581.3 km | 88.6 kN/m | 126.6 kN/m | 202.3 |
| 7,598, casing at rest | 0 | 8,178 | 600.6 GN | 8.88 GPa | 650.0 km | 97.5 kN/m | 183.0 kN/m | 435.6 |
| 8,111, ground speed (phase 2 begins) | 514 E | 8,207 | 682.3 GN | 10.09 GPa | 666.7 km | 96.4 kN/m | 193.2 kN/m | 501.1 |
| 8,114, design point (end) | 516 E | 8,270 | 845.6 GN | 12.50 GPa | 700.0 km | 95.0 kN/m | 214.5 kN/m | 647.1 |

Daily, every 50th day (the file has every day):

| Day | Altitude km | Tension GN | Stress GPa | Strain % | v_cable m/s | v_casing m/s (W positive) | v_rel m/s | Phase | LIM thrust total MN | Anchor reaction total MN |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 527.36 | 0.00 | 0.000 | 0.0000 | 7,597.5 | 7,597.5 | 0.0 | 1 | 0.0 | 0 |
| 50 | 536.49 | 44.70 | 0.661 | 0.1322 | 7,781.2 | 5,558.0 | 2,223.2 | 1 | 188.8 | 0 |
| 100 | 550.70 | 114.33 | 1.690 | 0.3380 | 7,881.9 | 4,324.9 | 3,557.0 | 1 | 127.1 | 0 |
| 150 | 565.56 | 187.13 | 2.766 | 0.5532 | 7,952.5 | 3,399.7 | 4,552.8 | 1 | 102.0 | 0 |
| 200 | 580.60 | 260.79 | 3.855 | 0.7710 | 8,008.0 | 2,630.8 | 5,377.2 | 1 | 87.4 | 0 |
| 250 | 595.68 | 334.63 | 4.947 | 0.9893 | 8,054.2 | 1,960.2 | 6,093.9 | 1 | 77.6 | 0 |
| 300 | 610.74 | 408.40 | 6.037 | 1.2074 | 8,093.9 | 1,358.5 | 6,735.3 | 1 | 70.5 | 0 |
| 350 | 625.75 | 481.95 | 7.124 | 1.4248 | 8,128.7 | 808.6 | 7,320.2 | 1 | 65.0 | 0 |
| 400 | 640.26 | 553.00 | 8.174 | 1.6349 | 8,159.0 | 314.2 | 7,844.8 | 1 | 56.9 | 0 |
| 450 | 653.75 | 619.08 | 9.151 | 1.8303 | 8,184.6 | -119.1 | 8,303.8 | 1 | 50.5 | 0 |
| 500 | 666.40 | 681.04 | 10.067 | 2.0134 | 8,206.9 | -506.0 | 8,712.9 | 1 | 45.5 | 0 |
| 550 | 677.91 | 737.43 | 10.901 | 2.1801 | 8,228.4 | -514.5 | 8,743.0 | 2 | 45.1 | 44.9 |
| 600 | 689.32 | 793.33 | 11.727 | 2.3454 | 8,249.9 | -515.4 | 8,765.3 | 2 | 44.9 | 44.7 |
| 647.1 | 700.00 | 845.63 | 12.500 | 2.5000 | 8,270.0 | -516.1 | 8,786.2 | 2 | 44.6 (log) | 44.4 |

The ring rises 0.30 km per day through the middle of the run and 0.23 km per day in phase 2. The casing passes through rest on day 435.6 at 650.0 km with the cable at 601 GN (fixed radius: day 433.2, 740 GN), and reaches the ground-synchronous speed, 513.7 m/s eastward at that radius, on day 501.1 at 666.66 km with the cable at 8,207.3 m/s and 682.3 GN, 10.09 GPa: the state the plan computed by hand. The day marks the book quotes: after 1 day 0.01 GN and 0.63 kN/m (slowdown 23.1 m/s); after 10 days 0.89 GN and 7.60 kN/m; after 30 days 17.3 GN and 31.3 kN/m (slowdown 1,265 m/s); after 100 days 114.3 GN and 68.7 kN/m. The casing's net weight reaches half its final 95 kN/m on day 50 (fixed radius: day 55) and peaks at 97.8 kN/m on day 394 before easing to 95.0 as the ring rises. The peak rise of the tension is 1.48 GN per day on day 222 (fixed: 1.85 GN per day on day 350).

### Phase 2

| Quantity | Value |
|---|---|
| Begins (casing at ground-synchronous speed) | day 501.1, 666.66 km, cable 8,207.3 m/s, 682.3 GN, 10.09 GPa |
| Ends (design tension and speed) | day 647.1, 700.00 km, 8,270.0 m/s, 845.6 GN |
| Duration | 146.1 days |
| Rise, cable gain, tension gain in phase 2 | 33.34 km, 62.7 m/s, 163.3 GN |
| Angular momentum from the Earth | 3.991e21 kg m^2/s |
| Its share of the cable's final angular momentum (3.272e23) | 1.220 % |
| Its share of the ring's at the design point | 1.242 % |
| Anchor reaction, total around the ring: peak | 45.16 MN on day 501.1 |
| Anchor reaction, total: mean over phase 2 | 44.79 MN |
| Per station (889 stations): peak, mean | 50.8 kN, 50.4 kN |
| LIM thrust, total, at the end of the run | 44.6 MN (250.9 N per LIM, 177,892 LIMs) |
| The Earth's work through the anchors | -0.291 EJ (from the ring to the Earth) |

The anchor reaction is the LIM thrust, almost exactly: the LIMs push the cable westward and the casing eastward, the casing cannot take the eastward push because the anchors hold it at ground speed, so the anchors carry the whole LIM reaction to the ground, 45 MN around the ring, 50 kN per station, for 146 days. The small difference between the two (0.2 MN) is the casing's own angular momentum growing with r^2 as the ring rises. Phase 2 is long because the cable's angular momentum has to grow by 4.0e21 kg m^2/s at a torque of 45 MN x 7,060 km, and only part of that momentum appears as speed: of the 101.5 m/s of L_c / (M_r r) added in phase 2, 62.7 m/s is the cable's speed gain and 38.8 m/s is taken up by the radius rising 33 km (v_c dr / r). The brief's 10 to 30 days with 650 to 220 MN assumed a much larger push than the LIMs can give at the end of the deployment, where the 100 kV coil voltage limits them to 250 N each.

### The energy budget

| Term | EJ |
|---|---|
| LIM work, thrust x relative speed, summed in the loop | 23.864092 |
| Cable kinetic energy change (161.33 to 191.16) | +29.824756 |
| Casing kinetic energy change (15.40 to 0.07) | -15.331467 |
| Potential energy change, -GM M / r, the ring lifted 172.64 km | +8.621142 |
| Elastic energy of the cable, 1/2 E A eps^2 C0 | +0.458631 |
| Total mechanical energy change | 23.573062 |
| The Earth's work through the anchors in phase 2 | -0.291036 |
| LIM work required, delta E minus the Earth's work | 23.864098 |
| Difference between the two ways | -0.000006 (6e-6 EJ, as the fixed-radius run's two energies differ) |
| Fixed-radius deployment energy (`ring_altitude.E_DEPLOY`, and the fixed run's 18.80385) | 18.803860 |

The LIM work is 23.86 EJ, 26.9 % more than the 18.8 EJ of the fixed radius: 23.57 EJ into the ring's mechanical energy (14.49 EJ net kinetic, 8.62 EJ of lift, 0.46 EJ elastic) and 0.29 EJ through the anchors into the Earth. The plan's first estimate of 23.6 EJ for the mechanical change was right; the energy through the anchors is on top of it. The energy per site is 268.3 TJ against 211.4 TJ.

## 5. The LIM operating history, fixed radius against rising radius

The LIM model sees only v_rel, and v_rel through the first year is within a few percent of the fixed-radius run's, so every peak and limit is the same to the precision the book quotes. What moves is the timing of the voltage limit and the length of the run.

| Quantity | Fixed radius (`lim_simulation.txt`, 2026-10-27) | Rising radius (`rising_radius_lim_simulation.txt`) |
|---|---|---|
| Deployment time | 499.7 days, 1.37 years, 16.4 months | 647.1 days, 1.77 years, 21.3 months (phase 1: 501.1 days) |
| Total LIM energy | 18.804 EJ | 23.864 EJ |
| Energy per site | 211.4 TJ | 268.3 TJ |
| Energy balance of the two sums | 6.0e-6 EJ | 6.0e-6 EJ |
| Peak coil current | 183.1 A, day 1 | 183.1 A, day 1 |
| Peak slip speed | 191.2 m/s, day 1 | 191.2 m/s, day 1 |
| Peak thrust per LIM | 1,937.5 N at 27.6 days | 1,937.6 N at 27.6 days |
| Site power reaches 8 MW | 27.6 days | 27.6 days |
| Coil voltage reaches 100 kV | 11.66 months | 11.82 months |
| Peak thrust power per LIM | 2.68 MW at 11.66 months | 2.68 MW at 11.82 months |
| Peak LIM power | 2.83 MW at 11.66 months | 2.83 MW at 11.82 months |
| Site power at the end | 6.81 MW | 6.81 MW (6.86 MW through phase 2) |
| Thrust at the end, per LIM | 250.9 N | 250.9 N (254.7 N at the start of phase 2) |
| Supply frequency at the end | 44.15 Hz | 44.15 Hz |
| Slip at the end | 0.5 % | 0.5 % |
| Peak eddy loss per LIM | 156.9 kW, day 1 | 156.9 kW, day 1 |
| Peak hysteresis per LIM | 5.09 kW at 1.75 months | 5.09 kW at 1.75 months |
| Peak cryogenic power per site | 2.20 MW at 1.75 months | 2.20 MW at 1.75 months |
| Peak cold-side load, radiator width | 26.4 kW, 4.28 m at 1.75 months | 26.4 kW, 4.28 m at 1.75 months |
| Peak cable, wall temperature | 382.2 K, 371.4 K, day 1 | 382.2 K, 371.4 K, day 1 |
| Cable temperature after the first months | about 259 K, 257 K at month 15 | about 259 K, 255 K through phase 2 |
| Peak warm-loop heat load | 728.2 W/m, day 1 | 728.2 W/m, day 1 |
| Cryogenic power at the end | 1.465 MW | 1.465 MW |
| Month 1: current, voltage, thrust, site power | 176.4 A, 42.6 kV, 1,768 N, 8.0 MW | 176.4 A, 42.6 kV, 1,769 N, 8.0 MW |
| Month 12: current, voltage, thrust, site power | 77.8 A, 100 kV, 352 N, 7.89 MW | 78.4 A, 100 kV, 357 N, 7.94 MW |
| Months at 100 kV | 4.8 (11.7 to 16.4) | 9.4 (11.8 to 21.3) |
| Months of phase 2, 6.86 to 6.81 MW per site, thrust 255 to 251 N | none | 4.8 (16.5 to 21.3) |

The 16 MW per site case (427.3 days at the fixed radius) was not rerun; the brief did not ask for it. By the same mechanism its phase 2 would be governed by the thrust at its 100 kV limit.

## 6. Every number the book quotes from this simulation, fixed against rising

Listed from Chapters 3, 8, 9 and 14 of the manuscript (`Vol_III/book`; L numbers are lines of the chapter files). "Same" means the rising run gives the same value to the precision the book quotes. The design point itself (8,270, 516, 8,786 m/s, 846 GN, 12.5 GPa, 44 Hz, 115,005 kg/m, 67.65 m^2) is the same by construction: the run ends at it.

| Chapter, place | The book (fixed radius) | Rising radius |
|---|---|---|
| 3 intro, 3.1, 3.2, Table 3.1, 3 summary; 9 intro, 9.8, 9.9, 9.10: 8,270 m/s, 516 m/s, 8,786 m/s, 846 GN, 12.5 GPa, 44 Hz | design point | same |
| 3.1 (L13), 3.2 (L52, L64, L65, L78, L126, L137, L150), Table 3.1, 9 intro (L5), 9.10 (L259): orbital velocity 7,504 (7,504.3) m/s at the start | 7,504 m/s at 700 km | 7,598 m/s at the build orbit (527 km); 7,504 m/s is the circular velocity at the design radius, which the ring reaches only at the end |
| 3.1 (L13), 3.2 (L60, L124, L156), Table 3.1 (L164), 9 intro: casing slowdown 8,020 (8,020.4) m/s | 8,020 m/s | 8,114 m/s (7,598 westward to 516 eastward) |
| 3.1 (L13), 3.2 (L137): cable speed gain 766 (765.7) m/s | 766 m/s | 672 m/s (7,598 to 8,270); the cable's angular momentum gains more than its speed shows because the radius rises |
| 3.1 (L17), 3.2 (L146), Table 3.1 (L178); 8 intro (L3); 9 Table 9.3 (L167), 9.6 (L176), 9.10 (L261): deployment energy 18.8 EJ | 18.8 EJ | 23.9 EJ of LIM work (23.6 EJ into the ring, 0.3 EJ through the anchors into the Earth) |
| 3.2 (L146): "integrates the LIM power over the 500 days", "better than one part per million" | 500 days, 6e-6 of 18.8 EJ | 647 days; the two energies differ by 6e-6 EJ, 2.5e-7 of the total |
| 3.2 (L39): 2.5 % strain, 1,100 km, "about 180 km" | 1,100 km, 180 km | 2.5 %, 1,085 km (eps C0), 172.6 km of rise from a build orbit at 527.4 km |
| 3.2 (L62): "21 GN after 30 days" | 20.94 GN | 17.3 GN after 30 days |
| 3.2 (L62), Figure 3.2 caption (L65): "half of its final 95 kN/m by day 55" | day 55 | day 50 (and a peak of 97.8 kN/m on day 394) |
| 3.2 (L62): "passes through rest on day 433 with the cable at 740 GN" | day 433, 740 GN, 10.94 GPa | day 436, 601 GN, 8.88 GPa, at 650 km |
| 3.2 (L62), Figure 3.2 (L64, L65): "845.6 GN and 12.5 GPa only at the end, on day 500" | day 500 | day 647; the casing is at ground speed from day 501 at 667 km, 682 GN, and the last 163 GN come in 146 days against the Earth |
| Figure 3.2 (L64): 95 kN/m, 214.5 kN/m, 516 m/s eastward, "over the 500 days" | | 95.0 kN/m, 214.5 kN/m, 516 m/s: same; over 647 days |
| 3.2 Table 3.1 (L176), L182: cable speed past the casing as a multiple of the wave speed, 3.4 | 3.4 | same (the design point) |
| 3.2 (L37, L154, L156), 3.1 (L11, L19): the stress figures 22.4 / 25.2 GPa, the wave-speed bands, 12.1 / 13.3 x 10^6 m^2/s^2, the prograde 6,988, 839, 8,343, 7,827 m/s, 865 GN | from the quadratic at 700 km | not changed by this run: they are design-point quantities (the quadratic at the design radius); the prograde deployment was not rerun with the rising radius |
| 8 intro (L3), 8.1 (L24, L26, L31, L33, L58), 8.2 (L68, L84, L94), 8.3 (L110), Table 8.2 (L136, L145): 88,946 sites at 8 MW, 712 GW | 712 GW | same peak; the demand runs 21 months instead of 16: 8 MW per site from day 28 to month 11.8, then 7.9 down to 6.8 MW to month 21.3 |
| 8.2 (L100): 44 Hz at the end, 8,786 m/s | | same |
| 8.3 (L106, L111), Table 8.2 (L137, L139): 1.5 to 1.53 MW of cryogenic power after the deployment, 136 GW; 1.7 MW steady, 2.5 MW peak | from Table 14.2 | the run's cryogenic power is 2.20 MW at the peak and 1.465 MW at the end in both runs; Table 14.2's 2.49 and 1.53 MW are figured from the cold-side load with the constants of Chapter 14 and do not come from the log's cryo line (an existing difference, not one this run introduces) |
| 9.2 (L28), Table 9.2 (L133, L134), Eq 9.1 (L31), 9.6 (L178): slip 5 to 200 m/s, 100 % at the start, 2 % or less at the end, under 1 Hz | 191.2 m/s, 0.5 %, 0.96 Hz | same |
| Table 9.2 (L126, L128, L131, L132, L135, L140, L141): 183 A, 100 kV, 44 Hz, about 1,900 N, 8 MW | | same (1,937.6 N) |
| Table 9.2 (L142), Table 9.3 (L160), 9.5, 9.8 (L212), 9.10 (L257), 9 summary (L265): 500 days (16.4 months), "sixteen months", "14 to 17 months" | 499.7 days | 647.1 days (21.3 months); 427 days at 16 MW not rerun |
| Table 9.3 (L161, L163), 9.6 (L172), 9 summary (L265): the 100 kV limit from 11.7 months, "the last five months" | 11.66 months | 11.82 months, the last 9.4 months |
| Table 9.3 (L162): 183 A on the first day | | same |
| Table 9.3 (L164), 9.6 (L172): 8 MW from day 28, easing to 6.9 MW in the last months | 27.6 d, 6.81 MW | same (6.86 MW through phase 2, 6.81 at the end) |
| Table 9.3 (L165), 9.6 (L172): 1,940 N at four weeks | 1,937.5 N | 1,937.6 N |
| Table 9.3 (L166): 2.7 MW at 11.7 months, thrust about 350 N | 2.68 MW, 352 N | 2.68 MW at 11.8 months, 357 N |
| Table 9.3 (L168, L169, L170), 9.6 (L178), 9.7 (L208): 157 kW on the first day, 5.1 kW at seven weeks, 382 K settling near 257 K | | same (the cable settles near 259 K, 255 K through phase 2) |
| 9.6 (L172): thrust "falls to about 260 N by the end", "177,892 LIMs for well over a year" | 250.9 N, 1.37 years | 250.9 N; 1.77 years |
| 9.7 (L204, L206): 208 turns reaches 100 kV after 2.5 months at 16 MW; 150 turns 26 and 77 days slower | variant runs | not rerun |
| 14 intro (L3): 88,946 sites at 500 m spacing | 500 m | 487.8 m at the build, 500 m at the design point (the sites are fixed in number; the casing's joints take the 2.5 %) |
| 14.1 (L18, L20), Table 14.1 (L26 to L32), 14.2 (L38 to L43), Figures 14.1 and 14.2, 14.3 (L49, L53, L54), Table 14.2, 14.4 (L81, L83), Table 14.3, 14.6: 314 kW, 10 kW, 382 K, 371 K, 257 K, 250 K, 26.4 kW at seven weeks, 4.28 m, 728 W/m, 18 kW and 2.9 m in months 12 to 15 | | same: every thermal peak is in the first seven weeks and the LIM model is unchanged; the steady-state window (months 12 to 15) has the same values; "over 17 months" on the figures becomes 21 months, and the cold loop runs 4.8 months longer at 16 kW and 2.85 m |
| 14.3 (L102): "trim its speed to hold the 12.5 GPa of tension" | | same |

Not quoted in the book but in the fixed run's log and changed: the tension marks at 1,000 / 3,000 / 5,000 m/s of casing slowdown (13.1 / 118.3 / 328.6 GN become 10.8 / 96.2 / 264.2 GN, at 530 / 547 / 581 km); the tension after 1, 10 and 100 days (0.01 / 1.07 / 139.1 GN become 0.01 / 0.89 / 114.3 GN); the peak tension rise (1.85 GN per day on day 350 becomes 1.48 on day 222).

## 7. What hangs: the array and the arms

From the daily altitude profile (`tools/rising_radius_hanging.py`, output `rising_radius_hanging.txt` and `rising_radius_hanging.csv`). A point d below the casing is at h(t) - d:

| Day | Ring | 20 km below | 50 km below | 100 km below | Spreader arm, 550 km below | Longest hang above 150 km |
|---|---|---|---|---|---|---|
| 0 | 527.4 | 507.4 | 477.4 | 427.4 | -22.6 (below the ground) | 377.4 |
| 30 | 530.9 | 510.9 | 480.9 | 430.9 | -19.1 | 380.9 |
| 100 | 550.7 | 530.7 | 500.7 | 450.7 | 0.7 | 400.7 |
| 200 | 580.6 | 560.6 | 530.6 | 480.6 | 30.6 | 430.6 |
| 300 | 610.7 | 590.7 | 560.7 | 510.7 | 60.7 | 460.7 |
| 400 | 640.3 | 620.3 | 590.3 | 540.3 | 90.3 | 490.3 |
| 500 | 666.4 | 646.4 | 616.4 | 566.4 | 116.4 | 516.4 |
| 600 | 689.3 | 669.3 | 639.3 | 589.3 | 139.3 | 539.3 |
| 647.1 | 700.0 | 680.0 | 650.0 | 600.0 | 150.0 | 550.0 |

- 20, 50 and 100 km below the casing: never below 150 km. The array on its 100 km support lines would start at 427 km and rise to 600 km with the ring; it is above the 150 km limit on every day, so the limit does not set its schedule. What does is where the array should be while the ring is low (427 to 600 km through 647 days: the air, the debris, the shadow), which Paul decides.
- The spreader arm as Section 12.11 hangs it, 550 km below the casing: below 150 km on every day of the run, and below the ground for the first 97 days. Hung from the first day it would have to be paid out with the ring: the longest hang that stays above 150 km is 377 km on day 0, 400 km on day 98, 450 km on day 265, 500 km on day 436, 550 km only at the end. A 600 km hang never clears 150 km.
- The anchor lines are attached when the casing comes to rest over the ground on day 501 at 667 km, and they pay out 33 km in the 146 days of phase 2 while carrying 50 kN per station.

## 8. Sanity check: the default run is unchanged

The fixed-radius run with the flag off, made with the edited code, reproduces the reference log of 2026-10-27 (`runs/lim_simulation.txt`, the log the book was proofed against) line for line apart from the timestamp, and `lim_deployment_tension.csv` byte for byte (4,076 rows). The unedited code was also run first, as a control, and gave the same two files. `python -m pytest tests` passes (26 fast tests; the slow `lim_simulation` log test passes on the new default log).

The run with the flag on and `--quick` (10 times the time steps) gives 647.1 days and 23.864 EJ as well, with an energy balance of 6e-5 EJ; the full run's is 6e-6 EJ.

## 9. Files and figures

Under `runs/` (the output folder, not tracked by git; a copy is in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/`):

- `rising_radius_lim_simulation.txt`: the log (the LIM history as in the fixed run, then the rising-radius report).
- `rising_radius_daily.csv`: day by day, altitude, radius, tension, stress, strain, speeds, casing weight, cable excess, phase, LIM thrust, anchor reaction total and per station, the three angular momenta, the energies. (The closing row's thrust column reads 0: the run was started before the fix that carries the last thrust into it; the final thrust is in the log's "final" row, 250.9 N per LIM, 44.6 MN in total.)
- `rising_radius_deployment_tension.csv`: the Figure 3.2 series (the same seven columns as `lim_deployment_tension.csv`, plus altitude, radius, strain, phase and anchor force), 4,053 rows.
- `rising_radius_summary.json`: the headline numbers.
- `rising_radius_hanging.txt`, `rising_radius_hanging.csv`: Section 7.
- `figures_rising_radius/rising_radius_altitude.svg|png`: the ring's altitude, the point 100 km below, the longest hang above 150 km.
- `figures_rising_radius/rising_radius_tension.svg|png`: tension and stress against time, rising against fixed.
- `figures_rising_radius/rising_radius_speeds.svg|png`: cable and casing speeds, both runs.
- `figures_rising_radius/rising_radius_anchor_phase2.svg|png`: the anchor reaction and the LIM thrust through phase 2.
- `figures_rising_radius/rising_radius_energy.svg|png`: the energy budget accumulating.
- `figures_rising_radius/fig03_02_deployment_tension_700_rising.svg|png`: Figure 3.2 redrawn from the rising series (`tools/fig03_02_deployment_tension.py --rising`). Its top and middle panels are against the casing's speed, which is constant through phase 2, so the last 163 GN appear as a vertical step at 8,114 m/s; a Figure 3.2 for the rising radius wants the cable's speed or time on that axis, which is a drawing decision for the book.

Reproduce: `python lim_simulation.py --rising-radius > runs/rising_radius_lim_simulation.txt` (about 3 hours when run beside two other deployment runs), then `python tools/rising_radius_hanging.py`, `python tools/fig_rising_radius.py`, `python tools/fig03_02_deployment_tension.py --rising --out=runs/figures_rising_radius`.

## 10. Round 2 (2026-11-06): the reruns the book quotes and the deployment study

Asked for in `Vol_III/Claude_Code_instructions_deployment_radius_round2_2026-11-06.md`. The code gained two things: `rising_radius.py` handles the prograde cable (the casing ends at +omega_E r in the launch frame, the anchors hold it there, the Earth's work takes the sign of the casing's rotation) and logs when each LIM limit binds (the sampled coil voltage, current and site power; the first day from which the voltage stays at 100 kV, the days at the current cap and at the power cap), and `tools/rising_radius_study.py` runs and tabulates the study. The default fixed-radius run is untouched: the edits are inside the `if rising` blocks of `lim_simulation.py` (the diff shows no other line) and the 26 fast reference tests pass. The 8 MW rising run of Sections 3 to 9 stays the reference.

### 10.1 The 16 MW per site rerun

`lim_config.configure(MAX_SITE_POWER=16e6)` with everything else as the 8 MW run, on the full time steps (1 s on the first day, then 10 s). Log and files in `runs/rising_radius_16MW/`.

| Quantity | Fixed radius, 16 MW (`lim_simulation_16MW.txt`, the book's Table 9.3 column) | Rising radius, 16 MW |
|---|---|---|
| Deployment time | 427.3 days (14.0 months) | 573.9 days (18.9 months, 1.57 years): phase 1 427.9 days, phase 2 146.1 days |
| The casing at ground speed, phase 2 begins | day 427.3 (the end) | day 427.9 at 666.66 km, cable 8,207.3 m/s, 682.3 GN |
| Phase 2: duration, anchor reaction, the Earth's share | none | 146.1 days, 44.8 MN mean (50.4 kN per station), 45.2 MN peak; 3.99e21 kg m^2/s, 1.22 % of the cable's angular momentum |
| LIM work | 18.8 EJ | 23.864 EJ (268.3 TJ per site); the balance of the two sums 8e-6 EJ |
| The 100 kV coil voltage arrives | 2.5 months | day 60.0 (2.0 months) with 1,709 N per LIM; at it for 514 of the 574 days |
| The current cap, 183.1 A | the first six weeks | 60.1 days |
| The 16 MW power cap | from 1.8 months, falls back to 6.9 MW | 8.1 days only (days 52.6 to 60.7); 16.00 MW peak at 1.80 months, 6.81 MW at the end |
| Peak thrust per LIM | 1,950 N at 1.7 months | 1,948.8 N at 49.5 days (1.65 months) |
| Peak thrust power, LIM power per LIM | 5.5 MW at 1.8 months | 5.50 MW at 1.82 months, 5.86 MW |
| Peak hysteresis per LIM | 8.6 kW at two months | 8.64 kW at 2.00 months |
| Peak cryogenic power per site, cold-side load, radiator width | | 2.79 MW, 33.5 kW, 5.43 m at 2.00 months |
| Cryogenic power at the end | | 1.465 MW |
| Peak cable and wall temperature | 382 K | 382.2 K, 371.4 K, day 1 |
| Thrust, site power at the end | | 250.9 N, 6.81 MW (the same end state as the 8 MW run) |

Phase 1 of the rising 16 MW run, 427.9 days, is within a day of the fixed-radius run's whole deployment, a coincidence of the two changes (the longer slowdown and the lower starting speed) against the stronger push. Phase 2 is the same 146.1 days as at 8 MW: once the voltage limit holds, the thrust per LIM is 250 N whatever the site is allowed to draw, so the extra power buys nothing after day 60. The quick-step run of the same case in the study gives 573.9 days and 23.864 EJ, the same figures. The book's "14 to 17 months" (9.8) becomes 18.9 to 21.3 months; the "427 days" of Table 9.2, Table 9.3, 9.5 and 9.7 becomes 574.

The hang table for the 16 MW run (`runs/rising_radius_16MW/rising_radius_hanging.txt` and `.csv`): the ring rises faster, so every date moves up by 70 to 100 days.

| Day | Ring km | 100 km below | Spreader arm, 550 km below | Longest hang above 150 km |
|---|---|---|---|---|
| 0 | 527.4 | 427.4 | -22.6 | 377.4 |
| 60 | 546.2 | 446.2 | -3.8 | 396.2 |
| 100 | 565.4 | 465.4 | 15.4 | 415.4 |
| 200 | 600.8 | 500.8 | 50.8 | 450.8 |
| 300 | 632.6 | 532.6 | 82.6 | 482.6 |
| 400 | 659.7 | 559.7 | 109.7 | 509.7 |
| 500 | 683.2 | 583.2 | 133.2 | 533.2 |
| 573.9 | 700.0 | 600.0 | 150.0 | 550.0 |

A 400 km hang clears 150 km from day 67 (day 98 at 8 MW), 450 km from day 198 (265), 500 km from day 363 (436), 550 km only at the end; the 550 km arm hung from day 0 would be below the ground for the first 66 days (97). The points 20, 50 and 100 km below the casing are never below 150 km.

### 10.2 The prograde deployment

`ORBITAL_RING_CABLE=prograde python lim_simulation.py --rising-radius`, full time steps, 8 MW per site, 208 turns; log and files in `runs/rising_radius_prograde/`. The ring is `ring_altitude.py`'s prograde case: cable 89,266 kg/m structural, 99,952 with hardware, 52.51 m^2, design speed 8,343.3 m/s eastward, design tension 656.4 GN (12.5 GPa on that area), relative speed 7,827 m/s. The brief named 865 GN for the design point; that figure is not the prograde cable's tension in `ring_altitude.py` (Chapter 3, L156, uses it for a cable of the prograde mass asked to absorb the retrograde momentum), and the run uses the code's 656.4 GN. The build orbit is the same 527.36 km, since the strain is the same 2.5 %: the ring is assembled prograde at 7,597.5 m/s and its casing slows eastward to 516 m/s while the cable goes to 8,343 m/s.

| Quantity | Prograde, rising radius | Retrograde, rising radius (Section 4) |
|---|---|---|
| Deployment time | 476.9 days (15.7 months): phase 1 385.1 days, phase 2 91.8 days | 647.1 days: 501.1 and 146.1 |
| Casing at ground speed | day 385.1 at 665.40 km, cable 8,281.2 m/s, 524.8 GN, 9.99 GPa | day 501.1 at 666.66 km, 8,207.3 m/s, 682.3 GN, 10.09 GPa |
| Phase 2: rise, cable gain, tension gain | 34.60 km, 62.1 m/s, 131.5 GN | 33.34 km, 62.7 m/s, 163.3 GN |
| Angular momentum from the Earth | 3.245e21 kg m^2/s, 1.236 % of the cable's final 2.625e23 (1.242 % of the ring's) | 3.991e21, 1.220 % |
| Anchor reaction, phase 2 | 57.9 MN mean, 65 kN per station; 58.4 MN peak | 44.8 MN, 50 kN |
| LIM work | 18.219 EJ (204.8 TJ per site) | 23.864 EJ |
| Cable kinetic, casing kinetic, potential, elastic | +26.422, -15.331, +7.010, +0.356 EJ: 18.456 EJ | +29.825, -15.331, +8.621, +0.459: 23.573 |
| The Earth's work through the anchors | +0.237 EJ (positive: the Earth pushes the casing eastward as it moves eastward) | -0.291 EJ |
| LIM work required, delta E minus the Earth's work; the balance | 18.219 EJ; 6e-6 EJ | 23.864; 6e-6 |
| Fixed-radius deployment energy (`ring_altitude`) | 14.595 EJ | 18.804 EJ |
| The 100 kV limit arrives; thrust then | day 346.7 (11.4 months), 364 N | day 353.2, 364 N |
| Thrust per LIM at the end; site power at the end | 321 N, 7.59 MW (the prograde relative speed is lower, 7,827 m/s, so the voltage limit leaves more thrust) | 251 N, 6.81 MW |
| Supply frequency at the end | 39.3 Hz | 44.2 Hz |
| Peak thrust, hysteresis, cryogenic power, cable temperature | 1,933 N at day 27.2; 5.05 kW; 2.20 MW; 382.2 K | 1,938 N at 27.6; 5.09 kW; 2.20 MW; 382.2 K |

The closed-form budget from the design point and the build orbit (Section 4's method with `ring_altitude.design(retrograde=False)`) gives the same 18.456 EJ of mechanical change (26.422 + 7.010 + 0.356 - 15.331), so the run and the closed form agree, and the book's 14.6 EJ against 18.8 EJ becomes 18.2 EJ against 23.9 EJ: the prograde deployment still takes a quarter less energy, and 170 fewer days, 477 against 647. The energy through the anchors changes sign between the two: the retrograde ring does 0.29 EJ of work on the Earth in its final tensioning, the prograde ring receives 0.24 EJ from it.

### 10.3 The turns variants at 16 MW

From the study's quick-step runs (Section 10.4): with 208 turns the 16 MW site reaches 100 kV on day 60.0 and deploys in 573.9 days; with 150 turns it reaches 100 kV on day 229.4, sits at the current cap for 233 days instead, and deploys in 689.0 days, 115 days slower (the book's fixed-radius figure is 77 days slower). At 8 MW the 150-turn coil never reaches 100 kV and deploys in 650.8 days, 3.7 days slower than 208 turns (the book: 26 days slower); 260 turns reaches 100 kV on day 119.7 and takes 866.9 days.

### 10.4 The deployment study

`tools/rising_radius_study.py --run` runs the cases (quick steps, five at a time, 2 to 9 minutes each), `--table` makes the table from their summaries: `runs/rising_radius_study/rising_radius_study.md`, `.csv`, `.svg|png`, one folder per case with its log, daily table, series and summary. The coil voltage limit is 100 kV in every case. The baseline is the 8 MW, 208-turn, two-LIM site of the reference run. "100 kV from day" is the first day from which the coil voltage stays at its limit; "days at I cap" and "days at P cap" are the days the current cap (183.1 A, the worst tape at half its measured limit) and the site power cap bind. "Plates" is the LIM reaction plate mass on the cable (936 kg/m per plate; the baseline's two plates are 1,872 kg/m, 936 t per 500 m site); "cable" is the cable with hardware after the quadratic has resized it for the plates; "tape" is the HTS tape of the site's LIMs (2 layers of 3 mm tape, 382 km and 1.63 t per LIM at 208 turns, `lim_hysteresis.tape_mass_per_lim`); "coil area" is the coil cryostat surface per site. The cryogenic power is per site: its peak and its value at the end of the run.

| case | MW/site | turns | LIMs/site | days | phase 1 d | phase 2 d | 100 kV from day | thrust/LIM at 100 kV N | days at I cap | days at P cap | LIM work EJ | peak site MW | peak plate K | cryo peak MW | cryo end MW | plates kg/m cable | cable kg/m | tape km/site | tape t/site | coil area m2/site |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| base_8MW | 8 | 208 | 2 | 647.1 | 501.0 | 146.1 | 353.2 | 364 | 28.0 | 343.3 | 23.864 | 8.00 | 382.2 | 2.20 | 1.465 | 1872 | 125,691 | 764 | 3.3 | 3101 |
| P12MW | 12 | 208 | 2 | 572.7 | 426.6 | 146.1 | 90.6 | 999 | 45.2 | 56.1 | 23.864 | 12.00 | 382.2 | 2.57 | 1.465 | 1872 | 125,691 | 764 | 3.3 | 3101 |
| P16MW | 16 | 208 | 2 | 573.9 | 427.9 | 146.1 | 60.0 | 1709 | 60.1 | 8.1 | 23.864 | 16.00 | 382.2 | 2.79 | 1.465 | 1872 | 125,691 | 764 | 3.3 | 3101 |
| P24MW | 24 | 208 | 2 | 575.9 | 429.9 | 146.1 | 58.8 | 1902 | 59.3 | 0.0 | 23.864 | 17.37 | 382.2 | 2.82 | 1.465 | 1872 | 125,691 | 764 | 3.3 | 3101 |
| N150 | 8 | 150 | 2 | 650.8 | 528.8 | 122.0 | never |  | 115.9 | 542.8 | 23.864 | 8.00 | 331.2 | 2.24 | 1.585 | 1872 | 125,691 | 551 | 2.4 | 2236 |
| N260 | 8 | 260 | 2 | 866.9 | 638.3 | 228.6 | 119.7 | 627 | 15.5 | 115.4 | 23.864 | 8.00 | 425.9 | 2.10 | 1.509 | 1872 | 125,691 | 955 | 4.1 | 3876 |
| L4 | 8 | 208 | 4 | 642.5 | 510.8 | 131.8 | never |  | 10.4 | 632.4 | 23.899 | 8.00 | 453.5 | 2.62 | 1.947 | 3744 | 127,427 | 1528 | 6.5 | 6202 |
| L6 | 8 | 208 | 6 | 703.0 | 556.1 | 146.9 | never |  | 0.0 | 697.6 | 23.935 | 8.00 | 492.6 | 3.16 | 2.435 | 5616 | 129,166 | 2291 | 9.8 | 9303 |
| N150_16MW | 16 | 150 | 2 | 689.0 | 507.4 | 181.6 | 229.4 | 617 | 232.8 | 0.0 | 23.864 | 12.49 | 331.2 | 3.23 | 2.047 | 1872 | 125,691 | 551 | 2.4 | 2236 |
| P16MW_L4 | 16 | 208 | 4 | 315.9 | 242.5 | 73.5 | 148.8 | 416 | 18.4 | 137.9 | 23.899 | 16.00 | 453.5 | 3.28 | 2.097 | 3744 | 127,427 | 1528 | 6.5 | 6202 |
| P24MW_L4 | 24 | 208 | 4 | 288.6 | 215.1 | 73.5 | 43.0 | 1107 | 25.0 | 19.7 | 23.899 | 24.00 | 453.5 | 3.94 | 2.097 | 3744 | 127,427 | 1528 | 6.5 | 6202 |
| P24MW_L6 | 24 | 208 | 6 | 210.3 | 160.9 | 49.3 | 94.4 | 435 | 5.3 | 84.6 | 23.935 | 24.00 | 492.6 | 4.72 | 2.729 | 5616 | 129,166 | 2291 | 9.8 | 9303 |
| P16MW_N150_L4 | 16 | 150 | 4 | 310.4 | 252.5 | 57.8 | never |  | 54.6 | 263.4 | 23.899 | 16.00 | 385.8 | 3.62 | 2.423 | 3744 | 127,427 | 1102 | 4.7 | 4472 |
| P24MW_N150_L6 | 24 | 150 | 6 | 205.8 | 167.7 | 38.2 | never |  | 35.3 | 171.7 | 23.935 | 24.00 | 425.8 | 4.73 | 3.259 | 5616 | 129,166 | 1652 | 7.1 | 6709 |

The extra hardware against the baseline, per site and for the ring:

| Case | Plates on the cable | Cable with hardware | Tape per site | Cryogenic power at the end | Radiator width, peak | Anchor reaction in phase 2 (mean, per station) |
|---|---|---|---|---|---|---|
| baseline, 2 LIMs, 208 turns | 1,872 kg/m, 936 t per site | 125,691 kg/m | 3.3 t (764 km) | 1.47 MW | 4.28 m | 44.8 MN, 50 kN |
| 4 LIMs, 208 turns | 3,744 kg/m, 1,872 t per site | 127,427 kg/m, +1,735 kg/m, +77 Mt of cable | 6.5 t (1,528 km) | 1.95 MW at 8 MW, 2.10 at 16 and 24 MW | 5.09 m (8 MW) | 90.2 MN, 101 kN (16 MW) |
| 6 LIMs, 208 turns | 5,616 kg/m, 2,808 t per site | 129,166 kg/m, +3,474 kg/m, +155 Mt of cable | 9.8 t (2,291 km) | 2.44 MW at 8 MW, 2.73 at 24 MW | 6.14 m (8 MW), 9.18 m (24 MW) | 136.0 MN, 153 kN (24 MW) |
| 4 LIMs, 150 turns | 3,744 kg/m | 127,427 kg/m | 4.7 t (1,102 km) | 2.42 MW (16 MW) | | |
| 6 LIMs, 150 turns | 5,616 kg/m | 129,166 kg/m | 7.1 t (1,652 km) | 3.26 MW (24 MW) | | |

Reading the table:

- **Site power alone** saves 73 days whatever the power: 12, 16 and 24 MW all give 573 to 576 days, because the run is at the 100 kV limit from day 60 to 90 on, and from there the thrust is set by the voltage, not the power. At 24 MW the site never draws more than 17.4 MW. Phase 2 is the same 146 days in every two-LIM, 208-turn case: it runs at the voltage limit with 250 N per LIM regardless of the power. Per unit of hardware, the 12 MW step is the best single change in the table, 18.6 days per MW of site power, and it adds nothing to the ring; the 16 MW step buys 9.1 days per MW, the 24 MW step 4.4.
- **Turns alone** do not help. 150 turns never reaches 100 kV but loses thrust per ampere and sits at the current cap for 116 days: 651 days, 4 days slower than 208, with 0.9 t less tape per site. 260 turns reaches 100 kV on day 120 with 627 N and takes 867 days.
- **LIMs per site alone**, at 8 MW, do not help either: four LIMs give 643 days (the 8 MW cap binds for 632 of them; the voltage limit is never reached), six LIMs 703 days, because the extra LIMs' cryogenic load (2.4 MW at the end against 1.5) comes out of the same 8 MW and the cable is 3,474 kg/m heavier. Their first-day plate temperature is 453 K and 493 K against the 500 K plate limit (the baseline's 382 K).
- **Power and LIMs together** is where the days are. 16 MW with four LIMs: 316 days (243 of phase 1, 74 of phase 2), for 3.3 t more tape per site, 1,735 kg/m more cable (77 Mt), 0.6 MW more cryogenic power per site at the end, and an anchor reaction of 90 MN (101 kN per station) in phase 2. 24 MW with four LIMs: 289 days. 24 MW with six LIMs: 210 days (161 and 49), for 6.5 t more tape, 155 Mt more cable, 1.3 MW more cryogenic power per site, a 9.2 m radiator, and 136 MN on the anchors (153 kN per station). With 150 turns the same two cases give 310 and 206 days with 1.8 and 2.7 t less tape than at 208 turns, and the shortest phase 2 of the table, 58 and 38 days, because the lower-turn coils are not at the voltage limit in phase 2 and push with 300 to 340 N per LIM instead of 250. Per tonne of tape per site these combined cases buy 100 to 230 days.
- **Phase 2 under a month**: no case in the table. Phase 2 scales with the thrust at the voltage limit times the LIMs per site: 146 days with two LIMs, 73.5 with four, 49.3 with six at 208 turns; 58 and 38 days with four and six at 150 turns. Under 30 days would take about ten LIMs per site at 208 turns or eight at 150 turns, at 100 kV, or a higher coil voltage. The anchor reaction in phase 2 is the LIM thrust, so it grows with the LIM count: 45 MN with two LIMs, 90 with four, 136 with six.
- The LIM work is the same 23.864 EJ in every two-LIM case, as it must be (the mechanical energy change is set by the design point and the build orbit), and 23.90 and 23.94 EJ with the heavier four- and six-plate cables.

Two cautions on the multi-LIM cases. The model puts the extra LIMs on extra reaction plates per side of the cable (`N_PLATES_PER_SIDE`, the knob the brief named), each 200 mm by 1.2 m of gamma-TiAl; where a second and third plate sit on the 8.2 m cable, and how two LIMs share a side of a 500 m site when each is 300 m long, are layout questions the model does not ask. And the first-day plate temperature of the six-LIM site, 493 K, is 7 K under the 500 K at which the controller cuts the current, so that case has no thermal margin on the first day.

The study ran on the quick time steps (10 s on the first day, 100 s after; the reference run uses 1 and 10 s). The 8 MW baseline on the quick steps gives 647.1 days and 23.864 EJ, the same as the full-step run; the 16 MW case on the full steps is in Section 10.1 for the comparison.

### 10.5 Files and figures of round 2

Under `runs/` (copies in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/`):

- `rising_radius_16MW/`: `rising_radius_lim_simulation.txt`, `rising_radius_daily.csv`, `rising_radius_deployment_tension.csv`, `rising_radius_summary.json`, `rising_radius_hanging.txt|csv`, `figures_rising_radius/` (altitude, tension, speeds, anchor reaction, energy, the Figure 3.2 redraw).
- `rising_radius_prograde/`: the same set for the prograde run.
- `rising_radius_study/`: `rising_radius_study.md|csv|svg|png` (the table and the bar chart of phase 1 and phase 2 days with the tape per site), `_run_console.txt`, and one folder per case with its log, daily table, series and summary.

Reproduce: `python -c "import lim_config as c; c.configure(MAX_SITE_POWER=16e6); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; runpy.run_path('lim_simulation.py', run_name='__main__')"` with `ORBITAL_RING_OUT=runs/rising_radius_16MW`; `ORBITAL_RING_CABLE=prograde python lim_simulation.py --rising-radius` with `ORBITAL_RING_OUT=runs/rising_radius_prograde`; `python tools/rising_radius_study.py --run --jobs=5` then `--table`; `python tools/rising_radius_hanging.py runs/rising_radius_16MW/rising_radius_daily.csv`; the figures with `ORBITAL_RING_OUT=<folder> python tools/fig_rising_radius.py` and `tools/fig03_02_deployment_tension.py --rising --out=<folder>/figures_rising_radius`.

### 10.6 Method notes of round 2

- The prograde case needed no structural change: the launch frame keeps the cable's direction positive, and the casing's end state, the phase 2 constraint and the Earth's work carry a sign, -1 for the retrograde cable and +1 for the prograde.
- "100 kV from day" is the first sampled day from which the coil voltage is at 99.9 % of the limit for 95 % of the remaining samples; "days at the current cap" counts samples within 1 % of 183.1 A; "days at the power cap" samples within 2 % of the site's limit. The log's `volts_max` time (the first sample at the rounded maximum) can differ by a month from the first.
- The study ran on the quick time steps to keep thirteen cases to half an hour; the one case run both ways (16 MW) agrees to 0.1 day and 1e-6 EJ.
- The cases with more LIMs per site set `ORBITAL_RING_HW_KGM` to the plate mass of that case so that `ring_altitude.py` resizes the cable for it; the design point of such a run is that heavier ring's (its cable speed and tension differ in the fourth digit), and the LIM work rises by the lift and spin of the extra plates, 0.035 EJ per plate pair.

## 11. Round 3 (2026-11-06): the deployment options table at full time steps

Asked for in `Vol_III/Claude_Code_instructions_deployment_options_table_2026-11-06.md`. The design stays at 8 MW per site, 208 turns and two LIMs per site with the rising radius (647 days); Chapter 9 gets a section on shortening the deployment with a table of the options, and that table needs the full time steps (1 s on the first day, then 10 s), not the quick steps of the Section 10.4 study. `tools/rising_radius_options.py --run` ran the eight cases of the brief (seven runs, 8 to 29 minutes each, seven at a time; the 16 MW case is the full-step run of Section 10.1 copied in), `--table` made the table from their summaries and daily tables: `runs/rising_radius_options/rising_radius_options.md`, `.csv`, one folder per case. The cases are the study's, with the study's names and overrides, so the cases with two and three reaction plates per side set `ORBITAL_RING_HW_KGM` as before and the quadratic sizes the cable for them. The default fixed-radius run is untouched (the code change of this round is the new tool and one line of `tools/rising_radius_study.py`, the study folder from the environment); the fast reference tests pass (28 passed, 11 skipped).

The full-step base case reproduces the Section 4 run to every digit (647.1141 days, 23.86409 EJ, 44.787 MN), which is also the check that the runs used the full steps: the quick-step study gave 647.1088 days. The other cases agree with their quick-step runs to 0.04 days or better, so the Section 10.4 table stands; the numbers below are the full-step ones.

### 11.1 The table

Case numbers as in the brief. Months are days over 30.44. "100 kV from day" is the first day from which the coil voltage stays at the limit (Section 10.6); "never" means the site power cap binds to the end instead. The ring total is the peak site power times the 88,946 sites; the panel is the deployment array that supplies it at Chapter 8's 291 W/m^2 averaged around the ring, 27,500 m^2 and 55 m of width per site at 8 MW, scaled with the peak site power.

**Deployment and power**

| # | Case | Days | Months | Phase 1 d | Phase 2 d | 100 kV from day | LIM work EJ | Peak site MW | Ring GW | Panel m^2 per site | Panel width m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 8 MW, 208 turns, 2 LIMs (the design) | 647.1 | 21.3 | 501.1 | 146.1 | 353.2 | 23.864 | 8.00 | 712 | 27,500 | 55 |
| 2 | 12 MW, 208, 2 | 572.7 | 18.8 | 426.6 | 146.1 | 90.6 | 23.864 | 12.00 | 1,067 | 41,250 | 83 |
| 3 | 16 MW, 208, 2 | 573.9 | 18.9 | 427.9 | 146.1 | 60.0 | 23.864 | 16.00 | 1,423 | 55,000 | 110 |
| 4 | 8 MW, 208, 4 | 642.6 | 21.1 | 510.8 | 131.8 | never | 23.899 | 8.00 | 712 | 27,500 | 55 |
| 5 | 16 MW, 208, 4 | 315.9 | 10.4 | 242.5 | 73.5 | 148.8 | 23.899 | 16.00 | 1,423 | 55,000 | 110 |
| 6 | 16 MW, 150, 4 | 310.4 | 10.2 | 252.5 | 57.8 | never | 23.899 | 16.00 | 1,423 | 55,000 | 110 |
| 7 | 24 MW, 208, 6 | 210.3 | 6.9 | 161.0 | 49.3 | 94.4 | 23.935 | 24.00 | 2,135 | 82,500 | 165 |
| 8 | 24 MW, 150, 6 | 205.9 | 6.8 | 167.7 | 38.2 | never | 23.935 | 24.00 | 2,135 | 82,500 | 165 |

**Thermal, and the hardware on the cable with the quadratic's design point for it**

The plate temperature is the peak of the LIM reaction plate against the 500 K limit at which the controller cuts the current. The cryogenic power is per site, its peak and its value at the end. The levitation plates are the same 8,814 kg/m in every case (four bearing plates); the LIM reaction plates are 936 kg/m each, two per plate pair. "Structural" is the cable's structural mass from `ring_altitude.cable_mass_structural` with that hardware, "cable" the cable with hardware, and the total is that times the 44,476 km ring against the design's 5,590 Mt. The design-point speed and tension are `ring_altitude.design(m_hw=...)`; each run ends at exactly these (the run's final cable speed and tension agree with them to all printed digits).

| # | Case | Plate K | Cryo peak MW | Cryo end MW | LIM plates kg/m | Hardware vs design kg/m | Structural kg/m | Cable kg/m | Cable Mt | vs design Mt | Cable area m^2 | Design v_cable m/s | Relative speed m/s | Design tension GN |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 8 MW, 208, 2 | 382.2 | 2.20 | 1.46 | 1,872 | 0 | 115,005 | 125,691 | 5,590 | 0 | 67.65 | 8,270.0 | 8,786.2 | 845.6 |
| 2 | 12 MW, 208, 2 | 382.2 | 2.57 | 1.46 | 1,872 | 0 | 115,005 | 125,691 | 5,590 | 0 | 67.65 | 8,270.0 | 8,786.2 | 845.6 |
| 3 | 16 MW, 208, 2 | 382.2 | 2.79 | 1.46 | 1,872 | 0 | 115,005 | 125,691 | 5,590 | 0 | 67.65 | 8,270.0 | 8,786.2 | 845.6 |
| 4 | 8 MW, 208, 4 | 453.5 | 2.62 | 1.95 | 3,744 | +1,872 | 114,869 | 127,427 | 5,667 | +77 | 67.57 | 8,259.6 | 8,775.7 | 844.6 |
| 5 | 16 MW, 208, 4 | 453.5 | 3.28 | 2.10 | 3,744 | +1,872 | 114,869 | 127,427 | 5,667 | +77 | 67.57 | 8,259.6 | 8,775.7 | 844.6 |
| 6 | 16 MW, 150, 4 | 385.8 | 3.62 | 2.42 | 3,744 | +1,872 | 114,869 | 127,427 | 5,667 | +77 | 67.57 | 8,259.6 | 8,775.7 | 844.6 |
| 7 | 24 MW, 208, 6 | 492.3 | 4.72 | 2.73 | 5,616 | +3,744 | 114,735 | 129,166 | 5,744 | +155 | 67.49 | 8,249.4 | 8,765.6 | 843.6 |
| 8 | 24 MW, 150, 6 | 425.8 | 4.73 | 3.26 | 5,616 | +3,744 | 114,735 | 129,166 | 5,744 | +155 | 67.49 | 8,249.4 | 8,765.6 | 843.6 |

The design point moves with the hardware, through the quadratic. Each extra plate pair, 1,872 kg/m, takes 136 kg/m off the structural cable (115,005 to 114,869 to 114,735 kg/m) and puts 1,735 kg/m on the cable with hardware (125,691 to 127,427 to 129,166), 77 Mt per pair on the ring. The cable speed falls 10.4 m/s per pair (8,270.0 to 8,259.6 to 8,249.4 m/s; the relative speed to the casing the same, 8,786 to 8,776 to 8,766 m/s) and the design tension 1.0 GN per pair (845.6 to 844.6 to 843.6 GN), since a heavier cable carries the casing's momentum deficit with less excess speed and the cross-section follows the structural mass (67.65 to 67.57 to 67.49 m^2). Phase 2 begins at 666.66 km and 682.3 GN with two LIMs, 666.32 km and 679.8 GN with four, 665.98 km and 677.4 GN with six. The LIM work rises 0.035 EJ per pair with the lift and spin of the plates: 23.864, 23.899, 23.935 EJ.

**Hardware on the casing, the anchors in phase 2, the hang schedule**

The tape is the HTS tape of the site's LIMs (two layers of 3 mm tape, `lim_hysteresis.tape_mass_per_lim`: 382 km and 1.63 t per LIM at 208 turns, 275 km and 1.18 t at 150). The model has no LIM stator mass: it carries the coil's tape length and the cryostat surface per site (the coil area), not the stator's structure, so that column is not in the table. The cryogenic plant per site is sized by the peak cryogenic power of the second column group. The anchor reaction is the mean over phase 2 (the peak is within 1 % of it) and per station is over the 889 stations at 50 km. The hang columns are from the daily tables: the first day on which a hang of that length, hung from the casing, is above the 150 km floor of Section 7; a 550 km hang (the spreader arm as hung during the deployment) clears 150 km only at the end of the deployment in every case, by construction, since it is the arm's length from the 700 km ring.

| # | Case | Tape per site t (km) | Coil area m^2 per site | Cryo plant MW per site | Anchor reaction MN (peak) | Per station kN | Thrust per LIM in phase 2 N | 450 km hang clears day | 550 km hang clears day |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 8 MW, 208, 2 | 3.3 (764) | 3,101 | 2.20 | 44.8 (45.2) | 50 | 251 | 265 | 647 |
| 2 | 12 MW, 208, 2 | 3.3 (764) | 3,101 | 2.57 | 44.8 (45.2) | 50 | 251 | 197 | 573 |
| 3 | 16 MW, 208, 2 | 3.3 (764) | 3,101 | 2.79 | 44.8 (45.2) | 50 | 251 | 198 | 574 |
| 4 | 8 MW, 208, 4 | 6.5 (1,528) | 6,202 | 2.62 | 50.3 (50.5) | 57 | 141 | 271 | 643 |
| 5 | 16 MW, 208, 4 | 6.5 (1,528) | 6,202 | 3.28 | 90.2 (90.9) | 101 | 253 | 128 | 316 |
| 6 | 16 MW, 150, 4 | 4.7 (1,102) | 4,472 | 3.62 | 114.5 (115.0) | 129 | 322 | 147 | 310 |
| 7 | 24 MW, 208, 6 | 9.8 (2,291) | 9,303 | 4.72 | 136.0 (137.1) | 153 | 254 | 85 | 210 |
| 8 | 24 MW, 150, 6 | 7.1 (1,652) | 6,709 | 4.73 | 175.6 (176.3) | 198 | 329 | 99 | 206 |

### 11.2 The hang schedule of cases 2 and 5, for Chapter 12

From `tools/rising_radius_hanging.py` on each case's daily table (`rising_radius_hanging.txt` and `.csv` in the case folder), beside the design's schedule of Section 7:

| Hang below the casing | Design, 8 MW (647 days) | Case 2, 12 MW (573 days) | Case 5, 16 MW and four LIMs (316 days) |
|---|---|---|---|
| 400 km clears 150 km from day | 98 | 74 | 49 |
| 450 km | 265 | 197 | 128 |
| 500 km | 436 | 362 | 211 |
| 550 km, the spreader arm | 647, the end | 573, the end | 316, the end |
| 600 km | never | never | never |

### 11.3 Reading, for Chapter 9

- **Site power alone** buys 74 days and nothing more. At 12 MW the run reaches the 100 kV coil limit on day 91 with 999 N per LIM and deploys in 573 days; at 16 MW the limit comes on day 60 with 1,709 N, and the deployment is 574 days, a day longer, because the push comes earlier but the voltage limit sets the thrust from there on. Phase 2 is the same 146 days in cases 1 to 3: the thrust at 100 kV is 251 N per LIM whatever the site may draw. The design's 8 MW site is at its power cap for 343 days (day 27 to 370) and reaches 100 kV only on day 353. The 12 MW case adds nothing to the ring; it widens the deployment panel from 55 to 83 m and the ring's draw from 712 to 1,067 GW.
- **More LIMs alone**, at 8 MW, buy nothing: four LIMs give 643 days, because the 8 MW cap binds from day 10 to the end (the coil voltage is 75 kV at the end, never at the limit) and the extra LIMs' cryogenic load, 1.95 MW at the end against 1.46, comes out of the same 8 MW; the 141 N per LIM at the end is less than the design's 251. It costs 77 Mt of cable and 3.3 t more tape per site, and the first-day plate temperature is 453 K.
- **Power and LIMs together** are where the days are. 16 MW with four LIMs, case 5: 316 days (243 of phase 1, 73.5 of phase 2), 100 kV from day 149; for 3.3 t more tape and 3,100 m^2 more coil cryostat per site, 77 Mt more cable, 1.1 MW more cryogenic plant per site (3.3 MW peak), and an anchor reaction of 90 MN in phase 2, 101 kN per station against the design's 50. With 150 turns instead, case 6: 310 days (252.5 and 57.8), never at 100 kV (81 kV at the end), the 16 MW cap from day 47 to the end; the lower-turn coil is at the power cap in phase 2 and pushes with 322 N per LIM instead of 253, so phase 2 is 16 days shorter and the anchors carry 114.5 MN, 129 kN per station; 1.8 t less tape than case 5 and a cooler plate, 386 K. 24 MW with six LIMs, case 7: 210 days (161 and 49), 100 kV from day 94; 6.5 t more tape and 6,200 m^2 more coil area per site than the design, 155 Mt more cable, 4.7 MW of cryogenic plant per site, a 9.2 m radiator, 136 MN on the anchors (153 kN per station), and a first-day plate temperature of 492 K, 8 K under the 500 K limit. With 150 turns, case 8: 206 days (168 and 38), never at 100 kV (82 kV at the end), the 24 MW cap from day 34 to the end, 329 N per LIM in phase 2, 176 MN on the anchors (198 kN per station), 426 K on the plate, 7.1 t of tape per site.
- **Where the 100 kV limit sits.** Day 353 in the design, 91 at 12 MW, 60 at 16 MW, 149 with 16 MW and four LIMs, 94 with 24 MW and six; never in the 150-turn cases and in the four-LIM 8 MW case, where the power cap binds to the end. The 208-turn cases end phase 2 at the voltage limit with 251 to 254 N per LIM; the 150-turn cases end it at the power cap with 322 to 329 N.
- **Phase 2 under a month.** No case in the table. At 208 turns and 100 kV the phase 2 thrust is 251 to 254 N per LIM, so phase 2 scales with the LIM count: 146 days with two, 73.5 with four, 49.3 with six. By that scaling ten LIMs per site would give 29 days, with 16 t of tape per site and some 230 MN on the anchors, 260 kN per station; this is the scaling of the table, not a run. At 150 turns the thrust at the power cap, 322 to 329 N, gives 58 days with four LIMs and 38 with six, so eight LIMs at about 32 MW per site would be the other route, 29 days by the same scaling. The other lever is the coil voltage: a higher limit than 100 kV raises the thrust per LIM at the limit and shortens phase 2 for every case that ends there, but no run of this round changed it. The anchor reaction in phase 2 is the LIM thrust, so every route to a short phase 2 puts a proportionally larger force through the anchor stations.
- The LIM work is 23.86 EJ in every two-LIM case, as it must be, and 23.90 and 23.94 EJ with the heavier four- and six-plate cables.

The two cautions of Section 10.4 hold: the model puts the extra LIMs on extra reaction plates per side of the 8.2 m cable without asking where they sit or how two 300 m LIMs share a side of a 500 m site, and the six-LIM site has 8 K of thermal margin on the first day (case 7) or 74 K (case 8).

### 11.4 Files

Under `runs/rising_radius_options/` (copies in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/rising_radius_options/`): `rising_radius_options.md` and `.csv` (the one-row-per-case table with every column above), `_run_console.txt`, and one folder per case with `rising_radius_lim_simulation.txt`, `rising_radius_daily.csv`, `rising_radius_deployment_tension.csv`, `rising_radius_summary.json`; `rising_radius_hanging.txt|csv` in the folders of cases 2 and 5. Reproduce: `python tools/rising_radius_options.py --run --jobs=7` then `--table`; the hang schedule with `python tools/rising_radius_hanging.py runs/rising_radius_options/<case>/rising_radius_daily.csv`.

## 12. The figures redrawn from the rising-radius run (2026-11-06)

Asked for in `Vol_III/Claude_Code_instructions_figures_rising_radius_2026-11-06.md`: Figures 3.2, 14.1 and 14.2 of the book, drawn until now from the fixed-radius run, redrawn from the 8 MW rising-radius run (Section 4) with each figure's layout, fonts and size kept as the current file in `Vol_III/images/figures/`. The originals are in `Vol_III/images/figures_backup_2026-11-06/`. No book file was edited.

- **Figure 3.2** (`fig03_02_deployment_tension_700.svg` and `.png`): drawn by `Vol_III/figure_scripts_2026-11-06/gen_fig03_02_rising.py`, which is the v14 generator (`figure_scripts_2026-11-03/gen_fig03_02.py`, the layout Sonnet set in proof v14) with the series of the rising run: the orbital speed is the series' own first row (the build orbit, 7,598 m/s), the tension over the radius uses the radius of each row, the prograde cable is the prograde rising run (Section 10.2) dashed at the same slowdown, and the bottom panel marks the casing at rest (day 436), at ground speed with the anchors attached (day 501) and the design point (day 647), with the anchored phase shaded. The top panel's marks are now 7,081 m/s (prograde deployment complete, 656 GN), 7,598 m/s (casing at rest, 601 GN) and 8,114 m/s (retrograde deployment complete, 846 GN); the anchored phase is the vertical segment at the end of each curve, 682 to 846 GN at one casing speed, and is labelled. Middle panel: 214.5 kN/m, 95.0 kN/m and 119.5 kN/m at the end, "casing at rest, day 436"; the prograde curves stop at their anchored phase. The caption (Ch03 L65) and the alt text (L64) agree with the figure: 7,598 m/s, 516 m/s, 845.6 GN and 12.5 GPa, 95 and 214.5 kN/m, 647 days, day 436, day 501. Series: `runs/rising_radius_deployment_tension.csv` (4,053 rows) and `runs/rising_radius_prograde/rising_radius_deployment_tension.csv`.
- **Figures 14.1 and 14.2** (`fig14_01_cable_wall_temperature_700.png`, `fig14_02_cold_load_radiator_700.png`): `tools/ch14_thermal_series.py` gained `--rising` (the same run through `lim_simulation.py --rising-radius`, full time steps, 20 minutes) and wrote `runs/rising_radius_thermal/ch14_thermal_series_rising.csv` (279,554 rows at 200 s); `tools/fig14_figures.py` drew the two figures from it unchanged, so the layout, fonts and size are the Oct 27 ones and the month axis runs to 21.6 months (ticks to 22). The values: the cable peaks at 382.6 K on day 1, is 312 K at one month, below 260 K from month 1.6, 256.5 K at two months, 259 K through months 5 to 12, and falls to 255.3 K in the anchored phase (months 16.7 to 21.6); the wall is 7.3 K below it, 248 K at the end. The cold-side load starts at 16.2 kW, peaks at 26.4 kW (4.28 m) at 1.8 months, is 18.4 kW at month 12, and settles at 17.6 kW and 2.85 m from month 16.7 to the end. The figures' own labels say 383 K on day 1, 255 K and 248 K at the end; 26.4 kW, 4.28 m at 1.8 months; 17.6 kW, 2.85 m at the end.
- **Caption checks for Chapter 14** (not edited; for Paul): the caption of Figure 14.1 (L39) says the cable settles at 257 K for the rest of the 21 months including the last five against the Earth, and the alt text (L38) says it stays near 257 K to the end; the run has 259 K through the middle of the deployment and 255 K through the anchored phase, and the figure's end label reads 255 K. The caption of Figure 14.2 (L54) and the alt text (L53) say 16 kW and a 2.85 m radiator through the last months; the run has 17.6 kW with the 2.85 m (the radiator width is proportional to the load, 0.162 m per kW, so 16 kW would be 2.59 m), the figure's end label reads 17.6 kW, and 16 kW is the load on day 0. The alt text's "width axis runs from 0 to 5 m" is 0 to 7.3 m in the figure (the right axis follows the 0 to 45 kW left axis), as it was in the fixed-radius figure.
- Files: `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/rising_radius_thermal/` holds the console, the summary, the figures and the series decimated to one row per 10,000 s (`ch14_thermal_series_rising_1in50.csv`); the full series stays in `runs/rising_radius_thermal/`. Reproduce: `ORBITAL_RING_OUT=runs/rising_radius_thermal python tools/ch14_thermal_series.py runs/rising_radius_thermal/ch14_thermal_series_rising.csv --rising`, then `python tools/fig14_figures.py runs/rising_radius_thermal/ch14_thermal_series_rising.csv <out>`; Figure 3.2 with `python Vol_III/figure_scripts_2026-11-06/gen_fig03_02_rising.py runs/rising_radius_deployment_tension.csv runs/rising_radius_prograde/rising_radius_deployment_tension.csv --out=<out>`.
