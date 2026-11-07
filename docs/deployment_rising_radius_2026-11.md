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

**Superseded (2026-11-06, round 4):** this table was made with the October cable; the cable of the rising quadratic and the site power sweep that replaces this table are in Section 13. Kept for the record.

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

## 13. Round 4 (2026-11-06): the rising quadratic, the verification run and the site power sweep

Asked for in `Vol_III/Claude_Code_instructions_rising_quadratic_and_power_2026-11-06.md`. Paul's rulings: the cable is re-sized so that the deployment from the build orbit ends at 700 km and 12.5 GPa on internal momentum alone, with no anchored phase against the Earth; then the best site power is found for two LIMs per site and for four. This section supersedes the options table of Section 11 (its cable was the October one, which needs the 146-day anchored phase) and the study of Section 10.4; both stay for the record. The October design point and every result of Sections 1 to 12 come back with `ORBITAL_RING_QUADRATIC=fixed`.

### 13.1 The quadratic with the build orbit

`ring_altitude.py` sizes the cable by the force balance at the design radius R under a momentum constraint. The October form took everything to start at the circular speed of R, v_orbit = sqrt(GM / R) = 7,504.3 m/s. The real deployment starts at the build orbit r0 = R / (1 + sigma / E) = 6,905.5 km (527.36 km, 172.64 km below the ring), where the ring is assembled unstressed on the circular orbit at 7,597.5 m/s, and what is conserved as the ring rises is angular momentum. The speed that stands for the ring's angular momentum at R is the build orbit's angular momentum per unit mass brought to R, v_mom = sqrt(GM r0) / R = 7,412.2 m/s, and the momentum constraint becomes m_c v_c - m_load v_geo = (m_c + m_load) v_mom (the retrograde ring; the prograde with the sign of v_geo reversed, as the code had it). The force balance at R and the stress constraint are unchanged. With v_mom in place of v_orbit the cable's own momentum term no longer cancels, D = v_mom^2 - v_orbit^2 = -1.373e6 m^2/s^2, and the quadratic is

    (S - D) m_s^2 + (S m_hw - K - 2 D m_hw) m_s - (K m_hw + P + D m_hw^2) = 0,
    S = sigma / rho,  K = 2 m_load dv v_mom - m_load g_net R,  P = (m_load dv)^2,  dv = v_mom + v_geo,

which is the October form when D = 0 (`cable_mass_structural` computes the same floating-point expressions then, so the October numbers come back to the last digit: 115,005.0919 kg/m, 8,270.0145 m/s). `ORBITAL_RING_QUADRATIC` selects `rising` (the default since this round) or `fixed`; `ring_altitude.py` also carries `R_BUILD`, `V_BUILD`, `V_MOM`, `DELTA_V_CASING`, `cable_velocity()` and `deployment_energy()` (the rising form: the kinetic energy change of cable and casing from the build orbit, the potential energy of the lift and the elastic energy of the cable, which with no anchored phase is the LIM work), and `design()` takes the mode. The reference tests run under the fixed quadratic (`tests/conftest.py`) and a new test holds the rising design point; 29 pass.

| Quantity | Retrograde, rising quadratic (the design) | Retrograde, October (fixed) | Prograde, rising quadratic | Prograde, October |
|---|---|---|---|---|
| Momentum velocity v_mom | 7,412.2 m/s | 7,504.3 (v_orbit) | 7,412.2 | 7,504.3 |
| Casing term dv in the constraint | 7,928.4 m/s | 8,020.4 | 6,896.1 | 6,988.1 |
| Casing speed change, build orbit to ground speed | 8,113.7 m/s | 8,020.4 | 7,081.4 | 6,988.1 |
| Cable structural | 92,883.3 kg/m | 115,005.1 | 71,386.0 | 89,266 |
| Cable with hardware (10,686 kg/m) | 103,569.3 kg/m | 125,691.5 | 82,072.0 | 99,952 |
| Cable speed | 8,330.8 m/s | 8,270.0 | 8,420.5 | 8,343.3 |
| Speed past the casing | 8,847.0 m/s | 8,786.2 | 7,904.4 | 7,827 |
| Cable tension, stress | 683.0 GN, 12.5 GPa | 845.6 GN | 524.9 GN | 656.4 GN |
| Cable area, side | 54.64 m^2, 7.39 m | 67.65, 8.22 | 41.99, 6.48 | 52.51, 7.25 |
| Cable mass with hardware; ring mass | 4,606 Mt; 5,140 Mt | 5,590; 6,124 | 3,650; 4,184 | 4,445; 4,979 |
| Deployment energy (closed form) | 19.175 EJ: kinetic 11.569, lift 7.236, elastic 0.370 | 18.804 | 14.902: 8.727, 5.890, 0.285 | 14.595 |
| Casing weight | 95.0 kN/m | same | same | same |

These are the brief's check values (92,883 / 103,569 / 8,330.8 / 8,847 / 683.0 / 4,606 Mt) to the digit. The cable is 19 % lighter than October's (22,122 kg/m less structural mass, 984 Mt less on the ring), because the ring arrives at 700 km with less angular momentum than the October constraint assumed, so less cable mass is needed to carry the casing's momentum deficit at the same stress, and the cable runs 61 m/s faster. The deployment energy is 2 % higher than October's figure for the heavier cable, since the lift of the whole ring, 7.2 EJ, is now in it.

**The hardware model.** The code's cable hardware does not scale with the cable's side: `lim_config.py` makes the four bearing plates 0.20 m by 1.4 m of iron (8,814 kg/m) and the two LIM reaction plates 0.20 m by 1.2 m of gamma-TiAl (1,872 kg/m), widths set by the coils, and `bearing_mass.py` and `final_bearing.py` size the bearing from the casing's weight (95.0 kN/m, unchanged) and check that the plates fit on the cable face. So m_hw stays 10,686 kg/m and there is nothing to iterate; the only question is whether the plates still fit on the 7.39 m face (they took 8.22 m before), which the bearing reruns below answer.

**The scripts that follow.** Every script that imports `ring_altitude.py` sees the new design point. `run_all.py --no-lim` was rerun into `runs/rising_quadratic/all/` (the October logs stay in `runs/` and `tests/reference/`); the equation sources of `tools/generate_equations_latex.py` were written under both quadratics (no LaTeX on this machine, so `.tex` only; the book's SVGs are untouched). The outputs that change, from a diff of the headline lines against the October logs:

| Output (log of `run_all.py`, or tool) | Changes under the rising quadratic | What changed |
|---|---|---|
| `ring_altitude_700`, `_700_prograde`, `_250`, `_800` | yes, 16 of 31 lines | the cable's mass, speed, area, tension, ring mass and deployment energy at every altitude and direction (the table above for 700 km; 800 km: structural 92,128 kg/m, 8,280.6 m/s, 677.4 GN against 113,809, 8,220.6, 836.8; 250 km: 96,543 kg/m, 8,571.0 m/s, 709.9 GN against 120,863, 8,506.3, 888.7) |
| `anchor_700_i0` to `_i30` (`anchor_line_analysis.py`) | yes, 42 of 135 lines | the cable's speed and mass in the lateral force of the inclined ring: the total force per anchor 5.317 GN against 5.331 at 20 degrees, and the line angles, lengths, stresses and tapers that follow |
| `anchor_lines_deployment` (`tools/anchor_lines_deployment.py` on the verification run's series) | yes (new series) | the casing's rate of slowing and net weight day by day, the pull on a stowed line and the tip speed through the air, from the 510-day deployment |
| `bearing_design`, `bearing_mass`, `final_bearing` | yes, 84, 85 and 120 lines | the ring's row of each (cable 103,569 kg/m, speed past the casing 8,847 m/s, cable face 7.39 m) and everything that uses the face width: the bearing's plates still fit, 4.33 m of the 7.39 m underside and 1.90 m on top (`bearing_mass`), the plates' edges at 3.70 m from the centerline (`final_bearing`); the casing's weight, 95.0 kN/m, and so the track ratings, tape and the 114 dynamics cases are unchanged (smallest gap 50.2 mm in `final_bearing`) |
| `bearing_fem_check`, `compensation_fields` | timing lines only | the solve times; every number the same |
| `compensation_tracks` | yes, 156 of 480 lines | the ring's row and the cable's moment of inertia (`CABLE_SIDE`), so the track forces and timings of the lateral compensation |
| `j2_retrograde_i5` to `_i30`, `j2_prograde_i5` to `_i30` | yes, 121 of 320 lines | the cable's mass per meter and speed in the J2 torque and the precession figures (the free precession rate at the J2 torque 6.58 deg/day against 6.50 at 20 degrees) |
| `levitation` | yes, 37 of 204 lines | the eddy heating against the gap at the new relative speed, 8,847 m/s (the plate's service temperature reached at a 51.4 mm gap and melting at 33.4 mm, against 51.3 and 33.3 in October) |
| `sled_passage`, `sled_passage_limits`, `sled_passage_pulsed`, `sled_passage_spread` | yes, 1,370, 29, 559 and 320 lines | the cable's mass, side and speed in the sled passage: the design-point line and the tables that follow |
| `casing_6000` (`casing_6000.py --bearing` and `--table`) | yes | the Chapter 3 quadratic for the three hardware cases of the 6,000 kg/m casing (the ring rows); the deployment logs behind Table 16.2 are the October fixed-radius runs of `run_700_casing_6000.sh` and were not rerun |
| `tools/generate_equations_latex.py` (`runs/rising_quadratic/all/equations/*.tex`) | yes, 11 of 26 equation sources | 2.5 (the lateral kinetic energy, through the ring mass), 3.6, 3.7 (now with v_mom, 7,412.2 m/s, in place of v_orbit; the new `eq03_05b_vmom_700` gives v_mom), 3.8 (683.0 GN), 3.9 and 3.10 (the second moment and L_char of the 7.39 m cable), 4.3 to 4.5 (the cable's centrifugal and net outward acceleration, 9.805 and 1.849 m/s^2, the gap time), 4.8 (the apogee of the released cable); no LaTeX on this machine, so `.tex` only and the book's SVGs untouched |
| `power_deployment`, `power_ops`, `power_deployment_ring_mounted_8MW` and `_16MW` | no | the power model uses the ring's circumference and the array's altitude, not the cable |
| `debris_orbits`, `debris_simulation`, `altitude_trade`, `md_simulation`, `md_trade_study`, `lsm_*`, `tether_stress`, `freefall`, `rendezvous`, `climber_descent`, `atmospheric_loading`, `ch4_numbers` | no | they take the altitude, the orbital speed or the circumference from `ring_altitude.py`, not the cable (`md_trade_study` differs by one path line) |

The two LIM deployment runs of `run_all.py` (`lim_simulation`, `lim_simulation_16MW`, the fixed-radius ones) were not rerun: the deployment of the new cable is the rising-radius run of 13.2.

### 13.2 The verification run: the deployment with the new cable

`ORBITAL_RING_OUT=runs/rising_quadratic/verify python tools/ch14_thermal_series.py <csv> --rising` (the same run as `lim_simulation.py --rising-radius`, full time steps, with the thermal series written): 8 MW per site, 208 turns, two LIMs per site, the retrograde cable of the rising quadratic. Log, daily table, Figure 3.2 series, summary, hang table, figures and the thermal series in `runs/rising_quadratic/verify/`.

**The end state.** The run ends on the design tension on day 510.23 with the three conditions met within the solver's precision: the ring at 700.000 km (14 m above the design radius), the casing at 516.146 m/s eastward (the ground-synchronous speed to the step), the cable at 8,330.817 m/s (0.003 m/s under the design 8,330.820) and the tension at 682.966 GN, 12.500 GPa (0.00005 GN over the design). The casing reached ground speed on day 510.2257 and the tension its design value on day 510.2263, one 10 s time step later, so the anchored phase lasted one step: the ring rose 0.00 km in it, the Earth supplied 1.8e16 kg m^2/s of angular momentum (7e-8 of the ring's 2.70e23) and did -1.3e-6 EJ of work. With the October cable the same run needed 146 days of anchored phase and 3.99e21 kg m^2/s from the Earth (Section 4). The total angular momentum at the build equals that at the design point to all printed digits (2.69656e23, +0.000 %), which is the quadratic's constraint.

**The first-list set of the round-1 brief, for Chapters 3, 9, 11, 12 and 14** (October cable's rising-radius value in brackets):

| Quantity | Value |
|---|---|
| Deployment time | 510.2 days, 16.8 months, 1.40 years (647.1 days) |
| Phase 1, phase 2 | 510.2 days; none (one time step) (501.1 and 146.1 days) |
| The casing at rest | day 442 at 679.2 km, 600.7 GN (day 436) |
| The 100 kV coil limit arrives; thrust then | day 347.8, 364 N per LIM (day 353.2) |
| The current cap, 183.1 A | the first 27.7 days (28.0) |
| The 8 MW power cap | 338 days, day 26.8 to 364.7 (343 days, 27.1 to 370.4) |
| Peak thrust per LIM | 1,934 N on day 27.3 (1,938 N on day 27.6) |
| LIM work | 19.175 EJ, 215.6 TJ per site (23.864 EJ, 268.3 TJ) |
| Its parts: cable kinetic, casing kinetic, lift, elastic | +26.900, -15.331, +7.236, +0.370 EJ; sum 19.175; the balance of the loop and the mechanical change 6e-6 EJ (+29.825, -15.331, +8.621, +0.459, minus the Earth's 0.291) |
| The closed form, `ring_altitude.E_DEPLOY` | 19.1753 EJ (the run 19.1753; the 1e-5 EJ between them is `lim_config`'s 10,686.4 kg/m of hardware against `ring_altitude`'s 10,686) |
| Altitude profile | 527.36 km at the build; 531.8 on day 30, 542.0 on day 60, 556.1 on day 100, 592.8 on day 200, 629.6 on day 300, 665.4 on day 400, 681.7 on day 450, 697.0 on day 500, 700.0 on day 510 |
| Tension profile | 17 GN on day 30, 58 on day 60, 114 on day 100, 259 on day 200, 404 on day 300, 546 on day 400 (10.0 GPa), 601 at rest (day 442), 683 at the end |
| The casing's net weight | half of its final 95.0 kN/m on day 51; peak 96.7 kN/m on day 391; 95.0 kN/m at the end; the cable's excess force 191.5 kN/m at the end (214.5 with the October cable) |
| Thrust and site power at the end | 247 N per LIM, 6.77 MW, 44.5 Hz (251 N, 6.81 MW) |
| Anchor reaction in the deployment | none: the LIM thrust at the end, 43.7 MN over the ring, is the cable's, not the anchors'; the anchors are attached at the end at 700 km (44.8 MN mean over 146 days) |

**The hang table** (`tools/rising_radius_hanging.py --depths=100,377,550`; the 150 km floor of Section 11.1):

| Day | Ring km | 100 km below | 377 km below | 550 km below (the spreader arm) | Longest hang above 150 km |
|---|---|---|---|---|---|
| 0 | 527.36 | 427.4 | 150.4 | -22.6 (below the ground) | 377 |
| 30 | 531.8 | 431.8 | 154.8 | -18.2 | 382 |
| 60 | 542.0 | 442.0 | 165.0 | -8.0 | 392 |
| 100 | 556.1 | 456.1 | 179.1 | 6.1 | 406 |
| 200 | 592.8 | 492.8 | 215.8 | 42.8 | 443 |
| 300 | 629.6 | 529.6 | 252.6 | 79.6 | 480 |
| 400 | 665.4 | 565.4 | 288.4 | 115.4 | 515 |
| 500 | 697.0 | 597.0 | 320.0 | 147.0 | 547 |
| 510.2 | 700.0 | 600.0 | 323.0 | 150.0 | 550 |

A hang of 100 km is never below 427 km; 377 km is the longest hang that clears 150 km from day 0 (150.4 km); a 400 km hang clears it from day 83, 450 km from day 220, 500 km from day 356, and the 550 km spreader arm only at the end (it would be below the ground for the first 82 days).

**The thermal history** (`ch14_thermal_series.csv`, 220,000 rows at 200 s; the figures drawn from it by `tools/fig14_figures.py` are in `figures_rising_radius/` for the record, not for the book yet). It is the October cable's history (the LIM model and its losses do not depend on the cable mass; the deployment is 137 days shorter): the cable peaks at 382.6 K on day 1, is 367.7 K at half a month and 309.5 K at one month, below 260 K from month 1.6, 256.5 K at two months, 257.8 at three, 258.7 at five, and 255.1 K at the end (month 16.8); the wall 7.3 K below it, 247.9 K at the end. The cold-side load starts at 16.2 kW, peaks at 26.35 kW (4.27 m of radiator) at 1.73 months, is 20.0 kW at five months and 17.5 kW, 2.84 m, at the end. The cryogenic power per site peaks at 2.20 MW on day 52 and ends at 1.46 MW; the peak hysteresis 5.06 kW per LIM.

### 13.3 The site power sweep

`tools/rising_radius_power_sweep.py --run` ran the eighteen cases of the brief at the full time steps (1 s on the first day, then 10 s; 4 to 45 minutes each, five at a time), each `lim_simulation.py --rising-radius` with the cable of the rising quadratic, 100 kV, and for the four-LIM cases `ORBITAL_RING_HW_KGM` set to the two-plates-per-side hardware so that the quadratic sizes that cable (12,558 kg/m of hardware: structural 92,451 kg/m, 105,010 kg/m with hardware, 4,670 Mt, 8,318.2 m/s, 679.8 GN; the LIM work of that ring is 19.148 EJ, the closed form again). `--table` made the table, the knees and the figure: `runs/rising_radius_sweep/rising_radius_power_sweep.md`, `.csv`, `.svg|png`, one folder per case. The 8 MW two-LIM case is the verification run of 13.2 to the digit (510.23 days). No case has an anchored phase: every run ends on the design tension at 700 km with the casing at ground speed, as the quadratic requires. "P cap last day" is the last day the site power cap binds, "100 kV from day" the first day from which the coil voltage stays at the limit (the definitions of Section 10.6); "never" in the first means the cap never binds (the voltage limit holds the draw under it), in the second that the limit is never reached. The panel width is 55 m at 8 MW scaled with the peak site power (Chapter 8, 291 W/m^2).

| LIMs | Turns | MW/site | Days | Months | P cap last day | 100 kV from day | LIM work EJ | Peak site MW | Panel width m | Plate K | Cryo peak MW | Cryo end MW | Tape t/site | Thrust at the end N/LIM | Cable kg/m, Mt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 208 | 6 | 695.3 | 22.8 | 695.3 (the end) | never | 19.175 | 6.00 | 41 | 382.2 | 1.98 | 1.43 | 3.3 | 212 | 103,570; 4,606 |
| 2 | 208 | 8 | 510.2 | 16.8 | 364.7 | 347.8 | 19.175 | 8.00 | 55 | 382.2 | 2.20 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 10 | 453.8 | 14.9 | 174.9 | 165.6 | 19.175 | 10.00 | 69 | 382.2 | 2.39 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 12 | 437.0 | 14.4 | 94.7 | 89.2 | 19.175 | 12.00 | 83 | 382.2 | 2.57 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 14 | 437.7 | 14.4 | 63.6 | 62.3 | 19.175 | 14.00 | 96 | 382.2 | 2.73 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 16 | 438.3 | 14.4 | 59.9 | 59.1 | 19.175 | 16.00 | 110 | 382.2 | 2.79 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 20 | 440.3 | 14.5 | never | 57.9 | 19.175 | 17.40 | 120 | 382.2 | 2.82 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 2 | 208 | 24 | 440.3 | 14.5 | never | 57.9 | 19.175 | 17.37 | 119 | 382.2 | 2.82 | 1.46 | 3.3 | 247 | 103,570; 4,606 |
| 4 | 208 | 8 | 517.8 | 17.0 | 517.8 (the end) | never | 19.148 | 8.00 | 55 | 453.5 | 2.62 | 1.95 | 6.5 | 141 | 105,010; 4,670 |
| 4 | 208 | 12 | 318.1 | 10.5 | 318.1 (the end) | never | 19.148 | 12.00 | 83 | 453.5 | 2.97 | 2.07 | 6.5 | 234 | 105,010; 4,670 |
| 4 | 208 | 16 | 247.0 | 8.1 | 153.6 | 146.5 | 19.148 | 16.00 | 110 | 453.5 | 3.29 | 2.09 | 6.5 | 249 | 105,010; 4,670 |
| 4 | 208 | 20 | 225.8 | 7.4 | 75.9 | 71.8 | 19.148 | 20.00 | 138 | 453.5 | 3.61 | 2.09 | 6.5 | 249 | 105,010; 4,670 |
| 4 | 208 | 24 | 220.2 | 7.2 | 43.7 | 42.5 | 19.148 | 24.00 | 165 | 453.5 | 3.94 | 2.09 | 6.5 | 249 | 105,010; 4,670 |
| 4 | 208 | 32 | 221.2 | 7.3 | 35.0 | 34.5 | 19.148 | 32.00 | 220 | 453.5 | 4.63 | 2.09 | 6.5 | 249 | 105,010; 4,670 |
| 4 | 150 | 12 | 334.1 | 11.0 | 334.1 (the end) | never | 19.148 | 12.00 | 83 | 385.8 | 3.16 | 2.06 | 4.7 | 234 | 105,010; 4,670 |
| 4 | 150 | 16 | 255.2 | 8.4 | 255.2 (the end) | never | 19.148 | 16.00 | 110 | 385.8 | 3.61 | 2.42 | 4.7 | 320 | 105,010; 4,670 |
| 4 | 150 | 24 | 187.8 | 6.2 | 187.8 (the end) | 184.5 | 19.148 | 24.00 | 165 | 385.8 | 4.29 | 3.30 | 4.7 | 482 | 105,010; 4,670 |
| 4 | 150 | 32 | 250.6 | 8.2 | never | 110.8 | 19.148 | 25.91 | 178 | 385.8 | 5.63 | 3.24 | 4.7 | 208 | 105,010; 4,670 |

**The curves and the knees** (the knee: the lowest power within 2 % of the curve's floor):

- **Two LIMs per site, 208 turns**: 695, 510, 454, 437, 438, 438, 440, 440 days at 6, 8, 10, 12, 14, 16, 20, 24 MW; the floor is 437.0 days at 12 MW. Days saved per MW: 92.5 from 6 to 8 MW, 28.2 from 8 to 10, 8.4 from 10 to 12, then nothing (the curve rises 0.3 to 0.5 days per MW to 24 MW: the earlier arrival at the voltage limit costs a little, since the voltage-limited thrust is lower than the power-limited one). The knee is **12 MW, 437.0 days (14.4 months)**: 100 kV from day 89, the power cap binding to day 95, LIM work 19.175 EJ, peak site 12.00 MW, panel 83 m, plate 382 K, cryogenic 2.57 MW peak and 1.46 MW at the end, 3.3 t of tape per site, 247 N per LIM at the end. 10 MW is 3.8 % above the floor (454 days) and misses the 2 % band; it saves 56 days over 8 MW for 2 MW. Above 16 MW the site never draws more than 17.4 MW.
- **Four LIMs per site, 208 turns**: 518, 318, 247, 226, 220, 221 days at 8, 12, 16, 20, 24, 32 MW; floor 220.2 days at 24 MW. Days saved per MW: 49.9 from 8 to 12, 17.8 from 12 to 16, 5.3 from 16 to 20, 1.4 from 20 to 24, none to 32. The knee is **24 MW, 220.2 days (7.2 months)**: 100 kV from day 42.5, the power cap to day 43.7, LIM work 19.148 EJ, panel 165 m, plate 453.5 K on the first day, cryogenic 3.94 MW peak and 2.09 MW at the end, 6.5 t of tape per site, 249 N per LIM at the end; the cable 105,010 kg/m with hardware, 4,670 Mt (64 Mt more than the design's), 8,318.2 m/s, 679.8 GN. 20 MW is 2.5 % above the floor (225.8 days), just outside the band, with a 138 m panel and 3.61 MW of cryogenic plant. The knee is not above 24 MW; the 32 MW case buys nothing (221.2 days) and costs a 220 m panel and 4.63 MW of cryogenic plant per site.
- **Four LIMs per site, 150 turns**: 334, 255, 188, 251 days at 12, 16, 24, 32 MW; floor 187.8 days at 24 MW, the shortest deployment of the sweep. Days saved per MW: 19.7 from 12 to 16, 8.4 from 16 to 24, and 7.8 lost per MW from 24 to 32. The knee is **24 MW, 187.8 days (6.2 months)**: the power cap binds from day 70 to the end and the coil reaches 100 kV only in the last 3 days, so the thrust per LIM stays high, 482 N at the end at 23.7 MW; LIM work 19.148 EJ, panel 165 m, plate 385.8 K, cryogenic 4.29 MW peak and 3.30 MW at the end, 4.7 t of tape per site. The 32 MW case is slower, 250.6 days: with the higher cap the coil draws its full 183.1 A for 112 days, reaches 100 kV on day 111 at 682 N, and from there runs voltage-limited, which for the 150-turn coil means 531 N at 173 A in month 4 and 211 N at 127 A in month 8, against 650 N at 148 A and 501 N at 129 A in the power-limited 24 MW run at the same months; the site draws 22 MW falling to 12.4 MW. The voltage-limited operating point of this coil pushes less per ampere than its power-limited one (the controller's slip at the limit), so that case is a control-strategy result, not a physical optimum, and the 150-turn curve has a true minimum at 24 MW only in the sense that more cap lets the controller reach the limit sooner. The same effect is the 0.3 to 0.5 days per MW rise of the two-LIM curve above 12 MW.

**The two candidate optima**, as the brief asks, with Paul to pick:

| | Two LIMs per site, 208 turns, 12 MW | Four LIMs per site, 208 turns, 24 MW | (Four LIMs, 150 turns, 24 MW) |
|---|---|---|---|
| Deployment | 437.0 days, 14.4 months | 220.2 days, 7.2 months | 187.8 days, 6.2 months |
| 100 kV from day; power cap to day | 89; 95 | 42.5; 43.7 | 184.5; the end |
| LIM work | 19.175 EJ (215.6 TJ per site) | 19.148 EJ | 19.148 EJ |
| Peak site power; ring; panel width | 12.0 MW; 1,067 GW; 83 m | 24.0 MW; 2,135 GW; 165 m | 24.0 MW; 2,135 GW; 165 m |
| Plate temperature, first day | 382.2 K | 453.5 K | 385.8 K |
| Cryogenic plant per site (peak); at the end | 2.57 MW; 1.46 | 3.94 MW; 2.09 | 4.29 MW; 3.30 |
| HTS tape per site | 3.3 t (764 km) | 6.5 t (1,528 km) | 4.7 t (1,102 km) |
| Cable hardware; structural; with hardware; ring total | 10,686; 92,883; 103,570 kg/m; 4,606 Mt | 12,558; 92,451; 105,010 kg/m; 4,670 Mt (+64) | the same |
| Design cable speed; tension | 8,330.8 m/s; 683.0 GN | 8,318.2 m/s; 679.8 GN | the same |
| Thrust per LIM at the end; site power at the end | 247 N; 6.8 MW | 249 N; 12.7 MW | 482 N; 23.7 MW |
| Anchored phase | none | none | none |

Against the design of Section 13.2 (8 MW, two LIMs, 510 days), the two-LIM knee saves 73 days for 4 MW more per site and nothing on the ring; the four-LIM knee saves 290 days for 16 MW more per site, 3.3 t more tape and 3,100 m^2 more coil cryostat per site, 1.4 MW more cryogenic plant per site, a 165 m panel, 64 Mt more cable, and the 453 K first-day plate (47 K under the limit) and the layout questions of Section 10.4. The 150-turn four-LIM case at 24 MW is 32 days faster than the 208-turn one with 1.8 t less tape and a cooler plate, at 0.35 MW more cryogenic plant and a higher end-of-run draw.

### 13.4 Files, and what this round supersedes

- The options table of Section 11 and the study of Section 10.4 were made with the October cable and are superseded by the sweep of 13.3; they stay for the record. The deployment of Section 4 (647 days with 146 anchored) and the figures of Section 12 are the October cable's; the book's figures are not regenerated in this round (the brief: not yet).
- Under `runs/` (copies in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/`): `rising_quadratic/verify/` (the verification run: log, daily table, Figure 3.2 series, summary, hang table, the figures of `tools/fig_rising_radius.py`, `tools/fig03_02_deployment_tension.py --rising` and `tools/fig14_figures.py`, and the thermal series, in the copy decimated to one row per 10,000 s); `rising_quadratic/all/` (every log of `run_all.py --no-lim` under the rising quadratic, with the figures and CSV files the scripts write, and `equations/` with the `.tex` sources; `rising_quadratic/fixed_equations/` the same sources under the October quadratic for the diff); `rising_radius_sweep/` (the table, the figure, one folder per case).
- Reproduce: `python ring_altitude.py` (the design point; `ORBITAL_RING_CABLE=prograde` for the prograde cable; `ORBITAL_RING_QUADRATIC=fixed` for October's); `ORBITAL_RING_OUT=runs/rising_quadratic/verify python tools/ch14_thermal_series.py runs/rising_quadratic/verify/ch14_thermal_series.csv --rising` (the verification run with the thermal series; `python lim_simulation.py --rising-radius` is the same run without it); `python tools/rising_radius_hanging.py <daily.csv> --depths=100,377,550`; `python tools/rising_radius_power_sweep.py --run --jobs=3` then `--table`; `ORBITAL_RING_OUT=runs/rising_quadratic/all python run_all.py --no-lim`; `python tools/generate_equations_latex.py`.
- The October results: `ORBITAL_RING_QUADRATIC=fixed` before any script, including `run_all.py` and the two LIM runs; `python -m pytest tests` runs under it by itself.

## 14. Round 5 (2026-11-07): the number of turns at 12 MW

Asked for in `Vol_III/Claude_Code_instructions_turns_12MW_2026-11-06.md`. Paul's ruling after round 4: the design is two LIMs per site at 12 MW with the cable of the rising quadratic. The 208 turns per phase coil were chosen at 8 MW with the October cable (Section 9.7 of the book: 100, 150, 208 and 300 turns run on the deployment simulation, the fastest kept, hysteresis the second check); this round redoes that choice at 12 MW with the new cable. `tools/rising_radius_turns_sweep.py --run` ran the cases at the full time steps (1 s on the first day, then 10 s; 7 to 30 minutes each, four at a time), `--table` made the table, the curve and the figure: `runs/rising_radius_turns/rising_radius_turns.md`, `.csv`, `.svg|png`, one folder per case. The 208-turn case is the 12 MW run of the round-4 sweep copied in (437.0 days).

### 14.1 What `lim_config.configure()` does with N

For every N the measured tape basis of 2026-10-27 stands: the flat winding of two 3 mm tape layers (`HTS_TAPE_LAYERS_MEASURED = 2`, `TAPE_LIMIT_MODE = measured`), the worst tape's critical current 183.1 A, and the coil current cap `I_TARGET = I_PEAK = I_COIL_LIMIT_MEASURED = 183.1 A`, the current at which that tape is at half its limit. `configure(N_TURNS=N)` recomputes the winding width ((3 + 1.06) mm x N: 0.49 m at 120 turns, 0.61 at 150, 0.84 at 208, 1.14 at 280), the tape length per LIM (`L_HTS_LIM` = 2 x (w_coil + tau_p) x N x 3 phases x 3 pitches, per layer), the coil cryostat surface, and the hysteresis loss table (`lim_hysteresis.loss_table` at 1.05 x the cap) for that N; the layer count and the cap it does not change, since those are measured figures for the 208-turn winding. The cap is therefore a choice at the other N: the worst tape of a flat winding sits in the field of the other turns, which grows with the winding width, so `lim_hysteresis.coil_limit("flat", N, 2, 3 mm)` gives each winding its own cap: 185.3 A at 120 turns, 184.8 at 135, 184.3 at 150, 184.0 at 165, 183.7 at 180, 183.1 at 208, 182.6 at 240, 182.1 at 280 (the worst tape's limit 182.1 to 185.3 A at 0.068 to 0.074 T and 7.3 to 8.1 degrees). The brief's six N were run both ways: with the 183.1 A cap as `configure()` leaves it, and with the winding's own cap (`configure(I_COIL_LIMIT_MEASURED=...)`, the `_own` cases). The two differ by 2 days or less at every N but 120 (16 days, the one case that is current-capped to the end), so the choice of cap does not move the choice of N. The grid was extended with 135 and 165 turns (the 183.1 A cap) to place the minimum, and the chosen N was run at 8 and 16 MW for Section 9.8.

### 14.2 The runs

The columns: the day the current cap stops binding; the days the 12 MW power cap binds (first and last); the day the 100 kV coil limit arrives (the definitions of Section 10.6); the peak thrust per LIM and its day; the LIM work; the peak hysteresis per LIM; the peak and end cryogenic power per site; the peak plate temperature; the HTS tape per site (two LIMs, two layers of 3 mm tape); the coil's inductance at the end of the run as the model implies it, L = V / (omega I) from the induced voltage omega N B A (the model has no separate coil inductance; L grows as N^2, 2.83 H at 150 turns against 5.44 H at 208); the voltage, current and thrust per LIM at the end.

| Turns | Cap A | Days (months) | I cap to day | P cap, days | 100 kV from day | Thrust max N/LIM (day) | LIM work EJ | Hyst peak kW/LIM | Cryo peak / end MW | Plate K | Tape t/site (km) | Coil width m | L at the end H | V, I, thrust at the end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 120 | 183.1 | 650.1 (21.4) | 650.1 (the end) | never | never | 649 (117.5) | 19.175 | 14.18 | 3.50 / 3.50 | 304.5 | 1.88 (441) | 0.487 | 1.81 | 93.9 kV, 183.1 A, 269 N |
| 120 | 185.3 (own) | 633.9 (20.8) | 633.9 | never | never | 664 (115.0) | 19.175 | 14.73 | 3.59 / 3.59 | 305.8 | 1.88 (441) | 0.487 | 1.81 | 95.1 kV, 185.3 A, 276 N |
| 135 | 183.1 | 448.4 (14.7) | 278.1 | 277 to 448 (the end) | never | 821 (95.3) | 19.175 | 10.76 | 2.96 / 2.20 | 317.8 | 2.12 (496) | 0.548 | 2.29 | 88.2 kV, 137.9 A, 459 N |
| **150** | **183.1** | **381.1 (12.5)** | 132.9 | 130 to 381 (the end) | never | 1,013 (79.6) | 19.175 | 8.15 | 2.57 / 2.03 | 331.2 | 2.35 (551) | 0.609 | 2.83 | 98.9 kV, 125.3 A, 468 N |
| 150 | 184.3 (own) | 379.1 (12.5) | 129.3 | 126 to 379 | never | 1,027 (78.7) | 19.175 | 8.27 | 2.59 / 2.03 | 332.1 | 2.35 (551) | 0.609 | 2.83 | 98.9 kV, 125.3 A, 468 N |
| 165 | 183.1 | 401.8 (13.2) | 106.3 | 88 to 280 | 270.3 | 1,226 (85.6) | 19.175 | 8.05 | 2.59 / 1.72 | 344.5 | 2.59 (606) | 0.670 | 3.42 | 100 kV, 104.1 A, 275 N |
| 180 | 183.1 | 400.8 (13.2) | 81.9 | 66 to 196 | 187.9 | 1,459 (31.6) | 19.175 | 7.83 | 2.59 / 1.56 | 357.8 | 2.82 (661) | 0.731 | 4.07 | 100 kV, 87.9 A, 331 N |
| 180 | 183.7 (own) | 399.3 (13.1) | 81.1 | 65 to 196 | 187.6 | 1,468 (31.6) | 19.175 | 7.87 | 2.60 / 1.56 | 358.3 | 2.82 (661) | 0.731 | 4.07 | 100 kV, 87.9 A, 331 N |
| 208 | 183.1 | 437.0 (14.4) | 44.0 | 40 to 95 | 89.2 | 1,949 (36.1) | 19.175 | 7.29 | 2.57 / 1.46 | 382.2 | 3.26 (764) | 0.844 | 5.44 | 100 kV, 65.8 A, 247 N |
| 240 | 183.1 | 541.1 (17.8) | 25.9 | 25 to 44 | 43.5 | 2,535 (25.7) | 19.175 | 6.18 | 2.46 / 1.48 | 409.4 | 3.76 (881) | 0.974 | 7.24 | 100 kV, 49.5 A, 185 N |
| 240 | 182.6 (own) | 541.2 (17.8) | 26.1 | 26 to 45 | 43.6 | 2,525 (25.9) | 19.175 | 6.16 | 2.46 / 1.48 | 408.9 | 3.76 (881) | 0.974 | 7.24 | 100 kV, 49.5 A, 185 N |
| 280 | 183.1 | 719.4 (23.6) | 17.4 | 17 to 29 | 28.2 | 3,033 (17.3) | 19.175 | 4.91 | 2.35 / 1.55 | 442.2 | 4.39 (1,028) | 1.137 | 9.85 | 100 kV, 36.3 A, 136 N |
| 280 | 182.1 (own) | 719.4 (23.6) | 17.6 | 17 to 29 | 28.4 | 3,017 (17.6) | 19.175 | 4.86 | 2.34 / 1.55 | 440.9 | 4.39 (1,028) | 1.137 | 9.85 | 100 kV, 36.3 A, 136 N |

The LIM work is 19.175 EJ in every case, as it must be with the same cable and no anchored phase; every run ends on the design tension at 700 km with the casing at ground speed.

### 14.3 The curve, the minimum and its width

Days against turns with the 183.1 A cap: 650, 448, 381, 402, 401, 437, 541, 719 days at 120, 135, 150, 165, 180, 208, 240, 280 turns; with each winding's own cap 634, 379, 399, 437, 541, 719 at 120, 150, 180, 208, 240, 280. **The minimum is at 150 turns, 381.1 days (12.5 months)**, 56 days (13 %) faster than the 208-turn design at the same 12 MW. The minimum is narrow: no other N of the grid is within 2 % of it (165 and 180 turns are 5 % above, 135 turns 18 %, 208 turns 15 %), so the band of N within 2 % is about 145 to 155 turns by the slopes on either side (4.5 days per turn below 150, 1.4 days per turn above). The figure `rising_radius_turns.png` shows the two curves.

Why 150 wins. The thrust per LIM at a given current grows with N (the field is N x I), so at the current cap more turns push harder: 649 N at 120 turns, 1,013 at 150, 1,949 at 208, 3,033 at 280 at the peak. But the voltage grows as N^2 at a given current (L = k N^2) and as N at a given field, and the 100 kV limit arrives sooner the more turns there are: day 270 at 165 turns, 188 at 180, 89 at 208, 44 at 240, 28 at 280, and never at 150 or fewer. Once the voltage limit holds, the current must fall as the relative speed and the supply frequency rise, and the thrust with it: at the end 136 N per LIM at 280 turns, 185 at 240, 247 at 208, 331 at 180, 275 at 165. The 150-turn coil is the one that uses all of the 12 MW to the end without reaching 100 kV: it sits at the power cap from day 130 to day 381 at 98.9 kV with 125 A and 468 N per LIM at the end, twice the 208-turn coil's 247 N. Below 150 turns the current cap binds instead: 135 turns reaches 12 MW only on day 277 (821 N peak) and 120 turns never (649 N, 183.1 A to the end, 94 kV).

The second check of Section 9.7, hysteresis and tape: the 150-turn coil's peak hysteresis is 8.15 kW per LIM against 7.29 at 208 turns (12 % more, because it runs more current for longer; the hysteresis per LIM falls with N at the cap, 14.2 kW at 120 turns to 4.9 at 280), the peak cryogenic power per site is the same 2.57 MW, and at the end 2.03 MW against 1.46 (the 150-turn coil is still at full power then). The tape per site is 2.35 t and 551 km against 3.26 t and 764 km, 28 % less, and the coil cryostat surface falls with the winding width (0.61 m against 0.84). The plate's first-day peak is 331 K against 382 K. So 150 turns is faster, lighter in tape and cooler, and costs 0.9 kW of hysteresis per LIM at the peak and 0.6 MW of cryogenic power per site at the end of the run. The one case where the cap choice matters, 120 turns, is 2 % faster with its own 185.3 A cap and is out of the running either way.

### 14.4 The chosen N at 8, 12 and 16 MW, for Section 9.8

| 150 turns, two LIMs | 8 MW | 12 MW | 16 MW |
|---|---|---|---|
| Deployment | 534.5 days (17.6 months) | 381.1 days (12.5) | 520.3 days (17.1) |
| Current cap to day; power cap, days; 100 kV from day | 113.6; 105 to 534 (the end); never | 132.9; 130 to 381 (the end); never | 228.8; never; 225.5 |
| Thrust max N/LIM; at the end | 1,013; 300 | 1,013; 468 | 1,013; 201 |
| Hysteresis peak kW/LIM; cryogenic peak / end MW per site | 6.15; 2.23 / 1.58 | 8.15; 2.57 / 2.03 | 12.13; 3.23 / 2.03 |
| Voltage, current, slip at the end | 79.4 kV, 100.5 A, (power-capped) | 98.9 kV, 125.3 A, 0.49 % | 100 kV, 124.9 A, 1.88 % |
| Against 208 turns at the same power (Section 13.3) | 510.2 days | 437.0 days | 438.3 days |

At 8 MW the 150-turn coil is 24 days slower than 208 turns (the book's Section 9.7 finding at 8 MW with the October cable was the same way round, 26 days); at 12 MW it is 56 days faster; at 16 MW it is 82 days slower than 208 turns and 139 days slower than itself at 12 MW. The 16 MW run is the same controller effect as the 150-turn 32 MW case of Section 13.3: with the higher cap the coil draws more current early, reaches 100 kV on day 225 while its slip is still high, and from there runs voltage-limited at 1.9 % slip with 201 N per LIM, whereas the 12 MW run, power-capped, has its slip brought down to 0.5 % and pushes 468 N at the same 125 A and 99 kV. In `lim_simulation.py` the current is scaled down at the voltage limit and at the power cap alike (by V_max / V and by sqrt(P_max / P)); the slip is set by its own rule, and the end states show that rule treating the two limits differently, 0.49 % slip when power-capped and 1.88 % when voltage-limited at the same 125 A, so a site that meets the voltage limit early keeps a high slip and a low thrust per ampere for the rest of the run. Which branch of the slip rule does this was not traced in this round. A 16 MW site could always run the 12 MW schedule, so the 16 MW figure is pessimistic by the controller, not by the physics; the 12 MW figure stands on its own, since that run never meets the voltage limit. The controller's slip rule at the voltage limit is worth a look before Section 9.8 quotes the 16 MW column; it also means that the "12 MW knee" of Section 13.3 for 208 turns (437 days, with the rise of 0.3 to 0.5 days per MW above it) carries the same caveat above 12 MW.

### 14.5 Files

Under `runs/rising_radius_turns/` (copies in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/rising_radius_turns/`): `rising_radius_turns.md`, `.csv`, `.svg|png`, `_run_console.txt`, `_run_console2.txt`, and one folder per case (`N120` to `N280`, the `_own` cases, `N150_P8MW`, `N150_P16MW`) with the log, daily table, Figure 3.2 series and summary. Reproduce: `python tools/rising_radius_turns_sweep.py --run --jobs=4` then `--table`.

## 15. Round 6a (2026-11-07): the slip rule at the voltage limit, and the reruns Section 9.8 needs

Asked for in `Vol_III/Claude_Code_instructions_slip_rule_and_sheet_2026-11-07.md`. Paul's design after round 5: two LIMs per site, 12 MW, 150 turns, the cable of the rising quadratic, 381.1 days. Section 14.4 noted that a coil held at the 100 kV limit kept about 1.9 % slip where the same coil power-capped ran at 0.5 %, so every run that met the limit early was slower than the physics requires.

### 15.1 The slip rule, traced

`lim_simulation.py` sets the slip target from the site's power margin, 1 - P / P_max: at a margin of 15 % or more the target is `SLIP_RATIO_NORMAL` (2 %), at 5 % or less `SLIP_RATIO_REDUCED` (0.5 %), between them a blend; a current under 70 % of the cap also pulls the target toward the reduced value. The slip command then follows the target with a 15-day time constant. The three limits act on the current only: the voltage limit scales it by V_max / V, the power cap by sqrt(P_max / P), the plate temperature by 0.95 per iteration. So a coil that meets the voltage limit has its current cut until V = V_max, its site power falls well under the cap (12.4 MW of 32 at the end of the 150-turn 32 MW run of Section 13.3, 6.5 of 16 in the 150-turn 16 MW run of 14.4), the power margin is wide, and the rule holds the slip at 2 %: the voltage-limited coil sits at the high slip the rule reserves for a site with power to spare. At the same current and voltage the 0.5 % slip gives twice the thrust (468 N against 201 N at 125 A and 99 to 100 kV in the two 150-turn runs of 14.4), since in the plate's eddy-current model the thrust per ampere falls with the slip over this range, and the slip adds only its own small share to the supply frequency, so the voltage at a given current barely moves with it.

What the physics allows: the induced voltage is omega N B A with omega from the supply frequency, f = (v_rel + v_slip) / (2 tau_p); lowering the slip lowers the frequency a little (v_slip is 0.5 to 2 % of v_rel) and raises the thrust per ampere, so at the voltage limit the controller can hold V = V_max with the reduced slip and the thrust of the power-capped case. The model only lacked the rule.

The change (`lim_config.SLIP_RULE`, environment `ORBITAL_RING_SLIP_RULE`): with `both`, the default since this round, the margin that sets the slip target is the smaller of the power margin and the voltage margin, 1 - V / V_max, so a coil at the voltage limit goes to the reduced slip as a site at its power cap does; `power` is the October rule, the power cap alone. Nothing else in the controller changed. Under `power` the code's expressions are the October ones, and the fixed-radius quick run (log and tension series) is byte-identical to the run of the code before the change; `tests/conftest.py` sets `power` for the reference checks, with `fixed` for the quadratic, and the 29 tests pass. The design run, 150 turns at 12 MW, is unchanged by the rule to the digit (381.1 days: it never meets the voltage limit, and its slip was already at the reduced value from the power cap); the 150-turn 16 MW run goes from 520.3 days to 354.9 and ends at 0.5 % slip with 477 N per LIM (quick steps).

### 15.2 The reruns Section 9.8 and Table 9.3 quote

`tools/rising_radius_design_power.py --run` (eleven full-step runs with the corrected rule, four at a time, 12 to 28 minutes each; `--table` for the table, the curves and the figure `runs/rising_radius_design/rising_radius_design_power.md|csv|svg|png`, one folder per case). The columns are those of Section 14.2, with the slip at the end of the run added. Every run ends on the design tension at 700 km with the casing at ground speed and no anchored phase.

| LIMs | Turns | MW/site | Days (months) | I cap to day | P cap, days | 100 kV from day | Thrust max N/LIM (day) | LIM work EJ | Hyst peak kW/LIM | Cryo peak / end MW | Plate K | Tape t/site | V, I, thrust, slip at the end | Cable kg/m, Mt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 150 | 8 | 534.5 (17.6) | 113.6 | 105 to the end | never | 1,013 (29.6) | 19.175 | 6.15 | 2.23 / 1.58 | 331.2 | 2.35 | 79.4 kV, 100.5 A, 300 N, 0.50 % | 103,570; 4,606 |
| 2 | 150 | **12** | **381.1 (12.5)** | 132.9 | 130 to the end | never | 1,013 (79.6) | 19.175 | 8.15 | 2.57 / 2.03 | 331.2 | 2.35 | 98.9 kV, 125.3 A, 468 N, 0.49 % | 103,570; 4,606 |
| 2 | 150 | 16 | 354.9 (11.7) | 210.2 | 209 to 256 | 249.6 | 1,013 (79.6) | 19.175 | 11.56 | 3.14 / 2.06 | 331.2 | 2.35 | 100 kV, 126.6 A, 477 N, 0.50 % | 103,570; 4,606 |
| 2 | 150 | 20 | 352.6 (11.6) | 222.2 | never | 219.9 | 1,013 (79.6) | 19.175 | 12.13 | 3.23 / 2.06 | 331.2 | 2.35 | 100 kV, 126.6 A, 477 N, 0.50 % | 103,570; 4,606 |
| 2 | 150 | 24 | 352.6 (11.6) | 222.2 | never | 219.9 | 1,013 (79.6) | 19.175 | 12.13 | 3.23 / 2.06 | 331.2 | 2.35 | 100 kV, 126.6 A, 477 N, 0.50 % | 103,570; 4,606 |
| 2 | 208 | 12 | 436.0 (14.3) | 44.0 | 40 to 95 | 89.2 | 1,949 (36.1) | 19.175 | 7.29 | 2.57 / 1.46 | 382.2 | 3.26 | 100 kV, 65.8 A, 247 N, 0.50 % | 103,570; 4,606 |
| 2 | 208 | 16 | 428.1 (14.1) | 59.2 | 52 to 60 | 59.1 | 1,949 (49.0) | 19.175 | 8.61 | 2.79 / 1.46 | 382.2 | 3.26 | 100 kV, 65.8 A, 247 N, 0.50 % | 103,570; 4,606 |
| 4 | 150 | 16 | 255.2 (8.4) | 50.7 | 46 to the end | never | 1,013 (42.9) | 19.148 | 6.09 | 3.61 / 2.42 | 385.8 | 4.70 | 81.6 kV, 103.5 A, 320 N, 0.49 % | 105,010; 4,670 |
| 4 | 150 | 24 | 187.8 (6.2) | 71.5 | 70 to the end | 184.5 | 1,013 (38.6) | 19.148 | 8.10 | 4.29 / 3.30 | 385.8 | 4.70 | 100 kV, 126.8 A, 482 N, 0.50 % | 105,010; 4,670 |
| 4 | 208 | 16 | 247.0 (8.1) | 18.2 | 18 to 154 | 146.5 | 1,705 (18.2) | 19.148 | 4.25 | 3.29 / 2.09 | 453.5 | 6.52 | 100 kV, 65.9 A, 249 N, 0.50 % | 105,010; 4,670 |
| 4 | 208 | 24 | 220.1 (7.2) | 24.7 | 24 to 44 | 42.5 | 1,856 (24.6) | 19.148 | 6.21 | 3.94 / 2.09 | 453.5 | 6.52 | 100 kV, 65.9 A, 249 N, 0.50 % | 105,010; 4,670 |

What the corrected rule changed, against the runs of Sections 13.3 and 14 (the October rule): every voltage-limited run now ends at 0.5 % slip with the thrust of its power-capped equivalent, 477 N per LIM at 150 turns and 247 to 249 N at 208 (the 208-turn coil's thrust at 100 kV was already its reduced-slip value in the old runs, since those met the limit late enough for the current-ratio rule to have pulled the slip down). The design, 150 turns at 12 MW, is 381.1 days as before. 150 turns at 16 MW: 354.9 days against 520.3. 208 turns at 12 MW: 436.0 against 437.0; at 16 MW: 428.1 against 438.3. The four-LIM cases do not move (255.2, 187.8, 247.0, 220.1 against 255.2, 187.8, 247.0, 220.2): they never met the limit, or met it with the slip already reduced.

**The design's power curve, 150 turns and two LIMs.** 534.5, 381.1, 354.9, 352.6, 352.6 days at 8, 12, 16, 20, 24 MW: 38.3 days per MW from 8 to 12, 6.6 from 12 to 16, 0.6 from 16 to 20, none beyond. The floor is 352.6 days and the knee by the 2 % rule is **16 MW, 354.9 days (11.7 months)**; 12 MW is 8 % above the floor. The coil meets 100 kV on day 250 at 16 MW, on day 220 at 20 and 24 MW, and from there the site never draws its cap (the 20 and 24 MW runs are the same run: the cap never binds). Against 12 MW, 16 MW buys 26 days for a 110 m panel in place of 83 m (1,423 GW against 1,067), 0.57 MW more cryogenic plant at the peak (3.14 against 2.57 MW) and 3.4 kW more hysteresis per LIM at the peak (11.56 against 8.15), with the same tape, plate temperature and end-of-run cryogenic power. Paul's 12 MW stands as the design; the curve is what Section 9.8 can print.

**The 208-turn coil, for 9.7 and 9.8.** At 12 MW the 150-turn coil deploys in 381.1 days against 436.0 for 208 turns (55 days, 13 %), at 16 MW in 354.9 against 428.1 (73 days). The 208-turn coil's own power curve is flat from 12 MW up (436, 428 days) because it meets 100 kV on day 89 at 12 MW and day 59 at 16 MW and runs voltage-limited at 247 N per LIM from there; the 150-turn coil holds 12 MW to the end without meeting the limit and ends with 468 N.

**The fast options, four LIMs per site.** 150 turns: 255.2 days at 16 MW, 187.8 at 24 MW; 208 turns: 247.0 and 220.1 days. The 150-turn four-LIM site at 24 MW is the shortest deployment in the set, 6.2 months, with 4.70 t of tape per site against 6.52 for 208 turns, a 386 K plate against 454 K, and 4.29 MW of cryogenic plant at the peak against 3.94.

### 15.3 The thermal series of the design, for Figures 14.1, 14.2 and Table 14.2

`tools/ch14_thermal_series.py runs/rising_radius_design/thermal/ch14_thermal_series_design.csv --rising --set N_TURNS=150 --set MAX_SITE_POWER=12e6` (the tool takes `lim_config.configure` overrides since this round; the same run as the design's, full steps, 164,600 rows at 200 s), the figures by `tools/fig14_figures.py` in `runs/rising_radius_design/thermal/figures/` and installed as the book's Figures 14.1 and 14.2 (the originals in `Vol_III/images/figures_backup_2026-11-06b/`).

The history is different in shape from the 208-turn, 8 MW one the book describes. The 150-turn coil runs at its 183 A cap for 133 days and then at the 12 MW cap to the end, so its losses do not fall away after the first weeks: the cable peaks at 331 K on day 1 (not 382), falls to 278 K at 1.5 months, rises again to 313 K at four months while the coil is at its current cap, and from month 5 holds 270 K to the end (not 257 K), the wall 7 K below it; the end-of-run eddy loss is 20.6 kW per LIM at 468 N and 0.5 % slip. The cold-side load starts at 14.5 kW, is 17.5 kW at one month, 25.0 at three, peaks at 30.8 kW at 5.0 months (the hysteresis peak, 16.3 kW for the two LIMs, is at the same time, not at seven weeks), is 26.3 kW at eight months and 24.4 kW at the end (12.7 months), with the radiator 4.99 m at the peak and 3.95 m at the end.

Table 14.2 of the book against the design (the compressor power by the book's rule, the cold-side load over the COP of 0.0106; the simulation's own `p_cryo` in brackets):

| Row | Book (208 turns, 8 MW, October cable): peak at seven weeks / steady / after | Design: peak at 5.0 months / steady state, months 6 to 10 / after deployment |
|---|---|---|
| Hysteresis in the two LIMs | 10.2 / 1.8 / 0 kW | 16.3 / 12.0 / 0 kW |
| Leak through the cryostat insulation | 6.2 kW | 4.5 kW (the 150-turn coil's cryostat is 2,236 m^2 per site against 3,101) |
| Heat onto the levitation coils | 10.0 kW | 10.0 kW |
| Cold-side load | 26.4 / 18.0 / 16.2 kW | 30.8 / 26.5 / 14.5 kW |
| Compressor power | 2.49 / 1.70 / 1.53 MW | 2.90 (2.57) / 2.50 / 1.37 (1.21) MW |
| Heat rejected at 400 K | 2.52 / 1.72 / 1.55 MW | 2.93 / 2.53 / 1.38 MW |
| Cryogenic radiator area | 1,925 / 1,314 / 1,183 m^2 | 2,247 / 1,935 / 1,058 m^2 |
| Radiator width over 450 m | 4.28 / 2.92 / 2.63 m | 4.99 / 4.30 / 2.35 m |

Every row changes: the peak is 17 % higher and five months in, the steady state 47 % higher (the coil is at full power to the end), and the standing load after the deployment 10 % lower, 14.5 kW and 1.37 MW per site, 122 GW around the ring (Chapter 8's largest standing load, 136 GW in the book). The steady-state column is the mean of months 6 to 10, the plateau of the power-capped run; the book's "months 12 to 15" window no longer exists in a 12.5-month deployment. Table 14.3 (the warm loop) moves with the cable's temperature, 270 K in the steady state against 257.

### 15.4 The numbers sheet

Part 4 of the brief, `docs/design_point_numbers_2026-11.md` (copied to `Vol_III/`): one table per chapter, where in the book, the October value, the new value, what the chapter carries today, and the script. The hard-coded October cable found and fixed on the way: `debris_simulation.ring_state` (its own copy of the fixed quadratic, so the Chapter 4 orbits, lifetimes and soot had not followed the design point), `radiation/sweetspot.py` and `radiation/fig05_04_sky_under_cable.py` (the 8.20 and 8.22 m cable), `tools/ch03_direction_numbers.py` (the cable's direction from `DELTA_V_LOAD`), `tools/ch08_numbers.py` (8,786.2 m/s typed in; a 12 MW case added). The figures that carry the cable were redrawn from the new runs in their current layouts (1.1, 3.2, 4.1 to 4.4, 4.6, 4.7, 5.4, 7.4, 14.1, 14.2; originals in `images/figures_backup_2026-11-06b/`); the chapter files were not edited.

### 15.5 Files

Under `runs/` (copies in `Vol_III/runs_700km_hw10686_2026-10-12/rising_radius_2026-11-06/`): `rising_radius_design/` (the eleven reruns, the table and figure, `thermal/` with the design's series and figures, `prograde/` with the prograde design run for Figure 3.2's dashed curve: 311.8 days, 14.902 EJ, no anchored phase); `rising_quadratic/all/` refreshed (the Chapter 4 numbers and figures, `debris_orbits`, `altitude_trade`, `debris_simulation`, `casing_6000`, `ch08_numbers`, `ch03_direction_numbers`, the radiation runs `radiation/` and Figure 5.4, the bearing figure); `design_point_sheet/` (the tools run under both quadratics). Reproduce: `python tools/rising_radius_design_power.py --run --jobs=4` then `--table`; the thermal series as above; `ORBITAL_RING_SLIP_RULE=power` for the October controller.
