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
- `rising_radius_deployment_tension.csv`: the Figure 3.2 series (the same seven columns as `lim_deployment_tension.csv`, plus altitude, radius, strain, phase and anchor force), 4,052 rows.
- `rising_radius_summary.json`: the headline numbers.
- `rising_radius_hanging.txt`, `rising_radius_hanging.csv`: Section 7.
- `figures_rising_radius/rising_radius_altitude.svg|png`: the ring's altitude, the point 100 km below, the longest hang above 150 km.
- `figures_rising_radius/rising_radius_tension.svg|png`: tension and stress against time, rising against fixed.
- `figures_rising_radius/rising_radius_speeds.svg|png`: cable and casing speeds, both runs.
- `figures_rising_radius/rising_radius_anchor_phase2.svg|png`: the anchor reaction and the LIM thrust through phase 2.
- `figures_rising_radius/rising_radius_energy.svg|png`: the energy budget accumulating.
- `figures_rising_radius/fig03_02_deployment_tension_700_rising.svg|png`: Figure 3.2 redrawn from the rising series (`tools/fig03_02_deployment_tension.py --rising`). Its top and middle panels are against the casing's speed, which is constant through phase 2, so the last 163 GN appear as a vertical step at 8,114 m/s; a Figure 3.2 for the rising radius wants the cable's speed or time on that axis, which is a drawing decision for the book.

Reproduce: `python lim_simulation.py --rising-radius > runs/rising_radius_lim_simulation.txt` (about 3 hours when run beside two other deployment runs), then `python tools/rising_radius_hanging.py`, `python tools/fig_rising_radius.py`, `python tools/fig03_02_deployment_tension.py --rising --out=runs/figures_rising_radius`.
