#!/usr/bin/env python3
"""
Section 3b of round 8's numbers sheet: the rows the Chapter 7 rewrite needs (thirteen items), as markdown lines, from runs/round8_logs/round8_rows.json
(bearing_round8.py --part=rows) and round 7's logs.  Called by tools/round8_sheet.py.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (Sections 7.5 to 7.9, Tables 7.3, 7.4 and 7.6), Chapters 6, 8 and 16.
"""
import math

SRC = "`bearing_round8.py --part=solve-sup`, `--part=rows`"


def kA(x): return f"{x/1e3:,.1f} kA"
def kN(x): return f"{x/1e3:.1f} kN/m"
def pc(x): return f"{x*100:.0f}%"


def rows_3b(R, ctx):
    """R: round8_rows.json; ctx: dict(half, fast, J12, DP, dump, sup7, cpl, hold) from the sheet."""
    L = []; half = ctx["half"]; fast = ctx["fast"]; J12 = ctx["J12"]; DP = ctx["DP"]
    T = R["tapes_leg"]; n_t = R["n_t"]; t73 = R["t73"]; cap_NI = half["cap"]["NI"]
    L.append("\n## Section 3b: the rows for the rewrite of Sections 7.5 to 7.9\n")
    L.append(f"All of these are for the round 7 design as ruled: two support tracks, each of two independent 0.4 m loops 50 mm apart at the center on a 1.2 m by 100 mm plate; two load tracks and two lift tracks of fast coils, "
             f"each one 0.6 m loop on a 1.0 m by 200 mm plate; six plates. The support solves are at 2.5 mm cells ({SRC}; logs `round8_solve_sup_*.txt`, `round8_rows.txt`). "
             f"The support winding is built at 50 percent of the worst tape's limit at the cap: {T:,.0f} tapes of 12 mm per leg, {n_t:.1f} turns of ten soldered layers per loop and coil, ten staggered coils of 1 km per loop, {half['cap']['tape']:.1f} kg/m of tape per track.\n")

    # 1
    L.append("### 1. Table 7.3, how hard to run the worst tape of a support track\n")
    L.append(f"The worst 12 mm tape at the cap of {kA(cap_NI)} per leg sits at {t73['B']:.2f} T, {t73['th']:.0f} degrees from its normal. Its limit there is {t73['ic70']:.0f} A at 70 K and {t73['ic77']:.0f} A at 77.5 K (the two measured surfaces at that field and angle), "
             f"a fall of {(t73['ic70']-t73['ic77'])/7.5:.1f} A per kelvin, {(t73['ic70']-t73['ic77'])/7.5/t73['ic70']*100:.1f} percent of the 70 K value per kelvin, taken linear between them as in the present table. Normal operation is {t73['frac_normal']*100:.0f} percent of the cap's ampere-turns (the present table's 44 percent was on the old plate).\n")
    L.append("| Share of the 70 K limit carried at the cap | Tape current at the cap | Warming the tape can take | Tape current in normal operation | Tapes per leg | Tape per track |\n|---|---|---|---|---|---|")
    for r in t73["rows"]:
        L.append(f"| {r['share']*100:.0f}% | {r['I_cap']:.0f} A | {r['dT']:.1f} K | {r['I_normal']:.0f} A | {r['tapes']:,.0f} | {r['tape']:.0f} kg/m |")
    L.append(f"\nThe design is built on the first row and rated on the third. In normal operation the field at the worst tape is lower and its limit higher, so the tape then runs at {t73['share_normal']*100:.0f} percent of its limit. Book today: 234 A at 1.64 T and 53 degrees; 117, 140, 164, 187 A; 5.6, 4.5, 3.3, 2.2 K; 51, 61, 72, 82 A.\n")

    # 2
    L.append("### 2. Table 7.4, the pull against the ampere-turns\n")
    L.append("One split support track, both loops unless the row says otherwise. A coil carries ten times the tape current, because a turn is ten layers.\n")
    L.append("| State | Ampere-turns per leg | Current per coil | Current per tape | Share of the worst tape's limit | Pull per meter of track | Plate saturated |\n|---|---|---|---|---|---|---|")
    for r in R["t74_support"]:
        L.append(f"| {r['state'][0].upper() + r['state'][1:]} | {r['NI']:,.0f} | {r['I_coil']:,.0f} A | {r['I_tape']:.1f} A | {r['share']*100:.0f}% | {r['F']/1e3:.1f} kN/m | {r['sat']*100:.0f}% |")
    pts = R["support_points"]
    L.append(f"\nThe solved points behind the table (both loops): " + "; ".join(f"{p['NI']/1e3:.0f} kA, {p['F']/1e3:.0f} kN/m, {p['sat']*100:.0f}% saturated, {p['dev']*100:.0f}% of the square law" for p in pts) + ". "
             f"The square law is the pull of the same track on a plate that does not saturate, {R['k_lin']*1e6:.3f} N/m per kA² per leg. 'Where the pull stops growing as the square' is given as the points where the solved pull is 5 and 10 percent under it; that 5 percent is a choice, the chapter may want another. "
             f"**The 100 mm support plate is under the square law over its whole working range: {(1-pts[0]['dev'])*100:.0f} percent under in normal operation and {(1-next(q for q in pts if abs(q['k']-1.0) < 1e-6)['dev'])*100:.0f} percent at the cap.** So the chapter's 'Below saturation the pull grows with the square of the ampere-turns' holds for a fast track (1 percent under at its requirement) and not for the support track; from normal operation to the cap the support's pull grows as about the {math.log(next(q for q in pts if abs(q['k']-1.0) < 1e-6)['F']/pts[0]['F'])/math.log(next(q for q in pts if abs(q['k']-1.0) < 1e-6)['NI']/pts[0]['NI']):.1f} power of the ampere-turns. "
             f"'Plate saturated' in these tables is the share of the plate's section whose permeability has fallen below a tenth of its unsaturated value.\n")
    L.append("The same for one fast track (lift or load: the same plate, loop and bundle). A coil has 10 turns and 20 coils lie over every point, so a coil carries the leg's ampere-turns over 200.\n")
    L.append("| State | Ampere-turns per leg | Current per coil | Current per 3 mm tape | Share of the worst tape's limit | Pull per meter of track | Plate saturated | Of the square law |\n|---|---|---|---|---|---|---|---|")
    for r in R["t74_fast"]:
        L.append(f"| {r['state'][0].upper() + r['state'][1:]} | {r['NI']:,.0f} | {r['I_coil']:,.0f} A | {r['I_tape']:.1f} A | {r['share']*100:.0f}% | {r['F']/1e3:.1f} kN/m | {r['sat']*100:.1f}% | {r['dev']*100:.0f}% |")
    L.append(f"\nSource: `bearing_round7.py --part=solve-fast` (the levels at 1, 1.25, 1.5 and 2 times the requirement's ampere-turns: " + "; ".join(f"{l['k']:g}: {l['F']/1e3:.0f} kN/m, {l['sat']*100:.1f}% saturated, tape at {l['share']*100:.0f}%" for l in R["fast_levels"]) + f"). The 198 kN/m row is interpolated between the first two levels. Book today: 141,000 / 287,000 / 322,000 / 403,000 / 493,000 / 610,000 ampere-turns for 48 / 198 / 245 / 330 / 444 / 601 kN/m, with the two windings on the same loops; the last two rows have no counterpart, because nothing shares the support's loops now.\n")

    # 3
    P = R["parked"]; t210 = P["times"]["210"]; t1k = P["times"]["1000"]
    L.append("### 3. Table 7.6, the operating points of a support track\n")
    L.append("| Case | Current per coil | Ampere-turns per leg | Current per tape | Share of the worst tape's limit | Pull |\n|---|---|---|---|---|---|")
    for o in R["ops"]:
        L.append(f"| {o['case'][0].upper() + o['case'][1:]} | {o['I']:,.0f} A | {o['NI']:,.0f} | {o['I']/10:.1f} A | {o['share']*100:.0f}% | {o['F']} |")
    ch = R["change"]
    L.append(f"\n- **Normal operation.** Each loop at {P['I_n']:.0f} A per coil.")
    L.append(f"- **A parked 5,000 tonne spacecraft on 100 m of sled**, 396 kN/m, 198 kN/m to each side: each support track goes to its cap of {P['F']/1e3:.1f} kN/m, {P['I_cap']:,.0f} A per coil. The loops' flux has to rise by {P['lam']:,.0f} volt-seconds per coil (ten coils of a loop together, both loops of the track rising), which the supply gives in **{t210['T']:.0f} s at the 210 V charging voltage and {t1k['T']:.1f} s at 1 kV**. "
             f"What carries the load until then is the load track of each side. If the whole load were on from the first moment, the integral of that track's pull over its full pull would be {t210['duty']:.1f} s at 210 V and {t1k['duty']:.1f} s at 1 kV, against the fast leads' {P['budget20']:.2f} s with all 20 families connected. "
             f"**So the load tracks cannot hold a parked spacecraft through the takeover on their pulse-sized leads.** Either the load comes on no faster than the supply raises the loops (a spacecraft is set down or fueled over minutes, and the support follows it), or the load tracks' leads are sized for it; no script sizes them for that. This is a finding of this round and is for Paul to rule on.")
    L.append(f"- **Not modeled, and to be said so:** the support coils are 1 km long in ten staggered families and the spacecraft covers 100 m. Raising the coils that cross the 100 m raises the support, in steps of a tenth, over a kilometer either side of it, where there is no load. No round has run the six-plate design through that case; the earlier rounds' takeover times (Stage B, 2.3 minutes at 1 kV on the old coil) have the same gap.")
    L.append(f"- **One loop lost.** After the flux jump the loop beside it carries {ch['jump']:.3f} times its current by itself ({ch['outer']['I0'][2]:.0f} A) and its side gives {ch['F_jump']/1e3:.1f} kN/m. After the trim at 1 kV: an outer loop lost, the lost side's inner loop at {ch['outer']['I1'][2]:.0f} A and the sound side's two at {ch['outer']['I1'][0]:.0f} A, in {ch['outer']['T_hi']:.1f} s, using {ch['outer']['int_hi']:.2f} s of the fast leads' {P['budget10']:.2f} s (ten families); an inner loop lost, {ch['inner']['I1'][2]:.0f} A and {ch['inner']['I1'][0]:.0f} A, {ch['inner']['T_hi']:.1f} s and {ch['inner']['int_hi']:.2f} s (`bearing_round7.py --part=support`).")
    L.append(f"- **One loop out for a long repair.** The same state as after the trim: it is steady, and no fast coil carries anything.")
    for s in R["standing"]:
        L.append(f"- **One loop out with a standing load on that stretch, {s['case']}, {s['limit']}:** the loop gives {s['F1']/1e3:.0f} kN/m with its worst tape at {s['share']*100:.0f} percent, so the standing load may be **{s['w']/1e3:.0f} kN/m, {s['tonnes']:,.0f} tonnes on 100 m**; the sound track carries {s['F_other']/1e3:.0f} kN/m of its {P['F']/1e3:.1f}.")
    ol = R["one_loop"]
    L.append(f"\nOne loop alone was solved at " + "; ".join(f"{p['NI']/1e3:.0f} kA: {p['F']/1e3:.0f} kN/m, tape at {p['share']*100:.0f}%, plate {p['sat']*100:.0f}% saturated" for p in ol["points"]) + ". The standing-load rows balance heave and roll with the lost side's remaining loop at its own position and the sound track whole. The present chapter's rule is 349 kN/m (4,410 tonnes on 100 m) with a track out.\n")

    # 4
    c = R["capacity"]
    L.append("### 4. The capacity sentences\n")
    L.append(f"- One support track (both loops) is built to carry {P['F']/1e3:.1f} kN/m, **{c['track_over_casing']:.1f} times** the casing's 95.0 kN/m (the chapter's 2.6 stands).")
    L.append(f"- The pair can carry {c['pair']/1e3:.0f} kN/m, {c['pair_over']:.1f} times the casing's weight (the chapter's 491 kN/m and 'about five times' stand).")
    L.append(f"- The rating, the worst tape at 70 percent, is **{c['rating']/1e3:.0f} kN/m, {(c['rating_over_req']-1)*100:.0f} percent above the requirement** (the chapter: 330 kN/m, 34 percent). The 100 mm plate is {R['t74_support'][3]['sat']*100:.0f} percent saturated there.")
    L.append(f"- One loop alone gives {c['loop_cap']/1e3:.0f} kN/m at the cap's coil current, **{c['loop_over_casing']:.2f} times the casing's weight**" + (f", and {c['loop70']/1e3:.0f} kN/m at its own 70 percent point, {c['loop70_over_casing']:.2f} times." if c.get("loop70") else ".") + " So one loop cannot carry the whole casing by itself; it carries its track's half with its tape at 21 percent of the limit.")
    L.append(f"- With one loop out, the other three carry the casing with the lost side's loop at {ch['outer']['I1'][2]/P['I_n']:.2f} times its normal current.\n")

    # 5
    f = R["failures"]
    L.append("### 5. The failures the chapter quotes\n")
    L.append("| Sentence in 7.7 | Book | New | Script |\n|---|---|---|---|")
    L.append(f"| One lift track lost during a launch at 30 km/s: the casing's movement; the smallest gap | 29 mm; above 70 mm | {f['lift_lost']['exc']:.1f} mm; {f['lift_lost']['gap']:.1f} mm | `bearing_round7.py --part=cases` (`bearing_round7_casing_10ms_d20cm.txt`, case B.5.3) |")
    L.append(f"| With one lift track out, the fastest launch of 5,000 tonnes inside 20 mm | 27 km/s | {f['fastest']} km/s (19.0 mm at 27 km/s, 23.7 at 28) | the same, the speed sweep |")
    L.append(f"| With one load track out, the heaviest spacecraft that can be parked and launched | about 2,480 tonnes | parked: not limited by a load track any more, because a standing load is carried by the support tracks at their cap. Launched: the full 5,000 tonnes passes 1 km/s on one load track with the casing moving {f['one_load_1kms']['exc']:.1f} mm and {f['one_load_1kms']['gap']:.1f} mm of gap left; inside one track's full pull of 198 kN/m (the pair's {2*f['pull_1kms_track']/1e3:.0f} kN/m at 1 km/s scaled by the mass) it is {f['mass_by_pull']:,.0f} tonnes | the same, F.1's one-track lane; the mass by arithmetic on `bearing_round4.pull_track_of_u` |")
    L.append(f"| The number of cases and the smallest gap among them | 114 cases; 50 mm, a pull 50 ms late | {f['n_cases']} cases; {f['gap_single']:.1f} mm, a pull 50 ms late | the same, B.6 |")
    L.append(f"\nA sentence the book can quote: \"The design was run through {f['n_cases']} cases in simulation. No single failure leads to contact, and the smallest gap in any single failure is {f['gap_single']:.0f} mm, when a pull fires 50 milliseconds late.\" "
             f"The cases are round 4's, rerun in round 7 with the 10 ms rise, the support frozen, the 10 mm protection and the load track's coupling into the split support (11.7 percent). In them a lost support track is half the support gone in heave, which bounds a lost loop; the roll of a lost loop is round 7's separate 2-DOF run (97.2 mm of gap left, 0.38 mrad). A sled passing during a change-over is not a case.\n")

    # 6
    L.append("### 6. The stray field\n")
    L.append("| State | 1 m below a support track | 1 m below a load track | At a LIM stator's center | Largest over the stator |\n|---|---|---|---|---|")
    for k, v in R["stray"].items():
        L.append(f"| {k[0].upper() + k[1:]} | {v['1 m below a support track']*1e3:.1f} mT | {v['1 m below a load track']*1e3:.1f} mT | {v['LIM stator, center']*1e3:.2f} mT | {v['largest over the stator']*1e3:.2f} mT |")
    L.append(f"\nAll six tracks in one 30 by 30 m box at 50 mm cells with linear plates, the method of `final_bearing.py` D.2 on the round 7 section (`bearing_round8_parts.stray`). The LIM stator is taken where `final_bearing.py` puts it, 0.5 m outside the cable's side and centered 1 m above the cable's centerline [UNVERIFIED position; the chapter says only 'offset above the centerline']. The plates saturate at the cap and above, which a linear plate does not carry, so the cap and rating rows are low by an amount this model cannot give. Book today: 15 mT one meter below a track and 0.3 mT at the LIM stators in normal operation, 13 mT at the stators at a track's rating (that rating had the load winding and the lift track on).\n")

    # 7
    L.append("### 7. The coil ends\n")
    fc = R["face"]
    L.append("| Coils | Heat, per meter of ring | Around the ring (heat, and the drag's power) |\n|---|---|---|")
    for e in R["ends"]:
        L.append(f"| {e['L']/1e3:.0f} km, ten per bundle, {e['lab']} | {e['P']:.2f} W/m | {e['MW']:,.0f} MW |")
    L.append(f"\nThe design's row is the first: 1 km coils, ten staggered families, two loops per track, **{R['ends'][0]['P']:.1f} W per meter of ring and {R['ends'][0]['MW']:.0f} MW around the ring**. "
             f"The model is `final_bearing.py` C.3's with the field put in summed across the plate, as round 7 did for the plate heat (`bearing_round8.py --part=solve-face`, `--part=rows`): a coil end takes a tenth of its own loop's field away (ten staggered coils) over 300 mm, twice the distance from the leg to the plate [UNVERIFIED, as before; a coil end in three dimensions is beyond the 2-D solver]; "
             f"the step's eddy heat is half a rise-and-fall cycle of the 1-D diffusion into the 100 mm plate, {fc['E1']:.1f} J/m² for a step of 1 T, as the step squared; it is summed across the 1.2 m with the face field of one loop alone at the normal current "
             f"(peak {fc['one_peak']:.2f} T, root mean square {fc['one_rms']:.2f} T, 5 mm cells); there are 2 ends times 2 loops times 10 coils per kilometer of track, and two tracks. The two loops of a track are separate coils now, which does not change the count. "
             f"Putting the plate's single peak-cell field of {R['B_plate_normal']:.2f} T into the same formula would give {R['ends'][0]['P_pk']:.0f} W/m, eight times more; that is the overstatement round 7 found in the plate heat and is not the figure to quote. "
             f"Book today: 0.1 to 4 W/m and 5 to 183 MW 'depending on the coil length' (4 W/m and 183 MW were the 1 km coils on the old plate).\n")

    # 8
    m = R["mass"]
    L.append("### 8. The bearing's mass on the casing\n")
    L.append("| Item | Per meter of casing | From |\n|---|---|---|")
    L.append(f"| Tape, two support tracks | {m['tape_sup']:.0f} kg/m | {half['cap']['tape']:.1f} kg/m per track, the 2.5 mm solve |")
    L.append(f"| Tape, four fast tracks | {m['tape_fast']:.0f} kg/m | {fast['tape']:.1f} kg/m per track, `bearing_round7.py --part=fast` |")
    L.append(f"| Formers, cryostats and supports | {m['formers']:.0f} kg/m | taken at one tape mass, the rule of Part 3b that the chapter uses; no design gives it, and the four support cryostats in place of two are not weighed |")
    L.append(f"| The stores, film capacitors | {m['store']:.0f} kg/m ({m['store_shared']:.0f} if a site's lift and load tracks share them) | {fast['E_coil']/1e6:.3f} MJ per coil, 3.2 coils per meter, 1.5 kJ/kg (round 5's figure for film capacitors) |")
    L.append(f"| The switches | {m['switches']:.0f} kg/m | two switches and two diodes per coil, 5 modules of 4.3 kg in series per position at {fast['V']/1e3:.1f} kV (round 5's rule) |")
    L.append(f"| The fast coils' leads | {m['leads_fast']:.1f} kg/m | two leads of 21 mm² by 1 m of copper per coil |")
    L.append(f"| **Total** | **{m['total']:,.0f} kg/m** ({m['total_shared']:,.0f} with the stores shared) | {m['total']/12000*100:.0f} percent of the casing's 12,000 kg/m |")
    L.append(f"\nIt is inside the 12,000 kg/m in the sense that it is smaller than it; Table 16.1's list of weighed items is what the claim has to be set against, and the stores and switches are new entries there (the chapter's 'about 1,030 kg' was tape and formers only). Not weighed in any script: the support coils' strips, diodes, dump resistors and retractable leads, the supplies, the sensors, the capacitor housings and the busbars.\n")

    # 9
    hd = ctx["hold"]
    L.append("### 9. Two terms, in a sentence each\n")
    L.append(f"- **\"The 10 mm protection in every lane.\"** Each side of the cable has its own gap sensors, and if the measured gap on a side is off by more than 10 mm, the fast track on that side that pulls the other way (the load track below the cable or the lift track above it) is switched into a position loop of 5 Hz after a sensor delay of 1 ms, and stays in it until the gap is back (`compensation_tracks.prot`: trigger 10 mm, delay 1 ms, 5 Hz; 'in every lane' means it is on in every simulated case, on both sides). It makes no standing heat; it acts only in an event.")
    L.append(f"- **\"The 0.1 mm regulation\" of the deployment.** Through the deployment ten of the 20 families of fast coils stay connected, and they hold the gap against the slow drift that the support coils do not follow between their weekly trims: the mean offset they work against is 0.1 mm of gap (33 N/m of pull per track at the constant-flux stiffness of 0.33 kN/m per millimeter) plus half a week of the casing's growing weight (44 N/m per track at 12.7 N/m a day). The signal is the same gap sensors. Its Joule heat in the leads is {DP['joule']:,.0f} W per site, beside {DP['cond']:,.0f} W of conduction through the ten families' leads: {DP['extra']:,.0f} W per site in all (`tools/round8_design_point.py --list`; the rule is round 6a's Section 4).\n")

    # 10
    L.append("### 10. The flux pumps\n")
    L.append("Confirmed: the round 7 design has no flux pump anywhere in the bearing. Stage A found that the published pumps give 0.05 to 1.2 mV, far too little for these coils, and round 4 replaced Stage B's transformer charger with leads. A support coil is charged and trimmed from a warm supply through a pair of copper leads that are withdrawn (or opened cold) when the coil is persistent. A strip of bare superconductor across the coil is the persistent switch: heated open to charge or trim, cooled closed to hold. A quench is dumped through the free-wheel diode with the dump resistor in series, the dump switch across the resistor opening to bring it in (round 6b, Section 1.2). The fast coils are not persistent at all; they are driven from their stores through their switches. So every 'flux pump' in Chapters 6, 7 and 8 (the definition box in 6.6, 'its own flux pump' in 7.5, Figure 7.3's caption, the pumps in 7.7 and Table 7.6, 'gap control, flux pumps' in Table 8.2 and 8.3) describes something the design does not have.\n")

    # 11
    L.append("### 11. The LIMs' lateral guidance\n")
    L.append("No script and no earlier report gives a force between a LIM stator and its reaction plate across the gap. `lim_physics.py`, `lim_simulation.py` and `lim_config.py` model the thrust along the track and the losses; a search of all three finds no normal, attraction, repulsion, lateral or centering term. `levitation_physics.py` has lateral terms, but for the old passive limiter plates, not for the LIMs. The bearing reports from Stage A to round 7 carry no lateral force either (round 6a notes that 'the lateral centering the support loops give the casing' is in no model). "
             f"For scale only: the LIMs' thrust is {J12['thrust_max_N']:,.0f} N each at the peak of the deployment, {2*J12['thrust_max_N']/500:.1f} N per meter of ring for the two of a site. The lateral loads of Chapter 2 are 0 for the equatorial ring and 10.9 kN/m at 5 degrees of inclination (Table 2.1), and the lost support loop's sideways force with inclined faces was 46 kN/m (round 7). "
             "Whether a LIM of this kind has a force across its gap far larger than its thrust is not something the model can say. So the book can state the lateral guidance as a requirement on the LIMs (they must center the cable and take the lateral loads), not as a result.\n")

    # 13
    L.append("### 13. Chapter 16's lighter casing (Section 16.6, Table 16.2)\n")
    L.append("**The three 6,000 kg/m cases were not rerun.** Each needs the bearing sized again for the lighter casing (the support at 23.75 + 198 kN/m, the plate trade, the fast tracks) and its own deployment runs, which this round did not do. They rest on the runs of 28 October (`Claude_Code_instructions_casing_6000_2026-10-28.md`): the old bearing of four 200 mm plates and the 150-turn coil. The design column alone, on the round 7 bearing and the 165-turn coil:\n")
    L.append("| Quantity | Book, the design case | New, the design case | Script |\n|---|---|---|---|")
    S9 = "`tools/round8_design_point.py`"
    L.append(f"| Hardware on the cable | 10,686 kg/m | 10,686 kg/m as the allowance; the bearing's plates and the LIM plates use {DP['hw_bearing']:,.0f} | Paul's ruling; `tools/round8_design_point.py --list` |")
    L.append("| Structural cable; cable with hardware; cable velocity; tension; deployment energy | 92,883 kg/m; 103,569 kg/m; 8,331 m/s; 683.0 GN; 19.2 EJ | unchanged | |")
    L.append(f"| Deployment at 12 MW per site | 381 days, 150 turns | {J12['t_end_days']:.0f} days, 165 turns | {S9} |")
    L.append(f"| Peak thrust per LIM | 1.01 kN | {J12['thrust_max_N']/1e3:.2f} kN | {S9} |")
    L.append(f"| First-day eddy heating per LIM; peak cable temperature | 77 kW; 331 K | {J12['p_eddy_max_kW']:.0f} kW; {J12['cable_temp_max_K']:.0f} K | {S9} |")
    L.append(f"| Levitation plates | Four 200 mm, 8,814 kg/m | Two support plates 1.2 m by 100 mm and four fast plates 1.0 m by 200 mm, {2*DP['m_s']+4*DP['m_l']:,.0f} kg/m | `tools/round8_design_point.py --list` |")
    L.append("| LIM reaction plates | Two 200 mm, 1,872 kg/m | unchanged | |")
    L.append(f"| Bearing on the casing, tape with formers and cryostats | About 1,030 kg/m | {m['tape']+m['formers']:,.0f} kg/m ({m['tape']:.0f} of tape); {m['total']:,.0f} with the stores, switches and leads | item 8 above |")
    L.append(f"| Standing cryogenic load per site (the cable's radiation onto the coils) | 10 kW | 10 kW (4.0 m of cold width), and {DP['cond']/1e3:.1f} kW through ten families' leads | {S9} |")
    L.append("| Anchor attachments per trolley line; per spreader arm | 8; 81 | unchanged (they follow the casing's weight) | |")
    L.append(f"\nWhat the text of 16.6 quotes from the old bearing and has no rerun: 'the bearing's tape falls only from 504 to 486 kg/m', '355 kg/m', 'the support plates stay at 200 mm', 'the lift plates can be 175 mm in every case ... a 550 kg/m remark', '10,136 kg/m', 'the bearing falls to 971 kg/m'. On the design case the tape is now {m['tape']:.0f} kg/m, the support plates are 100 mm and the fast plates 200 mm. The 193 days of the lighter cases would be near half of {J12['t_end_days']:.0f} with the 165-turn coil, but that is not run.\n")

    # 12
    L.append("### 12. Three plain descriptions to check the book's wording against\n")
    L.append("**What happens when a sled passes a place on the ring.** "
             f"The launch is planned, so the bearing knows when the spacecraft will reach each place. Under the casing floor and on its ceiling lie the four fast tracks: two load tracks below the cable, which pull the casing up, and two lift tracks above it, which pull it down. Each track is a row of coils 25 m long, in 20 staggered families, so that a new coil starts every 1.25 m and 20 coils lie over every point. Below 6,988 m/s the spacecraft presses the casing down and the load tracks act; above it the spacecraft lifts the casing and the lift tracks act. "
             f"As the spacecraft approaches, each coil in the zone around it closes its switches. Its store, a bank of film capacitors beside the coil, drives the current up through a pair of copper leads in 10 milliseconds, to {fast['I']:,.0f} A at full pull and {fast['V']/1e3:.1f} kV at the terminals. The pulled zone moves with the sled: {hd[1][1]:.0f} m long at 1 km/s and {hd[30][1]:,.0f} m at 30 km/s, so a place on the ring is pulled for about {hd[1][1]/1e3*1e3:.0f} ms and {hd[30][1]/30e3*1e3:.0f} ms. "
             f"The pull is {hd[1][0]/1e3:.0f} kN/m on each load track at 1 km/s and 196.5 kN/m on each lift track at 30 km/s. Its impulse cancels the spacecraft's kick, and the casing moves 12 to 14 mm and comes back. When the zone has passed, the switches turn the coil's energy back into its store, less the losses. "
             "Each side of the cable has its own sensors and is regulated separately, and if a gap is off by more than 10 mm the fast track on that side is switched into the gap protection. The support coils do nothing in the pass. They are persistent and keep their flux, which costs them about 12 percent of their pull while the load track beside them is at full pull; they get it back when the pull ends. For a launch all 20 families are connected. Between launches ten stay connected.\n")
    dn = ctx["dump"]
    L.append("**How a support coil is charged, held, trimmed and dumped.** "
             f"A support coil is one loop, 1 km long, wound of {n_t:.0f} turns, each turn ten layers of 12 mm tape soldered together. Ten such coils, staggered, make a loop, and two loops side by side make a track. Across the coil's two ends lies a strip of bare superconductor, and across the same ends a cold diode with a dump resistor in series, the resistor bridged by a switch. "
             f"To charge, a pair of copper leads is brought into contact with the coil's ends, the strip is heated until it is resistive, which takes about 0.7 s, and the warm supply drives the current up at about 210 V. The ten coils of a loop are charged together and share the loop's flux, so each sees {ctx['L_loop']:.1f} H, and the normal current of {P['I_n']:.0f} A takes about {ctx['L_loop']*P['I_n']/210:.0f} s from nothing ({(ctx['L_loop']+ctx['M_loop'])*P['I_n']/210:.0f} s if the loop beside it is charged at the same time). "
             "To hold, the strip is cooled until it is superconducting again. The current then runs around through the strip with no loss, the supply is brought to zero and the leads are withdrawn. Nothing is supplied while holding. "
             f"To trim, the leads are engaged, the strip is opened and the supply moves the current: 1 percent takes {ctx['L_loop']*P['I_n']*0.01/210:.1f} s at 210 V. After the loss of a loop the supply gives 1 kV and the change-over takes {ch['outer']['T_hi']:.1f} s. A parked spacecraft takes the track to its cap, {P['I_cap']:,.0f} A, in {t1k['T']:.0f} s at 1 kV. "
             f"To dump after a quench, the switch across the dump resistor opens, and the coil's current runs down through the diode and the resistor at about 210 V, in about 0.3 s. The nine other coils of the loop keep the loop's flux, so the dumped coil gives up only its leakage energy. With 0.8 s to detect the quench and open the strip, the hot spot reaches {dn[0]:.0f} K at the normal current and {dn[1]:.0f} K at the cap. "
             "(The dump's resistor, time and temperatures are round 7's, computed on the figure-eight coil at 455 A before the split; they were not repeated for the single loop at 371 A, where the spot is cooler.)\n")
    L.append("**What happens when a support loop is lost.** "
             f"The casing hangs from four loops, two under each side of the cable, each carrying a quarter of its 95.0 kN/m. Suppose a fault empties one loop, an outer one. Three things follow, on three time scales. "
             f"At once, with no one doing anything, the loop beside the lost one answers. The two loops of a track share flux, and a superconducting loop keeps its flux, so its current rises by itself by 21 percent, from {P['I_n']:.0f} to {ch['outer']['I0'][2]:.0f} A. That side now gives {ch['F_jump']/1e3:.1f} kN/m where it gave 47.5: {(W_side(ch) )/1e3:.1f} kN/m is missing, not 23.75. "
             f"Within milliseconds the fast coils take up the rest. The load track on the lost side pulls {ctx['first'][0]/1e3:.1f} kN/m and the one on the other side {ctx['first'][1]/1e3:.1f} kN/m, which holds the casing level and stops it rolling. The casing sags less than 3 mm and rolls 0.4 milliradian. The fast coils' leads are sized for pulses, so they cannot hold this for long: with ten families connected the budget is {P['budget10']:.2f} s of full pull. "
             f"Within about a second the supplies of the three sound loops are connected and their strips opened, and at 1 kV they move the currents to where no fast coil is needed: the lost side's remaining loop to {ch['outer']['I1'][2]:.0f} A, 1.70 times normal, and the other side's two loops down a little, to {ch['outer']['I1'][0]:.0f} A. That takes {ch['outer']['T_hi']:.1f} s and uses {ch['outer']['int_hi']:.2f} s of the budget. The strips are closed, the leads withdrawn, and the fast coils go back to rest. "
             f"The ring can stay like that for as long as the repair takes. The remaining loop's worst tape is at about {ctx['share_trim']*100:.0f} percent of its limit. If the lost loop is an inner one the numbers are a little smaller ({ch['inner']['I1'][2]:.0f} A, {ch['inner']['int_hi']:.2f} s).\n")
    return L


def W_side(ch):
    return 47500.0 - ch["F_jump"]
