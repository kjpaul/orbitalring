#!/usr/bin/env python3
"""
The numbers sheet of round 8 (2026-11-15): round 7's sheet on the book's cable.  The cable keeps its hardware allowance of 10,686 kg/m (Paul's
ruling), so there are no cable rows; the deployment rows come from the runs of tools/round8_design_point.py, the Chapter 7 rows from round 7's
solves, and Section 3b (the rows for the rewrite of Sections 7.5 to 7.9) from bearing_round8.py through tools/round8_sheet_3b.py.
Book: Orbital Ring Engineering, Volume III: Chapters 1, 3, 4, 7 to 12, 14 to 16, 19 and the appendix.
Run: python tools/round8_sheet.py [--out=FILE]
Reads runs/round8_design/ (tools/round8_design_point.py), runs/round7_logs/ (bearing_round7.py), runs/round8_logs/ (bearing_round8.py) and runs
tools/anchor_lines_deployment.py and tools/rising_radius_hanging.py.  The book's values are typed from the chapter files as they stand (they are
quotations, not results).  Environment: ORBITAL_RING_OUT (default runs).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import csv
import json
import math
import os
import re
import subprocess
import sys

import ring_altitude as ring
import round5_sheet as s5
import round8_sheet_3b as b3
import rising_radius_design_165 as d

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
R7 = os.path.join(HERE, OUT, "round8_design"); LOGS = os.path.join(HERE, OUT, "round7_logs"); LOGS8 = os.path.join(HERE, OUT, "round8_logs")
SHEET = os.path.join(HERE, "..", "..", "Vol_III", "design_point_numbers_round8_2026-11-15.md")
HW_BOOK = 10686.0
COP = 0.0106; RAD_W_M2 = 1306.0; SITES = 88946


def jload(*p):
    path = os.path.join(*p)
    return json.load(open(path, encoding="utf-8")) if os.path.exists(path) else None


def tool(args, env=None):
    e = dict(os.environ, PYTHONUTF8="1"); e.update(env or {})
    return subprocess.run([sys.executable] + args, cwd=HERE, env=e, capture_output=True, text=True, encoding="utf-8", errors="replace").stdout


def series(path):
    rows = list(csv.DictReader(open(path)))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def main():
    out = next((a[6:] for a in sys.argv if a.startswith("--out=")), SHEET)
    DP = jload(R7, "design_point.json")
    hw = DP["hw"]; hw_b = DP["hw_bearing"]
    R8 = jload(LOGS8, "round8_rows.json")
    d.BASE = R7; d.TURNS_REF = os.path.join(R7, "slip_optimal")
    A_book = {}; A_new = {}
    rb = ring.design(m_hw=HW_BOOK); rn = rb
    fast = jload(LOGS, "round7_fast_track.json"); heat = jload(LOGS, "round7_heat.json"); sup = jload(LOGS, "round7_support.json")
    half = jload(LOGS, "round7_half_2p5.json") or (jload(LOGS, "round7_half.json") or {}).get("0.050")
    half_h = "2.5 mm" if jload(LOGS, "round7_half_2p5.json") else "5 mm"
    cpl = jload(LOGS, "round7_coupling.json"); dump = jload(LOGS, "round7_dump.json")
    J12 = d.load("slip", "optimal", "N165_P12MW"); J8 = d.load("ch9", "optimal", "L2_N165_P8MW")
    L = []
    def table(title, head="| Where | Book | New | Script |\n|---|---|---|---|"):
        L.append(f"\n## {title}\n"); L.append(head)
    def row(where, book, new, script):
        L.append(f"| {where} | {book} | {new} | {script} |")

    L.append("# The numbers sheet for round 8: the design point on the book's cable (2026-11-15)\n")
    L.append(f"Asked for in `Vol_III/Claude_Code_instructions_round8_book_cable_2026-11-15.md`. **The cable is unchanged.** By Paul's ruling it keeps its hardware allowance of 10,686 kg/m, so every cable number of the book stands "
             f"(92,883 kg/m structural, 103,569 kg/m with hardware, 8,330.8 and 8,847.0 m/s, 683.0 GN, 54.64 m², 7.39 m, 4,606 Mt of cable, 5,140 Mt of ring) and this sheet has no cable rows. "
             f"The bearing of round 7 uses {hw_b:,.0f} kg/m of the allowance; the {hw - hw_b:,.0f} kg/m between them is the hardware's margin.\n")
    L.append(f"The design point: the ring of `ring_altitude.py` (700 km, the rising quadratic, the casing at 12,000 kg/m, the cable's hardware 10,686 kg/m); the LIM coil at 165 turns, 12 MW per site, two LIMs, the `optimal` slip rule; "
             f"the bearing as ruled after round 7: two support tracks, each of two independent 0.4 m loops 50 mm apart at the center on a 1.2 m by 100 mm plate, and four fast tracks (two load, two lift), each a 0.6 m loop with 204 by 50 mm bundles on a 1.0 m by 200 mm plate; "
             f"4.0 m of cold width; ten of the 20 fast families connected through the deployment and between launches, all 20 for a launch; the 21 mm² lead. "
             f"The extra cold load of the ten connected families is **{DP['extra']:,.0f} W per site**: {DP['cond']:,.0f} W of conduction and {DP['joule']:,.0f} W of the 0.1 mm regulation's Joule heat, at the accepted fast track's own {DP['I_full']:,.0f} A (round 7's runs used 4,406 W). "
             f"'Book' is the chapter text as it stands today: Chapter 9 is at 165 turns with no bearing heat (361.9 days), every other chapter at 150 turns and 381 days. 'New' is the run on this design point. "
             f"Every new number comes from a script run, named in the last column; where no script gives a number the row says so. No book chapter is edited.\n")
    L.append("@@CH134@@")

    # ------------------------------------------------------------------ Chapter 7
    table("Chapter 7: the bearing")
    r = fast; lv_note = "`bearing_round7.py --part=solve-fast`, `--part=fast` (2.5 mm cells)"
    row("Table 7.1 and 7.5: the cable's hardware", "10,686 kg/m; 'The four plates together weigh approximately 8,800 kg per meter of cable'", f"10,686 kg/m stays as the allowance. The bearing uses {hw_b:,.0f} kg/m: LIM plates 1,872, two support plates of {DP['m_s']:,.0f}, four fast plates of {DP['m_l']:,.0f} ({2*DP['m_s']+4*DP['m_l']:,.0f} for the six). The {hw - hw_b:,.0f} kg/m between them, {(hw - hw_b)/hw*100:.0f} percent, is the hardware's margin", "Paul's ruling (round 8, Section 1); `tools/round8_design_point.py --list`")
    row("7.5, the plates: lift and load", "1.4 m wide, 200 mm thick, four plates, 'about 8,800 kg per meter of cable'", f"1.0 m by 200 mm, {DP['m_l']:,.0f} kg/m each, four (two lift, two load)", "Paul's ruling (round 7, Section 0); round 6b Section 2.1")
    row("7.5, the plates: support", "1.4 m wide, 200 mm thick (the same four plates)", f"1.2 m by 100 mm, {DP['m_s']:,.0f} kg/m each, two", "the same")
    row("7.5: the lift track's loop and bundles", "'a single loop 0.8 m wide, with two bundles 200 mm wide and 50 mm tall'", "a single loop 0.6 m wide between leg centers, two bundles 204 mm wide and 50 mm tall, ten cells across and two up, the end and outer cells with their tapes' faces vertical", lv_note)
    row("7.5: the fast coils' length", "50 m on the lift tracks, 25 m on the load windings", "25 m on all four fast tracks, 20 staggered families", "round 4's ruling")
    row("7.7: the fast track's ampere-turns per leg at 196.5 kN/m", "421,000", f"{r['NI_req']:,.0f}", lv_note)
    row("7.7: the worst 3 mm tape (30 A at half)", "60 A limit", f"{r['ic']:.1f} A at {r['B']:.2f} T, {r['th']:.0f} degrees from its normal; {r['ic']/2:.1f} A at half", lv_note)
    row("7.6: the fast track's tape", "60 kg of 3 mm tape per meter", f"{r['tape']:.1f} kg/m, {r['tapes']:,.0f} tapes per leg", lv_note)
    row("7.6: the fast track's rating", "288 kN/m against 196.5", f"{r['F70']/1e3:.0f} kN/m, {r['margin']:.2f} times the requirement", lv_note)
    row("the fast plate at the requirement", "0.1 percent saturated, 1.53 T (the book's plate)", f"{r['sat']*100:.1f} percent saturated, peak {r['B_plate']:.2f} T", lv_note)
    row("7.7: a fast coil: turns, current, voltage", "not quoted (MVA only)", f"10 turns of {r['layers_whole']} layers, {r['I']:,.0f} A, {r['V']/1e3:.1f} kV at the terminals ({r['V_turn']:,.0f} V per turn)", lv_note)
    row("7.7: the rise", "2 ms (the protection), 20 ms (the planned ramp)", "10 ms", "Paul's ruling (round 5)")
    row("7.7, Equation 7.5: a coil's stored energy and supply", "50 m lift coil 0.30 MJ, 27 MVA, 301 MVA at 2 ms; 25 m load coil 0.15 MJ, 13 and 146 MVA", f"25 m coil {r['E_coil']/1e6:.3f} MJ, {2*r['E_coil']/10e-3/1e6:.0f} MVA at the 10 ms rise", lv_note)
    row("7.7: the film between turns; the cell; the bundle's height", "not quoted", f"turns 1 to 5 below, 6 to 10 above; {r['t1']*1e3:.2f} mm between neighbors, {5*r['t1']*1e3:.2f} mm between the rows; walls {r['t_fam']*1e3:.2f} mm; the rows need {r['pack']['rows in order']['widths'][0]*1e3:.1f} and {r['pack']['rows in order']['widths'][1]*1e3:.1f} mm of the 204 and {r['pack']['rows in order']['height']*1e3:.1f} mm of the 50", "`bearing_round7.py --part=fast`; round 6b Section 1.3")
    row("7.7 Losses: the tape's pulse heat per launch", "605 GJ, 57 TJ at the compressor, 1.2 percent; 207 GJ at 20 degrees, 53 GJ at 5 degrees", f"{(r['Ep']+r['Es'])/1e9:.0f} GJ ({r['Ep']/1e9:.0f} through the faces, {r['Es']/1e9:.0f} the soldered stacks), {(r['Ep']+r['Es'])/COP/1e12:.0f} TJ at the compressor, {(r['Ep']+r['Es'])/COP/4.6e15*100:.2f} percent of the launch's 4.6 PJ; {(r['Ep']+r['Es'])*7/86400/ring.L_RING:.2f} W/m at 7 a day, {(r['Ep']+r['Es'])*14/86400/ring.L_RING:.2f} at 14", lv_note)
    h = heat
    row("7.7 Losses, 7.8: the plate heat of a launch", "1,090 J per meter of plate per pass at 30 km/s, 238 GJ, 0.005 percent; 5.4 kJ per meter; 31 W/m", f"{h['rows'][-1]['hyst']+h['rows'][-1]['eddy']:.0f} J per meter of plate per pass at 30 km/s; {h['E']/1e9:.0f} GJ per launch ({h['E_mean']/1e9:.0f} with each column's mean field), {h['E']/4.6e15*100:.4f} percent; {h['E']/ring.L_RING/1e3:.1f} kJ per meter of ring; {h['q7']:.2f} W/m at 7 a day, {h['q14']:.2f} at 14; the cable {h['T14']-h['T0']:.2f} K warmer at 14 a day", "`bearing_round7.py --part=heat`")
    hc, hn, ha = half["cap"], half["normal"], half["alone"]
    n_t = hc["tapes_leg"] / 100.0
    row("7.5: the support track's loops and bundles", "'two loops side by side, each 0.5 m wide', bundles 150 by 67 mm, the center pair touching", f"two separate loops, each 0.4 m between leg centers, each its own coil in its own cryostat, 50 mm between the center bundles; bundles 150 by 67 mm with centers at -500, -100, +100 and +500 mm", f"`bearing_round7.py --part=solve-half` ({half_h} cells), `--part=support`")
    row("Equation 7.2, Tables 7.4 and 7.6, 7.6: the support's cap, ampere-turns per leg at 245 kN/m", "322,000", f"{hc['NI']:,.0f}", "the same")
    row("7.6, Table 7.3: the worst 12 mm tape at the cap", "1.64 T at 53 degrees, 234 A, 117 A at half", f"{hc['B']:.2f} T at {hc['th']:.0f} degrees, {hc['ic']:.0f} A, {hc['ic']/2:.0f} A at half", "the same")
    row("7.5, 7.6, Equation 7.3: tapes per leg; tape per track", "2,753; 94 kg/m", f"{hc['tapes_leg']:,.0f}; {hc['tape']:.1f} kg/m", "the same")
    row("7.6: the plate at the cap", "11 percent saturated", f"{hc['sat']*100:.0f} percent saturated, peak {hc['B_plate']:.2f} T", "the same")
    row("Tables 7.3, 7.4 and 7.6, 7.6: normal operation, ampere-turns per leg and tape current", "141,000; 51 A, about 15 percent of the limit", f"{hn['NI']:,.0f} (solved at 47.5 kN/m, not scaled from the cap); {hn['NI']/hc['tapes_leg']:.1f} A per tape, {hn['share']*100:.0f} percent of the limit there", "the same")
    row("7.5: a support coil", "275 turns per loop per coil (Stage B), then 55 turns of ten layers (round 4)", f"{n_t:.1f} turns of ten layers per loop, one coil per loop; {hn['NI']/(10*n_t):.0f} A normal, {hc['NI']/(10*n_t):,.0f} A at the cap", "the same")
    row("Table 7.6, 'one track lost' and 'one track out for repair'", "the other track at 199,000 ampere-turns, 72 A, 95 kN/m", f"one loop lost: its neighbor's current rises 21 percent by itself, then the three others are trimmed (the lost side's loop to {ha['NI']/hn['NI']*math.sqrt(51.6/47.5):.2f} times normal, the other side's two to 0.96); one loop alone carrying its track's 47.5 kN/m needs {ha['NI']/hn['NI']:.2f} times the normal current, its worst tape at {ha['share']*100:.0f} percent of its limit", "`bearing_round7.py --part=support`, `--part=dyn`")
    row("the load track's pull on the support beside it", "not quoted", f"{cpl['c_split']*100:.1f} percent of the support's normal force at the load track's full pull, at the 0.2 m spacing (round 6b's 9.2 percent was on a scaled normal current)", "`bearing_round7.py --part=coupling`")
    row("Table 7.5: the casing's movement at 1 / 3 / 5 / 10 / 20 / 30 km/s", "7.6 / 2.6 / 1.3 / 1.7 / 6.0 / 9.7 mm", "12.0 / 10.6 / 10.3 / 10.2 / 10.8 / 13.8 mm (the 10 ms rise, the 10 mm protection in every lane, the 11.7 percent coupling)", "`bearing_round7.py --part=cases` (round 4's cases)")
    row("7.7 Tolerances: the smallest gap in any single failure", "50 mm (a pull 50 ms late)", "24.3 mm (a pull 50 ms late), at the 10 ms rise", "the same")
    dn = next(x for x in dump["rows"] if x.get("held") == 0.8 and abs(x["I"] - dump["coil"]["I0"]) < 1 and x["R"] < 1)
    dc = next(x for x in dump["rows"] if x.get("held") == 0.8 and x["I"] > 1000 and x["R"] < 0.6)
    row("7.5, the support coil's protection", "flux pumps (the chapter); the dump (round 4)", f"the dump resistor in series with the diode, its switch across it; a dump at the charging voltage ({dn['V']:.0f} V, {dn['R']:.2f} ohm, {dn['tau']*1e3:.0f} ms): hot spot {dn['T']:.0f} K at the normal current, {dc['T']:.0f} K at the cap, with 0.8 s to detect and open the strip", "`bearing_round7.py --part=dump`; round 6b Section 1.2")
    tape_casing = 2 * hc["tape"] + 4 * r["tape"]
    row("7.6: tape per meter of casing; the bearing's weight", "94 + 104 kg/m per support track, 60 per lift track; 'about 1,030 kg for each meter of casing'", f"two support tracks {hc['tape']:.1f}, four fast tracks {r['tape']:.1f}: {tape_casing:.0f} kg/m of tape, {2*tape_casing:,.0f} kg/m with formers, cryostats and supports at one tape mass (Part 3b's rule); {tape_casing*ring.L_RING/1e9:.2f} Mt of tape around the ring", "arithmetic on the two solves")

    # ------------------------------------------------------------------ Section 3b
    import bearing_round4 as r4, bearing_round5 as r5
    hold = {u: (r4.pull_track_of_u(u * 1e3)[0], r5.hold_time(u * 1e3)[1]) for u in (1, 30)}
    dn0 = next(x for x in dump["rows"] if x.get("held") == 0.8 and abs(x["I"] - dump["coil"]["I0"]) < 1 and x["R"] < 1)
    dc0 = next(x for x in dump["rows"] if x.get("held") == 0.8 and x["I"] > 1000 and x["R"] < 0.6)
    hm = sup["half"]["0.050"]; tr_o = next(x for x in sup["rows"] if "50 mm" in x["name"] and x["loss"] == "one loop, others trimmed")
    share_trim = next(o["share"] for o in R8["ops"] if o["case"].startswith("an outer loop lost"))
    L3b = b3.rows_3b(R8, dict(half=half, fast=fast, J12=J12, DP=DP, dump=(dn0["T"], dc0["T"]), L_loop=hm["L_loop"], M_loop=hm["M_loop"], first=tr_o["first"], share_trim=share_trim, hold=hold))

    # ------------------------------------------------------------------ Chapter 8
    table("Chapter 8: the stores, the switches and the standing loads")
    W = r["W"]; E_ring = 4 * W; cap_kJ = 1.5
    n_s = math.ceil(r["V"] / 2800.0); m_sw = 4 * n_s * 4.3
    row("8.1: the store for the sled's passage", "'SMES ... is the natural fit', no number", f"film capacitors, warm, at every coil: {r['E_coil']/1e6:.3f} MJ per coil, {E_ring/1e3:.0f} kJ per meter of ring, {E_ring/cap_kJ/1e3:.0f} kg/m at 1.5 kJ/kg; shared by a site's four tracks {E_ring/2/cap_kJ/1e3:.0f} kg/m", "`bearing_round7.py --part=fast`; round 5 Section 3's rule")
    row("8.1: the pulse power of the pulled zone at 30 km/s", "not quoted", f"{W*30e3/1e9:.1f} GW per lift track at the front, {2*W*30e3/1e9:.1f} for the pair", "the same (the stored energy per meter times the speed)")
    row("the switches per coil", "not quoted", f"{n_s} StakPak IGBTs in series per position at {r['V']/1e3:.1f} kV, two switches and two diodes: {m_sw:.0f} kg per coil, {m_sw*3.2:.0f} kg/m of casing", "round 5 Section 4's rule on the new voltage")
    ser = series(os.path.join(R7, "ch14", "ch14_thermal_series.csv"))
    q0 = ser["q_cryo_cold_W"][1] - ser["p_hyst_site_W"][1]                # the standing load at the start: leak, levitation coils, the ten families' leads
    # after the deployment the casing's weight no longer grows, so of the regulation's mean pull (33 N/m for the 0.1 mm, 44 N/m for half a week of
    # weight drift) only the first part is left; its Joule heat goes as one over the number of connected families
    import bearing_round4 as r4s
    p_reg = r4s.K_FLUX * 1e-3 * 0.1; p_drift = 7 * 12.7 / 2
    j10 = DP["joule"] * p_reg / (p_reg + p_drift); j20 = j10 / 2
    q_leads10 = DP["cond"] + j10; q_leads20 = 2 * DP["cond"] + j20
    q_stand_dep = q0 - DP["extra"] + q_leads10                            # ten families connected: the book's standing figure
    q_stand_op = q0 - DP["extra"] + q_leads20                             # all 20 connected: a launch window
    qser = ser["q_cryo_cold_W"]; q_pk = max(qser)
    tser = [x / (30 * 86400.0) for x in ser["t_s"]]; mid8 = [i for i in range(len(tser)) if 6 <= tser[i] <= 10]; q_mid = sum(qser[i] for i in mid8) / len(mid8)
    def after(qw):
        P = qw / COP * SITES / 1e9; left = (1067.0 - P) * 0.95
        return f"{qw/1e3:.1f} kW cold and {qw/COP/1e6:.2f} MW per site, {P:.0f} GW around the ring; left for export {left:.0f} GW ((1,067 less {P:.0f}) less 5 percent for transmission, the rule behind Table 8.2's 'about 900 GW')"
    row("8.3 (L106), Table 8.2 'Operations, Cryogenic system', 12.10, 19.1, 19.2: the standing plant after the deployment, ten fast families connected (the book's standing figure)", "about 1.4 MW per site (1.37), 122 GW over 88,946 sites", after(q_stand_dep) + f". The site's cold load at the start of the run less the hysteresis is {q0/1e3:.2f} kW: leak, the cable's radiation onto 4.0 m of cold width, {DP['cond']:,.0f} W of conduction in ten families' leads and {DP['joule']:,.0f} W of the regulation's Joule heat. After the deployment the casing's weight no longer grows, so the regulation's heat falls to {j10:,.0f} W (the 0.1 mm part alone)", "`tools/ch14_thermal_series.py` on the design point (`tools/round8_design_point.py --part=ch14`); the regulation's split by round 6a's rule")
    row("the same with all 20 families connected, a launch window (before the launch's own heat)", "not in the book", after(q_stand_op), f"the same; {DP['cond']:,.0f} W more conduction per site and the regulation's heat halved")
    row("Table 8.2 'Export, 83 m panels'", "About 900 GW", f"about {(1067.0 - q_stand_dep/COP*SITES/1e9)*0.95:.0f} GW with ten families connected; about {(1067.0 - q_stand_op/COP*SITES/1e9)*0.95:.0f} GW in a launch window", "arithmetic on the two rows above")
    row("Table 8.2 'Deployment, Cryogenic cooling'", "2.5 MW of the 12 MW at each site through the second half of the deployment, 2.9 MW at the peak (Table 14.2)", f"{q_mid/COP/1e6:.2f} MW through months 6 to 10, {q_pk/COP/1e6:.2f} MW at the peak (Chapter 14 below)", "`tools/ch14_thermal_series.py` on the design point")
    row("Table 8.2 'Operations, Levitation'; 8.3: 'gap control, flux pumps, and the coil-end drag of 5 to 183 MW'", "flux pumps; 5 to 183 MW", f"no flux pumps (Section 3b, item 10); the coil-end drag {R8['ends'][0]['MW']:.0f} MW for the design's 1 km staggered coils (Section 3b, item 7)", "`bearing_round8.py --part=rows`")

    # ------------------------------------------------------------------ Chapter 9
    table("Chapter 9: Table 9.2", "| Row | Book | New | Script |\n|---|---|---|---|")
    S9 = "`tools/round8_design_point.py` (`lim_simulation.py --rising-radius`)"
    import lim_hysteresis as lh
    cl = lh.coil_limit("flat", 165, d.LAYERS, d.W_TAPE)
    row("Turns in each coil", "165", "165 (the ruling); see Table 9.4 for the fastest", S9)
    row("Peak coil current", "183 A", f"{J12['I_TARGET']:.0f} A", S9)
    row("Current limit of one tape in the field of the winding", "183 A, in about 0.07 T", f"{cl[0]:.0f} A, in {cl[1]['B']:.2f} T", "`lim_hysteresis.coil_limit`")
    row("Tape in each LIM", "1,293 kg (2,585 kg per site), 303 km", f"{J12['tape_kg_lim']:,.0f} kg ({2*J12['tape_kg_lim']:,.0f} kg per site), {J12['tape_km_lim']:,.0f} km", S9)
    row("Slip at the end", "0.4 %", f"{J12['final_slip_pct']:.2f} %", S9)
    row("Supply frequency at the end", "44 Hz", f"{J12['final_f_supply_Hz']:.1f} Hz", S9)
    row("Peak thrust of each LIM", "about 1,230 N", f"{J12['thrust_max_N']:,.0f} N", S9)
    row("Deployment time", "362 days at 12 MW, 512 days at 8 MW", f"{J12['t_end_days']:.0f} days at 12 MW" + (f", {J8['t_end_days']:.0f} days at 8 MW" if J8 else ", (8 MW not yet run)"), S9)

    table("Chapter 9: Table 9.3", "| Quantity | Book, 12 MW | New, 12 MW | Book, 8 MW | New, 8 MW |\n|---|---|---|---|---|")
    def q3(lab, b12, f, b8):
        L.append(f"| {lab} | {b12} | {f(J12)} | {b8} | {f(J8) if J8 else '(not run)'} |")
    def lim_txt(J):
        s = f"current to day {J['current_cap_last_day']:.0f}"
        if J.get("power_cap_first_day") is not None: s += f"; site power from day {J['power_cap_first_day']:.0f} to day {J['power_cap_last_day']:.0f}"
        s += f"; 100 kV from day {J['volt_limit_day']:.0f} to the end" if J.get("volt_limit_day") is not None else "; never 100 kV"
        return s
    def T_txt(J):
        ms = J.get("T_cable_months", []); late = [T for m, T in ms if m >= 6]
        return f"{J['cable_temp_max_K']:.0f} K on the first day" + (f", settling near {late[-1]:.0f} K" if late else "")
    q3("Deployment time", "362 days (11.9 months)", lambda J: f"{J['t_end_days']:.0f} days ({J['t_end_days']/30.4375:.1f} months)", "512 days (16.8 months)")
    q3("Governing limit, in order", "current to day 92; site power to day 279; 100 kV from day 269", lim_txt, "current to day 60; site power to the end; never 100 kV")
    q3("Coil voltage at the end", "100 kV, 105 A", lambda J: f"{J['final_volts_V']/1e3:.0f} kV, {J['final_current_A']:.0f} A", "87 kV, 91 A")
    q3("Site power at the end", "falling to 10.3 MW", lambda J: f"{J['final_site_power_MW']:.1f} MW", "8 MW")
    q3("Peak thrust per LIM; at the end", "1,230 N at one month; 401 N", lambda J: f"{J['thrust_max_N']:,.0f} N at {J['thrust_max_day']/30.4375:.1f} months; {J['final_thrust_N']:.0f} N", "1,230 N; 303 N")
    q3("Peak thrust power per LIM", "4.2 MW at nine months", lambda J: f"{J['p_thrust_max'][0]:.1f} MW at {J['p_thrust_max'][2]/30.4375:.1f} months" if "p_thrust_max" in J else f"{J['final_thrust_N']*J['final_v_rel_m_s']/1e6:.1f} MW at the end (the thrust times the speed past the casing)", "2.7 MW at the end")
    q3("Mechanical energy delivered", "19.2 EJ", lambda J: f"{J['E_lim_loop_EJ']:.2f} EJ", "19.2 EJ")
    q3("Peak eddy current heating per LIM", "95 kW on the first day", lambda J: f"{J['p_eddy_max_kW']:.0f} kW", "95 kW")
    q3("Peak hysteresis heating per LIM", "7.1 kW at three months", lambda J: f"{J['p_hyst_max_kW']:.1f} kW" + (f" at {J['p_hyst_max'][2]/30.4375:.1f} months" if "p_hyst_max" in J else ""), "4.4 kW at two months")
    q3("Peak cable temperature", "345 K, settling near 262 K", T_txt, "345 K, settling near 255 K")
    q3("Peak cryogenic plant per site (the book's COP)", "2.75 MW (9.8)", lambda J: f"{J['cryo_peak_cop_MW']:.2f} MW ({J['q_cryo_cold_max_kW']:.1f} kW cold)", "")

    table("Chapter 9: Table 9.4, the turns at 12 MW (and the 8 MW cases of the text)", "| Turns | Book: days; 100 kV from; peak thrust; thrust at the end; tape per site | New | Script |\n|---|---|---|---|")
    book94 = {120: "628; never; 650 N; 448 N; 1.88 t", 135: "449; never; 820; 460; 2.12", 150: "381; never; 1,010; 469; 2.35", 160: "364; day 305; 1,150; 427; 2.51", 165: "362; day 269; 1,230; 401; 2.59",
              170: "363; day 239; 1,300; 378; 2.66", 180: "371; day 190; 1,460; 337; 2.82", 208: "429; day 108; 1,950; 253; 3.26", 240: "533; day 64; 2,510; 190; 3.76", 280: "701; day 37; 3,010; 139; 4.39",
              140: "not in the table", 145: "not in the table", 155: "not in the table", 175: "not in the table (363 days in Section 18)", 190: "not in the table"}
    pts = {12: [], 8: []}
    for name, (lims, N, mw) in sorted(d.SLIP_CASES.items(), key=lambda kv: (-kv[1][2], kv[1][1])):
        J = d.load("slip", "optimal", name)
        if not J:
            L.append(f"| {N} at {mw} MW | | (not yet run) | {S9} |"); continue
        pts[mw].append((N, J["t_end_days"]))
        vl = J.get("volt_limit_day")
        new = f"{J['t_end_days']:.0f} days ({J['t_end_days']:.1f}); {'never' if vl is None else 'day %.0f' % vl}; {J['thrust_max_N']:,.0f} N; {J['final_thrust_N']:.0f} N; {J['tape_t_site']:.2f} t"
        book = book94.get(N, "") if mw == 12 else {150: "529 days (Section 18; not in the chapter)", 190: "496 days (9.7)", 208: "505 days (9.7)"}.get(N, "not quoted")
        L.append(f"| {N}{'' if mw == 12 else ' at 8 MW'} | {book} | {new} | {S9} |")
    for mw in (12, 8):
        if pts[mw]:
            fl = min(x for _, x in pts[mw]); best = [N for N, x in pts[mw] if x == fl][0]; within = [N for N, x in pts[mw] if x <= 1.02 * fl]
            L.append(f"\nAt {mw} MW the fastest coil of the grid has {best} turns, {fl:.1f} days; within 2 percent: {', '.join(str(n) for n in sorted(within))} turns."
                     + (f" The 165-turn design is {next(x for N, x in pts[mw] if N == 165):.1f} days, {(next(x for N, x in pts[mw] if N == 165)/fl-1)*100:.2f} percent above the fastest." if mw == 12 and any(N == 165 for N, _ in pts[mw]) else ""))

    table("Chapter 9: Table 9.5, the site power and the number of LIMs", "| Row (site power, turns, LIMs) | Book: days; 100 kV from; panel; tape per site; cable | New | Script |\n|---|---|---|---|")
    book95 = {"L2_N165_P8MW": "512 (16.8 months); never; 55 m; 2.59 t; 4,606 Mt", "L2_N165_P12MW": "362 (11.9); day 269; 83 m; 2.59 t; 4,606 Mt", "L2_N165_P16MW": "338 (11.1); day 163; 110 m; 2.59 t; 4,606 Mt",
              "L2_N165_P20MW": "337 (11.1); day 143; 118 m (20 or 24 MW)", "L2_N165_P24MW": "the same row", "L2_N150_P12MW": "381 (12.5); never; 83 m; 2.35 t", "L2_N150_P16MW": "355 (11.7); day 251; 110 m; 2.35 t",
              "L2_N208_P12MW": "429 (14.1); day 108; 82 m; 3.26 t", "L2_N208_P16MW": "417 (13.7); day 65; 110 m; 3.26 t", "L4_N165_P8MW": "510 (16.8); never; 55 m; 5.17 t; 4,670 Mt",
              "L4_N165_P16MW": "245 (8.0); never; 110 m; 5.17 t; 4,670 Mt", "L4_N165_P24MW": "181 (5.9); day 127; 165 m; 5.17 t; 4,670 Mt", "L6_N165_P24MW": "162 (5.3); never; 165 m; 7.76 t; 4,734 Mt"}
    rows95 = {}
    for name in ["L2_N165_P8MW", "L2_N165_P12MW", "L2_N165_P16MW", "L2_N165_P20MW", "L2_N165_P24MW", "L2_N150_P12MW", "L2_N150_P16MW", "L2_N208_P12MW", "L2_N208_P16MW", "L4_N165_P8MW", "L4_N165_P16MW", "L4_N165_P24MW", "L6_N165_P24MW"]:
        J = d.load("ch9", "optimal", name)
        lims, N, mw = d.CH9_CASES[name]
        if not J:
            L.append(f"| {mw} MW, {N}, {lims} | {book95[name]} | (not yet run) | {S9} |"); continue
        rows95[name] = J; vl = J.get("volt_limit_day")
        L.append(f"| {mw} MW, {N}, {lims} | {book95[name]} | {J['t_end_days']:.0f} days ({J['t_end_days']/30.4375:.1f} months); {'never' if vl is None else 'day %.0f' % vl}; {J['panel_m']:.0f} m; {J['tape_t_site']:.2f} t; "
                 f"{J['cable_Mt']:,.0f} Mt; cryogenic plant {J['cryo_peak_cop_MW']:.2f} MW; plate {J['plate_t_max_K']:.0f} K; hysteresis {J['p_hyst_max_kW']:.1f} kW per LIM; the site draws at most {J['site_power_max_MW']:.1f} MW | {S9} |")
    pc = sorted([(J["P_MW"], J["t_end_days"]) for n, J in rows95.items() if n.startswith("L2_N165")])
    table("Chapter 9: the text of Sections 9.6 to 9.8")
    row("9.6: 'on day 269 it reaches 100 kV', 'falls to 105 A on the last day', '10.3 MW at the end'", "day 269; 105 A; 10.3 MW", f"day {J12['volt_limit_day']:.0f}; {J12['final_current_A']:.0f} A; {J12['final_site_power_MW']:.1f} MW", S9)
    row("9.6: the thrust, 'about 1,230 N', '401 N at the end', '218 MN around the ring'", "1,230 N; 401 N; 218 MN", f"{J12['thrust_max_N']:,.0f} N; {J12['final_thrust_N']:.0f} N; {J12['thrust_max_N']*2*SITES/1e6:.0f} MN", S9)
    row("9.6: 'On day 362 ... the cable reaches 8,331 m/s and 683 GN'", "day 362; 8,331 m/s; 683 GN", f"day {J12['t_end_days']:.0f}; {J12['final_v_cable_m_s']:,.0f} m/s; {J12['final_tension_GN']:.0f} GN", S9)
    if J8:
        row("9.6: '8 MW per site takes 512 days ... first 60 days ... ends the run at 87 kV. The extra 150 days'", "512 days; 60 days; 87 kV; 150 days", f"{J8['t_end_days']:.0f} days; {J8['current_cap_last_day']:.0f} days; {J8['final_volts_V']/1e3:.0f} kV; {J8['t_end_days']-J12['t_end_days']:.0f} days", S9)
    row("9.6: the energy, '19.2 EJ ... 11.6 EJ is kinetic ... 7.2 EJ ... 0.4 EJ'", "19.2; 11.6; 7.2; 0.4 EJ", f"{J12['E_lim_loop_EJ']:.2f}; {J12['dKE_cable_EJ']+J12['dKE_casing_EJ']:.1f}; {J12['dPE_EJ']:.1f}; {J12['elastic_EJ']:.1f} EJ", S9)
    row("9.6: the heating, 'about 95 kW for each LIM', 'the slip is 0.4%'", "95 kW; 0.4 %", f"{J12['p_eddy_max_kW']:.0f} kW; {J12['final_slip_pct']:.2f} %", S9)
    row("9.6: the cable, '345 K ... about 262 K at the end'", "345 K; 262 K", T_txt(J12), S9)
    J150 = d.load("slip", "optimal", "N150_P12MW"); J208 = d.load("slip", "optimal", "N208_P12MW")
    if J150:
        row("9.7: 'At 12 MW that coil has 150 turns. It ends the run at 98 kV and takes 381 days'", "150 turns; 98 kV; 381 days", f"150 turns; {J150['final_volts_V']/1e3:.0f} kV; {J150['t_end_days']:.0f} days", S9)
        row("9.7: '401 N against 469 N ... 19 days faster'", "401 against 469 N; 19 days", f"{J12['final_thrust_N']:.0f} against {J150['final_thrust_N']:.0f} N; {J150['t_end_days']-J12['t_end_days']:.0f} days", S9)
    if J208:
        row("9.7: 'a 208-turn coil pushes 1,950 N ... 253 N', '429 days, 67 days longer'", "1,950 N; 253 N; 429 days, 67 longer", f"{J208['thrust_max_N']:,.0f} N; {J208['final_thrust_N']:.0f} N; {J208['t_end_days']:.0f} days, {J208['t_end_days']-J12['t_end_days']:.0f} longer", S9)
    J190 = d.load("slip", "optimal", "N190_P8MW"); J208_8 = d.load("slip", "optimal", "N208_P8MW")
    if J190 and J208_8:
        row("9.7, Figure 9.3: 'The fastest coil has 165 turns at 12 MW and 190 turns at 8 MW'; '190 turns ... 496 days ... 208 turns ... 505 days'", "165 and 190 turns; 496; 505 days", f"see the two sentences under Table 9.4 above; 190 turns {J190['t_end_days']:.0f} days, 208 turns {J208_8['t_end_days']:.0f} days at 8 MW", S9)
    row("9.7: '2.59 t for each site against 2.35 t ... 345 K against 331 K'", "2.59 against 2.35 t; 345 against 331 K", f"{J12['tape_t_site']:.2f} against " + (f"{J150['tape_t_site']:.2f} t; {J12['plate_t_max_K']:.0f} against {J150['plate_t_max_K']:.0f} K" if J150 else "(150 turns not run)"), S9)
    row("9.7: 'loses about 7 kW per LIM at the peak'", "7 kW", f"{J12['p_hyst_max_kW']:.1f} kW", S9)
    if len(pc) > 1:
        fl = min(x for _, x in pc); d12 = next(x for p, x in pc if p == 12)
        steps = "; ".join(f"{(x0-x1)/(p1-p0):.1f} days per MW from {p0:.0f} to {p1:.0f}" for (p0, x0), (p1, x1) in zip(pc, pc[1:]))
        row("9.8: 'saves 150 days, about 37 days for each megawatt ... 24 more, 6 days ... The floor is 337 days ... The 12 MW design sits 7% above'", "37 and 6 days per MW; floor 337 days; 7 %", f"{', '.join(f'{x:.0f} days at {p:.0f} MW' for p, x in pc)}; {steps}; the floor {fl:.0f} days; 12 MW sits {(d12/fl-1)*100:.1f} percent above it; the knee by the 2 percent rule {min(p for p, x in pc if x <= 1.02*fl):.0f} MW", S9)
    for nm, bk in (("L4_N165_P8MW", "9.8: '510 days against 512 with two'"), ("L4_N165_P16MW", "9.8: '245 days at 16 MW'"), ("L4_N165_P24MW", "9.8: '181 days at 24 MW'"), ("L6_N165_P24MW", "9.8: '162 days'")):
        if nm in rows95:
            row(bk, bk.split("'")[1], f"{rows95[nm]['t_end_days']:.0f} days", S9)
    if "L4_N165_P16MW" in rows95 and "L6_N165_P24MW" in rows95:
        row("9.8: 'cryogenic plant of 3.7 to 4.7 MW at the peak, against 2.75 MW for the design ... 404 K ... 447 K'", "3.7 to 4.7 MW against 2.75; 404 K; 447 K", f"{min(rows95[n]['cryo_peak_cop_MW'] for n in rows95 if n.startswith(('L4', 'L6'))):.1f} to {max(rows95[n]['cryo_peak_cop_MW'] for n in rows95 if n.startswith(('L4', 'L6'))):.1f} MW against {J12['cryo_peak_cop_MW']:.2f}; {rows95['L4_N165_P16MW']['plate_t_max_K']:.0f} K; {rows95['L6_N165_P24MW']['plate_t_max_K']:.0f} K", S9)
        d4 = ring.design(m_hw=hw + 1872.0)
        row("9.8: '4,670 Mt with four LIMs per site against 4,606'; '12.6 m/s less per pair ... 3.2 GN less tension'", "4,670 against 4,606 Mt; 12.6 m/s; 3.2 GN", f"{rows95['L4_N165_P16MW']['cable_Mt']:,.0f} against {J12['cable_Mt']:,.0f} Mt; {rn['v_cable']-d4['v_cable']:.1f} m/s less; {abs(d4['tension']-rn['tension'])/1e9:.1f} GN less", "`ring_altitude.design`")

    # ------------------------------------------------------------------ Chapter 10
    table("Chapter 10")
    row("the whole chapter", "the cable at 683.0 GN and 103,569 kg/m, 8,847 m/s past the casing, the wave speed 2,568 m/s, the circuit time 17,319 s, the 190 minute launch", "nothing changes: the cable is the book's, and `lsm_simulation.py` takes only the ring from it", "no run needed")

    # ------------------------------------------------------------------ Chapter 11
    table("Chapter 11")
    S11 = "`tools/anchor_lines_deployment.py` on the design run's tension series"
    tfile = os.path.join(R7, "slip_optimal", "N165_P12MW", "rising_radius_deployment_tension.csv")
    at = tool(["tools/anchor_lines_deployment.py", tfile])
    hg = tool(["tools/rising_radius_hanging.py", os.path.join(R7, "slip_optimal", "N165_P12MW", "rising_radius_daily.csv")])
    def g(pat, txt, cast=float):
        m = re.search(pat, txt); return cast(m.group(1)) if m else None
    T = series(tfile)
    k0 = next(i for i, v in enumerate(T["v_casing_m_s"]) if v <= 0.0)
    row("11.1: 'by 23 m/s on the first day and 16 m/s a day over the whole deployment'; '2.7e-4 m/s2 on the first day'", "23 m/s; 16 m/s a day; 2.7e-4 m/s2 (on the chapter's own 150-turn run the tool gives 11.4 m/s on the first day, 1.3e-4 m/s2, and 21.3 a day over the whole deployment: the sentence's figures match neither and are to be checked)", f"{g(r'over the first 1 day.s.: ([\d.]+) m/s per day', at):.1f} m/s on the first day ({g(r'over the first 1 day.s.: [\d.]+ m/s per day = ([\d.e+-]+) m/s2', at, str)} m/s2), {g(r'over the first 10 day.s.: ([\d.]+) m/s per day', at):.1f} a day over the first ten, {g(r'average ([\d.]+) m/s per day', at):.1f} a day over the whole deployment", S11)
    row("11.1: the casing's net weight per kilogram, '0.03 m/s2 after one day, 0.2 after a week, a tenth of its final value on day 20, 1.4 after 30 days, 4.9 after 100 days'", "0.03; 0.2; day 20; 1.4; 4.9", f"{g(r'day   1: ([\d.]+)', at):.2f}; {g(r'day   7: ([\d.]+)', at):.2f}; day {g(r'on day ([\d.]+) .text: day 30', at):.0f}; {g(r'day  30: ([\d.]+)', at):.1f}; {g(r'day 100: ([\d.]+)', at):.1f} m/s2", S11)
    hang = {int(m.group(1)): int(m.group(2)) for m in re.finditer(r"a hang of (\d+) km stays above 150 km from day (\d+)", hg)}
    row("11.1: '400 km by day 108, 450 km by day 203, 500 km by day 292, and the full 550 km only at the end'", "days 108, 203, 292, the end", f"days {hang.get(400)}, {hang.get(450)}, {hang.get(500)}, {hang.get(550)}", "`tools/rising_radius_hanging.py` on the design run's daily table")
    row("11.1: 'which it is from about day 342, two days before it passes through rest on day 344'", "day 342; day 344", f"day {g(r'that happens on day ([\d.]+)', at):.0f}; day {T['time_days'][k0]:.0f}", S11 + "; the first row of the series with the casing at or past rest")
    row("11.1: 'The casing reaches ground speed on day 381'", "day 381", f"day {J12['t_end_days']:.0f}", S9)
    row("11.1: the lowering spool (550 km): from the first day; after the deployment", "not quoted in days; 8.0 MWh and 177 kN, 59.6 MWh and 894 kN in the 2026-11 sheet", f"{g(r'from day 0 .ends day [\d.]+, casing weight then [\d.]+ m/s2.: brake energy ([\d,]+) kWh', at, str)} kWh and {g(r'from day 0 .*?peak drum tension ([\d.]+) kN', at):.0f} kN from day 0; {g(r'after the deployment .*?= ([\d.]+) MWh', at):.1f} MWh and {g(r'after the deployment.*?peak drum tension ([\d.]+) kN', at):.0f} kN after", S11)
    row("11.2: 'well inside the 396 kN/m that the support tracks allow'", "396 kN/m", "the load tracks' pair, 2 x 198 kN/m, unchanged as a requirement (the load compensation is on its own tracks since Stage B)", "`final_bearing.py`")
    row("the anchor lines' sizes", "sized from the casing's weight and the lateral force, not from the cable", "not rerun: `anchor_line_analysis.py` takes neither the hardware nor the cable's mass", "no change expected; not run")

    # ------------------------------------------------------------------ Chapter 12
    table("Chapter 12")
    row("12.11 (L261): '377 km of arm can hang from the first day ..., 400 km by day 108, 450 km by day 203, 500 km by day 292, and the full 550 km only at the end'", "377 km; days 108, 203, 292; the end (day 381)", f"{J12['build_altitude_km'] - 150:.0f} km; days {hang.get(400)}, {hang.get(450)}, {hang.get(500)}; the end (day {J12['t_end_days']:.0f})", "`tools/rising_radius_hanging.py` on the design run's daily table")
    row("12.10 (L227): 'The cryogenic plant alone takes about 122 GW (Chapter 8)'", "122 GW", f"{q_stand_dep/COP*SITES/1e9:.0f} GW with ten fast families connected ({q_stand_op/COP*SITES/1e9:.0f} GW in a launch window)", "`tools/ch14_thermal_series.py` on the design point")

    # ------------------------------------------------------------------ Chapter 14
    S14 = "`tools/ch14_thermal_series.py --rising` on the design point; `tools/fig14_figures.py`"
    t = [x / (30 * 86400.0) for x in ser["t_s"]]; q = ser["q_cryo_cold_W"]
    ip = max(range(len(q)), key=lambda i: q[i])
    def at_m(m): return min(range(len(t)), key=lambda i: abs(t[i] - m))
    mid = [i for i in range(len(t)) if 6 <= t[i] <= 10]
    def mean(key): return sum(ser[key][i] for i in mid) / len(mid)
    def plant(qw): return qw / COP, qw + qw / COP, (qw + qw / COP) / RAD_W_M2, (qw + qw / COP) / RAD_W_M2 / 450.0
    lev = 5.0 * DP["cold_width"] * 500.0
    leak = q0 - lev - DP["extra"]
    table("Chapter 14: Table 14.2 and the cold loop", "| Row | Book: peak (five months); steady state; after deployment | New: peak (" + f"{t[ip]:.1f} months); steady state (months 6 to 10); after deployment (ten families' leads; all 20) | Script |\n|---|---|---|---|")
    hp = ser["p_hyst_site_W"][ip]; hm = mean("p_hyst_site_W")
    def r14(lab, book, a, b, c1, c2, f="{:.1f} kW", k=1e3):
        L.append(f"| {lab} | {book} | {f.format(a/k)}; {f.format(b/k)}; {f.format(c1/k)}; {f.format(c2/k)} | {S14} |")
    r14("Hysteresis in the two LIMs", "16.3; 12.0; 0 kW", hp, hm, 0.0, 0.0)
    r14("Leak through the cryostat insulation", "4.5; 4.5; 4.5 kW", leak, leak, leak, leak)
    r14("Heat onto the levitation coils", "10.0; 10.0; 10.0 kW", lev, lev, lev, lev)
    r14("The fast coils' leads, conduction and the regulation's Joule heat (new row)", "not in the table", DP["extra"], DP["extra"], q_leads10, q_leads20, "{:.2f} kW")
    qm = mean("q_cryo_cold_W")
    r14("Cold-side load", "30.8; 26.5; 14.5 kW", q[ip], qm, q_stand_dep, q_stand_op)
    r14("Compressor power", "2.90; 2.50; 1.37 MW", plant(q[ip])[0], plant(qm)[0], plant(q_stand_dep)[0], plant(q_stand_op)[0], "{:.2f} MW", 1e6)
    r14("Heat rejected at 400 K", "2.93; 2.53; 1.38 MW", plant(q[ip])[1], plant(qm)[1], plant(q_stand_dep)[1], plant(q_stand_op)[1], "{:.2f} MW", 1e6)
    r14("Cryogenic radiator area", "2,247; 1,935; 1,058 m2", plant(q[ip])[2], plant(qm)[2], plant(q_stand_dep)[2], plant(q_stand_op)[2], "{:,.0f} m2", 1.0)
    r14("Cryogenic radiator width over 450 m", "4.99; 4.30; 2.35 m", plant(q[ip])[3], plant(qm)[3], plant(q_stand_dep)[3], plant(q_stand_op)[3], "{:.2f} m", 1.0)
    table("Chapter 14: the text, Tables 14.1 and 14.3, Figures 14.1 and 14.2")
    Tc = ser["cable_temp_K"]; Tw = ser["wall_temp_K"]; pe = ser["p_eddy_site_W"]
    i1 = max(range(min(len(Tc), at_m(0.5))), key=lambda i: Tc[i])
    row("Table 14.1: eddy currents in the reaction plates, per site: peak; steady state", "155 kW; 41 kW", f"{max(pe)/1e3:.0f} kW; {pe[-1]/1e3:.0f} kW (the end of the run)", S14)
    row("Table 14.1: hysteresis, per site: peak; steady state", "16 kW; 12 kW", f"{max(ser['p_hyst_site_W'])/1e3:.1f} kW; {hm/1e3:.1f} kW", S14)
    row("Table 14.3: eddy heat from the LIMs, W/m: steady state; first-day peak", "82 W/m; 310 W/m", f"{pe[-1]/500:.0f} W/m; {max(pe)/500:.0f} W/m", S14)
    row("Table 14.3, 14.2, 14.4, 14.6: the casing wall: steady state; first-day peak", "263 K; 322 K", f"{Tw[-1]:.0f} K; {max(Tw):.0f} K", S14)
    row("14.2, Figure 14.1, 14.6: the cable's temperature", "331 K on the first day, 278 K, 313 K at four months, near 271 K, the wall 263 K; 12.7 months", f"{max(Tc):.0f} K on the first day, {Tc[at_m(1.5)]:.0f} K at 1.5 months, {max(Tc[at_m(2.0):at_m(4.0)]):.0f} K at {t[max(range(at_m(2.0), at_m(4.0)), key=lambda i: Tc[i])]:.1f} months, {Tc[at_m(6.0)]:.0f} K at six months, {Tc[-1]:.0f} K at the end with the wall at {Tw[-1]:.0f} K; {t[-1]*30/30.4375:.1f} months", S14)
    row("14.2: 'cools as soon as the power limit takes over on day 133'", "day 133 (the 150-turn coil)", f"the site reaches 12 MW on day {J12['power_cap_first_day']:.0f}; the current cap ends on day {J12['current_cap_last_day']:.0f}", S9)
    row("14.3, Figure 14.2, 14.6: the cold-side load", "near 15 kW at the start, 30.8 kW at 5.0 months, 24.4 kW at the end; 4.99 m and 3.95 m", f"{q[1]/1e3:.1f} kW at the start, {q[ip]/1e3:.1f} kW at {t[ip]:.1f} months, {q[-1]/1e3:.1f} kW at the end; {plant(q[ip])[3]:.2f} m and {plant(q[-1])[3]:.2f} m", S14)
    row("14.3: '2,000 m2 and 10 kW at a 500 m site' (the levitation coils' cold width)", "4.0 m of cold width, 10 kW", f"{DP['cold_width']:.1f} m (two 0.4 m loops per support track, one 0.6 m loop per fast track), {lev/1e3:.0f} kW", "`tools/round8_design_point.py --list`")
    row("14.3, 14.6: 'about 1.4 MW at every site' after deployment; 'about 14.5 kW'", "14.5 kW, 1.4 MW", f"{q_stand_dep/1e3:.1f} kW and {plant(q_stand_dep)[0]/1e6:.2f} MW with ten families' leads connected; {q_stand_op/1e3:.1f} kW and {plant(q_stand_op)[0]/1e6:.2f} MW with all 20", S14)
    row("14.4: the radiator, '5.0 m wide ... about 4.3 m ... after the deployment 2.4 m'", "5.0; 4.3; 2.4 m", f"{plant(q[ip])[3]:.1f}; {plant(qm)[3]:.1f}; {plant(q_stand_dep)[3]:.1f} m with ten families connected ({plant(q_stand_op)[3]:.1f} m with all 20)", S14)
    row("14.6: '2.9 MW of compressor power'", "2.9 MW (2.75 MW in Chapter 9)", f"{plant(q[ip])[0]/1e6:.2f} MW", S14)
    cl = R8["cold_launch"]; E_c = cl["E_tape"] + cl["E_lead"]; spare = q[ip] - q_stand_op
    row("14.5: the sled's pulse in the cold, '5.74 kW on the loss factor of the 2025 book, about 2 kW on the strip model'", "5.74 kW; 2 kW", f"not rerun; the heaviest launch's cold heat is in the next rows", "")
    row("14.5 (new): the cold heat the heaviest launch leaves, in all", "605 GJ as the upper bound for the tape (7.7)", f"{cl['E_tape']/1e9:.0f} GJ in the fast coils' tape and {cl['E_lead']/1e9:.1f} GJ of Joule heat in their leads, {E_c/1e9:.0f} GJ; {E_c/COP/1e12:.0f} TJ at the compressor, {E_c/COP/4.6e15*100:.1f} percent of the launch's 4.6 PJ", "`bearing_round7.py --part=fast`; `bearing_round8.py --part=rows` (the leads: `bearing_round5.launch_integrals`)")
    row("14.5 (new): the same per meter of ring along the run, and for one 500 m site", "not in the book", f"the run is {cl['x_run']/1e3:,.0f} km in {cl['t_run']/60:.0f} minutes, {cl['laps']:.1f} laps of the ring, so every meter is passed {cl['laps']:.1f} times; summed over the passes, as a mean around the ring: {cl['per_m_tape']/1e3:.1f} kJ/m in the tape and {cl['per_m_lead']/1e3:.2f} kJ/m in the leads; {E_c/ring.L_RING*500/1e6:.1f} MJ per 500 m site", "the same; how the heat divides between the laps (each at its own speed) is inside the launch integral and is not tabulated")
    row("14.5 (new): how long a site's plant takes to remove it at its spare capacity", "not in the book", f"the plant is sized for the deployment's peak of {q[ip]/1e3:.1f} kW cold and carries {q_stand_op/1e3:.1f} kW in a launch window (all 20 families connected), so {spare/1e3:.1f} kW is spare: {E_c/ring.L_RING*500/spare:.0f} s, {E_c/ring.L_RING*500/spare/60:.1f} minutes, for one launch's heat; at 14 of the heaviest a day the launches add {E_c*14/86400/ring.L_RING*500/1e3:.2f} kW per site on average, {E_c*14/86400/ring.L_RING*500/spare*100:.0f} percent of the spare. The rise of the nitrogen's temperature while the heat waits is in no script", "arithmetic on the rows above and Table 14.2")

    # ------------------------------------------------------------------ Chapter 15 and Table 1.1
    table("Chapters 15 and 16")
    row("15.5 (L98, L106): 'about 380 days at 12 MW per LIM site (Chapter 9)'; 'about 380 days each way'", "about 380 days", f"about {round(J12['t_end_days'], -1):.0f} days ({J12['t_end_days']:.1f})", S9)
    row("15.5: 'about 19 EJ'", "about 19 EJ", f"{J12['E_lim_loop_EJ']:.2f} EJ, unchanged", S9)
    row("16.6 (L121): 'the deployment at 12 MW per site to 193 days against 381'; Table 16.2's 'Deployment at 12 MW per site, 150 turns: 381 days'; its footnote 'the 150-turn, 12 MW, two-LIM design'", "381 days; 193 days; 150 turns", f"{J12['t_end_days']:.0f} days at 165 turns for the design case; the three lighter cases are not rerun (Section 3b, item 13)", S9)
    row("Table 16.2: 'Peak thrust per LIM 1.01 kN'; 'First-day eddy heating per LIM; peak cable temperature: 77 kW; 331 K'", "1.01 kN; 77 kW; 331 K", f"{J12['thrust_max_N']/1e3:.2f} kN; {J12['p_eddy_max_kW']:.0f} kW; {J12['cable_temp_max_K']:.0f} K (the design case)", S9)
    row("Table 16.2: 'Levitation plates'; 'Bearing on the casing'; 'Standing cryogenic load per site'", "four 200 mm, 8,814 kg/m; about 1,030 kg/m; 10 kW", "see Section 3b, item 13", "")

    # ------------------------------------------------------------------ Chapter 19 and the appendix
    table("Chapter 19 and the appendix: the numbers they repeat")
    lt = open(os.path.join(R7, "slip_optimal", "N165_P12MW", "rising_radius_lim_simulation.txt"), encoding="utf-8", errors="replace").read()
    ke = re.search(r"cable kinetic energy change\s+([\d.-]+)\s+\(([\d.]+) to ([\d.]+)\)", lt); kc = re.search(r"casing kinetic energy change\s+([\d.-]+)\s+\(([\d.]+) to ([\d.]+)\)", lt)
    row("19.1, Table 19.1, 19.2: '381 days at 12 MW per site with a 150-turn coil that ends the run just under 100 kV'; 'Deployment time (381 days)'; '19.2 EJ over 381 days'", "381 days; 150 turns; 19.2 EJ", f"{J12['t_end_days']:.0f} days; 165 turns, at 100 kV from day {J12['volt_limit_day']:.0f}; {J12['E_lim_loop_EJ']:.1f} EJ", S9)
    row("19.1, 19.2: '1.37 MW per site and 122 GW in all'; 'the 122 GW of the cryogenic plant'", "1.37 MW; 122 GW", f"{plant(q_stand_dep)[0]/1e6:.2f} MW; {plant(q_stand_dep)[0]*SITES/1e9:.0f} GW with ten fast families connected ({plant(q_stand_op)[0]/1e6:.2f} MW and {plant(q_stand_op)[0]*SITES/1e9:.0f} GW in a launch window)", S14)
    row("Table 19.2: 'thinner than 200 mm, and whether the lift plates can go to 175 mm (550 kg/m off the cable's hardware)'", "200 mm plates; 175 mm; 550 kg/m", f"the fast plates are 1.0 m by 200 mm and the support plates 1.2 m by 100 mm; the bearing uses {hw_b:,.0f} kg/m of the 10,686 kg/m allowance, so the 550 kg/m remark goes", "round 6b Section 2.1; Paul's ruling")
    if ke and kc:
        row("Appendix A.2, the energy balance: the start, 'the cable and the casing both at 7,598 m/s, is 148.34 EJ: 132.94 EJ in the cable and 15.40 EJ in the casing'", "148.34; 132.94; 15.40 EJ", f"{float(ke.group(2))+float(kc.group(2)):.2f}; {float(ke.group(2)):.2f}; {float(kc.group(2)):.2f} EJ", S9 + " (the log's energy budget)")
        row("Appendix A.2: the end, 'the cable has 159.84 EJ and the casing 0.07 EJ'", "159.84; 0.07 EJ", f"{float(ke.group(3)):.2f}; {float(kc.group(3)):.2f} EJ", "the same")
        row("Appendix A.2: '11.57 EJ; the potential energy of the 173 km rise adds 7.24 EJ and the elastic energy ... 0.37 EJ, 19.18 EJ in all'", "11.57; 7.24; 0.37; 19.18 EJ", f"{float(ke.group(1))+float(kc.group(1)):.2f}; {J12['dPE_EJ']:.2f}; {J12['elastic_EJ']:.2f}; {J12['E_lim_loop_EJ']:.2f} EJ", "the same")
    row("Appendix (L23) and 3.x (L148): 'integrates the power over the 381 days'", "381 days", f"{J12['t_end_days']:.0f} days", S9)

    # ------------------------------------------------------------------ the figures
    table("The figures that need redrawing", "| Figure | What changes | Redrawn | Drawn by |\n|---|---|---|---|")
    FW = "`Vol_III/figure_work_2026-11-15/`"
    L.append(f"| 3.2, the deployment tension | the design run on the book's cable: {J12['t_end_days']:.0f} days, {rn['tension']/1e9:.0f} GN, the casing at rest on day {T['time_days'][k0]:.0f}; the prograde dashed curve from the prograde run of the same design point | yes, {FW}`fig03_02_deployment_tension_700` | `tools/round8_book_figures.py` on `Vol_III/figure_scripts_2026-11-07/gen_fig03_02_rising.py` |")
    L.append(f"| 9.3, the turns | the sweep at 12 and 8 MW on the book's cable | yes, {FW}`fig09_03_turns_sweep` | `tools/round8_book_figures.py` on `rising_radius_design_165.figure93` |")
    L.append(f"| 14.1, the cable and the wall | {max(Tc):.0f} K on the first day, {Tc[-1]:.0f} and {Tw[-1]:.0f} K at the end, {t[-1]:.1f} months of 30 days | yes, {FW}`fig14_01_cable_wall_temperature_700` | `tools/round8_book_figures.py` on `tools/fig14_figures.py` |")
    L.append(f"| 14.2, the cold-side load | {q[ip]/1e3:.1f} kW and {plant(q[ip])[3]:.2f} m at {t[ip]:.1f} months, {q[-1]/1e3:.1f} kW and {plant(q[-1])[3]:.2f} m at the end | yes, {FW}`fig14_02_cold_load_radiator_700` | the same |")
    L.append(f"| 7.4 and the Chapter 7 candidates | the plates as ruled, the split support, the inset | yes, {FW}`fig07_candidate_1_cross_section` | `tools/round7_figures.py` |")
    L.append(f"| 1.1, 4.1, 4.6, 5.4, 8.1 | nothing: the cable and the 12 MW per site are the book's | no | |")
    L.append("\nRound 7's versions of Figures 3.2, 9.3, 14.1 and 14.2 (on the 10,057 kg/m cable) are in `Vol_III/figure_work_2026-11-15/round7_cable_10057/`.")
    C = []
    def crow(where, book, new, script):
        C.append(f"| {where} | {book} | {new} | {script} |")
    def at_day(dd, key):
        i = min(range(len(T["time_days"])), key=lambda i: abs(T["time_days"][i] - dd)); return T[key][i]
    wts = T["casing_net_weight_kN_m"]; w_end = wts[-1]
    i_half = next(i for i, v in enumerate(wts) if v >= w_end / 2); i_pk = max(range(len(wts)), key=lambda i: wts[i])
    C.append("\n## Chapter 1\n"); C.append("| Where | Book | New | Script |\n|---|---|---|---|")
    crow("Table 1.1 (L169): Deployment", "381 days at 12 MW per LIM site, 1,067 GW in all, 19.2 EJ", f"{J12['t_end_days']:.0f} days at 12 MW per LIM site, 1,067 GW in all, {J12['E_lim_loop_EJ']:.1f} EJ", S9)
    crow("1.1 (L27): 'takes about 380 days at 12 MW per LIM site (Chapter 9)'", "about 380 days", f"about {round(J12['t_end_days'], -1):.0f} days ({J12['t_end_days']:.1f})", S9)
    crow("Table 1.1: Deployment altitude", "527 km, rising 173 km to 700 km", f"{J12['build_altitude_km']:.0f} km, rising {J12['rise_km']:.0f} km: unchanged", S9)
    C.append("\n## Chapter 3: Section 3.2, Figure 3.2 and the summary\n"); C.append("| Where | Book | New | Script |\n|---|---|---|---|")
    S3 = "the design run's tension series (`rising_radius_deployment_tension.csv`)"
    crow("3.2 (L62): 'at 4.6 GN after 30 days'", "4.6 GN", f"{at_day(30.0, 'tension_GN'):.1f} GN", S3)
    crow("3.2: 'while the LIMs push at their hardest, 180 MN around the ring'", "180 MN (the 150-turn coil's 1,010 N)", f"{J12['thrust_max_N']*2*SITES/1e6:.0f} MN ({J12['thrust_max_N']:,.0f} N per LIM, at {J12['thrust_max_day']:.0f} days)", S9)
    crow("3.2, Figure 3.2's caption: 'to half of its final 95 kN/m by day 79'", "day 79", f"day {T['time_days'][i_half]:.0f}", S3)
    crow("3.2: 'peaks near 96.5 kN/m around day 314 before it eases to its final value'", "96.5 kN/m; day 314", f"{wts[i_pk]:.1f} kN/m; day {T['time_days'][i_pk]:.0f}", S3)
    crow("3.2, the caption: 'The casing passes through rest on day 344, at 679 km with the cable at 601 GN and 11.0 GPa'", "day 344; 679 km; 601 GN; 11.0 GPa", f"day {T['time_days'][k0]:.0f}; {T['altitude_km'][k0]:.0f} km; {T['tension_GN'][k0]:.0f} GN; {T['stress_GPa'][k0]:.1f} GPa", S3)
    crow("3.2, the caption, the alt text: 'reaches ground speed on day 381'; 'over the 381 days'; 'The deployment ends on day 381, twelve and a half months after it began'", "day 381; twelve and a half months", f"day {J12['t_end_days']:.0f}; {J12['t_end_days']/30.4375:.1f} months (twelve months)", S9)
    crow("3.2: 'with the cable at 8,331 m/s and 683 GN, which is 12.5 GPa'", "8,331 m/s; 683 GN; 12.5 GPa", f"{J12['final_v_cable_m_s']:,.0f} m/s; {J12['final_tension_GN']:.0f} GN; {T['stress_GPa'][-1]:.1f} GPa: unchanged", S9)
    crow("3.x (L148): 'integrates the LIM power over the 381 days'; the summary (L262): 'the deployment takes 381 days and 19.2 EJ at 12 MW per LIM site'", "381 days; 19.2 EJ", f"{J12['t_end_days']:.0f} days; {J12['E_lim_loop_EJ']:.1f} EJ", S9)
    crow("Figure 3.2", "the 150-turn run, 381 days", "redrawn from the design run and a prograde run on the book's cable (see the figures at the end)", "`tools/round8_book_figures.py`")
    C.append("\n## Chapter 4\n"); C.append("| Where | Book | New | Script |\n|---|---|---|---|")
    crow("the place that quotes 381", "L414: 'a ring that operates at 600 km is constructed at 430 km and one at 550 km at 381 km'", "nothing to change: a search of the chapter finds this one '381', and it is a build altitude in kilometers, not the deployment's days. The chapter does not quote the deployment's days", "a search of `Vol_III_Ch04.md` (read only)")
    k = L.index("@@CH134@@"); L[k:k + 1] = C
    k = next(i for i, x in enumerate(L) if x.startswith("\n## Chapter 8"))
    L[k:k] = L3b
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {out}: {len(L)} lines")


if __name__ == "__main__":
    main()
