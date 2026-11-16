#!/usr/bin/env python3
"""
The numbers sheet of round 9 (2026-11-16, Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md, item 2): round 8's sheet, row for row,
on the corrected model.  'Book' is round 8's value (the chapters were re-keyed to round 8's sheet on 15 November); 'New' is round 9's value, or
'unchanged' where the two are the same, so that the book session need not check that row.
Book: Orbital Ring Engineering, Volume III: Chapters 1, 3, 7 to 12, 14 to 16, 19 and the appendix.
Run: python tools/round9_sheet.py [--out=FILE]
How: tools/round8_sheet.py is run on round 9's runs (runs/round9_design/, the bearing's cold load from each run's own bearing_cold_load.csv) into
runs/round9_logs/round9_sheet_raw.md; the rows that round 9 computes differently (the cases, the bearing's cold-load rows, Chapter 11 in full)
are replaced or added here; then each row is set against the same row of Vol_III/design_point_numbers_round8_2026-11-15.md.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import csv
import json
import os
import re
import sys

import ring_altitude as ring
import round8_sheet as s8

HERE = s8.HERE; OUT = s8.OUT
R9 = os.path.join(HERE, OUT, "round9_design"); L9 = os.path.join(HERE, OUT, "round9_logs")
SHEET8 = os.path.join(HERE, "..", "..", "Vol_III", "design_point_numbers_round8_2026-11-15.md")
SHEET9 = os.path.join(HERE, "..", "..", "Vol_III", "design_point_numbers_round9_2026-11-16.md")
RAW = os.path.join(L9, "round9_sheet_raw.md")
COP = s8.COP; SITES = s8.SITES


def bcl_series(folder):
    rows = list(csv.DictReader(open(os.path.join(folder, "bearing_cold_load.csv"), encoding="utf-8")))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def raw():
    """Round 8's sheet generator on round 9's runs."""
    B = bcl_series(os.path.join(R9, "ch14")); heat = json.load(open(os.path.join(L9, "round9_heat.json"), encoding="utf-8"))
    j0 = s8.jload

    def jload(*p):
        v = j0(*p)
        if p[-1] == "design_point.json":
            v = dict(v); v["extra"] = B["extra_W"][0]            # the bearing's cold load at the start of the run (the standing rows subtract it)
        if p[-1] == "round8_rows.json":
            v = dict(v); v["cold_launch"] = dict(v["cold_launch"], E_lead=heat["E_lead"], per_m_lead=heat["per_m_lead"], lead_s=heat["lead_s"])
        return v
    s8.jload = jload; s8.R7 = R9
    argv = sys.argv; sys.argv = [argv[0], f"--out={RAW}"]
    try:
        s8.main()
    finally:
        sys.argv = argv; s8.jload = j0
    return B


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def parse(path):
    """{(section, first cell): cells} for every table row, and the lines."""
    lines = open(path, encoding="utf-8").read().split("\n"); sec = ""; rows = {}
    for ln in lines:
        if ln.startswith("## "):
            sec = re.sub(r"\(.*?\)", "", ln[3:]).strip()
        elif ln.startswith("|") and not ln.startswith("|---"):
            c = cells(ln); rows[(sec, c[0])] = c
    return rows, lines


def fmt_row(c):
    return "| " + " | ".join(c) + " |"


def pump_section():
    """Addendum A2: the default pump beside the pump that carries the whole rise, and the dynamo."""
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "runs", "round9_design")
    cols = []
    for n in ("ch14", "ch14_whole", "ch14_dynamo"):
        S = json.load(open(os.path.join(base, n, "rising_radius_summary.json"), encoding="utf-8"))
        B = json.load(open(os.path.join(base, n, "bearing_cold_load.json"), encoding="utf-8"))
        T = json.load(open(B["table"], encoding="utf-8"))
        cols.append((S, B, T))
    row = lambda lab, f: fmt_row([lab] + [f(*c) for c in cols])
    L = ["", "## The flux pump of the support coils: the default beside a pump that carries the whole rise (addendum A2)", "",
         "`tools/round9_design_point.py --part=ch14` (the default), `--part=ch14 --pump=whole`, `--part=ch14 --pump=dynamo`; the pump tables from `tools/round9_cold_load.py`. "
         "The second column's voltage is 2 percent over the largest daily rise of the support current after the first day of the design run (9.85 A a day, from day 1 to day 2). "
         "The first day itself rises 22.5 A, which would take 5.6 mV into the two loops of a coil (4.6 mV into one loop), so the leads stay on for the first two days in that column. "
         "Its loss per pump is kept at the default's 10 W [UNVERIFIED for a pump of ten times the voltage]. The brief's estimate was about 3.3 mV.", "",
         fmt_row(["Quantity", "Default: 0.267 mV pump", "A pump that carries the whole rise", "Dynamo, 0.197 mV (second case)"]), "|---|---|---|---|",
         row("Pump voltage into a support coil (21.54 H)", lambda S, B, T: f"{T['pump_V']*1e3:.3g} mV"),
         row("Current the pump moves", lambda S, B, T: f"{B['pump_rate_A_per_day']:.2f} A a day"),
         row("The support's leads come off", lambda S, B, T: f"day {B['leads_off_day']:.1f}"),
         row("Deployment", lambda S, B, T: f"{S['t_end_days']:.1f} days"),
         row("Peak cold load per site", lambda S, B, T: f"{S['q_cryo_cold_max_kW']:.1f} kW on day {S['p_cryo_max_day']:.0f}"),
         row("Radiator width at the peak (the run's own figure)", lambda S, B, T: f"{S['radiator_width_max_m']:.2f} m"),
         row("The bearing's leads and pumps, mean over the deployment", lambda S, B, T: f"{B['extra_mean_W']:,.0f} W"),
         row("The same, largest", lambda S, B, T: f"{B['extra_max_W']:,.0f} W"),
         row("The same, at the end of the deployment (pumps still running)", lambda S, B, T: f"{B['extra_end_W']:,.0f} W"),
         "", "So a stronger pump takes 0.8 kW off the peak cold load and does not change the days (367.9 against 368.0). The support's leads matter only to the first hundred days' cold load."]
    return L


def main():
    out = next((a[6:] for a in sys.argv if a.startswith("--out=")), SHEET9)
    B = raw()
    r8, _ = parse(SHEET8); r9, lines = parse(RAW)
    ser = s8.series(os.path.join(R9, "ch14", "ch14_thermal_series.csv"))
    t = [x / (30 * 86400.0) for x in ser["t_s"]]; q = ser["q_cryo_cold_W"]; ip = max(range(len(q)), key=lambda i: q[i])
    DP = json.load(open(os.path.join(R9, "design_point.json"), encoding="utf-8"))
    bj = json.load(open(os.path.join(R9, "ch14", "bearing_cold_load.json"), encoding="utf-8"))
    J12 = json.load(open(os.path.join(R9, "slip_optimal", "N165_P12MW", "rising_radius_summary.json")))
    pj = json.load(open(os.path.join(HERE, OUT, "gap_trade_logs", "gap_trade_pump.json"), encoding="utf-8")); j_reg = pj["j_reg"]
    T75 = json.load(open(os.path.join(L9, "round9_table75.json"), encoding="utf-8"))["rows"]; FA = json.load(open(os.path.join(L9, "round9_failures.json"), encoding="utf-8"))

    def at_day(key, day):
        i = min(range(len(B["day"])), key=lambda i: abs(B["day"][i] - day)); return B[key][i]

    def mean(key, d0, d1):
        v = [B[key][i] for i in range(len(B["day"])) if d0 <= B["day"][i] <= d1]; return sum(v) / len(v)
    d_pk = t[ip] * 30.0
    over = {}            # first-cell substring -> new third cell
    over["Table 7.5: the casing's movement"] = (" / ".join(f"{r['exc']*1e3:.1f}" for r in T75) + " mm, with the planned pull of the pair " + " / ".join(f"{r['planned']/1e3:.0f}" for r in T75)
                                               + " kN/m (the support at constant flux; see the report, item 3, for the table that is true)")
    w = FA["worst"][0]
    over["7.7 Tolerances: the smallest gap in any single failure"] = "see the report, item 3: the point and beam lanes' smallest gap is " + next(f"{x['gap']:.1f} mm ({x['label']})" for x in FA["worst"] if not x["label"].startswith("roll")) + "; the roll lanes are open"
    extra_rows = {}      # first-cell substring -> rows to put after it
    parts = [("the fast coils' leads, conduction (ten families; all 20 in a launch window)", "fast_conduction_W", DP["cond"], 2 * DP["cond"]),
             ("the fast coils' regulation (the holding current's Joule heat, the support following the weight)", "regulation_W", j_reg, j_reg / 2),
             ("the support coils' leads while the supply carries the rise", "support_leads_W", 0.0, 0.0), ("the flux pumps while they run", "pumps_W", 0.0, 0.0),
             ("the bearing's leads and pumps together", "extra_W", DP["cond"] + j_reg, 2 * DP["cond"] + j_reg / 2)]
    new14 = []
    for lab, key, a10, a20 in parts:
        new14.append(["Table 14.2 (new row): " + lab, "not in round 8's sheet as a row of its own (round 8: 4.41 kW flat through the deployment; 3.63 and 6.41 kW after it)",
                      f"{at_day(key, d_pk)/1e3:.2f} kW at the peak (day {d_pk:.0f}); {mean(key, 180, 300)/1e3:.2f} kW as the mean of months 6 to 10; {a10/1e3:.2f} kW after the deployment between launches; {a20/1e3:.2f} kW in a launch window",
                      "`bearing_cold_load.py` in the design run (`tools/round9_design_point.py --part=ch14`), `bearing_cold_load.csv`"])
    new14.append(["Table 14.2 (new): the day the support's leads come off", "not in round 8 (no support lead in its cold load)",
                  f"day {bj['leads_off_day']:.0f}, from the run's own current series (the pump moves {bj['pump_rate_A_per_day']:.2f} A a day); the leads bring {bj['lead_heat_mean_W']:,.0f} W per site as a mean while attached ({bj['lead_heat_max_W']:,.0f} W at most); the coil current ends at {bj['I_end_A']:.0f} A",
                  "the same, `bearing_cold_load.json`"])
    # Chapter 11 in full
    tf = os.path.join(R9, "slip_optimal", "N165_P12MW", "rising_radius_deployment_tension.csv"); tf8 = tf.replace("round9_design", "round8_design")

    def ch11(tfile):
        a = s8.tool(["tools/anchor_lines_deployment.py", tfile]); b = s8.tool(["tools/anchor_lines_deployment.py", tfile, "--length=377e3"])
        g = lambda pat, txt=a: (re.search(pat, txt, re.S) or [None, "?"])[1]
        return [("11.1: the pull on a 148 t stowed line, first day; first ten days", f"{g(r'first 1 day.s.: .*?stowed line ([\d.]+) N')} N; {g(r'first 10 day.s.: .*?stowed line ([\d.]+) N')} N"),
                ("11.1: the lean of a loose line on day 1 at 100 m, 1 km and 10 km below", f"{g(r'0.1 km below: .*?lean ([\d.]+) deg')}; {g(r' 1.0 km below: .*?lean ([\d.]+) deg')}; {g(r'10.0 km below: .*?lean ([\d.]+) deg')} degrees"),
                ("11.1: half the final weight per kilogram", f"day {g(r'half on day ([\d.]+)')}"),
                ("11.1: the lowering spool for 377 km: from day 0; after the deployment", f"{g(r'from day 0 .*?brake energy ([\d,]+) kWh', b)} kWh and {g(r'from day 0 .*?peak drum tension ([\d.]+) kN', b)} kN, ending on day {g(r'from day 0 .ends day ([\d.]+)', b)}; {g(r'after the deployment .*?= ([\d.]+) MWh', b)} MWh and {g(r'after the deployment.*?peak drum tension ([\d.]+) kN', b)} kN"),
                ("11.1: the 550 km payout from day 0 ends on; the tenth of the weight is passed on; the payout speed that would finish by then", f"day {g(r'from day 0 .ends day ([\d.]+)')}; day {g(r'the tenth is passed on day ([\d.]+)')}; {g(r'the payout would have to run at ([\d.]+) m/s')} m/s"),
                ("11.1: the tip at 100 km is below 500 m/s through the air from; days of the deployment left", f"day {g(r'that happens on day ([\d.]+)')}; {g(r'on day [\d.]+; ([\d.]+) days of the deployment remain')} days")]
    c9 = ch11(tf); c8 = ch11(tf8)
    S11 = "`tools/anchor_lines_deployment.py` (and `--length=377e3`) on the design run's tension series"
    new11 = [[a[0], b[1] + " (round 8's run through the same tool)", a[1] if a[1] != b[1] else "unchanged", S11] for a, b in zip(c9, c8)]
    over["Table 8.2 'Operations, Levitation'"] = "flux pumps for the support's trims (Paul's ruling of 16 November), a supply through leads for charging and for a change in a hurry; the coil-end drag as in round 8 (229 MW), see the report, item 5b"
    over["14.5: the sled's pulse in the cold"] = "see the report, item 5b"
    over["Table 16.2: 'Levitation plates'"] = "see the report, item 5 (Table 16.2 run again on the present design)"
    over["16.6 (L121)"] = f"{J12['t_end_days']:.1f} days at 165 turns for the design case; the three lighter cases are in the report, item 5"
    L = ["# The numbers sheet for round 9: the reruns on the corrected model (2026-11-16)\n",
         "Asked for in `Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md`, item 2. The report is `Vol_III/round9_reruns_2026-11-16.md`. "
         "The design is unchanged: the book's cable (hardware allowance 10,686 kg/m, of which the plates use 10,057), 165 turns, 12 MW, two LIMs, the `optimal` slip rule, the six-track bearing at the 100 mm gap, "
         "ten of the 20 fast families connected in the deployment and between launches and all 20 for a launch.\n",
         f"**What changed in the model.** (1) The bearing's cold load in the deployment follows the run's own casing weight (`bearing_cold_load.py`): {DP['cond']:,.0f} W of conduction in ten families' fast leads, {j_reg:,.0f} W of regulation heat, "
         f"the support coils' own leads from day 0 until the 0.267 mV pump can follow the current (day {bj['leads_off_day']:.0f} in the design run), and 200 W of pumps after that. Round 8 used a flat 4,411 W. "
         "(2) A fast coil's on-time is the energized length over the sled's speed. (3) The support's pull follows its gap at constant flux in every dynamics model.\n",
         "**How to read the sheet.** 'Book' is round 8's value: the chapters were re-keyed to round 8's sheet on 15 November, and each row here is the same row of that sheet, made by the same code "
         "(`tools/round8_sheet.py`) on round 9's runs. The chapter files were not read again row by row, so where a chapter differs from round 8's sheet the difference is not seen here. "
         "'New' is round 9's value. **'unchanged' means the row's text came out the same to the last printed digit**, and the book session need not check it. Every new number comes from a script run, named in the last column. No book chapter was edited.\n"]
    sec = ""; skip = False; n_same = n_new = 0
    for ln in lines[1:]:
        if ln.startswith("## "):
            sec = re.sub(r"\(.*?\)", "", ln[3:]).strip(); skip = sec.startswith("Section 3b")
            if skip:
                L.append("\n## Section 3b: the rows for the rewrite of Sections 7.5 to 7.9\n")
                L.append("Round 8's Section 3b stands except where this list says otherwise. The new values are in the report (`Vol_III/round9_reruns_2026-11-16.md`).\n")
                L.append("| Item of round 8's Section 3b | In round 9 |\n|---|---|")
                for a, b in (("1. Table 7.3, the worst tape of a support track", "unchanged (the solves are round 8's)"), ("2. Table 7.4, the pull against the ampere-turns", "unchanged"),
                             ("3. Table 7.6, the operating points; the takeover of a parked spacecraft; a loop out with a standing load", "unchanged: none of these uses a pass's on-time (the report, item 3, the lost loop)"),
                             ("4. The capacity sentences", "unchanged"), ("5. The failures the chapter quotes", "**changed**: the report, item 3"),
                             ("6. The stray field", "run again with saturating plates: the report, item 5b"), ("7. The coil ends", "the report, item 5b"),
                             ("8. The bearing's mass on the casing", "unchanged (1,687 kg/m); the 6,000 kg/m casing is in the report, item 5"),
                             ("9. Two terms: the 10 mm protection; the 0.1 mm regulation", "the 0.1 mm regulation's heat is 576 W per site, not 1,352 W, and the support is trimmed continuously, not weekly: the report, items 1 and 2"),
                             ("10. The flux pumps", "**changed**: the design has a flux pump on every support coil for the trims (Paul's ruling of 16 November); the gap trade report, item 4, and this round's item 2"),
                             ("11. The LIMs' lateral guidance", "unchanged: no script gives it"), ("12. Three plain descriptions", "the charging description gains the pump for the trims; otherwise unchanged"),
                             ("13. Chapter 16's lighter casing", "**run again**: the report, item 5")):
                    L.append(f"| {a} | {b} |")
                continue
            L.append("\n" + ln if not ln.startswith("\n") else ln); continue
        if skip:
            continue
        if ln.startswith("|---"):
            L.append(ln); continue
        if ln.startswith("|"):
            c = cells(ln); key = (sec, c[0]); b = r8.get(key)
            if c[0] in ("Where", "Row", "Quantity", "Turns", "Figure") or c[0].startswith("Row ("):
                L.append(fmt_row([x.replace("Book", "Book (round 8)") if i else x for i, x in enumerate(c)])); continue
            if "The fast coils' leads, conduction and the regulation" in c[0]:
                for r in new14:
                    L.append(fmt_row(r))
                n_new += len(new14); continue
            ov = next((v for k, v in over.items() if k in c[0]), None)
            if "8.3 (L106), Table 8.2 'Operations, Cryogenic system'" in c[0]:
                ov = c[2].split(". The site's cold load")[0] + (f". After the deployment the bearing brings {DP['cond']:,.0f} W of conduction in ten families' fast leads and {j_reg:,.0f} W of regulation heat; "
                                                                 "the support's leads are off and the pumps run only for a trim (their loss while they run is not counted here)")
            if "Table 7.5: the casing's movement" in c[0] or "7.7 Tolerances: the smallest gap" in c[0]:
                c[3] = "`bearing_round9.py --part=cases`, `--part=table75`, `--part=failures`"
            if sec.startswith("The figures"):
                L.append(fmt_row([x.replace("figure_work_2026-11-15", "figure_work_2026-11-16").replace("round8_book_figures", "round9_book_figures") for x in c])); continue
            idx = (2, 4) if len(c) == 5 else (2,)
            new = list(c)
            for i in idx:
                old = b[i] if b and len(b) > i else "(not in round 8's sheet)"
                val = ov if (ov is not None and i == 2) else c[i]
                new[i - 1] = old
                if val == old:
                    new[i] = "unchanged"; n_same += 1
                else:
                    new[i] = val; n_new += 1
            new[-1] = new[-1].replace("round8_design_point", "round9_design_point") if len(c) == 4 else new[-1]
            L.append(fmt_row(new))
            if "'The casing reaches ground speed on day" in c[0]:
                for r in new11:
                    L.append(fmt_row(r))
            continue
        L.append(ln)
    L += pump_section()
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"wrote {out}: {len(L)} lines; {n_same} values unchanged, {n_new} new")


if __name__ == "__main__":
    main()
