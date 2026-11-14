#!/usr/bin/env python3
"""
The numbers sheet of round 5 (2026-11-13), in the form of Vol_III/design_point_numbers_2026-11.md: every quantity the book quotes that the optimum
plates, the 10 ms rise and the rating margin 1.35 move, old value (the book's design point, 10,686 kg/m of hardware) beside new.
Book: Orbital Ring Engineering, Volume III: Chapters 1, 3, 7, 8, 12, 14, 15, 16, 17.
Run: python tools/round5_sheet.py [--out=FILE]
     reads runs/bearing_round5_final.json, runs/bearing_round5_coldload.json, runs/bearing_round5_tapeloss.json, the deployment logs in
     runs/round5_hw/<tag>/lim_design/ and runs/round4_hw/book_4.0/lim_design/, and runs tools/design_point_numbers.py under both hardware masses.
Environment: ORBITAL_RING_OUT (default runs).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json
import math
import os
import re
import subprocess
import sys

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
SHEET = os.path.join(HERE, "..", "..", "Vol_III", "design_point_numbers_round5_2026-11-13.md")
HW_BOOK = 10686.0


def arithmetic(hw):
    env = dict(os.environ, ORBITAL_RING_HW_KGM=f"{hw:.1f}", PYTHONUTF8="1")
    txt = subprocess.run([sys.executable, "tools/design_point_numbers.py"], cwd=HERE, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace").stdout
    rows = {}; block = None
    for line in txt.splitlines():
        if "RETROGRADE" in line: block = "retro"
        elif "PROGRADE" in line: block = "pro"
        m = re.match(r"\s+(.+?)\s{2,}([\d,.\-]+)\s+([\d,.\-]+)\s*(\S*)\s*$", line)
        if m and block == "retro":
            rows[m.group(1).strip()] = (m.group(3), m.group(4))
    return rows


def deployment(path):
    if not os.path.exists(path):
        return None
    txt = open(path, encoding="utf-8", errors="replace").read()
    def grab(pat, cast=float):
        m = re.search(pat, txt); return cast(m.group(1)) if m else None
    return dict(days=grab(r"End of the run at ([\d.]+) days"), p_cryo=grab(r"p_cryo_max\s+([\d.]+) MW"), q_cold=grab(r"q_cryo_cold_max\s+([\d.]+) kW"),
                EJ=grab(r"Total energy delivered: ([\d.]+) EJ"))


def main():
    out = SHEET
    for a in sys.argv[1:]:
        if a.startswith("--out="): out = a[6:]
    fin = json.load(open(os.path.join(OUT, "bearing_round5_final.json"), encoding="utf-8"))
    cold = json.load(open(os.path.join(OUT, "bearing_round5_coldload.json"), encoding="utf-8"))
    tl = json.load(open(os.path.join(OUT, "bearing_round5_tapeloss.json"), encoding="utf-8"))
    rs = fin["results"]["support:opt135"]; rl = fin["results"]["fast:opt135"]
    hw = fin["hw"]
    A_book = arithmetic(HW_BOOK); A_new = arithmetic(hw)
    dep_book = deployment(os.path.join(OUT, "round4_hw", "book_4.0", "lim_design", "rising_radius_lim_simulation.txt"))
    dep = {tag: deployment(os.path.join(OUT, "round5_hw", tag, "lim_design", "rising_radius_lim_simulation.txt")) for tag in ("opt_0", "opt_10", "opt_14")}
    L = []
    L.append("# The numbers sheet for round 5 of the bearing (2026-11-13)\n")
    L.append(f"Asked for in Section 9 of `Vol_III/Claude_Code_instructions_bearing_round5_2026-11-13.md`. The design point: the ring of `ring_altitude.py` (700 km, the rising quadratic, the casing held at 12,000 kg/m) with the optimum plates of round 5 in place of the book's four 1.4 m by 200 mm plates: two support plates of {rs['Wp']:.2f} m by {rs['T']*1e3:.0f} mm with {rs['Wl']:.2f} m loops, four lift and load plates of {rl['Wp']:.2f} m by {rl['T']*1e3:.0f} mm with a {rl['Wl']:.2f} m loop, the cable hardware {hw:,.0f} kg/m against 10,686; the fast coils' rise 10 ms; the rating margin 1.35; the fast coils of 25 m in 20 families with 10 or 14 series turns. 'Book' is the chapter file's present value where the chapter quotes one, else the run of record (the 10,686 kg/m hardware); 'Round 5' is the run under the new hardware. Every number comes from a script run (`tools/design_point_numbers.py` under `ORBITAL_RING_HW_KGM`, `lim_simulation.py --rising-radius` through `tools/round5_deployment.py`, `bearing_round5.py`), or from the round 5 report's own tables, as the script column says. The book chapters are not edited.\n")
    L.append("## Chapters 1, 3, 12, 16 and 17: the cable under the new hardware\n")
    L.append("| Where | Book | Round 5 | Script |\n|---|---|---|---|")
    keys = [("Chapter 1 Table 1.1, Chapter 3 Table 3.1: structural cable", "Table 1.1, Table 3.1: structural cable"), ("cable with hardware", "Table 1.1, Table 3.1: cable with hardware"),
            ("ring (cable, hardware, casing)", "Table 3.1: ring mass per meter"), ("cable speed", "Table 1.1, Table 3.1: cable velocity"), ("speed past the casing", "Table 1.1, Table 3.1: relative velocity"),
            ("tension", "Table 1.1, Table 3.1, Equation 3.8: tension"), ("cross-section", "Table 1.1, Table 3.1: cross-section"), ("side", "Table 1.1, Figure 1.1, Table 3.1: side"),
            ("cable mass, structural, around the ring", "Chapter 3: structural cable around the ring"), ("cable mass with hardware, around the ring", "Chapter 3: cable with hardware around the ring"),
            ("ring mass", "Table 3.1: ring mass"), ("deployment energy (closed form)", "Chapters 3 and 9: deployment energy, closed form"),
            ("mass ratio, structural cable per kg of casing (Table 3.1)", "Table 3.1: mass ratio structural / casing"), ("mass ratio, cable with hardware per kg of casing (Ch 12, 16)", "Chapters 12 and 16: cable with hardware per kg of casing"),
            ("Chapter 3: hoop stress with no casing load, (m_c/A)(v_c^2 - v_orb^2)", "3.2: hoop stress with no casing load"), ("added load at which the hoop stress reaches zero", "3.2: added load at zero hoop stress"),
            ("second moment of area, side^4/12 (Eq 3.9)", "Equation 3.9: second moment"), ("static bending length sqrt(EI/T) (Eq 3.10)", "Equation 3.10: bending length"),
            ("cable centrifugal acceleration v_c^2/r (Eq 4.3)", "Equation 4.3: centrifugal acceleration"), ("net outward acceleration (Eq 4.4)", "Equation 4.4: net outward acceleration"),
            ("Chapter 12: cable claimed by the 1,200 kg/m array, at the with-hardware ratio", "Chapter 12: cable claimed by the array"), ("the same around the ring", "Chapter 12: the same around the ring"),
            ("Chapter 17: the Earth ring's tension quoted in the scaling", "Chapter 17: the Earth ring's tension"), ("the Earth ring's cable with hardware quoted in the scaling", "Chapter 17: the Earth ring's cable with hardware")]
    for k, where in keys:
        if k in A_book and k in A_new:
            b, u = A_book[k]; n, _ = A_new[k]
            L.append(f"| {where} | {b} {u} | {n} {u} | `tools/design_point_numbers.py` (arithmetic on `ring_altitude.py`) |")
    L.append("\n## Chapter 7: the bearing\n")
    L.append("| Where | Book | Round 5 | Script |\n|---|---|---|---|")
    tape_load = fin["tape_load"]; tape_casing = fin["tape_casing"]
    c10 = cold["leads"]["10"]; c14 = cold["leads"]["14"]
    rows7 = [
        ("7.3, the plates ('each plate is 1.4 m wide ... and is 200 mm thick'): support", "1.4 m x 200 mm, 2,204 kg/m each", f"{rs['Wp']:.2f} m x {rs['T']*1e3:.0f} mm, {rs['plate_m']:,.0f} kg/m each", "`bearing_round5.py --part=final`"),
        ("7.3 and 7.7, the plates: lift and load", "1.4 m x 200 mm (the book's four plates; six plates since Stage B)", f"{rl['Wp']:.2f} m x {rl['T']*1e3:.0f} mm, {rl['plate_m']:,.0f} kg/m each, four", "the same"),
        ("7.3, the support loops", "two of 0.5 m", f"two of {rs['Wl']:.2f} m", "the same"),
        ("7.7, the lift and load loops", "one of 0.8 m (Part 3b) / 1.0 m (the book's lift track)", f"one of {rl['Wl']:.2f} m", "the same"),
        ("Table 7.1, Chapters 1 and 3: cable hardware", "10,686 kg/m", f"{hw:,.0f} kg/m", "the same"),
        ("Equation 7.2, Tables 7.3, 7.4 and 7.6: the support winding's cap, ampere-turns per leg", "322,000", f"{rs['NI_req']:,.0f}", "the same"),
        ("Table 7.3: the worst 12 mm tape's limit (117 A at half)", "234 A at 1.64 T, 54 deg", f"{rs['ic_req']:.0f} A at {rs['B']:.2f} T, {rs['th']:.0f} deg", "the same"),
        ("Equation 7.2, 7.3, Table 7.6: tapes per leg; tape per track", "2,753; 94 kg/m", f"{rs['turns_leg']:,.0f}; {rs['tape']:.1f} kg/m", "the same"),
        ("support track: rating (70 percent) and cap (100 percent)", "336 and 476 kN/m (round 4 on the book's plate)", f"{rs['F70']/1e3:.0f} and {rs['F100']/1e3:.0f} kN/m ({rs['margin']:.2f} x the requirement)", "the same"),
        ("support plate: saturation and peak field at the requirement", "12.6 percent, 2.10 T", f"{rs['sat']*100:.1f} percent, {rs['B_plate']:.2f} T", "the same"),
        ("7.7: the lift track's ampere-turns per leg at 196.5 kN/m", "421,000", f"{rl['NI_req']:,.0f}", "the same"),
        ("7.7: the worst 3 mm tape's limit (30 A at half)", "60 A at 1.51 T, 60 deg", f"{rl['ic_req']:.0f} A at {rl['B']:.2f} T, {rl['th']:.0f} deg", "the same"),
        ("the lift track's tapes per leg; tape per track (Part 3b)", "14,002; 59.7 kg/m", f"{rl['turns_leg']:,.0f}; {rl['tape']:.1f} kg/m", "the same"),
        ("lift track: rating and cap", "288 and 415 kN/m", f"{rl['F70']/1e3:.0f} and {rl['F100']/1e3:.0f} kN/m ({rl['margin']:.2f} x)", "the same"),
        ("lift plate: saturation and peak field at the requirement", "0.1 percent, 1.53 T", f"{rl['sat']*100:.1f} percent, {rl['B_plate']:.2f} T", "the same"),
        ("Tables 7.4 and 7.6: the load winding (12,197 turns of 3 mm on the support loops in the book; the load track of the lift design since Stage B): tape per track", "60 kg/m (Stage B)", f"{tape_load:.1f} kg/m", "the same"),
        ("tape per meter of casing, six tracks; the bearing's casing side (x 2)", "424 kg/m; 933 kg/m (Stage B 6.4, four tracks)", f"{tape_casing:.0f} kg/m; {2*tape_casing:,.0f} kg/m", "the same"),
        ("tape around the ring", "", f"{tape_casing*44473.246e3/1e9:.2f} Mt", "the same"),
        ("the load track's coupling to the support track beside it", "5.2 percent of the support's force at 0.5 m (Stage B 3; its 0.9 percent was the flux)", f"{fin['coupling']['dF']*100:.2f} percent of the support's normal force at {fin['coupling']['d_edge']:.1f} m ({abs(fin['coupling']['dF'])*47.5:.1f} kN/m per track)", "the same (corrected 2026-11-14)"),
        ("7.7: the fast coils' rise", "2 ms (the protection)", "10 ms (Paul's ruling)", "the brief"),
        ("7.7: a 25 m coil's stored energy and supply rating", "0.15 MJ, 146 MVA at 2 ms", f"{rl['W_mag']*25/20/1e6:.3f} MJ, {2*rl['W_mag']*25/20/10e-3/1e6:.0f} MVA at 10 ms", "`bearing_round5.py --part=switches`"),
        ("7.7: the series turns and the terminal voltage", "'two or three turns' (Stage B; withdrawn in round 4)", f"10 turns of {rl['NI_req']/20/10/(0.5*rl['ic_req']):.0f} layers at {rl['NI_req']/20/10:,.0f} A and {10*2*rl['W_mag']/rl['NI_req']*25/10e-3/1e3:.1f} kV, or 14 of {rl['NI_req']/20/14/(0.5*rl['ic_req']):.0f} at {rl['NI_req']/20/14:,.0f} A and {14*2*rl['W_mag']/rl['NI_req']*25/10e-3/1e3:.1f} kV", "the same"),
        ("7.7: the plate heat of a launch", "238 GJ (the book), 231 (Stage B P4 on the book's plates)", f"{fin['plate_heat']['E_launch']/1e9:.0f} GJ on the optimum plates", "`bearing_round5.py --part=final`"),
        ("7.7: the tape's pulse heat per launch (the 605 GJ bound and the 53 GJ 5 degree case)", "53 to 605 GJ", f"{sum(tl['launch']['10']['v'])/1e9:.0f} GJ with the tapes' faces vertical ({tl['launch']['10']['v'][0]/1e9:.0f} through the face, {tl['launch']['10']['v'][1]/1e9:.0f} the soldered slab) at 10 turns; {sum(tl['launch']['14']['v'])/1e9:.0f} at 14", "`bearing_round5.py --part=tapeloss`"),
        ("Table 7.5: the casing's movement at 1 / 3 / 5 / 10 / 20 / 30 km/s", "7.6 / 2.6 / 1.3 / 1.7 / 6.0 / 9.7 mm", "12.6 / 11.0 / 10.2 / 10.2 / 10.8 / 13.8 mm (round 4's 10 ms run, the 10 mm protection in every lane)", "`bearing_round4.py --part=casing`"),
        ("the smallest single-failure gap, 10 ms", "", "24.3 mm (a pull 50 ms late)", "the same"),
        ("the e-folding time with the support frozen", "", "135 ms", "the same"),
    ]
    for where, b, n, s in rows7:
        L.append(f"| {where} | {b} | {n} | {s} |")
    L.append("\n## Chapter 8: the stores and the switches\n")
    L.append("| Where | Book | Round 5 | Script |\n|---|---|---|---|")
    L.append(f"| 8.1: the store for the sled's passage | SMES named; 'the storage on the ring is there for the fast transients' | film capacitors, warm, at every coil (0.155 MJ, 330 kg/m) or shared by a site's four tracks (half, 165 kg/m); SMES excluded on its cold loss | `bearing_round5.py --part=stores` |")
    L.append(f"| 8.1: the pulse power of the pulled zone at 30 km/s | | 3.7 GW per lift track at the front, 7.4 for the pair, returned 2.25 km behind; 619 times a site's 12 MW for 17 ms (the log `bearing_round5_all.txt`, S3.1) | the same |")
    L.append(f"| the section bus if it carried the pulses (way b) | 11,600 A at 100 kV, 18 tapes each way | 74 kA, 114 tapes each way, 6.4 times the bus's rating | the same |")
    L.append(f"| the switches per coil | | 10 turns: 5 StakPaks in series per position, 86 kg per coil, 275 kg/m; 14 turns: 7, 120 kg, 385 kg/m (the log, S4) | `bearing_round5.py --part=switches` |")
    L.append("\n## Chapter 14: the cold loop\n")
    L.append("| Where | Book | Round 5 | Script |\n|---|---|---|---|")
    L.append(f"| Table 14.2, after the deployment: cold-side load per site | 14.5 kW | {14.5 + c10['idle_site']/1e3:.1f} kW at 10 turns, {14.5 + c14['idle_site']/1e3:.1f} at 14 (the fast coils' leads idle) | `bearing_round5.py --part=coldload` |")
    L.append(f"| in operation at 500 launches a day, per site | | {cold['tot10']*500/1e3:.1f} kW at 10 turns, {cold['tot14']*500/1e3:.1f} at 14 (leads, tape; faces vertical) | the same |")
    L.append(f"| Table 14.2, compressor power after the deployment | 1.37 MW | {(14.5 + c10['idle_site']/1e3)*94.3/1e3:.2f} MW idle at 10 turns; {cold['tot10']*500*94.3/1e6:.1f} MW in operation | the same, x 94.3 |")
    L.append("\n## Chapter 15 and Table 1.1: the deployment\n")
    L.append("| Where | Book | Round 5 | Script |\n|---|---|---|---|")
    def dd(d, k, f):
        return f.format(d[k]) if d and d.get(k) is not None else "(not run)"
    L.append(f"| days at 12 MW, 150 turns | {dd(dep_book, 'days', '{:.1f}')} | {dd(dep['opt_0'], 'days', '{:.1f}')} (plates alone, the loops' 4.0 m of cold width), {dd(dep['opt_10'], 'days', '{:.1f}')} with the leads at 10 turns, {dd(dep['opt_14'], 'days', '{:.1f}')} at 14 | `tools/round5_deployment.py` |")
    L.append(f"| cryogenic plant, peak | {dd(dep_book, 'p_cryo', '{:.2f} MW')} | {dd(dep['opt_0'], 'p_cryo', '{:.2f}')} / {dd(dep['opt_10'], 'p_cryo', '{:.2f}')} / {dd(dep['opt_14'], 'p_cryo', '{:.2f} MW')} | the same |")
    L.append(f"| cold load, peak | {dd(dep_book, 'q_cold', '{:.1f} kW')} | {dd(dep['opt_0'], 'q_cold', '{:.1f}')} / {dd(dep['opt_10'], 'q_cold', '{:.1f}')} / {dd(dep['opt_14'], 'q_cold', '{:.1f} kW')} | the same |")
    L.append(f"| LIM work | {dd(dep_book, 'EJ', '{:.3f} EJ')} | {dd(dep['opt_0'], 'EJ', '{:.3f}')} / {dd(dep['opt_10'], 'EJ', '{:.3f}')} / {dd(dep['opt_14'], 'EJ', '{:.3f} EJ')} | the same |")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {out}: {len(L)} lines")


if __name__ == "__main__":
    main()
