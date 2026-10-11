#!/usr/bin/env python3
"""
List every place the Volume III chapter files quote one of the old ring numbers (cable hardware 4,076 kg/m), grouped by quantity,
as chapter, line and sentence fragment, for the renumbering pass after the change to 10,686 kg/m (2026-10-12).
Searches Vol_III/book/Vol_III_Ch*.md, not the _Paul or .bak files.  Hits that could be a different quantity with the same digits
are flagged.  Usage: python book_value_hits.py > log
"""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
BOOK = os.path.join(HERE, "..", "..", "Vol_III", "book")

# quantity -> (old forms as they appear in prose, new value text)
QUANTITIES = [
    ("cable with hardware 119,592 kg/m", [r"119,592", r"119592", r"119\.6\s*(thousand|t)"], "125,691 kg/m"),
    ("ring per meter 131,592 kg/m", [r"131,592", r"131592"], "137,691 kg/m"),
    ("cable structural 115,516 kg/m", [r"115,516", r"115516"], "115,005 kg/m"),
    ("hardware 4,076 kg/m", [r"4,076", r"4076"], "10,686 kg/m"),
    ("levitation plates 2,204 kg/m, two 100 mm plates", [r"2,204", r"2204", r"\b10 cm thick", r"100 mm thick", r"0\.10 × 1\.4"], "8,814 kg/m, four 200 mm plates"),
    ("cable speed 8,309 m/s", [r"8,309", r"8309"], "8,270 m/s"),
    ("relative speed 8,825 m/s", [r"8,825", r"8825"], "8,786 m/s"),
    ("tension 849.4 GN", [r"849\.4", r"\b849 GN", r"\b849\b"], "845.6 GN"),
    ("cable area 67.95 m2", [r"67\.95", r"\b67\.9\b", r"\b68 m²", r"\b68 m\^2"], "67.65 m2"),
    ("cable side 8.24 m", [r"8\.24 m", r"8\.243", r"8\.24\b"], "8.22 m (8.225)"),
    ("wave speed 2,712 / 2,665 m/s", [r"2,712", r"2712", r"2,665", r"2665"], "2,712 / 2,594 m/s"),
    ("cable net outward acceleration 1.798 m/s2", [r"1\.798", r"1\.80 m/s", r"1\.8 m/s²"], "1.706 m/s2"),
    ("cable centrifugal acceleration 9.754 m/s2", [r"9\.754", r"9\.75 m/s"], "9.663 m/s2"),
    ("deployment demand 711.6 GW", [r"711\.6", r"711,6", r"\b712 GW", r"711 GW"], "711.6 GW (unchanged: sites x 8 MW)"),
    ("deployment time 567 days", [r"\b567", r"1\.55 years"], "from the rerun lim_simulation"),
    ("deployment energy 18.9 EJ", [r"18\.9 EJ", r"18\.89", r"18\.887"], "18.804 EJ"),
    ("LIM supply frequency 44 Hz at 8,825 m/s", [r"44\.1 Hz", r"44 Hz", r"44\.35", r"44\.13"], "43.9 Hz at 8,786 m/s, pole pitch 100 m"),
    ("cable mass ratio 9.97", [r"9\.97", r"9\.63"], "10.47"),
    ("structural cable 5,137 Mt, with hardware 5,319 Mt, ring 5,852 Mt", [r"5,137", r"5,319", r"5,852", r"5\.85 Gt", r"5\.32 Gt", r"5\.14 Gt"], "5,115 / 5,590 / 6,124 Mt"),
    ("net momentum p 987.5 MN s/m", [r"987\.5", r"987 MN"], "from the rerun hand numbers"),
    ("Table 2.1 lateral forces 12.6 / 25.0 / 49.3 / 72.0 kN/m", [r"12\.6 kN", r"25\.0 kN", r"49\.3 kN", r"72\.0 kN", r"\b72 kN"], "from the rerun hand numbers"),
    ("casing weight 95.0 kN/m", [r"95\.0 kN", r"95 kN/m", r"95\.02"], "95.0 kN/m (unchanged)"),
    ("EML capacity 218 / 871 / 573 kPa / 776 kN/m (old sheet-formula bearing)", [r"\b218 kN", r"\b871", r"\b573 k", r"\b776"], "replaced by the air-core design values (final_bearing_2026-10-12.md)"),
    ("second moment 385 m4, L_char 15.0 m", [r"\b385 m", r"\b15\.0 m\b", r"15\.0 m\b"], "from the rerun hand numbers / equation log"),
    ("cable net outward 1.798 in the gap-closing time 0.143 s", [r"0\.143", r"1\.40 m/s"], "from the rerun hand numbers"),
    ("Halbach / limiter plate masses (molybdenum, copper) on the cable", [r"molybdenum", r"OFHC"], "limiters removed (2026-10-11)"),
]

files = sorted(f for f in glob.glob(os.path.join(BOOK, "Vol_III_Ch*.md")) if "_Paul" not in f and not f.endswith(".bak"))
print(f"Chapter files searched ({len(files)}): " + ", ".join(os.path.basename(f) for f in files))
print("Each hit: chapter, line, fragment.  '(check: same digits, maybe another quantity)' marks forms that can be something else.")
total = 0
for q, pats, new in QUANTITIES:
    hits = []
    for f in files:
        with open(f, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh, 1):
                for p in pats:
                    for m in re.finditer(p, line):
                        a = max(0, m.start() - 60); b = min(len(line), m.end() + 60)
                        frag = line[a:b].strip().replace("|", " ")
                        hits.append((os.path.basename(f).replace("Vol_III_", "").replace(".md", ""), i, p, frag))
    print()
    print(f"## {q}  ->  {new}   ({len(hits)} hits)")
    seen = set()
    for ch, i, p, frag in hits:
        key = (ch, i, frag)
        if key in seen:
            continue
        seen.add(key)
        flag = "  (check: same digits, maybe another quantity)" if (p.startswith(r"\b") or p in (r"8\.24\b", r"\b849\b", r"\b567", r"\b385 m", r"15\.0 m\b", r"\b72 kN", r"\b871", r"\b776", r"\b573 k", r"\b218 kN")) else ""
        print(f"  {ch} line {i}: ...{frag}...{flag}")
        total += 1
print()
print(f"{total} hits in all.")
