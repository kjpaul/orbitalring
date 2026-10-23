"""
Extracts the 70 K critical-current angle dependence of the Theva TPL4622-AP tape from the Robinson Research Institute
dataset (Pantoja and Strickland 2024, doi 10.6084/m9.figshare.25893079.v2; sample THV202, measured April and May 2024)
into data/THEVA_Robinson_70K_angle_dependence.csv, which bearing_optimize.py reads when ORBITAL_RING_TAPE_ANGLE=full.
Book: The Orbital Ring, Volume IV: Chapter 6 (Table 6.2) and Chapter 8 (the tape's worst field angle, coil optimization Part 3b).
Run: python tools/extract_theva_robinson_angles.py <folder with "70K angle dependences.xlsx">
     (the folder is not in the repository; the CSV it writes is)
Environment: none.

Columns written: Set field (T), Set angle (deg), Ic/w (A/cm). The angle is the angle between the sample's normal and the
applied field (the dataset's Description.txt): 0 = field through the tape face, 90 = field along it; the sheets run 0 to
240 degrees in 5 degree steps (finer near 120). Needs openpyxl.
"""
import os, sys
import openpyxl

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "data", "THEVA_Robinson_70K_angle_dependence.csv")


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    src = os.path.join(sys.argv[1], "70K angle dependences.xlsx")
    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    rows = []
    for ws in wb.worksheets:
        it = ws.iter_rows(values_only=True); head = next(it)
        iB = head.index("Set field (T)"); iA = head.index("Set angle (°)"); iI = head.index("Ic/w (A/cm)")
        for r in it:
            if r[iI] is None: continue
            rows.append((float(r[iB]), float(r[iA]), float(r[iI])))
    rows.sort()
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("Set field (T),Set angle (deg),Ic/w (A/cm)\n")
        for B, a, i in rows:
            f.write(f"{B:g},{a:g},{i:.4f}\n")
    fields = sorted(set(r[0] for r in rows))
    print(f"wrote {len(rows)} rows, fields {fields} T, to {OUT}")


if __name__ == "__main__":
    main()
