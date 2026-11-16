#!/usr/bin/env python3
"""
Figures 3.2, 9.3, 14.1 and 14.2 redrawn from round 9's runs (the bearing's cold load from the run's own weight series), by round 7's rules.
Book: Orbital Ring Engineering, Volume III: Chapters 3, 9 and 14.
Run: python tools/round9_book_figures.py [--only=3.2,9.3,14] [--out=FOLDER]     (default ../../Vol_III/figure_work_2026-11-16)
It is tools/round7_book_figures.py reading runs/round9_design/ (tools/round9_design_point.py), with one change: Figure 14.2's radiator width is
computed as Table 14.2 computes it (the cold load and the compressor's work at the book's COP of 0.0106, rejected at 1,306 W/m2 over 450 m), so
that the label and the table agree to the digit.  Round 8's figure took the simulation's own width column, which differed in the last digit.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import round7_book_figures as f7

WIDTH_OLD = 'w.append(float(row["radiator_width_m"]))'
WIDTH_NEW = 'w.append(float(row["q_cryo_cold_W"]) * (1 + 1 / 0.0106) / 1306.0 / 450.0)'

if __name__ == "__main__":
    f7.R7 = os.path.join(f7.HERE, "runs", "round9_design")
    _lifted = f7.lifted

    def lifted(src, *a):
        if "fig14" in src or "radiator_width_m" in src:
            assert WIDTH_OLD in src
            src = src.replace(WIDTH_OLD, WIDTH_NEW)
        return _lifted(src, *a)
    f7.lifted = lifted
    if not any(a.startswith("--out=") for a in sys.argv):
        sys.argv.append("--out=" + os.path.join(f7.VOL, "figure_work_2026-11-16"))
    f7.main()
