#!/usr/bin/env python3
"""
Figures 3.2, 9.3, 14.1 and 14.2 redrawn from round 8's runs on the book's cable (hardware 10,686 kg/m), by round 7's rules: the installed layouts,
the smallest text at least 8 points at 4.5 inches (Section 4 of Vol_III/Claude_Code_instructions_round8_book_cable_2026-11-15.md).
Book: Orbital Ring Engineering, Volume III: Chapters 3, 9 and 14.
Run: python tools/round8_book_figures.py [--only=3.2,9.3,14] [--out=FOLDER]     (default ../../Vol_III/figure_work_2026-11-15)
It is tools/round7_book_figures.py reading runs/round8_design/ (tools/round8_design_point.py).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os

import round7_book_figures as f7

if __name__ == "__main__":
    f7.R7 = os.path.join(f7.HERE, "runs", "round8_design")
    f7.main()
