#!/usr/bin/env python3
"""
The book figures that carry round 7's design point, redrawn from the new runs with their present layouts and the text raised to at least 8 points
at the width a figure prints (4.5 inches): Figures 3.2, 9.3, 14.1 and 14.2 (Section 3 of Vol_III/Claude_Code_instructions_bearing_round7_2026-11-15.md).
Book: Orbital Ring Engineering, Volume III: Chapters 3, 9 and 14.
Run: python tools/round7_book_figures.py [--only=3.2,9.3,14] [--out=FOLDER]     (default ../../Vol_III/figure_work_2026-11-15)
It runs the layout scripts the installed figures were made with (Vol_III/figure_scripts_2026-11-07/gen_fig03_02_rising.py, tools/fig14_figures.py,
tools/rising_radius_design_165.figure93) with two changes only: the figure's size, and every font size raised to the least that prints at 8 points.
Reads runs/round7_design/ (tools/round7_design_point.py, tools/ch14_thermal_series.py).  Prints each figure's smallest text at 4.5 inches.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VOL = os.path.join(HERE, "..", "..", "Vol_III")
R7 = os.path.join(HERE, "runs", "round7_design")
PRINT_W = 4.5; MIN_PT = 8.0


def lifted(src, fig_w, figsize_old, figsize_new):
    """The script's source with its figure size changed and every font size raised to what prints at MIN_PT."""
    need = MIN_PT * fig_w / PRINT_W + 0.05
    src = src.replace(f"figsize={figsize_old}", f"figsize={figsize_new}")
    src = re.sub(r"fontsize=([\d.]+)", lambda m: f"fontsize={max(float(m.group(1)), need):.2f}", src)
    src = re.sub(r'"font.size": ([\d.]+)', lambda m: f'"font.size": {max(float(m.group(1)), need):.2f}', src)
    return src


def smallest(path_note, figs):
    for name, fig in figs:
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
        s = min(t.get_fontsize() for t in texts); w = fig.get_figwidth()
        bbs = [t.get_window_extent(r) for t in texts]; fb = fig.bbox
        out = [t.get_text()[:30] for t, b in zip(texts, bbs) if b.x0 < fb.x0 - 1 or b.x1 > fb.x1 + 1 or b.y0 < fb.y0 - 1 or b.y1 > fb.y1 + 1]
        print(f"  {name}: {w:g} inches wide as drawn; smallest text {s:.1f} pt, {s * PRINT_W / w:.2f} pt at {PRINT_W} inches; texts outside the figure: {out or 'none'}")


def run_script(path, argv, fig_w, figsize_old, figsize_new, names):
    src = lifted(open(path, encoding="utf-8").read(), fig_w, figsize_old, figsize_new)
    keep = []
    _close = plt.close
    plt.close = lambda *a, **k: None                                  # keep the figures for the text check
    old_argv = sys.argv
    try:
        sys.argv = argv
        g = {"__name__": "__main__", "__file__": path}
        exec(compile(src, path, "exec"), g)
    finally:
        sys.argv = old_argv; plt.close = _close
    figs = [plt.figure(n) for n in plt.get_fignums()]
    smallest(path, list(zip(names, figs)))
    for f in figs:
        plt.close(f)


def main():
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), "3.2,9.3,14").split(",")
    out = os.path.abspath(next((a[6:] for a in sys.argv if a.startswith("--out=")), os.path.join(VOL, "figure_work_2026-11-15")))
    os.makedirs(out, exist_ok=True)
    if "14" in only:
        series = os.path.join(R7, "ch14", "ch14_thermal_series.csv")
        run_script(os.path.join(HERE, "tools", "fig14_figures.py"), ["fig14_figures.py", series, out], 6.2, "(10, 5)", "(6.2, 3.7)", ["Figure 14.1", "Figure 14.2"])
    if "3.2" in only:
        retro = os.path.join(R7, "slip_optimal", "N165_P12MW", "rising_radius_deployment_tension.csv")
        pro = os.path.join(R7, "prograde_N165_P12MW", "rising_radius_deployment_tension.csv")
        argv = ["gen_fig03_02_rising.py", retro] + ([pro] if os.path.exists(pro) else []) + [f"--out={out}"]
        run_script(os.path.join(VOL, "figure_scripts_2026-11-07", "gen_fig03_02_rising.py"), argv, 6.8, "(6.8, 9.6)", "(6.8, 11.4)", ["Figure 3.2"])
    if "9.3" in only:
        import rising_radius_design_165 as d
        import inspect
        d.BASE = R7; d.TURNS_REF = os.path.join(R7, "slip_optimal")
        src = inspect.getsource(d.figure93)
        src = lifted(src, 6.0, "(6.5, 4.4)", "(6.0, 4.4)").replace('out = os.path.join(BASE, f"fig09_03_turns_{rule}")', f'out = os.path.join(r"{out}", "fig09_03_turns_sweep")').replace("plt.close(fig)", "")
        g = dict(vars(d)); exec(src, g); g["figure93"]("optimal")
        smallest("", [("Figure 9.3", plt.figure(plt.get_fignums()[-1]))])
    print("  written to", out)


if __name__ == "__main__":
    main()
