#!/usr/bin/env python3
"""
Two figures of the gap trade round (Vol_III/Claude_Code_instructions_gap_trade_2026-11-16.md).
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run (from Python/orbitalring): python tools/gap_trade_figures.py [output folder]     (default ../../Vol_III/figure_work_2026-11-16)
  fig07_candidate_field_at_plate   item 3, a candidate for the book: the field across the face of a support plate (normal operation and the cap) and of a
                                   compensation plate (196.5 kN/m), in tesla against position, the bundles marked under the axis to scale, the even field
                                   with the same pull as a horizontal line.  From runs/gap_trade_logs/gap_trade_face.json (bearing_gap_trade.py --part=solve-face).
  fig07_candidate_2_two_coils      item 5, for Volume III: the two kinds of coil seen from above, the support track's two loops separate, the coils of one
                                   loop all the same width, one behind the other and offset along the track.
4.5 inches wide, the width it prints at, so a font size here is the size in print; the smallest text is 8 points or more.  SVG and PNG.  The font,
colors and line weights are those of tools/round6b_figures.py (the book's Figures 7.1 and 7.3).  Each PNG is to be opened and looked at.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

import round6b_figures as f6

C = f6.C; W_IN = f6.W_IN; LW, LW2 = f6.LW, f6.LW2
FS = 8.5
LOGS = os.path.join(os.environ.get("ORBITAL_RING_OUT", "runs"), "gap_trade_logs")


def check(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    bbs = [(t, t.get_window_extent(r)) for t in texts]; fb = fig.bbox
    out = [t.get_text() for t, b in bbs if b.x0 < fb.x0 - 0.5 or b.x1 > fb.x1 + 0.5 or b.y0 < fb.y0 - 0.5 or b.y1 > fb.y1 + 0.5]
    over = [(a.get_text(), c.get_text()) for i, (a, ba) in enumerate(bbs) for c, bc in bbs[i + 1:] if ba.x0 < bc.x1 and bc.x0 < ba.x1 and ba.y0 < bc.y1 and bc.y0 < ba.y1]
    small = min(t.get_fontsize() for t in texts)
    print(f"  {name}: {len(texts)} texts; smallest {small:g} pt at {W_IN} inches wide; outside the figure: {out or 'none'}; overlapping: {over or 'none'}; {fig.get_figheight():.2f} inches tall")
    return small, out, over


def save(fig, folder, name):
    res = check(fig, name)
    fig.savefig(os.path.join(folder, name + ".png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(folder, name + ".svg"), facecolor="white")
    plt.close(fig)
    return res


# ------------------------------------------------------------------------------------------------ item 3: the field across the face of the plate
def fig_field(folder):
    d = json.load(open(os.path.join(LOGS, "gap_trade_face.json"), encoding="utf-8"))
    fig, axes = plt.subplots(2, 1, figsize=(W_IN, 5.0))
    fig.subplots_adjust(left=0.125, right=0.975, top=0.985, bottom=0.085, hspace=0.30)
    panels = [(axes[0], "Support plate", [("support cap", C["support"], LW2, "the cap"), ("support normal", C["support"], LW2 * 0.6, "normal operation")], 1.45),
              (axes[1], "Compensation plate", [("fast 196.5", C["load"], LW2, "196.5 kN/m")], 1.45)]
    for ax, title, curves, ymax in panels:
        s0 = d[curves[0][0]]; Wp = s0["Wp"]
        yb = -0.20 * ymax / 1.45                                        # the band under the axis where the bundles are drawn
        for key, col, lw, lab in curves:
            s = d[key]; x = np.array(s["line"]["x"]); B = np.hypot(np.array(s["line"]["By"]), np.array(s["line"]["Bx"]))
            ax.plot(x * 1e3, B, color=col, lw=lw, zorder=4, solid_capstyle="round")
            ev = s["line"]["even"]
            ax.plot([-Wp / 2 * 1e3, Wp / 2 * 1e3], [ev, ev], color=col, lw=LW, ls=(0, (4, 2.5)), zorder=3)
            s["_even"] = ev
        ax.axhline(0, color=C["ink"], lw=LW, zorder=2)
        # the plate's width and the bundles under the axis, to the scale of the x axis (the heights to the same scale in millimeters)
        mm_per_T = 1.0                                                  # set below from the axes' aspect so that the bundles are drawn to scale
        ax.set_xlim(-Wp / 2 * 1e3 - 30, Wp / 2 * 1e3 + 255); ax.set_ylim(yb * 1.25, ymax)
        fig.canvas.draw()
        bb = ax.get_window_extent(); x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        T_per_mm = (y1 - y0) / bb.height * bb.width / (x1 - x0)          # tesla of the y axis that one millimeter of the x axis spans on the page
        h_b = s0["bh"] * 1e3 * T_per_mm
        for xl, sg in s0["legs"]:
            ax.add_patch(Rectangle(((xl - s0["bw"] / 2) * 1e3, -0.035 - h_b), s0["bw"] * 1e3, h_b, fc=curves[0][1], ec="none", zorder=4))
        ax.set_ylabel("Field at the face (T)", fontsize=FS, color=C["ink"])
        ax.set_yticks([0, 0.5, 1.0])
        ax.tick_params(labelsize=FS, colors=C["ink"], width=LW, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(LW); ax.spines[sp].set_color(C["ink"])
        ax.spines["left"].set_bounds(0, ymax); ax.spines["bottom"].set_bounds(-Wp / 2 * 1e3, Wp / 2 * 1e3)
        ax.set_xticks([t for t in range(-600, 601, 200) if abs(t) <= Wp / 2 * 1e3 + 1])
        xr = Wp / 2 * 1e3 + 18
        for key, col, lw, lab in curves:
            ax.text(xr, d[key]["_even"], f"even field\n{d[key]['_even']:.2f} T", fontsize=FS, color=col, ha="left", va="center", linespacing=1.05)
        ax.text(xr, -0.035 - h_b / 2, "bundles", fontsize=FS, color=curves[0][1], ha="left", va="center")
        ax.text(x0 + 8, ymax * 0.995, title, fontsize=FS, fontweight="bold", color=C["ink"], ha="left", va="top")
    axes[1].set_xlabel("Position across the plate (mm)", fontsize=FS, color=C["ink"], x=0.40)
    # the labels: each curve and its even field, placed where the curve leaves room
    a0, a1 = axes
    sc, sn, sf = d["support cap"], d["support normal"], d["fast 196.5"]
    a0.text(-300, 1.17, "the cap", fontsize=FS, color=C["support"], ha="center", va="bottom")
    a0.text(300, 0.475, "normal\noperation", fontsize=FS, color=C["support"], ha="center", va="bottom", linespacing=1.05)
    a1.text(0, 0.95, "196.5 kN/m", fontsize=FS, color=C["load"], ha="center", va="bottom")
    return save(fig, folder, "fig07_candidate_field_at_plate")


# ------------------------------------------------------------------------------------------------ item 5: the two kinds of coil
def fig_two_coils(folder):
    H = 3.5
    fig, ax = f6.new_fig(H, (0, 100), (0, 100 * H / W_IN))
    def arrow(x, y, dr, col):
        ax.annotate("", xy=(x + 2.4 * dr, y), xytext=(x - 2.4 * dr, y), arrowprops=dict(arrowstyle="-|>", lw=0, color=col, mutation_scale=9), zorder=6)
    def coils(n, L, st, x_first, y_front, hh, dy, hi, col, lw_faint):
        """n coils of the same width hh and length L, each st further along the track and dy further back (up the page); coil hi is drawn in color."""
        for k in range(n - 1, -1, -1):
            x0 = x_first + k * st; yc = y_front + k * dy
            xs, ys = f6.loop_path(x0, x0 + L, yc - hh / 2, yc + hh / 2, 0.9)
            if k == hi:
                continue
            ax.plot(xs, ys, color=C["faint"], lw=lw_faint, zorder=2, solid_capstyle="round")
        x0 = x_first + hi * st; yc = y_front + hi * dy
        xs, ys = f6.loop_path(x0, x0 + L, yc - hh / 2, yc + hh / 2, 0.9)
        ax.plot(xs, ys, color=col, lw=LW2, zorder=4)
        return x0, yc
    # the support track: two separate loops side by side, each of ten coils of one width, one behind the other and offset along the track
    L, st, n, hh, dy, hi = 55.0, 3.6, 10, 6.2, 0.62, 4
    xA, yA = coils(n, L, st, 3.0, 57.5, hh, dy, hi, C["support"], LW)
    xB, yB = coils(n, L, st, 3.0, 43.5, hh, dy, hi, C["support"], LW)
    xm = xA + L / 2 + 6
    arrow(xm, yA + hh / 2, +1, C["support"]); arrow(xm, yA - hh / 2, -1, C["support"])
    arrow(xm, yB + hh / 2, -1, C["support"]); arrow(xm, yB - hh / 2, +1, C["support"])
    f6.label(ax, (50, 73.2), "support track: two separate loops,\nten coils in each, one behind the other", color=C["support"])
    # the fast track: twenty coils of one width
    Lf, stf, nf, hf, dyf, hif = 27.0, 1.5, 20, 9.0, 0.36, 9
    xf = (100 - (Lf + (nf - 1) * stf)) / 2
    xF, yF = coils(nf, Lf, stf, xf, 13.5, hf, dyf, hif, C["load"], LW * 0.7)
    xm = xF + Lf / 2
    arrow(xm, yF + hf / 2, +1, C["load"]); arrow(xm, yF - hf / 2, -1, C["load"])
    f6.label(ax, (50, 31.5), "fast track: one loop, twenty coils", color=C["load"])
    f6.label(ax, (50, 3.2), "The lengths and the offsets are not to scale.", color=C["muted"])
    return save(fig, folder, "fig07_candidate_2_two_coils")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    folder = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(here, "..", "..", "..", "Vol_III", "figure_work_2026-11-16")
    os.makedirs(folder, exist_ok=True)
    which = [a[2:] for a in sys.argv[1:] if a.startswith("--")] or ["coils", "field"]
    if "coils" in which:
        fig_two_coils(folder)
    if "field" in which:
        fig_field(folder)


if __name__ == "__main__":
    main()
