#!/usr/bin/env python3
"""
Candidate 1 for Chapter 7 redrawn (the bearing, round 7, Section 4 of Vol_III/Claude_Code_instructions_bearing_round7_2026-11-15.md): the bearing in
cross-section with the plates as ruled (fast 1.0 m x 200 mm, support 1.2 m x 100 mm), the 204 x 50 mm fast bundles, each support track as two
separate loops (round 7, Section 1; provisional), and an enlarged inset of one lower corner.
Book: Orbital Ring Engineering, Volume III: Chapter 7.  A candidate, not a working drawing; the rules are round 6b's (tools/round6b_figures.py).
Run: python tools/round7_figures.py [output folder]     (default ../../Vol_III/figure_work_2026-11-15)
The proportions are true in the main view and in the inset.  The figure is made 4.5 inches wide, so a font size is the size in print.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import round6b_figures as f6

C = f6.C; FS = f6.FS; LW, LW2 = f6.LW, f6.LW2; W_IN = f6.W_IN
SIDE = 7.39; GAP = 0.10; WALL = 0.01
GC = 0.05                                                    # between the two loops of a support track at the center (each loop its own cryostat)
SUP = dict(Wp=1.2, T=0.10, legs=(-0.475 - GC / 2, -0.075 - GC / 2, 0.075 + GC / 2, 0.475 + GC / 2), bw=0.15, bh=0.067)
FAST = dict(Wp=1.0, T=0.20, legs=(-0.30, 0.30), bw=0.204, bh=0.05)
X_FAST = SIDE / 2 - FAST["Wp"] / 2; X_SUP = SIDE / 2 - FAST["Wp"] - 0.2 - SUP["Wp"] / 2


def draw(ax, lw_scale=1.0, half=5.0):
    s = SIDE / 2
    y_in = s + FAST["T"] + GAP; x_in = s + 0.2 + 0.2
    ax.add_patch(Rectangle((-half, -half), 2 * half, 2 * half, fc=C["casing"], ec=C["casing_edge"], lw=LW2 * lw_scale, zorder=1))
    ax.add_patch(Rectangle((-x_in, -y_in), 2 * x_in, 2 * y_in, fc="white", ec=C["casing_edge"], lw=LW * lw_scale, zorder=2))
    for sx in (-1, 1):                                       # under the thinner support plates the casing's face stands higher
        ax.add_patch(Rectangle((sx * X_SUP - SUP["Wp"] / 2 - 0.05, -y_in), SUP["Wp"] + 0.1, FAST["T"] - SUP["T"], fc=C["casing"], ec="none", zorder=2.5))
        ax.plot([sx * X_SUP - SUP["Wp"] / 2 - 0.05, sx * X_SUP - SUP["Wp"] / 2 - 0.05, sx * X_SUP + SUP["Wp"] / 2 + 0.05, sx * X_SUP + SUP["Wp"] / 2 + 0.05],
                [-y_in, -y_in + FAST["T"] - SUP["T"], -y_in + FAST["T"] - SUP["T"], -y_in], color=C["casing_edge"], lw=LW * lw_scale, zorder=2.6)
    ax.add_patch(Rectangle((-s, -s), SIDE, SIDE, fc=C["cable"], ec=C["cable_edge"], lw=LW2 * lw_scale, zorder=3))
    def plate(xc, y_face, T, Wp, up):
        ax.add_patch(Rectangle((xc - Wp / 2, y_face if up else y_face - T), Wp, T, fc=C["plate"], ec="none", zorder=4))
    def track(xc, y_plate_face, d, col, up):
        for x in d["legs"]:
            y0 = y_plate_face + GAP + WALL if up else y_plate_face - GAP - WALL - d["bh"]
            ax.add_patch(Rectangle((xc + x - d["bw"] / 2, y0), d["bw"], d["bh"], fc=col, ec="none", zorder=4))
    for sx in (-1, 1):
        plate(sx * X_FAST, s, FAST["T"], FAST["Wp"], True); track(sx * X_FAST, s + FAST["T"], FAST, C["lift"], True)
        plate(sx * X_FAST, -s, FAST["T"], FAST["Wp"], False); track(sx * X_FAST, -s - FAST["T"], FAST, C["load"], False)
        plate(sx * X_SUP, -s, SUP["T"], SUP["Wp"], False); track(sx * X_SUP, -s - SUP["T"], SUP, C["support"], False)
        ax.add_patch(Rectangle((sx * s if sx > 0 else -s - 0.2, 0.3), 0.2, 1.2, fc=C["lim_plate"], ec="none", zorder=4))
        ax.add_patch(Rectangle((sx * x_in if sx > 0 else -x_in - 0.3, 0.4), 0.3, 1.0, fc=C["lim_stator"], ec="none", zorder=4))


def fig_cross_section(folder):
    s = SIDE / 2
    main_h = W_IN; inset_h = 1.75; H = main_h + inset_h
    fig = plt.figure(figsize=(W_IN, H))
    ax = fig.add_axes([0, inset_h / H, 1, main_h / H]); ax.set_xlim(-5.15, 5.15); ax.set_ylim(-5.15, 5.15); ax.set_aspect("equal"); ax.axis("off")
    draw(ax)
    x_in = s + 0.4
    f6.label(ax, (0, 0.0), "cable", color=C["cable_edge"], weight="bold")
    f6.label(ax, (0, 4.52), "casing")
    f6.label(ax, (-1.5, 2.55), "lift track", to=(-X_FAST + 0.3, s + FAST["T"] + GAP + WALL + 0.01), color=C["lift"], gap=11)
    f6.label(ax, (1.5, 2.55), "plate", to=(X_FAST - 0.3, s + FAST["T"] / 2), color="#6b5f2a", gap=9)
    f6.label(ax, (-1.75, 0.9), "LIM plate", to=(-s - 0.1, 0.9), gap=24)
    f6.label(ax, (1.7, 0.9), "LIM stator", to=(x_in + 0.15, 0.9), gap=26)
    f6.label(ax, (-0.9, -2.45), "support track", to=(-X_SUP + 0.1, -s - SUP["T"] - GAP - WALL - 0.01), color=C["support"], gap=10)
    f6.label(ax, (-2.6, -1.45), "load track", to=(-X_FAST + 0.3, -s - FAST["T"] - GAP - WALL - 0.01), color=C["load"], gap=10)
    # the corner that the inset enlarges
    bx0, bx1, by0, by1 = 1.05, 3.95, -4.22, -3.42
    ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fc="none", ec=C["ink"], lw=LW, ls=(0, (3, 2)), zorder=7))
    # the inset: the same drawing, the same proportions, enlarged
    m = 0.045
    w_ax = 1 - 2 * m; h_ax = w_ax * W_IN * (by1 - by0) / (bx1 - bx0) / H
    y_ax = 0.36 / H
    axi = fig.add_axes([m, y_ax, w_ax, h_ax]); axi.set_xlim(bx0, bx1); axi.set_ylim(by0, by1); axi.set_aspect("equal")
    draw(axi, lw_scale=1.0)
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_linewidth(LW); sp.set_linestyle((0, (3, 2))); sp.set_color(C["ink"])
    # the two lines that tie the corner to the inset
    for xf in (bx0, bx1):
        p_main = fig.transFigure.inverted().transform(ax.transData.transform((xf, by0)))
        p_in = (m if xf == bx0 else m + w_ax, y_ax + h_ax)
        fig.add_artist(plt.Line2D([p_main[0], p_in[0]], [p_main[1], p_in[1]], color=C["muted"], lw=LW * 0.8, ls=(0, (3, 2))))
    # the inset's labels, under it, each with a leader up to its part
    yl = by0 - 0.16
    def ilabel(x, text, to, col=None):
        axi.annotate(text, xy=to, xytext=(x, yl), fontsize=FS, color=col or C["ink"], ha="center", va="top", annotation_clip=False,
                     arrowprops=dict(arrowstyle="-", lw=LW * 0.8, color=C["muted"], shrinkA=1.5, shrinkB=0), zorder=8)
    y_sup_b = -s - SUP["T"] - GAP - WALL - SUP["bh"]; y_load_b = -s - FAST["T"] - GAP - WALL - FAST["bh"]
    ilabel(1.50, "support track", (X_SUP - 0.5, y_sup_b + 0.005), C["support"])
    ilabel(2.25, "plate", (X_SUP + 0.45, -s - SUP["T"] / 2), "#6b5f2a")
    ilabel(2.85, "casing", (2.60, by0 + 0.05))
    ilabel(3.50, "load track", (X_FAST + 0.3, y_load_b + 0.01), C["load"])
    return f6.save(fig, folder, "fig07_candidate_1_cross_section")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "Vol_III", "figure_work_2026-11-15")
    os.makedirs(folder, exist_ok=True)
    fig_cross_section(folder)
    print("  saved to", os.path.abspath(folder))


if __name__ == "__main__":
    main()
