#!/usr/bin/env python3
"""
Four candidate figures for Chapter 7 (the bearing, round 6 part B, Section 3 of Vol_III/Claude_Code_instructions_bearing_round6b_2026-11-14.md).
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Candidates beside Figures 7.3 and 7.4; not working drawings.
Run: python tools/round6b_figures.py [output folder]     (default ../../Vol_III/figure_work_2026-11-14)
  fig07_candidate_1_cross_section    the bearing in cross-section, the proportions true
  fig07_candidate_2_two_coils        the two kinds of coil, seen from above (lengths not to scale)
  fig07_candidate_3_support_circuit  the support coil's circuit (the dump resistor in series with the diode, its switch across it)
  fig07_candidate_4_fast_circuit     the fast coil's circuit (the two-switch, two-diode bridge)
Rules (the brief): one idea per figure, at most ten labels, names of parts only, no numbers, no title inside; the figure is made 4.5 inches wide, the
width it prints at, so a font size here is the size in print; the smallest text is 8.5 points.  The script checks that no two texts overlap, that
none leaves the figure, and prints the smallest text; each PNG is to be opened and looked at.
The font, colors and line weights follow the book's Figures 7.1 and 7.3 (Segoe UI or Arial; blue, green and orange-red parts, a brown plate, gray ink).
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Arc, Polygon, Circle

plt.rcParams["font.family"] = ["Segoe UI", "Arial", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
W_IN = 4.5
FS = 8.5
LW, LW2 = 0.7, 1.4                                 # points in print: the book's 1.5 and 3 px on its 960 px figures, a little heavier for the small width
C = dict(ink="#222222", muted="#777777", wire="#444444", support="#cc4400", load="#2266aa", lift="#228833", plate="#9a8a4a", cable="#cfe0f0", cable_edge="#1b4a73",
         casing="#ececec", casing_edge="#444444", lim_plate="#b9b9b9", lim_stator="#666666", cold="#1b4a73", faint="#b8c4d0", strip="#cfe0f0")

# the design drawn (round 5's optimum plates, not yet accepted; the proportions, not the numbers, are what the figure shows)
SIDE = 7.39; GAP = 0.10
SUP = dict(Wp=1.2, T=0.10, legs=(-0.475, -0.075, 0.075, 0.475), bw=0.15, bh=0.067)
FAST = dict(Wp=1.0, T=0.15, legs=(-0.30, 0.30), bw=0.20, bh=0.05)
X_FAST = SIDE / 2 - FAST["Wp"] / 2; X_SUP = SIDE / 2 - FAST["Wp"] - 0.2 - SUP["Wp"] / 2
WALL = 0.01


def new_fig(h_in, xlim, ylim, equal=True):
    fig = plt.figure(figsize=(W_IN, h_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    if equal:
        ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def label(ax, xy, text, to=None, ha="center", va="center", color=None, weight="normal", gap=0.0):
    t = ax.text(xy[0], xy[1], text, fontsize=FS, ha=ha, va=va, color=color or C["ink"], fontweight=weight, linespacing=1.05, zorder=6)
    if to is not None:
        ax.annotate("", xy=to, xytext=(xy[0], xy[1]), arrowprops=dict(arrowstyle="-", lw=LW * 0.8, color=C["muted"], shrinkA=gap, shrinkB=0), zorder=5)
    return t


def check(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    bbs = [(t, t.get_window_extent(r)) for t in texts]; fb = fig.bbox
    out = [t.get_text() for t, b in bbs if b.x0 < fb.x0 or b.x1 > fb.x1 or b.y0 < fb.y0 or b.y1 > fb.y1]
    over = [(a.get_text(), c.get_text()) for i, (a, ba) in enumerate(bbs) for c, bc in bbs[i + 1:] if ba.x0 < bc.x1 and bc.x0 < ba.x1 and ba.y0 < bc.y1 and bc.y0 < ba.y1]
    small = min(t.get_fontsize() for t in texts)
    print(f"  {name}: {len(texts)} labels; smallest text {small:g} pt at {W_IN} inches wide; outside the figure: {out or 'none'}; overlapping: {over or 'none'}; {fig.get_figheight():.2f} inches tall")
    return len(texts), small, out, over


def save(fig, folder, name):
    res = check(fig, name)
    fig.savefig(os.path.join(folder, name + ".png"), dpi=300, facecolor="white")
    fig.savefig(os.path.join(folder, name + ".svg"), facecolor="white")
    plt.close(fig)
    return res


# ---------------------------------------------------------------------------------------------------------------- 1: the cross-section
def fig_cross_section(folder):
    half = 5.0
    fig, ax = new_fig(W_IN * 10.3 / 10.3, (-5.15, 5.15), (-5.15, 5.15))
    s = SIDE / 2
    y_in = s + FAST["T"] + GAP                       # the casing's inner face above and below the fast plates
    x_in = s + 0.2 + 0.2                             # the LIM plate and its gap
    # the casing: a frame around the cable
    ax.add_patch(Rectangle((-half, -half), 2 * half, 2 * half, fc=C["casing"], ec=C["casing_edge"], lw=LW2, zorder=1))
    ax.add_patch(Rectangle((-x_in, -y_in), 2 * x_in, 2 * y_in, fc="white", ec=C["casing_edge"], lw=LW, zorder=2))
    # under the support plates the casing's face stands higher, so that every track is one gap from its plate
    for sx in (-1, 1):
        ax.add_patch(Rectangle((sx * X_SUP - SUP["Wp"] / 2 - 0.05, -y_in), SUP["Wp"] + 0.1, FAST["T"] - SUP["T"], fc=C["casing"], ec="none", zorder=2.5))
    # the cable
    ax.add_patch(Rectangle((-s, -s), SIDE, SIDE, fc=C["cable"], ec=C["cable_edge"], lw=LW2, zorder=3))
    # plates and tracks
    def plate(xc, y_face, T, Wp, up):
        ax.add_patch(Rectangle((xc - Wp / 2, y_face if up else y_face - T), Wp, T, fc=C["plate"], ec="none", zorder=4))
    def track(xc, y_plate_face, d, col, up):
        for x in d["legs"]:
            y0 = y_plate_face + GAP + WALL if up else y_plate_face - GAP - WALL - d["bh"]
            ax.add_patch(Rectangle((xc + x - d["bw"] / 2, y0), d["bw"], d["bh"], fc=col, ec="none", zorder=4))
    for sx in (-1, 1):
        plate(sx * X_FAST, s, FAST["T"], FAST["Wp"], True); track(sx * X_FAST, s + FAST["T"], FAST, C["lift"], True)                # lift: on top
        plate(sx * X_FAST, -s, FAST["T"], FAST["Wp"], False); track(sx * X_FAST, -s - FAST["T"], FAST, C["load"], False)            # load: underneath, at the edges
        plate(sx * X_SUP, -s, SUP["T"], SUP["Wp"], False); track(sx * X_SUP, -s - SUP["T"], SUP, C["support"], False)               # support: underneath, inboard
        # the LIM: a reaction plate on the cable's side, a stator in the casing across the gap, above the center
        ax.add_patch(Rectangle((sx * s if sx > 0 else -s - 0.2, 0.3), 0.2, 1.2, fc=C["lim_plate"], ec="none", zorder=4))
        ax.add_patch(Rectangle((sx * x_in if sx > 0 else -x_in - 0.3, 0.4), 0.3, 1.0, fc=C["lim_stator"], ec="none", zorder=4))
    # labels, inside the cable, each with a leader to its part
    label(ax, (0, 0.0), "cable", color=C["cable_edge"], weight="bold")
    label(ax, (0, 4.52), "casing")
    label(ax, (-1.5, 2.55), "lift track", to=(-X_FAST + 0.3, s + FAST["T"] + GAP + WALL + 0.01), color=C["lift"], gap=11)
    label(ax, (1.5, 2.55), "plate", to=(X_FAST - 0.3, s + FAST["T"] / 2), color="#6b5f2a", gap=9)
    label(ax, (-1.75, 0.9), "LIM plate", to=(-s - 0.1, 0.9), gap=24)
    label(ax, (1.7, 0.9), "LIM stator", to=(x_in + 0.15, 0.9), gap=26)
    label(ax, (0, -2.45), "support track", to=(-X_SUP, -s - SUP["T"] - GAP - WALL - 0.01), color=C["support"], gap=10)
    ax.annotate("", xy=(X_SUP, -s - SUP["T"] - GAP - WALL - 0.01), xytext=(0, -2.45), arrowprops=dict(arrowstyle="-", lw=LW * 0.8, color=C["muted"], shrinkA=10, shrinkB=0), zorder=5)
    label(ax, (2.35, -1.45), "load track", to=(X_FAST - 0.3, -s - FAST["T"] - GAP - WALL - 0.01), color=C["load"], gap=10)
    return save(fig, folder, "fig07_candidate_1_cross_section")


# ---------------------------------------------------------------------------------------------------------------- 2: the two kinds of coil
def loop_path(x0, x1, y0, y1, r):
    """A rounded rectangle as points."""
    pts = []
    for cx, cy, a0 in ((x1 - r, y1 - r, 0), (x0 + r, y1 - r, 90), (x0 + r, y0 + r, 180), (x1 - r, y0 + r, 270)):
        for k in range(13):
            a = math.radians(a0 + 90 * k / 12)
            pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    pts.append(pts[0])
    return [p[0] for p in pts], [p[1] for p in pts]


def fig_two_coils(folder):
    fig, ax = new_fig(3.1, (0, 100), (0, 100 * 3.1 / W_IN))
    # the support track: ten coils, each a figure-eight of two loops, staggered along the track
    L, st, n = 58.0, 4.0, 10; x_first = 3.0
    yA, yB, hh = 55.5, 45.5, 8.6                                        # the two loops' center lines and height
    hi = 4
    for k in list(range(n)):
        if k == hi:
            continue
        ins = 0.37 * k
        for yc in (yA, yB):
            xs, ys = loop_path(x_first + k * st, x_first + k * st + L, yc - hh / 2 + ins, yc + hh / 2 - ins, 0.9)
            ax.plot(xs, ys, color=C["faint"], lw=LW, zorder=2, solid_capstyle="round")
    ins = 0.37 * hi; xa, xb = x_first + hi * st, x_first + hi * st + L
    for yc in (yA, yB):
        xs, ys = loop_path(xa, xb, yc - hh / 2 + ins, yc + hh / 2 - ins, 0.9)
        ax.plot(xs, ys, color=C["support"], lw=LW2, zorder=4)
    # the sense of the current: opposite in the two loops, the same in the two legs at the center
    def arrow(x, y, d, col):
        ax.annotate("", xy=(x + 2.4 * d, y), xytext=(x - 2.4 * d, y), arrowprops=dict(arrowstyle="-|>", lw=0, color=col, mutation_scale=9), zorder=5)
    xm = (xa + xb) / 2 + 8
    arrow(xm, yA + hh / 2 - ins, +1, C["support"]); arrow(xm, yA - hh / 2 + ins, -1, C["support"])
    arrow(xm, yB + hh / 2 - ins, -1, C["support"]); arrow(xm, yB - hh / 2 + ins, +1, C["support"])
    # the crossover between the loops and the two ends, at one end of the coil
    ax.plot([xa, xa - 1.6, xa - 1.6, xa], [yA - 1.2, yA - 1.2, yB + 1.2, yB + 1.2], color=C["support"], lw=LW2, zorder=4)
    label(ax, (50, 66.0), "support coil, one of ten that share the track", color=C["support"])
    # the fast track: twenty coils, each a single loop
    Lf, stf, nf = 29.0, 1.45, 20; xf = (100 - (Lf + (nf - 1) * stf)) / 2
    yF, hf = 18.5, 14.0; hif = 9
    for k in range(nf):
        if k == hif:
            continue
        ins = 0.3 * k
        xs, ys = loop_path(xf + k * stf, xf + k * stf + Lf, yF - hf / 2 + ins, yF + hf / 2 - ins, 0.9)
        ax.plot(xs, ys, color=C["faint"], lw=LW * 0.7, zorder=2)
    ins = 0.3 * hif
    xs, ys = loop_path(xf + hif * stf, xf + hif * stf + Lf, yF - hf / 2 + ins, yF + hf / 2 - ins, 0.9)
    ax.plot(xs, ys, color=C["load"], lw=LW2, zorder=4)
    xm = xf + hif * stf + Lf / 2
    arrow(xm, yF + hf / 2 - ins, +1, C["load"]); arrow(xm, yF - hf / 2 + ins, -1, C["load"])
    label(ax, (50, 30.5), "fast coil, one of twenty that share the track", color=C["load"])
    label(ax, (50, 5.2), "The lengths are not to scale.", color=C["muted"])
    return save(fig, folder, "fig07_candidate_2_two_coils")


# ---------------------------------------------------------------------------------------------------------------- circuit symbols
def wire(ax, pts, col=None, lw=LW2 * 0.8):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col or C["wire"], lw=lw, solid_joinstyle="miter", solid_capstyle="butt", zorder=3)


def dot(ax, x, y, r=0.9):
    ax.add_patch(Circle((x, y), r, fc=C["wire"], ec="none", zorder=4))


def diode_v(ax, x, y0, y1, col=None, s=2.6):
    """A diode in a vertical wire from y0 up to y1, conducting upward."""
    col = col or C["wire"]; ym = (y0 + y1) / 2
    wire(ax, [(x, y0), (x, ym - s)]); wire(ax, [(x, ym + s), (x, y1)])
    ax.add_patch(Polygon([(x - s, ym - s), (x + s, ym - s), (x, ym + s)], closed=True, fc=col, ec=col, lw=LW, zorder=4))
    ax.plot([x - s, x + s], [ym + s, ym + s], color=col, lw=LW2, zorder=4)


def resistor_v(ax, x, y0, y1, col=None, w=1.8, n=6):
    col = col or C["wire"]; pad = (y1 - y0) * 0.14
    ya, yb = y0 + pad, y1 - pad
    xs = [x, x]; ys = [y0, ya]
    for k in range(n):
        xs.append(x + (w if k % 2 == 0 else -w)); ys.append(ya + (yb - ya) * (k + 0.5) / n)
    xs += [x, x]; ys += [yb, y1]
    ax.plot(xs, ys, color=col, lw=LW2 * 0.8, solid_joinstyle="miter", zorder=3)


def switch_v(ax, x, y0, y1, col=None, closed=False):
    col = col or C["wire"]; ym = (y0 + y1) / 2; g = 3.2
    wire(ax, [(x, y0), (x, ym - g)]); wire(ax, [(x, ym + g), (x, y1)])
    ax.plot([x, x + (0 if closed else 2.6)], [ym - g, ym + g * 0.75], color=col, lw=LW2, zorder=4)
    for yy in (ym - g, ym + g):
        ax.add_patch(Circle((x, yy), 0.75, fc="white", ec=col, lw=LW, zorder=5))


def switch_h(ax, x0, x1, y, col=None):
    col = col or C["wire"]; xm = (x0 + x1) / 2; g = 3.2
    wire(ax, [(x0, y), (xm - g, y)]); wire(ax, [(xm + g, y), (x1, y)])
    ax.plot([xm - g, xm + g * 0.75], [y, y + 2.6], color=col, lw=LW2, zorder=4)
    for xx in (xm - g, xm + g):
        ax.add_patch(Circle((xx, y), 0.75, fc="white", ec=col, lw=LW, zorder=5))


def inductor_v(ax, x, y0, y1, col, n=4):
    pad = (y1 - y0) * 0.1; ya, yb = y0 + pad, y1 - pad; r = (yb - ya) / (2 * n)
    wire(ax, [(x, y0), (x, ya)]); wire(ax, [(x, yb), (x, y1)])
    for k in range(n):
        ax.add_patch(Arc((x, ya + r * (2 * k + 1)), 2 * r, 2 * r, theta1=-90, theta2=90, color=col, lw=LW2 * 1.2, zorder=4))


def inductor_h(ax, x0, x1, y, col, n=4):
    pad = (x1 - x0) * 0.1; xa, xb = x0 + pad, x1 - pad; r = (xb - xa) / (2 * n)
    wire(ax, [(x0, y), (xa, y)]); wire(ax, [(xb, y), (x1, y)])
    for k in range(n):
        ax.add_patch(Arc((xa + r * (2 * k + 1), y), 2 * r, 2 * r, theta1=0, theta2=180, color=col, lw=LW2 * 1.2, zorder=4))


def capacitor_v(ax, x, y0, y1, col=None, w=3.2, g=1.1):
    col = col or C["wire"]; ym = (y0 + y1) / 2
    wire(ax, [(x, y0), (x, ym - g)]); wire(ax, [(x, ym + g), (x, y1)])
    for yy in (ym - g, ym + g):
        ax.plot([x - w, x + w], [yy, yy], color=col, lw=LW2 * 1.3, zorder=4)


# ---------------------------------------------------------------------------------------------------------------- 3: the support coil's circuit
def fig_support_circuit(folder):
    H = 3.0
    fig, ax = new_fig(H, (0, 100), (0, 100 * H / W_IN))
    yt, yb = 50.0, 14.0
    xw = 24.0                                                            # the line between warm and cold
    ax.plot([xw, xw], [5, 63], color=C["cold"], lw=LW, ls=(0, (4, 3)), zorder=2)
    label(ax, (xw - 2.0, 61.0), "warm", ha="right", color=C["muted"]); label(ax, (xw + 2.0, 61.0), "cold", ha="left", color=C["cold"])
    # the supply
    ax.add_patch(FancyBboxPatch((5.5, 24.5), 11, 15, boxstyle="round,pad=0,rounding_size=1.2", fc="#f4f4f4", ec=C["wire"], lw=LW2 * 0.8, zorder=3))
    label(ax, (11, 32), "supply")
    wire(ax, [(11, 39.5), (11, yt), (xw - 4, yt)]); wire(ax, [(11, 24.5), (11, yb), (xw - 4, yb)])
    # the leads, each with its contact, across the line
    switch_h(ax, xw - 4, xw + 8, yt); switch_h(ax, xw - 4, xw + 8, yb)
    label(ax, (xw + 3.2, 42.5), "leads", ha="left")
    x_d, x_s, x_st, x_h, x_c = 45.0, 53.0, 70.0, 75.0, 92.0
    wire(ax, [(xw + 8, yt), (x_c, yt)]); wire(ax, [(xw + 8, yb), (x_c, yb)])
    # the diode, and in series with it the dump resistor with the dump switch across it
    y_n = 28.0; y_r1 = 46.0
    diode_v(ax, x_d, yb, y_n)
    resistor_v(ax, x_d, y_n, y_r1); wire(ax, [(x_d, y_r1), (x_d, yt)])
    wire(ax, [(x_d, y_n + 1.0), (x_s, y_n + 1.0)]); wire(ax, [(x_d, y_r1 - 1.0), (x_s, y_r1 - 1.0)])
    switch_v(ax, x_s, y_n + 1.0, y_r1 - 1.0, closed=True)
    for p in ((x_d, y_n + 1.0), (x_d, y_r1 - 1.0), (x_d, yt), (x_d, yb), (x_st, yt), (x_st, yb)):
        dot(ax, *p)
    label(ax, (x_d - 3.6, 21.0), "diode", ha="right")
    label(ax, (x_d - 3.6, 37.0), "dump\nresistor", ha="right")
    label(ax, (x_s + 2.8, 37.0), "dump\nswitch", ha="left")
    # the strip with its heater
    wire(ax, [(x_st, yb), (x_st, 22)]); wire(ax, [(x_st, 42), (x_st, yt)])
    ax.add_patch(Rectangle((x_st - 1.1, 22), 2.2, 20, fc=C["strip"], ec=C["cold"], lw=LW2 * 0.8, zorder=4))
    resistor_v(ax, x_h, 25, 39, col=C["support"], w=1.3, n=8)
    label(ax, (x_st - 2.8, 25.5), "strip", ha="right")
    label(ax, (x_h + 2.6, 32), "heater", ha="left")
    # the coil
    inductor_v(ax, x_c, yb, yt, C["support"])
    label(ax, (x_c - 3.0, 7.5), "support coil", color=C["support"])
    return save(fig, folder, "fig07_candidate_3_support_circuit")


# ---------------------------------------------------------------------------------------------------------------- 4: the fast coil's circuit
def fig_fast_circuit(folder):
    H = 2.6
    fig, ax = new_fig(H, (0, 100), (4.0, 4.0 + 100 * H / W_IN))
    yt, yb, ym = 56.0, 10.0, 33.0
    x_c, x1, x2 = 9.0, 30.0, 82.0
    wire(ax, [(x_c, yt), (x2, yt)]); wire(ax, [(x_c, yb), (x2, yb)])
    capacitor_v(ax, x_c, yb, yt)
    label(ax, (x_c + 4.6, ym), "store", ha="left")
    # the left leg: a switch from the upper rail, a diode up from the lower rail; the right leg the other way round
    switch_v(ax, x1, ym, yt); diode_v(ax, x1, yb, ym)
    diode_v(ax, x2, ym, yt); switch_v(ax, x2, yb, ym)
    for p in ((x1, yt), (x1, yb), (x1, ym), (x2, ym)):
        dot(ax, *p)
    label(ax, (x1 + 4.4, (ym + yt) / 2), "switch", ha="left")
    label(ax, (x1 + 4.4, (ym + yb) / 2), "diode", ha="left")
    # the coil in the cold, between the two legs; the leads cross the line
    bx0, bx1, by0, by1 = 44.5, 69.5, 20.5, 45.5
    ax.add_patch(FancyBboxPatch((bx0, by0), bx1 - bx0, by1 - by0, boxstyle="round,pad=0,rounding_size=2.0", fc="none", ec=C["cold"], lw=LW, ls=(0, (4, 3)), zorder=2))
    wire(ax, [(x1, ym), (bx0 + 4, ym)]); wire(ax, [(bx1 - 4, ym), (x2, ym)])
    inductor_h(ax, bx0 + 4, bx1 - 4, ym, C["load"])
    label(ax, ((bx0 + bx1) / 2, 25.5), "fast coil", color=C["load"])
    label(ax, ((bx0 + bx1) / 2, 41.5), "cold", color=C["cold"])
    label(ax, ((bx0 + bx1) / 2, 50.5), "warm", color=C["muted"])
    label(ax, (bx1 + 1.6, ym - 4.4), "leads", ha="left")
    return save(fig, folder, "fig07_candidate_4_fast_circuit")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "Vol_III", "figure_work_2026-11-14")
    os.makedirs(folder, exist_ok=True)
    for f in (fig_cross_section, fig_two_coils, fig_support_circuit, fig_fast_circuit):
        f(folder)
    print("  saved to", os.path.abspath(folder))


if __name__ == "__main__":
    main()
