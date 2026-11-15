#!/usr/bin/env python3
"""
Draws Figure 9.1: three pulsating phase fields add up to one traveling wave, with a strip of the three phase currents over one cycle above the graphs.
Book: Orbital Ring Engineering, Volume I: Chapter 9, Figure 9.1.
Run: python tools/fig09_01_traveling_wave.py [output name without extension] [--text=print8|present]
       (default <output folder>/figures_ch09/fig09_01_traveling_wave, as .png and .svg)
       --text=print8   (default) every text at least 8 points when the figure prints 4.5 inches wide
       --text=present  the text sizes of the figure before the strip was added (10 to 13 points on the 8 inch figure, 5.6 points and up in print)
Environment: ORBITAL_RING_OUT (output folder, default runs)

Phase k: current i_k(t) = cos(wt - 2 pi k/3); field b_k(x,t) = i_k(t) * cos(pi x/tau - 2 pi k/3).  Sum = 1.5 cos(wt - pi x/tau).
The strip (2026-11-15) shows i_k over one cycle with marks at the three moments of the graphs: the height of a phase's curve in a graph is
the current in that phase at that moment.  The script checks the nine currents against the curves it draws and prints them, checks that no
two texts overlap, and prints the smallest text size at 4.5 inches."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith("--")]
TEXT = "print8"
for a in sys.argv[1:]:
    if a.startswith("--text="): TEXT = a[7:]
PRINT_W = 4.5                                   # inches, the width the figure prints at on the 6 by 9 inch page
FIG_W, FIG_H = 8.0, 10.9
FS = dict(tick=14.5, label=14.5, axis=14.5, title=15.5, legend=14.5) if TEXT == "print8" else dict(tick=10, label=11, axis=12, title=13, legend=10)

COL = ["#0072B2", "#D55E00", "#009E73"]
INK, MUTED = "#222222", "#666666"
xs = [i * 0.01 for i in range(401)]            # position in pole pitches
def current(k, wt): return math.cos(wt - 2 * math.pi * k / 3)
def phase(k, x, wt): return current(k, wt) * math.cos(math.pi * x - 2 * math.pi * k / 3)
steps = [0, 60, 120]

# the nine currents, checked against the curves of the graphs: the largest height of phase k's curve at a moment is |i_k|, and its sign is the curve's sign at its own axis (x = 2k/3)
print("phase currents (peak = 1):")
for deg in steps:
    wt = math.radians(deg); row = []
    for k in range(3):
        i = current(k, wt)
        amp = max(abs(phase(k, x, wt)) for x in xs); at_axis = phase(k, 2 * k / 3, wt)
        assert abs(amp - abs(i)) < 1e-3 and abs(at_axis - i) < 1e-9, (deg, k, i, amp, at_axis)
        row.append(f"{'ABC'[k]} {i:+.1f}")
    print(f"  {deg:3d} deg: " + ", ".join(row))

fig = plt.figure(figsize=(FIG_W, FIG_H))
def box(left, bottom, width, height): return fig.add_axes([left / FIG_W, bottom / FIG_H, width / FIG_W, height / FIG_H])
L, W = 1.20, 6.45                               # inches: the left edge and the width of every plot
GH, GAP, B0 = 2.10, 0.25, 0.85                  # a graph's height, the gap between graphs, the bottom graph's lower edge
axes = [box(L, B0 + (2 - j) * (GH + GAP), W, GH) for j in range(3)]
SH = 0.74                                       # the strip: about a third of a graph's height
strip = box(L, B0 + 3 * GH + 2 * GAP + 0.90, W, SH)

# the strip: the three phase currents over one cycle
degs = [d * 0.5 for d in range(721)]
for k in range(3):
    strip.plot(degs, [current(k, math.radians(d)) for d in degs], color=COL[k], lw=1.3)
strip.axhline(0, color=MUTED, lw=0.6)
for deg in steps:
    strip.plot([deg, deg], [-1.25, 1.92], color=INK, lw=1.0, ls=(0, (4, 2)), zorder=4, clip_on=False)
strip.set_xlim(0, 360); strip.set_ylim(-1.25, 1.25)
strip.set_xticks(range(0, 361, 60)); strip.set_xticklabels([f"{d}°" for d in range(0, 361, 60)], fontsize=FS["tick"])
strip.set_yticks([-1, 0, 1]); strip.set_yticklabels(["−1", "0", "1"], fontsize=FS["tick"])
strip.grid(True, axis="x", color="#DDDDDD", lw=0.6); strip.set_axisbelow(True)
for s in ("top", "right"): strip.spines[s].set_visible(False)
strip.set_xlabel("Time, in degrees of the cycle", fontsize=FS["axis"])
strip.set_ylabel("Current\n(peak = 1)", fontsize=FS["label"])
mark_y = 2.00                                   # the labels sit clear of the y label; the marks run up to them
strip.text(0, mark_y, "Start of\nthe cycle", fontsize=FS["label"], color=INK, ha="center", va="bottom", linespacing=1.05)
strip.text(60, mark_y, "60°", fontsize=FS["label"], color=INK, ha="center", va="bottom")
strip.text(120, mark_y, "120°", fontsize=FS["label"], color=INK, ha="center", va="bottom")
strip.text(360, mark_y, "The current in each phase\nover one cycle", fontsize=FS["label"], color=INK, ha="right", va="bottom", linespacing=1.05)

# the three graphs (unchanged: the curves, the colors, the labels and the dot)
for ax, deg in zip(axes, steps):
    wt = math.radians(deg)
    for k in range(3):
        ax.plot(xs, [phase(k, x, wt) for x in xs], color=COL[k], lw=1.3, label=f"phase {'ABC'[k]}" if deg == 0 else None)
    tot = [sum(phase(k, x, wt) for k in range(3)) for x in xs]
    ax.plot(xs, tot, color=INK, lw=2.6, label="sum of the three phases" if deg == 0 else None)
    xp = deg / 180.0 + 2.0                      # one crest of the sum, followed from panel to panel
    ax.plot([xp], [1.5], "o", ms=10, color="#CC3311", markeredgecolor="white", markeredgewidth=2, zorder=5)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_ylim(-1.9, 2.5); ax.set_xlim(0, 4)
    ax.set_yticks([-1.5, 0, 1.5]); ax.set_yticklabels(["−1.5", "0", "1.5"], fontsize=FS["tick"])
    ax.set_xticks([0, 1, 2, 3, 4]); ax.tick_params(labelbottom=ax is axes[2], labelsize=FS["tick"])
    ax.grid(True, axis="x", color="#DDDDDD", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    lab = {0: "Start of the cycle", 60: "60° into the cycle: the crest has moved 1/3 of a pole pitch",
           120: "120° into the cycle: the crest has moved 2/3 of a pole pitch"}[deg]
    ax.text(0.04, 2.2, lab, fontsize=FS["label"], color=INK, va="center")
axes[1].set_ylabel("Field strength (one phase at its peak = 1)", fontsize=FS["label"])
axes[2].set_xlabel("Position along the stator, in pole pitches", fontsize=FS["axis"])
fig.suptitle("Three pulsating fields add up to one wave that travels", fontsize=FS["title"], color=INK, y=1 - 0.20 / FIG_H)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1 - 0.55 / FIG_H), ncol=4, frameon=False, fontsize=FS["legend"],
           handlelength=1.3, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)

# checks: every text inside the figure, no two texts overlapping, the smallest text at the print width
fig.canvas.draw(); rend = fig.canvas.get_renderer()
texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
boxes = [(t, t.get_window_extent(rend)) for t in texts]
fb = fig.bbox
bad = [t.get_text() for t, b in boxes if b.x0 < fb.x0 or b.x1 > fb.x1 or b.y0 < fb.y0 or b.y1 > fb.y1]
over = [(a.get_text(), c.get_text()) for i, (a, ba) in enumerate(boxes) for c, bc in boxes[i + 1:]
        if ba.x0 < bc.x1 and bc.x0 < ba.x1 and ba.y0 < bc.y1 and bc.y0 < ba.y1]
small = min(t.get_fontsize() for t in texts)
print(f"texts: {len(texts)}; outside the figure: {bad or 'none'}; overlapping pairs: {over or 'none'}")
print(f"smallest text: {small:g} pt on the {FIG_W:g} inch figure = {small * PRINT_W / FIG_W:.2f} pt at {PRINT_W} inches wide ({FIG_H * PRINT_W / FIG_W:.2f} inches tall)")
print(f"strip height / graph height = {SH / GH:.2f}")

import ring_altitude as ring
out = args[0] if args else _os.path.join(ring.out_dir("figures_ch09"), "fig09_01_traveling_wave")   # under ORBITAL_RING_OUT (default runs)
fig.savefig(out + ".png", dpi=300, facecolor="white"); fig.savefig(out + ".svg", facecolor="white"); print("saved", out)
