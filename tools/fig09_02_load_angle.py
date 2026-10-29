#!/usr/bin/env python3
"""
Draws the load angle figure of a linear synchronous motor (Figure 9.2 in the chapter; the file keeps the name it was written under): stator wave poles and sled magnets lagging by 60 degrees, and thrust against load angle.
Book: Orbital Ring Engineering, Volume I: Chapter 9, Figure 9.2.
Run: python tools/fig09_02_load_angle.py [output name without extension]   (default <output folder>/figures_ch09/fig09_02_load_angle, as .png and .svg)
Environment: ORBITAL_RING_OUT (output folder, default runs)

Figure 9.2: the load angle of a linear synchronous motor.
Top: stator wave poles and sled magnets, sled lagging by a load angle of 60 degrees.
Bottom: thrust as a fraction of the pull-out thrust, sin(delta)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

BLUE, ORANGE, INK, MUTED = "#0072B2", "#D55E00", "#222222", "#666666"
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7.2), gridspec_kw={"height_ratios": [1, 1.35]})

# --- top: poles ---
tau = 1.0                      # one pole pitch
delta = 60.0
lag = delta / 180.0 * tau
def row(ax, x0, y, label, color, n=5):
    for i in range(n):
        x = x0 + i * tau
        pole = "N" if i % 2 == 0 else "S"
        ax.add_patch(Rectangle((x + 0.04, y), tau - 0.08, 0.32, facecolor=color if pole == "N" else "white",
                               edgecolor=color, linewidth=1.6))
        ax.text(x + tau / 2, y + 0.16, pole, ha="center", va="center", fontsize=13,
                color="white" if pole == "N" else color, fontweight="bold")
    ax.text(x0 - 0.12, y + 0.16, label, ha="right", va="center", fontsize=11, color=INK)
row(ax1, 0.0, 0.95, "Stator field\n(traveling wave)", BLUE)
# sled poles are opposite to the stator pole they are attracted to, and lag behind it
for i in range(4):
    x = -lag + i * tau
    pole = "S" if i % 2 == 0 else "N"
    ax1.add_patch(Rectangle((x + 0.04, 0.28), tau - 0.08, 0.32, facecolor=ORANGE if pole == "N" else "white",
                            edgecolor=ORANGE, linewidth=1.6))
    ax1.text(x + tau / 2, 0.44, pole, ha="center", va="center", fontsize=13,
             color="white" if pole == "N" else ORANGE, fontweight="bold")
ax1.text(-lag - 0.12, 0.44, "Sled magnets", ha="right", va="center", fontsize=11, color=INK)
# pole centre markers and lag dimension
xc_st, xc_sl = tau / 2, tau / 2 - lag
ax1.plot([xc_st, xc_st], [0.60, 0.95], color=MUTED, lw=1, ls=":")
ax1.plot([xc_sl, xc_sl], [0.60, 0.80], color=MUTED, lw=1, ls=":")
ax1.annotate("", xy=(xc_st, 0.74), xytext=(xc_sl, 0.74), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
ax1.text((xc_st + xc_sl) / 2, 0.80, "lag", ha="center", va="bottom", fontsize=10, color=INK)
ax1.text(3.78, 0.44, "Load angle 60°:\nthe sled poles trail\nthe stator poles by one\nthird of a pole pitch.",
         ha="left", va="center", fontsize=10, color=INK)
ax1.add_patch(FancyArrowPatch((3.2, 1.45), (4.6, 1.45), arrowstyle="-|>", mutation_scale=16, color=BLUE, lw=1.6))
ax1.text(3.9, 1.50, "wave and sled move this way", ha="center", va="bottom", fontsize=10, color=INK)
ax1.set_xlim(-1.75, 5.3); ax1.set_ylim(0.15, 1.75); ax1.axis("off")
ax1.set_title("The sled trails the wave, and the pull between unlike poles drives it forward", fontsize=12, color=INK)

# --- bottom: thrust against load angle ---
xs = [i * 0.5 for i in range(361)]
ys = [math.sin(math.radians(x)) * 100 for x in xs]
ax2.axvspan(0, 60, color=BLUE, alpha=0.10, lw=0)
ax2.axvspan(90, 180, color=MUTED, alpha=0.10, lw=0)
ax2.plot(xs, ys, color=BLUE, lw=2)
for x, lab, dx, dy, ha in ((30, "30°: 50%", 6, -9, "left"), (60, "60°: 87%\noperating limit", -6, 3, "right"),
                           (90, "90°: 100%\npull-out limit", 0, 5, "center")):
    y = math.sin(math.radians(x)) * 100
    ax2.plot([x], [y], "o", ms=8, color=BLUE, markeredgecolor="white", markeredgewidth=2)
    ax2.text(x + dx, y + dy, lab, ha=ha, va="bottom" if dy > 0 else "top", fontsize=10, color=INK)
ax2.text(30, 8, "operating range", ha="center", fontsize=10, color=INK)
ax2.text(135, 8, "past pull-out: thrust falls\nand the sled loses step", ha="center", fontsize=10, color=INK)
ax2.set_xlim(0, 180); ax2.set_ylim(0, 125)
ax2.set_xticks(range(0, 181, 30)); ax2.set_xticklabels([f"{t}°" for t in range(0, 181, 30)])
ax2.set_yticks(range(0, 101, 25)); ax2.set_yticklabels([f"{t}%" for t in range(0, 101, 25)])
ax2.set_xlabel("Load angle δ", fontsize=12); ax2.set_ylabel("Thrust, share of pull-out thrust", fontsize=12)
ax2.grid(True, color="#DDDDDD", lw=0.6); ax2.set_axisbelow(True)
for s in ("top", "right"): ax2.spines[s].set_visible(False)
plt.tight_layout()
import sys
import ring_altitude as ring
out = sys.argv[1] if len(sys.argv) > 1 else _os.path.join(ring.out_dir("figures_ch09"), "fig09_02_load_angle")   # under ORBITAL_RING_OUT (default runs)
fig.savefig(out + ".png", dpi=300, facecolor="white"); fig.savefig(out + ".svg", facecolor="white")
print("saved", out)
