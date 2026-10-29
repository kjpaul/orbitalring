#!/usr/bin/env python3
"""
Draws Figure 9.1: three pulsating phase fields add up to one traveling wave.
Book: Orbital Ring Engineering, Volume I: Chapter 9, Figure 9.1.
Run: python tools/fig09_01_traveling_wave.py [output name without extension]   (default <output folder>/figures_ch09/fig09_01_traveling_wave, as .png and .svg)
Environment: ORBITAL_RING_OUT (output folder, default runs)

Figure 9.1: three pulsating phase fields add up to one traveling wave.
Phase k: b_k(x,t) = cos(wt - 2 pi k/3) * cos(pi x/tau - 2 pi k/3).  Sum = 1.5 cos(wt - pi x/tau)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = ["#0072B2", "#D55E00", "#009E73"]
INK, MUTED = "#222222", "#666666"
xs = [i * 0.01 for i in range(401)]            # position in pole pitches
def phase(k, x, wt): return math.cos(wt - 2 * math.pi * k / 3) * math.cos(math.pi * x - 2 * math.pi * k / 3)
steps = [0, 60, 120]
fig, axes = plt.subplots(3, 1, figsize=(8, 7.6), sharex=True)
peaks = []
for ax, deg in zip(axes, steps):
    wt = math.radians(deg)
    for k in range(3):
        ax.plot(xs, [phase(k, x, wt) for x in xs], color=COL[k], lw=1.3, label=f"phase {'ABC'[k]}" if deg == 0 else None)
    tot = [sum(phase(k, x, wt) for k in range(3)) for x in xs]
    ax.plot(xs, tot, color=INK, lw=2.6, label="sum of the three phases" if deg == 0 else None)
    xp = deg / 180.0 + 2.0                      # one crest of the sum, followed from panel to panel
    ax.plot([xp], [1.5], "o", ms=10, color="#CC3311", markeredgecolor="white", markeredgewidth=2, zorder=5)
    peaks.append(xp)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_ylim(-1.9, 2.5); ax.set_xlim(0, 4)
    ax.set_yticks([-1.5, 0, 1.5]); ax.set_yticklabels(["−1.5", "0", "1.5"])
    ax.grid(True, axis="x", color="#DDDDDD", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    lab = {0: "Start of the cycle", 60: "60° into the cycle: the crest has moved 1/3 of a pole pitch",
           120: "120° into the cycle: the crest has moved 2/3 of a pole pitch"}[deg]
    ax.text(0.04, 2.2, lab, fontsize=11, color=INK, va="center")
axes[1].set_ylabel("Field strength (one phase at its peak = 1)", fontsize=11)
axes[2].set_xticks([0, 1, 2, 3, 4]); axes[2].set_xlabel("Position along the stator, in pole pitches", fontsize=12)
fig.suptitle("Three pulsating fields add up to one wave that travels", fontsize=13, color=INK)
fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4, frameon=False, fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95))
import ring_altitude as ring
out = sys.argv[1] if len(sys.argv) > 1 else _os.path.join(ring.out_dir("figures_ch09"), "fig09_01_traveling_wave")   # under ORBITAL_RING_OUT (default runs)
fig.savefig(out + ".png", dpi=300, facecolor="white"); fig.savefig(out + ".svg", facecolor="white"); print("saved")
