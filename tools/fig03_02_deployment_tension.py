#!/usr/bin/env python3
"""
Draws Figure 3.2: the cable tension and stress, the casing's net weight and the cable's excess centrifugal force through
the deployment, from the time series lim_simulation.py writes (lim_deployment_tension.csv in the output folder).
Book: Orbital Ring Engineering, Volume I: Chapter 3, Figure 3.2 (and the deployment tension table of Section 3.2).
Run: python tools/fig03_02_deployment_tension.py [retrograde.csv] [prograde.csv] [--out=FOLDER]
     (defaults: <output folder>/lim_deployment_tension.csv, no prograde curve, figures to <output folder>/figures_ch03)
Environment: ORBITAL_RING_OUT (output folder, default runs)

Top panel: tension against the casing's delta-v from orbital speed (0 to 8,020 m/s), the same curve read as stress on the
right axis, the casing speed on an upper axis, and the three marks: prograde deployment complete (6,988 m/s), casing at
rest (7,504 m/s), retrograde deployment complete (8,020 m/s).  Middle panel: what the cable's excess force is spent on:
the casing's net weight (the bearing load) and the cable's excess force against the casing speed, their difference
(the tension over the orbital radius) shaded and drawn.  Bottom panel: tension and net weight against time in days.
Style as fig13_figures.py (DejaVu Sans, 11 pt, the house blue and red).  Revised 2026-10-23 after Paul's review.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ring_altitude as ring

BLUE, RED, GREEN, GRAY = "#2b6ca3", "#b5451b", "#2e8b57", "gray"
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})


def read_series(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = next((a[6:] for a in sys.argv[1:] if a.startswith("--out=")), ring.out_dir("figures_ch03"))
    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__); return
    retro_path = args[0] if args else os.path.join(ring.out_dir(), "lim_deployment_tension.csv")
    pro_path = args[1] if len(args) > 1 else None
    os.makedirs(out, exist_ok=True)
    R = read_series(retro_path)
    P = read_series(pro_path) if pro_path else None
    v_o = ring.V_ORBIT; v_g = ring.V_GROUND_SYNC; r_orb = ring.R_ORBIT
    vL = R["v_casing_m_s"]; dv = v_o - vL                     # the casing's delta-v from orbital speed
    T = R["tension_GN"]; S = R["stress_GPa"]; W = R["casing_net_weight_kN_m"]; F = R["cable_excess_kN_m"]
    T_over_r = T * 1e9 / r_orb / 1e3                           # kN/m, the difference between the excess force and the weight
    k_zero = int(np.argmin(np.abs(vL)))
    dv_end = dv[-1]; dv_rest = v_o
    dv_pro_end = (v_o - P["v_casing_m_s"][-1]) if P is not None else None
    T_pro_end = P["tension_GN"][-1] if P is not None else None

    fig = plt.figure(figsize=(10, 12.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.05, 0.85], hspace=0.62)
    ax1 = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1]); ax3 = fig.add_subplot(gs[2])

    # ---- top: tension against the casing's delta-v; the same curve read as stress on the right axis
    ax1.plot(dv, T, color=BLUE, lw=2.2, label="Cable tension, retrograde cable")
    if P is not None:
        ax1.plot(v_o - P["v_casing_m_s"], P["tension_GN"], color=BLUE, lw=1.4, ls="--", label="Cable tension, prograde cable")
    ax1.set_xlim(-150, dv_end + 150)
    ax1.set_ylim(0, T[-1] * 1.12)
    ax1.set_xlabel("Casing delta-v from orbital speed (m/s)")
    ax1.set_ylabel("Cable tension (GN)", color=BLUE)
    ax1b = ax1.twinx()
    ax1b.set_ylim(0, T[-1] * 1.12 * (S[-1] / T[-1]))
    ax1b.set_ylabel("Cable stress (GPa),\nthe same curve read on the right axis", color=RED)
    marks = [(dv_rest, f"casing at rest,\n{dv_rest:,.0f} m/s, {T[k_zero]:.0f} GN"), (dv_end, f"retrograde deployment\ncomplete, {dv_end:,.0f} m/s, {T[-1]:.0f} GN")]
    if dv_pro_end is not None:
        marks.insert(0, (dv_pro_end, f"prograde deployment\ncomplete, {dv_pro_end:,.0f} m/s, {T_pro_end:.0f} GN"))
    y_text = [0.10, 0.27, 0.44] if len(marks) == 3 else [0.15, 0.40]     # below the curve, each label left of its line
    for (x, lab), yt in zip(marks, y_text):
        ax1.axvline(x, color=GRAY, lw=0.9, ls="--")
        ax1.text(x - 90, T[-1] * 1.12 * yt, lab, color=GRAY, fontsize=8.5, ha="right", va="center",
                 bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    ticks = list(range(0, 7000, 1000)) + [m[0] for m in marks]
    ax1.set_xticks(ticks)
    ax1.set_xticklabels([f"{t:,.0f}" for t in ticks], rotation=35, ha="right", fontsize=9)
    # the casing speed on an upper axis: westward on the left, through rest, eastward on the right
    ax1t = ax1.secondary_xaxis("top", functions=(lambda d: v_o - d, lambda v: v_o - v))
    v_ticks = [7000, 6000, 5000, 4000, 3000, 2000, 1000, 0, -500]
    ax1t.set_xticks(v_ticks)
    ax1t.set_xticklabels([f"{v:,} W" if v > 0 else ("0" if v == 0 else f"{-v:,} E") for v in v_ticks], fontsize=9)
    ax1t.set_xlabel("Casing speed (m/s), westward (W) through rest to eastward (E)", fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper left", fontsize=9)
    ax1.set_title("Cable tension through the deployment of the 700 km ring", fontsize=12, pad=38)

    # ---- middle: what the cable's excess force is spent on
    x2 = vL
    ax2.fill_between(x2, W, F, color=BLUE, alpha=0.13, lw=0)
    ax2.plot(x2, F, color=RED, lw=2.2, label="Cable excess centrifugal force")
    ax2.plot(x2, W, color=GREEN, lw=2.2, label="Casing net weight (the bearing load)")
    ax2.plot(x2, T_over_r, color=BLUE, lw=1.1, label="Their difference, tension / r")
    if P is not None:
        xp = P["v_casing_m_s"]
        ax2.plot(xp, P["cable_excess_kN_m"], color=RED, lw=1.2, ls="--")
        ax2.plot(xp, P["casing_net_weight_kN_m"], color=GREEN, lw=1.2, ls="--")
    ax2.set_xlim(v_o + 150, -(v_g + 150))
    ax2.set_ylim(0, F[-1] * 1.18)
    ax2.set_xticks(v_ticks)
    ax2.set_xticklabels([f"{v:,} W" if v > 0 else ("0" if v == 0 else f"{-v:,} E") for v in v_ticks], fontsize=9)
    ax2.set_xlabel("Casing speed (m/s), westward through rest to eastward")
    ax2.set_ylabel("Force per meter of ring (kN/m)")
    ax2.axvline(0, color=GRAY, lw=0.9, ls="--")
    ax2.axhline(0, color=GRAY, lw=0.8)
    k13 = int(np.argmin(np.abs(vL - 1300)))
    ax2.annotate(f"difference = tension / r\n(T / r = {T_over_r[-1]:.1f} kN/m at the end)", xy=(1300, (W[k13] + F[k13]) / 2),
                 xytext=(3300, F[-1] * 0.96), color=BLUE, fontsize=9.5, ha="center", arrowprops=dict(arrowstyle="->", color=BLUE))
    ax2.text(-300, F[-1] + 4, f"{F[-1]:.1f} kN/m", color=RED, fontsize=9, ha="right")
    ax2.text(-300, W[-1] - 17, f"{W[-1]:.1f} kN/m", color=GREEN, fontsize=9, ha="right")
    ax2.text(-300, T_over_r[-1] + 4, f"{T_over_r[-1]:.1f} kN/m", color=BLUE, fontsize=9, ha="right")
    ax2.text(150, F[-1] * 1.08, f"casing at rest, day {R['time_days'][k_zero]:.0f}", color=GRAY, fontsize=8.5, ha="right")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=9)
    ax2.set_title("What the cable's excess force is spent on:\ncarrying the casing, and tensioning the cable", fontsize=12)

    # ---- bottom: against time
    t = R["time_days"]
    ax3.plot(t, T, color=BLUE, lw=2.2, label="Cable tension (GN)")
    ax3.set_ylabel("Cable tension (GN)", color=BLUE)
    ax3b = ax3.twinx()
    ax3b.plot(t, W, color=GREEN, lw=2.2, label="Casing net weight (kN/m)")
    ax3b.set_ylabel("Casing net weight\n(kN/m)", color=GREEN)
    ax3.set_xlabel("Time from the start of the deployment (days)")
    ax3.set_xlim(0, t[-1] * 1.02); ax3.set_ylim(0, T[-1] * 1.08); ax3b.set_ylim(0, 120)
    ax3.axvline(t[k_zero], color=GRAY, lw=0.9, ls="--")
    ax3.text(t[k_zero] - 8, T[-1] * 0.08, "casing at rest", color=GRAY, fontsize=9, ha="right")
    ax3.grid(True, alpha=0.3)
    l1 = ax3.get_legend_handles_labels(); l2 = ax3b.get_legend_handles_labels()
    ax3.legend(l1[0] + l2[0], l1[1] + l2[1], loc="upper left", fontsize=9)
    ax3.set_title("The same two quantities against time", fontsize=12)

    fig.subplots_adjust(left=0.10, right=0.86, top=0.90, bottom=0.06)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig03_02_deployment_tension_700.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig03_02_deployment_tension_700.svg and .png in {out}")


if __name__ == "__main__":
    main()
