#!/usr/bin/env python3
"""
Draws Figure 3.2: the cable tension, the cable stress, the casing's net weight and the cable's excess centrifugal force
through the deployment, from the time series lim_simulation.py writes (lim_deployment_tension.csv in the output folder).
Book: The Orbital Ring, Volume IV: Chapter 3, Figure 3.2 (and the deployment tension table of Section 3.4).
Run: python tools/fig03_02_deployment_tension.py [retrograde.csv] [prograde.csv] [--out=FOLDER]
     (defaults: <output folder>/lim_deployment_tension.csv, no prograde curve, figures to <output folder>/figures_ch03)
Environment: ORBITAL_RING_OUT (output folder, default runs)

Top panel: tension (GN) and stress (GPa) against the casing's speed, from orbital speed westward through zero to the
ground-synchronous speed eastward; middle: the casing's net weight and the cable's excess force (kN/m) on the same axis;
bottom: tension and net weight against time in days.  Style as fig17_figures.py (DejaVu Sans, 11 pt, the house blue and red).
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
    d = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return d


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
    v_o = ring.V_ORBIT; v_g = ring.V_GROUND_SYNC
    vL = R["v_casing_m_s"]
    k_zero = int(np.argmin(np.abs(vL)))

    fig = plt.figure(figsize=(10, 11))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.0, 0.85], hspace=0.32)
    ax1 = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1], sharex=ax1); ax3 = fig.add_subplot(gs[2])

    # top: tension and stress against the casing speed (westward positive on the left, eastward on the right)
    ax1.plot(vL, R["tension_GN"], color=BLUE, lw=2.0, label="Cable tension, retrograde cable")
    if P is not None:
        ax1.plot(P["v_casing_m_s"], P["tension_GN"], color=BLUE, lw=1.4, ls="--", label="Cable tension, prograde cable")
    ax1.set_ylabel("Cable tension (GN)", color=BLUE)
    ax1b = ax1.twinx()
    ax1b.plot(vL, R["stress_GPa"], color=RED, lw=1.2, ls=":", label="Cable stress")
    ax1b.set_ylabel("Cable stress (GPa)", color=RED)
    ax1b.set_ylim(0, 14); ax1.set_ylim(0, 14 * (R["tension_GN"][-1] / R["stress_GPa"][-1]))
    ax1.axhline(0, color=GRAY, lw=0.8)
    ax1.annotate(f"{R['tension_GN'][-1]:.1f} GN, {R['stress_GPa'][-1]:.1f} GPa at the end", xy=(vL[-1], R["tension_GN"][-1]),
                 xytext=(3300, R["tension_GN"][-1] * 1.02), color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
    ax1.annotate(f"{R['tension_GN'][k_zero]:.0f} GN with the casing at rest", xy=(0, R["tension_GN"][k_zero]),
                 xytext=(3900, R["tension_GN"][k_zero] * 0.78), color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
    ax1.set_title("Cable tension and bearing load through the deployment of the 700 km ring", fontsize=12)
    lines = ax1.get_legend_handles_labels(); lines2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines[0] + lines2[0], lines[1] + lines2[1], loc="upper left", fontsize=9)

    # middle: casing net weight and cable excess force
    ax2.plot(vL, R["casing_net_weight_kN_m"], color=GREEN, lw=2.0, label="Casing net weight (the bearing load)")
    ax2.plot(vL, R["cable_excess_kN_m"], color=RED, lw=2.0, label="Cable excess centrifugal force")
    if P is not None:
        ax2.plot(P["v_casing_m_s"], P["casing_net_weight_kN_m"], color=GREEN, lw=1.4, ls="--")
        ax2.plot(P["v_casing_m_s"], P["cable_excess_kN_m"], color=RED, lw=1.4, ls="--")
    ax2.set_ylabel("Force per meter of ring (kN/m)")
    ax2.axhline(0, color=GRAY, lw=0.8)
    ax2.annotate(f"{R['casing_net_weight_kN_m'][-1]:.1f} kN/m at the end", xy=(vL[-1], R["casing_net_weight_kN_m"][-1]),
                 xytext=(1900, R["casing_net_weight_kN_m"][-1] * 0.62), color=GREEN, arrowprops=dict(arrowstyle="->", color=GREEN))
    ax2.annotate(f"{R['cable_excess_kN_m'][-1]:.1f} kN/m", xy=(vL[-1], R["cable_excess_kN_m"][-1]),
                 xytext=(1500, R["cable_excess_kN_m"][-1] * 1.02), color=RED, arrowprops=dict(arrowstyle="->", color=RED))
    ax2.legend(loc="upper left", fontsize=9)
    ax2.set_xlabel("Casing speed in the direction the ring was launched (m/s): for the retrograde cable westward (W), through rest, to eastward (E); "
                   "the dashed prograde curves are drawn at the same casing slowdown (that casing stays eastward)")

    for ax in (ax1, ax2):
        ax.axvline(0, color=GRAY, lw=0.8, ls="--")
        ax.grid(True, alpha=0.3)
    ax1.text(120, ax1.get_ylim()[1] * 0.04, f"casing at rest, day {R['time_days'][k_zero]:.0f}", color=GRAY, fontsize=9, ha="right")
    ax1.set_xlim(v_o + 150, -(v_g + 150))
    ticks = [7000, 6000, 5000, 4000, 3000, 2000, 1000, 0, -500]
    ax2.set_xticks(ticks); ax2.set_xticklabels([f"{t:,} W" if t > 0 else ("0" if t == 0 else f"{-t:,} E") for t in ticks])
    plt.setp(ax1.get_xticklabels(), visible=False)

    # bottom: against time
    t = R["time_days"]
    ax3.plot(t, R["tension_GN"], color=BLUE, lw=2.0, label="Cable tension (GN)")
    ax3.set_ylabel("Cable tension (GN)", color=BLUE)
    ax3b = ax3.twinx()
    ax3b.plot(t, R["casing_net_weight_kN_m"], color=GREEN, lw=2.0, label="Casing net weight (kN/m)")
    ax3b.set_ylabel("Casing net weight (kN/m)", color=GREEN)
    ax3.set_xlabel("Time from the start of the deployment (days)")
    ax3.set_xlim(0, t[-1] * 1.02); ax3.set_ylim(0, R["tension_GN"][-1] * 1.08); ax3b.set_ylim(0, 120)
    ax3.axvline(R["time_days"][k_zero], color=GRAY, lw=0.8, ls="--")
    ax3.text(R["time_days"][k_zero] - 8, R["tension_GN"][-1] * 0.08, "casing at rest", color=GRAY, fontsize=9, ha="right")
    ax3.grid(True, alpha=0.3)
    l1 = ax3.get_legend_handles_labels(); l2 = ax3b.get_legend_handles_labels()
    ax3.legend(l1[0] + l2[0], l1[1] + l2[1], loc="upper left", fontsize=9)

    fig.subplots_adjust(left=0.09, right=0.91, top=0.96, bottom=0.06)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig03_02_deployment_tension_700.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig03_02_deployment_tension_700.svg and .png in {out}")


if __name__ == "__main__":
    main()
