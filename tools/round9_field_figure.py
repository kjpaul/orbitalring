#!/usr/bin/env python3
"""
The figure of the field at the plate, redrawn (round 9, item 5b; addendum A4): the size of the field THROUGH the face only (the part that pulls),
the plate's two edges marked, the bundles under the axis to scale, the even field as a dashed line.  The old version (the size of the whole field)
stays beside it as fig07_candidate_field_at_plate.
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run (from Python/orbitalring): python tools/round9_field_figure.py [output folder]     (default ../../Vol_III/figure_work_2026-11-19)
Round 22 (2026-11-19): the compensation panel is the state "fast 199" (a lift compensation track in the design launch at 30 km/s, the book's
199 kN/m) when the face solve holds it, else "fast 196.5" as on 16 November; its label is the state's own pull.  The states of the figure are
printed beside the book's fields, and the log goes to runs/round22_logs/round22_field_figure.json and .txt as well.
Reads runs/gap_trade_logs/gap_trade_face.json (bearing_gap_trade.py --part=solve-face, 2.5 mm cells).  Prints, for each state, the pull from the
integral of (Bn^2 - Bt^2) / (2 mu0) across the face against the solver's pull, and the share of it from the strips outside the outer bundles;
writes runs/round9_logs/round9_field_figure.json.
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

import gap_trade_figures as g

C = g.C; FS = g.FS; LW, LW2 = g.LW, g.LW2
MU0 = 4e-7 * math.pi


def shares(d):
    out = {}
    print(f"  {'state':>16} {'solver pull kN/m':>16} {'integral over the face':>22} {'difference %':>12} {'Bn alone':>9} {'Bt takes away %':>15} {'outside the outer bundles %':>27} {'outside their centers %':>23} {'Bn at the edge T':>16} {'rms Bn T':>8} {'even T':>7}")
    for k, s in d.items():
        x = np.array(s["line"]["x"]); Bn = np.array(s["line"]["By"]); Bt = np.array(s["line"]["Bx"]); h = x[1] - x[0]
        p = (Bn**2 - Bt**2) / (2 * MU0); F = p.sum() * h; Fn = (Bn**2 / (2 * MU0)).sum() * h
        xo = max(abs(l[0]) for l in s["legs"]); edge = xo + s["bw"] / 2
        xc = x + h / 2
        f_out = p[np.abs(xc) > edge].sum() * h / F; f_cen = p[np.abs(xc) > xo].sum() * h / F
        out[k] = dict(F=s["F"], F_face=F, diff=F / s["F"] - 1, F_n=Fn, Bt_takes=1 - F / Fn, share_outside=f_out, share_outside_centers=f_cen, Bn_edge=float(max(abs(Bn[0]), abs(Bn[-1]))),
                      rms_n=float(np.sqrt((Bn**2).mean())), even=s["line"]["even"], strip_mm=(s["Wp"] / 2 - edge) * 1e3)
        o = out[k]
        print(f"  {k:>16} {s['F']/1e3:16.1f} {F/1e3:22.1f} {o['diff']*100:12.1f} {Fn/1e3:9.1f} {o['Bt_takes']*100:15.1f} {f_out*100:27.1f} {f_cen*100:23.1f} {o['Bn_edge']:16.2f} {o['rms_n']:8.3f} {o['even']:7.3f}")
    return out


BOOK = {"support cap": ("0.72 T", "1.13 T"), "support normal": ("0.32 T", "not printed"), "fast 199": ("0.71 T", "not printed")}


def fig(folder, d):
    comp = "fast 199" if "fast 199" in d else "fast 196.5"
    comp_label = f"{d[comp]['F']/1e3:.0f} kN/m" if comp == "fast 199" else "196.5 kN/m"
    f, axes = plt.subplots(2, 1, figsize=(g.W_IN, 5.0))
    f.subplots_adjust(left=0.125, right=0.975, top=0.985, bottom=0.085, hspace=0.30)
    panels = [(axes[0], "Support plate", [("support cap", C["support"], LW2), ("support normal", C["support"], LW2 * 0.6)], 1.45),
              (axes[1], "Compensation plate", [(comp, C["load"], LW2)], 1.45)]
    for ax, title, curves, ymax in panels:
        s0 = d[curves[0][0]]; Wp = s0["Wp"]
        yb = -0.20 * ymax / 1.45
        for key, col, lw in curves:
            s = d[key]; x = np.array(s["line"]["x"]); B = np.abs(np.array(s["line"]["By"]))
            ax.plot(x * 1e3, B, color=col, lw=lw, zorder=4, solid_capstyle="round")
            ev = s["line"]["even"]
            ax.plot([-Wp / 2 * 1e3, Wp / 2 * 1e3], [ev, ev], color=col, lw=LW, ls=(0, (4, 2.5)), zorder=3)
        ax.axhline(0, color=C["ink"], lw=LW, zorder=2)
        ax.set_xlim(-Wp / 2 * 1e3 - 30, Wp / 2 * 1e3 + 255); ax.set_ylim(yb * 1.25, ymax)
        for xe in (-Wp / 2 * 1e3, Wp / 2 * 1e3):
            ax.plot([xe, xe], [0, ymax * 0.86], color=C["muted"], lw=LW * 0.8, zorder=1)
        f.canvas.draw()
        bb = ax.get_window_extent(); x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        T_per_mm = (y1 - y0) / bb.height * bb.width / (x1 - x0)
        h_b = s0["bh"] * 1e3 * T_per_mm
        for xl, sg in s0["legs"]:
            ax.add_patch(Rectangle(((xl - s0["bw"] / 2) * 1e3, -0.035 - h_b), s0["bw"] * 1e3, h_b, fc=curves[0][1], ec="none", zorder=4))
        ax.set_ylabel("Field through the face (T)", fontsize=FS, color=C["ink"])
        ax.set_yticks([0, 0.5, 1.0])
        ax.tick_params(labelsize=FS, colors=C["ink"], width=LW, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(LW); ax.spines[sp].set_color(C["ink"])
        ax.spines["left"].set_bounds(0, ymax); ax.spines["bottom"].set_bounds(-Wp / 2 * 1e3, Wp / 2 * 1e3)
        ax.set_xticks([t for t in range(-600, 601, 200) if abs(t) <= Wp / 2 * 1e3 + 1])
        xr = Wp / 2 * 1e3 + 18
        for key, col, lw in curves:
            ax.text(xr, d[key]["line"]["even"], f"even field\n{d[key]['line']['even']:.2f} T", fontsize=FS, color=col, ha="left", va="center", linespacing=1.05)
        ax.text(xr, -0.035 - h_b / 2, "bundles", fontsize=FS, color=curves[0][1], ha="left", va="center")
        ax.text(xr, ymax * 0.80, "edge of\nthe plate", fontsize=FS, color=C["muted"], ha="left", va="center", linespacing=1.05)
        ax.text(x0 + 8, ymax * 0.995, title, fontsize=FS, fontweight="bold", color=C["ink"], ha="left", va="top")
    axes[1].set_xlabel("Position across the plate (mm)", fontsize=FS, color=C["ink"], x=0.40)
    axes[0].text(-300, 1.20, "the cap", fontsize=FS, color=C["support"], ha="center", va="bottom")
    axes[0].text(300, 0.47, "normal\noperation", fontsize=FS, color=C["support"], ha="center", va="bottom", linespacing=1.05)
    axes[1].text(0, 0.98, comp_label, fontsize=FS, color=C["load"], ha="center", va="bottom")
    return g.save(f, folder, "fig07_candidate_field_through_face")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "..", "Vol_III", "figure_work_2026-11-19")
    d = json.load(open(os.path.join(g.LOGS, "gap_trade_face.json"), encoding="utf-8"))
    out = shares(d)
    os.makedirs(os.path.join(os.path.dirname(g.LOGS), "round9_logs"), exist_ok=True)
    json.dump(out, open(os.path.join(os.path.dirname(g.LOGS), "round9_logs", "round9_field_figure.json"), "w", encoding="utf-8"), indent=1)
    comp = "fast 199" if "fast 199" in d else "fast 196.5"
    print()
    print("  The states the figure shows, beside the book (Table 7.3 and the paragraph under it):")
    print(f"  {'state':>16} {'pull kN/m':>10} {'ampere-turns kA':>15} {'plate m':>8} {'even field T':>12} {'book':>8} {'peak through the face T':>24} {'peak of the whole field T':>26} {'book':>12} {'in the plate T':>14}")
    for k in ("support cap", "support normal", comp):
        s = d[k]; b = BOOK.get(k, ("", ""))
        out[k].update(NI=s["NI"], Wp=s["Wp"], peak_n=s["line"]["peak_n"], peak=s["line"]["peak"], B_plate=s["B_plate"], legs=s["legs"], bw=s["bw"], bh=s["bh"], y_line=s["line"]["y"], h=s["h"])
        print(f"  {k:>16} {s['F']/1e3:10.2f} {s['NI']/1e3:15.1f} {s['Wp']:8.1f} {s['line']['even']:12.4f} {b[0]:>8} {s['line']['peak_n']:24.4f} {s['line']['peak']:26.4f} {b[1]:>12} {s['B_plate']:14.2f}")
    r22 = os.path.join(os.path.dirname(g.LOGS), "round22_logs"); os.makedirs(r22, exist_ok=True)
    json.dump(out, open(os.path.join(r22, "round22_field_figure.json"), "w", encoding="utf-8"), indent=1)
    os.makedirs(folder, exist_ok=True)
    fig(folder, d)


if __name__ == "__main__":
    main()
