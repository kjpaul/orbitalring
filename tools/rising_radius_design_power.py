#!/usr/bin/env python3
"""
The reruns Section 9.8 and Table 9.3 of the book quote, with the corrected slip rule (round 6a of 2026-11-07, SLIP_RULE
"both": the controller lowers the slip at the coil voltage limit as it does at the site power cap): the design's own power
curve, 150 turns and two LIMs per site at 8, 12, 16, 20 and 24 MW; the book's old coil, 208 turns and two LIMs, at 12 and
16 MW; and the fast options, four LIMs per site with 150 and with 208 turns at 16 and 24 MW.  Every case is a full-step
lim_simulation.py --rising-radius run with the cable of the rising quadratic, 100 kV, through tools/rising_radius_study.py's
runner (the four-LIM cases with ORBITAL_RING_HW_KGM set to their plate mass), into runs/rising_radius_design/<case>/.
The table has the columns of Section 14.2 (tools/rising_radius_turns_sweep.py), and the days-against-power curve of the
design with its knee (the lowest power within 2 % of the curve's floor).
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Sections 9.7, 9.8, Table 9.3), Chapter 8, Chapter 14.
Run: python tools/rising_radius_design_power.py --run [--jobs=4] [--only=a,b] [--quick]
     python tools/rising_radius_design_power.py --table
Environment: ORBITAL_RING_OUT (output folder, default runs); ORBITAL_RING_SLIP_RULE (both, the default, or power)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_os.environ.setdefault("ORBITAL_RING_STUDY_DIR", "rising_radius_design")
import concurrent.futures
import csv
import os
import sys

import lim_config as cfg
import rising_radius_study as study
import rising_radius_turns_sweep as turns

GRID = {"L2_N150": (2, 150, [8, 12, 16, 20, 24]), "L2_N208": (2, 208, [12, 16]), "L4_N150": (4, 150, [16, 24]), "L4_N208": (4, 208, [16, 24])}
KNEE_TOL = 0.02

CASES = {}
for key, (lims, N, powers) in GRID.items():
    for mw in powers:
        over = {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": N}
        env = {}
        if lims == 4:
            over["N_PLATES_PER_SIDE"] = 2
            env["ORBITAL_RING_HW_KGM"] = f"{study.hw_kgm(2):.1f}"
        CASES[f"{key}_P{mw}MW"] = (over, env)
study.CASES.update(CASES)


def design_dir():
    return study.study_dir()


def load(name):
    p = os.path.join(design_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    import json, math
    import lim_hysteresis as lh
    J = json.load(open(p)); J["name"] = name
    N = J["N_TURNS"]
    J["tape_t_site"] = lh.tape_mass_per_lim(N, 2, 3e-3) * J["LIMS_PER_SITE"] / 1e3
    J["coil_width_m"] = lh.winding_width("flat", N, 2, 3e-3)
    f, I, V = J["final_f_supply_Hz"], J["final_current_A"], J["final_volts_V"]
    J["L_end_H"] = V / (2 * math.pi * f * I) if f > 0 and I > 0 else float("nan")
    return J


HDR = ["case", "LIMs", "turns", "MW/site"] + turns.HDR[3:] + ["slip at the end %", "slip rule", "cable kg/m", "cable Mt"]


def row(J):
    r = turns.row(J)
    import ring_altitude as ring
    d = ring.design(m_hw=J["M_HARDWARE_kg_m"])
    return [J["name"], J["LIMS_PER_SITE"], J["N_TURNS"], f"{J['MAX_SITE_POWER'] / 1e6:.0f}"] + r[3:] + [f"{J['final_slip_pct']:.2f}", J.get("SLIP_RULE", "power"), f"{d['m_cable']:,.0f}", f"{d['m_cable'] * d['L_ring'] / 1e9:,.0f}"]


def table():
    rows = [load(n) for n in CASES]
    rows = [r for r in rows if r]
    out = [row(J) for J in rows]
    md = "| " + " | ".join(HDR) + " |\n|" + "---|" * len(HDR) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out) + "\n"
    print(md)
    lines = []
    for key, (lims, N, powers) in GRID.items():
        pts = [(mw, J) for mw in powers for J in [next((r for r in rows if r["name"] == f"{key}_P{mw}MW"), None)] if J]
        if len(pts) < 2:
            continue
        floor = min(J["t_end_days"] for _, J in pts)
        knee = next((mw, J) for mw, J in pts if J["t_end_days"] <= floor * (1 + KNEE_TOL))
        lines.append(f"{lims} LIMs, {N} turns: " + "; ".join(f"{mw} MW {J['t_end_days']:.1f} d" for mw, J in pts)
                     + f"; floor {floor:.1f} d; knee {knee[0]} MW, {knee[1]['t_end_days']:.1f} d; days saved per MW: "
                     + "; ".join(f"{a} to {b} MW {(Ja['t_end_days'] - Jb['t_end_days']) / (b - a):+.1f}" for (a, Ja), (b, Jb) in zip(pts, pts[1:])))
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(design_dir(), "rising_radius_design_power.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n" + txt + "\n")
    with open(os.path.join(design_dir(), "rising_radius_design_power.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(HDR); wr.writerows(out)
    # the figure: the design's curve, with the 208-turn and four-LIM points
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for key, col, mk in (("L2_N150", "#2b6ca3", "o-"), ("L2_N208", "#888888", "s--"), ("L4_N150", "#2e8b57", "^-"), ("L4_N208", "#b5451b", "v--")):
        lims, N, powers = GRID[key]
        pts = [(mw, J) for mw in powers for J in [next((r for r in rows if r["name"] == f"{key}_P{mw}MW"), None)] if J]
        if pts:
            ax.plot([mw for mw, _ in pts], [J["t_end_days"] for _, J in pts], mk, color=col, lw=1.8, label=f"{lims} LIMs per site, {N} turns")
    ax.set_xlabel("Site power cap (MW)"); ax.set_ylabel("Deployment (days)")
    ax.set_title("Deployment days against site power with the corrected slip rule; the cable of the rising quadratic, 100 kV", fontsize=9.5)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(design_dir(), f"rising_radius_design_power.{ext}"), dpi=200 if ext == "png" else None)
    plt.close(fig)
    print(f"written rising_radius_design_power.md, .csv, .svg, .png in {design_dir()}")


def main():
    names = list(CASES)
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    if only:
        names = [n for n in only.split(",") if n in CASES]
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "4"))
    full = "--quick" not in sys.argv
    if "--run" in sys.argv:
        print(f"running {len(names)} cases, {jobs} at a time" + (" (full steps)" if full else " (quick steps)") + f", slip rule {cfg.SLIP_RULE}, into {design_dir()}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: study.run_case(n, full), names):
                print(f"  {name:16} rc={rc} {dt / 60:.1f} min", flush=True)
    if "--table" in sys.argv:
        table()
    if "--run" not in sys.argv and "--table" not in sys.argv:
        print(__doc__)


if __name__ == "__main__":
    main()
