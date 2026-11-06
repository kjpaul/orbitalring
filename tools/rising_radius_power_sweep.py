#!/usr/bin/env python3
"""
The site power sweep on the rising-radius deployment with the cable of the rising quadratic (round 4 of 2026-11-06):
two LIMs per site at 208 turns for 6 to 24 MW per site, four LIMs at 208 turns for 8 to 32 MW, and four LIMs at 150 turns
for 12 to 32 MW, every case lim_simulation.py --rising-radius at the full time steps (1 s on the first day, then 10 s)
through tools/rising_radius_study.py's runner (the parameters through lim_config.configure, the four-LIM cases with
ORBITAL_RING_HW_KGM set to their plate mass so that ring_altitude's quadratic sizes the cable for them), into
runs/rising_radius_sweep/<case>/.  The table (--table) gives, per case, the deployment days, the day the site power cap
stops binding and the day the 100 kV coil limit arrives, the LIM work, the peak plate temperature, the peak and end
cryogenic power per site, the HTS tape per site, the deployment panel width (55 m at 8 MW, scaled with the peak site
power) and the cable of the quadratic; then, per LIM count, the curve of days against site power with the knee (the
lowest power within 2 % of the curve's floor) and the days saved per MW at each step.
Book: Orbital Ring Engineering, Volume I: Chapter 9 (the section on shortening the deployment), Chapter 8 (the panel), Chapter 14.
Run: python tools/rising_radius_power_sweep.py --run [--jobs=3] [--only=a,b] [--quick]
     python tools/rising_radius_power_sweep.py --table
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_os.environ.setdefault("ORBITAL_RING_STUDY_DIR", "rising_radius_sweep")
import concurrent.futures
import csv
import json
import os
import sys

import ring_altitude as ring
import lim_hysteresis as lh
import rising_radius_study as study

SWEEP = {2: (208, [6, 8, 10, 12, 14, 16, 20, 24]), 4: (208, [8, 12, 16, 20, 24, 32]), "4N150": (150, [12, 16, 24, 32])}
PANEL_WIDTH_8MW = 55.0           # m, the deployment panel per site at 8 MW (Chapter 8: 291 W/m^2 averaged around the ring)
KNEE_TOL = 0.02                  # the knee: the lowest power within 2 % of the curve's floor


def case_name(lims, turns, mw):
    return f"L{lims}_N{turns}_P{mw}MW"


CASES = {}
for key, (turns, powers) in SWEEP.items():
    lims = 4 if key == "4N150" else key
    for mw in powers:
        over = {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": turns}
        env = {}
        if lims == 4:
            over["N_PLATES_PER_SIDE"] = 2
            env["ORBITAL_RING_HW_KGM"] = f"{study.hw_kgm(2):.1f}"
        CASES[case_name(lims, turns, mw)] = (over, env)
study.CASES.update(CASES)        # the study's runner and loader know the sweep's cases


def sweep_dir():
    return study.study_dir()


def load(name):
    p = os.path.join(sweep_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p))
    J["name"] = name
    J["tape_t_site"] = lh.tape_mass_per_lim(J["N_TURNS"], 2, 3e-3) * J["LIMS_PER_SITE"] / 1e3
    J["design"] = ring.design(m_hw=J["M_HARDWARE_kg_m"])
    return J


def row(J):
    f = lambda x, fmt: "never" if x is None else format(x, fmt)
    d = J["design"]; ps = J["site_power_max_MW"]
    return [J["name"], J["LIMS_PER_SITE"], J["N_TURNS"], f"{J['MAX_SITE_POWER'] / 1e6:.0f}", f"{J['t_end_days']:.1f}", f"{J['t_end_days'] / 30.4375:.1f}",
            f"{(J.get('phase2_duration_days') or 0):.1f}", f(J.get("power_cap_last_day"), ".1f"), f(J.get("volt_limit_day"), ".1f"),
            f"{J['E_lim_loop_EJ']:.3f}", f"{ps:.2f}", f"{PANEL_WIDTH_8MW * ps / 8:.0f}", f"{J['plate_t_max_K']:.1f}", f"{J['p_cryo_max_MW']:.2f}",
            f"{J['final_p_cryo_MW']:.2f}", f"{J['tape_t_site']:.1f}", f"{J['M_HARDWARE_kg_m']:,.0f}", f"{d['m_struct']:,.0f}", f"{d['m_cable']:,.0f}",
            f"{d['m_cable'] * d['L_ring'] / 1e9:,.0f}", f"{d['v_cable']:,.1f}", f"{d['tension'] / 1e9:.1f}", f"{J['final_thrust_N']:.0f}", f"{J['final_site_power_MW']:.2f}"]


HDR = ["case", "LIMs/site", "turns", "MW/site", "days", "months", "phase 2 d", "P cap last day", "100 kV from day", "LIM work EJ", "peak site MW",
       "panel width m", "plate K", "cryo peak MW", "cryo end MW", "tape t/site", "hardware kg/m", "structural kg/m", "cable kg/m", "cable Mt",
       "v_cable m/s", "tension GN", "end thrust N/LIM", "end site MW"]


def knees(rows):
    """Per LIM count and turns: the days against power, the days saved per MW at each step, the knee."""
    out = {}
    for key, (turns, powers) in SWEEP.items():
        lims = 4 if key == "4N150" else key
        pts = [(mw, J) for mw in powers for J in [next((r for r in rows if r["name"] == case_name(lims, turns, mw)), None)] if J]
        if not pts:
            continue
        floor = min(J["t_end_days"] for _, J in pts)
        knee = next(((mw, J) for mw, J in pts if J["t_end_days"] <= floor * (1 + KNEE_TOL)), None)
        steps = []
        for (p0, J0), (p1, J1) in zip(pts, pts[1:]):
            steps.append((p0, p1, (J0["t_end_days"] - J1["t_end_days"]) / (p1 - p0)))
        out[key] = dict(lims=lims, turns=turns, pts=pts, floor=floor, knee=knee, steps=steps)
    return out


def figure(rows, K):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for key, col in ((2, "#2b6ca3"), (4, "#b5451b"), ("4N150", "#2e8b57")):
        if key not in K:
            continue
        k = K[key]
        x = [mw for mw, _ in k["pts"]]; y = [J["t_end_days"] for _, J in k["pts"]]
        ax.plot(x, y, "o-", color=col, lw=1.8, label=f"{k['lims']} LIMs per site, {k['turns']} turns")
        if k["knee"]:
            mw, J = k["knee"]
            ax.plot([mw], [J["t_end_days"]], "s", color=col, ms=11, mfc="none", mew=2)
            ax.annotate(f"knee {mw} MW, {J['t_end_days']:.0f} d", (mw, J["t_end_days"]), xytext=(8, 10), textcoords="offset points", color=col, fontsize=9)
    ax.set_xlabel("Site power cap (MW)"); ax.set_ylabel("Deployment (days)")
    ax.set_title("Deployment days against site power, rising radius, 100 kV; the cable of the rising quadratic", fontsize=10.5)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(sweep_dir(), f"rising_radius_power_sweep.{ext}"), dpi=200 if ext == "png" else None)
    plt.close(fig)


def table():
    rows = [load(n) for n in CASES]
    rows = [r for r in rows if r]
    out = [row(J) for J in rows]
    md = "| " + " | ".join(HDR) + " |\n|" + "---|" * len(HDR) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out) + "\n"
    print(md)
    K = knees(rows)
    lines = []
    for key, k in K.items():
        lines.append(f"{k['lims']} LIMs per site, {k['turns']} turns: floor {k['floor']:.1f} days; "
                     + "; ".join(f"{mw} MW {J['t_end_days']:.1f} d" for mw, J in k["pts"]))
        lines.append("  days saved per MW at each step: " + "; ".join(f"{p0} to {p1} MW {s:+.1f} d/MW" for p0, p1, s in k["steps"]))
        if k["knee"]:
            mw, J = k["knee"]
            lines.append(f"  knee (lowest power within {KNEE_TOL * 100:.0f} % of the floor): {mw} MW, {J['t_end_days']:.1f} days, {J['t_end_days'] / k['floor'] * 100 - 100:+.2f} % above the floor")
    txt = "\n".join(lines)
    print(txt)
    figure(rows, K)
    with open(os.path.join(sweep_dir(), "rising_radius_power_sweep.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n" + txt + "\n")
    with open(os.path.join(sweep_dir(), "rising_radius_power_sweep.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(HDR); wr.writerows(out)
    print(f"written rising_radius_power_sweep.md, .csv, .svg, .png in {sweep_dir()}")


def main():
    names = list(CASES)
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    if only:
        names = [n for n in only.split(",") if n in CASES]
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "3"))
    full = "--quick" not in sys.argv
    if "--run" in sys.argv:
        print(f"running {len(names)} cases, {jobs} at a time" + (" (full steps)" if full else " (quick steps)") + f", into {sweep_dir()}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: study.run_case(n, full), names):
                print(f"  {name:16} rc={rc} {dt / 60:.1f} min", flush=True)
    if "--table" in sys.argv:
        table()
    if "--run" not in sys.argv and "--table" not in sys.argv:
        print(__doc__)


if __name__ == "__main__":
    main()
