#!/usr/bin/env python3
"""
The optimum deployment: a small study on the rising-radius deployment (round 2 of 2026-11-06).  Varies, one at a time from the
8 MW, 208-turn, two-LIMs-per-site baseline, the site power (8, 12, 16, 24 MW), the turns (150, 208, 260) and the LIMs per site
(2, 4, 6: one, two and three reaction plates per side of the cable), then a few combined cases; the coil voltage limit stays at
100 kV.  Each case is lim_simulation.py --rising-radius --quick (the quick time steps: 10 s on the first day, then 100 s; the
8 MW baseline gives the same 647.1 days and 23.864 EJ as the full-step run) in its own output folder under
runs/rising_radius_study/, with the parameters set through lim_config.configure().  More plates per side put more mass on
the cable, so those cases also set ORBITAL_RING_HW_KGM so that ring_altitude's quadratic sizes the cable for that hardware
(the design point the run ends at is then that ring's).
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Sections 9.5 to 9.7, Table 9.3), Chapter 3 (Section 3.2), Chapter 14.
Run: python tools/rising_radius_study.py --run [--jobs=5] [--only=a,b] [--full]     (runs the cases; --full for the 10 s steps)
     python tools/rising_radius_study.py --table                                     (the table from the summaries, to the console,
                                                                                      rising_radius_study.md and .csv in the study folder)
Environment: ORBITAL_RING_OUT (output folder, default runs); ORBITAL_RING_STUDY_DIR (the study's subfolder, default rising_radius_study)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import concurrent.futures
import csv
import json
import os
import subprocess
import sys
import time

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import ring_altitude as ring
import lim_hysteresis as lh

M_LEV = 4 * 7870.0 * 0.20 * 1.4          # the four bearing plates, kg/m (lim_config: N_LEV_PLATES x DENSITY_IRON x D_LEV x W_LEV)
M_LIM_PLATE = 3900.0 * 0.20 * 1.2        # one LIM reaction plate, kg/m (PLATE_DENSITY x T_PLATE x W_PLATE, gamma-TiAl)


def hw_kgm(n_plates_per_side):
    return M_LEV + 2 * n_plates_per_side * M_LIM_PLATE


# name: (configure overrides, environment overrides)
CASES = {
    "base_8MW":        ({}, {}),
    "P12MW":           ({"MAX_SITE_POWER": 12e6}, {}),
    "P16MW":           ({"MAX_SITE_POWER": 16e6}, {}),
    "P24MW":           ({"MAX_SITE_POWER": 24e6}, {}),
    "N150":            ({"N_TURNS": 150}, {}),
    "N260":            ({"N_TURNS": 260}, {}),
    "L4":              ({"N_PLATES_PER_SIDE": 2}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(2):.1f}"}),
    "L6":              ({"N_PLATES_PER_SIDE": 3}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(3):.1f}"}),
    "N150_16MW":       ({"N_TURNS": 150, "MAX_SITE_POWER": 16e6}, {}),
    "P16MW_L4":        ({"MAX_SITE_POWER": 16e6, "N_PLATES_PER_SIDE": 2}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(2):.1f}"}),
    "P24MW_L4":        ({"MAX_SITE_POWER": 24e6, "N_PLATES_PER_SIDE": 2}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(2):.1f}"}),
    "P24MW_L6":        ({"MAX_SITE_POWER": 24e6, "N_PLATES_PER_SIDE": 3}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(3):.1f}"}),
    "P16MW_N150_L4":   ({"MAX_SITE_POWER": 16e6, "N_TURNS": 150, "N_PLATES_PER_SIDE": 2}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(2):.1f}"}),
    "P24MW_N150_L6":   ({"MAX_SITE_POWER": 24e6, "N_TURNS": 150, "N_PLATES_PER_SIDE": 3}, {"ORBITAL_RING_HW_KGM": f"{hw_kgm(3):.1f}"}),
}
ORDER = list(CASES)


def study_dir():
    return ring.out_dir(os.environ.get("ORBITAL_RING_STUDY_DIR", "rising_radius_study"))


def run_case(name, full=False):
    cfg_over, env_over = CASES[name]
    out = os.path.join(study_dir(), name)
    os.makedirs(out, exist_ok=True)
    args = "['lim_simulation.py', '--rising-radius'" + ("" if full else ", '--quick'") + "]"
    code = (f"import lim_config as c; c.configure(**{cfg_over!r}); import runpy, sys; sys.argv={args}; "
            f"runpy.run_path('lim_simulation.py', run_name='__main__')")
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_OUT=out, **env_over)
    log = os.path.join(out, "rising_radius_lim_simulation.txt")
    t0 = time.time()
    with open(log, "w", encoding="utf-8", errors="replace") as lf:
        r = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    return name, r.returncode, time.time() - t0


def load(name):
    p = os.path.join(study_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p))
    J["name"] = name
    cfg_over, env_over = CASES[name]
    J["tape_mass_per_site_t"] = lh.tape_mass_per_lim(J["N_TURNS"], 2, 3e-3) * J["LIMS_PER_SITE"] / 1e3
    return J


def table():
    base = load("base_8MW")
    rows = [load(n) for n in ORDER]
    rows = [r for r in rows if r]
    hdr = ["case", "MW/site", "turns", "LIMs/site", "days", "phase 1 d", "phase 2 d", "100 kV from day", "thrust/LIM at 100 kV N", "days at I cap", "days at P cap",
           "LIM work EJ", "peak site MW", "peak plate K", "cryo peak MW", "cryo end MW", "plates kg/m cable", "cable kg/m", "tape km/site", "tape t/site", "coil area m2/site"]
    out_rows = []
    for J in rows:
        p2 = J.get("phase2_duration_days") or 0.0
        out_rows.append([J["name"], f"{J['MAX_SITE_POWER'] / 1e6:.0f}", J["N_TURNS"], J["LIMS_PER_SITE"], f"{J['t_end_days']:.1f}", f"{(J.get('t_phase2_days') or J['t_end_days']):.1f}",
                         f"{p2:.1f}", ("never" if J.get("volt_limit_day") is None else f"{J['volt_limit_day']:.1f}"),
                         ("" if J.get("thrust_at_volt_limit_N") is None else f"{J['thrust_at_volt_limit_N']:.0f}"),
                         f"{J.get('days_at_current_cap', 0):.1f}", f"{J.get('days_at_power_cap', 0):.1f}", f"{J['E_lim_loop_EJ']:.3f}", f"{J['site_power_max_MW']:.2f}",
                         f"{J['plate_t_max_K']:.1f}", f"{J['p_cryo_max_MW']:.2f}", f"{J['final_p_cryo_MW']:.3f}", f"{J['M_HARDWARE_kg_m'] - M_LEV:.0f}",
                         f"{J['M_CABLE_M_kg_m']:,.0f}", f"{J['L_HTS_SITE_m'] / 1e3:.0f}", f"{J['tape_mass_per_site_t']:.1f}", f"{J['COIL_SURFACE_AREA_PER_SITE_m2']:.0f}"])
    md = "| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out_rows) + "\n"
    print(md)
    if base:
        print(f"Baseline: {base['t_end_days']:.1f} days, {base['E_lim_loop_EJ']:.3f} EJ, cable {base['M_CABLE_M_kg_m']:,.0f} kg/m, tape {base['tape_mass_per_site_t']:.1f} t per site.")
        print("Days saved against the baseline, per unit of hardware:")
        for J in rows:
            if J["name"] == "base_8MW":
                continue
            saved = base["t_end_days"] - J["t_end_days"]
            d_tape = J["tape_mass_per_site_t"] - base["tape_mass_per_site_t"]
            d_cable = J["M_CABLE_M_kg_m"] - base["M_CABLE_M_kg_m"]
            d_power = (J["MAX_SITE_POWER"] - base["MAX_SITE_POWER"]) / 1e6
            print(f"  {J['name']:16} {saved:7.1f} days saved; site power {d_power:+.0f} MW, tape {d_tape:+.1f} t/site, cable {d_cable:+,.0f} kg/m, LIM energy {J['E_lim_loop_EJ'] - base['E_lim_loop_EJ']:+.3f} EJ"
                  + (f"; {saved / d_power:.1f} days per MW" if d_power else "") + (f"; {saved / d_tape:.1f} days per tonne of tape per site" if d_tape else ""))
    figure(rows, base)
    with open(os.path.join(study_dir(), "rising_radius_study.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    with open(os.path.join(study_dir(), "rising_radius_study.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(hdr); wr.writerows(out_rows)
    print(f"written rising_radius_study.md and .csv in {study_dir()}")


def figure(rows, base):
    """rising_radius_study.png/svg: the deployment days of every case split into phase 1 and phase 2, the tape per site beside them."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
    names = [J["name"] for J in rows]
    p1 = np.array([(J.get("t_phase2_days") or J["t_end_days"]) for J in rows]); p2 = np.array([J.get("phase2_duration_days") or 0.0 for J in rows])
    tape = np.array([J["tape_mass_per_site_t"] for J in rows]); x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(11, 5.6))
    ax.bar(x, p1, color="#2b6ca3", label="Phase 1: the casing to ground speed")
    ax.bar(x, p2, bottom=p1, color="#b5451b", label="Phase 2: the final tensioning against the Earth")
    for k, J in enumerate(rows):
        ax.text(k, p1[k] + p2[k] + 8, f"{p1[k] + p2[k]:.0f} d", ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Deployment time (days)"); ax.set_ylim(0, (p1 + p2).max() * 1.15)
    axb = ax.twinx(); axb.plot(x, tape, "o--", color="#2e8b57", lw=1.2, ms=5, label="HTS tape per site (t)"); axb.set_ylabel("HTS tape per site (t)", color="#2e8b57"); axb.set_ylim(0, tape.max() * 1.6)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_title("The deployment study: site power, turns and LIMs per site, one at a time and combined (100 kV coil limit throughout)", fontsize=11)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(study_dir(), f"rising_radius_study.{ext}"), **({"dpi": 200} if ext == "png" else {}))
    plt.close(fig)
    print("written rising_radius_study.svg and .png")


def main():
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "5"))
    only = next((a[7:].split(",") for a in sys.argv if a.startswith("--only=")), None)
    full = "--full" in sys.argv
    if "--run" in sys.argv:
        names = [n for n in ORDER if only is None or n in only]
        print(f"running {len(names)} cases, {jobs} at a time" + (" (full steps)" if full else " (quick steps)"), flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: run_case(n, full), names):
                print(f"  {name:16} exit {rc} in {dt / 60:.1f} min", flush=True)
    if "--table" in sys.argv or "--run" in sys.argv:
        table()
    if not any(a in sys.argv for a in ("--run", "--table")):
        print(__doc__)


if __name__ == "__main__":
    main()
