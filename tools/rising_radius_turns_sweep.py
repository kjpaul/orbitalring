#!/usr/bin/env python3
"""
The number of turns at 12 MW (round 5 of 2026-11-06): the choice of Section 9.7 of the book (100, 150, 208 and 300 turns run on
the deployment simulation at 8 MW with the October cable, the fastest kept, hysteresis the second check) redone for the design
of round 4, two LIMs per site at 12 MW with the cable of the rising quadratic and the rising radius, at the full time steps.
N_TURNS = 120, 150, 180, 208, 240, 280, each a lim_simulation.py --rising-radius run through tools/rising_radius_study.py's
runner with lim_config.configure(N_TURNS=N, MAX_SITE_POWER=12e6), into runs/rising_radius_turns/<case>/.

lim_config.configure() keeps, for every N, the measured tape basis of 2026-10-27: the flat winding of two 3 mm tape layers,
the worst tape at 183.1 A, and the coil current cap I_TARGET = I_PEAK = 183.1 A at which that tape is at half its limit, a
figure measured for the 208-turn winding.  The worst tape's limit moves a little with N (the field at the tape grows with the
winding width): lim_hysteresis.coil_limit gives 185.3 A at 120 turns down to 182.1 A at 280.  Since the cap is a choice at the
other N, each is also run with its own cap (the "_own" cases, configure(I_COIL_LIMIT_MEASURED=that coil limit)).
The 208-turn case is the 12 MW run of the power sweep (runs/rising_radius_sweep/L2_N208_P12MW) copied in when it exists.
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Sections 9.7 and 9.8), Chapter 14.
The grid was extended with 135 and 165 turns to place the minimum, and the chosen N is also run at 8 and 16 MW.
Run: python tools/rising_radius_turns_sweep.py --run [--jobs=3] [--only=a,b] [--quick]
     python tools/rising_radius_turns_sweep.py --table
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_os.environ.setdefault("ORBITAL_RING_STUDY_DIR", "rising_radius_turns")
import concurrent.futures
import csv
import json
import math
import os
import shutil
import sys

import ring_altitude as ring
import lim_hysteresis as lh
import rising_radius_study as study

TURNS = [120, 135, 150, 165, 180, 208, 240, 280]      # the brief's six, plus 135 and 165 to place the minimum between 120 and 180
OWN_TURNS = [120, 150, 180, 240, 280]                  # the brief's six less 208, run also with their own cap
POWER_CHECK = {"N150_P8MW": (150, 8e6), "N150_P16MW": (150, 16e6)}   # the chosen N at 8 and 16 MW, for Section 9.8's comparison
POWER = 12e6
LAYERS, W_TAPE = 2, 3e-3
COPY_FROM = {"N208": ("rising_radius_sweep", "L2_N208_P12MW")}
FILES = ["rising_radius_lim_simulation.txt", "rising_radius_daily.csv", "rising_radius_deployment_tension.csv", "rising_radius_summary.json"]
KNEE_TOL = 0.02

OWN_CAP = {N: lh.coil_limit("flat", N, LAYERS, W_TAPE)[0] for N in TURNS}   # each winding's own coil current limit, A

CASES = {}
for N in TURNS:
    CASES[f"N{N}"] = ({"MAX_SITE_POWER": POWER, "N_TURNS": N}, {})
    if N in OWN_TURNS:
        CASES[f"N{N}_own"] = ({"MAX_SITE_POWER": POWER, "N_TURNS": N, "I_COIL_LIMIT_MEASURED": OWN_CAP[N]}, {})
for name, (N, P) in POWER_CHECK.items():
    CASES[name] = ({"MAX_SITE_POWER": P, "N_TURNS": N}, {})
study.CASES.update(CASES)


def turns_dir():
    return study.study_dir()


def run_case(name, full=True):
    src = COPY_FROM.get(name)
    if src and os.path.exists(os.path.join(ring.out_dir(), src[0], src[1], "rising_radius_summary.json")):
        out = os.path.join(turns_dir(), name)
        os.makedirs(out, exist_ok=True)
        for f in FILES:
            shutil.copy2(os.path.join(ring.out_dir(), src[0], src[1], f), os.path.join(out, f))
        return name, 0, 0.0
    return study.run_case(name, full)


def load(name):
    p = os.path.join(turns_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p))
    J["name"] = name
    N = J["N_TURNS"]
    J["tape_t_site"] = lh.tape_mass_per_lim(N, LAYERS, W_TAPE) * J["LIMS_PER_SITE"] / 1e3
    J["coil_width_m"] = lh.winding_width("flat", N, LAYERS, W_TAPE)
    J["own_cap_A"] = OWN_CAP[N]
    # the coil's inductance implied by the model's flux linkage at the end of the run, L = V / (omega I)
    f, I, V = J["final_f_supply_Hz"], J["final_current_A"], J["final_volts_V"]
    J["L_end_H"] = V / (2 * math.pi * f * I) if f > 0 and I > 0 else float("nan")
    return J


HDR = ["case", "turns", "cap A", "days", "months", "I cap last day", "P cap first day", "P cap last day", "100 kV from day", "thrust max N/LIM (day)",
       "LIM work EJ", "hyst peak kW/LIM", "cryo peak MW", "cryo end MW", "plate K", "tape t/site", "tape km/site", "coil width m",
       "L at the end mH", "volts at the end", "amps at the end", "Hz at the end", "thrust at the end N"]


def row(J):
    f = lambda x, fmt: "never" if x is None else format(x, fmt)
    return [J["name"], J["N_TURNS"], f"{J['I_TARGET']:.1f}", f"{J['t_end_days']:.1f}", f"{J['t_end_days'] / 30.4375:.1f}",
            f(J.get("current_cap_last_day"), ".1f"), f(J.get("power_cap_first_day"), ".1f"), f(J.get("power_cap_last_day"), ".1f"), f(J.get("volt_limit_day"), ".1f"),
            f"{J['thrust_max_N']:,.0f} ({J['thrust_max_day']:.1f})", f"{J['E_lim_loop_EJ']:.3f}", f"{J['p_hyst_max_kW']:.2f}", f"{J['p_cryo_max_MW']:.2f}",
            f"{J['final_p_cryo_MW']:.2f}", f"{J['plate_t_max_K']:.1f}", f"{J['tape_t_site']:.2f}", f"{J['L_HTS_SITE_m'] / 1e3:,.0f}", f"{J['coil_width_m']:.3f}",
            f"{J['L_end_H'] * 1e3:.1f}", f"{J['final_volts_V']:,.0f}", f"{J['final_current_A']:.1f}", f"{J['final_f_supply_Hz']:.1f}", f"{J['final_thrust_N']:.0f}"]


def curve(rows, own):
    pts = []
    for N in TURNS:
        name = f"N{N}" if (N == 208 or not own) else f"N{N}_own"
        J = next((r for r in rows if r["name"] == name), None)
        if J:
            pts.append((N, J))
    if not pts:
        return None
    floor = min(J["t_end_days"] for _, J in pts)
    best = next((N, J) for N, J in pts if J["t_end_days"] == floor)
    within = [N for N, J in pts if J["t_end_days"] <= floor * (1 + KNEE_TOL)]
    steps = [(a, b, (Ja["t_end_days"] - Jb["t_end_days"]) / (b - a)) for (a, Ja), (b, Jb) in zip(pts, pts[1:])]
    return dict(pts=pts, floor=floor, best=best, within=within, steps=steps)


def figure(C_fixed, C_own):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for C, col, lab in ((C_fixed, "#2b6ca3", "the 183.1 A cap of the 208-turn winding"), (C_own, "#b5451b", "each winding's own cap (182.1 to 185.3 A)")):
        if not C:
            continue
        ax.plot([N for N, _ in C["pts"]], [J["t_end_days"] for _, J in C["pts"]], "o-", color=col, lw=1.8, label=lab)
        N, J = C["best"]
        ax.plot([N], [J["t_end_days"]], "s", color=col, ms=11, mfc="none", mew=2)
    ax.set_xlabel("Turns per phase coil"); ax.set_ylabel("Deployment (days)")
    ax.set_title("Deployment days against the number of turns: 12 MW per site, two LIMs, 100 kV, rising radius", fontsize=10)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(turns_dir(), f"rising_radius_turns.{ext}"), dpi=200 if ext == "png" else None)
    plt.close(fig)


def table():
    rows = [load(n) for n in CASES]
    rows = [r for r in rows if r]
    out = [row(J) for J in rows if J["name"] not in POWER_CHECK]
    pc = [row(J) for J in rows if J["name"] in POWER_CHECK or J["name"] == "N150"]   # the chosen N at 8, 12 and 16 MW
    md = "| " + " | ".join(HDR) + " |\n|" + "---|" * len(HDR) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out) + "\n"
    if len(pc) > 1:
        md += chr(10) + "The chosen N at 8, 12 and 16 MW:" + chr(10) * 2 + "| " + " | ".join(HDR) + " |" + chr(10) + "|" + "---|" * len(HDR) + chr(10) \
              + chr(10).join("| " + " | ".join(str(c) for c in r) + " |" for r in sorted(pc, key=lambda r: float(r[3]))) + chr(10)
    print(md)
    lines = ["Each winding's own coil current limit (lim_hysteresis.coil_limit, flat, 2 layers of 3 mm tape, the worst tape at half its limit): "
             + ", ".join(f"{N} turns {OWN_CAP[N]:.1f} A" for N in TURNS)]
    C_fixed, C_own = curve(rows, False), curve(rows, True)
    for C, lab in ((C_fixed, "the 183.1 A cap"), (C_own, "each winding's own cap")):
        if not C:
            continue
        N, J = C["best"]
        lines.append(f"Days against turns with {lab}: " + "; ".join(f"{N} turns {J['t_end_days']:.1f} d" for N, J in C["pts"]))
        lines.append(f"  minimum {J['t_end_days']:.1f} days at {N} turns; within {KNEE_TOL * 100:.0f} % of it: {', '.join(str(n) for n in C['within'])} turns")
        lines.append("  days saved per turn at each step: " + "; ".join(f"{a} to {b}: {s:+.2f} d/turn" for a, b, s in C["steps"]))
    txt = "\n".join(lines)
    print(txt)
    figure(C_fixed, C_own)
    with open(os.path.join(turns_dir(), "rising_radius_turns.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n" + txt + "\n")
    with open(os.path.join(turns_dir(), "rising_radius_turns.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(HDR); wr.writerows(out)
    print(f"written rising_radius_turns.md, .csv, .svg, .png in {turns_dir()}")


def main():
    names = list(CASES)
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    if only:
        names = [n for n in only.split(",") if n in CASES]
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "3"))
    full = "--quick" not in sys.argv
    if "--run" in sys.argv:
        print(f"running {len(names)} cases, {jobs} at a time" + (" (full steps)" if full else " (quick steps)") + f", into {turns_dir()}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: run_case(n, full), names):
                print(f"  {name:12} rc={rc} {dt / 60:.1f} min" + ("  (copied)" if dt == 0.0 else ""), flush=True)
    if "--table" in sys.argv:
        table()
    if "--run" not in sys.argv and "--table" not in sys.argv:
        print(__doc__)


if __name__ == "__main__":
    main()
