#!/usr/bin/env python3
"""
The turns sweep rerun with the corrected slip rule (2026-11-14, Vol_III/Claude_Code_instructions_turns_sweep_rerun_2026-11-14.md).
Round 5's sweep (docs/deployment_rising_radius_2026-11.md, Section 14) ran before round 6a corrected the slip rule (SLIP_RULE = both,
Section 15.1): a coil at the 100 kV limit now goes to the reduced slip as a site at its power cap does.  This reruns the sweep under
the corrected rule, at the full time steps, two LIMs per site, 100 kV, the rising quadratic and the rising radius, the 183.1 A cap as
configure() leaves it: 12 MW at N = 120 to 280, 8 MW at N = 150 to 250, 16 MW at N = 135, 150, 165; the cases that already exist
under the corrected rule (tools/rising_radius_design_power.py: 150 turns at 8, 12 and 16 MW, 208 at 12 MW) are copied in.  One more
run gives Section 9.7's hysteresis number: the design coil (150 turns, 12 MW) with HYSTERESIS_MODEL = "book", the loss factor of the
2025 book, in place of the Brandt-Indenbom model.
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Section 9.7, Table 9.3 and 9.4, a candidate Figure 9.3), Chapter 14.
Run: python tools/rising_radius_turns_rerun.py --run [--jobs=6] [--only=a,b] [--quick]
     python tools/rising_radius_turns_rerun.py --table          (the table, the curves, the checks, the figure, into the study folder)
Environment: ORBITAL_RING_OUT (default runs); ORBITAL_RING_SLIP_RULE must be 'both' (the default) for the runs to mean anything.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_os.environ.setdefault("ORBITAL_RING_STUDY_DIR", "rising_radius_turns_rerun")
import concurrent.futures
import csv
import json
import math
import os
import shutil
import sys

import ring_altitude as ring
import lim_config as cfg
import lim_hysteresis as lh
import lim_physics as phys
import rising_radius_study as study

TURNS_12 = [120, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 190, 208, 240, 280]
TURNS_8 = [150, 170, 190, 208, 230, 250]
TURNS_16 = [135, 150, 165]
LAYERS, W_TAPE = 2, 3e-3
COP_BOOK = 0.0106                                  # Chapter 14's COP: cryogenic plant = cold load / COP
KNEE_TOL = 0.02
FILES = ["rising_radius_lim_simulation.txt", "rising_radius_daily.csv", "rising_radius_deployment_tension.csv", "rising_radius_summary.json"]
COPY_FROM = {"N150_P12MW": ("rising_radius_design", "L2_N150_P12MW"), "N208_P12MW": ("rising_radius_design", "L2_N208_P12MW"),
             "N150_P8MW": ("rising_radius_design", "L2_N150_P8MW"), "N150_P16MW": ("rising_radius_design", "L2_N150_P16MW")}

CASES = {}
for P, Ns in ((12, TURNS_12), (8, TURNS_8), (16, TURNS_16)):
    for N in Ns:
        CASES[f"N{N}_P{P}MW"] = ({"MAX_SITE_POWER": P * 1e6, "N_TURNS": N}, {})
CASES["N150_P12MW_bookhyst"] = ({"MAX_SITE_POWER": 12e6, "N_TURNS": 150, "HYSTERESIS_MODEL": "book"}, {})
study.CASES.update(CASES)


def out_dir():
    return study.study_dir()


def run_case(name, full=True):
    src = COPY_FROM.get(name)
    if src and os.path.exists(os.path.join(ring.out_dir(), src[0], src[1], "rising_radius_summary.json")):
        out = os.path.join(out_dir(), name)
        os.makedirs(out, exist_ok=True)
        for f in FILES:
            shutil.copy2(os.path.join(ring.out_dir(), src[0], src[1], f), os.path.join(out, f))
        return name, 0, 0.0
    return study.run_case(name, full)


def load(name):
    p = os.path.join(out_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p)); J["name"] = name
    N = J["N_TURNS"]
    J["tape_t_site"] = lh.tape_mass_per_lim(N, LAYERS, W_TAPE) * J["LIMS_PER_SITE"] / 1e3
    J["P_MW"] = J["MAX_SITE_POWER"] / 1e6
    J["copied"] = name in COPY_FROM
    vl = J.get("volt_limit_day")
    J["days_at_V"] = (J["t_end_days"] - vl) if vl is not None else 0.0
    pf, pl = J.get("power_cap_first_day"), J.get("power_cap_last_day")
    J["days_P"] = (pl - pf) if (pf is not None and pl is not None) else 0.0
    J["cryo_peak_cop_MW"] = J["q_cryo_cold_max_kW"] / COP_BOOK / 1e3 if J.get("q_cryo_cold_max_kW") else float("nan")
    # the cold load at the end: the simulation's cryogenic power times its own COP (lim_physics.calc_cryo_power: p_cryo = q_cold / COP_sim, COP_sim =
    # calc_cryo_cop(77.4 K, 400 K, CRYO_EFF 0.05)), then over the book's COP
    J["cop_sim"] = phys.calc_cryo_cop(cfg.T_LN2_BOIL, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF)
    J["cryo_end_cop_MW"] = J["final_p_cryo_MW"] * J["cop_sim"] / COP_BOOK
    return J


HDR = ["case", "turns", "MW", "days (months)", "I cap to day", "P cap days (first to last)", "100 kV from day (days at it)", "thrust max N/LIM (day)",
       "end: thrust N, A, kV, slip %", "hyst peak kW/LIM", "cryo peak MW (sim / COP)", "cryo end MW (sim / COP)", "plate K", "tape t/site", "LIM work EJ", "note"]


def fmt(x, f=".1f"):
    return "never" if x is None else format(x, f)


def row(J):
    return [J["name"], J["N_TURNS"], f"{J['P_MW']:.0f}", f"{J['t_end_days']:.1f} ({J['t_end_days']/30.4375:.1f})", fmt(J.get("current_cap_last_day")),
            (f"{J['days_P']:.0f} ({J['power_cap_first_day']:.0f} to {J['power_cap_last_day']:.0f})" if J["days_P"] else "never"),
            (f"{J['volt_limit_day']:.1f} ({J['days_at_V']:.0f})" if J.get("volt_limit_day") is not None else "never"),
            f"{J['thrust_max_N']:,.0f} ({J['thrust_max_day']:.1f})",
            f"{J['final_thrust_N']:.0f}, {J['final_current_A']:.1f}, {J['final_volts_V']/1e3:.1f}, {J['final_slip_pct']:.2f}",
            f"{J['p_hyst_max_kW']:.2f}", f"{J['p_cryo_max_MW']:.2f} / {J['cryo_peak_cop_MW']:.2f}", f"{J['final_p_cryo_MW']:.2f} / {J['cryo_end_cop_MW']:.2f}",
            f"{J['plate_t_max_K']:.1f}", f"{J['tape_t_site']:.2f}", f"{J['E_lim_loop_EJ']:.3f}",
            ("copied from rising_radius_design (corrected rule)" if J["copied"] else "") + (" HYSTERESIS_MODEL = book" if "bookhyst" in J["name"] else "")]


def curve(rows, P):
    pts = sorted([(J["N_TURNS"], J) for J in rows if abs(J["P_MW"] - P) < 1e-9 and "bookhyst" not in J["name"]], key=lambda t: t[0])
    if not pts:
        return None
    floor = min(J["t_end_days"] for _, J in pts)
    best = next((N, J) for N, J in pts if J["t_end_days"] == floor)
    within = [N for N, J in pts if J["t_end_days"] <= floor * (1 + KNEE_TOL)]
    vlim = [N for N, J in pts if J.get("volt_limit_day") is not None]
    never = [N for N, J in pts if J.get("volt_limit_day") is None]
    return dict(P=P, pts=pts, floor=floor, best=best, within=within, vlim=vlim, never=never, largest_never=max(never) if never else None)


def figure(curves, path):
    """A candidate Figure 9.3 for the book: days against turns at 8 and 12 MW, one idea, few labels, large type for a 6 by 9 inch page.
    Open markers: the coil never reaches 100 kV; filled: it does (the voltage-limited part of the curve)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, ORANGE, INK = "#0072B2", "#D55E00", "#222222"
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for C, col, lab in ((curves.get(12), BLUE, "12 MW per site"), (curves.get(8), ORANGE, "8 MW per site")):
        if not C:
            continue
        Ns = [N for N, _ in C["pts"]]; ds = [J["t_end_days"] for _, J in C["pts"]]
        ax.plot(Ns, ds, "-", color=col, lw=2.2, label=lab, zorder=2)
        for N, J in C["pts"]:
            filled = J.get("volt_limit_day") is not None
            ax.plot([N], [J["t_end_days"]], "o", ms=7, color=col, mfc=col if filled else "white", mew=1.8, zorder=3)
    ax.set_xlabel("Turns per phase coil", fontsize=13, color=INK); ax.set_ylabel("Deployment, days", fontsize=13, color=INK)
    ax.tick_params(labelsize=12)
    ax.grid(True, color="#DDDDDD", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    leg = ax.legend(fontsize=12, frameon=False, loc="upper right")
    ax.text(0.02, 0.97, "filled: the coil reaches 100 kV\nopen: it never does", transform=ax.transAxes, fontsize=11, color=INK, ha="left", va="top")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{path}.{ext}", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def table():
    rows = [load(n) for n in CASES]
    rows = [r for r in rows if r]
    out = [row(J) for J in rows]
    md = "| " + " | ".join(HDR) + " |\n|" + "---|" * len(HDR) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out) + "\n"
    lines = []
    curves = {P: curve(rows, P) for P in (12, 8, 16)}
    for P, C in curves.items():
        if not C:
            continue
        N, J = C["best"]
        lines.append(f"{P} MW: days against turns: " + "; ".join(f"{n} turns {j['t_end_days']:.1f} d" for n, j in C["pts"]))
        lines.append(f"  the minimum {J['t_end_days']:.1f} days at {N} turns; within {KNEE_TOL*100:.0f} % of it: {', '.join(str(n) for n in C['within'])} turns;"
                     f" the coil reaches 100 kV at: {', '.join(str(n) for n in C['vlim']) or 'none'}; never at: {', '.join(str(n) for n in C['never']) or 'none'};"
                     f" the largest N that never reaches 100 kV: {C['largest_never']}")
    # the voltage check: V_end against N sqrt(P) for the runs that end power-capped (never at 100 kV)
    lines.append("The voltage at the end against N sqrt(P) for the runs that never reach 100 kV (V_end / (N sqrt(P / MW)) in volts):")
    for J in sorted(rows, key=lambda j: (j["P_MW"], j["N_TURNS"])):
        if J.get("volt_limit_day") is None and "bookhyst" not in J["name"]:
            lines.append(f"  {J['N_TURNS']} turns, {J['P_MW']:.0f} MW: V_end {J['final_volts_V']/1e3:.1f} kV, V / (N sqrt(P)) = {J['final_volts_V'] / (J['N_TURNS'] * math.sqrt(J['P_MW'])):.1f} V"
                         f" (I_end {J['final_current_A']:.1f} A, f_end {J['final_f_supply_Hz']:.1f} Hz, slip {J['final_slip_pct']:.2f} %)")
    # the hysteresis of the design coil on the book's loss factor
    Jb = next((J for J in rows if J["name"] == "N150_P12MW_bookhyst"), None); Jd = next((J for J in rows if J["name"] == "N150_P12MW"), None)
    if Jb and Jd:
        lines.append(f"The design coil's hysteresis peak: {Jd['p_hyst_max_kW']:.2f} kW per LIM on the Brandt-Indenbom model (the design run), {Jb['p_hyst_max_kW']:.2f} kW per LIM"
                     f" on the 2025 book's loss factor (HYSTERESIS_MODEL = book, the same run otherwise: {Jb['t_end_days']:.1f} days, cryo peak {Jb['p_cryo_max_MW']:.2f} MW).")
    txt = "\n".join(lines)
    print(md); print(txt)
    figure(curves, os.path.join(out_dir(), "fig09_03_turns_candidate"))
    with open(os.path.join(out_dir(), "rising_radius_turns_rerun.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n" + txt + "\n")
    with open(os.path.join(out_dir(), "rising_radius_turns_rerun.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(HDR); wr.writerows(out)
    with open(os.path.join(out_dir(), "rising_radius_turns_rerun.json"), "w", encoding="utf-8") as fh:
        json.dump({J["name"]: {k: v for k, v in J.items() if not isinstance(v, (list, dict))} for J in rows}, fh, indent=1, default=str)
    print(f"written rising_radius_turns_rerun.md, .csv, .json and fig09_03_turns_candidate.png|svg in {out_dir()}")


def main():
    names = list(CASES)
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    if only:
        names = [n for n in only.split(",") if n in CASES]
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "6"))
    full = "--quick" not in sys.argv
    if "--run" in sys.argv:
        if cfg.SLIP_RULE != "both":
            raise SystemExit("ORBITAL_RING_SLIP_RULE must be 'both' for this rerun")
        print(f"running {len(names)} cases, {jobs} at a time" + (" (full steps)" if full else " (quick steps)") + f", into {out_dir()}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: run_case(n, full), names):
                print(f"  {name:22} rc={rc} {dt / 60:.1f} min" + ("  (copied)" if dt == 0.0 else ""), flush=True)
    if "--table" in sys.argv:
        table()
    if "--run" not in sys.argv and "--table" not in sys.argv:
        print(__doc__)


if __name__ == "__main__":
    main()
