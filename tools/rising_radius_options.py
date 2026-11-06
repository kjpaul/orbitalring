#!/usr/bin/env python3
"""
The deployment options table (round 3 of 2026-11-06): the eight cases of the brief, each a full-time-step run of
lim_simulation.py --rising-radius (1 s steps on the first day, then 10 s), for the table of Chapter 9's section on
shortening the deployment.  The design stays at 8 MW per site, 208 turns and two LIMs per site; the table shows what
the other cases buy and cost.  The cases are those of tools/rising_radius_study.py (same names, same overrides, same
ORBITAL_RING_HW_KGM for the cases with more reaction plates), run into runs/rising_radius_options/<case>/.  The 16 MW
case is the full-step run of runs/rising_radius_16MW/ (Section 10.1 of the report) copied in, when that folder exists.
Book: Orbital Ring Engineering, Volume I: Chapter 9 (the new section on shortening the deployment), Chapter 12 (the hang
schedule of the options), Chapter 3 (the quadratic's design point with the extra hardware).
Run: python tools/rising_radius_options.py --run [--jobs=7] [--only=a,b]     (the full-step runs, hours each)
     python tools/rising_radius_options.py --table                           (the table from the summaries and the daily tables,
                                                                               to the console, rising_radius_options.md and .csv)
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_os.environ.setdefault("ORBITAL_RING_STUDY_DIR", "rising_radius_options")
import concurrent.futures
import csv
import json
import os
import shutil
import sys
import time

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import ring_altitude as ring
import lim_hysteresis as lh
import rising_radius_study as study

# the brief's order: 1 the design, 2 12 MW, 3 16 MW, 4 four LIMs at 8 MW, 5 four LIMs at 16 MW, 6 the same with 150 turns,
# 7 six LIMs at 24 MW, 8 the same with 150 turns
CASES = ["base_8MW", "P12MW", "P16MW", "L4", "P16MW_L4", "P16MW_N150_L4", "P24MW_L6", "P24MW_N150_L6"]
COPY_FROM = {"P16MW": "rising_radius_16MW"}     # full-step runs that exist already, copied in rather than rerun
FILES = ["rising_radius_lim_simulation.txt", "rising_radius_daily.csv", "rising_radius_deployment_tension.csv", "rising_radius_summary.json"]
FLOOR_KM = 150.0                                 # the hang: the lines let down to 150 km, no lower (Section 7 of the report)
HANGS_KM = [450.0, 550.0]                        # the hang depths Chapter 12 wants the schedule for
PANEL_AREA_8MW, PANEL_WIDTH_8MW = 27_500.0, 55.0 # the deployment panel at 8 MW per site, Chapter 8: 291 W/m^2 averaged around the ring
BASE_CABLE_MT = None                             # filled from the baseline run


def options_dir():
    return study.study_dir()


def run_case(name):
    src = COPY_FROM.get(name)
    if src and os.path.exists(os.path.join(ring.out_dir(), src, "rising_radius_summary.json")):
        out = os.path.join(options_dir(), name)
        os.makedirs(out, exist_ok=True)
        for f in FILES:
            shutil.copy2(os.path.join(ring.out_dir(), src, f), os.path.join(out, f))
        return name, 0, 0.0
    return study.run_case(name, full=True)


def hang_days(name):
    """From the daily table: the first day on which a hang of each depth, hung from the casing, is above the 150 km floor."""
    p = os.path.join(options_dir(), name, "rising_radius_daily.csv")
    if not os.path.exists(p):
        return {}
    rows = list(csv.DictReader(open(p)))
    out = {}
    for h in HANGS_KM:
        day = next((float(r["day"]) for r in rows if float(r["altitude_km"]) - h >= FLOOR_KM - 1e-6), None)
        out[h] = day
    return out


def load(name):
    p = os.path.join(options_dir(), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p))
    J["name"] = name
    J["tape_t_site"] = lh.tape_mass_per_lim(J["N_TURNS"], 2, 3e-3) * J["LIMS_PER_SITE"] / 1e3
    d = ring.design(m_hw=J["M_HARDWARE_kg_m"])                   # the quadratic with this case's hardware
    J["design"] = d
    J["cable_Mt"] = J["M_CABLE_M_kg_m"] * ring.L_RING / 1e9
    J["hang"] = hang_days(name)
    return J


def table():
    rows = [load(n) for n in CASES]
    rows = [r for r in rows if r]
    base = rows[0] if rows and rows[0]["name"] == "base_8MW" else None
    hdr = ["case", "MW/site", "turns", "LIMs/site", "days", "months", "phase 1 d", "phase 2 d", "100 kV from day", "LIM work EJ",
           "peak site MW", "ring GW", "panel m2/site", "panel width m", "plate K", "cryo peak MW", "cryo end MW",
           "LIM plates kg/m", "lev plates kg/m", "hardware vs base kg/m", "structural kg/m", "cable kg/m", "cable Mt", "cable vs base Mt",
           "v_cable m/s", "tension GN", "tape t/site", "stator t/site", "cryo plant MW/site", "anchor MN", "anchor kN/station",
           "450 km hang clears day", "550 km hang clears day"]
    out = []
    for J in rows:
        d = J["design"]; m_lim = J["M_LIM_PLATE_kg_m"] * 2 * J["N_PLATES_PER_SIDE"]; m_lev = J["M_HARDWARE_kg_m"] - m_lim
        ps = J["site_power_max_MW"]
        fmt = lambda x, f: "never" if x is None else format(x, f)
        out.append([J["name"], f"{J['MAX_SITE_POWER'] / 1e6:.0f}", J["N_TURNS"], J["LIMS_PER_SITE"], f"{J['t_end_days']:.1f}", f"{J['t_end_days'] / 30.4375:.1f}",
                    f"{(J.get('t_phase2_days') or J['t_end_days']):.1f}", f"{(J.get('phase2_duration_days') or 0):.1f}", fmt(J.get("volt_limit_day"), ".1f"),
                    f"{J['E_lim_loop_EJ']:.3f}", f"{ps:.2f}", f"{ps * J['lim_sites'] / 1e3:,.0f}", f"{PANEL_AREA_8MW * ps / 8:,.0f}", f"{PANEL_WIDTH_8MW * ps / 8:.0f}",
                    f"{J['plate_t_max_K']:.1f}", f"{J['p_cryo_max_MW']:.2f}", f"{J['final_p_cryo_MW']:.2f}",
                    f"{m_lim:,.0f}", f"{m_lev:,.0f}", f"{J['M_HARDWARE_kg_m'] - (base['M_HARDWARE_kg_m'] if base else 0):+,.0f}",
                    f"{J['M_CABLE_STRUCTURAL_kg_m']:,.0f}", f"{J['M_CABLE_M_kg_m']:,.0f}", f"{J['cable_Mt']:,.0f}", f"{J['cable_Mt'] - (base['cable_Mt'] if base else 0):+,.0f}",
                    f"{d['v_cable']:,.1f}", f"{d['tension'] / 1e9:.1f}", f"{J['tape_t_site']:.1f}", "not in model", f"{J['p_cryo_max_MW']:.2f}",
                    f"{J['anchor_mean_MN']:.1f}", f"{J['anchor_mean_MN'] * 1e3 / J['n_stations']:.0f}",
                    fmt(J["hang"].get(450.0), ".0f"), fmt(J["hang"].get(550.0), ".0f")])
    md = "| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n" + "\n".join("| " + " | ".join(str(c) for c in r) + " |" for r in out) + "\n"
    print(md)
    print("The quadratic's design point with each case's hardware (ring_altitude.design) against the run's end state:")
    for J in rows:
        d = J["design"]
        print(f"  {J['name']:16} hardware {J['M_HARDWARE_kg_m']:,.1f} kg/m: structural {d['m_struct']:,.1f}, with hardware {d['m_cable']:,.1f} kg/m,"
              f" v_cable {d['v_cable']:,.2f} m/s, v_rel {d['v_rel']:,.2f} m/s, tension {d['tension'] / 1e9:.2f} GN, area {d['area']:.2f} m^2, cable {d['m_cable'] * d['L_ring'] / 1e9:,.1f} Mt;"
              f" the run ends at {J['final_v_cable_m_s']:,.2f} m/s, {J['final_tension_GN']:.2f} GN, {J['M_CABLE_STRUCTURAL_kg_m']:,.1f} kg/m structural")
    with open(os.path.join(options_dir(), "rising_radius_options.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    with open(os.path.join(options_dir(), "rising_radius_options.csv"), "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh); wr.writerow(hdr); wr.writerows(out)
    print(f"written rising_radius_options.md and .csv in {options_dir()}")


def main():
    names = CASES
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    if only:
        names = [n for n in only.split(",") if n in CASES]
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "7"))
    if "--run" in sys.argv:
        print(f"running {len(names)} cases at full steps, {jobs} at a time, into {options_dir()}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(run_case, names):
                print(f"  {name:16} rc={rc} {dt / 60:.1f} min" + ("  (copied)" if dt == 0.0 else ""), flush=True)
    if "--table" in sys.argv:
        table()
    if "--run" not in sys.argv and "--table" not in sys.argv:
        print(__doc__)


if __name__ == "__main__":
    main()
