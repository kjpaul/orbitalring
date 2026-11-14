#!/usr/bin/env python3
"""
The design moves to 165 turns (2026-11-14, Vol_III/Claude_Code_instructions_design_165_turns_2026-11-14.md): the slip-rule check and the
Chapter 9 set at 165 turns, two LIMs per site, 12 MW, the rising quadratic, the rising radius, 100 kV, full time steps.
Part 1, the slip rule (--part=slip): the turns of the brief rerun under SLIP_RULE "optimal" (lim_config: the slip at a binding limit goes to
  the plate's thrust-optimal slip speed instead of the fixed 0.5 percent): 150 to 180 turns at 12 MW, 170, 190, 208 at 8 MW, into
  runs/rising_radius_design_165/slip_optimal/<case>/; the table sets them beside Section 17's runs under "both".
Part 2, the Chapter 9 set (--part=ch9): under the rule Part 1 settles on (--rule=both|optimal): two LIMs with 165 turns at 8, 12, 16, 20, 24 MW;
  two LIMs with 150 and 208 turns at 12 and 16 MW; four LIMs with 165 turns at 8, 16, 24 MW; six LIMs with 165 turns at 24 MW; into
  runs/rising_radius_design_165/ch9_<rule>/<case>/.  Cases that exist under the same rule elsewhere (Section 15.2's and 17's) are copied in.
  --table writes the tables of the brief (9.2's values, 9.3's two columns, 9.4's rows, the power curve) and redraws Figure 9.3.
Run: python tools/rising_radius_design_165.py --part=slip --run [--jobs=6]      then --part=slip --table
     python tools/rising_radius_design_165.py --part=ch9 --rule=both --run       then --part=ch9 --rule=both --table
Environment: ORBITAL_RING_OUT (default runs).  The runs set ORBITAL_RING_SLIP_RULE themselves.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import concurrent.futures
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import ring_altitude as ring
import lim_config as cfg
import lim_hysteresis as lh
import lim_physics as phys
import rising_radius_study as study

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
BASE = os.path.join(ring.out_dir(), "rising_radius_design_165")
TURNS_REF = os.path.join(ring.out_dir(), "rising_radius_turns_rerun")          # Section 17, SLIP_RULE both
DESIGN_REF = os.path.join(ring.out_dir(), "rising_radius_design")               # Section 15.2, SLIP_RULE both
FILES = ["rising_radius_lim_simulation.txt", "rising_radius_daily.csv", "rising_radius_deployment_tension.csv", "rising_radius_summary.json"]
COP_BOOK = 0.0106
LAYERS, W_TAPE = 2, 3e-3
PANEL_M_PER_MW = 55.0 / 8.0                                                      # the deployment panel: 55 m at 8 MW, scaled with the site power (Section 13.3)

SLIP_CASES = {f"N{N}_P12MW": (2, N, 12) for N in (150, 155, 160, 165, 170, 175, 180)}
SLIP_CASES.update({f"N{N}_P8MW": (2, N, 8) for N in (170, 190, 208)})
# the rest of Section 17's grid, so that Figure 9.3 is drawn from one rule
SLIP_CASES.update({f"N{N}_P12MW": (2, N, 12) for N in (120, 135, 140, 145, 190, 208, 240, 280)})
SLIP_CASES.update({f"N{N}_P8MW": (2, N, 8) for N in (150, 230, 250)})
CH9_CASES = {}
for mw in (8, 12, 16, 20, 24):
    CH9_CASES[f"L2_N165_P{mw}MW"] = (2, 165, mw)
for N in (150, 208):
    for mw in (12, 16):
        CH9_CASES[f"L2_N{N}_P{mw}MW"] = (2, N, mw)
for mw in (8, 16, 24):
    CH9_CASES[f"L4_N165_P{mw}MW"] = (4, 165, mw)
CH9_CASES["L6_N165_P24MW"] = (6, 165, 24)
# where a case already exists under SLIP_RULE both
REF_BOTH = {"L2_N150_P12MW": (DESIGN_REF, "L2_N150_P12MW"), "L2_N150_P16MW": (DESIGN_REF, "L2_N150_P16MW"), "L2_N208_P12MW": (DESIGN_REF, "L2_N208_P12MW"),
            "L2_N208_P16MW": (DESIGN_REF, "L2_N208_P16MW"), "L2_N165_P12MW": (TURNS_REF, "N165_P12MW"), "L2_N165_P16MW": (TURNS_REF, "N165_P16MW")}


def case_dir(part, rule):
    return os.path.join(BASE, "slip_optimal" if part == "slip" else f"ch9_{rule}")


def run_case(part, rule, name, spec, full=True):
    lims, N, mw = spec
    out = os.path.join(case_dir(part, rule), name)
    os.makedirs(out, exist_ok=True)
    if part == "ch9" and rule == "both" and name in REF_BOTH and os.path.exists(os.path.join(*REF_BOTH[name], "rising_radius_summary.json")):
        for f in FILES:
            shutil.copy2(os.path.join(*REF_BOTH[name], f), os.path.join(out, f))
        return name, 0, 0.0
    over = {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": N}
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_OUT=out, ORBITAL_RING_SLIP_RULE=rule)
    if lims > 2:
        over["N_PLATES_PER_SIDE"] = lims // 2
        env["ORBITAL_RING_HW_KGM"] = f"{study.hw_kgm(lims // 2):.1f}"
    args = "['lim_simulation.py', '--rising-radius'" + ("" if full else ", '--quick'") + "]"
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv={args}; runpy.run_path('lim_simulation.py', run_name='__main__')")
    t0 = time.time()
    with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
        r = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    return name, r.returncode, time.time() - t0


def parse_log(path):
    """The peaks and their times from the log's max/min block, the cable temperature per month, the limits line."""
    txt = open(path, encoding="utf-8", errors="replace").read()
    out = {}
    for key in ("p_eddy_max", "p_hyst_max", "p_thrust_max", "thrust_max", "cable_temp_max", "plate_t_max", "current_max", "volts_max", "site_power_max", "lim_power_max"):
        m = re.search(rf"^\s*{key}\s+([\d.]+)\s+(\S+)\s+([\d.]+)\s+(days|months)", txt, re.M)
        if m:
            v, unit, t, tu = float(m.group(1)), m.group(2), float(m.group(3)), m.group(4)
            out[key] = (v, unit, t * (30.4375 if tu == "months" else 1.0))
    m = re.search(r"^Limits: (.*)$", txt, re.M)
    out["limits"] = m.group(1) if m else ""
    # the monthly table: Month | Progress | Voltage | Current | Thrust | Site Power | Cryo Rad W | T_cable | ...
    months = []
    for line in txt.splitlines():
        if re.match(r"^\s*\d+(\.\d+)?\s*\|", line) and line.count("|") >= 10:
            cells = [c.strip() for c in line.split("|")]
            try:
                months.append((float(cells[0]), float(cells[7].split()[0])))
            except (ValueError, IndexError):
                pass
    out["T_cable_months"] = months
    return out


def load(part, rule, name):
    p = os.path.join(case_dir(part, rule), name, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p)); J["name"] = name
    N = J["N_TURNS"]; J["P_MW"] = J["MAX_SITE_POWER"] / 1e6
    J["tape_t_site"] = lh.tape_mass_per_lim(N, LAYERS, W_TAPE) * J["LIMS_PER_SITE"] / 1e3
    J["tape_kg_lim"] = lh.tape_mass_per_lim(N, LAYERS, W_TAPE)
    J["tape_km_lim"] = J["L_HTS_LIM_m"] / 1e3
    J["coil_width_m"] = lh.winding_width("flat", N, LAYERS, W_TAPE)
    J["cryo_peak_cop_MW"] = J["q_cryo_cold_max_kW"] / COP_BOOK / 1e3
    J["cop_sim"] = phys.calc_cryo_cop(cfg.T_LN2_BOIL, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF)
    J["cryo_end_cop_MW"] = J["final_p_cryo_MW"] * J["cop_sim"] / COP_BOOK
    vl = J.get("volt_limit_day"); J["days_at_V"] = (J["t_end_days"] - vl) if vl is not None else 0.0
    J["panel_m"] = PANEL_M_PER_MW * J["site_power_max_MW"]
    J["ring_GW"] = J["site_power_max_MW"] * J["lim_sites"] / 1e3
    d = ring.design(m_hw=J["M_HARDWARE_kg_m"]); J["cable_kg_m"] = d["m_cable"]; J["cable_Mt"] = d["m_cable"] * d["L_ring"] / 1e9
    J.update(parse_log(os.path.join(case_dir(part, rule), name, "rising_radius_lim_simulation.txt")))
    return J


def ref_both(name):
    """Section 17's or 15.2's run of the same case under 'both', for the slip table."""
    p = os.path.join(TURNS_REF, name, "rising_radius_summary.json")
    return json.load(open(p)) if os.path.exists(p) else None


def fmt(x, f=".1f"):
    return "never" if x is None else format(x, f)


def table_slip():
    rows = []
    for name, spec in SLIP_CASES.items():
        J = load("slip", "optimal", name); R = ref_both(name)
        if J:
            rows.append((name, spec, J, R))
    hdr = "| Turns | MW | Days, both (Section 17) | Days, optimal | Change | 100 kV from day, both / optimal | Days at 100 kV, both / optimal | P cap days, both / optimal | End: N/LIM, A, slip %, both | End, optimal | Hyst peak kW/LIM, both / optimal | Cryo peak MW (sim), both / optimal |"
    lines = [hdr, "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, (lims, N, mw), J, R in rows:
        lines.append(f"| {N} | {mw} | {R['t_end_days']:.1f} | {J['t_end_days']:.1f} | {J['t_end_days']-R['t_end_days']:+.1f} d ({(J['t_end_days']/R['t_end_days']-1)*100:+.2f} %) | "
                     f"{fmt(R.get('volt_limit_day'))} / {fmt(J.get('volt_limit_day'))} | {R.get('days_at_volt_limit', 0):.0f} / {J.get('days_at_volt_limit', 0):.0f} | "
                     f"{R.get('days_at_power_cap', 0):.0f} / {J.get('days_at_power_cap', 0):.0f} | {R['final_thrust_N']:.0f}, {R['final_current_A']:.1f}, {R['final_slip_pct']:.2f} | "
                     f"{J['final_thrust_N']:.0f}, {J['final_current_A']:.1f}, {J['final_slip_pct']:.2f} | {R['p_hyst_max_kW']:.2f} / {J['p_hyst_max_kW']:.2f} | {R['p_cryo_max_MW']:.2f} / {J['p_cryo_max_MW']:.2f} |")
    pts12 = sorted([(N, J["t_end_days"]) for name, (l, N, mw), J, R in rows if mw == 12]); pts8 = sorted([(N, J["t_end_days"]) for name, (l, N, mw), J, R in rows if mw == 8])
    for lab, pts in (("12 MW", pts12), ("8 MW", pts8)):
        if pts:
            floor = min(d for _, d in pts); best = next(N for N, d in pts if d == floor)
            within = [N for N, d in pts if d <= 1.02 * floor]
            lines.append(f"\n{lab} under optimal: " + "; ".join(f"{N} turns {d:.1f} d" for N, d in pts) + f". The minimum {floor:.1f} days at {best} turns; within 2 percent: {', '.join(map(str, within))}.")
    r165 = next((J for name, (l, N, mw), J, R in rows if N == 165 and mw == 12), None); b165 = ref_both("N165_P12MW")
    if r165 and b165:
        ch = (r165["t_end_days"] / b165["t_end_days"] - 1) * 100
        lines.append(f"\nThe 165-turn design run moves by {ch:+.2f} percent under the optimal rule ({r165['t_end_days']:.1f} against {b165['t_end_days']:.1f} days): "
                     + ("under half a percent, so 'both' stays the default and Section 17's runs stand." if abs(ch) < 0.5 else "more than half a percent, so the optimal rule becomes the default and the Chapter 9 set is run under it."))
    md = "\n".join(lines)
    print(md)
    os.makedirs(case_dir("slip", "optimal"), exist_ok=True)
    with open(os.path.join(case_dir("slip", "optimal"), "slip_optimal_table.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n")
    return md


def table_ch9(rule):
    rows = {}
    for name in CH9_CASES:
        J = load("ch9", rule, name)
        if J:
            rows[name] = J
    L = []
    d12, d8 = rows.get("L2_N165_P12MW"), rows.get("L2_N165_P8MW")
    if d12:
        N = d12["N_TURNS"]; cl = lh.coil_limit("flat", N, LAYERS, W_TAPE)
        L.append("### Table 9.2's values for the 165-turn coil\n")
        L.append("| Parameter | Value |\n|---|---|")
        L.append(f"| Turns in each coil | {N} |")
        L.append(f"| Peak coil current | {d12['I_TARGET']:.1f} A (the cap as configure() leaves it; the winding's own limit {cl[0]:.1f} A) |")
        L.append(f"| Current limit of one tape in the field of the winding | {cl[0]:.1f} A at the worst tape, at the edge of the sheet, in {cl[1]['B']:.3f} T at {cl[1]['theta']:.1f} degrees from its face's normal |")
        L.append(f"| Layers of tape in each turn | 2, the turns laid flat, side by side, in a sheet {d12['coil_width_m']*1e3:.0f} mm wide |")
        L.append(f"| Tape in each LIM | {d12['tape_kg_lim']:,.0f} kg ({2*d12['tape_kg_lim']:,.0f} kg per site), {d12['tape_km_lim']:,.0f} km of 3 mm tape ({2*d12['tape_km_lim']:,.0f} km per site) |")
        L.append(f"| Current in each tape | {d12['I_TARGET']/2:.1f} A, half of its limit |")
        L.append(f"| Supply frequency at the end of deployment | {d12['final_f_supply_Hz']:.1f} Hz |")
        L.append(f"| Peak thrust of each LIM | about {d12['thrust_max_N']:,.0f} N |")
        L.append(f"| Deployment time | {d12['t_end_days']:.0f} days at 12 MW" + (f", {d8['t_end_days']:.0f} days at 8 MW |" if d8 else " |"))
        L.append("")
    if d12 and d8:
        L.append("### Table 9.3, both columns, 165 turns\n")
        L.append("| Quantity | 12 MW per site, the design | 8 MW per site |\n|---|---|---|")
        def gov(J):
            parts = [f"coil current for the first {J['current_cap_last_day']:.0f} days"]
            if J.get("power_cap_first_day") is not None:
                parts.append(f"then site power from day {J['power_cap_first_day']:.0f} to day {J['power_cap_last_day']:.0f}")
            if J.get("volt_limit_day") is not None:
                parts.append(f"then the coil voltage at 100 kV from day {J['volt_limit_day']:.0f} to the end ({J['days_at_V']:.0f} days)")
            else:
                parts.append("the coil voltage never reaches 100 kV")
            return "; ".join(parts)
        def T_settle(J):
            ms = J.get("T_cable_months", [])
            if not ms:
                return "n/a"
            T1 = J["cable_temp_max"][0] if "cable_temp_max" in J else float("nan")
            late = [T for m, T in ms if m >= 6]
            return f"{T1:.0f} K on the first day, settling near {late[-1]:.0f} K" if late else f"{T1:.0f} K on the first day"
        R = [("Deployment time", lambda J: f"{J['t_end_days']:.0f} days ({J['t_end_days']/30.4375:.1f} months)"),
             ("Governing limit, in order", gov),
             ("Peak coil current", lambda J: f"{J['current_max_A']:.0f} A, from the first day to day {J['current_cap_last_day']:.0f}"),
             ("Coil voltage at the end", lambda J: f"{J['final_volts_V']/1e3:.0f} kV, with {J['final_current_A']:.0f} A in the coil"),
             ("Site power", lambda J: f"{J['P_MW']:.0f} MW from day {J['power_cap_first_day']:.0f} to day {J['power_cap_last_day']:.0f} ({J['days_at_power_cap']:.0f} days)" if J.get("power_cap_first_day") is not None else "the cap never binds"),
             ("Peak thrust per LIM", lambda J: f"{J['thrust_max_N']:,.0f} N, at {J['thrust_max_day']/30.4375:.1f} months (day {J['thrust_max_day']:.0f}); {J['final_thrust_N']:.0f} N at the end"),
             ("Peak thrust power per LIM", lambda J: f"{J['p_thrust_max'][0]:.1f} MW, at {J['p_thrust_max'][2]/30.4375:.1f} months" if "p_thrust_max" in J else "n/a"),
             ("Mechanical energy delivered", lambda J: f"{J['E_lim_loop_EJ']:.1f} EJ ({J['E_lim_loop_EJ']:.3f})"),
             ("Peak eddy current heating per LIM", lambda J: f"{J['p_eddy_max'][0]:.0f} kW, on day {J['p_eddy_max'][2]:.0f}" if "p_eddy_max" in J else "n/a"),
             ("Peak hysteresis heating per LIM", lambda J: f"{J['p_hyst_max_kW']:.1f} kW, at {J['p_hyst_max'][2]/30.4375:.1f} months" if "p_hyst_max" in J else f"{J['p_hyst_max_kW']:.1f} kW"),
             ("Peak cable temperature", T_settle),
             ("Peak cryogenic plant per site (by the book's COP; the simulation's own)", lambda J: f"{J['cryo_peak_cop_MW']:.2f} MW ({J['p_cryo_max_MW']:.2f}), at the end {J['cryo_end_cop_MW']:.2f} ({J['final_p_cryo_MW']:.2f})")]
        for lab, f in R:
            L.append(f"| {lab} | {f(d12)} | {f(d8)} |")
        L.append("")
    L.append("### Table 9.4's rows\n")
    L.append("| Site power | Turns | LIMs per site | Deployment | 100 kV from day | Panel width | HTS tape per site | Cable mass | Peak cryogenic plant (book COP) | Plate on day one | Peak hysteresis per LIM |\n|---|---|---|---|---|---|---|---|---|---|---|")
    order = ["L2_N165_P8MW", "L2_N165_P12MW", "L2_N165_P16MW", "L2_N165_P20MW", "L2_N165_P24MW", "L2_N150_P12MW", "L2_N150_P16MW", "L2_N208_P12MW", "L2_N208_P16MW",
             "L4_N165_P8MW", "L4_N165_P16MW", "L4_N165_P24MW", "L6_N165_P24MW"]
    for name in order:
        J = rows.get(name)
        if not J:
            L.append(f"| {name} | (not run) |"); continue
        L.append(f"| {J['P_MW']:.0f} MW | {J['N_TURNS']} | {J['LIMS_PER_SITE']} | {J['t_end_days']:.0f} days ({J['t_end_days']/30.4375:.1f} months) | {fmt(J.get('volt_limit_day'), '.0f')} | {J['panel_m']:.0f} m | "
                 f"{J['tape_t_site']:.2f} t | {J['cable_kg_m']:,.0f} kg/m, {J['cable_Mt']:,.0f} Mt | {J['cryo_peak_cop_MW']:.2f} MW | {J['plate_t_max_K']:.0f} K | {J['p_hyst_max_kW']:.1f} kW |")
    L.append("")
    pts = sorted([(J["P_MW"], J["t_end_days"]) for n, J in rows.items() if n.startswith("L2_N165")])
    if len(pts) > 1:
        floor = min(d for _, d in pts)
        L.append("### The power curve of the design (two LIMs, 165 turns)\n")
        L.append(", ".join(f"{d:.1f} days at {p:.0f} MW" for p, d in pts) + ". Days per megawatt between the steps: " +
                 "; ".join(f"{(d0-d1)/(p1-p0):.1f} from {p0:.0f} to {p1:.0f} MW" for (p0, d0), (p1, d1) in zip(pts, pts[1:])) +
                 f". The floor is {floor:.1f} days; 12 MW sits {(next(d for p, d in pts if p == 12)/floor-1)*100:.1f} percent above it; the knee by the 2 percent rule is "
                 f"{min(p for p, d in pts if d <= 1.02*floor):.0f} MW.")
    md = "\n".join(L)
    print(md)
    with open(os.path.join(case_dir("ch9", rule), "ch9_tables.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n")
    with open(os.path.join(case_dir("ch9", rule), "ch9_rows.json"), "w", encoding="utf-8") as fh:
        json.dump({n: {k: v for k, v in J.items() if not isinstance(v, (list, dict, tuple))} for n, J in rows.items()}, fh, indent=1, default=str)
    return md


def figure93(rule):
    """Figure 9.3 redrawn: Section 17's curves under 'both' (or the optimal runs where they exist under 'optimal'), the legend clear of the 280-turn point."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pts = {12: [], 8: []}
    for P in (12, 8):
        for d in sorted(os.listdir(TURNS_REF)):
            m = re.match(rf"N(\d+)_P{P}MW$", d)
            if not m:
                continue
            J = json.load(open(os.path.join(TURNS_REF, d, "rising_radius_summary.json")))
            if rule == "optimal":
                Jo = load("slip", "optimal", d)
                if Jo:
                    J = Jo
            pts[P].append((int(m.group(1)), J["t_end_days"], J.get("volt_limit_day") is not None))
    BLUE, ORANGE, INK = "#0072B2", "#D55E00", "#222222"
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for P, col, lab in ((12, BLUE, "12 MW per site"), (8, ORANGE, "8 MW per site")):
        Ns = [n for n, _, _ in pts[P]]; ds = [d for _, d, _ in pts[P]]
        ax.plot(Ns, ds, "-", color=col, lw=2.2, label=lab, zorder=2)
        for n, d, vl in pts[P]:
            ax.plot([n], [d], "o", ms=7, color=col, mfc=col if vl else "white", mew=1.8, zorder=3)
    ax.set_xlabel("Turns per phase coil", fontsize=13, color=INK); ax.set_ylabel("Deployment, days", fontsize=13, color=INK)
    ax.tick_params(labelsize=12); ax.grid(True, color="#DDDDDD", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(fontsize=12, frameon=False, loc="upper left", bbox_to_anchor=(0.47, 0.99))
    ax.text(0.02, 0.97, "filled: the coil reaches 100 kV\nopen: it never does", transform=ax.transAxes, fontsize=11, color=INK, ha="left", va="top")
    fig.tight_layout()
    out = os.path.join(BASE, f"fig09_03_turns_{rule}")
    for ext in ("png", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"figure written: {out}.png|svg")


def main():
    part = next((a[7:] for a in sys.argv if a.startswith("--part=")), "slip")
    rule = next((a[7:] for a in sys.argv if a.startswith("--rule=")), "optimal" if part == "slip" else "both")
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "6"))
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    cases = SLIP_CASES if part == "slip" else CH9_CASES
    names = [n for n in cases if not only or n in only.split(",")]
    full = "--quick" not in sys.argv
    if "--run" in sys.argv:
        print(f"running {len(names)} cases under SLIP_RULE {rule}, {jobs} at a time, into {case_dir(part, rule)}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: run_case(part, rule, n, cases[n], full), names):
                print(f"  {name:16} rc={rc} {dt/60:.1f} min" + ("  (copied)" if dt == 0.0 else ""), flush=True)
    if "--table" in sys.argv:
        if part == "slip":
            table_slip()
        else:
            table_ch9(rule)
    if "--figure" in sys.argv:
        figure93(rule)
    if not any(a in sys.argv for a in ("--run", "--table", "--figure")):
        print(__doc__)


if __name__ == "__main__":
    main()
