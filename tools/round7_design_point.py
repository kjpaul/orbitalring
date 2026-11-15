#!/usr/bin/env python3
"""
The deployment runs of round 7's design point (2026-11-15, Section 3 of Vol_III/Claude_Code_instructions_bearing_round7_2026-11-15.md): the
turns sweep and the Chapter 9 set rerun with the plates as ruled (four fast plates 1.0 m x 200 mm, two support plates 1.2 m x 100 mm), 4.0 m of
cold width and the cold load of ten connected fast families with the 0.1 mm regulation's heat on the 21 mm2 lead.
Book: Orbital Ring Engineering, Volume III: Chapter 9 (Tables 9.2 to 9.4, Figure 9.3), Chapter 14, Chapter 15.
Run: python tools/round7_design_point.py --list
     python tools/round7_design_point.py --part=slip --run [--jobs=5]     the sweep: 120 to 280 turns at 12 MW, 150 to 250 at 8 MW
     python tools/round7_design_point.py --part=ch9 --run                 the Chapter 9 set (cases the sweep already has are copied)
     python tools/round7_design_point.py --part=ch9 --table --figure
It drives tools/rising_radius_design_165.py's cases, tables and figure with the design point's overrides; the runs go to runs/round7_design/.
Environment: ORBITAL_RING_OUT (default runs).  Everything runs under SLIP_RULE optimal.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import time

import ring_altitude as ring
import rising_radius_design_165 as d

HERE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
BASE = os.path.join(ring.out_dir(), "round7_design")
RHO_PLATE = 7870.0; M_LIM_PLATES = 1872.0
SUP = (1.2, 0.10); FAST = (1.0, 0.20)                                # plate width and thickness, m
COLD_WIDTH = 2 * (2 * 0.4) + 4 * 0.6                                 # the loops' count: two 0.4 m loops per support track, one 0.6 m loop per fast track


def design_point():
    """The hardware, the plates' mean width and thickness for the cable's radiation, the cold width and the extra cold load per site."""
    import bearing_round5 as r5, bearing_round6b as r6, bearing_electrical as be, bearing_round4 as r4, final_bearing as fb
    m_s = RHO_PLATE * SUP[0] * SUP[1]; m_l = RHO_PLATE * FAST[0] * FAST[1]
    hw = M_LIM_PLATES + 2 * m_s + 4 * m_l
    w_mean = (2 * SUP[0] + 4 * FAST[0]) / 6; t_eff = (2 * m_s + 4 * m_l) / (6 * RHO_PLATE * w_mean)
    rs_opt, rl_opt, fin = r6.load_final()
    A = r5.lead_rows(rl_opt)[0][10]["A_m"]                             # the 21 mm2 lead, kept as ruled
    q_lead = be.lead_heat(0.0, A_over_l=A)[0]
    p = os.path.join(ring.out_dir(), "round7_logs", "round7_fast_track.json")
    if os.path.exists(p):
        I_full = json.load(open(p, encoding="utf-8"))["I"]; src = "the accepted fast track's solve"
    else:
        I_full = r5.fast_basics(r6.cached("fast", FAST[0], FAST[1], 0.6), 10)["I"]; src = "the 1.0 m x 200 mm plate's solve with the 200 x 50 mm bundle"
    P_full = I_full**2 * be.RHO_CU_77 * 1.0 / A
    n = 10; coils_site = 1600
    cond = coils_site * n / 20 * 2 * q_lead
    dW_day = 12.7                                                      # N/m per day, the casing's weight falling as the ring rises (round 6a, S4.1)
    p_mean = r4.K_FLUX * 1e-3 * 0.1 / fb.F_LOAD_TRACK + 7 * dW_day / 2 / fb.F_LOAD_TRACK
    joule = coils_site * n / 20 * 2 * p_mean * (20 / n)**2 * P_full
    return dict(hw=hw, m_s=m_s, m_l=m_l, w_mean=w_mean, t_eff=t_eff, cold_width=COLD_WIDTH, extra=cond + joule, cond=cond, joule=joule, I_full=I_full, src=src, A_lead=A)


DP = None


def run_case(part, rule, name, spec, full=True):
    lims, N, mw = spec
    out = os.path.join(d.case_dir(part, rule), name)
    os.makedirs(out, exist_ok=True)
    if part == "ch9":                                                  # the cases the sweep already ran
        twin = f"N{N}_P{mw}MW"
        src = os.path.join(d.case_dir("slip", rule), twin)
        if lims == 2 and os.path.exists(os.path.join(src, "rising_radius_summary.json")):
            for f in d.FILES:
                shutil.copy2(os.path.join(src, f), os.path.join(out, f))
            return name, 0, 0.0
    over = {"MAX_SITE_POWER": mw * 1e6, "N_TURNS": N, "N_LEV_PLATES": 6, "D_LEV": DP["t_eff"], "W_LEV": DP["w_mean"], "LEV_COIL_WIDTH": DP["cold_width"], "Q_COLD_EXTRA_PER_SITE": DP["extra"]}
    hw = DP["hw"]
    if lims > 2:
        over["N_PLATES_PER_SIDE"] = lims // 2
        hw += M_LIM_PLATES * (lims // 2 - 1)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_OUT=out, ORBITAL_RING_SLIP_RULE=rule, ORBITAL_RING_HW_KGM=f"{hw:.1f}")
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius']; runpy.run_path('lim_simulation.py', run_name='__main__')")
    t0 = time.time()
    with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
        r = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, stdout=lf, stderr=subprocess.STDOUT)
    return name, r.returncode, time.time() - t0


def main():
    global DP
    part = next((a[7:] for a in sys.argv if a.startswith("--part=")), "slip")
    jobs = int(next((a[7:] for a in sys.argv if a.startswith("--jobs=")), "5"))
    only = next((a[7:] for a in sys.argv if a.startswith("--only=")), None)
    rule = "optimal"
    d.BASE = BASE
    d.TURNS_REF = os.path.join(BASE, "slip_optimal")                    # the figure and the slip table read this round's own sweep
    d.ref_both = lambda name: (json.load(open(os.path.join(ring.out_dir(), "rising_radius_design_165", "slip_optimal", name, "rising_radius_summary.json")))
                               if os.path.exists(os.path.join(ring.out_dir(), "rising_radius_design_165", "slip_optimal", name, "rising_radius_summary.json")) else None)
    os.makedirs(BASE, exist_ok=True)
    pfile = os.path.join(BASE, "design_point.json")
    if os.path.exists(pfile) and "--list" not in sys.argv:
        DP = json.load(open(pfile, encoding="utf-8"))
    else:
        DP = design_point()
        if "--list" not in sys.argv:
            json.dump(DP, open(pfile, "w", encoding="utf-8"), indent=1)
    if "--list" in sys.argv:
        for k, v in DP.items():
            print(f"  {k}: {v}")
        return
    cases = d.SLIP_CASES if part == "slip" else d.CH9_CASES
    names = [n for n in cases if not only or n in only.split(",")]
    if "--run" in sys.argv:
        todo = [n for n in names if not os.path.exists(os.path.join(d.case_dir(part, rule), n, "rising_radius_summary.json"))]
        print(f"design point: hardware {DP['hw']:,.1f} kg/m, cold width {DP['cold_width']:.1f} m, extra cold load {DP['extra']:,.0f} W per site ({DP['cond']:,.0f} conduction + {DP['joule']:,.0f} Joule)", flush=True)
        print(f"running {len(todo)} of {len(names)} cases, {jobs} at a time, into {d.case_dir(part, rule)}", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            for name, rc, dt in ex.map(lambda n: run_case(part, rule, n, cases[n]), todo):
                print(f"  {name:16} rc={rc} {dt/60:.1f} min" + ("  (copied)" if dt == 0.0 else ""), flush=True)
    if "--table" in sys.argv:
        if part == "slip":
            d.table_slip()
        else:
            d.table_ch9(rule)
    if "--figure" in sys.argv:
        d.figure93(rule)


if __name__ == "__main__":
    main()
