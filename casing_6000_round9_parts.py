#!/usr/bin/env python3
"""
The lighter casing on the present design (round 9): the cold load tables, the deployment runs and Table 16.2 (run through
casing_6000_round9.py --part=tables|run|report).  Book: Orbital Ring Engineering, Volume III: Chapter 16.
"""
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_electrical as be
import bearing_round7 as r7
import bearing_round9 as r9
import casing_6000_round9 as c6

sys.path.insert(0, os.path.join(c6.HERE, "tools"))
HW_BOOK = 10686.0; MARGIN_HW = 629.2          # kg/m: the cable's allowance, and what the design's bearing leaves of it (10,686 - 10,057)
M_FAST = 4 * 7870.0 * 1.0 * 0.20              # the four compensation plates as ruled, 1.0 m x 200 mm
TURNS = (150, 160, 165, 170, 180, 190)
CASES = ("same", "kept", "halved")
L_ATTACH, F_TROLLEY, F_ARM = 30.0, 2.24e6, 23e6


def picked():
    """The picked support plate of each sized case at 2.5 mm cells (full solves)."""
    done = c6.cache(); out = {}
    for case in ("kept", "halved"):
        best, P = c6.pick(done, case, 0.0025)
        out[case] = best
    return out


def hardware(case, pk=None):
    """(the bearing's plates and the LIM plates, the cable's hardware allowance)."""
    if case in ("design", "same"):
        return HW_BOOK - MARGIN_HW, HW_BOOK
    pk = pk or picked()
    used = 2 * pk[case]["plate_m"] + M_FAST + c6.M_LIM
    return used, used + MARGIN_HW


def table_for(case, pk):
    """bearing_cold_load's table for a 6,000 kg/m casing: the design's support (same) or the sized one."""
    import round9_cold_load as cl
    t = cl.table("mallett")
    wn = c6.M_NEW * c6.G / 2
    t["reg_W"] = t["reg_W"] * wn / t["F_normal_N_m"]              # the regulation's heat goes as the pull the 0.1 mm offset asks for, which goes as the support's pull
    if case == "same":
        t["note"] = "the design's support track under a 6,000 kg/m casing"
        return t
    v = pk[case]; T = v["cap"]["tapes_leg"]; n_t = T / 100.0
    k = v["k_lin"]; NI_n = v["normal"]["NI"]; F_n = v["normal"]["F"]; d_n = F_n / (k * NI_n**2); NI_c = v["cap"]["NI"]; F_c = v["cap"]["F"]

    def F_low(NI):
        return k * NI**2 * (1 - (1 - d_n) * NI / NI_n)

    def NI_of_F(F):
        if F <= 0:
            return 0.0
        if F >= F_n:
            return math.exp(math.log(NI_n) + (math.log(F) - math.log(F_n)) * (math.log(NI_c) - math.log(NI_n)) / (math.log(F_c) - math.log(F_n)))
        lo, hi = 0.0, NI_n
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if F_low(mid) < F: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)
    F = np.concatenate([np.linspace(0.0, 1000.0, 81)[:-1], np.linspace(1000.0, 1.3 * F_n, 400)])
    I_cap = 10 * NI_c / T
    t.update(L_coil_H=n_t**2 * 1e4 * (v["L11"] + v["L21"]), lead_aol_m=be.lead_heat(I_cap)[1], F_track_N_m=[float(x) for x in F], I_coil_A=[float(10 * NI_of_F(x) / T) for x in F],
             I_normal_A=10 * NI_n / T, F_normal_N_m=F_n, note=f"the support sized for 6,000 kg/m, the standing load {case}: plate {v['T']*1e3:.0f} mm, {n_t:.1f} turns of ten layers, cap {I_cap:,.0f} A")
    return t


def part_tables(jobs=None):
    pk = picked()
    for case in CASES:
        t = table_for(case, pk); p = os.path.join(r9.LOGS, f"bearing_cold_load_casing6000_{case}.json")
        json.dump(t, open(p, "w", encoding="utf-8"), indent=1)
        print(f"  {case}: {t['note']}; coil {t['L_coil_H']:.2f} H, normal {t['I_normal_A']:.0f} A, pump {t['pump_V']/t['L_coil_H']*86400:.2f} A a day, regulation {t['reg_W']:.0f} W, hardware {hardware(case, pk)[1]:,.0f} kg/m")


def run_dir(case, N, lim_mm=200):
    return os.path.join(c6.BASE, f"{case}_N{N}" + ("" if lim_mm == 200 else f"_lim{lim_mm}"))


def run_one(a):
    case, N, lim_mm, hw = a
    out = run_dir(case, N, lim_mm); os.makedirs(out, exist_ok=True)
    if os.path.exists(os.path.join(out, "rising_radius_summary.json")):
        return case, N, lim_mm, 0, 0.0
    tpath = os.path.abspath(os.path.join(r9.LOGS, f"bearing_cold_load_casing6000_{case}.json")).replace("\\", "/")
    hw_run = hw - (1872.0 * (1 - lim_mm / 200.0))
    over = {"MAX_SITE_POWER": 12e6, "N_TURNS": N, "LEV_COIL_WIDTH": 4.0, "Q_COLD_EXTRA_PER_SITE": 0.0, "BEARING_COLD_LOAD": tpath,
            "D_LEV": (hw - 1872.0) / (4 * 7870.0 * 1.4), "T_PLATE": lim_mm / 1e3}
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUTF8="1", ORBITAL_RING_OUT=out, ORBITAL_RING_SLIP_RULE="optimal", ORBITAL_RING_HW_KGM=f"{hw_run:.1f}")
    code = (f"import lim_config as c; c.configure(**{over!r}); import runpy, sys; sys.argv=['lim_simulation.py', '--rising-radius', '--m_load=6000']; runpy.run_path('lim_simulation.py', run_name='__main__')")
    t0 = time.time()
    with open(os.path.join(out, "rising_radius_lim_simulation.txt"), "w", encoding="utf-8", errors="replace") as lf:
        rc = subprocess.run([sys.executable, "-c", code], cwd=c6.HERE, env=env, stdout=lf, stderr=subprocess.STDOUT).returncode
    return case, N, lim_mm, rc, time.time() - t0


def part_run(jobs=6):
    import concurrent.futures
    pk = picked()
    tasks = [(c, 165, 200, hardware(c, pk)[1]) for c in CASES] + [("kept", 165, 100, hardware("kept", pk)[1])] + [(c, N, 200, hardware(c, pk)[1]) for N in TURNS if N != 165 for c in CASES]
    print(f"  {len(tasks)} runs, {jobs} at a time, into {c6.BASE}", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        for case, N, lim, rc, dt in ex.map(run_one, tasks):
            print(f"    {case} N{N} LIM plate {lim} mm: rc={rc} {dt/60:.1f} min", flush=True)


def load_run(case, N=165, lim_mm=200, folder=None):
    import rising_radius_design_165 as d
    f = folder or run_dir(case, N, lim_mm)
    p = os.path.join(f, "rising_radius_summary.json")
    if not os.path.exists(p):
        return None
    J = json.load(open(p)); J.update(d.parse_log(os.path.join(f, "rising_radius_lim_simulation.txt")))
    b = os.path.join(f, "bearing_cold_load.json")
    if os.path.exists(b):
        J["bcl"] = json.load(open(b))
    return J


def part_report(jobs=None):
    pk = picked(); done = c6.cache()
    c6.show(done, 0.005); c6.show(done, 0.0025)
    half = json.load(open(os.path.join(r7.OUT, "round7_logs", "round7_half_2p5.json"), encoding="utf-8")); fast = r7.load("round7_fast_track.json")
    R8 = json.load(open(os.path.join(r7.OUT, "round8_logs", "round8_rows.json"), encoding="utf-8")); mass8 = R8["mass"]
    des = load_run("design", folder=os.path.join(r7.OUT, "round9_design", "slip_optimal", "N165_P12MW"))
    cols = ["design"] + list(CASES); rows = []; out = {}

    def add(lab, f):
        vals = [f(c) for c in cols]; rows.append((lab, vals)); out[lab] = vals
    m_load = lambda c: 12000.0 if c == "design" else 6000.0
    D = {c: ring.design(m_load=m_load(c), m_hw=hardware(c, pk)[1]) for c in cols}
    J = {c: (des if c == "design" else load_run(c)) for c in cols}
    sup = {"design": dict(T=0.10, plate_m=7870 * 1.2 * 0.10, sat=half["cap"]["sat"], tape=half["cap"]["tape"], NI=half["cap"]["NI"], F=half["cap"]["F"], margin=341.1 / 245.5)}
    sup["same"] = sup["design"]
    for c in ("kept", "halved"):
        v = pk[c]; sup[c] = dict(T=v["T"], plate_m=v["plate_m"], sat=v["cap"]["sat"], tape=v["cap"]["tape"], NI=v["cap"]["NI"], F=v["cap"]["F"], margin=v["rating"]["margin"])
    tape = lambda c: 2 * sup[c]["tape"] + 4 * fast["tape"]
    others = mass8["total"] - 2 * mass8["tape"]                       # the stores, switches and fast leads of round 8's weighing
    add("Standing load allowed per support track", lambda c: "99 kN/m" if c == "halved" else "198 kN/m")
    add("The support track's requirement", lambda c: f"{sup[c]['F']/1e3:.1f} kN/m")
    add("Hardware on the cable (the allowance)", lambda c: f"{hardware(c, pk)[1]:,.0f} kg/m")
    add("The sized hardware itself (the plates of the bearing and the LIM plates)", lambda c: f"{hardware(c, pk)[0]:,.0f} kg/m")
    add("Structural cable", lambda c: f"{D[c]['m_struct']:,.0f} kg/m")
    add("Cable with hardware", lambda c: f"{D[c]['m_cable']:,.0f} kg/m")
    add("Cable velocity", lambda c: f"{D[c]['v_cable']:,.0f} m/s")
    add("Tension", lambda c: f"{D[c]['tension']/1e9:.1f} GN")
    add("Deployment energy", lambda c: f"{D[c]['E_deploy']/1e18:.2f} EJ")
    ok = lambda c: J[c] is not None
    add("Deployment at 12 MW per site, 165 turns", lambda c: f"{J[c]['t_end_days']:.1f} days" if ok(c) else "not run")
    add("The day 100 kV is reached", lambda c: ("never" if J[c].get("volt_limit_day") is None else f"day {J[c]['volt_limit_day']:.0f}") if ok(c) else "not run")
    add("Peak thrust per LIM", lambda c: f"{J[c]['thrust_max_N']/1e3:.2f} kN" if ok(c) else "not run")
    add("First-day eddy heating per LIM; peak cable temperature", lambda c: f"{J[c]['p_eddy_max_kW']:.0f} kW; {J[c]['cable_temp_max_K']:.0f} K" if ok(c) else "not run")
    add("Peak cold load per site; the support's leads come off", lambda c: (f"{J[c]['q_cryo_cold_max_kW']:.1f} kW; " + (f"day {J[c]['bcl']['leads_off_day']:.0f}" if J[c].get("bcl") and J[c]["bcl"]["leads_off_day"] is not None else "not in the run")) if ok(c) else "not run")

    def turns(c):
        if c == "design":
            return "yes (round 9's sweep, see item 2)"
        pts = [(N, load_run(c, N)) for N in TURNS]; pts = [(N, x["t_end_days"]) for N, x in pts if x]
        if len(pts) < 3:
            return "not run"
        fl = min(x for _, x in pts); nb = [N for N, x in pts if x == fl][0]; d165 = dict(pts).get(165)
        return f"{'yes' if d165 <= 1.02 * fl else 'NO'}: fastest {nb} turns, {fl:.1f} days; 165 turns {(d165/fl-1)*100:+.1f} percent (" + ", ".join(f"{N}: {x:.1f}" for N, x in pts) + ")"
    add("Is 165 turns within 2 percent of the fastest", turns)
    add("Support plates (two), 1.2 m wide", lambda c: f"{sup[c]['T']*1e3:.0f} mm, {2*sup[c]['plate_m']:,.0f} kg/m, {sup[c]['sat']*100:.0f} percent saturated at the cap, rated {sup[c]['margin']:.2f} times")
    add("Compensation plates (four), 1.0 m wide", lambda c: f"200 mm, {M_FAST:,.0f} kg/m, {fast['sat']*100:.1f} percent saturated at 196.5 kN/m")
    add("LIM reaction plates (two)", lambda c: "200 mm, 1,872 kg/m")
    add("Bearing on the casing: tape", lambda c: f"{tape(c):.0f} kg/m")
    add("with formers and cryostats (one tape mass more)", lambda c: f"{2*tape(c):,.0f} kg/m")
    add("with stores, switches and fast leads", lambda c: f"{2*tape(c)+others:,.0f} kg/m")
    att = lambda c: f"{math.ceil(F_TROLLEY/(0.1*m_load(c)*ring.G_NET*L_ATTACH))}; {F_ARM/(0.1*m_load(c)*ring.G_NET*L_ATTACH):.0f}"
    add("Anchor attachments per trolley line; per spreader arm", att)
    add("Standing cryogenic load per site: the cable's heat onto 4.0 m of coils", lambda c: "10.0 kW")
    print("\n  TABLE 16.2 ON THE PRESENT DESIGN")
    print("  | Quantity | 12,000 kg/m, the design case | 6,000 kg/m, plates unchanged | 6,000 kg/m, plates sized, standing load kept | 6,000 kg/m, plates sized, standing load halved |")
    print("  |---|---|---|---|---|")
    for lab, vals in rows:
        print(f"  | {lab} | " + " | ".join(vals) + " |")
    # addendum A3: what a kilogram per meter of casing costs
    print("\n  WHAT A KILOGRAM OF CASING COSTS (each case against the design, 6,000 kg/m taken off the casing)")
    print("  | Case | Cable with hardware saved | per kg/m of casing | Structural cable saved | Hardware saved | Days saved (165 turns) | per kg/m of casing | Deployment energy saved | per kg/m of casing |")
    print("  |---|---|---|---|---|---|---|---|---|")
    a3 = {}
    names = dict(same="6,000 kg/m, plates unchanged", kept="6,000 kg/m, plates sized, standing load kept", halved="6,000 kg/m, plates sized, standing load halved")
    for c in CASES:
        dm = D["design"]["m_cable"] - D[c]["m_cable"]; ds = D["design"]["m_struct"] - D[c]["m_struct"]; dh = hardware("design", pk)[1] - hardware(c, pk)[1]
        dE = D["design"]["E_deploy"] - D[c]["E_deploy"]; dd = (J["design"]["t_end_days"] - J[c]["t_end_days"]) if ok(c) and J["design"] else None
        a3[c] = dict(cable=dm, struct=ds, hw=dh, E=dE, days=dd)
        print(f"  | {names[c]} | {dm:,.0f} kg/m | {dm/6000:.2f} kg/m | {ds:,.0f} kg/m | {dh:,.0f} kg/m | " + (f"{dd:.1f} | {dd/6000*1000:.1f} days per tonne per meter" if dd is not None else "not run | not run")
              + f" | {dE/1e18:.2f} EJ | {dE/6000/1e15:.2f} PJ |")
    out["a3"] = a3
    # the sentences of 16.6
    Jl = load_run("kept", 165, 100); Jk = J["kept"]
    sent = {}
    if Jl and Jk:
        sent["lim100"] = dict(eddy=Jl["p_eddy_max_kW"], eddy200=Jk["p_eddy_max_kW"], T=Jl["cable_temp_max_K"], T200=Jk["cable_temp_max_K"], v=Jl.get("volt_limit_day"), v200=Jk.get("volt_limit_day"),
                              days=Jl["t_end_days"], days200=Jk["t_end_days"], saved=936.0, cable=D["kept"]["m_cable"])
        print(f"\n  The LIM plate halved to 100 mm (on the sized case, standing load kept): first-day eddy heating {Jl['p_eddy_max_kW']:.0f} kW per LIM against {Jk['p_eddy_max_kW']:.0f}; the cable {Jl['cable_temp_max_K']:.0f} K against {Jk['cable_temp_max_K']:.0f};")
        print(f"  100 kV on day {Jl.get('volt_limit_day')} against {Jk.get('volt_limit_day')}; {Jl['t_end_days']:.1f} days against {Jk['t_end_days']:.1f} ({Jl['t_end_days']-Jk['t_end_days']:+.1f}); 936 kg/m saved on a cable of {D['kept']['m_cable']:,.0f} kg/m.")
    print(f"\n  The support plates: {sup['design']['T']*1e3:.0f} mm in the design, {sup['kept']['T']*1e3:.0f} mm with the standing load kept, {sup['halved']['T']*1e3:.0f} mm with it halved.  The bearing's tape: {tape('design'):.0f}, {tape('kept'):.0f} and {tape('halved'):.0f} kg/m of casing")
    print(f"  ({tape('same'):.0f} with the plates unchanged, since the design's winding is kept).  The standing load is {198e3/sup['kept']['F']*100:.0f} percent of the sized track's requirement.")
    print(f"  The sized hardware: {hardware('kept', pk)[0]:,.0f} kg/m (kept) and {hardware('halved', pk)[0]:,.0f} kg/m (halved) against the design's {hardware('design')[0]:,.0f}: {(hardware('kept', pk)[0]/hardware('design')[0]-1)*100:+.1f} and {(hardware('halved', pk)[0]/hardware('design')[0]-1)*100:+.1f} percent.")
    known = lambda c: 2 * tape(c) + others + 1200 + 62 + 65 + 5.2 + 16
    print(f"  Table 16.1's known items at 6,000 kg/m: the bearing {2*tape('kept')+others:,.0f} kg/m (kept; {2*tape('halved')+others:,.0f} halved; {2*tape('design')+others:,.0f} in the design), the array 1,200, the panels 41 to 83 (62 taken),")
    print(f"  the ring lines 65, the LIM tape 5.2, the stowed lines 16: about {known('kept'):,.0f} kg/m known, {6000-known('kept'):,.0f} left for the unweighed items (design: {known('design'):,.0f} known, {12000-known('design'):,.0f} left).")
    r9.save("casing6000_report.json", dict(rows=out, cols=cols, sup=sup, tape={c: tape(c) for c in cols}, others=others, hardware={c: hardware(c, pk) for c in cols}, sent=sent,
                                             known={c: known(c) for c in cols}))
