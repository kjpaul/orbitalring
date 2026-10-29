#!/usr/bin/env python3
"""
The LIM stator hysteresis on the Brandt-Indenbom strip model, compared with the book's factor and Norris, and the width,
winding and turns trade on it (Claude_Code_instructions_lim_hysteresis_model_2026-10-27.md, items 2 to 4).
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Sections 9.5 to 9.7), Chapter 14 (Table 14.2).
Run: python tools/lim_hysteresis_model.py            (the checks, the comparison at the design point, the trade table; a few minutes)
     python tools/lim_hysteresis_model.py --tables   (Tables 9.3 and 14.2 from the deployment logs of the candidates)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT

Reference operating conditions, from the 208-turn deployment on the measured tape limit (lim_simulation_tapelimit.txt): the
hysteresis peak at 1.75 months, where the coil current is near its cap and the supply frequency about 10.6 Hz (the book factor's
14.4 kW per LIM over its 1,361 W/Hz at 163 A); the steady state of months 12 to 15 at 72 A and about 36 Hz (progress 76 to
86 percent, v_rel / 2 tau_p); the end of deployment at 66 A and 44 Hz.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, os, sys
import numpy as np
import lim_config as cfg
import lim_physics as phys
import lim_hysteresis as lh

LOGDIR = os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "..", "..", "Vol_III", "runs_700km_hw10686_2026-10-12")
F_PEAK, F_STEADY, F_END = 10.6, 36.0, 44.15
I_STEADY, I_END = 72.0, 66.0
SHARE = 0.5


def per_lim_per_hz(q_per_m_turn, N):
    """W per hertz per LIM from a loss per metre of turn per cycle"""
    return q_per_m_turn * 2 * (cfg.W_COIL + cfg.TAU_P) * N * cfg.LIM_PHASES * cfg.PITCH_COUNT


def book_q(I, N, w):
    return phys.q_hts_loss_factor(I, N, cfg.W_COIL, w, cfg.ALPHA_TAPE)


def layers_by_rule(kind, N, w, I0=163.0):
    for L in (1, 2, 3, 4, 6, 8):
        wt = lh.worst_tape(kind, N, L, w, I0)
        if I0 / L / wt["ic"] <= SHARE:
            return L, wt
    return L, wt


def comparison():
    N, w = 208, 3e-3
    print("3.  THE THREE MODELS AT THE DESIGN POINT: 208 turns, 3 mm tape, per LIM")
    for kind, L in (("pancake", 3), ("flat", 2)):
        I_cap, wtl = lh.coil_limit(kind, N, L, w)
        print(f"  {kind} winding, {L} layers (the 50 percent rule: cap {I_cap:.1f} A, worst tape {wtl['ic']:.1f} A at {wtl['B']:.3f} T, {wtl['theta']:.1f} deg from its normal)")
        print(f"  {'coil current':>13} {'book W/Hz':>10} {'norris W/Hz':>11} {'brandt perp W/Hz':>16} {'+ norris':>9} {'brandt / book':>13} {'end tapes share':>15}")
        for I in (163.0, I_cap, 184.6, I_STEADY, I_END):
            qp, qn, d = lh.winding_loss(kind, N, L, w, I)
            Ic_w = lh.worst_tape(kind, N, L, w, I)["ic"]
            qb = book_q(I, N, w); qnor = L * phys.q_hysteresis_norris(I / L, Ic_w)
            n = d["n"]; ends = np.r_[0:10 * L, n - 10 * L:n]
            share = d["qp"][ends].sum() / max(d["qp"].sum(), 1e-30)
            print(f"  {I:11.1f} A {per_lim_per_hz(qb, N):10.0f} {per_lim_per_hz(qnor, N):11.1f} {per_lim_per_hz(qp, N):16.0f} {per_lim_per_hz(qp + qn, N):9.0f} {(qp + qn) / qb:13.2f} {share:15.0%}")
        for label, I, f in (("the hysteresis peak, 1.75 months", min(I_cap, 184.6), F_PEAK), ("months 12 to 15", I_STEADY, F_STEADY), ("the end of deployment", I_END, F_END)):
            qp, qn, d = lh.winding_loss(kind, N, L, w, I)
            qb = book_q(I, N, w)
            print(f"    at {label}: {I:.0f} A, {f:.1f} Hz: book {per_lim_per_hz(qb, N) * f / 1e3:.1f} kW, brandt {per_lim_per_hz(qp + qn, N) * f / 1e3:.1f} kW per LIM "
                  f"(perpendicular {per_lim_per_hz(qp, N) * f / 1e3:.1f}, Norris {per_lim_per_hz(qn, N) * f / 1e3:.2f}); hottest tape {d['qp'].max() * f:.2f} W/m")
        print("  (end tapes share: the loss in the ten tapes at each end of the bundle, of the perpendicular-field loss of all the tapes;")
        print("   the Norris column is the self-field transport loss alone, the brandt column the strip loss in the winding's field with Norris added as a bound)")
        print()


def trade():
    cop = phys.calc_cryo_cop(cfg.T_LN2_SUPPLY, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF)
    p_m2 = cfg.EM_HEATSINK * cfg.STEFAN_BOLTZMANN * (cfg.T_RADIATOR_HOT**4 - cfg.T_SPACE**4)
    print("4.  THE TRADE ON THE BRANDT MODEL (3 mm at 100, 150, 208, 300 turns; 4, 6 and 12 mm at 150 and 208; flat and pancake; layers by the 50 percent rule)")
    print(f"  cold-side load of a site = hysteresis of the two LIMs + the cryostat leak of that winding's width + 10 kW on the levitation coils;")
    print(f"  peak = the cap current at {F_PEAK} Hz, steady = {I_STEADY:.0f} A at {F_STEADY:.0f} Hz (the design run's conditions, item 3); compressor = load / COP {cop:.4f}")
    print(f"  {'w':>5} {'winding':<8} {'N':>4} {'L':>2} {'cap A':>6} {'worst tape':>26} {'tape kg/LIM':>11} {'band mm':>8} {'leak kW':>7} {'book W/Hz':>10} {'brandt W/Hz at cap':>18} {'brandt W/Hz at 72 A':>19} {'peak kW':>8} {'steady kW':>9} {'MW peak':>8} {'MW steady':>9}")
    rows = []
    for wmm in (3, 4, 6, 12):
        w = wmm * 1e-3
        for kind in ("flat", "pancake"):
            for N in ((100, 150, 208, 300) if wmm == 3 else (150, 208)):
                L, wt = layers_by_rule(kind, N, w)
                I_cap, wtl = lh.coil_limit(kind, N, L, w)
                qp, qn, _ = lh.winding_loss(kind, N, L, w, I_cap)
                qps, qns, _ = lh.winding_loss(kind, N, L, w, I_STEADY)
                qb = book_q(I_cap, N, w)
                band = lh.winding_width(kind, N, L, w)
                leak = (cfg.Q_SUN + cfg.Q_EARTH_ALBEDO) * 2 * (cfg.TAU_P + cfg.W_COIL) * band * 3 * cfg.PITCH_COUNT * 2 * cfg.N_PLATES_PER_SIDE * cfg.COIL_MLI_EFFECTIVENESS
                hp = per_lim_per_hz(qp + qn, N) * F_PEAK * cfg.LIMS_PER_SITE; hs = per_lim_per_hz(qps + qns, N) * F_STEADY * cfg.LIMS_PER_SITE
                q_peak = hp + leak + 10e3; q_steady = hs + leak + 10e3
                m = lh.tape_mass_per_lim(N, L, w, cfg.W_COIL, cfg.TAU_P, cfg.LIM_PHASES, cfg.PITCH_COUNT)
                rows.append(dict(w=wmm, kind=kind, N=N, L=L, cap=I_cap, m=m, band=band, leak=leak, hp=hp, hs=hs, q_peak=q_peak, q_steady=q_steady, wtl=wtl))
                print(f"  {wmm:3d} mm {kind:<8} {N:4d} {L:2d} {I_cap:6.1f} {wtl['ic']:6.1f} A at {wtl['B']:.3f} T, {wtl['theta']:4.1f} deg {m:11,.0f} {band*1e3:8.0f} {leak/1e3:7.2f} {per_lim_per_hz(qb, N):10.0f} {per_lim_per_hz(qp + qn, N):18.0f} {per_lim_per_hz(qps + qns, N):19.1f} {q_peak/1e3:8.1f} {q_steady/1e3:9.1f} {q_peak/cop/1e6:8.2f} {q_steady/cop/1e6:9.2f}")
    print("  (cap = the coil current at which the worst tape reaches 50 percent of its measured limit; band = the winding's width across the coil's surface, the cryostat band;")
    print("   the ampere-turns at the cap, N x cap, set the thrust per coil and the voltage; the deployment runs of the candidates give the times)")
    return rows


def parse_log(path):
    import re
    s = open(path, encoding="utf-8", errors="replace").read()
    g = lambda pat, cast=float: cast(re.search(pat, s).group(1)) if re.search(pat, s) else None
    months = [tuple(float(v) for v in m.groups()) for m in re.finditer(r"^\s*(\d+) \|\s*([\d.]+)% \|\s*(\d+) V \|\s*([\d.]+) A \|\s*(\d+) N \|\s*([\d.]+) MW \|\s*([\d.]+) m\s+\|\s*([\d.]+) K", s, re.M)]
    mm = {m.group(1): (float(m.group(2)), m.group(3), float(m.group(4)), m.group(5)) for m in re.finditer(r"^\s*(\w+)\s+([\d.]+) (\S+)\s+([\d.]+) (\w+)$", s, re.M)}
    return dict(days=g(r"Deployment time: ([\d.]+) days"), E_total=g(r"Total energy delivered: ([\d.]+) EJ"), E_diff=g(r"These values differ by: ([\d.]+) EJ"),
                I_target=g(r"I_TARGET\s+A\s+([\d.]+)"), layers=g(r"HTS_TAPE_LAYERS\s+#\s+([\d.]+)"), width=g(r"HTS_TAPE_WIDTH\s+mm\s+([\d.]+)"), power=g(r"MAX_SITE_POWER\s+MW\s+([\d.]+)"),
                n_turns=g(r"N_TURNS\s+#\s+([\d.]+)"), model=g(r"HYSTERESIS_MODEL\s+(\w+)", str), T_cable_max=g(r"Max cable temperature: ([\d.]+) K"), q_cold_max=g(r"Max cryo cold-side load: ([\d.]+) kW"), months=months, mm=mm)


def when(mm, key):
    v = mm.get(key)
    return f"{v[2]:.2f} {v[3]}" if v else "?"


def tables():
    import glob
    cop = phys.calc_cryo_cop(cfg.T_LN2_SUPPLY, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF)
    p_m2 = cfg.EM_HEATSINK * cfg.STEFAN_BOLTZMANN * (cfg.T_RADIATOR_HOT**4 - cfg.T_SPACE**4)
    files = sorted(glob.glob(os.path.join(LOGDIR, "lim_simulation_brandt_*.txt")))
    logs = [(os.path.basename(f)[len("lim_simulation_brandt_"):-4], parse_log(f)) for f in files]
    print("TABLE 11.3 FORM, THE CANDIDATES ON THE BRANDT MODEL\n")
    print("| Quantity | " + " | ".join(l for l, _ in logs) + " |")
    print("|---|" + "---|" * len(logs))
    rows = [("Winding", lambda d: f"{d['width']:.0f} mm, {int(d['n_turns'])} turns, {int(d['layers'])} layer(s), cap {d['I_target']:.0f} A, {d['power']:.0f} MW"),
            ("Deployment time", lambda d: f"{d['days']:.1f} days ({d['days']/30.44:.1f} months)"),
            ("Peak coil current", lambda d: f"{d['mm']['current_max'][0]:.1f} A at {when(d['mm'], 'current_max')}"),
            ("Coil voltage", lambda d: f"peak {d['mm']['volts_max'][0]/1e3:.1f} kV at {when(d['mm'], 'volts_max')}"),
            ("Site power", lambda d: f"peak {d['mm']['site_power_max'][0]:.2f} MW at {when(d['mm'], 'site_power_max')}; last month {d['months'][-1][5]:.2f} MW"),
            ("Peak thrust per LIM", lambda d: f"{d['mm']['thrust_max'][0]:,.0f} N at {when(d['mm'], 'thrust_max')}; last month {d['months'][-1][4]:.0f} N"),
            ("Peak thrust power per LIM", lambda d: f"{d['mm']['p_thrust_max'][0]:.2f} MW at {when(d['mm'], 'p_thrust_max')}"),
            ("Mechanical energy", lambda d: f"{d['E_total']:.2f} EJ (balance {d['E_diff']:.1e} EJ)"),
            ("Peak eddy heating per LIM", lambda d: f"{d['mm']['p_eddy_max'][0]:.1f} kW at {when(d['mm'], 'p_eddy_max')}"),
            ("Peak hysteresis per LIM", lambda d: f"{d['mm']['p_hyst_max'][0]:.2f} kW at {when(d['mm'], 'p_hyst_max')}"),
            ("Peak cold-side load, site", lambda d: f"{d['q_cold_max']:.1f} kW at {when(d['mm'], 'q_cryo_cold_max')}; compressor {d['q_cold_max']*1e3/cop/1e6:.2f} MW"),
            ("Peak cryo compressor (run)", lambda d: f"{d['mm']['p_cryo_max'][0]:.2f} MW at {when(d['mm'], 'p_cryo_max')}"),
            ("Steady radiator width, months 12 to 15", lambda d: f"{np.mean([m[6] for m in d['months'] if 12 <= m[0] <= 15]):.2f} m" if len(d['months']) > 15 else f"{d['months'][-1][6]:.2f} m (last month)"),
            ("Peak cable temperature", lambda d: f"{d['T_cable_max']:.1f} K")]
    for name, f in rows:
        print(f"| {name} | " + " | ".join(f(d) for _, d in logs) + " |")
    print("\nTABLE 18.2 FORM (peak from the run, steady state from the run's radiator width in months 12 to 15, levitation coils 10.0 kW)\n")
    for l, d in logs:
        N, L, w = int(d["n_turns"]), int(d["layers"]), d["width"] * 1e-3
        kind = "flat" if N * (w + 1.06e-3) * 1e3 < 1e9 and ("flat" in l) else "pancake"
        band = lh.winding_width(kind, N, L, w)
        leak = (cfg.Q_SUN + cfg.Q_EARTH_ALBEDO) * 2 * (cfg.TAU_P + cfg.W_COIL) * band * 3 * cfg.PITCH_COUNT * 2 * cfg.N_PLATES_PER_SIDE * cfg.COIL_MLI_EFFECTIVENESS
        hyst_peak = d["mm"]["p_hyst_max"][0] * 1e3 * cfg.LIMS_PER_SITE; q_peak = d["q_cold_max"] * 1e3; q_lev = q_peak - hyst_peak - leak
        w_steady = np.mean([m[6] for m in d["months"] if 12 <= m[0] <= 15]) if len(d["months"]) > 15 else d["months"][-1][6]
        q_steady = w_steady * cfg.HEATSINK_LENGTH * p_m2 / (1 + 1 / cop)
        print(f"{l}: leak {leak/1e3:.2f} kW ({kind}, band {band*1e3:.0f} mm), levitation coils {q_lev/1e3:.1f} kW")
        print("| Quantity | Peak, at %s | Steady state, months 12 to 15 | After deployment |" % when(d["mm"], "p_hyst_max"))
        print("|---|---|---|---|")
        for name, vals in (("Hysteresis in the two LIMs", (hyst_peak, q_steady - leak - q_lev, 0.0)), ("Leak through the cryostat insulation", (leak,) * 3), ("Heat onto the levitation coils", (q_lev,) * 3), ("Cold-side load", (q_peak, q_steady, leak + q_lev))):
            print(f"| {name} | " + " | ".join(f"{v/1e3:.1f} kW" for v in vals) + " |")
        for name, fn in (("Compressor power", lambda q: f"{q/cop/1e6:.2f} MW"), ("Heat rejected at 400 K", lambda q: f"{q*(1+1/cop)/1e6:.2f} MW"), ("Cryogenic radiator area", lambda q: f"{q*(1+1/cop)/p_m2:,.0f} m2"), ("Cryogenic radiator width over 450 m", lambda q: f"{q*(1+1/cop)/p_m2/cfg.HEATSINK_LENGTH:.2f} m")):
            print(f"| {name} | " + " | ".join(fn(q) for q in (q_peak, q_steady, leak + q_lev)) + " |")
        print()


if __name__ == "__main__":
    if "--tables" in sys.argv:
        tables()
    else:
        print("THE LIM STATOR HYSTERESIS ON THE BRANDT-INDENBOM STRIP MODEL (2026-10-27)\n")
        lh.checks(); print()
        comparison()
        trade()
