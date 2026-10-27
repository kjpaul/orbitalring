#!/usr/bin/env python3
"""
The deployment LIM coil on the measured tape limit: the field at every tape of the coil's own winding, the worst tape's
critical current at its own field and angle at 70 K from the Robinson dataset, the layers of 3 mm tape per turn that keep
every tape at or below 50 percent of its limit, the coil current limit that follows, the hysteresis loss on that basis,
and the same for the number-of-turns trade (Claude_Code_instructions_lim_tape_limit_2026-10-26.md, items 1 to 4 and 6).
Book: The Orbital Ring, Volume IV: Chapter 11 (Table 11.2, Sections 11.5 to 11.7), Chapter 18 (Table 18.2).
Run: python tools/lim_tape_limit.py            (the field model and the tables; about a minute)
     python tools/lim_tape_limit.py --tables   (Tables 11.3 and 18.2 from the deployment logs in Vol_III/runs_700km_hw10686_2026-10-12/)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT

Corrected 2026-10-27: the pancake's tapes are on edge, face to face along the coil's surface, one turn beside the next (the
run of 2026-10-26 had them laid the other way with the wrong normal; its pancake figures are superseded, lim_hysteresis.py).

The winding.  lim_config.py gives the coil 208 turns, 2 m wide, of 3 mm tape, and budgets each turn at 1 mm of Kapton,
0.06 mm of adhesive and 0.12 mm per layer of tape (COIL_WINDING_WIDTH).  Two layouts are modelled with that budget:
  flat      the tapes lie flat, edge to edge, one turn beside the next across the coil's surface: a sheet N x (w + 1.06 mm)
            wide and 0.08 mm x layers thick (the conductor; the 0.12 mm budget per layer leaves room for the solder);
  pancake   the tapes stand on edge, face to face, one turn beside the next: a band N x (1.06 + 0.12 x layers) mm wide and w tall.
A third, the compact pancake of tools/lim_tape_current_check.py (0.05 mm of insulation per turn, the book's 109 A figure),
is carried for comparison only.  The field is the 2-D field of the straight bundle with uniform current density in each
turn, by the closed-form field of a rectangle (exact for 2-D magnetostatics with no iron; a magnetostatic_2d grid solve
of the pancake is printed beside it as a check).  The plate is gamma titanium aluminide, non-magnetic, so it adds no field
at DC; its eddy currents at the supply frequency oppose the winding's field and are neglected (conservative).  The other
coil side of the loop (100 m away) and the neighbouring phase belts (33 m away) add under 0.3 mT and are neglected.
The field is taken at the AC peak coil current.  Every tape (each layer of each turn) is sampled at its centre and both edges,
with the field of every other conductor but not its own (the measured limit already holds the tape's self-field); the
angle is that between the field and the tape's normal, folded to the lower of theta and 180 - theta, and the critical
current is read from data/THEVA_Robinson_70K_angle_dependence.csv as bearing_optimize.ic_angle_full does.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, os, re, sys
import numpy as np
import ring_altitude as ring
import lim_config as cfg
import lim_physics as phys
from lim_hysteresis import ic_tape, robinson, TAPE_T, TAPE_RHO, KAPTON, ADHESIVE, T_LAYER_BUDGET
T_COND = TAPE_T

MU0 = 4e-7 * math.pi
I_BOOK = 163.0            # A, the peak coil current of Table 11.2
SHARE = 0.5               # the bearing rule: the worst tape at or below 50 percent of its measured limit
LOGDIR = os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "..", "..", "Vol_III", "runs_700km_hw10686_2026-10-12")


# the field model, the tape data and the limits live in lim_hysteresis.py since 2026-10-27 (shared with the loss model)
from lim_hysteresis import layout, rect_field, worst_tape, coil_limit, winding_width, ROBINSON_CSV


def tape_mass(N, L, w):
    """kg of tape per LIM (3 phases x 3 pole pitches) and per site (2 LIMs)"""
    per_coil = 2 * (cfg.W_COIL + cfg.TAU_P) * N * L * w * T_COND * TAPE_RHO
    per_lim = per_coil * cfg.LIM_PHASES * cfg.PITCH_COUNT
    return per_lim, per_lim * cfg.LIMS_PER_SITE


def grid_check(N, L, w, I_coil, h=1.0e-3):
    """the pancake solved on a magnetostatic_2d grid (linear, A = 0 on a far boundary): the largest field over the bundle"""
    import magnetostatic_2d as ms
    rects, normal, bw, bh = layout("pancake", N, L, w)
    Lx, Ly = 2.2 * bw + 0.4, 16 * bh + 0.4
    g = ms.Grid(Lx, Ly, h, x0=-Lx / 2 + bw / 2, y0=-Ly / 2 + bh / 2)
    for (x1, x2, y1, y2) in rects:
        g.add_current(x1, x2, y1, y2, I_coil)
    g.solve()
    f = g.region_field(0.0, bw, 0.0, bh)
    return f, g.nx * g.ny


def hysteresis_per_hz(N, L, w, I_coil, ic_tape):
    """hysteresis loss per LIM per hertz of supply frequency, W/Hz: the book's loss factor and the Norris formula per tape"""
    params = dict(n_turns=N, w_coil=cfg.W_COIL, w_tape=w, alpha_tape=cfg.ALPHA_TAPE, l_hts_coil=2 * (cfg.W_COIL + cfg.TAU_P) * N,
                  lim_phases=cfg.LIM_PHASES, pitch_count=cfg.PITCH_COUNT, i_c=ic_tape, layers=L)
    q_book = phys.q_hts_loss_factor(I_coil, N, cfg.W_COIL, w, cfg.ALPHA_TAPE)
    q_norris = L * phys.q_hysteresis_norris(I_coil / L, ic_tape)
    scale = params["l_hts_coil"] * cfg.LIM_PHASES * cfg.PITCH_COUNT
    return q_book * scale, q_norris * scale


def field_model():
    N = cfg.N_TURNS; w = cfg.W_TAPE
    print("THE LIM COIL ON THE MEASURED TAPE LIMIT (2026-10-26)")
    print(f"  Coil: {N} turns, {cfg.W_COIL:.0f} m wide, pole pitch {cfg.TAU_P:.0f} m, {w*1e3:.0f} mm Theva AP tape at 70 K; the book's peak coil current {I_BOOK:.0f} A.")
    print(f"  Per-turn budget from lim_config.COIL_WINDING_WIDTH: {KAPTON*1e3:.2f} mm Kapton + {ADHESIVE*1e3:.2f} mm adhesive + {T_LAYER_BUDGET*1e3:.2f} mm per layer; conductor {T_COND*1e3:.2f} mm per layer.")
    ic0 = float(ic_tape(np.array([0.0]), np.array([0.0]), w)[0])
    print(f"  Tape with no field: {ic0:.0f} A per {w*1e3:.0f} mm at 70 K (Robinson); lim_config's legacy rule {cfg.IC_PER_MM_PER_LAYER:.1f} A/mm x {w*1e3:.0f} mm = {cfg.IC_PER_MM_PER_LAYER*w*1e3:.0f} A per layer.")
    print(f"  Far parts neglected: the other coil side {cfg.TAU_P:.0f} m away adds {MU0*N*I_BOOK/(2*math.pi*cfg.TAU_P)*1e3:.3f} mT, the nearest other phase belt {cfg.TAU_P/3:.0f} m away {MU0*N*I_BOOK/(2*math.pi*cfg.TAU_P/3)*1e3:.2f} mT.")
    print()
    print("1. and 2.  The field at the tape and the layers, at the 163 A coil current")
    print(f"  {'layout':<16} {'layers':>6} {'bundle':>22} {'worst tape field':>16} {'angle':>6} {'limit':>7} {'perp. limit':>11} {'per tape':>9} {'share':>6} {'coil limit at 50 %':>18} {'tape kg per LIM / site':>22}")
    chosen = {}
    for kind in ("flat", "pancake", "compact"):
        for L in (1, 2, 3, 4):
            wt = worst_tape(kind, N, L, w, I_BOOK)
            I_lim, wtl = coil_limit(kind, N, L, w)
            m_lim, m_site = tape_mass(N, L, w)
            share = I_BOOK / L / wt["ic"]
            print(f"  {kind:<16} {L:6d} {wt['width']*1e3:8.0f} x {wt['height']*1e3:6.2f} mm {wt['B']:14.4f} T {wt['theta']:6.1f} {wt['ic']:7.1f} {wt['ic_perp']:11.1f} {I_BOOK/L:9.1f} {share:6.0%} {I_lim:12.1f} A ({wtl['ic']:.0f} A, {wtl['B']:.3f} T) {m_lim:10.0f} / {m_site:.0f}")
            if kind not in chosen and share <= SHARE:
                chosen[kind] = dict(L=L, I_lim=I_lim, wt=wt, wtl=wtl, m_lim=m_lim, m_site=m_site)
    print("  (bundle = width x height of the winding cross-section; the worst tape field and angle at 163 A; the coil limit is the current at which")
    print("   the worst tape reaches 50 percent of its limit, with the limit read at that current's field; tape mass with 0.08 mm conductor at 8,890 kg/m3)")
    print()
    print("3.  The chosen layer count and the coil current limit")
    for kind, c in chosen.items():
        print(f"  {kind}: {c['L']} layer(s); at 163 A the worst tape sees {c['wt']['B']:.3f} T at {c['wt']['theta']:.1f} degrees, limit {c['wt']['ic']:.1f} A, carries {I_BOOK/c['L']:.1f} A ({I_BOOK/c['L']/c['wt']['ic']:.0%});")
        print(f"         coil current limit {c['I_lim']:.1f} A (worst tape {c['wtl']['ic']:.1f} A at {c['wtl']['B']:.3f} T, {c['wtl']['theta']:.1f} degrees, at {SHARE:.0%}); tape {c['m_lim']:,.0f} kg per LIM, {c['m_site']:,.0f} kg per site.")
    # the grid check on the pancake
    cp = chosen["pancake"]
    f, ncell = grid_check(N, cp["L"], w, I_BOOK)
    wt = worst_tape("pancake", N, cp["L"], w, I_BOOK)
    print(f"  Grid check (magnetostatic_2d, {ncell:,} cells of 1 mm, pancake with {cp['L']} layers at 163 A): largest cell-averaged field over the bundle {f['B_max']:.4f} T,")
    print(f"  mean |Bx| {f['Bx_mean']:.4f} T; the closed form's largest field at a tape {wt['Bmax']:.4f} T (a 1 mm cell cannot resolve the {cp['L']*T_COND*1e3:.2f} mm tapes, so the grid gives the")
    print(f"  bundle's smoothed field and the closed form the field each tape sees from its neighbours; the sheet limit mu0 N I / (2 width) is {MU0*N*I_BOOK/(2*wt['width'])*1e3:.1f} mT)")
    print()
    print("4.  Hysteresis on this basis (loss per LIM per hertz of supply frequency, at the 163 A coil current; the run's peak follows from its frequency)")
    print(f"  {'tape width':>10} {'layout':<9} {'layers':>6} {'worst limit':>11} {'per tape':>8} {'share':>6} {'book factor W/Hz':>16} {'Norris W/Hz':>12} {'tape kg/LIM':>11} {'coil limit':>10}")
    for wmm in (3, 4, 6, 12):
        ww = wmm * 1e-3
        for kind in ("flat", "pancake"):
            for L in (1, 2, 3, 4, 6, 8):
                wt = worst_tape(kind, N, L, ww, I_BOOK)
                if I_BOOK / L / wt["ic"] <= SHARE:
                    break
            I_lim, wtl = coil_limit(kind, N, L, ww)
            hb, hn = hysteresis_per_hz(N, L, ww, I_BOOK, wt["ic"])
            print(f"  {wmm:7d} mm {kind:<9} {L:6d} {wt['ic']:9.1f} A {I_BOOK/L:8.1f} {I_BOOK/L/wt['ic']:6.0%} {hb:16.1f} {hn:12.2f} {tape_mass(N, L, ww)[0]:11.0f} {I_lim:8.1f} A")
    print("  (book factor: B_coil x I x w x sin(20 deg) per metre of turn, independent of the layers; Norris: layers x Q(I/layers, Ic of the worst tape),")
    print("   which falls as the fourth power of the per-tape share, so it drops with the layers and with the width)")
    print()
    print("6.  The turns trade: the layers, the coil current limit and the tape mass at each N (flat layout, 3 mm tape, 50 percent rule)")
    print(f"  {'N':>5} {'layers':>6} {'coil limit':>10} {'worst tape at the limit':>24} {'ampere-turns':>12} {'tape kg/LIM':>11} {'tape kg/site':>12}")
    trade = {}
    for Nn in (100, 150, 208, 300):
        for L in (1, 2, 3, 4):
            wt = worst_tape("flat", Nn, L, w, I_BOOK)
            if I_BOOK / L / wt["ic"] <= SHARE:
                break
        I_lim, wtl = coil_limit("flat", Nn, L, w)
        m_lim, m_site = tape_mass(Nn, L, w)
        trade[Nn] = dict(L=L, I_lim=I_lim)
        print(f"  {Nn:5d} {L:6d} {I_lim:8.1f} A {wtl['ic']:8.1f} A at {wtl['B']:.3f} T, {wtl['theta']:4.1f} deg {Nn*I_lim/1e3:9.1f} kA {m_lim:11.0f} {m_site:12.0f}")
    print("  (the sheet's field is set by the current per unit width, I / pitch, so the limit hardly depends on N; the ampere-turns and the thrust do)")
    return chosen, trade


# ---------------------------------------------------------------- the tables from the deployment logs
def parse_log(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    g = lambda pat, cast=float: cast(re.search(pat, s).group(1)) if re.search(pat, s) else None
    months = []
    for m in re.finditer(r"^\s*(\d+) \|\s*([\d.]+)% \|\s*(\d+) V \|\s*([\d.]+) A \|\s*(\d+) N \|\s*([\d.]+) MW \|\s*([\d.]+) m\s+\|\s*([\d.]+) K", s, re.M):
        months.append(tuple(float(v) for v in m.groups()))
    mm = {}
    for m in re.finditer(r"^\s*(\w+)\s+([\d.]+) (\S+)\s+([\d.]+) (\w+)$", s, re.M):
        mm[m.group(1)] = (float(m.group(2)), m.group(3), float(m.group(4)), m.group(5))
    return dict(days=g(r"Deployment time: ([\d.]+) days"), E_total=g(r"Total energy delivered: ([\d.]+) EJ"), E_ke=g(r"Total energy from 1/2 m v\^2: ([\d.]+) EJ"),
                E_diff=g(r"These values differ by: ([\d.]+) EJ"), I_target=g(r"I_TARGET\s+A\s+([\d.]+)"), layers=g(r"HTS_TAPE_LAYERS\s+#\s+([\d.]+)"),
                I_c=g(r"I_C\s+A\s+([\d.]+)"), power=g(r"MAX_SITE_POWER\s+MW\s+([\d.]+)"), n_turns=g(r"N_TURNS\s+#\s+([\d.]+)"),
                T_cable_max=g(r"Max cable temperature: ([\d.]+) K"), q_cold_max=g(r"Max cryo cold-side load: ([\d.]+) kW"), months=months, mm=mm)


def when(mm, key):
    v = mm.get(key)
    return f"{v[2]:.2f} {v[3]}" if v else "?"


def table_11_3(logs):
    rows = [("Deployment time", lambda d: f"{d['days']:.1f} days ({d['days']/30.44:.1f} months)"),
            ("Peak coil current", lambda d: f"{d['mm']['current_max'][0]:.1f} A, at {when(d['mm'], 'current_max')} (cap {d['I_target']:.0f} A)"),
            ("Coil voltage at the end", lambda d: f"{d['months'][-1][2]/1e3:.0f} kV in the last month printed; peak {d['mm']['volts_max'][0]/1e3:.1f} kV at {when(d['mm'], 'volts_max')}"),
            ("Site power", lambda d: f"peak {d['mm']['site_power_max'][0]:.2f} MW at {when(d['mm'], 'site_power_max')}; by month: " + ", ".join(f"{int(m[0])}: {m[5]:.2f}" for m in d['months'][:8]) + " ..."),
            ("Peak thrust per LIM", lambda d: f"{d['mm']['thrust_max'][0]:,.0f} N at {when(d['mm'], 'thrust_max')}; last month {d['months'][-1][4]:.0f} N"),
            ("Peak thrust power per LIM", lambda d: f"{d['mm']['p_thrust_max'][0]:.2f} MW at {when(d['mm'], 'p_thrust_max')}"),
            ("Mechanical energy delivered", lambda d: f"{d['E_total']:.2f} EJ (from the kinetic energies {d['E_ke']:.2f} EJ, difference {d['E_diff']:.1e} EJ)"),
            ("Peak eddy current heating per LIM", lambda d: f"{d['mm']['p_eddy_max'][0]:.1f} kW at {when(d['mm'], 'p_eddy_max')}"),
            ("Peak hysteresis heating per LIM", lambda d: f"{d['mm']['p_hyst_max'][0]:.2f} kW at {when(d['mm'], 'p_hyst_max')}"),
            ("Peak cable temperature", lambda d: f"{d['T_cable_max']:.1f} K at {when(d['mm'], 'cable_temp_max')}, settling near {d['months'][-1][7]:.1f} K"),
            ("Peak supply frequency", lambda d: f"{d['mm']['f_supply_max'][0]:.1f} Hz")]
    print("| Quantity | " + " | ".join(l for l, _ in logs) + " |")
    print("|---|" + "---|" * len(logs))
    for name, f in rows:
        print(f"| {name} | " + " | ".join(f(d) for _, d in logs) + " |")


def table_18_2(label, d):
    """the cold loop at one site in the form of Table 18.2: the peak from the run's own peak (hysteresis 2 x p_hyst_max, the leak of
    that run's winding, the levitation coils as the remainder of the run's peak cold-side load), the steady state from the run's
    radiator width in months 12 to 15 as the chapter does, and after deployment the leak and the levitation coils alone"""
    cop = phys.calc_cryo_cop(cfg.T_LN2_SUPPLY, cfg.T_RADIATOR_HOT, cfg.CRYO_EFF)
    p_m2 = cfg.EM_HEATSINK * cfg.STEFAN_BOLTZMANN * (cfg.T_RADIATOR_HOT**4 - cfg.T_SPACE**4)
    N, L = int(d["n_turns"]), int(d["layers"])
    width = (1 + 0.06 + 0.12 * L) * N / 1000                                   # that run's COIL_WINDING_WIDTH
    leak = (cfg.Q_SUN + cfg.Q_EARTH_ALBEDO) * 2 * (cfg.TAU_P + cfg.W_COIL) * width * 3 * cfg.PITCH_COUNT * 2 * cfg.N_PLATES_PER_SIDE * cfg.COIL_MLI_EFFECTIVENESS
    hyst_peak = d["mm"]["p_hyst_max"][0] * 1e3 * cfg.LIMS_PER_SITE
    q_peak = d["q_cold_max"] * 1e3
    q_lev = q_peak - hyst_peak - leak
    w_steady = np.mean([m[6] for m in d["months"] if 12 <= m[0] <= 15]) if len(d["months"]) > 15 else d["months"][-1][6]
    q_steady = w_steady * cfg.HEATSINK_LENGTH * p_m2 / (1 + 1 / cop)
    cols = [("Peak, at %s" % when(d["mm"], "p_hyst_max"), hyst_peak, q_peak), ("Steady state during deployment (months 12 to 15)", q_steady - leak - q_lev, q_steady), ("After deployment", 0.0, leak + q_lev)]
    print(f"{label}: COP {cop:.4f} at {cfg.T_LN2_SUPPLY:.0f} K cold and {cfg.T_RADIATOR_HOT:.0f} K hot, radiator {p_m2:.0f} W/m2 over {cfg.HEATSINK_LENGTH:.0f} m; leak {leak/1e3:.2f} kW (winding {width*1e3:.0f} mm, {L} layer(s), {N} turns), levitation coils {q_lev/1e3:.1f} kW as the remainder of the run's peak cold-side load {q_peak/1e3:.1f} kW")
    print("| Quantity | " + " | ".join(c[0] for c in cols) + " |")
    print("|---|---|---|---|")
    print("| Hysteresis in the two LIMs | " + " | ".join(f"{c[1]/1e3:.1f} kW" for c in cols) + " |")
    print("| Leak through the cryostat insulation | " + " | ".join(f"{leak/1e3:.1f} kW" for c in cols) + " |")
    print("| Heat onto the levitation coils | " + " | ".join(f"{q_lev/1e3:.1f} kW" for c in cols) + " |")
    print("| Cold-side load | " + " | ".join(f"{c[2]/1e3:.1f} kW" for c in cols) + " |")
    print("| Compressor power | " + " | ".join(f"{c[2]/cop/1e6:.2f} MW" for c in cols) + " |")
    print("| Heat rejected at 400 K | " + " | ".join(f"{c[2]*(1+1/cop)/1e6:.2f} MW" for c in cols) + " |")
    print("| Cryogenic radiator area | " + " | ".join(f"{c[2]*(1+1/cop)/p_m2:,.0f} m2" for c in cols) + " |")
    print("| Cryogenic radiator width over 450 m | " + " | ".join(f"{c[2]*(1+1/cop)/p_m2/cfg.HEATSINK_LENGTH:.2f} m" for c in cols) + " |")


def tables():
    names = [("old basis, 8 MW", "lim_simulation.txt"), ("old basis, 16 MW", "lim_simulation_16MW.txt"),
             ("tape limit, 8 MW", "lim_simulation_tapelimit.txt"), ("tape limit, 16 MW", "lim_simulation_tapelimit_16MW.txt")]
    logs = [(l, parse_log(os.path.join(LOGDIR, f))) for l, f in names if os.path.exists(os.path.join(LOGDIR, f))]
    print("TABLE 11.3 FORM, FROM THE LOGS\n")
    table_11_3(logs)
    print("\nTABLE 18.2 FORM\n")
    for l, d in logs:
        table_18_2(l, d); print()
    print("THE TURNS TRADE (item 6), from the logs lim_simulation_tapelimit_N*.txt\n")
    print("| N | layers | coil limit | 100 kV reached at 16 MW | deployment at 8 MW | at 16 MW | peak thrust per LIM | tape per LIM |")
    print("|---|---|---|---|---|---|---|---|")
    for Nn in (100, 150, 208, 300):
        f8 = os.path.join(LOGDIR, f"lim_simulation_tapelimit_N{Nn}_8MW.txt" if Nn != 208 else "lim_simulation_tapelimit.txt")
        f16 = os.path.join(LOGDIR, f"lim_simulation_tapelimit_N{Nn}_16MW.txt" if Nn != 208 else "lim_simulation_tapelimit_16MW.txt")
        if not (os.path.exists(f8) and os.path.exists(f16)):
            continue
        d8, d16 = parse_log(f8), parse_log(f16)
        v = d16["mm"].get("volts_max")
        kv = f"{v[2]:.2f} {v[3]}" if v and v[0] >= 99.99e3 else "not reached"
        m_lim, _ = tape_mass(Nn, int(d8["layers"]), cfg.W_TAPE)
        print(f"| {Nn} | {int(d8['layers'])} | {d8['I_target']:.0f} A | {kv} | {d8['days']:.1f} days | {d16['days']:.1f} days | {d8['mm']['thrust_max'][0]:,.0f} N / {d16['mm']['thrust_max'][0]:,.0f} N | {m_lim:,.0f} kg |")


if __name__ == "__main__":
    if "--tables" in sys.argv:
        tables()
    else:
        field_model()
