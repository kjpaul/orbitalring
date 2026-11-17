#!/usr/bin/env python3
"""
Round 9, item 5b, part A: the sled's pulse in the mass driver's stator coils (Section 14.5).
Book: Orbital Ring Engineering, Volume III: Chapter 14, Section 14.5 ("5.74 kW on the loss factor of the 2025 book, and about 2 kW on the strip model").
Brief: Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md, Section 5b.
Run: python tools/round9_stator_pulse.py          (about a minute; writes runs/round9_logs/round9_stator_pulse.txt and .json)

What it does:
  1. runs lsm_simulation.run_lsm_simulation() as it stands (lsm_config.py's stator) and takes the peak hysteresis of the launch, where it falls, and
     the energy of the whole launch (the 2025 book's loss factor, lsm_physics.calc_hysteresis_power);
  2. rebuilds that peak by hand from the same formula, so every factor is on the page;
  3. puts the same stator winding (50 turns, five layers of 12 mm tape) through the Brandt and Indenbom strip model of lim_hysteresis.py
     (winding_loss: each tape at its own field and its measured 70 K limit, plus Norris's transport loss as an upper bound), flat and pancake;
  4. states the figures per 500 m site and the energy one pass leaves there.
It changes no file but its own two logs.
"""
import contextlib
import io
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the repository root
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # round 13: this study was run on the old 130 m stator of lsm_config.py (its launch and its
                                                   # constants); the switch keeps it reproducing now that the 6 m stator is the default
import matplotlib
matplotlib.use("Agg")
import numpy as np

import lsm_config as cfg
import lsm_physics as phys
import lsm_simulation as sim
import lim_hysteresis as lh
import lim_physics as limphys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round9_logs")
SITE = 500.0                                   # m, the spacing of the ring's sites (lsm_config.P_HVDC_MAX and tools/lsm_chapter_metrics.py take 500 m)
_lines = []


def out(s=""):
    print(s)
    _lines.append(s)


def main():
    os.makedirs(LOGS, exist_ok=True)
    res = {}
    out("=" * 140)
    out("  ROUND 9, ITEM 5b, PART A: THE SLED'S PULSE IN THE STATOR COILS (SECTION 14.5)")
    out("=" * 140)

    # ---- 1. the launch as lsm_simulation runs it today
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok = sim.run_lsm_simulation(quick_mode=False)
    P = np.array(sim.data_P_hysteresis); v = np.array(sim.data_v); I = np.array(sim.data_I_stator); f = np.array(sim.data_f_supply); t = np.array(sim.data_t)
    k = int(P.argmax()); E_tot = float(sim.data_E_hyst_cumulative[-1]); x_run = float(sim.data_x[-1]) + float(v[-1]) * cfg.DT
    n_sites_sled = cfg.L_SLED / SITE
    res["run"] = dict(ok=bool(ok), t_end_min=float(t[-1] + cfg.DT) / 60, x_run_km=x_run / 1e3, laps=x_run / cfg.L_RING, P_peak_W=float(P[k]), v_at_peak=float(v[k]), I_at_peak=float(I[k]), f_at_peak=float(f[k]),
                      t_at_peak_min=float(t[k]) / 60, E_total_J=E_tot, P_last_W=float(P[-1]), v_last=float(v[-1]), I_last=float(I[-1]), f_last=float(f[-1]))
    out(f"  1. lsm_simulation.run_lsm_simulation() as it stands (dt {cfg.DT} s): reached {'the' if ok else 'NOT the'} target, {(t[-1]+cfg.DT)/60:.1f} minutes, {x_run/1e3:,.1f} km, {x_run/cfg.L_RING:.2f} laps.")
    out(f"     Peak hysteresis over the launch (the whole 10 km sled's worth of stator): {P[k]/1e3:.2f} kW, at the last step: v = {v[k]:,.1f} m/s, coil current {I[k]:.1f} A, supply {f[k]:.2f} Hz.")
    out(f"     Hysteresis energy of the whole launch: {E_tot:.4e} J ({E_tot/3.6e6:.1f} kWh).")
    out(f"     tools/lsm_chapter_metrics.py holds these two as constants (114.70 kW and 8.488e8 J) and prints the peak over {n_sites_sled:.0f} sites: {P[k]/n_sites_sled/1e3:.2f} kW per 500 m site.")

    # ---- 2. the same peak by hand
    N = cfg.N_STATOR; Lyr = cfg.N_HTS_LAYERS; w = cfg.TAPE_WIDTH; I_pk = float(I[k]); f_pk = float(f[k])
    l_coil = cfg.TAU_P / 3.0
    B_coil = cfg.MU0 * N * I_pk / l_coil
    depth = N * 0.15e-3; perim = 2.0 * (cfg.W_COIL + depth); L_layer = N * perim
    Q = B_coil * cfg.I_C_PER_LAYER * w * cfg.SIN_ALPHA                 # J per meter of one tape per cycle
    P_coil = Lyr * L_layer * Q * f_pk
    n_coils = cfg.LIM_PHASES * cfg.N_LIM_SIDES * cfg.N_UNITS
    P_hand = P_coil * n_coils
    P_fn = phys.calc_hysteresis_power(f_pk, I_pk, cfg.get_physics_params(), cfg.N_UNITS)
    B_gap = phys.calc_B_stator(I_pk, N, cfg.W_COIL, cfg.G_GAP)
    t_pass = cfg.L_SLED / cfg.V_LAUNCH
    res["stator"] = dict(turns=N, layers=Lyr, tape_width_m=w, tau_p=cfg.TAU_P, w_coil=cfg.W_COIL, gap=cfg.G_GAP, I_c_per_layer=cfg.I_C_PER_LAYER, I_critical=cfg.I_CRITICAL, I_target=cfg.I_TARGET, I_peak=cfg.I_PEAK,
                         sin_alpha=cfg.SIN_ALPHA, l_coil=l_coil, B_coil_T=B_coil, B_gap_T=B_gap, turn_perimeter_m=perim, tape_per_layer_m=L_layer, tape_per_coil_m=Lyr * L_layer, Q_J_per_m_tape_cycle=Q,
                         q_J_per_m_turn_cycle=Lyr * Q, P_coil_W=P_coil, n_units=cfg.N_UNITS, L_unit=cfg.L_UNIT, n_coils=n_coils, P_hand_W=P_hand, P_function_W=P_fn, f_30kms=cfg.V_LAUNCH / (2 * cfg.TAU_P), T_stator_K=cfg.T_STATOR)
    out()
    out(f"  2. The same figure by hand (lsm_physics.calc_hysteresis_power, lsm_config.py):")
    out(f"     stator coil: {N} turns, {Lyr} layers of {w*1e3:.0f} mm tape in parallel, pole pitch {cfg.TAU_P:.0f} m, coil height {cfg.W_COIL:.1f} m, gap {cfg.G_GAP*1e3:.0f} mm, {cfg.T_STATOR:.0f} K in lsm_config.")
    out(f"     tape rating in lsm_config: 66.7 A/mm x 12 mm x (1 - sin 20 deg) = {cfg.I_C_PER_LAYER:.1f} A per layer, {cfg.I_CRITICAL:.0f} A per coil; I_TARGET 70 percent = {cfg.I_TARGET:.0f} A, I_PEAK 80 percent = {cfg.I_PEAK:.1f} A.")
    out(f"     coil current at the peak {I_pk:.1f} A ({I_pk/Lyr:.1f} A per tape); supply frequency v / (2 tau_p) = {cfg.V_LAUNCH/(2*cfg.TAU_P):.2f} Hz at 30 km/s ({f_pk:.2f} Hz at the last step).")
    out(f"     field the formula uses: B_coil = mu0 N I / l_coil with l_coil = tau_p / 3 = {l_coil:.2f} m: {B_coil*1e3:.3f} mT.  (The stator's field at the gap, calc_B_stator: {B_gap*1e3:.1f} mT.)")
    out(f"     loss per meter of one tape per cycle: Q = B_coil x Ic_layer x w x sin 20 deg = {Q*1e3:.3f} mJ/m; per meter of turn (5 tapes): {Lyr*Q*1e3:.2f} mJ/m.")
    out(f"     tape length the formula counts: turn perimeter 2 x (w_coil + N x 0.15 mm) = {perim:.3f} m; x {N} turns = {L_layer:.2f} m per layer, {Lyr*L_layer:.1f} m of tape per coil.")
    out(f"     per coil: {Lyr} x {L_layer:.2f} m x Q x f = {P_coil:.1f} W.  Coils counted under the sled: 3 phases x 2 sides x {cfg.N_UNITS} units (unit = 3 pole pitches + 2 m = {cfg.L_UNIT:.0f} m) = {n_coils}.")
    out(f"     total {P_hand/1e3:.2f} kW (the function itself: {P_fn/1e3:.2f} kW); over {n_sites_sled:.0f} sites of 500 m: {P_hand/n_sites_sled/1e3:.2f} kW per site.")

    # ---- 3. the strip model on the same winding
    out()
    out(f"  3. The strip model (lim_hysteresis.winding_loss: Brandt and Indenbom's perpendicular-field loss of every tape at its own field and its measured 70 K limit,")
    out(f"     Robinson data, plus Norris's transport loss per tape as an upper bound), on the same {N} turns of {Lyr} x {w*1e3:.0f} mm tape at {I_pk:.1f} A:")
    rows = []
    out(f"     {'winding':>8} {'width mm':>8} {'q perp mJ/m':>11} {'q Norris mJ/m':>13} {'q sum mJ/m':>10} {'x book factor':>13} {'worst tape Ic A':>15} {'at T':>6} {'deg':>5} {'share %':>7} {'largest field T':>15}")
    for kind in ("flat", "pancake"):
        qp, qn, det = lh.winding_loss(kind, N, Lyr, w, I_pk)
        wt = lh.worst_tape(kind, N, Lyr, w, I_pk)
        r = dict(kind=kind, width_m=wt["width"], q_perp=qp, q_norris=qn, q=qp + qn, ratio=(qp + qn) / (Lyr * Q), ic_worst=wt["ic"], B_worst=wt["B"], theta_worst=wt["theta"], share=I_pk / Lyr / wt["ic"], B_max=wt["Bmax"],
                 Bperp_max=float(det["Bperp"].max()))
        rows.append(r)
        out(f"     {kind:>8} {r['width_m']*1e3:8.0f} {qp*1e3:11.3f} {qn*1e3:13.3f} {r['q']*1e3:10.3f} {r['ratio']:13.3f} {r['ic_worst']:15.1f} {r['B_worst']:6.3f} {r['theta_worst']:5.1f} {r['share']*100:7.1f} {r['B_max']:15.3f}")
    out(f"     (q per meter of turn per cycle, the five tapes of the turn together.  'x book factor' is against the {Lyr*Q*1e3:.2f} mJ/m of part 2.  'share' is the tape's {I_pk/Lyr:.1f} A over the worst tape's limit.)")
    q_lim_book = limphys.q_hts_loss_factor(I_pk, N, cfg.W_COIL, w, math.radians(cfg.ALPHA_TAPE_FIELD))
    out(f"     For information: the deployment LIM's form of the book factor (lim_physics.q_hts_loss_factor, B = mu0 N I / w_coil = {cfg.MU0*N*I_pk/cfg.W_COIL*1e3:.1f} mT, times the coil current) gives {q_lim_book*1e3:.1f} mJ/m of turn,")
    out(f"     {q_lim_book/(Lyr*Q):.1f} times the mass driver script's form.  The two scripts do not code the 2025 book's factor the same way.")
    res["strip"] = rows; res["q_lim_book_form"] = q_lim_book

    # ---- 4. per site, and the energy of one pass
    out()
    out(f"  4. Per 500 m site at 30 km/s, and the energy one pass leaves.  The sled is {cfg.L_SLED/1e3:.0f} km long, so {n_sites_sled:.0f} sites are under it at once and each meter of stator is under it for")
    out(f"     L_sled / v = {t_pass:.4f} s.  Power per site = (power of all the coils under the sled) / {n_sites_sled:.0f}; it is the power while the sled is over the site, not a mean over the lap.")
    out(f"     Energy per site per pass = that power x {t_pass:.4f} s.  Time to the next pass at that speed: (L_ring - L_sled) / v = {(cfg.L_RING - cfg.L_SLED)/cfg.V_LAUNCH/60:.1f} minutes.")
    turn_emf = 2.0 * (cfg.TAU_P + cfg.W_COIL)
    cases = [("as lsm_physics counts the tape", perim, n_coils, ""),
             ("a turn around a coil one pole pitch long", turn_emf, n_coils, "[UNVERIFIED geometry]"),
             ("the same, 3 phases x 3 pole pitches per unit", turn_emf, 3 * n_coils, "[UNVERIFIED geometry]")]
    out(f"     {'tape length counted':>46} {'m per turn':>10} {'coils':>6} {'tape under the sled km':>22} | {'book factor':>24} | {'strip, flat':>24} | {'strip, pancake':>24}")
    out(f"     {'':>46} {'':>10} {'':>6} {'':>22} | {'kW per site':>11} {'J per pass':>12} | {'kW per site':>11} {'J per pass':>12} | {'kW per site':>11} {'J per pass':>12}")
    table = []
    for lab, per, nc, flag in cases:
        turn_m = per * N * nc                                                   # meters of turn under the sled
        row = dict(case=lab, m_per_turn=per, coils=nc, turn_m=turn_m, tape_m=turn_m * Lyr, flag=flag)
        s = f"     {lab:>46} {per:10.3f} {nc:6d} {turn_m*Lyr/1e3:22.1f} |"
        for name, q in (("book", Lyr * Q), ("flat", rows[0]["q"]), ("pancake", rows[1]["q"])):
            Pt = q * turn_m * f_pk; Ps = Pt / n_sites_sled; Es = Ps * t_pass
            row[name] = dict(P_total_W=Pt, P_site_W=Ps, E_site_J=Es)
            s += f" {Ps/1e3:11.3f} {Es:12.1f} |"
        table.append(row)
        out(s + (" " + flag if flag else ""))
    res["per_site"] = table; res["t_pass_s"] = t_pass; res["t_next_min"] = (cfg.L_RING - cfg.L_SLED) / cfg.V_LAUNCH / 60
    out(f"     Row 1 is the script's own count and is the like-for-like answer to the book's sentence.  Rows 2 and 3 are not from any script in the repository: lsm_physics.calc_EMF's")
    out(f"     own derivation has each stator coil span one pole pitch with N turns of height w_coil, which makes a turn 2 x (tau_p + w_coil) = {turn_emf:.0f} m long, and lim_config.py counts the")
    out(f"     deployment LIM's tape that way (L_HTS_COIL = 2 (W_COIL + TAU_P) N_TURNS).  calc_hysteresis_power counts {perim:.3f} m per turn, {turn_emf/perim:.1f} times less.  Which stator the mass")
    out(f"     driver has (how long a turn is, how many coils lie in a unit) is not settled by any script, so rows 2 and 3 are marked unverified and no book figure should be taken from them")
    out(f"     until the stator's winding is fixed.")
    Es_mean = E_tot / (cfg.L_RING / SITE)
    out()
    out(f"     The whole launch on the book factor: {E_tot:.4e} J over {x_run/cfg.L_RING:.2f} laps; as a mean over the ring's {cfg.L_RING/SITE:,.0f} sites, {Es_mean/1e3:.2f} kJ per site per launch (row 1's count).")
    res["E_site_mean_launch_J"] = Es_mean
    tab = lh.loss_table("flat", N, Lyr, w, 1.05 * cfg.I_PEAK)
    E_strip = float(np.sum(np.interp(I, tab["I"], tab["q"]) * f) * cfg.DT) * perim * N * n_coils
    out(f"     The whole launch on the strip model, flat winding (lim_hysteresis.loss_table at each step's coil current and frequency, row 1's count): {E_strip:.4e} J, {E_strip/E_tot:.1f} times the book factor's;")
    out(f"     as a mean over the ring's sites, {E_strip/(cfg.L_RING/SITE)/1e3:.1f} kJ per site per launch.  The coil current runs between {I.min():.0f} and {I.max():.0f} A through the launch.")
    res["E_launch_strip_flat_J"] = E_strip; res["E_site_mean_launch_strip_flat_J"] = E_strip / (cfg.L_RING / SITE); res["I_min"] = float(I.min()); res["I_max"] = float(I.max())
    out()
    out(f"  Limits of this run: the strip model's tape limits are the 70 K Robinson measurements (lsm_config has the stator at {cfg.T_STATOR:.0f} K; no 77 K angle data is in lim_hysteresis.py);")
    out(f"  the field at each tape is the winding's own two-dimensional field (the sled's {cfg.B_SLED_NOMINAL*1e3:.0f} mT field at the stator, which also cycles as the sled's poles pass, is not added);")
    out(f"  eddy current loss in the tape's metal layers and coupling loss between the five soldered tapes are not in either model.")

    with open(os.path.join(LOGS, "round9_stator_pulse.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, "round9_stator_pulse.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  written: {os.path.join(LOGS, 'round9_stator_pulse.txt')} and .json")


if __name__ == "__main__":
    main()
