#!/usr/bin/env python3
"""
Round 12, item 4: the stator the round ends on, in full, for the book's rows.
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 9.3, 9.10, 9.11, 10.1, 14.5.
Run: python tools/round12_stator_final.py     (about a minute; writes runs/round12_logs/round12_stator_final.txt and .json)

The stator: 6 m pole pitch, active sides on edge, end runs flat in a 1 m band set back 1 m, 3 mm copper-free tape, 70 K, two
sides, on a sled whose field coil is a stepped winding of 13 loops per pole; above 5 km/s the sled's field as s^0.4 and the
coil current as s^0.6.  Beside it the same stator at s^0.6 / s^0.4 (the smallest inverter).
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

import lsm_config as cfg
import ring_altitude as ring
import round11_stator_model as m
import round12_stator_model as m12
import round12_stator as r12

LOGS = r12.LOGS
_lines = []


def out(s=""):
    print(s, flush=True); _lines.append(s)


def build(alpha, m_sc=cfg.M_SPACECRAFT, F=r12.F_BOOK, P=r12.P_BOOK, L_sled=10000.0, v_f=cfg.V_LAUNCH):
    st = m.Stator(6.0, "edge", 0.003, 70.0, 1.0, 1.0)
    m12.apply_sled(st, m12.sled_ratios(6.0, 13, 0.2))
    return m12.evaluate12(st, F, P, L_sled=L_sled, m_sc=m_sc, v_f=v_f, alpha=alpha)


def main():
    res = {}
    speeds = (1000.0, 5000.0, 10000.0, 20000.0, 30000.0)
    q = m12.sled_ratios(6.0, 13, 0.2); q130 = m12.sled_ratios(130.0, 1, 0.2)
    out("=" * 150)
    out("  ROUND 12, ITEM 4: THE STATOR THE ROUND ENDS ON")
    out("=" * 150)
    for alpha in (0.4, 0.6):
        r = build(alpha); b = build(alpha, m_sc=0.0)
        E, E_low = m12.total(r); Eb, _ = m12.total(b)
        out()
        out(f"  alpha = {alpha}: pole pitch 6 m, {r['N']} turns, {r['layers']} layers of 3 mm tape ({r['d_stack']*1e3:.2f} mm stack), tape current {r['i_t']:.1f} A, coil current {r['I_pk']:.1f} A peak,")
        out(f"  {r['q']:.2f} coils per phase per pole pitch per side, {r['n_coils']:,.0f} coils facing the sled, K1 {r['K1']:,.0f} A/m (its own field mu0 K1 / 2 = {0.5*m.MU0*r['K1']*1e3:.1f} mT, mu0 K1 = {m.MU0*r['K1']*1e3:.1f} mT),")
        out(f"  tape {r['tape_kg_per_m']:.2f} kg per meter of ring, {r['tape_ring_t']/1e3:,.0f} kt; thrust {r['F_pk']/1e6:.2f} MN, power {r['P']/1e9:.2f} GW, launch {r['t']/60:.1f} min, {r['x']/1e3:,.0f} km, {r['laps']:.2f} laps.")
        out(f"      {'km/s':>5} {'f Hz':>7} {'sled field mT':>13} {'coil A':>8} {'back-EMF kV':>12} {'inductive kV':>13} {'terminal kV':>12} {'if flux adds kV':>15} {'coil MW':>8} {'coil MVA':>9} {'power factor':>12}")
        ph = m12.phasors(r, speeds)
        for p in ph:
            out(f"      {p['v']/1e3:5.0f} {p['f']:7.1f} {p['B0']*1e3:13.1f} {p['I']:8.1f} {p['E']/1e3:12.2f} {p['X']/1e3:13.2f} {p['V']/1e3:12.2f} {p['V_mag']/1e3:15.2f} {p['P']/1e6:8.2f} {p['S']/1e6:9.2f} {p['pf']:12.3f}")
        out(f"  largest terminal {r['V_max']/1e3:.1f} kV, largest back-EMF {r['E_max']/1e3:.1f} kV, largest apparent power {r['S_max']/1e9:.1f} GVA, inverter rating {r['rating_coil']/1e6:.1f} MVA a coil, {r['rating']/1e9:,.0f} GVA, {r['rating']/r['P']:.2f} x the power.")
        e = r["eddy"]
        out(f"  heat per launch: hysteresis {r['E_hyst']:.3e}, eddy (no copper) {e['no copper']['E']:.3e}, soldered layers {r['E_sold']:.3e}, total {E:.3e} J, {E/r['KE']*100:.4f} % of {r['KE']/1e15:.3f} PJ;")
        out(f"  with 20 um copper {m12.total(r, '20 um copper')[0]:.3e} J, with 40 um {m12.total(r, '40 um copper')[0]:.3e} J; braking the empty sled {Eb:.3e} J ({b['t']/60:.1f} min).")
        site = E / m.N_SITES; siteb = (E + Eb) / m.N_SITES; wpw = m.work_per_watt(70.0)
        out(f"  per site: mean {site/1e6:.2f} MJ per launch, {siteb/1e6:.2f} MJ with braking; first {r['x_c']/1e3:,.0f} km: {E_low/r['n_start']/1e6:.1f} MJ per site; {r['site_kw_c']:,.0f} kW while the sled is overhead at 5 km/s, {r['site_kw_f']:,.0f} kW at 30 km/s ({r['site_kw_f']*cfg.L_SLED/cfg.V_LAUNCH:,.0f} kJ per pass).")
        out(f"  compressor work {E*wpw:.3e} J ({E*wpw/r['KE']*100:.2f} %); launches a day on 11.9 kW: stator {r12.SPARE_W*86400/site:.1f}, with the bearing {r12.SPARE_W*86400/(site+r12.BEARING_SITE_J):.1f}, with braking and the bearing {r12.SPARE_W*86400/(siteb+r12.BEARING_SITE_J):.1f};")
        out(f"  fourteen a day: {14*site/86400/1e3:.2f} kW per site, {14*siteb/86400/1e3:.2f} kW with braking, {14*(siteb+r12.BEARING_SITE_J)/86400/1e3:.2f} kW with braking and the bearing; over the 11.9 kW spare by {max(0.0, 14*(siteb+r12.BEARING_SITE_J)/86400-r12.SPARE_W)/1e3:.2f} kW,")
        over = max(0.0, 14 * (siteb + r12.BEARING_SITE_J) / 86400 - r12.SPARE_W)
        out(f"  which is {over*wpw/1e6:.3f} MW more compressor power per site and {over*wpw*m.N_SITES/1e9:.1f} GW over the ring.")
        light = build(alpha, m_sc=500e3, F=r12.F_BOOK / 10, P=r12.P_BOOK / 10, L_sled=1000.0, v_f=12000.0)
        out(f"  light launch (500 t to 12 km/s on a 1 km sled): {light['t']/60:.1f} min, {light['x']/1e3:,.0f} km, heat {m12.total(light)[0]:.3e} J.")
        res[f"alpha_{alpha}"] = dict(design={k: v for k, v in r.items()}, brake_E=Eb, brake_t=b["t"], E=E, E_low=E_low, phasors=ph, light_E=m12.total(light)[0], light_t=light["t"])
    out()
    v1 = cfg.V_ORBIT - ring.V_GROUND_SYNC; v2 = 9800.0 - ring.V_GROUND_SYNC
    e = m.kinematics(r12.F_BOOK, r12.P_BOOK, cfg.M_SLED_HARDWARE, v2, marks=(v1, v2))
    out(f"  the empty sled at {r12.F_BOOK/1e6:.2f} MN: {r12.F_BOOK/cfg.M_SLED_HARDWARE/cfg.G_0:.3f} g; to 7,504 m/s in space {e['marks'][v1][0]/60:.1f} min, {e['marks'][v1][1]/1e3:,.0f} km; to 9,800 m/s {e['marks'][v2][0]/60:.1f} min, {e['marks'][v2][1]/1e3:,.0f} km.")
    out()
    out("  The sled's field coils (round12_stator_model.sled_ratios; one side of the sled):")
    n_new = cfg.L_SLED * m.ACTIVE / 6.0; n_old = cfg.N_UNITS * cfg.N_POLES_PER_UNIT
    out(f"      present: {n_old} poles of 130 m, one racetrack each: {q130['At_for_0p1T']/1e3:.1f} kA-turns per pole for 0.10 T at the stator, turn {q130['turn_m']:.0f} m, {q130['Am_per_m']/1e3:.1f} kA m of conductor per meter of sled.")
    out(f"      6 m poles: {n_new:,.0f} poles, 13 nested loops each: {q['At_for_0p1T']/1e3:.1f} kA-turns per pole, mean turn {q['turn_m']:.1f} m, {q['Am_per_m']/1e3:.1f} kA m of conductor per meter of sled ({q['Am_per_m']/q130['Am_per_m']:.3f} of the present).")
    out(f"      loops in all: {n_new*13:,.0f} per side against {n_old}.  The sled's budget gives its field coils 30 kg/m of the 525 kg/m (lsm_config.py); at the same conductor per ampere-meter that part becomes {30*q['Am_per_m']/q130['Am_per_m']:.1f} kg/m.")
    res["sled"] = dict(new=q, old=q130, poles_new=n_new, poles_old=n_old)
    with open(os.path.join(LOGS, "round12_stator_final.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, "round12_stator_final.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
