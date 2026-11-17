#!/usr/bin/env python3
"""
Round 13, item 1: the loss table of the ruled stator, for lsm_config.py and lsm_simulation.py.
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 10.1 and 14.5.
Run: python tools/round13_stator_table.py     (about a minute; writes data/lsm_stator_ruled.json and
                                              runs/round13_logs/round13_stator_table.txt)

The stator Paul ruled on (17 November 2026) is round 12's: 6 m pole pitch, 17 turns of 30 soldered layers of 3 mm
copper-free tape on edge, end runs flat in a 1 m band set back 1 m, two sides, 70 K, on a sled whose field coil is a stepped
winding of 13 loops per pole; above 5 km/s the sled's field goes as s^0.4 and the coil current as s^0.6 (s = 5 km/s / v);
load angle 60 degrees.  Its tape current, its coil's inductance and its losses come from the tape-by-tape model of rounds 11
and 12 (tools/round11_stator_model.py, tools/round12_stator_model.py), which is too heavy to run inside the launch
simulation.  This script runs that model once along the launch's path and writes what the simulation needs:

  per meter of track of one side, per cycle, at 24 values of s from 1 (5 km/s and below) down to 1/8:
    q_hyst   J    hysteresis (strip model for the field through each tape's face, and Norris for the transport current)
    q_sold   J    the soldered layers of a turn in the field along the tapes' faces (the bearing's slab model)
    b2       m T2 tape length x mean square field through the face, for the eddy loss in the tape's metal:
                  P = pi^2 f^2 w^3 (t / rho) b2 / 6 per meter of track
  and the tape current, the layers, the coil current, the tape per meter, the coil's cyclic inductance.

The values of s follow the ruled split only (alpha = 0.4) at the design's peak current sheet.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # rounds 11 and 12's models are written on the old lsm_config.py; see their headers
import numpy as np

import lsm_config as cfg
import round11_stator_model as m
import round12_stator_model as m12

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_DESIGN = 83_294_324.34775715          # N, the design launch's peak thrust (the 83.29 MN of the book's launch tables)
ALPHA = 0.4; TAU = 6.0; N_TURNS = 17; V_CORNER = 5000.0


def main():
    st = m.Stator(TAU, "edge", 0.003, 70.0, 1.0, 1.0)
    sled = m12.sled_ratios(TAU, 13, 0.2)
    m12.apply_sled(st, sled)
    W = m.W; L_act = cfg.L_SLED * m.ACTIVE
    K1 = 2 * F_DESIGN / (st.sides * m.B0N * math.sin(m.DELTA) * W * L_act); Kb = K1 / m.KD
    f3 = m12._fields(st, Kb, m.B0N)
    ic_min = min(m12._region(a, b, Kb, c, d, 1.0, st.w_t, st.T, 1e-3)[0] for a, b, c, d in f3)
    i_t = m.SHARE * ic_min; n_t = Kb / i_t
    layers = max(1, int(round(2106.6 / i_t))); I_pk = layers * i_t; d_stack = layers * m12.PITCH
    sg = np.geomspace(1.0 / 8.0, 1.0, 24)
    lens = (W, TAU, 2 * st.d)
    rows = []
    for s in sg:
        K = Kb * s**(1 - ALPHA); B0 = m.B0N * s**ALPHA
        hy = sl = b2 = 0.0
        for (a, b, c, d), ln in zip(m12._fields(st, K, B0), lens):
            _, qh, qs, qf, Bn2, Bp, Btm, Bnm = m12._region(a, b, K, c, d, i_t * s**(1 - ALPHA), st.w_t, st.T, d_stack)
            hy += n_t * ln * qh; sl += n_t * ln * qs; b2 += n_t * ln * Bn2
        rows.append((float(s), hy, sl, b2))
    q = Kb * TAU / 6.0 / (N_TURNS * I_pk)
    tape_m, tape_kg = st.tape_per_m(n_t)
    L_coil = N_TURNS * st.phi_unit * Kb / I_pk
    out = dict(source="tools/round13_stator_table.py on tools/round12_stator_model.py; the stator ruled on 2026-11-17",
               tau_p=TAU, n_turns=N_TURNS, alpha=ALPHA, v_corner=V_CORNER, delta_deg=math.degrees(m.DELTA), k_w=m.KD, sides=st.sides, w_coil=W,
               active_fraction=m.ACTIVE, F_design=F_DESIGN, K1=K1, K_belt=Kb, tape_width=st.w_t, tape_current=i_t, tapes_per_m=n_t, layers=layers, I_peak=I_pk,
               coils_per_phase_per_pole_pitch=q, coils_facing_10km=st.sides * 3 * q * L_act / TAU, turn_length=2 * (W + TAU + 2 * st.d),
               tape_m_per_m_one_side=tape_m, tape_kg_per_m_ring=st.sides * tape_kg, L_coil_cyclic=L_coil, flux_per_turn_per_K=st.phi_unit,
               eddy_sheet_conductance=m12.G_BARE, sled_loops_per_pole=13, sled_At_per_pole=sled["At_for_0p1T"],
               s=[r[0] for r in rows], q_hyst=[r[1] for r in rows], q_sold=[r[2] for r in rows], b2=[r[3] for r in rows])
    with open(os.path.join(HERE, "data", "lsm_stator_ruled.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    lines = ["ROUND 13: THE LOSS TABLE OF THE RULED STATOR (data/lsm_stator_ruled.json)",
             f"  pole pitch {TAU:g} m, {N_TURNS} turns, {layers} layers, tape current {i_t:.2f} A, coil current {I_pk:.1f} A, {n_t:.1f} tapes per meter of track in an active side,",
             f"  K1 {K1:,.1f} A/m, belt current density {Kb:,.1f} A/m, {q:.4f} coils per phase per pole pitch per side, {out['coils_facing_10km']:,.0f} coils facing a 10 km sled,",
             f"  turn {out['turn_length']:.0f} m, tape {out['tape_kg_per_m_ring']:.2f} kg per meter of ring, coil cyclic inductance {L_coil*1e3:.3f} mH.",
             f"  {'s':>8} {'speed km/s':>10} {'hysteresis J/m':>15} {'soldered J/m':>13} {'b2 m T2':>12}"]
    for s, hy, sl, b2 in rows:
        lines.append(f"  {s:8.4f} {V_CORNER/s/1e3:10.2f} {hy:15.4f} {sl:13.5f} {b2:12.5f}")
    os.makedirs(os.path.join(HERE, "runs", "round13_logs"), exist_ok=True)
    with open(os.path.join(HERE, "runs", "round13_logs", "round13_stator_table.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
