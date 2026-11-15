#!/usr/bin/env python3
"""
The bearing, round 6, part B (2026-11-14, run 2026-11-15): Vol_III/Claude_Code_instructions_bearing_round6b_2026-11-14.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 14 (the cold loop), Chapter 15 (the deployment).
Run (from Python/orbitalring):
  python bearing_round6b.py --part=legs      S1.1  where the legs are in the plate solve, the coupling run and the roll run; the coupling rerun on the solve's legs
  python bearing_round6b.py --part=dump      S1.2  the support coil's dump with D2 and the strip across the coil: the numbers of the corrected circuit
  python bearing_round6b.py --part=film      S1.3  the order of the ten turns in a fast coil's cell, the voltage across every film, the cell and the bundle
  python bearing_round6b.py --part=plates    S2.1  the trade table of plates (the cached 2.5 mm solves, the tape's pulse loss and the plate heat of a launch on each)
  python bearing_round6b.py --part=bundle    S2.2  one attempt at the bundle: shapes and stack orientations at 5 mm cells  [--jobs=4]
  python bearing_round6b.py --part=lost      S2.3  the lost support track with ten families connected
  python bearing_round6b.py --part=all
Reads runs/bearing_round5_plates.json (the plate solves) and runs/bearing_round5_final.json.  Writes runs/round6b_logs/*.json.
Environment: ORBITAL_RING_OUT (default runs).  US spelling in the printed text.
"""
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round5 as r5
import bearing_stageB as sb
import compensation_fields as cf
import final_bearing as fb
import lim_hysteresis as lh
import magnetostatic_2d as ms

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
LOGS = os.path.join(OUT, "round6b_logs")
MU0 = 4e-7 * math.pi
W_HALF = be.W_HALF
K_FLUX = r4.K_FLUX
TAU = 10e-3
N_TURNS = 10                                   # Paul's ruling: 10 series turns
N_FAM, L_FAST = 20, 25.0
E_INS = 5e6                                    # V/m, the insulation's design stress (5 kV/mm, round 5)
TAPE_PITCH = r5.TAPE_PITCH
E_LAUNCH = 4.6e15                              # J, the heaviest launch (5,000 tonnes to 30 km/s, Table 9.5; Paul's figure)
COP = 0.0106                                   # Chapter 14's compressor COP
LAUNCH_RATES = (7, 14)                         # of the heaviest a day
SIDE = ring.CABLE_SIDE

# the candidates of the trade table: (plate width m, thickness m, loop width m, label)
FAST_CAND = [(1.4, 0.20, 0.8, "the book's plate, the book's 0.8 m loop"), (1.4, 0.20, 0.6, "the book's plate, 0.6 m loop"), (1.2, 0.15, 0.6, ""), (1.0, 0.20, 0.6, "out of saturation"),
             (0.8, 0.20, 0.6, "out of saturation, narrow"), (1.0, 0.15, 0.6, "round 5's optimum"), (1.2, 0.10, 0.6, "thin")]
SUP_CAND = [(1.4, 0.20, 0.5, "the book's plate, the book's 0.5 m loops"), (1.4, 0.20, 0.4, "the book's plate, 0.4 m loops"), (1.2, 0.15, 0.4, ""), (1.2, 0.125, 0.4, ""),
            (1.4, 0.10, 0.4, ""), (1.2, 0.10, 0.4, "round 5's optimum"), (1.2, 0.075, 0.4, "thin")]


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))


def load_final():
    fin = json.load(open(os.path.join(OUT, "bearing_round5_final.json"), encoding="utf-8"))
    return fin["results"]["support:opt135"], fin["results"]["fast:opt135"], fin


def cached(kind, Wp, T, Wl, h=0.0025):
    c = r5.load_cache().get(r5.key_of(kind, Wp, T, Wl, h))
    return c if c and c.get("ok") else None


def cable_of(hw):
    """The cable's mass per meter with hardware and its side for a hardware mass (the ring's quadratic, in a subprocess: ring_altitude reads the
    hardware from the environment when it is imported)."""
    env = dict(os.environ, ORBITAL_RING_HW_KGM=f"{hw:.1f}", PYTHONUTF8="1")
    out = subprocess.run([sys.executable, "-c", "import ring_altitude as r; print(r.M_CABLE_M, r.CABLE_SIDE, r.V_REL)"], cwd=os.path.dirname(os.path.abspath(__file__)),
                         env=env, capture_output=True, text=True).stdout.split()
    return float(out[0]), float(out[1]), float(out[2])


# =============================================================================
# S1.1  where the legs are
# =============================================================================

def track_legs(kind, Wl):
    geom = dict(bo.BEST_SUPPORT if kind == "support" else bo.BEST_LIFT); geom["W"] = Wl
    tr = be.Track(geom, 0.01, w_tape=0.012 if kind == "support" else 0.003, I_ref=bo.I_REF_12 if kind == "support" else bo.I_REF_3, NI_ref=300e3)
    return tr


def two_track_legs(d_edge, rs, rl, NI_s=0.0, NI_l=0.0, h=0.01):
    """round 5's two_track with the legs of the plate solve: the support's four 150 x 67 mm bundles at the solve's positions (the center pair
    touching), the load's two 200 x 50 mm bundles, the bundles' tops 10 mm under the casing's face as in the solve; linear plates."""
    g = ms.Grid(9.0, 3.0, h, x0=-3.0, y0=-1.0)
    mu = cf.PLATE_MU; y_top = -bo.WALL_DESIGN
    g.add_iron(-rs["Wp"] / 2, rs["Wp"] / 2, bd.GAP0, bd.GAP0 + rs["T"], mu)
    xl = rs["Wp"] / 2 + d_edge + rl["Wp"] / 2
    g.add_iron(xl - rl["Wp"] / 2, xl + rl["Wp"] / 2, bd.GAP0, bd.GAP0 + rl["T"], mu)
    ts = track_legs("support", rs["Wl"]); tl = track_legs("fast", rl["Wl"])
    legs_s = [(x, s * NI_s, ts.wb, ts.hb) for x, s in ts.legs]
    legs_l = [(xl + x, s * NI_l, tl.wb, tl.hb) for x, s in tl.legs]
    for x, a, wb, hb in legs_s + legs_l:
        if a:
            g.add_current(x - wb / 2, x + wb / 2, y_top - hb, y_top, a)
    g.solve()
    def mean_A(x, wb, hb):
        return float(g.A[g.box(x - wb / 2, x + wb / 2, y_top - hb, y_top)].mean())
    A = [mean_A(x, wb, hb) for x, a, wb, hb in legs_s]
    sg = [s for x, s in ts.legs]
    flux = sum(s * a for s, a in zip(sg, A))                             # the figure-eight's flux per meter in its own sense: sum of sign x A over its four legs
    return g, flux, legs_s, legs_l, xl


def coupling_legs(rs, rl, d_edge):
    g, phi_s, legs_s, _, xl = two_track_legs(d_edge, rs, rl, NI_s=1e3)
    L_s = phi_s / 1e3
    g, phi_m, _, _, _ = two_track_legs(d_edge, rs, rl, NI_l=1e3)
    M = phi_m / 1e3
    NI_l_full = rl["NI_req"] * math.sqrt(fb.F_LOAD_TRACK / rl["F_req"]); NI_s_op = rs["NI_req"] * math.sqrt(W_HALF / rs["F_req"])
    dNI = -M / L_s * NI_l_full
    return dict(d_edge=d_edge, xl=xl, L_s=L_s, M=M, k=M / L_s, dF=2 * dNI / NI_s_op, NI_l_full=NI_l_full, NI_s_op=NI_s_op)


def part_legs():
    rs, rl, fin = load_final()
    hr("S1.1  WHERE THE LEGS ARE: THE PLATE SOLVE, THE COUPLING RUN, THE ROLL RUN")
    ts = track_legs("support", rs["Wl"]); tl = track_legs("fast", rl["Wl"])
    print(f"  The plate solve (bearing_round5.solve_point, bearing_electrical.Track on bearing_optimize's geometry), the support track on its {rs['Wp']:.1f} m plate:")
    print(f"  n = 2 loops of W = {rs['Wl']:.2f} m between leg centers, bundles '{ts.shape}' ({ts.wb*1e3:.0f} x {ts.hb*1e3:.1f} mm), the two center bundles touching (gc = 0, their centers {ts.wb*1e3:.0f} mm apart):")
    for x, s in ts.legs:
        print(f"    leg at {x*1e3:+7.1f} mm from the plate's center, sign {s:+d}; the bundle from {(x-ts.wb/2)*1e3:+7.1f} to {(x+ts.wb/2)*1e3:+7.1f} mm")
    span = ts.legs[-1][0] - ts.legs[0][0] + ts.wb
    print(f"  Four bundles, symmetric about the plate's center, two equal loops of {rs['Wl']*1e3:.0f} mm between leg centers; the windings span {span*1e3:.0f} mm under the {rs['Wp']*1e3:.0f} mm plate,")
    print(f"  {(rs['Wp']-span)/2*1e3:.0f} mm of plate past each outer bundle's edge.  The fast track on its {rl['Wp']:.1f} m plate: one loop of {rl['Wl']:.2f} m, bundles '{tl.shape}':")
    for x, s in tl.legs:
        print(f"    leg at {x*1e3:+7.1f} mm, sign {s:+d}; the bundle from {(x-tl.wb/2)*1e3:+7.1f} to {(x+tl.wb/2)*1e3:+7.1f} mm")
    span_l = tl.legs[-1][0] - tl.legs[0][0] + tl.wb
    print(f"  spanning {span_l*1e3:.0f} mm under the {rl['Wp']*1e3:.0f} mm plate ({(rl['Wp']-span_l)/2*1e3:.0f} mm past each bundle's edge).")
    print(f"  The coupling run (round 5's two_track, Stage B's model, used in rounds 5 and 6a): the support as THREE uniform 100 x 100 mm bundles at {-rs['Wl']*1e3:+.0f}, 0 and {rs['Wl']*1e3:+.0f} mm")
    print(f"  (the center pair merged into one bundle with twice the ampere-turns), the load as two 100 x 100 mm bundles at +-{rl['Wl']/2*1e3:.0f} mm from its plate's center.  Two equal loops,")
    print(f"  symmetric, but the outer legs {abs(ts.legs[0][0])*1e3 - rs['Wl']*1e3:.0f} mm further in than the solve's and the bundles not the solve's shape.  The roll run (compensation_tracks.roll_model): no legs at all;")
    print(f"  each track is a force at its plate's center (support +-{SIDE/2 - rl['Wp'] - 0.2 - rs['Wp']/2:.3f} m at 0.2 m spacing, lift and load +-{SIDE/2 - rl['Wp']/2:.3f} m).  A symmetric track's force acts at its center, so the")
    print(f"  roll run does not depend on the legs.  Drawing 1b's bundles at -400, 0, +150, +400 mm were none of these: a drawing error.")
    print()
    print(f"  The coupling rerun on the solve's legs (four {ts.wb*1e3:.0f} x {ts.hb*1e3:.0f} mm support bundles, two {tl.wb*1e3:.0f} x {tl.hb*1e3:.0f} mm load bundles, tops 10 mm under the casing's face; linear plates, 10 mm cells):")
    print(f"  {'d m':>6} | {'three-bundle model (rounds 5, 6a): M/L %':>40} {'dF/F %':>7} | {'the solve s legs: M/L %':>24} {'dF/F %':>7} {'kN/m per track':>14} {'kN/m pair':>9} {'as gap offset mm':>16}")
    rows = []
    for d in (0.2, 0.6, 1.0):
        old = r5.coupling(rs, rl, d); new = coupling_legs(rs, rl, d)
        dF = abs(new["dF"]); per = dF * W_HALF
        rows.append(dict(d=d, k_old=old["k"], dF_old=abs(old["dF"]), k=new["k"], dF=dF, per_track=per, pair=2 * per, eq_mm=2 * per / (2 * K_FLUX) * 1e3, L_s=new["L_s"], L_s_old=old["L_s"], M=new["M"], M_old=old["M"]))
        print(f"  {d:6.2f} | {old['k']*100:40.2f} {abs(old['dF'])*100:7.2f} | {new['k']*100:24.2f} {dF*100:7.2f} {per/1e3:14.2f} {2*per/1e3:9.2f} {2*per/(2*K_FLUX)*1e3:16.1f}", flush=True)
    r0 = rows[0]
    print(f"  At 0.2 m: the support's own flux per ampere-turn {r0['L_s_old']*1e6:.3f} -> {r0['L_s']*1e6:.3f} uWb/m per At, the load's flux through it {r0['M_old']*1e6:.4f} -> {r0['M']*1e6:.4f}.")
    save("round6b_legs.json", dict(support_legs=[(x, s) for x, s in ts.legs], fast_legs=[(x, s) for x, s in tl.legs], wb_s=ts.wb, hb_s=ts.hb, wb_l=tl.wb, hb_l=tl.hb, rows=rows))
    return rows


# =============================================================================
# S1.2  the dump
# =============================================================================

def support_coil(rs):
    """The optimum plate's support coil as whole turns: 81 (41 + 40), its currents, and its inductances two ways: round 4's 55-turn values (the
    book's plate) scaled as the turns squared, and from the optimum plate's own stored energy."""
    turns = 81
    NI0 = rs["NI_req"] * math.sqrt(W_HALF / rs["F_req"])                 # per leg at the normal force
    n_leg = rs["turns_leg"] / (r5.N_SUP * 10)                            # the solve's turns of one coil per leg (40.7 of ten layers; 41 and 40 as whole turns), as rounds 5 and 6a
    I0 = NI0 / (r5.N_SUP * n_leg); I_cap = rs["NI_req"] / (r5.N_SUP * n_leg)
    Lc_scaled = 21.3 * (turns / 55)**2; Ll_scaled = 0.059 * (turns / 55)**2
    W0 = rs["W_mag"] * (W_HALF / rs["F_req"])                            # J per m of track at the normal force (the energy as the force, both as NI^2)
    Lc_energy = 2 * W0 * r5.L_SUP / (r5.N_SUP * I0**2)                   # 10 coils x 1/2 L I^2 = W0 x 1000 m
    return dict(turns=turns, I0=I0, I_cap=I_cap, NI0=NI0, Lc_scaled=Lc_scaled, Ll_scaled=Ll_scaled, Lc_energy=Lc_energy, Ll_energy=Lc_energy * 0.059 / 21.3, W0=W0)


def part_dump():
    rs, rl, fin = load_final()
    hr("S1.2  THE DUMP OF A SUPPORT COIL WITH D2 AND THE STRIP ACROSS THE COIL")
    sc = support_coil(rs); I0, Ic = sc["I0"], sc["I_cap"]
    print(f"  The coil: {sc['turns']} turns, {I0:.0f} A normal, {Ic:,.0f} A at the cap.  Inductance with the ten coils of the bundle together: {sc['Lc_scaled']:.1f} H (round 4's 21.3 H at 55 turns on the book's")
    print(f"  plate, scaled as the turns squared; round 6a's 47 H) or {sc['Lc_energy']:.1f} H (from the optimum plate's own stored energy, {sc['W0']/1e3:.1f} kJ per meter of track at the normal force:")
    print(f"  ten coils x 1/2 L I^2 = W x 1,000 m).  One coil alone against nine persistent neighbors (its leakage, 0.28 percent of the common mode): {sc['Ll_scaled']*1e3:.0f} or {sc['Ll_energy']*1e3:.0f} mH.")
    print(f"  The two differ because the optimum plate is thinner and narrower than the book's and its loops are 0.4 m, not 0.5 m; the energy-based value is the optimum plate's own.")
    print()
    print(f"  As drawn (6a) and as written (round 4, 2.1; Stage B 2.4): the dump branch (a switch into a resistor) across the coil, D2 across the coil, the strip across D2.")
    print(f"  That cannot dump: with the strip heated open the coil's current free-wheels through D2 at 1.0 V, and a resistor R in a parallel branch takes 1.0 V / R:")
    for R in (2.2, 0.85, 20.0):
        print(f"    R = {R:5.2f} ohm: {1.0/R:.2f} A of the coil's {I0:.0f} A; the coil's voltage stays at 1 V and its current decays at 1 V / L_leak = {1.0/sc['Ll_scaled']:.1f} A/s (one alone).")
    print(f"  The fault is in the design as written, not only in the drawing.  The correction: the resistor goes in SERIES with D2, and the dump switch across the resistor,")
    print(f"  normally closed.  Across the coil's terminals then sit (i) the strip, (ii) D2 in series with [the resistor, bypassed by the closed switch], (iii) the supply through")
    print(f"  the lead contacts.  Charging: the strip open, the supply's voltage across the coil, D2 reverse-biased and blocking it; nothing in branch (ii) whatever the switch does.")
    print(f"  Free-wheel (the supply removed, the strip still open): through D2 and the closed switch, at D2's 1.0 V plus the switch's on-drop.  Persistent: the strip carries")
    print(f"  the current at zero volts.  Dump: the strip heated open and the switch opened; the current moves from the strip into D2 and the resistor, and the coil sees I R.")
    print(f"  {'case':>44} {'I A':>6} {'R ohm':>6} {'V at the start':>14} {'L_leak / R ms':>13} {'energy kJ (leakage)':>19} {'in the open strip %':>19}")
    out = []
    for lab, I, V in (("normal current, 1 kV", I0, 1000.0), ("the cap, 1 kV", Ic, 1000.0), ("normal current, R as for the cap", I0, None), ("normal current, the charging voltage", I0, "charge")):
        for Ll, tag in ((sc["Ll_scaled"], "scaled"), (sc["Ll_energy"], "energy")):
            V_ch = (sc["Lc_scaled"] if tag == "scaled" else sc["Lc_energy"]) * 0.01 * I0
            R_strip = V_ch / (0.01 * I0)
            if V is None: R = 1000.0 / Ic; Vs = I * R
            elif V == "charge": R = V_ch / I; Vs = V_ch
            else: R = V / I; Vs = V
            E = 0.5 * Ll * I**2
            out.append(dict(case=lab, L=tag, I=I, R=R, V=Vs, tau=Ll / R, E=E, R_strip=R_strip, strip_share=R / (R + R_strip)))
            print(f"  {lab + ' (' + tag + ' L)':>44} {I:6.0f} {R:6.2f} {Vs:12.0f} V {Ll/R*1e3:13.0f} {E/1e3:19.1f} {R/(R+R_strip)*100:19.1f}")
    V_ch = sc["Lc_scaled"] * 0.01 * I0
    print(f"  The charging voltage for 1 percent per second, ten together: {V_ch:.0f} V (scaled L) or {sc['Lc_energy']*0.01*I0:.0f} V (energy L); the strip by Stage B's rule (1 percent of the coil current")
    print(f"  through it while charging) {V_ch/(0.01*I0):.0f} ohm.  With the strip open across the coil during a 1 kV dump it carries 1 kV / R_strip = {1000/(V_ch/(0.01*I0)):.0f} A and {1e6/(V_ch/(0.01*I0))/1e3:.0f} kW")
    print(f"  for the dump's tens of milliseconds ({1e6/(V_ch/(0.01*I0))*sc['Ll_scaled']/2.2/2/1e3:.2f} kJ over a decay at {sc['Ll_scaled']/2.2*1e3:.0f} ms); the strip was sized for {V_ch**2/(V_ch/(0.01*I0))/1e3:.2f} kW while charging.  A dump at the charging voltage")
    print(f"  instead ({V_ch:.0f} V, R = {V_ch/I0:.2f} ohm) takes {sc['Ll_scaled']/(V_ch/I0)*1e3:.0f} ms, the 0.3 s of the 551-turn coil's 1 kV dump, and keeps the strip at its charging load.")
    print(f"  The switch: closed it carries the coil's current only while free-wheeling (seconds, in procedure (a) of round 5 Section 8); open it blocks the dump voltage.")
    print(f"  Its on-drop adds to D2's 1.0 V in the free-wheel: at 2 V in all the ten coils together lose {2.0*60/sc['Lc_scaled']:.1f} A in a minute ({2.0*60/sc['Lc_scaled']/I0*100:.2f} percent), a lone coil {2.0/sc['Ll_scaled']:.0f} A/s.")
    save("round6b_dump.json", dict(coil=sc, rows=out))
    return out


# =============================================================================
# S1.3  the film between turns in a fast coil's cell
# =============================================================================

def cell_orders(n_x, n_y, N):
    """Orders of N turns in an n_x by n_y grid of stacks: name -> list of (col, row) for turn 1..N."""
    rowmajor = [(k % n_x, k // n_x) for k in range(N)]
    snake = [((k % n_x) if (k // n_x) % 2 == 0 else n_x - 1 - (k % n_x), k // n_x) for k in range(N)]
    colmajor = [(k // n_y, k % n_y) for k in range(N)]
    colsnake = [(k // n_y, (k % n_y) if (k // n_y) % 2 == 0 else n_y - 1 - (k % n_y)) for k in range(N)]
    return {"rows, each left to right": rowmajor, "rows, back and forth": snake, "columns, each bottom to top": colmajor, "columns, up and down": colsnake}


def film_layout(order, n_x, n_y, w_s, h_s, V_turn, E=E_INS, film_min=0.05e-3):
    """Films sized for the turns they separate: every film between two columns as thick as the largest voltage across it needs, the same for
    the films between rows; then the cell's size.  Returns the films' voltages and the cell."""
    pos = {p: k + 1 for k, p in enumerate(order)}
    v_col = [0.0] * (n_x - 1); v_row = [0.0] * (n_y - 1); pairs = []; v_diag = 0.0
    for (c, r), k in pos.items():
        if (c + 1, r) in pos:
            dv = abs(pos[(c + 1, r)] - k) * V_turn; v_col[c] = max(v_col[c], dv); pairs.append(("side", k, pos[(c + 1, r)], dv))
        if (c, r + 1) in pos:
            dv = abs(pos[(c, r + 1)] - k) * V_turn; v_row[r] = max(v_row[r], dv); pairs.append(("above", k, pos[(c, r + 1)], dv))
        for dc in (-1, 1):
            if (c + dc, r + 1) in pos:
                dv = abs(pos[(c + dc, r + 1)] - k) * V_turn; v_diag = max(v_diag, dv); pairs.append(("diagonal", k, pos[(c + dc, r + 1)], dv))
    t_col = [max(v / E, film_min) for v in v_col]; t_row = [max(v / E, film_min) for v in v_row]
    w_c = n_x * w_s + sum(t_col); h_c = n_y * h_s + sum(t_row)
    # the diagonal: two stacks that meet corner to corner through the crossing of a column film and a row film
    worst_diag = 0.0
    for kind, a, b, dv in pairs:
        if kind == "diagonal":
            (ca, ra), (cb, rb) = order[a - 1], order[b - 1]
            dist = math.hypot(t_col[min(ca, cb)], t_row[min(ra, rb)])
            worst_diag = max(worst_diag, dv / dist)
    return dict(t_col=t_col, t_row=t_row, v_col=v_col, v_row=v_row, w_c=w_c, h_c=h_c, pairs=pairs, E_diag=worst_diag)


def part_film():
    rs, rl, fin = load_final()
    hr("S1.3  THE FILM BETWEEN TURNS IN A FAST COIL'S CELL (10 turns, 5 stacks across and 2 up)")
    fbk = r5.fast_basics(rl, N_TURNS)
    layers = math.ceil(fbk["layers"]); V_turn = fbk["V_turn"]; V = fbk["V"]
    w_s, h_s = 0.003, layers * TAPE_PITCH
    t_fam = V / E_INS
    w_cell = (0.200 - 11 * t_fam) / 10
    print(f"  The coil on round 5's optimum plate: {N_TURNS} turns of {layers} layers, {fbk['I']:,.0f} A, {V_turn:,.0f} V per turn, {V/1e3:.1f} kV at the terminals at the 10 ms rise.  Stack 3.00 x {h_s*1e3:.2f} mm")
    print(f"  (faces horizontal).  Walls between cells and to ground for the terminal voltage at {E_INS/1e6:.0f} kV/mm: {t_fam*1e3:.2f} mm; ten cells across 200 mm leave {w_cell*1e3:.2f} mm inside a cell.")
    print(f"  A film is sized for the largest voltage between the turns it separates, at {E_INS/1e6:.0f} kV/mm.  Four orders of the ten turns in the 5 by 2 grid:")
    out = {}
    for name, order in cell_orders(5, 2, N_TURNS).items():
        fl = film_layout(order, 5, 2, w_s, h_s, V_turn)
        h_bundle = 2 * fl["h_c"] + 3 * t_fam
        fits_w = fl["w_c"] <= w_cell + 1e-9
        grid = [[0] * 5 for _ in range(2)]
        for k, (c, r) in enumerate(order):
            grid[r][c] = k + 1
        out[name] = dict(order=order, **{k: v for k, v in fl.items() if k != "pairs"}, h_bundle=h_bundle, fits_w=fits_w, grid=grid)
        print(f"\n  {name}:   top row {grid[1]}   bottom row {grid[0]}")
        print(f"    films between columns: {', '.join(f'{v/1e3:.1f} kV -> {t*1e3:.2f} mm' for v, t in zip(fl['v_col'], fl['t_col']))}")
        print(f"    film between the rows: {', '.join(f'{v/1e3:.1f} kV -> {t*1e3:.2f} mm' for v, t in zip(fl['v_row'], fl['t_row']))}")
        ab = [p for p in fl["pairs"] if p[0] == "above"]
        print(f"    the pairs facing across the rows: " + ", ".join(f"{a}/{b} {dv/1e3:.1f} kV" for _, a, b, dv in sorted(ab, key=lambda p: p[1])))
        print(f"    cell {fl['w_c']*1e3:.2f} x {fl['h_c']*1e3:.2f} mm ({'fits' if fits_w else 'DOES NOT FIT'} the {w_cell*1e3:.2f} mm inside width); bundle height 2 cells + 3 walls = {h_bundle*1e3:.1f} mm of the 50;"
              f" corner to corner {fl['E_diag']/1e6:.1f} kV/mm")
    best = out["rows, each left to right"]
    # the row film thickened so that the corner-to-corner path is at the design stress too
    order = best["order"]; fl = film_layout(order, 5, 2, w_s, h_s, V_turn)
    v_diag = max(dv for kind, a, b, dv in fl["pairs"] if kind == "diagonal")
    t_need = math.sqrt(max((v_diag / E_INS)**2 - fl["t_col"][0]**2, 0.0))
    h_c2 = 2 * h_s + max(t_need, fl["t_row"][0])
    print(f"\n  The order carried: rows, each left to right (turns 1 to 5 in the lower row, 6 to 10 above them).  Eight films of {fl['t_col'][0]*1e3:.2f} mm between neighbors in a row ({V_turn/1e3:.2f} kV),")
    print(f"  one film of {fl['t_row'][0]*1e3:.2f} mm between the rows (five turns, {5*V_turn/1e3:.1f} kV).  Cell {fl['w_c']*1e3:.2f} x {fl['h_c']*1e3:.2f} mm, bundle {best['h_bundle']*1e3:.1f} mm of the 50 (round 6a's 43.1 had the one-turn film")
    print(f"  between the rows).  The stacks that meet corner to corner across the film crossing are six turns apart ({v_diag/1e3:.1f} kV over {math.hypot(fl['t_col'][0], fl['t_row'][0])*1e3:.2f} mm, {fl['E_diag']/1e6:.1f} kV/mm);")
    print(f"  holding that path to {E_INS/1e6:.0f} kV/mm too takes a row film of {t_need*1e3:.2f} mm: cell height {h_c2*1e3:.2f} mm, bundle {(2*h_c2 + 3*t_fam)*1e3:.1f} mm of the 50.")
    out["carried"] = dict(layers=layers, V_turn=V_turn, V=V, t_fam=t_fam, w_cell=w_cell, h_s=h_s, t_row_diag=t_need, h_bundle_diag=2 * h_c2 + 3 * t_fam)
    save("round6b_film.json", out)
    return out


# =============================================================================
# the tape's pulse loss of a launch for a plate, a bundle shape and a stack orientation
# =============================================================================

def shape_name(w, h):
    name = f"r6b {w*1e3:.0f}x{h*1e3:.0f}"
    bo.SHAPES[name] = dict(kind="rect", w=w, h=h)
    return name


def leg_fields(Wp, T, Wl, NI, h, shape=None):
    """The fast track at NI ampere-turns per leg on the plate Wp x T: the field at every model tape of both legs."""
    geom = dict(bo.BEST_LIFT); geom["W"] = Wl
    if shape:
        geom["shape"] = shape
    with r4.plate_patch(Wp, T):
        tr = be.Track(geom, h, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=NI)
        g, s = tr.solve(tr.NI_design(NI))
        fields = [bo.field_at_points(g, xs, ys) for xs, ys in tr.lattices]
    return fields, s, tr


def face_components(Bx, By, face):
    """The field through the face (Bn) and along it (Bt) for stacks lying faces horizontal (h), faces vertical (v), each stack the better of the
    two (best), or every stack turned along its local field (aligned, the floor)."""
    ax, ay = np.abs(Bx), np.abs(By)
    if face == "h": return ay, ax
    if face == "v": return ax, ay
    if face == "best": return np.minimum(ax, ay), np.maximum(ax, ay)
    return np.zeros_like(ax), np.hypot(ax, ay)


SPEEDS = [1e3, 2e3, 3e3, 5e3, 7e3, 10e3, 15e3, 20e3, 25e3, 30e3]
_LI = {}


def launch_profile(W_mag, F_req):
    key = (round(W_mag), round(F_req))
    if key not in _LI:
        _LI[key] = r5.launch_integrals(W_mag, F_req)
    return _LI[key]


def tape_loss_launch(fields, face, tapes_leg, layers, li):
    """Per launch, the pair of pulsing tracks: the Brandt strip loss of the field through the faces and the Bean slab loss of the field along
    them (the soldered stack one slab), round 5 Section 5's model with the orientation free.  J."""
    w = 0.003; d = layers * TAPE_PITCH
    n_model = sum(len(f[0]) for f in fields) / len(fields); scale = tapes_leg / n_model
    per_u = {}
    for u in SPEEDS:
        F, kind = r4.pull_track_of_u(u); k = math.sqrt(F / li["F_req"]) if F > 0 else 0.0
        qp = qs = 0.0
        if k > 0:
            for (Bx, By) in fields:
                Bn, Bt = face_components(Bx * k, By * k, face)
                Bm = np.hypot(Bx, By) * k; th = np.degrees(np.arctan2(Bt, Bn)); Ic = lh.ic_tape(Bm, th, w)
                qp += float(lh.q_strip_perp(Ic, w, Bn / 2.0).sum()) * scale
                Bp = MU0 * (Ic / (w * TAPE_PITCH)) * d / 2.0
                qs += float(np.array([r5.q_bean_slab(b / 2.0, bp) for b, bp in zip(Bt, Bp)]).sum()) * w * TAPE_PITCH * scale
        per_u[u] = (2 * qp, 2 * qs)
    Ep = Es = 0.0; t_p, v_p = li["t_p"], li["v_p"]
    for i in range(1, len(t_p)):
        vv = float(v_p[i]); dtt = float(t_p[i] - t_p[i - 1])
        if vv < SPEEDS[0]:
            continue
        Ep += float(np.interp(vv, SPEEDS, [per_u[u][0] for u in SPEEDS])) * vv * dtt
        Es += float(np.interp(vv, SPEEDS, [per_u[u][1] for u in SPEEDS])) * vv * dtt
    return Ep, Es


def worst_ic(fields, face, theva):
    """The worst tape's limit for the orientation (the THEVA and Robinson limit of the plate solves), with its field and angle from the normal."""
    best = None
    for (Bx, By) in fields:
        Bn, Bt = face_components(Bx, By, face); Bm = np.hypot(Bx, By); th = np.degrees(np.arctan2(Bt, Bn))
        ic = bo.ic_angle(theva, Bm, th, 0.003); k = int(ic.argmin())
        if best is None or ic[k] < best[0]:
            best = (float(ic[k]), float(Bm[k]), float(th[k]))
    return best


def _plate_loss_worker(args):
    Wp, T, Wl, lab = args
    res = cached("fast", Wp, T, Wl)
    if res is None:
        return dict(Wp=Wp, T=T, Wl=Wl, ok=False)
    t0 = time.time()
    fields, s, tr = leg_fields(Wp, T, Wl, res["NI_req"], 0.005)
    li = dict(launch_profile(res["W_mag"], res["F_req"]), F_req=res["F_req"])
    fbk = r5.fast_basics(res, N_TURNS)
    Ep, Es = tape_loss_launch(fields, "h", res["turns_leg"], fbk["layers"], li)
    ph = r5.plate_heat_launch(res)
    return dict(Wp=Wp, T=T, Wl=Wl, ok=True, Ep=Ep, Es=Es, plate_heat=ph["E_launch"], Bp30=ph["rows"][-1]["Bp"], F5=s["F"], seconds=time.time() - t0, I=fbk["I"], V=fbk["V"], layers=fbk["layers"],
                E_cycle=li["E_cycle"], lead_s=li["lead_s"], x_run=li["x_run"])


def part_plates(jobs=4):
    import multiprocessing as mp
    rs_opt, rl_opt, fin = load_final()
    hr("S2.1  A TRADE TABLE OF PLATES (2.5 mm cells, the worst tape at half its limit at the requirement, the convergence step; 10 turns, faces horizontal)")
    print(f"  The launch: the heaviest (5,000 tonnes to 30 km/s), {E_LAUNCH/1e15:.1f} PJ (Paul's figure, Table 9.5).  Cold heat is charged at the compressor's COP of {COP} (Chapter 14): a joule in the")
    print(f"  cold costs {1/COP:.0f} J at the compressor.  The tape's pulse loss is round 5 Section 5's model (the strip loss of the field through the faces plus the slab loss of the")
    print(f"  field along them, the track resolved at 5 mm at the 2.5 mm solve's ampere-turns); the plate heat round 5's (hysteresis and eddy, the solved peak field of each pass).")
    cands = [c for c in FAST_CAND if cached("fast", *c[:3])]
    missing = [c for c in FAST_CAND if not cached("fast", *c[:3])]
    if missing:
        print(f"  NOT YET SOLVED at 2.5 mm (left out): " + "; ".join(f"{c[0]:.1f} x {c[1]*1e3:.0f}, loop {c[2]:.1f}" for c in missing))
    with mp.Pool(jobs) as pool:
        losses = pool.map(_plate_loss_worker, cands)
    rows = []
    print(f"\n  The fast (lift and load) plate; the hardware and the cable with round 5's optimum support plates ({rs_opt['plate_m']:,.0f} kg/m each):")
    print(f"  {'plate W x T, loop':>24} {'plate kg/m':>10} {'hardware kg/m':>13} {'cable + hw kg/m':>15} {'sat %':>6} {'peak T':>6} {'kA':>4} {'tape kg/m':>9} {'margin':>6} {'tape loss GJ':>12} {'(strip + slab)':>14} {'plate heat GJ':>13} {'cold, % of launch':>17} {'plate, %':>8}  note")
    for (Wp, T, Wl, lab), L in zip(cands, losses):
        res = cached("fast", Wp, T, Wl)
        hw = fb.M_LIM_PLATES + 2 * rs_opt["plate_m"] + 4 * res["plate_m"]
        m_c, side, v_rel = cable_of(hw)
        Et = L["Ep"] + L["Es"]
        row = dict(Wp=Wp, T=T, Wl=Wl, label=lab, plate_m=res["plate_m"], hw=hw, cable=m_c, side=side, sat=res["sat"], B_plate=res["B_plate"], NI=res["NI_req"], tape=res["tape"], margin=res["margin"],
                   ic=res["ic_req"], B=res["B"], th=res["th"], E_tape=Et, Ep=L["Ep"], Es=L["Es"], E_plate=L["plate_heat"], Bp30=L["Bp30"], cold_share=Et / COP / E_LAUNCH, plate_share=L["plate_heat"] / E_LAUNCH,
                   I=L["I"], V=L["V"], layers=L["layers"], E_cycle=L["E_cycle"], lead_s=L["lead_s"], x_run=L["x_run"], turns_leg=res["turns_leg"])
        rows.append(row)
        print(f"  {f'{Wp:.1f} m x {T*1e3:.0f} mm, {Wl:.1f} m':>24} {res['plate_m']:10,.0f} {hw:13,.0f} {m_c:15,.0f} {res['sat']*100:6.1f} {res['B_plate']:6.2f} {res['NI_req']/1e3:4.0f} {res['tape']:9.1f} {res['margin']:6.2f} {Et/1e9:12.0f} "
              f"{f'{L['Ep']/1e9:.0f} + {L['Es']/1e9:.0f}':>14} {L['plate_heat']/1e9:13.0f} {Et/COP/E_LAUNCH*100:17.2f} {L['plate_heat']/E_LAUNCH*100:8.3f}  {lab}", flush=True)
    print(f"\n  The support plate (a DC track: no pulse loss in its tape and no cyclic heat in its plate; it is chosen on mass, tape, margin and saturation); the hardware and")
    print(f"  the cable with round 5's optimum fast plates ({rl_opt['plate_m']:,.0f} kg/m each):")
    print(f"  {'plate W x T, loops':>24} {'plate kg/m':>10} {'hardware kg/m':>13} {'cable + hw kg/m':>15} {'sat % at the cap':>16} {'peak T':>6} {'kA':>4} {'tape kg/m':>9} {'margin':>6} {'worst tape':>22}  note")
    srows = []
    for Wp, T, Wl, lab in SUP_CAND:
        res = cached("support", Wp, T, Wl)
        if res is None:
            print(f"  {f'{Wp:.1f} m x {T*1e3:.0f} mm, {Wl:.1f} m':>24}  NOT YET SOLVED at 2.5 mm"); continue
        hw = fb.M_LIM_PLATES + 2 * res["plate_m"] + 4 * rl_opt["plate_m"]
        m_c, side, v_rel = cable_of(hw)
        srows.append(dict(Wp=Wp, T=T, Wl=Wl, label=lab, plate_m=res["plate_m"], hw=hw, cable=m_c, sat=res["sat"], B_plate=res["B_plate"], NI=res["NI_req"], tape=res["tape"], margin=res["margin"], ic=res["ic_req"], B=res["B"], th=res["th"]))
        print(f"  {f'{Wp:.1f} m x {T*1e3:.0f} mm, {Wl:.1f} m':>24} {res['plate_m']:10,.0f} {hw:13,.0f} {m_c:15,.0f} {res['sat']*100:16.1f} {res['B_plate']:6.2f} {res['NI_req']/1e3:4.0f} {res['tape']:9.1f} {res['margin']:6.2f} "
              f"{f'{res['B']:.2f} T at {res['th']:.0f} deg, {res['ic_req']:.0f} A':>22}  {lab}")
    # every heat of the bearing per launch, on round 5's optimum and on each candidate: the leads' pulse heat, the converter's loss
    print(f"\n  Every heat of the bearing per launch (the heaviest), as a share of its {E_LAUNCH/1e15:.1f} PJ, and as a rate at {LAUNCH_RATES[0]} and {LAUNCH_RATES[1]} of the heaviest a day around the {ring.L_RING/1e3:,.0f} km ring:")
    heats = {}
    for row in rows:
        A_p = r5.lead_rows(cached("fast", row["Wp"], row["T"], row["Wl"]))[0][N_TURNS]["A_m"]; R_p = be.RHO_CU_77 * 1.0 / A_p      # round 5's lead
        E_leads = 2 * row["I"]**2 * R_p * row["lead_s"]
        E_conv = 0.03 * row["E_cycle"]
        heats[f"{row['Wp']:.1f}x{row['T']*1e3:.0f}/{row['Wl']:.1f}"] = dict(E_tape=row["E_tape"], E_leads=E_leads, E_plate=row["E_plate"], E_conv=E_conv, E_cycle=row["E_cycle"], A_lead=A_p)
        row.update(E_leads=E_leads, E_conv=E_conv)
    print(f"  {'plate':>24} | {'item':>44} {'GJ per launch':>13} {'% of the launch':>15} {'at the compressor, %':>20} {'W/m at 7 a day':>14} {'at 14':>7} {'MW at 14':>9}")
    for row in rows:
        items = [("the tape's pulse loss (cold)", row["E_tape"], True), ("the leads' Joule heat in the pulses (cold)", row["E_leads"], True), ("the plates' hysteresis and eddy heat (warm, on the cable)", row["E_plate"], False),
                 ("the converters at 3 percent of the energy cycled (warm) [assumed]", row["E_conv"], False)]
        for lab, E, cold in items:
            print(f"  {f'{row['Wp']:.1f} x {row['T']*1e3:.0f}, {row['Wl']:.1f}':>24} | {lab:>44.44} {E/1e9:13.1f} {E/E_LAUNCH*100:15.4f} {(f'{E/COP/E_LAUNCH*100:.3f}' if cold else '-'):>20} "
                  f"{E*LAUNCH_RATES[0]/86400/ring.L_RING:14.3f} {E*LAUNCH_RATES[1]/86400/ring.L_RING:7.3f} {E*LAUNCH_RATES[1]/86400/1e6:9.1f}")
        cold = row["E_tape"] + row["E_leads"]
        print(f"  {'':>24} | {'cold in all':>44} {cold/1e9:13.1f} {cold/E_LAUNCH*100:15.4f} {cold/COP/E_LAUNCH*100:20.3f} {cold*LAUNCH_RATES[0]/86400/ring.L_RING:14.3f} {cold*LAUNCH_RATES[1]/86400/ring.L_RING:7.3f} {cold*LAUNCH_RATES[1]/86400/1e6:9.1f}")
    print(f"  For scale: the standing cold load is {r5.Q_STANDING_W_PER_M:.0f} W per meter of ring (Chapter 14), the fast coils' leads conduct 12.2 W/m idle with all 20 families connected (round 5).")
    save("round6b_plates.json", dict(fast=rows, support=srows, heats=heats))
    return rows, srows


# =============================================================================
# S2.2  one attempt at the bundle
# =============================================================================

def pack(bundle_w, face, layers, V, V_turn, N=N_TURNS, n_fam=N_FAM):
    """Stack-by-stack packing of n_fam cells of N stacks in a bundle of width bundle_w: for every cell grid (cells across x cells up) and every
    arrangement of the stacks in a cell (columns x rows, row order or column order, the films sized for the turns they separate), the bundle
    height needed.  Returns the lowest."""
    t_s = layers * TAPE_PITCH
    w_s, h_s = (0.003, t_s) if face == "h" else (t_s, 0.003)
    t_fam = V / E_INS; t1 = max(V_turn / E_INS, 0.05e-3)
    best = None
    for n_cx in (1, 2, 4, 5, 10, 20):
        n_cy = n_fam // n_cx
        w_c = (bundle_w - (n_cx + 1) * t_fam) / n_cx
        if w_c <= 0:
            continue
        for n_x in range(1, N + 1):
            n_y = math.ceil(N / n_x)
            for order in ("row", "col"):
                t_col = t1 if order == "row" else n_y * t1               # row order: neighbors in a row one turn apart, rows n_x turns apart
                t_row = n_x * t1 if order == "row" else t1
                w_need = n_x * w_s + (n_x - 1) * t_col
                if w_need > w_c + 1e-9:
                    continue
                h_c = n_y * h_s + (n_y - 1) * t_row
                h_need = n_cy * h_c + (n_cy + 1) * t_fam
                if best is None or h_need < best["h_need"]:
                    best = dict(n_cx=n_cx, n_cy=n_cy, n_x=n_x, n_y=n_y, order=order, w_c=w_c, w_need=w_need, h_c=h_c, h_need=h_need, t_fam=t_fam, t_col=t_col, t_row=t_row, w_s=w_s, h_s=h_s)
    return best


def _shape_worker(args):
    Wp, T, Wl, w, hgt, h, tag = args
    t0 = time.time()
    theva, theva77 = be.load_theva()
    name = shape_name(w, hgt)
    r_sup, r_lift = be.load_records(); F_req = r_lift["F_req"]
    geom = dict(bo.BEST_LIFT); geom["W"] = Wl; geom["shape"] = name
    base = cached("fast", Wp, T, Wl) or cached("fast", Wp, T, Wl, 0.005)
    NI1 = base["NI_req"] if base else r_lift["NI_req"]
    with r4.plate_patch(Wp, T):
        tr = be.Track(geom, h, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=NI1)
        def solve(NI):
            g, s = tr.solve(tr.NI_design(NI))
            return s, [bo.field_at_points(g, xs, ys) for xs, ys in tr.lattices]
        s1, f1 = solve(NI1)
        NIa, Fa = NI1, s1["F"]; NIb = NI1 * math.sqrt(F_req / s1["F"]); sb_, fb_ = solve(NIb); Fb = sb_["F"]; n_it = 0
        while abs(Fb - F_req) / F_req > r5.CONVERGE and n_it < 4:
            NIn = math.exp(math.log(NIb) + (math.log(F_req) - math.log(Fb)) * (math.log(NIb) - math.log(NIa)) / (math.log(Fb) - math.log(Fa)))
            NIa, Fa = NIb, Fb; sb_, fb_ = solve(NIn); NIb, Fb = NIn, sb_["F"]; n_it += 1
        NI_req, s_req, f_req = NIb, sb_, fb_
        curve = [(NI_req, s_req, f_req)]
        for k in (1.25, 1.5, 2.0):
            sk, fk = solve(k * NI_req); curve.append((k * NI_req, sk, fk))
    li = dict(launch_profile(s_req["W"], F_req), F_req=F_req)
    out = dict(tag=tag, Wp=Wp, T=T, Wl=Wl, w=w, hgt=hgt, h=h, NI_req=NI_req, F=s_req["F"], sat=s_req["sat"], B_plate=s_req["B_plate"], W_mag=s_req["W"], faces={}, n_it=n_it)
    NIs = np.array([c[0] for c in curve]); Fs = np.array([c[1]["F"] for c in curve])
    Phi = 2 * s_req["W"] / NI_req; V_turn = Phi * L_FAST / TAU; V = N_TURNS * V_turn
    for face in ("h", "v", "best", "aligned"):
        wi = [worst_ic(c[2], face, theva) for c in curve]
        ic_req, B, th = wi[0]
        tapes_leg = NI_req / (0.5 * ic_req); tape = 2 * tapes_leg * 0.003 * bd.TAPE_T * bd.TAPE_RHO
        share = (NIs / tapes_leg) / np.array([x[0] for x in wi])
        if share.max() >= 0.7:
            NI70 = float(np.interp(0.7, share, NIs)); F70 = math.exp(float(np.interp(math.log(NI70), np.log(NIs), np.log(Fs))))
        else:
            F70 = float("nan")
        layers = (NI_req / N_FAM) / (0.5 * ic_req) / N_TURNS
        Ep, Es = tape_loss_launch(f_req, face, tapes_leg, layers, li)
        pk = pack(w, face, math.ceil(layers), V, V_turn) if face in ("h", "v") else None
        Bn, Bt = zip(*[face_components(Bx, By, face) for Bx, By in f_req])
        out["faces"][face] = dict(ic=ic_req, B=B, th=th, tapes_leg=tapes_leg, tape=tape, margin=F70 / F_req, layers=layers, Ep=Ep, Es=Es, pack=pk,
                                  Bn_mean=float(np.concatenate(Bn).mean()), Bt_mean=float(np.concatenate(Bt).mean()))
    out.update(V=V, V_turn=V_turn, I=NI_req / N_FAM / N_TURNS, seconds=time.time() - t0)
    return out


SHAPES_TRY = [(0.200, 0.050, "the 200 x 50 mm bundle"), (0.150, 0.067, "taller, narrower"), (0.100, 0.100, "square"), (0.067, 0.150, "tall"), (0.300, 0.034, "wider, flatter"), (0.240, 0.042, "a little wider")]


def part_bundle(jobs=4, h=0.005, plates=((1.0, 0.15, 0.6), (1.0, 0.20, 0.6))):
    import multiprocessing as mp
    hr(f"S2.2  ONE ATTEMPT AT THE BUNDLE: SHAPES AND STACK ORIENTATIONS (3 mm tape, 10 turns, 20 overlapping coils; {h*1e3:g} mm cells)")
    print(f"  Each shape is solved for the requirement on the plate (the convergence step), then at 1.25, 1.5 and 2 times its ampere-turns for the rating.  For each stack")
    print(f"  orientation: the worst tape at the requirement (its field, the angle from its normal, its limit), the tape at half that limit, the rating margin (the worst tape")
    print(f"  at 70 percent), the pulse loss of the heaviest launch, and whether the 20 cells of 10 stacks pack in the shape's width (the films sized for the turns they")
    print(f"  separate, the walls for the terminal voltage at {E_INS/1e6:.0f} kV/mm).  Orientations: faces horizontal (h), faces vertical (v), each stack the better of the two")
    print(f"  ('along the field', as two orientations allow), and every stack turned to its local field (the floor: no field through any face).")
    tasks = [(Wp, T, Wl, w, hg, h, lab) for (Wp, T, Wl) in plates for (w, hg, lab) in SHAPES_TRY]
    p = os.path.join(LOGS, "round6b_bundle.json")
    done = json.load(open(p, encoding="utf-8")) if os.path.exists(p) else {}
    def key(t): return f"{t[0]:.2f}|{t[1]:.3f}|{t[2]:.2f}|{t[3]:.3f}|{t[4]:.3f}|{t[5]:.4f}"
    todo = [t for t in tasks if key(t) not in done]
    if todo:
        with mp.Pool(jobs) as pool:
            for t, res in zip(todo, pool.imap(_shape_worker, todo)):
                done[key(t)] = res; save("round6b_bundle.json", done)
                print(f"    solved {t[6]} on {t[0]:.1f} x {t[1]*1e3:.0f}: {res['seconds']:.0f} s", flush=True)
    for (Wp, T, Wl) in plates:
        print(f"\n  On the {Wp:.1f} m x {T*1e3:.0f} mm plate with a {Wl:.1f} m loop (between leg centers):")
        print(f"  {'bundle W x H mm':>16} {'kA':>4} {'plate sat %':>11} {'cryostat width m':>16} | {'faces':>8} {'worst tape: T, deg from normal, A':>34} {'tape kg/m':>9} {'margin':>6} {'loss GJ (strip + slab)':>24} {'B through / along, mean T':>26} | packs in the width as; height needed mm")
        base = None
        for (w, hg, lab) in SHAPES_TRY:
            res = done[key((Wp, T, Wl, w, hg, h, lab))]
            for face in ("h", "v", "best", "aligned"):
                f = res["faces"][face]; pk = f["pack"]
                if base is None and face == "h":
                    base = f["Ep"] + f["Es"]
                ptxt = "-" if face in ("best", "aligned") else ("does not pack" if pk is None else
                        f"{pk['n_cx']} x {pk['n_cy']} cells, {pk['n_x']} x {pk['n_y']} stacks, {pk['order']} order: {pk['h_need']*1e3:.1f} of {hg*1e3:.0f} ({'fits' if pk['h_need'] <= hg + 1e-9 else 'DOES NOT FIT'})")
                print(f"  {(f'{w*1e3:.0f} x {hg*1e3:.0f}' if face == 'h' else ''):>16} {(f'{res['NI_req']/1e3:.0f}' if face == 'h' else ''):>4} {(f'{res['sat']*100:.1f}' if face == 'h' else ''):>11} {(f'{Wl + w:.2f}' if face == 'h' else ''):>16} | "
                      f"{face:>8} {f'{f['B']:.2f} T, {f['th']:.0f} deg, {f['ic']:.0f} A':>34} {f['tape']:9.1f} {f['margin']:6.2f} {f'{(f['Ep']+f['Es'])/1e9:.0f} ({f['Ep']/1e9:.0f} + {f['Es']/1e9:.0f})':>24} "
                      f"{f'{f['Bn_mean']:.2f} / {f['Bt_mean']:.2f}':>26} | {ptxt}")
    return done


# =============================================================================
# S2.3  the lost support track with ten families connected
# =============================================================================

def part_lost():
    rs, rl, fin = load_final()
    hr("S2.3  THE LOST SUPPORT TRACK WITH TEN FAMILIES CONNECTED")
    sc = support_coil(rs); I0 = sc["I0"]
    fbk = r5.fast_basics(rl, N_TURNS); I_full = fbk["I"]
    A_p = r5.lead_rows(rl)[0][N_TURNS]["A_m"]; R_lead = be.RHO_CU_77 * 1.0 / A_p; P_full = I_full**2 * R_lead      # round 5's lead (21 mm2), as round 6a
    m_lead = 8960.0 * A_p * 1.0; E20 = 20 * 400.0 * m_lead                   # J for a 20 K adiabatic rise of one lead (round 6a)
    P_steady = 20 / (1 / (8 * 400 * A_p))
    W = 2 * W_HALF; F_full = fb.F_LOAD_TRACK
    x_s = SIDE / 2 - rl["Wp"] - 0.2 - rs["Wp"] / 2; x_c = SIDE / 2 - rl["Wp"] / 2
    print(f"  The numbers: the casing's weight {W/1e3:.1f} kN/m, {W_HALF/1e3:.1f} per support track; the support tracks at +-{x_s:.3f} m from the centerline, the load and lift tracks at +-{x_c:.3f} m;")
    print(f"  a load track's full pull {F_full/1e3:.0f} kN/m at {I_full:,.0f} A per coil; the pulse-sized lead {A_p*1e6:.0f} mm2, {P_full:.0f} W of Joule heat at full current, {E20:.0f} J for a 20 K rise,")
    print(f"  {P_steady:.2f} W steadily.  With n of the 20 families connected a pull p (of full) heats each connected lead at p (20/n)^2 x {P_full:.0f} W, so the 20 K budget is")
    print(f"  an integral: p dt <= {E20/P_full/1:.2f} s x (n/20)^2: {E20/P_full:.2f} s with 20 connected, {E20/P_full/4:.3f} s with ten.")
    # (A) the heave-only picture of the earlier rounds
    print()
    print(f"  (A) Heave only, as rounds 4 to 6a carried the case (one meter of casing, no roll).  The two load tracks share the lost track's {W_HALF/1e3:.1f} kN/m: {W_HALF/2/1e3:.1f} each,")
    p0 = W_HALF / 2 / F_full
    print(f"  p = {p0*100:.1f} percent, held for {E20/P_full/p0:.1f} s with 20 connected and {E20/P_full/4/p0:.1f} s with ten (round 6a's 21 and 5.4 s).  The surviving track must go from {W_HALF/1e3:.1f} to {W/1e3:.1f} kN/m.")
    print(f"  The force goes as the ampere-turns squared, so its current rises by the square root of two, {I0:.0f} to {I0*math.sqrt(2):.0f} A ({rs['NI_req']*math.sqrt(W_HALF/rs['F_req'])/1e3:.0f} to {rs['NI_req']*math.sqrt(W/rs['F_req'])/1e3:.0f} kA; Stage B's 141 to 199 kA),")
    print(f"  not by a factor of two: round 6a's '21 s to double its current' overstated the step by 2.4.  While it rises the load tracks' pull falls, and the lead budget")
    print(f"  is the integral of p: with the current rising linearly over T after a delay t_c (strips opened, contacts closed), integral p dt = p0 (t_c + {1 - (0.828/2 + 0.1716/3):.3f} T).")
    dI = I0 * (math.sqrt(2) - 1); k_int = 1 - (2 * (math.sqrt(2) - 1) / 2 + (math.sqrt(2) - 1)**2 / 3)
    print(f"  {'supply V':>9} {'L H':>6} {'T s':>6} | {'20 connected: t_c may be s':>26} | {'ten connected: t_c may be s':>27} | {'peak power per coil kW':>22} {'energy per coil MJ':>18}")
    rowsA = []
    for Lc, tag in ((sc["Lc_scaled"], "scaled"), (sc["Lc_energy"], "energy")):
        for V in (sc["Lc_scaled"] * 0.01 * I0 if tag == "scaled" else sc["Lc_energy"] * 0.01 * I0, 500.0, 1000.0, 2000.0, 3000.0):
            T = Lc * dI / V
            tc20 = E20 / P_full / p0 - k_int * T; tc10 = E20 / P_full / 4 / p0 - k_int * T
            rowsA.append(dict(L=tag, Lc=Lc, V=V, T=T, tc20=tc20, tc10=tc10))
            print(f"  {V:9.0f} {Lc:6.1f} {T:6.1f} | {(f'{tc20:.1f}' if tc20 > 0 else 'does not close'):>26} | {(f'{tc10:.1f}' if tc10 > 0 else 'does not close'):>27} | {V*I0*math.sqrt(2)/1e3:22.0f} {0.5*Lc*((I0*math.sqrt(2))**2 - I0**2)/1e6:18.2f}")
    V10 = lambda Lc, tc: Lc * dI * k_int / (E20 / P_full / 4 / p0 - tc)
    print(f"  With ten connected and 1 s to open the strips and close the contacts the supply must give {V10(sc['Lc_scaled'], 1.0):,.0f} V (scaled L) or {V10(sc['Lc_energy'], 1.0):,.0f} V (energy L); with 2 s, {V10(sc['Lc_scaled'], 2.0):,.0f} or {V10(sc['Lc_energy'], 2.0):,.0f} V.")
    # (B) with the roll
    print()
    print(f"  (B) With the roll, which no round has carried for this case (every 'one support track lost' case is the point model; the roll model has only the lift-alone case).")
    print(f"  One support track at +{x_s:.3f} m cannot hold a casing whose weight acts at the centerline: the moment S x {x_s:.3f} m has to be balanced.  The only track that can is")
    print(f"  the load track on the lost side (it pulls up at -{x_c:.3f} m); the lift tracks pull down and the other load track is on the wrong side.  Statics, per meter:")
    S_now = W_HALF
    # immediately: S = W_HALF (frozen flux), unknown load pulls a (lost side) and b (surviving side): a + b = W - S; a x_c = S x_s + b x_c
    a0 = ((W - S_now) + S_now * x_s / x_c) / 2; b0 = (W - S_now) - a0
    # settled: b = 0: S + a = W, a x_c = S x_s
    S1 = W / (1 + x_s / x_c); a1 = W - S1
    print(f"    at the loss (the surviving track still at {S_now/1e3:.1f} kN/m): the lost side's load track {a0/1e3:.1f} kN/m ({a0/F_full*100:.1f} percent of full), the other load track {b0/1e3:.1f} ({b0/F_full*100:.1f} percent);")
    print(f"    settled (the surviving track raised until the other load track is idle): the surviving support track {S1/1e3:.1f} kN/m ({S1/W_HALF:.2f} x its normal force, {math.sqrt(S1/W_HALF):.3f} x its")
    print(f"    current, {I0:.0f} to {I0*math.sqrt(S1/W_HALF):.0f} A) and the lost side's load track {a1/1e3:.1f} kN/m, {a1/F_full*100:.1f} percent of its full pull, for as long as the support track is out.")
    print(f"    (The alternative, the surviving side's lift track pulling down: S = {W/(1 - x_s/x_c)/1e3:.0f} kN/m against the cap of {rs['F_req']/1e3:.0f}, the lift track {W/(1 - x_s/x_c)*x_s/x_c/1e3:.0f} kN/m: worse.)")
    print(f"  {'n connected on that load track':>30} {'current per coil, % of full':>27} {'A':>6} {'Joule per lead W':>16} {'20 K in s':>9} {'against the steady':>18} {P_steady:.2f} W")
    rowsB = []
    for n in (20, 10):
        for lab, a in (("at the loss", a0), ("settled", a1)):
            p = a / F_full; f = math.sqrt(p) * 20 / n; PJ = p * (20 / n)**2 * P_full
            rowsB.append(dict(n=n, when=lab, a=a, p=p, f=f, PJ=PJ, t20=E20 / PJ))
            print(f"  {f'{n}, {lab}':>30} {f*100:27.0f} {f*I_full:6.0f} {PJ:16.0f} {E20/PJ:9.1f} {PJ/P_steady:18.0f} x")
    dI_B = I0 * (math.sqrt(S1 / W_HALF) - 1)
    print(f"  The surviving track's step is small ({dI_B:.0f} A: {sc['Lc_scaled']*dI_B/1000:.1f} s at 1 kV, {sc['Lc_scaled']*dI_B/(sc['Lc_scaled']*0.01*I0):.0f} s at the 1 percent per second supply), but nothing the support's supply does removes the")
    print(f"  {a1/1e3:.0f} kN/m from the lost side's load track.  Its pulse-sized leads hold that for {E20/(a1/F_full*P_full):.0f} s with all 20 families connected and {E20/(a1/F_full*4*P_full):.1f} s with ten.")
    # a lead that carries the settled pull steadily
    print(f"  A lead that carries it steadily (the optimum conduction-cooled lead for the current, bearing_electrical.lead_heat):")
    print(f"  {'m families carry it':>19} {'A per coil':>10} {'lead mm2 (1 m)':>14} {'W per lead, carrying':>20} {'W per lead, idle':>16} {'against the pulse lead s':>24} {be.lead_heat(0.0, A_over_l=A_p)[0]:.2f} W")
    rowsL = []
    for m in (20, 10, 5):
        f = math.sqrt(a1 / F_full) * 20 / m; I = f * I_full
        if f > 1.21:
            print(f"  {m:19d} {I:10,.0f}  over the coil's rating (1.21 x the full current)"); continue
        q, aol = be.lead_heat(I); q0 = be.lead_heat(0.0, A_over_l=aol)[0]
        rowsL.append(dict(m=m, I=I, A=aol, q=q, q0=q0))
        print(f"  {m:19d} {I:10,.0f} {aol*1e6:14.0f} {q:20.1f} {q0:16.1f}")
    save("round6b_lost.json", dict(coil=sc, x_s=x_s, x_c=x_c, P_full=P_full, E20=E20, P_steady=P_steady, A=rowsA, B=rowsB, leads=rowsL, a0=a0, b0=b0, S1=S1, a1=a1, dI=dI, k_int=k_int, p0=p0))
    return rowsA, rowsB


def main():
    part = "all"; jobs = 4
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    print("=" * 150 + "\n  THE BEARING, ROUND 6 PART B\n" + "=" * 150)
    print(f"  Ring at 700 km, rising quadratic: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m with hardware {ring.M_HARDWARE_M:,.0f}, side {SIDE:.2f} m; the rise {TAU*1e3:.0f} ms; {N_TURNS} turns.")
    t0 = time.time()
    if part in ("legs", "all"): part_legs()
    if part in ("dump", "all"): part_dump()
    if part in ("film", "all"): part_film()
    if part in ("plates", "all"): part_plates(jobs)
    if part in ("bundle", "all"): part_bundle(jobs)
    if part in ("lost", "all"): part_lost()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
