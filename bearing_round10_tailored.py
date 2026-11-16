#!/usr/bin/env python3
"""
The tailored winding for a local load (round 10, item 2, 2026-11-17: Vol_III/Claude_Code_instructions_round10_2026-11-17.md, Section 2).
Book: Orbital Ring Engineering, Volume III: Chapter 7 (Section 7.7), Chapter 16.

Round 9 (bearing_round9_local.py) found that the staggered 1 km support coils cannot carry a parked 5,000 tonne spacecraft on 100 m.  This script
sizes the winding that was named there and not run: a winding only as long as the load, ending at its edges.

The kind taken: extra turns INSIDE the four existing 150 x 67 mm bundles of the split support track, spread evenly among the base coils' turns, on the
same 1.2 m x 100 mm plate.  Reason: in the 2D section this is the same bundle box with more ampere-turns in it, so the solved 2.5 mm curve of the
split support (round 8, bearing_round8.py --part=solve-sup: pull, plate saturation, loop flux and worst tape against the ampere-turns per leg) is its
exact 2D answer.  A winding beside the bundles would be a new section and a new set of solves; none was run.  No new field solve is made here.

Cases: A, 198 kN/m per track more on 100 m; B, the same on 1 km (a factory).  The base 1 km coils stay at their normal current under the winding.
  part 1  the ampere-turns, the tape, the plate, the room in the bundle, and what the plate allows (item 1 and item 5)
  part 2  persistence and charging: the volt-seconds from the solved loop flux (item 2), and what the charging does to the persistent base coils
  part 3  the leftover along the track and the gap's movement by bearing_round9_local.movement (item 3)
  part 4  one sled pass, the point model as bearing_round9_local.part_pass (item 4)   [--pass]
Not modeled: the winding's ends (the field of a coil end is 3D; the 2D curve is applied cell by cell along the track).
Run (from Python/orbitalring): python bearing_round10_tailored.py [--pass]     writes runs/round10_logs/round10_tailored.txt and .json
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq, least_squares

import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round7 as r7
import bearing_round8 as r8
import bearing_round9 as r9
import bearing_round9_local as r9l

LOGS = os.path.join(r7.OUT, "round10_logs")
W_T = r9l.W_T; F_LOAD = r9l.F_LOAD
T_PER_M = 50.0                                       # tonnes per meter of ring that 2 x F_LOAD stands for (5,000 t on 100 m)
VOLTS = (210.0, 1000.0)
SHARE_CAP, SHARE_RATE = 0.5, 0.7


class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8"); self.o = sys.stdout
    def write(self, s):
        self.o.write(s); self.f.write(s)
    def flush(self):
        self.o.flush(); self.f.flush()


def curves():
    """The solved curve of the split support at 2.5 mm cells (round 8's points): pull, loop flux per meter and turn, worst 12 mm tape at 70 K, plate saturation."""
    F_of, n_t, I_n, I_cap, NIs, Fs = r9l.curve()
    rows = r8.load("round8_rows.json"); pts = rows["support_points"]
    phis = np.array([p["phi"] for p in pts]); ics = np.array([p["ic"] for p in pts]); sats = np.array([p["sat"] for p in pts])
    half = r7.load("round7_half_2p5.json"); Lin = half["L11"] + half["L21"]; NI_n = NIs[0]; c = 1 - phis[0] / (Lin * NI_n)

    def phi_of(NI):
        NI = np.asarray(NI, dtype=float); a = np.abs(NI)
        low = Lin * a * (1 - c * a / NI_n)                                   # the linear plate's flux, bent to meet the first solved point
        hi = np.exp(np.interp(np.log(np.maximum(a, 1.0)), np.log(NIs), np.log(phis)))
        return np.sign(NI) * np.where(a < NI_n, low, hi)
    ic_of = lambda NI: np.interp(NI, NIs, ics)                               # linear between the solved points, as round 8's rating
    sat_of = lambda NI: np.interp(NI, NIs, sats)
    return dict(F_of=lambda NI: F_of(np.abs(NI)), phi_of=phi_of, ic_of=ic_of, sat_of=sat_of, n_t=n_t, I_n=I_n, I_cap=I_cap, NIs=NIs, Fs=Fs, NI_n=NI_n,
                tapes_base=rows["tapes_leg"], rating=rows["rating_support"], I_tape_n=pts[0]["I_tape"], Lin=Lin)


def kg_per_m(tapes_leg):
    return 4 * tapes_leg * 0.012 * bd.TAPE_T * bd.TAPE_RHO                   # four legs to a track (two loops)


def part1(c):
    print("\n  PART 1.  THE WINDING: extra turns in the four bundles of a support track, only as long as the load.  Per support track; the base coils at their normal current.")
    F_of, ic_of, sat_of, NI_n, NIs = c["F_of"], c["ic_of"], c["sat_of"], c["NI_n"], c["NIs"]
    A_box = bo.shape_area(bo.BEST_SUPPORT["shape"]); a_tape = 0.012 * bd.TAPE_T
    fill = lambda tapes: tapes * a_tape / be.F_CABLE / A_box                 # the tape stacks' area (tape thickness over pitch F_CABLE) over the bundle's box

    def state(NI, tapes_t=None, share=SHARE_CAP):
        ic = float(ic_of(NI)); t = (NI - NI_n) / (share * ic) if tapes_t is None else tapes_t
        return dict(NI=NI, NI_tail=NI - NI_n, F=float(F_of(NI)), load=float(F_of(NI)) - W_T, t_per_m=(float(F_of(NI)) - W_T) / F_LOAD * T_PER_M, sat=float(sat_of(NI)), ic=ic,
                    tapes_tail=t, share_tail=(NI - NI_n) / t / ic, share_base=c["I_tape_n"] / ic, kg=kg_per_m(t), fill=fill(c["tapes_base"] + t), fill_base=fill(c["tapes_base"]))
    NI_req = brentq(lambda NI: F_of(NI) - (W_T + F_LOAD), NIs[0], NIs[-1])
    req = state(NI_req)
    # the rating of that winding: its worst tape at 70 percent, the base coils still at their normal current
    NI_r = brentq(lambda NI: (NI - NI_n) / req["tapes_tail"] - SHARE_RATE * ic_of(NI), NI_req, NIs[-1])
    rate = state(NI_r, req["tapes_tail"])
    # the tape that keeps the track's rating of round 8 (341 kN/m with the worst tape at 70 percent)
    NI_341 = c["rating"]["NI"]; keep = state(NI_341, share=SHARE_RATE)
    # what the room in the bundle allows with the tape at 50 percent: the stacks filling the box
    NI_room = brentq(lambda NI: state(NI)["fill"] - 1.0, NI_req, NIs[-1]) if state(NIs[-1])["fill"] > 1.0 else NIs[-1]
    room = state(NI_room); top = state(NIs[-1])
    print(f"  The solved curve (bearing_round8.py --part=solve-sup, 2.5 mm cells): from {NIs[0]/1e3:.1f} to {NIs[-1]/1e3:.1f} kA per leg.  Base coils: {c['tapes_base']:,.0f} tapes per leg, {NI_n/1e3:.1f} kA per leg, {c['I_tape_n']:.1f} A per tape.")
    print(f"  Bundle box {A_box*1e6:,.0f} mm2; a 12 mm tape {a_tape*1e6:.2f} mm2 ({bd.TAPE_T*1e3:.2f} mm thick [UNVERIFIED]); stacks at tape thickness over pitch {be.F_CABLE:.2f} [assumption of bearing_electrical].")
    hdr = f"  {'state':72} {'kA per leg, all':>15} {'the winding kA':>14} {'pull kN/m':>9} {'standing t/m':>12} {'plate sat %':>11} {'worst tape A':>12} {'winding tapes/leg':>17} {'its share %':>11} {'base share %':>12} {'tape kg/m':>9} {'box filled %':>12}"
    print(hdr)
    rows = [("the load: 198 kN/m more per track, the winding's tape at 50 percent", req), ("the same winding at its rating (its worst tape at 70 percent)", rate),
            ("the winding sized to keep the track's 341 kN/m rating at 70 percent", keep), ("the bundle's box full of tape stacks, the winding at 50 percent", room),
            ("the top of the solved curve, the winding at 50 percent (does not fit)", top)]
    for lab, s in rows:
        print(f"  {lab:72} {s['NI']/1e3:15.1f} {s['NI_tail']/1e3:14.1f} {s['F']/1e3:9.1f} {s['t_per_m']:12.1f} {s['sat']*100:11.1f} {s['ic']:12.1f} {s['tapes_tail']:17,.0f} {s['share_tail']*100:11.1f} {s['share_base']*100:12.1f} {s['kg']:9.1f} {s['fill']*100:12.1f}")
    print(f"  (the base coils alone fill {req['fill_base']*100:.1f} percent of the box by the same count; the base tape is {kg_per_m(c['tapes_base']):.1f} kg/m per track.  The fill counts tape stacks only: no formers,")
    print(f"  insulation between coils or cooling passages, which no script sizes [UNVERIFIED that the turns fit].  The worst tape is the bundle's worst place at the solved state; the base")
    print(f"  coils' tapes there carry their normal current, so their share is lower than the winding's.)")
    return dict(req=req, rate=rate, keep=keep, room=room, top=top, A_box=A_box)


def part2(c, p1):
    print("\n  PART 2.  PERSISTENCE AND CHARGING.  The loop flux per meter and turn from the same solves; both loops of the track and all the winding's coils rising together.")
    phi_of, n_t, NI_n = c["phi_of"], c["n_t"], c["NI_n"]; req = p1["req"]
    dphi = float(phi_of(req["NI"]) - phi_of(NI_n)); I_coil = 10 * SHARE_CAP * req["ic"]          # a coil made as the support coils: n_t turns of ten tapes
    turns = req["tapes_tail"] / 10.0; n_coils = turns / n_t
    print(f"  The winding per loop: {req['tapes_tail']:,.0f} tapes per leg, as turns of ten tapes {turns:.1f} turns at {I_coil:,.0f} A; made as the support coils ({n_t:.2f} turns each) that is {n_coils:.2f} coils per loop.")
    print(f"  The loop flux goes from {float(phi_of(NI_n)):.4f} to {float(phi_of(req['NI'])):.4f} Wb per meter and turn ({dphi:.4f}) as the winding goes from zero to full with the base coils held at normal.")
    print(f"  {'case':28} {'length m':>8} {'V s per coil of ' + format(n_t, '.2f') + ' turns':>30} {'effective H per coil':>20} {'s at 210 V':>10} {'s at 1,000 V':>12} {'V s, a loop winding in series':>29} {'s at 210 V':>10} {'s at 1,000 V':>12} {'V s put back into each base coil':>32} {'s at 210 V':>10}")
    out = {}
    for tag, L in (("A, 100 m", 100.0), ("B, 1 km factory", 1000.0)):
        lam = n_t * L * dphi; lam_all = turns * L * dphi
        lam_base = n_t * min(L, 1000.0) * dphi                                                   # a base coil lying wholly over the winding (B) or covering all of it (A)
        out[tag] = dict(L=L, lam=lam, L_eff=lam / I_coil, t=[lam / V for V in VOLTS], lam_all=lam_all, t_all=[lam_all / V for V in VOLTS], lam_base=lam_base, t_base=lam_base / VOLTS[0])
        print(f"  {tag:28} {L:8.0f} {lam:30,.0f} {lam / I_coil:20.2f} {lam / VOLTS[0]:10.1f} {lam / VOLTS[1]:12.2f} {lam_all:29,.0f} {lam_all / VOLTS[0]:10.1f} {lam_all / VOLTS[1]:12.1f} {lam_base:32,.0f} {lam_base / VOLTS[0]:10.1f}")
    lam_own = n_t * 1000.0 * float(phi_of(NI_n))
    print(f"  (a base coil's own flux at normal current is {lam_own:,.0f} V s.  The base coils are closed superconducting loops, so the winding's flux through them lowers their current unless")
    print(f"  their supplies put the same volt-seconds back: the last two columns, for a base coil that covers the whole winding (A) or lies wholly over it (B).  Part 3 runs both ways.")
    print(f"  The ends' inductance is not in the 2D flux.  The switch, the joints and the leads of such a coil are in no script.)")
    out.update(dphi=dphi, I_coil=I_coil, turns=turns, n_coils=n_coils, lam_own=lam_own)
    return out


def profile(c, NI_tail, L, s, window, shift=0.0, ramp=0.0):
    """Cells of s on a periodic window; the load on |x| < L/2; the winding on |x - shift| < L/2, its ampere-turns rising linearly over 'ramp' centered on each end."""
    N = int(round(window / s)); x = (np.arange(N) - N / 2 + 0.5) * s
    w = np.where(np.abs(x) < L / 2, W_T + F_LOAD, W_T)
    d = L / 2 - np.abs(x - shift)                                            # distance inside the winding's end
    cov = np.where(d > 0, 1.0, 0.0) if ramp <= 0 else np.clip(d / ramp + 0.5, 0.0, 1.0)
    F = c["F_of"](c["NI_n"] + NI_tail * cov)
    return dict(N=N, s=s, x=x, w=w, F=F, r=F - w)


def flux_frozen(c, NI_tail, L, window):
    """The base coils left as closed loops while the winding is charged (not put back on their supplies): each keeps its flux linkage.  Cells of 100 m,
    ten staggered 1 km coils over every cell, as bearing_round9_local.solve."""
    s = 100.0; N = int(round(window / s)); m = int(round(L / s)); c0 = N // 2 - m // 2; n_t = c["n_t"]
    A = np.zeros((N, N))
    for j in range(N):
        for q in range(10):
            A[(j + q) % N, j] = 1.0
    tail = np.zeros(N); tail[c0:c0 + m] = NI_tail
    phi_n = float(c["phi_of"](c["NI_n"]))
    res = least_squares(lambda I: (A.T @ c["phi_of"](n_t * (A @ I) + tail)) / (10 * phi_n) - 1.0, np.full(N, c["I_n"]), xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=2000)
    I = res.x; NI = n_t * (A @ I) + tail; F = c["F_of"](NI)
    w = np.full(N, W_T); w[c0:c0 + m] += F_LOAD
    x = (np.arange(N) - c0 - m / 2 + 0.5) * s
    return dict(N=N, s=s, x=x, w=w, F=F, r=F - w, I=I, NI=NI, resid=float(np.abs(res.fun).max()), sl=slice(c0, c0 + m))


def row_of(lab, sol):
    mv, par = r9l.movement(sol); fh = mv["fast tracks holding"]; rel = fh["rel"]; hold = fh["hold"] / 2; r = sol["r"]
    big = np.abs(r) > 1e3
    return dict(lab=lab, r_min=float(r.min()), r_max=float(r.max()), len_over=float(big.sum() * sol["s"]), net=float(r.sum() * sol["s"]), rel=float(np.abs(rel).max()), hold=float(np.abs(hold).max()), K=par["K"], k_reg=par["k_reg"])


def part3(c, p1):
    print("\n  PART 3.  THE LEFTOVER ALONG THE TRACK AND THE GAP'S MOVEMENT (bearing_round9_local.movement: the casing a Timoshenko beam, the cable a string, both support tracks'")
    print("  constant-flux pull, the fast tracks holding by the protection's 5 Hz law [UNVERIFIED as the standing regulation's law]).  Per support track.  The 2D curve cell by cell;")
    print("  the winding's ends are not modeled: a sharp end is a step of the ampere-turns, and the 'end spread' rows put a linear rise in its place [UNVERIFIED, 3D].")
    NI_tail = p1["req"]["NI_tail"]; out = {}
    old = r9.load("round9_local.json")["rows"][0]
    print(f"  {'case':10} {'what':86} {'leftover, lowest kN/m':>21} {'highest kN/m':>12} {'on m of track':>13} {'net kN':>8} {'gap moves mm':>12} {'a load track holds kN/m':>23}")
    print(f"  {'A':10} {'round 9, no winding: the 1 km coils at their best currents (round9_local.json)':86} {-old['short']/1e3:21.1f} {old['over']/1e3:12.1f} {'about 2,000':>13} {old['net']/1e3:8.0f} {old['rel_max']*1e3:12.2f} {old['hold_max']/1e3:23.2f}")
    for tag, L, window in (("A", 100.0, 3000.0), ("B", 1000.0, 6000.0)):
        rows = []
        rows.append(row_of("the winding ends at the load's edges, base coils at normal (1 m cells)", profile(c, NI_tail, L, 1.0, window)))
        for sh in (1.0, 5.0, 10.0):
            rows.append(row_of(f"the load parked {sh:.0f} m off the winding (1 m cells)", profile(c, NI_tail, L, 1.0, window, shift=sh)))
        for rp in (0.3, 1.0):
            rows.append(row_of(f"end spread: the winding's ampere-turns rise over {rp:.1f} m at each end (0.05 m cells)", profile(c, NI_tail, L, 0.05, window, ramp=rp)))
        ff = flux_frozen(c, NI_tail, L, window); rw = row_of("the base coils NOT put back on their supplies: constant flux (100 m cells)", ff)
        nb = ff["NI"].copy(); nb[ff["sl"]] -= NI_tail                                             # the base coils' own ampere-turns per leg, cell by cell
        rw.update(NIb_min=float(nb.min()), NIb_max=float(nb.max()), NIb_load=float(nb[ff["sl"]].mean()), F_load=float(ff["F"][ff["sl"]].mean()), resid=ff["resid"]); rows.append(rw)
        for r in rows:
            print(f"  {tag:10} {r['lab']:86} {r['r_min']/1e3:21.1f} {r['r_max']/1e3:12.1f} {r['len_over']:13.1f} {r['net']/1e3:8.0f} {r['rel']*1e3:12.3f} {r['hold']/1e3:23.2f}")
        print(f"  {'':10} (constant flux: the base coils' ampere-turns per leg go from {c['NI_n']/1e3:.1f} kA to {rw['NIb_load']/1e3:.1f} kA under the load and between {rw['NIb_min']/1e3:.1f} and {rw['NIb_max']/1e3:.1f} kA along the window; the pull under the load is {rw['F_load']/1e3:.1f} of {(W_T+F_LOAD)/1e3:.1f} kN/m; flux residual {rw['resid']:.1e})")
        out[tag] = rows
    k = out["A"][0]
    print(f"  (stiffness per meter of ring in the movement: the support's constant-flux pull {k['K']/1e6:.2f} MN/m2 as a mean over the window, the fast tracks' law {k['k_reg']/1e6:.1f} MN/m2.  With the base coils")
    print(f"  at their normal current everywhere the stagger of the 1 km coils puts nothing into the leftover: every cell has the same ten coils' ampere-turns.  'On m of track' counts cells")
    print(f"  with more than 1 kN/m of leftover.  A load track's leads carry {r9.load('round9_local.json')['p_steady']/1e3:.2f} kN/m steadily (round 9).)")
    return out


def part4(p3):
    """One sled pass, as bearing_round9_local.part_pass: the point model, the pair, the 10 mm protection, the casing's mass and weight at the place with the
    parked load (the same per meter in case A and case B).  The leftover is put in as a steady support that differs from the weight (F_STEADY = W + 2 r)."""
    import sled_passage_pulsed as pz
    import bearing_round9_parts as p9
    pz_, sp_, ct, tabs, lane, cc = p9._setup()
    old = pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING
    m_load = T_PER_M * 1e3
    left = [("no leftover: the winding ends at the load's edges, base coils at normal", 0.0),
            ("a point under the load with no winding at all, taken alone (no bending, no neighbors)", -F_LOAD)]
    rows = []
    try:
        for u in (30e3, 1e3):
            base, _ = pz.point_pulsed([lane(u, True, label="bare")], u, tabs)
            rows.append(dict(u=u, lab="bare casing, no load", r=0.0, exc=base[0]["max_exc_all"], gap=min(base[0]["min_gap_b"], base[0]["min_gap_t"]), contact=bool(base[0]["contact"]), pk_up=base[0].get("pk_hu"), pk_dn=base[0].get("pk_hd")))
            for lab, r in left:
                pz.MU = old[0] + m_load; pz.W = old[1] + 2 * F_LOAD; pz.F_STEADY = pz.W + 2 * r; pz.MATCH_RING = pz.MU_C / (pz.MU + pz.MU_C)
                res, _ = pz.point_pulsed([lane(u, True, label=lab)], u, tabs)
                pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING = old
                rows.append(dict(u=u, lab=lab, r=r, exc=res[0]["max_exc_all"], gap=min(res[0]["min_gap_b"], res[0]["min_gap_t"]), contact=bool(res[0]["contact"]), pk_up=res[0].get("pk_hu"), pk_dn=res[0].get("pk_hd")))
    finally:
        pz.MU, pz.W, pz.F_STEADY, pz.MATCH_RING = old
    print(f"\n  PART 4.  ONE SLED PASS OVER THE PLACE (sled_passage_pulsed.point_pulsed as bearing_round9_local.part_pass: the pair, the heaviest sled, the 10 mm protection; the loaded place")
    print(f"  {old[0]+m_load:,.0f} kg/m and {(old[1]+2*F_LOAD)/1e3:.0f} kN/m, the support at constant flux; the same per meter for case A and case B, so one run stands for both):")
    print(f"  {'u km/s':>7} {'the place':84} {'leftover per track kN/m':>23} {'movement mm':>11} {'smallest gap mm':>15} {'largest pull up kN/m':>20} {'down':>8}")
    for r in rows:
        f = lambda v: f"{v/1e3:.0f}" if v is not None else "n/a"
        print(f"  {r['u']/1e3:7.0f} {r['lab']:84} {r['r']/1e3:23.1f} {r['exc']*1e3:11.2f} {r['gap']*1e3:15.1f}{' CONTACT' if r['contact'] else ''} {f(r['pk_up']):>20} {f(r['pk_dn']):>8}")
    print("  (not in it: the casing's bending between the loaded stretch and the bare casing beside it; the winding's ends; any standing regulator, so a point with a leftover falls")
    print("  until the 10 mm protection takes it; the protection's pull is capped at the pair of load tracks' 396 kN/m, and whether their leads can hold it steadily is round 9's answer: no.)")
    return rows


def main():
    os.makedirs(LOGS, exist_ok=True)
    sys.stdout = Tee(os.path.join(LOGS, "round10_tailored.txt"))
    print("=" * 150 + "\n  ROUND 10, ITEM 2: THE TAILORED WINDING FOR A LOCAL LOAD (bearing_round10_tailored.py)\n" + "=" * 150)
    c = curves()
    print(f"  One support track: {W_T/1e3:.1f} kN/m of casing everywhere and {F_LOAD/1e3:.0f} kN/m more under the load ({T_PER_M:.0f} tonnes per meter of ring on the two support tracks).")
    out = {}
    out["part1"] = p1 = part1(c)
    out["part2"] = part2(c, p1)
    out["part3"] = part3(c, p1)
    if "--pass" in sys.argv:
        out["part4"] = part4(out["part3"])
    json.dump(out, open(os.path.join(LOGS, "round10_tailored.json"), "w", encoding="utf-8"), indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else (o.tolist() if isinstance(o, np.ndarray) else float(o)))
    print("\n  DONE")


if __name__ == "__main__":
    main()
