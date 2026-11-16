#!/usr/bin/env python3
"""
Round 10 (2026-11-17), item 1: the gap regulator in the compensation coils, with the roll carried.
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Brief: Vol_III/Claude_Code_instructions_round10_2026-11-17.md.

    python bearing_round10.py --part=law        the law and its gains
    python bearing_round10.py --part=table      the launch table (1, 3, 5, 10, 20, 30 km/s), point and cross-section models
    python bearing_round10.py --part=cases      the 114 cases with the regulator in every lane that had the protection
    python bearing_round10.py --part=roll       the cases that can roll, in the cross-section model (gap_regulator.section_run)
    python bearing_round10.py --part=cost       what the regulator adds to the coils' currents and heat

Writes runs/round10_logs/round10_*.json; the tables go to stdout (redirect to runs/round10_logs/round10_*.txt).
"""
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_round7 as r7
import final_bearing as fb
import gap_regulator as gr

LOGS = os.path.join(r7.OUT, "round10_logs")
SPEEDS = (1, 3, 5, 10, 20, 30)
MU0 = 4e-7 * math.pi


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    json.dump(obj, open(os.path.join(LOGS, name), "w", encoding="utf-8"), indent=1, default=float)


def load(name):
    return json.load(open(os.path.join(LOGS, name), encoding="utf-8"))


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def setup():
    """Round 9's lanes (bearing_round9_parts._setup) with the regulator in place of the protection, and the cross-section's layout."""
    import bearing_round9_parts as p9
    import bearing_round7_parts as p7
    pz, sp, ct, tabs, lane9, c = p9._setup()
    fast = r7.load("round7_fast_track.json")
    k = fast["F70"] / fast["F_req"]                                   # the fast track's rating over its full pull
    lift1, load1 = c["lift1"], c["load1"]

    def reg_pair(err=0.0, n_lift=2, n_load=2):
        return pz.regulator(pz.MU, 2 * lift1, k * n_load * load1, k * n_lift * lift1, sensor_err=err)

    def lane(u, **kw):
        kw.setdefault("reg", reg_pair())
        return lane9(u, False, **kw)

    def lay(halves=False, split_sup=True):
        return gr.layout(p7.X_C, p7.X_S, lift1, load1, k * lift1, k * load1, pz.W, halves=halves, split_sup=split_sup)

    reg_side = pz.regulator(pz.MU / 2, lift1, k * load1, k * lift1)
    return dict(pz=pz, sp=sp, ct=ct, tabs=tabs, lane=lane, lane9=lane9, c=c, k=k, lay=lay, reg_side=reg_side, reg_pair=reg_pair, X_C=p7.X_C, X_S=p7.X_S, fast=fast)


def lose(fast, kind, side=-1, which=None, t=None):
    """Mark the track(s) of a kind on a side as lost (from the start unless t is given).  which: index among that side's tracks of the kind (halves)."""
    out = [dict(f) for f in fast]; n = 0
    for f in out:
        if f["kind"] == kind and (f["x"] > 0) == (side > 0):
            if which is None or n == which:
                f["lost_t"] = "start" if t is None else t
            n += 1
    return out


def section(S, u, fast=None, sup=None, plan="kept", side_err=(0.0, 0.0), reg_off=False, t_settle=3.0, keep_hist=False, **lane_kw):
    """One pass in the cross-section model: the point model's lane gives the sled's load, the expected load and the planned command."""
    pz = S["pz"]
    ln = S["lane"](u, **lane_kw)
    res, h = pz.point_pulsed([ln], u, S["tabs"], t_settle=t_settle)
    f0, s0 = S["lay"]()
    reg = dict(S["reg_side"])
    if reg_off:
        reg.update(om=0.0, ki=0.0)
    out = gr.section_run(h["t"], h["q"][0], h["q_pred"][0], h["F_cmd"][0], fast or f0, sup or s0, reg, tau=ln["tau"], side_err=side_err, couple=S["c"]["cp"], plan=plan,
                         keep_hist=keep_hist)
    return out, res[0], h


def part_law():
    hr("THE REGULATOR: THE LAW AND ITS GAINS")
    S = setup(); pz = S["pz"]; rp = S["reg_pair"](); rs = S["reg_side"]
    print(f"  On each side of the cable, every {rs['sample']*1e3:.0f} ms (the gap sensors' rate), from the gap of that side read {rs['delay']*1e3:.0f} ms ago:")
    print(f"      e    = the side's gap error: the measured movement of the casing against the cable at that side, less the movement the plan predicts")
    print(f"      F    = - m (om^2 e + 2 om de/dt) - k_i integral(e)        de/dt from two readings {rs['window']*1e3:.0f} ms apart")
    print(f"  F is added to that side's share of the planned pull.  The sum goes to the load track of the side if it is up and to the lift track if it is down.")
    print(f"  The coils' current (the root of the pull) moves no faster than zero to full pull ({S['c']['lift1']/1e3:.1f} kN/m per track) in {rs['rise']*1e3:.0f} ms, the supply's voltage,")
    print(f"  and no track pulls more than its rating ({S['k']*S['c']['lift1']/1e3:.1f} kN/m, {S['k']:.3f} times its full pull).  There is no dead band.")
    print(f"  The predicted movement: two free masses (casing {pz.MU:,.0f} kg/m, cable {pz.MU_C:,.0f} kg/m) under the planned pull as the coils deliver it and the sled load")
    print(f"  the plan expects, led back to zero a decade slower than the loop (om / 10).  The planned pull moves the casing on purpose before the spacecraft arrives,")
    print(f"  so a regulator that held the gap at 100 mm would fight the plan (tried: it moved the casing 39 mm at 30 km/s, against {12.6:.1f} mm with the reference).")
    print(f"  The gains, by rule: critically damped; om set for {pz.REG_MARGIN_DEG:.0f} degrees of phase margin against the loop's own dead time T_d = delay + sample/2 + window/2")
    print(f"  = {rs['T_d']*1e3:.1f} ms (a mass under this law crosses over at 2.058 om with 76.3 degrees of lead, so om = (76.3 - {pz.REG_MARGIN_DEG:.0f}) degrees / (2.058 T_d));")
    print(f"  the integral a decade slower, k_i = m om^3 / 10; m the mass per meter one side moves, half the casing's.")
    print(f"      om = {rs['om']:.1f} rad/s ({rs['hz']:.2f} Hz)     per side: m = {rs['mass']:,.0f} kg/m, m om^2 = {rs['mass']*rs['om']**2/1e6:.2f} MN/m per m (kN/m per mm), 2 m om = {2*rs['mass']*rs['om']/1e3:.0f} kN s/m per m,"
          f" k_i = {rs['ki']/1e6:.0f} MN/m s per m")
    print(f"  Against it the support at constant flux pushes the wrong way with 0.35 kN/m per mm per track (round 9), a fiftieth of the loop's stiffness per side.")
    print(f"  In the point and beam models the two sides act as one: m = {rp['mass']:,.0f} kg/m, the same om.  The beam model has no sample and hold (it acts every step).")
    save("round10_law.json", dict(side=rs, pair=rp, rating=S["k"]))


def part_table():
    hr("THE LAUNCH TABLE WITH THE REGULATOR (the 10 km sled, 5,000 t on 100 m, the pair of tracks)")
    S = setup(); pz = S["pz"]
    print(f"  {'u km/s':>7} {'planned pull, pair kN/m':>23} {'track':>6} {'even field T':>12} {'plan trim':>9} {'largest pull, pair kN/m':>23} {'of it the regulator':>19} {'casing moves mm':>15} "
          f"{'the plan predicts mm':>20} {'off the plan mm':>15} {'roll mrad':>9} {'section: moves mm':>17} {'gap left mm':>11} {'round 9 mm':>10}")
    r9 = {1: 12.14, 3: 11.05, 5: 11.15, 10: 11.14, 20: 11.59, 30: 13.90}
    rows = []
    for uk in SPEEDS:
        u = uk * 1e3
        ln = S["lane"](u, label="pair")
        res, h = pz.point_pulsed([ln], u, S["tabs"], t_settle=3.0)
        r = res[0]; t = h["t"]
        rel = h["w"][0] - h["wc"][0]
        planned = float(np.abs(h["F_ff"][0]).max()); kind = "lift" if h["F_ff"][0].min() < -h["F_ff"][0].max() else "load"
        peak = max(r["pk_hd"], r["pk_hu"]); B = math.sqrt(2 * MU0 * (planned / 2) / 1.0)
        err = float(np.abs(rel - h["z_ref"][0]).max()); ref = float(np.abs(h["z_ref"][0]).max())
        sec, _, _ = section(S, u, label="pair")
        trim = ln.get("_plan_trim_down" if kind == "lift" else "_plan_trim_up", 1.0)
        rows.append(dict(u=u, planned=planned, kind=kind, B=B, trim=trim, peak=peak, pk_fb=r["pk_fb"], move=r["max_exc_all"], ref=ref, err=err, roll=sec["roll"], move_sec=sec["move_edge"],
                         gap=min(r["min_gap_b"], r["min_gap_t"]), contact=r["contact"], J_abs=r["J_abs"], J_ff=r["J_ff"], sec_pk=sec["pk"], sec_names=sec["names"]))
        print(f"  {uk:7d} {planned/1e3:23.1f} {kind:>6} {B:12.3f} {trim:9.3f} {peak/1e3:23.1f} {r['pk_fb']/1e3:19.1f} {r['max_exc_all']*1e3:15.2f} {ref*1e3:20.2f} {err*1e3:15.2f} {sec['roll']*1e3:9.3f} "
              f"{sec['move_edge']*1e3:17.2f} {min(r['min_gap_b'], r['min_gap_t'])*1e3:11.1f} {r9[uk]:10.2f}")
    print(f"  (planned pull: the largest planned pull the coils deliver, after the trim.  Even field: sqrt(2 mu0 F / w), F half the pair's planned pull, w = 1.0 m.  Section: the same")
    print(f"  pass in the cross-section model with the two sides regulated separately; a pass with nothing lost has no roll in it.)")
    save("round10_table.json", dict(rows=rows))


def roll_for_cases(S):
    """A stand-in for compensation_tracks.roll_model in the 114 cases: the cross-section model with one lift track lost from the start and the
    plan moved to the track that is left (the rows are 'one lift track alone')."""
    def f(hist, i_lane, u, with_prot=True, support_lost=False, cap=None, **kw):
        fast, sup = S["lay"]()
        fast = lose(fast, "lift")
        if support_lost:
            sup = [dict(p, lost_t=-1e9) for p in sup]
        reg = dict(S["reg_side"])
        if not with_prot:
            reg.update(om=0.0, ki=0.0)
        o = gr.section_run(hist["t"], hist["q"][i_lane], hist["q_pred"][i_lane], hist["F_cmd"][i_lane], fast, sup, reg, tau=10e-3, couple=S["c"]["cp"], plan="moved")
        return dict(peak_roll=o["roll"], worst_edge=o["move_edge"], worst_cable_edge=o["move_edge"], min_gap_b=max(o["gap_b"], 0.0), min_gap_t=max(o["gap_t"], 0.0), peak_M=0.0, t_peak=0.0,
                    contact=o["contact"], t_end=o["t_contact"])
    return f


def part_cases():
    """Round 7's 114 cases (as round 9 ran them: the support at constant flux in every model) with the regulator in every lane that had the 10 mm
    protection, and the roll rows in the cross-section model."""
    import bearing_layout as bl
    import bearing_round6a as a6
    S = setup(); cp = S["c"]["cp"]
    A2, theva, sup, lift, tabs = a6.casing_setup()
    t0 = time.time()
    gaps, records, txt = bl.run_cases(a6.TAU, True, A2, theva, sup, lift, tabs, os.path.join("round10_logs", "round10_cases.txt"), prot_all=True,
                                      head=f"[round 10: the continuous regulator in place of the 10 mm protection, ratings {S['k']:.3f} times the full pull, the support at constant flux, "
                                           f"coupling {cp*100:.2f} percent; the roll rows in gap_regulator.section_run]", couple=cp, regulator=dict(rating=S["k"]), roll_model=roll_for_cases(S))
    rows = [dict(label=g[0], model=str(g[1]), gap_b=float(g[2]), gap_t=float(g[3]), contact=bool(g[4]), single=bool(bl.is_single(g[0]))) for g in gaps]
    recs = [dict(label=r[0], u=float(r[1]), **{k: float(v) for k, v in r[4].items()}) for r in records]
    old = json.load(open(os.path.join(r7.OUT, "round9_logs", "round9_cases_flux.json"), encoding="utf-8"))
    O = {r["label"]: r for r in old["rows"]}
    mg = lambda r: min(r["gap_b"], r["gap_t"]) * 1e3
    hr("THE 114 CASES WITH THE REGULATOR (smallest gap in mm; round 9 is the 10 mm protection)")
    print(f"  {'case':>76} {'model':>6} {'round 9':>8} {'round 10':>8} {'movement mm':>11} {'single failure':>14}")
    for r in rows:
        o = O.get(r["label"])
        print(f"  {r['label'][-76:]:>76} {r['model']:>6} {mg(o) if o else float('nan'):8.1f} {mg(r):8.1f} {100 - mg(r):11.1f} {'yes' if r['single'] else 'no':>14}" + ("   CONTACT" if r["contact"] else ""))
    singles = [r for r in rows if r["single"]]
    worst = sorted(singles, key=mg)[:8]
    print(f"\n  {len(rows)} cases.  The smallest gaps in any single failure: " + "; ".join(f"{r['label']} ({r['model']}) {mg(r):.1f} mm" for r in worst))
    print(f"  Contacts among the single failures: {[r['label'] for r in singles if r['contact']] or 'none'}.")
    print(f"  All contacts: {[r['label'] for r in rows if r['contact']]}")
    tol = {}
    for r in rows:
        if r["label"].startswith("tolerance "):
            key = r["label"][10:].rsplit(" at ", 1)[0]; tol.setdefault(key, []).append((100 - mg(r), r["label"].rsplit(" at ", 1)[1]))
    kinds = [("timing 5 ms", ("timing -5 ms", "timing +5 ms")), ("mass 10 percent", ("mass -10 %", "mass +10 %")), ("speed 5 percent", ("speed -5 %", "speed +5 %")),
             ("gap measurement 5 mm", ("sensor -5 mm", "sensor +5 mm")), ("the rise halved", ("rise -50 %",)), ("the rise half again as long", ("rise +50 %",))]
    print(f"\n  The tolerances (the pair), the largest movement in mm over 1, 3, 10 and 30 km/s (the movement includes the plan's own, see the launch table):")
    tout = []
    for name, keys in kinds:
        allv = [(m, f"{k} at {at}") for k in keys for m, at in tol.get(k, [])]
        best = max(allv, default=(float("nan"), ""))
        per = "; ".join(f"{k}: " + ", ".join(f"{m:.1f}" for m, at in tol.get(k, [])) for k in keys)
        tout.append(dict(name=name, mm=best[0], case=best[1])); print(f"    {name:>28}: {best[0]:6.2f} mm ({best[1]})    [{per}]")
    save("round10_cases.json", dict(rows=rows, records=recs, tolerances=tout, worst=[dict(label=r["label"], model=r["model"], gap=mg(r)) for r in worst], seconds=time.time() - t0))
    print(f"\n  ({time.time()-t0:.0f} s)")


def main():
    part = "law"
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a in ("-h", "--help"): print(__doc__); return
    t0 = time.time()
    if part in ("law", "table", "cases"):
        globals()["part_" + part]()
    else:
        import bearing_round10_parts as p10
        getattr(p10, "part_" + part)()
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
