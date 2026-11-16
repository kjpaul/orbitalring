#!/usr/bin/env python3
"""
Round 10 (2026-11-17), item 1: the cases that can roll, and what the regulator costs (run through bearing_round10.py --part=roll|cost).
"""
import math

import numpy as np

import bearing_round10 as r10
import gap_regulator as gr
from bearing_round10 import hr, save, setup, section, lose


def show(lab, o, extra=""):
    pk = max(o["pk"]) / 1e3
    print(f"  {lab:<92} {('CONTACT at %.0f ms' % (o['t_contact']*1e3)) if o['contact'] else 'holds':>18} {o['gap']*1e3:8.1f} {o['move_c']*1e3:9.1f} {o['roll']*1e3:9.2f} {pk:10.1f} {max(o['pk_reg'])/1e3:10.1f}{extra}")
    return dict(label=lab, contact=o["contact"], t_contact=o["t_contact"], gap=o["gap"], move_c=o["move_c"], roll=o["roll"], pk=o["pk"], names=o["names"], heat=o["heat"], pk_reg=o["pk_reg"])


def head(t):
    hr(t)
    print(f"  {'case':<92} {'':>18} {'gap mm':>8} {'center mm':>9} {'roll mrad':>9} {'peak kN/m':>10} {'regulator':>10}")
    print(f"  (gap: the smallest gap at either edge of the cable's face; center: the casing's largest movement on the centerline; peak: the largest pull of any one track or half;")
    print(f"  regulator: the largest pull one side's regulator asks for)")


def part_roll():
    S = setup(); pz = S["pz"]; rows = []
    f0, s0 = S["lay"](); fh, _ = S["lay"](halves=True)
    I_cas = pz.MU * gr.R_GYR**2
    hr("WHY A LOST FAST TRACK CANNOT BE HELD BY ANY REGULATOR IN THIS CROSS-SECTION (an argument, then the runs)")
    print(f"  With one lift track gone, the side that has lost it has nothing that can pull the casing down.  Call the other side A (it has a lift and a load track, so any pull up or")
    print(f"  down at +{S['X_C']:.2f} m) and the side of the loss B (a load track only, so only a pull up at -{S['X_C']:.2f} m).  Above the neutral speed the sled lifts the casing on the centerline.")
    print(f"  For the casing to end the pass at rest, the pulls must take out the sled's whole impulse J with no turning impulse left: the down impulse on A and on B must each be J / 2.")
    print(f"  B cannot pull down at all.  Whatever A does alone leaves either the heave or the roll moving, and the support at constant flux does not bring either back (it pushes the")
    print(f"  wrong way).  The same holds below the neutral speed with one load track gone.  This is the momentum of two rigid bodies per meter; it does not depend on the law or its gains.")
    print(f"  What is not in it: the casing's torsional stiffness along the track (no script has it), and the support coils, which are frozen through an event.")
    head("ONE LIFT TRACK LOST, THE ROLL CARRIED (the cross-section model; the regulator on each side; three ways of handling the plan)")
    for uk in (10, 15, 20, 25, 27, 30):
        for plan, txt in (("kept", "the track left keeps its own half of the plan"), ("moved", "the whole plan to the track left, up to its rating"), ("dropped", "no planned pull")):
            o, _, _ = section(S, uk * 1e3, fast=lose(f0, "lift"), plan=plan, label="x")
            rows.append(dict(show(f"one lift track lost at {uk} km/s: {txt}", o), group="lift lost", u=uk, plan=plan))
    head("ONE LOAD TRACK LOST AT 1 km/s")
    for plan, txt in (("kept", "the track left keeps its own half of the plan"), ("moved", "the whole plan to the track left, up to its rating"), ("dropped", "no planned pull")):
        o, _, _ = section(S, 1e3, fast=lose(f0, "load"), plan=plan, label="x")
        rows.append(dict(show(f"one load track lost at 1 km/s: {txt}", o), group="load lost", u=1, plan=plan))
    head("THE REMEDY: EACH FAST TRACK AS TWO INDEPENDENT HALVES (half the pull and half the rating each, 0.25 m either side of the track's center [the halves are not solved])")
    for uk in SPEEDS_R:
        o, _, _ = section(S, uk * 1e3, fast=fh, label="x")
        rows.append(dict(show(f"halves, nothing lost, {uk} km/s", o), group="halves normal", u=uk, plan="kept"))
    for uk in (10, 15, 20, 25, 27, 30):
        for plan, txt in (("moved", "the lost half's share to the three halves left, up to their ratings"), ("balanced", "each side gives the same, what the weaker side can"), ("kept", "nothing reassigned")):
            o, _, _ = section(S, uk * 1e3, fast=lose(fh, "lift", which=0), plan=plan, label="x")
            rows.append(dict(show(f"one half of a lift track lost at {uk} km/s: {txt}", o), group="lift half lost", u=uk, plan=plan))
    for uk in (1, 3):
        for plan, txt in (("moved", "the lost half's share to the three halves left, up to their ratings"), ("balanced", "each side gives the same, what the weaker side can"), ("kept", "nothing reassigned")):
            o, _, _ = section(S, uk * 1e3, fast=lose(fh, "load", which=0), plan=plan, label="x")
            rows.append(dict(show(f"one half of a load track lost at {uk} km/s: {txt}", o), group="load half lost", u=uk, plan=plan))
    o, _, _ = section(S, 30e3, fast=lose(fh, "lift", which=0, t=0.5 * 100.0 / 30e3), plan="kept", label="x")
    rows.append(dict(show("one half of a lift track lost in the middle of the pull at 30 km/s, nothing reassigned", o), group="lift half lost", u=30, plan="kept mid"))
    head("OTHER CASES THAT CAN ROLL (the six tracks as designed)")
    for uk in (30, 1):
        for e in (0.005, -0.005):
            o, _, _ = section(S, uk * 1e3, side_err=(e, 0.0), label="x")
            rows.append(dict(show(f"the gap sensors of one side read {e*1e3:+.0f} mm wrong, pass at {uk} km/s", o), group="sensor", u=uk))
    o, _, _ = section(S, 30e3, side_err=(0.005, 0.0), no_sled=True, cancel=False, label="x")
    rows.append(dict(show("the gap sensors of one side read +5 mm wrong, no sled", o), group="sensor", u=0))
    o, _, _ = section(S, 30e3, side_err=(0.005, -0.005), label="x")
    rows.append(dict(show("the two sides' sensors read 5 mm wrong in opposite senses, pass at 30 km/s", o), group="sensor", u=30))
    for lab, kw in (("as the zone arrives at 30 km/s", dict()), ("at rest, no sled", dict(no_sled=True, cancel=False))):
        sup = [dict(p) for p in s0]; sup[0]["lost_t"] = 0.0               # the outer loop of the + side
        o, _, _ = section(S, 30e3, sup=sup, label="x", **kw)
        rows.append(dict(show(f"one loop of a support track lost {lab} (the fast tracks hold; the other loops are not trimmed in the run)", o), group="support loop", u=30))
    sup = [dict(p) for p in s0]; sup[0]["lost_t"] = 0.0; sup[1]["lost_t"] = 0.0
    o, _, _ = section(S, 30e3, sup=sup, no_sled=True, cancel=False, label="x")
    rows.append(dict(show("both loops of one support track lost at rest (a double failure)", o), group="support loop", u=0))
    for lab, kw in (("the pull 50 ms late at 30 km/s", dict(dt_err=50e-3)), ("the pull 50 ms early at 30 km/s", dict(dt_err=-50e-3)), ("the pull fired with no sled at 30 km/s", dict(no_sled=True)),
                    ("the pull 10 ms late at 30 km/s", dict(dt_err=10e-3))):
        o, _, _ = section(S, 30e3, label="x", **kw)
        rows.append(dict(show(lab + " (nothing lost; no roll)", o), group="timing", u=30))
    for lab, kw in (("the pull 50 ms late at 30 km/s", dict(dt_err=50e-3)), ("the pull fired with no sled at 30 km/s", dict(no_sled=True))):
        o, _, _ = section(S, 30e3, fast=lose(fh, "lift", which=0), plan="moved", label="x", **kw)
        rows.append(dict(show(lab + ", halves, one lift half already lost (a double failure)", o), group="timing", u=30))
    save("round10_roll.json", dict(rows=rows, X_C=S["X_C"], X_S=S["X_S"], I_cas=I_cas))


SPEEDS_R = (1, 30)


def part_cost():
    hr("WHAT THE REGULATOR COSTS")
    S = setup(); pz = S["pz"]; fast = S["fast"]
    import bearing_round7_parts as p7
    import bearing_round9 as r9
    ld = p7.lead(); heat9 = r9.load("round9_heat.json") if hasattr(r9, "load") else None
    print(f"  Through a pass (the point model, the pair, nothing wrong): the pull-seconds of the planned pull alone against the pull-seconds with the regulator's part in, both as")
    print(f"  the integral of |pull| over the pass at one place.  The leads' Joule heat and the tape's pulse loss go with the current squared, which is the pull, so the ratio is")
    print(f"  the factor on a pass's lead heat at that speed.")
    print(f"  {'u km/s':>7} {'planned, kN s/m':>15} {'with the regulator':>18} {'ratio':>6} {'regulator peak, pair kN/m':>25} {'as a share of full pull':>23} {'equal current, share of full':>28}")
    rows = []
    for uk in r10.SPEEDS:
        u = uk * 1e3
        res, h = pz.point_pulsed([S["lane"](u, label="pair"), S["lane"](u, couple=None, label="no coupling")], u, S["tabs"], t_settle=3.0)
        r = res[0]; rn = res[1]
        full = 2 * S["c"]["lift1"]
        rows.append(dict(u=u, J_ff=r["J_ff"], J_abs=r["J_abs"], ratio=r["J_abs"] / r["J_ff"], pk_fb=r["pk_fb"], ratio_nocouple=rn["J_abs"] / rn["J_ff"]))
        print(f"  {uk:7d} {r['J_ff']/1e3:15.2f} {r['J_abs']/1e3:18.2f} {r['J_abs']/r['J_ff']:6.3f} {r['pk_fb']/1e3:25.1f} {r['pk_fb']/full*100:22.1f}% {math.sqrt(r['pk_fb']/full)*100:27.1f}%"
              f"   (with no load-support coupling the ratio is {rn['J_abs']/rn['J_ff']:.3f})")
    print(f"  Where the ratio is under 1 the regulator has taken pull off the plan (the plan's trim and the lag leave it a little high there).  Where it is over 1 at 1 to 5 km/s, most")
    print(f"  of it is the load tracks making up what their own coupling takes from the support beside them (round 7's 12.4 percent), which the 10 mm protection also had to give.")
    # between launches
    res, h = pz.point_pulsed([S["lane"](30e3, no_sled=True, cancel=False, label="at rest")], 30e3, S["tabs"], t_settle=3.0)
    r = res[0]
    print(f"\n  Between launches (the same lane with no sled and no pull, 5 s): the casing moves {r['max_exc_all']*1e6:.3f} micrometers and the regulator's largest pull is {r['pk_fb']:.3f} N/m.")
    print(f"  With nothing disturbing the casing the model's regulator draws nothing.  The model has no sensor noise, no vibration of the cable and no drift of the support, so it")
    print(f"  cannot give the standing regulation current.  The 576 W per site of regulation heat is round 6a's allowance (4.8 percent of full current held steadily) and stays as it is:")
    print(f"  nothing run here confirms it or moves it.  [the regulator's standing current with real sensor noise: not known]")
    save("round10_cost.json", dict(rows=rows, rest=dict(move=r["max_exc_all"], pk_fb=r["pk_fb"])))


def part_launch():
    """The regulator's factor on a pass (part_cost) carried along the heaviest launch: the leads' full-current seconds and the cycled energy."""
    import bearing_round5 as r5
    import bearing_round4 as r4
    import sled_passage_pulsed as pz
    hr("THE REGULATOR'S PART IN THE HEAT OF A LAUNCH")
    cost = r10.load("round10_cost.json")["rows"]; fast = r7_load()
    us = np.array([r["u"] for r in cost]); k = np.array([r["ratio"] for r in cost])
    t_p, v_p, P_p = pz.lsm_profile()
    a = b = c = d = 0.0
    for i in range(1, len(t_p)):
        u = float(v_p[i]); dt = float(t_p[i] - t_p[i - 1]); dx = u * dt
        if u <= 0:
            continue
        F, kind = r4.pull_track_of_u(u); f = float(np.interp(u, us, k))
        ls = (2 * r5.N_FAM / r5.L_FAST) * (F / fast["F_req"]) * r5.hold_time(u)[0] * dx
        a += ls; b += ls * f; c += (F / fast["F_req"]) * dx; d += (F / fast["F_req"]) * dx * f
    print(f"  The factor of part_cost (1.24 at 1 km/s, 1.30 at 3, 1.49 at 5, 1.00 from 10 km/s up, straight lines between, the 1 km/s value below it) weighted along the launch:")
    print(f"  the leads' full-current seconds {a:.4e} planned, {b:.4e} with the regulator, {b/a:.4f} times; the cycled pull-meters {d/c:.4f} times.")
    E_lead = 21.1e9; E_tape = 587e9
    print(f"  On round 9's figures: the leads' Joule heat {E_lead/1e9:.1f} GJ becomes {E_lead*b/a/1e9:.1f} GJ; if the tape's pulse loss went with the same pull-meters (it is a loss per cycle of")
    print(f"  the planned pull, and the regulator's part is a change of amplitude, not more cycles [an upper figure]) {E_tape/1e9:.0f} GJ would become {E_tape*d/c/1e9:.0f} GJ.")
    print(f"  The cold heat of a launch: 608 GJ becomes {(E_lead*b/a + E_tape)/1e9:.0f} GJ with the leads alone, at most {(E_lead*b/a + E_tape*d/c)/1e9:.0f} GJ.")
    save("round10_launch.json", dict(lead_factor=b / a, cycle_factor=d / c, E_lead=E_lead * b / a, E_cold_lo=E_lead * b / a + E_tape, E_cold_hi=E_lead * b / a + E_tape * d / c))


def r7_load():
    import bearing_round7 as r7
    return r7.load("round7_fast_track.json")
