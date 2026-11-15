#!/usr/bin/env python3
"""
The bearing, round 7: the parts that follow the solves (run through bearing_round7.py --part=support|dyn|heat|fast|dump).
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Brief: Vol_III/Claude_Code_instructions_bearing_round7_2026-11-15.md.
"""
import json
import math
import os

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_layout as bl
import bearing_optimize as bo
import bearing_round4 as r4
import bearing_round5 as r5
import bearing_round6b as r6
import final_bearing as fb
import levitation_config as lc
import lim_hysteresis as lh
from bearing_round7 import (LOGS, OUT, MU0, W_HALF, W_CASING, F_FULL, K_FLUX, SUP, FAST, E_LAUNCH, COP, hr, save, load)

SIDE = ring.CABLE_SIDE
X_C = SIDE / 2 - FAST["Wp"] / 2                                  # the load and lift tracks' centers
X_S = SIDE / 2 - FAST["Wp"] - 0.2 - SUP["Wp"] / 2                # the two support tracks' centers at the 0.2 m spacing
def half_solves():
    """The support track's solves by the gap between its center bundles: 5 mm cells, with the 2.5 mm solves in their place where they were run."""
    half = dict(load("round7_half.json") or {})
    for name, key in (("round7_half_2p5_gc0.json", "0.000"), ("round7_half_2p5.json", "0.050")):
        h = load(name)
        if h:
            half[key] = h
    return half


_h0 = half_solves().get("0.000")
if _h0:                                                          # the present track from its own 5 mm solve: the normal current and the inductance solved, not scaled
    _nt = _h0["cap"]["tapes_leg"] / 100.0
    I_REF = _h0["normal"]["NI"] / (10 * _nt)                     # A per coil at the normal 47.5 kN/m (rounds 5 to 6b scaled 455 A down from the partly saturated cap)
    L_COIL = 2 * 10 * _nt**2 * 1000.0 * (_h0["L11"] + _h0["L21"])  # H per coil, ten together (rounds 6a and 6b: 46.2 H scaled from round 4, 38.0 H from the stored energy)
else:
    I_REF, L_COIL = 455.0, 46.2
V_CH, V_HI = 210.0, 1000.0
T_C = 1.0                                                        # s: the strips heated open (0.7 s, Stage B) and the contacts closed [the contacts' time is in no source]


def lead():
    rs, rl, fin = r6.load_final()
    fbk = r5.fast_basics(rl, 10); I_full = fbk["I"]
    A = r5.lead_rows(rl)[0][10]["A_m"]; R = be.RHO_CU_77 * 1.0 / A; P_full = I_full**2 * R
    E20 = 20 * 400.0 * 8960.0 * A
    return dict(A=A, P_full=P_full, E20=E20, budget20=E20 / P_full, budget10=E20 / P_full / 4, P_steady=20 * 8 * 400 * A)


def fast_pulls(forces):
    """The load tracks' pulls that close the heave and the roll for the support forces [(x, F)]: P at +X_C, Q at -X_C; a negative value is the lift
    track of that side pulling down instead."""
    S = sum(F for x, F in forces); M = sum(x * F for x, F in forces)
    tot = W_CASING - S                                           # P + Q
    dif = -M / X_C                                               # P - Q
    return (tot + dif) / 2, (tot - dif) / 2


def changeover(I_start, I_end, force_fn, V, L, tc=T_C, dt=0.01):
    """Currents held at I_start until tc, then each driven linearly to I_end at V / L.  Returns the two fast tracks' integrals of |pull| / full,
    the time the last current arrives, the pulls at the start and at the end."""
    I_start = np.array(I_start, float); I_end = np.array(I_end, float)
    T = L * np.abs(I_end - I_start) / V; T_all = float(T.max()) if len(T) else 0.0
    t = 0.0; iP = iQ = 0.0; first = None
    while t <= tc + T_all + dt:
        s = np.clip((t - tc) / np.maximum(T, 1e-9), 0.0, 1.0) if t > tc else np.zeros_like(T)
        I = I_start + (I_end - I_start) * s
        P, Q = fast_pulls(force_fn(I))
        if first is None:
            first = (P, Q)
        iP += abs(P) / F_FULL * dt; iQ += abs(Q) / F_FULL * dt
        t += dt
    Pe, Qe = fast_pulls(force_fn(I_end))
    # the part before tc and the part of the ramp, apart: the delay allowed for a budget b is (b - ramp) / hold
    hold = max(abs(first[0]), abs(first[1])) / F_FULL
    ramp = max(iP, iQ) - hold * tc
    return dict(iP=iP, iQ=iQ, worst=max(iP, iQ), T=T_all, first=first, last=(Pe, Qe), hold=hold, ramp=ramp)


def full_tracks(xs):
    def fn(I):
        return [(x, W_HALF * (i / I_REF)**2) for x, i in zip(xs, I)]
    return fn


def I_of(F):
    return I_REF * math.sqrt(max(F, 0.0) / W_HALF)


def half_model(hs):
    """A track of two separate loops: the side's force as a quadratic form of its two loops' currents, from the two solved states (both loops at the
    normal ampere-turns, one loop alone carrying the track's normal force): F = ks (Ii^2 + Io^2) + 2 km Ii Io.  Currents in units of the normal one."""
    r = hs["alone"]["NI"] / hs["normal"]["NI"]
    ks = W_HALF / r**2; km = (W_HALF - 2 * ks) / 2
    return ks, km, r


def part_support():
    hr("S1  THE LOST SUPPORT TRACK, BY GEOMETRY")
    ld = lead()
    print(f"  The casing's weight {W_CASING/1e3:.1f} kN/m.  Support tracks at +-{X_S:.3f} m, load and lift tracks at +-{X_C:.3f} m (the {SIDE:.2f} m side).  A full support track gives {W_HALF/1e3:.1f} kN/m at {I_REF:.0f} A per coil,")
    print(f"  its force as the current squared; {L_COIL:.1f} H per coil with the ten coils of a bundle together.  A change-over: the lost part's force is gone at once; the other supplies")
    print(f"  start after {T_C:.1f} s (the strips heated open, the contacts closed) and drive their currents linearly at V / L, up or down, to the settled values; the load tracks (or, where")
    print(f"  the sign asks for it, the lift track of that side) close the heave and the roll meanwhile.  The fast leads' budget with ten families connected: the integral of pull / full")
    print(f"  over time {ld['budget10']:.2f} s ({ld['budget20']:.2f} s with all 20).  A supply that can drive its coil down as fast as up is assumed; with a dump resistor instead the fall is slower.")
    rows = []

    def report(name, loss, I0, I1, fn, L=L_COIL, note=""):
        out = {}
        for V in (V_CH, V_HI):
            c = changeover(I0, I1, fn, V, L)
            out[V] = c
        c1 = out[V_HI]; c2 = out[V_CH]
        Pe, Qe = c1["last"]
        ok_settled = abs(Pe) < 50.0 and abs(Qe) < 50.0
        tc_ok = (ld["budget10"] - c1["ramp"]) / c1["hold"] if c1["hold"] > 0 else float("inf")
        V_min = None
        if ok_settled and c1["hold"] > 0:
            lo, hi = 50.0, 20000.0
            if changeover(I0, I1, fn, hi, L)["worst"] <= ld["budget10"]:
                for _ in range(30):
                    mid = math.sqrt(lo * hi)
                    if changeover(I0, I1, fn, mid, L, dt=0.02)["worst"] <= ld["budget10"]: hi = mid
                    else: lo = mid
                V_min = hi
        elif ok_settled:
            V_min = 0.0
        rows.append(dict(name=name, loss=loss, first=c1["first"], last=c1["last"], T_ch=c2["T"], T_hi=c1["T"], int_ch=c2["worst"], int_hi=c1["worst"], ramp_hi=c1["ramp"], hold=c1["hold"], tc_ok=tc_ok,
                         settled=ok_settled, I0=list(map(float, I0)), I1=list(map(float, I1)), note=note, V_min=V_min))
        print(f"  {name:>44} | {loss:>22} | {c1['first'][0]/1e3:6.1f} {c1['first'][1]/1e3:6.1f} | {Pe/1e3:6.1f} {Qe/1e3:6.1f} | " + " ".join(f"{a:4.0f}>{b:4.0f}" for a, b in zip(I0, I1)).ljust(40)
              + f" | {c2['T']:6.1f} {c2['worst']:6.2f} | {c1['T']:6.1f} {c1['worst']:6.2f} | {('none' if tc_ok <= 0 else ('any' if tc_ok > 1e3 else '%.1f' % tc_ok)):>6} | {'yes' if ok_settled else 'NO':>3} | {('-' if V_min is None else f'{V_min:,.0f}'):>6}")

    print()
    print(f"  {'arrangement':>44} | {'the single loss':>22} | {'fast at the loss':>13} | {'settled':>13} | {'currents A, each remaining part, from > to':<40} | {'210 V: s, int':>13} | {'1 kV: s, int':>13} | {'t_c ok':>6} | steady with no fast coil | least supply V for 0.64 s")
    print(f"  {'':>44} | {'':>22} | {'lost side, other kN/m':>13} | {'kN/m':>13} | {'':<40} | {'':>13} | {'':>13} | {'s':>6} |")
    # 0. the present two tracks (round 6b's 2.3 B)
    S1 = W_CASING / (1 + X_S / X_C)
    report("0. two tracks (round 6b)", "one track", [I_REF], [I_of(S1)], full_tracks([-X_S]))
    # 1. a third track on the centerline
    for lab, c in (("equal thirds", W_CASING / 3), ("the center idle", 0.0), ("the center at a tenth", 0.1 * W_CASING), ("the center carries half", W_CASING / 2), ("the center carries it all", W_CASING),
                   ("the center carries 0.6", 0.6 * W_CASING), ("the center carries 0.7", 0.7 * W_CASING)):
        a = (W_CASING - c) / 2
        report(f"1. three tracks, {lab}", "the center track", [I_of(a), I_of(a)], [I_of(W_HALF), I_of(W_HALF)], full_tracks([-X_S, X_S]))
        report(f"1. three tracks, {lab}", "an outer track", [I_of(a), I_of(c)], [0.0, I_of(W_CASING)], full_tracks([-X_S, 0.0]))
    # 2. four half-tracks
    half = half_solves()
    out_half = {}
    for key in sorted(half):
        hs = half[key]; gc = hs["gc"]
        ks, km, r = half_model(hs)
        d = (SUP["Wl"] + 0.15 + gc) / 2                                    # a loop's center from its track's center
        xi, xo = X_S - d, X_S + d
        n_t = hs["cap"]["tapes_leg"] / 100.0                               # turns of one coil in one loop (ten coils, ten layers)
        I0 = hs["normal"]["NI"] / (10 * n_t)
        L_loop = 10 * n_t**2 * 1000.0 * hs["L11"]; M_loop = 10 * n_t**2 * 1000.0 * hs["L21"]
        L_fig8 = 2 * (L_loop + M_loop)
        jump = 1 + hs["L21"] / hs["L11"]
        out_half[key] = dict(gc=gc, ks=ks, km=km, r=r, xi=xi, xo=xo, n_t=n_t, I0=I0, L_loop=L_loop, M_loop=M_loop, L_fig8=L_fig8, jump=jump, tape=hs["cap"]["tape"], NI_cap=hs["cap"]["NI"], sat_cap=hs["cap"]["sat"],
                             share_alone=hs["alone"]["share"], sat_alone=hs["alone"]["sat"], F_jump=hs["jump"]["F"])

        def side(sign, Ii, Io):                                            # currents in amperes; the force of a side at its loops' centers
            ui, uo = Ii / I0, Io / I0
            Fi = ks * ui * ui + km * ui * uo; Fo = ks * uo * uo + km * ui * uo
            return [(sign * xi, Fi), (sign * xo, Fo)]
        # 2a. independent: an outer loop is lost on the + side; its inner neighbor jumps by flux conservation; the mirror outer loop is driven to zero, the two inner loops to r
        def fn_a(I):                                                       # I = [left outer, left inner, right inner]
            return side(-1, I[1], I[0]) + side(+1, I[2], 0.0)
        report(f"2. four half-tracks, {gc*1e3:.0f} mm between", "one loop, mirror emptied", [I0, I0, I0 * jump], [0.0, I0 * r, I0 * r], fn_a, L=L_loop, note="independent")
        # 2a'. the same loss with the mirror kept: the sound side's two loops trimmed together, the lost side's inner loop raised until no fast coil pulls
        F_r = W_CASING * X_S / (xi + X_S); F_l = W_CASING - F_r
        report(f"2. four half-tracks, {gc*1e3:.0f} mm between", "one loop, others trimmed", [I0, I0, I0 * jump], [I0 * math.sqrt(F_l / W_HALF)] * 2 + [I0 * r * math.sqrt(F_r / W_HALF)], fn_a, L=L_loop, note="independent, trimmed")
        # the inner loop lost instead: the outer loop of that side carries, at its own position
        def fn_c(I):                                                       # I = [left outer, left inner, right outer]
            return side(-1, I[1], I[0]) + side(+1, 0.0, I[2])
        F_r2 = W_CASING * X_S / (xo + X_S); F_l2 = W_CASING - F_r2
        report(f"2. four half-tracks, {gc*1e3:.0f} mm between", "an inner loop, trimmed", [I0, I0, I0 * jump], [I0 * math.sqrt(F_l2 / W_HALF)] * 2 + [I0 * r * math.sqrt(F_r2 / W_HALF)], fn_c, L=L_loop, note="independent, trimmed")
        # 2b. the two loops at +-x_out in one circuit and the two at +-x_in in another: a fault in one loop takes its whole circuit, both sides at once
        def fn_b(I):                                                       # I = [the surviving circuit]
            return side(-1, I[0], 0.0) + side(+1, I[0], 0.0)
        report(f"4. two symmetric circuits, {gc*1e3:.0f} mm between", "one circuit", [I0 * jump], [I0 * r], fn_b, L=2 * L_loop, note="symmetric")
    save("round7_support.json", dict(rows=rows, half=out_half, lead=ld, X_S=X_S, X_C=X_C, L_COIL=L_COIL))
    print(f"  ('fast at the loss' and 'settled': the load track on the lost side and the other one; a negative value is that side's lift track.  't_c ok': the delay the ten connected")
    print(f"  families' leads allow at 1 kV.  int: the larger of the two fast tracks' integrals with the {T_C:.0f} s delay, to set against {ld['budget10']:.2f} s.)")
    if out_half:
        print()
        print(f"  The support track as two separate loops (2.5 mm cells for 0 and 50 mm, 5 mm for 100 mm; the plate {SUP['Wp']:.1f} m x {SUP['T']*1e3:.0f} mm, the loops {SUP['Wl']:.1f} m between leg centers, the cap 245.5 kN/m with both loops):")
        print(f"  {'between the center bundles':>26} {'cells mm':>8} {'cap kA':>7} {'sat %':>6} {'tape kg/m':>9} {'turns per loop':>14} {'normal A':>8} {'one loop alone for 47.5 kN/m: x current':>40} {'its worst tape share':>20} {'sat %':>6} {'L loop H':>8} {'M H':>6} {'figure-eight H':>14} {'jump':>6} {'force after the jump kN/m':>26}")
        for key in sorted(out_half):
            o = out_half[key]; hs = half[key]
            print(f"  {o['gc']*1e3:23.0f} mm {hs['h']*1e3:8.1f} {o['NI_cap']/1e3:7.0f} {o['sat_cap']*100:6.1f} {o['tape']:9.1f} {o['n_t']:14.1f} {o['I0']:8.0f} {o['r']:40.3f} {o['share_alone']:20.2f} {o['sat_alone']*100:6.1f} {o['L_loop']:8.1f} {o['M_loop']:6.1f} {o['L_fig8']:14.1f} {o['jump']:6.3f} {o['F_jump']/1e3:26.1f}")
    # 3. inclined faces
    print()
    y_track = SIDE / 2 + SUP["T"] + bd.GAP0
    th = math.atan2(X_S, y_track)
    print(f"  3. Inclined faces.  With the casing's center of weight taken at the cable's center [the casing's own center of mass is in no script], a support track at {X_S:.3f} m out and")
    print(f"  {y_track:.2f} m below it must pull along a line {math.degrees(th):.1f} degrees from the vertical to point at it.  Normal: two tracks at +-{math.degrees(th):.1f} degrees, each {W_HALF/math.cos(th)/1e3:.1f} kN/m, their sideways parts")
    print(f"  ({W_HALF*math.tan(th)/1e3:.1f} kN/m each) cancelling.  After a loss: the other track alone, {W_CASING/math.cos(th)/1e3:.1f} kN/m, no moment, and a sideways force of {W_CASING*math.tan(th)/1e3:.1f} kN/m on every meter that nothing in the")
    print(f"  bearing takes: the tracks act along their own normals only and the LIM model has no force across its gap (1.5).  It trades the roll for a sideways hold.")
    # 1.5 the LIMs
    print()
    print(f"  1.5  The LIMs.  lim_physics and lim_simulation model the thrust along the track and the losses; neither has a force between the stator and its plate across the gap")
    print(f"  (no normal, attraction or repulsion term in either file).  For scale only, the thrust: 1,226 N per LIM at the peak of the design run (Section 18 of the deployment doc),")
    print(f"  two LIMs per 500 m site, {2*1226/500:.1f} N per meter of ring; a force of that size acting across the gap at the cable's side, {SIDE/2 + 0.2:.2f} m out, would give {2*1226/500*(SIDE/2+0.2):.0f} N m per meter")
    a1 = W_CASING - W_CASING / (1 + X_S / X_C)
    print(f"  against the {a1*X_C/1e3:.0f} kN m per meter of the lost track: a factor of {a1*X_C/(2*1226/500*(SIDE/2+0.2)):,.0f} short.")
    # 1.6 the casing's torsional stiffness
    print()
    m0 = a1 * X_C; a = 1000.0
    print(f"  1.6  The casing's own stiffness.  A fault at one place takes the coils that pass through it: the support force falls linearly from sound at 1 km either side to zero at")
    print(f"  the fault, so the moment per meter is m0 (1 - |x| / a) with m0 = {m0/1e3:.0f} kN m/m and a = 1,000 m, {m0*a/1e6:.0f} MN m in all, half carried to each side.  With the casing a shaft of")
    print(f"  torsional stiffness GJ held against twist at the two ends of the lost length, the torque at x is m0 (x - x^2 / 2a) and the twist at the fault m0 a^2 / (3 GJ).")
    for z in (0.010, 0.050, 0.100):
        phi = z / X_C
        print(f"    twist allowed {z*1e3:.0f} mm at the load tracks ({phi*1e3:.1f} mrad): GJ >= {m0*a*a/(3*phi):.2e} N m2")
    print(f"  No script or chapter gives the casing's torsional stiffness, so only the need is stated.  The sound tracks either side would also have to take the {m0*a/2/1e6:.0f} MN m handed to them.")
    return rows


# =============================================================================
# S1  the dynamics: heave and roll of casing and cable through the loss of a support part
# =============================================================================

def roll_run(force_fn, I_of_t, n_conn, t_end, dt=2e-4, hz=5.0, delay=1e-3, tau=10e-3):
    """Two bodies per meter, each free in heave and roll (the cable with no torsional stiffness, the casing a 10 m square frame with a radius of
    gyration of 5.77 m [UNVERIFIED, as bearing_design F.3]).  The support parts pull at their positions with their commanded force plus the
    constant-flux stiffness on their own gap; a floor (load) track and a ceiling (lift) track on each side hold their side's gap with a 5 Hz loop
    (proportional, rate and reset), the current lagging by tau, capped at what the connected families can pull.  Returns the smallest gaps, the
    peak roll, the fast tracks' peak pulls and the integrals of pull / full."""
    MU, MU_C = ring.M_LOAD_M, ring.M_CABLE_M
    r_g = 10.0 / math.sqrt(3.0); I_cas = MU * r_g**2; I_cab = MU_C * (2 * SIDE**2) / 12.0
    om = 2 * math.pi * hz; k_p = (MU / 2) * om**2; k_d = (MU / 2) * 2 * om; k_i = k_p * om / 5.0
    cap = F_FULL * min(1.0, (n_conn * 1.21 / 20.0)**2)
    w = v = th = thd = wc = vc = thc = thcd = 0.0
    sides = (+X_C, -X_C); integ = [0.0, 0.0]; cur = [0.0, 0.0]; heat = [0.0, 0.0]; peak = [0.0, 0.0]
    buf = []; nd = max(1, int(round(delay / dt))); lag = 1 - math.exp(-dt / tau)
    gmin_b = gmin_t = bd.GAP0; roll_pk = 0.0; t = 0.0; n = int(t_end / dt)
    edge = SIDE / 2; rel = rth = 0.0
    for j in range(n):
        t = j * dt
        rel = w - wc; reld = v - vc; rth = th - thc; rthd = thd - thcd
        F_up = 0.0; M = 0.0
        for x, F in force_fn(I_of_t(t)):
            if F <= 0:
                continue
            z = rel + rth * x
            Fx = F * (1 + (K_FLUX / W_HALF) * z)                       # constant flux: the pull grows as its gap closes (K_FLUX per full track at its normal force)
            F_up += Fx; M += Fx * x
        buf.append((rel, rth, reld, rthd))
        if len(buf) > nd + 1:
            buf.pop(0)
        rm, rthm, rdm, rthdm = buf[0]
        for k, x in enumerate(sides):
            z = rm + rthm * x; zd = rdm + rthdm * x
            integ[k] += z * dt
            F_des = -(k_p * z + k_d * zd + k_i * integ[k])                # positive: pull the casing up (the floor track); negative: down (the ceiling track)
            if abs(F_des) >= cap:
                integ[k] -= z * dt                                        # no wind-up while the track is at its cap
            F_des = max(-cap, min(cap, F_des))
            sgn = 1.0 if F_des >= 0 else -1.0
            target = sgn * math.sqrt(abs(F_des) / F_FULL)
            cur[k] += (target - cur[k]) * lag                             # the signed current fraction of that side's active track
            Fk = (1.0 if cur[k] >= 0 else -1.0) * F_FULL * cur[k]**2
            F_up += Fk; M += Fk * x
            heat[k] += abs(Fk) / F_FULL * dt; peak[k] = max(peak[k], abs(Fk))
        a = (F_up - W_CASING) / MU; al = M / I_cas
        ac = -(F_up - W_CASING) / MU_C; alc = -M / I_cab                 # the cable feels the reactions; its own orbit carried the normal load
        v += a * dt; w += v * dt; thd += al * dt; th += thd * dt
        vc += ac * dt; wc += vc * dt; thcd += alc * dt; thc += thcd * dt
        rel = w - wc; rth = th - thc
        roll_pk = max(roll_pk, abs(rth))
        for y in (-edge, edge):
            r_y = rel + rth * y
            gmin_b = min(gmin_b, bd.GAP0 - r_y); gmin_t = min(gmin_t, bd.GAP0 + r_y)
        if gmin_b <= 0 or gmin_t <= 0:
            break
    return dict(gap_b=gmin_b, gap_t=gmin_t, roll=roll_pk, heat=heat, peak=peak, cap=cap, contact=bool(gmin_b <= 0 or gmin_t <= 0), t_end=t,
                last=[(1.0 if c >= 0 else -1.0) * F_FULL * c * c for c in cur], rel_end=rel, rth_end=rth)


def ramp(I_start, I_end, V, L, tc=T_C):
    I_start = np.array(I_start, float); I_end = np.array(I_end, float)
    T = L * np.abs(I_end - I_start) / V
    def f(t):
        if t <= tc:
            return I_start
        s = np.clip((t - tc) / np.maximum(T, 1e-9), 0.0, 1.0)
        return I_start + (I_end - I_start) * s
    return f, float(T.max()) + tc


def part_dyn():
    hr("S1  THE DYNAMICS: THE LOSS OF A SUPPORT PART WITH THE ROLL CARRIED")
    ld = lead(); sup = load("round7_support.json")
    print(f"  Casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m, both free in heave and roll.  The lost part's force goes to zero at t = 0; the others start after {T_C:.0f} s.  Each side's")
    print(f"  fast pair (a load track below, a lift track above, at +-{X_C:.3f} m) holds its own side's gap with a 5 Hz loop, a 1 ms sensor delay and the 10 ms rise.  Deployment: ten families")
    print(f"  connected (a track can pull {min(1.0, (10*1.21/20)**2)*F_FULL/1e3:.0f} kN/m at the coils' rating, the leads' budget {ld['budget10']:.2f} s); operation: all 20 ({F_FULL/1e3:.0f} kN/m, {ld['budget20']:.2f} s).  No sled in these runs.")
    print(f"  {'arrangement, the loss':>58} {'supply':>7} {'families':>8} | {'smallest gap mm: floor':>22} {'ceiling':>8} {'peak roll mrad':>14} | {'peak pull kN/m: lost side':>25} {'other':>6} | {'integral s: lost side':>21} {'other':>6} {'budget':>6} {'within':>6} | {'pull left at the end kN/m':>25}")
    rows = []

    def run(name, fn, I0, I1, L, V, n):
        f, T = ramp(I0, I1, V, L)
        r = roll_run(fn, f, n, T + 6.0)
        bud = ld["budget10"] if n == 10 else ld["budget20"]
        ok = (max(r["heat"]) <= bud) and not r["contact"]
        rows.append(dict(name=name, V=V, n=n, **{k: v for k, v in r.items()}, budget=bud, ok=ok, T=T))
        print(f"  {name:>58} {V:6.0f}V {n:8d} | {r['gap_b']*1e3:22.1f} {r['gap_t']*1e3:8.1f} {r['roll']*1e3:14.2f} | {r['peak'][0]/1e3:25.1f} {r['peak'][1]/1e3:6.1f} | {r['heat'][0]:21.2f} {r['heat'][1]:6.2f} {bud:6.2f} {'yes' if ok else 'NO':>6} |"
              f" {r['last'][0]/1e3:12.1f} {r['last'][1]/1e3:6.1f}{'  CONTACT' if r['contact'] else ''}", flush=True)

    # the case to beat: two tracks, one lost (the + side)
    S1 = W_CASING / (1 + X_S / X_C)
    for V in (V_CH, V_HI):
        for n in (10, 20):
            run("0. two tracks, one lost (round 6b's case)", full_tracks([-X_S]), [I_REF], [I_of(S1)], L_COIL, V, n)
    # three tracks, the center carrying half
    c = W_CASING / 2; a = (W_CASING - c) / 2
    for V in (V_CH, V_HI):
        for n in (10, 20):
            run("1. three tracks, center half: the center lost", full_tracks([-X_S, X_S]), [I_of(a), I_of(a)], [I_of(W_HALF), I_of(W_HALF)], L_COIL, V, n)
            run("1. three tracks, center half: an outer lost", full_tracks([-X_S, 0.0]), [I_of(a), I_of(c)], [0.0, I_of(W_CASING)], L_COIL, V, n)
    # the half-track arrangements, from the solves
    half = (sup or {}).get("half", {})
    for key in sorted(half):
        o = half[key]
        if abs(o["gc"] - 0.05) > 1e-6 and len(half) > 1:
            continue
        ks, km, I0, r, xi, xo = o["ks"], o["km"], o["I0"], o["r"], o["xi"], o["xo"]
        def side(sign, Ii, Io, ks=ks, km=km, I0=I0, xi=xi, xo=xo):
            ui, uo = Ii / I0, Io / I0
            return [(sign * xi, ks * ui * ui + km * ui * uo), (sign * xo, ks * uo * uo + km * ui * uo)]
        fn_a = lambda I: side(-1, I[1], I[0]) + side(+1, I[2], 0.0)          # the + side's outer loop lost: I = [left outer, left inner, right inner]
        fn_c = lambda I: side(-1, I[1], I[0]) + side(+1, 0.0, I[2])          # the + side's inner loop lost: I = [left outer, left inner, right outer]
        fn_b_out = lambda I: side(-1, I[0], 0.0) + side(+1, I[0], 0.0)      # the outer circuit lost, the inner one carries
        F_r = W_CASING * X_S / (xi + X_S); F_l = W_CASING - F_r
        F_r2 = W_CASING * X_S / (xo + X_S); F_l2 = W_CASING - F_r2
        tgt_a = [I0 * math.sqrt(F_l / W_HALF)] * 2 + [I0 * r * math.sqrt(F_r / W_HALF)]
        tgt_c = [I0 * math.sqrt(F_l2 / W_HALF)] * 2 + [I0 * r * math.sqrt(F_r2 / W_HALF)]
        for V in (V_CH, 300.0, V_HI):
            for n in (10, 20):
                run(f"2. four half-tracks ({o['gc']*1e3:.0f} mm): an outer loop lost", fn_a, [I0, I0, I0 * o["jump"]], tgt_a, o["L_loop"], V, n)
                run(f"2. four half-tracks ({o['gc']*1e3:.0f} mm): an inner loop lost", fn_c, [I0, I0, I0 * o["jump"]], tgt_c, o["L_loop"], V, n)
        for V in (V_CH, V_HI):
            for n in (10, 20):
                run(f"4. two symmetric circuits ({o['gc']*1e3:.0f} mm): one circuit lost", fn_b_out, [I0 * o["jump"]], [I0 * r], 2 * o["L_loop"], V, n)
        # the supplies never act: the fast coils alone, to show what the change-over is for
        run(f"2. four half-tracks ({o['gc']*1e3:.0f} mm): an outer loop lost, no supply acts", fn_a, [I0, I0, I0 * o["jump"]], [I0, I0, I0 * o["jump"]], o["L_loop"], V_HI, 10)
    save("round7_dyn.json", rows)
    return rows


# =============================================================================
# S2.3  the hot spot of a dump at the charging voltage
# =============================================================================

def part_dump():
    hr("S2.3  A SUPPORT COIL QUENCHED AND DUMPED AT THE CHARGING VOLTAGE: THE HOT SPOT")
    rs, rl, fin = r6.load_final(); sc = r6.support_coil(rs)
    print(f"  Method: bearing_layout.hot_spot_T, the adiabatic heating of a normal length of one 12 mm tape (its copper and Hastelloy in parallel, the copper's heat capacity from the")
    print(f"  NIST fit, the Hastelloy taken with copper's specific heat [UNVERIFIED, as round 4]), from 70 K, carrying the tape's share of the coil's current: a turn is ten soldered layers,")
    print(f"  so a tape carries a tenth.  The current stays at its value for the detection time and the strip's opening, then falls as exp(-t / tau), tau = L_leak / R.  No heat leaves the spot.")
    print(f"  {'case':>44} {'coil A':>7} {'tape A':>7} {'R ohm':>6} {'V':>6} {'tau ms':>7} {'held before the fall s':>22} {'hot spot K':>10} {'time to 300 K at constant current s':>36}")
    out = []
    Ll = sc["Ll_scaled"]; V_ch = sc["Lc_scaled"] * 0.01 * sc["I0"]
    cases = []
    for lab, I in (("normal current", sc["I0"]), ("the cap", sc["I_cap"])):
        R_ch = V_ch / sc["I0"]                                         # one resistor: the charging voltage at the normal current
        for held in (0.1, 0.8, 2.0):
            cases.append((f"{lab}, the charging-voltage resistor", I, R_ch, held))
        cases.append((f"{lab}, 1 kV at this current", I, 1000.0 / I, 0.8))
    for lab, I, R, held in cases:
        tau = Ll / R; It = I / 10.0
        f = lambda t, It=It, held=held, tau=tau: It if t < held else It * math.exp(-(t - held) / tau)
        T = bl.hot_spot_T(f, held + 10 * tau)
        t300 = bl.hot_spot_time(It)
        out.append(dict(case=lab, I=I, I_tape=It, R=R, V=I * R, tau=tau, held=held, T=T, t300=t300))
        print(f"  {lab:>44} {I:7.0f} {It:7.1f} {R:6.2f} {I*R:6.0f} {tau*1e3:7.0f} {held:22.1f} {T:10.1f} {t300:36.1f}")
    print(f"  The slowest dump that keeps the spot under a limit (0.8 s held before the fall):")
    for lab, I in (("normal current", sc["I0"]), ("the cap", sc["I_cap"])):
        It = I / 10.0
        for T_lim in (150.0, 300.0):
            lo, hi = 0.01, 500.0
            for _ in range(30):
                mid = math.sqrt(lo * hi)
                f = lambda t, It=It, mid=mid: It if t < 0.8 else It * math.exp(-(t - 0.8) / mid)
                if bl.hot_spot_T(f, 0.8 + 10 * mid, dt=max(1e-3, mid / 500)) > T_lim: hi = mid
                else: lo = mid
            R = Ll / lo
            out.append(dict(case=f"slowest, {lab}, {T_lim:.0f} K", I=I, tau=lo, R=R, V=I * R))
            print(f"    {lab}, {T_lim:.0f} K: tau {lo:.1f} s, a resistor of {R*1e3:.1f} mohm, {I*R:.1f} V at the start")
    save("round7_dump.json", dict(coil=sc, rows=out, V_ch=V_ch))
    return out


# =============================================================================
# the saved solve of the accepted fast track
# =============================================================================

def fast_solve():
    meta = load("round7_fast_solve.json")
    if meta is None:
        raise SystemExit("runs/round7_logs/round7_fast_solve.json missing: run 'python bearing_round7.py --part=solve-fast' first")
    z = np.load(os.path.join(LOGS, "round7_fast_fields.npz"))
    return meta, z


def level_fields(z, k):
    tag = f"{k:.2f}"
    return [(z[f"Bx{i}_{tag}"], z[f"By{i}_{tag}"]) for i in range(len(z["leg_x"]))]


# the cells turned as ruled (round 6b, Section 2.2): per leg, the upper and the lower row of ten cells; the leg at +x is the mirror
PATTERN_NEG = ("vvhhhhvvvv", "vhhhhhhhvv")                              # the leg at -300 mm: upper row, lower row, from -x to +x
PATTERN_POS = ("vvvvhhhhvv", "vvhhhhhhhv")


def cell_index(z, i, ncx=10, ncy=2):
    xs, ys = z[f"xs{i}"], z[f"ys{i}"]; x = float(z["leg_x"][i])
    ix = np.clip(((xs - (x - FAST["bw"] / 2)) / FAST["bw"] * ncx).astype(int), 0, ncx - 1)
    iy = np.clip(((ys - ys.min()) / (ys.max() - ys.min() + 1e-12) * ncy).astype(int), 0, ncy - 1)
    return ix, iy


def turned(z, fields, pattern=None):
    """The fields with the two components exchanged in the cells that lie faces vertical, so that the faces-horizontal formulas give the turned
    cells' loss and limit.  pattern None: each cell the better way, found from the fields; returns the fields and the pattern used."""
    out = []; pats = []
    for i, (Bx, By) in enumerate(fields):
        ix, iy = cell_index(z, i)
        x = float(z["leg_x"][i])
        if pattern is None:
            rows = [["h"] * 10 for _ in range(2)]
            for cx in range(10):
                for cy in range(2):
                    m = (ix == cx) & (iy == cy)
                    if (Bx[m]**2).mean() < (By[m]**2).mean():
                        rows[cy][cx] = "v"
            pat = ("".join(rows[1]), "".join(rows[0]))
        else:
            pat = pattern[0] if x < 0 else pattern[1]
        Bx2, By2 = Bx.copy(), By.copy()
        for cx in range(10):
            for cy in range(2):
                if pat[1 - cy][cx] == "v":                               # pat[0] is the upper row (cy = 1)
                    m = (ix == cx) & (iy == cy)
                    Bx2[m], By2[m] = By[m], Bx[m]
        out.append((Bx2, By2)); pats.append(pat)
    return out, pats


def part_fast():
    hr("S2.2  THE ACCEPTED FAST TRACK AS ONE SOLVE (2.5 mm cells)")
    meta, z = fast_solve(); theva, theva77 = be.load_theva()
    F_req = meta["F_req"]; NI_req = meta["NI_req"]
    lv = {round(l["k"], 2): l for l in meta["levels"]}
    print(f"  The 1.0 m x 200 mm plate, the 0.6 m loop, the 204 x 50 mm bundle (its top 10 mm under the casing's face), solved in {meta['seconds']/60:.0f} min: {NI_req/1e3:.1f} kA per leg for the")
    print(f"  requirement of {F_req/1e3:.1f} kN/m (solved {lv[1.0]['F']/1e3:.1f}); the plate {lv[1.0]['sat']*100:.1f} percent saturated, its peak field {lv[1.0]['B_plate']:.2f} T; stored energy {lv[1.0]['W']/1e3:.1f} kJ per meter of track.")
    out = dict(NI_req=NI_req, F_req=F_req)
    ks = (1.0, 1.25, 1.5, 2.0)
    NIs = np.array([lv[k]["NI"] for k in ks]); Fs = np.array([lv[k]["F"] for k in ks])
    li = dict(r6.launch_profile(lv[1.0]["W"], F_req), F_req=F_req)
    Phi = 2 * lv[1.0]["W"] / NI_req; V_turn = Phi * 25.0 / 10e-3; V = 10 * V_turn; I = NI_req / 20 / 10
    f_req = level_fields(z, 1.0)
    best_f, best_pat = turned(z, f_req, None)
    print(f"  {'the cells':>42} {'worst tape: T, deg from its normal, A':>38} {'tapes per leg':>13} {'layers per turn':>15} {'tape kg/m':>9} {'margin':>6} {'loss per launch GJ (strip + slab)':>34}")
    res = {}
    for lab, pattern in (("all faces horizontal", ("hhhhhhhhhh", "hhhhhhhhhh")), ("turned as ruled (round 6b's pattern)", "ruled"), ("each cell the better way at 2.5 mm", None)):
        pat = (PATTERN_NEG, PATTERN_POS) if pattern == "ruled" else pattern
        if isinstance(pat, tuple) and isinstance(pat[0], str):
            pat = (pat, pat)                                              # the same rows for both legs (all horizontal)
        wi = []
        for k in ks:
            fk, pats = turned(z, level_fields(z, k), pat)
            wi.append(r6.worst_ic(fk, "h", theva))
        ic, B, th = wi[0]
        tapes = NI_req / (0.5 * ic); tape = 2 * tapes * 0.003 * bd.TAPE_T * bd.TAPE_RHO
        share = (NIs / tapes) / np.array([w[0] for w in wi])
        NI70 = float(np.interp(0.7, share, NIs)); F70 = math.exp(float(np.interp(math.log(NI70), np.log(NIs), np.log(Fs))))
        NI100 = float(np.interp(1.0, share, NIs)) if share.max() >= 1.0 else float("nan")
        layers = (NI_req / 20) / (0.5 * ic) / 10
        fm, pats = turned(z, f_req, pat)
        Ep, Es = r6.tape_loss_launch(fm, "h", tapes, layers, li)
        res[lab] = dict(ic=ic, B=B, th=th, tapes=tapes, tape=tape, margin=F70 / F_req, F70=F70, layers=layers, Ep=Ep, Es=Es, pats=pats, I70=NI70 / NI_req)
        print(f"  {lab:>42} {f'{B:.2f} T, {th:.0f} deg, {ic:.1f} A':>38} {tapes:13,.0f} {layers:15.1f} {tape:9.1f} {F70/F_req:6.2f} {f'{(Ep+Es)/1e9:.0f} ({Ep/1e9:.0f} + {Es/1e9:.0f})':>34}")
    r = res["turned as ruled (round 6b's pattern)"]
    print(f"  The pattern found at 2.5 mm, per leg, upper row / lower row: " + "; ".join(f"{p[0]} / {p[1]}" for p in res['each cell the better way at 2.5 mm']['pats']) + f"; as ruled: {PATTERN_NEG[0]} / {PATTERN_NEG[1]} and its mirror.")
    layers = math.ceil(r["layers"]); t_s = layers * r5.TAPE_PITCH
    print(f"  The coil, cells turned as ruled: {I:,.0f} A per turn at the requirement ({NI_req/20/1e3:.2f} kA per coil, 10 turns), {r['layers']:.1f} layers per turn ({layers} as whole layers, a stack {t_s*1e3:.2f} mm thick),")
    print(f"  {V_turn:,.0f} V per turn and {V/1e3:.2f} kV at the terminals for the 10 ms rise; {lv[1.0]['W']*25/20/1e6:.3f} MJ per coil.  Tape per track {r['tape']:.1f} kg/m; the rating (the worst tape at 70 percent)")
    print(f"  {r['F70']/1e3:.0f} kN/m, {r['margin']:.2f} times the requirement, at {r['I70']:.2f} times the current.")
    # the packing with Section 1.3's films
    t_fam = V / r6.E_INS; t1 = max(V_turn / r6.E_INS, 0.05e-3)
    h_cell = dict(w=5 * 0.003 + 4 * t1, h=2 * t_s + 5 * t1)               # faces horizontal: 5 across, 2 up, rows left to right (a five-turn film between the rows)
    v_row = dict(w=2 * t_s + t1, h=5 * 0.003 + 4 * 2 * t1)                # faces vertical, two across and five up, rows in order: one turn across, two turns between rows
    v_col = dict(w=2 * t_s + 5 * t1, h=5 * 0.003 + 4 * t1)                # faces vertical, columns in order: one turn up a column, five turns between the columns
    print(f"  Packing, films sized for the turns they separate at {r6.E_INS/1e6:.0f} kV/mm: walls {t_fam*1e3:.2f} mm, one-turn film {t1*1e3:.2f} mm.  A cell faces horizontal {h_cell['w']*1e3:.2f} x {h_cell['h']*1e3:.2f} mm;")
    print(f"  a cell faces vertical {v_row['w']*1e3:.2f} x {v_row['h']*1e3:.2f} mm (rows in order) or {v_col['w']*1e3:.2f} x {v_col['h']*1e3:.2f} mm (columns in order).")
    pack = {}
    for nm, vc in (("rows in order", v_row), ("columns in order", v_col)):
        widths = []
        for row in PATTERN_NEG:
            nv = row.count("v"); widths.append(nv * vc["w"] + (10 - nv) * h_cell["w"] + 11 * t_fam)
        height = 2 * max(h_cell["h"], vc["h"]) + 3 * t_fam
        pack[nm] = dict(widths=widths, height=height)
        print(f"    vertical cells with {nm}: the upper row needs {widths[0]*1e3:.1f} mm and the lower {widths[1]*1e3:.1f} mm of the {FAST['bw']*1e3:.0f}; the height {height*1e3:.1f} mm of the {FAST['bh']*1e3:.0f}"
              f" ({'fits' if max(widths) <= FAST['bw'] + 1e-9 and height <= FAST['bh'] + 1e-9 else 'DOES NOT FIT'})")
    out.update(res=res, I=I, V=V, V_turn=V_turn, layers=r["layers"], layers_whole=layers, t_s=t_s, E_coil=lv[1.0]["W"] * 25 / 20, W=lv[1.0]["W"], sat=lv[1.0]["sat"], B_plate=lv[1.0]["B_plate"], pack=pack,
               t_fam=t_fam, t1=t1, li={k: v for k, v in li.items() if k not in ("t_p", "v_p")})
    save("round7_fast_track.json", dict(out, **{k: v for k, v in r.items() if k != "pats"}))
    return out


# =============================================================================
# S2.1  the plate heat as a sum over the plate's section
# =============================================================================

def plate_grid(z, k):
    tag = f"{k:.2f}"
    x, y = z[f"plate_x_{tag}"], z[f"plate_y_{tag}"]
    B = np.hypot(z[f"plate_Bx_{tag}"], z[f"plate_By_{tag}"])
    xs = np.unique(np.round(x, 6)); ys = np.unique(np.round(y, 6))
    G = np.zeros((len(xs), len(ys)))
    G[np.searchsorted(xs, np.round(x, 6)), np.searchsorted(ys, np.round(y, 6))] = B
    return xs, ys, G


def part_heat():
    hr("S2.1  THE PLATE HEAT AS A SUM OVER THE PLATE'S SECTION (the 1.0 m x 200 mm fast plate, the accepted track's 2.5 mm solve)")
    meta, z = fast_solve()
    F_req = meta["F_req"]; NI_req = meta["NI_req"]; h = meta["h"]
    lv = {round(l["k"], 2): l for l in meta["levels"]}
    xs, ys, G1 = plate_grid(z, 1.0); _, _, G05 = plate_grid(z, 0.5)
    lin = float((G05.sum() / G1.sum()) / 0.5)
    print(f"  The plate's section in {G1.shape[0]} x {G1.shape[1]} cells of {h*1e3:g} mm.  At the requirement ({NI_req/1e3:.0f} kA): the field's peak in the plate {G1.max():.2f} T, its mean {G1.mean():.2f} T,")
    print(f"  its root mean square {math.sqrt((G1**2).mean()):.2f} T; half the cells are under {np.median(G1):.2f} T.  At half the ampere-turns the plate's field is {lin:.3f} times half of that: the plate is linear")
    print(f"  to {abs(lin-1)*100:.1f} percent, so each cell's swing in a pass is its field at the requirement times the pass's ampere-turns over the requirement's.")
    print(f"  Hysteresis: the book's {lc.FERRO_HYSTERESIS_COEFF} J/m3 per cycle at 1.2 T, as (B / 1.2)^1.6, summed over the cells, each with its own swing.  Eddy: the book's 1-D diffusion of the field into")
    print(f"  the plate from its face (bearing_design.eddy_cycle, resistivity {lc.FERRO_RESISTIVITY:g} ohm m, relative permeability {lc.FERRO_MU_R}), for each column of cells across the width with the swing")
    print(f"  of its own face cell; the eddy energy goes as the field squared, so the diffusion is run once per speed.  Both formulas are the book's; what changes is the field put in.")
    hyst_unit = lc.FERRO_HYSTERESIS_COEFF * ((G1 / 1.2)**1.6).sum() * h * h         # J per m of plate per cycle at the requirement
    face = G1[:, 0]                                                                 # the cells at the gap's face
    face2 = float((face**2).sum() * h); mean2 = float((G1.mean(axis=1)**2).sum() * h)
    rows = []; per_u = {}; per_u_pk = {}; per_u_mean = {}
    print(f"  {'u km/s':>7} {'winding':>7} {'pull/track kN/m':>15} {'kA':>5} {'of the requirement':>18} {'hold ms':>7} | {'summed: hyst J/m':>16} {'eddy J/m':>9} | {'one field, the peak: hyst':>25} {'eddy':>7} | {'per m of ring, the pair, summed J/m':>36}")
    for u_km in (1, 3, 5, 10, 20, 30):
        u = u_km * 1e3; F, kind = r4.pull_track_of_u(u)
        s = math.sqrt(F / F_req)
        L_on = r5.hold_time(u)[1]; v_sw = u + bd.V_REL
        e1 = bd.eddy_cycle(1.0, bd.GAP0 / v_sw, L_on / v_sw, d_p=FAST["T"], mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY)       # J/m2 at 1 T
        hyst = hyst_unit * s**1.6; eddy = e1 * face2 * s * s; eddy_mean = e1 * mean2 * s * s
        Bp = G1.max() * s
        hyst_pk = lc.FERRO_HYSTERESIS_COEFF * (Bp / 1.2)**1.6 * FAST["T"] * FAST["Wp"]; eddy_pk = e1 * Bp * Bp * FAST["Wp"]
        rows.append(dict(u=u, kind=kind, F=F, s=s, hold=L_on / v_sw, hyst=hyst, eddy=eddy, eddy_mean=eddy_mean, hyst_pk=hyst_pk, eddy_pk=eddy_pk))
        per_u[u] = 2 * (hyst + eddy); per_u_pk[u] = 2 * (hyst_pk + eddy_pk); per_u_mean[u] = 2 * (hyst + eddy_mean)
        print(f"  {u_km:7.0f} {kind:>7} {F/1e3:15.0f} {NI_req*s/1e3:5.0f} {s:18.2f} {L_on/v_sw*1e3:7.1f} | {hyst:16.1f} {eddy:9.1f} | {hyst_pk:25.1f} {eddy_pk:7.1f} | {2*(hyst+eddy):36.1f}")
    li = r6.launch_profile(lv[1.0]["W"], F_req); t_p, v_p = li["t_p"], li["v_p"]
    def launch(per):
        us = sorted(per); E = 0.0
        for i in range(1, len(t_p)):
            vv = float(v_p[i]); dtt = float(t_p[i] - t_p[i - 1])
            if vv < us[0]:
                continue
            E += float(np.interp(vv, us, [per[u] for u in us])) * (vv + bd.V_REL) * dtt
        return E
    E = launch(per_u); E_pk = launch(per_u_pk); E_mean = launch(per_u_mean)
    E_h = launch({u: 2 * r["hyst"] for u, r in zip(sorted(per_u), rows)})
    print(f"  Per launch (the heaviest), the pair of plates along the run: {E/1e9:.0f} GJ summed over the section ({E_h/1e9:.0f} hysteresis, {(E-E_h)/1e9:.0f} eddy), {E/E_LAUNCH*100:.4f} percent of the launch's 4.6 PJ;")
    print(f"  {E*7/86400/ring.L_RING:.2f} W per meter of ring at 7 a day and {E*14/86400/ring.L_RING:.2f} at 14; {E*14/86400/1e6:.0f} MW around the ring at 14.  With one field for the whole plate, the solve's peak")
    print(f"  ({G1.max():.2f} T), the same formulas give {E_pk/1e9:.0f} GJ (round 6b's 819 GJ was this with the 200 x 50 mm bundle); with each column's mean field instead of its face cell, {E_mean/1e9:.0f} GJ.")
    # the cable's heat balance: lim_physics has the cable's radiative equilibrium; the plate heat set beside the LIM heat it is sized for
    import lim_config as cfg, lim_physics as phys
    q14 = E * 14 / 86400 / ring.L_RING; q7 = E * 7 / 86400 / ring.L_RING
    out = dict(E=E, E_h=E_h, E_pk=E_pk, E_mean=E_mean, rows=rows, peak=float(G1.max()), mean=float(G1.mean()), rms=float(math.sqrt((G1**2).mean())), lin=lin, q7=q7, q14=q14)
    cfg.configure()
    def T_cable(q_per_m):                                                # the LIM model's two-node balance with q W per meter of ring put into the cable
        p_equiv = q_per_m * cfg.LIM_SPACING / cfg.LIMS_PER_SITE          # the same heat per meter as a per-LIM eddy power
        return phys.calc_two_node_thermal(p_equiv, cfg.LIMS_PER_SITE, cfg.LIM_SITES, cfg.CABLE_EMISSIVITY, cfg.CABLE_RADIATING_AREA_PER_M, cfg.CASING_WALL_EMISSIVITY,
                                          cfg.WARM_RADIATOR_AREA_PER_M, cfg.WARM_RADIATOR_EMISSIVITY, cfg.Q_ABSORBED_PER_M)
    T0c, T0w = T_cable(0.0); T7c, _ = T_cable(q7); T14c, T14w = T_cable(q14)
    out.update(T0=T0c, T7=T7c, T14=T14c, Tw0=T0w, Tw14=T14w)
    print(f"  How it leaves, and the rise: lim_physics.calc_two_node_thermal is the cable's heat balance of the deployment runs (the cable radiates to the casing's wall from")
    print(f"  {cfg.CABLE_RADIATING_AREA_PER_M:g} m2 per meter at an emissivity of {cfg.CABLE_EMISSIVITY:g}, the wall to space through the warm radiator).  With no LIM heat (operation) it gives the cable {T0c:.1f} K; with the")
    print(f"  plate heat spread around the ring, {q7:.2f} W/m at 7 a day and {q14:.2f} at 14, it gives {T7c:.1f} and {T14c:.1f} K: a rise of {T7c-T0c:.2f} and {T14c-T0c:.2f} K.  That is the mean around the ring;")
    print(f"  the heat is laid down along the sled's run, lap after lap, and a model of the cable's temperature along the run through a launch is in no script.")
    save("round7_heat.json", out)
    return out


# =============================================================================
# S1  the coupling of a load track into a support track of two separate loops, and round 4's cases with it
# =============================================================================

def part_coupling(gc=0.05, d_edge=0.2, h=0.01):
    import magnetostatic_2d as ms
    import compensation_fields as cf
    hr("S1  THE LOAD TRACK'S COUPLING INTO A SUPPORT TRACK OF TWO SEPARATE LOOPS")
    half = half_solves()
    hs = half[f"{gc:.3f}"]
    hs0 = half["0.000"]
    r_sup, r_lift = be.load_records()
    geom = dict(bo.BEST_SUPPORT); geom["W"] = SUP["Wl"]; geom["gc"] = gc
    ts = be.Track(geom, h, w_tape=0.012, I_ref=bo.I_REF_12, NI_ref=300e3)
    xl = SUP["Wp"] / 2 + d_edge + FAST["Wp"] / 2
    y_top = -bo.WALL_DESIGN
    load_legs = [(xl - FAST["Wl"] / 2, +1), (xl + FAST["Wl"] / 2, -1)]

    def solve(NI_sup, NI_load):
        g = ms.Grid(9.0, 3.0, h, x0=-3.0, y0=-1.0)
        g.add_iron(-SUP["Wp"] / 2, SUP["Wp"] / 2, bd.GAP0, bd.GAP0 + SUP["T"], cf.PLATE_MU)
        g.add_iron(xl - FAST["Wp"] / 2, xl + FAST["Wp"] / 2, bd.GAP0, bd.GAP0 + FAST["T"], cf.PLATE_MU)
        for (x, sg), a in zip(ts.legs, NI_sup):
            if a:
                g.add_current(x - ts.wb / 2, x + ts.wb / 2, y_top - ts.hb, y_top, a)
        for x, sg in load_legs:
            if NI_load:
                g.add_current(x - FAST["bw"] / 2, x + FAST["bw"] / 2, y_top - FAST["bh"], y_top, sg * NI_load)
        g.solve()
        A = [float(g.A[g.box(x - ts.wb / 2, x + ts.wb / 2, y_top - ts.hb, y_top)].mean()) for x, sg in ts.legs]
        return [sum(s * A[i] for i, s in loop) for loop in ts.loops]

    L = np.zeros((2, 2))
    for k in range(2):
        L[:, k] = np.array(solve(ts.NI_loop(k, 1e3), 0.0)) / 1e3
    M = np.array(solve([0.0] * len(ts.legs), 1e3)) / 1e3
    NI_l = r_lift["NI_req"] * math.sqrt(F_FULL / r_lift["F_req"])
    fast = load("round7_fast_track.json")
    if fast:
        NI_l = fast["NI_req"] * math.sqrt(F_FULL / fast["F_req"])
    NI0 = hs["normal"]["NI"]
    dNI = -np.linalg.solve(L, M) * NI_l                                  # each loop keeps its own flux
    dF_split = float(dNI.sum() / NI0)
    Ls = L.sum(); Ms = M.sum(); dNI_series = -Ms / Ls * NI_l             # the two loops in series: the sum of their fluxes kept (the present figure-eight)
    dF_series = float(2 * dNI_series / NI0)
    print(f"  Support plate {SUP['Wp']:.1f} m x {SUP['T']*1e3:.0f} mm with two 0.4 m loops {gc*1e3:.0f} mm apart at the center, the load track (1.0 m x 200 mm plate, 0.6 m loop, 204 x 50 mm bundles) {d_edge:.1f} m")
    print(f"  from it, plate edge to plate edge; linear plates, {h*1e3:.0f} mm cells.  Loop 1 is the far loop, loop 2 the one next to the load track.")
    print(f"  Flux per ampere-turn, uWb/m: loop 1 own {L[0,0]*1e6:.3f}, loop 2 own {L[1,1]*1e6:.3f}, between them {L[0,1]*1e6:.3f}; from the load track into loop 1 {M[0]*1e6:.4f}, into loop 2 {M[1]*1e6:.4f}.")
    print(f"  The load track at its full {F_FULL/1e3:.0f} kN/m, {NI_l/1e3:.0f} kA; the support at its solved normal {NI0/1e3:.1f} kA per leg (rounds 5 to 6b took {r6.load_final()[0]['NI_req']*math.sqrt(W_HALF/r6.load_final()[0]['F_req'])/1e3:.0f} kA, scaled")
    print(f"  as the square root of the force from the partly saturated cap).  Each loop keeping its own flux: loop 1 {dNI[0]/1e3:+.2f} kA, loop 2 {dNI[1]/1e3:+.2f} kA; the track's force {dF_split*100:+.2f} percent")
    print(f"  ({abs(dF_split)*W_HALF/1e3:.2f} kN/m per track, {2*abs(dF_split)*W_HALF/1e3:.2f} for the pair).  The two loops in series, as the present figure-eight, on the same solve: {dF_series*100:+.2f} percent; round 6b's")
    print(f"  9.21 percent was that figure-eight (touching bundles) on the scaled 185 kA, and is {9.21*185.0/ (hs0['normal']['NI']/1e3):.1f} percent on its solved {hs0['normal']['NI']/1e3:.1f} kA.")
    out = dict(gc=gc, L=L.tolist(), M=M.tolist(), NI_l=NI_l, NI0=NI0, dNI=dNI.tolist(), c_split=abs(dF_split), c_series=abs(dF_series))
    save("round7_coupling.json", out)
    return out


def part_cases():
    """Round 4's 10 ms cases with the split support's coupling carried (as round 6b's rerun)."""
    import bearing_round6a as a6
    hr("S1  ROUND 4'S CASES WITH THE SPLIT SUPPORT'S COUPLING")
    cp = load("round7_coupling.json") or part_coupling()
    a6.LOGS = LOGS
    _rc = bl.run_cases
    def rc(tau, prot, A2, theva, sup, lift, tabs, fname, **kw):
        return _rc(tau, prot, A2, theva, sup, lift, tabs, fname.replace("round6a", "round7"), **kw)
    bl.run_cases = rc
    a6.part_cases((0.2,), spacing_rows=[dict(d=0.2, dF=cp["c_split"])])
