#!/usr/bin/env python3
"""
The bearing, round 8: the rows of Section 3b (run through bearing_round8.py --part=rows).
Book: Orbital Ring Engineering, Volume III: Chapter 7 (Tables 7.3, 7.4, 7.6; Sections 7.5 to 7.9).
Brief: Vol_III/Claude_Code_instructions_round8_book_cable_2026-11-15.md.  Reads runs/round8_logs/round8_sup_*.json (bearing_round8.py
--part=solve-sup), runs/round7_logs/ (bearing_round7.py) and runs/bearing_round7_casing_10ms_d20cm.txt.  Writes runs/round8_logs/round8_rows.json.
"""
import json
import math
import os
import re

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round5 as r5
import bearing_round6b as r6
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_round8 as r8
import compensation_fields as cf
import final_bearing as fb
import levitation_config as lc
import magnetostatic_2d as ms

W_HALF = r7.W_HALF; W_CASING = r7.W_CASING; F_FULL = r7.F_FULL; SUP = r7.SUP; FAST = r7.FAST
X_S, X_C, SIDE = p7.X_S, p7.X_C, p7.SIDE
G_700 = W_CASING / ring.M_LOAD_M                                   # N per kg of casing at 700 km (95.0 kN/m over 12,000 kg/m)


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def sup_rows():
    rows = []
    for s in ("norm", "both", "high", "one"):
        d = r8.load(f"round8_sup_{s}.json")
        if d:
            rows += d["rows"]
    both = sorted([r for r in rows if r["kind"] == "both"], key=lambda r: r["NI"])
    one = sorted([r for r in rows if r["kind"] == "one"], key=lambda r: r["NI"])
    lin = next((r for r in rows if r["kind"] == "lin"), None)
    return both, one, lin


def interp_log(x, xs, ys):
    return math.exp(float(np.interp(math.log(x), np.log(xs), np.log(ys))))


def stray(NI_sup, NI_load, NI_lift, h=0.05):
    """All six tracks in a 30 x 30 m linear box (final_bearing D.2's method on the round 7 section): the cable SIDE square from y = 0 (underside) to SIDE;
    support plates 1.2 m x 100 mm and load plates 1.0 m x 200 mm let into the underside, lift plates into the top; the bundles' tops 110 mm from their
    plates (the 100 mm gap and 10 mm of cryostat).  Returns the field (T) 1 m below a support track's bundles and at the LIM stator points."""
    g = ms.Grid(30.0, 30.0, h, x0=-15.0, y0=-10.0)
    d = (SUP["Wl"] + 0.15 + 0.05) / 2
    y_b = -(bd.GAP0 + 0.010)
    for sx in (-1, 1):
        xs = sx * X_S; xc = sx * X_C
        g.add_iron(xs - SUP["Wp"] / 2, xs + SUP["Wp"] / 2, 0.0, SUP["T"], cf.PLATE_MU)
        g.add_iron(xc - FAST["Wp"] / 2, xc + FAST["Wp"] / 2, 0.0, FAST["T"], cf.PLATE_MU)
        g.add_iron(xc - FAST["Wp"] / 2, xc + FAST["Wp"] / 2, SIDE - FAST["T"], SIDE, cf.PLATE_MU)
        if NI_sup:
            for dx, a in ((-d - 0.2, NI_sup), (-d + 0.2, -NI_sup), (d - 0.2, -NI_sup), (d + 0.2, NI_sup)):
                g.add_current(xs + dx - 0.075, xs + dx + 0.075, y_b - 0.067, y_b, a)
        for NI, y0, y1 in ((NI_load, y_b - FAST["bh"], y_b), (NI_lift, SIDE - y_b, SIDE - y_b + FAST["bh"])):
            if NI:
                for dx, a in ((-0.3, NI), (0.3, -NI)):
                    g.add_current(xc + dx - FAST["bw"] / 2, xc + dx + FAST["bw"] / 2, y0, y1, a)
    g.solve()
    out = {"1 m below a support track": math.hypot(*g.field_at(X_S, y_b - 0.067 - 1.0)), "1 m below a load track": math.hypot(*g.field_at(X_C, y_b - FAST["bh"] - 1.0)),
           "LIM stator, center": math.hypot(*g.field_at(fb.STATOR_X, fb.STATOR_Y))}
    m_st = g.box(fb.STATOR_X - 0.1, fb.STATOR_X + 0.1, fb.STATOR_Y - 1.0, fb.STATOR_Y + 1.0)
    out["largest over the stator"] = float(np.sqrt(g.Bx[m_st]**2 + g.By[m_st]**2).max())
    return out


def part_rows():
    hr("ROUND 8, SECTION 3b: THE ROWS FOR THE CHAPTER 7 REWRITE")
    half = r7.load("round7_half_2p5.json"); sup7 = r7.load("round7_support.json"); fast = r7.load("round7_fast_track.json"); fsol = r7.load("round7_fast_solve.json")
    both, one, lin = sup_rows()
    theva, theva77 = be.load_theva()
    out = {}
    T = half["cap"]["tapes_leg"]; n_t = T / 100.0; NI_cap = half["cap"]["NI"]; NI_n = half["normal"]["NI"]
    cap = min(both, key=lambda r: abs(r["k"] - 1.0)); nrm = min(both, key=lambda r: abs(r["NI"] - NI_n))
    print(f"  The split support track (two 0.4 m loops, 50 mm between the center bundles, the 1.2 m x 100 mm plate), 2.5 mm cells.  Built at 50 percent: {T:,.0f} tapes per leg,")
    print(f"  {n_t:.2f} turns of ten layers per loop and coil, ten coils per loop, {half['cap']['tape']:.1f} kg/m of tape per track.  Solved points, both loops:")
    k_lin = lin["F"] / lin["NI"]**2 if lin else None
    print(f"  {'kA per leg':>10} {'x cap':>6} {'pull kN/m':>9} {'sat %':>6} {'plate peak T':>12} {'of the square law':>17} {'worst tape A':>12} {'T':>5} {'deg':>4} {'A per tape':>10} {'share %':>7} {'A per coil':>10} {'flux per loop Wb/m':>18}")
    pts = []
    for r in both:
        dev = r["F"] / (k_lin * r["NI"]**2) if k_lin else float("nan")
        It = r["NI"] / T
        pts.append(dict(NI=r["NI"], k=r["k"], F=r["F"], sat=r["sat"], B_plate=r["B_plate"], dev=dev, ic=r["ic"], B=r["B"], th=r["th"], I_tape=It, share=It / r["ic"], I_coil=10 * It, phi=float(np.mean(np.abs(r["phi"]))), ic77=r.get("ic77")))
        p = pts[-1]
        print(f"  {p['NI']/1e3:10.1f} {p['k']:6.3f} {p['F']/1e3:9.1f} {p['sat']*100:6.1f} {p['B_plate']:12.2f} {dev:17.3f} {p['ic']:12.1f} {p['B']:5.2f} {p['th']:4.0f} {It:10.1f} {It/p['ic']*100:7.1f} {10*It:10.0f} {p['phi']:18.4f}")
    NIs = [p["NI"] for p in pts]
    out["support_points"] = pts; out["k_lin"] = k_lin; out["tapes_leg"] = T; out["n_t"] = n_t

    # ---- 1. Table 7.3
    hr("1. TABLE 7.3: HOW HARD TO RUN THE WORST TAPE OF A SUPPORT TRACK")
    ic70 = cap["ic"]; ic77 = cap["ic77"]
    print(f"  The worst 12 mm tape at the cap ({NI_cap/1e3:.1f} kA per leg): {ic70:.1f} A at 70 K and {ic77:.1f} A at 77.5 K, at {cap['B']:.2f} T and {cap['th']:.0f} degrees from its normal; the limit falls")
    print(f"  {(ic70-ic77)/7.5:.1f} A per kelvin, {(ic70-ic77)/7.5/ic70*100:.1f} percent of the 70 K value per kelvin (linear between the two measurements, as the present table).  Normal operation is {NI_n/1e3:.1f} kA,")
    print(f"  {NI_n/NI_cap*100:.1f} percent of the cap.")
    print(f"  {'share at the cap':>16} {'tape A at the cap':>17} {'warming K':>9} {'tape A, normal':>14} {'tapes per leg':>13} {'tape kg/m per track':>19}")
    t73 = []
    for s in (0.5, 0.6, 0.7, 0.8):
        I = s * ic70; tapes = NI_cap / I; dT = 7.5 * (1 - s) * ic70 / (ic70 - ic77); In = NI_n / tapes
        tape = 4 * tapes * 0.012 * bd.TAPE_T * bd.TAPE_RHO
        t73.append(dict(share=s, I_cap=I, dT=dT, I_normal=In, tapes=tapes, tape=tape))
        print(f"  {s*100:15.0f}% {I:17.1f} {dT:9.1f} {In:14.1f} {tapes:13,.0f} {tape:19.1f}")
    out["t73"] = dict(ic70=ic70, ic77=ic77, B=cap["B"], th=cap["th"], rows=t73, frac_normal=NI_n / NI_cap, share_normal=half["normal"]["share"])

    # ---- 2. Table 7.4
    hr("2. TABLE 7.4: THE PULL AGAINST THE AMPERE-TURNS")
    shares = [p["share"] for p in pts]
    NI70 = float(np.interp(0.7, shares, NIs)); F70 = interp_log(NI70, NIs, [p["F"] for p in pts]); sat70 = float(np.interp(NI70, NIs, [p["sat"] for p in pts])); ic70r = float(np.interp(NI70, NIs, [p["ic"] for p in pts]))
    devs = [p["dev"] for p in pts]
    def NI_dev(x):
        for a, b in zip([dict(NI=0.0, dev=1.0)] + pts, pts):            # the square law holds as the ampere-turns go to zero
            if a["dev"] >= x >= b["dev"]:
                return a["NI"] + (a["dev"] - x) / (a["dev"] - b["dev"]) * (b["NI"] - a["NI"])
        return None
    NI95 = NI_dev(0.95); NI90 = NI_dev(0.90)
    al = half["alone"]
    t74 = [dict(state="normal operation, both loops", NI=NI_n, F=nrm["F"], I_tape=NI_n / T, share=half["normal"]["share"], sat=half["normal"]["sat"]),
           dict(state="one loop alone carrying its track's 47.5 kN/m", NI=al["NI"], F=al["F"], I_tape=al["NI"] / T, share=al["share"], sat=al["sat"]),
           dict(state="the cap, the requirement", NI=NI_cap, F=cap["F"], I_tape=NI_cap / T, share=NI_cap / T / cap["ic"], sat=cap["sat"]),
           dict(state="the rating, the worst tape at 70 percent", NI=NI70, F=F70, I_tape=NI70 / T, share=0.7, sat=sat70)]
    for lab, NIx in (("5 percent under the square law", NI95), ("10 percent under the square law", NI90)):
        if NIx:
            below = NIx < NIs[0]
            t74.append(dict(state=lab + (" (between zero and the normal point; no solve there)" if below else ""), NI=NIx, F=(k_lin * NIx**2 * float(lab.split()[0]) / 100 * 0 + k_lin * NIx**2 * (1 - float(lab.split()[0]) / 100)) if below else interp_log(NIx, NIs, [p["F"] for p in pts]),
                            I_tape=NIx / T, share=float(np.interp(NIx, [0.0] + NIs, [0.0] + shares)), sat=float(np.interp(NIx, [0.0] + NIs, [0.0] + [p["sat"] for p in pts]))))
    print(f"  One split support track.  The square law is the linear plate's: {k_lin*1e6:.4f} N/m per kA^2 per leg (both loops).  The pull is under it by {(1-pts[0]['dev'])*100:.1f} percent in normal operation and")
    print(f"  {(1-min(pts, key=lambda q: abs(q['k']-1.0))['dev'])*100:.1f} percent at the cap.")
    print(f"  {'state':>48} {'kA per leg':>10} {'A per coil':>10} {'A per tape':>10} {'share %':>7} {'pull kN/m':>9} {'sat %':>6}")
    for r in t74:
        r["I_coil"] = 10 * r["I_tape"]
        print(f"  {r['state']:>48} {r['NI']/1e3:10.1f} {r['I_coil']:10.0f} {r['I_tape']:10.1f} {r['share']*100:7.1f} {r['F']/1e3:9.1f} {r['sat']*100:6.1f}")
    out["t74_support"] = t74; out["rating_support"] = dict(NI=NI70, F=F70, sat=sat70, ic=ic70r, margin=F70 / cap["F"])
    # the fast track, from round 7's solve and its saved tape fields
    meta, z = p7.fast_solve(); lv = {round(l["k"], 2): l for l in meta["levels"]}
    ks = (1.0, 1.25, 1.5, 2.0); pat = (p7.PATTERN_NEG, p7.PATTERN_POS)
    fNI = np.array([lv[k]["NI"] for k in ks]); fF = np.array([lv[k]["F"] for k in ks]); fsat = np.array([lv[k]["sat"] for k in ks])
    fic = np.array([r6.worst_ic(p7.turned(z, p7.level_fields(z, k), pat)[0], "h", theva)[0] for k in ks])
    Tf = fast["tapes"]; fshare = (fNI / Tf) / fic
    kf = lv[0.5]["F"] / lv[0.5]["NI"]**2
    fdev = fF / (kf * fNI**2)
    def frow(state, NI):
        return dict(state=state, NI=NI, F=interp_log(NI, fNI, fF), I_coil=NI / 200.0, I_tape=NI / Tf, share=float(np.interp(NI, fNI, fshare)), sat=float(np.interp(NI, fNI, fsat)), dev=float(np.interp(NI, fNI, fdev)))
    NI198 = interp_log(F_FULL, fF, fNI)
    f74 = [frow("lift, 196.5 kN/m", meta["NI_req"]), frow("load, 198 kN/m", NI198), frow("the rating, the worst tape at 70 percent", float(np.interp(0.7, fshare, fNI)))]
    fx = None
    for i in range(len(ks) - 1):
        if fdev[i] >= 0.95 >= fdev[i + 1]:
            fx = fNI[i] + (fdev[i] - 0.95) / (fdev[i] - fdev[i + 1]) * (fNI[i + 1] - fNI[i])
    if fx:
        f74.append(frow("5 percent under the square law", float(fx)))
    print(f"\n  One fast track (the 1.0 m x 200 mm plate, the 0.6 m loop, {Tf:,.0f} tapes of 3 mm per leg, 20 coils of 10 turns over every point).  The square law from the solve at half the")
    print(f"  requirement's ampere-turns.  Solved levels: " + ", ".join(f"{k:g} x: {F/1e3:.0f} kN/m, sat {s*100:.1f} %, {d:.3f} of the square law, tape at {sh*100:.0f} %" for k, F, s, d, sh in zip(ks, fF, fsat, fdev, fshare)))
    print(f"  {'state':>48} {'kA per leg':>10} {'A per coil':>10} {'A per tape':>10} {'share %':>7} {'pull kN/m':>9} {'sat %':>6} {'of the square law':>17}")
    for r in f74:
        print(f"  {r['state']:>48} {r['NI']/1e3:10.1f} {r['I_coil']:10.0f} {r['I_tape']:10.1f} {r['share']*100:7.1f} {r['F']/1e3:9.1f} {r['sat']*100:6.1f} {r['dev']:17.3f}")
    out["t74_fast"] = f74; out["fast_levels"] = [dict(k=k, NI=float(a), F=float(b), sat=float(c), dev=float(d), share=float(e), ic=float(f)) for k, a, b, c, d, e, f in zip(ks, fNI, fF, fsat, fdev, fshare, fic)]

    # ---- 3. Table 7.6 and 4. the capacity sentences
    hr("3. TABLE 7.6: THE OPERATING POINTS OF A SUPPORT TRACK;  4. THE CAPACITY SENTENCES")
    ld = p7.lead(); hm = sup7["half"]["0.050"]
    phi_n = float(np.mean(np.abs(half["normal"]["phi"]))); phi_c = pts[[p["NI"] for p in pts].index(cap["NI"])]["phi"]
    lam = n_t * 1000.0 * (phi_c - phi_n)                                  # volt-seconds per coil, ten coils of a loop together, from the normal state to the cap
    L_eff = lam / (10 * (NI_cap - NI_n) / T)
    # the load tracks' duty if the standing load is on from the start: the pull needed falls as the support rises
    phis = [p["phi"] for p in pts if p["NI"] <= NI_cap * 1.001]; Fs = [p["F"] for p in pts if p["NI"] <= NI_cap * 1.001]
    def duty(V, n=2000):
        Tt = lam / V; s = 0.0
        for i in range(n):
            ph = phi_n + (phi_c - phi_n) * (i + 0.5) / n
            Fsup = float(np.interp(ph, phis, Fs))
            s += max(0.0, (cap["F"] - Fsup)) / F_FULL * Tt / n
        return Tt, s
    park = {}
    for V in (210.0, 1000.0):
        Tt, s = duty(V); park[V] = dict(T=Tt, duty=s)
    print(f"  Parked 5,000 t on 100 m: {5e6*G_700/100/1e3:.0f} kN/m, 198 kN/m to each side, so each support track goes to its cap of {cap['F']/1e3:.1f} kN/m.  The loops' flux goes from {phi_n:.4f} to")
    print(f"  {phi_c:.4f} Wb per meter and turn, {lam:,.0f} volt-seconds per coil with the ten coils of a loop together (an effective {L_eff:.1f} H over the step, both loops of the track rising; the")
    print(f"  small-signal value is {hm['L_loop']+hm['M_loop']:.1f} H).  The supply takes " + " and ".join(f"{park[V]['T']:.1f} s at {V:,.0f} V" for V in park) + f", the coil current from {10*NI_n/T:.0f} to {10*NI_cap/T:,.0f} A.")
    print(f"  If the whole load were on from the start, the load track of each side would hold what the support does not yet give: the integral of its pull over its full pull is " + " and ".join(f"{park[V]['duty']:.1f} s at {V:,.0f} V" for V in park) + ",")
    print(f"  against the fast leads' budget of {ld['budget20']:.2f} s with all 20 families connected.  So the pulse-sized leads cannot hold a parked spacecraft through the takeover: the load has to")
    print(f"  come on no faster than the supply raises the loops (the support follows the loading), or the load tracks' leads must be sized for it.  No script sizes them for that.")
    out["parked"] = dict(F=cap["F"], lam=lam, L_eff=L_eff, I_n=10 * NI_n / T, I_cap=10 * NI_cap / T, times={str(int(V)): park[V] for V in park}, budget20=ld["budget20"], budget10=ld["budget10"])
    # one loop alone: its pull curve
    opts = [dict(NI=al["NI"], F=al["F"], share=al["share"], sat=al["sat"])] + [dict(NI=r["NI"], F=r["F"], share=r["NI"] / T / r["ic"], sat=r["sat"]) for r in one if r["NI"] > al["NI"] * 1.05]
    oNI = [o["NI"] for o in opts]; oF = [o["F"] for o in opts]; osh = [o["share"] for o in opts]; osat = [o["sat"] for o in opts]
    print(f"  One loop alone (the other loop of its track empty): " + "; ".join(f"{o['NI']/1e3:.0f} kA, {o['F']/1e3:.1f} kN/m, tape at {o['share']*100:.0f} %, sat {o['sat']*100:.0f} %" for o in opts))
    F1_cap = interp_log(NI_cap, oNI, oF); sh1_cap = float(np.interp(NI_cap, oNI, osh))
    NI1_70 = float(np.interp(0.7, osh, oNI)) if max(osh) >= 0.7 else None
    F1_70 = interp_log(NI1_70, oNI, oF) if NI1_70 else None
    # the loss of a loop: the jump, the trim
    rows_h = {r["loss"]: r for r in sup7["rows"] if "50 mm" in r["name"] and r["name"].startswith("2.")}
    tr_o = rows_h["one loop, others trimmed"]; tr_i = rows_h["an inner loop, trimmed"]
    def share_one(I):
        return float(np.interp(I * T / 10, oNI, osh))
    ops = []
    ops.append(dict(case="normal operation, each loop", I=10 * NI_n / T, NI=NI_n, F="47.5 kN/m per track", share=half["normal"]["share"]))
    ops.append(dict(case="parked 5,000 t, the support at its cap", I=10 * NI_cap / T, NI=NI_cap, F=f"{cap['F']/1e3:.1f} kN/m per track", share=NI_cap / T / cap["ic"]))
    ops.append(dict(case="one loop lost, after the flux jump: the loop beside it", I=tr_o["I0"][2], NI=tr_o["I0"][2] * T / 10, F=f"{half['jump']['F']/1e3:.1f} kN/m on that side", share=share_one(tr_o["I0"][2])))
    ops.append(dict(case="an outer loop lost, after the trim: the lost side's inner loop", I=tr_o["I1"][2], NI=tr_o["I1"][2] * T / 10, F=f"{W_CASING*X_S/(hm['xi']+X_S)/1e3:.1f} kN/m", share=share_one(tr_o["I1"][2])))
    ops.append(dict(case="the same: the sound side's two loops", I=tr_o["I1"][0], NI=tr_o["I1"][0] * T / 10, F=f"{(W_CASING - W_CASING*X_S/(hm['xi']+X_S))/1e3:.1f} kN/m", share=float(np.interp(tr_o["I1"][0] * T / 10, NIs, shares))))
    ops.append(dict(case="an inner loop lost, after the trim: the lost side's outer loop", I=tr_i["I1"][2], NI=tr_i["I1"][2] * T / 10, F=f"{W_CASING*X_S/(hm['xo']+X_S)/1e3:.1f} kN/m", share=share_one(tr_i["I1"][2])))
    ops.append(dict(case="the same: the sound side's two loops", I=tr_i["I1"][0], NI=tr_i["I1"][0] * T / 10, F=f"{(W_CASING - W_CASING*X_S/(hm['xo']+X_S))/1e3:.1f} kN/m", share=float(np.interp(tr_i["I1"][0] * T / 10, NIs, shares))))
    print(f"  {'case':>66} {'A per coil':>10} {'kA per leg':>10} {'A per tape':>10} {'share %':>7}  pull")
    for o in ops:
        print(f"  {o['case']:>66} {o['I']:10.0f} {o['NI']/1e3:10.1f} {o['I']/10:10.1f} {o['share']*100:7.1f}  {o['F']}")
    # one loop out with a standing load w per meter: F_r (x_r + X_S) = (W + w) X_S
    stand = []
    for lab, x_r in (("an outer loop out (the inner loop of that side carries)", hm["xi"]), ("an inner loop out (the outer loop carries)", hm["xo"])):
        for lim, F1, sh in (("the loop at the cap's coil current", F1_cap, sh1_cap), ("the loop at its rating (70 percent)", F1_70, 0.7)):
            if F1 is None:
                continue
            w = F1 * (x_r + X_S) / X_S - W_CASING; F_l = W_CASING + w - F1
            stand.append(dict(case=lab, limit=lim, F1=F1, share=sh, w=w, tonnes=w * 100 / G_700 / 1e3, F_other=F_l))
            print(f"  One loop out with a standing load, {lab}, {lim}: the loop gives {F1/1e3:.1f} kN/m (tape at {sh*100:.0f} %), the standing load may be {w/1e3:.0f} kN/m,")
            print(f"    {w*100/G_700/1e3:,.0f} tonnes on 100 m; the sound track then carries {F_l/1e3:.1f} kN/m of its {cap['F']/1e3:.1f} cap.")
    out["ops"] = ops; out["standing"] = stand; out["one_loop"] = dict(points=opts, F_cap=F1_cap, share_cap=sh1_cap, NI70=NI1_70, F70=F1_70)
    out["change"] = dict(jump=hm["jump"], F_jump=half["jump"]["F"], outer=dict(I0=tr_o["I0"], I1=tr_o["I1"], V_min=tr_o["V_min"], int_hi=tr_o["int_hi"], T_hi=tr_o["T_hi"]), inner=dict(I0=tr_i["I0"], I1=tr_i["I1"], V_min=tr_i["V_min"], int_hi=tr_i["int_hi"], T_hi=tr_i["T_hi"]))
    capac = dict(track_over_casing=cap["F"] / W_CASING, pair=2 * cap["F"], pair_over=2 * cap["F"] / W_CASING, rating=F70, rating_over_req=F70 / cap["F"], loop_cap=F1_cap, loop_over_casing=F1_cap / W_CASING,
                 loop70=F1_70, loop70_over_casing=(F1_70 / W_CASING if F1_70 else None))
    print(f"\n  Capacity: one track carries {cap['F']/1e3:.1f} kN/m, {capac['track_over_casing']:.2f} times the casing's {W_CASING/1e3:.1f}; the pair {capac['pair']/1e3:.0f} kN/m, {capac['pair_over']:.2f} times; the rating {F70/1e3:.0f} kN/m,")
    print(f"  {(F70/cap['F']-1)*100:.0f} percent over the requirement.  One loop alone gives {F1_cap/1e3:.1f} kN/m at the cap's coil current, {capac['loop_over_casing']:.2f} times the casing"
          + (f", and {F1_70/1e3:.1f} kN/m at its own rating, {F1_70/W_CASING:.2f} times." if F1_70 else "; its own 70 percent point is above the solved range."))
    out["capacity"] = capac

    # ---- 5. failures, from round 7's rerun of round 4's cases
    hr("5. THE FAILURES THE CHAPTER QUOTES")
    p = os.path.join(r7.OUT, "bearing_round7_casing_10ms_d20cm.txt"); txt = open(p, encoding="utf-8", errors="replace").read()
    m3 = re.search(r"3\. one lift track lost at 30 km/s.*?\)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", txt)
    mf = re.search(r"fastest launch inside 20 mm with one lift track out: (\d+) km/s; at 30 km/s the casing moves ([\d.]+) mm \(gap ([\d.]+) mm\)", txt)
    mn = re.search(r"B\.6\s+(\d+) cases in all", txt); mg = re.search(r"Smallest gap in any single failure: ([\d.]+) mm in '([^']+)'\.", txt)
    cases7 = open(os.path.join(r7.LOGS, "round7_cases.txt"), encoding="utf-8", errors="replace").read()
    m1 = re.search(r"d = 0\.20 m, c = [\d.]+ %\s+([\d.]+) / ([\d.]+), gap ([\d.]+)", cases7)
    pull_1 = r5.r4.pull_track_of_u(1000.0)[0] if hasattr(r5, "r4") else None
    import bearing_round4 as r4
    F1k = r4.pull_track_of_u(1000.0)[0]
    fail = dict(lift_lost=dict(exc=float(m3.group(1)), gap=float(m3.group(2))), fastest=int(mf.group(1)), n_cases=int(mn.group(1)), gap_single=float(mg.group(1)), which=mg.group(2),
                one_load_1kms=dict(exc=float(m1.group(2)), gap=float(m1.group(3))), pull_1kms_track=F1k, mass_by_pull=5000.0 * F_FULL / (2 * F1k), parked_by_pull=5000.0 * F_FULL / (2 * F_FULL))
    print(f"  From {p} (round 7's rerun of the cases with the split support's coupling, the 10 ms rise, the 10 mm protection, the support frozen):")
    print(f"  one lift track lost at 30 km/s: the casing moves {fail['lift_lost']['exc']:.1f} mm, the smallest gap {fail['lift_lost']['gap']:.1f} mm; the fastest launch of 5,000 t inside 20 mm with one lift track out: {fail['fastest']} km/s;")
    print(f"  {fail['n_cases']} cases, no single failure in contact, the smallest single-failure gap {fail['gap_single']:.1f} mm ({fail['which']}).  One load track out: the other alone at 1 km/s with the full 5,000 t moves")
    print(f"  the casing {fail['one_load_1kms']['exc']:.1f} mm, gap {fail['one_load_1kms']['gap']:.1f} mm.  By the pull alone (the pass at 1 km/s asks {2*F1k/1e3:.0f} kN/m of the pair, a track gives {F_FULL/1e3:.0f}): {fail['mass_by_pull']:,.0f} t launched;")
    print(f"  parking no longer rests on the load tracks (the support carries a standing load at its cap), so a load track out does not limit the parked mass.")
    out["failures"] = fail

    # ---- 6. the stray field
    hr("6. THE STRAY FIELD")
    NI_lift = meta["NI_req"]; NI_lift70 = f74[2]["NI"]
    sf = {"normal operation (support at 47.5 kN/m, fast tracks off)": stray(NI_n, 0.0, 0.0), "the support tracks at their cap (245.5 kN/m), fast off": stray(NI_cap, 0.0, 0.0),
          "the support tracks at their rating, fast off": stray(NI70, 0.0, 0.0), "a launch at 30 km/s: support normal, the lift tracks at 196.5 kN/m": stray(NI_n, 0.0, NI_lift),
          "the lift tracks at their rating, support normal": stray(NI_n, 0.0, NI_lift70), "the load tracks at 198 kN/m, support normal": stray(NI_n, NI198, 0.0)}
    print(f"  All six tracks in one 30 x 30 m box, 50 mm cells, the plates linear (final_bearing D.2's method on the round 7 section).  LIM stators at x = +-{fb.STATOR_X:.2f} m, 2 m tall, centered")
    print(f"  {fb.STATOR_Y - SIDE/2:.1f} m above the cable's centerline [UNVERIFIED position, as final_bearing].  The plates saturate at the cap and the rating (35 percent and more), which a linear plate does not carry.")
    for k, v in sf.items():
        print(f"  {k:>72}: " + ", ".join(f"{a} {b*1e3:.2f} mT" for a, b in v.items()))
    out["stray"] = sf

    # ---- 7. the coil ends
    hr("7. THE COIL ENDS")
    Bp = nrm["B_plate"]; ell = 2 * (bd.GAP0 + 0.05)
    face = r8.load("round8_face.json"); hf = face["h"]
    Bf = np.array(face["both"]["B_face"]); B1 = np.array(face["one"]["B_face"])
    E1 = 0.5 * bd.eddy_cycle(1.0, ell / bd.V_REL, 0.0, d_p=SUP["T"], mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY)      # J/m2 for a step of 1 T; the energy goes as the step squared
    f2_one = float((B1**2).sum() * hf)                                    # T2 m: one loop's own face field squared, across the plate
    f2_pk = Bp**2 * SUP["Wp"] / 2
    print(f"  The plate's field in normal operation, 5 mm cells: at the gap's face it peaks at {Bf.max():.2f} T with a root mean square of {math.sqrt((Bf**2).mean()):.2f} T across the {SUP['Wp']:.1f} m; the peak")
    print(f"  anywhere in the plate is {face['both']['B_max']:.2f} T ({Bp:.2f} T at 2.5 mm), in a few cells above the bundles.  One loop alone at the same current: face peak {B1.max():.2f} T, rms {math.sqrt((B1**2).mean()):.2f} T.")
    print(f"  A coil end takes a tenth of its own loop's field away over {ell*1e3:.0f} mm; the step's eddy heat is half a rise-and-fall cycle of the 1-D diffusion, {E1:.1f} J/m2 for 1 T, as the step squared,")
    print(f"  summed across the plate with the one-loop face field: {f2_one:.3f} T2 m.  (With one field for the half plate, the peak cell's {Bp:.2f} T: {f2_pk:.3f} T2 m, {f2_pk/f2_one:.1f} times more.)")
    ends = []
    for Lc in (1000.0, 2000.0, 10000.0):
        for lab, nn in (("staggered, a tenth of the loop's field per end", 10), ("ending together, the whole field", 1)):
            steps = 2 * 2 * nn / Lc                                       # per meter of one track: 2 ends x 2 loops x n coils / L
            P = 2 * steps * E1 * f2_one / nn**2 * bd.V_REL                # two tracks
            P_pk = 2 * steps * E1 * f2_pk / nn**2 * bd.V_REL
            ends.append(dict(L=Lc, n=10, lab=lab, dB=float(B1.max()) / nn, P=P, MW=P * ring.L_RING / 1e6, P_pk=P_pk))
            print(f"  {Lc/1e3:.0f} km coils, ten per bundle, {lab}: {steps:.3f} steps per meter of track, {P:.2f} W per meter of ring, {P*ring.L_RING/1e6:.0f} MW around the ring (with the peak cell's field: {P_pk:.1f} W/m)")
    print(f"  (final_bearing C.3's model, a coil end stepping the field over twice the distance from the leg to the plate [UNVERIFIED]; the 100 mm plate; the heat is also the drag's power, which the")
    print(f"  LIMs make up.  The two loops of a track are separate coils now; the count of ends per meter is the same.)")
    out["face"] = dict(face_peak=float(Bf.max()), face_rms=float(math.sqrt((Bf**2).mean())), one_peak=float(B1.max()), one_rms=float(math.sqrt((B1**2).mean())), f2_one=f2_one, E1=E1)
    out["ends"] = ends; out["B_plate_normal"] = Bp

    # ---- 8. the bearing's mass on the casing
    hr("8. THE BEARING'S MASS ON THE CASING")
    tape = 2 * half["cap"]["tape"] + 4 * fast["tape"]
    E_m = 4 * fast["W"]; m_store = E_m / 1.5e3; n_s = math.ceil(fast["V"] / 2800.0); m_sw_coil = 4 * n_s * 4.3; m_sw = m_sw_coil * 3.2
    A = ld["A"]; m_lead_fast = 2 * 8960.0 * A * 1.0 * 3.2
    mass = dict(tape=tape, tape_sup=2 * half["cap"]["tape"], tape_fast=4 * fast["tape"], formers=tape, store=m_store, store_shared=m_store / 2, switches=m_sw, leads_fast=m_lead_fast)
    mass["total"] = tape + tape + m_store + m_sw + m_lead_fast; mass["total_shared"] = mass["total"] - m_store / 2
    print(f"  tape {tape:.0f} kg/m (two support tracks {2*half['cap']['tape']:.0f}, four fast tracks {4*fast['tape']:.0f}); formers, cryostats and supports at one tape mass (Part 3b's rule) {tape:.0f};")
    print(f"  the stores, film capacitors at 1.5 kJ/kg for {E_m/1e3:.0f} kJ per meter of ring: {m_store:.0f} kg/m ({m_store/2:.0f} if a site's lift and load tracks share them, since they never pull together);")
    print(f"  the switches, {n_s} IGBT modules in series per position at {fast['V']/1e3:.1f} kV, two switches and two diodes of 4.3 kg a module: {m_sw_coil:.0f} kg per coil, {m_sw:.0f} kg/m at 3.2 coils per meter;")
    print(f"  the fast coils' leads, two of {A*1e6:.0f} mm2 by 1 m of copper per coil: {m_lead_fast:.1f} kg/m.  Total {mass['total']:,.0f} kg/m ({mass['total_shared']:,.0f} with the stores shared), {mass['total']/ring.M_LOAD_M*100:.0f} percent of the 12,000 kg/m.")
    print(f"  Not weighed in any script: the support coils' strips, diodes, dump resistors and retractable leads, the supplies, the sensors, the capacitors' housings and the busbars.")
    out["mass"] = mass

    # ---- 14. the cold heat of the heaviest launch
    hr("CHAPTER 14: THE COLD HEAT OF THE HEAVIEST LAUNCH")
    li = r6.launch_profile(fast["W"], fast["F_req"]); P_full = fast["I"]**2 * be.RHO_CU_77 * 1.0 / A
    E_tape = fast["Ep"] + fast["Es"]; E_lead = li["lead_s"] * 2 * P_full
    laps = li["x_run"] / ring.L_RING
    cold = dict(E_tape=E_tape, E_lead=E_lead, x_run=li["x_run"], t_run=li["t_run"], laps=laps, per_m_tape=E_tape / ring.L_RING, per_m_lead=E_lead / ring.L_RING, P_full=P_full, lead_s=li["lead_s"])
    print(f"  The launch profile (sled_passage_pulsed.lsm_profile): {li['t_run']/60:.0f} minutes, {li['x_run']/1e3:,.0f} km of run, {laps:.2f} laps of the ring.  The tape's pulse loss {E_tape/1e9:.0f} GJ and the leads' Joule heat")
    print(f"  {E_lead/1e9:.1f} GJ ({li['lead_s']:.3e} full-current coil-seconds, two leads per coil at {P_full:.0f} W).  The run covers the ring {laps:.1f} times, so as a mean per meter of ring: {E_tape/ring.L_RING/1e3:.1f} kJ/m of tape loss")
    print(f"  and {E_lead/ring.L_RING/1e3:.2f} kJ/m of lead heat; per 500 m site {(E_tape+E_lead)/ring.L_RING*500/1e6:.2f} MJ.  The split of the heat between the laps (each at its own speed) is in the run's integral, not in a table.")
    out["cold_launch"] = cold
    r8.save("round8_rows.json", out)
    return out
