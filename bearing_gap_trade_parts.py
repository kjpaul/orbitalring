#!/usr/bin/env python3
"""
The gap trade (2026-11-16): the tables (run through bearing_gap_trade.py --part=report|pump|start).
Brief: Vol_III/Claude_Code_instructions_gap_trade_2026-11-16.md.  Reads runs/gap_trade_logs/ (bearing_gap_trade.py), runs/round7_logs/,
runs/round8_logs/ and the round 8 design run's series (runs/round8_design/ch14/).  Writes runs/gap_trade_logs/gap_trade_*.json and .txt.
"""
import csv
import json
import math
import os

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_round4 as r4
import bearing_round5 as r5
import bearing_round6b as r6
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_round8 as r8
import bearing_gap_trade as gt
import final_bearing as fb

MU0 = gt.MU0
W_HALF = r7.W_HALF; W_CASING = 2 * W_HALF; F_FULL = r7.F_FULL; SUP = r7.SUP; FAST = r7.FAST
G_700 = W_CASING / ring.M_LOAD_M
PUMP_HEAT = ((0.15, "a dynamo: 135 mW of ac loss in the stator tape at 25 Hz (Ghabeli et al. 2021, a model, read in full for Stage A)"),
             (10.0, "a self-rectifier: about 1 J per cycle in the bridge at 10 Hz (Mallett et al. 2023, read from its Fig. 8 for Stage A, approximate)"))
D_PER_KW = 2.0          # days of deployment per kW of cold load per site: round 6a's runs (8.8 days for 4.5 kW; round 5: 12.3 days for 6.1 kW)


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


# =============================================================================
# item 4: the flux pump for the trims
# =============================================================================

def support_current():
    """One support track's coil current against its pull, from round 8's 2.5 mm solves of the split support: the solved points from normal operation
    up, and below it the linear plate's square law with the solved shortfall at the normal point brought in linearly with the ampere-turns."""
    rows = r8.load("round8_rows.json"); pts = rows["support_points"]; k = rows["k_lin"]; T = rows["tapes_leg"]
    NIs = np.array([p["NI"] for p in pts]); Fs = np.array([p["F"] for p in pts]); NI_n = NIs[0]; d_n = pts[0]["dev"]
    def F_low(NI):
        return k * NI**2 * (1 - (1 - d_n) * NI / NI_n)
    def NI_of_F(F):
        if F <= 0:
            return 0.0
        if F >= Fs[0]:
            return math.exp(float(np.interp(math.log(F), np.log(Fs), np.log(NIs))))
        lo, hi = 0.0, NI_n
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if F_low(mid) < F: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)
    return (lambda F: 10 * NI_of_F(F) / T), dict(T=T, NI_n=NI_n, F_n=Fs[0], d_n=d_n, k=k, I_n=10 * NI_n / T)


def weight_series():
    p = os.path.join(r7.OUT, "round8_design", "ch14", "rising_radius_deployment_tension.csv")
    t, w = [], []
    for r in csv.DictReader(open(p, encoding="utf-8")):
        t.append(float(r["time_days"])); w.append(float(r["casing_net_weight_kN_m"]) * 1e3)
    return np.array(t), np.array(w), p


def part_pump():
    hr("ITEM 4: THE FLUX PUMP FOR THE TRIMS")
    I_of_F, c = support_current()
    hm = r7.load("round7_support.json")["half"]["0.050"]; L1 = hm["L_loop"]; M = hm["M_loop"]; L2 = L1 + M
    DP = json.load(open(os.path.join(r7.OUT, "round8_design", "design_point.json"), encoding="utf-8"))
    out = dict(L_loop=L1, M_loop=M, I_n=c["I_n"])
    day = 86400.0
    print(f"  One loop of a support track: ten staggered 1 km coils of {c['T']/100:.2f} turns; each coil sees {L1:.2f} H with the ten of its loop moving together and {L2:.2f} H when the loop")
    print(f"  beside it on the same plate moves with it (the mutual {M:.2f} H), which is the case in a trim of the track (round 7, `round7_support.json`).  Normal current {c['I_n']:.1f} A.")
    # 4.1 the proof-read arithmetic
    chk = [("1 mV moves, one loop alone", 1e-3 / L1 * day, "4.9 A a day"), ("1 mV moves, both loops of a track", 1e-3 / L2 * day, "not given"),
           ("0.05 mV moves, one loop alone", 0.05e-3 / L1 * day, "0.24 A a day"), ("0.05 mV moves, both loops", 0.05e-3 / L2 * day, "not given"),
           ("the first charge of 371 A at 1 mV, one loop alone (days)", c["I_n"] * L1 / 1e-3 / day, "76 days"), ("the same, both loops of a track together (days)", c["I_n"] * L2 / 1e-3 / day, "not given")]
    print(f"\n  4.1  The proof-read session's arithmetic:")
    for lab, v, said in chk:
        print(f"    {lab:>62}: {v:8.2f}   (the session: {said})")
    # the local slope of the current against the pull at the normal point
    dF = 50.0; s_n = (I_of_F(c["F_n"] + dF) - I_of_F(c["F_n"] - dF)) / (2 * dF)
    e_n = s_n * c["F_n"] / c["I_n"]
    d127 = 12.7 / 2 * s_n
    print(f"    At the normal point the current moves {s_n*1e3:.3f} A per kN/m of the track's pull, {e_n:.2f} percent per percent (the square law would be 0.50).  The session's 0.05 A a day put")
    print(f"    12.7 N/m a day on one track; 12.7 N/m a day is round 6a's figure for the whole casing (the fall of gravity as the ring rises 173 km), 6.35 per track, which is {d127:.3f} A a day.")
    out.update(check=[dict(label=a, value=b, said=s) for a, b, s in chk], slope_A_per_N=s_n, exponent=e_n, drift_127=d127)
    # 4.2 the drift through the deployment from the weight series
    t, w, p = weight_series()
    td = np.arange(0.0, math.floor(t[-1]) + 1e-9, 1.0); wd = np.interp(td, t, w); Id = np.array([I_of_F(x / 2) for x in wd])
    dI = np.gradient(Id, td); dw = np.gradient(wd, td)
    print(f"\n  4.2  The drift through the deployment, from the round 8 design run's series ({p}, the column casing_net_weight_kN_m): the casing's weight on the")
    print(f"  bearing is not round 6a's 12.7 N/m a day.  The casing starts in orbit with no weight on the bearing and is slowed to rest, so the whole {wd.max()/1e3:.1f} kN/m comes on during the deployment.")
    print(f"  Each track carries half.  The coil current from the solved pull curve; the voltage a coil's pump must give is L dI/dt with both loops of the track moving ({L2:.1f} H).")
    print(f"  {'day':>5} {'weight kN/m':>11} {'change N/m per day':>18} {'A per coil':>10} {'A per day':>9} {'pump uV (21.5 H)':>16} {'one loop alone uV (17.8 H)':>26}")
    rows = []
    for d in (1, 2, 5, 10, 20, 30, 50, 66, 100, 150, 200, 250, 295, 327, 350, int(td[-1])):
        i = int(d)
        rows.append(dict(day=d, w=wd[i], dw=dw[i], I=Id[i], dI=dI[i], V2=L2 * dI[i] / day, V1=L1 * dI[i] / day))
        print(f"  {d:5d} {wd[i]/1e3:11.2f} {dw[i]:18.1f} {Id[i]:10.1f} {dI[i]:9.3f} {L2*dI[i]/day*1e6:16.1f} {L1*dI[i]/day*1e6:26.1f}")
    i_mx = 1 + int(np.argmax(np.abs(dI[1:]))); i_w = int(np.argmax(np.abs(dw)))
    print(f"  The weight changes fastest on day {td[i_w]:.0f}, {dw[i_w]:.0f} N/m a day ({dw[i_w]/12.7:.0f} times round 6a's figure).  The current changes fastest at the start, where it goes as the square root")
    print(f"  of the weight: {dI[1]:.1f} A a day on day 1, {dI[10]:.2f} on day 10, {dI[30]:.2f} on day 30.  At the end (day {td[-1]:.0f}) the weight is falling {abs(dw[-1]):.0f} N/m a day and the current {abs(dI[-1]):.3f} A a day.")
    print(f"  The mean over the deployment is {Id[-1]/td[-1]:.2f} A a day ({Id[-1]:.0f} A in {td[-1]:.0f} days), {L2*Id[-1]/td[-1]/day*1e6:.0f} uV per coil.")
    out.update(series=rows, dw_max=float(dw[i_w]), dw_max_day=float(td[i_w]), dI_day1=float(dI[1]), dI_end=float(dI[-1]), dw_end=float(dw[-1]), I_end=float(Id[-1]), I_peak=float(Id.max()))
    # which pump keeps up
    print(f"\n  Which pump keeps up (one pump per coil, {L2:.1f} H; the pump runs from day 0 and the supply is not used):")
    print(f"  {'pump mV':>8} {'A per day':>9} {'keeps up with the rate from day':>31} {'has caught up with the current by day':>37} {'largest shortfall A':>19} {'source of the voltage':<70}")
    pumps = []
    for V, src in ((0.046e-3, "Geng et al. 2019, transformer-rectifier, 1150 A, time-averaged"), (0.197e-3, "Bumby et al. 2016, dynamo, 30 A, time-averaged"),
                   (0.267e-3, "Mallett et al. 2023, self-rectifier, 320 A, time-averaged"), (1.0e-3, "the session's round figure"), (1.2e-3, "Jiang et al. (in Zhou et al. 2023), dynamo at 600 Hz, 58 A, saturated")):
        rate = V / L2 * day
        ok = np.abs(dI) <= rate
        d_rate = next((int(td[i]) for i in range(len(td)) if ok[i:].all()), None)
        Ip = np.minimum(rate * td, 1e9); short = Id - Ip
        # the pump charging flat out until it reaches the demand, then following it
        d_catch = next((int(td[i]) for i in range(1, len(td)) if short[i] <= 0), None)
        pumps.append(dict(V=V, rate=rate, d_rate=d_rate, d_catch=d_catch, short=float(short.max()), src=src))
        print(f"  {V*1e3:8.3f} {rate:9.2f} {str(d_rate):>31} {str(d_catch):>37} {max(short.max(), 0.0):19.1f} {src:<70}")
    print(f"  (the voltages are the time-averaged outputs read for Stage A, `runs/bearing_layout_pumps.txt`; a shortfall is current the support does not yet carry, so weight the fast tracks")
    print(f"  or the supply through the leads must carry meanwhile.  No pump read has given 1 mV at the support coil's {c['I_n']:.0f} A: the 1.2 mV dynamo carried 58 A, the 320 A and")
    print(f"  1150 A pumps averaged 0.27 and 0.05 mV [the combination is UNVERIFIED].)")
    out["pumps"] = pumps
    # 4.3 the heat
    n_site = 4 * 10 * 500.0 / 1000.0
    print(f"\n  4.3  The pump's heat in the cold.  One pump per coil: 4 loops x 10 coils per kilometer, {n_site:.0f} pumps per 500 m site, {n_site*2*ring.L_RING/1e3/1e6/ (500/500):.2f} million around the ring.")
    P_out = L2 * abs(dI).max() / day * Id.max()
    print(f"  What a pump delivers is small: at the mean rate {L2*Id[-1]/td[-1]/day*c['I_n']*1e3:.0f} mW at {c['I_n']:.0f} A.  Its loss is a standing loss of its moving field while it runs, not a share of its output")
    print(f"  (Stage A, A.4):")
    heat = []
    for q, lab in PUMP_HEAT:
        heat.append(dict(q=q, site=q * n_site, label=lab))
        print(f"    {q:5.2f} W per pump, {q*n_site:6.1f} W per site, {q*n_site/500*ring.L_RING/1e6:6.1f} MW around the ring if every pump runs all the time: {lab}")
    print(f"  Neither figure is for a pump carrying {c['I_n']:.0f} A into 21.5 H [UNVERIFIED for this duty].  After the deployment the weight does not drift and the pumps run only for a trim.")
    out["pump_heat"] = heat; out["pumps_per_site"] = n_site
    # 4.4 the fast coils' regulation heat
    p_reg = r4.K_FLUX * 1e-3 * 0.1; p_d6 = 7 * 12.7 / 2
    j_reg = DP["joule"] * p_reg / (p_reg + p_d6)
    def weekly(dW): return DP["joule"] * (p_reg + 7 * abs(dW) / 2) / (p_reg + p_d6)
    dw_mean = float(np.abs(dw).mean())
    print(f"\n  4.4  The fast coils' regulation heat.  Round 6a's rule: with ten families connected the Joule heat in their leads goes as the mean pull a load track holds, {p_reg:.0f} N/m for a")
    print(f"  0.1 mm regulation ({r4.K_FLUX/1e6:.2f} kN/m per mm) plus the weight a weekly trim of the support leaves, 7 x 12.7 / 2 = {p_d6:.1f} N/m: {DP['joule']:,.0f} W per site, of which {j_reg:,.0f} W is the regulation.")
    print(f"  With the weight series of the run in place of 12.7 N/m a day the same rule gives, for weekly trims: {weekly(dw_mean)/1e3:.1f} kW per site at the deployment's mean drift ({dw_mean:.0f} N/m a day)")
    print(f"  and {weekly(dw[i_w])/1e3:.1f} kW at its largest ({dw[i_w]:.0f}).  So weekly trims were never enough for the deployment, and the {DP['joule']:,.0f} W was low for that scheme.")
    print(f"  With the support trimmed continuously (pumps, the loop closed) the drift's part goes and the regulation's {j_reg:,.0f} W is left.")
    print(f"  {'the support trimmed by':>44} {'regulation W':>12} {'pumps W':>8} {'conduction W':>12} {'extra cold load per site W':>26} {'against 4,411 W':>15} {'days, at 2.0 per kW':>19} {'deployment days':>15}")
    rows4 = []
    for lab, q in (("weekly trims through the leads (round 8's run)", None), ("pumps, dynamo loss", PUMP_HEAT[0][0] * n_site), ("pumps, self-rectifier loss", PUMP_HEAT[1][0] * n_site)):
        if q is None:
            extra = DP["extra"]; jj = DP["joule"]; qq = 0.0
        else:
            jj = j_reg; qq = q; extra = DP["cond"] + jj + qq
        dd = (extra - DP["extra"]) / 1e3 * D_PER_KW
        rows4.append(dict(label=lab, joule=jj, pumps=qq, cond=DP["cond"], extra=extra, d=dd, days=368.9 + dd))
        print(f"  {lab:>44} {jj:12,.0f} {qq:8.0f} {DP['cond']:12,.0f} {extra:26,.0f} {extra-DP['extra']:+15,.0f} {dd:+19.1f} {368.9+dd:15.1f}")
    print(f"  (2.0 days per kW of cold load per site: round 6a's deployment runs, 8.8 days for 4.5 kW, and round 5's 12.3 days for 6.1 kW.  The deployment was not rerun.)")
    # 4.5 the supply through the leads while the pumps cannot follow
    q_cap, aol = be.lead_heat(rows_cap := r8.load("round8_rows.json")["parked"]["I_cap"])
    q_n = be.lead_heat(c["I_n"], A_over_l=aol)[0]; q_0 = be.lead_heat(0.0, A_over_l=aol)[0]
    n_lead = 2 * n_site
    print(f"\n  4.5  While the pumps cannot follow, the supply carries the current through the leads.  A conduction-cooled copper lead sized for the cap's {rows_cap:,.0f} A (be.lead_heat, Part 1c)")
    print(f"  brings {q_0:.1f} W into the cold with no current, {q_n:.1f} W at {c['I_n']:.0f} A and {q_cap:.1f} W at the cap: {n_lead:.0f} leads per site, {q_0*n_lead:,.0f} W attached and idle, {q_n*n_lead:,.0f} W at the normal current,")
    print(f"  {q_n*n_lead/1e3*D_PER_KW:.1f} days of deployment if they stayed attached for all of it.  Round 8's run has no support lead in its cold load (the leads withdrawn).")
    out["sup_leads"] = dict(q0=q_0, qn=q_n, qcap=q_cap, n=n_lead, site_idle=q_0 * n_lead, site_normal=q_n * n_lead)
    out.update(j_reg=j_reg, weekly_mean=weekly(dw_mean), weekly_max=weekly(dw[i_w]), dw_mean=dw_mean, rows4=rows4, cond=DP["cond"], joule=DP["joule"], extra=DP["extra"])
    gt.save("gap_trade_pump.json", out)
    return out


# =============================================================================
# item 5: the start of a launch
# =============================================================================

def part_start():
    hr("ITEM 5: THE START OF A LAUNCH, FROM REST TO 1 km/s")
    rows = r8.load("round8_rows.json"); park = rows["parked"]; ld = p7.lead(); fast = r7.load("round7_fast_track.json")
    a = 0.83 * 9.80665
    t1 = 1000.0 / a; x1 = 0.5 * a * t1**2
    L_sc = bd.L_SC; L_c = 25.0; L_sup = 1000.0
    out = dict(a=a, t1=t1, x1=x1)
    print(f"  The spacecraft at 0.83 g ({a:.2f} m/s2): 1 km/s after {t1:.0f} s and {x1/1e3:.0f} km.  Parked, 5,000 t on {L_sc:.0f} m is {5e6*G_700/L_sc/1e3:.0f} kN/m, and each support track under it is at its cap of")
    print(f"  {park['F']/1e3:.1f} kN/m.  No script runs this stretch: the cases start at 1 km/s.  What follows is arithmetic on round 8's numbers, not a run.")
    # the support following through its leads
    T1k = park["times"]["1000"]["T"]; T210 = park["times"]["210"]["T"]
    u_sup = L_sc / T1k; u_sup210 = L_sc / T210
    print(f"\n  (1) The support following through its leads.  A stretch of support goes from normal to the cap in {T1k:.1f} s at 1 kV ({T210:.0f} s at 210 V): {park['lam']:,.0f} volt-seconds per coil.")
    print(f"  To follow the load, the support ahead must rise and the support behind fall by the time the load has moved its own {L_sc:.0f} m: {u_sup:.1f} m/s at 1 kV, {u_sup210:.1f} m/s at 210 V.")
    print(f"  At 0.83 g the spacecraft passes {u_sup:.1f} m/s after {u_sup/a:.1f} s and {0.5*u_sup**2/a:.1f} m.  And the support coils are {L_sup/1e3:.0f} km long in ten staggered families, so the support cannot be")
    print(f"  shaped finer than the coils: the spacecraft leaves the kilometer its raised coils cover after {math.sqrt(2*L_sup/a):.0f} s, at {math.sqrt(2*L_sup*a):.0f} m/s.")
    # the load tracks inside their leads' budget
    b20 = ld["budget20"]
    t_on = lambda u: (L_sc + L_c) / u + 2 * 10e-3
    u_load = (L_sc + L_c) / (b20 - 2 * 10e-3)
    print(f"\n  (2) The load tracks carrying the weight.  A 25 m load coil on the casing is on while the {L_sc:.0f} m spacecraft is over it, ({L_sc:.0f} + {L_c:.0f}) m / u, at its full current if the load tracks")
    print(f"  carry the whole weight (198 kN/m each).  The fast leads take {b20:.2f} s of full current for a 20 K rise with all 20 families connected (round 7, `p7.lead`).  So a pass is inside the")
    print(f"  budget above {u_load:.0f} m/s, which the spacecraft reaches after {u_load/a:.1f} s and {0.5*u_load**2/a:.0f} m.  Below that speed a coil under the spacecraft is on for longer than its leads allow.")
    print(f"  {'u m/s':>6} {'time from the start s':>20} {'distance m':>10} {'coil on s':>9} {'of the budget':>13}")
    tab = []
    for u in (1.0, u_sup, 20.0, u_load, 100.0, 300.0, 1000.0):
        tab.append(dict(u=u, t=u / a, x=0.5 * u * u / a, on=t_on(u), share=t_on(u) / b20))
        print(f"  {u:6.1f} {u/a:20.1f} {0.5*u*u/a:10.1f} {t_on(u):9.2f} {t_on(u)/b20:13.2f}")
    print(f"\n  (3) The band between.  From {u_sup:.0f} to {u_load:.0f} m/s ({u_sup/a:.0f} to {u_load/a:.0f} s, the first {0.5*u_load**2/a:.0f} m) neither holds by these two rules: the support cannot be moved that fast through")
    print(f"  its leads and a load coil would be on too long.  Three things soften it, none of them computed: the raised support coils reach a kilometer, so for the first {math.sqrt(2*L_sup/a):.0f} s the")
    print(f"  spacecraft is still over support that is at its cap (how the cap is confined to 100 m with 1 km coils is the open item of round 8); a load coil does not need its full current")
    print(f"  if the support under it still carries part; and the leads' budget is a 20 K rise from one pass, not a limit of the tape.  A run of this stretch would need the support's coils")
    print(f"  resolved along the track, which no script has.")
    out.update(u_sup=u_sup, u_sup210=u_sup210, u_load=u_load, budget20=b20, table=tab, t_leave=math.sqrt(2 * L_sup / a), u_leave=math.sqrt(2 * L_sup * a))
    # a note on the coil's on-time in the launch integrals
    li = r6.launch_profile(fast["W"], fast["F_req"]); t_p, v_p = li["t_p"], li["v_p"]
    s_old = s_new = 0.0
    for i in range(1, len(t_p)):
        u = float(v_p[i]); dt = float(t_p[i] - t_p[i - 1])
        if u <= 0:
            continue
        F, kind = r4.pull_track_of_u(u)
        T_old, L_on = r5.hold_time(u)
        s_old += (2 * 20 / 25.0) * (F / fast["F_req"]) * T_old * u * dt
        s_new += (2 * 20 / 25.0) * (F / fast["F_req"]) * (L_on / u) * u * dt
    print(f"\n  A thing to check, not changed here.  bearing_round4.hold_time and bearing_round5.hold_time give a coil's on-time in a pass as the energized length over (u + V_rel), final_bearing")
    print(f"  C.2's time for a point of the PLATE, which moves with the cable.  The coils are on the casing, where the load moves at u, so a coil's on-time is the energized length over u:")
    print(f"  {r5.hold_time(1000.0)[1]/1000.0*1e3:.0f} ms against {r5.hold_time(1000.0)[0]*1e3:.0f} ms at 1 km/s, {r5.hold_time(30000.0)[1]/30000.0*1e3:.0f} against {r5.hold_time(30000.0)[0]*1e3:.0f} ms at 30 km/s.  Over the heaviest launch the leads' full-current coil-seconds")
    print(f"  are {s_new:.3e} by the casing's frame against {s_old:.3e} as computed ({s_new/s_old:.2f} times), so the leads' pulse heat of 15.2 GJ per launch would be {15.2*s_new/s_old:.0f} GJ and the lead section sized")
    print(f"  from the day's mean heat would grow.  The plate's dwell time (u + V_rel) is right for the plate heat.  [Flagged for a ruling; no earlier number is changed in this round.]")
    out.update(lead_s_old=s_old, lead_s_casing=s_new)
    gt.save("gap_trade_start.json", out)
    return out
