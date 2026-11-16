"""The bearing, round 4 (2026-11-13): Vol_III/Claude_Code_instructions_bearing_round4_2026-11-13.md.
What is left after the Stage B report: the slow coils charged from a DC supply through leads (no transformer), the fast coils'
series turns, how the fast coils are fed and what the supplies and stores weigh, the bias, the plates at the design's resolution,
the cases on the book's 12,000 kg/m casing, the cold width and the deployment.

Run: python bearing_round4.py --part=slow|turns|feed|bias|plates|casing|width|all [--quick]
     (plates: eight plate cases at 2.5 mm, about two hours; casing: six case runs on the 12,000 kg/m casing, about half an hour;
      the deployment runs of the plate candidates are tools/round4_deployment.py; the rest is seconds)

S1  slow    the DC circuit (supply, leads, D2, strip, dump), the lead heat at 551 and 55 turns, the strip's substrate, the
            charging procedure's allowance, the failure table
S2  turns   the series turns of the 25 m fast coil (1 to 20) at 2, 5 and 10 ms: supply voltage and current, insulation, layers,
            lead heat; Paul's 1,500 A case; 12 mm against 3 mm tape
S3  feed    where the supply and the store sit: the heat into the cold per metre of ring idle and over 500 launches a day, the mass
            of the supplies and stores, the family count (the cases are in part casing), the 120 MVA per metre invariant
S4  bias    what the bias costs for each feed and turn count, what it buys (the Stage B compare), a way that costs nothing standing
S5  plates  the candidates (a) to (d) solved as Part 3b solved the book's tracks (2.5 mm cells, the real bundles, the worst tape)
S6  casing  the 25 m coils at 2, 5 and 10 ms and the family sweep (2, 5, 10, 20) on the 12,000 kg/m casing
S7  width   how the cold width is counted; the deployment days of the candidates from tools/round4_deployment.py's logs

Nothing here changes a default of the existing scripts.  Outputs go to runs/ (ORBITAL_RING_OUT).
"""
import contextlib
import io
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_optimize as bo
import bearing_electrical as be
import bearing_design as bd
import bearing_mass as bm
import bearing_layout as bl
import bearing_stageB as sb
import final_bearing as fb
import sled_passage_pulsed as pz
import sled_passage_simulation as sp
import compensation_tracks as ct
import lim_config

MU0 = 4e-7 * math.pi
N_SUP = bl.N_SUP_COILS                         # 10 staggered support coils per track
L_SUP = 1000.0                                 # m, the support coil
L_FAST = 25.0                                  # m, the chosen lift and load coil (Stage B)
TAU_FAST = 2e-3                                # s, the chosen rise
N_FAM = fb.N_FAM                               # 20 staggered families
SERIES_TURNS = (1, 2, 5, 10, 14, 20)
RISES = (2e-3, 5e-3, 10e-3)
FAMILIES = (2, 5, 10, 20)
I_SUPPLY_PAUL = 1500.0                         # A, Paul's case for the fast coils' supply
LAUNCHES_PER_DAY = 500
V_D = sb.V_DIODE_COLD                          # 1.0 V, the free-wheel diode at 70 K (Liao et al. 2024)
T_CLOSE = (1.0, 10.0, 60.0)                    # s, the strip's closing time: its own thermal link (Stage B 2.3), MRI switches (Ogawa et al. 2014)
Q_SITE_BOOK = 14.5e3                           # W, Chapter 14's standing cold load per 500 m site after the deployment (Table 14.2)
Q_SITE_SIX = 17.5e3                            # W, the same with the six tracks' 5.2 m of cold width (Stage B 7.3)
T_MGO = 3e-6                                   # m, the MgO buffer (Chapter 6, Table 6.1: "about 3 um"; the THEVA ISD-MgO paper read gives no thickness)
E_MGO_FILM = (4.4e8, 6.0e8)                    # V/m, the breakdown field of a 200 nm CVD MgO film, two capacitor dots (Talukdar et al. 2018, read in full)
E_MGO_BULK = 0.92e8                            # V/m, pure MgO ceramic (RSC supplementary data, search summary only)
RHO_HAST = bl.RHO_HAST                         # 1.24 uohm m, Hastelloy C-276
T_HAST = 100e-6
K_KAPTON_80 = 0.12                             # W/(m K), polyimide near 80 K [UNVERIFIED: handbook order of magnitude]
CAP_J_PER_KG = 1500.0                          # J/kg, metallized polypropylene pulse capacitors for low pulse duty, 1.5 J/cc (Kerrigan, CDE, CARTS 2007, read in full)
CAP_J_PER_CC = (0.6, 1.3, 1.5, 2.4)            # J/cc: ARL's two 12.5 kJ units 0.64 and 1.3; CDE 1.2 to 1.5; 2.4 for ms-discharge units above 50 kJ (search summaries)
IGBT_MODULE = dict(V=6.5e3, I_pulse=1500.0, kg=1.5, use=0.6)   # a 6.5 kV / 750 A module, its 1 ms surge rating and mass [UNVERIFIED: typical values from memory, not a source read]
J_PULSE_LEAD = math.sqrt(8960 * 200 * 20 / (be.RHO_CU_77 * 0.1))   # A/m2, Part 1c's pulse-sized lead: 20 K rise in a 0.1 s pulse [assumptions as Part 1c]
Q_LEAD_OPT = math.sqrt(be.L0 * (be.T_WARM**2 - be.T_COLD**2))    # W/A, the optimized steady lead, 45.6 W/kA
K_FLUX = 0.33e6                                # N/m per m per track, Part 1's constant-flux stiffness (the e-folding time)
PLATE_CANDIDATES = {                           # support (W, T), lift and load (W, T), metres
    "book": dict(sup=(1.4, 0.20), lift=(1.4, 0.20), note="the book's four plates; Part 3b's lift design is on the same 1.4 m plate"),
    "a": dict(sup=(1.2, 0.15), lift=(1.2, 0.15), note="six plates of 1.2 m by 150 mm (Stage B's choice)"),
    "b": dict(sup=(1.2, 0.10), lift=(1.2, 0.20), note="fast plates 1.2 m by 200 mm, support plates 1.2 m by 100 mm"),
    "c": dict(sup=(1.2, 0.15), lift=(1.2, 0.20), note="fast plates 1.2 m by 200 mm, support plates 1.2 m by 150 mm"),
    "d": dict(sup=(1.2, None), lift=(1.2, None), note="six equal plates with the hardware held at the book's 10,686 kg/m"),
}
T_D = (ring.M_HARDWARE_M - fb.M_LIM_PLATES) / (6 * bd.RHO_PLATE * 1.2)     # 155.6 mm: six equal 1.2 m plates at the book's hardware
PLATE_CANDIDATES["d"]["sup"] = (1.2, T_D); PLATE_CANDIDATES["d"]["lift"] = (1.2, T_D)
COLD_WIDTHS = dict(book=4.0, four_loops=2 * 1.0 + 2 * 0.8, six_loops=2 * 1.0 + 4 * 0.8)


def hr(t=""):
    print(); print("=" * 150)
    if t:
        print("  " + t); print("=" * 150)


def fmt_t(s):
    return bl.fmt_t(s)


def hw_mass(cand):
    (ws, ts), (wl, tl) = cand["sup"], cand["lift"]
    m_s = bd.RHO_PLATE * ws * ts; m_l = bd.RHO_PLATE * wl * tl
    return fb.M_LIM_PLATES + 2 * m_s + 4 * m_l, m_s, m_l


# =============================================================================
# S1  THE SLOW COILS ON DC
# =============================================================================

def turn_variants(c, sc):
    """The 551-turn (one layer) and 55-turn (ten layers) support coil: current, inductances, the strip, at the same ampere-turns."""
    out = []
    for N_tot, layers in ((551, 1), (110, 5), (55, 10)):
        nt = N_tot / 2; k = nt / c["n_t"]
        I0 = c["NI_s0"] / (N_SUP * nt); I_cap = c["NI_cap"] / (N_SUP * nt)
        out.append(dict(N=N_tot, layers=layers, nt=nt, I0=I0, I_cap=I_cap, L_common=sc["L_common"] * k**2, L_leak=sc["L_leak"] * k**2,
                        L_self=sc["L1"] * k**2, V_1pct_s=sc["L_common"] * k**2 * 0.01 * I0))
    return out


def part_slow(quick=False):
    hr("S1  THE SLOW COILS CHARGED FROM A DC SUPPLY THROUGH LEADS (the transformer and D1 gone)")
    r_sup, r_lift = be.load_records(); c = bl.part1_constants(r_sup); sc = bl.switched_circuit(c, L_SUP)
    out = {}
    coils_per_m = 2 * N_SUP / L_SUP                                    # both tracks, 10 staggered 1 km coils each: 0.02 per metre of ring
    print(f"  The circuit now: a DC supply fed from the ring's bus, a pair of leads into the cold, the free-wheel diode D2 across the coil, the strip of bare")
    print(f"  superconductor across D2, and the dump branch (a cold switch into a 20 ohm resistor) across the coil.  Charging: the strip is heated open, the supply")
    print(f"  drives the coil at V less the lead drops; D2 is reverse-biased and idle.  Holding: the strip is cooled closed and carries the whole current; the")
    print(f"  supply's current is brought to zero and the leads carry nothing, or are withdrawn.  D2 is the path the coil's current takes if the supply drops")
    print(f"  while the strip is still open: it free-wheels at the diode drop ({V_D:.1f} V) and drains the coil at V_d / L.")
    # S1.1 lead heat
    print()
    print(f"  S1.1  The heat through the leads.  A conduction-cooled copper lead optimized for a steady current brings {Q_LEAD_OPT*1e3:.1f} W/kA into 70 K (Part 1c, Wiedemann-Franz);")
    print(f"  the same lead with no current leaks {be.lead_heat(0.0, A_over_l=be.lead_heat(1000.0)[1])[0]/be.lead_heat(1000.0)[0]*100:.0f} percent of that.  The lead is sized for the largest current it must carry steadily, the cap (a parked")
    print(f"  spacecraft's takeover); its heat is given at the normal current, at the cap, attached and idle (zero current), and withdrawn (no metal path: nothing).")
    print(f"  Coils: {N_SUP} per track per kilometre, both tracks, {coils_per_m:.2f} per metre of ring, {coils_per_m*500:.0f} per 500 m site.  Chapter 14's standing cold load per site after the deployment:")
    print(f"  {Q_SITE_BOOK/1e3:.1f} kW (Table 14.2); {Q_SITE_SIX/1e3:.1f} kW with the six tracks' cold width (Stage B).")
    print(f"  {'turns':>5} {'layers':>6} {'I normal A':>10} {'I cap A':>7} {'lead pair W: at normal':>22} {'at the cap':>10} {'idle, attached':>14} {'withdrawn':>9} {'per m of ring idle W':>20} {'per site idle W':>15} {'of 14.5 kW':>10} {'charging per site W':>19}")
    variants = turn_variants(c, sc)
    for v in variants:
        q_cap, aol = be.lead_heat(v["I_cap"])
        q_norm = be.lead_heat(v["I0"], A_over_l=aol)[0]; q_idle = be.lead_heat(0.0, A_over_l=aol)[0]
        v.update(q_cap=2 * q_cap, q_norm=2 * q_norm, q_idle=2 * q_idle, aol=aol)
        print(f"  {v['N']:5d} {v['layers']:6d} {v['I0']:10.0f} {v['I_cap']:7.0f} {2*q_norm:22.1f} {2*q_cap:10.1f} {2*q_idle:14.1f} {0.0:9.1f} {2*q_idle*coils_per_m:20.2f} {2*q_idle*coils_per_m*500:15.0f} {2*q_idle*coils_per_m*500/Q_SITE_BOOK*100:9.1f} % {2*q_cap*coils_per_m*500:19.0f}")
    print(f"  (per site idle: every coil's leads attached and carrying nothing; charging per site: all ten coils of the site charging at the cap at once, the")
    print(f"  worst case.  The 110-turn, five-layer coil is added as the middle.)")
    out["leads"] = variants
    v551, v55 = variants[0], variants[-1]
    print(f"  What the leads say about the turn count: at 551 turns the idle leads of a site are {v551['q_idle']*coils_per_m*500:.0f} W, {v551['q_idle']*coils_per_m*500/Q_SITE_BOOK*100:.1f} percent of the standing load; at 55")
    print(f"  turns {v55['q_idle']*coils_per_m*500:.0f} W, {v55['q_idle']*coils_per_m*500/Q_SITE_BOOK*100:.1f} percent, ten times more, because the lead is sized for ten times the current.  With the leads withdrawn after")
    print(f"  each charge neither costs anything between charges, and the charge itself costs {v55['q_cap']*fmt_dur(sc, c)[0]/1e3:.0f} kJ of lead heat at 55 turns for a takeover ({v551['q_cap']*fmt_dur(sc, c)[0]/1e3:.1f} at 551)")
    print(f"  against the {sb_charge_heat(sc, c):.0f} kJ the strip and the diode cost in Stage B's count.  So 55 turns of ten layers stays the balanced choice if the leads are withdrawn or")
    print(f"  disconnected inside the cold between charges (a cold contactor, or the retractable leads of conduction-cooled MRI magnets: GE's US 9,182,464, read in part);")
    print(f"  with leads left attached the 55-turn coil's {v55['q_idle']*coils_per_m*500:.0f} W per site is a cost the 551-turn coil does not have, and 110 turns of five layers ({variants[1]['q_idle']*coils_per_m*500:.0f} W per site)")
    print(f"  keeps most of the defect tolerance at half the lead.")
    # S1.2 the strip's substrate
    print()
    print(f"  S1.2  The strip's substrate.  The strip is a bare superconducting layer in Kapton (Stage B 2.3): {sb.T_REBCO*1e6:.1f} um of REBCO at about 100 uohm cm, no silver, no copper.")
    print(f"  What is under it: the tape's {T_MGO*1e6:.0f} um MgO buffer (Chapter 6, Table 6.1) on a {T_HAST*1e6:.0f} um Hastelloy substrate.  With the strip open and V along it, a conducting")
    print(f"  substrate floats at the strip's mid potential and the buffer at the strip's ends sees V / 2 across {T_MGO*1e6:.0f} um.  If the buffer fails at two places the substrate")
    print(f"  shunts the layer between them: Hastelloy at {RHO_HAST*1e6:.2f} uohm m and {T_HAST*1e6:.0f} um is {RHO_HAST/(T_HAST*0.012):.2f} ohm/m at 12 mm against the layer's {sb.strip_R_per_m(0.012):.0f} ohm/m, {sb.strip_R_per_m(0.012)/(RHO_HAST/(T_HAST*0.012)):.0f} times less, so a")
    print(f"  shunted strip bypasses up to {sb.strip_R_per_m(0.012)/(RHO_HAST/(T_HAST*0.012)):.0f} times the design's 1 percent of the charge and puts that heat into the Hastelloy.  Breakdown fields: MgO thin films")
    print(f"  {E_MGO_FILM[0]/1e8:.1f} to {E_MGO_FILM[1]/1e8:.1f} MV/cm (a 200 nm CVD film, two capacitor dots, the lower 'may be due to a local defect', Talukdar et al. 2018, read in full; 'compare well with")
    print(f"  literature values for MgO films grown by ALD'); bulk MgO ceramic about {E_MGO_BULK/1e8:.1f} MV/cm (search summary).  The ISD-MgO of the tape is a {T_MGO*1e6:.0f} um polycrystalline")
    print(f"  template, not a capacitor dielectric; its own breakdown field was not found [UNVERIFIED], so the film values are the optimistic end.")
    print(f"  {'V charge':>8} {'strip R ohm':>11} {'12 mm strip m':>13} {'V along V/m':>11} {'buffer sees V':>13} {'buffer MV/cm':>12} {'margin to 4.4 MV/cm':>19} {'to bulk 0.9':>11} {'verdict':<60}")
    rows = []
    for V in sb.V_CHARGE:
        R_s = 100.0 * V / sc["I0"]; L12 = R_s / sb.strip_R_per_m(0.012)
        E_b = (V / 2) / T_MGO
        m_film = E_MGO_FILM[0] / E_b; m_bulk = E_MGO_BULK / E_b
        verdict = ("non-conducting substrate or thicker buffer needed" if m_film < 10 else ("ordinary tape without silver and copper will do" if m_bulk > 3 else "ordinary tape, with the buffer tested"))
        rows.append(dict(V=V, R_s=R_s, L12=L12, E_along=V / L12, V_buf=V / 2, E_buf=E_b, m_film=m_film, m_bulk=m_bulk, verdict=verdict))
        print(f"  {V:8g} {R_s:11.0f} {L12:13.1f} {V/L12:11.0f} {V/2:13.0f} {E_b/1e8:12.2f} {m_film:19.1f} {m_bulk:11.1f} {verdict:<60}")
    out["substrate"] = rows
    print(f"  Paul's hand estimate holds: at 1 kV the buffer sees 500 V, {rows[0]['E_buf']/1e8:.1f} MV/cm, a factor {rows[0]['m_film']:.1f} to {E_MGO_FILM[1]/rows[0]['E_buf']:.1f} under the film values and above the bulk ceramic's, on a")
    print(f"  {rows[0]['L12']:.0f} m strip whose buffer must have no pinhole pair: not a design margin.  The way out is Paul's: a non-conducting, non-magnetic substrate (the Kapton")
    print(f"  carrier itself, or a ceramic), or a buffer ten times thicker ({10*T_MGO*1e6:.0f} um brings 1 kV to {rows[0]['E_buf']/10/1e8:.2f} MV/cm), or charging at 100 V, where the buffer sees")
    print(f"  {rows[1]['V_buf']:.0f} V and {rows[1]['E_buf']/1e8:.2f} MV/cm and the tape as made, stripped of its silver and copper, will do.  The 55-turn coil charges at {v55['V_1pct_s']:.0f} V for the same 1 percent")
    print(f"  per second (Stage B 2.5), so at 55 turns the question does not arise.  Kapton in place of the silver and copper: the open strip's {sb.Q_STRIP_MAX*0 + 0.41:.2f} W/cm2 per face")
    print(f"  crosses 25 um of polyimide with a {0.41e4*25e-6/K_KAPTON_80:.1f} K drop [k = {K_KAPTON_80} W/(m K) UNVERIFIED]; what the strip loses with the metal is any stabilizer of its own: a strip that")
    R_s1k = 100.0 * 1000.0 / sc["I0"]
    print(f"  quenches by itself (a bare spot at {sc['I_cap']:.0f} A) puts the coil's current through its rising resistance, {sc['I0']*R_s1k/1e3:.0f} kV across the whole {R_s1k:,.0f} ohm strip at {sc['I0']:.0f} A, so the dump")
    print(f"  switch must close on a volt-level strip voltage within microseconds (a semiconductor switch does), and the strip is sized at twice the cap current in a low field.")
    # S1.3 the procedure and its allowance
    print()
    print(f"  S1.3  The procedure: the supply runs until the coil is just above the current wanted, the power is removed, the strip is cooled.  While the strip cools the")
    print(f"  coil free-wheels through D2 at {V_D:.1f} V and loses V_d t / L.  Ten coils charged together see the common-mode inductance; one coil trimmed alone against nine")
    print(f"  persistent neighbours sees only its leakage inductance (Stage B 2.1), {sc['L_leak']/sc['L_common']*100:.2f} percent of it.  The allowance (the overshoot to charge to) for the closing times")
    print(f"  of Stage B 2.3 (about a second with the strip's own thermal link; tens of seconds to a minute for MRI switches, Ogawa et al. 2014):")
    print(f"  {'turns':>5} {'I normal A':>10} {'L ten together H':>16} {'L one alone H':>13} " + " ".join(f"{'ten: ' + f'{t:g}' + ' s':>13} {'one: ' + f'{t:g}' + ' s':>13}" for t in T_CLOSE))
    allow = []
    for v in variants:
        cells = []
        rr = dict(N=v["N"], ten=[], one=[])
        for t in T_CLOSE:
            d_ten = V_D * t / v["L_common"]; d_one = V_D * t / v["L_leak"]
            rr["ten"].append(d_ten); rr["one"].append(d_one)
            cells.append(f"{d_ten:8.3f} A {d_ten/v['I0']*100:4.2f}%")
            cells.append(f"{d_one:8.2f} A {min(d_one/v['I0']*100, 999):4.0f}%")
        allow.append(rr)
        print(f"  {v['N']:5d} {v['I0']:10.0f} {v['L_common']:16.1f} {v['L_leak']:13.3f} " + " ".join(f"{cc:>13}" for cc in cells))
    out["allowance"] = allow
    print(f"  Ten coils together: the allowance is a fraction of a percent at any closing time, {allow[0]['ten'][-1]/v551['I0']*100:.2f} percent at 551 turns and a minute.  One coil alone: the drain is")
    print(f"  {allow[0]['one'][0]:.2f} A per second at 551 turns ({allow[0]['one'][0]/v551['I0']*100:.1f} percent of its current per second) and {allow[-1]['one'][0]:.0f} A per second at 55 turns ({allow[-1]['one'][0]/v55['I0']*100:.1f} percent per second):")
    print(f"  a lone trim needs the strip closed within a second, or the other order of the standard persistent-magnet procedure: the supply holds the current")
    print(f"  (the coil's voltage is then zero and D2 carries nothing) while the strip cools and closes, and is ramped down afterwards, which moves the")
    print(f"  current into the strip by flux conservation with no allowance at all; the cost is the supply on and the leads carrying the current for the closing time.")
    # the lone trim's voltage and the stuck supply
    print()
    print(f"  S1.4  The lone trim and the stuck supply.  A 1 percent trim of one coil alone is L_leak x 0.01 I = {v551['L_leak']*0.01*v551['I0']:.1f} V s at 551 turns: {v551['L_leak']*0.01*v551['I0']/1000*1e3:.0f} ms at 1 kV, {v551['L_leak']*0.01*v551['I0']/10:.2f} s at 10 V.")
    print(f"  A supply stuck at its full voltage drives a lone coil at V / L_leak: at 1 kV {1000/v551['L_leak']:.0f} A/s at 551 turns (from {v551['I0']:.0f} A to the tape's {bl.IC_12MM:.0f} A limit in {(bl.IC_12MM - v551['I0'])*v551['L_leak']/1000:.1f} s) and")
    print(f"  {1000/v55['L_leak']:,.0f} A/s at 55 turns (to the ten layers' {10*bl.IC_12MM:,.0f} A in {(10*bl.IC_12MM - v55['I0'])*v55['L_leak']/1000:.2f} s); with all ten strips open the rate is V / L_common, {1000/v551['L_common']:.2f} A/s at 551 turns (Stage B's 112 s charge).")
    print(f"  The 1 kV is for the bundle's charge with all ten coils together; a lone trim is done at a few volts, or the supply is current-limited and the dump is")
    print(f"  armed on overcurrent.")
    out["lone"] = dict(Vs_551=v551["L_leak"] * 0.01 * v551["I0"], rate_551=1000 / v551["L_leak"], rate_55=1000 / v55["L_leak"])
    # S1.5 failures
    print()
    print(f"  S1.5  What fails when any one part fails (the DC circuit; the dynamics cases are Stage B's and S6's):")
    items = [("the supply fails off (open)", f"the coil holds through the strip if closed; if the strip is open the coil free-wheels through D2, {V_D/v551['L_common']*1e3:.2f} mA/s with the bundle, {V_D/v551['L_leak']:.2f} A/s alone: close the strip or restore the supply within seconds"),
             ("the supply stuck on at full voltage", f"strip closed: the supply's current goes into the strip until its limit, then the strip quenches and the dump takes the coil; strip open: V / L (S1.4), the supply's current limit and the dump; the interlock: full voltage only with the ten strips open"),
             ("the supply's output shorted", "a short across the coil through the leads; the strip (zero resistance) keeps the current, the leads carry none; the coil cannot be charged until repaired"),
             ("a lead breaks (open)", "no charging; the coil holds through the strip"),
             ("a lead's cold end quenches or warms while charging", f"the lead's drop comes off the charging voltage; at 55 turns a lead at 1 V is {v55['I_cap']:.0f} W into the cold; the supply's voltage and current limits stop it; the coil holds"),
             ("a lead left attached and idle", f"{v551['q_idle']:.1f} W per coil at 551 turns, {v55['q_idle']:.0f} W at 55 (S1.1); withdrawn: nothing"),
             ("D2 open", "the coil has no free-wheel path while the strip is open: the supply must hold until the strip closes, or the dump is closed; the supply's own free-wheel path is the backup"),
             ("D2 shorted", "a short across the coil through the diode's metal (Stage B 6.2): the coil holds, cannot be charged, the dump is bypassed; a quench needs the heaters; repair"),
             ("the strip stuck closed (heater fails)", "the coil stays persistent and cannot be trimmed or dumped; its nine neighbours are trimmed around it (B.1 (f)); on a quench the heaters of A.5"),
             ("the strip stuck open (will not cool)", f"the coil free-wheels through D2: 1 percent in {0.01*v551['I0']*v551['L_common']/V_D/60:.0f} minutes at 551 turns, {0.01*v55['I0']*v55['L_common']/V_D/60:.1f} minutes at 55; the supply holds it at {V_D*v55['I_cap']:.0f} W of lead and diode heat at 55 turns, or the dump moves it to the neighbours"),
             ("the strip quenches by itself", "the coil's current meets the strip's rising resistance: the dump switch closes on a volt-level strip voltage within microseconds (S1.2), 0.3 s to the neighbours; the strip is sized at twice the cap current in a low field"),
             ("the dump switch fails open", "a quench falls back on the heaters of A.5 (the strip cannot be the dump, Stage B 2.4)"),
             ("the dump switch fails closed", f"the 20 ohm resistor sits across the coil: {sc['I0']**2*20/1e3:.0f} kW at first, the coil's current decays into its neighbours in L_leak / R = {sc['L_leak']/20*1e3:.0f} ms and the resistor takes the leakage energy, {0.5*sc['L_leak']*sc['I0']**2/1e3:.1f} kJ; the coil is lost until repair, the bundle keeps its flux"),
             ("the strip's heater fails", "as the strip stuck closed"),
             ("a support coil quenches (whole turn or more)", "the dump at 1 kV into the nine neighbours in 0.3 s at 51 A, 0.7 s at 117 A (Stage B 2.4, B.1)"),
             ("one of the ten coils open", "the nine keep the bundle's flux (A.5): 0.004 percent of force"),
             ("ring power lost", "the coils hold for years through the strips; no supply is needed to hold; the gap is the fast loops' business (S6)")]
    for part_, what in items:
        print(f"    {part_:>52}: {what}")
    out["failures"] = items
    return out


def fmt_dur(sc, c):
    """The takeover's duration at 1 kV (seconds) and the charge from zero, for the lead-heat sentence."""
    t_take = c["n_t"] * 2 * L_SUP * bl.dphi(c, c["NI_s0"], c["NI_cap"]) / 1000.0
    t_zero = c["n_t"] * 2 * L_SUP * bl.dphi(c, 0.0, c["NI_s0"]) / 1000.0
    return t_take, t_zero


def sb_charge_heat(sc, c):
    """Stage B 2.7's cold heat of a takeover at 1 kV, kJ (strip 1 percent bypass plus one diode), for comparison."""
    t_take = fmt_dur(sc, c)[0]
    return (0.01 * 1000 * sc["I_cap"] + V_D * sc["I_cap"]) * t_take / 1e3


# =============================================================================
# S2  THE FAST COILS: SERIES TURNS
# =============================================================================

def fast_coil_basics(r_lift, L=L_FAST, n_fam=N_FAM):
    Phi = 2 * r_lift["W_req"] / r_lift["NI_req"]                       # Wb/m per turn at the requirement (one loop: 2W/NI)
    NI_coil = r_lift["NI_req"] / n_fam                                  # ampere-turns one coil carries (one of n_fam families)
    E_coil = r_lift["W_req"] * L / n_fam
    I_half = 0.5 * r_lift["ic3"]                                        # the worst 3 mm tape at half its limit
    tapes = NI_coil / I_half                                            # tapes through one coil's cross-section
    return dict(Phi=Phi, NI_coil=NI_coil, E_coil=E_coil, I_half=I_half, tapes=tapes, A_fam=tapes * 0.003 * bd.TAPE_T)


def part_turns(quick=False):
    hr("S2  THE FAST COILS: THE SERIES TURNS OF THE 25 M COIL")
    r_sup, r_lift = be.load_records(); theva, theva77 = be.load_theva()
    fc = fast_coil_basics(r_lift)
    out = dict(basics=fc)
    print(f"  The words: 'turns' are the series turns of one coil; the {r_lift['turns3']/2:,.0f} tapes per leg of the lift track (Part 3b, the worst 3 mm tape at {fc['I_half']:.0f} A) are shared among")
    print(f"  the {N_FAM} coils that overlap there, {fc['tapes']:.0f} tapes through one coil's cross-section, carrying {fc['NI_coil']/1e3:.1f} kA of the track's {r_lift['NI_req']/1e3:.0f} kA per leg.  With N series turns")
    print(f"  each turn is a stack of {fc['tapes']:.0f} / N soldered layers, carries {fc['NI_coil']/1e3:.1f} kA / N, and links the whole loop flux: its voltage is Phi' L / tau whatever N, and")
    print(f"  the terminal voltage is N times it.  Stage B's Section 5.3 sentence 'two or three turns in series, 2.4 to 3.6 kV each' divided where it should have")
    print(f"  multiplied: two turns are 7.1 kV each and 14.2 kV at the terminals.  Loop flux at the requirement {fc['Phi']:.4f} Wb/m per turn (2 W / NI, Part 1a), {fc['Phi']*L_FAST:.2f} Wb per 25 m coil;")
    print(f"  stored energy {fc['E_coil']/1e6:.3f} MJ per coil; the supply's rating 2 E / tau whatever N.")
    print()
    print(f"  {'N':>3} {'rise ms':>7} {'I per turn A':>12} {'layers/turn':>11} {'V per turn':>10} {'V supply':>10} {'MVA':>6} {'turn-turn':>9} {'to ground (one end / mid)':>26} {'film at 5 kV/mm: turn / ground mm':>33} {'fits 200x50 at 5 / 10 / 20 kV/mm':>32} {'steady lead pair W':>18} {'pulse lead pair idle W':>22}")
    rows = []
    for tau in RISES:
        V_turn = fc["Phi"] * L_FAST / tau
        for N in SERIES_TURNS:
            I = fc["NI_coil"] / N; V = N * V_turn; VA = V * I
            layers = fc["tapes"] / N
            fits = []
            for ek, E in be.E_DESIGN.items():
                L = be.bundle_layout(0.20, 0.05, N_FAM, fc["A_fam"], N, V_turn, V, E, 10, 2)
                fits.append("yes" if (L and L["fits"]) else "NO")
            q_dc = 2 * be.lead_heat(I)[0]
            A_p = I / J_PULSE_LEAD; q_p0 = 2 * be.lead_heat(0.0, A_over_l=A_p)[0]
            rows.append(dict(N=N, tau=tau, I=I, layers=layers, V_turn=V_turn, V=V, VA=VA, V_tt=V_turn if N > 1 else 0.0, V_g=V, fits=dict(zip(be.E_DESIGN, fits)), q_dc=q_dc, q_p0=q_p0, A_p=A_p))
            print(f"  {N:3d} {tau*1e3:7g} {I:12,.0f} {layers:11.0f} {be.fmt_v(V_turn):>10} {be.fmt_v(V):>10} {VA/1e6:6.0f} {be.fmt_v(V_turn if N > 1 else 0):>9} {be.fmt_v(V) + ' / ' + be.fmt_v(V/2):>26} {(V_turn/5e6*1e3 if N > 1 else 0):14.2f} / {V/5e6*1e3:14.2f} {' / '.join(fits):>32} {q_dc:18,.0f} {q_p0:22.2f}")
    out["rows"] = rows
    print(f"  (turn-turn: neighbouring series turns differ by one turn voltage; to ground: the terminal voltage with one end grounded, half with the supply's midpoint")
    print(f"  grounded; the film thickness is the voltage over the design stress, 5 kV/mm the gas-void case of Part 1b; the fit is Part 1b's bundle layout, 20 family")
    print(f"  cells in 10 columns of 2 with walls of V / E between cells and turns and to ground, the tape at a fill of {be.F_CABLE:.2f} and {be.F_COOL:.0%} coolant; steady lead pair: two")
    print(f"  optimized leads at the turn current, {Q_LEAD_OPT*1e3:.1f} W/kA, the cost if the coil is fed warm and holds a standing current; pulse lead pair idle: Part 1c's lead sized for")
    print(f"  the pulse's adiabatic heating, 1 m long, carrying nothing.)")
    # Paul's 1,500 A case
    N_paul = round(fc["NI_coil"] / I_SUPPLY_PAUL)
    print()
    print(f"  Paul's {I_SUPPLY_PAUL:.0f} A supply: N = {fc['NI_coil']:.0f} / {I_SUPPLY_PAUL:.0f} = {fc['NI_coil']/I_SUPPLY_PAUL:.1f}, so {N_paul} turns of {fc['tapes']/N_paul:.0f} soldered layers at {fc['NI_coil']/N_paul:,.0f} A: "
          + ", ".join(f"{N_paul*fc['Phi']*L_FAST/tau/1e3:.0f} kV at {tau*1e3:g} ms" for tau in RISES) + f" (the hand estimate: 100, 40, 20 kV).")
    print(f"  The LIM coil of Chapter 9 came down to 150 turns because its supply is a continuous inverter at 100 kV; the fast coil's supply is a pulse of 2 E / tau whatever")
    print(f"  the turns, and the turns only trade current for voltage: {N_paul} turns at 2 ms is a 100 kV, 1.5 kA pulse converter, one turn a 7 kV, 21 kA one.  At 100 kV the ground")
    print(f"  wall is {N_paul*fc['Phi']*L_FAST/2e-3/5e6*1e3:.0f} mm at 5 kV/mm and the 200 x 50 mm bundle does not hold it (the table); at 10 ms (20 kV) it does.")
    out["paul"] = dict(N=N_paul, I=fc["NI_coil"] / N_paul, V=[N_paul * fc["Phi"] * L_FAST / tau for tau in RISES])
    # 12 mm against 3 mm tape
    print()
    lim3 = r_lift["lim3"]; d = min(lim3, key=lambda q: min(q["ic_v"], q["ic_h"]))
    face = "v" if d["ic_v"] <= d["ic_h"] else "h"
    ic12 = float(bo.ic_angle(theva, np.array([d["B_" + face]]), np.array([d["th_" + face]]), 0.012)[0])
    tapes12 = fc["NI_coil"] / (0.5 * ic12)
    m3 = fc["tapes"] * 0.003 * bd.TAPE_T * bd.TAPE_RHO; m12 = tapes12 * 0.012 * bd.TAPE_T * bd.TAPE_RHO
    ratio = (0.012 * ic12) / (0.003 * r_lift["ic3"]) * (fc["tapes"] / tapes12) / (fc["tapes"] / fc["tapes"])
    print(f"  12 mm tape in these coils: the worst position of the lift bundle ({d['B_' + face]:.2f} T at {d['th_' + face]:.0f} degrees) gives a 12 mm tape a limit of {ic12:.0f} A against the 3 mm tape's {r_lift['ic3']:.0f} A,")
    print(f"  so one coil needs {tapes12:.0f} layers of 12 mm instead of {fc['tapes']:.0f} of 3 mm, {m12:.2f} against {m3:.2f} kg per metre of coil cross-section (the same tape within {abs(m12/m3-1)*100:.0f} percent).  The pulse loss:")
    print(f"  the magnetization loss of a tape in a field change through its face goes as its width times its limit times the field change (the strip model of Section 7.7),")
    print(f"  so per ampere-turn it goes as the width: {0.012*ic12/(0.003*r_lift['ic3'])*(r_lift['ic3']/ic12):.1f} times for 12 mm against 3 mm here, and {4.5:.1f} times in Part 1d's solve of the support bundles")
    print(f"  (11,350 against 2,530 J per metre of track per pulse, 4.6).  The book's bound for the driven tapes, 605 GJ per launch with the whole field through the face and")
    print(f"  53 GJ with it along the face (Section 7.7), becomes {605*4:,.0f} and {53*4:.0f} GJ with 12 mm tape: 12 mm tape costs four times the pulse loss for no saving in mass, and the")
    print(f"  3 mm tape stays.  What 12 mm would buy is a quarter of the layer count per turn ({tapes12/N_paul:.0f} against {fc['tapes']/N_paul:.0f} at {N_paul} turns) and so a quarter of the solder joints.")
    out["tape12"] = dict(ic12=ic12, tapes12=tapes12, m12=m12, m3=m3, loss_ratio=4.0)
    return out


# =============================================================================
# S3  HOW THE FAST COILS ARE FED, AND WHAT IT WEIGHS
# =============================================================================

def pull_track_of_u(u):
    """One track's pull at sled speed u (N/m), from Part 3b's launch (PASS_PULLS, the pair's pull): the load pair below orbital velocity, the lift
    pair above; linear between the points, zero at the orbital speed, flat below 1 km/s."""
    pts = sorted([(uk * 1e3, F * 1e3 / 2, kind) for uk, F, kind in bo.PASS_PULLS])
    us = [p[0] for p in pts]; Fs = [p[1] for p in pts]
    u_orb = be.V_ORBIT
    xs = [0.0] + [p for p in us if p < u_orb] + [u_orb] + [p for p in us if p > u_orb]
    ys = [Fs[0]] + [F for uu, F in zip(us, Fs) if uu < u_orb] + [0.0] + [F for uu, F in zip(us, Fs) if uu > u_orb]
    return float(np.interp(u, xs, ys)), ("load" if u < u_orb else "lift")


ON_TIME_FRAME = "casing"      # round 9 (2026-11-16): the coils are on the casing, where the load moves at u.  "plate" restores rounds 4 to 8


def hold_time(u, L_coil, frame=None):
    """The time one coil is energized in a pass at speed u: the energized length (final_bearing C.2) over the sled's speed.  The coil is on the
    casing, so the zone passes it at u.  Rounds 4 to 8 divided by u + V_rel, which is the time a point of the PLATE spends in the zone (the plate
    moves with the cable); frame="plate" gives that."""
    L_on = u * max(0.05, bd.L_SC / u) + bd.L_SC + 2 * L_coil + 2 * u * bd.TAU_HW
    return L_on / (u if (frame or ON_TIME_FRAME) == "casing" else u + bd.V_REL)


def launch_integrals(r_lift, L_coil=L_FAST, n_fam=N_FAM, frame=None):
    """Per launch, along the sled's run: the magnetic energy cycled through the pulsing pair's coils (J), and the full-current-equivalent
    lead seconds per metre of ring (sum over the pulsing coils of (I/I_full)^2 x hold time, per metre), from lsm_simulation's profile."""
    t_p, v_p, P_p = pz.lsm_profile()
    F_req = r_lift["F_req"]; W_req = r_lift["W_req"]
    E_cycle = 0.0; lead_s = 0.0; x_run = 0.0
    per_u = {}
    for i in range(1, len(t_p)):
        u = float(v_p[i]); dt = float(t_p[i] - t_p[i - 1]); dx = u * dt
        if u <= 0:
            continue
        F, kind = pull_track_of_u(u)
        E_cycle += 2 * W_req * (F / F_req) * dx                     # both tracks of the pulsing pair, per metre of run
        lead_s += (2 * n_fam / L_coil) * (F / F_req) * hold_time(u, L_coil, frame) * dx   # I^2 ~ F: full-current seconds of the coils per metre
        x_run += dx
    for uk in (1, 3, 5, 10, 20, 30):
        u = uk * 1e3; F, kind = pull_track_of_u(u)
        per_u[uk] = dict(F=F, kind=kind, hold=hold_time(u, L_coil, frame), frac=F / F_req)
    return dict(E_cycle=E_cycle, lead_s=lead_s, x_run=x_run, per_u=per_u)


def part_feed(quick=False):
    hr("S3  HOW THE FAST COILS ARE FED, AND WHAT IT WEIGHS")
    r_sup, r_lift = be.load_records()
    fc = fast_coil_basics(r_lift)
    L_ring = ring.L_RING
    n_m = 4 * N_FAM / L_FAST                                            # coils per metre of ring, four tracks
    n_tot = n_m * L_ring
    E_coil = fc["E_coil"]; VA = 2 * E_coil / TAU_FAST
    out = dict(n_per_m=n_m, n_total=n_tot, E_coil=E_coil, VA=VA)
    print(f"  Four tracks, {N_FAM} staggered families of {L_FAST:.0f} m coils: a coil, a supply and a store every {L_FAST/N_FAM:.2f} m of each track, {n_m:.1f} per metre of ring, {n_tot/1e6:.1f} million")
    print(f"  around the {L_ring/1e3:,.0f} km ring, each {E_coil/1e6:.3f} MJ at full pull and {VA/1e6:.0f} MVA for the {TAU_FAST*1e3:.0f} ms rise.")
    # S3.0 the invariant
    print()
    print(f"  S3.0  The invariant.  The supply rating per metre of track is 2 W / tau with W the track's stored energy per metre at full pull ({r_lift['W_req']/1e3:.1f} kJ/m): {2*r_lift['W_req']/TAU_FAST/1e6:.0f} MVA per")
    print(f"  metre of track at 2 ms, {4*2*r_lift['W_req']/TAU_FAST/1e6:.0f} per metre of ring for four tracks, whatever the coil length, the family count or the series turns; and the store")
    print(f"  per metre of track is W, {r_lift['W_req']/1e3:.1f} kJ/m, {4*r_lift['W_req']/1e3:.0f} kJ per metre of ring.  Check over length, families and turns (MVA per coil x coils per metre):")
    print(f"  {'coil m':>6} {'families':>8} {'turns':>5} {'coils per m of track':>20} {'E per coil MJ':>13} {'MVA per coil':>12} {'MVA per m of track':>18} {'A per coil':>10} {'kV per coil':>11}")
    inv = []
    for L in (12.5, 25.0, 50.0):
        for nf in FAMILIES:
            for N in (1, 14):
                E = r_lift["W_req"] * L / nf; va = 2 * E / TAU_FAST; per_m = nf / L
                I = r_lift["NI_req"] / nf / N; V = N * fc["Phi"] * L / TAU_FAST
                inv.append(dict(L=L, nf=nf, N=N, E=E, va=va, per_m=per_m, I=I, V=V))
                if N == 1 or (L == 25.0 and nf == 20):
                    print(f"  {L:6g} {nf:8d} {N:5d} {per_m:20.2f} {E/1e6:13.3f} {va/1e6:12.0f} {va*per_m/1e6:18.1f} {I:10,.0f} {V/1e3:11.1f}")
    out["invariant"] = inv
    print(f"  Confirmed: {2*r_lift['W_req']/TAU_FAST/1e6:.1f} MVA per metre of track in every row.  What the families and the length do move: the size of one unit (its energy, current")
    print(f"  and rating) and the count, in inverse proportion; what the turns move: the unit's voltage against its current.")
    # S3.1 where the supply and the store sit
    li = launch_integrals(r_lift)
    out["launch"] = {k: v for k, v in li.items()}
    print()
    print(f"  S3.1  Where the supply and the store sit.  The launch (lsm_simulation's profile: {li['x_run']/1e3:,.0f} km of run, {li['x_run']/L_ring:.1f} laps of the ring to 30 km/s, as Part 4's plate heat")
    print(f"  used it) cycles {li['E_cycle']/1e12:.1f} TJ of magnetic energy through the pulsing pair's coils ({li['E_cycle']/(2*r_lift['W_req']*li['x_run'])*100:.0f} percent of the full-pull energy over the run), and the coils'")
    print(f"  full-current-equivalent on-time summed over the ring is {li['lead_s']/1e6:,.1f} million coil-seconds per launch ({li['lead_s']/n_tot*1e3:.0f} ms per coil averaged over all {n_tot/1e6:.0f} million, each coil")
    print(f"  passed once per lap; a lift coil at 30 km/s is on for {li['per_u'][30]['hold']*1e3:.0f} ms at full current).  Per day, {LAUNCHES_PER_DAY} launches.")
    print(f"  Three ways to feed a coil, the heat into the cold per metre of ring, idle and averaged over the day:")
    print(f"  {'turns':>5} {'I full A':>8} | {'(i) warm, steady leads: idle W/m':>32} {'per ring GW':>11} {'+ pulses W/m':>12} | {'(ii) pulse-sized lead: idle W/m':>31} {'+ pulses W/m':>12} {'lead dT K':>9} {'resized for 20 K: mm2':>21} {'idle + pulses W/m':>17} | {'(iii) cold: per % of pulse energy W/m':>37}")
    rows = []
    E_day_m = li["E_cycle"] * LAUNCHES_PER_DAY / 86400 / L_ring       # W per metre of ring of energy cycled, averaged
    duty = li["lead_s"] * LAUNCHES_PER_DAY / 86400 / (n_m * L_ring)   # full-current-equivalent seconds per second, per coil, averaged
    for N in SERIES_TURNS:
        I = fc["NI_coil"] / N
        q_dc, aol = be.lead_heat(I); q_dc_idle = be.lead_heat(0.0, A_over_l=aol)[0]
        idle_i = 2 * q_dc_idle * n_m                                     # W/m, every coil's steady lead pair leaking with no current
        pulse_i = 2 * (q_dc - q_dc_idle) * duty * n_m                     # the extra while carrying current, averaged
        A_p = I / J_PULSE_LEAD; R_p = be.RHO_CU_77 * 1.0 / A_p
        q_p0 = be.lead_heat(0.0, A_over_l=A_p)[0]
        idle_ii = 2 * q_p0 * n_m
        pulse_ii = 2 * I**2 * R_p * duty * n_m                            # the Joule heat of the pulses, all of it to the cold [conservative]
        P_lead_mean = I**2 * R_p * duty                                   # mean Joule power per lead, W
        dT_lead = P_lead_mean * 1.0 / (2 * 400.0 * A_p)                   # a 1 m lead's mid-point rise from its mean Joule heating [1-D, k 400 W/(m K)]
        # the lead resized so that its mean Joule heating lifts its middle 20 K at most: dT ~ 1/A^2 at fixed duty
        A_m = A_p * math.sqrt(dT_lead / 20.0) if dT_lead > 20.0 else A_p
        R_m = be.RHO_CU_77 * 1.0 / A_m
        idle_m = 2 * be.lead_heat(0.0, A_over_l=A_m)[0] * n_m; pulse_m = 2 * I**2 * R_m * duty * n_m
        cold_per_pct = 0.01 * E_day_m
        rows.append(dict(N=N, I=I, idle_i=idle_i, pulse_i=pulse_i, idle_ii=idle_ii, pulse_ii=pulse_ii, dT_lead=dT_lead, A_p=A_p, A_m=A_m, idle_m=idle_m, pulse_m=pulse_m, cold_per_pct=cold_per_pct))
        print(f"  {N:5d} {I:8,.0f} | {idle_i:32,.0f} {idle_i*L_ring/1e9:11.1f} {pulse_i:12.2f} | {idle_ii:31.1f} {pulse_ii:12.1f} {dT_lead:9.0f} {A_m*1e6:21.0f} {idle_m + pulse_m:17.1f} | {cold_per_pct:37.1f}")
    out["feeds"] = rows; out["duty"] = duty
    q_site = Q_SITE_BOOK / 500.0
    print(f"  (Chapter 14's standing cold load is {q_site:.0f} W per metre of ring.  (i): a lead pair optimized for the coil's full current, leaking with no current all day, plus the")
    print(f"  extra while it carries the pulses.  (ii): Part 1c's lead sized for a 20 K adiabatic rise in a 0.1 s pulse, 1 m long, {rows[0]['A_p']*1e6:.0f} mm2 at one turn: its conduction leak, the")
    print(f"  Joule heat of the pulses taken all into the cold, and the rise of its middle from that heat averaged over the day (a copper lead does not cool between")
    print(f"  launches {86400/LAUNCHES_PER_DAY:.0f} s apart); then the same lead made thick enough to hold that rise to 20 K, and its heat.  (iii): the converter and the store at 70 K,")
    print(f"  the loss per pulse a share of the energy cycled [the share is UNVERIFIED: 'a few percent' in Part 1c], refilled through a lead carrying milliamperes at the bus")
    print(f"  voltage, nothing else standing.  The duty: {duty*1e3:.2f} ms of full-current-equivalent per second per coil, {duty*100:.2f} percent.)")
    print(f"  What it says: (i) is ruinous at one turn, {rows[0]['idle_i']*L_ring/1e9:.0f} GW of 70 K heat for the leads alone, {rows[0]['idle_i']/q_site:,.0f} times the standing load, and still {rows[-1]['idle_i']/q_site:.0f} times it at 20 turns;")
    print(f"  (ii) as Part 1c sized it would run {rows[0]['dT_lead']:.0f} K hot at one turn; resized for the day's mean heat it costs {rows[0]['idle_m']+rows[0]['pulse_m']:.0f} W/m at one turn ({(rows[0]['idle_m']+rows[0]['pulse_m'])/q_site:.1f} times the standing")
    print(f"  load), {rows[3]['idle_m']+rows[3]['pulse_m']:.1f} W/m at 10 turns and {rows[-1]['idle_m']+rows[-1]['pulse_m']:.1f} at 20 ({(rows[-1]['idle_m']+rows[-1]['pulse_m'])/q_site*100:.0f} percent), and it cannot hold a standing current (the bias, S4); (iii) costs {rows[0]['cold_per_pct']:.0f} W/m per percent")
    print(f"  of pulse energy lost in the cold, {3*rows[0]['cold_per_pct']:.0f} W/m at 3 percent, {3*rows[0]['cold_per_pct']/q_site:.1f} times the standing load, the same at every turn count, because the loss is on the")
    print(f"  energy cycled and the launch cycles {li['E_cycle']/1e12:.0f} TJ through the coils.  The feed is therefore warm through pulse-sized leads at ten turns and more ((ii)), with")
    print(f"  the converter and the store warm, where their loss costs nothing cold; a cold converter ((iii)) must lose well under a percent per pulse to compete; a")
    print(f"  standing current through warm leads is out at any turn count.")
    # S3.2 the mass
    print()
    E_m = 4 * r_lift["W_req"]                                           # J per metre of ring stored at full pull, four tracks
    print(f"  S3.2  The mass of the supplies and the stores per metre of casing.  The store must hold the coil's full-pull energy and give it back: {E_m/1e3:.0f} kJ per metre of ring")
    print(f"  ({E_coil/1e6:.3f} MJ per coil), whatever the families or the length.  Pulse capacitors: metallized polypropylene for low pulse duty reaches 1.2 to 1.5 J/cc and about")
    print(f"  {CAP_J_PER_KG/1e3:.1f} kJ/kg (Kerrigan, CDE, CARTS 2007, read in full: a resin-sealed flat element at 75 percent of the film's breakdown stress, for hundreds to thousands of")
    print(f"  pulses); long-life units sit lower, 0.6 to 1.3 J/cc for two 12.5 kJ units (ARL, search summary) and up to 2.4 J/cc for millisecond-discharge units above 50 kJ")
    print(f"  (IEEE PPC 2015, search summary); the ring's coils fire {LAUNCHES_PER_DAY} times a day, {LAUNCHES_PER_DAY*365/1e6:.2f} million times a year, so the long-life end applies [the kg per litre of those")
    print(f"  units is UNVERIFIED; CDE's 1.5 J/cc at 1.5 kJ/kg is 1 kg/L and is used to convert].")
    print(f"  {'store':>44} {'J/kg':>6} {'kg/m of ring':>12} {'of the casing (12,000)':>22}")
    mass_rows = []
    for lab, jcc in (("long-life PP film, 0.6 J/cc", 0.6), ("PP film, 1.3 J/cc (ARL's GA unit)", 1.3), ("low-duty PP film, 1.5 J/cc (CDE, read in full)", 1.5), ("ms-discharge >50 kJ, 2.4 J/cc", 2.4)):
        jkg = jcc * 1e3
        m = E_m / jkg
        mass_rows.append(dict(label=lab, J_per_kg=jkg, kg_per_m=m))
        print(f"  {lab:>44} {jkg:6.0f} {m:12.0f} {m/ring.M_LOAD_M*100:21.1f} %")
    out["store_mass"] = mass_rows
    # the switches
    print(f"  The switches: a pulse converter per coil that blocks the terminal voltage and passes the turn current for 2 ms; two switch positions at least (a half")
    print(f"  bridge) with modules in series for the voltage and in parallel for the current.  With a {IGBT_MODULE['V']/1e3:.1f} kV module used at {IGBT_MODULE['use']:.0%} of its voltage, {IGBT_MODULE['I_pulse']:.0f} A of 1 ms surge and")
    print(f"  {IGBT_MODULE['kg']:.1f} kg [UNVERIFIED: typical IGBT module ratings from memory, not a source read]:")
    print(f"  {'turns':>5} {'V supply kV':>11} {'I A':>8} {'series':>6} {'parallel':>8} {'modules per coil':>16} {'kg per coil':>11} {'kg/m of ring':>12}")
    sw_rows = []
    for N in SERIES_TURNS:
        I = fc["NI_coil"] / N; V = N * fc["Phi"] * L_FAST / TAU_FAST
        n_s = math.ceil(V / (IGBT_MODULE["use"] * IGBT_MODULE["V"])); n_p = math.ceil(I / IGBT_MODULE["I_pulse"])
        n_mod = 2 * n_s * n_p; m = n_mod * IGBT_MODULE["kg"]
        sw_rows.append(dict(N=N, V=V, I=I, n_s=n_s, n_p=n_p, n_mod=n_mod, kg=m, kg_per_m=m * n_m))
        print(f"  {N:5d} {V/1e3:11.1f} {I:8,.0f} {n_s:6d} {n_p:8d} {n_mod:16d} {m:11.0f} {m*n_m:12.0f}")
    out["switch_mass"] = sw_rows
    print(f"  So the supplies and stores weigh of the order of {min(r['kg_per_m'] for r in mass_rows) + min(r['kg_per_m'] for r in sw_rows):,.0f} to {max(r['kg_per_m'] for r in mass_rows) + max(r['kg_per_m'] for r in sw_rows):,.0f} kg per metre of casing, the store the larger part at the long-life end: comparable")
    print(f"  to the bearing's own {933:.0f} kg/m, which Stage B's Section 6.4 left them out of, and {min(r['kg_per_m'] for r in mass_rows)/ring.M_LOAD_M*100:.0f} to {(max(r['kg_per_m'] for r in mass_rows) + max(r['kg_per_m'] for r in sw_rows))/ring.M_LOAD_M*100:.0f} percent of the 12,000 kg/m.  The mass does not depend on the family count or the coil")
    print(f"  length (the invariant); it depends on the store's specific energy and on nothing else the layout can change.  A store that is cold adds its cryostat; a store shared")
    print(f"  by a site's coils ({E_m*500/1e6:.0f} MJ per 500 m site) weighs the same and needs a 2 ms distribution to {n_m*500:.0f} coils.")
    # S3.3 the families
    print()
    print(f"  S3.3  The number of families at {L_FAST:.0f} m.  Per supply: the current (one turn), the energy, the rating; the count; a dead coil as a share of the track's pull.  The")
    print(f"  cases are run in S6 (part casing) on the 12,000 kg/m casing.")
    print(f"  {'families':>8} {'coils per m of track':>20} {'A per coil (one turn)':>21} {'E per coil MJ':>13} {'MVA per coil':>12} {'supplies per track':>18} {'supplies, four tracks':>21} {'dead coil, share of pull':>24}")
    fam_rows = []
    for nf in FAMILIES:
        I = r_lift["NI_req"] / nf; E = r_lift["W_req"] * L_FAST / nf; va = 2 * E / TAU_FAST; n_tr = nf * L_ring / L_FAST
        fam_rows.append(dict(nf=nf, I=I, E=E, va=va, n_track=n_tr, n_all=4 * n_tr, dead=1.0 / nf))
        print(f"  {nf:8d} {nf/L_FAST:20.2f} {I:21,.0f} {E/1e6:13.3f} {va/1e6:12.0f} {n_tr/1e6:15.1f} million {4*n_tr/1e6:17.1f} million {1/nf*100:23.0f} %")
    out["families"] = fam_rows
    return out


# =============================================================================
# S4  THE BIAS
# =============================================================================

def part_bias(quick=False):
    hr("S4  THE BIAS: WHAT IT COSTS FOR EACH FEED AND TURN COUNT, WHAT IT BUYS, AND A WAY THAT COSTS NOTHING STANDING")
    r_sup, r_lift = be.load_records(); fc = fast_coil_basics(r_lift)
    L_ring = ring.L_RING; n_m = 4 * N_FAM / L_FAST
    NI_bias = r_lift["NI_req"] * math.sqrt(sb.BIAS / r_lift["F_req"])            # 95 kA per leg for 10 kN/m per face
    out = dict(NI_bias=NI_bias)
    print(f"  The {sb.BIAS/1e3:.0f} kN/m bias is {NI_bias/1e3:.1f} kA per leg in every fast track, {NI_bias/N_FAM/1e3:.2f} kA per coil at one turn, standing.  Its cost by the feed of S3:")
    print(f"  {'turns':>5} {'bias A per coil':>15} {'(i) steady leads W per coil':>27} {'W/m':>8} {'per ring GW':>11} {'(ii) pulse lead':>15} {'(iii) cold switch 1 mohm W/coil':>31} {'W/m':>8} {'cold IGBT 1.5 V W/coil':>22} {'W/m':>8} {'persistent strip':>16}")
    rows = []
    for N in SERIES_TURNS:
        I_b = NI_bias / N_FAM / N
        q_i = 2 * be.lead_heat(I_b)[0]
        q_sw = I_b**2 * 1e-3; q_igbt = 1.5 * I_b
        rows.append(dict(N=N, I_b=I_b, q_leads=q_i, q_sw=q_sw, q_igbt=q_igbt))
        print(f"  {N:5d} {I_b:15,.0f} {q_i:27,.0f} {q_i*n_m:8,.0f} {q_i*n_m*L_ring/1e9:11.1f} {'cannot hold it':>15} {q_sw:31.2f} {q_sw*n_m:8.2f} {q_igbt:22,.0f} {q_igbt*n_m:8,.0f} {'0':>16}")
    out["cost"] = rows
    q_site = Q_SITE_BOOK / 500.0
    print(f"  (Paul's hand estimate at one turn, 0.43 kW per coil and 61 GW, is the first row; the cold switch at 1 milliohm is the on-resistance of parallel silicon dies in")
    print(f"  liquid nitrogen (Stage B of the layout report); 1.5 V is an IGBT's on-state drop [UNVERIFIED: typical]; the persistent strip is a superconducting switch across")
    print(f"  the coil, as the slow coils have.)  Against the {q_site:.0f} W/m standing load: warm leads carry the bias at {rows[0]['q_leads']*n_m/q_site:.0f} times the standing load at one turn and {rows[-1]['q_leads']*n_m/q_site:.1f} times")
    print(f"  at 20 turns; a cold IGBT free-wheeling it costs {rows[0]['q_igbt']*n_m/q_site:.0f} times the standing load at one turn and {rows[-1]['q_igbt']*n_m/q_site:.0f} times at 20; a cold 1 milliohm switch {rows[0]['q_sw']*n_m/q_site:,.0f} times")
    print(f"  at one turn (I^2 R at 4.75 kA is {rows[0]['q_sw']/1e3:.1f} kW per coil) and {rows[-1]['q_sw']*n_m/q_site:.1f} times at 20; only a superconducting path across the coil holds it for nothing.")
    # what the bias buys
    print()
    print(f"  S4.2  What the bias buys.  The point and beam models already command the current that gives the pull they want (I = sqrt(F / F_ref), the square-root")
    print(f"  command of the brief, through the first-order current lag of the rise): the 'no bias' rows of Stage B are that case.  Stage B's rerun of the 114 cases with")
    print(f"  the bias against without (both with the support frozen, 25 m coils at 2 ms, the on-top casing) moved 4 cases by more than 0.5 mm: the pull 50 ms late +0.6 mm")
    print(f"  (35.7 to 36.3, the smallest single-failure gap), the spread pull with no sled (beam) +0.6, one lift track alone at 18 km/s +1.9, the false full lift pull")
    print(f"  of 50 ms -0.7; the planned passes within 0.1 mm.  The bias buys 0.6 mm on the worst single failure and nothing on the passes: with a square-root command")
    print(f"  a small pull from zero is not slow (the current for 1 percent of the pull is 10 percent of the full current, reached in a tenth of the rise), and the")
    print(f"  authority dF/dNI = 2F/NI being zero at zero is a property of the plant the controller inverts, not a limit on it.")
    out["buys"] = dict(worst_gain_mm=0.6, pass_change_mm=0.1, moved=4)
    # a way that costs nothing standing
    print()
    tape_bias = 2 * NI_bias / (0.5 * bl.IC_12MM) * 0.012 * bd.TAPE_T * bd.TAPE_RHO    # a persistent 12 mm winding at half its limit, two legs
    print(f"  S4.3  If a bias is wanted, a way that costs nothing standing: a separate persistent winding on the same loop, charged once like the slow coils (DC supply,")
    print(f"  withdrawn leads, a strip across it), {NI_bias/1e3:.0f} kA per leg of 12 mm tape at half its limit: {tape_bias:.1f} kg/m of tape per track, {NI_bias/(0.5*bl.IC_12MM):,.0f} turns in its two legs,")
    print(f"  in ten staggered 1 km coils.  The driven turn then works from a field that is already there, with the authority of the total ampere-turns.  Its cost is the")
    print(f"  coupling Part 1d found for a persistent winding under a driven one: the 1 km persistent coil keeps its flux integrated over its length, so a 25 m")
    print(f"  pulse of dNI in the driven coil takes {L_FAST/L_SUP*100:.1f} percent of dNI out of the persistent coil along its whole kilometre, and gives it back after the pulse;")
    print(f"  the pull lost under the pulse is {L_FAST/L_SUP*100:.1f} percent of the drive's and the bias dips by the same share along the kilometre.  [SPECULATIVE] A second way:")
    print(f"  a strip across each driven coil with a limit just above the bias current, so that a command pushes the strip past its limit and opens it in")
    print(f"  microseconds, the converter drives the coil, and the strip closes again seconds after the event; the standing cost is zero but two events within")
    print(f"  the strip's closing time find it open, and a self-triggered switch of this kind was not found in the literature read.")
    print(f"  The recommendation: no bias.  The square-root command makes it unnecessary, the models show it buys 0.6 mm, and every way of carrying it that is not")
    print(f"  superconducting costs more than the whole cryogenic plant's standing load.")
    out["persistent_bias"] = dict(tape_per_track=tape_bias, coupling_share=L_FAST / L_SUP)
    return out


# =============================================================================
# S5  THE PLATES AT THE DESIGN'S RESOLUTION
# =============================================================================

@contextlib.contextmanager
def plate_patch(W, T):
    """Point bearing_electrical's Track (and bearing_optimize's plate box) at a plate of width W and thickness T."""
    old = be.PLATE_W, be.PLATE_T, bo.PLATE_W, bo.PLATE_T
    be.PLATE_W, be.PLATE_T, bo.PLATE_W, bo.PLATE_T = W, T, W, T
    try:
        yield
    finally:
        be.PLATE_W, be.PLATE_T, bo.PLATE_W, bo.PLATE_T = old


def solve_plate_case(kind, W, T, h, theva, theva77, r, F_req, w_tape, ks=(1.25, 1.5, 2.0), verbose=True):
    """One track on a plate W x T, Part 3b's method: the ampere-turns for the requirement (two nonlinear solves and a log-log step), the worst
    tape there, the turns at the 50 percent share and the tape mass; then the worst tape at ks x NI_req for the rating (the share of the worst
    tape's limit at the design turn count reaching 70 percent) and the cap (100 percent), and the pull curve to 2 x NI_req."""
    t0 = time.time()
    geom = bo.BEST_SUPPORT if kind == "support" else bo.BEST_LIFT
    I_ref = bo.I_REF_12 if w_tape == 0.012 else bo.I_REF_3
    with plate_patch(W, T):
        tr = be.Track(geom, h, w_tape=w_tape, I_ref=I_ref, NI_ref=r["NI_req"])
        NI1 = r["NI_req"]
        g1, s1 = tr.solve(tr.NI_design(NI1))
        NI2 = NI1 * math.sqrt(F_req / s1["F"])
        g2, s2 = tr.solve(tr.NI_design(NI2))
        NI_req = math.exp(math.log(NI1) + (math.log(F_req) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI1)) / (math.log(s2["F"]) - math.log(s1["F"])))
        g3, s3 = tr.solve(tr.NI_design(NI_req))
        wt, lim = be.worst_tape(g3, tr, tr.NI_design(NI_req), w_tape, theva, theva77)
        ic_req = wt[0]; I50 = 0.5 * ic_req
        n_legs = len(tr.legs); turns_leg = NI_req / I50
        tape = n_legs * turns_leg * w_tape * bd.TAPE_T * bd.TAPE_RHO
        curve = [(NI_req, s3["F"], s3["sat"], s3["B_plate"], ic_req)]
        for k in ks:
            NIk = k * NI_req
            gk, sk = tr.solve(tr.NI_design(NIk))
            wtk, _ = be.worst_tape(gk, tr, tr.NI_design(NIk), w_tape, theva, theva77)
            curve.append((NIk, sk["F"], sk["sat"], sk["B_plate"], wtk[0]))
    NIs = np.array([c[0] for c in curve]); Fs = np.array([c[1] for c in curve]); ics = np.array([c[4] for c in curve])
    share = (NIs / turns_leg) / ics                                     # the worst tape's share of its limit at each level
    def at_share(s_target):
        if share.max() < s_target:
            return float("nan"), float("nan")
        NI_s = float(np.interp(s_target, share, NIs))
        F_s = math.exp(float(np.interp(math.log(NI_s), np.log(NIs), np.log(Fs))))
        return NI_s, F_s
    NI70, F70 = at_share(0.70); NI100, F100 = at_share(1.0)
    res = dict(kind=kind, W=W, T=T, h=h, NI_req=NI_req, F_req=F_req, F_solved=s3["F"], ic_req=ic_req, B=wt[1], th=wt[2], turns_leg=turns_leg, n_legs=n_legs, tape=tape,
               sat=s3["sat"], B_plate=s3["B_plate"], W_mag=s3["W"], plate_m=bd.RHO_PLATE * W * T, curve=[list(c) for c in curve], share=share.tolist(),
               NI70=NI70, F70=F70, NI100=NI100, F100=F100, F_2x=float(Fs[-1]), sat_2x=float(curve[-1][2]), seconds=time.time() - t0)
    if verbose:
        print(f"    {kind:>7} {W:.1f} m x {T*1e3:.0f} mm: {NI_req/1e3:.0f} kA for {F_req/1e3:.1f} kN/m (solved {s3['F']/1e3:.1f}); worst {w_tape*1e3:.0f} mm tape {wt[1]:.2f} T at {wt[2]:.0f} deg, {ic_req:.0f} A; {turns_leg:,.0f} turns per leg,")
        print(f"            {tape:.1f} kg/m; plate {s3['sat']*100:.1f} % saturated, {s3['B_plate']:.2f} T inside, {bd.RHO_PLATE*W*T:,.0f} kg/m; rating (70 %) {NI70/1e3:.0f} kA, {F70/1e3:.0f} kN/m; cap (100 %) {NI100/1e3:.0f} kA, {F100/1e3:.0f} kN/m;")
        print(f"            at 2 x NI_req {Fs[-1]/1e3:.0f} kN/m, {curve[-1][2]*100:.0f} % saturated ({time.time()-t0:.0f} s)", flush=True)
    return res


def part_plates(quick=False):
    hr("S5  THE PLATES AT THE DESIGN'S RESOLUTION (Part 3b's method: 2.5 mm cells, the real bundles, the worst tape)")
    r_sup, r_lift = be.load_records(); theva, theva77 = be.load_theva()
    h = 0.005 if quick else 0.0025
    print(f"  Cells {h*1e3:g} mm.  Support track: Part 3b's two 0.5 m loops in four 150 x 67 mm bundles, 12 mm tape, {r_sup['F_req']/1e3:.1f} kN/m; lift and load tracks: one 0.8 m loop in two")
    print(f"  200 x 50 mm bundles, 3 mm tape, {r_lift['F_req']/1e3:.1f} kN/m.  Part 3b's records on the 1.4 m x 200 mm plate: support {r_sup['NI_req']/1e3:.0f} kA, worst tape {r_sup['ic12']:.0f} A, {r_sup['tape12']:.1f} kg/m,")
    print(f"  rating (70 %) {bo.rated_tape(r_sup, 0.7)/1e3:.0f} kN/m at {r_sup['rate_tape']['70'][0]/1e3:.0f} kA; lift {r_lift['NI_req']/1e3:.0f} kA, {r_lift['ic3']:.0f} A, {r_lift['tape3']:.1f} kg/m, rating {bo.rated_tape(r_lift, 0.7)/1e3:.0f} kN/m at {r_lift['rate_tape']['70'][0]/1e3:.0f} kA (the book's 288), the 1 T rule's")
    print(f"  {r_lift['F_rated']/1e3:.0f} kN/m at {r_lift['NI_rated']/1e3:.0f} kA (the brief's 242).  Note: Part 3b solved the lift track on the fixed 1.4 m plate (rule 7 of the October brief); the Stage B report")
    print(f"  called the lift design's plate 1.2 m, which was the air-core table's row, not Part 3b's.  Both are solved here.")
    cases = {}
    for key, cand in PLATE_CANDIDATES.items():
        for kind, (W, T) in (("support", cand["sup"]), ("lift", cand["lift"])):
            cases.setdefault((kind, round(W, 3), round(T, 4)), [])
            cases[(kind, round(W, 3), round(T, 4))].append(key)
    if quick:
        cases = {k: v for k, v in cases.items() if "a" in v or "book" in v}
    print(f"  {len(cases)} plate cases to solve (each about {6 if not quick else 2} solves):")
    results = {}
    for (kind, W, T), keys in cases.items():
        r, F_req, w_tape = (r_sup, r_sup["F_req"], 0.012) if kind == "support" else (r_lift, r_lift["F_req"], 0.003)
        res = solve_plate_case(kind, W, T, h, theva, theva77, r, F_req, w_tape, ks=((1.5,) if quick else (1.25, 1.5, 2.0)))
        res["candidates"] = keys
        results[f"{kind}|{W}|{T}"] = res
    # the candidates' table
    print()
    print(f"  The candidates (plate and hardware mass with the LIM plates' {fb.M_LIM_PLATES:,.0f} kg/m; tape per track at the 50 percent share):")
    print(f"  {'cand':>4} {'support plate':>14} {'kA':>4} {'tape kg/m':>9} {'sat %':>5} {'rating kN/m':>11} {'cap kN/m':>8} {'at 2x kN/m':>10} | {'fast plate':>14} {'kA':>4} {'tape kg/m':>9} {'sat %':>5} {'rating':>7} {'cap':>5} {'at 2x':>6} | {'plates kg/m':>11} {'hardware kg/m':>13} {'casing tape kg/m':>16}")
    table = []
    for key, cand in PLATE_CANDIDATES.items():
        rs = results.get(f"support|{round(cand['sup'][0],3)}|{round(cand['sup'][1],4)}"); rl = results.get(f"lift|{round(cand['lift'][0],3)}|{round(cand['lift'][1],4)}")
        if rs is None or rl is None:
            continue
        hw, m_s, m_l = hw_mass(cand)
        tape_casing = 2 * rs["tape"] + 2 * rl["tape"] * math.sqrt(fb.F_LOAD_TRACK / r_lift["F_req"]) + 2 * rl["tape"]
        table.append(dict(cand=key, hw=hw, plates=hw - fb.M_LIM_PLATES, m_sup=m_s, m_lift=m_l, sup=rs, lift=rl, tape_casing=tape_casing))
        print(f"  {key:>4} {f'{cand['sup'][0]:.1f} m x {cand['sup'][1]*1e3:.0f} mm':>14} {rs['NI_req']/1e3:4.0f} {rs['tape']:9.1f} {rs['sat']*100:5.1f} {rs['F70']/1e3:11.0f} {rs['F100']/1e3:8.0f} {rs['F_2x']/1e3:10.0f} | {f'{cand['lift'][0]:.1f} m x {cand['lift'][1]*1e3:.0f} mm':>14} {rl['NI_req']/1e3:4.0f} {rl['tape']:9.1f} {rl['sat']*100:5.1f} {rl['F70']/1e3:7.0f} {rl['F100']/1e3:5.0f} {rl['F_2x']/1e3:6.0f} | {hw - fb.M_LIM_PLATES:11,.0f} {hw:13,.0f} {tape_casing:16.0f}")
    print(f"  (rating: the pull at which the worst tape of the design's turn count reaches 70 percent of its limit, Part 3b's rule; cap: 100 percent; at 2x: the pull at")
    print(f"  twice the requirement's ampere-turns, the plate's plateau; casing tape: two support, two load (the lift design at 198 kN/m) and two lift tracks.)")
    ra = [t for t in table if t["cand"] == "a"]
    if ra:
        rl = ra[0]["lift"]
        print(f"  What (a) loses on the fast tracks: on the 1.2 m x 150 mm plate the lift track reaches its {r_lift['F_req']/1e3:.1f} kN/m at {rl['NI_req']/1e3:.0f} kA with the plate {rl['sat']*100:.0f} percent saturated; its")
        print(f"  rating by the tape rule is {rl['F70']/1e3:.0f} kN/m and its cap {rl['F100']/1e3:.0f}, against {bo.rated_tape(r_lift, 0.7)/1e3:.0f} and more on the book's plate; at twice the ampere-turns the plate gives {rl['F_2x']/1e3:.0f} kN/m.  The cases")
        print(f"  cap a lift track at its {fb.F_LIFT_TRACK/1e3:.1f} kN/m requirement, so nothing in the 114 cases moves; what is lost is the headroom above it: a heavier or faster launch, a")
        print(f"  track carrying more than its half of the pair.  The air-core table's 220 kN/m (Stage B 7.1) was that model's plateau, not the tape rating.")
    return dict(results=results, table=table, h=h)


# =============================================================================
# S6  THE CASES ON THE BOOK'S CASING
# =============================================================================

def inside_A():
    """final_bearing's Part A with the bearing counted inside the 12,000 kg/m casing (the book's), the cable of the present design point."""
    theva = bd.read_theva(); sup, lift = fb.load_air()
    with contextlib.redirect_stdout(io.StringIO()):
        A = fb.part_A(theva, sup, lift)
    A2 = dict(A); A2["top"] = A["inside"]; A2["m_top"] = ring.M_LOAD_M; A2["dm_top"] = A["dm_in"]
    return A2, theva, sup, lift


def part_casing(quick=False):
    hr("S6  THE CASES ON THE BOOK'S CASING: 12,000 kg/m, THE BEARING INSIDE IT")
    A2, theva, sup, lift = inside_A()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    d = A2["top"]
    print(f"  Casing {A2['m_top']:,.0f} kg/m (the book: the bearing's {A2['dm_top']['casing_side']:.0f} kg/m counted inside it), weight {A2['m_top']*fb.G_NET/1e3:.1f} kN/m, cable {d['m_cable']:,.0f} kg/m with the book's")
    print(f"  {ring.M_HARDWARE_M:,.0f} kg/m of hardware, tension {d['tension']/1e9:.1f} GN, support requirement {A2['dm_top']['F_s']/1e3:.1f} kN/m per track.  Stage B's runs used final_bearing's on-top casing:")
    print(f"  12,681 kg/m, cable 108,980 kg/m, tension 723 GN, requirement 248 kN/m (the Stage B report's Section 5.2 said 103,328 for those runs; the logs say 108,980).")
    print(f"  The six-plate hardware (10,372 kg/m) would make the cable 103,328 kg/m; the cable mass enters the dynamics through the cable's share of the")
    print(f"  pull and the ring matching (0.2 percent): not rerun.  Lift and load coils {L_FAST:.0f} m, the support frozen, the protection in every lane.")
    L0 = fb.L_LIFT, fb.L_LOAD; nf0 = fb.N_FAM, bd.N_FAM_F
    runs = {}
    configs = [(f"{tau*1e3:g} ms, 20 families", tau, 20) for tau in RISES] + [(f"2 ms, {nf} families", TAU_FAST, nf) for nf in FAMILIES if nf != 20]
    if quick:
        configs = [configs[0], configs[-1]]
    for name, tau, nf in configs:
        fb.L_LIFT, fb.L_LOAD = L_FAST, L_FAST; fb.N_FAM = nf; bd.N_FAM_F = nf
        t0 = time.time()
        fname = "bearing_round4_casing_" + name.replace(" ms, ", "ms_").replace(" families", "fam").replace(" ", "") + ".txt"
        try:
            gaps, records, txt = bl.run_cases(tau, True, A2, theva, sup, lift, tabs, fname, prot_all=True, head=f"[round 4, the 12,000 kg/m casing: {name}]")
        finally:
            fb.L_LIFT, fb.L_LOAD = L0; fb.N_FAM, bd.N_FAM_F = nf0
        runs[name] = dict(gaps=gaps, summary=bl.parse_summary(txt), txt=txt, tau=tau, nf=nf, seconds=time.time() - t0, log=fname)
        s = runs[name]["summary"]
        print(f"  ran {name}: pair {s['e1']} / {s['e30']} mm at 1 / 30 km/s, one track {s['o30']}, smallest single-failure gap {s['gap']} mm ({s['case'][:44]}), {time.time()-t0:.0f} s -> {fname}", flush=True)
    print()
    print(f"  {'run':<22} {'1 km/s pair mm':>14} {'30 km/s pair mm':>15} {'30 km/s one mm':>14} {'smallest single-failure gap mm':>30} {'in':<40} {'single contacts':>15} {'band margin mm':>14} {'dead coil, beam mm':>18} {'dead coil, one lift track mm':>28}")
    out = dict(casing=A2["m_top"], cable=d["m_cable"], tension=d["tension"], F_s=A2["dm_top"]["F_s"], runs={})
    def gap_of(gaps, sub):
        hits = [g for g in gaps if sub in g[0]]
        return min(min(g[2], g[3]) for g in hits) * 1e3 if hits else float("nan")
    for name, r in runs.items():
        s = r["summary"]; singles = [g for g in r["gaps"] if bl.is_single(g[0])]; cs = [g[0] for g in singles if g[4]]
        try:
            e1 = float(s["e1"]); e30 = float(s["e30"]); o30 = float(s["o30"]); gap = float(s["gap"])
        except ValueError:
            e1 = e30 = o30 = gap = float("nan")
        dead = gap_of(r["gaps"], "one coil dead, beam"); dead1 = gap_of(r["gaps"], "one coil dead, one lift track alone")
        out["runs"][name] = dict(e1=e1, e30=e30, o30=o30, gap=gap, case=s["case"], contacts=cs, margin=20 - max(e1, e30), dead=dead, dead_one=dead1, tau=r["tau"], nf=r["nf"], log=r["log"], seconds=r["seconds"])
        print(f"  {name:<22} {e1:14.2f} {e30:15.2f} {o30:14.2f} {gap:30.1f} {s['case'][:40]:<40} {len(cs):15d} {20-max(e1,e30):14.1f} {dead:18.1f} {dead1:28.1f}")
    # Table 7.5's rows and the failure list from the 2 ms, 20 families run
    base_name = configs[0][0]; base = runs[base_name]
    print()
    print(f"  Table 7.5's rows ({base_name}): the F.1 block of the log:")
    lines = base["txt"].splitlines()
    i0 = next(i for i, l in enumerate(lines) if "F.1  Sled passage" in l)
    for l in lines[i0:i0 + 9]:
        print("  " + l)
    i1 = next(i for i, l in enumerate(lines) if "fastest launch inside 20 mm" in l)
    print("  " + lines[i1].strip())
    out["table75"] = lines[i0:i0 + 9]; out["fastest_one_lift"] = lines[i1].strip()
    print()
    print(f"  The failure list ({base_name}; the smallest gap of the case, mm; C = contact):")
    gaps = base["gaps"]
    def g_of(sub):
        hits = [g for g in gaps if sub in g[0]]
        return "; ".join(f"{min(g[2], g[3])*1e3:.1f}{' C' if g[4] else ''}" for g in hits[:3]) if hits else "n/a"
    items = [("a driven coil's supply (one of the 20 families): one coil dead, beam / one lift track alone", g_of("one coil dead, beam") + " / " + g_of("one coil dead, one lift track alone")),
             ("one lift track lost at 30 km/s", g_of("one lift track lost at 30")), ("one load track lost: the one-track pass at 1 km/s", g_of("sled passage 1 km/s, one track")),
             ("one support track lost, both windings, at rest / as the zone arrives", g_of("1. one support track loses both windings at rest") + " / " + g_of("1. one support track loses both windings as the zone")),
             ("both support tracks lost as the zone arrives / at rest (load pair catches)", g_of("both support tracks lost as the zone arrives") + " / " + g_of("both support tracks lost, load pair catches")),
             ("a gap sensor +5 / -5 mm", g_of("sensor error +5 mm") + " / " + g_of("sensor error -5 mm")),
             ("a false full lift pull 50 ms / a false load pull 50 ms", g_of("false full lift pull 50 ms") + " / " + g_of("false full load pull 50 ms")),
             ("a pull 50 ms late / early", g_of("pull 50 ms late") + " / " + g_of("pull 50 ms early")), ("the spread pull with no sled, beam", g_of("spread pull with no sled, beam")),
             ("the roll with one lift track alone (per-side regulation)", g_of("roll: one lift track alone, per-side"))]
    for part_, what in items:
        print(f"    {part_:>88}: {what}")
    out["failures"] = items
    tau_e = math.sqrt(A2["m_top"] / (2 * K_FLUX))
    print()
    print(f"  The e-folding time of the attractive bearing with the support frozen at constant flux: sqrt(m / (2 k)) with k = {K_FLUX/1e6:.2f} kN/m per mm per track (Part 1's")
    print(f"  constant-flux stiffness) and m = {A2['m_top']:,.0f} kg/m: {tau_e*1e3:.0f} ms (139 ms on the 12,681 kg/m casing).")
    out["efold_ms"] = tau_e * 1e3
    return out


# =============================================================================
# S7  THE COLD WIDTH AND THE DEPLOYMENT
# =============================================================================

def part_width(quick=False):
    hr("S7  THE COLD WIDTH: HOW IT IS COUNTED, AND THE DEPLOYMENT DAYS OF THE CANDIDATES")
    r_sup, r_lift = be.load_records()
    print(f"  lim_config counts the cold width the cable radiates onto as LEV_COIL_WIDTH = D_LEV_COIL x N_LEV_PLATES = {lim_config.D_LEV_COIL:.1f} m x {lim_config.N_LEV_PLATES} = {lim_config.LEV_COIL_WIDTH:.1f} m: one metre of coil per")
    print(f"  track, four tracks, and the heat onto the coils is {lim_config.LEV_COIL_MLI_HEAT_FLUX_REF:.0f} W/m2 of that width at the reference cable temperature, scaled with T^4 (lim_physics.calc_levitation_coil_heat):")
    print(f"  {lim_config.LEV_COIL_MLI_HEAT_FLUX_REF*lim_config.LEV_COIL_WIDTH*500/1e3:.0f} kW per 500 m site at the reference, Table 14.2's 10.0 kW.  Other ways to count the same tracks:")
    w_sup_loop = 2 * bo.BEST_SUPPORT["W"]; w_lift_loop = bo.BEST_LIFT["W"]
    wb_sup = bo.shape_width(bo.BEST_SUPPORT["shape"]); wb_lift = bo.shape_width(bo.BEST_LIFT["shape"])
    rows = [("the book's count: 1 m of coil per track", 4 * 1.0, 6 * 1.0),
            ("the loops' widths (support 2 x 0.5 m, lift and load 0.8 m)", 2 * w_sup_loop + 2 * w_lift_loop, 2 * w_sup_loop + 4 * w_lift_loop),
            ("the tracks' widths, loop plus one bundle (support 1.3 m, lift 1.0 m)", 2 * r_sup["track_width"] + 2 * r_lift["track_width"], 2 * r_sup["track_width"] + 4 * r_lift["track_width"]),
            ("the bundles' faces alone (support 4 x 150 mm, lift 2 x 200 mm)", 2 * 4 * wb_sup + 2 * 2 * wb_lift, 2 * 4 * wb_sup + 4 * 2 * wb_lift),
            ("the plates (book 1.4 m; candidates 1.2 m)", 4 * 1.4, 6 * 1.2)]
    print(f"  {'count':>70} {'four tracks m':>13} {'six tracks m':>12} {'kW per site, four':>17} {'six':>6}")
    for lab, w4, w6 in rows:
        print(f"  {lab:>70} {w4:13.1f} {w6:12.1f} {lim_config.LEV_COIL_MLI_HEAT_FLUX_REF*w4*500/1e3:17.1f} {lim_config.LEV_COIL_MLI_HEAT_FLUX_REF*w6*500/1e3:6.1f}")
    print(f"  Stage B's 5.2 m is the loops' count for six tracks ({w_sup_loop:.1f} + {w_sup_loop:.1f} + 4 x {w_lift_loop:.1f}); the book's four tracks counted the same way are {2*w_sup_loop + 2*w_lift_loop:.1f} m, not 4.0, so 0.4 m of the")
    print(f"  5.2 - 4.0 = 1.2 m step is the counting and 0.8 m is the two load tracks.  What the cable sees is the cryostats' top faces, which are wider than the bundles")
    print(f"  (the loop's inside is former, cold if it is inside the cryostat) and narrower than the plates; the loops' count is kept as the design's measure, and the")
    print(f"  table shows what the others would give.  The deployment runs (tools/round4_deployment.py, lim_config's design run with each candidate's hardware and the")
    print(f"  loops' cold width) are read below if their logs exist.")
    out = dict(counts=[dict(label=l, four=a, six=b) for l, a, b in rows])
    base = os.path.join(ring.out_dir(), "round4_hw")
    print()
    print(f"  {'run':>34} {'hardware kg/m':>13} {'cable kg/m':>10} {'cold width m':>12} {'days':>6} {'cryo peak MW':>12} {'cold load peak kW':>17} {'LIM work EJ':>11}")
    deps = {}
    for tag in ("book_4.0", "book_3.6", "a", "b", "c", "d"):
        p = os.path.join(base, tag, "lim_design", "rising_radius_lim_simulation.txt")
        meta_p = os.path.join(base, tag, "meta.json")
        if not os.path.exists(p):
            continue
        txt = open(p, encoding="utf-8", errors="replace").read()
        meta = json.load(open(meta_p, encoding="utf-8")) if os.path.exists(meta_p) else {}
        def grab(pat, conv=float):
            import re
            m = re.search(pat, txt)
            return conv(m.group(1).replace(",", "")) if m else float("nan")
        days = grab(r"End of the run at ([\d.,]+) days"); cryo = grab(r"p_cryo_max\s+([\d.]+) MW"); qc = grab(r"q_cryo_cold_max\s+([\d.]+) kW")
        cable = grab(r"Cable mass M_CABLE_M=([\d.]+) kg/m"); work = grab(r"([\d.]+) EJ")
        deps[tag] = dict(days=days, cryo_peak_MW=cryo, q_cold_peak_kW=qc, cable=cable, hw=meta.get("hw"), cold=meta.get("cold_width"), work_EJ=work)
        print(f"  {tag:>34} {meta.get('hw', float('nan')):13,.0f} {cable:10,.0f} {meta.get('cold_width', float('nan')):12.1f} {days:6.1f} {cryo:12.2f} {qc:17.1f} {work:11.3f}")
    if not deps:
        print(f"  (no deployment logs yet under {base})")
    out["deployment"] = deps
    return out


def main():
    part = "all"; quick = False
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--part="):
            part = arg.split("=")[1]
        elif arg == "--quick":
            quick = True
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    hr("THE BEARING, ROUND 4 (2026-11-13)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km, {ring.QUADRATIC} quadratic: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m with hardware {ring.M_HARDWARE_M:,.0f}, side {ring.CABLE_SIDE:.2f} m{' (quick)' if quick else ''}.")
    res = {}
    cache = os.path.join(ring.out_dir(), "bearing_round4.json" if part == "all" else f"bearing_round4_{part}.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            res = json.load(f)
    def save():
        with open(cache, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (o.tolist() if isinstance(o, np.ndarray) else (bool(o) if isinstance(o, np.bool_) else str(o))))
    for name, fn in PARTS.items():
        if part in (name, "all"):
            res[name] = fn(quick); save()
    hr(f"DONE ({time.time()-t_all:.0f} s)")


PARTS = {"slow": part_slow, "turns": part_turns, "feed": part_feed, "bias": part_bias, "plates": part_plates, "casing": part_casing, "width": part_width}

if __name__ == "__main__":
    main()
