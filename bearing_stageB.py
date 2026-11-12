"""The bearing's coils, Stage B (2026-11-12, evening): Vol_III/Claude_Code_instructions_bearing_stageB_2026-11-12.md.
Slow coils (the persistent support coils charged through a transformer, a cold diode and a superconducting switch strip, no
flux pumps and no leads), fast coils (the lift and load compensation coils, their length and rise, a standing bias), the layout
with the load compensation coils on their own tracks and plates, the whole bearing's cases and failures, the plates and the
design point.

Run: python bearing_stageB.py --part=slow|fast|bias|layout|plates|whole|gap|table74|heat|all [--quick]
     (the default quadratic is the rising one, the present cable; the fast-coil sweep and the whole-bearing cases take an hour
      and more each: run them in the background)

S1  slow     the support coil's charging circuit: settling times L I / V, the takeover and the recovery, the switch strip's
             resistance, heat, opening and closing, the series-turn choice (551 / 55 / 6 turns), the protection by the soldered
             layers and by the strip and the dump, the losses in the cold, what the slow coils can do, the coil length
S3  fast     final_bearing's 114 cases with the support frozen and the protection in every lane for lift and load coils of 12.5,
             25 and 50 m at rises of 2, 5, 10 and 20 ms; the turn voltage and the supply rating of each coil length
S3  bias     a standing current in the load windings against a standing pull of the lift tracks, and the chosen combination rerun
S2  layout   layout (A): the load compensation coils on their own tracks and plates beside the support tracks: the coupling
             against the spacing, the load track's design, the plates and the fit on the 7.39 m underside; layout (B) far enough
             to compare: the single-loop load winding on the outer bundles solved with the plate, the roll moment, the residual
             coupling with an imperfect symmetry
S6  plates   the lightest plate for each kind of track (overhang, thickness, a thickness that varies), the pull, the worst tape and
             the tape mass, the launch heat in the narrower plates
S4  whole    the 114 cases of the chosen combination on the present cable, the support frozen and the slow coils' switches able to
             act between events, the failures, the hold with no power, the mass on the casing
P2  gap      the support track at 75 and 50 mm (Parts 2 of the electrical brief)
P3  table74  the fields in tesla for Table 7.4 (Part 3)
P4  heat     the plate heat of a launch on the present tracks (Part 4)

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
import final_bearing as fb
import sled_passage_pulsed as pz
import sled_passage_simulation as sp
import compensation_tracks as ct
import compensation_fields as cf
import magnetostatic_2d as ms
import levitation_config as lc

MU0 = 4e-7 * math.pi
N_SUP = bl.N_SUP_COILS
LENGTHS_FAST = ((50.0, 25.0), (50.0, 50.0), (25.0, 25.0), (12.5, 12.5))     # (lift, load) coil lengths, m; the first is the book's
RISES = (2e-3, 5e-3, 10e-3, 20e-3)
V_CHARGE = (1000.0, 100.0, 10.0)                                          # V, the charging voltages the brief asks about
V_DIODE_RT = 0.7                                                          # V, a silicon diode's forward drop at room temperature (the brief's figure)
V_DIODE_COLD = 1.0                                                        # V, at 70 K: 10 to 80 percent higher than at 300 K (Liao et al., iEnergy 2024, read in full)
RHO_N_REBCO = (81e-8, 100e-8, 177e-8)                                     # ohm m, the normal-state resistivity of the superconducting layer near 100 K:
                                                                          # 81 to 177 uohm cm on six films (Lyatti et al., Materials 17 (2024) 558, as the brief
                                                                          # cites it; the arXiv text gives 208 uohm cm for a YBCO film and 307 to 419 for the Tb-doped
                                                                          # ones), about 100 uohm cm at 100 K (Cantoni et al., cond-mat/0101384, read in full)
T_REBCO = 3.5e-6                                                          # m, the superconducting layer (Chapter 6: 1 to 3.5 um)
T_KAPTON = 25e-6; RHO_KAPTON = 1420.0; CP_KAPTON_80 = 500.0               # the strip's coating [cp of polyimide near 80 K, UNVERIFIED]
RHO_REBCO_D = 6300.0; CP_REBCO_100 = 200.0                                # kg/m3 and J/(kg K) near 100 K [UNVERIFIED: handbook order of magnitude]
T_C = 92.0; T_OP = 70.0
Q_STRIP_MAX = 2.5                                                         # W/cm2 per face, the transient heat flux Cantoni et al. design a stabilizer for,
                                                                          # "well below the critical heat flux of liquid nitrogen"
R_CONTACT_SOLDER = 30e-9 * 1e-4                                           # ohm m2, a soldered lap joint (Stage A's figure, abstract only)
R_CONTACT_SOLDER_HI = 450e-9 * 1e-4                                       # ohm m2, a fully solder-potted coil (Fleiter and Ballarino, cited by Tokamak Energy, WAMHTS-5 slides)
F_LIFT_ONE = fb.F_LIFT_TRACK; F_LOAD_ONE = fb.F_LOAD_TRACK
Y_PLATE_CENTER = fb.Y_TRACK                                               # the support plates' centres from the cable's centreline (the plates at the cable edges)


def hr(t=""):
    print()
    print("=" * 150)
    if t:
        print("  " + t)
        print("=" * 150)


def fmt_t(s):
    return bl.fmt_t(s)


# =============================================================================
# S1  THE SLOW COILS
# =============================================================================

def strip_R_per_m(w, rho=RHO_N_REBCO[1]):
    """Normal-state resistance per metre of a bare superconducting layer of width w (no silver, no copper)."""
    return rho / (T_REBCO * w)


def strip_heat_capacity_per_m(w):
    """J/(m K) of a bare layer of width w with Kapton on both faces, near 80 to 100 K."""
    return w * (T_REBCO * RHO_REBCO_D * CP_REBCO_100 + 2 * T_KAPTON * RHO_KAPTON * CP_KAPTON_80)


def part_slow(quick=False):
    hr("S1  THE SLOW COILS: CHARGING, HOLDING AND PROTECTING A SUPPORT COIL WITHOUT PUMPS OR LEADS")
    r_sup, r_lift = be.load_records(); c = bl.part1_constants(r_sup)
    sc = bl.switched_circuit(c, 1000.0); L_c = sc["L_common"]; L_leak = sc["L_leak"]; n_t = c["n_t"]
    I0 = sc["I0"]; I_cap = sc["I_cap"]
    NI_one = bo.NI_of_F(r_sup, pz.W)                                   # one track carrying the whole casing (the other lost): 95 kN/m
    I_one = NI_one / (N_SUP * n_t)
    E0 = 0.5 * L_c * I0**2
    out = dict(L_common=L_c, L_leak=L_leak, I0=I0, I_cap=I_cap, I_one=I_one, NI_one=NI_one, E0=E0)
    print(f"  The 1 km coil of the book: {2*n_t:.0f} turns (both loops in series), ten coils per bundle charged together, so each sees the common-mode inductance")
    print(f"  {L_c:,.0f} H (the brief's 2,190 H counts the self inductance ten times; the leakage share of A.5 takes 2.4 percent off), {E0/1e6:.2f} MJ at {I0:.1f} A (141 kA per track);")
    print(f"  {I_cap:.1f} A at the 322 kA cap under a parked spacecraft, {I_one:.1f} A at {NI_one/1e3:.0f} kA when one track carries the whole casing ({pz.W/1e3:.1f} kN/m).")
    # (1) settling times
    print()
    print(f"  S1.1  Settling time L I / V (the saturating flux per ampere-turn of Part 1d integrated; the linear figure in brackets):")
    print(f"  {'change':>44} " + " ".join(f"{'at ' + f'{V:g}' + ' V':>18}" for V in V_CHARGE))
    rows = []
    for lab, NI1, NI2 in (("charge from zero to 141 kA (51 A)", 0.0, c["NI_s0"]), ("takeover of a parked spacecraft, 141 to 322 kA", c["NI_s0"], c["NI_cap"]),
                          (f"recovery after a lost track, 141 to {NI_one/1e3:.0f} kA", c["NI_s0"], NI_one), ("a 1 percent trim", c["NI_s0"], 1.01 * c["NI_s0"])):
        Vs = n_t * 2 * 1000.0 * bl.dphi(c, NI1, NI2)                  # volt-seconds per coil (ten coils share the bundle's flux change: each coil's n_t turns see its own share... see note)
        Vs_lin = L_c * (NI2 - NI1) / (N_SUP * n_t)
        ts = [Vs / V for V in V_CHARGE]
        rows.append(dict(label=lab, Vs=Vs, Vs_lin=Vs_lin, t=ts))
        print(f"  {lab:>44} " + " ".join(f"{fmt_t(t):>10} ({fmt_t(Vs_lin / V)})" for t, V in zip(ts, V_CHARGE)) + f"   {Vs:,.0f} V s per coil")
    print(f"  (The volt-seconds are per coil: all ten coils of a bundle are charged together, each from its own transformer; a coil charged alone against nine")
    print(f"  persistent neighbours sees only the leakage inductance, {L_leak:.2f} H, and takes {L_leak*I0/1000:.2f} s to 51 A at 1 kV, as B.1 found.)")
    out["settling"] = rows
    # (2) the diode and the transformer
    print()
    print(f"  S1.2  The circuit.  A single diode in series with the coil cannot rectify: the coil's current is continuous, so in the half-cycle the diode")
    print(f"  blocks, the current has to flow somewhere, and if that somewhere is the strip (hot, resistive) or the source itself, the reverse half-cycle's")
    print(f"  volt-seconds come straight out of the coil: no net gain.  The circuit that works is the half-wave rectifier with a free-wheel path:")
    print(f"    D1 in series with the transformer's secondary turns (lets the current in on the forward half-cycle);")
    print(f"    D2 across the coil (the free-wheel: on the reverse half-cycle the coil's current circulates through it at the diode drop);")
    print(f"    the strip across the coil, bridging D2: superconducting, it shorts D2 and the coil is persistent with no drop; heated, it is a resistor")
    print(f"    much larger than V / I and the coil charges through D1 at the secondary's whole voltage less the diode drop.")
    print(f"  The diode drop at 70 K: {V_DIODE_COLD:.1f} V taken here (Liao et al., iEnergy 2024, read in full: 10 to 80 percent above the room-temperature value for Si")
    print(f"  and SiC diodes; the brief's 0.7 V is the warm figure).  With the strip still open and the source off, D2 drains the coil at V_d / L = {V_DIODE_COLD/L_c*1e3:.2f} mA/s")
    print(f"  ({V_DIODE_RT/L_c*1e3:.2f} at 0.7 V): 1 percent ({0.01*I0:.2f} A) in {fmt_t(0.01*I0*L_c/V_DIODE_COLD)}, so the strip must close within minutes of the end of a charge.")
    print(f"  The transformer: the secondary is a few turns of the coil's own tape passing through a core whose primary is warm (the core crosses the")
    print(f"  cryostat wall); V_rms = 4.44 f N_s B A.  Core area for N_s = 5 turns at 1.5 T:")
    print(f"  {'f Hz':>6} " + " ".join(f"{f'{V:g} V':>12}" for V in V_CHARGE))
    for f in (50.0, 400.0, 1000.0):
        print(f"  {f:6.0f} " + " ".join(f"{V/(4.44*f*5*1.5)*1e4:9.0f} cm2" for V in V_CHARGE))
    # AC loss in the secondary turns: the current is a half-wave pulse of the coil's current at f; Norris strip loss per cycle
    def norris_strip(I_pk, Ic):
        i = min(I_pk / Ic, 0.999)
        return MU0 * Ic**2 / math.pi * ((1 - i) * math.log(1 - i) + (1 + i) * math.log(1 + i) - i**2)
    Ic_tape = bl.IC_12MM
    print(f"  AC loss in the secondary turns (the coil's own tape, Ic {Ic_tape:.0f} A, carrying the coil current as half-wave pulses; Norris strip loss per cycle, an upper")
    print(f"  bound for a stack): per turn of {2*1000:.0f} m, 5 turns:")
    for f in (50.0, 400.0):
        for I in (I0, I_cap):
            q = norris_strip(I, Ic_tape)
            print(f"    {f:4.0f} Hz at {I:5.1f} A: {q*1e6:7.2f} uJ/m per cycle, {q*f*2000*5:8.1f} W in the five turns while charging")
    # (3) the strip
    print()
    print(f"  S1.3  The switch strip: a bare superconducting layer ({T_REBCO*1e6:.1f} um, no silver or copper) in Kapton, across the coil.")
    print(f"  Normal-state resistivity near 100 K: 81 to 177 uohm cm (six films, Lyatti et al. 2024, as the brief cites it; the group's arXiv text gives 208 uohm cm for a YBCO")
    print(f"  film and 307 to 419 for Tb-doped ones), about 100 uohm cm (Cantoni et al. 2001, read in full); nothing for THEVA's GdBCO.  Resistance per metre:")
    print(f"  {'width':>8} " + " ".join(f"{f'{r*1e8:.0f} uohm cm':>14}" for r in RHO_N_REBCO))
    for w in (0.003, 0.012):
        print(f"  {w*1e3:5.0f} mm " + " ".join(f"{strip_R_per_m(w, r):10.0f} ohm/m" for r in RHO_N_REBCO))
    print(f"  The strip's resistance when open must be much more than V / I, or it bypasses the charge: for a 1 percent bypass R_s = 100 V / I.  Its heat is")
    print(f"  then 1 percent of the charging power whatever the voltage.  What the voltage sets is the strip's length and the field along it.")
    print(f"  {'V charge':>9} {'R_s ohm':>8} {'3 mm strip m':>12} {'12 mm strip m':>13} {'heat W':>7} {'E V/m (3 mm)':>12} {'W/cm2 per face':>14} {'open energy J (3 mm)':>20} {'current I_cap':>13}")
    strips = []
    for V in V_CHARGE:
        R_s = 100.0 * V / I0
        L3 = R_s / strip_R_per_m(0.003); L12 = R_s / strip_R_per_m(0.012)
        P = V**2 / R_s; i_b = V / R_s
        E_field = V / L3; q = E_field * i_b / (2 * 0.3 * 100)          # W per cm2 per face: E (V/m) x i (A) over 2 faces x 0.3 cm x 100 cm/m
        E_open = strip_heat_capacity_per_m(0.003) * L3 * (100.0 - T_OP)
        strips.append(dict(V=V, R_s=R_s, L3=L3, L12=L12, P=P, E_field=E_field, q=q, E_open=E_open))
        print(f"  {V:9g} {R_s:8.0f} {L3:12.1f} {L12:13.1f} {P:7.0f} {E_field:12.0f} {q:14.2f} {E_open:20.0f} {I_cap:13.0f} A")
    print(f"  (open energy: the layer and two 25 um Kapton faces from 70 to 100 K, no former [cp values UNVERIFIED]; the brief's 600 J includes a plate.  The field along")
    print(f"  the strip is what the brief's 100 V/m is: at a 1 percent bypass the heat flux is {strips[0]['q']:.2f} W/cm2 per face, below the {Q_STRIP_MAX:.1f} W/cm2 Cantoni et al. take as safe")
    print(f"  in liquid nitrogen.)  Closed, the strip carries the coil's whole current: {I_cap:.0f} A at the cap, so a 3 mm strip (about 60 A in the lift")
    print(f"  track's field, {r_lift['ic3']:.0f} A in the records; more in a low field) is not enough: a 12 mm strip, or four 3 mm strips in parallel, placed in a low field.")
    print(f"  Opening: a heater raises the strip above {T_C:.0f} K; with {strips[0]['E_open']:.0f} J of heat capacity a 100 W heater opens it in {strips[0]['E_open']/100:.1f} s, and its own dissipation")
    print(f"  ({strips[0]['P']:.0f} W at 1 kV) keeps it open.  Closing: the strip must cool below the temperature at which it carries the persistent current with")
    print(f"  margin (about {T_OP + 5:.0f} K at half its limit); the time is its heat capacity over the thermal link to the bath.  A link that carries the open-state")
    print(f"  dissipation at a 30 K rise is {strips[0]['P']/30:.0f} W/K, and cools the {strips[0]['E_open']:.0f} J in about {strips[0]['E_open']/(strips[0]['P']/30)/30*math.log(30/5):.1f} s; MRI persistent switches close in tens of")
    print(f"  seconds to a minute (the Mini-RT REBCO floating coil's heater switch: 'a switching time of around 1 min', Ogawa et al. 2014, read in full).")
    out["strips"] = strips
    # the strip as the quench dump
    print()
    print(f"  S1.4  The strip as the dump.  Opening the strip on a quench puts the coil's current through its open resistance: V = I R_s, and the leakage")
    print(f"  inductance against the nine persistent neighbours gives the transfer time L_leak / R_s:")
    print(f"  {'R_s ohm':>8} {'V at 51 A':>10} {'V at 117 A':>11} {'L_leak / R_s':>13}")
    for st in strips:
        print(f"  {st['R_s']:8.0f} {I0*st['R_s']/1e3:8.1f} kV {I_cap*st['R_s']/1e3:9.1f} kV {L_leak/st['R_s']*1e3:10.1f} ms")
    R_d = 1000.0 / I0; R_d_cap = 1000.0 / I_cap
    print(f"  A strip sized for a 1 kV charge dumps at 51 kV: not a dump but a flashover.  A dump at the 1 kV the insulation is designed for needs {R_d:.0f} ohm")
    print(f"  at 51 A ({R_d_cap:.1f} ohm at 117 A) and takes {L_leak*I0/1000*1e3:.0f} ms (B.1 (d)), and a {R_d:.0f} ohm strip would bypass a 1 kV charge entirely.  The two jobs want")
    print(f"  resistances a hundred times apart.  Three ways out:")
    print(f"    (i) charge at 10 V through a {R_d:.0f} ohm strip (a 0.5 A bypass, 5 W; {fmt_t(rows[0]['t'][2])} to charge, {fmt_t(rows[1]['t'][2])} for the takeover) and dump at 1 kV through the same strip;")
    print(f"    (ii) charge at 1 kV through a 2 kohm strip and dump through a separate path: a cold semiconductor switch (B.1's) into a {R_d:.0f} ohm resistor, closed on detection;")
    print(f"    (iii) charge at 100 V through a 200 ohm strip, dump through it at 10 kV: the Part 1b insulation at 5 kV/mm and 2 mm between turns would take it, but")
    print(f"         the strip then sees 10 kV along 2 m of 3 mm layer, 5 kV/m, and the energy of the dump in a layer of micrograms.")
    m_strip_hast = 20.0 * 0.012 * 100e-6 * 8890.0
    print(f"  The dump's energy is the leakage energy, {0.5*L_leak*I0**2/1e3:.1f} kJ at 51 A and {0.5*L_leak*I_cap**2/1e3:.1f} kJ at 117 A, and a bare layer cannot take it: a {R_d:.0f} ohm resistor of")
    print(f"  REBCO on its {100:.0f} um Hastelloy substrate (1.0 ohm/m at 12 mm, {R_d/1.03:.0f} m, {m_strip_hast*1e3:.0f} g) would rise {0.5*L_leak*I0**2/(m_strip_hast*250):.0f} K at 51 A and {0.5*L_leak*I_cap**2/(m_strip_hast*250):.0f} K at 117 A [cp 250 J/(kg K),")
    print(f"  UNVERIFIED]; the resistor is a metal one of a few kilograms.  So: the strip is the persistent switch, the dump is a switched resistor, and (ii) is the")
    print(f"  circuit: transformer, D1, D2, the strip, and a dump branch; the strip is heated open at the same moment the dump closes (its {strips[0]['E_open']/100:.1f} s at 100 W is")
    print(f"  inside the {bl.hot_spot_time(I0):.0f} s a bare spot has at 51 A, {bl.hot_spot_time(I_cap):.1f} s at 117 A), and the dump's {L_leak*I0/1000*1e3:.0f} ms does the rest.")
    out["dump"] = dict(R_d=R_d, t_dump=L_leak * I0 / 1000.0, E_leak=0.5 * L_leak * I0**2)
    # (4) the series turns
    print()
    print(f"  S1.5  The number of series turns (the power for a 1 percent per second trim is the same, {L_c*0.01*I0*I0/1e3:.0f} kW per coil, for any count):")
    print(f"  {'turns per coil':>14} {'layers per turn':>15} {'tape A':>7} {'L common H':>10} {'V for 1 %/s':>11} {'diodes W':>9} {'strip R_s (1 %)':>15} {'strip Ic needed':>15} {'dump ohm':>8} {'turn V at 1 kV':>14} {'spot bypass':>11}")
    turns_rows = []
    for N_tot, layers in ((551, 1), (55, 10), (6, 100)):
        nt = N_tot / 2; I = c["NI_s0"] / (N_SUP * nt); Ic = c["NI_cap"] / (N_SUP * nt)
        Lc = L_c * (nt / n_t)**2; V1 = Lc * 0.01 * I; Pd = 2 * V_DIODE_COLD * I
        R_s = 100.0 * V1 / I; R_dump = 1000.0 / Ic
        turns_rows.append(dict(N=N_tot, layers=layers, I=I, L=Lc, V1=V1, Pd=Pd, R_s=R_s, R_dump=R_dump))
        print(f"  {N_tot:14d} {layers:15d} {I:7.0f} {Lc:10.1f} {V1:11.0f} {Pd:9.0f} {R_s:15.1f} {Ic:12.0f} A {R_dump:8.2f} {1000/N_tot:14.1f} {'none' if layers == 1 else 'yes':>11}")
    print(f"  (tape A: the current per turn at 141 kA; strip Ic needed: the persistent current at the cap; diodes: the two drops at the charging current, only while")
    print(f"  charging.  The turn voltage at a 1 kV charge is 1 kV over the turns.)")
    print(f"  What changes with fewer turns: the diode loss grows in proportion to the current (but only while charging, kilojoules per charge); the strip and")
    print(f"  the dump must carry it (at 5,100 A, 22 widths of 12 mm tape in parallel); the transformer's secondary carries it; the charging voltage falls to the")
    print(f"  tens of volts, where the insulation is nothing; and a turn of many soldered layers carries its current around a defect in one layer (S1.6),")
    print(f"  which a single-tape turn cannot.  The strip's film area and heat are the same at any count.  The balanced choice is the middle row: about 55")
    print(f"  turns of 10 soldered layers at {turns_rows[1]['I']:.0f} A, {turns_rows[1]['V1']:.0f} V for a 1 percent per second trim, {turns_rows[1]['Pd']:.0f} W in the diodes while charging, a {turns_rows[1]['R_s']:.0f} ohm strip and a")
    print(f"  {turns_rows[1]['R_dump']:.1f} ohm dump; the 551-turn coil has no defect tolerance and the 6-turn coil asks the strip and the diodes for kiloamperes.")
    out["turns"] = turns_rows
    # (5) protection
    print()
    r70 = bl.tape_resistance_per_m(T_OP)
    print(f"  S1.6  Protection.  Current transfer between soldered layers: a normal spot in one layer sheds its current into the layers against it over the")
    print(f"  transfer length sqrt(rho_c / (w r)), with rho_c the contact resistivity, w the width and r the tape's normal resistance per metre ({r70*1e3:.1f} mohm/m at {T_OP:.0f} K):")
    for lab, rc in (("soldered lap joint, 30 nohm cm2", R_CONTACT_SOLDER), ("fully solder-potted coil, 450 nohm cm2 (Fleiter and Ballarino, via Tokamak Energy)", R_CONTACT_SOLDER_HI)):
        lam = math.sqrt(rc / (0.012 * r70))
        print(f"    {lab:>88}: {lam*1e3:.2f} mm")
    print(f"  A defect in one layer of a ten-layer turn is bypassed within a millimetre; the nine others carry {10/9:.2f} times their share there, and the spot, carrying")
    print(f"  no current, makes no heat.  This is the defect tolerance Tokamak Energy reports for its solder-potted coils ('the electrical current can almost")
    print(f"  freely bypass a hotspot', MT27 slides; 'parallel tapes for mitigating defects', WAMHTS-5 slides; slide-level claims, not a paper) and the 'preferentially")
    print(f"  shared within a turn to bypass defects' of their patent US 11,101,060.  Their protection is a breaker that forces the current through a partial")
    print(f"  insulation between the turns (a distributed heater driven by the coil's own energy); here the turns are insulated and the equivalent is the dump.")
    print(f"  What is left for the dump is a quench of a whole turn or more (a loss of cooling, a mechanical event): the dump at 1 kV moves the coil's current")
    print(f"  into its nine neighbours in {L_leak*I0/1000*1e3:.0f} ms at 51 A and {L_leak*I_cap/1000*1e3:.0f} ms at 117 A, against the {bl.hot_spot_time(I0):.0f} and {bl.hot_spot_time(I_cap):.1f} s a bare spot has to 300 K (B.1 (d):")
    print(f"  the spot ends at 72 to 89 K after a 100 ms detection).  Heaters are not needed: the soldered layers protect against a local defect, the dump")
    print(f"  against a quench, and the heaters of A.5 remain only as a backup for a dump that fails to close.")
    # (6) losses
    print()
    P_strip = 0.01 * 1000.0 * I0; P_diode = V_DIODE_COLD * I0
    t_ch = rows[0]["t"][0]; t_trim = rows[3]["t"][0]
    print(f"  S1.7  Losses in the cold (551 turns, 1 kV): while charging, the diode {P_diode:.0f} W (D1; D2 only on the reverse half-cycles), the strip {P_strip:.0f} W (the 1 percent")
    print(f"  bypass), the secondary turns' AC loss (S1.2, tens of watts), the core warm: about {P_strip + P_diode:.0f} W for the {fmt_t(t_ch)} of a charge from zero, {(P_strip+P_diode)*t_ch/1e3:.0f} kJ; while")
    print(f"  holding, nothing (the strip superconducting, the diodes at zero volts); a 1 percent trim at 1 kV takes {t_trim:.1f} s and {(P_strip+P_diode)*t_trim:.0f} J.  The takeover of a parked")
    print(f"  spacecraft: {fmt_t(rows[1]['t'][0])} at up to {1000*I_cap/1e3:.0f} kW per coil from the bus, {(P_strip+P_diode)*rows[1]['t'][0]/1e3:.0f} kJ of cold heat per coil.  At 55 turns and 100 V the strip's share is")
    print(f"  the same and the diodes' ten times more, {V_DIODE_COLD*turns_rows[1]['I']:.0f} W.")
    # (7) what the slow coils can do, and the length
    print()
    print(f"  S1.8  What the slow coils can do: nothing inside an event (the strip must be heated open first, seconds, and the rate is then 1 percent per")
    print(f"  second at 1 kV, B.2 showed that is invisible in a 2 to 3 second event); between events, set, trim and take over a standing load in minutes")
    print(f"  (S1.1).  The length: with no pumps and no leads the parts per coil are a transformer, two diodes, a strip, a dump switch and a resistor,")
    print(f"  twenty sets per kilometre of casing at 1 km; the settling time at a given voltage is proportional to the length; the coil-end eddy heat is")
    print(f"  inversely proportional to it (A.6: 4.1 W/m at 1 km, 39 at 100 m, 386 at 10 m with step ends); the dump energy per coil is proportional to it.")
    print(f"  Nothing now asks for a coil longer than the book's 1 km, and a shorter one buys only faster charging at the price of plate heat and parts.  The")
    print(f"  book's 1 km staggered coils stay.")
    return out


# =============================================================================
# S2  LAYOUT (A): THE LOAD COMPENSATION COILS ON THEIR OWN TRACKS AND PLATES
# =============================================================================

PLATE_T = fb.PLATE_T; PLATE_W_SUP = fb.PLATE_W
W_LIFT_LOOP = bo.BEST_LIFT["W"]; PLATE_W_LIFT = W_LIFT_LOOP + 0.4           # 0.8 m loop, 1.2 m plate (the lift track's design, Part 3b)
BUNDLE = 0.10
SPACINGS = (0.1, 0.3, 0.5, 0.7, 1.0)                                       # m, plate edge to plate edge between a support track and the load track beside it


def two_track_linear(d_edge, h=0.01, NI_s=0.0, NI_l=0.0, plate_shift=0.0, gap_step=0.0, Lx=9.0, Ly=3.0):
    """The support track (two 0.5 m loops, 1.4 m plate) centred at x = 0 and a load track of the lift design (one 0.8 m loop, 1.2 m plate)
    beside it, plate edge to plate edge d_edge, both at the 100 mm gap, uniform 0.1 m bundles, linear plates (mu_r of the plate), solved
    once; returns the grid and the leg positions.  plate_shift moves the support plate sideways; gap_step raises the right half of it."""
    g = ms.Grid(Lx, Ly, h, x0=-3.0, y0=-1.0)
    mu = cf.PLATE_MU
    if gap_step:
        g.add_iron(-PLATE_W_SUP / 2 + plate_shift, plate_shift, bd.GAP0, bd.GAP0 + PLATE_T, mu)
        g.add_iron(plate_shift, PLATE_W_SUP / 2 + plate_shift, bd.GAP0 + gap_step, bd.GAP0 + gap_step + PLATE_T, mu)
    else:
        g.add_iron(-PLATE_W_SUP / 2 + plate_shift, PLATE_W_SUP / 2 + plate_shift, bd.GAP0, bd.GAP0 + PLATE_T, mu)
    xl = PLATE_W_SUP / 2 + d_edge + PLATE_W_LIFT / 2                       # the load plate's centre
    g.add_iron(xl - PLATE_W_LIFT / 2, xl + PLATE_W_LIFT / 2, bd.GAP0, bd.GAP0 + PLATE_T, mu)
    legs_s = [(-0.5, NI_s), (0.0, -2 * NI_s), (0.5, NI_s)]
    legs_l = [(xl - W_LIFT_LOOP / 2, NI_l), (xl + W_LIFT_LOOP / 2, -NI_l)]
    for x, a in legs_s + legs_l:
        if a:
            g.add_current(x - BUNDLE / 2, x + BUNDLE / 2, -BUNDLE, 0.0, a)
    g.solve()
    return g, legs_s, legs_l, xl


def mean_A(g, x):
    m = g.box(x - BUNDLE / 2, x + BUNDLE / 2, -BUNDLE, 0.0)
    return float(g.A[m].mean())


def fig8_flux(g, legs_s):
    """Flux per metre through the support track's two loops in their own senses, summed (the figure-eight circuit)."""
    A = [mean_A(g, x) for x, _ in legs_s]
    phi1 = A[0] - A[1]; phi2 = -(A[1] - A[2])                               # loop 1 (left, +), loop 2 (right, -) in own senses
    return phi1 + phi2, phi1, phi2


def part_layout(quick=False):
    hr("S2  LAYOUT (A): THE LOAD COMPENSATION COILS ON THEIR OWN TRACKS AND PLATES BESIDE THE SUPPORT TRACKS")
    r_sup, r_lift = be.load_records(); theva4 = bd.read_theva(); theva77 = bo.read_theva_file(bo.THEVA_77)
    rob = bo.read_robinson_angles(); theva = dict(theva4); theva["full"] = rob
    out = {}
    side = ring.CABLE_SIDE
    NI_l_full = r_lift["NI_req"] * math.sqrt(F_LOAD_ONE / r_lift["F_req"])
    print(f"  The load track is a copy of the lift track's design (Part 3b: one {W_LIFT_LOOP:.1f} m loop, two 200 x 50 mm bundles, a {PLATE_W_LIFT:.1f} m plate; {r_lift['NI_req']/1e3:.0f} kA per leg for")
    print(f"  {r_lift['F_req']/1e3:.1f} kN/m with the worst 3 mm tape at {r_lift['ic3']:.0f} A, {r_lift['tape3']:.1f} kg/m of tape at the 50 percent share): {NI_l_full/1e3:.0f} kA for its {F_LOAD_ONE/1e3:.0f} kN/m, {r_lift['tape3']*NI_l_full/r_lift['NI_req']:.1f} kg/m,")
    print(f"  against the {r_sup['tape3_446']:.0f} kg/m of the book's load winding in the support bundles (3 mm tape in the support's field at the 446 kN/m state).")
    # S2.1 the coupling against the spacing
    print()
    print(f"  S2.1  The coupling between a load track and the support track beside it (both under the cable, linear plates, {10:.0f} mm cells, uniform bundles):")
    print(f"  the flux the load loop puts through the support's figure-eight circuit per ampere-turn, against the support's own flux per ampere-turn, and")
    print(f"  what flux conservation then does to the support current when the load track pulls its {F_LOAD_ONE/1e3:.0f} kN/m ({NI_l_full/1e3:.0f} kA) beside a support at 141 kA.")
    g_s, legs_s, legs_l, xl = two_track_linear(0.5, NI_s=1e3)
    L_s, phi1_s, phi2_s = fig8_flux(g_s, legs_s); L_s /= 1e3
    print(f"  The support's figure-eight flux per ampere-turn {L_s*1e6:.3f} uWb/m, both loops (Part 1d's solver with the modelled winding: 2 x {bl.part1_constants(r_sup)['Lp_lin']*1e6:.3f} = {2*bl.part1_constants(r_sup)['Lp_lin']*1e6:.3f}).")
    print(f"  {'edge spacing m':>14} {'load plate centre m':>19} {'M uWb/m per At':>14} {'loop 1 / loop 2':>24} {'M / L_s':>9} {'dNI_s per kA of load':>20} {'force change at 198 kN/m':>24} {'fits 7.39 m':>11}")
    rows = []
    for d in SPACINGS:
        g, legs_s, legs_l, xl = two_track_linear(d, NI_l=1e3)
        M, p1, p2 = fig8_flux(g, legs_s); M /= 1e3; p1 /= 1e3; p2 /= 1e3
        dNI = -M / L_s                                                      # support ampere-turns per ampere-turn of load (flux conserved)
        dF = 2 * dNI * NI_l_full / 141e3                                    # fraction of the support force
        x_outer = xl + PLATE_W_LIFT / 2                                     # the load plate's outer edge from the support's centre
        fits = (Y_PLATE_CENTER + x_outer - PLATE_W_SUP / 2 * 0) <= side / 2  # see the fit below
        rows.append(dict(d=d, xl=xl, M=M, p1=p1, p2=p2, k=M / L_s, dNI_per_kA=dNI * 1e3, dF=dF))
        print(f"  {d:14.1f} {xl:19.2f} {M*1e6:14.4f} {p1*1e6:11.4f} / {p2*1e6:9.4f} {M/L_s:9.5f} {dNI*1e3:17.2f} At {dF*100:22.2f} %")
    out["coupling"] = rows
    print(f"  (M / L_s is the share of the support's flux the load loop can move; the figure-eight links only the difference between the two loops, so the")
    print(f"  coupling falls with the distance faster than a single loop's would.  The lift tracks on the ceiling couple to the support through the cable's")
    print(f"  two plates: Part 1d D.9.)")
    # S2.2 the plates and the fit
    print()
    y_sup = side / 2 - PLATE_W_SUP / 2
    print(f"  S2.2  The plates on the {side:.2f} m underside: two support plates of {PLATE_W_SUP:.1f} m at the edges (centres at {y_sup:.2f} m) and two load plates of {PLATE_W_LIFT:.1f} m inside them:")
    for d in SPACINGS:
        y_l = y_sup - PLATE_W_SUP / 2 - d - PLATE_W_LIFT / 2
        inner = y_l - PLATE_W_LIFT / 2
        print(f"    spacing {d:.1f} m: load plate centres at {y_l:.2f} m, inner edges {2*inner:.2f} m apart{' (overlap: does not fit)' if inner < 0 else ''}; {2*PLATE_W_SUP + 2*PLATE_W_LIFT + 2*d:.2f} m of plate and gaps on {side:.2f} m")
    m_plate_sup = PLATE_T * PLATE_W_SUP * bd.RHO_PLATE; m_plate_lift = PLATE_T * PLATE_W_LIFT * bd.RHO_PLATE
    print(f"  Plate masses at 200 mm: {m_plate_sup:,.0f} kg/m for a {PLATE_W_SUP:.1f} m plate, {m_plate_lift:,.0f} kg/m for a {PLATE_W_LIFT:.1f} m one; six plates (two support, two load, two lift) {2*m_plate_sup + 4*m_plate_lift:,.0f} kg/m against")
    print(f"  the four of {4*m_plate_sup:,.0f}; the cable hardware {fb.M_LIM_PLATES + 2*m_plate_sup + 4*m_plate_lift:,.0f} kg/m against 10,686 (the lift plates of the book are 1.4 m, {m_plate_sup:,.0f} kg/m; at the lift")
    print(f"  design's 1.2 m they are {m_plate_lift:,.0f}).  S6 gives the lighter plates; the design point is rerun there.")
    out["plates"] = dict(m_sup=m_plate_sup, m_lift=m_plate_lift, six=2 * m_plate_sup + 4 * m_plate_lift, hw_six=fb.M_LIM_PLATES + 2 * m_plate_sup + 4 * m_plate_lift)
    # S2.3 what it removes
    print()
    print(f"  S2.3  What layout (A) removes: the pulsed field in the support bundles (Part 1's 12 mm persistent tapes lost 4.5 times the 3 mm tapes' strip loss")
    print(f"  per pulse; now none), the combined-field limit (Table 7.6's rule that support plus load may not exceed 493,000 ampere-turns, and the 3 mm")
    print(f"  tape's 2.38 T in the design state: the load tape now sees only its own track's {r_lift['B_mean_req']:.2f} T mean, worst tape {r_lift['ic3']:.0f} A), the dependence on symmetry")
    print(f"  (the coupling is a geometric mutual, S2.1, not a cancellation), and the shared cryostat (each track in its own).  What it costs: two more plates")
    print(f"  on the cable, {4*m_plate_lift:,.0f} kg/m at 200 mm and 1.2 m (S6 for less), and the roll moment of a load track acting alone is that of its")
    print(f"  own offset from the centreline, the same as a lift track's.")
    # S2.4 layout (B) far enough to compare
    print()
    print(f"  S2.4  Layout (B), the single wide load loop on the outer bundles of the support track (Part 1 D.8(a)), with the saturating plate at 5 mm cells:")
    h = 0.005
    tr = be.Track(bo.BEST_SUPPORT, h)
    NI_s0 = bo.NI_of_F(r_sup, be.W_HALF); NI_cap = r_sup["NI_req"]
    def outer(NIl):
        v = [0.0] * len(tr.legs); v[0] += NIl; v[-1] -= NIl; return v
    # the single loop's own pull curve: two solves, log-log
    NIa = 468e3
    g_a, s_a = tr.solve(outer(NIa))
    NIb = NIa * math.sqrt(F_LOAD_ONE / s_a["F"])
    g_b, s_b = tr.solve(outer(NIb))
    NI_lB = math.exp(math.log(NIa) + (math.log(F_LOAD_ONE) - math.log(s_a["F"])) * (math.log(NIb) - math.log(NIa)) / (math.log(s_b["F"]) - math.log(s_a["F"])))
    g_l, s_l = tr.solve(outer(NI_lB))
    print(f"  The single loop alone: {s_a['F']/1e3:.1f} kN/m at {NIa/1e3:.0f} kA, {s_b['F']/1e3:.1f} at {NIb/1e3:.0f} kA; {F_LOAD_ONE/1e3:.0f} kN/m at {NI_lB/1e3:.0f} kA ({s_l['F']/1e3:.1f} solved), plate {s_l['sat']*100:.1f} % saturated.")
    wt3, lim3 = be.worst_tape(g_l, tr, outer(NI_lB), 0.003, theva, theva77)
    turns3 = 2 * NI_lB / (0.5 * wt3[0]); tape3 = turns3 * 0.003 * bo.TAPE_T * bo.TAPE_RHO
    print(f"  Its worst 3 mm tape alone: {wt3[1]:.2f} T at {wt3[2]:.0f} degrees, {wt3[0]:.1f} A; at the 50 percent share {turns3:,.0f} turns in its two legs, {tape3:.1f} kg/m of tape (the figure-eight winding of the book: {r_sup['tape3_446']:.0f};")
    print(f"  the lift-type track of layout (A): {r_lift['tape3']*NI_l_full/r_lift['NI_req']:.1f}).")
    rowsB = []
    for lab, NIs in (("support at 141 kA", NI_s0), ("support at its 322 kA cap", NI_cap)):
        legs_c = [a + b for a, b in zip(tr.NI_design(NIs), outer(NI_lB))]
        g_c, s_c = tr.solve(legs_c); g_s, s_s = tr.solve(tr.NI_design(NIs))
        wt12, _ = be.worst_tape(g_c, tr, legs_c, 0.012, theva, theva77); wt3c, _ = be.worst_tape(g_c, tr, legs_c, 0.003, theva, theva77)
        wt12s, _ = be.worst_tape(g_s, tr, tr.NI_design(NIs), 0.012, theva, theva77)
        # the force distribution along the plate: Maxwell stress on a line just under the plate; the centre of pull and the moment about the track centre
        y_line = bd.GAP0 - 2 * h
        xs, By, Bx = g_c.line_By(y_line, -PLATE_W_SUP / 2 - 0.1, PLATE_W_SUP / 2 + 0.1)
        f_x = (By**2 - Bx**2) / (2 * MU0)
        F_line = float(np.trapezoid(f_x, xs)); x_cp = float(np.trapezoid(xs * f_x, xs)) / F_line
        I_s = NIs / (N_SUP * bl.part1_constants(r_sup)["n_t"])
        rowsB.append(dict(label=lab, F=s_c["F"], F_sup=s_s["F"], F_sum=s_s["F"] + s_l["F"], sat=s_c["sat"], B12=wt12[1], th12=wt12[2], ic12=wt12[0], B12_alone=wt12s[1], ic12_alone=wt12s[0],
                          B3=wt3c[1], th3=wt3c[2], ic3=wt3c[0], x_cp=x_cp, M=F_line * x_cp))
        print(f"  {lab}: both windings {s_c['F']/1e3:.1f} kN/m (support alone {s_s['F']/1e3:.1f} + load alone {s_l['F']/1e3:.1f} = {(s_s['F']+s_l['F'])/1e3:.1f}: the cross term {(s_c['F']-s_s['F']-s_l['F'])/1e3:+.1f}), plate {s_c['sat']*100:.0f} % saturated;")
        print(f"    worst 12 mm tape {wt12[1]:.2f} T at {wt12[2]:.0f} deg, limit {wt12[0]:.0f} A (alone {wt12s[1]:.2f} T, {wt12s[0]:.0f} A) against its {I_s:.0f} A ({I_s/wt12[0]*100:.0f} % of the limit); worst 3 mm tape {wt3c[1]:.2f} T, {wt3c[0]:.0f} A;")
        print(f"    the pull's centre {x_cp*1e3:+.0f} mm from the track's centre: a moment of {F_line*x_cp/1e3:.1f} kN m per metre about it, and {F_line*(Y_PLATE_CENTER + x_cp)/1e3:.0f} kN m per metre about the casing's")
        print(f"    centreline when this track's load winding acts alone (the lift-type track of (A) acting alone: {F_LOAD_ONE*(side/2 - PLATE_W_SUP - 0.5 - PLATE_W_LIFT/2)/1e3:.0f} kN m per metre at a 0.5 m spacing).")
    out["layoutB"] = dict(NI_l=NI_lB, tape3=tape3, rows=rowsB)
    # the residual coupling with an imperfect symmetry (linear, the two-track grid with the load loop on the support's outer legs)
    print()
    print(f"  The residual coupling of (B) with an imperfect symmetry (linear plate, 10 mm cells, the load loop on the outer legs at 1 kA, the figure-eight flux")
    print(f"  it links, against the support's own {L_s*1e6:.3f} uWb/m per ampere-turn):")
    def fig8_from_outer(plate_shift=0.0, gap_step=0.0):
        g = ms.Grid(6.0, 3.0, 0.01, x0=-3.0, y0=-1.0); mu = cf.PLATE_MU
        if gap_step:
            g.add_iron(-PLATE_W_SUP / 2 + plate_shift, plate_shift, bd.GAP0, bd.GAP0 + PLATE_T, mu)
            g.add_iron(plate_shift, PLATE_W_SUP / 2 + plate_shift, bd.GAP0 + gap_step, bd.GAP0 + gap_step + PLATE_T, mu)
        else:
            g.add_iron(-PLATE_W_SUP / 2 + plate_shift, PLATE_W_SUP / 2 + plate_shift, bd.GAP0, bd.GAP0 + PLATE_T, mu)
        for x, a in ((-0.5, 1e3), (0.5, -1e3)):
            g.add_current(x - BUNDLE / 2, x + BUNDLE / 2, -BUNDLE, 0.0, a)
        g.solve()
        M, p1, p2 = fig8_flux(g, [(-0.5, 0), (0.0, 0), (0.5, 0)])
        return M / 1e3
    resid = []
    for lab, kw, scale in (("perfect symmetry", {}, 1.0), ("the casing 10 mm off-centre (the plate shifted 10 mm)", dict(plate_shift=0.01), 1.0),
                           ("one loop's gap 5 mm larger (the plate stepped at the centre)", dict(gap_step=0.005), 1.0),
                           ("the casing rolled 1 milliradian (the step result scaled to 0.6 mm across a loop)", dict(gap_step=0.005), 0.6 / 5.0)):
        M = fig8_from_outer(**kw) * scale
        resid.append(dict(label=lab, M=M, k=M / L_s, dF=2 * (-M / L_s) * NI_lB / 141e3))
        print(f"    {lab:>78}: M {M*1e9:+8.3f} nWb/m per At, M / L_s {M/L_s:+.5f}, support force {2*(-M/L_s)*NI_lB/141e3*100:+.3f} % at the parked load")
    print(f"    {'one of the ten support coils dead':>78}: no asymmetry (the figure-eight carries both loops in series, so a dead coil takes a tenth from both)")
    out["residual"] = resid
    return out


# =============================================================================
# P2  THE SUPPORT TRACK AT 75 AND 50 MM
# =============================================================================

def part_gap(quick=False):
    hr("P2  THE SAME BEARING AT A 75 AND A 50 MM GAP (Part 2 of the electrical brief)")
    r_sup, r_lift = be.load_records(); theva4 = bd.read_theva(); theva77 = bo.read_theva_file(bo.THEVA_77)
    rob = bo.read_robinson_angles(); theva = dict(theva4); theva["full"] = rob
    h = 0.005; out = {}
    F_s = r_sup["F_req"]; F_half = be.W_HALF
    print(f"  Support track (Part 3b geometry, {h*1e3:g} mm cells), the plate saturating; at each gap the ampere-turns for {F_s/1e3:.1f} kN/m and {F_half/1e3:.1f} kN/m by two nonlinear solves")
    print(f"  and a log-log interpolation, the worst 12 mm tape there, the turns at the 50 percent share, the tape mass (Equation 7.3), the plate saturation.")
    print(f"  {'gap mm':>6} {'NI for 245.5 kA':>15} {'NI for 47.5 kA':>14} {'worst tape T, deg, A':>22} {'turns per leg':>13} {'tape kg/m':>9} {'sat %':>5} {'Lp uWb/m per At':>15} {'k_I kN/m per mm':>15} {'k_phi':>7}")
    gaps = (0.05, 0.075, 0.10) if quick else (0.05, 0.06, 0.075, 0.09, 0.10)
    sweep = be.gap_sweep(bo.BEST_SUPPORT, h, gaps, NI_lin=1e3)
    rows = []
    for gp in (0.05, 0.075, 0.10):
        tr = be.Track(bo.BEST_SUPPORT, h, gap=gp)
        NI1 = 322e3 * (gp / 0.10)**0.9
        g1, s1 = tr.solve(tr.NI_design(NI1)); NI2 = NI1 * math.sqrt(F_s / s1["F"]); g2, s2 = tr.solve(tr.NI_design(NI2))
        NI_req = math.exp(math.log(NI1) + (math.log(F_s) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI1)) / (math.log(s2["F"]) - math.log(s1["F"])))
        g3, s3 = tr.solve(tr.NI_design(NI_req))
        wt, _ = be.worst_tape(g3, tr, tr.NI_design(NI_req), 0.012, theva, theva77)
        I50 = 0.5 * wt[0]; turns = NI_req / I50; tape = 4 * turns * 0.012 * bo.TAPE_T * bo.TAPE_RHO
        NI_half = NI_req * math.sqrt(F_half / F_s)                          # near-linear below the cap
        k_I, k_phi, _, _ = be.stiffness(sweep, gp, NI_req)
        Lp = float(np.interp(gp, [r["gap"] for r in sweep], [r["Lp"] for r in sweep]))
        rows.append(dict(gap=gp, NI_req=NI_req, NI_half=NI_half, B=wt[1], th=wt[2], ic=wt[0], turns=turns, tape=tape, sat=s3["sat"], Lp=Lp, k_I=k_I, k_phi=k_phi, F=s3["F"]))
        print(f"  {gp*1e3:6.0f} {NI_req/1e3:12.0f} kA {NI_half/1e3:11.0f} kA {wt[1]:7.2f} {wt[2]:5.0f} {wt[0]:6.0f} {turns:13,.0f} {tape:9.1f} {s3['sat']*100:5.1f} {Lp*1e6:15.3f} {k_I/1e6:15.2f} {k_phi/1e6:7.2f}")
    out["support"] = rows
    r100 = rows[-1]; r50 = rows[0]
    print(f"  At 50 mm the track needs {r50['NI_req']/r100['NI_req']*100:.0f} percent of the ampere-turns of the 100 mm design for the same pull (the hand estimate 73), the worst tape sits in")
    print(f"  {r50['B']:.2f} T against {r100['B']:.2f}, and the tape is {r50['tape']/r100['tape']*100:.0f} percent of the 100 mm mass; the negative stiffness at constant current is {r50['k_I']/r100['k_I']:.1f} times larger.")
    # the lift track at 196.5 kN/m
    print()
    print(f"  Lift track (and the load track of layout (A)), {F_LIFT_ONE/1e3:.1f} kN/m, 3 mm tape:")
    print(f"  {'gap mm':>6} {'NI kA':>8} {'worst tape T, deg, A':>22} {'turns per leg':>13} {'tape kg/m':>9} {'sat %':>5} {'E per 50 m coil kJ':>18}")
    rowsL = []
    for gp in (0.05, 0.075, 0.10):
        tr = be.Track(bo.BEST_LIFT, h, gap=gp, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=421e3)
        NI1 = 421e3 * (gp / 0.10)**0.9
        g1, s1 = tr.solve(tr.NI_design(NI1)); NI2 = NI1 * math.sqrt(F_LIFT_ONE / s1["F"]); g2, s2 = tr.solve(tr.NI_design(NI2))
        NI_req = math.exp(math.log(NI1) + (math.log(F_LIFT_ONE) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI1)) / (math.log(s2["F"]) - math.log(s1["F"])))
        g3, s3 = tr.solve(tr.NI_design(NI_req))
        wt, _ = be.worst_tape(g3, tr, tr.NI_design(NI_req), 0.003, theva, theva77)
        I50 = 0.5 * wt[0]; turns = NI_req / I50; tape = 2 * turns * 0.003 * bo.TAPE_T * bo.TAPE_RHO
        rowsL.append(dict(gap=gp, NI_req=NI_req, B=wt[1], th=wt[2], ic=wt[0], turns=turns, tape=tape, sat=s3["sat"], W=s3["W"]))
        print(f"  {gp*1e3:6.0f} {NI_req/1e3:8.0f} {wt[1]:7.2f} {wt[2]:5.0f} {wt[0]:6.0f} {turns:13,.0f} {tape:9.1f} {s3['sat']*100:5.1f} {s3['W']*50/fb.N_FAM/1e3:18.1f}")
    out["lift"] = rowsL
    # a re-optimization at 50 mm: narrower loops
    if not quick:
        print()
        print(f"  One re-optimization at 50 mm: the support's loop width free (the bundles as they are), {F_s/1e3:.1f} kN/m:")
        print(f"  {'loop W m':>8} {'NI kA':>8} {'worst tape T, A':>16} {'turns per leg':>13} {'tape kg/m':>9} {'sat %':>5} {'track width m':>13}")
        for W in (0.3, 0.4, 0.5, 0.6):
            geom = dict(bo.BEST_SUPPORT, W=W)
            tr = be.Track(geom, h, gap=0.05)
            NI1 = r50["NI_req"]; g1, s1 = tr.solve(tr.NI_design(NI1)); NI2 = NI1 * math.sqrt(F_s / s1["F"]); g2, s2 = tr.solve(tr.NI_design(NI2))
            NI_req = math.exp(math.log(NI1) + (math.log(F_s) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI1)) / (math.log(s2["F"]) - math.log(s1["F"])))
            g3, s3 = tr.solve(tr.NI_design(NI_req))
            wt, _ = be.worst_tape(g3, tr, tr.NI_design(NI_req), 0.012, theva, theva77)
            turns = NI_req / (0.5 * wt[0]); tape = 4 * turns * 0.012 * bo.TAPE_T * bo.TAPE_RHO
            print(f"  {W:8.1f} {NI_req/1e3:8.0f} {wt[1]:9.2f} {wt[0]:6.0f} {turns:13,.0f} {tape:9.1f} {s3['sat']*100:5.1f} {2*W + tr.wb:13.2f}")
    # what the smaller gap costs in the cases: the excursions that exceed the gap
    print()
    print(f"  What the smaller gap costs in the dynamics: final_bearing's lanes have gap_dep = False (the driven windings' current is set for the measured")
    print(f"  gap), so the casing's movement in a case does not depend on the nominal gap and a case touches when its excursion exceeds the gap.  S4 counts")
    print(f"  the single-failure cases of the chosen combination whose excursion exceeds 50 and 25 mm, which would touch at a 50 and a 75 mm gap.")
    return out


def parse_case_gaps(path):
    """(label, smallest gap) of every case in a Part B log's F.4-style listing is not printed; read the per-case tables instead."""
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.rstrip()
            parts = s.split()
            # F.2 / B.5 rows: label ... exc gap_b gap_t [engaged|contact]
            if len(parts) >= 4 and s.startswith("  ") and not s.startswith("  F.") and not s.startswith("  ("):
                try:
                    nums = [float(x) for x in parts[-4:-1]]
                except ValueError:
                    try:
                        nums = [float(x) for x in parts[-3:]]
                    except ValueError:
                        continue
                if len(nums) == 3 and all(0 <= v <= 100.0001 for v in nums):
                    lab = " ".join(parts[:-4]) if len(parts) >= 5 else " ".join(parts[:-3])
                    if not any(ch.isalpha() for ch in lab) or lab.startswith("u km/s"):
                        continue                                             # the F.1 speed rows and headers are not cases
                    out.append((lab, min(nums[1], nums[2]) * 1e-3))
    return out


# =============================================================================
# P3  THE FIELDS FOR TABLE 7.4
# =============================================================================

def part_table74(quick=False):
    hr("P3  THE FIELDS IN TESLA FOR TABLE 7.4 (Part 3 of the electrical brief)")
    r_sup, r_lift = be.load_records(); theva4 = bd.read_theva(); theva77 = bo.read_theva_file(bo.THEVA_77)
    rob = bo.read_robinson_angles(); theva = dict(theva4); theva["full"] = rob
    h = 0.005 if quick else 0.0025
    tr = be.Track(bo.BEST_SUPPORT, h)
    print(f"  Support track at {h*1e3:g} mm cells, the plate saturating; both windings in the same bundles (the book's layout), so the rows are the total ampere-turns per leg.")
    print(f"  {'kA per leg':>10} {'pull kN/m':>9} {'book':>5} {'worst 12 mm tape T, deg, A':>27} {'worst 3 mm tape T, deg, A':>26} {'plate face: peak, mean, rms T':>29} {'inside the plate peak T':>23} {'book column T':>13} {'sat %':>5}")
    rows = []
    book = {141e3: (48, 0.29), 287e3: (198, 0.60), 322e3: (245, 0.66), 403e3: (330, 0.77), 493e3: (444, 0.89), 610e3: (601, 1.04)}
    for NI in (141e3, 287e3, 322e3, 403e3, 493e3, 610e3):
        if quick and NI not in (141e3, 322e3, 493e3):
            continue
        g, s_ = tr.solve(tr.NI_design(NI))
        wt12, _ = be.worst_tape(g, tr, tr.NI_design(NI), 0.012, theva, theva77); wt3, _ = be.worst_tape(g, tr, tr.NI_design(NI), 0.003, theva, theva77)
        xs, By, Bx = g.line_By(bd.GAP0 - 2 * h, -PLATE_W_SUP / 2, PLATE_W_SUP / 2)
        m = g.box(-PLATE_W_SUP / 2, PLATE_W_SUP / 2, bd.GAP0, bd.GAP0 + PLATE_T)
        B_in = float(np.sqrt(g.Bx[m]**2 + g.By[m]**2).max())
        B_rms = float(np.sqrt(np.mean(By**2))); B_mean = float(np.mean(np.abs(By))); B_pk = float(np.abs(By).max())
        B_even = math.sqrt(2 * MU0 * s_["F"] / PLATE_W_SUP)
        rows.append(dict(NI=NI, F=s_["F"], B12=wt12[1], th12=wt12[2], ic12=wt12[0], B3=wt3[1], th3=wt3[2], ic3=wt3[0], B_pk=B_pk, B_mean=B_mean, B_rms=B_rms, B_in=B_in, B_even=B_even, sat=s_["sat"]))
        print(f"  {NI/1e3:10.0f} {s_['F']/1e3:9.1f} {book[NI][0]:5.0f} {wt12[1]:9.2f} {wt12[2]:5.0f} {wt12[0]:6.0f} {wt3[1]:9.2f} {wt3[2]:5.0f} {wt3[0]:6.0f} {B_pk:9.2f} {B_mean:6.2f} {B_rms:6.2f} {B_in:23.2f} {B_even:6.2f} ({book[NI][1]:.2f}) {s_['sat']*100:5.0f}")
    print(f"  (the book column is sqrt(2 mu0 F / w) with w = 1.4 m from the solved pull, the book's own value in brackets; the rms normal field on the face is what that")
    print(f"  formula stands for.  The book says 'the strongest field any tape sees is 1.64 T' at the cap: that is the 12 mm tape at 322 kA; the 3 mm tape in the")
    print(f"  design state (493 kA) sits in the field of the row above it.)")
    rowsL = []
    trl = be.Track(bo.BEST_LIFT, h, w_tape=0.003, I_ref=bo.I_REF_3, NI_ref=421e3)
    print(f"  Lift track:")
    for NI, lab in ((r_lift["NI_req"], "196.5 kN/m"), (r_lift["NI_rated"], "the 1 T rating")):
        g, s_ = trl.solve(trl.NI_design(NI))
        wt3, _ = be.worst_tape(g, trl, trl.NI_design(NI), 0.003, theva, theva77)
        xs, By, Bx = g.line_By(bd.GAP0 - 2 * h, -PLATE_W_LIFT / 2, PLATE_W_LIFT / 2)
        m = g.box(-PLATE_W_LIFT / 2, PLATE_W_LIFT / 2, bd.GAP0, bd.GAP0 + PLATE_T)
        B_in = float(np.sqrt(g.Bx[m]**2 + g.By[m]**2).max())
        rowsL.append(dict(NI=NI, F=s_["F"], B3=wt3[1], th3=wt3[2], ic3=wt3[0], B_pk=float(np.abs(By).max()), B_mean=float(np.mean(np.abs(By))), B_rms=float(np.sqrt(np.mean(By**2))), B_in=B_in, sat=s_["sat"]))
        print(f"    {lab:>14}: {NI/1e3:.0f} kA, {s_['F']/1e3:.1f} kN/m; worst 3 mm tape {wt3[1]:.2f} T at {wt3[2]:.0f} deg, {wt3[0]:.0f} A; plate face peak {rowsL[-1]['B_pk']:.2f}, mean {rowsL[-1]['B_mean']:.2f}, rms {rowsL[-1]['B_rms']:.2f} T; inside {B_in:.2f} T; {s_['sat']*100:.0f} % saturated")
    return dict(support=rows, lift=rowsL, h=h)


# =============================================================================
# P4  THE PLATE HEAT OF A LAUNCH ON THE PRESENT TRACKS
# =============================================================================

def plate_heat_pass(Bp, u, L_coil, plate_t, plate_w):
    """Hysteresis and eddy heat per metre of one plate for one pass at speed u: final_bearing C.2's model with the plate's dimensions free."""
    L_on = u * max(0.05, bd.L_SC / u) + bd.L_SC + 2 * L_coil + 2 * u * bd.TAU_HW
    v_sw = u + bd.V_REL
    hyst = lc.FERRO_HYSTERESIS_COEFF * (Bp / 1.2) ** 1.6 * plate_t * plate_w
    eddy = bd.eddy_cycle(Bp, bd.GAP0 / v_sw, L_on / v_sw, d_p=plate_t, mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY) * plate_w
    return hyst, eddy, L_on / v_sw


def part_heat(quick=False, plate_t_lift=PLATE_T, plate_w_lift=PLATE_W_SUP, plate_t_load=PLATE_T, plate_w_load=PLATE_W_SUP, L_lift=None, L_load=None, label="the book's plates (1.4 m, 200 mm)"):
    hr(f"P4  THE PLATE HEAT OF A LAUNCH ON THE PRESENT TRACKS, {label} (Part 4 of the electrical brief)")
    theva = bd.read_theva(); sup, lift = fb.load_air()
    L_lift = L_lift or fb.L_LIFT; L_load = L_load or fb.L_LOAD
    t_p, v_p, P_p = pz.lsm_profile(); dt_p = np.diff(t_p)
    KE = 0.5 * (bd.M_SC0 + bd.M_SLED_M * bd.L_SLED) * 30000**2
    print(f"  Pulls per track from Part 3b's launch (bearing_optimize.PASS_PULLS: the pair's pull at each speed), the plate's peak field from the air-core curve of")
    print(f"  the track at that pull (the lift track's curve for the load track too in layout (A)); hysteresis 50 J/m3 per cycle at 1.2 T scaled as (B/1.2)^1.6")
    print(f"  [UNVERIFIED], eddy by 1-D diffusion; hold = the energized length over u + V_rel.  Per plate (one of the pair) and per metre of ring (two plates).")
    print(f"  {'u km/s':>7} {'winding':>7} {'pull per track kN/m':>19} {'plate peak T':>12} {'hold ms':>7} {'hyst J/m':>8} {'eddy J/m':>8} {'per plate J/m':>13} {'per m of ring J/m':>17}")
    per_u = {}
    rows = []
    for u_km, F_pair, kind in bo.PASS_PULLS:
        u = u_km * 1e3; Fu = F_pair * 1e3 / 2
        rec = lift                                                          # layout (A): the load track is the lift design; the book's load winding sits on the support loops
        NI, Bp = fb.curve_at(rec, Fu, "B_peak")
        pt, pw, L = (plate_t_lift, plate_w_lift, L_lift) if kind == "lift" else (plate_t_load, plate_w_load, L_load)
        hyst, eddy, hold = plate_heat_pass(Bp, u, L, pt, pw)
        per_u[u] = 2 * (hyst + eddy)
        rows.append(dict(u=u, kind=kind, F=Fu, Bp=Bp, hyst=hyst, eddy=eddy, hold=hold))
        print(f"  {u_km:7.0f} {kind:>7} {Fu/1e3:19.0f} {Bp:12.2f} {hold*1e3:7.1f} {hyst:8.1f} {eddy:8.1f} {hyst+eddy:13.1f} {2*(hyst+eddy):17.1f}")
    us = sorted(per_u)
    E_pl = 0.0
    for i, dtt in enumerate(dt_p):
        vv = v_p[i + 1]
        if vv < us[0]:
            continue
        E_pl += float(np.interp(vv, us, [per_u[u] for u in us])) * (vv + bd.V_REL) * dtt
    print(f"  Per launch {E_pl/1e9:.1f} GJ in the plates ({E_pl/KE*100:.4f} % of the launch energy; the book: 238 GJ and 1,090 J per metre of plate at 30 km/s).")
    return dict(rows=rows, E_launch=E_pl)


# =============================================================================
# S6  THE PLATES: THE LIGHTEST PLATE FOR EACH KIND OF TRACK
# =============================================================================

PLATE_OVERS = (0.0, 0.1, 0.2)
PLATE_TS = (0.05, 0.10, 0.15, 0.20)
AIR_LEVELS_B = (100e3, 200e3, 300e3, 400e3, 500e3)


def plate_case(kind, geom, theva, F_req, w_tape, levels=AIR_LEVELS_B):
    """One plate geometry of the air-core model (uniform 0.1 m bundles, bearing_mass.build_air): the pull curve, the ampere-turns for the
    requirement, the plate's saturation and peak field there, the winding's mean field and the tape at the 50 percent share of the tape's
    limit at that mean field (bearing_mass.air_tape's rule), the plate mass."""
    t0 = time.time()
    curve, meta = bm.air_curve(levels=levels, **geom)
    rows = [c for c in curve if c.get("nonlinear")]
    Fs = [c["F"] for c in rows]; NIs = [c["NI"] for c in rows]
    rec = dict(curve=curve, meta={k: v for k, v in meta.items() if k != "legs"})
    out = dict(kind=kind, geom=geom, plate_w=meta["plate_w"], plate_m=meta["plate_m2"] * bd.RHO_PLATE, F_max=max(Fs), seconds=time.time() - t0)
    if F_req <= max(Fs):
        tp = bm.air_tape(rec, F_req, theva, w_tape)
        NI = tp["NI"]
        for k in ("B_peak", "sat_frac", "B_plate_max", "B_plate_mid", "W_mag"):
            out[k] = float(np.interp(NI, NIs, [c[k] for c in rows]))
        out.update(NI=NI, Bw=tp["Bw"], Ic=tp["Ic"], I_op=tp["I_op"], tape=tp["mass"], reached=True)
    else:
        out.update(reached=False, NI=float("nan"), tape=float("nan"), sat_frac=float("nan"), B_peak=float("nan"), B_plate_max=float("nan"), Bw=float("nan"), Ic=float("nan"))
    return out


def part_plates(quick=False):
    hr("S6  THE PLATES: THE LIGHTEST PLATE FOR EACH KIND OF TRACK (air-core model, 10 mm cells, uniform 0.1 m bundles)")
    theva = bd.read_theva()
    r_sup, r_lift = be.load_records()
    F_s = r_sup["F_req"]
    print(f"  The model is bearing_mass's (the one final_bearing sizes its masses with): each leg a 0.1 m square bundle of uniform current, the plate saturating at")
    print(f"  {cf.PLATE_B_SAT} T, the tape at half of its limit at the mean field over the bundles (the Part 3b design uses the worst tape at 2.5 mm cells and gives more")
    print(f"  tape: 94 against 67 kg/m for the support; the ratios between plates are what this table is for).  Support: two 0.5 m loops, {F_s/1e3:.1f} kN/m, 12 mm tape.")
    print(f"  Lift and load (layout (A)): one 0.8 m loop, {F_LIFT_ONE/1e3:.1f} kN/m, 3 mm tape.  Overhang = plate beyond the loop's outer legs each side; 'shaped' = thin over")
    print(f"  the legs and thick between them.")
    kinds = (("support", dict(n_loops=2, W=0.5), F_s, 0.012), ("lift / load", dict(n_loops=1, W=0.8), F_LIFT_ONE, 0.003))
    out = {}
    for kind, base, F_req, w_tape in kinds:
        print()
        print(f"  {kind}, {F_req/1e3:.1f} kN/m per track:")
        print(f"  {'overhang m':>10} {'thickness mm':>12} {'plate w m':>9} {'plate kg/m':>10} {'curve max kN/m':>14} {'NI kA per leg':>13} {'sat %':>5} {'plate peak T':>12} {'inside T':>8} {'winding T':>9} {'Ic A':>5} {'tape kg/m':>9} {'plate + 8.6 x tape kg/m':>23} {'s':>4}")
        cases = []
        for po in PLATE_OVERS:
            for pt in PLATE_TS:
                if quick and (po, pt) not in ((0.2, 0.2), (0.0, 0.1), (0.1, 0.1)):
                    continue
                cases.append(dict(base, plate_t=pt, plate_over=po))
        if not quick:
            cases.append(dict(base, plate_t=0.10, plate_over=0.1, plate_shape=0.20))
            cases.append(dict(base, plate_t=0.05, plate_over=0.1, plate_shape=0.15))
        rows = []
        for geom in cases:
            r = plate_case(kind, geom, theva, F_req, w_tape)
            rows.append(r)
            po = geom["plate_over"]; pt = geom["plate_t"]; sh = geom.get("plate_shape")
            lab_t = f"{pt*1e3:.0f}" + (f" ({sh*1e3:.0f} between)" if sh else "")
            if r["reached"]:
                print(f"  {po:10.1f} {lab_t:>12} {r['plate_w']:9.2f} {r['plate_m']:10.0f} {r['F_max']/1e3:14.0f} {r['NI']/1e3:13.0f} {r['sat_frac']*100:5.0f} {r['B_peak']:12.2f} {r['B_plate_max']:8.2f} {r['Bw']:9.2f} {r['Ic']:5.0f} {r['tape']:9.1f} {r['plate_m'] + 8.6*r['tape']:23.0f} {r['seconds']:4.0f}", flush=True)
            else:
                print(f"  {po:10.1f} {lab_t:>12} {r['plate_w']:9.2f} {r['plate_m']:10.0f} {r['F_max']/1e3:14.0f} {'not reached':>13} {r['seconds']:>70.0f}", flush=True)
        out[kind] = rows
        ok = [r for r in rows if r["reached"]]
        if ok:
            light = min(ok, key=lambda r: r["plate_m"]); best = min(ok, key=lambda r: r["plate_m"] + 8.6 * r["tape"])
            print(f"  Lightest plate that reaches the pull: overhang {light['geom']['plate_over']:.1f} m, {light['geom']['plate_t']*1e3:.0f} mm, {light['plate_m']:.0f} kg/m ({light['sat_frac']*100:.0f} % saturated, {light['tape']:.1f} kg/m of tape);")
            print(f"  lightest for the cable (plate plus 8.6 times the tape, the cable's cost of casing mass): overhang {best['geom']['plate_over']:.1f} m, {best['geom']['plate_t']*1e3:.0f} mm{' (' + f'{best['geom']['plate_shape']*1e3:.0f} between)' if best['geom'].get('plate_shape') else ''}, {best['plate_m']:.0f} kg/m, {best['tape']:.1f} kg/m of tape.")
    # the launch heat in the narrower plates: the lift and load plates at their candidates
    print()
    print(f"  The launch heat in the lift and load plates (P4's model) for the plate candidates:")
    theva2 = theva; sup, lift = fb.load_air()
    heats = []
    for po, pt in ((0.2, 0.20), (0.1, 0.20), (0.1, 0.10), (0.0, 0.10), (0.1, 0.05)):
        pw = 0.8 + 2 * po
        with contextlib.redirect_stdout(io.StringIO()):
            hres = part_heat(plate_t_lift=pt, plate_w_lift=pw, plate_t_load=pt, plate_w_load=pw, label=f"{pw:.1f} m, {pt*1e3:.0f} mm")
        r30 = [r for r in hres["rows"] if r["u"] == 30000.0][0]; r1 = [r for r in hres["rows"] if r["u"] == 1000.0][0]
        heats.append(dict(po=po, pt=pt, pw=pw, E=hres["E_launch"], q30=2 * (r30["hyst"] + r30["eddy"]), q1=2 * (r1["hyst"] + r1["eddy"]), hold30=r30["hold"]))
        print(f"    {pw:.1f} m, {pt*1e3:3.0f} mm: {hres['E_launch']/1e9:6.1f} GJ per launch; {2*(r30['hyst']+r30['eddy']):6.0f} J per metre of ring per pass at 30 km/s, {2*(r1['hyst']+r1['eddy']):5.0f} at 1 km/s")
    print(f"  (the eddy heat of a 1-D diffusion falls with the plate's thickness as the field's penetration does; the plate's own field at a given pull comes from")
    print(f"  the air-core curve of the 1.4 m, 200 mm plate, so the narrower plates' higher peak field is not in these rows: see the table above for it)")
    out["heat"] = heats
    return out


def part_fast(quick=False, lengths=LENGTHS_FAST, rises=RISES):
    """final_bearing's 114 cases with the support frozen and the protection in every lane, for each (lift, load) coil length and rise."""
    hr("S3  THE FAST COILS: COIL LENGTH AND RISE, THE SUPPORT FROZEN, THE PROTECTION IN EVERY LANE")
    theva = bd.read_theva(); sup, lift = fb.load_air()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    with contextlib.redirect_stdout(io.StringIO()):
        A = fb.part_A(theva, sup, lift)
    print(f"  Casing {A['m_top']:,.0f} kg/m, cable {A['top']['m_cable']:,.0f} kg/m ({ring.QUADRATIC} quadratic).")
    L0 = fb.L_LIFT, fb.L_LOAD
    runs = {}
    for (L_lift, L_load) in lengths:
        for tau in rises:
            if quick and not (tau in (2e-3, 20e-3) and (L_lift, L_load) in (lengths[0], lengths[-1])):
                continue
            fb.L_LIFT, fb.L_LOAD = L_lift, L_load
            key = f"lift {L_lift:g} m, load {L_load:g} m, {tau*1e3:g} ms"
            t0 = time.time()
            fname = f"bearing_stageB_fast_L{L_lift:g}_{L_load:g}_tau{tau*1e3:g}.txt"
            try:
                gaps, records, txt = bl.run_cases(tau, True, A, theva, sup, lift, tabs, fname, prot_all=True, head=f"[Stage B fast coils: {key}]")
            finally:
                fb.L_LIFT, fb.L_LOAD = L0
            runs[key] = dict(gaps=gaps, summary=bl.parse_summary(txt), tau=tau, L_lift=L_lift, L_load=L_load, seconds=time.time() - t0)
            s = runs[key]["summary"]
            print(f"  ran {key}: pair {s['e1']} / {s['e30']} mm at 1 / 30 km/s, one track {s['o30']}, smallest single-failure gap {s['gap']} mm ({s['case'][:40]}), {time.time()-t0:.0f} s", flush=True)
    print()
    print(f"  {'lift / load m':>14} {'rise ms':>7} {'1 km/s pair mm':>14} {'30 km/s pair mm':>15} {'30 km/s one mm':>14} {'smallest single-failure gap mm':>30} {'in':<44} {'single contacts':>15} {'band margin mm':>14}")
    rows = []
    for key, r in runs.items():
        s = r["summary"]
        singles = [g for g in r["gaps"] if bl.is_single(g[0])]
        cs = [g[0] for g in singles if g[4]]
        try:
            e1 = float(s["e1"]); e30 = float(s["e30"]); o30 = float(s["o30"]); gap = float(s["gap"])
        except ValueError:
            e1 = e30 = o30 = gap = float("nan")
        margin = 20.0 - max(e1, e30)
        rows.append(dict(L_lift=r["L_lift"], L_load=r["L_load"], tau=r["tau"], e1=e1, e30=e30, o30=o30, gap=gap, contacts=cs, margin=margin, seconds=r["seconds"]))
        print(f"  {r['L_lift']:>6g} / {r['L_load']:<5g} {r['tau']*1e3:7g} {e1:14.2f} {e30:15.2f} {o30:14.2f} {gap:30.1f} {s['case'][:44]:<44} {len(cs):15d} {margin:14.1f}")
    print(f"  (band margin = 20 mm less the larger of the pair's movements at 1 and 30 km/s; single contacts = single-failure cases in contact)")
    ok = [r for r in rows if not r["contacts"] and r["margin"] > 0]
    if ok:
        best = max(ok, key=lambda r: r["gap"] + r["margin"])
        print(f"  Best margins among the combinations with no single-failure contact (the smallest single-failure gap plus the band margin): lift {best['L_lift']:g} m, load {best['L_load']:g} m at {best['tau']*1e3:g} ms "
              f"(smallest single-failure gap {best['gap']:.1f} mm, band margin {best['margin']:.1f} mm).")
    return dict(rows=rows)


# =============================================================================
# S3b  THE BIAS, AND  S4  THE WHOLE BEARING
# =============================================================================

BIAS = 10e3                                   # N/m per face: the standing pull the load windings and the lift tracks hold against each other
CHOSEN = dict(L_lift=50.0, L_load=25.0, tau=20e-3)   # replaced by the fast sweep's best when its json is present


def chosen_combination():
    rows = []
    for name in ("bearing_stageB_fast.json", "bearing_stageB.json"):
        p = os.path.join(ring.out_dir(), name)
        if os.path.exists(p) and not rows:
            with open(p, encoding="utf-8") as f:
                rows = json.load(f).get("fast", {}).get("rows", [])
    if rows:
        ok = [r for r in rows if not r["contacts"] and r["margin"] > 0]
        if ok:
            best = max(ok, key=lambda r: r["gap"] + r["margin"])
            return dict(L_lift=best["L_lift"], L_load=best["L_load"], tau=best["tau"], source="the fast sweep")
    return dict(CHOSEN, source="the default (the book's lengths at 20 ms)")


def part_fasttable(quick=False):
    """The turn voltage and the supply rating of each coil length and rise (Part 1a's method: one series turn, the loop flux at the track's
    full ampere-turns over the coil's length and the rise; the stored energy per coil is the track's per metre times the length over the
    20 families), and the sweep's choice."""
    hr("S3  THE FAST COILS: THE TURN VOLTAGE AND THE SUPPLY RATING OF EACH LENGTH AND RISE, AND THE CHOICE")
    r_sup, r_lift = be.load_records()
    rows = []
    print(f"  Lift track at {r_lift['F_req']/1e3:.1f} kN/m ({r_lift['NI_req']/1e3:.0f} kA per leg, {r_lift['W_req']/1e3:.1f} kJ per metre of track stored); the load track of layout (A) at 198 kN/m is the same within a percent.")
    print(f"  {'coil m':>6} {'E per coil MJ':>13} " + " ".join(f"{'V turn, MVA at ' + f'{t*1e3:g}' + ' ms':>22}" for t in RISES))
    Phi = 2 * r_lift["W_req"] / r_lift["NI_req"]
    for L in (12.5, 25.0, 50.0):
        E = r_lift["W_req"] * L / fb.N_FAM
        cells = []
        for t in RISES:
            V = Phi * L / t; VA = 2 * E / t
            rows.append(dict(L=L, tau=t, E=E, V=V, VA=VA)); cells.append(f"{V/1e3:8.1f} kV {VA/1e6:7.0f} MVA")
        print(f"  {L:6g} {E/1e6:13.3f} " + " ".join(f"{c:>22}" for c in cells))
    print(f"  (Part 1a: 14.2 kV per turn of a 50 m lift coil at 2 ms; the book's Equation 7.5: 301 MVA for the 50 m coil at 2 ms, 27 at the 22 ms planned ramp.)")
    ch = chosen_combination()
    print(f"  The sweep's choice: lift {ch['L_lift']:g} m, load {ch['L_load']:g} m at {ch['tau']*1e3:g} ms ({ch['source']}).")
    return dict(rows=rows, chosen=ch)


def part_bias(quick=False):
    hr("S3b  THE BIAS: A STANDING CURRENT IN THE LOAD WINDINGS AGAINST A STANDING PULL OF THE LIFT TRACKS")
    r_sup, r_lift = be.load_records()
    F_l = F_LOAD_ONE; NI_l = r_lift["NI_req"] * math.sqrt(F_l / r_lift["F_req"])        # the load track of layout (A), the lift design
    print(f"  A driven winding's pull goes as the square of its own current: dF/dNI = 2 F / NI, zero at zero.  The load track of layout (A) ({NI_l/1e3:.0f} kA for {F_l/1e3:.0f} kN/m):")
    print(f"  {'standing pull per face':>22} {'standing kA':>11} {'share of cap':>12} {'kA for +1 % of W (0.95 kN/m)':>28} {'kA for +5 %':>11} {'dF/dNI N/m per At':>17}")
    for B in (0.0, 2e3, 5e3, 10e3, 20e3, 40e3):
        NI_B = NI_l * math.sqrt(B / F_l) if B else 0.0
        d1 = NI_l * math.sqrt((B + 0.01 * pz.W) / F_l) - NI_B; d5 = NI_l * math.sqrt((B + 0.05 * pz.W) / F_l) - NI_B
        slope = 2 * B / NI_B if B else 0.0
        print(f"  {B/1e3:19.0f} kN/m {NI_B/1e3:11.1f} {NI_B/NI_l*100:11.0f} % {d1/1e3:28.1f} {d5/1e3:11.1f} {slope:17.2f}")
    print(f"  (the hand estimate's 23 kA for 1 percent of a track's load is the first row's kind of number: from zero, a small pull takes a large current; from a")
    print(f"  standing pull of {BIAS/1e3:.0f} kN/m per face it takes {NI_l*math.sqrt((BIAS + 0.01*pz.W)/F_l)/1e3 - NI_l*math.sqrt(BIAS/F_l)/1e3:.1f} kA more, and the slope is {2*BIAS/(NI_l*math.sqrt(BIAS/F_l)):.2f} N/m per ampere-turn.)")
    print(f"  The bias costs {BIAS/1e3:.0f} kN/m of each track's cap, a standing {NI_l*math.sqrt(BIAS/F_l)/1e3:.0f} kA in the load and the lift windings (no loss in the tape; the plates see a")
    print(f"  steady field), and {2*BIAS/1e3:.0f} kN/m of standing pull on the cable from each pair, equal and opposite.  In the dynamics models (`bias` in the point and beam")
    print(f"  models) both faces hold the bias and a net command moves them to bias + c/2 and bias - c/2, so a small correction starts from a working current.")
    print(f"  The planned pulses are shaped for the response from zero (the feed-forward's impulse matching), so a bias makes their first-order response faster")
    print(f"  than the shaping expects: the planned passes move a little more and the unplanned corrections a little less (S4's runs).")
    return dict(bias=BIAS, NI_bias=NI_l * math.sqrt(BIAS / F_l))


def part_whole(quick=False):
    ch = chosen_combination()
    hr(f"S4  THE WHOLE BEARING: THE 114 CASES OF THE CHOSEN COMBINATION (lift {ch['L_lift']:g} m, load {ch['L_load']:g} m, {ch['tau']*1e3:g} ms, from {ch['source']})")
    theva = bd.read_theva(); sup, lift = fb.load_air()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    with contextlib.redirect_stdout(io.StringIO()):
        A = fb.part_A(theva, sup, lift)
    r_sup, r_lift = be.load_records(); c = bl.part1_constants(r_sup)
    print(f"  Casing {A['m_top']:,.0f} kg/m, cable {A['top']['m_cable']:,.0f} kg/m ({ring.QUADRATIC} quadratic, the present cable).")
    dFdNI = (bo.F_of_NI(r_sup, c["NI_s0"] + 100) - bo.F_of_NI(r_sup, c["NI_s0"] - 100)) / 200.0
    rate = 2 * dFdNI * bl.switched_circuit(c, 1000.0)["rate_NI"]            # the slow coils at 1 kV: 1.7 kN/m per s on the pair
    L0 = fb.L_LIFT, fb.L_LOAD
    configs = [("support frozen", dict(freeze=True, prot_all=True)),
               ("support frozen, bias", dict(freeze=True, prot_all=True, bias=BIAS)),
               ("slow coils acting between events (1 kV)", dict(freeze=False, prot_all=True, rate=rate))]
    runs = {}
    for name, kw in configs:
        if quick and "bias" in name:
            continue
        fb.L_LIFT, fb.L_LOAD = ch["L_lift"], ch["L_load"]
        key = name; t0 = time.time()
        fname = "bearing_stageB_whole_" + name.split(" (")[0].replace(" ", "").replace(",", "_") + ".txt"
        try:
            gaps, records, txt = bl.run_cases(ch["tau"], A=A, theva=theva, sup=sup, lift=lift, tabs=tabs, log_name=fname, head=f"[Stage B whole bearing: {name}]", **kw)
        finally:
            fb.L_LIFT, fb.L_LOAD = L0
        runs[key] = dict(gaps=gaps, summary=bl.parse_summary(txt), seconds=time.time() - t0)
        print(f"  ran {key}: {len(gaps)} cases, {time.time()-t0:.0f} s -> {fname}", flush=True)
    print()
    print(f"  {'run':<42} {'1 km/s pair mm':>14} {'30 km/s pair mm':>15} {'30 km/s one mm':>14} {'smallest single-failure gap mm':>30} {'in':<44} {'single contacts':>15}")
    out = dict(chosen=ch, rate=rate, summaries={})
    for key, r in runs.items():
        s_ = r["summary"]; singles = [g for g in r["gaps"] if bl.is_single(g[0])]; cs = [g[0] for g in singles if g[4]]
        out["summaries"][key] = dict(s_, contacts_single=cs)
        print(f"  {key:<42} {s_['e1']:>14} {s_['e30']:>15} {s_['o30']:>14} {s_['gap']:>30} {s_['case'][:44]:<44} {len(cs):15d}")
    keys = list(runs)
    if len(keys) > 1:
        for k in keys[1:]:
            print(); print(f"  {keys[0]} against {k}:")
            out[f"compare|{k}"] = bl.compare_cases(runs[keys[0]]["gaps"], runs[k]["gaps"], "frozen", k.split(" (")[0][:12])
    # the smaller gaps of P2: the single-failure cases whose excursion exceeds 50 and 25 mm
    base = runs[keys[0]]["gaps"]
    singles = [g for g in base if bl.is_single(g[0])]
    o50 = [g for g in singles if min(g[2], g[3]) < 0.050]; o75 = [g for g in singles if min(g[2], g[3]) < 0.075]
    print()
    print(f"  S4.1b  For P2 (the 75 and 50 mm gap): of the {len(singles)} single-failure cases, {len(o75)} would touch at a 75 mm gap and {len(o50)} at a 50 mm gap")
    print(f"    (at 50 mm: {'; '.join(f'{g[0][:44]} {min(g[2], g[3])*1e3:.0f} mm' for g in o50) or 'none'}).")
    out["smaller_gaps"] = dict(n_single=len(singles), touch_75=[g[0] for g in o75], touch_50=[g[0] for g in o50])
    # the failure list
    print()
    print(f"  S4.2  What fails when any one part fails (the case gaps at the chosen combination, support frozen; the circuit items from S1 and B.1):")
    base = runs[keys[0]]["gaps"]
    def gap_of(sub):
        hits = [g for g in base if sub in g[0]]
        return "; ".join(f"{min(g[2], g[3])*1e3:.1f} mm{' C' if g[4] else ''}" for g in hits[:3]) if hits else "n/a"
    items = [("a support coil (one of ten): quench or open", "the nine neighbours keep the bundle's flux (A.5): 0.004 % of force; the dump moves its current in 0.3 s (S1.4)"),
             ("a support coil's switch strip stuck closed (will not open)", "the coil stays persistent and cannot be charged or dumped: its nine neighbours are trimmed around it (B.1 (f): a shorted turn, 1.6 % per step); on a quench the heaters of A.5"),
             ("a strip stuck open (will not close)", "the coil free-wheels through D2 at 1 V: 0.47 mA/s, 1 % in 18 min; the dump takes it into the neighbours and the coil waits for repair"),
             ("a diode: D1 open", "the coil cannot be charged; it holds; its neighbours take over its trim"), ("a diode: D1 shorted", "the secondary drives the coil both ways: no charge, the coil's current rings at the mains frequency through the strip: open the strip and dump"),
             ("a diode: D2 open", "the strip must close before the source stops, else the coil has no free-wheel path and the reverse half-cycle goes into the coil"), ("a diode: D2 shorted", "a short across the coil through the diode's metal: the coil holds (a milliohm short decays it over weeks, B.1 (b)), cannot be charged, and the dump branch is bypassed: a quench needs the heaters; the coil waits for repair"),
             ("a transformer (primary or core)", "no charge; the coil holds; the bundle's other nine carry its trim; the coil-end heat unchanged"),
             ("a driven coil's supply (one of the 20 families)", gap_of("one coil dead")),
             ("one lift track lost at 30 km/s", gap_of("one lift track lost at 30")), ("one load track lost (layout (A): the other load track and the support)", "the book's rule: the heaviest parked spacecraft about 2,480 t; the passages below orbital velocity on one load track: " + gap_of("sled passage 1 km/s, one track")),
             ("one support track lost, both windings, at rest / as the zone arrives", gap_of("1. one support track loses both windings at rest") + " / " + gap_of("1. one support track loses both windings as the zone")),
             ("both support tracks lost", gap_of("both support tracks lost as the zone arrives")), ("a gap sensor +5 mm / -5 mm", gap_of("sensor error +5 mm") + " / " + gap_of("sensor error -5 mm")),
             ("a false full lift pull 50 ms / a false load pull 50 ms", gap_of("false full lift pull 50 ms") + " / " + gap_of("false full load pull 50 ms")), ("a pull 50 ms late / early", gap_of("pull 50 ms late") + " / " + gap_of("pull 50 ms early")),
             ("the spread pull with no sled, beam", gap_of("spread pull with no sled, beam")), ("ring power lost", "the support coils hold for years (S4.3); the gap is lost in about 135 ms without a running gap loop (B.1): the loops' stores decide")]
    for part_, what in items:
        print(f"    {part_:>72}: {what}")
    # no power
    print()
    L_c = bl.switched_circuit(c, 1000.0)["L_common"]
    print(f"  S4.3  How long the bearing holds with no power at all: the support coils are persistent through the superconducting strip, with no diode drop and")
    print(f"  no switch resistance (a superconducting joint of 1 nohm would give {L_c/1e-9/3.156e7:,.0f} years); the driven windings' currents last as long as their supplies'")
    print(f"  free-wheel paths; nothing holds the gap without the fast loops: the attractive track's e-folding time is {math.sqrt(pz.MU/(2*0.33e6))*1e3:.0f} ms (Part 1's constant-flux stiffness),")
    print(f"  so the answer is the time the loops' cold stores can run the gap control, Parts 2 to 4's number, not a property of the coils.")
    # the mass on the casing
    print()
    tape_sup = r_sup["tape12"]; tape_load_A = r_lift["tape3"] * r_lift["NI_req"] * math.sqrt(F_LOAD_ONE / r_lift["F_req"]) / r_lift["NI_req"]; tape_lift = r_lift["tape3"]
    tape_book = 2 * (r_sup["tape12"] + r_sup["tape3_446"]) + 2 * r_lift["tape3"]
    tape_A = 2 * tape_sup + 2 * tape_load_A + 2 * tape_lift
    print(f"  S4.4  The mass on the casing (the Part 3b basis of the book's 1,030 kg/m: tape times two for formers, cryostats and supports):")
    print(f"    the book: support {tape_sup:.0f} + load {r_sup['tape3_446']:.0f} + lift {tape_lift:.0f} kg/m per side, {tape_book:.0f} kg/m of tape, {2*tape_book:,.0f} kg/m;")
    print(f"    layout (A): support {tape_sup:.0f} + load {tape_load_A:.0f} (the lift design at {F_LOAD_ONE/1e3:.0f} kN/m) + lift {tape_lift:.0f}, {tape_A:.0f} kg/m of tape, {2*tape_A:,.0f} kg/m; the chargers, diodes, strips and")
    print(f"    dump resistors of the slow coils (twenty sets per km, a few kilograms each: under 0.1 kg/m), the fast coils' supplies and stores (Parts 2 to 4);")
    print(f"    the casing side falls by {2*(tape_book - tape_A):.0f} kg/m, {2*(tape_book - tape_A)*8.6:,.0f} kg/m of cable at 8.6 per kilogram.  The cable's hardware rises by the two load plates (S6).")
    out["mass"] = dict(tape_book=tape_book, tape_A=tape_A, casing_book=2 * tape_book, casing_A=2 * tape_A)
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
    hr("THE BEARING'S COILS, STAGE B (2026-11-12)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km, {ring.QUADRATIC} quadratic: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m with hardware, side {ring.CABLE_SIDE:.2f} m{' (quick)' if quick else ''}.")
    res = {}
    cache = os.path.join(ring.out_dir(), "bearing_stageB.json" if part == "all" else f"bearing_stageB_{part}.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            res = json.load(f)
    def save():
        with open(cache, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (o.tolist() if isinstance(o, np.ndarray) else (bool(o) if isinstance(o, np.bool_) else str(o))))
    parts = {}
    for name, fn in PARTS.items():
        if part in (name, "all"):
            res[name] = fn(quick); save()
    hr(f"DONE ({time.time()-t_all:.0f} s)")


PARTS = {"slow": part_slow, "fast": part_fast, "layout": part_layout, "gap": part_gap, "table74": part_table74, "heat": part_heat, "plates": part_plates, "bias": part_bias, "whole": part_whole, "fasttable": part_fasttable}

if __name__ == "__main__":
    main()
