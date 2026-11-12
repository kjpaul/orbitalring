"""The support and load windings' layout, round 2 (2026-11-12): what the flux pumps are asked to do, what they can do, and
whether a soldered persistent coil can be charged or protected at all.  Stage A of Vol_III/Claude_Code_instructions_bearing_layout_2026-11-12.md.

Run: python bearing_layout.py [--part=demand|frozen|passages|pumps|solder|length|all] [--quick]
     (ORBITAL_RING_QUADRATIC=fixed for the October cable the Part 1 report used)

A.1  demand    the volt-seconds, peak voltage and power the 1 Hz support loop asks of the flux pumps in each of final_bearing's
               114 cases (the support force history of every point-model lane, through the support track's pull curve and the
               incremental loop inductance of bearing_electrical Part 1d) and in the Part 1d passages
A.2  frozen    final_bearing's cases rerun with the support current frozen (steady_loop_hz = 0, the coils keeping their flux) at the
               2 ms rise of the book and the 20 ms rise of this round, then the longest rise that still keeps the 20 mm band and every
               single-failure gap with the support frozen (30, 50, 75, 100 ms)
A.3  passages  the Part 1d coupled passages (1 km staggered, 10 km staggered, 10 km ending together) with the loop on and frozen
A.4  pumps     N flux pumps in series per support coil: the time to trim 1 percent, to take over a parked spacecraft, and the N a
               1.5 percent trim in one second needs, at the output voltages the literature has reached
A.5  solder    a soldered support coil charged through its ends: the turn-to-turn resistance, the L / R of the spiral, the layer a
               one-day charge would need; the quench: the stored energy against the winding's heat capacity, the whole-coil copper
               decay, the hot spot's time to room temperature, the bypass a resistive layer gives, the decay of a coil through a
               10 cm normal zone on its own and on its leakage inductance to the other nine coils of the bundle, what the heaters,
               the cascade, the induction and a pump acting on the whole stack can do
A.6  length    the support coil length as a free choice (10 m to 10 km, ten per bundle, staggered and ending together): the coupling
               with a shared-loop load winding, the coil-end eddy heat in the plates (final_bearing Part C.3's model), the stored energy
               per coil, coils, pumps and joints per kilometre, the charging time through a layer that still bypasses a 10 cm normal zone,
               and the decay of a coil through such a zone into its neighbours

Nothing here changes a default of the existing scripts: the point models only record the support force (hist['S'], passage_point's
S) and passage_point takes loop_hz = 1.0 as before.  Outputs go to runs/ (ORBITAL_RING_OUT).
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
import final_bearing as fb
import sled_passage_pulsed as pz
import sled_passage_simulation as sp

MU0 = 4e-7 * math.pi
N_SUP_COILS = be.N_SUP_COILS                  # 10 staggered persistent coils per bundle
L_SUP_COILS = be.L_SUP_COILS                  # (1 km, 10 km)
W_HALF = be.W_HALF                            # one support track's share of the casing weight
GAP = be.GAP
RISE_SWEEP = (30e-3, 50e-3, 75e-3, 100e-3)    # the rises tried with the support frozen, after the 2 and 20 ms runs
DEMAND_WINDOW = 0.01                          # s: the support command is read every 10 ms for the rate (the loop's own time scale is 1 s)
W_TAPE = 0.012; T_CU = 20e-6; T_HAST = 100e-6; T_AG = 2e-6; T_REBCO = 3e-6     # the 12 mm support tape's stack [copper 20 um in all as Part 1;
                                                                            # Hastelloy 100 um, silver and REBCO as THEVA tapes: Barusco et al. 2024]
RHO_HAST = 1.24e-6                            # ohm m, Hastelloy C-276, flat from 4 to 300 K (Lu, Choi, Zhou 2008, J. Appl. Phys. 103 064908, via the supermagnet fit plot)
DENS_CU, DENS_HAST = 8960.0, 8890.0           # kg/m^3 [handbook]
RRR_CU = 50.0
# the ideal (phonon) resistivity of copper against temperature, nohm m [handbook table, UNVERIFIED], the residual part from RRR
RHO_IDEAL_T = np.array([10, 20, 30, 40, 50, 60, 70, 77, 80, 100, 150, 200, 250, 300.0])
RHO_IDEAL = np.array([0.0002, 0.008, 0.06, 0.22, 0.52, 0.97, 1.52, 1.94, 2.10, 3.4, 6.9, 10.3, 13.5, 16.8]) * 1e-9
# NIST cryogenic material properties, OFHC copper specific heat fit (4 to 300 K), log10 cp = sum a_k (log10 T)^k
CP_CU_COEF = (-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0.0)
T_OP, T_SHARE, T_C = 70.0, 85.0, 92.0         # operating, a current-sharing temperature for the half-Ic tape [UNVERIFIED], Tc of REBCO

# the flux pumps of the literature: (label, time-averaged dc output V, the current reached A, the load inductance H, how it was read)
PUMPS = [("Geng et al. 2019, transformer-rectifier, 77 K: 1150 A into 7-9 uH in 200 s", 8e-6 * 1150 / 200, 1150.0, 8e-6, "read in full (arXiv 1901.07983)"),
         ("Mallett et al. 2023, self-rectifier, 77 K: 320 A into 2.5 uH in under 3 s, 60 mV peak per cycle", 2.5e-6 * 320 / 3, 320.0, 2.5e-6, "read in full (arXiv 2210.12999)"),
         ("Bumby et al. 2016, external-rotor dynamo, 77 K: > 30 A into 1.97 mH over the first 500 s, joints 3.2 uohm", 1.97e-3 * 30 / 300, 30.0, 1.97e-3, "read in full (MT-24 preprint)"),
         ("Jiang et al. (in Zhou et al. 2023), dynamo at 600 Hz, 77 K: 58 A into 2.4 mH, measured saturated voltage 1.2 mV", 1.2e-3, 58.0, 2.4e-3, "read in full (KCL accepted manuscript)"),
         ("Iftikhar et al. 2025, conduction-cooled transformer-rectifier: 160 A into 32 uH at 77 K, 32 mV on the bridge while it pumps; the charging time is not stated", float("nan"), 160.0, 32e-6, "read in full")]
V_PUMP_TABLE = ((0.27e-3, "0.27 mV: the best time-averaged output read (Mallett et al. 2023)"), (1.2e-3, "1.2 mV: the dynamo's measured saturated voltage (Zhou et al. 2023)"),
                (5e-3, "5 mV: the brief's hand assumption"), (60e-3, "60 mV: the highest peak output read, as if it were held all the time (an upper bound no device reaches on average)"))
# turn-to-turn contact resistivities, ohm m^2, with their sources (read = the paper was read in full here)
R_CONTACTS = [("soldered lap joint, Cu-stabilized tape, 30 nohm cm2 (Zhang & Li, IEEE TASC; abstract only)", 30e-9 * 1e-4),
              ("soldered joint, Ag-dispersed, 9 to 18 nohm cm2 (IEEE TASC 2022; abstract only)", 15e-9 * 1e-4),
              ("the brief's soldered figure, 1 uohm cm2", 1e-6 * 1e-4),
              ("dry copper-to-copper, NI coil, 6.3 uohm cm2 (Sohn, Ha, Kim 2018, read)", 6.3e-6 * 1e-4),
              ("dry copper-to-copper, NI coil, 71 uohm cm2 (Wang, Hahn et al. 2013, read)", 71e-6 * 1e-4),
              ("the brief's dry figure, 100 uohm cm2", 100e-6 * 1e-4),
              ("stainless co-wind, MI coil, 725 uohm cm2 (Sohn, Ha, Kim 2018, read)", 724.8e-6 * 1e-4),
              ("oxidized copper surfaces, 35 mohm cm2 at 4.2 K (Lu et al., read)", 35e-3 * 1e-4),
              ("as-received stainless co-wind, native oxide, about 100,000 uohm cm2 at 77 K, worn to under 100 at 4.2 K by pressure cycling (Lu et al., read)", 0.1 * 1e-4),
              ("the brief's tailored layer, 11 ohm cm2", 11.0 * 1e-4)]


def hr(t=""):
    print(); print("=" * 150); print(f"  {t}"); print("=" * 150)


# =============================================================================
# the support track's electrical constants from Part 1
# =============================================================================

def part1_constants(r_sup):
    """The incremental loop inductances of bearing_electrical Part 1d (runs/bearing_electrical_5mm.json) and the support coil's turns."""
    p = os.path.join(ring.out_dir(), "bearing_electrical_5mm.json")
    if not os.path.exists(p):
        raise SystemExit(f"{p} missing: run 'python bearing_electrical.py --part=1d' first (its D.1 inductances are used here)")
    with open(p, encoding="utf-8") as f:
        d = json.load(f)["1d"]
    L_lin = np.array(d["L_lin"]); L_ops = {float(k) * 1e3: np.array(v) for k, v in d["L_ops"].items()}
    NI_pts = sorted(L_ops); Lp = {NI: float(L_ops[NI][0, 0] + L_ops[NI][0, 1]) for NI in NI_pts}     # flux per loop-turn per ampere-turn, both loops at NI
    NI_s0 = bo.NI_of_F(r_sup, W_HALF); NI_cap = r_sup["NI_req"]
    n_t = NI_cap / r_sup["I12"] / N_SUP_COILS                                                       # turns per loop per coil
    return dict(L_lin=L_lin, L_ops=L_ops, Lp=Lp, Lp_lin=float(L_lin[0, 0] + L_lin[0, 1]), NI_s0=NI_s0, NI_cap=NI_cap, n_t=n_t, I_tape0=NI_s0 / NI_cap * r_sup["I12"])


def Lp_of(c, NI):
    """Incremental flux per loop-turn per ampere-turn at NI, both loops moving together: the linear value at zero, Part 1d's incremental
    values at 141, 322 and 493 kA, straight lines between."""
    xs = [0.0] + sorted(c["Lp"]); ys = [c["Lp_lin"]] + [c["Lp"][k] for k in sorted(c["Lp"])]
    return float(np.interp(NI, xs, ys))


def dphi(c, NI1, NI2, n=200):
    """Flux change per loop-turn per metre between two ampere-turn levels, Wb/m."""
    x = np.linspace(NI1, NI2, n); return float(np.trapezoid([Lp_of(c, v) for v in x], x))


def energy_per_m(c, NI):
    """Stored energy per metre of track at NI in both loops: 2 x integral of NI' L'(NI') dNI'."""
    x = np.linspace(0.0, NI, 400); return 2 * float(np.trapezoid([v * Lp_of(c, v) for v in x], x))


# =============================================================================
# A.1  the demand on the pumps from the support force histories
# =============================================================================

class SupportCurve:
    """Vectorized inverse of the support track's pull curve (bearing_optimize's log-log interpolation with the quadratic tail)."""

    def __init__(self, r_sup):
        self.x, self.y = bo.pull_curve(r_sup)

    def NI_of_F(self, F):
        F = np.asarray(F, float)
        out = np.where(F <= self.y[0], self.x[0] * np.sqrt(np.maximum(F, 0.0) / self.y[0]), 0.0)
        hi = F > self.y[0]
        if hi.any():
            out[hi] = np.exp(np.interp(np.log(np.minimum(F[hi], self.y[-1])), np.log(self.y), np.log(self.x)))
        return out


def pump_demand(t, S, curve, c, window=DEMAND_WINDOW):
    """What one lane's support force history asks of the pumps: the ampere-turns per track from S / 2 through the pull curve, read
    every `window` seconds; the total variation (the volt-seconds the pumps must put in, both directions), the peak rate, the extremes."""
    dt = float(t[1] - t[0]); k = max(1, int(round(window / dt)))
    NI = curve.NI_of_F(np.asarray(S, float)[::k] / 2)
    d = np.diff(NI)
    rate = d / (k * dt)
    j = int(np.argmax(np.abs(rate))) if len(rate) else 0
    Lp = Lp_of(c, float(NI[0]))
    return dict(NI0=float(NI[0]), NI_min=float(NI.min()), NI_max=float(NI.max()), TV=float(np.abs(d).sum()), rate_pk=float(rate[j]) if len(rate) else 0.0,
                NI_at_pk=float(NI[j]) if len(rate) else float(NI[0]), Lp=Lp,
                Vs_turn_m=Lp * float(np.abs(d).sum()),                       # volt-seconds per loop-turn per metre of coil
                V_turn_m=Lp * (abs(float(rate[j])) if len(rate) else 0.0),   # peak volts per loop-turn per metre of coil
                P_m=2 * Lp * (float(NI[j]) * abs(float(rate[j])) if len(rate) else 0.0))   # peak W per metre of track (both loops)


def demand_table(records, curve, c, title, n_show=None):
    """Print the demand of every recorded lane and return the rows."""
    n_t = c["n_t"]
    print(f"  {title}")
    print(f"  Per loop-turn and per metre of coil: Vs = L' x (total variation of the support ampere-turns), V = L' x peak rate, L' the incremental flux per")
    print(f"  ampere-turn at the lane's level; per coil (both loops in series, {n_t:,.0f} turns per loop): x 2 x {n_t:,.0f} x the coil length; power per metre of")
    print(f"  track: 2 L' NI dNI/dt (both loops) at the peak; the rate is read every {DEMAND_WINDOW*1e3:.0f} ms.")
    print(f"  {'case':<70} {'u km/s':>6} {'NI kA min':>9} {'max':>7} {'swing kA':>8} {'rate kA/s':>9} {'Vs/turn.m':>9} {'V/turn.m':>8} {'1 km coil: Vs':>13} {'V peak':>8} {'10 km: Vs':>10} {'V peak':>8} {'W/m pk':>8}")
    rows = []
    for lab, u, t, S, extra in records:
        d = pump_demand(t, S, curve, c)
        Vs1 = d["Vs_turn_m"] * 1e3 * 2 * n_t; V1 = d["V_turn_m"] * 1e3 * 2 * n_t
        rows.append(dict(case=lab, u=u, **d, Vs_coil_1km=Vs1, V_coil_1km=V1, Vs_coil_10km=10 * Vs1, V_coil_10km=10 * V1, **extra))
    rows_show = rows if n_show is None else sorted(rows, key=lambda r: -r["TV"])[:n_show]
    for r in rows_show:
        print(f"  {r['case'][:70]:<70} {r['u']/1e3:6.0f} {r['NI_min']/1e3:9.1f} {r['NI_max']/1e3:7.1f} {r['TV']/1e3:8.1f} {r['rate_pk']/1e3:9.2f} {r['Vs_turn_m']:9.4f} {r['V_turn_m']:8.4f} "
              f"{r['Vs_coil_1km']:13.0f} {r['V_coil_1km']:8.0f} {r['Vs_coil_10km']:10.0f} {r['V_coil_10km']:8.0f} {r['P_m']:8.0f}")
    return rows


# =============================================================================
# A.2  final_bearing's cases with the support loop on and frozen
# =============================================================================

SINGLE_EXCLUDE = ("no per-side", "no protection", "support lost too")


def is_single(label):
    return not any(s in label for s in SINGLE_EXCLUDE) and not label.startswith("2.")


def run_cases(tau, freeze, A, theva, sup, lift, tabs, log_name, prot_all=False):
    """final_bearing's Part B at the rise tau, the support loop on or frozen in every lane (point and beam), recording every point
    lane's support force history.  prot_all: the lanes that final_bearing runs without the 10 mm gap protection (the F.1 passages and
    tolerances, which lean on the 1 Hz support loop to arrest the drift after the pass) get it, at their own track rating, as the F.2
    cases have it; the 'no protection' what-ifs keep none.  Returns (all_gaps, records, text)."""
    import compensation_tracks as ct
    orig_point, orig_beam = pz.point_pulsed, pz.beam_pulsed
    records = []

    def point_wrapped(lanes, u, tabs_, **kw):
        if freeze:
            for ln in lanes:
                ln["steady_loop_hz"] = 0.0
        if prot_all:
            for ln in lanes:
                if ln.get("prot") is None and "no protection" not in (ln.get("label") or ""):
                    cap_up = ln["hu_width"] * bd.P12 if ln.get("hu_width") else 2 * fb.F_LOAD_TRACK
                    cap_down = ln["hd_width"] * bd.P12 if ln.get("hd_width") else 2 * fb.F_LIFT_TRACK
                    ln["prot"] = ct.prot(0.010, 1e-3, cap_up=cap_up, cap_down=cap_down)
        res, hist = orig_point(lanes, u, tabs_, **kw)
        for i, (ln, r) in enumerate(zip(lanes, res)):
            records.append((ln.get("label", "") or f"lane {i}", u, hist["t"], hist["S"][i].copy(), dict(exc=r["max_exc_all"], gap_b=r["min_gap_b"], gap_t=r["min_gap_t"], contact=r["contact"])))
        return res, hist

    def beam_wrapped(*args, **kw):
        if freeze:
            kw["loop_hz"] = 0.0
        return orig_beam(*args, **kw)

    tau0 = bd.TAU_HW, fb.TAU_HW
    pz.point_pulsed, pz.beam_pulsed = point_wrapped, beam_wrapped
    bd.TAU_HW = tau; fb.TAU_HW = tau
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            gaps = fb.part_B(A, theva, sup, lift, tabs)
    finally:
        pz.point_pulsed, pz.beam_pulsed = orig_point, orig_beam
        bd.TAU_HW, fb.TAU_HW = tau0
    txt = buf.getvalue()
    with open(os.path.join(ring.out_dir(), log_name), "w", encoding="utf-8") as f:
        f.write(f"# final_bearing Part B, rise {tau*1e3:g} ms, support loop {'frozen (steady_loop_hz = 0, beam loop_hz = 0)' if freeze else 'on (1 Hz)'}{', the 10 mm gap protection in every lane' if prot_all else ''}\n")
        f.write(txt)
    return gaps, records, txt


def parse_summary(txt):
    """The headline numbers of a Part B log: the pair excursions at 1 and 30 km/s, the one-track excursion at 30 km/s, the smallest
    single-failure gap and its case, the contact cases."""
    out = dict(e1="?", e30="?", o30="?", gap="?", case="?", contacts="?")
    for line in txt.splitlines():
        s = line.strip()
        if s.startswith("Smallest gap in any single failure:"):
            out["gap"] = s.split(":")[1].split("mm")[0].strip(); out["case"] = line.split("'")[1] if "'" in line else "?"
        if "cases in all.  Contact in:" in line:
            out["contacts"] = line.split("Contact in:")[1].strip().rstrip(".")
        parts = s.split()
        if len(parts) >= 6 and parts[0] in ("1", "30") and line.startswith("   ") and "sled" not in line:
            try:
                float(parts[1]); float(parts[2]); float(parts[3])
                if parts[0] == "1" and out["e1"] == "?": out["e1"] = parts[2]
                if parts[0] == "30" and out["e30"] == "?": out["e30"], out["o30"] = parts[2], parts[3]
            except ValueError:
                pass
    return out


def compare_cases(base, other, label_b, label_o):
    """Case-by-case comparison of two Part B gap lists (same order); prints the rows that move by more than 0.5 mm or change contact."""
    print(f"  {'case':<80} {label_b + ' gap mm':>18} {label_o + ' gap mm':>18} {'change':>7} {'contact':>12}")
    moved = 0
    for (lab, model, gb, gt, cb), (lab2, model2, gb2, gt2, cb2) in zip(base, other):
        g1 = min(gb, gt); g2 = min(gb2, gt2)
        if abs(g2 - g1) > 0.5e-3 or cb != cb2:
            moved += 1
            print(f"  {lab[:80]:<80} {g1*1e3:18.1f} {g2*1e3:18.1f} {(g2-g1)*1e3:+7.1f} {('yes' if cb else 'no') + ' -> ' + ('yes' if cb2 else 'no'):>12}")
    print(f"  ({moved} of {len(base)} cases move by more than 0.5 mm or change contact; the rest are unchanged)")
    s1 = [g for g in base if is_single(g[0])]; s2 = [g for g in other if is_single(g[0])]
    w1 = min(s1, key=lambda g: min(g[2], g[3])); w2 = min(s2, key=lambda g: min(g[2], g[3]))
    c1 = [g[0] for g in base if g[4]]; c2 = [g[0] for g in other if g[4]]
    print(f"  Smallest single-failure gap: {label_b} {min(w1[2], w1[3])*1e3:.1f} mm ({w1[0]}); {label_o} {min(w2[2], w2[3])*1e3:.1f} mm ({w2[0]}).")
    print(f"  Contact cases: {label_b} {len(c1)} ({'; '.join(c1) or 'none'}); {label_o} {len(c2)} ({'; '.join(c2) or 'none'}).")
    return dict(worst_b=(w1[0], min(w1[2], w1[3])), worst_o=(w2[0], min(w2[2], w2[3])), contacts_b=c1, contacts_o=c2, moved=moved)


def part_A2(curve, c, quick=False):
    hr("A.1 / A.2  FINAL_BEARING'S CASES: THE DEMAND ON THE PUMPS, THEN THE SUPPORT FROZEN, THEN THE LONGEST RISE WITH IT FROZEN")
    theva = bd.read_theva(); sup, lift = fb.load_air()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    with contextlib.redirect_stdout(io.StringIO()):
        A = fb.part_A(theva, sup, lift)
    print(f"  Casing {A['m_top']:,.0f} kg/m, cable {A['top']['m_cable']:,.0f} kg/m (final_bearing's on-top version, as Part 1 ran it).")
    out = {}
    configs = [(2e-3, False, False), (2e-3, True, False), (2e-3, True, True), (20e-3, False, False), (20e-3, True, False), (20e-3, True, True)] + ([] if quick else [(tau, True, True) for tau in RISE_SWEEP])
    runs = {}
    for tau, freeze, prot_all in configs:
        key = f"{tau*1e3:g} ms, {'frozen + prot' if (freeze and prot_all) else ('frozen' if freeze else 'loop on')}"
        t0 = time.time()
        gaps, records, txt = run_cases(tau, freeze, A, theva, sup, lift, tabs, f"bearing_layout_cases_tau{tau*1e3:g}_{'frozen_prot' if (freeze and prot_all) else ('frozen' if freeze else 'on')}.txt", prot_all=prot_all)
        runs[key] = dict(gaps=gaps, records=records, summary=parse_summary(txt), tau=tau, freeze=freeze, prot_all=prot_all, seconds=time.time() - t0)
        print(f"  ran {key}: {len(gaps)} cases, {len(records)} point lanes recorded, {time.time()-t0:.0f} s", flush=True)
    # A.1 the demand, from the loop-on runs
    for tau in (2e-3, 20e-3):
        key = f"{tau*1e3:g} ms, loop on"
        print()
        rows = demand_table(runs[key]["records"], curve, c, f"A.1  Demand on the pumps, rise {tau*1e3:g} ms, every point lane of final_bearing's Part B (the beam-model lanes are not recorded):")
        tv = np.array([r["TV"] for r in rows]); pk = np.array([abs(r["rate_pk"]) for r in rows]); vs1 = np.array([r["Vs_coil_1km"] for r in rows]); v1 = np.array([r["V_coil_1km"] for r in rows])
        big = max(rows, key=lambda r: r["TV"]); fast = max(rows, key=lambda r: abs(r["rate_pk"]))
        print(f"  Over the {len(rows)} lanes: swing of the support ampere-turns median {np.median(tv)/1e3:.1f} kA, largest {tv.max()/1e3:.1f} kA ({big['case']}); peak rate median {np.median(pk)/1e3:.2f} kA/s,")
        print(f"  largest {pk.max()/1e3:.1f} kA/s ({fast['case']}).  Per 1 km coil (both loops, {c['n_t']:,.0f} turns per loop): volt-seconds median {np.median(vs1):.0f}, largest {vs1.max():.0f} V s; peak")
        print(f"  voltage median {np.median(v1):.0f} V, largest {v1.max():.0f} V; ten times those for a 10 km coil.  Peak power per metre of track largest {max(r['P_m'] for r in rows):,.0f} W/m")
        print(f"  ({max(r['P_m'] for r in rows)*1e3/1e6:.1f} MW per km of track).")
        out[f"demand_{tau*1e3:g}ms"] = dict(rows=[{k: v for k, v in r.items()} for r in rows])
    # A.2 frozen against loop on
    for tau in (2e-3, 20e-3):
        print()
        print(f"  A.2  Rise {tau*1e3:g} ms: the support loop on (the book's model) against the support frozen, nothing else changed:")
        out[f"compare_{tau*1e3:g}ms"] = compare_cases(runs[f"{tau*1e3:g} ms, loop on"]["gaps"], runs[f"{tau*1e3:g} ms, frozen"]["gaps"], "loop on", "frozen")
        print()
        print(f"  A.2  Rise {tau*1e3:g} ms: the support loop on (the book's model) against the support frozen with the 10 mm gap protection in every lane:")
        out[f"compare_prot_{tau*1e3:g}ms"] = compare_cases(runs[f"{tau*1e3:g} ms, loop on"]["gaps"], runs[f"{tau*1e3:g} ms, frozen + prot"]["gaps"], "loop on", "frozen+prot")
    print()
    print(f"  A.2  Summary of every run (the pair's casing movement at 1 and 30 km/s, one track at 30 km/s, the smallest single-failure gap, contacts):")
    print(f"  {'run':<20} {'1 km/s pair mm':>14} {'30 km/s pair mm':>15} {'30 km/s one mm':>14} {'smallest single-failure gap mm':>30} {'in':<52} {'contacts':<70} {'s':>5}")
    for key, r in runs.items():
        s = r["summary"]
        print(f"  {key:<20} {s['e1']:>14} {s['e30']:>15} {s['o30']:>14} {s['gap']:>30} {s['case'][:52]:<52} {s['contacts'][:70]:<70} {r['seconds']:5.0f}")
    out["summaries"] = {k: dict(r["summary"], seconds=r["seconds"]) for k, r in runs.items()}
    if not quick:
        # the longest rise with the support frozen that keeps the 20 mm band at 30 km/s and every single-failure gap open
        ok = []
        for key, r in runs.items():
            if not (r["freeze"] and r["prot_all"]):
                continue
            s = r["summary"]
            try:
                e30 = float(s["e30"]); o30 = float(s["o30"]); gap = float(s["gap"])
            except ValueError:
                continue
            singles = [g for g in r["gaps"] if is_single(g[0])]
            contact_single = [g[0] for g in singles if g[4]]
            ok.append((r["tau"], e30, o30, gap, contact_single))
        print()
        print(f"  Longest rise with the support frozen and the protection in every lane (band 20 mm on the pair's movement at 30 km/s, no single-failure contact):")
        for tau, e30, o30, gap, cs in sorted(ok):
            print(f"    {tau*1e3:5.0f} ms: pair {e30:5.1f} mm, one track {o30:5.1f} mm, smallest single-failure gap {gap:5.1f} mm, single-failure contacts {len(cs)}{' (' + '; '.join(cs) + ')' if cs else ''}"
                  f"  -> {'keeps the band' if e30 <= 20.0 and not cs else 'FAILS'}")
        out["rise_sweep"] = [dict(tau=t, e30=e, o30=o, gap=g, contacts=cs) for t, e, o, g, cs in sorted(ok)]
    return out


# =============================================================================
# A.3  the Part 1d passages, loop on and frozen
# =============================================================================

def part_A3(r_sup, curve, c, quick=False):
    hr("A.3  THE PART 1D COUPLED PASSAGES WITH THE SUPPORT LOOP ON AND FROZEN")
    L_lin = c["L_lin"]
    systems = {}
    for L_c in L_SUP_COILS:
        systems[f"{L_c/1e3:.0f} km, {N_SUP_COILS} staggered"] = be.CoilSystem(L_c, L_lin, n_fam=N_SUP_COILS, cells_per_coil=100 if quick else 200)
    systems["10 km, ending together"] = be.CoilSystem(max(L_SUP_COILS), L_lin, n_fam=1, cells_per_coil=100 if quick else 200)
    kernels = {}
    for key, cs in systems.items():
        xc = cs.x[cs.n_cells // 2]
        dNI_s, delta = cs.solve(cs.stretch(xc, cs.dx, 1.0))
        xi = cs.x - xc; m = np.abs(xi) <= 3 * cs.L_c
        kernels[key] = (xi[m], (dNI_s[0] / cs.dx)[m])
    print(f"  The point model of bearing_electrical Part 1d (D.6) at 1, 3 and 5 km/s, the load pair's pull from bearing_optimize.PASS_PULLS, the 20 ms")
    print(f"  command ramp and the {be.TAU_HW*1e3:.0f} ms lag, protection as there; the support's 1 Hz loop on, then frozen (loop_hz = 0: the support pull follows")
    print(f"  only the coupled coil currents).  The demand columns are the loop-on run's.")
    print(f"  {'u km/s':>6} {'coils':<24} {'protection':>12} {'loop':>7} {'casing mm':>9} {'floor gap':>9} {'ceiling':>8} {'dip kA':>7} {'contact':>8} {'swing kA':>8} {'rate kA/s':>9} {'1 km coil Vs':>12} {'V pk':>6} {'10 km: Vs':>9} {'V pk':>6}")
    rows = []
    for u_kms, F_pair, kind in bo.PASS_PULLS:
        if kind != "load":
            continue
        u = u_kms * 1e3
        for lh in (1.0, 0.0):
            base = be.passage_point(u, be.L_LOAD, F_pair * 1e3, kind, r_sup=r_sup, coupled=False, loop_hz=lh)
            d = pump_demand(base["t"], base["S"], curve, c) if lh else None
            rows.append(dict(u=u, coils="-", prot="on", loop=lh, **{k: base[k] for k in ("exc", "gap_b", "gap_t", "contact", "t_contact")}, dNI_min=0.0, demand=d))
            dd = (f"{d['TV']/1e3:8.1f} {d['rate_pk']/1e3:9.2f} {d['Vs_turn_m']*1e3*2*c['n_t']:12.0f} {d['V_turn_m']*1e3*2*c['n_t']:6.0f} {d['Vs_turn_m']*1e4*2*c['n_t']:9.0f} {d['V_turn_m']*1e4*2*c['n_t']:6.0f}" if d else "")
            print(f"  {u_kms:6.0f} {'uncoupled':<24} {'on':>12} {('1 Hz' if lh else 'frozen'):>7} {base['exc']*1e3:9.1f} {base['gap_b']*1e3:9.1f} {base['gap_t']*1e3:8.1f} {0.0:7.1f} {('yes' if base['contact'] else 'no'):>8} {dd}")
        for key, cs in systems.items():
            L_c = cs.L_c; stagger = L_c / cs.n_fam
            for plabel, pon, pf in ((f"local x {max(0.0, 1 - 2 * be.L_LOAD / stagger):.2f}", True, max(0.0, 1 - 2 * be.L_LOAD / stagger)), ("off", False, 0.0)):
                for lh in (1.0, 0.0):
                    r = be.passage_point(u, be.L_LOAD, F_pair * 1e3, kind, kernel=kernels[key], L_c=L_c, coupled=True, r_sup=r_sup, prot=pon, prot_factor=pf, t_settle=1.0 + L_c / u, loop_hz=lh)
                    d = pump_demand(r["t"], r["S"], curve, c) if lh else None
                    rows.append(dict(u=u, coils=key, prot=plabel, loop=lh, **{k: r[k] for k in ("exc", "gap_b", "gap_t", "contact", "t_contact", "dNI_min")}, demand=d))
                    dd = (f"{d['TV']/1e3:8.1f} {d['rate_pk']/1e3:9.2f} {d['Vs_turn_m']*1e3*2*c['n_t']:12.0f} {d['V_turn_m']*1e3*2*c['n_t']:6.0f} {d['Vs_turn_m']*1e4*2*c['n_t']:9.0f} {d['V_turn_m']*1e4*2*c['n_t']:6.0f}" if d else "")
                    print(f"  {u_kms:6.0f} {key:<24} {plabel:>12} {('1 Hz' if lh else 'frozen'):>7} {r['exc']*1e3:9.1f} {r['gap_b']*1e3:9.1f} {r['gap_t']*1e3:8.1f} {r['dNI_min']/1e3:7.1f} {('yes' if r['contact'] else 'no'):>8} {dd}", flush=True)
    print(f"  (dip = the largest loss of support ampere-turns at the point from the coupling; the loop-on rows' 'swing' is what the 1 Hz loop commanded on top of")
    print(f"  that, which is what the pumps would have to deliver within the passage)")
    return rows


# =============================================================================
# A.4  N pumps in series per coil
# =============================================================================

def part_A4(c):
    hr("A.4  N FLUX PUMPS IN SERIES IN EACH SUPPORT COIL")
    n_t = c["n_t"]; NI_s0 = c["NI_s0"]; NI_cap = c["NI_cap"]
    print(f"  What the literature has delivered as a time-averaged dc voltage into a coil (L dI/dt over the charge, or the measured value):")
    for lab, V, I, L, how in PUMPS:
        print(f"    {(f'{V*1e3:7.3f} mV' if V == V else '  (n/a)  '):>10}  {I:6.0f} A  {lab}  [{how}]")
    print(f"  The transformer-rectifier pumps put tens of millivolts on the bridge during their pumping half-cycle (Iftikhar 32 mV, Mallett 60 to 65 mV peak)")
    print(f"  and the dynamos about a millivolt; averaged over the cycle and the charge, no device read here exceeds about 1 mV into its load.")
    print(f"  Any pump in series with a support coil carries the coil's current, {c['I_tape0']:.0f} A per tape at the normal share and {NI_cap/N_SUP_COILS/n_t:.0f} A at the cap:")
    print(f"  within what the pumps above have carried.")
    Vs_1pct = 0.01 * NI_s0 * Lp_of(c, NI_s0)                       # per loop-turn per metre
    Vs_take = dphi(c, NI_s0, NI_cap)
    Vs_15 = 0.015 * NI_s0 * Lp_of(c, NI_s0)
    print()
    print(f"  Volt-seconds per loop-turn per metre of coil: 1 percent of {NI_s0/1e3:.0f} kA {Vs_1pct:.5f}; the takeover 141 to {NI_cap/1e3:.0f} kA {Vs_take:.4f} (the flux curve through")
    print(f"  the saturating plate); 1.5 percent {Vs_15:.5f}.  Per coil (both loops, {n_t:,.0f} turns per loop): x {2*n_t:,.0f} x the coil length.")
    print(f"  {'coil':>6} {'pump mV':>8} {'N':>6} {'1 percent':>14} {'takeover 141 to 322 kA':>24} {'N for 1.5 percent in 1 s':>24} {'pumps per km of track (both tracks)':>36}")
    rows = []
    for L_c in L_SUP_COILS:
        vs1 = Vs_1pct * L_c * 2 * n_t; vst = Vs_take * L_c * 2 * n_t; vs15 = Vs_15 * L_c * 2 * n_t
        coils_per_km = N_SUP_COILS * 1000.0 / L_c * 2            # staggered coils starting per km, both tracks
        for V, vlab in V_PUMP_TABLE:
            for N in (1, 100, 10000):
                t1 = vs1 / (N * V); tt = vst / (N * V); N15 = vs15 / (V * 1.0)
                rows.append(dict(L_c=L_c, V=V, N=N, t_1pct=t1, t_take=tt, N_15=N15, Vs_1pct=vs1, Vs_take=vst))
                print(f"  {L_c/1e3:4.0f} km {V*1e3:8.2f} {N:6d} {fmt_t(t1):>14} {fmt_t(tt):>24} {N15:24,.0f} {N*coils_per_km:36,.0f}")
    print(f"  (the N for the one-second trim is the same for every N column: it is the N the trim needs)")
    print()
    print(f"  The pumps' own cold loss.  Mallett et al. 2023 (read): the bridge of the self-rectifier dissipates of the order of 1 J per cycle at 10 Hz")
    print(f"  while putting about 0.1 W of dc power into the load (320 A at 0.27 mV): about 100 W of heat per watt delivered [read from the paper's")
    print(f"  Fig. 8 axis, approximate].  Ghabeli et al. 2021 (read, a model of a dynamo): 135 mW of ac loss in the stator tape at 25 Hz for a device")
    print(f"  whose output is of the order of a millivolt at tens of amperes: heat of the order of the delivered power or more.  A pump's loss is a")
    print(f"  standing loss of its moving field (f x hysteresis), not a share of its output: it runs whether or not the coil needs volt-seconds.")
    for L_c in L_SUP_COILS:
        coils_per_km = N_SUP_COILS * 1000.0 / L_c * 2
        for N in (1, 100, 10000):
            for q_pump, qlab in ((0.15, "0.15 W per pump, a dynamo [Ghabeli's ac loss]"), (10.0, "10 W per pump, a self-rectifier [Mallett's bridge loss]")):
                print(f"    {L_c/1e3:.0f} km coils, N = {N:6d}: {N*coils_per_km*q_pump:14,.0f} W of cold heat per km of track at {qlab}")
    print(f"  For scale: Part 1c's lead heat for a driven coil's 21 kA pair was 1.9 kW per coil, 93 GW around the ring, and was judged unaffordable.")
    return rows


def fmt_t(s):
    if s < 120: return f"{s:.1f} s"
    if s < 7200: return f"{s/60:.1f} min"
    if s < 2 * 86400: return f"{s/3600:.1f} h"
    if s < 2 * 365.25 * 86400: return f"{s/86400:.1f} d"
    return f"{s/365.25/86400:,.0f} yr"


# =============================================================================
# A.5  the soldered coil: charging through the turns, and the quench
# =============================================================================

def rho_cu(T):
    """Copper resistivity, ohm m: the ideal part from the table and the residual from RRR."""
    return float(np.interp(T, RHO_IDEAL_T, RHO_IDEAL)) + RHO_IDEAL[-1] / RRR_CU


def cp_cu(T):
    """NIST OFHC copper specific heat fit, J/(kg K), 4 to 300 K."""
    lt = math.log10(T); return 10 ** sum(a * lt**k for k, a in enumerate(CP_CU_COEF))


def tape_heat_capacity(T):
    """Heat capacity per metre of 12 mm tape, J/(m K): copper from the NIST fit, Hastelloy taken with copper's specific heat per kg
    [UNVERIFIED: Lu, Choi, Zhou 2008 is the source for Hastelloy; its value was not read]; silver and REBCO neglected (5 um of 125)."""
    return W_TAPE * (T_CU * DENS_CU + T_HAST * DENS_HAST) * cp_cu(T)


def tape_resistance_per_m(T):
    """Normal-state resistance per metre of the tape: copper and Hastelloy in parallel (the REBCO and silver neglected)."""
    g = T_CU * W_TAPE / rho_cu(T) + T_HAST * W_TAPE / RHO_HAST
    return 1.0 / g


def hot_spot_time(I, T_end=300.0, T0=T_OP):
    """Adiabatic heating of a normal length of tape carrying I: the time from T0 to T_end."""
    Ts = np.linspace(T0, T_end, 2000); t = 0.0
    for a, b in zip(Ts[:-1], Ts[1:]):
        Tm = 0.5 * (a + b); t += tape_heat_capacity(Tm) * (b - a) / (I**2 * tape_resistance_per_m(Tm))
    return t


def leakage_inductances(tr, c, stride=3):
    """The ten coils of one bundle as groups of tapes in the bundle lattice: their per-metre inductance matrix from free-space filaments
    (one filament per tape with the strip GMD, the tape-to-tape term through sqrt(d^2 + gmd^2) so that stacked tapes stay positive
    definite [an approximation of the strip-to-strip GMD]), the return leg mirrored and the plate as an image; grouped by layers
    (block: a coil is one tenth of the bundle's height) and interleaved (every tenth tape).  Returns the self and leakage inductance
    per metre per coil for the two groupings, with the tapes subsampled by `stride` (each sample standing for stride tapes)."""
    xs, ys = tr.lattices[0]
    x_ret = tr.legs[1][0]; y_plate = tr.gap
    order = np.argsort(ys); xs = xs[order][::stride]; ys = ys[order][::stride]
    n = len(xs); gmd = 0.2235 * W_TAPE
    k = (be.cf.PLATE_MU - 1) / (be.cf.PLATE_MU + 1)
    def lnd(x1, y1, x2, y2):
        d2 = (x1[:, None] - x2[None, :])**2 + (y1[:, None] - y2[None, :])**2
        return 0.5 * np.log(d2 + gmd**2)
    Xr = 2 * x_ret - xs; Yi = 2 * y_plate - ys
    M = lnd(xs, ys, Xr, ys) - lnd(xs, ys, xs, ys) + lnd(Xr, ys, xs, ys) - lnd(Xr, ys, Xr, ys)
    M += k * (lnd(xs, ys, Xr, Yi) - lnd(xs, ys, xs, Yi)) + k * (lnd(Xr, ys, xs, Yi) - lnd(Xr, ys, Xr, Yi))
    M *= MU0 / (2 * math.pi) * stride**2
    out = {}
    for name, groups in (("block", np.minimum((np.arange(n) * N_SUP_COILS) // n, N_SUP_COILS - 1)), ("interleaved", np.arange(n) % N_SUP_COILS)):
        G = np.zeros((n, N_SUP_COILS)); G[np.arange(n), groups] = 1.0
        Lc = G.T @ M @ G                                        # per metre, each coil's tapes in series
        leak = []
        for i in range(N_SUP_COILS):
            o = [j for j in range(N_SUP_COILS) if j != i]
            leak.append(Lc[i, i] - Lc[i, o] @ np.linalg.solve(Lc[np.ix_(o, o)], Lc[o, i]))
        out[name] = dict(L_self=float(np.mean(np.diag(Lc))), L_leak=float(np.mean(leak)), L_sum=float(Lc.sum()))
    return out


def part_A5(r_sup, c, h=0.005):
    hr("A.5  THE SOLDERED SUPPORT COIL: CAN IT BE CHARGED, AND WHAT HAPPENS IN A QUENCH")
    n_t = c["n_t"]; NI_s0 = c["NI_s0"]; NI_cap = c["NI_cap"]; I0 = c["I_tape0"]; I_cap = NI_cap / N_SUP_COILS / n_t
    out = {}
    # 1. the coil as a circuit
    print(f"  Support coil: {n_t:,.0f} turns per loop, both loops in series, length L_c, the 12 mm tape at {I0:.0f} A ({NI_s0/1e3:.0f} kA per loop over {N_SUP_COILS} coils) and {I_cap:.0f} A at the cap.")
    print(f"  Flux per loop-turn per ampere-turn (both loops moving together): linear {c['Lp_lin']*1e6:.4f} uWb/m, incremental {Lp_of(c, NI_s0)*1e6:.4f} at {NI_s0/1e3:.0f} kA, {Lp_of(c, NI_cap)*1e6:.4f} at the cap.")
    print(f"  {'coil':>6} {'L alone H':>10} {'all ten together H':>19} {'contact area per turn pair m2':>30} {'R end to end at 1 uohm cm2':>28}")
    for L_c in L_SUP_COILS:
        L1 = 2 * n_t**2 * L_c * Lp_of(c, NI_s0); L10 = N_SUP_COILS * L1
        A_pair = 4 * L_c * W_TAPE                                  # four legs of tape face per turn pair
        R = (n_t - 1) * 1e-10 / A_pair
        print(f"  {L_c/1e3:4.0f} km {L1:10.0f} {L10:19.0f} {A_pair:30.0f} {R*1e9:25.3f} nohm")
        out[f"L_{L_c:.0f}"] = dict(L1=L1, L10=L10, A_pair=A_pair)
    print(f"  (the brief's 109 H and 1,090 H count one loop of a 1 km coil; a coil has both loops in series, so its inductance is twice that and its")
    print(f"  contact area per turn pair is four tape faces, 48 m2 for 1 km; L / R is the same either way)")
    print()
    print(f"  A pump or supply at the coil's ends drives the spiral in parallel with the stack of contacts: the current settles into the spiral with")
    print(f"  tau = L / R_c, R_c the contact stack end to end (Wang, Hahn et al. 2013, the no-insulation coil's equivalent circuit), with L the")
    print(f"  inductance the pump sees: all ten coils of the bundle charged together, each sees ten times its own.  tau grows as L_c^2.")
    print(f"  {'turn-to-turn contact':<120} {'1 km coil tau':>14} {'10 km coil tau':>15}")
    for lab, r in R_CONTACTS:
        taus = []
        for L_c in L_SUP_COILS:
            L10 = N_SUP_COILS * 2 * n_t**2 * L_c * Lp_of(c, NI_s0)
            R = (n_t - 1) * r / (4 * L_c * W_TAPE)
            taus.append(L10 / R)
        print(f"  {lab:<120} {fmt_t(taus[0]):>14} {fmt_t(taus[1]):>15}")
    for L_c in L_SUP_COILS:
        L10 = N_SUP_COILS * 2 * n_t**2 * L_c * Lp_of(c, NI_s0)
        for tau, tl in ((86400.0, "one day"), (3600.0, "one hour")):
            r = tau * L10 / (n_t - 1) * 4 * L_c * W_TAPE / L10 * (L10 / L10)
            r = (L10 / tau) * 4 * L_c * W_TAPE / (n_t - 1)
            print(f"    {L_c/1e3:.0f} km coil, tau = {tl}: needs {r*1e4:,.1f} ohm cm2 between turns")
    print(f"  Laboratory check of the scaling: Wang et al.'s NI60 (60 turns, 4 mm tape, L 432 uH, R_c 534 uohm, tau 0.81 s, 71 uohm cm2) has turns of")
    print(f"  about 0.3 m; the time grows with the turn count, the square of the turn length over the tape width, and the contact resistivity.")
    # 2. the quench: energy and heat capacity
    print()
    E0 = energy_per_m(c, NI_s0); Ecap = energy_per_m(c, NI_cap)
    m_tape = r_sup["tape12"]
    cpT = cp_cu(T_OP)
    print(f"  The quench.  Stored energy per metre of track {E0/1e3:.1f} kJ/m at {NI_s0/1e3:.0f} kA and {Ecap/1e3:.1f} kJ/m at the cap (the flux curve), into {m_tape:.0f} kg/m of")
    print(f"  12 mm tape: {E0/m_tape/cpT:.1f} K and {Ecap/m_tape/cpT:.1f} K of warming at {cpT:.0f} J/(kg K) [copper's specific heat at {T_OP:.0f} K from the NIST fit, applied to the")
    print(f"  whole tape; UNVERIFIED for Hastelloy].  Spread over the whole coil a dump is harmless.")
    for L_c in L_SUP_COILS:
        L1 = 2 * n_t**2 * L_c * Lp_of(c, NI_s0)
        R_all = tape_resistance_per_m(T_OP) * 4 * L_c * n_t
        print(f"    {L_c/1e3:.0f} km coil wholly normal at {T_OP:.0f} K: {R_all/1e3:,.1f} kohm of copper (the Hastelloy in parallel).  All ten coils of the bundle normal together: each sees")
        print(f"    ten times its own inductance through the others, L / R = {N_SUP_COILS*L1/R_all*1e3:.0f} ms, and the bundle's {E0*L_c/1e6:.1f} MJ per coil length goes into {m_tape*L_c/1e3:.0f} t of tape, {E0/m_tape/cpT:.1f} K.")
        print(f"    One coil normal with nine persistent: its current leaves on the leakage inductance (below), in a millisecond, and only the leakage energy is dissipated.")
    # 3. the hot spot
    print()
    J = I0 / (T_CU * W_TAPE)
    print(f"  The hot spot.  A normal zone's current runs in the tape's copper ({T_CU*1e6:.0f} um over {W_TAPE*1e3:.0f} mm): {J/1e6:.0f} A/mm2 at {I0:.0f} A, {I_cap/(T_CU*W_TAPE)/1e6:.0f} A/mm2 at {I_cap:.0f} A.")
    print(f"  Copper resistivity {rho_cu(T_OP)*1e9:.2f} nohm m at {T_OP:.0f} K (RRR {RRR_CU:.0f}), {rho_cu(300)*1e9:.1f} at 300 K [handbook table, UNVERIFIED]; heat capacity of the tape")
    print(f"  {tape_heat_capacity(T_OP):.3f} J/(m K) at {T_OP:.0f} K, {tape_heat_capacity(300):.3f} at 300 K.  Adiabatic time for a normal zone to reach room temperature:")
    for I in (I0, I_cap, 2 * I_cap):
        print(f"    {I:6.0f} A: {hot_spot_time(I):6.2f} s to 300 K, {hot_spot_time(I, 150.0):6.2f} s to 150 K")
    out["hot_spot"] = {f"{I:.0f}": hot_spot_time(I) for I in (I0, I_cap)}
    # 4. the bypass through a resistive layer
    print()
    print(f"  A layer of resistivity r (ohm m2) between the turns lets the zone's current cross into the neighbouring turns: the current leaves the")
    print(f"  normal tape over a transfer length lambda = sqrt(r / R_s), R_s = rho / t_cu the copper's sheet resistance (1-D current transfer, the")
    print(f"  neighbours superconducting).  A zone longer than 2 lambda carries its full current only over lambda at each end; a zone shorter than")
    print(f"  that is bypassed in part only.  Half the current around a 10 cm zone needs lambda of about 5 cm.")
    Rs = rho_cu(T_OP) / T_CU
    print(f"  {'layer r':>12} {'lambda at 70 K':>14} {'at 300 K':>9} {'share of a 10 cm zone bypassed [1 - lambda / 5 cm, floored]':>60} {'1 km coil charge tau':>20}")
    for r_cm2 in (1e-6, 1e-4, 5e-3, 1e-1, 11.0, 300.0):
        r = r_cm2 * 1e-4
        lam = math.sqrt(r / Rs); lam300 = math.sqrt(r / (rho_cu(300) / T_CU))
        L10 = N_SUP_COILS * 2 * n_t**2 * 1000.0 * Lp_of(c, NI_s0); R = (n_t - 1) * r / (4 * 1000.0 * W_TAPE)
        print(f"  {r_cm2:9.0e} ohm cm2 {lam*1e2:12.2f} cm {lam300*1e2:7.2f} cm {max(0.0, 1 - lam / 0.05):60.2f} {fmt_t(L10/R):>20}")
    r_half = (0.05)**2 * Rs
    print(f"  The layer that bypasses half of a 10 cm zone at {T_OP:.0f} K is r = {r_half*1e4:.1e} ohm cm2; it charges a 1 km coil in {fmt_t(N_SUP_COILS * 2 * n_t**2 * 1000.0 * Lp_of(c, NI_s0) / ((n_t - 1) * r_half / (4 * 1000.0 * W_TAPE)))}.")
    print(f"  The 11 ohm cm2 layer of a one-day charge has lambda {math.sqrt(11e-4/Rs):.1f} m: it bypasses nothing until the zone is metres long, and by then (the")
    print(f"  propagation below) the spot has long reached room temperature.")
    # 5. decay through a 10 cm zone: own inductance and leakage to the nine others
    print()
    tr = be.Track(bo.BEST_SUPPORT, h)
    lk = leakage_inductances(tr, c)
    print(f"  A coil's current decaying through a 10 cm normal zone (the turns insulated or soldered alike: the zone is in series with the spiral).")
    print(f"  Zone resistance {tape_resistance_per_m(T_OP)*0.1*1e3:.2f} mohm at {T_OP:.0f} K, {tape_resistance_per_m(300)*0.1*1e3:.2f} mohm at 300 K.  The other nine coils of the bundle keep their flux, so the")
    print(f"  decaying coil sees only its leakage inductance to them: from the bundle lattice ({len(tr.lattices[0][0]):,} tapes per leg, every 3rd taken), the coil's")
    print(f"  own inductance per metre and its leakage when the coils are blocks of layers and when they are interleaved:")
    for name, v in lk.items():
        print(f"    {name:<12}: self {v['L_self']*1e3:.3f} mH/m, leakage {v['L_leak']*1e6:.2f} uH/m ({v['L_leak']/v['L_self']*100:.3f} percent of self); the bundle's ten coils in series {v['L_sum']*1e3:.2f} mH/m")
    L11 = float(c["L_lin"][0, 0])
    print(f"    (check: the ten coils in series are one loop of the bundle; the 2-D solver's L11 x (10 x {n_t:,.0f})^2 = {L11*(N_SUP_COILS*n_t)**2*1e3:,.0f} mH/m; the filament sum is")
    print(f"    {lk['block']['L_sum']/(L11*(N_SUP_COILS*n_t)**2)*100:.0f} percent of it, the plate as an image only; the leakage is taken as the filament model's share of the solver's own inductance)")
    out["leakage"] = lk
    print(f"  {'coil':>6} {'grouping':<12} {'own L (both loops) H':>20} {'leakage H':>10} {'tau through 10 cm at 70 K, own':>31} {'leakage':>10} {'at 300 K, own':>14} {'leakage':>10}")
    for L_c in L_SUP_COILS:
        for name, v in lk.items():
            Ls = 2 * n_t**2 * L_c * Lp_of(c, NI_s0); Ll = Ls * v["L_leak"] / v["L_self"]
            R70 = tape_resistance_per_m(T_OP) * 0.1; R300 = tape_resistance_per_m(300) * 0.1
            print(f"  {L_c/1e3:4.0f} km {name:<12} {Ls:20.0f} {Ll:10.2f} {fmt_t(Ls/R70):>31} {fmt_t(Ll/R70):>10} {fmt_t(Ls/R300):>14} {fmt_t(Ll/R300):>10}")
    print(f"  Against the hot spot's {hot_spot_time(I0):.1f} s to room temperature at {I0:.0f} A: the current does not leave the zone in time on either inductance; the")
    print(f"  spot burns unless something else removes the current, and a coil that has lost its current to its nine neighbours has not lost the")
    print(f"  bundle's ampere-turns (the flux is conserved: the pull stays), so a dead coil costs the track nothing until its neighbours are at their limit.")
    # 6. propagation, heaters, the load winding, a pump on the stack
    print()
    print(f"  Normal zone propagation.  Bonura and Senatore 2016 (read): in copper-stabilized REBCO tape the longitudinal velocity is set by the")
    print(f"  operating current almost regardless of field and temperature, 1 to 100 cm/s over 100 to 1,000 A in their Fig. 3 (SuperPower tape at 4.5 to")
    print(f"  40 K), confirmed against measurements to within 25 percent above 100 A; at {I0:.0f} A per tape it is at or below the bottom of that range")
    print(f"  [extrapolated, UNVERIFIED: below their 0.2 Ic].  Barusco et al. 2024 (read): pristine THEVA tape at 77 K, with a current-flow-diverter")
    print(f"  coating 5 to 8 times faster.  At 5 cm/s a zone crosses a 1 km coil in {1000/0.05/3600:.1f} hours; the cascade cannot spread a quench.")
    print()
    print(f"  What can take a coil from a spot to wholly normal, or take its current away, within a second:")
    E_h = lambda L_c, f, dT: m_tape * L_c / N_SUP_COILS * f * cp_cu(T_OP + dT / 2) * dT
    ratio = lk["block"]["L_leak"] / lk["block"]["L_self"]
    for L_c in L_SUP_COILS:
        L1 = 2 * n_t**2 * L_c * Lp_of(c, NI_s0); E_coil = E0 * L_c / N_SUP_COILS; L_leak = ratio * L1; E_leak = 0.5 * L_leak * I0**2
        print(f"    {L_c/1e3:.0f} km coil ({m_tape*L_c/N_SUP_COILS/1e3:.0f} t of tape, a tenth of the bundle's {E0*L_c/1e6:.1f} MJ; leakage inductance {L_leak:.1f} H holding {E_leak/1e3:.1f} kJ at {I0:.0f} A):")
        for f in (1.0, 0.1):
            R = tape_resistance_per_m(T_OP) * 4 * L_c * n_t * f
            dT_dump = E_leak / (m_tape * L_c / N_SUP_COILS * f) / cp_cu(T_OP)
            print(f"      heaters over {f:.0%} of its length to {T_C:.0f} K: {E_h(L_c, f, T_C - T_OP)/1e6:7.1f} MJ of heater energy ({E_h(L_c, f, T_C - T_OP)/E_coil:5.1f} x the coil's share of the stored energy), {E_h(L_c, f, T_C - T_OP)/1e6:.0f} MW for one")
            print(f"        second; the coil's current then leaves on its leakage inductance in L_leak / R = {L_leak/R*1e3:.2f} ms into the nine neighbours, which keep the bundle's ampere-turns; the")
            print(f"        leakage energy warms the heated length by {dT_dump:.2f} K more.  The quenched coil is then recharged; the bundle has lost nothing.")
        print(f"      a persistent-current switch and dump resistor: to bring the whole bundle down in tau needs V = L_10 I / tau at the switch:")
        for tau in (1.0, 60.0, 3600.0):
            V = N_SUP_COILS * L1 * I0 / tau
            print(f"        tau {tau:6.0f} s: {V/1e3:8.1f} kV across each coil's switch (ten coils together)")
    print(f"    the load compensation winding's pulse: Part 1 D.7 found a 287 kA pulse over 50 m swings the 1 km staggered coils' tapes by up to 137 A")
    print(f"    in the pulse, against a limit of {r_sup['ic12']:.0f} A at the design field: it drives current into the support coils but does not quench them, and")
    print(f"    what it drives in is taken out of the same coils elsewhere (D.2): it changes no coil's total current.")
    N_pen = 2 * W_TAPE * 0.5 / (MU0 * r_sup["ic12"])
    print(f"    a dynamo acting on the whole stack: a stack of N tapes screens an applied field up to about mu0 N Ic / (2 w) [slab bound]: {MU0*n_t*r_sup['ic12']/(2*W_TAPE):.1f} T for")
    print(f"    the {n_t:,.0f} turns of one coil, so a 0.5 T rotor field reaches about {N_pen:.0f} tapes; the published dynamos use one stator tape or a few in parallel")
    print(f"    (Bumby et al. 2016, Geng et al. 2019), none a soldered stack.  The turns the rotor does not reach see the pumped turns' volts across the")
    print(f"    solder: the same leakage as a pump at the ends.  [No published device of this kind was found; the bound is a hand estimate.]")
    return out


# =============================================================================
# A.6  the coil length as a free choice
# =============================================================================

COIL_LENGTHS = (10.0, 30.0, 100.0, 1000.0, 10000.0)


def coil_end_heat(L_c, n_steps, sup, split=True):
    """final_bearing Part C.3's coil-end ripple model: each coil end steps the plate field by 1/n of the loop's field (n coils per
    bundle, staggered) or by the whole of it (coils ending together), over 2 x (gap + 50 mm) at the relative speed; the eddy energy of a
    step is half a rise-and-fall cycle of the 1-D diffusion; steps per metre 2 ends x 2 loops x n / L_c.  W per metre of ring, two tracks."""
    NI_w, Bp_w = fb.curve_at(sup, W_HALF, "B_peak")
    ell = 2 * (bd.GAP0 + 0.05)
    n = n_steps if split else 1
    dB = Bp_w / n
    E_step = 0.5 * bd.eddy_cycle(dB, ell / bd.V_REL, 0.0, d_p=fb.PLATE_T, mu_r=fb.lc.FERRO_MU_R, rho_p=fb.lc.FERRO_RESISTIVITY)
    steps = 2 * 2 * n / L_c
    return 2 * E_step * (0.5 + 0.4) * bd.V_REL * steps, Bp_w


def part_A6(r_sup, c, quick=False):
    hr("A.6  WHY ARE THE SUPPORT COILS SO LONG: THE COIL LENGTH AS A FREE CHOICE")
    L_lin = c["L_lin"]; n_t = c["n_t"]; NI_s0 = c["NI_s0"]; NI_cap = c["NI_cap"]; I0 = c["I_tape0"]
    NI_l = r_sup["NI_load"]; NI_lost = bo.NI_of_F(r_sup, 2 * W_HALF) - NI_s0
    sup, lift = fb.load_air()
    out = {}
    # (a) the coupling with a shared-loop load winding
    print(f"  (a) The coupling with a shared-loop load winding (Part 1d's flux-conserving coil system): the share of the load winding's added pull that")
    print(f"  survives under its stretch after the persistent coils have kept their flux, and the support ampere-turns the stretch's neighbourhood is left")
    print(f"  with.  The decoupled layout (the figure-eight support, the load winding on the outer legs) links no flux into the support at any length:")
    print(f"  its coupling is zero by symmetry (Part 1 D.8 found no cross term in the field solve), so only the shared-loop column depends on the length.")
    rows = [("parked 5,000 t: 181 kA over 125 m", NI_cap - NI_s0, 125.0), ("lost track at rest: 58 kA over 100 m", NI_lost, 100.0), ("protection pulse: 287 kA over 50 m", NI_l, 50.0)]
    F_s0 = bo.F_of_NI(r_sup, NI_s0)
    print(f"  {'coil':>6} {'layout':<16} " + " ".join(f"{name + ': pull kN/m, share kept, min support kA':>62}" for name, _, _ in rows))
    coupling = []
    for L_c in COIL_LENGTHS:
        for n_fam, lay in ((N_SUP_COILS, "10 staggered"), (1, "ending together")):
            n_span = max(16, int(math.ceil(1200.0 / L_c)))
            cs = be.CoilSystem(L_c, L_lin, n_fam=n_fam, n_span=n_span, cells_per_coil=(100 if L_c < 100 or quick else 200))
            xc = cs.x[cs.n_cells // 2]
            cells = []
            for name, NIl, l in rows:
                dNI_s, delta = cs.solve(cs.stretch(xc, l, NIl))
                prof = dNI_s[0]; inside = np.abs(cs.x - xc) < l / 2
                s_under = NI_s0 + float(prof[inside].mean()); F_under = bo.F_of_NI(r_sup, s_under + NIl); F_book = bo.F_of_NI(r_sup, NI_s0 + NIl)
                s_min = NI_s0 + float(prof.min())
                share = (F_under - F_s0) / (F_book - F_s0)
                cells.append((F_under, share, s_min))
                coupling.append(dict(L_c=L_c, layout=lay, case=name, F_under=F_under, share=share, s_min=s_min, F_book=F_book))
            print(f"  {L_c:6.0f} {lay:<16} " + " ".join(f"{F/1e3:34.1f} {sh:12.2f} {sm/1e3:14.1f}" for F, sh, sm in cells))
    out["coupling"] = coupling
    print(f"  (share kept = (pull under the stretch - the support's own {F_s0/1e3:.0f} kN/m) / (the book's pull at the same ampere-turns - {F_s0/1e3:.0f}); coils shorter than the")
    print(f"  stretch lie wholly inside it and cancel all of it, so with shared loops a short coil leaves the load winding nothing to do)")
    # (b) coil-end eddy heat, (c) energy, (d) counts, (e) charging through a protecting layer, (f) the decay into the neighbours
    print()
    P1, Bp = coil_end_heat(1000.0, N_SUP_COILS, sup)
    print(f"  (b) to (f) against the length.  Coil-end eddy heat from final_bearing Part C.3's model (plate peak field {Bp:.3f} T at the weight share; a staggered")
    print(f"  bundle steps a tenth of it at every coil end, a bundle ending together all of it at once; {P1:.1f} W/m for the book's 1 km staggered coils), against")
    print(f"  the 82 W/m of LIM heat the brief quotes from Table 14.3.  Stored energy: the bundle's {energy_per_m(c, NI_s0)/1e3:.1f} kJ/m x L_c / 10 as a coil's share, and a coil's")
    print(f"  own L I^2 / 2.  Coils per km per track {N_SUP_COILS} x 1000 / L_c in either layout; one pump and two joints per coil at N = 1.  Charging through a layer")
    print(f"  that bypasses half of a 10 cm zone at 70 K (r = 2.3 mohm cm2, A.5), the brief's 5 mohm cm2 and Lu's oxidized copper (35 mohm cm2): tau = L_10 / R_c,")
    print(f"  falling as L_c^2.  The decay of a coil's current through a 10 cm zone at 70 K into its nine flux-conserving neighbours (the block leakage share")
    print(f"  of A.5), against the hot spot's {hot_spot_time(I0):.0f} s to room temperature at {I0:.0f} A.")
    tr = be.Track(bo.BEST_SUPPORT, 0.005)
    lk = leakage_inductances(tr, c)
    ratio = lk["block"]["L_leak"] / lk["block"]["L_self"]; ratio_i = lk["interleaved"]["L_leak"] / lk["interleaved"]["L_self"]
    R70 = tape_resistance_per_m(T_OP) * 0.1
    Lp0 = Lp_of(c, NI_s0)
    print(f"  {'coil':>6} {'eddy W/m stag.':>14} {'ending tog.':>11} {'share of 82':>11} {'E share MJ':>10} {'own LI2/2 kJ':>12} {'coils/km':>8} {'joints/km':>9} {'coils ring':>10} {'tau 2.3 mohm cm2':>16} {'5 mohm':>10} {'35 mohm':>10} {'1 pump 1.2 mV, 1 pct':>20} {'zone decay block':>16} {'interleaved':>11}")
    table = []
    for L_c in COIL_LENGTHS:
        Ps, _ = coil_end_heat(L_c, N_SUP_COILS, sup); Pe, _ = coil_end_heat(L_c, N_SUP_COILS, sup, split=False)
        E_share = energy_per_m(c, NI_s0) * L_c / N_SUP_COILS
        L1 = 2 * n_t**2 * L_c * Lp0; L10 = N_SUP_COILS * L1
        E_own = 0.5 * L1 * I0**2
        coils_km = N_SUP_COILS * 1000.0 / L_c
        taus = [L10 / ((n_t - 1) * r / (4 * L_c * W_TAPE)) for r in (2.3e-7, 5e-7, 35e-7)]
        t_pump = 0.01 * NI_s0 * Lp0 * L_c * 2 * n_t / 1.2e-3
        tz_b = ratio * L1 / R70; tz_i = ratio_i * L1 / R70
        table.append(dict(L_c=L_c, P_stag=Ps, P_end=Pe, E_share=E_share, E_own=E_own, coils_km=coils_km, tau_23=taus[0], tau_5=taus[1], tau_35=taus[2], t_pump_1pct=t_pump, tz_block=tz_b, tz_inter=tz_i))
        print(f"  {L_c:6.0f} {Ps:14.1f} {Pe:11.1f} {Ps/82:11.2f} {E_share/1e6:10.2f} {E_own/1e3:12.0f} {coils_km:8.0f} {2*coils_km:9.0f} {coils_km*bd.L_RING/1e3:10,.0f} {fmt_t(taus[0]):>16} {fmt_t(taus[1]):>10} {fmt_t(taus[2]):>10} {fmt_t(t_pump):>20} {fmt_t(tz_b):>16} {fmt_t(tz_i):>11}")
    out["table"] = table
    Vs_km = 0.01 * NI_s0 * Lp0 * 1000.0 * 2 * n_t * N_SUP_COILS
    print(f"  Pumps per kilometre: the volt-seconds a 1 percent trim needs per km of track are the same at every length, {Vs_km:,.0f} V s (2 loops x {n_t:,.0f} turns x")
    print(f"  10 coils x 0.00203 Wb/m x 1,000 m), so a trim of 1 percent per day needs {Vs_km/(1.2e-3*86400):.0f} pumps per km per track at 1.2 mV whatever the coil length;")
    print(f"  the length only sets how those pumps are shared out ({N_SUP_COILS*1000/COIL_LENGTHS[0]:.0f} coils per km at 10 m, 1 at 10 km).  The charging time through the turns' layer and")
    print(f"  the decay into the neighbours fall as L_c^2 and L_c: they are what the length decides.")
    return out


# =============================================================================

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
    hr("THE SUPPORT AND LOAD WINDINGS' LAYOUT, STAGE A (2026-11-12)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km: casing {ring.M_LOAD_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m{' (quick)' if quick else ''}.")
    r_sup, r_lift = be.load_records()
    c = part1_constants(r_sup); curve = SupportCurve(r_sup)
    print(f"  Part 1d constants: flux per loop-turn per ampere-turn {c['Lp_lin']*1e6:.4f} uWb/m linear, {Lp_of(c, c['NI_s0'])*1e6:.4f} incremental at {c['NI_s0']/1e3:.0f} kA; {c['n_t']:,.0f} turns per loop per coil.")
    res = {}
    cache = os.path.join(ring.out_dir(), "bearing_layout_A.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            res = json.load(f)
    def save():
        with open(cache, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (o.tolist() if isinstance(o, np.ndarray) else (bool(o) if isinstance(o, np.bool_) else str(o))))
    if part in ("pumps", "all"):
        res["A4"] = part_A4(c); save()
    if part in ("solder", "all"):
        res["A5"] = part_A5(r_sup, c); save()
    if part in ("length", "all"):
        res["A6"] = part_A6(r_sup, c, quick); save()
    if part in ("passages", "all"):
        res["A3"] = part_A3(r_sup, curve, c, quick); save()
    if part in ("demand", "frozen", "all"):
        res["A2"] = part_A2(curve, c, quick); save()
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
