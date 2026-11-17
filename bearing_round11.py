#!/usr/bin/env python3
"""
Round 11 (2026-11-17), item 1: the failure cases as the design has them, in a model with length along the track (gap_regulator_track.py).
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Brief: Vol_III/Claude_Code_instructions_round11_2026-11-17.md.

    python bearing_round11.py --part=model      the stiffnesses, the lengths that follow, the model against round 10 with nothing lost, the mesh
    python bearing_round11.py --part=cases      each unit of failure lost inside the pull's window and known before the pass, 1, 10, 20, 30 km/s
    python bearing_round11.py --part=outage     the longest outage of one track a pass can cross (5,000 t, 500 t, the empty sled)
    python bearing_round11.py --part=outage1    the same for the load track at 1 km/s
    python bearing_round11.py --part=early      the pull 50 ms early: the point model with one signed current and with a current per track

    --jobs=N  processes (default 4).   --soft  also the runs with a tenth of the torsional stiffness (cases, outage).

Writes runs/round11_logs/round11_*.json; the tables go to stdout (redirect to runs/round11_logs/round11_*.txt).
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

import bearing_round7 as r7

LOGS = os.path.join(r7.OUT, "round11_logs")
_S = {}


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    json.dump(obj, open(os.path.join(LOGS, name), "w", encoding="utf-8"), indent=1, default=float)


def hr(t=""):
    print("\n" + "=" * 150 + ("\n  " + t if t else "") + "\n" + "=" * 150)


def S_():
    if not _S:
        import bearing_round10 as r10
        _S.update(r10.setup()); _S["h"] = {}
    return _S


def hist(u, M_sc=None, **kw):
    """The point lane's histories at one place: the staggered sets' command (the reference path), one coil's command, the sled's load."""
    S = S_(); key = (u, M_sc, tuple(sorted(kw.items())))
    if key not in S["h"]:
        pz = S["pz"]
        if M_sc is not None:
            kw = dict(kw, M_sc=M_sc)
        res, h = pz.point_pulsed([S["lane"](u, label="sets", **kw), S["lane"](u, label="coil", n_coils=1, **kw)], u, S["tabs"], t_settle=3.0)
        F = h["F_cmd"][0]; big = np.abs(F) > 0.01 * np.abs(F).max() if np.abs(F).max() > 0 else np.zeros(len(F), dtype=bool)
        S["h"][key] = dict(t=h["t"], q=h["q"][0], q_pred=h["q_pred"][0], F_cmd=F, F_coil=h["F_cmd"][1], t_c0=float(h["t"][big][0]) if big.any() else 0.0, point=res[0])
    return S["h"][key]


def one(c):
    """One run.  c: u (m/s), M_sc (kg or None), track ('face': the track the plan uses, on the - side), sets (None: all), Lz (m, 0: the coils over
    one point), t_rel ('start' of the pull, 'half', 'arrival', 'mid' of the spacecraft, or None: from before the pass), known, plan, scale (torsion),
    L, dx, settle."""
    import gap_regulator_track as gt
    S = S_(); u = c["u"]; h = hist(u, c.get("M_sc"))
    st = gt.stiffness(torsion_scale=c.get("scale", 1.0))
    Lz = c.get("Lz", 0.0); L = c.get("L") or max(3000.0, Lz + 2250.0); dx = c.get("dx", 6.25)
    face = 2 if (h["F_cmd"].min() < -h["F_cmd"].max()) else 3                       # lift - or load -
    tr = {"face": face, "lift": 2, "load": 3}[c.get("track", "face")]
    tl = None
    if c.get("unit", "none") != "none":
        g = gt.grid(L, dx)
        when = {None: -np.inf, "start": h["t_c0"], "half": h["t_c0"] / 2, "arrival": 0.0, "mid": 50.0 / u}[c.get("t_rel")]
        tl = gt.outage(g, tr, L / 2 - Lz / 2, L / 2 + Lz / 2, sets=c.get("sets"), t=when)
    t0 = time.time()
    o = gt.run(h["t"], h["q"], h["q_pred"], h["F_cmd"], h["F_coil"], u, S["reg_side"], st, S["c"]["lift1"], S["k"] * S["c"]["lift1"], S["X_C"], S["X_S"], L=L, dx=dx, t_lost=tl,
               known=c.get("known", False), plan=c.get("plan", "moved"), couple=S["c"]["cp"], t0=min(h["t_c0"], 0.0) - 0.15, t1=100.0 / u + c.get("settle", 1.5))
    o.pop("snap", None)
    o.update(case=c, track=gt.TRACKS[tr], seconds=time.time() - t0, planned=float(np.abs(h["F_cmd"]).max()), t_c0=h["t_c0"])
    return o


def run_all(cases, jobs):
    t0 = time.time()
    if jobs > 1:
        with Pool(jobs) as p:
            out = p.map(one, cases, chunksize=1)
    else:
        out = [one(c) for c in cases]
    print(f"  ({len(cases)} runs, {time.time()-t0:.0f} s)")
    return out


def fmt(o):
    pk = max(o["pk_coil_share"])
    return (f"{('CONTACT %.0f ms' % (o['t_contact']*1e3)) if o['contact'] else 'holds':>15} {o['gap']*1e3:7.1f} {o['move_c']*1e3:9.1f} {o['roll']*1e3:9.2f} {o['x_gap']:8.0f} {o['t_gap']*1e3:8.0f} "
            f"{pk*100:7.0f}%")


HEAD = f"{'':>15} {'gap mm':>7} {'center mm':>9} {'roll mrad':>9} {'where m':>8} {'when ms':>8} {'coil':>8}"
NOTE = ("  (gap: the smallest gap at either edge of the cable's face, anywhere in the stretch; center: the casing's largest movement on the centerline; where, when: the place of the\n"
        "  smallest gap from the middle of the outage and its time after the spacecraft's front passes that place; coil: the hardest worked coil's pull as a share of its rating)")


def part_model(jobs):
    import gap_regulator_track as gt
    S = S_(); pz = S["pz"]; st = gt.stiffness(); rs = S["reg_side"]
    hr("THE MODEL ALONG THE TRACK: THE STIFFNESSES AND WHERE THEY COME FROM")
    print(f"  Casing in bending: the beam lanes' truss, {st['truss']} (10 m deep, chords with a tenth of the casing's {pz.MU:,.0f} kg/m, diagonals with a twentieth):")
    print(f"      EI = {st['EI']:.3e} N m^2, GA = {st['GA']:.3e} N, shear wave {st['c_s']:,.0f} m/s, bar wave {st['c_b']:,.0f} m/s.")
    print(f"  Casing in torsion, a first figure: a thin-walled square tube of the casing's outline, {st['side']:.0f} m by {st['side']:.0f} m [the outline is the 10 m of the design rule taken as")
    print(f"  height and width, UNVERIFIED as before], each of the four walls a lattice with half the shear stiffness of the beam lanes' web: G t = (GA / 2) / {st['side']:.0f} m = {st['Gt']:.3e} N/m")
    print(f"  (the same as {st['Gt']/ (500e9/2.6) *1e3:.1f} mm of a solid wall with G = E / 2.6 = 192 GPa [a comparison only]).  Bredt: GJ = 4 A^2 / (perimeter / G t) = side^3 G t = {st['GJ']:.3e} N m^2.")
    print(f"  Casing's roll inertia {st['I_cas']:,.0f} kg m^2 per m (the mass on the perimeter of the 10 m square, radius of gyration {gt.gr.R_GYR:.2f} m [UNVERIFIED]); torsion wave {math.sqrt(st['GJ']/st['I_cas']):,.0f} m/s.")
    print(f"  Cable in heave: a string under its tension {st['T']:.3e} N, wave {math.sqrt(st['T']/pz.MU_C):,.0f} m/s.  Cable in roll: the stiffness its tension gives, T a^2 / 6 = {st['GJ_cab']:.3e} N m^2")
    print(f"  for the {gt.ring.CABLE_SIDE:.2f} m square (derivation in gap_regulator_track.py); inertia {st['I_cab']:,.0f} kg m^2 per m; its own shear stiffness left out.")
    k_h = 2 * rs["mass"] * rs["om"]**2; k_r = 2 * rs["mass"] * rs["om"]**2 * S["X_C"]**2
    print(f"\n  The lengths that follow, with the regulator's stiffness (per side m om^2 = {rs['mass']*rs['om']**2/1e6:.2f} MN/m per m, at +-{S['X_C']:.2f} m):")
    print(f"      roll: the regulators' stiffness in roll is 2 m om^2 X^2 = {k_r/1e6:.0f} MN m per rad per m.  A twist dies away along the track over sqrt(GJ / k) = {math.sqrt(st['GJ']/k_r):.0f} m")
    print(f"      (a tenth of the stiffness: {math.sqrt(st['GJ']/10/k_r):.0f} m).  Where one side has no coil that can pull the right way, the regulators give no stiffness in roll, and the twist has")
    print(f"      to reach live coils at the ends of the outage: a stretch of length L with its ends held and a moment m per m along it twists m L^2 / (8 GJ) at the middle.")
    print(f"      heave: the shear length sqrt(GA / k) = {math.sqrt(st['GA']/k_h):.0f} m and the bending length (4 EI / k)^(1/4) = {(4*st['EI']/k_h)**0.25:.0f} m.")
    cases = [dict(u=uk * 1e3, label=f"nothing lost, {uk} km/s") for uk in (1, 3, 5, 10, 20, 30)]
    cases += [dict(u=30e3, dx=3.125, label="nothing lost, 30 km/s, cells 3.125 m"), dict(u=30e3, unit="track", Lz=500.0, t_rel="start", label="one lift track over 500 m lost at the pull's start, 30 km/s"),
              dict(u=30e3, unit="track", Lz=500.0, t_rel="start", dx=3.125, label="the same, cells 3.125 m"),
              dict(u=30e3, unit="track", Lz=100.0, t_rel="start", label="one lift track over 100 m lost at the pull's start, 30 km/s"),
              dict(u=30e3, unit="track", Lz=100.0, t_rel="start", dx=3.125, label="the same, cells 3.125 m")]
    out = run_all(cases, jobs)
    r10 = {1: 10.9, 3: 3.6, 5: 1.6, 10: 2.1, 20: 7.7, 30: 12.6}
    hr("THE MODEL AGAINST ROUND 10 WITH NOTHING LOST, AND THE MESH (cells 6.25 m, four to a coil)")
    print(f"  {'case':<62} {HEAD} {'off the plan mm':>15} {'round 10 moves mm':>17} {'step ms':>8}")
    for c, o in zip(cases, out):
        print(f"  {c['label']:<62} {fmt(o)} {o['err']*1e3:15.2f} {r10.get(int(c['u']/1e3), float('nan')) if 'unit' not in c and 'dx' not in c else float('nan'):17.1f} {o['dt']*1e3:8.3f}")
    print(NOTE)
    print(f"  At 1 to 5 km/s the sled is slower than the casing's shear wave ({st['c_s']:,.0f} m/s), so the casing's own stiffness carries part of the spacecraft's load to the casing ahead and behind.")
    print(f"  The plan's predicted path is that of free masses and does not have that in it, which is why the casing is further off the plan there than in round 10's models.")
    save("round11_model.json", dict(st={k: v for k, v in st.items()}, k_roll=k_r, k_heave=k_h, rows=out))


UNITS = (("one coil (25 m, one set of one track)", dict(unit="coil", sets=[0], Lz=0.0)),
         ("one set of one track over 500 m", dict(unit="set", sets=[0], Lz=500.0)),
         ("one whole track over 500 m", dict(unit="track", Lz=500.0)),
         ("one whole track over 2.25 km", dict(unit="track", Lz=2250.0)))


def part_cases(jobs, soft):
    cases = []
    for scale in ((1.0, 0.1) if soft else (1.0,)):
        for name, un in UNITS:
            for uk in (1, 10, 20, 30):
                for when in ("start", "half", "arrival", "mid"):
                    cases.append(dict(un, u=uk * 1e3, t_rel=when, known=False, scale=scale, name=name, row="lost inside the window, the plan as it was", when=when))
                for pl in (("moved", "kept", "dropped") if un["unit"] == "track" else ("moved",)):
                    cases.append(dict(un, u=uk * 1e3, t_rel=None, known=True, plan=pl, scale=scale, name=name, row="known before the pass", when=pl))
    out = run_all(cases, jobs)
    hr("EACH UNIT OF FAILURE AT 1, 10, 20 AND 30 km/s (the load track at 1 km/s, the lift track above the neutral speed; the 5,000 t spacecraft)")
    print(f"  {'unit':<38} {'torsion':>7} {'km/s':>5} {'case':<44} {'moment or plan':>14} {HEAD}")
    for c, o in zip(cases, out):
        print(f"  {c['name']:<38} {c['scale']:7.1f} {c['u']/1e3:5.0f} {c['row']:<44} {c['when']:>14} {fmt(o)}")
    print(NOTE)
    print("  (moment: start = as the planned pull starts at the place; half = half way to the spacecraft's arrival; arrival; mid = with the spacecraft half over the place.")
    print("  plan, for a whole track known to be out: moved = the other side's track takes the whole planned pull over the outage, up to its rating; kept = it keeps its own half;")
    print("  dropped = no planned pull over the outage.  For a coil or a set known to be out, the live coils over it share what the place asks for.)")
    save("round11_cases.json", dict(rows=out))


LENGTHS = (25.0, 50.0, 100.0, 200.0, 300.0, 500.0, 750.0, 1000.0, 1500.0, 2250.0)


def part_outage(jobs, soft):
    cases = []
    for scale in ((1.0, 0.1) if soft else (1.0,)):
        for M in (5000e3, 500e3, 0.0):
            for uk in (10, 20, 30):
                for Lz in LENGTHS:
                    cases.append(dict(unit="track", Lz=Lz, u=uk * 1e3, M_sc=M, t_rel="start", known=False, scale=scale, how="lost as the pull starts"))
                    for pl in ("moved", "kept"):
                        cases.append(dict(unit="track", Lz=Lz, u=uk * 1e3, M_sc=M, t_rel=None, known=True, plan=pl, scale=scale, how="known, plan " + pl))
    out = run_all(cases, jobs)
    hr("THE LONGEST OUTAGE OF ONE LIFT TRACK THAT A PASS CAN CROSS (the smallest gap in mm; C = contact)")
    hows = ("lost as the pull starts", "known, plan moved", "known, plan kept")
    res = {}
    for c, o in zip(cases, out):
        res[(c["scale"], c["M_sc"], c["u"], c["how"], c["Lz"])] = o
    summ = []
    for scale in sorted({c["scale"] for c in cases}, reverse=True):
        print(f"\n  Torsional stiffness x {scale:g}")
        print(f"  {'spacecraft':>10} {'km/s':>5} {'case':<24} " + " ".join(f"{Lz:>7.0f}" for Lz in LENGTHS) + "  m of outage" + f"   {'longest with no contact':>24} {'with 80 mm left':>15}")
        for M in (5000e3, 500e3, 0.0):
            for uk in (10, 20, 30):
                for how in hows:
                    row = [res[(scale, M, uk * 1e3, how, Lz)] for Lz in LENGTHS]
                    ok = [Lz for Lz, o in zip(LENGTHS, row) if not o["contact"]]; ok80 = [Lz for Lz, o in zip(LENGTHS, row) if o["gap"] >= 0.080]
                    # the longest length below which every shorter length also holds
                    def longest(good):
                        best = 0.0
                        for Lz in LENGTHS:
                            if Lz in good: best = Lz
                            else: break
                        return best
                    lo, l80 = longest(ok), longest(ok80)
                    print(f"  {M/1e3:8,.0f} t {uk:5d} {how:<24} " + " ".join((f"{'C':>7}" if o["contact"] else f"{o['gap']*1e3:7.1f}") for o in row) + f"   {'':>11} {lo:>24.0f} {l80:>15.0f}")
                    summ.append(dict(scale=scale, M_sc=M, u=uk * 1e3, how=how, longest=lo, longest80=l80, gaps=[o["gap"] for o in row], contact=[o["contact"] for o in row],
                                     roll=[o["roll"] for o in row], coil=[max(o["pk_coil_share"]) for o in row]))
    print(f"\n  (0 t is the empty sled.  The longest lengths are those of the list tried, {', '.join(f'{Lz:.0f}' for Lz in LENGTHS)} m: the true limit lies between the last that holds and the next.)")
    save("round11_outage.json", dict(lengths=LENGTHS, summary=summ, rows=out))


def part_outage1(jobs):
    """The same scan for the load track at 1 km/s (the 5,000 t spacecraft), both torsion figures."""
    Ls = (25.0, 50.0, 100.0, 200.0, 300.0)
    cases = []
    for scale in (1.0, 0.1):
        for Lz in Ls:
            cases.append(dict(unit="track", Lz=Lz, u=1e3, t_rel="start", known=False, scale=scale, how="lost as the pull starts", L=2000.0))
            cases.append(dict(unit="track", Lz=Lz, u=1e3, t_rel=None, known=True, plan="moved", scale=scale, how="known, plan moved", L=2000.0))
    out = run_all(cases, jobs)
    hr("THE LONGEST OUTAGE OF ONE LOAD TRACK THAT A PASS AT 1 km/s CAN CROSS (5,000 t; the smallest gap in mm, the roll in mrad; C = contact)")
    print(f"  {'torsion':>7} {'case':<24} " + " ".join(f"{Lz:>14.0f}" for Lz in Ls) + "  m of outage")
    for scale in (1.0, 0.1):
        for how in ("lost as the pull starts", "known, plan moved"):
            row = [o for c, o in zip(cases, out) if c["scale"] == scale and c["how"] == how]
            print(f"  {scale:7.1f} {how:<24} " + " ".join((f"{'C':>14}" if o["contact"] else f"{o['gap']*1e3:7.1f} {o['roll']*1e3:6.2f}") for o in row))
    save("round11_outage1.json", dict(lengths=Ls, rows=out))


def part_early(jobs):
    import bearing_round10 as r10
    S = S_(); pz = S["pz"]
    hr("THE PULL 50 ms EARLY: WHY THE POINT MODEL GAVE 10.4 mm AND THE CROSS-SECTION MODEL 51.5 mm")
    print(f"  The same lane in the point model run two ways: with one signed current for the pair (round 10: the lift tracks must fall to zero before the load tracks can start), and with")
    print(f"  a current for the lift tracks and one for the load tracks (pz.REG_SEPARATE, as the hardware and the cross-section model have).")
    print(f"  {'case':<40} {'km/s':>5} {'point, one signed current':>26} {'point, a current per track':>26} {'cross-section':>14}   smallest gap, mm")
    rows = []
    for lab, kw in (("nothing wrong", {}), ("the pull 50 ms early", dict(dt_err=-50e-3)), ("the pull 50 ms late", dict(dt_err=50e-3)), ("the pull 10 ms late", dict(dt_err=10e-3)),
                    ("the pull fired with no sled", dict(no_sled=True)), ("the pull 5 ms early", dict(dt_err=-5e-3))):
        for uk in (30, 1):
            g = []
            for sep in (False, True):
                pz.REG_SEPARATE = sep
                res, _ = pz.point_pulsed([S["lane"](uk * 1e3, label="x", **kw)], uk * 1e3, S["tabs"], t_settle=3.0)
                g.append(min(res[0]["min_gap_b"], res[0]["min_gap_t"]))
            pz.REG_SEPARATE = True
            o, _, _ = r10.section(S, uk * 1e3, label="x", **kw)
            rows.append(dict(case=lab, u=uk * 1e3, signed=g[0], separate=g[1], section=o["gap"]))
            print(f"  {lab:<40} {uk:5d} {g[0]*1e3:26.1f} {g[1]*1e3:26.1f} {o['gap']*1e3:14.1f}")
    save("round11_early.json", dict(rows=rows))


def main():
    part = "model"; jobs = 4; soft = False
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a == "--soft": soft = True
        elif a in ("-h", "--help"): print(__doc__); return
    t0 = time.time()
    if part == "model": part_model(jobs)
    elif part == "cases": part_cases(jobs, soft)
    elif part == "outage": part_outage(jobs, soft)
    elif part == "early": part_early(jobs)
    elif part == "outage1": part_outage1(jobs)
    print(f"\n  DONE ({time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
