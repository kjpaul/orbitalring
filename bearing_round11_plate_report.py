#!/usr/bin/env python3
"""
Round 11, item 1b: the tables of bearing_round11_plate.py, from its logs (runs/round11_logs/round11_plate_*.json).
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run (from Python/orbitalring): python bearing_round11_plate.py --part=report      writes runs/round11_logs/round11_plate_report.txt and .json
Nothing is solved here; every number is read from a log or is arithmetic on logged numbers that the text shows.
"""
import math

import numpy as np

import bearing_round6b as r6
import bearing_round7 as r7
import bearing_round7_parts as p7
import bearing_gap_trade as gt
import bearing_round11_plate as p
from bearing_round11_plate import Plate, RHO, SUP, FAST, HALF_FACE, M_LIM_PLATES, HW_ALLOW, BEARING_R8

BF = 1.2                       # T, the fill field of the plates carried through items 2 and 3
FACES = (7.39, 10.0, 12.0)
SIZED = "sized on the face at 1.2 T"


def allrecs():
    d = {}
    for f in ("round11_plate_tracks.json", "round11_plate_width.json", "round11_plate_width2.json", "round11_plate_width3.json", "round11_plate_alone.json", "round11_plate_design.json"):
        d.update(p.load(f) or {})
    return {k: v for k, v in d.items() if v.get("ok")}


def mass_rule(tape_sup, fasts):
    """Round 8's rule for the bearing on the casing (tools/gap_trade_report.mass_rule, bearing_round8_parts item 8): the tape of two support tracks
    and of the compensation tracks; formers, cryostats and supports at one tape mass; film capacitors at 1.5 kJ/kg for the compensation tracks'
    stored energy; 4.3 kg switch modules, one for each 2.8 kV, two switches and two diodes per coil, 0.8 coils per meter of each track (20 sets of
    25 m coils); the leads.  fasts: (how many tracks, tape kg/m, stored energy J/m, volts per coil, coil current, loops per track); a track of
    n loops has n coils where the design has one."""
    ld = dict(p7.lead(), I=r7.load("round7_fast_track.json")["I"])
    tape = 2 * tape_sup + sum(c * t for c, t, W, V, I, n in fasts)
    store = sum(c * W for c, t, W, V, I, n in fasts) / 1.5e3
    sw = sum(c * 4 * math.ceil(V / 2800.0) * 4.3 * 0.8 * n for c, t, W, V, I, n in fasts)
    leads = sum(c * 2 * 8960.0 * ld["A"] * (I / ld["I"]) * 1.0 * 0.8 * n for c, t, W, V, I, n in fasts)
    return dict(tape=tape, formers=tape, store=store, n_s=max(math.ceil(V / 2800.0) for c, t, W, V, I, n in fasts), switches=sw, leads=leads, total=2 * tape + store + sw + leads)


def fast4(r, n=1, count=4):
    return (count, r["tape"], r["W"], r["V"] / n, r["I"], n)


def get(R, kind, Wl, plate, n=None, F=None, bundle=None, h=0.005):
    if kind == "support":
        a = p.sup_task(Wl, plate, h)
    else:
        a = p.fast_task(Wl, plate, h, F=F or p.F_FAST, n=n or 1, bundle=bundle or (FAST["bw"], FAST["bh"]))
    return R.get(p.tkey(a))


def tl(b=BF):
    return dict(mode="tailor", Bf=b, tmin=p.T_MIN)


def width_of(r):
    a, b = r["plate_extent"]
    return b - a


def main():
    R = allrecs()
    out = {}
    pr = print
    base_s = get(R, "support", SUP["Wl"], dict(mode="rect", Wp=SUP["Wp"], T=SUP["T"])); base_f = get(R, "fast", FAST["Wl"], dict(mode="rect", Wp=FAST["Wp"], T=FAST["T"]))
    A = gt.load("gap_trade_a.json")
    gs = A[gt.key("support", 0.10, SUP["Wl"], gt.SUP_BUNDLE[0], gt.SUP_BUNDLE[1], 0.005)]; gf = A[gt.key("fast", 0.10, FAST["Wl"], FAST["bw"], FAST["bh"], 0.005)]
    half7 = r7.load("round7_half_2p5.json"); fast7 = r7.load("round7_fast_track.json")

    pr("=" * 150); pr("0. THE MODEL AGAINST THE EARLIER ROUNDS (5 mm cells)"); pr("=" * 150)
    pr(f"  Support track on its 1.2 m x 100 mm plate at the cap: gap trade {gs['cap']['NI']/1e3:.1f} kA for {gs['cap']['F']/1e3:.1f} kN/m, {gs['cap']['tape']:.1f} kg/m, rating {gs['rating']['margin']:.2f};")
    pr(f"    this model {base_s['NI']/1e3:.1f} kA for {base_s['F']/1e3:.1f} kN/m, {base_s['tape']:.1f} kg/m, rating {base_s['rating']['margin']:.2f}; normal operation {base_s['normal']['NI']/1e3:.1f} kA (gap trade {gs['normal']['NI']/1e3:.1f}).")
    pr(f"    The gap trade's ampere-turns moved to 245.5 kN/m by the square law: {gs['cap']['NI']*math.sqrt(base_s['F']/gs['cap']['F'])/1e3:.1f} kA.")
    pr(f"    The pull on the box around the bundles {base_s.get('F_solved', base_s['F'])/1e3:.2f} kN/m, on the box around the plate {base_s['F_plate']/1e3:.2f} kN/m, in the same solve.")
    pr(f"  Compensation track on its 1.0 m x 200 mm plate at 196.5 kN/m: gap trade {gf['NI']/1e3:.1f} kA, {gf['ways']['each stack the better way']['tape']:.1f} kg/m, rating {gf['ways']['each stack the better way']['margin']:.2f};")
    pr(f"    this model {base_f['NI']/1e3:.1f} kA, {base_f['tape']:.1f} kg/m, rating {base_f['rating']['margin']:.2f}; stored energy {base_f['W']/1e3:.1f} kJ/m (gap trade {gf['W']/1e3:.1f}).")
    pr(f"  At 2.5 mm cells (rounds 7 and 8, the book's rows): support {half7['cap']['NI']/1e3:.1f} kA and {half7['cap']['tape']:.1f} kg/m, compensation {fast7['NI_req']/1e3:.1f} kA and {fast7['tape']:.1f} kg/m.")
    pr(f"  The sideways flux in the plate, at its largest: support {base_s['flux_max']:.3f} Wb/m in 100 mm, {base_s['flux_max']/0.1:.2f} T if it were even over the section; compensation {base_f['flux_max']:.3f} Wb/m in 200 mm, {base_f['flux_max']/0.2:.2f} T.")
    m8 = mass_rule(base_s["tape"], [fast4(base_f)])
    m8b = mass_rule(half7["cap"]["tape"], [(4, fast7["tape"], fast7["W"], fast7["V"], fast7["I"], 1)])
    pr(f"  The bearing on the casing by round 8's rule: on this model's 5 mm rows tape {m8['tape']:.0f}, stores {m8['store']:.0f}, switches {m8['switches']:.0f}, leads {m8['leads']:.1f}, total {m8['total']:.0f} kg/m;")
    pr(f"    on round 7's 2.5 mm rows tape {m8b['tape']:.0f}, stores {m8b['store']:.0f}, switches {m8b['switches']:.0f}, leads {m8b['leads']:.1f}, total {m8b['total']:.0f} kg/m (round 8: 1,687; the gap trade's 5 mm row: 1,653).")
    plates0 = 2 * base_s["plate_mass"] + 4 * base_f["plate_mass"]
    pr(f"  The six plates {plates0:.0f} kg/m ({2*base_s['plate_mass']+2*base_f['plate_mass']:.0f} under the cable, {2*base_f['plate_mass']:.0f} over it), with the LIM plates {plates0+M_LIM_PLATES:.0f} kg/m of the {HW_ALLOW:.0f}.")
    out["base"] = dict(sup=dict(NI=base_s["NI"], tape=base_s["tape"], rating=base_s["rating"]["margin"], NI_n=base_s["normal"]["NI"], sat=base_s["sat"], plate=base_s["plate_mass"]),
                       fast=dict(NI=base_f["NI"], tape=base_f["tape"], rating=base_f["rating"]["margin"], W=base_f["W"], V=base_f["V"], I=base_f["I"], plate=base_f["plate_mass"]), m8=m8, m8_25=m8b, plates=plates0,
                       gt=dict(sup_NI=gs["cap"]["NI"], sup_tape=gs["cap"]["tape"], fast_NI=gf["NI"], fast_tape=gf["ways"]["each stack the better way"]["tape"]),
                       r7=dict(sup_NI=half7["cap"]["NI"], sup_tape=half7["cap"]["tape"], fast_NI=fast7["NI_req"], fast_tape=fast7["tape"]), flux=dict(sup=base_s["flux_max"], fast=base_f["flux_max"]),
                       F_box=(base_s.get("F_solved", base_s["F"]), base_s["F_plate"]))

    SK = p.load("round11_plate_skirt.json")
    if SK:
        pr(); pr("=" * 150); pr("1a. HOW FAR THE PLATE REACHES PAST THE OUTER BUNDLES (each track alone at its own design ampere-turns, a rectangular plate; 5 mm cells)"); pr("=" * 150)
        pr(f"  {'track':<8} {'plate mm':>8} {'past the outer bundles m':>24} {'plate width m':>13} {'pull kN/m':>10} {'against the design plate':>24} {'saturated %':>11}")
        for r in SK:
            ref = [q for q in SK if q["kind"] == r["kind"] and q["skirt"] == 0.025 and q["T"] == (SUP["T"] if r["kind"] == "support" else FAST["T"])][0]
            pr(f"  {r['kind']:<8} {r['T']*1e3:8.0f} {r['skirt']:24.3f} {r['width']:13.2f} {r['F']/1e3:10.1f} {(r['F']/ref['F']-1)*100:+22.0f} % {r['sat']*100:11.0f}")
        out["skirt"] = SK

    pr(); pr("=" * 150); pr("1b. A PLATE SIZED FROM ITS FLUX, EACH TRACK ALONE, THE DESIGN'S LOOPS (support 0.4 m loops, compensation 0.6 m loop), 5 mm cells"); pr("=" * 150)
    pr("  The plate's thickness at each 5 mm column is the sideways flux there (from the track on a 300 mm plate that reaches 0.6 m past its outer bundles, at its")
    pr("  requirement) over the fill field, rounded up to 5 mm, 10 mm at least; the plate ends where 10 mm is enough.")
    for kind, base, Wl in (("support", base_s, SUP["Wl"]), ("fast", base_f, FAST["Wl"])):
        nm = "support" if kind == "support" else "compensation"
        pr(f"  {'the ' + nm + ' track':<24} {'plate width m':>13} {'thickest mm':>11} {'mean mm':>8} {'plate kg/m':>10} {'kA per leg':>10} {'worst tape T':>12} {'limit A':>8} {'tape kg/m':>9} {'rating':>7} {'saturated %':>11}")
        rows = []
        Wp, T = (SUP["Wp"], SUP["T"]) if kind == "support" else (FAST["Wp"], FAST["T"])
        pr(f"  {'today: ' + f'{Wp:.1f} m x {T*1e3:.0f} mm':<24} {Wp:13.2f} {T*1e3:11.0f} {T*1e3:8.0f} {base['plate_mass']:10.0f} {base['NI']/1e3:10.1f} {base['B']:12.2f} {base['ic']:8.1f} {base['tape']:9.1f} {base['rating']['margin']:7.2f} {base['sat']*100:11.0f}")
        for b in p.BF_LIST:
            r = get(R, kind, Wl, tl(b))
            if not r:
                continue
            w = width_of(r)
            rows.append(dict(Bf=b, width=w, Tmax=r["plate_Tmax"], mean=r["plate_mass"] / RHO / w, mass=r["plate_mass"], NI=r["NI"], B=r["B"], ic=r["ic"], tape=r["tape"], rating=r["rating"]["margin"], sat=r["sat"]))
            pr(f"  {'fill field ' + f'{b:.1f} T':<24} {w:13.2f} {r['plate_Tmax']*1e3:11.0f} {r['plate_mass']/RHO/w*1e3:8.0f} {r['plate_mass']:10.0f} {r['NI']/1e3:10.1f} {r['B']:12.2f} {r['ic']:8.1f} {r['tape']:9.1f} {r['rating']['margin']:7.2f} {r['sat']*100:11.0f}")
        if rows:
            th = get(R, kind, Wl, tl(rows[0]["Bf"]))["thick"]
            pr(f"    on the 300 mm plate: {th['NI']/1e3:.1f} kA per leg, the sideways flux {th['flux_max']:.3f} Wb/m at its largest.")
            out[f"thick_{kind}"] = dict(NI=th["NI"], flux_max=th["flux_max"])
            ms_ = [r["mass"] for r in rows]; ts = [r["tape"] for r in rows]
            if min(ts) <= base["tape"] <= max(ts):
                o = np.argsort(ts); m_eq = float(np.interp(base["tape"], np.array(ts)[o], np.array(ms_)[o])); b_eq = float(np.interp(base["tape"], np.array(ts)[o], np.array([r["Bf"] for r in rows])[o]))
                pr(f"    The same tape as today's plate ({base['tape']:.1f} kg/m) at a fill field of {b_eq:.2f} T: {m_eq:.0f} kg/m of plate against {base['plate_mass']:.0f} (interpolated between the solved rows).")
                out[f"equal_tape_{kind}"] = dict(Bf=b_eq, mass=m_eq)
            if min(ms_) <= base["plate_mass"] <= max(ms_):
                o = np.argsort(ms_); t_eq = float(np.interp(base["plate_mass"], np.array(ms_)[o], np.array(ts)[o]))
                pr(f"    The same plate mass as today's ({base['plate_mass']:.0f} kg/m): {t_eq:.1f} kg/m of tape against {base['tape']:.1f} (interpolated).")
                out[f"equal_mass_{kind}"] = dict(tape=t_eq)
        out[f"tailor_{kind}"] = rows

    # ---- the rectangular plates of round 5's grid, for the same loops
    import bearing_round5 as r5
    c5 = r5.load_cache()
    pr(); pr("  RECTANGULAR PLATES OF ROUND 5'S GRID for the same loops (bearing_round5.py --part=sweep, runs/bearing_round5_plates.json, 5 mm cells; the support there has its two")
    pr("  center bundles touching, so its rows are compared with its own 1.2 m x 100 mm row):")
    rect = {"support": [], "fast": []}
    for kind, Wl, ref in (("support", SUP["Wl"], (1.2, 0.1)), ("fast", FAST["Wl"], (1.0, 0.2))):
        rows = [v for v in c5.values() if v.get("ok") and v["kind"] == kind and abs(v.get("h", 0) - 0.005) < 1e-9 and abs(v.get("Wl", 0) - Wl) < 1e-9 and v["sat"] < 0.5]
        r0 = [v for v in rows if abs(v["Wp"] - ref[0]) < 1e-9 and abs(v["T"] - ref[1]) < 1e-9]
        if not r0:
            continue
        r0 = r0[0]
        pr(f"  {('support' if kind == 'support' else 'compensation') + ' track: plate':<28} {'plate kg/m':>10} {'against its own design row':>26} {'kA per leg':>10} {'tape kg/m':>9} {'against its own design row':>26} {'rating':>7} {'saturated %':>11}")
        for v in sorted(rows, key=lambda v: v["plate_m"]):
            rt = v.get("F70", 0) / v["F_req"]
            rect[kind].append(dict(Wp=v["Wp"], T=v["T"], mass=v["plate_m"], NI=v["NI_req"], tape=v["tape"], rating=rt, sat=v["sat"], d_mass=v["plate_m"] / r0["plate_m"] - 1, d_tape=v["tape"] / r0["tape"] - 1))
            lab = "%.1f m x %.0f mm" % (v["Wp"], v["T"] * 1e3)
            pr(f"  {lab:<28} {v['plate_m']:10.0f} {(v['plate_m']/r0['plate_m']-1)*100:+24.0f} % {v['NI_req']/1e3:10.1f} {v['tape']:9.1f} {(v['tape']/r0['tape']-1)*100:+24.0f} % {rt:7.2f} {v['sat']*100:11.0f}")
    out["rect5"] = rect

    # ---- which plates fit the cable's allowance
    pr(); pr(f"  WHICH PLATES FIT THE {HW_ALLOW:,.0f} kg/m (the design's loops; two support plates, four compensation plates, the LIM plates {M_LIM_PLATES:,.0f} kg/m; a rating of at least 1.35):")
    opts_s = [("today's 1.2 m x 100 mm", base_s["plate_mass"], base_s["tape"], base_s["rating"]["margin"], None)] + [(f"sized at {r['Bf']:.1f} T", r["mass"], r["tape"], r["rating"], r["Bf"]) for r in out.get("tailor_support", [])]
    opts_f = [("today's 1.0 m x 200 mm", base_f["plate_mass"], base_f["tape"], base_f["rating"]["margin"], None)] + [(f"sized at {r['Bf']:.1f} T", r["mass"], r["tape"], r["rating"], r["Bf"]) for r in out.get("tailor_fast", [])]
    pr(f"  {'support plate':<24} {'compensation plate':<24} {'six plates kg/m':>15} {'hardware kg/m':>13} {'of the allowance':>16} {'tape of the six kg/m':>20} {'against the design':>18} {'bearing on the casing kg/m':>26} {'fits':>5}")
    combos = []
    for a in opts_s:
        for b in opts_f:
            if a[3] < p.gt.MARGIN or b[3] < p.gt.MARGIN:
                continue
            plates = 2 * a[1] + 4 * b[1]; hw = plates + M_LIM_PLATES; tape6 = 2 * a[2] + 4 * b[2]
            fr_ = base_f if b[4] is None else get(R, "fast", FAST["Wl"], tl(b[4]))
            m = mass_rule(a[2], [fast4(fr_)])
            combos.append(dict(sup=a[0], fast=b[0], plates=plates, hw=hw, tape=tape6, bearing=m["total"], fits=hw <= HW_ALLOW))
    for c in sorted(combos, key=lambda c: c["tape"]):
        pr(f"  {c['sup']:<24} {c['fast']:<24} {c['plates']:15.0f} {c['hw']:13.0f} {c['hw']/HW_ALLOW*100:14.0f} % {c['tape']:20.0f} {(c['tape']/m8['tape']-1)*100:+16.0f} % {c['bearing']:26.0f} {'yes' if c['fits'] else 'no':>5}")
    out["combos"] = sorted(combos, key=lambda c: c["tape"])

    # ---- the half face
    F = p.load("round11_plate_face.json") or {}
    C = p.load("round11_plate_coupling.json") or {}
    pr(); pr("=" * 150); pr("1c. THE 7.39 m FACE: THE SUPPORT AND THE LOAD TRACK UNDER ONE PLATE (half face, centerline at x = 0; 5 mm cells, coarse box inside 10 x 8 m)"); pr("=" * 150)
    names = sorted([v["name"] for v in F.values() if v.get("ok")], key=lambda n: (0 if n == "four plates" else 1 if "web" in n else 2 if "everywhere" in n else 3, n))
    pr("  Each track alone at its requirement on the plate; the tape at 50 percent of the worst tape's limit; the rating where the worst tape reaches 70 percent.")
    pr(f"  {'plate':<30} {'lower face plates kg/m':>22} {'thickest mm':>11} | {'support: kA':>11} {'tape kg/m':>9} {'rating':>6} {'sat %':>5} {'normal kA':>9} | {'load: kA':>8} {'tape kg/m':>9} {'rating':>6} {'sat %':>5} {'kJ/m':>6} | "
       f"{'cap with the load full: support, load kN/m':>42} {'worst tape share %: support, load':>33}")
    face_rows = {}
    for nm in names:
        v = F[f"{nm}|pol +1|0.0050"]; s, l = v["sup"], v["load"]; d = v["states"]["support at its cap, load track full"]
        face_rows[nm] = dict(mass=2 * v["plate_mass_half"], Tmax=v["plate_Tmax"], sup=dict(NI=s["NI"], tape=s["tape"], rating=s["rating"]["margin"], sat=s["sat"], NI_n=s["normal"]["NI"], B=s["B"], ic=s["ic"]),
                             load=dict(NI=l["NI"], tape=l["tape"], rating=l["rating"]["margin"], sat=l["sat"], W=l["W"], V=l["V"], I=l["I"], B=l["B"], ic=l["ic"]),
                             both=dict(F_sup=d["F_sup"], F_load=d["F_load"], share_sup=d["sup_share"], share_load=d["load_share"]), thick=v.get("thick"))
        pr(f"  {nm:<30} {2*v['plate_mass_half']:22.0f} {v['plate_Tmax']*1e3:11.0f} | {s['NI']/1e3:11.1f} {s['tape']:9.1f} {s['rating']['margin']:6.2f} {s['sat']*100:5.0f} {s['normal']['NI']/1e3:9.1f} | {l['NI']/1e3:8.1f} {l['tape']:9.1f} "
           f"{l['rating']['margin']:6.2f} {l['sat']*100:5.0f} {l['W']/1e3:6.1f} | {d['F_sup']/1e3:20.1f} {d['F_load']/1e3:21.1f} {d['sup_share']*100:16.1f} {d['load_share']*100:16.1f}")
    out["face"] = face_rows

    pr(); pr("  THE COUPLING: the load track at its full 198 kN/m beside a support track in normal operation (47.5 kN/m), each support loop keeping its flux.")
    pr("  Round 7's method: linear plates, the flux per ampere-turn of each loop and of the load track through each loop, round 7's ampere-turns (453.4 and 156.7 kA).")
    pr("  Direct: the plate saturating, this plate's own ampere-turns, the pull on the support's bundles before and after, the loops' ampere-turns found so that each keeps its flux.")
    pr(f"  {'plate':<30} | {'load flux through loop 1, loop 2, of its own':>44} | {'round 7 method: 10 mm':>21} {'5 mm':>8} | {'direct: support pull':>20} {'kN/m lost':>9} {'loop 1 kA':>9} {'loop 2 kA':>9} {'load pull kN/m':>14} | {'at constant current':>19}")
    cp_rows = {}
    for nm in list(dict.fromkeys(names + [v["name"] for v in C.values()])):
        c10 = C.get(f"{nm}|pol +1|0.0100"); c5 = C.get(f"{nm}|pol +1|0.0050"); v = F.get(f"{nm}|pol +1|0.0050")
        row = dict()
        s1 = f"{c5['share'][0]*100:+6.2f} %, {c5['share'][1]*100:+6.2f} %" if c5 else ""
        s2 = f"{c10['signed']*100:+19.2f} %" if c10 else " " * 21
        s3 = f"{c5['signed']*100:+6.2f} %" if c5 else " " * 8
        if c5:
            row.update(lin10=c10["signed"] if c10 else None, lin5=c5["signed"], share=c5["share"], dNI=c5["dNI"])
        s4 = ""
        if v and v.get("ok"):
            c = v["coupling"]
            row.update(direct=-c["loss"], lost=c["loss"] * c["F0"], x=c["x"], F0=c["F0"], NI0=c["NI0"], const=c["F_sup_const_current"] / c["F0"] - 1, resid=c["resid"], lin_face=c["signed_linear"], F_load=c["F_load"], n_it=len(c["iters"]))
            s4 = f"{-c['loss']*100:+18.2f} % {c['loss']*c['F0']/1e3:9.2f} {c['x'][0]/1e3:9.1f} {c['x'][1]/1e3:9.1f} {c['F_load']/1e3:14.1f} | {(c['F_sup_const_current']/c['F0']-1)*100:+17.2f} %"
        cp_rows[nm] = row
        pr(f"  {nm:<30} | {s1:>44} | {s2} {s3} | {s4}")
    cp7 = r7.load("round7_coupling.json")
    pr(f"  Round 7 (its own box, 9 x 3 m, 10 mm cells, four plates): {-cp7['c_split']*100:+.2f} %.  The load track's polarity is round 7's; turned over, every sign of the linear columns turns.")
    out["coupling"] = cp_rows; out["c_round7"] = cp7["c_split"]
    M = p.load("round11_plate_mesh.json") or {}
    if M:
        pr(); pr("  THE BOX AND THE CELLS (round 7's method, linear plates):")
        for k, v in M.items():
            if k.startswith("box|"):
                pr(f"    {v['lab']:<72} {v['name']:<22} {v['signed']*100:+.2f} %")
        for k, v in M.items():
            if not k.startswith("box|"):
                pr(f"    {v['name']:<22} {v['h']*1e3:4.1f} mm cells, window from {v['win'][0]:.1f} m {'(coarse box only)' if v['h'] == 0.01 else '':<18} {v['signed']*100:+.2f} %")
        out["mesh"] = {k: v["signed"] for k, v in M.items()}
    S = p.load("round11_plate_spacing.json") or {}
    if S:
        pr(); pr("  THE COUPLING AGAINST THE SPACING (a 12 m face; the design's rectangular plates under each track, a web between them; round 7's method, 5 mm cells):")
        ds = sorted({v["d"] for v in S.values()}); ws = sorted({v["web"] for v in S.values()})
        pr(f"    {'bundle to bundle':>16} | " + " ".join(f"{('separate plates' if w == 0 else f'web {w*1e3:.0f} mm'):>16}" for w in ws))
        sp = {}
        for d in ds:
            row = [S.get(f"{d:.3f}|{w:.3f}|+1|0.0050") for w in ws]
            sp[f"{d:.3f}"] = [r["signed"] if r else None for r in row]
            pr(f"    {d*1e3:>13.0f} mm | " + " ".join(f"{(r['signed']*100 if r else float('nan')):>+14.2f} %" for r in row))
        out["spacing"] = dict(d=ds, web=ws, rows=sp)
    FF = p.load("round11_plate_fullface.json") or {}
    if FF:
        pr(); pr("  THE WHOLE FACE IN ONE GRID (10 mm cells, no centerline wall): the change of one support track's pull from each source")
        for k, v in FF.items():
            nm, mir = k.split("|")
            pr(f"    {nm:<24} the far half's currents {'as the near half, mirrored' if mir == '+1' else 'turned over':<27}: near load track {v['load']*100:+.2f} %, far load track {v['loadL']*100:+.3f} %, far support track {v['supL']*100:+.3f} %")
        out["fullface"] = FF

    # ---- item 2
    pr(); pr("=" * 150); pr(f"2. WIDER LOOPS (each track alone on a plate sized from its flux at {BF} T, reaching up to 0.6 m past its outer bundles; 5 mm cells)"); pr("=" * 150)
    wrows = {"support": [], "fast": []}
    for kind, Wls in (("support", p.SUP_WL), ("fast", p.FAST_WL)):
        nm = "support" if kind == "support" else "compensation"
        pr(f"  {'the ' + nm + ' track: loop m':<30} {'bundles span m':>14} {'kA per leg':>10} {'worst tape T':>12} {'limit A':>8} {'tape kg/m':>9} {'rating':>7} {'fills %':>7} {'plate width m':>13} {'thickest mm':>11} {'plate kg/m':>10} {'sat %':>6}" +
           ("" if kind == "support" else f" {'energy kJ/m':>11} {'kV':>6} {'coil A':>7}"))
        for Wl in Wls:
            r = get(R, kind, Wl, tl())
            if not r:
                continue
            wrows[kind].append(r)
            pr(f"  {Wl:<30.2f} {r['span']:14.2f} {r['NI']/1e3:10.1f} {r['B']:12.2f} {r['ic']:8.1f} {r['tape']:9.1f} {r['rating']['margin']:7.2f} {r['fill']*100:7.0f} {width_of(r):13.2f} {r['plate_Tmax']*1e3:11.0f} {r['plate_mass']:10.0f} {r['sat']*100:6.0f}" +
               ("" if kind == "support" else f" {r['W']/1e3:11.1f} {r['V']/1e3:6.1f} {r['I']:7.0f}"))
    out["width"] = {k: [dict(Wl=r["Wl"], span=r["span"], NI=r["NI"], tape=r["tape"], rating=r["rating"]["margin"], fill=r["fill"], width=width_of(r), Tmax=r["plate_Tmax"], mass=r["plate_mass"], sat=r["sat"],
                             W=r.get("W"), V=r.get("V"), I=r.get("I"), B=r["B"], ic=r["ic"]) for r in v] for k, v in wrows.items()}

    def face_table(rule, pick, margin_min=p.gt.MARGIN):
        pr(f"  {rule}")
        pr(f"  {'face m':>6} | {'support loops m':>15} {'tape kg/m':>9} {'rating':>6} | {'load and lift loop m':>20} {'tape kg/m':>9} {'rating':>6} | {'tape of the six kg/m':>20} {'against the design':>18} | {'support plate':>13} {'comp. plate':>11} {'six plates kg/m':>15} {'hardware kg/m':>13} {'of 10,686':>9} | "
           f"{'bearing on the casing kg/m':>26} {'against the design':>18} | {'face left over per half m':>25}")
        res = {}
        for face in FACES:
            best = None
            for rs in wrows["support"]:
                for rf in wrows["fast"]:
                    room = face / 2 - width_of(rs) - width_of(rf)
                    if room < -1e-6 or rs["rating"]["margin"] < margin_min or rf["rating"]["margin"] < margin_min:
                        continue
                    c = pick(rs, rf)
                    if best is None or c < best[0]:
                        best = (c, rs, rf, room)
            if not best:
                pr(f"  {face:6.2f} | no pair of table 2 fits side by side with its plates whole"); continue
            c, rs, rf, room = best
            m = mass_rule(rs["tape"], [fast4(rf)])
            plates = 2 * rs["plate_mass"] + 4 * rf["plate_mass"]; hw = plates + M_LIM_PLATES
            m0 = out["base"]["m8"]
            res[f"{face:.2f}"] = dict(Wls=rs["Wl"], Wlf=rf["Wl"], tape_s=rs["tape"], tape_f=rf["tape"], rating_s=rs["rating"]["margin"], rating_f=rf["rating"]["margin"], tape=m["tape"], plate_s=rs["plate_mass"], plate_f=rf["plate_mass"],
                                      plates=plates, hw=hw, bearing=m["total"], m=m, room=room, width_s=width_of(rs), width_f=width_of(rf), span_s=rs["span"], span_f=rf["span"])
            pr(f"  {face:6.2f} | {rs['Wl']:15.2f} {rs['tape']:9.1f} {rs['rating']['margin']:6.2f} | {rf['Wl']:20.2f} {rf['tape']:9.1f} {rf['rating']['margin']:6.2f} | {m['tape']:20.0f} {(m['tape']/m0['tape']-1)*100:+16.0f} % | {rs['plate_mass']:13.0f} {rf['plate_mass']:11.0f} {plates:15.0f} {hw:13.0f} {hw/HW_ALLOW*100:7.0f} % | "
               f"{m['total']:26.0f} {(m['total']/m0['total']-1)*100:+16.0f} % | {room:25.2f}")
        return res

    if wrows["support"] and wrows["fast"]:
        pr()
        pr("  The fit: on each half of the lower face a support track and a load track side by side, each with its whole plate (so the plates touch at most and the load track's plate ends at")
        pr("  the face's edge); the lift tracks have the upper face to themselves and take the load track's loop; a rating of at least 1.35 (round 5's ruling).")
        out["faces_tape"] = face_table("(a) The loops that need the least tape:", lambda rs, rf: 2 * rs["tape"] + 4 * rf["tape"])
        out["faces_ring"] = face_table("(b) The loops that put the least on the ring, a kilogram on the casing counted as nine on the cable (Chapter 7), the casing's share by round 8's rule:",
                                       lambda rs, rf: 9 * mass_rule(rs["tape"], [fast4(rf)])["total"] + 2 * rs["plate_mass"] + 4 * rf["plate_mass"])
        out["faces_plate"] = face_table("(c) The loops that need the least plate:", lambda rs, rf: 2 * rs["plate_mass"] + 4 * rf["plate_mass"])
    if SIZED in face_rows:
        fr = face_rows[SIZED]; lift = get(R, "fast", FAST["Wl"], tl())
        pr()
        pr("  The 7.39 m face with the design's loops and one plate sized on the face (table 1c), for the same columns:")
        for lab, lf in (("the lift tracks on today's 1.0 m x 200 mm plates", base_f), ("the lift tracks on plates sized at 1.2 T (table 1b)", lift)):
            if not lf:
                continue
            m = mass_rule(fr["sup"]["tape"], [(2, fr["load"]["tape"], fr["load"]["W"], fr["load"]["V"], fr["load"]["I"], 1), fast4(lf, count=2)])
            plates = fr["mass"] + 2 * lf["plate_mass"]; hw = plates + M_LIM_PLATES; m0 = out["base"]["m8"]
            pr(f"    {lab}: tape of the six {m['tape']:.0f} kg/m ({(m['tape']/m0['tape']-1)*100:+.0f} %), plates {fr['mass']:.0f} + {2*lf['plate_mass']:.0f} = {plates:.0f} kg/m, hardware {hw:.0f} kg/m ({hw/HW_ALLOW*100:.0f} % of 10,686), "
               f"the bearing on the casing {m['total']:.0f} kg/m ({(m['total']/m0['total']-1)*100:+.0f} %)")
            out.setdefault("face739", {})[lab] = dict(tape=m["tape"], plates=plates, hw=hw, bearing=m["total"], m=m, lift_tape=lf["tape"], lift_plate=lf["plate_mass"])

    # ---- item 3
    pr(); pr("=" * 150); pr(f"3. A COMPENSATION TRACK THAT DOES THE PAIR'S JOB ALONE: {p.F_ALONE/1e3:.0f} kN/m (each track alone on a plate sized at {BF} T; 5 mm cells)"); pr("=" * 150)
    pr(f"  {'loops':>5} {'loop m':>6} {'bundle mm':>10} {'bundles span m':>14} {'kA per leg':>10} {'worst tape T':>12} {'limit A':>8} {'tape kg/m':>9} {'rating':>7} {'fills %':>7} {'bundle height needed mm':>23} {'plate width m':>13} {'thickest mm':>11} {'plate kg/m':>10} {'energy kJ/m':>11} {'kV per coil':>11} {'coil A':>7}")
    arows = []
    for n, Wl, b in p.ALONE:
        r = get(R, "fast", Wl, tl(), n=n, F=p.F_ALONE, bundle=b)
        if not r:
            continue
        layers = (r["NI"] / 20) / (0.5 * r["ic"]) / 10; V = r["V"] / n
        pk = [r6.pack(b[0], f, layers, V, V / 10) for f in ("h", "v")]; pk = [q for q in pk if q]
        h_need = min(q["h_need"] for q in pk) if pk else float("nan")
        arows.append(dict(n=n, Wl=Wl, bundle=b, span=r["span"], NI=r["NI"], tape=r["tape"], rating=r["rating"]["margin"], fill=r["fill"], h_need=h_need, width=width_of(r), Tmax=r["plate_Tmax"], mass=r["plate_mass"],
                          W=r["W"], V=V, I=r["I"], B=r["B"], ic=r["ic"], sat=r["sat"]))
        pr(f"  {n:5d} {Wl:6.2f} {b[0]*1e3:4.0f} x {b[1]*1e3:3.0f} {r['span']:14.2f} {r['NI']/1e3:10.1f} {r['B']:12.2f} {r['ic']:8.1f} {r['tape']:9.1f} {r['rating']['margin']:7.2f} {r['fill']*100:7.0f} {h_need*1e3:23.0f} {width_of(r):13.2f} {r['plate_Tmax']*1e3:11.0f} "
           f"{r['plate_mass']:10.0f} {r['W']/1e3:11.1f} {V/1e3:11.1f} {r['I']:7.0f}")
    layers = (base_f["NI"] / 20) / (0.5 * base_f["ic"]) / 10
    pk = [r6.pack(FAST["bw"], f, layers, base_f["V"], base_f["V"] / 10) for f in ("h", "v")]; pk = [q for q in pk if q]
    pr(f"  The design's track by the same packing (bearing_round6b.pack): {layers:.1f} tapes per stack, {base_f['V']/1e3:.1f} kV, the 204 mm bundle needs {min(q['h_need'] for q in pk)*1e3:.0f} mm of height and has 50; the tape fills {base_f['fill']*100:.0f} % of it.")
    out["pack_design"] = dict(h_need=min(q["h_need"] for q in pk), fill=base_f["fill"])
    out["alone"] = arows
    if arows and wrows["support"]:
        pr()
        rs = [r for r in wrows["support"] if abs(r["Wl"] - SUP["Wl"]) < 1e-9][0]
        pr(f"  On a 12 m face, beside the support track with the design's loops on its sized plate ({rs['Wl']:.1f} m loops, plate {width_of(rs):.2f} m wide), each with its whole plate:")
        m0 = out["base"]["m8"]
        pr(f"  {'loops':>5} {'loop m':>6} {'bundle mm':>10} | {'lower face left over per half m':>31} {'upper face left over m':>22} | {'tape of the six kg/m':>20} {'six plates kg/m':>15} {'hardware kg/m':>13} {'of 10,686':>9} | {'bearing on the casing kg/m':>26} {'against the design':>18}")
        for a in arows:
            m = mass_rule(rs["tape"], [(4, a["tape"], a["W"], a["V"], a["I"], a["n"])])
            plates = 2 * rs["plate_mass"] + 4 * a["mass"]; hw = plates + M_LIM_PLATES
            a.update(room_bottom=6.0 - width_of(rs) - a["width"], room_top=12.0 - 2 * a["width"], tape6=m["tape"], plates=plates, hw=hw, bearing=m["total"], m=m, sup_Wl=rs["Wl"])
            pr(f"  {a['n']:5d} {a['Wl']:6.2f} {a['bundle'][0]*1e3:4.0f} x {a['bundle'][1]*1e3:3.0f} | {a['room_bottom']:31.2f} {a['room_top']:22.2f} | {m['tape']:20.0f} {plates:15.0f} {hw:13.0f} {hw/HW_ALLOW*100:7.0f} % | {m['total']:26.0f} {(m['total']/m0['total']-1)*100:+16.0f} %")

    CV = p.load("round11_plate_conv.json") or {}
    if CV:
        pr(); pr("=" * 150); pr("THE CELL SIZE AND THE BOX FOR THE SIZED PLATES"); pr("=" * 150)
        out["conv"] = []
        for k, v in CV.items():
            if not v.get("ok"):
                continue
            ref = get(R, v["kind"], v["Wl"], v["plate"])
            big = k.endswith("|outer")
            pr(f"    {v['kind']:<8} loop {v['Wl']:.1f} m, plate sized at {v['plate']['Bf']} T, {v['h']*1e3:.1f} mm cells{', the coarse box inside a larger one' if big else ''}: {v['NI']/1e3:.1f} kA, tape {v['tape']:.1f} kg/m, plate {v['plate_mass']:.0f} kg/m, rating {v['rating']['margin']:.2f}"
               + (f"   (5 mm, the usual box: {ref['NI']/1e3:.1f} kA, {ref['tape']:.1f} kg/m, {ref['plate_mass']:.0f} kg/m, rating {ref['rating']['margin']:.2f})" if ref else ""))
            out["conv"].append(dict(kind=v["kind"], h=v["h"], outer=big, picard=v.get("outer"), NI=v["NI"], tape=v["tape"], mass=v["plate_mass"], rating=v["rating"]["margin"],
                                    ref=dict(NI=ref["NI"], tape=ref["tape"], mass=ref["plate_mass"], rating=ref["rating"]["margin"]) if ref else None))
    p.save("round11_plate_report.json", out)
    return out


if __name__ == "__main__":
    main()
