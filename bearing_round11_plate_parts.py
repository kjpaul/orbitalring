#!/usr/bin/env python3
"""
Round 11, item 1b: the parts of bearing_round11_plate.py that put two tracks under one plate (the half face of the cable), the load coil's flux
through the support loops, and the mesh checks.
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run: through bearing_round11_plate.py --part=face|coupling|spacing|mesh (see its header).  Writes runs/round11_logs/round11_plate_*.json.
"""
import math
import time

import numpy as np

import bearing_electrical as be
import bearing_design as bd
import bearing_round7 as r7
import bearing_round11_plate as p
from bearing_round11_plate import Plate, Trk, Face, GAP, HALF_FACE, SUP, FAST, SUP_BUNDLE, F_FAST, F_FULL, W_HALF

X_S = HALF_FACE - FAST["Wp"] - 0.2 - SUP["Wp"] / 2          # 1.895 m: the support track's center on the 7.39 m face (bearing_round7_parts.X_S)
X_L = HALF_FACE - FAST["Wp"] / 2                            # 3.195 m: the load track's center
BOX = (0.0, 5.7, -1.5, 1.5)
WIN = (0.0, 4.3, -0.45, 0.65)
OUTER = (0.0, 10.0, -4.0, 4.0)                              # a 20 mm box around the coarse box: the coupling needs it (see part_mesh)


def tracks_739(pol=1, xs=X_S, xl=X_L):
    return [Trk("sup", "support", xs, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1]), Trk("load", "fast", xl, FAST["Wl"], FAST["bw"], FAST["bh"], pol=pol)]


def rec(key):
    d = {}
    for f in ("round11_plate_tracks.json", "round11_plate_width.json", "round11_plate_width2.json", "round11_plate_width3.json", "round11_plate_alone.json", "round11_plate_design.json"):
        d.update(p.load(f) or {})
    r = d.get(key)
    if not r or not r.get("ok"):
        raise SystemExit(f"missing solve: {key}")
    return r


def base_recs(h):
    s = rec(p.tkey(p.tasks_base(h)[0])); f = rec(p.tkey(p.tasks_base(h)[1]))
    return s, f


def tail_recs(h, bfs, bff, Wls=None, Wlf=None):
    s = rec(p.tkey(p.sup_task(Wls or SUP["Wl"], dict(mode="tailor", Bf=bfs, tmin=p.T_MIN), h)))
    f = rec(p.tkey(p.fast_task(Wlf or FAST["Wl"], dict(mode="tailor", Bf=bff, tmin=p.T_MIN), h)))
    return s, f


def plates_739(h=0.005):
    """The plates of the half face, 0 to 3.695 m: the four plates of the design (here two of them) and the same with a web between them."""
    four = [Plate.rect(X_S - SUP["Wp"] / 2, X_S + SUP["Wp"] / 2, SUP["T"]), Plate.rect(X_L - FAST["Wp"] / 2, X_L + FAST["Wp"] / 2, FAST["T"])]
    out = {"four plates": Plate.combine(four, 0.0, HALF_FACE)}
    for w in (0.010, 0.025, 0.050, 0.100):
        out[f"one plate, web {w*1e3:.0f} mm"] = Plate.combine(four, 0.0, HALF_FACE, floor=w)
    for T in (0.100, 0.200):
        out[f"one plate, {T*1e3:.0f} mm everywhere"] = Plate.rect(0.0, HALF_FACE, T)
    return out


def sized_plates():
    """The plates sized on the face itself by job_face, read back from its log (for the coupling table)."""
    F = p.load("round11_plate_face.json") or {}
    return {v["name"]: Plate.from_json(v["plate"]["json"]) for v in F.values() if v.get("ok") and v["name"].startswith("sized")}


def loop_NI(t, x):
    a = np.array(t.loop(0, x[0])) + np.array(t.loop(1, x[1]))
    return a.tolist()


def job_face(a):
    """One plate of the 7.39 m half face with the support and the load track under it.
    plate: {mode 'given', json} or {mode 'sized', Bf, tmin}.  Sized: both tracks under a 300 mm plate over the whole half face, each at its
    requirement alone; the plate's thickness at each 5 mm column is then the sum of the two sideways fluxes there over Bf, rounded up to 5 mm,
    tmin at least, over the whole half face (one plate).
    Then on the plate: each track alone at its requirement (the ampere-turns, the worst tape, the tape, the rating), the support in normal
    operation, the support at its cap with the load track at its full 198 kN/m, and the coupling by the direct way: the support's two loops keep
    the flux they had in normal operation while the load track comes to its full pull (Newton steps on the two loops' ampere-turns with the
    linear plate's inductances), and the pull on the support's own bundles before and after."""
    t0 = time.time()
    theva, theva77 = be.load_theva()
    name, h, pol, g0 = a["name"], a["h"], a["pol"], a["guess"]
    out = dict(a, key=a["key"])
    try:
        tr = tracks_739(pol); sup, load = tr
        F_cap = p.gt.F_cap()
        if a["plate"]["mode"] == "sized":
            thick = Face(tr, Plate.rect(0.0, HALF_FACE, 0.30), BOX, WIN, h=h, outer=OUTER)
            env = None; th = {}
            for t, F_req, gs in ((sup, F_cap, g0["cap"]), (load, F_FAST, g0["fast"])):
                def solveF(NI, t=t):
                    g, s = thick.solve({t.name: t.design(NI)})
                    return s["tracks"][t.name]["F"], (g, s)
                NI0, NI0s, F0s, (g, s) = p.find_NI(solveF, F_req, gs)
                x, phi = p.plate_flux(g, thick.plate); phi = np.abs(phi) * (NI0 / NI0s)
                env = phi if env is None else env + phi
                th[t.name] = dict(NI=NI0, flux_max=float(phi.max()))
            hp = p.HP; n = int(round(HALF_FACE / hp)); c = (np.arange(n) + 0.5) * hp
            need = np.interp(c, x, env) / a["plate"]["Bf"]
            plate = Plate(0.0, np.minimum(np.maximum(np.ceil(need / hp - 1e-6) * hp, a["plate"]["tmin"]), 0.40))
            out["thick"] = dict(th, solves=thick.solves)
            g0 = dict(cap=th["sup"]["NI"], fast=th["load"]["NI"])
        else:
            plate = Plate.from_json(a["plate"]["json"])
        face = Face(tr, plate, BOX, WIN, h=h, outer=OUTER)
        rs, _ = p.evaluate(face, sup, F_cap, g0["cap"], theva, theva77, normal=True)
        rl, _ = p.evaluate(face, load, F_FAST, g0["fast"], theva, theva77)
        out["sup"] = rs; out["load"] = rl
        NI_full = rl["NI"] * math.sqrt(F_FULL / F_FAST); NI_n = rs["normal"]["NI"]
        out["states"] = {}
        for lab, ns, nl in (("support normal", NI_n, 0.0), ("support normal, load track full", NI_n, NI_full), ("support at its cap, load track full", rs["NI"], NI_full)):
            ni = {"sup": sup.design(ns), "load": load.design(nl)}
            g, s = face.solve(ni)
            d = dict(NI_sup=ns, NI_load=nl, F_sup=s["tracks"]["sup"]["F"], F_load=s["tracks"]["load"]["F"], F_plate=s["F_plate"], sat_sup=s["tracks"]["sup"]["sat"], sat_load=s["tracks"]["load"]["sat"],
                     phi=s["tracks"]["sup"]["phi"], outer=s["outer"])
            ic = p.worst(g, sup, ni["sup"], theva, theva77); d["sup_ic"] = ic[0]; d["sup_share"] = (ns / rs["tapes_leg"]) / ic[0]
            if nl:
                ic = p.worst(g, load, ni["load"], theva, theva77); d["load_ic"] = ic[0]; d["load_share"] = (nl / rl["tapes_leg"]) / ic[0]
            out["states"][lab] = d
        s0 = out["states"]["support normal"]; phi0 = np.array(s0["phi"]); F0 = s0["F_sup"]
        L = np.zeros((2, 2))
        for k in range(2):
            g, s = face.solve({"sup": sup.loop(k, 1e3)}, linear=True); L[:, k] = np.array(s["tracks"]["sup"]["phi"]) / 1e3
        g, s = face.solve({"load": load.design(1e3)}, linear=True); M = np.array(s["tracks"]["sup"]["phi"]) / 1e3
        x = np.array([NI_n, NI_n], dtype=float); it = []
        for n in range(5):
            g, s = face.solve({"sup": loop_NI(sup, x), "load": load.design(NI_full)})
            phi = np.array(s["tracks"]["sup"]["phi"]); r = phi - phi0
            it.append(dict(x=x.tolist(), phi=phi.tolist(), F_sup=s["tracks"]["sup"]["F"], F_load=s["tracks"]["load"]["F"], F_plate=s["F_plate"]))
            if np.abs(r).max() / np.abs(phi0).max() < 2e-3:
                break
            x = x - np.linalg.solve(L, r)
        dNI = -np.linalg.solve(L, M) * NI_full
        out["coupling"] = dict(L=L.tolist(), M=M.tolist(), NI_full=NI_full, NI0=NI_n, F0=F0, x=it[-1]["x"], F_sup=it[-1]["F_sup"], F_load=it[-1]["F_load"],
                               F_plate=it[-1]["F_plate"], loss=1 - it[-1]["F_sup"] / F0, iters=it, resid=float(np.abs(np.array(it[-1]["phi"]) - phi0).max() / np.abs(phi0).max()),
                               c_linear=float(abs(dNI.sum()) / NI_n), signed_linear=float(dNI.sum() / NI_n), dNI_linear=dNI.tolist(),
                               F_sup_const_current=out["states"]["support normal, load track full"]["F_sup"])
        out["plate"] = dict(a["plate"], json=plate.to_json())
        out["plate_mass_half"] = plate.mass(); out["plate_Tmax"] = float(plate.T.max()); out["solves"] = face.solves; out["ok"] = True
    except Exception as e:
        import traceback
        out["ok"] = False; out["error"] = repr(e) + traceback.format_exc()
    out["seconds"] = time.time() - t0
    return out


def part_face(h=0.005, jobs=2, bf="1.2,1.8", only=None, **kw):
    sb, fb = base_recs(h)
    guess = dict(cap=sb["NI"], fast=fb["NI"])
    tasks = []
    for name, pl in plates_739(h).items():
        if not only and ("50 mm" in name or "web 100" in name or "200 mm everywhere" in name):        # left out to save time; --only= runs any one of them
            continue
        tasks.append(dict(key=f"{name}|pol +1|{h:.4f}", name=name, h=h, pol=1, guess=guess, plate=dict(mode="given", json=pl.to_json())))
    for b in [float(x) for x in bf.split(",") if x]:
        name = f"sized on the face at {b:.1f} T"
        tasks.append(dict(key=f"{name}|pol +1|{h:.4f}", name=name, h=h, pol=1, guess=guess, plate=dict(mode="sized", Bf=b, tmin=p.T_MIN)))
    if only:
        tasks = [t for t in tasks if only in t["name"]]
    print(f"  The half face of 7.39 m: the support track at {X_S:.3f} m and the load track at {X_L:.3f} m from the cable's centerline; {len(tasks)} cases", flush=True)
    import multiprocessing as mp
    done = p.load("round11_plate_face.json") or {}
    todo = [t for t in tasks if not (done.get(t["key"]) or {}).get("ok")]
    if todo:
        with mp.Pool(min(jobs, len(todo))) as pool:
            for res in pool.imap_unordered(job_face, todo):
                done = p.load("round11_plate_face.json") or {}; done[res["key"]] = res; p.save("round11_plate_face.json", done)
                if res["ok"]:
                    c = res["coupling"]
                    print(f"    {res['key']}: plate {2*res['plate_mass_half']:.0f} kg/m for the face; support {res['sup']['NI']/1e3:.1f} kA, {res['sup']['tape']:.1f} kg/m, rating {res['sup']['rating']['margin']:.2f}; "
                          f"load {res['load']['NI']/1e3:.1f} kA, {res['load']['tape']:.1f} kg/m, rating {res['load']['rating']['margin']:.2f}; "
                          f"coupling {-c['loss']*100:+.1f} % direct, {c['signed_linear']*100:+.1f} % by the linear fluxes ({res['seconds']:.0f} s)", flush=True)
                else:
                    print(f"    {res['key']}: FAILED {res['error'][:600]}", flush=True)
    return done


def job_skirt(a):
    kind, skirt, T, NI, h = a
    t = Trk("t", kind, 0.0, SUP["Wl"] if kind == "support" else FAST["Wl"], *(SUP_BUNDLE if kind == "support" else (FAST["bw"], FAST["bh"])))
    plate = Plate.rect(-(t.span / 2 + skirt), t.span / 2 + skirt, T)
    face = p.single_face(t, plate, h, outer=True)
    g, s = face.solve({"t": t.design(NI)})
    return dict(kind=kind, skirt=skirt, T=T, NI=NI, F=s["tracks"]["t"]["F"], sat=s["tracks"]["t"]["sat"], width=t.span + 2 * skirt, mass=plate.mass())


def part_skirt(h=0.005, jobs=2, **kw):
    """How far the plate should reach past the outer bundles: each track at its own ampere-turns on rectangular plates of its own thickness and of
    300 mm, the plate's edge 25 mm to 1.5 m past the outer bundles."""
    sb, fb = base_recs(h)
    tasks = []
    for kind, NI, T0 in (("support", sb["NI"], SUP["T"]), ("fast", fb["NI"], FAST["T"])):
        for T in (T0, 0.30):
            for sk in (0.025, 0.2, 0.6, 1.2):
                tasks.append((kind, sk, T, NI, h))
    import multiprocessing as mp
    with mp.Pool(jobs) as pool:
        rows = pool.map(job_skirt, tasks)
    p.save("round11_plate_skirt.json", rows)
    print(f"  {'track':<8} {'plate mm':>8} {'past the outer bundles m':>24} {'plate width m':>13} {'pull kN/m':>10} {'saturated %':>11}")
    for r in rows:
        print(f"  {r['kind']:<8} {r['T']*1e3:8.0f} {r['skirt']:24.3f} {r['width']:13.2f} {r['F']/1e3:10.1f} {r['sat']*100:11.0f}")
    return rows


# =============================================================================
# the load coil's flux through the support loops, round 7's method (linear plates), for any plate and spacing
# =============================================================================

def coupling_linear(tr, plate, box, win, h, coarse_only=False, NI_l=None, NI0=None, outer=None):
    """bearing_round7_parts.part_coupling's method on this model: the two support loops' own and mutual flux per ampere-turn, the load track's flux
    through each, and the support's ampere-turns that keep each loop's flux when the load track comes to its full pull."""
    cp = r7.load("round7_coupling.json")
    NI_l = NI_l or cp["NI_l"]; NI0 = NI0 or cp["NI0"]
    sup, load = tr[0], tr[1]
    face = Face(tr, plate, box, win, h=h, outer=outer)
    L = np.zeros((2, 2))
    for k in range(2):
        g, s = face.solve({"sup": sup.loop(k, 1e3)}, linear=True, coarse_only=coarse_only); L[:, k] = np.array(s["tracks"]["sup"]["phi"]) / 1e3
    g, s = face.solve({"load": load.design(1e3)}, linear=True, coarse_only=coarse_only); M = np.array(s["tracks"]["sup"]["phi"]) / 1e3
    Ll = s["tracks"]["load"]["phi"][0] / 1e3
    dNI = -np.linalg.solve(L, M) * NI_l
    return dict(L=L.tolist(), M=M.tolist(), L_load=Ll, dNI=dNI.tolist(), c=float(abs(dNI.sum()) / NI0), signed=float(dNI.sum() / NI0), NI_l=NI_l, NI0=NI0,
                share=[float(M[0] / Ll), float(M[1] / Ll)])


def part_coupling(h=0.005, jobs=1, **kw):
    """The 7.39 m face: every plate at 10 mm (the coarse box alone, round 7's cell) and at h."""
    plates = dict(plates_739(h)); plates.update(sized_plates())
    out = p.load("round11_plate_coupling.json") or {}
    cp = r7.load("round7_coupling.json")
    print(f"  Round 7 (bearing_round7.py --part=coupling, its own grid, 10 mm cells, two plates 0.2 m apart): {cp['c_split']*100:.2f} percent.")
    print(f"  The load track at its full pull {cp['NI_l']/1e3:.1f} kA per leg, the support in normal operation {cp['NI0']/1e3:.1f} kA per leg (round 7's values, 2.5 mm solves).")
    print(f"  {'plate':<34} {'cells':>6} | {'load flux through loop 1, 2 (share of its own)':>46} | {'loop 1 kA':>9} {'loop 2 kA':>9} | {'support pull':>12}")
    for name, pl in plates.items():
        for pol in (1,):                                              # the plates are linear here, so turning the load track over only turns the sign
            for hh, co in ((0.01, True), (h, False)):
                k = f"{name}|pol {pol:+d}|{hh:.4f}"
                if k not in out:
                    out[k] = dict(coupling_linear(tracks_739(pol), pl, BOX, WIN, hh if not co else h, coarse_only=co, outer=OUTER), name=name, pol=pol, h=hh, mass_face=2 * pl.mass())
                    p.save("round11_plate_coupling.json", out)
                c = out[k]
                print(f"  {name:<34} {hh*1e3:>4.0f}mm | {c['share'][0]*100:>21.2f} % {c['share'][1]*100:>21.2f} % | {c['dNI'][0]/1e3:>+9.2f} {c['dNI'][1]/1e3:>+9.2f} | {c['signed']*100:>+10.2f} %", flush=True)
    return out


def part_spacing(h=0.005, jobs=1, **kw):
    """A 12 m face (half 6 m): the load track at the edge as today, the support track moved in so that the clear distance between the two tracks'
    bundles is d; rectangular plates under each track as today's, with a web between them of 0 (separate plates), 10, 25 and 100 mm."""
    out = p.load("round11_plate_spacing.json") or {}
    half = 6.0; xl = half - FAST["Wp"] / 2
    box = (0.0, 8.0, -1.5, 1.5); win = (0.0, 6.6, -0.45, 0.65)
    d0 = (X_L - FAST["Wl"] / 2 - FAST["bw"] / 2) - (X_S + SUP["Wl"] + (SUP_BUNDLE[0] + p.GC) / 2 + SUP_BUNDLE[0] / 2)
    print(f"  The clear distance between the support's outer bundle and the load track's inner bundle is {d0*1e3:.0f} mm in the design (plates 0.2 m apart).")
    print(f"  {'bundle to bundle':>16} {'web':>8} {'load polarity':>13} | {'support pull':>12}")
    for d in (d0, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        xs = xl - FAST["Wl"] / 2 - FAST["bw"] / 2 - d - SUP_BUNDLE[0] / 2 - SUP["Wl"] - (SUP_BUNDLE[0] + p.GC) / 2
        two = [Plate.rect(xs - SUP["Wp"] / 2, xs + SUP["Wp"] / 2, SUP["T"]), Plate.rect(xl - FAST["Wp"] / 2, xl + FAST["Wp"] / 2, FAST["T"])]
        for w in (0.0, 0.010, 0.025, 0.100):
            for pol in (1,):
                k = f"{d:.3f}|{w:.3f}|{pol:+d}|{h:.4f}"
                if k not in out:
                    out[k] = dict(coupling_linear(tracks_739(pol, xs, xl), Plate.combine(two, 0.0, half, floor=w), box, win, h, outer=(0.0, 12.0, -4.0, 4.0)), d=d, web=w, pol=pol, h=h, xs=xs)
                    p.save("round11_plate_spacing.json", out)
                print(f"  {d*1e3:>13.0f} mm {w*1e3:>5.0f} mm {pol:>+13d} | {out[k]['signed']*100:>+10.2f} %", flush=True)
    return out


def part_fullface(h=0.01, jobs=1, **kw):
    """The whole 7.39 m face in one 10 mm grid (no symmetry plane), one plate with a 25 mm web and the four plates: the flux of the far load track
    through the near support loops, and whether the two halves' polarities matter."""
    out = {}
    for name, floor in (("four plates", 0.0), ("one plate, web 25 mm", 0.025), ("one plate, web 100 mm", 0.100)):
        for mirror in (1, -1):
            tr = tracks_739(1) + [Trk("supL", "support", -X_S, SUP["Wl"], SUP_BUNDLE[0], SUP_BUNDLE[1], pol=mirror), Trk("loadL", "fast", -X_L, FAST["Wl"], FAST["bw"], FAST["bh"], pol=mirror)]
            four = [Plate.rect(X_S - 0.6, X_S + 0.6, 0.1), Plate.rect(X_L - 0.5, X_L + 0.5, 0.2), Plate.rect(-X_S - 0.6, -X_S + 0.6, 0.1), Plate.rect(-X_L - 0.5, -X_L + 0.5, 0.2)]
            pl = Plate.combine(four, -HALF_FACE, HALF_FACE, floor=floor)
            face = Face(tr, pl, (-5.7, 5.7, -1.5, 1.5), (-4.3, 4.3, -0.45, 0.65), h=h, outer=(-10.0, 10.0, -4.0, 4.0))
            cp = r7.load("round7_coupling.json")
            L = np.zeros((2, 2)); sup = tr[0]
            for k in range(2):
                g, s = face.solve({"sup": sup.loop(k, 1e3)}, linear=True, coarse_only=True); L[:, k] = np.array(s["tracks"]["sup"]["phi"]) / 1e3
            res = {}
            for src in ("load", "loadL", "supL"):
                t = [x for x in tr if x.name == src][0]
                g, s = face.solve({src: t.design(1e3)}, linear=True, coarse_only=True); M = np.array(s["tracks"]["sup"]["phi"]) / 1e3
                NI = cp["NI_l"] if src != "supL" else cp["NI0"]
                dNI = -np.linalg.solve(L, M) * NI
                res[src] = float(dNI.sum() / cp["NI0"])
            out[f"{name}|{mirror:+d}"] = res
            print(f"  {name:<24} far half's currents {'as the near half' if mirror == 1 else 'turned over':<17}: near load track {res['load']*100:+.2f} %, far load track {res['loadL']*100:+.3f} %, "
                  f"far support track at its normal current {res['supL']*100:+.3f} % of the support's pull", flush=True)
    p.save("round11_plate_fullface.json", out)
    return out


def part_mesh(h=0.005, jobs=1, **kw):
    """The box and the cell size for the coupling (linear plates): round 7's box, the coarse box alone, the coarse box inside a larger one, and
    cells of 10, 5 and 2.5 mm (the window narrowed for 2.5 mm)."""
    out = p.load("round11_plate_mesh.json") or {}
    plates = plates_739(h)
    xs7, xl7 = 0.0, SUP["Wp"] / 2 + 0.2 + FAST["Wp"] / 2
    two7 = [Plate.rect(-0.6, 0.6, SUP["T"]), Plate.rect(xl7 - 0.5, xl7 + 0.5, FAST["T"])]
    print("  The box (10 mm cells, the load track's polarity as round 7's):")
    for lab, box, outer in (("round 7's box, 9 x 3 m, the support alone at its center, no centerline", (-3.0, 6.0, -1.0, 2.0), None),
                            ("the same in a box 15 x 8 m", (-3.0, 6.0, -1.0, 2.0), (-6.0, 9.0, -4.0, 4.0))):
        for name, fl in (("four plates", 0.0), ("one plate, web 25 mm", 0.025)):
            k = f"box|{lab}|{name}"
            if k not in out:
                out[k] = dict(coupling_linear(tracks_739(1, xs7, xl7), Plate.combine(two7, -0.6, xl7 + 0.5, floor=fl), box, (-1.5, 2.5, -0.45, 0.65), 0.005, coarse_only=True, outer=outer), name=name, lab=lab)
                p.save("round11_plate_mesh.json", out)
            print(f"    {lab:<72} {name:<22} {out[k]['signed']*100:+.2f} %", flush=True)
    for lab, outer in (("the half face, coarse box 5.7 x 3 m alone", None), ("the half face, coarse box inside 10 x 8 m", OUTER), ("the half face, coarse box inside 14 x 12 m", (0.0, 14.0, -6.0, 6.0))):
        for name in ("four plates", "one plate, web 25 mm"):
            k = f"box|{lab}|{name}"
            if k not in out:
                out[k] = dict(coupling_linear(tracks_739(1), plates[name], BOX, WIN, 0.005, coarse_only=True, outer=outer), name=name, lab=lab)
                p.save("round11_plate_mesh.json", out)
            print(f"    {lab:<72} {name:<22} {out[k]['signed']*100:+.2f} %", flush=True)
    print("  The cells (the half face, the coarse box inside 10 x 8 m):")
    for name in ("four plates", "one plate, web 25 mm"):
        for hh, co, win in ((0.01, True, WIN), (0.005, False, WIN), (0.005, False, (1.0, 4.0, -0.35, 0.45)), (0.0025, False, (1.0, 4.0, -0.35, 0.45))):
            k = f"{name}|{hh:.4f}|{win[0]:.1f}"
            if k not in out:
                out[k] = dict(coupling_linear(tracks_739(1), plates[name], BOX, win, hh if not co else 0.005, coarse_only=co, outer=OUTER), name=name, h=hh, win=win)
                p.save("round11_plate_mesh.json", out)
            print(f"    {name:<24} {hh*1e3:4.1f} mm cells, window from {win[0]:.1f} m: the support's pull {out[k]['signed']*100:+.2f} %", flush=True)
    return out


def part_report(**kw):
    import bearing_round11_plate_report as rp
    rp.main()
