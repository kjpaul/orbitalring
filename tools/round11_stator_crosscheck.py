#!/usr/bin/env python3
"""
Round 11, item 2: a second model for the loss of the chosen stators' active sides, and the loss split by cause.
Book: Orbital Ring Engineering, Volume III: Section 14.5.
Run: python tools/round11_stator_crosscheck.py     (seconds; needs runs/round11_logs/round11_stator_study.json from
     tools/round11_stator.py; writes runs/round11_logs/round11_stator_crosscheck.txt and .json)

The study's loss is the strip model: each tape alone in the field the others make (Brandt and Indenbom) plus Norris.
The active sides of the chosen stators are tapes on edge, face to face, closer together than they are wide.  Such a stack
screens like a slab (Bean's critical state for a slab; the infinite-stack limit reviewed by Grilli et al., IEEE Trans.
Appl. Supercond. 24, 8200433, 2014).  For a slab of thickness w_t (the tape's width), sheet critical current
K_c = n_t Ic (A/m) and so J_c = K_c / w_t:
  transport alone, sheet current K (field + and - K/2 at the two faces):   Q = mu0 K^3 w_t / (6 K_c)          J per m2 per cycle
  a field H_a along the slab's faces (through the tape faces), H_a below the full penetration field K_c / 2:
                                                                           Q = 4 mu0 H_a^3 w_t / (3 K_c)      J per m2 per cycle
Both are per square meter of the sheet (per meter of track times per meter of coil width).
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

import round11_stator_model as m

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round11_logs")
_lines = []


def out(s=""):
    print(s)
    _lines.append(s)


# ------------------------------------------------------------------------------------------------ the sled's field from racetrack coils
def seg_field(P, A, B):
    """Biot-Savart field (T per A) at points P (n, 3) of a straight filament from A to B"""
    r1 = P - A; r2 = P - B
    n1 = np.linalg.norm(r1, axis=1); n2 = np.linalg.norm(r2, axis=1)
    cr = np.cross(r1, r2); den = n1 * n2 * (n1 * n2 + np.einsum("ij,ij->i", r1, r2))
    return 1e-7 * cr * ((n1 + n2) / np.maximum(den, 1e-30))[:, None]


def racetrack_sled(tau, pack, z, ys, n_poles=6, nx=240, n_f=5):
    """The field of a row of rectangular sled coils, one per pole (tau long, 2 m wide over the conductor centres, alternate
    polarity), each a winding pack `pack` wide laid flat (n_f filaments across it; pack = 0 is one filament), at the plane z
    below the coils' plane, at the lateral positions ys, along one period in x.  Returns Bx, By, Bz arrays (len(ys), nx) per
    ampere-turn of a coil."""
    xs = (np.arange(nx) + 0.5) / nx * 2 * tau
    X, Y = np.meshgrid(xs, ys); P = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z)], axis=1)
    offs = [0.0] if pack <= 0 else list((np.arange(n_f) + 0.5) / n_f * pack - pack / 2)
    B = np.zeros_like(P)
    for j in range(-n_poles, n_poles + 2):
        sgn = 1.0 if j % 2 == 0 else -1.0
        for o in offs:                                     # the filament at offset o from the pack's centre line (outward positive)
            x1 = j * tau + pack / 2 - o; x2 = (j + 1) * tau - pack / 2 + o; a = 1.0 + o
            c = [np.array([x1, -a, 0.0]), np.array([x2, -a, 0.0]), np.array([x2, a, 0.0]), np.array([x1, a, 0.0])]
            for k in range(4):
                B += sgn / len(offs) * seg_field(P, c[k], c[(k + 1) % 4])
    B = B.reshape(len(ys), nx, 3)
    return B[:, :, 0], B[:, :, 1], B[:, :, 2]


def sled_ratios(tau, pack, g=m.G):
    """the racetrack sled's field at the stator as ratios to B0, B0 being the fundamental of its field through the face
    averaged over the 2 m width (the quantity Equation 9.3 uses)"""
    ys = np.linspace(-0.98, 0.98, 29)
    Bx, By, Bz = racetrack_sled(tau, pack, -g, ys)
    fund = lambda A: 2 * np.abs(np.fft.rfft(A, axis=1)[:, 1]) / A.shape[1]
    B0 = fund(Bz.mean(axis=0, keepdims=True))[0]
    yb = np.linspace(0.5, 1.5, 21)                                     # the end band, 1 m wide, centred on the coil's edge
    _, By0, Bz0 = racetrack_sled(tau, pack, -g, yb)
    _, By1, Bz1 = racetrack_sled(tau, pack, -g - 1.0, yb)
    return dict(tau=tau, pack=pack, rz_peak=float(np.abs(Bz).max(axis=1).mean() / B0), rx_fund=float(fund(Bx).mean() / B0), rx_peak=float(np.abs(Bx).max(axis=1).mean() / B0),
                cz_face=float(np.abs(Bz0).max(axis=1).mean() / B0), cz_back=float(np.abs(Bz1).max(axis=1).mean() / B0), cy_face=float(np.abs(By0).max(axis=1).mean() / B0),
                cz_face_fund=float(fund(Bz0).mean() / B0), cz_back_fund=float(fund(Bz1).mean() / B0), cy_face_fund=float(fund(By0).mean() / B0),
                ampere_turns_for_0p1T=float(0.1 / B0))


def main():
    with open(os.path.join(LOGS, "round11_stator_study.json"), encoding="utf-8") as fh:
        study = json.load(fh)
    res = {}
    out("=" * 150)
    out("  ROUND 11, ITEM 2: THE ACTIVE SIDES ON A SECOND MODEL (BEAN SLAB), AND THE LOSS PER CYCLE BY CAUSE, AT THE PEAK CURRENT AND THE FULL SLED FIELD")
    out("=" * 150)
    out(f"  {'launch':<8} {'tau m':>5} {'tapes/m':>8} {'pitch mm':>9} {'tape mm':>8} {'K_belt A/m':>10} {'sled field through the faces mT':>32} | {'strip model J/m: own field only':>32} {'sled field only':>16} {'both':>8} |"
        f" {'slab J/m: transport':>20} {'sled field':>11} {'sum':>8} | {'slab / strip':>12}")
    for name, c in study["chosen"].items():
        r = c["design"]
        st = m.Stator(r["tau"], r["lay"], r["w_t"], r["T"], r["b_e"], r["d"], r["sides"], r["c_end"], r["end_lay"])
        Kb, i_t, n_t, wt = r["K_belt"], r["i_t"], r["n_t"], r["w_t"]
        a_own = st.state(Kb, 0.0, i_t)[0]; a_sled = st.state(1e-9, m.B0N, 1e-9)[0]; a_both = st.state(Kb, m.B0N, i_t)[0]
        q_own = n_t * m.W * (a_own["qp"] + a_own["qn"]); q_sled = n_t * m.W * (a_sled["qp"] + a_sled["qn"]); q_both = n_t * m.W * (a_both["qp"] + a_both["qn"])
        Kc = n_t * a_both["ic_mean"]
        Ha = m.B0N * st.rx / m.MU0
        s_tr = m.MU0 * Kb**3 * wt / (6 * Kc) * m.W
        s_f = (4 * m.MU0 * Ha**3 * wt / (3 * Kc) if Ha < Kc / 2 else m.MU0 * Kc * wt * (Ha - Kc / 3)) * m.W
        res[name] = dict(tau=r["tau"], n_t=n_t, pitch=1 / n_t, w_t=wt, K_belt=Kb, B_face=m.B0N * st.rx, strip_own=q_own, strip_sled=q_sled, strip_both=q_both, slab_transport=s_tr, slab_field=s_f,
                         slab_sum=s_tr + s_f, ratio=(s_tr + s_f) / q_both, lay=r["lay"])
        out(f"  {name:<8} {r['tau']:5.0f} {n_t:8.0f} {1e3/n_t:9.2f} {wt*1e3:8.0f} {Kb:10.0f} {m.B0N*st.rx*1e3:32.1f} | {q_own:32.3f} {q_sled:16.3f} {q_both:8.3f} | {s_tr:20.3f} {s_f:11.3f} {s_tr+s_f:8.3f} | {(s_tr+s_f)/q_both:12.2f}")
    out("  J/m is per meter of track of one side per cycle, the active sides only (2 m of tape per tape).  The slab's field term above full penetration is mu0 K_c w_t (H_a - K_c / 3).")
    out("  The slab model applies where the tapes stand on edge (lay 'edge'): " + ", ".join(f"{k} {v['lay']}" for k, v in res.items()) + ".")
    out()
    out("  BRAKING THE EMPTY SLED FROM 30 km/s ON THE SAME STATOR (the launch of the sled alone run backward: the same thrust and power, the same states)")
    out(f"  {'launch':<8} {'time min':>9} {'distance km':>12} {'heat J':>10} {'against the launch':>19} {'launch and braking J':>21} {'MJ per site':>12}")
    for name, c in study["chosen"].items():
        r = c["design"]
        st = m.Stator(r["tau"], r["lay"], r["w_t"], r["T"], r["b_e"], r["d"], r["sides"], r["c_end"], r["end_lay"])
        b = m.evaluate(st, r["F_pk"], r["P"], r["L_sled"], 0.0, r["v_f"], r["alpha"])
        res[name]["brake"] = dict(t=b["t"], x=b["x"], E_hyst=b["E_hyst"], ratio=b["E_hyst"] / r["E_hyst"], both=b["E_hyst"] + r["E_hyst"], both_site=(b["E_hyst"] + r["E_hyst"]) / m.N_SITES)
        out(f"  {name:<8} {b['t']/60:9.1f} {b['x']/1e3:12,.0f} {b['E_hyst']:10.3e} {b['E_hyst']/r['E_hyst']:19.2f} {b['E_hyst']+r['E_hyst']:21.3e} {(b['E_hyst']+r['E_hyst'])/m.N_SITES/1e6:12.2f}")
    out()
    out("  THE SLED'S FIELD FROM RACETRACK COILS (Biot-Savart over straight filaments): one rectangular coil per pole, 2 m wide over its conductor centres, 100 mm above the")
    out("  stator's face, alternate polarity; its winding pack laid flat and 0, 0.2 or 0.5 m wide.  B0 is the fundamental of the field through the face averaged over the 2 m")
    out("  width, the quantity of Equation 9.3.  Each ratio is the peak over a cycle, averaged over the tape positions, divided by B0.")
    out(f"  {'tau m':>5} {'pack m':>7} {'kA-turns per pole for B0 = 0.1 T':>33} | {'through the face (flat tapes)':>30} {'along the track (on-edge tapes): fundamental':>45} {'peak':>6} {'the study takes':>16} |"
        f" {'end band flat, at the face: peak':>33} {'fundamental':>12} {'set back 1 m: peak':>19} {'fundamental':>12} {'the study takes':>16} {'and':>6}")
    sled = []
    for tau in (4.0, 6.0, 10.0, 20.0, 130.0):
        for pack in (0.0, 0.2, 0.5):
            if pack >= tau / 2:
                continue
            q = sled_ratios(tau, pack); sled.append(q)
            st = m.Stator(tau, "edge", 0.003, 70.0, 1.0, 1.0)
            out(f"  {tau:5.0f} {pack:7.1f} {q['ampere_turns_for_0p1T']/1e3:33.1f} | {q['rz_peak']:30.2f} {q['rx_fund']:45.2f} {q['rx_peak']:6.2f} {st.rx:16.2f} | {q['cz_face']:33.2f} {q['cz_face_fund']:12.2f} {q['cz_back']:19.2f} {q['cz_back_fund']:12.2f} {0.5:16.2f} {st.c_run:6.2f}")
    res["sled_racetrack"] = sled
    out()
    out("  THE BOOK'S LAUNCH WITH THE RACETRACK SLED'S FIELDS AT THE TAPES IN PLACE OF THE STUDY'S (3 mm tape, 70 K; the best lay, band and alpha at each pole pitch)")
    out("  'fundamental': only the fundamental of the racetrack sled's field along the track is taken (a sled wound so that its field is a smooth wave).  'peak': the peak of")
    out("  the field over a cycle is taken, which is what a tape sees under the ends of plain racetrack coils.")
    out(f"  {'sled field':>12} {'pack m':>7} {'tau m':>5} {'lay, active / end runs':<26} {'band':>4} {'alpha':>5} {'A/tape':>7} {'kg/m':>8} {'N':>4} {'heat J':>10} {'% of KE':>8} {'MJ/site':>9} {'against the study at this pitch':>31}")
    b = study["chosen"]["book"]["design"]
    alt = []
    for mode, pack in (("fundamental", 0.2), ("peak", 0.5), ("peak", 0.2), ("peak", 0.0)):
        for tau in (4.0, 6.0, 10.0, 20.0):
            q = next(x for x in sled if x["tau"] == tau and x["pack"] == pack)
            sfx = "_fund" if mode == "fundamental" else ""
            best = None
            for lay, end_lay, d in (("edge", "edge", 0.0), ("edge", "flat", 0.0), ("edge", "flat", 1.0), ("flat", "flat", 0.0), ("flat", "flat", 1.0)):
                for be in (0.5, 1.0):
                    st = m.Stator(tau, lay, 0.003, 70.0, be, d, 2, 0.5, end_lay)
                    st.rx = q["rx_fund"] if mode == "fundamental" else q["rx_peak"]
                    st.rz = 1.0 if mode == "fundamental" else q["rz_peak"]
                    if end_lay == "flat":
                        st.c_run = q["cz_back" + sfx] if d > 0 else q["cz_face" + sfx]
                    else:
                        st.c_run = q["cy_face" + sfx]
                    st.c_end = q["cz_face" + sfx]
                    for al in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
                        r = m.evaluate(st, b["F_pk"], b["P"], b["L_sled"], 5e6, b["v_f"], al)
                        if r["N"] >= 1 and (best is None or r["E_hyst"] < best["E_hyst"]):
                            best = r
            ref = study["search"]["book"]["3mm_70K"]["by_tau"].get(str(tau)) or dict(E_hyst=float("nan"))
            alt.append(dict(mode=mode, pack=pack, tau=tau, lay=best["lay"], end_lay=best["end_lay"], d=best["d"], b_e=best["b_e"], alpha=best["alpha"], i_t=best["i_t"], kg=best["tape_kg_per_m"], N=best["N"],
                            E_hyst=best["E_hyst"], share=best["share"], E_site=best["E_site"], ratio=best["E_hyst"] / ref["E_hyst"]))
            a_ = "on edge" if best["lay"] == "edge" else "flat"; e_ = "on edge" if best["end_lay"] == "edge" else "flat"
            out(f"  {mode:>12} {pack:7.1f} {tau:5.0f} {a_ + ' / ' + e_ + (', set back 1 m' if best['d'] > 0 else ''):<26} {best['b_e']:4.1f} {best['alpha']:5.1f} {best['i_t']:7.1f} {best['tape_kg_per_m']:8.2f} {best['N']:4d} {best['E_hyst']:10.3e} {best['share']*100:8.4f} {best['E_site']/1e6:9.2f} {best['E_hyst']/ref['E_hyst']:31.1f}")
    res["book_on_racetrack_sled"] = alt
    with open(os.path.join(LOGS, "round11_stator_crosscheck.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_lines) + "\n")
    with open(os.path.join(LOGS, "round11_stator_crosscheck.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
