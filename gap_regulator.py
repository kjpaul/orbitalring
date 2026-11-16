#!/usr/bin/env python3
"""
The gap regulator with the roll carried (round 10, 2026-11-17).
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Brief: Vol_III/Claude_Code_instructions_round10_2026-11-17.md, item 1.

One cross-section of the ring, per meter of track: the casing and the cable, each free in heave and roll (the cable with no torsional stiffness, the
casing a 10 m square frame, radius of gyration 5.77 m [UNVERIFIED, as bearing_design F.3]; nothing along the track, so no torsion of the casing
between a loaded stretch and the bare casing beside it).  The cable's orbit carries the casing's weight.

  supports   each part (a whole track or one loop of the split track) pulls at its own place with its frozen pull times the constant-flux factor
             on its own gap (sled_passage_pulsed.SUPPORT_OF_GAP), less the coupling to the load track beside it.
  fast       load tracks (below, they pull the casing up) and lift tracks (above, they pull it down), each at its own place, each only able to pull,
             each with its own rating and its own current that moves no faster than full pull in one rise time (the supply's voltage).
  sled       its load q(t) on the centerline, from the point model's lane.
  plan       the planned pull of the launch (the point model's command, trimmed there so its delivered impulse is the planned one), shared among the
             tracks of the face it uses.
  regulator  sled_passage_pulsed.regulator() on each side of the cable separately: the side's gap error is its measured gap against the path the plan
             predicts, read every sample with the sensors' delay; the pull it asks for is added to that side's share of the plan and goes to the load
             track(s) of the side if it is up and the lift track(s) if it is down.

The point model (sled_passage_pulsed.point_pulsed with reg=) is the same law with the two sides acting as one.
"""
import math

import numpy as np

import ring_altitude as ring
import sled_passage_pulsed as pz

R_GYR = 10.0 / math.sqrt(3.0)


def track(kind, x, F_full, cap, share, lost_t=None, name=""):
    """A fast track or a part of one.  kind: 'load' (pulls the casing up) or 'lift' (down).  share: its share of the planned pull of its face.
    lost_t: the time it is lost (None: never; -inf: from the start)."""
    return dict(kind=kind, x=x, F_full=F_full, cap=cap, share=share, lost_t=lost_t, name=name or f"{kind} {x:+.2f}")


def support(x, S0, lost_t=None, name=""):
    return dict(x=x, S0=S0, lost_t=lost_t, name=name or f"support {x:+.2f}")


def layout(X_C, X_S, F_full_lift, F_full_load, cap_lift, cap_load, W, halves=False, split_sup=True, half_off=0.25, sup_off=0.225):
    """The six tracks.  halves: each fast track as two independent halves (half the pull and half the rating each, half_off either side of the
    track's center [the halves' own solves are not made: an assumption]).  split_sup: each support track as its two loops."""
    fast = []
    for sgn in (+1, -1):
        for kind, Ff, cap in (("lift", F_full_lift, cap_lift), ("load", F_full_load, cap_load)):
            if halves:
                for o in (-half_off, half_off):
                    fast.append(track(kind, sgn * X_C + o, Ff / 2, cap / 2, 0.25))
            else:
                fast.append(track(kind, sgn * X_C, Ff, cap, 0.5))
    sup = []
    for sgn in (+1, -1):
        if split_sup:
            for o in (-sup_off, sup_off):
                sup.append(support(sgn * X_S + o, W / 4))
        else:
            sup.append(support(sgn * X_S, W / 2))
    return fast, sup


def section_run(t, q, q_pred, F_cmd, fast, sup, reg, tau=10e-3, W=None, MU=None, MU_C=None, X_side=None, side_err=(0.0, 0.0), couple=0.0, plan="kept",
                edge=None, keep_hist=False):
    """Integrate one cross-section.  t, q, q_pred, F_cmd: from a point lane (hist of point_pulsed with reg=).  reg: a regulator() dict made for one
    side (mass = MU / 2).  side_err: the gap sensors' error on the + side and the - side.  couple: the load track's coupling to the support of its
    side (dF/F of that support at the load track's full pull).  plan, when a track of the planned face is lost: 'kept' (each track left keeps its
    own share), 'moved' (the lost track's share goes to the tracks left on its face, up to their ratings), 'dropped' (no planned pull at all),
    'balanced' (each side gives the same planned pull, the smaller of what the two sides can give)."""
    W = pz.W if W is None else W; MU = pz.MU if MU is None else MU; MU_C = pz.MU_C if MU_C is None else MU_C
    edge = ring.CABLE_SIDE / 2 if edge is None else edge
    X_side = max(abs(f["x"]) for f in fast) if X_side is None else X_side
    I_cas = MU * R_GYR**2; I_cab = MU_C * (2 * ring.CABLE_SIDE**2) / 12.0
    dt = float(t[1] - t[0]); nt = len(t); sog = pz.SUPPORT_OF_GAP; G0 = pz.GAP0
    nd = int(round(reg["delay"] / dt)); ks = max(int(round(reg["sample"] / dt)), 1); kw = max(int(round(reg["window"] / dt)), 1)
    om = reg["om"]; m_side = reg["mass"]; ki = reg["ki"]; omr = om / 10.0
    lag = 1.0 - math.exp(-dt / tau)
    nf = len(fast)
    fx = np.array([f["x"] for f in fast]); fsgn = np.array([1.0 if f["kind"] == "load" else -1.0 for f in fast])
    ffull = np.array([f["F_full"] for f in fast]); fcap = np.array([f["cap"] for f in fast]); fshare = np.array([f["share"] for f in fast])
    flost = np.array([(-np.inf if f["lost_t"] == "start" else (np.inf if f["lost_t"] is None else f["lost_t"])) for f in fast])
    fside = np.where(fx > 0, 0, 1)                       # 0: the + side, 1: the - side
    rate = np.sqrt(ffull) / reg["rise"]
    s = np.zeros(nf)                                     # the root of each track's pull
    sx = np.array([p["x"] for p in sup]); s0 = np.array([p["S0"] for p in sup]); slost = np.array([(np.inf if p["lost_t"] is None else p["lost_t"]) for p in sup])
    sside = np.where(sx > 0, 0, 1)
    xs = np.array([+X_side, -X_side])
    w = v = th = thd = wc = vc = thc = thcd = 0.0
    zr = vr = 0.0; I_pu = I_pd = 0.0
    ref_u = max(float(F_cmd.max()), 1.0); ref_d = max(float(-F_cmd.min()), 1.0)
    ebuf = np.zeros((nd + kw + 2, 2)); integ = np.zeros(2); F_reg = np.zeros(2)
    out = dict(gap_b=G0, gap_t=G0, roll=0.0, move_c=0.0, err=0.0, contact=False, t_contact=float("nan"), pk=np.zeros(nf), heat=np.zeros(nf), pk_reg=np.zeros(2),
               heat_plan=0.0)
    hist = dict(t=t, rel=np.zeros(nt), rth=np.zeros(nt), F=np.zeros((nt, nf)), zr=np.zeros(nt)) if keep_hist else None
    side_err = np.array(side_err, dtype=float)
    for j in range(nt):
        tj = t[j]
        rel = w - wc; rth = th - thc
        # the plan as the coils deliver it, and the reference path
        I_pu += (math.sqrt(max(F_cmd[j], 0.0) / ref_u) - I_pu) * lag; I_pd += (math.sqrt(max(-F_cmd[j], 0.0) / ref_d) - I_pd) * lag
        F_plan = ref_u * I_pu**2 - ref_d * I_pd**2
        vr += ((-q_pred[j] + F_plan) / MU + F_plan / MU_C - 2 * omr * vr - omr**2 * zr) * dt; zr += vr * dt
        ebuf = np.roll(ebuf, 1, axis=0); ebuf[0] = rel + rth * xs - zr
        alive = tj < flost
        # the regulator, each side on its own gap
        if j % ks == 0 and j >= nd + kw:
            z = ebuf[nd] + side_err; zd = (ebuf[nd] - ebuf[nd + kw]) / (kw * dt)
            integ += z * ks * dt
            F_reg = -m_side * (om**2 * z + 2 * om * zd) - ki * integ
        # the plan's share of each side
        face = fsgn * F_plan > 0                          # the tracks of the face the plan is using now
        plan_side = np.zeros(2)
        if plan != "dropped" or alive[face].sum() == face.sum():
            if plan == "moved" and alive[face].sum() < face.sum() and alive[face].any():
                # the lost share goes to the tracks left, in proportion to their own shares, up to their ratings
                tot = fshare[face & alive].sum()
                for k in (0, 1):
                    m = face & alive & (fside == k)
                    plan_side[k] = np.sign(F_plan) * min(abs(F_plan) * fshare[m].sum() / tot, fcap[m].sum()) if m.any() else 0.0
            elif plan == "balanced" and alive[face].sum() < face.sum():
                capk = [fcap[face & alive & (fside == k)].sum() for k in (0, 1)]
                each = min(abs(F_plan) / 2, capk[0], capk[1])
                plan_side[:] = np.sign(F_plan) * each
            else:
                for k in (0, 1):
                    plan_side[k] = F_plan * fshare[face & alive & (fside == k)].sum()
        # what each side wants in all, to the tracks of the right face on that side
        tgt = np.zeros(nf)
        for k in (0, 1):
            c = plan_side[k] + F_reg[k]
            m = alive & (fside == k) & (fsgn * c > 0)
            if m.any():
                capm = fcap[m].sum(); cc = min(abs(c), capm)
                if abs(c) >= capm:
                    integ[k] -= (ebuf[nd, k] + side_err[k]) * dt        # no wind-up while the side is at its rating
                tgt[m] = cc * fcap[m] / capm
            out["pk_reg"][k] = max(out["pk_reg"][k], abs(F_reg[k]))
        s_t = np.sqrt(tgt)
        s += np.clip(s_t - s, -rate * dt, rate * dt)
        s = np.where(alive, s, 0.0)
        # the forces
        z_f = rel + rth * fx
        gap_f = np.where(fsgn > 0, G0 - z_f, G0 + z_f)
        err_f = np.where(fsgn > 0, -side_err[fside], side_err[fside])         # the regulator sets the current for the gap it measures
        Ff = np.minimum(s**2 * ((gap_f + err_f) / np.maximum(gap_f, 1e-3))**2, fcap)
        F_up = float((fsgn * Ff).sum()); M = float((fsgn * Ff * fx).sum())
        z_s = rel + rth * sx
        Ss = np.where(tj < slost, s0, 0.0)
        if sog is not None:
            Ss = Ss * np.interp(G0 - z_s, sog[0], sog[1])
        if couple:
            for k in (0, 1):
                m = (fside == k) & (fsgn > 0)
                Ss = np.where(sside == k, Ss * (1.0 - couple * math.sqrt(min(Ff[m].sum() / max(ffull[m].sum(), 1.0), 4.0))), Ss)
        F_up += float(Ss.sum()); M += float((Ss * sx).sum())
        a = (F_up - W - q[j]) / MU; al = M / I_cas
        ac = -(F_up - W) / MU_C; alc = -M / I_cab
        v += a * dt; w += v * dt; thd += al * dt; th += thd * dt
        vc += ac * dt; wc += vc * dt; thcd += alc * dt; thc += thcd * dt
        rel = w - wc; rth = th - thc
        out["roll"] = max(out["roll"], abs(rth)); out["move_c"] = max(out["move_c"], abs(rel)); out["err"] = max(out["err"], abs(rel - zr))
        np.maximum(out["pk"], Ff, out=out["pk"]); out["heat"] += Ff / ffull * dt; out["heat_plan"] += abs(F_plan) / max(ffull[face].sum(), 1.0) * dt
        for y in (-edge, edge):
            r_y = rel + rth * y
            out["gap_b"] = min(out["gap_b"], G0 - r_y); out["gap_t"] = min(out["gap_t"], G0 + r_y)
        if keep_hist:
            hist["rel"][j] = rel; hist["rth"][j] = rth; hist["F"][j] = fsgn * Ff; hist["zr"][j] = zr
        if out["gap_b"] <= 0 or out["gap_t"] <= 0:
            out["contact"] = True; out["t_contact"] = float(tj); break
    out["gap"] = max(min(out["gap_b"], out["gap_t"]), 0.0)
    out["move_edge"] = G0 - out["gap"]
    out["pk"] = out["pk"].tolist(); out["heat"] = out["heat"].tolist(); out["pk_reg"] = out["pk_reg"].tolist()
    out["names"] = [f["name"] for f in fast]
    return (out, hist) if keep_hist else out
