#!/usr/bin/env python3
"""
The gap regulator along the track, heave and roll together (round 11, 2026-11-17).
Book: Orbital Ring Engineering, Volume III: Chapter 7.  Brief: Vol_III/Claude_Code_instructions_round11_2026-11-17.md, item 1.

gap_regulator.py is one cross-section per meter: nothing along the track.  This is the same cross-section at every node of a stretch of ring,
with what joins the nodes:

  casing   a beam in bending (the Timoshenko truss of the beam lanes, sled_passage_simulation.truss(10, 0.1)) and in torsion (GJ, stiffness()).
  cable    a string under its tension in heave, and in roll the torsional stiffness its tension gives it, T x (I_p / A) = T a^2 / 6 for a square of
           side a (a bundle of parallel fibers under tension: a fiber at radius r is tilted r dtheta/dx by a twist, and its tension then has a
           component r dtheta/dx across the radius, a torque sigma dA r^2 dtheta/dx; summed over the section, T I_p / A).  The cable's own shear
           stiffness is left out.
  coils    each load and lift track is N_SETS sets of coils L_COIL long, a new coil every L_COIL / N_SETS, each coil with its own current, its own
           rating (the track's over N_SETS) and its own slew limit.  A coil can be dead from a given time.
  plan     the planned pull of the launch, taken from the point model's lane as a history at one place and moved along the track at the sled's
           speed (so are the sled's load and the load the plan expects).
  regulator  sled_passage_pulsed.regulator() at every node on each side of the cable, as in gap_regulator.section_run.  A coil is given the mean,
           over its own length, of what its side asks for, a twentieth of it, or more when the controller knows that coils over it are dead
           (known=True: the live coils over a node share what the node asks for, each inside its rating).
  supports the four support loops at constant flux, each on its own gap, with the load tracks' coupling.

Free ends.  The scores are taken over the middle of the stretch (margin).
"""
import math

import numpy as np
import scipy.sparse as sps

import ring_altitude as ring
import sled_passage_pulsed as pz
import sled_passage_simulation as sp
import gap_regulator as gr

N_SETS = 20
L_COIL = 25.0
TRACKS = ("lift +", "load +", "lift -", "load -")          # index 0..3; the + side is at +X_C
FSGN = np.array([-1.0, 1.0, -1.0, 1.0])                    # a load track pulls the casing up, a lift track down
FSIDE = np.array([0, 0, 1, 1])


def stiffness(truss=None, side=10.0, torsion_scale=1.0):
    """The casing's and the cable's stiffnesses and where they come from.
    Bending: the beam lanes' truss (10 m deep, chords with a tenth of the casing's mass, CNT).
    Torsion, a first figure: a thin-walled square tube of the casing's outline (side x side, the 10 m of the design rule taken as height and width
    [UNVERIFIED, as bearing_design F.3]), each of its four walls a lattice with half the shear stiffness of the beam lanes' web, G t h = GA / 2.
    Bredt: GJ = 4 A^2 / (perimeter / (G t)) = side^3 G t."""
    tr = truss or sp.truss(10.0, 0.1)
    Gt = (tr["GA"] / 2.0) / side
    GJ = side**3 * Gt * torsion_scale
    a = ring.CABLE_SIDE
    return dict(EI=tr["EI"], GA=tr["GA"], J=tr["J"], c_s=tr["c_s"], c_b=tr["c_b"], Gt=Gt, side=side, GJ=GJ, I_cas=pz.MU * gr.R_GYR**2,
                GJ_cab=sp.T_CABLE * a * a / 6.0, I_cab=pz.MU_C * a * a / 6.0, T=sp.T_CABLE, torsion_scale=torsion_scale, truss=tr["label"])


def grid(L, dx):
    """Nodes and coils.  Coil (s, j) of a track covers [s L_COIL / N_SETS + (j - 1) L_COIL, + L_COIL)."""
    N = int(round(L / dx)) + 1
    x = np.arange(N) * dx
    n_sec = int(math.ceil(L / L_COIL)) + 1
    rows = []; cols = []
    for s in range(N_SETS):
        j = np.floor((x - s * L_COIL / N_SETS) / L_COIL).astype(int) + 1
        rows.append(s * n_sec + j); cols.append(np.arange(N))
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    P = sps.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N_SETS * n_sec, N))
    cnt = np.asarray(P.sum(axis=1)).ravel()
    A = sps.diags(1.0 / np.maximum(cnt, 1.0)) @ P
    x0 = (np.arange(N_SETS)[:, None] * L_COIL / N_SETS + (np.arange(n_sec)[None, :] - 1) * L_COIL).ravel()
    return dict(N=N, x=x, dx=dx, n_sec=n_sec, n_coil=N_SETS * n_sec, P=P, PT=P.T.tocsr(), A=A.tocsr(), x0=x0, x1=x0 + L_COIL, xmid=x0 + L_COIL / 2, set=np.repeat(np.arange(N_SETS), n_sec))


def outage(g, track, x_a, x_b, sets=None, t=-np.inf, t_lost=None):
    """Mark coils of a track dead from time t: the coils that overlap [x_a, x_b] (x_a = x_b: the coils over that point), of the given sets
    (None: all).  Returns the (4, n_coil) array of the times of loss."""
    if t_lost is None:
        t_lost = np.full((4, g["n_coil"]), np.inf)
    m = (g["x1"] > x_a) & (g["x0"] <= x_b)
    if sets is not None:
        m &= np.isin(g["set"], sets)
    t_lost[track, m] = t
    return t_lost


def run(h_t, q, q_pred, F_cmd, F_coil, u, reg, st, F_full, cap, X_C, X_S, L=3000.0, dx=6.25, t_lost=None, known=False, plan="moved", tau=10e-3, couple=0.0,
        t0=None, t1=None, margin=0.15, sup_off=0.225, W=None, snap_t=None):
    """One pass along a stretch L of ring.  h_t, q, q_pred, F_cmd: the point lane's histories at one place (the sled's load, the load the plan
    expects, the planned command, positive up).  F_coil: the planned command of one coil (times N_SETS), in the time of the coil's middle: the
    point lane run with one coil in place of the staggered sets.  F_cmd is the sum of the staggered coils at a place, and is what the reference
    path is made from.  F_full, cap: one track's full pull and rating, N/m.  reg: a regulator() dict for one side.
    t_lost: outage() array, in the time of the place x = L / 2 (0: the spacecraft's front arrives there).  known: the controller knows which coils
    are dead from the moment they are.  plan, where a node has no live coil of the planned face on one side and the controller knows it: 'moved'
    (the other side's track takes the whole plan up to its rating), 'kept' (the other side keeps its half), 'dropped' (no planned pull there)."""
    W = pz.W if W is None else W; MU = pz.MU; MU_C = pz.MU_C; G0 = pz.GAP0; sog = pz.SUPPORT_OF_GAP
    g = grid(L, dx); N = g["N"]; x = g["x"]; A = g["A"]; PT = g["PT"]; nc = g["n_coil"]
    xm = L / 2.0
    c_max = max(st["c_b"], st["c_s"], math.sqrt(st["T"] / MU_C), math.sqrt(st["GJ"] / st["I_cas"]), math.sqrt(st["GJ_cab"] / st["I_cab"]))
    ks = int(math.ceil(reg["sample"] / (0.45 * dx / c_max))); dt = reg["sample"] / ks
    nd = int(round(reg["delay"] / dt)); kw = max(int(round(reg["window"] / dt)), 1); nb = nd + kw + 2
    t0 = float(h_t[0]) if t0 is None else t0; t1 = float(h_t[-1]) if t1 is None else t1
    ta = t0 + 0.0; tb = t1 + L / u                    # in the time of the place x = 0
    nt = int(math.ceil((tb - ta) / dt)) + 1
    lag = 1.0 - math.exp(-dt / tau)
    om = reg["om"]; m_side = reg["mass"]; ki = reg["ki"]; omr = om / 10.0
    ref_u = max(float(F_cmd.max()), 1.0); ref_d = max(float(-F_cmd.min()), 1.0)
    rate = math.sqrt(F_full / N_SETS) / reg["rise"]; cap_coil = cap / N_SETS
    t_lost = np.full((4, nc), np.inf) if t_lost is None else t_lost + xm / u          # to the time of x = 0
    edge = ring.CABLE_SIDE / 2.0
    sx = np.array([X_S - sup_off, X_S + sup_off, -X_S - sup_off, -X_S + sup_off]); sside = np.array([0, 0, 1, 1])
    EI, GA, J = st["EI"], st["GA"], st["J"]
    w = np.zeros(N); v = np.zeros(N); phi = np.zeros(N - 1); phid = np.zeros(N - 1); th = np.zeros(N); thd = np.zeros(N)
    wc = np.zeros(N); vc = np.zeros(N); thc = np.zeros(N); thcd = np.zeros(N)
    Vp = np.zeros(N + 1); Mp = np.zeros(N); Tq = np.zeros(N + 1); Tqc = np.zeros(N + 1)
    Ic_u = np.zeros(nc); Ic_d = np.zeros(nc); xmid = g["xmid"]
    I_pu = np.zeros(N); I_pd = np.zeros(N); zr = np.zeros(N); vr = np.zeros(N)
    buf = np.zeros((nb, 2, N)); integ = np.zeros((2, N)); F_reg = np.zeros((2, N)); sat = np.zeros((2, N), dtype=bool)
    s = np.zeros((4, nc))
    alive = np.ones((4, nc), dtype=bool); n_al = np.full((4, N), float(N_SETS)); n_dead = -1
    i0 = int(margin * N); i1 = N - i0; reg_sl = slice(i0, i1)
    out = dict(gap=G0, roll=0.0, move_c=0.0, err=0.0, contact=False, t_contact=float("nan"), x_gap=float("nan"), t_gap=float("nan"), pk_coil=np.zeros(4), pk_node=np.zeros(4),
               dt=dt, N=N, nt=nt, dx=dx)
    snap = None
    for j in range(nt):
        t = ta + j * dt
        tn = t - x / u
        Fc = np.interp(tn, h_t, F_cmd); qn = np.interp(tn, h_t, q); qp = np.interp(tn, h_t, q_pred)
        I_pu += (np.sqrt(np.maximum(Fc, 0.0) / ref_u) - I_pu) * lag; I_pd += (np.sqrt(np.maximum(-Fc, 0.0) / ref_d) - I_pd) * lag
        F_plan = ref_u * I_pu**2 - ref_d * I_pd**2
        vr += ((-qp + F_plan) / MU + F_plan / MU_C - 2 * omr * vr - omr**2 * zr) * dt; zr += vr * dt
        rel = w - wc; rth = th - thc
        buf[j % nb, 0] = rel + rth * X_C - zr; buf[j % nb, 1] = rel - rth * X_C - zr
        al_now = t < t_lost
        nd_now = int(al_now.size - al_now.sum())
        if nd_now != n_dead:
            alive = al_now; n_dead = nd_now
            n_al = (PT @ alive.T.astype(float)).T
        if j % ks == 0 and j >= nd + kw:
            z = buf[(j - nd) % nb]; zd = (z - buf[(j - nd - kw) % nb]) / (kw * dt)
            integ += np.where(sat, 0.0, z * ks * dt)
            F_reg = -m_side * (om**2 * z + 2 * om * zd) - ki * integ
        # each coil's own planned command, through its rise, and the regulator's part as the mean over the coil's length
        Fcc = np.interp(t - xmid / u, h_t, F_coil)
        Ic_u += (np.sqrt(np.maximum(Fcc, 0.0) / ref_u) - Ic_u) * lag; Ic_d += (np.sqrt(np.maximum(-Fcc, 0.0) / ref_d) - Ic_d) * lag
        Fp = ref_u * Ic_u**2 - ref_d * Ic_d**2                                           # (n_coil,), the plan of the pair of tracks at this coil
        Rc = (A @ F_reg.T).T                                                             # (2, n_coil)
        if known and n_dead:
            inv = (A @ (N_SETS / np.maximum(n_al.T, 1.0))).T                             # (4, n_coil): how much more a live coil gives where coils over it are dead
            av = (A @ (n_al.T > 0).astype(float)).T > 0.5                                # (4, n_coil): the track has live coils here
            up = Fp > 0
            av0 = np.where(up, av[1], av[0]); av1 = np.where(up, av[3], av[2]); both = av0 & av1
            if plan == "moved":
                one = np.sign(Fp) * np.minimum(np.abs(Fp), cap)
                p0 = np.where(both, Fp / 2, np.where(av0, one, 0.0)); p1 = np.where(both, Fp / 2, np.where(av1, one, 0.0))
            elif plan == "kept":
                p0 = np.where(av0, Fp / 2, 0.0); p1 = np.where(av1, Fp / 2, 0.0)
            else:
                p0 = np.where(both, Fp / 2, 0.0); p1 = p0
        else:
            inv = 1.0; p0 = p1 = Fp / 2
        c0 = p0 + Rc[0]; c1 = p1 + Rc[1]
        tgt = np.stack([np.maximum(-c0, 0.0), np.maximum(c0, 0.0), np.maximum(-c1, 0.0), np.maximum(c1, 0.0)]) * inv / N_SETS     # (4, n_coil)
        hit_c = (tgt >= cap_coil).astype(float)
        tgt = np.minimum(tgt, cap_coil)
        hit = (PT @ hit_c.T) >= np.maximum(n_al.T, 1.0) - 0.5                            # (N, 4): every live coil over the node is at its rating
        sat[0] = hit[:, 0] | hit[:, 1] | ((F_reg[0] < 0) & (n_al[0] == 0)) | ((F_reg[0] > 0) & (n_al[1] == 0))        # no wind-up while a side is at its
        sat[1] = hit[:, 2] | hit[:, 3] | ((F_reg[1] < 0) & (n_al[2] == 0)) | ((F_reg[1] > 0) & (n_al[3] == 0))        # rating or has no coil for the pull
        s += np.clip(np.sqrt(tgt) - s, -rate * dt, rate * dt)
        s *= alive
        s2 = s * s
        Ff = (PT @ s2.T).T                                                               # (4, N), each track's pull at each node
        F_up = Ff[1] + Ff[3] - Ff[0] - Ff[2]
        M = X_C * (Ff[1] - Ff[0]) - X_C * (Ff[3] - Ff[2])
        for i in range(4):
            Ss = W / 4.0
            if sog is not None:
                Ss = Ss * np.interp(G0 - (rel + rth * sx[i]), sog[0], sog[1])
            if couple:
                Fl = Ff[1] if sside[i] == 0 else Ff[3]
                Ss = Ss * (1.0 - couple * np.sqrt(np.minimum(Fl / F_full, 4.0)))
            F_up = F_up + Ss; M = M + Ss * sx[i]
        # the casing: bending and torsion
        V = GA * ((w[1:] - w[:-1]) / dx - phi)
        Vp[1:-1] = V
        Mp[1:-1] = EI * (phi[1:] - phi[:-1]) / dx
        a_phi = ((Mp[1:] - Mp[:-1]) / dx + V) / J
        Tq[1:-1] = st["GJ"] * (th[1:] - th[:-1]) / dx
        a = ((Vp[1:] - Vp[:-1]) / dx + F_up - W - qn) / MU
        al = ((Tq[1:] - Tq[:-1]) / dx + M) / st["I_cas"]
        # the cable: a string, and its tension's torsional stiffness
        f_cab = np.zeros(N); f_cab[1:-1] = st["T"] * (wc[2:] - 2 * wc[1:-1] + wc[:-2]) / dx**2
        Tqc[1:-1] = st["GJ_cab"] * (thc[1:] - thc[:-1]) / dx
        ac = (f_cab - (F_up - W)) / MU_C
        alc = ((Tqc[1:] - Tqc[:-1]) / dx - M) / st["I_cab"]
        v += a * dt; w += v * dt; phid += a_phi * dt; phi += phid * dt; thd += al * dt; th += thd * dt
        vc += ac * dt; wc += vc * dt; thcd += alc * dt; thc += thcd * dt
        rel = (w - wc)[reg_sl]; rth = (th - thc)[reg_sl]
        r_e = np.maximum(np.abs(rel + rth * edge), np.abs(rel - rth * edge))
        k = int(np.argmax(r_e)); gp = G0 - float(r_e[k])
        if gp < out["gap"]:
            out["gap"] = gp; out["x_gap"] = float(x[i0 + k] - xm); out["t_gap"] = float(t - x[i0 + k] / u)
        out["roll"] = max(out["roll"], float(np.abs(rth).max())); out["move_c"] = max(out["move_c"], float(np.abs(rel).max()))
        out["err"] = max(out["err"], float(np.abs(rel - zr[reg_sl]).max()))
        np.maximum(out["pk_coil"], s2.max(axis=1), out=out["pk_coil"]); np.maximum(out["pk_node"], Ff[:, reg_sl].max(axis=1), out=out["pk_node"])
        if snap_t is not None and snap is None and t - xm / u >= snap_t:
            snap = dict(x=x - xm, rel=(w - wc).copy(), rth=(th - thc).copy(), F=Ff.copy())
        if gp <= 0:
            out["contact"] = True; out["t_contact"] = out["t_gap"]; break
    out["gap"] = max(out["gap"], 0.0)
    out["pk_coil_share"] = (out["pk_coil"] / cap_coil).tolist()              # the hardest worked coil of each track, as a share of its rating
    out["pk_coil"] = out["pk_coil"].tolist(); out["pk_node"] = out["pk_node"].tolist()
    if snap is not None:
        out["snap"] = snap
    return out
