#!/usr/bin/env python3
"""
Round 12, items 1 to 4: helpers for tools/round12_stator.py (the 6 m stator at 2,500 Hz).
Book: Orbital Ring Engineering, Volume III: Table 9.6, Sections 9.3, 9.10, 9.11, 10.1, 14.5.
Run: imported by tools/round12_stator.py; `python tools/round12_stator_model.py` runs its self-checks.

It builds on tools/round11_stator_model.py (the winding, the fields at the tapes, the strip loss) and adds:

  inductance   the mutual inductance of two loops in parallel planes a distance D apart (the other side of the stator, the
               sled's coils), by the same Neumann integral over straight filaments;
  the sled     the field of a row of sled coils by Biot-Savart: a plain racetrack per pole, or a stepped winding of n nested
               loops per pole whose cross conductors are placed so that the current sheet is a sine wave in n equal steps;
  the launch   evaluate12: round 11's launch integral again (it reproduces round11_stator_model.evaluate), with three more
               losses along the launch and the inverter's rating:
                 eddy      in the tape's own metal for the field through the face: pi^2 f^2 B^2 w^3 (t / rho) / 6 per meter of
                           tape (a strip of width w in a uniform alternating field; the induced field omega B x integrated over
                           the width), for a stated layer build;
                 soldered  the layers of a turn soldered into one stack, in the field along their faces: the bearing's model
                           (bearing_round5.py, Section 5): the stack as one Bean slab of the stack's thickness with the stack's
                           engineering current density, at the pitch of 0.107 mm per layer;
                 film      the same field on insulated layers: the film's own slab loss (3.5 micrometers);
               and the turns, the terminal voltage, the real and apparent power and the inverter rating (largest terminal
               voltage times largest current) of a coil.
"""
import functools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["ORBITAL_RING_LSM_STATOR"] = "table96"   # round 13: this study was run on the old 130 m stator of lsm_config.py (its launch and its
                                                   # constants); the switch keeps it reproducing now that the 6 m stator is the default
import numpy as np

import lsm_config as cfg
import lim_hysteresis as lh
import round11_stator_model as m

MU0 = m.MU0
PITCH = lh.TAPE_T / 0.75            # m per soldered layer: bearing_round5.TAPE_PITCH = bearing_design.TAPE_T / bearing_electrical.F_CABLE
T_FILM = 3.5e-6                     # m, the superconducting film (bearing_round5.py's figure)
RHO_CU = 2.0e-9                     # ohm m, copper near 70 to 77 K at RRR about 100 [UNVERIFIED: handbook value from memory]
G_BARE = m.EDDY_G                   # siemens: 80 um Hastelloy at 1.24 uohm m and 2 um silver at 3 nohm m (round 11's build)
BUILDS = (("no copper", 0.0), ("20 um copper", 20e-6), ("40 um copper", 40e-6))


def sheet_conductance(t_cu):
    return G_BARE + t_cu / RHO_CU


def q_bean_slab(Ba, Bp):
    """Bean slab in a parallel field of amplitude Ba, penetration field Bp: loss per cycle per m3 (bearing_round5.q_bean_slab, vectorized)"""
    Ba = np.asarray(Ba, float); Bp = np.asarray(Bp, float)
    return np.where(Ba <= Bp, (2.0 / 3.0) * Ba**3 / (MU0 * Bp), (2.0 / MU0) * (Ba * Bp - 2.0 * Bp**2 / 3.0))


# ------------------------------------------------------------------------------------------------ inductance between parallel planes
def m_loops_z(X, tau, w, g_act, g_end, D, tau2=None):
    """mutual inductance (H) of two one-turn rectangular loops in parallel planes a distance D apart: the first tau x w, the
    second tau2 x w (default tau) shifted X along the track.  With D = 0 it is round11_stator_model.m_loops."""
    X = np.asarray(X, float); tau2 = tau if tau2 is None else tau2
    runs = 2 * m.m_par(0.0, tau, X, X + tau2, math.hypot(g_end, D)) - 2 * m.m_par(0.0, tau, X, X + tau2, math.hypot(w, D))
    def side(dx):
        return m.m_par(-w / 2, w / 2, -w / 2, w / 2, np.sqrt(np.maximum(np.abs(dx), g_act)**2 + D * D))
    return runs + side(X) + side(X + tau2 - tau) - side(X + tau2) - side(X - tau)


def cyclic(tau, w, q, g_act, g_end, D=0.0, n_pp=12):
    """per-turn-squared cyclic inductance of a coil with the coils of a winding in a plane D away (D = 0: its own side)"""
    nq = max(1, int(round(3 * q))); j = np.arange(-nq * n_pp, nq * n_pp + 1)
    M = m_loops_z(j * tau / nq, tau, w, g_act, g_end, D)
    return float(np.sum(M * np.cos(math.pi * j / nq)))


def sled_coupling(tau, w, q_m, g_act, g_end, gap, g_sled=0.05, n_pp=12):
    """The stator winding and a row of one-turn racetrack sled coils (tau x w, alternate polarity, a distance gap away):
    M_sf  the flux through one stator turn per ampere in the sled's coils (coil aligned with a sled pole);
    M_fs  the flux through one sled coil per ampere-turn in each stator coil, the stator's wave aligned with the sled (d axis);
    L_f   the cyclic inductance of a sled coil in its row."""
    p = np.arange(-n_pp, n_pp + 1); sgn = np.where(p % 2 == 0, 1.0, -1.0)
    M_sf = float(np.sum(sgn * m_loops_z(p * tau, tau, w, g_act, g_end, gap)))
    nq = max(1, int(round(3 * q_m))); j = np.arange(-nq * n_pp, nq * n_pp + 1)
    M_fs = float(np.sum(m_loops_z(j * tau / nq, tau, w, g_act, g_end, gap) * np.cos(math.pi * j / nq)))
    L_f = float(np.sum(sgn * m.m_loops(p * tau, tau, w, g_sled, g_sled)))
    return M_sf, M_fs, L_f


# ------------------------------------------------------------------------------------------------ the sled's coils (Biot-Savart)
def seg_field(P, A, B):
    r1 = P - A; r2 = P - B
    n1 = np.linalg.norm(r1, axis=1); n2 = np.linalg.norm(r2, axis=1)
    cr = np.cross(r1, r2); den = n1 * n2 * (n1 * n2 + np.einsum("ij,ij->i", r1, r2))
    return 1e-7 * cr * ((n1 + n2) / np.maximum(den, 1e-30))[:, None]


def sled_loops(tau, n_steps, pack=0.2):
    """the loops of one pole centred on x = 0: (half length along the track, half width, share of the pole's ampere-turns).
    n_steps = 1: a plain racetrack, its pack laid flat and `pack` wide (five filaments).  n_steps > 1: n nested loops of equal
    current whose cross conductors sit where a sine-wave current sheet has carried (j - 1/2) / n of its current, counted
    from the pole's centre: a = (tau / pi) acos(1 - (j - 1/2) / n); their long sides are spread over `pack` about y = 1 m."""
    if n_steps == 1:
        offs = (np.arange(5) + 0.5) / 5 * pack - pack / 2 if pack > 0 else np.array([0.0])
        return [(tau / 2 - pack / 2 + o, 1.0 + o, 1.0 / len(offs)) for o in offs]
    out = []
    for j in range(1, n_steps + 1):
        a = tau / math.pi * math.acos(1.0 - (j - 0.5) / n_steps)
        out.append((a, 1.0 - pack / 2 + pack * (j - 0.5) / n_steps, 1.0 / n_steps))
    return out


def sled_field(tau, loops, z, ys, n_poles=6, nx=None):
    """Bx, By, Bz (len(ys), nx) per ampere-turn of a pole at the plane z below the coils, over one period along the track"""
    nx = min(2400, max(240, int(2 * tau / 0.025))) if nx is None else nx
    xs = (np.arange(nx) + 0.5) / nx * 2 * tau
    X, Y = np.meshgrid(xs, ys); P = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z)], axis=1)
    B = np.zeros_like(P)
    for j in range(-n_poles, n_poles + 2):
        sgn = 1.0 if j % 2 == 0 else -1.0; xc = (j + 0.5) * tau
        for a, hw, share in loops:
            c = [np.array([xc - a, -hw, 0.0]), np.array([xc + a, -hw, 0.0]), np.array([xc + a, hw, 0.0]), np.array([xc - a, hw, 0.0])]
            for k in range(4):
                B += sgn * share * seg_field(P, c[k], c[(k + 1) % 4])
    B = B.reshape(len(ys), nx, 3)
    return B[:, :, 0], B[:, :, 1], B[:, :, 2]


@functools.lru_cache(maxsize=None)
def sled_ratios(tau, n_steps, pack=0.2, g=m.G):
    """the sled's field at the stator as ratios to B0 (the fundamental of the field through the face averaged over the 2 m
    width): peaks over a cycle averaged over the tape positions, the fundamentals, the low harmonics, and the conductor"""
    loops = sled_loops(tau, n_steps, pack)
    ys = np.linspace(-0.95, 0.95, 13)
    Bx, By, Bz = sled_field(tau, loops, -g, ys)
    def harm(A, k):
        return 2 * np.abs(np.fft.rfft(A, axis=1)[:, k]) / A.shape[1]
    B0 = harm(Bz.mean(axis=0, keepdims=True), 1)[0]
    yb = np.linspace(0.5, 1.5, 11)
    _, By0, Bz0 = sled_field(tau, loops, -g, yb)
    Bx1, By1, Bz1 = sled_field(tau, loops, -g - 1.0, yb)
    ya = np.array([-1.0, 1.0])
    Bxl = np.concatenate([sled_field(tau, loops, -g - zz, ya)[0] for zz in (0.1, 0.3, 0.5, 0.7, 0.9)])     # along-track field down a leg
    out = dict(tau=tau, n_steps=n_steps, pack=pack, B0_per_At=float(B0), At_for_0p1T=float(0.1 / B0),
               rz_peak=float(np.abs(Bz).max(axis=1).mean() / B0), rx_fund=float(harm(Bx, 1).mean() / B0), rx_peak=float(np.abs(Bx).max(axis=1).mean() / B0),
               cz_face=float(np.abs(Bz0).max(axis=1).mean() / B0), cz_back=float(np.abs(Bz1).max(axis=1).mean() / B0),
               cy_face=float(np.abs(By0).max(axis=1).mean() / B0), cy_back=float(np.abs(By1).max(axis=1).mean() / B0),
               leg_x=float(np.abs(Bxl).max(axis=1).mean() / B0),
               hz=[float(harm(Bz.mean(axis=0, keepdims=True), k)[0] / B0) for k in (3, 5, 7, 9, 11)],
               hx=[float(harm(Bx, k).mean() / B0) for k in (3, 5, 7, 9, 11)],
               turn_m=float(sum(share * 2 * (2 * a + 2 * hw) for a, hw, share in loops)))
    out["Am_per_m"] = out["At_for_0p1T"] * out["turn_m"] / tau          # ampere-meters of sled conductor per meter of sled, one side
    return out


def apply_sled(st, q):
    """put a Biot-Savart sled's peak fields at the tapes of the stator st in place of the sine wave's"""
    st.rx = q["rx_peak"]; st.rz = q["rz_peak"]
    if st.end_lay == "flat":
        st.c_run = q["cz_back"] if st.d > 0 else q["cz_face"]
    else:
        st.c_run = q["cy_face"]
    st.c_end = q["cz_face"]
    st.leg_x = q["leg_x"]
    return st


# ------------------------------------------------------------------------------------------------ the launch with everything in it
def _region(BnK, BtK, K, Bn_s, Bt_s, i_t, w_t, T, d_stack):
    Bn = BnK * K + Bn_s; Bt = BtK * K + Bt_s
    B = np.hypot(Bn, Bt); th = np.degrees(np.arctan2(Bt, Bn))
    ic_min = float(m.ic(B.ravel(), th.ravel(), w_t, T).min())
    Bm = B.mean(axis=1); Bnm = Bn.mean(axis=1); Btm = Bt.mean(axis=1)
    thm = np.degrees(np.arccos(np.clip(Bnm / np.maximum(Bm, 1e-30), 0, 1)))
    Ic = m.ic(Bm, thm, w_t, T)
    qh = lh.q_strip_perp(Ic, w_t, Bnm) + lh.q_norris_strip(i_t, Ic)
    Bp = MU0 * (Ic / (w_t * PITCH)) * d_stack / 2.0
    qs = q_bean_slab(Btm, Bp) * (w_t * PITCH)
    Bpf = MU0 * (Ic / (w_t * T_FILM)) * T_FILM / 2.0
    qf = q_bean_slab(Btm, Bpf) * (w_t * T_FILM)
    return ic_min, float(qh.mean()), float(qs.mean()), float(qf.mean()), float((Bnm**2).mean()), float(Bp.mean()), float(Btm.mean()), float(Bnm.mean())


def _fields(st, K, B0):
    """(Bn_own per K, Bt_own per K, Bn_sled, Bt_sled) of the three regions"""
    own = np.maximum(st.aBt, st.B_avg) if st.lay == "edge" else np.maximum(st.aBn, st.B_avg)
    act = (st.aBn, own, B0 * st.rx, B0 * st.rz) if st.lay == "edge" else (own, st.aBt, B0 * st.rz, B0 * st.rx)
    end = (st.eBn, st.eBt, B0 * st.c_run, B0 * st.c_run)
    leg_x = getattr(st, "leg_x", None)
    leg = (st.eBn, st.eBt, B0 * (st.rx * st.leg_decay if leg_x is None else leg_x), B0 * st.c_end * st.leg_decay) if st.d > 0 else end
    return act, end, leg


def evaluate12(st, F_pk, P, L_sled=10000.0, m_sc=cfg.M_SPACECRAFT, v_f=cfg.V_LAUNCH, alpha=0.2, delta=m.DELTA, phi_scale=1.0, I_coil=cfg.I_PEAK, n_s=8, n_v=3000):
    """One stator on one launch (the arguments of round11_stator_model.evaluate; delta the load angle; phi_scale multiplies
    the stator's own flux per turn, for a different inductance).  Returns round 11's quantities and the added ones."""
    tau = st.tau; W = m.W; KD = m.KD; B0N = m.B0N
    L_act = L_sled * m.ACTIVE; mass = cfg.M_SLED_PER_M * L_sled + m_sc
    K1 = 2 * F_pk / (st.sides * B0N * math.sin(delta) * W * L_act); Kb = K1 / KD
    f3 = _fields(st, Kb, B0N)
    ic_min = min(_region(a, b, Kb, c, d_, 1.0, st.w_t, st.T, 1e-3)[0] for a, b, c, d_ in f3)
    i_t = m.SHARE * ic_min; n_t = Kb / i_t
    layers = max(1, int(round(I_coil / i_t))); I_pk = layers * i_t; d_stack = layers * PITCH
    v_c = P / F_pk; s_min = min(1.0, v_c / v_f)
    sg = np.geomspace(s_min, 1.0, n_s) if s_min < 1 else np.array([1.0, 1.0])
    lens = (W, tau, 2 * st.d)
    rows = []
    for s in sg:
        K = Kb * s**(1 - alpha); B0 = B0N * s**alpha
        hy = sl = fi = b2 = 0.0; parts = []
        for (a, b, c, d_), ln in zip(_fields(st, K, B0), lens):
            _, qh, qs, qf, Bn2, Bp, Btm, Bnm = _region(a, b, K, c, d_, i_t * s**(1 - alpha), st.w_t, st.T, d_stack)
            hy += n_t * ln * qh; sl += n_t * ln * qs; fi += n_t * ln * qf; b2 += n_t * ln * Bn2; parts.append(n_t * ln * qh)
        rows.append([hy, sl, fi, b2] + parts)
    rows = np.array(rows)
    ve = np.linspace(0.0, v_f, n_v + 1); v = 0.5 * (ve[1:] + ve[:-1]); dv = ve[1] - ve[0]
    Fv = np.minimum(F_pk, P / v); dt = mass * dv / Fv
    s = np.minimum(1.0, v_c / v); f = v / (2 * tau)
    def along(col):
        y = rows[:, col]
        if s_min >= 1:
            return np.full_like(v, y[-1])
        return np.exp(np.interp(np.log(s), np.log(sg), np.log(np.maximum(y, 1e-300))))
    scale = st.sides * L_act
    qh, qs, qf, b2 = along(0), along(1), along(2), along(3)
    low = v < v_c; x_c = float(np.sum((v * dt)[low])); n_start = max(x_c / m.SITE, 1.0)
    def tot(p):                                  # p: W per meter of track of one side along the launch
        return scale * float(np.sum(p * dt)), scale * float(np.sum((p * dt)[low]))
    E_h, E_h_low = tot(qh * f); E_s, E_s_low = tot(qs * f); E_f, _ = tot(qf * f)
    p_e1 = (math.pi**2 / 6.0) * f**2 * st.w_t**3 * b2                     # per siemens of sheet conductance
    E_e1, E_e1_low = tot(p_e1)
    eddy = {name: dict(E=E_e1 * sheet_conductance(t), E_low=E_e1_low * sheet_conductance(t)) for name, t in BUILDS}
    # the voltage of a turn along the launch
    B0v = B0N * s**alpha; Kv = Kb * s**(1 - alpha)
    e_t = 2 * KD * v * W * B0v
    x_t = (math.pi * v / tau) * st.phi_unit * phi_scale * Kv
    cd = math.cos(delta)
    v_dem = np.sqrt(e_t**2 + x_t**2 - 2 * e_t * x_t * cd); v_mag = np.sqrt(e_t**2 + x_t**2 + 2 * e_t * x_t * cd)
    N = int(m.VLIM // float(max(v_dem.max(), e_t.max())))
    q = Kb * tau / 6.0 / (max(N, 1) * I_pk)
    n_coils = st.sides * 3 * q * L_act / tau
    Iv = I_pk * s**(1 - alpha)
    Pc = 0.5 * N * e_t * Iv * math.sin(delta); Sc = 0.5 * N * v_dem * Iv
    kin_t = float(dt.sum()); kin_x = float((v * dt).sum())
    tape_m, tape_kg = st.tape_per_m(n_t)
    KE = 0.5 * mass * v_f**2
    k_c = int(np.searchsorted(v, v_c)) if v_c < v_f else n_v - 1
    def site_kw(k, G):
        return st.sides * L_act * float((qh[k] + qs[k]) * f[k] + p_e1[k] * G) / (L_sled / m.SITE) / 1e3
    return dict(tau=tau, lay=st.lay, end_lay=st.end_lay, d=st.d, b_e=st.b_e, w_t=st.w_t, T=st.T, sides=st.sides, alpha=alpha, delta_deg=math.degrees(delta),
                F_pk=F_pk, P=P, v_c=v_c, K1=K1, K_belt=Kb, i_t=i_t, n_t=n_t, layers=layers, I_pk=I_pk, d_stack=d_stack, N=N, q=q, n_coils=n_coils,
                f_end=v_f / (2 * tau), tape_kg_per_m=st.sides * tape_kg, tape_ring_t=st.sides * tape_kg * cfg.L_RING / 1e3,
                t=kin_t, x=kin_x, laps=kin_x / cfg.L_RING, KE=KE, x_c=x_c,
                E_hyst=E_h, E_hyst_low=E_h_low, E_sold=E_s, E_sold_low=E_s_low, E_film=E_f, eddy=eddy,
                E_act=scale * float(np.sum(along(4) * f * dt)), E_end=scale * float(np.sum(along(5) * f * dt)), E_leg=scale * float(np.sum(along(6) * f * dt)),
                V_max=N * float(v_dem.max()), v_at_Vmax=float(v[int(v_dem.argmax())]), E_max=N * float(e_t.max()), Vmag_max=N * float(v_mag.max()),
                S_max=float(Sc.max()) * n_coils, P_check=float(Pc[-1]) * n_coils, pf_min=float((Pc / np.maximum(Sc, 1e-30))[low.sum():].min()) if low.sum() < n_v else 1.0,
                rating_coil=0.5 * N * float(v_dem.max()) * I_pk, rating=0.5 * N * float(v_dem.max()) * I_pk * n_coils,
                site_kw_c=site_kw(k_c, G_BARE), site_kw_f=site_kw(n_v - 1, G_BARE), n_start=n_start,
                phi_unit=st.phi_unit * phi_scale, B_avg_pk=st.B_avg * Kb, rx=st.rx, rz=st.rz, c_run=st.c_run)


def total(r, build="no copper", soldered=True):
    """the heat of a launch with everything in it, and the part of it left below the corner speed"""
    E = r["E_hyst"] + r["eddy"][build]["E"] + (r["E_sold"] if soldered else r["E_film"])
    E_low = r["E_hyst_low"] + r["eddy"][build]["E_low"] + (r["E_sold_low"] if soldered else 0.0)
    return E, E_low


def phasors(r, speeds, st_tau=None, delta=None):
    """the phasor budget of a coil of the result r at the given speeds (m/s)"""
    tau = r["tau"]; delta = math.radians(r["delta_deg"]) if delta is None else delta
    out = []
    for v in speeds:
        s = min(1.0, r["v_c"] / v); a = r["alpha"]
        B0 = m.B0N * s**a; I = r["I_pk"] * s**(1 - a); f = v / (2 * tau)
        E = r["N"] * 2 * m.KD * v * m.W * B0
        X = r["N"] * (math.pi * v / tau) * r["phi_unit"] * r["K_belt"] * s**(1 - a)
        V = math.sqrt(E * E + X * X - 2 * E * X * math.cos(delta)); Vm = math.sqrt(E * E + X * X + 2 * E * X * math.cos(delta))
        Pc = 0.5 * E * I * math.sin(delta); Sc = 0.5 * V * I
        Qc = 0.5 * I * (X - E * math.cos(delta))                         # reactive power the supply gives: the coil's own less what the sled's flux returns
        phi_v = math.degrees(math.atan2(Qc, Pc))
        out.append(dict(v=v, f=f, B0=B0, I=I, E=E, X=X, V=V, V_mag=Vm, P=Pc, S=Sc, Q=Qc, pf=Pc / Sc, angle_V_I=phi_v, L_eff=X / (2 * math.pi * f * I)))
    return out


def selfcheck():
    print("SELF-CHECKS OF tools/round12_stator_model.py")
    st = m.Stator(6.0, "edge", 0.003, 70.0, 1.0, 1.0)
    F = 83294324.34775715; P = F * 5000.0
    a = m.evaluate(st, F, P, 10000.0, cfg.M_SPACECRAFT, cfg.V_LAUNCH, 0.2); b = evaluate12(st, F, P, alpha=0.2)
    print(f"  the launch integral: round 11 hysteresis {a['E_hyst']:.6e} J, here {b['E_hyst']:.6e} J; eddy (round 11's build) {a['E_eddy']:.4e}, here {b['eddy']['no copper']['E']:.4e} J;")
    print(f"  turns {a['N']}, {b['N']}; largest terminal voltage per turn {a['V_turn_max']:.2f}, here {b['V_max']/b['N']:.2f} V; tape {a['tape_kg_per_m']:.3f}, {b['tape_kg_per_m']:.3f} kg/m")
    print(f"  m_loops_z at D = 0 against round 11's m_loops: {float(m_loops_z(2.0, 6.0, 2.0, 0.1, 0.2, 0.0)):.6e}, {float(m.m_loops(2.0, 6.0, 2.0, 0.1, 0.2)):.6e} H")
    try:
        import bearing_round5 as r5
        print(f"  the slab formula against bearing_round5.q_bean_slab: {float(q_bean_slab(0.1, 0.9)):.6e}, {r5.q_bean_slab(0.1, 0.9):.6e}; above penetration {float(q_bean_slab(1.5, 0.9)):.6e}, {r5.q_bean_slab(1.5, 0.9):.6e} J/m3")
        print(f"  the layer pitch: here {PITCH*1e3:.4f} mm, bearing_round5.TAPE_PITCH {r5.TAPE_PITCH*1e3:.4f} mm")
    except Exception as e:                                                # the bearing module is heavy to import; the formula is two lines
        print(f"  bearing_round5 not imported ({type(e).__name__}); the slab formula is copied from it")
    q1 = sled_ratios(6.0, 1, 0.2)
    print(f"  the racetrack sled at 6 m, 0.2 m pack: peak along the track {q1['rx_peak']:.2f} B0, through the face {q1['rz_peak']:.2f} B0 (round 11's cross-check: 2.44 and 1.80 on its grid)")


if __name__ == "__main__":
    selfcheck()
