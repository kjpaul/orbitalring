#!/usr/bin/env python3
"""
Design history: the pulsed-track step of the sled passage work that led to the adopted bearing in `final_bearing.py`.
Pulsed electromagnetic tracks for the sled passage: staggered coils, the protection law, block sections, failure cases.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (library and the pulsed-track step).
Run: python sled_passage_pulsed.py [see --help]
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

Pulsed EML tracks for the sled passage.

Tests the design of 2026-10-10: the two steady EML tracks under the cable
carry the casing weight; a pulsed hold-down track on the casing ceiling
(above the cable) and a pulsed hold-up track on the casing floor cancel the
spacecraft's load while the sled is overhead; the four Halbach limiters are a
failure backstop only and must not engage in a normal launch.

Models:
  * Point model: one meter of casing (12,000 kg/m) and one meter of cable
    (119,592 kg/m), both free to move, with weight, the steady EML, the moving
    sled load, the pulsed tracks (section quantization, current rise lag,
    lead, timing and mass errors, gap dependence of the attractive force with
    plate saturation), the limiters and contact.  The 2026-10-09 run showed
    one meter of casing on its own is exact above about 5 km/s.
  * Beam model: the casing as a Timoshenko truss and the cable as a tensioned
    string on the same grid, with pulsed sections along the ring that can be
    made to fail, fire late or fire spuriously, and the parked sled.
  * Closed form: forces and widths (A), energy and power (E), plate heat (F),
    the push-down option (J), and the ring's radial equation driven by the
    launch profile of lsm_simulation.py (G).

Usage:
    python sled_passage_pulsed.py                 # everything, graphs on
    python sled_passage_pulsed.py --no-graphs
    python sled_passage_pulsed.py --l_sc=100      # one spacecraft load length (m)
    python sled_passage_pulsed.py --b_sat=1.6     # plate saturation (T) for the sizing (1.2 and 1.6 are always tabulated)
    python sled_passage_pulsed.py --tau=5         # baseline current rise time (ms)
    python sled_passage_pulsed.py --parts=ABK     # run a subset of parts A B C D E F G H I J K
    python sled_passage_pulsed.py --help

Conventions: sled speed u relative to the casing; inertial speed u + 516 m/s;
radial load from the casing frame (12.6 g at 30 km/s).  Forces per meter of
ring.  Up is positive for casing and cable motion.  The spacecraft zone is at
the front of the sled unless stated.

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import contextlib
import io
import math
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ring_altitude as ring
import lsm_config as lsm
import levitation_config as lc
import levitation_physics as lp
import sled_passage_simulation as sp

# =============================================================================
# PARAMETERS
# =============================================================================

MU = sp.MU                      # casing, kg/m
MU_C = sp.MU_CABLE              # cable with hardware, kg/m
W = sp.W_CASING                 # casing weight, N/m
G_NET = sp.G_NET
G0 = sp.G0_STD
GAP0 = sp.GAP0                  # 0.100 m
GAP_FORBIDDEN = sp.GAP_FORBIDDEN  # 0.020 m
BAND = 0.020                    # the casing must stay within 100 +/- 20 mm of the cable
BAND_STRICT = 0.010
L_SLED = sp.L_SLED
M_SLED_M = sp.M_SLED_M
M_SC0 = sp.M_SC                 # 5,000 t
V_REL = sp.V_REL
CABLE_FACE = ring.CABLE_SIDE    # 8.24 m, the cable face width available to a track
MU0 = lc.MU0
B_SAT_CASES = (1.2, 1.6)        # T
P_SAT = {b: b**2 / (2 * MU0) for b in B_SAT_CASES}   # Pa, pull per square meter of track at saturation
F_STEADY = W                    # steady tracks hold the casing weight
F_STEADY_CAP = sp.F_EML_OP      # 218 kN/m, two 1 m tracks at 400 A
L_SC_CASES = (50.0, 100.0, 200.0)
POSITIONS = ("front", "middle", "rear")
TAU_CASES = (1e-3, 2e-3, 5e-3, 10e-3)
L_SECTION_CASES = (10.0, 25.0, 50.0, 100.0)
DT_ERR_CASES = (0.0, 0.5e-3, -0.5e-3, 1e-3, -1e-3, 2e-3, -2e-3, 5e-3, -5e-3)
ERR_CASES = (0.01, -0.01, 0.05, -0.05)
LEAD_TAUS = 3.0                 # the command starts this many rise times early
PROT_V_WINDOW = 2e-3            # s, the gap protection estimates the closing speed over this window of sensor samples (at least the delay)
PROT_X_WINDOW = 50.0            # m, the beam model's protection acts on the offset averaged over this length (sensor and coil resolution)
FRINGE_FACTOR = 1.0             # gap field energy multiplier for fringing (1.0 = none assumed)
RECOVERY_EFF = 0.90             # fraction of a section's field energy recovered on release (assumed)
COIL_LOSS_FRACTION = 0.02       # resistive and AC loss per pulse as a fraction of stored energy (assumed)
PALETTE = sp.PALETTE
MATCH_RING = MU_C / (MU + MU_C)   # fraction of the sled impulse the pulse cancels so casing and cable end up moving together


def a_down(u):
    return sp.a_down(u)


def q_zone(u, L_sc, M_sc=M_SC0):
    """Sled load per meter under the spacecraft zone, N/m, positive down."""
    return (M_SLED_M + M_sc / L_sc) * a_down(u)


def q_bare(u):
    return M_SLED_M * a_down(u)


def hr(t=""):
    print("\n" + "=" * 78)
    if t:
        print(t)
        print("=" * 78)


# =============================================================================
# POINT MODEL
# =============================================================================

def lane(**kw):
    """A parameter set for the point model, with defaults."""
    d = dict(L_sc=100.0, pos="front", tau=2e-3, L_section=50.0, pos_in_section=0.5, dt_err=0.0, mass_err=0.0,
             speed_err=0.0, tau_err=0.0, steady="hold", steady_loop_hz=0.0, match="ring", hd_width=None, hu_width=None,
             B_sat=1.2, T_spread=None, w0=0.0, gap_dep=True, cable_free=True, push_share=0.0, push_k=0.0, M_sc=M_SC0,
             cancel=True, gap_err=0.0, label="",
             # 2026-10-11 layout options (defaults reproduce the earlier runs)
             n_coils=1, L_coil=None,          # staggered coils: n_coils overlapping coils of length L_coil, a new one every L_coil/n_coils
             no_sled=False,                   # no sled load at all (spurious pull, protection tests at rest)
             support_fail_t=None,             # time (s, sled front at 0) at which the casing support tracks drop
             support_fail_frac=1.0,           # the fraction of the support force lost then (0.5 = one track of the pair)
             prot=None,                       # gap protection on the compensation tracks: dict(trigger, delay, hz, cap_up, cap_down, sensor_err)
             false_pull=None,                 # (t0, duration, force): a pull commanded with no cause, force up positive
             pull_scale=None,                 # [(t0, t1, factor), ...]: the commanded pull is multiplied by factor between t0 and t1
             cmd_ramp=0.0,                    # s: the planned pull command is smoothed with a moving average of this width (a slow ramp by command)
             steady_cap=None)                 # N/m: ceiling of the casing support force for this lane (one track out), default F_STEADY_CAP
    d.update(kw)
    return d


def box_corrections(T, tau):
    """A coil whose current follows the command with a first-order lag of time constant tau, and whose
    force goes as the current squared, is commanded with a box.  Returns (d_a, d_b): shift the command's
    start and end by these amounts from the target box [a, a + T] so that the delivered force has the
    target impulse (level x T) and the same centroid.  Closed form for the lagged box response."""
    if T <= 0:
        return 0.0, 0.0

    def J_and_m(Tp):
        e1 = math.exp(-Tp / tau); e2 = math.exp(-2 * Tp / tau)
        Ib = 1 - e1
        J = Tp - 2 * tau * (1 - e1) + 0.5 * tau * (1 - e2) + Ib**2 * tau / 2
        m1 = Tp**2 / 2 - 2 * (tau**2 - tau * e1 * (Tp + tau)) + (tau**2 / 4 - 0.5 * tau * e2 * (Tp + tau / 2))
        m_tail = Ib**2 * (Tp * tau / 2 + tau**2 / 4)
        return J, (m1 + m_tail) / J          # impulse (unit level) and centroid measured from the command start

    lo, hi = T, T + 4 * tau
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if J_and_m(mid)[0] < T:
            lo = mid
        else:
            hi = mid
    Tp = 0.5 * (lo + hi)
    _, m = J_and_m(Tp)
    d_a = T / 2 - m
    d_b = Tp - T + d_a
    return d_a, d_b


def delivered_impulse(boxes, tau):
    """Impulse (in units of the largest level) delivered through the current lag by a piecewise-constant
    command made of boxes (a, b, level) of one sign, with the force going as the current squared.
    Where boxes overlap the larger level is commanded.  Closed form segment by segment."""
    F_ref = max(abs(b[2]) for b in boxes)
    pts = sorted(set([a for a, _, _ in boxes] + [b for _, b, _ in boxes]))
    I = 0.0
    J = 0.0
    for p0, p1 in zip(pts[:-1], pts[1:]):
        lev = sum(abs(l) for a, b, l in boxes if a <= p0 and b >= p1)     # one coil: overlapping boxes add
        c = math.sqrt(lev / F_ref)
        T = p1 - p0
        e1 = math.exp(-T / tau); e2 = math.exp(-2 * T / tau)
        J += c * c * T + 2 * c * (I - c) * tau * (1 - e1) + (I - c) ** 2 * (tau / 2) * (1 - e2)
        I = c + (I - c) * e1
    J += I * I * tau / 2
    return J * F_ref


def shape_boxes(boxes, tau):
    """Shift the command boxes so that the lagged, squared current delivers the target impulse with the
    right centroid for each box, then trim the end of the largest box to correct for the coupling between
    adjacent boxes (the current does not fall to zero between them)."""
    if not boxes:
        return boxes
    J_tgt = sum(abs(lev) * (b - a) for a, b, lev in boxes)
    shaped = []
    for a, b, lev in boxes:
        d_a, d_b = box_corrections(b - a, tau)
        shaped.append([a + d_a, b + d_b, lev])
    k = max(range(len(shaped)), key=lambda i: abs(shaped[i][2]))
    for _ in range(6):
        J = delivered_impulse([tuple(s) for s in shaped], tau)
        shaped[k][1] += (J_tgt - J) / abs(shaped[k][2])
    return [tuple(s) for s in shaped]


def point_pulsed(lanes, u, tabs, t_settle=1.5, dt=None):
    """Integrate casing and cable through one pass at speed u for every lane.
    Returns a list of result dicts (one per lane) and the time histories."""
    n = len(lanes)
    taus = [ln["tau"] for ln in lanes]
    Lmin = min(ln["L_sc"] for ln in lanes)
    if dt is None:
        dt = min(min(taus) / 20.0, Lmin / (24.0 * u), 0.25e-3)
    t_start = -(max(LEAD_TAUS * t for t in taus) + 10e-3 + 1.2)   # steady-track ramp for case (b) needs 1 s
    t_end = L_SLED / u + max(ln["T_spread"] or 0 for ln in lanes) + t_settle
    nt = int(math.ceil((t_end - t_start) / dt)) + 1
    t = t_start + np.arange(nt) * dt

    # per lane: sled load at the point, commanded pulsed force (up positive), steady force
    q = np.zeros((n, nt))
    F_cmd = np.zeros((n, nt))
    S = np.full((n, nt), F_STEADY)
    F_ref_hd = np.ones(n)
    F_ref_hu = np.ones(n)
    cap_hd = np.full(n, np.inf)
    cap_hu = np.full(n, np.inf)
    J_sled = np.zeros(n)
    for i, ln in enumerate(lanes):
        L_sc = ln["L_sc"]
        u_assumed = u * (1 + ln["speed_err"])
        d_pos = {"front": 0.0, "middle": 0.5 * (L_SLED - L_sc), "rear": L_SLED - L_sc}[ln["pos"]]
        tz0, tz1 = d_pos / u, (d_pos + L_sc) / u          # zone over the point
        cov_zone = np.clip((np.minimum(tz1, t + dt) - np.maximum(tz0, t)) / dt, 0.0, 1.0)      # fraction of each step under the zone
        cov_sled = np.clip((np.minimum(L_SLED / u, t + dt) - np.maximum(0.0, t)) / dt, 0.0, 1.0)
        qz = float(q_zone(u, L_sc, ln["M_sc"]))
        qb = float(q_bare(u))
        q[i] = qz * cov_zone + qb * (cov_sled - cov_zone)
        if ln["no_sled"]:
            q[i] = 0.0
        if ln["support_fail_t"] is not None:
            S[i, t >= ln["support_fail_t"]] *= (1.0 - ln["support_fail_frac"])
        J_sled[i] = qz * L_sc / u + qb * (L_SLED - L_sc) / u
        f_match = {"ring": MATCH_RING, "casing": 1.0}.get(ln["match"], ln["match"])
        if not ln["cable_free"] and ln["match"] == "ring":
            f_match = 1.0
        # what the feed-forward thinks the loads are (mass and speed errors)
        qz_ff = float(q_zone(u_assumed, L_sc, ln["M_sc"] * (1 + ln["mass_err"])))
        qb_ff = float(q_bare(u_assumed))
        if ln["steady"] == "reduced":
            # steady track set to W - |q_bare| (clipped to its range) from 1 s before the sled to 1 s after, ramped
            S_low = min(max(W - abs(qb_ff), 0.0), F_STEADY_CAP) if qb_ff < 0 else min(W + qb_ff, F_STEADY_CAP)
            ramp_in = np.clip((t + 1.0) / 1.0, 0, 1)
            ramp_out = np.clip((t - L_SLED / u) / 1.0, 0, 1)
            S[i] = W + (S_low - W) * np.where(t < L_SLED / u, ramp_in, 1 - ramp_out)
            qb_left = 0.0                                  # the bare sled is carried by the steady track
            qz_cmd = qz_ff - qb_ff                         # the pulse carries only the zone's excess
        else:
            qb_left = qb_ff
            qz_cmd = qz_ff
        if not ln["cancel"]:
            continue
        # target boxes (start, end, level); the feed-forward assumes a rise time tau_assumed and shapes each box for it
        boxes = []
        if ln["T_spread"]:
            # the pull window of length u T sweeps the ring; a section switches when the window reaches it, so a point
            # sees the window shifted by its place in the section (jitter of up to L_section / u)
            T = ln["T_spread"]
            tc = 0.5 * (tz0 + tz1) + (ln["pos_in_section"] - 0.5) * ln["L_section"] / u
            boxes.append((tc - T / 2, tc + T / 2, f_match * qz_cmd * (L_sc / u) / T))
            s0, s1 = 0.0, 0.0
            z0, z1 = tz0, tz1
        else:
            Ls, p = ln["L_section"], ln["pos_in_section"]
            s0, s1 = -p * Ls, (1 - p) * Ls
            t_on = (s0 + d_pos) / u
            t_off = (s1 + L_sc + d_pos) / u
            boxes.append((t_on, t_off, f_match * qz_cmd * L_sc / (L_sc + Ls)))
            z0, z1 = t_on, t_off
        if qb_left:
            tb0, tb1 = s0 / u, (s1 + L_SLED) / u
            if tb0 < z0:
                boxes.append((tb0, z0, f_match * qb_left))
            if z1 < tb1:
                boxes.append((z1, tb1, f_match * qb_left))
        tau_assumed = ln["tau"] / (1 + ln["tau_err"])
        shaped = shape_boxes(boxes, tau_assumed)
        if ln["n_coils"] > 1 or ln["L_coil"] is not None:
            # staggered coils: the point is covered by n_coils coils of length L_coil, starting every L_coil/n_coils behind it;
            # each coil is on while the pull window overlaps it and carries 1/n of the level, impulse matched over its own on-time.
            # n_coils = 1 with L_coil set is a section switched as a block: on for the window time plus L_coil / u.
            Lc, nc = ln["L_coil"] or L_SLED, ln["n_coils"]
            sub = []
            for (a, b, lev) in shaped:
                for k in range(nc):
                    d_k = (k + 0.5) * Lc / nc
                    sub.append((a - d_k / u, b + (Lc - d_k) / u, lev * (b - a) / ((b - a) + Lc / u) / nc))
            shaped = sub
        cmd = np.zeros(nt)
        for (a, b, lev) in shaped:    # area-exact edges; overlapping boxes add (one coil carries both loads)
            cov = np.clip((np.minimum(b + ln["dt_err"], t + dt) - np.maximum(a + ln["dt_err"], t)) / dt, 0.0, 1.0)
            cmd += lev * cov
        if ln["pull_scale"]:
            for (ps0, ps1, fac) in ln["pull_scale"]:
                cmd[(t >= ps0) & (t < ps1)] *= fac
        if ln["cmd_ramp"] and ln["cmd_ramp"] > dt:
            k_r = int(round(ln["cmd_ramp"] / dt))
            cmd = np.convolve(cmd, np.ones(k_r) / k_r, mode="same")       # impulse preserved, edges ramped over cmd_ramp
        F_cmd[i] = cmd
        lanes[i]["_lead"] = boxes[0][0] - shaped[0][0]
        F_ref_hd[i] = max(-cmd.min(), 1.0)      # hold-down reference (down)
        F_ref_hu[i] = max(cmd.max(), 1.0)       # hold-up reference (up)
        ps = P_SAT[ln["B_sat"]]
        wd_hd = ln["hd_width"] if ln["hd_width"] is not None else F_ref_hd[i] / ps
        wd_hu = ln["hu_width"] if ln["hu_width"] is not None else F_ref_hu[i] / ps
        cap_hd[i] = wd_hd * ps
        cap_hu[i] = wd_hu * ps
        lanes[i]["_hd_width"] = wd_hd
        lanes[i]["_hu_width"] = wd_hu

    # gap protection and false pulls use the compensation tracks at their track rating
    prot_on = np.array([ln["prot"] is not None for ln in lanes])
    prot_trig = np.array([(ln["prot"] or {}).get("trigger", 0.01) for ln in lanes])
    prot_delay = np.array([int(round((ln["prot"] or {}).get("delay", 0.0) / dt)) for ln in lanes])
    prot_om = 2 * math.pi * np.array([(ln["prot"] or {}).get("hz", 5.0) for ln in lanes])
    prot_err = np.array([(ln["prot"] or {}).get("sensor_err", 0.0) for ln in lanes])
    for i, ln in enumerate(lanes):
        if ln["prot"] is not None:
            cap_hu[i] = ln["prot"].get("cap_up", cap_hu[i])
            cap_hd[i] = ln["prot"].get("cap_down", cap_hd[i])
        if ln["prot"] is not None or ln["false_pull"] is not None:
            F_ref_hu[i] = max(F_ref_hu[i], cap_hu[i])
            F_ref_hd[i] = max(F_ref_hd[i], cap_hd[i])
    false_pull = [ln["false_pull"] for ln in lanes]
    engaged = np.zeros(n, dtype=bool)
    t_engaged = np.full(n, np.nan)
    k_win = max(int(round(PROT_V_WINDOW / dt)), int(prot_delay.max()), 1)   # velocity estimate over at least the sensor delay
    n_buf = int(prot_delay.max()) + k_win + 2
    rel_buf = np.zeros((n, n_buf))
    rel_meas_prev = np.zeros(n)

    tau_arr = np.array(taus)
    lag_step = 1.0 - np.exp(-dt / tau_arr)     # exact first-order lag for a piecewise-constant command
    loop_om = 2 * math.pi * np.array([ln["steady_loop_hz"] for ln in lanes])
    steady_cap = np.array([ln["steady_cap"] if ln["steady_cap"] is not None else F_STEADY_CAP for ln in lanes])
    S = np.minimum(S, steady_cap[:, None])
    t_pass_end = L_SLED / u + max(ln["T_spread"] or 0 for ln in lanes) + 0.1
    gap_dep = np.array([ln["gap_dep"] for ln in lanes])
    gap_err = np.array([ln["gap_err"] for ln in lanes])
    cable_free = np.array([ln["cable_free"] for ln in lanes])
    push_share = np.array([ln["push_share"] for ln in lanes])
    push_k = np.array([ln["push_k"] for ln in lanes])
    F0, k2 = tabs["F0"], 2 * tabs["k"]

    w = np.array([ln["w0"] for ln in lanes], dtype=float)   # casing position relative to its equilibrium
    v = np.zeros(n)
    wc = np.zeros(n)                                        # cable position
    vc = np.zeros(n)
    I_hd = np.zeros(n)
    I_hu = np.zeros(n)
    min_gb = np.full(n, GAP0)
    min_gt = np.full(n, GAP0)
    max_exc = np.zeros(n)
    max_exc_pass = np.zeros(n)
    pk_lim = np.zeros(n)
    pk_hd = np.zeros(n)
    pk_hu = np.zeros(n)
    J_hd = np.zeros(n)
    J_hu = np.zeros(n)
    contact = np.zeros(n, dtype=bool)
    hist = dict(t=t, w=np.zeros((n, nt)), wc=np.zeros((n, nt)), F=np.zeros((n, nt)))

    for j in range(nt):
        cmd = F_cmd[:, j].copy()
        if prot_on.any():
            # the sensors see the casing-to-cable offset with a delay and an error; once a face's gap has closed by the
            # trigger the compensation tracks run a position loop (critically damped at prot hz) plus the weight the
            # support tracks are not carrying, clipped to the track ratings, on whichever face can pull the right way
            rel_now = w - wc
            rel_buf = np.roll(rel_buf, 1, axis=1); rel_buf[:, 0] = rel_now
            rel_meas = rel_buf[np.arange(n), prot_delay] + prot_err
            v_meas = (rel_meas - rel_buf[np.arange(n), prot_delay + k_win] - prot_err) / (k_win * dt) if j >= k_win else np.zeros(n)
            newly = prot_on & (~engaged) & (np.abs(rel_meas) > prot_trig)
            t_engaged[newly] = t[j]
            engaged |= newly
            S_now = S[:, j]
            F_des = np.where(engaged, (W - S_now) - MU * (2 * prot_om * v_meas + prot_om**2 * rel_meas), 0.0)
            cmd = cmd + np.clip(F_des, -cap_hd, cap_hu)
        for i_fp, fp in enumerate(false_pull):
            if fp is not None and fp[0] <= t[j] < fp[0] + fp[1]:
                cmd[i_fp] += fp[2]
        # commanded currents (force goes as current squared), first-order lag with time constant tau
        I_cmd_hd = np.sqrt(np.clip(-cmd, 0, None) / F_ref_hd)
        I_cmd_hu = np.sqrt(np.clip(cmd, 0, None) / F_ref_hu)
        I_hd += (I_cmd_hd - I_hd) * lag_step
        I_hu += (I_cmd_hu - I_hu) * lag_step
        rel = w - wc
        gap_t = GAP0 + rel
        gap_b = GAP0 - rel
        # open loop: force goes as 1/gap^2; regulated: the current is set for the measured gap (true gap plus sensor error)
        fac_t = np.where(gap_dep, (GAP0 / np.maximum(gap_t, 1e-3)) ** 2, ((gap_t + gap_err) / np.maximum(gap_t, 1e-3)) ** 2)
        fac_b = np.where(gap_dep, (GAP0 / np.maximum(gap_b, 1e-3)) ** 2, ((gap_b + gap_err) / np.maximum(gap_b, 1e-3)) ** 2)
        F_hd = np.minimum(F_ref_hd * I_hd**2 * (1 - push_share) * fac_t, cap_hd)      # pulls casing down, cable up
        F_pd = F_ref_hd * I_hd**2 * push_share * np.exp(np.minimum(push_k * rel, 5.0))  # push-down, stable with the gap
        F_hu = np.minimum(F_ref_hu * I_hu**2 * fac_b, cap_hu)                          # pulls casing up, cable down
        F_lb = F0 * np.exp(-k2 * np.maximum(gap_b, 0.0))
        F_lt = F0 * np.exp(-k2 * np.maximum(gap_t, 0.0))
        Sj = S[:, j]
        if loop_om.any():
            Sj = np.clip(Sj - MU * (2 * loop_om * (v - vc) + loop_om**2 * rel), 0.0, steady_cap)
        a_cas = (-W + Sj - q[:, j] - F_hd - F_pd + F_hu - F_lb + F_lt) / MU
        a_cab = np.where(cable_free, (W - Sj + F_hd + F_pd - F_hu + F_lb - F_lt) / MU_C, 0.0)
        v += a_cas * dt
        w += v * dt
        vc += a_cab * dt
        wc += vc * dt
        rel = w - wc
        hb = rel >= GAP0
        if hb.any():
            w[hb] = wc[hb] + GAP0
            v[hb] = np.minimum(v[hb], vc[hb])
            contact |= hb
        ht = rel <= -GAP0
        if ht.any():
            w[ht] = wc[ht] - GAP0
            v[ht] = np.maximum(v[ht], vc[ht])
            contact |= ht
        rel = w - wc
        np.minimum(min_gb, GAP0 - rel, out=min_gb)
        np.minimum(min_gt, GAP0 + rel, out=min_gt)
        np.maximum(max_exc, np.abs(rel), out=max_exc)
        if t[j] <= t_pass_end:
            np.maximum(max_exc_pass, np.abs(rel), out=max_exc_pass)
        np.maximum(pk_lim, np.maximum(F_lb, F_lt), out=pk_lim)
        np.maximum(pk_hd, F_hd + F_pd, out=pk_hd)
        np.maximum(pk_hu, F_hu, out=pk_hu)
        J_hd += (F_hd + F_pd) * dt
        J_hu += F_hu * dt
        hist["w"][:, j] = w
        hist["wc"][:, j] = wc
        hist["F"][:, j] = F_hu - F_hd - F_pd

    out = []
    for i, ln in enumerate(lanes):
        out.append(dict(label=ln["label"], u=u, L_sc=ln["L_sc"], pos=ln["pos"], tau=ln["tau"], L_section=ln["L_section"],
                        dt_err=ln["dt_err"], mass_err=ln["mass_err"], speed_err=ln["speed_err"], steady=ln["steady"],
                        T_spread=ln["T_spread"], hd_width=ln.get("_hd_width", 0.0), hu_width=ln.get("_hu_width", 0.0),
                        min_gap_b=float(min_gb[i]), min_gap_t=float(min_gt[i]), max_exc=float(max_exc_pass[i]),
                        max_exc_all=float(max_exc[i]), lead=ln.get("_lead", 0.0),
                        v_rel_end=float(v[i] - vc[i]), v_cas_end=float(v[i]), v_cab_end=float(vc[i]),
                        t_20=(BAND / abs(v[i] - vc[i]) if abs(v[i] - vc[i]) > 1e-9 else float('inf')),
                        pk_lim=float(pk_lim[i]), pk_hd=float(pk_hd[i]), pk_hu=float(pk_hu[i]),
                        J_hd=float(J_hd[i]), J_hu=float(J_hu[i]), J_sled=float(J_sled[i]), contact=bool(contact[i]),
                        t_engaged=float(t_engaged[i]), engaged=bool(engaged[i]),
                        ok_band=bool(max_exc_pass[i] <= BAND and pk_lim[i] < 1e3), ok_strict=bool(max_exc_pass[i] <= BAND_STRICT),
                        dt=dt))
    return out, hist


def point_table(rows, title, cols=("u", "L_sc", "pos", "tau", "L_section", "dt_err", "max_exc", "min_gap_b", "min_gap_t",
                                   "v_rel_end", "pk_lim", "pk_hd", "hd_width", "ok")):
    print(f"\n{title}")
    hdr = {"u": "u km/s", "L_sc": "L_sc", "pos": "pos", "tau": "tau ms", "L_section": "sect m", "dt_err": "dt ms",
           "max_exc": "exc mm", "min_gap_b": "gap_b", "min_gap_t": "gap_t", "v_rel_end": "v_end", "pk_lim": "F_lim",
           "pk_hd": "F_hd", "hd_width": "w_hd m", "ok": "band", "label": "case", "T_spread": "T ms", "mass_err": "m err",
           "steady": "steady", "pk_hu": "F_hu", "hu_width": "w_hu m", "J_hd": "J_hd", "J_sled": "J_sled", "t_20": "t20 s",
           "max_exc_all": "exc all"}
    print("  " + " ".join(f"{hdr[c]:>9}" for c in cols))
    for r in rows:
        cells = []
        for c in cols:
            if c == "u":
                cells.append(f"{r['u']/1e3:9.1f}")
            elif c in ("tau", "dt_err", "T_spread"):
                val = r[c]
                cells.append(f"{(val or 0)*1e3:9.2f}")
            elif c in ("max_exc", "min_gap_b", "min_gap_t", "max_exc_all"):
                cells.append(f"{r[c]*1e3:9.1f}")
            elif c == "t_20":
                cells.append(f"{min(r[c], 999.0):9.2f}")
            elif c in ("pk_lim", "pk_hd", "pk_hu"):
                cells.append(f"{r[c]/1e3:9.1f}")
            elif c in ("J_hd", "J_sled"):
                cells.append(f"{r[c]/1e3:9.2f}")
            elif c == "ok":
                cells.append(f"{('ok' if r['ok_band'] else 'FAIL') + ('/10' if r['ok_strict'] else ''):>9}")
            elif c in ("hd_width", "hu_width", "v_rel_end", "mass_err"):
                cells.append(f"{r[c]:9.3f}")
            else:
                cells.append(f"{str(r[c])[:9]:>9}")
        print("  " + " ".join(cells))


# =============================================================================
# BEAM MODEL (casing truss + cable string) WITH PULSED SECTIONS
# =============================================================================

def beam_pulsed(u, L_sc, trusses, cases, tabs, L_section=50.0, tau=2e-3, L_win=6000.0, t_settle=1.0,
                match=MATCH_RING, dx=None, parked=False, F_steady_cap=None, loop_hz=1.0, T_spread=None, snap=False,
                n_coils=1, L_coil=None, track_width=None, prot=None):
    """2026-10-11 options: n_coils staggered coils of length L_coil (a new one every L_coil / n_coils) instead of
    sections; track_width (m per face) sets the force cap at 1.2 T; prot = dict(trigger, delay, hz, cap_up, cap_down,
    sensor_err) switches on the gap protection on every node (a latched position loop on the compensation tracks
    after the measured offset passes the trigger).  Case keys added: no_sled, support_fail_t, false_pull=(t0, dur, F).
    """
    if n_coils > 1:
        L_section = L_coil or L_section
    """Casing (Timoshenko truss lanes) and cable (tensioned string) on one grid.
    cases: list of dicts with keys
        fail_x: (x0, x1) sections whose supply has failed (no pulse) or None
        open_x: (x0, x1) sections that have lost their gap regulation (open loop) or None
        regulated: False for the whole pulse open loop
        late: schedule delay (s), negative = early
        spurious: True to fire the full pulse schedule with no sled
        no_pulse: True for no pulsed track at all
        no_loop: True to switch the steady track's slow loop off
    Lanes = trusses x cases.  Steady EML at W with a loop_hz gap loop unless parked, where the
    steady track is an ideal controller with ceiling F_steady_cap.
    T_spread: the zone's pull is spread over a window of length u T centered on the zone; a section
    is on while the window overlaps it, at the impulse-matched level.
    snap: also return the casing-to-cable offset along the window when the zone center passes its middle.
    parked: the zone load sits at the window center, ramped on over 1 s; u is ignored for motion."""
    n_s, n_c = len(trusses), len(cases)
    n = n_s * n_c
    if dx is None:
        dx = L_sc / 8.0
    N = int(round(L_win / dx)) + 1
    x = np.arange(N) * dx
    lo_e, hi_e = x - dx / 2, x + dx / 2
    i_lo, i_hi = int(0.4 * N), int(0.6 * N)
    region = slice(i_lo, i_hi + 1)
    EI = np.array([trusses[i % n_s]["EI"] for i in range(n)])[:, None]
    GA = np.array([trusses[i % n_s]["GA"] for i in range(n)])[:, None]
    J = np.array([trusses[i % n_s]["J"] for i in range(n)])[:, None]
    c_max = max(max(max(s["c_s"], s["c_b"]) for s in trusses), math.sqrt(sp.T_CABLE / MU_C))
    dt = min(0.5 * dx / c_max, tau / 10.0, (L_sc / u / 12.0) if not parked else 1e-3)
    if parked:
        t0, t_end = 0.0, 4.0
    else:
        t0 = -LEAD_TAUS * tau - 20e-3
        t_end = (x[i_hi] + L_SLED) / u + t_settle
    nt = int(math.ceil((t_end - t0) / dt)) + 1
    qz = float(q_zone(u, L_sc)) if not parked else float(q_zone(0.0, L_sc))
    qb = float(q_bare(u)) if not parked else 0.0
    # coil families: family k has coils of length L_section starting at k * L_section / n_coils + m * L_section
    n_sec = int(math.ceil(L_win / L_section)) + 1
    stagger = L_section / n_coils
    fam_off = np.arange(n_coils) * stagger
    sec_x0 = fam_off[:, None] + (np.arange(n_sec) - 1)[None, :] * L_section       # (n_coils, n_sec)
    sec_x1 = sec_x0 + L_section
    sec_of_node = np.clip(((x[None, :] - fam_off[:, None]) // L_section).astype(int) + 1, 0, n_sec - 1)   # (n_coils, N)
    # steady track at W with a 1 Hz loop; the sections cancel the zone (impulse matched over the time the zone overlaps a
    # section, shaped for the current lag; a reference section at x = 0, others shifted in time) and the bare sled at its
    # own level while any bare part of the sled is over the section
    if parked:
        level = 0.0
        ref_box = (0.0, 0.0, 0.0)
    elif T_spread:
        T_on = T_spread + L_section / u                       # window of length u T plus the section
        level = match * abs(qz) * (L_sc / u) / T_on
        ref_box = shape_boxes([(-T_spread / 2 + L_sc / (2 * u), T_on - T_spread / 2 + L_sc / (2 * u), level)], tau)[0]
    else:
        level = match * abs(qz) * L_sc / (L_sc + L_section)
        ref_box = shape_boxes([(0.0, (L_section + L_sc) / u, level)], tau)[0] if level != 0 else (0.0, 0.0, 0.0)
    level_bare = match * abs(qb) if not parked else 0.0
    t_sec0 = sec_x0 / u if not parked else sec_x0                  # the zone front reaches the section start
    t_bare0 = (sec_x0 + L_sc) / u if not parked else sec_x0        # the bare sled reaches the section start
    t_bare1 = (sec_x1 + L_SLED) / u if not parked else sec_x1      # the sled tail leaves the section end
    t_snap = (x[N // 2] + L_sc / 2) / u if not parked else 0.0     # zone center at the window middle
    # per case: coil command schedule
    fail_mask = np.ones((n_c, n_coils, n_sec))
    reg_mask = np.ones((n_c, n_coils, n_sec))
    late = np.zeros(n_c)
    spurious = np.zeros(n_c, dtype=bool)
    no_pulse = np.zeros(n_c, dtype=bool)
    loop_on = np.ones(n_c)
    t_fail = np.full(n_c, np.inf)
    false_pull = [None] * n_c
    for ic, c in enumerate(cases):
        if c.get("fail_x"):
            x0, x1 = c["fail_x"]
            fail_mask[ic][(sec_x1 > x0) & (sec_x0 < x1)] = 0.0
        if c.get("fail_coil") is not None:          # one coil dead: family k, the coil covering x
            k, xx = c["fail_coil"]
            fail_mask[ic, k, (sec_x1[k] > xx) & (sec_x0[k] <= xx)] = 0.0
        if c.get("open_x"):
            x0, x1 = c["open_x"]
            reg_mask[ic][(sec_x1 > x0) & (sec_x0 < x1)] = 0.0
        if not c.get("regulated", True):
            reg_mask[ic] = 0.0
        late[ic] = c.get("late", 0.0)
        spurious[ic] = c.get("spurious", False) or c.get("no_sled", False)
        no_pulse[ic] = c.get("no_pulse", False)
        loop_on[ic] = 0.0 if c.get("no_loop", False) else 1.0
        if c.get("support_fail_t") is not None:
            t_fail[ic] = c["support_fail_t"]
        false_pull[ic] = c.get("false_pull")
    # a node is regulated if the coils covering it are (fraction of regulated coils over the node)
    reg_node_c = np.mean(np.stack([reg_mask[:, k, :][:, sec_of_node[k]] for k in range(n_coils)], axis=0), axis=0)  # (n_c, N)
    reg_node = np.repeat(reg_node_c, n_s, axis=0)
    om_loop = 2 * math.pi * loop_hz * np.repeat(loop_on, n_s)[:, None]
    t_fail_lane = np.repeat(t_fail, n_s)[:, None]
    I_sec = np.zeros((n, n_coils, n_sec))
    lag_step = 1.0 - math.exp(-dt / tau)
    F_ref = abs(level) if level != 0 else 1.0
    ps = P_SAT[1.2]
    width = track_width if track_width is not None else F_ref / ps
    cap = width * ps
    # gap protection (per node, latched), and false pulls, through their own current lag
    prot_on = prot is not None
    if prot_on:
        p_trig = prot.get("trigger", 0.01); p_delay = int(round(prot.get("delay", 0.0) / dt)); p_om = 2 * math.pi * prot.get("hz", 5.0)
        p_cap_up = prot.get("cap_up", cap); p_cap_dn = prot.get("cap_down", cap); p_err = prot.get("sensor_err", 0.0)
    else:
        p_delay = 0
    n_buf = p_delay + 2
    k_win = max(int(round(PROT_V_WINDOW / dt)), p_delay, 1)       # the closing speed is estimated over at least the sensor delay
    n_buf = n_buf + k_win
    m_avg = max(int(round(PROT_X_WINDOW / (2 * dx))), 0)          # the measured offset is averaged over +- m_avg nodes (about 50 m)
    ker = np.ones(2 * m_avg + 1) / (2 * m_avg + 1)
    rel_buf = np.zeros((n_buf, n, N))
    engaged = np.zeros((n, N), dtype=bool)
    I_prot = np.zeros((n, N))
    F_ref_prot = max(cap, (prot or {}).get("cap_up", 0.0), (prot or {}).get("cap_down", 0.0), 1.0)
    t_eng = np.full(n, np.nan)
    rel_meas_prev = np.zeros((n, N))

    T_c = sp.T_CABLE
    snapshot = None
    w = np.zeros((n, N)); v = np.zeros((n, N))
    phi = np.zeros((n, N - 1)); phid = np.zeros((n, N - 1))
    wc = np.zeros((n, N)); vcb = np.zeros((n, N))
    Vp = np.zeros((n, N + 1)); Mp = np.zeros((n, N))
    F0, k2 = tabs["F0"], 2 * tabs["k"]
    om = 2 * math.pi * sp.PD_FREQ_HZ
    min_gb = np.full(n, GAP0); min_gt = np.full(n, GAP0); max_exc = np.zeros(n)
    pk_lim = np.zeros(n); pk_steady = np.zeros(n); contact = np.zeros(n, dtype=bool)
    v_pk = np.zeros(n)
    profile = None
    for j in range(nt):
        t = t0 + j * dt
        if parked:
            xf = x[N // 2] + L_sc / 2
            ramp = min(max(t / 1.0, 0.0), 1.0)
            fz = np.clip((np.minimum(hi_e, xf) - np.maximum(lo_e, xf - L_sc)) / dx, 0, 1) * ramp
            fb = np.zeros(N)
        else:
            xf = u * t
            fz = np.clip((np.minimum(hi_e, xf) - np.maximum(lo_e, xf - L_sc)) / dx, 0, 1)
            fb = np.clip((np.minimum(hi_e, xf - L_sc) - np.maximum(lo_e, xf - L_SLED)) / dx, 0, 1)
        q = qz * fz + qb * fb
        # coil commands: each coil is on from t_sec0 + a' to t_sec0 + b' (shaped box), per case and family
        F_node = np.zeros((n, N))
        if not parked:
            for ic in range(n_c):
                if no_pulse[ic]:
                    continue
                tt = t - late[ic]
                on = ((tt >= t_sec0 + ref_box[0]) & (tt < t_sec0 + ref_box[1])).astype(float) * (level / n_coils)   # (n_coils, n_sec)
                bare = ((tt >= t_bare0) & (tt < t_bare1)).astype(float) * (level_bare / n_coils)
                I_cmd = np.sqrt((on + bare) / F_ref) * fail_mask[ic]      # one coil carries both loads
                sl = slice(ic * n_s, (ic + 1) * n_s)
                I_sec[sl] += (I_cmd[None, :, :] - I_sec[sl]) * lag_step
                for k in range(n_coils):
                    F_node[sl] += F_ref * I_sec[sl, k, :][:, sec_of_node[k]] ** 2
        rel = w - wc
        gap_t = GAP0 + rel
        gap_b = GAP0 - rel
        sign_down = 1.0 if qz < 0 else -1.0          # above the neutral speed the pulse pulls down
        if sign_down > 0:
            fac = reg_node + (1.0 - reg_node) * (GAP0 / np.maximum(gap_t, 1e-3)) ** 2
            F_p = np.minimum(F_node * fac, cap)      # lift compensation, attractive to the top plate
            F_on_casing = -F_p
        else:
            fac = reg_node + (1.0 - reg_node) * (GAP0 / np.maximum(gap_b, 1e-3)) ** 2
            F_p = np.minimum(F_node * fac, cap)      # load compensation
            F_on_casing = F_p
        # gap protection and false pulls (positive = floor tracks pull up)
        S_fail = (t >= t_fail_lane)                  # (n, 1) lanes whose support tracks are gone
        if prot_on or any(fp is not None for fp in false_pull):
            F_extra = np.zeros((n, N))
            if prot_on:
                rel_buf[j % n_buf] = rel
                rel_d = rel_buf[(j - p_delay) % n_buf]
                rel_dw = rel_buf[(j - p_delay - k_win) % n_buf]
                if m_avg > 0:
                    rel_d = np.array([np.convolve(np.pad(r, m_avg, mode="edge"), ker, mode="valid") for r in rel_d])
                    rel_dw = np.array([np.convolve(np.pad(r, m_avg, mode="edge"), ker, mode="valid") for r in rel_dw])
                rel_meas = rel_d + p_err
                v_meas = (rel_d - rel_dw) / (k_win * dt) if j >= p_delay + k_win else np.zeros((n, N))
                rel_meas_prev = rel_meas.copy()
                newly = (~engaged) & (np.abs(rel_meas) > p_trig)
                for il in range(n):
                    if np.isnan(t_eng[il]) and newly[il].any():
                        t_eng[il] = t
                engaged |= newly
                W_missing = np.where(S_fail, W, 0.0)
                F_des = np.where(engaged, W_missing - MU * (2 * p_om * v_meas + p_om**2 * rel_meas), 0.0)
                F_extra += np.clip(F_des, -p_cap_dn, p_cap_up)
            for ic, fp in enumerate(false_pull):
                if fp is not None and fp[0] <= t < fp[0] + fp[1]:
                    F_extra[ic * n_s:(ic + 1) * n_s] += fp[2]
            I_cmd_p = np.sign(F_extra) * np.sqrt(np.abs(F_extra) / F_ref_prot)
            I_prot += (I_cmd_p - I_prot) * lag_step
            F_on_casing = F_on_casing + np.sign(I_prot) * F_ref_prot * I_prot**2
        if snap and (snapshot is None) and t >= t_snap:
            snapshot = dict(x=x - x[N // 2], rel=(w - wc).copy(), F=F_on_casing.copy())
        # beam internal
        V = GA * ((w[:, 1:] - w[:, :-1]) / dx - phi)
        Vp[:, 1:-1] = V
        f_beam = (Vp[:, 1:] - Vp[:, :-1]) / dx
        Mp[:, 1:-1] = EI * (phi[:, 1:] - phi[:, :-1]) / dx
        a_phi = ((Mp[:, 1:] - Mp[:, :-1]) / dx + V) / J
        # cable string
        f_cab = np.zeros((n, N))
        f_cab[:, 1:-1] = T_c * (wc[:, 2:] - 2 * wc[:, 1:-1] + wc[:, :-2]) / dx**2
        F_lb = F0 * np.exp(-k2 * np.maximum(gap_b, 0.0))
        F_lt = F0 * np.exp(-k2 * np.maximum(gap_t, 0.0))
        if parked:
            S = np.clip(W + q[None, :] - MU * (2 * om * (v - vcb) + om * om * rel), 0.0, F_steady_cap)
        else:
            S = np.clip(W - MU * (2 * om_loop * (v - vcb) + om_loop * om_loop * rel), 0.0, F_STEADY_CAP)
            S = np.where(S_fail, 0.0, S)
        q_cas = q[None, :] * np.ones((n, 1))
        if spurious.any():
            for ic in range(n_c):
                if spurious[ic]:
                    q_cas[ic * n_s:(ic + 1) * n_s] = 0.0
        a = (f_beam - W + S - q_cas + F_on_casing - F_lb + F_lt) / MU
        a_c = (f_cab + W - S - F_on_casing + F_lb - F_lt) / MU_C
        v += a * dt; w += v * dt
        phid += a_phi * dt; phi += phid * dt
        vcb += a_c * dt; wc += vcb * dt
        rel = w - wc
        hb = rel >= GAP0
        if hb.any():
            w[hb] = wc[hb] + GAP0; v[hb] = np.minimum(v[hb], vcb[hb]); contact |= hb[:, region].any(axis=1)
        ht = rel <= -GAP0
        if ht.any():
            w[ht] = wc[ht] - GAP0; v[ht] = np.maximum(v[ht], vcb[ht]); contact |= ht[:, region].any(axis=1)
        rel = (w - wc)[:, region]
        np.minimum(min_gb, GAP0 - rel.max(axis=1), out=min_gb)
        np.minimum(min_gt, GAP0 + rel.min(axis=1), out=min_gt)
        np.maximum(max_exc, np.abs(rel).max(axis=1), out=max_exc)
        np.maximum(pk_lim, np.maximum(F_lb[:, region], F_lt[:, region]).max(axis=1), out=pk_lim)
        np.maximum(pk_steady, S[:, region].max(axis=1), out=pk_steady)
        np.maximum(v_pk, np.abs(v - vcb)[:, region].max(axis=1), out=v_pk)
    if parked:
        profile = dict(x=x - x[N // 2], S=S.copy(), rel=(w - wc).copy())
    elif snap:
        profile = snapshot
    out = []
    for il in range(n):
        out.append(dict(truss=trusses[il % n_s]["label"], case=cases[il // n_s].get("label", ""),
                        min_gap_b=float(min_gb[il]), min_gap_t=float(min_gt[il]), max_exc=float(max_exc[il]),
                        pk_lim=float(pk_lim[il]), pk_steady=float(pk_steady[il]), v_pk=float(v_pk[il]),
                        contact=bool(contact[il]), width=width, dt=dt, nt=nt, t_engaged=float(t_eng[il])))
    return out, profile


# =============================================================================
# PART G: RING RADIAL EQUATION
# =============================================================================

def lsm_profile():
    """Run lsm_simulation in quick mode, silently, and return t, v (relative to the casing), P_mech."""
    import lsm_simulation as ls
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ls.run_lsm_simulation(quick_mode=True)
    return np.array(ls.data_t), np.array(ls.data_v), np.array(ls.data_P_mech)


def ring_radial(t_prof, v_prof, with_anchors):
    """Uniform radial (breathing) mode of the whole ring driven by the sled's net outward force,
    which is spread around the ring by its motion.  Returns t, dr, and the mode parameters.

      mu_r dr'' = -k_eff dr + F_total(t) / L_ring
      k_eff = (E A - T) / R^2  [hoop]  + mu_c (3 V_c^2/R^2 - 2 g/R) + mu_l (3 V_gs^2/R^2 - 2 g/R)  [centrifugal and gravity gradient]
              + k_anchor  [anchor lines, uniform-area estimate]
    """
    R = ring.R_ORBIT
    mu_r = MU + MU_C
    EA = sp.E_CNT * ring.CABLE_AREA
    k_hoop = (EA - sp.T_CABLE) / R**2
    k_grav = (MU_C * (3 * ring.V_CABLE**2 / R**2 - 2 * sp.G_LOCAL / R) + MU * (3 * sp.V_GS**2 / R**2 - 2 * sp.G_LOCAL / R))
    spacing = 50e3
    n_lines = ring.L_RING / spacing
    A_line = MU * G_NET * spacing / ring.SIGMA_CABLE
    k_line = sp.E_CNT * A_line / ring.ALTITUDE
    k_anchor = k_line / spacing if with_anchors else 0.0
    k_eff = k_hoop + k_grav + k_anchor
    om = math.sqrt(k_eff / mu_r)
    v_i = v_prof + sp.V_GS
    F_tot = (M_SC0 + M_SLED_M * L_SLED) * (v_i**2 / R - sp.G_LOCAL)    # outward positive
    f = F_tot / ring.L_RING
    dt = 1.0
    nt = int(t_prof[-1] / dt) + 1
    tt = np.arange(nt) * dt
    ff = np.interp(tt, t_prof, f)
    dr = np.zeros(nt); vr = 0.0
    for j in range(1, nt):
        a = (-k_eff * dr[j - 1] + ff[j - 1]) / mu_r
        vr += a * dt
        dr[j] = dr[j - 1] + vr * dt
    return dict(t=tt, dr=dr, f=ff, F_tot=np.interp(tt, t_prof, F_tot), k_hoop=k_hoop, k_grav=k_grav, k_anchor=k_anchor,
                k_eff=k_eff, om=om, period=2 * math.pi / om, A_line=A_line, k_line=k_line, n_lines=n_lines,
                dr_static_end=f[-1] / k_eff)


# =============================================================================
# PART J: PUSH-DOWN OPTION (closed form)
# =============================================================================

def push_down(v_slip, wavelength, B_surface, gap, material="molybdenum", plate_thickness=lc.SAFETY_PLATE_THICKNESS):
    mat = lp.get_material(material)
    rho = mat["resistivity"]
    k = 2 * math.pi / wavelength
    om = k * v_slip
    delta = math.sqrt(2 * rho / (MU0 * om))
    kd = k * delta
    B_p = B_surface * math.exp(-k * gap)
    thick = plate_thickness >= 3 * delta
    P = B_p**2 / (2 * MU0) / (1 + kd**2)           # image-current limit for k delta << 1, falling beyond  [UNVERIFIED form]
    drag_lift = kd
    heat_per_N = v_slip * kd                        # W per N of push
    f_drive = (V_REL - v_slip) / wavelength
    K = math.sqrt(2) * B_surface / MU0              # Halbach sheet current for B_surface
    NI_pole = K * wavelength / 2
    E_field = B_surface**2 / (2 * MU0) / k          # J per m^2 of track in the field above the coils
    Q_reactive = 2 * math.pi * abs(f_drive) * E_field   # VA per m^2 of track
    return dict(delta=delta, kd=kd, B_p=B_p, P=P, drag_lift=drag_lift, heat_per_N=heat_per_N, f_drive=f_drive,
                K=K, NI_pole=NI_pole, E_field=E_field, Q_reactive=Q_reactive, thick=thick)


# =============================================================================
# GRAPH HELPERS
# =============================================================================

def _save(fig, name, out_dir):
    sp._save(fig, name, out_dir)


# =============================================================================
# MAIN
# =============================================================================

def main():
    L_scs = list(L_SC_CASES)
    b_sat = 1.2
    tau0 = 2e-3
    graphs = True
    parts = "ABCDEFGHIJK"
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--l_sc="):
            L_scs = [float(arg.split("=")[1])]
        elif arg.startswith("--b_sat="):
            b_sat = float(arg.split("=")[1])
        elif arg.startswith("--tau="):
            tau0 = float(arg.split("=")[1]) * 1e-3
        elif arg.startswith("--parts="):
            parts = arg.split("=")[1].upper()
        elif arg == "--no-graphs":
            graphs = False
        else:
            print(f"Unknown option {arg}"); return
    if b_sat not in P_SAT:
        P_SAT[b_sat] = b_sat**2 / (2 * MU0)
    out_dir = ring.out_dir("graphs_sled_passage")   # under ORBITAL_RING_OUT (default runs)
    t_all = time.time()
    tabs = sp.limiter_tables()
    ps = P_SAT[b_sat]
    L_main = 100.0 if 100.0 in L_scs else L_scs[0]

    hr("PULSED EML TRACKS FOR THE SLED PASSAGE")
    print(f"  Design point {ring.ALTITUDE_KM:.0f} km.  Casing {MU:,.0f} kg/m, weight {W/1e3:.1f} kN/m, cable {MU_C:,.0f} kg/m, "
          f"relative speed {V_REL:,.0f} m/s, cable face {CABLE_FACE:.2f} m.")
    print(f"  Spacecraft {M_SC0/1e3:,.0f} t on {L_scs} m; sled {L_SLED:,.0f} m at {M_SLED_M:.0f} kg/m.")
    print(f"  Steady EML held at {F_STEADY/1e3:.1f} kN/m (case a) unless stated; pulsed tracks sized at B_sat = {b_sat} T "
          f"({ps/1e3:.0f} kPa).  Baseline rise time {tau0*1e3:.0f} ms, lead {LEAD_TAUS:.0f} rise times, sections 50 m, point mid-section.")
    print(f"  Pulse target: cancel {MATCH_RING*100:.1f} % of the sled impulse (the cable's mass share), so casing and cable end "
          f"up moving outward together; 'casing' target cancels 100 % and leaves the cable rising relative to the casing.")
    print(f"  Limiter law from levitation_physics: {tabs['F0']/1e3:.1f} kN/m at contact, k = {tabs['k']:.1f} /m.")

    # ------------------------------------------------------------------ A
    if "A" in parts:
        hr("PART A: FORCE AND WIDTH NEEDED OVER THE LAUNCH (closed form)")
        print(f"  Pull per square meter: {P_SAT[1.2]/1e3:.0f} kPa at 1.2 T, {P_SAT[1.6]/1e3:.0f} kPa at 1.6 T.")
        print(f"  Hold-down force with the steady track at 95 kN/m = |q_zone|; with it reduced to zero = |q_zone| - 95; with it")
        print(f"  reduced to 95 - |q_bare| (case b as run) = |q_zone| - |q_bare|.  Hold-up force below the neutral speed = q_zone.")
        us = [0, 1000, 2000, 3000, 5000, 6988.0, 8000, 10000, 15000, 20000, 25000, 30000]
        for L in L_scs:
            print(f"\n  L_sc = {L:g} m, 5,000 t.  Forces in kN/m, widths in m (total track width, at 1.2 T and 1.6 T).")
            print(f"  {'u km/s':>7} {'q_zone':>8} {'q_bare':>7} {'hold-up':>8} {'w1.2':>6} {'w1.6':>6} {'hd a':>8} {'w1.2':>6} {'w1.6':>6} "
                  f"{'hd b0':>8} {'w1.2':>6} {'hd b':>8} {'w1.2':>6} {'>face':>5}")
            for u in us:
                qz, qb = float(q_zone(u, L)), float(q_bare(u))
                hu = max(qz, 0.0)
                hd_a = max(-qz, 0.0)
                hd_b0 = max(-qz - W, 0.0)
                hd_b = max(-qz + qb, 0.0) if qz < 0 else 0.0
                flag = "yes" if hd_a / P_SAT[1.2] > CABLE_FACE else ""
                print(f"  {u/1e3:7.2f} {qz/1e3:8.1f} {qb/1e3:7.1f} {hu/1e3:8.1f} {hu/P_SAT[1.2]:6.2f} {hu/P_SAT[1.6]:6.2f} "
                      f"{hd_a/1e3:8.1f} {hd_a/P_SAT[1.2]:6.2f} {hd_a/P_SAT[1.6]:6.2f} {hd_b0/1e3:8.1f} {hd_b0/P_SAT[1.2]:6.2f} "
                      f"{hd_b/1e3:8.1f} {hd_b/P_SAT[1.2]:6.2f} {flag:>5}")
        print(f"\n  The job passes from the hold-up track to the hold-down track at the circular velocity, u = {sp.V_ORBIT - sp.V_GS:,.0f} m/s,")
        print(f"  where the sled load is zero.  Width flag '>face': the case (a) hold-down width exceeds the {CABLE_FACE:.2f} m cable face.")
        print(f"  Largest spacecraft mass (t) a fixed hold-down width cancels in full (case a, L_sc = {L_main:g} m, 1.2 T):")
        print(f"  {'u km/s':>7} " + " ".join(f"{'%g m' % wd:>8}" for wd in (2, 4, 6, 8)) + f"   {'8.24 m face':>12}")
        for u in (10000, 15000, 20000, 25000, 30000):
            ad = -float(a_down(u))
            line = f"  {u/1e3:7.0f} "
            for wd in (2, 4, 6, 8, CABLE_FACE):
                m = L_main * (wd * P_SAT[1.2] / ad - M_SLED_M)
                line += f"{m/1e3:8,.0f} " if wd != CABLE_FACE else f"  {m/1e3:12,.0f}"
            print(line)
        qb30 = float(q_bare(30000.0))
        print(f"\n  Bare sled: {-qb30/1e3:.1f} kN/m upward at 30 km/s, under the casing weight, but with the steady track held at 95 kN/m")
        print(f"  the casing's net force under it is {-qb30/1e3:.1f} kN/m up for {L_SLED/30000:.2f} s, which gives {-qb30*L_SLED/30000/MU:.2f} m/s.")
        print(f"  So in case (a) the hold-down track must also cover the bare sled, at {-qb30/1e3:.1f} kN/m over its whole 10 km")
        print(f"  ({-qb30/P_SAT[1.2]:.2f} m of track at 1.2 T), or the steady track must be reduced by that much ahead of the sled (case b).")

    # ------------------------------------------------------------------ B
    rows_B = []
    sweep_rows = []
    if "B" in parts:
        hr("PART B: TIME DOMAIN WITH REAL PULSES (point model, casing and cable free)")
        print(f"  exc = largest |casing - cable| offset from the 100 mm gap; gap_b, gap_t = smallest gaps; v_end = casing speed relative")
        print(f"  to the cable after the pass; F_lim = peak limiter force; F_hd = peak hold-down (or hold-up) force; w_hd = its track width")
        print(f"  at {b_sat} T; band = ok if exc <= 20 mm and the limiter force stays under 1 kN/m, '/10' if also within 10 mm.")
        # baseline sweep over speed
        base = dict(tau=tau0, L_section=50.0, B_sat=b_sat)
        t0 = time.time()
        for L in L_scs:
            for u in [1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 8000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0]:
                lanes = [lane(L_sc=L, label="open", **base),
                         lane(L_sc=L, label="regulated", gap_dep=False, **base),
                         lane(L_sc=L, label="open+loop", steady_loop_hz=1.0, **base),
                         lane(L_sc=L, label="reg+loop", steady_loop_hz=1.0, gap_dep=False, **base),
                         lane(L_sc=L, label="reg+loop casing tgt", steady_loop_hz=1.0, gap_dep=False, match="casing", **base),
                         lane(L_sc=L, label="b reduced ahead", steady="reduced", gap_dep=False, **base),
                         lane(L_sc=L, label="no pulse", cancel=False, **base),
                         lane(L_sc=L, label="no pulse+loop", cancel=False, steady_loop_hz=1.0, **base)]
                res, _ = point_pulsed(lanes, u, tabs)
                sweep_rows.extend(res)
        print(f"  (sweep {time.time()-t0:.0f} s)")
        print(f"  Cases.  Steady track at 95 kN/m and the pulsed track cancels zone and bare sled (case a) unless stated.  open = pulse")
        print(f"  current set in advance, force varies with the gap (attractive instability); regulated = pulse current corrected for")
        print(f"  the measured gap so the force is the commanded one (a local loop on the pulsed coil, not a casing position loop).")
        print(f"  +loop = the steady track also has a slow gap loop, critically damped at 1 Hz, clipped to 0 to 218 kN/m.")
        print(f"  'b reduced ahead' = steady track reduced by the bare sled's load from 1 s before the sled (case b as written):")
        print(f"  the casing falls while the track is reduced and the sled is not yet there, so this case cannot work and is kept")
        print(f"  only to show it.  'casing tgt' cancels 100 % of the impulse instead of the cable's share.")
        for L in L_scs:
            point_table([r for r in sweep_rows if r["L_sc"] == L],
                        f"Speed sweep, L_sc = {L:g} m, tau {tau0*1e3:.0f} ms, 50 m sections, no timing error",
                        cols=("u", "label", "max_exc", "min_gap_b", "min_gap_t", "v_rel_end", "t_20", "pk_lim", "pk_hd", "hd_width", "pk_hu", "hu_width", "ok"))
        print(f"  t20 = time for the residual relative speed to drift 20 mm, what the slow steady loop has to correct it.")
        # tolerances at 30 km/s, case (b), pulse regulated against the gap and open loop
        t0 = time.time()
        REG_CASES = (("regulated", False), ("open loop", True))   # both with a 1 Hz loop on the steady track
        for reg_label, gd in REG_CASES:
            print(f"\n  Tolerance scans at 30 km/s (ring target, steady track with a 1 Hz loop, pulse {reg_label}):")
            for L in L_scs:
                lanes = []
                kw = dict(L_sc=L, B_sat=b_sat, steady="hold", steady_loop_hz=1.0, gap_dep=gd)
                for tau in TAU_CASES:
                    for de in DT_ERR_CASES:
                        lanes.append(lane(tau=tau, dt_err=de, L_section=50.0, label=f"tau{tau*1e3:g} dt{de*1e3:+g}", **kw))
                for Ls in L_SECTION_CASES:
                    for pis in (0.0, 0.5, 1.0):
                        lanes.append(lane(tau=tau0, L_section=Ls, pos_in_section=pis, label=f"sect{Ls:g} p{pis:g}", **kw))
                for e in ERR_CASES:
                    lanes.append(lane(tau=tau0, mass_err=e, label=f"mass{e*100:+g}%", **kw))
                    lanes.append(lane(tau=tau0, speed_err=e, label=f"speed{e*100:+g}%", **kw))
                    lanes.append(lane(tau=tau0, tau_err=4 * e, label=f"rise{e*400:+g}%", **kw))
                for pos in POSITIONS:
                    lanes.append(lane(tau=tau0, pos=pos, label=f"pos {pos}", **kw))
                res, _ = point_pulsed(lanes, 30000.0, tabs)
                rows_B.extend([dict(r, reg=reg_label) for r in res])
                print(f"\n  L_sc = {L:g} m: rise time and timing error (exc mm; '*' = outside 20 mm; '+' = outside 10 mm):")
                print(f"  {'tau ms':>7} " + " ".join(f"{'dt%+g' % (de*1e3):>8}" for de in DT_ERR_CASES))
                for tau in TAU_CASES:
                    line = f"  {tau*1e3:7.0f} "
                    for de in DT_ERR_CASES:
                        r = next(r for r in res if r["label"] == f"tau{tau*1e3:g} dt{de*1e3:+g}")
                        mark = "*" if not r["ok_band"] else ("+" if not r["ok_strict"] else " ")
                        line += f"{r['max_exc']*1e3:7.1f}{mark} "
                    print(line)
                print(f"  Section length and point position in the section (exc mm):")
                print(f"  {'sect m':>7} {'start':>8} {'middle':>8} {'end':>8}")
                for Ls in L_SECTION_CASES:
                    vals = [next(r for r in res if r["label"] == f"sect{Ls:g} p{pis:g}")["max_exc"] * 1e3 for pis in (0.0, 0.5, 1.0)]
                    print(f"  {Ls:7.0f} " + " ".join(f"{vv:8.1f}" for vv in vals))
                print(f"  Mass, speed and rise-time error, spacecraft position (exc mm, v_end m/s, t20 s):")
                for r in res:
                    if r["label"].startswith(("mass", "speed", "rise", "pos")):
                        print(f"    {r['label']:12} exc {r['max_exc']*1e3:6.1f} mm, v_end {r['v_rel_end']:+.3f} m/s, t20 {min(r['t_20'], 999):6.1f} s, "
                              f"gap_b {r['min_gap_b']*1e3:5.1f}, gap_t {r['min_gap_t']*1e3:5.1f}, F_lim {r['pk_lim']/1e3:.2f} kN/m")
            print(f"\n  Section length at low speed (hold-up track, pulse {reg_label}), L_sc = {L_main:g} m, exc mm:")
            print(f"  {'sect m':>7} " + " ".join(f"{'%g km/s' % (u/1e3):>9}" for u in (1000.0, 2000.0, 3000.0, 5000.0)))
            for Ls in L_SECTION_CASES:
                line = f"  {Ls:7.0f} "
                for u in (1000.0, 2000.0, 3000.0, 5000.0):
                    res, _ = point_pulsed([lane(L_sc=L_main, tau=tau0, L_section=Ls, B_sat=b_sat, steady="hold", steady_loop_hz=1.0, gap_dep=gd)], u, tabs)
                    line += f"{res[0]['max_exc']*1e3:8.1f}{'*' if not res[0]['ok_band'] else ' '}"
                print(line)
        print(f"  (tolerance scans {time.time()-t0:.0f} s)")
        # tolerance over the whole speed range: tau x dt_err at 5, 10, 20, 30 km/s
        for reg_label, gd in REG_CASES:
            print(f"\n  Timing tolerance over the launch (largest |dt_err| that keeps exc <= 20 mm at 5, 10, 20 and 30 km/s), 1 Hz steady loop, pulse {reg_label}, L_sc = {L_main:g} m:")
            for tau in TAU_CASES:
                ok_de = []
                for de_abs in (0.5e-3, 1e-3, 2e-3, 5e-3):
                    good = True
                    for u in (5000.0, 10000.0, 20000.0, 30000.0):
                        lanes = [lane(L_sc=L_main, tau=tau, dt_err=s_ * de_abs, B_sat=b_sat, steady="hold", steady_loop_hz=1.0, gap_dep=gd) for s_ in (1, -1)]
                        res, _ = point_pulsed(lanes, u, tabs)
                        if not all(r["ok_band"] for r in res):
                            good = False; break
                    if good:
                        ok_de.append(de_abs)
                print(f"    tau {tau*1e3:3.0f} ms: {'none' if not ok_de else f'+/- {max(ok_de)*1e3:g} ms'}")

    # ------------------------------------------------------------------ C
    rows_C = []
    if "C" in parts:
        hr("PART C: PARTIAL CANCELLING, smallest hold-down width at 30 km/s (point model)")
        print(f"  The hold-down pulse is matched to the zone in time; its force is capped by the width at 1.2 T, so a narrow track")
        print(f"  cancels only part of the kick.  Steady track at 95 kN/m with no loop (nothing restores a residual speed but the")
        print(f"  limiters) and with the 1 Hz loop.  Pulse current regulated against the gap.")
        widths = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 8.24, 10, 12, 16, 20]
        for L in L_scs:
            for pos in POSITIONS:
                for loop_label, lhz in (("no loop", 0.0), ("1 Hz loop", 1.0)):
                    lanes = [lane(L_sc=L, pos=pos, steady_loop_hz=lhz, hd_width=wd, hu_width=wd, tau=tau0, B_sat=1.2, gap_dep=False,
                                  label=f"w{wd:g}") for wd in widths]
                    res, _ = point_pulsed(lanes, 30000.0, tabs)
                    rows_C.extend([dict(r, reg=loop_label) for r in res])
                    need_w = next((wd for wd, r in zip(widths, res) if r["ok_band"]), None)
                    need_w10 = next((wd for wd, r in zip(widths, res) if r["ok_strict"]), None)
                    nc = next((wd for wd, r in zip(widths, res) if not r["contact"]), None)
                    print(f"  L_sc {L:4g} m, {pos:6}, {loop_label:9}: width within 20 mm {('none' if need_w is None else f'{need_w:g} m'):>7}, "
                          f"within 10 mm {('none' if need_w10 is None else f'{need_w10:g} m'):>7}, no contact {('none' if nc is None else f'{nc:g} m'):>7}; "
                          f"exc at 4 m {next(r for wd, r in zip(widths, res) if wd == 4)['max_exc']*1e3:5.1f}, at 8 m "
                          f"{next(r for wd, r in zip(widths, res) if wd == 8)['max_exc']*1e3:5.1f} mm")

    # ------------------------------------------------------------------ D
    if "D" in parts:
        hr("PART D: STABILITY OF THE ATTRACTIVE PULSE IN THE GAP")
        print(f"  Force on an iron plate across a gap g goes as 1/g^2 (B = mu0 NI / 2g), so the stiffness is -2F/g: the hold-down pull")
        print(f"  weakens as the casing rises away from the cable's top face, the hold-up pull strengthens as the casing rises toward")
        print(f"  the underside.  Growth time = sqrt(g mu / 2F).  (With F going as 1/g the stiffness would be F/g and the time 1.41 times")
        print(f"  longer; that is the 14 ms of the hand values.  The magnetic circuit here gives 1/g^2 and 10 ms at 6,000 kN/m.)")
        print(f"  {'u km/s':>7} " + " ".join(f"{'F kN/m L='+str(int(L)):>12} {'growth ms':>9} {'pulse ms':>8} {'ratio':>6}" for L in L_scs))
        for u in (1000, 2000, 3000, 5000, 8000, 10000, 15000, 20000, 25000, 30000):
            line = f"  {u/1e3:7.0f} "
            for L in L_scs:
                F = abs(float(q_zone(u, L)))
                tg = math.sqrt(GAP0 * MU / (2 * F)) if F > 0 else float('inf')
                tp = L / u
                line += f"{F/1e3:12.0f} {tg*1e3:9.1f} {tp*1e3:8.1f} {tp/tg:6.2f} "
            print(line)
        print(f"  Above the neutral speed the ratio pulse/growth is nearly constant (force and pulse length both scale with speed);")
        print(f"  below it, where the hold-up track carries a near-constant weight for longer, the ratio grows as the speed falls.")
        print(f"\n  Open-loop run with a 5 mm initial gap error (point model, ring target, 1 Hz steady loop, tau {tau0*1e3:.0f} ms), L_sc = {L_main:g} m:")
        print(f"  {'u km/s':>7} {'w0 mm':>6} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'v_end':>7} {'F_lim':>7} {'band':>6}   {'same, pulse regulated':>24}")
        for u in (1000.0, 2000.0, 3000.0, 5000.0, 10000.0, 20000.0, 30000.0):
            for w0 in (0.005, -0.005):
                lanes = [lane(L_sc=L_main, tau=tau0, w0=w0, B_sat=b_sat, steady_loop_hz=1.0), lane(L_sc=L_main, tau=tau0, w0=w0, B_sat=b_sat, steady_loop_hz=1.0, gap_dep=False)]
                res, _ = point_pulsed(lanes, u, tabs)
                r, r2 = res
                print(f"  {u/1e3:7.0f} {w0*1e3:+6.0f} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['v_rel_end']:+7.3f} "
                      f"{r['pk_lim']/1e3:7.2f} {('ok' if r['ok_band'] else 'FAIL'):>6}   exc {r2['max_exc']*1e3:6.1f} mm, v_end {r2['v_rel_end']:+.3f}")

    # ------------------------------------------------------------------ E
    if "E" in parts:
        hr("PART E: POWER AND ENERGY (hold-down track sized for case a at 1.2 T, fringe factor %.1f, core energy neglected)" % FRINGE_FACTOR)
        print(f"  Gap field energy per meter of track = P_sat x gap x width x fringe factor (the width counts both pole faces, so one")
        print(f"  gap length per unit of width covers the two gaps of a pole pair).  Energized length = zone + rise and fall lengths u tau.")
        print(f"  Power into new sections = energy per meter x u.  Net electrical energy per pass = stored x (1 - recovery {RECOVERY_EFF:.0%}) + "
              f"coil loss {COIL_LOSS_FRACTION:.0%} of stored.")
        NI = 1.2 * 2 * GAP0 / MU0
        print(f"  Ampere-turns per pole pair for 1.2 T across two {GAP0*1e3:.0f} mm gaps: {NI:,.0f} A-turns.")
        for L in L_scs:
            print(f"\n  L_sc = {L:g} m, tau {tau0*1e3:.0f} ms, sections {50:.0f} m.")
            print(f"  {'u km/s':>7} {'F_hd kN/m':>9} {'width m':>8} {'E J/m':>10} {'L_on m':>8} {'E_on MJ':>8} {'P GW':>7} {'P/P_lsm':>8} {'P/712GW':>8} {'E_sec MJ':>8} {'E_net TJ/lap':>12}")
            for u in (8000, 10000, 15000, 20000, 25000, 30000):
                F = abs(float(q_zone(u, L)))
                wd = F / P_SAT[1.2]
                E_m = P_SAT[1.2] * GAP0 * wd * FRINGE_FACTOR
                L_on = L + 2 * u * tau0 + 50.0
                E_on = E_m * L_on
                P = E_m * u
                E_sec = E_m * 50.0
                n_sections_per_pass = ring.L_RING / 50.0
                E_net = E_m * ring.L_RING * ((1 - RECOVERY_EFF) + COIL_LOSS_FRACTION)   # per lap, every meter energized once
                print(f"  {u/1e3:7.0f} {F/1e3:9.0f} {wd:8.2f} {E_m:10,.0f} {L_on:8.0f} {E_on/1e6:8.2f} {P/1e9:7.2f} {P/416.47e9:8.3f} {P/lsm.P_HVDC_MAX:8.3f} "
                      f"{E_sec/1e6:8.2f} {E_net/1e12:12.2f}")
            print(f"  E_net is per lap (every meter of ring energized once), not per section.  P_lsm = 416.47 GW peak from lsm_simulation.txt;")
            print(f"  712 GW = lsm_config.P_HVDC_MAX.  The bare sled's cancellation adds {abs(float(q_bare(30000.0)))*GAP0*1.0:,.0f} J/m over 10 km at 30 km/s "
                  f"for a 1 m wide track at its lower field ({math.sqrt(2*MU0*abs(float(q_bare(30000.0)))/1.0):.2f} T).")

    # ------------------------------------------------------------------ F
    if "F" in parts:
        hr("PART F: HEAT IN THE CABLE'S PLATES PER PASS [UNVERIFIED plate model]")
        mu_r = lc.FERRO_MU_R; rho_p = lc.FERRO_RESISTIVITY; d_p = lc.FERRO_PLATE_THICKNESS
        wh = lc.FERRO_HYSTERESIS_COEFF
        print(f"  Plate: sintered iron in resin, mu_r {mu_r}, resistivity {rho_p:.0e} ohm m, {d_p*1e3:.0f} mm thick, hysteresis "
              f"{wh} J/m^3 per cycle at saturation (levitation_config.py), density {lc.FERRO_DENSITY} kg/m^3, specific heat 450 J/kg K assumed.")
        print(f"  Each piece of plate sees one field cycle per pass: the energized region sweeps over it at u + {V_REL:,.0f} m/s.  The field")
        print(f"  edge is set by the pole fringing, about one gap (0.1 m), so the rise takes 0.1 / (u + V_rel), a few microseconds.")
        for L in L_scs:
            print(f"\n  L_sc = {L:g} m (hold-down track at 1.2 T, width from case a)")
            print(f"  {'u km/s':>7} {'width m':>7} {'edge us':>8} {'skin mm':>8} {'hyst J/m':>9} {'eddy J/m':>9} {'dT hyst uK':>10} {'dT eddy mK':>10}")
            for u in (10000, 20000, 30000):
                wd = abs(float(q_zone(u, L))) / P_SAT[1.2]
                v_sweep = u + V_REL
                t_edge = 0.1 / v_sweep
                om = math.pi / (2 * t_edge)                 # quarter-period rise
                skin = math.sqrt(2 * rho_p / (om * MU0 * mu_r))
                hyst = wh * d_p * wd                        # J per meter of plate, one full cycle at saturation
                B = 1.2
                eddy = 2 * (B**2 / (2 * MU0 * mu_r)) * min(skin, d_p) * wd   # field energy in the skin layer dissipated at each edge  [UNVERIFIED]
                cp = 450.0
                mass_plate = lc.FERRO_DENSITY * d_p * wd
                mass_skin = lc.FERRO_DENSITY * min(skin, d_p) * wd
                print(f"  {u/1e3:7.0f} {wd:7.2f} {t_edge*1e6:8.2f} {skin*1e3:8.1f} {hyst:9.1f} {eddy:9.1f} {hyst/(mass_plate*cp)*1e6:10.2f} "
                      f"{eddy/(mass_skin*cp)*1e3:10.2f}")
        print(f"  Missing: the plate's real high-frequency behavior (the 1e-4 ohm m effective resistivity and mu_r 200 are design guesses),")
        print(f"  the field edge shape, and the limiter-zone fields.  The hysteresis number is a per-cycle bulk figure, fine for a millisecond")
        print(f"  pulse if the loop closes; the eddy estimate takes the field energy in one skin depth as dissipated at each edge.")

    # ------------------------------------------------------------------ G
    ringres = {}
    if "G" in parts:
        hr("PART G: REACTION ON THE CABLE AND THE RING")
        for L in L_scs:
            for u in (10000, 20000, 30000):
                Jz = abs(float(q_zone(u, L))) * L / u
                Jb = abs(float(q_bare(u))) * (L_SLED - L) / u
                print(f"  L_sc {L:4g} m, u {u/1e3:3.0f} km/s: impulse on the cable per pass {Jz/1e3:6.2f} (zone) + {Jb/1e3:5.2f} (bare) kN s/m "
                      f"= {(Jz+Jb)/1e3:6.2f}; cable gains {(Jz+Jb)/MU_C:.3f} m/s; casing and cable together {(Jz+Jb)/(MU+MU_C):.3f} m/s; "
                      f"per lap the same, delivered once to every meter.")
        t_p, v_p, P_p = lsm_profile()
        print(f"\n  Launch profile from lsm_simulation (quick mode): {t_p[-1]/60:.1f} min to {v_p[-1]/1e3:.1f} km/s, peak LSM power {P_p.max()/1e9:.1f} GW.")
        for wa in (False, True):
            rr = ring_radial(t_p, v_p, wa)
            ringres[wa] = rr
            print(f"\n  {'With' if wa else 'Without'} anchor lines: k_hoop {rr['k_hoop']:.3e}, k_grav {rr['k_grav']:.3e}, k_anchor {rr['k_anchor']:.3e} N/m per m "
                  f"(sum {rr['k_eff']:.3e}); breathing period {rr['period']/60:.1f} min.")
            print(f"    Ring-average outward force at the end of the launch {rr['f'][-1]:.1f} N/m (sled total {rr['F_tot'][-1]/1e6:,.0f} MN); "
                  f"static radial expansion {rr['dr_static_end']:.1f} m; peak excursion in the run {rr['dr'].max():.1f} m, lowest {rr['dr'].min():.1f} m, "
                  f"at the end {rr['dr'][-1]:.1f} m.")
            if wa:
                print(f"    Anchor estimate: {rr['n_lines']:.0f} lines at 50 km spacing, uniform area {rr['A_line']:.3f} m^2 from the casing weight share at "
                      f"{ring.SIGMA_CABLE/1e9:.1f} GPa, {ring.ALTITUDE/1e3:.0f} km long, E 500 GPa: {rr['k_line']/1e3:.0f} kN/m per line.  [UNVERIFIED] The ring")
                print(f"    between anchors is treated as moving with them; the real motion bulges between anchor points.")
        print(f"  The 0.32 m/s per lap and 141 m excursion of the first report were the hoop term alone; the gravity gradient and the")
        print(f"  cable's centrifugal term stiffen the mode, and the drive is continuous around the ring, not a kick per lap.")

    # ------------------------------------------------------------------ J
    if "J" in parts:
        hr("PART J: THE PUSH-DOWN OPTION (driven traveling-wave coils under the cable, repelling a conducting plate)")
        print(f"  Thick-plate model: skin depth delta = sqrt(2 rho / (mu0 k v_slip)), drag/lift = k delta, heat per newton = v_slip k delta,")
        print(f"  push P = B_plate^2/(2 mu0) / (1 + (k delta)^2) [UNVERIFIED interpolation between the image limit and the weak regime],")
        print(f"  B_plate = B_surface exp(-k gap), drive frequency (V_rel - v_slip)/lambda, Halbach sheet current K = 1.414 B_surface/mu0,")
        print(f"  field energy above the coils B_surface^2/(2 mu0 k) per m^2, reactive power = 2 pi f x that.  Molybdenum plate, 20 mm, gap 100 mm.")
        print(f"  {'lambda m':>8} {'slip m/s':>8} {'delta mm':>8} {'k delta':>7} {'thick?':>6} {'heat W/N':>8} {'f kHz':>8} {'B_p/B_s':>7} "
              f"{'P kPa@1T':>8} {'K MA/m':>7} {'NI/pole':>8} {'Q MVA/m2':>9}")
        for wl in (0.2, 0.5, 1.0, 2.0):
            for vs in (10, 100, 1000, V_REL):
                pd = push_down(vs, wl, 1.0, GAP0)
                print(f"  {wl:8.1f} {vs:8.0f} {pd['delta']*1e3:8.2f} {pd['kd']:7.3f} {('yes' if pd['thick'] else 'no'):>6} {pd['heat_per_N']:8.1f} "
                      f"{pd['f_drive']/1e3:8.2f} {pd['B_p']:7.3f} {pd['P']/1e3:8.1f} {pd['K']/1e6:7.2f} {pd['NI_pole']:8,.0f} {pd['Q_reactive']/1e6:9.1f}")
        print(f"\n  Stability split.  Hold-down stiffness -2 F_hd / g; push-down stiffness +2 k F_pd (force grows as exp(-2 k gap)).")
        print(f"  Neutral when F_pd / F_hd = 1 / (k g) = lambda / (2 pi g).  For the case (a) force at 30 km/s, L_sc = {L_main:g} m:")
        F_tot = abs(float(q_zone(30000.0, L_main)))
        print(f"  {'lambda m':>8} {'F_pd/F_hd':>9} {'F_hd kN/m':>9} {'F_pd kN/m':>9} {'w_hd m':>7} {'B_s T':>6} {'w_pd m':>7} {'sum m':>6} {'heat kJ/m/pass':>15} {'dT K':>6} {'Q GVA':>7}")
        for wl in (0.2, 0.5, 1.0, 2.0):
            ratio = wl / (2 * math.pi * GAP0)
            F_hd = F_tot / (1 + ratio)
            F_pd = F_tot - F_hd
            w_hd = F_hd / P_SAT[1.2]
            for B_s in (1.0, 2.0):
                pd = push_down(100.0, wl, B_s, GAP0)
                w_pd = F_pd / pd["P"]
                mat = lp.get_material("molybdenum")
                heat_ring = F_pd * pd["heat_per_N"] * L_main / 30000.0           # J per meter of ring per pass
                heat_plate = F_pd * pd["heat_per_N"] / (30000.0 + V_REL)       # J per meter of plate per pass
                m_skin = mat["density"] * pd["delta"] * w_pd
                dT = heat_plate / (m_skin * mat["specific_heat"]) if m_skin > 0 else float('inf')
                Q = pd["Q_reactive"] * w_pd * (L_main + 2 * 30000.0 * tau0)
                print(f"  {wl:8.1f} {ratio:9.2f} {F_hd/1e3:9.0f} {F_pd/1e3:9.0f} {w_hd:7.2f} {B_s:6.1f} {w_pd:7.1f} {w_hd+w_pd:6.1f} {heat_plate/1e3:15.2f} {dT:6.1f} {Q/1e9:7.0f}")
        print(f"  Slip 100 m/s in the split table.  At zero slip the loss vanishes but so does the push; at full slip (static coils, the")
        print(f"  limiter case) the heat per newton is the {push_down(V_REL, 1.0, 1.0, GAP0)['heat_per_N']:.0f} W/N of the first row block.")
        used = 2 * lc.EML_PLATE_WIDTH + 2 * lc.SAFETY_PLATE_WIDTH
        print(f"  Underside width already used: two steady EML plates {2*lc.EML_PLATE_WIDTH:.1f} m, two limiter plates {2*lc.SAFETY_PLATE_WIDTH:.1f} m, "
              f"= {used:.1f} m of {CABLE_FACE:.2f} m; the hold-up plate and a push-down plate must fit in the remaining {CABLE_FACE-used:.2f} m.")
        print(f"  The limiter plates are molybdenum, 20 mm, 1 m wide each; a push-down coil row could face them, but their width is 2 m in total")
        print(f"  and their gap is the limiter gap, so a push-down track of more than 2 m needs its own plate.")

    # ------------------------------------------------------------------ K
    rows_K = []
    if "K" in parts:
        hr("PART K: ACTING AHEAD OF THE SLED")
        print(f"  Lead for the rise time: each section starts {LEAD_TAUS:.0f} tau early.  Lead distance u tau:")
        for tau in TAU_CASES:
            print(f"    tau {tau*1e3:3.0f} ms: " + ", ".join(f"{u/1e3:.0f} km/s -> {u*tau:,.0f} m" for u in (5000, 10000, 20000, 30000)))
        print(f"\n  Spread pull: the hold-down impulse J spread evenly over T centered on the pass.  Closed form: excursion J T / (8 mu), force J / T.")
        print(f"  Point model with the gap dependence (instability) on and off, ring target, tau {tau0*1e3:.0f} ms, 5,000 t:")
        Ts = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
        for L in L_scs:
            for u in (10000.0, 20000.0, 30000.0):
                Jz = abs(float(q_zone(u, L))) * L / u * MATCH_RING
                print(f"\n  L_sc = {L:g} m, u = {u/1e3:.0f} km/s: J (ring target) = {Jz/1e3:.2f} kN s/m, free casing speed {Jz/MU:.3f} m/s.")
                print(f"  {'T ms':>6} {'J T/8mu mm':>10} {'F kN/m':>8} {'width m':>8} {'L_on m':>8} {'P GW':>7} {'exc mm':>7} {'band':>6} {'v_end':>7} "
                      f"{'exc no-inst':>11} {'growth ms':>9} {'T/growth':>8}")
                lanes = []
                for T in Ts:
                    lanes.append(lane(L_sc=L, tau=tau0, T_spread=(T or None), B_sat=b_sat, steady_loop_hz=1.0, label=f"T{T}"))
                    lanes.append(lane(L_sc=L, tau=tau0, T_spread=(T or None), B_sat=b_sat, steady_loop_hz=1.0, gap_dep=False, label=f"T{T}n"))
                res, _ = point_pulsed(lanes, u, tabs)
                for T in Ts:
                    r = next(r for r in res if r["label"] == f"T{T}")
                    r2 = next(r for r in res if r["label"] == f"T{T}n")
                    rows_K.append(r)
                    Teff = T if T else L / u
                    F = Jz / Teff
                    wd = F / P_SAT[1.2]
                    L_on = u * Teff + L + 2 * u * tau0
                    P = P_SAT[1.2] * GAP0 * wd * u
                    tg = math.sqrt(GAP0 * MU / (2 * F))
                    print(f"  {Teff*1e3:6.1f} {Jz*Teff/(8*MU)*1e3:10.1f} {F/1e3:8.0f} {wd:8.2f} {L_on:8.0f} {P/1e9:7.2f} {r['max_exc']*1e3:7.1f} "
                          f"{('ok' if r['ok_band'] else 'FAIL'):>6} {r['v_rel_end']:+7.3f} {r2['max_exc']*1e3:11.1f} {tg*1e3:9.1f} {Teff/tg:8.2f}")
        print(f"  The longest safe T is the last row with 'ok' and a small v_end; past it the attractive instability runs away during the pull.")
        print(f"  With the push-down track taking its neutral share (Part J), the instability term is removed; the no-inst column is that case.")

    # ------------------------------------------------------------------ H
    rows_H = []
    if "H" in parts:
        hr("PART H: FAILURE CASES (beam model: 10 m truss and zero stiffness, cable as a string, 50 m sections, tau %.0f ms)" % (tau0 * 1e3))
        trusses = [sp.truss(10.0, 0.1), sp.ZERO_STIFFNESS]
        for u in (10000.0, 30000.0):
            cases = [dict(label="normal, regulated"),
                     dict(label="normal, open loop", regulated=False),
                     dict(label="one 50 m section out", fail_x=(3000.0, 3050.0)),
                     dict(label="200 m out", fail_x=(2900.0, 3100.0)),
                     dict(label="1 km supply out", fail_x=(2500.0, 3500.0)),
                     dict(label="late 10 ms", late=10e-3),
                     dict(label="late 2 ms", late=2e-3),
                     dict(label="spurious, no sled", spurious=True),
                     dict(label="no pulse", no_pulse=True)]
            t0 = time.time()
            res, _ = beam_pulsed(u, L_main, trusses, cases, tabs, L_section=50.0, tau=tau0)
            rows_H.extend([dict(r, u=u) for r in res])
            print(f"\n  u = {u/1e3:.0f} km/s, L_sc = {L_main:g} m ({time.time()-t0:.0f} s; dt {res[0]['dt']*1e6:.0f} us, {res[0]['nt']:,} steps):")
            print(f"  {'case':>22} {'truss':>22} {'exc mm':>7} {'gap_b':>6} {'gap_t':>6} {'v_pk':>6} {'F_lim kN/m':>10} {'contact':>7}")
            for r in res:
                print(f"  {r['case']:>22} {r['truss']:>22} {r['max_exc']*1e3:7.1f} {r['min_gap_b']*1e3:6.1f} {r['min_gap_t']*1e3:6.1f} {r['v_pk']:6.3f} "
                      f"{r['pk_lim']/1e3:10.1f} {('yes' if r['contact'] else 'no'):>7}")
        print(f"\n  A spurious full pulse with no sled (point model): the hold-down delivers the whole impulse downward.")
        for u in (10000.0, 30000.0):
            lanes = [lane(L_sc=L_main, tau=tau0, B_sat=b_sat, M_sc=M_SC0, label="spurious")]
            # emulate: no sled load but the pulse fires: use mass_err to keep the command, and remove the sled by M_sc=0 in q only
            res, hist = point_pulsed(lanes, u, tabs)
            # point_pulsed cannot zero q alone; compute the spurious case analytically from J_hd
            J = res[0]["J_hd"]
            v_k = J / MU
            U_top = tabs["F0"] / (2 * tabs["k"])
            print(f"    u {u/1e3:3.0f} km/s: pulse impulse {J/1e3:.2f} kN s/m -> casing {v_k:.2f} m/s toward the cable's top face; the top limiter stores "
                  f"{U_top/1e3:.2f} kJ/m, enough for {math.sqrt(2*U_top/MU):.2f} m/s, so it {'holds' if 0.5*MU*v_k**2 < U_top else 'does not hold'}; "
                  f"within 20 mm it stores {sp.limiter_U(GAP_FORBIDDEN, tabs)/1e3:.2f} kJ/m, {math.sqrt(2*sp.limiter_U(GAP_FORBIDDEN, tabs)/MU):.2f} m/s.")
        print(f"  Largest fault the limiters catch: a kick of {math.sqrt(2*tabs['F0']/(2*tabs['k'])/MU):.2f} m/s either way before contact "
              f"({math.sqrt(2*tabs['F0']/(2*tabs['k'])/MU)*MU/1e3:.1f} kN s/m, {math.sqrt(2*tabs['F0']/(2*tabs['k'])/MU)*MU/abs(float(q_zone(30000.0, L_main)))/L_main*30000.0*100:.0f} % "
              f"of a 30 km/s pass on {L_main:g} m), and {math.sqrt(2*sp.limiter_U(GAP_FORBIDDEN, tabs)/MU):.2f} m/s inside the 20 mm zone.")

    # ------------------------------------------------------------------ I
    rows_I = []
    if "I" in parts:
        hr("PART I: THE PARKED SLED (beam model, static zone load, ideal steady controller with a force ceiling)")
        trusses = [sp.truss(10.0, 0.1), sp.truss(20.0, 0.1), sp.truss(30.0, 0.1), sp.ZERO_STIFFNESS]
        for L in L_scs:
            F_need = W + float(q_zone(0.0, L))
            for cap in (F_STEADY_CAP, F_need):
                cases = [dict(label=f"cap {cap/1e3:.0f} kN/m")]
                res, prof = beam_pulsed(0.0, L, trusses, cases, tabs, parked=True, F_steady_cap=cap, L_win=4000.0)
                for il, r in enumerate(res):
                    S = prof["S"][il]; rel = prof["rel"][il]; x = prof["x"]
                    spread = x[S > W + 1e3]
                    print(f"  L_sc {L:4g} m, needs {F_need/1e3:5.0f} kN/m, ceiling {cap/1e3:5.0f} kN/m, {r['truss']:>22}: top gap {r['min_gap_t']*1e3:5.1f} mm, "
                          f"load spread over {(spread.max()-spread.min()) if spread.size else 0:6.0f} m, peak steady force {r['pk_steady']/1e3:5.0f} kN/m, "
                          f"{'contact' if r['contact'] else 'no contact'}")
                rows_I.extend([dict(r, L_sc=L, cap=cap) for r in res])
        print(f"  Steady track width needed at a loading station (400 A, {sp.F_EML_OP/2e3:.0f} kN/m per meter of track): " +
              ", ".join(f"{L:g} m: {(W + float(q_zone(0.0, L)))/(sp.F_EML_OP/2.0):.1f} m" for L in L_scs))

    # ------------------------------------------------------------------ hand values
    hr("HAND VALUES OF THE INSTRUCTIONS (section 4)")
    checks = []
    checks.append(("Pull per square meter at 1.2 T (kPa)", 573, P_SAT[1.2] / 1e3))
    for L, hv_F, hv_w in ((200.0, 3060, 5.3), (100.0, 6150, 10.7), (50.0, 12330, 21.5)):
        F = -float(q_zone(30000.0, L)) - W
        checks.append((f"Hold-down force, {L:g} m zone, steady at zero (kN/m)", hv_F, F / 1e3))
        checks.append((f"  its width at 1.2 T (m)", hv_w, F / P_SAT[1.2]))
    F100 = -float(q_zone(30000.0, 100.0)) - W
    checks.append(("Field across the 8.24 m face for the 100 m zone (T)", 1.4, math.sqrt(2 * MU0 * F100 / CABLE_FACE)))
    checks.append(("Instability growth time at 6,000 kN/m (ms)", 14, math.sqrt(GAP0 * MU / (2 * 6e6)) * 1e3))
    checks.append(("Pulse length at 30 km/s, 100 m (ms)", 3.3, 100 / 30000 * 1e3))
    checks.append(("Pulse length at 30 km/s, 200 m (ms)", 6.7, 200 / 30000 * 1e3))
    checks.append(("Gap field energy, 5.3 m track at 1.2 T (MJ/m)", 0.30, P_SAT[1.2] * GAP0 * 5.3 / 1e6))
    checks.append(("Power into new sections at 30 km/s (GW)", 9, P_SAT[1.2] * GAP0 * 5.3 * 30000 / 1e9))
    checks.append(("Ampere-turns per pole pair, 1.2 T, two 100 mm gaps", 191000, 1.2 * 2 * GAP0 / MU0))
    checks.append(("Timing error 1 ms at 30 km/s (m)", 30, 30000 * 1e-3))
    J30 = abs(float(q_zone(30000.0, 100.0))) * 100 / 30000
    checks.append(("Impulse on the cable per pass (kN s/m)", 20, J30 / 1e3))
    checks.append(("Cable speed gained (m/s)", 0.17, J30 / MU_C))
    checks.append(("Push-down heat per N, Mo, 1 m, slip 8,825 (W/N)", 69, push_down(V_REL, 1.0, 1.0, GAP0)["heat_per_N"]))
    checks.append(("Push-down heat per N, slip 100 (W/N)", 7.3, push_down(100.0, 1.0, 1.0, GAP0)["heat_per_N"]))
    checks.append(("Skin depth, slip 100 (mm)", 11.6, push_down(100.0, 1.0, 1.0, GAP0)["delta"] * 1e3))
    checks.append(("Drive frequency, 1 m, slip 100 (kHz)", 8.7, push_down(100.0, 1.0, 1.0, GAP0)["f_drive"] / 1e3))
    print(f"  {'Quantity':58} {'hand':>9} {'script':>10} {'check':>6}")
    for name, hv, sv in checks:
        rel = abs(sv - hv) / abs(hv)
        print(f"  {name:58} {hv:9.4g} {sv:10.4g} {('pass' if rel <= 0.05 else 'FAIL'):>6}")
    print(f"  The growth time hand value of 14 ms assumes a stiffness of F/g; the 1/g^2 force law of the magnetic circuit gives 2F/g and 10 ms.")

    # ------------------------------------------------------------------ graphs
    if graphs:
        hr("GRAPHS")
        u_f = np.linspace(0, 30000, 601)
        # 11: force required from each pulsed track, with capacities
        fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
        for L, col in zip(L_scs, PALETTE):
            qz = q_zone(u_f, L)
            ax.plot(u_f / 1e3, np.maximum(-qz, 0) / 1e3, color=col, lw=1.5, label=f"hold-down, L_sc = {L:g} m")
            ax.plot(u_f / 1e3, np.maximum(qz, 0) / 1e3, color=col, lw=1.5, ls="--", label=f"hold-up, L_sc = {L:g} m")
        for wd, ls in ((2, ":"), (4, "-."), (8.24, (0, (5, 2)))):
            ax.axhline(wd * P_SAT[1.2] / 1e3, color="gray", ls=ls, lw=1, label=f"{wd:g} m of track at 1.2 T")
        ax.set_yscale("log"); ax.set_ylim(10, 2e4)
        ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel("Pulsed track force needed, case a (kN/m)")
        ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7, ncol=2)
        _save(fig, "fig_sled_11_pulsed_force.png", out_dir)
        # 12: width needed
        fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
        for L, col in zip(L_scs, PALETTE):
            qz = np.abs(q_zone(u_f, L))
            ax.plot(u_f / 1e3, qz / P_SAT[1.2], color=col, lw=1.5, label=f"L_sc = {L:g} m, 1.2 T")
            ax.plot(u_f / 1e3, qz / P_SAT[1.6], color=col, lw=1.2, ls="--", label=f"L_sc = {L:g} m, 1.6 T")
        ax.axhline(CABLE_FACE, color="black", lw=1, label=f"cable face {CABLE_FACE:.2f} m")
        ax.set_ylim(0, 25); ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel("Pulsed track width needed, case a (m)")
        ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2)
        _save(fig, "fig_sled_12_pulsed_width.png", out_dir)
        # 13: smallest gap with the pulsed tracks (baseline)
        if sweep_rows:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            for L, col in zip(L_scs, PALETTE):
                for lab, ls in (("reg+loop", "-"), ("regulated", "--"), ("open+loop", "-."), ("no pulse+loop", ":")):
                    sel = sorted([r for r in sweep_rows if r["L_sc"] == L and r["label"] == lab], key=lambda r: r["u"])
                    ax.plot([r["u"] / 1e3 for r in sel], [min(r["min_gap_b"], r["min_gap_t"]) * 1e3 for r in sel], color=col, ls=ls, lw=1.4,
                            label=f"L_sc = {L:g} m, {lab}")
            ax.axhspan(80, 120, color="gray", alpha=0.12); ax.axhline(GAP_FORBIDDEN * 1e3, color="gray", ls="--", lw=1)
            ax.set_ylim(-3, 105); ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel("Smallest gap, either side (mm)")
            ax.grid(alpha=0.3); ax.legend(fontsize=6.5, ncol=2)
            _save(fig, "fig_sled_13_pulsed_min_gap.png", out_dir)
        # 14: excursion at 30 km/s against timing error per rise time
        if rows_B:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            for tau, col in zip(TAU_CASES, PALETTE):
                des = sorted(DT_ERR_CASES)
                ys = [next(r for r in rows_B if r["L_sc"] == L_main and r.get("reg") == "regulated" and r["label"] == f"tau{tau*1e3:g} dt{de*1e3:+g}")["max_exc"] * 1e3 for de in des]
                ax.plot([d * 1e3 for d in des], ys, color=col, marker="o", ms=3, lw=1.4, label=f"rise time {tau*1e3:g} ms")
            ax.axhline(BAND * 1e3, color="gray", ls="--", lw=1); ax.axhline(BAND_STRICT * 1e3, color="gray", ls=":", lw=1)
            ax.set_xlabel("Timing error of the pulse (ms, positive = late)"); ax.set_ylabel(f"Largest gap excursion at 30 km/s, L_sc = {L_main:g} m (mm)")
            ax.grid(alpha=0.3); ax.legend(fontsize=8)
            _save(fig, "fig_sled_14_timing_error.png", out_dir)
        # 15: displacement against time, with and without the pulse
        fig, axes = plt.subplots(2, 1, figsize=(sp.GRAPH_W, sp.GRAPH_H * 1.4))
        for u, col in zip((5000.0, 10000.0, 20000.0, 30000.0), PALETTE):
            lanes = [lane(L_sc=L_main, tau=tau0, B_sat=b_sat, steady_loop_hz=1.0, gap_dep=False), lane(L_sc=L_main, tau=tau0, B_sat=b_sat, steady_loop_hz=1.0, cancel=False)]
            res, h = point_pulsed(lanes, u, tabs)
            tm = h["t"] * 1e3
            axes[0].plot(tm, (h["w"][0] - h["wc"][0]) * 1e3, color=col, lw=1.3, label=f"{u/1e3:g} km/s, pulsed tracks")
            axes[0].plot(tm, (h["w"][1] - h["wc"][1]) * 1e3, color=col, lw=1.1, ls=":", label=f"{u/1e3:g} km/s, no pulse")
            axes[1].plot(tm, h["F"][0] / 1e3, color=col, lw=1.3)
        axes[0].axhspan(-20, 20, color="gray", alpha=0.12); axes[0].set_ylim(-105, 105); axes[0].set_xlim(-20, 400)
        axes[0].set_ylabel("Casing offset from the cable (mm, up positive)"); axes[0].legend(fontsize=7, ncol=2)
        axes[1].set_xlim(-20, 400); axes[1].set_ylabel("Pulsed track force on the casing (kN/m, up positive)")
        axes[1].set_xlabel(f"Time after the sled front passes the point (ms), L_sc = {L_main:g} m")
        for a_ in axes: a_.grid(alpha=0.3)
        _save(fig, "fig_sled_15_pulsed_time_history.png", out_dir)
        # 16: power and energized length
        fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H)); ax2 = ax.twinx()
        for L, col in zip(L_scs, PALETTE):
            F = np.abs(q_zone(u_f, L)); wd = F / P_SAT[1.2]
            ax.plot(u_f / 1e3, P_SAT[1.2] * GAP0 * wd * u_f / 1e9, color=col, lw=1.5, label=f"power, L_sc = {L:g} m")
            ax2.plot(u_f / 1e3, L + 2 * u_f * tau0 + 50, color=col, lw=1.2, ls="--", label=f"energized length, L_sc = {L:g} m")
        ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel("Power into newly energized sections (GW)")
        ax2.set_ylabel(f"Energized length, zone plus rise and fall at tau = {tau0*1e3:g} ms plus one section (m)")
        ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left"); ax2.legend(fontsize=7, loc="lower right")
        _save(fig, "fig_sled_16_pulsed_power.png", out_dir)
        # 17: ring radial motion
        if ringres:
            fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
            for wa, col in ((False, PALETTE[0]), (True, PALETTE[1])):
                rr = ringres[wa]
                ax.plot(rr["t"] / 60, rr["dr"], color=col, lw=1.4, label=f"{'with' if wa else 'without'} anchor lines, period {rr['period']/60:.0f} min")
            ax.set_xlabel("Time from the start of the launch (min)"); ax.set_ylabel("Radial expansion of the whole ring (m)")
            ax.grid(alpha=0.3); ax.legend(fontsize=8)
            _save(fig, "fig_sled_17_ring_radial.png", out_dir)
        # 18: mass against speed for fixed hold-down widths (full cancel, closed form)
        fig, ax = plt.subplots(figsize=(sp.GRAPH_W, sp.GRAPH_H))
        us = np.linspace(8000, 30000, 200)
        for wd, col in zip((2, 4, 6, 8), PALETTE):
            m = L_main * (wd * P_SAT[1.2] / (-a_down(us)) - M_SLED_M)
            ax.plot(us / 1e3, m / 1e3, color=col, lw=1.5, label=f"{wd:g} m hold-down at 1.2 T")
        m = L_main * (CABLE_FACE * P_SAT[1.2] / (-a_down(us)) - M_SLED_M)
        ax.plot(us / 1e3, m / 1e3, color="black", lw=1.2, ls="--", label=f"{CABLE_FACE:.2f} m, the whole cable face")
        ax.axhline(M_SC0 / 1e3, color="gray", ls=":", lw=1, label="5,000 t")
        ax.set_yscale("log"); ax.set_ylim(300, 50000)
        ax.set_xlabel("Sled speed relative to the casing (km/s)"); ax.set_ylabel(f"Largest spacecraft mass fully cancelled, L_sc = {L_main:g} m (t)")
        ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
        _save(fig, "fig_sled_18_pulsed_mass_envelope.png", out_dir)

    hr(f"DONE in {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
