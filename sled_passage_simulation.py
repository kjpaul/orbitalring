#!/usr/bin/env python3
"""
Sled passage load on the levitation bearing.

What the casing and its levitation bearing do while a mass driver sled passes,
for sled speeds from zero to 30 km/s relative to the casing, at the design
altitude set in ring_altitude.py (700 km, retrograde cable).

Three parts:

  A. Closed form. Load per meter on the casing against sled speed, under the
     bare sled and under the spacecraft zone, the speed at which each cancels
     the casing weight, passage times, impulse per meter and the speed that
     impulse gives a free meter of casing.

  B. The hypersonic excess velocity equations (Mach number, effective tension,
     characteristic length) over the whole run, for the cable and for the
     casing.  The casing is a truss beam whose depth is the design variable.
     Chord material, depth and chord mass fraction are printed as a table.

  C. Time domain.  The casing is integrated directly as a Timoshenko truss
     beam (and in the zero stiffness limit, each meter on its own) with its
     weight, the moving sled load, the EML pull under several control cases,
     and the passive Halbach limiter forces from levitation_physics.py.
     Nothing is assumed about steady state.  This is the fact check.

Usage:
    python sled_passage_simulation.py                   # full run at the design point, graphs on
    python sled_passage_simulation.py --no-graphs
    python sled_passage_simulation.py --l_sc=100        # one spacecraft load length (m) instead of 50, 100, 200
    python sled_passage_simulation.py --truss_depth=20  # one truss depth (m) instead of 10, 20, 30
    python sled_passage_simulation.py --chord_fraction=0.2
    python sled_passage_simulation.py --web_fraction=0.1
    python sled_passage_simulation.py --material=steel  # chord material for the time domain runs (CNT or steel)
    python sled_passage_simulation.py --workers=4       # parallel processes for the sweep (default: cores - 2, at most 6)
    python sled_passage_simulation.py --plate=molybdenum  # limiter reaction plate material (levitation_config default: copper_OFHC)
    python sled_passage_simulation.py --quick           # coarser speed steps, for a fast look
    python sled_passage_simulation.py --help

Frame convention (Paul, 2026-10-09): the sled speed u is relative to the
casing, because the LSM stator is on the casing.  The sled's inertial speed is
u + V_GROUND_SYNC.  Note that lsm_simulation.py starts the sled at v = 0 and
uses that v directly in v^2/R, so its radial load treats v as inertial; its
12.2 g at 30 km/s is 12.6 g in this convention.  Both are printed.

Mass driver rails are taken on top of the casing (Chapter 8 and 14 layout).
The spacecraft zone is at the front of the sled.  The vertical force balance
does not depend on either choice, but the order in which a point of casing
sees the zone and the bare sled does.

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

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

# =============================================================================
# RING, SLED, BEARING VALUES (nothing hard coded that ring_altitude.py gives)
# =============================================================================

R_ORBIT = ring.R_ORBIT
V_ORBIT = ring.V_ORBIT
V_GS = ring.V_GROUND_SYNC          # casing speed, eastward, inertial (m/s)
G_LOCAL = ring.G_LOCAL
G_NET = ring.G_NET
MU = ring.M_LOAD_M                 # casing mass per meter (kg/m)
V_REL = ring.V_REL                 # cable speed relative to the casing (m/s)
V_CABLE = ring.V_CABLE
T_CABLE = ring.CABLE_TENSION       # N
MU_CABLE = ring.M_CABLE_M          # cable with hardware (kg/m)
MU_CABLE_STRUCT = ring.M_CABLE_STRUCTURAL
CABLE_SIDE = ring.CABLE_SIDE
RHO_CNT = ring.RHO_CNT
E_CNT = 500e9                      # Pa, Chapter 3 (secant modulus at 5 % break strain)
W_CASING = MU * G_NET              # casing weight per meter (N/m)

L_SLED = float(lsm.L_SLED)         # m
M_SLED_M = float(lsm.M_SLED_PER_M) # kg/m sled hardware
M_SC = float(lsm.M_SPACECRAFT)     # kg spacecraft

# EML bearing (Chapter 8, Equations 8.1 and 8.2).  Two tracks under the cable,
# pull only.  The Halbach model in levitation_physics.py gives larger numbers
# (409 and 1,604 kN/m in the 2026-10-04 run); the chapter values are used here
# and the difference is reported.
F_EML_OP = 218e3                   # N/m at 400 A
F_EML_MAX = 871e3                  # N/m at 800 A (stated maximum)
GAP0 = lc.EML_NOMINAL_GAP          # 0.100 m nominal gap, contact at 0
GAP_FORBIDDEN = lc.FORBIDDEN_ZONE_INNER   # 0.020 m
U_MAX = 30_000.0                   # top sled speed relative to the casing (m/s)

G0_STD = 9.80665

# Time domain settings
PD_FREQ_HZ = 10.0                  # ideal controller return to the nominal gap: critically damped at this frequency
L_WIN_SUBSONIC = 20_000.0          # window length (m) for u <= U_WINDOW_SWITCH, where the free ends can matter
L_WIN_SUPERSONIC = 6_000.0         # window length (m) above it
U_WINDOW_SWITCH = 5_000.0          # m/s; the truss shear wave speeds are 1.9 to 3.3 km/s (CNT) and below 1 km/s (steel)
T_SETTLE = 1.5                     # s after the sled tail leaves the measurement region

# Color-blind safe palette (same family as anchor_line_analysis.py)
PALETTE = ['#0072B2', '#009E73', '#E69F00', '#CC0000', '#56B4E9', '#D55E00', '#CC79A7', '#000000']
GRAPH_DPI = 300
GRAPH_W, GRAPH_H = 10, 6

# =============================================================================
# PART A: CLOSED FORM LOADS
# =============================================================================

def a_down(u):
    """Net downward acceleration felt by a mass riding the casing at speed u (m/s^2).
    Positive down.  u is relative to the casing; inertial speed is u + V_GS."""
    v_i = np.asarray(u, dtype=float) + V_GS
    return G_LOCAL - v_i**2 / R_ORBIT


def a_down_inertial(v_i):
    """Same, with v_i already inertial (the lsm_simulation.py convention)."""
    return G_LOCAL - np.asarray(v_i, dtype=float)**2 / R_ORBIT


def mass_zone(L_sc):
    """Mass per meter on the casing in the spacecraft zone (kg/m)."""
    return M_SLED_M + M_SC / L_sc


def load_bare(u):
    """Load per meter under the bare sled, N/m, positive down."""
    return M_SLED_M * a_down(u)


def load_zone(u, L_sc):
    """Load per meter under the spacecraft zone, N/m, positive down."""
    return mass_zone(L_sc) * a_down(u)


def neutral_speed(m_per_m):
    """Inertial speed at which m_per_m kg/m pulls up with exactly the casing weight.
    Returns (v_inertial, u_relative)."""
    a_up_needed = W_CASING / m_per_m
    v_i = math.sqrt((G_LOCAL + a_up_needed) * R_ORBIT)
    return v_i, v_i - V_GS


def impulse_zone(u, L_sc):
    """Impulse per meter of ring from the spacecraft zone, N s/m, positive down.
    Each meter sees the zone load for L_sc/u seconds."""
    u = np.asarray(u, dtype=float)
    with np.errstate(divide='ignore'):
        return np.where(u > 0, load_zone(u, L_sc) * L_sc / u, np.inf)


def impulse_bare(u):
    """Impulse per meter from the bare sled over the whole sled length, N s/m, positive down."""
    u = np.asarray(u, dtype=float)
    with np.errstate(divide='ignore'):
        return np.where(u > 0, load_bare(u) * L_SLED / u, np.inf)


# =============================================================================
# PART B: TRUSS CASING AND THE HYPERSONIC EQUATIONS
# =============================================================================

MATERIALS = {
    "CNT":   dict(E=E_CNT, rho=RHO_CNT, label="CNT, E = 500 GPa, rho = 1,700 kg/m^3"),
    "steel": dict(E=200e9, rho=7850.0, label="structural steel, E = 200 GPa, rho = 7,850 kg/m^3"),
}


def truss(h, f, material="CNT", web_fraction=None, theta_deg=45.0):
    """Truss beam properties for a casing of depth h (m) whose two chords carry
    a fraction f of the 12,000 kg/m, with diagonals at theta to the chords
    carrying web_fraction of the mass (default half the chord fraction).

      EI  = E * A_chords * (h/2)^2              bending stiffness (N m^2)
      J   = f * mu * (h/2)^2                    rotational inertia per meter (kg m)
      GA  = E * A_web * sin^2(theta) cos^2(theta)   effective shear stiffness (N)
      c_s = sqrt(GA / mu)                       shear wave speed, the top speed of the bending branch
      c_b = sqrt(EI / J) = sqrt(E / rho)        rotational branch, the bar wave in the chords

    GA follows from a Warren truss: the shear V is carried by the diagonal
    force V / sin(theta), the diagonal strain is gamma sin(theta) cos(theta)
    for a shear angle gamma, and the diagonal mass per meter of truss is
    rho A_d / cos(theta).  Chord shear (Vierendeel action) is ignored.
    """
    m = MATERIALS[material]
    E, rho = m["E"], m["rho"]
    if web_fraction is None:
        web_fraction = 0.5 * f
    A_ch = f * MU / rho
    A_web = web_fraction * MU / rho
    EI = E * A_ch * (h / 2.0) ** 2
    J = f * MU * (h / 2.0) ** 2
    th = math.radians(theta_deg)
    GA = E * A_web * math.sin(th) ** 2 * math.cos(th) ** 2
    return dict(label=f"h={h:g} m, f={f:g}, {material}", h=h, f=f, material=material,
                web_fraction=web_fraction, EI=EI, GA=GA, J=J,
                c_s=math.sqrt(GA / MU), c_b=math.sqrt(EI / J))


ZERO_STIFFNESS = dict(label="zero stiffness", h=0.0, f=0.0, material="none", web_fraction=0.0,
                      EI=0.0, GA=0.0, J=1.0, c_s=0.0, c_b=0.0)


def cable_properties():
    """Cable bending stiffness and wave speed, as in Chapter 3 (square section)."""
    I = CABLE_SIDE ** 4 / 12.0
    EI = E_CNT * I
    v_wave = math.sqrt(T_CABLE / MU_CABLE)                # transverse wave speed, full mass
    v_wave_struct = math.sqrt(T_CABLE / MU_CABLE_STRUCT)  # structural mass only (anchor script form)
    return dict(I=I, EI=EI, v_wave=v_wave, v_wave_struct=v_wave_struct)


def hypersonic_cable(u, cab):
    """Mach number, excess tension and characteristic length for the cable, for
    a load moving at u + V_REL along it.

    In the frame moving with the load, a tensioned beam of mass per meter mu
    obeys  EI W'''' + (mu v^2 - T) W'' = q.  The coefficient mu v^2 - T =
    T (M^2 - 1) is the excess term.  It has the sign of a COMPRESSIVE axial
    force, so it does not stiffen the beam against the load.  Ahead of the
    load nothing moves (M > 1); behind it the term produces a trailing wake
    whose mean slope leaves each point moving at F/(mu v) * M^2/(M^2 - 1).
    L_char = sqrt(EI / T_eff) is the wake length scale (Chapter 3, Eq. 3.10 form).
    """
    v = np.asarray(u, dtype=float) + V_REL
    M = v / cab["v_wave"]
    T_eff = np.where(M > 1, T_CABLE * (M**2 - 1), T_CABLE)
    L_char = np.sqrt(cab["EI"] / T_eff)
    L_char_anchor_form = (cab["EI"] / T_eff) ** 0.25   # the form in anchor_line_analysis.py, dimensionally wrong
    return dict(v=v, M=M, T_eff=T_eff, L_char=L_char, L_char_anchor_form=L_char_anchor_form)


def hypersonic_casing(u, tr):
    """Same for the casing, which has no pretension.  In the moving frame
    EI W'''' + mu u^2 W'' = q, so mu u^2 plays the part of the excess tension,
    again with the compressive sign.  L_char = sqrt(EI / (mu u^2)) is 1/k of
    the trailing bending wake (the wave whose phase speed equals u).
    The Mach number uses the shear wave speed, the top speed of the bending
    branch of a Timoshenko beam; above it the casing ahead of the load gets
    no warning at all."""
    u = np.asarray(u, dtype=float)
    M = np.where(tr["c_s"] > 0, u / max(tr["c_s"], 1e-9), np.inf)
    with np.errstate(divide='ignore'):
        T_eff = MU * u**2
        L_char = np.where(u > 0, np.sqrt(tr["EI"] / np.maximum(T_eff, 1e-30)), np.inf)
    return dict(M=M, T_eff=T_eff, L_char=L_char)


def string_term_deflection(q, L_sc, T_eff):
    """Deflection at the rear edge of a load q (N/m) of length L_sc moving at
    supersonic speed, from the string term alone: q L_sc^2 / (2 T_eff).
    For the casing (T_eff = mu u^2) this is the same as 1/2 a t^2 with
    a = q/mu and t = L_sc/u, the hand estimate."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.abs(q) * L_sc**2 / (2.0 * T_eff)


# =============================================================================
# LIMITER LAW FROM levitation_physics.py
# =============================================================================

def limiter_tables(n=2501, g_max=0.25):
    """Tabulate the passive Halbach limiter force (N/m) and reaction plate
    eddy heating (W/m) against gap, using the functions in
    levitation_physics.py with the settings in levitation_config.py and the
    true cable-casing relative speed ring_altitude.V_REL.  At this speed the
    force is in the image-current limit, F = F0 exp(-2 k gap) exactly, and the
    thick-plate heating has the same gap dependence; the time domain solver
    uses that form with F0 and Q0 taken from the tabulated functions, and the
    largest deviation from the table is reported."""
    wavelength = lc.SAFETY_WAVELENGTH
    coils_per_m = 4.0 / wavelength
    B_surf = lp.halbach_surface_field(lc.SAFETY_TURNS_PER_COIL, lc.SAFETY_I_PERSISTENT,
                                      lc.SAFETY_COIL_WIDTH, wavelength, coils_per_meter=coils_per_m)
    mat = lp.get_material(lc.SAFETY_PLATE_MATERIAL)
    sigma = 1.0 / mat["resistivity"]
    gaps = np.linspace(0.0, g_max, n)
    F = np.array([lp.safety_repulsive_force(B_surf, g, wavelength, sigma, lc.SAFETY_PLATE_THICKNESS,
                                            lc.SAFETY_PLATE_WIDTH, V_REL, n_arrays=2) for g in gaps])
    Q = np.array([lp.eddy_current_heating_halbach(B_surf, g, wavelength, V_REL, mat,
                                                  lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH,
                                                  n_arrays=2) for g in gaps])
    k = 2 * math.pi / wavelength
    F0, Q0 = float(F[0]), float(Q[0])
    sel = gaps < 0.15
    err_F = float(np.max(np.abs(F0 * np.exp(-2 * k * gaps[sel]) - F[sel]) / np.maximum(F[sel], 1e-12)))
    err_Q = float(np.max(np.abs(Q0 * np.exp(-2 * k * gaps[sel]) - Q[sel]) / np.maximum(Q[sel], 1e-12)))
    # energy stored against the limiter from gap g to the nominal gap: U(g) = int_g^GAP0 F dg'
    U = F0 / (2 * k) * (np.exp(-2 * k * gaps) - math.exp(-2 * k * GAP0))
    U = np.maximum(U, 0.0)
    f_eddy = V_REL / wavelength
    skin = math.sqrt(2 * mat["resistivity"] / (2 * math.pi * f_eddy * lc.MU0))
    return dict(gaps=gaps, F=F, Q=Q, U=U, F0=F0, Q0=Q0, k=k, err_F=err_F, err_Q=err_Q,
                B_surf=B_surf, mat=mat, f_eddy=f_eddy, skin=skin, wavelength=wavelength)


def limiter_U(gap, tabs):
    return tabs["F0"] / (2 * tabs["k"]) * max(0.0, math.exp(-2 * tabs["k"] * max(gap, 0.0)) - math.exp(-2 * tabs["k"] * GAP0))


# =============================================================================
# PART C: TIME DOMAIN
# =============================================================================

CONTROL_CASES = [
    ("ideal_871", "ideal feed-forward, ceiling 871 kN/m"),
    ("ideal_218", "ideal feed-forward, ceiling 218 kN/m"),
    ("hold_95",   "EML held at 95.0 kN/m"),
    ("off_sled",  "EML off under the whole sled"),
    ("off_zone",  "EML off under the spacecraft zone only"),
]
CONTROL_NAMES = dict(CONTROL_CASES)


def simulate(u, L_sc, stiff_cases, control_keys, tabs, dx=None, L_win=None, margin_frac=0.1,
             t_settle=T_SETTLE, pd_freq_hz=PD_FREQ_HZ, pd_zeta=1.0, probe=False, record_dt=None,
             dt_override=None, steps_per_pass=12):
    """Integrate the casing through one sled pass at sled speed u (m/s, relative
    to the casing) with the spacecraft spread over L_sc (m).

    Lanes: every (control case, stiffness case) pair is integrated side by side
    on the same grid.  Lane index = i_ctrl * n_stiff + i_stiff.

    Model, per meter of casing, up positive:
      mu w_tt = dV/dx - W_CASING - q_sled(x,t) + F_eml - F_lim,bottom(gap_b) + F_lim,top(gap_t)
      J phi_tt = dM/dx + V
      V = GA (w_x - phi) at half nodes, M = EI phi_x at nodes, free ends.
    Explicit staggered-grid leapfrog.  dt is the smaller of 0.5 dx / c_max and
    L_sc / (steps_per_pass u), so the zone pass is resolved by at least
    steps_per_pass steps.

    Control cases:
      ideal_*  : F_eml = weight + sled load (feed-forward of the known loads) plus a
                 critically damped return to the nominal gap at pd_freq_hz, clipped
                 to [0, ceiling].  Bandwidth is unlimited and the return gain is high
                 enough that the EML uses its full authority whenever the gap or its
                 rate is off; this is the best any controller can do with that force range.
                 It does not fight the limiters or the beam.
      hold_95  : F_eml = W_CASING always.
      off_sled : F_eml = 0 wherever any part of the sled is over the node, W_CASING elsewhere.
      off_zone : F_eml = 0 under the spacecraft zone only.

    Contact: if the bottom (top) gap reaches zero the node stops there with its
    velocity zeroed (dead stop) and the lane is flagged.

    Results are taken over the central part of the window (measurement region),
    away from the free ends.
    """
    n_s, n_c = len(stiff_cases), len(control_keys)
    n_lane = n_s * n_c
    if dx is None:
        dx = L_sc / 8.0
    if L_win is None:
        L_win = L_WIN_SUBSONIC if u <= U_WINDOW_SWITCH else L_WIN_SUPERSONIC
    N = int(round(L_win / dx)) + 1
    x = np.arange(N) * dx
    lo_edge, hi_edge = x - dx / 2, x + dx / 2
    i_lo = int(round(N * (0.5 - margin_frac)))
    i_hi = int(round(N * (0.5 + margin_frac)))
    i_probe = N // 2
    region = slice(i_lo, i_hi + 1)

    EI = np.array([stiff_cases[i % n_s]["EI"] for i in range(n_lane)])[:, None]
    GA = np.array([stiff_cases[i % n_s]["GA"] for i in range(n_lane)])[:, None]
    J = np.array([stiff_cases[i % n_s]["J"] for i in range(n_lane)])[:, None]
    c_max = max(max(s["c_s"], s["c_b"]) for s in stiff_cases)
    ctrl_slices = [(ck, slice(ic * n_s, (ic + 1) * n_s)) for ic, ck in enumerate(control_keys)]

    ad = float(a_down(u))
    q_zone = mass_zone(L_sc) * ad          # N/m, down positive
    q_bare = M_SLED_M * ad

    dt_cfl = 0.5 * dx / c_max if c_max > 0 else np.inf
    dt_pass = L_sc / u / steps_per_pass
    dt = min(dt_cfl, dt_pass, 2e-3)
    if dt_override is not None:
        dt = dt_override
    t_end = (x[i_hi] + L_SLED) / u + t_settle
    n_steps = int(math.ceil(t_end / dt)) + 1

    om = 2 * math.pi * pd_freq_hz
    two_zeta_om = 2 * pd_zeta * om
    om2 = om * om
    F0, k2, Q0 = tabs["F0"], 2 * tabs["k"], tabs["Q0"]
    QF = Q0 / F0

    w = np.zeros((n_lane, N))
    v = np.zeros((n_lane, N))
    phi = np.zeros((n_lane, N - 1))
    phid = np.zeros((n_lane, N - 1))
    Vp = np.zeros((n_lane, N + 1))     # shear with zero padding at the free ends
    Mp = np.zeros((n_lane, N))         # moment with zero at the free ends
    F_eml = np.zeros((n_lane, N))

    min_gap_b = np.full(n_lane, GAP0)
    min_gap_t = np.full(n_lane, GAP0)
    v_peak = np.zeros(n_lane)
    F_lim_peak = np.zeros(n_lane)
    F_eml_peak = np.zeros(n_lane)
    Q_acc = np.zeros((n_lane, i_hi + 1 - i_lo))
    contact_b = np.zeros(n_lane, dtype=bool)
    contact_t = np.zeros(n_lane, dtype=bool)
    active = np.ones(n_lane, dtype=bool)      # records are frozen at a lane's first contact
    w_zone_end = np.full(n_lane, np.nan)
    v_zone_end = np.full(n_lane, np.nan)
    zone_end_done = False
    x_zone_end_probe = x[i_probe] + dx / 2 + L_sc    # sled front position when the zone has fully passed the probe node

    hist_t, hist_w, hist_v, hist_eml = [], [], [], []
    if probe:
        if record_dt is None:
            record_dt = max(dt, 1e-3)
        rec_every = max(1, int(round(record_dt / dt)))

    t0 = time.time()
    for step in range(n_steps):
        t = step * dt
        xf = u * t
        fz = np.clip((np.minimum(hi_edge, xf) - np.maximum(lo_edge, xf - L_sc)) / dx, 0.0, 1.0)
        fb = np.clip((np.minimum(hi_edge, xf - L_sc) - np.maximum(lo_edge, xf - L_SLED)) / dx, 0.0, 1.0)
        f_known = -W_CASING - (q_zone * fz + q_bare * fb)        # weight and sled load, N/m up

        # beam internal forces
        V = GA * ((w[:, 1:] - w[:, :-1]) / dx - phi)
        Vp[:, 1:-1] = V
        f_beam = (Vp[:, 1:] - Vp[:, :-1]) / dx
        Mp[:, 1:-1] = EI * (phi[:, 1:] - phi[:, :-1]) / dx
        a_phi = ((Mp[:, 1:] - Mp[:, :-1]) / dx + V) / J

        # limiters (image-current limit, F0 exp(-2k gap); gap clamped at contact)
        F_lb = F0 * np.exp(-k2 * np.maximum(GAP0 - w, 0.0))
        F_lt = F0 * np.exp(-k2 * np.maximum(GAP0 + w, 0.0))

        # EML
        for ck, sl in ctrl_slices:
            if ck == "ideal_871":
                np.clip(-f_known - MU * (two_zeta_om * v[sl] + om2 * w[sl]), 0.0, F_EML_MAX, out=F_eml[sl])
            elif ck == "ideal_218":
                np.clip(-f_known - MU * (two_zeta_om * v[sl] + om2 * w[sl]), 0.0, F_EML_OP, out=F_eml[sl])
            elif ck == "hold_95":
                F_eml[sl] = W_CASING
            elif ck == "off_sled":
                F_eml[sl] = W_CASING * ((fz + fb) <= 0.0)
            elif ck == "off_zone":
                F_eml[sl] = W_CASING * (fz <= 0.0)

        a = (f_beam + f_known + F_lt - F_lb + F_eml) / MU
        v += a * dt
        w += v * dt
        phid += a_phi * dt
        phi += phid * dt

        # contact
        hb = w >= GAP0
        if hb.any():
            w[hb] = GAP0
            v[hb] = np.minimum(v[hb], 0.0)
            contact_b |= hb[:, region].any(axis=1)
        ht = w <= -GAP0
        if ht.any():
            w[ht] = -GAP0
            v[ht] = np.maximum(v[ht], 0.0)
            contact_t |= ht[:, region].any(axis=1)

        # records over the measurement region, up to each lane's first contact
        wr = w[:, region]
        np.minimum(min_gap_b, GAP0 - wr.max(axis=1), out=min_gap_b)
        np.minimum(min_gap_t, GAP0 + wr.min(axis=1), out=min_gap_t)
        np.maximum(v_peak, active * np.abs(v[:, region]).max(axis=1), out=v_peak)
        np.maximum(F_lim_peak, active * np.maximum(F_lb[:, region], F_lt[:, region]).max(axis=1), out=F_lim_peak)
        np.maximum(F_eml_peak, active * F_eml[:, region].max(axis=1), out=F_eml_peak)
        Q_acc += (QF * dt) * F_lb[:, region] * active[:, None]
        active &= ~(contact_b | contact_t)

        if not zone_end_done and xf >= x_zone_end_probe:
            w_zone_end[:] = w[:, i_probe]
            v_zone_end[:] = v[:, i_probe]
            zone_end_done = True

        if probe and step % rec_every == 0:
            hist_t.append(t)
            hist_w.append(w[:, i_probe].copy())
            hist_v.append(v[:, i_probe].copy())
            hist_eml.append(F_eml[:, i_probe].copy())

    wall = time.time() - t0
    out = []
    for i in range(n_lane):
        s = stiff_cases[i % n_s]
        out.append(dict(
            u=u, L_sc=L_sc, stiff=s["label"], h=s["h"], control=control_keys[i // n_s],
            min_gap_b=float(min_gap_b[i]), min_gap_t=float(min_gap_t[i]),
            v_peak=float(v_peak[i]), F_lim_peak=float(F_lim_peak[i]), F_eml_peak=float(F_eml_peak[i]),
            U_lim=limiter_U(float(min_gap_b[i]), tabs),
            Q_heat=float(Q_acc[i].max()),
            contact_b=bool(contact_b[i]), contact_t=bool(contact_t[i]),
            w_zone_end=float(w_zone_end[i]), v_zone_end=float(v_zone_end[i]),
            dx=dx, dt=dt, n_steps=n_steps, N=N, L_win=L_win, wall=wall,
        ))
    hist = None
    if probe:
        hist = dict(t=np.array(hist_t), w=np.array(hist_w), v=np.array(hist_v), eml=np.array(hist_eml),
                    t_front=x[i_probe] / u, t_zone_end=(x[i_probe] + L_sc) / u,
                    t_tail=(x[i_probe] + L_SLED) / u)
    return out, hist


def _task(args):
    u, L, stiff_cases, ctrl_keys, tabs, kw = args
    res, _ = simulate(u, L, stiff_cases, ctrl_keys, tabs, **kw)
    return res


def run_many(tasks, workers):
    """Run simulate() over a list of (u, L_sc, stiff_cases, ctrl_keys, tabs, kwargs) tasks,
    in parallel processes if workers > 1.  Results come back in task order."""
    if workers <= 1:
        return [_task(t) for t in tasks]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(_task, tasks))


# =============================================================================
# PRINTING HELPERS
# =============================================================================

def hr(title=""):
    print("\n" + "=" * 78)
    if title:
        print(title)
        print("=" * 78)


def fmt_gap(g):
    return f"{g*1000:7.1f}"


def print_results_table(rows, title):
    print(f"\n{title}")
    print(f"  {'u km/s':>7} {'L_sc':>5} {'stiffness':>22} {'control':>10} {'gap_b mm':>9} {'gap_t mm':>9} "
          f"{'v_pk m/s':>9} {'F_lim kN/m':>10} {'F_eml kN/m':>10} {'U kJ/m':>8} {'Q kJ/m':>8} {'hit':>5}")
    for r in rows:
        hit = ("B" if r["contact_b"] else "") + ("T" if r["contact_t"] else "")
        print(f"  {r['u']/1000:7.2f} {r['L_sc']:5.0f} {r['stiff']:>22} {r['control']:>10} "
              f"{fmt_gap(r['min_gap_b']):>9} {fmt_gap(r['min_gap_t']):>9} {r['v_peak']:9.3f} "
              f"{r['F_lim_peak']/1e3:10.1f} {r['F_eml_peak']/1e3:10.1f} {r['U_lim']/1e3:8.2f} {r['Q_heat']/1e3:8.1f} {hit:>5}")


# =============================================================================
# GRAPHS
# =============================================================================

def _save(fig, name, out_dir):
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=GRAPH_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def plot_loads(u_arr, L_scs, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(GRAPH_W, GRAPH_H * 1.4))
    for ax, (umax, ylim) in zip(axes, [(U_MAX, None), (12_000.0, (-900, 900))]):
        sel = u_arr <= umax
        ax.plot(u_arr[sel] / 1e3, -load_bare(u_arr[sel]) / 1e3, color=PALETTE[7], lw=1.5, label="bare sled, 525 kg/m")
        for L, col in zip(L_scs, PALETTE):
            ax.plot(u_arr[sel] / 1e3, -load_zone(u_arr[sel], L) / 1e3, color=col, lw=1.5,
                    label=f"spacecraft zone, L_sc = {L:g} m")
            vn, un = neutral_speed(mass_zone(L))
            if un <= umax:
                ax.axvline(un / 1e3, color=col, ls=":", lw=1)
        ax.axhline(W_CASING / 1e3, color="gray", ls="--", lw=1, label=f"casing weight {W_CASING/1e3:.1f} kN/m")
        ax.axhline(0, color="gray", lw=0.5)
        ax.axvline((V_ORBIT - V_GS) / 1e3, color="gray", ls="-.", lw=1, label="circular velocity")
        ax.set_xlabel("Sled speed relative to the casing (km/s)")
        ax.set_ylabel("Load on the casing, upward (kN/m)")
        if ylim:
            ax.set_ylim(*ylim)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper left")
    _save(fig, "fig_sled_01_load_vs_speed.png", out_dir)


def plot_impulse(u_arr, L_scs, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(GRAPH_W, GRAPH_H * 1.4), sharex=True)
    sel = u_arr >= 500
    for L, col in zip(L_scs, PALETTE):
        I = -impulse_zone(u_arr[sel], L)
        axes[0].plot(u_arr[sel] / 1e3, I / 1e3, color=col, lw=1.5, label=f"spacecraft zone, L_sc = {L:g} m")
        axes[1].plot(u_arr[sel] / 1e3, I / MU, color=col, lw=1.5)
    Ib = -impulse_bare(u_arr[sel])
    axes[0].plot(u_arr[sel] / 1e3, Ib / 1e3, color=PALETTE[7], lw=1.5, label="bare sled, 10 km")
    axes[1].plot(u_arr[sel] / 1e3, Ib / MU, color=PALETTE[7], lw=1.5)
    axes[0].set_ylabel("Impulse per meter of ring, upward (kN s/m)")
    axes[1].set_ylabel("Speed given to a free 12,000 kg/m casing (m/s)")
    axes[1].set_xlabel("Sled speed relative to the casing (km/s)")
    axes[0].set_ylim(-60, 40)
    axes[1].set_ylim(-5, 3.5)
    for ax in axes:
        ax.axhline(0, color="gray", lw=0.5)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    _save(fig, "fig_sled_02_impulse.png", out_dir)


def plot_hypersonic(u_arr, cab, trusses, L_sc, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(GRAPH_W * 1.2, GRAPH_H * 1.4))
    hc = hypersonic_cable(u_arr, cab)
    sel = u_arr > 0
    axes[0, 0].plot(u_arr / 1e3, hc["M"], color=PALETTE[3], lw=1.5, label="cable (load at u + V_rel)")
    axes[0, 1].plot(u_arr / 1e3, hc["T_eff"] / 1e9, color=PALETTE[3], lw=1.5, label="cable, T (M^2 - 1)")
    axes[1, 0].plot(u_arr[sel] / 1e3, hc["L_char"][sel], color=PALETTE[3], lw=1.5, label="cable")
    q = load_zone(u_arr, L_sc)
    axes[1, 1].plot(u_arr[sel] / 1e3, string_term_deflection(q[sel], L_sc, hc["T_eff"][sel]) * 1e3,
                    color=PALETTE[3], lw=1.5, label="cable, if the full zone load reached it")
    for tr, col in zip(trusses, PALETTE):
        hs = hypersonic_casing(u_arr, tr)
        axes[0, 0].plot(u_arr / 1e3, hs["M"], color=col, lw=1.5, label=f"casing, {tr['label']} (u / c_s)")
        axes[1, 0].plot(u_arr[sel] / 1e3, hs["L_char"][sel], color=col, lw=1.5, label=f"casing, {tr['label']}")
    axes[0, 1].plot(u_arr / 1e3, MU * u_arr**2 / 1e9, color=PALETTE[0], lw=1.5, label="casing, mu u^2 (any depth)")
    axes[1, 1].plot(u_arr[sel] / 1e3, string_term_deflection(q[sel], L_sc, MU * u_arr[sel]**2) * 1e3,
                    color=PALETTE[0], lw=1.5, label="casing, any stiffness (string term)")
    axes[0, 0].set_ylabel("Mach number")
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_ylabel("Excess tension term (GN)")
    axes[0, 1].set_yscale("log")
    axes[1, 0].set_ylabel("Characteristic length sqrt(EI / T_eff) (m)")
    axes[1, 0].set_yscale("log")
    axes[1, 1].set_ylabel(f"Deflection at rear of zone, L_sc = {L_sc:g} m (mm)")
    axes[1, 1].set_yscale("log")
    for ax in axes.flat:
        ax.set_xlabel("Sled speed relative to the casing (km/s)")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
    _save(fig, "fig_sled_03_hypersonic.png", out_dir)


LINESTYLES = {"ideal_871": "-", "ideal_218": "--", "hold_95": "-.", "off_sled": ":", "off_zone": (0, (3, 1, 1, 1, 1, 1))}


def plot_min_gap(rows, L_scs, stiff_label, tag, out_dir):
    fig, ax = plt.subplots(figsize=(GRAPH_W, GRAPH_H))
    for L, col in zip(L_scs, PALETTE):
        for ck, _ in CONTROL_CASES:
            sel = sorted([r for r in rows if r["L_sc"] == L and r["stiff"] == stiff_label and r["control"] == ck],
                         key=lambda r: r["u"])
            if not sel:
                continue
            ax.plot([r["u"] / 1e3 for r in sel], [r["min_gap_b"] * 1e3 for r in sel], color=col,
                    ls=LINESTYLES[ck], lw=1.4, label=f"L_sc = {L:g} m, {CONTROL_NAMES[ck]}")
    ax.axhline(GAP_FORBIDDEN * 1e3, color="gray", ls="--", lw=1)
    ax.axhline(0, color="black", lw=1)
    ax.text(0.3, GAP_FORBIDDEN * 1e3 + 1.5, "forbidden zone boundary, 20 mm", fontsize=8, color="gray")
    ax.text(0.3, 1.5, "contact with the cable", fontsize=8)
    ax.set_xlabel("Sled speed relative to the casing (km/s)")
    ax.set_ylabel(f"Smallest bottom gap during the pass (mm), {stiff_label}")
    ax.set_ylim(-3, 105)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.5, ncol=2, loc="upper right")
    _save(fig, f"fig_sled_04_min_gap_{tag}.png", out_dir)


def plot_required_depth(scan_rows, L_scs, depths, out_dir, u_show):
    fig, axes = plt.subplots(1, len(L_scs), figsize=(GRAPH_W * 1.4, GRAPH_H), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, L in zip(axes, L_scs):
        for ck, _ in CONTROL_CASES:
            for u, col in zip(u_show, PALETTE):
                pts = sorted([r for r in scan_rows if r["L_sc"] == L and r["control"] == ck and abs(r["u"] - u) < 1
                              and r["h"] > 0], key=lambda r: r["h"])
                if pts:
                    ax.plot([r["h"] for r in pts], [r["min_gap_b"] * 1e3 for r in pts], color=col, ls=LINESTYLES[ck],
                            marker="o", ms=3, lw=1.2,
                            label=f"{u/1e3:g} km/s, {CONTROL_NAMES[ck]}" if ax is axes[0] else None)
        ax.axhline(GAP_FORBIDDEN * 1e3, color="gray", ls="--", lw=1)
        ax.set_xlabel(f"Truss depth (m), L_sc = {L:g} m")
        ax.set_xscale("log")
        ax.set_xticks(depths)
        ax.set_xticklabels([f"{d:g}" for d in depths])
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("Smallest bottom gap during the pass (mm)")
    axes[0].set_ylim(-3, 105)
    axes[0].legend(fontsize=6, ncol=1, loc="center left")
    _save(fig, "fig_sled_07_gap_vs_truss_depth.png", out_dir)


def plot_time_history(hists, L_sc, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(GRAPH_W, GRAPH_H * 1.4), sharex=False)
    for (u, h, lane_labels), col in zip(hists, PALETTE):
        t_ms = (h["t"] - h["t_front"]) * 1e3
        for (lane, lab, ls) in lane_labels:
            axes[0].plot(t_ms, h["w"][:, lane] * 1e3, color=col, ls=ls, lw=1.3, label=f"{u/1e3:g} km/s, {lab}")
            axes[1].plot(t_ms, h["v"][:, lane], color=col, ls=ls, lw=1.3)
    for ax in axes:
        ax.axhline(GAP0 * 1e3 if ax is axes[0] else 0, color="black", lw=1)
        ax.grid(alpha=0.3)
        ax.set_xlim(-20, 400)
    axes[0].axhline((GAP0 - GAP_FORBIDDEN) * 1e3, color="gray", ls="--", lw=1)
    axes[0].text(2, (GAP0 - GAP_FORBIDDEN) * 1e3 + 2, "bottom forbidden zone boundary", fontsize=8, color="gray")
    axes[0].text(2, GAP0 * 1e3 + 2, "contact with the cable (casing has risen 100 mm)", fontsize=8)
    axes[0].set_ylabel("Casing rise at one point (mm)")
    axes[1].set_ylabel("Casing speed at that point, upward (m/s)")
    axes[1].set_xlabel(f"Time after the sled front passes the point (ms), L_sc = {L_sc:g} m")
    axes[0].legend(fontsize=7, ncol=2)
    _save(fig, "fig_sled_05_time_history.png", out_dir)


def plot_limiter_energy(rows, L_scs, stiff_label, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(GRAPH_W, GRAPH_H * 1.4), sharex=True)
    for L, col in zip(L_scs, PALETTE):
        for ck in ("ideal_871", "hold_95", "off_zone"):
            sel = sorted([r for r in rows if r["L_sc"] == L and r["stiff"] == stiff_label and r["control"] == ck],
                         key=lambda r: r["u"])
            if not sel:
                continue
            axes[0].plot([r["u"] / 1e3 for r in sel], [r["U_lim"] / 1e3 for r in sel], color=col, ls=LINESTYLES[ck],
                         lw=1.3, label=f"L_sc = {L:g} m, {CONTROL_NAMES[ck]}")
            axes[1].plot([r["u"] / 1e3 for r in sel], [r["Q_heat"] / 1e3 for r in sel], color=col, ls=LINESTYLES[ck], lw=1.3)
    axes[0].set_ylabel("Energy stored against the bottom limiter (kJ/m)")
    axes[1].set_ylabel("Reaction plate eddy heat per pass (kJ/m)")
    axes[1].set_xlabel("Sled speed relative to the casing (km/s)")
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)
    _save(fig, "fig_sled_06_limiter_energy.png", out_dir)


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    opts = dict(l_sc=[50.0, 100.0, 200.0], depths=[10.0, 20.0, 30.0], f=0.1, material="CNT",
                graphs=True, quick=False, web_fraction=None, plate=lc.SAFETY_PLATE_MATERIAL,
                workers=max(1, min(6, (os.cpu_count() or 2) - 2)))
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__)
            sys.exit(0)
        elif arg.startswith("--l_sc="):
            opts["l_sc"] = [float(arg.split("=")[1])]
        elif arg.startswith("--truss_depth="):
            opts["depths"] = [float(arg.split("=")[1])]
        elif arg.startswith("--chord_fraction="):
            opts["f"] = float(arg.split("=")[1])
        elif arg.startswith("--web_fraction="):
            opts["web_fraction"] = float(arg.split("=")[1])
        elif arg.startswith("--material="):
            opts["material"] = arg.split("=")[1]
        elif arg.startswith("--workers="):
            opts["workers"] = int(arg.split("=")[1])
        elif arg.startswith("--plate="):
            opts["plate"] = arg.split("=")[1]
        elif arg == "--no-graphs":
            opts["graphs"] = False
        elif arg == "--quick":
            opts["quick"] = True
        else:
            print(f"Unknown option {arg}; see --help")
            sys.exit(1)
    return opts


def main():
    global MU, W_CASING, M_SC
    opts = parse_args()
    L_scs = opts["l_sc"]
    f = opts["f"]
    material = opts["material"]
    web = opts["web_fraction"]
    workers = opts["workers"]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "sled_passage")
    os.makedirs(out_dir, exist_ok=True)
    t_all = time.time()

    hr("SLED PASSAGE LOAD ON THE LEVITATION BEARING")
    print(f"Design point: {ring.ALTITUDE_KM:.0f} km, {'retrograde' if ring.RETROGRADE else 'prograde'} cable")
    print(f"  R_orbit = {R_ORBIT/1e3:,.1f} km   V_orbit = {V_ORBIT:,.1f} m/s   V_ground_sync = {V_GS:,.1f} m/s")
    print(f"  g_local = {G_LOCAL:.3f} m/s^2   g_net on casing = {G_NET:.3f} m/s^2   casing weight = {W_CASING/1e3:.1f} kN/m")
    print(f"  cable: {MU_CABLE:,.0f} kg/m with hardware, tension {T_CABLE/1e9:,.1f} GN, side {CABLE_SIDE:.2f} m, "
          f"V_cable = {V_CABLE:,.1f} m/s, V_rel = {V_REL:,.1f} m/s")
    print(f"  sled: {L_SLED:,.0f} m at {M_SLED_M:.0f} kg/m, spacecraft {M_SC/1e3:,.0f} t, zone lengths {L_scs} m")
    print(f"  EML: {F_EML_OP/1e3:.0f} kN/m at 400 A, {F_EML_MAX/1e3:.0f} kN/m stated maximum, pull only, two tracks under the cable")
    print(f"  Sled speed u is relative to the casing; inertial speed is u + {V_GS:.0f} m/s.  Top speed {U_MAX/1e3:.0f} km/s.")
    print(f"  Rails on top of the casing; spacecraft zone at the front of the sled.")
    print(f"  Chord material for the time domain: {MATERIALS[material]['label']}; chord fraction f = {f:g}; "
          f"web fraction = {('f/2' if web is None else str(web))}")
    print(f"Output directory: {out_dir}")
    print(f"Workers: {workers}")

    # ------------------------------------------------------------------ Part A
    hr("PART A: LOAD ON THE CASING AGAINST SLED SPEED (closed form)")
    u_tab = np.array([0, 1000, 2000, 3000, 4000, 5000, 6000, V_ORBIT - V_GS, 7000, 8000, 9000, 10000, 12000,
                      15000, 20000, 25000, 30000])
    print(f"\n  Load per meter on the casing, positive down (kN/m).  v_i = u + {V_GS:.0f} m/s.")
    print(f"  {'u km/s':>7} {'v_i km/s':>9} {'a_down':>8} {'g':>7} {'bare':>9} " +
          " ".join(f"{'zone '+str(int(L))+' m':>12}" for L in L_scs))
    for u in u_tab:
        ad = float(a_down(u))
        print(f"  {u/1e3:7.3f} {(u+V_GS)/1e3:9.3f} {ad:8.3f} {-ad/G0_STD:7.2f} {load_bare(u)/1e3:9.1f} " +
              " ".join(f"{load_zone(u, L)/1e3:12.1f}" for L in L_scs))
    print(f"\n  Mass per meter in the spacecraft zone: " +
          ", ".join(f"L_sc = {L:g} m: {mass_zone(L):,.0f} kg/m" for L in L_scs))
    print(f"  Load is zero at the circular velocity: u = {V_ORBIT - V_GS:,.0f} m/s (v_i = {V_ORBIT:,.0f} m/s).")
    print(f"\n  Speed at which the load cancels the casing weight ({W_CASING/1e3:.1f} kN/m), net bearing load reverses:")
    vn_b, un_b = neutral_speed(M_SLED_M)
    print(f"    bare sled alone:  v_i = {vn_b:,.0f} m/s, u = {un_b:,.0f} m/s  (above the {U_MAX/1e3:.0f} km/s run)")
    neutral = {}
    for L in L_scs:
        vn, un = neutral_speed(mass_zone(L))
        neutral[L] = un
        print(f"    zone, L_sc = {L:5g} m:  v_i = {vn:,.0f} m/s, u = {un:,.0f} m/s")

    print(f"\n  Passage times at a fixed point and impulse per meter of ring (zone and bare sled), positive down:")
    print(f"  {'u km/s':>7} " + " ".join(f"{'t_zone '+str(int(L)):>11}" for L in L_scs) + f" {'t_sled s':>9} " +
          " ".join(f"{'I_zone '+str(int(L)):>11}" for L in L_scs) + f" {'I_bare':>9} " +
          " ".join(f"{'v_free '+str(int(L)):>11}" for L in L_scs) + f" {'v_bare':>8}")
    print(f"  {'':>7} " + " ".join(f"{'(ms)':>11}" for L in L_scs) + f" {'':>9} " +
          " ".join(f"{'(kN s/m)':>11}" for L in L_scs) + f" {'(kN s/m)':>9} " +
          " ".join(f"{'(m/s)':>11}" for L in L_scs) + f" {'(m/s)':>8}")
    for u in u_tab[1:]:
        Iz = [float(impulse_zone(u, L)) for L in L_scs]
        Ib = float(impulse_bare(u))
        print(f"  {u/1e3:7.3f} " + " ".join(f"{L/u*1e3:11.2f}" for L in L_scs) + f" {L_SLED/u:9.3f} " +
              " ".join(f"{I/1e3:11.2f}" for I in Iz) + f" {Ib/1e3:9.2f} " +
              " ".join(f"{-I/MU:11.3f}" for I in Iz) + f" {-Ib/MU:8.3f}")
    print(f"\n  The zone impulse per meter is M_spacecraft * a / u plus the sled hardware share, so it hardly")
    print(f"  depends on L_sc.  The free casing speed after the zone passes is that impulse over {MU:,.0f} kg/m.")

    # ------------------------------------------------------------ hand checks
    hr("HAND VALUES OF THE INSTRUCTIONS (L_sc = 200 m, u = 30 km/s)")
    L_h, u_h = 200.0, 30_000.0
    a_frame = -float(a_down(u_h))                   # upward, casing-frame convention
    a_lsm = -float(a_down_inertial(u_h))            # upward, lsm_simulation convention (30 km/s inertial)
    tp = L_h / u_h
    acc_frame = (mass_zone(L_h) * a_frame - W_CASING) / MU
    acc_lsm = (mass_zone(L_h) * a_lsm - W_CASING) / MU
    hand = [
        ("Mass per meter in the zone (kg/m)", 25_525, mass_zone(L_h), mass_zone(L_h)),
        ("Peak radial load (g)", 12.2, a_frame / G0_STD, a_lsm / G0_STD),
        ("Load under the zone, upward (kN/m)", 3050, mass_zone(L_h) * a_frame / 1e3, mass_zone(L_h) * a_lsm / 1e3),
        ("Load under the bare sled, upward (kN/m)", 62.8, M_SLED_M * a_frame / 1e3, M_SLED_M * a_lsm / 1e3),
        ("Neutral speed, zone, inertial (km/s)", 9.1, neutral_speed(mass_zone(L_h))[0] / 1e3, neutral_speed(mass_zone(L_h))[0] / 1e3),
        ("Zone passage time (ms)", 6.7, tp * 1e3, tp * 1e3),
        ("Impulse per meter, zone (kN s/m)", 19.9, mass_zone(L_h) * a_frame * tp / 1e3, mass_zone(L_h) * a_lsm * tp / 1e3),
        ("Casing displacement during the pass, zero stiffness, EML off (mm)", 5.5, 0.5 * acc_frame * tp**2 * 1e3, 0.5 * acc_lsm * tp**2 * 1e3),
        ("Casing speed after the pass (m/s)", 1.6, acc_frame * tp, acc_lsm * tp),
        ("Casing kinetic energy per meter (kJ/m)", 16, 0.5 * MU * (acc_frame * tp)**2 / 1e3, 0.5 * MU * (acc_lsm * tp)**2 / 1e3),
    ]
    print(f"\n  {'Quantity':64} {'hand':>9} {'script':>9} {'lsm conv':>9} {'check':>6}")
    print(f"  {'':64} {'':>9} {'(u rel.)':>9} {'(v inert)':>9}")
    for name, hv, sv, lv in hand:
        rel = abs(sv - hv) / abs(hv)
        rel_l = abs(lv - hv) / abs(hv)
        chk = "pass" if rel <= 0.05 else ("lsm" if rel_l <= 0.05 else "FAIL")
        print(f"  {name:64} {hv:9.4g} {sv:9.4g} {lv:9.4g} {chk:>6}")
    print(f"\n  'script' uses the casing-frame convention (v_i = {u_h + V_GS:,.0f} m/s).  'lsm conv' uses")
    print(f"  v_i = {u_h:,.0f} m/s, which is what lsm_simulation.py does: it starts the sled at v = 0 and")
    print(f"  puts that v straight into v^2/R, so its radial load treats the launch speed as inertial.")
    print(f"  The hand values follow the lsm number (12.2 g).  In the casing-frame convention the peak")
    print(f"  radial load is {a_frame/G0_STD:.2f} g and every load-based hand value is {100*(a_frame/a_lsm-1):.1f} % higher.")
    print(f"  'pass' = within 5 % of the script value; 'lsm' = within 5 % only under the lsm convention.")
    print(f"  The displacement, speed and energy lines include the casing weight acting during the pass (EML off).")

    # ------------------------------------------------------------------ Part B
    hr("PART B: HYPERSONIC EXCESS VELOCITY EQUATIONS")
    cab = cable_properties()
    print(f"\n  Cable: side {CABLE_SIDE:.2f} m, I = {cab['I']:.1f} m^4, EI = {cab['EI']:.3e} N m^2 (E = {E_CNT/1e9:.0f} GPa)")
    print(f"  Transverse wave speed sqrt(T/mu): {cab['v_wave']:,.0f} m/s with hardware mass "
          f"({cab['v_wave_struct']:,.0f} m/s on structural mass only, the anchor_line_analysis.py form).")
    print(f"  Static L_char = sqrt(EI/T) = {math.sqrt(cab['EI']/T_CABLE):.1f} m (Chapter 3, Eq. 3.10).")
    print(f"  Note: anchor_line_analysis.py computes L_char as (EI/T_eff)^(1/4), which is not dimensionally a")
    print(f"  length, and uses E = 25 GPa there; this script uses sqrt(EI/T_eff) with E = 500 GPa as Chapter 3 does.")
    print(f"\n  Frame moving with the load: EI W'''' + (mu v^2 - T) W'' = q.  The excess term mu v^2 - T = T (M^2 - 1)")
    print(f"  carries the sign of a compressive axial force.  It does not stiffen the beam against the load.")
    print(f"  For M > 1 nothing moves ahead of the load; behind it the term produces a trailing wake, and")
    print(f"  the mean of that wake leaves each point moving at F/(mu v) * M^2/(M^2 - 1), which is the")
    print(f"  zero-stiffness impulse result with a small tension correction.  Stiffness sets the wake's")
    print(f"  wavelength 2 pi L_char, not the momentum left behind.")

    u_fine = np.linspace(0, U_MAX, 601)
    trusses_B = [truss(h, f, material, web) for h in opts["depths"]]
    print(f"\n  {'u km/s':>7} | {'cable: v_load':>13} {'Mach':>6} {'T_eff GN':>10} {'L_char m':>9} | " +
          " ".join(f"{'casing h='+str(int(tr['h'])):>14} {'L_char':>7}" for tr in trusses_B))
    for u in [0, 1000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000]:
        hc = hypersonic_cable(u, cab)
        line = f"  {u/1e3:7.1f} | {float(hc['v']):13,.0f} {float(hc['M']):6.2f} {float(hc['T_eff'])/1e9:10,.0f} {float(hc['L_char']):9.2f} | "
        for tr in trusses_B:
            hs = hypersonic_casing(u, tr)
            line += f"{'M='+format(float(hs['M']), '.2f'):>14} {float(hs['L_char']):7.2f} "
        print(line)

    print(f"\n  Truss casing table.  EI = E A_chords (h/2)^2, A_chords = f * {MU:,.0f} / rho; shear from diagonals at 45 deg")
    print(f"  carrying half the chord mass fraction (web fraction = f/2); c_s = sqrt(GA/mu), c_b = sqrt(EI/J) = sqrt(E/rho).")
    print(f"  {'material':>8} {'h m':>5} {'f':>5} {'EI N m^2':>11} {'GA N':>11} {'J kg m':>10} {'c_s m/s':>9} {'c_b m/s':>9} "
          f"{'sqrt(EI/mu)':>12} {'L_char@30km/s':>14}")
    for mat in ("CNT", "steel"):
        for h in (10.0, 20.0, 30.0):
            for ff in (0.1, 0.2, 0.3):
                tr = truss(h, ff, mat, None if web is None else web)
                print(f"  {mat:>8} {h:5.0f} {ff:5.2f} {tr['EI']:11.3e} {tr['GA']:11.3e} {tr['J']:10.3e} {tr['c_s']:9,.0f} "
                      f"{tr['c_b']:9,.0f} {math.sqrt(tr['EI']/MU):12,.0f} {math.sqrt(tr['EI']/MU)/U_MAX:14.2f}")
    print(f"\n  Every shear and bar wave speed in the table is below the sled speed over most of the run, so")
    print(f"  the casing ahead of the sled gets no warning.  What stiffness can do is spread the impulse")
    print(f"  over casing behind and beside the load.  Part C measures that instead of assuming it.")

    for L in L_scs:
        print(f"\n  String-term deflection at the rear of the zone, L_sc = {L:g} m: q L_sc^2 / (2 T_eff)")
        print(f"  {'u km/s':>7} {'q_zone kN/m':>12} {'cable mm':>9} {'casing mm':>10}   (cable: if the full zone load reached it)")
        for u in [1000, 5000, 10000, 20000, 30000]:
            hc = hypersonic_cable(u, cab)
            q = float(load_zone(u, L))
            print(f"  {u/1e3:7.1f} {q/1e3:12.1f} {float(string_term_deflection(q, L, hc['T_eff']))*1e3:9.3f} "
                  f"{float(string_term_deflection(q, L, MU*u**2))*1e3:10.3f}")

    # ------------------------------------------------------------------ limiters
    hr("LIMITER FORCE LAW (levitation_physics.safety_repulsive_force, V_rel from ring_altitude)")
    lc.SAFETY_PLATE_MATERIAL = opts["plate"]
    tabs = limiter_tables()
    print(f"  Plate: {tabs['mat']['name']}, {lc.SAFETY_PLATE_THICKNESS*1e3:.0f} mm thick, {lc.SAFETY_PLATE_WIDTH:.1f} m wide, 2 arrays per side")
    print(f"  Halbach surface field {tabs['B_surf']:.3f} T at {lc.SAFETY_I_PERSISTENT} A, wavelength {tabs['wavelength']*1e3:.0f} mm")
    print(f"  Eddy frequency V_rel / wavelength = {tabs['f_eddy']/1e3:.1f} kHz, skin depth {tabs['skin']*1e3:.2f} mm")
    print(f"  Image-current limit holds; F = {tabs['F0']/1e3:.1f} kN/m * exp(-2 k gap), Q = {tabs['Q0']/1e6:.2f} MW/m * exp(-2 k gap),")
    print(f"  k = {tabs['k']:.1f} /m.  Largest deviation from the tabulated functions: F {tabs['err_F']:.1e}, Q {tabs['err_Q']:.1e}.")
    print(f"  Note: levitation_simulation.py hard codes 8,620 m/s for V_rel in four places (lines 273, 510, 649, 785);")
    print(f"  this script uses ring_altitude.V_REL = {V_REL:,.0f} m/s.  The force is in the image-current limit either way.")
    print(f"  {'gap mm':>7} {'F kN/m':>8} {'U kJ/m':>8} {'Q MW/m':>8}")
    for g in (0.0, 0.002, 0.005, 0.010, 0.020, 0.030, 0.050, 0.080, 0.100):
        print(f"  {g*1e3:7.1f} {np.interp(g, tabs['gaps'], tabs['F'])/1e3:8.1f} {limiter_U(g, tabs)/1e3:8.2f} "
              f"{np.interp(g, tabs['gaps'], tabs['Q'])/1e6:8.3f}")
    print(f"  U is the energy stored pushing the casing from the nominal gap down to that gap (int F dg).")
    print(f"  The whole limiter travel to contact stores {limiter_U(0.0, tabs)/1e3:.2f} kJ/m.")
    print("")
    print(f"  Reaction plate materials at contact (gap 0), same force law; Chapter 8 (8.4a) names molybdenum, levitation_config.py")
    print(f"  defaults to copper, and the time domain above uses {tabs['mat']['name']}:")
    print(f"  {'material':>28} {'F(0) kN/m':>10} {'Q(0) MW/m':>10} {'skin mm':>8} {'skin layer kg/m':>16} {'cp J/kgK':>9}")
    for name in ("copper_OFHC", "molybdenum", "aluminum_6061", "gamma_TiAl"):
        mat = lp.get_material(name)
        sig = 1.0 / mat["resistivity"]
        Fm = lp.safety_repulsive_force(tabs['B_surf'], 0.0, tabs['wavelength'], sig, lc.SAFETY_PLATE_THICKNESS,
                                       lc.SAFETY_PLATE_WIDTH, V_REL, n_arrays=2)
        Qm = lp.eddy_current_heating_halbach(tabs['B_surf'], 0.0, tabs['wavelength'], V_REL, mat,
                                             lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, n_arrays=2)
        sk = math.sqrt(2 * mat["resistivity"] / (2 * math.pi * tabs['f_eddy'] * lc.MU0))
        print(f"  {mat['name']:>28} {Fm/1e3:10.1f} {Qm/1e6:10.2f} {sk*1e3:8.2f} {2*lc.SAFETY_PLATE_WIDTH*sk*mat['density']:16.2f} {mat['specific_heat']:9.0f}")
    print(f"  The skin layer is the plate volume within one skin depth on both bottom arrays, the mass that takes the heat")
    print(f"  of one pass before conduction spreads it.")

    # ------------------------------------------------------------------ Part C
    hr("PART C: TIME DOMAIN")
    trusses = [truss(h, f, material, web) for h in opts["depths"]]
    stiff_cases = trusses + [ZERO_STIFFNESS]
    ctrl_keys = [k for k, _ in CONTROL_CASES]
    print(f"  Stiffness cases: " + "; ".join(s["label"] for s in stiff_cases))
    print(f"  Control cases: " + "; ".join(f"{k} = {n}" for k, n in CONTROL_CASES))
    print(f"  Ideal controller: feed-forward of weight and sled load plus a critically damped return to the")
    print(f"  nominal gap at {PD_FREQ_HZ:.0f} Hz, clipped to the force range.  The control law is not designed; this is")
    print(f"  the best any controller can do with the stated force range and no bandwidth limit.")
    print(f"  Window {L_WIN_SUBSONIC/1e3:.0f} km for u <= {U_WINDOW_SWITCH/1e3:.0f} km/s, {L_WIN_SUPERSONIC/1e3:.0f} km above; free ends; "
          f"results from the central 20 %.  Cable fixed.  Contact modeled as a dead stop.")
    print(f"  dx = L_sc / 8; dt = min(0.5 dx / c_max, L_sc / (12 u)).")

    speeds = list(range(1000, int(U_MAX) + 1, 3000 if opts["quick"] else 1000))
    speeds = sorted(set([float(s) for s in speeds] + [round(neutral[L], 1) for L in L_scs]))
    tasks = [(u, L, stiff_cases, ctrl_keys, tabs, {}) for L in L_scs for u in speeds]
    tasks.sort(key=lambda t: (t[0], -t[1]))     # slowest (low u, short zone) first for load balance
    t_start = time.time()
    results = run_many(tasks, workers)
    rows = [r for res in results for r in res]
    for res in results:
        r0 = res[0]
        print(f"  ran u = {r0['u']/1e3:6.2f} km/s, L_sc = {r0['L_sc']:4g} m: window {r0['L_win']/1e3:.0f} km, dx = {r0['dx']:.2f} m, "
              f"dt = {r0['dt']*1e6:.0f} us, {r0['n_steps']:,} steps, {r0['wall']:.1f} s")
    print(f"  Sweep wall time {time.time()-t_start:.0f} s")

    for s in stiff_cases:
        print_results_table(sorted([r for r in rows if r["stiff"] == s["label"]], key=lambda r: (r["L_sc"], r["u"], ctrl_keys.index(r["control"]))),
                            f"Results, {s['label']} (gap_b, gap_t: smallest bottom and top gap; v_pk: peak casing speed; "
                            f"F_lim: peak limiter force; F_eml: peak EML force; U: energy stored against the bottom limiter at its "
                            f"smallest gap; Q: plate eddy heat per pass; hit B/T: contact at bottom/top; peaks and Q are "
                            f"taken up to the first contact)")

    # hand check of the zero-stiffness, EML-off-under-zone lane
    hr("TIME DOMAIN CHECK AGAINST THE CLOSED FORM (zero stiffness, EML off under the zone)")
    print(f"  Rise and upward speed of the probe point at the instant the zone has passed it, against 1/2 a t^2 and a t")
    print(f"  with a = (q_zone - weight)/mu and t = L_sc/u.")
    for L in L_scs:
        for u in (10_000.0, 20_000.0, 30_000.0):
            r = [r for r in rows if r["L_sc"] == L and abs(r["u"] - u) < 1 and r["stiff"] == "zero stiffness" and r["control"] == "off_zone"]
            if r:
                r = r[0]
                q = float(load_zone(u, L))
                acc = (-q - W_CASING) / MU
                tp = L / u
                print(f"  L_sc = {L:4g} m, u = {u/1e3:4.0f} km/s: rise {r['w_zone_end']*1e3:7.3f} mm (closed form {0.5*acc*tp**2*1e3:7.3f}), "
                      f"speed {r['v_zone_end']:6.3f} m/s (closed form {acc*tp:6.3f}), smallest bottom gap {r['min_gap_b']*1e3:6.1f} mm, "
                      f"{'contact' if r['contact_b'] else 'no contact'}")

    # summary: first failing speed per (L_sc, stiffness, control)
    hr("SUMMARY: SLED SPEED AT WHICH THE BOTTOM GAP FIRST DROPS BELOW 20 mm, AND AT WHICH THE CASING FIRST TOUCHES")
    print(f"  {'L_sc':>5} {'stiffness':>22} {'control':>10} {'<20 mm at':>10} {'contact at':>11} {'gap_b 30 km/s':>14} {'worst gap_b':>12} {'at u':>7}  "
          f"{'top <20 mm at':>13} {'top hit at':>11}")
    for L in L_scs:
        for s in stiff_cases:
            for ck in ctrl_keys:
                sel = sorted([r for r in rows if r["L_sc"] == L and r["stiff"] == s["label"] and r["control"] == ck], key=lambda r: r["u"])
                u_20 = next((r["u"] for r in sel if r["min_gap_b"] < GAP_FORBIDDEN), None)
                u_hit = next((r["u"] for r in sel if r["contact_b"]), None)
                u_t20 = [r["u"] for r in sel if r["min_gap_t"] < GAP_FORBIDDEN]
                u_thit = [r["u"] for r in sel if r["contact_t"]]
                r30 = next((r for r in sel if abs(r["u"] - U_MAX) < 1), None)
                worst = min(sel, key=lambda r: r["min_gap_b"])
                g30 = f"{r30['min_gap_b']*1e3:.1f} mm" if r30 else "n/a"
                rng = lambda lst: "none" if not lst else (f"{min(lst)/1e3:.0f} to {max(lst)/1e3:.0f} km/s" if len(lst) > 1 else f"{lst[0]/1e3:.0f} km/s")
                print(f"  {L:5g} {s['label']:>22} {ck:>10} {('none' if u_20 is None else f'{u_20/1e3:.0f} km/s'):>10} "
                      f"{('none' if u_hit is None else f'{u_hit/1e3:.0f} km/s'):>11} {g30:>14} {worst['min_gap_b']*1e3:9.1f} mm {worst['u']/1e3:5.0f}  "
                      f"{rng(u_t20):>13} {rng(u_thit):>11}")

    # ------------------------------------------------------------ depth scan
    hr("REQUIRED TRUSS DEPTH: SCAN OF DEPTH AT SELECTED SPEEDS (%s chords, f = %.2f)" % (material, f))
    depths_scan = [10.0, 20.0, 30.0, 50.0, 100.0]
    u_scan = [5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0] if not opts["quick"] else [10000.0, 20000.0, 30000.0]
    scan_stiff = [truss(h, f, material, web) for h in depths_scan] + [ZERO_STIFFNESS]
    tasks = [(u, L, scan_stiff, ctrl_keys, tabs, {}) for L in L_scs for u in u_scan]
    tasks.sort(key=lambda t: (t[0], -t[1]))
    t_start = time.time()
    scan_rows = [r for res in run_many(tasks, workers) for r in res]
    print(f"  Scan wall time {time.time()-t_start:.0f} s")
    print(f"\n  Smallest bottom gap (mm) against truss depth; 'hit' marks contact.")
    for L in L_scs:
        for ck in ctrl_keys:
            print(f"\n  L_sc = {L:g} m, {CONTROL_NAMES[ck]}")
            print(f"  {'u km/s':>7} " + " ".join(f"{s['label'][:14]:>15}" for s in scan_stiff))
            for u in u_scan:
                line = f"  {u/1e3:7.0f} "
                for s in scan_stiff:
                    r = next(r for r in scan_rows if r["L_sc"] == L and abs(r["u"] - u) < 1 and r["stiff"] == s["label"] and r["control"] == ck)
                    line += f"{r['min_gap_b']*1e3:11.1f}{' hit' if r['contact_b'] else '    '}"
                print(line)
    print(f"\n  Smallest depth that keeps the bottom gap above 20 mm at every scanned speed:")
    for L in L_scs:
        for ck in ctrl_keys:
            req = None
            for h in depths_scan:
                ok = all(r["min_gap_b"] >= GAP_FORBIDDEN for r in scan_rows
                         if r["L_sc"] == L and r["control"] == ck and r["h"] == h)
                if ok:
                    req = h
                    break
            g100 = next(r["min_gap_b"] for r in scan_rows if r["L_sc"] == L and r["control"] == ck and r["h"] == 100.0 and abs(r["u"] - U_MAX) < 1)
            g0 = next(r["min_gap_b"] for r in scan_rows if r["L_sc"] == L and r["control"] == ck and r["h"] == 0.0 and abs(r["u"] - U_MAX) < 1)
            print(f"    L_sc = {L:4g} m, {CONTROL_NAMES[ck]:42}: "
                  f"{('none up to 100 m' if req is None else f'{req:.0f} m')};  at 30 km/s the gap is {g100*1e3:.1f} mm with a 100 m truss "
                  f"and {g0*1e3:.1f} mm with zero stiffness")
    print(f"\n  If the 100 m truss and zero stiffness give the same gap, the beam is not the lever: behind a load")
    print(f"  that outruns every wave in the casing, every meter has received the same impulse, the motion is")
    print(f"  uniform along the ring, and a uniform rise carries no shear or bending for the truss to resist.")
    print(f"  What limits the gap then is the force balance per meter: {W_CASING/1e3:.1f} kN/m of weight against the")
    print(f"  kick, with the EML unable to push and the limiters storing {limiter_U(0.0, tabs)/1e3:.2f} kJ/m at most.")

    # ------------------------------------------------------------ convergence
    Lc = L_scs[min(1, len(L_scs) - 1)]
    hr("CONVERGENCE CHECK (u = 30 km/s, L_sc = %g m, %s, ideal_871 and off_zone)" % (Lc, trusses[0]["label"]))
    conv_stiff = [trusses[0], ZERO_STIFFNESS]
    print(f"  {'dx m':>7} {'dt us':>7} {'steps':>8} {'stiffness':>22} {'control':>10} {'gap_b mm':>9} {'gap_t mm':>9} {'v_pk m/s':>9} "
          f"{'v_zone_end':>10} {'U kJ/m':>8} {'Q kJ/m':>8}")
    for div in (4, 8, 16, 32):
        res, _ = simulate(U_MAX, Lc, conv_stiff, ["ideal_871", "off_zone"], tabs, dx=Lc / div)
        for r in res:
            print(f"  {r['dx']:7.2f} {r['dt']*1e6:7.1f} {r['n_steps']:8,} {r['stiff']:>22} {r['control']:>10} {fmt_gap(r['min_gap_b']):>9} "
                  f"{fmt_gap(r['min_gap_t']):>9} {r['v_peak']:9.3f} {r['v_zone_end']:10.3f} {r['U_lim']/1e3:8.2f} {r['Q_heat']/1e3:8.1f}")
    print(f"  Time step check at dx = L_sc/8 (dt halved and quartered):")
    for fac in (1.0, 0.5, 0.25):
        res0, _ = simulate(U_MAX, Lc, conv_stiff, ["ideal_871", "off_zone"], tabs)
        dt0 = res0[0]["dt"]
        res, _ = simulate(U_MAX, Lc, conv_stiff, ["ideal_871", "off_zone"], tabs, dt_override=dt0 * fac)
        for r in res:
            print(f"  {r['dx']:7.2f} {r['dt']*1e6:7.1f} {r['n_steps']:8,} {r['stiff']:>22} {r['control']:>10} {fmt_gap(r['min_gap_b']):>9} "
                  f"{fmt_gap(r['min_gap_t']):>9} {r['v_peak']:9.3f} {r['v_zone_end']:10.3f} {r['U_lim']/1e3:8.2f} {r['Q_heat']/1e3:8.1f}")
    print(f"  Window length check (free ends), L_sc = {Lc:g} m, {trusses[0]['label']}:")
    for uw in (3000.0, 10000.0):
        for Lw in (6e3, 12e3, 20e3):
            res, _ = simulate(uw, Lc, [trusses[0]], ["ideal_218", "hold_95"], tabs, L_win=Lw)
            for r in res:
                print(f"    u = {uw/1e3:4.0f} km/s, window {Lw/1e3:4.0f} km, {r['control']:>10}: gap_b {fmt_gap(r['min_gap_b'])} mm, "
                      f"gap_t {fmt_gap(r['min_gap_t'])} mm, v_pk {r['v_peak']:.3f} m/s")

    # ------------------------------------------------------------ probe histories for graph 5
    hists = []
    if opts["graphs"]:
        Lh = 100.0 if 100.0 in L_scs else L_scs[0]
        for u in (5000.0, 10000.0, 20000.0, 30000.0):
            res, h = simulate(u, Lh, [trusses[0], ZERO_STIFFNESS], ["ideal_871", "off_zone"], tabs, probe=True,
                              record_dt=2e-4, t_settle=1.0)
            lanes = [(0, f"{trusses[0]['label']}, ideal 871", "-"), (3, "zero stiffness, EML off under zone", ":")]
            hists.append((u, h, lanes))

    # ------------------------------------------------------------ what would change the result
    Lo = 100.0 if 100.0 in L_scs else L_scs[0]
    hr("WHAT WOULD CHANGE THE RESULT (bisection at u = 30 km/s, L_sc = %g m, ideal_871, zero stiffness)" % Lo)
    print(f"  Zero stiffness is used because the depth scan shows the truss depth does not change the gap at this speed.")
    print(f"  Each line is the value at which the smallest bottom gap just reaches the target.")

    def bisect(fn, lo, hi, target, increasing, n=14):
        for _ in range(n):
            mid = 0.5 * (lo + hi)
            ok = fn(mid) >= target
            if increasing == ok:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def gap_at(u, L, tabs_=None):
        r, _ = simulate(u, L, [ZERO_STIFFNESS], ["ideal_871"], tabs_ or tabs)
        return r[0]["min_gap_b"]

    t_start = time.time()
    for target, name in ((GAP_FORBIDDEN, "gap >= 20 mm"), (0.0005, "no contact")):
        u_lim = bisect(lambda u: gap_at(u, Lo), 10_000.0, U_MAX, target, increasing=False)
        print(f"  Top sled speed, design as is, {name:12}: {u_lim/1e3:5.1f} km/s")
    for target, name in ((GAP_FORBIDDEN, "gap >= 20 mm"), (0.0005, "no contact")):
        def with_F0(m):
            t = dict(tabs)
            t["F0"] = tabs["F0"] * m
            t["Q0"] = tabs["Q0"] * m
            return gap_at(U_MAX, Lo, t)
        m = bisect(with_F0, 1.0, 400.0, target, increasing=True)
        print(f"  Limiter force at 30 km/s, {name:12}: {m:6.1f} times the present law, {tabs['F0']*m/1e3:,.0f} kN/m at contact, "
              f"surface field x{math.sqrt(m):.2f} ({tabs['B_surf']*math.sqrt(m):.2f} T), persistent current x{math.sqrt(m):.2f} "
              f"({lc.SAFETY_I_PERSISTENT*math.sqrt(m):,.0f} A at the present turns)")
    MU0_, W0_ = MU, W_CASING
    for target, name in ((GAP_FORBIDDEN, "gap >= 20 mm"), (0.0005, "no contact")):
        def with_mu(mu):
            global MU, W_CASING
            MU, W_CASING = mu, mu * G_NET
            try:
                return gap_at(U_MAX, Lo)
            finally:
                MU, W_CASING = MU0_, W0_
        mu_req = bisect(with_mu, MU0_, 10 * MU0_, target, increasing=True)
        print(f"  Casing mass along the mass driver at 30 km/s, {name:12}: {mu_req:,.0f} kg/m ({mu_req/MU0_:.2f} times), "
              f"weight {mu_req*G_NET/1e3:.0f} kN/m for the EML to hold")
    MSC0_ = M_SC
    for target, name in ((GAP_FORBIDDEN, "gap >= 20 mm"), (0.0005, "no contact")):
        def with_msc(m):
            global M_SC
            M_SC = m
            try:
                return gap_at(U_MAX, Lo)
            finally:
                M_SC = MSC0_
        msc = bisect(with_msc, 1e5, MSC0_, target, increasing=False)
        print(f"  Spacecraft mass at 30 km/s, {name:12}: {msc/1e3:,.0f} t (of {MSC0_/1e3:,.0f} t)")
    print(f"  Options wall time {time.time()-t_start:.0f} s")
    print(f"  Not run: limiter wavelength (a longer wavelength reaches further at the cost of drag and heating at the")
    print(f"  nominal gap), an EML that can push (upper tracks are out of the design, Paul 2026-10-09), and a longer")
    print(f"  spacecraft load length, which the sweep shows does not change the outcome above 10 km/s.")

    # ------------------------------------------------------------ global ring note
    hr("SIDE NOTE: OUTWARD MOMENTUM GIVEN TO THE WHOLE RING PER LAP [UNVERIFIED, outside this task]")
    print(f"  Whatever the bearing does, the sled's net outward force ends up in the ring.  Per lap at speed u the")
    print(f"  ring (cable + casing, {MU_CABLE + MU:,.0f} kg/m) receives (F_zone + F_bare)/u per meter.")
    EA = E_CNT * ring.CABLE_AREA
    om_breath = math.sqrt(EA / (MU_CABLE + MU)) / R_ORBIT
    print(f"  Hoop breathing mode from EA/(mu R^2) alone (no gravity gradient term): period {2*math.pi/om_breath/60:.1f} min.")
    print(f"  {'u km/s':>7} {'F_total MN':>11} {'impulse kN s/m':>15} {'dv ring m/s':>12} {'excursion m':>12}")
    for u in (10000.0, 15000.0, 20000.0, 25000.0, 30000.0):
        F_tot = -(float(load_zone(u, L_scs[0])) * L_scs[0] + float(load_bare(u)) * (L_SLED - L_scs[0]))
        I = F_tot / u
        dv = I / (MU_CABLE + MU)
        print(f"  {u/1e3:7.0f} {F_tot/1e6:11.1f} {I/1e3:15.2f} {dv:12.3f} {dv/om_breath:12.1f}")

    # ------------------------------------------------------------ graphs
    if opts["graphs"]:
        hr("GRAPHS")
        plot_loads(u_fine, L_scs, out_dir)
        plot_impulse(u_fine, L_scs, out_dir)
        plot_hypersonic(u_fine, cab, trusses, 100.0 if 100.0 in L_scs else L_scs[0], out_dir)
        for s in stiff_cases:
            tag = f"h{s['h']:.0f}" if s["h"] > 0 else "zero"
            plot_min_gap(rows, L_scs, s["label"], tag, out_dir)
        plot_required_depth(scan_rows, L_scs, depths_scan, out_dir, [10000.0, 20000.0, 30000.0])
        plot_time_history(hists, 100.0 if 100.0 in L_scs else L_scs[0], out_dir)
        plot_limiter_energy(rows, L_scs, trusses[0]["label"], out_dir)

    # ------------------------------------------------------------ CSV
    csv_path = os.path.join(out_dir, "sled_passage_results.csv")
    with open(csv_path, "w") as fh:
        fh.write("u_m_s,L_sc_m,stiffness,truss_depth_m,control,min_gap_bottom_m,min_gap_top_m,peak_casing_speed_m_s,"
                 "peak_limiter_force_N_m,peak_eml_force_N_m,limiter_energy_J_m,plate_heat_J_m,contact_bottom,contact_top,"
                 "rise_at_zone_end_m,speed_at_zone_end_m_s\n")
        for r in rows + scan_rows:
            fh.write(f"{r['u']:.1f},{r['L_sc']:.0f},{r['stiff']},{r['h']:.0f},{r['control']},{r['min_gap_b']:.6f},{r['min_gap_t']:.6f},"
                     f"{r['v_peak']:.4f},{r['F_lim_peak']:.1f},{r['F_eml_peak']:.1f},{r['U_lim']:.1f},{r['Q_heat']:.1f},"
                     f"{int(r['contact_b'])},{int(r['contact_t'])},{r['w_zone_end']:.6f},{r['v_zone_end']:.4f}\n")
    print(f"\nSaved: {csv_path}")
    hr(f"DONE in {time.time()-t_all:.0f} s")


if __name__ == "__main__":
    main()
