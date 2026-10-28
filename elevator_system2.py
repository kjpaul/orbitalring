#!/usr/bin/env python3
"""
The elevator system, second pass: the track's width and conductance, the track as hot as it can go, where the cruise energy
goes, 4,000 m/s and the stop acceleration, the magnet arrangement, throughput against pod spacing, and the safety case.
Book: The Orbital Ring, Volume IV: Chapter 17.
Run: python elevator_system2.py > elevator_system2.txt    (after radiation/elevator_dose_rates.py; imports elevator_system.py)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring); ORBITAL_RING_OUT (output folder, default runs)

Paul's decisions since the first pass (task file of 2026-10-28, second pass): the brake track is a consumable carbon part whose
conductivity and thickness are chosen for eddy braking, run as hot as it can be; the structural yarn never runs hot; faster is
better; the objective is throughput per tonne of tether; a pod that cannot continue is ejected; the tether is a flat ribbon.
Sections follow the task file's items 1 to 7. Every ASSUMPTION is printed as such.
"""
import os, sys, math, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import elevator_system as e1
import freefall_trajectory as ff
import tether_config as tcfg
import ring_altitude as ring

OUT = e1.OUT; FIG = os.path.join(OUT, "elevator_system2"); os.makedirs(FIG, exist_ok=True)
BLUE, RED, GREEN, GREY = e1.BLUE, e1.RED, e1.GREEN, e1.GREY
GM, RE, W, r_r, r_geo, g0, MU0, SB, RHO, SIG = e1.GM, e1.RE, e1.W, e1.r_r, e1.r_geo, e1.g0, e1.MU0, e1.SB, e1.RHO, e1.SIG
anet, Phi = e1.anet, e1.Phi
section = e1.section

# ---- JANAF graphite extended to 2,500 K (janaf.nist.gov/tables/C-002.html): T, H - H(298.15) kJ/mol ----
e1.JANAF_T = np.array([298.15, 300, 500, 800, 1000, 1200, 1500, 1600, 1700, 1800, 1900, 2000, 2500])
e1.JANAF_H = np.array([0.0, 0.016, 2.365, 7.637, 11.795, 16.240, 23.253, 25.658, 28.090, 30.547, 33.026, 35.525, 48.289]) * 1e3 / 12.011
T0 = e1.T_TRACK0
DH = {T: e1.h_graphite(T) - e1.h_graphite(T0) for T in (1000.0, 1500.0, 2000.0)}
EPS_TRACK = e1.EPS_TRACK

# ---- the track (item 1) ----
W_MAG = 1.0                             # ribbon width = magnet width (Paul: 1 m, or what the carriage needs)
SIG_MIN, SIG_MAX = 0.05e6, 2.0e6        # what carbon offers for the ribbon's conductivity (Paul: 0.05 to 2 MS/m)
V_HAND = e1.V_HAND

# ---- the magnets (item 5) ----
BR, RHO_M = 1.35, 7500.0                # NdFeB; SmCo 1.10 T, 8,400 kg/m3 as the alternative
BR_SMCO, RHO_SMCO = 1.10, 8400.0
GAP = e1.GAP; T_RIBBON_NOM = 0.001      # the two arrays are 2 x 20 mm plus the ribbon apart
B_SAT_FE, RHO_FE = 1.6, 7870.0          # ASSUMPTION: iron yoke worked to 1.6 T

# ---- heat into the carriage (item 2) ----
EPS_SHIELD = (0.03, 0.1)                # multilayer reflective shield on the carriage face, effective emissivity (Paul: 0.03 to 0.1)
T_W_LIMIT, T_W0, C_W = e1.T_W_LIMIT, e1.T_W0, e1.C_W
L_BOIL = 2.256e6                        # latent heat of water, J/kg
Q_RAD = 0.9 * SB * (373.0**4 - 250.0**4)  # ASSUMPTION: condenser radiator at 100 C, emissivity 0.9, sink 250 K (earthshine and sun not resolved)

# ---- regeneration (item 3) ----
SIG_CU = 60e6                           # CNT at the conductivity of copper (Chapter 10)
F_LOSS = 0.10                           # I^2R loss in the winding (task)
V_LINE = (100e3, 1e6)                   # bus voltage on the ribbon's conductor: tether_config V_BUS and the 1 MV of Chapter 15's lines
TAPE_PER_M2 = 1630.0 / (3 * 100.0 * 0.844)   # LIM stator tape per m2 of band: 1,630 kg per LIM over 3 pole pitches of 100 m and the 0.844 m band (Chapter 11, 2026-10-27)
CRYO_PER_M2 = 5.0                       # ASSUMPTION: a vacuum-insulated cryostat skin 5 kg/m2
LEAK_PER_M2 = 6200.0 / (450.0 * 0.844)  # W/m2: the site's 6.2 kW leak over its 450 m band (Chapter 18)
COP = 0.0106                            # cryocooler COP at 70 K (Chapter 18)

# ---- throughput (item 6) ----
SPACINGS = (300.0, 900.0, 1800.0, 3600.0)      # 5 min added beyond the task's 15, 30, 60

# ---- safety (item 7) ----
ISP = (230.0, 300.0)
BUDGETS = (0.05, 0.10)

def fmt(x, nd=0): return f"{x:,.{nd}f}"

# ======================================================================================================================
# item 5: magnet arrangements (needed first: the pods are sized on the lightest)
# ======================================================================================================================
def arrangements():
    rows = []
    D = 2 * GAP + T_RIBBON_NOM
    for tau in (0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3):
        for h in (0.025, 0.05, 0.075, 0.1, 0.15, 0.2):
            B1 = e1.halbach_B(BR, tau, h, GAP)                   # one array at the ribbon
            # (a) single-sided pair, fields counted separately (first pass): drag per m2 of ribbon 2 x B1^2/(4 mu0)
            rows.append(dict(kind="single-sided pair", tau=tau, h=h, B=B1, p=2 * B1**2 / (4 * MU0), m_mag=2 * RHO_M * h, m_fe=0.0))
            # (b) facing pair, flux of one closing through the other: the fields add at the ribbon
            B2 = 2 * e1.halbach_B(BR, tau, h, D / 2)
            rows.append(dict(kind="facing pair", tau=tau, h=h, B=B2, p=B2**2 / (4 * MU0), m_mag=2 * RHO_M * h, m_fe=0.0))
            # (c) iron-backed C core: magnets h thick on both poles in series with the gap D, yokes carry B tau / 2 at B_SAT
            Bc = BR * 2 * h / (2 * h + D)
            t_y = Bc * tau / (2 * B_SAT_FE)
            rows.append(dict(kind="iron C core", tau=tau, h=h, B=Bc, p=Bc**2 / (4 * MU0), m_mag=2 * RHO_M * h, m_fe=2 * RHO_FE * t_y))
    for r in rows:
        r["kg_per_N"] = (r["m_mag"] + r["m_fe"]) / r["p"]; r["N_per_m"] = r["p"] * W_MAG
    return rows

# ======================================================================================================================
# item 1: the track
# ======================================================================================================================
def track2(m, v_c, dH, p_peak, a_stop=e1.A_STOP, rg=None, w=W_MAG, free_fall=False, v_entry=None, rule="heat", L=None):
    """ribbon of width w. Thickness from the heat of a pass under the magnets, t = F/(rho w dH); conductivity to keep the
    characteristic speed at the pod speed, sigma = 2/(mu0 t v), clipped to what carbon offers (SIG_MIN to SIG_MAX).
    rule="strict": where sigma would exceed SIG_MAX the thickness follows the conductance at SIG_MAX (the ribbon is then
    thicker than the heat needs and runs cooler). rule="heat": the thickness stays at the heat value, the conductivity at
    SIG_MAX, the characteristic speed is then above the pod speed (u = v/w_char < 1) and the drag is 2u/(1+u^2) of its peak;
    that is allowed wherever the force needed, F, is within that fraction of the carriage's peak force p_peak w L, and
    where it is not the thickness rises to the value the force needs. The carriage is sized for the force at the ring, so
    the stop's bottom always needs u = 1 and is conductance-set; the cruise needs a third of that force or less.
    free_fall: no braking outside the stop, which starts from v_entry."""
    rg = e1.r_grid() if rg is None else rg
    s = rg - r_r; a = anet(rg)
    if free_fall:
        s_stop = v_entry**2 / (2 * a_stop); in_stop = rg < r_r + s_stop
        v = np.where(in_stop, np.sqrt(np.maximum(2 * a_stop * s, 0.0)), 0.0); v = np.maximum(v, V_HAND)
        F = np.where(in_stop, m * (a + a_stop), 0.0)
    else:
        in_stop = rg < r_r + v_c**2 / (2 * a_stop)
        v = np.where(in_stop, np.sqrt(np.maximum(2 * a_stop * s, 0.0)), v_c); v = np.maximum(v, V_HAND)
        F = np.where(in_stop, m * (a + a_stop), m * a)
    F_max = m * (anet(r_r) + a_stop); L = F_max / (p_peak * w) if L is None else L
    t_heat = F / (RHO * w * dH)
    t_cond = 2.0 / (MU0 * SIG_MAX * v)                                   # the thickness for u = 1 at SIG_MAX
    f = np.clip(F / (p_peak * w * L), 0.0, 1.0)                          # share of the carriage's peak force needed here
    u_min = np.where(f > 1e-9, (1 - np.sqrt(np.maximum(1 - f**2, 0.0))) / np.maximum(f, 1e-9), 0.0)   # 2u/(1+u^2) = f, small root
    t_force = u_min * t_cond
    if rule == "strict":
        sig_req = np.where(t_heat > 0, 2.0 / (MU0 * np.maximum(t_heat, 1e-12) * v), np.inf)
        t = np.where(sig_req > SIG_MAX, t_cond, t_heat)
    else:
        t = np.maximum(t_heat, t_force)
        sig_req = np.where(t > 0, 2.0 / (MU0 * np.maximum(t, 1e-12) * v), np.inf)
    t = np.where(F > 0, t, 0.0)
    sig = np.clip(sig_req, SIG_MIN, SIG_MAX)
    w_char = 2.0 / (MU0 * sig * np.maximum(t, 1e-12)); u = np.where(t > 0, v / w_char, 0.0)
    drag_frac = 2 * u / (1 + u**2)
    mu = RHO * t * w; dH_local = np.where(mu > 0, F / np.maximum(mu, 1e-12), 0.0)
    w_force = F / (p_peak * np.maximum(drag_frac, 1e-9))
    short = (F > 0) & (drag_frac * p_peak * w * L < F * 0.999)            # the force cannot be met (only where sigma is clipped at SIG_MIN)
    return dict(r=rg, alt=rg - RE, s=s, v=v, F=F, t=t, t_heat=t_heat, t_force=t_force, t_cond=t_cond, sig_req=sig_req, sig=sig, u=u, drag_frac=drag_frac, mu=mu,
                dH_local=dH_local, in_stop=in_stop, w=np.full(len(rg), w), w_force=w_force, below_min=(sig_req < SIG_MIN) & (F > 0), above_max=(sig_req > SIG_MAX) & (F > 0), short=short, f=f, L=L)

def integrate(tr, mask=None):
    mu = tr["mu"] if mask is None else np.where(mask, tr["mu"], 0.0)
    return float(np.trapezoid(mu, tr["r"]))

def cut_track(tr, res):
    own = RHO * res["A"]; above = np.where(own >= tr["mu"] / e1.TETHER_SHARE)[0]
    return tr["r"][above[0]] if len(above) else r_geo

def size_with_track(d, tr, rg):
    """tether with the pods on the line and the track, the track cut where the tether's own mass is 10 x the track's (two passes)"""
    _, _, E_tot, Ew_tot = e1.pod_envelope(d, rg); E = np.full(len(rg), E_tot)
    res = e1.size_tether(rg, tr["mu"], E)
    for _ in range(2):
        r_cut = cut_track(tr, res); mu_cut = np.where(rg < r_cut, tr["mu"], 0.0); res = e1.size_tether(rg, mu_cut, E)
    res.update(r_cut=r_cut, m_track=float(np.trapezoid(mu_cut, rg)), E_tot=E_tot, Ew_tot=Ew_tot, mu_cut=mu_cut)
    return res

# ======================================================================================================================
# pods on the facing pair
# ======================================================================================================================
def size_pod2(N, x_shield, mag, a_stop=e1.A_STOP):
    m = 15e3 * N / 25
    for _ in range(80):
        F_max = m * (a_stop + anet(r_r)); m_mag = F_max * mag["kg_per_N"]
        L = max(e1.L_CAR, F_max / (mag["p"] * W_MAG))
        g = e1.pod_geometry(N, x_shield, m_mag, W_MAG, mag["h"]); m_new = g["m"]
        if abs(m_new - m) < 0.1: break
        m = m_new
    g.update(F_max=F_max, A_face=2 * L * W_MAG, L_car=L, a_stop=a_stop)
    return g

# ======================================================================================================================
# the free-fall-then-stop profile (item 3, last option)
# ======================================================================================================================
def freefall_profile(m, a_stop=e1.A_STOP, dt=2.0, v0=10.0):
    """released just below the platform with v0 (ASSUMPTION: a push to 10 m/s), no braking until the 1.5 g stop that ends at the ring"""
    # speed at radius r in free fall from GEO: v^2 = v0^2 + 2 (Phi(r_geo) - Phi(r)); the stop needs v^2 = 2 a_stop (r - r_r)
    lo, hi = r_r, r_geo
    for _ in range(200):
        mid = 0.5 * (lo + hi); vff2 = v0**2 + 2 * (Phi(r_geo) - Phi(mid)); vst2 = 2 * a_stop * (mid - r_r)
        if vff2 > vst2: lo = mid
        else: hi = mid
    r_stop = hi; v_entry = math.sqrt(2 * a_stop * (r_stop - r_r))
    T_, R_, V_, F_, P_, PH = [], [], [], [], [], []
    r = r_geo; v = v0; t = 0.0
    while r > r_stop:
        v += anet(r) * dt; r -= v * dt; t += dt
        T_.append(t); R_.append(r); V_.append(v); F_.append(0.0); P_.append(0.0); PH.append(1)
    t_fall = t
    while r > r_r and v > 0:
        dts = min(dt, 0.5); v -= a_stop * dts; r -= v * dts; t += dts; F = m * (anet(r) + a_stop)
        T_.append(t); R_.append(r); V_.append(max(v, 0.0)); F_.append(F); P_.append(F * max(v, 0.0)); PH.append(2)
    d = {k: np.array(x) for k, x in dict(t=T_, r=R_, v=V_, F=F_, P=P_, phase=PH).items()}
    d["alt"] = d["r"] - RE; d["s_km"] = (r_geo - d["r"]) / 1e3; d["dt"] = None; d["m"] = m; d["v_c"] = v_entry
    d["r_stop"] = r_stop; d["s_stop"] = r_stop - r_r; d["total"] = t; d["t_fall"] = t_fall; d["v_entry"] = v_entry
    d["P_stop_peak"] = float(d["P"].max()); d["F_max"] = m * (anet(r_r) + a_stop)
    d["t_below_6000"] = float(np.sum(np.diff(np.concatenate([[0.0], d["t"]]))[d["alt"] < 6000e3]))
    d["n_on_line"] = int(math.ceil(t / e1.HEADWAY))
    return d

def dose_profile(rates, d, key="max_+0"):
    """dose per descent for a profile with variable steps"""
    dts = np.diff(np.concatenate([[0.0], d["t"]]))
    sec, _ = np.histogram(d["alt"], bins=rates["edges"], weights=dts)
    p = (rates[f"Hp_{key}"] * sec[None, :]).sum(1) * 1e3; e = (rates[f"De_{key}"] * sec[None, :]).sum(1) * 1e3
    return p, e, p + e

# ======================================================================================================================
# item 2: heat into the carriage and the water, with the shield
# ======================================================================================================================
def heat_into_water(d, tr, A_face, eps_face):
    e1.EPS_MAG = eps_face
    Q, q_peak = e1.magnet_heat(d, tr, A_face)
    return Q, q_peak

def cooling_curve(t_mm, T_hot, seconds):
    """ribbon temperature after cooling for the given seconds from T_hot, both faces radiating at EPS_TRACK to 300 K surroundings"""
    T = T_hot; tt = 0.0; dt = 0.02
    while tt < seconds and T > 301.0:
        cp = (e1.h_graphite(T + 1) - e1.h_graphite(T - 1)) / 2.0
        T -= 2 * EPS_TRACK * SB * (T**4 - 300.0**4) / (RHO * t_mm * 1e-3 * cp) * dt; tt += dt
    return T

def sublimation_per_pass(t_mm, T_hot):
    """kg/m2 lost from both faces while the ribbon cools from T_hot, on the graphite rate fit 4.6e30 exp(-81,200/T) atoms/cm2/s
    (Guhathakurta and Draine 1989, as quoted; UNVERIFIED here), unit sticking coefficient"""
    T = T_hot; tt = 0.0; dt = 0.02; lost = 0.0
    while T > 1200.0:
        flux = 4.6e30 * math.exp(-81200.0 / T) * 1e4          # atoms/m2/s
        lost += 2 * flux * 12.011 * 1.6605e-27 * dt
        cp = (e1.h_graphite(T + 1) - e1.h_graphite(T - 1)) / 2.0
        T -= 2 * EPS_TRACK * SB * (T**4 - 300.0**4) / (RHO * t_mm * 1e-3 * cp) * dt; tt += dt
    return lost, tt

# ======================================================================================================================
# item 3: regenerating stator
# ======================================================================================================================
def regen_masses(d, tr, r_top, V, dH, superconducting=False):
    """mass per metre of a regenerating stator on the ribbon from the ring up to r_top. Resistive: the winding takes F_LOSS of
    the brake energy of each pass as heat within the rise dH (so its mass is F_LOSS x the passive heat-sized track), and the
    ribbon's own CNT conductor carries the power to the ring at V with a loss of F_LOSS over the distance, sized for the largest
    power that passes each point. Superconducting: tape and cryostat per m2 of the LIM's band, the same conductor, and the
    cryogenic plant's power as a note. Above r_top the passive ribbon of tr."""
    rg = tr["r"]; in_sec = rg < r_top
    P_at = np.interp(rg, d["r"][::-1], d["P"][::-1])                     # brake power when the pod is at r
    P_env = np.maximum.accumulate(P_at[::-1])[::-1]                     # the largest power generated at or above r passes r
    A_cond = P_env * (rg - r_r) / (SIG_CU * F_LOSS * V**2)              # conductor section for F_LOSS over the run to the ring
    mu_cond = np.where(in_sec, RHO * A_cond, 0.0)
    if superconducting:
        mu_wind = np.where(in_sec, (TAPE_PER_M2 + CRYO_PER_M2) * W_MAG, 0.0)
        leak_W = float(np.trapezoid(np.where(in_sec, LEAK_PER_M2 * W_MAG, 0.0), rg))
    else:
        mu_wind = np.where(in_sec, F_LOSS * tr["F"] / dH, 0.0); leak_W = 0.0
    mu_pass = np.where(in_sec, 0.0, tr["mu"])
    return dict(mu=mu_wind + mu_cond + mu_pass, mu_wind=mu_wind, mu_cond=mu_cond, mu_pass=mu_pass, leak_W=leak_W, P_peak=float(P_env.max()), A_cond_max=float(A_cond.max()))

# ======================================================================================================================
# item 7: free fall with an initial velocity, and the thrust
# ======================================================================================================================
def fall(r0, vr0, vt0, R_nose, beta):
    state = np.array([r0, 0.0, vr0, vt0]); t = 0.0
    eps = 0.5 * (vr0**2 + vt0**2) - GM / r0; hh = r0 * vt0
    if eps < 0:
        a = -GM / (2 * eps); e = math.sqrt(max(1 - hh**2 / (GM * a), 0.0)); r_p = a * (1 - e)
        if r_p > RE + 100e3: return dict(orbit=True, perigee=r_p - RE, apogee=a * (1 + e) - RE)
    q_tot = 0.0; q_peak = 0.0; g_peak = 0.0; v_entry = None; gam = None; t20 = 0.0; t40 = 0.0; t_peak = 0.0
    for _ in range(4_000_000):
        alt = state[0] - RE
        dt = 2.0 if alt > 400e3 else (0.2 if alt > 150e3 else 0.02)
        state = ff.rk4_step(state, dt, beta); t += dt
        r, th, vr, vt = state; alt = max(r - RE, 0.0)
        vt_rel = vt - W * r; v_rel = math.sqrt(vr**2 + vt_rel**2)
        rho = ff.atmosphere_density(alt)
        q = tcfg.K_SUTTON_GRAVES * math.sqrt(rho / R_nose) * v_rel**3; q_tot += q * dt; q_peak = max(q_peak, q)
        gg = 0.5 * rho * v_rel**2 / beta / g0
        if gg > g_peak: g_peak = gg; t_peak = t
        t20 += dt * (gg > 20); t40 += dt * (gg > 40)
        if v_entry is None and alt <= 100e3: v_entry = v_rel; gam = math.degrees(math.atan2(-vr, abs(vt_rel)))
        if alt <= 0: break
    return dict(orbit=False, v_entry=v_entry, gamma=gam, q_peak=q_peak, q_total=q_tot, g_peak=g_peak, t20=t20, t40=t40, t_total=t, v_ground=v_rel)

def perigee_of(r, vr, vt):
    eps = 0.5 * (vr**2 + vt**2) - GM / r; hh = r * vt
    if eps >= 0: return float("inf")
    a = -GM / (2 * eps); e = math.sqrt(max(1 - hh**2 / (GM * a), 0.0)); return a * (1 - e)

def min_dv(r, vr, vt, r_p_target):
    """smallest single impulse that lifts the perigee to r_p_target: search over direction and magnitude"""
    if perigee_of(r, vr, vt) >= r_p_target: return 0.0, (0.0, 0.0)
    best = (float("inf"), None)
    for ang in np.linspace(0, 2 * math.pi, 181):
        ux, uy = math.cos(ang), math.sin(ang)                              # (radial, tangential) direction
        lo, hi = 0.0, 20000.0
        if perigee_of(r, vr + hi * ux, vt + hi * uy) < r_p_target: continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if perigee_of(r, vr + mid * ux, vt + mid * uy) >= r_p_target: hi = mid
            else: lo = mid
        if hi < best[0]: best = (hi, (ux, uy))
    return best

def scale_height():
    _, _, rho30 = tcfg.atmosphere(30e3); _, _, rho60 = tcfg.atmosphere(60e3)
    return 30e3 / math.log(rho30 / rho60)

# ======================================================================================================================
# main
# ======================================================================================================================
def main():
    print("elevator_system2.py: second pass on the elevator system")
    print(f"JANAF graphite enthalpy from {T0:.0f} K: to 1,000 K {DH[1000.0]/1e3:.0f} kJ/kg, to 1,500 K {DH[1500.0]/1e3:,.0f}, to 2,000 K {DH[2000.0]/1e3:,.0f} kJ/kg ; ribbon emissivity {EPS_TRACK} [UNVERIFIED at these temperatures]")
    rates = e1.load_rates(); T_sh = rates["T"]; rg = e1.r_grid()

    # ---------------- item 5: magnets ----------------
    section("5. MAGNET ARRANGEMENTS ON THE RIBBON (first, because the pods are sized on the lightest)")
    print(f"ribbon under a millimetre; the two arrays {2*GAP*1e3+T_RIBBON_NOM*1e3:.0f} mm apart. NdFeB {BR} T, {RHO_M:.0f} kg/m3. Drag per m2 of ribbon at peak B^2/(4 mu0) with B the field in the ribbon.")
    print("(a) single-sided pair as run: each array's field counted alone, drag 2 x B1^2/(4 mu0). (b) facing pair with the flux of one closing through the other: the fields add in the")
    print(f"ribbon, B = 2 B1 at the midplane, drag (2 B1)^2/(4 mu0). (c) iron C core: magnets h thick on both poles in series with the {(2*GAP+T_RIBBON_NOM)*1e3:.0f} mm gap, B = Br 2h/(2h + gap),")
    print(f"yokes on both sides carrying B tau/2 at {B_SAT_FE} T (ASSUMPTION), iron {RHO_FE:.0f} kg/m3, leakage ignored (favours the core).")
    rows = arrangements()
    print("\n  arrangement | pitch mm | magnet mm | B in ribbon T | peak drag kPa | magnet kg/m2 | iron kg/m2 | kg per MN | kN per metre of carriage at 1 m width")
    for kind in ("single-sided pair", "facing pair", "iron C core"):
        sub = [r for r in rows if r["kind"] == kind]
        for r in sub:
            if r["h"] in (0.05, 0.1) and r["tau"] in (0.1, 0.125, 0.2):
                print(f"  {kind:18s} | {r['tau']*1e3:8.0f} | {r['h']*1e3:9.0f} | {r['B']:13.3f} | {r['p']/1e3:13.1f} | {r['m_mag']:12.0f} | {r['m_fe']:10.0f} | {r['kg_per_N']*1e6:9,.0f} | {r['N_per_m']/1e3:.0f}")
        b = min(sub, key=lambda r: r["kg_per_N"])
        print(f"  lightest {kind}: pitch {b['tau']*1e3:.0f} mm, magnet {b['h']*1e3:.0f} mm, {b['B']:.3f} T, {b['p']/1e3:.1f} kPa, {b['kg_per_N']*1e6:,.0f} kg per MN, {b['N_per_m']/1e3:.0f} kN per metre of carriage")
    best = {k: min([r for r in rows if r["kind"] == k], key=lambda r: r["kg_per_N"]) for k in ("single-sided pair", "facing pair", "iron C core")}
    mag = dict(best["facing pair"]); mag["name"] = "NdFeB facing pair"
    smco = dict(mag); smco["B"] = mag["B"] * BR_SMCO / BR; smco["p"] = smco["B"]**2 / (4 * MU0); smco["m_mag"] = 2 * RHO_SMCO * mag["h"]; smco["kg_per_N"] = smco["m_mag"] / smco["p"]
    print(f"\ncarried forward: the facing pair, {mag['B']:.3f} T in the ribbon, {mag['p']/1e3:.0f} kPa, {mag['kg_per_N']*1e6:,.0f} kg per MN ({best['single-sided pair']['kg_per_N']/mag['kg_per_N']:.2f} x lighter than the single-sided pair, "
          f"{best['iron C core']['kg_per_N']/mag['kg_per_N']:.2f} x lighter than the iron core). SmCo facing pair: {smco['B']:.3f} T, {smco['p']/1e3:.0f} kPa, {smco['kg_per_N']*1e6:,.0f} kg per MN.")

    # ---------------- descents, doses, pods for the speeds ----------------
    section("DESCENTS, DOSES AND PODS AT 2,000, 3,000 AND 4,000 m/s (1 g drive, 1.5 g stop; 2 mSv shield; facing-pair magnets)")
    SPEEDS = (2000.0, 3000.0, 4000.0)
    prof = {v: e1.descent(v, 1000.0) for v in SPEEDS}; shield = {}; dose = {}
    print("  m/s | trip h | below 6,000 km min | drive min, km | stop min, km | peak brake MW per tonne | dose mSv at 10 / 20 / 30 / 50 g/cm2 | shield for 5 mSv g/cm2 | for 2 mSv")
    for v in SPEEDS:
        d = prof[v]; p, e, tot = e1.dose_per_descent(rates, d); dose[v] = tot
        for target in (5.0, 2.0): shield[(v, target)] = e1.shield_for(T_sh, tot, target)
        print(f"  {v:5,.0f} | {d['total']/3600:6.2f} | {d['t_below_6000']/60:18.0f} | {d['t_drive']/60:.1f}, {d['s_drive']/1e3:,.0f} | {d['t_stop']/60:.1f}, {d['s_stop']/1e3:,.0f} | {d['P_stop_peak']/1e6:23.0f} | "
              f"{tot[6]:.1f} / {tot[8]:.1f} / {tot[9]:.2f} / {tot[11]:.2f} | {shield[(v,5.0)]:22.1f} | {shield[(v,2.0)]:.1f}")
    pods = {}
    print("\n  m/s | seats | pod t | shield t | magnets t | carriage t | outer dia m | brake force at ring kN | carriage m | face m2 | seats per tonne")
    for v in SPEEDS:
        for N in (100, 400):
            g = size_pod2(N, shield[(v, 2.0)], mag); pods[(N, v)] = g
            print(f"  {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {g['m_sh']/1e3:8.1f} | {g['m_mag']/1e3:9.1f} | {g['m_car']/1e3:10.1f} | {2*g['R_out']:11.2f} | {g['F_max']/1e3:22,.0f} | {g['L_car']:10.1f} | {g['A_face']:7.1f} | {N/(g['m']/1e3):.2f}")
    print("(first pass, single-sided pair: 100 seats 137.5 t and 400 seats 304.5 t at 2,000 m/s)")

    # ---------------- item 1: the track ----------------
    section("1. THE TRACK: width fixed at the magnets, thickness from the heat of a pass, conductivity to keep w = v")
    print(f"ribbon {W_MAG:.1f} m wide; t = F/(rho w dH); sigma = 2/(mu0 t v), carbon offers {SIG_MIN/1e6:.2f} to {SIG_MAX/1e6:.0f} MS/m.")
    print("STRICT rule (the task's): where sigma would have to exceed 2 MS/m the thickness follows the conductance at 2 MS/m and the ribbon runs cooler and heavier.")
    print("HEAT rule (used from here on): the thickness stays at the heat value with sigma at 2 MS/m, so the characteristic speed is above the pod speed and the drag is")
    print("2u/(1+u^2) of its peak; that is allowed where the force needed is within that fraction of the carriage's peak force (the cruise needs a third or less), and")
    print("where it is not (the bottom of the stop, where the whole force is needed) the thickness rises to what the force needs, which is the conductance value.")
    tracks = {}
    g = pods[(100, 2000.0)]
    for rule in ("strict", "heat"):
        for Tmax in (1000.0, 2000.0):
            tr = track2(g["m"], 2000.0, DH[Tmax], mag["p"], rg=rg, rule=rule)
            print(f"\n{rule.upper()} rule, 100 seats ({g['m']/1e3:.1f} t), 2,000 m/s, ribbon to {Tmax:.0f} K: km above ring | speed m/s | force kN | share of peak force | heat thickness mm | thickness used mm | sigma asked MS/m | sigma used | u = v/w_char | drag share | rise K | kg/m")
            for s_km in (0.34, 1, 3, 10, 30, 100, 136, 137, 500, 1000, 2000, 3000, 5000, 10000, 20000):
                j = min(np.searchsorted(rg, r_r + s_km * 1e3), len(rg) - 1)
                print(f"  {s_km:8,.2f} | {tr['v'][j]:7,.0f} | {tr['F'][j]/1e3:8,.0f} | {tr['f'][j]:6.3f} | {tr['t_heat'][j]*1e3:8.3f} | {tr['t'][j]*1e3:8.3f} | {min(tr['sig_req'][j],1e9)/1e6:9.3f} | {tr['sig'][j]/1e6:6.3f} | {tr['u'][j]:5.3f} | {tr['drag_frac'][j]:5.3f} | {e1.T_from_dH(tr['dH_local'][j])-T0:5.0f} | {tr['mu'][j]:.4f}")
            s_hi = rg[tr["above_max"] & ~tr["in_stop"]]; s_lo_stop = rg[tr["above_max"] & tr["in_stop"]]
            print(f"  sigma asked > 2 MS/m: the bottom {(s_lo_stop.max()-r_r)/1e3 if len(s_lo_stop) else 0:.1f} km of the stop, and the cruise above {(s_hi.min()-r_r)/1e3 if len(s_hi) else float('nan'):,.0f} km; "
                  f"below 0.05 MS/m: {int(tr['below_min'].sum())} grid points; force not met: {int(tr['short'].sum())} grid points; "
                  f"track to the platform {integrate(tr)/1e3:,.0f} t (stop {integrate(tr, tr['in_stop'])/1e3:,.0f}, cruise {integrate(tr, ~tr['in_stop'])/1e3:,.0f})")
    print("\ntrack masses on the new rule, cut where the tether's own mass is 10 x the track's (the tether of item 4 below):")
    print("  m/s | seats | pod t | ribbon K | stop section t | cruise section to the platform t | total t | x pod | below the cut t (cut km) | heat lower bound m dPhi/dH t")
    results = {}
    for v in SPEEDS:
        for N in (100, 400):
            g = pods[(N, v)]; d = e1.descent(v, g["m"])
            for Tmax in (1000.0, 1500.0, 2000.0):
                tr = track2(g["m"], v, DH[Tmax], mag["p"], rg=rg); tracks[(N, v, Tmax)] = tr
                res = size_with_track(d, tr, rg); res.update(m_pod=g["m"], P_peak=d["P_stop_peak"], n_on_line=d["n_on_line"], d=d); results[(N, v, Tmax)] = res
                m_stop = integrate(tr, tr["in_stop"]); m_cr = integrate(tr, ~tr["in_stop"])
                print(f"  {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {Tmax:8.0f} | {m_stop/1e3:14,.0f} | {m_cr/1e3:32,.0f} | {(m_stop+m_cr)/1e3:7,.0f} | {(m_stop+m_cr)/g['m']:5.1f} | {res['m_track']/1e3:,.0f} ({(res['r_cut']-RE)/1e3-ring.ALTITUDE_KM:,.0f}) | {g['m']*(Phi(r_geo)-Phi(r_r))/DH[Tmax]/1e3:,.0f}")
    print("(first pass, 100 seats at 2,000 m/s, 137.5 t: 6,221 t to the platform at 1,000 K, 2,624 t below the cut)")

    # ---------------- item 2: as hot as it can go ----------------
    section("2. THE TRACK AS HOT AS IT CAN GO: heat into the carriage and the water with a reflective shield")
    print(f"black-body flux of the ribbon at its peak: 1,000 K {SB*1000**4/1e3:.0f} kW/m2, 1,500 K {SB*1500**4/1e3:.0f}, 2,000 K {SB*2000**4/1e3:.0f} kW/m2; the carriage face sees it across the gap through")
    print(f"a reflective shield of effective emissivity {EPS_SHIELD[0]} or {EPS_SHIELD[1]} (Paul), the magnet pole faces plated to the same figure (ASSUMPTION); the ribbon at {EPS_TRACK}.")
    print(f"the water starts at {T_W0:.0f} C, limit {T_W_LIMIT:.0f} C; {e1.Q_PERSON:.0f} W per passenger; condenser radiator at 100 C, {Q_RAD:.0f} W/m2 net (ASSUMPTION: sink 250 K, emissivity 0.9).")
    print("\n  m/s | seats | ribbon K | shield eps | heat into the carriage GJ | passengers GJ | water t | rise K (no boiling) | wall C | water for 35 C t | water boiled if latent t | condenser m2 over the trip | peak face flux kW/m2")
    heat = {}
    for v in (2000.0, 3000.0):
        for N in (100, 400):
            g = pods[(N, v)]; d = results[(N, v, 1000.0)]["d"]
            for Tmax in (1000.0, 1500.0, 2000.0):
                tr = tracks[(N, v, Tmax)]
                for eps in EPS_SHIELD:
                    Q, q_peak = heat_into_water(d, tr, g["A_face"], eps); Qp = e1.Q_PERSON * N * d["total"]
                    dT = (Q + Qp) / (g["m_sh"] * C_W); m_w35 = (Q + Qp) / (C_W * (T_W_LIMIT - T_W0)); boiled = max(Q + Qp - g["m_sh"] * C_W * (100 - T_W0), 0.0) / L_BOIL
                    A_rad = Q / (d["total"] * Q_RAD); heat[(N, v, Tmax, eps)] = (Q, dT)
                    print(f"  {v:5,.0f} | {N:5d} | {Tmax:8.0f} | {eps:10.2f} | {Q/1e9:25.2f} | {Qp/1e9:13.2f} | {g['m_sh']/1e3:7.1f} | {dT:19.1f} | {T_W0+dT:6.1f} | {m_w35/1e3:16.1f} | {boiled/1e3:24.1f} | {A_rad:26.0f} | {q_peak/1e3:.1f}")
    g = pods[(100, 2000.0)]; d = results[(100, 2000.0, 1000.0)]["d"]; tr = tracks[(100, 2000.0, 2000.0)]
    print("\nshield emissivity that holds the wall under 35 C at 2,000 K (100 seats, 2,000 m/s):")
    for eps in (0.01, 0.02, 0.03, 0.05, 0.1):
        Q, _ = heat_into_water(d, tr, g["A_face"], eps); dT = (Q + e1.Q_PERSON * 100 * d["total"]) / (g["m_sh"] * C_W)
        print(f"  eps {eps:.2f}: {Q/1e9:.2f} GJ, wall {T_W0+dT:.1f} C {'(under 35 C)' if T_W0+dT < T_W_LIMIT else '(over)'}")
    print("\ncooling of the ribbon after a pass and the sublimation while it is hot (graphite rate fit, UNVERIFIED):")
    for Tmax in (1000.0, 1500.0, 2000.0):
        for t_mm, where in ((0.3, "cruise ribbon"), (8.0, "8 mm bottom")):
            lost, th = sublimation_per_pass(t_mm, Tmax)
            yr = 24 * 365 * lost                                     # one pass an hour
            print(f"  {Tmax:.0f} K, {where} {t_mm} mm: above 1,200 K for {th:.1f} s per pass, sublimation {lost*1e6:.3g} mg/m2 per pass, {yr*1e3:.3g} g/m2 per year at one pod an hour "
                  f"= {yr/(RHO*t_mm*1e-3)*100:.3g}% of the ribbon per year; to 350 K in {e1.track_cooling(t_mm, Tmax):.0f} s")
    print("  the mechanical wear of a non-contact brake cannot be stated; the thermal cycling of a thin carbon sheet to 2,000 K every pass is the real life question [QUESTION FOR AUTHOR].")

    # ---------------- item 3: where the cruise energy goes ----------------
    section("3. WHERE THE CRUISE ENERGY GOES: passive ribbon, regenerating stator over the stop, over the whole line, and free fall then one stop")
    print(f"resistive stator: CNT winding at {SIG_CU/1e6:.0f} MS/m with {F_LOSS*100:.0f}% I^2R loss: that loss is heat in the winding at every pass, so the winding weighs {F_LOSS*100:.0f}% of the heat-sized passive ribbon")
    print(f"at the same rise; the ribbon's own conductor carries the power to the ring at {V_LINE[0]/1e3:.0f} kV or {V_LINE[1]/1e6:.0f} MV with {F_LOSS*100:.0f}% loss over the run, sized for the largest power passing each point.")
    print(f"superconducting stator: {TAPE_PER_M2:.1f} kg/m2 of tape (the LIM's 1,630 kg over 3 x 100 m and 0.844 m) and {CRYO_PER_M2:.0f} kg/m2 of cryostat (ASSUMPTION), the same conductor; its leak {LEAK_PER_M2:.1f} W/m2 at COP {COP}.")
    print("the power electronics that put each stator segment's current onto the line are not counted [UNVERIFIED]. The pods and the tether as in item 1.")
    print("\n  m/s | ribbon K | option | track + stator below the cut t | conductor t | winding or tape t | passive ribbon t | tether whole kt | per passenger per day t | note")
    g = pods[(100, 2000.0)]
    for v in (2000.0, 3000.0):
        g = pods[(100, v)]; d = e1.descent(v, g["m"])
        for Tmax in (1000.0, 2000.0):
            tr = tracks[(100, v, Tmax)]; base = results[(100, v, Tmax)]
            print(f"  {v:5,.0f} | {Tmax:8.0f} | passive ribbon | {base['m_track']/1e3:30,.0f} | {0:11.0f} | {0:17.0f} | {base['m_track']/1e3:16,.0f} | {base['total']/1e6:15,.3f} | {base['total']/1e3/(100*24):23,.1f} |")
            for label, r_top in (("stator over the stop", d["r_stop"]), ("stator over the line below the cut", base["r_cut"])):
                for V in V_LINE:
                    rm = regen_masses(d, tr, r_top, V, DH[Tmax])
                    mu = np.where(rg < base["r_cut"], rm["mu"], 0.0); _, _, E_tot, _ = e1.pod_envelope(d, rg); res = e1.size_tether(rg, mu, np.full(len(rg), E_tot))
                    print(f"  {v:5,.0f} | {Tmax:8.0f} | {label} at {V/1e3:,.0f} kV | {np.trapezoid(mu, rg)/1e3:30,.0f} | {np.trapezoid(np.where(rg<base['r_cut'], rm['mu_cond'], 0), rg)/1e3:11,.0f} | "
                          f"{np.trapezoid(np.where(rg<base['r_cut'], rm['mu_wind'], 0), rg)/1e3:17,.0f} | {np.trapezoid(np.where(rg<base['r_cut'], rm['mu_pass'], 0), rg)/1e3:16,.0f} | {res['total']/1e6:15,.3f} | {res['total']/1e3/(100*24):23,.1f} | conductor up to {rm['A_cond_max']*1e4:,.0f} cm2, {rm['P_peak']/1e9:.2f} GW")
            rm = regen_masses(d, tr, base["r_cut"], V_LINE[1], DH[Tmax], superconducting=True)
            mu = np.where(rg < base["r_cut"], rm["mu"], 0.0); res = e1.size_tether(rg, mu, np.full(len(rg), E_tot))
            print(f"  {v:5,.0f} | {Tmax:8.0f} | superconducting stator below the cut at 1 MV | {np.trapezoid(mu, rg)/1e3:30,.0f} | {np.trapezoid(np.where(rg<base['r_cut'], rm['mu_cond'], 0), rg)/1e3:11,.0f} | "
                  f"{np.trapezoid(np.where(rg<base['r_cut'], rm['mu_wind'], 0), rg)/1e3:17,.0f} | 0 | {res['total']/1e6:15,.3f} | {res['total']/1e3/(100*24):23,.1f} | leak {rm['leak_W']/1e6:.1f} MW cold, {rm['leak_W']/COP/1e9:.1f} GW of compressors along the line")
    # free fall then one stop
    print("\nno cruise braking: the pod falls from just below the platform and takes the whole energy in a 1.5 g stop:")
    for v_c_label in ("100 seats at 2,000 m/s pod", "100 seats at 3,000 m/s pod"):
        g = pods[(100, 2000.0 if "2,000" in v_c_label else 3000.0)]
        dff = freefall_profile(g["m"]); d1 = freefall_profile(g["m"], v0=1.0)
        print(f"  (ASSUMPTION: a push to 10 m/s at the platform; with 1 m/s the fall takes {d1['t_fall']/3600:.1f} h, which is the chapter's 24 h)")
        p, e, tot = dose_profile(rates, dff)
        x5 = e1.shield_for(T_sh, tot, 5.0); x2 = e1.shield_for(T_sh, tot, 2.0)
        print(f"  {v_c_label} ({g['m']/1e3:.1f} t): free fall {dff['t_fall']/3600:.1f} h, enters the stop at {dff['v_entry']:,.0f} m/s {dff['s_stop']/1e3:,.0f} km above the ring, stop {(dff['total']-dff['t_fall'])/60:.1f} min, whole {dff['total']/3600:.1f} h; "
              f"below 6,000 km {dff['t_below_6000']/60:.0f} min; peak brake {dff['P_stop_peak']/1e9:.1f} GW; pods on the line at one an hour {dff['n_on_line']}")
        print(f"    dose at 10 / 20 / 30 / 50 g/cm2: {tot[6]:.1f} / {tot[8]:.1f} / {tot[9]:.2f} / {tot[11]:.2f} mSv (electrons {e[8]:.2f} mSv at 20 g/cm2 against {e1.dose_per_descent(rates, e1.descent(2000.0, 1.0))[1][8]:.2f} at 2,000 m/s); shield for 5 mSv {x5:.1f} g/cm2, for 2 mSv {x2:.1f} g/cm2 (2,000 m/s: {shield[(2000.0,5.0)]:.1f} and {shield[(2000.0,2.0)]:.1f})")
        for Tmax in (1000.0, 2000.0):
            tr = track2(g["m"], None, DH[Tmax], mag["p"], rg=rg, free_fall=True, v_entry=dff["v_entry"])
            res = size_with_track(dff, tr, rg)
            print(f"    ribbon to {Tmax:.0f} K: track {integrate(tr)/1e3:,.0f} t all in the stop ({res['m_track']/1e3:,.0f} below the cut at {(res['r_cut']-r_r)/1e3:,.0f} km), sigma asked {tr['sig_req'][tr['in_stop']].min()/1e6:.3f} to {min(tr['sig_req'][tr['in_stop']].max(),1e9)/1e6:.1f} MS/m, "
                  f"largest pod load {res['E_tot']/(g['m']*anet(r_r)):.2f} pod weights, tether whole {res['total']/1e6:.3f} Mt, {res['total']/1e3/(100*24):,.1f} t per passenger per day (with the 2 mSv shield of the cruise pod; the free-fall shield is {x2:.0f} g/cm2)")

    # ---------------- item 4: speed and the stop ----------------
    section("4. SPEED AND THE STOP: 4,000 m/s, and the stop acceleration")
    print("the speeds (2 mSv shield, ribbon to 2,000 K; 1,000 K in brackets):")
    print("  m/s | seats | pod t | belt time min | dose at the design shield mSv | shield g/cm2 | peak brake GW | track below cut t | largest pod load, pod weights | tether whole kt | per passenger per day t")
    for v in SPEEDS:
        for N in (100, 400):
            g = pods[(N, v)]; r2 = results[(N, v, 2000.0)]; r1 = results[(N, v, 1000.0)]; d = r2["d"]
            print(f"  {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {d['t_below_6000']/60:13.0f} | 2.0 | {shield[(v,2.0)]:12.1f} | {d['P_stop_peak']/1e9:13.2f} | {r2['m_track']/1e3:,.0f} ({r1['m_track']/1e3:,.0f}) | {r2['E_tot']/(g['m']*anet(r_r)):29.2f} | "
                  f"{r2['total']/1e3:,.0f} ({r1['total']/1e3:,.0f}) | {r2['total']/1e3/(N*24):,.1f} ({r1['total']/1e3/(N*24):,.1f})")
    print("\nthe stop (100 seats, 2 mSv): tether and track against the stop acceleration, track's stop and cruise sections separately")
    print("  m/s | stop g | stop km | belt min | shield g/cm2 | pod t | brake force at ring, pod weights | largest pod load, pod weights | ribbon K | track stop section t | cruise section t | below cut t | tether whole kt | per passenger per day t")
    stopg = {}
    for v in (2000.0, 3000.0):
        for ag in (1 / 3, 1.0, 1.5, 2.0):
            a_stop = ag * g0; d1 = e1.descent(v, 1000.0, a_stop=a_stop); p, e, tot = e1.dose_per_descent(rates, d1); x = e1.shield_for(T_sh, tot, 2.0)
            g = size_pod2(100, x, mag, a_stop=a_stop); d = e1.descent(v, g["m"], a_stop=a_stop)
            for Tmax in (1000.0, 2000.0):
                tr = track2(g["m"], v, DH[Tmax], mag["p"], a_stop=a_stop, rg=rg); res = size_with_track(d, tr, rg)
                stopg[(v, ag, Tmax)] = dict(res=res, m_stop=integrate(tr, tr["in_stop"]), m_cr=integrate(tr, ~tr["in_stop"]), m_pod=g["m"])
                print(f"  {v:5,.0f} | {ag:6.2f} | {d['s_stop']/1e3:7,.0f} | {d['t_below_6000']/60:8.0f} | {x:12.1f} | {g['m']/1e3:5.1f} | {(a_stop+anet(r_r))/anet(r_r):31.2f} | {res['E_tot']/(g['m']*anet(r_r)):29.2f} | {Tmax:8.0f} | "
                      f"{integrate(tr, tr['in_stop'])/1e3:20,.0f} | {integrate(tr, ~tr['in_stop'])/1e3:16,.0f} | {res['m_track']/1e3:11,.0f} | {res['total']/1e3:15,.0f} | {res['total']/1e3/2400:.1f}")

    # ---------------- item 6: throughput ----------------
    section("6. THROUGHPUT: pod spacing of 15, 30 and 60 minutes")
    print("pods on the line = trip / spacing; the largest pod load over the phase of the train; passengers per day = seats x 1440 / spacing; the ribbon's temperature when the next pod")
    print("arrives, at the 8 mm bottom and in the cruise ribbon, from the cooling curve; and the peak the next pass then reaches at the bottom (its rise is fixed by its mass).")
    print("\n  spacing min | m/s | seats | pod t | ribbon K | pods on line | largest pod load, pod weights | tether whole kt | passengers per day | tether per passenger per day t | bottom K at next pod | next peak at bottom K | cruise ribbon K at next pod")
    thr = []
    for sp in SPACINGS:
        e1.HEADWAY = sp
        for v in SPEEDS:
            for N in (100, 400):
                g = pods[(N, v)]
                for Tmax in (1000.0, 2000.0):
                    d = e1.descent(v, g["m"]); tr = tracks[(N, v, Tmax)]; res = size_with_track(d, tr, rg)
                    t_bot = float(np.interp(r_r + V_HAND**2 / (2 * e1.A_STOP), rg, tr["t"])) * 1e3; t_cr = float(np.interp(d["r_stop"] + 1e3, rg, tr["t"])) * 1e3
                    T_bot_peak = e1.T_from_dH(float(np.interp(r_r + V_HAND**2 / (2 * e1.A_STOP), rg, tr["dH_local"])))
                    Tb = cooling_curve(t_bot, T_bot_peak, sp); Tc = cooling_curve(t_cr, Tmax, sp)
                    T_next = e1.T_from_dH(e1.h_graphite(Tb) - e1.h_graphite(T0) + float(np.interp(r_r + V_HAND**2 / (2 * e1.A_STOP), rg, tr["dH_local"])))
                    pax = N * 86400 / sp; row = dict(sp=sp, v=v, N=N, Tmax=Tmax, pax=pax, tpp=res["total"] / 1e3 / pax, total=res["total"], n=d["n_on_line"], E=res["E_tot"] / (g["m"] * anet(r_r)), Tb=Tb, Tn=T_next, Tc=Tc, d=d)
                    thr.append(row)
                    print(f"  {sp/60:11.0f} | {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {Tmax:8.0f} | {d['n_on_line']:12d} | {row['E']:29.2f} | {res['total']/1e3:15,.0f} | {pax:18,.0f} | {row['tpp']:30,.1f} | {Tb:20.0f} | {T_next:21.0f} | {Tc:.0f}")
    e1.HEADWAY = 3600.0
    # figures: tether per passenger per day against passengers per day; dose, belt time and peak brake power against speed
    fig, ax = plt.subplots(figsize=(8, 5))
    mk = {100: "o", 400: "s"}; col = {2000.0: BLUE, 3000.0: GREEN, 4000.0: RED}
    for Tmax, alpha in ((2000.0, 1.0), (1000.0, 0.35)):
        for v in SPEEDS:
            for N in (100, 400):
                pts = sorted([(r["pax"], r["tpp"]) for r in thr if r["v"] == v and r["N"] == N and r["Tmax"] == Tmax])
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col[v], marker=mk[N], ms=5, lw=1.2, alpha=alpha, label=f"{v:,.0f} m/s, {N} seats, {Tmax:.0f} K" if Tmax == 2000.0 else None)
    ax.set_xscale("log"); ax.set_xlabel("passengers per day"); ax.set_ylabel("tether per passenger per day (t)"); ax.set_ylim(0, None); ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, title="solid: 2,000 K ribbon; faint: 1,000 K; points: 60, 30, 15, 5 min spacing", title_fontsize=7)
    ax.set_title("whole tether per passenger per day, 2 mSv shield, 1.5 g stop", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "throughput.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
    for N, ls in ((100, "-"), (400, "--")):
        axs[0].plot(SPEEDS, [dose[v][np.searchsorted(T_sh, 50)] for v in SPEEDS], color=BLUE, ls="-", marker="o", ms=4, label="at 50 g/cm2" if N == 100 else None)
        axs[0].plot(SPEEDS, [dose[v][np.searchsorted(T_sh, 30)] for v in SPEEDS], color=RED, ls="-", marker="o", ms=4, label="at 30 g/cm2" if N == 100 else None)
        axs[1].plot(SPEEDS, [results[(N, v, 2000.0)]["d"]["t_below_6000"] / 60 for v in SPEEDS], color=GREEN, ls=ls, marker="o", ms=4, label=f"{N} seats")
        axs[2].plot(SPEEDS, [results[(N, v, 2000.0)]["d"]["P_stop_peak"] / 1e9 for v in SPEEDS], color=RED, ls=ls, marker="o", ms=4, label=f"{N} seats")
    axs[0].set_ylabel("dose per descent (mSv)"); axs[1].set_ylabel("time below 6,000 km (min)"); axs[2].set_ylabel("peak brake power (GW)")
    for a in axs: a.set_xlabel("cruise speed (m/s)"); a.set_xticks(SPEEDS); a.grid(alpha=0.3); a.set_ylim(0, None); a.legend(fontsize=7)
    fig.suptitle("dose, belt time and peak brake power against cruise speed (1.5 g stop, 2 mSv pods)", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "speed_tradeoffs.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    fig, ax = plt.subplots(figsize=(8, 4.8)); ax2 = ax.twinx()
    for v, c in ((2000.0, BLUE), (3000.0, RED)):
        gs = [1 / 3, 1.0, 1.5, 2.0]
        ax.plot(gs, [stopg[(v, ag, 2000.0)]["res"]["total"] / 1e3 for ag in gs], color=c, marker="o", ms=4, label=f"tether, {v:,.0f} m/s")
        ax2.plot(gs, [stopg[(v, ag, 2000.0)]["res"]["m_track"] / 1e3 for ag in gs], color=c, ls="--", marker="s", ms=4, label=f"track below the cut, {v:,.0f} m/s")
        ax2.plot(gs, [stopg[(v, ag, 2000.0)]["m_stop"] / 1e3 for ag in gs], color=c, ls=":", marker="^", ms=4, label=f"track, stop section, {v:,.0f} m/s")
    ax.set_xlabel("stop deceleration (g)"); ax.set_ylabel("whole tether (kt)"); ax2.set_ylabel("track (t)"); ax.set_ylim(0, None); ax2.set_ylim(0, None); ax.grid(alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=7)
    ax.set_title("100-seat pod, 2 mSv shield, ribbon to 2,000 K: the stop's acceleration", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "stop_g.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    with open(os.path.join(FIG, "throughput.csv"), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["spacing_min", "speed_m_s", "seats", "ribbon_K", "pods_on_line", "largest_load_pod_weights", "tether_whole_t", "passengers_per_day", "tether_per_passenger_per_day_t", "bottom_K_at_next_pod", "next_peak_bottom_K", "cruise_ribbon_K_at_next_pod"])
        for r in thr: wr.writerow([f"{r['sp']/60:.0f}", f"{r['v']:.0f}", r["N"], f"{r['Tmax']:.0f}", r["n"], f"{r['E']:.2f}", f"{r['total']/1e3:.0f}", f"{r['pax']:.0f}", f"{r['tpp']:.1f}", f"{r['Tb']:.0f}", f"{r['Tn']:.0f}", f"{r['Tc']:.0f}"])

    # ---------------- item 7: safety ----------------
    section("7. THE SAFETY CASE")
    H = scale_height()
    print(f"Allen and Eggers (NACA Report 1381, 1958): for a steep ballistic entry in an exponential atmosphere the peak deceleration is v_E^2 sin(gamma) / (2 e H), with H the scale height")
    print(f"and gamma the entry angle; it does not contain the ballistic coefficient, which only moves the altitude of the peak, H ln(rho0 H / (beta sin gamma)). Scale height of the model atmosphere between 30 and 60 km: {H/1e3:.2f} km.")
    g = pods[(100, 2000.0)]
    print(f"\n100 seats ({g['m']/1e3:.1f} t, radius {g['R_out']:.2f} m, beta {g['beta']:,.0f} kg/m2) released at rest on the tether; the same pod with a ballute (beta 300 kg/m2) and a heavy pod (beta 10,000):")
    print("  release km | entry speed m/s | entry angle deg | Allen-Eggers peak g | simulated peak g (beta of the pod) | ballute 300 kg/m2 | heavy 10,000 kg/m2 | seconds above 20 g (pod)")
    for alt0 in (700e3, 2000e3, 5000e3, 10000e3, 20000e3):
        r0 = RE + alt0; f = fall(r0, 0.0, W * r0, g["R_out"], g["beta"]); fb = fall(r0, 0.0, W * r0, g["R_out"], 300.0); fh = fall(r0, 0.0, W * r0, g["R_out"], 10000.0)
        ae = f["v_entry"]**2 * math.sin(math.radians(f["gamma"])) / (2 * math.e * H) / g0
        print(f"  {alt0/1e3:10,.0f} | {f['v_entry']:15,.0f} | {f['gamma']:15.1f} | {ae:19.0f} | {f['g_peak']:33.0f} | {fb['g_peak']:17.0f} | {fh['g_peak']:18.0f} | {f['t20']:.1f}")
    print("so neither a ballute nor a heavier pod lowers the peak. Lift would (a lifting body flies a shallower path and spreads the deceleration), or a shallower entry, which a pod falling")
    print("from a tether does not have: its only sideways speed is the tether's own, W r, under 3 km/s, against a fall of 3 to 9 km/s.")

    # break below the pod
    print("\nA break below the pod: the pod rides the upper piece. The piece from the break to the tether's end was held by the tension at the break; after the break that pull is unbalanced")
    print("and the piece accelerates outward at A = T_break / M_piece. The tension at a point of the piece is the sum over the part below it of (a(r) + A) dm, the inward pull of that part")
    print("plus what it takes to drag it outward at A; it is zero at the free end and positive wherever that sum is. A pod on the piece hangs on it while a(r_pod) + A > 0, which is")
    print("everywhere below GEO; above GEO the pod is flung outward with the piece.")
    print("The piece's centre of mass then moves as a body with the tether's rotation: its energy says whether it stays bound. (Rigid-line treatment at the instant of the break [UNVERIFIED beyond that].)")
    base = results[(100, 2000.0, 2000.0)]; A_prof = base["A"]; mu_pass = RHO * A_prof + base["mu_cut"]; d = base["d"]; m_pod = g["m"]
    T_geo = base["T_geo"]; L_up = base["L_up"]; m_up = base["m_up_per_m"]; r_end = base["r_up_end"]
    print(f"  design case 100 seats, 2,000 m/s, 2,000 K: lower tether {base['m_lower']/1e3:,.0f} t, upper {base['m_upper']/1e3:,.0f} t to {(r_end-RE)/1e3:,.0f} km, tension at the platform {T_geo/1e6:.0f} MN")
    print("  break km | pod km | piece mass kt | tension at the break MN | outward acceleration m/s2 | pull at the break a(r) m/s2 | piece slack anywhere | piece CM km, speed m/s (escape speed there) | bound? perigee / apogee km | the pod")
    rg_up = np.linspace(r_geo, r_end, 4000)
    for br_alt in (2000e3, 10000e3, 25000e3):
        rb = RE + br_alt
        for pod_alt in (br_alt + 100e3, 0.5 * (br_alt + (r_geo - RE)), (r_geo - RE) - 500e3):
            rp = RE + pod_alt
            # mass and first moment of the piece: lower part rb..r_geo, upper r_geo..r_end, the pod at rp
            sel = rg >= rb; m_low = float(np.trapezoid(mu_pass[sel], rg[sel])); M = m_low + m_up * L_up + m_pod
            rcm = (float(np.trapezoid(mu_pass[sel] * rg[sel], rg[sel])) + m_up * 0.5 * (r_end**2 - r_geo**2) + m_pod * rp) / M
            # actual tension at the break with this one pod on the line: T_geo minus the passive weight above the break minus the pod's weight if above
            W_above = float(np.trapezoid((mu_pass * anet(rg))[sel], rg[sel])); T_b = T_geo - W_above - (m_pod * anet(rp) if rp > rb else 0.0)
            A = T_b / M
            # tension along the piece after the break: integral of (a - A) dm from the break upward, pod included
            Tp = []; acc = 0.0; prev_r = rb; slack = False
            for k in np.where(sel)[0][::10]:
                r_k = rg[k]; acc += (anet(r_k) + A) * mu_pass[k] * (r_k - prev_r); prev_r = r_k
                if rp <= r_k and not any(t[0] >= rp for t in Tp): acc += (anet(rp) + A) * m_pod
                Tp.append((r_k, acc))
            T_min = min(t[1] for t in Tp)
            for r_k in rg_up[::20]:
                acc += (anet(r_k) + A) * m_up * (rg_up[1] - rg_up[0]) * 20; T_min = min(T_min, acc)
            v_cm = W * rcm; eps = 0.5 * v_cm**2 - GM / rcm
            if eps < 0:
                a_ = -GM / (2 * eps); e_ = math.sqrt(max(1 - (rcm * v_cm)**2 / (GM * a_), 0.0)); orb = f"bound, {(a_*(1-e_)-RE)/1e3:,.0f} / {(a_*(1+e_)-RE)/1e3:,.0f}"
            else: orb = "unbound, escapes"
            pod_note = ("hangs on the piece and is carried outward" if T_min >= 0 else "the piece goes slack below the pod") if anet(rp) + A > 0 else "above GEO, flung outward with the piece"
            print(f"  {br_alt/1e3:8,.0f} | {pod_alt/1e3:6,.0f} | {M/1e6:12.2f} | {T_b/1e6:22.0f} | {A:24.3f} | {anet(rb):27.3f} | {'yes' if T_min < 0 else 'no':19s} | {(rcm-RE)/1e3:,.0f}, {v_cm:,.0f} ({math.sqrt(2*GM/rcm):,.0f}) | {orb:28s} | {pod_note}")
    print("  a break above the pod: the pod and the ribbon below it fall together (the ribbon can no longer hold the pod up), the pod free-falls as in the table above, and the piece below the")
    print("  break comes down on the ring; a break at the ring itself frees the whole elevator, whose centre of mass is above GEO and whose speed W r exceeds the circular speed there: it escapes.")

    # the survivable set
    print("\nThe survivable set. By heat (first pass, 2 to 3% of the kinetic energy into the vehicle) a release below about 2,000 km boils a quarter of the water, below 5,000 km half.")
    print("By deceleration: Eiband (NASA Memorandum 5-19-59E, 1959) reports spineward (eyeballs-in) plateaus endured without injury of 45 g for 0.044 s, 40 g for 0.1 s and 25 g for 0.2 s by")
    print("well-restrained subjects; the plateaus here are 7 to 10 s above 20 g from every release altitude, far beyond the longest durations in his figures [the 1 to 10 s boundary UNVERIFIED here].")
    print("So X, the altitude below which a break above the pod is survivable without thrust, is zero: no ballistic release from the tether is survivable by a passenger, and the band of loss is the")
    print(f"whole descent below 23,512 km. Share of the trip spent below 23,512 km (and above the ring):")
    for v in SPEEDS:
        d = prof[v]; frac = float(np.sum(d["alt"] < 23512e3)) * d["dt"] / d["total"]
        print(f"  {v:,.0f} m/s: {frac*100:.0f}% of {d['total']/3600:.2f} h = {frac*d['total']/3600:.2f} h")

    # thrust
    print("\nThrust: the smallest single impulse at release that lifts the perigee above 100 km (entry avoided) or above 300 km (a rescue orbit), and the propellant fraction at Isp 230 and 300 s.")
    print("  release km | moving down m/s | dv for perigee 100 km m/s (propellant at 230 / 300 s) | dv for 300 km (propellant) | with 5% of the pod at 300 s: entry speed, peak g | with 10%: entry speed, peak g")
    thrust = {}
    for alt0 in (700e3, 2000e3, 5000e3, 10000e3, 15000e3, 20000e3, 21000e3, 22000e3, 22500e3, 23000e3):
        r0 = RE + alt0
        for vd in (0.0, 2000.0, 3000.0):
            vr0, vt0 = -vd, W * r0
            dv1, u1 = min_dv(r0, vr0, vt0, RE + 100e3); dv3, u3 = min_dv(r0, vr0, vt0, RE + 300e3)
            def pf(dv, isp): return 1 - math.exp(-dv / (isp * g0))
            outs = []
            for fb in BUDGETS:
                dvb = 300.0 * g0 * (-math.log(1 - fb))
                if dvb >= dv3: outs.append("rescue orbit reached")
                else:
                    ux, uy = u3 if u3 else (1.0, 0.0)
                    f = fall(r0, vr0 + dvb * ux, vt0 + dvb * uy, g["R_out"], g["beta"])
                    outs.append("orbit" if f["orbit"] else f"{f['v_entry']:,.0f} m/s, {f['g_peak']:.0f} g")
            thrust[(alt0, vd)] = (dv1, dv3)
            print(f"  {alt0/1e3:10,.0f} | {vd:14,.0f} | {dv1:8,.0f} ({pf(dv1,230)*100:.0f}% / {pf(dv1,300)*100:.0f}%) | {dv3:8,.0f} ({pf(dv3,230)*100:.0f}% / {pf(dv3,300)*100:.0f}%) | {outs[0]:28s} | {outs[1]}")
    for isp in ISP:
        for fb in BUDGETS:
            dvb = isp * g0 * (-math.log(1 - fb))
            for vd in (0.0, 2000.0, 3000.0):
                ok = [alt for (alt, v_) in sorted(thrust) if v_ == vd and thrust[(alt, v_)][1] <= dvb]
                print(f"  Isp {isp:.0f} s, {fb*100:.0f}% of the pod = {dvb:.0f} m/s: a rescue orbit from {ok[0]/1e3:,.0f} km up at {vd:,.0f} m/s down" if ok else f"  Isp {isp:.0f} s, {fb*100:.0f}% of the pod = {dvb:.0f} m/s: no rescue orbit from any altitude run at {vd:,.0f} m/s down")
    f0 = fall(RE + 700e3, 0.0, W * (RE + 700e3), g["R_out"], g["beta"])
    print(f"  for scale: from 700 km at rest the pod enters at {f0['v_entry']:,.0f} m/s and sees {f0['g_peak']:.0f} g; from the top of the stop while moving at 2,000 m/s the first pass gave 51 g.")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    alts = sorted(set(a for a, _ in thrust))
    for vd, c in ((0.0, BLUE), (2000.0, GREEN), (3000.0, RED)):
        ax.plot([a / 1e3 for a in alts], [thrust[(a, vd)][1] for a in alts], color=c, marker="o", ms=4, label=f"rescue orbit (perigee 300 km), moving down at {vd:,.0f} m/s")
        ax.plot([a / 1e3 for a in alts], [thrust[(a, vd)][0] for a in alts], color=c, ls="--", marker="s", ms=3, label=f"perigee 100 km, {vd:,.0f} m/s")
    for isp, fb, c in ((300.0, 0.05, GREY), (300.0, 0.10, GREY), (230.0, 0.10, "#bbbbbb")):
        dvb = isp * g0 * (-math.log(1 - fb)); ax.axhline(dvb, color=c, ls=":", lw=1); ax.text(alts[0] / 1e3, dvb + 20, f"{fb*100:.0f}% of the pod at Isp {isp:.0f} s: {dvb:.0f} m/s", color=c, fontsize=7)
    ax.set_xlabel("release altitude (km)"); ax.set_ylabel("impulse needed (m/s)"); ax.set_ylim(0, None); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("impulse at release that keeps the pod out of the air", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "rescue_dv.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    print("done")

if __name__ == "__main__":
    main()
