#!/usr/bin/env python3
"""
The elevator pod, its shield, the graded brake track, the carriage magnets and the tether sized as one system.
Book: Orbital Ring Engineering, Volume I: Chapter 13 (pods of 25 to 400 seats at 2,000 and 3,000 m/s with a 1.5 g stop).
Run: python elevator_system.py > elevator_system.txt      (reads runs/elevator_dose_rates.npz from radiation/elevator_dose_rates.py)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py);
             ORBITAL_RING_OUT (folder of the dose-rate file and of the figures and CSVs, default runs)

What it does, in the order of the report (Vol_III/runs_700km_hw10686_2026-10-12/elevator_system_report_2026-10-28.md):
  1. descent profiles: 1 g drive from the platform, cruise at v_c, 1.5 g stop at the ring (felt 2.3 g), wheels below 100 m/s
  2. dose per descent behind water shields from the dose rates along the tether, and the shield for 5 and 2 mSv per descent
  3. permanent-magnet Halbach arrays on the carriage: field at a 20 mm gap, peak eddy drag B^2/(4 mu0) per m^2 of face
     (thin-sheet model of Reitz 1970 and Rossing and Hull 1991), magnet mass per newton of brake force
  4. the pod: a sphere with a tube through it; seats, water-polyethylene shield, aluminum hull, magnets, carriage
  5. the brake track: a CNT ribbon whose conductance sigma*t tracks the pod speed (characteristic speed 2/(mu0 sigma t) = v)
     and whose mass takes the heat of each pass (energy per metre = brake force) within a 300 K to 1,000 K or 1,500 K rise
  6. heat radiated from the track onto the magnets and into the water, and the water temperature
  7. the lower tether with the pods on the line and the track as distributed loads, the upper tether to balance, at 18.75 GPa
  8. loss of the tether: free fall of the pod from 700 to 35,000 km, entry speed, heat, water boiled, peak g
  9. cargo drop pods: release from the ring and from a trolley line at 300, 200 and 100 km; faster cargo descents
Every assumption is printed with the word ASSUMPTION so the report can list them.
"""
import os, sys, math, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring_altitude as ring
import tether_config as tcfg
import freefall_trajectory as ff

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
FIG = os.path.join(OUT, "elevator_system"); os.makedirs(FIG, exist_ok=True)
plt.rcParams["font.family"] = "DejaVu Sans"
BLUE, RED, GREEN, GREY = "#2b6ca3", "#b5451b", "#2e8b57", "#777777"

# ----------------------------------------------------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------------------------------------------------
GM = ring.GM; RE = ring.R_EARTH_EQ; W = ring.OMEGA_SIDEREAL
r_r = ring.R_ORBIT                      # the ring
r_geo = (GM / W**2)**(1.0 / 3.0)        # the platform (GEO), as ch17_options.py
g0 = 9.80665; MU0 = 4e-7 * math.pi; SB = 5.67e-8
RHO = 1700.0                            # CNT yarn, kg/m3 (Chapter 13)
SIG = 18.75e9                           # working stress, Pa (Chapter 13, Table 13.3)

# the pod (Paul, 2026-10-27/28)
V_SEAT = 0.56                           # m3 per seat, 0.5 x 0.8 x 1.4 m, reclined
M_SEAT = 130.0                          # kg per seat with seat and supplies
HULL_T_AL = 0.005                       # ASSUMPTION: aluminum hull 5 mm
RHO_AL = 2700.0
X_AL = HULL_T_AL * 100 * 2.7            # g/cm2 of aluminum in the hull (1.35)
AL_AS_WATER = 0.8                       # aluminum counts at 0.8 of water per g/cm2 (Paul)
RHO_SHIELD = 1000.0                     # ASSUMPTION: the water-polyethylene layer counted as water, 1,000 kg/m3
SKIN = 1.0                              # ASSUMPTION: reflective skin with self-healing layer, 1 kg/m2
M_FIXED = 2000.0                        # ASSUMPTION: coolant loop, pumps, wheels, avionics, 2 t per pod
CARRIAGE_FACTOR = 1.0                   # ASSUMPTION: carriage structure (frames, hinges, back plates) equal to the magnet mass
TUBE_CLEAR = 0.15                       # ASSUMPTION: 150 mm of structure between the magnet backs and the tube wall
C_D = tcfg.C_D_POD                      # 1.2, blunt body (tether_config)
SEATS = [25, 50, 100, 200, 400]
HARDWARE_START = (8000.0, 0.05)         # the starting allowance of the task: 8 t plus 5% of the rest

# the descent
A_DRIVE = g0; A_STOP = 1.5 * g0
V_HAND = 100.0                          # ASSUMPTION: below 100 m/s the wheels (and the lowest coils) finish the stop
HEADWAY = 3600.0                        # one pod per hour
PODS_PER_DAY = 24
SPEEDS = [2000.0, 3000.0]
DT = 0.5

# the track
SIG_CNT = 2e6                           # S/m, plain CNT ribbon (Paul)
T_TRACK0 = 300.0                        # ASSUMPTION: the track is back at 300 K before each pod (it cools in seconds, see the log)
# NIST-JANAF graphite, C1(ref): T (K), H - H(298.15) (kJ/mol); Cp (J/mol/K) 8.517 at 298, 14.623 at 500, 19.827 at 800,
# 21.610 at 1000, 22.766 at 1200, 23.904 at 1500 K. janaf.nist.gov/tables/C-002.html
JANAF_T = np.array([298.15, 300, 500, 800, 1000, 1200, 1500])
JANAF_H = np.array([0.0, 0.016, 2.365, 7.637, 11.795, 16.240, 23.253]) * 1e3 / 12.011     # J/g
def h_graphite(T): return np.interp(T, JANAF_T, JANAF_H) * 1e3                           # J/kg above 298.15 K
DH = {1000.0: h_graphite(1000.0) - h_graphite(T_TRACK0), 1500.0: h_graphite(1500.0) - h_graphite(T_TRACK0)}   # J/kg per pass
EPS_TRACK = tcfg.EPSILON_CNT_IR         # 0.85 (tether_config); unverified at 1,000 to 1,500 K
EPS_MAG = 0.10                          # ASSUMPTION: polished or plated magnet face, emissivity 0.1
TETHER_SHARE = 0.1                      # ASSUMPTION: the track ends where the tether's own mass per metre is 10 x the track's,
                                        #   so the bare tether warms by a tenth of the track's rise per pass

# the magnets
GAP = 0.020
M_SEG = 4                               # segments per wavelength of the Halbach array
MAGNETS = {"NdFeB": dict(Br=1.35, rho=7500.0, T_max="80 to 150 C by grade"),
           "SmCo": dict(Br=1.10, rho=8400.0, T_max="about 300 C")}
PITCHES = [0.050, 0.075, 0.100, 0.125, 0.150, 0.200]
THICKNESSES = [0.025, 0.050, 0.075, 0.100, 0.150]
L_CARRIAGES = [5.0, 10.0, 20.0]
L_CAR = 10.0                            # the shortest carriage used for the pod and track sizing
W_MAX = 1.0                             # ASSUMPTION: the ribbon is at most 1 m wide; a heavier pod gets a longer carriage, A_face/(2 W_MAX)
K_MAG = 9.0                             # ASSUMPTION: thermal conductivity of sintered NdFeB, 9 W/m/K

# the water
C_W = 4186.0; T_W0 = 20.0; T_W_LIMIT = 35.0   # ASSUMPTION: water from 20 C; wall limit 35 C (see the report)
Q_PERSON = 100.0                        # ASSUMPTION: 100 W of metabolic heat per passenger

# free fall and cargo
F_HEAT = (0.02, 0.03)                   # share of the kinetic energy reaching the vehicle behind an ablative shield (Paul; unverified)
H_BOIL = C_W * 80.0 + 2.256e6           # J/kg to boil water from 20 C
RELEASE_ALTS = [700e3, 2000e3, 5000e3, 10000e3, 20000e3, 35000e3]
CARGO_RHO = 500.0                       # ASSUMPTION: cargo packs at 500 kg/m3
CARGO_STRUCT = 0.10                     # ASSUMPTION: cargo pod structure 10% of the cargo
PARACHUTE_KG_PER_T = 25.0               # derived below from an 8 m/s splashdown; ASSUMPTION on the canopy areal mass
TROLLEY_LOAD = 100e3                    # a trolley line carries 100 t at its lower end (Chapter 11)

def anet(r): return GM / r**2 - W * W * r
def Phi(r): return -GM / r - 0.5 * W * W * r * r
def pull_above(r2): return 0.5 * W * W * (r2 * r2 - r_geo * r_geo) + GM * (1 / r2 - 1 / r_geo)   # J/kg, constant-section line from GEO

def fmt(x, nd=0): return f"{x:,.{nd}f}"
def section(title): print("\n" + "=" * 110 + "\n" + title + "\n" + "=" * 110)

# ----------------------------------------------------------------------------------------------------------------------
# 1. descent profiles
# ----------------------------------------------------------------------------------------------------------------------
def descent(v_c, m, a_drive=A_DRIVE, a_stop=A_STOP, dt=DT):
    """time series of a descent from the platform to the ring: drive at a_drive (gravity helping), cruise at v_c, stop at a_stop.
    Arrays: t, r, v, F (force the brakes put on the track, N: 0 in the drive, m*anet in the cruise, m*(anet + a_stop) in the stop),
    P (brake power, W), phase (0 drive, 1 cruise, 2 stop). The stop is integrated upward from the ring and reversed."""
    T_, R_, V_, F_, P_, PH = [], [], [], [], [], []
    r = r_geo; v = 0.0; t = 0.0; P_drive = 0.0
    while v < v_c:
        v += (a_drive + anet(r)) * dt; r -= v * dt; t += dt; P_drive = max(P_drive, m * a_drive * v)
        T_.append(t); R_.append(r); V_.append(v); F_.append(0.0); P_.append(0.0); PH.append(0)
    s_drive = r_geo - r; t_drive = t
    r3 = r_r; v3 = 0.0; up = []
    while v3 < v_c:
        v3 += a_stop * dt; r3 += v3 * dt; up.append((r3, min(v3, v_c)))
    r_stop = r3
    while r > r_stop:
        r -= v_c * dt; t += dt; F = m * anet(r)
        T_.append(t); R_.append(r); V_.append(v_c); F_.append(F); P_.append(F * v_c); PH.append(1)
    t_cruise = t - t_drive
    for (rr, vv) in reversed(up):
        t += dt; F = m * (anet(rr) + a_stop)
        T_.append(t); R_.append(rr); V_.append(vv); F_.append(F); P_.append(F * vv); PH.append(2)
    t_stop = t - t_drive - t_cruise
    d = {k: np.array(v) for k, v in dict(t=T_, r=R_, v=V_, F=F_, P=P_, phase=PH).items()}
    d["alt"] = d["r"] - RE; d["s_km"] = (r_geo - d["r"]) / 1e3; d["dt"] = dt; d["m"] = m; d["v_c"] = v_c
    d["r_stop"] = r_stop; d["s_stop"] = r_stop - r_r; d["t_drive"] = t_drive; d["s_drive"] = s_drive; d["P_drive"] = P_drive
    d["E_drive"] = 0.5 * m * v_c**2; d["t_cruise"] = t_cruise; d["t_stop"] = t_stop; d["total"] = t
    d["P_cruise_end"] = m * anet(r_stop) * v_c; d["P_stop_peak"] = d["P"][d["phase"] == 2].max()
    d["F_max"] = m * (anet(r_r) + a_stop); d["F_stop_start"] = m * (anet(r_stop) + a_stop)
    d["t_below_6000"] = dt * np.sum(d["alt"] < 6000e3)
    d["n_on_line"] = int(math.ceil(t / HEADWAY))
    return d

def print_descent(d, label):
    print(f"\n{label}: pod {d['m']/1e3:,.1f} t, cruise {d['v_c']:,.0f} m/s ({d['v_c']*3.6:,.0f} km/h), drive {A_DRIVE/g0:.1f} g, stop {A_STOP/g0:.1f} g "
          f"(felt {(A_STOP+anet(r_r))/g0:.2f} g at the ring)")
    print(f"  drive: {d['t_drive']/60:.1f} min over {d['s_drive']/1e3:,.0f} km, peak drive power {d['P_drive']/1e6:,.0f} MW, energy {d['E_drive']/1e9:,.0f} GJ ({d['E_drive']/3.6e9:,.1f} MWh)")
    print(f"  cruise: {d['t_cruise']/3600:.2f} h, brake force at the end of the cruise {d['P_cruise_end']/d['v_c']/1e3:,.0f} kN, brake power {d['P_cruise_end']/1e6:,.0f} MW")
    print(f"  stop: {d['t_stop']/60:.1f} min over {d['s_stop']/1e3:,.0f} km, peak brake power {d['P_stop_peak']/1e6:,.0f} MW at the start of the stop, "
          f"brake force {d['F_stop_start']/1e3:,.0f} kN there and {d['F_max']/1e3:,.0f} kN at the ring ({d['F_max']/1e3:,.0f} kJ per metre of track)")
    print(f"  whole trip {d['total']/3600:.2f} h ; {d['t_below_6000']/60:.0f} min below 6,000 km ; pods on the line at one per hour: {d['n_on_line']}")

# ----------------------------------------------------------------------------------------------------------------------
# 2. dose and shield
# ----------------------------------------------------------------------------------------------------------------------
def load_rates():
    f = os.path.join(OUT, "elevator_dose_rates.npz")
    if not os.path.exists(f):
        sys.exit(f"missing {f}: run radiation/elevator_dose_rates.py first (needs aep8)")
    return np.load(f)

def dose_per_descent(rates, d, key="max_+0"):
    """mSv per descent (protons, electrons, total) at each shield depth rates['T'] for the descent d"""
    sec, _ = np.histogram(d["alt"], bins=rates["edges"]); sec = sec.astype(float) * d["dt"]
    p = (rates[f"Hp_{key}"] * sec[None, :]).sum(1) * 1e3; e = (rates[f"De_{key}"] * sec[None, :]).sum(1) * 1e3
    return p, e, p + e

def shield_for(T, tot, target):
    """areal density (g/cm2 of water) for a dose target, by log interpolation on the decreasing total"""
    i = np.searchsorted(-tot, -target)
    if i == 0 or i >= len(T): return float("nan")
    return float(np.exp(np.interp(np.log(target), np.log(tot[[i, i - 1]]), np.log(T[[i, i - 1]]))))

# ----------------------------------------------------------------------------------------------------------------------
# 3. magnets
# ----------------------------------------------------------------------------------------------------------------------
def halbach_B(Br, tau, h, g=GAP, M=M_SEG):
    """peak field of a planar Halbach array at distance g from its face: Br (1 - e^-kh) e^-kg sin(pi/M)/(pi/M), k = pi/tau
    (Halbach 1980, Nucl. Instrum. Methods 169, 1; the segment factor for M blocks per wavelength)"""
    k = math.pi / tau
    return Br * (1 - math.exp(-k * h)) * math.exp(-k * g) * math.sin(math.pi / M) / (math.pi / M)

def magnet_table(name):
    mat = MAGNETS[name]; rows = []
    for tau in PITCHES:
        for h in THICKNESSES:
            B = halbach_B(mat["Br"], tau, h); p = B * B / (4 * MU0); kgm2 = mat["rho"] * h
            rows.append(dict(tau=tau, h=h, B=B, p_peak=p, kgm2=kgm2, kg_per_N=kgm2 / p))
    return rows

# ----------------------------------------------------------------------------------------------------------------------
# 4. the pod
# ----------------------------------------------------------------------------------------------------------------------
def pod_geometry(N, x_shield, m_magnet, w_track, h_mag):
    """sphere with a tube through it. x_shield: g/cm2 of water-equivalent shield; the aluminum hull counts at 0.8."""
    t_sh = max(x_shield - AL_AS_WATER * X_AL, 0.0) * 10.0 / RHO_SHIELD           # m of the water-polyethylene layer
    r_t = w_track / 2 + GAP + h_mag + TUBE_CLEAR                                   # bare tube radius (track, gaps, magnets, structure)
    r_ti = r_t + t_sh                                                              # tube with its shield sleeve
    V_c = V_SEAT * N
    R_in = math.sqrt((3 * V_c / (4 * math.pi))**(2.0 / 3.0) + r_ti**2)             # sphere with a cylindrical hole: V = 4/3 pi (R^2 - r^2)^1.5
    R_out = R_in + t_sh
    V_sh = 4.0 / 3.0 * math.pi * (R_out**2 - r_t**2)**1.5 - V_c
    m_sh = V_sh * RHO_SHIELD
    A_hull = 4 * math.pi * R_out**2 + 2 * math.pi * r_t * 2 * R_out
    m_hull = A_hull * (HULL_T_AL * RHO_AL + SKIN)
    m_car = m_magnet * (1 + CARRIAGE_FACTOR)
    m = M_SEAT * N + m_sh + m_hull + m_car + M_FIXED
    return dict(N=N, x=x_shield, t_sh=t_sh, r_t=r_t, R_in=R_in, R_out=R_out, V_c=V_c, m_sh=m_sh, A_hull=A_hull, m_hull=m_hull,
                m_mag=m_magnet, m_car=m_car, m_pax=M_SEAT * N, m=m, beta=m / (C_D * math.pi * R_out**2))

def size_pod(N, x_shield, mag):
    """iterate the pod mass with the magnets sized for the 1.5 g stop at the ring (force m (a_stop + anet) on a carriage L_CAR long)"""
    m = 20e3 * N / 25
    for _ in range(60):
        F_max = m * (A_STOP + anet(r_r))
        m_mag = F_max * mag["kg_per_N"]; L = max(L_CAR, F_max / (2 * W_MAX * mag["p_peak"])); w = F_max / (2 * L * mag["p_peak"])
        g = pod_geometry(N, x_shield, m_mag, w, mag["h"]); m_new = g["m"]
        if abs(m_new - m) < 0.1: break
        m = m_new
    g["F_max"] = F_max; g["A_face"] = F_max / mag["p_peak"]; g["w_track_force"] = w; g["L_car"] = L
    return g

# ----------------------------------------------------------------------------------------------------------------------
# 5. the track
# ----------------------------------------------------------------------------------------------------------------------
def r_grid(n=12000):
    x = np.linspace(0, 1, n)
    return r_r + (r_geo - r_r) * x**2                                              # dense near the ring

def track_profile(m, v_c, dH, p_peak, L=L_CAR, rg=None):
    """track per metre along the tether for one pod mass and cruise speed.
    Thickness from the characteristic speed: t = 2/(mu0 sigma v). Width: the larger of what the heat needs
    (F/(rho t dH), so one pass raises the ribbon by dH) and what the force needs at peak drag (F/(2 L p_peak))."""
    rg = r_grid() if rg is None else rg
    s = rg - r_r; a = anet(rg)
    in_stop = rg < r_r + v_c**2 / (2 * A_STOP)
    v = np.where(in_stop, np.sqrt(np.maximum(2 * A_STOP * s, 0.0)), v_c); v = np.maximum(v, V_HAND)
    F = np.where(in_stop, m * (a + A_STOP), m * a)
    t = 2.0 / (MU0 * SIG_CNT * v)
    w_heat = F / (RHO * t * dH); w_force = F / (2 * L * p_peak)
    w = np.maximum(w_heat, w_force); mu = RHO * t * w
    dH_local = F / mu                                                              # J/kg actually put into the ribbon per pass
    return dict(r=rg, alt=rg - RE, s=s, v=v, F=F, t=t, w=w, w_heat=w_heat, w_force=w_force, mu=mu, dH_local=dH_local, in_stop=in_stop)

def integrate_mu(tr, mask=None):
    mu = tr["mu"] if mask is None else np.where(mask, tr["mu"], 0.0)
    return float(np.trapezoid(mu, tr["r"]))

def track_cooling(t_mm, T_hot):
    """seconds for a ribbon t_mm thick radiating from both faces to fall from T_hot to 350 K (emissivity EPS_TRACK, 300 K surroundings)"""
    T = T_hot; tt = 0.0; dt = 0.01
    while T > 350.0:
        cp = (h_graphite(T + 1) - h_graphite(T - 1)) / 2.0
        T -= 2 * EPS_TRACK * SB * (T**4 - 300.0**4) / (RHO * t_mm * 1e-3 * cp) * dt; tt += dt
        if tt > 1e5: break
    return tt

# ----------------------------------------------------------------------------------------------------------------------
# 6. heat onto the magnets and into the water
# ----------------------------------------------------------------------------------------------------------------------
def T_from_dH(dh):
    """ribbon temperature after absorbing dh J/kg from T_TRACK0"""
    H = h_graphite(T_TRACK0) + dh
    if H > JANAF_H[-1] * 1e3 * (1 + 1e-6): return float("nan")
    return float(np.interp(min(H, JANAF_H[-1] * 1e3), JANAF_H * 1e3, JANAF_T))

def mean_flux(dh, T_mag=300.0):
    """mean radiant flux (W/m2) from the ribbon onto the magnet face over the carriage length, the ribbon warming linearly in
    enthalpy from T_TRACK0 at the front to T(dh) at the back; two large parallel faces, emissivities EPS_TRACK and EPS_MAG"""
    eps = 1.0 / (1.0 / EPS_TRACK + 1.0 / EPS_MAG - 1.0)
    f = np.linspace(0, 1, 41); T = np.array([T_from_dH(x * dh) for x in f])
    return float(eps * SB * np.mean(T**4 - T_mag**4)), float(eps * SB * (T[-1]**4 - T_mag**4))

def magnet_heat(d, tr, A_face):
    """heat into the magnet faces over one descent (J): the flux at each moment from the local ribbon enthalpy rise"""
    dh = np.interp(d["r"], tr["r"], tr["dH_local"]); dh = np.where(d["phase"] == 0, 0.0, dh)
    q = np.array([mean_flux(x)[0] if x > 0 else 0.0 for x in dh[::20]])            # every 20th step, then scaled
    Q = float(np.sum(q) * d["dt"] * 20 * A_face)
    q_peak = max(mean_flux(float(dh.max()))[1], 0.0)
    return Q, q_peak

# ----------------------------------------------------------------------------------------------------------------------
# 7. the tether
# ----------------------------------------------------------------------------------------------------------------------
def pod_envelope(d, rg, n_phase=120):
    """E(r): the largest sum, over the phase of the pod train, of the forces of the pods below radius r.
    Pods leave one per hour; each pulls m*anet in the cruise and m*(anet + a_stop) in the stop; in the drive the coils push it
    down and the reaction on the tether is upward, counted as zero. Also returns the weight-only envelope (Paul's 1.39 pods)."""
    E = np.zeros(len(rg)); Ew = np.zeros(len(rg)); m = d["m"]; E_tot = 0.0; Ew_tot = 0.0
    for ph in np.linspace(0, HEADWAY, n_phase, endpoint=False):
        e = np.zeros(len(rg)); ew = np.zeros(len(rg)); et = 0.0; ewt = 0.0
        for i in range(d["n_on_line"] + 1):
            ti = ph + HEADWAY * i
            if ti > d["total"]: break
            j = np.searchsorted(d["t"], ti); j = min(j, len(d["t"]) - 1)
            ri = d["r"][j]; F = d["F"][j]; Fw = m * anet(ri) if d["phase"][j] > 0 else 0.0
            below = rg > ri
            e[below] += F; ew[below] += Fw; et += F; ewt += Fw
        E = np.maximum(E, e); Ew = np.maximum(Ew, ew); E_tot = max(E_tot, et); Ew_tot = max(Ew_tot, ewt)
    return E, Ew, E_tot, Ew_tot

def size_tether(rg, mu_track, E):
    """march the tension up from the ring: dT/dr = (rho T/sigma + mu_track) anet + dE ; A = T/sigma.
    The pull of the upper tether at the platform is fixed, so the section at r carries the most when no pod is above it:
    T(r) = (largest total pod load on the line) + weight of tether and track below r. E is that total, constant in r;
    the ring holds the difference as pretension and must keep it above zero with every pod on the line.
    Returns the profile and the upper tether (constant section A_geo) long enough to balance T_geo with no end mass."""
    n = len(rg); T = np.zeros(n); T[0] = E[0]
    for k in range(n - 1):
        dr = rg[k + 1] - rg[k]
        T[k + 1] = T[k] + (RHO * T[k] / SIG + mu_track[k]) * anet(rg[k]) * dr + (E[k + 1] - E[k])
    A = T / SIG; m_lower = float(np.trapezoid(RHO * A, rg))
    m_up = RHO * A[-1]
    lo, hi = r_geo + 1e3, r_geo + 3e8
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if m_up * pull_above(mid) < T[-1]: lo = mid
        else: hi = mid
    L_up = hi - r_geo
    return dict(T=T, A=A, m_lower=m_lower, T_geo=T[-1], A_geo=A[-1], m_up_per_m=m_up, L_up=L_up, r_up_end=hi,
                m_upper=m_up * L_up, total=m_lower + m_up * L_up, T_ring=T[0], d_ring=2 * math.sqrt(A[0] / math.pi), d_geo=2 * math.sqrt(A[-1] / math.pi))

# ----------------------------------------------------------------------------------------------------------------------
# 8. free fall
# ----------------------------------------------------------------------------------------------------------------------
def freefall(alt0, R_nose, beta, vr0=0.0):
    """pod released from the tether at altitude alt0 (at rest on the tether, so with the tether's own tangential speed W r0,
    and a downward speed vr0 if it was moving). RK4 of freefall_trajectory.rk4_step with a step that shortens in the air."""
    r0 = RE + alt0; state = np.array([r0, 0.0, -vr0, W * r0]); t = 0.0
    eps = 0.5 * (vr0**2 + (W * r0)**2) - GM / r0; hh = r0 * W * r0
    if eps < 0:
        a = -GM / (2 * eps); e = math.sqrt(max(1 - hh**2 / (GM * a), 0.0)); r_p = a * (1 - e)
        if r_p > RE + 100e3:
            return dict(alt0=alt0, v_entry=None, perigee=r_p - RE, apogee=a * (1 + e) - RE, period=2 * math.pi * math.sqrt(a**3 / GM))
    q_tot = 0.0; q_peak = 0.0; g_peak = 0.0; v_entry = None; t_entry = None; v_peak = 0.0; alt_q = 0.0; t20 = 0.0; t40 = 0.0
    for _ in range(4_000_000):
        alt = state[0] - RE
        dt = 2.0 if alt > 400e3 else (0.2 if alt > 150e3 else 0.02)
        state = ff.rk4_step(state, dt, beta); t += dt
        r, th, vr, vt = state; alt = max(r - RE, 0.0)
        v_rel = math.sqrt(vr**2 + (vt - W * r)**2); v_peak = max(v_peak, v_rel)
        rho = ff.atmosphere_density(alt)
        q = tcfg.K_SUTTON_GRAVES * math.sqrt(rho / R_nose) * v_rel**3; q_tot += q * dt
        if q > q_peak: q_peak = q; alt_q = alt
        gg = 0.5 * rho * v_rel**2 / beta / g0; g_peak = max(g_peak, gg); t20 += dt * (gg > 20); t40 += dt * (gg > 40)
        if v_entry is None and alt <= 100e3: v_entry = v_rel; t_entry = t
        if alt <= 0: break
    return dict(alt0=alt0, v_entry=v_entry, t_entry=t_entry, q_peak=q_peak, alt_q=alt_q, q_total=q_tot, g_peak=g_peak, t20=t20, t40=t40,
                t_total=t, v_ground=v_rel, v_peak=v_peak, KE_entry=0.5 * v_entry**2)

# ----------------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------------
def main():
    print("elevator_system.py: the pod, shield, brake track, magnets and tether as one system")
    print(f"ring at {ring.ALTITUDE_KM:.0f} km (r = {r_r:,.0f} m), platform at GEO r = {r_geo:,.0f} m (altitude {(r_geo-RE)/1e3:,.0f} km), "
          f"net pull at the ring {anet(r_r):.3f} m/s2 ; CNT {RHO:.0f} kg/m3 at {SIG/1e9:.2f} GPa ; potential ring to platform {(Phi(r_geo)-Phi(r_r))/1e6:.2f} MJ/kg")
    print(f"the upper tether of constant section balances when its outward pull per kg/m reaches sigma/rho = {SIG/RHO/1e6:.2f} MJ/kg")
    rates = load_rates(); T_sh = rates["T"]

    # ---- 1. descent profiles (the accelerations do not depend on the pod mass; shown per tonne and for the pods below) ----
    section("1. DESCENT PROFILES: 1 g drive, cruise, 1.5 g stop")
    prof = {v: descent(v, 1000.0) for v in SPEEDS}
    for v, d in prof.items():
        print_descent(d, f"per tonne of pod at {v:,.0f} m/s")
    print("\nnet pull along the line, for the pod-train weights (one pod an hour at 2,000 m/s puts them 7,200 km apart):")
    for v in SPEEDS:
        d = prof[v]; tot = 0.0; parts = []
        for i in range(d["n_on_line"]):
            j = np.searchsorted(d["t"], d["total"] - HEADWAY * i); j = min(j, len(d["t"]) - 1)
            parts.append(f"{d['alt'][j]/1e3:,.0f} km: {anet(d['r'][j]):.3f}"); tot += anet(d["r"][j])
        print(f"  {v:,.0f} m/s, one pod arriving: pods at " + " ; ".join(parts) + f" m/s2 ; sum {tot:.2f} = {tot/anet(r_r):.2f} pod weights at the ring's pull")

    # ---- 2. dose and shield ----
    section("2. DOSE PER DESCENT AND THE SHIELD (AP-8 and AE-8 at solar max, longitude 0, folded with the Chapter 5 responses)")
    shield = {}
    for v in SPEEDS:
        d = prof[v]; p, e, tot = dose_per_descent(rates, d)
        print(f"\ncruise {v:,.0f} m/s: {d['total']/3600:.2f} h on the tether, {d['t_below_6000']/60:.0f} min below 6,000 km")
        print("  shield g/cm2 | water cm | protons mSv | electrons mSv | total mSv | other cases (total mSv): lon -30 ; lon +120 ; solar min lon 0")
        others = [dose_per_descent(rates, d, k)[2] for k in ("max_-30", "max_+120", "min_+0")]
        for i in range(len(T_sh)):
            if T_sh[i] > 150: break
            print(f"  {T_sh[i]:11g} | {T_sh[i]:8g} | {p[i]:11.3g} | {e[i]:13.3g} | {tot[i]:9.3g} | " + " ; ".join(f"{o[i]:.3g}" for o in others))
        for target in (20.0, 10.0, 5.0, 2.0, 1.0):
            x = shield_for(T_sh, tot, target); shield[(v, target)] = x
            print(f"  shield for {target:g} mSv per descent: {x:.1f} g/cm2 of water-equivalent = {x*10:.0f} kg/m2 ; "
                  f"as water-polyethylene behind the {HULL_T_AL*1e3:.0f} mm aluminum hull ({X_AL:.2f} g/cm2, counted as {AL_AS_WATER*X_AL:.2f}): {(x-AL_AS_WATER*X_AL):.1f} cm")
    print("\nNote: the dose is at the centre of a solid water sphere of the given radius (the Chapter 5 geometry); the pod's shield is a shell of that")
    print("thickness around the compartment, with the tube through the pod carrying the same thickness; the open ends of the tube are not modelled [UNVERIFIED].")

    # ---- 3. magnets ----
    section("3. HALBACH MAGNETS ON THE CARRIAGE: field at a 20 mm gap and the mass per newton of brake force")
    print("model: B = Br (1 - e^-kh) e^-kg sin(pi/M)/(pi/M), k = pi/pitch, h magnet thickness, g gap, M = 4 blocks per wavelength (Halbach 1980).")
    print("eddy drag per m2 of face at its peak B^2/(4 mu0), reached when the ribbon's characteristic speed 2/(mu0 sigma t) equals the pod speed;")
    print("the repulsion on each array is then also B^2/(4 mu0) per m2 (B^2/(2 mu0) in the high-speed limit), balanced by the array on the other side")
    print("(thin-sheet model: Reitz 1970, J. Appl. Phys. 41, 2067; Rossing and Hull 1991, The Physics Teacher 29, 552). Face area = both sides together.")
    best = {}
    for name in MAGNETS:
        rows = magnet_table(name); mat = MAGNETS[name]
        print(f"\n{name}: remanence {mat['Br']:.2f} T, {mat['rho']:.0f} kg/m3, limit {mat['T_max']}")
        print("  pitch mm | thickness mm | B at the track T | peak drag kPa | magnet kg/m2 | kg per MN of brake force | face m2 per MN")
        for r in rows:
            print(f"  {r['tau']*1e3:8.0f} | {r['h']*1e3:12.0f} | {r['B']:16.3f} | {r['p_peak']/1e3:13.1f} | {r['kgm2']:12.0f} | {r['kg_per_N']*1e6:24,.0f} | {1e6/r['p_peak']:14.1f}")
        b = min(rows, key=lambda r: r["kg_per_N"]); best[name] = b
        print(f"  lightest: pitch {b['tau']*1e3:.0f} mm, thickness {b['h']*1e3:.0f} mm, {b['B']:.3f} T, {b['p_peak']/1e3:.1f} kPa, {b['kg_per_N']*1e6:,.0f} kg per MN")
    k_opt = 1 / (2 * GAP); x_opt = 1.26
    print(f"\nanalytic optimum of h/B^2 for the ideal array: k = 1/(2g) = {k_opt:.1f} /m (pitch {math.pi/k_opt*1e3:.0f} mm), kh = {x_opt} (h = {x_opt/k_opt*1e3:.0f} mm); the sweep above brackets it")
    paul = dict(B=0.7, p_peak=0.7**2 / (4 * MU0))
    print(f"Paul's estimate point, 0.7 T: peak drag {paul['p_peak']/1e3:.1f} kPa ; 25 t pod at 1.5 g: force {25e3*(A_STOP+anet(r_r))/1e3:,.0f} kN, face {25e3*(A_STOP+anet(r_r))/paul['p_peak']:.1f} m2, "
          f"width on a 10 m carriage {25e3*(A_STOP+anet(r_r))/paul['p_peak']/20:.2f} m, peak power at 2,000 m/s {25e3*(A_STOP+anet(r_r+2000**2/(2*A_STOP)))*2000/1e9:.2f} GW")
    mag = dict(best["NdFeB"]); mag["name"] = "NdFeB"
    print(f"\ndesign magnets: {mag['name']}, pitch {mag['tau']*1e3:.0f} mm, thickness {mag['h']*1e3:.0f} mm, {mag['B']:.3f} T at the track, {mag['p_peak']/1e3:.1f} kPa peak drag, {mag['kg_per_N']*1e6:,.0f} kg per MN")
    smco = best["SmCo"]
    print(f"SmCo alternative: pitch {smco['tau']*1e3:.0f} mm, thickness {smco['h']*1e3:.0f} mm, {smco['B']:.3f} T, {smco['p_peak']/1e3:.1f} kPa, {smco['kg_per_N']*1e6:,.0f} kg per MN ({smco['kg_per_N']/mag['kg_per_N']:.2f} x the NdFeB mass)")
    print("thin-sheet check: skin depth at the pod speed must exceed the ribbon thickness; printed with the track below.")

    # ---- 4. pods ----
    section("4. THE PODS: sphere with a tube, seats, shield, hull, magnets, carriage (iterated with the magnets at the 1.5 g stop)")
    print(f"ASSUMPTIONS: hull {HULL_T_AL*1e3:.0f} mm aluminum ({HULL_T_AL*RHO_AL:.1f} kg/m2) ; skin {SKIN:.0f} kg/m2 ; shield layer as water at {RHO_SHIELD:.0f} kg/m3 ; "
          f"carriage structure = magnet mass ; {M_FIXED/1e3:.0f} t of loop, pumps, wheels, avionics ; tube clearance {TUBE_CLEAR*1e3:.0f} mm beyond the magnets ; carriage {L_CAR:.0f} m")
    pods = {}
    for target in (5.0, 2.0):
        for v in SPEEDS:
            x = shield[(v, target)]
            print(f"\nshield for {target:g} mSv at {v:,.0f} m/s: {x:.1f} g/cm2 ; water-polyethylene layer {(x-AL_AS_WATER*X_AL):.1f} cm")
            print("  seats | pod t | passengers+seats t | shield t | hull+skin t | magnets t | carriage t | outer dia m | compartment dia m | tube dia m | shield m2 | brake force at ring kN | face m2 | carriage m | track width m (5/10/20 m carriage) | start allowance pod t")
            for N in SEATS:
                g = size_pod(N, x, mag); pods[(N, v, target)] = g
                widths = "/".join(f"{g['F_max']/(2*L*mag['p_peak']):.2f}" for L in L_CARRIAGES)
                # the task's starting allowance: 8 t + 5% of the rest, with the same shield geometry
                g0_ = pod_geometry(N, x, 0.0, g["w_track_force"], mag["h"]); rest = g0_["m_pax"] + g0_["m_sh"] + g0_["m_hull"]
                m_start = rest + HARDWARE_START[0] + HARDWARE_START[1] * rest
                print(f"  {N:5d} | {g['m']/1e3:5.1f} | {g['m_pax']/1e3:18.1f} | {g['m_sh']/1e3:8.1f} | {g['m_hull']/1e3:11.1f} | {g['m_mag']/1e3:9.1f} | {g['m_car']/1e3:10.1f} | {2*g['R_out']:11.2f} | "
                      f"{2*g['R_in']:17.2f} | {2*g['r_t']:10.2f} | {g['A_hull']:9.0f} | {g['F_max']/1e3:22,.0f} | {g['A_face']:7.1f} | {g['L_car']:10.0f} | {widths:34s} | {m_start/1e3:.1f}")
    print("\nseats per tonne of pod and shield per seat:")
    for (N, v, target), g in sorted(pods.items(), key=lambda kv: (kv[0][2], kv[0][1], kv[0][0])):
        print(f"  {target:g} mSv, {v:,.0f} m/s, {N:3d} seats: {g['m']/1e3:6.1f} t, {N/(g['m']/1e3):.2f} seats per tonne, shield {g['m_sh']/N:.0f} kg per seat, magnets {g['m_mag']/g['m']*100:.0f}% of the pod")

    # ---- 5. the track ----
    section("5. THE BRAKE TRACK: a CNT ribbon graded so 2/(mu0 sigma t) = v, sized by the heat of each pass")
    print(f"heat per pass: energy per metre of track = brake force. Carbon enthalpy from NIST-JANAF graphite: {T_TRACK0:.0f} to 1,000 K {DH[1000.0]/1e3:.0f} kJ/kg, "
          f"{T_TRACK0:.0f} to 1,500 K {DH[1500.0]/1e3:,.0f} kJ/kg (Cp 8.5 J/mol/K at 298 K, 21.6 at 1,000 K, 23.9 at 1,500 K); CNT yarn taken as graphite [UNVERIFIED].")
    for v in SPEEDS + [4000.0]:
        print(f"  conductance for the characteristic speed {v:,.0f} m/s: sigma t = {2/(MU0*v):,.0f} S = {2/(MU0*v)/SIG_CNT*1e3:.2f} mm of {SIG_CNT/1e6:.0f} MS/m CNT")
    print(f"  at {V_HAND:.0f} m/s: {2/(MU0*V_HAND):,.0f} S = {2/(MU0*V_HAND)/SIG_CNT*1e3:.1f} mm ; the last {V_HAND**2/(2*A_STOP):.0f} m below {V_HAND:.0f} m/s are the wheels' (ASSUMPTION)")
    for v in SPEEDS:
        for tau in (0.1, 0.2):
            om = math.pi * v / tau; delta = math.sqrt(2 / (MU0 * SIG_CNT * om))
            print(f"  skin depth at {v:,.0f} m/s, pitch {tau*1e3:.0f} mm: {delta*1e3:.1f} mm vs ribbon {2/(MU0*v)/SIG_CNT*1e3:.2f} mm ; at {V_HAND:.0f} m/s: {math.sqrt(2/(MU0*SIG_CNT*math.pi*V_HAND/tau))*1e3:.0f} mm vs {2/(MU0*V_HAND)/SIG_CNT*1e3:.1f} mm")
    print(f"  cooling of the ribbon after a pass, both faces radiating at emissivity {EPS_TRACK} [UNVERIFIED at these temperatures]: "
          f"0.4 mm from 1,000 K to 350 K in {track_cooling(0.4, 1000.0):.1f} s, from 1,500 K in {track_cooling(0.4, 1500.0):.1f} s ; 8 mm from 1,000 K in {track_cooling(8.0, 1000.0):.0f} s, from 1,500 K in {track_cooling(8.0, 1500.0):.0f} s")
    print("  so the hour between pods is never a constraint; the ribbon's mass is set by the rise in one pass.")
    rg = r_grid()
    tracks = {}
    print("\ntrack mass: the lower bound is the descent energy over the enthalpy rise, m x 42.3 MJ/kg / dH; the width rule adds mass where the force needs more width than the heat")
    print("  mSv | m/s | seats | pod t | rise K | stop section t (last v^2/3g) | cruise section t (to the platform) | total t | x pod mass | lower bound m dPhi/dH t | widest point: m at km above the ring | width at the ring m | thickness at the ring mm | track kg/m at the ring, at 1,000 km, 5,000 km, 10,000 km, 20,000 km")
    for target in (2.0, 5.0):
        for v in SPEEDS:
            for N in SEATS:
                g = pods[(N, v, target)]
                for Tmax in (1000.0, 1500.0):
                    tr = track_profile(g["m"], v, DH[Tmax], mag["p_peak"], L=g["L_car"], rg=rg); tracks[(N, v, target, Tmax)] = tr
                    m_stop = integrate_mu(tr, tr["in_stop"]); m_cruise = integrate_mu(tr, ~tr["in_stop"]); m_tot = m_stop + m_cruise
                    lb = g["m"] * (Phi(r_geo) - Phi(r_r)) / DH[Tmax]
                    iw = int(np.argmax(tr["w"])); vals = [np.interp(RE + a, tr["r"], tr["mu"]) for a in (700e3 + 500.0, 1000e3, 5000e3, 10000e3, 20000e3)]
                    print(f"  {target:3g} | {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {Tmax-T_TRACK0:6.0f} | {m_stop/1e3:28,.1f} | {m_cruise/1e3:34,.0f} | {m_tot/1e3:7,.0f} | {m_tot/g['m']:10.1f} | {lb/1e3:24,.0f} | "
                          f"{tr['w'][iw]:.2f} at {tr['alt'][iw]/1e3-ring.ALTITUDE_KM:,.0f} | {np.interp(r_r+V_HAND**2/(2*A_STOP), tr['r'], tr['w']):.2f} | {np.interp(r_r+V_HAND**2/(2*A_STOP), tr['r'], tr['t'])*1e3:.1f} | "
                          + ", ".join(f"{x:.3f}" for x in vals))
    # profile of the 100-seat, 2 mSv, 2,000 m/s, 1,000 K track in detail
    g = pods[(100, 2000.0, 2.0)]; tr = tracks[(100, 2000.0, 2.0, 1000.0)]
    print(f"\ntrack profile, 100 seats ({g['m']/1e3:.1f} t), 2 mSv shield, 2,000 m/s, 1,000 K: altitude km | km above ring | pod speed m/s | brake force kN | thickness mm | width m (heat / force) | kg/m | rise K per pass")
    for s_km in [0.34, 1, 3, 10, 30, 60, 100, 136, 137, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, 35000]:
        rr = r_r + s_km * 1e3; j = np.searchsorted(tr["r"], rr); j = min(j, len(rg) - 1)
        print(f"  {tr['alt'][j]/1e3:9,.0f} | {s_km:7,.2f} | {tr['v'][j]:7,.0f} | {tr['F'][j]/1e3:8,.0f} | {tr['t'][j]*1e3:6.2f} | {tr['w'][j]:.3f} ({tr['w_heat'][j]:.3f} / {tr['w_force'][j]:.3f}) | {tr['mu'][j]:.4f} | {T_from_dH(tr['dH_local'][j])-T_TRACK0:.0f}")

    # ---- 6. magnet heat and water ----
    section("6. HEAT FROM THE TRACK ONTO THE MAGNETS AND INTO THE WATER")
    print(f"the ribbon under the carriage warms from {T_TRACK0:.0f} K at the front to its rise at the back; the face sees it across the {GAP*1e3:.0f} mm gap.")
    print(f"radiant exchange between two large parallel faces: eps_eff = 1/(1/{EPS_TRACK} + 1/{EPS_MAG} - 1) = {1/(1/EPS_TRACK+1/EPS_MAG-1):.3f} (ASSUMPTION: magnet face emissivity {EPS_MAG}; black track at 1 would see {1/(1/EPS_TRACK+1/1.0-1):.2f}).")
    for Tmax in (1000.0, 1500.0):
        qm, qb = mean_flux(DH[Tmax])
        print(f"  ribbon at {Tmax:.0f} K at the back: black-body {SB*Tmax**4/1e3:.0f} kW/m2 ; onto a {EPS_MAG} face: {qb/1e3:.1f} kW/m2 at the back of the carriage, {qm/1e3:.1f} kW/m2 averaged over its length")
        print(f"    conduction across a {mag['h']*1e3:.0f} mm NdFeB magnet at the back-of-carriage flux: {qb*mag['h']/K_MAG:.0f} K (ASSUMPTION k = {K_MAG} W/m/K); a bare black face would see {mean_flux(DH[Tmax])[1]/ (1/(1/EPS_TRACK+1/EPS_MAG-1)) * (1/(1/EPS_TRACK+1/1.0-1)) * mag['h'] / K_MAG:.0f} K")
    print("\n  mSv | m/s | seats | track K | heat into the magnets per descent MJ | passengers' heat MJ | water t | water rise K | wall (water) C | limit C | peak face flux kW/m2")
    heat = {}
    for target in (2.0, 5.0):
        for v in SPEEDS:
            for N in SEATS:
                g = pods[(N, v, target)]; d = descent(v, g["m"])
                for Tmax in (1000.0, 1500.0):
                    tr = tracks[(N, v, target, Tmax)]
                    Q, q_peak = magnet_heat(d, tr, g["A_face"]); Qp = Q_PERSON * N * d["total"]
                    dT = (Q + Qp) / (g["m_sh"] * C_W); heat[(N, v, target, Tmax)] = (Q, Qp, dT, q_peak)
                    print(f"  {target:3g} | {v:5,.0f} | {N:5d} | {Tmax:7.0f} | {Q/1e6:36.0f} | {Qp/1e6:19.0f} | {g['m_sh']/1e3:7.1f} | {dT:12.2f} | {T_W0+dT:14.1f} | {T_W_LIMIT:7.0f} | {q_peak/1e3:.1f}")
    print(f"\nASSUMPTIONS: water starts at {T_W0:.0f} C; the compartment wall is the water, and the limit is {T_W_LIMIT:.0f} C for air in a dry cabin (the critical wet-bulb limit of young adults")
    print("is 30.5 C, Vecellio et al. 2022, J. Appl. Physiol. 132, 340; a cabin dehumidified below that wet bulb keeps a 35 C wall tolerable, not comfortable). The magnets' own")
    print(f"temperature depends on the coolant path through the carriage, not modelled; the face flux and the conduction figures above bound it.")

    # ---- 7. the tether ----
    section("7. THE TETHER WITH THE PODS ON THE LINE AND THE TRACK AS DISTRIBUTED LOADS, AT 18.75 GPa")
    # check against the chapter: tip load only
    E_chk = np.full(len(rg), 50e3 * anet(r_r)); chk = size_tether(rg, np.zeros(len(rg)), E_chk)
    print(f"check, Table 13.3 method (50 t at the tip, no track, no pods on the line): tension at the ring {chk['T_ring']/1e3:.0f} kN, diameter {chk['d_ring']*1e3:.1f} mm, "
          f"at the platform {chk['T_geo']/1e6:.2f} MN and {chk['d_geo']*1e3:.1f} mm, lower tether {chk['m_lower']/1e3:,.0f} t ({chk['m_lower']/50e3:.0f} t per tonne), "
          f"upper {chk['m_up_per_m']:.3f} kg/m to {(chk['r_up_end']-RE)/1e3:,.0f} km ({chk['L_up']/1e3:,.0f} km, {chk['m_upper']/1e3:,.0f} t) ; whole {chk['total']/1e3:,.0f} t")
    print("  (Table 13.3: 396 kN, 5.2 mm, 18.3 MN, 35 mm, 36,700 t, to 82,000 km, 46,000 km, 76,300 t)")
    print(f"\nthe upper tether always reaches {(chk['r_up_end']-RE)/1e3:,.0f} km altitude ({chk['L_up']/1e3:,.0f} km long): the balance condition sigma/rho does not depend on the load.")
    print("\nloads: the largest sum over the phase of the pod train of the brake reactions of all the pods on the line (weight in the cruise, m (anet + 1.5 g) in the stop),")
    print("carried by every section because the upper tether's pull is fixed and the ring holds the slack as pretension; the track's weight; the tether's own weight.")
    print("The track above the crossover, where the tether's own mass per metre is 10 x the track's, is dropped (ASSUMPTION).")
    results = {}
    hdr = ("  mSv | m/s | seats | pod t | track K | pods on line | largest pod load on the line, pod weights | weight-only (Paul) | track t (cut at km) | pretension at ring kN | dia ring mm | "
           "tension at platform MN | dia platform mm | lower tether t | upper tether t | whole t | per tonne of pod | per seat on the line t | per passenger per day t | peak brake GW")
    for Tmax in (1000.0, 1500.0):
        for target in (2.0, 5.0):
            print(f"\n-- track to {Tmax:.0f} K, shield for {target:g} mSv --"); print(hdr)
            for v in SPEEDS:
                for N in SEATS:
                    g = pods[(N, v, target)]; d = descent(v, g["m"]); tr = tracks[(N, v, target, Tmax)]
                    Eb, Ewb, E_tot, Ew_tot = pod_envelope(d, rg); E = np.full(len(rg), E_tot)
                    res = size_tether(rg, tr["mu"], E)
                    # crossover: tether's own mass per metre >= 10 x track's; cut the track there and resize once
                    own = RHO * res["A"]; above = np.where(own >= tr["mu"] / TETHER_SHARE)[0]
                    r_cut = rg[above[0]] if len(above) else r_geo
                    mu_cut = np.where(rg < r_cut, tr["mu"], 0.0); m_track = float(np.trapezoid(mu_cut, rg))
                    res = size_tether(rg, mu_cut, E)
                    own = RHO * res["A"]; above = np.where(own >= tr["mu"] / TETHER_SHARE)[0]; r_cut = rg[above[0]] if len(above) else r_geo
                    mu_cut = np.where(rg < r_cut, tr["mu"], 0.0); m_track = float(np.trapezoid(mu_cut, rg)); res = size_tether(rg, mu_cut, E)
                    res.update(m_track=m_track, r_cut=r_cut, E_tot=E_tot, Ew_tot=Ew_tot, m_pod=g["m"], P_peak=d["P_stop_peak"], n_on_line=d["n_on_line"])
                    results[(N, v, target, Tmax)] = res
                    print(f"  {target:3g} | {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {Tmax:7.0f} | {d['n_on_line']:12d} | {E_tot/(g['m']*anet(r_r)):33.2f} | {Ew_tot/(g['m']*anet(r_r)):18.2f} | "
                          f"{m_track/1e3:6,.0f} ({(r_cut-RE)/1e3:,.0f}) | {res['T_ring']/1e3:18,.0f} | {res['d_ring']*1e3:11.1f} | {res['T_geo']/1e6:22.1f} | {res['d_geo']*1e3:15.0f} | "
                          f"{res['m_lower']/1e3:14,.0f} | {res['m_upper']/1e3:14,.0f} | {res['total']/1e3:7,.0f} | {res['total']/g['m']:17,.0f} | "
                          f"{res['total']/1e3/(N*d['n_on_line']):22,.1f} | {res['total']/1e3/(N*PODS_PER_DAY):22,.1f} | {d['P_stop_peak']/1e9:.2f}")
    # the same without the track (pods only) and without pods (track only), to show what each costs, for the 2 mSv, 1,000 K case
    print("\nwhat each load costs (2 mSv shield, 1,000 K track): whole tether in t for the pods only, the track only (with one pod's weight at the tip), both; 1,500 K track for comparison")
    print("  m/s | seats | pod t | pods only (in brackets: if each section carried only the pods below it) | track only (1,000 K) | both (1,000 K) | both (1,500 K) | both / pods only")
    for v in SPEEDS:
        for N in SEATS:
            g = pods[(N, v, 2.0)]; d = descent(v, g["m"]); Eb, _, E_tot, _ = pod_envelope(d, rg); E = np.full(len(rg), E_tot)
            a = size_tether(rg, np.zeros(len(rg)), E); a2 = size_tether(rg, np.zeros(len(rg)), Eb)
            res = results[(N, v, 2.0, 1000.0)]; mu_cut = np.where(rg < res["r_cut"], tracks[(N, v, 2.0, 1000.0)]["mu"], 0.0)
            b = size_tether(rg, mu_cut, np.full(len(rg), g["m"] * anet(r_r)))
            print(f"  {v:5,.0f} | {N:5d} | {g['m']/1e3:5.1f} | {a['total']/1e3:9,.0f} ({a2['total']/1e3:,.0f}) | {b['total']/1e3:20,.0f} | {res['total']/1e3:14,.0f} | {results[(N, v, 2.0, 1500.0)]['total']/1e3:14,.0f} | {res['total']/a['total']:.1f}")
    print("\ntether per passenger per day against pod size (whole tether / (24 pods x seats)), t:")
    print("  seats | " + " | ".join(f"{target:g} mSv {v:,.0f} m/s {Tmax:.0f} K" for Tmax in (1000.0, 1500.0) for target in (2.0, 5.0) for v in SPEEDS))
    for N in SEATS:
        print(f"  {N:5d} | " + " | ".join(f"{results[(N, v, target, Tmax)]['total']/1e3/(N*PODS_PER_DAY):22,.1f}" for Tmax in (1000.0, 1500.0) for target in (2.0, 5.0) for v in SPEEDS))
    print("fall from one pod size to the next (per passenger per day), 2 mSv, 1,000 K:")
    for v in SPEEDS:
        vals = [results[(N, v, 2.0, 1000.0)]["total"] / N for N in SEATS]
        print(f"  {v:,.0f} m/s: " + " ; ".join(f"{SEATS[i]}->{SEATS[i+1]} seats {100*(1-vals[i+1]/vals[i]):.0f}%" for i in range(len(SEATS) - 1)))

    # ---- 8. the slot ----
    section("8. THE SLOT: numbers for the note")
    g = pods[(100, 2000.0, 2.0)]
    print(f"100 seats, 2 mSv, 2,000 m/s: brake force at the ring {g['F_max']/1e3:,.0f} kN on {g['A_face']:.1f} m2 of face; repulsion on each side at peak drag {mag['p_peak']/1e3:.1f} kPa x {g['A_face']/2:.1f} m2 = {mag['p_peak']*g['A_face']/2/1e3:,.0f} kN,")
    print(f"  {2*mag['p_peak']*g['A_face']/2/1e3:,.0f} kN in the high-speed limit; a C-shaped carriage carries that as a bending moment across its opening, two hinged halves carry it as tension in the hinge and latch.")
    print(f"  the ribbon at the ring is {np.interp(r_r+V_HAND**2/(2*A_STOP), tracks[(100,2000.0,2.0,1000.0)]['r'], tracks[(100,2000.0,2.0,1000.0)]['w']):.2f} m wide for the force; the eddy loops close within about a pole pitch ({mag['tau']*1e3:.0f} mm) of the magnet edge, "
          f"so the ribbon needs that margin on the free edge and the tether bonded along the other edge [UNVERIFIED margin].")

    # ---- 9. loss of the tether ----
    section("9. LOSS OF THE TETHER: free fall of the pod, released at rest on the tether (tangential speed W r), sphere of the pod's size")
    print(f"heat flux: Sutton-Graves with the pod's radius as the nose radius; heat into the vehicle {F_HEAT[0]*100:.0f} to {F_HEAT[1]*100:.0f}% of the kinetic energy at 100 km (Paul, unverified);")
    print(f"water boiled at {H_BOIL/1e6:.2f} MJ/kg from 20 C. Ballistic coefficient m/(C_D A) with C_D {C_D}. Peak g is the drag deceleration (no lift).")
    surv = {}
    for (N, v, target) in ((25, 2000.0, 2.0), (100, 2000.0, 2.0), (400, 2000.0, 2.0), (100, 2000.0, 5.0)):
        g = pods[(N, v, target)]
        print(f"\n{N} seats, {target:g} mSv shield, pod {g['m']/1e3:.1f} t, radius {g['R_out']:.2f} m, water {g['m_sh']/1e3:.1f} t, beta {g['beta']:,.0f} kg/m2")
        print("  release km | fall to 100 km | entry speed at 100 km m/s | peak heat flux MW/m2 (at km) | heat load at the stagnation point MJ/m2 | kinetic energy at entry GJ | water boiled at 2% / 3% t | share of the water | peak g (s above 20 g, above 40 g) | speed at the ground m/s (no parachute)")
        for alt0 in RELEASE_ALTS:
            f = freefall(alt0, g["R_out"], g["beta"]); surv[(N, target, alt0)] = f
            if f["v_entry"] is None:
                print(f"  {alt0/1e3:10,.0f} | stays in orbit: perigee {f['perigee']/1e3:,.0f} km, apogee {f['apogee']/1e3:,.0f} km, period {f['period']/3600:.1f} h ; it never reaches the air"); continue
            KE = f["KE_entry"] * g["m"]; boiled = [F * KE / H_BOIL for F in F_HEAT]
            print(f"  {alt0/1e3:10,.0f} | {f['t_entry']/60:11.1f} min | {f['v_entry']:24,.0f} | {f['q_peak']/1e6:8.2f} ({f['alt_q']/1e3:.0f}) | {f['q_total']/1e6:38.0f} | {KE/1e9:26,.0f} | {boiled[0]/1e3:.1f} / {boiled[1]/1e3:.1f} | "
                  f"{boiled[0]/g['m_sh']*100:.0f} to {boiled[1]/g['m_sh']*100:.0f}% | {f['g_peak']:6.1f} ({f['t20']:.1f}, {f['t40']:.1f}) | {f['v_ground']:.0f}")
    lo, hi = 700e3, 35000e3
    for _ in range(60):
        mid = 0.5 * (lo + hi); r0 = RE + mid; eps = 0.5 * (W * r0)**2 - GM / r0; a = -GM / (2 * eps); e = math.sqrt(max(1 - (r0 * W * r0)**2 / (GM * a), 0.0))
        if a * (1 - e) > RE + 100e3: hi = mid
        else: lo = mid
    print(f"\nreleased at rest on the tether above {hi/1e3:,.0f} km the pod's perigee is above 100 km and it stays in orbit; below that it comes down.")
    g = pods[(100, 2000.0, 2.0)]
    f = freefall(700e3, g["R_out"], g["beta"], vr0=2000.0); fs = freefall(700e3 + 136e3, g["R_out"], g["beta"], vr0=2000.0)
    print(f"\n100 seats released while still moving at 2,000 m/s down: from 700 km entry {f['v_entry']:,.0f} m/s, {f['g_peak']:.0f} g, water boiled {F_HEAT[1]*f['KE_entry']*g['m']/H_BOIL/g['m_sh']*100:.0f}% at 3%; "
          f"from the top of the stop ({(700e3+136e3)/1e3:.0f} km) entry {fs['v_entry']:,.0f} m/s, {fs['g_peak']:.0f} g, {F_HEAT[1]*fs['KE_entry']*g['m']/H_BOIL/g['m_sh']*100:.0f}%")
    print("for the parachute and the splashdown the pod must also survive the peak g above. The reference for restrained passengers is Eiband's survey (NASA Memo 5-19-59E, 1959):")
    print("its curves allow a few tens of g only for fractions of a second [exact values UNVERIFIED here]; the peaks above last seconds, so the release altitudes that")
    print("the heat allows are limited by the deceleration first. A steep entry from a stationary release has no lift to spread it [QUESTION FOR AUTHOR: a lifting shape?].")

    # ---- 10. cargo drop pods ----
    section("10. CARGO DROP PODS: unshielded, dropped from the ring or from a trolley line, heat shield and parachute")
    print(f"ASSUMPTIONS: cargo at {CARGO_RHO:.0f} kg/m3 in a sphere; structure {CARGO_STRUCT*100:.0f}% of the cargo; heat-shield mass fraction from the correlation of Laub and Venkatapathy (2004),")
    print("TPS % of entry mass = 0.091 x (stagnation heat load in J/cm2)^0.51575 [UNVERIFIED correlation, fitted to orbital and planetary capsules];")
    A_chute = 2 * g0 / (1.225 * 0.8 * 8.0**2)
    print(f"parachute for an 8 m/s splashdown at C_D 0.8: {A_chute:.3f} m2 per kg of landed mass = {A_chute*1e3:.0f} m2 per tonne; at 60 g/m2 of canopy and lines and a 30% rigging allowance {A_chute*1e3*0.06*1.3:.0f} kg per tonne; {PARACHUTE_KG_PER_T:.0f} kg per tonne used.")
    print("\n  cargo t | radius m | beta kg/m2 | release | entry speed at 100 km m/s | peak heat flux MW/m2 | heat load MJ/m2 (J/cm2) | peak g | TPS % | TPS t | parachute t | empty pod t | empty pod per tonne of cargo | speed at the ground m/s")
    for M_c in (10e3, 50e3):
        R = (3 * M_c / CARGO_RHO / (4 * math.pi))**(1.0 / 3.0)
        for label, alt0 in (("ring 700 km", 700e3), ("trolley 300 km", 300e3), ("trolley 200 km", 200e3), ("trolley 100 km", 100e3)):
            m_entry = M_c * (1 + CARGO_STRUCT) * 1.15; beta = m_entry / (C_D * math.pi * R**2)
            f = freefall(alt0, R, beta)
            tps = 0.091 * (f["q_total"] / 1e4)**0.51575 / 100 if f["q_total"] > 0 else 0.0
            m_tps = tps * m_entry; m_chute = PARACHUTE_KG_PER_T * (M_c * (1 + CARGO_STRUCT) + m_tps) / 1e3
            m_empty = M_c * CARGO_STRUCT + m_tps + m_chute
            print(f"  {M_c/1e3:7.0f} | {R:8.2f} | {beta:10,.0f} | {label:14s} | {f['v_entry'] if f['v_entry'] else 0:24,.0f} | {f['q_peak']/1e6:19.2f} | {f['q_total']/1e6:8.1f} ({f['q_total']/1e4:,.0f}) | {f['g_peak']:6.1f} | {tps*100:5.1f} | {m_tps/1e3:5.1f} | {m_chute/1e3:11.2f} | {m_empty/1e3:11.1f} | {m_empty/M_c:27.2f} | {f['v_ground']:.0f}")
    print(f"  a trolley line carries {TROLLEY_LOAD/1e3:.0f} t at its lower end (Chapter 11), so the empty pods above ride up on a trolley car.")
    print("\ncargo descents at 3,000 and 4,000 m/s with the 1.5 g stop (an unshielded 50 t cargo pod of 55 t with its structure; the drive and the track scale with mass):")
    for v in (3000.0, 4000.0):
        d = descent(v, 55e3); print_descent(d, f"cargo pod at {v:,.0f} m/s")
    print("\ntrack consequences: the ribbon's conductance is graded for the passenger pods' speeds. On it a faster pod runs at u = v / v_char > 1, and the thin-sheet drag")
    print("falls to 2u/(1+u^2) of its peak, so the carriage needs more face (or a smaller gap) by the inverse:")
    for v in (3000.0, 4000.0):
        for v_pass in SPEEDS:
            s_stop_fast = v**2 / (2 * A_STOP); s_stop_pass = v_pass**2 / (2 * A_STOP)
            u_top = v / v_pass; frac = 2 * u_top / (1 + u_top**2)
            print(f"  cargo at {v:,.0f} m/s on the track graded for {v_pass:,.0f} m/s: its stop starts {s_stop_fast/1e3:,.0f} km up, the passenger stop {s_stop_pass/1e3:,.0f} km; in the cruise and down to {s_stop_pass/1e3:,.0f} km the cargo pod sees u = {u_top:.2f}, "
                  f"{frac*100:.0f}% of peak drag ({1/frac:.2f} x the face); below {s_stop_pass/1e3:,.0f} km both pods follow v = sqrt(3 g s) and the grading matches.")
    print("  the heat per metre is the brake force whatever the speed, so a 55 t cargo pod leaves 55/pod-mass of a passenger pod's energy per metre; the ribbon width must cover the heavier pod.")

    # ---- figures and CSVs ----
    section("FIGURES AND FILES")
    # descent profiles, Figure 13.2 form, for the 25-seat 2 mSv pod
    for v in SPEEDS:
        g = pods[(25, v, 2.0)]; d = descent(v, g["m"])
        fig, ax = plt.subplots(figsize=(8, 4.8)); ax2 = ax.twinx()
        ax.plot(d["s_km"], d["P"] / 1e6, color=RED, lw=1.6, label="brake power"); ax2.plot(d["s_km"], d["v"] * 3.6, color=BLUE, lw=1.6, label="speed")
        ax.set_xlabel("distance down the tether from the platform (km)"); ax.set_ylabel("brake power (MW)", color=RED); ax2.set_ylabel("speed (km/h)", color=BLUE)
        ax.set_xlim(0, (r_geo - r_r) / 1e3); ax.set_ylim(0, None); ax2.set_ylim(0, v * 3.6 * 1.15); ax.grid(alpha=0.3)
        j = int(np.argmax(d["P"])); ax.annotate(f"{d['P'][j]/1e6:,.0f} MW at the start of the stop", xy=(d["s_km"][j], d["P"][j] / 1e6), xytext=(0.45 * (r_geo - r_r) / 1e3, 0.85 * d["P"][j] / 1e6), color=RED, arrowprops=dict(arrowstyle="->", color=RED))
        ax.set_title(f"{g['m']/1e3:.0f} t pod (25 seats, 2 mSv shield): {v:,.0f} m/s cruise, 1 g drive, 1.5 g stop over {d['s_stop']/1e3:.0f} km", fontsize=10)
        fig.tight_layout(); fn = os.path.join(FIG, f"descent_profile_{v:.0f}.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
        with open(os.path.join(FIG, f"descent_profile_{v:.0f}.csv"), "w", newline="") as fh:
            wr = csv.writer(fh); wr.writerow(["t_s", "distance_from_platform_km", "altitude_km", "speed_m_s", "brake_force_kN", "brake_power_MW", "phase"])
            for k in range(0, len(d["t"]), 20): wr.writerow([f"{d['t'][k]:.1f}", f"{d['s_km'][k]:.1f}", f"{d['alt'][k]/1e3:.1f}", f"{d['v'][k]:.1f}", f"{d['F'][k]/1e3:.2f}", f"{d['P'][k]/1e6:.2f}", d["phase"][k]])
    # track mass profiles
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for N, c in zip((25, 100, 400), (BLUE, GREEN, RED)):
        for Tmax, ls in ((1000.0, "-"), (1500.0, "--")):
            tr = tracks[(N, 2000.0, 2.0, Tmax)]; res = results[(N, 2000.0, 2.0, Tmax)]
            m_ = np.where(tr["r"] < res["r_cut"], tr["mu"], np.nan)
            ax.plot((tr["r"] - r_r) / 1e3, m_, color=c, ls=ls, lw=1.4, label=f"{N} seats, {Tmax:.0f} K")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 5000); ax.set_xlabel("distance above the ring (km)"); ax.set_ylabel("track mass (kg/m)"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
    ax.set_title("brake track along the tether, 2,000 m/s cruise, 2 mSv shield", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "track_mass_profiles.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    with open(os.path.join(FIG, "track_profiles_2000.csv"), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["km_above_ring"] + [f"{N}seats_{Tmax:.0f}K_kg_per_m" for N in SEATS for Tmax in (1000.0, 1500.0)] + ["thickness_mm", "speed_m_s"])
        tr0 = tracks[(25, 2000.0, 2.0, 1000.0)]
        for k in range(0, len(rg), 40):
            wr.writerow([f"{(rg[k]-r_r)/1e3:.2f}"] + [f"{tracks[(N, 2000.0, 2.0, Tmax)]['mu'][k]:.5f}" for N in SEATS for Tmax in (1000.0, 1500.0)] + [f"{tr0['t'][k]*1e3:.3f}", f"{tr0['v'][k]:.0f}"])
    # tether per passenger per day vs seats
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for v, c in zip(SPEEDS, (BLUE, RED)):
        for Tmax, ls in ((1000.0, "-"), (1500.0, "--")):
            ax.plot(SEATS, [results[(N, v, 2.0, Tmax)]["total"] / 1e3 / (N * PODS_PER_DAY) for N in SEATS], color=c, ls=ls, marker="o", ms=4, lw=1.4, label=f"{v:,.0f} m/s, track to {Tmax:.0f} K")
    ax.set_xscale("log"); ax.set_xticks(SEATS); ax.set_xticklabels([str(s) for s in SEATS]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlabel("seats per pod"); ax.set_ylabel("tether per passenger per day (t)"); ax.set_ylim(0, None)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8); ax.set_title("whole tether per passenger per day, 2 mSv shield, 24 pods a day", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "tether_per_passenger.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    # survivability
    fig, ax = plt.subplots(figsize=(8, 4.8)); ax2 = ax.twinx()
    g = pods[(100, 2000.0, 2.0)]
    RELEASE_ALTS[:] = [a for a in RELEASE_ALTS if surv[(100, 2.0, a)]["v_entry"] is not None]
    alts = [a / 1e3 for a in RELEASE_ALTS]
    lo = [F_HEAT[0] * surv[(100, 2.0, a)]["KE_entry"] * g["m"] / H_BOIL / g["m_sh"] * 100 for a in RELEASE_ALTS]
    hi = [F_HEAT[1] * surv[(100, 2.0, a)]["KE_entry"] * g["m"] / H_BOIL / g["m_sh"] * 100 for a in RELEASE_ALTS]
    ax.fill_between(alts, lo, hi, color=BLUE, alpha=0.25, label="water boiled, 2 to 3% of the kinetic energy"); ax.plot(alts, hi, color=BLUE, lw=1.2); ax.plot(alts, lo, color=BLUE, lw=1.2)
    ax.axhline(100, color=GREY, ls=":", lw=1); ax.text(alts[0], 103, "all the water gone", color=GREY, fontsize=8)
    ax2.plot(alts, [surv[(100, 2.0, a)]["g_peak"] for a in RELEASE_ALTS], color=RED, marker="o", ms=4, lw=1.4, label="peak deceleration")
    ax.set_xscale("log"); ax.set_xlabel("release altitude (km)"); ax.set_ylabel("share of the pod's water boiled (%)", color=BLUE); ax2.set_ylabel("peak deceleration (g)", color=RED)
    ax.set_ylim(0, None); ax2.set_ylim(0, None); ax.grid(alpha=0.3, which="both"); ax.set_title(f"100-seat pod ({g['m']/1e3:.0f} t, {g['m_sh']/1e3:.0f} t of water) released from the tether; above 23,500 km it stays in orbit", fontsize=10)
    fig.tight_layout(); fn = os.path.join(FIG, "survivability.png"); fig.savefig(fn, dpi=200); plt.close(fig); print("wrote", fn)
    with open(os.path.join(FIG, "tether_vs_seats.csv"), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["mSv", "cruise_m_s", "track_K", "seats", "pod_t", "track_t", "track_cut_km", "lower_t", "upper_t", "whole_t", "per_seat_on_line_t", "per_passenger_per_day_t", "peak_brake_GW"])
        for (N, v, target, Tmax), res in sorted(results.items()):
            wr.writerow([target, v, Tmax, N, f"{res['m_pod']/1e3:.1f}", f"{res['m_track']/1e3:.0f}", f"{(res['r_cut']-RE)/1e3:.0f}", f"{res['m_lower']/1e3:.0f}", f"{res['m_upper']/1e3:.0f}", f"{res['total']/1e3:.0f}",
                         f"{res['total']/1e3/(N*res['n_on_line']):.1f}", f"{res['total']/1e3/(N*PODS_PER_DAY):.1f}", f"{res['P_peak']/1e9:.2f}"])
    print("done")

if __name__ == "__main__":
    main()
