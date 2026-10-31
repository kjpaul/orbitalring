"""
The anchor lines during the deployment (Chapter 11, Section 11.1): the casing's rate of slowing and its net weight day by day from
the deployment time series, the pull on a stowed line, the gradient term below the casing, the brake energy of paying out a 550 km
trolley line early against after the deployment (the line model written out below), and the day a line tip at 100 km first moves
through the air at less than 500 m/s.
Book: Orbital Ring Engineering, Volume I: Chapter 11, Section 11.1 (the deployment paragraphs); Section 11.3 (the line's taper).
Run: python tools/anchor_lines_deployment.py <deployment_tension csv>   (default ../../Vol_III/runs_700km_hw10686_2026-10-12/deployment_tension_700_brandt_flat.csv)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

The line model.  A trolley line for a 100 t load at its tip, the constant-stress taper of Section 11.3: the cross-section at height z
above the ground is A(z) = A_tip exp((rho / sigma) * integral_0^z a_net(z') dz') with a_net(z') = GM / r'^2 - W^2 r' at r' = R_E + z'
(the line at rest against the turning Earth), A_tip = 100,000 kg x 9.80665 m/s2 / 12.5 GPa, rho = 1,700 kg/m3; mass per metre
mu(z) = rho A(z).  The line is stowed on the casing and paid out tip first at u = 0.3 m/s, so when a length L is out the point at
depth d below the drum is the line's point z = L - d.  At the drum the tension is the weight of what is out:
    T(t) = integral_0^L mu(L - d) [a_c(t) + g'(t) d] dd  =  a_c(t) M(L) + g'(t) [L M(L) - S(L)],
    M(L) = integral_0^L mu(z) dz,  S(L) = integral_0^L z mu(z) dz,
with a_c(t) the casing's net weight per kilogram from the run (casing_net_weight_kN_m x 1000 / 12,000) and g'(t) = 2 GM / r^3 + omega(t)^2
the growth of that weight per metre of depth, omega(t) = |v_casing(t)| / r (the difference between gravity and the centripetal demand at
the casing's angular rate, one metre lower; at the start omega^2 = GM / r^3 and g' = 3 GM / r^3; at the end omega = W).  The casing's
slowing adds a horizontal pull a_dec(t) M(L); the drum sees the vector sum.  The brake absorbs E = integral T u dt = integral T dL.
After the deployment a_c = G_NET and omega = W with no horizontal pull.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys, csv, math
import numpy as np
import ring_altitude as ring

GM, RE, W, r_r, G_NET = ring.GM, ring.R_EARTH_EQ, ring.OMEGA_SIDEREAL, ring.R_ORBIT, ring.G_NET
RHO, SIGMA, G0 = 1700.0, 12.5e9, 9.80665
M_LOAD, M_TIP, L_LINE = 12000.0, 100e3, r_r - RE
U_PAYOUT, L_PAYOUT, M_LINE_BOOK = 0.3, 550e3, 148e3
DAY = 86400.0
src = sys.argv[1] if len(sys.argv) > 1 else _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "..", "Vol_III", "runs_700km_hw10686_2026-10-12", "deployment_tension_700_brandt_flat.csv")

# ---- the run ----
rows = list(csv.DictReader(open(src)))
t = np.array([float(r["time_days"]) for r in rows]); v = np.array([float(r["v_casing_m_s"]) for r in rows])   # launch frame, westward positive
w_kN = np.array([float(r["casing_net_weight_kN_m"]) for r in rows])
a_c = w_kN * 1e3 / M_LOAD                                   # m/s2, the casing's net weight per kilogram
def at(day, y): return float(np.interp(day, t, y))
t_end, v_end = float(t[-1]), float(v[-1])
print(f"Run: {_os.path.basename(src)} ({len(rows)} rows), the default deployment of 2026-10-27 (flat 3 mm winding, Brandt-Indenbom hysteresis, 8 MW per site; Python commit e758e24).")
print(f"End of run: day {t_end:.1f}, casing {-v_end:.1f} m/s eastward (ground-synchronous {ring.V_GROUND_SYNC:.1f}), total slowdown {v[0]-v_end:,.1f} m/s, average {(v[0]-v_end)/t_end:.1f} m/s per day.")
print()
print("1. THE CASING'S RATE OF SLOWING AND THE PULL ON A STOWED LINE")
for d in (1.0, 10.0):
    rate = (v[0] - at(d, v)) / d
    print(f"  over the first {d:.0f} day(s): {rate:.1f} m/s per day = {rate/DAY:.2e} m/s2; on a {M_LINE_BOOK/1e3:.0f} t stowed line {M_LINE_BOOK*rate/DAY:.0f} N")
rate1 = (v[0] - at(1.0, v)) / 1.0; a_dec1 = rate1 / DAY
print(f"  (text: 13 m/s per day, 1.5e-4 m/s2, 23 N)")
print()
print("2. THE GRADIENT TERM BELOW THE CASING")
gp0 = 3 * GM / r_r**3; gp_end = 2 * GM / r_r**3 + W**2
print(f"  a_net(r) = GM/r^2 - omega^2 r; one metre lower it is larger by -da/dr = 2 GM/r^3 + omega^2.")
print(f"  at the start omega = v_orbit/r, omega^2 = GM/r^3: g' = 3 GM/r^3 = 3 x {GM:.6g} / {r_r:.0f}^3 = {gp0:.3e} m/s2 per metre (text 3.4e-6)")
print(f"  at the end omega = W: g' = 2 GM/r^3 + W^2 = {2*GM/r_r**3:.3e} + {W**2:.3e} = {gp_end:.3e} m/s2 per metre")
print(f"  lean of a loose line ahead of the vertical on day 1 (tan = a_dec / (g' d), a_dec = {a_dec1:.2e}, g' at day 1 = {2*GM/r_r**3 + (at(1,v)/r_r)**2:.3e}):")
gp1 = 2 * GM / r_r**3 + (at(1.0, v) / r_r)**2
for dd in (100.0, 1e3, 10e3, 100e3):
    print(f"    {dd/1e3:>5.1f} km below: gradient pull {gp1*dd:.3e} m/s2, lean {math.degrees(math.atan2(a_dec1, gp1*dd)):.1f} deg (text: 25, 2.6, 0.25 deg at 100 m, 1 km, 10 km; 0.34 m/s2 at 100 km)")
print()
print("3. THE CASING'S NET WEIGHT PER KILOGRAM FROM THE RUN")
for d in (1, 7, 30, 100):
    print(f"  day {d:>3}: {at(d, a_c):.3f} m/s2 (text: 0.03, 0.2, 0.8, 2.6)")
d_tenth = float(np.interp(0.1 * G_NET, a_c, t))
d_half = float(np.interp(0.5 * G_NET, a_c, t))
print(f"  a tenth of the final {G_NET:.3f} m/s2 ({0.1*G_NET:.3f}) on day {d_tenth:.1f} (text: day 30); half on day {d_half:.1f}")
print()

# ---- the line ----
z = np.linspace(0.0, L_LINE, 70001)                          # 10 m steps from the ground to the ring
a_net_z = GM / (RE + z)**2 - W**2 * (RE + z)
I = np.concatenate([[0.0], np.cumsum(0.5 * (a_net_z[1:] + a_net_z[:-1]) * np.diff(z))])
A = (M_TIP * G0 / SIGMA) * np.exp(RHO / SIGMA * I)
mu = RHO * A
M = np.concatenate([[0.0], np.cumsum(0.5 * (mu[1:] + mu[:-1]) * np.diff(z))])
S = np.concatenate([[0.0], np.cumsum(0.5 * (mu[1:] * z[1:] + mu[:-1] * z[:-1]) * np.diff(z))])
print("4. THE LINE (Section 11.3 reproduced)")
print(f"  integral of a_net from the ground to the ring {I[-1]:.4e} m2/s2 (text 6.155e6); exponent rho/sigma x that = {RHO/SIGMA*I[-1]:.3f} (text 0.837); taper ratio {A[-1]/A[0]:.3f} (text 2.31)")
print(f"  A_tip = {M_TIP:.0f} x {G0} / {SIGMA:.3g} = {A[0]:.3e} m2, {2e3*math.sqrt(A[0]/math.pi):.1f} mm across (text 10 mm); at the ring {A[-1]:.3e} m2, {2e3*math.sqrt(A[-1]/math.pi):.1f} mm (text 15 mm)")
print(f"  mass {M[-1]/1e3:.1f} t (text 148 t); mass per metre {mu[0]:.3f} kg/m at the tip, {mu[-1]:.3f} at the ring; the first 550 km from the tip {np.interp(L_PAYOUT, z, M)/1e3:.1f} t")
print()

# ---- the payout ----
def ML(L): return float(np.interp(L, z, M))
def SL(L): return float(np.interp(L, z, S))
def tension(L, ac, gp, adec=0.0):
    Tv = ac * ML(L) + gp * (L * ML(L) - SL(L)); Th = adec * ML(L)
    return math.hypot(Tv, Th), Tv, Th
def payout(t0=None, n=20000):
    """brake energy and peak drum tension for paying out L_PAYOUT at U_PAYOUT from day t0 (None: after the deployment)"""
    Ls = np.linspace(0.0, L_PAYOUT, n + 1); E = 0.0; Tprev = None; Tmax = 0.0; parts = np.zeros(2)
    for k, L in enumerate(Ls):
        if t0 is None:
            ac, gp, adec = G_NET, gp_end, 0.0
        else:
            day = t0 + L / U_PAYOUT / DAY
            ac = at(day, a_c); om = at(day, v) / r_r; gp = 2 * GM / r_r**3 + om * om
            adec = (at(day - 0.05, v) - at(day + 0.05, v)) / (0.1 * DAY)
        T, Tv, Th = tension(L, ac, gp, adec); Tmax = max(Tmax, T)
        if Tprev is not None:
            E += 0.5 * (T + Tprev) * (Ls[k] - Ls[k - 1])
            parts += 0.5 * np.array([ac * ML(L) + acp * ML(Ls[k-1]), gp * (L * ML(L) - SL(L)) + gpp * (Ls[k-1] * ML(Ls[k-1]) - SL(Ls[k-1]))]) * (Ls[k] - Ls[k - 1])
        Tprev, acp, gpp = T, ac, gp
    return E, Tmax, parts, (t0 + L_PAYOUT / U_PAYOUT / DAY) if t0 is not None else None

dur = L_PAYOUT / U_PAYOUT / DAY
print("5. PAYING OUT 550 km OF THE LINE AT 0.3 m/s")
print(f"  payout takes {L_PAYOUT/U_PAYOUT/DAY:.1f} days (text: three weeks); the tension at the drum is T = a_c M(L) + g' [L M(L) - S(L)] (vertical) with the horizontal a_dec M(L) added as a vector.")
E_after, T_after, p_after, _ = payout(None)
print(f"  after the deployment (a_c = {G_NET:.3f} m/s2, g' = {gp_end:.3e}): brake energy {E_after/3.6e6:,.0f} kWh = {E_after/3.6e9:.1f} MWh = {E_after/1e9:.1f} GJ (text 78 MWh); of which the casing-weight term {p_after[0]/3.6e9:.1f} MWh and the gradient term {p_after[1]/3.6e9:.2f} MWh; peak drum tension {T_after/1e3:.1f} kN at the end of the payout (the line out weighs {ML(L_PAYOUT)/1e3:.1f} t)")
for t0 in (0.0, 1.0):
    E, Tm, p, tend = payout(t0)
    print(f"  from day {t0:.0f} (ends day {tend:.1f}, casing weight then {at(tend, a_c):.3f} m/s2): brake energy {E/3.6e6:,.0f} kWh = {E/3.6e9:.2f} MWh = {E/1e9:.2f} GJ (text 12,000 kWh); casing-weight term {p[0]/3.6e6:,.0f} kWh, gradient term {p[1]/3.6e6:,.0f} kWh; peak drum tension {Tm/1e3:.2f} kN; ratio after/early {E_after/E:.1f} (text: about seven)")
print("  start-day sweep (energy of the whole payout; the payout must end before the casing weight passes a tenth of its final value):")
best = None
for t0 in (0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 200, 300, 400, 478):
    E, Tm, p, tend = payout(t0, n=4000)
    ok = tend <= d_tenth
    print(f"    start day {t0:>3}: ends day {tend:6.1f}, energy {E/3.6e6:>9,.0f} kWh, peak tension {Tm/1e3:7.2f} kN{'  (ends before the tenth)' if ok else ''}")
    if ok and (best is None or E < best[1]): best = (t0, E)
print(f"  the tenth is passed on day {d_tenth:.1f} and the payout takes {dur:.1f} days, so {'the latest start that ends before it is day ' + str(best[0]) if best else 'no start day ends before it: a payout at 0.3 m/s cannot finish while the weight is still under a tenth of its final value'}.")
print(f"  the weight only grows, so the smallest energy is the earliest start: day 0 of the deployment.")
t0 = None
# the speed that would finish before the tenth
u_need = L_PAYOUT / (d_tenth * DAY)
print(f"  to finish 550 km by day {d_tenth:.1f} the payout would have to run at {u_need:.2f} m/s.")
print()

# ---- the air ----
print("6. THE AIR")
r100 = RE + 100e3; r150 = RE + 150e3
v_air100 = W * r100; v_air150 = W * r150
print(f"  air at 100 km co-rotating with the ground: W x {r100:.0f} = {v_air100:.1f} m/s eastward (465.1 x 6478/6378 = {465.1*6478/6378:.1f}); at 150 km {v_air150:.1f} m/s")
# tip at the casing's angular rate: westward speed v * r_tip / r; relative to the eastward air: v_tip + v_air
v_rel100 = v * r100 / r_r + v_air100
v_rel150 = v * r150 / r_r + v_air150
day500 = float(np.interp(-500.0, -v_rel100, t)) if v_rel100[0] > 500 else 0.0
v_cas_500 = (500.0 - v_air100) * r_r / r100
print(f"  a tip at 100 km moving at the casing's angular rate: v_tip = v_casing x {r100:.0f}/{r_r:.0f}; relative to the air v_tip + {v_air100:.1f}.")
print(f"  below 500 m/s when the casing is slower than {v_cas_500:.1f} m/s westward in the inertial frame, which is within {v_cas_500 - v_end:.0f} m/s of its final {-v_end:.1f} m/s eastward (text: about 550 m/s).")
print(f"  that happens on day {day500:.1f}; {t_end - day500:.1f} days of the deployment remain (text: the last 41 days; the brief's reading of the CSV: day 429, 70 days).")
for dd in (400, 429, 450, 470, 490, 499):
    print(f"    day {dd}: casing {at(dd, v):8.1f} m/s westward, tip at 100 km {at(dd, v_rel100):7.1f} m/s through the air, at 150 km {at(dd, v_rel150):7.1f} m/s")
print(f"  at the end of the deployment the tip moves with the air at every altitude ({at(t_end, v_rel100):.1f} m/s at 100 km, {at(t_end, v_rel150):.1f} at 150 km), so the relative speed falls monotonically to zero and the least heating is at the end.")
print("  No model of the heating of a line tip in the air exists in the repository, so the day of least heating at 150 km cannot be computed beyond that.")
