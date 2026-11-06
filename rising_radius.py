#!/usr/bin/env python3
"""
The ring's radius as a state variable of the LIM deployment: the rising-radius mode of lim_simulation.py (2026-11-06).
Book: Orbital Ring Engineering, Volume I: Chapter 3 (Section 3.2, the strain and the build orbit; Figure 3.2), Chapter 9
(the deployment run), Chapter 11 (Section 11.1, the anchor lines in the final tensioning), Chapter 12 (Section 12.11).
Run: python lim_simulation.py --rising-radius      (or ORBITAL_RING_RISING=1; the default run is unchanged)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT

The physics.  A closed loop is longer by its strain, C = C0 (1 + eps), eps = T / (E A).  The cable can only be assembled
unstressed, so a ring that operates at r = R_ORBIT with eps_d = sigma / E is built at r0 = R_ORBIT / (1 + eps_d) and rises as
the tension builds.  The LIM sites are fixed in number; they move angular momentum between the cable and the casing, which
is conserved (not the linear momentum per meter of the fixed-radius run):

    dL_c = +F r dt,   dL_l = -F r dt            (F the total LIM thrust around the ring, westward positive for the retrograde cable)
    v_c = L_c / (M_r r),   v_l = L_l / (M_l r)     (M_r, M_l the fixed total masses of the cable with hardware and of the casing)
    T = m_r (v_c^2 - GM / r) - m_l (GM / r - v_l^2),   m = M / (2 pi r)     (radial equilibrium at r; the casing's own centrifugal term)
    eps = T / (E A),   r = r0 (1 + eps)            (A the structural CNT area at the design point)

With L_c and L_l given, T(r) = a / r^3 - b / r^2 with a = (L_c^2 / M_r + L_l^2 / M_l) / (2 pi) and b = GM (M_r + M_l) / (2 pi), and
r = r0 (1 + T(r) / (E A)) is solved by Newton's method each step (quasi-static: the radial period is of the order of the orbit,
the deployment hundreds of days).  Phase 2 begins when the casing reaches the ground-synchronous speed, v_l = -omega_E r.  The
anchor lines then hold it there, L_l = -M_l omega_E r^2, the LIMs keep pushing the cable, and the difference between what the
LIM reaction takes from the casing and what the casing's momentum requires is supplied by the Earth: the external torque
tau = dL_E / dt, the anchor reaction tau / r in total and over round(L_RING / 50 km) stations.  Phase 2 ends when the cable
reaches the design speed or the tension reaches the design tension, whichever comes first.

Energy.  E = 1/2 M_r v_c^2 + 1/2 M_l v_l^2 - GM (M_r + M_l) / r + 1/2 E A eps^2 C0 (the elastic energy of a linear loop
stretched by eps C0, which equals 1/2 T eps C0).  Because the radius is at equilibrium, dE/dr = 0 at fixed momenta, so the
balance dE = F (v_c - v_l) dt + tau Omega_l dt holds exactly in continuous time: the LIM work (thrust times relative speed, as
in the fixed-radius run) plus the Earth's work through the anchors, tau x (-omega_E), the casing turning eastward with the
ground while the Earth pulls it westward.  The run prints both sides and their difference.
"""
import csv
import json
import math
import os

import ring_altitude as ring

E_CNT = 500e9                  # Pa, the cable's modulus (Chapter 3: 25 GPa at 5 % break strain; sled_passage_simulation.E_CNT)
ANCHOR_STATION_SPACING = 50e3  # m of ring per anchor station (the brief of 2026-11-06: round(L_RING / 50e3) stations)
DAY = 86400.0
TWO_PI = 2.0 * math.pi


class RisingRing:
    """The ring's angular momenta and radius through the deployment; one instance per run."""

    def __init__(self, m_rotor, m_load, m_struct, lim_sites, lims_per_site, sigma, v_cable_design=None, tension_design=None, E=E_CNT):
        self.GM, self.omega = ring.GM, ring.OMEGA_SIDEREAL
        self.R_design, self.L_design = ring.R_ORBIT, ring.L_RING
        self.M_r = m_rotor * self.L_design                 # total masses, fixed through the run
        self.M_l = m_load * self.L_design
        self.m_rotor, self.m_load = m_rotor, m_load         # per meter at the design point (for the record)
        self.A = m_struct / ring.RHO_CNT                    # structural CNT area at the design point
        self.E, self.EA = E, E * self.A
        self.sigma = sigma
        self.T_design = tension_design if tension_design is not None else sigma * self.A
        self.eps_design = self.T_design / self.EA
        if v_cable_design is None:                           # momentum conservation at the design radius, as ring_altitude.V_CABLE
            v_cable_design = ring.V_ORBIT + m_load * ring.DELTA_V_LOAD / m_rotor
        self.v_cable_design = v_cable_design
        self.r0 = self.R_design / (1.0 + self.eps_design)
        self.C0 = TWO_PI * self.r0
        self.v0 = math.sqrt(self.GM / self.r0)
        self.lim_sites, self.lims_per_site = lim_sites, lims_per_site
        self.n_stations = round(self.L_design / ANCHOR_STATION_SPACING)
        # state
        self.L_c = self.M_r * self.r0 * self.v0
        self.L_l = self.M_l * self.r0 * self.v0
        self.r, self.T, self.eps = self.r0, 0.0, 0.0
        self.v_c = self.v_l = self.v0
        self.phase = 1
        self.newton_iters_max = 0
        # the Earth, in phase 2
        self.L_earth = 0.0; self.W_earth = 0.0
        self.F_anchor = 0.0; self.F_anchor_peak = (0.0, 0.0); self.F_anchor_int = 0.0; self.t_phase2_dur = 0.0
        self.t_phase2 = None; self.r_phase2 = None; self.v_c_phase2 = None; self.T_phase2 = None; self.L_c_phase2 = None
        self.t_v_design = None; self.t_T_design = None; self.end_trigger = None
        self.t_end = None
        self.F_last = 0.0
        # bookkeeping
        self.W_lim = 0.0
        self.E0 = sum(self.energy())
        self.rows = []; self.last_day = -1
        self._log_day(0.0, 0.0, 0.0)
        # the sampled series in step with lim_simulation's data lists
        self.s_alt = [self.altitude_km]; self.s_strain = [0.0]; self.s_phase = [1]; self.s_anchor = [0.0]; self.s_r = [self.r]

    # ------------------------------------------------------------------ state
    @property
    def altitude_km(self):
        return (self.r - ring.R_EARTH_EQ) / 1e3

    def running(self):
        return self.phase < 3

    def progress(self):
        """Percent of the cable's speed gain from the build orbit to the design speed (a display figure)."""
        return 100.0 * (self.v_c - self.v0) / (self.v_cable_design - self.v0)

    def _solve(self):
        """Radial equilibrium at the current momenta: r = r0 (1 + T(r) / EA) by Newton's method; sets r, T, eps, v_c, v_l."""
        a_c = self.L_c ** 2 / self.M_r
        b = self.GM * (self.M_r + self.M_l) / TWO_PI
        if self.phase == 1:
            a = (a_c + self.L_l ** 2 / self.M_l) / TWO_PI; c = 0.0
        else:                                                # the casing held at omega_E: L_l^2 / M_l = M_l omega^2 r^4
            a = a_c / TWO_PI; c = self.M_l * self.omega ** 2 / TWO_PI
        r = self.r
        for i in range(100):
            T = a / r ** 3 + c * r - b / r ** 2
            dT = -3.0 * a / r ** 4 + c + 2.0 * b / r ** 3
            f = r - self.r0 * (1.0 + T / self.EA)
            df = 1.0 - self.r0 * dT / self.EA
            dr = f / df
            r -= dr
            if abs(dr) < 1e-6:                               # metres; the round-off floor of r is about 1e-9 m
                break
        self.newton_iters_max = max(self.newton_iters_max, i + 1)
        self.r = r
        self.T = a / r ** 3 + c * r - b / r ** 2
        self.eps = self.T / self.EA
        if self.phase == 2:
            self.L_l = -self.M_l * self.omega * r ** 2
        self.v_c = self.L_c / (self.M_r * r)
        self.v_l = self.L_l / (self.M_l * r)

    def _earth(self, dL_E, dt, t_new):
        self.L_earth += dL_E
        self.W_earth += -self.omega * dL_E                   # torque x Omega_l, Omega_l = -omega_E in the launch frame
        self.F_anchor = dL_E / (self.r * dt)
        if self.F_anchor > self.F_anchor_peak[0]:
            self.F_anchor_peak = (self.F_anchor, t_new)
        self.F_anchor_int += self.F_anchor * dt
        self.t_phase2_dur += dt

    def step(self, F_total, dt, t_new, v_rel):
        """Apply the LIM thrust (total around the ring, N) over dt; returns (v_c, v_l) in the launch frame."""
        dL = F_total * self.r * dt
        self.F_last = F_total
        self.W_lim += F_total * v_rel * dt
        self.L_c += dL
        if self.phase == 1:
            self.L_l -= dL
            self._solve()
            if self.v_l <= -self.omega * self.r:             # the casing at the ground-synchronous speed: the anchors take over
                self.phase = 2
                self.t_phase2, self.r_phase2, self.v_c_phase2, self.T_phase2, self.L_c_phase2 = t_new, self.r, self.v_c, self.T, self.L_c
                L_l_free = self.L_l
                self._solve()
                self._earth(self.L_l - L_l_free, dt, t_new)
        else:
            L_l_free = self.L_l - dL
            self._solve()
            self._earth(self.L_l - L_l_free, dt, t_new)
            if self.t_v_design is None and self.v_c >= self.v_cable_design:
                self.t_v_design = t_new
            if self.t_T_design is None and self.T >= self.T_design:
                self.t_T_design = t_new
            if self.t_v_design is not None or self.t_T_design is not None:
                self.phase = 3
                self.t_end = t_new
                both = self.t_v_design is not None and self.t_T_design is not None
                self.end_trigger = ("both in the same step" if both else
                                    ("the design cable speed" if self.t_v_design is not None else "the design tension"))
        self._log_day(t_new, F_total, dt)
        return self.v_c, self.v_l

    def forces(self):
        """(tension N, casing net weight N/m, cable excess force N/m, stress Pa) at the current radius; T / r = excess - weight."""
        r = self.r; g = self.GM / r ** 2
        m_r = self.M_r / (TWO_PI * r); m_l = self.M_l / (TWO_PI * r)
        W = m_l * (g - self.v_l ** 2 / r)
        F = m_r * (self.v_c ** 2 / r - g)
        return self.T, W, F, self.T / self.A

    def energy(self):
        """(cable kinetic, casing kinetic, potential, elastic) J, totals for the ring."""
        return (0.5 * self.M_r * self.v_c ** 2, 0.5 * self.M_l * self.v_l ** 2,
                -self.GM * (self.M_r + self.M_l) / self.r, 0.5 * self.EA * self.eps ** 2 * self.C0)

    def lim_work_required(self):
        """The LIM work the mechanical energy change implies: delta E minus the Earth's work."""
        return sum(self.energy()) - self.E0 - self.W_earth

    def sample(self):
        self.s_alt.append(self.altitude_km); self.s_strain.append(self.eps); self.s_phase.append(self.phase)
        self.s_anchor.append(self.F_anchor if self.phase >= 2 else 0.0); self.s_r.append(self.r)

    def _log_day(self, t, F_total, dt):
        day = int(t // DAY)
        if day <= self.last_day and t > 0.0:
            return
        self.last_day = day
        T, W, F_exc, S = self.forces()
        KE_c, KE_l, PE, U = self.energy()
        self.rows.append(dict(day=t / DAY, r_m=self.r, altitude_km=self.altitude_km, tension_GN=T / 1e9, stress_GPa=S / 1e9,
                              strain=self.eps, v_cable_m_s=self.v_c, v_casing_m_s=self.v_l, v_rel_m_s=self.v_c - self.v_l,
                              casing_net_weight_kN_m=W / 1e3, cable_excess_kN_m=F_exc / 1e3, phase=min(self.phase, 2),
                              lim_thrust_total_MN=F_total / 1e6, anchor_force_total_MN=(self.F_anchor / 1e6 if self.phase >= 2 else 0.0),
                              anchor_force_per_station_MN=(self.F_anchor / 1e6 / self.n_stations if self.phase >= 2 else 0.0),
                              L_cable=self.L_c, L_casing=self.L_l, L_earth=self.L_earth,
                              KE_cable_EJ=KE_c / 1e18, KE_casing_EJ=KE_l / 1e18, PE_EJ=PE / 1e18, elastic_EJ=U / 1e18,
                              lim_work_EJ=self.W_lim / 1e18, earth_work_EJ=self.W_earth / 1e18))

    # ------------------------------------------------------------------ output
    def print_build(self):
        print("\n" + "=" * 70)
        print("RISING RADIUS: the ring is built unstressed and rises with the cable's strain")
        print("=" * 70)
        print(f"Design point r = {self.R_design:,.0f} m ({(self.R_design - ring.R_EARTH_EQ) / 1e3:,.1f} km), design tension {self.T_design / 1e9:,.1f} GN on"
              f" A = {self.A:.2f} m^2 at E = {self.E / 1e9:.0f} GPa: strain {self.eps_design * 100:.3f} %.")
        print(f"Build orbit r0 = r / (1 + eps) = {self.r0:,.1f} m, altitude {(self.r0 - ring.R_EARTH_EQ) / 1e3:,.2f} km, circular velocity"
              f" {self.v0:,.1f} m/s, unstressed circumference C0 = {self.C0 / 1e3:,.1f} km; rise to the design point {(self.R_design - self.r0) / 1e3:,.2f} km.")
        print(f"{self.lim_sites:,} LIM sites: spacing {self.C0 / self.lim_sites:,.1f} m at the build, {self.L_design / self.lim_sites:,.1f} m at the design point;"
              f" {self.n_stations} anchor stations (one per {ANCHOR_STATION_SPACING / 1e3:.0f} km of ring).")
        print(f"Total masses: cable with hardware {self.M_r / 1e9:,.1f} Mt, casing {self.M_l / 1e9:,.1f} Mt; angular momentum at the build"
              f" {self.L_c + self.L_l:.5e} kg m^2/s, at the design point {self.M_r * self.R_design * self.v_cable_design - self.M_l * self.R_design * ring.V_GROUND_SYNC:.5e}"
              f" ({(self.M_r * self.R_design * self.v_cable_design - self.M_l * self.R_design * ring.V_GROUND_SYNC) / (self.L_c + self.L_l) * 100 - 100:+.3f} %, from the Earth in phase 2).")
        print("=" * 70 + "\n")

    def report(self, E_total_loop, t_final, series, csv_stride_rows=4000):
        """The rising-radius report at the end of the run, and the files rising_radius_*.csv / .json in the output folder.
        series: (time_days, v_cable, v_casing, tension, casing_weight, cable_excess, stress) lists from lim_simulation, with
        the final state already appended (this object's sampled lists are in step with them)."""
        import numpy as np
        out = ring.out_dir()
        t_days, v_cab, v_cas, T_s, W_s, F_s, S_s = series
        T_end, W_end, F_end, S_end = self.forces()
        KE_c, KE_l, PE, U = self.energy()
        KE_c0 = 0.5 * self.M_r * self.v0 ** 2; KE_l0 = 0.5 * self.M_l * self.v0 ** 2; PE0 = -self.GM * (self.M_r + self.M_l) / self.r0
        dE = sum(self.energy()) - self.E0
        rows = self.rows
        if rows[-1]["day"] < t_final / DAY - 1e-9:
            self.last_day = -1; self._log_day(t_final, self.F_last, 0.0)
        h = np.array([r["altitude_km"] for r in rows]); dd = np.array([r["day"] for r in rows])
        print("\n" + "=" * 70)
        print("RISING RADIUS: the deployment with the ring's radius as a state variable")
        print("=" * 70)
        print(f"Build orbit: r0 {self.r0:,.1f} m, altitude {(self.r0 - ring.R_EARTH_EQ) / 1e3:,.2f} km, circular velocity {self.v0:,.1f} m/s,"
              f" C0 {self.C0 / 1e3:,.1f} km, site spacing {self.C0 / self.lim_sites:,.1f} m ({self.lim_sites:,} sites).")
        print(f"End of the run at {t_final / DAY:,.1f} days: altitude {self.altitude_km:,.2f} km (r {self.r:,.1f} m), tension {T_end / 1e9:,.2f} GN,"
              f" stress {S_end / 1e9:.3f} GPa, strain {self.eps * 100:.4f} %, cable {self.v_c:,.2f} m/s, casing {-self.v_l:,.2f} m/s eastward,"
              f" relative {self.v_c - self.v_l:,.2f} m/s.")
        print(f"Phase 2 ended on {self.end_trigger}: cable speed {self.v_c - self.v_cable_design:+.4f} m/s from the design {self.v_cable_design:,.2f} m/s,"
              f" tension {(T_end - self.T_design) / 1e9:+.5f} GN from the design {self.T_design / 1e9:,.2f} GN"
              + (f"; the design speed at {self.t_v_design / DAY:,.2f} d" if self.t_v_design else "; the design speed not reached")
              + (f", the design tension at {self.t_T_design / DAY:,.2f} d." if self.t_T_design else ", the design tension not reached."))
        print(f"Newton solve of the radius: at most {self.newton_iters_max} iterations in a step.")
        print()
        print("Marks (interpolated on the sampled series; casing slowdown x = v0 - v_L in the launch frame):")
        print(f"{'x':>10} {'casing v_L':>12} {'cable v_c':>11} {'tension':>10} {'stress':>9} {'altitude':>10} {'weight':>11} {'excess':>11} {'time':>9}")
        xs = np.array([self.v0 - v for v in v_cas]); alt = np.array(self.s_alt)
        marks = [0.0, 1000.0, 3000.0, 5000.0, self.v0]
        if self.t_phase2 is not None:
            marks.append(self.v0 + self.omega * self.r_phase2)
        marks.append(xs[-1])
        for xm in marks:
            if xm > xs[-1] + 0.5:
                continue
            k = min(len(xs) - 1, int(np.searchsorted(xs, min(xm, xs[-1]))))      # xs is monotone in phase 1 and almost flat in phase 2
            vL = v_cas[k]; direction = "W" if vL > 0.5 else ("E" if vL < -0.5 else "")
            print(f"{xs[k]:8,.0f} m/s {abs(vL):8,.0f} m/s {direction:1} {v_cab[k]:8,.0f} m/s {T_s[k] / 1e9:7.1f} GN {S_s[k] / 1e9:6.2f} GPa {alt[k]:7.1f} km"
                  f" {W_s[k] / 1e3:7.1f} kN/m {F_s[k] / 1e3:7.1f} kN/m {t_days[k]:7.1f} d")
        print()
        if self.t_phase2 is not None:
            print(f"Phase 1 (no external torque) ends at {self.t_phase2 / DAY:,.1f} days: the casing at the ground-synchronous speed"
                  f" {self.omega * self.r_phase2:,.1f} m/s eastward, altitude {(self.r_phase2 - ring.R_EARTH_EQ) / 1e3:,.2f} km, cable {self.v_c_phase2:,.1f} m/s,"
                  f" tension {self.T_phase2 / 1e9:,.1f} GN, stress {self.T_phase2 / self.A / 1e9:.2f} GPa.")
            dur = (self.t_end if self.t_end else t_final) - self.t_phase2
            F_mean = self.F_anchor_int / max(self.t_phase2_dur, 1e-9)
            L_c_end = self.L_c
            print(f"Phase 2 (the casing held by the anchor lines, the LIMs pushing the cable) lasts {dur / DAY:,.1f} days, to {(self.t_end or t_final) / DAY:,.1f} days;"
                  f" the ring rises {(self.r - self.r_phase2) / 1e3:,.2f} km, the cable gains {self.v_c - self.v_c_phase2:,.1f} m/s, the tension {(T_end - self.T_phase2) / 1e9:,.1f} GN.")
            print(f"Angular momentum from the Earth {self.L_earth:.4e} kg m^2/s = {self.L_earth / L_c_end * 100:.3f} % of the cable's final {L_c_end:.4e}"
                  f" ({self.L_earth / (self.L_c_phase2 + self.L_l - 0) * 100:.3f} % of the ring's), {(self.L_c - self.L_c_phase2):.4e} to the cable in phase 2.")
            print(f"Anchor reaction (the Earth's torque over r): peak {self.F_anchor_peak[0] / 1e6:,.2f} MN total at {self.F_anchor_peak[1] / DAY:,.1f} d,"
                  f" {self.F_anchor_peak[0] / 1e6 / self.n_stations:.4f} MN per station; mean {F_mean / 1e6:,.2f} MN total, {F_mean / 1e6 / self.n_stations:.4f} MN per station"
                  f" ({self.n_stations} stations); the LIM thrust at the end {self.s_anchor[-1] / 1e6:,.2f} MN total.")
        else:
            print("The casing did not reach the ground-synchronous speed: no phase 2.")
        print()
        print("Energy budget (EJ):")
        print(f"  LIM work, thrust x relative speed (loop)        {E_total_loop / 1e18:12.6f}")
        print(f"  LIM work, the same sum in this module           {self.W_lim / 1e18:12.6f}")
        print(f"  cable kinetic energy change                     {(KE_c - KE_c0) / 1e18:12.6f}   ({KE_c0 / 1e18:.4f} to {KE_c / 1e18:.4f})")
        print(f"  casing kinetic energy change                    {(KE_l - KE_l0) / 1e18:12.6f}   ({KE_l0 / 1e18:.4f} to {KE_l / 1e18:.4f})")
        print(f"  potential energy change, -GM M / r              {(PE - PE0) / 1e18:12.6f}   (the ring lifted {(self.r - self.r0) / 1e3:,.2f} km)")
        print(f"  elastic energy, 1/2 E A eps^2 C0                {U / 1e18:12.6f}")
        print(f"  total mechanical energy change                  {dE / 1e18:12.6f}")
        print(f"  the Earth's work through the anchors (phase 2)  {self.W_earth / 1e18:12.6f}   (negative: the ring works on the Earth)")
        print(f"  LIM work required, delta E - W_earth            {(dE - self.W_earth) / 1e18:12.6f}")
        print(f"  difference, LIM work (loop) - required          {(E_total_loop - (dE - self.W_earth)) / 1e18:12.8f}")
        print(f"  fixed-radius deployment energy, ring_altitude   {ring.E_DEPLOY / 1e18:12.6f}")
        print()
        # ---- the daily table (every 25 days in the log; every day in the CSV)
        print("Daily (every 25th day in the log; all days in rising_radius_daily.csv):")
        print(f"{'day':>6} {'alt km':>8} {'T GN':>8} {'GPa':>6} {'strain%':>8} {'v_c':>8} {'v_L':>8} {'v_rel':>8} {'ph':>2} {'thrust MN':>9} {'anchor MN':>9}")
        for r_ in rows:
            if abs(r_["day"] - round(r_["day"])) < 1e-6 and round(r_["day"]) % 25 == 0 or r_ is rows[-1]:
                print(f"{r_['day']:6.1f} {r_['altitude_km']:8.2f} {r_['tension_GN']:8.2f} {r_['stress_GPa']:6.3f} {r_['strain'] * 100:8.4f} {r_['v_cable_m_s']:8.1f}"
                      f" {r_['v_casing_m_s']:8.1f} {r_['v_rel_m_s']:8.1f} {r_['phase']:2d} {r_['lim_thrust_total_MN']:9.2f} {r_['anchor_force_total_MN']:9.2f}")
        # ---- files
        daily = os.path.join(out, "rising_radius_daily.csv")
        with open(daily, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            for r_ in rows:
                wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) and k not in ("r_m", "L_cable", "L_casing", "L_earth") else
                                 (f"{v:.10e}" if k in ("L_cable", "L_casing", "L_earth") else (f"{v:.3f}" if k == "r_m" else v))) for k, v in r_.items()})
        ser = os.path.join(out, "rising_radius_deployment_tension.csv")
        n = len(t_days)
        with open(ser, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["time_days", "v_cable_m_s", "v_casing_m_s", "tension_GN", "casing_net_weight_kN_m", "cable_excess_kN_m", "stress_GPa",
                         "altitude_km", "radius_m", "strain", "phase", "anchor_force_total_MN"])
            stride = max(1, n // csv_stride_rows)
            idx = list(range(0, n, stride)) + ([n - 1] if (n - 1) % stride else [])
            for k in idx:
                wr.writerow([f"{t_days[k]:.5f}", f"{v_cab[k]:.3f}", f"{v_cas[k]:.3f}", f"{T_s[k] / 1e9:.6f}", f"{W_s[k] / 1e3:.5f}", f"{F_s[k] / 1e3:.5f}",
                             f"{S_s[k] / 1e9:.6f}", f"{self.s_alt[k]:.4f}", f"{self.s_r[k]:.2f}", f"{self.s_strain[k]:.7f}", self.s_phase[k], f"{self.s_anchor[k] / 1e6:.4f}"])
        print(f"Time series ({n} points; every {stride}th written with the initial and final states) in {os.path.basename(ser)}; the daily table in {os.path.basename(daily)}.")
        summary = dict(
            r0_m=self.r0, build_altitude_km=(self.r0 - ring.R_EARTH_EQ) / 1e3, v0_m_s=self.v0, C0_km=self.C0 / 1e3,
            site_spacing_build_m=self.C0 / self.lim_sites, site_spacing_design_m=self.L_design / self.lim_sites, lim_sites=self.lim_sites,
            rise_km=(self.R_design - self.r0) / 1e3, strain_design=self.eps_design, E_Pa=self.E, A_m2=self.A,
            t_end_days=t_final / DAY, end_trigger=self.end_trigger,
            t_v_design_days=(self.t_v_design / DAY if self.t_v_design else None), t_T_design_days=(self.t_T_design / DAY if self.t_T_design else None),
            t_phase2_days=(self.t_phase2 / DAY if self.t_phase2 else None), altitude_phase2_km=((self.r_phase2 - ring.R_EARTH_EQ) / 1e3 if self.r_phase2 else None),
            v_cable_phase2_m_s=self.v_c_phase2, tension_phase2_GN=(self.T_phase2 / 1e9 if self.T_phase2 else None),
            phase2_duration_days=(((self.t_end or t_final) - self.t_phase2) / DAY if self.t_phase2 else None),
            final_altitude_km=self.altitude_km, final_tension_GN=T_end / 1e9, final_stress_GPa=S_end / 1e9, final_strain=self.eps,
            final_v_cable_m_s=self.v_c, final_v_casing_m_s=self.v_l, final_v_rel_m_s=self.v_c - self.v_l,
            final_casing_weight_kN_m=W_end / 1e3, final_cable_excess_kN_m=F_end / 1e3,
            L_earth=self.L_earth, L_cable_final=self.L_c, L_earth_share_of_cable=self.L_earth / self.L_c,
            anchor_peak_MN=self.F_anchor_peak[0] / 1e6, anchor_peak_day=self.F_anchor_peak[1] / DAY,
            anchor_mean_MN=self.F_anchor_int / max(self.t_phase2_dur, 1e-9) / 1e6, n_stations=self.n_stations,
            E_lim_loop_EJ=E_total_loop / 1e18, E_lim_module_EJ=self.W_lim / 1e18, dKE_cable_EJ=(KE_c - KE_c0) / 1e18, dKE_casing_EJ=(KE_l - KE_l0) / 1e18,
            dPE_EJ=(PE - PE0) / 1e18, elastic_EJ=U / 1e18, dE_mech_EJ=dE / 1e18, W_earth_EJ=self.W_earth / 1e18,
            balance_EJ=(E_total_loop - (dE - self.W_earth)) / 1e18, E_deploy_fixed_EJ=ring.E_DEPLOY / 1e18, newton_iters_max=self.newton_iters_max)
        with open(os.path.join(out, "rising_radius_summary.json"), "w") as fh:
            json.dump(summary, fh, indent=1)
        print("Summary in rising_radius_summary.json.")
        print("=" * 70)
        return summary
