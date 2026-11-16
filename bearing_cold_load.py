#!/usr/bin/env python3
"""
The bearing's cold load through the deployment, as the design makes it (round 9, 2026-11-16, correction 3 of
Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md).  It replaces the flat lim_config.Q_COLD_EXTRA_PER_SITE of rounds 5 to 8.
Book: Orbital Ring Engineering, Volume III: Chapters 7, 9 and 14.

Per LIM site (500 m of ring), four parts:
  1. the fast coils' leads' conduction with the connected families (constant);
  2. the regulation's Joule heat in those leads with the support following the casing's weight continuously (constant);
  3. the support coils' own leads while the supply carries the rise of the current: the conduction-cooled lead's heat at the coil current
     of that moment (the Wiedemann-Franz lead of bearing_electrical.lead_heat, the section fixed);
  4. the flux pumps' loss after the leads come off, while the pumps run.
The coil current follows the casing's weight on the bearing through the solved pull curve of a support track (each track carries half).
The leads are attached from day 0.  They come off when the rate the current must change at, taken over the last WINDOW days of the run's own
series, has fallen to what the pump's voltage gives into the coil's inductance (V / L), and they go back on if the rate ever exceeds it again.
Nothing here is put in by hand: the day comes out of the run.

Use: lim_config.BEARING_COLD_LOAD = <path to a JSON table> (tools/round9_cold_load.py writes it).  lim_simulation.py calls update() each step
and write() at the end (bearing_cold_load.csv and bearing_cold_load.json in the run's folder).
"""
import json
import math
import os

L0 = 2.44e-8             # W ohm / K2, the Lorenz number
DAY = 86400.0
WINDOW = 1.0             # days over which the current's rate is taken


def lead_heat(I, aol, T_h, T_c, k_cu):
    """The heat one conduction-cooled copper lead of section over length aol (m) brings to its cold end at current I (bearing_electrical.lead_heat)."""
    w = I * math.sqrt(L0) / (k_cu * aol)
    if w < 1e-9:
        return k_cu * aol * (T_h - T_c)
    if w >= math.pi or abs(math.sin(w)) < 1e-9:
        return float("inf")
    b = (T_c - T_h * math.cos(w)) / math.sin(w)
    return -k_cu * aol * (-T_h * w * math.sin(w) + b * w * math.cos(w))


class BearingColdLoad:
    def __init__(self, path):
        t = json.load(open(path, encoding="utf-8"))
        self.t = t; self.path = path
        self.F = t["F_track_N_m"]; self.I = t["I_coil_A"]
        self.rate_pump = t["pump_V"] / t["L_coil_H"] * DAY            # A per day the pump can move the coil's current
        self.attached = True
        self.hist = [(0.0, 0.0)]                                      # (day, coil current) samples, at least 0.05 day apart
        self.rate = float("inf")
        self.events = [(0.0, "leads attached")]
        self.rows = []; self._last_row = -1.0
        self.q = self.parts(0.0)[0]

    def current(self, w_casing):
        """The support coil current (A) with the casing's weight on the bearing w_casing (N/m); each of the two tracks carries half."""
        f = max(w_casing, 0.0) / 2.0
        F, I = self.F, self.I
        if f <= F[0]:
            return I[0]
        if f >= F[-1]:
            return I[-1]
        lo, hi = 0, len(F) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if F[m] <= f: lo = m
            else: hi = m
        return I[lo] + (I[hi] - I[lo]) * (f - F[lo]) / (F[hi] - F[lo])

    def parts(self, I):
        t = self.t
        q_lead = t["n_leads_site"] * lead_heat(I, t["lead_aol_m"], t["T_h"], t["T_c"], t["k_cu"]) if self.attached else 0.0
        q_pump = 0.0 if self.attached else t["n_pumps_site"] * t["pump_q_W"]
        return t["cond_W"] + t["reg_W"] + q_lead + q_pump, t["cond_W"], t["reg_W"], q_lead, q_pump

    def update(self, t_s, w_casing):
        """The extra cold load per site (W) at time t_s (s) with the casing's weight on the bearing w_casing (N/m)."""
        d = t_s / DAY; I = self.current(w_casing)
        if d - self.hist[-1][0] >= 0.05:
            self.hist.append((d, I))
            while len(self.hist) > 2 and d - self.hist[1][0] >= WINDOW:
                self.hist.pop(0)
            d0, I0 = self.hist[0]
            if d - d0 >= 0.999 * WINDOW:
                self.rate = (I - I0) / (d - d0)
                if self.attached and abs(self.rate) <= self.rate_pump:
                    self.attached = False; self.events.append((d, "leads off, the pumps follow"))
                elif not self.attached and abs(self.rate) > self.rate_pump:
                    self.attached = True; self.events.append((d, "leads back on"))
        p = self.parts(I); self.q = p[0]
        if d - self._last_row >= 1.0 - 1e-9 or not self.rows:
            self._last_row = math.floor(d + 1e-9)
            self.rows.append((d, w_casing, I, self.rate if math.isfinite(self.rate) else float("nan"), int(self.attached)) + p)
        return self.q

    def write(self, folder):
        with open(os.path.join(folder, "bearing_cold_load.csv"), "w", encoding="utf-8") as f:
            f.write("day,casing_weight_N_m,coil_current_A,rate_A_per_day,leads_attached,extra_W,fast_conduction_W,regulation_W,support_leads_W,pumps_W\n")
            for r in self.rows:
                f.write(",".join(f"{x:.6g}" for x in r) + "\n")
        off = [e[0] for e in self.events if e[1].startswith("leads off")]
        att = [r for r in self.rows if r[4]]
        json.dump(dict(table=self.path, pump_rate_A_per_day=self.rate_pump, events=self.events, leads_off_day=off[-1] if off and not self.attached else None,
                       leads_on_days=len(att), lead_heat_mean_W=sum(r[8] for r in att) / max(len(att), 1), lead_heat_max_W=max((r[8] for r in att), default=0.0),
                       extra_mean_W=sum(r[5] for r in self.rows) / max(len(self.rows), 1), extra_max_W=max(r[5] for r in self.rows), extra_end_W=self.rows[-1][5],
                       I_max_A=max(r[2] for r in self.rows), I_end_A=self.rows[-1][2]),
                  open(os.path.join(folder, "bearing_cold_load.json"), "w", encoding="utf-8"), indent=1)
