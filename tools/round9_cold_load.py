#!/usr/bin/env python3
"""
The table bearing_cold_load.py reads (round 9, 2026-11-16, correction 3): the bearing's cold load per site through the deployment.
Book: Orbital Ring Engineering, Volume III: Chapters 7, 9 and 14.
Run (from Python/orbitalring): python tools/round9_cold_load.py [--pump=mallett|dynamo]     writes runs/round9_logs/bearing_cold_load_<pump>.json
Sources: the support track's solved pull curve and the coil's inductance (rounds 7 and 8, bearing_gap_trade_parts.support_current), the fast leads'
conduction with ten families (round 8's design point), the regulation's Joule heat with the support following the weight continuously and the
support's leads (the gap trade, item 4), the pumps read for Stage A (runs/bearing_layout_pumps.txt).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json
import os
import sys

import numpy as np

import ring_altitude as ring
import bearing_electrical as be
import bearing_round7 as r7
import bearing_round8 as r8
import bearing_gap_trade_parts as gp

PUMPS = {"mallett": dict(V=0.267e-3, q=10.0, src="Mallett et al. 2023, self-rectifier, 320 A, 0.267 mV time-averaged; about 1 J per cycle in the bridge at 10 Hz (read for Stage A, approximate)"),
         "dynamo": dict(V=0.197e-3, q=0.15, src="Bumby et al. 2016, dynamo, 30 A, 0.197 mV time-averaged; 135 mW of ac loss in the stator tape (Ghabeli et al. 2021, a model)")}


def table(pump="mallett", F_max=None):
    I_of_F, c = gp.support_current()
    hm = r7.load("round7_support.json")["half"]["0.050"]; L2 = hm["L_loop"] + hm["M_loop"]
    pj = json.load(open(os.path.join(r7.OUT, "gap_trade_logs", "gap_trade_pump.json"), encoding="utf-8"))
    I_cap = r8.load("round8_rows.json")["parked"]["I_cap"]
    aol = be.lead_heat(I_cap)[1]
    F = np.concatenate([np.linspace(0.0, 2000.0, 81)[:-1], np.linspace(2000.0, F_max or 1.3 * c["F_n"], 400)])
    P = PUMPS[pump]
    return dict(pump=pump, pump_source=P["src"], pump_V=P["V"], pump_q_W=P["q"], L_coil_H=L2, n_pumps_site=pj["pumps_per_site"], n_leads_site=pj["sup_leads"]["n"],
                lead_aol_m=aol, T_h=be.T_WARM, T_c=be.T_COLD, k_cu=400.0, cond_W=pj["cond"], reg_W=pj["j_reg"],
                F_track_N_m=[float(x) for x in F], I_coil_A=[float(I_of_F(x)) for x in F], I_normal_A=c["I_n"], F_normal_N_m=c["F_n"])


def main():
    pump = next((a[7:] for a in sys.argv if a.startswith("--pump=")), "mallett")
    t = table(pump)
    out = os.path.join(r7.OUT, "round9_logs"); os.makedirs(out, exist_ok=True)
    p = os.path.join(out, f"bearing_cold_load_{pump}.json")
    json.dump(t, open(p, "w", encoding="utf-8"), indent=1)
    import bearing_cold_load as b
    o = b.BearingColdLoad(p)
    print(f"wrote {p}: pump {t['pump_V']*1e3:.3f} mV into {t['L_coil_H']:.2f} H = {o.rate_pump:.3f} A a day; conduction {t['cond_W']:,.0f} W, regulation {t['reg_W']:,.0f} W,")
    print(f"  support leads {t['n_leads_site']:.0f} per site, A/l {t['lead_aol_m']:.3e} m: {b.lead_heat(0, t['lead_aol_m'], t['T_h'], t['T_c'], 400.0)*t['n_leads_site']:,.0f} W idle, "
          f"{b.lead_heat(t['I_normal_A'], t['lead_aol_m'], t['T_h'], t['T_c'], 400.0)*t['n_leads_site']:,.0f} W at {t['I_normal_A']:.0f} A; pumps {t['n_pumps_site']*t['pump_q_W']:,.0f} W per site while they run")


if __name__ == "__main__":
    main()
