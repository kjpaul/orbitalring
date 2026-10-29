"""
Runs the LIM deployment simulation at its defaults (8 MW per site, the cable warm, the winding and loss model of lim_config.py)
and writes the thermal time series Figures 14.1 and 14.2 are drawn from: cable and wall temperature, the cold-side load and the
radiator width, the compressor power, the hysteresis and eddy loss of the two LIMs of a site, the warm-loop load, the coil current
and the slip frequency, one row per SKIP steps.
Book: Orbital Ring Engineering, Volume I: Chapter 14 (Figures 14.1 and 14.2, Table 14.1).
Run: python tools/ch14_thermal_series.py <out.csv>     (about ten minutes, the same run as lim_simulation.py)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT
Copied from Vol_III/renumber_scripts_2026-10-12/ch18_thermal_series.py (written under the old numbering) on 2026-10-27.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import sys, os, csv
import lim_config as cfg
import lim_simulation as sim

out = sys.argv[1]
v_slip = cfg.V_SLIP_MAX
i_peak = cfg.I_MIN
param_str, success, total_time = sim.run_deployment_simulation(v_slip, i_peak, cfg.THRUST_MODEL, cfg.EDDY_HEAT_TO_CABLE, False, 0.0)
print("success", success, "total_time_days", total_time/cfg.DAY, "samples", len(sim.data_cable_temp), "skip_s", cfg.SKIP)
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s","cable_temp_K","wall_temp_K","q_cryo_cold_W","radiator_width_m","p_cryo_W","p_hyst_site_W","p_eddy_site_W","warm_W_per_m","current_A","f_slip_Hz"])   # the site columns are the two LIMs together
    for i in range(len(sim.data_cable_temp)):
        w.writerow([i*cfg.SKIP, sim.data_cable_temp[i], sim.data_wall_temp[i], sim.data_q_cryo_cold[i], sim.data_radiator_width[i], sim.data_p_cryo[i], sim.data_p_hyst[i], sim.data_p_eddy[i], sim.data_warm_rad_width[i], sim.data_current[i], sim.data_f_slip[i]])
print("wrote", out)
