#!/usr/bin/env python3
"""
The Chapter 8 numbers: ring-mounted panel width for 8 and 16 MW per site, noon output, power flow and conductor mass around the ring, the expanded 400 m case, section buses; from power_config.py and power_physics.py.
Book: Orbital Ring Engineering, Volume I: Chapter 8 (Table 8.1 and the text).
Run: ORBITAL_RING_ARRAY_OFFSET_KM=0 python tools/ch08_numbers.py   (the panels on the ring, as the chapter assumes)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_ARRAY_OFFSET_KM (must be 0 for the chapter's numbers)

Chapter 8 numbers for the 700 km ring (ring-mounted panels, array offset 0).
Uses power_config.py and power_physics.py unchanged."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, os, sys
import power_config as cfg, power_physics as phys

print("Ring circumference km", cfg.L_RING/1e3, " sites", cfg.LIM_SITES)
print("Array altitude km", cfg.ALTITUDE/1e3)
half = math.degrees(math.asin(cfg.R_EARTH/cfg.R_ORBIT))
print("Earth angular radius seen from ring deg %.2f ; sunlit half-angle %.2f ; sunlit fraction %.4f" % (half, 180-half, (180-half)/180))
print("Albedo peak W/m2 %.1f" % cfg.ALBEDO_PEAK)
prof = phys.compute_ring_power_profile(1.0)
avg, peak = prof['p_avg_per_m2'], prof['p_peak_per_m2']
print("avg %.2f peak %.2f W/m2  peak/avg %.3f" % (avg, peak, peak/avg))
print("output at 90 deg (dawn/dusk) %.1f ; at midnight %.2f W/m2" % (phys.panel_electrical_output(math.pi/2), phys.panel_electrical_output(math.pi)))
n = cfg.N_POINTS; dphi = 2*math.pi/n

def case(label, site_power):
    demand_total = cfg.LIM_SITES*site_power
    width = demand_total/(avg*cfg.L_ARRAY)
    area_site = site_power/avg
    p = phys.compute_ring_power_profile(width)
    dem_m = demand_total/cfg.L_RING
    flow = phys.compute_power_flow(p['p_gen_per_m'], dem_m, dphi)
    # crossover angle: where generation per m = demand per m
    cross = None
    for i in range(n//2, n-1):
        if p['p_gen_per_m'][i] >= dem_m > p['p_gen_per_m'][i+1]:
            cross = math.degrees(p['phi'][i]+dphi/2); break
    c = phys.size_conductor(flow['p_flow_peak'])
    print("\n== %s: %.2f MW per site ==" % (label, site_power/1e6))
    print(" total %.2f GW ; matched width %.2f m (%.2f per side) ; panel area per site %.0f m2 ; per m of ring %.1f kW" % (demand_total/1e9, width, width/2, area_site, width*avg/1e3))
    print(" noon output per site %.2f MW ; surplus/deficit crossover at +/- %.1f deg" % (peak*width*cfg.LIM_SPACING/1e6, cross))
    print(" peak flow %.2f GW ; current %.0f A ; conductor %.1f cm2 total (%.4f m2) ; %d cables of %.0f cm ; %.1f kg/m ; %.2f Mt ; peak loss %.2f GW" % (
        flow['p_flow_peak']/1e9, c['i_peak'], c['a_total']*1e4, c['a_total'], c['n_cables_total'], cfg.CABLE_DIAMETER*100, c['m_per_m'], c['m_total']/1e9, c['p_loss_peak']/1e9))
    print(" share of 12,000 kg/m casing: conductor %.2f %%" % (c['m_per_m']/12000*100))
    ops = cfg.P_DEMAND_OPS
    print(" after deployment: generation %.2f - operations %.2f - loss %.2f = export %.2f GW (even export around the ring gives the same flow as above)" % (
        demand_total/1e9, ops/1e9, c['p_loss_peak']/1e9, (demand_total-ops-c['p_loss_peak'])/1e9))
    print(" per 10 km tap section (20 sites): %.1f MW ; at 100 kV %.2f kA" % (20*site_power/1e6, 20*site_power/100e3/1e3))
    return width

case("8 MW", 8e6)
case("16 MW", 16e6)
# expanded 400 m
w = 400.0
ps = avg*w*cfg.LIM_SPACING
print("\nExpanded 400 m: %.1f kW per m ; %.2f MW per site ; total %.3f TW" % (avg*w/1e3, ps/1e6, ps*cfg.LIM_SITES/1e12))
case("400 m expanded", ps)
print("\nPanel mass at 0.5 kg/m2: 110 m -> %.0f kg/m ; 400 m -> %.0f kg/m (%.2f %% of 12,000)" % (0.5*110, 0.5*400, 0.5*400/120))
print("Operations demand %.2f GW ; ratio 8 MW / 0.5 MW = %.0f" % (cfg.P_DEMAND_OPS/1e9, 8/0.5))
taps = cfg.L_RING/10e3
print("Tap points at 10 km: %.0f ; anchor stations at 50 km: %.0f" % (taps, cfg.L_RING/50e3))
print("Specific conductivity CNT %.0f S m2/kg ; copper %.0f ; ratio %.2f" % (cfg.SPECIFIC_CONDUCTIVITY, 5.96e7/8960, cfg.SPECIFIC_CONDUCTIVITY/(5.96e7/8960)))
# section bus: resistive versus HTS, expanded case
for label, P in (("16 MW", 20*16e6), ("400 m", 20*ps)):
    I = P/100e3
    # resistive: feed from tap at section centre, longest run 5 km, out and back = 10 km of conductor, 1 % loss at full current over full run
    R = 0.01*100e3/I
    A = 10e3/(cfg.SIGMA_CONDUCTOR*R)
    m_res = 2*A*cfg.RHO_CONDUCTOR   # two conductors of area A each? A is per conductor (10 km loop total), mass per m of ring = 2 conductors
    tapes = math.ceil(I/650.0)
    m_tape = 12e-3*0.1e-3*8900      # 12 mm x 0.1 mm, mostly Hastelloy/copper, 8,900 kg/m3
    print("Section bus %s: %.0f MW, %.2f kA ; resistive at 1%% loss: %.1f cm2 each way, %.1f kg/m ; HTS: %d tapes each way at 650 A, %.2f kg/m ; ratio %.0f" % (
        label, P/1e6, I/1e3, A*1e4, m_res, tapes, 2*tapes*m_tape, m_res/(2*tapes*m_tape)))
print("Export per anchor station (889): 8 MW case %.0f MW ; 16 MW case %.0f MW" % (658.52e3/889, 1361.51e3/889))
print("Sled 10 km long at 30 km/s passes a point in %.2f s" % (10/30.0))
print("LIM supply frequency: v_rel/(2 tau) = 8786.2/200 = %.1f Hz (lim_simulation.txt gives f_supply_max 44.15 Hz with slip)" % (8786.2/200))
