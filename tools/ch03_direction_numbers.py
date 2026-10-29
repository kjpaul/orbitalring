#!/usr/bin/env python3
"""
The Chapter 3 numbers for one cable direction: masses, speeds, tension, the static wave speed and the Mach number of the cable past the casing, deployment energy, LIM frequency, the launch speed window.
Book: Orbital Ring Engineering, Volume I: Chapter 3 (Table 3.1).
Run: python tools/ch03_direction_numbers.py, then ORBITAL_RING_CABLE=prograde python tools/ch03_direction_numbers.py   (one column each)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

Chapter 3 numbers: retrograde against prograde cable at 700 km.  Run once with
ORBITAL_RING_CABLE unset (retrograde) and once with ORBITAL_RING_CABLE=prograde."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import math, os, sys
import ring_altitude as R
m_c = R.M_CABLE_M; m_l = 12000.0
print("cable", os.environ.get("ORBITAL_RING_CABLE", "retrograde"))
print(" structural %.0f  with hw %.0f  ring %.0f kg/m ; structural:load %.2f" % (R.M_CABLE_STRUCTURAL, m_c, R.M_RING_M, R.M_CABLE_STRUCTURAL/m_l))
print(" v_orbit %.1f v_cable %.1f v_geo %.1f dv_load %.1f dv_cable %.1f" % (R.V_ORBIT, R.V_CABLE, R.V_GROUND_SYNC, R.DELTA_V_LOAD, R.V_CABLE-R.V_ORBIT))
vrel = getattr(R, "V_REL", None)
sign = -1 if abs(R.DELTA_V_LOAD-(R.V_ORBIT-R.V_GROUND_SYNC)) < 1 else 1
vrel = R.V_CABLE + sign*R.V_GROUND_SYNC
c = math.sqrt(R.CABLE_TENSION/m_c)
print(" v_rel %.1f  tension %.1f GN  area %.2f m2  side %.2f m" % (vrel, R.CABLE_TENSION/1e9, R.CABLE_AREA, R.CABLE_SIDE))
print(" static wave speed sqrt(T/m) %.0f m/s ; Mach of cable past casing %.2f" % (c, vrel/c))
v_cas0 = R.V_ORBIT; v_cas1 = R.V_GROUND_SYNC
e = 0.5*m_c*(R.V_CABLE**2-R.V_ORBIT**2) + 0.5*m_l*(v_cas1**2-v_cas0**2)
print(" deployment energy %.3e J/m  -> %.2f EJ" % (e, e*R.L_RING/1e18))
print(" LIM frequency v_rel/200 = %.1f Hz" % (vrel/200))
print(" westward launch on this cable: sled matches the cable at %.0f m/s past the casing; subsonic between %.0f and %.0f m/s" % (vrel, vrel-c, vrel+c))
print(" g %.4f centrifugal %.4f g_net %.4f ; bearing load %.1f against %.1f kN/m" % (R.G_LOCAL, R.V_GROUND_SYNC**2/R.R_ORBIT, R.G_NET, 12*R.G_NET, 12*R.G_LOCAL))
