#!/usr/bin/env python3
"""
Single source for the ring's design altitude and everything that follows
from it.  Every simulation in this folder imports its altitude-dependent
constants from here, so the design point is changed in one place.

Design case (2026-10-04): 700 km, retrograde cable, suspended array at 600 km.
(2026-10-02 to 2026-10-04 the design case was 800 km with the array 50 km below.)
Override for comparison runs without editing any file:

    ORBITAL_RING_ALT_KM=2000 python lim_simulation.py
    ORBITAL_RING_CABLE=prograde python j2_simulation.py

Run this file to print the reference table for the current altitude.

Reference: "The Orbital Ring" by Paul G de Jong
"""
import math
import os

GM = 3.986004418e14            # Earth gravitational parameter (m^3/s^2)
R_EARTH_EQ = 6_378_137.0       # Earth equatorial radius (m)
OMEGA_SIDEREAL = 7.2921159e-5  # Earth sidereal rotation rate (rad/s)

ALTITUDE_KM = float(os.environ.get("ORBITAL_RING_ALT_KM", "700"))
RETROGRADE = os.environ.get("ORBITAL_RING_CABLE", "retrograde").lower() != "prograde"

M_LOAD_M = 12_000.0            # casing mass per meter (kg/m)
M_HARDWARE_M = 4_076.0         # cable hardware per meter (kg/m)
RHO_CNT = 1700.0               # kg/m^3
SIGMA_CABLE = 12.5e9           # Pa, operating stress (safety factor 2 on 25 GPa)

ALTITUDE = ALTITUDE_KM * 1000.0
R_ORBIT = R_EARTH_EQ + ALTITUDE
L_RING = 2 * math.pi * R_ORBIT
V_ORBIT = math.sqrt(GM / R_ORBIT)                  # circular velocity
V_ESCAPE = math.sqrt(2 * GM / R_ORBIT)
T_ORBIT = 2 * math.pi * R_ORBIT / V_ORBIT
G_LOCAL = GM / R_ORBIT**2                          # gravity at the ring
V_GROUND_SYNC = OMEGA_SIDEREAL * R_ORBIT           # casing speed, eastward, inertial
G_NET = G_LOCAL - V_GROUND_SYNC**2 / R_ORBIT       # net downward accel. on the casing

# Velocity change of the casing during deployment.
DELTA_V_LOAD = V_ORBIT + V_GROUND_SYNC if RETROGRADE else V_ORBIT - V_GROUND_SYNC

# Casing end-of-deployment velocity measured in the launch direction of the
# ring (westward positive for a retrograde ring).  The deployment loop runs
# the casing from V_ORBIT down to this value; for the retrograde ring it
# passes through zero and ends at -V_GROUND_SYNC.
V_CASING_FINAL_LAUNCH_FRAME = -V_GROUND_SYNC if RETROGRADE else V_GROUND_SYNC


def cable_mass_structural(m_load=M_LOAD_M, m_hw=M_HARDWARE_M, delta_v=DELTA_V_LOAD,
                          v_orbit=V_ORBIT, g_net=G_NET, r_orbit=R_ORBIT,
                          rho=RHO_CNT, sigma=SIGMA_CABLE):
    """Positive root of the cable sizing quadratic (same as lim_physics)."""
    S = sigma / rho
    K = m_load * delta_v * 2 * v_orbit - m_load * g_net * r_orbit
    P = (m_load * delta_v) ** 2
    a, b, c = S, S * m_hw - K, -(K * m_hw + P)
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


M_CABLE_STRUCTURAL = cable_mass_structural()
M_CABLE_M = M_CABLE_STRUCTURAL + M_HARDWARE_M      # cable with hardware
M_RING_M = M_CABLE_M + M_LOAD_M
# Cable speed (magnitude) from momentum conservation.
V_CABLE = V_ORBIT + M_LOAD_M * DELTA_V_LOAD / M_CABLE_M
V_REL = V_CABLE + V_GROUND_SYNC if RETROGRADE else V_CABLE - V_GROUND_SYNC
CABLE_AREA = M_CABLE_STRUCTURAL / RHO_CNT
CABLE_SIDE = math.sqrt(CABLE_AREA)
CABLE_TENSION = SIGMA_CABLE * CABLE_AREA
E_DEPLOY = 0.5 * L_RING * (M_CABLE_M * (V_CABLE**2 - V_ORBIT**2)
                           + M_LOAD_M * (V_GROUND_SYNC**2 - V_ORBIT**2))


def design(altitude_km=ALTITUDE_KM, retrograde=RETROGRADE, m_load=M_LOAD_M, m_hw=M_HARDWARE_M):
    """Ring design at any altitude and cable direction (dict), for comparison tables."""
    r = R_EARTH_EQ + altitude_km * 1000.0
    v_orb = math.sqrt(GM / r)
    v_g = OMEGA_SIDEREAL * r
    g = GM / r**2
    g_net = g - v_g**2 / r
    dv = v_orb + v_g if retrograde else v_orb - v_g
    m_s = cable_mass_structural(m_load, m_hw, dv, v_orb, g_net, r)
    m_c = m_s + m_hw
    v_c = v_orb + m_load * dv / m_c
    L = 2 * math.pi * r
    return dict(altitude_km=altitude_km, retrograde=retrograde, r=r, L_ring=L, v_orbit=v_orb,
                v_ground=v_g, g=g, g_net=g_net, delta_v_load=dv, m_struct=m_s, m_cable=m_c,
                m_ring=m_c + m_load, v_cable=v_c,
                v_rel=v_c + v_g if retrograde else v_c - v_g,
                area=m_s / RHO_CNT, side=math.sqrt(m_s / RHO_CNT), tension=SIGMA_CABLE * m_s / RHO_CNT,
                M_ring=(m_c + m_load) * L,
                E_deploy=0.5 * L * (m_c * (v_c**2 - v_orb**2) + m_load * (v_g**2 - v_orb**2)))


# Solar array: hangs below the ring, spread by the anchor lines (Paul, 2026-10-04: 100 km, array at 600 km).
ARRAY_OFFSET = float(os.environ.get("ORBITAL_RING_ARRAY_OFFSET_KM", "100")) * 1000.0  # 0 = panels mounted on the ring
ARRAY_ALTITUDE = ALTITUDE - ARRAY_OFFSET
R_ARRAY = R_EARTH_EQ + ARRAY_ALTITUDE
L_ARRAY = 2 * math.pi * R_ARRAY


def report():
    d = "retrograde" if RETROGRADE else "prograde"
    rows = [
        ("Altitude", ALTITUDE_KM, "km"), ("Cable direction", d, ""),
        ("Radius", R_ORBIT / 1e3, "km"), ("Circumference", L_RING / 1e3, "km"),
        ("g at ring", G_LOCAL, "m/s^2"), ("g net on casing", G_NET, "m/s^2"),
        ("Circular velocity", V_ORBIT, "m/s"), ("Escape velocity", V_ESCAPE, "m/s"),
        ("Orbital period", T_ORBIT / 60, "min"),
        ("Ground-sync velocity", V_GROUND_SYNC, "m/s"),
        ("Delta-v load", DELTA_V_LOAD, "m/s"),
        ("Cable structural", M_CABLE_STRUCTURAL, "kg/m"),
        ("Cable with hardware", M_CABLE_M, "kg/m"), ("Ring", M_RING_M, "kg/m"),
        ("Cable velocity", V_CABLE, "m/s"), ("Cable-casing relative", V_REL, "m/s"),
        ("Cable area", CABLE_AREA, "m^2"), ("Cable side", CABLE_SIDE, "m"),
        ("Cable tension", CABLE_TENSION / 1e9, "GN"),
        ("Cable centrifugal accel.", V_CABLE**2 / R_ORBIT, "m/s^2"),
        ("Cable net outward accel.", V_CABLE**2 / R_ORBIT - G_LOCAL, "m/s^2"),
        ("Casing weight", M_LOAD_M * G_NET / 1e3, "kN/m"),
        ("Cable mass ratio (with hw / casing)", M_CABLE_M / M_LOAD_M, ""),
        ("Cable mass, structural", M_CABLE_STRUCTURAL * L_RING / 1e9, "Mt"),
        ("Cable mass, with hardware", M_CABLE_M * L_RING / 1e9, "Mt"),
        ("Casing mass", M_LOAD_M * L_RING / 1e9, "Mt"),
        ("Ring mass", M_RING_M * L_RING / 1e9, "Mt"),
        ("Deployment energy", E_DEPLOY / 1e18, "EJ"),
        ("LIM sites at 500 m", round(L_RING / 500.0), ""),
        ("Array altitude", ARRAY_ALTITUDE / 1e3, "km"),
        ("Array circumference", L_ARRAY / 1e3, "km"),
    ]
    for name, val, unit in rows:
        if isinstance(val, str):
            print(f"  {name:38} {val:>14} {unit}")
        else:
            print(f"  {name:38} {val:>14,.3f} {unit}")


if __name__ == "__main__":
    report()
