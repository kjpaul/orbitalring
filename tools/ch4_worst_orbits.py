#!/usr/bin/env python3
"""
Table 4.7 of the book (the released cable's lifetimes on its worst-case orbits) and the two lifetime rows of Table 4.12, under
the design point of ring_altitude.py: the five orbits (the slow half with the full radial kick at 700, 600 and 550 km, and the
700 km slow half that has swept up a 1,200 or a 3,000 kg/m casing roof, with the kick), each for the full cross-section (the
ballistic coefficient of the cable's side), a 0.8 m piece (600 kg/m^2) and an 8 cm piece (60 kg/m^2), by the rules the
repository has: (A) debris_orbits.lifetime, the orbit-averaged decay on the Vallado mean atmosphere; (B) the same divided by
the factor the book's footnote quotes for the cycle-averaged atmosphere, 1.5 and 2, and by the factor the October book
implied against the October run of rule A (1.67 at 700 km, 1.44 at 600 km); (C) debris_orbits.lifetime_quick, the circular
orbit at the perigee (the pessimistic bound).  Paul rules on the rule.
Book: Orbital Ring Engineering, Volume I: Chapter 4, Tables 4.7 and 4.12, the floor sentence of 4.5.
Run: python tools/ch4_worst_orbits.py
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM, ORBITAL_RING_QUADRATIC (the ring, through ring_altitude.py)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import math
import debris_orbits as do
import ring_altitude as ring

V = do.V_RECOIL
side = math.sqrt(ring.M_CABLE_STRUCTURAL / do.RHO_CNT)
B_ROD = do.ballistic_coefficient(side)
BS = [("full section", B_ROD), ("0.8 m piece", 600.0), ("8 cm piece", 60.0)]
OCT_FACTOR = {700: 65950.48 / 39400.0, 600: 15086.13 / 10500.0}     # the October book over the October run of rule A


def orbit_kick(h_ring, m_swept=0.0):
    """the slow half (-V along the orbit), having swept up m_swept kg/m of the roof, with the full radial kick V"""
    s = do.ring_state(h_ring)
    v = do.swept_cable_speed(h_ring, m_swept, -V) if m_swept else abs(s['v_cable']) - V
    h_p, h_a, a, e = do.orbit_from_state(s['r'], v, V)
    return v, h_p, h_a


cases = [("700 km, slow half with the radial kick", 700e3, 0.0), ("700 km, with the 1,200 kg/m roof swept up", 700e3, 1200.0),
         ("700 km, with a 3,000 kg/m roof swept up", 700e3, 3000.0), ("600 km, slow half with the radial kick", 600e3, 0.0),
         ("550 km, slow half with the radial kick", 550e3, 0.0)]

print(f"Design point: cable {ring.M_CABLE_STRUCTURAL:,.0f} kg/m structural, side {side:.2f} m, full-section B = {B_ROD:,.0f} kg/m^2 (C_D 2.2); quadratic {ring.QUADRATIC}.")
print("Years until the perigee falls below 100 km.  A: orbit-averaged decay on the Vallado mean atmosphere (debris_orbits.lifetime).")
print("B: A over the cycle-averaged factor, 1.5 and 2 (the footnote of Table 4.7) and over the factor the October book implied (1.67 at 700 km, 1.44 at 600).")
print("C: the circular orbit at the perigee (debris_orbits.lifetime_quick), the pessimistic bound.\n")
hdr = f"{'case':46} {'orbit km':>16} {'piece':>13} {'A years':>10} {'A/1.5':>9} {'A/2':>9} {'A/Oct':>9} {'C years':>10}"
print(hdr)
for label, h, m_swept in cases:
    v, h_p, h_a = orbit_kick(h, m_swept)
    f_oct = OCT_FACTOR.get(round(h / 1e3))
    for name, B in BS:
        A = do.lifetime(h_p, h_a, B)
        C = do.lifetime_quick(h_p, B)
        fmt = lambda x: "over 1e5" if x == float('inf') else (f"{x:,.0f}" if x >= 100 else f"{x:,.1f}")
        print(f"{label:46} {h_p / 1e3:7.0f} x {h_a / 1e3:6.0f} {name:>13} {fmt(A):>10} {fmt(A / 1.5):>9} {fmt(A / 2):>9} {(fmt(A / f_oct) if f_oct else '-'):>9} {fmt(C):>10}")
    print(f"{'':46} (tangential speed {v:,.0f} m/s, radial {V:.0f} m/s)")
print("\nTable 4.12's lifetime rows: the full section on the worst orbit at 250 and 700 km, and the 8 cm pieces:")
for h in (250e3, 700e3):
    v, h_p, h_a = orbit_kick(h)
    print(f"  {h / 1e3:.0f} km: orbit {h_p / 1e3:.0f} x {h_a / 1e3:.0f} km; full section A {do.lifetime(h_p, h_a, B_ROD):,.1f} y, 8 cm A {do.lifetime(h_p, h_a, 60.0):,.2f} y; C {do.lifetime_quick(h_p, B_ROD):,.1f} and {do.lifetime_quick(h_p, 60.0):,.2f} y")
