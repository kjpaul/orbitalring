"""
debris_orbits.py -- Where the released cable goes after a catastrophic break.

Companion to debris_simulation.py (Volume IV, Chapter 4).  That script's
Population 2 model started the large cable fragments from the circular
orbital velocity at the ring altitude.  The cable does not move at circular
velocity.  It moves faster than circular (that excess is what carries the
casing's weight and the hoop tension), so when the cable is set free the
break point is the PERIGEE of the fragment's orbit, not the apogee, and the
fragments go up, not down.  This module does the orbital mechanics properly
and replaces the "critical delta-v from circular" criterion for Population 2
with a perigee-and-drag-lifetime criterion.

Contents
  1. Released-cable orbit: perigee and apogee for a given ring altitude and
     recoil components (tangential and radial).
  2. Momentum bookkeeping: how much casing mass the cable would have to
     sweep up before its perigee reaches the drag zone.
  3. Atmosphere model (Vallado exponential table, extended above 1,000 km).
  4. Orbit-averaged drag decay and lifetime for a fragment of given
     ballistic coefficient on an elliptical orbit.
  5. Tables for Chapter 4.

Run:  python3 debris_orbits.py            (prints the Chapter 4 tables)
      python3 debris_orbits.py --save     (also writes the SVG figures)
"""

import math
import sys
import numpy as np

import debris_simulation as ds
from debris_simulation import GM, R_EARTH, OMEGA_SIDEREAL, ring_state, V_RECOIL, RHO_CNT

YEAR = 365.25 * 86400.0

# -----------------------------------------------------------------------------
# 1. Released-cable orbit
# -----------------------------------------------------------------------------

def orbit_from_state(r, v_t, v_r=0.0):
    """Perigee and apogee altitude (m) of the orbit through radius r with
    tangential speed v_t and radial speed v_r.  Returns (h_p, h_a, a, e);
    h_a is inf for an unbound orbit."""
    v2 = v_t * v_t + v_r * v_r
    eps = v2 / 2.0 - GM / r
    h = r * v_t
    if eps >= 0:
        a = float('inf'); e = 1.0
        p = h * h / GM
        return p / (1 + e) - R_EARTH, float('inf'), a, e
    a = -GM / (2.0 * eps)
    e = math.sqrt(max(0.0, 1.0 - h * h / (GM * a)))
    return a * (1 - e) - R_EARTH, a * (1 + e) - R_EARTH, a, e


def released_cable(h_ring, dv_tan=0.0, dv_rad=0.0, retrograde=True):
    """Orbit of a cable fragment released at ring altitude h_ring (m) with
    an extra tangential speed dv_tan (positive = faster) and radial speed
    dv_rad (m/s).  Returns dict with v_cable, v_circ, excess, h_p, h_a."""
    s = ring_state(h_ring, retrograde=retrograde)
    r = s['r']
    v_c = abs(s['v_cable'])
    v_circ = s['v_orbit']
    v_t = v_c + dv_tan
    h_p, h_a, a, e = orbit_from_state(r, v_t, dv_rad)
    v_esc = v_circ * math.sqrt(2.0)
    return {'v_cable': v_c, 'v_circ': v_circ, 'v_esc': v_esc,
            'excess': v_c - v_circ, 'v_frag': v_t, 'h_p': h_p, 'h_a': h_a,
            'a': a, 'e': e, 'state': s}


def slowdown_to_perigee(h_ring, h_target, retrograde=True):
    """Tangential speed loss (m/s) a cable fragment would need, from its
    actual cable speed, to bring its perigee down to h_target (m)."""
    s = ring_state(h_ring, retrograde=retrograde)
    r = s['r']; rp = R_EARTH + h_target
    if rp >= r:
        return 0.0
    a = 0.5 * (r + rp)
    v_apo = math.sqrt(GM * (2.0 / r - 1.0 / a))
    return abs(s['v_cable']) - v_apo


# -----------------------------------------------------------------------------
# 2. Momentum bookkeeping
# -----------------------------------------------------------------------------

def swept_cable_speed(h_ring, m_swept, dv_tan=0.0, retrograde=True):
    """Speed of the cable (m/s, inertial, magnitude) after it has swept up
    m_swept kg/m of ground-synchronous casing and carried it along.  The
    cable's own tangential recoil dv_tan is applied before the sweep-up."""
    s = ring_state(h_ring, retrograde=retrograde)
    m_c = s['m_cable_total']
    v_c = abs(s['v_cable']) + dv_tan
    v_g = s['v_ground']
    if retrograde:
        # cable westward, casing eastward: opposite directions
        return (m_c * v_c - m_swept * v_g) / (m_c + m_swept)
    return (m_c * v_c + m_swept * v_g) / (m_c + m_swept)


def casing_mass_to_reach_perigee(h_ring, h_target, dv_tan=0.0, retrograde=True):
    """kg/m of casing the cable must sweep up (fully, to its own speed)
    before its perigee drops to h_target.  None if it never does."""
    s = ring_state(h_ring, retrograde=retrograde)
    r = s['r']; rp = R_EARTH + h_target
    a = 0.5 * (r + rp)
    v_need = math.sqrt(GM * (2.0 / r - 1.0 / a))
    m_c = s['m_cable_total']; v_c = abs(s['v_cable']) + dv_tan; v_g = s['v_ground']
    if retrograde:
        # (m_c v_c - m v_g)/(m_c + m) = v_need  ->  m = m_c (v_c - v_need)/(v_need + v_g)
        m = m_c * (v_c - v_need) / (v_need + v_g)
    else:
        m = m_c * (v_c - v_need) / (v_need - v_g)
    return m if m > 0 else None


# -----------------------------------------------------------------------------
# 3. Atmosphere (Vallado, Fundamentals of Astrodynamics, Table 8-4; mean solar
#    conditions).  Above 1,000 km the last scale height is extended.  The real
#    density above 1,000 km varies by an order of magnitude with solar
#    activity, so lifetimes computed up there are order-of-magnitude only.
# -----------------------------------------------------------------------------

_ATM = [  # (base altitude km, density kg/m^3, scale height km)
    (0, 1.225, 7.249), (25, 3.899e-2, 6.349), (30, 1.774e-2, 6.682),
    (40, 3.972e-3, 7.554), (50, 1.057e-3, 8.382), (60, 3.206e-4, 7.714),
    (70, 8.770e-5, 6.549), (80, 1.905e-5, 5.799), (90, 3.396e-6, 5.382),
    (100, 5.297e-7, 5.877), (110, 9.661e-8, 7.263), (120, 2.438e-8, 9.473),
    (130, 8.484e-9, 12.636), (140, 3.845e-9, 16.149), (150, 2.070e-9, 22.523),
    (180, 5.464e-10, 29.740), (200, 2.789e-10, 37.105), (250, 7.248e-11, 45.546),
    (300, 2.418e-11, 53.628), (350, 9.518e-12, 53.298), (400, 3.725e-12, 58.515),
    (450, 1.585e-12, 60.828), (500, 6.967e-13, 63.822), (600, 1.454e-13, 71.835),
    (700, 3.614e-14, 88.667), (800, 1.170e-14, 124.64), (900, 5.245e-15, 181.05),
    (1000, 3.019e-15, 268.00),
]

def rho_atm(h_m):
    """Atmospheric density (kg/m^3) at altitude h_m (m)."""
    h = h_m / 1e3
    if h <= 0:
        return _ATM[0][1]
    row = _ATM[0]
    for r in _ATM:
        if h >= r[0]:
            row = r
    h0, rho0, H = row
    return rho0 * math.exp(-(h - h0) / H)

def scale_height(h_m):
    h = h_m / 1e3
    row = _ATM[0]
    for r in _ATM:
        if h >= r[0]:
            row = r
    return row[2] * 1e3


# -----------------------------------------------------------------------------
# 4. Drag decay and lifetime
# -----------------------------------------------------------------------------

def ballistic_coefficient(side_m, cd=2.2):
    """m/(C_D A) for a long square-section CNT rod seen broadside (kg/m^2).
    A rod of side s and length L has mass rho s^2 L and broadside area s L."""
    return RHO_CNT * side_m / cd


def decay_per_orbit(a, e, B, n_steps=180):
    """Change in semi-major axis and eccentricity over one orbit from drag,
    by quadrature over eccentric anomaly (King-Hele orbit averaging).
    Drag acceleration f = rho v^2 / (2B) opposite to the velocity;
    Gauss planetary equations for a tangential perturbation:
      da/dt = 2 a^2 v f / GM,   de/dt = 2 (e + cos nu) f / v."""
    E = 2 * np.pi * (np.arange(n_steps) + 0.5) / n_steps
    r = a * (1 - e * np.cos(E))
    h = r - R_EARTH
    if np.any(h < 0):
        return -a, 0.0
    rho = np.array([rho_atm(x) for x in h])
    v = np.sqrt(GM * (2.0 / r - 1.0 / a))
    n_mean = math.sqrt(GM / a**3)
    dt = (r / (n_mean * a)) * (2 * np.pi / n_steps)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    f = rho * v * v / (2.0 * B)
    da = float(np.sum(-(2.0 * a * a * v / GM) * f * dt))
    de = float(np.sum(-(2.0 * (e + np.cos(nu)) / v) * f * dt))
    return da, de


def lifetime(h_p, h_a, B, h_reentry=100e3, max_years=1e5, verbose=False):
    """Years until perigee falls below h_reentry for a fragment with ballistic
    coefficient B (kg/m^2) starting on the orbit h_p x h_a (m).  Orbit-averaged
    decay with adaptive batching of orbits.  Returns inf if longer than
    max_years."""
    rp = R_EARTH + h_p; ra = R_EARTH + h_a
    a = 0.5 * (rp + ra); e = (ra - rp) / (ra + rp)
    t = 0.0
    while True:
        hp_now = a * (1 - e) - R_EARTH
        if hp_now <= h_reentry:
            return t / YEAR
        if t > max_years * YEAR:
            return float('inf')
        T = 2 * math.pi * math.sqrt(a**3 / GM)
        da, de = decay_per_orbit(a, e, B)
        if da >= 0:
            da = -1e-12
        dhp = abs(da * (1 - e) - a * de)
        # quick exit: if even the current (slowest) rate needs > max_years to drop 50 km
        if dhp > 0 and (50e3 / dhp) * T > max_years * YEAR * 3:
            return float('inf')
        n = max(1, int(min(0.01 * a / abs(da), 2000.0 / max(dhp, 1e-12))))
        n = min(n, 5_000_000)
        a += n * da; e = max(0.0, e + n * de); t += n * T
        if verbose:
            print(f"t={t/YEAR:9.2f} yr  hp={(a*(1-e)-R_EARTH)/1e3:8.1f} km  ha={(a*(1+e)-R_EARTH)/1e3:8.1f} km  n={n}")


def lifetime_quick(h, B):
    """Closed-form estimate for a circular orbit: tau = B H / (rho a v)."""
    r = R_EARTH + h
    v = math.sqrt(GM / r)
    return B * scale_height(h) / (rho_atm(h) * r * v) / YEAR


# -----------------------------------------------------------------------------
# 5. Chapter 4 tables
# -----------------------------------------------------------------------------

def print_tables():
    print("=" * 78)
    print("RELEASED CABLE: WHERE IT GOES  (retrograde cable, V_RECOIL = %.0f m/s)" % V_RECOIL)
    print("=" * 78)
    print(f"\n{'h_ring':>7} {'v_cable':>8} {'v_circ':>7} {'excess':>7} {'v_esc':>7} "
          f"{'apogee':>8} {'apo -rec':>8} {'apo +rec':>8} {'min peri':>8} {'to 500km':>8}")
    print(f"{'(km)':>7} {'(m/s)':>8} {'(m/s)':>7} {'(m/s)':>7} {'(m/s)':>7} "
          f"{'(km)':>8} {'(km)':>8} {'(km)':>8} {'(km)':>8} {'(m/s)':>8}")
    for h_km in [250, 500, 600, 800, 1000, 1500, 2000, 2500, 3000, 5000]:
        h = h_km * 1e3
        n = released_cable(h)
        lo = released_cable(h, -V_RECOIL); hi = released_cable(h, +V_RECOIL)
        worst = released_cable(h, -V_RECOIL, V_RECOIL)
        sd = slowdown_to_perigee(h, 500e3)
        print(f"{h_km:>7,} {n['v_cable']:>8,.0f} {n['v_circ']:>7,.0f} {n['excess']:>7,.0f} {n['v_esc']:>7,.0f} "
              f"{n['h_a']/1e3:>8,.0f} {lo['h_a']/1e3:>8,.0f} {hi['h_a']/1e3:>8,.0f} {worst['h_p']/1e3:>8,.0f} {sd:>8,.0f}")

    print("\nMomentum: casing mass (kg/m) the cable must sweep up to bring its perigee to 500 km")
    for h_km in [600, 800, 1000, 1500, 2000]:
        h = h_km * 1e3
        m0 = casing_mass_to_reach_perigee(h, 500e3)
        m1 = casing_mass_to_reach_perigee(h, 500e3, -V_RECOIL)
        print(f"   {h_km:>5,} km: nominal {m0 if m0 is None else round(m0):>7}   slow half (-recoil) {m1 if m1 is None else round(m1):>7}   (casing total 12,000; roof ~3,000)")

    import ring_altitude as _ring
    H_DESIGN = _ring.ALTITUDE
    print(f"\nCable speed after sweeping up the casing roof (3,000 kg/m) and the whole casing (12,000 kg/m), {H_DESIGN/1e3:,.0f} km:")
    for m in [0, 3000, 6000, 12000]:
        for dv in [-V_RECOIL, 0, V_RECOIL]:
            v = swept_cable_speed(H_DESIGN, m, dv)
            hp, ha, _, _ = orbit_from_state(R_EARTH + H_DESIGN, v)
            print(f"   swept {m:>6,} kg/m, recoil {dv:+5.0f}: v = {v:,.0f} m/s, perigee {hp/1e3:,.0f} km, apogee {ha/1e3:,.0f} km")

    print("\n" + "=" * 78)
    print("DRAG LIFETIME of a dense cable fragment (B = m/(C_D A))")
    print("=" * 78)
    s2000 = ring_state(H_DESIGN); side = math.sqrt(s2000['m_struct'] / RHO_CNT)
    B_rod = ballistic_coefficient(side)
    print(f"  {side:.2f} m square CNT rod broadside: m/A = {RHO_CNT*side:,.0f} kg/m^2, B = {B_rod:,.0f} kg/m^2 (C_D 2.2)")
    print(f"  {'perigee':>8} {'circular':>12} {'ellipse to 2x':>14} {'B=600 (0.8 m)':>14} {'B=60 (8 cm)':>12}")
    for hp_km in [150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000]:
        hp = hp_km * 1e3
        t_circ = lifetime(hp, hp, B_rod)
        t_ell = lifetime(hp, 2 * hp + 1500e3, B_rod)
        t_mid = lifetime(hp, hp, 600.0)
        t_small = lifetime(hp, hp, 60.0)
        def f(t):
            return "> 100 kyr" if t == float('inf') else (f"{t*365.25:.0f} d" if t < 1 else f"{t:,.0f} yr")
        print(f"  {hp_km:>6,} km {f(t_circ):>12} {f(t_ell):>14} {f(t_mid):>14} {f(t_small):>12}")

    print("\nThe 250 km case: released-cable orbits and their lifetimes (B = rod)")
    for dv_t, dv_r, label in [(0, 0, "nominal"), (-V_RECOIL, 0, "slow half"), (V_RECOIL, 0, "fast half"),
                              (0, V_RECOIL, "radial recoil"), (-V_RECOIL, V_RECOIL, "slow + radial")]:
        o = released_cable(250e3, dv_t, dv_r)
        t = lifetime(max(o['h_p'], 0), o['h_a'], B_rod) if o['h_p'] > 100e3 else 0.0
        print(f"   {label:>14}: perigee {o['h_p']/1e3:6.0f} km, apogee {o['h_a']/1e3:6.0f} km, lifetime {t:.2f} yr (rod); "
              f"{lifetime(max(o['h_p'],0), o['h_a'], 600.0) if o['h_p']>100e3 else 0:.3f} yr (B=600)")

    for H_CASE in sorted({H_DESIGN, 2000e3}):
      print(f"\nThe {H_CASE/1e3:,.0f} km case:")
      for dv_t, dv_r, label in [(0, 0, "nominal"), (-V_RECOIL, 0, "slow half"), (V_RECOIL, 0, "fast half"), (0, V_RECOIL, "radial recoil"), (-V_RECOIL, V_RECOIL, "slow + radial"), (V_RECOIL, V_RECOIL, "fast + radial")]:
        o = released_cable(H_CASE, dv_t, dv_r)
        tl = [lifetime(max(o['h_p'], 1), o['h_a'], B) for B in (B_rod, 600.0, 60.0)]
        print(f"   {label:>14}: perigee {o['h_p']/1e3:6.0f} km, apogee {o['h_a']/1e3:6.0f} km, e = {o['e']:.3f}, lifetime yr rod/0.8 m/8 cm: " + " / ".join("inf" if t == float("inf") else f"{t:,.0f}" for t in tl))
    print("\n(2,000 km short form, as before)")
    for dv_t, dv_r, label in [(0, 0, "nominal"), (-V_RECOIL, 0, "slow half"), (V_RECOIL, 0, "fast half"),
                              (-V_RECOIL, V_RECOIL, "slow + radial")]:
        o = released_cable(2000e3, dv_t, dv_r)
        print(f"   {label:>14}: perigee {o['h_p']/1e3:6.0f} km, apogee {o['h_a']/1e3:6.0f} km, e = {o['e']:.3f}")

    print("\n" + "=" * 78)
    print("POPULATION 1 (collision ejecta) prompt reentry particulate, 30% soot, full casing participation")
    print("=" * 78)
    print(f"  {'h_ring':>7} " + " ".join(f"{'beta='+str(b):>10}" for b in [1.0, 1.2, 1.5, 2.0, 3.0]))
    for h_km in [250, 400, 500, 600, 700, 800, 1000, 1500, 2000]:
        h = h_km * 1e3
        row = [ds.particulate_from_pop1(h, beta=b) / 1e9 for b in [1.0, 1.2, 1.5, 2.0, 3.0]]
        print(f"  {h_km:>5,} km " + " ".join(f"{x:>10.1f}" for x in row) + "  Mt")
    print("  (quarter casing participation, beta 1.5):",
          ", ".join(f"{h}: {ds.particulate_from_pop1(h*1e3, beta=1.5, casing_frac=0.25)/1e9:.1f} Mt" for h in [250, 600, 800, 2000]))


if __name__ == "__main__":
    print_tables()
