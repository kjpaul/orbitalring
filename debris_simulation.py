#!/usr/bin/env python3
"""
Minimum Safe Altitude -- Orbital Ring Catastrophic Failure Analysis

Models the debris hazard from a catastrophic cable-casing collision on the
orbital ring, accounting for the asymmetric collision geometry (cable scrapes
along the upper casing wall) and three distinct debris populations:

  Population 1: Direct collision ejecta (ablated/spalled surface material)
    - Small mass fraction (1-5% of cable + effective casing)
    - Scattered at high velocity by the collision
    - High particulate conversion on reentry (10-50%)

  Population 2: Large cable fragments (hoop-stress breakup)
    - Dominant mass (~95% of cable)
    - Velocity near original cable velocity (delta_v ~ recoil speed)
    - Low particulate conversion on reentry (1-10%, size-dependent)

  Population 3: Ballistic casing debris (non-participating lower structure)
    - Falls from ring altitude at casing velocity + gravitational accel
    - Entry velocity too low for significant stratospheric particulate
    - Zero contribution to nuclear winter threshold

The minimum safe altitude is where the total stratospheric particulate
from all populations equals the 50 Mt nuclear winter threshold.

Usage:
    python debris_simulation.py [keywords] [options]

    Keywords:
      sweep        Mass vs. altitude for multiple beta values
      safe         Minimum safe altitude threshold crossing detail
      delta_v      Critical delta_v vs. altitude
      parametric   Safe altitude vs. beta exponent
      contour      2D contour: safe altitude vs. casing fraction & particulate eff
      populations  Breakdown of debris populations vs. altitude
      all          Generate all plots

    Options:
      --save           Save plots to graphs_debris/ (default)
      --show           Display plots interactively
      --beta=N         Override power-law exponent (e.g. --beta=1.8)
      --casing=N       Override casing participation fraction (0.0-1.0)
      --particulate=N  Override large-fragment particulate efficiency (0.0-1.0)
      --simple         Use original single-population model (Ch 13 baseline)
      --help, -h       Show this help

References:
  Housen & Holsapple (2011). Ejecta from impact craters.
    Icarus 211(1), 856-875.
  Grady & Kipp (1985). Geometric statistics and dynamic fragmentation.
    J. Appl. Phys. 58(3), 1210-1222.
  Melosh (1989). Impact Cratering: A Geologic Process.
    Oxford Univ. Press. Ch. 7 (ejecta scaling laws).
"""

import sys
import os
import math
import csv

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

import lim_physics as phys


# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

GM = 3.986004418e14          # Earth gravitational parameter (m^3/s^2)
R_EARTH = 6_378_137          # Earth equatorial radius (m)
OMEGA_SIDEREAL = 7.2921159e-5  # Earth sidereal rotation rate (rad/s)


# =============================================================================
# RING DESIGN PARAMETERS (configurable)
#
# These are the only fundamental parameters of the ring.  Everything else
# (cable mass, cable velocity, ring circumference, relative velocity at
# the collision, collision energy) is derived from these as a function of
# the ring altitude.
# =============================================================================

M_CASING_PER_M = 12_000.0    # Casing/load mass per meter (kg/m)
M_HARDWARE_PER_M = 4_076.0   # Cable hardware mass per meter (kg/m, fixed)
RETROGRADE = True            # True = retrograde cable, False = prograde
H_BASELINE = 250_000.0       # Baseline ring altitude (m), for reference

H_DRAG = 500_000.0           # Atmospheric drag threshold altitude (m)
M_THRESHOLD = 5.0e10         # Nuclear winter threshold (kg) = 50 Mt


# =============================================================================
# COLLISION GEOMETRY PARAMETERS
# =============================================================================

# What fraction of the casing mass participates in the collision?
# The cable scrapes along the upper wall. Most casing mass is below.
CASING_PARTICIPATION = 1.0   # Fraction (0.25 = quarter, 0.5 = half, 1.0 = full)

# CNT yarn material properties
#
# We treat the cable as a homogeneous CNT yarn defined by three measurable
# properties:  density, ultimate tensile strength, and strain at failure.
# Young's modulus, axial wave speed, and operating stress are then DERIVED
# so the chain (rho, sigma_ult, eps_break, safety) -> (E, c_axial, v_recoil)
# is fully transparent.
#
# The values here represent a future high-quality CNT yarn:
#   25 GPa ultimate strength at 5% failure strain implies E ~= 500 GPa.
# State-of-the-art lab CNT fibers (Behabtu 2013, Bai 2018) currently achieve
# 5-10 GPa at E ~= 200-300 GPa; the values below extrapolate to a yarn that
# is roughly 3x stronger and stiffer, consistent with the book's assumption
# of mature CNT manufacturing.
RHO_CNT          = 1_700.0   # Yarn density (kg/m^3)
SIGMA_ULTIMATE   = 25.0e9    # Ultimate tensile strength (Pa)
EPS_BREAK        = 0.05      # Strain at failure (5%)
SAFETY_FACTOR    = 2.0       # Operating stress = sigma_ult / safety_factor

# Derived elastic properties
E_CNT            = SIGMA_ULTIMATE / EPS_BREAK              # Young's modulus (Pa)
SIGMA_OPERATING  = SIGMA_ULTIMATE / SAFETY_FACTOR          # Operating stress (Pa)
C_AXIAL_CNT      = math.sqrt(E_CNT / RHO_CNT)              # Thin-rod axial wave speed (m/s)

# Transverse (shear) wave speed in the yarn.  Bundle slip and weak inter-tube
# coupling make G << E in CNT yarns, so this is treated as an empirical
# parameter rather than computed from E and Poisson's ratio.
C_TRANSVERSE_CNT = 4_000.0   # Transverse sound speed (m/s)

# Hoop-stress recoil velocity when the cable suddenly relaxes:
#   v_recoil = c_axial * eps_operating = sigma_op / (rho * c_axial)
# This is the particle velocity behind the relief wave that propagates
# inward from a fresh cable break.
V_RECOIL = SIGMA_OPERATING / (RHO_CNT * C_AXIAL_CNT)


# =============================================================================
# EJECTA DISTRIBUTION PARAMETERS
#
# The power-law M(>v) = M_retro * (v/v_min)^(-beta) describes the
# cumulative mass of collision ejecta with retrograde delta_v > v.
#
# To model MORE violent fragmentation (more mass at high delta_v):
#   - Decrease BETA toward 1.0
#
# To model LESS violent fragmentation:
#   - Increase BETA toward 2.5 or higher
# =============================================================================

BETA = 1.5                   # Power-law exponent (typical: 1.2-2.0)
V_MIN = 10.0                 # Minimum fragment delta_v (m/s)
F_RETROGRADE = 0.5           # Fraction of ejecta mass scattered retrograde


# =============================================================================
# MULTI-POPULATION MODEL PARAMETERS
# =============================================================================

# Population 1: Direct collision ejecta
EJECTA_MASS_FRACTION = 0.05  # 5% of collision mass becomes fine ejecta
EJECTA_PARTICULATE_EFF = 0.30  # 30% of reentry ejecta becomes stratospheric soot

# Population 2: Large cable fragments
LARGE_FRAG_PARTICULATE_EFF = 0.03  # 3% baseline (key uncertainty: 1-10%)

# Population 3: Ballistic casing (zero stratospheric contribution)
CASING_PARTICULATE_EFF = 0.0


# =============================================================================
# OUTPUT CONTROL
# =============================================================================

GRAPH_DPI = 200
GRAPH_DIR = "graphs_debris"
SAVE_GRAPHS = True


# =============================================================================
# SECTION 1: ORBITAL MECHANICS
# =============================================================================

def v_orbital(h):
    """Keplerian circular orbital velocity at altitude h (m) above surface."""
    return math.sqrt(GM / (R_EARTH + h))


def delta_v_critical(h_ring, h_drag=H_DRAG):
    """Minimum retrograde delta_v to drop a fragment's perigee to h_drag.

    For an orbit with apogee at r_ring and perigee at r_drag:
      a = (r_ring + r_drag) / 2
      v_apogee = sqrt(GM * (2/r_ring - 1/a))
      delta_v = v_orbit - v_apogee
    """
    if h_ring <= h_drag:
        return 0.0

    r_ring = R_EARTH + h_ring
    r_drag = R_EARTH + h_drag

    a = (r_ring + r_drag) / 2.0
    v_orb = math.sqrt(GM / r_ring)
    v_apogee = math.sqrt(GM * (2.0 / r_ring - 1.0 / a))

    return v_orb - v_apogee


# =============================================================================
# SECTION 2: RING DESIGN AT ALTITUDE
#
# The cable mass per meter, cable velocity, ring circumference, and
# relative collision velocity are not constants -- they all depend on
# the ring altitude (and the casing/load mass).  At each altitude we
# solve the quadratic cable sizing equation from lim_physics to get the
# self-consistent ring design.
# =============================================================================

def ring_state(h_ring, m_load=M_CASING_PER_M, m_hw=M_HARDWARE_PER_M,
               retrograde=RETROGRADE):
    """Self-consistent ring design at altitude h_ring.

    Returns a dict with all the altitude-dependent ring parameters:
      r          radius from Earth's center (m)
      L_ring     ring circumference (m)
      v_orbit    Keplerian circular orbital velocity (m/s)
      v_ground   ground-synchronous velocity at this radius (m/s)
      g_net      net downward acceleration on a stationary mass (m/s^2)
      m_struct   structural cable mass per meter (kg/m, from quadratic)
      m_rotor    total rotor mass = m_struct + m_hw (kg/m)
      m_load     casing/load mass per meter (kg/m)
      v_cable    cable velocity in inertial frame, signed (m/s)
                   negative = retrograde, positive = prograde
      v_rel      cable-casing relative velocity magnitude (m/s)
      M_total    total ring mass (kg) = (m_rotor + m_load) * L_ring
      m_cable_total  alias for m_rotor (kg/m)
      m_casing_total casing total mass (kg) = m_load * L_ring

    Cable velocity is set by momentum conservation.  The orbital ring's
    centre of mass is in a Keplerian circular orbit at v_orbit.  The
    casing rides at the slower (or opposite) ground-synchronous speed,
    so the cable must compensate:

      m_rotor * v_cable + m_load * v_ground = (m_rotor + m_load) * v_com

    For a retrograde cable, v_com = -v_orbit and v_cable is negative.
    For a prograde cable, v_com = +v_orbit.
    """
    r = R_EARTH + h_ring
    L = 2 * math.pi * r
    v_orb = math.sqrt(GM / r)
    v_ground = OMEGA_SIDEREAL * r
    g_net = GM / r**2 - OMEGA_SIDEREAL**2 * r

    if retrograde:
        delta_v_load = v_orb + v_ground
    else:
        delta_v_load = v_orb - v_ground

    m_struct = phys.calc_cable_mass(
        load_mass=m_load, m_hw=m_hw, delta_v=delta_v_load,
        v_orbit=v_orb, g_net=g_net, r_orbit=r,
        rho_cnt=RHO_CNT, sigma_cable=SIGMA_OPERATING)
    m_rotor = m_struct + m_hw

    if retrograde:
        # System COM at -v_orbit, casing at +v_ground, cable at -|v_cable|
        v_cable_mag = ((m_rotor + m_load) * v_orb + m_load * v_ground) / m_rotor
        v_cable = -v_cable_mag
        v_rel = v_cable_mag + v_ground
    else:
        # System COM at +v_orbit, casing and cable both prograde
        v_cable = ((m_rotor + m_load) * v_orb - m_load * v_ground) / m_rotor
        v_rel = v_cable - v_ground

    return {
        'r': r,
        'L_ring': L,
        'v_orbit': v_orb,
        'v_ground': v_ground,
        'g_net': g_net,
        'm_struct': m_struct,
        'm_rotor': m_rotor,
        'm_cable_total': m_rotor,
        'm_load': m_load,
        'v_cable': v_cable,
        'v_rel': v_rel,
        'M_total': (m_rotor + m_load) * L,
        'm_casing_total': m_load * L,
    }


# =============================================================================
# SECTION 2b: COLLISION ENERGY AND EJECTA
# =============================================================================

def reduced_mass_per_m(state, casing_frac=CASING_PARTICIPATION):
    """Reduced mass per meter for cable-casing collision (kg/m)."""
    m_eff = state['m_load'] * casing_frac
    m_cable = state['m_cable_total']
    return (m_cable * m_eff) / (m_cable + m_eff)


def collision_energy_per_m(state, casing_frac=CASING_PARTICIPATION):
    """Collision energy per meter (J/m) at this altitude."""
    mu = reduced_mass_per_m(state, casing_frac)
    return 0.5 * mu * state['v_rel']**2


def v_max_ejecta(state, casing_frac=CASING_PARTICIPATION):
    """Maximum ejecta velocity in the centre-of-mass frame (m/s).

    The lighter body (casing) takes the larger share of the relative
    velocity:  v_casing_cm = m_cable / (m_cable + m_casing_eff) * v_rel
    """
    m_eff = state['m_load'] * casing_frac
    m_cable = state['m_cable_total']
    return m_cable / (m_cable + m_eff) * state['v_rel']


def contact_dynamics(state, wall_thickness=0.005):
    """Contact dynamics for cable scraping along casing wall.

    A point on the cable surface is in contact for t = wall_thickness / v_rel.
    The transverse shock wave in the CNT propagates ~c_transverse * t into
    the cable surface during that interval, spalling a thin layer.
    """
    v_rel = state['v_rel']
    m_cable = state['m_cable_total']

    t_contact = wall_thickness / v_rel
    shock_depth = C_TRANSVERSE_CNT * t_contact
    a_cable = m_cable / RHO_CNT
    cable_side = math.sqrt(a_cable)
    m_shocked_cable = cable_side * shock_depth * RHO_CNT

    return {
        'shock_depth_m': shock_depth,
        't_contact_s': t_contact,
        'cable_side_m': cable_side,
        'm_shocked_cable_per_m': m_shocked_cable,
        'ablated_cable_fraction': m_shocked_cable / m_cable,
    }


# =============================================================================
# SECTION 3: SINGLE-POPULATION MODEL (original Ch 13 baseline)
# =============================================================================

def mass_below_drag_zone_simple(h_ring, beta=BETA, v_min=V_MIN, v_max=None,
                                f_retro=F_RETROGRADE, m_total=None,
                                state=None):
    """Total debris mass with perigee below H_DRAG (single-population model).

    Power-law ejecta: M(>v) = M_retro * (v/v_min)^(-beta)
    All debris reaching the drag zone is assumed 100% particulate.

    If m_total or v_max are not provided, they are derived from the
    self-consistent ring state at h_ring.
    """
    if m_total is None or v_max is None:
        if state is None:
            state = ring_state(h_ring)
        if m_total is None:
            m_total = state['M_total']
        if v_max is None:
            v_max = v_max_ejecta(state)

    if h_ring <= H_DRAG:
        return f_retro * m_total

    dv = delta_v_critical(h_ring)

    if dv <= v_min:
        return f_retro * m_total
    if dv >= v_max:
        return 0.0

    return f_retro * m_total * (dv / v_min) ** (-beta)


# =============================================================================
# SECTION 4: MULTI-POPULATION MODEL (revised Ch 12)
# =============================================================================

def particulate_from_pop1(h_ring, beta=BETA, casing_frac=CASING_PARTICIPATION,
                          ejecta_particulate_eff=EJECTA_PARTICULATE_EFF):
    """Population 1: Stratospheric particulate from direct collision ejecta."""
    state = ring_state(h_ring)
    m_casing_eff = state['m_load'] * casing_frac
    m_collision = (state['m_cable_total'] + m_casing_eff) * state['L_ring']
    m_ejecta_total = m_collision * EJECTA_MASS_FRACTION
    v_max = v_max_ejecta(state, casing_frac)

    m_reentry = mass_below_drag_zone_simple(
        h_ring, beta=beta, v_max=v_max, m_total=m_ejecta_total, state=state)

    return m_reentry * ejecta_particulate_eff


def particulate_from_pop2(h_ring, large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Population 2: Stratospheric particulate from large cable fragments.

    ~95% of cable mass stays in large fragments with delta_v ~ V_RECOIL
    (~408 m/s) from hoop-stress breakup. For isotropic recoil, the fraction
    with retrograde component > dv_crit determines reentry mass.
    """
    state = ring_state(h_ring)
    cd = contact_dynamics(state)
    m_large = (state['m_cable_total'] - cd['m_shocked_cable_per_m']) * state['L_ring']

    dv_crit = delta_v_critical(h_ring)

    if dv_crit <= 0:
        m_reentry = m_large
    elif dv_crit >= V_RECOIL * 3:
        m_reentry = 0.0
    elif dv_crit <= V_RECOIL:
        # Isotropic 3D recoil: fraction with 1D component > v
        frac = 0.5 * (1.0 - dv_crit / V_RECOIL)
        m_reentry = m_large * frac
    else:
        # Tail beyond V_RECOIL
        frac = 0.5 * (V_RECOIL / dv_crit) ** 2 * 0.1
        m_reentry = m_large * frac

    return m_reentry * large_frag_eff


def particulate_total(h_ring, beta=BETA, casing_frac=CASING_PARTICIPATION,
                      ejecta_particulate_eff=EJECTA_PARTICULATE_EFF,
                      large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Total stratospheric particulate from all populations (kg)."""
    p1 = particulate_from_pop1(h_ring, beta, casing_frac, ejecta_particulate_eff)
    p2 = particulate_from_pop2(h_ring, large_frag_eff)
    return p1 + p2


# =============================================================================
# SECTION 5: ROOT FINDING
# =============================================================================

def find_safe_altitude_simple(target_mass=M_THRESHOLD, beta=BETA,
                              h_low=501e3, h_high=20_000e3):
    """Find altitude where debris mass = target (single-population model)."""
    def objective(h):
        return mass_below_drag_zone_simple(h, beta=beta) - target_mass

    f_low, f_high = objective(h_low), objective(h_high)
    if f_low <= 0:
        return h_low
    if f_high >= 0:
        return h_high
    return brentq(objective, h_low, h_high, rtol=1e-12)


def find_safe_altitude_multi(target_mass=M_THRESHOLD, beta=BETA,
                             casing_frac=CASING_PARTICIPATION,
                             ejecta_particulate_eff=EJECTA_PARTICULATE_EFF,
                             large_frag_eff=LARGE_FRAG_PARTICULATE_EFF,
                             h_low=501e3, h_high=50_000e3):
    """Find altitude where total particulate = target (multi-population model)."""
    def objective(h):
        return particulate_total(h, beta, casing_frac,
                                 ejecta_particulate_eff, large_frag_eff) - target_mass

    f_low, f_high = objective(h_low), objective(h_high)
    if f_low <= 0:
        return h_low
    if f_high >= 0:
        return h_high
    return brentq(objective, h_low, h_high, rtol=1e-12)


def find_safe_altitude(beta=BETA, simple=False, **kwargs):
    """Find minimum safe altitude using selected model."""
    if simple:
        return find_safe_altitude_simple(beta=beta)
    return find_safe_altitude_multi(beta=beta, **kwargs)


# =============================================================================
# SECTION 6: CONSOLE OUTPUT
# =============================================================================

def print_parameters(beta=BETA, simple=False, casing_frac=CASING_PARTICIPATION,
                     large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Print simulation parameters and key results."""
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')

    sep = "=" * 70
    model = "SINGLE-POPULATION (Ch 13)" if simple else "MULTI-POPULATION (Ch 12)"

    print(sep)
    print(f"MINIMUM SAFE ALTITUDE -- {model}")
    print(sep)

    # Compute the minimum safe altitude up front so the ring parameters
    # can be displayed for the self-consistent design at that altitude.
    try:
        h_safe = find_safe_altitude(beta=beta, simple=simple,
                                    casing_frac=casing_frac,
                                    large_frag_eff=large_frag_eff)
    except (ValueError, RuntimeError):
        h_safe = None

    h_display = h_safe if h_safe is not None else H_BASELINE
    state = ring_state(h_display)

    print(f"\n  Ring Parameters (at h = {h_display/1e3:.0f} km):")
    print(f"    Cable mass/m        {state['m_cable_total']:>12,.0f} kg/m")
    print(f"    Casing mass/m       {state['m_load']:>12,.0f} kg/m")
    print(f"    Cable velocity      {abs(state['v_cable']):>12,.0f} m/s "
          f"({'retrograde' if state['v_cable'] < 0 else 'prograde'})")
    print(f"    Casing velocity     {state['v_ground']:>12,.0f} m/s (prograde)")
    print(f"    Relative velocity   {state['v_rel']:>12,.0f} m/s")
    print(f"    Ring circumference  {state['L_ring']:>12,.0f} m ({state['L_ring']/1e3:.0f} km)")
    print(f"    Total ring mass     {state['M_total']:>12.3e} kg")
    print(f"    Drag threshold      {H_DRAG/1e3:>12.0f} km")
    print(f"    Nuclear winter thr  {M_THRESHOLD:>12.2e} kg ({M_THRESHOLD/1e9:.0f} Mt)")

    if not simple:
        m_eff = state['m_load'] * casing_frac
        print(f"\n  Collision Geometry:")
        print(f"    Casing participation {casing_frac:>12.0%}")
        print(f"    Effective casing     {m_eff:>12,.0f} kg/m")
        print(f"    Reduced mass         {reduced_mass_per_m(state, casing_frac):>12,.0f} kg/m")
        print(f"    Collision energy     {collision_energy_per_m(state, casing_frac)/1e9:>12.1f} GJ/m")
        print(f"    V_max (CM frame)     {v_max_ejecta(state, casing_frac):>12.0f} m/s")

        cd = contact_dynamics(state)
        print(f"\n  Contact Dynamics (5 mm Al wall):")
        print(f"    Contact time/point   {cd['t_contact_s']*1e6:>12.2f} us")
        print(f"    Shock depth          {cd['shock_depth_m']*1000:>12.1f} mm")
        print(f"    Cable side           {cd['cable_side_m']:>12.1f} m")
        print(f"    Shocked cable frac   {cd['ablated_cable_fraction']*100:>12.2f} %")
        print(f"    V_recoil (hoop)      {V_RECOIL:>12.0f} m/s")

        print(f"\n  Population Parameters:")
        print(f"    Pop 1 ejecta frac    {EJECTA_MASS_FRACTION*100:>12.1f} %")
        print(f"    Pop 1 particulate    {EJECTA_PARTICULATE_EFF*100:>12.0f} %")
        print(f"    Pop 2 particulate    {large_frag_eff*100:>12.0f} %")
        print(f"    Pop 3 particulate    {CASING_PARTICULATE_EFF*100:>12.0f} %")

    print(f"\n  Ejecta Distribution:")
    print(f"    Beta (power law)    {beta:>12.2f}")
    print(f"    V_min               {V_MIN:>12.0f} m/s")
    print(f"    F_retrograde        {F_RETROGRADE:>12.2f}")

    # Altitude sweep table
    print(f"\n  {'Altitude':>10}  {'dv_crit':>10}  ", end="")
    if simple:
        print(f"{'Mass < 500km':>14}")
        print(f"  {'(km)':>10}  {'(m/s)':>10}  {'(Mt)':>14}")
        print(f"  {'-'*10}  {'-'*10}  {'-'*14}")
    else:
        print(f"{'Pop 1':>10}  {'Pop 2':>10}  {'Total':>10}")
        print(f"  {'(km)':>10}  {'(m/s)':>10}  {'(Mt)':>10}  {'(Mt)':>10}  {'(Mt)':>10}")
        print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for h_km in [600, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000]:
        h = h_km * 1e3
        dv = delta_v_critical(h)
        if simple:
            m = mass_below_drag_zone_simple(h, beta=beta)
            flag = " <--" if abs(m - M_THRESHOLD) / M_THRESHOLD < 0.5 else ""
            print(f"  {h_km:>10,}  {dv:>10.1f}  {m/1e9:>14.1f}{flag}")
        else:
            p1 = particulate_from_pop1(h, beta, casing_frac)
            p2 = particulate_from_pop2(h, large_frag_eff)
            total = p1 + p2
            flag = " <--" if abs(total - M_THRESHOLD) / M_THRESHOLD < 0.5 else ""
            print(f"  {h_km:>10,}  {dv:>10.1f}  {p1/1e9:>10.1f}"
                  f"  {p2/1e9:>10.1f}  {total/1e9:>10.1f}{flag}")

    if h_safe is not None:
        dv_safe = delta_v_critical(h_safe)
        print(f"\n  RESULT:")
        print(f"    Minimum safe altitude:  {h_safe/1e3:.1f} km")
        print(f"    Critical delta_v:       {dv_safe:.1f} m/s")
        print(f"    v_orbital at that alt:  {v_orbital(h_safe):.1f} m/s")
    else:
        print(f"\n  RESULT: No crossing found in search range")

    # Beta sensitivity
    print(f"\n  Sensitivity to beta:")
    print(f"    {'beta':>6}  {'Safe alt (km)':>13}  {'dv_crit (m/s)':>13}")
    print(f"    {'-'*6}  {'-'*13}  {'-'*13}")
    for b in [1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0]:
        try:
            hs = find_safe_altitude(beta=b, simple=simple,
                                    casing_frac=casing_frac,
                                    large_frag_eff=large_frag_eff)
            print(f"    {b:>6.1f}  {hs/1e3:>13.1f}  {delta_v_critical(hs):>13.1f}")
        except (ValueError, RuntimeError):
            print(f"    {b:>6.1f}  {'> 50,000':>13}  {'N/A':>13}")

    print()
    return h_safe


# =============================================================================
# SECTION 7: PLOTTING
# =============================================================================

def save_or_show(fig, filename, save=True, show=False):
    """Save figure to GRAPH_DIR or show interactively."""
    if save:
        os.makedirs(GRAPH_DIR, exist_ok=True)
        path = os.path.join(GRAPH_DIR, filename)
        fig.savefig(path, dpi=GRAPH_DPI, bbox_inches='tight')
        print(f"  Saved: {path}")
    if show:
        plt.show()
    plt.close(fig)


def plot_mass_vs_altitude(save=True, show=False, beta=BETA, simple=False,
                          casing_frac=CASING_PARTICIPATION,
                          large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Plot 1: Mass/particulate vs altitude for multiple beta values."""
    print("Generating: mass vs. altitude sweep...")

    h_km = np.linspace(550, 10000, 2000)
    betas = [1.2, 1.5, 1.8, 2.0]
    if beta not in betas:
        betas.append(beta)
        betas.sort()

    fig, ax = plt.subplots(figsize=(10, 6))

    for b in betas:
        if simple:
            mass = np.array([mass_below_drag_zone_simple(h*1e3, beta=b) for h in h_km])
        else:
            mass = np.array([particulate_total(h*1e3, beta=b, casing_frac=casing_frac,
                             large_frag_eff=large_frag_eff) for h in h_km])
        lw = 2.5 if abs(b - beta) < 0.01 else 1.2
        ax.plot(h_km, mass / 1e9, linewidth=lw, label=f"beta = {b:.1f}")

    ax.axhline(M_THRESHOLD / 1e9, color='red', linestyle='--', linewidth=1.5,
               label=f"Extinction threshold ({M_THRESHOLD/1e9:.0f} Mt)")

    try:
        h_safe = find_safe_altitude(beta=beta, simple=simple,
                                    casing_frac=casing_frac,
                                    large_frag_eff=large_frag_eff)
        ax.axvline(h_safe/1e3, color='gray', linestyle=':', linewidth=1.0, alpha=0.7)
        ax.plot(h_safe/1e3, M_THRESHOLD/1e9, 'ko', markersize=8, zorder=5)
        ax.annotate(f"  {h_safe/1e3:.0f} km (beta={beta})",
                    xy=(h_safe/1e3, M_THRESHOLD/1e9), fontsize=10, va='bottom')
    except (ValueError, RuntimeError):
        pass

    ax.set_yscale('log')
    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ylabel = "Debris Mass (Mt)" if simple else "Stratospheric Particulate (Mt)"
    ax.set_ylabel(ylabel, fontsize=12)
    title = "Single-Population" if simple else "Multi-Population"
    ax.set_title(f"{title}: Debris vs. Ring Altitude", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(500, 10000)
    ax.set_ylim(1e-1, 5e3)

    save_or_show(fig, "01-mass_vs_altitude.png", save, show)


def plot_safe_altitude(save=True, show=False, beta=BETA, simple=False,
                       casing_frac=CASING_PARTICIPATION,
                       large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Plot 2: Zoomed view of threshold crossing."""
    print("Generating: safe altitude detail...")

    try:
        h_safe = find_safe_altitude(beta=beta, simple=simple,
                                    casing_frac=casing_frac,
                                    large_frag_eff=large_frag_eff)
    except (ValueError, RuntimeError):
        print("  Skipping: no threshold crossing found.")
        return

    h_center = h_safe / 1e3
    h_lo = max(550, h_center - 500)
    h_hi = h_center + 500

    h_km = np.linspace(h_lo, h_hi, 1000)
    if simple:
        mass = np.array([mass_below_drag_zone_simple(h*1e3, beta=beta) for h in h_km])
    else:
        mass = np.array([particulate_total(h*1e3, beta=beta, casing_frac=casing_frac,
                         large_frag_eff=large_frag_eff) for h in h_km])
    mass_mt = mass / 1e9

    fig, ax = plt.subplots(figsize=(10, 6))
    label = "simple" if simple else "multi-pop"
    ax.plot(h_km, mass_mt, 'b-', linewidth=2.5, label=f"Debris ({label}, beta={beta})")
    ax.axhline(M_THRESHOLD/1e9, color='red', linestyle='--', linewidth=1.5,
               label=f"Extinction threshold ({M_THRESHOLD/1e9:.0f} Mt)")
    ax.axvline(h_safe/1e3, color='green', linestyle='-', linewidth=2.0, alpha=0.7,
               label=f"Min safe altitude = {h_safe/1e3:.0f} km")
    ax.plot(h_safe/1e3, M_THRESHOLD/1e9, 'r*', markersize=15, zorder=5)

    ax.set_yscale('log')
    y_lo = min(mass_mt) * 0.5
    y_hi = max(mass_mt) * 2.0
    ax.set_ylim(y_lo, y_hi)

    ax.axvspan(h_lo, h_safe/1e3, alpha=0.08, color='red')
    ax.axvspan(h_safe/1e3, h_hi, alpha=0.08, color='green')
    ax.text(h_lo + 20, y_lo * 1.5,
            "DANGER", fontsize=14, color='red', alpha=0.4, va='bottom')
    ax.text(h_hi - 20, y_lo * 1.5,
            "SAFE", fontsize=14, color='green', alpha=0.4, va='bottom', ha='right')

    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ylabel = "Debris Mass (Mt)" if simple else "Stratospheric Particulate (Mt)"
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title("Minimum Safe Altitude -- Threshold Crossing Detail", fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    save_or_show(fig, "02-safe_altitude.png", save, show)


def plot_delta_v(save=True, show=False):
    """Plot 3: Critical delta_v vs altitude with V_MAX and V_RECOIL lines."""
    print("Generating: delta_v vs. altitude...")

    h_km = np.linspace(550, 10000, 1000)
    dv = np.array([delta_v_critical(h * 1e3) for h in h_km])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(h_km, dv, 'b-', linewidth=2.0, label="delta_v to reach 500 km perigee")

    # V_max varies (weakly) with altitude; show curves rather than horizontal
    # lines so the altitude dependence is honest.
    for frac, color, ls in [(1.0, 'red', '--'), (0.5, 'orange', '-.'),
                            (0.25, 'green', ':')]:
        vm = np.array([v_max_ejecta(ring_state(h * 1e3), frac) for h in h_km])
        ax.plot(h_km, vm, color=color, linestyle=ls, linewidth=1.2,
                label=f"V_max ({frac:.0%} casing)")

    ax.axhline(V_RECOIL, color='purple', linestyle=':', linewidth=1.0,
               label=f"V_recoil (hoop stress) = {V_RECOIL:.0f} m/s")

    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ax.set_ylabel("Retrograde delta_v to Reach 500 km Perigee (m/s)", fontsize=12)
    ax.set_title("Critical Delta-v for Atmospheric Re-entry", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    save_or_show(fig, "03-delta_v_vs_altitude.png", save, show)


def plot_parametric(save=True, show=False, simple=False,
                    casing_frac=CASING_PARTICIPATION,
                    large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Plot 4: Minimum safe altitude vs beta exponent."""
    print("Generating: parametric beta sweep...")

    betas = np.linspace(1.0, 3.0, 200)
    h_safe_km = []
    for b in betas:
        try:
            hs = find_safe_altitude(beta=b, simple=simple, casing_frac=casing_frac,
                                    large_frag_eff=large_frag_eff)
            h_safe_km.append(hs / 1e3)
        except (ValueError, RuntimeError):
            h_safe_km.append(np.nan)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(betas, h_safe_km, 'b-', linewidth=2.5)

    try:
        h_def = find_safe_altitude(beta=BETA, simple=simple, casing_frac=casing_frac,
                                   large_frag_eff=large_frag_eff) / 1e3
        ax.plot(BETA, h_def, 'ro', markersize=10, zorder=5,
                label=f"Default: beta={BETA}, H={h_def:.0f} km")
    except (ValueError, RuntimeError):
        pass

    ax.axvspan(1.2, 2.0, alpha=0.08, color='blue',
               label="Typical hypervelocity range (1.2-2.0)")
    ax.set_xlabel("Power-Law Exponent (beta)", fontsize=12)
    ax.set_ylabel("Minimum Safe Altitude (km)", fontsize=12)
    title = "Single-Population" if simple else "Multi-Population"
    ax.set_title(f"{title}: Safe Altitude Sensitivity to Beta", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    save_or_show(fig, "04-parametric_beta.png", save, show)


def plot_contour(save=True, show=False, beta=BETA):
    """Plot 5: 2D contour of safe altitude vs casing fraction & particulate eff."""
    print("Generating: contour plot...")

    casing_fracs = np.linspace(0.1, 1.0, 30)
    part_effs = np.linspace(0.01, 0.10, 30)
    CF, PE = np.meshgrid(casing_fracs, part_effs)
    H_SAFE = np.zeros_like(CF)

    for i in range(CF.shape[0]):
        for j in range(CF.shape[1]):
            try:
                H_SAFE[i, j] = find_safe_altitude_multi(
                    beta=beta, casing_frac=CF[i, j],
                    large_frag_eff=PE[i, j]) / 1e3
            except (ValueError, RuntimeError):
                H_SAFE[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(10, 7))
    levels = [500, 750, 1000, 1250, 1500, 2000, 3000, 5000]
    cs = ax.contourf(CF * 100, PE * 100, H_SAFE, levels=20, cmap='RdYlGn')
    cl = ax.contour(CF * 100, PE * 100, H_SAFE, levels=levels,
                    colors='black', linewidths=0.8)
    ax.clabel(cl, fmt='%.0f km', fontsize=9)
    plt.colorbar(cs, ax=ax, label='Minimum Safe Altitude (km)')

    ax.plot(CASING_PARTICIPATION * 100, LARGE_FRAG_PARTICULATE_EFF * 100,
            'w*', markersize=15, markeredgecolor='black', markeredgewidth=1.0,
            label='Baseline parameters', zorder=5)

    ax.set_xlabel("Casing Participation (%)", fontsize=12)
    ax.set_ylabel("Large Fragment Particulate Efficiency (%)", fontsize=12)
    ax.set_title(f"Minimum Safe Altitude (beta={beta})", fontsize=13)
    ax.legend(fontsize=10, loc='upper right')

    save_or_show(fig, "05-contour_safe_altitude.png", save, show)


def plot_populations(save=True, show=False, beta=BETA,
                     casing_frac=CASING_PARTICIPATION,
                     large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Plot 6: Population breakdown vs altitude."""
    print("Generating: population breakdown...")

    h_km = np.linspace(550, 10000, 500)
    p1 = np.array([particulate_from_pop1(h*1e3, beta, casing_frac) for h in h_km])
    p2 = np.array([particulate_from_pop2(h*1e3, large_frag_eff) for h in h_km])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(h_km, (p1 + p2) / 1e9, 'k-', linewidth=2.5, label='Total particulate')
    ax.plot(h_km, p1 / 1e9, 'b--', linewidth=1.5, label='Pop 1: Collision ejecta')
    ax.plot(h_km, p2 / 1e9, 'r-.', linewidth=1.5, label='Pop 2: Large fragments')
    ax.axhline(M_THRESHOLD / 1e9, color='red', linestyle=':', linewidth=1.0,
               label=f"Threshold ({M_THRESHOLD/1e9:.0f} Mt)")

    ax.set_yscale('log')
    ax.set_xlabel("Ring Altitude (km)", fontsize=12)
    ax.set_ylabel("Stratospheric Particulate (Mt)", fontsize=12)
    ax.set_title("Debris Population Breakdown vs. Ring Altitude", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(500, 10000)
    ax.set_ylim(1e-2, 5e3)

    save_or_show(fig, "06-populations.png", save, show)


# =============================================================================
# SECTION 8: CSV EXPORT
# =============================================================================

def write_csv(beta=BETA, simple=False, casing_frac=CASING_PARTICIPATION,
              large_frag_eff=LARGE_FRAG_PARTICULATE_EFF):
    """Write results to CSV files."""
    os.makedirs(GRAPH_DIR, exist_ok=True)

    path = os.path.join(GRAPH_DIR, "debris_altitude_sweep.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        if simple:
            w.writerow(["altitude_km", "delta_v_crit_m_s", "mass_below_500km_Mt"])
        else:
            w.writerow(["altitude_km", "delta_v_crit_m_s",
                         "pop1_Mt", "pop2_Mt", "total_Mt"])
        for h_km in range(550, 10001, 10):
            h = h_km * 1e3
            dv = delta_v_critical(h)
            if simple:
                m = mass_below_drag_zone_simple(h, beta=beta)
                w.writerow([h_km, f"{dv:.2f}", f"{m/1e9:.4f}"])
            else:
                p1 = particulate_from_pop1(h, beta, casing_frac)
                p2 = particulate_from_pop2(h, large_frag_eff)
                w.writerow([h_km, f"{dv:.2f}", f"{p1/1e9:.4f}",
                            f"{p2/1e9:.4f}", f"{(p1+p2)/1e9:.4f}"])
    print(f"  Saved: {path}")

    path = os.path.join(GRAPH_DIR, "debris_beta_sensitivity.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["beta", "safe_altitude_km", "delta_v_crit_m_s"])
        for b100 in range(100, 301, 5):
            b = b100 / 100.0
            try:
                hs = find_safe_altitude(beta=b, simple=simple, casing_frac=casing_frac,
                                        large_frag_eff=large_frag_eff)
                w.writerow([f"{b:.2f}", f"{hs/1e3:.1f}", f"{delta_v_critical(hs):.1f}"])
            except (ValueError, RuntimeError):
                w.writerow([f"{b:.2f}", "N/A", "N/A"])
    print(f"  Saved: {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')

    show_graphs = []
    save = SAVE_GRAPHS
    show = False
    beta = BETA
    simple = False
    casing_frac = CASING_PARTICIPATION
    large_frag_eff = LARGE_FRAG_PARTICULATE_EFF

    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg == "--save":
                save = True
            elif arg == "--show":
                show = True
                save = False
            elif arg == "--simple":
                simple = True
            elif arg.startswith("--beta="):
                beta = float(arg.split("=")[1])
            elif arg.startswith("--casing="):
                casing_frac = float(arg.split("=")[1])
            elif arg.startswith("--particulate="):
                large_frag_eff = float(arg.split("=")[1])
            else:
                show_graphs.append(arg.lower())

    h_safe = print_parameters(beta=beta, simple=simple, casing_frac=casing_frac,
                              large_frag_eff=large_frag_eff)

    if not show_graphs:
        return

    keywords = set(show_graphs)
    do_all = "all" in keywords
    print("Generating graphs...")

    if do_all or "sweep" in keywords:
        plot_mass_vs_altitude(save, show, beta, simple, casing_frac, large_frag_eff)
    if do_all or "safe" in keywords:
        plot_safe_altitude(save, show, beta, simple, casing_frac, large_frag_eff)
    if do_all or "delta_v" in keywords:
        plot_delta_v(save, show)
    if do_all or "parametric" in keywords:
        plot_parametric(save, show, simple, casing_frac, large_frag_eff)
    if (do_all or "contour" in keywords) and not simple:
        plot_contour(save, show, beta)
    if (do_all or "populations" in keywords) and not simple:
        plot_populations(save, show, beta, casing_frac, large_frag_eff)

    if save:
        write_csv(beta, simple, casing_frac, large_frag_eff)

    print("Done.")


if __name__ == "__main__":
    main()
