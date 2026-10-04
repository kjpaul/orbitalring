#!/usr/bin/env python3
"""
Power Physics Module - Solar generation and HVDC distribution calculations

Pure physics functions for computing solar power generation around the
orbital ring and sizing the circumferential HVDC transmission system.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

import power_config as cfg


# =============================================================================
# SOLAR POWER GENERATION
# =============================================================================
#
# Physics model restored from Chapter 7 "Powering the Orbital Ring" (Vol III).
# This replaces an earlier oversimplified model that dropped two real effects:
#
#   1. Extended limb albedo. Even when a ring segment is past the geometric
#      day-night terminator (|phi| > pi/2), the Earth-facing back surface can
#      still see the illuminated portion of Earth's limb. This continues to
#      deliver albedo flux all the way around to true midnight.
#
#   2. Anti-reflection nanostructure coating. Standard AR coatings capture
#      ~89% of the geometric limit; advanced moth-eye and inverted-pyramid
#      nanostructures push this to ~96%, a net improvement of about 7.85%.
#      The multiplier is applied after the cell efficiency.
#
# The reference numbers from Ch7 are recovered exactly with this model:
#
#     albedo peak (noon):       378 W/m^2   = A * S0 * (R_E/r)^2
#     peak flux (noon):        1739 W/m^2   = S0 + 378
#     geometric collection:    45.7% (front-only + bifacial albedo)
#     AR-boosted collection:   49.3% = 45.7% * 1.0785
#     average electrical:       302 W/m^2   (text rounds to "about 300")
#     peak electrical:          844 W/m^2   = 1739 * 0.45 * 1.0785
#


def direct_solar_flux(phi):
    """Direct solar flux on the outward (radially-outward) panel face.

    A bifacial panel on the orbital ring has its front surface facing
    radially outward and its back surface facing Earth. At angle phi
    from solar noon, the outward face sees the sun at an incidence angle
    of phi. The cosine of the incidence angle gives the projected area
    factor. The front face only collects sunlight while |phi| < pi/2
    (the cosine is positive). Past pi/2 the outward face is shadowed
    by its own substrate; the back face takes over through the albedo
    term and, within the narrow sunlit limb beyond 90 deg, by direct
    sun. That back-face direct contribution is small and is rolled
    into the empirical 49% collection factor from Ch7.

    Args:
        phi: Angle from subsolar point (rad), range [-pi, pi]

    Returns:
        Direct solar flux (W/m^2) on the outward-facing surface.
    """
    cos_phi = math.cos(phi)
    if cos_phi <= 0.0:
        return 0.0
    return cfg.SOLAR_CONSTANT * cos_phi


def albedo_flux(phi):
    """Albedo flux on the inward (Earth-facing) panel face.

    At altitude h = 250 km, Earth fills almost the entire downward
    hemisphere. The maximum albedo flux, received when the panel sits
    directly between the sun and Earth and looks down at Earth's fully
    illuminated dayside, is

        f_albedo_peak = A * S0 * (R_E / r_orbit)^2
                      = 0.30 * 1361 * (6371/6621)^2
                      = 378 W/m^2

    As the panel moves around the ring, the fraction of Earth's visible
    disk that is illuminated follows the standard phase-angle formula
    for a disk partially lit from one side. At phase angle alpha the
    illuminated fraction is (1 + cos(alpha))/2:

        phi =   0 deg (noon)     -> full disk visible   ->  378 W/m^2
        phi = +/-90 deg (dawn)   -> half-lit disk       ->  189 W/m^2
        phi = 180 deg (midnight) -> dark side           ->    0 W/m^2

    This formula also gives the "extended limb albedo" the text
    describes: even inside Earth's geometric shadow, the back surface
    sees some illuminated limb, so the albedo goes smoothly to zero at
    true midnight rather than abruptly at the shadow boundary (~106 deg).

    Args:
        phi: Angle from subsolar point (rad)

    Returns:
        Albedo flux on the inward-facing surface (W/m^2)
    """
    phase_factor = 0.5 * (1.0 + math.cos(phi))   # 1 at noon, 0 at midnight
    if phase_factor <= 0.0:
        return 0.0
    return cfg.ALBEDO_PEAK * phase_factor


def panel_electrical_output(phi):
    """Total electrical output per m^2 of a bifacial panel at angle phi.

    Combines direct (outward face) and albedo (inward face) flux through
    the cell efficiency and the anti-reflection coating boost factor.

    The packing factor (cell area / panel area) is applied if set to
    anything other than 1.0; the Ch7 derivation treats "per square meter
    of solar array" as the full array footprint, so the restored model
    uses PANEL_PACKING = 1.0 by default.

    Args:
        phi: Angle from subsolar point (rad)

    Returns:
        Electrical power per m^2 of panel area (W/m^2)
    """
    f_direct = direct_solar_flux(phi)
    f_albedo = albedo_flux(phi)
    f_total = f_direct + f_albedo
    return f_total * cfg.CELL_EFFICIENCY * cfg.AR_BOOST * cfg.PANEL_PACKING


def compute_ring_power_profile(panel_width, n_points=None):
    """Compute power generation at each point around the ring.

    Args:
        panel_width: Total solar panel width (m)
        n_points: Number of angular sample points (default: cfg.N_POINTS)

    Returns:
        dict with:
            phi: list of angles (rad)
            p_gen_per_m: list of generation per metre of ring (W/m)
            p_avg_per_m2: orbit-averaged electrical output (W/m^2)
            p_peak_per_m2: peak electrical output (W/m^2)
            p_total: total generation (W)
    """
    if n_points is None:
        n_points = cfg.N_POINTS

    dphi = 2 * math.pi / n_points
    phi_list = []
    p_gen_per_m = []

    total_flux = 0.0
    peak_flux = 0.0

    for i in range(n_points):
        phi = -math.pi + (i + 0.5) * dphi
        flux = panel_electrical_output(phi)
        gen_per_m = flux * panel_width  # W per metre of ring

        phi_list.append(phi)
        p_gen_per_m.append(gen_per_m)
        total_flux += flux
        if flux > peak_flux:
            peak_flux = flux

    avg_flux = total_flux / n_points
    total_gen = avg_flux * panel_width * cfg.L_ARRAY

    return {
        'phi': phi_list,
        'p_gen_per_m': p_gen_per_m,
        'p_avg_per_m2': avg_flux,
        'p_peak_per_m2': peak_flux,
        'p_total': total_gen,
    }


# =============================================================================
# POWER DEMAND AND NET FLOW
# =============================================================================

def compute_power_flow(p_gen_per_m, demand_per_m, dphi):
    """Compute net power surplus/deficit and cumulative circumferential flow.

    Power flows from surplus regions (dayside) to deficit regions (nightside).
    The cumulative flow at each point is the integral of net surplus from
    phi = -pi to that point. The flow is then shifted so that the total
    integral is zero (power is conserved — what goes in one direction must
    come back the other way).

    Args:
        p_gen_per_m: list of generation per metre (W/m) at each angular point
        demand_per_m: uniform demand per metre of ring (W/m)
        dphi: angular step (rad)

    Returns:
        dict with:
            net_per_m: list of net surplus per metre (W/m), positive = surplus
            p_flow: list of cumulative power flow (W) at each cross-section
            p_flow_peak: peak absolute power flow (W)
    """
    n = len(p_gen_per_m)
    arc_step = cfg.R_ORBIT * dphi  # metres of ring per angular step

    net_per_m = []
    for i in range(n):
        net = p_gen_per_m[i] - demand_per_m
        net_per_m.append(net)

    # Cumulative power flow: integrate net surplus around the ring.
    # P_flow(phi) = integral from -pi to phi of net(phi') * ds
    # where ds = R_orbit * dphi is the arc length step.
    p_flow = []
    cumulative = 0.0
    for i in range(n):
        cumulative += net_per_m[i] * arc_step
        p_flow.append(cumulative)

    # Shift so the mean flow is zero (the ring is a closed loop —
    # the net integral must vanish, and we choose the reference so
    # the flow is symmetric about zero).
    mean_flow = sum(p_flow) / n
    p_flow = [f - mean_flow for f in p_flow]

    # Peak absolute flow
    p_flow_peak = max(abs(f) for f in p_flow)

    return {
        'net_per_m': net_per_m,
        'p_flow': p_flow,
        'p_flow_peak': p_flow_peak,
    }


# =============================================================================
# HVDC CONDUCTOR SIZING
# =============================================================================

def size_conductor(p_flow_peak, v_hvdc=None, loss_budget=None):
    """Size the HVDC conductor cross-section for a given peak flow and loss budget.

    The conductor sizing formula:
        A_total = (rho_e * L_trans * P_trans) / (eta_loss * V_pole^2)

    where:
        rho_e = 1/sigma = electrical resistivity
        L_trans = average transmission distance (quarter circumference)
        P_trans = peak transmission power
        eta_loss = fractional loss budget
        V_pole = half of pole-to-pole voltage

    Args:
        p_flow_peak: Peak power flow through any cross-section (W)
        v_hvdc: Pole-to-pole voltage (V), default cfg.V_HVDC
        loss_budget: Fractional loss target, default cfg.LOSS_BUDGET

    Returns:
        dict with conductor sizing results
    """
    if v_hvdc is None:
        v_hvdc = cfg.V_HVDC
    if loss_budget is None:
        loss_budget = cfg.LOSS_BUDGET

    v_pole = v_hvdc / 2.0
    rho_e = cfg.RHO_ELEC
    l_trans = cfg.L_RING / 4.0  # Average transmission distance (quarter ring)

    # Total conductor cross-section (both poles combined)
    a_total = (rho_e * l_trans * p_flow_peak) / (loss_budget * v_pole ** 2)

    # Per-pole cross-section
    a_per_pole = a_total / 2.0

    # Number of cables: each pole needs ceil(A_per_pole / A_cable)
    # Bipolar system requires equal cables per pole
    a_cable = math.pi * cfg.CABLE_DIAMETER ** 2 / 4.0
    n_cables_per_pole = math.ceil(a_per_pole / a_cable)
    n_cables_total = 2 * n_cables_per_pole

    # Conductor mass per metre of ring
    m_per_m = a_total * cfg.RHO_CONDUCTOR

    # Peak current per pole
    i_peak = p_flow_peak / (2.0 * v_pole)  # Each pole carries half the power

    # Actual resistance per metre (both poles in series for the loop)
    r_per_m = rho_e / a_per_pole  # ohm/m per pole

    # Actual loss at peak flow
    # Loss = I^2 * R * L_trans for each pole, times 2 poles
    p_loss_peak = 2.0 * i_peak ** 2 * r_per_m * l_trans
    loss_fraction = p_loss_peak / p_flow_peak if p_flow_peak > 0 else 0.0

    # Current density
    j_peak = i_peak / a_per_pole if a_per_pole > 0 else 0.0

    return {
        'a_total': a_total,
        'a_per_pole': a_per_pole,
        'a_cable': a_cable,
        'n_cables_per_pole': n_cables_per_pole,
        'n_cables_total': n_cables_total,
        'm_per_m': m_per_m,
        'm_total': m_per_m * cfg.L_RING,
        'i_peak': i_peak,
        'j_peak': j_peak,
        'r_per_m': r_per_m,
        'p_loss_peak': p_loss_peak,
        'loss_fraction': loss_fraction,
        'v_hvdc': v_hvdc,
        'v_pole': v_pole,
        'l_trans': l_trans,
        'loss_budget': loss_budget,
    }
