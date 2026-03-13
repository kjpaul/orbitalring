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

def direct_solar_flux(phi):
    """Direct solar flux on a radially-outward panel at angle phi from subsolar point.

    Args:
        phi: Angle from subsolar point (rad), range [-pi, pi]

    Returns:
        Solar flux (W/m^2) on the outward-facing surface.
        Zero when panel faces away from the sun (|phi| > pi/2)
        or when in Earth's shadow.
    """
    # In Earth's shadow?
    if abs(phi) > cfg.SHADOW_HALF_ANGLE:
        return 0.0
    # Facing away from sun?
    cos_phi = math.cos(phi)
    if cos_phi <= 0:
        return 0.0
    return cfg.SOLAR_CONSTANT * cos_phi


def albedo_flux(phi):
    """Albedo flux on the inward-facing (Earth-facing) panel surface.

    The albedo contribution depends on:
    1. Earth's albedo coefficient
    2. View factor (Earth fills most of the downward hemisphere at 250 km)
    3. Whether the Earth surface below is illuminated

    At 250 km altitude, the Earth subtends a large solid angle.
    The view factor from an infinite flat plate at 250 km to Earth is:
        F = (R_E / r_orbit)^2 ~ 0.926

    The illuminated fraction of the visible Earth surface below depends
    on the panel's position relative to the subsolar point.

    Args:
        phi: Angle from subsolar point (rad)

    Returns:
        Albedo flux on the inward-facing surface (W/m^2)
    """
    # In Earth's shadow — no albedo either (Earth surface below is dark)
    if abs(phi) > cfg.SHADOW_HALF_ANGLE:
        return 0.0

    # View factor: fraction of hemisphere filled by Earth
    view_factor = (cfg.R_EARTH / cfg.R_ORBIT) ** 2

    # The Earth surface below the panel receives solar flux proportional
    # to the cosine of the solar zenith angle at that point.
    # For a point at angle phi on the ring, the Earth surface directly
    # below sees the sun at zenith angle phi.
    cos_phi = math.cos(phi)
    if cos_phi <= 0:
        # Earth surface below is in darkness — minimal albedo
        return 0.0

    # Albedo flux = solar constant * albedo * view factor * cos(zenith)
    return cfg.SOLAR_CONSTANT * cfg.EARTH_ALBEDO * view_factor * cos_phi


def panel_electrical_output(phi):
    """Total electrical output per m^2 of panel area at angle phi.

    Combines direct (outward) and albedo (inward) flux through cell efficiency.

    Args:
        phi: Angle from subsolar point (rad)

    Returns:
        Electrical power per m^2 of panel area (W/m^2)
    """
    f_direct = direct_solar_flux(phi)
    f_albedo = albedo_flux(phi)
    return (f_direct + f_albedo) * cfg.CELL_EFFICIENCY * cfg.PANEL_PACKING


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
    total_gen = avg_flux * panel_width * cfg.L_RING

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
