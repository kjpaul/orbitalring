#!/usr/bin/env python3
"""
Levitation Physics Module - Magnetic field, force, and thermal calculations

This module contains the physics for:
  - Halbach array magnetic field profile (exponential decay with gap)
  - Solenoid magnetic field profile (for comparison)
  - EML attractive force between Halbach array and ferromagnetic plate
  - Repulsive eddy current force from safety Halbach arrays
  - Thermal analysis: eddy current heating and radiative cooling
  - Material property lookups

The key physics:

1. HALBACH ARRAY field:
   A Halbach array concentrates the magnetic field on one side (the "front")
   and cancels it on the other side (the "back"). The field on the front side
   decays exponentially with distance from the array surface:

     B(z) = B_surface * exp(-k * z)

   where k = 2*pi / lambda is the spatial wavenumber and lambda is the
   Halbach wavelength (the repeat distance of the magnet pattern). Short
   wavelengths give stronger surface fields but faster decay. This is
   exactly what we want for the safety backstop: strong force at close
   range, negligible field at the nominal gap (so no eddy current drag
   during normal operations).

2. EML ATTRACTIVE FORCE:
   The attractive force between a DC magnetic field source and a ferromagnetic
   plate follows the Maxwell stress tensor:

     P = B_gap^2 / (2 * mu_0)    [Pa, pressure]

   where B_gap is the field in the air gap. For an EML bearing with a
   Halbach array source, B_gap depends on the gap distance, the coil
   current, and the array geometry.

3. EDDY CURRENT REPULSION (safety):
   When a conducting plate moves relative to a Halbach array (or when the
   gap changes rapidly), eddy currents are induced. For a stationary plate
   in a time-varying field (safety engagement transient), the repulsive
   force depends on the rate of field change and the plate conductivity.

   For the safety backstop, the cable is NOT moving laterally through the
   Halbach array. Instead, the cable is drifting vertically (closing the gap)
   and the increasing field induces eddy currents that oppose the motion.
   The force is:

     F_repulsive proportional to sigma * (dB/dt)^2 * volume

   But we can also model it as a quasi-static image current repulsion
   for a superconducting Halbach array near a conducting plate, which gives
   a force that increases rapidly as the gap closes.

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math
from levitation_config import (
    MU0, STEFAN_BOLTZMANN, T_SPACE, T_CASING_WALL,
    REACTION_PLATE_MATERIALS
)


# =============================================================================
# HALBACH ARRAY MAGNETIC FIELD
# =============================================================================

def halbach_surface_field(n_turns, i_current, coil_width, wavelength,
                          coils_per_meter=None):
    """
    Peak magnetic field at the surface of a Halbach array.

    For an array of HTS coils arranged in a Halbach pattern, the surface
    field is determined by the sheet current density (total ampere-turns
    per meter of array length). The Halbach arrangement concentrates the
    field on one side, giving approximately sqrt(2) enhancement over a
    simple solenoid arrangement.

    The key formula: for a current sheet with linear current density K (A/m),
    the field at the surface is:

      B_surface = mu_0 * K * halbach_factor

    where K = coils_per_meter * n_turns * i_current.

    For a 4-segment Halbach (90 deg rotation per segment), the enhancement
    factor is sqrt(2) ≈ 1.414. The field then decays exponentially with
    distance from the surface:

      B(z) = B_surface * exp(-2*pi*z / lambda)

    Parameters
    ----------
    n_turns : int
        Turns per coil
    i_current : float
        Operating current (A)
    coil_width : float
        Coil width (m), used as fallback for single-coil estimate
    wavelength : float
        Halbach wavelength (m), the repeat distance
    coils_per_meter : float or None
        Number of coils per meter of array length. If None, estimates
        as 4/wavelength (4 coils per Halbach wavelength).

    Returns
    -------
    float
        Peak surface field B_surface (T)
    """
    if coils_per_meter is None:
        # Default: 4 coils per Halbach wavelength (90 deg segments)
        coils_per_meter = 4.0 / wavelength

    # Sheet current density (A/m)
    K = coils_per_meter * n_turns * i_current

    # Halbach enhancement factor (4-segment array)
    halbach_factor = math.sqrt(2)

    # Surface field from a current sheet
    # For an infinite current sheet, B = mu_0 * K / 2 on each side.
    # The Halbach arrangement concentrates it to one side, roughly doubling
    # the one-side field, giving B ≈ mu_0 * K on the front side.
    # With the sqrt(2) Halbach factor for a 4-magnet pattern:
    B_surface = MU0 * K * halbach_factor / 2

    # Note: this can exceed physical limits for very high K. In practice,
    # the field is limited by the ferromagnetic plate saturation (1.2 T)
    # for the EML bearing, or by the HTS critical current for both systems.
    return B_surface


def halbach_field_at_gap(B_surface, gap, wavelength):
    """
    Halbach array field at distance 'gap' from the surface.

    The field decays exponentially:
      B(z) = B_surface * exp(-2*pi*z / lambda)

    This is the fundamental property that makes Halbach arrays useful
    for the safety backstop: short wavelength means rapid decay, so
    the field is strong close up but negligible at the nominal gap.

    Parameters
    ----------
    B_surface : float
        Surface field (T)
    gap : float
        Distance from array surface (m)
    wavelength : float
        Halbach wavelength (m)

    Returns
    -------
    float
        Field at the gap (T)
    """
    k = 2 * math.pi / wavelength
    return B_surface * math.exp(-k * gap)


def halbach_field_gradient(B_surface, gap, wavelength):
    """
    dB/dz for a Halbach array at distance z = gap.

    dB/dz = -k * B_surface * exp(-k * z) = -k * B(z)

    Returns magnitude (positive value).
    """
    k = 2 * math.pi / wavelength
    B_at_gap = halbach_field_at_gap(B_surface, gap, wavelength)
    return k * B_at_gap


# =============================================================================
# SOLENOID MAGNETIC FIELD (for comparison)
# =============================================================================

def solenoid_surface_field(n_turns, i_current, coil_width, coils_per_meter=None):
    """
    Approximate field at the surface of a solenoid coil array.

    For a solenoid (non-Halbach) arrangement of the same coils, the
    surface field is determined by the same sheet current density but
    without the Halbach concentration. The field is split equally
    between both sides of the array.

      B_surface = mu_0 * K / 2

    where K = coils_per_meter * n_turns * i_current.

    This is lower than the Halbach case by the factor sqrt(2).

    The field profile with distance is NOT exponential. For a solenoid,
    the field falls off more slowly (roughly as 1/z for an extended array).
    """
    if coils_per_meter is None:
        coils_per_meter = 40  # default matching EML config

    K = coils_per_meter * n_turns * i_current
    return MU0 * K / 2


def solenoid_field_at_gap(B_surface, gap, coil_width):
    """
    Approximate solenoid field at distance 'gap' from the surface.

    For a wide, flat solenoid coil (width >> gap), the field is roughly
    constant near the surface and falls off approximately as:

      B(z) = B_surface / (1 + (z / w)^2)^(3/2)

    for a single dipole, but for an extended array of coils, the decay
    is slower. We use a semi-empirical model that captures the key
    difference from Halbach: slower decay with distance.

    For a continuous array of solenoid coils with width w, the field
    at distance z can be approximated as:

      B(z) = B_surface * w / sqrt(w^2 + (pi*z)^2)

    This gives B_surface at z=0 and falls off as ~1/z for z >> w.
    """
    return B_surface * coil_width / math.sqrt(coil_width**2 + (math.pi * gap)**2)


# =============================================================================
# EML ATTRACTIVE FORCE
# =============================================================================

def magnetic_pressure(B_gap):
    """
    Maxwell stress tensor: magnetic pressure between a field source
    and a ferromagnetic plate.

      P = B^2 / (2 * mu_0)    [Pa]

    This is the attractive pressure pulling the plate toward the coil.
    """
    return B_gap**2 / (2 * MU0)


def eml_attractive_force_halbach(n_turns, i_current, coil_width, wavelength,
                                  gap, plate_width, n_tracks=2,
                                  coils_per_meter=None):
    """
    Total attractive force per meter of ring from the EML Halbach array
    levitation bearing.

    Parameters
    ----------
    n_turns, i_current, coil_width, wavelength : see halbach_surface_field
    gap : float
        Air gap (m)
    plate_width : float
        Ferromagnetic plate width (m)
    n_tracks : int
        Number of levitation tracks (default 2)
    coils_per_meter : float or None
        Coils per meter of array. See halbach_surface_field.

    Returns
    -------
    float
        Force per meter of ring (N/m), attractive (positive = pulling casing toward cable)
    """
    B_surface = halbach_surface_field(n_turns, i_current, coil_width, wavelength,
                                       coils_per_meter=coils_per_meter)
    B_gap = halbach_field_at_gap(B_surface, gap, wavelength)

    # Limit to saturation
    if B_gap > 1.2:  # Sintered iron saturation
        B_gap = 1.2

    P = magnetic_pressure(B_gap)
    # Force per meter = pressure * plate_width * n_tracks
    return P * plate_width * n_tracks


def eml_attractive_force_solenoid(n_turns, i_current, coil_width,
                                   gap, plate_width, n_tracks=2,
                                   coils_per_meter=None):
    """
    Total attractive force per meter from a solenoid-type EML bearing
    (non-Halbach arrangement). For comparison.
    """
    B_surface = solenoid_surface_field(n_turns, i_current, coil_width,
                                        coils_per_meter=coils_per_meter)
    B_gap = solenoid_field_at_gap(B_surface, gap, coil_width)

    if B_gap > 1.2:
        B_gap = 1.2

    P = magnetic_pressure(B_gap)
    return P * plate_width * n_tracks


def eml_force_from_chapter5(i_current, i_c, gap, nominal_gap=0.100,
                             n_tracks=2, plate_width=1.4):
    """
    EML force using the simplified model from Chapter 5.

    The Chapter 5 analysis derived the operating point:
      - At I = 400 A (Ic/2), gap = 100 mm: B_gap = 0.523 T per track
      - Pressure = 109 kPa per track
      - Force = 218 kN/m total (2 tracks)
      - Maximum at Ic = 800 A: 4x force = 871 kN/m

    The force scales as I^2 (through B^2) and approximately as
    (nominal_gap / gap)^3 for small deviations from nominal gap
    (inverse cube law for a magnetic circuit with air gap dominating).

    Parameters
    ----------
    i_current : float
        Operating current (A)
    i_c : float
        Critical current (A) for Ic limiting
    gap : float
        Air gap (m)
    nominal_gap : float
        Nominal gap (m), default 0.100
    n_tracks : int
        Number of tracks
    plate_width : float
        Plate width (m)

    Returns
    -------
    float
        Force per meter (N/m), attractive
    """
    if i_current > i_c:
        i_current = i_c

    # Operating point from Chapter 5
    B_nominal = 0.523  # T per track at 400 A, 100 mm gap
    I_nominal = 400.0

    # Scale field with current
    B_at_current = B_nominal * (i_current / I_nominal)

    # Scale with gap: inverse-cube approximation for attractive EML
    # F ~ 1/gap^2 for magnetic pressure (since B ~ 1/gap for a magnetic circuit)
    # More precisely, B ~ 1/gap for a simple magnetic circuit, so P ~ 1/gap^2
    gap_ratio = nominal_gap / max(gap, 0.001)
    B_at_gap = B_at_current * gap_ratio

    # Enforce saturation limit
    B_sat = 1.2
    if B_at_gap > B_sat:
        B_at_gap = B_sat

    P = magnetic_pressure(B_at_gap)
    return P * plate_width * n_tracks


# =============================================================================
# SAFETY BACKSTOP: REPULSIVE FORCE FROM HALBACH ARRAY
# =============================================================================

def safety_repulsive_force(B_surface, gap_to_plate, wavelength,
                            plate_conductivity, plate_thickness,
                            plate_width, cable_velocity_lateral=0,
                            n_arrays=2):
    """
    Repulsive force from a superconducting Halbach safety array
    acting on a conductive reaction plate.

    There are two contributions to the repulsive force:

    1. IMAGE CURRENT REPULSION (quasi-static):
       When a conducting plate is near a Halbach array, the time-varying
       field component (from cable vibration or gap oscillation) induces
       eddy currents. For a superconducting Halbach array with persistent
       current, the field is DC, and eddy currents are only induced when
       the gap is CHANGING (i.e., during a safety event).

       However, there is a subtlety: the cable is moving along the ring
       at ~8,600 m/s relative to the casing. If the safety Halbach array
       has any spatial variation in the direction of travel (it does: the
       wavelength lambda creates a periodic field pattern), then the cable
       motion through this pattern induces eddy currents continuously.
       This is the mechanism that provides the repulsive force.

       The repulsive force per unit area from a Halbach array on a
       conducting plate moving at velocity v through the periodic field is:

         P_repulsive = (sigma * delta * v * k * B_surface^2 * exp(-2k*z)) /
                       (1 + (sigma * delta * mu_0 * v * k)^2)

       where sigma is conductivity, delta is plate thickness, v is the
       velocity component along the array, k = 2*pi/lambda, z is gap.

       At high velocity (sigma*delta*mu_0*v*k >> 1), this saturates to:

         P_repulsive = B_surface^2 * exp(-2k*z) / (2 * mu_0)

       which is the full image current limit (the plate acts like a
       perfect conductor and the repulsive pressure equals the full
       Maxwell stress).

    2. At the orbital ring's cable velocity of ~8,600 m/s and with a
       short Halbach wavelength of 100 mm, the skin depth parameter
       sigma*delta*mu_0*v*k is enormous, so we are firmly in the
       image-current-limited regime. The force is essentially:

         P = B^2(z) / (2*mu_0) per array

    Parameters
    ----------
    B_surface : float
        Halbach surface field (T)
    gap_to_plate : float
        Distance from array face to reaction plate (m)
    wavelength : float
        Halbach wavelength (m)
    plate_conductivity : float
        Plate conductivity (S/m)
    plate_thickness : float
        Plate thickness (m)
    plate_width : float
        Plate width (m)
    cable_velocity_lateral : float
        Cable velocity along the array direction (m/s).
        For the orbital ring, this is ~8,600 m/s.
    n_arrays : int
        Number of safety arrays on this side (top or bottom)

    Returns
    -------
    float
        Repulsive force per meter of ring (N/m)
    """
    k = 2 * math.pi / wavelength
    B_at_gap = B_surface * math.exp(-k * gap_to_plate)

    # Check if we are in the image-current-limited regime
    sigma = plate_conductivity
    delta = plate_thickness
    v = abs(cable_velocity_lateral)

    if v > 0:
        skin_param = sigma * delta * MU0 * v * k
    else:
        skin_param = 0

    if skin_param > 10:
        # Image current limit: full repulsion
        P_repulsive = B_at_gap**2 / (2 * MU0)
    elif skin_param > 0:
        # Intermediate regime
        P_repulsive = (B_at_gap**2 / (2 * MU0)) * skin_param**2 / (1 + skin_param**2)
    else:
        # No motion: no eddy current repulsion from DC field
        # But during a safety event, the gap is changing, so dB/dt != 0
        # We model this separately in the transient simulation
        P_repulsive = 0

    return P_repulsive * plate_width * n_arrays


def safety_drag_force(B_surface, gap_to_plate, wavelength,
                       plate_conductivity, plate_thickness,
                       plate_width, cable_velocity_lateral,
                       n_arrays=2):
    """
    Drag force from eddy currents in the safety reaction plate.

    This is the parasitic drag that exists even when the safety is
    NOT engaged (i.e., at the nominal gap). The Halbach wavelength
    is chosen short enough that the field is negligible at the
    nominal gap, making this drag very small.

    The drag force is:
      F_drag = P_repulsive * (1/skin_param) for skin_param >> 1

    Actually, for a Halbach array and a plate in the image current limit:
      Drag/Lift = 1/(sigma * delta * mu_0 * v * k)

    This ratio is extremely small at orbital velocities.
    """
    k = 2 * math.pi / wavelength
    B_at_gap = B_surface * math.exp(-k * gap_to_plate)
    sigma = plate_conductivity
    delta = plate_thickness
    v = abs(cable_velocity_lateral)

    if v < 1:
        return 0

    skin_param = sigma * delta * MU0 * v * k
    P_repulsive = B_at_gap**2 / (2 * MU0)

    if skin_param > 0.01:
        drag_to_lift = 1.0 / skin_param
    else:
        drag_to_lift = 1.0

    return P_repulsive * plate_width * n_arrays * drag_to_lift


# =============================================================================
# THERMAL ANALYSIS
# =============================================================================

def eddy_current_heating_rate(B_field, frequency, sigma, volume):
    """
    Eddy current power dissipation in a conducting plate.

    For a plate in a sinusoidal field of amplitude B and frequency f:
      P = (pi^2 * B^2 * f^2 * d^2 * sigma) / 6 * volume

    where d is the plate thickness (for thin plates where skin depth >> d).

    For the safety backstop, the "frequency" is determined by the cable
    velocity and the Halbach wavelength:
      f = v_cable / wavelength

    Parameters
    ----------
    B_field : float
        Field amplitude at the plate (T)
    frequency : float
        Frequency of field variation (Hz)
    sigma : float
        Electrical conductivity (S/m)
    volume : float
        Plate volume per meter of ring (m^3/m = m^2)

    Returns
    -------
    float
        Heating power per meter of ring (W/m)
    """
    # For thin plate (thickness << skin depth)
    # P/V = (pi * f * B)^2 * sigma / 6  ... no, let me use the standard formula
    #
    # Actually, the standard eddy current loss formula for a plate of thickness d
    # in a uniform alternating field B*sin(2*pi*f*t) is:
    #   P/V = (pi * B * f * d)^2 / (6 * rho)
    # where rho is resistivity.
    #
    # But this is for a plate much thinner than the skin depth. If the plate
    # is thicker than the skin depth, only the skin layer heats.
    #
    # We'll compute both cases and use the appropriate one.

    rho = 1.0 / sigma if sigma > 0 else 1e10
    omega = 2 * math.pi * frequency

    # Skin depth
    if frequency > 0 and sigma > 0:
        skin_depth = math.sqrt(2 * rho / (omega * MU0))
    else:
        return 0

    # The heating per unit volume in the skin layer is approximately:
    # P/V_skin = sigma * (omega * B)^2 / 2 in the skin layer
    # Total power per unit area = sigma * (omega * B)^2 * skin_depth / 2
    #
    # More precisely, for a semi-infinite conductor:
    # P/area = (omega * B)^2 / (2 * sigma * skin_depth_factor)
    # = B^2 * omega / (2 * mu_0 * skin_depth)  [approximate]
    #
    # Standard result: P per unit area = B^2 * sqrt(omega/(2*mu_0*sigma)) * omega/(2*mu_0)
    # Let's use: P/A = B^2 * omega * skin_depth * sigma / 4
    # which is equivalent to B^2 * omega / (4 * mu_0 / skin_depth)

    # Cleaner: Power dissipated per unit area of a conducting half-space
    # in an AC field of amplitude B_0 at the surface:
    #   P/A = B_0^2 * omega / (4 * mu_0 * sigma * skin_depth^(-1))
    #   Actually: P/A = B_0^2 / (2 * mu_0) * sqrt(omega * rho / (2 * mu_0))
    #   = B_0^2 / (2 * mu_0) * (1 / skin_depth) * (rho / mu_0)  ... getting circular.
    #
    # Standard textbook result (Griffiths, Jackson):
    #   P/A = B_0^2 * omega * delta / (4 * mu_0)
    # where delta is the skin depth = sqrt(2*rho/(omega*mu_0))

    P_per_area = B_field**2 * omega * skin_depth / (4 * MU0)

    return P_per_area * volume  # volume here is width * n_arrays (area per meter)


def eddy_current_heating_halbach(B_surface, gap, wavelength, v_cable,
                                  mat_props, plate_thickness, plate_width,
                                  n_arrays=2):
    """
    Eddy current heating in the safety reaction plate from the cable
    moving through the periodic Halbach field.

    The frequency is f = v_cable / wavelength.
    The field at the plate is B(gap) = B_surface * exp(-k * gap).
    """
    k = 2 * math.pi / wavelength
    B_at_plate = B_surface * math.exp(-k * gap)
    f = abs(v_cable) / wavelength
    sigma = 1.0 / mat_props["resistivity"]

    # For the Halbach case, the field penetrates the plate.
    # Use the standard formula for power per unit area.
    omega = 2 * math.pi * f
    rho = mat_props["resistivity"]

    if omega < 1e-6:
        return 0

    skin_depth = math.sqrt(2 * rho / (omega * MU0))

    # If plate is thinner than skin depth, use thin-plate formula
    if plate_thickness < skin_depth:
        # Thin plate: P/A = (pi * B * f)^2 * plate_thickness^2 / (6 * rho)
        P_per_area = (math.pi * B_at_plate * f)**2 * plate_thickness**2 / (6 * rho)
    else:
        # Thick plate: P/A = B^2 * omega * skin_depth / (4 * mu_0)
        P_per_area = B_at_plate**2 * omega * skin_depth / (4 * MU0)

    return P_per_area * plate_width * n_arrays


def radiative_cooling_rate(T_plate, emissivity, plate_width, n_arrays=2,
                           T_environment=None):
    """
    Radiative cooling power per meter of ring from the safety reaction plate.

    The plate radiates to the casing interior. In steady state, the casing
    wall is at ~247 K. The plate radiates from its exposed surface area.

    P/length = epsilon * sigma_SB * (T_plate^4 - T_env^4) * width * n_arrays

    Parameters
    ----------
    T_plate : float
        Plate temperature (K)
    emissivity : float
        Plate surface emissivity
    plate_width : float
        Width of exposed surface (m)
    n_arrays : int
        Number of plates radiating
    T_environment : float or None
        Surrounding temperature (K). If None, uses T_CASING_WALL.

    Returns
    -------
    float
        Cooling power per meter (W/m)
    """
    if T_environment is None:
        T_environment = T_CASING_WALL

    return (emissivity * STEFAN_BOLTZMANN *
            (T_plate**4 - T_environment**4) *
            plate_width * n_arrays)


def equilibrium_temperature(heating_rate, emissivity, plate_width, n_arrays=2,
                            T_environment=None):
    """
    Find the equilibrium temperature where radiative cooling equals heating.

    P_heat = epsilon * sigma_SB * (T^4 - T_env^4) * A

    Solve for T:
      T = (P_heat / (epsilon * sigma_SB * A) + T_env^4)^(1/4)
    """
    if T_environment is None:
        T_environment = T_CASING_WALL

    A = plate_width * n_arrays
    if A < 1e-10 or emissivity < 1e-10:
        return float('inf')

    T4 = heating_rate / (emissivity * STEFAN_BOLTZMANN * A) + T_environment**4
    if T4 < 0:
        return T_environment
    return T4**0.25


# =============================================================================
# MATERIAL PROPERTY LOOKUP
# =============================================================================

def get_material(name):
    """Get material properties dictionary by name."""
    if name in REACTION_PLATE_MATERIALS:
        return REACTION_PLATE_MATERIALS[name]
    raise ValueError(f"Unknown material: {name}. Available: {list(REACTION_PLATE_MATERIALS.keys())}")


def resistivity_at_temp(mat_props, T):
    """
    Resistivity at temperature T, accounting for temperature coefficient.

    rho(T) = rho_293 * (1 + alpha * (T - 293))

    For materials with negative temp coefficient (CNT), resistivity decreases
    with temperature.
    """
    rho_293 = mat_props["resistivity"]
    alpha = mat_props["temp_coeff_rho"]
    return rho_293 * (1 + alpha * (T - 293))


def conductivity_at_temp(mat_props, T):
    """Conductivity at temperature T."""
    rho = resistivity_at_temp(mat_props, T)
    return 1.0 / rho if rho > 0 else 0
