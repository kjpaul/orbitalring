#!/usr/bin/env python3
"""
Levitation System Simulation

Compares Halbach array vs solenoid bearing designs, models the EML
attractive bearing and the passive superconducting Halbach safety
backstop, and performs thermal analysis of the safety reaction plate.

Usage:
    py levitation_simulation.py                    # Run all analyses, show all graphs
    py levitation_simulation.py force_comparison   # Halbach vs solenoid force profile
    py levitation_simulation.py safety_crossover   # Safety vs EML force crossover
    py levitation_simulation.py thermal            # Thermal analysis of safety plate
    py levitation_simulation.py materials          # Material trade study
    py levitation_simulation.py drag               # Parasitic drag at nominal gap
    py levitation_simulation.py all                # All graphs
    py levitation_simulation.py --material=copper_OFHC   # Set safety plate material
    py levitation_simulation.py --wavelength=0.05        # Set Halbach wavelength (m)

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import math
import os

# Try to import matplotlib; if not available, skip plotting
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not found. Graphs will not be generated.")

from levitation_config import *
from levitation_physics import *


# =============================================================================
# COMMAND LINE PARSING
# =============================================================================

def parse_args():
    """Parse command line arguments."""
    args = {
        "graphs": [],
        "material": SAFETY_PLATE_MATERIAL,
        "wavelength": SAFETY_WAVELENGTH,
        "show_all": False,
    }

    for arg in sys.argv[1:]:
        if arg.startswith("--material="):
            args["material"] = arg.split("=", 1)[1]
        elif arg.startswith("--wavelength="):
            args["wavelength"] = float(arg.split("=", 1)[1])
        elif arg == "all":
            args["show_all"] = True
        else:
            args["graphs"].append(arg)

    if not args["graphs"] and not args["show_all"]:
        args["show_all"] = True

    return args


# =============================================================================
# ANALYSIS 1: HALBACH vs SOLENOID FORCE COMPARISON
# =============================================================================

def run_force_comparison():
    """
    Compare the force-vs-gap profile of a Halbach array EML bearing
    vs a solenoid EML bearing, both using the same HTS tape and current.

    The key difference:
    - Halbach concentrates field on one side, giving ~40% more surface field
    - Halbach field decays exponentially with characteristic length lambda/(2*pi)
    - Solenoid field decays more slowly (roughly as 1/z)

    For the EML bearing (large wavelength, ~50 m pitch), both give similar
    results at the 100 mm operating gap. The Halbach advantage is moderate
    for the EML bearing but dramatic for the short-wavelength safety arrays.
    """
    print("\n" + "="*72)
    print("ANALYSIS 1: HALBACH vs SOLENOID EML BEARING COMPARISON")
    print("="*72)

    # EML bearing parameters (from Chapter 5 / Section 10.3.1)
    n_turns = EML_TURNS_PER_COIL
    i_op = EML_I_OPERATING
    i_max = EML_HTS_IC
    coil_width = EML_COIL_WIDTH
    plate_width = EML_PLATE_WIDTH
    n_tracks = EML_N_TRACKS

    # The EML Halbach array uses long coils (~50 m) arranged end-to-end
    # with staggered joints. The field is quasi-uniform in the direction
    # of cable travel (by design, to minimize losses). The relevant field
    # decay with gap distance is set by the TRANSVERSE coil geometry.
    #
    # For a wide, flat array (width >> gap), the field at the gap is
    # determined by the magnetic circuit: the coil produces a field,
    # some of which crosses the air gap to reach the ferromagnetic plate.
    # The effective "wavelength" for the exponential decay model is
    # determined by the transverse coil width: lambda_eff ~ pi * width.
    # For the 1 m wide EML array: lambda_eff ~ 3.14 m.
    #
    # This gives a decay length of lambda/(2*pi) = width/2 = 0.5 m,
    # meaning at 100 mm gap the field is exp(-0.1/0.5) = 0.82 of surface.
    eml_wavelength = math.pi * coil_width  # ~3.14 m for 1 m wide coils
    eml_coils_per_m = EML_COILS_PER_METER  # 40 coils/m

    # Generate gap sweep
    gaps = [GAP_SWEEP_MIN + i * (GAP_SWEEP_MAX - GAP_SWEEP_MIN) / GAP_SWEEP_POINTS
            for i in range(GAP_SWEEP_POINTS + 1)]

    # Compute forces at operating current
    f_halbach_op = []
    f_solenoid_op = []
    f_ch5_op = []

    # Compute forces at maximum current
    f_halbach_max = []
    f_solenoid_max = []
    f_ch5_max = []

    for g in gaps:
        # Halbach EML at operating current
        fh = eml_attractive_force_halbach(n_turns, i_op, coil_width, eml_wavelength,
                                           g, plate_width, n_tracks,
                                           coils_per_meter=eml_coils_per_m)
        f_halbach_op.append(fh / 1000)  # kN/m

        # Solenoid EML at operating current
        fs = eml_attractive_force_solenoid(n_turns, i_op, coil_width,
                                            g, plate_width, n_tracks,
                                            coils_per_meter=eml_coils_per_m)
        f_solenoid_op.append(fs / 1000)

        # Chapter 5 model at operating current
        fc = eml_force_from_chapter5(i_op, i_max, g, EML_NOMINAL_GAP,
                                      n_tracks, plate_width)
        f_ch5_op.append(fc / 1000)

        # At maximum current
        fh_max = eml_attractive_force_halbach(n_turns, i_max, coil_width, eml_wavelength,
                                               g, plate_width, n_tracks,
                                               coils_per_meter=eml_coils_per_m)
        f_halbach_max.append(fh_max / 1000)

        fs_max = eml_attractive_force_solenoid(n_turns, i_max, coil_width,
                                                g, plate_width, n_tracks,
                                                coils_per_meter=eml_coils_per_m)
        f_solenoid_max.append(fs_max / 1000)

        fc_max = eml_force_from_chapter5(i_max, i_max, g, EML_NOMINAL_GAP,
                                          n_tracks, plate_width)
        f_ch5_max.append(fc_max / 1000)

    # Print key values
    nominal_load = M_LOAD_M * G_NET / 1000  # kN/m
    print(f"\nNominal casing load: {nominal_load:.1f} kN/m")
    print(f"  ({M_LOAD_M:,} kg/m × {G_NET} m/s²)")

    # Find forces at nominal gap
    idx_nom = min(range(len(gaps)), key=lambda i: abs(gaps[i] - EML_NOMINAL_GAP))
    print(f"\nAt nominal gap ({EML_NOMINAL_GAP*1000:.0f} mm), operating current ({i_op} A):")
    print(f"  Halbach EML:   {f_halbach_op[idx_nom]:.1f} kN/m")
    print(f"  Solenoid EML:  {f_solenoid_op[idx_nom]:.1f} kN/m")
    print(f"  Chapter 5 model: {f_ch5_op[idx_nom]:.1f} kN/m")
    print(f"  Halbach advantage: {f_halbach_op[idx_nom]/f_solenoid_op[idx_nom]:.2f}x")

    print(f"\nAt nominal gap, maximum current ({i_max} A):")
    print(f"  Halbach EML:   {f_halbach_max[idx_nom]:.1f} kN/m")
    print(f"  Solenoid EML:  {f_solenoid_max[idx_nom]:.1f} kN/m")
    print(f"  Chapter 5 model: {f_ch5_max[idx_nom]:.1f} kN/m")

    # Margin
    margin_halbach = f_halbach_max[idx_nom] - nominal_load
    margin_solenoid = f_solenoid_max[idx_nom] - nominal_load
    print(f"\nMargin above nominal load (at Ic):")
    print(f"  Halbach: {margin_halbach:.1f} kN/m  ({margin_halbach/nominal_load*100:.0f}% of load)")
    print(f"  Solenoid: {margin_solenoid:.1f} kN/m  ({margin_solenoid/nominal_load*100:.0f}% of load)")

    # Field comparison
    B_h_surface = halbach_surface_field(n_turns, i_op, coil_width, eml_wavelength,
                                         coils_per_meter=eml_coils_per_m)
    B_s_surface = solenoid_surface_field(n_turns, i_op, coil_width,
                                          coils_per_meter=eml_coils_per_m)
    B_h_gap = halbach_field_at_gap(B_h_surface, EML_NOMINAL_GAP, eml_wavelength)
    B_s_gap = solenoid_field_at_gap(B_s_surface, EML_NOMINAL_GAP, coil_width)

    print(f"\nField comparison at {i_op} A:")
    print(f"  Halbach surface field: {B_h_surface:.3f} T")
    print(f"  Solenoid surface field: {B_s_surface:.3f} T")
    print(f"  Halbach field at {EML_NOMINAL_GAP*1000:.0f} mm gap: {B_h_gap:.3f} T")
    print(f"  Solenoid field at {EML_NOMINAL_GAP*1000:.0f} mm gap: {B_s_gap:.3f} T")

    return {
        "gaps_mm": [g * 1000 for g in gaps],
        "f_halbach_op": f_halbach_op,
        "f_solenoid_op": f_solenoid_op,
        "f_ch5_op": f_ch5_op,
        "f_halbach_max": f_halbach_max,
        "f_solenoid_max": f_solenoid_max,
        "f_ch5_max": f_ch5_max,
        "nominal_load": nominal_load,
    }


# =============================================================================
# ANALYSIS 2: SAFETY vs EML FORCE CROSSOVER
# =============================================================================

def run_safety_crossover(material_name, wavelength):
    """
    Find the gap at which the safety backstop repulsive force exceeds
    the EML attractive force. This is the critical design requirement:
    the safety must win at close range to prevent cable-casing contact.

    The EML attractive force increases as gap closes (inverse square).
    The safety repulsive force also increases as gap closes (exponential).
    Because the safety uses a short Halbach wavelength, its force grows
    FASTER than the EML force at close range.
    """
    print("\n" + "="*72)
    print("ANALYSIS 2: SAFETY BACKSTOP vs EML FORCE CROSSOVER")
    print("="*72)

    mat = get_material(material_name)
    print(f"\nSafety reaction plate material: {mat['name']}")
    print(f"  Resistivity: {mat['resistivity']:.2e} ohm-m")
    print(f"  Density: {mat['density']} kg/m^3")
    print(f"  Melting point: {mat['melting_point']} K ({mat['melting_point']-273:.0f} C)")
    print(f"  Max service temp: {mat['max_service_temp']} K ({mat['max_service_temp']-273:.0f} C)")
    print(f"\nSafety Halbach wavelength: {wavelength*1000:.0f} mm")

    # Safety array parameters
    # The safety Halbach arrays have 4 coils per wavelength (one per 90 deg
    # segment). With wavelength = 100 mm, coil pitch = 25 mm, so 40 coils/m.
    # Each coil has fewer turns than the EML coils (smaller cross-section),
    # but operates at higher persistent current.
    safety_coils_per_m = 4.0 / wavelength  # 4 coils per Halbach period

    # Safety coils: 50 mm thick assembly, 25 mm coil pitch, 3 mm tape
    # Can fit 8 turns/layer x 3 layers = 24 turns per coil
    safety_turns = SAFETY_TURNS_PER_COIL

    B_safety_surface = halbach_surface_field(
        safety_turns,
        SAFETY_I_PERSISTENT,
        SAFETY_COIL_WIDTH,
        wavelength,
        coils_per_meter=safety_coils_per_m
    )
    print(f"Safety Halbach surface field: {B_safety_surface:.3f} T")
    print(f"  at persistent current {SAFETY_I_PERSISTENT} A")
    print(f"  {safety_turns} turns/coil, {safety_coils_per_m:.0f} coils/m")

    # Decay length
    decay_length = wavelength / (2 * math.pi)
    print(f"  Decay length (lambda/2pi): {decay_length*1000:.1f} mm")
    print(f"  Field at 20 mm: {halbach_field_at_gap(B_safety_surface, 0.020, wavelength):.4f} T")
    print(f"  Field at 50 mm: {halbach_field_at_gap(B_safety_surface, 0.050, wavelength):.6f} T")
    print(f"  Field at 100 mm (nominal gap): {halbach_field_at_gap(B_safety_surface, 0.100, wavelength):.2e} T")

    # Cable velocity (for eddy current calculation)
    v_cable = 8620  # m/s (cable velocity relative to casing)

    sigma = 1.0 / mat["resistivity"]

    # ---------------------------------------------------------------
    # FAILURE MODE 1: Total EML power failure
    # ---------------------------------------------------------------
    # The EML loses power. The cable is in orbit (weightless). The casing
    # falls at g_net relative to the cable. From the casing's perspective,
    # the cable rises toward the casing ceiling. The top safety arrays
    # engage. The repulsive force acts on BOTH cable (pushed down) and
    # casing (pushed up, opposing fall).
    #
    # The safety must arrest the casing's fall before the cable contacts
    # the casing ceiling. Required: F_safety(gap) > m_load * g_net at
    # some gap, and the integrated impulse must be enough to stop the
    # casing within the available travel distance (~80 mm from forbidden
    # zone boundary to contact).
    #
    # ---------------------------------------------------------------
    # FAILURE MODE 2: EML overcurrent (gap control malfunction)
    # ---------------------------------------------------------------
    # The EML current surges to Ic, pulling the cable toward the casing
    # bottom with more force than gravity requires. The excess force
    # accelerates the cable downward. The bottom safety engages.
    # Required: F_safety > F_EML(gap) - m_load*g_net (the excess pull).
    # ---------------------------------------------------------------

    nominal_load = M_LOAD_M * G_NET  # N/m

    # Sweep gap from 1 mm to 100 mm (distance from safety array to plate)
    gaps = [0.001 + i * 0.099 / GAP_SWEEP_POINTS for i in range(GAP_SWEEP_POINTS + 1)]

    f_gravity = []    # Gravitational load (constant, what safety must overcome in Mode 1)
    f_safety = []     # Safety repulsive force (kN/m)
    f_eml_excess = [] # Excess EML force in Mode 2 (EML at Ic minus gravity)
    f_net_mode1 = []  # Net force, Mode 1 (safety vs gravity)
    f_net_mode2 = []  # Net force, Mode 2 (safety vs excess EML)

    for g in gaps:
        # Safety repulsive force at gap g (2 arrays on the engaged side)
        fs = safety_repulsive_force(
            B_safety_surface, g, wavelength,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=2
        )
        f_safety.append(fs / 1000)
        f_gravity.append(nominal_load / 1000)

        # Mode 1: safety vs gravity
        f_net_mode1.append((fs - nominal_load) / 1000)

        # Mode 2: EML at Ic with gap closing
        # As cable moves toward bottom wall, EML gap to bottom decreases.
        # If safety gap = g, the cable has moved (nominal_safety_gap - g)
        # from its nominal position. The EML gap has decreased by the same amount.
        # EML gap = nominal_EML_gap - (nominal_safety_gap - g)
        # For simplicity, assume nominal safety gap ≈ nominal EML gap = 100 mm
        eml_gap = g  # approximately: when safety gap = g, EML gap ≈ g too
        if eml_gap < 0.002:
            eml_gap = 0.002
        fe_at_ic = eml_force_from_chapter5(EML_HTS_IC, EML_HTS_IC, eml_gap,
                                            EML_NOMINAL_GAP, EML_N_TRACKS, EML_PLATE_WIDTH)
        excess = fe_at_ic - nominal_load
        f_eml_excess.append(max(excess, 0) / 1000)

        # Mode 2 net: safety vs excess EML
        f_net_mode2.append((fs - max(excess, 0)) / 1000)

    # Find crossover points
    crossover_mode1 = None
    for i in range(len(gaps) - 1):
        if f_net_mode1[i+1] >= 0 and f_net_mode1[i] < 0:
            frac = -f_net_mode1[i] / (f_net_mode1[i+1] - f_net_mode1[i])
            crossover_mode1 = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break

    crossover_mode2 = None
    for i in range(len(gaps) - 1):
        if f_net_mode2[i+1] >= 0 and f_net_mode2[i] < 0:
            frac = -f_net_mode2[i] / (f_net_mode2[i+1] - f_net_mode2[i])
            crossover_mode2 = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break

    print(f"\n--- FAILURE MODE 1: Total EML Power Loss ---")
    print(f"  Casing weight to arrest: {nominal_load/1000:.1f} kN/m")
    if crossover_mode1 is not None:
        print(f"  Safety exceeds gravity at gap = {crossover_mode1*1000:.1f} mm")
    elif f_net_mode1[0] > 0:
        print(f"  Safety exceeds gravity at ALL gaps (strong backstop)")
    else:
        print(f"  WARNING: Safety never exceeds gravity in sweep range!")

    print(f"\n--- FAILURE MODE 2: EML Overcurrent (stuck at Ic) ---")
    if crossover_mode2 is not None:
        print(f"  Safety exceeds excess EML at gap = {crossover_mode2*1000:.1f} mm")
    elif f_net_mode2[0] > 0:
        print(f"  Safety exceeds excess EML at ALL gaps")
    else:
        print(f"  WARNING: Safety never overcomes excess EML force!")

    # Print force at key gaps
    for g_check in [0.002, 0.005, 0.010, 0.020, 0.050, 0.080, 0.100]:
        idx = min(range(len(gaps)), key=lambda i: abs(gaps[i] - g_check))
        print(f"\n  At safety gap = {g_check*1000:.0f} mm:")
        print(f"    Safety repulsive:    {f_safety[idx]:>10.1f} kN/m")
        print(f"    Gravity load:        {f_gravity[idx]:>10.1f} kN/m")
        print(f"    Excess EML (at Ic):  {f_eml_excess[idx]:>10.1f} kN/m")
        print(f"    Net Mode 1 (+ safe): {f_net_mode1[idx]:>10.1f} kN/m")
        print(f"    Net Mode 2 (+ safe): {f_net_mode2[idx]:>10.1f} kN/m")

    # Transient analysis: can the safety arrest the casing in Mode 1?
    print(f"\n--- Transient: Can safety arrest casing fall (Mode 1)? ---")
    # When EML fails, the cable stays in orbit (weightless). The casing
    # falls at g_net relative to the cable. From the cable's reference
    # frame, the casing falls downward, and the cable appears to rise
    # toward the casing ceiling. The top safety arrays see a closing gap.
    #
    # Convention: v_fall = positive = gap closing (casing falling)
    # F_safety opposes the fall (pushes casing back up).
    #
    # Equation of motion:
    #   m * dv/dt = m * g_net - F_safety(gap)
    # where m = M_LOAD_M (per meter), g_net = 9.038 m/s^2
    #
    # Starting condition: gap = 80 mm (just entering forbidden zone),
    # with some velocity from free-fall over the 20 mm between nominal
    # (100 mm) and forbidden zone boundary (80 mm).
    # v_initial = sqrt(2 * g_net * 0.020) = sqrt(0.361) = 0.601 m/s

    dt = 0.00001  # 10 microseconds
    gap_initial = 0.080
    v_initial = math.sqrt(2 * G_NET * (EML_NOMINAL_GAP - gap_initial))
    gap_current = gap_initial
    v_fall = v_initial
    t = 0
    min_gap = gap_current
    arrested = False
    max_t = 2.0

    # Track trajectory for printout
    checkpoints = {0.001: False, 0.005: False, 0.010: False, 0.050: False}

    while t < max_t and gap_current > 0.0005:
        # Safety force at current gap (per meter of ring)
        fs = safety_repulsive_force(
            B_safety_surface, max(gap_current, 0.001), wavelength,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable, n_arrays=2)

        # Net acceleration: gravity pulls casing down, safety pushes up
        a_net = G_NET - fs / M_LOAD_M

        # Update velocity and position
        v_fall += a_net * dt
        gap_current -= v_fall * dt

        if gap_current < min_gap:
            min_gap = gap_current

        # Check if arrested (velocity reversed while being decelerated)
        if v_fall <= 0:
            arrested = True
            break

        # Print checkpoints
        for ck in list(checkpoints.keys()):
            if not checkpoints[ck] and gap_current <= ck + 0.0005:
                checkpoints[ck] = True
                print(f"  t = {t*1000:.2f} ms: gap = {gap_current*1000:.1f} mm, "
                      f"v = {v_fall*1000:.1f} mm/s, F_safety = {fs/1000:.1f} kN/m")

        t += dt

    if arrested:
        print(f"\n  Casing fall ARRESTED at t = {t*1000:.2f} ms")
        print(f"  Minimum gap reached: {min_gap*1000:.2f} mm")
        print(f"  Entry velocity was: {v_initial*1000:.1f} mm/s")
        print(f"  Safety force at min gap: {safety_repulsive_force(B_safety_surface, max(min_gap,0.001), wavelength, sigma, SAFETY_PLATE_THICKNESS, SAFETY_PLATE_WIDTH, v_cable, 2)/1000:.1f} kN/m")
    else:
        print(f"\n  WARNING: Casing NOT arrested in {max_t} s")
        print(f"  Gap reached: {gap_current*1000:.2f} mm")
        print(f"  Velocity at end: {v_fall*1000:.1f} mm/s")

    return {
        "gaps_mm": [g * 1000 for g in gaps],
        "f_gravity": f_gravity,
        "f_safety": f_safety,
        "f_eml_excess": f_eml_excess,
        "f_net_mode1": f_net_mode1,
        "f_net_mode2": f_net_mode2,
        "crossover_mode1_mm": crossover_mode1 * 1000 if crossover_mode1 else None,
        "crossover_mode2_mm": crossover_mode2 * 1000 if crossover_mode2 else None,
        "material": mat,
    }


# =============================================================================
# ANALYSIS 3: THERMAL ANALYSIS OF SAFETY REACTION PLATE
# =============================================================================

def run_thermal_analysis(material_name, wavelength):
    """
    Determine at what gap the safety reaction plate reaches its thermal limit.

    The cable moves through the safety Halbach array's periodic field at
    ~8,600 m/s, inducing eddy currents in the reaction plate. These
    eddy currents provide the repulsive force but also heat the plate.
    The plate cools by radiation to the casing interior.

    At the nominal gap (100 mm), the Halbach field is negligible and
    heating is essentially zero. As the gap closes during a safety event,
    heating increases exponentially. We need to know:
    1. Equilibrium temperature at each gap
    2. Gap at which the plate reaches its service temperature limit
    3. Gap at which the plate reaches its melting point
    4. How long the plate can survive at various gaps (transient analysis)
    """
    print("\n" + "="*72)
    print("ANALYSIS 3: SAFETY REACTION PLATE THERMAL ANALYSIS")
    print("="*72)

    mat = get_material(material_name)
    print(f"\nMaterial: {mat['name']}")
    print(f"  Melting point: {mat['melting_point']} K ({mat['melting_point']-273:.0f} C)")
    print(f"  Max service temp: {mat['max_service_temp']} K ({mat['max_service_temp']-273:.0f} C)")
    print(f"  Emissivity: {mat['emissivity']}")

    safety_coils_per_m = 4.0 / wavelength
    safety_turns = SAFETY_TURNS_PER_COIL
    B_safety_surface = halbach_surface_field(
        safety_turns, SAFETY_I_PERSISTENT,
        SAFETY_COIL_WIDTH, wavelength,
        coils_per_meter=safety_coils_per_m
    )

    v_cable = 8620  # m/s

    # Sweep gap
    gaps = [0.001 + i * 0.199 / GAP_SWEEP_POINTS for i in range(GAP_SWEEP_POINTS + 1)]

    heating = []
    cooling_at_service = []
    T_eq = []

    for g in gaps:
        # Heating rate
        Q_heat = eddy_current_heating_halbach(
            B_safety_surface, g, wavelength, v_cable,
            mat, SAFETY_PLATE_THICKNESS, SAFETY_PLATE_WIDTH,
            n_arrays=2
        )
        heating.append(Q_heat)

        # Cooling at service temperature
        Q_cool = radiative_cooling_rate(
            mat["max_service_temp"], mat["emissivity"],
            SAFETY_PLATE_WIDTH, n_arrays=2
        )
        cooling_at_service.append(Q_cool)

        # Equilibrium temperature
        T = equilibrium_temperature(
            Q_heat, mat["emissivity"],
            SAFETY_PLATE_WIDTH, n_arrays=2
        )
        T_eq.append(T)

    # Find gap at which T_eq = service temp
    service_gap = None
    for i in range(len(gaps) - 1):
        if T_eq[i] > mat["max_service_temp"] and T_eq[i+1] <= mat["max_service_temp"]:
            frac = (T_eq[i] - mat["max_service_temp"]) / (T_eq[i] - T_eq[i+1])
            service_gap = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break
        elif T_eq[i] <= mat["max_service_temp"] and T_eq[i+1] > mat["max_service_temp"]:
            frac = (mat["max_service_temp"] - T_eq[i]) / (T_eq[i+1] - T_eq[i])
            service_gap = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break

    # Find gap at which T_eq = melting point
    melt_gap = None
    for i in range(len(gaps) - 1):
        if T_eq[i] > mat["melting_point"] and T_eq[i+1] <= mat["melting_point"]:
            frac = (T_eq[i] - mat["melting_point"]) / (T_eq[i] - T_eq[i+1])
            melt_gap = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break
        elif T_eq[i] <= mat["melting_point"] and T_eq[i+1] > mat["melting_point"]:
            frac = (mat["melting_point"] - T_eq[i]) / (T_eq[i+1] - T_eq[i])
            melt_gap = gaps[i] + frac * (gaps[i+1] - gaps[i])
            break

    print(f"\n--- Steady-State Thermal Limits ---")
    if service_gap:
        print(f"Service temp ({mat['max_service_temp']} K) reached at gap = {service_gap*1000:.1f} mm")
    else:
        # Check if always above or below
        if T_eq[0] < mat["max_service_temp"]:
            print(f"Equilibrium temperature stays below service limit at all gaps.")
        else:
            print(f"Equilibrium temperature exceeds service limit at all gaps in range.")

    if melt_gap:
        print(f"Melting point ({mat['melting_point']} K) reached at gap = {melt_gap*1000:.1f} mm")
    else:
        if T_eq[0] < mat["melting_point"]:
            print(f"Equilibrium temperature stays below melting point at all gaps.")
        else:
            print(f"Equilibrium temperature exceeds melting point at all gaps in range.")

    # Print at key gaps
    for g_check in [0.005, 0.010, 0.020, 0.050, 0.100]:
        idx = min(range(len(gaps)), key=lambda i: abs(gaps[i] - g_check))
        print(f"\n  At gap = {g_check*1000:.0f} mm:")
        print(f"    Heating rate: {heating[idx]:.1f} W/m")
        print(f"    Equilibrium temp: {T_eq[idx]:.0f} K ({T_eq[idx]-273:.0f} C)")

    # Transient analysis: how long can the plate survive at a close gap?
    print(f"\n--- Transient Survival Time ---")
    print(f"  (Time to reach service temperature from cable steady-state {T_CABLE_STEADY} K)")

    for g_check in [0.005, 0.010, 0.020]:
        idx = min(range(len(gaps)), key=lambda i: abs(gaps[i] - g_check))
        Q = heating[idx]
        if Q <= 0:
            continue
        # Thermal mass per meter: density * thickness * width * n_arrays
        thermal_mass = mat["density"] * SAFETY_PLATE_THICKNESS * SAFETY_PLATE_WIDTH * 2
        cp = mat["specific_heat"]
        delta_T = mat["max_service_temp"] - T_CABLE_STEADY

        if delta_T <= 0:
            continue

        # Simple estimate: t = m*cp*dT / Q (neglecting radiation, conservative)
        t_service = thermal_mass * cp * delta_T / Q

        delta_T_melt = mat["melting_point"] - T_CABLE_STEADY
        t_melt = thermal_mass * cp * delta_T_melt / Q

        print(f"\n  At gap = {g_check*1000:.0f} mm:")
        print(f"    Heating rate: {Q:.1f} W/m")
        print(f"    Thermal mass: {thermal_mass:.1f} kg/m")
        print(f"    Time to service temp: {t_service:.1f} s")
        print(f"    Time to melting point: {t_melt:.1f} s")

    return {
        "gaps_mm": [g * 1000 for g in gaps],
        "heating": heating,
        "T_eq": T_eq,
        "service_gap_mm": service_gap * 1000 if service_gap else None,
        "melt_gap_mm": melt_gap * 1000 if melt_gap else None,
        "material": mat,
    }


# =============================================================================
# ANALYSIS 4: MATERIAL TRADE STUDY
# =============================================================================

def run_material_study(wavelength):
    """
    Compare all candidate materials for the safety reaction plate.

    For each material, compute:
    - Repulsive force at key gaps
    - Parasitic drag at nominal gap
    - Equilibrium temperature at key gaps
    - Thermal survival time
    - Mass penalty
    """
    print("\n" + "="*72)
    print("ANALYSIS 4: SAFETY REACTION PLATE MATERIAL TRADE STUDY")
    print("="*72)

    v_cable = 8620  # m/s

    safety_coils_per_m = 4.0 / wavelength
    safety_turns = SAFETY_TURNS_PER_COIL
    B_safety_surface = halbach_surface_field(
        safety_turns, SAFETY_I_PERSISTENT,
        SAFETY_COIL_WIDTH, wavelength,
        coils_per_meter=safety_coils_per_m
    )

    results = {}

    for mat_name, mat in REACTION_PLATE_MATERIALS.items():
        sigma = 1.0 / mat["resistivity"]

        # Force at 20 mm gap (forbidden zone boundary)
        F_20mm = safety_repulsive_force(
            B_safety_surface, 0.020, wavelength,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=2
        )

        # Force at 5 mm gap (very close)
        F_5mm = safety_repulsive_force(
            B_safety_surface, 0.005, wavelength,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=2
        )

        # Drag at nominal gap (100 mm)
        drag_100mm = safety_drag_force(
            B_safety_surface, 0.100, wavelength,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=2
        )

        # Drag power at nominal gap
        drag_power_100mm = drag_100mm * v_cable

        # Heating at 20 mm gap
        Q_20mm = eddy_current_heating_halbach(
            B_safety_surface, 0.020, wavelength, v_cable,
            mat, SAFETY_PLATE_THICKNESS, SAFETY_PLATE_WIDTH,
            n_arrays=2
        )

        # Equilibrium temp at 20 mm gap
        T_eq_20mm = equilibrium_temperature(
            Q_20mm, mat["emissivity"],
            SAFETY_PLATE_WIDTH, n_arrays=2
        )

        # Thermal survival time at 10 mm gap
        Q_10mm = eddy_current_heating_halbach(
            B_safety_surface, 0.010, wavelength, v_cable,
            mat, SAFETY_PLATE_THICKNESS, SAFETY_PLATE_WIDTH,
            n_arrays=2
        )
        thermal_mass = mat["density"] * SAFETY_PLATE_THICKNESS * SAFETY_PLATE_WIDTH * 2
        cp = mat["specific_heat"]
        delta_T_melt = mat["melting_point"] - T_CABLE_STEADY
        t_melt_10mm = thermal_mass * cp * delta_T_melt / Q_10mm if Q_10mm > 0 else float('inf')

        # Mass per meter for safety plates (4 plates total: 2 top, 2 bottom sides)
        mass_per_m = mat["density"] * SAFETY_PLATE_THICKNESS * SAFETY_PLATE_WIDTH * 4

        results[mat_name] = {
            "name": mat["name"],
            "F_20mm_kN": F_20mm / 1000,
            "F_5mm_kN": F_5mm / 1000,
            "drag_100mm_N": drag_100mm,
            "drag_power_W": drag_power_100mm,
            "Q_20mm_W": Q_20mm,
            "T_eq_20mm_K": T_eq_20mm,
            "t_melt_10mm_s": t_melt_10mm,
            "mass_per_m_kg": mass_per_m,
            "melting_point_K": mat["melting_point"],
            "max_service_K": mat["max_service_temp"],
        }

    # Print comparison table
    print(f"\nHalbach wavelength: {wavelength*1000:.0f} mm")
    print(f"Safety Halbach surface field: {B_safety_surface:.3f} T")
    print(f"Cable velocity: {v_cable} m/s")
    print(f"Plate thickness: {SAFETY_PLATE_THICKNESS*1000:.0f} mm")
    print(f"Plate width: {SAFETY_PLATE_WIDTH*1000:.0f} mm")

    print(f"\n{'Material':<30} {'F@20mm':>8} {'F@5mm':>8} {'Drag@100':>10} {'DragPwr':>10} {'T_eq@20':>8} {'t_melt@10':>10} {'Mass':>8}")
    print(f"{'':30} {'kN/m':>8} {'kN/m':>8} {'N/m':>10} {'W/m':>10} {'K':>8} {'sec':>10} {'kg/m':>8}")
    print("-" * 122)

    for mat_name, r in results.items():
        t_melt_str = f"{r['t_melt_10mm_s']:.1f}" if r['t_melt_10mm_s'] < 1e6 else "inf"
        print(f"{r['name']:<30} {r['F_20mm_kN']:>8.1f} {r['F_5mm_kN']:>8.1f} "
              f"{r['drag_100mm_N']:>10.2e} {r['drag_power_W']:>10.2e} "
              f"{r['T_eq_20mm_K']:>8.0f} {t_melt_str:>10} {r['mass_per_m_kg']:>8.1f}")

    # Print recommendations
    print(f"\n--- Recommendations ---")
    print(f"\nKey trade-offs:")
    print(f"  - Copper OFHC: Strongest repulsive force (highest conductivity), but heaviest.")
    print(f"    Best if mass is not the primary constraint.")
    print(f"  - Aluminum 6061: Good force, light, but lowest melting point.")
    print(f"    Risk of melting during prolonged safety engagement.")
    print(f"  - gamma-TiAl: Already on the cable for LIM plates. High service temp.")
    print(f"    Moderate force. A practical choice that avoids adding a new material.")
    print(f"  - Molybdenum: Highest melting point, good conductivity, extremely heavy.")
    print(f"    Best for absolute worst-case thermal survival.")
    print(f"  - Graphene/CNT composite: Lightest, highest temp tolerance in vacuum.")
    print(f"    Speculative conductivity. Best if projections hold.")

    return results


# =============================================================================
# ANALYSIS 5: PARASITIC DRAG AT NOMINAL GAP
# =============================================================================

def run_drag_analysis(material_name, wavelength):
    """
    Compute the parasitic drag from the safety Halbach arrays during
    normal operations (at the nominal 100 mm gap).

    This drag must be small compared to other power losses in the system.
    The Halbach wavelength is the primary tuning parameter: shorter
    wavelength means the field decays faster with distance, giving less
    drag at the nominal gap.
    """
    print("\n" + "="*72)
    print("ANALYSIS 5: PARASITIC DRAG FROM SAFETY ARRAYS AT NOMINAL GAP")
    print("="*72)

    mat = get_material(material_name)
    v_cable = 8620  # m/s
    sigma = 1.0 / mat["resistivity"]

    safety_coils_per_m = 4.0 / wavelength
    safety_turns = SAFETY_TURNS_PER_COIL
    B_safety_surface = halbach_surface_field(
        safety_turns, SAFETY_I_PERSISTENT,
        SAFETY_COIL_WIDTH, wavelength,
        coils_per_meter=safety_coils_per_m
    )

    # Sweep wavelength to show sensitivity
    wavelengths = [0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50]

    print(f"\nMaterial: {mat['name']}")
    print(f"Cable velocity: {v_cable} m/s")
    print(f"Nominal gap: {EML_NOMINAL_GAP*1000:.0f} mm")

    print(f"\n{'Wavelength':>12} {'B_surface':>10} {'B@100mm':>12} {'Drag':>12} {'Drag Power':>12} {'B@20mm':>12} {'F_safety@20':>12}")
    print(f"{'(mm)':>12} {'(T)':>10} {'(T)':>12} {'(N/m)':>12} {'(W/m)':>12} {'(T)':>12} {'(kN/m)':>12}")
    print("-" * 92)

    for wl in wavelengths:
        B_surf = halbach_surface_field(
            EML_TURNS_PER_COIL, SAFETY_I_PERSISTENT,
            SAFETY_COIL_WIDTH, wl
        )
        B_at_100mm = halbach_field_at_gap(B_surf, 0.100, wl)
        B_at_20mm = halbach_field_at_gap(B_surf, 0.020, wl)

        drag = safety_drag_force(
            B_surf, 0.100, wl,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=4  # All 4 safety arrays contribute drag
        )
        drag_power = drag * v_cable

        F_at_20 = safety_repulsive_force(
            B_surf, 0.020, wl,
            sigma, SAFETY_PLATE_THICKNESS,
            SAFETY_PLATE_WIDTH, v_cable,
            n_arrays=2
        )

        print(f"{wl*1000:>12.0f} {B_surf:>10.3f} {B_at_100mm:>12.2e} "
              f"{drag:>12.2e} {drag_power:>12.2e} {B_at_20mm:>12.4f} {F_at_20/1000:>12.1f}")

    # Total drag power for the full ring at the selected wavelength
    B_at_nom = halbach_field_at_gap(B_safety_surface, EML_NOMINAL_GAP, wavelength)
    total_drag = safety_drag_force(
        B_safety_surface, EML_NOMINAL_GAP, wavelength,
        sigma, SAFETY_PLATE_THICKNESS,
        SAFETY_PLATE_WIDTH, v_cable,
        n_arrays=4
    )
    total_drag_power = total_drag * v_cable
    ring_drag_power = total_drag_power * L_RING

    print(f"\nAt selected wavelength = {wavelength*1000:.0f} mm:")
    print(f"  B at nominal gap: {B_at_nom:.2e} T")
    print(f"  Drag per meter: {total_drag:.2e} N/m")
    print(f"  Drag power per meter: {total_drag_power:.2e} W/m")
    print(f"  Total ring drag power: {ring_drag_power:.2e} W ({ring_drag_power/1e6:.3f} MW)")

    # Compare to system power
    print(f"\n  For context:")
    print(f"    LIM deployment power: 666,000 MW")
    print(f"    Post-deployment operations: 42,000 MW")
    print(f"    Safety drag: {ring_drag_power/1e6:.3f} MW")
    if ring_drag_power > 0:
        print(f"    Fraction of operations budget: {ring_drag_power/(42e9)*100:.6f}%")

    return {
        "wavelengths": wavelengths,
        "selected_wavelength": wavelength,
        "total_drag_power_MW": ring_drag_power / 1e6,
    }


# =============================================================================
# PLOTTING
# =============================================================================

def plot_force_comparison(data, output_dir):
    """Plot Halbach vs solenoid force-vs-gap comparison."""
    if not HAS_MATPLOTLIB:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Operating current
    ax1.plot(data["gaps_mm"], data["f_halbach_op"], 'b-', linewidth=2, label='Halbach array')
    ax1.plot(data["gaps_mm"], data["f_solenoid_op"], 'r--', linewidth=2, label='Solenoid')
    ax1.plot(data["gaps_mm"], data["f_ch5_op"], 'g:', linewidth=2, label='Chapter 5 model')
    ax1.axhline(y=data["nominal_load"], color='k', linestyle='-.', alpha=0.5, label=f'Required load ({data["nominal_load"]:.0f} kN/m)')
    ax1.axvline(x=EML_NOMINAL_GAP*1000, color='gray', linestyle=':', alpha=0.3)
    ax1.set_xlabel('Gap (mm)')
    ax1.set_ylabel('Attractive Force (kN/m)')
    ax1.set_title(f'EML Force at Operating Current ({EML_I_OPERATING} A)')
    ax1.legend(fontsize=9)
    ax1.set_xlim(0, 300)
    ax1.grid(True, alpha=0.3)

    # Maximum current
    ax2.plot(data["gaps_mm"], data["f_halbach_max"], 'b-', linewidth=2, label='Halbach array')
    ax2.plot(data["gaps_mm"], data["f_solenoid_max"], 'r--', linewidth=2, label='Solenoid')
    ax2.plot(data["gaps_mm"], data["f_ch5_max"], 'g:', linewidth=2, label='Chapter 5 model')
    ax2.axhline(y=data["nominal_load"], color='k', linestyle='-.', alpha=0.5, label=f'Required load ({data["nominal_load"]:.0f} kN/m)')
    ax2.axvline(x=EML_NOMINAL_GAP*1000, color='gray', linestyle=':', alpha=0.3)
    ax2.set_xlabel('Gap (mm)')
    ax2.set_ylabel('Attractive Force (kN/m)')
    ax2.set_title(f'EML Force at Maximum Current ({EML_HTS_IC} A)')
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, 300)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(output_dir, "fig_levitation_halbach_vs_solenoid.png")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {path}")


def plot_safety_crossover(data, output_dir):
    """Plot safety backstop force vs failure mode loads."""
    if not HAS_MATPLOTLIB:
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Mode 1: Safety vs gravity (EML total failure)
    ax1.semilogy(data["gaps_mm"], [max(f, 1e-3) for f in data["f_safety"]],
                  'b-', linewidth=2.5, label='Safety repulsive force')
    ax1.semilogy(data["gaps_mm"], [max(f, 1e-3) for f in data["f_gravity"]],
                  'r--', linewidth=2, label=f'Gravity load ({data["f_gravity"][0]:.0f} kN/m)')
    ax1.semilogy(data["gaps_mm"], [max(f, 1e-3) for f in data["f_eml_excess"]],
                  'm:', linewidth=2, label='Excess EML (overcurrent at Ic)')
    if data["crossover_mode1_mm"]:
        ax1.axvline(x=data["crossover_mode1_mm"], color='green', linestyle='--', alpha=0.7,
                     label=f'Mode 1 crossover: {data["crossover_mode1_mm"]:.1f} mm')
    ax1.set_xlabel('Gap to Safety Array (mm)')
    ax1.set_ylabel('Force (kN/m)')
    ax1.set_title(f'Safety Force vs Failure Loads ({data["material"]["name"]})')
    ax1.legend(fontsize=9)
    ax1.set_xlim(0, 100)
    ax1.grid(True, alpha=0.3, which='both')

    # Net forces for both modes
    ax2.plot(data["gaps_mm"], data["f_net_mode1"], 'b-', linewidth=2,
              label='Mode 1: Safety - Gravity')
    ax2.plot(data["gaps_mm"], data["f_net_mode2"], 'm--', linewidth=2,
              label='Mode 2: Safety - Excess EML')
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.5)
    ax2.fill_between(data["gaps_mm"], data["f_net_mode1"], 0,
                      where=[f > 0 for f in data["f_net_mode1"]],
                      alpha=0.15, color='green', label='Mode 1: Safety wins')
    ax2.fill_between(data["gaps_mm"], data["f_net_mode1"], 0,
                      where=[f < 0 for f in data["f_net_mode1"]],
                      alpha=0.15, color='red')
    ax2.set_xlabel('Gap to Safety Array (mm)')
    ax2.set_ylabel('Net Force (kN/m)')
    ax2.set_title('Net Force by Failure Mode (positive = safe)')
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, 100)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(output_dir, "fig_levitation_safety_crossover.png")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_thermal(data, output_dir):
    """Plot thermal analysis results."""
    if not HAS_MATPLOTLIB:
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    mat = data["material"]

    # Equilibrium temperature vs gap
    # Cap at 5000 K for display
    T_capped = [min(T, 5000) for T in data["T_eq"]]
    ax1.plot(data["gaps_mm"], T_capped, 'r-', linewidth=2, label='Equilibrium temperature')
    ax1.axhline(y=mat["max_service_temp"], color='orange', linestyle='--',
                 label=f'Service limit ({mat["max_service_temp"]} K)')
    ax1.axhline(y=mat["melting_point"], color='red', linestyle='--',
                 label=f'Melting point ({mat["melting_point"]} K)')
    ax1.axhline(y=T_CABLE_STEADY, color='blue', linestyle=':', alpha=0.5,
                 label=f'Cable steady state ({T_CABLE_STEADY} K)')
    if data["service_gap_mm"]:
        ax1.axvline(x=data["service_gap_mm"], color='orange', linestyle=':', alpha=0.5)
    if data["melt_gap_mm"]:
        ax1.axvline(x=data["melt_gap_mm"], color='red', linestyle=':', alpha=0.5)
    ax1.set_xlabel('Gap to Safety Array (mm)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title(f'Equilibrium Temperature vs Gap ({mat["name"]})')
    ax1.legend(fontsize=9)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, min(max(T_capped[:100]) * 1.2, 5000))
    ax1.grid(True, alpha=0.3)

    # Heating rate vs gap (log scale)
    ax2.semilogy(data["gaps_mm"], [max(h, 1e-10) for h in data["heating"]],
                  'r-', linewidth=2, label='Eddy current heating')
    ax2.set_xlabel('Gap to Safety Array (mm)')
    ax2.set_ylabel('Heating Rate (W/m)')
    ax2.set_title('Eddy Current Heating Rate vs Gap')
    ax2.legend(fontsize=9)
    ax2.set_xlim(0, 100)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    path = os.path.join(output_dir, "fig_levitation_thermal.png")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_material_comparison(results, output_dir):
    """Plot material trade study comparison."""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    names = [r["name"] for r in results.values()]
    short_names = [n.split('(')[0].strip()[:15] for n in names]

    # Force at 20 mm
    ax = axes[0, 0]
    forces = [r["F_20mm_kN"] for r in results.values()]
    bars = ax.barh(short_names, forces, color='steelblue')
    ax.set_xlabel('Repulsive Force at 20 mm (kN/m)')
    ax.set_title('Safety Force at Forbidden Zone Boundary')
    ax.grid(True, alpha=0.3, axis='x')

    # Thermal survival at 10 mm
    ax = axes[0, 1]
    times = [min(r["t_melt_10mm_s"], 1000) for r in results.values()]
    colors = ['green' if t > 60 else 'orange' if t > 10 else 'red' for t in times]
    bars = ax.barh(short_names, times, color=colors)
    ax.set_xlabel('Time to Melt at 10 mm Gap (s)')
    ax.set_title('Thermal Survival Time')
    ax.grid(True, alpha=0.3, axis='x')

    # Mass per meter
    ax = axes[1, 0]
    masses = [r["mass_per_m_kg"] for r in results.values()]
    bars = ax.barh(short_names, masses, color='coral')
    ax.set_xlabel('Mass (kg/m)')
    ax.set_title('Safety Plate Mass Penalty (4 plates)')
    ax.grid(True, alpha=0.3, axis='x')

    # Melting point
    ax = axes[1, 1]
    melt_pts = [r["melting_point_K"] for r in results.values()]
    bars = ax.barh(short_names, melt_pts, color='orange')
    ax.axvline(x=1000, color='r', linestyle=':', alpha=0.5, label='1000 K reference')
    ax.set_xlabel('Melting Point (K)')
    ax.set_title('Melting Point')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')

    plt.suptitle('Safety Reaction Plate Material Trade Study', fontsize=14, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(output_dir, "fig_levitation_material_trade.png")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    material = args["material"]
    wavelength = args["wavelength"]

    # Output directory
    output_dir = os.path.dirname(os.path.abspath(__file__))
    fig_dir = os.path.join(os.path.dirname(output_dir), "localbrain", "figures")
    if os.path.isdir(fig_dir):
        output_dir = fig_dir
    else:
        output_dir = output_dir + "\\output\\levitation\\"

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 72)
    print("ORBITAL RING LEVITATION SYSTEM SIMULATION")
    print("=" * 72)
    print(f"Safety reaction plate material: {material}")
    print(f"Safety Halbach wavelength: {wavelength*1000:.0f} mm")
    print(f"Output directory: {output_dir}")

    show = args["show_all"]
    graphs = args["graphs"]

    # Run analyses
    if show or "force_comparison" in graphs:
        data1 = run_force_comparison()
        plot_force_comparison(data1, output_dir)

    if show or "safety_crossover" in graphs:
        data2 = run_safety_crossover(material, wavelength)
        plot_safety_crossover(data2, output_dir)

    if show or "thermal" in graphs:
        data3 = run_thermal_analysis(material, wavelength)
        plot_thermal(data3, output_dir)

    if show or "materials" in graphs:
        data4 = run_material_study(wavelength)
        plot_material_comparison(data4, output_dir)

    if show or "drag" in graphs:
        data5 = run_drag_analysis(material, wavelength)

    print("\n" + "=" * 72)
    print("SIMULATION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
