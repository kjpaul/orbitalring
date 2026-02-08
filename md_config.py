#!/usr/bin/env python3
"""
Mass Driver Configuration Module - Parameters and material properties

This module contains all configurable parameters for the orbital ring
mass driver launch simulation. The mass driver accelerates a sled carrying
a spacecraft along the orbital ring using three sequential LIM stages with
increasing pole pitches.

Architecture follows the orbital ring LIM simulator (lim_config.py).

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import math

# =============================================================================
# SECTION 1: PHYSICAL CONSTANTS
# =============================================================================

MU0 = 4 * math.pi * 1e-7          # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.670374e-8     # Stefan-Boltzmann constant (W/(m²·K⁴))
G_ACCEL = 9.80665                  # Standard gravity (m/s²)
T_SPACE = 2.7                      # Background space temperature (K)
T_AMBIENT = 293.15                 # Ambient temperature for material props (K)

# =============================================================================
# SECTION 2: ORBITAL RING PARAMETERS
# =============================================================================

R_EARTH_EQ = 6_378_137.0           # Earth equatorial radius (m)
ALT = 250_000.0                    # Orbital altitude (m)
R_ORBIT = R_EARTH_EQ + ALT         # Orbital radius (m)
L_RING = 2 * math.pi * R_ORBIT     # Ring circumference (m) ≈ 41,645,813 m
V_ORBIT = 7_754.866                # Orbital velocity at 250 km (m/s)
G_LOCAL = 9.073                    # Local gravity at 250 km (m/s²)
V_GROUND_SYNC = 483.331            # Ground-synchronous velocity at 250 km (m/s)

# Orbital ring structural limits (at 1000 m mass distribution)
# These are breaking limits — operational limit is ~50%
RING_LIMIT_10G = 1_674_000         # Max launch mass at 10g (kg) — unmanned cargo
RING_LIMIT_3G  = 5_581_000         # Max launch mass at 3g (kg) — human-rated
RING_LIMIT_1G  = 16_744_000        # Max launch mass at 1g (kg) — heavy cargo
RING_SAFETY_FACTOR = 0.50          # Operating at 50% of breaking limit

# =============================================================================
# SECTION 3: LIM STAGE DEFINITIONS
#
# Three LIM stages with increasing pole pitch, placed back-to-back:
#   Stage 1 (LV): Low velocity,  τ_p = 10 m,  handles 0 → ~1 km/s
#   Stage 2 (MV): Mid velocity,  τ_p = 60 m,  handles ~1 → ~5 km/s
#   Stage 3 (HV): High velocity, τ_p = 200 m, handles ~5 → 15 km/s
#
# L_LIM(τ_p, pitches) = pitches × τ_p + (2/3) × τ_p
# L_all_LIMs = L_LIM(10,3) + L_LIM(60,3) + L_LIM(200,3)
#            = 36.67 + 206.67 + 733.33 = 976.67 m ≈ 990 m (w/ gaps)
# =============================================================================

class LIMStageConfig:
    """Configuration for one LIM stage."""
    def __init__(self, name, tau_p, pitches, n_turns, hts_width_mm,
                 hts_layers, w_coil, gap, f_handoff_hz):
        self.name = name
        self.tau_p = tau_p              # Pole pitch (m)
        self.pitches = pitches          # Number of pole pitches
        self.n_turns = n_turns          # Turns per phase coil
        self.hts_width_mm = hts_width_mm  # HTS tape width (mm)
        self.hts_layers = hts_layers    # HTS tape layers in parallel
        self.w_coil = w_coil            # Coil width (m)
        self.gap = gap                  # Air gap coil-to-plate (m)
        self.f_handoff_hz = f_handoff_hz  # Handoff frequency (Hz)

        # Derived: active length of this stage
        self.L_active = pitches * tau_p + (2.0/3.0) * tau_p

        # HTS current limits (same model as orbital ring code)
        self.Ic_per_mm_per_layer = 66.7  # A/mm per layer at 77K, self-field
        self.I_c = self.Ic_per_mm_per_layer * hts_width_mm * hts_layers
        self.I_peak = 0.90 * self.I_c    # 90% of critical for safety
        self.I_target = 0.80 * self.I_c   # 80% target operating point

    def __repr__(self):
        return (f"LIMStage({self.name}: τ_p={self.tau_p}m, "
                f"L={self.L_active:.1f}m, I_c={self.I_c:.0f}A)")


# Stage 1: Low velocity (0 → ~1 km/s)
# Many turns for high B at low frequency; V_back stays manageable at low v_sled
STAGE_LV = LIMStageConfig(
    name="LV",
    tau_p=10.0,
    pitches=3,
    n_turns=360,
    hts_width_mm=12,
    hts_layers=3,
    w_coil=0.5,
    gap=0.100,          # 100 mm gap
    f_handoff_hz=100.0  # Hand off when f_supply reaches 100 Hz
)

# Stage 2: Mid velocity (~1 → ~5 km/s)
# With PFC, voltage scales with v_slip (~50 m/s), not v_wave (~3000 m/s)
# So we can maintain high turn count for strong B-field
STAGE_MV = LIMStageConfig(
    name="MV",
    tau_p=60.0,
    pitches=3,
    n_turns=360,
    hts_width_mm=12,
    hts_layers=3,
    w_coil=0.5,
    gap=0.100,
    f_handoff_hz=60.0   # Hand off when f_supply reaches 60 Hz
)

# Stage 3: High velocity (~5 → 15 km/s)
# With PFC and v_slip ~30 m/s: V = 360 × 2.17 × 0.5 × 30 × 0.64 ≈ 7.5 kV
# Well within 100 kV limit
STAGE_HV = LIMStageConfig(
    name="HV",
    tau_p=200.0,
    pitches=3,
    n_turns=360,
    hts_width_mm=12,
    hts_layers=3,
    w_coil=0.5,
    gap=0.100,
    f_handoff_hz=None   # Final stage — runs to target velocity
)

# Ordered list of stages (sequential operation)
LIM_STAGES = [STAGE_LV, STAGE_MV, STAGE_HV]

# Total length of one repeating unit (all three LIMs back-to-back + gaps)
L_REPEATING_UNIT = sum(s.L_active for s in LIM_STAGES)
L_REPEATING_UNIT_WITH_GAPS = 990.0  # Including structural transitions (m)

# =============================================================================
# SECTION 4: SLED AND REACTION PLATE PARAMETERS
# =============================================================================

# Reaction plate geometry
PLATE_HEIGHT = 0.5          # Plate height (vertical dimension, m) — was "width" in earlier convos
PLATE_THICKNESS = 0.200     # Plate thickness (EM interaction depth, m)
SLED_LENGTH = 5000.0        # Sled length (m) — 5 km

# Reaction plate material (default: γ-TiAl)
PLATE_MATERIAL = "gamma_tial"

# Sled structural mass (frame, ribs, couplers, etc.)
# Estimated as 15% of plate mass
SLED_STRUCTURE_FRACTION = 0.15

# =============================================================================
# SECTION 5: MATERIAL PROPERTIES DATABASE
#
# Same material database as the orbital ring code for consistency.
# Properties at 293 K unless noted. The simulation uses temperature-dependent
# resistivity where available.
# =============================================================================

MATERIALS = {
    "aluminum": {
        "name": "Aluminum 6061-T6",
        "rho_e": 3.99e-8,          # Electrical resistivity (Ω·m) at 293K
        "density": 2700,            # kg/m³
        "Cp": 896,                  # Specific heat (J/(kg·K))
        "T_max": 750,               # Max service temperature (K)
        "alpha_rho": 0.0040,        # Temp coefficient of resistivity (1/K)
        "emissivity": 0.25,         # Thermal emissivity
        "UTS": 310e6,               # Ultimate tensile strength (Pa)
    },
    "gamma_tial": {
        "name": "γ-TiAl (Ti-48Al-2Cr-2Nb)",
        "rho_e": 1.55e-6,          # Electrical resistivity (Ω·m) at 293K
        "density": 3900,            # kg/m³ (lighter than steel, heavier than Al)
        "Cp": 620,                  # Specific heat (J/(kg·K))
        "T_max": 1200,              # Max service temperature (K) — outstanding high-T
        "alpha_rho": 0.0015,        # Temp coefficient of resistivity (1/K)
        "emissivity": 0.35,         # Thermal emissivity (oxidized surface)
        "UTS": 450e6,               # Ultimate tensile strength (Pa)
    },
    "cuni7030": {
        "name": "CuNi 70/30 (Monel-like)",
        "rho_e": 3.75e-7,          # Electrical resistivity (Ω·m) at 293K
        "density": 8900,            # kg/m³
        "Cp": 377,                  # Specific heat (J/(kg·K))
        "T_max": 800,               # Max service temperature (K)
        "alpha_rho": 0.00020,       # Very low temp coefficient
        "emissivity": 0.30,
        "UTS": 540e6,
    },
    "inconel718": {
        "name": "Inconel 718",
        "rho_e": 1.25e-6,          # Electrical resistivity (Ω·m) at 293K
        "density": 8190,            # kg/m³
        "Cp": 435,                  # Specific heat (J/(kg·K))
        "T_max": 1200,              # Max service temperature (K)
        "alpha_rho": 0.0010,
        "emissivity": 0.35,
        "UTS": 1240e6,
    },
    "copper": {
        "name": "Copper C10200 (OFHC)",
        "rho_e": 1.72e-8,          # Electrical resistivity (Ω·m) at 293K
        "density": 8960,
        "Cp": 385,
        "T_max": 600,
        "alpha_rho": 0.0039,
        "emissivity": 0.15,
        "UTS": 220e6,
    },
    "stainless316": {
        "name": "Stainless Steel 316L",
        "rho_e": 7.40e-7,
        "density": 7990,
        "Cp": 500,
        "T_max": 1100,
        "alpha_rho": 0.0010,
        "emissivity": 0.40,
        "UTS": 485e6,
    },
}

# =============================================================================
# SECTION 6: OPERATING LIMITS
# =============================================================================

VOLTS_MAX = 100_000.0       # Maximum voltage (V) — 100 kV limit
# Power is NOT limited — sled draws from entire ring grid
# Current limits are per-stage from HTS tape (see LIMStageConfig)

# =============================================================================
# SECTION 7: LAUNCH MISSION PARAMETERS
# =============================================================================

V_LAUNCH = 15_000.0         # Target launch velocity (m/s) — Mars transfer
LAUNCH_CLASS = "3g"          # "10g", "3g", or "1g" — sets ring loading limit
MAX_ACCEL_G = 0.5            # Maximum acceleration (g) — human-rated comfort

# Spacecraft mass (what the sled carries, not including the sled itself)
M_SPACECRAFT = 500_000.0    # Spacecraft mass (kg) — fits within 3g ring limit

# =============================================================================
# SECTION 8: SIMULATION PARAMETERS
# =============================================================================

DT = 0.1                    # Time step (s)
DT_QUICK = 1.0              # Quick-run time step (s)
V_INITIAL = 0.0             # Initial sled velocity (m/s) — stationary on ring

THRUST_MODEL = 2             # Default thrust model: 1=eddy, 2=goodness, 3=slip×pressure

# =============================================================================
# SECTION 9: DERIVED PARAMETERS
# =============================================================================

def get_material():
    """Return the material properties dict for the configured plate material."""
    return MATERIALS[PLATE_MATERIAL]

def calc_derived():
    """Calculate all derived parameters and return as a dict."""
    mat = get_material()

    # Plate cross-section and mass
    plate_cross_section = PLATE_HEIGHT * PLATE_THICKNESS  # m²
    plate_linear_mass = plate_cross_section * mat["density"]  # kg/m
    plate_mass = plate_linear_mass * SLED_LENGTH  # kg

    # Sled structure mass (ribs, frame, couplers)
    sled_structure_mass = SLED_STRUCTURE_FRACTION * plate_mass  # kg

    # Total sled mass (plate + structure, no spacecraft)
    sled_mass = plate_mass + sled_structure_mass  # kg

    # Total launch mass
    total_mass = sled_mass + M_SPACECRAFT  # kg

    # Ring loading limit
    if LAUNCH_CLASS == "10g":
        ring_limit = RING_LIMIT_10G * RING_SAFETY_FACTOR
    elif LAUNCH_CLASS == "3g":
        ring_limit = RING_LIMIT_3G * RING_SAFETY_FACTOR
    else:
        ring_limit = RING_LIMIT_1G * RING_SAFETY_FACTOR

    # Number of active repeating units on the sled at any time
    n_active_units = SLED_LENGTH / L_REPEATING_UNIT_WITH_GAPS

    # Thermal capacity of the plate
    thermal_capacity = plate_mass * mat["Cp"] * (mat["T_max"] - T_AMBIENT)  # J

    # Kinetic energy at launch velocity
    KE_launch = 0.5 * total_mass * V_LAUNCH**2  # J

    # Maximum acceleration force
    F_max_accel = total_mass * MAX_ACCEL_G * G_ACCEL  # N

    # Launch time estimate (constant acceleration)
    t_launch_est = V_LAUNCH / (MAX_ACCEL_G * G_ACCEL)  # s

    # Launch distance estimate
    d_launch_est = 0.5 * MAX_ACCEL_G * G_ACCEL * t_launch_est**2  # m

    return {
        "mat": mat,
        "plate_cross_section": plate_cross_section,
        "plate_linear_mass": plate_linear_mass,
        "plate_mass": plate_mass,
        "sled_structure_mass": sled_structure_mass,
        "sled_mass": sled_mass,
        "total_mass": total_mass,
        "ring_limit": ring_limit,
        "n_active_units": n_active_units,
        "thermal_capacity": thermal_capacity,
        "KE_launch": KE_launch,
        "F_max_accel": F_max_accel,
        "t_launch_est": t_launch_est,
        "d_launch_est": d_launch_est,
    }


# =============================================================================
# SECTION 10: DISPLAY
# =============================================================================

def print_config():
    """Print all configuration parameters in a formatted table."""
    d = calc_derived()
    mat = d["mat"]

    sep = "=" * 72
    print(f"\n{sep}")
    print("MASS DRIVER LAUNCH CONFIGURATION")
    print(sep)

    print(f"\n  MISSION")
    print(f"    Target velocity          {V_LAUNCH/1000:>10.1f} km/s")
    print(f"    Launch class                    {LAUNCH_CLASS:>6}")
    print(f"    Max acceleration         {MAX_ACCEL_G:>10.2f} g")

    print(f"\n  REACTION PLATE ({mat['name']})")
    print(f"    Height                   {PLATE_HEIGHT*1000:>10.0f} mm")
    print(f"    Thickness                {PLATE_THICKNESS*1000:>10.0f} mm")
    print(f"    Sled length              {SLED_LENGTH/1000:>10.1f} km")
    print(f"    Cross-section            {d['plate_cross_section']:>10.4f} m²")
    print(f"    Linear mass              {d['plate_linear_mass']:>10.1f} kg/m")
    print(f"    Plate mass               {d['plate_mass']/1000:>10.0f} tonne")
    print(f"    Density                  {mat['density']:>10.0f} kg/m³")
    print(f"    Cp                       {mat['Cp']:>10.0f} J/(kg·K)")
    print(f"    T_max                    {mat['T_max']:>10.0f} K")
    print(f"    ρ_e (293K)               {mat['rho_e']:>12.2e} Ω·m")

    print(f"\n  MASS BUDGET")
    print(f"    Reaction plate           {d['plate_mass']/1000:>10.0f} tonne")
    print(f"    Sled structure           {d['sled_structure_mass']/1000:>10.0f} tonne")
    print(f"    Spacecraft               {M_SPACECRAFT/1000:>10.0f} tonne")
    print(f"    TOTAL                    {d['total_mass']/1000:>10.0f} tonne")
    print(f"    Ring limit ({LAUNCH_CLASS})         {d['ring_limit']/1000:>10.0f} tonne")
    margin = d['ring_limit'] - d['total_mass']
    if margin < 0:
        print(f"    *** OVER LIMIT BY        {-margin/1000:>10.0f} tonne ***")
    else:
        print(f"    Margin                   {margin/1000:>10.0f} tonne")

    print(f"\n  LIM STAGES")
    for s in LIM_STAGES:
        print(f"    {s.name}: τ_p={s.tau_p:>6.0f} m, L={s.L_active:>7.1f} m, "
              f"N={s.n_turns:>4}, HTS {s.hts_width_mm}mm×{s.hts_layers}L, "
              f"I_c={s.I_c:>7.0f} A, gap={s.gap*1000:.0f} mm"
              f"{f', handoff @ {s.f_handoff_hz:.0f} Hz' if s.f_handoff_hz else ''}")
    print(f"    Repeating unit           {L_REPEATING_UNIT_WITH_GAPS:>10.0f} m")
    print(f"    Active units on sled     {d['n_active_units']:>10.1f}")

    print(f"\n  ESTIMATES (constant {MAX_ACCEL_G}g)")
    print(f"    Launch time              {d['t_launch_est']:>10.0f} s  ({d['t_launch_est']/60:>.1f} min)")
    print(f"    Launch distance          {d['d_launch_est']/1000:>10.0f} km  ({d['d_launch_est']/L_RING*100:.2f}% of ring)")
    print(f"    Kinetic energy           {d['KE_launch']:>12.3e} J  ({d['KE_launch']/3.6e12:.1f} GWh)")
    print(f"    Thermal capacity         {d['thermal_capacity']:>12.3e} J  ({d['thermal_capacity']/d['KE_launch']*100:.1f}% of KE)")
    print(f"    Required thrust          {d['F_max_accel']/1e6:>10.2f} MN")
    print(sep)


if __name__ == "__main__":
    print_config()
