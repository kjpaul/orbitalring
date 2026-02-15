#!/usr/bin/env python3
"""
Mass Driver Configuration Module

Orbital ring mass driver launch simulation. Four sequential LIM stages
accelerate a sled from rest to target velocity (e.g. 30 km/s).

Key design features:
  - Stator coils are FIXED ON THE RING at 77 K (cryo-cooled)
  - Reaction plate is ON THE SLED (γ-TiAl, starts at ambient)
  - Four LIM stages with increasing pole pitch
  - Power unlimited (sled draws from entire ring grid)
  - Cryo system on ring in coast mode — barely used, available for cooling
  - Each ring section has a full orbit to recover after sled passes
  - LN2 storage tanks provide thermal buffer

Architecture matches the orbital ring LIM simulator (lim_config.py).
Reference: "Orbital Ring Engineering" by Paul G de Jong
"""
import math

# ── Physical constants ───────────────────────────────────────────────
MU0 = 4 * math.pi * 1e-7
STEFAN_BOLTZMANN = 5.670374e-8
G_ACCEL = 9.80665
T_SPACE = 2.7
T_AMBIENT = 293.15              # Initial plate temperature (K)
T_STATOR = 77.0                 # Stator temperature (K) — HTS cryo LN2

# ── Plate radiation coolant ───────────────────────────────────────────
# Thermal shield between hot plate and cryo coils absorbs plate radiation.
# "water" = water-cooled shield at ~350 K (large liquid range, easy to reject)
# "ln2"   = plate radiates directly to 77 K cryo surfaces
PLATE_COOLANT = "water"

# Water properties (plate radiation coolant)
T_WATER_SUPPLY = 300.0         # K — supply temperature (~27°C)
T_WATER_BOIL = 373.0           # K — boiling point at 1 atm
C_P_WATER = 4186               # J/(kg·K) — specific heat
L_V_WATER = 2_260_000          # J/kg — latent heat of vaporization
RHO_WATER = 1000.0             # kg/m³
T_SHIELD_WATER = 350.0         # K — thermal shield operating temperature

# LN2 properties (HTS cryo system; also plate coolant if PLATE_COOLANT="ln2")
T_LN2_BOIL = 77.4              # K — boiling point at 1 atm
T_LN2_SUPPLY = 70.0            # K — subcooled supply temperature
C_P_LN2 = 2040                 # J/(kg·K) — specific heat (liquid)
L_V_LN2 = 199_000              # J/kg — latent heat of vaporization
RHO_LN2 = 807.0                # kg/m³ — liquid density at 77 K

# ── HTS tape and Gömöry hysteresis parameters ────────────────────
ALPHA_TAPE_FIELD = 20.0                      # degrees — Gömöry perpendicular field angle
SIN_ALPHA = math.sin(math.radians(ALPHA_TAPE_FIELD))  # 0.34202
I_C_PER_MM_LAYER = 66.7                      # A/mm-width/layer (critical current density)

# ── Orbital ring parameters ─────────────────────────────────────────
R_EARTH_EQ = 6_378_137.0
ALT = 250_000.0
R_ORBIT = R_EARTH_EQ + ALT
L_RING = 2 * math.pi * R_ORBIT  # ≈ 41,645,813 m
V_ORBIT = 7_754.866
G_LOCAL = 9.073
V_GROUND_SYNC = 483.331

RING_LIMIT_10G = 1_674_000
RING_LIMIT_3G  = 5_581_000
RING_LIMIT_1G  = 16_744_000
RING_SAFETY_FACTOR = 0.50

# ── LIM stage definitions ───────────────────────────────────────────

class LIMStageConfig:
    """Configuration for one LIM stage.

    Per-stage settable parameters:
      tau_p        — pole pitch (m)
      pitches      — number of pole pitches per stage section
      n_turns      — turns per phase coil
      hts_width_mm — HTS tape width (3, 4, 6, or 12 mm)
      hts_layers   — number of HTS tape layers (1 = no de-rating)
      w_coil       — coil width (m), max 2.0
      gap          — mechanical air gap (m)
      f_handoff_hz — supply frequency at which to hand off to next stage

    Derived quantities:
      de_rating    — Ic reduction factor for multi-layer stacks:
                     1.0 for 1 layer, (1 - sin(α)) per layer for >1
      I_c_per_layer — effective critical current per layer (A)
      I_c          — total critical current (A)
      I_peak       — max operating current, 90% of I_c (A)
      I_target     — nominal operating current, 80% of I_c (A)
      L_active     — active stator length per section (m)
    """
    def __init__(self, name, tau_p, pitches, n_turns, hts_width_mm,
                 hts_layers, w_coil, gap, f_handoff_hz):
        self.name = name
        self.tau_p = tau_p
        self.pitches = pitches
        self.n_turns = n_turns
        self.hts_width_mm = hts_width_mm
        self.hts_layers = hts_layers
        self.w_coil = w_coil
        self.gap = gap
        self.f_handoff_hz = f_handoff_hz

        # Derived geometry
        self.L_active = pitches * tau_p + (2.0 / 3.0) * tau_p

        # HTS critical current with de-rating for multi-layer stacks
        # Single layer: full Ic.  Multiple layers: perpendicular field
        # from adjacent layers reduces Ic by sin(α) per layer.
        if hts_layers <= 1:
            self.de_rating = 1.0
        else:
            self.de_rating = 1.0 - SIN_ALPHA       # 0.658

        I_c_per_layer_nominal = I_C_PER_MM_LAYER * hts_width_mm   # A
        self.I_c_per_layer = I_c_per_layer_nominal * self.de_rating
        self.I_c = self.I_c_per_layer * hts_layers
        self.I_peak = 0.90 * self.I_c
        self.I_target = 0.80 * self.I_c

# Stage 1: 0 → ~0.6 km/s, τ_p=32m, f=0–10 Hz
# Large τ_p keeps HTS hysteresis under 500 kW; per-stage tiling
# gives n_active ≈ 35 on the 10 km sled for strong thrust.
STAGE_S1 = LIMStageConfig("S1", tau_p=32.0, pitches=8, n_turns=360,
    hts_width_mm=3.0, hts_layers=16, w_coil=2.0, gap=0.100, f_handoff_hz=10.0)

# Stage 2: ~0.6 → ~3.2 km/s, τ_p=56m, f=5–30 Hz
STAGE_S2 = LIMStageConfig("S2", tau_p=56.0, pitches=8, n_turns=360,
    hts_width_mm=3.0, hts_layers=16, w_coil=2.0, gap=0.100, f_handoff_hz=30.0)

# Stage 3: ~3.2 → ~10.2 km/s, τ_p=82m, f=20–66 Hz
STAGE_S3 = LIMStageConfig("S3", tau_p=82.0, pitches=8, n_turns=360,
    hts_width_mm=3.0, hts_layers=16, w_coil=2.0, gap=0.100, f_handoff_hz=66.0)

# Stage 4: ~10.2 → 30 km/s, τ_p=130m, f=40–124 Hz
STAGE_S4 = LIMStageConfig("S4", tau_p=130.0, pitches=8, n_turns=360,
    hts_width_mm=3.0, hts_layers=16, w_coil=2.0, gap=0.100, f_handoff_hz=None)

LIM_STAGES = [STAGE_S1, STAGE_S2, STAGE_S3, STAGE_S4]

L_REPEATING_UNIT = sum(s.L_active for s in LIM_STAGES)
L_STAGE_GAP = 5.0
L_REPEATING_UNIT_WITH_GAPS = L_REPEATING_UNIT + L_STAGE_GAP * len(LIM_STAGES)

# ── Sled and reaction plate ─────────────────────────────────────────
PLATE_HEIGHT = 2.0          # 2000 mm — optimized for air-core B
PLATE_THICKNESS = 0.100     # 100 mm
SLED_LENGTH = 10_000.0
PLATE_MATERIAL = "molybdenum"
SLED_STRUCTURE_FRACTION = 0.15
N_LIM_SIDES = 2            # LIMs on both sides of the plate

# ── Materials database ───────────────────────────────────────────────
MATERIALS = {
    "aluminum": {
        "name": "Aluminum 6061-T6",
        "rho_e": 3.99e-8, "density": 2700, "Cp": 896, "T_max": 750,
        "alpha_rho": 0.0040, "emissivity": 0.25, "UTS": 310e6,
    },
    "gamma_tial": {
        "name": "γ-TiAl (Ti-48Al-2Cr-2Nb)",
        "rho_e": 1.55e-6, "density": 3900, "Cp": 620, "T_max": 1200,
        "alpha_rho": 0.0015, "emissivity": 0.35, "UTS": 450e6,
    },
    "cuni7030": {
        "name": "CuNi 70/30",
        "rho_e": 3.75e-7, "density": 8900, "Cp": 377, "T_max": 800,
        "alpha_rho": 0.00020, "emissivity": 0.30, "UTS": 540e6,
    },
    "inconel718": {
        "name": "Inconel 718",
        "rho_e": 1.25e-6, "density": 8190, "Cp": 435, "T_max": 1200,
        "alpha_rho": 0.0010, "emissivity": 0.35, "UTS": 1240e6,
    },
    "copper": {
        "name": "Copper C10200 (OFHC)",
        "rho_e": 1.72e-8, "density": 8960, "Cp": 385, "T_max": 600,
        "alpha_rho": 0.0039, "emissivity": 0.15, "UTS": 220e6,
    },
    "stainless316": {
        "name": "Stainless Steel 316L",
        "rho_e": 7.40e-7, "density": 7990, "Cp": 500, "T_max": 1100,
        "alpha_rho": 0.0010, "emissivity": 0.40, "UTS": 485e6,
    },
    "molybdenum": {
        "name": "Molybdenum (Mo)",
        "rho_e": 5.3e-8, "density": 10220, "Cp": 251, "T_max": 2000,
        "alpha_rho": 0.0047, "emissivity": 0.30, "UTS": 585e6,
    },
    "mo_tzm": {
        "name": "Mo-TZM (Mo-0.5Ti-0.1Zr)",
        "rho_e": 5.7e-8, "density": 10160, "Cp": 255, "T_max": 1900,
        "alpha_rho": 0.0047, "emissivity": 0.35, "UTS": 860e6,
    },
    "niobium": {
        "name": "Niobium (Nb)",
        "rho_e": 1.52e-7, "density": 8570, "Cp": 265, "T_max": 2200,
        "alpha_rho": 0.0039, "emissivity": 0.35, "UTS": 275e6,
    },
    "nb_c103": {
        "name": "Nb-C103 (Nb-10Hf-1Ti)",
        "rho_e": 1.36e-7, "density": 8850, "Cp": 272, "T_max": 1920,
        "alpha_rho": 0.0038, "emissivity": 0.35, "UTS": 390e6,
    },
}

# ── Operating limits ─────────────────────────────────────────────────
VOLTS_COIL_MAX = 100_000.0      # Max coil insulation voltage (V) — cryo limit

# ── Launch mission ───────────────────────────────────────────────────
V_LAUNCH = 30_000.0
LAUNCH_CLASS = "3g"
MAX_ACCEL_G = 0.5
M_SPACECRAFT = 5_000_000.0

# ── Simulation ───────────────────────────────────────────────────────
DT = 0.1
DT_QUICK = 1.0
V_INITIAL = 0.0
THRUST_MODEL = 1                     # Model 1 only — eddy current, narrow plate
GRAPH_DIR = "graphs_lim_md_3"         # Output directory for plots and CSV

# ── Derived ──────────────────────────────────────────────────────────
def get_material():
    return MATERIALS[PLATE_MATERIAL]

def calc_derived():
    mat = get_material()
    plate_xs = PLATE_HEIGHT * PLATE_THICKNESS
    plm = plate_xs * mat["density"]
    pm = plm * SLED_LENGTH
    ssm = SLED_STRUCTURE_FRACTION * pm
    sm = pm + ssm
    tm = sm + M_SPACECRAFT
    lim = {"10g": RING_LIMIT_10G, "3g": RING_LIMIT_3G, "1g": RING_LIMIT_1G}
    rl = lim.get(LAUNCH_CLASS, RING_LIMIT_3G) * RING_SAFETY_FACTOR
    na = SLED_LENGTH / L_REPEATING_UNIT_WITH_GAPS
    tc = pm * mat["Cp"] * (mat["T_max"] - T_AMBIENT)
    ke = 0.5 * tm * V_LAUNCH**2
    fa = tm * MAX_ACCEL_G * G_ACCEL
    tl = V_LAUNCH / (MAX_ACCEL_G * G_ACCEL)
    dl = 0.5 * MAX_ACCEL_G * G_ACCEL * tl**2
    return dict(mat=mat, plate_xs=plate_xs, plate_linear_mass=plm,
                plate_mass=pm, sled_structure_mass=ssm, sled_mass=sm,
                total_mass=tm, ring_limit=rl, n_active=na,
                thermal_capacity=tc, KE_launch=ke, F_max_accel=fa,
                t_launch_est=tl, d_launch_est=dl)

def print_config():
    d = calc_derived()
    mat = d["mat"]
    sep = "=" * 76
    print(f"\n{sep}")
    print("MASS DRIVER LAUNCH CONFIGURATION")
    print(sep)
    print(f"  Target velocity        {V_LAUNCH/1000:>8.1f} km/s    Max accel     {MAX_ACCEL_G:>8.2f} g")
    print(f"  Launch class           {LAUNCH_CLASS:>8}        Thrust model  {THRUST_MODEL:>8} (eddy)")
    print(f"  Stator temp            {T_STATOR:>8.0f} K       Plate T_init  {T_AMBIENT:>8.1f} K")
    print(f"  Stator core            {'air-core':>8}       LIM sides     {N_LIM_SIDES:>8}")
    print(f"\n  PLATE: {mat['name']}, {PLATE_HEIGHT*1000:.0f}×{PLATE_THICKNESS*1000:.0f} mm, {SLED_LENGTH/1000:.0f} km, LIMs on {N_LIM_SIDES} sides")
    print(f"    Linear mass {d['plate_linear_mass']:.0f} kg/m, Total {d['plate_mass']/1000:.0f} t, "
          f"ρ_e={mat['rho_e']:.2e} Ω·m, T_max={mat['T_max']} K")
    print(f"\n  MASS: plate {d['plate_mass']/1000:.0f}t + structure {d['sled_structure_mass']/1000:.0f}t "
          f"+ spacecraft {M_SPACECRAFT/1000:.0f}t = {d['total_mass']/1000:.0f}t "
          f"(limit {d['ring_limit']/1000:.0f}t, "
          f"{'OK +'+str(int((d['ring_limit']-d['total_mass'])/1000))+'t' if d['ring_limit']>=d['total_mass'] else 'OVER '+str(int((d['total_mass']-d['ring_limit'])/1000))+'t'})")
    print(f"\n  STAGES:")
    for s in LIM_STAGES:
        hoff = f"→{s.f_handoff_hz:.0f}Hz" if s.f_handoff_hz else "→END"
        na_stage = SLED_LENGTH / (s.L_active + L_STAGE_GAP)
        print(f"    {s.name}: τ_p={s.tau_p:>5.0f}m  L={s.L_active:>6.1f}m  N={s.n_turns:>4}  "
              f"{s.hts_layers}×{s.hts_width_mm:.0f}mm  "
              f"I_c={s.I_c:>6.0f}A  gap={s.gap*1000:.0f}mm  "
              f"n_a={na_stage:>5.1f}  {hoff}")
    print(f"    Repeating unit {L_REPEATING_UNIT_WITH_GAPS:.1f} m, "
          f"{d['n_active']:.1f} active on {SLED_LENGTH/1000:.0f} km sled")
    print(f"\n  EST: {d['t_launch_est']:.0f}s ({d['t_launch_est']/60:.1f}min), "
          f"{d['d_launch_est']/1000:.0f}km, {d['F_max_accel']/1e6:.1f}MN, "
          f"KE={d['KE_launch']/3.6e12:.1f}GWh")
    print(sep)

if __name__ == "__main__":
    print_config()
