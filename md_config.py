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
T_STATOR = 77.0                 # Stator temperature (K) — cryo LN2

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
        self.L_active = pitches * tau_p + (2.0/3.0) * tau_p
        self.Ic_per_mm_per_layer = 66.7
        self.I_c = self.Ic_per_mm_per_layer * hts_width_mm * hts_layers
        self.I_peak = 0.90 * self.I_c
        self.I_target = 0.80 * self.I_c

# Stage 1: 0 → ~0.5 km/s, τ_p=3m, f=0–83 Hz, excellent EM coupling
STAGE_S1 = LIMStageConfig("S1", tau_p=3.0, pitches=5, n_turns=360,
    hts_width_mm=12, hts_layers=3, w_coil=2.0, gap=0.100, f_handoff_hz=100.0)

# Stage 2: ~0.5 → ~3 km/s, τ_p=5m, f=50–300 Hz, good coupling
STAGE_S2 = LIMStageConfig("S2", tau_p=5.0, pitches=5, n_turns=360,
    hts_width_mm=12, hts_layers=3, w_coil=2.0, gap=0.100, f_handoff_hz=300.0)

# Stage 3: ~3 → ~8 km/s, τ_p=10m, f=150–400 Hz, moderate coupling
STAGE_S3 = LIMStageConfig("S3", tau_p=10.0, pitches=3, n_turns=360,
    hts_width_mm=12, hts_layers=3, w_coil=2.0, gap=0.100, f_handoff_hz=400.0)

# Stage 4: ~8 → 15 km/s, τ_p=20m, f=200–375 Hz, fair coupling
STAGE_S4 = LIMStageConfig("S4", tau_p=20.0, pitches=3, n_turns=360,
    hts_width_mm=12, hts_layers=3, w_coil=2.0, gap=0.100, f_handoff_hz=None)

LIM_STAGES = [STAGE_S1, STAGE_S2, STAGE_S3, STAGE_S4]

L_REPEATING_UNIT = sum(s.L_active for s in LIM_STAGES)
L_STAGE_GAP = 5.0
L_REPEATING_UNIT_WITH_GAPS = L_REPEATING_UNIT + L_STAGE_GAP * len(LIM_STAGES)

# ── Sled and reaction plate ─────────────────────────────────────────
PLATE_HEIGHT = 2.0          # 2000 mm — optimized for air-core B
PLATE_THICKNESS = 0.100     # 100 mm
SLED_LENGTH = 5000.0
PLATE_MATERIAL = "gamma_tial"
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
}

# ── Operating limits ─────────────────────────────────────────────────
VOLTS_MAX = 100_000.0       # Max supply voltage with PFC (V)

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
GRAPH_DIR = "graphs_lim_md"         # Output directory for plots and CSV

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
        print(f"    {s.name}: τ_p={s.tau_p:>5.0f}m  L={s.L_active:>6.1f}m  N={s.n_turns:>4}  "
              f"I_c={s.I_c:>6.0f}A  gap={s.gap*1000:.0f}mm  {hoff}")
    print(f"    Repeating unit {L_REPEATING_UNIT_WITH_GAPS:.1f} m, "
          f"{d['n_active']:.1f} active on {SLED_LENGTH/1000:.0f} km sled")
    print(f"\n  EST: {d['t_launch_est']:.0f}s ({d['t_launch_est']/60:.1f}min), "
          f"{d['d_launch_est']/1000:.0f}km, {d['F_max_accel']/1e6:.1f}MN, "
          f"KE={d['KE_launch']/3.6e12:.1f}GWh")
    print(sep)

if __name__ == "__main__":
    print_config()
