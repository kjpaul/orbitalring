import math

# --- 1. MATERIAL DATABASE ---
MATERIALS = {
    "gamma_titanium": {
        "name": "Gamma TiAl (Ti-48Al-2Cr-2Nb)",
        "rho_293": 7.5e-7,  # High resistivity
        "alpha": 0.0012,    # Temp coeff
        "density": 3900,
        "cp": 600,          # Specific heat
        "max_temp": 1100    # Softening point (K)
    },
    "aluminum": {
        "name": "Aluminum 6061-T6",
        "rho_293": 2.65e-8,
        "alpha": 0.004,
        "density": 2700,
        "cp": 900,
        "max_temp": 500     # Melts low
    },
    "copper": {
        "name": "Copper (OFHC)",
        "rho_293": 1.68e-8,
        "alpha": 0.0039,
        "density": 8960,
        "cp": 385,
        "max_temp": 600
    },
    "cuni7030": {
        "name": "CuNi 70/30",
        "rho_293": 3.75e-7,
        "alpha": 0.0004,
        "density": 8900,
        "cp": 377,
        "max_temp": 700
    }
}

# --- 2. HTS TAPE CONFIGURATION (Theva TLP AP) ---
class HTSTape:
    def __init__(self, width_mm=12.0, layers=2):
        self.width_mm = width_mm
        self.layers = layers
        # Theva TLP AP Spec @ 70K
        self.Ic_per_mm = 66.7 
        # Angular de-rating (Sin(20 deg) penalty approx 0.34)
        self.derating = 1.0 - math.sin(math.radians(20))
        
    @property
    def max_safe_current(self):
        # 80% Safety Margin
        total_Ic = self.Ic_per_mm * self.width_mm * self.layers * self.derating
        return 0.80 * total_Ic

# --- 3. LIM STAGE DEFINITIONS ---
class LIMStage:
    def __init__(self, name, tau_p, pitches, turns, tape_cfg, handoff_freq_hz):
        self.name = name
        self.tau_p = tau_p
        self.pitches = pitches
        self.turns = turns
        self.tape = tape_cfg
        self.handoff_freq = handoff_freq_hz
        
        # Geometry
        # L_LIM formula: (pitches * tau) + (2/3 * tau)
        self.length_active = (pitches * tau_p) + ((2/3) * tau_p)
        
        # Coil Width (slightly wider than plate to capture flux)
        self.w_coil = 0.7  # m (Plate is 0.5m)
        self.gap = 0.15    # m

# Define the 3 Stages
STAGES = [
    LIMStage(
        name="Low Velocity", 
        tau_p=20.0, 
        pitches=3, 
        turns=150, 
        tape_cfg=HTSTape(width_mm=12, layers=2),
        handoff_freq_hz=100.0
    ),
    LIMStage(
        name="Medium Velocity", 
        tau_p=60.0, 
        pitches=3, 
        turns=50, 
        tape_cfg=HTSTape(width_mm=12, layers=2),
        handoff_freq_hz=100.0
    ),
    LIMStage(
        name="High Velocity", 
        tau_p=200.0, 
        pitches=3, 
        turns=50, 
        tape_cfg=HTSTape(width_mm=12, layers=2),
        handoff_freq_hz=60.0 # Hand off at max velocity
    )
]

# --- 4. SLED & PLATE CONFIGURATION ---
class SledConfig:
    def __init__(self, material_key="gamma_titanium"):
        self.mat_props = MATERIALS[material_key]
        self.length = 5000.0   # 5 km
        self.height = 0.5      # Height in gap
        self.thickness = 0.2   # 200 mm
        self.mass_total = 8_372_000 # 8,372 Tonnes
        
        # Calculate derived mass of the plate itself (for thermal)
        vol_plate = self.length * self.height * self.thickness
        self.mass_plate = vol_plate * self.mat_props['density']

# --- 5. SYSTEM LIMITS ---
VOLTS_MAX = 100_000.0  # 100 kV
CRYO_PENALTY_RATIO = 50.0  # W_wall / W_cryo

# Calculate repeating unit length
L_UNIT_TOTAL = sum([s.length_active for s in STAGES])