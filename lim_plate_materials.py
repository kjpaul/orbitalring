# Reaction plate material properties
# Options: "aluminum", "cuni7030", "titanium", "alpha_titanium" or "gamma_titanium"
PLATE_MATERIAL = "titanium"

# Reaction plate materials
# Aluminum properties
RHO_ALU_293K = 2.65e-8          # Electrical resistivity at 293 K (Ω·m)
RHO_ALU_70K = 4.853e-9          # Electrical resistivity at 70 K (Ω·m)
DENSITY_ALU = 2700              # Mass density (kg/m³)
C_P_ALU = 900                   # Specific heat capacity (J/kg·K)
K_ALU = 205                     # Thermal conductivity (W/m·K)
EM_ALU = 0.85                   # Emissivity (black anodized)
ALPHA_ALU = 3.663e-3            # Temperature coefficient of resistivity (1/K)

# CuNi 70/30 (Cupronickel) properties
RHO_CUNI_293K = 38e-8           # Electrical resistivity at 293 K (Ω·m) - nearly constant!
DENSITY_CUNI = 8900             # Mass density (kg/m³)
C_P_CUNI = 377                  # Specific heat capacity (J/kg·K)
K_CUNI = 29                     # Thermal conductivity (W/m·K)
EM_CUNI = 0.65                  # Emissivity (oxidized)
ALPHA_CUNI = 0.0004             # Temperature coefficient - very small!

# Pure Titanium properties (lunar-available from ilmenite)
RHO_TI_293K = 42e-8             # Electrical resistivity at 293 K (Ω·m)
DENSITY_TI = 4500               # Mass density (kg/m³)
C_P_TI = 520                    # Specific heat capacity (J/kg·K)
K_TI = 22                       # Thermal conductivity (W/m·K)
EM_TI = 0.60                    # Emissivity (oxidized)
ALPHA_TI = 0.0035               # Temperature coefficient (1/K)

# Alpha-2 titanium aluminide, α₂-Ti₃Al properties
RHO_aTI_293K = 50e-8             # Electrical resistivity at 293 K (Ω·m)
DENSITY_aTI = 4200               # Mass density (kg/m³)
C_P_aTI = 550                    # Specific heat capacity (J/kg·K)
K_aTI = 17                       # Thermal conductivity (W/m·K)
EM_aTI = 0.60                    # Emissivity (oxidized)
ALPHA_aTI = 0.0015               # Temperature coefficient (1/K)

# Gamma titanium aluminide, Ti-48Al-2Cr-2Nb properties
RHO_gTI_293K = 75e-8             # Electrical resistivity at 293 K (Ω·m)
DENSITY_gTI = 3900               # Mass density (kg/m³)
C_P_gTI = 570                    # Specific heat capacity (J/kg·K)
K_gTI = 15                       # Thermal conductivity (W/m·K)
EM_gTI = 0.60                    # Emissivity (oxidized)
ALPHA_gTI = 0.0012               # Temperature coefficient (1/K)

# Iron properties (for levitation plates)
DENSITY_IRON = 7870             # Mass density (kg/m³)


# Material-dependent properties
if PLATE_MATERIAL == "cuni7030":
    PLATE_DENSITY = DENSITY_CUNI
    PLATE_RHO_293K = RHO_CUNI_293K
    PLATE_ALPHA = ALPHA_CUNI
    PLATE_CP = C_P_CUNI
    PLATE_K = K_CUNI
    PLATE_EM = EM_CUNI
elif PLATE_MATERIAL == "titanium":
    PLATE_DENSITY = DENSITY_TI
    PLATE_RHO_293K = RHO_TI_293K
    PLATE_ALPHA = ALPHA_TI
    PLATE_CP = C_P_TI
    PLATE_K = K_TI
    PLATE_EM = EM_TI
elif PLATE_MATERIAL == "alpha_titanium":
    PLATE_DENSITY = DENSITY_aTI
    PLATE_RHO_293K = RHO_aTI_293K
    PLATE_ALPHA = ALPHA_aTI
    PLATE_CP = C_P_aTI
    PLATE_K = K_aTI
    PLATE_EM = EM_aTI
elif PLATE_MATERIAL == "gamma_titanium":
    PLATE_DENSITY = DENSITY_gTI
    PLATE_RHO_293K = RHO_gTI_293K
    PLATE_ALPHA = ALPHA_gTI
    PLATE_CP = C_P_gTI
    PLATE_K = K_gTI
    PLATE_EM = EM_gTI
else:  # aluminum
    PLATE_DENSITY = DENSITY_ALU
    PLATE_RHO_293K = RHO_ALU_293K
    PLATE_ALPHA = ALPHA_ALU
    PLATE_CP = C_P_ALU
    PLATE_K = K_ALU
    PLATE_EM = EM_ALU

