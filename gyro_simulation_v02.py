#!/usr/bin/env python3
"""
Orbital Ring Aerodynamic Braking & Winch Payout Simulation

Models the lateral atmospheric drag on tethers and parachutes dropped from 
an orbital ring at 250 km altitude drifting laterally at 241 m/s.
Calculates the aerodynamic braking force and the kinetic energy injected 
into the global weather system.
"""

import math
import matplotlib.pyplot as plt

# =============================================================================
# SYSTEM CONSTANTS
# =============================================================================
ALTITUDE_RING = 250_000.0   # m
V_LATERAL_DRIFT = 241.0     # m/s (Max lateral drift at equator for 30 deg inclination)
NUM_TETHERS = 800           # Number of anchor lines deployed globally

# Tether Specs (Assuming Kevlar/Zylon macro-cables)
TETHER_DIAMETER = 0.5       # meters (thick enough to survive and hold load)
C_D_CYLINDER = 1.2          # Drag coefficient of a long cylinder

# Parachute / Drogue Specs
DROGUE_DIAMETER = 50.0      # meters (massive supersonic drogue chutes)
C_D_PARACHUTE = 1.5         # Drag coefficient of parachute

# =============================================================================
# ATMOSPHERIC MODEL (Simplified Exponential)
# =============================================================================
def get_air_density(altitude_m):
    """Standard atmosphere exponential density model."""
    rho_0 = 1.225 # kg/m^3 at sea level
    scale_height = 8500.0 # meters
    if altitude_m > 150_000:
        return 0.0 # Negligible drag above 150km
    return rho_0 * math.exp(-altitude_m / scale_height)

# =============================================================================
# SIMULATION
# =============================================================================
def run_aerobraking_sim():
    print("=== ORBITAL RING AERODYNAMIC BRAKING SIMULATION ===")
    
    # Divide the 250km tether into 100-meter segments to integrate drag
    segment_length = 100.0
    altitudes = [h for h in range(0, int(ALTITUDE_RING), int(segment_length))]
    
    total_tether_drag_N = 0.0
    drag_profile = []
    
    print("Integrating atmospheric drag along 250km tether...")
    
    # 1. Calculate drag on the tether cable itself
    for h in altitudes:
        rho = get_air_density(h)
        # Drag equation: F = 0.5 * rho * v^2 * Cd * A
        area = TETHER_DIAMETER * segment_length
        drag = 0.5 * rho * (V_LATERAL_DRIFT**2) * C_D_CYLINDER * area
        total_tether_drag_N += drag
        if h % 5000 == 0:
            drag_profile.append((h/1000, drag))
            
    # 2. Calculate drag from the Drogue Parachute (deployed at 2,000m altitude)
    rho_drogue = get_air_density(2000.0)
    area_drogue = math.pi * (DROGUE_DIAMETER / 2.0)**2
    parachute_drag_N = 0.5 * rho_drogue * (V_LATERAL_DRIFT**2) * C_D_PARACHUTE * area_drogue
    
    total_force_per_station_N = total_tether_drag_N + parachute_drag_N
    global_aerobraking_force = total_force_per_station_N * NUM_TETHERS
    
    # 3. Weather Impact (Power Dissipation = Force * Velocity)
    power_per_station_W = total_force_per_station_N * V_LATERAL_DRIFT
    global_weather_power_GW = (global_aerobraking_force * V_LATERAL_DRIFT) / 1e9
    
    print("\n--- SINGLE TETHER RESULTS ---")
    print(f"Cable Drag Force:     {total_tether_drag_N/1000:,.1f} kN")
    print(f"Parachute Drag Force: {parachute_drag_N/1000:,.1f} kN")
    print(f"Total Braking Force:  {total_force_per_station_N/1000:,.1f} kN per station")
    print(f"Heat/Turbulence Gen:  {power_per_station_W/1e6:,.1f} MW per station")
    
    print("\n--- GLOBAL SYSTEM RESULTS ---")
    print(f"Total Braking Force:  {global_aerobraking_force/1e6:,.1f} MN (MegaNewtons)")
    print(f"Weather Injection:    {global_weather_power_GW:,.1f} GW (GigaWatts) of atmospheric turbulence")
    
    # Plotting the Drag Profile
    alts, drags = zip(*drag_profile)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(drags, alts, 'b-', linewidth=2)
    ax.set_title("Aerodynamic Drag Profile of a Single 250km Tether")
    ax.set_xlabel("Drag Force per 100m segment (Newtons)")
    ax.set_ylabel("Altitude (km)")
    ax.grid(True, alpha=0.3)
    ax.fill_betweenx(alts, 0, drags, color='blue', alpha=0.1)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_aerobraking_sim()