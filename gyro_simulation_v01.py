#!/usr/bin/env python3
"""
Orbital Ring Rocket-to-Tether Handover Simulation

Simulates the fuel mass required to actively hold the orbital ring in Earth-synchronous 
precession using lateral rocket thrusters, followed by a smooth mechanical handover 
to the ground anchor tethers.

Evaluates a 1-meter slice of the ring at the equatorial crossing (maximum lateral load).
"""

import math
import matplotlib.pyplot as plt

# =============================================================================
# SYSTEM CONSTANTS & PARAMETERS
# =============================================================================
G0 = 9.80665
OMEGA_EARTH = 7.2921159e-5  # Earth rotation rate (rad/s)

# Ring Configuration (from prior momentum balance)
ALTITUDE = 250_000.0        # m
INCLINATION = 30.0          # degrees
V_CABLE = 8776.0            # m/s (Retrograde cable speed)
M_DRY_PER_METER = 100_000.0 # kg/m (Dry mass of casing + cable)

# Rocket Engine Specs (Assuming highly efficient Hydrolox engines)
I_SP = 450.0                # Specific Impulse (seconds)
V_EXHAUST = I_SP * G0       # Exhaust velocity (m/s)

# Handover Timeline (User Configurable)
T_HOLD_HOURS = 24.0         # Time rockets hold the ring while tethers drop & attach
T_TRANS_HOURS = 48.0        # Time over which thrust transfers to tethers

# =============================================================================
# PHYSICS CALCULATIONS
# =============================================================================
def run_handover_simulation():
    print("=== ORBITAL RING ROCKET HANDOVER SIMULATION ===")
    
    # 1. Peak Lateral Acceleration (Coriolis Precession Torque)
    # Required to force the gyroscope to precess with Earth's rotation
    a_lat_peak = 2.0 * OMEGA_EARTH * V_CABLE * math.sin(math.radians(INCLINATION))
    print(f"Required Lateral Accel: {a_lat_peak:.4f} m/s^2")
    
    # Timeline setup
    t_hold = T_HOLD_HOURS * 3600
    t_trans = T_TRANS_HOURS * 3600
    t_total = t_hold + t_trans
    dt = 60.0 # 1-minute time steps
    
    # Data logging
    history = {
        "time_hr": [], "mass_kg": [], "thrust_kN": [], 
        "tether_tension_kN": [], "rocket_fraction": []
    }
    
    # 2. Backward Integration of the Rocket Equation
    # We must integrate backwards from the END of the maneuver, because 
    # we know the final mass is exactly M_DRY.
    
    current_mass = M_DRY_PER_METER
    time_steps = int(t_total / dt)
    
    # Pre-allocate arrays for reverse-time calculation
    times = [i * dt for i in range(time_steps)]
    
    print(f"\nSimulating {T_HOLD_HOURS}h Hold + {T_TRANS_HOURS}h Transition...")
    
    for t in reversed(times):
        # Determine how much of the lateral load the rocket is taking
        if t <= t_hold:
            rocket_fraction = 1.0  # Rocket does 100% of the work
        else:
            # Linear ramp down during transition phase
            time_in_trans = t - t_hold
            rocket_fraction = 1.0 - (time_in_trans / t_trans)
            
        tether_fraction = 1.0 - rocket_fraction
        
        # Calculate forces based on CURRENT mass
        # F = m * a
        total_lateral_force = current_mass * a_lat_peak
        rocket_thrust = total_lateral_force * rocket_fraction
        tether_tension = total_lateral_force * tether_fraction
        
        # Calculate fuel burned in this time step
        # dm = (F / v_e) * dt
        dm = (rocket_thrust / V_EXHAUST) * dt
        
        # Because we are moving BACKWARDS in time, the mass at the previous 
        # moment must be LARGER by the amount of fuel burned.
        current_mass += dm
        
        # Log data (Pre-pending so it ends up in chronological order)
        history["time_hr"].insert(0, t / 3600)
        history["mass_kg"].insert(0, current_mass)
        history["thrust_kN"].insert(0, rocket_thrust / 1000)
        history["tether_tension_kN"].insert(0, tether_tension / 1000)
        history["rocket_fraction"].insert(0, rocket_fraction)

    # 3. Results Output
    initial_mass = history["mass_kg"][0]
    fuel_required = initial_mass - M_DRY_PER_METER
    mass_ratio = initial_mass / M_DRY_PER_METER
    
    print("\n--- RESULTS PER METER OF RING (At Equator) ---")
    print(f"Dry Mass of Ring:      {M_DRY_PER_METER:,.0f} kg/m")
    print(f"Fuel Required:         {fuel_required:,.0f} kg/m")
    print(f"Total Initial Mass:    {initial_mass:,.0f} kg/m")
    print(f"Gear Ratio (Wet/Dry):  {mass_ratio:.2f}")
    print(f"Max Rocket Thrust:     {history['thrust_kN'][0]:,.1f} kN/m")
    print(f"Final Tether Tension:  {history['tether_tension_kN'][-1]:,.1f} kN/m")
    
    plot_handover(history, T_HOLD_HOURS)

def plot_handover(data, t_hold_hr):
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Orbital Ring: Rocket-to-Tether Load Handover', fontsize=14)
    
    # Top Plot: Forces (Thrust vs Tension)
    axs[0].plot(data["time_hr"], data["thrust_kN"], 'r-', linewidth=2, label="Rocket Thrust")
    axs[0].plot(data["time_hr"], data["tether_tension_kN"], 'b-', linewidth=2, label="Tether Tension")
    axs[0].axvline(t_hold_hr, color='k', linestyle=':', label="Handover Begins")
    axs[0].set_ylabel("Force (kN per meter)")
    axs[0].set_title("Force Distribution Over Time")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()
    
    # Bottom Plot: Total Mass (The tyranny of the rocket equation)
    axs[1].plot(data["time_hr"], [m/1000 for m in data["mass_kg"]], 'g-', linewidth=2, label="Total Mass (Ring + Fuel)")
    axs[1].axhline(M_DRY_PER_METER/1000, color='k', linestyle='--', label="Dry Mass Limit")
    axs[1].axvline(t_hold_hr, color='k', linestyle=':')
    axs[1].set_xlabel("Time (Hours)")
    axs[1].set_ylabel("Mass (Tonnes per meter)")
    axs[1].set_title("System Mass Depletion")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_handover_simulation()