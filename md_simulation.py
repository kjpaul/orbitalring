"""
Mass Driver Simulation Main Loop
Simulates the launch of an 8,372-tonne sled.
"""
import md_config as cfg
import md_physics as phys
import math

def run_launch_simulation():
    print(f"--- MASS DRIVER LAUNCH SIMULATION ---")
    print(f"Payload Class: {cfg.SLED_MASS_TOTAL/1000} tonnes")
    print(f"Target Velocity: {cfg.V_LAUNCH_TARGET} m/s")
    print(f"Accel Limit: {cfg.G_FORCE_LIMIT:.2f} m/s^2")
    print(f"Material: {cfg.PLATE_MATERIAL} (High Resistivity)\n")
    
    t = 0
    x = 0
    v = 0
    temp_plate = 300.0 # Kelvin
    
    # Header
    print(f"{'Time(s)':<8} {'Dist(km)':<10} {'Vel(m/s)':<10} {'Accel(g)':<10} {'Thrust(MN)':<12} {'Slip(m/s)':<10} {'Temp(K)':<8} {'G-Factor':<8}")
    
    while v < cfg.V_LAUNCH_TARGET:
        # 1. Control Loop: Determine desired Thrust
        # We want constant acceleration of G_FORCE_LIMIT
        req_force = cfg.SLED_MASS_TOTAL * cfg.G_FORCE_LIMIT
        
        # 2. Physics Step
        # Optimize slip for max thrust if needed, or regulate for efficiency
        # At low speed, we need high slip to get any induction in TiAl
        if v < 100:
            current_slip = 100.0
        else:
            current_slip = v * cfg.TARGET_SLIP_RATIO
            current_slip = max(cfg.V_SLIP_MIN, min(cfg.V_SLIP_MAX, current_slip))
            
        avail_thrust, G, B_field = phys.calc_thrust_force(v, current_slip, temp_plate)
        
        # 3. Apply Limits (Throttle)
        if avail_thrust > req_force:
            applied_thrust = req_force
            throttle = req_force / avail_thrust
        else:
            applied_thrust = avail_thrust
            throttle = 1.0
            
        accel = applied_thrust / cfg.SLED_MASS_TOTAL
        
        # 4. Thermal Update
        # Heat = Power * Slip_Fraction roughly, or F * v_slip
        # Heat goes into the TiAl plates
        q_heat = applied_thrust * current_slip * throttle
        plate_volume = cfg.SLED_LENGTH * cfg.PLATE_WIDTH * cfg.PLATE_THICKNESS
        plate_mass = plate_volume * 3900 # Density of TiAl approx 3900 kg/m3
        c_p = 600 # Specific heat approx 600 J/kgK
        
        temp_rise = (q_heat * cfg.DT) / (plate_mass * c_p)
        
        # Radiative cooling (Stefan-Boltzmann)
        area_rad = cfg.SLED_LENGTH * cfg.PLATE_WIDTH * 2 # Top and bottom
        q_rad = phys.STEFAN_BOLTZMANN * 0.8 * area_rad * (temp_plate**4 - 200**4)
        temp_drop = (q_rad * cfg.DT) / (plate_mass * c_p)
        
        temp_plate += (temp_rise - temp_drop)
        
        # 5. Integration
        v += accel * cfg.DT
        x += v * cfg.DT
        t += cfg.DT
        
        # Log every 10 seconds or significant events
        if int(t) % 10 == 0 and (t - int(t) < cfg.DT):
             print(f"{t:<8.1f} {x/1000:<10.1f} {v:<10.1f} {accel/9.81:<10.2f} {applied_thrust/1e6:<12.1f} {current_slip:<10.1f} {int(temp_plate):<8} {G:<8.2f}")

        # Safety Cutoff
        if temp_plate > 1200:
            print("!!! ABORT: PLATE OVERHEATING !!!")
            break
            
    print(f"\nLaunch Complete.")
    print(f"Final Velocity: {v:.2f} m/s")
    print(f"Distance Covered: {x/1000:.2f} km")
    
if __name__ == "__main__":
    run_launch_simulation()