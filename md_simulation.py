"""
Mass Driver Simulation - Constraint Limited with Visualization
"""
import md_config as cfg
import md_physics as phys
import matplotlib.pyplot as plt
import numpy as np

def run_simulation():
    print("=== MASS DRIVER SIMULATION: CONSTRAINED ENGINEERING ===")
    print(f"Sled Mass: {cfg.SLED_MASS/1000:,.0f} tonnes")
    print(f"HTS Limit: {cfg.I_TARGET:.0f} A (80% of Ic)")
    print(f"Power Limit: {cfg.MAX_SITE_POWER/1e6} MW/site")
    print(f"Cryo Penalty: {cfg.CRYO_COP_PENALTY}:1")
    
    # --- Data Containers for Plotting ---
    history = {
        "time": [],
        "velocity": [],
        "accel_g": [],
        "thrust_mn": [],
        "current_ka": [],
        "power_total_mw": [],
        "power_mech_mw": [],
        "power_cryo_mw": [],
        "voltage_kv": [],
        "temp_plate": []
    }
    
    t = 0.0
    x = 0.0
    v = 0.0
    temp_plate = 300.0
    
    # Log interval (simulation steps between data capture)
    log_interval = 10 
    print_interval = 500 # Print less often than we plot
    
    print(f"{'T(h)':<8} {'Vel(km/s)':<10} {'Accel(g)':<10} {'Current(A)':<10} {'Cryo(MW)':<10} {'Volt(kV)':<8}")
    
    while v < cfg.V_LAUNCH_TARGET:
        
        # --- CONTROLLER LOGIC ---
        # 1. Target Slip Optimization (Model 1)
        # At low speed, we need high slip to induce current. 
        # As we speed up, we optimize for efficiency.
        v_slip_target = 40.0 
        
        # 2. Binary Search for Max Permissible Current
        i_min = 0.0
        i_max = cfg.I_TARGET
        i_op = 0.0
        
        # Placeholders for best valid state
        state_op = {}

        # Iteratively solve for max current that satisfies constraints
        for _ in range(12): 
            i_test = (i_min + i_max) / 2
            
            F, f_sup, I_ed, B = phys.calc_thrust_model1(v, v_slip_target, i_test, temp_plate)
            P_tot, P_cry, P_hys = phys.calc_power_balance(v, v_slip_target, i_test, F, f_sup)
            V_ind = phys.calc_voltage(i_test, f_sup)
            
            # Constraints: Max Power OR Max Voltage
            if P_tot > cfg.MAX_SITE_POWER or V_ind > cfg.VOLTS_MAX:
                i_max = i_test 
            else:
                i_min = i_test
                i_op = i_test
                
                # Store valid state
                state_op = {
                    "thrust": F,
                    "p_total": P_tot,
                    "p_cryo": P_cry,
                    "p_mech": F * v, # Mechanical power delivered to sled
                    "voltage": V_ind,
                    "f_supply": f_sup
                }

        # Fallback if no valid state found (should not happen with i_min=0)
        if not state_op:
            state_op = {"thrust": 0, "p_total": 0, "p_cryo": 0, "p_mech": 0, "voltage": 0}

        # --- UPDATE PHYSICAL STATE ---
        thrust = state_op["thrust"]
        accel = thrust / cfg.SLED_MASS
        
        # Thermal Update (Sled Plates)
        # Heat = Thrust * Slip Velocity
        q_plate = thrust * v_slip_target
        plate_vol = cfg.SLED_LENGTH * cfg.W_PLATE * cfg.T_PLATE
        m_plate = plate_vol * 3900 # Density TiAl
        
        # Temp Rise
        temp_plate += (q_plate * cfg.DT) / (m_plate * 600)
        # Temp Drop (Radiative)
        q_rad = 5.67e-8 * 0.8 * (cfg.SLED_LENGTH*cfg.W_PLATE*2) * (temp_plate**4 - 250**4)
        temp_plate -= (q_rad * cfg.DT) / (m_plate * 600)
        
        # Kinematics
        v += accel * cfg.DT
        x += v * cfg.DT
        t += cfg.DT
        
        # --- DATA LOGGING ---
        if int(t) % log_interval == 0:
            history["time"].append(t / 3600.0) # Hours
            history["velocity"].append(v / 1000.0) # km/s
            history["accel_g"].append(accel / 9.81)
            history["thrust_mn"].append(thrust / 1e6)
            history["current_ka"].append(i_op / 1000.0)
            history["power_total_mw"].append(state_op["p_total"] / 1e6)
            history["power_mech_mw"].append(state_op["p_mech"] / 1e6)
            history["power_cryo_mw"].append(state_op["p_cryo"] / 1e6)
            history["voltage_kv"].append(state_op["voltage"] / 1e3)
            history["temp_plate"].append(temp_plate)

        # Terminal Output
        if int(t) % print_interval == 0:
             print(f"{t/3600:<8.2f} {v/1000:<10.2f} {accel/9.81:<10.3f} {i_op:<10.0f} {state_op['p_cryo']/1e6:<10.2f} {state_op['voltage']/1e3:<8.1f}")

        # Safety & Timeout
        if t > 3600 * 24: # 24 hour timeout
            print("!!! TIMEOUT: Launch taking too long !!!") 
            break
            
    print("--- Launch Complete ---")
    print(f"Total Time: {t/3600:.2f} hours")
    
    # --- VISUALIZATION ---
    plot_results(history)

def plot_results(data):
    """
    Generates a 2x2 dashboard of launch metrics.
    """
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Mass Driver Launch Analysis (Constraint-Limited)', fontsize=16)
    
    # 1. VELOCITY & ACCELERATION vs TIME
    ax1 = axs[0, 0]
    ax1.set_title("Launch Profile")
    ax1.plot(data["time"], data["velocity"], 'b-', label="Velocity (km/s)")
    ax1.set_xlabel("Time (Hours)")
    ax1.set_ylabel("Velocity (km/s)", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True, alpha=0.3)
    
    ax2 = ax1.twinx()
    ax2.plot(data["time"], data["accel_g"], 'r--', label="Acceleration (g)")
    ax2.set_ylabel("Acceleration (g)", color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # 2. POWER BREAKDOWN vs VELOCITY (The "Cryo Wall")
    ax3 = axs[0, 1]
    ax3.set_title("Power Distribution vs. Velocity")
    # Stack plot or just overlapping lines
    ax3.plot(data["velocity"], data["power_total_mw"], 'k-', linewidth=2, label="Total Power Limit")
    ax3.plot(data["velocity"], data["power_cryo_mw"], 'r-', label="Cryo Penalty (Hysteresis)")
    ax3.plot(data["velocity"], data["power_mech_mw"], 'g-', label="Useful Mechanical Power")
    ax3.set_xlabel("Velocity (km/s)")
    ax3.set_ylabel("Power (MW per Site)")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.fill_between(data["velocity"], 0, data["power_cryo_mw"], color='red', alpha=0.1)

    # 3. CONSTRAINTS (Current & Voltage)
    ax4 = axs[1, 0]
    ax4.set_title("Engineering Constraints")
    ax4.plot(data["velocity"], data["current_ka"], 'm-', label="Current (kA)")
    ax4.set_xlabel("Velocity (km/s)")
    ax4.set_ylabel("Current (kA)", color='m')
    ax4.tick_params(axis='y', labelcolor='m')
    ax4.set_ylim(0, cfg.I_TARGET/1000 * 1.1)
    ax4.grid(True, alpha=0.3)
    
    ax5 = ax4.twinx()
    ax5.plot(data["velocity"], data["voltage_kv"], 'c--', label="Induced Voltage (kV)")
    ax5.set_ylabel("Voltage (kV)", color='c')
    ax5.tick_params(axis='y', labelcolor='c')
    ax5.axhline(cfg.VOLTS_MAX/1000, color='c', linestyle=':', label="Max Voltage")

    # 4. THERMAL (Plate Temp)
    ax6 = axs[1, 1]
    ax6.set_title("Sled Reaction Plate Temperature")
    ax6.plot(data["time"], data["temp_plate"], 'tab:orange', label="Plate Temp (K)")
    ax6.set_xlabel("Time (Hours)")
    ax6.set_ylabel("Temperature (K)")
    ax6.axhline(1100, color='red', linestyle='--', label="TiAl Softening Limit")
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_simulation()