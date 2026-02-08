"""
Mass Driver Simulation - 3g Human Rated Profile
"""
import md_config as cfg
import md_physics as phys
import matplotlib.pyplot as plt

def run_simulation():
    print(f"=== ORBITAL RING MASS DRIVER: {cfg.LAUNCH_PROFILE.upper()} PROFILE ===")
    print(f"Ring Spare Capacity: {cfg.F_MAX_LIFT_M/1e3:.1f} kN/m")
    print(f"Sled Mass Limit: {cfg.SLED_MASS_TOTAL/1000:,.0f} tonnes")
    print(f"V_Launch Target: {cfg.V_LAUNCH_TARGET:.0f} m/s")
    print(f"Power Limit: {cfg.MAX_SITE_POWER/1e6} MW/site")
    
    # Data logging setup (same as before)
    history = {"time": [], "velocity": [], "accel_g": [], "power_cryo_mw": [], "current_ka": []}
    
    t, x, v = 0.0, 0.0, 0.0
    temp_plate = 300.0
    
    print(f"{'T(h)':<8} {'Vel(km/s)':<10} {'Accel(g)':<10} {'Current(A)':<10} {'Cryo(MW)':<10} {'Volt(kV)':<8}")
    
    while v < cfg.V_LAUNCH_TARGET:
        v_slip_target = 40.0 
        
        # Constraint Solver
        i_min, i_max = 0.0, cfg.I_TARGET
        i_op = 0.0
        state_op = {}

        for _ in range(12): 
            i_test = (i_min + i_max) / 2
            F, f_sup, I_ed, B = phys.calc_thrust_model1(v, v_slip_target, i_test, temp_plate)
            P_tot, P_cry, P_hys = phys.calc_power_balance(v, v_slip_target, i_test, F, f_sup)
            V_ind = phys.calc_voltage(i_test, f_sup)
            
            if P_tot > cfg.MAX_SITE_POWER or V_ind > cfg.VOLTS_MAX:
                i_max = i_test 
            else:
                i_min = i_test
                i_op = i_test
                state_op = {"thrust": F, "p_cryo": P_cry, "voltage": V_ind}

        if not state_op: state_op = {"thrust": 0, "p_cryo": 0, "voltage": 0}

        thrust = state_op["thrust"]
        accel = thrust / cfg.SLED_MASS_TOTAL
        
        # Thermal & Kinematics
        q_plate = thrust * v_slip_target
        m_plate = (cfg.SLED_LENGTH * cfg.W_PLATE * cfg.T_PLATE) * 3900
        temp_plate += (q_plate * cfg.DT) / (m_plate * 600)
        q_rad = 5.67e-8 * 0.8 * (cfg.SLED_LENGTH*cfg.W_PLATE*2) * (temp_plate**4 - 250**4)
        temp_plate -= (q_rad * cfg.DT) / (m_plate * 600)
        
        v += accel * cfg.DT
        x += v * cfg.DT
        t += cfg.DT
        
        if int(t) % 20 == 0:
            history["time"].append(t/3600)
            history["velocity"].append(v/1000)
            history["accel_g"].append(accel/9.81)
            history["power_cryo_mw"].append(state_op["p_cryo"]/1e6)
            history["current_ka"].append(i_op/1000)

        if int(t) % 1000 == 0:
             print(f"{t/3600:<8.2f} {v/1000:<10.2f} {accel/9.81:<10.3f} {i_op:<10.0f} {state_op['p_cryo']/1e6:<10.2f} {state_op['voltage']/1e3:<8.1f}")
            
        if t > 3600 * 48: break

    print(f"Launch Complete. Time: {t/3600:.2f} h. Dist: {x/1000:,.0f} km")
    
    # Basic Plot
    fig, ax1 = plt.subplots()
    ax1.plot(history["time"], history["velocity"], 'b-', label="Velocity")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("Velocity (km/s)", color='b')
    ax2 = ax1.twinx()
    ax2.plot(history["time"], history["power_cryo_mw"], 'r--', label="Cryo Load")
    ax2.set_ylabel("Cryo Load (MW)", color='r')
    plt.title(f"3g Profile Launch: {cfg.SLED_MASS_TOTAL/1000:.0f}t Sled")
    plt.show()

if __name__ == "__main__":
    run_simulation()