"""
Mass Driver Simulation - Constraint Limited
"""
import md_config as cfg
import md_physics as phys

def run_simulation():
    print("=== MASS DRIVER SIMULATION: CONSTRAINED ENGINEERING ===")
    print(f"Sled Mass: {cfg.SLED_MASS/1000} tonnes")
    print(f"HTS Limit: {cfg.I_TARGET:.1f} A (80% of Ic)")
    print(f"Power Limit: {cfg.MAX_SITE_POWER/1e6} MW/site")
    print(f"Cryo Penalty: {cfg.CRYO_COP_PENALTY}:1")
    
    t = 0.0
    x = 0.0
    v = 0.0
    temp_plate = 300.0
    
    # We step every 10s for output log
    log_interval = 100 
    
    print(f"{'T(s)':<8} {'Vel(m/s)':<10} {'Accel(g)':<10} {'Current(A)':<10} {'Thrust(MN)':<10} {'Power(MW)':<10} {'Cryo(MW)':<10} {'Volt(kV)':<8}")
    
    while v < cfg.V_LAUNCH_TARGET:
        
        # --- CONTROLLER LOGIC ---
        # We try to push max current, but respect constraints.
        
        # 1. Determine optimal slip (Model 1 implies slip ~ R/X peak)
        # For Model 1, peak thrust often occurs at specific f_slip.
        # We'll use a constant v_slip target that is efficient.
        v_slip_target = 40.0 # m/s (Tuned for this Tau_p/Resistance)
        
        # 2. Binary Search for Max Permissible Current
        # We need to find the highest I_peak that satisfies V < Vmax and P < Pmax
        i_min = 0.0
        i_max = cfg.I_TARGET
        i_op = 0.0
        
        # Physics results placeholders
        thrust_op = 0
        p_total_op = 0
        p_cryo_op = 0
        v_induced_op = 0
        
        # Iteratively solve for max current
        for _ in range(10): 
            i_test = (i_min + i_max) / 2
            
            F, f_sup, I_ed, B = phys.calc_thrust_model1(v, v_slip_target, i_test, temp_plate)
            P_tot, P_cry, P_hys = phys.calc_power_balance(v, v_slip_target, i_test, F, f_sup)
            V_ind = phys.calc_voltage(i_test, f_sup)
            
            # Check constraints
            if P_tot > cfg.MAX_SITE_POWER or V_ind > cfg.VOLTS_MAX:
                i_max = i_test # Too high, lower ceiling
            else:
                i_min = i_test # Safe, try raising floor
                i_op = i_test
                thrust_op = F
                p_total_op = P_tot
                p_cryo_op = P_cry
                v_induced_op = V_ind

        # --- UPDATE STATE ---
        accel = thrust_op / cfg.SLED_MASS
        
        # Thermal update (Sled plates)
        # Heat into plate = Slip Power = Thrust * v_slip
        q_plate = thrust_op * v_slip_target
        plate_vol = cfg.SLED_LENGTH * cfg.W_PLATE * cfg.T_PLATE
        m_plate = plate_vol * 3900 # TiAl density
        temp_plate += (q_plate * cfg.DT) / (m_plate * 600) # cp approx
        
        # Radiative cooling
        q_rad = 5.67e-8 * 0.8 * (cfg.SLED_LENGTH*cfg.W_PLATE*2) * (temp_plate**4 - 250**4)
        temp_plate -= (q_rad * cfg.DT) / (m_plate * 600)
        
        v += accel * cfg.DT
        x += v * cfg.DT
        t += cfg.DT
        
        if int(t) % log_interval == 0:
             print(f"{t:<8.0f} {v:<10.1f} {accel/9.81:<10.3f} {i_op:<10.0f} {thrust_op/1e6:<10.2f} {p_total_op/1e6:<10.2f} {p_cryo_op/1e6:<10.2f} {v_induced_op/1e3:<8.1f}")

        # Safety
        if t > 50000: # Timeout
            print("timeout") 
            break
            
    print("--- Launch Complete ---")
    print(f"Total Time: {t/3600:.2f} hours")
    print(f"Final Distance: {x/1000:.0f} km")

if __name__ == "__main__":
    run_simulation()