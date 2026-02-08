import numpy as np
import matplotlib.pyplot as plt
import md_config as cfg
import md_physics as phys

def run_simulation():
    # --- SETUP ---
    sled = cfg.SledConfig(material_key="gamma_titanium") # Change material here
    
    # Calculate how many of each stage fit under the 5km plate
    # Total Unit Length = L_low + L_mid + L_high
    # Number of active units = Plate_Length / Unit_Length
    n_active_units = sled.length / cfg.L_UNIT_TOTAL
    
    print(f"=== MASS DRIVER SIMULATION ===")
    print(f"Plate Material: {sled.mat_props['name']}")
    print(f"Plate Length: {sled.length} m")
    print(f"Repeating Unit Length: {cfg.L_UNIT_TOTAL:.2f} m")
    print(f"Active Modules under Sled: {n_active_units:.2f}")
    
    # State Variables
    t = 0
    v = 0
    x = 0
    temp = 300.0 # K
    dt = 1.0 # s
    
    # Data Logs
    history = {
        "time": [], "velocity": [], "dist": [], "temp": [],
        "thrust_total": [], "power_total": [],
        # Per Stage Logs
        "s1_current": [], "s1_volts": [], "s1_thrust": [],
        "s2_current": [], "s2_volts": [], "s2_thrust": [],
        "s3_current": [], "s3_volts": [], "s3_thrust": []
    }
    
    target_v = 12000 # m/s (example target)
    
    print(f"{'Time':<9} {'Vel(m/s)':<10} {'Amps (A)':<10} {'Volts (kV)':<10} {'Thrust(N)':<12} {'Temp(K)':<8} {'Stage':<10}")
    
    while v < target_v:
        # Determine active stage based on velocity/frequency efficiency
        # We simulate all 3, but the physics engine returns 0 if f > 60Hz
        
        # Optimal Slip Strategy:
        # At low speed, we need slip ~ Tau/TimeConstant. 
        # Simple heuristic: v_slip = 10% of v_sync or capped
        v_slip_target = 50.0 # m/s constant slip target usually works well for LIMs
        
        total_thrust = 0
        total_power = 0
        stage_data = [] # To store [I, V, F] for each stage
        
        # --- CALCULATE FOR EACH STAGE TYPE ---
        for stage in cfg.STAGES:
            F, I, V, P_wall, P_mech, P_cryo = phys.solve_stage_performance(
                stage, v, v_slip_target, temp, sled
            )
            
            # Scale by number of active units under the sled
            F_total_stage = F * n_active_units
            P_wall_total_stage = P_wall * n_active_units
            
            total_thrust += F_total_stage
            total_power += P_wall_total_stage
            
            stage_data.append((I, V, F_total_stage))
            
        # --- KINEMATICS ---
        accel = total_thrust / sled.mass_total
        v += accel * dt
        x += v * dt
        t += dt
        
        # --- THERMAL ---
        # Heat = Total Slip Power = Thrust * Slip_Velocity
        # Note: This assumes all slip energy goes into the plate (conservative)
        Q_heat = total_thrust * v_slip_target * dt
        
        # Adiabatic rise
        temp_rise = Q_heat / (sled.mass_plate * sled.mat_props['cp'])
        
        # Radiative Cooling
        # Area = 2 * Length * Height (Sides) + 2 * Length * Thickness (Top/Bot)
        area_rad = 2 * sled.length * (sled.height + sled.thickness)
        Q_rad = 5.67e-8 * 0.85 * area_rad * (temp**4 - 250**4) * dt
        temp_drop = Q_rad / (sled.mass_plate * sled.mat_props['cp'])
        
        temp += (temp_rise - temp_drop)
        
        # --- LOGGING ---
        if t % 5 < dt: # Log every ~5 sim seconds
            history["time"].append(t)
            history["velocity"].append(v)
            history["dist"].append(x/1000)
            history["temp"].append(temp)
            history["thrust_total"].append(total_thrust/1e6)
            history["power_total"].append(total_power/1e6)
            
            history["s1_current"].append(stage_data[0][0])
            history["s1_volts"].append(stage_data[0][1]/1000)
            history["s1_thrust"].append(stage_data[0][2]/1e6)
            
            history["s2_current"].append(stage_data[1][0])
            history["s2_volts"].append(stage_data[1][1]/1000)
            history["s2_thrust"].append(stage_data[1][2]/1e6)
            
            history["s3_current"].append(stage_data[2][0])
            history["s3_volts"].append(stage_data[2][1]/1000)
            history["s3_thrust"].append(stage_data[2][2]/1e6)
            
        if t % 100000 < dt:
            dom_stage = "None"
            if stage_data[0][2] > 100: dom_stage = "Low"
            elif stage_data[1][2] > 100: dom_stage = "Mid"
            elif stage_data[2][2] > 100: dom_stage = "High"
            
            print(f"{t/3600:<6.1f} hr {v:<6.1f} m/s {I:<8.1f} A {V/1000:<7.1f} kV {total_thrust:<10.2f} N {temp:<6.1f} K {dom_stage}")
            
        if temp > sled.mat_props['max_temp']:
            print("!!! PLATE MELTED !!!")
            break
            
    plot_dashboard(history)

def plot_dashboard(data):
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Segmented Mass Driver Performance', fontsize=16)
    
    # 1. Velocity & Thrust
    ax1 = axs[0,0]
    ax1.set_title("Kinematics")
    ax1.plot(data['time'], data['velocity'], 'k-', label='Velocity (m/s)')
    ax1_r = ax1.twinx()
    ax1_r.plot(data['time'], data['thrust_total'], 'r-', label='Total Thrust (MN)', alpha=0.6)
    ax1.set_ylabel("Velocity (m/s)")
    ax1_r.set_ylabel("Thrust (MN)")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_r.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc=2)
    
    # 2. Stage Handoff (Thrust per Stage)
    ax2 = axs[0,1]
    ax2.set_title("Stage Handoff (Thrust Contribution)")
    ax2.plot(data['velocity'], data['s1_thrust'], label='Low (Tau=3m)')
    ax2.plot(data['velocity'], data['s2_thrust'], label='Mid (Tau=30m)')
    ax2.plot(data['velocity'], data['s3_thrust'], label='High (Tau=200m)')
    ax2.set_xlabel("Sled Velocity (m/s)")
    ax2.set_ylabel("Thrust (MN)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Electrical (Voltage)
    ax3 = axs[1,0]
    ax3.set_title("Stator Voltage per Stage")
    ax3.plot(data['velocity'], data['s1_volts'], '--', label='Low V')
    ax3.plot(data['velocity'], data['s2_volts'], '--', label='Mid V')
    ax3.plot(data['velocity'], data['s3_volts'], '--', label='High V')
    ax3.axhline(cfg.VOLTS_MAX/1000, color='r', linestyle=':', label='Limit (100kV)')
    ax3.set_xlabel("Sled Velocity (m/s)")
    ax3.set_ylabel("Voltage (kV)")
    ax3.legend()
    
    # 4. Thermal
    ax4 = axs[1,1]
    ax4.set_title("Reaction Plate Temperature")
    ax4.plot(data['time'], data['temp'], 'r-')
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Temp (K)")
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_simulation()