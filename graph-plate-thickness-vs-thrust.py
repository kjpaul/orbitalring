import math
import matplotlib.pyplot as plt
import numpy as np

# Parameters from user's data and prompt
MU0 = 4 * math.pi * 1e-7
TAU_P = 100.0  # m
W_COIL = 1.0   # m
L_ACTIVE = 3*TAU_P # m (assumed for scaling)
F_SLIP = 0.5  # Hz
B = 0.03104        # Tesla (assumed for scaling)
RHO = 420.0e-9 # Ohm-m (Al at 150K)
ALPHA = 0.0035 
ALU  = [26.5, 0.0039, "Pure Aluminum", 2700] # [resistivity in nano ohm meters, alpha in 1/K, name, density in kg/m^3]
TI   = [420, 0.0035, "Pure Titanium", 4500]
TI_ALPHA = [500, 0.0015, "Alpha-2 titanium aluminide, α₂-Ti₃Al", 4200]
TI_GAMMA = [750, 0.0012, "Gamma titanium aluminide, Ti-48Al-2Cr-2Nb", 3900]
TI_G5 = [1710, 0.0010, "Grade 5 Titanium", 4430]
CuNi = [380, 0.0004, "Cupronickel, C71500, Marine alloy, 70/30 alloy", 8900]
CuNi55 = [380, 0.00001, "Eureka, CuNi 55/45", 8900]
MG   = [44, 0.0040, "Magnesium", 1740]
CA   = [34, 0.0035, "Calcium", 1526]

TEMP = 150
METALS = [ALU, TI, TI_ALPHA, TI_GAMMA, TI_G5, CuNi, CuNi55, MG, CA]
#METALS = [ALU]

def get_rho (metal=ALU, temp=TEMP):
    if len(metal) < 2:
        print ("Invalid metal.")
        return
    rho20 = metal[0] * 1e-9
    alpha = metal[1]
    return rho20 * (1 + alpha*(temp - 293))

def calc_loop_inductance(tau_p, w_coil):
    return (MU0 * tau_p / math.pi) * math.log(tau_p / w_coil)

def model1_thrust(d_plate, rho, f_slip, b, tau_p, w_coil, l_active):
    v_slip = f_slip * 2 * tau_p
    L = calc_loop_inductance(tau_p, w_coil)
    omega = 2 * math.pi * f_slip
    X = omega * L
    
    # Resistance of the loop
    R = rho * tau_p / (d_plate * w_coil)
    Z_sq = R**2 + X**2
    
    EMF = (4.0 / math.pi) * v_slip * b * w_coil
    I_loop = EMF / math.sqrt(Z_sq)
    N_loops = l_active / tau_p
    P = N_loops * (I_loop**2) * R
    F = P / v_slip
    return F, R, X

def makeGraph(metal):
    if len(metal) < 2:
        print ("Invalid metal.")
        return

    d_range = np.linspace(0.001, 0.2, 200) # 1mm to 200mm
    rho = get_rho(metal, TEMP)
    thrusts = [model1_thrust(d, rho, F_SLIP, B, TAU_P, W_COIL, L_ACTIVE)[0] for d in d_range]
    _, R_100, X = model1_thrust(0.15, rho, F_SLIP, B, TAU_P, W_COIL, L_ACTIVE)

    xr = X/R_100
    thr_ratio = 2*xr/(1 + xr**2)

    i = 0
    max_thrust = [0, 0]
    for force in thrusts:
        i = i + 1
        if force > max_thrust[0]:
            max_thrust[0] = force
            max_thrust[1] = i
    

    print("-"*70)
    print(f"max thrust: {max_thrust[0]} at {max_thrust[1]}")
    print(f"{metal[2]} at 100mm: R = {R_100:.2e}, X = {X:.2e}, PF = {thr_ratio*100:.2f}%")
    print(f"Ratio X/R: {X/R_100:.2f}")

    plt.rc('font', size=22)
    plt.plot(d_range * 1000, np.array(thrusts)/1000)
    plt.axvline(100, color='r', linestyle='--', label='Plate thickness = 100mm')
    plt.axvline(max_thrust[1], color='b', linestyle='--', label=f'Max point    = {max_thrust[1]}mm\nThrust ratio = {thr_ratio*100:.2f}%')
    plt.xlabel('Plate Thickness (mm)')
    plt.ylabel('Thrust (kN)')
    plt.title(f'Thrust vs Thickness (Model 1) for {metal[2]}')
    plt.legend()
    plt.grid(True)
    plt.savefig('thrust_vs_thickness.png')
    plt.show()

#ax.annotate('Critical Point', (3, 4.5), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14, color='red')
#ax.annotate('Low', (3, 0.5), textcoords="offset points", xytext=(0,0), ha='center', fontsize=10, color='blue')


for metal in METALS:
    makeGraph(metal)

print("-"*70)