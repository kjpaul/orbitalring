import math

# Constants and setup
g0 = 9.81                   # m/s^2 at Earth's surface
Re = 6371000.0             # Earth radius (m)
rho0 = 1.225               # kg/m^3 at sea level
H = 8400.0                 # scale height (m)
Cd = 1.3
A = 707.0                  # m^2 r = 15 m
m = 10000.0                # kg
beta = m / (Cd*A)

# Initial conditions
h = 300000.0               # starting altitude in m (300 km)
v = 0.0                    # initial downward velocity (m/s)
t = 0.0                    # time (s)
dt = 0.1                   # time step (s)

max_g_load = 0.0
max_v = 0.0
data = []
max_g_h = 0.0
max_v_h = 0.0
max_drag = 0.0
max_drag_h = 0.0
max_drag_g_local = 0.0
max_drag_a = 0.0
max_drag_v = 0.0


while h > 0.0:
    # Calculate local gravity
    r = Re + h
    g_local = g0 * (Re / r)**2
    
    # Calculate local air density
    rho = rho0 * math.exp(-h / H) if h < 1.0e6 else 0.0  # no atmosphere above ~1,000 km
    
    # Drag acceleration
    # If v>0, drag is upward => subtract from gravity
    drag_acc = 0.5*rho*(v**2)/beta
    
    # Net acceleration (downward = +)
    a = g_local - drag_acc
    
    # Update velocity and altitude
    v_new = v + a*dt
    h_new = h - 0.5*(v + v_new)*dt  # average velocity * dt
    
    # Get max velocity
    if v_new > max_v:
        max_v = v_new
        max_v_h = h_new

    # Get max drag
    if drag_acc > max_drag:
        max_drag = drag_acc
        max_drag_h = h_new
        max_drag_g_local = g_local
        max_drag_a = a
        max_drag_v = v_new

    
    # Log g-load in multiples of g0
    g_load = a / g0
    if abs(g_load) > abs(max_g_load):
        max_g_load = g_load
        max_g_h = h_new
    
    # Update states
    v = v_new
    h = h_new
    t += dt
    
    data.append((t, h, v, a, g_load))
    

# Print summary
print(f"Time to reach ground: {t/60:.1f} minutes")
print(f"Impact velocity (m/s): {v:.1f}")
print(f"Peak g-load: {max_g_load:.2f} g")
print(f"Peak g-load h: {max_g_h/1000:.2f} km")
print(f"Max velocity: {max_v:.2f} m/s")
print(f"Max velocity h: {max_v_h/1000:.2f} km")
print(f"Max drag: {max_drag:.2f} m/s^2")
print(f"Max drag h: {max_drag_h/1000:.2f} km")
print(f"Max drag g_local: {max_drag_g_local:.2f} m/s^2")
print(f"Max drag a: {max_drag_a:.2f} m/s^2")
print(f"Max drag v: {max_drag_v:.2f} m/s")
