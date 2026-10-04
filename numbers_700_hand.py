"""Hand calculations for the design point in ring_altitude.py (the values the 800 km sheet did by hand)."""
import math, ring_altitude as R
w = R.OMEGA_SIDEREAL; Re = R.R_EARTH_EQ; r = R.R_ORBIT
vg = R.V_GROUND_SYNC
print(f"altitude {R.ALTITUDE_KM:.0f} km, array {R.ARRAY_ALTITUDE/1e3:.0f} km")
a_rel = R.V_CABLE**2/r; t_gap = math.sqrt(0.2/a_rel)
print(f"relative radial acceleration at a break: {R.V_CABLE**2/r - R.G_LOCAL:.3f} + {R.G_LOCAL:.3f} = {a_rel:.3f} m/s2; gap time sqrt(0.2/{a_rel:.3f}) = {t_gap:.3f} s; contact speed {a_rel*t_gap:.2f} m/s")
I = R.CABLE_SIDE**4/12; print(f"second moment {R.CABLE_SIDE:.2f}^4/12 = {I:.0f} m4; L_char = sqrt(500e9 x {I:.0f} / {R.CABLE_TENSION:.4g}) = {math.sqrt(500e9*I/R.CABLE_TENSION):.1f} m")
for i in (5, 10, 30): print(f"drift at {i} deg: {vg:.1f} x sin = {vg*math.sin(math.radians(i)):.0f} m/s")
M = R.M_RING_M*R.L_RING; print(f"lateral KE at 30 deg: 0.5 x {M:.3e} x {vg*0.5:.1f}^2 x 0.5 = {0.5*M*(vg*0.5)**2*0.5:.3e} J")
print(f"rotation radius at 30 deg apex {r*math.cos(math.radians(30))/1e3:.0f} km, casing speed there {vg*math.cos(math.radians(30)):.1f} m/s")
p = R.M_CABLE_M*R.V_CABLE - R.M_LOAD_M*vg; W = R.M_LOAD_M*R.G_NET
print(f"net momentum per metre p = {R.M_CABLE_M:.0f} x {R.V_CABLE:.1f} - 12000 x {vg:.1f} = {p/1e6:.1f} MN s/m; casing weight {W/1e3:.1f} kN/m")
for i in (5, 10, 20, 30):
    F = 2*w*p*math.sin(math.radians(i)); print(f"  Table 2.1, {i} deg: F = 2 w p sin i = {F/1e3:.1f} kN/m = {100*F/W:.0f}% of vertical load")
ang = math.asin(r*math.sin(math.radians(45))/Re); reach = math.degrees(ang) - 45; print(f"latitude reach of a 45 deg line: +-{reach:.1f} deg, length {Re*math.sin(math.radians(reach))/math.sin(math.radians(45))/1e3:.0f} km")
for h in (R.ALTITUDE, R.ARRAY_ALTITUDE): print(f"horizon from {h/1e3:.0f} km: {math.degrees(math.acos(Re/(Re+h))):.1f} deg")
t = 2*math.sqrt(R.ALTITUDE/0.981); print(f"anchor transit at 0.1 g: 2 sqrt({R.ALTITUDE:.0f}/0.981) = {t:.0f} s = {t/60:.1f} min, peak {0.981*t/2:.0f} m/s; at 100 to 130 km/h: {R.ALTITUDE_KM/130:.1f} to {R.ALTITUDE_KM/100:.1f} h")
ra = R.R_ARRAY; print(f"array sunlit share at equinox: 1 - asin(Re/r)/pi = {100*(1-math.asin(Re/ra)/math.pi):.1f}%; suspension {R.ARRAY_OFFSET/1e3:.0f} km, self-weight stress rho g L = {1700*R.G_LOCAL*R.ARRAY_OFFSET/1e9:.2f} GPa (using g at the ring)")
print(f"sunlit share at the ring: {100*(1-math.asin(Re/r)/math.pi):.1f}%")
for W_ in (55.0+10, 91+10, 110.1+10):
    x = W_/(R.ARRAY_OFFSET*0.0093); print(f"shadow of {W_:.0f} m (panels + 10 m ring) on the array: band {R.ARRAY_OFFSET*0.0093:.0f} m, peak dimming {100*2/math.pi*(math.asin(x)+x*math.sqrt(1-x*x)):.1f}%")
d = math.degrees(math.atan(5950/R.ARRAY_OFFSET)); print(f"shadow on an 11.9 km array when |declination| < {d:.2f} deg: {365.25*2/math.pi*math.asin(math.sin(math.radians(d))/math.sin(math.radians(23.44))):.0f} days per year")
print(f"LIM sites {round(R.L_RING/500)}; x 8 MW = {round(R.L_RING/500)*8/1e3:.1f} GW; x 16 MW = {round(R.L_RING/500)*16/1e3:.1f} GW")
print(f"mass driver: LSM launch 222,015.5 km = {222015.5e3/R.L_RING:.2f} laps; LIM stages 370,712 km = {370712e3/R.L_RING:.2f} laps")
