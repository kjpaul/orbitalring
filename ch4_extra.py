"""Chapter 4 numbers that the chapter text needs and ch4_numbers.py does not print:
collision energy and Table 4.1, gap crossing, contact numbers, force balance, the worked
radial-kick orbit of Section 4.4, roof-sweep cases, a 550/700 km row for Tables 4.2 and 4.5,
and the casing fall. Altitude from ring_altitude.py (ORBITAL_RING_ALT_KM)."""
import math
import ring_altitude as R
import debris_orbits as do, debris_simulation as ds
from debris_orbits import released_cable, swept_cable_speed, orbit_from_state, slowdown_to_perigee, casing_mass_to_reach_perigee, R_EARTH, V_RECOIL, GM
V=V_RECOIL; h=R.ALTITUDE; r=R.R_ORBIT; L=R.L_RING
mc=R.M_CABLE_M; vrel=R.V_REL
print(f"altitude {h/1e3:.0f} km  v_cable {R.V_CABLE:.1f}  v_ground {R.V_GROUND_SYNC:.1f}  v_rel {vrel:.1f}  v_orbit {R.V_ORBIT:.1f}  L {L/1e3:.0f} km")
print("--- Table 4.1")
for m in (12000,6000,3000,1200):
    mu=mc*m/(mc+m); E=0.5*mu*vrel**2
    print(f"  casing {m:6d}: reduced mass {mu:,.0f} kg/m  E {E/1e9:.1f} GJ/m  per kg {E/(mc+m)/1e6:.2f} MJ/kg  ring total {E*L:.3e} J = {E*L/4.184e18:.2f} Gt TNT = {E*L/4.184e18/0.05:.0f} x 50 Mt")
print(f"  fracture fronts: half circumference {L/2e3:.0f} km / 17,150 m/s = {L/2/17150:.0f} s = {L/2/17150/60:.1f} min")
d=R.design(h/1e3,False); print(f"  prograde: v_cable {d['v_cable']:.1f} v_rel {d['v_rel']:.1f} ({(1-d['v_rel']/vrel)*100:.0f}% lower); a_cent {d['v_cable']**2/r:.3f} net {d['v_cable']**2/r-R.G_LOCAL:.3f} rel {d['v_cable']**2/r:.3f}; side {d['side']:.2f}")
print("--- Section 4.2")
ac=R.V_CABLE**2/r; an=ac-R.G_LOCAL; t=math.sqrt(0.2/ac)
print(f"  a_cent {ac:.3f}  g {R.G_LOCAL:.3f}  a_net {an:.3f}  a_rel {ac:.3f}  t_gap {t:.4f} s  d_tan {vrel*t:.0f} m  v_radial {ac*t:.2f} m/s")
print(f"  Mach wedge half-angle {math.degrees(math.asin(4000/vrel)):.1f} deg; ratio 4000/v_rel {4000/vrel:.3f}; 5 cm fragment {0.05/vrel*1e6:.2f} us, depth {4000*0.05/vrel*100:.2f} cm, share of side {4000*0.05/vrel/R.CABLE_SIDE*100:.2f}%")
F=an*mc; print(f"  outward force {F/1e3:.1f} kN/m = casing weight {12000*R.G_NET/1e3:.1f} + T/r {R.CABLE_TENSION/r/1e3:.1f} = {(12000*R.G_NET+R.CABLE_TENSION/r)/1e3:.1f}; pressure on face {F/R.CABLE_SIDE/1e3:.1f} kPa")
print(f"  ejecta: 0.3% of cable {0.003*mc:.0f} kg/m; 5% {0.05*mc:.0f}; +1,200 = {0.05*mc+1200:.0f} ({(0.05*mc+1200)/R.M_RING_M*100:.1f}% of ring); +3,000 = {0.05*mc+3000:.0f} ({(0.05*mc+3000)/R.M_RING_M*100:.1f}%)")
print(f"  cable total {mc*L/1e9:.0f} Mt; casing total {12000*L/1e9:.0f} Mt")
print("--- Section 4.4")
two=2*GM/r; print(f"  2GM/r = {two:.4e}; excess {R.V_CABLE-R.V_ORBIT:.1f}")
for lab,dv in (("fast",V),("slow",-V)):
    v=R.V_CABLE+dv; den=two/v**2-1; print(f"  {lab} half v {v:.0f}: denominator {den:.4f}, r_apo {r/den/1e3:.0f} km, altitude {r/den/1e3-R_EARTH/1e3:.0f} km")
vt=R.V_CABLE-V; eps=(vt**2+V**2)/2-GM/r; a=GM/(2*abs(eps)); hh=r*vt; e=math.sqrt(1-hh*hh/(GM*a))
print(f"  worst case: v_t {vt:.1f}; KE {(vt**2+V**2)/2:.4e}; GM/r {GM/r:.4e}; eps {eps:.4e}; a {a/1e3:.0f} km; h {hh:.4e}; h^2/(GM a) {hh*hh/(GM*a):.5f}; e {e:.4f}; perigee r {a*(1-e)/1e3:.0f} = {a*(1-e)/1e3-R_EARTH/1e3:.0f} km alt; apogee r {a*(1+e)/1e3:.0f} = {a*(1+e)/1e3-R_EARTH/1e3:.0f} km alt")
print(f"  slowdown to 500 km perigee {slowdown_to_perigee(h,500e3):.0f} m/s; dv_crit from circular {ds.delta_v_critical(h):.1f} m/s")
mneed=casing_mass_to_reach_perigee(h,500e3,-V); print(f"  casing the slow half must sweep to reach 500 km: {mneed:.0f} kg/m = {mneed/120:.0f}% of casing")
for m in (1200,3000):
    v=swept_cable_speed(h,m,-V); hp,ha,_,_=orbit_from_state(R_EARTH+h,v); hp2,ha2,_,_=orbit_from_state(R_EARTH+h,v,V)
    print(f"  slow half + {m} kg/m roof: v {v:.0f} m/s, orbit {hp/1e3:.0f} x {ha/1e3:.0f} km; with radial kick {hp2/1e3:.0f} x {ha2/1e3:.0f} km")
print("--- Table 4.2 rows")
for hk in (250,500,550,600,650,700,750,800,1000,1500):
    H=hk*1e3; n=released_cable(H); lo=released_cable(H,-V); hi=released_cable(H,V); w=released_cable(H,-V,V)
    print(f"  {hk:5d} | {n['v_cable']:,.0f} | {n['v_circ']:,.0f} | {n['excess']:.0f} | {n['h_a']/1e3:,.0f} | {lo['h_a']/1e3:,.0f} to {hi['h_a']/1e3:,.0f} | {w['h_p']/1e3:.0f} | {slowdown_to_perigee(H,500e3):.0f} | v_esc {n['v_esc']:,.0f}")
print("--- Table 4.5 rows (5% ejecta, 30% soot, 10% casing), Mt, beta 1.0 1.2 1.5 2.0 3.0")
for hk in (250,500,550,600,650,700,750,800,1000,1500):
    print(f"  {hk:5d} | dv_crit {ds.delta_v_critical(hk*1e3):5.1f} | "+" | ".join(f"{ds.particulate_from_pop1(hk*1e3,beta=b,casing_frac=0.1)/1e9:.1f}" for b in (1.0,1.2,1.5,2.0,3.0)))
print("--- Casing fall (no drag, inverse-square gravity), from ring altitude to the surface")
def fall(hk):
    r0=R_EARTH+hk*1e3; v0=R.OMEGA_SIDEREAL*r0
    x,y,vx,vy=r0,0.0,0.0,v0; tt=0.0; dt=0.05
    while math.hypot(x,y)>R_EARTH:
        rr=math.hypot(x,y); ax=-GM*x/rr**3; ay=-GM*y/rr**3
        vx+=ax*dt; vy+=ay*dt; x+=vx*dt; y+=vy*dt; tt+=dt
    th=math.atan2(y,x); east=(th-R.OMEGA_SIDEREAL*tt)*R_EARTH; v=math.hypot(vx,vy)
    Lr=2*math.pi*r0
    print(f"  from {hk} km: time {tt:.0f} s ({tt/60:.1f} min), speed at surface {v:.0f} m/s, lands {east/1e3:.1f} km east; casing {12000*Lr/1e9:.0f} Mt; {0.5*12000*v*v/4.184e9:.1f} t TNT per m; total {0.5*12000*v*v*Lr/4.184e18:.2f} Gt TNT")
for hk in (250,700,800): fall(hk)
