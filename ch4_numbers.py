import math, json
import numpy as np
import debris_orbits as do
import debris_simulation as ds
from debris_orbits import *
V=V_RECOIL
import ring_altitude as ring
H_DESIGN=ring.ALTITUDE
out={}
# Table 4.2 rows
rows=[]
for h_km in [250,500,600,800,1000,1500,2000,2500,3000]:
    h=h_km*1e3; n=released_cable(h); lo=released_cable(h,-V); hi=released_cable(h,V); w=released_cable(h,-V,V)
    rows.append(dict(h=h_km, v_cable=n['v_cable'], v_circ=n['v_circ'], excess=n['excess'], v_esc=n['v_esc'],
                     apo=n['h_a']/1e3, apo_lo=lo['h_a']/1e3, apo_hi=hi['h_a']/1e3, peri_min=w['h_p']/1e3,
                     slow500=slowdown_to_perigee(h,500e3), dvcrit_circ=ds.delta_v_critical(h)))
out['t42']=rows
for r in rows: print({k:(round(v) if isinstance(v,float) else v) for k,v in r.items()})
# lifetimes
s=ring_state(H_DESIGN); side=math.sqrt(s['m_struct']/RHO_CNT); Brod=ballistic_coefficient(side)
print("side",side,"Brod",Brod, "m/A", RHO_CNT*side)
life=[]
for hp_km in [150,200,250,300,400,500,600,700,800,1000,1500,2000]:
    hp=hp_km*1e3
    d=dict(hp=hp_km, circ_rod=lifetime(hp,hp,Brod), circ_600=lifetime(hp,hp,600.), circ_60=lifetime(hp,hp,60.))
    life.append(d); print(d)
out['life']=life
# released cable at ring altitude: worst-case (slow+radial) perigee and nominal orbit lifetimes for rod and B=600
rel=[]
for h_km in [250,400,500,600,700,800,1000,1500,2000]:
    h=h_km*1e3; n=released_cable(h); w=released_cable(h,-V,V)
    d=dict(h=h_km, peri_nom=n['h_p']/1e3, apo_nom=n['h_a']/1e3, peri_w=w['h_p']/1e3, apo_w=w['h_a']/1e3,
           life_nom_rod=lifetime(max(n['h_p'],1),n['h_a'],Brod), life_w_rod=lifetime(max(w['h_p'],1),w['h_a'],Brod),
           life_w_600=lifetime(max(w['h_p'],1),w['h_a'],600.), life_w_60=lifetime(max(w['h_p'],1),w['h_a'],60.))
    rel.append(d); print({k:(round(v,2) if isinstance(v,float) else v) for k,v in d.items()})
out['rel']=rel
# 250 km cases
c250=[]
for dv_t,dv_r,label in [(0,0,"nominal"),(-V,0,"slow half"),(V,0,"fast half"),(0,V,"radial recoil"),(-V,V,"slow half plus radial"),(V,V,"fast half plus radial")]:
    o=released_cable(250e3,dv_t,dv_r)
    d=dict(label=label, peri=o['h_p']/1e3, apo=o['h_a']/1e3, rod=lifetime(o['h_p'],o['h_a'],Brod), b600=lifetime(o['h_p'],o['h_a'],600.), b60=lifetime(o['h_p'],o['h_a'],60.))
    c250.append(d); print({k:(round(v,2) if isinstance(v,float) else v) for k,v in d.items()})
out['c250']=c250
# Pop 1
p1=[]
for h_km in [600,700,800,1000,1500,2000]:
    h=h_km*1e3
    d=dict(h=h_km, **{f"b{b}":ds.particulate_from_pop1(h,beta=b)/1e9 for b in [1.0,1.2,1.5,2.0,3.0]})
    p1.append(d); print({k:(round(v,2) if isinstance(v,float) else v) for k,v in d.items()})
out['p1']=p1
# 250 km prompt ejecta
s250=ring_state(250e3); L=s250['L_ring']
ej5=0.05*s250['m_cable_total']*L; print("250 km ejecta 5% cable only Mt", ej5/1e9, "soot 30%", 0.3*ej5/1e9)
ej_shock=0.003*s250['m_cable_total']*L; print("250 km ejecta 0.3% (shock layer) Mt", ej_shock/1e9, "soot", 0.3*ej_shock/1e9)
print("250 km: m_cable_total", s250['m_cable_total'], "L km", L/1e3, "cable total Mt", s250['m_cable_total']*L/1e9, "casing Mt", 12000*L/1e9, "ring Mt", s250['M_total']/1e9, "v_ground", s250['v_ground'])
s2=ring_state(2000e3); print("2000 km: cable total Mt", s2['m_cable_total']*s2['L_ring']/1e9, "casing Mt", 12000*s2['L_ring']/1e9, "v_ground", s2['v_ground'], "struct", s2['m_struct'])
sd=ring_state(H_DESIGN); print(f'{H_DESIGN/1e3:.0f} km: cable total Mt', sd['m_cable_total']*sd['L_ring']/1e9, 'casing Mt', 12000*sd['L_ring']/1e9, 'ring Mt', sd['M_total']/1e9, 'v_ground', sd['v_ground'], 'struct', sd['m_struct'], 'v_cable', sd['v_cable'], 'v_rel', sd['v_rel'])
for m in [3000]:
    for hh in sorted({H_DESIGN}):
        v=swept_cable_speed(hh,m,-V); hp,ha,_,_=orbit_from_state(R_EARTH+hh,v); print(f'slow half + roof 3000 at {hh/1e3:.0f} km:',v,hp/1e3,ha/1e3)
        print(f'casing needed slow half {hh/1e3:.0f}:', casing_mass_to_reach_perigee(hh,500e3,-V), 'nominal', casing_mass_to_reach_perigee(hh,500e3))
        r_=R_EARTH+hh; vg_=OMEGA_SIDEREAL*r_; print('casing tangential at ground',vg_*r_/R_EARTH,'east excess',vg_*r_/R_EARTH-OMEGA_SIDEREAL*R_EARTH)
# roof sweep at 2000: slow half
for m in [3000]:
    v=swept_cable_speed(2000e3,m,-V); hp,ha,_,_=orbit_from_state(R_EARTH+2000e3,v); print("slow half + roof 3000:",v,hp/1e3,ha/1e3)
print("casing needed slow half 2000:", casing_mass_to_reach_perigee(2000e3,500e3,-V))
# atmosphere densities
for h in [250,500,600,700,800,1000,2000]: print("rho",h,rho_atm(h*1e3), "H", scale_height(h*1e3)/1e3, "quick circ rod yr", lifetime_quick(h*1e3,Brod))
# recoil elastic energy
print("elastic energy per kg", (12.5e9)**2/(2*1700*500e9), "v", math.sqrt((12.5e9)**2/(1700*500e9)))
# casing fall: angular momentum eastward speed at ground
om=OMEGA_SIDEREAL; r=R_EARTH+2000e3; vg=om*r; L_=vg*r; v_tan_ground=L_/R_EARTH; print("casing tangential at ground",v_tan_ground,"ground speed",om*R_EARTH, "east excess", v_tan_ground-om*R_EARTH)
# casing KE at impact 5500 m/s
print("casing KE per m MJ", 0.5*12000*5500**2/1e6, "t TNT per m", 0.5*12000*5500**2/4.184e9, "total Gt TNT", 0.5*12000*5500**2*s2['L_ring']/4.184e18)
json.dump(out, open("ch4_numbers.json","w"), indent=1)
