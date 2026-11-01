"""
Book: Orbital Ring Engineering, Volume I: Chapter 13, Space Elevators (the pod descending the space elevator through the Van Allen belts).
Run from the radiation folder: python elevator_dose.py > elevator_dose.txt   (needs aep8, astropy, numpy; uses fold_trapped.py)

Trapped-proton dose to a pod riding the tether from GEO down to the 700 km ring, on the equator at a chosen longitude.
AP-8 integral spectra (solar min and max) at each altitude, folded with the Geant4 water depth-dose responses of
fold_trapped.py, weighted by the time the pod spends at that altitude on the descent profiles of Chapter 13 (Table 13.6)."""
import numpy as np, warnings, sys, os
warnings.filterwarnings("ignore")
import astropy.units as u, aep8
from astropy.coordinates import EarthLocation
from astropy.time import Time
import fold_trapped as ft

GM=3.986004418e14; RE=6378137.0; W=7.2921159e-5; g0=9.80665
r_geo=(GM/W**2)**(1/3); r_r=RE+700e3
def anet(r): return GM/r**2 - W*W*r
def profile(v_c, a_drive, a_stop, dt=1.0):
    """altitude (m) at each second of the descent"""
    alts=[]; r=r_geo; v=0.0
    while v<v_c:
        v+=(a_drive+anet(r))*dt; r-=v*dt; alts.append(r-RE)
    r3=r_r; v3=0.0
    while v3<v_c: v3+=a_stop*dt; r3+=v3*dt
    r_stop=r3
    while r>r_stop:
        r-=v_c*dt; alts.append(r-RE)
    r3=r_r; v3=0.0; stop=[]
    while v3<v_c: v3+=a_stop*dt; r3+=v3*dt; stop.append(r3-RE)
    alts.extend(reversed(stop))
    return np.array(alts)

cases={"2,000 m/s, 1 g drive, 3 m/s2 stop": (2000.0, g0, 3.0),
       "1,000 m/s, 1 g drive, 1 m/s2 stop": (1000.0, g0, 1.0)}
edges=np.arange(700e3, 35786e3+1, 100e3); edges[-1]=35786e3+1
mid=0.5*(edges[1:]+edges[:-1])
t=Time("2026-01-01"); Eg=ft.np.array([10,15,20,30,40,50,60,80,100,125,150,175,200,250,300,350,400.])
T=ft.T; YR=ft.YR
print("shield depths (g/cm2 of water, centre of a sphere):", list(T))
print("in aluminum at 2.7 g/cm3 the same g/cm2 is 1 g/cm2 = 3.7 mm;  10 g/cm2 = 3.7 cm;  20 g/cm2 = 7.4 cm")
spectra={}
for solar in ("min","max"):
    for lon in (-30.0, 0.0, 60.0, 120.0, 180.0):
        loc=EarthLocation.from_geodetic(lon*u.deg, 0*u.deg, (mid/1e3)*u.km)
        I=np.array([aep8.flux(loc, t, E*u.MeV, kind="integral", solar=solar, particle="p").to_value(u.cm**-2/u.s) for E in Eg])  # (nE, nalt)
        spectra[(solar,lon)]=I
        k=np.nanargmax(np.nan_to_num(I[0]))
        print(f"solar {solar}, lon {lon:+.0f}: >10 MeV flux peaks at {mid[k]/1e3:,.0f} km with {I[0,k]:,.0f} /cm2/s ; >100 MeV there {I[8,k]:,.0f} /cm2/s")
# dose rate per altitude (Gy/yr, Sv/yr) at each shield depth, for the worst case
rates={}
for key,I in spectra.items():
    D=np.zeros((len(T),len(mid))); H=np.zeros((len(T),len(mid)))
    for j in range(len(mid)):
        d,h=ft.dose(Eg, I[:,j], tail="exp"); D[:,j]=d; H[:,j]=h
    rates[key]=(D,H)
print("\n-- dose rate along the tether, solar max, lon -30 (worst), Sv per hour --")
D,H=rates[("max",-30.0)]
for alt_km in (1000,1500,2000,2500,3000,3500,4000,5000,6000,8000,10000,15000,20000):
    j=np.searchsorted(mid, alt_km*1e3)
    print(f"  {alt_km:>6,} km: " + " ; ".join(f"{T[i]:g} g/cm2 {H[i,j]/YR*3600*1e3:.3g} mSv/h" for i in (0,2,4,6,8,9,11,13,15)))
print("\n-- dose per descent (mSv, dose equivalent; Gy absorbed in brackets), by case, solar activity, longitude --")
for name,(v_c,ad,as_) in cases.items():
    alts=profile(v_c,ad,as_); hist,_=np.histogram(alts, bins=edges); sec=hist.astype(float)
    print(f"\n{name}: {sec.sum()/3600:.2f} h on the tether ; {sec[mid<6000e3].sum()/60:.0f} min below 6,000 km")
    for solar in ("min","max"):
        for lon in (-30.0, 0.0, 60.0, 120.0, 180.0):
            D,H=rates[(solar,lon)]
            perT_H=(H*sec[None,:]).sum(1)/YR*1e3; perT_D=(D*sec[None,:]).sum(1)/YR*1e3
            print(f"  solar {solar}, lon {lon:+4.0f}: " + " ; ".join(f"{T[i]:g} g/cm2 {perT_H[i]:.3g} mSv ({perT_D[i]:.3g} mGy)" for i in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,15)))

# ---------------- trapped electrons (AE-8), outer belt mostly, with bremsstrahlung ----------------
# Geant4 electron depth-dose runs in water: out_elec/ (QBBC, 60 cm slab, 0.2 cm bins, 100,000 events per energy;
# DEPTH_CM=60 NBIN=300 python g4_depthdose.py QBBC e- out_elec 100000:0.3 ... 100000:7.0). All deposition is electrons and photons, Q = 1.
import glob, re
eruns=[]
for f in sorted(glob.glob("out_elec/*.npz")):
    d=np.load(f); names=list(d["names"]); dat=dict(zip(names, d["data"]/float(d["nev"])))
    E=float(re.search(r"_([\d.]+)MeV", f).group(1)); w=0.2
    eruns.append(dict(E=E, z=(np.arange(300)+0.5)*w, R=dat["tot"]/w))   # MeV/g per (e/cm2) at depth z (g/cm2)
eruns.sort(key=lambda r: r["E"]); EE=np.array([r["E"] for r in eruns])
def eresp(E, t):
    if E<=EE[0]: return np.interp(t, eruns[0]["z"], eruns[0]["R"], right=0.0)*(E/EE[0])
    i=np.clip(np.searchsorted(EE,E)-1, 0, len(eruns)-2); a,b=eruns[i],eruns[i+1]
    f=(np.log(E)-np.log(a["E"]))/(np.log(b["E"])-np.log(a["E"]))
    return (1-f)*np.interp(t,a["z"],a["R"],right=0.0)+f*np.interp(t,b["z"],b["R"],right=0.0)
Ege=np.array([0.04,0.1,0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0])
EFe=np.geomspace(0.3, 7.0, 200); EFm=0.5*(EFe[1:]+EFe[:-1])
RESPe=np.array([[eresp(E,t) for E in EFm] for t in T])
def edose_rate(I):
    """Gy/s at the centre of water spheres of radius T from an AE-8 integral spectrum I(Ege)"""
    I=np.nan_to_num(np.asarray(I,float)); ok=I>0
    if ok.sum()<2: return np.zeros(len(T))
    li=np.interp(np.log(EFe), np.log(Ege[ok]), np.log(I[ok]), right=-np.inf); If=np.exp(li); n=-np.diff(If); n[n<0]=0
    return (RESPe*n[None,:]).sum(1)*ft.GY
print("\n-- trapped electrons (AE-8): flux along the tether and dose per descent (Gy = Sv for electrons and photons) --")
for solar in ("min","max"):
    for lon in (-30.0, 0.0, 120.0):
        loc=EarthLocation.from_geodetic(lon*u.deg, 0*u.deg, (mid/1e3)*u.km)
        Ie=np.array([aep8.flux(loc, t, E*u.MeV, kind="integral", solar=solar, particle="e").to_value(u.cm**-2/u.s) for E in Ege])
        k=np.nanargmax(np.nan_to_num(Ie[3])); 
        print(f" solar {solar}, lon {lon:+.0f}: >0.5 MeV electrons peak at {mid[k]/1e3:,.0f} km with {Ie[3,k]:.3g} /cm2/s ; >2 MeV there {Ie[7,k]:.3g} ; >0.5 MeV at 20,000 km {Ie[3,np.searchsorted(mid,20000e3)]:.3g}, at 30,000 km {Ie[3,np.searchsorted(mid,30000e3)]:.3g}")
        De=np.array([edose_rate(Ie[:,j]) for j in range(len(mid))]).T   # (nT, nalt) Gy/s
        for name,(v_c,ad,as_) in cases.items():
            alts=profile(v_c,ad,as_); hist,_=np.histogram(alts, bins=edges); sec=hist.astype(float)
            per=(De*sec[None,:]).sum(1)*1e3
            print(f"   {name}: " + " ; ".join(f"{T[i]:g} g/cm2 {per[i]:.3g} mSv" for i in (0,1,2,3,4,6,8,9,11)))

# ---------------- combined table for Chapter 13, Space Elevators (Table 13.6): worst longitude (0 deg), solar max ----------------
print("\n-- Table for Chapter 13 (Table 13.6): dose per descent behind water shields, solar max, longitude 0 (the worst of the five), protons + electrons --")
loc=EarthLocation.from_geodetic(0*u.deg, 0*u.deg, (mid/1e3)*u.km)
Ie=np.array([aep8.flux(loc, t, E*u.MeV, kind="integral", solar="max", particle="e").to_value(u.cm**-2/u.s) for E in Ege])
De=np.array([edose_rate(Ie[:,j]) for j in range(len(mid))]).T
Dp,Hp=rates[("max",0.0)]
print(" shield g/cm2 | Al mm | PE cm |  2,000 m/s: protons mSv, electrons mSv, total mSv |  1,000 m/s: total mSv")
res={}
for name,(v_c,ad,as_) in cases.items():
    alts=profile(v_c,ad,as_); hist,_=np.histogram(alts, bins=edges); sec=hist.astype(float)
    res[name]=((Hp*sec[None,:]).sum(1)/YR*1e3, (De*sec[None,:]).sum(1)*1e3)
f,s=list(cases)
for i in range(len(T)):
    if T[i]>150: break
    print(f"  {T[i]:6g}   | {T[i]/2.7*10:5.0f} | {T[i]/0.95:5.1f} | {res[f][0][i]:8.3g} {res[f][1][i]:8.3g} {res[f][0][i]+res[f][1][i]:8.3g} | {res[s][0][i]+res[s][1][i]:8.3g}")
tot=res[f][0]+res[f][1]
for target in (20.0, 10.0, 5.0, 2.0, 1.0):
    i=np.searchsorted(-tot, -target)   # tot is decreasing
    if i==0 or i>=len(T): print(f"  {target} mSv: outside the table"); continue
    # log interpolation in T
    x=np.interp(np.log(target), np.log(tot[[i,i-1]]), np.log(T[[i,i-1]])); tw=np.exp(x)
    print(f"  shield for {target:g} mSv per descent at 2,000 m/s: {tw:.0f} g/cm2 = {tw:.0f} cm of water, {tw/2.7:.1f} cm of aluminum, {tw/0.95:.0f} cm of polyethylene ; {tw*10:.0f} kg/m2")
# shield mass for a compartment and for a whole pod
for label,area in (("compartment 2.5 m across, 3 m long", np.pi*2.5*3+2*np.pi*1.25**2), ("pod cabin 3 m across, 6 m long", np.pi*3*6+2*np.pi*1.5**2)):
    print(f"  {label}: {area:.0f} m2 ; at 35 g/cm2 {area*350/1e3:.1f} t ; at 22 g/cm2 {area*220/1e3:.1f} t")
