# Trapped-proton dose at the centre of a spherical water-equivalent shield, CSDA slowing only.
# Range-energy (Bragg-Kleeman fit to NIST PSTAR, water): R = A*E**P g/cm2, A=0.0022, P=1.77
# No nuclear attenuation, no secondaries (neutrons). Absorbed dose in water (Gy); not Sv.
import numpy as np, astropy.units as u, warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import EarthLocation
from astropy.time import Time
import aep8
A,P=0.0022,1.77
t=Time("2026-01-01"); lons=np.arange(-180,180,2.0)
Eg=np.array([10,15,20,30,40,50,60,80,100,125,150,200,250,300.])
def integ(alt,lo,solar="min"):
    loc=EarthLocation.from_geodetic(lo*u.deg,0*u.deg,alt*u.km)
    return np.array([aep8.flux(loc,t,E*u.MeV,kind="integral",solar=solar,particle="p").to_value(u.cm**-2/u.s) for E in Eg])
def dose_gy_yr(I,tk):
    # integral spectrum interpolated log-linearly onto a fine grid; above 300 MeV the 250-300 MeV
    # exponential slope is continued to 400 MeV (AP-8 upper limit). Flux above 400 MeV set to zero.
    if I[0]<=0: return 0.0
    ok=I>0; Ef=np.linspace(10,400,3901)
    li=np.interp(Ef,Eg[ok],np.log(I[ok]))
    if ok[-1] and ok[-2]:
        k=(np.log(I[-1])-np.log(I[-2]))/(Eg[-1]-Eg[-2]); m=Ef>Eg[-1]; li[m]=np.log(I[-1])+k*(Ef[m]-Eg[-1])
    If=np.exp(li); n=-np.diff(If); n=np.append(n,0.0)[:len(Ef)-1]; Em=0.5*(Ef[1:]+Ef[:-1])
    r=A*Em**P-tk; m=r>0
    Er=np.maximum((r[m]/A)**(1/P),1.0)            # residual energy floored at 1 MeV
    S=Er**(1-P)/(A*P)
    return float((n[m]*S).sum()*1.602e-13*1e3*3.156e7)
T=[1,5,10,20,30,50,75,100]
for alt in (700,750,800):
    I=[integ(alt,lo) for lo in lons]
    pk=max(I,key=lambda x:x[0])
    print(f"\n{alt} km, solar min. Worst-longitude integral flux p/cm2/s at E>",dict(zip(Eg.astype(int),np.round(pk,1))))
    print(" shield g/cm2 | worst longitude Gy/yr | ring mean Gy/yr | share of ring > 0.02 Gy/yr | > 0.001 Gy/yr")
    for tk in T:
        D=np.array([dose_gy_yr(x,tk) for x in I])
        print(f"   {tk:5d}      | {D.max():10.3f}            | {D.mean():8.4f}        | {100*(D>0.02).mean():5.1f}%   | {100*(D>0.001).mean():5.1f}%")
