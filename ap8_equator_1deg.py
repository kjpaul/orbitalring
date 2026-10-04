import numpy as np, astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
import aep8
t=Time("2026-01-01")
lons=np.arange(-180,180,1.0)
def ring(alt,E,solar,lat=0.0):
    out=[]
    for lo in lons:
        loc=EarthLocation.from_geodetic(lo*u.deg,lat*u.deg,alt*u.km)
        out.append(aep8.flux(loc,t,E*u.MeV,kind="integral",solar=solar,particle="p").to_value(u.cm**-2/u.s))
    return np.array(out)
print("geodetic latitude 0.000, 1 deg longitude grid, integral proton flux p/cm2/s")
for E in (10,30,100):
  for solar in ("min","max"):
    for alt in (500,550,600,650,700,750,800,1000,1250,1500,2000):
        f=ring(alt,E,solar); nz=lons[f>0]
        print(f">{E}MeV {solar} {alt:5d} km peak {f.max():10.1f} at lon {lons[f.argmax()]:6.0f} mean {f.mean():10.2f} nonzero {100*(f>0).mean():5.1f}% span {nz.min() if len(nz) else 0:.0f}..{nz.max() if len(nz) else 0:.0f}")
print("lat sensitivity at 800 km, >10 MeV min: lat, peak, mean")
for lat in (-2,-1,0,1,2):
    f=ring(800,10,"min",lat); print(lat,round(f.max(),1),round(f.mean(),1))
