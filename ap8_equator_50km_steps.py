import numpy as np, astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
import aep8
t=Time("2026-01-01"); lons=np.arange(-180,180,1.0)
def ring(alt,E,solar,p="p"):
    return np.array([aep8.flux(EarthLocation.from_geodetic(lo*u.deg,0*u.deg,alt*u.km),t,E*u.MeV,kind="integral",solar=solar,particle=p).to_value(u.cm**-2/u.s) for lo in lons])
def span(f):
    nz=lons[f>0]
    if len(nz)==0: return "none"
    if len(nz)==360: return "all"
    # find wrap-aware contiguous start/end
    z=np.where(f<=0)[0]; start=lons[(z.max()+1)%360]; end=lons[(z.min()-1)%360]
    return f"{start:.0f}..{end:.0f}"
print("AP-8 protons, geodetic latitude 0, 1 deg longitude grid. Flux in p/cm2/s. 'clean' = share of ring with zero flux.")
for solar in ("min","max"):
  print(f"\nsolar {solar}")
  print(" alt_km | >10MeV: peak  lon  mean  min  clean%  sector | >30MeV peak mean | >100MeV peak mean | share>100 share>1000 (>10MeV)")
  for alt in range(600,1551,50):
    f=ring(alt,10,solar); g=ring(alt,30,solar); h=ring(alt,100,solar)
    print(f"{alt:6d} | {f.max():9.1f} {lons[f.argmax()]:5.0f} {f.mean():9.1f} {f.min():8.1f} {100*(f<=0).mean():5.1f} {span(f):>10} | {g.max():9.1f} {g.mean():9.1f} | {h.max():9.1f} {h.mean():9.1f} | {100*(f>100).mean():5.1f} {100*(f>1000).mean():5.1f}")
print("\nAE-8 electrons >0.5 MeV and >1 MeV, solar max (e/cm2/s): alt, peak, mean, clean% (>0.5) | peak, mean (>1)")
for alt in range(600,1551,50):
    f=ring(alt,0.5,"max","e"); g=ring(alt,1.0,"max","e")
    print(f"{alt:6d} {f.max():11.1f} {f.mean():11.1f} {100*(f<=0).mean():5.1f} | {g.max():11.1f} {g.mean():11.1f}")
