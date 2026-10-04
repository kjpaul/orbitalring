import numpy as np, astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
import aep8
t=Time("2026-01-01"); lons=np.arange(-180,180,1.0)
def ring(alt,E,solar,p="p"):
    return np.array([aep8.flux(EarthLocation.from_geodetic(lo*u.deg,0*u.deg,alt*u.km),t,E*u.MeV,kind="integral",solar=solar,particle=p).to_value(u.cm**-2/u.s) for lo in lons])
def sect(f,thr):
    m=f>thr
    if not m.any(): return "none"
    if m.all(): return "all"
    z=np.where(~m)[0]; return f"{lons[(z.max()+1)%360]:.0f}..{lons[(z.min()-1)%360]:.0f} ({m.sum()} deg)"
YR=3.156e7
print("AP-8, geodetic latitude 0, 1 deg grid. p/cm2/s. Sector longitudes are in the AP-8 field epoch (1960s); today shift about 20 deg west.")
for solar in ("min","max"):
  print(f"\nPROTONS solar {solar}\n alt | >10MeV peak lon mean | >30 peak mean | >100 peak mean | sector>0 | sector>10 | sector>100 | panel 10% loss yr: worst, mean")
  for alt in range(500,1001,25):
    f=ring(alt,10,solar); g=ring(alt,30,solar); h=ring(alt,100,solar)
    w=3e10/(f.max()*YR) if f.max()>0 else float('inf'); m=3e10/(f.mean()*YR) if f.mean()>0 else float('inf')
    print(f"{alt:5d} | {f.max():8.1f} {lons[f.argmax()]:4.0f} {f.mean():8.2f} | {g.max():8.1f} {g.mean():8.2f} | {h.max():8.1f} {h.mean():8.2f} | {sect(f,0)} | {sect(f,10)} | {sect(f,100)} | {w:8.1f} {m:9.1f}")
print("\nELECTRONS AE-8 max\n alt | >0.5MeV peak lon mean sector>0 sector>100 | >1MeV peak mean")
for alt in range(500,1001,25):
    f=ring(alt,0.5,"max","e"); g=ring(alt,1.0,"max","e")
    print(f"{alt:5d} | {f.max():9.1f} {lons[f.argmax()]:4.0f} {f.mean():9.2f} {sect(f,0)} {sect(f,100)} | {g.max():9.1f} {g.mean():9.2f}")
