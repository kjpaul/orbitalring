#!/usr/bin/env python3
"""AP-8 trapped proton flux around an equatorial orbital ring, by longitude and altitude.
Uses the NASA AP-8 model as packaged in `aep8` (IRBEM).  pip install aep8
The ring casing is ground-stationary, so each longitude is a fixed place on the ring."""
import numpy as np, astropy.units as u, sys, warnings
from astropy.coordinates import EarthLocation
from astropy.time import Time
import aep8
warnings.filterwarnings("ignore")
t = Time("2020-01-01")
lons = np.arange(-180, 180, 5.0)
def ring(alt_km, E, solar="min"):
    loc = EarthLocation.from_geodetic(lons*u.deg, 0*u.deg, alt_km*u.km)
    f = aep8.flux(loc, t, E*u.MeV, kind="integral", solar=solar, particle="p").to_value(1/(u.cm**2*u.s))
    return np.nan_to_num(f)
if __name__ == "__main__":
    print("AP-8 integral proton flux on the geographic equator, protons/cm^2/s, field epoch", t.iso[:10])
    for solar in ("min","max"):
        print(f"\n=== AP-8 {solar.upper()} ===")
        print(f"{'alt km':>7} {'E MeV':>6} {'ring mean':>11} {'peak':>11} {'peak lon':>9} {'min':>10} {'frac of ring >1% of peak':>26} {'lon range >10% of peak'}")
        for alt in (250,400,500,600,700,750,800,900,1000,1250,1500,2000,2500,3000):
            for E in (10,30,100):
                f = ring(alt,E,solar); pk=f.max(); 
                hi = lons[f>0.1*pk] if pk>0 else []
                rng = f"{hi.min():.0f} to {hi.max():.0f}" if len(hi) else "-"
                print(f"{alt:>7} {E:>6} {f.mean():>11.3g} {pk:>11.3g} {lons[f.argmax()]:>9.0f} {f.min():>10.3g} {np.mean(f>0.01*pk) if pk>0 else 0:>26.2f} {rng}")
    for alt in (750,800,2000):
        print(f"\nLongitude profile at {alt} km, >10 MeV, solar min:")
        f = ring(alt,10)
        print(" ".join(f"{lo:.0f}:{v:.3g}" for lo,v in zip(lons,f)))
