"""AP-8 integral proton spectra on the geographic equator, 1 deg longitude steps -> ap8_spectra.npz"""
import numpy as np, astropy.units as u, warnings, aep8
warnings.filterwarnings("ignore")
from astropy.coordinates import EarthLocation
from astropy.time import Time
t = Time("2026-01-01")   # aep8 uses the model's own field epoch; the date does not change the result (checked below)
Eg = np.array([10,15,20,30,40,50,60,80,100,125,150,175,200,250,300,350,400.])
lons = np.arange(-180, 180, 1.0)
alts = [550,600,650,700,750,800,850,900,925,950,1000]
out = {}
for solar in ("min", "max"):
    mdl = aep8.model("p", solar)
    for alt in alts:
        loc = EarthLocation.from_geodetic(lons*u.deg, 0*u.deg, alt*u.km)
        L, BB0 = mdl.geomagnetic_coordinates(loc, t)
        I = np.array([mdl.integral_flux_for_geomagnetic_coordinates(L, BB0, E*u.MeV).to_value(u.cm**-2/u.s) for E in Eg])
        out[f"{solar}_{alt}"] = I; out[f"L_{solar}_{alt}"] = np.asarray(L); out[f"BB0_{solar}_{alt}"] = np.asarray(BB0)
np.savez("ap8_spectra.npz", Eg=Eg, lons=lons, alts=alts, **out)
for alt in (600,700,750,800):
    I = out[f"min_{alt}"]; k = I[0].argmax()
    nz = lons[I[0] > 0]
    print(f"{alt} km solar min: peak >10 MeV {I[0,k]:.1f} at lon {lons[k]:.0f}; nonzero {nz.min() if len(nz) else None}..{nz.max() if len(nz) else None} ({len(nz)} deg); L={out[f'L_min_{alt}'][k]:.3f} B/B0={out[f'BB0_min_{alt}'][k]:.3f}")
    print("   spectrum at peak:", dict(zip(Eg.astype(int), np.round(I[:,k],1))))
loc = EarthLocation.from_geodetic(-27*u.deg, 0*u.deg, 800*u.km)
print("date check 1965 vs 2026:", [aep8.flux(loc, Time(d), 10*u.MeV, kind="integral", solar="min", particle="p") for d in ("1965-01-01","2026-01-01")])
