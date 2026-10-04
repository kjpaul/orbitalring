"""Measured proton flux on the geographic equator at about 820 km, MetOp-C MEPED omnidirectional detectors
(NOAA NCEI processed files, August-September 2026), compared with AP-8 at the same altitude.
The MEPED 'omni' flux is differential (per cm2 s sr MeV) at 25, 50 and 100 MeV, from a dome detector that looks at the zenith."""
import numpy as np, netCDF4 as nc, glob, astropy.units as u, aep8, warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import EarthLocation; from astropy.time import Time
lat, lon, alt, f1, f2, f3 = [], [], [], [], [], []
files = sorted(glob.glob("poes/poes_m03_2026*_proc.nc"))
for f in files:
    try: d = nc.Dataset(f)
    except Exception as e: print("skip", f, e); continue
    la = d["lat"][:].filled(np.nan); m = np.abs(la) < 1.0
    lat.append(la[m]); lo = d["lon"][:].filled(np.nan)[m]; lon.append(((lo+180) % 360)-180); alt.append(d["alt"][:].filled(np.nan)[m])
    for arr, k in ((f1,"mep_omni_flux_p1"),(f2,"mep_omni_flux_p2"),(f3,"mep_omni_flux_p3")): arr.append(d[k][:].filled(np.nan)[m])
lat, lon, alt, f1, f2, f3 = map(np.concatenate, (lat, lon, alt, f1, f2, f3))
print(f"{len(files)} daily files, {len(lat)} samples within 1 deg of the equator, altitude {np.nanmin(alt):.0f} to {np.nanmax(alt):.0f} km (mean {np.nanmean(alt):.0f})")
for name, f in (("25 MeV", f1), ("50 MeV", f2), ("100 MeV", f3)):
    ok = np.isfinite(f) & (f >= 0)
    print(name, "valid", ok.sum(), "median outside sector (lon 60..180):", np.nanmedian(f[ok & (lon > 60)]), "max", np.nanmax(f[ok]))
bins = np.arange(-180, 181, 5); c = 0.5*(bins[1:]+bins[:-1])
mdl = aep8.model("p", "min"); A = float(np.nanmean(alt))
loc = EarthLocation.from_geodetic(c*u.deg, 0*u.deg, A*u.km); t = Time("2026-01-01")
ap = {E: np.nan_to_num(mdl.differential_flux(loc, t, E*u.MeV).to_value(1/(u.cm**2*u.s*u.MeV)))/(4*np.pi) for E in (25, 50, 100)}
print(f"\nlon bin | MEPED median flux at 25 / 50 / 100 MeV (per cm2 s sr MeV) | n | AP-8 MIN (1960s field) omnidirectional/4pi at {A:.0f} km, same energies")
rows = []
for i in range(len(c)):
    m = (lon >= bins[i]) & (lon < bins[i+1])
    v = [np.nanmedian(f[m & np.isfinite(f)]) if (m & np.isfinite(f)).sum() else np.nan for f in (f1, f2, f3)]
    rows.append([c[i], *v, m.sum(), ap[25][i], ap[50][i], ap[100][i]])
    print(f"{bins[i]:5d}..{bins[i+1]:4d} | {v[0]:9.3f} {v[1]:9.3f} {v[2]:9.3f} | {m.sum():4d} | {ap[25][i]:8.3f} {ap[50][i]:8.3f} {ap[100][i]:8.3f}")
np.save("poes_equator_rows.npy", np.array(rows))
# sector limits at fixed fractions of the peak, 2-degree bins, 25 MeV and 100 MeV channels
b2 = np.arange(-180, 181, 2); c2 = 0.5*(b2[1:]+b2[:-1])
for name, f in (("25 MeV", f1), ("100 MeV", f3)):
    med = np.array([np.nanmedian(f[(lon >= b2[i]) & (lon < b2[i+1])]) for i in range(len(c2))])
    bg = np.nanmedian(med[c2 > 60]); pk = np.nanmax(med)
    print(f"\n{name}: peak {pk:.3f} at lon {c2[np.nanargmax(med)]:.0f}, background (lon > 60) {bg:.4f}")
    for frac in (0.5, 0.1, 0.01):
        m = med - bg > frac*(pk - bg); print(f"   above {100*frac:.0f}% of peak: {c2[m].min()-1:.0f} to {c2[m].max()+1:.0f} ({m.sum()*2} deg)")
