"""Low-energy trapped protons (0.1 to 1 MeV) at the worst longitude, with the IGRF 2026 field. 2026-10-06.
The installed aep8 (0.1.0) only takes a location and uses AP-8's own 1960s field. Method: take the (L, B/B0) of the
2026 worst point from ap8_igrf2026_spectra.npz, find the equatorial point (longitude, altitude) with the same
(L, B/B0) in the IGRF field of 1960 (grids made with igrf_sector.L_B, 1 degree by 25 km, interpolated), and read
AP-8 there. Check: the flux above 10 MeV there must reproduce the stored 2026 value."""
import numpy as np, glob, os, astropy.units as u, aep8, warnings
from astropy.coordinates import EarthLocation
from astropy.time import Time
from scipy.interpolate import RegularGridInterpolator
warnings.filterwarnings("ignore")
YR = 365.25*86400; t = Time("1965-01-01")
fs = sorted(glob.glob(os.path.expanduser('~/fine_*.npz')), key=lambda f: int(os.path.basename(f)[5:-4]))
alts = np.array([int(os.path.basename(f)[5:-4]) for f in fs], float); g0 = np.load(fs[0]); lons = g0['lons']
Lg = np.array([np.load(f)['L'] for f in fs]); bg = np.array([np.load(f)['b'] for f in fs])
fL = RegularGridInterpolator((alts, lons), Lg); fb = RegularGridInterpolator((alts, lons), bg)
A, LO = np.meshgrid(np.arange(alts[0], alts[-1]+0.01, 1.0), np.arange(lons[0], lons[-1]+0.001, 0.05), indexing='ij')
P = np.stack([A.ravel(), LO.ravel()], -1); LL = fL(P); BB = fb(P)
d = np.load('ap8_igrf2026_spectra.npz'); Eg = list(d['Eg'])
for alt in (600, 700):
    for solar in ('min', 'max'):
        J10 = d[f'{solar}_{alt}'][Eg.index(10.0)]; i = int(np.argmax(d[f'min_{alt}'][Eg.index(10.0)]))
        Ls, bs = d[f'L_{alt}'][i], d[f'b_{alt}'][i]
        k = int(np.argmin(((LL-Ls)/0.01)**2 + ((BB-bs)/0.005)**2)); a, lo = P[k]
        loc = EarthLocation.from_geodetic(lo*u.deg, 0*u.deg, a*u.km)
        J = {E: float(np.nan_to_num(aep8.flux(loc, t, E*u.MeV, kind="integral", solar=solar, particle="p").value)) for E in (0.1, 1.0, 10.0)}
        band = J[0.1]-J[1.0]
        print(f"{alt} km solar {solar}: 2026 worst lon {d['lons'][i]:.0f}, L {Ls:.4f}, B/B0 {bs:.4f}, stored J>10 {J10[i]:.1f}; "
              f"matched 1960 point lon {lo:.2f}, alt {a:.0f} km (L {LL[k]:.4f}, B/B0 {BB[k]:.4f}); AP-8 there: J>10 {J[10.0]:.1f} (check, ratio {J[10.0]/J10[i]:.2f}), "
              f"J>1 {J[1.0]:.1f}, J>0.1 {J[0.1]:.1f} /cm2/s; 0.1-1 MeV {band:.1f} /cm2/s = {band*YR*30:.3g} /cm2 in 30 yr")
