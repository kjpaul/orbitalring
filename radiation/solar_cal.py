"""Solar cell degradation calibration.
Datasheet (SolAero/Rocket Lab ZTJ triple junction): remaining power after 3 MeV protons 0.97 @4.5e9, 0.96 @2.3e10, 0.91 @8.8e10, 0.81 @3.4e11 p/cm2;
and 0.975 after 3 years in a 650 km, 30 deg orbit behind 4 mil coverglass.
Here: AP-8 MIN fluence above 10 MeV for that orbit, to express the on-orbit point as a >10 MeV trapped-proton fluence."""
import numpy as np, astropy.units as u, aep8, warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import EarthLocation; from astropy.time import Time
t = Time("2026-01-01"); rng = np.random.default_rng(1)
n = 40000; lon = rng.uniform(-180, 180, n); uu = rng.uniform(0, 2*np.pi, n)
lat = np.degrees(np.arcsin(np.sin(np.radians(30))*np.sin(uu)))       # uniform in argument of latitude; longitude uniform after many orbits
loc = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg, 650*u.km)
for solar in ("min", "max"):
    for E in (4, 10, 30):
        J = np.nan_to_num(aep8.flux(loc, t, E*u.MeV, kind="integral", solar=solar, particle="p").to_value(1/(u.cm**2*u.s)))
        print(f"650 km x 30 deg, AP-8 {solar.upper()}: orbit-average flux >{E} MeV {J.mean():.1f} p/cm2/s, 3-year fluence {J.mean()*3*3.156e7:.2e} p/cm2")
# fit of the 3 MeV datasheet points: P/P0 = 1 - C log10(1 + F/Fx)
from scipy.optimize import curve_fit
F = np.array([4.5e9, 2.3e10, 8.8e10, 3.4e11]); P = np.array([0.97, 0.96, 0.91, 0.81])
f = lambda F, C, Fx: 1 - C*np.log10(1 + F/Fx)
(C, Fx), _ = curve_fit(f, F, P, p0=(0.15, 2e10)); print(f"fit: C={C:.3f}, Fx={Fx:.2e}; residuals {np.round(f(F,C,Fx)-P,3)}")
from scipy.optimize import brentq
for loss in (0.025, 0.10, 0.20, 0.25): print(f"3 MeV proton fluence for {100*loss:.1f}% loss: {brentq(lambda x: f(x,C,Fx)-(1-loss), 1e8, 1e13):.2e} p/cm2")
