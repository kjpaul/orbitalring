"""
Book: Orbital Ring Engineering, Volume I: Chapter 13 (the pod descending the space elevator through the Van Allen belts).
Run from the radiation folder: python elevator_dose_rates.py > elevator_dose_rates.txt   (needs aep8, astropy, numpy; uses fold_trapped.py and out_elec/)
Environment: ORBITAL_RING_OUT (output folder for the rate file, default ../runs)

Dose rate along the tether, from the ring at 700 km to the platform at GEO, at the centre of water spheres of radius T (g/cm2),
for trapped protons (AP-8, folded with the Geant4 depth-dose responses of fold_trapped.py) and trapped electrons (AE-8, folded
with the Geant4 electron runs in out_elec/). The rates are written to elevator_dose_rates.npz so that elevator_system.py can
weight them by the time any descent profile spends at each altitude without needing the aep8 package. The method is that of
elevator_dose.py (Table 13.6); this script only separates the rates from the profile.

Cases: solar max at longitudes 0 (the worst of the five longitudes checked for Table 13.6), -30 and +120 deg, and solar min at 0 deg.
The check at the end refolds the Table 13.6 profile (2,000 m/s, 1 g drive, 3 m/s2 stop) and prints it beside the table's values."""
import numpy as np, warnings, sys, os, glob, re
warnings.filterwarnings("ignore")
import astropy.units as u, aep8
from astropy.coordinates import EarthLocation
from astropy.time import Time
import fold_trapped as ft

OUT = os.environ.get("ORBITAL_RING_OUT", os.path.join("..", "runs"))
os.makedirs(OUT, exist_ok=True)
GM = 3.986004418e14; RE = 6378137.0; W = 7.2921159e-5; g0 = 9.80665
r_geo = (GM / W**2)**(1 / 3); r_r = RE + 700e3
T = ft.T; YR = ft.YR
t = Time("2026-01-01")
edges = np.arange(700e3, 35786e3 + 1, 100e3); edges[-1] = 35786e3 + 1
mid = 0.5 * (edges[1:] + edges[:-1])
Eg = np.array([10, 15, 20, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350, 400.])

# ---- electron response from the Geant4 runs in out_elec/ (as in elevator_dose.py) ----
eruns = []
for f in sorted(glob.glob("out_elec/*.npz")):
    d = np.load(f); names = list(d["names"]); dat = dict(zip(names, d["data"] / float(d["nev"])))
    E = float(re.search(r"_([\d.]+)MeV", f).group(1)); w = 0.2
    eruns.append(dict(E=E, z=(np.arange(300) + 0.5) * w, R=dat["tot"] / w))   # MeV/g per (e/cm2) at depth z (g/cm2)
eruns.sort(key=lambda r: r["E"]); EE = np.array([r["E"] for r in eruns])
def eresp(E, tt):
    if E <= EE[0]: return np.interp(tt, eruns[0]["z"], eruns[0]["R"], right=0.0) * (E / EE[0])
    i = np.clip(np.searchsorted(EE, E) - 1, 0, len(eruns) - 2); a, b = eruns[i], eruns[i + 1]
    f = (np.log(E) - np.log(a["E"])) / (np.log(b["E"]) - np.log(a["E"]))
    return (1 - f) * np.interp(tt, a["z"], a["R"], right=0.0) + f * np.interp(tt, b["z"], b["R"], right=0.0)
Ege = np.array([0.04, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])
EFe = np.geomspace(0.3, 7.0, 200); EFm = 0.5 * (EFe[1:] + EFe[:-1])
RESPe = np.array([[eresp(E, tt) for E in EFm] for tt in T])
def edose_rate(I):
    """Gy/s at the centre of water spheres of radius T from an AE-8 integral spectrum I(Ege)"""
    I = np.nan_to_num(np.asarray(I, float)); ok = I > 0
    if ok.sum() < 2: return np.zeros(len(T))
    li = np.interp(np.log(EFe), np.log(Ege[ok]), np.log(I[ok]), right=-np.inf); If = np.exp(li); n = -np.diff(If); n[n < 0] = 0
    return (RESPe * n[None, :]).sum(1) * ft.GY

print("elevator_dose_rates.py: dose rates along the tether, 100 km bins from 700 km to GEO")
print("shield depths T (g/cm2 of water, centre of a sphere):", list(T))
print("electron runs:", len(eruns), "energies", EE[0], "to", EE[-1], "MeV")
out = dict(mid=mid, edges=edges, T=T)
cases = [("max", 0.0), ("max", -30.0), ("max", 120.0), ("min", 0.0)]
for solar, lon in cases:
    loc = EarthLocation.from_geodetic(lon * u.deg, 0 * u.deg, (mid / 1e3) * u.km)
    Ip = np.array([aep8.flux(loc, t, E * u.MeV, kind="integral", solar=solar, particle="p").to_value(u.cm**-2 / u.s) for E in Eg])
    Ie = np.array([aep8.flux(loc, t, E * u.MeV, kind="integral", solar=solar, particle="e").to_value(u.cm**-2 / u.s) for E in Ege])
    Dp = np.zeros((len(T), len(mid))); Hp = np.zeros((len(T), len(mid))); De = np.zeros((len(T), len(mid)))
    for j in range(len(mid)):
        d, h = ft.dose(Eg, Ip[:, j], tail="exp"); Dp[:, j] = d / YR; Hp[:, j] = h / YR    # Gy/s, Sv/s
        De[:, j] = edose_rate(Ie[:, j])                                                    # Gy/s (= Sv/s for electrons)
    key = f"{solar}_{int(lon):+d}"
    out[f"Hp_{key}"] = Hp; out[f"Dp_{key}"] = Dp; out[f"De_{key}"] = De
    kp = np.nanargmax(np.nan_to_num(Ip[0])); ke = np.nanargmax(np.nan_to_num(Ie[3]))
    print(f"\nsolar {solar}, lon {lon:+.0f}: >10 MeV protons peak at {mid[kp]/1e3:,.0f} km with {Ip[0,kp]:,.0f} /cm2/s ; "
          f">0.5 MeV electrons peak at {mid[ke]/1e3:,.0f} km with {Ie[3,ke]:.3g} /cm2/s")
    print("  altitude   | proton dose equivalent rate (mSv/h) at 5, 20, 50 g/cm2 | electron dose rate (mSv/h) at 2, 5, 20 g/cm2")
    for alt_km in (1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 25000, 30000, 35000):
        j = np.searchsorted(mid, alt_km * 1e3)
        print(f"  {alt_km:>7,} km | " + " ; ".join(f"{Hp[i,j]*3600e3:.3g}" for i in (4, 8, 11)) + " | " + " ; ".join(f"{De[i,j]*3600e3:.3g}" for i in (2, 4, 8)))
np.savez(os.path.join(OUT, "elevator_dose_rates.npz"), **out)
print(f"\nwrote {os.path.join(OUT, 'elevator_dose_rates.npz')}")

# ---- check: refold the Table 13.6 profile (2,000 m/s, 1 g drive, 3 m/s2 stop; solar max, lon 0) ----
def anet(r): return GM / r**2 - W * W * r
def profile(v_c, a_drive, a_stop, dt=1.0):
    alts = []; r = r_geo; v = 0.0
    while v < v_c:
        v += (a_drive + anet(r)) * dt; r -= v * dt; alts.append(r - RE)
    r3 = r_r; v3 = 0.0
    while v3 < v_c: v3 += a_stop * dt; r3 += v3 * dt
    r_stop = r3
    while r > r_stop:
        r -= v_c * dt; alts.append(r - RE)
    r3 = r_r; v3 = 0.0; stop = []
    while v3 < v_c: v3 += a_stop * dt; r3 += v3 * dt; stop.append(r3 - RE)
    alts.extend(reversed(stop))
    return np.array(alts)
alts = profile(2000.0, g0, 3.0); sec, _ = np.histogram(alts, bins=edges); sec = sec.astype(float)
Hp = out["Hp_max_+0"]; De = out["De_max_+0"]
p = (Hp * sec[None, :]).sum(1) * 1e3; e = (De * sec[None, :]).sum(1) * 1e3
table_17_6 = {2: (103, 15.8, 118), 5: (43.2, 2.84, 46), 10: (22.5, 1.55, 24), 20: (10.7, 0.691, 11.4), 30: (6.14, 0.312, 6.45), 50: (2.61, 0.0565, 2.66)}
print(f"\ncheck against Table 13.6 (2,000 m/s, 1 g drive, 3 m/s2 stop, {sec.sum()/3600:.2f} h on the tether; solar max, lon 0):")
print("  T g/cm2 | this run: protons, electrons, total mSv | elevator_dose.txt: protons, electrons, total mSv")
for i in range(len(T)):
    tv = int(T[i]) if T[i] in table_17_6 else None
    ref = table_17_6.get(tv)
    print(f"  {T[i]:6g}  | {p[i]:8.3g} {e[i]:8.3g} {p[i]+e[i]:8.3g}" + (f" | {ref[0]:8.3g} {ref[1]:8.3g} {ref[2]:8.3g}" if ref else ""))
