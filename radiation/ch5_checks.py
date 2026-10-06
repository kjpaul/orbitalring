"""Chapter 5 checks, 2026-10-06: (1) low-energy trapped protons at 600 km, (2) atomic oxygen fluence,
(3) dose for one trolley trip on an anchor line. Run from the radiation folder. Needs aep8, pymsis, astropy."""
import numpy as np, math, re, astropy.units as u, aep8, pymsis, warnings
warnings.filterwarnings("ignore")
YR = 365.25*86400
d = np.load('ap8_igrf2026_spectra.npz'); lons = d['lons']
try:
    print("1. Trapped protons 0.1 to 1 MeV, AP-8 with the IGRF 2026 field, equator")
    for alt in (600, 700):
        L, b = d[f'L_{alt}'], d[f'b_{alt}']
        for solar in ('min', 'max'):
            m = aep8.model('p', solar)
            J = {E: np.nan_to_num(m.integral_flux_for_geomagnetic_coordinates(L, b, E*u.MeV).value) for E in (0.1, 1.0, 10.0)}
            band = J[0.1] - J[1.0]; i = int(np.argmax(band))
            print(f"  {alt} km solar {solar}: worst lon {lons[i]:.0f}, J>0.1 {J[0.1][i]:.3g}, J>1 {J[1.0][i]:.3g}, J>10 {J[10.0][i]:.3g} /cm2/s; "
                  f"0.1-1 MeV {band[i]:.3g} /cm2/s = {band[i]*YR*30:.3g} /cm2 in 30 yr; ring mean {band.mean()*YR*30:.3g}")
except AttributeError as e:
    print("  skipped: this needs the aep8 build with aep8.model (L, B/B0 interface); installed version lacks it:", e)
print("2. Atomic oxygen, NRLMSIS 2.1, equator, F10.7 cycle 70 to 200 (same cycle as altitude_trade.py)")
kB = 1.380649e-23; mO = 16*1.66054e-27; GM = 3.986004418e14; RE = 6378137.0; W = 7.2921e-5
def msis(alt, F):
    n = []; T = []
    for hr in range(0, 24, 3):
        o = pymsis.calculate(np.array([np.datetime64('2020-03-21T00:00') + np.timedelta64(hr, 'h')]), 0., 0., alt, F, F, [[15]*7])
        o = np.squeeze(o); n.append(o[pymsis.Variable.O]); T.append(o[pymsis.Variable.TEMPERATURE])
    return np.mean(n), np.mean(T)
phase = np.linspace(0, math.pi, 23)[:-1]
for alt in (400, 600, 700):
    r = RE + alt*1e3; vorb = math.sqrt(GM/r) - W*r
    res = {}
    for name, Fs in (("cycle average", [70 + 130*math.sin(p)**2 for p in phase]), ("solar min F=70", [70]), ("solar max F=200", [200])):
        nT = [msis(alt, F) for F in Fs]
        n = np.mean([x[0] for x in nT]); ram = np.mean([x[0] for x in nT])*vorb
        th = np.mean([x[0]*math.sqrt(8*kB*x[1]/(math.pi*mO))/4 for x in nT])
        print(f"  {alt} km {name}: n(O) {n:.3g} /m3; ram flux on an orbiting surface {ram*1e-4:.3g} /cm2/s = {ram*1e-4*YR:.3g} /cm2/yr; "
              f"thermal flux on a surface at rest in the air {th*1e-4:.3g} /cm2/s = {th*1e-4*YR:.3g} /cm2/yr")
print("3. Trolley trip on an anchor line (trapped protons at the worst longitude, cosmic rays anywhere)")
alts = []; pk = []
for line in open('igrf_altscan_clean.txt'):
    m = re.match(r"\s*(\d+) km solar min: >10 MeV peak\s+([\d.]+)", line)
    if m and int(m.group(1)) <= 700: alts.append(int(m.group(1))); pk.append(float(m.group(2)))
alts = [450] + alts; pk = [0.0] + pk     # the model gives 114 /cm2/s at 500 km; taken as zero at 450 km
eq_km = np.trapz(pk, alts)/pk[-1]
print(f"  peak flux by altitude {dict(zip(alts, pk))}; integral = {eq_km:.1f} km at the 700 km flux")
RATE_OPEN = 65.0   # Sv/yr at 700 km, peak longitude, walls equal to 5 cm of water (Table 5.8)
GCR = (32.0, 40.0) # mSv/yr at 700 km in an ordinary room (Table 5.7), used for the whole trip as an upper bound
for v in (100.0, 130.0):
    hrs = 700/v; eq_h = eq_km/v
    print(f"  {v:.0f} km/h: trip {hrs:.1f} h; trapped protons, worst longitude {RATE_OPEN*1e3/8766*eq_h:.1f} mSv per trip "
          f"({eq_h:.2f} h equivalent at the 700 km rate of {RATE_OPEN*1e3/8766:.1f} mSv/h); cosmic rays under {GCR[0]/8766*hrs*1e3:.0f} to {GCR[1]/8766*hrs*1e3:.0f} microSv per trip")
for v_top in (300.0, 500.0):
    eq_h = eq_km/v_top
    print(f"  if the top 250 km is covered at {v_top:.0f} km/h: trapped protons at the worst longitude {RATE_OPEN*1e3/8766*eq_h:.1f} mSv per trip")
