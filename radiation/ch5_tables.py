"""Chapter 5 tables from the present-field AP-8 spectra: electronics dose behind aluminium, solar panel life."""
import numpy as np
from fold_trapped import dose, T
d = np.load("ap8_igrf2026_spectra.npz"); Eg = d["Eg"]; lons = d["lons"]; YR = 3.156e7
AL_TO_WATER = 7.718/10.01      # PSTAR CSDA range ratio water/aluminium at 100 MeV
SI_OVER_WATER = 0.80           # NIST PSTAR stopping power silicon/water: 0.783 at 30 MeV, 0.801 at 100 MeV, 0.810 at 300 MeV
print("Proton dose behind aluminium (trapped protons above 10 MeV only; trapped electrons NOT included), krad(Si) per year")
print(" altitude | solar | 3 mm Al worst / mean | 10 mm worst / mean | 30 mm worst / mean | 100 mm worst / mean")
for alt in (600, 700, 750, 800):
    for solar in ("min", "max"):
        I = np.nan_to_num(d[f"{solar}_{alt}"])
        D = np.array([dose(Eg, I[:,k], "exp")[0] for k in range(len(lons))])
        row = []
        for mm in (3, 10, 30, 100):
            tw = mm*0.1*2.70*AL_TO_WATER; v = np.array([np.exp(np.interp(tw, T, np.log(np.maximum(x, 1e-12)))) for x in D])*100*SI_OVER_WATER/1000
            row.append(f"{v.max():8.2f} / {v.mean():7.3f}")
        print(f"  {alt} km | {solar} | " + " | ".join(row))
print("\nSolar panels: annual fluence above 10 MeV and time to 10% and 20% power loss.")
print("Thresholds 9e10 and 4e11 p/cm2 (triple junction, 4 mil coverglass; from the ZTJ datasheet, see solar_cal.py)")
print(" altitude | solar | worst longitude: fluence/yr, yr to 10%, yr to 20% | ring mean: fluence/yr, yr to 10%, yr to 20% | share of ring reaching 10% loss within 30 yr")
for alt in (500, 550, 600, 650, 700, 750, 800):
    for solar in ("min", "max"):
        J = np.nan_to_num(d[f"{solar}_{alt}"][0]); F = J*YR
        print(f"  {alt} km | {solar} | {F.max():.2e}, {9e10/F.max():6.1f}, {4e11/F.max():6.1f} | {F.mean():.2e}, {9e10/F.mean():6.1f}, {4e11/F.mean():6.1f} | {100*(F*30 > 9e10).mean():.0f}%")
