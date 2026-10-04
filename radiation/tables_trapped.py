import numpy as np, sys
from fold_trapped import dose, T
src = sys.argv[1] if len(sys.argv) > 1 else "epoch"
d = np.load("ap8_spectra.npz" if src == "epoch" else "ap8_igrf2026_spectra.npz"); Eg = d["Eg"][:15]; lons = d["lons"]
show = [1, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300]
idx = [int(np.where(T == s)[0][0]) for s in show]
print(f"Field: {'AP-8 own field (1960s epoch)' if src=='epoch' else 'IGRF 2026'}; Geant4 QGSP_BIC_HP response; solid water sphere; isotropic flux.")
for alt in (600, 700, 750, 800):
    for solar in ("min",):
        I = np.nan_to_num(d[f"{solar}_{alt}"][:15])
        for tail in ("cut", "exp"):
            res = np.array([dose(Eg, I[:,k], tail) for k in range(len(lons))])     # lon, (D,H), T
            D, H = res[:,0,:], res[:,1,:]; k = I[0].argmax()
            print(f"\n{alt} km, solar {solar}, spectrum above 300 MeV: {'none (AP-8 limit, flux above 300 MeV placed at 300 to 400 MeV)' if tail=='cut' else 'exponential continuation to 2 GeV'}; flux peak at longitude {lons[k]:.0f}")
            print(" shield g/cm2 | worst: Gy/yr  Sv/yr | ring mean: Gy/yr  Sv/yr | share of ring above 20 mSv/yr | above 1 mSv/yr")
            for s, i in zip(show, idx):
                print(f"   {s:5d}      | {D[:,i].max():9.4f} {H[:,i].max():9.4f} | {D[:,i].mean():9.5f} {H[:,i].mean():9.5f} | {100*(H[:,i]>0.02).mean():5.1f}% | {100*(H[:,i]>0.001).mean():5.1f}%")
