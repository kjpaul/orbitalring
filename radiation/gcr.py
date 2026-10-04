"""Galactic cosmic ray dose behind shielding at a given geomagnetic cutoff, from Geant4 depth-dose runs (FTFP_BERT)
folded with the Particle Data Group primary spectrum:
  I_N(E) = 1.8e4 E^-2.7 nucleons/(m2 s sr GeV), E = total energy per nucleon in GeV; 74% of nucleons are free protons,
  70% of the rest are in helium; heavier nuclei split by the PDG abundance table at 10.6 GeV/nucleon (F x A weights).
Dose is per steradian of open sky; the geometry module multiplies by the open solid angle and the direction-dependent cutoff."""
import numpy as np
from response import load
MP = 0.938; GY = 1.602e-10; YR = 3.156e7
SPECIES = {   # name: (file glob, A, Z, share of primary nucleons, Q of unidentified ion fragments)
 "H":  ("out_gcr/proton_*.npz",    1,  1, 0.74,              20.0),
 "He": ("out_gcr/alpha_*.npz",     4,  2, 0.26*0.70,         20.0),
 "O":  ("out_gcr/ion_8_16_*.npz", 16,  8, 0.26*0.30*0.551,   10.0),
 "Si": ("out_gcr/ion_14_28_*.npz",28, 14, 0.26*0.30*0.282,   10.0),
 "Fe": ("out_gcr/ion_26_56_*.npz",56, 26, 0.26*0.30*0.167,   10.0),
}
Z_DEPTH = (np.arange(300)+0.5)*2.0
def table(name, key):
    g, A, Z, share, qres = SPECIES[name]
    runs = load(g, q_resid=qres, prim_Z=(Z if Z > 2 else None), A=A)
    if not runs: return None
    runs.sort(key=lambda r: r["E"])
    En = np.array([r["E"] for r in runs])/A/1000.0           # kinetic energy per nucleon, GeV
    R = np.array([r[key]/2.0 for r in runs])                # MeV/g per (particle/cm2), [energy, depth]
    return En, R
def dose_per_sr(Rc, key="H", species=None, smooth=5):
    """dose (Gy/yr or Sv/yr) per steradian of sky versus depth in water, for rigidity cutoff Rc (GV)"""
    tot = np.zeros(300); parts = {}
    for name in (species or SPECIES):
        tb = table(name, key)
        if tb is None: continue
        En, R = tb; g, A, Z, share, _ = SPECIES[name]
        pn = Rc*Z/A; Ecut = np.sqrt(pn**2+MP**2) - MP           # kinetic energy per nucleon at the cutoff
        Ef = np.geomspace(max(Ecut, 0.5), 2e4, 400); Em = np.sqrt(Ef[1:]*Ef[:-1]); dE = np.diff(Ef)
        flux = share*1.8e4*(Em+MP)**-2.7/A*1e-4                 # nuclei /(cm2 s sr GeV per nucleon)
        logR = np.log(np.maximum(R, 1e-12))
        # log-log interpolation in energy at each depth; power-law continuation outside the simulated energies
        le = np.log(En); x = np.log(Em); i = np.clip(np.searchsorted(le, x)-1, 0, len(le)-2)
        f = (x-le[i])/(le[i+1]-le[i])
        r = np.exp(logR[i]*(1-f)[:,None] + logR[i+1]*f[:,None])
        d = (flux[:,None]*dE[:,None]*r).sum(0)*GY*YR
        parts[name] = d; tot += d
    if smooth > 1:
        k = np.ones(smooth)/smooth; tot = np.convolve(np.pad(tot, smooth//2, mode="edge"), k, mode="valid")
    return tot, parts
if __name__ == "__main__":
    for Rc in (8.1, 11.8, 14.9, 47.0):
        H, ph = dose_per_sr(Rc, "H"); D, pd = dose_per_sr(Rc, "D")
        print(f"\ncutoff {Rc} GV: per steradian. depth g/cm2: dose equivalent mSv/yr/sr, absorbed mGy/yr/sr, share by species (H-eq)")
        for t in (2, 10, 20, 40, 60, 100, 150, 200, 300, 400, 500):
            j = int(t/2)
            print(f"  {t:4d}: {1e3*H[j]:7.2f} {1e3*D[j]:7.2f}   " + " ".join(f"{n} {100*ph[n][j]/sum(p[j] for p in ph.values()):.0f}%" for n in ph))
