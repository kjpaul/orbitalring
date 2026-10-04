"""Turns the Geant4 depth-dose runs into response functions and folds them with spectra."""
import numpy as np, glob, re
PS = np.loadtxt("pstar_water.txt")           # NIST PSTAR liquid water: E, S_el, S_nuc, S_tot (MeV cm2/g), CSDA, projected, detour
def S_p(E):   return np.exp(np.interp(np.log(E), np.log(PS[:,0]), np.log(PS[:,3])))
def R_p(E):   return np.exp(np.interp(np.log(E), np.log(PS[:,0]), np.log(PS[:,4])))
def Q_of_L(L):                                # ICRP Publication 60 quality factor, L in keV/um in water
    L = np.asarray(L, float)
    return np.where(L < 10, 1.0, np.where(L <= 100, 0.32*L-2.2, 300/np.sqrt(np.maximum(L,1e-9))))
PBANDS = [0, 0.3, 1, 3, 10, 30, 100, 300, 1e9]; ABANDS = [0, 2, 8, 30, 100, 400, 1e9]
def band_Q(lo, hi, zeff2=1.0, mass=1.0, cap=None):
    hi = min(hi, 1e4); E = np.linspace(max(lo, 1e-3*mass), hi, 4000)      # energy-averaged Q of a particle slowing through the band
    L = zeff2*S_p(E/mass)*0.1
    if cap: L = np.minimum(L, cap)
    return float(Q_of_L(L).mean())
QP = [band_Q(PBANDS[i], PBANDS[i+1]) for i in range(8)]
QA = [band_Q(ABANDS[i], ABANDS[i+1], 4.0, 3.97, cap=226.0) for i in range(6)]
Q_DT, Q_RECOIL = 2.0, 20.0                     # assumptions: deuterons/tritons; heavy target recoils (ICRP 60 maximum is 30)
def load(dirglob, q_resid=Q_RECOIL, prim_Z=None, A=1):
    runs = []
    for f in sorted(glob.glob(dirglob)):
        d = np.load(f); names = list(d["names"]); dat = dict(zip(names, d["data"]/float(d["nev"])))
        nb = len(dat["tot"]); m = re.search(r"_([\d.]+)MeV", f); E = float(m.group(1))
        depth = float(d["depth"]) if "depth" in d else None
        known = sum(dat[n] for n in names if n not in ("tot", "prim"))
        prim = dat.get("prim", 0*dat["tot"]); resid = np.maximum(dat["tot"] - known - prim, 0)
        H = dat["em"] + dat["mupi"] + sum(QP[i]*dat[f"p{i}"] for i in range(8)) + Q_DT*dat["dt"] \
            + sum(QA[i]*dat[f"a{i}"] for i in range(6)) + q_resid*resid
        if prim_Z:
            Lp = prim_Z**2*S_p(min(E/A, 1e4))*0.1; H = H + float(Q_of_L(Lp))*prim
        runs.append(dict(E=E, D=dat["tot"], H=H, nb=nb, f=f))
    return runs
