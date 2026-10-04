"""Trapped-proton dose behind a water shield: AP-8 spectra folded with the Geant4 depth-dose runs.
Dose at the centre of a solid water sphere of radius t under isotropic flux (same geometry as SHIELDOSE's solid sphere).
Includes nuclear interactions and all secondaries (neutrons, recoil protons, alphas, heavy recoils, gammas)."""
import numpy as np
from response import load, R_p
GY = 1.602e-10; YR = 3.156e7
low = load("out_trap_lowE/*.npz"); main = load("out_trap/*.npz")
runs = []
for r in low:                       # thin slabs: depth = 1.25 x range, 250 bins
    w = 1.25*float(R_p(r["E"]))/250; z = (np.arange(250)+0.5)*w
    runs.append(dict(E=r["E"], z=z, D=r["D"]/w, H=r["H"]/w, R=float(R_p(r["E"])), zmax=250*w))
Elow_max = max(r["E"] for r in runs)
for r in main:                      # 600 cm slab, 2 cm bins
    if r["E"] <= Elow_max: continue
    z = (np.arange(300)+0.5)*2.0
    runs.append(dict(E=r["E"], z=z, D=r["D"]/2.0, H=r["H"]/2.0, R=float(R_p(r["E"])), zmax=600.0))
runs.sort(key=lambda r: r["E"]); ER = np.array([r["E"] for r in runs])
def prof(r, key, z):                # MeV/g per (proton/cm2) at depth z (g/cm2), 0 outside the scored slab
    return np.interp(z, r["z"], r[key], left=r[key][0], right=0.0) * (z <= r["zmax"])
def resp(E, t, key):
    """response at energy E and depth t by interpolating between the two neighbouring runs at the same depth/range ratio"""
    i = np.clip(np.searchsorted(ER, E)-1, 0, len(runs)-2); a, b = runs[i], runs[i+1]
    f = (np.log(E)-np.log(a["E"]))/(np.log(b["E"])-np.log(a["E"])); R = float(R_p(E)); x = t/R
    if x <= 3.0: va, vb = prof(a, key, x*a["R"]), prof(b, key, x*b["R"])
    else:        va, vb = prof(a, key, t), prof(b, key, t)
    return (1-f)*va + f*vb
EF = np.geomspace(10.0, 2000.0, 900)
T = np.array([0.5,1,2,3,5,7.5,10,15,20,30,40,50,75,100,150,200,300,400,500])
RESP = {k: np.array([[resp(E, t, k) if E >= ER[0] else 0.0 for E in EF] for t in T]) for k in ("D","H")}
def spectrum(Eg, I, tail):
    """differential flux on EF from the integral spectrum I(Eg). tail: 'cut' = all flux above the last AP-8 point lies within 4/3 of it (300 to 400 MeV);
    'exp' = continue the last e-folding to 2 GeV."""
    I = np.nan_to_num(np.asarray(I, float)); ok = I > 0
    if ok.sum() < 2: return np.zeros(len(EF)-1), 0.5*(EF[1:]+EF[:-1])
    Eg, I = Eg[ok], I[ok]
    li = np.interp(EF, Eg, np.log(I), right=-np.inf)
    if tail == "exp":
        k = (np.log(I[-1])-np.log(I[-2]))/(Eg[-1]-Eg[-2]); m = EF > Eg[-1]; li[m] = np.log(I[-1]) + k*(EF[m]-Eg[-1])
    If = np.exp(li); n = -np.diff(If); n[n < 0] = 0
    if tail == "cut":                # AP-8 ends at 400 MeV: all flux above the last tabulated point is put between that
        j = np.searchsorted(EF, Eg[-1])-1; hi = np.searchsorted(EF, Eg[-1]*4/3)   # point and 4/3 of it (300 -> 400 MeV)
        k = (np.log(I[-1])-np.log(I[-2]))/(Eg[-1]-Eg[-2]); wgt = np.exp(k*(EF[j:hi]-Eg[-1])); n[j:] = 0; n[j:hi] = If[j]*wgt/wgt.sum()
    return n, None
def dose(Eg, I, tail="exp"):
    """returns absorbed dose (Gy/yr) and dose equivalent (Sv/yr) at the centre of water spheres of radius T (g/cm2)"""
    n, _ = spectrum(Eg, I, tail)
    out = []
    for k in ("D", "H"):
        r = 0.5*(RESP[k][:,1:]+RESP[k][:,:-1]); out.append((r*n[None,:]).sum(1)*GY*YR)
    return out
if __name__ == "__main__":
    print("runs:", len(runs), "energies", ER[0], "to", ER[-1])
    # check 1: entrance dose per unit fluence vs PSTAR stopping power; check 2: neutron/secondary tail beyond the range
    from response import S_p
    for E in (50, 100, 200, 400, 1000):
        print(f"E={E}: response at t=0.5 g/cm2 {resp(E,0.5,'D'):.2f} MeV/g per p/cm2 (PSTAR S={float(S_p(E)):.2f}); at 1.2 x range {resp(E,1.2*float(R_p(E)),'D'):.4f} abs, {resp(E,1.2*float(R_p(E)),'H'):.4f} eq; H/D at half range {resp(E,0.5*float(R_p(E)),'H')/resp(E,0.5*float(R_p(E)),'D'):.2f}")
