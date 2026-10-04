"""Trapped-proton dose at a point under the cable, by direction, for longitudes in the anomaly sector.
Angular distribution: from AP-8 itself (pitch.py), symmetric about the plane perpendicular to the field.
Field direction: IGRF 2026. No east-west asymmetry is applied (it changes east against west, not above against below)."""
import numpy as np, datetime, ppigrf, warnings, sys
warnings.filterwarnings("ignore")
from sweetspot import U, thickness
from pitch import pitch_pdf
from fold_trapped import dose, T
d = np.load("ap8_igrf2026_spectra.npz"); Eg = d["Eg"]; lons = d["lons"]
def weights(alt, k):
    L = float(d[f"L_{alt}"][k]); b = float(d[f"b_{alt}"][k])
    r = pitch_pdf(L, b, 30.0)
    if r is None: return None
    mu_e, j, _ = r
    Be, Bn, Bu = [float(np.ravel(v)[0]) for v in ppigrf.igrf(float(lons[k]), 0.0, float(alt), datetime.datetime(2026,1,1))]
    bh = np.array([Be, Bn, Bu]); bh /= np.linalg.norm(bh)
    mu = np.abs(U @ bh); idx = np.searchsorted(mu_e, mu)-1
    w = np.where((idx >= 0) & (idx < len(j)), j[np.clip(idx, 0, len(j)-1)], 0.0)
    return w/w.sum(), bh
def point_dose(alt, k, t, tail="exp"):
    wb = weights(alt, k)
    if wb is None: return 0.0, 0.0
    w, bh = wb
    D, H = dose(Eg, np.nan_to_num(d[f"min_{alt}"][:, k]), tail)       # sphere-centre dose versus T
    lt = np.log(np.maximum(H, 1e-12)); Ht = np.exp(np.interp(np.clip(t, T[0], None), T, lt, right=-30))
    # beyond 500 g/cm2: continue the 300-500 slope
    sl = (lt[-1]-lt[-3])/(T[-1]-T[-3]); Ht = np.where(t > T[-1], np.exp(lt[-1] + sl*(t-T[-1])), Ht)
    return float((w*Ht).sum()), float(H[np.where(T == 5)[0][0]])
if __name__ == "__main__":
    for alt in (700, 800):
        J = d[f"min_{alt}"][0]; kpk = int(J.argmax())
        print(f"\nRing at {alt} km, IGRF 2026 field, AP-8 MIN, spectrum continued above 300 MeV. Trapped-proton dose equivalent, Sv/yr")
        print(" longitude | flux >10 MeV | open, 5 g/cm2 sphere | under cable h=1 m: walls+floor 5 | 100 | 300 g/cm2 | no cable, same room, 100 | cable + floor 300, walls 100")
        for lon in sorted(set([int(lons[kpk])] + list(range(-110, 21, 10)))):
            k = int(np.where(lons == lon)[0][0])
            if J[k] <= 0: continue
            row = []
            for tw, tf, cab in ((5,5,True), (100,100,True), (300,300,True), (100,100,False), (100,300,True)):
                t, hit = thickness(1.0, t_wall=tw, t_floor=tf, cable=cab)
                if not cab:                       # no cable: a ceiling slab of the wall thickness instead
                    ux, uy, uz = U.T; up = (uz > 0) & (np.abs(1.0*uy/np.where(uz > 0, uz, 1)) < 2.0); t = np.where(up, tw/np.maximum(uz, 1e-9), t)
                v, open5 = point_dose(alt, k, t); row.append(v)
            print(f"   {lon:5d}   | {J[k]:8.1f}    | {open5:9.3f}   | " + " | ".join(f"{v:8.4f}" for v in row))
        w, bh = weights(alt, kpk); t, hit = thickness(1.0)
        print(f"  at the flux peak: share of the trapped flux arriving through the cable {100*w[hit].sum():.0f}%, from below the horizontal {100*w[U[:,2]<0].sum():.0f}%, field direction (east, north, up) {np.round(bh,3)}")
