import numpy as np, sys, datetime, ppigrf, warnings
warnings.filterwarnings("ignore")
from sweetspot import U, thickness, earth_blocked, stormer
import gcr
ZD = gcr.Z_DEPTH; dOm = 4*np.pi/len(U)
RC = np.geomspace(5, 60, 14)
def gcr_tables(key):
    tab = np.array([gcr.dose_per_sr(r, key)[0] for r in RC])       # [Rc, depth]
    return tab
TAB = {k: gcr_tables(k) for k in ("H", "D")}
def gcr_at(t, rc, key="H"):
    """dose per sr for water-equivalent thickness t and cutoff rc (arrays). Beyond 590 g/cm2 the 400-590 g/cm2 slope is continued."""
    tab = TAB[key]; lr = np.log(np.clip(rc, RC[0], RC[-1])); i = np.clip(np.searchsorted(np.log(RC), lr)-1, 0, len(RC)-2)
    f = (lr-np.log(RC[i]))/(np.log(RC[i+1])-np.log(RC[i]))
    tt = np.clip(t, 1.0, 590.0); j = np.clip((tt/2.0-0.5).astype(int), 0, 298); g = tt/2.0-0.5-j
    v = lambda ii: tab[ii, j]*(1-g) + tab[ii, j+1]*g
    val = v(i)*(1-f) + v(i+1)*f
    # exponential continuation past 590 g/cm2 using the 400 to 590 slope at this cutoff
    a = tab[i, 200]*(1-f)+tab[i+1, 200]*f; b = tab[i, 295]*(1-f)+tab[i+1, 295]*f
    lam = np.maximum(190.0/np.log(np.maximum(a/b, 1.0001)), 50.0)
    return np.where(t > 590, val*np.exp(-(t-590)/lam), val)
def gcr_point(alt, t, decl=0.0, key="H", blocked_extra=None):
    rc = stormer(alt, decl); open_ = ~earth_blocked(alt)
    d = gcr_at(t, rc, key)*open_
    return d.sum()*dOm
if __name__ == "__main__":
    print("Species available:", [n for n in gcr.SPECIES if gcr.table(n, 'H') is not None])
    # --- validation: ISS-like point, 400 km, uniform shield, vertical cutoff from the same Stormer formula
    for tsh in (10, 20, 40, 60):
        t = np.full(len(U), float(tsh))
        print(f"ISS check, 400 km equator, uniform {tsh} g/cm2: absorbed {gcr_point(400, t, key='D')*1e3/365.25:.3f} mGy/day in water, dose equivalent {gcr_point(400, t)*1e3/365.25:.3f} mSv/day  (DOSTEL at cutoff > 10 GV: about 0.05 mGy/day in silicon, 0.06 in water)")
    # --- validation 2: the atmosphere at the equator (FAA DOT/FAA/AM-00/33 Table 1, 0 deg, 20 E): slab of X g/cm2 overhead
    uz = U[:,2]
    for km, X, faa in ((12.2, 188, 3.0), (9.1, 307, 1.6), (6.1, 472, 0.54)):
        t = np.where(uz > 0.02, X/np.maximum(uz, 0.02), 1e6); rc = stormer(km, 0.0)
        v = (gcr_at(t, rc)*(uz > 0.02)).sum()*dOm
        print(f"Atmosphere check, {km} km ({X} g/cm2 overhead): model {v*1e6/8766:.2f} microSv/h, FAA table {faa} microSv/h")
    for alt in (700, 800):
        print(f"\nRing at {alt} km, clean arc (no trapped protons). Galactic cosmic ray dose equivalent, mSv/yr")
        for tsh in (5, 20, 50, 100, 200, 300, 500):
            t = np.full(len(U), float(tsh)); print(f"  no cable, uniform shield {tsh:4d} g/cm2 water: {gcr_point(alt, t)*1e3:6.1f}  (absorbed {gcr_point(alt, t, key='D')*1e3:5.1f} mGy/yr)")
        print("  under the cable (room half-width 2 m, point 1 m above the floor):")
        print("   h below cable | sky: cable / Earth / open | thin walls 5 g/cm2 | walls 100 | walls 300 | walls+floor 500")
        for h in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
            row = []
            for tw in (5, 100, 300, 500):
                t, hit = thickness(h, t_wall=tw, t_floor=tw); row.append(gcr_point(alt, t)*1e3)
            eb = earth_blocked(alt)
            print(f"     {h:4.1f} m     |  {100*hit.mean():4.1f}% / {100*eb.mean():4.1f}% / {100*(~hit & ~eb).mean():4.1f}%  | " + " | ".join(f"{v:6.1f}" for v in row))
        for y0 in (1.5,):
            t, hit = thickness(1.0, t_wall=5, t_floor=5, y0=y0); print(f"   point 1.5 m off the centre line, h = 1 m, thin walls: {gcr_point(alt, t)*1e3:.1f}")
