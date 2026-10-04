"""Present-day anomaly sector on the equator from a current field model (IGRF via ppigrf).
1. Field direction (declination, dip) and strength along the geographic equator.
2. McIlwain L and B/B0 by field-line tracing in IGRF, checked against the AP-8 package's own field (Jensen-Cain 1960).
3. AP-8 MIN flux evaluated with present-day (L, B/B0): an UPPER BOUND (ECSS/Heynderickx warn that AP-8 with a
   modern field overestimates), next to the standard 0.3 deg/yr westward shift."""
import numpy as np, ppigrf, datetime, astropy.units as u, aep8, sys
RE = 6371.2; M = 0.311653  # Gauss Re^3, the dipole moment constant used with AP-8 for B0 = M/L^3
def bvec(xyz, date):
    r = np.linalg.norm(xyz, axis=-1); th = np.degrees(np.arccos(xyz[...,2]/r)); ph = np.degrees(np.arctan2(xyz[...,1], xyz[...,0]))
    Br, Bt, Bp = [np.asarray(a).reshape(r.shape) for a in ppigrf.igrf_gc(r, th, ph, date)]
    t = np.radians(th); p = np.radians(ph)
    er = np.stack([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)], -1)
    et = np.stack([np.cos(t)*np.cos(p), np.cos(t)*np.sin(p), -np.sin(t)], -1)
    ep = np.stack([-np.sin(p), np.cos(p), 0*p], -1)
    return (Br[...,None]*er + Bt[...,None]*et + Bp[...,None]*ep)*1e-5   # Gauss
def L_B(lons, alt, date, ds=15.0, nmax=1500):
    a = 6378.137; p = np.radians(lons)
    x0 = np.stack([(a+alt)*np.cos(p), (a+alt)*np.sin(p), 0*p], -1)
    B0v = bvec(x0, date); Bm = np.linalg.norm(B0v, axis=-1)
    I = np.zeros(len(lons))
    for sgn in (+1, -1):
        x = x0.copy(); alive = np.ones(len(lons), bool); Bprev = Bm.copy()
        for _ in range(nmax):
            def f(y):
                b = bvec(y, date); return sgn*b/np.linalg.norm(b, axis=-1, keepdims=True), np.linalg.norm(b, axis=-1)
            k1, B1 = f(x); k2, _ = f(x+0.5*ds*k1); k3, _ = f(x+0.5*ds*k2); k4, _ = f(x+ds*k3)
            xn = x + ds*(k1+2*k2+2*k3+k4)/6; Bn = np.linalg.norm(bvec(xn, date), axis=-1)
            inside = alive & (Bn < Bm)          # still between the mirror points
            # first step: decide if this direction goes toward lower B at all
            cross = alive & ~inside
            # partial last step: fraction until B reaches Bm
            frac = np.where(cross & (Bn != Bprev), np.clip((Bm-Bprev)/(Bn-Bprev), 0, 1), 0.0)
            g_prev = np.sqrt(np.clip(1-Bprev/Bm, 0, None)); g_new = np.sqrt(np.clip(1-Bn/Bm, 0, None))
            I += np.where(inside, 0.5*(g_prev+g_new)*ds, 0.0) + np.where(cross, 0.5*g_prev*frac*ds, 0.0)
            x = np.where(inside[:,None], xn, x); Bprev = np.where(inside, Bn, Bprev); alive = inside
            if not alive.any(): break
    I = I/RE
    X = np.log(np.maximum(I**3*Bm/M, 1e-300))
    # Hilton (1971) approximation
    X3 = (I**3*Bm/M)
    L = ((1 + 1.35047*X3**(1/3) + 0.465376*X3**(2/3) + 0.0475455*X3) * M/Bm)**(1/3)
    return L, Bm/(M/L**3), Bm
if __name__ == "__main__":
    lons = np.arange(-180, 180, 1.0); t = aep8.model("p","min"); out = {}
    from astropy.coordinates import EarthLocation; from astropy.time import Time
    for alt in (600, 700, 750, 800):
        loc = EarthLocation.from_geodetic(lons*u.deg, 0*u.deg, alt*u.km)
        Lp, bp = t.geomagnetic_coordinates(loc, Time("2026-01-01")); Lp = np.asarray(Lp); bp = np.asarray(bp)
        J_pkg = np.nan_to_num(t.integral_flux_for_geomagnetic_coordinates(Lp, bp, 10*u.MeV).value)
        for yr in (1960, 2026):
            L, b, Bm = L_B(lons, alt, datetime.datetime(yr,1,1))
            J = np.nan_to_num(t.integral_flux_for_geomagnetic_coordinates(L, b, 10*u.MeV).value)
            nz = lons[J > 0]; out[(alt,yr)] = (L,b,Bm,J)
            print(f"{alt} km IGRF {yr}: min |B| {Bm.min():.4f} G at lon {lons[Bm.argmin()]:.0f}; AP-8 MIN >10 MeV with this field: peak {J.max():.0f} at {lons[J.argmax()]:.0f}, nonzero {nz.min():.0f}..{nz.max():.0f} ({len(nz)} deg), ring mean {J.mean():.1f}")
            if yr == 1960:
                print(f"     check against package field (Jensen-Cain 1960): max |dL| {np.abs(L-Lp).max():.4f}, max |d(B/B0)| {np.abs(b-bp).max():.4f}; package: peak {J_pkg.max():.0f} at {lons[J_pkg.argmax()]:.0f}, nonzero {lons[J_pkg>0].min():.0f}..{lons[J_pkg>0].max():.0f}")
        # shift that best maps the 1960 profile onto the 2026 one
        J0 = out[(alt,1960)][3]; J1 = out[(alt,2026)][3]
        if J0.max() > 0:
            c0 = (lons*J0).sum()/J0.sum(); c1 = (lons*J1).sum()/J1.sum()
            print(f"     flux-weighted centre longitude 1960 {c0:.1f}, 2026 {c1:.1f}: shift {c1-c0:.1f} deg in 66 yr = {(c1-c0)/66:.2f} deg/yr")
    print("\nField direction on the equator at 800 km, IGRF 2026 (declination D east of north, dip I positive down):")
    for lon in range(-95, 25, 10):
        Be, Bn, Bu = [float(np.ravel(v)[0]) for v in ppigrf.igrf(lon, 0.0, 800.0, datetime.datetime(2026,1,1))]
        H = np.hypot(Be, Bn); print(f"  lon {lon:4d}: |B| {np.sqrt(H*H+Bu*Bu)*1e-5:.4f} G, D {np.degrees(np.arctan2(Be,Bn)):6.1f}, I {np.degrees(np.arctan2(-Bu,H)):6.1f}")
    np.savez("igrf_LB.npz", lons=lons, **{f"{k}_{a}_{y}": v[i] for (a,y),v in out.items() for i,k in enumerate(("L","b","B","J"))})
