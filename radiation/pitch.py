"""Local pitch-angle distribution of trapped protons from AP-8 itself.
Along a field line (fixed L) AP-8 gives the omnidirectional flux J(b), b = B/B0. With Liouville's theorem the
unidirectional flux at pitch angle a where the field is b equals the 90-degree flux at the mirror point
b_m = b / sin^2(a). So J(b) = 4 pi * integral_0^1 j90(b/(1-mu^2)) dmu, which is inverted here for j90."""
import numpy as np, astropy.units as u, aep8
from scipy.optimize import nnls
def pitch_pdf(L, b, E=30.0, solar="min", n=60):
    mdl = aep8.model("p", solar)
    J = lambda bb: float(np.nan_to_num(mdl.integral_flux_for_geomagnetic_coordinates(L, bb, E*u.MeV).to_value(u.cm**-2/u.s)))
    if J(b) <= 0: return None
    hi = b
    while J(hi) > 0 and hi < 50*b: hi *= 1.02          # find cutoff b_c
    edges = np.linspace(b, hi, n+1); mid = 0.5*(edges[1:]+edges[:-1])
    # j90 piecewise constant on mirror-field cells; row i: J at b_i = edges[i]
    A = np.zeros((n, n)); y = np.zeros(n)
    for i in range(n):
        bi = edges[i]; y[i] = J(bi)
        for k in range(i, n):
            mu_lo = np.sqrt(1 - bi/edges[k]) ; mu_hi = np.sqrt(1 - bi/edges[k+1])
            A[i, k] = 4*np.pi*(mu_hi - mu_lo)
    j90, _ = nnls(A, y)
    # local distribution at b: j(mu) = j90(b/(1-mu^2)); cells in mu
    mu_edges = np.sqrt(1 - b/edges)
    return mu_edges, j90, hi/b      # j constant between mu_edges[k], mu_edges[k+1]
if __name__ == "__main__":
    d = np.load("ap8_spectra.npz")
    for alt in (700, 800):
        I = d[f"min_{alt}"][0]
        for lon in (-60,-45,-27,-10,5):
            k = int(np.where(d["lons"]==lon)[0][0]); L = d[f"L_min_{alt}"][k]; b = d[f"BB0_min_{alt}"][k]
            for E in (30.,100.):
                r = pitch_pdf(L, b, E)
                if r is None: print(alt, lon, E, "zero flux"); continue
                mu, j, ratio = r
                w = j*np.diff(mu); w /= w.sum()
                cum = np.cumsum(w)
                a50 = 90-np.degrees(np.arccos(np.interp(0.5,cum,mu[1:]))) if False else None
                m50 = np.interp(0.5, cum, mu[1:]); m90 = np.interp(0.9, cum, mu[1:]); mmax = mu[1:][w>1e-6].max()
                print(f"{alt} km lon {lon:4d} E>{E:.0f}: L={L:.3f} B/B0={b:.3f} cutoff/b={ratio:.2f}; half of flux within {np.degrees(np.arcsin(m50)):.0f} deg of 90, 90% within {np.degrees(np.arcsin(m90)):.0f} deg, edge {np.degrees(np.arcsin(mmax)):.0f} deg")
