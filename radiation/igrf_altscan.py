"""AP-8 on the equator with the present-day field (IGRF 2026): L and B/B0 by field-line tracing, then the AP-8 maps.
Checked against MetOp-C MEPED at 823 km (poes_equator.py, igrf_823.py)."""
import numpy as np, datetime, astropy.units as u, aep8, warnings
warnings.filterwarnings("ignore")
from igrf_sector import L_B
lons = np.arange(-180, 180, 1.0); Eg = np.array([10,15,20,30,40,50,60,80,100,125,150,175,200,250,300.])
out = {}
for alt in (500,550,600,650,700,750,800,850,900,925,950,1000):
    L, b, Bm = L_B(lons, float(alt), datetime.datetime(2026,1,1))
    out[f"L_{alt}"] = L; out[f"b_{alt}"] = b
    for solar in ("min","max"):
        m = aep8.model("p", solar)
        I = np.array([np.nan_to_num(m.integral_flux_for_geomagnetic_coordinates(L, b, E*u.MeV).value) for E in Eg]); out[f"{solar}_{alt}"] = I
        J = I[0]; nz = lons[J > 0]
        e = aep8.model("e", "max"); Je = np.nan_to_num(e.integral_flux_for_geomagnetic_coordinates(L, b, 0.5*u.MeV).value) if solar == "max" else None
        print(f"{alt:5d} km solar {solar}: >10 MeV peak {J.max():8.1f} at {lons[J.argmax()]:5.0f}, mean {J.mean():7.1f}, nonzero {len(nz):3d} deg ({(nz.min() if len(nz) else 0):.0f}..{(nz.max() if len(nz) else 0):.0f}), clean {100*(J==0).mean():5.1f}%, >10/cm2/s over {(J>10).sum():3d} deg, >100 over {(J>100).sum():3d} deg, >1000 over {(J>1000).sum():3d} deg; >100 MeV peak {I[8].max():7.1f}" + (f"; electrons >0.5 MeV (AE-8 MAX) peak {Je.max():.0f}, mean {Je.mean():.0f}, nonzero {(Je>0).sum()} deg" if Je is not None else ""), flush=True)
np.savez("ap8_igrf2026_spectra.npz", Eg=Eg, lons=lons, **out)
