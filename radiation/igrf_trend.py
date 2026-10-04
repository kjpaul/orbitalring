"""How the anomaly on the equator has changed with the field: AP-8 MIN read with IGRF at four epochs."""
import numpy as np, datetime, astropy.units as u, aep8, warnings
warnings.filterwarnings("ignore")
from igrf_sector import L_B
lons = np.arange(-180, 180, 2.0); m = aep8.model("p", "min")
for alt in (600, 650, 700, 800):
    for yr in (1960, 1980, 2000, 2026):
        L, b, Bm = L_B(lons, float(alt), datetime.datetime(yr, 1, 1))
        J = np.nan_to_num(m.integral_flux_for_geomagnetic_coordinates(L, b, 10*u.MeV).value)
        c = (lons*J).sum()/J.sum() if J.sum() > 0 else float('nan')
        print(f"{alt} km, field of {yr}: peak {J.max():7.0f}, ring mean {J.mean():6.1f}, clean {100*(J==0).mean():5.1f}%, centre {c:6.1f}, min B {Bm.min()*1e5:.0f} nT", flush=True)
