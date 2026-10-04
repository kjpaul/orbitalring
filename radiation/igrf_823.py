import numpy as np, datetime, astropy.units as u, aep8, warnings
warnings.filterwarnings("ignore")
from igrf_sector import L_B
lons = np.arange(-177.5, 180, 5.0); t = aep8.model("p","min"); tm = aep8.model("p","max")
L, b, Bm = L_B(lons, 823.0, datetime.datetime(2026,9,1))
rows = np.load("poes_equator_rows.npy")
print("lon | AP-8 MIN with IGRF 2026 field, differential/4pi at 25, 50, 100 MeV | AP-8 MAX same | MEPED measured Aug-Sep 2026")
for i, lo in enumerate(lons):
    v = [np.nan_to_num(m.differential_flux_for_geomagnetic_coordinates(L[i], b[i], E*u.MeV).value)/(4*np.pi) for m in (t, tm) for E in (25,50,100)]
    print(f"{lo:7.1f} | {v[0]:7.3f} {v[1]:7.3f} {v[2]:7.3f} | {v[3]:7.3f} {v[4]:7.3f} {v[5]:7.3f} | {rows[i,1]:7.3f} {rows[i,2]:7.3f} {rows[i,3]:7.3f}")
