"""Round 13, item 6: the ambient air's pressure at 250 km (and, for comparison, at 700 km), from NRLMSIS as altitude_trade.py runs it.

    python tools/round13_air_250km.py

The atmosphere is called exactly as altitude_trade.msis() calls it: pymsis (NRLMSIS 2.1), the equator at longitude 0 on 2020-03-21, eight times
through the day (every 3 hours) and averaged, Ap = 15, F10.7 and its 81-day mean set equal.  Three cases: solar minimum (F10.7 = 70), the
cycle's mean (the average over 22 phases of F10.7 = 70 + 130 sin^2, from 70 to 200, mean 135), solar maximum (F10.7 = 230).
altitude_trade.py keeps only the mass density; this keeps the number densities and the temperature too, and gives p = n k T with n the sum of
the species' number densities (N2, O2, O, He, H, Ar, N, anomalous O, NO).  Each quantity is averaged over the day (and the cycle) by itself;
the pressure is averaged as n k T of each call, not as the product of the averages.
"""
import json
import math
import os

import numpy as np
import pymsis
from pymsis import Variable as V

K_B = 1.380649e-23
ALT = np.array([250.0, 700.0])
SPECIES = (V.N2, V.O2, V.O, V.HE, V.H, V.AR, V.N, V.ANOMALOUS_O, V.NO)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs", "round13_logs")


def call(F):
    """Day mean at F10.7 = F: rows of (pressure Pa, number density 1/m^3, temperature K, mass density kg/m^3) for each altitude."""
    r = []
    for hr in range(0, 24, 3):
        o = np.squeeze(pymsis.calculate(np.array([np.datetime64("2020-03-21T00:00") + np.timedelta64(hr, "h")]), 0.0, 0.0, ALT, F, F, [[15] * 7]))
        n = np.nansum(np.stack([o[..., int(s)] for s in SPECIES]), axis=0)
        T = o[..., int(V.TEMPERATURE)]
        r.append(np.stack([n * K_B * T, n, T, o[..., int(V.MASS_DENSITY)]], axis=-1))
    return np.mean(r, axis=0)


def main():
    phase = (np.arange(22) + 0.5) / 22 * np.pi
    cases = (("solar minimum, F10.7 = 70", call(70.0)), ("the cycle's mean, F10.7 from 70 to 200", np.mean([call(70 + 130 * math.sin(p) ** 2) for p in phase], axis=0)),
             ("solar maximum, F10.7 = 230", call(230.0)))
    print(f"pymsis {pymsis.__version__} (NRLMSIS 2.1), equator, 2020-03-21, mean of eight times of day, Ap = 15")
    print(f"{'altitude':>9} {'case':<42} {'pressure n k T, Pa':>19} {'number density, 1/m^3':>22} {'temperature, K':>15} {'mass density, kg/m^3':>21}")
    rows = []
    for i, h in enumerate(ALT):
        for name, c in cases:
            p, n, T, rho = (float(x) for x in c[i])
            rows.append(dict(alt_km=float(h), case=name, p_Pa=p, n_m3=n, T_K=T, rho_kg_m3=rho))
            print(f"{h:7.0f} km {name:<42} {p:19.3e} {n:22.3e} {T:15.0f} {rho:21.3e}")
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(pymsis=pymsis.__version__, rows=rows), open(os.path.join(OUT, "round13_air_250km.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
