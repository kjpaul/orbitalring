"""Ring altitude trade, debris side: released-cable orbits and drag lifetimes for ring altitudes 500 to 800 km,
with three atmospheres: the Vallado mean table used in Chapter 4, NRLMSIS 2.1 averaged over a solar cycle,
and NRLMSIS held at solar maximum (a bound). Prompt soot from the collision ejecta from debris_simulation.py."""
import math, numpy as np, pymsis, warnings
warnings.filterwarnings("ignore")
import debris_orbits as do, debris_simulation as ds
from debris_orbits import released_cable, swept_cable_speed, orbit_from_state, ballistic_coefficient, R_EARTH, V_RECOIL, RHO_CNT, ring_state
ALT = np.arange(100, 1001, 25.0)
def msis(F):
    r = []
    for hr in range(0, 24, 3):
        o = pymsis.calculate(np.array([np.datetime64('2020-03-21T00:00') + np.timedelta64(hr, 'h')]), 0., 0., ALT, F, F, [[15]*7])
        r.append(np.squeeze(o)[..., 0])
    return np.mean(r, axis=0)
phase = (np.arange(22)+0.5)/22*np.pi
cyc = np.mean([msis(70 + 130*math.sin(p)**2) for p in phase], axis=0)     # F10.7 from 70 to 200, mean 135
smax = msis(230.0); smin = msis(70.0)
VALLADO = do.rho_atm
def make(tab):
    lt = np.log(tab); H_top = -(ALT[-1]-ALT[-5])/(lt[-1]-lt[-5])
    def rho(h_m):
        h = h_m/1e3
        if h < ALT[0]: return VALLADO(h_m)
        if h > ALT[-1]: return float(tab[-1]*math.exp(-(h-ALT[-1])/H_top))
        return float(math.exp(np.interp(h, ALT, lt)))
    return rho
ATM = {"Vallado mean (Chapter 4)": VALLADO, "MSIS, solar-cycle average": make(cyc), "MSIS, permanent solar maximum": make(smax)}
print("Density kg/m3 at 500 / 600 / 700 / 800 km:")
for n, f in ATM.items(): print(f"  {n:32s}", " ".join(f"{f(h*1e3):.2e}" for h in (500, 600, 700, 800)))
print(f"  {'MSIS, solar minimum':32s}", " ".join(f"{make(smin)(h*1e3):.2e}" for h in (500, 600, 700, 800)))
V = V_RECOIL
def fmt(t): return ">100,000" if t == float("inf") else (f"{t:,.0f}" if t >= 10 else f"{t:.1f}")
for name, f in ATM.items():
    do.rho_atm = f
    print(f"\n=== Atmosphere: {name} ===")
    print(" ring | case                      | perigee x apogee km | lifetime yr: full section / 0.8 m / 8 cm")
    for hk in (500, 550, 600, 650, 700, 750, 800):
        h = hk*1e3; s = ring_state(h); Brod = ballistic_coefficient(math.sqrt(s['m_struct']/RHO_CNT))
        cases = []
        o = released_cable(h); cases.append(("nominal release", o['h_p'], o['h_a']))
        o = released_cable(h, -V, V); cases.append(("slow half + radial recoil", o['h_p'], o['h_a']))
        for m, lab in ((1200, "same + 10% of casing swept"), (3000, "same + 25% of casing swept")):
            v = swept_cable_speed(h, m, -V); hp, ha, _, _ = orbit_from_state(R_EARTH+h, v, V); cases.append((lab, hp, ha))
        for lab, hp, ha in cases:
            L = [do.lifetime(max(hp, 1), ha, B) for B in (Brod, 600.0, 60.0)]
            print(f"  {hk} | {lab:26s}| {hp/1e3:5.0f} x {ha/1e3:6.0f}     | " + " / ".join(fmt(t) for t in L))
print("\nPrompt soot from collision ejecta (5% ejecta bound, 30% soot, 10% casing), Mt, beta = 1.0 / 1.5 / 2.0; threshold 50 Mt")
for hk in (250, 400, 500, 550, 600, 650, 700, 750, 800):
    print(f"  {hk} km: " + " / ".join(f"{ds.particulate_from_pop1(hk*1e3, beta=b, casing_frac=0.1)/1e9:.1f}" for b in (1.0, 1.5, 2.0)))
