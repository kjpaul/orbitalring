"""
Solar particle event dose on the ring at 700 km (Chapter 5, Section 5.9).
Book: Orbital Ring Engineering, Volume I: Chapter 5, Section 5.9.
Run from the radiation folder: python spe_700.py [--log=FILE]
Environment: none; reads the Geant4 proton runs in out_gcr/ and out_trap*/ (through gcr.py and fold_trapped.py), the
PSTAR table, and sweetspot.py for the room geometry; writes the log to FILE (default spe_700.txt here).

What is computed
1. Reference events. The event-integrated proton fluence spectra of
   (a) the August 1972 event in King's form, J(>E) = 7.9e9 exp((30 - E)/26.5) protons/cm2, E in MeV (King 1974,
       J. Spacecraft Rockets 11, 401, as given in the SPENVIS help on solar proton models), the classic design event:
       the largest fluence of the space age below 100 MeV and a soft spectrum; and
   (b) the ground level enhancements with the hardest spectra, in the modified Band function of Koldobskiy, Raukunen,
       Vainio, Kovaltsov and Usoskin 2021 (A&A 647, A132, doi 10.1051/0004-6361/202040058, Table 2): F(>R) = J1 R^-g1
       exp(-R/R1) below Rb, J2 R^-g2 exp(-R/R2) above, R in GV, omnidirectional fluence per cm2, fitted to spacecraft
       data above 30 MeV and to the neutron monitor network above 1 GV. The events used: GLE 5 (23 Feb 1956, the hardest
       spectrum on record), GLE 42 (29 Sep 1989), GLE 43 to 45 (October 1989), GLE 69 (20 Jan 2005), GLE 24 (4 Aug 1972,
       the same event as (a) in this fit), GLE 59 (14 Jul 2000) and GLE 65 (28 Oct 2003).
2. Quiet-time cutoff: Stormer's formula for the field of 2026 in every direction of the sky (sweetspot.stormer, the
   same function the cosmic ray dose uses: 11.6 GV from overhead, 7.9 GV from the western horizon, 46 GV from the east).
   The share of each event's fluence above the cutoff.
3. Storm-time cutoff: the cutoff drop found for the equator in the 20 November 2003 superstorm (Dst -422 nT), 0.5 to
   0.6 GV (Danilova, Demina, Ptitsyna and Tyasto 2019, Geomagn. Aeron. 59, 147, doi 10.1134/S0016793219020051,
   computed on a 5 x 15 degree world grid in a model magnetosphere at the storm maximum; the largest drop anywhere was
   1.8 GV, at middle and low latitudes). The drop is applied as the same subtraction in every direction of arrival,
   which is an assumption: the paper gives the vertical cutoff. The 1.8 GV drop is run as a bound.
4. Dose per event: the part of the spectrum above the cutoff of each direction, folded with the Geant4 depth-dose
   response of the cosmic ray proton runs (out_gcr/, 3 GeV to 1 TeV, the response gcr.py uses), direction by direction
   over the sky as sweetspot_run.py does, with the Earth blocking the directions below its limb, for a room with walls
   equal to 5, 20 and 100 cm of water all round and for the room under the cable (sweetspot.thickness, h = 1 m, thin
   walls and floor). The event fluence is taken as isotropic outside the magnetosphere, which is right for the
   event-integrated fluence (the anisotropic phase of a GLE lasts an hour or so).
   For scale, the same events in free space (no cutoff, no Earth) behind the three walls, with the trapped-proton
   response (fold_trapped, 12 MeV to 1.2 GeV) for the part of the spectrum below 2 GeV and the cosmic ray response above.
Every number printed is from this run or from the cited documents.
"""
import os, sys
import numpy as np
from sweetspot import U, thickness, earth_blocked, stormer
import gcr
import fold_trapped as ft

MP = 0.938272                      # GeV, the proton rest energy
GY = 1.602e-10                     # Gy per MeV/g
dOm = 4 * np.pi / len(U)
ALT = 700.0
NASA_SPE_LIMIT_MSV = 250.0         # NASA-STD-3001, the chapter's reference [23]

# ---- the reference event spectra
KING = dict(J0=7.9e9, E0=26.5, Eref=30.0)                      # King 1974, August 1972: J(>E) = J0 exp((Eref - E)/E0), E in MeV
MBF = {                                                        # Koldobskiy et al. 2021, Table 2: g1, R1 (GV), J2 (cm-2), g2, R2 (GV), Rb (GV)
    "GLE 5, 23 Feb 1956":  (1.59, 0.770, 1.63e8, 4.84, 8.614, 2.748),
    "GLE 24, 4 Aug 1972":  (2.02, 0.152, 4.35e6, 11.96, np.inf, 1.511),
    "GLE 42, 29 Sep 1989": (-1.95, 0.069, 9.89e6, 3.39, 5.935, 0.373),
    "GLE 43, 19 Oct 1989": (-1.21, 0.109, 6.59e6, 4.44, 2.960, 0.639),
    "GLE 44, 22 Oct 1989": (0.35, 0.122, 1.84e6, 6.07, 5.807, 0.713),
    "GLE 45, 24 Oct 1989": (-1.77, 0.077, 1.19e7, 2.95, 1.131, 0.390),
    "GLE 59, 14 Jul 2000": (2.80, 0.276, 3.46e6, 5.85, 2.191, 0.963),
    "GLE 65, 28 Oct 2003": (6.49, -0.525, 5.31e6, 0.01, 0.447, 1.565),
    "GLE 69, 20 Jan 2005": (1.16, 0.295, 8.03e6, 4.73, 2.392, 1.201),
}


def R_of_E(E_GeV):
    return np.sqrt(E_GeV * (E_GeV + 2 * MP))


def E_of_R(R_GV):
    return np.sqrt(R_GV**2 + MP**2) - MP


def king_F(E_MeV):
    """integral fluence above E (MeV), protons/cm2"""
    return KING["J0"] * np.exp((KING["Eref"] - np.asarray(E_MeV, float)) / KING["E0"])


def mbf_F(R, p):
    """integral fluence above rigidity R (GV), protons/cm2, modified Band function (Koldobskiy et al. 2021, Eqs 1 to 3)"""
    g1, R1, J2, g2, R2, Rb = p
    R = np.asarray(R, float)
    g0 = g2 - g1
    J1 = J2 * Rb**(-g0) * np.exp(g0)
    lo = J1 * R**(-g1) * np.exp(-R / R1)
    hi = J2 * R**(-g2) * (np.exp(-R / R2) if np.isfinite(R2) else 1.0)
    return np.where(R < Rb, lo, hi)


def F_of_E(name, E_MeV):
    """integral fluence above kinetic energy E (MeV) for the named event"""
    if name == "King 1972":
        return king_F(E_MeV)
    return mbf_F(R_of_E(np.asarray(E_MeV, float) / 1e3), MBF[name])


EVENTS = ["King 1972"] + list(MBF)

# ---- the proton depth-dose responses
_En_g, _R_g = {}, {}
for key in ("H", "D"):
    En, R = gcr.table("H", key)                # En in GeV (kinetic), R[energy, depth bin of 2 g/cm2], MeV/g per proton/cm2
    _En_g[key], _R_g[key] = En, R


def resp_gcr(E_GeV, t, key="H"):
    """response at energy E (GeV, array) and slant depth t (g/cm2, array, same length), from the cosmic ray proton runs:
    log-log interpolation in energy (power-law continuation outside 3 GeV to 1 TeV), linear in depth, exponential
    continuation beyond 590 g/cm2 with the 400 to 590 slope (as sweetspot_run.gcr_at)"""
    En, R = _En_g[key], _R_g[key]
    logR = np.log(np.maximum(R, 1e-12)); le = np.log(En); x = np.log(np.clip(E_GeV, 1e-3, 1e5))
    i = np.clip(np.searchsorted(le, x) - 1, 0, len(le) - 2); f = (x - le[i]) / (le[i + 1] - le[i])
    tt = np.clip(t, 1.0, 590.0); j = np.clip((tt / 2.0 - 0.5).astype(int), 0, 298); g = tt / 2.0 - 0.5 - j
    def at(jj): return np.exp(logR[i, jj] * (1 - f) + logR[i + 1, jj] * f)
    val = at(j) * (1 - g) + at(j + 1) * g
    a = at(np.full_like(j, 200)); b = at(np.full_like(j, 295))
    lam = np.maximum(190.0 / np.log(np.maximum(a / b, 1.0001)), 50.0)
    return np.where(t > 590, val * np.exp(-(t - 590) / lam), val)


def dose_above_cutoff(name, rc, t, key="H", n_e=160):
    """dose (Sv or Gy) at the point from the part of the event fluence above the cutoff rc (GV, per direction) through
    slant depth t (g/cm2, per direction), isotropic fluence, summed over the open directions; rc may be inf (blocked)"""
    out = np.zeros(len(U))
    ok = np.isfinite(rc)
    if not ok.any(): return 0.0
    # energy grid from the lowest cutoff to 100 GV, shared; each direction integrates from its own cutoff
    Emin = E_of_R(np.min(rc[ok])); Emax = E_of_R(100.0)
    Eg = np.geomspace(max(Emin, 1e-3), Emax, n_e); Em = np.sqrt(Eg[1:] * Eg[:-1])
    F = F_of_E(name, Eg * 1e3); n = -np.diff(F); n[n < 0] = 0        # protons/cm2 in each energy bin
    Ecut = E_of_R(rc[ok]); tt = t[ok]
    for b in range(len(Em)):
        if n[b] <= 0: continue
        w = np.where(Eg[b + 1] <= Ecut, 0.0, np.where(Eg[b] >= Ecut, 1.0, (Eg[b + 1] - Ecut) / (Eg[b + 1] - Eg[b])))   # share of the bin above the cutoff
        out[ok] += n[b] * w * resp_gcr(np.full(len(tt), Em[b]), tt, key)
    return float((out * dOm / (4 * np.pi)).sum() * GY)


def dose_no_cutoff(name, t, key="H"):
    """bounding case, the field switched off: the whole isotropic fluence through the open sky (the Earth still blocks),
    direction by direction. Below 2 GeV the trapped-proton response of fold_trapped (12 MeV to 2 GeV on its EF grid,
    interpolated in depth over its T grid in the logarithm, the 300 to 500 g/cm2 slope continued beyond 500 as
    sweetspot_trapped does); above 2 GeV the cosmic ray proton response."""
    eb = earth_blocked(ALT); ok = ~eb; tt = t[ok]
    EF = ft.EF; F = F_of_E(name, EF); n = -np.diff(F); n[n < 0] = 0
    T = ft.T; lt = np.log(np.maximum(ft.RESP[key], 1e-30))             # [T, EF]
    r_mid = 0.5 * (lt[:, 1:] + lt[:, :-1])
    sl = np.minimum((r_mid[-1] - r_mid[-3]) / (T[-1] - T[-3]), 0.0)         # the response never grows with depth past 5 m
    acc = np.zeros(len(tt))
    for b in range(len(n)):
        if n[b] <= 0 or r_mid[:, b].max() < np.log(1e-20): continue
        v = np.interp(np.clip(tt, T[0], T[-1]), T, r_mid[:, b])
        v = np.where(tt > T[-1], r_mid[-1, b] + sl[b] * (tt - T[-1]), v)
        v = np.where((tt > T[-1]) & (r_mid[-1, b] < np.log(1e-20)), -np.inf, v)   # nothing reaches 5 m at this energy: nothing deeper either
        acc += n[b] * np.exp(v)
    low = float((acc * dOm / (4 * np.pi)).sum() * GY)
    rc = np.where(eb, np.inf, R_of_E(2.0))                              # the part above 2 GeV, no cutoff
    high = dose_above_cutoff(name, rc, t, key)
    return low + high


def dose_free_space(name, t_uniform, key="H"):
    """dose at the centre of a water sphere of radius t under the whole isotropic event fluence, no cutoff, no Earth:
    trapped-proton response (12 MeV to 2 GeV grid of fold_trapped) plus the cosmic ray response above 2 GeV"""
    EF = ft.EF; F = F_of_E(name, EF); n = -np.diff(F); n[n < 0] = 0
    k = int(np.where(ft.T == t_uniform)[0][0])
    r = 0.5 * (ft.RESP[key][k, 1:] + ft.RESP[key][k, :-1])
    low = float((r * n).sum() * GY)
    Eg = np.geomspace(2.0, E_of_R(100.0), 200); Em = np.sqrt(Eg[1:] * Eg[:-1])
    Fh = F_of_E(name, Eg * 1e3); nh = -np.diff(Fh); nh[nh < 0] = 0
    high = float((nh * resp_gcr(Em, np.full(len(Em), float(t_uniform)), key)).sum() * GY)
    return low + high, low, high


def main():
    log = next((a[6:] for a in sys.argv[1:] if a.startswith("--log=")), "spe_700.txt")
    lines = []
    def P(s=""):
        print(s); lines.append(s)
    rc0 = stormer(ALT)                                  # quiet-time cutoff by direction, GV
    eb = earth_blocked(ALT)
    uz = U[:, 2]; ux = U[:, 0]
    k_up = int(np.argmax(uz)); k_w = int(np.argmin(ux + 50 * np.abs(uz)))
    # the three reference directions by the formula itself (the grid has no point exactly on the horizon)
    r_e = (6371.2 + ALT) / 6371.2
    rc_up_f, rc_w_f, rc_e_f = (57.0 / (r_e * r_e * (1 + np.sqrt(1 - cg))**2) for cg in (0.0, -1.0, 1.0))
    P("SOLAR PARTICLE EVENT DOSE ON THE RING AT 700 KM (Chapter 5, Section 5.9)")
    P("Inputs: King 1974 August 1972 spectrum (SPENVIS form), Koldobskiy et al. 2021 Table 2 modified Band fits for the GLEs,")
    P("Stormer cutoff of sweetspot.py (field of 2026), storm-time drop from Danilova et al. 2019 (20 Nov 2003, Dst -422 nT).")
    P()
    P("1. Reference events: event-integrated omnidirectional fluence, protons/cm2 above the energy")
    P(f"   {'event':24s} | {'>10 MeV':>9s} | {'>30 MeV':>9s} | {'>100 MeV':>9s} | {'>1 GeV':>9s} | {'>7.0 GeV (7.9 GV)':>17s} | {'>10.7 GeV (11.6 GV)':>19s}")
    for name in EVENTS:
        f = [F_of_E(name, e) for e in (10, 30, 100, 1000, E_of_R(7.9) * 1e3, E_of_R(11.6) * 1e3)]
        P(f"   {name:24s} | {f[0]:9.2e} | {f[1]:9.2e} | {f[2]:9.2e} | {f[3]:9.2e} | {f[4]:17.2e} | {f[5]:19.2e}")
    P("   (the Band fits are not meant below 30 MeV; the >10 MeV column of those rows is the fit extrapolated, for the record only)")
    P("   King 1972 above 7.0 GeV: exp((30 - 7000)/26.5) = exp(-263) of the >30 MeV fluence: zero to every purpose.")
    P()
    P("2. Quiet-time cutoff, Stormer formula with the 2026 dipole (sweetspot.stormer, K = 57.0 GV):")
    P(f"   from overhead {rc_up_f:.2f} GV ({E_of_R(rc_up_f):.2f} GeV), from the western horizon {rc_w_f:.2f} GV ({E_of_R(rc_w_f):.2f} GeV), "
      f"from the eastern horizon {rc_e_f:.1f} GV; the chapter gives 11.6, 7.9 and 46 GV (Equation 5.1 uses C = 56.9 GV and R_E = 6,378 km, 0.2 percent apart).")
    P(f"   On the direction grid used below: overhead {rc0[k_up]:.2f} GV, the most westerly grid direction {rc0[k_w]:.2f} GV.")
    P(f"   Earth blocks {100 * eb.mean():.1f} percent of the sky; the lowest cutoff over the open sky is {rc0[~eb].min():.2f} GV.")
    P()

    def fractions(dR):
        rc = np.where(eb, np.inf, np.maximum(rc0 - dR, 0.0))
        rows = []
        for name in EVENTS:
            F30 = F_of_E(name, 30.0)
            f_up = F_of_E(name, E_of_R(max(rc0[k_up] - dR, 0)) * 1e3) / F30
            f_w = F_of_E(name, E_of_R(max(rc0[k_w] - dR, 0)) * 1e3) / F30
            # the share of the isotropic fluence that reaches the point through the open sky, all directions
            Fdir = np.where(np.isfinite(rc), F_of_E(name, E_of_R(np.where(np.isfinite(rc), rc, 1.0)) * 1e3), 0.0)
            f_sky = float((Fdir * dOm).sum() / (4 * np.pi)) / F30
            rows.append((name, f_up, f_w, f_sky))
        return rc, rows

    cases = [("quiet", 0.0), ("storm, equatorial drop 0.6 GV (20 Nov 2003 class, Dst -422 nT)", 0.6),
             ("bound, the largest drop seen anywhere in that storm, 1.8 GV", 1.8)]
    P("3. Share of each event's fluence above 30 MeV that passes the cutoff (quiet, and with the storm-time drop):")
    for label, dR in cases:
        rc, rows = fractions(dR)
        P(f"   {label}: cutoff overhead {max(rc0[k_up] - dR, 0):.2f} GV, western horizon {max(rc0[k_w] - dR, 0):.2f} GV")
        P(f"     {'event':24s} | {'overhead':>10s} | {'from the west':>13s} | {'whole open sky':>14s}")
        for name, f_up, f_w, f_sky in rows:
            P(f"     {name:24s} | {f_up:10.2e} | {f_w:13.2e} | {f_sky:14.2e}")
    P()

    P("4. Dose per event, mSv (dose equivalent, ICRP 60 quality factors as in the chapter), at the point in the room:")
    rooms = [("walls 5 cm water all round", lambda: (np.full(len(U), 5.0), None)),
             ("walls 20 cm water all round", lambda: (np.full(len(U), 20.0), None)),
             ("walls 100 cm water all round", lambda: (np.full(len(U), 100.0), None)),
             ("under the cable, h = 1 m, thin walls and floor (5 cm)", lambda: thickness(1.0, t_wall=5.0, t_floor=5.0))]
    for label, dR in cases + [("no cutoff at all (bounding case, the field switched off; the Earth still blocks 28 percent of the sky)", None)]:
        rc = None if dR is None else np.where(eb, np.inf, np.maximum(rc0 - dR, 0.0))
        P(f"   {label}")
        P(f"     {'event':24s} | " + " | ".join(f"{r[0][:28]:>28s}" for r in rooms))
        for name in EVENTS:
            vals = []
            for _, geom in rooms:
                t, _ = geom()
                vals.append((dose_no_cutoff(name, t) if rc is None else dose_above_cutoff(name, rc, t)) * 1e3)
            P(f"     {name:24s} | " + " | ".join(f"{v:28.3g}" for v in vals))
    P(f"   NASA's design limit for one solar particle event: {NASA_SPE_LIMIT_MSV:.0f} mSv (NASA-STD-3001, reference [23] of the chapter).")
    P()
    P("   For scale, the same events in free space (no magnetic cutoff, no Earth), dose at the centre of a water sphere, mSv:")
    P(f"     {'event':24s} | {'5 cm':>10s} | {'20 cm':>10s} | {'100 cm':>10s}   (of which above 2 GeV, 100 cm)")
    for name in EVENTS:
        d = [dose_free_space(name, t) for t in (5.0, 20.0, 100.0)]
        P(f"     {name:24s} | {d[0][0] * 1e3:10.3g} | {d[1][0] * 1e3:10.3g} | {d[2][0] * 1e3:10.3g}   ({d[2][2] * 1e3:.3g})")
    P("   (the free-space figures use the trapped-proton response below 2 GeV, which has no entry below 12 MeV; a response")
    P("    that stops at 2 GeV for the King spectrum changes nothing, its fluence above 2 GeV is exp(-74) of the total)")
    P()
    P("5. Under the cable, the directions that pass the storm-time cutoff: share of the open sky and the cable's cover")
    rc = np.where(eb, np.inf, np.maximum(rc0 - 0.6, 0.0)); t, hit = thickness(1.0, t_wall=5.0, t_floor=5.0)
    open_ = ~eb
    P(f"   open sky {100 * open_.mean():.1f} percent of all directions; of it, through the cable {100 * (open_ & hit).mean():.1f} percent, past the cable {100 * (open_ & ~hit).mean():.1f} percent")
    P(f"   the lowest storm-time cutoff over the open sky {rc[open_].min():.2f} GV, the mean over the open sky {rc[open_].mean():.2f} GV")
    with open(log, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nlog written to {log}")


if __name__ == "__main__":
    main()
