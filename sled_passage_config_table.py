#!/usr/bin/env python3
"""
Prints the configuration table of the sled passage problem for review.
Book: The Orbital Ring, Volume IV: Not quoted by a chapter (a review tool).
Run: python sled_passage_config_table.py
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

Configuration table for the sled passage problem (2026-10-09, for Paul's review).

Uses the single-point casing model of sled_passage_limits.py (one meter of casing, no
stiffness, ideal feed-forward controller, pull-only EML), with three things made
variable that the limits script holds fixed:
  - the nominal gap G (the script's 100 mm),
  - the smallest gap allowed (the "no close calls" margin; the script's 20 mm),
  - where the spacecraft sits on the sled (front or rear).
The limiter and the EML are assumed to share the same gap, as in the other scripts.
Reaction plate: molybdenum.
"""
import math, sys
import numpy as np
import levitation_config as lc
import levitation_physics as lp
import sled_passage_simulation as sp

MU, W = sp.MU, sp.W_CASING
L_SLED, M_SLED_M = sp.L_SLED, sp.M_SLED_M
PLATE = "molybdenum"

def limiter(k_mult, wavelength, gap_nom):
    coils_per_m = 4.0 / wavelength
    K0 = (4.0 / lc.SAFETY_WAVELENGTH) * lc.SAFETY_TURNS_PER_COIL * lc.SAFETY_I_PERSISTENT
    K = k_mult * K0
    B = lp.halbach_surface_field(1, K / coils_per_m, lc.SAFETY_COIL_WIDTH, wavelength, coils_per_meter=coils_per_m)
    mat = lp.get_material(PLATE); sig = 1.0 / mat["resistivity"]
    a = (lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, sp.V_REL)
    F0 = lp.safety_repulsive_force(B, 0.0, wavelength, sig, *a, n_arrays=2)
    Fn = lp.safety_repulsive_force(B, gap_nom, wavelength, sig, *a, n_arrays=2)
    Qn = lp.eddy_current_heating_halbach(B, gap_nom, wavelength, sp.V_REL, mat, lc.SAFETY_PLATE_THICKNESS, lc.SAFETY_PLATE_WIDTH, n_arrays=2)
    return dict(F0=F0, k=2 * math.pi / wavelength, B=B, K=K, F_nom=Fn, Q_nom=Qn)

def response(M_sc, u, L_sc, lim, G, rear, F_max=sp.F_EML_MAX, f_hz=sp.PD_FREQ_HZ, t_settle=2.0, steps=12):
    """Smallest bottom gap and smallest top gap, arrays (n_mass, n_speed)."""
    M = np.asarray(M_sc, float)[:, None]; U = np.asarray(u, float)[None, :]
    ad = sp.a_down(U)
    q_zone = (M_SLED_M + M / L_sc) * ad; q_bare = M_SLED_M * ad
    t_zone = L_sc / U; t_sled = L_SLED / U
    dt = min(float(L_sc / np.max(U) / steps), 2e-3)
    n = int(math.ceil((float(np.max(t_sled)) + t_settle) / dt)) + 1
    om = 2 * math.pi * f_hz; k2, F0 = 2 * lim["k"], lim["F0"]
    w = np.zeros(np.broadcast(M, U).shape); v = np.zeros_like(w)
    gb = np.full_like(w, G); gt = np.full_like(w, G)
    for i in range(n):
        t = i * dt
        if rear:
            in_bare = t < (t_sled - t_zone); in_zone = (t >= t_sled - t_zone) & (t < t_sled)
        else:
            in_zone = t < t_zone; in_bare = (t >= t_zone) & (t < t_sled)
        q = np.where(in_zone, q_zone, 0.0) + np.where(in_bare, q_bare, 0.0)
        fk = -W - q
        Fe = np.clip(-fk - MU * (2 * om * v + om * om * w), 0.0, F_max)
        Flb = F0 * np.exp(-k2 * np.maximum(G - w, 0.0)); Flt = F0 * np.exp(-k2 * np.maximum(G + w, 0.0))
        v += (fk + Fe + Flt - Flb) / MU * dt; w += v * dt
        hb = w >= G; w[hb] = G; v[hb] = np.minimum(v[hb], 0.0)
        ht = w <= -G; w[ht] = -G; v[ht] = np.maximum(v[ht], 0.0)
        np.minimum(gb, G - w, out=gb); np.minimum(gt, G + w, out=gt)
    return gb, gt

MASSES = np.arange(250.0, 20001.0, 250.0) * 1000.0
SPEEDS = np.array([15e3, 20e3, 25e3, 30e3])
SCAN = np.arange(8e3, 30001.0, 250.0)

# (label, nominal gap m, smallest gap allowed m, K multiple, wavelength m, spacecraft at rear)
CONFIGS = [
    ("A  present design",                      0.100, 0.020, 1, 0.100, False),
    ("B  present, spacecraft at rear",         0.100, 0.020, 1, 0.100, True),
    ("C  K x3 limiter",                        0.100, 0.020, 3, 0.100, False),
    ("D  gap 150, margin 50",                  0.150, 0.050, 1, 0.100, False),
    ("E  gap 150, margin 50, K x3",            0.150, 0.050, 3, 0.100, False),
    ("F  gap 200, margin 100",                 0.200, 0.100, 1, 0.100, False),
    ("G  gap 200, margin 100, K x3",           0.200, 0.100, 3, 0.100, False),
    ("H  gap 200, margin 100, K x3, 200 mm",   0.200, 0.100, 3, 0.200, False),
    ("I  gap 200, margin 100, K x3, 300 mm",   0.200, 0.100, 3, 0.300, False),
    ("J  gap 200, margin 100, rear",           0.200, 0.100, 1, 0.100, True),
    ("K  gap 200, margin 100, K x3, rear",     0.200, 0.100, 3, 0.100, True),
    ("L  gap 300, margin 100, rear",           0.300, 0.100, 1, 0.100, True),
    ("M  gap 300, margin 100, K x3, 300 mm, rear", 0.300, 0.100, 3, 0.300, True),
    ("N  gap 300, margin 100, K x3, 300 mm",   0.300, 0.100, 3, 0.300, False),
]

def main():
    L_sc = 100.0
    print("SLED PASSAGE CONFIGURATION TABLE (single-point model, ideal controller, EML ceiling 871 kN/m, molybdenum plate)")
    print(f"Spacecraft zone {L_sc:.0f} m. V_REL {sp.V_REL:,.1f} m/s. Casing {MU:,.0f} kg/m, weight {W/1e3:.1f} kN/m.")
    print("Largest spacecraft mass (t) that keeps both gaps at or above the margin, by sled speed; top speed for 5,000 t;")
    print("limiter steady force and plate heating at the nominal gap; smallest bottom gap for 5,000 t at 30 km/s.\n")
    hdr = f"{'config':44s} {'gap':>4s} {'marg':>4s} {'K':>2s} {'wl':>4s} {'15km/s':>7s} {'20km/s':>7s} {'25km/s':>7s} {'30km/s':>7s} {'vmax5kt':>8s} {'Fnom':>8s} {'Qnom':>9s} {'gmin5kt30':>9s} {'Bsurf':>6s}"
    print(hdr); print("-" * len(hdr))
    for lab, G, marg, km, wl, rear in CONFIGS:
        lim = limiter(km, wl, G)
        gb, gt = response(MASSES, SPEEDS, L_sc, lim, G, rear)
        ok = (gb >= marg) & (gt >= marg)
        mm = []
        for j in range(len(SPEEDS)):
            idx = np.where(ok[:, j])[0]
            # largest mass below the first failure (envelope must be contiguous from light masses)
            bad = np.where(~ok[:, j])[0]
            top = (bad.min() - 1) if bad.size else len(MASSES) - 1
            mm.append("none" if top < 0 else (">20000" if not bad.size else f"{MASSES[top]/1e3:,.0f}"))
        g5, t5 = response([5.0e6], SCAN, L_sc, lim, G, rear)
        ok5 = (g5[0] >= marg) & (t5[0] >= marg)
        bad5 = np.where(~ok5)[0]
        vmax = SCAN[-1] if not bad5.size else (SCAN[bad5.min() - 1] if bad5.min() > 0 else float('nan'))
        print(f"{lab:44s} {G*1e3:4.0f} {marg*1e3:4.0f} {km:2d} {wl*1e3:4.0f} {mm[0]:>7s} {mm[1]:>7s} {mm[2]:>7s} {mm[3]:>7s} {vmax/1e3:8.2f} {lim['F_nom']/1e3:8.3f} {lim['Q_nom']:9.2f} {g5[0][-1]*1e3:9.1f} {lim['B']:6.2f}")
    print("\nUnits: gap, margin, wavelength, gmin in mm; masses in tonnes; vmax5kt in km/s (250 m/s steps); Fnom in kN/m; Qnom in W/m; Bsurf in T.")
    print("Free casing speed after a 5,000 t pass at 30 km/s: %.3f m/s. Travel for the casing weight alone to stop it: %.0f mm." % (
        (5e6 * -sp.a_down(30e3) / 30e3) / MU, 1e3 * ((5e6 * -sp.a_down(30e3) / 30e3) / MU) ** 2 / (2 * W / MU)))

if __name__ == "__main__":
    main()
