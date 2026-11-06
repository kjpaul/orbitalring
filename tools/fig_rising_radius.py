#!/usr/bin/env python3
"""
Figures of the rising-radius deployment (lim_simulation.py --rising-radius): the ring's altitude against time with the
hang limits, the tension and stress against time beside the fixed-radius run, the cable and casing speeds, the anchor
reaction of phase 2, and the energy budget as it accumulates.  Reads rising_radius_daily.csv, rising_radius_deployment_tension.csv
and rising_radius_summary.json from the output folder, and lim_deployment_tension.csv (the fixed-radius run) when it is there.
Book: Orbital Ring Engineering, Volume I: Chapter 3 (Section 3.2, Figure 3.2), Chapter 9, Chapter 11 (Section 11.1), Chapter 12 (Section 12.11).
Run: python tools/fig_rising_radius.py [--out=FOLDER]     (default <output folder>/figures_rising_radius; then
     python tools/fig03_02_deployment_tension.py --rising --out=<the same folder>   for the Figure 3.2 rerun)
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ring_altitude as ring

BLUE, RED, GREEN, GRAY, ORANGE = "#2b6ca3", "#b5451b", "#2e8b57", "gray", "#c98a1a"
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
FLOOR_KM = 150.0


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def save(fig, out, name):
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), **kw)
    plt.close(fig)
    print(f"saved {name}.svg and .png")


def main():
    src = ring.out_dir()
    out = next((a[6:] for a in sys.argv[1:] if a.startswith("--out=")), ring.out_dir("figures_rising_radius"))
    os.makedirs(out, exist_ok=True)
    D = read_csv(os.path.join(src, "rising_radius_daily.csv"))
    S = read_csv(os.path.join(src, "rising_radius_deployment_tension.csv"))
    J = json.load(open(os.path.join(src, "rising_radius_summary.json")))
    fixed_path = os.path.join(src, "lim_deployment_tension.csv")
    F = read_csv(fixed_path) if os.path.exists(fixed_path) else None
    t2 = J.get("t_phase2_days"); t_end = J["t_end_days"]
    day = D["day"]; h = D["altitude_km"]

    # 1. altitude against time, with what can hang
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.plot(day, h, color=BLUE, lw=2.4, label="Ring altitude (the casing)")
    ax.plot(day, h - 100.0, color=GREEN, lw=1.6, ls="--", label="100 km below the casing (the array support lines)")
    ax.plot(day, h - FLOOR_KM, color=ORANGE, lw=1.6, ls=":", label=f"Longest hang that stays above {FLOOR_KM:.0f} km")
    ax.axhline(FLOOR_KM, color=GRAY, lw=0.9, ls="--"); ax.text(3, FLOOR_KM + 6, f"{FLOOR_KM:.0f} km, the lower limit for what hangs (Section 11.1)", color=GRAY, fontsize=9)
    if t2:
        ax.axvline(t2, color=GRAY, lw=0.9, ls="--"); ax.text(t2 - 4, 40, f"casing at ground speed, day {t2:.0f}\n({J['altitude_phase2_km']:.0f} km): anchors attached", color=GRAY, fontsize=9, ha="right")
    ax.text(3, h[0] + 8, f"build orbit {h[0]:.0f} km", color=BLUE, fontsize=9); ax.text(t_end - 3, h[-1] + 8, f"{h[-1]:.0f} km at day {t_end:.0f}", color=BLUE, fontsize=9, ha="right")
    ax.set_xlim(0, t_end * 1.02); ax.set_ylim(0, 760); ax.set_xlabel("Time from the start of the deployment (days)"); ax.set_ylabel("Altitude (km)")
    ax.grid(True, alpha=0.3); ax.legend(loc="center right", fontsize=9)
    ax.set_title("The ring rises from the build orbit as the cable's strain builds", fontsize=12)
    fig.tight_layout(); save(fig, out, "rising_radius_altitude")

    # 2. tension and stress against time, both runs
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.plot(S["time_days"], S["tension_GN"], color=BLUE, lw=2.4, label="Cable tension, rising radius")
    if F is not None:
        ax.plot(F["time_days"], F["tension_GN"], color=BLUE, lw=1.4, ls="--", label="Cable tension, fixed 700 km radius")
    ax.set_xlabel("Time from the start of the deployment (days)"); ax.set_ylabel("Cable tension (GN)", color=BLUE)
    axb = ax.twinx(); axb.set_ylim(0, S["tension_GN"][-1] * 1.08 * S["stress_GPa"][-1] / S["tension_GN"][-1]); axb.set_ylabel("Cable stress (GPa), the same curves", color=RED)
    ax.set_ylim(0, S["tension_GN"][-1] * 1.08); ax.set_xlim(0, t_end * 1.02)
    if t2:
        ax.axvline(t2, color=GRAY, lw=0.9, ls="--"); ax.text(t2 - 4, S["tension_GN"][-1] * 0.05, f"casing at ground speed, day {t2:.0f}\n{J['tension_phase2_GN']:.0f} GN, phase 2 begins", color=GRAY, fontsize=9, ha="right")
    k = int(np.argmin(np.abs(S["v_casing_m_s"]))); ax.plot(S["time_days"][k], S["tension_GN"][k], "o", color=BLUE, ms=5)
    ax.text(S["time_days"][k] - 4, S["tension_GN"][k] + 25, f"casing at rest, day {S['time_days'][k]:.0f}, {S['tension_GN'][k]:.0f} GN", color=BLUE, fontsize=9, ha="right")
    ax.grid(True, alpha=0.3); ax.legend(loc="upper left", fontsize=9)
    ax.set_title("Cable tension through the deployment: the rising radius against the fixed radius", fontsize=12)
    fig.tight_layout(); save(fig, out, "rising_radius_tension")

    # 3. speeds
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.plot(S["time_days"], S["v_cable_m_s"], color=RED, lw=2.2, label="Cable speed westward, rising radius")
    ax.plot(S["time_days"], S["v_casing_m_s"], color=GREEN, lw=2.2, label="Casing speed westward (negative: eastward), rising radius")
    if F is not None:
        ax.plot(F["time_days"], F["v_cable_m_s"], color=RED, lw=1.2, ls="--", label="Cable speed, fixed radius")
        ax.plot(F["time_days"], F["v_casing_m_s"], color=GREEN, lw=1.2, ls="--", label="Casing speed, fixed radius")
    ax.axhline(0, color=GRAY, lw=0.8)
    if t2:
        ax.axvline(t2, color=GRAY, lw=0.9, ls="--")
    ax.set_xlim(0, t_end * 1.02); ax.set_xlabel("Time from the start of the deployment (days)"); ax.set_ylabel("Speed in the launch frame (m/s)")
    ax.grid(True, alpha=0.3); ax.legend(loc="center right", fontsize=9)
    ax.set_title("Cable and casing speeds: the build orbit's 7,598 m/s against 7,504 m/s at 700 km", fontsize=12)
    fig.tight_layout(); save(fig, out, "rising_radius_speeds")

    # 4. phase 2: the anchor reaction
    if t2:
        m = D["phase"] >= 2
        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.plot(day[m], D["anchor_force_total_MN"][m], color=RED, lw=2.2, label="Anchor reaction, total around the ring (the Earth's torque / r)")
        mt = m & (D["lim_thrust_total_MN"] > 0)             # the closing row of a run made before 2026-11-06's F_last fix carries no thrust
        ax.plot(day[mt], D["lim_thrust_total_MN"][mt], color=BLUE, lw=1.4, ls="--", label="LIM thrust, total around the ring")
        ax.set_xlabel("Time from the start of the deployment (days)"); ax.set_ylabel("Force (MN)")
        axb = ax.twinx(); axb.set_ylabel(f"Per station, {J['n_stations']} stations (kN)")
        lo, hi = ax.get_ylim(); ax.set_ylim(0, hi * 1.15); axb.set_ylim(0, hi * 1.15 / J["n_stations"] * 1e3)
        ax.text(day[m][0] + 2, hi * 1.05, f"phase 2: {J['phase2_duration_days']:.0f} days, mean {J['anchor_mean_MN']:.1f} MN total, {J['anchor_mean_MN'] / J['n_stations'] * 1e3:.1f} kN per station;"
                f"\nthe Earth supplies {J['L_earth']:.3g} kg m$^2$/s, {J['L_earth_share_of_cable'] * 100:.2f} % of the cable's angular momentum", fontsize=9, color=GRAY)
        ax.grid(True, alpha=0.3); ax.legend(loc="lower left", fontsize=9)
        ax.set_title("Phase 2: the casing held by the anchor lines while the LIMs finish the cable", fontsize=12)
        fig.tight_layout(); save(fig, out, "rising_radius_anchor_phase2")

    # 5. energy
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.plot(day, D["lim_work_EJ"], color="black", lw=2.4, label="LIM work (thrust x relative speed)")
    ax.plot(day, D["KE_cable_EJ"] - D["KE_cable_EJ"][0], color=RED, lw=1.8, label="Cable kinetic energy change")
    ax.plot(day, D["KE_casing_EJ"] - D["KE_casing_EJ"][0], color=GREEN, lw=1.8, label="Casing kinetic energy change")
    ax.plot(day, D["PE_EJ"] - D["PE_EJ"][0], color=BLUE, lw=1.8, label="Potential energy change (the ring lifted)")
    ax.plot(day, D["elastic_EJ"], color=ORANGE, lw=1.8, label="Elastic energy of the cable")
    ax.plot(day, D["earth_work_EJ"], color=GRAY, lw=1.8, label="The Earth's work through the anchors")
    ax.axhline(ring.E_DEPLOY / 1e18, color="black", lw=0.9, ls=":"); ax.text(3, ring.E_DEPLOY / 1e18 + 0.4, f"fixed-radius deployment energy {ring.E_DEPLOY / 1e18:.1f} EJ", fontsize=9)
    ax.axhline(0, color=GRAY, lw=0.8)
    if t2:
        ax.axvline(t2, color=GRAY, lw=0.9, ls="--")
    ax.set_xlim(0, t_end * 1.02); ax.set_xlabel("Time from the start of the deployment (days)"); ax.set_ylabel("Energy (EJ)")
    ax.grid(True, alpha=0.3); ax.legend(loc="lower right", fontsize=9)
    ax.set_title(f"Energy budget of the rising-radius deployment: {J['E_lim_loop_EJ']:.1f} EJ of LIM work", fontsize=12)
    fig.tight_layout(); save(fig, out, "rising_radius_energy")
    print(f"figures in {out}")


if __name__ == "__main__":
    main()
