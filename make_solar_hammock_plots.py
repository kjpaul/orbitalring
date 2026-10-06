#!/usr/bin/env python3
"""
Regenerate the solar hammock 24-hour power profile plots at 250, 750, and
1000 km altitudes using the restored Chapter 7 physics from the orbital
ring power module.

The "solar hammock" is a flat horizontal bifacial panel platform at a
given altitude above the equator. As Earth rotates (or equivalently, as
the hammock orbits), the sun's position relative to the platform cycles
through all angles once per rotation. The 24-hour average is therefore
the same as the phi-averaged electrical output from power_physics.

Altitude dependence comes from two competing effects:
  1. The Earth view factor (R_E/r)^2 shrinks with altitude, which
     reduces the peak albedo.
  2. The shadow half-angle grows with altitude (the Earth blocks less
     of the sky), which extends the sunlit arc.

The two effects partially cancel, leaving a weak altitude dependence.
Higher altitudes gain a longer sunlit period but lose some albedo flux.

Outputs:
  solar_hammock_250km_power.png    (orbital ring reference)
  solar_hammock_750km_power.png
  solar_hammock_1000km_power.png
  solar_hammock_24hr_comparison.png
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import power_config as cfg
import power_physics as pp


def ring_physics_at_altitude(h_m):
    """Compute albedo peak, view factor, and shadow half-angle for a
    flat bifacial panel at altitude h_m above Earth.

    Returns a dict with:
        r           orbital radius (m)
        F_view      Earth view factor (R_E/r)^2
        alb_peak    peak Earth-albedo flux on back face (W/m^2)
        shadow_ha   shadow half-angle measured from solar noon (rad)
    """
    r = cfg.R_EARTH + h_m
    F_view = (cfg.R_EARTH / r) ** 2
    alb_peak = cfg.EARTH_ALBEDO * cfg.SOLAR_CONSTANT * F_view
    shadow_ha = math.pi - math.asin(cfg.R_EARTH / r)
    return {
        "h_m": h_m,
        "r": r,
        "F_view": F_view,
        "alb_peak": alb_peak,
        "shadow_ha": shadow_ha,
    }


def direct_flux(phi):
    """S0 * max(0, cos(phi)). Front-face direct solar flux."""
    c = math.cos(phi)
    return cfg.SOLAR_CONSTANT * c if c > 0 else 0.0


def albedo_flux(phi, alb_peak):
    """A * S0 * F * max(0, (1+cos(phi))/2). Back-face Earth albedo."""
    phase = 0.5 * (1.0 + math.cos(phi))
    return alb_peak * phase if phase > 0 else 0.0


def hammock_profile(h_km, n=2000):
    """Compute a 24-hour power profile at altitude h_km (kilometers)."""
    phys = ring_physics_at_altitude(h_km * 1000.0)

    phi = np.linspace(-math.pi, math.pi, n)
    eff = cfg.CELL_EFFICIENCY * cfg.AR_BOOST * cfg.PANEL_PACKING

    p_direct = np.array([direct_flux(p) for p in phi]) * eff
    p_albedo = np.array([albedo_flux(p, phys["alb_peak"]) for p in phi]) * eff
    p_total = p_direct + p_albedo

    # Fine-grid 24-hour average (independent of plot resolution).
    n_fine = 200_000
    phi_fine = np.linspace(-math.pi, math.pi, n_fine, endpoint=False)
    p_fine = np.array([
        (direct_flux(p) + albedo_flux(p, phys["alb_peak"])) * eff
        for p in phi_fine
    ])
    p_avg = float(np.mean(p_fine))
    p_peak = float(np.max(p_total))

    # Local solar time from phi. phi = 0 -> noon (12 h); phi = -pi -> midnight (0 h).
    t_hours = 12.0 + phi * 12.0 / math.pi

    # Shadow boundary in hours: sunlit when |phi| < shadow_ha.
    t_shadow_morning = 12.0 - phys["shadow_ha"] * 12.0 / math.pi
    t_shadow_evening = 12.0 + phys["shadow_ha"] * 12.0 / math.pi

    return {
        "h_km": h_km,
        "phys": phys,
        "t_hours": t_hours,
        "p_direct": p_direct,
        "p_albedo": p_albedo,
        "p_total": p_total,
        "p_avg": p_avg,
        "p_peak": p_peak,
        "t_shadow_morning": t_shadow_morning,
        "t_shadow_evening": t_shadow_evening,
    }


def plot_single_altitude(profile, out_path):
    """Draw a single-altitude 24-hour power chart in the Ch7 style."""
    h_km = profile["h_km"]
    t = profile["t_hours"]
    p_direct = profile["p_direct"]
    p_total = profile["p_total"]

    fig, ax = plt.subplots(figsize=(10, 5))

    direct_color = "#f4d03f"   # warm yellow
    albedo_color = "#5dade2"   # ocean blue
    shadow_color = "#d0d0d0"

    # Shadow background
    ax.axvspan(0.0, profile["t_shadow_morning"],
               color=shadow_color, alpha=0.55, zorder=1)
    ax.axvspan(profile["t_shadow_evening"], 24.0,
               color=shadow_color, alpha=0.55, zorder=1)

    # Stacked fills
    ax.fill_between(t, 0.0, p_direct,
                    facecolor=direct_color, edgecolor="#c8a93a",
                    linewidth=0.7, label="Direct sunlight (front surface)",
                    zorder=3)
    ax.fill_between(t, p_direct, p_total,
                    facecolor=albedo_color, edgecolor="#3f7ca8",
                    linewidth=0.7, label="Earth albedo (back surface)",
                    zorder=3)

    ax.plot(t, p_total, color="#1a2e4a", linewidth=1.1,
            label="Total electrical output", zorder=4)

    ax.axhline(profile["p_avg"], color="#c0392b", linestyle="--",
               linewidth=1.6,
               label=f"24-hour average: {profile['p_avg']:.0f} W/m\u00b2",
               zorder=5)

    # Shadow labels
    ax.text(profile["t_shadow_morning"] / 2.0, 150.0, "Earth's\nShadow",
            ha="center", va="center", color="#555555",
            fontsize=10, fontstyle="italic", zorder=6)
    ax.text((profile["t_shadow_evening"] + 24.0) / 2.0, 150.0,
            "Earth's\nShadow",
            ha="center", va="center", color="#555555",
            fontsize=10, fontstyle="italic", zorder=6)

    ax.set_xlim(0.0, 24.0)
    ax.set_ylim(0.0, 900.0)
    ax.set_xticks(np.arange(0, 25, 3))
    ax.set_xticklabels([
        "00:00\n(midnight)", "03:00", "06:00", "09:00", "12:00\n(noon)",
        "15:00", "18:00", "21:00", "24:00\n(midnight)",
    ])
    ax.set_yticks(np.arange(0, 901, 100))
    ax.set_xlabel("Local Solar Time (hours)")
    ax.set_ylabel("Electrical Power Output (W/m\u00b2)")
    ax.set_title(
        f"24-Hour Power Output of 1 m\u00b2 Bifacial Solar Panel "
        f"on Solar Hammock\n"
        f"({h_km} km altitude, 45% efficient cells with nanostructured AR "
        f"coating)",
        fontsize=11,
    )
    ax.grid(True, which="major", linestyle="-", linewidth=0.4,
            color="#dddddd", zorder=0)
    ax.legend(loc="upper right", framealpha=0.95, fontsize=9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  saved: {out_path}")
    plt.close(fig)


def plot_comparison(profiles, out_path):
    """Three-panel comparison plot for 250, 750, and 1000 km."""
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

    direct_color = "#f4d03f"
    albedo_color = "#5dade2"
    shadow_color = "#d0d0d0"

    for ax, profile in zip(axes, profiles):
        t = profile["t_hours"]
        p_direct = profile["p_direct"]
        p_total = profile["p_total"]

        ax.axvspan(0.0, profile["t_shadow_morning"],
                   color=shadow_color, alpha=0.55, zorder=1)
        ax.axvspan(profile["t_shadow_evening"], 24.0,
                   color=shadow_color, alpha=0.55, zorder=1)

        ax.fill_between(t, 0.0, p_direct,
                        facecolor=direct_color, edgecolor="#c8a93a",
                        linewidth=0.7,
                        label="Direct sunlight (front surface)", zorder=3)
        ax.fill_between(t, p_direct, p_total,
                        facecolor=albedo_color, edgecolor="#3f7ca8",
                        linewidth=0.7,
                        label="Earth albedo (back surface)", zorder=3)
        ax.plot(t, p_total, color="#1a2e4a", linewidth=1.1, zorder=4)

        ax.axhline(profile["p_avg"], color="#c0392b", linestyle="--",
                   linewidth=1.5,
                   label=f"24-hr avg: {profile['p_avg']:.0f} W/m\u00b2",
                   zorder=5)

        ax.set_title(
            f"{profile['h_km']} km altitude   "
            f"(peak {profile['p_peak']:.0f} W/m\u00b2, "
            f"avg {profile['p_avg']:.0f} W/m\u00b2, "
            f"view factor {profile['phys']['F_view']:.3f}, "
            f"shadow \u00b1{math.degrees(profile['phys']['shadow_ha']):.1f}\u00b0)",
            fontsize=11,
        )
        ax.set_ylabel("Electrical power (W/m\u00b2)")
        ax.set_xlim(0.0, 24.0)
        ax.set_ylim(0.0, 900.0)
        ax.set_yticks(np.arange(0, 901, 200))
        ax.grid(True, linestyle="-", linewidth=0.4, color="#dddddd",
                zorder=0)
        if profile is profiles[0]:
            ax.legend(loc="upper right", framealpha=0.95, fontsize=9)

    axes[-1].set_xlabel("Local Solar Time (hours)")
    axes[-1].set_xticks(np.arange(0, 25, 3))
    axes[-1].set_xticklabels([
        "00:00", "03:00", "06:00", "09:00", "12:00",
        "15:00", "18:00", "21:00", "24:00",
    ])

    fig.suptitle(
        "Solar Hammock 24-hour Power Profile at Three Altitudes\n"
        "(flat horizontal bifacial panel, 45% cells with AR coating, "
        "A = 0.30)",
        fontsize=12, y=0.995,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  saved: {out_path}")
    plt.close(fig)


def main():
    altitudes_km = [250, 750, 1000]
    profiles = [hammock_profile(h) for h in altitudes_km]

    # Summary table
    print()
    print("=" * 78)
    print("  SOLAR HAMMOCK 24-HOUR AVERAGE POWER  "
          f"(eta = {cfg.CELL_EFFICIENCY*100:.1f}%, AR = {cfg.AR_BOOST:.4f})")
    print("=" * 78)
    header = f"  {'Altitude (km)':<22}"
    for h in altitudes_km:
        header += f"{h:>14}"
    print(header)
    print("  " + "-" * 22 + ("-" * 14) * len(altitudes_km))

    def row(label, key, fmt):
        line = f"  {label:<22}"
        for p in profiles:
            if key == "shadow_deg":
                v = math.degrees(p["phys"]["shadow_ha"])
            elif key == "view_factor":
                v = p["phys"]["F_view"]
            else:
                v = p[key]
            line += f"{v:>14{fmt}}"
        return line

    print(row("View factor (R_E/r)^2", "view_factor", ".4f"))
    print(row("Shadow half-angle", "shadow_deg", ".2f"))
    print(row("Peak (W/m^2)", "p_peak", ".1f"))
    print(row("24-hour avg (W/m^2)", "p_avg", ".1f"))
    print("=" * 78)
    print()

    # Save outputs to multiple locations so the book build and the
    # original "Python" scratch folder both pick up the refreshed plots.
    out_dirs = [
        "/sessions/confident-amazing-lamport/mnt/localbrain/Python/orbitalring/graphs_power",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Python",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Vol_I/images",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Vol_II/images",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Vol_III/images",
    ]
    for d in out_dirs:
        if not os.path.isdir(d):
            os.makedirs(d, exist_ok=True)

    for profile in profiles:
        filename = f"solar_hammock_{profile['h_km']}km_power.png"
        for d in out_dirs:
            plot_single_altitude(profile, os.path.join(d, filename))

    for d in out_dirs:
        plot_comparison(profiles, os.path.join(d, "solar_hammock_24hr_comparison.png"))


if __name__ == "__main__":
    main()
