#!/usr/bin/env python3
"""
Regenerate the bifacial solar panel 24-hour output graph for Chapter 7
"Powering the Orbital Ring" (now in Volume III).

This script produces the "bifacial_solar_24hr_output.png" figure that
shows the electrical output of a 1 m^2 bifacial solar panel riding on
the orbital ring over a 24-hour local solar day.

Physics model: restored from Chapter 7 of "Orbital Ring Engineering"
by Paul G. de Jong. Uses:
  - Front-face direct flux: S0 * max(0, cos(phi))
  - Back-face albedo flux:  A * S0 * (R_E/r)^2 * (1+cos(phi))/2
  - 45% efficient multi-junction cells
  - Anti-reflection nanostructure boost of 1.0785x
  - Peak at noon:          844 W/m^2
  - 24-hour average:       302 W/m^2

The x-axis is local solar time (0 to 24 hours). The y-axis is the
instantaneous electrical output in W/m^2. Yellow fills the direct
contribution, blue stacks the albedo contribution on top, the red
dashed line marks the 24-hour average, and the shadow regions are
shaded gray with "Earth's Shadow" callouts.
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import power_config as cfg
import power_physics as pp


def local_time_from_phi(phi_rad):
    """Convert angle from solar noon (rad) to local solar time (hours)."""
    # phi = 0 -> noon (12 h); phi = +/-pi -> midnight (0 h or 24 h)
    # Hours past midnight = 12 + phi * 12/pi
    return 12.0 + phi_rad * 12.0 / math.pi


def main():
    # Sample the full day.
    # phi ranges from -pi to +pi, corresponding to 0:00 to 24:00 local time.
    n = 2000
    phi = np.linspace(-math.pi, math.pi, n)

    # Compute direct and albedo electrical outputs separately so they
    # can be plotted as stacked filled regions.
    eff = cfg.CELL_EFFICIENCY * cfg.AR_BOOST * cfg.PANEL_PACKING
    p_direct = np.array([pp.direct_solar_flux(p) for p in phi]) * eff
    p_albedo = np.array([pp.albedo_flux(p) for p in phi]) * eff
    p_total = p_direct + p_albedo

    t_hours = np.array([local_time_from_phi(p) for p in phi])

    # 24-hour average. Recompute with a fine uniform sum so the result
    # is independent of the plot resolution.
    n_avg = 360_000
    phi_fine = np.linspace(-math.pi, math.pi, n_avg, endpoint=False)
    p_avg = np.mean([pp.panel_electrical_output(p) for p in phi_fine])
    p_peak = float(np.max(p_total))

    print(f"Peak electrical output: {p_peak:.2f} W/m^2")
    print(f"24-hour average:        {p_avg:.2f} W/m^2")

    # Shadow region boundaries in local solar time.
    # Sunlit region spans phi in [-shadow_half, +shadow_half] where
    # shadow_half = 105.8 deg; outside that arc the panel is in shadow.
    t_shadow_start_morning = local_time_from_phi(-cfg.SHADOW_HALF_ANGLE)
    t_shadow_start_evening = local_time_from_phi(+cfg.SHADOW_HALF_ANGLE)

    # --- Plot -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 5))

    # Stacked filled areas
    direct_color = "#f4d03f"   # warm yellow for front-face direct sun
    albedo_color = "#5dade2"   # ocean blue for back-face Earth albedo

    ax.fill_between(t_hours, 0.0, p_direct,
                    facecolor=direct_color, edgecolor="#c8a93a",
                    linewidth=0.7, label="Direct sunlight (front surface)",
                    zorder=3)
    ax.fill_between(t_hours, p_direct, p_total,
                    facecolor=albedo_color, edgecolor="#3f7ca8",
                    linewidth=0.7, label="Earth albedo (back surface)",
                    zorder=3)

    # Total curve (outline on top of the stack)
    ax.plot(t_hours, p_total, color="#1a2e4a", linewidth=1.1,
            label="Total electrical output", zorder=4)

    # 24-hour average line
    ax.axhline(p_avg, color="#c0392b", linestyle="--", linewidth=1.6,
               label=f"24-hour average: {p_avg:.0f} W/m\u00b2", zorder=5)

    # Shadow shading
    shadow_color = "#d0d0d0"
    ax.axvspan(0.0, t_shadow_start_morning,
               color=shadow_color, alpha=0.55, zorder=1)
    ax.axvspan(t_shadow_start_evening, 24.0,
               color=shadow_color, alpha=0.55, zorder=1)

    # "Earth's Shadow" text callouts
    ax.text(t_shadow_start_morning / 2.0, 150.0, "Earth's\nShadow",
            ha="center", va="center", color="#555555",
            fontsize=10, fontstyle="italic", zorder=6)
    ax.text((t_shadow_start_evening + 24.0) / 2.0, 150.0, "Earth's\nShadow",
            ha="center", va="center", color="#555555",
            fontsize=10, fontstyle="italic", zorder=6)

    # Axes and labels
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
        "24-Hour Power Output of 1 m\u00b2 Bifacial Solar Panel on Orbital Ring\n"
        "(250 km altitude, 45% efficient cells with nanostructured AR coating)",
        fontsize=11,
    )

    ax.grid(True, which="major", linestyle="-", linewidth=0.4,
            color="#dddddd", zorder=0)
    ax.legend(loc="upper right", framealpha=0.95, fontsize=9)

    plt.tight_layout()

    # Save to the local graphs_power folder AND to the canonical
    # localbrain image locations used by the book build.
    out_dirs = [
        "/sessions/confident-amazing-lamport/mnt/localbrain/Python/orbitalring/graphs_power",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Vol_II/images",
        "/sessions/confident-amazing-lamport/mnt/localbrain/Vol_III/images",
    ]
    filename = "bifacial_solar_24hr_output.png"
    for d in out_dirs:
        if not os.path.isdir(d):
            os.makedirs(d, exist_ok=True)
        path = os.path.join(d, filename)
        fig.savefig(path, dpi=300, bbox_inches="tight",
                    facecolor="white")
        print(f"  saved: {path}")

    plt.close(fig)


if __name__ == "__main__":
    main()
