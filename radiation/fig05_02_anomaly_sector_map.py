"""
Figure 5.4: the anomaly sector on the equator, as a strip map.
Book: Orbital Ring Engineering, Volume I: Chapter 5, Section 5.2 (Figure 5.4).
Run from the radiation folder: python fig05_02_anomaly_sector_map.py [--out=FOLDER]
Environment: none; reads ap8_igrf2026_spectra.npz (AP-8 MIN read with the IGRF 2026 field, made by igrf_altscan.py; the
same file igrf_sector.txt and ch5_figure.txt are printed from). Writes fig05_02_anomaly_sector_map.svg and .png.

A strip of the equator from 180°W to 180°E. No coastline package (cartopy) is installed here, so the land along the
equator is drawn as a labelled bar from the longitudes at which the equator enters and leaves each country or island,
as listed in the Wikipedia article "Equator" (table "Starting at the Prime Meridian and heading eastwards, the equator
passes through", revision of 15 October 2026); the page gives no longitude at which the equator leaves Sumatra and the
Riau and Lingga islands, so that stretch is drawn as islands and straits up to the Borneo crossing. Shaded: the sector
where the trapped proton flux above 10 MeV is nonzero at 700 km at solar minimum, darker where it exceeds 100 and
1,000 protons/cm2/s, with the flux peak marked; a thin line gives the 800 km sector. All longitudes of the sector are
from the run, not drawn by eye.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

BLUE, RED, GREEN, GRAY = "#2b6ca3", "#b5451b", "#2e8b57", "gray"
plt.rcParams.update({"font.size": 12, "font.family": "DejaVu Sans"})

# land on the equator: (longitude entered, longitude left, label), decimal degrees east, from the Wikipedia table
LAND = [
    (-91.583, -91.217, "Galápagos"),
    (-80.117, -49.350, "South America: Ecuador, Colombia, Brazil"),
    (6.519, 6.527, ""),                       # Ilhéu das Rolas, São Tomé and Príncipe
    (9.350, 42.883, "Africa: Gabon to Somalia"),
    (98.217, 109.167, "Sumatra and islands"),     # Sumatra, Riau, Lingga: exit longitude not listed; drawn to the Borneo crossing
    (109.167, 117.517, "Borneo"),
    (119.667, 120.100, ""),                   # Sulawesi
    (127.417, 127.883, ""),                   # Kayoa and Halmahera
    (129.350, 130.183, ""),                   # Gebe and Kawe
]
WATER_IN_LAND = [(32.217, 34.000, "Lake Victoria")]
ISLANDS_AND_STRAITS = [(98.217, 109.167)]      # drawn hatched: the page lists no exit for Sumatra


def sector(d, alt, threshold=0.0):
    """west and east ends, the peak longitude and flux, the number of 1-degree longitudes above the threshold, and those longitudes"""
    J = d[f"min_{alt}"][0]; lon = d["lons"]
    m = J > threshold
    return float(lon[m].min()), float(lon[m].max()), float(lon[J.argmax()]), float(J.max()), int(m.sum()), lon[m]


def shade(ax, lons, **kw):
    """one cell per 1-degree longitude above the threshold, so that any gap the model gives inside the sector shows"""
    for x in lons:
        ax.add_patch(Rectangle((x - 0.5, 0), 1.0, 1, lw=0, **kw))


def lon_label(t):
    return f"{abs(t):.0f}°{'W' if t < 0 else ('E' if t > 0 else '')}"


def main():
    out = next((a[6:] for a in sys.argv[1:] if a.startswith("--out=")), ".")
    os.makedirs(out, exist_ok=True)
    d = np.load("ap8_igrf2026_spectra.npz")
    assert abs(d["Eg"][0] - 10.0) < 1e-9
    w0, e0, pk, jpk, n0, l0 = sector(d, 700, 0.0)
    w100, e100, _, _, n100, l100 = sector(d, 700, 100.0)
    w1000, e1000, _, _, n1000, l1000 = sector(d, 700, 1000.0)
    w8, e8, pk8, _, n8, _ = sector(d, 800, 0.0)
    gaps = [int(x) for x in range(int(w0), int(e0) + 1) if x not in set(l0.astype(int))]
    print(f"700 km, AP-8 MIN, IGRF 2026, flux > 10 MeV: nonzero {w0:.0f} to {e0:.0f} ({n0} degrees of longitude; zero inside the span at {gaps[0]} to {gaps[-1]}), "
          f"above 100: {w100:.0f} to {e100:.0f} ({n100} deg), above 1000: {w1000:.0f} to {e1000:.0f} ({n1000} deg), peak {jpk:.0f} at {pk:.0f}")
    print(f"800 km: nonzero {w8:.0f} to {e8:.0f} ({n8} deg), peak at {pk8:.0f}")

    fig, ax = plt.subplots(figsize=(10, 5.0))
    # the sector shading, three levels, one cell per degree of longitude from the run
    shade(ax, l0, fc=RED, alpha=0.14)
    shade(ax, l100, fc=RED, alpha=0.20)
    shade(ax, l1000, fc=RED, alpha=0.32)
    ax.axvline(pk, color=RED, lw=1.2, ls="--")
    # the 800 km sector as a thin line near the top
    ax.plot([w8, e8], [0.93, 0.93], color=RED, lw=1.4)
    for x in (w8, e8): ax.plot([x, x], [0.91, 0.95], color=RED, lw=1.4)
    ax.text(-178, 0.955, f"800 km sector: {lon_label(w8)} to {lon_label(e8)}", color=RED, fontsize=10, ha="left", va="bottom")
    # the equator and the land along it
    y0, hh = 0.50, 0.10
    ax.plot([-180, 180], [y0, y0], color="#4a7fb5", lw=1.0)
    for a, b, lab in LAND:
        hatched = (a, b) in ISLANDS_AND_STRAITS
        ax.add_patch(Rectangle((a, y0 - hh / 2), b - a, hh, fc="none" if hatched else "#8c8c6e", ec="#5a5a45", lw=0.8,
                               hatch="////" if hatched else None))
    for a, b, lab in WATER_IN_LAND:
        ax.add_patch(Rectangle((a, y0 - hh / 2), b - a, hh, fc="#4a7fb5", ec="none"))
    ax.text(-64.7, y0 - hh / 2 - 0.03, "South America:\nEcuador, Colombia, Brazil", fontsize=10.5, ha="center", va="top")
    ax.text(26.1, y0 - hh / 2 - 0.03, "Africa: Gabon to Somalia\n(Lake Victoria in blue)", fontsize=10.5, ha="center", va="top")
    ax.text(114.2, y0 - hh / 2 - 0.03, "Indonesia: Sumatra,\nBorneo, Sulawesi, Halmahera", fontsize=10.5, ha="center", va="top")
    ax.text(-91.4, y0 + hh / 2 + 0.02, "Galápagos", fontsize=10, ha="center", va="bottom")
    ax.text(6.5, y0 + hh / 2 + 0.02, "São Tomé", fontsize=10, ha="center", va="bottom")
    ax.text(-135, y0 + hh / 2 + 0.02, "Pacific Ocean", fontsize=10.5, ha="center", va="bottom", color="#4a7fb5")
    ax.text(160, y0 + hh / 2 + 0.02, "Pacific Ocean", fontsize=10.5, ha="center", va="bottom", color="#4a7fb5")
    ax.text(70, y0 + hh / 2 + 0.02, "Indian Ocean", fontsize=10.5, ha="center", va="bottom", color="#4a7fb5")
    ax.text(-25, y0 + hh / 2 + 0.02, "Atlantic", fontsize=10.5, ha="center", va="bottom", color="#4a7fb5")
    # the sector labels
    ax.text(pk + 2, 0.88, f"flux peak, {lon_label(pk)}", color=RED, fontsize=10, ha="left", va="top")
    ax.text(-177, 0.86, f"anomaly sector at 700 km,\n{lon_label(w0)} to {lon_label(e0)}:\neastern Pacific, Ecuador,\nColombia, Brazil, Atlantic\nto the Gulf of Guinea", color=RED, fontsize=10, ha="left", va="top")
    ax.text((w1000 + e1000) / 2, 0.16, f"above 1,000\nprotons/cm²/s\n{lon_label(w1000)} to {lon_label(e1000)}", color="white", fontsize=10, ha="center", va="center",
            bbox=dict(fc=RED, ec="none", alpha=0.75, pad=2))
    ax.text(w100 - 2, 0.07, f"above 100: {lon_label(w100)} to {lon_label(e100)}", color=RED, fontsize=10, ha="right", va="center")
    ax.text(96, 0.26, "no trapped protons:\nAfrica, Indian Ocean, Indonesia, the Pacific", color=GREEN, fontsize=10.5, ha="center", va="center")
    ax.set_xlim(-180, 180); ax.set_ylim(0, 1)
    ax.set_xticks(range(-180, 181, 30)); ax.set_xticklabels([lon_label(t) for t in range(-180, 181, 30)])
    ax.set_yticks([])
    ax.set_xlabel("Longitude on the equator")
    ax.set_title("The anomaly sector on the equator: where the ring meets trapped protons (AP-8, field of 2026, solar minimum)", fontsize=12.5)
    fig.subplots_adjust(left=0.03, right=0.99, top=0.91, bottom=0.14)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig05_02_anomaly_sector_map.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig05_02_anomaly_sector_map.svg and .png in {out}")


if __name__ == "__main__":
    main()
