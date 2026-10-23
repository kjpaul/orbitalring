"""
Figure 5.2: the sky seen from the room under the cable.
Book: The Orbital Ring, Volume IV: Chapter 5, Section 5.7 (Figure 5.2).
Run from the radiation folder: python fig05_02_sky_under_cable.py [--out=FOLDER] [--h=1.0]
Environment: none; reads sweetspot_gcr.txt and sweetspot_trapped.txt in the radiation folder (the logs of
sweetspot_run.py and sweetspot_trapped.py), writes fig05_02_sky_under_the_cable.svg and .png to FOLDER (default: here).
--h=2.0 draws the version with the ceiling 1 m below the cable (the EML track and its cryostat between them) and the
dose point 2 m below the cable, from the h = 2.0 m row of the log, as fig05_02_sky_under_the_cable_h2.svg; the
trapped-proton share blocked by the cable is in the log for h = 1 m only, so that line is left off the h = 2 m version.

Left: a cross-section looking along the ring: the 8.22 m square cable, the 4 m wide room under it with its ceiling at the
cable's bottom face, and the dose point 1 m above the floor (the geometry of sweetspot.py, h = 1 m below the cable).
Right: the sky around the dose point as a circle, split into the three parts the run gives for that point: covered by
the cable, covered by the Earth, open. The shares are read from the log, not drawn by eye. Beside it, the cosmic ray dose
under the cable against the same room with no cable, and the share of the trapped-proton flux the cable blocks at the
flux peak. Style as tools/fig03_02_deployment_tension.py (DejaVu Sans 11 pt, the house blue and red).
"""
import math, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Circle

BLUE, RED, GREEN, GRAY = "#2b6ca3", "#b5451b", "#2e8b57", "gray"
plt.rcParams.update({"font.size": 12, "font.family": "DejaVu Sans"})
W_CABLE = 8.22          # m, the cable's side (Chapter 3; sweetspot.py uses 8.20 m in the direction sums)
ROOM_W, ROOM_H, Z_POINT = 4.0, 2.0, 1.0   # m: room width, height, dose point above the floor; the ceiling is h - 1 m below the cable


def read_shares(path="sweetspot_gcr.txt", alt=700, h=1.0):
    """the sky shares and the doses of the h row of the 'under the cable' table of the given altitude block"""
    lines = open(path).read().splitlines()
    k0 = next(i for i, l in enumerate(lines) if l.startswith(f"Ring at {alt} km"))
    block = lines[k0:k0 + 20]
    no_cable_5 = next(float(l.split(":")[1].split("(")[0]) for l in block if "uniform shield    5 g/cm2" in l)
    row = next(l for l in block if re.match(rf"\s+{h:.1f} m\s+\|", l))
    m = re.search(r"\|\s*([\d.]+)% / ([\d.]+)% / ([\d.]+)%\s*\|\s*([\d.]+)", row)
    cable, earth, open_, dose_thin = (float(m.group(i)) for i in (1, 2, 3, 4))
    return cable, earth, open_, dose_thin, no_cable_5


def read_trapped_share(path="sweetspot_trapped.txt", alt=700):
    lines = open(path).read().splitlines()
    k0 = next(i for i, l in enumerate(lines) if l.startswith(f"Ring at {alt} km"))
    line = next(l for l in lines[k0:] if "share of the trapped flux arriving through the cable" in l)
    return float(re.search(r"through the cable (\d+)%", line).group(1))


def main():
    out = next((a[6:] for a in sys.argv[1:] if a.startswith("--out=")), ".")
    h = float(next((a[4:] for a in sys.argv[1:] if a.startswith("--h=")), "1.0"))
    os.makedirs(out, exist_ok=True)
    cable, earth, open_, dose_under, dose_no_cable = read_shares(h=h)
    trapped_blocked = read_trapped_share() if h == 1.0 else None
    suffix = "" if h == 1.0 else f"_h{h:.0f}"
    print(f"sky at {h:.0f} m below the cable: cable {cable}%, Earth {earth}%, open {open_}% (sweetspot_gcr.txt)")
    print(f"cosmic ray dose, thin walls: {dose_under} mSv/yr under the cable, {dose_no_cable} mSv/yr with no cable")
    print(f"trapped flux through the cable at the flux peak: {trapped_blocked}% (sweetspot_trapped.txt, h = 1 m only)" if trapped_blocked else "trapped-proton share: not in the log for this h, line left off")

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 5.5), gridspec_kw=dict(width_ratios=[1.0, 1.15], wspace=0.08))

    # ---- left: the cross-section, looking along the ring (east into the page); x is north-south, y up
    ax.set_aspect("equal")
    y_bot = 0.0                                             # the cable's bottom face
    y_ceil = y_bot - (h - Z_POINT)                          # the room's ceiling: at the cable for h = 1 m, 1 m lower for h = 2 m
    y_floor = y_ceil - ROOM_H
    ax.add_patch(Rectangle((-W_CABLE / 2, y_bot), W_CABLE, W_CABLE, fc="#d9d9d9", ec="black", lw=1.2))
    ax.text(0, y_bot + W_CABLE / 2, f"cable\n{W_CABLE} m square\nof carbon", ha="center", va="center", fontsize=11)
    if y_ceil < y_bot - 1e-9:
        ax.add_patch(Rectangle((-ROOM_W / 2 - 1.0, y_ceil), ROOM_W + 2.0, y_bot - y_ceil, fc="#eeeeee", ec=GRAY, lw=0.8, ls="--"))
        ax.text(-ROOM_W / 2 - 1.2, (y_bot + y_ceil) / 2, "EML track and\nits cryostat", color=GRAY, fontsize=10, ha="right", va="center")
    ax.add_patch(Rectangle((-ROOM_W / 2, y_floor), ROOM_W, ROOM_H, fc="white", ec=BLUE, lw=1.8))
    ax.plot([-ROOM_W / 2 - 3, ROOM_W / 2 + 3], [y_floor, y_floor], color=BLUE, lw=1.0)   # the floor line
    ax.text(-ROOM_W / 2 - 0.7, y_floor - 0.35, "room, 4 m wide,\nthin walls and floor", ha="right", va="top", fontsize=10, color=BLUE)
    ax.plot(0, y_floor + Z_POINT, "o", color=RED, ms=7, zorder=5)
    ax.annotate(f"dose point,\n1 m above the floor,\n{h:.0f} m below the cable", xy=(0, y_floor + Z_POINT), xytext=(5.3, y_floor + Z_POINT - 0.6),
                fontsize=10, color=RED, ha="left", va="center", arrowprops=dict(arrowstyle="->", color=RED))
    # dimension marks
    ax.annotate("", xy=(W_CABLE / 2 + 0.6, y_bot), xytext=(W_CABLE / 2 + 0.6, y_bot + W_CABLE), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(W_CABLE / 2 + 0.85, y_bot + W_CABLE / 2, f"{W_CABLE} m", color=GRAY, fontsize=10, va="center")
    ax.annotate("", xy=(-ROOM_W / 2, y_floor - 0.5), xytext=(ROOM_W / 2, y_floor - 0.5), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(0, y_floor - 0.95, "4 m", color=GRAY, fontsize=10, ha="center")
    ax.annotate("", xy=(ROOM_W / 2 + 0.5, y_floor), xytext=(ROOM_W / 2 + 0.5, y_ceil), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(ROOM_W / 2 + 0.75, (y_floor + y_ceil) / 2, "2 m", color=GRAY, fontsize=10, ha="left", va="center")
    ax.text(0, y_bot + W_CABLE + 0.6, "up (away from the Earth)", ha="center", fontsize=10, color=GRAY)
    ax.text(0, y_floor - 2.3, "Earth, 700 km below", ha="center", fontsize=10, color=GRAY)
    ax.set_xlim(-7.2, 9.6); ax.set_ylim(y_floor - 2.8, W_CABLE + 1.4)
    ax.axis("off")
    ax.set_title("Cross-section, looking along the ring", fontsize=12)

    # ---- right: the sky as a circle, split by solid angle: cable on top, Earth below, open sky at the sides
    bx.set_aspect("equal"); bx.axis("off")
    r = 1.0
    deg_cable, deg_earth = 360 * cable / 100, 360 * earth / 100
    deg_side = (360 - deg_cable - deg_earth) / 2
    a0 = 90 - deg_cable / 2                                 # the cable wedge is centred on straight up
    bx.add_patch(Wedge((0, 0), r, a0, a0 + deg_cable, fc="#d9d9d9", ec="black", lw=1.0))
    bx.add_patch(Wedge((0, 0), r, a0 + deg_cable, a0 + deg_cable + deg_side, fc="#dbe7f3", ec="black", lw=1.0))     # open, one side
    bx.add_patch(Wedge((0, 0), r, 270 - deg_earth / 2, 270 + deg_earth / 2, fc="#c9b8a8", ec="black", lw=1.0))     # Earth, centred on down
    bx.add_patch(Wedge((0, 0), r, 270 + deg_earth / 2, 270 + deg_earth / 2 + deg_side, fc="#dbe7f3", ec="black", lw=1.0))   # open, other side
    bx.add_patch(Circle((0, 0), r, fc="none", ec="black", lw=1.2))
    bx.plot(0, 0, "o", color=RED, ms=6, zorder=5)
    bx.text(0, 0.60, f"cable\n{cable:.1f}% of the sky", ha="center", va="center", fontsize=11)
    bx.text(0, -0.60, f"Earth\n{earth:.1f}%", ha="center", va="center", fontsize=11)
    # the two open-sky parts: small labels outside the circle with leader lines, so that nothing crosses a wedge line
    a_side = math.radians(a0 + deg_cable + deg_side / 2)      # the middle of the left open wedge
    for sgn, lab in ((-1, f"open sky\n{open_ / 2:.1f}%\n(up and to\nthe sides)"), (1, f"open sky\n{open_ / 2:.1f}%")):
        xm, ym = sgn * abs(math.cos(a_side)) * 0.72, math.sin(a_side) * 0.72
        bx.annotate(lab, xy=(xm, ym), xytext=(sgn * 1.32, 0.30), fontsize=10, color=BLUE, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.8))
    bx.text(0, 1.12, f"open sky in total {open_:.1f}%", ha="center", fontsize=10, color=BLUE)
    bx.set_xlim(-1.75, 1.75); bx.set_ylim(-1.95, 1.32)
    bx.set_title("The sky from the dose point, by solid angle", fontsize=12)
    # the numbers beside the sketch
    tbl = (f"Cosmic ray dose, thin walls (model):\n"
           f"   under the cable          {dose_under:.0f} mSv/yr\n"
           f"   same room, no cable   {dose_no_cable:.0f} mSv/yr")
    if trapped_blocked is not None:
        tbl += (f"\nTrapped protons at the flux peak:\n"
                f"   blocked by the cable   {trapped_blocked:.0f}% of the flux")
    bx.text(-1.5, -1.25, tbl, ha="left", va="top", fontsize=10.5, family="DejaVu Sans Mono",
            bbox=dict(fc="white", ec=GRAY, lw=0.8, pad=6))

    fig.suptitle(f"What the room under the cable sees of the sky, at 700 km (dose point {h:.0f} m below the cable)", fontsize=12.5, y=0.99)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.04)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig05_02_sky_under_the_cable{suffix}.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig05_02_sky_under_the_cable{suffix}.svg and .png in {out}")


if __name__ == "__main__":
    main()
