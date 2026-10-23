"""
Figure 5.2: the sky seen from the room under the cable.
Book: The Orbital Ring, Volume IV: Chapter 5, Section 5.7 (Figure 5.2).
Run from the radiation folder: python fig05_02_sky_under_cable.py [--out=FOLDER]
Environment: none; reads sweetspot_gcr.txt and sweetspot_trapped.txt in the radiation folder (the logs of
sweetspot_run.py and sweetspot_trapped.py), writes fig05_02_sky_under_the_cable.svg and .png to FOLDER (default: here).

Left: a cross-section looking along the ring: the 8.22 m square cable, the 4 m wide room under it with its ceiling at the
cable's bottom face, and the dose point 1 m above the floor (the geometry of sweetspot.py, h = 1 m below the cable).
Right: the sky around the dose point as a circle, split into the three parts the run gives for that point: covered by
the cable, covered by the Earth, open. The shares are read from the log, not drawn by eye. Beside it, the cosmic ray dose
under the cable against the same room with no cable, and the share of the trapped-proton flux the cable blocks at the
flux peak. Style as tools/fig03_02_deployment_tension.py (DejaVu Sans 11 pt, the house blue and red).
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Circle

BLUE, RED, GREEN, GRAY = "#2b6ca3", "#b5451b", "#2e8b57", "gray"
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
W_CABLE = 8.22          # m, the cable's side (Chapter 3; sweetspot.py uses 8.20 m in the direction sums)
ROOM_W, ROOM_H, Z_POINT = 4.0, 2.0, 1.0   # m: room width, height (ceiling at the cable's bottom face), dose point above the floor


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
    os.makedirs(out, exist_ok=True)
    cable, earth, open_, dose_under, dose_no_cable = read_shares()
    trapped_blocked = read_trapped_share()
    print(f"sky at 1 m below the cable: cable {cable}%, Earth {earth}%, open {open_}% (sweetspot_gcr.txt)")
    print(f"cosmic ray dose, thin walls: {dose_under} mSv/yr under the cable, {dose_no_cable} mSv/yr with no cable")
    print(f"trapped flux through the cable at the flux peak: {trapped_blocked}% (sweetspot_trapped.txt)")

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 5.4), gridspec_kw=dict(width_ratios=[1.0, 1.15], wspace=0.08))

    # ---- left: the cross-section, looking along the ring (east into the page); x is north-south, y up
    ax.set_aspect("equal")
    y_bot = 0.0                                             # the cable's bottom face
    ax.add_patch(Rectangle((-W_CABLE / 2, y_bot), W_CABLE, W_CABLE, fc="#d9d9d9", ec="black", lw=1.2))
    ax.text(0, y_bot + W_CABLE / 2, f"cable\n{W_CABLE} m square\nof carbon", ha="center", va="center", fontsize=10)
    ax.add_patch(Rectangle((-ROOM_W / 2, y_bot - ROOM_H), ROOM_W, ROOM_H, fc="white", ec=BLUE, lw=1.8))
    ax.plot([-ROOM_W / 2 - 3, ROOM_W / 2 + 3], [y_bot - ROOM_H, y_bot - ROOM_H], color=BLUE, lw=1.0)   # the floor line
    ax.text(-ROOM_W / 2 - 0.7, y_bot - ROOM_H - 0.35, "room, 4 m wide,\nthin walls and floor", ha="right", va="top", fontsize=9, color=BLUE)
    ax.plot(0, y_bot - ROOM_H + Z_POINT, "o", color=RED, ms=7, zorder=5)
    ax.annotate("dose point,\n1 m above the floor,\n1 m below the cable", xy=(0, y_bot - ROOM_H + Z_POINT), xytext=(5.3, -1.6),
                fontsize=9, color=RED, ha="left", va="center", arrowprops=dict(arrowstyle="->", color=RED))
    # dimension marks
    ax.annotate("", xy=(W_CABLE / 2 + 0.6, y_bot), xytext=(W_CABLE / 2 + 0.6, y_bot + W_CABLE), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(W_CABLE / 2 + 0.85, y_bot + W_CABLE / 2, f"{W_CABLE} m", color=GRAY, fontsize=9, va="center")
    ax.annotate("", xy=(-ROOM_W / 2, y_bot - ROOM_H - 0.5), xytext=(ROOM_W / 2, y_bot - ROOM_H - 0.5), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(0, y_bot - ROOM_H - 0.95, "4 m", color=GRAY, fontsize=9, ha="center")
    ax.annotate("", xy=(ROOM_W / 2 + 0.5, y_bot - ROOM_H), xytext=(ROOM_W / 2 + 0.5, y_bot), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.8))
    ax.text(ROOM_W / 2 + 0.75, y_bot - ROOM_H / 2, "2 m", color=GRAY, fontsize=9, ha="left", va="center")
    ax.text(0, y_bot + W_CABLE + 0.6, "up (away from the Earth)", ha="center", fontsize=9, color=GRAY)
    ax.text(0, y_bot - ROOM_H - 2.3, "Earth, 700 km below", ha="center", fontsize=9, color=GRAY)
    ax.set_xlim(-7.2, 9.2); ax.set_ylim(-4.6, W_CABLE + 1.4)
    ax.axis("off")
    ax.set_title("Cross-section, looking along the ring", fontsize=11)

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
    bx.text(0, 0.62, f"cable\n{cable:.1f}% of the sky", ha="center", va="center", fontsize=10)
    bx.text(0, -0.62, f"Earth\n{earth:.1f}%", ha="center", va="center", fontsize=10)
    bx.text(-0.66, 0.0, f"open sky\n{open_ / 2:.1f}%\n(up and to\nthe sides)", ha="center", va="center", fontsize=9, color=BLUE)
    bx.text(0.66, 0.0, f"open sky\n{open_ / 2:.1f}%", ha="center", va="center", fontsize=9, color=BLUE)
    bx.text(0, 1.12, f"open sky in total {open_:.1f}%", ha="center", fontsize=9, color=BLUE)
    bx.set_xlim(-1.35, 1.35); bx.set_ylim(-2.05, 1.32)
    bx.set_title("The sky from the dose point, by solid angle", fontsize=11)
    # the numbers beside the sketch
    tbl = (f"Cosmic ray dose, thin walls (model):\n"
           f"   under the cable          {dose_under:.0f} mSv/yr\n"
           f"   same room, no cable   {dose_no_cable:.0f} mSv/yr\n"
           f"Trapped protons at the flux peak:\n"
           f"   blocked by the cable   {trapped_blocked:.0f}% of the flux")
    bx.text(-1.3, -1.25, tbl, ha="left", va="top", fontsize=9.5, family="DejaVu Sans Mono",
            bbox=dict(fc="white", ec=GRAY, lw=0.8, pad=6))

    fig.suptitle("What the room under the cable sees of the sky, at 700 km", fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.04)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig05_02_sky_under_the_cable.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig05_02_sky_under_the_cable.svg and .png in {out}")


if __name__ == "__main__":
    main()
