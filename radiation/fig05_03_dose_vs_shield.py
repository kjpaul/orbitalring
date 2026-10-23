"""
Figure 5.3: dose against shield thickness, trapped protons and cosmic rays, at 700 km.
Book: The Orbital Ring, Volume IV: Chapter 5, Sections 5.3 and 5.6 (Figure 5.3).
Run from the radiation folder: python fig05_03_dose_vs_shield.py [--out=FOLDER]
Environment: none; reads tables_trapped_igrf2026.txt (the log of tables_trapped.py: the trapped-proton dose at the worst
longitude, spectrum continued above 300 MeV, the block Table 5.6 of the chapter is taken from) and sweetspot_gcr.txt (the
log of sweetspot_run.py: the cosmic ray dose in the clean part of the ring behind a uniform water shield, Table 5.9).
Writes fig05_03_dose_vs_shield.svg and .png to FOLDER (default: here).

One panel, both axes logarithmic: the dose equivalent in mSv per year against the shield thickness in centimetres of
water (1 g/cm2 = 1 cm). The trapped-proton curve falls at every step; the cosmic ray curve rises to a peak near 2 m of
water, where the shower of secondary particles is largest, before it falls. All thickness points in the logs are plotted.
Style as tools/fig03_02_deployment_tension.py (DejaVu Sans 11 pt, the house blue and red).
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, RED, GREEN, GRAY = "#2b6ca3", "#b5451b", "#2e8b57", "gray"
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
TARGET = 20.0           # mSv/yr, the crew target (Section 5.8)


def read_trapped(path="tables_trapped_igrf2026.txt", alt=700):
    """worst-longitude dose equivalent (Sv/yr) by shield g/cm2, from the 'exponential continuation to 2 GeV' block"""
    lines = open(path).read().splitlines()
    k0 = next(i for i, l in enumerate(lines) if l.startswith(f"{alt} km") and "exponential continuation" in l)
    t, sv = [], []
    for l in lines[k0 + 2:]:
        m = re.match(r"\s+(\d+)\s+\|\s+([\d.]+)\s+([\d.]+)\s+\|", l)
        if not m: break
        t.append(float(m.group(1))); sv.append(float(m.group(3)))
    return np.array(t), np.array(sv)


def read_gcr(path="sweetspot_gcr.txt", alt=700):
    """cosmic ray dose equivalent (mSv/yr) behind a uniform water shield, no cable, clean part of the ring"""
    lines = open(path).read().splitlines()
    k0 = next(i for i, l in enumerate(lines) if l.startswith(f"Ring at {alt} km"))
    t, msv = [], []
    for l in lines[k0 + 1:]:
        m = re.match(r"\s+no cable, uniform shield\s+(\d+) g/cm2 water:\s+([\d.]+)", l)
        if not m:
            if t: break
            continue
        t.append(float(m.group(1))); msv.append(float(m.group(2)))
    return np.array(t), np.array(msv)


def main():
    out = next((a[6:] for a in sys.argv[1:] if a.startswith("--out=")), ".")
    os.makedirs(out, exist_ok=True)
    tt, trapped_sv = read_trapped(); tg, gcr = read_gcr()
    trapped = trapped_sv * 1e3
    print("trapped protons, worst longitude, 700 km (tables_trapped_igrf2026.txt): " + ", ".join(f"{a:.0f} cm {b:.3g} Sv/yr" for a, b in zip(tt, trapped_sv)))
    print("cosmic rays, clean part, 700 km (sweetspot_gcr.txt): " + ", ".join(f"{a:.0f} cm {b:.1f} mSv/yr" for a, b in zip(tg, gcr)))
    k_pk = int(np.argmax(gcr))

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.plot(tt, trapped, color=RED, lw=2.2, marker="o", ms=5, label="Trapped protons, worst longitude in the anomaly sector (56°W)")
    ax.plot(tg, gcr, color=BLUE, lw=2.2, marker="o", ms=5, label="Cosmic rays, clean part of the ring, no cable overhead")
    ax.axhline(TARGET, color=GREEN, lw=1.2, ls="--")
    ax.text(1.05, TARGET * 1.12, f"crew target, {TARGET:.0f} mSv/yr", color=GREEN, fontsize=9.5, ha="left", va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1, 500); ax.set_ylim(10, 2e5)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 200, 500])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "50", "100", "200", "500"])
    ax.set_yticks([10, 100, 1e3, 1e4, 1e5])
    ax.set_yticklabels(["10", "100", "1,000", "10,000", "100,000"])
    ax.set_xlabel("Shield thickness, uniform, in centimetres of water (1 cm = 1 g/cm²)")
    ax.set_ylabel("Dose equivalent (mSv per year)")
    ax.grid(True, which="major", alpha=0.3); ax.grid(True, which="minor", alpha=0.12)
    ax.annotate("the shower peaks at about 2 m of water,\nas it does at 15 to 20 km in the atmosphere",
                xy=(tg[k_pk], gcr[k_pk]), xytext=(25, 260), color=BLUE, fontsize=9.5, ha="center",
                arrowprops=dict(arrowstyle="->", color=BLUE))
    ax.annotate(f"{trapped[-1]:.0f} mSv/yr behind 3 m", xy=(tt[-1], trapped[-1]), xytext=(120, 150), color=RED, fontsize=9.5, ha="center",
                arrowprops=dict(arrowstyle="->", color=RED))
    ax.text(1.05, trapped[0] * 1.25, f"{trapped[0] / 1e3:.0f} Sv/yr behind 1 cm", color=RED, fontsize=9.5, ha="left")
    ax.text(5.3, gcr[0] * 0.70, f"{gcr[0]:.0f} mSv/yr", color=BLUE, fontsize=9, ha="left")
    ax.text(400, 105, f"{gcr[-1]:.0f} mSv/yr\nbehind 5 m", color=BLUE, fontsize=9, ha="center", va="bottom")
    ax.legend(loc="upper right", fontsize=9.5)
    ax.set_title("Dose against shield thickness at 700 km: one curve falls, the other rises before it falls", fontsize=12)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.11)
    for ext, kw in (("svg", {}), ("png", dict(dpi=200))):
        fig.savefig(os.path.join(out, f"fig05_03_dose_vs_shield.{ext}"), **kw)
    plt.close(fig)
    print(f"saved fig05_03_dose_vs_shield.svg and .png in {out}")


if __name__ == "__main__":
    main()
