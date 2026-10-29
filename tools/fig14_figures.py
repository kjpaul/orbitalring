"""
Draws Figures 14.1 (cable and casing wall temperature through the deployment) and 14.2 (the cold-side heat load of a site with
the cryogenic radiator width on the right axis) from the thermal time series of tools/ch14_thermal_series.py, and prints the
values the chapter quotes (day-one peaks, the first months, the settled values, the cold-side peak and its month).
Book: Orbital Ring Engineering, Volume I: Chapter 14, Figures 14.1 and 14.2.
Run: python tools/fig14_figures.py <series.csv> <output folder>
Environment: none.  Copied from Vol_III/renumber_scripts_2026-10-12/fig18_figures.py (written under the old numbering) on 2026-10-27; the month ticks follow the run's length.
"""
import csv, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "ch14_thermal_series.csv"
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
MONTH = 30 * 86400.0

t, Tc, Tw, q, w = [], [], [], [], []
with open(src) as f:
    for row in csv.DictReader(f):
        t.append(float(row["t_s"]) / MONTH)
        Tc.append(float(row["cable_temp_K"]))
        Tw.append(float(row["wall_temp_K"]))
        q.append(float(row["q_cryo_cold_W"]) / 1e3)
        w.append(float(row["radiator_width_m"]))

CABLE, WALL, COLD = "#b5451b", "#2b6ca3", "#2b6ca3"
plt.rcParams.update({"font.size": 11})

# Figure 14.1: cable and wall temperature
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, Tc, color=CABLE, lw=1.6, label="Cable")
ax.plot(t, Tw, color=WALL, lw=1.6, label="Casing wall")
ax.set_xlabel("Time from the start of deployment (months)")
ax.set_ylabel("Temperature (K)")
ax.set_xlim(0, max(t))
ax.set_ylim(200, 390)          # the day-one peak (383 K) must lie inside the axes for its annotation to draw
ax.grid(True, alpha=0.3)
i = Tc.index(max(Tc))
ax.annotate(f"{max(Tc):.0f} K on day 1", xy=(t[i], Tc[i]), xytext=(2.2, 366), color=CABLE,
            arrowprops=dict(arrowstyle="->", color=CABLE))
ax.annotate(f"{Tc[-1]:.0f} K", xy=(t[-1], Tc[-1]), xytext=(t[-1] - 2.5, Tc[-1] + 14), color=CABLE,
            arrowprops=dict(arrowstyle="->", color=CABLE))
ax.annotate(f"{Tw[-1]:.0f} K", xy=(t[-1], Tw[-1]), xytext=(t[-1] - 2.5, Tw[-1] - 22), color=WALL,
            arrowprops=dict(arrowstyle="->", color=WALL))
ax.axhline(273.15, color="gray", lw=0.8, ls=":")
ax.text(0.3, 275, "0°C", color="gray", fontsize=9)
ax.legend(loc="upper right")
ax.set_xticks(range(0, int(max(t)) + 2, 2))
fig.tight_layout()
fig.savefig(os.path.join(outdir, "fig14_01_cable_wall_temperature_700.png"), dpi=200)
plt.close(fig)

# Figure 14.2: cold-side load with radiator width on the right axis
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, q, color=COLD, lw=1.6)
ax.set_xlabel("Time from the start of deployment (months)")
ax.set_ylabel("Cold-side heat load per site (kW)")
ax.set_xlim(0, max(t))
ax.set_ylim(0, 45)
ax.grid(True, alpha=0.3)
j = q.index(max(q))
ax.annotate(f"{max(q):.1f} kW, {w[j]:.2f} m at {t[j]:.1f} months", xy=(t[j], q[j]), xytext=(t[j] + 1.5, q[j] + 2),
            color=COLD, arrowprops=dict(arrowstyle="->", color=COLD))
ax.annotate(f"{q[-1]:.1f} kW, {w[-1]:.2f} m", xy=(t[-1], q[-1]), xytext=(t[-1] - 4.5, q[-1] - 8),
            color=COLD, arrowprops=dict(arrowstyle="->", color=COLD))
# right axis: width is proportional to load
k = w[j] / q[j]
ax2 = ax.twinx()
ax2.set_ylim(0, 45 * k)
ax2.set_ylabel("Cryogenic radiator width over 450 m (m)")
ax.set_xticks(range(0, int(max(t)) + 2, 2))
fig.tight_layout()
fig.savefig(os.path.join(outdir, "fig14_02_cold_load_radiator_700.png"), dpi=200)
plt.close(fig)

print("peak cable", max(Tc), "at month", t[i], "final", Tc[-1], Tw[-1])
print("peak cold", max(q), "width", w[j], "at month", t[j], "final", q[-1], w[-1], "k", k)
# a few checkpoints for the text
for m in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
    idx = min(range(len(t)), key=lambda n: abs(t[n] - m))
    print(f"month {m}: cable {Tc[idx]:.1f} K wall {Tw[idx]:.1f} K cold {q[idx]:.1f} kW width {w[idx]:.2f} m")
# when does the cable first fall below 260 K after the peak?
for n in range(i, len(t)):
    if Tc[n] < 260:
        print("below 260 K at month", round(t[n], 2)); break
print("start", Tc[0], Tw[0], q[0], w[0])
