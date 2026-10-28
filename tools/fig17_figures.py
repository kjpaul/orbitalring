"""
Draws Figures 17.2 and 17.3 of the rewritten Chapter 17 (2026-10-28) and prints the values the captions quote.
  Figure 17.2  fig17_02_descent_profile_700.png: speed (right axis, km/h) and brake power (left axis, MW) of the 100-seat
               pod at 4,000 m/s against distance down the tether, 1 g drive from the coils, cruise, 1.5 g stop. The pod is
               sized as in elevator_system2.py (facing-pair magnets, the shield for 2 mSv from the Chapter 5 dose rates).
  Figure 17.3  fig17_03_throughput_700.png: whole elevator per passenger per day against passengers per day (log axis) from
               elevator_system2's throughput table, six lines (100 and 400 seats at 2,000, 3,000 and 4,000 m/s) with points
               at 60, 30, 15 and 5 minute spacing on the 2,000 K ribbon, fainter lines for the 1,000 K ribbon.
Figure 17.1 (the net pull) is unchanged and stays with Vol_III/renumber_scripts_2026-10-12/fig17_figures.py.
Book: The Orbital Ring, Volume IV: Chapter 17, Figures 17.2 and 17.3.
Run from Python/orbitalring:  python tools/fig17_figures.py <throughput.csv> <output folder>
Environment: ORBITAL_RING_OUT (folder holding elevator_dose_rates.npz; default runs). Same size, fonts and colors as the
figures of 2026-10-21 (10 x 5 in at 200 dpi, DejaVu Sans 11 pt, blue #2b6ca3, red #b5451b, green #2e8b57).
"""
import os, sys, csv, math
sys.path.insert(0, os.getcwd())
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import elevator_system as e1
import elevator_system2 as e2

src = sys.argv[1] if len(sys.argv) > 1 else "runs/elevator_system2/throughput.csv"
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
BLUE, RED, GREEN = e1.BLUE, e1.RED, e1.GREEN
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
r_r, r_geo, g0 = e1.r_r, e1.r_geo, e1.g0
L_LINE = (r_geo - r_r) / 1e3          # km, platform to ring

# ----------------------------------------------------------------------------------------------------------------------
# Figure 17.2: the 100-seat pod at 4,000 m/s
# ----------------------------------------------------------------------------------------------------------------------
V_C, N_SEATS = 4000.0, 100
rows = e2.arrangements()
mag = dict(min([r for r in rows if r["kind"] == "facing pair"], key=lambda r: r["kg_per_N"])); mag["name"] = "NdFeB facing pair"
rates = e1.load_rates(); T_sh = rates["T"]
d1 = e1.descent(V_C, 1000.0); p_, e_, tot = e1.dose_per_descent(rates, d1)
x_shield = e1.shield_for(T_sh, tot, 2.0)
pod = e2.size_pod2(N_SEATS, x_shield, mag)
d = e1.descent(V_C, pod["m"])
S, V, P = d["s_km"], d["v"], d["P"]
j = int(np.argmax(P))
def at(s_km):
    k = int(np.argmin(np.abs(S - s_km))); return P[k] / 1e6, V[k] * 3.6
print("-- Figure 17.2: the 100-seat pod at 4,000 m/s, 1 g drive, 1.5 g stop --")
print(f"  pod {pod['m']/1e3:.1f} t with the shield for 2 mSv ({x_shield:.1f} g/cm2 of water), facing-pair magnets {mag['kg_per_N']*1e6:,.0f} kg per MN")
print(f"  line from the platform to the ring {L_LINE:,.0f} km")
print(f"  drive: {d['t_drive']/60:.1f} min over {d['s_drive']/1e3:,.0f} km to {V_C*3.6:,.0f} km/h ; peak drive power {d['P_drive']/1e6:,.0f} MW ; energy {d['E_drive']/1e9:,.0f} GJ")
s_cruise = L_LINE - d["s_drive"] / 1e3 - d["s_stop"] / 1e3
print(f"  cruise: {d['t_cruise']/3600:.2f} h over {s_cruise:,.0f} km ; brake power at the end of the cruise {d['P_cruise_end']/1e6:,.0f} MW")
print(f"  stop: {d['t_stop']/60:.1f} min over {d['s_stop']/1e3:,.0f} km ; brake force {d['F_stop_start']/1e3:,.0f} kN at the start and {d['F_max']/1e3:,.0f} kN at the ring")
print(f"  peak brake power {P[j]/1e6:,.0f} MW at {S[j]:,.0f} km from the platform (the start of the stop) ; whole trip {d['total']/3600:.2f} h")
for s in (1000, 5000, 10000, 20000, 30000, 34000):
    pw, sp = at(s); print(f"  at {s:>6,} km: brake power {pw:7.0f} MW, speed {sp:,.0f} km/h")
print(f"  brake power at the top (first cruise point) {P[d['phase'] == 1][0]/1e6:.1f} MW ; at the ring {P[-1]/1e6:.1f} MW ; speed at the ring {V[-1]*3.6:.0f} km/h")

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(S, P / 1e6, color=RED, lw=1.8, label="Brake power")
ax.set_xlabel("Distance down the tether from the platform (km)"); ax.set_ylabel("Brake power (MW)", color=RED)
ax.set_xlim(0, L_LINE); ax.set_ylim(0, 8000)
ax.set_xticks(range(0, 35001, 5000)); ax.set_xticklabels([f"{x:,}" for x in range(0, 35001, 5000)])
ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
ax2.plot(S, V * 3.6, color=BLUE, lw=1.8, ls="--", label="Speed")
ax2.set_ylabel("Speed (km/h)", color=BLUE); ax2.set_ylim(0, 18000)
ax2.yaxis.set_major_formatter(mt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.annotate(f"{round(P[j]/1e6, -2):,.0f} MW at the start of the stop", xy=(S[j], P[j] / 1e6), xytext=(18500, 7300), color=RED, arrowprops=dict(arrowstyle="->", color=RED))
ax.text(L_LINE - 300, 80, "Ring", color="gray", ha="right", fontsize=9)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left")
fig.tight_layout(); fn = os.path.join(outdir, "fig17_02_descent_profile_700.png"); fig.savefig(fn, dpi=200); plt.close(fig)
print("  wrote", fn)

# ----------------------------------------------------------------------------------------------------------------------
# Figure 17.3: throughput against spacing, from elevator_system2's table
# ----------------------------------------------------------------------------------------------------------------------
T = []
with open(src) as f:
    for row in csv.DictReader(f):
        T.append({k: float(v) for k, v in row.items()})
SPEEDS, SEATS, SPACINGS = (2000, 3000, 4000), (100, 400), (60, 30, 15, 5)
COLOR = {2000: BLUE, 3000: GREEN, 4000: RED}; MARK = {100: "o", 400: "s"}
def series(v, n, K):
    pts = [next(r for r in T if r["speed_m_s"] == v and r["seats"] == n and r["ribbon_K"] == K and r["spacing_min"] == sp) for sp in SPACINGS]
    return [r["passengers_per_day"] for r in pts], [r["tether_per_passenger_per_day_t"] for r in pts], pts

print("\n-- Figure 17.3: whole elevator per passenger per day (t) at 60 / 30 / 15 / 5 min spacing, 2,000 K ribbon (1,000 K in brackets) --")
fig, ax = plt.subplots(figsize=(10, 5))
for v in SPEEDS:
    for n in SEATS:
        x1, y1, p1 = series(v, n, 1000); x2, y2, p2 = series(v, n, 2000)
        ax.plot(x1, y1, color=COLOR[v], marker=MARK[n], ms=5, lw=1.4, alpha=0.25)
        ax.plot(x2, y2, color=COLOR[v], marker=MARK[n], ms=5, lw=1.8, label=f"{v:,} m/s, {n} seats")
        print(f"  {v:,} m/s, {n:>3} seats: passengers/day " + " / ".join(f"{x:,.0f}" for x in x2) + " ; t per passenger per day "
              + " / ".join(f"{b:.0f} ({a:.0f})" for a, b in zip(y1, y2)) + " ; pods on the line " + " / ".join(f"{r['pods_on_line']:.0f}" for r in p2))
ax.plot([], [], color="gray", alpha=0.3, lw=1.4, label="Same with a 1,000 K ribbon")
# the spacings, written along the heaviest solid line
x0, y0, _ = series(2000, 100, 2000)
for sp, x, y in zip(SPACINGS, x0, y0):
    ax.annotate(f"{sp} min" if sp in (60, 5) else f"{sp}", xy=(x, y), xytext=(6, 5) if sp != 5 else (6, -12), textcoords="offset points", color=BLUE, fontsize=9)
ax.set_xscale("log"); ax.set_xlim(2000, 130000); ax.set_ylim(0, 1400)
ax.set_xticks([2000, 5000, 10000, 20000, 50000, 100000]); ax.set_xticklabels([f"{x:,}" for x in (2000, 5000, 10000, 20000, 50000, 100000)])
ax.xaxis.set_minor_formatter(mt.NullFormatter())
ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.set_xlabel("Passengers per day"); ax.set_ylabel("Whole elevator per passenger per day (t)")
ax.grid(True, alpha=0.3, which="major")
ax.legend(loc="upper right", fontsize=10)
fig.tight_layout(); fn = os.path.join(outdir, "fig17_03_throughput_700.png"); fig.savefig(fn, dpi=200); plt.close(fig)
print("  wrote", fn)
r400 = [next(r for r in T if r["speed_m_s"] == 4000 and r["seats"] == 400 and r["ribbon_K"] == 2000 and r["spacing_min"] == sp) for sp in SPACINGS]
print("  Table 17.9 line (400 seats, 4,000 m/s, 2,000 K): " + " ; ".join(
    f"{r['spacing_min']:.0f} min: {r['pods_on_line']:.0f} pods, {r['largest_load_pod_weights']:.2f} pod weights, {r['tether_whole_t']:,.0f} t, {r['passengers_per_day']:,.0f} pax/day, {r['tether_per_passenger_per_day_t']:.1f} t"
    for r in r400))
print(f"  bottom of the ribbon at the next pod {r400[-1]['bottom_K_at_next_pod']:.0f} K, next peak {r400[-1]['next_peak_bottom_K']:.0f} K at 5 min spacing")
xs = [r["passengers_per_day"] for r in T]; ys = [r["tether_per_passenger_per_day_t"] for r in T]
print(f"  data range: passengers per day {min(xs):,.0f} to {max(xs):,.0f} ; t per passenger per day {min(ys):.1f} to {max(ys):.1f}")
