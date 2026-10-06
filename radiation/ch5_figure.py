"""Figure 5.1 and the Table 5.2 rows for Chapter 5: trapped proton flux on the equator, AP-8 read with the IGRF 2026
field (ap8_igrf2026_spectra.npz, made by ap8_spectra.py / igrf_altscan.py). Prints the table and writes the figure."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=np.load("ap8_igrf2026_spectra.npz"); E=d["Eg"]; lon=d["lons"]
print("energies MeV:", E)
i10=int(np.argmin(abs(E-10))); i30=int(np.argmin(abs(E-30))); i100=int(np.argmin(abs(E-100)))
print("Table 5.2: altitude | peak >10 | peak >30 | peak >100 | ring mean >10 | clean share solar min | solar max | annual fluence >10 MeV at the worst longitude")
for alt in (500,550,600,650,700,750,800,850,900,1000):
    a=d[f"min_{alt}"]; m=d[f"max_{alt}"]
    print(f"  {alt:5d} | {a[i10].max():7.0f} | {a[i30].max():7.0f} | {a[i100].max():7.0f} | {a[i10].mean():7.1f} | {(a[i10]<=0).mean()*100:5.1f}% | {(m[i10]<=0).mean()*100:5.1f}% | {a[i10].max()*3.156e7:.2e}")
a=d["min_700"][i10]; m=d["max_700"][i10]
print("700 km, solar min, flux >10 MeV by longitude (every 10 deg):")
for L in range(-130,21,10):
    k=int(np.where(lon==L)[0][0]); print(f"  {L:5d}: min {a[k]:7.1f}  max {m[k]:7.1f}")
plt.rcParams.update({"font.size":10,"axes.titlesize":13,"axes.labelsize":12,"legend.fontsize":10,"svg.fonttype":"path"})
order=np.argsort(lon); x=lon[order]
fig,ax=plt.subplots(figsize=(8.5,5.0))
for alt,col in ((800,"red"),(700,"black"),(600,"blue")):
    ax.plot(x,d[f"min_{alt}"][i10][order],color=col,linewidth=2.5 if alt==700 else 1.8,label=f"{alt} km"+(" (the ring)" if alt==700 else (" (the suspended array)" if alt==600 else "")))
ax.plot(x,d["max_700"][i10][order],color="black",linewidth=1.5,linestyle="--",label="700 km at solar maximum")
ax.set_xlim(-180,180); ax.set_ylim(0,3200); ax.set_xticks(range(-180,181,30))
ax.set_xticklabels([f"{abs(t)}°{'W' if t<0 else ('E' if t>0 else '')}" for t in range(-180,181,30)])
ax.set_xlabel("Longitude on the Equator"); ax.set_ylabel("Trapped Protons Above 10 MeV (per cm² per second)")
ax.set_title("Trapped Proton Flux Around the Ring (Magnetic Field of 2026)")
ax.grid(alpha=0.3); ax.legend(loc="upper right")
ax.text(95,1500,"No trapped protons\non this part of the ring",ha="center",fontsize=10,color="darkgreen")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
fig.savefig("fig05_01_proton_flux_around_ring.svg",bbox_inches="tight"); fig.savefig("fig05_01_proton_flux_around_ring.png",dpi=110,bbox_inches="tight")
