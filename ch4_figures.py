import math, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import debris_orbits as do
import debris_simulation as ds
from debris_orbits import *
V=V_RECOIL
import ring_altitude as ring
H_DESIGN=ring.ALTITUDE; HD=H_DESIGN/1e3
plt.rcParams.update({"font.size":10, "axes.titlesize":13, "axes.labelsize":12, "legend.fontsize":10, "svg.fonttype":"path"})
def save(fig,name):
    fig.savefig(name+".svg", bbox_inches="tight"); fig.savefig(name+".png", dpi=110, bbox_inches="tight"); plt.close(fig)

# Fig 4.1 released cable orbits vs ring altitude
h=np.linspace(200,3000,141)
apo=[];apo_lo=[];apo_hi=[];peri_w=[]
for hk in h:
    x=hk*1e3; n=released_cable(x); lo=released_cable(x,-V); hi=released_cable(x,V); w=released_cable(x,-V,V)
    apo.append(n['h_a']/1e3); apo_lo.append(lo['h_a']/1e3); apo_hi.append(hi['h_a']/1e3); peri_w.append(w['h_p']/1e3)
fig,ax=plt.subplots(figsize=(7.5,5.5))
ax.fill_between(h, apo_lo, apo_hi, color="blue", alpha=0.12, label="Apogee range, fast to slow half of the ring (±429 m/s)")
ax.plot(h, apo, color="blue", linewidth=2.5, label="Apogee of the released cable (no recoil)")
ax.plot(h, h, color="green", linewidth=2.5, label="Perigee of the released cable (= ring altitude)")
ax.plot(h, peri_w, color="green", linewidth=1.5, linestyle="--", label="Lowest perigee with a 429 m/s radial kick")
ax.axhspan(0, 500, color="red", alpha=0.10)
ax.text(2950, 250, "Atmospheric drag region (below 500 km)", color="red", ha="right", va="center", fontsize=10)
ax.axvline(HD, color="gray", linestyle=":", linewidth=1.0)
ax.text(HD+30, 8600, f"{HD:,.0f} km\ndesign altitude", fontsize=9, color="gray", va="bottom")
ax.axvline(2000, color="gray", linestyle=":", linewidth=1.0)
ax.text(2030, 3600, "2,000 km\nhigh option", fontsize=9, color="gray", va="bottom")
ax.set_xlabel("Ring Altitude (km)"); ax.set_ylabel("Altitude of the Released Cable's Orbit (km)")
ax.set_title("Where the Cable Goes When It Is Set Free")
ax.set_xlim(200,3000); ax.set_ylim(0,14000); ax.grid(alpha=0.3); ax.legend(loc="upper left")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}")); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_01_released_cable_orbits")

# Fig 4.2 lifetime vs perigee (circular orbits), three ballistic coefficients
s=ring_state(H_DESIGN); side=math.sqrt(s['m_struct']/RHO_CNT); Brod=ballistic_coefficient(side)
hp=np.arange(150,1301,25)
fig,ax=plt.subplots(figsize=(7.5,5.5))
for B,lab,col in [(Brod,f"Full cable cross-section, {side:.1f} m (B ≈ {round(Brod,-2):,.0f} kg/m²)","blue"),(600.,"0.8 m piece (B ≈ 600 kg/m²)","purple"),(60.,"8 cm piece (B ≈ 60 kg/m²)","orange")]:
    ys=[]
    for x in hp:
        t=lifetime(x*1e3,x*1e3,B, max_years=1e6)
        ys.append(t if t!=float('inf') else np.nan)
    ax.plot(hp, ys, color=col, linewidth=2.5, label=lab)
for y,lab in [(1,"1 year"),(100,"1 century"),(1e4,"10,000 years")]:
    ax.axhline(y, color="gray", linestyle=":", linewidth=1.0); ax.text(1290, y*1.15, lab, ha="right", fontsize=9, color="gray")
ax.set_yscale("log"); ax.set_ylim(1e-3, 1e6); ax.set_xlim(150,1300)
ax.set_xlabel("Perigee Altitude (km)"); ax.set_ylabel("Time Until Reentry (years)")
ax.set_title("Orbital Lifetime of a Dense CNT Fragment")
ax.grid(alpha=0.3, which="both"); ax.legend(loc="lower right")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_02_fragment_lifetime")

# Fig 4.3 Pop 1 prompt soot vs ring altitude
hh=np.linspace(550,2500,100)
fig,ax=plt.subplots(figsize=(7.5,5.5))
for b,col in [(1.0,"red"),(1.5,"blue"),(2.0,"green"),(3.0,"purple")]:
    ys=[ds.particulate_from_pop1(x*1e3,beta=b,casing_frac=0.10)/1e9 for x in hh]
    ax.plot(hh, ys, color=col, linewidth=2.5, label=f"β = {b}" + (" (default)" if b==1.5 else ""))
ax.axhline(50, color="red", linestyle="--", linewidth=1.5, label="Nuclear winter threshold, 50 Mt")
s250=ring_state(250e3); soot250=0.3*0.05*s250['m_cable_total']*s250['L_ring']/1e9
ax.plot([250],[soot250], marker="o", color="black", markersize=7, linestyle="none", label=f"250 km: every fragment's perigee is in the atmosphere ({soot250:.0f} Mt)")
ax.axvspan(0,500,color="red",alpha=0.10)
ax.set_xlim(200,2500); ax.set_ylim(0,90)
ax.set_xlabel("Ring Altitude (km)"); ax.set_ylabel("Prompt Stratospheric Soot from Collision Ejecta (Mt)")
ax.set_title("Population 1: The Only Prompt Reentry")
ax.grid(alpha=0.3); ax.legend(loc="upper right")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_03_ejecta_soot_vs_altitude")

# Fig 4.6 design space
fig,ax=plt.subplots(figsize=(10,4.4))
ax.axvspan(0,600,color="red",alpha=0.18); ax.axvspan(600,800,color="orange",alpha=0.18); ax.axvspan(800,2500,color="green",alpha=0.15); ax.axvspan(2500,3500,color="red",alpha=0.18)
ax.text(300,0.80,"Released cable\nreaches the\natmosphere\nwithin centuries",ha="center",va="center",fontsize=9,color="darkred")
ax.text(700,0.55,"Margin for solar-maximum atmosphere",ha="center",va="center",fontsize=8.5,color="saddlebrown",rotation=90)
ax.text(1650,0.82,"Buildable band",ha="center",va="center",fontsize=11,color="darkgreen",fontweight="bold")
ax.text(1650,0.68,"Released cable stays in orbit for millennia\nRadiation manageable with shielding",ha="center",va="center",fontsize=9.5,color="darkgreen")
ax.text(3000,0.80,"Van Allen ceiling\nSolar panels degrade\n10% in weeks",ha="center",va="center",fontsize=9,color="darkred")
ax.axvline(HD,color="black",linewidth=2.5); ax.text(HD+40,0.36,f"{HD:,.0f} km\ndesign altitude",ha="left",fontsize=10,fontweight="bold")
ax.axvline(2000,color="black",linewidth=1.5,linestyle="--"); ax.text(2000,0.22,"2,000 km\n(high option, at the\nedge of the inner belt)",ha="center",fontsize=8.5)
ax.axvline(250,color="black",linewidth=1.5,linestyle="--"); ax.text(250,0.22,"250 km\n(baseline of the\nearlier volumes)",ha="center",fontsize=8.5)
ax.set_xlim(0,3500); ax.set_ylim(0,1); ax.set_yticks([]); ax.set_xlabel("Ring Altitude (km)")
ax.set_title("The Altitude Design Space"); ax.grid(axis="x",alpha=0.3)
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_06_design_space")
# Fig 4.5 ring mass across the altitude trade space
alt=np.linspace(250,3000,111)
D=[ring.design(x,True) for x in alt]
fig,axs=plt.subplots(3,1,figsize=(7.5,8.0),sharex=True)
axs[0].plot(alt,[d['m_cable']/1e3 for d in D],color="blue",linewidth=2.5); axs[0].set_ylabel("Cable mass with\nhardware (tonnes/m)")
axs[1].plot(alt,[d['L_ring']/1e6 for d in D],color="green",linewidth=2.5); axs[1].set_ylabel("Ring circumference\n(thousand km)")
axs[2].plot(alt,[d['M_ring']/1e9 for d in D],color="red",linewidth=2.5); axs[2].set_ylabel("Total ring mass (Mt)")
d8=ring.design(HD,True); d20=ring.design(2000,True)
for ax,fmt8,fmt20 in [(axs[0],f"{d8['m_cable']/1e3:.0f} t/m",f"{d20['m_cable']/1e3:.0f} t/m"),(axs[1],f"{d8['L_ring']/1e3:,.0f} km",f"{d20['L_ring']/1e3:,.0f} km"),(axs[2],f"{d8['M_ring']/1e9:,.0f} Mt",f"{d20['M_ring']/1e9:,.0f} Mt")]:
    ax.axvline(HD,color="black",linewidth=1.8); ax.axvline(2000,color="gray",linestyle="--",linewidth=1.2); ax.grid(alpha=0.3)
    yl=ax.get_ylim(); ym=yl[0]+0.12*(yl[1]-yl[0])
    ax.text(HD+40,ym,fmt8+f" at {HD:,.0f} km",fontsize=9); ax.text(2040,ym,fmt20+" at 2,000 km",fontsize=9,color="gray")
axs[2].set_xlabel("Ring Altitude (km)"); axs[0].set_title("Ring Mass Across the Altitude Trade Space (Retrograde Cable)")
axs[2].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}")); axs[2].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_05_ring_mass_vs_altitude")
print("done", soot250)
