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
h=np.linspace(200,1500,131)
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
ax.text(1480, 250, "Atmospheric drag region (below 500 km)", color="red", ha="right", va="center", fontsize=10)
ax.axvline(HD, color="gray", linestyle=":", linewidth=1.0)
ax.text(HD+15, 1000, f"{HD:,.0f} km\ndesign altitude", fontsize=9, color="gray", va="bottom")
ax.set_xlabel("Ring Altitude (km)"); ax.set_ylabel("Altitude of the Released Cable's Orbit (km)")
ax.set_title("Where the Cable Goes When It Is Set Free")
ax.set_xlim(200,1500); ax.set_ylim(0,12000); ax.grid(alpha=0.3); ax.legend(loc="upper left")
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
hh=np.linspace(550,1500,96)
fig,ax=plt.subplots(figsize=(7.5,5.5))
for b,col in [(1.0,"red"),(1.5,"blue"),(2.0,"green"),(3.0,"purple")]:
    ys=[ds.particulate_from_pop1(x*1e3,beta=b,casing_frac=0.10)/1e9 for x in hh]
    ax.plot(hh, ys, color=col, linewidth=2.5, label=f"β = {b}" + (" (default)" if b==1.5 else ""))
ax.axhline(50, color="red", linestyle="--", linewidth=1.5, label="Nuclear winter threshold, 50 Mt")
s250=ring_state(250e3); soot250=0.3*0.05*s250['m_cable_total']*s250['L_ring']/1e9
ax.plot([250],[soot250], marker="o", color="black", markersize=7, linestyle="none", label=f"250 km: every fragment's perigee is in the atmosphere ({soot250:.0f} Mt)")
ax.axvspan(0,500,color="red",alpha=0.10)
ax.set_xlim(200,1500); ax.set_ylim(0,90)
ax.set_xlabel("Ring Altitude (km)"); ax.set_ylabel("Prompt Stratospheric Soot from Collision Ejecta (Mt)")
ax.set_title("Population 1: The Only Prompt Reentry")
ax.grid(alpha=0.3); ax.legend(loc="upper right")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_03_ejecta_soot_vs_altitude")

# Fig 4.6 design space
fig,ax=plt.subplots(figsize=(10,4.4))
ax.axvspan(0,600,color="red",alpha=0.18); ax.axvspan(600,800,color="green",alpha=0.15); ax.axvspan(800,1200,color="orange",alpha=0.18)
ax.text(300,0.80,"Below the debris floor\nPrompt soot rises steeply\nSmall fragments reenter\nwithin a century",ha="center",va="center",fontsize=9,color="darkred")
ax.text(700,0.88,"Buildable\nband",ha="center",va="center",fontsize=10,color="darkgreen",fontweight="bold")
ax.text(1000,0.80,"Above the radiation ceiling\nLess than half of the ring\nis free of trapped protons\n(2026 magnetic field)",ha="center",va="center",fontsize=9,color="saddlebrown")
ax.axvline(600,color="darkred",linewidth=1.5,linestyle="--"); ax.text(590,0.30,"600 km\ndebris floor",ha="right",fontsize=9,color="darkred")
ax.axvline(HD,color="black",linewidth=2.5); ax.text(HD-8,0.55,f"{HD:,.0f} km\ndesign\naltitude",ha="right",fontsize=10,fontweight="bold")
ax.axvline(800,color="saddlebrown",linewidth=1.5,linestyle="--"); ax.text(810,0.30,"about 800 km\nradiation ceiling",ha="left",fontsize=9,color="saddlebrown")
ax.axvline(250,color="black",linewidth=1.5,linestyle="--"); ax.text(250,0.22,"250 km\n(baseline of the\nearlier volumes)",ha="center",fontsize=8.5)
ax.set_xlim(0,1200); ax.set_ylim(0,1); ax.set_yticks([]); ax.set_xlabel("Ring Altitude (km)")
ax.set_title("The Altitude Design Space"); ax.grid(axis="x",alpha=0.3)
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_06_design_space")
# Fig 4.5 ring mass across the altitude trade space
alt=np.linspace(250,1500,101)
D=[ring.design(x,True) for x in alt]
fig,axs=plt.subplots(3,1,figsize=(7.5,8.0),sharex=True)
axs[0].plot(alt,[d['m_cable']/1e3 for d in D],color="blue",linewidth=2.5); axs[0].set_ylabel("Cable mass with\nhardware (tonnes/m)")
axs[1].plot(alt,[d['L_ring']/1e6 for d in D],color="green",linewidth=2.5); axs[1].set_ylabel("Ring circumference\n(thousand km)")
axs[2].plot(alt,[d['M_ring']/1e9 for d in D],color="red",linewidth=2.5); axs[2].set_ylabel("Total ring mass (Mt)")
d8=ring.design(HD,True); d25=ring.design(250,True)
for ax,fmt8,fmt25 in [(axs[0],f"{d8['m_cable']/1e3:.0f} t/m",f"{d25['m_cable']/1e3:.0f} t/m"),(axs[1],f"{d8['L_ring']/1e3:,.0f} km",f"{d25['L_ring']/1e3:,.0f} km"),(axs[2],f"{d8['M_ring']/1e9:,.0f} Mt",f"{d25['M_ring']/1e9:,.0f} Mt")]:
    ax.axvline(HD,color="black",linewidth=1.8); ax.grid(alpha=0.3)
    yl=ax.get_ylim(); ym=yl[0]+0.12*(yl[1]-yl[0]); yt=yl[0]+0.80*(yl[1]-yl[0])
    ax.text(HD+20,ym,fmt8+f" at {HD:,.0f} km",fontsize=9); ax.text(265,yt,fmt25+" at 250 km",fontsize=9,color="gray")
axs[2].set_xlabel("Ring Altitude (km)"); axs[0].set_title("Ring Mass Across the Altitude Trade Space (Retrograde Cable)")
axs[2].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}")); axs[2].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_05_ring_mass_vs_altitude")
# Fig 4.4 share of the ring free of trapped protons (AP-8 read with the IGRF 2026 field; radiation/igrf_altscan.py)
import re, os
scan=os.environ.get("IGRF_ALTSCAN", os.path.join(os.path.dirname(os.path.abspath(__file__)), "radiation", "igrf_altscan_clean.txt"))
mn={}; mx={}
for line in open(scan):
    m=re.match(r"\s*(\d+) km solar (min|max):.*clean\s+([\d.]+)%", line)
    if m: (mn if m.group(2)=="min" else mx)[int(m.group(1))]=float(m.group(3))
ks=sorted(mn)
fig,ax=plt.subplots(figsize=(7.5,5.0))
ax.plot(ks,[mn[k] for k in ks],color="blue",linewidth=2.5,marker="o",label="Solar minimum")
ax.plot(ks,[mx[k] for k in ks],color="orange",linewidth=2.5,marker="o",label="Solar maximum")
ax.axhline(50,color="gray",linestyle="--",linewidth=1.2); ax.text(505,51.5,"Half of the ring",fontsize=9,color="gray")
ax.axvline(HD,color="black",linewidth=1.8); ax.text(HD+8,92,f"{HD:,.0f} km design altitude\n{mn[int(HD)]:.0f}% (solar minimum), {mx[int(HD)]:.0f}% (solar maximum)",fontsize=9,va="top")
ax.set_xlim(500,1000); ax.set_ylim(30,100); ax.grid(alpha=0.3); ax.legend(loc="lower left")
ax.set_xlabel("Ring Altitude (km)"); ax.set_ylabel("Share of the Ring With No Trapped Protons (%)")
ax.set_title("The Radiation Ceiling (Magnetic Field of 2026)")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,p: f"{v:,.0f}"))
save(fig,"fig04_04_proton_free_share")
# Kepler orbit diagram for Section 4.4: the released cable's ellipse, drawn to scale at the design altitude
rp=R_EARTH+H_DESIGN; n0=released_cable(H_DESIGN); ra=R_EARTH+n0['h_a']; a_=0.5*(rp+ra); e_=(ra-rp)/(ra+rp); b_=a_*math.sqrt(1-e_**2)
th=np.linspace(0,2*np.pi,400); k=1e-3
fig,ax=plt.subplots(figsize=(8.5,6.2))
ax.add_patch(plt.Circle((0,0),R_EARTH*k,color="#9ec9e8",zorder=2)); ax.text(-1000,3200,"Earth",ha="center",va="center",fontsize=12,zorder=3); ax.text(0,-700,"focus at\nEarth's center",ha="center",va="top",fontsize=9,zorder=3)
ax.plot(rp*k*np.cos(th),rp*k*np.sin(th),color="gray",linestyle="--",linewidth=1.5,label=f"The ring: circular orbit at {HD:,.0f} km")
xc=-(a_-rp)*k   # ellipse centre; Earth at the right-hand focus, perigee on the +x side
ax.plot(xc+a_*k*np.cos(th),b_*k*np.sin(th),color="blue",linewidth=2.5,label="The released cable: elliptical orbit")
ax.plot([rp*k],[0],"o",color="green",markersize=8,zorder=5); ax.plot([-ra*k],[0],"o",color="red",markersize=8,zorder=5)
ax.plot([0],[0],"k+",markersize=10,zorder=4); ax.plot([2*xc],[0],"k+",markersize=8); ax.text(2*xc,-700,"empty focus",ha="center",va="top",fontsize=9)
ax.annotate("",xy=(rp*k,2600),xytext=(rp*k,0),arrowprops=dict(arrowstyle="-|>",color="green",lw=2)); ax.text(rp*k+250,1500,"v$_{cable}$\n(tangential)",color="green",fontsize=10,va="center")
ax.text(rp*k+250,-900,f"Perigee: the break point\nr$_{{perigee}}$ = {rp*k:,.0f} km\n({HD:,.0f} km altitude)",color="green",fontsize=10,va="top")
ax.text(-ra*k-250,900,f"Apogee\nr$_{{apogee}}$ = {ra*k:,.0f} km\n({n0['h_a']*k:,.0f} km altitude)",color="red",fontsize=10,ha="right",va="bottom")
y0=-b_*k-1500
ax.annotate("",xy=(rp*k,y0),xytext=(-ra*k,y0),arrowprops=dict(arrowstyle="<->",color="black",lw=1.2)); ax.text(xc,y0-250,f"2a = r$_{{perigee}}$ + r$_{{apogee}}$ = {2*a_*k:,.0f} km",ha="center",va="top",fontsize=10)
ax.plot([0,rp*k],[0,0],color="green",linewidth=1.2); ax.plot([0,-ra*k],[0,0],color="red",linewidth=1.2)
ax.set_aspect("equal"); ax.set_xlim(-ra*k-4200,rp*k+4300); ax.set_ylim(y0-1500,b_*k+2600); ax.axis("off")
ax.legend(loc="upper center",ncol=2,frameon=False); ax.set_title("The Orbit of the Released Cable")
save(fig,"fig04_kepler_orbit")
print("kepler orbit: rp",rp*k,"ra",ra*k,"a",a_*k,"e",e_,"period h",2*math.pi*math.sqrt(a_**3/GM)/3600)
print("done", soot250)
