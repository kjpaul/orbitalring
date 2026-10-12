"""
Writes the LaTeX source of the worked equations of Chapters 1 to 4 and 10 with every number taken from
ring_altitude.py (and power_config.py for Chapter 10), and renders them to SVG when latex and dvisvgm are installed.
Book: The Orbital Ring, Volume IV: Equations 1.3 to 1.8, 2.2, 3.6 to 3.10, 4.3 to 4.5, 4.7 to 4.9, 4.11, 10.1, 10.2 (the equation images).
Run: python tools/generate_equations_latex.py   (writes <output folder>/equations/<name>.tex, and .svg when LaTeX is present)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import subprocess, os, shutil
PRE = r"""\documentclass[12pt]{article}
\usepackage{amsmath}
\usepackage[active,tightpage]{preview}
\setlength{\PreviewBorder}{0pt}
\begin{document}
\begin{preview}
\parbox{305.54pt}{\fontsize{14.4pt}{16.5pt}\selectfont\setlength{\jot}{1.5pt}\setlength{\abovedisplayskip}{4pt}\setlength{\belowdisplayskip}{0pt}%
\begin{align*}
%s
\end{align*}}
\end{preview}
\end{document}
"""
# Numbers come from ring_altitude.py, never typed by hand.  The design point is the one set there (700 km, retrograde cable);
# ORBITAL_RING_ALT_KM=2000 reproduces the 2,000 km set of the first edition.
import sys, math
import ring_altitude as R

def c(x, nd=0):
    """Number with LaTeX thousands separators."""
    return f"{x:,.{nd}f}".replace(",", "{,}")

def sci(x, nd=2):
    m, e = f"{x:.{nd}e}".split("e"); return rf"{m} \times 10^{{{int(e)}}}"

ALT = f"{R.ALTITUDE_KM:.0f}"
a_cent = R.V_CABLE**2 / R.R_ORBIT
a_net = a_cent - R.G_LOCAL
a_rel = a_net + R.G_LOCAL
t_gap = math.sqrt(2 * 0.1 / a_rel)
I_sec = R.CABLE_SIDE**4 / 12
E_CNT = 500e9
L_char = math.sqrt(E_CNT * I_sec / R.CABLE_TENSION)
M_ring = R.M_RING_M * R.L_RING
v_lat = R.V_GROUND_SYNC * math.sin(math.radians(30))
KE_lat = 0.5 * M_ring * v_lat**2 * 0.5
r_apo = R.R_ORBIT / (2 * R.GM / (R.R_ORBIT * R.V_CABLE**2) - 1)
r_drag = R.R_EARTH_EQ + 500e3
v_apo_500 = math.sqrt(R.GM * (2 / R.R_ORBIT - 2 / (R.R_ORBIT + r_drag)))  # apogee speed of an ellipse from the ring down to 500 km
TAU_49 = 6400 * 45500 / (7.2e-11 * 6.63e6 * 7755)  # Equation 4.9: full-section fragment, circular orbit at 250 km
# Chapter 10: ring-mounted panels (at the ring's altitude, array offset 0), from power_config.py and power_physics.py
import power_config as PC, power_physics as PP
PC.R_ORBIT = PC.R_EARTH + R.ALTITUDE_KM * 1e3
PC.ALBEDO_PEAK = PC.EARTH_ALBEDO * PC.SOLAR_CONSTANT * (PC.R_EARTH / PC.R_ORBIT) ** 2
P_AVG_M2 = PP.compute_ring_power_profile(1.0)['p_avg_per_m2']   # W/m^2, 24 hour average
N_SITES = round(R.L_RING / 500.0)
P_DEPLOY = N_SITES * 8e6
P_SITE_400 = P_AVG_M2 * 400.0 * 500.0
P_EXPANDED = N_SITES * P_SITE_400
VALUES = dict(P_AVG_M2=P_AVG_M2, N_SITES=N_SITES, P_DEPLOY_GW=P_DEPLOY/1e9, P_SITE_400_MW=P_SITE_400/1e6, P_EXPANDED_TW=P_EXPANDED/1e12, a_cent=a_cent, a_net=a_net, a_rel=a_rel, t_gap=t_gap, v_radial=a_rel*t_gap, I_sec=I_sec, L_char=L_char,
              M_ring=M_ring, v_lat=v_lat, KE_lat=KE_lat, r_apo_km=r_apo/1e3, h_apo_km=(r_apo-R.R_EARTH_EQ)/1e3)

EQS = {
 f"eq03_vorbit_{ALT}": rf"""v_{{\mathrm{{orbit}}}} &= \sqrt{{\frac{{GM}}{{r}}}}\\
&= \sqrt{{\frac{{3.986 \times 10^{{14}}}}{{{c(R.R_ORBIT)}}}}}\\
&= {c(R.V_ORBIT,1)}\ \mathrm{{m/s}}""",
 f"eq04_g_{ALT}": rf"""g(r) &= \frac{{GM}}{{r^{{2}}}}\\
&= \frac{{\mu}}{{r^{{2}}}}\\
&= \frac{{3.986004418 \times 10^{{14}}}}{{({c(R.R_ORBIT)})^{{2}}}}\\
&= {R.G_LOCAL:.4f}\ \mathrm{{m/s}}^{{2}}""",
 f"eq05_vgeo_{ALT}": rf"""v_{{\mathrm{{geo}}}} &= \omega_{{\mathrm{{sidereal}}}} \times r\\
&= 7.2921 \times 10^{{-5}} \times {c(R.R_ORBIT)}\\
&= {c(R.V_GROUND_SYNC,1)}\ \mathrm{{m/s}}""",
 f"eq06_gnet_{ALT}": rf"""g_{{\mathrm{{net}}}} &= g(r) - \omega_{{\mathrm{{sidereal}}}}^{{2}} \cdot r\\
&= {R.G_LOCAL:.3f} - {R.V_GROUND_SYNC**2/R.R_ORBIT:.3f}\\
&= {R.G_NET:.3f}\ \mathrm{{m/s}}^{{2}}""",
 f"eq08_circumference_{ALT}": rf"""C &= 2\pi r\\
&= 2\pi \times {c(R.R_ORBIT)}\\
&= {c(R.L_RING/1e3)}\ \mathrm{{km}}""",
 "ch2_eq03_a_centrifugal": rf"""a_{{\mathrm{{centrifugal}}}} &= \frac{{v_{{\mathrm{{cable}}}}^{{2}}}}{{r}}\\
&= \frac{{({c(R.V_CABLE)})^{{2}}}}{{{c(R.R_ORBIT)}}}\\
&= {a_cent:.3f}\ \mathrm{{m/s}}^{{2}}""",
 "ch2_eq04_a_net_outward": rf"""a_{{\mathrm{{net}}}} &= a_{{\mathrm{{centrifugal}}}} - g(r)\\
&= {a_cent:.3f} - {R.G_LOCAL:.3f}\\
&= {round(a_cent,3) - round(R.G_LOCAL,3):.3f}\ \mathrm{{m/s}}^{{2}}""",
 "ch2_eq05_t_gap": rf"""t &= \sqrt{{\frac{{2d}}{{a_{{\mathrm{{rel}}}}}}}}\\
&= \sqrt{{\frac{{2 \times 0.1}}{{{a_rel:.2f}}}}}\\
&\approx {t_gap:.3f}\ \mathrm{{s}}""",
 "eq09_dv_retro_eq": rf"""\Delta v_{{\mathrm{{load,retro}}}} &= v_{{\mathrm{{orbit}}}} + v_{{\mathrm{{geo}}}}\\
&= {c(R.V_ORBIT,1)} + {c(R.V_GROUND_SYNC,1)}\\
&= {c(R.DELTA_V_LOAD,1)}\ \mathrm{{m/s}}""",
 "eq15_vcable_retro_eq": rf"""v_{{\mathrm{{cable,retro}}}} &= {c(R.V_ORBIT,1)} + \frac{{12{{,}}000 \times {c(R.DELTA_V_LOAD,1)}}}{{{c(R.M_CABLE_M)}}}\\
&= {c(R.V_CABLE,1)}\ \mathrm{{m/s}}""",
 "eq16_tension_retro_eq": rf"""T &= \sigma_{{\mathrm{{cable}}}} \times A\\
&= 12.5 \times 10^{{9}} \times {R.CABLE_AREA:.2f}\\
&= {R.CABLE_TENSION/1e9:.1f}\ \mathrm{{GN}}""",
 "eq17_second_moment": rf"""I &= \frac{{a^4}}{{12}}\\
&= \frac{{({R.CABLE_SIDE:.2f})^4}}{{12}}\\
&= {I_sec:.0f}\ \mathrm{{m}}^4""",
 "eq18_Lchar": rf"""L_{{\mathrm{{char}}}} &= \sqrt{{\frac{{EI}}{{T}}}}\\
&= \sqrt{{\frac{{500 \times 10^{{9}} \times {I_sec:.0f}}}{{{R.CABLE_TENSION/1e9:.1f} \times 10^{{9}}}}}}\\
&= {L_char:.1f}\ \mathrm{{m}}""",
 "eq20_lateral_KE": rf"""KE_{{\mathrm{{lat}}}} &= \frac{{1}}{{2}} m_{{\mathrm{{total}}}}\, v_{{\mathrm{{max}}}}^{{2}}\, \langle \cos^{{2}}(u) \rangle\\
&\approx \frac{{1}}{{2}} \times {sci(M_ring)} \times {v_lat:.0f}^{{2}} \times 0.5\\
&\approx {sci(KE_lat)}\ \mathrm{{J}}""",
 f"ch2_eq09_dv_crit_{ALT}": rf"""\Delta v_{{\mathrm{{crit}}}}({c(R.ALTITUDE_KM)}\ \mathrm{{km}}) &= {c(R.V_ORBIT,1)} - {c(v_apo_500,1)}\\
&= {R.V_ORBIT - v_apo_500:.1f}\ \mathrm{{m/s}}""",
 "eq04_07_apogee": r"""r_{\mathrm{apogee}} &= \frac{r_{\mathrm{perigee}}}{\dfrac{2GM}{r_{\mathrm{perigee}}\, v_{\mathrm{cable}}^{2}} - 1}""",
 f"eq04_08_apogee_{ALT}": rf"""r_{{\mathrm{{apogee}}}} &= \frac{{{c(R.R_ORBIT)}}}{{\dfrac{{2 \times 3.986 \times 10^{{14}}}}{{{c(R.R_ORBIT)} \times ({c(R.V_CABLE)})^{{2}}}} - 1}}\\
&= {c(r_apo/1e3)}\ \mathrm{{km}}\quad ({c((r_apo-R.R_EARTH_EQ)/1e3)}\ \mathrm{{km\ altitude}})""",
 "eq04_09_lifetime": rf"""\tau &\approx \frac{{B\, H}}{{\rho\, r\, v}}\\
&= \frac{{6{{,}}400 \times 45{{,}}500}}{{7.2 \times 10^{{-11}} \times 6.63 \times 10^{{6}} \times 7{{,}}755}}\\
&\approx {sci(TAU_49,1)}\ \mathrm{{s}} \approx {TAU_49/3.156e7:.1f}\ \mathrm{{years}}""",
 "eq11_01_slip_frequency": r"""f_{\mathrm{slip}} &= \frac{v_{\mathrm{wave}} - v_{\mathrm{plate}}}{2\, \tau_{p}}""",
 "eq11_02_lsm_thrust": r"""F &= \frac{B_{\mathrm{stator}}\, B_{\mathrm{sled}}\, \sin\delta}{2\, \mu_{0}}\; A_{\mathrm{active}}""",
 "eq11_03_back_emf": r"""\mathrm{EMF}_{\mathrm{peak}} &= 2\, N\, v\, B_{\mathrm{sled}}\, w_{\mathrm{coil}}""",
 f"eq10_01_deploy_power_{ALT}": rf"""P_{{\mathrm{{deploy}}}} &= N_{{\mathrm{{sites}}}} \times P_{{\mathrm{{site}}}}\\
&= {c(N_SITES)} \times 8\ \mathrm{{MW}}\\
&\approx {P_DEPLOY/1e9:.0f}\ \mathrm{{GW}}""",
 f"eq10_02_expanded_power_{ALT}": rf"""P_{{\mathrm{{expanded}}}} &= N_{{\mathrm{{sites}}}} \times P_{{\mathrm{{site,400}}}}\\
&= {c(N_SITES)} \times {P_SITE_400/1e6:.1f}\ \mathrm{{MW}}\\
&\approx {P_EXPANDED/1e12:.2f}\ \mathrm{{TW}}""",
}
for k, v in VALUES.items(): print(f"  {k:10} {v:,.6g}")
OUT = R.out_dir("equations")   # under ORBITAL_RING_OUT (default runs)
for name, body in EQS.items():
    open(os.path.join(OUT, f"{name}.tex"), "w").write(PRE.replace("%s", body))
    try:
        r = subprocess.run(["latex","-interaction=nonstopmode",f"{name}.tex"],capture_output=True,text=True,cwd=OUT)
    except FileNotFoundError:
        print(name, "NO LATEX on this machine; .tex written"); continue
    if r.returncode: print(name, "LATEX FAIL"); print(r.stdout[-1500:]); continue
    r = subprocess.run(["dvisvgm","--no-fonts","--bbox=preview",f"{name}.dvi","-o",f"{name}.svg"],capture_output=True,text=True,cwd=OUT)
    hdr = open(os.path.join(OUT, f"{name}.svg")).read()
    import re; m=re.search(r"width='([^']*)' height='([^']*)'", hdr); print(name, m.group(1), m.group(2))
