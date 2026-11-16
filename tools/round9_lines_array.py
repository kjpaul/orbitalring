#!/usr/bin/env python3
"""
Round 9, item 5c: the lines and the array (Chapters 11 and 12 of Volume III).  Three checks:
  1. the tent of the suspended array: the hanging sheet as a catenary, against the two parabola readings the chapter scripts mix;
  2. the pull at the top of the array support lines (4.7 against 4.8 kN/m);
  3. the reaction lines as hanging lines (constant stress, gravity that changes with height, the centrifugal term, a round Earth),
     against the straight line at 30 degrees that the chapter script sizes.
Book: Orbital Ring Engineering, Volume III: Chapter 11 (Table 11.1, Table 11.2, Sections 11.3 and 11.6), Chapter 12 (Sections 12.3, 12.4, 12.8, Table 12.1).
Run from Python/orbitalring:  python tools/round9_lines_array.py
Log: runs/round9_logs/round9_lines_array.txt and .json
Imports ring_altitude.py, power_config.py and power_physics.py unchanged.  The formulas of the chapter scripts (Vol_III/renumber_scripts_2026-10-12/
ch16_numbers.py sections P, Q and R, ch15_arm_choices.py, ch15_numbers.py line(), spreader_cases.py), which sit outside this repository and run on import,
are written out again here and first checked against their logged numbers, so that a difference below is a difference of the model and not of the inputs.
"""
import io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root
import ring_altitude as R, power_config as cfg, power_physics as phys

LOG_DIR = os.path.join(ROOT, "runs", "round9_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n")
OUT = {}

mu = R.GM; RE = R.R_EARTH_EQ; W_E = R.OMEGA_SIDEREAL
RHO = R.RHO_CNT; S_OP = R.SIGMA_CABLE; S_BR = 25e9; S75 = 0.75*S_BR; E_CNT = R.E_CNT
h_r = R.ALTITUDE; h_a = R.ARRAY_ALTITUDE; L = R.L_RING; La = R.L_ARRAY
r_top = RE + h_r

def geff(h):
    """Effective gravity of a body that turns with the Earth, in the equatorial plane (gravity less the centrifugal term), m/s2."""
    r = RE + h; return mu/r**2 - W_E*W_E*r
def pot(h1, h2):
    """Potential difference from altitude h1 up to h2 for a body that turns with the Earth, J/kg."""
    r1, r2 = RE + h1, RE + h2
    return mu*(1/r1 - 1/r2) - 0.5*W_E*W_E*(r2*r2 - r1*r1)
def dia(T, s): return 2e3*math.sqrt(T/s/math.pi)   # mm
g_a = geff(h_a); g_r = geff(h_r)

P("Round 9 item 5c: the lines and the array.  ring %.0f km, array %.0f km ; g_eff at the ring %.4f, at the array %.4f, at the ground %.4f m/s2" % (h_r/1e3, h_a/1e3, g_r, g_a, geff(0)))
P("CNT: density %.0f kg/m3 ; anchor lines and arms %.2f GPa ; array support lines and retention lines %.2f GPa ; modulus %.0f GPa" % (RHO, S_OP/1e9, S75/1e9, E_CNT/1e9))

# =====================================================================================================================
# PART 1.  THE TENT
# =====================================================================================================================
P("\n" + "="*110 + "\nPART 1.  The tent: which width is which\n" + "="*110)

def catenary(theta_deg):
    """Uniform sheet (weight per arc length), flat at the low edge, slope theta at the ridge, arc length s per side.
    a = s/tan(theta) ; half-span = a asinh(tan theta) ; drop = a (sec theta - 1).  Returns half-span/s, drop/s."""
    t = math.tan(math.radians(theta_deg))
    return math.asinh(t)/t, (math.sqrt(1 + t*t) - 1)/t
def catenary_numeric(theta_deg, n=200000):
    """The same by adding up the sheet: constant horizontal pull H, vertical pull q s' at arc s' from the edge, so tan(psi) = (s'/s) tan(theta).
    half-span = sum cos(psi) ds', drop = sum sin(psi) ds'.  Independent of the closed form."""
    t = math.tan(math.radians(theta_deg)); u = (np.arange(n) + 0.5)/n; psi = np.arctan(u*t)
    return float(np.cos(psi).mean()), float(np.sin(psi).mean())
def parabola(theta_deg):
    """Sheet loaded per horizontal meter (a parabola), flat at the edge, slope theta at the ridge, half-span x0.
    drop = x0 tan(theta)/2 ; arc = x0 (t sqrt(1+t2) + asinh t)/(2 t).  Returns arc/x0, drop/x0."""
    t = math.tan(math.radians(theta_deg))
    return (t*math.sqrt(1 + t*t) + math.asinh(t))/(2*t), t/2

P("\n1a. The catenary factors, closed form against the sum (for a sheet of arc length s per side)")
P("  ridge slope | half-span / s (closed, summed) | drop / s (closed, summed) | sheet / footprint | parabola: sheet / footprint, drop / half-span")
OUT["catenary"] = {}
for th in (30, 45, 60):
    hs, dr = catenary(th); hs2, dr2 = catenary_numeric(th); fa, pd = parabola(th)
    P("  %2d deg      | %.4f, %.4f                 | %.4f, %.4f             | %.4f            | %.4f, %.4f" % (th, hs, hs2, dr, dr2, 1/hs, fa, pd))
    OUT["catenary"][th] = dict(half_span_over_s=hs, drop_over_s=dr, sheet_over_footprint=1/hs, parabola_sheet_over_footprint=fa)

# ---- the design case of ch16_numbers.py section Q (v2), written out again ----
dP_arr = pot(h_a, h_r)
kc = RHO*(h_r - h_a)*g_a/(S75 - RHO*dP_arr)       # ch16_numbers.py section Q: lines of one thickness
S_SP = 10e3
def v2(sA, M=1200.0, S=S_SP, lam=45.0):
    """ch16_numbers.py section Q, function v2, unchanged: the strip that fits the budget M (kg per m of ring)."""
    q = sA*g_a; TR = q*S*S/2; mR = 2*RHO*TR/S75
    per_w = sA*La/L*(1 + kc) + q/(g_r*math.tan(math.radians(lam)))
    Wd = (M - mR)/per_w; H = q*Wd/2; FR = H*S
    Tlow = FR/math.sin(math.radians(lam)); Tup = FR/math.tan(math.radians(lam))
    return dict(q=q, TR=TR, mR=mR, W=Wd, H=H, FR=FR, Tlow=Tlow, Tup=Tup, pan=sA*La/L*(1 + kc)*Wd, arms=q*Wd/(g_r*math.tan(math.radians(lam))))
d5 = v2(0.5)
P("\n1b. The design case again (ch16_numbers.py section Q, v2): width %.1f m ; panels and lines %.1f + arms %.1f + retention %.1f = %.1f kg/m ; edge pull %.3f kN/m ; sideways %.2f MN ; arm below the kink %.2f MN"
  % (d5['W'], d5['pan'], d5['arms'], d5['mR'], d5['pan'] + d5['arms'] + d5['mR'], d5['H']/1e3, d5['FR']/1e6, d5['Tlow']/1e6))
P("    logged by ch16_numbers.py: width 1.11 km ; 591 + 572 + 37 = 1200 kg/m ; edge pull 2.26 kN/m ; sideways 22.6 MN ; arm below 32 MN  -> reproduced")
P("    In v2 the width W enters three ways:")
P("      mass:  panels = sA x W x La/L, so W is the width of the SHEET (arc length across the tent).")
P("      pull:  H = q W/2 with q = sA x g at the array, the weight of half the sheet for a 45 deg ridge, so W is again the width of the SHEET.")
P("      drop:  printed as W/4, which is the parabola's drop for a horizontal SPAN of W.")
P("    Section R then takes the footprint as W / 1.1478 (parabola with slope 0 to 1), so there W is the SHEET.  Section R's power is p_tent x W x La, W the SHEET.")

Wd = d5['W']; s_half = Wd/2
fa45, _ = parabola(45)
P("\n1c. Three readings of the design case, sheet %.1f m wide where the reading keeps it (0.5 kg/m2, ridge at 45 deg)" % Wd)
P("  reading                                             | sheet (arc) | footprint (span) | edge below the ridge | pull at each edge")
rows = []
rows.append(("A. W is the span, parabola (the drop of 12.3, Fig 12.3)", fa45*Wd, Wd, Wd/4, 0.5*g_a*fa45*Wd/2))
rows.append(("B. W is the sheet, parabola (footprint of 12.8)", Wd, Wd/fa45, Wd/fa45/4, d5['H']))
hs45, dr45 = catenary(45)
rows.append(("C. W is the sheet, catenary (the hanging sheet)", Wd, 2*hs45*s_half, dr45*s_half, d5['H']))
for nm, arc, span, drop, H in rows:
    P("  %-52s| %7.1f m   | %7.1f m        | %6.1f m             | %.3f kN/m" % (nm, arc, span, drop, H/1e3))
P("  The pull H = (weight of half the sheet) / tan(ridge slope) holds for any shape, because the edge is flat and the ridge carries the whole weight of the half.")
P("  Reading A would need a sheet %.0f m wide, with %.1f %% more mass and pull than the budget has.  The scripts' mass, pull and power are reading B or C; only the drop is reading A." % (fa45*Wd, (fa45 - 1)*100))
P("  drop as a share: catenary %.4f of the sheet width, %.4f of the footprint (the book says a quarter of the width)" % (dr45/2, dr45/(2*hs45)))
OUT["tent_design"] = dict(sheet_m=Wd, footprint_catenary_m=2*hs45*s_half, drop_catenary_m=dr45*s_half, footprint_script_m=Wd/fa45, drop_script_m=Wd/4,
                          drop_parabola_arc_m=Wd/fa45/4, edge_pull_kN_m=d5['H']/1e3)

# ---- Table 11.1 (ch15_arm_choices.py), with the catenary's drop and footprint ----
P("\n1d. Table 11.1 again (ch15_arm_choices.py, function case: W = 1,110 m, S = 10 km, arm at 45 deg and 12.5 GPa, 65.5 t per MN, retention line at 18.75 GPa)")
P("    The pull, the arm and the retention line are the script's and do not change.  The last columns are the edge below the ridge: script (W/4 x tan theta), and the catenary; then the footprint.")
ARM_T_PER_MN = 65.5; W11 = 1110.0
def case11(name, areal, theta, bow_frac):
    H = (areal*g_a*W11/2)/math.tan(math.radians(theta)); pull = H*S_SP; T = pull/math.sin(math.radians(45)); m = T/1e6*ARM_T_PER_MN
    f = bow_frac*W11; TR = H*S_SP**2/(8*f); mR = 2*RHO*TR/S75
    hs, dr = catenary(theta); drop_script = W11/4*math.tan(math.radians(theta))
    P("  %-34s| pull %5.1f MN | arm %2.0f mm, %5.0f t | retention %3.0f MN, %3.0f mm, %2.0f kg/m | edge: script %3.0f m, catenary %5.1f m | footprint %6.1f m | sheet/footprint %.4f"
      % (name, pull/1e6, dia(T, S_OP), m, TR/1e6, dia(TR, S75), mR, drop_script, dr*W11/2, 2*hs*W11/2, 1/hs))
    return dict(pull_MN=pull/1e6, arm_mm=dia(T, S_OP), arm_t=m, TR_MN=TR/1e6, ret_mm=dia(TR, S75), ret_kgm=mR, drop_script_m=drop_script, drop_catenary_m=dr*W11/2, footprint_m=hs*W11)
OUT["table_11_1"] = {}
for nm, a, th, b in (("Design case (0.5, 45 deg, W/8)", 0.5, 45, 1/8), ("Lighter panels (0.2, 45 deg, W/8)", 0.2, 45, 1/8), ("Steeper sheet (0.5, 60 deg, W/8)", 0.5, 60, 1/8),
                     ("Deeper bow (0.5, 45 deg, W/4)", 0.5, 45, 1/4), ("All three (0.2, 60 deg, W/4)", 0.2, 60, 1/4), ("not in the book: 0.5, 30 deg, W/8", 0.5, 30, 1/8)):
    OUT["table_11_1"][nm] = case11(nm, a, th, b)

# ---- Section 12.3, the lean of the arms: 30, 45, 60 deg, the strip rebalanced inside 1,200 kg/m, ridge at 45 deg ----
P("\n1e. Section 12.3, the lean of the arm (v2 with lam = 30, 45, 60 deg, the strip rebalanced inside 1,200 kg/m, ridge slope 45 deg in all three)")
OUT["lean"] = {}
for lam in (30, 45, 60):
    d = v2(0.5, lam=lam)
    P("  lean %2d deg: sheet %.3f km ; sideways pull %.1f MN ; casing load %.0f t per spreader ; edge below the ridge: script %.0f m, catenary %.0f m ; footprint: script %.3f km, catenary %.3f km"
      % (lam, d['W']/1e3, d['FR']/1e6, 2*d['Tup']/g_r/1e3, d['W']/4, dr45*d['W']/2, d['W']/fa45/1e3, hs45*d['W']/1e3))
    OUT["lean"][lam] = dict(sheet_km=d['W']/1e3, drop_script_m=d['W']/4, drop_catenary_m=dr45*d['W']/2, footprint_script_km=d['W']/fa45/1e3, footprint_catenary_km=hs45*d['W']/1e3, casing_t=2*d['Tup']/g_r/1e3)
P("  (book: 0.82 km and 7,300 t at 30 deg ; 1.11 km and 5,700 t at 45 deg ; 1.40 km and 4,200 t at 60 deg ; these come from spreader_cases.py, which used the tapered line fraction 0.0767)")

# ---- Section R: the light the tent loses, the power, the shadow ----
P("\n1f. What the shape does to the light, the power and the shadow (ch16_numbers.py section R, with the sheet/footprint ratio as the only change)")
prof = phys.compute_ring_power_profile(1.0); avg = prof['p_avg_per_m2']
t = np.linspace(0, 1, 20001); dec = np.arcsin(math.sin(math.radians(23.44))*np.sin(2*math.pi*t)); cosd = float(np.cos(dec).mean())
phi = np.linspace(-math.pi, math.pi, 1441); eta = cfg.CELL_EFFICIENCY*cfg.AR_BOOST
top = np.array([phys.direct_solar_flux(abs(x)) for x in phi])*eta; bot = np.array([phys.albedo_flux(abs(x)) for x in phi])*eta
dir_frac = top.mean()/(top.mean() + bot.mean()); f_season = dir_frac*cosd + (1 - dir_frac)
sun = 0.533; n_sp = L/S_SP; n_st = L/50e3
def days(dmax): return 365.25*2*math.asin(min(1, dmax/23.44))/math.pi
def anchor(P_total, loss, sig=5.8e7, V=1.0e6):
    """ch16_numbers.py section F, function anchor, unchanged: anchor conductors, kg per m of ring."""
    Pst = P_total/n_st; I = Pst/V; Rl = loss*V/I; A = 2*h_r/(sig*Rl); m = 2*RHO*A*h_r
    return m*n_st/L
P("  flat panel %.2f W/m2 ; season factor %.4f (mean cos of the declination %.4f, direct share %.3f)" % (avg, f_season, cosd, dir_frac))
OUT["power"] = {}
for nm, f_arc in (("parabola (script)", fa45), ("catenary", 1/hs45)):
    p_tent = avg/f_arc*f_season
    P("  -- %s: sheet / footprint %.4f ; a square meter of panel gives %.1f %% of flat ; with the season %.4f ; p_tent %.1f W/m2" % (nm, f_arc, 100/f_arc, f_season/f_arc, p_tent))
    OUT["power"][nm] = dict(f_arc=f_arc, pct_of_flat=100/f_arc, p_tent=p_tent, rows={})
    for sA in (1.0, 0.5, 0.2):
        d = v2(sA); Pw = p_tent*d['W']*La; Wfp = d['W']/f_arc
        ang = math.degrees(Wfp/h_a); x = min(1.0, ang/sun); frac = (2/math.pi)*(math.asin(x) + x*math.sqrt(1 - x*x))
        row = dict(sheet_km=d['W']/1e3, footprint_km=Wfp/1e3, drop_m=(d['W']/4 if nm.startswith("par") else dr45*d['W']/2), P_TW=Pw/1e12, GW_per_spreader=Pw/1e9/n_sp,
                   kgm_per_TW=1200/(Pw/1e12), m_per_TW=d['W']/(Pw/1e12), MW_per_km=Pw/1e6/(La/1e3), anchor_1MV=anchor(Pw, 0.02), anchor_2MV=anchor(Pw, 0.02, V=2e6),
                   md_pct=416e9/Pw*100, ang_deg=ang, cover_pct=frac*100, fade_km=(Wfp + h_a*math.radians(sun))/2e3, band_days=days(math.degrees(math.atan(Wfp/2/(h_r - h_a)))),
                   intercepted=Wfp*2*(RE + h_a)/(math.pi*RE**2))
        OUT["power"][nm]["rows"][sA] = row
        P("    %.1f kg/m2: sheet %.3f km, footprint %.3f km ; power %.2f TW ; %.2f GW per spreader ; %.1f kg/m per TW ; %.1f m of sheet per TW ; %.0f MW per km of strip ; anchor lines all by wire 1 MV %.0f, 2 MV %.0f kg/m ; design launch %.1f %%"
          % (sA, row['sheet_km'], row['footprint_km'], row['P_TW'], row['GW_per_spreader'], row['kgm_per_TW'], row['m_per_TW'], row['MW_per_km'], row['anchor_1MV'], row['anchor_2MV'], row['md_pct']))
        P("        shadow: footprint %.3f deg against the Sun's %.3f ; covers at most %.1f %% ; fades at %.2f km ; ring-panel band %.1f days a year ; share of the Earth's sunlight intercepted %.2e"
          % (row['ang_deg'], sun, row['cover_pct'], row['fade_km'], row['band_days'], row['intercepted']))
P("  logged by ch16_numbers.py section R: 87.1 % of flat, p_tent 248 ; 1.0: 0.54 km, 0.47 km, 5.8 TW, 1.31 GW, 205, anchor 377/94, 7.1 %, 11 %, 3.0 km ; 0.5: 1.11, 0.97, 12.1 TW, 2.71 GW, 99, 275 MW/km, 780/195, 3.4 %, 22 %, 3.3 km ; 0.2: 2.83, 2.47, 30.8 TW, 6.92 GW, 39")

# =====================================================================================================================
# PART 2.  THE PULL AT THE TOP OF THE ARRAY SUPPORT LINES
# =====================================================================================================================
P("\n" + "="*110 + "\nPART 2.  The pull at the top of the array support lines\n" + "="*110)
OUT["support"] = {}
for sA in (1.0, 0.5, 0.2):
    d = v2(sA); m_pan = sA*d['W']*La/L                        # panels, kg per m of ring
    x = RHO*dP_arr/S75
    A_line = m_pan*g_a/(S75 - RHO*dP_arr)                     # m2 of line per m of ring: S75 A = m_pan g_a + rho A dPhi
    m_line = RHO*A_line*(h_r - h_a)
    T_top = S75*A_line                                         # N per m of ring at the top of the lines
    # the same by adding up the weight of the line piece by piece
    n = 20000; zz = h_a + (np.arange(n) + 0.5)*(h_r - h_a)/n; w_line = float((RHO*A_line*np.array([geff(z) for z in zz])).sum()*(h_r - h_a)/n)
    T_sum = m_pan*g_a + w_line
    script = d['pan']*g_r
    OUT["support"][sA] = dict(panels_kgm=m_pan, lines_kgm=m_line, total_kgm=m_pan + m_line, factor=1/(1 - x), T_top_kN_m=T_top/1e3, T_sum_kN_m=T_sum/1e3, script_kN_m=script/1e3,
                              equiv_kgm_at_ring=T_top/g_r, per_km_MN=T_top*1e3/1e6, d_mm=dia(T_top*1e3, S75), per_50km_MN=T_top*50e3/1e6, hanger_200m_MN=T_top*200/1e6, hanger_200m_mm=dia(T_top*200, S75),
                              panels_only_kN_m=m_pan*g_a/1e3, hanger_panels_200m_MN=m_pan*g_a*200/1e6, hanger_panels_200m_mm=dia(m_pan*g_a*200, S75), share_of_casing=T_top/(R.M_LOAD_M*g_r))
    P("  %.1f kg/m2: panels %.1f kg/m of ring + lines %.1f = %.1f kg/m (line fraction %.4f)" % (sA, m_pan, m_line, m_pan + m_line, m_line/m_pan))
    P("     script (section Q): (panels + lines) x g at the ring = %.1f x %.4f = %.3f kN/m" % (d['pan'], g_r, script/1e3))
    P("     top of the lines:   panels x g at the array / (1 - rho dPhi / sigma) = %.1f x %.4f x %.4f = %.3f kN/m ; by adding up the line's weight %.3f kN/m ; that is %.1f kg/m at the ring's %.3f m/s2"
      % (m_pan, g_a, 1/(1 - x), T_top/1e3, T_sum/1e3, T_top/g_r, g_r))
    P("     a line every km %.2f MN, %.1f mm at 18.75 GPa ; on 50 km %.1f MN ; share of the casing's %.1f kN/m: %.2f %% ; hanger every 200 m: %.2f MN, %.1f mm (panels alone %.2f MN, %.1f mm)"
      % (T_top*1e3/1e6, dia(T_top*1e3, S75), T_top*50e3/1e6, R.M_LOAD_M*g_r/1e3, T_top/(R.M_LOAD_M*g_r)*100, T_top*200/1e6, dia(T_top*200, S75), m_pan*g_a*200/1e6, dia(m_pan*g_a*200, S75)))
s5 = OUT["support"][0.5]
F_line = 83.29e6/n_st
P("  rho dPhi / sigma = %.5f ; dPhi (array to ring) = %.4e J/kg ; mean g over the line %.4f m/s2" % (RHO*dP_arr/S75, dP_arr, dP_arr/(h_r - h_a)))
P("  Section 11.6, the lean of the hangers: %.1f kN / %.1f MN = %.2e of their length ; hangers of 500 m lean %.1f cm   (book: 93.6 kN / 235 MN = 4 x 10^-4, 20 cm)"
  % (F_line/1e3, s5['per_50km_MN'], F_line/(s5['per_50km_MN']*1e6), F_line/(s5['per_50km_MN']*1e6)*500*100))
# the budget counted by weight on the casing
def v2w(sA, M=1200.0, S=S_SP, lam=45.0):
    q = sA*g_a; TR = q*S*S/2; mR = 2*RHO*TR/S75
    per_w = sA*La/L*g_a/(1 - RHO*dP_arr/S75)/g_r + q/(g_r*math.tan(math.radians(lam)))
    Wd = (M - mR)/per_w
    return Wd, sA*La/L*g_a/(1 - RHO*dP_arr/S75)/g_r*Wd, q*Wd/(g_r*math.tan(math.radians(lam))), mR
Ww, panw, armw, mRw = v2w(0.5)
P("  The budget: the arms are counted as their pull on the casing divided by g at the ring (572 kg/m), the panels and lines as their true mass (591 kg/m).")
P("  Counted the same way as the arms, the panels and lines are %.1f kg/m and the design case is %.1f kg/m, not 1,200." % (s5['equiv_kgm_at_ring'], s5['equiv_kgm_at_ring'] + d5['arms'] + d5['mR']))
P("  If the 1,200 kg/m is held by weight on the casing instead: sheet %.1f m (%.2f %% narrower) ; panels and lines %.1f + arms %.1f + retention %.1f kg/m of weight at the ring" % (Ww, (1 - Ww/d5['W'])*100, panw, armw, mRw))
OUT["support"]["budget_by_weight"] = dict(sheet_m=Ww, pan=panw, arms=armw, ret=mRw, total_as_is=s5['equiv_kgm_at_ring'] + d5['arms'] + d5['mR'])

# =====================================================================================================================
# PART 3.  THE REACTION LINES AS HANGING LINES
# =====================================================================================================================
P("\n" + "="*110 + "\nPART 3.  The reaction lines as hanging lines\n" + "="*110)
P("""
Equations.  The line lies in the plane of the equator and turns with the Earth, so the only body force on it is radial:
    g(r) = GM / r^2 - w^2 r                      (gravity less the centrifugal term), with potential Phi(r), dPhi/dr = g.
Along the line, s is the length, T the tension, A the cross-section, beta the angle between the line and the local vertical.
    along the line:    dT/ds = rho A g cos(beta) = rho A g dr/ds
    constant stress:   T = sigma A                ->  T(r) = T_ground exp( rho (Phi(r) - Phi(R_E)) / sigma )
    across the line:   the force is central, so the moment about the Earth's center is the same all along the line:
                       T r sin(beta) = K           ->  sin(beta(r)) = K / (r T(r))
    shape:             d(phi)/dr = tan(beta) / r ,  ds/dr = 1 / cos(beta) ,  mass = sum of rho (T/sigma) ds
The force along the ring that the line delivers at the casing is F = T_top sin(beta_top) = K / r_top, so K = F r_top.
On a flat Earth K/r would be the same top and bottom (the horizontal pull is the same all along the line).  On the round Earth the
horizontal pull at the ground is F r_top / R_E, 11 percent more than at the casing.
The taper T_top / T_ground = exp(rho dPhi / sigma) depends on the height only, whatever the path.""")

dP_g = pot(0, h_r); TAPER = math.exp(RHO*dP_g/S_OP)
F_launch = 83.29e6                       # N, peak reaction of the design launch (ch15_numbers.py)
F1 = F_launch/n_st                       # N per reaction line, Equation 11.4
P("ground to ring: dPhi %.4e J/kg ; area taper %.4f ; lines %.2f at 50 km ; force along the ring per line %.2f kN" % (dP_g, TAPER, n_st, F1/1e3))

# ---- the book's row first: ch15_numbers.py line() ----
def line_script(F_tip, theta_deg=0.0, n=20000):
    """ch15_numbers.py, function line, unchanged: constant-stress straight line; the tip tension is put at the ground."""
    A_tip = F_tip/S_OP; m = 0.0; c = math.cos(math.radians(theta_deg))
    for i in range(n):
        z = (i + 0.5)*h_r/n
        m += RHO*A_tip*math.exp(RHO*pot(0, z)/S_OP)*(h_r/n)/c
    return A_tip, A_tip*TAPER, m, S_OP*A_tip*TAPER
Tc = F1/math.sin(math.radians(30)); A_tip, A_top, m_b, Tt_b = line_script(Tc, 30.0)
P("\n3a. The book's row (ch15_numbers.py, function line): straight at 30 deg, T = F / sin 30 = %.1f kN put at the ground ; %.2f mm at the ground, %.2f mm at the ring ; length %.1f km ; mass %.2f t ; all lines %.0f t ; %.3f kg per m of ring"
  % (Tc/1e3, dia(Tc, S_OP), dia(Tt_b, S_OP), h_r/1e3/math.cos(math.radians(30)), m_b/1e3, m_b*n_st/1e3, m_b/50e3))
P("    logged: tension 187 kN ; bottom 4.4 mm ; top 6.6 mm ; length 808 km ; mass 32.7 t ; all 29065 t  -> reproduced")
P("    What this row leaves out: with 187 kN at the ground the tension at the ring is %.1f kN.  A straight line at 30 deg would then pull along the ring with %.1f kN at the casing, %.2f times the %.1f kN asked for."
  % (Tt_b/1e3, Tt_b*0.5/1e3, Tt_b*0.5/F1, F1/1e3))
P("    A hanging line cannot do that: its horizontal pull is (nearly) the same top and bottom, so it steepens toward the ring.")
OUT["reaction"] = dict(book=dict(T_ground_kN=Tc/1e3, d_ground_mm=dia(Tc, S_OP), d_top_mm=dia(Tt_b, S_OP), length_km=h_r/1e3/math.cos(math.radians(30)), mass_t=m_b/1e3, all_t=m_b*n_st/1e3, kgm=m_b/50e3))

NQ = 400000
r_mid = RE + (np.arange(NQ) + 0.5)*h_r/NQ; dr = h_r/NQ
r_edge = RE + np.arange(NQ + 1)*h_r/NQ
def potv(r): return mu*(1/RE - 1/r) - 0.5*W_E*W_E*(r*r - RE*RE)
EXPO_mid = np.exp(RHO*potv(r_mid)/S_OP); EXPO_edge = np.exp(RHO*potv(r_edge)/S_OP)

def hang(F, beta_g_deg, stress=S_OP, want_path=False):
    """Constant-stress hanging line that delivers F along the ring at the casing and leaves the ground at beta_g from the vertical."""
    if stress != S_OP: em = np.exp(RHO*potv(r_mid)/stress); ee = np.exp(RHO*potv(r_edge)/stress)
    else: em, ee = EXPO_mid, EXPO_edge
    K = F*r_top; T_g = K/(RE*math.sin(math.radians(beta_g_deg)))
    sb = K/(r_mid*T_g*em)
    if sb.max() >= 1.0: return None
    cb = np.sqrt(1 - sb*sb)
    dphi = (sb/cb)/r_mid*dr; ds = dr/cb; dm = RHO*(T_g*em/stress)*ds
    phi_top = float(dphi.sum()); T_top = T_g*ee[-1]; sbt = K/(r_top*T_top)
    res = dict(F=F, beta_g=beta_g_deg, T_g=T_g, T_top=T_top, beta_top=math.degrees(math.asin(sbt)), phi_top=phi_top, ground_km=RE*phi_top/1e3,
               length_km=float(ds.sum())/1e3, mass_t=float(dm.sum())/1e3, d_g=dia(T_g, stress), d_top=dia(T_top, stress), V_top_kN=T_top*math.sqrt(1 - sbt*sbt)/1e3,
               H_ground_kN=K/RE/1e3, K=K)
    # the chord and the sag
    ph_e = np.concatenate(([0.0], np.cumsum(dphi))); x = r_edge*np.sin(ph_e); y = r_edge*np.cos(ph_e)
    x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]; cl = math.hypot(x1 - x0, y1 - y0)
    dist = ((x - x0)*(y1 - y0) - (y - y0)*(x1 - x0))/cl          # perpendicular distance from the chord
    i = int(np.argmax(np.abs(dist)))
    # angle of the chord from the vertical at the casing and at the ground
    ux, uy = (x1 - x0)/cl, (y1 - y0)/cl
    cas_v = (math.sin(phi_top), math.cos(phi_top)); gr_v = (0.0, 1.0)
    res.update(chord_km=cl/1e3, sag_km=abs(float(dist[i]))/1e3, sag_alt_km=(r_edge[i] - RE)/1e3, chord_from_vertical_at_casing=math.degrees(math.acos(ux*cas_v[0] + uy*cas_v[1])),
               chord_from_vertical_at_ground=math.degrees(math.acos(ux*gr_v[0] + uy*gr_v[1])), slack_km=float(ds.sum())/1e3 - cl/1e3)
    if want_path: res["path"] = (r_edge, ph_e, np.concatenate(([0.0], np.cumsum(ds))), T_g*ee)
    return res

# ---- independent check: integrate the vector equation d(T t)/ds = rho A g r_hat step by step in x, y ----
def hang_xy(F, beta_g_deg, ds=50.0):
    K = F*r_top; T_g = K/(RE*math.sin(math.radians(beta_g_deg))); b = math.radians(beta_g_deg)
    pos = np.array([0.0, RE]); Pv = T_g*np.array([math.sin(b), math.cos(b)])     # at the ground: r_hat = (0,1), along the ring = (1,0)
    def f(y):
        p = y[:2]; q = y[2:]; r = math.hypot(p[0], p[1]); T = math.hypot(q[0], q[1]); g = mu/r**2 - W_E*W_E*r
        return np.concatenate((q/T, RHO*(T/S_OP)*g*p/r))
    y = np.concatenate((pos, Pv)); s = 0.0; m = 0.0
    while True:
        k1 = f(y); k2 = f(y + 0.5*ds*k1); k3 = f(y + 0.5*ds*k2); k4 = f(y + ds*k3); yn = y + ds/6*(k1 + 2*k2 + 2*k3 + k4)
        rn = math.hypot(yn[0], yn[1]); ro = math.hypot(y[0], y[1])
        if rn >= r_top:
            fr = (r_top - ro)/(rn - ro); m += RHO*math.hypot(y[2], y[3])/S_OP*ds*fr; y = y + fr*(yn - y); s += fr*ds; break
        m += RHO*0.5*(math.hypot(y[2], y[3]) + math.hypot(yn[2], yn[3]))/S_OP*ds; y = yn; s += ds
    T = math.hypot(y[2], y[3]); r = math.hypot(y[0], y[1]); ph = math.atan2(y[0], y[1])
    rh = np.array([y[0], y[1]])/r; cosb = float(np.dot(y[2:], rh))/T
    along = float(np.dot(y[2:], np.array([rh[1], -rh[0]])))
    return dict(T_top=T, beta_top=math.degrees(math.acos(cosb)), ground_km=RE*ph/1e3, length_km=s/1e3, mass_t=m/1e3, F_along_kN=along/1e3)

def show(nm, c):
    P("  %-44s| ground %5.2f deg, %6.1f kN, %.2f mm | casing %5.2f deg, %6.1f kN, %.2f mm | anchorage %6.1f km | length %6.1f km | mass %5.2f t | sag %5.1f km at %3.0f km up | chord %4.1f deg at the casing | down-pull %5.1f kN"
      % (nm, c['beta_g'], c['T_g']/1e3, c['d_g'], c['beta_top'], c['T_top']/1e3, c['d_top'], c['ground_km'], c['length_km'], c['mass_t'], c['sag_km'], c['sag_alt_km'], c['chord_from_vertical_at_casing'], c['V_top_kN']))
def trim(c): return {k: v for k, v in c.items() if k != "path"}

P("\n3b. The hanging line, for the same force along the ring at the casing (%.2f kN) and the same 12.5 GPa everywhere" % (F1/1e3))
cases = {}
cA = hang(F1, 30.0); cases["A. 30 deg at the ground"] = cA
# B: the anchorage where the book's straight line puts it on a flat Earth, 700 km x tan 30 along the ground
gB = h_r*math.tan(math.radians(30))/1e3
bB = brentq(lambda b: hang(F1, b)['ground_km'] - gB, 15.0, 89.0); cB = hang(F1, bB); cases["B. anchorage %.0f km away (700 tan 30)" % gB] = cB
# B2: the chord at 30 deg from the casing's own vertical
bB2 = brentq(lambda b: hang(F1, b)['chord_from_vertical_at_casing'] - 30.0, 15.0, 89.0); cB2 = hang(F1, bB2); cases["B2. chord 30 deg from the casing's vertical"] = cB2
# C: the angle at the casing given
sb_max = 1.0/(TAPER*r_top/RE)
P("  The angle at the casing cannot be chosen freely: sin(beta_top) = sin(beta_ground) x R_E / (r_top x taper) <= 1 / (%.4f x %.4f), so beta_top <= %.2f deg, with the line lying flat at the ground."
  % (r_top/RE, TAPER, math.degrees(math.asin(sb_max))))
P("  A line at 30 deg at the casing does not exist at 12.5 GPa: sin(30 deg) x %.4f x %.4f = %.3f > 1." % (TAPER, r_top/RE, 0.5*TAPER*r_top/RE))
for bt in (10.0, 15.0, 20.0):
    sbg = math.sin(math.radians(bt))*TAPER*r_top/RE; cases["C. %2.0f deg at the casing" % bt] = hang(F1, math.degrees(math.asin(sbg)))
for nm, c in cases.items(): show(nm, c)
P("  %-44s| ground 30.00 deg, %6.1f kN, %.2f mm | casing 30.00 deg, %6.1f kN, %.2f mm | anchorage %6.1f km | length %6.1f km | mass %5.2f t | sag   0.0 km"
  % ("book (straight line, ch15_numbers.py)", Tc/1e3, dia(Tc, S_OP), Tt_b/1e3, dia(Tt_b, S_OP), gB, h_r/1e3/math.cos(math.radians(30)), m_b/1e3))
OUT["reaction"]["cases"] = {k: trim(v) for k, v in cases.items()}

chk = hang_xy(F1, 30.0)
P("  check of case A by stepping the vector equation (50 m steps): casing %.3f deg, %.2f kN ; anchorage %.2f km ; length %.2f km ; mass %.3f t ; force along the ring at the casing %.3f kN  (asked: %.3f)"
  % (chk['beta_top'], chk['T_top']/1e3, chk['ground_km'], chk['length_km'], chk['mass_t'], chk['F_along_kN'], F1/1e3))
OUT["reaction"]["check_xy_caseA"] = chk

# flat Earth, uniform field of the same potential drop, for the size of the round-Earth and centrifugal terms
P("\n3c. What each term does (case A, 30 deg at the ground)")
def hang_flat(F, beta_g_deg, use_centrifugal=True, n=200000):
    """Flat Earth: the horizontal pull F is the same top and bottom; T(z) = T_g exp(rho dPhi(z)/sigma)."""
    T_g = F/math.sin(math.radians(beta_g_deg)); z = (np.arange(n) + 0.5)*h_r/n; r = RE + z
    pv = mu*(1/RE - 1/r) - (0.5*W_E*W_E*(r*r - RE*RE) if use_centrifugal else 0.0); e = np.exp(RHO*pv/S_OP)
    sb = F/(T_g*e); cb = np.sqrt(1 - sb*sb); dz = h_r/n
    pt = mu*(1/RE - 1/r_top) - (0.5*W_E*W_E*(r_top**2 - RE**2) if use_centrifugal else 0.0); Tt = T_g*math.exp(RHO*pt/S_OP)
    return dict(T_g=T_g, T_top=Tt, beta_top=math.degrees(math.asin(F/Tt)), ground_km=float((sb/cb).sum()*dz)/1e3, length_km=float((1/cb).sum()*dz)/1e3, mass_t=float((RHO*T_g*e/S_OP/cb).sum()*dz)/1e3)
fl = hang_flat(F1, 30.0); fl0 = hang_flat(F1, 30.0, use_centrifugal=False)
P("  flat Earth, gravity falling with height, with the centrifugal term:    ground %.1f kN ; casing %.2f deg, %.1f kN ; anchorage %.1f km ; length %.1f km ; mass %.2f t" % (fl['T_g']/1e3, fl['beta_top'], fl['T_top']/1e3, fl['ground_km'], fl['length_km'], fl['mass_t']))
P("  flat Earth, gravity falling with height, without the centrifugal term: ground %.1f kN ; casing %.2f deg, %.1f kN ; anchorage %.1f km ; length %.1f km ; mass %.2f t" % (fl0['T_g']/1e3, fl0['beta_top'], fl0['T_top']/1e3, fl0['ground_km'], fl0['length_km'], fl0['mass_t']))
P("  round Earth (case A above):                                            ground %.1f kN ; casing %.2f deg, %.1f kN ; anchorage %.1f km ; length %.1f km ; mass %.2f t" % (cA['T_g']/1e3, cA['beta_top'], cA['T_top']/1e3, cA['ground_km'], cA['length_km'], cA['mass_t']))
OUT["reaction"]["flat"] = dict(with_centrifugal=fl, without=fl0)

# ---- every figure of Chapter 11 that quotes a reaction line, for A, B and B2 ----
P("\n3d. The figures of Chapter 11 that quote a reaction line")
c_sound = math.sqrt(E_CNT/RHO)
def figures(nm, c):
    A_g = c['T_g']/S_OP
    f = dict(T_g_kN=c['T_g']/1e3, d_g=c['d_g'], d_top=c['d_top'], length_km=c['length_km'], mass_t=c['mass_t'], kgm=c['mass_t']*1e3/50e3, all_t=c['mass_t']*n_st,
             area_ground_cm2=A_g*n_st*1e4, area_along_cm2=F_launch/S_OP*1e4, N_per_N=c['T_g']/F1, N_per_N_top=c['T_top']/F1,
             d_25km=c['d_g']/math.sqrt(2), d_100km=c['d_g']*math.sqrt(2), d_four=c['d_g']*2, sound_s=c['length_km']*1e3/c_sound, stretch_km=S_OP/E_CNT*c['length_km'],
             longer_pct=(c['length_km']*1e3/h_r - 1)*100, anchorage_km=c['ground_km'], beta_top=c['beta_top'], V_top_kN=c['V_top_kN'])
    P("  %s" % nm)
    P("     tension at the ground %.1f kN (%.2f N for each N along the ring ; at the ring %.1f kN, %.2f N for each N) ; %.2f mm at the ground, %.2f mm at the ring ; length %.1f km (%.1f %% longer than vertical) ; mass %.2f t ; %.3f kg per m of ring ; all %.0f lines %.0f t"
      % (f['T_g_kN'], f['N_per_N'], c['T_top']/1e3, f['N_per_N_top'], f['d_g'], f['d_top'], f['length_km'], f['longer_pct'], f['mass_t'], f['kgm'], n_st, f['all_t']))
    P("     combined cross-section at the ground %.1f cm2 (along the ring %.1f cm2) ; at 25 km spacing %.2f mm, at 100 km %.2f mm ; four launches %.2f mm ; sound takes %.1f s ; stretch at 2.5 %% %.1f km"
      % (f['area_ground_cm2'], f['area_along_cm2'], f['d_25km'], f['d_100km'], f['d_four'], f['sound_s'], f['stretch_km']))
    return f
OUT["reaction"]["figures"] = {}
bookc = dict(T_g=Tc, T_top=Tt_b, d_g=dia(Tc, S_OP), d_top=dia(Tt_b, S_OP), length_km=h_r/1e3/math.cos(math.radians(30)), mass_t=m_b/1e3, ground_km=gB, beta_top=30.0, V_top_kN=Tt_b*math.cos(math.radians(30))/1e3)
OUT["reaction"]["figures"]["book"] = figures("book (straight line, ch15_numbers.py)", bookc)
for nm in list(cases)[:3]: OUT["reaction"]["figures"][nm] = figures(nm, cases[nm])
P("  speed of sound in the line sqrt(E/rho) = %.2f km/s" % (c_sound/1e3))

# ---- the spring: an elastic hanging line, built for the design shape, with the casing moved along the ring ----
P("\n3e. The line as a spring (elastic hanging line; the cross-section is the one built for the design shape, E = %.0f GPa)" % (E_CNT/1e9))
P("""    Along the unstretched length s0, with stretch 1 + eps and eps = T / (E A_built(s0)):
        dr/ds0 = (1 + eps) cos(beta) ,  d(phi)/ds0 = (1 + eps) sin(beta) / r ,  dT/ds0 = m0(s0) g(r) cos(beta) ,  sin(beta) = K / (r T)
    m0 is the mass per unstretched meter (the built line's mass does not change).  The two unknowns T_ground and K are found so that the top of the
    line is at the casing's radius and at the casing's new place, x meters along the ring from the design place (x > 0 away from the anchorage).
    The force along the ring at the casing is K / r_top.  The book's spring is the straight line's, F / (stretch / sin 30) = 2.32 N/m (ch15_numbers.py).""")
def spring(c_design, xs, n=4000, say=True):
    r_e, ph_e, s_e, T_e = c_design["path"]
    lam_d = 1 + S_OP/E_CNT; s0_e = s_e/lam_d; L0 = s0_e[-1]
    Td = lambda s0: np.interp(s0, s0_e, T_e)                   # design tension at s0 ; built area = Td / S_OP ; m0 = rho (Td/S_OP) lam_d
    h = L0/n
    def shoot(Tg, K):
        y = np.array([RE, 0.0, Tg]); ok = True
        def f(s0, y):
            r, ph, T = y; td = float(Td(s0)); sb = K/(r*T)
            if abs(sb) >= 1.0: sb = math.copysign(0.999999, sb)
            cb = math.sqrt(1 - sb*sb); lam = 1 + (S_OP/E_CNT)*T/td
            return np.array([lam*cb, lam*sb/r, RHO*(td/S_OP)*lam_d*(mu/r**2 - W_E*W_E*r)*cb])
        s0 = 0.0
        for i in range(n):
            k1 = f(s0, y); k2 = f(s0 + h/2, y + h/2*k1); k3 = f(s0 + h/2, y + h/2*k2); k4 = f(s0 + h, y + h*k3)
            y = y + h/6*(k1 + 2*k2 + 2*k3 + k4); s0 += h
        return y
    Tg0 = c_design['T_g']; K0 = c_design['K']; phi0 = c_design['phi_top']
    out = []; guess = np.array([1.0, 1.0])
    y0 = shoot(Tg0, K0)
    if say: P("    built length (unstretched) %.2f km ; check at x = 0: top at %.1f m from the casing's radius and %.1f m along the ring from the design place" % (L0/1e3, y0[0] - r_top, (y0[1] - phi0)*r_top))
    for x in xs:
        def res(u):
            y = shoot(u[0]*Tg0, u[1]*K0); return [(y[0] - r_top)/1e3, (y[1] - phi0 - x/r_top)*r_top/1e3]
        sol = root(res, guess, method="hybr", tol=1e-10)
        if not sol.success or max(abs(v) for v in res(sol.x)) > 1e-4:
            out.append((x, None, None, None)); continue
        Tg = sol.x[0]*Tg0; K = sol.x[1]*K0; sbg = K/(RE*Tg)
        y = shoot(Tg, K)
        out.append((x, K/r_top, Tg, math.degrees(math.asin(min(1.0, sbg))) if sbg < 1 else 90.0, y[2]))
        guess = sol.x
    return out, L0
OUT["reaction"]["spring"] = {}
for nm, c0, bg in (("A. 30 deg at the ground", cA, 30.0), ("B. anchorage %.0f km away" % gB, cB, bB)):
    cd = hang(F1, bg, want_path=True)
    xs_up = [0.0, 100.0, 800.0, 5e3, 10e3]; xs_dn = [-100.0, -800.0, -5e3, -10e3, -20e3, -40e3, -80e3]
    P("  %s (built to give %.1f kN along the ring at x = 0):" % (nm, F1/1e3))
    up, L0 = spring(cd, xs_up); dn, _ = spring(cd, xs_dn, say=False)
    tab = {}
    for row in sorted(up[1:] + dn + [up[0]], key=lambda t: t[0]):
        if row[1] is None: P("     x = %8.1f km: no solution found (the line would lie on the ground or the solver did not converge)" % (row[0]/1e3)); continue
        x, H, Tg, bgx, Ttop = row
        P("     x = %8.2f km: force along the ring %7.2f kN (%5.1f %% of design) ; ground %6.1f kN at %5.2f deg ; top %6.1f kN = %5.1f %% of the working stress" % (x/1e3, H/1e3, H/F1*100, Tg/1e3, bgx, Ttop/1e3, Ttop/c0['T_top']*100))
        tab[x] = dict(H_kN=H/1e3, T_g_kN=Tg/1e3, beta_g=bgx, T_top_kN=Ttop/1e3)
    k = (tab[100.0]['H_kN'] - tab[-100.0]['H_kN'])*1e3/200.0
    P("     spring at the design load: %.3f N for each meter the casing moves (book: 2.3 N/m) ; force change at the 800 m swing: +%.2f / %.2f kN" % (k, tab[800.0]['H_kN'] - F1/1e3, tab[-800.0]['H_kN'] - F1/1e3))
    OUT["reaction"]["spring"][nm] = dict(k_N_per_m=k, built_length_km=L0/1e3, table={str(kk): vv for kk, vv in tab.items()})
k_arm = 1/(1/(23e6/100e3) + 1/(32e6/893e3))
P("  The ring's hold (ch15_spreader_update.py: 2 x %.1f arms at %.1f N/m + %.1f reaction lines at 2.32 N/m = %.0f N/m):" % (n_sp, k_arm, n_st, 2*n_sp*k_arm + n_st*2.32))
for nm, v in OUT["reaction"]["spring"].items():
    kt = 2*n_sp*k_arm + n_st*v['k_N_per_m']; M_cas = 12000.0*L; om = math.sqrt(kt/M_cas)
    P("     with %s: %.0f N/m ; period %.0f s (%.2f h) ; swing after a launch %.0f m, with the sled braked %.0f m ; force change in a reaction line at that swing %.2f kN"
      % (nm, kt, 2*math.pi/om, 2*math.pi/om/3600, 0.576/om, 0.281/om, v['k_N_per_m']*0.576/om/1e3))
    v.update(k_ring=kt, period_s=2*math.pi/om, swing_m=0.576/om, swing_braked_m=0.281/om)
P("  The standing pull: a line built for the design shape pulls along the ring with its design force whenever the casing is at its design place, launch or no launch.")
P("  All %.0f lines lean the same way, so together that is %.1f MN along the ring at all times, to be held by something else or let off by building the lines longer." % (n_st, F_launch/1e6))

with open(os.path.join(LOG_DIR, "round9_lines_array.txt"), "w", encoding="utf-8") as f: f.write(_buf.getvalue())
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
with open(os.path.join(LOG_DIR, "round9_lines_array.json"), "w", encoding="utf-8") as f: json.dump(clean(OUT), f, indent=1)
print("\nlog written: runs/round9_logs/round9_lines_array.txt and .json")
