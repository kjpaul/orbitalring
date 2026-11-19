#!/usr/bin/env python3
"""
Round 24: Paul's layout of the lines to the ground (layout F), the guide line on the round Earth, the tide and the Earth's own field.
Brief: Vol_III/Claude_Code_instructions_round24_2026-11-19.md.   Book: Orbital Ring Engineering, Volume III: Chapter 11, Chapter 2 (Section 2.3), Chapter 7.

  1. Layout F: reaction lines in opposed pairs every 25 km, a vertical 100 t guide line every 100 km or 200 km.  F1: the pairs sized for the
     chapter's 5,096 N/m.  F2: the pairs sized for the chapter's line mass.  Layouts A (the chapter) and D beside them.  Every guide line at its
     round-Earth pendulum stiffness.
  2. The sideways loads: the lunar and solar tide ; the Earth's own field from EGM2008 to degree 12 ; what carries them ; the cable's bow.
  3. The conductors on the guide lines.

The lines are round 9's model through round 23's log (runs/round23_logs/round23_reaction_vs_guide.json: the design line's size and stiffness, its
pull and top tension against the casing's place, the guide line on the round Earth).  The launch, the braking and the casing's integrator are
round 20's (tools/round20_brake_felt_pods.py), taken as they stand.

Sources of the constants that are not the repository's own:
  EGM2008 coefficients: data/egm2008_to_degree_12.gfc, the head of ICGEM's EGM2008.gfc (https://icgem.gfz.de/tom_longtime, entry 104), fetched
      19 November 2026.  Fully normalized, tide free, GM 3.986004415e14 m3/s2, radius 6,378,136.3 m.
  GM of the Sun 1.32712440041279419e20 m3/s2, astronomical unit 149,597,870,700 m, obliquity 84,381.412 arcsec: JPL Solar System Dynamics,
      Astrodynamic Parameters (https://ssd.jpl.nasa.gov/astro_par.html), opened 19 November 2026.
  GM of the Moon 4.9048695e12 m3/s2 and its mean distance 384,400 km: the repository's lunar/lunar_config.py (MU_MOON, A_MOON_ORBIT).
  Moon's perigee 363,300 km and its greatest declination 28 deg 36': Wikipedia, "Orbit of the Moon" (opened 19 November 2026).
  Earth's perihelion 147.10 million km: Wikipedia, "Earth's orbit" (opened 19 November 2026).

Run from Python/orbitalring:  python tools/round24_layout_f.py          (about 1 minute)
Log: runs/round24_logs/round24_layout_f.txt and .json
"""
import contextlib, io, json, math, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
import ring_altitude as R

T_START = time.time()
LOG_DIR = os.path.join(ROOT, "runs", "round24_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
OUT = {}
L23 = json.load(open(os.path.join(ROOT, "runs", "round23_logs", "round23_reaction_vs_guide.json"), encoding="utf-8"))
MU = R.GM; RE = R.R_EARTH_EQ; W_E = R.OMEGA_SIDEREAL; RHO = R.RHO_CNT; S_OP = R.SIGMA_CABLE; h_r = R.ALTITUDE; Lr = R.L_RING; r_top = RE + h_r
def potv(r): return MU*(1/RE - 1/r) - 0.5*W_E*W_E*(r*r - RE*RE)          # round 9's potential for a body that turns with the Earth
G0 = MU/RE**2 - W_E*W_E*RE
F_LAUNCH = 83.29e6; n50 = Lr/50e3; F1 = F_LAUNCH/n50
row404 = next(r for r in L23["by_distance"] if abs(r["d_km"] - 404.145188432738) < 1e-6)
kr = row404["k"]; m_r = row404["mass_t"]; kns_r = row404["k_ns"]; dg_r = row404["d_g"]; dt_r = row404["d_top"]; Tg_r = row404["T_g_kN"]; Tt_r = row404["T_top_kN"]; Vt_r = row404["V_top_kN"]
G = L23["guide"]; k_g = G["k_p_round"]; k_g_flat = G["k_p"]; m_g = G["mass_t"]; T_g_top = G["T_top_kN"]*1e3; T_g_gnd = G["T_g_kN"]*1e3
TAB = {float(k): v for k, v in L23["both"]["table"].items()}; xs_t = np.array(sorted(TAB)); H_t = np.array([TAB[x]["H"] for x in xs_t]); Tt_t = np.array([TAB[x]["T_top"] for x in xs_t]); k_t = np.array([TAB[x]["k"] for x in xs_t])
def stress_at(x):
    """Top tension of a reaction line over its tension at the design place, with the casing x m away from its anchorage (round 23's table, interpolated)."""
    return float(np.interp(x, xs_t, Tt_t))/float(np.interp(0.0, xs_t, Tt_t))
RULE = 0.1*95.0e3*30.0

P("Round 24.  Layout F, the guide line on the round Earth, the tide and the Earth's own field.")
P("From round 23's log: the design reaction line (anchorage %.2f km, %.3f kN along the ring) %.4f N/m, %.3f t, %.3f and %.3f mm ; the 100 t guide line on the round Earth %.4f N/m (flat Earth %.4f), %.2f t."
  % (row404["d_km"], F1/1e3, kr, m_r, dg_r, dt_r, k_g, k_g_flat, m_g))

# ---------------------------------------------------------------------------------------------------------------------
# round 20's launch, braking and integrator
src20 = open(os.path.join(HERE, "round20_brake_felt_pods.py"), encoding="utf-8").read(); cut = src20.index('P("Round 20.'); assert cut > 0
g20 = {"__name__": "round20_prefix", "__file__": os.path.join(HERE, "round20_brake_felt_pods.py")}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src20[:cut], "round20_brake_felt_pods.py (prefix)", "exec"), g20)
    g20["LOG_DIR"] = LOG_DIR
    D = g20["lsm_run"]("design_brake")
M_CAS = g20["M_CAS"]; F_LIM = g20["F_LIM"]; rebin = g20["rebin"]; swing_once = g20["swing_once"]; swing_daily = g20["swing_daily"]
L20 = json.load(open(os.path.join(ROOT, "runs", "round20_logs", "round20_brake_felt_pods.json"), encoding="utf-8"))
assert abs(D["dur"] - L20["design"]["dur"]) < 1e-6 and abs(D["impulse"] - L20["design"]["impulse"]) < 1.0
Fl = D["F"]; Fb = D["FB"]; F_nb = rebin(Fl, 0.1, 10.0); F_br = rebin(np.concatenate((Fl, -Fb)), 0.1, 10.0); J_br = float(Fl.sum()*0.1 - Fb.sum()*0.1)
M_CAS_M = M_CAS/Lr

# =====================================================================================================================
P("\n" + "="*118 + "\n1. The layouts\n" + "="*118)
P("The design launch and its braking: lsm_simulation, run here, the same as round 20's log (%.1f min, %.4e N s ; braking %.1f min, %.4e N s).  Casing %.4e kg, %.0f kg/m.  LIMs %.1f MN."
  % (D["dur"]/60, D["impulse"], D["brake"]["duration_s"]/60, D["brake"]["impulse"], M_CAS, M_CAS_M, F_LIM/1e6))
K_BOOK = n50*(kr + 2*k_g_flat); KGM_BOOK = (m_r + 2*m_g)*1e3/50e3
P("The chapter's printed hold: %.1f x (%.4f + 2 x %.4f) = %.1f N/m (the guide line on a flat Earth).  The chapter's lines: %.3f kg for each meter of ring." % (n50, kr, k_g_flat, K_BOOK, KGM_BOOK))

def layout(key, name, s, react_spacing, n_react, guide_spacing, n_guide, paired):
    npt = Lr/react_spacing; ng = (Lr/guide_spacing)*n_guide if n_guide else 0.0; nl = npt*n_react
    Y = dict(key=key, name=name, s=s, react_spacing_km=react_spacing/1e3, guide_spacing_km=(guide_spacing/1e3 if n_guide else None), paired=paired, n_react_lines=nl, n_guide_lines=ng,
             F_kN=s*F1/1e3, d_g=dg_r*math.sqrt(s), d_top=dt_r*math.sqrt(s), T_g_kN=s*Tg_r, T_top_kN=s*Tt_r, V_top_kN=s*Vt_r, mass_t=s*m_r, k_line=s*kr)
    k = nl*s*kr + ng*k_g; Y["hold"] = k; Y["hold_react"] = nl*s*kr; Y["hold_guide"] = ng*k_g
    Y["period_s"] = 2*math.pi*math.sqrt(M_CAS/k); c = 2*math.sqrt(k*M_CAS); Y["c_crit"] = c
    Y["once"] = {nm: swing_once(F, 10.0, c, k=k) for nm, F in (("not braked", F_nb), ("braked", F_br))}
    Y["daily"] = {N: swing_daily(F_br, 10.0, N, c, False, k=k, J_exact=J_br) for N in (7, 14)}
    Y["all_t"] = nl*s*m_r + ng*m_g; Y["kgm"] = Y["all_t"]*1e3/Lr; Y["kgm_react"] = nl*s*m_r*1e3/Lr; Y["kgm_guide"] = ng*m_g*1e3/Lr
    Y["anchorages"] = nl + ng; Y["hold_ns"] = nl*s*kns_r + ng*k_g
    Y["line_dF_kN"] = {nm: s*kr*Y["once"][nm]["x_max"]/1e3 for nm in Y["once"]}
    if paired:
        Y["rest_m"] = 0.0; Y["rest_pull"] = 1.0; Y["standing_MN"] = 0.0
        Y["stress"] = dict(once=stress_at(Y["once"]["not braked"]["x_max"]), **{"d%d" % N: stress_at(Y["daily"][N]["x_mean"] + Y["daily"][N]["swing"]) for N in (7, 14)})
        Y["stress_slack"] = dict(once=stress_at(-Y["once"]["not braked"]["x_max"]), **{"d%d" % N: stress_at(-(Y["daily"][N]["x_mean"] + Y["daily"][N]["swing"])) for N in (7, 14)})
    else:
        # one-way lines built for the design shape, held by the guide lines as pendulums: where the casing rests
        from scipy.optimize import brentq
        gpp = ng/npt                                                         # guide lines for each reaction line
        xa = brentq(lambda x: float(np.interp(x, xs_t, H_t)) + gpp*k_g*x/s, xs_t[0], 0.0)
        Y["rest_m"] = xa; Y["rest_pull"] = float(np.interp(xa, xs_t, H_t))/F1; Y["k_at_rest"] = s*float(np.interp(xa, xs_t, k_t)); Y["hold_at_rest"] = nl*Y["k_at_rest"] + ng*k_g
        Y["standing_MN"] = nl*s*F1/1e6
        Y["stress"] = dict(once=stress_at(Y["once"]["not braked"]["x_max"]), **{"d%d" % N: stress_at(Y["daily"][N]["x_mean"] + Y["daily"][N]["swing"]) for N in (7, 14)})
    n_par = 1 if s <= 1.0 + 1e-9 else math.ceil(s - 1e-9)
    Y["parallel"] = n_par; Y["d_each"] = dg_r*math.sqrt(s/n_par)
    return Y

nP = Lr/25e3
def s_hold(gs): return (K_BOOK - (Lr/gs)*k_g)/(2*nP*kr)
def s_mass(gs): return (KGM_BOOK*Lr/1e3 - (Lr/gs)*m_g)/(2*nP*m_r)
sD = (m_r + 2*m_g)/(2*m_r)
LAY = [layout("A", "A. The chapter: one reaction line and two guide lines every 50 km", 1.0, 50e3, 1, 50e3, 2, False),
       layout("D", "D. Opposed pairs every 50 km with the mass of layout A's lines, no guide lines (round 23)", sD, 50e3, 2, None, 0, True),
       layout("F1-100", "F1, guide line every 100 km: pairs every 25 km sized for the chapter's %.0f N/m" % K_BOOK, s_hold(100e3), 25e3, 2, 100e3, 1, True),
       layout("F1-200", "F1, guide line every 200 km: pairs every 25 km sized for the chapter's %.0f N/m" % K_BOOK, s_hold(200e3), 25e3, 2, 200e3, 1, True),
       layout("F2-100", "F2, guide line every 100 km: pairs every 25 km sized for the chapter's %.2f kg/m" % KGM_BOOK, s_mass(100e3), 25e3, 2, 100e3, 1, True),
       layout("F2-200", "F2, guide line every 200 km: pairs every 25 km sized for the chapter's %.2f kg/m" % KGM_BOOK, s_mass(200e3), 25e3, 2, 200e3, 1, True)]
for Y in LAY:
    P("\n%s" % Y["name"])
    P("   reaction lines: %.1f (%s every %.0f km) ; each %.4f times the design line, built for %.2f kN along the ring ; %.3f mm at the ground, %.3f mm at the top ; %.1f kN at the ground, %.1f kN at the top ; %.3f t ; %.4f N/m"
      % (Y["n_react_lines"], "an opposed pair" if Y["paired"] else "one", Y["react_spacing_km"], Y["s"], Y["F_kN"], Y["d_g"], Y["d_top"], Y["T_g_kN"], Y["T_top_kN"], Y["mass_t"], Y["k_line"]))
    P("   by the chapter's rule (the design line is %.2f mm ; four launches at once are %.1f mm, carried as several lines): %s" % (dg_r, 2*dg_r, "one line, thinner than the chapter's" if Y["s"] < 1 - 1e-9 else ("one line" if Y["parallel"] == 1 else "%d lines of %.2f mm in parallel" % (Y["parallel"], Y["d_each"]))))
    P("   guide lines: %.1f%s, each %.4f N/m on the round Earth, %.2f t" % (Y["n_guide_lines"], " (one every %.0f km)" % Y["guide_spacing_km"] if Y["key"].startswith("F") else (" (two every 50 km)" if Y["n_guide_lines"] else ""), k_g, m_g))
    P("   hold along the ring %.1f N/m (reaction lines %.1f, guide lines %.1f) ; period %.0f s, %.2f h ; critical damping %.3e N s/m" % (Y["hold"], Y["hold_react"], Y["hold_guide"], Y["period_s"], Y["period_s"]/3600, Y["c_crit"]))
    for nm in ("not braked", "braked"):
        o = Y["once"][nm]
        P("   one design launch, LIMs damping, sled %-10s: farthest %.0f m, %.2f h after the launch begins ; within 1 %% after %.1f h ; LIMs' peak push %.2f MN ; a reaction line's pull changes by %.2f kN"
          % (nm, o["x_max"], o["peak_h"], o["die_h"], o["F_lim_MN"], Y["line_dF_kN"][nm]))
    for N in (7, 14):
        o = Y["daily"][N]
        P("   %2d braked design launches a day (mean push %.3f MN): the casing sits at %.0f m and swings %.0f m about it ; LIMs' peak push %.1f MN" % (N, o["mean_push_MN"], o["x_mean"], o["swing"], o["F_lim_peak_MN"]))
    P("   lines for each meter of ring %.3f kg/m (reaction lines %.3f, guide lines %.3f) ; all %.0f t ; anchorages on the ground %.1f" % (Y["kgm"], Y["kgm_react"], Y["kgm_guide"], Y["all_t"], Y["anchorages"]))
    if Y["paired"]:
        P("   with no traffic the casing rests at the design place: the pairs' standing pulls cancel (0 m, by symmetry) ; each line stands at its design pull")
        P("   top tension of the tighter line, of the working stress: %.2f %% at the swing of one launch ; %.2f %% at 7 a day ; %.2f %% at 14 a day (the slacker line %.2f, %.2f, %.2f %%)"
          % (100*Y["stress"]["once"], 100*Y["stress"]["d7"], 100*Y["stress"]["d14"], 100*Y["stress_slack"]["once"], 100*Y["stress_slack"]["d7"], 100*Y["stress_slack"]["d14"]))
    else:
        P("   with no traffic, the lines built for the design shape: the casing rests %.1f km toward the anchorages, the reaction lines at %.1f %% of their design pull ; hold there %.0f N/m.  Standing pull at the design place %.1f MN."
          % (-Y["rest_m"]/1e3, 100*Y["rest_pull"], Y["hold_at_rest"], Y["standing_MN"]))
        P("   top tension of a reaction line counted from the design place, of the working stress: %.2f %% at the swing of one launch ; %.2f %% at 7 a day ; %.2f %% at 14 a day" % (100*Y["stress"]["once"], 100*Y["stress"]["d7"], 100*Y["stress"]["d14"]))
    P("   at the top: a reaction line %.1f kN, %.2f times the rule's %.0f kN (it hangs from the force spreader cable)%s" % (Y["T_top_kN"], Y["T_top_kN"]*1e3/RULE, RULE/1e3, " ; a guide line %.2f MN, eight branches on the casing" % (T_g_top/1e6) if Y["n_guide_lines"] else ""))
    P("   north-south pendulum hold of the lines [estimate for the reaction lines, round 23]: %.0f N/m for each meter the casing moves sideways (guide lines %.0f)" % (Y["hold_ns"], Y["hold_guide"]))
OUT["layouts"] = LAY; OUT["K_BOOK"] = K_BOOK; OUT["KGM_BOOK"] = KGM_BOOK
P("  (%.0f s so far)" % (time.time() - T_START))

# =====================================================================================================================
P("\n" + "="*118 + "\n2. The sideways loads\n" + "="*118)
cfgtxt = open(os.path.join(ROOT, "lunar", "lunar_config.py"), encoding="utf-8").read()
def cfgval(nm): return float(eval(re.search(r"^%s\s*=\s*([0-9_.eE+\-]+)" % nm, cfgtxt, re.M).group(1).replace("_", "")))
GM_MOON = cfgval("MU_MOON"); D_MOON = cfgval("A_MOON_ORBIT")
GM_SUN = 1.32712440041279419e20; AU = 149597870700.0; OBLIQ = 84381.412/3600.0
D_MOON_PERI = 363.3e6; DEC_MOON = 28 + 36/60.0; D_SUN_PERI = 147.10e9
M_CABLE_M = R.M_CABLE_M; T_CABLE = R.CABLE_TENSION
P("(a) The tide.  A body of gravitational parameter GM at the distance d and the declination delta pulls a point of the ring, at r from the Earth's center, out of the")
P("    plane of the equator with the acceleration (3 GM r / d^3) sin(delta) cos(delta) cos(H), H the body's hour angle at that point.  The largest value is 1.5 (GM / d^3) r sin(2 delta).")
P("    r = %.0f km.  Casing %.0f kg/m ; cable %.0f kg/m (ring_altitude.M_CABLE_M)." % (r_top/1e3, M_CAS_M, M_CABLE_M))
tide = {}
P("    case                                              | acceleration   | on the casing | on the cable and casing")
for nm, GMb, d, dec in (("Moon, mean distance %.0f km, %.2f deg" % (D_MOON/1e3, DEC_MOON), GM_MOON, D_MOON, DEC_MOON), ("Moon, perigee %.0f km, %.2f deg" % (D_MOON_PERI/1e3, DEC_MOON), GM_MOON, D_MOON_PERI, DEC_MOON),
                        ("Sun, 1 au, %.3f deg" % OBLIQ, GM_SUN, AU, OBLIQ), ("Sun, perihelion %.2f million km, %.3f deg" % (D_SUN_PERI/1e9, OBLIQ), GM_SUN, D_SUN_PERI, OBLIQ)):
    a = 1.5*(GMb/d**3)*r_top*math.sin(math.radians(2*dec)); tide[nm] = dict(a=a, casing=a*M_CAS_M, both=a*(M_CAS_M + M_CABLE_M))
    P("    %-49s | %.4e m/s2 | %.5f N/m   | %.5f N/m" % (nm, a, a*M_CAS_M, a*(M_CAS_M + M_CABLE_M)))
ks = list(tide); a_mean = tide[ks[0]]["a"] + tide[ks[2]]["a"]; a_worst = tide[ks[1]]["a"] + tide[ks[3]]["a"]
P("    Moon and Sun together, both at their greatest declination on the same side and the same hour angle: mean distances %.4e m/s2, %.5f N/m on the casing, %.5f N/m on the cable and casing ;"
  % (a_mean, a_mean*M_CAS_M, a_mean*(M_CAS_M + M_CABLE_M)))
P("    at perigee and perihelion %.4e m/s2, %.5f N/m on the casing, %.5f N/m on the cable and casing." % (a_worst, a_worst*M_CAS_M, a_worst*(M_CAS_M + M_CABLE_M)))
P("    The old chapter script (spreader_loads_check.py) has the same two accelerations at the mean distances (7.70e-07 and 3.07e-07 m/s2) and puts them on 22,686 kg/m, 0.0245 N/m.  That mass is the cable's hardware (%.0f kg/m) and the casing,"
  % R.M_HARDWARE_M)
P("    not the cable: the script looks for a name ring_altitude does not have and falls back.  On the cable and casing as they are (%.0f kg/m) the same accelerations give %.4f N/m." % (M_CAS_M + M_CABLE_M, a_mean*(M_CAS_M + M_CABLE_M)))
P("    The book's 'on the order of 0.1 N/m' (Chapter 2, Section 2.3) is %.1f times the worst case on the cable and casing together, and %.1f times the worst case on the casing." % (0.1/(a_worst*(M_CAS_M + M_CABLE_M)), 0.1/(a_worst*M_CAS_M)))
P("    The tide's force goes as cos(H) around the ring: one way on the side toward the body and on the side away from it, the other way a quarter of the ring from them.")
OUT["tide"] = dict(cases=tide, a_mean=a_mean, a_worst=a_worst, casing_mean=a_mean*M_CAS_M, casing_worst=a_worst*M_CAS_M, both_mean=a_mean*(M_CAS_M + M_CABLE_M), both_worst=a_worst*(M_CAS_M + M_CABLE_M))

# ---- (b) the Earth's own field
P("\n(b) The Earth's own field, EGM2008 (data/egm2008_to_degree_12.gfc).")
C = {}; Sx = {}; GM_E = A_E = None
for l in open(os.path.join(ROOT, "data", "egm2008_to_degree_12.gfc"), encoding="utf-8"):
    p = l.split()
    if not p: continue
    if p[0] == "earth_gravity_constant": GM_E = float(p[1])
    elif p[0] == "radius": A_E = float(p[1])
    elif p[0] == "gfc":
        n, m = int(p[1]), int(p[2]); C[(n, m)] = float(p[3].replace("d", "e")); Sx[(n, m)] = float(p[4].replace("d", "e"))
NMAX_FILE = max(n for n, m in C)
P("    Model constants: GM %.10e m3/s2, radius %.1f m ; fully normalized, tide free ; degrees 2 to %d read.  C20 %.15e, C30 %.15e, C22 %.15e, S22 %.15e." % (GM_E, A_E, NMAX_FILE, C[(2, 0)], C[(3, 0)], C[(2, 2)], Sx[(2, 2)]))
def pbar(nmax, x):
    """Fully normalized associated Legendre functions P[n][m](x), x = sin(latitude), by the standard forward recursions."""
    Pm = [[0.0]*(nmax + 2) for _ in range(nmax + 2)]; u = math.sqrt(1 - x*x)
    Pm[0][0] = 1.0; Pm[1][1] = math.sqrt(3.0)*u
    for n in range(2, nmax + 1): Pm[n][n] = u*math.sqrt((2*n + 1)/(2.0*n))*Pm[n - 1][n - 1]
    for m in range(0, nmax + 1):
        for n in range(m + 1, nmax + 1):
            a = math.sqrt((4.0*n*n - 1)/(n*n - m*m)); b = math.sqrt(((n - 1.0)**2 - m*m)/(4.0*(n - 1)**2 - 1)) if n - 1 > m else 0.0
            Pm[n][m] = a*(x*Pm[n - 1][m] - b*(Pm[n - 2][m] if n >= 2 else 0.0))
    return Pm
def field(r, lam, nmax, nmin=2, dphi=1e-6):
    """Potential terms of degrees nmin..nmax at latitude 0: returns g_north = (1/r) dU/dphi, g_east = (1/r) dU/dlambda, g_down = -dU/dr (arrays over lam), m/s2."""
    Pp = pbar(nmax, math.sin(dphi)); Pn = pbar(nmax, math.sin(-dphi)); P0 = pbar(nmax, 0.0)
    gN = np.zeros_like(lam); gE = np.zeros_like(lam); gD = np.zeros_like(lam)
    for n in range(nmin, nmax + 1):
        f = (GM_E/r)*(A_E/r)**n
        for m in range(0, n + 1):
            cs = C[(n, m)]*np.cos(m*lam) + Sx[(n, m)]*np.sin(m*lam); dcs = m*(-C[(n, m)]*np.sin(m*lam) + Sx[(n, m)]*np.cos(m*lam))
            gN += f*(Pp[n][m] - Pn[n][m])/(2*dphi)*cs/r; gE += f*P0[n][m]*dcs/r; gD += (n + 1)*f/r*P0[n][m]*cs
    return gN, gE, gD
lam = np.radians(np.arange(0.0, 360.0, 0.1))
_, _, gd = field(A_E, lam, NMAX_FILE); g_surf = GM_E/A_E**2 + gd.mean()
P("    Check of the sum: gravity at the model's radius on the equator, mean around it, to degree %d: %.5f m/s2 (the accepted equatorial gravitation, 9.7803 + 0.0339 of rotation, is 9.8142)." % (NMAX_FILE, g_surf))
assert abs(g_surf - 9.8142) < 2e-3
P("    At r = %.0f km on the equator.  North is positive ; east is positive ; longitude east of Greenwich.  Force on the casing = %.0f kg/m x the acceleration." % (r_top/1e3, M_CAS_M))
P("    to degree | north-south: mean around the ring | largest (where)              | smallest (where)             | along the ring: largest (where)      | rms")
earth = {}
for nmax in (2, 3, 4, 6, 8, 10, 12):
    gN, gE, _ = field(r_top, lam, nmax); fN = gN*M_CAS_M; fE = gE*M_CAS_M
    i1 = int(np.argmax(fN)); i2 = int(np.argmin(fN)); j = int(np.argmax(np.abs(fE)))
    earth[nmax] = dict(gN_mean=float(gN.mean()), fN_mean=float(fN.mean()), fN_max=float(fN[i1]), lon_max=float(np.degrees(lam[i1])), fN_min=float(fN[i2]), lon_min=float(np.degrees(lam[i2])), fN_absmax=float(np.abs(fN).max()),
                       fE_max=float(fE[j]), lon_E=float(np.degrees(lam[j])), fE_rms=float(np.sqrt((fE**2).mean())), gN_absmax=float(np.abs(gN).max()), gE_absmax=float(np.abs(gE).max()))
    e = earth[nmax]
    P("    %9d | %+.5f N/m (%+.4e m/s2)     | %+.5f N/m (%5.1f deg E)     | %+.5f N/m (%5.1f deg E)     | %+.5f N/m (%5.1f deg E)             | %.5f N/m" % (nmax, e["fN_mean"], e["gN_mean"], e["fN_max"], e["lon_max"], e["fN_min"], e["lon_min"], e["fE_max"], e["lon_E"], e["fE_rms"]))
gN, gE, _ = field(r_top, lam, NMAX_FILE); fN = gN*M_CAS_M; fE = gE*M_CAS_M
P("    The parts of the north-south force at degree %d: the same all the way around (the zonal terms of odd degree, J3 first) %+.5f N/m ; the part that changes with the longitude, largest %.5f N/m, rms %.5f N/m."
  % (NMAX_FILE, fN.mean(), float(np.abs(fN - fN.mean()).max()), float(np.sqrt(((fN - fN.mean())**2).mean()))))
P("    Every 30 degrees of longitude (degree %d): north-south N/m, along the ring N/m" % NMAX_FILE)
P("      " + " ; ".join("%3.0f E: %+.3f, %+.3f" % (np.degrees(lam[i]), fN[i], fE[i]) for i in range(0, len(lam), 300)))
ds = Lr/len(lam); cum = np.cumsum(fE - fE.mean())*ds
P("    Along the ring the force adds to nothing around the ring (mean %+.2e N/m).  Carried by the casing along its own length it is a tension that swings by %.2f MN, peak to peak, around the ring." % (fE.mean(), (cum.max() - cum.min())/1e6))
z0 = fN.mean()/M_CAS_M*r_top**3/MU
P("    Not converged: the terms fall only as (%.4f)^n at this height, and the sum is given to degree %d because that is as far as the file kept here goes.  The largest north-south force moved from %.3f N/m at degree 4 to %.3f at degree 8 and %.3f at degree 12."
  % (A_E/r_top, NMAX_FILE, earth[4]["fN_absmax"], earth[8]["fN_absmax"], earth[12]["fN_absmax"]))
P("    The same acceleration acts on the cable: on its %.0f kg/m the largest is %.3f N/m.  The cable passes through the pattern at its own speed, so it feels the part that changes with longitude as a load that changes in time ; that is not computed here."
  % (M_CABLE_M, earth[NMAX_FILE]["gN_absmax"]*M_CABLE_M))
P("    [Statics only, not simulated] A sideways acceleration that is the same all the way around moves the plane of a ring that circles the Earth off the equator by a r^3 / GM = %.1f m, where gravity's own pull back toward the equator balances it." % z0)
OUT["earth"] = dict(by_degree=earth, nmax=NMAX_FILE, fN_mean=float(fN.mean()), fN_var_max=float(np.abs(fN - fN.mean()).max()), fE_tension_pp_MN=float((cum.max() - cum.min())/1e6), z0_m=z0,
                    every30=[dict(lon=float(np.degrees(lam[i])), fN=float(fN[i]), fE=float(fE[i])) for i in range(0, len(lam), 300)], check_g_surface=g_surf)

# ---- (c) and (d)
fT = a_worst*M_CAS_M; fEm = earth[NMAX_FILE]["fN_absmax"]; fEz = abs(float(fN.mean())); fEv = float(np.abs(fN - fN.mean()).max())
P("\n(c) What carries them in layout F.  Sideways loads on the casing: the tide at its worst %.5f N/m ; the Earth's field %.4f N/m at its largest (%.4f the same all around, up to %.4f changing with longitude)." % (fT, fEm, fEz, fEv))
P("    layout | guide lines' pendulum hold  | reaction pairs' sideways hold [estimate] | both, for each meter of ring | the casing's sideways shift for the lines to carry: the tide | the Earth's field")
carry = {}
for Y in LAY:
    if not Y["key"].startswith("F"): continue
    kg_m = Y["hold_guide"]/Lr; kr_m = (Y["hold_ns"] - Y["hold_guide"])/Lr; kt = Y["hold_ns"]/Lr
    carry[Y["key"]] = dict(guide_per_m=kg_m, react_per_m=kr_m, both_per_m=kt, shift_tide_m=fT/kt, shift_earth_m=fEm/kt)
    P("    %-6s | %.1f N/m, %.3e /m     | %.1f N/m, %.3e /m                  | %.3e N/m for each m     | %.0f m                                                      | %.0f m" % (Y["key"], Y["hold_guide"], kg_m, Y["hold_ns"] - Y["hold_guide"], kr_m, kt, fT/kt, fEm/kt))
P("    The lines are far too soft to hold the casing in place sideways: for them to carry these loads by the casing's shift, the casing would have to move hundreds of meters to kilometers.")
P("    A line can carry a steady sideways load another way, by where its foot is put.  The Earth's field does not change in the Earth's frame, and the casing does not move in it.")
for gs in (100e3, 200e3):
    Hn = fEm*gs; off = Hn/k_g
    P("    A guide line every %.0f km that takes the whole of the largest load on its stretch: %.1f kN sideways, %.2f %% of its %.0f kN at the ground ; its foot %.1f km north or south of the point under its top, a lean of %.2f deg.  For the tide's %.5f N/m: %.2f kN, %.2f km, changing with the Moon."
      % (gs/1e3, Hn/1e3, 100*Hn/T_g_gnd, T_g_gnd/1e3, off/1e3, math.degrees(math.atan(off/h_r)), fT, fT*gs/1e3, fT*gs/k_g/1e3))
    carry["foot_%d" % (gs/1e3)] = dict(H_kN=Hn/1e3, offset_km=off/1e3, lean_deg=math.degrees(math.atan(off/h_r)), tide_kN=fT*gs/1e3, tide_km=fT*gs/k_g/1e3)
P("    What is left for the lateral guidance: between the lines' hold points the casing's sideways load has nowhere to go but through the guidance to the cable, the whole %.4f N/m at the worst place and %.5f N/m of tide,"
  % (fEm, fT))
P("    unless the casing is built to carry it along its own length to the hold points (not computed).  With no sideways hold from the lines at all, the guidance carries the same.")
OUT["carry"] = carry
P("\n(d) The track's flatness.  A load f for each meter on the cable, held at points L apart, bows it by f L^2 / (8 T).  T = %.1f GN (ring_altitude.CABLE_TENSION ; the book prints 683 GN)." % (T_CABLE/1e9))
tol = 41e-6; f_tol = 8*T_CABLE*tol/100.0**2
P("    Chapter 2's tolerance, %.0f micrometers over 100 m, is a curvature.  The sideways load that bows the cable by that much over 100 m is 8 T x 41 um / (100 m)^2 = %.0f N/m." % (tol*1e6, f_tol))
bow = {}
P("    load on the cable                                         | N/m      | bow over 100 m     | of the 41 um | bow over 100 km | over 200 km")
for nm, f in (("tide on the casing, worst", fT), ("tide on the cable and casing, worst", a_worst*(M_CAS_M + M_CABLE_M)), ("Earth's field on the casing, largest", fEm), ("Earth's field on the cable and casing, largest", earth[NMAX_FILE]["gN_absmax"]*(M_CAS_M + M_CABLE_M)), ("the book's 0.1 N/m", 0.1)):
    b100 = f*100.0**2/(8*T_CABLE); bow[nm] = dict(f=f, bow_100m=b100, share=b100/tol, bow_100km=f*100e3**2/(8*T_CABLE), bow_200km=f*200e3**2/(8*T_CABLE))
    P("    %-57s | %.5f | %.3e m (%.2e um) | %.1e   | %.2f mm         | %.2f mm" % (nm, f, b100, b100*1e6, b100/tol, bow[nm]["bow_100km"]*1e3, bow[nm]["bow_200km"]*1e3))
P("    The bow over the full spacing is the cable's shift between hold points, a long gentle curve.  Over any 100 m of it the track departs from a straight line by the 100 m figure, whatever the spacing.")
v_rel = R.V_CABLE + R.V_GROUND_SYNC; mv2 = M_CABLE_M*v_rel**2
P("    A caution on the formula.  It is for a cable at rest.  The cable moves past the casing at %.0f m/s, and its mass for each meter times that speed squared is %.2e N, %.1f times its tension." % (v_rel, mv2, mv2/T_CABLE))
P("    Its real path under a sideways load that stands still in the Earth's frame is set by its motion as much as by its tension.  That path is not computed here ; the figures above are the brief's formula and no more.")
OUT["bow_caution"] = dict(v_rel=v_rel, m_v2=mv2, ratio=mv2/T_CABLE)
P("    The lateral guidance: Chapter 7 (Section 7.5) states it as a requirement on the LIMs and gives no figure ; 'It has not computed it' and 'how large a force a LIM of this kind can make across its gap is an open item'.")
P("    What it has to carry to keep the casing centered on the cable is the casing's sideways load: up to %.4f N/m from the Earth's field (degree %d) and %.5f N/m of tide, against the casing's weight of 95.0 kN/m (%.1e of it)." % (fEm, NMAX_FILE, fT, (fEm + fT)/95.0e3))
OUT["bow"] = bow; OUT["f_tolerance"] = f_tol

# =====================================================================================================================
P("\n" + "="*118 + "\n3. The conductors on the guide lines\n" + "="*118)
SIG = 5.8e7; V = 1e6; LOSS = 0.02; PW = 900e9
L16 = json.load(open(os.path.join(ROOT, "runs", "round16_logs", "round16_number_sheet.json"), encoding="utf-8")); c16 = L16["J"]["ring panels' export, 900 GW at 1 MV"]
NQ = 400000; r_mid = RE + (np.arange(NQ) + 0.5)*h_r/NQ; dr = h_r/NQ; dP = potv(r_top); phi_mid = potv(r_mid)
P("900 GW at 1 MV, a loss of 2 %% down and back, CNT at copper's conductivity (%.0f MS/m), 700 km down and 700 km back.  The guide lines are the poles, one pole on each, alternating along the ring," % (SIG/1e6))
P("so a pole pair is two neighboring guide lines.  The run along the ring between them is not counted (round 16 did not count it either).")
P("  guide line every | pole pairs | at each pair | current in each line | each conductor        | one conductor | for each meter of ring | uniform conductor hanging alone: stress at the top | its top tension | of the rule's %.0f kN | with a 100 t tip load: pendulum" % (RULE/1e3))
cond = {}
for gs in (50e3 / 2, 100e3, 200e3):
    npairs = Lr/gs/2; Pp = PW/npairs; I = Pp/V; A_ = 2*h_r*I/(SIG*LOSS*V); m1 = RHO*A_*h_r; kgm = 2*m1*npairs/Lr
    s_top = RHO*dP; T_top = A_*s_top; T_tip = 100e3*G0; Tm = T_tip + RHO*A_*phi_mid; kpend = 1.0/(r_top*r_top*float((1/(r_mid*r_mid*Tm)).sum()*dr))
    cond[gs] = dict(spacing_km=gs/1e3, pairs=npairs, P_MW=Pp/1e6, I_A=I, A_cm2=A_*1e4, d_mm=2e3*math.sqrt(A_/math.pi), t_each=m1/1e3, kgm=kgm, stress_top_GPa=s_top/1e9, T_top_MN=T_top/1e6, rule=T_top/RULE,
                    tip_spare_t=(S_OP - s_top)*A_/G0/1e3, k_pend=kpend, hold_ring=(Lr/gs)*kpend)
    c = cond[gs]
    P("  %6.0f km        | %7.1f    | %7.1f MW   | %8.1f A           | %6.2f cm2, %6.2f mm | %7.0f t      | %.2f kg/m             | %.2f GPa of 12.5                                    | %.1f MN         | %.0f times           | %.1f N/m"
      % (gs/1e3, npairs, Pp/1e6, I, c["A_cm2"], c["d_mm"], c["t_each"], kgm, s_top/1e9, c["T_top_MN"], c["rule"], kpend))
assert abs(cond[25e3]["kgm"] - c16["kgm"]) < 1e-6 and abs(cond[25e3]["A_cm2"] - c16["A_cm2"]) < 1e-9
P("  (The first row is round 16's, a pole pair at every 50 km anchor point: the same %.2f kg/m.  The mass for each meter of ring does not depend on the spacing.)" % c16["kgm"])
P("  A conductor of one cross-section all the way up, hanging by itself, is at rho x (the potential from the ground to the ring) = %.2f GPa at its top, %.0f %% of the working stress, whatever its size." % (RHO*dP/1e9, 100*RHO*dP/S_OP))
P("  So it can hang as a line of its own at 12.5 GPa if the conducting CNT is as strong as the lines' CNT.  The repository has no strength for CNT made to conduct like copper ; that is not known here.")
for gs in (100e3, 200e3):
    c = cond[gs]
    P("  At %.0f km: the conductor has %.1f times the cross-section of the guide line at the ring (%.2f cm2), and %.0f t against the guide line's %.0f t.  It could carry a tip load of %.0f t beside its own weight within 12.5 GPa, so it would be the guide line."
      % (gs/1e3, c["A_cm2"]/(T_g_top/S_OP*1e4), T_g_top/S_OP*1e4, c["t_each"], m_g, c["tip_spare_t"]))
    P("     As a pendulum with a 100 t tip load it gives %.1f N/m in place of %.2f, and all of them %.0f N/m along the ring in place of %.0f.  Its %.1f MN at the top is %.0f attachments by the rule of a tenth, against 8 for a guide line."
      % (c["k_pend"], k_g, c["hold_ring"], (Lr/gs)*k_g, c["T_top_MN"], math.ceil(c["rule"])))
P("  Not run: layout F's swing and offsets with the conductors' pendulum hold added.")
OUT["conductors"] = {str(int(k/1e3)): v for k, v in cond.items()}

def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items() if not isinstance(v, np.ndarray)}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
P("\n(%.0f s)" % (time.time() - T_START))
with open(os.path.join(LOG_DIR, "round24_layout_f.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
with open(os.path.join(LOG_DIR, "round24_layout_f.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
print("log written: runs/round24_logs/round24_layout_f.txt and .json")
