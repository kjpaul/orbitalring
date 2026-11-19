#!/usr/bin/env python3
"""
Round 23: the hold along the ring from reaction lines against guide lines.
Brief: Vol_III/Claude_Code_instructions_round23_2026-11-19.md.   Book: Orbital Ring Engineering, Volume III: Chapter 11 (Sections 11.5 to 11.7).

The line model is round 9's (tools/round9_lines_array.py): the constant-stress hanging line with gravity that changes with height and the
centrifugal term on a round Earth (function hang), and the same line as an elastic line with the casing moved along the ring (function spring).
Those functions and their constants are taken from that file as they stand (by ast, so that its Parts 1 and 2 are not run again).
The casing's integrator, the simulated design launch and its braking are round 20's (tools/round20_brake_felt_pods.py), taken as they stand.

  1. the line against the distance of its anchorage along the ring: size, angles, stiffness of its top along the ring, hold for each tonne
  2. both directions: the pull and the stiffness against the casing's place, toward the anchorage until the line lies on the ground
  3. five layouts: hold, period, swing, offset under traffic, mass, anchorages, attachments
  4. the tide: only what the repository holds (nothing) and what the chapter script outside it prints
  5. the conductors on the lines of layouts C and E

Run from Python/orbitalring:  python tools/round23_reaction_vs_guide.py          (about 10 minutes)
Log: runs/round23_logs/round23_reaction_vs_guide.txt and .json
"""
import ast, contextlib, io, json, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root
import ring_altitude as R

T_START = time.time()
LOG_DIR = os.path.join(ROOT, "runs", "round23_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
OUT = {}

# ---------------------------------------------------------------------------------------------------------------------
# round 9's line model, as it stands
SRC9 = open(os.path.join(HERE, "round9_lines_array.py"), encoding="utf-8").read()
WANT_F = {"geff", "pot", "dia", "potv", "hang", "spring"}
WANT_A = {"mu", "RE", "W_E", "RHO", "S_OP", "S_BR", "S75", "E_CNT", "h_r", "h_a", "L", "La", "r_top", "dP_g", "TAPER", "NQ", "r_mid", "dr", "r_edge", "EXPO_mid", "EXPO_edge"}
ns = dict(np=np, math=math, R=R, root=root, brentq=brentq, P=lambda s="": None)
got = set()
for node in ast.parse(SRC9).body:
    if isinstance(node, ast.FunctionDef) and node.name in WANT_F:
        exec(compile(ast.Module(body=[node], type_ignores=[]), "round9_lines_array.py", "exec"), ns); got.add(node.name)
    elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id in WANT_A for t in node.targets):
        exec(compile(ast.Module(body=[node], type_ignores=[]), "round9_lines_array.py", "exec"), ns); got.update(t.id for t in node.targets)
assert got == WANT_F | WANT_A, (WANT_F | WANT_A) - got
hang = ns["hang"]; spring = ns["spring"]; pot = ns["pot"]; geff = ns["geff"]; dia = ns["dia"]
RE = ns["RE"]; RHO = ns["RHO"]; S_OP = ns["S_OP"]; E_CNT = ns["E_CNT"]; h_r = ns["h_r"]; Lr = ns["L"]; r_top = ns["r_top"]; TAPER = ns["TAPER"]
r_mid = ns["r_mid"]; dr = ns["dr"]; EXPO_mid = ns["EXPO_mid"]
L9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))["reaction"]
L16 = json.load(open(os.path.join(ROOT, "runs", "round16_logs", "round16_number_sheet.json"), encoding="utf-8"))
n_st = Lr/50e3; F_LAUNCH = 83.29e6; F1 = F_LAUNCH/n_st; g0 = geff(0)
B9 = L9["cases"]["B. anchorage 404 km away (700 tan 30)"]; K9 = L9["spring"]["B. anchorage 404 km away"]["k_N_per_m"]
assert abs(F1 - B9["F"]) < 1e-6

P("Round 23.  The hold along the ring from reaction lines against guide lines.")
P("CNT %.0f kg/m3, working stress %.1f GPa, modulus %.0f GPa ; ring at %.0f km ; taper %.4f ; anchor points %.2f at 50 km ; push of the design launch %.2f MN, %.3f kN for each anchor point (Equation 11.4)."
  % (RHO, S_OP/1e9, E_CNT/1e9, h_r/1e3, TAPER, n_st, F_LAUNCH/1e6, F1/1e3))

def solve_d(F, d_km):
    """Round 9's hanging line that delivers F along the ring at the casing and reaches the ground d_km along the ring."""
    f = lambda b: hang(F, b)["ground_km"] - d_km
    return brentq(f, 0.5, 89.99, xtol=1e-10)

def walk(cd, xs):
    """Round 9's spring on a list of casing places, kept in the order given (each solve starts from the one before).  Returns {x: (H, T_ground, beta_ground, T_top)}."""
    out, L0 = spring(cd, xs, say=False)
    return {r[0]: (r[1:] if r[1] is not None else None) for r in out}, L0

def k_at(cd, x=0.0, h=100.0):
    """The stiffness of the top along the ring at the casing's place x: the change of the force along the ring over 2 h, both solves of round 9's spring."""
    t, L0 = walk(cd, [x, x + h, x - h] if x >= 0 else [0.0] + list(np.arange(-5e3, x, -5e3)) + [x, x + h, x - h])
    if t[x] is None or t[x + h] is None or t[x - h] is None: return None
    return dict(H=t[x][0], k=(t[x + h][0] - t[x - h][0])/(2*h), T_g=t[x][1], beta_g=t[x][2], T_top=t[x][3], L0_km=L0/1e3)

# =====================================================================================================================
P("\n" + "="*118 + "\n1. The line against the distance of its anchorage along the ring\n" + "="*118)
P("""How the stiffness is computed.  Round 9's function spring: the line is built (its unstretched length and its cross-section all the way up)
for its design shape at 12.5 GPa.  It is then an elastic hanging line (modulus 500 GPa) whose top is put at the casing's radius, x meters
along the ring from the design place.  The two unknowns, the tension at the ground and the line's moment about the Earth's center, are found
for each x.  The force along the ring at the casing is read at x = +100 m and x = -100 m, and the stiffness is the difference over 200 m.
The stretch of the line and the change of its hanging shape are both in it.""")
cd9 = hang(F1, B9["beta_g"], want_path=True)
k9 = k_at(cd9)
P("Check against round 9: anchorage %.2f km, %.3f kN: stiffness %.4f N/m (round 9's log: %.4f) ; mass %.2f t ; %.2f mm at the ground." % (cd9["ground_km"], F1/1e3, k9["k"], K9, cd9["mass_t"], cd9["d_g"]))
assert abs(k9["k"] - K9) < 1e-6

# the vertical guide line, 100 t tip load: round 16's pendulum, and the same line in round 9's elastic model with a small lean
A_tip = 100e3*g0/S_OP; T_g100 = 100e3*g0
m100 = float((RHO*A_tip*EXPO_mid).sum()*dr); T100 = T_g100*TAPER; k_p = 1.0/float((1/(T_g100*EXPO_mid)).sum()*dr)
assert abs(k_p - L16["H"]["k_p"]) < 2e-4 and abs(m100/1e3 - L16["H"]["m100_t"]) < 0.02, (k_p, L16["H"]["k_p"], m100)
P("\nThe vertical line under a 100 t tip load (the guide line): %.2f mm at the ground, %.2f mm at the ring, %.2f t, %.1f kN at the ground, %.1f kN at the top." % (dia(T_g100, S_OP), dia(T100, S_OP), m100/1e3, T_g100/1e3, T100/1e3))
P("  As a pendulum, 1 / (the sum of ds / T up the line): %.4f N/m (round 16: %.4f).  A vertical line does not stretch when its top moves sideways." % (k_p, L16["H"]["k_p"]))
pend_chk = {}
for lean_km in (1.0, 10.0):
    b = brentq(lambda b: hang(T_g100*RE*math.sin(math.radians(b))/r_top, b)["ground_km"] - lean_km, 0.001, 5.0, xtol=1e-12)
    c = hang(T_g100*RE*math.sin(math.radians(b))/r_top, b, want_path=True); kk = k_at(c)
    pend_chk[lean_km] = dict(k=kk["k"], F=c["F"], mass_t=c["mass_t"])
    P("  The same line in round 9's elastic model, its foot %.0f km along the ring from under its top: %.4f N/m (it then pulls %.1f N along the ring)." % (lean_km, kk["k"], c["F"]))
k_p_round = 1.0/(r_top*r_top*float((1/(r_mid*r_mid*T_g100*EXPO_mid)).sum()*dr))
P("  On the round Earth the same pendulum is 1 / (r_top^2 x the sum of dr / (r^2 T)): %.4f N/m.  The line's moment about the Earth's center is the same all the way up," % k_p_round)
P("  so its sideways pull at the casing is R_E / r_top = %.4f of its pull at the ground.  Round 9's elastic model agrees with this figure, not with the %.4f N/m." % (RE/r_top, k_p))
P("  The chapter script's %.2f N/m (ch15_numbers.py, a flat Earth) is therefore %.1f %% higher than the guide line gives in the model the reaction line is computed in." % (k_p, 100*(k_p/k_p_round - 1)))
OUT["guide"] = dict(k_p_round=k_p_round, per_tonne_round=k_p_round/(m100/1e3),mass_t=m100/1e3, d_g=dia(T_g100, S_OP), d_top=dia(T100, S_OP), T_g_kN=T_g100/1e3, T_top_kN=T100/1e3, k_p=k_p, per_tonne=k_p/(m100/1e3), check_elastic=pend_chk)

b_max = hang(F1, 89.99); d_max = b_max["ground_km"]
P("\nThe farthest anchorage a line at 12.5 GPa can have (the line lying flat at the ground): %.1f km along the ring ; the line is then %.1f deg from the vertical at the casing." % (d_max, b_max["beta_top"]))
P("\nLeaning lines sized for %.3f kN along the ring at the casing (Equation 11.4), 12.5 GPa everywhere at that pull:" % (F1/1e3))
P("  anchorage | at the ground        | at the casing        | across (ground, ring) | length   | mass     | down-pull at the casing | stiffness along the ring | hold for each tonne | against the guide line | north-south pendulum [estimate]")
rows = []
row0 = dict(d_km=0.0, beta_g=0.0, beta_top=0.0, T_g_kN=T_g100/1e3, T_top_kN=T100/1e3, d_g=dia(T_g100, S_OP), d_top=dia(T100, S_OP), length_km=h_r/1e3, mass_t=m100/1e3, V_top_kN=T100/1e3, k=k_p, per_tonne=k_p/(m100/1e3), k_ns=k_p, sized="100 t tip load")
P("  %6.1f km | %5.2f deg, %7.1f kN | %5.2f deg, %7.1f kN | %5.2f, %5.2f mm        | %6.1f km | %7.2f t | %7.1f kN              | %7.4f N/m              | %.5f N/m per t    | %5.2f                  | %.4f N/m      (the guide line, 100 t tip load)"
  % (0.0, 0.0, T_g100/1e3, 0.0, T100/1e3, row0["d_g"], row0["d_top"], h_r/1e3, m100/1e3, T100/1e3, k_p, row0["per_tonne"], 1.0, k_p))
rows.append(row0); DES = {}
for d in (100.0, 200.0, 300.0, B9["ground_km"], 500.0, 600.0):
    if d >= d_max:
        P("  %6.1f km | no such line at 12.5 GPa: the farthest anchorage is %.1f km" % (d, d_max)); rows.append(dict(d_km=d, exists=False)); continue
    b = solve_d(F1, d); c = hang(F1, b, want_path=True); kk = k_at(c)
    em = EXPO_mid; sb = c["K"]/(r_mid*c["T_g"]*em); cb = np.sqrt(1 - sb*sb); k_ns = 1.0/float((1/(c["T_g"]*em)/cb).sum()*dr)      # 1 / sum of ds / T along the hanging line
    row = dict(d_km=d, exists=True, beta_g=b, beta_top=c["beta_top"], T_g_kN=c["T_g"]/1e3, T_top_kN=c["T_top"]/1e3, d_g=c["d_g"], d_top=c["d_top"], length_km=c["length_km"], mass_t=c["mass_t"], V_top_kN=c["V_top_kN"],
               sag_km=c["sag_km"], k=kk["k"], per_tonne=kk["k"]/c["mass_t"], k_ns=k_ns, sized="%.3f kN along the ring" % (F1/1e3), built_km=kk["L0_km"])
    rows.append(row); DES[d] = c
    P("  %6.1f km | %5.2f deg, %7.1f kN | %5.2f deg, %7.1f kN | %5.2f, %5.2f mm        | %6.1f km | %7.2f t | %7.1f kN              | %7.4f N/m              | %.5f N/m per t    | %5.2f                  | %.4f N/m"
      % (d, b, row["T_g_kN"], row["beta_top"], row["T_top_kN"], row["d_g"], row["d_top"], row["length_km"], row["mass_t"], row["V_top_kN"], row["k"], row["per_tonne"], row["per_tonne"]/row0["per_tonne"], k_ns))
OUT["by_distance"] = rows; OUT["d_max_km"] = d_max
P("  Against the guide line in the same round-Earth model (%.4f N/m, %.5f N/m for each tonne): %s."
  % (k_p_round, k_p_round/(m100/1e3), " ; ".join("%.0f km %.2f times" % (r_["d_km"], r_["per_tonne"]/(k_p_round/(m100/1e3))) for r_ in rows[1:] if r_.get("exists"))))
P("  (Angles are from the vertical.  North-south pendulum: 1 / the sum of ds / T along the line, the pull back for each meter the top is moved out of the line's plane.")
P("   It is an estimate: round 9's model is in the plane of the ring and does not carry a sideways shift.)")

# what the stiffness is made of, and how it goes with the size
c404 = DES[B9["ground_km"]]
E_keep = ns["E_CNT"]; ns["E_CNT"] = E_keep*1e4
c_stiff = hang(F1, B9["beta_g"], want_path=True); k_inext = k_at(c_stiff)["k"]
ns["E_CNT"] = E_keep
th = math.radians(30); k_straight = F1/((S_OP/E_CNT)*(h_r/math.cos(th))/math.sin(th))
k_el = 1/(1/k9["k"] - 1/k_inext)
P("\nWhat the stiffness is made of (the design line, anchorage %.1f km, %.3f kN):" % (B9["ground_km"], F1/1e3))
P("  with the modulus 10,000 times higher (a line that cannot stretch, shape only): %.4f N/m" % k_inext)
P("  with the real modulus: %.4f N/m.  As two springs in a row, the stretch alone is then worth %.4f N/m (1/k = 1/k_shape + 1/k_stretch)." % (k9["k"], k_el))
P("  the book's old straight line at 30 deg that only stretches: %.4f N/m." % k_straight)
OUT["parts"] = dict(k=k9["k"], k_shape_only=k_inext, k_stretch_in_series=k_el, k_straight_30=k_straight)
c6 = hang(6*F1, B9["beta_g"], want_path=True); k6 = k_at(c6)
c2 = hang(2*F1, B9["beta_g"], want_path=True); k2 = k_at(c2)
P("\nThe size at a fixed anchorage and working stress.  The same line for 2 and 6 times the pull:")
for nm, c, kk in (("1 x", cd9, k9), ("2 x", c2, k2), ("6 x", c6, k6)):
    P("  %s: %.2f kN along the ring ; %.2f mm at the ground, %.2f mm at the ring ; %.2f t ; anchorage %.2f km ; %.5f N/m ; %.6f N/m for each tonne" % (nm, c["F"]/1e3, c["d_g"], c["d_top"], c["mass_t"], c["ground_km"], kk["k"], kk["k"]/c["mass_t"]))
P("  The stiffness is %.5f and %.5f times the design line's for 2 and 6 times the cross-section: linear in the area, as expected." % (k2["k"]/k9["k"], k6["k"]/k9["k"]))
OUT["scaling"] = dict(k1=k9["k"], k2=k2["k"], k6=k6["k"], m1=cd9["mass_t"], m2=c2["mass_t"], m6=c6["mass_t"])
P("  (%.0f s so far)" % (time.time() - T_START))

# =====================================================================================================================
P("\n" + "="*118 + "\n2. Both directions\n" + "="*118)
P("""x is the casing's place along the ring from the line's design place, positive away from the anchorage (the way a launch pushes).
The line is built for its design shape and is not changed.  Its pull and its stiffness at each place, from round 9's spring:""")
xs_pos = [0.0, 100.0, 2000.0, 2100.0, 1900.0, 5000.0, 5100.0, 4900.0, 10000.0, 10100.0, 9900.0, 20000.0, 20100.0, 19900.0]
tp, _ = walk(cd9, xs_pos)
neg_pts = [-2e3, -5e3, -10e3, -20e3, -40e3, -60e3, -80e3, -100e3, -120e3, -140e3, -160e3, -180e3, -200e3, -220e3, -240e3, -260e3, -280e3, -300e3]
xs_neg = [0.0, -100.0]
for x in neg_pts: xs_neg += [x + 100.0, x, x - 100.0]
tn, _ = walk(cd9, xs_neg)
tab = {}
def put(x, t):
    if t.get(x) is None or t.get(x + 100.0) is None or t.get(x - 100.0) is None: return False
    tab[x] = dict(H=t[x][0], k=(t[x + 100.0][0] - t[x - 100.0][0])/200.0, T_g=t[x][1], beta_g=t[x][2], T_top=t[x][3]); return True
tall = dict(tn); tall.update(tp)
for x in (20e3, 10e3, 5e3, 2e3, 0.0): put(x, tall)
last_ok = 0.0; first_bad = None
for x in neg_pts:
    if put(x, tn): last_ok = x
    else:
        first_bad = x; break
P("  casing's place | pull along the ring (design line) | of design | stiffness | of design | at the ground      | top tension, of the working stress | the line six times heavier: pull, stiffness")
for x in sorted(tab, reverse=True):
    r = tab[x]
    P("  %+9.1f km   | %8.2f kN                       | %5.1f %%   | %.4f N/m | %5.1f %%   | %5.2f deg, %6.1f kN | %6.1f kN, %5.1f %%                  | %7.1f kN, %.3f N/m"
      % (x/1e3, r["H"]/1e3, 100*r["H"]/F1, r["k"], 100*r["k"]/k9["k"], r["beta_g"], r["T_g"]/1e3, r["T_top"]/1e3, 100*r["T_top"]/cd9["T_top"], 6*r["H"]/1e3, 6*r["k"]))
# where the line lies down at the ground: refine between the last place solved and the first that was not
slack = None
if first_bad is not None:
    lo, hi = last_ok, first_bad
    for _ in range(6):
        mid = 0.5*(lo + hi); t, _ = walk(cd9, [0.0] + [x for x in neg_pts if x > mid] + [mid])
        if t[mid] is not None and t[mid][2] < 89.9: lo = mid; last = t[mid]
        else: hi = mid
    t, _ = walk(cd9, [0.0] + [x for x in neg_pts if x > lo] + [lo]); last = t[lo]
    slack = dict(x_last=lo, x_fail=hi, H=last[0], beta_g=last[2], T_g=last[1])
    P("  The last place solved toward the anchorage is %.1f km ; at %.1f km no hanging shape is found.  There the line is %.1f deg from the vertical at the ground and pulls %.2f kN (%.1f %% of design), %.2f kN at the ground."
      % (lo/1e3, hi/1e3, last[2], last[0]/1e3, 100*last[0]/F1, last[1]/1e3))
else:
    xl = min(tab); P("  The line is still hanging clear of the ground at %.0f km, the last place tried (%.1f deg at the ground, %.1f %% of the design pull)." % (xl/1e3, tab[xl]["beta_g"], 100*tab[xl]["H"]/F1))
P("  The line six times heavier has the same shape at every place: six times the pull and six times the stiffness, at the same distances (Section 1).")
OUT["both"] = dict(table={str(k): v for k, v in tab.items()}, slack=slack)
xs_t = np.array(sorted(tab)); Hs_t = np.array([tab[x]["H"] for x in xs_t]); ks_t = np.array([tab[x]["k"] for x in xs_t])
P("\n  The stiffness against the standing pull (read from the table by interpolation):")
frac = {}
for f in (1.0, 0.75, 0.5, 0.25):
    xx = float(np.interp(f*F1, Hs_t, xs_t)); kk = float(np.interp(xx, xs_t, ks_t)); frac[f] = dict(x=xx, k=kk)
    P("    standing at %3.0f %% of the design pull (%.1f kN): the casing %.1f km toward the anchorage from the design place ; stiffness %.3f N/m, %.0f %% of the %.3f N/m at the full pull" % (100*f, f*F1/1e3, -xx/1e3, kk, 100*kk/k9["k"], k9["k"]))
OUT["both"]["by_standing_pull"] = {str(k): v for k, v in frac.items()}
# layout A as built: one reaction line against two pendulums
xa = brentq(lambda x: float(np.interp(x, xs_t, Hs_t)) + 2*k_p*x, xs_t[0], 0.0); Ha = float(np.interp(xa, xs_t, Hs_t)); ka = float(np.interp(xa, xs_t, ks_t))
P("\n  The chapter's layout, if the reaction line is built for its design shape: its standing pull is held by the two guide lines as pendulums.")
P("    They balance with the casing %.1f km toward the anchorages: reaction line %.1f kN (%.0f %% of design), stiffness %.3f N/m ; guide lines leaning %.2f deg.  Hold there %.0f N/m, against %.0f N/m at the design place."
  % (-xa/1e3, Ha/1e3, 100*Ha/F1, ka, math.degrees(math.atan(-xa/h_r)), n_st*(ka + 2*k_p), n_st*(k9["k"] + 2*k_p)))
OUT["both"]["layout_A_rest"] = dict(x=xa, H=Ha, k=ka, hold=n_st*(ka + 2*k_p))
P("""  Opposed pairs.  Two lines of one size at an anchor point, one leaning each way, each built for its design shape: at the design place their
    standing pulls cancel, so nothing else has to hold them, and the casing rests there.  When the casing moves x, one line gains k x and the
    other loses k x, so the pair holds with 2 k in both directions, the same either way.  The hold for each tonne is that of the single line.""")
sw = 2000.0; r2 = tab[2e3]; rm2 = tab[-2e3]
P("    At a swing of 2 km: the line pulled tighter is at %.1f %% of its working stress at the top, the other at %.1f %% ; the pair's net pull %.2f kN, 2 k x = %.2f kN."
  % (100*r2["T_top"]/cd9["T_top"], 100*rm2["T_top"]/cd9["T_top"], (r2["H"] - rm2["H"])/1e3, 2*k9["k"]*sw/1e3))
r5 = tab[5e3]; rm5 = tab[-5e3]
P("    At 5 km: %.1f %% and %.1f %% ; net pull %.2f kN, 2 k x = %.2f kN." % (100*r5["T_top"]/cd9["T_top"], 100*rm5["T_top"]/cd9["T_top"], (r5["H"] - rm5["H"])/1e3, 2*k9["k"]*5e3/1e3))
OUT["both"]["pair"] = dict(over_at_2km=r2["T_top"]/cd9["T_top"], over_at_5km=r5["T_top"]/cd9["T_top"], net_2km=r2["H"] - rm2["H"], net_5km=r5["H"] - rm5["H"])
P("  (%.0f s so far)" % (time.time() - T_START))

# =====================================================================================================================
P("\n" + "="*118 + "\n3. The layouts\n" + "="*118)
src20 = open(os.path.join(HERE, "round20_brake_felt_pods.py"), encoding="utf-8").read(); cut = src20.index('P("Round 20.'); assert cut > 0
g20 = {"__name__": "round20_prefix", "__file__": os.path.join(HERE, "round20_brake_felt_pods.py")}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src20[:cut], "round20_brake_felt_pods.py (prefix)", "exec"), g20)
g20["LOG_DIR"] = LOG_DIR
M_CAS = g20["M_CAS"]; F_LIM = g20["F_LIM"]; rebin = g20["rebin"]; swing_once = g20["swing_once"]; swing_daily = g20["swing_daily"]
with contextlib.redirect_stdout(io.StringIO()):
    D = g20["lsm_run"]("design_brake")
L20 = json.load(open(os.path.join(ROOT, "runs", "round20_logs", "round20_brake_felt_pods.json"), encoding="utf-8"))
assert abs(D["dur"] - L20["design"]["dur"]) < 1e-6 and abs(D["impulse"] - L20["design"]["impulse"]) < 1.0
Fl = D["F"]; Fb = D["FB"]; F_nb = rebin(Fl, 0.1, 10.0); F_br = rebin(np.concatenate((Fl, -Fb)), 0.1, 10.0); J_br = float(Fl.sum()*0.1 - Fb.sum()*0.1)
P("The design launch and its braking are lsm_simulation's (run here, the same as round 20's log): %.1f min, impulse %.4e N s ; braking %.1f min, %.4e N s.  Casing %.4e kg.  LIMs %.1f MN." % (D["dur"]/60, D["impulse"], D["brake"]["duration_s"]/60, D["brake"]["impulse"], M_CAS, F_LIM/1e6))
RULE = 0.1*95.0e3*30.0
P("Section 11.2's rule: no attachment adds more than a tenth of the casing's weight over 30 m, 0.1 x 95.0 kN/m x 30 m = %.0f kN." % (RULE/1e3))
m_r = cd9["mass_t"]; kr = k9["k"]; k_A = n_st*(kr + 2*k_p); mass_A = m_r + 2*m100/1e3
def rl(s):
    """The design reaction line with s times the cross-section (Section 1: everything of it goes with s)."""
    return dict(s=s, F_kN=s*F1/1e3, d_g=cd9["d_g"]*math.sqrt(s), d_top=cd9["d_top"]*math.sqrt(s), mass_t=s*m_r, T_g_kN=s*cd9["T_g"]/1e3, T_top_kN=s*cd9["T_top"]/1e3, V_top_kN=s*cd9["V_top_kN"], k=s*kr, k_ns=s*next(r_["k_ns"] for r_ in rows if abs(r_["d_km"] - B9["ground_km"]) < 1e-6))
sB = k_A/(n_st*kr); sC = k_A/(2*n_st*kr); sD = mass_A/(2*m_r)
LAY = [
 dict(key="A", name="A. The chapter: one reaction line and two guide lines", react=rl(1.0), n_react=1, n_guide=2, paired=False),
 dict(key="B", name="B. Reaction lines only, one every 50 km, for the same hold", react=rl(sB), n_react=1, n_guide=0, paired=False),
 dict(key="C", name="C. Reaction lines in opposed pairs, for the same hold", react=rl(sC), n_react=2, n_guide=0, paired=True),
 dict(key="D", name="D. Opposed pairs with the mass of layout A's lines", react=rl(sD), n_react=2, n_guide=0, paired=True),
 dict(key="E", name="E. As C, with one guide line kept at each anchor point", react=rl(sC), n_react=2, n_guide=1, paired=True),
]
for Y in LAY:
    r = Y["react"]; k = n_st*(Y["n_react"]*r["k"] + Y["n_guide"]*k_p); Y["hold"] = k
    Y["period_s"] = 2*math.pi*math.sqrt(M_CAS/k); c = 2*math.sqrt(k*M_CAS)
    Y["once"] = {nm: swing_once(F, 10.0, c, k=k) for nm, F in (("not braked", F_nb), ("braked", F_br))}
    Y["daily"] = {N: swing_daily(F_br, 10.0, N, c, False, k=k, J_exact=J_br) for N in (7, 14)}
    Y["mass_pt_t"] = Y["n_react"]*r["mass_t"] + Y["n_guide"]*m100/1e3; Y["kgm"] = Y["mass_pt_t"]*1e3/50e3; Y["all_t"] = Y["mass_pt_t"]*n_st
    Y["anchorages"] = n_st*(Y["n_react"] + Y["n_guide"]); Y["reach_km"] = (2 if Y["paired"] else 1)*cd9["ground_km"]
    Y["att_react"] = r["T_top_kN"]*1e3/RULE; Y["att_guide"] = T100/RULE if Y["n_guide"] else 0.0
    Y["hold_ns"] = n_st*(Y["n_react"]*r["k_ns"] + Y["n_guide"]*k_p)
    Y["standing_MN"] = 0.0 if Y["paired"] else n_st*r["F_kN"]/1e3
    xo = Y["once"]["not braked"]["x_max"]; Y["line_change_kN"] = r["k"]*xo/1e3; Y["line_change_share"] = r["k"]*xo/(r["F_kN"]*1e3)
    P("\n%s" % Y["name"])
    P("   reaction line: %.3f times the design line ; built for %.1f kN along the ring ; %.2f mm at the ground, %.2f mm at the ring ; %.1f kN at the ground, %.1f kN at the top ; %.2f t ; %.4f N/m each ; anchorage %.1f km away, length %.1f km"
      % (r["s"], r["F_kN"], r["d_g"], r["d_top"], r["T_g_kN"], r["T_top_kN"], r["mass_t"], r["k"], cd9["ground_km"], cd9["length_km"]))
    P("   by the chapter's rule for a busy ring (four design lines' cross-section is %.1f mm at the ground, carried as several lines in parallel): %.2f design lines' worth ; %s"
      % (2*cd9["d_g"], r["s"], "one line" if r["s"] <= 1.0 + 1e-9 else "%d lines of %.2f mm in parallel, or more thinner ones" % (math.ceil(r["s"] - 1e-9), cd9["d_g"]*math.sqrt(r["s"]/math.ceil(r["s"] - 1e-9)))))
    P("   hold along the ring %.1f N/m (%.1f x (%d x %.4f + %d x %.4f)) ; period %.0f s, %.2f h ; critical damping %.3e N s/m" % (k, n_st, Y["n_react"], r["k"], Y["n_guide"], k_p, Y["period_s"], Y["period_s"]/3600, c))
    for nm in ("not braked", "braked"):
        o = Y["once"][nm]
        P("   one design launch, LIMs damping, sled %-10s: farthest %.0f m, %.2f h after the launch begins ; within 1 %% after %.1f h ; LIM push %.2f MN" % (nm, o["x_max"], o["peak_h"], o["die_h"], o["F_lim_MN"]))
    for N in (7, 14):
        o = Y["daily"][N]
        P("   %2d braked design launches a day (mean push %.3f MN): the casing sits at %.0f m and swings %.0f m about it ; LIM push at peak %.1f MN" % (N, o["mean_push_MN"], o["x_mean"], o["swing"], o["F_lim_peak_MN"]))
    P("   a reaction line's pull changes by %.2f kN at the swing of one launch, %.1f %% of the pull it is built for" % (Y["line_change_kN"], 100*Y["line_change_share"]))
    P("   lines for each anchor point %.1f t ; for each meter of ring %.3f kg/m ; all %.0f t ; anchorages on the ground %.0f ; ground reached along the ring %.0f km%s"
      % (Y["mass_pt_t"], Y["kgm"], Y["all_t"], Y["anchorages"], Y["reach_km"], " (%.0f km each way)" % cd9["ground_km"] if Y["paired"] else " (one way)"))
    P("   at the top: a reaction line %.1f kN, %.2f attachments' worth by the rule (it hangs from the force spreader cable, down-pull %.1f kN)%s"
      % (r["T_top_kN"], Y["att_react"], r["V_top_kN"], " ; a guide line %.2f MN, %.2f, eight branches on the casing" % (T100/1e6, Y["att_guide"]) if Y["n_guide"] else " ; no line meets the casing directly"))
    P("   standing pull along the ring that nothing in the layout cancels, at the design place: %.1f MN" % Y["standing_MN"])
    P("   north-south pendulum hold of these lines [estimate, see Section 1]: %.0f N/m for each meter the casing moves sideways" % Y["hold_ns"])
# layout A against round 20 and the book
a20 = L20["design"]["casing"]; A = LAY[0]
assert abs(A["once"]["not braked"]["x_max"] - a20["not braked"]["x_max"]) < 0.5 and abs(A["once"]["braked"]["x_max"] - a20["braked"]["x_max"]) < 0.5
d20 = {(r_["N"], r_["sled"], r_["lims"]): r_ for r_ in a20["daily"]}
assert abs(A["daily"][7]["x_mean"] - d20[(7, "braked", "damping, critical")]["x_mean"]) < 0.5 and abs(A["daily"][14]["x_mean"] - d20[(14, "braked", "damping, critical")]["x_mean"]) < 0.5
k_A_round = n_st*(kr + 2*k_p_round); OUT["hold_A_round"] = k_A_round
P("\nWith the guide line's pendulum taken on the round Earth (%.4f N/m), layout A's hold is %.1f x (%.4f + 2 x %.4f) = %.0f N/m and its period %.2f h.  The layouts above keep the chapter's %.4f N/m, as the brief states it."
  % (k_p_round, n_st, kr, k_p_round, k_A_round, 2*math.pi*math.sqrt(M_CAS/k_A_round)/3600, k_p))
P("\nLayout A is round 20's row to the meter (%.0f m, %.0f m braked ; %.0f and %.0f m under traffic)." % (A["once"]["not braked"]["x_max"], A["once"]["braked"]["x_max"], A["daily"][7]["x_mean"], A["daily"][14]["x_mean"]))
# layout B: what happens to one-way lines with nothing to pull against
B = LAY[1]; sBv = B["react"]["s"]
P("\nLayout B as it would stand.  Its lines all lean one way and nothing pulls back, so the hold above is only there while each line is pulled with %.1f kN, %.1f MN in all." % (B["react"]["F_kN"], B["standing_MN"]))
B["rest"] = {}
for nm, Fm in (("7 braked launches a day", LAY[0]["daily"][7]["mean_push_MN"]*1e6), ("14 braked launches a day", LAY[0]["daily"][14]["mean_push_MN"]*1e6), ("the LIMs' whole %.1f MN" % (F_LIM/1e6), F_LIM)):
    Hneed = Fm/n_st/sBv
    if Hneed < Hs_t[0]:
        P("   held by %s (%.2f MN): each line would stand at %.1f %% of its design pull.  That is below the %.1f %% at which the line comes to lie on the ground at its anchorage (the casing %.0f km toward it), where the model ends."
          % (nm, Fm/1e6, 100*Hneed/F1, 100*slack["H"]/F1, -slack["x_last"]/1e3))
        B["rest"][nm] = dict(share=Hneed/F1, x=None)
    else:
        xx = float(np.interp(Hneed, Hs_t, xs_t)); kk = sBv*float(np.interp(xx, xs_t, ks_t))
        P("   held by %s (%.2f MN): each line stands at %.1f %% of its design pull, the casing %.0f km toward the anchorages from the design place ; hold there %.0f N/m, period %.1f h."
          % (nm, Fm/1e6, 100*Hneed/F1, -xx/1e3, n_st*kk, 2*math.pi*math.sqrt(M_CAS/(n_st*kk))/3600))
        B["rest"][nm] = dict(share=Hneed/F1, x=xx, hold=n_st*kk)
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items() if not isinstance(v, np.ndarray)}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
OUT["layouts"] = clean(LAY); OUT["rule_kN"] = RULE/1e3

# =====================================================================================================================
P("\n" + "="*118 + "\n4. The north-south load\n" + "="*118)
hits = []
for dp, dn, fn in os.walk(ROOT):
    if any(s in dp for s in (os.sep + "runs", os.sep + ".git", "__pycache__", os.sep + "lunar")): continue
    for f in fn:
        if f.endswith(".py") and f != os.path.basename(__file__):
            try: txt = open(os.path.join(dp, f), encoding="utf-8", errors="replace").read().lower()
            except OSError: continue
            if "tidal" in txt or " tide" in txt: hits.append(os.path.relpath(os.path.join(dp, f), ROOT))
P("Scripts of the repository that hold a tidal term (searched for 'tidal' and 'tide', the lunar mass driver's folder left out): %s" % (", ".join(hits) if hits else "none"))
CH = os.path.normpath(os.path.join(ROOT, "..", "..", "Vol_III", "renumber_scripts_2026-10-12", "spreader_loads_check.py")); tide_lines = []
if os.path.exists(CH):
    o = subprocess.run([sys.executable, CH], cwd=ROOT, capture_output=True, text=True, timeout=300).stdout
    tide_lines = [l.strip() for l in o.splitlines() if "tidal" in l]
    P("Outside the repository, the chapter script Vol_III/renumber_scripts_2026-10-12/spreader_loads_check.py (run as it stands, not edited) prints:")
    for l in tide_lines: P("   " + l)
P("The book's 'on the order of 0.1 N/m' is not reproduced from the repository.  Stopped here on this item: no sideways force is taken and no lean is computed.")
OUT["tide"] = dict(repository_scripts=hits, chapter_script_lines=tide_lines)

# =====================================================================================================================
P("\n" + "="*118 + "\n5. The conductors (as round 16, section J)\n" + "="*118)
SIG = 5.8e7
def cond(P_total, L1, L2, V=1e6, loss=0.02):
    """Round 16's cond (ch16_numbers.py, function anchor), with the two conductors' lengths given: one thickness for both, the current down L1 and back on L2."""
    Pp = P_total/n_st; I = Pp/V; A_ = (L1 + L2)*I/(SIG*loss*V); m = RHO*A_*(L1 + L2)
    return dict(P_MW=Pp/1e6, I_A=I, A_cm2=A_*1e4, d_mm=2e3*math.sqrt(A_/math.pi), t_pair=m/1e3, kgm=m*n_st/Lr, L1_km=L1/1e3, L2_km=L2/1e3)
c16 = L16["J"]["ring panels' export, 900 GW at 1 MV"]; ca = cond(900e9, h_r, h_r)
assert abs(ca["kgm"] - c16["kgm"]) < 1e-6 and abs(ca["A_cm2"] - c16["A_cm2"]) < 1e-9
Lre = cd9["length_km"]*1e3; OUT["conductors"] = {}
P("900 GW over %.1f anchor points at 1 MV, a loss of 2 %% down and back, CNT at copper's conductivity (%.0f MS/m).  The mass goes with the square of the path's length." % (n_st, SIG/1e6))
P("  layout and path                                              | at each point | current | each conductor      | pair    | for each meter of ring")
for nm, L1, L2 in (("A (the chapter): down one guide line, back on the other", h_r, h_r), ("C: down one reaction line of the pair, back on the other", Lre, Lre),
                   ("E: down the guide line, back on one reaction line", h_r, Lre), ("E: down one reaction line, back on the other (as C)", Lre, Lre)):
    c = cond(900e9, L1, L2); OUT["conductors"][nm] = c
    P("  %-60s | %6.1f MW     | %6.1f A | %5.2f cm2, %5.2f mm | %6.1f t | %.2f kg/m" % (nm, c["P_MW"], c["I_A"], c["A_cm2"], c["d_mm"], c["t_pair"], c["kgm"]))
rC = LAY[2]["react"]
P("  The CNT of the lines for their structural job: a reaction line of layout C %.3f cm2 at the ground and %.3f cm2 at the ring ; a guide line %.3f and %.3f cm2."
  % (rC["T_g_kN"]*1e3/S_OP*1e4, rC["T_top_kN"]*1e3/S_OP*1e4, A_tip*1e4, A_tip*TAPER*1e4))
cC = OUT["conductors"]["C: down one reaction line of the pair, back on the other"]
P("  In layout C a conductor is %.0f times the cross-section of its line at the ring, and a pair of conductors is %.1f t against %.1f t for the pair of lines." % (cC["A_cm2"]/(rC["T_top_kN"]*1e3/S_OP*1e4), cC["t_pair"], 2*rC["mass_t"]))
P("  Not run: the hanging shape and the stiffness of a reaction line that carries such a conductor.  The conductor's weight is far more than the line's own,")
P("  so the sizes and the hold of layouts C and E in Section 3 are for lines without it.")
OUT["conductors"]["structural_cm2"] = dict(react_C_ground=rC["T_g_kN"]*1e3/S_OP*1e4, react_C_top=rC["T_top_kN"]*1e3/S_OP*1e4, guide_ground=A_tip*1e4, guide_top=A_tip*TAPER*1e4)

P("\n(%.0f s)" % (time.time() - T_START))
with open(os.path.join(LOG_DIR, "round23_reaction_vs_guide.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
with open(os.path.join(LOG_DIR, "round23_reaction_vs_guide.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
print("log written: runs/round23_logs/round23_reaction_vs_guide.txt and .json")
