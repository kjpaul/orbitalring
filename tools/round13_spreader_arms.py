#!/usr/bin/env python3
"""
Round 13, item 5: the spreader arms below the kink as hanging lines (Chapters 11 and 12 of Volume III).
Round 9 (tools/round9_lines_array.py, part 3) sized the reaction lines as hanging lines.  The spreader arms below the kink were left as
straight lines: 893 km at 45 degrees, 32 MN at the kink, 2,100 t at 65.5 t per MN (Vol_III/renumber_scripts_2026-10-12/spreader_cases.py and
ch15_arm_choices.py).  This script sizes them the same way as the reaction lines: constant stress 12.5 GPa, hanging under their own weight,
gravity that changes with height, a round Earth.

The physics is round 9's and is not written out again.  round9_lines_array.py runs its whole study on import (about 4 minutes, and it rewrites
its own log), so its functions and constants are loaded from its source file by name (ast), without running the rest:
    geff, pot, dia, potv, hang, hang_xy, catenary, v2, case11 and the constants they use.
For the arm the top of the line is the kink at the array's altitude (600 km), not the ring (700 km); that is the only thing changed, and it is
changed by setting h_r and r_top before the part-3 constants of round 9 are evaluated.

Checks, in this order:
  0. the loaded functions reproduce round 9's logged rows A and B for the reaction lines (runs/round9_logs/round9_lines_array.json);
  1. the limit: a very strong line (and so a line whose weight does not matter) gives the straight line's numbers (45 deg, 632 km, 893 km, 32 MN);
  2. an independent step-by-step integration in the plane of the meridian, with the centrifugal term as it really is off the equator
     (round 9's hang() treats the plane as the equator's; the arm lies in a meridian plane).

Run from Python/orbitalring:  python tools/round13_spreader_arms.py
Log: runs/round13_logs/round13_spreader_arms.txt and .json
"""
import ast, io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root
import ring_altitude as R, power_config as cfg, power_physics as phys

LOG_DIR = os.path.join(ROOT, "runs", "round13_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n")
OUT = {}

# ---------------------------------------------------------------------------------------------------------------------
# Loading round 9's method without running round 9's study
# ---------------------------------------------------------------------------------------------------------------------
R9 = os.path.join(HERE, "round9_lines_array.py")
STAGE1 = ["mu", "RE", "W_E", "RHO", "S_OP", "S_BR", "S75", "E_CNT", "h_r", "h_a", "L", "La", "r_top", "g_a", "g_r", "dP_arr", "kc", "S_SP",
          "ARM_T_PER_MN", "W11", "sun", "n_sp", "n_st"]
STAGE2 = ["dP_g", "TAPER", "F_launch", "F1", "NQ", "r_mid", "dr", "r_edge", "EXPO_mid", "EXPO_edge"]
FUNCS = ["geff", "pot", "dia", "catenary", "parabola", "v2", "case11", "potv", "hang", "hang_xy"]

def load_round9(top_is_array=False, no_spin=False):
    """Round 9's functions and constants in a namespace of their own.  top_is_array: the line ends at the array's altitude (the kink) and
    not at the ring.  no_spin: the centrifugal term off (for the check of the size of that term)."""
    tree = ast.parse(open(R9, encoding="utf-8").read())
    ns = dict(np=np, math=math, brentq=brentq, root=root, R=R, cfg=cfg, phys=phys, P=lambda *a: None)
    def names(node):
        if isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) for t in node.targets): return [t.id for t in node.targets]
        return []
    def run(node): exec(compile(ast.Module(body=[node], type_ignores=[]), R9, "exec"), ns)
    stage2 = False
    for node in tree.body:                                   # in the order of the source, so that default arguments find their constants
        if isinstance(node, ast.FunctionDef) and node.name in FUNCS: run(node); continue
        nm = names(node)
        if nm and all(n in STAGE1 for n in nm):
            run(node)
            if no_spin and "W_E" in nm: ns["W_E"] = 0.0
        elif nm and all(n in STAGE2 for n in nm):
            if not stage2:                                   # the first constant of round 9's part 3: set the top of the line here
                stage2 = True; ns["h_ring"] = ns["h_r"]
                if top_is_array: ns["h_r"] = ns["h_a"]; ns["r_top"] = ns["RE"] + ns["h_a"]
            run(node)
    missing = [n for n in STAGE1 + STAGE2 + FUNCS if n not in ns]
    if missing: raise RuntimeError("not found in round9_lines_array.py: %s" % missing)
    return ns

ring = load_round9()                         # the reaction lines, as round 9 ran them
arm = load_round9(top_is_array=True)         # the arm below the kink: top at 600 km
arm0 = load_round9(top_is_array=True, no_spin=True)
hang = arm["hang"]; dia = arm["dia"]; RE = arm["RE"]; mu = arm["mu"]; W_E = arm["W_E"]; RHO = arm["RHO"]; S_OP = arm["S_OP"]; S75 = arm["S75"]
h_a = arm["h_a"]; h_ring = arm["h_ring"]; g_r = arm["g_r"]; g_a = arm["g_a"]; L = arm["L"]; La = arm["La"]; S_SP = arm["S_SP"]; r_k = arm["r_top"]
n_sp = arm["n_sp"]; n_st = arm["n_st"]

P("Round 13 item 5: the spreader arms below the kink as hanging lines.")
P("kink at %.0f km, ring at %.0f km ; CNT %.0f kg/m3 at %.2f GPa ; spreaders every %.0f km, %.1f pairs" % (h_a/1e3, h_ring/1e3, RHO, S_OP/1e9, S_SP/1e3, n_sp))
P("area taper, ground to the kink: %.4f (ground to the ring, round 9: %.4f) ; radius ratio at the kink %.4f" % (arm["TAPER"], ring["TAPER"], r_k/RE))

# ---------------------------------------------------------------------------------------------------------------------
# CHECK 0.  The loaded functions give round 9's logged reaction-line rows
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nCHECK 0.  The loaded method against round 9's log (reaction lines, rows A and B)\n" + "="*110)
log9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))["reaction"]["cases"]
F1 = ring["F1"]; gB = ring["h_r"]*math.tan(math.radians(30))/1e3
cA = ring["hang"](F1, 30.0)
bB = brentq(lambda b: ring["hang"](F1, b)["ground_km"] - gB, 15.0, 89.0); cB = ring["hang"](F1, bB)
OUT["check0"] = {}; worst = 0.0
for nm, c in (("A. 30 deg at the ground", cA), ("B. anchorage 404 km away (700 tan 30)", cB)):
    lg = log9[nm]
    for k in ("beta_g", "beta_top", "T_g", "T_top", "ground_km", "length_km", "mass_t", "sag_km", "d_g", "d_top"):
        worst = max(worst, abs(c[k]/lg[k] - 1))
    P("  %-40s here: ground %.2f deg, %.2f kN, %.3f mm ; casing %.2f deg ; anchorage %.2f km ; length %.2f km ; mass %.3f t ; sag %.2f km"
      % (nm, c["beta_g"], c["T_g"]/1e3, c["d_g"], c["beta_top"], c["ground_km"], c["length_km"], c["mass_t"], c["sag_km"]))
    P("  %-40s log:  ground %.2f deg, %.2f kN, %.3f mm ; casing %.2f deg ; anchorage %.2f km ; length %.2f km ; mass %.3f t ; sag %.2f km"
      % ("", lg["beta_g"], lg["T_g"]/1e3, lg["d_g"], lg["beta_top"], lg["ground_km"], lg["length_km"], lg["mass_t"], lg["sag_km"]))
    OUT["check0"][nm] = {k: c[k] for k in ("beta_g", "beta_top", "T_g", "T_top", "ground_km", "length_km", "mass_t", "sag_km")}
P("  largest relative difference over ten quantities and both rows: %.2e   -> %s" % (worst, "reproduced" if worst < 1e-9 else "NOT reproduced"))
OUT["check0"]["worst_rel_diff"] = worst

# ---------------------------------------------------------------------------------------------------------------------
# The book's arm: the straight line of spreader_cases.py and ch15_arm_choices.py
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nThe book's arm below the kink (straight line)\n" + "="*110)
d5 = arm["v2"](0.5)                           # ch16_numbers.py section Q as round 9 wrote it out: the design case
FR = d5["FR"]                                 # N, sideways pull on one arm at the kink
def ground_straight(lam_deg):
    """spreader_cases.py, function ground: straight line from the kink leaning lam from the local vertical; distance along the ground from the
    point under the kink (km), length (km), and the angle from the local vertical where it meets the ground (deg)."""
    Pk = (r_k, 0.0); d = (-math.cos(math.radians(lam_deg)), math.sin(math.radians(lam_deg)))
    b = 2*(Pk[0]*d[0] + Pk[1]*d[1]); c = Pk[0]**2 + Pk[1]**2 - RE**2
    t = (-b - math.sqrt(b*b - 4*c))/2; x = Pk[0] + t*d[0]; y = Pk[1] + t*d[1]
    return RE*math.atan2(y, x)/1e3, t/1e3, math.degrees(math.asin(math.sin(math.radians(lam_deg))*r_k/RE))
ARM_T_PER_MN = arm["ARM_T_PER_MN"]
book = {}
for lam in (30.0, 45.0, 60.0):
    gd, Ls, bg = ground_straight(lam); book[lam] = dict(anchorage_km=gd, length_km=Ls, beta_g=bg)
    P("  straight at %2.0f deg at the kink: anchorage %.1f km from the equator ; length %.1f km ; meets the ground at %.2f deg from the vertical" % (lam, gd, Ls, bg))
gd45, Ls45, bg45 = ground_straight(45.0)
T_low_b = FR/math.sin(math.radians(45)); T_up_b = FR/math.tan(math.radians(45))
A_g_b = (T_low_b/S_OP)/arm["TAPER"]
P("  design case (v2, 0.5 kg/m2): sheet %.1f m ; sideways pull on an arm %.2f MN ; arm below the kink %.2f MN, %.1f mm at the kink ; above the kink %.2f MN, %.1f mm"
  % (d5["W"], FR/1e6, T_low_b/1e6, dia(T_low_b, S_OP), T_up_b/1e6, dia(T_up_b, S_OP)))
P("  mass %.1f t per MN at the kink x %.2f MN = %.0f t ; pair every 10 km %.0f kg/m ; casing load 2 x %.2f MN / g at the ring = %.0f t per spreader = %.1f kg/m"
  % (ARM_T_PER_MN, T_low_b/1e6, T_low_b/1e6*ARM_T_PER_MN, 2*T_low_b/1e6*ARM_T_PER_MN*1e3/S_SP, T_up_b/1e6, 2*T_up_b/g_r/1e3, 2*T_up_b/g_r/S_SP))
P("  at the ground, if the 32 MN at the kink is tapered down by the height: %.1f mm (the book says about 40 mm)" % dia(A_g_b*S_OP, S_OP))
P("  What the straight line leaves out, as for the reaction lines: with %.1f MN at the kink a constant-stress line has %.1f MN at the ground." % (T_low_b/1e6, T_low_b/arm["TAPER"]/1e6))
P("  Its sideways pull at the ground would be %.1f MN x sin(%.2f deg) = %.1f MN, against %.1f MN x %.4f = %.1f MN that the moment about the Earth's center asks for."
  % (T_low_b/arm["TAPER"]/1e6, bg45, T_low_b/arm["TAPER"]*math.sin(math.radians(bg45))/1e6, FR/1e6, r_k/RE, FR*r_k/RE/1e6))
OUT["book"] = dict(FR_MN=FR/1e6, sheet_m=d5["W"], T_low_MN=T_low_b/1e6, T_up_MN=T_up_b/1e6, d_kink_mm=dia(T_low_b, S_OP), d_up_mm=dia(T_up_b, S_OP),
                   mass_t=T_low_b/1e6*ARM_T_PER_MN, pair_kgm=2*T_low_b/1e6*ARM_T_PER_MN*1e3/S_SP, casing_t=2*T_up_b/g_r/1e3, casing_kgm=2*T_up_b/g_r/S_SP,
                   straight={str(k): v for k, v in book.items()}, d_ground_tapered_mm=dia(A_g_b*S_OP, S_OP))

# ---------------------------------------------------------------------------------------------------------------------
# CHECK 1.  The limit: a very strong line gives the straight line
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nCHECK 1.  The limit of a very strong line (anchorage held at the straight line's %.1f km)\n" % gd45 + "="*110)
P("  stress      | ground angle | kink angle | tension at the kink | length    | sag     | mass      | mass / (rho x T_kink/stress x length)")
OUT["check1"] = {}
for s in (12.5e9, 25e9, 50e9, 125e9, 1.25e12, 1.25e14, 1.25e16):
    f = lambda b: hang(FR, b, stress=s)["ground_km"] - gd45
    b = brentq(f, 20.0, 89.5, xtol=1e-12); c = hang(FR, b, stress=s)
    ratio = c["mass_t"]*1e3/(RHO*c["T_top"]/s*c["length_km"]*1e3)
    P("  %9.4g GPa | %7.3f deg  | %6.3f deg | %8.3f MN         | %7.2f km | %6.3f km | %9.4g t | %.4f" % (s/1e9, c["beta_g"], c["beta_top"], c["T_top"]/1e6, c["length_km"], c["sag_km"], c["mass_t"], ratio))
    OUT["check1"]["%g" % (s/1e9)] = dict(beta_g=c["beta_g"], beta_top=c["beta_top"], T_top_MN=c["T_top"]/1e6, length_km=c["length_km"], sag_km=c["sag_km"], mass_t=c["mass_t"], mass_ratio=ratio)
lim = OUT["check1"]["%g" % 1.25e7]
P("  straight line (spreader_cases.py): ground %.3f deg ; kink 45.000 deg ; tension at the kink %.3f MN ; length %.2f km ; sag 0 ; anchorage %.1f km" % (bg45, T_low_b/1e6, Ls45, gd45))
ok1 = abs(lim["beta_top"] - 45) < 1e-3 and abs(lim["T_top_MN"]/(T_low_b/1e6) - 1) < 1e-4 and abs(lim["length_km"]/Ls45 - 1) < 1e-4 and abs(lim["beta_g"] - bg45) < 1e-3 and abs(lim["mass_ratio"] - 1) < 1e-3
P("  -> %s: as the stress grows the hanging line goes to 45 deg at the kink, %.2f MN, %.1f km, no sag, and its mass to density x area x length." % ("limit confirmed" if ok1 else "LIMIT NOT CONFIRMED", lim["T_top_MN"], lim["length_km"]))
OUT["check1"]["confirmed"] = bool(ok1)

# ---------------------------------------------------------------------------------------------------------------------
# The hanging arm
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nThe hanging arm at 12.5 GPa, for the same sideways pull at the kink (%.2f MN)\n" % (FR/1e6) + "="*110)
sb_max = 1.0/(arm["TAPER"]*r_k/RE); bt_max = math.degrees(math.asin(sb_max))
P("  The angle at the kink cannot be chosen freely: sin(kink angle) = sin(ground angle) x R_E / (r_kink x taper) <= 1 / (%.4f x %.4f), so the arm leans at most %.2f deg at the kink,"
  % (r_k/RE, arm["TAPER"], bt_max))
P("  and then it lies flat at the ground.  An arm at 45 deg at the kink does not exist at 12.5 GPa: sin(45 deg) x %.4f x %.4f = %.3f > 1." % (arm["TAPER"], r_k/RE, math.sin(math.radians(45))*arm["TAPER"]*r_k/RE))
s45 = brentq(lambda s: math.sin(math.radians(45))*math.exp(RHO*arm["potv"](r_k)/s)*r_k/RE - 1.0, 12.5e9, 1e12)
P("  For 45 deg at the kink to exist at all (arm flat at the ground) the working stress would have to be above %.1f GPa." % (s45/1e9))
OUT["limit_angle"] = dict(max_kink_angle_deg=bt_max, stress_for_45_GPa=s45/1e9)

def row(c):
    """The figures of one hanging arm, for the sideways pull c['F']."""
    V = c["V_top_kN"]*1e3
    return dict(beta_g=c["beta_g"], beta_top=c["beta_top"], T_g_MN=c["T_g"]/1e6, T_top_MN=c["T_top"]/1e6, d_g=c["d_g"], d_top=c["d_top"], anchorage_km=c["ground_km"],
                length_km=c["length_km"], mass_t=c["mass_t"], sag_km=c["sag_km"], sag_alt_km=c["sag_alt_km"], chord_deg=c["chord_from_vertical_at_casing"],
                V_MN=V/1e6, d_up=dia(V, S_OP), V_over_F=V/c["F"], T_over_F=c["T_top"]/c["F"], t_per_MN_kink=c["mass_t"]/(c["T_top"]/1e6), t_per_MN_pull=c["mass_t"]/(c["F"]/1e6),
                pair_kgm=2*c["mass_t"]*1e3/S_SP, casing_t=2*V/g_r/1e3, casing_kgm=2*V/g_r/S_SP)
def show(nm, r):
    P("  %-46s| ground %5.2f deg, %6.1f MN, %5.1f mm | kink %5.2f deg, %6.1f MN, %5.1f mm | above the kink %6.1f MN, %5.1f mm | anchorage %7.1f km | length %7.1f km | sag %5.1f km at %3.0f km up | mass %6.0f t (%5.1f t per MN at the kink) | pair %5.0f kg/m | casing %6.0f t per spreader, %6.1f kg/m"
      % (nm, r["beta_g"], r["T_g_MN"], r["d_g"], r["beta_top"], r["T_top_MN"], r["d_top"], r["V_MN"], r["d_up"], r["anchorage_km"], r["length_km"], r["sag_km"], r["sag_alt_km"], r["mass_t"], r["t_per_MN_kink"], r["pair_kgm"], r["casing_t"], r["casing_kgm"]))
def by_anchorage(km): return brentq(lambda b: hang(FR, b)["ground_km"] - km, 5.0, 89.9, xtol=1e-12)
cases = {}
cases["A. 45 deg at the ground"] = hang(FR, 45.0)
cases["A2. %.1f deg at the ground (the straight line's)" % bg45] = hang(FR, bg45)
bBarm = by_anchorage(gd45); cases["B. anchorage %.0f km from the equator" % gd45] = hang(FR, bBarm)
bB2 = brentq(lambda b: hang(FR, b)["chord_from_vertical_at_casing"] - 45.0, 20.0, 89.9, xtol=1e-12); cases["B2. chord 45 deg from the kink's vertical"] = hang(FR, bB2)
far = hang(FR, 89.99)["ground_km"]
P("  The farthest anchorage a hanging arm can reach at 12.5 GPa (leaving the ground at 89.99 deg): %.0f km from the equator." % far)
OUT["limit_angle"]["farthest_anchorage_km"] = far
lean_names = {}
for lam in (30.0, 60.0):
    km = book[lam]["anchorage_km"]
    if km < far:
        lean_names[lam] = "anchorage %.0f km (the straight line's at %.0f deg)" % (km, lam); cases[lean_names[lam]] = hang(FR, by_anchorage(km))
    else:
        P("  The straight line's anchorage at %.0f deg, %.0f km, is out of reach: no hanging arm at 12.5 GPa gets there." % (lam, km))
for bt in (15.0, 20.0, 25.0):
    cases["%2.0f deg at the kink" % bt] = hang(FR, math.degrees(math.asin(math.sin(math.radians(bt))*arm["TAPER"]*r_k/RE)))
for bg in (60.0, 70.0, 80.0, 89.0):
    cases["%2.0f deg at the ground" % bg] = hang(FR, bg)
OUT["cases"] = {}
for nm, c in cases.items():
    OUT["cases"][nm] = row(c); show(nm, OUT["cases"][nm])
P("  %-46s| ground %5.2f deg, %6.1f MN, %5.1f mm | kink 45.00 deg, %6.1f MN, %5.1f mm | above the kink %6.1f MN, %5.1f mm | anchorage %7.1f km | length %7.1f km | sag   0.0 km               | mass %6.0f t (%5.1f t per MN at the kink) | pair %5.0f kg/m | casing %6.0f t per spreader, %6.1f kg/m"
  % ("book (straight line)", bg45, T_low_b/arm["TAPER"]/1e6, dia(A_g_b*S_OP, S_OP), T_low_b/1e6, dia(T_low_b, S_OP), T_up_b/1e6, dia(T_up_b, S_OP), gd45, Ls45, T_low_b/1e6*ARM_T_PER_MN, ARM_T_PER_MN,
     OUT["book"]["pair_kgm"], OUT["book"]["casing_t"], OUT["book"]["casing_kgm"]))
keyB = "B. anchorage %.0f km from the equator" % gd45; rB = OUT["cases"][keyB]; cBarm = cases[keyB]
P("\n  Row B against the straight line, same sideways pull %.2f MN and same anchorage %.0f km:" % (FR/1e6, gd45))
P("    tension at the kink %.1f MN, not %.1f (x %.2f) ; across at the kink %.0f mm, not %.0f ; pull carried up to the casing %.1f MN, not %.1f (x %.2f) ; mass %.0f t, not %.0f (x %.2f)"
  % (rB["T_top_MN"], T_low_b/1e6, rB["T_top_MN"]/(T_low_b/1e6), rB["d_top"], dia(T_low_b, S_OP), rB["V_MN"], T_up_b/1e6, rB["V_MN"]/(T_up_b/1e6), rB["mass_t"], T_low_b/1e6*ARM_T_PER_MN, rB["mass_t"]/(T_low_b/1e6*ARM_T_PER_MN)))

# ---- how the casing load and the mass go with the anchorage ----
P("\n  The arm for every anchorage (sideways pull %.2f MN):" % (FR/1e6))
P("  ground angle | anchorage  | kink angle | N at the kink per N sideways | N on the casing per N sideways | mass     | casing load per spreader")
OUT["scan"] = []
for bg in (30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 88, 89.5):
    c = hang(FR, float(bg)); r = row(c)
    P("  %5.1f deg    | %7.1f km | %5.2f deg  | %.3f                        | %.3f                          | %6.0f t | %6.0f t" % (bg, r["anchorage_km"], r["beta_top"], r["T_over_F"], r["V_over_F"], r["mass_t"], r["casing_t"]))
    OUT["scan"].append(dict(beta_g=bg, anchorage_km=r["anchorage_km"], beta_top=r["beta_top"], T_over_F=r["T_over_F"], V_over_F=r["V_over_F"], mass_t=r["mass_t"], casing_t=r["casing_t"]))
P("  The straight line has 1.414 N at the kink and 1.000 N on the casing for each N sideways.  No hanging arm at 12.5 GPa gets below %.3f and %.3f." % (1/sb_max, math.sqrt(1 - sb_max**2)/sb_max))

# ---- the same arm at other working stresses ----
P("\n  The arm to the %.0f km anchorage at other working stresses (sideways pull %.2f MN):" % (gd45, FR/1e6))
OUT["stress"] = {}
for s in (12.5e9, 15e9, 18.75e9, 25e9):
    b = brentq(lambda bb: hang(FR, bb, stress=s)["ground_km"] - gd45, 20.0, 89.5, xtol=1e-12); c = hang(FR, b, stress=s); V = c["V_top_kN"]*1e3
    P("    %.2f GPa: ground %.2f deg ; kink %.2f deg, %.1f MN, %.1f mm ; on the casing %.1f MN (%.0f t per spreader, %.0f kg/m) ; mass %.0f t ; sag %.1f km"
      % (s/1e9, c["beta_g"], c["beta_top"], c["T_top"]/1e6, c["d_top"], V/1e6, 2*V/g_r/1e3, 2*V/g_r/S_SP, c["mass_t"], c["sag_km"]))
    OUT["stress"]["%g" % (s/1e9)] = dict(beta_g=c["beta_g"], beta_top=c["beta_top"], T_top_MN=c["T_top"]/1e6, d_top=c["d_top"], V_MN=V/1e6, casing_t=2*V/g_r/1e3, casing_kgm=2*V/g_r/S_SP, mass_t=c["mass_t"], sag_km=c["sag_km"])
P("    (18.75 GPa is the stress the book gives the retention lines and the array support lines ; 25 GPa is the yarn's strength, no margin.)")

# ---------------------------------------------------------------------------------------------------------------------
# CHECK 2.  Independent integration in the meridian plane
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nCHECK 2.  Row B stepped through in the plane of the meridian (50 m steps)\n" + "="*110)
def step_meridian(F, beta_g_deg, phi_g, w, ds=50.0):
    """Independent of hang(): the vector equation d(T t)/ds = -rho A b stepped from the ground anchorage to the kink's radius.
    Axes: y from the Earth's center through the point on the equator under the kink, x to the north (along the spin axis).
    Body force per kg: b = -GM p/r^3 + w^2 (0, y), the centrifugal term pointing away from the spin axis (not along the radius, off the equator).
    The anchorage is at latitude phi_g ; the line leaves it at beta_g from the local vertical, leaning toward the equator."""
    K = F*r_k; T_g = K/(RE*math.sin(math.radians(beta_g_deg))); b = math.radians(beta_g_deg)
    pos = np.array([RE*math.sin(phi_g), RE*math.cos(phi_g)]); up = pos/RE; toeq = np.array([-up[1], up[0]])
    q = T_g*(math.cos(b)*up + math.sin(b)*toeq)
    def f(y):
        p = y[:2]; qq = y[2:]; r = math.hypot(p[0], p[1]); T = math.hypot(qq[0], qq[1])
        body = -mu*p/r**3 + np.array([0.0, w*w*p[1]])
        return np.concatenate((qq/T, -RHO*(T/S_OP)*body))
    y = np.concatenate((pos, q)); s = 0.0; m = 0.0
    while True:
        k1 = f(y); k2 = f(y + 0.5*ds*k1); k3 = f(y + 0.5*ds*k2); k4 = f(y + ds*k3); yn = y + ds/6*(k1 + 2*k2 + 2*k3 + k4)
        rn = math.hypot(yn[0], yn[1]); ro = math.hypot(y[0], y[1])
        if rn >= r_k:
            fr = (r_k - ro)/(rn - ro); m += RHO*math.hypot(y[2], y[3])/S_OP*ds*fr; y = y + fr*(yn - y); s += fr*ds; break
        m += RHO*0.5*(math.hypot(y[2], y[3]) + math.hypot(yn[2], yn[3]))/S_OP*ds; y = yn; s += ds
    T = math.hypot(y[2], y[3]); r = math.hypot(y[0], y[1]); rh = y[:2]/r
    cosb = float(np.dot(y[2:], rh))/T; side = float(np.dot(y[2:], np.array([-rh[1], rh[0]])))
    return dict(T_top_MN=T/1e6, beta_top=math.degrees(math.acos(cosb)), miss_km=RE*math.atan2(y[0], y[1])/1e3, length_km=s/1e3, mass_t=m/1e3, F_side_MN=side/1e6, V_MN=T*cosb/1e6)
c0 = arm0["hang"](FR, bBarm)                                       # hang() with the spin off, same ground angle
st0 = step_meridian(FR, bBarm, c0["phi_top"], 0.0)
st1 = step_meridian(FR, bBarm, cBarm["phi_top"], W_E)
P("  no spin:   hang()   kink %.4f deg, %.4f MN ; length %.3f km ; mass %.2f t ; sideways %.4f MN ; anchorage %.2f km" % (c0["beta_top"], c0["T_top"]/1e6, c0["length_km"], c0["mass_t"], FR/1e6, c0["ground_km"]))
P("             stepped  kink %.4f deg, %.4f MN ; length %.3f km ; mass %.2f t ; sideways %.4f MN ; ends %.3f km from the kink's place" % (st0["beta_top"], st0["T_top_MN"], st0["length_km"], st0["mass_t"], st0["F_side_MN"], st0["miss_km"]))
P("  with spin: hang()   kink %.4f deg, %.4f MN ; length %.3f km ; mass %.2f t ; sideways %.4f MN   (the plane taken as the equator's)" % (cBarm["beta_top"], cBarm["T_top"]/1e6, cBarm["length_km"], cBarm["mass_t"], FR/1e6))
P("             stepped  kink %.4f deg, %.4f MN ; length %.3f km ; mass %.2f t ; sideways %.4f MN ; ends %.3f km from the kink's place   (the centrifugal term as it is in the meridian plane)"
  % (st1["beta_top"], st1["T_top_MN"], st1["length_km"], st1["mass_t"], st1["F_side_MN"], st1["miss_km"]))
d0 = max(abs(st0["T_top_MN"]/(c0["T_top"]/1e6) - 1), abs(st0["mass_t"]/c0["mass_t"] - 1), abs(st0["length_km"]/c0["length_km"] - 1), abs(st0["F_side_MN"]/(FR/1e6) - 1))
d1 = max(abs(st1["T_top_MN"]/(cBarm["T_top"]/1e6) - 1), abs(st1["mass_t"]/cBarm["mass_t"] - 1), abs(st1["length_km"]/cBarm["length_km"] - 1), abs(st1["F_side_MN"]/(FR/1e6) - 1))
P("  largest relative difference, no spin: %.1e (the method) ; with spin: %.1e (what taking the meridian plane for the equator's plane costs)" % (d0, d1))
P("  the spin itself: tension at the kink %.4f MN with it, %.4f MN without (%.2f %%)" % (cBarm["T_top"]/1e6, c0["T_top"]/1e6, (cBarm["T_top"]/c0["T_top"] - 1)*100))
OUT["check2"] = dict(no_spin_hang=dict(beta_top=c0["beta_top"], T_top_MN=c0["T_top"]/1e6, length_km=c0["length_km"], mass_t=c0["mass_t"]), no_spin_stepped=st0, spin_stepped=st1, diff_no_spin=d0, diff_spin=d1)

# ---------------------------------------------------------------------------------------------------------------------
# What follows in Tables 11.1 and 11.2 and in the budget
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nTable 11.1 with the hanging arm of row B (anchorage %.0f km)\n" % gd45 + "="*110)
kT = rB["T_over_F"]; kV = rB["V_over_F"]; kM = rB["t_per_MN_pull"]
P("  Every figure of the arm is in proportion to the sideways pull.  For each MN of sideways pull, row B: %.3f MN at the kink (straight: 1.414), %.3f MN on the casing (straight: 1.000), %.1f t of arm (straight: 92.6)."
  % (kT, kV, kM))
P("  case                               | pull     | arm at the kink: book, hanging | arm mass: book, hanging | above the kink: book, hanging | ground: hanging")
OUT["table_11_1"] = {}
for nm, a, th, bw in (("Design case (0.5, 45 deg, W/8)", 0.5, 45, 1/8), ("Lighter panels (0.2, 45 deg, W/8)", 0.2, 45, 1/8), ("Steeper sheet (0.5, 60 deg, W/8)", 0.5, 60, 1/8),
                      ("Deeper bow (0.5, 45 deg, W/4)", 0.5, 45, 1/4), ("All three (0.2, 60 deg, W/4)", 0.2, 60, 1/4)):
    k = arm["case11"](nm, a, th, bw); pull = k["pull_MN"]*1e6
    r = dict(pull_MN=pull/1e6, book_mm=k["arm_mm"], book_t=k["arm_t"], hang_mm=dia(kT*pull, S_OP), hang_t=kM*pull/1e6, book_up_mm=dia(pull, S_OP), hang_up_mm=dia(kV*pull, S_OP),
             hang_ground_mm=dia(kT*pull/arm["TAPER"], S_OP), hang_kink_MN=kT*pull/1e6, hang_up_MN=kV*pull/1e6, ret_MN=k["TR_MN"], ret_mm=k["ret_mm"], ret_kgm=k["ret_kgm"])
    OUT["table_11_1"][nm] = r
    P("  %-35s| %5.2f MN | %4.1f mm, %4.1f mm (%.1f MN)     | %5.0f t, %5.0f t        | %4.1f mm, %4.1f mm (%.1f MN)    | %4.1f mm"
      % (nm, r["pull_MN"], r["book_mm"], r["hang_mm"], r["hang_kink_MN"], r["book_t"], r["hang_t"], r["book_up_mm"], r["hang_up_mm"], r["hang_up_MN"], r["hang_ground_mm"]))
P("  (Table 11.1 uses a sheet of 1,110 m ; the retention line columns and the last column do not depend on the arm.)")

P("\n" + "="*110 + "\nThe 1,200 kg/m budget\n" + "="*110)
P("  The book counts the arms on the casing as the pull they carry up to it, divided by g at the ring: 2 x V / (g x 10 km).  V is the downward pull of the lower arm at the kink.")
P("  The weight of the 100 km of arm above the kink is not in that count, here or in the book.")
tot_keep = d5["pan"] + rB["casing_kgm"] + d5["mR"]
P("  (i) the sheet kept at %.1f m: panels and lines %.1f + arms %.1f + retention %.1f = %.1f kg/m (the book: %.1f + %.1f + %.1f = %.1f)"
  % (d5["W"], d5["pan"], rB["casing_kgm"], d5["mR"], tot_keep, d5["pan"], d5["arms"], d5["mR"], d5["pan"] + d5["arms"] + d5["mR"]))
OUT["budget"] = dict(keep_sheet=dict(sheet_m=d5["W"], pan=d5["pan"], arms=rB["casing_kgm"], ret=d5["mR"], total=tot_keep))
p9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))["power"]["catenary"]["rows"]["0.5"]
P("  (ii) the 1,200 kg/m held, the strip made narrower.  v2 of ch16_numbers.py (as round 9 wrote it out) takes the arm's lean only as tan and sin of the angle at the kink,")
P("       so it is run with the hanging arm's angle at the kink in place of 45 deg.  The anchorage is held, so the angle at the kink is the same whatever the pull.")
OUT["budget"]["hold_1200"] = {}
for nm in [keyB] + list(lean_names.values()):
    r = OUT["cases"][nm]; d = arm["v2"](0.5, lam=r["beta_top"]); c = hang(d["FR"], r["beta_g"]); rr = row(c)
    Pw = p9["P_TW"]*d["W"]/p9["sheet_km"]/1e3
    P("       %-46s kink %.2f deg: sheet %.3f km ; panels and lines %.1f + arms %.1f + retention %.1f = %.1f kg/m ; sideways pull %.2f MN ; kink %.1f MN, %.1f mm ; above %.1f MN, %.1f mm ; ground %.1f mm ; arm %.0f t, pair %.0f kg/m ; casing %.0f t per spreader ; power %.2f TW"
      % (nm, r["beta_top"], d["W"]/1e3, d["pan"], d["arms"], d["mR"], d["pan"] + d["arms"] + d["mR"], d["FR"]/1e6, rr["T_top_MN"], rr["d_top"], rr["V_MN"], rr["d_up"], rr["d_g"], rr["mass_t"], rr["pair_kgm"], rr["casing_t"], Pw))
    OUT["budget"]["hold_1200"][nm] = dict(beta_top=r["beta_top"], sheet_km=d["W"]/1e3, pan=d["pan"], arms=d["arms"], ret=d["mR"], FR_MN=d["FR"]/1e6, T_kink_MN=rr["T_top_MN"], d_kink=rr["d_top"],
                                          V_MN=rr["V_MN"], d_up=rr["d_up"], d_g=rr["d_g"], mass_t=rr["mass_t"], pair_kgm=rr["pair_kgm"], casing_t=rr["casing_t"], casing_kgm=rr["casing_kgm"], P_TW=Pw,
                                          anchorage_km=rr["anchorage_km"], length_km=rr["length_km"])
P("       power: round 9's catenary row (%.2f TW for a sheet of %.3f km, runs/round9_logs/round9_lines_array.json) times the ratio of the sheet widths ; the power is in proportion to the width in ch16_numbers.py section R."
  % (p9["P_TW"], p9["sheet_km"]))
P("       book, straight arms: 0.82 km and 7,300 t at 30 deg (352 km) ; 1.11 km and 5,700 t at 45 deg (632 km) ; 1.40 km and 4,200 t at 60 deg (1,264 km).")

# ---- the swing along the ring, with the book's own formula ----
P("\n" + "="*110 + "\nThe swing along the ring (Section 11.6), the book's formula with the hanging arm's tensions\n" + "="*110)
P("  ch15_spreader_update.py: each arm is two pendulums in series, tension / length: the 100 km above the kink and the arm below it.  That formula is kept.")
P("  [UNVERIFIED] for a sagging arm: tension / length is the straight line's sideways stiffness, and the tension is taken at the kink.")
M_cas = 12000.0*L; k_re = 1.55            # N/m, round 9 row B, runs/round9_logs/round9_lines_array.json
def swing(T_up, L_up, T_low, L_low):
    k_arm = 1/(1/(T_up/L_up) + 1/(T_low/L_low)); kt = 2*n_sp*k_arm + n_st*k_re; om = math.sqrt(kt/M_cas)
    return dict(k_up=T_up/L_up, k_low=T_low/L_low, k_arm=k_arm, k_ring=kt, period_s=2*math.pi/om, swing_m=0.576/om, swing_braked_m=0.281/om)
OUT["swing"] = {}
for nm, a in (("book (23 MN, 100 km ; 32 MN, 893 km)", (23e6, 100e3, 32e6, 893e3)),
              ("hanging arm, row B, sheet kept", (rB["V_MN"]*1e6, 100e3, rB["T_top_MN"]*1e6, rB["length_km"]*1e3)),
              ("hanging arm, row B, 1,200 kg/m held", (OUT["budget"]["hold_1200"][keyB]["V_MN"]*1e6, 100e3, OUT["budget"]["hold_1200"][keyB]["T_kink_MN"]*1e6, OUT["budget"]["hold_1200"][keyB]["length_km"]*1e3))):
    s = swing(*a); OUT["swing"][nm] = s
    P("  %-38s upper arm %.0f N/m, lower arm %.1f N/m, in series %.1f N/m ; ring held by %.0f N/m ; period %.0f s (%.2f h) ; swing after a launch %.0f m, with the sled braked %.0f m"
      % (nm, s["k_up"], s["k_low"], s["k_arm"], s["k_ring"], s["period_s"], s["period_s"]/3600, s["swing_m"], s["swing_braked_m"]))

# ---- the path of row B and of the reaction line, for the figures ----
def path_xy(ns, F, bg, n=60):
    c = ns["hang"](F, bg, want_path=True); r_e, ph_e, s_e, T_e = c["path"]; idx = np.linspace(0, len(r_e) - 1, n).astype(int)
    return dict(alt_km=[float((r_e[i] - ns["RE"])/1e3) for i in idx], along_km=[float(ns["RE"]*(c["phi_top"] - ph_e[i])/1e3) for i in idx])
OUT["path_arm_B"] = path_xy(arm, FR, bBarm); OUT["path_reaction_B"] = path_xy(ring, F1, bB)
OUT["reaction_B"] = dict(beta_g=cB["beta_g"], beta_top=cB["beta_top"], anchorage_km=cB["ground_km"], sag_km=cB["sag_km"], sag_alt_km=cB["sag_alt_km"], length_km=cB["length_km"])
P("\npaths of the reaction line (row B of round 9) and of the arm (row B here) saved in the json for the figures: altitude against distance along the ground from the point under the top.")

with open(os.path.join(LOG_DIR, "round13_spreader_arms.txt"), "w", encoding="utf-8") as f: f.write(_buf.getvalue())
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    return o
with open(os.path.join(LOG_DIR, "round13_spreader_arms.json"), "w", encoding="utf-8") as f: json.dump(clean(OUT), f, indent=1)
print("\nlog written: runs/round13_logs/round13_spreader_arms.txt and .json")
