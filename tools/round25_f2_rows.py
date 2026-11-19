#!/usr/bin/env python3
"""
Round 25: the rows of Section 11.6 that round 24 did not run for layout F2 (guide line every 100 km).
Brief: Vol_III/Claude_Code_instructions_round25_2026-11-19.md.   Book: Orbital Ring Engineering, Volume III, Chapter 11, Section 11.6.

  1. The design launch with nothing to damp it, sled not braked and braked.
  2. The crewed launch at the 3 g limit and the return of ten empty pods, damped, each with and without braking.
  3. Traffic: seven (and fourteen) unbraked design launches a day ; the pull in the tighter and the slacker line of a pair at the offset.
  4. The lean of the lines and the rotation of the casing after the design launch.

Layout F2 is read from round 24's log (runs/round24_logs/round24_layout_f.json) and not changed.  The launches, the braking and the casing's
integrator are round 20's (tools/round20_brake_felt_pods.py, its prefix run as it stands) ; the crewed launch is round 21's call.  The lines
are round 9's model (tools/round9_lines_array.py: hang and spring), taken the way round 23 takes them.

Run from Python/orbitalring:  python tools/round25_f2_rows.py
Log: runs/round25_logs/round25_f2_rows.txt and .json ; the simulation's own printing for each run: runs/round25_logs/lsm_<case>.txt
"""
import ast, contextlib, io, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root
import ring_altitude as R

T_START = time.time()
LOG_DIR = os.path.join(ROOT, "runs", "round25_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
OUT = {}
def jl(*p): return json.load(open(os.path.join(ROOT, "runs", *p), encoding="utf-8"))
L20 = jl("round20_logs", "round20_brake_felt_pods.json"); L21 = jl("round21_logs", "round21_crew_limit.json"); L24 = jl("round24_logs", "round24_layout_f.json"); L23 = jl("round23_logs", "round23_reaction_vs_guide.json")

# ---------------------------------------------------------------------------------------------------------------------
# layout F2, guide line every 100 km, from round 24's log
Y = next(y for y in L24["layouts"] if y["key"] == "F2-100")
K = Y["hold"]; s = Y["s"]; k_line = Y["k_line"]; F_BUILT = Y["F_kN"]*1e3; n_lines = Y["n_react_lines"]; n_pairs = n_lines/2; n_guide = Y["n_guide_lines"]; k_g = L23["guide"]["k_p_round"]
assert abs(n_lines*k_line + n_guide*k_g - K) < 1e-6 and abs(Y["kgm"] - L24["KGM_BOOK"]) < 1e-9

# round 20's launch, braking and integrator
src20 = open(os.path.join(HERE, "round20_brake_felt_pods.py"), encoding="utf-8").read(); cut = src20.index('P("Round 20.'); assert cut > 0
g20 = {"__name__": "round20_prefix", "__file__": os.path.join(HERE, "round20_brake_felt_pods.py")}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src20[:cut], "round20_brake_felt_pods.py (prefix)", "exec"), g20)
g20["LOG_DIR"] = LOG_DIR
lsm_run = g20["lsm_run"]; cfg = g20["cfg"]; e1 = g20["e1"]; e2 = g20["e2"]
M_CAS = g20["M_CAS"]; F_LIM = g20["F_LIM"]; Lr = g20["Lr"]; rebin = g20["rebin"]; swing_once = g20["swing_once"]; swing_daily = g20["swing_daily"]; undamped = g20["undamped"]
K_OLD = g20["K"]
C = 2*math.sqrt(K*M_CAS); OM = math.sqrt(K/M_CAS); PERIOD = 2*math.pi/OM
assert abs(C - Y["c_crit"]) < 1e-3*C and abs(PERIOD - Y["period_s"]) < 1e-6

P("Round 25.  The rows of Section 11.6 that round 24 did not run for layout F2.")
P("Layout F2 from round 24's log: %.1f reaction lines in opposed pairs every 25 km, each %.4f times the design line, built for %.2f kN along the ring, %.4f N/m ; %.1f guide lines at %.4f N/m ;"
  % (n_lines, s, F_BUILT/1e3, k_line, n_guide, k_g))
P("hold %.1f N/m ; casing %.4e kg ; period %.0f s (%.2f h) ; LIMs damping at the critical rate %.3e N s/m, limited to %.1f MN.  (The old layout's hold in round 20: %.1f N/m.)" % (K, M_CAS, PERIOD, PERIOD/3600, C, F_LIM/1e6, K_OLD))

def runs_of(c):
    """The two pushes of one launch on the casing, in 10 s steps: sled not braked, and braked (the braking pushes the casing the other way)."""
    Fl = c["F"]; Fb = c["FB"]
    return {"not braked": (rebin(Fl, 0.1, 10.0), float(Fl.sum()*0.1), len(Fl)*0.1), "braked": (rebin(np.concatenate((Fl, -Fb)), 0.1, 10.0), float(Fl.sum()*0.1 - Fb.sum()*0.1), (len(Fl) + len(Fb))*0.1)}

def once_rows(c, tag):
    res = {}
    for nm, (F10, J, dur) in runs_of(c).items():
        o = swing_once(F10, 10.0, C, k=K); u = undamped(F10, 10.0, k=K)
        o.update(push_s=dur, impulse=J, undamped_max=u["x_max"], undamped_after=u["amp_after"], undamped_during=u["x_max_during"], sudden=J/M_CAS/OM, line_dF_kN=k_line*o["x_max"]/1e3,
                 lean_deg=math.degrees(math.atan(o["x_max"]/R.ALTITUDE)), rot_deg=360*o["x_max"]/Lr)
        res[nm] = o
        P("   %-10s: push lasts %.0f s (%.1f min), net impulse %.4e N s, %.4f m/s on the casing.  LIMs damping: farthest %.0f m, %.2f h after the launch begins ; within 1 %% of that after %.1f h ; LIMs' peak push %.2f MN of %.1f."
          % (nm, dur, dur/60, J, J/M_CAS, o["x_max"], o["peak_h"], o["die_h"], o["F_lim_MN"], F_LIM/1e6))
        P("   %s  No damping: farthest %.0f m, left swinging %.0f m (a sudden push of the same impulse: %.0f m).  A reaction line's pull changes by %.2f kN, %.2f %% of %.1f kN."
          % (" "*10, u["x_max"], u["amp_after"], o["sudden"], o["line_dF_kN"], 100*o["line_dF_kN"]*1e3/F_BUILT, F_BUILT/1e3))
    OUT[tag] = res
    return res

# =====================================================================================================================
P("\n" + "="*118 + "\n1. The design launch (lsm_simulation.py --brake), with nothing to damp it and with the LIMs damping\n" + "="*118)
with contextlib.redirect_stdout(io.StringIO()): D = lsm_run("design_brake")
d20 = L20["design"]
assert abs(D["dur"] - d20["dur"]) < 1e-6 and abs(D["impulse"] - d20["impulse"]) < 1.0 and abs(D["brake"]["duration_s"] - d20["brake"]["duration_s"]) < 1e-6 and abs(D["brake"]["impulse"] - d20["brake"]["impulse"]) < 1.0
P("   The launch: %.1f minutes, peak %.2f MN, impulse %.4e N s ; the braking %.1f minutes from the release, impulse given back %.4e N s.  The same as round 20's log to the last digit."
  % (D["dur"]/60, D["F_peak"]/1e6, D["impulse"], D["brake"]["duration_s"]/60, D["brake"]["impulse"]))
rD = once_rows(D, "design")
o24 = Y["once"]
assert abs(rD["not braked"]["x_max"] - o24["not braked"]["x_max"]) < 1e-6 and abs(rD["braked"]["x_max"] - o24["braked"]["x_max"]) < 1e-6
P("   The damped figures are round 24's (%.0f m and %.0f m), reproduced." % (o24["not braked"]["x_max"], o24["braked"]["x_max"]))
P("   Sudden-push arithmetic as the chapter prints it: %.4f m/s x %.0f s / 2 pi = %.0f m." % (D["impulse"]/M_CAS, PERIOD, D["impulse"]/M_CAS*PERIOD/(2*math.pi)))
c20 = d20["casing"]
P("   The old layout (round 20, hold %.1f N/m): no damping %.0f m not braked, %.0f m braked." % (K_OLD, c20["not braked"]["undamped_after"], c20["braked"]["undamped_after"]))

# =====================================================================================================================
P("\n" + "="*118 + "\n2. The crewed launch at the 3 g limit and the return of ten empty pods, LIMs damping\n" + "="*118)
Rr = cfg.R_ORBIT; v_rot = cfg.V_CASING_SPACE
v3_space = math.sqrt((3*cfg.G_0 + cfg.G_250)*Rr); v3_cas = v3_space - v_rot
with contextlib.redirect_stdout(io.StringIO()): Cw = lsm_run("crewed", V_LAUNCH=v3_cas, G_LOAD_MAX=3.0)
c21 = L21["crewed"]
assert abs(Cw["dur"] - c21["dur"]) < 1e-6 and abs(Cw["impulse"] - c21["impulse"]) < 1.0 and abs(Cw["brake"]["impulse"] - c21["brake"]["impulse"]) < 1.0
P("   The crewed launch (round 21's call: --v_launch=%.1f --g-limit=3 --brake): %.2f minutes, release %.1f m/s in space, net outward %.4f g, peak thrust %.2f MN, impulse %.4e N s ;"
  % (v3_cas, Cw["dur"]/60, Cw["v_space"], Cw["a_out_g"], Cw["F_peak"]/1e6, Cw["impulse"]))
P("   braking of its sled %.2f minutes, impulse given back %.4e N s.  The same as round 21's log." % (Cw["brake"]["duration_s"]/60, Cw["brake"]["impulse"]))
rC = once_rows(Cw, "crewed")
k21 = c21["casing"]
P("   The old layout (round 21): farthest %.0f m not braked, %.0f m braked." % (k21["not braked"]["x_max"], k21["braked"]["x_max"]))

# the empty pod, as round 20 has it
rates = e1.load_rates(); T_sh = rates["T"]
d4 = e1.descent(4000.0, 1000.0); _p, _e, tot = e1.dose_per_descent(rates, d4); x_sh = e1.shield_for(T_sh, tot, 2.0)
rows = e2.arrangements(); mag = dict(min([r for r in rows if r["kind"] == "facing pair"], key=lambda r: r["kg_per_N"]))
gp = e2.size_pod2(100, x_sh, mag); m_empty = gp["m"] - gp["m_pax"]
assert abs(m_empty - L20["pod"]["m_empty"]) < 1e-6
r_a = e1.r_geo; v_h = math.sqrt(e1.GM*(2/Rr - 2/(Rr + r_a))); v_h_cas = v_h - v_rot
with contextlib.redirect_stdout(io.StringIO()): Pd = lsm_run("pod_10", M_SPACECRAFT=10*m_empty, V_LAUNCH=v_h_cas, A_MAX_G=2.0, A_BRAKE_MAX_G=None)
p20 = L20["pod_10"]
assert abs(Pd["dur"] - p20["dur"]) < 1e-6 and abs(Pd["impulse"] - p20["impulse"]) < 1.0 and abs(Pd["brake"]["impulse"] - p20["brake"]["impulse"]) < 1.0
P("\n   Ten empty pods (round 20's call: %.2f t each, %.1f t on the 5,250 t sled, --v_launch=%.1f --accel=2 --brake): %.2f minutes, peak thrust %.2f MN, impulse %.4e N s ;"
  % (m_empty/1e3, 10*m_empty/1e3, v_h_cas, Pd["dur"]/60, Pd["F_peak"]/1e6, Pd["impulse"]))
P("   braking of the sled %.2f minutes, impulse given back %.4e N s.  The same as round 20's log." % (Pd["brake"]["duration_s"]/60, Pd["brake"]["impulse"]))
rP = once_rows(Pd, "pod_10")
kp = p20["casing"]
P("   The old layout (round 20): farthest %.0f m not braked, %.0f m braked." % (kp["not braked"]["x_max"], kp["braked"]["x_max"]))
OUT["launches"] = dict(crewed=dict(dur=Cw["dur"], v_space=Cw["v_space"], a_out_g=Cw["a_out_g"], impulse=Cw["impulse"], brake_s=Cw["brake"]["duration_s"], brake_impulse=Cw["brake"]["impulse"], v3_cas=v3_cas),
                       pod_10=dict(dur=Pd["dur"], m_empty=m_empty, impulse=Pd["impulse"], brake_s=Pd["brake"]["duration_s"], brake_impulse=Pd["brake"]["impulse"], v_h_cas=v_h_cas),
                       design=dict(dur=D["dur"], impulse=D["impulse"], brake_s=D["brake"]["duration_s"], brake_impulse=D["brake"]["impulse"]))
P("  (%.0f s so far)" % (time.time() - T_START))

# =====================================================================================================================
P("\n" + "="*118 + "\n3. Traffic\n" + "="*118)
RD = runs_of(D); daily = {}
P("   launches | sled       | mean push  | under way at once | the casing sits at | swings about it | farthest | LIMs' push: peak, mean")
for N in (7, 14):
    for nm in ("braked", "not braked"):
        F10, J, dur = RD[nm]
        o = swing_daily(F10, 10.0, N, C, False, k=K, J_exact=J); daily[(N, nm)] = o
        P("   %2d a day | %-10s | %6.3f MN  | %4.2f              | %7.1f m          | %6.1f m        | %6.0f m | %5.1f, %5.1f MN" % (N, nm, o["mean_push_MN"], o["overlap"], o["x_mean"], o["swing"], o["x_max"], o["F_lim_peak_MN"], o["F_lim_mean_MN"]))
        assert abs(o["x_mean"] - o["mean_push_MN"]*1e6/K) < 0.01*o["x_mean"]
for N in (7, 14):
    assert abs(daily[(N, "braked")]["x_mean"] - Y["daily"][str(N)]["x_mean"]) < 1e-6
P("   The braked rows are round 24's, reproduced.  The offset is the mean push divided by the hold.")
od = [q for q in c20["daily"] if q["N"] == 7 and q["sled"] == "not braked" and q["lims"].startswith("damping, critical")][0]
P("   The old layout (round 20): seven unbraked a day, the casing sits at %.0f m." % od["x_mean"])
OUT["daily"] = {"%d %s" % k: v for k, v in daily.items()}

# ---- the lines of a pair at the offset: round 9's elastic line solved at the casing's place
SRC9 = open(os.path.join(HERE, "round9_lines_array.py"), encoding="utf-8").read()
WANT_F = {"geff", "pot", "dia", "potv", "hang", "spring"}
WANT_A = {"mu", "RE", "W_E", "RHO", "S_OP", "S_BR", "S75", "E_CNT", "h_r", "h_a", "L", "La", "r_top", "dP_g", "TAPER", "NQ", "r_mid", "dr", "r_edge", "EXPO_mid", "EXPO_edge"}
ns = dict(np=np, math=math, R=R, root=root, brentq=brentq, P=lambda s="": None); got = set()
for node in ast.parse(SRC9).body:
    if isinstance(node, ast.FunctionDef) and node.name in WANT_F:
        exec(compile(ast.Module(body=[node], type_ignores=[]), "round9_lines_array.py", "exec"), ns); got.add(node.name)
    elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id in WANT_A for t in node.targets):
        exec(compile(ast.Module(body=[node], type_ignores=[]), "round9_lines_array.py", "exec"), ns); got.update(t.id for t in node.targets)
assert got == WANT_F | WANT_A
hang = ns["hang"]; spring = ns["spring"]
L9 = jl("round9_logs", "round9_lines_array.json")["reaction"]; B9 = L9["cases"]["B. anchorage 404 km away (700 tan 30)"]
F1 = 83.29e6/(Lr/50e3); assert abs(F1 - B9["F"]) < 1e-6 and abs(s*F1 - F_BUILT) < 1e-6
cd9 = hang(F1, B9["beta_g"], want_path=True)
places = [("one design launch, at the farthest", rD["not braked"]["x_max"])]
for N in (7, 14):
    o = daily[(N, "braked")]; places += [("%d braked a day, at the offset" % N, o["x_mean"]), ("%d braked a day, offset and swing" % N, o["x_mean"] + o["swing"])]
for N in (7, 14):
    o = daily[(N, "not braked")]; places += [("%d unbraked a day, at the offset" % N, o["x_mean"]), ("%d unbraked a day, offset and swing" % N, o["x_mean"] + o["swing"])]
xs_all = sorted(set(x for _, x in places))
tp, _ = spring(cd9, [0.0] + xs_all, say=False); tn, _ = spring(cd9, [0.0] + [-x for x in xs_all], say=False)
tab = {r[0]: r[1:] for r in tp}; tab.update({r[0]: r[1:] for r in tn})
assert all(v[0] is not None for v in tab.values())
H0, Tg0, bg0, Tt0 = tab[0.0]
assert abs(H0 - F1) < 1e-3*F1
slack = L23["both"]["slack"]
P("\n   The two lines of a pair, round 9's elastic line solved at the casing's place and scaled by %.4f (pull and stiffness go as the cross-section, round 23).  With no traffic each line pulls %.2f kN along the ring," % (s, s*H0/1e3))
P("   %.1f kN at its top and %.1f kN at the ground, %.2f deg from the vertical at the ground.  'Tighter' is the line whose anchorage the casing has moved away from." % (s*Tt0/1e3, s*Tg0/1e3, bg0))
P("   the casing's place                    | m       | tighter line: pull along the ring | of %.1f kN | top tension, of working | slacker line: pull | of %.1f kN | top tension, of working | at the ground      | the pair's net pull | all pairs and guide lines | the mean push" % (F_BUILT/1e3, F_BUILT/1e3))
pulls = {}
for nm, x in places:
    Hp, Tgp, bgp, Ttp = tab[x]; Hm, Tgm, bgm, Ttm = tab[-x]
    net = s*(Hp - Hm); tot = n_pairs*net + n_guide*k_g*x
    key = nm.split(",")[0]; mp = None
    for N in (7, 14):
        for sl in ("braked", "not braked"):
            if nm.startswith("%d %s a day" % (N, "unbraked" if sl == "not braked" else "braked")): mp = daily[(N, sl)]["mean_push_MN"]
    pulls[nm] = dict(x=x, tight_kN=s*Hp/1e3, tight_share=Hp/H0, tight_top=Ttp/Tt0, slack_kN=s*Hm/1e3, slack_share=Hm/H0, slack_top=Ttm/Tt0, slack_beta_g=bgm, slack_Tg_kN=s*Tgm/1e3, tight_beta_g=bgp, net_kN=net/1e3, all_MN=tot/1e6, mean_push_MN=mp,
                     linear_kN=k_line*x/1e3)
    P("   %-37s | %7.1f | %7.2f kN (+%.2f)                | %6.2f %%    | %6.2f %%                | %7.2f kN (%.2f) | %6.2f %%    | %6.2f %%                | %5.2f deg, %5.1f kN | %6.2f kN           | %6.3f MN                 | %s"
      % (nm, x, s*Hp/1e3, s*(Hp - H0)/1e3, 100*Hp/H0, 100*Ttp/Tt0, s*Hm/1e3, s*(Hm - H0)/1e3, 100*Hm/H0, 100*Ttm/Tt0, bgm, s*Tgm/1e3, net/1e3, tot/1e6, ("%.3f MN" % mp) if (mp and "swing" not in nm) else ""))
P("   The integrator holds the casing with the straight-line hold %.1f N/m.  The lines solved at the offset give back the mean push to within %.1f %% (last two columns), so the offsets stand."
  % (K, max(100*abs(v["all_MN"]/v["mean_push_MN"] - 1) for nm, v in pulls.items() if v["mean_push_MN"] and "swing" not in nm)))
xw = max(x for _, x in places)
P("   Slack: round 23 found that a line of this shape lies down at the ground when the casing has moved %.1f km toward its anchorage, where it still pulls %.1f %% of its design pull.  The largest place here is %.2f km." % (-slack["x_last"]/1e3, 100*slack["H"]/F1, xw/1e3))
P("   At that place the slacker line still pulls %.1f %% of its %.1f kN.  No line of a pair comes near slack." % (100*tab[-xw][0]/H0, F_BUILT/1e3))
OUT["pulls"] = pulls; OUT["no_traffic"] = dict(H_kN=s*H0/1e3, T_top_kN=s*Tt0/1e3, T_g_kN=s*Tg0/1e3, beta_g=bg0); OUT["slack"] = dict(x_km=-slack["x_last"]/1e3, share=slack["H"]/F1)

# =====================================================================================================================
P("\n" + "="*118 + "\n4. The lean and the rotation after the design launch\n" + "="*118)
h_r = R.ALTITUDE
for nm in ("not braked", "braked"):
    o = rD[nm]
    P("   sled %-10s: the casing's farthest place %.1f m.  A line 700 km high leans by atan(%.1f / 700,000) = %.4f deg.  The casing as a rigid body turns 360 x %.1f / %.0f = %.5f deg." % (nm, o["x_max"], o["x_max"], o["lean_deg"], o["x_max"], Lr, o["rot_deg"]))
x1 = rD["not braked"]["x_max"]
P("   The guide line is the 700 km line, vertical at rest.  A reaction line already leans: at the ground it is %.3f deg from the vertical at rest, %.3f deg in the tighter line and %.3f deg in the slacker line at %.0f m (round 9's elastic line)."
  % (bg0, tab[x1][2], tab[-x1][2], x1))
P("   Its pull changes by %.2f kN in the tighter line and %.2f kN in the slacker (the straight-line figure is %.2f kN)." % (s*(tab[x1][0] - H0)/1e3, s*(tab[-x1][0] - H0)/1e3, k_line*x1/1e3))
OUT["lean"] = dict(x=x1, lean_deg=rD["not braked"]["lean_deg"], rot_deg=rD["not braked"]["rot_deg"], braked_x=rD["braked"]["x_max"], braked_lean=rD["braked"]["lean_deg"], braked_rot=rD["braked"]["rot_deg"],
                   beta_g_rest=bg0, beta_g_tight=tab[x1][2], beta_g_slack=tab[-x1][2], dH_tight_kN=s*(tab[x1][0] - H0)/1e3, dH_slack_kN=s*(tab[-x1][0] - H0)/1e3)
OUT["layout"] = dict(K=K, s=s, k_line=k_line, F_built_kN=F_BUILT/1e3, n_lines=n_lines, n_guide=n_guide, k_g=k_g, period_s=PERIOD, c_crit=C, M_CAS=M_CAS, F_LIM_MN=F_LIM/1e6)

def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items() if not isinstance(v, np.ndarray)}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
P("\n(%.0f s)" % (time.time() - T_START))
with open(os.path.join(LOG_DIR, "round25_f2_rows.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
with open(os.path.join(LOG_DIR, "round25_f2_rows.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
print("log written: runs/round25_logs/round25_f2_rows.txt and .json")
