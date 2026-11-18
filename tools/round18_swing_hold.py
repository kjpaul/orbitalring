#!/usr/bin/env python3
"""
Round 18: the hold of the casing along the ring with the reaction line as round 9's hanging line (1.55 N/m), printed 5,100 N/m,
and every swing number of Section 11.6 of Volume III run again on it.
Brief: Vol_III/Claude_Code_instructions_round18_2026-11-18.md.

No new model.  The code that integrates the swing is the code of the earlier rounds, taken from their files as it stands:
  one_launch  from tools/round16_number_sheet.py  (section I: one design launch, LIMs damping, time step 10 s, 12 days)
  sim         from tools/round15_truss_array.py   (item 5: seven and fourteen launches a day, time step 20 s, 40 days, last 5 read)
Only the hold K_RING changes, from 5,800 N/m to 5,100 N/m, and with it the critical damping rate 2 sqrt(k m).
The guide line's pendulum stiffness is computed as in round 16 section H (the same lines of code) and checked against round 16's log.
The reaction line's spring is read from round 9's log (tools/round9_lines_array.py, part 3e, case B, anchorage 404 km away).

Run from Python/orbitalring:  python tools/round18_swing_hold.py
Log: runs/round18_logs/round18_swing_hold.txt and .json
"""
import ast, contextlib, io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
import numpy as np
R15 = os.path.join(HERE, "round15_truss_array.py"); R16 = os.path.join(HERE, "round16_number_sheet.py")
SRC15 = open(R15, encoding="utf-8").read(); SRC16 = open(R16, encoding="utf-8").read()
MARK = "# ---- checks against round 14"
assert MARK in SRC15
ns = {"__name__": "round15_prefix", "__file__": R15}
_argv = sys.argv; sys.argv = [R15]
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(SRC15.split(MARK)[0], R15, "exec"), ns)
sys.argv = _argv
def take(src, path, name, into):
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            exec(compile(ast.Module(body=[node], type_ignores=[]), path, "exec"), into); return
    raise KeyError(name)
pot = ns["pot"]; geff = ns["geff"]; RHO = ns["RHO"]; S_OP = ns["S_OP"]; BOOKENV = ns["BOOKENV"]; Lr = BOOKENV.Lr
LOG9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))
LOG15 = json.load(open(os.path.join(ROOT, "runs", "round15_logs", "round15_truss_array.json"), encoding="utf-8"))
LOG16 = json.load(open(os.path.join(ROOT, "runs", "round16_logs", "round16_number_sheet.json"), encoding="utf-8"))

LOG_DIR = os.path.join(ROOT, "runs", "round18_logs"); os.makedirs(LOG_DIR, exist_ok=True)
buf = io.StringIO()
def P(s=""):
    print(s); buf.write(s + "\n"); sys.stdout.flush()
OUT = {}

# ---------------------------------------------------------------------------------------------------------------------
P("Round 18, item 1.  The hold along the ring.")
# the guide line as a pendulum: the lines of round 16 section H
A_tip = 100e3*geff(0)/S_OP; nq = 20000; zz = (np.arange(nq) + 0.5)*700e3/nq; ex = np.exp(RHO*np.array([pot(0, z) for z in zz])/S_OP)
k_p = 1.0/float((1/(S_OP*A_tip*ex)).sum()*700e3/nq)
assert abs(k_p - LOG16["I"]["k_p"]) < 1e-9, (k_p, LOG16["I"]["k_p"])
n_st = Lr/50e3
SP = LOG9["reaction"]["spring"]["B. anchorage 404 km away"]; k_r = SP["k_N_per_m"]; k_r_book = LOG16["I"]["k_r_book"]
assert abs(k_r - LOG16["I"]["k_r_hang"]) < 1e-12
k_hang = n_st*(k_r + 2*k_p); K_PRINT = 5100.0
P("  Anchor points: %.2f (the ring's %.1f km over 50 km), printed 889." % (n_st, Lr/1e3))
P("  Reaction line, hanging, anchorage 404 km away (round 9, part 3e, case B): %.4f N for each meter the casing moves." % k_r)
P("  Guide line, a vertical 100 t line as a pendulum: %.4f N/m ; two of them %.4f N/m." % (k_p, 2*k_p))
P("  Ring: %.2f x (%.4f + 2 x %.4f) = %.1f N/m, printed 5,100.  The guide lines give %.1f %% of it." % (n_st, k_r, k_p, k_hang, 100*2*k_p/(k_r + 2*k_p)))
P("  (The chapter script's straight line at 30 deg: %.4f N/m, ring %.1f N/m, printed 5,800.)" % (k_r_book, n_st*(k_r_book + 2*k_p)))
P("  How round 9's figure depends on the pull.  It is the slope of the force along the ring at the design place, taken between 100 m either side,")
P("  for a line built to carry the design force of one launch (%.2f kN along the ring) at that place.  The line gets softer as it slackens:" % SP["table"]["0.0"]["H_kN"])
P("    casing moved toward the anchorage by | force along the ring | of the design force | mean slope from the design place")
H0 = SP["table"]["0.0"]["H_kN"]; dep = []
for key in sorted(SP["table"], key=float):
    x = float(key); H = SP["table"][key]["H_kN"]
    if x == 0: continue
    sec = (H - H0)*1e3/x; dep.append(dict(x_m=x, H_kN=H, share=H/H0, secant_N_m=sec))
    P("    %+9.1f km                           | %6.2f kN            | %5.1f %%             | %.3f N/m" % (-x/1e3, H, 100*H/H0, sec))
P("  (A minus sign in the first column is the casing moved away from the anchorage, which tightens the line.)")
xs_ = sorted(float(k_) for k_ in SP["table"]); loc = []
for xa, xb in zip(xs_[:-1], xs_[1:]):
    ka = "%.1f" % xa; kb = "%.1f" % xb; loc.append(dict(x_from=xa, x_to=xb, slope=(SP["table"][kb]["H_kN"] - SP["table"][ka]["H_kN"])*1e3/(xb - xa)))
P("  Slope between neighboring rows (the local spring): " + " ; ".join("%.0f to %.0f km: %.2f" % (-q["x_to"]/1e3, -q["x_from"]/1e3, q["slope"]) for q in loc[:4]) + " N/m.")
k_guides = n_st*2*k_p; T_g = 2*math.pi*math.sqrt(12000.0*Lr/k_guides)
P("  If the reaction lines gave nothing, the guide lines alone would hold the ring with %.0f N/m (period %.1f h).  That is the floor." % (k_guides, T_g/3600))
OUT["hold"] = dict(n_st=n_st, k_r_hang=k_r, k_p=k_p, k_ring=k_hang, k_print=K_PRINT, k_r_book=k_r_book, k_ring_book=n_st*(k_r_book + 2*k_p), guide_share=2*k_p/(k_r + 2*k_p), H0_kN=H0, dependence=dep, local=loc, k_guides_only=k_guides, period_guides_only_h=T_g/3600)

# ---------------------------------------------------------------------------------------------------------------------
P("\nRound 18, item 2.  The swing rows on the new hold.")
M_CAS = 12000.0*Lr; N_SITES = 88946; F_LIM = 2*401.2*N_SITES; V_REL = 8846.96; F_LINE = 83.29e6/n_st
def run(K, full=True):
    """All the rows at one hold K.  Each integrator gets its own namespace with the constants of the round it comes from."""
    g16 = dict(np=np, math=math, K_RING=K, M_CAS=M_CAS, F_LIM=F_LIM); take(SRC16, R16, "one_launch", g16)
    J16 = (10.25e6*30000.0, 5.0e6*30000.0)                                     # round 16 section I
    om = math.sqrt(K/M_CAS); C = 2*math.sqrt(K*M_CAS); T = 2*math.pi/om
    o = dict(K=K, period_s=T, period_h=T/3600, c_crit=C, swing_m=J16[0]/M_CAS/om, swing_braked_m=J16[1]/M_CAS/om, v_full=J16[0]/M_CAS, v_net=J16[1]/M_CAS, one={}, rows=[])
    for nm, J in (("not braked", J16[0]), ("braked", J16[1])):
        xm, fm, tdie, tpk = g16["one_launch"](J, C)
        o["one"][nm] = dict(x_max=xm, F_lim_MN=fm/1e6, die_h=tdie, peak_h=tpk, line_kN=k_r*xm/1e3, GW=fm*V_REL/1e9, lean_deg=math.degrees(math.atan(xm/700e3)), rot_deg=360*xm/Lr)
    if not full: return o
    J15 = (0.281*M_CAS, 0.576*M_CAS)                                           # round 15 item 5
    g15 = dict(np=np, math=math, K_RING=K, M_CAS=M_CAS, F_LIM=F_LIM, QUICK=False); take(SRC15, R15, "sim", g15)
    for Nn in (7, 14):
        for nm, J in (("braked", J15[0]), ("not braked", J15[1])):
            Fm = Nn*J/86400.0
            for cn, c, hold in (("none", 0.0, False), ("damping, critical", C, False), ("damping and holding the mean push", C, True)):
                if c: r_ = g15["sim"](Nn, J, c, hold)
                else:
                    xm = Fm/K; r_ = dict(x_max=2*xm + J/(M_CAS*om), x_mean=xm, swing=xm + J/(M_CAS*om), F_lim_peak_MN=0.0, F_lim_mean_MN=0.0)      # the bound of round 15, from rest
                o["rows"].append(dict(N=Nn, sled=nm, mean_push_MN=Fm/1e6, lims=cn, **r_, line_kN=k_r*r_["x_max"]/1e3, lean_deg=math.degrees(math.atan(r_["x_max"]/700e3))))
    return o
def show(o):
    P("  Hold %.1f N/m ; casing %.4e kg ; period %.0f s (%.2f h) ; critical damping %.3e N s/m" % (o["K"], M_CAS, o["period_s"], o["period_h"], o["c_crit"]))
    P("  One design launch leaves the casing moving at %.4f m/s, %.4f m/s with the sled braked.  No damping: swing %.0f m, %.0f m braked." % (o["v_full"], o["v_net"], o["swing_m"], o["swing_braked_m"]))
    for nm in ("not braked", "braked"):
        q = o["one"][nm]
        P("    LIMs damping, sled %-10s: farthest %.0f m, %.2f h after the start of the launch ; within 1 %% of that after %.1f h ; LIM push at most %.2f MN of %.1f MN, %.0f GW ;"
          % (nm, q["x_max"], q["peak_h"], q["die_h"], q["F_lim_MN"], F_LIM/1e6, q["GW"]))
        P("    %s reaction line %.2f kN (%.1f %% of %.1f kN) ; lines lean %.3f deg ; casing turned %.4f deg" % (" "*26, q["line_kN"], 100*q["line_kN"]*1e3/F_LINE, F_LINE/1e3, q["lean_deg"], q["rot_deg"]))
    if o["rows"]:
        P("    launches | sled       | mean push | LIMs                              | casing sits at | swings about it | farthest | reaction line there | lean     | LIM push: peak, mean")
        for r_ in o["rows"]:
            P("    %2d a day | %-10s | %5.2f MN  | %-33s | %7.0f m      | %7.0f m       | %6.0f m | %5.2f kN (%4.1f %%)    | %.3f deg | %5.1f, %5.1f MN%s"
              % (r_["N"], r_["sled"], r_["mean_push_MN"], r_["lims"], r_["x_mean"], r_["swing"], r_["x_max"], r_["line_kN"], 100*r_["line_kN"]*1e3/F_LINE, r_["lean_deg"], r_["F_lim_peak_MN"], r_["F_lim_mean_MN"],
                 "  (bound, from rest)" if r_["lims"] == "none" else ""))
P("  The launch is a steady push over 3,690 s [ASSUMED shape, as rounds 15 and 16].  LIMs: %.1f MN in all (round 15)." % (F_LIM/1e6))
new = run(K_PRINT); show(new)
P("\n  For comparison, the same code at the 5,800 N/m the book has now:")
old = run(5800.0); show(old)
P("\n  And the single launch at the exact %.1f N/m, to see that printing 5,100 changes no printed digit:" % k_hang)
exact = run(k_hang, full=False); show(exact)
OUT["new"] = new; OUT["old_5800"] = old; OUT["exact"] = exact
# checks against the earlier logs at 5,800
o16 = LOG16["I"]["one"]
for nm in ("not braked", "braked"):
    assert abs(old["one"][nm]["x_max"] - o16[nm]["x_max"]) < 1e-6 and abs(old["one"][nm]["F_lim_MN"] - o16[nm]["F_lim_MN"]) < 1e-9, nm
for a, b in zip(old["rows"], LOG15["item5"]["rows"]):
    assert a["N"] == b["N"] and a["sled"] == b["sled"] and abs(a["x_mean"] - b["x_mean"]) < 1e-6 and abs(a["swing"] - b["swing"]) < 1e-6 and abs(a["F_lim_peak_MN"] - b["F_lim_peak_MN"]) < 1e-9, (a, b)
P("\n  Check: at 5,800 N/m this script gives round 16's single-launch numbers and round 15's twelve traffic rows again, to the last digit logged.")
dv = LOG15["item5"]["dv_cable_per_launch"]; dv2 = 0.281*M_CAS/(103569.71*Lr)
assert abs(dv - dv2) < 1e-12
P("  The cable's gain of speed if the LIMs hold the mean push does not contain the hold: %.5f m/s for each braked launch, %.3f m/s a day and %.1f m/s a year at seven a day.  Unchanged."
  % (dv, 7*dv, 7*dv*365.25))
OUT["dv_cable"] = dict(per_launch=dv, per_year_7=7*dv*365.25)
OUT["const"] = dict(M_CAS=M_CAS, F_LIM_MN=F_LIM/1e6, F_LINE_kN=F_LINE/1e3, Lr=Lr)
with open(os.path.join(LOG_DIR, "round18_swing_hold.txt"), "w", encoding="utf-8") as fh: fh.write(buf.getvalue())
with open(os.path.join(LOG_DIR, "round18_swing_hold.json"), "w", encoding="utf-8") as fh: json.dump(OUT, fh, indent=1, default=float)
print("ok")
