#!/usr/bin/env python3
"""
Round 16: the number sheet for Chapters 11 and 12 of Volume III on the ruled truss frame, and Paul's deployment sequence of the array.
Brief: Vol_III/Claude_Code_instructions_round16_2026-11-18.md (as it stood at 11:08 on 18 November 2026).

The ruling: no ground spreader arms ; a truss frame under the sheet every 10 km ; sheet 1.11 km, 0.5 kg/m2 ; two suspension lines on each
leg (wing) at 2/5 and 4/5 of the leg from the hinge, of fixed length ; the center array support line to the hinge, on a spool ; lines 100 km.

The frame model, the wall rule (NASA SP-8007, Equation 44) and the budget are round 15's.  The first part of tools/round15_truss_array.py
(up to its checks against round 14) is run here as it stands, with its printing muted, so that the functions are the same code.
New here: the lines as springs (a CNT line stretches 3.75 percent at 18.75 GPa, which is 3.75 km in 100 km), the opening with two lines,
the drawing out of the sheet with the ridge let down, the tables at a finer spacing, the ring's own lines, the swing, the conductors,
and the deployment sequence day by day.

Run from Python/orbitalring:  python tools/round16_number_sheet.py
Log: runs/round16_logs/round16_number_sheet.txt and .json ; the number sheet: runs/round16_logs/round16_numbers.md
"""
import ast, contextlib, csv, io, json, math, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq

T0 = time.time()
TRAPZ = getattr(np, "trapezoid", None) or np.trapz
R15 = os.path.join(HERE, "round15_truss_array.py")
SRC15 = open(R15, encoding="utf-8").read()
MARK = "# ---- checks against round 14"
assert MARK in SRC15
ns = {"__name__": "round15_prefix", "__file__": R15}
_argv = sys.argv; sys.argv = [R15]
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(SRC15.split(MARK)[0], R15, "exec"), ns)
    for node in ast.parse(SRC15).body:
        if isinstance(node, ast.FunctionDef) and node.name == "sheet_at":
            exec(compile(ast.Module(body=[node], type_ignores=[]), R15, "exec"), ns)
sys.argv = _argv
g = ns.get
frame_system = ns["frame_system"]; truss = ns["truss"]; layout = ns["layout"]; lp_bending = ns["lp_bending"]; Env = ns["Env"]; BOOKENV = ns["BOOKENV"]
catenary = ns["catenary"]; dia = ns["dia"]; geff = ns["geff"]; pot = ns["pot"]; sheet_at = ns["sheet_at"]; lit = ns["lit"]
MU = ns["MU"]; RE = ns["RE"]; W_E = ns["W_E"]; RHO = ns["RHO"]; S_OP = ns["S_OP"]; S75 = ns["S75"]; E_CNT = ns["E_CNT"]; RHO_CF = ns["RHO_CF"]
HS45 = ns["HS45"]; DR45 = ns["DR45"]; d5 = ns["d5"]; P9 = ns["P9"]; LOG9 = ns["LOG9"]; LOG14 = ns["LOG14"]; TW_PER_KM_600 = ns["TW_PER_KM_600"]
FOS = ns["FOS"]; BRACE = ns["BRACE"]; E_CF = ns["E_CF"]; SC_CF = ns["SC_CF"]; K_SRC = ns["K_SRC"]; WALL_NAME = ns["WALL_NAME"]; JOINTS = ns["JOINTS"]
LOG15 = json.load(open(os.path.join(ROOT, "runs", "round15_logs", "round15_truss_array.json"), encoding="utf-8"))
LOG13T = open(os.path.join(ROOT, "runs", "round13_logs", "round13_spreader_arms.txt"), encoding="utf-8").read()

LOG_DIR = os.path.join(ROOT, "runs", "round16_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
def head(s): P("\n" + "="*118 + "\n" + s + "\n" + "="*118)
OUT = {}; SHEET = []
SCRIPT = "tools/round16_number_sheet.py"
def N(sec, what, show, exact, goes, src=None, key=None):
    """One row of the number sheet: section letter, what it is, the value as the book should print it, the exact value, where it goes."""
    SHEET.append(dict(sec=sec, what=what, show=show, exact=exact, src=src or ("%s, log section %s" % (SCRIPT, key or sec.upper())), goes=goes))
def fx(v, n=4):
    if isinstance(v, str): return v
    return ("%%.%dg" % (n + 2)) % v

ga = BOOKENV.g_a; gr = BOOKENV.g_r; kc = BOOKENV.kc; D = BOOKENV.D; Lr = BOOKENV.Lr; La = BOOKENV.La
W5 = d5["W"]; S10 = 10e3; s_half = W5/2; qs = 0.5*ga; b_half = HS45*s_half; drop = DR45*s_half
TOPF = 1/(1 - kc/(1 + kc))                         # pull at the top of a line of one thickness over the pull at its lower end
EPS75 = S75/E_CNT
P("Round 16: the number sheet on the ruled truss frame, and Paul's deployment sequence.")
P("ring 700 km, array 600 km ; g at the array %.4f, at the ring %.4f m/s2 ; lines %.2f GPa, modulus %.0f GPa, so a line at its working stress is %.2f %% longer than slack: %.2f km in 100 km"
  % (ga, gr, S75/1e9, E_CNT/1e9, 100*EPS75, EPS75*100))
P("sheet %.1f m ; half footprint %.1f m ; edge %.1f m below the ridge ; line fraction kc %.4f ; top pull / lower pull %.4f" % (W5, b_half, drop, kc, TOPF))

# =====================================================================================================================
head("A.  The frame of the design case: 0.5 kg/m2, sheet 1.11 km, a frame every 10 km, two lines on each leg at 2/5 and 4/5")
# =====================================================================================================================
lay = layout(2, 0.8); layf = layout(2, 0.8, True)
fr = frame_system(0.5, W5, S10, n=2, f=0.8); frf = frame_system(0.5, W5, S10, n=2, f=0.8, hinge_free=True)
t = fr["truss"]; tf = frf["truss"]; leg = fr["leg"]; Lf = fr["Lf"]; F = fr["F"]; Mfr = t["mass"]; Wl = Mfr*ga/2; lever = 0.866*t["d"]
P("  Line pulls set for the least bending from the frame's own weight (a linear program, round 15):")
P("    hinge carrying moment:  cM %.5f, cV %.3f ; shares of a leg's weight: inner %.3f, outer %.3f, center line %.3f for each leg ; moment at the hinge %+.4f w l^2"
  % (lay["cM"], lay["cV"], lay["R"][0], lay["R"][1], lay["Rc"], lay["M0"]))
P("    no moment at the hinge: cM %.5f, cV %.3f ; shares: inner %.3f, outer %.3f, center line %.3f for each leg" % (layf["cM"], layf["cV"], layf["R"][0], layf["R"][1], layf["Rc"]))
def show_t(nm, x):
    tt = x["truss"]
    P("    %-28s %7.1f t ; %.1f m long, %.2f m deep, %d bays of %.2f m ; chords %.0f mm across, wall %.2f mm, %.0f MPa at work (%s) ; battens %.0f mm, wall %.2f mm ; diagonals %.1f mm"
      % (nm, tt["mass"]/1e3, x["Lf"], tt["d"], tt["bays"], tt["lb"], 2e3*tt["r_chord"], 1e3*tt["t_chord"], tt["stress_chord_MPa"], "whole strut governs" if tt["euler_governs"] else "chord in its bay governs",
         2e3*tt["r_batten"], 1e3*tt["t_batten"], tt["d_diag_mm"]))
    P("    %-28s chords %.1f t + battens %.1f t + joints and hinge %.1f t + diagonals %.2f t ; chord force at work %.2f MN ; bending from its own weight %.1f MN m"
      % ("", tt["m_chords"]/1e3, tt["m_battens"]/1e3, tt["m_joints"]/1e3, tt["m_diag"]/1e3, tt["N_chord_MN"], tt["M_weight_MNm"]))
show_t("hinge carrying moment", fr); show_t("no moment at the hinge", frf)
M0_MNm = lay["M0"]*Wl*leg/1e6
P("  The sign of the hinge moment matters.  The stop (pins on the lower chords, upper chords butting) takes a moment that tries to lift the tips.  The least-bending")
P("  layout asks for %.1f MN m the other way (the tips heavy), which is %.2f MN of pull across the upper chords.  A butt joint cannot pull.  But once the sheet is out the"
  % (-M0_MNm, -M0_MNm/lever))
P("  push of %.2f MN runs through the joint, and the line of thrust only moves %.2f m toward the pins (the joint stays shut as long as that is under %.1f m)."
  % (F/1e6, -M0_MNm*1e6/F, 0.2887*t["d"]))
P("  So the lighter frame is right for the working state, and the frame with no moment at the hinge is the one that needs no push to hold its joint.")
P("  Laminate: %s ; axial modulus %.1f GPa, compressive strength %.0f MPa, wall constant K %.1f GPa (SP-8007 Eq. 44), density %.0f kg/m3." % (WALL_NAME, E_CF/1e9, SC_CF/1e6, K_SRC/1e9, RHO_CF))
P("  Push %.3f MN ; weight of one leg %.3f MN (%.1f t) ; chord stress at work %.0f MPa against %.0f MPa (%.1f %%) ; bays %.0f chord radii long."
  % (F/1e6, Wl/1e6, Mfr/2e3, t["stress_chord_MPa"], SC_CF/1e6, 100*t["stress_chord_MPa"]*1e6/SC_CF, t["lb"]/t["r_chord"]))
r15_1 = LOG15["item1"]["0.800"]["frame_t"]
P("  Round 15 for comparison: one line at 4/5, %.0f t ; one line at 2/3, %.0f t ; round 14 (isotropic wall rule) 452 t." % (r15_1, LOG15["item1"]["0.667"]["frame_t"]))
OUT["A"] = dict(hinge_moment=dict(mass_t=Mfr/1e3, **{k: (float(v) if not isinstance(v, bool) else v) for k, v in t.items()}), hinge_free=dict(mass_t=tf["mass"]/1e3, **{k: (float(v) if not isinstance(v, bool) else v) for k, v in tf.items()}),
                Lf=Lf, leg=leg, push_MN=F/1e6, Wl_MN=Wl/1e6, cM=lay["cM"], cV=lay["cV"], R=lay["R"], Rc=lay["Rc"], M0=lay["M0"], cM_free=layf["cM"], R_free=layf["R"], Rc_free=layf["Rc"],
                M0_MNm=M0_MNm, tie_MN=-M0_MNm/lever, thrust_shift_m=-M0_MNm*1e6/F, kern_m=0.2887*t["d"], lever_m=lever)
G11 = "Ch11 line 24, Table 11.2 ; Ch12 Section 12.3 (lines 52, 54)"
N("a", "Frame mass, hinge carrying moment (the working state)", "%.0f t" % (Mfr/1e3), Mfr/1e3, G11)
N("a", "Frame mass, no moment at the hinge", "%.0f t" % (tf["mass"]/1e3), tf["mass"]/1e3, G11)
N("a", "Frame length, tip to tip", "%.0f m" % Lf, Lf, G11); N("a", "Length of one leg (wing)", "%.0f m" % leg, leg, G11)
N("a", "Frame depth (side of the triangle)", "%.0f m" % t["d"], t["d"], G11)
N("a", "Chords: three tubes, across and wall", "%.2f m, %.1f mm" % (2*t["r_chord"], 1e3*t["t_chord"]), "%.4f m, %.3f mm" % (2*t["r_chord"], 1e3*t["t_chord"]), G11)
N("a", "Battens: across and wall", "%.2f m, %.1f mm" % (2*t["r_batten"], 1e3*t["t_batten"]), "%.4f m, %.3f mm" % (2*t["r_batten"], 1e3*t["t_batten"]), G11)
N("a", "Bays: number and length", "%d of %.1f m" % (t["bays"], t["lb"]), "%d, %.3f m" % (t["bays"], t["lb"]), G11)
N("a", "Diagonals (CNT, crossed, in each face)", "%.0f mm" % t["d_diag_mm"], t["d_diag_mm"], G11)
N("a", "Mass split: chords, battens, joints and hinge, diagonals", "%.0f, %.0f, %.0f, %.0f t" % (t["m_chords"]/1e3, t["m_battens"]/1e3, t["m_joints"]/1e3, t["m_diag"]/1e3),
  "%.2f, %.2f, %.2f, %.2f t" % (t["m_chords"]/1e3, t["m_battens"]/1e3, t["m_joints"]/1e3, t["m_diag"]/1e3), G11)
N("a", "Laminate", "IM7 / 8552, 75 % of plies at 0 deg, 25 % at +-45 deg", "E %.1f GPa, strength %.0f MPa, K %.2f GPa, %.0f kg/m3" % (E_CF/1e9, SC_CF/1e6, K_SRC/1e9, RHO_CF), G11, src="tools/round15_truss_array.py, log ITEM 7 (NASA SP-8007 Eq. 44 ; Hexcel HexPly 8552 data sheet)")
N("a", "Push on a frame", "%.1f MN" % (F/1e6), F/1e6, "Ch12 line 52 (was the 23 MN into each arm) ; Table 11.1")
N("a", "Weight of one leg", "%.1f MN (%.0f t)" % (Wl/1e6, Mfr/2e3), Wl/1e6, G11)
N("a", "Chord stress at work, against the laminate's strength", "%.0f MPa of %.0f MPa" % (t["stress_chord_MPa"], SC_CF/1e6), "%.2f MPa, %.1f MPa" % (t["stress_chord_MPa"], SC_CF/1e6), G11)
N("a", "Bending from the frame's own weight, share of w l^2 (hinge carrying moment ; none)", "%.3f ; %.3f" % (lay["cM"], layf["cM"]), "%.5f ; %.5f" % (lay["cM"], layf["cM"]), G11)

# =====================================================================================================================
head("B.  The lines of one frame, as springs: the opening, the open frame, the drawing out of the sheet")
# =====================================================================================================================
xo = 0.8*leg; xi = 0.4*leg; Zt = D + drop                      # the open frame is at the height of the edges, Zt below the top of the lines
P("  The suspension lines are of fixed length.  The center array support line is on a spool.  A CNT line is a soft spring: E A / L.  At its working stress it is")
P("  %.2f km longer than slack.  So every pull below is found with the lines as springs and the leg as a stiff bar hinged at the center:" % (EPS75*D/1e3))
P("  for each depth of the hinge below the top of the lines, the leg's angle is the one of least energy, the spool's pull is what is left of the weight.")
def state(h, A, L0, phi_lo=math.radians(1.0)):
    """Hinge h below the top of the lines.  A, L0: areas and slack lengths of the (inner, outer) suspension line.  The leg's angle phi from the vertical
    (0 folded, 90 deg flat) minimises weight + spring energy on [phi_lo, 90 deg].  Returns phi, the two line pulls, the center line's pull (both legs), the hinge moment at flat."""
    ph = np.linspace(phi_lo, math.pi/2, 6001); U = -Wl*(leg/2)*np.cos(ph); Ts = []
    for x, a, l0 in zip((xi, xo), A, L0):
        dist = np.hypot(x*np.sin(ph), h + x*np.cos(ph)); ext = np.maximum(dist - l0, 0.0); k = E_CNT*a/l0
        U = U + 0.5*k*ext*ext; Ts.append((k*ext, (h + x*np.cos(ph))/dist))
    i = int(np.argmin(U))
    T = [float(v[0][i]) for v in Ts]; cz = [float(v[1][i]) for v in Ts]
    Tc = 2*Wl - 2*sum(a*b for a, b in zip(T, cz))
    Mh = (T[0]*cz[0]*xi + T[1]*cz[1]*xo - Wl*leg/2) if i == len(ph) - 1 else 0.0
    return dict(h=h, phi=math.degrees(ph[i]), Ti=T[0], To=T[1], Tc=Tc, Mh=Mh)
def design(RiT, RoT, ri_size, ro_size):
    """Slack lengths cut so that, flat at the depth Zt, the lines carry RiT and RoT of a leg's weight ; areas for ri_size, ro_size of a leg's weight at 18.75 GPa at the top."""
    A = (ri_size*Wl*TOPF/S75, ro_size*Wl*TOPF/S75)
    L0 = (math.hypot(xi, Zt)/(1 + RiT*Wl/(E_CNT*A[0])), math.hypot(xo, Zt)/(1 + RoT*Wl/(E_CNT*A[1])))
    return A, L0
def sequence(A, L0):
    hs = np.concatenate((np.arange(Zt - 9000.0, Zt - 600.0, 50.0), np.arange(Zt - 600.0, Zt + 1500.0, 2.0)))
    return [state(h, A, L0) for h in hs]
ri = ro = 0.5
for it in range(6):
    A, L0 = design(lay["R"][0], lay["R"][1], ri, ro)
    seq = sequence(A, L0)
    flat = next(s for s in seq if s["phi"] > 89.999)                       # the first depth at which the frame is flat with no push on it
    upto = [s for s in seq if s["h"] <= flat["h"] + 1e-9]
    ri_n = max(max(s["Ti"] for s in upto)/Wl, lay["R"][0]); ro_n = max(max(s["To"] for s in upto)/Wl, lay["R"][1])
    if abs(ri_n - ri) + abs(ro_n - ro) < 1e-4: break
    ri, ro = ri_n, ro_n
taut = next(s for s in seq if s["Ti"] + s["To"] > 1.0)                      # first pull in a suspension line
start = next(s for s in seq if s["phi"] > 1.0 + 0.05)                       # the legs begin to turn
P("  Suspension lines sized for the largest pull they see up to the flat frame, at 18.75 GPa at the top: inner %.3f, outer %.3f of a leg's weight (%.2f, %.2f MN)." % (ri, ro, ri*Wl/1e6, ro*Wl/1e6))
P("  Slack lengths: inner %.3f km, outer %.3f km (cut so that the open frame at work, %.3f km below the top of the lines, has the least-bending pulls)." % (L0[0]/1e3, L0[1]/1e3, Zt/1e3))
P("\n    hinge below its working place | leg from the vertical | inner line | outer line | center line (both legs) | of the frame's weight on the spool")
shown = set()
for s in seq:
    dh = s["h"] - Zt
    key = round(dh/250)*250 if dh < -700 else (round(dh/50)*50 if abs(dh) > 60 else round(dh/10)*10)
    if key in shown or abs(dh - key) > 1.01: continue
    if dh < taut["h"] - Zt - 300 or dh > flat["h"] - Zt + 120: continue
    shown.add(key)
    P("    %+8.0f m                    | %6.2f deg            | %6.3f MN  | %6.3f MN  | %6.3f MN               | %5.1f %%%s"
      % (dh, s["phi"], s["Ti"]/1e6, s["To"]/1e6, s["Tc"]/1e6, 50*s["Tc"]/Wl, "   stop shut, %.0f MN m" % (s["Mh"]/1e6) if s["Mh"] > 0 else ""))
jump = max(abs(b["phi"] - a["phi"]) for a, b in zip(seq[:-1], seq[1:]) if b["h"] - a["h"] < 3)
E_brake = sum(0.5*(a["Tc"] + b["Tc"])*(b["h"] - a["h"]) for a, b in zip(seq[:-1], seq[1:]) if a["h"] >= taut["h"] and b["h"] <= flat["h"])
tc_min = min(s["Tc"] for s in upto)
work = state(Zt, A, L0)                                                         # the working state needs the push to hold the joint: computed with the stop tied
Ti_w, To_w = lay["R"][0]*Wl, lay["R"][1]*Wl; Tc_w = 2*lay["Rc"]*Wl
P("\n  The suspension lines first pull with the hinge %.0f m above its working place.  The legs begin to turn at %.0f m above it and are flat at %+.0f m." % (Zt - taut["h"], Zt - start["h"], flat["h"] - Zt))
P("  So the take-up and the opening use %.2f km of the spool's pay-out, not the %.0f m of the stiff-line geometry (0.8 of the leg)." % ((flat["h"] - taut["h"])/1e3, xo))
P("  The largest change of the leg's angle in a 2 m step of the spool is %.2f deg: it opens gradually, with no snap.  The center line's least pull on the way is %.3f MN (it stays taut)." % (jump, tc_min/1e6))
P("  The spool's brake takes %.2f GJ from first pull to flat (round 15, stiff line at 4/5: 0.71 GJ)." % (E_brake/1e9))
P("  Flat with no push, the joint is just shut.  With the sheet out, the push holds the joint and the spool takes up %.0f m to bring the frame to its working pulls:" % (flat["h"] - Zt))
P("    working pulls: inner %.3f MN, outer %.3f MN, center line %.3f MN.   At flat with no push: inner %.3f, outer %.3f, center %.3f MN." % (Ti_w/1e6, To_w/1e6, Tc_w/1e6, flat["Ti"]/1e6, flat["To"]/1e6, flat["Tc"]/1e6))
# center line: phase 2 carries both legs
Tc_size = 2*Wl
A_c = Tc_size*TOPF/S75
def line(nm, T_size, T_open, T_opening, lean_x):
    Aa = T_size*TOPF/S75; m = RHO*Aa*D
    return dict(name=nm, size_MN=T_size/1e6, open_MN=T_open/1e6, opening_MN=T_opening/1e6, top_MN=T_size*TOPF/1e6, d_mm=dia(T_size*TOPF, S75), mass_t=m/1e3, area_cm2=Aa*1e4,
                stress_open_GPa=T_open*TOPF/Aa/1e9, stretch_open_km=T_open/(E_CNT*Aa)*D/1e3, k_N_m=E_CNT*Aa/D, lean_deg=math.degrees(math.atan(lean_x/Zt)))
LN = [line("center array support line", Tc_size, Tc_w, max(s["Tc"] for s in upto if s["h"] >= taut["h"]), 0.0),
      line("inner suspension line (2/5)", ri*Wl, Ti_w, max(s["Ti"] for s in upto), xi),
      line("outer suspension line (4/5)", ro*Wl, To_w, max(s["To"] for s in upto), xo)]
P("\n    line                         | open, at work | largest while opening | sized for          | across  | mass    | stress at work | longer than slack at work | spring     | lean")
for l_ in LN:
    P("    %-28s | %6.3f MN     | %6.3f MN             | %6.3f MN%s | %5.1f mm | %5.1f t | %5.2f GPa      | %5.2f km                   | %5.0f N/m  | %.3f deg"
      % (l_["name"], l_["open_MN"], l_["opening_MN"], l_["size_MN"], " (both legs, folded)" if l_ is LN[0] else "                    ", l_["d_mm"], l_["mass_t"], l_["stress_open_GPa"], l_["stretch_open_km"], l_["k_N_m"], l_["lean_deg"]))
m_lines = LN[0]["mass_t"] + 2*LN[1]["mass_t"] + 2*LN[2]["mass_t"]
C_LINES = m_lines*1e3/(kc*Mfr)
P("  The five lines: %.1f t, %.2f kg for each meter of ring at 10 km, %.1f %% of the frame's mass (round 15 counted %.1f %%: the open frame only)." % (m_lines, m_lines*1e3/S10, 100*m_lines*1e3/Mfr, 100*kc))
P("  Working stress: 18.75 GPa at the top, as the array support lines of Section 12.4, because the pull is steady once the frame is open and the larger pulls come once, slowly.")
P("  At the 12.5 GPa of the anchor lines every line here would be 1.5 times heavier: %.1f t." % (1.5*m_lines*(1 + 0.0)))
alt_settle = None
# ---- the stop ----
M_stop_flat = flat["Mh"]; M_slack = (ro*0 + 0)  # placeholder, replaced below
sl = state(flat["h"] + 400.0, A, L0)
P("\n  The stop.  Flat with no push, any more pay-out of the spool puts more of the weight on the suspension lines and shuts the joint harder:")
for dh in (0.0, 100.0, 200.0, 400.0, 800.0):
    s = state(flat["h"] + dh, A, L0)
    P("    spool paid out %4.0f m past flat: center line %6.3f MN ; joint held shut with %6.1f MN m (%.2f MN across the upper chords) ; lines %.3f and %.3f MN (%.0f and %.0f %% of what they are sized for)"
      % (dh, s["Tc"]/1e6, s["Mh"]/1e6, s["Mh"]/lever/1e6, s["Ti"]/1e6, s["To"]/1e6, 100*s["Ti"]/(ri*Wl), 100*s["To"]/(ro*Wl)))
P("  At work (sheet out, spool taken up): the least-bending pulls leave %.1f MN m that would open the joint ; the push of %.2f MN holds it (thrust line %.2f m off center, limit %.1f m)."
  % (-M0_MNm, F/1e6, -M0_MNm*1e6/F, 0.2887*t["d"]))
# ---- drawing the sheet out ----
P("\n  Drawing the sheet out, first with the ridge held where it is (round 15's way), two lines on each leg, the center line slack, no tie at the joint:")
P("    tether  | push     | weight on each tip | inner line | outer line | joint          | chord force / sized | the outer line would stretch by")
rows = []
for tl in (585, 570, 540, 500, 450, 400, 300, 200, 100, 50, 20, 5, 0.5):
    def g_(al):
        sh = sheet_at(b_half - tl*math.cos(al), -drop - tl*math.sin(al), s_half, qs)
        return None if sh is None else (-sh[1]) - sh[0]*math.tan(al)
    al_s = np.linspace(1e-4, math.pi/2 - 1e-4, 2000); vals = [g_(a_) for a_ in al_s]; al0 = None
    for i in range(len(al_s) - 1):
        if vals[i] is not None and vals[i + 1] is not None and vals[i]*vals[i + 1] <= 0: al0 = brentq(g_, al_s[i], al_s[i + 1]); break
    if al0 is None: continue
    sh = sheet_at(b_half - tl*math.cos(al0), -drop - tl*math.sin(al0), s_half, qs); Fp = sh[0]*S10; Pt = -sh[1]*S10
    Tin = 0.75*Wl - 0.5*Pt; Tout = 0.25*Wl + 1.5*Pt; Mh = 0.0
    if Tin < 0: Tin = 0.0; Tout = Wl + Pt; Mh = Tout*xo - Wl*leg/2 - Pt*leg
    xs = np.linspace(0, leg, 2001); Mx = -Pt*(leg - xs) - (Wl/leg)*(leg - xs)**2/2 + Tin*np.maximum(xi - xs, 0) + Tout*np.maximum(xo - xs, 0)
    Nch = Fp/3 + float(np.abs(Mx).max())/lever
    rows.append(dict(tether_m=tl, push_MN=Fp/1e6, tip_MN=Pt/1e6, Ti_MN=Tin/1e6, To_MN=Tout/1e6, Mh_MNm=Mh/1e6, chord_ratio=Nch/1e6/t["N_chord_MN"], Mmax_MNm=float(np.abs(Mx).max())/1e6))
tom = max(r_["To_MN"] for r_ in rows)*1e6; A_draw = tom*TOPF/S75
for r_ in rows:
    P("    %5.1f m | %5.2f MN | %5.2f MN           | %5.2f MN   | %5.2f MN   | %-14s | %.2f                | %.2f km"
      % (r_["tether_m"], r_["push_MN"], r_["tip_MN"], r_["Ti_MN"], r_["To_MN"], ("tie %.2f MN" % (-r_["Mh_MNm"]/lever*1e6/1e6)) if r_["Mh_MNm"] < 0 else "free", r_["chord_ratio"], (r_["To_MN"]*1e6 - To_w)/(E_CNT*A_draw)*D/1e3))
ptm = max(r_["tip_MN"] for r_ in rows); crm = max(r_["chord_ratio"] for r_ in rows); mhm = min(r_["Mh_MNm"] for r_ in rows); tim = max(r_["Ti_MN"] for r_ in rows)
sink = (tom - To_w)/(E_CNT*A_draw)*D
m_draw = 2*RHO*A_draw*D
P("  Peak: %.2f MN on each tip ; outer line %.2f MN (%.2f at work) ; the joint pulled open with %.0f MN m, a tie of %.2f MN ; chord force %.2f times the sized working force (factor of safety %.1f)."
  % (ptm, tom/1e6, To_w/1e6, -mhm, -mhm/lever, crm, FOS))
P("  The statics pass.  The stretch does not.  An outer line sized for that peak (%.0f mm, %.0f t the pair) is %.2f km longer at the peak than at work, and the line is of fixed length:"
  % (dia(tom*TOPF, S75), m_draw/1e3, sink/1e3))
P("  the tip would sink kilometers under a load that is only there because the tip is %.0f m above the hanging edge.  Round 15 did not look at the stretch.  This way of drawing out" % (s_half - drop))
P("  does not work with suspension lines of fixed length.")
P("\n  Drawing the sheet out with the ridge let down instead.  The ridge's array support lines are on spools.  Start with the ridge high, so that the hanging edge is level with")
P("  the tips (%.0f m above the frame, not %.0f m).  Haul the tethers in level and pay the ridge out so that the edge stays level with the tips.  The tether is then horizontal" % (s_half, drop))
P("  all the way, the tip carries no weight, and the frame sees only the push as it grows.  The sheet at each stage is flat at its edge and slopes theta at the ridge:")
P("    slope at the ridge | edge from the center | ridge above the frame | push on the frame | of the final push")
lv = []
for th in (89.0, 85, 80, 75, 70, 65, 60, 55, 50, 45):
    tt_ = math.tan(math.radians(th)); y = s_half*math.asinh(tt_)/tt_; z = s_half*(math.sqrt(1 + tt_*tt_) - 1)/tt_; H = qs*s_half/tt_
    lv.append((th, y, z, H*S10))
    P("    %4.0f deg           | %6.1f m             | %6.1f m              | %6.2f MN         | %5.1f %%" % (th, y, z, H*S10/1e6, 100*H*S10/F))
thg = np.radians(np.linspace(89.999, 45.0, 40001)); tg = np.tan(thg); yg = s_half*np.arcsinh(tg)/tg; zg = s_half*(np.sqrt(1 + tg*tg) - 1)/tg; Hg = qs*s_half/tg
W_teth = float(np.sum(0.5*(Hg[1:] + Hg[:-1])*np.diff(yg)))*S10
def zc_above_edge(th):
    tt_ = math.tan(th); u = (np.arange(200000) + 0.5)/200000; psi = np.arctan(u*tt_); dz = np.sin(psi)*s_half/200000
    zz = np.cumsum(dz) - dz/2; return float(zz.mean())
dPE = qs*s_half*S10*(zc_above_edge(math.radians(45.0)) - s_half/2)                # the frame stays ; the sheet's center of mass comes down
E_ridge = W_teth - dPE
P("  Work of the tether winches: %.2f GJ for each half-sheet for each 10 km (round 14 and 15, ridge held: 2.94 GJ).  The ridge comes down %.0f m, and the ridge spools' brakes take %.2f GJ."
  % (W_teth/1e9, s_half - drop, E_ridge/1e9))
P("  No line and no chord sees more than its working load.  [This sequence is my proposal ; it is quasi-static, with the retention line taken as moving with the tips.]")
OUT["B"] = dict(lines=LN, lines_t=m_lines, lines_kgm=m_lines*1e3/S10, C_LINES=C_LINES, slack_len_km=[v/1e3 for v in L0], Zt_km=Zt/1e3, ri=ri, ro=ro,
                taut_above_m=Zt - taut["h"], turn_above_m=Zt - start["h"], flat_below_m=flat["h"] - Zt, takeup_km=(flat["h"] - taut["h"])/1e3, jump_deg=jump, Tc_min_MN=tc_min/1e6,
                E_brake_GJ=E_brake/1e9, flat=flat, work=dict(Ti_MN=Ti_w/1e6, To_MN=To_w/1e6, Tc_MN=Tc_w/1e6), travel_out_m=xo, travel_in_m=xi,
                seq=[s for s in seq if taut["h"] - 300 <= s["h"] <= flat["h"] + 100][::10],
                draw_ridge_held=dict(rows=rows, tip_max_MN=ptm, To_max_MN=tom/1e6, Ti_max_MN=tim, tie_MNm=-mhm, tie_MN=-mhm/lever, chord_ratio=crm, sink_km=sink/1e3, outer_mm=dia(tom*TOPF, S75), pair_t=m_draw/1e3),
                draw_level=dict(rows=lv, W_tether_GJ=W_teth/1e9, ridge_down_m=s_half - drop, E_ridge_GJ=E_ridge/1e9))
GB = "Ch12 Section 12.2 and 12.3 (lines 29 to 54), Table 11.2 (new frame-lines row)"
for l_ in LN:
    N("b", l_["name"] + ": pull open at work ; largest while opening ; sized for", "%.2f ; %.2f ; %.2f MN" % (l_["open_MN"], l_["opening_MN"], l_["size_MN"]), "%.4f ; %.4f ; %.4f MN" % (l_["open_MN"], l_["opening_MN"], l_["size_MN"]), GB)
    N("b", l_["name"] + ": across ; mass ; lean", "%.0f mm ; %.0f t ; %.2f deg" % (l_["d_mm"], l_["mass_t"], l_["lean_deg"]), "%.2f mm ; %.2f t ; %.4f deg" % (l_["d_mm"], l_["mass_t"], l_["lean_deg"]), GB)
N("b", "Working stress of the frame's lines", "18.75 GPa at the top", "18.75 GPa", GB + " ; as Section 12.4")
N("b", "The five lines of one frame ; for each meter of ring", "%.0f t ; %.1f kg/m" % (m_lines, m_lines*1e3/S10), "%.2f t ; %.3f kg/m" % (m_lines, m_lines*1e3/S10), GB)
N("b", "Stretch of a CNT line at 18.75 GPa over 100 km", "%.2f km (%.2f %%)" % (EPS75*100, EPS75*100), EPS75*100, "Ch12 Section 12.4 or 12.11 (new)")
N("b", "Sideways travel of the line points while opening (inner ; outer)", "%.0f m ; %.0f m" % (xi, xo), "%.2f ; %.2f m" % (xi, xo), "Ch12 Section 12.11")
N("b", "Hinge travel, stiff-line geometry (0.8 of the leg)", "%.0f m" % xo, LOG15["item1"]["0.800"]["travel_m"], "Ch12 Section 12.11", src="tools/round15_truss_array.py, log ITEM 1, lines at 0.800")
N("b", "Spool pay-out from the first pull in the suspension lines to the flat frame, lines as springs", "%.1f km" % ((flat["h"] - taut["h"])/1e3), (flat["h"] - taut["h"])/1e3, "Ch12 Section 12.11")
N("b", "Least pull of the center line during the opening", "%.2f MN" % (tc_min/1e6), tc_min/1e6, "Ch12 Section 12.11")
N("b", "Energy to the center line's brake, first pull to flat", "%.1f GJ" % (E_brake/1e9), E_brake/1e9, "Ch12 Section 12.11")
N("b", "Stop: moment that would open the joint at work ; pull across the upper chords ; shift of the thrust line under the push ; its limit", "%.0f MN m ; %.2f MN ; %.1f m ; %.1f m" % (-M0_MNm, -M0_MNm/lever, -M0_MNm*1e6/F, 0.2887*t["d"]),
  "%.2f MN m ; %.3f MN ; %.3f m ; %.2f m" % (-M0_MNm, -M0_MNm/lever, -M0_MNm*1e6/F, 0.2887*t["d"]), "Ch12 Section 12.3", key="A and B")
N("b", "Drawing out with the ridge held: peak weight on a tip ; outer line ; tie at the joint ; chord force over sized", "%.1f MN ; %.1f MN ; %.2f MN ; %.2f" % (ptm, tom/1e6, -mhm/lever, crm), "%.3f ; %.3f ; %.3f ; %.3f" % (ptm, tom/1e6, -mhm/lever, crm), "not for the book's design ; Ch12 Section 12.11 if the reason is given")
N("b", "Drawing out with the ridge held: the fixed outer line would be longer by", "%.1f km" % (sink/1e3), sink/1e3, "Ch12 Section 12.11 (why the ridge is let down)")
N("b", "Drawing out with the ridge let down: tether winch work for each half-sheet for each 10 km ; ridge comes down ; ridge brakes", "%.2f GJ ; %.0f m ; %.2f GJ" % (W_teth/1e9, s_half - drop, E_ridge/1e9), "%.4f GJ ; %.2f m ; %.4f GJ" % (W_teth/1e9, s_half - drop, E_ridge/1e9), "Ch12 Section 12.11")
E_open = sum(0.5*(a["Tc"] + b["Tc"])*(b["h"] - a["h"]) for a, b in zip(seq[:-1], seq[1:]) if a["h"] >= start["h"] and b["h"] <= flat["h"])
P("  Of the %.2f GJ to the center line's brake, %.2f GJ is taken while the legs turn (the last %.0f m) ; the rest is the take-up before it, while the spool still holds most of the frame."
  % (E_brake/1e9, E_open/1e9, flat["h"] - start["h"]))
OUT["B"]["E_open_GJ"] = E_open/1e9; OUT["B"]["turn_m"] = flat["h"] - start["h"]
N("b", "Spool pay-out while the legs turn from folded to flat ; brake energy in that stretch", "%.0f m ; %.2f GJ" % (flat["h"] - start["h"], E_open/1e9), "%.1f m ; %.4f GJ" % (flat["h"] - start["h"], E_open/1e9), "Ch12 Section 12.11")
# where the array settles
P("  Where the array settles.  With the slack lengths above the open frame at work is %.3f km below the top of the lines and the ridge is at 600.0 km: the lines are cut for it." % (Zt/1e3))
P("  If the outer lines were cut to a slack length of exactly 100 km, they would be %.2f km long at work, and the ridge would settle at %.1f km." % (100 + 100*To_w/(E_CNT*A[1]), 700 - (100e3*(1 + To_w/(E_CNT*A[1])) - drop)/1e3))
P("  [The top of the lines is taken at the casing, as round 9 did ; the force spreader cable hangs a few hundred meters lower, which shifts everything by that much.]")
OUT["B"]["settle_km_if_100km_slack"] = 700 - (100e3*(1 + To_w/(E_CNT*A[1])) - drop)/1e3
N("b", "Slack (cut) length of the suspension lines for a ridge at 600 km (inner ; outer)", "%.1f km ; %.1f km" % (L0[0]/1e3, L0[1]/1e3), "%.3f ; %.3f km" % (L0[0]/1e3, L0[1]/1e3), "Ch12 Section 12.11")
N("b", "Altitude of the ridge if the outer lines are cut to a slack 100 km and the casing ends at 700 km", "%.1f km" % OUT["B"]["settle_km_if_100km_slack"], OUT["B"]["settle_km_if_100km_slack"], "Ch12 Section 12.11 (brief, line 15)")

# =====================================================================================================================
head("C.  Yaw: what holds the frame against turning")
# =====================================================================================================================
yaw_up = F/leg; yaw_sheet = F/b_half; yaw_lines = (To_w*0.64 + Ti_w*0.16)/D
stay_T = BRACE*F; stay_len = math.hypot(S10, Lf/2); stay_m = 4*stay_T/S75*RHO*stay_len; stay_k = 2*E_CNT*(stay_T/S75)/stay_len*(S10/stay_len)**2
P("  A tip moved 1 m along the ring: the push on the turned strut upsets with %.0f N ; the sheet, whose pull points back at the ridge, restores %.0f N ; the lines %.1f N." % (yaw_up, yaw_sheet, yaw_lines))
P("  Net: neutral, as round 14 found.  The sheet and the 100 km lines do not hold the frame in yaw.")
P("  What holds it: four CNT stays, from each tip to the hinges of the next frames on either side, %.0f m long, each for %.0f %% of the push (%.2f MN) at 18.75 GPa: %.1f mm, %.2f t the four." % (stay_len, 100*BRACE, stay_T/1e6, dia(stay_T, S75), stay_m/1e3))
P("  The push and the sheet cancel, so the stays' %.0f N for each meter a tip moves along the ring is what is left to hold it (%.0f times the lines' share).  The 2 %% is the rule of thumb for bracing a strut [ASSUMED, no source opened]." % (stay_k, stay_k/yaw_lines))
OUT["C"] = dict(upset_N_m=yaw_up, sheet_N_m=yaw_sheet, lines_N_m=yaw_lines, stay_MN=stay_T/1e6, stay_len_m=stay_len, stay_mm=dia(stay_T, S75), stays_t=stay_m/1e3, stay_k_N_m=stay_k)
N("c", "Yaw stays: four for each frame, tip to the next frames' hinges ; each sized for ; across ; mass of the four", "%.2f MN (2 %% of the push) ; %.0f mm ; %.1f t" % (stay_T/1e6, dia(stay_T, S75), stay_m/1e3),
  "%.4f MN ; %.2f mm ; %.3f t" % (stay_T/1e6, dia(stay_T, S75), stay_m/1e3), "Ch12 Section 12.3 ; the 2 % rule is ASSUMED")
N("c", "Yaw: upsetting and restoring force for each meter a tip moves along the ring (push ; sheet ; lines ; stays)", "%.0f ; %.0f ; %.0f ; %.0f N/m" % (yaw_up, yaw_sheet, yaw_lines, stay_k), "%.1f ; %.1f ; %.2f ; %.1f" % (yaw_up, yaw_sheet, yaw_lines, stay_k), "Ch12 Section 12.3")

# =====================================================================================================================
head("D.  What the casing carries")
# =====================================================================================================================
s5 = LOG9["support"]["0.5"]
m_frame_sys = Mfr + m_lines*1e3 + stay_m; cas_kgm = m_frame_sys/S10
dPhi = pot(600e3, 700e3)
top_frame = (LN[0]["open_MN"] + 2*LN[1]["open_MN"] + 2*LN[2]["open_MN"])*1e6 + sum(RHO*l_["area_cm2"]*1e-4*dPhi*n_ for l_, n_ in zip(LN, (1, 2, 2)))
w_fsc = s5["T_top_kN_m"]*1e3 + (top_frame + stay_m*ga)/S10
T_ridge = s5["per_km_MN"]*1e6; T_center_top = LN[0]["open_MN"]*1e6 + RHO*LN[0]["area_cm2"]*1e-4*dPhi
F_line = 83.29e6/(Lr/50e3); lean_h = F_line/(w_fsc*50e3)
hang200 = w_fsc*200
P("  Frame %.1f t + five lines %.1f t + four stays %.2f t = %.1f t every 10 km: %.2f kg for each meter of ring (round 15, one line at 2/3: 62 to 65 ; round 14: 48.9 ; the book's arms: 572)." % (Mfr/1e3, m_lines, stay_m/1e3, m_frame_sys/1e3, cas_kgm))
P("  Pull at the top of the five lines of a frame, at work: %.3f MN (the frame's weight %.3f MN and the lines' own).  In the folded state the center line alone brings %.2f MN." % (top_frame/1e6, Mfr*ga/1e6, LN[0]["top_MN"]))
P("  Weight on the force spreader cable: array %.3f kN/m + frames, lines and stays %.3f kN/m = %.3f kN for each meter ; %.1f MN on 50 km (the book: 4.8 kN/m, 241 MN)." % (s5["T_top_kN_m"], (top_frame + stay_m*ga)/S10/1e3, w_fsc/1e3, w_fsc*50e3/1e6))
P("  The hangers' lean that passes the %.1f kN of a reaction line: %.2e of their length ; hangers of 500 m lean %.1f cm (the book: 4 x 10^-4, 20 cm)." % (F_line/1e3, lean_h, lean_h*500*100))
P("  An array support line every km carries %.2f MN at the top (18 mm).  The center array support line carries %.2f MN at the top at work and %.2f MN folded." % (T_ridge/1e6, T_center_top/1e6, LN[0]["top_MN"]))
P("  If the ridge's support line at a frame and the center line are one line, it carries %.2f MN at work, %.2f times an ordinary support line." % ((T_ridge + T_center_top)/1e6, (T_ridge + T_center_top)/T_ridge))
P("  Heaviest single attachment on the casing itself: a hanger of the force spreader cable.  Every 200 m it carries %.2f MN (the book, panels alone: 0.9 MN)." % (hang200/1e6))
P("  A trolley line or a guide line brings 2.26 MN, but in eight branches of 0.28 MN (Section 11.2).  Nothing like the arm's 23 MN is left.")
OUT["D"] = dict(frame_t=Mfr/1e3, lines_t=m_lines, stays_t=stay_m/1e3, casing_kgm=cas_kgm, top_frame_MN=top_frame/1e6, w_fsc_kN_m=w_fsc/1e3, per50km_MN=w_fsc*50e3/1e6, hanger_lean=lean_h,
                ridge_line_MN=T_ridge/1e6, center_top_work_MN=T_center_top/1e6, center_top_folded_MN=LN[0]["top_MN"], combined_MN=(T_ridge + T_center_top)/1e6, hanger200_MN=hang200/1e6, F_line_kN=F_line/1e3)
N("d", "Frame, its five lines and its four stays, for each meter of ring", "%.0f kg/m" % cas_kgm, cas_kgm, "Ch11 Table 11.2, Table 11.3 ; Ch12 line 125")
N("d", "Weight on the force spreader cable for each meter ; on 50 km", "%.1f kN/m ; %.0f MN" % (w_fsc/1e3, w_fsc*50e3/1e6), "%.4f kN/m ; %.2f MN" % (w_fsc/1e3, w_fsc*50e3/1e6), "Ch11 line 242 ; Ch12 lines 79, 263, 305")
N("d", "Lean of the hangers that passes 93.6 kN ; for hangers 500 m long", "%.1f x 10^-4 ; %.0f cm" % (lean_h*1e4, lean_h*5e4), "%.4e ; %.2f cm" % (lean_h, lean_h*5e4), "Ch11 line 242")
N("d", "Array support line every km, pull at the top (unchanged)", "%.1f MN, 18 mm" % (T_ridge/1e6), T_ridge/1e6, "Ch12 line 79", src="tools/round9_lines_array.py, log PART 2")
N("d", "Center array support line, pull at the top: at work ; folded (both legs)", "%.2f MN ; %.1f MN" % (T_center_top/1e6, LN[0]["top_MN"]), "%.4f ; %.4f MN" % (T_center_top/1e6, LN[0]["top_MN"]), "Ch12 Section 12.4")
N("d", "What a frame brings into the force spreader cable, five lines, at work", "%.1f MN" % (top_frame/1e6), top_frame/1e6, "Ch11 line 138 (was 23 MN at each arm) ; Ch12 line 305")
N("d", "Heaviest single attachment on the casing: a hanger every 200 m", "%.1f MN" % (hang200/1e6), hang200/1e6, "Ch11 line 138 ; Ch12 line 81")

# =====================================================================================================================
head("E and F.  The budget of the design case, and what 1,200 kg/m buys (two lines, the outer at 4/5)")
# =====================================================================================================================
def fs16(sA, W, S):
    """Round 15's frame with two lines at 2/5 and 4/5, and the lines counted as section B sizes them (the same shares of the frame's weight)."""
    x = dict(frame_system(sA, W, S, n=2, f=0.8)); m = x["mass_t"]*1e3
    x["lines_t"] = kc*C_LINES*m/1e3; x["casing_kgm"] = (m + kc*C_LINES*m + x["stays_t"]*1e3)/S
    return x
PW = lambda sA: sA*La/Lr*(1 + kc)
def budget16(sA, S, M=1200.0):
    def tot(W):
        x = fs16(sA, W, S); return PW(sA)*W + x["ret_kgm"] + x["casing_kgm"] - M
    if tot(50.0) > 0: return None
    W = brentq(tot, 50.0, M/PW(sA), xtol=0.02); x = fs16(sA, W, S)
    return dict(S=S, W=W, pan=PW(sA)*W, ret=x["ret_kgm"], frames=x["casing_kgm"], P_TW=BOOKENV.TW_PER_KM*W/1e3, f=x)
x5 = fs16(0.5, W5, S10); pan5 = PW(0.5)*W5; tot5 = pan5 + x5["ret_kgm"] + x5["casing_kgm"]; P5 = BOOKENV.TW_PER_KM*W5/1e3
P("  Design case: panels and array support lines %.1f + retention lines %.1f (%.0f MN, %.0f mm) + frames with lines and stays %.1f = %.1f kg for each meter of ring."
  % (pan5, x5["ret_kgm"], x5["TR_MN"], x5["ret_mm"], x5["casing_kgm"], tot5))
P("  Left of the 1,200 kg/m: %.0f kg/m.  Power %.2f TW ; %.1f kg/m for each TW (the book with ground arms: 98).  The bow of the retention line is %.0f m." % (1200 - tot5, P5, tot5/P5, W5/8))
OUT["E"] = dict(pan=pan5, ret=x5["ret_kgm"], TR_MN=x5["TR_MN"], ret_mm=x5["ret_mm"], frames=x5["casing_kgm"], total=tot5, left=1200 - tot5, P_TW=P5, kgm_per_TW=tot5/P5)
GE = "Ch12 line 125 ; Ch11 line 196 and Table 11.3 last rows"
N("e", "Design case: panels and array support lines", "%.0f kg/m" % pan5, pan5, GE)
N("e", "Design case: retention lines (tension, across, both for each meter)", "%.0f MN, %.0f mm, %.0f kg/m" % (x5["TR_MN"], x5["ret_mm"], x5["ret_kgm"]), "%.2f MN, %.2f mm, %.3f kg/m" % (x5["TR_MN"], x5["ret_mm"], x5["ret_kgm"]), GE + " ; Ch12 line 50")
N("e", "Design case: frames with their lines and stays", "%.0f kg/m" % x5["casing_kgm"], x5["casing_kgm"], GE)
N("e", "Design case: total ; left of 1,200 kg/m", "%.0f kg/m ; %.0f kg/m" % (tot5, 1200 - tot5), "%.2f ; %.2f" % (tot5, 1200 - tot5), GE)
N("e", "Design case: power ; kg/m for each terawatt", "%.1f TW ; %.0f kg/m" % (P5, tot5/P5), "%.3f TW ; %.2f" % (P5, tot5/P5), GE)
OUT["F"] = {}
RNG = {1.0: (2e3, 14e3), 0.5: (4e3, 28e3), 0.2: (8e3, 70e3)}
BOOKROW = {1.0: (0.54, 5.9), 0.5: (1.11, 12.2), 0.2: (2.83, 31.1)}
P("\n    panels    | spacing            | sheet    | footprint | power    | frame: mass, length, depth    | push     | panels + retention + frames | GW for each frame | kg/m for each TW | against the best")
for sA in (1.0, 0.5, 0.2):
    grid = list(np.geomspace(RNG[sA][0], RNG[sA][1], 15)); bb = [b for b in (budget16(sA, S) for S in grid) if b]
    best = max(bb, key=lambda b: b["W"])
    for S in np.geomspace(best["S"]/1.14, best["S"]*1.14, 9):
        b = budget16(sA, S)
        if b and b["W"] > best["W"]: best = b
    b10 = budget16(sA, 10e3)
    OUT["F"][str(sA)] = {}
    for nm, b in (("best", best), ("10 km", b10)):
        x = b["f"]
        P("    %.1f kg/m2 | %-5s %5.1f km     | %.3f km | %.3f km  | %5.2f TW | %6.0f t, %5.0f m, %5.1f m     | %5.1f MN | %.0f + %.0f + %.0f            | %5.2f GW           | %5.1f            | %.1f %%"
          % (sA, nm, b["S"]/1e3, b["W"]/1e3, HS45*b["W"]/1e3, b["P_TW"], x["mass_t"], x["Lf"], x["truss"]["d"], x["F"]/1e6, b["pan"], b["ret"], b["frames"], b["P_TW"]*1e3/(Lr/b["S"]), 1200/b["P_TW"], 100*(b["P_TW"]/best["P_TW"] - 1)))
        OUT["F"][str(sA)][nm] = dict(S_km=b["S"]/1e3, sheet_km=b["W"]/1e3, footprint_km=HS45*b["W"]/1e3, P_TW=b["P_TW"], frame_t=x["mass_t"], Lf=x["Lf"], depth_m=x["truss"]["d"], push_MN=x["F"]/1e6, pan=b["pan"], ret=b["ret"],
                                     frames=b["frames"], GW_per_frame=b["P_TW"]*1e3/(Lr/b["S"]), kgm_per_TW=1200/b["P_TW"], drop_m=DR45*b["W"]/2, ret_MN=x["TR_MN"], ret_mm=x["ret_mm"], vs_best=b["P_TW"]/best["P_TW"] - 1)
        lab = "best spacing" if nm == "best" else "frames every 10 km"
        N("f", "Table 12.1, %.1f kg/m2, %s: spacing ; sheet ; footprint ; power" % (sA, lab), "%.0f km ; %.2f km ; %.2f km ; %.1f TW" % (b["S"]/1e3, b["W"]/1e3, HS45*b["W"]/1e3, b["P_TW"]),
          "%.2f km ; %.4f km ; %.4f km ; %.3f TW" % (b["S"]/1e3, b["W"]/1e3, HS45*b["W"]/1e3, b["P_TW"]), "Ch12 Table 12.1 (lines 117 to 125), lines 7, 283, 304 ; Ch8 line 46")
        N("f", "Table 12.1, %.1f kg/m2, %s: frame mass, length, depth ; power for each frame ; kg/m for each TW" % (sA, lab),
          "%.0f t, %.0f m, %.0f m ; %.1f GW ; %.0f kg/m" % (x["mass_t"], x["Lf"], x["truss"]["d"], b["P_TW"]*1e3/(Lr/b["S"]), 1200/b["P_TW"]), "%.1f t, %.1f m, %.2f m ; %.3f GW ; %.2f" % (x["mass_t"], x["Lf"], x["truss"]["d"], b["P_TW"]*1e3/(Lr/b["S"]), 1200/b["P_TW"]), "Ch12 Table 12.1, line 125")
    ok10 = bool(b10["P_TW"] > 0.98*best["P_TW"])
    P("      book, ground arms: %.2f km, %.1f TW.   10 km is %.1f %% under the best spacing: %s" % (BOOKROW[sA][0], BOOKROW[sA][1], 100*(1 - b10["P_TW"]/best["P_TW"]), "the book can say 10 km" if ok10 else "NOT within 2 %"))
    OUT["F"][str(sA)]["ten_within_2pct"] = ok10
    N("f", "Table 12.1, %.1f kg/m2: 10 km against the best spacing" % sA, "%.1f %% less power ; %s" % (100*(1 - b10["P_TW"]/best["P_TW"]), "within 2 %" if ok10 else "not within 2 %"), 1 - b10["P_TW"]/best["P_TW"], "Ch12 Table 12.1")

# =====================================================================================================================
head("G.  What sets the size of a frame (the new Table 11.1): sheet 1,110 m wide, a frame every 10 km")
# =====================================================================================================================
def frame_case(sA, W, S, theta, bow):
    hs, dr = catenary(theta); tt_ = math.tan(math.radians(theta)); q = sA*ga; H = q*W/2/tt_; Fp = H*S; lg = hs*W/2; Dz = D + dr*W/2
    m = 0.0; tr = None
    for _ in range(3):
        add = sum(Rk*(m*ga/2)*x_*lg/Dz for Rk, x_ in zip(lay["R"], lay["x"]))
        tr = truss(Fp + add, 2*lg, 2*lg, lay["cM"], lay["cV"], ga, None, detail=True)
        if abs(tr["mass"] - m) < 1e-3*tr["mass"]: m = tr["mass"]; break
        m = tr["mass"]
    TR = H*S*S/(8*bow*W)
    return dict(push_MN=Fp/1e6, frame_t=m/1e3, depth_m=tr["d"], Lf=2*lg, TR_MN=TR/1e6, ret_mm=dia(TR, S75), ret_kgm=2*RHO*TR/S75, drop_m=dr*W/2, footprint_m=2*lg, chord_mm=2e3*tr["r_chord"], wall_mm=1e3*tr["t_chord"])
OUT["G"] = {}
P("    case                                   | push on a frame | frame: mass, depth, length  | retention line     | both, for each meter | edge below the ridge | (the book's arm)")
BK = {"Design case": "22.6 MN, 2,100 t", "Lighter panels (0.2 kg/m2)": "9.0 MN, 840 t", "Steeper sheet (60 deg at the ridge)": "13.1 MN, 1,210 t", "Deeper bow (a quarter of the width)": "22.6 MN, 2,100 t", "All three together": "5.2 MN, 480 t"}
for nm, a, th, b in (("Design case", 0.5, 45, 1/8), ("Lighter panels (0.2 kg/m2)", 0.2, 45, 1/8), ("Steeper sheet (60 deg at the ridge)", 0.5, 60, 1/8), ("Deeper bow (a quarter of the width)", 0.5, 45, 1/4), ("All three together", 0.2, 60, 1/4)):
    c = frame_case(a, 1110.0, S10, th, b); OUT["G"][nm] = c
    P("    %-38s | %5.1f MN        | %5.0f t, %4.1f m, %4.0f m     | %4.0f MN, %4.0f mm    | %5.1f kg/m           | %4.0f m               | %s"
      % (nm, c["push_MN"], c["frame_t"], c["depth_m"], c["Lf"], c["TR_MN"], c["ret_mm"], c["ret_kgm"], c["drop_m"], BK[nm]))
    N("g", "Table 11.1, %s: push ; frame mass and depth ; retention line ; retention lines for each meter ; edge below the ridge" % nm,
      "%.1f MN ; %.0f t, %.0f m ; %.0f MN, %.0f mm ; %.0f kg/m ; %.0f m" % (c["push_MN"], c["frame_t"], c["depth_m"], c["TR_MN"], c["ret_mm"], c["ret_kgm"], c["drop_m"]),
      "%.3f MN ; %.1f t, %.2f m ; %.2f MN, %.2f mm ; %.3f kg/m ; %.1f m" % (c["push_MN"], c["frame_t"], c["depth_m"], c["TR_MN"], c["ret_mm"], c["ret_kgm"], c["drop_m"]), "Ch11 Table 11.1 (lines 30 to 40), line 28")
P("  The frame is as long as the footprint: %.0f m at 45 deg, %.0f m at 60 deg.  The bow does not change the frame, only the retention line, as the book says of the arm." % (OUT["G"]["Design case"]["Lf"], OUT["G"]["Steeper sheet (60 deg at the ridge)"]["Lf"]))
P("  (The table's sheet is the book's round 1,110 m ; the design case of section A is 1,111.1 m, which is why its frame is %.1f t there and %.1f t here.)" % (Mfr/1e3, OUT["G"]["Design case"]["frame_t"]))

# =====================================================================================================================
head("H.  The ring's own lines: the guide lines, the reaction lines, the frame lines (Tables 11.2 and 11.3 without arms)")
# =====================================================================================================================
CH15 = os.path.normpath(os.path.join(ROOT, "..", "..", "Vol_III", "renumber_scripts_2026-10-12", "ch15_numbers.py"))
ch15_out = ""
if os.path.exists(CH15):
    ch15_out = subprocess.run([sys.executable, CH15], cwd=ROOT, capture_output=True, text=True, timeout=600).stdout
    with open(os.path.join(LOG_DIR, "ch15_numbers_stdout.txt"), "w", encoding="utf-8") as fh: fh.write(ch15_out)
def grab(pat):
    m = re.search(pat, ch15_out); return m.group(0) if m else "(chapter script not found on this machine)"
P("  Table 11.3's 'without the array' column comes from Vol_III/renumber_scripts_2026-10-12/ch15_numbers.py (the chapter script, outside the repository).")
P("  It was run again as it stands ; its output is saved as runs/round16_logs/ch15_numbers_stdout.txt.  Its lines for that column:")
for pat in (r"100 t tip load:.*", r" two arms as 100 t lines.*", r"vertical 100 t line, sideways stiffness.*", r"swing period.*", r"reaction line: stretch.*"):
    P("    " + grab(pat)[:230])
# the same again here
dP_g = pot(0, 700e3); TAPER = math.exp(RHO*dP_g/S_OP); g0 = geff(0)
A_tip = 100e3*g0/S_OP; nq = 20000; zz = (np.arange(nq) + 0.5)*700e3/nq; ex = np.exp(RHO*np.array([pot(0, z) for z in zz])/S_OP)
m100 = float((RHO*A_tip*ex).sum()*700e3/nq); T100 = S_OP*A_tip*TAPER; k_p = 1.0/float((1/(S_OP*A_tip*ex)).sum()*700e3/nq)
n_st = Lr/50e3; n_fr = Lr/S10
th30 = math.radians(30); k_r_book = F_line/((S_OP/E_CNT)*(700e3/math.cos(th30))/math.sin(th30))
k_r_hang = LOG9["reaction"]["spring"]["B. anchorage 404 km away"]["k_N_per_m"]
RB = LOG9["reaction"]["figures"]["B. anchorage 404 km away (700 tan 30)"]
P("  Recomputed here: a 100 t line is %.1f mm at the ground and %.1f mm at the ring, %.1f t, %.2f MN at the top ; a pair every 50 km is %.2f kg/m ; %d guide lines and as many anchorages."
  % (dia(100e3*g0, S_OP), dia(T100, S_OP), m100/1e3, T100/1e6, 2*m100/50e3, 2*round(n_st)))
P("  Confirmed: 10 mm, 148 t, 1,778, 5.9 kg/m.  What the guide lines are sized for: the book gives no load.  Section 11.5 says only 'each a 100 tonne line like a trolley line'.")
P("  The chapter script takes the trolley line's 100 t tip load for them and nothing else.  They are not sized for a north or south force anywhere in the book or its scripts.")
P("  Anchorages: %d guide lines + %d reaction lines = %d on the ground (the book with arms: 8,894 + 889)." % (2*round(n_st), round(n_st), 3*round(n_st)))
P("  Reaction line (round 9, hanging, anchorage 404 km away): %.1f mm at the ground, %.1f mm at the ring, %.0f km, %.1f t, %.2f kg/m ; all %d: %.0f t." % (RB["d_g"], RB["d_top"], RB["length_km"], RB["mass_t"], RB["kgm"], round(n_st), RB["all_t"]))
P("  Frame lines: %d frames, five lines each, %d lines ; they go to the force spreader cable, not to the ground ; 100 km ; lean %.2f and %.2f deg ; %.0f, %.0f and %.0f mm ; %.1f, %.1f and %.1f t ; %.2f kg/m."
  % (round(n_fr), round(5*n_fr), LN[1]["lean_deg"], LN[2]["lean_deg"], LN[0]["d_mm"], LN[1]["d_mm"], LN[2]["d_mm"], LN[0]["mass_t"], LN[1]["mass_t"], LN[2]["mass_t"], m_lines*1e3/S10))
lines_ground = 2*m100/50e3 + RB["kgm"]
P("  The ring's lines to the ground: guide lines %.2f + reaction lines %.2f = %.2f kg/m, plus the trolley lines' 0.2 to 0.3 kg/m.  With the array: + frames, lines and stays %.1f kg/m, which hang from the casing and do not reach the ground." % (2*m100/50e3, RB["kgm"], lines_ground, cas_kgm))
OUT["H"] = dict(m100_t=m100/1e3, d_ground_mm=dia(100e3*g0, S_OP), d_ring_mm=dia(T100, S_OP), T100_MN=T100/1e6, guide_kgm=2*m100/50e3, n_guide=2*n_st, n_reaction=n_st, anchorages=3*n_st, n_frames=n_fr, n_frame_lines=5*n_fr,
                reaction=RB, k_p=k_p, k_r_book=k_r_book, k_r_hang=k_r_hang, lines_ground_kgm=lines_ground)
GH = "Ch11 Table 11.2 (lines 44 to 51), Table 11.3 (lines 196 to 214), line 116"
N("h", "Guide line (a pair every 50 km): at the ground ; at the ring ; length ; mass", "%.0f mm ; %.0f mm ; 700 km ; %.0f t" % (dia(100e3*g0, S_OP), dia(T100, S_OP), m100/1e3), "%.2f mm ; %.2f mm ; %.2f t" % (dia(100e3*g0, S_OP), dia(T100, S_OP), m100/1e3), GH, src="ch15_numbers.py (chapter script), run again, and " + SCRIPT + " section H")
N("h", "Guide lines: how many ; for each meter of ring", "%d ; %.1f kg/m" % (2*round(n_st), 2*m100/50e3), "%.1f ; %.3f kg/m" % (2*n_st, 2*m100/50e3), GH)
N("h", "Load the guide lines are sized for", "none stated in the book ; the script uses the trolley line's 100 t tip load", "n/a", GH)
N("h", "Anchorages on the ground: guide lines ; reaction lines ; all", "%d ; %d ; %d" % (2*round(n_st), round(n_st), 3*round(n_st)), "%.1f ; %.1f ; %.1f" % (2*n_st, n_st, 3*n_st), GH)
N("h", "Frame lines: how many ; lean (inner, outer) ; across (center, inner, outer) ; length ; mass (center, inner, outer) ; for each meter", "%d (five at each of %d frames) ; %.2f and %.2f deg ; %.0f, %.0f, %.0f mm ; 100 km ; %.0f, %.0f, %.0f t ; %.1f kg/m"
  % (round(5*n_fr), round(n_fr), LN[1]["lean_deg"], LN[2]["lean_deg"], LN[0]["d_mm"], LN[1]["d_mm"], LN[2]["d_mm"], LN[0]["mass_t"], LN[1]["mass_t"], LN[2]["mass_t"], m_lines*1e3/S10), "see section b", GH + " (new row)", key="B and H")
N("h", "Table 11.3, with the array: frames, their lines and stays ; array in all, design case", "%.0f kg/m ; %.0f kg/m" % (cas_kgm, tot5), "%.2f ; %.2f" % (cas_kgm, tot5), GH, key="D and E")

# =====================================================================================================================
head("I.  The swing along the ring (Section 11.6)")
# =====================================================================================================================
M_CAS = 12000.0*Lr; J_FULL = 10.25e6*30000.0; J_NET = 5.0e6*30000.0; N_SITES = 88946; F_LIM = 2*401.2*N_SITES; V_REL = 8846.96
k_script = n_st*(k_r_book + 2*k_p); k_hang = n_st*(k_r_hang + 2*k_p); K_RING = 5800.0
P("  The 5,800 N/m is from ch15_numbers.py (added 13 and 16 October): for each of the %.0f anchor points, one reaction line and two vertical 100 t lines." % n_st)
P("    reaction line, as a straight line at 30 deg that stretches: %.3f N for each meter the casing moves" % k_r_book)
P("    vertical 100 t line, as a pendulum (1 / the sum of ds / T up the line): %.3f N/m ; two of them %.3f N/m" % (k_p, 2*k_p))
P("    ring: %.1f x (%.3f + %.3f) = %.0f N/m, printed in the book as 5,800." % (n_st, k_r_book, 2*k_p, k_script))
P("    With round 9's hanging reaction line (%.3f N/m, the book's '1.6 N'): %.0f N/m.  The guide lines give %.0f %% of it." % (k_r_hang, k_hang, 100*2*k_p/(k_r_hang + 2*k_p)))
P("    The frames' lines add nothing along the ring to the casing: they hang from it and move with it.")
def swing(k): om = math.sqrt(k/M_CAS); return 2*math.pi/om, J_FULL/M_CAS/om, J_NET/M_CAS/om, om
for nm, k in (("5,800 N/m (the book's figure)", K_RING), ("%.0f N/m (round 9's reaction line)" % k_hang, k_hang)):
    T_, a1, a2, om = swing(k)
    P("  %-38s casing %.3e kg ; period %.0f s (%.2f h) ; swing after one design launch %.0f m, %.0f m with the sled braked ; a reaction line then %.2f kN of its %.1f kN"
      % (nm, M_CAS, T_, T_/3600, a1, a2, k_r_hang*a1/1e3, F_line/1e3))
T_, a1, a2, om0 = swing(K_RING); C_CRIT = 2*math.sqrt(K_RING*M_CAS)
def one_launch(J, c, pulse=3690.0, days=12, dt=10.0):
    x = v = 0.0; xm = 0.0; fm = 0.0; t_die = None; n = int(days*86400/dt); xs = []
    for i in range(n):
        tt_ = i*dt; Fl = J/pulse if tt_ < pulse else 0.0; Fd = max(-F_LIM, min(F_LIM, -c*v)); a = (Fl + Fd - K_RING*x)/M_CAS; v += a*dt; x += v*dt
        xs.append(x); fm = max(fm, abs(Fd))
    xs = np.array(xs); xm = float(np.abs(xs).max()); late = np.where(np.abs(xs) > 0.01*xm)[0]
    return xm, fm, (late[-1] + 1)*dt/3600.0, float(np.argmax(np.abs(xs)))*dt/3600.0
OUT["I"] = dict(k_script=k_script, k_hang=k_hang, k_r_book=k_r_book, k_r_hang=k_r_hang, k_p=k_p, period_h=T_/3600, swing_m=a1, swing_braked_m=a2, M_cas=M_CAS, c_crit=C_CRIT, one={})
P("  One design launch with the LIMs damping at the critical rate (%.2e N s/m), push spread over 3,690 s [ASSUMED shape, as round 15]:" % C_CRIT)
for nm, J in (("not braked", J_FULL), ("braked", J_NET)):
    xm, fm, tdie, tpk = one_launch(J, C_CRIT)
    P("    sled %-10s: farthest %.0f m, %.1f h after the start ; back within 1 %% of that after %.0f h ; LIM push at most %.1f MN of %.1f MN, %.0f GW ; reaction line %.2f kN"
      % (nm, xm, tpk, tdie, fm/1e6, F_LIM/1e6, fm*V_REL/1e9, k_r_hang*xm/1e3))
    OUT["I"]["one"][nm] = dict(x_max=xm, F_lim_MN=fm/1e6, die_h=tdie, peak_h=tpk, line_kN=k_r_hang*xm/1e3, GW=fm*V_REL/1e9)
P("  Seven and fourteen a day: round 15's rows (tools/round15_truss_array.py, log ITEM 5), copied:")
P("    launches | sled       | LIMs                               | casing sits at | swings about it | reaction line | LIM push peak | power")
R5 = [r_ for r_ in LOG15["item5"]["rows"]]
for r_ in R5:
    P("    %2d a day | %-10s | %-34s | %7.0f m      | %6.0f m        | %5.1f kN      | %5.1f MN      | %4.0f GW" % (r_["N"], r_["sled"], r_["lims"][:34], r_["x_mean"], r_["swing"], r_["line_kN"], r_["F_lim_peak_MN"], r_["F_lim_peak_MN"]*1e6*V_REL/1e9))
dvc = LOG15["item5"]["dv_cable_per_launch"]
P("  If the LIMs also hold the mean push, the cable gains %.4f m/s for each braked launch, %.0f m/s a year at seven a day (round 15)." % (dvc, 7*dvc*365.25))
OUT["I"]["round15_rows"] = R5; OUT["I"]["dv_cable"] = dvc
GI = "Ch11 Section 11.6 (lines 252 to 261), Table 11.3 last row, line 293"
N("i", "Hold along the ring with no arms: the terms", "889 x (%.2f + 2 x %.2f) = %.0f N/m, printed 5,800" % (k_r_book, k_p, k_script), "%.4f ; %.4f ; %.1f" % (k_r_book, k_p, k_script), GI, src="ch15_numbers.py (chapter script), run again ; recomputed in " + SCRIPT + " section I")
N("i", "The same with the hanging reaction line of round 9 (1.55 N/m, the book's 1.6)", "%.0f N/m" % k_hang, k_hang, GI)
N("i", "Casing mass ; period on 5,800 N/m", "%.2f x 10^11 kg ; %.1f h" % (M_CAS/1e11, T_/3600), "%.4e kg ; %.1f s" % (M_CAS, T_), GI)
N("i", "Swing after one design launch, no damping: not braked ; braked", "%.1f km ; %.1f km" % (a1/1e3, a2/1e3), "%.1f m ; %.1f m" % (a1, a2), GI)
for nm in ("not braked", "braked"):
    o = OUT["I"]["one"][nm]
    N("i", "One design launch, LIMs damping at the critical rate, sled %s: farthest ; time for the swing to die to 1 %% ; LIM push ; reaction line" % nm, "%.1f km ; %.0f h ; %.0f MN ; %.1f kN" % (o["x_max"]/1e3, o["die_h"], o["F_lim_MN"], o["line_kN"]),
      "%.1f m ; %.2f h ; %.2f MN ; %.3f kN" % (o["x_max"], o["die_h"], o["F_lim_MN"], o["line_kN"]), GI)
for r_ in R5:
    if r_["lims"].startswith("damping, critical"):
        N("i", "%d launches a day, sled %s, LIMs damping: casing sits at ; swings about it ; reaction line ; LIM push at peak" % (r_["N"], r_["sled"]), "%.1f km ; %.0f m ; %.1f kN ; %.0f MN" % (r_["x_mean"]/1e3, r_["swing"], r_["line_kN"], r_["F_lim_peak_MN"]),
          "%.1f m ; %.1f m ; %.3f kN ; %.2f MN" % (r_["x_mean"], r_["swing"], r_["line_kN"], r_["F_lim_peak_MN"]), GI, src="tools/round15_truss_array.py, log ITEM 5")
N("i", "LIMs' push available ; power at full push", "%.1f MN ; %.0f GW" % (F_LIM/1e6, F_LIM*V_REL/1e9), "%.2f MN ; %.1f GW" % (F_LIM/1e6, F_LIM*V_REL/1e9), GI, src="tools/round15_truss_array.py, log ITEM 5")
N("i", "Cable's gain of speed if the LIMs hold the mean push: for each braked launch ; a year at seven a day", "%.3f m/s ; %.0f m/s" % (dvc, 7*dvc*365.25), "%.5f ; %.2f" % (dvc, 7*dvc*365.25), GI, src="tools/round15_truss_array.py, log ITEM 5")

# =====================================================================================================================
head("J.  The conductors (Sections 11.7 and 12.10)")
# =====================================================================================================================
SIG = 5.8e7; h_r = 700e3
def cond(P_total, n_pts, V=1e6, loss=0.02):
    """ch16_numbers.py, function anchor: two conductors of one thickness for each anchor point, 700 km down and 700 km back, a loss of 'loss' of the power."""
    Pp = P_total/n_pts; I = Pp/V; A_ = 2*h_r*I/(SIG*loss*V); m = 2*RHO*A_*h_r
    return dict(P_MW=Pp/1e6, I_A=I, A_cm2=A_*1e4, d_mm=2e3*math.sqrt(A_/math.pi), t_pair=m/1e3, kgm=m*n_pts/Lr)
OUT["J"] = {}
A_guide = (A_tip*1e4, A_tip*TAPER*1e4)
P("  CNT at the conductivity of copper (%.0f MS/m) ; a loss of 2 %% over 700 km down and 700 km back ; one pole pair at each of the %.0f anchor points (the two guide lines)." % (SIG/1e6, n_st))
P("    power             | voltage | at each anchor point | current  | each conductor      | pair      | for each meter of ring | for each TW")
for nm, Pt_, V in (("ring panels' export, 900 GW", 900e9, 1e6), ("array, %.1f TW" % P5, P5*1e12, 1e6), ("array, %.1f TW" % P5, P5*1e12, 2e6)):
    c = cond(Pt_, n_st, V); OUT["J"]["%s at %.0f MV" % (nm, V/1e6)] = c
    P("    %-17s | %.0f MV    | %7.0f MW           | %6.0f A | %6.1f cm2, %5.1f mm | %6.0f t  | %6.1f kg/m            | %.1f kg/m" % (nm[:17], V/1e6, c["P_MW"], c["I_A"], c["A_cm2"], c["d_mm"], c["t_pair"], c["kgm"], c["kgm"]/(Pt_/1e12)))
c9 = OUT["J"]["ring panels' export, 900 GW at 1 MV"]; ca1 = OUT["J"]["array, %.1f TW at 1 MV" % P5]; ca2 = OUT["J"]["array, %.1f TW at 2 MV" % P5]
P("  A guide line has %.2f cm2 of CNT at the ground and %.2f cm2 at the ring for its structural job.  The export of the ring-mounted panels asks for %.1f cm2 in each: %.0f times more." % (A_guide[0], A_guide[1], c9["A_cm2"], c9["A_cm2"]/A_guide[1]))
P("  So nearly all of the %.0f kg/m is added mass (the guide lines themselves are %.1f kg/m).  'Adds no mass to the arms at all' (Section 11.7) no longer holds." % (c9["kgm"], 2*m100/50e3))
P("  The mass for each meter of ring does not depend on how many lines share the current: with a pair every 10 km in place of every 50 km the conductors are a fifth as thick and five times as many.")
P("  The book's rows: 789 and 197 kg/m, 65 kg/m for each TW (round 9's log: %.1f and %.1f).  Confirmed.  Each anchor point of a pole pair every 50 km handles %.1f GW of the %.1f TW ; the book's 2.7 GW was for each spreader." % (P9["0.5"]["anchor_1MV"], P9["0.5"]["anchor_2MV"], ca1["P_MW"]/1e3, P5))
P("  Chapter 11 itself is not consistent on the export: line 269 has 900 GW over 4,447 lines as 'about 200 MW', and line 279 then uses 135 MW ; the chapter script used 659 GW over 889 stations.  Chapter 8 line 124 has 900 GW over 889 stations, 1,010 MW.")
GJ = "Ch11 Section 11.7 (lines 269, 277, 279), line 293 ; Ch12 Section 12.10 (line 229), Table 12.2"
N("j", "Export of the ring-mounted panels, 900 GW at 1 MV: at each anchor point ; current ; each conductor ; added for each meter of ring", "%.0f MW ; %.0f A ; %.0f cm2, %.0f mm ; %.0f kg/m" % (c9["P_MW"], c9["I_A"], c9["A_cm2"], c9["d_mm"], c9["kgm"]),
  "%.1f MW ; %.1f A ; %.3f cm2, %.2f mm ; %.2f kg/m" % (c9["P_MW"], c9["I_A"], c9["A_cm2"], c9["d_mm"], c9["kgm"]), GJ)
N("j", "CNT in a guide line for its structural job: at the ground ; at the ring", "%.1f cm2 ; %.1f cm2" % A_guide, "%.3f ; %.3f cm2" % A_guide, GJ)
N("j", "Array, all %.1f TW by wire, anchor lines at 1 MV ; at 2 MV ; for each TW at 1 MV" % P5, "%.0f kg/m ; %.0f kg/m ; %.0f kg/m" % (ca1["kgm"], ca2["kgm"], ca1["kgm"]/P5), "%.2f ; %.2f ; %.3f" % (ca1["kgm"], ca2["kgm"], ca1["kgm"]/P5), GJ + " (confirmed, unchanged)")
N("j", "Array, all by wire at 1 MV, a pole pair every 50 km: at each anchor point ; each conductor", "%.1f GW ; %.0f mm" % (ca1["P_MW"]/1e3, ca1["d_mm"]), "%.3f GW ; %.2f mm" % (ca1["P_MW"]/1e3, ca1["d_mm"]), GJ + " (was 2.7 GW and 65 mm for each spreader)")

# =====================================================================================================================
head("K, L, M.  Why not from the ground ; the clearance of Section 12.7 ; Chapter 16 line 142")
# =====================================================================================================================
def find13(pat):
    for i, ln in enumerate(LOG13T.splitlines(), 1):
        if re.search(pat, ln): return i, ln.strip()
    return 0, "(not found)"
P("  K.  No rerun.  Round 13's log, runs/round13_logs/round13_spreader_arms.txt (tools/round13_spreader_arms.py):")
for pat in (r"leans at most 26\.19 deg", r"pull carried up to the casing 47\.3 MN", r"11949 t per spreader, 1195 kg/m", r"2\.317 MN at the kink"):
    i, ln = find13(pat); P("    line %d: %s" % (i, ln[:250]))
N("k", "A hanging arm at 12.5 GPa leans at most, at the kink", "26.2 deg", "26.19 deg", "Ch12 Section 12.3 (one paragraph)", src="tools/round13_spreader_arms.py, log round13_spreader_arms.txt line %d" % find13(r"leans at most 26\.19 deg")[0])
N("k", "Pull on the casing for each newton of sideways pull ; tension at the kink for each newton", "2.09 N ; 2.32 N", "2.090 ; 2.317", "Ch12 Section 12.3", src="tools/round13_spreader_arms.py, log line %d" % find13(r"2\.317 MN at the kink")[0])
N("k", "Hanging arms on the casing, design case", "1,195 kg/m", "1195 kg/m (11,949 t for each spreader)", "Ch12 Section 12.3", src="tools/round13_spreader_arms.py, log line %d" % find13(r"11949 t per spreader")[0])
SUN = math.radians(0.533); W_CAS = 10.0
P("\n  L.  The casing's width.  The book gives no settled figure: Chapter 3 (line 238) has 'approximately ... 9 to 11 m per side for the retrograde cable, but these numbers are not yet")
P("      finalized' ; Figure 1.1 prints 'about 9 to 11 m on a side' ; Chapter 14 (line 100) takes 'a 10 m wide casing' for its heat.  Round 15's 12 m was my assumption and is dropped.")
P("      Full shadow of a body of width b ends b / 0.533 deg below it: %.2f km for 9 m, %.2f km for 10 m, %.2f km for 11 m.  With the cryogenic radiators (5.4 m, Chapter 14) beside a 10 m casing: %.2f km."
  % (9/SUN/1e3, 10/SUN/1e3, 11/SUN/1e3, 15.4/SUN/1e3))
P("      The ring-mounted panels: 55, 83, 110 m give %.1f, %.1f, %.1f km (Section 12.7's 6 to 12 km) ; with the 10 m casing between two 41.5 m wings, %.1f km." % (55/SUN/1e3, 83/SUN/1e3, 110/SUN/1e3, 93/SUN/1e3))
OUT["L"] = dict(casing_m="9 to 11 (Ch3, not finalized) ; 10 (Ch14)", umbra_km={"9": 9/SUN/1e3, "10": 10/SUN/1e3, "11": 11/SUN/1e3, "83": 83/SUN/1e3, "93": 93/SUN/1e3})
N("l", "Casing width in the book", "9 to 11 m on a side (Chapter 3, not final) ; 10 m in Chapter 14", "n/a", "Ch12 Section 12.7 (the passage on doing without the ring-mounted panels)", src="Vol_III_Ch03.md line 238, Vol_III_Ch14.md line 100 (read, not computed)")
N("l", "Full shadow of the casing alone, 9 to 11 m wide ; 10 m", "1.0 to 1.2 km ; 1.1 km", "%.3f to %.3f km ; %.3f km" % (9/SUN/1e3, 11/SUN/1e3, 10/SUN/1e3), "Ch12 Section 12.7")
att = 0.1*95.0e3*30.0
P("\n  M.  Chapter 16 line 142.  The rule of Section 11.2: no attachment adds more than a tenth of the casing's weight over 30 m, %.0f kN.  A trolley line's 2.26 MN gives 8 ; the arm's 23 MN gave 81." % (att/1e3))
P("      The frame has no attachment of its own on the casing: its five lines go to the force spreader cable.  By the same rule the %.2f MN it brings would be %.0f attachments' worth" % (top_frame/1e6, top_frame/att))
P("      of hangers (%.0f for the 6,000 kg/m casing), spread by the force spreader cable.  A guide line is a 100 t line: 8 (16)." % (top_frame/(att/2)))
P("      One thing the rule shows up: a hanger every 200 m carrying %.2f MN is %.1f times the rule's %.0f kN.  The book's own 0.9 MN (Section 12.4) already was.  The hangers need to be about %.0f m apart, or their fittings longer than 30 m."
  % (hang200/1e6, hang200/att, att/1e3, att/w_fsc))
OUT["M"] = dict(rule_kN=att/1e3, frame_equiv=top_frame/att, frame_equiv_6000=top_frame/(att/2), hanger_spacing_m=att/w_fsc)
N("m", "Ch16 line 142, in place of 'per spreader arm': attachments for each guide line ; a frame", "8 (16 for the lighter casings) ; none of its own, its lines go to the force spreader cable (%.0f hangers' worth by the rule of a tenth ; %.0f)" % (top_frame/att, top_frame/(att/2)),
  "%.2f ; %.2f" % (top_frame/att, top_frame/(att/2)), "Ch16 line 142")
N("m", "Hanger spacing that keeps a hanger within a tenth of the casing's weight over 30 m", "%.0f m" % (att/w_fsc), att/w_fsc, "Ch12 line 81 (the book has 200 m and 0.9 MN)")

# =====================================================================================================================
head("ITEM 2.  Paul's deployment sequence of the array (Section 12.11)")
# =====================================================================================================================
RUN = os.path.join(ROOT, "runs", "round9_design", "ch9_optimal", "L2_N165_P12MW")
rows_ = list(csv.DictReader(open(os.path.join(RUN, "rising_radius_daily.csv"))))
day = np.array([float(r_["day"]) for r_ in rows_]); rr = np.array([float(r_["r_m"]) for r_ in rows_]); vc = np.array([float(r_["v_casing_m_s"]) for r_ in rows_])
txt = open(os.path.join(RUN, "rising_radius_lim_simulation.txt"), encoding="utf-8").read()
mon = [(float(m.group(1)), float(m.group(2)), float(m.group(3))) for m in re.finditer(r"^(\d+) mth\s+(-?[\d.]+) m/s.*?([\d.]+) MW\s*$", txt, re.M)]
fin = re.search(r"^final\s+(-?[\d.]+) m/s.*?([\d.]+) MW\s*$", txt, re.M)
md = [float(np.interp(-v, -vc, day)) for _, v, _ in mon] + [day[-1]]; mp = [p for _, _, p in mon] + [float(fin.group(2))]
DEND = float(day[-1]); R_END = float(rr[-1])
def cas(d):
    """Casing on day d: radius, angular rate (west positive, as the run's speed), rate of rise (m/s), rate of slowing (m/s2)."""
    r_c = float(np.interp(d, day, rr)); v = float(np.interp(d, day, vc)); e = 0.5
    d1, d2 = max(d - e, 0.0), min(d + e, DEND)
    rise = (float(np.interp(d2, day, rr)) - float(np.interp(d1, day, rr)))/((d2 - d1)*86400.0)
    dec = -(float(np.interp(d2, day, vc)) - float(np.interp(d1, day, vc)))/((d2 - d1)*86400.0)
    return r_c, v/r_c, rise, dec
def gz(d, z):
    """Weight (m/s2) of a body that turns with the casing, z below it, on day d."""
    r_c, om, _, _ = cas(d); r = r_c - z; return MU/r**2 - om*om*r
def day_at(h_km): return float(np.interp(RE + h_km*1e3, rr, day))
d600 = day_at(600.0); d601 = day_at(601.0); d602 = day_at(602.0)
takeup = Zt - taut["h"]; turn0 = Zt - start["h"]
d_taut = float(np.interp(R_END - takeup, rr, day)); d_turn = float(np.interp(R_END - turn0, rr, day))
P("  Run: runs/round9_design/ch9_optimal/L2_N165_P12MW ; %.0f days ; casing from %.2f to %.2f km." % (DEND, (rr[0] - RE)/1e3, (R_END - RE)/1e3))
P("  (a) The casing passes 600 km on day %.1f, 601 km on day %.1f, 602 km on day %.1f." % (d600, d601, d602))
P("      Phase 1, folded close under the casing: day 0 to day %.0f (ridge 1 km below) or %.0f (2 km)." % (d601, d602))
P("      Phase 2, folded and held at 600 km while the spools pay out: day %.0f to day %.1f, %.0f days." % (d601, d_taut, d_taut - d601))
P("      Phase 3, the suspension lines pull: from day %.1f, when the casing is %.2f km short of 700 km.  The legs begin to turn on day %.2f, %.0f m short, %.1f hours before the last day ends."
  % (d_taut, takeup/1e3, d_turn, turn0, (DEND - d_turn)*24))
ris = [cas(d)[2] for d in np.arange(math.ceil(d601), DEND, 1.0)]
r_last = cas(DEND - 0.3)[2]
P("      Rate of rise of the casing: %.2f to %.2f mm/s through phase 2 (mean %.2f mm/s, %.0f m a day) ; %.2f mm/s, %.0f m a day, over the last days.  391 m of rise takes %.1f hours ; the %.2f km of take-up, %.1f days."
  % (1e3*min(ris), 1e3*max(ris), 1e3*float(np.mean(ris)), 86400*float(np.mean(ris)), 1e3*r_last, 86400*r_last, 391/r_last/3600, takeup/1e3, takeup/r_last/86400))
P("      When the casing reaches 700 km with the hinge held at its working height, the legs are at %.0f deg.  The spool pays out %.0f m more to bring them flat (about %.0f hours at the same rate, or as fast as the spool allows),"
  % (state(Zt, A, L0)["phi"], flat["h"] - Zt, (flat["h"] - Zt)/r_last/3600))
P("      then the sheet is drawn out, then the spool takes the same %.0f m up again.  The frame does not depend on the casing's exact final height." % (flat["h"] - Zt))
OUT["item2"] = dict(a=dict(d600=d600, d601=d601, d602=d602, d_taut=d_taut, d_turn=d_turn, rise_mean_mm_s=1e3*float(np.mean(ris)), rise_min=1e3*min(ris), rise_max=1e3*max(ris), rise_last_mm_s=1e3*r_last,
                           h391=391/r_last/3600, takeup_days=takeup/r_last/86400, phi_at_700=state(Zt, A, L0)["phi"], more_m=flat["h"] - Zt))
# (b) phase 1
m_sheet_km = 0.5*W5*1e3*La/Lr                       # kg of panel for each km of ring
lam_o = RHO*A[1]; lam_i = RHO*A[0]
P("\n  (b) Phase 1.  Weight at the folded array, 1 km and 2 km below the casing (the gradient pull and the casing's own weight as it slows):")
P("    day  | casing  | weight 1 km below | 2 km below | of the final | line ahead of vertical (1 km) | center line (both legs) | a ridge support line (1 km of sheet) | bottom of the slack bights | air there, at the array")
tab_b = []
import pymsis, warnings
warnings.filterwarnings("ignore")
ALT = np.arange(400.0, 721.0, 10.0)
def msis(Fx):
    acc = []
    for hr in range(0, 24, 3):
        o = pymsis.calculate(np.array([np.datetime64("2020-03-21T00:00") + np.timedelta64(hr, "h")]), 0., 0., ALT, Fx, Fx, [[15]*7]); acc.append(np.squeeze(o))
    return np.mean(acc, axis=0)
phase_ = (np.arange(22) + 0.5)/22*np.pi
CYC = np.mean([msis(70 + 130*math.sin(p)**2) for p in phase_], axis=0)
KB = 1.380649e-23
def air(h_km):
    tab = CYC
    rho = math.exp(np.interp(h_km, ALT, np.log(tab[:, 0]))); nO = math.exp(np.interp(h_km, ALT, np.log(tab[:, 3]))); T = float(np.interp(h_km, ALT, tab[:, 10]))
    ntot = sum(math.exp(np.interp(h_km, ALT, np.log(np.maximum(tab[:, j], 1e-30)))) for j in (1, 2, 3, 4, 5, 6, 7)); mbar = rho/ntot
    return rho, nO, T, math.sqrt(8*KB*T/(math.pi*mbar)), math.sqrt(8*KB*T/(math.pi*16*1.66054e-27))
for d in (0, 1, 7, 30, 60, 90, 120, 150, d601):
    r_c, om, rise, dec = cas(d); g1 = gz(d, 1e3); g2 = gz(d, 2e3); hc = (r_c - RE)/1e3
    bot = (L0[1] + 1e3 + xo)/2; zz_ = np.linspace(1e3 + xo, bot, 400); wb = lam_o*float(TRAPZ([gz(d, z) for z in zz_], zz_))
    lean = math.degrees(math.atan2(dec, g1))
    row = dict(day=d, casing_km=hc, g1=g1, g2=g2, lean_deg=lean, Tc_MN=Mfr*g1/1e6, Tridge_MN=m_sheet_km*g1/1e6, bight_bottom_km=hc - bot/1e3, bight_on_leg_kN=wb/1e3, rho_bight=air(hc - bot/1e3)[0], rho_array=air(hc - 1.0)[0])
    tab_b.append(row)
    P("    %5.1f | %6.1f km | %8.5f m/s2    | %8.5f   | 1/%-8.1f   | %6.2f deg                    | %9.4f MN            | %9.4f MN                         | %5.1f km, %6.2f kN on a leg | %.1e, %.1e kg/m3"
      % (d, hc, g1, g2, ga/g1, lean, row["Tc_MN"], row["Tridge_MN"], row["bight_bottom_km"], row["bight_on_leg_kN"], row["rho_bight"], row["rho_array"]))
b0 = tab_b[0]; b1 = tab_b[1]
P("      Day 0: at 1 km the pull is %.2e m/s2, 1/%.0f of the final weight.  The frame's %.0f t pull %.1f kN on the center line ; a kilometer of sheet (%.0f t) pulls %.1f kN on its ridge line." % (b0["g1"], ga/b0["g1"], Mfr/1e3, b0["Tc_MN"]*1e3, m_sheet_km/1e3, b0["Tridge_MN"]*1e3))
P("      The folded array hangs %.1f deg ahead of the vertical on the first day (Section 11.1 has 2.7 deg at 1 km), %.2f deg after a week and %.2f deg after a month.  It needs no tending:" % (b0["lean_deg"], tab_b[2]["lean_deg"], tab_b[3]["lean_deg"]))
P("      the pull is small but never zero or upward, and it grows every day.  Parts closer than about 100 m to the casing lean 25 deg or more on the first day, as the book says of the lines.")
P("      Hinge, folded: each pin carries its leg, %.1f kN on day 0, %.2f MN when phase 1 ends (the leg's full weight is %.2f MN)." % (Mfr/2*b0["g1"]/1e3, Mfr/2*tab_b[-1]["g1"]/1e6, Wl/1e6))
e_c = 0.2887*t["d"]; tilt = math.degrees(math.atan(e_c/(leg/2)))
P("      The hanging legs do not leave a gap by themselves.  With the pins on the lower chords the two legs hang back to back, and each leg's center of mass is %.1f m outside its pin line" % e_c)
P("      (a third of the section's height [taking the pin line through two lower chords: ASSUMED, the joint is not designed]).  A leg hanging free would swing its tip %.1f deg INWARD, %.0f m," % (tilt, 2*e_c))
P("      so the tips press on each other with about %.0f kN at full weight.  That is the wrong way for the opening.  Two cures, both simple: fasten the suspension lines on the upper-chord side," % (Wl*e_c/leg/1e3))
P("      which is the outer face of the folded leg, %.1f m outside the pin line: the lines' pull then turns the tip outward once the two together carry more than %.2f of the leg's weight (they rise to %.2f) ;"
  % (0.866*t["d"], e_c/(0.866*t["d"]), (flat["To"] + flat["Ti"])/Wl))
P("      or hold the pins %.0f m apart on a short center piece, which gives Paul's gap.  [Hand statics on the section of the script's frame.]" % (4*e_c))
P("      The slack suspension lines hang in bights whose lowest point is %.1f km below the casing: %.0f km up on day 0, against %.0f km for the array.  The air there is %.1f times as dense."
  % ((L0[1] + 1e3 + xo)/2e3, b0["bight_bottom_km"], b0["casing_km"] - 1, b0["rho_bight"]/b0["rho_array"]))
P("      A bight pulls on its leg with %.2f kN on day 0 and %.1f kN when phase 1 ends ; the outer line weighs %.1f t.  The lines are CNT and are sheathed for the air in any case (Section 11.1)." % (b0["bight_on_leg_kN"], tab_b[-1]["bight_on_leg_kN"], LN[2]["mass_t"]))
OUT["item2"]["b"] = dict(table=tab_b, e_c_m=e_c, tilt_in_deg=tilt, pins_apart_m=4*e_c, tip_press_kN=Wl*e_c/leg/1e3, attach_ratio=e_c/(0.866*t["d"]), bight_depth_km=(L0[1] + 1e3 + xo)/2e3)
# (c) phase 2
A_sup = LOG9["support"]["0.5"]["T_top_kN_m"]*1e3/S75             # m2 of ridge support line for each meter of ring
m_arr = 0.5*W5*La/Lr + Mfr/S10                                   # panels and frames for each meter of ring
def spool_pull(d, ell, load, lamA):
    """Pull at a spool on the casing: the load 'ell' below it, and the line paid out so far (RHO x lamA for each meter)."""
    zz_ = np.linspace(0, ell, 60); return load*gz(d, ell) + RHO*lamA*float(TRAPZ([gz(d, z) for z in zz_], zz_))
dd = np.arange(d601, d_taut, 0.5); E2 = 0.0; prev = None; pk = 0.0; tcm = 0.0; trm = 0.0
for d in dd:
    r_c = cas(d)[0]; ell = r_c - (RE + 600e3); T_m = spool_pull(d, ell, m_arr, A_sup)
    tcm = max(tcm, spool_pull(d, ell, Mfr, A_c)); trm = max(trm, spool_pull(d, ell, m_sheet_km, A_sup*1e3))
    if prev is not None: E2 += 0.5*(T_m + prev[1])*(ell - prev[0])
    prev = (ell, T_m)
E2_ring = E2*Lr; P2_mean = E2_ring/((d_taut - d601)*86400)
r0 = cas(0)[0]; ell_s = np.linspace(0, 100e3, 400); T0s = [spool_pull(0.0, e, m_arr, A_sup) for e in ell_s]; E0 = float(TRAPZ(T0s, ell_s))*Lr
Tc0 = spool_pull(0.0, 100e3, Mfr, A_c); Tr0 = spool_pull(0.0, 100e3, m_sheet_km, A_sup*1e3)
P("\n  (c) Phase 2.  The spools pay out at the casing's rate of rise, %.2f mm/s on the mean, %.1f km in %.0f days.  The array is at 600 km and has nearly its full weight by then (%.2f m/s2 on day %.0f, %.2f at the end)."
  % (1e3*float(np.mean(ris)), (cas(d_taut)[0] - (RE + 601e3))/1e3 + 1.0, d_taut - d601, gz(d601, 1e3), d601, gz(d_taut, cas(d_taut)[0] - RE - 600e3)))
P("      Center array support line (both legs and its own weight): up to %.2f MN at the spool ; it is sized for %.2f MN at the top, against %.2f MN while the frame opens and %.2f MN at work." % (tcm/1e6, LN[0]["top_MN"], LN[0]["opening_MN"], LN[0]["open_MN"]))
P("      A ridge support line (1 km of sheet): up to %.2f MN at the spool (%.2f MN at work)." % (trm/1e6, T_ridge/1e6))
P("      Brakes around the whole ring: %.3f EJ over phase 2, %.2f GW on the mean, %.0f W for each meter of ring.  For scale: the deployment itself takes 19.2 EJ." % (E2_ring/1e18, P2_mean/1e9, P2_mean/Lr))
P("      Against Section 12.11 as it stands, the whole 100 km let down on the first day: %.4f EJ, %.0f times less ; center line %.2f MN and ridge line %.2f MN at the spool at the end of the lowering."
  % (E0/1e18, E2_ring/E0, Tc0/1e6, Tr0/1e6))
P("      So Paul's sequence costs the brakes %.3f EJ more, heat that the spools on the casing have to shed at %.0f W for each meter of ring, and spools that work at full line tension." % ((E2_ring - E0)/1e18, P2_mean/Lr))
P("      What it buys is in (d).  [The spools are on the casing here.  The book's belay spool rides down its own line ; at 100 km of 18 mm line that is a different machine, not looked at.]")
OUT["item2"]["c"] = dict(days=d_taut - d601, Tc_spool_max_MN=tcm/1e6, Tridge_spool_max_MN=trm/1e6, E_ring_EJ=E2_ring/1e18, P_mean_GW=P2_mean/1e9, W_per_m=P2_mean/Lr, E_day0_EJ=E0/1e18, ratio=E2_ring/E0, Tc_day0_MN=Tc0/1e6, Tr_day0_MN=Tr0/1e6)
# (d) the air
def path(depth):
    """Air along a path. depth(h_casing_km) -> altitude of the strip in km.  Faces along the motion (grazing), as round 15."""
    out = []; prev = None; flu = flr = imp = 0.0
    for d in list(np.arange(0.0, DEND, 1.0)) + [DEND]:
        r_c, om, _, _ = cas(d); h_s = depth((r_c - RE)/1e3); r_s = RE + h_s*1e3; v = (om + W_E)*r_s
        rho, nO, T, cb, cO = air(h_s); tau = rho*cb/4*abs(v); drag = 2*W5*tau*2*math.pi*r_s
        cur = dict(day=d, alt=h_s, v=v, rho=rho, drag_MN=drag/1e6, heat=rho*cb/4*v*v/2, f_face=nO*cO/4/1e4, f_ram=nO*abs(v)/1e4, eV=0.5*16*1.66054e-27*v*v/1.602177e-19)
        if prev is not None:
            dt = (d - prev["day"])*86400; flu += 0.5*(cur["f_face"] + prev["f_face"])*dt; flr += 0.5*(cur["f_ram"] + prev["f_ram"])*dt; imp += 0.5*(cur["drag_MN"] + prev["drag_MN"])*1e6*dt
        prev = cur; out.append(cur)
    return dict(rows=out, flu=flu, flr=flr, imp=imp, drag0=out[0]["drag_MN"], drag_max=max(c["drag_MN"] for c in out), heat_max=max(c["heat"] for c in out), alt0=out[0]["alt"],
                days_1eV=max([c["day"] for c in out if c["eV"] > 1.0] + [0]), eV0=out[0]["eV"], flu_above_1eV=sum(0.5*(a["f_face"] + b["f_face"])*(b["day"] - a["day"])*86400 for a, b in zip(out[:-1], out[1:]) if a["eV"] > 1.0))
PATHS = {"Paul's, 1 km under the casing, then 600 km": lambda hc: min(hc - 1.0, 600.0), "Paul's, 2 km under": lambda hc: min(hc - 2.0, 600.0), "100 km under from the first day (round 15, Section 12.11 now)": lambda hc: hc - 100.0}
PA = {k: path(f_) for k, f_ in PATHS.items()}
yr = 365.25*86400; r600 = air(600.0); f600 = r600[1]*r600[4]/4/1e4
P("\n  (d) The air (NRLMSIS 2.1 as round 15 ran it: equator, 21 March, eight local times, mean of a solar cycle).  The sheets lie along the motion, folded or spread, so each face is grazed.")
P("    path                                                          | first day: up, through the air | drag, first day and largest | heat on a face, largest | O atoms on a face over the deployment | of a year at rest at 600 km | on a forward face | struck above 1 eV for | at more than 1 eV, on a face")
for k, p in PA.items():
    P("    %-61s | %5.1f km, %5.0f m/s, %.2f eV    | %6.3f, %6.3f MN            | %.4f W/m2             | %.2e /cm2                        | %5.2f years                  | %.2e /cm2      | %3.0f days              | %.2e /cm2"
      % (k, p["alt0"], p["rows"][0]["v"], p["eV0"], p["drag0"], p["drag_max"], p["heat_max"], p["flu"], p["flu"]/(f600*yr), p["flr"], p["days_1eV"], p["flu_above_1eV"]))
P("    a year at rest at 600 km: %.2e O atoms/cm2 on a face, at about 0.1 eV." % (f600*yr))
pk1 = PA["Paul's, 1 km under the casing, then 600 km"]; pk15 = PA["100 km under from the first day (round 15, Section 12.11 now)"]
P("    The energy of the atoms by day on Paul's path: " + " ; ".join("day %d %.2f eV" % (c["day"], c["eV"]) for c in pk1["rows"] if c["day"] in (0, 30, 60, 90, 120, 150, 180, 210, 240)))
P("    Paul's path takes %.1f times fewer oxygen atoms than the path 100 km down, and %.1f times fewer of the fast ones.  The speed through the air is the same on both paths, within 1.5 %%: it is the casing's." % (pk15["flu"]/pk1["flu"], pk15["flu_above_1eV"]/pk1["flu_above_1eV"]))
P("    What changes is the density: %.1e kg/m3 at %.0f km against %.1e at %.0f km on the first day, a factor of %.1f." % (pk15["rows"][0]["rho"], pk15["alt0"], pk1["rows"][0]["rho"], pk1["alt0"], pk15["rows"][0]["rho"]/pk1["rows"][0]["rho"]))
P("    Folded, the inner faces of the two half-sheets shield each other in part ; that was not counted.  What a grazing 5 eV atom does to a coated film is still not known to me.")
OUT["item2"]["d"] = {k: {kk: vv for kk, vv in p.items() if kk != "rows"} for k, p in PA.items()}; OUT["item2"]["d"]["rest600_per_year"] = f600*yr
OUT["item2"]["d"]["eV_by_day"] = {int(c["day"]): c["eV"] for c in pk1["rows"] if c["day"] in (0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)}
# (e) phase 3
P_TETHER = 1.0e6
P("\n  (e) Phase 3.  When the suspension lines come taut: with the legs hanging DOWN from the hinge, a line point is below the hinge, so a suspension line pulls when")
P("      the center line's length reaches the suspension line's length MINUS the distance from the hinge to its point on the leg (%.0f m for the outer line), not plus.  The brief's 'plus' would be a frame" % xo)
P("      folded with its legs pointing up.  And the lengths are the stretched ones: slack lengths %.3f km (outer), %.3f km (inner) against %.3f km for the center line at work." % (L0[1]/1e3, L0[0]/1e3, Zt/1e3))
P("      Because the lines stretch %.1f km under their working pull, both suspension lines begin to pull %.2f km before the end and take the weight over gradually (section B) ; the inner one is not slack" % (LN[2]["stretch_open_km"], takeup/1e3))
P("      during the opening, as it would be with stiff lines.  The frame is at its full weight (%.2f m/s2 at the array): section B is this case." % ga)
P("      Round 15's opening numbers were for one stiff line at 4/5 and a 598 t frame.  Its pulls (1.52 and 1.83 MN) are replaced by section B's.  Its drawing-out numbers are replaced too (section B).")
P("      Time: the take-up is %.1f days of the casing's rise ; the legs turn in the last %.1f hours and the %.0f m after ; the spool can do the last part faster." % (takeup/r_last/86400, (DEND - d_turn)*24, flat["h"] - Zt))
P("      Drawing out: %.2f GJ for each half-sheet for each 10 km.  At a tether winch power of %.0f MW for each half-sheet for each 10 km [my choice: 100 W for each meter of edge] it takes %.0f minutes," % (W_teth/1e9, P_TETHER/1e6, W_teth/P_TETHER/60))
P("      and %.1f GW around the ring if every span is drawn out at once.  The ridge spools' brakes take %.2f GJ for each half-sheet for each 10 km in the same time, %.1f MW." % (2*n_fr*P_TETHER/1e9, E_ridge/1e9, E_ridge/(W_teth/P_TETHER)/1e6))
P("      Driven by the casing's rise or by the brake: the same forces.  Seen from the ground the hinge stays where it is and the suspension lines lift the line points as the casing rises.")
P("      The casing's rise does the work: it lifts the frame's center of mass by half a leg (%.2f GJ), and it pulls the center line off its spool against the brake (%.2f GJ while the legs turn)." % (2*Wl*leg/2/1e9, E_open/1e9))
P("      The brake's heat is the same as if the hinge were let down from a casing at rest.  For all %d frames the two together are %.1e of the deployment's 19.2 EJ." % (round(n_fr), n_fr*(2*Wl*leg/2 + E_open)/19.2e18))
OUT["item2"]["e"] = dict(takeup_km=takeup/1e3, turn_hours=(DEND - d_turn)*24, P_tether_MW=P_TETHER/1e6, draw_min=W_teth/P_TETHER/60, ring_GW=2*n_fr*P_TETHER/1e9, ridge_brake_MW=E_ridge/(W_teth/P_TETHER)/1e6, lift_GJ=2*Wl*leg/2/1e9)
# (f) the possibility
P("\n  (f) The possibility: no panels on the ring, the strip spread from the start as the deployment's power.")
g1d0 = gz(0.0, 1e3); g2d0 = gz(0.0, 2e3)
P("      Opened 1 km under the casing on day 0, at %.2e m/s2 (1/%.0f of the final weight): frame %.1f kN, push on a frame %.1f kN, tether work %.2f MJ for each half-sheet for each 10 km (2.94 GJ at full weight)."
  % (g1d0, ga/g1d0, Mfr*g1d0/1e3, F*g1d0/ga/1e3, W_teth*g1d0/ga/1e6))
P("      Even drawn out with the ridge held, the tip load would be %.1f kN in place of %.2f MN.  At 2 km: 1/%.0f.  Every pull then grows together as the casing slows, in the same proportion," % (ptm*1e6*g1d0/ga/1e3, ptm, ga/g2d0))
P("      so the lines stretch together and the shape holds.  Then the five lines of each frame and the ridge lines are spooled out equally (the suspension lines cannot be of fixed length in this version).")
def supply(h_km): return TW_PER_KM_600*(W5/1e3)*1e3*(RE + h_km*1e3)/(RE + 600e3)*lit(h_km*1e3)/lit(600e3)
def sched(depth):
    rat = []; sup = []
    for d in np.arange(0.0, DEND + 0.1, 1.0):
        hc = (cas(min(d, DEND))[0] - RE)/1e3; s_ = supply(depth(hc)); dem = float(np.interp(d, md, mp))*N_SITES/1e3; rat.append(s_/dem); sup.append(s_)
    return min(rat), min(sup), max(sup)
dem_pk = max(mp)*N_SITES/1e3
P("    schedule                                         | strip gives      | least, against the demand (peak %.0f GW) | drag first day | O atoms on a face  | years at rest at 600 km | above 1 eV" % dem_pk)
OUT["item2"]["f"] = dict(g_1km_day0=g1d0, frac=g1d0/ga, schedules={})
for nm, key in (("A. 100 km below from the start (427 km up)", "100 km under from the first day (round 15, Section 12.11 now)"), ("B. 1 km below until 600 km, then held at 600 km", "Paul's, 1 km under the casing, then 600 km")):
    rmin, smin, smax = sched(PATHS[key]); p = PA[key]
    P("    %-48s | %5.0f to %5.0f GW | %5.1f times                             | %.3f MN       | %.2e /cm2      | %5.2f                   | %3.0f days" % (nm, smin, smax, rmin, p["drag0"], p["flu"], p["flu"]/(f600*yr), p["days_1eV"]))
    OUT["item2"]["f"]["schedules"][nm] = dict(ratio_min=rmin, supply_min_GW=smin, supply_max_GW=smax, drag0_MN=p["drag0"], AO=p["flu"], AO_years=p["flu"]/(f600*yr), days_1eV=p["days_1eV"])
foot = 2*b_half
P("      The casing's own shadow on a strip spread close under it, casing %.0f m wide [the book: 9 to 11 m, Chapter 14 takes 10 m]:" % W_CAS)
for z in (1e3, 2e3):
    um = max(0.0, W_CAS - z*SUN); pen = W_CAS + z*SUN; x = math.degrees(math.atan(b_half/z))/23.44; dd_ = 365.25*(2/math.pi)*math.asin(min(1.0, x))
    P("        %.0f km under: full shadow %.1f m wide, part shadow %.0f m wide ; it takes at most the casing's width out of the light, %.1f %% of the %.0f m footprint ; it falls on the strip %.0f days a year."
      % (z/1e3, um, pen, 100*W_CAS/foot, foot, dd_))
    OUT["item2"]["f"]["shadow_%gkm" % (z/1e3)] = dict(umbra_m=um, penumbra_m=pen, loss=W_CAS/foot, days=dd_)
P("      With no panels on the ring the shadow costs about 1 %%, so it does not ask for 100 km.  With the ring-mounted panels there (83 m and the casing), a strip 1 or 2 km under them would lose up to %.1f %%," % (100*93/foot))
P("      inside their full shadow, which ends %.0f km down.  That is the reason for a long hang when both sets of panels are carried." % (93/SUN/1e3))
# (g)
P("\n  (g) A least hang distance in phase 1.  What the book gives: the force spreader cable hangs 'a few hundred meters below the casing' (Section 12.4) on hangers 'a few hundred meters apart' ;")
P("      Section 11.6 uses hangers 500 m long as its example.  The ring-mounted panels reach 41 m to each side of the casing (Chapter 8) and the cryogenic radiators are 5.4 m wide along its outside (Chapter 14) ;")
P("      neither is below the casing.  Section 12.7's 20 km is for the shadow on a working array, not for a folded one.  So the book sets nothing but the force spreader cable: the folded array has to hang below it,")
P("      and the 1 km of the brief does.  The folded frame is %.0f m from hinge to tip and the folded sheet %.0f m, so the bundle is about %.1f km tall whatever the distance.  No clearance was invented." % (leg, s_half, (xo - drop + s_half)/1e3))
GD = "Ch12 Section 12.11 (lines 257 to 263) ; Ch11 Section 11.1 (lines 100 to 106)"
N("2a", "Day the casing passes 600 km ; 601 km", "day %.0f ; day %.0f" % (d600, d601), "%.2f ; %.2f" % (d600, d601), GD, key="ITEM 2 (a)")
N("2a", "Phase 2 (held at 600 km, spools paying out): from, to, days", "day %.0f to day %.0f, %.0f days" % (d601, d_taut, d_taut - d601), "%.2f ; %.2f" % (d601, d_taut), GD, key="ITEM 2 (a)")
N("2a", "Phase 3 begins (suspension lines first pull): day ; casing short of 700 km", "day %.0f ; %.1f km" % (d_taut, takeup/1e3), "%.2f ; %.3f km" % (d_taut, takeup/1e3), GD, key="ITEM 2 (a)")
N("2a", "Legs begin to turn: casing short of 700 km ; hours before the end", "%.0f m ; %.0f hours" % (turn0, (DEND - d_turn)*24), "%.1f m ; %.2f h" % (turn0, (DEND - d_turn)*24), GD, key="ITEM 2 (a)")
N("2a", "Casing's rate of rise: mean of phase 2 ; over the last days", "%.1f mm/s (%.0f m a day) ; %.1f mm/s (%.0f m a day)" % (1e3*float(np.mean(ris)), 86400*float(np.mean(ris)), 1e3*r_last, 86400*r_last), "%.4f ; %.4f mm/s" % (1e3*float(np.mean(ris)), 1e3*r_last), GD, key="ITEM 2 (a)")
N("2b", "Weight 1 km below the casing on day 0 ; share of the final weight", "%.4f m/s2 ; 1/%.0f" % (b0["g1"], ga/b0["g1"]), "%.6f ; %.1f" % (b0["g1"], ga/b0["g1"]), GD, key="ITEM 2 (b)")
N("2b", "Weight 1 km below: day 7 ; day 30 ; day 90 ; when phase 1 ends", "%.2f ; %.2f ; %.2f ; %.2f m/s2" % (tab_b[2]["g1"], tab_b[3]["g1"], tab_b[5]["g1"], tab_b[-1]["g1"]), "%.4f ; %.4f ; %.4f ; %.4f" % (tab_b[2]["g1"], tab_b[3]["g1"], tab_b[5]["g1"], tab_b[-1]["g1"]), GD, key="ITEM 2 (b)")
N("2b", "Day 0, folded 1 km below: center line (frame) ; a ridge line (1 km of sheet) ; lean ahead of the vertical", "%.1f kN ; %.1f kN ; %.1f deg" % (b0["Tc_MN"]*1e3, b0["Tridge_MN"]*1e3, b0["lean_deg"]), "%.3f kN ; %.3f kN ; %.3f deg" % (b0["Tc_MN"]*1e3, b0["Tridge_MN"]*1e3, b0["lean_deg"]), GD, key="ITEM 2 (b)")
N("2b", "Lowest point of the slack suspension lines: below the casing ; altitude on day 0", "%.0f km ; %.0f km" % ((L0[1] + 1e3 + xo)/2e3, b0["bight_bottom_km"]), "%.2f ; %.2f km" % ((L0[1] + 1e3 + xo)/2e3, b0["bight_bottom_km"]), GD, key="ITEM 2 (b)")
N("2c", "Phase 2 brakes, whole ring: energy ; mean power ; for each meter of ring", "%.3f EJ ; %.1f GW ; %.0f W/m" % (E2_ring/1e18, P2_mean/1e9, P2_mean/Lr), "%.5f EJ ; %.4f GW ; %.2f W/m" % (E2_ring/1e18, P2_mean/1e9, P2_mean/Lr), GD, key="ITEM 2 (c)")
N("2c", "The same lowering done on the first day (Section 12.11 now): energy ; ratio", "%.4f EJ ; %.0f times less" % (E0/1e18, E2_ring/E0), "%.6f EJ ; %.2f" % (E0/1e18, E2_ring/E0), GD, key="ITEM 2 (c)")
N("2c", "Pull at the spool in phase 2: center line ; a ridge line", "%.1f MN ; %.1f MN" % (tcm/1e6, trm/1e6), "%.4f ; %.4f MN" % (tcm/1e6, trm/1e6), GD, key="ITEM 2 (c)")
for k, p in PA.items():
    N("2d", "Air, path '%s': altitude and speed through the air on day 0 ; drag on day 0 ; O atoms on a face over the deployment ; in years at rest at 600 km ; days above 1 eV" % k,
      "%.0f km, %.1f km/s ; %.2f MN ; %.1e /cm2 ; %.1f years ; %.0f days" % (p["alt0"], p["rows"][0]["v"]/1e3, p["drag0"], p["flu"], p["flu"]/(f600*yr), p["days_1eV"]),
      "%.2f km, %.1f m/s ; %.4f MN ; %.4e ; %.3f ; %.0f" % (p["alt0"], p["rows"][0]["v"], p["drag0"], p["flu"], p["flu"]/(f600*yr), p["days_1eV"]), GD, key="ITEM 2 (d)")
N("2d", "O atoms on a face in a year at rest at 600 km", "%.1e /cm2" % (f600*yr), f600*yr, "Ch12 Section 12.7, 12.11", key="ITEM 2 (d)")
N("2e", "Drawing out at 1 MW for each half-sheet for each 10 km [my choice]: time ; around the ring", "%.0f minutes ; %.1f GW" % (W_teth/P_TETHER/60, 2*n_fr*P_TETHER/1e9), "%.2f min ; %.3f GW" % (W_teth/P_TETHER/60, 2*n_fr*P_TETHER/1e9), GD, key="ITEM 2 (e)")
for nm, o in OUT["item2"]["f"]["schedules"].items():
    N("2f", "Possibility, schedule %s: strip's power ; least multiple of the demand ; O atoms in years at rest at 600 km" % nm, "%.1f to %.1f TW ; %.0f times ; %.1f years" % (o["supply_min_GW"]/1e3, o["supply_max_GW"]/1e3, o["ratio_min"], o["AO_years"]),
      "%.1f to %.1f GW ; %.2f ; %.3f" % (o["supply_min_GW"], o["supply_max_GW"], o["ratio_min"], o["AO_years"]), "Ch12 Section 12.11 (the possibility passage)", key="ITEM 2 (f)")
N("2f", "Casing's shadow on a strip 1 km ; 2 km under it: largest loss ; days a year", "%.1f %% ; %.0f and %.0f days" % (100*W_CAS/foot, OUT["item2"]["f"]["shadow_1km"]["days"], OUT["item2"]["f"]["shadow_2km"]["days"]), "%.4f" % (W_CAS/foot), "Ch12 Section 12.7 / 12.11", key="ITEM 2 (f)")

# ---------------------------------------------------------------------------------------------------------------------
with open(os.path.join(LOG_DIR, "round16_number_sheet.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.ndarray): return [clean(v) for v in o.tolist()]
    return o
with open(os.path.join(LOG_DIR, "round16_number_sheet.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
with open(os.path.join(LOG_DIR, "round16_numbers.md"), "w", encoding="utf-8") as fh:
    fh.write("| No. | Item | Number | Print in the book | Exact value | Script and log section | Goes to |\n|---|---|---|---|---|---|---|\n")
    for i, r_ in enumerate(SHEET, 1):
        fh.write("| %d | %s | %s | %s | %s | %s | %s |\n" % (i, r_["sec"], r_["what"].replace("|", "/"), r_["show"], fx(r_["exact"]), r_["src"], r_["goes"]))
print("\nlog written: runs/round16_logs/round16_number_sheet.txt and .json ; sheet rows: %d ; %.0f s" % (len(SHEET), time.time() - T0))
