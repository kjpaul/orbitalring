#!/usr/bin/env python3
"""
Round 17: Table 11.1 of Volume III at the width of the design case (1,111.1 m, not the round 1,110 m), and the folded state of the
truss frame at the moment its legs begin to turn, for the opening figure (Figure 12.4).
Brief: Vol_III/Claude_Code_instructions_round17_2026-11-18.md, items 5 and 1.

No new model.  Sections A and B of tools/round16_number_sheet.py are run here as they stand, with their printing muted, and the
function frame_case of its section G is taken from the same file, so that the code is the same code.

Run from Python/orbitalring:  python tools/round17_table_and_fold.py
Log: runs/round17_logs/round17_table_and_fold.txt and .json
"""
import ast, contextlib, io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
R16 = os.path.join(HERE, "round16_number_sheet.py")
SRC = open(R16, encoding="utf-8").read()
MARK = 'head("C.  Yaw'
assert MARK in SRC
ns = {"__name__": "round16_prefix", "__file__": R16}
_argv = sys.argv; sys.argv = [R16]
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(SRC.split(MARK)[0], R16, "exec"), ns)
    for node in ast.parse(SRC).body:
        if isinstance(node, ast.FunctionDef) and node.name == "frame_case":
            exec(compile(ast.Module(body=[node], type_ignores=[]), R16, "exec"), ns)
sys.argv = _argv
frame_case = ns["frame_case"]; W5 = ns["W5"]; S10 = ns["S10"]; Mfr = ns["Mfr"]; F = ns["F"]; t = ns["t"]; leg = ns["leg"]; drop = ns["drop"]
state = ns["state"]; A = ns["A"]; L0 = ns["L0"]; Zt = ns["Zt"]; seq = ns["seq"]; taut = ns["taut"]; start = ns["start"]; flat = ns["flat"]; Wl = ns["Wl"]
xi = ns["xi"]; xo = ns["xo"]

LOG_DIR = os.path.join(ROOT, "runs", "round17_logs"); os.makedirs(LOG_DIR, exist_ok=True)
buf = io.StringIO()
def P(s=""):
    print(s); buf.write(s + "\n")
OUT = {}

P("Round 17, item 5.  Table 11.1 at the width of the design case.")
P("  Sheet of the design case (round 9, d5): %.4f m.  Round 16 section G used 1,110.0 m." % W5)
CASES = (("Design case", 0.5, 45, 1/8), ("Lighter panels (0.2 kg/m2)", 0.2, 45, 1/8), ("Steeper sheet (60 deg at the ridge)", 0.5, 60, 1/8),
         ("Deeper bow (a quarter of the width)", 0.5, 45, 1/4), ("All three together", 0.2, 60, 1/4))
def printed(c):
    return ["%.1f MN" % c["push_MN"], "%.0f t" % c["frame_t"], "%.0f m" % c["depth_m"], "%.0f MN, %.0f mm" % (c["TR_MN"], c["ret_mm"]), "%.0f kg/m" % c["ret_kgm"], "%.0f m" % c["drop_m"], "%.0f m" % c["Lf"]]
COLS = ["push", "frame mass", "depth", "retention line", "retention lines for each meter", "edge below the ridge", "frame length"]
OUT["table_11_1"] = {}; moved = []
P("    case                                | push     | frame  | depth | retention line   | for each meter | edge below ridge | frame length |  exact: push, mass, depth, pull, across, kg/m, drop, length")
for nm, a, th, b in CASES:
    c = frame_case(a, W5, S10, th, b); o = frame_case(a, 1110.0, S10, th, b); pc, po = printed(c), printed(o)
    OUT["table_11_1"][nm] = dict(at_1111=c, at_1110=o, printed_1111=pc, printed_1110=po)
    P("    %-35s | %-8s | %-6s | %-5s | %-16s | %-14s | %-16s | %-12s |  %.3f MN, %.2f t, %.3f m, %.2f MN, %.2f mm, %.3f kg/m, %.2f m, %.2f m"
      % (nm, pc[0], pc[1], pc[2], pc[3], pc[4], pc[5], pc[6], c["push_MN"], c["frame_t"], c["depth_m"], c["TR_MN"], c["ret_mm"], c["ret_kgm"], c["drop_m"], c["Lf"]))
    for k, x, y in zip(COLS, po, pc):
        if x != y: moved.append((nm, k, x, y))
P("  Printed values that move from the 1,110 m table to the 1,111 m table:")
for m in moved: P("    %-35s %-32s %s  becomes  %s" % m)
if not moved: P("    none")
OUT["moved"] = moved
dc = OUT["table_11_1"]["Design case"]["at_1111"]
P("  Check against round 16 section A (frame_system, the design frame): %.3f t, push %.4f MN, depth %.3f m.  Here: %.3f t, %.4f MN, %.3f m.  Difference in mass %.4f t."
  % (Mfr/1e3, F/1e6, t["d"], dc["frame_t"], dc["push_MN"], dc["depth_m"], dc["frame_t"] - Mfr/1e3))
OUT["section_A"] = dict(frame_t=Mfr/1e3, push_MN=F/1e6, depth_m=t["d"])

P("\nRound 17, item 1.  The folded frame at the moment the legs begin to turn (round 16 section B, the lines as springs).")
P("  'Working place' is the hinge's depth with the frame open at work, %.3f km below the top of the lines.  Hinge above (+) or below (-) it:" % (Zt/1e3))
def row(nm, s):
    P("    %-46s hinge %+8.1f m ; legs %5.2f deg from the vertical ; inner line %.3f MN, outer %.3f MN, center line %.3f MN (%.1f %% of the frame's weight)"
      % (nm, Zt - s["h"], s["phi"], s["Ti"]/1e6, s["To"]/1e6, s["Tc"]/1e6, 50*s["Tc"]/Wl))
    return dict(hinge_above_work_m=Zt - s["h"], phi_deg=s["phi"], Ti_MN=s["Ti"]/1e6, To_MN=s["To"]/1e6, Tc_MN=s["Tc"]/1e6)
OUT["fold"] = dict(first_pull=row("suspension lines first pull", taut), begin_turn=row("legs begin to turn", start))
before = [s for s in seq if s["h"] < start["h"]][-1]
OUT["fold"]["just_before_turn"] = row("last step before the legs turn", before)
for ph in (30.0, 45.0, 60.0):
    s = next(s for s in seq if s["phi"] >= ph); OUT["fold"]["phi_%d" % ph] = row("legs at %.0f deg" % ph, s)
OUT["fold"]["flat_no_push"] = row("flat, no push from the sheet yet", flat)
bt = OUT["fold"]["begin_turn"]; fl = OUT["fold"]["flat_no_push"]
P("  From the legs beginning to turn to flat the spool pays out %.0f m: %.0f m down to the working place and %.0f m past it." % (bt["hinge_above_work_m"] - fl["hinge_above_work_m"], bt["hinge_above_work_m"], -fl["hinge_above_work_m"]))
P("  So the folded hinge in the figure is %.0f m above the open frame, which is drawn at its working place.  (The model holds the folded legs at 1 deg ; that is its lower limit, not a result.)" % bt["hinge_above_work_m"])
P("  With the hinge there, the line points on the folded legs are %.0f m and %.0f m below the hinge, so the outer point is %.0f m below the open frame and the tips %.0f m below it."
  % (xi*math.cos(math.radians(bt["phi_deg"])), xo*math.cos(math.radians(bt["phi_deg"])), xo*math.cos(math.radians(bt["phi_deg"])) - bt["hinge_above_work_m"], leg*math.cos(math.radians(bt["phi_deg"])) - bt["hinge_above_work_m"]))
OUT["fold"].update(leg_m=leg, drop_m=drop, xi_m=xi, xo_m=xo, spool_turn_to_flat_m=bt["hinge_above_work_m"] - fl["hinge_above_work_m"])
with open(os.path.join(LOG_DIR, "round17_table_and_fold.txt"), "w", encoding="utf-8") as fh: fh.write(buf.getvalue())
with open(os.path.join(LOG_DIR, "round17_table_and_fold.json"), "w", encoding="utf-8") as fh: json.dump(OUT, fh, indent=1, default=float)
print("ok")
