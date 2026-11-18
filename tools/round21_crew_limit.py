#!/usr/bin/env python3
"""
Round 21: the crew limit of 3 g on the net outward load alone.
Brief: Vol_III/Claude_Code_instructions_round21_2026-11-18.md.

lsm_simulation.py's g-load limit now acts on the net outward load alone and does not throttle the thrust.  This script runs the
crewed launch at the 3 g limit again, with the braking of its sled, and the casing's swing after it, and sets them against round 20.
It also runs the design launch and checks that nothing of it has moved.  The run function and the casing's integrator are round 20's,
taken from tools/round20_brake_felt_pods.py as they stand.

Run from Python/orbitalring:  python tools/round21_crew_limit.py
Log: runs/round21_logs/round21_crew_limit.txt and .json ; the simulation's own printing: runs/round21_logs/lsm_<case>.txt
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
src = open(os.path.join(HERE, "round20_brake_felt_pods.py"), encoding="utf-8").read()
cut = src.index('P("Round 20.'); assert cut > 0
__file__ = os.path.join(HERE, "round20_brake_felt_pods.py")
exec(compile(src[:cut], "round20_brake_felt_pods.py (prefix)", "exec"))          # lsm_run, casing_rows, the hold, P, OUT, cfg, sim, np
LOG_DIR = os.path.join(ROOT, "runs", "round21_logs"); os.makedirs(LOG_DIR, exist_ok=True)
L20 = json.load(open(os.path.join(ROOT, "runs", "round20_logs", "round20_brake_felt_pods.json"), encoding="utf-8"))

P("Round 21.  The crew limit of 3 g on the net outward load alone.")
P("\n0. The design launch, to see that nothing of it moved")
D = lsm_run("design_brake"); d20 = L20["design"]
P("   %.1f minutes ; %.2f laps ; %.2f MN ; %.2f GW ; impulse %.4e N s ; peak load %.3f g ; braking %.1f s, %.4e N s." % (D["dur"]/60, D["x_end"]/cfg.L_RING, D["F_peak"]/1e6, D["P_peak"]/1e9, D["impulse"], D["g_peak"], D["brake"]["duration_s"], D["brake"]["impulse"]))
assert abs(D["dur"] - d20["dur"]) < 1e-6 and abs(D["impulse"] - d20["impulse"]) < 1.0 and abs(D["brake"]["duration_s"] - d20["brake"]["duration_s"]) < 1e-6 and abs(D["brake"]["impulse"] - d20["brake"]["impulse"]) < 1.0 and abs(D["g_peak"] - d20["g_peak"]) < 1e-9
P("   The same as round 20 to the last digit.")
OUT["design"] = strip(D)

P("\n1. The crewed launch at the 3 g limit")
R = cfg.R_ORBIT; v_rot = cfg.V_CASING_SPACE
v3_space = math.sqrt((3*cfg.G_0 + cfg.G_250)*R); v3_cas = v3_space - v_rot
P("   A net outward load of 3 g is reached at %.1f m/s in space, %.1f m/s over the casing.  lsm_simulation.py --v_launch=%.1f --g-limit=3 --brake" % (v3_space, v3_cas, v3_cas))
C = lsm_run("crewed", V_LAUNCH=v3_cas, G_LOAD_MAX=3.0); c20 = L20["crewed"]; B = C["brake"]; b20 = c20["brake"]
v_esc = math.sqrt(2*e1.GM/R); v_inf = math.sqrt(C["v_space"]**2 - v_esc**2)
a_out_end = (C["v_space"]**2/R - cfg.G_250)/cfg.G_0; a_out_last = C["ar"][-1]/cfg.G_0
P("   %-44s | %-18s | %-18s" % ("", "round 20", "round 21"))
def row(what, a, b, fmt):
    P("   %-44s | %-18s | %-18s" % (what, fmt % a, fmt % b))
row("Release in space (m/s)", c20["v_space"], C["v_space"], "%.1f"); row("Release over the casing (m/s)", c20["v_end"], C["v_end"], "%.1f")
row("Duration (minutes)", c20["dur"]/60, C["dur"]/60, "%.2f"); row("Duration (s)", c20["dur"], C["dur"], "%.1f"); row("Distance (km)", c20["x_end"]/1e3, C["x_end"]/1e3, "%.1f"); row("Laps", c20["x_end"]/cfg.L_RING, C["x_end"]/cfg.L_RING, "%.3f")
row("Peak thrust (MN)", c20["F_peak"]/1e6, C["F_peak"]/1e6, "%.2f"); row("Thrust at the release (MN)", c20["F_end"]/1e6, C["F_end"]/1e6, "%.2f"); row("Impulse (N s)", c20["impulse"], C["impulse"], "%.4e")
row("Net outward load at the release (g)", c20["a_out_g"], C["a_out_g"], "%.4f"); row("Peak load, vector sum (g)", c20["g_peak"], C["g_peak"], "%.3f")
row("Braking of the sled (minutes)", b20["duration_s"]/60, B["duration_s"]/60, "%.2f"); row("Braking of the sled (s)", b20["duration_s"], B["duration_s"], "%.1f"); row("Impulse given back (N s)", b20["impulse"], B["impulse"], "%.4e"); row("Energy the sled gives up (TJ)", b20["E_mech"]/1e12, B["E_mech"]/1e12, "%.1f")
P("   The limit does not bind (G-load limited: %s): the outward load at the last step before the release is %.5f g, and the launch ends at the target speed.  The thrust is not throttled." % (C["g_limited"], a_out_last))
P("   The release is %.2f m/s past the target, one time step's gain, so the outward load at the release reads %.4f g." % (C["v_end"] - v3_cas, a_out_end))
P("   The vector sum of the outward load and the %.3f g along the track at the release is %.3f g.  It is printed for the record and limits nothing." % (C["F_end"]/C["m_total"]/cfg.G_0, C["g_peak"]))
P("   Hyperbolic excess: escape speed %.1f m/s, excess %.2f km/s.  Chapter 10 prints 16.3 km/s, 3.0 g and 12.3 km/s." % (v_esc, v_inf/1e3))
OUT["crewed"] = strip(C); OUT["crewed"].update(v3_space=v3_space, v3_casing=v3_cas, v_escape=v_esc, v_inf=v_inf, a_out_last_step_g=a_out_last)

P("\n2. The casing's swing after the crewed launch.  Hold %.1f N/m ; LIMs damping at the critical rate, limited to %.1f MN." % (K, F_LIM/1e6))
OUT["crewed"]["casing"] = casing_rows(C["F"], C["FB"], daily=False)
a = c20["casing"]; b = OUT["crewed"]["casing"]
for nm in ("not braked", "braked"):
    P("   %-10s: farthest %.1f m (round 20: %.1f) ; when %.2f h (%.2f) ; within 1 %% after %.1f h (%.1f) ; LIM push %.2f MN (%.2f) ; reaction line %.3f kN (%.3f)."
      % (nm, b[nm]["x_max"], a[nm]["x_max"], b[nm]["peak_h"], a[nm]["peak_h"], b[nm]["die_h"], a[nm]["die_h"], b[nm]["F_lim_MN"], a[nm]["F_lim_MN"], b[nm]["line_kN"], a[nm]["line_kN"]))
P("   The book's Section 11.6 prints 56 minutes and 1.1 km for this launch: %.2f minutes and %.0f m." % (C["dur"]/60, b["not braked"]["x_max"]))
with open(os.path.join(LOG_DIR, "round21_crew_limit.txt"), "w", encoding="utf-8") as fh: fh.write(buf.getvalue())
with open(os.path.join(LOG_DIR, "round21_crew_limit.json"), "w", encoding="utf-8") as fh: json.dump(OUT, fh, indent=1, default=float)
print("ok")
