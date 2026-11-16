#!/usr/bin/env python3
"""
Round 9, item 5b, part B: the stray field of the bearing with plates that saturate (Section 7.9).
Book: Orbital Ring Engineering, Volume III: Chapter 7, Section 7.9 (18 mT one meter below a support track, under 1 mT at the LIM stators, about 12 mT at a
stator in a launch, about 50 mT below a support track at its cap: round 8's figures, bearing_round8_parts.stray, linear plates).
Brief: Vol_III/Claude_Code_instructions_round9_reruns_2026-11-16.md, Section 5b.
Run: python tools/round9_stray_field.py                 (all three parts in turn; one to two hours on a loaded machine)
     python tools/round9_stray_field.py --part=main     (the six states, linear and saturating, 30 x 30 m box, 50 mm cells; the as-drawn layout, linear)
     python tools/round9_stray_field.py --part=grid     (the grid check: a 16 x 16 m box at 50 and 25 mm, linear and saturating, three states)
     python tools/round9_stray_field.py --part=ends     (part C, Section 7.8: how round 8's coil-end estimate moves with its one assumed length; no new model)
     python tools/round9_stray_field.py --part=report   (prints the tables from runs/round9_logs/round9_stray.json and writes round9_stray.txt)
Reads runs/round8_logs/round8_rows.json and runs/round7_logs/round7_half_2p5.json (the ampere-turns of the six states and round 8's linear figures).
Writes runs/round9_logs/round9_stray.json and round9_stray.txt and nothing else.

The method is round 8's (bearing_round8_parts.stray, which is final_bearing D.2's on the round 7 section): all six tracks in one box, A = 0 on its edge.
The saturating plates follow the rule of the bearing solves (bearing_electrical.Track._iron: Grid.add_iron(..., cf.PLATE_MU, B_sat=cf.PLATE_B_SAT), solved by
magnetostatic_2d.Grid.solve_nonlinear: B = mu0 H + B_sat tanh(mu0 (mu_r - 1) H / B_sat), mu_r 200, B_sat 1.2 T).  solve_nonlinear stops after four passes
when it is given no force box, so it is called again and again here until the probe fields and the permeability stop changing.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the repository root
import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_round7 as r7
import compensation_fields as cf
import final_bearing as fb
import magnetostatic_2d as ms

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS = os.path.join(HERE, "runs", "round9_logs")
JSON = os.path.join(LOGS, "round9_stray.json")
SUP, FAST = r7.SUP, r7.FAST
SIDE = ring.CABLE_SIDE
X_C = SIDE / 2 - FAST["Wp"] / 2                                   # as bearing_round7_parts.X_C, X_S
X_S = SIDE / 2 - FAST["Wp"] - 0.2 - SUP["Wp"] / 2
SUP_BW, SUP_BH = 0.15, 0.067                                      # the support bundle, as round 8's stray()
STANDOFF = bd.GAP0 + 0.010                                        # 110 mm from a plate's face to its bundles (the gap and 10 mm of cryostat)
BOX30 = (30.0, 30.0, -15.0, -10.0)                                # round 8's box: Lx, Ly, x0, y0 (the cable's underside at y = 0)
BOX16 = (16.0, 16.0, -8.0, -4.3)
# the LIM stator: round 8's assumption (final_bearing.STATOR_X, STATOR_Y: a strip 0.2 m wide and 2 m tall) and the body drawn in tools/round7_figures.py
# (draw(): Rectangle((x_in, 0.4), 0.3, 1.0) with x_in = SIDE/2 + 0.4, in coordinates from the cable's center)
R8_ST = dict(x0=fb.STATOR_X - 0.1, x1=fb.STATOR_X + 0.1, y0=fb.STATOR_Y - 1.0, y1=fb.STATOR_Y + 1.0, xc=fb.STATOR_X, yc=fb.STATOR_Y)
DR_ST = dict(x0=SIDE / 2 + 0.4, x1=SIDE / 2 + 0.7, y0=SIDE / 2 + 0.4, y1=SIDE / 2 + 1.4, xc=SIDE / 2 + 0.55, yc=SIDE / 2 + 0.9)


def states():
    d8 = json.load(open(os.path.join(HERE, "runs", "round8_logs", "round8_rows.json")))
    half = json.load(open(os.path.join(HERE, "runs", "round7_logs", "round7_half_2p5.json")))
    NI_n = half["normal"]["NI"]; NI_cap = half["cap"]["NI"]; NI70 = d8["rating_support"]["NI"]
    f74 = d8["t74_fast"]; NI_lift = f74[0]["NI"]; NI198 = f74[1]["NI"]; NI_lift70 = f74[2]["NI"]
    st = [("normal operation (support at 47.5 kN/m, fast tracks off)", NI_n, 0.0, 0.0),
          ("a launch at 30 km/s: support normal, the lift tracks at 196.5 kN/m", NI_n, 0.0, NI_lift),
          ("the lift tracks at their rating, support normal", NI_n, 0.0, NI_lift70),
          ("the load tracks at 198 kN/m, support normal", NI_n, NI198, 0.0),
          ("the support tracks at their cap (245.5 kN/m), fast off", NI_cap, 0.0, 0.0),
          ("the support tracks at their rating, fast off", NI70, 0.0, 0.0)]
    return st, d8["stray"]


class Section:
    """The six tracks in one box.  layout 'let in' is round 8's (the plates let into the cable's faces); 'as drawn' has the plates standing proud of
    the cable's faces, as tools/round7_figures.py draws them, with the bundles the same 110 mm from the plates."""

    def __init__(self, h, box=BOX30, sat=False, layout="let in", air=False):
        self.h, self.sat, self.layout = h, sat, layout
        g = self.g = ms.Grid(box[0], box[1], h, x0=box[2], y0=box[3])
        bs = cf.PLATE_B_SAT if sat else None
        proud = layout == "as drawn"
        ys0, ys1 = (-SUP["T"], 0.0) if proud else (0.0, SUP["T"])                 # support plate
        yl0, yl1 = (-FAST["T"], 0.0) if proud else (0.0, FAST["T"])               # load plate
        yt0, yt1 = (SIDE, SIDE + FAST["T"]) if proud else (SIDE - FAST["T"], SIDE)  # lift plate
        self.y_sup = ys0 - STANDOFF; self.y_load = yl0 - STANDOFF; self.y_lift = yt1 + STANDOFF      # the bundles' faces toward their plates
        self.plates = {}
        for sx, side in ((-1, "left"), (1, "right")):
            xs = sx * X_S; xc = sx * X_C
            for name, (a, b, c, d) in (("support", (xs - SUP["Wp"] / 2, xs + SUP["Wp"] / 2, ys0, ys1)), ("load", (xc - FAST["Wp"] / 2, xc + FAST["Wp"] / 2, yl0, yl1)),
                                       ("lift", (xc - FAST["Wp"] / 2, xc + FAST["Wp"] / 2, yt0, yt1))):
                g.add_iron(a, b, c, d, cf.PLATE_MU, B_sat=bs)
                self.plates[(name, side)] = (a, b, c, d)
        if air:                                                    # no plates at all: the bare coils, for the bracket
            g.mu_r[:] = 1.0

    def set_currents(self, NI_sup, NI_load, NI_lift):
        g = self.g; g.J[:] = 0.0
        d = (SUP["Wl"] + 0.15 + 0.05) / 2
        for y0, y1 in ((self.y_sup - SUP_BH, self.y_sup), (self.y_load - FAST["bh"], self.y_load), (self.y_lift, self.y_lift + FAST["bh"])):
            if not ((g.y >= y0) & (g.y < y1)).any():                # Grid.add_current drops a bundle that holds no node, without a word
                raise SystemExit(f"the grid of {self.h*1e3:g} mm has no node row inside the bundle from y = {y0:.3f} to {y1:.3f} m")
        for sx in (-1, 1):
            xs = sx * X_S; xc = sx * X_C
            if NI_sup:
                for dx, a in ((-d - 0.2, NI_sup), (-d + 0.2, -NI_sup), (d - 0.2, -NI_sup), (d + 0.2, NI_sup)):
                    g.add_current(xs + dx - SUP_BW / 2, xs + dx + SUP_BW / 2, self.y_sup - SUP_BH, self.y_sup, a)
            for NI, y0, y1 in ((NI_load, self.y_load - FAST["bh"], self.y_load), (NI_lift, self.y_lift, self.y_lift + FAST["bh"])):
                if NI:
                    for dx, a in ((-0.3, NI), (0.3, -NI)):
                        g.add_current(xc + dx - FAST["bw"] / 2, xc + dx + FAST["bw"] / 2, y0, y1, a)

    def probes(self):
        g = self.g; B = np.sqrt(g.Bx**2 + g.By**2)
        o = {"1 m below a support track": math.hypot(*g.field_at(X_S, self.y_sup - SUP_BH - 1.0)), "1 m below a load track": math.hypot(*g.field_at(X_C, self.y_load - FAST["bh"] - 1.0))}
        for tag, st in (("round 8 stator", R8_ST), ("drawn stator", DR_ST)):
            for sx, side in ((1, "right"), (-1, "left")):
                x0, x1 = (st["x0"], st["x1"]) if sx > 0 else (-st["x1"], -st["x0"])
                m = g.box(x0, x1 + 1e-9, st["y0"], st["y1"] + 1e-9)
                o[f"{tag}, {side}, center"] = math.hypot(*g.field_at(sx * st["xc"], st["yc"]))
                o[f"{tag}, {side}, largest"] = float(B[m].max())
        # round 8's own two stator figures, its own mask (x from X - 0.1 to X + 0.1, y over 2 m, upper edges open)
        o["LIM stator, center"] = math.hypot(*g.field_at(fb.STATOR_X, fb.STATOR_Y))
        o["largest over the stator"] = float(B[g.box(fb.STATOR_X - 0.1, fb.STATOR_X + 0.1, fb.STATOR_Y - 1.0, fb.STATOR_Y + 1.0)].max())
        return o

    def plate_state(self):
        g = self.g; h = self.h; o = {}
        for (name, side), (a, b, c, d) in self.plates.items():
            if side != "right":
                continue
            m = g.box(a, b, c, d)
            o[name] = dict(sat=float(np.mean(g.mu_r[m] < 0.1 * cf.PLATE_MU)), B_max=float(np.sqrt(g.Bx[m]**2 + g.By[m]**2).max()), mu_mean=float(g.mu_r[m].mean()),
                           pull=abs(float(g.force_on_box(a - 0.1, b + 0.1, c - 2 * h, d + 2 * h))))
        return o

    def solve(self, tol=1e-3, max_calls=14, verbose=True):
        g = self.g; t0 = time.time()
        if not self.sat:
            g.solve()
            return dict(solves=1, seconds=time.time() - t0)
        keys = ("1 m below a support track", "1 m below a load track", "round 8 stator, right, center", "round 8 stator, right, largest", "drawn stator, right, largest")
        prev = None; d = 1.0; n = 0
        for c in range(max_calls):
            n += g.solve_nonlinear(max_outer=6) + 1
            p = self.probes(); cur = np.array([p[k] for k in keys])
            if prev is not None:
                d = float(np.max(np.abs(cur - prev) / np.maximum(np.abs(cur), 1e-9)))
            prev = cur
            if verbose:
                print(f"      call {c+1}: permeability change {g.outer_change:.5f}, probe change {d:.5f}, below support {cur[0]*1e3:.3f} mT, stator {cur[3]*1e3:.3f} mT ({time.time()-t0:.0f} s)", flush=True)
            if d < tol and g.outer_change < 0.002:
                break
        return dict(solves=n, calls=c + 1, probe_change=d, mu_change=float(g.outer_change), seconds=time.time() - t0)


def merge(key, val):
    """Put one part's results under its key in the JSON.  The two parts may run as two processes, so the read is retried and the write checked."""
    os.makedirs(LOGS, exist_ok=True)
    for _ in range(5):
        d = {}
        for _ in range(10):
            try:
                d = json.load(open(JSON, encoding="utf-8")) if os.path.exists(JSON) else {}
                break
            except (json.JSONDecodeError, OSError):
                time.sleep(1.0)
        d[key] = val
        with open(JSON, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=1)
        time.sleep(0.5)
        try:
            if key in json.load(open(JSON, encoding="utf-8")):
                return
        except (json.JSONDecodeError, OSError):
            pass


def run_set(h, box, sat, layout, st, which=None, tag="", air=False):
    sec = Section(h, box, sat, layout, air=air); rows = {}
    for name, NI_s, NI_ld, NI_lf in st:
        if which is not None and name not in which:
            continue
        print(f"    {tag}{'no plates' if air else 'saturating' if sat else 'linear'}, {layout}, {box[0]:.0f} m box, {h*1e3:.1f} mm: {name}", flush=True)
        sec.set_currents(NI_s, NI_ld, NI_lf)
        info = sec.solve()
        rows[name] = dict(NI_sup=NI_s, NI_load=NI_ld, NI_lift=NI_lf, fields=sec.probes(), plates=sec.plate_state(), solve=info)
    return rows


def part_main():
    st, r8 = states()
    out = dict(geometry=dict(SIDE=SIDE, X_S=X_S, X_C=X_C, round8_stator=R8_ST, drawn_stator=DR_ST, box=BOX30, h=0.05, plate_mu=cf.PLATE_MU, plate_bsat=cf.PLATE_B_SAT))
    lin = run_set(0.05, BOX30, False, "let in", st)
    worst = 0.0
    for name, row in lin.items():
        for k, v in r8[name].items():
            worst = max(worst, abs(row["fields"][k] - v) / v)
    out["linear"] = lin; out["round8_repro_worst_rel_diff"] = worst
    print(f"    linear rows against runs/round8_logs/round8_rows.json: largest relative difference {worst:.2e}", flush=True)
    out["no_plates"] = run_set(0.05, BOX30, False, "let in", st, air=True)
    merge("main", out)
    out["saturating"] = run_set(0.05, BOX30, True, "let in", st)
    merge("main", out)
    out["as_drawn_linear"] = run_set(0.05, BOX30, False, "as drawn", st)
    merge("main", out)


def part_grid():
    st, _ = states()
    which = [st[0][0], st[4][0], st[1][0]]                           # normal, the cap, a launch
    out = dict(box=BOX16, states=which)
    for h in (0.05, 0.025):
        for sat in (False, True):
            out[f"{'sat' if sat else 'lin'}_{h*1e3:g}mm"] = run_set(h, BOX16, sat, "let in", st, which, tag="grid check: ")
            merge("grid", out)


def diffusion_half_cycle(B, t_edge, d_p, mu_r, rho_p, nz, t_after=15e-3):
    """bearing_design.eddy_cycle with no hold, halved, on a grid of nz cells through the plate (eddy_cycle itself has 100): the eddy energy of one
    rise or fall of the face field B over the time t_edge, J/m2, by the same explicit 1-D diffusion."""
    dz = d_p / nz; D = rho_p / (ms.MU0 * mu_r); dt = 0.4 * dz * dz / D
    n = int((2 * t_edge + t_after) / dt) + 1
    Bz = np.zeros(nz + 1); E = 0.0
    for j in range(n):
        t = j * dt
        Bz[0] = B * t / t_edge if t < t_edge else (B * (1 - (t - t_edge) / t_edge) if t < 2 * t_edge else 0.0); Bz[-1] = Bz[-2]
        Jc = (Bz[1:] - Bz[:-1]) / dz / (ms.MU0 * mu_r)
        E += rho_p * float(np.sum(Jc * Jc)) * dz * dt
        Bz[1:-1] += D * (Bz[2:] - 2 * Bz[1:-1] + Bz[:-2]) / dz**2 * dt
    return 0.5 * E


def part_ends():
    """Part C: round 8's coil-end estimate (bearing_round8_parts.part_rows, section 7) again, with the length of the step varied and the diffusion grid refined.
    The model is unchanged: a coil end takes a tenth of its own loop's face field away over a length ell, the heat is half a rise-and-fall cycle of the
    1-D diffusion, summed across the plate with the one-loop face field squared.  Round 8 takes ell = 2 x (gap + 50 mm) = 300 mm."""
    import levitation_config as lc
    d8 = json.load(open(os.path.join(HERE, "runs", "round8_logs", "round8_rows.json")))
    f2 = d8["face"]["f2_one"]; nn = 10; Lc = 1000.0; steps = 2 * 2 * nn / Lc
    rows = []
    for ell in (0.15, 0.30, 0.60):
        t_edge = ell / bd.V_REL
        E_bd = 0.5 * bd.eddy_cycle(1.0, t_edge, 0.0, d_p=SUP["T"], mu_r=lc.FERRO_MU_R, rho_p=lc.FERRO_RESISTIVITY)
        E_fine = diffusion_half_cycle(1.0, t_edge, SUP["T"], lc.FERRO_MU_R, lc.FERRO_RESISTIVITY, 800)
        skin = math.sqrt(lc.FERRO_RESISTIVITY * t_edge / (ms.MU0 * lc.FERRO_MU_R))
        r = dict(ell=ell, t_edge=t_edge, skin_m=skin, E1_100=E_bd, E1_800=E_fine)
        for k, E in (("100", E_bd), ("800", E_fine)):
            P = 2 * steps * E * f2 / nn**2 * bd.V_REL
            r[f"P_{k}"] = P; r[f"MW_{k}"] = P * ring.L_RING / 1e6
        rows.append(r)
        print(f"    ell {ell*1e3:.0f} mm: transit {t_edge*1e6:.1f} microseconds, sqrt(rho t / mu) {skin*1e3:.2f} mm; E1 {E_bd:.3f} J/m2 (100 cells), {E_fine:.3f} (800 cells); {r['P_100']:.2f} and {r['P_800']:.2f} W/m", flush=True)
    merge("ends_sensitivity", dict(rows=rows, f2_one=f2, round8_E1=d8["face"]["E1"], round8_P=d8["ends"][0]["P"], round8_MW=d8["ends"][0]["MW"], V_REL=bd.V_REL, mu_r=lc.FERRO_MU_R, rho=lc.FERRO_RESISTIVITY, plate_T=SUP["T"]))


def part_report():
    d = json.load(open(JSON)); lines = []
    def out(s=""):
        print(s); lines.append(s)
    m = d["main"]; g = m["geometry"]; lin, sat = m["linear"], m.get("saturating", {})
    out("=" * 150); out("  ROUND 9, ITEM 5b, PART B: THE STRAY FIELD WITH SATURATING PLATES (SECTION 7.9)"); out("=" * 150)
    out(f"  All six tracks in one {g['box'][0]:.0f} x {g['box'][1]:.0f} m box, {g['h']*1e3:.0f} mm cells, A = 0 on the edge (round 8's method, bearing_round8_parts.stray).  Plates: mu_r {g['plate_mu']:.0f}, saturating at")
    out(f"  {g['plate_bsat']} T by magnetostatic_2d.Grid.solve_nonlinear.  Cable side {g['SIDE']:.4f} m, its underside at y = 0; support tracks at x = +-{g['X_S']:.3f} m, load and lift tracks at +-{g['X_C']:.3f} m.")
    out(f"  The linear rows reproduce round 8's saved figures to a relative difference of {m['round8_repro_worst_rel_diff']:.1e}.")
    out()
    out("  Table B1 (round9_stray_field.py, part_main): linear against saturating plates, millitesla.  'stator' is round 8's strip, 0.2 m wide and 2 m tall, on the right.")
    out(f"  {'state':>68} | {'1 m below support':>19} | {'1 m below load':>19} | {'stator center':>19} | {'largest over stator':>19}")
    out(f"  {'':>68} | {'linear':>9} {'saturating':>9} | {'linear':>9} {'saturating':>9} | {'linear':>9} {'saturating':>9} | {'linear':>9} {'saturating':>9}")
    keys = ("1 m below a support track", "1 m below a load track", "LIM stator, center", "largest over the stator")
    for name in lin:
        s = f"  {name:>68} |"
        for k in keys:
            a = lin[name]["fields"][k] * 1e3; b = sat[name]["fields"][k] * 1e3 if name in sat else float("nan")
            s += f" {a:9.2f} {b:9.2f} |"
        out(s)
    out()
    out("  Table B2 (part_main, Section.plate_state): the plates in the saturating run, right-hand side.  'sat' is the share of the plate's cells under a tenth of mu_r; 'pull' is the Maxwell")
    out("  stress on a box around the plate on this 50 mm grid, to set against the fine solves' 47.5 kN/m (normal), 245.5 (cap), 196.5 (lift) and 198 (load).")
    out(f"  {'state':>68} | {'support: sat %':>14} {'peak T':>7} {'pull kN/m':>9} | {'load: sat %':>11} {'peak T':>7} {'pull kN/m':>9} | {'lift: sat %':>11} {'peak T':>7} {'pull kN/m':>9} | {'solves':>6} {'probe chg':>9}")
    for name in sat:
        p = sat[name]["plates"]; sv = sat[name]["solve"]
        out(f"  {name:>68} | {p['support']['sat']*100:14.1f} {p['support']['B_max']:7.2f} {p['support']['pull']/1e3:9.1f} | {p['load']['sat']*100:11.1f} {p['load']['B_max']:7.2f} {p['load']['pull']/1e3:9.1f} | "
            f"{p['lift']['sat']*100:11.1f} {p['lift']['B_max']:7.2f} {p['lift']['pull']/1e3:9.1f} | {sv['solves']:6d} {sv.get('probe_change', 0):9.1e}")
    out(f"  The same pulls with linear plates on this grid: " + "; ".join(f"{n.split(' (')[0].split(':')[0]}: support {lin[n]['plates']['support']['pull']/1e3:.1f}, load {lin[n]['plates']['load']['pull']/1e3:.1f}, lift {lin[n]['plates']['lift']['pull']/1e3:.1f}" for n in lin))
    out()
    r, dr = g["round8_stator"], g["drawn_stator"]; c = g["SIDE"] / 2
    out("  Table B3 (part_main, Section.probes): the field at the LIM stator, round 8's assumed position against the body drawn in tools/round7_figures.py, millitesla, saturating plates")
    out(f"  (linear in brackets).  Round 8: a strip from x = {r['x0']:.3f} to {r['x1']:.3f} m, from {r['y0']-c:+.2f} to {r['y1']-c:+.2f} m about the cable's centerline (center {r['xc']:.3f} m, {r['yc']-c:+.2f} m).")
    out(f"  As drawn: a body from x = {dr['x0']:.3f} to {dr['x1']:.3f} m, from {dr['y0']-c:+.2f} to {dr['y1']-c:+.2f} m about the centerline (center {dr['xc']:.3f} m, {dr['yc']-c:+.2f} m).  Both sides of the cable.")
    out(f"  {'state':>68} | {'round 8 strip, right':>22} {'round 8 strip, left':>22} | {'drawn body, right':>22} {'drawn body, left':>22} | {'drawn center, right':>19}")
    for name in lin:
        s = f"  {name:>68} |"
        for k in ("round 8 stator, right, largest", "round 8 stator, left, largest", "drawn stator, right, largest", "drawn stator, left, largest", "drawn stator, right, center"):
            a = lin[name]["fields"][k] * 1e3; b = sat[name]["fields"][k] * 1e3 if name in sat else float("nan")
            s += f" {b:10.2f} ({a:8.2f}) " + ("|" if k in ("round 8 stator, left, largest", "drawn stator, left, largest") else "")
        out(s)
    if "as_drawn_linear" in m:
        ad = m["as_drawn_linear"]
        out()
        out("  Table B4 (part_main, layout 'as drawn'): the same six states with the plates standing proud of the cable's faces, as tools/round7_figures.py draws them (round 8 lets them into")
        out("  the cable: the lift bundles are then 0.2 m nearer the cable's centerline, the support bundles 0.1 m and the load bundles 0.2 m).  Linear plates, millitesla.")
        out(f"  {'state':>68} | {'1 m below support':>17} | {'1 m below load':>14} | {'drawn body, right':>17} | {'drawn body, left':>16} | {'round 8 strip, right':>20}")
        for name in ad:
            f = ad[name]["fields"]
            out(f"  {name:>68} | {f['1 m below a support track']*1e3:17.2f} | {f['1 m below a load track']*1e3:14.2f} | {f['drawn stator, right, largest']*1e3:17.2f} | {f['drawn stator, left, largest']*1e3:16.2f} | {f['round 8 stator, right, largest']*1e3:20.2f}")
    if "grid" in d:
        gr = d["grid"]; out()
        out(f"  Table B5 (part_grid): the grid check.  A {gr['box'][0]:.0f} x {gr['box'][1]:.0f} m box at 50 and 25 mm cells, and the 30 m box at 50 mm beside it for the effect of the box.  Millitesla; pulls in kN/m.")
        out(f"  {'state':>68} {'plates':>10} {'box m':>5} {'cell mm':>7} | {'1 m below support':>17} {'1 m below load':>14} {'stator center':>13} {'largest over stator':>19} | {'support pull':>12} {'support sat %':>13} {'lift pull':>9}")
        for name in gr["states"]:
            for kind, src in (("linear", "lin"), ("saturating", "sat")):
                rows = [(30, 50, (lin if src == "lin" else sat).get(name))] + [(gr["box"][0], hh, gr.get(f"{src}_{hh:g}mm", {}).get(name)) for hh in (50, 25)]
                for bx, hh, row in rows:
                    if not row:
                        continue
                    f = row["fields"]; p = row["plates"]
                    out(f"  {name:>68} {kind:>10} {bx:5.0f} {hh:7.0f} | {f['1 m below a support track']*1e3:17.2f} {f['1 m below a load track']*1e3:14.2f} {f['LIM stator, center']*1e3:13.2f} {f['largest over the stator']*1e3:19.2f} | "
                        f"{p['support']['pull']/1e3:12.1f} {p['support']['sat']*100:13.1f} {p['lift']['pull']/1e3:9.1f}")
    if "ends_sensitivity" in d:
        e = d["ends_sensitivity"]; out()
        out("  PART C (SECTION 7.8).  Table C1 (part_ends): round 8's coil-end estimate with its one assumed length varied.  The model is round 8's and is not improved: a coil end takes a tenth of")
        out(f"  its loop's face field away over a length ell; half a rise-and-fall cycle of the 1-D diffusion (plate {e['plate_T']*1e3:.0f} mm, mu_r {e['mu_r']:.0f}, {e['rho']:.1e} ohm m); 1 km coils, ten per loop, staggered.")
        out(f"  Round 8: ell = 300 mm, {e['round8_E1']:.3f} J/m2 for 1 T, {e['round8_P']:.2f} W per meter of ring, {e['round8_MW']:.0f} MW.  '100 cells' is bearing_design.eddy_cycle as round 8 called it; '800 cells' the same")
        out("  diffusion on a finer grid through the plate.")
        out(f"  {'ell mm':>7} {'transit microseconds':>20} {'sqrt(rho t / mu) mm':>19} | {'J/m2 for 1 T, 100 cells':>23} {'W/m':>6} {'MW':>6} | {'J/m2 for 1 T, 800 cells':>23} {'W/m':>6} {'MW':>6}")
        for r in e["rows"]:
            out(f"  {r['ell']*1e3:7.0f} {r['t_edge']*1e6:20.1f} {r['skin_m']*1e3:19.2f} | {r['E1_100']:23.3f} {r['P_100']:6.2f} {r['MW_100']:6.0f} | {r['E1_800']:23.3f} {r['P_800']:6.2f} {r['MW_800']:6.0f}")
    with open(os.path.join(LOGS, "round9_stray.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    part = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--part=")), "all")
    if part in ("main", "all"):
        part_main()
    if part in ("grid", "all"):
        part_grid()
    if part in ("ends", "all"):
        part_ends()
    if part in ("report", "all"):
        part_report()


if __name__ == "__main__":
    main()
