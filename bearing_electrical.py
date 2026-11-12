#!/usr/bin/env python3
"""
The electrical design of the levitation bearing's driven coils: the voltage a turn sees, series turns against parallel
tapes and the insulation they need, current sharing, leads and supplies, and the coupling between the driven load
compensation winding and the persistent support winding that share the same loops (Part 1 of
Claude_Code_instructions_bearing_electrical_2026-11-12.md).
Book: Orbital Ring Engineering, Volume I: Chapter 7 (Sections 7.5 to 7.7: the driven windings, Table 7.6, the supply ratings).
Run: python bearing_electrical.py --part=1a|1b|1c|1d|all [--h=5|2.5] [--quick]
     python bearing_electrical.py --part=sweep        (Part 1c's protection-rise and lift-coil-length sweep of final_bearing's cases; slow)
Environment: ORBITAL_RING_OUT (output folder, default runs); the ring through ring_altitude.py as the other bearing scripts

Nothing here changes a default of bearing_optimize.py, final_bearing.py, bearing_design.py or magnetostatic_2d.py.  The best
tracks of the Part 3b optimization (two 0.5 m loops with four 150 x 67 mm bundles and a 10 mm wall; one 0.8 m lift loop with
two 200 x 50 mm bundles) are rebuilt with their lattices and solved again on the same windowed grids, with the currents of
each bundle free, so that one loop can be excited alone and the flux through each loop read from the vector potential.

Models
  Flux linkage.  The 2-D solver gives A_z; the flux per meter through a loop is the difference of the mean vector potential
  over the two bundles that form it (the mean over the tape lattice of the bundle).  For the whole track the sum of
  NI_leg x A_leg / 2 is the stored energy, so the two routes to the flux are the same quantity, and the report prints both.
  A turn of a driven coil of length L_coil sees the voltage L_coil x dPhi'/dt, where Phi' is the flux per meter and dPhi'/dt
  follows the ampere-turns of all the families in the bundle (they share the loop), so when every family rises together
  over the coil's length the turn sees the full change of the loop's flux.
  Coupled circuits (Part 1d).  A persistent support coil of length L_c with n turns is a closed superconducting loop: its
  flux linkage n x integral(Phi'(x) dx) over its length is constant on every time scale of the driven windings (its L/R with
  soldered joints is years).  With the two windings on the same loops, a change dNI_l(x) of the load winding over a
  stretch forces the support coils covering that stretch to change their current so that each coil's integral is kept.
  The ten staggered coils per bundle overlap, so the system is linear in the coil currents and is solved as a banded
  linear system along the ring; the pull then follows from the local total ampere-turns through the track's pull curve
  (bearing_optimize Part 3b).  The gap dependence of the loop inductance and of the pull (a gap sweep of linear solves)
  gives the track's stiffness at constant current and at constant flux.
  Sled passage (Part 1d).  A point model of one meter of casing and cable, as sled_passage_pulsed's point model with the
  2026-10-11 options (spread pull over max(50 ms, L_sc/u), 20 staggered coils, the 20 ms command ramp, the 2 ms lag, the
  10 mm protection), re-implemented here so that the support pull can follow the coupled coil currents (the load winding's
  moving footprint convolved with the coil system's response kernel).

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""
import contextlib
import io
import json
import math
import os
import sys
import time

import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_mass as bm
import bearing_optimize as bo
import compensation_fields as cf
import final_bearing as fb
import levitation_config as lc
import magnetostatic_2d as ms
import sled_passage_pulsed as pz

MU0 = ms.MU0
GAP = bd.GAP0                       # 0.100 m
PLATE_T, PLATE_W = fb.PLATE_T, fb.PLATE_W
TAU_HW = bd.TAU_HW                  # 2 ms, the protection rise
RAMP = fb.RAMP                      # 20 ms command ramp; the planned rise is RAMP + TAU_HW as final_bearing rates the supplies
T_PLAN = RAMP + TAU_HW
L_LIFT, L_LOAD, N_FAM = fb.L_LIFT, fb.L_LOAD, fb.N_FAM
N_SUP_COILS = 10                    # persistent support coils per bundle (Chapter 7)
L_SUP_COILS = (1000.0, 10000.0)     # support coil lengths to model
TAPE_T, TAPE_RHO = bo.TAPE_T, bo.TAPE_RHO
I_REF_12, I_REF_3 = bo.I_REF_12, bo.I_REF_3
COARSE_H = 0.01
W_HALF = ring.M_LOAD_M * ring.G_NET / 2          # 47.5 kN/m, half the casing weight: one support track in normal operation
F_LOAD = fb.F_LOAD_TRACK                          # 198 kN/m
F_LIFT = fb.F_LIFT_TRACK                          # 196.5 kN/m
PASS_PULLS = bo.PASS_PULLS                        # (km/s, kN/m on the pair, winding)
V_ORBIT = 6988.0                                  # m/s relative to the casing: the sled presses on nothing (Chapter 7)
L_SC = bd.L_SC                                    # 100 m spacecraft zone
# insulation design stresses, kV/mm, with their basis (Part 1b); the sources are listed in the report
E_DESIGN = {"5 (bubble-limited)": 5.0, "10 (PD-free impregnated)": 10.0, "20 (film, bubble-free)": 20.0}
F_CABLE = 0.75                                    # tape stack fill inside a cable (tape thickness over pitch) [assumption]
F_COOL = 0.15                                     # share of the bundle area left as coolant channels [assumption]
RHO_CU_77 = 2.2e-9                                # ohm m, copper at 77 K, RRR about 50 [UNVERIFIED: handbook value]
L0 = 2.44e-8                                      # W ohm / K^2, the Lorenz number
T_WARM, T_COLD = 300.0, 70.0


def hr(t=""):
    fb.hr(t)


# =============================================================================
# the chosen tracks
# =============================================================================

def load_records():
    """The Part 3b records of the best support and lift tracks (pull curves, tape limits, turns)."""
    out = {}
    for kind, geom in (("support", bo.BEST_SUPPORT), ("lift", bo.BEST_LIFT)):
        for folder in (ring.out_dir(), os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "Vol_III", "runs_coil_optimization_2026-10-13")):
            p = os.path.join(folder, f"bearing_optimize_part3b_{kind}_2.5mm.json")
            if os.path.exists(p):
                with open(p, encoding="utf-8") as f:
                    out[kind] = json.load(f)[bo.geom_key(geom)]
                break
        if kind not in out:
            raise SystemExit(f"the Part 3b cache for the {kind} track is missing (bearing_optimize_part3b_{kind}_2.5mm.json)")
    return out["support"], out["lift"]


def load_theva():
    theva = dict(bd.read_theva()); theva["full"] = bo.read_robinson_angles()
    theva77 = bo.read_theva_file(bo.THEVA_77)
    return theva, theva77


class Track:
    """A track of bearing_optimize's geometry with the current of every bundle free.  Legs (x, design sign) as
    build_geometry makes them; loops as (leg indices, own-sense signs): loop 1 is the left pair, loop 2 the right pair,
    the own sense being the design's polarity, so at the design NI_1 = NI_2 = NI and both loop fluxes are equal."""

    def __init__(self, geom, h, gap=GAP, w_tape=0.012, I_ref=I_REF_12, NI_ref=322e3, window=bo.WINDOW):
        self.geom, self.h, self.gap, self.w_tape, self.window = geom, h, gap, w_tape, window
        n, W, gc, shape, allow = geom["n"], geom["W"], geom["gc"], geom["shape"], geom["allow"]
        wb = bo.shape_width(shape); split = 0.0 if gc is None else wb + gc
        x0 = -(n * W + (n - 1) * split) / 2
        legs = []
        for k in range(n):
            sgn = 1 if k % 2 == 0 else -1
            xl = x0 + k * (W + split); xr = xl + W
            legs.append((xl, sgn)); legs.append((xr, -sgn))
        merged = {}
        for x, s in legs:
            merged[round(x, 6)] = merged.get(round(x, 6), 0) + s
        self.legs = sorted(merged.items())                    # (x, design sign)
        self.loops = []
        for k in range(n):
            xl = round(x0 + k * (W + split), 6); xr = round(xl + W, 6)
            il = [i for i, (x, s) in enumerate(self.legs) if x == xl][0]; ir = [i for i, (x, s) in enumerate(self.legs) if x == xr][0]
            self.loops.append(((il, self.legs[il][1]), (ir, self.legs[ir][1])))
        self.allow = allow; self.shape = shape; self.wb = wb; self.hb = bo.shape_height(shape)
        N = max(1, int(round(NI_ref / I_ref)))
        self.lattices = [bo.tape_lattice(shape, x, -allow, w_tape, N) for x, s in self.legs]
        self.solves = 0

    def NI_design(self, NI):
        return [s * NI for x, s in self.legs]

    def NI_loop(self, k, NI):
        """Loop k alone at NI ampere-turns in its own sense."""
        v = [0.0] * len(self.legs)
        for i, s in self.loops[k]:
            v[i] = s * NI
        return v

    def _deposit(self, g, NI_legs):
        for (xs, ys), a in zip(self.lattices, NI_legs):
            if a == 0.0:
                continue
            I = a / len(xs)
            i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
            np.add.at(g.J, (i, j), I / (g.h * g.h))

    def _iron(self, g, linear):
        if linear:
            g.add_iron(-PLATE_W / 2, PLATE_W / 2, self.gap, self.gap + PLATE_T, cf.PLATE_MU)
        else:
            g.add_iron(-PLATE_W / 2, PLATE_W / 2, self.gap, self.gap + PLATE_T, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)

    def solve(self, NI_legs, linear=False, h=None):
        """Coarse full box then the window at h (bearing_optimize's nesting).  Returns the window grid and the summary."""
        h = h or self.h
        g0 = ms.Grid(4.0, 3.0, COARSE_H)
        self._iron(g0, linear); self._deposit(g0, NI_legs)
        box0 = (-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, self.gap - 2 * COARSE_H, self.gap + PLATE_T + 2 * COARSE_H)
        if linear:
            g0.solve()
        else:
            g0.solve_nonlinear(force_box=box0)
        w = bo.Window(g0, h, window=self.window)
        self._iron(w, linear); self._deposit(w, NI_legs)
        box = (-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, self.gap - 2 * h, self.gap + PLATE_T + 2 * h)
        if linear:
            w.solve()
        else:
            w.solve_nonlinear(force_box=box)
        self.solves += 1
        return w, self.summary(w, NI_legs, box)

    def summary(self, g, NI_legs, box):
        A_legs = []
        for (xs, ys) in self.lattices:
            i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
            A_legs.append(float(g.A[i, j].mean()))
        phi = [sum(s * A_legs[i] for i, s in loop) for loop in self.loops]
        W = 0.5 * sum(a * A for a, A in zip(NI_legs, A_legs))
        m = g.box(-PLATE_W / 2, PLATE_W / 2, self.gap, self.gap + PLATE_T)
        sat = float(np.mean(g.mu_r[m] < 0.1 * cf.PLATE_MU)); B_plate = float(np.sqrt(g.Bx[m]**2 + g.By[m]**2).max())
        return dict(NI_legs=list(NI_legs), F=g.force_on_box(*box), A_legs=A_legs, phi=phi, W=W, W_cells=bo.stored_energy(g), sat=sat, B_plate=B_plate)

    def bundles_for_limits(self, NI_legs):
        return [dict(x=x, NI=a, xs=xs, ys=ys) for (x, s), a, (xs, ys) in zip(self.legs, NI_legs, self.lattices)]


def worst_tape(g, track, NI_legs, w_tape, theva, theva77):
    """The worst tape of the modelled winding at the solved state (faces vertical or horizontal, whichever is lower)."""
    lim = bo.tape_limits(g, track.bundles_for_limits(NI_legs), w_tape, theva, theva77)
    best = None
    for d in lim:
        for face in ("v", "h"):
            if best is None or d["ic_" + face] < best[0]:
                best = (d["ic_" + face], d["B_" + face], d["th_" + face], face, d["x"])
    return best, lim


# =============================================================================
# PART 1A: the voltage a turn sees
# =============================================================================

def flux_curve(track, NIs, linear_too=True):
    """Symmetric (design) excitation at each NI: F, loop flux, energy, saturation; and the linear flux per ampere-turn."""
    rows = []
    for NI in NIs:
        t0 = time.time()
        g, s = track.solve(track.NI_design(NI))
        s["NI"] = NI; s["t"] = time.time() - t0
        rows.append(s)
    lin = None
    if linear_too:
        g, s = track.solve(track.NI_design(1e3), linear=True)
        lin = s["phi"][0] / 1e3
    return rows, lin


def phi_at(rows, lin, NI):
    """Loop flux per meter at NI from the solved curve (NI squared-free: flux is near linear in NI; log-log interpolation
    between the points, the linear value below the first)."""
    x = np.array([r["NI"] for r in rows]); y = np.array([r["phi"][0] for r in rows])
    if NI <= 0:
        return 0.0
    if NI < x[0]:
        return float(NI * (lin if lin else y[0] / x[0]))
    return float(math.exp(np.interp(math.log(min(NI, x[-1])), np.log(x), np.log(y))))


def part_1a(h, theva, theva77, r_sup, r_lift, quick=False):
    hr("PART 1A: THE VOLTAGE A TURN SEES (flux linkage from the solver)")
    out = {}
    NI_s0 = bo.NI_of_F(r_sup, W_HALF)                       # the support winding in normal operation, 141 kA
    NI_cap = r_sup["NI_req"]; NI_l = r_sup["NI_load"]; NI_446 = r_sup["NI_446"]
    print(f"  Best support track ({bo.geom_key(bo.BEST_SUPPORT)}), {h*1e3:g} mm cells.  Support winding {NI_s0/1e3:.0f} kA per leg in normal operation,")
    print(f"  cap {NI_cap/1e3:.0f} kA; load winding {NI_l/1e3:.0f} kA for its own {F_LOAD/1e3:.0f} kN/m; design state {NI_446/1e3:.0f} kA.")
    tr = Track(bo.BEST_SUPPORT, h)
    NIs = sorted(set([NI_s0, 2 * NI_s0, NI_l, NI_cap, NI_s0 + (NI_cap - NI_s0), NI_s0 + NI_l, NI_446]))
    if quick:
        NIs = [NI_s0, NI_l, NI_cap, NI_446]
    rows, lin = flux_curve(tr, NIs)
    print(f"  {'NI kA':>6} {'pull kN/m':>9} {'flux/loop Wb/m':>14} {'W kJ/m':>7} {'W/NI':>7} {'sat %':>5} {'s':>5}")
    for r in rows:
        print(f"  {r['NI']/1e3:6.0f} {r['F']/1e3:9.1f} {r['phi'][0]:14.4f} {r['W']/1e3:7.1f} {r['W']/r['NI']:7.4f} {r['sat']*100:5.1f} {r['t']:5.0f}")
    print(f"  Linear flux per ampere-turn (unsaturated plate): {lin*1e6:.4f} uWb/m per A; at {NI_l/1e3:.0f} kA that is {lin*NI_l:.4f} Wb/m against the solved {phi_at(rows, lin, NI_l):.4f}.")
    print(f"  W/NI equals the flux per loop (the stored energy is NI x Phi' for the two loops), so the hand estimate from the energy and the flux")
    print(f"  linkage from A are the same number; the book's 0.408 Wb/m per loop at 287 kA is confirmed.")
    out["support"] = dict(rows=[{k: v for k, v in r.items()} for r in rows], lin=lin, NI_s0=NI_s0)
    # the transitions a load coil sees
    print()
    print(f"  Load winding, 25 m coils, every family rising together over the coil: the change of the loop flux between the start and the end state,")
    print(f"  the average turn voltage over the rise (2 ms protection rise, {T_PLAN*1e3:.0f} ms planned ramp), and the peak over average for a linear")
    print(f"  current ramp (the slope dPhi'/dNI at the start of the rise, from the linear solve, over the average slope).")
    cases = [("load winding alone, 0 to its cap (the book's basis)", 0.0, NI_l),
             ("parked 5,000 t before the pumps: load 0 to 182 kA on the support at 141", NI_s0, NI_cap),
             ("launch at 1 km/s: load 0 to 141 kA on the support at 141", NI_s0, 2 * NI_s0),
             ("protection pulse at the load cap on the support at 141", NI_s0, NI_s0 + NI_l),
             ("design state: load 0 to 171 kA on the support at its cap", NI_cap, NI_446)]
    print(f"  {'case':<78} {'dPhi Wb/m':>9} {'per coil Wb':>11} {'V/turn 2 ms':>11} {'V/turn 22 ms':>12} {'peak/avg':>8}")
    trans = {}
    for name, a, b in cases:
        dphi = phi_at(rows, lin, b) - phi_at(rows, lin, a)
        slope0 = lin if a == 0 else (phi_at(rows, lin, a * 1.02) - phi_at(rows, lin, a * 0.98)) / (0.04 * a)
        pk = slope0 / (dphi / (b - a))
        trans[name] = dict(dphi=dphi, per_coil=dphi * L_LOAD, V2=dphi * L_LOAD / TAU_HW, V22=dphi * L_LOAD / T_PLAN, peak=pk)
        print(f"  {name:<78} {dphi:9.4f} {dphi*L_LOAD:11.2f} {dphi*L_LOAD/TAU_HW/1e3:9.2f} kV {dphi*L_LOAD/T_PLAN/1e3:9.2f} kV {pk:8.2f}")
    out["load_transitions"] = trans
    # the lift track
    print()
    NI_lift = r_lift["NI_req"]; NI_lift_rated = bo.rated_tape(r_lift, 0.7)
    trl = Track(bo.BEST_LIFT, h, w_tape=0.003, I_ref=I_REF_3, NI_ref=NI_lift)
    NIs_l = [NI_lift] if quick else [0.5 * NI_lift, NI_lift, r_lift["rate_tape"]["70"][0]]
    rows_l, lin_l = flux_curve(trl, NIs_l)
    print(f"  Best lift track ({bo.geom_key(bo.BEST_LIFT)}): {NI_lift/1e3:.0f} kA per leg for {F_LIFT/1e3:.1f} kN/m, 50 m coils.")
    print(f"  {'NI kA':>6} {'pull kN/m':>9} {'flux Wb/m':>10} {'W kJ/m':>7} {'2W/NI':>7} {'sat %':>5}")
    for r in rows_l:
        print(f"  {r['NI']/1e3:6.0f} {r['F']/1e3:9.1f} {r['phi'][0]:10.4f} {r['W']/1e3:7.1f} {2*r['W']/r['NI']:7.4f} {r['sat']*100:5.1f}")
    dphi_l = phi_at(rows_l, lin_l, NI_lift)
    print(f"  One loop: 2W/NI is the flux per turn; the hand estimate 2 x 120.3 / 421 = {2*120.3/421:.3f} Wb/m against the solved {dphi_l:.4f}.")
    print(f"  Lift coil 0 to {NI_lift/1e3:.0f} kA: {dphi_l*L_LIFT:.2f} Wb per 50 m coil; {dphi_l*L_LIFT/TAU_HW/1e3:.2f} kV per turn at 2 ms, {dphi_l*L_LIFT/T_PLAN/1e3:.2f} kV at {T_PLAN*1e3:.0f} ms; "
          f"peak over average {lin_l/(dphi_l/NI_lift):.2f}.")
    out["lift"] = dict(rows=rows_l, lin=lin_l, dphi=dphi_l, V2=dphi_l * L_LIFT / TAU_HW, V22=dphi_l * L_LIFT / T_PLAN)
    # the real pulses of a sled passage
    print()
    print(f"  The sled-passage pulses (final_bearing's point model): the pull window of max(50 ms, L_sc/u) sweeps the ring at u; the command is ramped")
    print(f"  over {RAMP*1e3:.0f} ms and the current lags {TAU_HW*1e3:.0f} ms, so the ramp front is u x {T_PLAN*1e3:.0f} ms long along the ring.  A coil whose length L is longer than the")
    print(f"  front has only u x {T_PLAN*1e3:.0f} ms of its length ramping at any instant: the turn sees dPhi' x min(L / {T_PLAN*1e3:.0f} ms, u).  The protection rise (2 ms, every")
    print(f"  family over the stretch at once) is the same at every speed.  Load winding on the support at its normal {NI_s0/1e3:.0f} kA; lift alone.")
    print(f"  {'u km/s':>6} {'winding':>7} {'pull/track kN/m':>15} {'NI_l kA':>7} {'dPhi Wb/m':>9} {'front m':>7} {'V/turn planned':>14} {'V/turn 2 ms':>11} {'coil m':>6}")
    passes = []
    for u_kms, F_pair, kind in PASS_PULLS:
        u = u_kms * 1e3; F_t = F_pair * 1e3 / 2
        if kind == "load":
            NIl = bo.NI_of_F(r_sup, W_HALF + F_t) - NI_s0
            dphi = phi_at(rows, lin, NI_s0 + NIl) - phi_at(rows, lin, NI_s0); L = L_LOAD
        else:
            NIl = bo.NI_of_F(r_lift, F_t)
            dphi = phi_at(rows_l, lin_l, NIl); L = L_LIFT
        front = u * T_PLAN
        V_plan = dphi * min(L / T_PLAN, u); V_prot = dphi * L / TAU_HW
        passes.append(dict(u=u, kind=kind, F_track=F_t, NI_l=NIl, dphi=dphi, front=front, V_plan=V_plan, V_prot=V_prot, L=L))
        print(f"  {u_kms:6.0f} {kind:>7} {F_t/1e3:15.1f} {NIl/1e3:7.0f} {dphi:9.4f} {front:7.0f} {V_plan/1e3:11.2f} kV {V_prot/1e3:8.2f} kV {L:6.0f}")
    print(f"  Below 5 km/s the moving ramp front is shorter than a 25 m load coil, so the planned pull's turn voltage is dPhi' x u, hundreds of volts;")
    print(f"  above it the whole coil ramps together and the voltage is the {T_PLAN*1e3:.0f} ms figure.  The protection rise sets the insulation at every speed.")
    out["passes"] = passes
    out["trans_lift"] = dict(V2=dphi_l * L_LIFT / TAU_HW, V22=dphi_l * L_LIFT / T_PLAN, dphi=dphi_l)
    if not quick and h > 0.003:
        print()
        print("  Grid check at 2.5 mm (one solve each):")
        for name, geom, NI, wt, Ir, ref in (("support at the load cap", bo.BEST_SUPPORT, NI_l, 0.012, I_REF_12, phi_at(rows, lin, NI_l)), ("lift at its requirement", bo.BEST_LIFT, NI_lift, 0.003, I_REF_3, dphi_l)):
            t25 = Track(geom, 0.0025, w_tape=wt, I_ref=Ir, NI_ref=NI)
            t0 = time.time(); g, s = t25.solve(t25.NI_design(NI))
            print(f"    {name}: flux {s['phi'][0]:.4f} Wb/m at 2.5 mm against {ref:.4f} at {h*1e3:g} mm ({(s['phi'][0]/ref-1)*100:+.2f} percent), pull {s['F']/1e3:.1f} kN/m ({time.time()-t0:.0f} s)")
    return out


# =============================================================================
# PART 1B: series turns, parallel tapes, insulation, the bundle
# =============================================================================

def series_table(label, N_t, I_tape, V_turn2, V_turn22, E_coil):
    """Terminal voltage and current against the number of series turns s (each a cable of N_t/s tapes in parallel)."""
    print(f"  {label}: {N_t:,.0f} turns of tape per coil at {I_tape:.1f} A; one turn sees {V_turn2/1e3:.2f} kV in the 2 ms rise and {V_turn22/1e3:.2f} kV in the {T_PLAN*1e3:.0f} ms ramp;")
    print(f"  VA = 2 E / t_rise = {2*E_coil/TAU_HW/1e6:.0f} MVA at 2 ms and {2*E_coil/T_PLAN/1e6:.1f} MVA at {T_PLAN*1e3:.0f} ms whatever the connection.")
    print(f"  {'series turns s':>14} {'tapes/cable':>11} {'current A':>9} {'V 2 ms':>10} {'V 22 ms':>10} {'hold':>12} {'turn-turn':>10} {'family-family':>13} {'to ground':>10}")
    rows = []
    for s in (1, 2, 5, 10, 50, int(N_t)):
        n_par = N_t / s; I = n_par * I_tape; V2 = s * V_turn2; V22 = s * V_turn22
        tt = V_turn2 if s > 1 else 0.0
        rows.append(dict(s=s, n_par=n_par, I=I, V2=V2, V22=V22, V_tt=tt, V_ff=V2, V_g=V2))
        print(f"  {s:14,d} {n_par:11,.0f} {I:9,.0f} {fmt_v(V2):>10} {fmt_v(V22):>10} {'~0 V, ' + f'{I:,.0f} A':>12} {fmt_v(tt):>10} {fmt_v(V2):>13} {fmt_v(V2):>10}")
    print(f"  (hold: a superconducting coil at steady current drops only its joint and lead voltages; the current is the same as in the pulse.")
    print(f"  Turn-turn: neighbouring series turns differ by one turn voltage.  Family-family: the 20 coils of a bundle end at staggered places, so a tape")
    print(f"  at a coil's end lies beside a tape of a coil in its middle: up to the full terminal voltage.  To ground: the terminal voltage with one end")
    print(f"  grounded, half of it with the supply's midpoint grounded.)")
    return rows


def fmt_v(v):
    return f"{v/1e6:.2f} MV" if v >= 1e6 else (f"{v/1e3:.1f} kV" if v >= 1e3 else f"{v:.0f} V")


def bundle_layout(w_b, h_b, n_fam, A_fam, s, V_turn, V_fam, E_kVmm, n_x, n_y, extra_area=0.0, f_cable=F_CABLE, f_cool=F_COOL):
    """Height a bundle of width w_b needs for n_fam family cells (n_x by n_y) each holding s series turns (cables) of
    conductor area A_fam / s, with insulation at the design stress E (kV/mm): walls of V/E between cells, between the
    series turns in a cell, and to ground around the bundle; coolant channels as f_cool of the conductor's area; extra_area is
    conductor outside the cells (the persistent winding in a support bundle), with its own thin insulation and cooling."""
    E = E_kVmm * 1e6                                          # V/m
    t_fam = V_fam / E; t_turn = V_turn / E if s > 1 else 0.0; t_g = V_fam / E
    t_fam = max(t_fam, 50e-6); t_g = max(t_g, 50e-6)          # a film at least 50 um [assumption]
    c_w = (w_b - 2 * t_g - (n_x - 1) * t_fam) / n_x
    if c_w <= 0:
        return None
    c_h = (A_fam / f_cable) * (1 + f_cool) / c_w + (s - 1) * t_turn
    h_cells = n_y * c_h + (n_y - 1) * t_fam
    h_extra = (extra_area / f_cable) * (1 + f_cool) / (w_b - 2 * t_g) if extra_area else 0.0
    h_need = h_cells + h_extra + 2 * t_g
    return dict(t_fam=t_fam, t_turn=t_turn, t_g=t_g, c_w=c_w, c_h=c_h, h_need=h_need, fits=h_need <= h_b + 1e-9, h_b=h_b, w_b=w_b)


def solve_enlarged(label, geom_base, shape_key, w, hgt, allow, NI_req_guess, F_req, w_tape, theva, theva77, h, NI_other=None, kind="support"):
    """Register a bundle shape, solve the track at the requirement (two nonlinear solves) and give the pull, the worst tape,
    the turns and the tape mass.  For a support track the 3 mm load winding's worst tape is taken in the design state."""
    bo.SHAPES[shape_key] = dict(kind="rect", w=w, h=hgt)
    geom = dict(geom_base); geom["shape"] = shape_key; geom["allow"] = allow
    I_ref = I_REF_12 if w_tape == 0.012 else I_REF_3
    tr = Track(geom, h, w_tape=w_tape, I_ref=I_ref, NI_ref=NI_req_guess)
    g, s1 = tr.solve(tr.NI_design(NI_req_guess))
    NI2 = NI_req_guess * math.sqrt(F_req / s1["F"])
    g, s2 = tr.solve(tr.NI_design(NI2))
    # log-log correction to the requirement
    NI_req = math.exp(math.log(NI_req_guess) + (math.log(F_req) - math.log(s1["F"])) * (math.log(NI2) - math.log(NI_req_guess)) / (math.log(s2["F"]) - math.log(s1["F"])))
    wt, lim = worst_tape(g, tr, tr.NI_design(NI2), w_tape, theva, theva77)
    sc = NI_req / NI2
    ic = float(bo.ic_angle(theva, np.array([wt[1] * sc]), np.array([wt[2]]), w_tape)[0])
    I50 = 0.5 * ic; turns = len(tr.legs) * NI_req / I50; mass = turns * w_tape * TAPE_T * TAPE_RHO
    res = dict(label=label, w=w, h=hgt, allow=allow, NI_req=NI_req, F_at_guess=s1["F"], ic=ic, B=wt[1] * sc, th=wt[2], face=wt[3], I50=I50, turns_total=turns, mass=mass, sat=s2["sat"], W=s2["W"] * sc**2)
    if kind == "support":
        NI_446 = bo.NI_of_F(r_sup_global, F_req + F_LOAD) * NI_req / r_sup_global["NI_req"]
        g3, s3 = tr.solve(tr.NI_design(NI_446))
        wt3, _ = worst_tape(g3, tr, tr.NI_design(NI_446), 0.003, theva, theva77)
        I3 = 0.5 * wt3[0]; NI_l_own = bo.NI_of_F(r_sup_global, F_LOAD) * NI_req / r_sup_global["NI_req"]
        turns3 = len(tr.legs) * NI_l_own / I3; mass3 = turns3 * 0.003 * TAPE_T * TAPE_RHO
        res.update(NI_446=NI_446, F_446=s3["F"], ic3=wt3[0], I3=I3, turns3_total=turns3, mass3=mass3, NI_l_own=NI_l_own)
    return res


r_sup_global = None


def part_1b(h, theva, theva77, r_sup, r_lift, a1, quick=False):
    global r_sup_global
    r_sup_global = r_sup
    hr("PART 1B: SERIES TURNS, PARALLEL TAPES AND THE INSULATION")
    out = {}
    # turn counts per coil
    N_lift_leg = r_lift["turns3"] / 2; N_lift_coil = N_lift_leg / N_FAM; I_lift = r_lift["I3"]
    N_load_leg = r_sup["turns3_446"] / 4; N_load_coil = N_load_leg / N_FAM; I_load = r_sup["I3_446"]
    N_sup_leg = r_sup["turns12"] / 4; I_sup = r_sup["I12"]
    E_lift = r_lift["E_coil"]; E_load = r_sup["E_coil"]
    V_l2, V_l22 = a1["trans_lift"]["V2"], a1["trans_lift"]["V22"]
    tl = a1["load_transitions"]["load winding alone, 0 to its cap (the book's basis)"]
    V_d2, V_d22 = tl["V2"], tl["V22"]
    print(f"  Turn counts: lift winding {N_lift_leg:,.0f} turns of 3 mm tape per leg ({N_lift_coil:,.0f} per 50 m coil, 20 families) at {I_lift:.1f} A; load winding")
    print(f"  {N_load_leg:,.0f} per leg ({N_load_coil:,.0f} per 25 m coil) at {I_load:.1f} A; persistent support winding {N_sup_leg:,.0f} turns of 12 mm per leg at {I_sup:.1f} A.")
    print(f"  Stored energy per coil: lift {E_lift/1e6:.3f} MJ, load {E_load/1e6:.3f} MJ (the Part 3b records, the book's 0.301 and 0.146).")
    print()
    out["lift_series"] = series_table("Lift coil, 50 m", N_lift_coil, I_lift, V_l2, V_l22, E_lift)
    print()
    out["load_series"] = series_table("Load coil, 25 m (its own field, 0 to 287 kA)", N_load_coil, I_load, V_d2, V_d22, E_load)
    tp = a1["load_transitions"]["protection pulse at the load cap on the support at 141"]
    print(f"  On the support winding's 141 kA the same 0 to 287 kA pulse changes the loop flux by {tp['dphi']:.3f} Wb/m instead of {tl['dphi']:.3f}: the turn voltage")
    print(f"  is {tp['V2']/1e3:.1f} kV at 2 ms; the insulation is set by the larger of the two.")
    # the insulation stress
    print()
    print("  Design stress.  The coils sit in nitrogen at 0.4 atm and 70 K, at its boiling point, so gas is present in any void and the gas sets the")
    print("  limit, not the film.  The numbers taken from the literature (listed in the report):")
    print("    Kapton HN film, room temperature, short-time AC: 154 kV/mm at 125 um to 303 kV/mm at 25 um (DuPont data sheet).")
    print("    Kapton in liquid nitrogen, AC: the strength falls with the sheet thickness and levels off near 91 kV/mm (Nomex and Kapton tests in LN2).")
    print("    Liquid nitrogen without bubbles: 12 to 50 kV/mm at 200 to 500 kPa and 77 K; bubbles are field concentrations and set the strength (bubble reviews).")
    print("    Nitrogen gas: Paschen minimum 251 V at 0.67 Torr cm; at 0.4 atm and 70 K the gas density is that of 1.7 atm at room temperature,")
    print("    so a millimetre gap holds about 1.7 x 3 kV/mm, 5 kV/mm, and a 0.1 mm void more per millimetre but only about 1.2 kV in all.")
    print("    HTS cable practice: the LN2 / laminated-paper insulation is designed partial-discharge-free (Kim, Kwag et al. 2007); design stresses of")
    print("    the order of 10 kV/mm are used [UNVERIFIED: the figure is not read from the paper].")
    print(f"  Three design stresses are carried: {', '.join(E_DESIGN)} kV/mm.  The 5 kV/mm case assumes gas-filled voids of millimetre size in")
    print("  the insulation stack at 0.4 atm; 10 kV/mm a bubble-free impregnated multilayer; 20 kV/mm a Kapton film at a quarter of its LN2 strength with")
    print("  no gas in series with it.  Film thickness is never taken below 50 um.")
    # bundle sizing
    print()
    print(f"  Bundle cross-section.  Each family coil is a cell: its s series turns are cables of N_t / s tapes stacked at a fill of {F_CABLE:.2f} (tape thickness over")
    print(f"  pitch), coolant channels take {F_COOL:.0%} of the conductor area, walls of V / E between cells, between series turns and around the bundle.  The lift")
    print(f"  bundle (200 x 50 mm) has its 20 cells in 10 columns of 2; the support bundle (150 x 67 mm) holds the persistent winding's 12 mm tape (thin")
    print(f"  insulation: a closed persistent loop has no terminal voltage) below 20 load-winding cells in 5 columns of 4.")
    A_lift_fam = N_lift_coil * 0.003 * TAPE_T
    A_load_fam = N_load_coil * 0.003 * TAPE_T
    A_sup_bundle = N_sup_leg * 0.012 * TAPE_T                  # the persistent winding of one bundle
    print(f"  Conductor per family: lift {A_lift_fam*1e6:.0f} mm2 (bundle {20*A_lift_fam*1e6:,.0f} of {200*50:,} mm2, {20*A_lift_fam/0.01:.0%}); load {A_load_fam*1e6:.0f} mm2 (20 families {20*A_load_fam*1e6:,.0f} mm2) plus the")
    print(f"  support winding {A_sup_bundle*1e6:,.0f} mm2 in the 150 x 67 = {150*67:,} mm2 bundle ({(20*A_load_fam+A_sup_bundle)/(0.15*0.0667):.0%}).")
    lay = {}
    for tag, (w_b, h_b, n_x, n_y, A_fam, V2, extra) in {"lift": (0.20, 0.05, 10, 2, A_lift_fam, V_l2, 0.0), "load": (0.15, 0.0667, 5, 4, A_load_fam, max(V_d2, tp["V2"]), A_sup_bundle)}.items():
        print()
        print(f"  {tag} bundle, {w_b*1e3:.0f} x {h_b*1e3:.0f} mm envelope:")
        print(f"  {'E kV/mm':>24} {'s':>5} {'t turn mm':>9} {'t family mm':>11} {'t ground mm':>11} {'cell mm':>12} {'height needed mm':>16} {'fits':>5}")
        for ek, E in E_DESIGN.items():
            for s in (1, 2, 5, 10):
                L = bundle_layout(w_b, h_b, N_FAM, A_fam, s, V2, s * V2, E, n_x, n_y, extra_area=extra)
                if L is None:
                    print(f"  {ek:>24} {s:5d} {'-':>9} {'-':>11} {'-':>11} {'no room for the walls':>12}"); continue
                lay[(tag, ek, s)] = L
                print(f"  {ek:>24} {s:5d} {L['t_turn']*1e3:9.2f} {L['t_fam']*1e3:11.2f} {L['t_g']*1e3:11.2f} {L['c_w']*1e3:5.1f} x {L['c_h']*1e3:5.1f} {L['h_need']*1e3:16.1f} {'yes' if L['fits'] else 'NO':>5}")
        L0 = bundle_layout(w_b, h_b, N_FAM, A_fam, 1, 0.0, 0.0, 1e9, n_x, n_y, extra_area=extra)
        print(f"  (with no insulation at all, the cells and their coolant need {L0['h_need']*1e3:.1f} mm of the {h_b*1e3:.0f} mm height)")
    out["layouts"] = {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in lay.items()}
    # enlarged bundles solved
    print()
    print(f"  The enlarged bundles solved ({h*1e3:g} mm cells, saturating plate): the width kept, the height the layout needs, the 10 mm tape-to-gap allowance")
    print(f"  grown by the ground wall beyond the 5 mm insulation face; the pull at the present ampere-turns, the ampere-turns for the requirement, the worst tape")
    print(f"  (full angle data), the turns at 50 percent of its limit, and the tape mass per track (Equation 7.3: turns x width x thickness x density).")
    runs = []
    base_cases = [("lift", "design 200 x 50, wall 10", 0.20, 0.05, 0.010)]
    for ek, E in E_DESIGN.items():
        for s in (1, 5):
            L = lay.get(("lift", ek, s))
            if L and not L["fits"]:
                base_cases.append(("lift", f"E {ek.split()[0]} kV/mm, s = {s}: 200 x {L['h_need']*1e3:.0f}", 0.20, round(L["h_need"], 4), 0.010 + max(0.0, L["t_g"] - 0.005)))
    sup_cases = [("support", "design 150 x 67, wall 10", 0.15, 0.0667, 0.010)]
    for ek, E in E_DESIGN.items():
        for s in (1, 5):
            L = lay.get(("load", ek, s))
            if L and not L["fits"]:
                sup_cases.append(("support", f"E {ek.split()[0]} kV/mm, s = {s}: 150 x {L['h_need']*1e3:.0f}", 0.15, round(L["h_need"], 4), 0.010 + max(0.0, L["t_g"] - 0.005)))
    if quick:
        base_cases = base_cases[:2]; sup_cases = sup_cases[:2]
    seen = set()
    print(f"  {'track':>7} {'bundle':<36} {'allow mm':>8} {'F at design NI':>14} {'NI req kA':>9} {'worst tape':>22} {'I50 A':>6} {'turns/leg':>9} {'tape kg/m':>9} {'sat %':>5} {'load tape kg/m':>14}")
    for kind, label, w, hgt, allow in base_cases + sup_cases:
        key = (kind, round(w, 4), round(hgt, 4), round(allow, 4))
        if key in seen:
            continue
        seen.add(key)
        if kind == "lift":
            res = solve_enlarged(label, bo.BEST_LIFT, f"ins {w*1e3:.0f}x{hgt*1e3:.0f}", w, hgt, allow, r_lift["NI_req"], F_LIFT, 0.003, theva, theva77, h, kind="lift")
            print(f"  {kind:>7} {label:<36} {allow*1e3:8.0f} {res['F_at_guess']/1e3:14.1f} {res['NI_req']/1e3:9.0f} {res['B']:5.2f} T {res['th']:4.0f} deg {res['ic']:5.1f} A {res['I50']:6.1f} {res['turns_total']/2:9,.0f} {res['mass']:9.1f} {res['sat']*100:5.1f} {'-':>14}")
        else:
            res = solve_enlarged(label, bo.BEST_SUPPORT, f"ins {w*1e3:.0f}x{hgt*1e3:.0f}", w, hgt, allow, r_sup["NI_req"], r_sup["F_req"], 0.012, theva, theva77, h, kind="support")
            print(f"  {kind:>7} {label:<36} {allow*1e3:8.0f} {res['F_at_guess']/1e3:14.1f} {res['NI_req']/1e3:9.0f} {res['B']:5.2f} T {res['th']:4.0f} deg {res['ic']:5.1f} A {res['I50']:6.1f} {res['turns_total']/4:9,.0f} {res['mass']:9.1f} {res['sat']*100:5.1f} {res['mass3']:14.1f}")
        res["kind"] = kind; runs.append(res)
    print(f"  (the design rows are this script's own solve of the Part 3b geometries, to be read against the Part 3b records: support {r_sup['NI_req']/1e3:.0f} kA, {r_sup['ic12']:.0f} A,")
    print(f"  {r_sup['tape12']:.1f} kg/m + {r_sup['tape3_446']:.1f} load; lift {r_lift['NI_req']/1e3:.0f} kA, {r_lift['ic3']:.0f} A, {r_lift['tape3']:.1f} kg/m; the grid here is coarser)")
    out["enlarged"] = runs
    return out


# =============================================================================
# PART 1C: what a few-turn, high-current winding brings with it
# =============================================================================

def strip_lattice(w_b, h_b, xc, y_top, w_tape, N):
    """N tapes of width w_tape filling a w_b x h_b block whose top is at y_top: columns across, rows up (as tape_centers)."""
    xs, ys, n_col, n_row, pitch = bo.tape_centers(xc - w_b / 2, xc + w_b / 2, y_top - h_b, y_top, w_tape, N)
    return xs[:N], ys[:N]


def cable_inductance(xs, ys, w_tape, x_ret, y_plate, mu_plate=cf.PLATE_MU, n_sub=6):
    """Per-meter inductance matrix of parallel flat tapes at (xs, ys), width w_tape along x, with their returns mirrored at
    x -> 2 x_ret - x (the other bundle, opposite current), in free space with the plate as an image (strength (mu - 1) / (mu + 1))
    at y -> 2 y_plate - y.  Each tape is n_sub filaments across its width and the tape-to-tape term is the mean over the
    filament pairs (tapes stacked closer than their width need this; a centre-to-centre log would exceed the self term); a
    filament's own term uses its geometric mean distance 0.2235 x its width [textbook value for a thin strip]."""
    n = len(xs)
    X = np.asarray(xs, float); Y = np.asarray(ys, float)
    off = ((np.arange(n_sub) + 0.5) / n_sub - 0.5) * w_tape
    FX = (X[:, None] + off[None, :]).ravel(); FY = np.repeat(Y, n_sub)
    FXr = 2 * x_ret - FX; FYi = 2 * y_plate - FY
    k = (mu_plate - 1) / (mu_plate + 1)
    gmd = 0.2235 * w_tape / n_sub
    def lnd(x1, y1, x2, y2, self_term=False):
        d = np.hypot(x1[:, None] - x2[None, :], y1[:, None] - y2[None, :])
        if self_term:
            np.fill_diagonal(d, gmd)
        return np.log(np.maximum(d, 1e-12))
    # the flux between go i and return i per unit current in loop j (go j, return j), and the plate's images (same sign for go, opposite for return)
    M = lnd(FX, FY, FXr, FY) - lnd(FX, FY, FX, FY, True)
    M += lnd(FXr, FY, FX, FY) - lnd(FXr, FY, FXr, FY, True)
    M += k * (lnd(FX, FY, FXr, FYi) - lnd(FX, FY, FX, FYi))
    M += k * (lnd(FXr, FY, FX, FYi) - lnd(FXr, FY, FXr, FYi))
    M *= MU0 / (2 * math.pi)
    return M.reshape(n, n_sub, n, n_sub).mean(axis=(1, 3))


def current_sharing(L):
    """Parallel tapes fed by one voltage: L dI/dt = V 1 -> the current division is L^-1 1, normalized to a mean of one."""
    x = np.linalg.solve(L, np.ones(L.shape[0]))
    return x / x.mean()


def lead_heat(I, T_h=T_WARM, T_c=T_COLD, k_cu=400.0, A_over_l=None, n=200):
    """One conduction-cooled copper lead, Wiedemann-Franz (rho = L0 T / k, k constant [UNVERIFIED simplification]): the heat into the
    cold end for a given A/l at current I, and the optimum A/l when none is given.  Returns (Q_cold W, A_over_l m)."""
    def Q_cold(aol):
        # steady 1-D: k A T'' + I^2 rho / A = 0 with rho = L0 T / k; nondimensional in s = x/l: T'' = -(I^2 L0 T) / (k^2 (A/l)^2)
        # linear ODE T'' + c T = 0, c = (I L0^0.5 / (k A/l))^2 ; T(0) = T_h, T(1) = T_c ; Q_cold = -k (A/l) T'(1)
        c = (I * math.sqrt(L0) / (k_cu * aol))**2
        w = math.sqrt(c)
        if w < 1e-9:
            return k_cu * aol * (T_h - T_c)
        # T = a cos(w s) + b sin(w s); T(0) = T_h -> a = T_h; T(1) = T_c -> b = (T_c - T_h cos w) / sin w
        if abs(math.sin(w)) < 1e-9 or w >= math.pi:
            return float("inf")
        b = (T_c - T_h * math.cos(w)) / math.sin(w)
        dT1 = -T_h * w * math.sin(w) + b * w * math.cos(w)
        return -k_cu * aol * dT1
    if A_over_l is None:
        lo, hi = 1e-9, 1e-1
        best = None
        for aol in np.logspace(math.log10(lo), math.log10(hi), 400):
            q = Q_cold(aol)
            if best is None or q < best[0]:
                best = (q, aol)
        return best
    return Q_cold(A_over_l), A_over_l


def part_1c(h, theva, theva77, r_sup, r_lift, a1, b1, quick=False):
    hr("PART 1C: WHAT A FEW-TURN, HIGH-CURRENT WINDING BRINGS WITH IT")
    out = {}
    # C.1 current sharing in an untransposed cable
    print("  C.1  Current sharing among parallel tapes during the 2 ms rise.  The tapes of one family coil form one cable (s = 1); fed by one voltage,")
    print("  zero-resistance strands divide the current by their inductance matrix alone (the common loop flux cancels out of the division), and with no")
    print("  resistance the division set during the rise stays.  The matrix is the 2-D free-space one of thin strips at the lattice positions with their")
    print("  returns in the other bundle and the plate's image.  The family's tapes are taken as one block of the bundle (a design choice; a block is the")
    print("  worst case for an untransposed cable because its tapes span the bundle's field gradient).")
    N_lift_coil = r_lift["turns3"] / 2 / N_FAM; N_load_coil = r_sup["turns3_446"] / 4 / N_FAM
    for tag, geom, N_t, w_b, h_b, n_x, n_y, xl, xr in (("lift", bo.BEST_LIFT, N_lift_coil, 0.20, 0.05, 10, 2, -0.4, 0.4), ("load", bo.BEST_SUPPORT, N_load_coil, 0.15, 0.0667, 5, 4, -0.575, -0.075)):
        allow = geom["allow"]; y_plate = GAP
        c_w = w_b / n_x; c_h = h_b / n_y
        # block at the bundle's upper outer corner (nearest the plate and the loop's inside)
        xs, ys = strip_lattice(c_w, c_h, xl + w_b / 2 - c_w / 2, -allow, 0.003, int(N_t))
        if quick:
            xs, ys = xs[::4], ys[::4]
        L = cable_inductance(xs, ys, 0.003, (xl + xr) / 2, y_plate)
        f = current_sharing(L)
        # a Roebel-like transposed cable: every strand sees the same mean position -> uniform by construction
        print(f"    {tag} coil cable, {len(xs)} tapes in a {c_w*1e3:.0f} x {c_h*1e3:.1f} mm block: tape current from {f.min():.2f} to {f.max():.2f} of the mean (untransposed);")
        print(f"      {np.mean(f > 1.5):.0%} of the tapes carry more than 1.5 x the mean, {np.mean(f > 2.0):.0%} more than twice; the heaviest tape would reach the design current's")
        print(f"      share x {f.max():.2f} of the worst tape's limit: {0.5*f.max():.0%} at the 50 percent design (the limit is 100 percent).  A transposed cable (Roebel,")
        print(f"      twisted stack, CORC) makes every strand see the same inductance and carries the mean in each.")
        out[f"sharing_{tag}"] = dict(n=int(len(xs)), fmin=float(f.min()), fmax=float(f.max()), over15=float(np.mean(f > 1.5)), over2=float(np.mean(f > 2.0)))
        # the simplest check: two tapes side by side
        L2 = cable_inductance(np.array([xs[0], xs[-1]]), np.array([ys[0], ys[-1]]), 0.003, (xl + xr) / 2, y_plate)
        f2 = current_sharing(L2)
        print(f"      (two tapes at the block's extreme positions alone: {f2[0]:.2f} and {f2[1]:.2f} of their mean)")
    print("  Transposition: Roebel strands are punched from the tape into a meander about half the tape's width wide and assembled into a fully transposed")
    print("  stack (the width lost is tape lost; in-plane bending is limited) [Goldacker et al. 2014 review]; CORC winds 2 to 4 mm tapes helically on a round")
    print("  former a few millimetres across, well inside the 40 mm double-bend diameter given for this tape, so a CORC on a 40 mm former would be a very")
    print("  open cable [UNVERIFIED: fill factors of Roebel about 0.5 and CORC about 0.3 of a plain stack are not from a source read here]; a twisted stacked")
    print("  tape keeps the stack's fill and transposes by a slow twist (pitch of the order of 100 mm for a 12 mm stack), which is the one that fits the")
    print("  40 mm bend limit and the flat bundle.  Its cost is the twist's torsion strain on the tape, not fill [UNVERIFIED].")
    print("  Loss: the hysteresis loss of each tape in the pulse is set by the field change through its face, not by the cabling, so the strip estimate")
    print("  (605 GJ per launch as the upper bound, Section 7.7) is unchanged by transposition; what transposition removes is the current imbalance above,")
    print("  which in an untransposed cable drives the heaviest tapes to their limit on every pulse (flux-flow heat and the risk of a quench), and the")
    print("  coupling currents between strands, which in superconducting strands do not decay and so do not dissipate unless the strands touch.")
    # C.2 leads
    print()
    print(f"  C.2  Current leads from {T_WARM:.0f} K to {T_COLD:.0f} K.  A conduction-cooled copper lead optimized for a steady current brings in sqrt(L0 (T_h^2 - T_c^2)) per ampere")
    q_opt = math.sqrt(L0 * (T_WARM**2 - T_COLD**2))
    print(f"  = {q_opt*1e3:.1f} W/kA (Wiedemann-Franz, no vapour cooling; the 1-D lead model below reproduces it).  Pulsed duty: a lead sized only for the pulse's")
    print(f"  adiabatic heating (copper at {RHO_CU_77*1e9:.1f} nohm m, 20 K rise allowed in a 0.1 s pulse [assumptions]) has a far smaller cross-section and conduction leak,")
    print(f"  but cannot hold a standing load.  The load winding holds standing loads with no time limit, so its leads are the steady kind unless the")
    print(f"  electronics are cold.")
    q1, aol1 = lead_heat(1000.0)
    print(f"  1-D lead model check at 1 kA: {q1:.1f} W per lead at the optimum A/l = {aol1*1e6:.0f} mm2 per m; leak with no current through the same lead {lead_heat(0.0, A_over_l=aol1)[0]:.1f} W.")
    rows = []
    N_lift_coil = r_lift["turns3"] / 2 / N_FAM; N_load_coil = r_sup["turns3_446"] / 4 / N_FAM
    for tag, N_t, I_t, L_coil in (("lift", N_lift_coil, r_lift["I3"], L_LIFT), ("load", N_load_coil, r_sup["I3_446"], L_LOAD)):
        for s in (1, 5, 50, int(N_t)):
            I = N_t * I_t / s
            q_dc, aol = lead_heat(I)
            # a pulse-sized lead: adiabatic heating J^2 rho t = rho_m c dT -> J = sqrt(rho_m c dT / (rho t)); copper 8960 kg/m3, c 200 J/kg/K at 77 K [UNVERIFIED]
            J = math.sqrt(8960 * 200 * 20 / (RHO_CU_77 * 0.1)); A_p = I / J
            q_p0 = lead_heat(0.0, A_over_l=A_p / 1.0)[0]                       # 1 m lead, no current
            coils_site = 500.0 / L_coil * N_FAM                                 # one 500 m site, one track
            coils_parked = (L_SC + L_coil) / L_coil * N_FAM                     # the families under a 100 m parked load, one track
            rows.append((tag, s, I, 2 * q_dc, 2 * q_p0, coils_site, coils_parked))
    print(f"  {'coil':>5} {'s':>6} {'lead A':>8} {'steady pair W':>13} {'pulse-sized pair W (I = 0)':>26} {'per 500 m site, one track kW':>28} {'under a parked 5,000 t, both tracks kW':>38}")
    for tag, s, I, q_dc, q_p, cs, cp in rows:
        print(f"  {tag:>5} {s:6,d} {I:8,.0f} {q_dc:13,.0f} {q_p:26.2f} {q_dc*cs/1e3:28.1f} {(q_dc*cp*2/1e3 if tag == 'load' else 0.0):38.1f}")
    print(f"  (per site: the coils of one track whose leads reach one 500 m site, {500/L_LIFT*N_FAM:.0f} lift or {500/L_LOAD*N_FAM:.0f} load coils; the parked column is the load")
    print(f"  coils energized under a 100 m spacecraft on both tracks, {(L_SC+L_LOAD)/L_LOAD*N_FAM*2:.0f} coils, with steady leads.  At {N_FAM*ring.L_RING/L_LOAD/1e6:.1f} million load coils per")
    print(f"  track every steady lead pair is a continuous heat load whether the coil is on or not: the few-turn design with steady leads is {rows[4][3]*N_FAM*ring.L_RING/L_LOAD*2/1e9:.1f} GW of 70 K")
    print(f"  heat per ring for the load windings alone, and the all-series one {rows[7][3]*N_FAM*ring.L_RING/L_LOAD*2/1e6:.0f} MW, {rows[7][3]*N_FAM*ring.L_RING/L_LOAD*2/1e6*94/1e3:.1f} GW of electricity at 94 J per J.)")
    print(f"  Alternatives: cold power electronics (the converter at 70 K: its own loss, of the order of a few percent of the pulse energy per pulse")
    print(f"  [UNVERIFIED], lands in the cryostat but only when the coil pulses; the dc link comes down as a high-voltage, low-current feed whose leads are")
    print(f"  the all-series kind); a cold energy store (a superconducting store beside the coil, charged slowly through a small lead and switched into the")
    print(f"  coil by a cold converter, no high-current lead at all).  Either keeps the warm-to-cold path at the all-series lead's size or below.")
    out["leads"] = rows
    # C.3 supplies
    print()
    print(f"  C.3  Supplies.  One supply per coil: {N_FAM*ring.L_RING/L_LIFT/1e6:.1f} million per lift track ({L_LIFT:.0f} m coils), {N_FAM*ring.L_RING/L_LOAD/1e6:.1f} million per load winding ({L_LOAD:.0f} m).  The rating per")
    print(f"  supply is 2 E / t_rise whatever s: lift {2*r_lift['E_coil']/TAU_HW/1e6:.0f} MVA at 2 ms ({2*r_lift['E_coil']/T_PLAN/1e6:.1f} at {T_PLAN*1e3:.0f} ms), load {2*r_sup['E_coil']/TAU_HW/1e6:.0f} MVA ({2*r_sup['E_coil']/T_PLAN/1e6:.1f}).  The choice of s moves the")
    print(f"  rating between voltage and current (Part 1b's table): s = 1 is a {a1['trans_lift']['V2']/1e3:.0f} kV, {N_lift_coil*r_lift['I3']/1e3:.0f} kA pulse converter for a lift coil, the all-series case a")
    print(f"  {N_lift_coil*a1['trans_lift']['V2']/1e6:.0f} MV, {r_lift['I3']:.0f} A one; s of 5 to 10 is the range of pulsed-power practice (tens of kilovolts, kiloamperes) [UNVERIFIED as a design rule].")
    print(f"  Only the coils under a pass are on: at most about {(L_SC + 0.05*30000 + L_LIFT)/L_LIFT*N_FAM:.0f} lift coils per track at 30 km/s and {(L_SC+0.1*1000+L_LOAD)/L_LOAD*N_FAM:.0f} load coils at 1 km/s.")
    # C.4 protection rise and coil length: the sweep runs separately
    print()
    print("  C.4  The two requirements that set the voltage, the 2 ms protection rise and the coil length, are swept through final_bearing's cases by")
    print("  'python bearing_electrical.py --part=sweep' (slow); its summary is in the report.  From the flux: the turn voltage scales with the coil length")
    print("  and inversely with the rise time, so a 12.5 m lift coil at 2 ms sees a quarter of the 50 m coil's turn voltage and a 50 m coil at 8 ms the same.")
    for L in (12.5, 25.0, 50.0):
        E = r_lift["W_req"] * L / N_FAM
        print(f"    lift coil {L:5.1f} m: {E/1e6:.3f} MJ, {2*E/TAU_HW/1e6:4.0f} MVA at 2 ms, turn voltage {a1['trans_lift']['dphi']*L/TAU_HW/1e3:5.1f} kV at 2 ms, {N_FAM*ring.L_RING/L/1e6:5.1f} million supplies per track")
    return out


def sweep_cases(configs=None):
    """Part 1c: final_bearing's Part B cases (the 87 of bearing_design.part_F plus the added ones) with the protection rise and the
    lift coil length changed through the module constants the way final_bearing itself points the dynamics modules at its masses."""
    hr("PART 1C SWEEP: THE PROTECTION RISE AND THE LIFT COIL LENGTH THROUGH FINAL_BEARING'S CASES")
    import sled_passage_simulation as sp
    theva = bd.read_theva(); sup, lift = fb.load_air()
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    with contextlib.redirect_stdout(io.StringIO()):
        A = fb.part_A(theva, sup, lift)
    tau0, L0_ = bd.TAU_HW, fb.L_LIFT
    configs = configs or [(2e-3, 50.0), (4e-3, 50.0), (6e-3, 50.0), (10e-3, 50.0), (20e-3, 50.0), (2e-3, 25.0), (2e-3, 12.5)]
    print(f"  {'rise ms':>7} {'lift coil m':>11} {'1 km/s pair mm':>14} {'30 km/s pair mm':>15} {'30 km/s one mm':>14} {'smallest single-failure gap mm':>30} {'in':<48} {'contacts':<60} {'s':>5}")
    for tau, L in configs:
        bd.TAU_HW = tau; fb.TAU_HW = tau; fb.L_LIFT = L
        buf = io.StringIO(); t0 = time.time()
        with contextlib.redirect_stdout(buf):
            fb.part_B(A, theva, sup, lift, tabs)
        txt = buf.getvalue()
        gap = case = contacts = "?"; e1 = e30 = o30 = "?"
        for line in txt.splitlines():
            if line.strip().startswith("Smallest gap in any single failure:"):
                gap = line.split(":")[1].split("mm")[0].strip(); case = line.split("'")[1] if "'" in line else "?"
            if "cases in all.  Contact in:" in line:
                contacts = line.split("Contact in:")[1].strip().rstrip(".")
            parts = line.split()
            if len(parts) >= 6 and parts[0] in ("1", "30") and line.startswith("   ") and "sled" not in line:
                try:
                    float(parts[1]); float(parts[2]); float(parts[3])
                    if parts[0] == "1" and e1 == "?": e1 = parts[2]
                    if parts[0] == "30" and e30 == "?": e30, o30 = parts[2], parts[3]
                except ValueError:
                    pass
        print(f"  {tau*1e3:7.0f} {L:11.1f} {e1:>14} {e30:>15} {o30:>14} {gap:>30} {case:<48} {contacts:<60} {time.time()-t0:5.0f}", flush=True)
        with open(os.path.join(ring.out_dir(), f"bearing_electrical_sweep_tau{tau*1e3:g}_L{L:g}.txt"), "w", encoding="utf-8") as f:
            f.write(txt)
    bd.TAU_HW = tau0; fb.TAU_HW = tau0; fb.L_LIFT = L0_
    print("  (the full case tables of every configuration are in bearing_electrical_sweep_tau<ms>_L<m>.txt; the rise applies to the planned pull's current")
    print("  lag and to the protection alike, as final_bearing uses one hardware rise for both)")


# =============================================================================
# PART 1D: the driven load winding and the persistent support winding share their loops
# =============================================================================

def inductance_matrix(track, linear=True, NI_op=None):
    """Per-meter inductance matrix of the loops, L_km = flux per meter through loop k per ampere-turn in loop m (own senses):
    linear (unsaturated plate) by exciting each loop alone; or incremental about a symmetric operating point NI_op
    (two nonlinear solves per loop, the loop's current moved by 2 percent)."""
    n = len(track.loops)
    L = np.zeros((n, n))
    if linear:
        for m in range(n):
            g, s = track.solve(track.NI_loop(m, 1e3), linear=True)
            L[:, m] = np.array(s["phi"]) / 1e3
        return L
    for m in range(n):
        base = np.array(track.NI_design(NI_op), float)
        d = 0.02 * NI_op
        up = base.copy(); dn = base.copy()
        for i, sgn in track.loops[m]:
            up[i] += sgn * d; dn[i] -= sgn * d
        g1, s1 = track.solve(list(up)); g2, s2 = track.solve(list(dn))
        L[:, m] = (np.array(s1["phi"]) - np.array(s2["phi"])) / (2 * d)
    return L


class CoilSystem:
    """The persistent support coils of a track as flux-conserving loops along the ring.  n_fam staggered coils per bundle of
    length L_c, each with n_turns turns; both loops of the track, with the per-meter inductance matrix L (2 x 2, own senses).
    The ring is a segment of n_span coil lengths with cells of length dx; the coils are those that start in the segment."""

    def __init__(self, L_c, L_mat, n_fam=N_SUP_COILS, n_span=16, cells_per_coil=200, dx_max=12.5):
        self.L_c = L_c; self.L = np.asarray(L_mat, float); self.n_fam = n_fam
        self.dx = min(L_c / cells_per_coil, dx_max)               # cells fine enough for a 125 m stretch on 10 km coils
        self.x = np.arange(0.0, n_span * L_c, self.dx) + self.dx / 2
        self.n_cells = len(self.x)
        stagger = L_c / n_fam
        starts = np.arange(-L_c + stagger, n_span * L_c, stagger)
        self.starts = starts
        self.n_coils = len(starts)
        # coverage matrix C[j, i] = 1 where coil j covers cell i (same for both loops)
        self.C = ((self.x[None, :] >= starts[:, None]) & (self.x[None, :] < starts[:, None] + L_c)).astype(float)
        self.cover_count = self.C.sum(axis=0)                 # coils covering each cell (n_fam inside)
        # overlap lengths between coils: O[j, i] = integral over coil j of coil i's coverage
        self.O = self.C @ self.C.T * self.dx
        self.cells_per_coil = cells_per_coil

    def solve(self, dNI_l):
        """dNI_l: (2, n_cells) change of the driven ampere-turns in each loop along the ring.  Returns the change of the
        support ampere-turns (2, n_cells) and the coil current changes (2, n_coils) in ampere-turns per coil."""
        nL = 2; nc = self.n_coils
        # unknowns delta[k, j] (loop k, coil j); equation for (k, j): sum_m L[k, m] (sum_i delta[m, i] O[j, i] + int_coil_j dNI_l[m] dx) = 0
        M = np.zeros((nL * nc, nL * nc)); b = np.zeros(nL * nc)
        rhs_int = self.C @ dNI_l.T * self.dx                     # (n_coils, 2): integral of dNI_l over each coil
        for k in range(nL):
            for m in range(nL):
                M[k * nc:(k + 1) * nc, m * nc:(m + 1) * nc] = self.L[k, m] * self.O
                b[k * nc:(k + 1) * nc] -= self.L[k, m] * rhs_int[:, m]
        # Overlapping staggered coils have circulation modes that change no ampere-turns anywhere (and so no flux): M is singular
        # and those modes are not excited by a flux change, so the least-squares (minimum-norm) solution is the physical one.
        delta = np.linalg.lstsq(M, b, rcond=1e-10)[0].reshape(nL, nc)
        dNI_s = delta @ self.C                                     # (2, n_cells)
        return dNI_s, delta

    def stretch(self, center, length, NI, loops=(0, 1)):
        d = np.zeros((2, self.n_cells))
        m = (self.x >= center - length / 2) & (self.x < center + length / 2)
        for k in loops:
            d[k, m] = NI
        return d


def gap_sweep(geom, h, gaps, NI_lin=1e3, NI_nl=None):
    """Linear solves at each gap: the loop flux per ampere-turn and the pull per ampere-turn squared; nonlinear at NI_nl when given."""
    rows = []
    for gp in gaps:
        tr = Track(geom, h, gap=gp)
        g, s = tr.solve(tr.NI_design(NI_lin), linear=True)
        row = dict(gap=gp, Lp=s["phi"][0] / NI_lin, Fc=s["F"] / NI_lin**2)
        if NI_nl:
            g2, s2 = tr.solve(tr.NI_design(NI_nl))
            row.update(F_nl=s2["F"], phi_nl=s2["phi"][0])
        rows.append(row)
    return rows


def stiffness(rows, g0, NI):
    """At constant current and at constant flux from the gap sweep (central differences at g0).  Pull positive toward the
    plate; dF/dg negative = the pull falls as the gap opens (unstable)."""
    gs = np.array([r["gap"] for r in rows]); Lp = np.array([r["Lp"] for r in rows]); Fc = np.array([r["Fc"] for r in rows])
    F_I = Fc * NI**2
    phi0 = float(np.interp(g0, gs, Lp)) * NI
    F_phi = Fc * (phi0 / Lp)**2
    k_I = float(np.gradient(F_I, gs)[np.argmin(np.abs(gs - g0))]); k_phi = float(np.gradient(F_phi, gs)[np.argmin(np.abs(gs - g0))])
    return k_I, k_phi, F_I, F_phi


def passage_point(u, L_coil, F_pair, kind, kernel=None, L_c=None, coupled=False, prot=True, t_settle=2.0, support_scale=None,
                  gap_dep=None, cap_up=None, cap_down=None, prot_factor=1.0, r_sup=None, loop_hz=1.0):
    """Point model of a sled passage: one meter of casing (m_c) and cable (m_k) free, the support tracks at the casing weight,
    the spread pull of the driven pair over T = max(50 ms, L_sc/u) in 20 staggered coils of length L_coil with the 20 ms command
    ramp and the 2 ms lag (sled_passage_pulsed's point model, the impulse matched so that casing and cable end up moving together),
    the 10 mm protection with a 1 ms delay (proportional over the next 10 mm, 2 ms lag).  Coupled: the support ampere-turns at
    the point follow the load winding's moving footprint through the coil system's kernel (dNI_s per unit NI m at offset xi),
    and the protection's own pull is cut by prot_factor (its stretch over the coil length).  Returns the extremes."""
    m_c = ring.M_LOAD_M; m_k = ring.M_CABLE_M; W = m_c * ring.G_NET
    match = m_k / (m_c + m_k)
    dt = 1e-4
    T = max(0.05, L_SC / u)
    qz = float(pz.q_zone(u, L_SC)); qb = float(pz.q_bare(u))
    tz0, tz1 = 0.0, L_SC / u
    t0 = -(3 * TAU_HW + 0.05 + (L_c / u if (coupled and L_c) else 0.0) + 0.2)
    t1 = pz.L_SLED / u + T + t_settle + (L_c / u if (coupled and L_c) else 0.0)
    t = np.arange(t0, t1, dt); nt = len(t)
    cov_zone = np.clip((np.minimum(tz1, t + dt) - np.maximum(tz0, t)) / dt, 0, 1)
    cov_sled = np.clip((np.minimum(pz.L_SLED / u, t + dt) - np.maximum(0.0, t)) / dt, 0, 1)
    q = qz * cov_zone + qb * (cov_sled - cov_zone)
    # the planned command: the zone's impulse spread over T, the bare sled's load after it, staggered coils, ramp, lag
    tc = 0.5 * (tz0 + tz1)
    boxes = [(tc - T / 2, tc + T / 2, match * qz * (L_SC / u) / T)]
    z1 = tc + T / 2; tb1 = pz.L_SLED / u
    if z1 < tb1:
        boxes.append((z1, tb1, match * qb))
    shaped = pz.shape_boxes(boxes, TAU_HW)
    sub = []
    for (a, b, lev) in shaped:
        for k in range(N_FAM):
            d_k = (k + 0.5) * L_coil / N_FAM
            sub.append((a - d_k / u, b + (L_coil - d_k) / u, lev * (b - a) / ((b - a) + L_coil / u) / N_FAM))
    cmd = np.zeros(nt)
    for (a, b, lev) in sub:
        cmd += lev * np.clip((np.minimum(b, t + dt) - np.maximum(a, t)) / dt, 0, 1)
    k_r = int(round(RAMP / dt)); cmd = np.convolve(cmd, np.ones(k_r) / k_r, mode="same")
    F_pair_cap = F_pair
    cmd = np.clip(cmd, -F_pair_cap, F_pair_cap) if kind == "lift" else np.clip(cmd, 0.0, F_pair_cap)
    # the coupled support deficit: dNI_s(t) = sum_xi kernel(xi) NI_l(t - xi/u) dxi ; NI_l from the command through the pull curve
    NI_s0 = bo.NI_of_F(r_sup, W_HALF)
    dNI = np.zeros(nt)
    if coupled and kernel is not None and kind == "load":
        xi, R = kernel                                              # R: ampere-turns of support change per (ampere-turn x m) of load at offset xi
        NIl_t = np.array([max(bo.NI_of_F(r_sup, W_HALF + c / 2) - NI_s0, 0.0) if c > 0 else 0.0 for c in cmd])
        # convolution in time: c(t - xi/u) -> shift by xi/u
        dxi = xi[1] - xi[0]
        for x_off, Rv in zip(xi, R):
            sh = int(round(x_off / u / dt))
            if abs(sh) >= nt or Rv == 0.0:
                continue
            if sh >= 0:
                dNI[sh:] += Rv * dxi * NIl_t[:nt - sh]
            else:
                dNI[:nt + sh] += Rv * dxi * NIl_t[-sh:]
    # integrate, as sled_passage_pulsed.point_pulsed does: w = casing, wc = cable, up positive; rel = w - wc; the floor (support) gap
    # is GAP - rel and the ceiling (lift) gap GAP + rel; the protection, once a face has closed by 10 mm (sensed with a 1 ms delay),
    # runs a critically damped 5 Hz position loop plus the weight the support is not carrying, clipped to the caps; the support has
    # its own 1 Hz loop (the flux pumps' trim) and a cap; currents follow the commands with the 2 ms lag, the force as current squared.
    lag = 1.0 - math.exp(-dt / TAU_HW)
    cap_u = cap_up if cap_up is not None else 2 * F_LOAD; cap_d = cap_down if cap_down is not None else 2 * F_LIFT
    steady_cap = 2 * r_sup["F_req"]
    om_p = 2 * math.pi * 5.0; om_l = 2 * math.pi * loop_hz                # loop_hz = 0: the support current frozen (bearing_layout, 2026-11-12)
    delay_n = int(round(1e-3 / dt)); k_win = max(int(round(pz.PROT_V_WINDOW / dt)), delay_n, 1)
    buf = np.zeros(delay_n + k_win + 2)
    F_ref_hu = max(cmd.max(), cap_u, 1.0) if kind == "load" else max(cap_u, 1.0)
    F_ref_hd = max(-cmd.min(), cap_d, 1.0) if kind == "lift" else max(cap_d, 1.0)
    I_hu = I_hd = 0.0; I_pu = I_pd = 0.0
    w = v = wc = vc = 0.0
    exc = 0.0; gap_b_min = GAP; gap_t_min = GAP; contact = False; t_contact = float("nan"); engaged = False
    S_min = W; S_hist = np.zeros(nt)
    for i in range(nt):
        rel = w - wc
        buf = np.roll(buf, 1); buf[0] = rel
        rel_meas = buf[delay_n]; v_meas = (rel_meas - buf[delay_n + k_win]) / (k_win * dt) if i >= k_win else 0.0
        # the support pull: both tracks at the casing weight, through the coupled ampere-turns, then the 1 Hz loop and the cap
        NI_s = NI_s0 + dNI[i]
        S_c = 2 * bo.F_of_NI(r_sup, NI_s) if coupled else W
        if gap_dep is not None:
            S_c *= gap_dep(GAP - rel)
        Sj = min(max(S_c - m_c * (2 * om_l * (v - vc) + om_l**2 * rel), 0.0), steady_cap)
        S_min = min(S_min, Sj); S_hist[i] = Sj
        # the planned command through the lag
        c_plan = cmd[i]
        I_hu += (math.sqrt(max(c_plan, 0.0) / F_ref_hu) - I_hu) * lag; I_hd += (math.sqrt(max(-c_plan, 0.0) / F_ref_hd) - I_hd) * lag
        F_hu = min(F_ref_hu * I_hu**2, cap_u); F_hd = min(F_ref_hd * I_hd**2, cap_d)
        if coupled and kind == "load":
            # the load winding's delivered pull: the total ampere-turns through the pull curve, less the depleted support's own
            NIl = max(bo.NI_of_F(r_sup, W_HALF + F_hu / 2) - NI_s0, 0.0) if F_hu > 0 else 0.0
            F_hu = 2 * bo.F_of_NI(r_sup, NI_s + NIl) - S_c
        # the protection
        if prot and not engaged and abs(rel_meas) > 0.010:
            engaged = True
        F_des = ((W - Sj) - m_c * (2 * om_p * v_meas + om_p**2 * rel_meas)) if engaged else 0.0
        F_des = min(max(F_des, -cap_d), cap_u)
        I_pu += (math.sqrt(max(F_des, 0.0) / F_ref_hu) - I_pu) * lag; I_pd += (math.sqrt(max(-F_des, 0.0) / F_ref_hd) - I_pd) * lag
        F_pu = min(F_ref_hu * I_pu**2, cap_u) * (prot_factor if coupled else 1.0); F_pd = min(F_ref_hd * I_pd**2, cap_d)
        F_up = F_hu + F_pu; F_dn = F_hd + F_pd
        a_c = (-W + Sj - q[i] - F_dn + F_up) / m_c
        a_k = (W - Sj + F_dn - F_up) / m_k
        v += a_c * dt; w += v * dt; vc += a_k * dt; wc += vc * dt
        rel = w - wc
        if rel >= GAP:                                              # the casing floor on the cable's underside plate
            w = wc + GAP; v = min(v, vc)
            if not contact: contact = True; t_contact = t[i]
        if rel <= -GAP:                                             # the casing ceiling on the cable's top plate
            w = wc - GAP; v = max(v, vc)
            if not contact: contact = True; t_contact = t[i]
        rel = w - wc
        exc = max(exc, abs(rel)); gap_b_min = min(gap_b_min, GAP - rel); gap_t_min = min(gap_t_min, GAP + rel)
    return dict(u=u, exc=exc, gap_b=gap_b_min, gap_t=gap_t_min, contact=contact, t_contact=t_contact, dNI_min=float(dNI.min()), dNI_max=float(dNI.max()),
                S_min=S_min, engaged=engaged, t=t, cmd=cmd, dNI=dNI, S=S_hist)


def part_1d(h, theva, theva77, r_sup, r_lift, a1, quick=False):
    hr("PART 1D: THE DRIVEN LOAD WINDING AND THE PERSISTENT SUPPORT WINDING SHARE THEIR LOOPS")
    out = {}
    NI_s0 = bo.NI_of_F(r_sup, W_HALF); NI_cap = r_sup["NI_req"]; NI_l = r_sup["NI_load"]; NI_446 = r_sup["NI_446"]
    n_turns = NI_cap / r_sup["I12"] / N_SUP_COILS                  # turns per support coil (2,753 per leg over 10 coils)
    tr = Track(bo.BEST_SUPPORT, h)
    # D.1 inductances
    print(f"  D.1  Inductances per meter of the two loops (own senses: at the design both loops carry +NI and their fluxes are equal), from the")
    print(f"  solver: linear (unsaturated plate) by exciting one loop alone, and incremental about the operating points.")
    L_lin = inductance_matrix(tr, linear=True)
    print(f"    linear: L11 {L_lin[0,0]*1e6:.4f}, L12 {L_lin[0,1]*1e6:.4f}, L22 {L_lin[1,1]*1e6:.4f} uH/m per ampere-turn squared; L11 + L12 = {(L_lin[0,0]+L_lin[0,1])*1e6:.4f} (the design's flux per ampere-turn,")
    print(f"    Part 1a's linear {a1['support']['lin']*1e6:.4f}); coupling of the two loops L12 / L11 = {L_lin[0,1]/L_lin[0,0]:.3f}.")
    L_ops = {}
    for NI in ([NI_s0] if quick else [NI_s0, NI_cap, NI_446]):
        Li = inductance_matrix(tr, linear=False, NI_op=NI)
        L_ops[NI] = Li
        print(f"    incremental at {NI/1e3:.0f} kA: L11 {Li[0,0]*1e6:.4f}, L12 {Li[0,1]*1e6:.4f}, L21 {Li[1,0]*1e6:.4f}, L22 {Li[1,1]*1e6:.4f} uH/m")
    print(f"  The two windings of a track lie in the same bundles, so the mutual inductance between a support turn and a load turn of the same loop is the")
    print(f"  loop inductance itself less the leakage between their positions in the bundle.  Leakage check: a turn in the top half of the bundle against")
    # leakage between the two halves of a bundle: excite the top-half lattice only and read A at the bottom half
    xs0, ys0 = tr.lattices[0]; xs1, ys1 = tr.lattices[1]
    g, s = tr.solve(tr.NI_loop(0, 1e3), linear=True)
    def meanA(g, xs, ys):
        i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
        return float(g.A[i, j].mean())
    top0 = ys0 > ys0.mean(); top1 = ys1 > ys1.mean()
    phi_top = meanA(g, xs0[top0], ys0[top0]) - meanA(g, xs1[top1], ys1[top1]); phi_bot = meanA(g, xs0[~top0], ys0[~top0]) - meanA(g, xs1[~top1], ys1[~top1])
    print(f"  one in the bottom half links {min(phi_top, phi_bot)/max(phi_top, phi_bot):.3f} of the other's flux (the whole loop excited); the coupling coefficient between the windings is above 0.98")
    print(f"  wherever they sit in the bundle, and the model takes it as one.")
    out["L_lin"] = L_lin.tolist(); out["L_ops"] = {f"{k/1e3:.0f}": v.tolist() for k, v in L_ops.items()}
    # D.2 the coil system: the stretch result
    print()
    print(f"  D.2  Flux conservation along the ring.  {N_SUP_COILS} staggered support coils per bundle of {n_turns:,.0f} turns each; a load winding change dNI over a stretch l.")
    print(f"  Each coil keeps its flux: the integral of (dNI_support + dNI_load) over its length is zero.  With the load winding symmetric on both loops the")
    print(f"  inductances drop out and the result is geometric.  The coils of a bundle overlap, staggered by L_c / {N_SUP_COILS}; every point is covered by all")
    print(f"  {N_SUP_COILS}, and a coil that covers part of the stretch must keep its flux too.  The hand estimate (the loss spread as l / L_c over the whole coil)")
    print(f"  is the case of one coil per bundle, no stagger; it is printed beside the staggered design.")
    systems = {}
    print(f"  {'coil length':>11} {'coils/bundle':>12} {'stagger m':>9} {'stretch m':>9} {'load dNI kA':>11} {'support change under':>20} {'+-stagger':>10} {'+-L_c/2':>8} {'+-L_c':>8} {'sum kA m':>9} {'= -l dNI':>9} {'span m':>7}")
    for L_c in L_SUP_COILS:
        for n_fam in (N_SUP_COILS, 1):
            cs = CoilSystem(L_c, L_lin, n_fam=n_fam, cells_per_coil=100 if quick else 200)
            if n_fam == N_SUP_COILS:
                systems[f"{L_c/1e3:.0f} km, {N_SUP_COILS} staggered"] = cs
            elif L_c == max(L_SUP_COILS):
                systems[f"{L_c/1e3:.0f} km, ending together"] = cs
            xc = cs.x[cs.n_cells // 2]
            for l, dNI in ((125.0, NI_cap - NI_s0), (125.0, NI_l), (L_c / 10, 10e3), (2 * L_c, NI_l)):
                d = cs.stretch(xc, l, dNI)
                dNI_s, delta = cs.solve(d)
                prof = dNI_s[0]
                at = lambda off: float(np.interp(xc + off, cs.x, prof))
                tot = float(prof.sum() * cs.dx)
                span = float(np.sum(prof < -0.02 * dNI) * cs.dx)
                print(f"  {L_c/1e3:9.0f} km {n_fam:12d} {L_c/n_fam:9.0f} {l:9.0f} {dNI/1e3:11.0f} {at(0)/1e3:17.2f} kA {at(L_c/n_fam)/1e3:10.2f} {at(L_c/2)/1e3:8.2f} {at(L_c)/1e3:8.2f} {tot/1e3:9.0f} {-l*dNI/1e3:9.0f} {span:7.0f}")
    print(f"  With one coil per bundle the loss is l / L_c of the load's ampere-turns along the whole coil, as the hand estimate says.  With {N_SUP_COILS} staggered")
    print(f"  coils the loss is confined to the stagger cells the stretch touches and is l / stagger there: the stagger, not the coil length, is the")
    print(f"  dilution length.  A stretch longer than the stagger is cancelled almost completely where it lies, with the remainder pushed into the")
    print(f"  neighbouring cells.  The sum over the ring is minus l x dNI in every case: the driven winding cannot change the total ampere-turns of")
    print(f"  a stretch of persistent coils, it can only move them along the ring, and staggering keeps it from moving them far.")
    # D.3 the kernel for the passage model: response to a unit stretch (NI x m) at offset xi
    kernels = {}
    for key, cs in systems.items():
        L_c = cs.L_c
        xc = cs.x[cs.n_cells // 2]
        d = cs.stretch(xc, cs.dx, 1.0)
        dNI_s, delta = cs.solve(d)
        R = dNI_s[0] / cs.dx                                       # per (ampere-turn x m)
        xi = cs.x - xc
        m = np.abs(xi) <= 3 * L_c
        kernels[key] = (xi[m], R[m])
    # D.4 Table 7.6 rows with flux-conserving coils
    print()
    print(f"  D.4  Table 7.6 with flux-conserving support coils.  The load winding's stretch: a parked 100 m spacecraft energizes the 25 m coils that cover it,")
    print(f"  {L_SC + L_LOAD:.0f} m; a launch at 1 km/s the moving footprint of {1000*max(0.05, L_SC/1000) + L_LOAD:.0f} m (the row is the static picture; D.6 moves it); a lost track is taken as a")
    print(f"  lost stretch of 100 m, one support coil length and ten; the protection pulse over 50 m (a local drift) and over 2 km (a long one).  Support current")
    print(f"  = the persistent ampere-turns at the stretch after the coils have kept their flux; pull from the track's pull curve at the local total; the")
    print(f"  'rest of the coils' column is the pull where the depletion is largest outside the stretch (its edge) and the span it covers.")
    NI_lost = bo.NI_of_F(r_sup, 2 * W_HALF) - NI_s0
    rows76 = [("Parked 5,000 t before the pumps (141 + 182 over the load)", NI_s0, NI_cap - NI_s0, L_SC + L_LOAD, r_sup["F_req"], W_HALF),
              ("Launch at 1 km/s, static (141 + 141)", NI_s0, NI_s0, 1000 * max(0.05, L_SC / 1000) + L_LOAD, W_HALF + 285e3 / 2, W_HALF),
              ("One track lost at rest over 100 m: the other's load winding catches the other half (141 + 58)", NI_s0, NI_lost, 100.0, 2 * W_HALF, W_HALF),
              ("One track lost over one coil length, the same", NI_s0, NI_lost, None, 2 * W_HALF, W_HALF),
              ("One track lost over ten coil lengths, the same", NI_s0, NI_lost, "10", 2 * W_HALF, W_HALF),
              ("One track out and a standing load, at the limit (322 + 171 over the load)", NI_cap, NI_446 - NI_cap, L_SC + L_LOAD, r_sup["F_446"], r_sup["F_req"]),
              ("Protection pulse at the load cap over 50 m (141 + 287)", NI_s0, NI_l, 50.0, W_HALF, W_HALF),
              ("Protection pulse at the load cap over 2 km", NI_s0, NI_l, 2000.0, W_HALF, W_HALF)]
    print(f"  {'case':<92} {'coils':<24} {'support kA':>10} {'load kA':>7} {'pull under kN/m':>15} {'book':>6} {'needed':>7} {'rest: min support kA':>20} {'pull kN/m':>9} {'needed':>7} {'span m':>7}")
    table76 = []
    for name, NIs, NIl, l, F_need_under, need_out in rows76:
        for key, cs in systems.items():
            L_c = cs.L_c
            ll = L_c if l is None else (10 * L_c if l == "10" else l)
            if ll > 0.5 * cs.x[-1]:
                ll = min(ll, 0.5 * cs.x[-1])
            xc = cs.x[cs.n_cells // 2]
            dNI_s, delta = cs.solve(cs.stretch(xc, ll, NIl))
            prof = dNI_s[0]
            inside = np.abs(cs.x - xc) < ll / 2
            s_under = NIs + float(prof[inside].mean()); F_under = bo.F_of_NI(r_sup, s_under + NIl)
            outside = ~inside
            k = np.argmin(prof[outside]); s_min = NIs + float(prof[outside][k]); F_min = bo.F_of_NI(r_sup, s_min)
            span = float(np.sum(prof < -0.05 * NIl) * cs.dx) - ll if NIl else 0.0
            F_book = bo.F_of_NI(r_sup, NIs + NIl)
            print(f"  {name:<92} {key:<24} {s_under/1e3:10.1f} {NIl/1e3:7.0f} {F_under/1e3:15.1f} {F_book/1e3:6.0f} {F_need_under/1e3:7.1f} {s_min/1e3:20.1f} {F_min/1e3:9.1f} {need_out/1e3:7.1f} {max(span, 0):7.0f}")
            table76.append(dict(case=name, coils=key, L_c=L_c, support_under=s_under, load=NIl, F_under=F_under, F_book=F_book, F_need_under=F_need_under, support_min=s_min, F_min=F_min, F_need=need_out, span=span))
    out["table76"] = table76
    print(f"  Reading: with the book's 1 km coils in ten staggered families the stagger is 100 m, shorter than every load stretch, so the persistent")
    print(f"  coils cancel most of what the load winding adds where it adds it and take the rest out of the cells beside it; the pull under a parked")
    print(f"  spacecraft is far below the book's row and the track beside it loses most of its weight share.  With 10 km coils (1 km stagger) the loss")
    print(f"  is the hand estimate's share of the load over the neighbouring kilometre.  A lost track over a stretch longer than the stagger cannot be")
    print(f"  caught by the other track's load winding over that stretch at all; the load winding there cannot restore it (D.2), and the flux pumps must.")
    # D.5 the gap: stiffness at constant current and at constant flux
    print()
    print(f"  D.5  The gap that follows.  Gap sweep of the support track (linear solves, the plate unsaturated as it is at 141 kA): the loop inductance and")
    print(f"  the pull per ampere-turn squared against the gap; the pull at constant current and at constant flux (NI = Phi' / L'(g)); their slopes at 100 mm.")
    gaps = [0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13] if not quick else [0.08, 0.10, 0.12]
    sw = gap_sweep(bo.BEST_SUPPORT, h, gaps, NI_nl=None if quick else NI_cap)
    k_I, k_phi, F_I, F_phi = stiffness(sw, GAP, NI_s0)
    print(f"  {'gap mm':>6} {'L uWb/m per A':>13} {'F/NI^2 nN/A^2':>13} {'F at 141 kA, const I kN/m':>25} {'const flux kN/m':>15}" + ("" if quick else f" {'F at 322 kA kN/m':>16}"))
    for r, fi, fp in zip(sw, F_I, F_phi):
        print(f"  {r['gap']*1e3:6.0f} {r['Lp']*1e6:13.4f} {r['Fc']*1e9:13.3f} {fi/1e3:25.1f} {fp/1e3:15.1f}" + ("" if quick else f" {r['F_nl']/1e3:16.1f}"))
    print(f"  Slope at 100 mm: constant current {k_I/1e3/1e3:+.2f} kN/m per mm ({k_I/W_HALF*1e-3*100:+.1f} percent of the weight share per mm), constant flux {k_phi/1e3/1e3:+.2f} kN/m per mm")
    print(f"  ({k_phi/W_HALF*1e-3*100:+.1f} percent per mm).  A negative slope is the EMS instability (the pull falls as the gap opens).  At constant flux the current rises as")
    print(f"  the inductance falls with the gap, which {'turns the track stable: the pull rises as the gap opens' if k_phi > 0 else 'takes ' + f'{(1 - k_phi/k_I)*100:.0f}' + ' percent of the instability away but leaves the track unstable'}.")
    print(f"  Constant flux holds for a change that is uniform along a coil (the loss of pull over the depleted span is one); a local gap change over a")
    print(f"  stretch l sees only l / L_c of it.")
    out["gap_sweep"] = sw; out["k_I"] = k_I; out["k_phi"] = k_phi
    # the sag of a depleted span at constant flux: the gap at which the constant-flux pull equals the needed share
    gs = np.array([r["gap"] for r in sw]); i0 = int(np.argmin(np.abs(gs - GAP)))
    for name, key, s_min in [(t["case"], t["coils"], t["support_min"]) for t in table76 if t["case"].startswith("Parked") or t["case"].startswith("Launch")]:
        Fphi_dep = F_phi * (s_min / NI_s0)**2                        # the depleted span at constant (its new, lower) flux
        d = 1 - Fphi_dep[i0] / W_HALF
        t_fall = math.sqrt(2 * GAP / (d * W_HALF * 2 / ring.M_LOAD_M)) if d > 0 else float("inf")
        above = np.nonzero(Fphi_dep[i0:] >= W_HALF)[0]
        if d > 0 and len(above):
            j = i0 + above[0]
            g_eq = float(gs[j - 1] + (W_HALF - Fphi_dep[j - 1]) * (gs[j] - gs[j - 1]) / (Fphi_dep[j] - Fphi_dep[j - 1])) if j > i0 else GAP
            print(f"    {name}, {key} coils: the depleted span is {d:.0%} short of its weight share at 100 mm; at constant flux its pull is back to the share at a gap")
            print(f"    of {g_eq*1e3:.0f} mm, so the span settles {(g_eq-GAP)*1e3:.0f} mm lower (the ceiling gap {GAP*1e3 - (g_eq-GAP)*1e3:.0f} mm) until the pumps restore the current; the fall there takes about {t_fall:.2f} s.")
        elif d > 0:
            print(f"    {name}, {key} coils: the depleted span is {d:.0%} short of its weight share and the constant-flux pull does not recover it within the")
            print(f"    sweep; free fall under that deficit closes the 100 mm ceiling gap in {t_fall:.2f} s if nothing else acts (both tracks depleted alike).")
    # D.6 the sled passage with the coupled model
    print()
    print(f"  D.6  Sled passage below orbital velocity with the coupled model (point model as described in the header; the planned command is the")
    print(f"  uncoupled feed-forward, the support pull follows the coupled coil currents).  The protection is run two ways: as a local pulse over two")
    print(f"  load coil lengths, credited with 1 - 50 m / stagger of its pull (the rest is taken back by the cell), and off, which is the case of a")
    print(f"  sag over a whole stagger cell, where the load winding can add nothing.")
    print(f"  {'u km/s':>6} {'coils':<24} {'coupled':>7} {'protection':>12} {'casing movement mm':>18} {'floor gap mm':>12} {'ceiling gap mm':>14} {'support dip kA':>14} {'contact':>8} {'at s':>6}")
    pass_rows = []
    for u_kms, F_pair, kind in PASS_PULLS:
        if kind != "load":
            continue
        u = u_kms * 1e3
        base = passage_point(u, L_LOAD, F_pair * 1e3, kind, r_sup=r_sup, coupled=False)
        print(f"  {u_kms:6.0f} {'-':<24} {'no':>7} {'on':>12} {base['exc']*1e3:18.1f} {base['gap_b']*1e3:12.1f} {base['gap_t']*1e3:14.1f} {0.0:14.1f} {('yes' if base['contact'] else 'no'):>8} {base['t_contact']:6.2f}")
        pass_rows.append(dict(u=u, L_c=None, coupled=False, prot="on", **{k: base[k] for k in ("exc", "gap_b", "gap_t", "contact", "t_contact")}))
        for key, cs in systems.items():
            L_c = cs.L_c; stagger = L_c / cs.n_fam
            for plabel, pon, pf in ((f"local x {max(0.0, 1 - 2 * L_LOAD / stagger):.2f}", True, max(0.0, 1 - 2 * L_LOAD / stagger)), ("off", False, 0.0)):
                r = passage_point(u, L_LOAD, F_pair * 1e3, kind, kernel=kernels[key], L_c=L_c, coupled=True, r_sup=r_sup, prot=pon, prot_factor=pf, t_settle=1.0 + L_c / u)
                print(f"  {u_kms:6.0f} {key:<24} {'yes':>7} {plabel:>12} {r['exc']*1e3:18.1f} {r['gap_b']*1e3:12.1f} {r['gap_t']*1e3:14.1f} {r['dNI_min']/1e3:14.1f} {('yes' if r['contact'] else 'no'):>8} {r['t_contact']:6.2f}")
                pass_rows.append(dict(u=u, coils=key, L_c=L_c, coupled=True, prot=plabel, **{k: r[k] for k in ("exc", "gap_b", "gap_t", "contact", "t_contact", "dNI_min")}))
    out["passes"] = pass_rows
    print(f"  (the uncoupled rows are this model's own baseline against final_bearing's 7.6, 2.6 and 1.3 mm; the coupled rows carry the support's loss of")
    print(f"  current as the footprint passes the point's stagger cell: with 1 km coils the cell is 100 m and the loss is most of the footprint's")
    print(f"  ampere-turns, with 10 km coils the cell is 1 km and the loss is the footprint's share of it; the point then sees a pull deficit for as")
    print(f"  long as the footprint is in or beside its cell, and the response kernel is that of a footprint centred in a cell)")
    # D.7 the support winding: induced currents, voltages and heat
    print()
    print(f"  D.7  What the pulses do inside the support winding.  A load pulse of dNI over a stretch l in the 2 ms rise changes the covering coils'")
    print(f"  currents (D.2) by up to the amounts below per tape ({n_turns:,.0f} turns per coil); the net EMF around a closed superconducting turn is zero at")
    print(f"  every instant (the stretch's dPhi/dt is balanced by the coil's own L dI/dt along its whole length), so no voltage appears between the")
    print(f"  soldered turns beyond the joint drop of the current change; the heat is the joints' I^2 R of that change [joint resistance 10 nohm per")
    print(f"  joint, one joint per turn, UNVERIFIED], and the tape's own hysteresis in the pulse's field change, which is the same strip estimate as")
    print(f"  for the driven tapes with the 12 mm width and limit:")
    dNI_p = NI_l; l_p = 50.0
    for key, cs in systems.items():
        xc = cs.x[cs.n_cells // 2]
        dNI_s, delta = cs.solve(cs.stretch(xc, l_p, dNI_p))
        n_t = NI_cap / r_sup["I12"] / cs.n_fam
        dI = float(np.abs(delta).max()) / n_t
        E_j = 4 * cs.n_fam * n_t * dI**2 * 10e-9 * TAU_HW       # four bundles, the change flowing through every joint for the rise
        print(f"    {key} coils, protection pulse of {dNI_p/1e3:.0f} kA over {l_p:.0f} m: up to {dI:.2f} A per tape in 2 ms ({dI/TAU_HW:,.0f} A/s); joint heat {E_j*1e6:.3f} uJ per pulse per track.")
    # the hysteresis of the 12 mm tapes under the load pulse: the field change at the support lattice between 141 and 428 kA
    g1, s1 = tr.solve(tr.NI_design(NI_s0)); B1 = [bo.field_at_points(g1, xs, ys) for xs, ys in tr.lattices]
    g2, s2 = tr.solve(tr.NI_design(NI_s0 + NI_l)); B2 = [bo.field_at_points(g2, xs, ys) for xs, ys in tr.lattices]
    dB = np.concatenate([np.hypot(b2[0] - b1[0], b2[1] - b1[1]) for b1, b2 in zip(B1, B2)])
    Ic12 = r_sup["ic12"]; Ic3 = r_sup["ic3_446"]
    q12 = 0.5 * float(dB.mean()) * 0.012 * Ic12 * (r_sup["turns12"])          # J per meter of track per pulse, all the field through the face
    q3 = 0.5 * float(dB.mean()) * 0.003 * Ic3 * (r_sup["turns3_446"])
    print(f"    field change at the support tapes for the protection pulse (141 to 428 kA): mean {dB.mean():.2f} T, largest {dB.max():.2f} T.  Strip estimate (all of it through the")
    print(f"    face, the bound Section 7.7 uses): 12 mm support winding {q12:,.0f} J per meter of track per pulse against the 3 mm load winding's {q3:,.0f} J in the same")
    print(f"    bundles, {q12/q3:.1f} times, because the loss per tape goes with width x limit and the 12 mm tape has four times both at a quarter of the count.")
    print(f"    The book's 605 GJ bound counts only the driven tapes; the persistent tapes under them add {q12/q3:.1f} times that for the pulses below orbital velocity.")
    out["dB_support"] = dict(mean=float(dB.mean()), max=float(dB.max()), q12=q12, q3=q3)
    # D.8 decoupling: a single wide loop for the load winding on the figure-eight support
    print()
    print(f"  D.8  Removing the coupling.  (a) The support loops in series as one figure-eight circuit and the load winding wound the same way on both")
    print(f"  loops (not opposite): its center legs cancel and it is one 1.15 m loop on the outer bundles; its flux through loop 1 and loop 2 is equal and")
    print(f"  opposite in the figure-eight's sense, so the persistent circuit links none of it.  The solver with the support at 141 kA (figure-eight) and")
    print(f"  the load on the outer legs only:")
    def outer(NIl):
        v = tr.NI_design(0.0); v = [0.0] * len(tr.legs)
        v[0] += NIl; v[-1] -= NIl
        return v
    for NIl in ([NI_s0] if quick else [NI_s0, NI_l]):
        legs_c = [a + b for a, b in zip(tr.NI_design(NI_s0), outer(NIl))]
        g, s = tr.solve(legs_c)
        g_l, s_l = tr.solve(outer(NIl))
        g_s, s_s = tr.solve(tr.NI_design(NI_s0))
        wt, _ = worst_tape(g, tr, legs_c, 0.012, theva, theva77)
        # the mutual: flux of the figure-eight (phi1 + phi2 in own senses) from the single loop alone
        print(f"    load {NIl/1e3:.0f} kA on the outer legs: pull alone {s_l['F']/1e3:.1f} kN/m (the figure-eight load winding gives {bo.F_of_NI(r_sup, NIl)/1e3:.0f} at the same ampere-turns), with the")
        print(f"    support's 141 kA {s['F']/1e3:.1f} kN/m (support alone {s_s['F']/1e3:.1f}; the sum {(s_l['F']+s_s['F'])/1e3:.1f}: no cross term); figure-eight flux from the single loop {(s_l['phi'][0]+s_l['phi'][1]):+.4f} Wb/m")
        print(f"    (zero by symmetry); worst 12 mm tape {wt[1]:.2f} T at {wt[2]:.0f} degrees, {wt[0]:.0f} A (the outer bundle carries both windings' current on one side).")
    print(f"    Cost: the single loop needs about {math.sqrt(bo.F_of_NI(r_sup, NI_l)/max(s_l['F'],1.0))*NI_l/1e3:.0f} kA for the 198 kN/m the figure-eight gives at {NI_l/1e3:.0f} kA, "
          f"{(math.sqrt(bo.F_of_NI(r_sup, NI_l)/max(s_l['F'],1.0))-1)*100:.0f} percent more tape, and its field adds to the")
    print(f"    support's in one outer bundle and subtracts in the other, so the pull sits off-centre on the track (a roll moment on the casing).")
    print(f"  (b) Separate loops for the load winding, beside or below the support bundles, keep a mutual of the order of the geometric overlap; only the")
    print(f"  symmetry of (a) makes it zero.  (c) The stagger.  The dilution length is the stagger (D.2): coils of a bundle that end together (one coil")
    print(f"  per bundle, or ten that start and end at the same places) spread the loss as l / L_c over the coil, and 10 km coils of that kind take about")
    print(f"  {125/10000*100:.1f} percent of a parked spacecraft's load winding out of the support, which the constant-flux response (D.5) carries with a few")
    print(f"  millimetres of sag.  Ten coils ending together also keep each other's ampere-turns when one quenches (each keeps its flux, so the nine pick up")
    print(f"  the tenth's share, cell by cell in the staggered case too).  The price is the coil-end ripple of Section 7.8: one full field step per coil")
    print(f"  length instead of ten steps of a tenth, ten times the eddy heat of the staggered coil of the same length, so 10 km coils ending together")
    print(f"  heat the plates like the book's 1 km staggered ones, {4.1:.1f} W/m [final_bearing Part C.3's model]; the step itself and the bend of a whole")
    print(f"  bundle's 2,753 turns to their pumps at one place are a coil-end design the chapter does not have.  (d) A driven support winding (no")
    print(f"  persistent loops) removes the constraint and gives up the persistent mode's zero-lead, zero-supply holding of the casing.")
    # D.9 lift-support coupling through the cable
    print()
    print(f"  D.9  Coupling between the lift tracks on the ceiling and the support tracks through the cable's plates (one 20 x 20 m linear box, 50 mm cells,")
    print(f"  both plates, the lift loop excited alone, the flux through the support loops read from A).")
    side = ring.CABLE_SIDE
    gbig = ms.Grid(20.0, 20.0, 0.05, x0=-10.0, y0=-8.0)
    xc = fb.Y_TRACK
    gbig.add_iron(xc - 0.7, xc + 0.7, 0.0, PLATE_T, cf.PLATE_MU); gbig.add_iron(xc - 0.7, xc + 0.7, side - PLATE_T, side, cf.PLATE_MU)
    NIl_big = 1e3
    for dx, a in ((-0.4, NIl_big), (0.4, -NIl_big)):
        gbig.add_current(xc + dx - 0.1, xc + dx + 0.1, side + PLATE_T + 0.01, side + PLATE_T + 0.06, a)
    gbig.solve()
    def Abox(x_lo, x_hi, y_lo, y_hi):
        m = gbig.box(x_lo, x_hi, y_lo, y_hi); return float(gbig.A[m].mean())
    y_s1, y_s0 = -0.01, -0.0767
    A_s = [Abox(xc + x - 0.075, xc + x + 0.075, y_s0, y_s1) for x in (-0.575, -0.075, 0.075, 0.575)]
    phi_s1 = A_s[0] - A_s[1]; phi_s2 = A_s[3] - A_s[2]
    A_l = [Abox(xc + x - 0.1, xc + x + 0.1, side + PLATE_T + 0.01, side + PLATE_T + 0.06) for x in (-0.4, 0.4)]
    phi_ll = A_l[0] - A_l[1]
    print(f"    lift loop self flux {phi_ll/NIl_big*1e6:.3f} uWb/m per A (coarse); flux through support loop 1 {phi_s1/NIl_big*1e9:+.3f}, loop 2 {phi_s2/NIl_big*1e9:+.3f} nWb/m per A: {abs(phi_s1)/phi_ll*100:.3f} percent of")
    print(f"    the lift loop's own, {abs(phi_s1)/(L_lin[0,0]+L_lin[0,1])/NIl_big*100:.3f} percent of a support loop's; a full lift pulse of {r_lift['NI_req']/1e3:.0f} kA moves the support's flux as {abs(phi_s1)/NIl_big*r_lift['NI_req']/((L_lin[0,0]+L_lin[0,1])*NI_s0)*100:.2f} percent of its")
    print(f"    own at 141 kA would: small, as expected; the {side:.2f} m of cable and the two plates between them short the field.")
    out["lift_support_coupling"] = dict(phi_ll=phi_ll / NIl_big, phi_s1=phi_s1 / NIl_big, phi_s2=phi_s2 / NIl_big)
    return out


# =============================================================================

def main():
    part = "all"; h = 0.005; quick = False
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--part="):
            part = arg.split("=")[1]
        elif arg.startswith("--h="):
            h = float(arg.split("=")[1]) * 1e-3
        elif arg == "--quick":
            quick = True
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    hr("THE ELECTRICAL DESIGN OF THE DRIVEN COILS (2026-11-12)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km: casing {ring.M_LOAD_M:,.0f} kg/m ({ring.M_LOAD_M*ring.G_NET/1e3:.1f} kN/m), cable {ring.M_CABLE_M:,.0f} kg/m.  Grid {h*1e3:g} mm{' (quick)' if quick else ''}.")
    if part == "sweep":
        sweep_cases(); hr(f"DONE ({time.time()-t_all:.0f} s)"); return
    theva, theva77 = load_theva(); r_sup, r_lift = load_records()
    cache = os.path.join(ring.out_dir(), f"bearing_electrical_{h*1e3:g}mm.json")
    res = {}
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            res = json.load(f)
    def save():
        with open(cache, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (o.tolist() if isinstance(o, np.ndarray) else str(o)))
    if part in ("1a", "all") or "1a" not in res:
        res["1a"] = part_1a(h, theva, theva77, r_sup, r_lift, quick); save()
    if part in ("1b", "all"):
        res["1b"] = part_1b(h, theva, theva77, r_sup, r_lift, res["1a"], quick); save()
    if part in ("1c", "all"):
        res["1c"] = part_1c(h, theva, theva77, r_sup, r_lift, res["1a"], res.get("1b"), quick); save()
    if part in ("1d", "all"):
        res["1d"] = part_1d(h, theva, theva77, r_sup, r_lift, res["1a"], quick); save()
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
