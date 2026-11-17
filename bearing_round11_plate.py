#!/usr/bin/env python3
"""
Round 11, item 1b (2026-11-17): one plate under the cable, and a wider face.  Vol_III/Claude_Code_instructions_round11_2026-11-17.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing's plates and loop widths).

The design is the one of rounds 7 and 8 at the 100 mm gap: two support tracks (two separate 0.4 m loops 50 mm apart at the center, 150 x 67 mm
bundles of 12 mm tape, a 1.2 m x 100 mm plate), two load tracks under the cable and two lift tracks over it (one 0.6 m loop, 204 x 50 mm bundles
of 3 mm tape, a 1.0 m x 200 mm plate).  This script puts the tracks of one face in one field model, with a plate whose thickness is given column
by column, and reads the pull of each track on its own coils.

The model (magnetostatic_2d.Grid, bearing_optimize.Window's nesting, bearing_electrical.Track's legs and tape lattices):
  x across the face, y up; the bundles' tops at y = -10 mm, the plate's face at y = 100 mm, the plate's thickness T(x) above it (the face toward
  the coils is flat, so every track is one gap from the plate); sintered iron, mu_r 200, saturating at 1.2 T (compensation_fields);
  a coarse box at 10 mm cells, then a window at h (5 mm by default) with its edge values from the coarse solve;
  the pull of a track is the Maxwell stress around a box that holds its bundles and nothing else (so two tracks under one plate are told apart);
  a half face has the cable's centerline at x = 0 with A = 0 there (the tracks of the other half then carry the mirrored, opposite currents).

Run (from Python/orbitalring):
  python bearing_round11_plate.py --part=base              the design's tracks on their own plates in this model, against the gap trade's 5 mm rows
  python bearing_round11_plate.py --part=tailor            a plate sized from its flux, at the design's loop widths, for a list of fill fields
  python bearing_round11_plate.py --part=width             wider loops on a plate sized from its flux
  python bearing_round11_plate.py --part=alone             a compensation track that gives the pair's whole pull alone
  python bearing_round11_plate.py --part=face              the half face of 7.39 m: four plates, one plate with a web, the plate sized from its flux
  python bearing_round11_plate.py --part=coupling          the load coil's flux through the support loops (round 7's method) for each plate and spacing
  python bearing_round11_plate.py --part=mesh              the cell size and box checks
  python bearing_round11_plate.py --part=report            the tables (bearing_round11_plate_report.py)
Options: --h=0.005, --jobs=2.  Writes runs/round11_logs/round11_plate_*.json and .txt.  Environment: ORBITAL_RING_OUT (default runs).
"""
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import bearing_round6b as r6
import bearing_round7 as r7
import bearing_gap_trade as gt
import compensation_fields as cf
import magnetostatic_2d as ms

OUT = os.environ.get("ORBITAL_RING_OUT", "runs")
LOGS = os.path.join(OUT, "round11_logs")
MU0 = 4e-7 * math.pi
GAP = bd.GAP0                                   # 0.100 m
MU = cf.PLATE_MU                                # 200
BS = cf.PLATE_B_SAT                             # 1.2 T
RHO = bd.RHO_PLATE                              # 7,870 kg/m3
Y_TOP = -bo.WALL_DESIGN                         # the bundles' tops, 10 mm behind the casing's face
HP = 0.005                                      # m, the column width in which a plate's thickness is given
GC = 0.05                                       # m between the two center bundles of a two-loop track
SUP = r7.SUP; FAST = r7.FAST
SUP_BUNDLE = gt.SUP_BUNDLE
F_FAST = gt.F_FAST                              # 196.5 kN/m, a lift track's (and a load track's) requirement
F_FULL = r7.F_FULL                              # 198 kN/m, a load track's full pull (the coupling's state, as round 7)
W_HALF = r7.W_HALF                              # 47.5 kN/m, a support track in normal operation
F_ALONE = 396e3                                 # N/m, a compensation track doing the pair's job alone: the parked 5,000 t on 100 m (bearing_design.F_PARKED);
                                                # the lift pair's 30 km/s pull is 393 kN/m (bearing_optimize.PASS_PULLS), 1 percent under it
HALF_FACE = 3.695                               # m, half of the 7.39 m face
TOL = 0.003                                     # the solved pull within 0.3 percent of the requirement before the tape is counted
M_LIM_PLATES = 1872.0                           # kg/m (Table 7.5)
HW_ALLOW = 10686.0                              # kg/m, the hardware the cable is sized for
BEARING_R8 = 1687.0                             # kg/m, the bearing on the casing (round 8)


def save(name, obj):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: bool(o) if isinstance(o, np.bool_) else (o.tolist() if isinstance(o, np.ndarray) else float(o)))


def load(name):
    p = os.path.join(LOGS, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


class Tee:
    def __init__(self, path, mode="a"):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.f = open(path, mode, encoding="utf-8"); self.o = sys.stdout             # a part run again adds to its log

    def write(self, s):
        self.o.write(s); self.f.write(s)

    def flush(self):
        self.o.flush(); self.f.flush()


# =============================================================================
# the plate, the tracks, the model
# =============================================================================

class Plate:
    """A plate on the cable's face: its thickness in columns HP wide from x0 on (m).  The face toward the coils is at y = GAP."""

    def __init__(self, x0, T, hp=HP):
        self.x0 = float(x0); self.hp = hp; self.T = np.asarray(T, dtype=float).copy()

    @classmethod
    def rect(cls, a, b, T, hp=HP):
        return cls(a, np.full(int(round((b - a) / hp)), T), hp)

    @property
    def x1(self):
        return self.x0 + len(self.T) * self.hp

    def centers(self):
        return self.x0 + (np.arange(len(self.T)) + 0.5) * self.hp

    def T_at(self, X):
        k = np.floor((np.asarray(X) - self.x0) / self.hp + 1e-6).astype(int)
        ok = (k >= 0) & (k < len(self.T))
        out = np.zeros(np.shape(X)); out[ok] = self.T[k[ok]]
        return out

    def mask(self, g):
        return (g.Y >= GAP - 1e-9) & (g.Y < GAP + self.T_at(g.X) - 1e-9)

    def area(self, a=None, b=None):
        c = self.centers(); m = np.ones(len(c), bool) if a is None else (c >= a) & (c < b)
        return float(self.T[m].sum() * self.hp)

    def mass(self, a=None, b=None):
        return RHO * self.area(a, b)

    def extent(self, t_min=0.0):
        k = np.nonzero(self.T > t_min + 1e-9)[0]
        return (self.x0 + k[0] * self.hp, self.x0 + (k[-1] + 1) * self.hp) if len(k) else (self.x0, self.x0)

    def shifted(self, dx):
        return Plate(self.x0 + dx, self.T, self.hp)

    def to_json(self):
        return dict(x0=self.x0, hp=self.hp, T=[round(float(t), 6) for t in self.T])

    @classmethod
    def from_json(cls, d):
        return cls(d["x0"], d["T"], d["hp"])

    @classmethod
    def combine(cls, plates, a, b, floor=0.0, hp=HP):
        """The plates laid on one face from a to b: the thicknesses added where they overlap, and a web of at least floor everywhere."""
        n = int(round((b - a) / hp)); c = a + (np.arange(n) + 0.5) * hp
        T = np.zeros(n)
        for p in plates:
            T += p.T_at(c)
        return cls(a, np.maximum(T, floor), hp)


class Trk:
    """A track's bundles as bearing_electrical.Track lays them (legs, loops, the tape lattice of the field model), centered on xc.
    kind 'support': 12 mm tape, two loops by default; 'fast': 3 mm tape, one loop by default.  pol turns every current over."""

    def __init__(self, name, kind, xc, Wl, bw, bh, n=None, gc=GC, pol=1):
        self.name, self.kind, self.xc, self.Wl, self.bw, self.bh, self.pol = name, kind, xc, Wl, bw, bh, pol
        if kind == "support":
            geom = dict(bo.BEST_SUPPORT); geom.update(W=Wl, gc=gc, shape=r6.shape_name(bw, bh), n=n or 2)
            self.w_tape = 0.012; I_ref = bo.I_REF_12; NI_ref = 437e3
        else:
            geom = dict(bo.BEST_LIFT); geom.update(W=Wl, gc=gc, shape=r6.shape_name(bw, bh), n=n or 1)
            self.w_tape = 0.003; I_ref = bo.I_REF_3; NI_ref = 452e3
        t = be.Track(geom, 0.005, w_tape=self.w_tape, I_ref=I_ref, NI_ref=NI_ref)
        self.n = geom["n"]; self.gc = gc
        self.legs = [(x + xc, s * pol) for x, s in t.legs]
        self.loops = [tuple((i, s * pol) for i, s in loop) for loop in t.loops]
        self.lat = [(xs + xc, ys) for xs, ys in t.lattices]
        self.x_lo = self.legs[0][0] - bw / 2; self.x_hi = self.legs[-1][0] + bw / 2
        self.span = self.x_hi - self.x_lo

    def design(self, NI):
        return [s * NI for x, s in self.legs]

    def loop(self, k, NI):
        v = [0.0] * len(self.legs)
        for i, s in self.loops[k]:
            v[i] = s * NI
        return v

    def bundles_for_limits(self, NI_legs):
        return [dict(x=x, NI=a, xs=xs, ys=ys) for (x, s), a, (xs, ys) in zip(self.legs, NI_legs, self.lat)]


class Win(bo.Window):
    """bearing_optimize.Window (edge values from a coarser solve) with the sparse direct solve at every size."""

    def _solve_direct(self):
        nx, ny, h = self.nx, self.ny, self.h
        t0 = time.time()
        nu = 1.0 / self.mu_r
        nxw = np.zeros((nx, ny)); nxe = np.zeros((nx, ny)); nys = np.zeros((nx, ny)); nyn = np.zeros((nx, ny))
        nxw[1:, :] = 2 * nu[1:, :] * nu[:-1, :] / (nu[1:, :] + nu[:-1, :])
        nxe[:-1, :] = nxw[1:, :]
        nys[:, 1:] = 2 * nu[:, 1:] * nu[:, :-1] / (nu[:, 1:] + nu[:, :-1])
        nyn[:, :-1] = nys[:, 1:]
        diag = nxw + nxe + nys + nyn
        rhs = MU0 * self.J * h * h
        idx = -np.ones((nx, ny), dtype=np.int64)
        inner = ~self.edge
        n_in = int(inner.sum())
        idx[inner] = np.arange(n_in)
        ii, jj = np.nonzero(inner)
        rows = [idx[ii, jj]]; cols = [idx[ii, jj]]; vals = [diag[ii, jj]]
        b = rhs[ii, jj].copy()
        for di, dj, coef in ((-1, 0, nxw), (1, 0, nxe), (0, -1, nys), (0, 1, nyn)):
            i2, j2 = ii + di, jj + dj
            ok = inner[i2, j2]
            rows.append(idx[ii[ok], jj[ok]]); cols.append(idx[i2[ok], j2[ok]]); vals.append(-coef[ii[ok], jj[ok]])
            b[~ok] += coef[ii[~ok], jj[~ok]] * self.A_edge[i2[~ok], j2[~ok]]
        M = sps.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_in, n_in))
        A = np.array(self.A_edge, copy=True)
        A[ii, jj] = spla.spsolve(M, b)
        self.A = A
        self.Bx = np.zeros((nx, ny)); self.By = np.zeros((nx, ny))
        self.Bx[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
        self.By[1:-1, :] = -(A[2:, :] - A[:-2, :]) / (2 * h)
        self.iters = 1
        self.solve_time = time.time() - t0
        return A


class Face:
    """Tracks under one plate (or several: a Plate with zero thickness between them).  box: the coarse grid (x_lo, x_hi, y_lo, y_hi), A = 0 on its
    edge; win: the fine window."""

    def __init__(self, tracks, plate, box, win, h=0.005, h0=0.01, outer=None, h_outer=0.02):
        """outer: a larger box (x_lo, x_hi, y_lo, y_hi) at h_outer cells around the coarse box, which then takes its edge values from it."""
        self.tracks, self.plate, self.box, self.win, self.h, self.h0, self.outer, self.h_outer = tracks, plate, box, win, h, h0, outer, h_outer
        self.solves = 0

    def _build(self, g, NI, linear):
        m = self.plate.mask(g)
        g.mu_r[m] = MU
        if not linear:
            g.materials.append((m, MU, BS))
        for t in self.tracks:
            for (xs, ys), a in zip(t.lat, NI.get(t.name, [])):
                if a == 0.0:
                    continue
                i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
                np.add.at(g.J, (i, j), (a / len(xs)) / (g.h * g.h))
        return m

    def plate_box(self, g):
        a, b = self.plate.extent()
        return (max(a - 0.1, g.x[0]), min(b + 0.1, g.x[-1]), GAP - 2 * g.h, GAP + float(self.plate.T.max()) + 2 * g.h)

    def coil_box(self, t):
        """A box around a track's bundles and nothing else: to the middle of the space to the next track (0.15 m at most), from 0.15 m under the
        bundles to the middle of the gap."""
        lo = [o.x_hi for o in self.tracks if o.x_hi <= t.x_lo]; hi = [o.x_lo for o in self.tracks if o.x_lo >= t.x_hi]
        a = t.x_lo - min(0.15, (t.x_lo - max(lo)) / 2 if lo else 0.15, t.x_lo - self.win[0] - 2 * self.h)
        b = t.x_hi + min(0.15, (min(hi) - t.x_hi) / 2 if hi else 0.15)
        return (a, b, Y_TOP - t.bh - 0.15, GAP / 2)

    def solve(self, NI, linear=False, coarse_only=False):
        """NI: {track name: ampere-turns of each leg}.  Returns the window grid and the summary."""
        x_lo, x_hi, y_lo, y_hi = self.box
        if self.outer:
            a_, b_, c_, d_ = self.outer
            go = ms.Grid(b_ - a_, d_ - c_, self.h_outer, x0=a_, y0=c_)
            self._build(go, NI, linear)
            if linear:
                go.solve()
            else:
                go.solve_nonlinear(force_box=self.plate_box(go), max_outer=8)       # the outer box only gives the coarse box its edge values
            g0 = Win(go, self.h0, window=self.box)
        else:
            g0 = ms.Grid(x_hi - x_lo, y_hi - y_lo, self.h0, x0=x_lo, y0=y_lo)
        self._build(g0, NI, linear)
        if linear:
            g0.solve()
        else:
            g0.solve_nonlinear(force_box=self.plate_box(g0))
        if coarse_only:
            self.solves += 1
            return g0, self.summary(g0, NI)
        w = Win(g0, self.h, window=self.win)
        self._build(w, NI, linear)
        if linear:
            w.solve()
        else:
            w.solve_nonlinear(force_box=self.plate_box(w))
        self.solves += 1
        return w, self.summary(w, NI)

    def summary(self, g, NI):
        out = dict(tracks={})
        pm = self.plate.mask(g)
        Bm = np.hypot(g.Bx, g.By)
        for t in self.tracks:
            A_legs = []
            for (xs, ys) in t.lat:
                i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
                A_legs.append(float(g.A[i, j].mean()))
            ni = NI.get(t.name, [0.0] * len(t.legs))
            zone = pm & (g.X >= t.x_lo - 0.1) & (g.X < t.x_hi + 0.1)
            out["tracks"][t.name] = dict(F=-float(g.force_on_box(*self.coil_box(t))), A_legs=A_legs, phi=[sum(s * A_legs[i] for i, s in loop) for loop in t.loops],
                                         W=0.5 * sum(a * A for a, A in zip(ni, A_legs)),
                                         sat=float(np.mean(g.mu_r[zone] < 0.1 * MU)) if zone.any() else 0.0, B_plate=float(Bm[zone].max()) if zone.any() else 0.0)
        out["F_plate"] = float(g.force_on_box(*self.plate_box(g)))
        out["sat"] = float(np.mean(g.mu_r[pm] < 0.1 * MU)) if pm.any() else 0.0
        out["outer"] = (getattr(g, "outer_iters", 0), getattr(g, "outer_change", 0.0), getattr(g, "outer_dF", 0.0))
        return out


def plate_flux(g, plate):
    """The flux that runs sideways through the plate's section at each column of the grid, Wb per meter of ring: A above the plate less A on its face."""
    T = plate.T_at(g.x)
    j0 = int(round((GAP - g.y0) / g.h)); jt = np.minimum(j0 + np.round(T / g.h).astype(int), g.ny - 1)
    return g.x.copy(), g.A[np.arange(g.nx), jt] - g.A[:, j0]


def tailor(x, phi, B_fill, t_min, hp=HP, t_max=0.40):
    """A plate sized from the flux it carries: at each column the thickness that carries the sideways flux at B_fill, rounded up to a whole column
    step, never under t_min, and cut off outside the last columns that need more than t_min."""
    a = math.floor(x[0] / hp + 0.5) * hp; n = int(round((x[-1] - a) / hp))
    c = a + (np.arange(n) + 0.5) * hp
    need = np.abs(np.interp(c, x, phi)) / B_fill
    T = np.minimum(np.maximum(np.ceil(need / hp - 1e-6) * hp, t_min), t_max)
    k = np.nonzero(need > t_min)[0]
    if len(k) == 0:
        return Plate(a, T * 0.0, hp)
    T[:k[0]] = 0.0; T[k[-1] + 1:] = 0.0
    return Plate(a + k[0] * hp, T[k[0]:k[-1] + 1], hp)


# =============================================================================
# one track on its own plate: the requirement, the worst tape, the tape, the rating
# =============================================================================

def worst(g, t, NI_legs, theva, theva77):
    """The worst tape's limit, field and angle: the support as bearing_electrical.worst_tape, the fast track with each stack laid the better way
    (bearing_round6b.worst_ic, 'best')."""
    if t.kind == "support":
        wt, _ = be.worst_tape(g, t, NI_legs, 0.012, theva, theva77)
        return wt[0], wt[1], wt[2]
    return r6.worst_ic([bo.field_at_points(g, xs, ys) for xs, ys in t.lat], "best", theva)


def single_face(t, plate, h, y_hi=0.65, outer=False):
    wide = t.span > 1.6
    xb = max(2.0, math.ceil((t.span / 2 + 1.4) * 10) / 10); yb = 2.0 if wide else 1.5
    xw = max(1.3, math.ceil((max(t.span / 2, plate.x1, -plate.x0) + 0.3) * 10) / 10)
    return Face([t], plate, (-xb, xb, -yb, yb), (-xw, xw, -0.45, y_hi), h=h, outer=(-xb - 3.0, xb + 3.0, -4.0, 4.0) if outer else None)


def find_NI(solve_F, F_target, NI1, tol=0.01, n_max=5):
    """bearing_round7.find_NI (a secant in the logs) with one more step taken on paper: when a solve is within tol of the target, the ampere-turns
    are moved to the target along the local power law of the last two solves.  Returns the ampere-turns at the target, the ampere-turns and pull of
    the last solve (within tol of the target; its fields are used as they are) and that solve."""
    Fa, pa = solve_F(NI1); NIa = NI1
    NIb = NI1 * math.sqrt(F_target / Fa); Fb, pb = solve_F(NIb); n = 0
    while abs(Fb - F_target) / F_target > tol and n < n_max:
        NIn = math.exp(math.log(NIb) + (math.log(F_target) - math.log(Fb)) * (math.log(NIb) - math.log(NIa)) / (math.log(Fb) - math.log(Fa)))
        NIa, Fa = NIb, Fb; Fb, pb = solve_F(NIn); NIb = NIn; n += 1
    pw = (math.log(Fb) - math.log(Fa)) / (math.log(NIb) - math.log(NIa)) if abs(NIb / NIa - 1) > 1e-4 else 2.0
    pw = min(max(pw, 1.0), 2.2)
    return NIb * (F_target / Fb) ** (1 / pw), NIb, Fb, pb


def evaluate(face, t, F_req, NI_guess, theva, theva77, normal=False, levels=(1.25, 1.5)):
    """The track at its requirement: the ampere-turns per leg, the worst tape, the tape at 50 percent, the rating (the worst tape at 70 percent),
    the stored energy; for a support track also normal operation."""
    def solveF(NI):
        g, s = face.solve({t.name: t.design(NI)})
        return s["tracks"][t.name]["F"], (g, s)
    NI, NI_s, F_s, (g, s) = find_NI(solveF, F_req, NI_guess)
    st = s["tracks"][t.name]; F = F_req
    ic, B, th = worst(g, t, t.design(NI_s), theva, theva77)
    tapes = NI / (0.5 * ic); tape = len(t.legs) * tapes * t.w_tape * bd.TAPE_T * bd.TAPE_RHO
    x, phi = plate_flux(g, face.plate)
    st = dict(st, W=st["W"] * (NI / NI_s) ** 2)
    out = dict(NI=NI, F=F, NI_solved=NI_s, F_solved=F_s, F_plate=s["F_plate"], ic=ic, B=B, th=th, tapes_leg=tapes, tape=tape, W=st["W"], phi=st["phi"], sat=st["sat"], B_plate=st["B_plate"],
               fill=tapes * t.w_tape * bd.TAPE_T / (t.bw * t.bh), outer=s["outer"], flux_max=float(np.abs(phi).max()))
    NIs, Fs, ics = [NI], [F], [ic]
    for k in levels:
        g2, s2 = face.solve({t.name: t.design(k * NI)}); i2 = worst(g2, t, t.design(k * NI), theva, theva77)[0]
        NIs.append(k * NI); Fs.append(s2["tracks"][t.name]["F"]); ics.append(i2)
        if (k * NI / tapes) / i2 >= 0.7:
            break
    NI70, F70, share = gt.rating(NIs, Fs, ics, tapes)
    out["rating"] = dict(NI=NI70, F=F70, margin=F70 / F_req, levels=[dict(NI=n, F=f, ic=i) for n, f, i in zip(NIs, Fs, ics)])
    if t.kind != "support":
        Phi = 2 * st["W"] / NI
        out.update(V=10 * Phi * 25.0 / 10e-3, I=NI / 200.0, E_coil=st["W"] * 25.0 / 20)
    if normal:
        NIn, NIn_s, Fn_s, (g, s) = find_NI(solveF, W_HALF, NI * 0.36); Fn = W_HALF
        icn = worst(g, t, t.design(NIn_s), theva, theva77)[0]
        out["normal"] = dict(NI=NIn, F=Fn, sat=s["tracks"][t.name]["sat"], share=(NIn / tapes) / icn, phi=s["tracks"][t.name]["phi"])
    return out, (x, phi)


def tkey(a):
    p = a["plate"]
    ps = f"rect {p['Wp']:.3f}x{p['T']:.3f}" if p["mode"] == "rect" else f"tailor {p['Bf']:.2f} tmin {p['tmin']:.3f}"
    return f"{a['kind']}|n{a['n']}|{a['Wl']:.3f}|{a['bw']:.3f}x{a['bh']:.4f}|{a['F_req']/1e3:.1f}|{ps}|{a['h']:.4f}" + ("|outer" if a.get("outer") else "")


def gkey(a):
    return f"{a['kind']}|n{a['n']}|{a['Wl']:.3f}|{a['bw']:.3f}x{a['bh']:.4f}|{a['F_req']/1e3:.1f}|{a['h']:.4f}" + ("|outer" if a.get("outer") else "")


def job_group(tasks):
    """Tracks of one geometry alone, each on its own plate.  A task: kind, n, Wl, bw, bh, F_req, h, NI_guess, plate {mode 'rect', Wp, T} or
    {mode 'tailor', Bf, tmin}, do_normal.  For the tailored plates: first the track on a thick plate (300 mm, 0.6 m past the outer bundles) at its
    requirement, once; its sideways flux sizes each plate; then the track on each plate."""
    theva, theva77 = be.load_theva()
    outs = []; thick = None; t_th = 0.0
    for a in tasks:
        t0 = time.time()
        out = dict(a, key=tkey(a))
        try:
            t = Trk("t", a["kind"], 0.0, a["Wl"], a["bw"], a["bh"], n=a["n"])
            p = a["plate"]
            out.update(span=t.span, legs=[(x, s) for x, s in t.legs])
            if p["mode"] == "rect":
                plate = Plate.rect(-p["Wp"] / 2, p["Wp"] / 2, p["T"])
                face = single_face(t, plate, a["h"], y_hi=a.get("y_hi", 0.65), outer=a.get("outer", False))
                res, (x, phi) = evaluate(face, t, a["F_req"], a["NI_guess"], theva, theva77, normal=a.get("do_normal", False))
                out["solves"] = face.solves
            else:
                if thick is None:
                    ref = Plate.rect(-(t.span / 2 + 0.6), t.span / 2 + 0.6, 0.30)
                    fr = single_face(t, ref, a["h"], outer=a.get("outer", False))

                    def solveF(NI):
                        g, s = fr.solve({t.name: t.design(NI)})
                        return s["tracks"][t.name]["F"], (g, s)
                    NI0, NI0s, F0s, (g, s) = find_NI(solveF, a["F_req"], a["NI_guess"])
                    x, phi = plate_flux(g, ref); phi = phi * (NI0 / NI0s)
                    thick = dict(NI=NI0, flux_max=float(np.abs(phi).max()), ic=worst(g, t, t.design(NI0s), theva, theva77)[0], x=x, phi=phi, solves=fr.solves)
                    t_th = time.time() - t0
                out["thick"] = dict(NI=thick["NI"], flux_max=thick["flux_max"], ic=thick["ic"], solves=thick["solves"], seconds=t_th)
                plate = tailor(thick["x"], thick["phi"], p["Bf"], p["tmin"])
                face = single_face(t, plate, a["h"], outer=a.get("outer", False))
                res, (x, phi) = evaluate(face, t, a["F_req"], thick["NI"], theva, theva77, normal=a.get("do_normal", False))
                out["solves"] = face.solves
            out.update(res)
            out["plate_json"] = plate.to_json(); out["plate_mass"] = plate.mass(); out["plate_extent"] = plate.extent(); out["plate_Tmax"] = float(plate.T.max())
            sel = slice(None, None, max(1, len(x) // 400))
            out["flux_x"] = x[sel].tolist(); out["flux"] = phi[sel].tolist()
            out["ok"] = True
        except Exception as e:
            import traceback
            out["ok"] = False; out["error"] = repr(e) + traceback.format_exc()
        out["seconds"] = time.time() - t0
        outs.append(out)
    return outs


def run_pool(tasks, fname, jobs=2):
    import multiprocessing as mp
    done = load(fname) or {}
    todo = [t for t in tasks if not (done.get(tkey(t)) or {}).get("ok")]
    groups = {}
    for t in todo:
        groups.setdefault(gkey(t), []).append(t)
    groups = list(groups.values())
    print(f"  {len(tasks)} points, {len(todo)} to solve in {len(groups)} geometries, {jobs} processes", flush=True)
    t0 = time.time()
    if groups:
        with mp.Pool(min(jobs, len(groups))) as pool:
            for outs in pool.imap_unordered(job_group, groups):
                done = load(fname) or {}
                for res in outs:
                    done[res["key"]] = res
                    if res["ok"]:
                        v = f"{res['NI']/1e3:.1f} kA, tape {res['tape']:.1f} kg/m, rating {res['rating']['margin']:.2f}, plate {res['plate_mass']:.0f} kg/m, sat {res['sat']*100:.0f} %"
                    else:
                        v = "FAILED " + res["error"][:300]
                    print(f"    {res['key']}: {v} ({res['seconds']:.0f} s; {(time.time()-t0)/60:.1f} min)", flush=True)
                save(fname, done)
    return load(fname) or {}


def sup_task(Wl, plate, h, F=None, normal=True, bundle=SUP_BUNDLE, guess=437e3):
    return dict(kind="support", n=2, Wl=Wl, bw=bundle[0], bh=bundle[1], F_req=F or gt.F_cap(), h=h, NI_guess=guess, plate=plate, do_normal=normal)


def fast_task(Wl, plate, h, F=F_FAST, n=1, bundle=(FAST["bw"], FAST["bh"]), guess=452e3):
    return dict(kind="fast", n=n, Wl=Wl, bw=bundle[0], bh=bundle[1], F_req=F, h=h, NI_guess=guess, plate=plate, do_normal=False)


BF_LIST = (1.0, 1.2, 1.4, 1.6, 1.8, 2.2)
T_MIN = 0.010
SUP_WL = (0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
FAST_WL = (0.4, 0.6, 0.8, 1.0, 1.4, 2.0)


def tasks_base(h):
    return [dict(sup_task(SUP["Wl"], dict(mode="rect", Wp=SUP["Wp"], T=SUP["T"]), h), y_hi=0.55),
            dict(fast_task(FAST["Wl"], dict(mode="rect", Wp=FAST["Wp"], T=FAST["T"]), h), y_hi=0.55)]


def tasks_tailor(h, bfs=BF_LIST):
    return ([sup_task(SUP["Wl"], dict(mode="tailor", Bf=b, tmin=T_MIN), h, normal=False) for b in bfs] +
            [fast_task(FAST["Wl"], dict(mode="tailor", Bf=b, tmin=T_MIN), h) for b in bfs])


def tasks_width(h, bf_s, bf_f, wls=SUP_WL, wlf=FAST_WL):
    return ([sup_task(W, dict(mode="tailor", Bf=bf_s, tmin=T_MIN), h, normal=False) for W in wls] +
            [fast_task(W, dict(mode="tailor", Bf=bf_f, tmin=T_MIN), h) for W in wlf])


ALONE = [(1, 0.6, (0.204, 0.100)), (1, 1.0, (0.204, 0.100)), (1, 1.0, (0.300, 0.068)),
         (2, 0.6, (0.204, 0.050)), (2, 1.0, (0.204, 0.050)), (3, 0.6, (0.204, 0.050))]          # loops, loop width, bundle


def tasks_alone(h, bf_f):
    return [fast_task(W, dict(mode="tailor", Bf=bf_f, tmin=T_MIN), h, F=F_ALONE, n=n, bundle=b, guess=452e3 * (1.45 if n == 1 else 1.0 if n == 2 else 0.85)) for n, W, b in ALONE]


def main():
    part = "report"; h = 0.005; jobs = 2; extra = {}
    for a in sys.argv[1:]:
        if a.startswith("--part="): part = a[7:]
        elif a.startswith("--h="): h = float(a[4:])
        elif a.startswith("--jobs="): jobs = int(a[7:])
        elif a.startswith("--") and "=" in a: extra[a[2:a.index("=")]] = a[a.index("=") + 1:]
        elif a in ("-h", "--help"): print(__doc__); return
    os.makedirs(LOGS, exist_ok=True)
    tag = extra.get("tag", "")
    sys.stdout = Tee(os.path.join(LOGS, f"round11_plate_{part}{tag}.txt"), "w" if part == "report" else "a")
    t0 = time.time()
    print(f"bearing_round11_plate.py --part={part} --h={h} {extra}   ({time.strftime('%Y-%m-%d %H:%M')})", flush=True)
    if part == "base":
        run_pool(tasks_base(h), "round11_plate_tracks.json", jobs)
    elif part == "tailor":
        bfs = tuple(float(x) for x in extra["bf"].split(",")) if "bf" in extra else BF_LIST
        run_pool(tasks_tailor(h, bfs), "round11_plate_tracks.json", jobs)
    elif part == "width":
        lst = lambda k, d: tuple(float(x) for x in extra[k].split(",") if x) if k in extra else d       # --wls=0.5,0.6 --wlf=0.8,1.0 --file=...
        run_pool(tasks_width(h, float(extra["bfs"]), float(extra["bff"]), lst("wls", SUP_WL), lst("wlf", FAST_WL)), extra.get("file", "round11_plate_width.json"), jobs)
    elif part == "alone":
        run_pool(tasks_alone(h, float(extra["bff"])), "round11_plate_alone.json", jobs)
    elif part == "design":
        # the support track on its sized plate with normal operation solved too (the half face needs it); --wls=0.4,0.5 --bfs=1.2
        run_pool([sup_task(float(w), dict(mode="tailor", Bf=float(extra["bfs"]), tmin=T_MIN), h, normal=True) for w in extra.get("wls", "0.4").split(",")], "round11_plate_design.json", jobs)
    elif part == "conv":
        # the cell size and the box for one headline: the support and the fast track on plates sized at --bf, at h (2.5 mm) and in a larger box
        bf = float(extra.get("bf", "1.2"))
        tasks = []
        for hh, outer in ((0.0025, False), (0.005, True)):
            tasks.append(dict(sup_task(SUP["Wl"], dict(mode="tailor", Bf=bf, tmin=T_MIN), hh, normal=False), outer=outer))
            tasks.append(dict(fast_task(FAST["Wl"], dict(mode="tailor", Bf=bf, tmin=T_MIN), hh), outer=outer))
        run_pool(tasks, "round11_plate_conv.json", jobs)
    else:
        import bearing_round11_plate_parts as pp
        getattr(pp, "part_" + part)(h=h, jobs=jobs, **extra)
    print(f"\n  DONE ({time.time()-t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
