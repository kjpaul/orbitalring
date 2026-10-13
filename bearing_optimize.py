#!/usr/bin/env python3
"""
Optimization of the levitation bearing's coils and gaps (Claude_Code_instructions_coil_optimization_2026-10-12.md).
Part 1: make the field solution trustworthy: grid convergence of the present design, the tape current limit taken from the
field on every tape of a modelled winding instead of the largest cell field at a bundle corner, and an independent check
of the pull by the method of images.
Book: The Orbital Ring, Volume IV: Chapter 8 (Sections 8.1 to 8.5, Tables 8.2 and 8.3; the optimized numbers are for a later pass).
Run: python bearing_optimize.py --part=1 [--quick]   (about an hour; --quick skips the finest grid)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

Nothing here changes a default of final_bearing.py, bearing_mass.py or magnetostatic_2d.py; the present design is rebuilt
with their functions and solved again on finer grids.

Method (Part 1)
  Grids.  The full 4.0 x 3.0 m box of final_bearing.build_shaped is solved at 5 mm as before (saturating plate, Picard loop).
  A window 2.6 x 1.0 m around the track and the plate is then solved again at 5, 2.5 and 1.25 mm with the vector potential
  on the window's edge taken from the full-box solution (Dirichlet values, bilinear interpolation), its own Picard loop for
  the plate, and a conjugate-gradient solve with a Jacobi preconditioner (the direct LU does not fit in memory at 1.25 mm).
  The 5 mm window must reproduce the 5 mm full box; that is the check of the nesting.
  Winding.  Each bundle is the design's 100 mm x 100 mm envelope with the design's turns spread evenly over it: 12 mm tape
  in 8 columns (96 mm) and N/8 rows, 3 mm tape in 33 columns (99 mm).  The field at the center of every tape is read from
  the fine solution and split into the component through the tape face and the component along it, for both ways the
  tape can be wound (faces vertical, faces horizontal).  The tape's critical current at that field and angle comes from the
  THEVA 70 K data (angles 0, 30, 90 and 120 degrees, 0 = field perpendicular to the face; 120 is used for 60 by symmetry,
  linear interpolation in between [UNVERIFIED: THEVA's angle convention beyond 90 degrees]).  A series winding carries one
  current, so the worst tape sets the limit.  The two bounding cases (all of the field through the face, all along it) are
  reported beside as reference orientations (they are not bounds: at 1.6 T the THEVA 30 degree curve lies below the 0 degree one).
  Image check.  At 25 kA per leg the plate is far from saturation and the solver is linear; the pull is compared with the
  method of images for an unsaturated plate (image strength (mu - 1)/(mu + 1) with mu = 200, the plate's back face ignored),
  integrating over the bundle cross-sections.
"""
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

import ring_altitude as ring
import bearing_design as bd
import bearing_mass as bm
import compensation_fields as cf
import final_bearing as fb
import magnetostatic_2d as ms

MU0 = ms.MU0
TAPE_T = bd.TAPE_T                    # 0.08 mm conductor [UNVERIFIED], as bearing_design.py
TAPE_RHO = bd.TAPE_RHO
GAP = bd.GAP0                         # 0.100 m mechanical gap
PLATE_T = fb.PLATE_T                  # 0.20 m
BUNDLE = 0.10                         # m, the design's bundle envelope (square)
SPLIT = 0.10                          # m, center-to-center spacing of the two center bundles (they touch edge to edge)
WINDOW = (-1.3, 1.3, -0.45, 0.55)     # m, the fine window around the track and the plate
GRIDS = (0.005, 0.0025, 0.00125)


def hr(t=""):
    fb.hr(t)


# =============================================================================
# the present design, as final_bearing.py builds it
# =============================================================================

def design_NI(A_top=None):
    """The design ampere-turns per leg from the bearing_air cache, as final_bearing Part A and D use them."""
    sup, lift = fb.load_air()
    theva = bd.read_theva()
    # the casing with the bearing's own mass on top, iterated as final_bearing.part_A does ('on top': 12,681 kg/m, 248 kN/m)
    m = ring.M_LOAD_M
    for _ in range(30):
        dm = fb.design_masses(m * ring.G_NET, theva, sup, lift)
        m_new = ring.M_LOAD_M + dm["casing_side"]
        if abs(m_new - m) < 0.5:
            break
        m = m_new
    F_s = dm["F_s"]
    NI_s, _ = fb.curve_at(sup, F_s, "F")
    NI_l, _ = fb.curve_at(lift, fb.F_LIFT_TRACK, "F")
    return dict(F_s=F_s, NI_s=NI_s, NI_l=NI_l, I_s=dm["t_sup"]["I_op"], I_l=dm["t_lift"]["I_op"], m_casing=m, theva=theva, sup=sup, lift=lift)


CASES = {
    "support, center bundles split 100 mm": dict(W=0.5, n_loops=2, split=SPLIT, w_tape=0.012),
    "support, center bundles joined": dict(W=0.5, n_loops=2, split=0.0, w_tape=0.012),
    "lift, one 1 m loop": dict(W=1.0, n_loops=1, split=0.0, w_tape=0.003),
}


def coarse_solve(case, NI, h=0.005):
    """The full-box solve of final_bearing.part_D for this case (saturating plate)."""
    g, boxes, meta = fb.build_shaped("square", case["W"], case["n_loops"], NI, h=h, split=case["split"])
    box = (-meta["plate_w"] / 2 - 0.1, meta["plate_w"] / 2 + 0.1, GAP - 2 * h, GAP + PLATE_T + 2 * h)
    g.solve_nonlinear(force_box=box)
    return g, boxes, meta, box


# =============================================================================
# a window grid with Dirichlet edge values from a coarser solution
# =============================================================================

class Window(ms.Grid):
    """A magnetostatic_2d.Grid whose edge values of A are prescribed (from a coarser solution) instead of zero."""

    def __init__(self, parent, h, window=WINDOW):
        x_lo, x_hi, y_lo, y_hi = window
        super().__init__(x_hi - x_lo, y_hi - y_lo, h, x0=x_lo, y0=y_lo)
        from scipy.interpolate import RegularGridInterpolator
        f = RegularGridInterpolator((parent.x, parent.y), parent.A, bounds_error=False, fill_value=0.0)
        self.A_edge = np.zeros((self.nx, self.ny))
        edge = np.zeros((self.nx, self.ny), dtype=bool)
        edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
        pts = np.column_stack([self.X[edge], self.Y[edge]])
        self.A_edge[edge] = f(pts)
        self.edge = edge
        self.cg_info = []

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
            b[~ok] += coef[ii[~ok], jj[~ok]] * self.A_edge[i2[~ok], j2[~ok]]      # the prescribed neighbours go to the right-hand side
        M = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_in, n_in))
        A = np.array(self.A_edge, copy=True)
        if n_in <= 300_000:
            a = spla.spsolve(M.tocsc(), b)
            info = "direct"
        else:
            x0 = getattr(self, "A", None)
            x0 = x0[ii, jj] if x0 is not None else np.zeros(n_in)
            Dinv = 1.0 / M.diagonal()
            Pre = spla.LinearOperator(M.shape, matvec=lambda v: Dinv * v)
            count = [0]
            def cb(_): count[0] += 1
            a, flag = spla.cg(M, b, x0=x0, rtol=1e-9, maxiter=40000, M=Pre, callback=cb)
            info = f"cg {count[0]} iterations, flag {flag}"
        A[ii, jj] = a
        self.A = A
        self.Bx = np.zeros((nx, ny)); self.By = np.zeros((nx, ny))
        self.Bx[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
        self.By[1:-1, :] = -(A[2:, :] - A[:-2, :]) / (2 * h)
        self.iters = 1
        self.solve_time = time.time() - t0
        self.cg_info.append((n_in, info, self.solve_time))
        return A


def window_solve(parent, h, case, NI, linear=False):
    """Rebuild the case's conductors and plate on a window grid and solve (saturating plate unless linear)."""
    w = Window(parent, h)
    legs = []
    x0 = -case["n_loops"] * case["W"] / 2
    for k in range(case["n_loops"]):
        sgn = 1 if k % 2 == 0 else -1
        off = (k - (case["n_loops"] - 1) / 2) * case["split"]
        legs.append((x0 + k * case["W"] + off, sgn * NI)); legs.append((x0 + (k + 1) * case["W"] + off, -sgn * NI))
    merged = {}
    for x, a in legs:
        merged[round(x, 6)] = merged.get(round(x, 6), 0.0) + a
    legs = sorted(merged.items())
    plate_w = case["n_loops"] * case["W"] + 0.4 + (case["n_loops"] - 1) * case["split"]
    if linear:
        w.add_iron(-plate_w / 2, plate_w / 2, GAP, GAP + PLATE_T, cf.PLATE_MU)
    else:
        w.add_iron(-plate_w / 2, plate_w / 2, GAP, GAP + PLATE_T, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
    boxes = []
    for x, a in legs:
        w.add_current(x - BUNDLE / 2, x + BUNDLE / 2, -BUNDLE, 0.0, a); boxes.append((x - BUNDLE / 2, x + BUNDLE / 2, -BUNDLE, 0.0, a))
    box = (-plate_w / 2 - 0.1, plate_w / 2 + 0.1, GAP - 2 * h, GAP + PLATE_T + 2 * h)
    t0 = time.time()
    if linear:
        w.solve()
    else:
        w.solve_nonlinear(force_box=box)
    w.wall = time.time() - t0
    return w, boxes, box, legs, plate_w


# =============================================================================
# the field on each tape of a modelled winding
# =============================================================================

def tape_centers(x_lo, x_hi, y_lo, y_hi, w_tape, n_turns):
    """Centers of n_turns tapes of width w_tape spread evenly over the envelope: columns across x, rows up y."""
    n_col = int((x_hi - x_lo) // w_tape)
    n_row = int(math.ceil(n_turns / n_col))
    xs = x_lo + (x_hi - x_lo - n_col * w_tape) / 2 + w_tape * (np.arange(n_col) + 0.5)
    pitch = (y_hi - y_lo) / n_row
    ys = y_lo + pitch * (np.arange(n_row) + 0.5)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X.ravel(), Y.ravel(), n_col, n_row, pitch


def field_at_points(g, xs, ys):
    """Bilinear interpolation of Bx, By at points."""
    fx = (xs - g.x0) / g.h; fy = (ys - g.y0) / g.h
    i = np.clip(np.floor(fx).astype(int), 0, g.nx - 2); j = np.clip(np.floor(fy).astype(int), 0, g.ny - 2)
    tx = fx - i; ty = fy - j
    def lerp(F):
        return ((1 - tx) * (1 - ty) * F[i, j] + tx * (1 - ty) * F[i + 1, j] + (1 - tx) * ty * F[i, j + 1] + tx * ty * F[i + 1, j + 1])
    return lerp(g.Bx), lerp(g.By)


def ic_angle(theva, B, theta_deg, w_tape):
    """Critical current (A) of a tape of width w_tape at field B at angle theta (degrees) between the field and the tape's
    normal: 0 = through the face, 90 = along it.  THEVA curves at 0, 30, 90 and 120 (used for 60 by symmetry), linear in angle."""
    angles = [0.0, 30.0, 60.0, 90.0]
    curves = [theva[0], theva[30], theva[120], theva[90]]
    th = np.abs(np.asarray(theta_deg, float))
    th = np.where(th > 90.0, 180.0 - th, th)
    ic = np.array([np.interp(B, *c) for c in curves])      # 4 x n
    k = np.clip(np.searchsorted(angles, th, side="right") - 1, 0, 2)
    t = (th - np.array(angles)[k]) / 30.0
    rows = np.arange(len(np.atleast_1d(th)))
    return (ic[k, rows] * (1 - t) + ic[k + 1, rows] * t) * w_tape * 100.0


def winding_limits(g, boxes, w_tape, I_op, theva):
    """For every bundle: the field on each tape of the modelled winding and the critical current the worst tape allows,
    for faces vertical and faces horizontal, plus the two bounding cases."""
    out = []
    for (x_lo, x_hi, y_lo, y_hi, a) in boxes:
        n_turns = int(round(abs(a) / I_op))
        xs, ys, n_col, n_row, pitch = tape_centers(x_lo, x_hi, y_lo, y_hi, w_tape, n_turns)
        Bx, By = field_at_points(g, xs, ys)
        Bm = np.hypot(Bx, By)
        # faces vertical (tape width along y): through = Bx, along = By; faces horizontal: through = By
        th_v = np.degrees(np.arctan2(np.abs(By), np.abs(Bx)))
        th_h = np.degrees(np.arctan2(np.abs(Bx), np.abs(By)))
        ic_v = ic_angle(theva, Bm, th_v, w_tape); ic_h = ic_angle(theva, Bm, th_h, w_tape)
        ic_through = ic_angle(theva, Bm, np.zeros_like(Bm), w_tape); ic_along = ic_angle(theva, Bm, np.full_like(Bm, 90.0), w_tape)
        k_v = int(ic_v.argmin()); k_h = int(ic_h.argmin())
        out.append(dict(x=(x_lo + x_hi) / 2, NI=a, n_turns=n_turns, n_col=n_col, n_row=n_row, pitch=pitch,
                        fill=n_turns * w_tape * TAPE_T / ((x_hi - x_lo) * (y_hi - y_lo)),
                        B_mean=float(Bm.mean()), B_max=float(Bm.max()), Bx_max=float(np.abs(Bx).max()), By_max=float(np.abs(By).max()),
                        ic_v=float(ic_v.min()), ic_h=float(ic_h.min()), ic_through=float(ic_through.min()), ic_along=float(ic_along.min()),
                        ic_v_mean=float(ic_v.mean()), ic_h_mean=float(ic_h.mean()),
                        where_v=(float(xs[k_v]), float(ys[k_v]), float(Bm[k_v]), float(th_v[k_v])),
                        where_h=(float(xs[k_h]), float(ys[k_h]), float(Bm[k_h]), float(th_h[k_h]))))
    return out


def bundle_cells(g, boxes):
    """Mean and largest cell field over the conductor cells, as final_bearing Part D reads them."""
    Bmax = 0.0; s = 0.0; n = 0
    for b in boxes:
        m = g.box(*b[:4]) & (g.J != 0)
        Bm = np.sqrt(g.Bx[m]**2 + g.By[m]**2)
        Bmax = max(Bmax, float(Bm.max())); s += float(Bm.sum()); n += int(m.sum())
    return s / max(n, 1), Bmax


# =============================================================================
# the image method
# =============================================================================

def image_pull(legs, NI_scale=1.0, mu=cf.PLATE_MU, n_sub=20):
    """Pull per meter on bundles of uniform current density below an unsaturated plate, method of images: every conductor
    element attracts the image (strength (mu-1)/(mu+1)) of every element, mirrored in the plate's lower face at y = GAP.
    The plate's back face (0.2 m further) is ignored, which makes this an upper bound by a small amount."""
    k_img = (mu - 1.0) / (mu + 1.0)
    pts = []
    for x, a in legs:
        xs = x - BUNDLE / 2 + BUNDLE * (np.arange(n_sub) + 0.5) / n_sub
        ys = -BUNDLE + BUNDLE * (np.arange(n_sub) + 0.5) / n_sub
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        pts.append(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, a * NI_scale / X.size)]))
    P = np.vstack(pts)
    xi, yi, Ii = P[:, 0], P[:, 1], P[:, 2]
    yimg = 2 * GAP - yi                                      # mirror in y = GAP
    Fy = 0.0
    for n in range(0, len(P), 2000):
        x = xi[n:n + 2000, None]; y = yi[n:n + 2000, None]; I = Ii[n:n + 2000, None]
        dx = x - xi[None, :]; dy = y - yimg[None, :]
        r2 = dx * dx + dy * dy
        # force on element (x, y) from image current I_j at (xj, yimg_j): F = mu0 I Ij / (2 pi r) toward the image if parallel
        Fy += float(np.sum(MU0 * I * (k_img * Ii[None, :]) / (2 * math.pi) * dy / r2))
    return -Fy        # dy < 0 means the image is above: the force is upward (toward the plate); return as a positive pull


# =============================================================================
# Part 1
# =============================================================================

def part_1(quick=False):
    hr("PART 1: THE FIELD SOLUTION MADE TRUSTWORTHY")
    d = design_NI()
    theva = d["theva"]
    grids = GRIDS[:-1] if quick else GRIDS
    print(f"  Present design (final_bearing.py Part A and D, bearing_air.json): support track two 0.5 m loops, 100 mm square bundles, 12 mm tape,")
    print(f"  {d['NI_s']/1e3:.0f} kA per leg at {d['I_s']:.0f} A ({d['NI_s']/d['I_s']:,.0f} turns per bundle) for {d['F_s']/1e3:.0f} kN/m (casing {d['m_casing']:,.0f} kg/m with the bearing on top, as Part A);")
    print(f"  lift track one 1 m loop, 3 mm tape, {d['NI_l']/1e3:.0f} kA per leg")
    print(f"  at {d['I_l']:.0f} A ({d['NI_l']/d['I_l']:,.0f} turns).  Gap {GAP*1e3:.0f} mm, plate {PLATE_T*1e3:.0f} mm thick, saturating at {cf.PLATE_B_SAT} T (mu {cf.PLATE_MU:.0f}).")
    print(f"  Center bundles of the support track {SPLIT*1e3:.0f} mm apart center to center (the 100 mm bundles touch edge to edge).")
    print(f"  Winding model: tape {TAPE_T*1e3:.2f} mm thick [UNVERIFIED], turns spread evenly over the 100 x 100 mm envelope; the fill factor is the")
    print(f"  conductor area over the envelope area.  Window {WINDOW[0]} to {WINDOW[1]} m across, {WINDOW[2]} to {WINDOW[3]} m up, edge values from the 5 mm full box.")
    results = {}
    for label, case in CASES.items():
        NI = d["NI_s"] if label.startswith("support") else d["NI_l"]
        I_op = d["I_s"] if label.startswith("support") else d["I_l"]
        w_tape = case["w_tape"]
        print()
        print(f"  1.1  {label}: {NI/1e3:.0f} kA per leg, {w_tape*1e3:.0f} mm tape at {I_op:.0f} A")
        t0 = time.time()
        g0, boxes0, meta0, box0 = coarse_solve(case, NI)
        F0 = g0.force_on_box(*box0)
        Bmean0, Bmax0 = bundle_cells(g0, boxes0)
        print(f"       full box 4.0 x 3.0 m at 5 mm (as final_bearing): pull {F0/1e3:7.1f} kN/m, bundle mean {Bmean0:.3f} T, largest cell {Bmax0:.3f} T ({time.time()-t0:.0f} s)")
        print(f"       {'cells':>8} {'pull kN/m':>10} {'mean T':>7} {'cell max T':>10} {'Ic worst tape, faces vertical':>30} {'faces horizontal':>17} {'all through':>12} {'all along':>10} {'solve s':>8}")
        rows = []
        for h in grids:
            w, boxes, box, legs, plate_w = window_solve(g0, h, case, NI)
            F = w.force_on_box(*box)
            Bmean, Bmax = bundle_cells(w, boxes)
            lim = winding_limits(w, boxes, w_tape, I_op, theva)
            worst_v = min(l["ic_v"] for l in lim); worst_h = min(l["ic_h"] for l in lim)
            worst_t = min(l["ic_through"] for l in lim); worst_a = min(l["ic_along"] for l in lim)
            rows.append(dict(h=h, F=F, Bmean=Bmean, Bmax=Bmax, ic_v=worst_v, ic_h=worst_h, ic_t=worst_t, ic_a=worst_a, lim=lim, legs=legs, plate_w=plate_w, w=w, boxes=boxes, box=box))
            cg = [c for c in w.cg_info if c[1].startswith("cg")]
            cg_note = f", {len(cg)} CG solves of {cg[0][0]:,} unknowns, {sum(int(c[1].split()[1]) for c in cg)} iterations in all" if cg else f", {w.cg_info[0][0]:,} unknowns direct"
            print(f"       {h*1e3:5.2f} mm {F/1e3:10.1f} {Bmean:7.3f} {Bmax:10.3f} {worst_v:30.0f} {worst_h:17.0f} {worst_t:12.0f} {worst_a:10.0f} {w.wall:8.0f}   (Picard {w.outer_iters}, dF {w.outer_dF:.4f}{cg_note})", flush=True)
        # convergence: Richardson on the two finest
        if len(rows) >= 2:
            a, b = rows[-2], rows[-1]
            for key, name, unit, fmt in (("F", "pull", "kN/m", 1e3), ("Bmean", "mean field", "T", 1), ("Bmax", "largest cell field", "T", 1), ("ic_v", "worst-tape Ic, faces vertical", "A", 1), ("ic_h", "worst-tape Ic, faces horizontal", "A", 1)):
                ch = (b[key] - a[key]) / max(abs(b[key]), 1e-30) * 100
                ext = b[key] + (b[key] - a[key])           # first-order Richardson (cell edges: the force and the tape fields converge about linearly)
                print(f"       {name:34}: {a[key]/fmt:9.3f} -> {b[key]/fmt:9.3f} {unit} ({ch:+.2f} % from {a['h']*1e3:.2f} to {b['h']*1e3:.2f} mm; extrapolated {ext/fmt:.3f})")
        # the worst tape of the finest grid
        fine = rows[-1]
        print(f"       worst tapes on the {fine['h']*1e3:.2f} mm grid (x across the track, y below the plate face at 0):")
        for l in fine["lim"]:
            xv, yv, Bv, tv = l["where_v"]; xh, yh, Bh, th = l["where_h"]
            print(f"         bundle at x = {l['x']:+.3f} m, {abs(l['NI'])/1e3:.0f} kA, {l['n_turns']:,} turns in {l['n_col']} columns x {l['n_row']} rows (pitch {l['pitch']*1e3:.2f} mm, fill {l['fill']*100:.0f} %):"
                  f" mean {l['B_mean']:.3f} T, largest {l['B_max']:.3f} T;")
            print(f"           faces vertical:   worst Ic {l['ic_v']:5.0f} A (mean over the tapes {l['ic_v_mean']:5.0f}) at x = {xv:+.4f}, y = {yv:+.4f} m, {Bv:.3f} T at {tv:4.1f} deg from the normal")
            print(f"           faces horizontal: worst Ic {l['ic_h']:5.0f} A (mean over the tapes {l['ic_h_mean']:5.0f}) at x = {xh:+.4f}, y = {yh:+.4f} m, {Bh:.3f} T at {th:4.1f} deg from the normal")
        results[label] = dict(NI=NI, I_op=I_op, w_tape=w_tape, rows=rows, F0=F0, Bmean0=Bmean0, Bmax0=Bmax0)
    # 1.3 image check on the support track at 25 kA (linear)
    print()
    print("  1.3  Independent check of the pull: method of images against the linear solver at 25 kA per leg (plate unsaturated, mu 200).")
    case = CASES["support, center bundles split 100 mm"]
    NI_lin = 25e3
    g0, boxes0, meta0, box0 = coarse_solve(case, NI_lin)        # the Picard loop stays linear at this current
    for h in (0.005, 0.0025):
        w, boxes, box, legs, plate_w = window_solve(g0, h, case, NI_lin, linear=True)
        F = w.force_on_box(*box)
        F_img = image_pull(legs)
        print(f"       {h*1e3:5.2f} mm window, linear plate: solver {F:8.1f} N/m; images {F_img:8.1f} N/m; difference {(F/F_img-1)*100:+.2f} % (images ignore the plate's back face and its 1.4 m width)")
    # restated baseline
    print()
    hr("THE PRESENT DESIGN ON THE CONVERGED BASIS")
    for label, r in results.items():
        fine = r["rows"][-1]; coarse = r["rows"][0]
        I50_v = 0.5 * fine["ic_v"]; I50_h = 0.5 * fine["ic_h"]
        print(f"  {label}: {r['NI']/1e3:.0f} kA per leg")
        print(f"    pull {fine['F']/1e3:.1f} kN/m on the {fine['h']*1e3:.2f} mm grid ({r['F0']/1e3:.1f} on the 5 mm full box, {(fine['F']/r['F0']-1)*100:+.1f} %)")
        print(f"    bundle mean field {fine['Bmean']:.3f} T ({r['Bmean0']:.3f} at 5 mm); largest cell field {fine['Bmax']:.2f} T ({r['Bmax0']:.2f} at 5 mm): the cell maximum does not converge, the tape fields do")
        print(f"    tape limit from the worst tape: faces vertical {fine['ic_v']:.0f} A, faces horizontal {fine['ic_h']:.0f} A; with the field taken all through the face {fine['ic_t']:.0f} A, all along it {fine['ic_a']:.0f} A")
        print(f"    at 50 percent of the limit: {I50_v:.0f} A (vertical) or {I50_h:.0f} A (horizontal) against the design's {r['I_op']:.0f} A; the turns the design needs at {I50_v:.0f} A: {r['NI']/I50_v:,.0f} per bundle,")
        n_legs = 4 if label.startswith("support") else 2
        m_v = r["NI"] / I50_v * n_legs * TAPE_T * r["w_tape"] * TAPE_RHO
        m_h = r["NI"] / I50_h * n_legs * TAPE_T * r["w_tape"] * TAPE_RHO
        print(f"    tape per meter of track at 50 percent: {m_v:.1f} kg/m (faces vertical), {m_h:.1f} kg/m (faces horizontal), counting {n_legs} leg-equivalents of {r['w_tape']*1e3:.0f} mm tape")
    print()
    print("  Book figures of 2026-10-12 to compare with (Chapter 8): corner field 2.3 T joined and 1.7 T split; 245 A limit at the corner (split), 122 A at 50 percent;")
    print("  field 0.92 T average over the bundle at 322,000 ampere-turns; support track 250.6 kN/m at 322 kA per leg (bearing_air.json, 10 mm cells).")
    return results


def main():
    part = "1"; quick = False
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--part="):
            part = arg.split("=")[1]
        elif arg == "--quick":
            quick = True
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    hr("BEARING COIL AND GAP OPTIMIZATION (2026-10-13)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km: casing {ring.M_LOAD_M:,.0f} kg/m ({ring.M_LOAD_M*ring.G_NET/1e3:.1f} kN/m), cable hardware {ring.M_HARDWARE_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m, face {ring.CABLE_SIDE:.2f} m.")
    if "1" in part:
        part_1(quick)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
