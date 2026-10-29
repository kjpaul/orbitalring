#!/usr/bin/env python3
"""
Optimization of the levitation bearing's coils and gaps (Claude_Code_instructions_coil_optimization_2026-10-12.md).
Part 1: make the field solution trustworthy: grid convergence of the present design, the tape current limit taken from the
field on every tape of a modelled winding instead of the largest cell field at a bundle corner, and an independent check
of the pull by the method of images.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (Sections 7.5 to 7.9, Tables 7.2 and 7.3; the optimized numbers are for a later pass).
Run: python bearing_optimize.py --part=1 [--quick]   (about two hours; --quick skips the finest grid)
     python bearing_optimize.py --part=2 --stage=screen|confirm25|confirm125|report   (Part 2, see below)
     python bearing_optimize.py --part=2b --stage=rerate|rerate25|wall   (Part 2 redone on Paul's rules of 2026-10-21, and the cryostat wall)
     python bearing_optimize.py --part=3 --stage=lift|lift25|report3   (Part 3: lift tracks and the load winding, coil length)
     python bearing_optimize.py --part=3b --stage=rerun|report   (Part 3b, Paul's decisions of 2026-10-23: the Design Case requirement,
                                the full tape angle data, the 70 percent summary, the composite wall at 0.4 bar, the control plan)
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

Part 2 (support track search; results cached in <output folder>/bearing_optimize_part2_<cells>.json, so a stage can be rerun)
  --stage=screen      5 mm window: loops, widths and spacings with square 100 mm bundles and a 20 mm wall; then the bundle
                      cross-sections on the four best geometries and the present one; then the wall allowance; then the
                      one-at-a-time neighbours of the best.  About four hours.
  --stage=confirm25   the eight best of the screen and the present design at 2.5 mm (about four hours).
  --stage=confirm125  the two best at 1.25 mm, requirement point only (about three hours).
  --stage=report      the tables from the caches, without solving.
  Each case: a 10 mm full box solved once (saturating plate) for the window's edge values, which then scale with the
  ampere-turns; a linear solve for the linear pull; nonlinear solves at the requirement (two, interpolated in log-log),
  at 1.13 times the requirement (the rating check) and at the load winding's 198 kN/m.  The tape limit is the worst tape of
  the modelled winding (Part 1), faces vertical or horizontal, whichever is lower; one current per track.
"""
import json
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
import levitation_config as lc
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


def read_robinson_angles(path=None):
    """The full 70 K angle dependence of the Theva AP tape from the Robinson dataset (data/THEVA_Robinson_70K_angle_dependence.csv,
    made by tools/extract_theva_robinson_angles.py): a dict field (T) -> (angles deg, Ic/w A/cm), angles 0 to 240 degrees."""
    import csv
    path = path or os.path.join(ring.DATA_DIR, "THEVA_Robinson_70K_angle_dependence.csv")
    data = {}
    with open(path, encoding="utf-8-sig") as f:
        rd = csv.reader(f); next(rd)
        for B, a, i in rd:
            data.setdefault(float(B), []).append((float(a), float(i)))
    return {B: (np.array([q[0] for q in sorted(v)]), np.array([q[1] for q in sorted(v)])) for B, v in data.items()}


def ic_angle_full(rob, B, theta_deg, w_tape):
    """Critical current (A) from the full Robinson angle data: at each measured field the Ic/w at the tape's angle, taking the
    lower of theta and 180 - theta (the tape's turn in the bundle is not known, and this tape's Ic is not symmetric about
    90 degrees: at 1 T it is lowest near 45 degrees and highest near 120); linear between the measured fields."""
    Bs = np.array(sorted(rob)); th = np.abs(np.asarray(theta_deg, float)); th = np.where(th > 90.0, 180.0 - th, th)
    B = np.asarray(B, float)
    per_B = np.array([np.minimum(np.interp(th, *rob[b]), np.interp(180.0 - th, *rob[b])) for b in Bs])     # [fields, points]
    Bc = np.clip(B, Bs[0], Bs[-1]); k = np.clip(np.searchsorted(Bs, Bc) - 1, 0, len(Bs) - 2)
    t = (Bc - Bs[k]) / (Bs[k + 1] - Bs[k]); rows = np.arange(len(np.atleast_1d(th)))
    return (per_B[k, rows] * (1 - t) + per_B[k + 1, rows] * t) * w_tape * 100.0


def ic_angle(theva, B, theta_deg, w_tape):
    """Critical current (A) of a tape of width w_tape at field B at angle theta (degrees) between the field and the tape's
    normal: 0 = through the face, 90 = along it.  THEVA curves at 0, 30, 90 and 120 (used for 60 by symmetry), linear in angle.
    When the table carries the full Robinson angle data under the key "full" (Part 3b), that data is used instead."""
    if "full" in theva:
        return ic_angle_full(theva["full"], B, theta_deg, w_tape)
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
    print("  Book figures of 2026-10-12 to compare with (Chapter 7): corner field 2.3 T joined and 1.7 T split; 245 A limit at the corner (split), 122 A at 50 percent;")
    print("  field 0.92 T average over the bundle at 322,000 ampere-turns; support track 250.6 kN/m at 322 kA per leg (bearing_air.json, 10 mm cells).")
    return results

# =============================================================================
# Part 2: support track search
# =============================================================================

PLATE_W = fb.PLATE_W                  # 1.4 m, fixed (rule 7)
I_REF_12, I_REF_3 = 163.0, 43.0       # A, the present design's currents: they set the turn count of the field model only
CRYO = bm.CRYO_FACTOR                 # formers, cryostat and supports as a multiple of the tape mass [UNVERIFIED, as bearing_mass]
F_REQ = None                          # set from design_NI: 248 kN/m
F_RATE_MARGIN = 1.13
F_LOAD = fb.F_LOAD_TRACK              # 198 kN/m, the load compensation winding on the same loops
SAT_LIMIT = 0.12                      # the plate may not be saturated beyond the present 12 percent (task, Part 2); the present computes to 9.7 to 10 percent
COARSE_H = 0.01
SHAPES = {
    "square 100": dict(kind="rect", w=0.10, h=0.10),
    "wide 150x67": dict(kind="rect", w=0.15, h=0.0667),
    "tall 67x150": dict(kind="rect", w=0.0667, h=0.15),
    "wide 150x100": dict(kind="rect", w=0.15, h=0.10),
    "tall 100x150": dict(kind="rect", w=0.10, h=0.15),
    "square 141": dict(kind="rect", w=0.1414, h=0.1414),
    "flat 200x50": dict(kind="rect", w=0.20, h=0.05),
    "round 56": dict(kind="round", R=0.0564),
    "round 80": dict(kind="round", R=0.08),
    "graded 80": dict(kind="graded", R=0.08),
    "graded 100": dict(kind="graded", R=0.10),
}
THEVA_77 = os.path.join(ring.DATA_DIR, "THEVA_Pro-Line_Advanced_Pinning_2G_HTS_77.5K_Field_Dependence.csv")


def read_theva_file(path):
    import csv
    data = {}
    with open(path, encoding="utf-8-sig") as f:
        rd = csv.reader(f); next(rd)
        for row in rd:
            if len(row) < 4 or not row[1]:
                continue
            data.setdefault(int(round(float(row[2]))), []).append((float(row[1]), float(row[3])))
    return {a: (np.array([q[0] for q in sorted(v)]), np.array([q[1] for q in sorted(v)])) for a, v in data.items()}


def shape_width(shape):
    s = SHAPES[shape]
    return s["w"] if s["kind"] == "rect" else 2 * s["R"]


def shape_height(shape):
    s = SHAPES[shape]
    return s["h"] if s["kind"] == "rect" else 2 * s["R"]


def shape_area(shape):
    s = SHAPES[shape]
    return s["w"] * s["h"] if s["kind"] == "rect" else math.pi * s["R"]**2


def geom_key(g):
    sp = "joined" if g["gc"] is None else f"{g['gc']*1e3:.0f}"
    return f"{g['n']} x {g['W']:.2f} m, spacing {sp}, {g['shape']}, wall {g['allow']*1e3:.0f} mm"


def tape_lattice(shape, xc, y_top, w_tape, N):
    """N tape centers filling the bundle whose top face is at y_top, centered on xc.  rect: columns of w_tape, rows evenly
    spread; round: the same lattice clipped to the circle, row pitch set so that N tapes fit; graded: the clipped lattice
    thinned toward the rim with density 1 - (r/R)^2 (a quasi-random threshold), row pitch set so that N remain."""
    s = SHAPES[shape]
    if s["kind"] == "rect":
        xs, ys, n_col, n_row, pitch = tape_centers(xc - s["w"] / 2, xc + s["w"] / 2, y_top - s["h"], y_top, w_tape, N)
        return xs, ys
    R = s["R"]; yc = y_top - R
    n_col = max(1, int((2 * R) // w_tape))
    cols = xc + w_tape * (np.arange(n_col) - (n_col - 1) / 2)
    def lattice(pitch):
        n_row = int(2 * R / pitch)
        rows = yc + pitch * (np.arange(n_row) - (n_row - 1) / 2)
        X, Y = np.meshgrid(cols, rows, indexing="ij")
        I, J = np.meshgrid(np.arange(n_col), np.arange(n_row), indexing="ij")
        r = np.hypot(X - xc, Y - yc)
        keep = r <= R
        if s["kind"] == "graded":
            u = np.mod(0.7548776662 * I + 0.5698402910 * J, 1.0)
            keep &= (1 - (r / R)**2) >= u
        return X[keep], Y[keep]
    lo, hi = TAPE_T, 2 * R
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        n = len(lattice(mid)[0])
        if n > N: lo = mid
        else: hi = mid
    xs, ys = lattice(lo)
    return xs[:N], ys[:N]


def build_geometry(geom, h, NI, grid=None, window_parent=None, w_tape=0.012, I_ref=I_REF_12, linear=False):
    """The track of `geom` with NI ampere-turns per leg on a grid of cell h: either a full box (grid=None, window_parent=None)
    or a Window of window_parent.  The current is deposited tape by tape from the lattice of the field model (I_ref per tape),
    so the field is that of the modelled winding.  Returns the grid, the bundles (x, NI, tape lattice) and the plate box."""
    n, W, gc, shape, allow = geom["n"], geom["W"], geom["gc"], geom["shape"], geom["allow"]
    wb = shape_width(shape)
    split = 0.0 if gc is None else wb + gc                      # center-to-center distance of the two center bundles
    x0 = -(n * W + (n - 1) * split) / 2
    legs = []
    for k in range(n):
        sgn = 1 if k % 2 == 0 else -1
        xl = x0 + k * (W + split); xr = xl + W
        legs.append((xl, sgn * NI)); legs.append((xr, -sgn * NI))
    merged = {}
    for x, a in legs:
        merged[round(x, 6)] = merged.get(round(x, 6), 0.0) + a
    legs = sorted(merged.items())
    if grid is None:
        g = Window(window_parent, h) if window_parent is not None else ms.Grid(4.0, 3.0, h)
    else:
        g = grid
    if linear:
        g.add_iron(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T, cf.PLATE_MU)
    else:
        g.add_iron(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T, cf.PLATE_MU, B_sat=cf.PLATE_B_SAT)
    bundles = []
    y_top = -allow
    for x, a in legs:
        N = max(1, int(round(abs(a) / I_ref)))
        xs, ys = tape_lattice(shape, x, y_top, w_tape, N)
        I = a / len(xs)
        i = np.clip(np.round((xs - g.x0) / g.h).astype(int), 0, g.nx - 1); j = np.clip(np.round((ys - g.y0) / g.h).astype(int), 0, g.ny - 1)
        np.add.at(g.J, (i, j), I / (g.h * g.h))
        bundles.append(dict(x=x, NI=a, xs=xs, ys=ys))
    box = (-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, GAP - 2 * h, GAP + PLATE_T + 2 * h)
    return g, bundles, box, legs, split


def tape_limits(g, bundles, w_tape, theva, theva77=None, scale=1.0):
    """Worst tape of each bundle (faces vertical and horizontal), with the 77.5 K value at the same tape when given.
    scale: the solved ampere-turns over the build ampere-turns of the bundles."""
    out = []
    for b in bundles:
        Bx, By = field_at_points(g, b["xs"], b["ys"])
        Bm = np.hypot(Bx, By)
        th_v = np.degrees(np.arctan2(np.abs(By), np.abs(Bx))); th_h = np.degrees(np.arctan2(np.abs(Bx), np.abs(By)))
        ic_v = ic_angle(theva, Bm, th_v, w_tape); ic_h = ic_angle(theva, Bm, th_h, w_tape)
        k_v = int(ic_v.argmin()); k_h = int(ic_h.argmin())
        d = dict(x=b["x"], NI=b["NI"] * scale, n_model=len(b["xs"]), B_mean=float(Bm.mean()), B_max=float(Bm.max()),
                 ic_v=float(ic_v[k_v]), ic_h=float(ic_h[k_h]), B_v=float(Bm[k_v]), th_v=float(th_v[k_v]), B_h=float(Bm[k_h]), th_h=float(th_h[k_h]),
                 at_v=(float(b["xs"][k_v]), float(b["ys"][k_v])), at_h=(float(b["xs"][k_h]), float(b["ys"][k_h])))
        if theva77 is not None:
            d["ic77_v"] = float(ic_angle(theva77, np.array([Bm[k_v]]), np.array([th_v[k_v]]), w_tape)[0])
            d["ic77_h"] = float(ic_angle(theva77, np.array([Bm[k_h]]), np.array([th_h[k_h]]), w_tape)[0])
        out.append(d)
    return out


def plate_saturation(g):
    m = g.box(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T)
    return float(np.mean(g.mu_r[m] < 0.1 * cf.PLATE_MU)), float(np.sqrt(g.Bx[m]**2 + g.By[m]**2).max())


def stored_energy(g):
    return float(0.5 * np.sum(g.A * g.J) * g.h * g.h)


class Evaluator:
    """Everything Part 2 needs for one geometry at one cell size: the pull curve near the requirement, the rating check,
    the load winding, the tape limits and the masses."""

    def __init__(self, geom, h, theva, theva77, verbose=False):
        self.geom, self.h, self.theva, self.theva77, self.verbose = geom, h, theva, theva77, verbose
        self.solves = 0; self.t0 = time.time()
        # the coarse full box, nonlinear, at an estimate of the required ampere-turns (edge values scale with NI from it)
        self.NI_est = 322e3
        g0, _, box0, _, _ = build_geometry(geom, COARSE_H, self.NI_est)
        g0.solve_nonlinear(force_box=box0)
        self.w, self.bundles, self.box, self.legs, self.split = build_geometry(geom, h, self.NI_est, window_parent=g0)
        self.J_ref = self.w.J.copy(); self.A_edge_ref = self.w.A_edge.copy()
        self.wb = shape_width(geom["shape"]); self.hb = shape_height(geom["shape"])
        self.track_width = geom["n"] * geom["W"] + (geom["n"] - 1) * self.split + self.wb

    def run(self, NI, linear=False):
        s = NI / self.NI_est
        self.w.J = self.J_ref * s; self.w.A_edge = self.A_edge_ref * s
        if linear:
            mu_save = self.w.mu_r.copy(); mats = self.w.materials
            self.w.mu_r[:] = 1.0; self.w.mu_r[self.w.box(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T)] = cf.PLATE_MU
            self.w.materials = []
            self.w.solve()
            F = self.w.force_on_box(*self.box)
            self.w.mu_r = mu_save; self.w.materials = mats
            self.solves += 1
            return F
        self.w.solve_nonlinear(force_box=self.box)
        self.solves += 1
        return self.w.force_on_box(*self.box)

    def curve_point(self, NI):
        F = self.run(NI)
        sat, Bpl = plate_saturation(self.w)
        lim12 = tape_limits(self.w, self.bundles, 0.012, self.theva, self.theva77, scale=NI / self.NI_est)
        return dict(NI=NI, F=F, sat=sat, B_plate=Bpl, W_mag=stored_energy(self.w), eff=F / (self.F_lin1 * (NI / self.NI_est)**2),
                    B_mean=float(np.mean([d["B_mean"] for d in lim12])), lim12=lim12,
                    lim3=tape_limits(self.w, self.bundles, 0.003, self.theva, self.theva77, scale=NI / self.NI_est))

    @staticmethod
    def interp_NI(points, F):
        """NI for a pull F from the solved points: log-log interpolation of F(NI), extrapolated with the local slope."""
        pts = sorted(points, key=lambda p: p["NI"])
        x = np.log([p["NI"] for p in pts]); y = np.log([p["F"] for p in pts]); lf = math.log(F)
        if len(pts) == 1:
            return pts[0]["NI"] * math.sqrt(F / pts[0]["F"])
        if y[0] <= lf <= y[-1]:
            return float(np.exp(np.interp(lf, y, x)))
        a, b = (0, 1) if lf < y[0] else (len(pts) - 2, len(pts) - 1)      # extrapolate with the slope of the nearest pair
        return float(np.exp(x[a] + (lf - y[a]) * (x[b] - x[a]) / (y[b] - y[a])))

    def evaluate(self, F_req=None, s_share=0.5, light_from=None):
        """light_from: a result of a coarser grid; then only the requirement point is solved here (one linear, one nonlinear
        solve) and the rating and load quantities are taken from that result."""
        F_req = F_req or F_REQ
        self.F_lin1 = self.run(self.NI_est, linear=True)
        if light_from is not None:
            return self.evaluate_light(light_from, F_req, s_share)
        pts = []
        # the requirement: two nonlinear solves and a log-log interpolation
        NI1 = self.NI_est * math.sqrt(F_req / self.F_lin1 / 0.6)      # the saturating pull is about 0.6 of the linear one at the requirement
        pts.append(self.curve_point(NI1))
        NI2 = self.interp_NI(pts, F_req)
        if abs(NI2 / NI1 - 1) < 0.004:
            NI2 = NI1 * 1.03
        pts.append(self.curve_point(NI2))
        NI_req = self.interp_NI(pts, F_req)
        p_req = min(pts, key=lambda p: abs(p["NI"] - NI_req))
        if abs(p_req["F"] / F_req - 1) > 0.01:
            pts.append(self.curve_point(NI_req)); p_req = pts[-1]; NI_req = self.interp_NI(pts, F_req)
        # the rating, as bearing_mass rates an air-core track: the pull where the mean field over the bundles reaches
        # B_WIND_LIMIT (1.0 T); one solve at the ampere-turns that scaling predicts, then a local power-law correction
        NI_r = NI_req * bm.B_WIND_LIMIT / p_req["B_mean"]
        p_r = self.curve_point(NI_r); pts.append(p_r)
        NI_e = NI_r * bm.B_WIND_LIMIT / p_r["B_mean"]
        F_rated = math.exp(float(np.interp(math.log(NI_e), np.log([q["NI"] for q in sorted(pts, key=lambda q: q["NI"])]), np.log([q["F"] for q in sorted(pts, key=lambda q: q["NI"])]))))
        # the load winding at 198 kN/m
        NI_l = self.interp_NI(pts, F_LOAD)
        p_l = self.curve_point(NI_l); pts.append(p_l)
        NI_l = self.interp_NI(pts, F_LOAD)
        # scale the solved fields to the interpolated ampere-turns (a few percent at most) for the limits
        def limits_at(p, NI, key):
            out = []
            for d in p[key]:
                sc = NI / p["NI"]
                w_t = 0.012 if key == "lim12" else 0.003
                r = dict(d)
                for face in ("v", "h"):
                    B = d[f"B_{face}"] * sc; th = d[f"th_{face}"]
                    r[f"ic_{face}"] = float(ic_angle(self.theva, np.array([B]), np.array([th]), w_t)[0])
                    r[f"ic77_{face}"] = float(ic_angle(self.theva77, np.array([B]), np.array([th]), w_t)[0])
                    r[f"B_{face}"] = B
                out.append(r)
            return out
        lim12 = limits_at(p_req, NI_req, "lim12")
        lim3 = limits_at(p_l, NI_l, "lim3")
        # the tape: one current per track (its worst tape), at the share s of the limit
        ic12 = min(min(d["ic_v"], d["ic_h"]) for d in lim12); face12 = "v" if min(d["ic_v"] for d in lim12) <= min(d["ic_h"] for d in lim12) else "h"
        ic3 = min(min(d["ic_v"], d["ic_h"]) for d in lim3); face3 = "v" if min(d["ic_v"] for d in lim3) <= min(d["ic_h"] for d in lim3) else "h"
        I12 = s_share * ic12; I3 = s_share * ic3
        turns12 = sum(abs(d["NI"]) for d in lim12) * NI_req / p_req["NI"] / I12
        turns3 = sum(abs(d["NI"]) for d in lim3) * NI_l / p_l["NI"] / I3
        tape12 = turns12 * 0.012 * TAPE_T * TAPE_RHO; tape3 = turns3 * 0.003 * TAPE_T * TAPE_RHO
        fill = (turns12 * 0.012 + turns3 * 0.003) * TAPE_T / (shape_area(self.geom["shape"]) * len(self.bundles))
        W_req = p_req["W_mag"] * (NI_req / p_req["NI"])**2
        W_load = p_l["W_mag"] * (NI_l / p_l["NI"])**2
        E_coil = W_load * fb.L_LOAD / fb.N_FAM
        res = dict(key=geom_key(self.geom), geom=self.geom, h=self.h, NI_req=NI_req, F_req=F_req, F_lin_req=self.F_lin1 * (NI_req / self.NI_est)**2,
                   eff_req=p_req["eff"], sat=p_req["sat"], B_plate=p_req["B_plate"], NI_rate=NI_r, F_rate=p_r["F"], eff_rate=p_r["eff"],
                   rating_ok=F_rated >= F_RATE_MARGIN * F_req, NI_rated=NI_e, F_rated=F_rated, B_mean_req=p_req["B_mean"], NI_load=NI_l, W_req=W_req, W_load=W_load, E_coil=E_coil, VA=2 * E_coil / fb.TAU_HW,
                   ic12=ic12, face12=face12, ic3=ic3, face3=face3, I12=I12, I3=I3, turns12=turns12, turns3=turns3, tape12=tape12, tape3=tape3, fill=fill,
                   track_width=self.track_width, fits=self.track_width <= PLATE_W + 1e-9, lim12=lim12, lim3=lim3, solves=self.solves, seconds=time.time() - self.t0,
                   B_mean=float(np.mean([d["B_mean"] for d in lim12])), split=self.split)
        return res


def _evaluate_light(self, light_from, F_req, s_share):
    p = self.curve_point(light_from["NI_req"])
    NI_req = light_from["NI_req"] * math.sqrt(F_req / p["F"])           # one correction, linear scaling (a fraction of a percent)
    lim12 = []
    for d in p["lim12"]:
        r = dict(d); sc = NI_req / p["NI"]
        for face in ("v", "h"):
            B = d["B_" + face] * sc; th = d["th_" + face]
            r["ic_" + face] = float(ic_angle(self.theva, np.array([B]), np.array([th]), 0.012)[0])
            r["ic77_" + face] = float(ic_angle(self.theva77, np.array([B]), np.array([th]), 0.012)[0])
            r["B_" + face] = B
        lim12.append(r)
    ic12 = min(min(d["ic_v"], d["ic_h"]) for d in lim12); face12 = "v" if min(d["ic_v"] for d in lim12) <= min(d["ic_h"] for d in lim12) else "h"
    I12 = s_share * ic12
    turns12 = sum(abs(d["NI"]) for d in lim12) * NI_req / p["NI"] / I12
    tape12 = turns12 * 0.012 * TAPE_T * TAPE_RHO
    res = dict(light_from)
    res.update(h=self.h, NI_req=NI_req, F_at_coarse_NI=p["F"], B_mean_req=p["B_mean"], F_lin_req=self.F_lin1 * (NI_req / self.NI_est)**2, eff_req=p["eff"], sat=p["sat"], B_plate=p["B_plate"],
               W_req=p["W_mag"] * (NI_req / p["NI"])**2, ic12=ic12, face12=face12, I12=I12, turns12=turns12, tape12=tape12, lim12=lim12,
               fill=(turns12 * 0.012 + light_from["turns3"] * 0.003) * TAPE_T / (shape_area(self.geom["shape"]) * len(self.bundles)),
               B_mean=float(np.mean([d["B_mean"] for d in lim12])), solves=self.solves, seconds=time.time() - self.t0, light=True)
    return res


Evaluator.evaluate_light = _evaluate_light


def casing_mass(res, lift_tape):
    """Casing-side mass per meter for both support tracks (support and load windings) and the two lift tracks, with formers,
    cryostat and supports at CRYO times the tape (as bearing_mass.py)."""
    tape = 2 * (res["tape12"] + res["tape3"]) + 2 * lift_tape
    return tape, tape * (1 + CRYO)


def mass_factors():
    d0 = ring.design(); d1 = ring.design(m_load=ring.M_LOAD_M + 1.0); d2 = ring.design(m_hw=ring.M_HARDWARE_M + 1.0)
    return d1["m_cable"] - d0["m_cable"], d2["m_cable"] - d0["m_cable"]


def temperature_table(lim, face, shares=(0.5, 0.6, 0.7, 0.8), NI_normal_frac=None):
    """The tape current table at the worst tape: current at each share, the warming the tape can take before its limit falls
    to that current (critical current linear in temperature between the 70 K and 77.5 K data [UNVERIFIED]), and the current in
    normal operation (the ampere-turns of the weight share over those of the requirement)."""
    d = min(lim, key=lambda q: q[f"ic_{face}"])
    ic70 = d[f"ic_{face}"]; ic77 = d[f"ic77_{face}"]
    slope = (ic70 - ic77) / 7.5                       # A per kelvin
    rows = []
    for s in shares:
        I = s * ic70
        dT = (ic70 - I) / slope if slope > 0 else float("inf")
        rows.append((s, I, dT, I * NI_normal_frac if NI_normal_frac else None))
    return ic70, ic77, slope, rows


def stray_field_check(geom, NI_normal, NI_rating, NI_lift=0.0):
    """The stray field of all four tracks in the 30 x 30 m linear box of final_bearing Part D.2 (50 mm cells), with the
    candidate's legs under both support plates and the present lift loops above the top plates; returns the field (T) at the
    LIM stator points and 1, 3 and 10 m below, beside and above the ring, for the normal and the rating currents."""
    side = ring.CABLE_SIDE
    out = {}
    for mode, NI_s, NI_l in (("normal", NI_normal, 0.0), ("rating", NI_rating, NI_lift)):
        g = ms.Grid(30.0, 30.0, 0.05, x0=-15.0, y0=-10.0)
        for sx in (-1, 1):
            xc = sx * fb.Y_TRACK
            g.add_iron(xc - PLATE_W / 2, xc + PLATE_W / 2, 0.0, PLATE_T, cf.PLATE_MU)
            g.add_iron(xc - 0.7, xc + 0.7, side - PLATE_T, side, cf.PLATE_MU)
            # the candidate's legs (merged bundles), their tops at -allow below the underside plate's lower face
            n, W, gc, shape, allow = geom["n"], geom["W"], geom["gc"], geom["shape"], geom["allow"]
            wb = shape_width(shape); hb = shape_height(shape); split = 0.0 if gc is None else wb + gc
            x0 = -(n * W + (n - 1) * split) / 2
            legs = {}
            for k in range(n):
                sgn = 1 if k % 2 == 0 else -1
                xl = x0 + k * (W + split)
                for x, a in ((xl, sgn * NI_s), (xl + W, -sgn * NI_s)):
                    legs[round(x, 6)] = legs.get(round(x, 6), 0.0) + a
            for x, a in legs.items():
                g.add_current(xc + x - wb / 2, xc + x + wb / 2, -PLATE_T - allow - hb, -PLATE_T - allow, a)
            if NI_l:
                for dx, a in ((-0.5, NI_l), (0.5, -NI_l)):
                    g.add_current(xc + dx - 0.05, xc + dx + 0.05, side + 0.1, side + 0.2, a)
        g.solve()
        pts = {"LIM stator, center": (fb.STATOR_X, fb.STATOR_Y), "LIM stator, bottom": (fb.STATOR_X, fb.STATOR_Y - 1.0), "LIM stator, top": (fb.STATOR_X, fb.STATOR_Y + 1.0)}
        for d_ in (1.0, 3.0, 10.0):
            pts[f"{d_:.0f} m below the ring"] = (fb.Y_TRACK, -PLATE_T - geom["allow"] - shape_height(geom["shape"]) - d_)
            pts[f"{d_:.0f} m beside the ring"] = (side / 2 + d_, side / 2)
            pts[f"{d_:.0f} m above the ring"] = (fb.Y_TRACK, side + 0.2 + d_)
        res = {k: math.hypot(*g.field_at(x, y)) for k, (x, y) in pts.items()}
        m_st = g.box(fb.STATOR_X - 0.1, fb.STATOR_X + 0.1, fb.STATOR_Y - 1.0, fb.STATOR_Y + 1.0)
        res["largest over the stator"] = float(np.sqrt(g.Bx[m_st]**2 + g.By[m_st]**2).max())
        out[mode] = res
    return out


def part2_cache(h):
    return os.path.join(ring.out_dir(), f"bearing_optimize_part2_{h*1e3:g}mm.json")


def load_cache(h):
    p = part2_cache(h)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(h, cache):
    with open(part2_cache(h), "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=1, default=float)


def run_cases(geoms, h, theva, theva77, cache, light=None):
    for geom in geoms:
        k = geom_key(geom)
        if k in cache:
            continue
        if not geometry_fits(geom):
            cache[k] = dict(key=k, geom=geom, h=h, fits=False, track_width=track_width_of(geom)); save_cache(h, cache)
            print(f"    {k:60} does not fit the {PLATE_W:.1f} m plate ({track_width_of(geom):.2f} m)", flush=True)
            continue
        ev = Evaluator(geom, h, theva, theva77)
        r = ev.evaluate(light_from=(light.get(k) if light else None))
        cache[k] = r; save_cache(h, cache)
        print(f"    {k:60} NI {r['NI_req']/1e3:5.0f} kA (mean {r['B_mean_req']:.2f} T), rated {r['F_rated']/1e3:5.1f} kN/m at {r['NI_rated']/1e3:4.0f} kA {'ok' if r['rating_ok'] else 'NO'}, sat {r['sat']*100:4.1f} %, "
              f"Ic {r['ic12']:5.0f} A ({r['face12']}), tape {r['tape12']:5.1f} + load {r['tape3']:5.1f} kg/m, W {r['W_req']/1e3:5.1f} kJ/m, fill {r['fill']*100:3.0f} %  ({r['solves']} solves, {r['seconds']:.0f} s)", flush=True)
    return cache


def track_width_of(geom):
    wb = shape_width(geom["shape"]); split = 0.0 if geom["gc"] is None else wb + geom["gc"]
    return geom["n"] * geom["W"] + (geom["n"] - 1) * split + wb


def geometry_fits(geom):
    return track_width_of(geom) <= PLATE_W + 1e-9


def passes(r):
    """Rating, fit and plate saturation, all at once."""
    return "tape12" in r and r["rating_ok"] and r["fits"] and r["sat"] <= SAT_LIMIT + 1e-9


def objective(r, lift_tape, f_casing):
    if not r.get("fits", True) or "tape12" not in r:
        return float("inf")
    tape, casing = casing_mass(r, lift_tape)
    return casing * f_casing


def present_geom(allow=0.0, shape="square 100", gc=0.0):
    return dict(n=2, W=0.5, gc=gc, shape=shape, allow=allow)


def part_2(stage="screen"):
    global F_REQ
    d = design_NI()
    theva = d["theva"]; theva77 = read_theva_file(THEVA_77)
    F_REQ = d["F_s"]
    f_casing, f_cable = mass_factors()
    hr(f"PART 2: SUPPORT TRACK SEARCH ({stage})")
    print(f"  Requirement {F_REQ/1e3:.1f} kN/m per track; rating = the pull where the mean field over the bundles reaches {bm.B_WIND_LIMIT:.1f} T (bearing_mass.B_WIND_LIMIT,")
    print(f"  the rule behind the book's 280 kN/m at 346 kA), required to be at least {F_RATE_MARGIN:.2f} x the requirement ({F_RATE_MARGIN*F_REQ/1e3:.1f} kN/m); load winding {F_LOAD/1e3:.0f} kN/m")
    print(f"  on the same loops; plate {PLATE_W:.1f} m x {PLATE_T*1e3:.0f} mm fixed,")
    print(f"  gap {GAP*1e3:.0f} mm plus the cryostat wall allowance; the track must fit the plate.  Mass factors from ring_altitude.design(): {f_casing:.2f} kg of cable")
    print(f"  per kg on the casing, {f_cable:.2f} kg of cable per kg on the cable (plates fixed, so only the casing side moves).  Spacing of the center")
    print(f"  bundles is edge to edge ('joined' = one bundle holding both loops' turns).  Tape currents: the worst tape of the track at the share;")
    print(f"  the field model carries {I_REF_12:.0f} A per 12 mm tape and {I_REF_3:.0f} A per 3 mm tape (the present currents) to place its tapes.")
    print(f"  Coil ends: assumed bent down below the ring with a loop as large as needed, so they do not set the tape limit [UNVERIFIED].")
    # the lift tracks at the worst-tape current (Part 1: 60 A limit, 30 A at 50 percent), carried unchanged into every candidate's mass
    lift_tape = d["NI_l"] / (0.5 * 59.5) * 2 * 0.003 * TAPE_T * TAPE_RHO
    print(f"  Lift tracks carried unchanged: {d['NI_l']/1e3:.0f} kA per leg at 30 A (half of the Part 1 worst-tape limit 59.5 A), {lift_tape:.1f} kg/m of tape per track.")
    stages = dict(screen=0.005, confirm25=0.0025, confirm125=0.00125)
    if stage in stages:
        h = stages[stage]
        cache = load_cache(h)
        if stage == "screen":
            geoms = [present_geom()]
            for n in (1, 2, 3):
                for W in (0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2):
                    for gc in ((None, 0.0, 0.05, 0.10, 0.15, 0.20, 0.30) if n >= 2 else (0.0,)):
                        geoms.append(dict(n=n, W=W, gc=gc, shape="square 100", allow=0.02))
            print(f"\n  2.1  Loops, width and spacing at 5 mm cells, square 100 mm bundles, wall 20 mm ({len(geoms)} cases, the present design first at wall 0):")
            cache = run_cases(geoms, h, theva, theva77, cache)
            ranked = sorted((r for r in cache.values() if passes(r)), key=lambda r: objective(r, lift_tape, f_casing))
            best4 = [r["geom"] for r in ranked[:4]]
            geoms = []
            for g in best4 + [present_geom(allow=0.02)]:
                for shape in SHAPES:
                    geoms.append(dict(g, shape=shape))
            print(f"\n  2.2  Bundle cross-sections on the four best geometries and on the present one ({len(geoms)} cases):")
            cache = run_cases(geoms, h, theva, theva77, cache)
            ranked = sorted((r for r in cache.values() if passes(r)), key=lambda r: objective(r, lift_tape, f_casing))
            best2 = [r["geom"] for r in ranked[:2]]
            geoms = []
            for g in best2 + [present_geom()]:
                for allow in (0.0, 0.01, 0.02, 0.03, 0.04):
                    geoms.append(dict(g, allow=allow))
            print(f"\n  2.3  Cryostat wall allowance 0 to 40 mm on the two best and on the present ({len(geoms)} cases):")
            cache = run_cases(geoms, h, theva, theva77, cache)
            # one-at-a-time sensitivities around the best
            best = ranked[0]["geom"]
            geoms = []
            for W in (best["W"] - 0.1, best["W"] + 0.1):
                geoms.append(dict(best, W=round(W, 3)))
            if best["gc"] is not None:
                for gc in (max(0.0, best["gc"] - 0.05), best["gc"] + 0.05):
                    geoms.append(dict(best, gc=round(gc, 3)))
            print(f"\n  2.4  Sensitivities around the best ({len(geoms)} extra cases):")
            cache = run_cases(geoms, h, theva, theva77, cache)
            # 2.5 the combined grid: the three cross-sections that raised the limit most, with widths, spacings and walls
            geoms = []
            for shape in ("flat 200x50", "wide 150x67", "wide 150x100"):
                for W in (0.3, 0.4, 0.5):
                    for gc in (0.0, 0.05, 0.10):
                        for allow in (0.0, 0.01, 0.02):
                            geoms.append(dict(n=2, W=W, gc=gc, shape=shape, allow=allow))
                for allow in (0.0, 0.02):
                    geoms.append(dict(n=3, W=0.3, gc=None, shape=shape, allow=allow))
                    geoms.append(dict(n=1, W=1.0, gc=0.0, shape=shape, allow=allow))
            geoms = [g for g in geoms if geom_key(g) not in cache]
            print(f"\n  2.5  The combined grid: flat and wide bundles with two loops of 0.3 to 0.5 m, center spacing 0 to 100 mm, wall 0 to 20 mm,")
            print(f"       plus three joined 0.3 m loops and one 1 m loop with the same bundles ({len(geoms)} new cases):")
            cache = run_cases(geoms, h, theva, theva77, cache)
        else:
            cache5 = load_cache(0.005) if stage == "confirm25" else load_cache(0.0025)
            ranked = sorted((r for r in cache5.values() if passes(r)), key=lambda r: objective(r, lift_tape, f_casing))
            n_top = 8 if stage == "confirm25" else 2
            geoms = [r["geom"] for r in ranked[:n_top]] + ([present_geom()] if stage == "confirm25" else [])
            print(f"\n  Confirming the {n_top} best of the coarser stage{' and the present design' if stage == 'confirm25' else ' (the present design is in Part 1)'} at {h*1e3:g} mm:")
            cache = run_cases(geoms, h, theva, theva77, cache, light=(cache5 if stage == "confirm125" else None))
        report_tables(cache, lift_tape, f_casing, f_cable, d, h)
    elif stage == "report":
        for h in (0.005, 0.0025, 0.00125):
            cache = load_cache(h)
            if cache:
                report_tables(cache, lift_tape, f_casing, f_cable, d, h, full=True)


def report_tables(cache, lift_tape, f_casing, f_cable, d, h, full=False):
    rows = [r for r in cache.values() if "tape12" in r]
    ok = [r for r in rows if passes(r)]
    present_key = geom_key(present_geom())
    sat_ref = cache[present_key]["sat"] if present_key in cache and "sat" in cache[present_key] else None
    ranked = sorted(ok, key=lambda r: objective(r, lift_tape, f_casing))
    print()
    note = "" if sat_ref is None else f" (the present design computes to {sat_ref*100:.1f} % here)"
    print(f"  Ranking at {h*1e3:g} mm cells ({len(rows)} cases run, {len(ok)} meet rating, fit and plate saturation <= {SAT_LIMIT:.0%}{note}); cable-equivalent = casing side x {f_casing:.2f}:")
    print(f"  {'rank':>4} {'geometry':60} {'NI kA':>6} {'rated':>6} {'sat %':>5} {'Ic A':>5} {'I50 A':>5} {'tape kg/m':>9} {'load kg/m':>9} {'casing kg/m':>11} {'cable-eq kg/m':>13} {'W kJ/m':>7} {'MVA':>5} {'fill %':>6}")
    for i, r in enumerate(ranked[:12]):
        tape, casing = casing_mass(r, lift_tape)
        print(f"  {i+1:4d} {r['key']:60} {r['NI_req']/1e3:6.0f} {r['F_rated']/1e3:6.1f} {r['sat']*100:5.1f} {r['ic12']:5.0f} {r['I12']:5.0f} {r['tape12']:9.1f} {r['tape3']:9.1f} {casing:11.1f} {casing*f_casing:13.0f} {r['W_req']/1e3:7.1f} {r['VA']/1e6:5.0f} {r['fill']*100:6.0f}")
    if present_key in cache and "tape12" in cache[present_key]:
        r = cache[present_key]; tape, casing = casing_mass(r, lift_tape)
        print(f"  {'now':>4} {r['key']:60} {r['NI_req']/1e3:6.0f} {r['F_rated']/1e3:6.1f} {r['sat']*100:5.1f} {r['ic12']:5.0f} {r['I12']:5.0f} {r['tape12']:9.1f} {r['tape3']:9.1f} {casing:11.1f} {casing*f_casing:13.0f} {r['W_req']/1e3:7.1f} {r['VA']/1e6:5.0f} {r['fill']*100:6.0f}")
        print(f"       (the present design at wall 0; its tape at the design currents 163 A and 43 A was 2 x (67.3 + 55.8) + 2 x 47.1 = 340 kg/m, casing side 681 kg/m)")
    rejected = [r for r in rows if not passes(r)]
    if rejected:
        print(f"  Rejected ({len(rejected)}): " + "; ".join(f"{r['key']} ({'rating' if not r['rating_ok'] else 'fit' if not r['fits'] else 'saturation'})" for r in rejected[:30]))
    if not full or not ranked:
        return
    # shapes: how the worst-tape limit rises with the cross-section (the best geometry's family and the present one)
    print()
    print("  Bundle cross-section against the worst-tape limit (same loops, same ampere-turns for the requirement):")
    for base in {geom_key(dict(r["geom"], shape="square 100")) for r in ranked[:4]} | {geom_key(present_geom(allow=0.02, shape="square 100"))}:
        fam = [r for r in rows if geom_key(dict(r["geom"], shape="square 100")) == base]
        if len(fam) < 2:
            continue
        ref = next((r for r in fam if r["geom"]["shape"] == "square 100"), fam[0])
        print(f"    {base}:")
        for r in sorted(fam, key=lambda q: -q["ic12"]):
            print(f"      {r['geom']['shape']:13} Ic {r['ic12']:5.0f} A ({(r['ic12']/ref['ic12']-1)*100:+5.1f} %), NI {r['NI_req']/1e3:5.0f} kA ({(r['NI_req']/ref['NI_req']-1)*100:+5.1f} %), tape {r['tape12']:5.1f} kg/m ({(r['tape12']/ref['tape12']-1)*100:+5.1f} %), rating {'ok' if r['rating_ok'] else 'NO'}, fill {r['fill']*100:3.0f} %")
    # allowance
    print()
    print("  Cryostat wall allowance (bundle face to plate = 100 mm + wall):")
    for base in {geom_key(dict(r["geom"], allow=0.0)) for r in rows}:
        fam = sorted([r for r in rows if geom_key(dict(r["geom"], allow=0.0)) == base], key=lambda q: q["geom"]["allow"])
        if len(fam) < 3:
            continue
        print(f"    {base}: " + "; ".join(f"{r['geom']['allow']*1e3:.0f} mm: NI {r['NI_req']/1e3:.0f} kA, Ic {r['ic12']:.0f} A, tape {r['tape12']:.1f} kg/m" for r in fam))
    # sensitivities: one line per variable from the cases around the best
    best = ranked[0]
    print()
    print(f"  Sensitivities around the best ({best['key']}); pull at the best's ampere-turns is not solved separately, so each line gives the")
    print(f"  ampere-turns for the requirement and the tape mass of the neighbouring case:")
    g = best["geom"]
    for var, name in (("W", "loop width"), ("gc", "center spacing"), ("allow", "wall allowance"), ("n", "loops per track"), ("shape", "bundle shape")):
        parts = []
        for r in sorted(rows, key=lambda q: str(q["geom"][var])):
            if r is best or any(r["geom"][k] != g[k] for k in g if k != var) or r["geom"][var] == g[var]:
                continue
            v = r["geom"][var]
            if var == "gc":
                vs = "joined" if v is None else f"{v*1e3:.0f} mm"
            elif var == "allow":
                vs = f"{v*1e3:.0f} mm"
            elif var == "W":
                vs = f"{v:g} m"
            else:
                vs = str(v)
            parts.append(f"{vs}: NI {(r['NI_req']/best['NI_req']-1)*100:+.1f} %, Ic {(r['ic12']/best['ic12']-1)*100:+.1f} %, tape {(r['tape12']/best['tape12']-1)*100:+.1f} %, rated {r['F_rated']/1e3:.0f} kN/m{'' if passes(r) else ' (fails)'}")
        print(f"    {name:16}: " + ("; ".join(parts) if parts else "no neighbouring case run"))
    # the tape current table (Table 7.2) for the best and the present
    print()
    NI_normal = None
    weight_share = d["m_casing"] * ring.G_NET / 2
    for label, r in (("best", best), ("present", cache.get(present_key))):
        if not r or "lim12" not in r:
            continue
        NI_normal = r["NI_req"] * math.sqrt(weight_share / r["F_lin_req"])      # the weight share is far below saturation: linear scaling
        ic70, ic77, slope, trows = temperature_table(r["lim12"], r["face12"], NI_normal_frac=NI_normal / r["NI_req"])
        print(f"  Tape current table, {label} ({r['key']}): worst 12 mm tape {ic70:.0f} A at 70 K, {ic77:.0f} A at 77.5 K ({slope:.1f} A per kelvin, {slope/ic70*100:.1f} % per kelvin);")
        print(f"    full load {r['NI_req']/1e3:.0f} kA per leg, normal operation (the weight share, {weight_share/1e3:.1f} kN/m) {NI_normal/1e3:.0f} kA")
        print(f"    {'share':>6} {'full-load current A':>20} {'warming the tape can take K':>28} {'normal current A':>16} {'turns per bundle':>16} {'tape kg/m':>10}")
        for s_, I, dT, I_n in trows:
            turns = r["turns12"] * r["I12"] / I
            print(f"    {s_*100:5.0f}% {I:20.0f} {dT:28.1f} {I_n:16.0f} {turns/ (len(r['lim12'])):16,.0f} {r['tape12']*r['I12']/I:10.1f}")
    # stray field of the three best and the present, at the normal (weight share) and the rating currents
    print()
    print("  Stray field (30 x 30 m linear box, 50 mm cells, as final_bearing Part D.2): normal = support at the weight share, lift off; rating = support at the")
    print("  requirement plus the load winding (ampere-turns from the log-log fit of the solved points, an extrapolation), lift at its present 446 kA.")
    weight_share = d["m_casing"] * ring.G_NET / 2
    ref = None
    for label, r in [("present", cache.get(present_key))] + [(f"rank {i+1}", q) for i, q in enumerate(ranked[:3])]:
        if not r or "lim12" not in r:
            continue
        NI_normal = r["NI_req"] * math.sqrt(weight_share / r["F_lin_req"])
        pts = sorted([(r["NI_req"], r["F_req"]), (r["NI_rated"], r["F_rated"])])
        NI_rating = math.exp(math.log(pts[1][0]) + (math.log(r["F_req"] + F_LOAD) - math.log(pts[1][1])) * (math.log(pts[1][0]) - math.log(pts[0][0])) / (math.log(pts[1][1]) - math.log(pts[0][1])))
        sf = stray_field_check(r["geom"], NI_normal, NI_rating, NI_lift=d["NI_l"])
        if ref is None:
            ref = sf
        print(f"    {label} ({r['key']}): normal {NI_normal/1e3:.0f} kA, rating {NI_rating/1e3:.0f} kA per leg")
        for mode in ("normal", "rating"):
            print(f"      {mode:7}: " + ", ".join(f"{k} {v*1e3:.3f} mT" + (f" ({v/ref[mode][k]:.2f} x)" if ref is not sf and ref[mode][k] > 0 else "") for k, v in sf[mode].items() if k in ("LIM stator, center", "largest over the stator", "1 m below the ring")))
        r["stray"] = sf
    print()
    if present_key in cache and "tape12" in cache[present_key]:
        ob, op = objective(best, lift_tape, f_casing), objective(cache[present_key], lift_tape, f_casing)
        print(f"  Is the present design near the optimum?  Best cable-equivalent {ob:,.0f} kg/m against the present {op:,.0f} kg/m ({(ob/op-1)*100:+.1f} %).")
    else:
        print(f"  Best cable-equivalent {objective(best, lift_tape, f_casing):,.0f} kg/m (the present design was not run at this cell size here; Part 1 has it).")

# =============================================================================
# Paul's decisions of 2026-10-21: the design cryostat allowance is 10 mm (5 mm as a sensitivity); the 1.0 T mean-field rule is
# dropped and the rating is the pull at which the worst tape reaches a share of its critical current at its own field and
# angle (the share was left as "XX": 60, 70 and 80 percent are reported, [QUESTION FOR AUTHOR]).  Then the wall check and Part 3.
# =============================================================================

WALL_DESIGN, WALL_SENS = 0.010, 0.005
RATE_SHARES = (0.6, 0.7, 0.8)
LIFT_REQ = fb.F_LIFT_TRACK                     # 196.5 kN/m per lift track
RHO_TI, RHO_SS = 4430.0, 7900.0                # kg/m^3 [UNVERIFIED: handbook values]
RES_TI_77, RES_SS_77 = 1.5e-6, 0.5e-6          # ohm m at 77 K, Ti-6Al-4V and 304 stainless [UNVERIFIED: handbook values, not from a run]
CP_TI_77, CP_SS_77 = 215.0, 200.0              # J/kg/K at 77 K [UNVERIFIED]
YIELD_TI_RT, YIELD_TI_77 = 880e6, 1300e6       # Pa, Ti-6Al-4V [UNVERIFIED: typical handbook values]
WALL_T, LN2_T, MLI_T = 0.002, 0.003, 0.005     # m: 2 mm Ti wall, 3 mm liquid nitrogen, the rest of the 10 mm as insulation face (Paul)


def tape_rating(points, NI_req, I_req, w_tape, shares=RATE_SHARES):
    """The pull at which the worst tape of a winding built for the requirement at 50 percent reaches each share of its
    critical current, from the solved points (share interpolated linearly in log NI, pull in log-log)."""
    key = "lim12" if w_tape == 0.012 else "lim3"
    pts = sorted(points, key=lambda p: p["NI"])
    N_leg = NI_req / I_req                                         # turns per leg-equivalent, fixed by the requirement
    xs = []; ss = []; fs = []
    for p in pts:
        ic = min(min(d["ic_v"], d["ic_h"]) for d in p[key])
        xs.append(math.log(p["NI"])); ss.append((p["NI"] / N_leg) / ic); fs.append(math.log(p["F"]))
    out = {}
    for s in shares:
        if s < min(ss) or s > max(ss):
            out[s] = (float("nan"), float("nan")); continue
        x = float(np.interp(s, ss, xs))
        out[s] = (math.exp(x), math.exp(float(np.interp(x, xs, fs))))
    return out


def _evaluate_rerate(self, F_req=None, s_share=0.5, kind="support", margin=F_RATE_MARGIN):
    """The full evaluation on the new rules: requirement, two points above it (1.3 and 1.6 times) for the tape rating,
    the load winding at its own field and at the combined field (support kind only), the old 1.0 T rating by interpolation."""
    F_req = F_req or F_REQ
    w_main = 0.012 if kind == "support" else 0.003
    I_ref = I_REF_12 if kind == "support" else I_REF_3
    self.F_lin1 = self.run(self.NI_est, linear=True)
    pts = []
    NI1 = self.NI_est * math.sqrt(F_req / self.F_lin1 / 0.6)
    pts.append(self.curve_point(NI1))
    NI2 = self.interp_NI(pts, F_req)
    if abs(NI2 / NI1 - 1) < 0.004:
        NI2 = NI1 * 1.03
    pts.append(self.curve_point(NI2))
    NI_req = self.interp_NI(pts, F_req)
    p_req = min(pts, key=lambda p: abs(p["NI"] - NI_req))
    if abs(p_req["F"] / F_req - 1) > 0.01:
        pts.append(self.curve_point(NI_req)); p_req = pts[-1]; NI_req = self.interp_NI(pts, F_req)
    p13 = self.curve_point(1.3 * NI_req); pts.append(p13)
    if self.h >= 0.004 or kind == "lift":                          # at 2.5 mm the combined point below serves as the high point
        p16 = self.curve_point(1.6 * NI_req); pts.append(p16)
    else:
        p16 = p13
    if kind == "support":
        NI_l = self.interp_NI(pts, F_LOAD)
        p_l = self.curve_point(NI_l); pts.append(p_l)
        NI_l = self.interp_NI(pts, F_LOAD)
        p_c = self.curve_point(NI_req + NI_l); pts.append(p_c)      # both windings at their full ampere-turns: the upper bound for the 3 mm tape
        NI_446 = self.interp_NI(pts, F_req + F_LOAD)                 # the book's full-load state: requirement plus the load share, 446 kN/m together
        p_446 = self.curve_point(NI_446); pts.append(p_446)
        NI_446 = self.interp_NI(pts, F_req + F_LOAD)
    key = "lim12" if kind == "support" else "lim3"
    def limits_at(p, NI, k, w_t):
        out = []
        for d in p[k]:
            sc = NI / p["NI"]; r = dict(d)
            for face in ("v", "h"):
                B = d["B_" + face] * sc; th = d["th_" + face]
                r["ic_" + face] = float(ic_angle(self.theva, np.array([B]), np.array([th]), w_t)[0])
                r["ic77_" + face] = float(ic_angle(self.theva77, np.array([B]), np.array([th]), w_t)[0])
                r["B_" + face] = B
            out.append(r)
        return out
    lim = limits_at(p_req, NI_req, key, w_main)
    ic = min(min(d["ic_v"], d["ic_h"]) for d in lim); face = "v" if min(d["ic_v"] for d in lim) <= min(d["ic_h"] for d in lim) else "h"
    I_op = s_share * ic
    turns = sum(abs(d["NI"]) for d in lim) * NI_req / p_req["NI"] / I_op
    tape = turns * w_main * TAPE_T * TAPE_RHO
    rate = tape_rating(pts, NI_req, I_op, w_main)
    # the old rule for reference: mean field over the bundles reaching 1.0 T
    bm_pts = sorted(((p["NI"], p["B_mean"], p["F"]) for p in pts), key=lambda q: q[0])
    bms = [q[1] for q in bm_pts]
    if min(bms) <= bm.B_WIND_LIMIT <= max(bms):
        NI_1T = math.exp(float(np.interp(bm.B_WIND_LIMIT, bms, np.log([q[0] for q in bm_pts]))))
        F_1T = math.exp(float(np.interp(math.log(NI_1T), np.log([q[0] for q in bm_pts]), np.log([q[2] for q in bm_pts]))))
    else:
        NI_1T = F_1T = float("nan")
    W_req = p_req["W_mag"] * (NI_req / p_req["NI"])**2
    res = dict(key=geom_key(self.geom), geom=self.geom, h=self.h, kind=kind, NI_req=NI_req, F_req=F_req, F_lin_req=self.F_lin1 * (NI_req / self.NI_est)**2,
               eff_req=p_req["eff"], sat=p_req["sat"], B_plate=p_req["B_plate"], B_mean_req=p_req["B_mean"], W_req=W_req,
               F_rated=F_1T, NI_rated=NI_1T, rate_tape={str(int(s * 100)): v for s, v in rate.items()},
               track_width=self.track_width, fits=self.track_width <= PLATE_W + 1e-9, solves=self.solves, seconds=time.time() - self.t0, split=self.split,
               sat13=p13["sat"], sat16=p16["sat"], rules="tape")
    if kind == "support":
        lim3 = limits_at(p_l, NI_l, "lim3", 0.003)
        lim3c = limits_at(p_c, NI_req + NI_l, "lim3", 0.003)
        lim3f = limits_at(p_446, NI_446, "lim3", 0.003)
        ic3f = min(min(d["ic_v"], d["ic_h"]) for d in lim3f); I3f = s_share * ic3f
        turns3f = sum(abs(d["NI"]) for d in lim3) * NI_l / p_l["NI"] / I3f
        ic3 = min(min(d["ic_v"], d["ic_h"]) for d in lim3); ic3c = min(min(d["ic_v"], d["ic_h"]) for d in lim3c)
        I3 = s_share * ic3; I3c = s_share * ic3c
        turns3 = sum(abs(d["NI"]) for d in lim3) * NI_l / p_l["NI"] / I3; turns3c = sum(abs(d["NI"]) for d in lim3) * NI_l / p_l["NI"] / I3c
        W_load = p_l["W_mag"] * (NI_l / p_l["NI"])**2
        E_coil = W_load * fb.L_LOAD / fb.N_FAM
        res.update(ic12=ic, face12=face, I12=I_op, turns12=turns, tape12=tape, lim12=lim, NI_load=NI_l, W_load=W_load, E_coil=E_coil, VA=2 * E_coil / fb.TAU_HW,
                   ic3=ic3, I3=I3, turns3=turns3, tape3=turns3 * 0.003 * TAPE_T * TAPE_RHO, lim3=lim3,
                   ic3_comb=ic3c, I3_comb=I3c, turns3_comb=turns3c, tape3_comb=turns3c * 0.003 * TAPE_T * TAPE_RHO, lim3_comb=lim3c,
                   F_comb=p_c["F"], sat_comb=p_c["sat"], W_comb=p_c["W_mag"],
                   NI_446=NI_446, F_446=p_446["F"], sat_446=p_446["sat"], ic3_446=ic3f, I3_446=I3f, turns3_446=turns3f, tape3_446=turns3f * 0.003 * TAPE_T * TAPE_RHO, lim3_446=lim3f,
                   fill=(turns * 0.012 + turns3c * 0.003) * TAPE_T / (shape_area(self.geom["shape"]) * len(self.bundles)),
                   rating_ok=True)
    else:
        res.update(ic3=ic, face3=face, I3=I_op, turns3=turns, tape3=tape, lim3=lim, W_load=W_req, E_coil=W_req * fb.L_LIFT / fb.N_FAM,
                   VA=2 * W_req * fb.L_LIFT / fb.N_FAM / fb.TAU_HW, tape12=0.0, turns12=0.0, ic12=float("nan"), I12=float("nan"), face12="-",
                   fill=turns * 0.003 * TAPE_T / (shape_area(self.geom["shape"]) * len(self.bundles)), rating_ok=True, B_wall=None)
    return res


Evaluator.evaluate_rerate = _evaluate_rerate


def rated_tape(r, share):
    v = r.get("rate_tape", {}).get(str(int(share * 100)))
    return v[1] if v else float("nan")


def passes_tape(r, share, margin=F_RATE_MARGIN):
    return "rate_tape" in r and r["fits"] and r["sat"] <= SAT_LIMIT + 1e-9 and rated_tape(r, share) >= margin * r["F_req"]


def part3_cache(h, kind):
    return os.path.join(ring.out_dir(), f"bearing_optimize_part3_{kind}_{h*1e3:g}mm.json")


def load_json(p):
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_json(p, d):
    with open(p, "w", encoding="utf-8") as f:
        json.dump(d, f, indent=1, default=float)


def run_cases_rerate(geoms, h, theva, theva77, cache, path, kind="support", F_req=None):
    for geom in geoms:
        k = geom_key(geom)
        if k in cache:
            continue
        if not geometry_fits(geom):
            cache[k] = dict(key=k, geom=geom, h=h, fits=False, track_width=track_width_of(geom)); save_json(path, cache)
            print(f"    {k:60} does not fit the {PLATE_W:.1f} m plate ({track_width_of(geom):.2f} m)", flush=True)
            continue
        ev = Evaluator(geom, h, theva, theva77)
        r = ev.evaluate_rerate(F_req=F_req, kind=kind)
        cache[k] = r; save_json(path, cache)
        rt = " ".join(f"{s}:{rated_tape(r, s)/1e3:5.1f}" for s in RATE_SHARES)
        if kind == "support":
            print(f"    {k:60} NI {r['NI_req']/1e3:5.0f} kA, rated(tape 60/70/80 %) {rt} kN/m, 1 T {r['F_rated']/1e3:5.1f}, sat {r['sat']*100:4.1f} %, Ic {r['ic12']:5.0f} A, "
                  f"tape {r['tape12']:5.1f} + load {r['tape3']:5.1f} ({r['tape3_446']:5.1f} at 446 kN/m together, {r['tape3_comb']:5.1f} with both at full ampere-turns; Ic3 {r['ic3']:.0f}/{r['ic3_446']:.0f}/{r['ic3_comb']:.0f} A) kg/m, W {r['W_req']/1e3:5.1f} kJ/m ({r['solves']} solves, {r['seconds']:.0f} s)", flush=True)
        else:
            print(f"    {k:60} NI {r['NI_req']/1e3:5.0f} kA, rated(tape 60/70/80 %) {rt} kN/m, 1 T {r['F_rated']/1e3:5.1f}, sat {r['sat']*100:4.1f} %, Ic3 {r['ic3']:4.0f} A ({r['face3']}), "
                  f"tape {r['tape3']:5.1f} kg/m, W {r['W_req']/1e3:5.1f} kJ/m, {r['E_coil']/1e6:.3f} MJ per 50 m coil ({r['solves']} solves, {r['seconds']:.0f} s)", flush=True)
    return cache


def casing_mass_new(r_sup, r_lift):
    """Casing side from a support result and a lift result, both on the worst-tape basis, the load winding at the combined field."""
    tape = 2 * (r_sup["tape12"] + r_sup["tape3_comb"]) + 2 * r_lift["tape3"]
    return tape, tape * (1 + CRYO)


def ranking_tables(cache, title, kind, share_list=RATE_SHARES, lift_tape=None):
    rows = [r for r in cache.values() if "rate_tape" in r]
    for s in share_list:
        ok = sorted((r for r in rows if passes_tape(r, s)), key=lambda r: (r["tape12"] + r.get("tape3_comb", r.get("tape3", 0.0))) if kind == "support" else r["tape3"])
        print(f"\n  {title}: rating = the pull at which the worst tape reaches {s:.0%} of its limit, >= {F_RATE_MARGIN:.2f} x the requirement; {len(ok)} of {len(rows)} pass (with fit and saturation <= {SAT_LIMIT:.0%}):")
        if kind == "support":
            print(f"  {'rank':>4} {'geometry':60} {'NI kA':>6} {'rated kN/m':>10} {'1 T rule':>8} {'sat %':>5} {'Ic A':>5} {'I50 A':>5} {'tape':>6} {'load':>6} {'load*':>6} {'casing kg/m':>11} {'cable-eq':>9} {'W kJ/m':>7} {'MVA':>5} {'fill %':>6}")
            for i, r in enumerate(ok[:10]):
                tape = 2 * (r["tape12"] + r["tape3_comb"]) + 2 * (lift_tape or 0.0); casing = tape * (1 + CRYO)
                print(f"  {i+1:4d} {r['key']:60} {r['NI_req']/1e3:6.0f} {rated_tape(r, s)/1e3:10.1f} {r['F_rated']/1e3:8.1f} {r['sat']*100:5.1f} {r['ic12']:5.0f} {r['I12']:5.0f} {r['tape12']:6.1f} {r['tape3']:6.1f} {r['tape3_comb']:6.1f} {casing:11.1f} {casing*9.65:9.0f} {r['W_req']/1e3:7.1f} {r['VA']/1e6:5.0f} {r['fill']*100:6.0f}")
        else:
            print(f"  {'rank':>4} {'geometry':60} {'NI kA':>6} {'rated kN/m':>10} {'1 T rule':>8} {'sat %':>5} {'Ic3 A':>5} {'I50 A':>5} {'tape kg/m':>9} {'W kJ/m':>7} {'MJ/coil':>7} {'MVA 2 ms':>8} {'fill %':>6}")
            for i, r in enumerate(ok[:10]):
                print(f"  {i+1:4d} {r['key']:60} {r['NI_req']/1e3:6.0f} {rated_tape(r, s)/1e3:10.1f} {r['F_rated']/1e3:8.1f} {r['sat']*100:5.1f} {r['ic3']:5.0f} {r['I3']:5.0f} {r['tape3']:9.1f} {r['W_req']/1e3:7.1f} {r['E_coil']/1e6:7.3f} {r['VA']/1e6:8.0f} {r['fill']*100:6.0f}")
        rej = [r for r in rows if not passes_tape(r, s)]
        if rej:
            print("  rejected: " + "; ".join(f"{r['key']} ({'rating' if (r['fits'] and r['sat'] <= SAT_LIMIT) else 'fit' if not r['fits'] else 'saturation'}, {rated_tape(r, s)/1e3:.0f})" for r in rej[:20]))


def best_of(cache, share, kind="support", wall=None):
    rows = [r for r in cache.values() if "rate_tape" in r and passes_tape(r, share) and (wall is None or abs(r["geom"]["allow"] - wall) < 1e-9)]
    if not rows:
        return None
    return min(rows, key=lambda r: (r["tape12"] + r["tape3_comb"]) if kind == "support" else r["tape3"])


# ------------------------------------------------------------------ the cryostat wall
def wall_check(geom, h, theva, theva77, NI_req, NI_load, label, kind="support"):
    """The 2 mm titanium wall between the bundle face (10 mm below the plate's clearance) and the gap: bending stress against
    the coolant pressure for the span of one bundle width (clamped and simply supported strips), and the eddy loss in the
    wall when a driven winding pulses (its field at the wall plane, rising in the 2 ms protection rise or the 20 ms ramp)."""
    wb = shape_width(geom["shape"])
    print(f"  {label}: bundle {wb*1e3:.0f} mm wide; wall {WALL_T*1e3:.0f} mm Ti-6Al-4V, {LN2_T*1e3:.0f} mm liquid nitrogen behind it, {MLI_T*1e3:.0f} mm insulation face to the bundle (the 10 mm allowance).")
    print(f"  Coolant pressure at the coil: no cryogenic hydraulics code exists in the repository and Chapter 14 states no pressure [QUESTION FOR AUTHOR];")
    print(f"  the stress is tabulated against pressure.  A strip of width b and thickness t under pressure p: clamped edges sigma = p b^2 / (2 t^2),")
    print(f"  simply supported 3 p b^2 / (4 t^2) [UNVERIFIED: plate-strip formulas, Roark]; yield Ti-6Al-4V {YIELD_TI_RT/1e6:.0f} MPa at room temperature,")
    print(f"  {YIELD_TI_77/1e6:.0f} MPa at 77 K [UNVERIFIED: handbook values].  The outside of the cryostat is in vacuum, so the whole absolute pressure loads the wall.")
    print(f"  {'pressure':>9} {'clamped MPa':>11} {'margin RT':>9} {'margin 77 K':>11} {'supported MPa':>13} {'margin RT':>9} {'span for yield/2 at RT, clamped':>32}")
    for p_bar in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
        p = p_bar * 1e5
        s_cl = p * wb**2 / (2 * WALL_T**2); s_ss = 3 * p * wb**2 / (4 * WALL_T**2)
        b_half = math.sqrt(0.5 * YIELD_TI_RT * 2 * WALL_T**2 / p)
        print(f"  {p_bar:7.1f} bar {s_cl/1e6:11.0f} {YIELD_TI_RT/s_cl:9.2f} {YIELD_TI_77/s_cl:11.2f} {s_ss/1e6:13.0f} {YIELD_TI_RT/s_ss:9.2f} {b_half*1e3:29.0f} mm")
    # the eddy loss: the driven winding's field along the wall plane
    ev = Evaluator(geom, h, theva, theva77)
    for tag, NI in (("load winding pulse (its own field)", NI_load), ("support plus load (the full-load state)", NI_req + NI_load)) if kind == "support" else (("lift winding pulse", NI_req),):
        ev.run(NI)
        y_wall = -(MLI_T + LN2_T) - WALL_T / 2                      # the wall's mid-plane, 2 mm below the gap's lower edge... measured from the bundle face at -allow
        y_wall = -geom["allow"] + MLI_T + LN2_T + WALL_T / 2       # bundle face at -allow; insulation, nitrogen, then the wall nearest the gap
        xs = np.arange(-ev.track_width / 2 - 0.3, ev.track_width / 2 + 0.3, ev.w.h)
        Bx, By = field_at_points(ev.w, xs, np.full_like(xs, y_wall))
        # a long thin sheet in a normal field B_y(x, t): E_z(x) = -d/dt int B_y dx + C, C from zero net current along z
        dBdt_unit = By                                             # per unit of dB/dt scaling: B rises from 0 to By in tau
        Ez = -np.cumsum(dBdt_unit) * ev.w.h
        Ez -= Ez.mean()
        integral = float(np.sum(Ez**2) * ev.w.h)                   # (V/m)^2 m at unit rate (1/s)
        Bmax = float(np.abs(By).max())
        print(f"  {tag}: {NI/1e3:.0f} kA per leg; field through the wall plane up to {Bmax:.3f} T (the pulse's change at the wall).")
        print(f"  {'wall':>22} {'resistivity':>12} {'W/m of ring in the 2 ms rise':>28} {'J/m per pulse (rise+fall)':>26} {'at the 20 ms ramp':>18} {'per launch, ring':>16} {'wall dT per pulse':>17} {'slit into 20 mm strips':>22}")
        for name, rho, dens, cp in (("Ti-6Al-4V", RES_TI_77, RHO_TI, CP_TI_77), ("304 stainless", RES_SS_77, RHO_SS, CP_SS_77), ("composite (non-conducting)", float("inf"), 1800.0, 800.0)):
            sigma = 0.0 if rho == float("inf") else 1.0 / rho
            P2 = sigma * WALL_T * integral / fb.TAU_HW**2           # W per m of ring during the 2 ms rise
            E2 = 2 * P2 * fb.TAU_HW                                   # rise and fall
            E20 = 2 * sigma * WALL_T * integral / (fb.RAMP + fb.TAU_HW)
            m_wall = dens * WALL_T * (xs[-1] - xs[0])
            dT = E2 / (m_wall * cp) if sigma else 0.0
            n_strip = (xs[-1] - xs[0]) / 0.020
            print(f"  {name:>22} {('%.2f' % (rho*1e6)) if sigma else '-':>9} uOhm m {P2:28.1f} {E2:26.2f} {E20:18.2f} {E2*ring.L_RING/1e6:13.1f} MJ {dT*1e3:14.2f} mK {E2/n_strip**2:16.3f} J/m")
        print(f"  (a continuous sheet across the whole track width plus 0.3 m each side; the eddy currents run along the ring and return across the sheet's width,")
        print(f"  so the loss scales with the cube of the sheet width and the thin-sheet assumption (no self-shielding) overstates it once the sheet's own time")
        print(f"  constant, about mu0 sigma t w, nears the rise time [UNVERIFIED]; slitting the wall along the ring into strips 20 mm wide across the track divides")
        print(f"  the loss by the square of the strip count, the last column; a non-conducting composite wall has none)")


# ------------------------------------------------------------------ coil length
def coil_length_table(r, kind, label):
    """Stored energy per coil, supply rating at the two rise times, supplies per track and the coil-end ripple heating against
    the coil length, from the result's stored energy per meter and its plate peak field (final_bearing Part C.3 and C.4)."""
    Wm = r["W_load"] if kind == "support" else r["W_req"]
    I_op = r["I3"]
    Bp = r["B_plate"]
    print(f"  {label}: stored {Wm/1e3:.1f} kJ per m of track at the winding's working pull; plate peak field {Bp:.2f} T; {I_op:.0f} A per tape.")
    print(f"  {'coil length':>11} {'MJ per coil':>11} {'MVA at 2 ms':>11} {'MVA at 20 ms':>12} {'kV per tape at 2 ms':>19} {'supplies per track':>18} {'ripple W/m (10 per bundle)':>26} {'ring MW':>8}")
    ell = 2 * (GAP + 0.05)
    for L in (12.5, 25.0, 50.0, 100.0, 200.0, 400.0):
        E = Wm * L / fb.N_FAM
        n = 10
        dB = Bp / n
        E_step = 0.5 * bd.eddy_cycle(dB, ell / bd.V_REL, 0.0, d_p=PLATE_T, mu_r=cf.PLATE_MU, rho_p=lc.FERRO_RESISTIVITY)
        steps = 2 * 2 * n / L
        P = 2 * E_step * (0.5 + 0.4) * bd.V_REL * steps
        print(f"  {L:9.1f} m {E/1e6:11.3f} {2*E/fb.TAU_HW/1e6:11.0f} {2*E/(fb.RAMP+fb.TAU_HW)/1e6:12.1f} {2*E/fb.TAU_HW/I_op/1e3:19.0f} {fb.N_FAM*ring.L_RING/L/1e6:15.1f} M {P:26.2f} {P*ring.L_RING/1e6:8.1f}")
    print(f"  (ripple as final_bearing Part C.3: each coil end steps the plate field by 1/n of the winding's field over {ell*1e3:.0f} mm, n = 10 coils per bundle;")
    print(f"  it applies while the winding is energized, that is during a pass for the driven windings; the planned 20 ms ramp and the 2 ms protection rise are kept)")


# ------------------------------------------------------------------ the stages
def part_2b(stage):
    """Redo of the Part 2 ranking on the new rules, and the wall check."""
    global F_REQ
    d = design_NI(); theva = d["theva"]; theva77 = read_theva_file(THEVA_77); F_REQ = d["F_s"]
    lift_tape = d["NI_l"] / (0.5 * 59.5) * 2 * 0.003 * TAPE_T * TAPE_RHO
    hr(f"PART 2 REDONE ON PAUL'S RULES ({stage}): wall {WALL_DESIGN*1e3:.0f} mm design, {WALL_SENS*1e3:.0f} mm sensitivity; rating from the worst tape's share")
    print(f"  Requirement {F_REQ/1e3:.1f} kN/m; the winding is built for it at 50 percent of the worst tape's limit; the rating is the pull at which that tape")
    print(f"  reaches 60, 70 or 80 percent of its limit at its own field and angle (the share was left open [QUESTION FOR AUTHOR]), required >= {F_RATE_MARGIN:.2f} x the")
    print(f"  requirement.  The load winding is sized at its own field and at the combined field (support at the requirement plus load at {F_LOAD/1e3:.0f} kN/m).")
    print(f"  Plate saturation <= {SAT_LIMIT:.0%} at the requirement; fit in the {PLATE_W:.1f} m plate; lift tracks carried at {lift_tape:.1f} kg/m of tape per track until Part 3.")
    stages = dict(rerate=0.005, rerate25=0.0025, rerate125=0.00125)
    h = stages.get(stage, 0.0025)
    path = part3_cache(h, "support"); cache = load_json(path)
    if stage == "rerate":
        geoms = [present_geom()]
        for shape in ("square 100", "flat 200x50", "wide 150x67", "wide 150x100", "square 141"):
            for W in (0.4, 0.5):
                for gc in (0.0, 0.05):
                    for allow in (0.0, WALL_SENS, WALL_DESIGN, 0.02):
                        geoms.append(dict(n=2, W=W, gc=gc, shape=shape, allow=allow))
        for shape in ("flat 200x50", "wide 150x67"):
            for allow in (WALL_SENS, WALL_DESIGN):
                geoms.append(dict(n=2, W=0.3, gc=0.0, shape=shape, allow=allow))
                geoms.append(dict(n=1, W=1.0, gc=0.0, shape=shape, allow=allow))
        print(f"\n  2b.1  {len(geoms)} geometries at 5 mm cells (two loops of 0.4 and 0.5 m, five bundle shapes, spacing 0 and 50 mm, walls 0, 5, 10, 20 mm; plus 0.3 m and one loop):")
        cache = run_cases_rerate(geoms, h, theva, theva77, cache, path)
    elif stage in ("rerate25", "rerate125"):
        coarse = load_json(part3_cache(0.005 if stage == "rerate25" else 0.0025, "support"))
        geoms = []
        for s in RATE_SHARES:
            for wall in (WALL_DESIGN, WALL_SENS):
                b = best_of(coarse, s, "support", wall)
                if b and b["geom"] not in geoms:
                    geoms.append(b["geom"])
        if stage == "rerate25":
            geoms.append(present_geom()); geoms.append(present_geom(allow=WALL_DESIGN))
            ranked = sorted((r for r in coarse.values() if "rate_tape" in r and passes_tape(r, 0.7) and abs(r["geom"]["allow"] - WALL_DESIGN) < 1e-9), key=lambda r: r["tape12"] + r["tape3_comb"])
            for r in ranked[:6]:
                if r["geom"] not in geoms:
                    geoms.append(r["geom"])
        else:
            geoms = geoms[:2]
        print(f"\n  2b.2  {len(geoms)} geometries at {h*1e3:g} mm:")
        cache = run_cases_rerate(geoms, h, theva, theva77, cache, path)
    elif stage == "wall":
        c25 = load_json(part3_cache(0.0025, "support"))
        hr("THE CRYOSTAT WALL")
        for label, r in (("present loops at the design wall", c25.get(geom_key(present_geom(allow=WALL_DESIGN)))), ("best at the design wall (70 percent rule)", best_of(c25, 0.7, "support", WALL_DESIGN))):
            if r:
                wall_check(r["geom"], 0.005, theva, theva77, r["NI_req"], r["NI_load"], f"{label}: {r['key']}")
                coil_length_table(r, "support", f"load winding on {r['key']}")
        return
    ranking_tables(cache, f"Support track ranking at {h*1e3:g} mm", "support", lift_tape=lift_tape)
    for wall in (WALL_DESIGN, WALL_SENS, 0.0):
        line = []
        for s in RATE_SHARES:
            b = best_of(cache, s, "support", wall)
            line.append(f"{s:.0%}: {b['key'] if b else 'none'}" + (f" ({b['tape12']:.1f} + {b['tape3_comb']:.1f} kg/m, rated {rated_tape(b, s)/1e3:.0f})" if b else ""))
        print(f"  best at wall {wall*1e3:.0f} mm by share: " + "; ".join(line))
    pk = geom_key(present_geom())
    if pk in cache and "rate_tape" in cache[pk]:
        r = cache[pk]
        print(f"  present design (wall 0): rated(tape) " + ", ".join(f"{s:.0%} {rated_tape(r, s)/1e3:.0f} kN/m" for s in RATE_SHARES) + f"; 1 T rule {r['F_rated']/1e3:.0f}; load winding {r['tape3']:.1f} kg/m at its own field, {r['tape3_comb']:.1f} at the combined field (Ic3 {r['ic3']:.0f} against {r['ic3_comb']:.0f} A)")


def part_3(stage):
    """Part 3: the lift compensation tracks and the load winding, with the coil length."""
    global F_REQ
    d = design_NI(); theva = d["theva"]; theva77 = read_theva_file(THEVA_77); F_REQ = d["F_s"]
    hr(f"PART 3: LIFT COMPENSATION TRACKS AND THE LOAD WINDING ({stage})")
    print(f"  Lift track requirement {LIFT_REQ/1e3:.1f} kN/m per track (the pair shares the 30 km/s pull), 3 mm tape, driven, built at 50 percent of the worst tape;")
    print(f"  rating by the tape rule as for the support tracks, >= {F_RATE_MARGIN:.2f} x; plate {PLATE_W:.1f} m x {PLATE_T*1e3:.0f} mm on top of the cable; wall {WALL_DESIGN*1e3:.0f} mm")
    print(f"  ({WALL_SENS*1e3:.0f} mm as a sensitivity).  The field model places its tapes at {I_REF_3:.0f} A.  Coil length from the stored energy (Part C.4) and the ripple (C.3).")
    stages = dict(lift=0.005, lift25=0.0025)
    h = stages.get(stage, 0.0025)
    path = part3_cache(h, "lift"); cache = load_json(path)
    if stage == "lift":
        geoms = [dict(n=1, W=1.0, gc=0.0, shape="square 100", allow=0.0)]          # the present lift track
        for shape in ("square 100", "flat 200x50", "wide 150x67", "wide 150x100"):
            for allow in (WALL_DESIGN, WALL_SENS):
                for W in (0.5, 0.6, 0.8, 1.0, 1.2):
                    geoms.append(dict(n=1, W=W, gc=0.0, shape=shape, allow=allow))
                for W in (0.3, 0.4, 0.5):
                    for gc in (0.0, 0.05):
                        geoms.append(dict(n=2, W=W, gc=gc, shape=shape, allow=allow))
        print(f"\n  3.1  {len(geoms)} lift track geometries at 5 mm (the present first):")
        cache = run_cases_rerate(geoms, h, theva, theva77, cache, path, kind="lift", F_req=LIFT_REQ)
    elif stage == "lift25":
        coarse = load_json(part3_cache(0.005, "lift"))
        geoms = [dict(n=1, W=1.0, gc=0.0, shape="square 100", allow=0.0)]
        for s in RATE_SHARES:
            for wall in (WALL_DESIGN, WALL_SENS):
                b = best_of(coarse, s, "lift", wall)
                if b and b["geom"] not in geoms:
                    geoms.append(b["geom"])
        print(f"\n  3.2  {len(geoms)} lift geometries at 2.5 mm:")
        cache = run_cases_rerate(geoms, h, theva, theva77, cache, path, kind="lift", F_req=LIFT_REQ)
    ranking_tables(cache, f"Lift track ranking at {h*1e3:g} mm", "lift")
    for wall in (WALL_DESIGN, WALL_SENS, 0.0):
        line = []
        for s in RATE_SHARES:
            b = best_of(cache, s, "lift", wall)
            line.append(f"{s:.0%}: {b['key'] if b else 'none'}" + (f" ({b['tape3']:.1f} kg/m, rated {rated_tape(b, s)/1e3:.0f})" if b else ""))
        print(f"  best lift track at wall {wall*1e3:.0f} mm by share: " + "; ".join(line))
    if stage == "lift25" or stage == "report3":
        pk = geom_key(dict(n=1, W=1.0, gc=0.0, shape="square 100", allow=0.0))
        for label, r in (("present lift track", cache.get(pk)), ("best lift track at the design wall (70 percent rule)", best_of(cache, 0.7, "lift", WALL_DESIGN))):
            if r and "rate_tape" in r:
                print()
                coil_length_table(r, "lift", f"{label}: {r['key']}")
        # the load winding's coil length on the chosen support track
        c25 = load_json(part3_cache(0.0025, "support"))
        b = best_of(c25, 0.7, "support", WALL_DESIGN)
        if b:
            print()
            coil_length_table(b, "support", f"load winding on the best support track at the design wall: {b['key']}")
            print(f"  Load winding on that track: {b['NI_load']/1e3:.0f} kA per leg for {F_LOAD/1e3:.0f} kN/m; worst 3 mm tape {b['ic3']:.0f} A at its own field, {b['ic3_comb']:.0f} A with the support winding at the")
            print(f"  requirement ({b['F_comb']/1e3:.0f} kN/m together, plate saturation {b['sat_comb']*100:.1f} %); at 50 percent {b['I3']:.0f} and {b['I3_comb']:.0f} A; tape {b['tape3']:.1f} and {b['tape3_comb']:.1f} kg/m per track.")
        # the casing side with both
        bl = best_of(cache, 0.7, "lift", WALL_DESIGN)
        if b and bl:
            tape, casing = casing_mass_new(b, bl)
            print(f"\n  Casing side with the best support track and the best lift track at the design wall, 70 percent rule, all windings at 50 percent of their worst tape")
            print(f"  (load winding at the combined field): tape {tape:.1f} kg/m, casing side {casing:.0f} kg/m, cable-equivalent {casing*9.65:,.0f} kg/m.  Book: 681 kg/m at the mean-field currents.")


# =============================================================================
# Part 3b (2026-10-23 and 24): Paul's decisions on sections 2b and 3, and the control plan
#   1. rating share 70 percent;  2. the load winding's design state is the full-load state (support at its requirement plus
#   the load share, 446 kN/m together), the both-at-full-ampere-turns case is the bound the control plan must rule out;
#   3. coolant 0.4 bar absolute at the coil;  4. a non-conducting composite cryostat wall;  5. the requirement from the
#   Design Case numbers (12,000 kg/m, 7.918 m/s2, 5,000 t on 100 m: 245.5 kN/m, not the 248.2 of the "on top" loop);
#   7. the tape's critical current at the worst tape's own field angle from the full Robinson angle data.
#   The rerun stage solves the chosen geometries again on that basis; the report stage prints everything from the caches.
# =============================================================================
F_REQ_DESIGN_CASE = ring.M_LOAD_M * ring.G_NET / 2 + fb.F_LOAD_TRACK     # 245.5 kN/m per support track
SHARE_DESIGN = 0.7
P_COOLANT = 0.4e5                                                        # Pa absolute at the coil (Paul: 70 K loops at about 0.4 atm, Chapter 6)
# G-10CR glass-cloth epoxy, the usual non-conducting cryostat laminate: tensile strength warp 415 MPa and fill 257 MPa at
# room temperature, 825 and 459 MPa at 76 K (Kasen et al. 1980 as quoted in a fusion design document found by web search;
# the original table was not read) [UNVERIFIED]; the fill direction is the weaker and is used for the design.
COMP_STRENGTH = dict(fill_RT=257e6, warp_RT=415e6, fill_76K=459e6, warp_76K=825e6)
BEST_SUPPORT = dict(n=2, W=0.5, gc=0.0, shape="wide 150x67", allow=WALL_DESIGN)
BEST_LIFT = dict(n=1, W=0.8, gc=0.0, shape="flat 200x50", allow=WALL_DESIGN)
PRESENT_LIFT = dict(n=1, W=1.0, gc=0.0, shape="square 100", allow=0.0)
# the sled passage, from the runs reported in Vol_III/sled_passage_spread_2026-10-10.md and final_bearing_2026-10-12.md
SPREAD_T, REG_BW_HZ, PROT_TRIGGER, PROT_DELAY = 0.050, 270.0, 0.010, 0.001
PASS_PULLS = ((1, 285, "load"), (3, 166, "load"), (5, 63, "load"), (10, 68, "lift"), (20, 251, "lift"), (30, 393, "lift"))   # km/s, kN/m on the pair


def part3b_cache(h, kind):
    return os.path.join(ring.out_dir(), f"bearing_optimize_part3b_{kind}_{h*1e3:g}mm.json")


def pull_curve(r):
    """(NI, F) of a cached record, sorted: the requirement, the tape-share points, and for a support track the load winding,
    the 446 kN/m state and the both-at-full point; below the first point the pull is taken as NI squared (the linear region)."""
    pts = [(r["NI_req"], r["F_req"])] + [tuple(v) for v in r["rate_tape"].values() if not math.isnan(v[0])]
    if "NI_load" in r:
        pts += [(r["NI_load"], F_LOAD), (r["NI_446"], r["F_446"]), (r["NI_req"] + r["NI_load"], r["F_comb"])]
    pts = sorted(set((float(a), float(b)) for a, b in pts))
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def F_of_NI(r, NI):
    x, y = pull_curve(r)
    if NI <= x[0]:
        return y[0] * (NI / x[0])**2
    return math.exp(float(np.interp(math.log(min(NI, x[-1])), np.log(x), np.log(y))))


def NI_of_F(r, F):
    x, y = pull_curve(r)
    if F <= y[0]:
        return x[0] * math.sqrt(F / y[0])
    return math.exp(float(np.interp(math.log(min(F, y[-1])), np.log(y), np.log(x))))


def worst_ic_at(r, theva, NI_total, key, ref_NI, w_tape):
    """The worst tape's critical current when the bundles carry NI_total ampere-turns per leg in all: the stored worst-tape
    fields (at ref_NI) scaled with the total ampere-turns, the angle unchanged."""
    sc = NI_total / ref_NI; best = None
    for d in r[key]:
        for face in ("v", "h"):
            ic = float(ic_angle(theva, np.array([d["B_" + face] * sc]), np.array([d["th_" + face]]), w_tape)[0])
            if best is None or ic < best[0]:
                best = (ic, d["B_" + face] * sc, d["th_" + face], face, d["x"])
    return best


def casing_side_3b(r_sup, r_lift):
    tape = 2 * (r_sup["tape12"] + r_sup["tape3_446"]) + 2 * r_lift["tape3"]
    return tape, tape * (1 + CRYO)


def requirement_note(d):
    g = ring.G_NET; w_half = ring.M_LOAD_M * g / 2; l_half = fb.F_LOAD_TRACK
    print(f"  5. The requirement.  bearing_mass.py, final_bearing.py Part A and Parts 1 to 3 of this optimization used {d['F_s']/1e3:.1f} kN/m per")
    print(f"     support track: design_NI iterates final_bearing's 'on top' loop, the casing {ring.M_LOAD_M:,.0f} kg/m plus the bearing's own casing side")
    print(f"     ({d['m_casing']-ring.M_LOAD_M:.0f} kg/m at the mean-field currents), so {d['m_casing']:,.0f} x {g:.3f} / 2 = {d['m_casing']*g/2/1e3:.1f} plus the load share {l_half/1e3:.1f} = {d['F_s']/1e3:.1f} kN/m.")
    print(f"     The Design Case numbers give {ring.M_LOAD_M:,.0f} x {g:.3f} / 2 = {w_half/1e3:.2f} kN/m (half the casing weight) plus 5,000,000 x {g:.3f} / 100 / 2 = {l_half/1e3:.2f} kN/m")
    print(f"     (half the parked 5,000 t on 100 m) = {F_REQ_DESIGN_CASE/1e3:.1f} kN/m.  Chapter 7 and final_bearing's own decision (Part A: the bearing's casing side is")
    print(f"     inside the 12,000 kg/m) put the bearing inside the casing budget, so the Design Case value is the requirement from here; the 'on top'")
    print(f"     value with the Part 3 casing side (999 kg/m) would be {(ring.M_LOAD_M+999)*g/2/1e3+l_half/1e3:.1f} kN/m.  1.5 x the casing weight share ({1.5*ring.M_LOAD_M*g/1e3:.1f} kN/m")
    print(f"     for the whole casing, {1.5*ring.M_LOAD_M*g/2/1e3:.1f} per track) is below it, so the parked load sets the requirement.")


def angle_check(theva4, rob):
    """Item 7: the tape's critical current at the worst tape's own field angle, from the full angle data."""
    print("  7. The tape's worst field angle.  The optimization's limit used four measured curves (0, 30, 90 and 120 degrees between")
    print("     the field and the tape's normal), linear in angle, with the 120 degree curve standing in for 60 degrees.  The full")
    print("     Robinson dataset at 70 K (Pantoja and Strickland 2024, doi 10.6084/m9.figshare.25893079.v2, 5 degree steps, 0 to 8 T):")
    for B in (1.0, 1.641, 2.0, 3.0):
        Bs = np.array(sorted(rob)); k = int(np.argmin(np.abs(Bs - B))); a, i = rob[Bs[k]]
        if abs(Bs[k] - B) < 1e-6:
            kmin = int(i.argmin()); kmax = int(i.argmax())
            print(f"     at {B:.3f} T: perpendicular (0 deg) {i[0]*1.2:.0f} A per 12 mm, along the face (90 deg) {float(np.interp(90, a, i))*1.2:.0f} A, "
                  f"lowest {i[kmin]*1.2:.0f} A at {a[kmin]:.0f} deg, highest {i[kmax]*1.2:.0f} A at {a[kmax]:.0f} deg; at 30 deg {float(np.interp(30, a, i))*1.2:.0f} A, "
                  f"at 60 deg {float(np.interp(60, a, i))*1.2:.0f} A, at 120 deg {float(np.interp(120, a, i))*1.2:.0f} A")
        else:
            th = np.arange(0, 181, 1.0)
            ic = ic_angle_full(rob, np.full(len(th), B), th, 0.012)
            ic4 = ic_angle(theva4, np.full(len(th), B), th, 0.012)
            print(f"     at {B:.3f} T (between the 1.5 and 2 T sheets): lowest {ic.min():.0f} A at {th[ic.argmin()]:.0f} deg (folded), perpendicular {ic[0]:.0f} A, at 30 deg {ic[30]:.0f} A "
                  f"(four-curve model: {ic4[30]:.0f} A), at 45 deg {ic[45]:.0f} A (four-curve {ic4[45]:.0f} A), at 60 deg {ic[60]:.0f} A (four-curve {ic4[60]:.0f} A)")
    print("     The 120 degree curve is this tape's maximum, not its 60 degree value, so the four-curve model overrated every tape whose")
    print("     field lies between 30 and 90 degrees from its normal.  The rerun below uses the full data, taking at each angle the")
    print("     lower of theta and 180 - theta (which way the tape is turned in the bundle is not fixed by the design).")


def wall_04(label, wb):
    """Items 3 and 4: the cryostat wall at 0.4 bar absolute against vacuum, titanium and the composite."""
    t = WALL_T
    s_cl = P_COOLANT * wb**2 / (2 * t**2); s_ss = 3 * P_COOLANT * wb**2 / (4 * t**2)
    print(f"  {label}: bundle {wb*1e3:.0f} mm wide, wall {t*1e3:.0f} mm spanning it, {P_COOLANT/1e5:.1f} bar absolute inside against vacuum outside.")
    print(f"     bending stress: clamped edges {s_cl/1e6:.0f} MPa, simply supported {s_ss/1e6:.0f} MPa [UNVERIFIED: plate-strip formulas as in section 2b.5].")
    print(f"     Ti-6Al-4V 2 mm (the comparison): margin on yield {YIELD_TI_RT/s_cl:.1f} clamped / {YIELD_TI_RT/s_ss:.1f} supported at room temperature, {YIELD_TI_77/s_cl:.1f} / {YIELD_TI_77/s_ss:.1f} at 77 K.")
    for name, S in COMP_STRENGTH.items():
        print(f"     G-10CR composite 2 mm, {name.replace('_', ' ')} strength {S/1e6:.0f} MPa [UNVERIFIED]: margin {S/s_cl:.1f} clamped / {S/s_ss:.1f} supported; "
              f"thickness for a margin of 3: {wb*math.sqrt(3*P_COOLANT/(2*S))*1e3:.1f} mm clamped, {wb*math.sqrt(9*P_COOLANT/(4*S))*1e3:.1f} mm supported")
    print(f"     Eddy heating of a composite wall during the driven-winding pulses: zero (no conducting path); the titanium rows of section 2b.5 stand as the comparison.")


def control_plan(r, rl, theva):
    """Part 3b: the operating points of the persistent support winding and the driven load winding on one support track,
    the limiter that rules out the both-at-full bound, and the rise times."""
    F_req = r["F_req"]; w_half = ring.M_LOAD_M * ring.G_NET / 2; w_full = 2 * w_half
    NI_s_cap = r["NI_req"]; NI_l_own = r["NI_load"]; NI_446 = r["NI_446"]; F_446 = r["F_446"]
    I12 = r["I12"]; I3 = r["I3_446"]
    turns12 = NI_s_cap / I12                       # turns per leg (r["turns12"] counts all four bundles, for the tape mass)
    turns3 = NI_l_own / I3                         # the load winding is wound to give its own 198 kN/m alone at the design-state current
    x, y = pull_curve(r)
    print(f"  Pull curve of the track (ampere-turns per leg -> pull): " + ", ".join(f"{a/1e3:.0f} kA {b/1e3:.0f} kN/m" for a, b in zip(x, y)) + "; NI squared below the first point.")
    print(f"  Windings: persistent 12 mm support winding {turns12:,.0f} turns per leg at {I12:.1f} A = {NI_s_cap/1e3:.0f} kA, its cap, the requirement {F_req/1e3:.1f} kN/m;")
    print(f"  driven 3 mm load winding {turns3:,.0f} turns per leg at {I3:.1f} A = {NI_l_own/1e3:.0f} kA, which gives its {F_LOAD/1e3:.0f} kN/m alone; {I3:.1f} A is 50 percent of the worst 3 mm tape's")
    print(f"  limit in the design state ({NI_446/1e3:.0f} kA in all, support at its cap, {F_446/1e3:.0f} kN/m together), in which the load winding supplies {(NI_446-NI_s_cap)/1e3:.0f} kA at {(NI_446-NI_s_cap)/turns3:.1f} A per tape.")
    print(f"  The limiter: the driven winding's supply current is capped so that NI_support (measured) + NI_load <= {NI_446/1e3:.0f} kA per leg, the design state, and")
    print(f"  the winding's own cap is {NI_l_own/1e3:.0f} kA ({I3:.1f} A per tape); the both-at-full bound ({(NI_s_cap+NI_l_own)/1e3:.0f} kA, {r['F_comb']/1e3:.0f} kN/m, plate {r['sat_comb']*100:.0f} percent")
    print(f"  saturated) can then not be commanded, by construction.")
    print()
    cases = [
        ("normal, both tracks: half the casing weight", w_half, 0.0),
        ("parked 5,000 t, both tracks, before the flux pumps move: support at its weight share, load winding the rest", w_half, F_req - w_half),
        ("parked 5,000 t, both tracks, after the flux pumps: support at its cap, load winding off", F_req, 0.0),
        ("launch at 1 km/s, both tracks: the load pair's 285 kN/m spread over 50 ms, half on this track", w_half, 285e3 / 2 - 0.0),
        ("one track lost at rest: the other at its weight share, load winding catches the other half", w_half, w_half),
        ("one track out for repair, after the flux pumps: support at the whole casing weight", w_full, 0.0),
        ("one track out and a standing load: support at its cap, load winding to the limiter (the design state)", F_req, F_446 - F_req),
        ("protection pulse at the full load rating with the support at its weight share", w_half, None),
        ("the bound: support at its cap and the load winding at its own full ampere-turns", F_req, "bound"),
    ]
    print(f"  {'case':<110} {'NI_s kA':>8} {'NI_l kA':>8} {'total':>6} {'pull kN/m':>9} {'I12 A':>6} {'share':>6} {'I3 A':>6} {'share':>6} {'limiter':>8}")
    for name, F_s, F_l in cases:
        NI_s = NI_of_F(r, F_s)
        if F_l == "bound":
            NI_l = NI_l_own
        elif F_l is None:
            NI_l = min(NI_l_own, NI_446 - NI_s)
        else:
            NI_l = max(NI_of_F(r, F_s + F_l) - NI_s, 0.0)
        tot = NI_s + NI_l; F = F_of_NI(r, tot)
        ic12 = worst_ic_at(r, theva, tot, "lim12", r["NI_req"], 0.012)[0]
        ic3 = worst_ic_at(r, theva, tot, "lim3_446", r["NI_446"], 0.003)[0]
        i12 = NI_s / turns12; i3 = NI_l / turns3
        ok = "ok" if tot <= NI_446 * 1.001 else "EXCEEDS"
        print(f"  {name:<110} {NI_s/1e3:8.0f} {NI_l/1e3:8.0f} {tot/1e3:6.0f} {F/1e3:9.0f} {i12:6.1f} {i12/ic12:6.0%} {i3:6.1f} {i3/ic3:6.0%} {ok:>8}")
    F_stand = F_446 - w_full
    print(f"  (share = the tape's current over the worst tape's critical current at that total field: 12 mm tape {I12:.0f} A against {r['ic12']:.0f} A at the requirement;")
    print(f"  3 mm tape {I3:.0f} A against {r['ic3_446']:.0f} A in the design state.  The design state carries {F_446/1e3:.0f} kN/m on one track, so with one track out the standing")
    print(f"  load on that stretch may not exceed {F_stand/1e3:.0f} kN/m, {F_stand*100/ring.G_NET/1e3:,.0f} t on 100 m; a parked 5,000 t ({fb.F_LOAD_TRACK*2/1e3:.0f} kN/m) with one track out needs {(w_full+2*fb.F_LOAD_TRACK)/1e3:.0f} kN/m")
    print(f"  and is outside the design: a track is not taken out under a parked spacecraft, and a spacecraft is not parked over a track that is out.)")
    print()
    print(f"  Rise times.  Driven windings: the planned {fb.RAMP*1e3:.0f} ms ramp and the {fb.TAU_HW*1e3:.0f} ms protection rise (final_bearing.py); the sled's load is spread over")
    print(f"  {SPREAD_T*1e3:.0f} ms centred on the pass with the current regulated against the measured gap at {REG_BW_HZ:.0f} Hz, protection at {PROT_TRIGGER*1e3:.0f} mm of drift")
    print(f"  with a {PROT_DELAY*1e3:.0f} ms sensor delay, per face and per side (sled_passage_spread and final_bearing reports).  Persistent winding: its flux pumps")
    print(f"  change the current slowly; the rate is not in the repository or the chapter [QUESTION FOR AUTHOR: the flux pump's ramp rate, A per second per tape],")
    print(f"  so the plan gives the driven winding every standing load until the pumps have taken it over, with no time limit (a driven superconducting")
    print(f"  winding holds a steady current with no loss but its joints and supply).")
    print()
    print(f"  During a sled passage (the pair, kN/m): " + "; ".join(f"{u} km/s {F} {k}" for u, F, k in PASS_PULLS) + " (Table 7.5 of the chapter, final_bearing run).")
    print(f"  Below 6,988 m/s the load windings on the support tracks pull with the sled's weight, half each; above it the lift tracks on the ceiling pull")
    print(f"  the casing down, half each; the support winding does not move during a pass.  The two driven windings of a face never act together;")
    print(f"  the opposite face's windings act only as protection, when the gap drifts {PROT_TRIGGER*1e3:.0f} mm, and then within their own limiter.  The lift track's")
    print(f"  worst 3 mm tape runs at {rl['I3']:.0f} A against {rl['ic3']:.0f} A at its requirement ({rl['F_req']/1e3:.1f} kN/m, {rl['NI_req']/1e3:.0f} kA); it shares no bundle with a persistent winding,")
    print(f"  so it has no combined state; its cap is its own rating.")


def summary_tables(c_sup_new, c_lift_new, c_sup_old, c_lift_old, theva):
    """The 70 percent summary tables on the new basis, the old basis beside."""
    s = SHARE_DESIGN
    print(f"  Support tracks, rated at the {s:.0%} share (the pull at which the worst 12 mm tape reaches {s:.0%} of its limit), the winding built at 50 percent:")
    print(f"  {'design':<60} {'basis':<34} {'req':>6} {'NI kA':>6} {'rated':>6} {'x req':>5} {'Ic12':>5} {'I50':>5} {'tape12':>7} {'Ic3 own/446/both':>17} {'load tape own/446/both':>24} {'sat':>5} {'W kJ/m':>7}")
    rows = []
    for label, cache, basis in (("Design Case 245.5 kN/m, full angle data", c_sup_new, "new"), ("248.2 kN/m, four-curve angles", c_sup_old, "old")):
        for g, name in ((BEST_SUPPORT, "best: 2 x 0.50 m touching, wide 150 x 67, wall 10 mm"), (present_geom(), "present loops and bundle, wall 0")):
            r = cache.get(geom_key(g))
            if not r or "rate_tape" not in r: continue
            rated = rated_tape(r, s)
            print(f"  {name:<60} {label:<34} {r['F_req']/1e3:6.1f} {r['NI_req']/1e3:6.0f} {rated/1e3:6.1f} {rated/r['F_req']:5.2f} {r['ic12']:5.0f} {r['I12']:5.0f} {r['tape12']:7.1f} "
                  f"{r['ic3']:5.0f}/{r['ic3_446']:3.0f}/{r['ic3_comb']:3.0f} {r['tape3']:7.1f}/{r['tape3_446']:6.1f}/{r['tape3_comb']:6.1f} {r['sat']*100:4.1f}% {r['W_req']/1e3:7.1f}")
    print()
    print(f"  Lift tracks, rated at the {s:.0%} share:")
    print(f"  {'design':<60} {'basis':<34} {'req':>6} {'NI kA':>6} {'rated':>6} {'x req':>5} {'Ic3':>5} {'I50':>5} {'tape3':>7} {'sat':>5} {'W kJ/m':>7} {'MJ/50 m coil':>12} {'MVA 2 ms':>8}")
    for label, cache in (("Design Case, full angle data", c_lift_new), ("Part 3, four-curve angles", c_lift_old)):
        for g, name in ((BEST_LIFT, "best: 1 x 0.80 m, flat 200 x 50, wall 10 mm"), (PRESENT_LIFT, "present: 1 x 1.0 m, square 100, wall 0")):
            r = cache.get(geom_key(g))
            if not r or "rate_tape" not in r: continue
            rated = rated_tape(r, s)
            print(f"  {name:<60} {label:<34} {r['F_req']/1e3:6.1f} {r['NI_req']/1e3:6.0f} {rated/1e3:6.1f} {rated/r['F_req']:5.2f} {r['ic3']:5.0f} {r['I3']:5.0f} {r['tape3']:7.1f} {r['sat']*100:4.1f}% {r['W_req']/1e3:7.1f} {r['E_coil']/1e6:12.3f} {r['VA']/1e6:8.0f}")
    print()
    for label, cs, cl in (("new basis", c_sup_new, c_lift_new), ("old basis", c_sup_old, c_lift_old)):
        rs = cs.get(geom_key(BEST_SUPPORT)); rl = cl.get(geom_key(BEST_LIFT))
        if rs and rl and "rate_tape" in rs and "rate_tape" in rl:
            tape, casing = casing_side_3b(rs, rl)
            tape_b = 2 * (rs["tape12"] + rs["tape3_comb"]) + 2 * rl["tape3"]
            print(f"  Casing side, {label}: best support track ({rs['tape12']:.1f} support + {rs['tape3_446']:.1f} load tape at the design state) and best lift track ({rl['tape3']:.1f}), two of each,")
            print(f"  formers, cryostat and supports at one tape mass: tape {tape:.1f} kg/m, casing side {casing:.0f} kg/m, cable-equivalent {casing*9.65:,.0f} kg/m; at the bound {tape_b:.1f} and {tape_b*(1+CRYO):.0f} kg/m.")


def part_3b(stage):
    global F_REQ
    d = design_NI(); theva4 = d["theva"]; theva77 = read_theva_file(THEVA_77)
    rob = read_robinson_angles(); theva = dict(theva4); theva["full"] = rob
    F_REQ = F_REQ_DESIGN_CASE
    hr(f"PART 3B: PAUL'S DECISIONS OF 2026-10-23 AND THE CONTROL PLAN ({stage})")
    print(f"  Decisions: rating share {SHARE_DESIGN:.0%}; the load winding's design state is the full-load state (support at its requirement plus the load share),")
    print(f"  the both-at-full case the bound; coolant {P_COOLANT/1e5:.1f} bar absolute; a non-conducting composite wall; the requirement from the Design Case numbers;")
    print(f"  the tape's limit from the full Robinson angle data.  Requirement used here {F_REQ/1e3:.1f} kN/m per support track, lift {LIFT_REQ/1e3:.1f} kN/m per track.")
    print()
    requirement_note(d)
    print()
    angle_check(theva4, rob)
    print()
    h = 0.0025
    if stage == "rerun":
        for kind, geoms, F in (("support", [BEST_SUPPORT, present_geom()], F_REQ), ("lift", [BEST_LIFT, PRESENT_LIFT], LIFT_REQ)):
            path = part3b_cache(h, kind); cache = load_json(path)
            print(f"  3b.{kind}: {len(geoms)} geometries at {h*1e3:g} mm on the new basis:")
            run_cases_rerate(geoms, h, theva, theva77, cache, path, kind=kind, F_req=F)
            print()
    c_sup_new = load_json(part3b_cache(h, "support")); c_lift_new = load_json(part3b_cache(h, "lift"))
    c_sup_old = load_json(part3_cache(h, "support")); c_lift_old = load_json(part3_cache(h, "lift"))
    hr("THE WORST TAPES ON THE NEW BASIS")
    for kind, cache, g, key, ref in (("support, best", c_sup_new, BEST_SUPPORT, "lim12", "NI_req"), ("support, present", c_sup_new, present_geom(), "lim12", "NI_req"),
                                     ("lift, best", c_lift_new, BEST_LIFT, "lim3", "NI_req"), ("lift, present", c_lift_new, PRESENT_LIFT, "lim3", "NI_req")):
        r = cache.get(geom_key(g))
        if not r or key not in r: continue
        w = 0.012 if key == "lim12" else 0.003
        ic, B, th, face, xpos = worst_ic_at(r, theva, r[ref], key, r[ref], w)
        ic4 = worst_ic_at(r, theva4, r[ref], key, r[ref], w)[0]
        print(f"  {kind}: worst {w*1e3:.0f} mm tape at x = {xpos:+.3f} m, {B:.3f} T at {th:.1f} degrees from the normal (faces {'vertical' if face == 'v' else 'horizontal'}): {ic:.1f} A (the four-curve model at the same tape {ic4:.1f} A)")
        if key == "lim12":
            ic3, B3, th3, f3, x3 = worst_ic_at(r, theva, r["NI_446"], "lim3_446", r["NI_446"], 0.003)
            print(f"     3 mm load tape in the design state: {B3:.3f} T at {th3:.1f} degrees, {ic3:.1f} A; at the bound {worst_ic_at(r, theva, r['NI_req']+r['NI_load'], 'lim3_comb', r['NI_req']+r['NI_load'], 0.003)[0]:.1f} A")
    print()
    hr("SUMMARY AT THE 70 PERCENT SHARE")
    summary_tables(c_sup_new, c_lift_new, c_sup_old, c_lift_old, theva)
    print()
    hr("THE CRYOSTAT WALL AT 0.4 BAR, COMPOSITE")
    wall_04("best support track", shape_width(BEST_SUPPORT["shape"])); print()
    wall_04("best lift track", shape_width(BEST_LIFT["shape"])); print()
    wall_04("present 100 mm bundle", 0.10)
    print()
    hr("THE CONTROL PLAN (PART 3B)")
    rs = c_sup_new.get(geom_key(BEST_SUPPORT)); rl = c_lift_new.get(geom_key(BEST_LIFT))
    if rs and rl and "rate_tape" in rs and "rate_tape" in rl:
        control_plan(rs, rl, theva)
    else:
        print("  (the new-basis caches are not complete: run --stage=rerun first)")


def main():
    part = "1"; quick = False; stage = "screen"
    for arg in sys.argv[1:]:
        if arg in ("--help", "-h"):
            print(__doc__); return
        elif arg.startswith("--part="):
            part = arg.split("=")[1]
        elif arg.startswith("--stage="):
            stage = arg.split("=")[1]
        elif arg == "--quick":
            quick = True
        else:
            print(f"Unknown option {arg}"); return
    t_all = time.time()
    hr("BEARING COIL AND GAP OPTIMIZATION (2026-10-13)")
    print(f"  Ring at {ring.ALTITUDE_KM:.0f} km: casing {ring.M_LOAD_M:,.0f} kg/m ({ring.M_LOAD_M*ring.G_NET/1e3:.1f} kN/m), cable hardware {ring.M_HARDWARE_M:,.0f} kg/m, cable {ring.M_CABLE_M:,.0f} kg/m, face {ring.CABLE_SIDE:.2f} m.")
    if "1" in part:
        part_1(quick)
    if part == "2":
        part_2(stage)
    if part == "2b":
        part_2b(stage)
    if part == "3":
        part_3(stage)
    if part == "3b":
        part_3b(stage)
    hr(f"DONE ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
