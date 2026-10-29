#!/usr/bin/env python3
"""
Checks the 2-D magnetostatic solver against known answers and against scikit-fem, with grid and box convergence.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (the sentence that the solver was checked).
Run: python bearing_fem_check.py
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

Part A of the 2026-10-12 bearing task: verify the field solver.

1. Known answers: two parallel wires (Biot-Savart, exact), an infinite current sheet (mu0 K / 2), and a C-core with two
   gaps (reluctance formula, which ignores fringing, so the finite-element and finite-difference answers both sit a few
   percent above it; the two numerical answers are compared with each other).
2. Independent discretization: the same material and current layout solved with scikit-fem (first-order triangular finite
   elements, scipy sparse direct solve) and compared with magnetostatic_2d.py (five-point finite differences, sparse direct
   solve) for the H7 compensation geometry (linear and saturating) and the S4 support geometry.  scikit-fem is an
   established open finite element package from PyPI; FEMM is not installed on this machine and pyfemm needs it.
3. Grid and domain convergence of the finite-difference solver: 10, 5 and 2.5 mm cells; domain 3.0 x 1.6 (the 2026-10-11
   runs), 3.8 x 2.0 and 6 x 3 m.

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import math
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import skfem
from skfem import MeshTri, ElementTriP1, Basis, BilinearForm, LinearForm, asm, condense, solve
from skfem.helpers import dot, grad

import magnetostatic_2d as ms
import compensation_fields as cf

MU0 = ms.MU0


def hr(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


# ---------------------------------------------------------------- FEM on the same cell layout
def fem_solve(g):
    """Solve d/dx(nu dA/dx) + d/dy(nu dA/dy) = -mu0 J on the grid's cell layout with P1 triangles (two per cell), A = 0 on the
    boundary.  Returns A at the grid nodes (nx, ny) so the same post-processing can be used."""
    nx, ny, h = g.nx, g.ny, g.h
    xs = g.x; ys = g.y
    m = MeshTri.init_tensor(xs, ys)
    # per-element material: element centroid -> cell (the cell value is defined at the node (i, j) covering [x_i, x_i + h))
    cx = m.p[0][m.t].mean(axis=0); cy = m.p[1][m.t].mean(axis=0)
    ic = np.clip(np.floor((cx - g.x0) / h).astype(int), 0, nx - 1)
    jc = np.clip(np.floor((cy - g.y0) / h).astype(int), 0, ny - 1)
    nu_e = 1.0 / g.mu_r[ic, jc]
    J_e = g.J[ic, jc]
    e = ElementTriP1()
    basis = Basis(m, e)

    @BilinearForm
    def stiff(u, v, w):
        return w["nu"] * dot(grad(u), grad(v))

    @LinearForm
    def load(v, w):
        return MU0 * w["J"] * v

    nu_field = basis.zero_w() + nu_e[None, :].repeat(basis.zero_w().shape[0], axis=0) if False else None
    # element-wise constant coefficients: pass as (n_elements, n_quad) arrays
    nq = basis.X.shape[1]
    K = asm(stiff, basis, nu=np.repeat(nu_e[:, None], nq, axis=1))
    b = asm(load, basis, J=np.repeat(J_e[:, None], nq, axis=1))
    D = basis.get_dofs()                     # boundary dofs
    A = solve(*condense(K, b, D=D))
    # the P1 dofs are the mesh vertices, ordered as the tensor mesh's points; map to (nx, ny)
    Agrid = np.zeros((nx, ny))
    ix = np.clip(np.round((m.p[0] - g.x0) / h).astype(int), 0, nx - 1)
    iy = np.clip(np.round((m.p[1] - g.y0) / h).astype(int), 0, ny - 1)
    Agrid[ix, iy] = A
    return Agrid


def set_fields_from_A(g, A):
    g.A = A
    nx, ny, h = g.nx, g.ny, g.h
    g.Bx = np.zeros((nx, ny)); g.By = np.zeros((nx, ny))
    g.Bx[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
    g.By[1:-1, :] = -(A[2:, :] - A[:-2, :]) / (2 * h)


def fem_nonlinear(g, force_box, max_outer=60, relax=0.3):
    """Picard iteration on the cell permeability with the FEM solve inside (same material law as the FD solver)."""
    H_tab = np.concatenate([[0.0], np.logspace(0, 7.5, 400)])
    inv = {}
    for k, (mk, mu0r, Bs) in enumerate(g.materials):
        inv[k] = (MU0 * H_tab + Bs * np.tanh(MU0 * (mu0r - 1) * H_tab / Bs), H_tab)
    F_prev = None; n_still = 0
    for it in range(max_outer):
        A = fem_solve(g); set_fields_from_A(g, A)
        for k, (mk, mu0r, Bs) in enumerate(g.materials):
            Bm = np.sqrt(g.Bx[mk]**2 + g.By[mk]**2)
            B_tab, H_tab2 = inv[k]
            Hm = np.interp(Bm, B_tab, H_tab2)
            mu_new = np.clip(np.where(Hm > 1e-9, Bm / (MU0 * np.maximum(Hm, 1e-9)), mu0r), 1.0, mu0r)
            g.mu_r[mk] = np.exp((1 - relax) * np.log(g.mu_r[mk]) + relax * np.log(mu_new))
        F = g.force_on_box(*force_box)
        dF = abs(F - F_prev) / max(abs(F), 1.0) if F_prev is not None else 1.0
        F_prev = F
        n_still = n_still + 1 if dF < 0.003 else 0
        if n_still >= 2 and it > 2:
            break
    A = fem_solve(g); set_fields_from_A(g, A)
    return it + 1, dF


# ---------------------------------------------------------------- known answers
def check_wires(h=0.01, L=6.0):
    """Two long parallel wires +I and -I, 0.4 m apart, field on the midline between them from Biot-Savart; the A = 0 box at
    +-L/2 adds an image error that falls with the box size."""
    I = 1e5
    g = ms.Grid(L, L, h)
    g.add_current(-0.2 - h, -0.2 + h, -h, h, I)
    g.add_current(0.2 - h, 0.2 + h, -h, h, -I)
    g.solve()
    pts = (0.0, 0.05, 0.1)
    out = []
    for y in pts:
        bx, by = g.field_at(0.0, y)
        # exact: each wire gives mu0 I / (2 pi r), direction perpendicular to r; on the midline the x components add
        r = math.hypot(0.2, y)
        B_exact = 2 * MU0 * I / (2 * math.pi * r) * (0.2 / r)
        out.append((y, abs(by), B_exact))
    return out


def main():
    t_all = time.time()
    hr("PART A.1: KNOWN ANSWERS")
    print("  Two parallel wires (+-100 kA, 0.4 m apart), field on the midline, finite-difference solver, 10 mm cells:")
    for L in (3.0, 6.0, 12.0):
        res = check_wires(0.01, L)
        print(f"    box {L:4.0f} m: " + ", ".join(f"y = {y*1e3:3.0f} mm: {b:.5f} T against {be:.5f} ({(b-be)/be*100:+.2f} %)" for y, b, be in res))
    print("  The remaining error at 12 m is the A = 0 box (an image current at the wall) and the cell size at the wire.")
    print("  Infinite current sheet (B = mu0 K / 2 on each side):")
    for h in (0.02, 0.01, 0.005):
        bx, bref, it, t = ms.check_sheet(h)
        print(f"    {h*1e3:4.1f} mm cells: {abs(bx):.4f} T against {bref:.4f} ({(abs(bx)-bref)/bref*100:+.2f} %); the sheet is 4 m wide in a 6 m box, the error is the finite width")
    print("  C-core with two 100 mm gaps, 200 mm poles (B = mu0 NI / 2 g for infinite permeability and no fringing):")
    for h in (0.01, 0.005, 0.0025):
        bg, bref, it, t = ms.check_gap(h)
        print(f"    {h*1e3:4.1f} mm cells: {bg:.4f} T against {bref:.4f} ({(bg-bref)/bref*100:+.2f} %), {t:.1f} s")
    print("  The C-core sits above the formula because the formula has no fringing: flux also crosses the gap outside the pole face, and")
    print("  the mean over the face rises.  The numbers converge with the grid (the 2.5 mm value is the reference for the FEM below).")

    hr("PART A.2: INDEPENDENT DISCRETIZATION (scikit-fem P1 triangles against the five-point finite differences)")
    print(f"  scikit-fem {skfem.__version__}: the same cells carry the same permeability and current; two triangles per cell; A = 0 on the boundary.")
    # C-core, both methods
    for h in (0.01, 0.005):
        g = ms.Grid(2.0, 1.4, h)
        mu = 1000.0; gap = 0.1
        g.add_iron(-0.6, 0.6, -0.5, -0.3, mu); g.add_iron(-0.6, -0.4, -0.3, 0.0, mu); g.add_iron(0.4, 0.6, -0.3, 0.0, mu); g.add_iron(-0.7, 0.7, gap, gap + 0.2, mu)
        NI = 1.2 * 2 * gap / MU0
        g.add_current(-0.35, -0.05, -0.3, 0.0, NI); g.add_current(-0.35, -0.05, -0.8, -0.5, -NI)
        t0 = time.time(); g.solve(); t_fd = time.time() - t0
        xs, By, Bx = g.line_By(gap / 2, -0.6, -0.4); b_fd = float(np.abs(By).mean())
        t0 = time.time(); A = fem_solve(g); t_fem = time.time() - t0; set_fields_from_A(g, A)
        xs, By, Bx = g.line_By(gap / 2, -0.6, -0.4); b_fem = float(np.abs(By).mean())
        print(f"  C-core {h*1e3:.0f} mm: FD {b_fd:.4f} T ({t_fd:.1f} s), FEM {b_fem:.4f} T ({t_fem:.1f} s), difference {(b_fem-b_fd)/b_fd*100:+.2f} %")
    # H7 linear and S4 linear at 10 and 5 mm
    geoms = {"H7 compensation (1 m, 300 mm pole, 200 mm slots, 100 mm plate)": dict(W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10),
             "S4 support (1 m, 200 mm poles x3 + halves, 50 mm slots, 100 mm plate)": dict(W=1.0, n_full=3, pole_w=0.20, slot_w=0.05, plate_t=0.10)}
    for name, geom in geoms.items():
        print(f"  {name}, linear, 100 kA-turns per pole:")
        for h in (0.01, 0.005):
            g, faces, slots, meta = cf.build_track(h=h, **geom)
            t0 = time.time(); g.solve(); t_fd = time.time() - t0
            r_fd = cf.evaluate_track(g, faces, slots, meta)
            t0 = time.time(); A = fem_solve(g); t_fem = time.time() - t0; set_fields_from_A(g, A)
            r_fem = cf.evaluate_track(g, faces, slots, meta)
            print(f"    {h*1e3:4.1f} mm: force FD {r_fd['F']/1e3:7.1f} kN/m, FEM {r_fem['F']/1e3:7.1f} ({(r_fem['F']-r_fd['F'])/r_fd['F']*100:+.2f} %); track mean FD {r_fd['B_track']:.4f}, FEM {r_fem['B_track']:.4f} T; "
                  f"plate flux FD {r_fd['B_plate']:.3f}, FEM {r_fem['B_plate']:.3f} T; stored FD {r_fd['W_mag']/1e3:.1f}, FEM {r_fem['W_mag']/1e3:.1f} kJ/m; winding across FD {r_fd['Bx_wm']:.3f}, FEM {r_fem['Bx_wm']:.3f} T  "
                  f"({t_fd:.1f} / {t_fem:.1f} s)")
    # H7 saturating at 100 and 200 kA, 10 mm, FD against FEM Picard
    print("  H7 saturating (plate 1.2 T, yoke 1.8 T), 10 mm cells, Picard iteration judged on the force:")
    for NI in (100e3, 200e3):
        g, faces, slots, meta = cf.build_track(h=0.01, NI=NI, W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10)
        box = (-meta["plate_w"] / 2 - 0.1, meta["plate_w"] / 2 + 0.1, meta["gap"] - 0.02, meta["gap"] + meta["plate_t"] + 0.02)
        t0 = time.time(); g.solve_nonlinear(force_box=box, max_outer=60); t_fd = time.time() - t0
        r_fd = cf.evaluate_track(g, faces, slots, meta)
        g2, faces, slots, meta = cf.build_track(h=0.01, NI=NI, W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10)
        t0 = time.time(); n_out, dF = fem_nonlinear(g2, box); t_fem = time.time() - t0
        r_fem = cf.evaluate_track(g2, faces, slots, meta)
        print(f"    {NI/1e3:.0f} kA: force FD {r_fd['F']/1e3:.1f} kN/m, FEM {r_fem['F']/1e3:.1f} ({(r_fem['F']-r_fd['F'])/r_fd['F']*100:+.2f} %); track mean {r_fd['B_track']:.3f} / {r_fem['B_track']:.3f} T; "
              f"plate max {r_fd['B_plate_max']:.2f} / {r_fem['B_plate_max']:.2f} T; saturated {r_fd['plate_sat_frac']*100:.0f} / {r_fem['plate_sat_frac']*100:.0f} %  ({t_fd:.0f} / {t_fem:.0f} s, FEM {n_out} outer, last dF {dF:.4f})")

    hr("PART A.3: GRID AND DOMAIN CONVERGENCE OF THE FINITE-DIFFERENCE SOLVER")
    print("  H7 linear, 100 kA per pole.  The 2026-10-11 runs used a 3.0 x 1.6 m box; the A = 0 wall was 0.3 m under the yoke and 0.6 m over the plate.")
    for (Lx, Ly) in ((3.0, 1.6), (3.8, 2.0), (6.0, 3.0), (10.0, 5.0)):
        for h in ((0.01, 0.005) if Ly <= 3.0 else (0.01,)):
            g, faces, slots, meta = cf.build_track(h=h, W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10, domain_x=Lx, domain_y=Ly)
            t0 = time.time(); g.solve(); t = time.time() - t0
            r = cf.evaluate_track(g, faces, slots, meta)
            print(f"    box {Lx:4.1f} x {Ly:3.1f} m, {h*1e3:4.1f} mm: force {r['F']/1e3:6.1f} kN/m, track mean {r['B_track']:.4f} T, plate flux {r['B_plate']:.3f} T, stored {r['W_mag']/1e3:.1f} kJ/m, "
                  f"winding across {r['Bx_wm']:.3f} T ({t:.1f} s)")
    print("  H7 linear at 2.5 mm in the 3.8 x 2.0 box:")
    g, faces, slots, meta = cf.build_track(h=0.0025, W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10)
    t0 = time.time(); g.solve(); t = time.time() - t0
    r = cf.evaluate_track(g, faces, slots, meta)
    print(f"    2.5 mm: force {r['F']/1e3:.1f} kN/m, track mean {r['B_track']:.4f} T, plate flux {r['B_plate']:.3f} T, stored {r['W_mag']/1e3:.1f} kJ/m, winding across {r['Bx_wm']:.3f} T ({t:.0f} s)")
    print("  S4 support linear, 100 kA per pole (50 mm slots, five cells at 10 mm):")
    for h in (0.01, 0.005, 0.0025):
        g, faces, slots, meta = cf.build_track(h=h, W=1.0, n_full=3, pole_w=0.20, slot_w=0.05, plate_t=0.10)
        t0 = time.time(); g.solve(); t = time.time() - t0
        r = cf.evaluate_track(g, faces, slots, meta)
        print(f"    {h*1e3:4.1f} mm: force {r['F']/1e3:6.1f} kN/m, track mean {r['B_track']:.4f} T, plate flux {r['B_plate']:.3f} T, winding across {r['Bx_wm']:.3f} T ({t:.0f} s)")
    print("  H7 saturating force curve, 10 and 5 mm cells (3.8 x 2.0 m box):")
    for h in (0.01, 0.005):
        out = cf.track_curve((50e3, 100e3, 150e3, 200e3, 300e3), W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, plate_t=0.10, h=h)
        print(f"    {h*1e3:4.1f} mm: " + ", ".join(f"{r['NI']/1e3:.0f} kA {r['F']/1e3:.0f}" for r in out[1:]) + " kN/m")
    print(f"\nDone ({time.time()-t_all:.0f} s)")


if __name__ == "__main__":
    main()
