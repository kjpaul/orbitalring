#!/usr/bin/env python3
"""
Two-dimensional magnetostatics for a long track facing a plate: the cross-section
of a coil array, a possible iron yoke, the air gap and the ferromagnetic plate on
the cable.  Everything is uniform along the ring, so the vector potential has
one component A_z(x, y) and

    d/dx (1/mu dA/dx) + d/dy (1/mu dA/dy) = -J_z

on a uniform grid with A = 0 on the far boundary, solved by red-black successive
over-relaxation.  B_x = dA/dy, B_y = -dA/dx.  Iron is linear with a relative
permeability; saturation is reported by comparing the field in the iron with
the plate's 1.2 T, not modeled.

Used by compensation_tracks.py for Part B of the 2026-10-11 task.  Run this
file alone for the two checks: a wide current sheet (B = mu0 K / 2 on each
side) and a gapped iron circuit (B = mu0 NI / 2g).

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import math
import time

import numpy as np

try:
    import scipy.sparse as sps
    import scipy.sparse.linalg as spla
    HAVE_SCIPY = True
except ImportError:          # the relaxation solver still works without scipy
    HAVE_SCIPY = False

MU0 = 4e-7 * math.pi
SOLVER = "direct"            # "direct" (sparse LU through scipy, exact to round-off) or "sor" (red-black relaxation)


class Grid:
    def __init__(self, Lx, Ly, h, x0=None, y0=None):
        """Domain Lx by Ly (m), cell size h, origin at (x0, y0) (defaults: centered in x, y0 = -Ly/2)."""
        self.h = h
        self.nx = int(round(Lx / h)) + 1
        self.ny = int(round(Ly / h)) + 1
        self.x0 = -Lx / 2 if x0 is None else x0
        self.y0 = -Ly / 2 if y0 is None else y0
        self.x = self.x0 + np.arange(self.nx) * h
        self.y = self.y0 + np.arange(self.ny) * h
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing="ij")
        self.mu_r = np.ones((self.nx, self.ny))
        self.J = np.zeros((self.nx, self.ny))
        self.materials = []                                   # (mask, mu_r0, B_sat) for the saturating solve

    def box(self, x_lo, x_hi, y_lo, y_hi):
        return (self.X >= x_lo) & (self.X < x_hi) & (self.Y >= y_lo) & (self.Y < y_hi)

    def add_iron(self, x_lo, x_hi, y_lo, y_hi, mu_r, B_sat=None):
        m = self.box(x_lo, x_hi, y_lo, y_hi)
        self.mu_r[m] = mu_r
        if B_sat is not None:
            self.materials.append((m, mu_r, B_sat))

    def add_current(self, x_lo, x_hi, y_lo, y_hi, ampere_turns):
        """A rectangular conductor region carrying ampere_turns in total (sign = direction of z)."""
        m = self.box(x_lo, x_hi, y_lo, y_hi)
        n = m.sum()
        if n:
            self.J[m] += ampere_turns / (n * self.h * self.h)

    def solve(self, tol=1e-6, max_iter=50000, omega=None, verbose=False, A0=None, method=None):
        """Solve for A_z.  method 'direct' assembles the five-point system as a sparse matrix and solves it with a sparse LU
        (scipy); 'sor' is the red-black relaxation.  Both use the same face reluctivities, so they agree to the solver
        tolerance; the direct solve is exact to round-off and tens of times faster on fine grids."""
        method = method or (SOLVER if HAVE_SCIPY else "sor")
        if method == "direct":
            return self._solve_direct()
        nx, ny, h = self.nx, self.ny, self.h
        nu = 1.0 / self.mu_r                                  # reluctivity per cell
        # face reluctivities (harmonic mean of the two cells sharing the face)
        nxw = np.zeros((nx, ny)); nxe = np.zeros((nx, ny)); nys = np.zeros((nx, ny)); nyn = np.zeros((nx, ny))
        nxw[1:, :] = 2 * nu[1:, :] * nu[:-1, :] / (nu[1:, :] + nu[:-1, :])
        nxe[:-1, :] = nxw[1:, :]
        nys[:, 1:] = 2 * nu[:, 1:] * nu[:, :-1] / (nu[:, 1:] + nu[:, :-1])
        nyn[:, :-1] = nys[:, 1:]
        diag = nxw + nxe + nys + nyn
        rhs = MU0 * self.J * h * h
        A = np.zeros((nx, ny)) if A0 is None else A0.copy()
        if omega is None:
            omega = 2.0 / (1.0 + math.pi / max(nx, ny))
        i_idx = np.arange(nx)[:, None]; j_idx = np.arange(ny)[None, :]
        red = ((i_idx + j_idx) % 2 == 0)
        inner = np.zeros((nx, ny), dtype=bool); inner[1:-1, 1:-1] = True
        masks = [red & inner, (~red) & inner]
        t0 = time.time()
        for it in range(max_iter):
            delta_max = 0.0
            for m in masks:
                Anew = (nxw * np.roll(A, 1, axis=0) + nxe * np.roll(A, -1, axis=0)
                        + nys * np.roll(A, 1, axis=1) + nyn * np.roll(A, -1, axis=1) + rhs) / np.maximum(diag, 1e-30)
                d = omega * (Anew - A)
                A[m] += d[m]
                delta_max = max(delta_max, np.abs(d[m]).max() if m.any() else 0.0)
            if it % 200 == 0:
                scale = max(np.abs(A).max(), 1e-30)
                if verbose:
                    print(f"    iter {it}: max change {delta_max/scale:.2e}")
                if delta_max / scale < tol and it > 0:
                    break
        self.A = A
        self.Bx = np.zeros((nx, ny)); self.By = np.zeros((nx, ny))
        self.Bx[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
        self.By[1:-1, :] = -(A[2:, :] - A[:-2, :]) / (2 * h)
        self.iters = it
        self.solve_time = time.time() - t0
        return A

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
        # unknowns: the inner nodes; A = 0 on the boundary
        idx = -np.ones((nx, ny), dtype=np.int64)
        inner = np.zeros((nx, ny), dtype=bool); inner[1:-1, 1:-1] = True
        n_in = int(inner.sum())
        idx[inner] = np.arange(n_in)
        ii, jj = np.nonzero(inner)
        rows = [idx[ii, jj]]; cols = [idx[ii, jj]]; vals = [diag[ii, jj]]
        for di, dj, coef in ((-1, 0, nxw), (1, 0, nxe), (0, -1, nys), (0, 1, nyn)):
            i2, j2 = ii + di, jj + dj
            ok = inner[i2, j2]
            rows.append(idx[ii[ok], jj[ok]]); cols.append(idx[i2[ok], j2[ok]]); vals.append(-coef[ii[ok], jj[ok]])
        M = sps.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_in, n_in))
        b = rhs[ii, jj]
        a = spla.spsolve(M, b)
        A = np.zeros((nx, ny)); A[ii, jj] = a
        self.A = A
        self.Bx = np.zeros((nx, ny)); self.By = np.zeros((nx, ny))
        self.Bx[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
        self.By[1:-1, :] = -(A[2:, :] - A[:-2, :]) / (2 * h)
        self.iters = 1
        self.solve_time = time.time() - t0
        return A

    def solve_nonlinear(self, tol=1e-6, max_outer=25, relax=0.3, verbose=False, force_box=None):
        """Picard iteration with saturating iron: each material follows B(H) = mu0 H + B_sat tanh(mu0 (mu_r0 - 1) H / B_sat),
        so it starts at mu_r0 and saturates at about B_sat.  After each linear solve the cell permeability is reset from the
        cell's |B| and relaxed toward the new value.  The permeability distribution present at the call is the starting
        point, so a sequence of calls with rising current continues from the last solution.  Converged when the force
        through force_box (if given) changes by under 0.5 % and the mean relative permeability change is under 1 %.
        Returns the number of outer iterations; the fields are those of the last solve."""
        H_tab = np.concatenate([[0.0], np.logspace(0, 7.5, 400)])
        inv = {}
        for k, (m, mu0r, Bs) in enumerate(self.materials):
            B_tab = MU0 * H_tab + Bs * np.tanh(MU0 * (mu0r - 1) * H_tab / Bs)
            inv[k] = (B_tab, H_tab)
        A = getattr(self, "A", None)
        t0 = time.time()
        F_prev = None
        change = 1.0; dF = 1.0; n_still = 0
        for it in range(max_outer):
            A = self.solve(tol=10 * tol, A0=A)
            change = 0.0
            for k, (m, mu0r, Bs) in enumerate(self.materials):
                Bm = np.sqrt(self.Bx[m]**2 + self.By[m]**2)
                B_tab, H_tab = inv[k]
                Hm = np.interp(Bm, B_tab, H_tab)
                mu_new = np.where(Hm > 1e-9, Bm / (MU0 * np.maximum(Hm, 1e-9)), mu0r)
                mu_new = np.clip(mu_new, 1.0, mu0r)
                mu_old = self.mu_r[m]
                mu_rel = np.exp((1 - relax) * np.log(mu_old) + relax * np.log(mu_new))     # relax in log space
                change = max(change, float(np.mean(np.abs(mu_rel - mu_old) / mu_old)))
                self.mu_r[m] = mu_rel
            if force_box is not None:
                F = self.force_on_box(*force_box)
                dF = abs(F - F_prev) / max(abs(F), 1.0) if F_prev is not None else 1.0
                F_prev = F
            else:
                dF = 0.0
            if verbose:
                print(f"    outer {it}: mean permeability change {change:.4f}, force change {dF:.4f} ({time.time()-t0:.0f} s)")
            n_still = n_still + 1 if dF < 0.003 else 0
            if (change < 0.01 or n_still >= 2) and it > 2:
                break
        self.solve(tol=tol, A0=A)
        self.outer_iters = it + 1
        self.outer_change = change
        self.outer_dF = dF
        return it + 1

    def field_at(self, x, y):
        i = int(round((x - self.x0) / self.h)); j = int(round((y - self.y0) / self.h))
        i = min(max(i, 0), self.nx - 1); j = min(max(j, 0), self.ny - 1)
        return self.Bx[i, j], self.By[i, j]

    def line_By(self, y, x_lo, x_hi):
        """B_y along a horizontal line y from x_lo to x_hi: arrays x, By."""
        j = int(round((y - self.y0) / self.h))
        sel = (self.x >= x_lo) & (self.x <= x_hi)
        return self.x[sel], self.By[sel, j], self.Bx[sel, j]

    def force_on_line(self, y, x_lo, x_hi):
        """Maxwell stress through a horizontal line (per meter along z), N/m: positive = pull in +y on what is above."""
        xs, By, Bx = self.line_By(y, x_lo, x_hi)
        return np.sum((By**2 - Bx**2) / (2 * MU0)) * self.h

    def force_on_box(self, x_lo, x_hi, y_lo, y_hi):
        """Maxwell stress around a closed rectangle (per meter along z), N/m: the attraction of what is inside the box
        toward -y (positive = pulled down, toward a track below it)."""
        h = self.h
        i0 = int(round((x_lo - self.x0) / h)); i1 = int(round((x_hi - self.x0) / h))
        j0 = int(round((y_lo - self.y0) / h)); j1 = int(round((y_hi - self.y0) / h))
        Fy = 0.0
        for jj, sgn in ((j0, -1), (j1, +1)):
            Bx = self.Bx[i0:i1 + 1, jj]; By = self.By[i0:i1 + 1, jj]
            Fy += sgn * np.sum((By**2 - Bx**2) / (2 * MU0)) * h
        for ii, sgn in ((i0, -1), (i1, +1)):
            Bx = self.Bx[ii, j0:j1 + 1]; By = self.By[ii, j0:j1 + 1]
            Fy += sgn * np.sum(Bx * By / MU0) * h
        return -Fy

    def region_field(self, x_lo, x_hi, y_lo, y_hi):
        """Mean and max |Bx|, |By| over a region (for the winding field)."""
        m = self.box(x_lo, x_hi, y_lo, y_hi)
        return dict(Bx_mean=float(np.abs(self.Bx[m]).mean()), Bx_max=float(np.abs(self.Bx[m]).max()),
                    By_mean=float(np.abs(self.By[m]).mean()), By_max=float(np.abs(self.By[m]).max()),
                    B_max=float(np.sqrt(self.Bx[m]**2 + self.By[m]**2).max()))


def check_sheet(h=0.01):
    """A wide sheet of current K (A/m) in free space: B = mu0 K / 2 on each side, parallel to the sheet."""
    g = Grid(6.0, 3.0, h)
    K = 1e6
    g.add_current(-2.0, 2.0, -h / 2, h / 2, K * 4.0)
    g.solve(tol=1e-7)
    bx, by = g.field_at(0.0, 0.1)
    return bx, MU0 * K / 2, g.iters, g.solve_time


def check_gap(h=0.005):
    """Iron C-circuit with two 0.1 m gaps and NI ampere-turns: B_gap = mu0 NI / (2 g) for mu_r -> infinity."""
    g = Grid(2.0, 1.4, h)
    mu = 1000.0
    gap = 0.1
    # yoke: bottom bar, two legs, and a top bar (the "plate") with gaps between legs and plate
    g.add_iron(-0.6, 0.6, -0.5, -0.3, mu)                 # bottom bar
    g.add_iron(-0.6, -0.4, -0.3, 0.0, mu)                 # left leg (pole)
    g.add_iron(0.4, 0.6, -0.3, 0.0, mu)                   # right leg (pole)
    g.add_iron(-0.7, 0.7, gap, gap + 0.2, mu)             # plate above the gaps
    NI = 1.2 * 2 * gap / MU0
    g.add_current(-0.35, -0.05, -0.3, 0.0, NI)            # winding around the bottom bar: go conductor inside the window
    g.add_current(-0.35, -0.05, -0.8, -0.5, -NI)          # return conductor below the bar, outside the circuit
    g.solve(tol=1e-7)
    xs, By, Bx = g.line_By(gap / 2, -0.6, -0.4)
    return float(np.abs(By).mean()), MU0 * NI / (2 * gap), g.iters, g.solve_time


if __name__ == "__main__":
    bx, b_ref, it, t = check_sheet()
    print(f"Current sheet check: B parallel at 100 mm = {bx:.4f} T, mu0 K/2 = {b_ref:.4f} T ({it} iterations, {t:.1f} s)")
    bg, b_ref, it, t = check_gap()
    print(f"Gapped iron circuit check: B in the gap = {bg:.3f} T, mu0 NI/2g = {b_ref:.3f} T ({it} iterations, {t:.1f} s)")
