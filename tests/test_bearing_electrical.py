"""Checks of bearing_electrical.py's models that run in seconds: the flux linkage read from the vector potential against the
stored energy, the flux-conserving coil system, the current division of parallel tapes, the lead heat leak and the invariance
of the supply rating in the number of series turns.  The solves use a coarse grid; the numbers are checked for consistency,
not against the book."""
import math

import numpy as np
import pytest

import bearing_electrical as be
import bearing_optimize as bo


def test_flux_equals_energy_single_loop():
    """One loop (the lift track) on a coarse grid: the flux per meter from A between the legs is 2 W / NI."""
    tr = be.Track(bo.BEST_LIFT, 0.02, w_tape=0.003, I_ref=be.I_REF_3, NI_ref=421e3)
    g, s = tr.solve(tr.NI_design(100e3), linear=True)
    assert s["phi"][0] == pytest.approx(2 * s["W"] / 100e3, rel=1e-6)
    assert s["phi"][0] > 0 and s["F"] > 0


def test_flux_equals_energy_two_loops():
    """The support track: the stored energy is NI x Phi' with both loops equal by symmetry."""
    tr = be.Track(bo.BEST_SUPPORT, 0.02)
    g, s = tr.solve(tr.NI_design(100e3), linear=True)
    assert s["phi"][0] == pytest.approx(s["phi"][1], rel=1e-3)
    assert s["W"] == pytest.approx(100e3 * s["phi"][0], rel=1e-3)


def test_coil_system_conserves_flux_and_dilutes():
    """A stretch of l inside coils of length L_c: every coil keeps its flux (the sum of the support change over the ring is
    -l x dNI).  One coil per bundle spreads the loss as -(l / L_c) dNI along the coil; ten staggered coils confine it to the
    stagger cells the stretch touches (here a 100 m stretch across the boundary of two 100 m cells: -dNI / 2 in each)."""
    L = np.array([[1.0e-6, 0.4e-6], [0.4e-6, 1.0e-6]])
    cs1 = be.CoilSystem(1000.0, L, n_fam=1, n_span=12, cells_per_coil=100)
    xc = cs1.x[cs1.n_cells // 2] + 500.0                      # mid-coil, so the stretch lies inside one coil
    dNI_s, delta = cs1.solve(cs1.stretch(xc, 100.0, 1000.0))
    assert float(dNI_s[0].sum() * cs1.dx) == pytest.approx(-100.0 * 1000.0, rel=0.02)
    assert float(np.interp(xc, cs1.x, dNI_s[0])) == pytest.approx(-0.1 * 1000.0, rel=0.02)
    cs = be.CoilSystem(1000.0, L, n_fam=10, n_span=12, cells_per_coil=100)
    xc = cs.x[cs.n_cells // 2]
    d = cs.stretch(xc, 100.0, 1000.0)
    dNI_s, delta = cs.solve(d)
    assert float(dNI_s[0].sum() * cs.dx) == pytest.approx(-100.0 * 1000.0, rel=0.02)
    assert float(np.interp(xc, cs.x, dNI_s[0])) == pytest.approx(-0.5 * 1000.0, rel=0.05)
    assert abs(float(np.interp(xc + 250.0, cs.x, dNI_s[0]))) < 0.01 * 1000.0       # nothing beyond the touched cells
    # the flux of every coil is kept: for each coil the integral of (dNI_s + dNI_l) over its span is zero (loop 1, symmetric)
    tot = cs.C @ ((dNI_s[0] + d[0]) * cs.dx)
    assert np.abs(tot).max() < 1e-6 * 100.0 * 1000.0
    # a stretch of two coil lengths takes the full dNI out of the coils inside it
    dNI_s2, _ = cs.solve(cs.stretch(xc, 2000.0, 1000.0))
    assert float(np.interp(xc, cs.x, dNI_s2[0])) == pytest.approx(-1000.0, rel=0.05)


def test_current_sharing_symmetric_pair_is_equal():
    """Two tapes placed symmetrically about the loop's mid-plane in y carry the same current; a tape nearer the return carries more."""
    xs = np.array([-0.4, -0.4]); ys = np.array([-0.02, -0.04])
    L = be.cable_inductance(xs, ys, 0.003, 0.0, 0.1, mu_plate=1.0)
    f = be.current_sharing(L)
    assert f[0] == pytest.approx(f[1], rel=0.05)              # near-equal: same distance to the return, the plate absent
    xs = np.array([-0.41, -0.39]); ys = np.array([-0.03, -0.03])
    f = be.current_sharing(be.cable_inductance(xs, ys, 0.003, 0.0, 0.1, mu_plate=1.0))
    assert f[1] > f[0]                                          # the inner tape (nearer the return leg) links less flux and takes more


def test_lead_heat_wiedemann_franz():
    """The 1-D lead at its optimum reproduces sqrt(L0 (T_h^2 - T_c^2)) per ampere to a few percent."""
    q, aol = be.lead_heat(1000.0)
    assert q == pytest.approx(1000.0 * math.sqrt(be.L0 * (be.T_WARM**2 - be.T_COLD**2)), rel=0.03)


def test_series_turns_keep_the_va():
    """V x I of the coil is the same for every number of series turns."""
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        rows = be.series_table("test", 700, 30.0, 14.3e3, 1.3e3, 0.3e6)
    va = [r["V2"] * r["I"] for r in rows]
    assert max(va) == pytest.approx(min(va), rel=1e-9)
    assert rows[0]["V_tt"] == 0.0 and rows[1]["V_tt"] == pytest.approx(14.3e3)


def test_phi_at_interpolation():
    rows = [dict(NI=100e3, phi=[0.14]), dict(NI=200e3, phi=[0.27])]
    assert be.phi_at(rows, 1.4e-6, 50e3) == pytest.approx(0.07)
    assert be.phi_at(rows, 1.4e-6, 150e3) == pytest.approx(math.exp(np.interp(math.log(150e3), np.log([100e3, 200e3]), np.log([0.14, 0.27]))))
