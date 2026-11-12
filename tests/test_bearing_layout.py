"""Checks of bearing_layout.py's circuit pieces that run in seconds: the vectorized pull-curve inverse, the demand read from a support
force history, the frozen point model, the soldered coil's L / R scaling, the hot-spot integration and the leakage inductances; Stage B: the
rate-limited support loop of both dynamics models and the switched bundle's circuit modes."""
import math

import numpy as np
import pytest

import bearing_layout as bl
import bearing_electrical as be
import bearing_optimize as bo


@pytest.fixture(scope="module")
def rec():
    r_sup, r_lift = be.load_records()
    return r_sup


@pytest.fixture(scope="module")
def consts(rec):
    return bl.part1_constants(rec)


def test_support_curve_inverse_matches_scalar(rec):
    cv = bl.SupportCurve(rec)
    for F in (10e3, 47.5e3, 120e3, 245e3):
        assert float(cv.NI_of_F(np.array([F]))[0]) == pytest.approx(bo.NI_of_F(rec, F), rel=1e-9)
    assert bo.F_of_NI(rec, float(cv.NI_of_F(np.array([47.5e3]))[0])) == pytest.approx(47.5e3, rel=1e-6)


def test_pump_demand_of_a_ramp(rec, consts):
    """A support force ramped from W to 1.1 W over 1 s and back: the swing is twice the ampere-turn change and the rate is its slope."""
    cv = bl.SupportCurve(rec)
    W = 2 * bl.W_HALF
    t = np.arange(0.0, 3.0, 1e-4)
    S = W * (1 + 0.1 * np.clip(t, 0, 1) - 0.1 * np.clip(t - 1.5, 0, 1))
    d = bl.pump_demand(t, S, cv, consts)
    dNI = bo.NI_of_F(rec, 1.1 * bl.W_HALF) - bo.NI_of_F(rec, bl.W_HALF)
    assert d["TV"] == pytest.approx(2 * dNI, rel=1e-3)
    assert abs(d["rate_pk"]) == pytest.approx(dNI / 1.0, rel=0.05)
    assert d["V_turn_m"] == pytest.approx(d["Lp"] * abs(d["rate_pk"]))


def test_frozen_passage_keeps_the_support_constant(rec):
    r = be.passage_point(5000.0, be.L_LOAD, 63e3, "load", r_sup=rec, coupled=False, loop_hz=0.0, t_settle=0.2)
    assert np.allclose(r["S"], r["S"][0])
    r1 = be.passage_point(5000.0, be.L_LOAD, 63e3, "load", r_sup=rec, coupled=False, loop_hz=1.0, t_settle=0.2)
    assert r1["S"].std() > 0


def test_solder_time_constant_scales_with_length_squared(consts):
    n_t = consts["n_t"]; Lp = bl.Lp_of(consts, consts["NI_s0"]); r = 1e-10
    taus = []
    for L_c in (1000.0, 10000.0):
        L10 = bl.N_SUP_COILS * 2 * n_t**2 * L_c * Lp
        R = (n_t - 1) * r / (4 * L_c * bl.W_TAPE)
        taus.append(L10 / R)
    assert taus[1] / taus[0] == pytest.approx(100.0)
    assert taus[0] > 1000 * 365.25 * 86400           # a soldered 1 km coil: thousands of years


def test_hot_spot_time_decreases_with_current():
    t1 = bl.hot_spot_time(51.0); t2 = bl.hot_spot_time(117.0)
    assert 0.1 < t2 < t1 < 100.0
    assert t1 / t2 > 3.0                               # heating goes as I^2


def test_copper_fit_values():
    assert 150.0 < bl.cp_cu(77.0) < 230.0              # NIST fit near 77 K
    assert 370.0 < bl.cp_cu(300.0) < 400.0
    assert 1.5e-9 < bl.rho_cu(77.0) < 3e-9


def test_leakage_inductance_is_small_and_positive(consts):
    tr = be.Track(bo.BEST_SUPPORT, 0.02)
    lk = bl.leakage_inductances(tr, consts, stride=9)
    for name, v in lk.items():
        assert 0 < v["L_leak"] < 0.05 * v["L_self"]
    assert lk["interleaved"]["L_leak"] < lk["block"]["L_leak"]
    # the ten coils in series are one loop: within a factor of two of the solver's loop inductance per turn squared (filament model, subsampled)
    Lloop = bl.Lp_of(consts, consts["NI_s0"]) * (bl.N_SUP_COILS * consts["n_t"])**2 / 2
    assert 0.5 < lk["block"]["L_sum"] / Lloop < 2.0


def test_energy_matches_quadratic_in_the_linear_range(consts):
    NI = 20e3
    assert bl.energy_per_m(consts, NI) == pytest.approx(consts["Lp_lin"] * NI**2, rel=0.02)


def test_short_coils_cancel_a_shared_loop_stretch(consts):
    """Coils much shorter than a stretch of load-winding ampere-turns lie inside it and keep their flux: the support loses the whole of it."""
    cs = be.CoilSystem(10.0, consts["L_lin"], n_fam=10, n_span=60, cells_per_coil=50)
    xc = cs.x[cs.n_cells // 2]
    dNI_s, _ = cs.solve(cs.stretch(xc, 125.0, 1000.0))
    inside = np.abs(cs.x - xc) < 50.0
    assert float(dNI_s[0][inside].mean()) == pytest.approx(-1000.0, rel=0.02)


def test_coil_end_heat_scales_with_the_coil_count():
    import final_bearing as fb
    sup, lift = fb.load_air()
    P1, Bp = bl.coil_end_heat(1000.0, 10, sup)
    P2, _ = bl.coil_end_heat(100.0, 10, sup)
    assert P2 / P1 == pytest.approx(10.0, rel=1e-6)
    Pe, _ = bl.coil_end_heat(1000.0, 10, sup, split=False)
    assert Pe > P1                                   # one full step heats more than ten steps of a tenth


# ---- Stage B (2026-11-12): the switched support winding

def test_rate_limited_support_loop_slews(rec):
    """A lane at rest with the casing 5 mm off: the 1 Hz loop's correction is slewed at steady_rate, so the support force changes no faster than
    that; without the limit it jumps."""
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    lanes = [pz.lane(no_sled=True, cancel=False, steady_loop_hz=1.0, w0=0.005, steady_rate=2e3, label="slewed"),
             pz.lane(no_sled=True, cancel=False, steady_loop_hz=1.0, w0=0.005, label="free")]
    res, hist = pz.point_pulsed(lanes, 30000.0, tabs, t_settle=0.5)
    dt = float(hist["t"][1] - hist["t"][0])
    slope = np.abs(np.diff(hist["S"], axis=1)) / dt
    assert slope[0].max() <= 2e3 * 1.001
    assert slope[1].max() > 1e4
    assert res[0]["max_exc_all"] >= res[1]["max_exc_all"]      # the slewed loop corrects the offset no better than the free one


def test_rate_limited_beam_loop(rec):
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    import compensation_tracks as ct
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    cases = [dict(label="false pull", no_sled=True, no_pulse=True, false_pull=(0.1, 20e-3, -50e3))]
    r_free, _ = pz.beam_pulsed(30000.0, 100.0, [sp.truss(10.0, 0.1)], cases, tabs, tau=20e-3, L_win=4000.0, t_settle=0.5, prot=ct.prot(0.010, 1e-3, cap_up=1e5, cap_down=1e5))
    r_slew, _ = pz.beam_pulsed(30000.0, 100.0, [sp.truss(10.0, 0.1)], cases, tabs, tau=20e-3, L_win=4000.0, t_settle=0.5, prot=ct.prot(0.010, 1e-3, cap_up=1e5, cap_down=1e5), loop_rate=1e3)
    assert r_slew[0]["pk_steady"] <= r_free[0]["pk_steady"] + 1.0   # the slewed support cannot rise above the free one's peak
    assert r_slew[0]["pk_steady"] < pz.W + 1e3 * 0.7                 # a 1 kN/m per s slew over the run's time


def test_switched_circuit_modes(consts):
    """The eigen-decomposition reproduces the common and differential time constants, and a coil opened at the supply voltage hands its
    ampere-turns to the nine within the leakage share."""
    sc = bl.switched_circuit(consts, 1000.0)
    lam = np.linalg.eigvalsh(sc["L"])
    assert lam.max() / bl.R_SW == pytest.approx(sc["tau_common"], rel=1e-9)
    assert lam.min() / bl.R_SW == pytest.approx(sc["tau_diff"], rel=1e-9)
    assert sc["tau_common"] > 20 * 86400 and sc["tau_diff"] < 2 * 3600
    N = sc["N"]; I0 = sc["I0"]
    t = np.linspace(0, 1.0, 1001)
    I = bl.circuit_run(sc["L"], np.array([-bl.V_SUP] + [0.0] * (N - 1)), np.full(N, I0), t)
    j = int(np.argmax(I[:, 0] <= 0.0))
    assert 0.2 < t[j] < 0.5                                         # L_leak I / V, about 0.34 s for the 1 km coil
    kept = (N - 1) * I[j, 1] / (N * I0)
    assert 0.995 < kept <= 1.0005
    # no power: the common mode decays with its own constant
    I_d = bl.circuit_run(sc["L"], np.zeros(N), np.full(N, I0), np.array([sc["tau_common"]]))
    assert I_d[0, 0] == pytest.approx(I0 / math.e, rel=1e-6)


def test_takeover_time_matches_the_integral(consts):
    sc = bl.switched_circuit(consts, 1000.0)
    t_lin = sc["L_common"] * (sc["I_cap"] - sc["I0"]) / bl.V_SUP
    t_int = consts["n_t"] * 2 * 1000.0 * bl.dphi(consts, consts["NI_s0"], consts["NI_cap"]) / bl.V_SUP
    assert t_int == pytest.approx(t_lin, rel=0.05)
    assert 100 < t_int < 200


def test_hot_spot_history_is_consistent_with_the_time():
    T = bl.hot_spot_T(lambda t: 51.0, bl.hot_spot_time(51.0, T_end=150.0), dt=2e-3)
    assert T == pytest.approx(150.0, rel=0.02)
