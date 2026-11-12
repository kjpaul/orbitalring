"""Checks of bearing_stageB.py and the bias option of the dynamics models that run in seconds."""
import math

import numpy as np
import pytest

import bearing_stageB as sb
import bearing_layout as bl
import bearing_electrical as be


@pytest.fixture(scope="module")
def consts():
    r_sup, r_lift = be.load_records()
    return bl.part1_constants(r_sup)


def test_strip_resistance_per_metre():
    """A bare 3.5 um layer at 100 uohm cm: 95 ohm/m at 3 mm, 24 at 12 mm (the brief's 77 to 169 and 19 to 42 at 81 to 177 uohm cm)."""
    assert sb.strip_R_per_m(0.003) == pytest.approx(95.2, rel=0.01)
    assert sb.strip_R_per_m(0.012, 81e-8) == pytest.approx(19.3, rel=0.01)
    assert sb.strip_R_per_m(0.003, 177e-8) == pytest.approx(168.6, rel=0.01)


def test_settling_time_is_L_I_over_V(consts):
    """A 1 km coil, ten charged together: L_common x 51 A / 1 kV is about two minutes, as the brief estimates."""
    sc = bl.switched_circuit(consts, 1000.0)
    t = sc["L_common"] * sc["I0"] / 1000.0
    assert 100 < t < 125
    assert sc["L_common"] * sc["I0"] / 10.0 == pytest.approx(100 * t)


def test_bias_point_model_net_force_and_standing_currents():
    """With a bias both faces hold it at rest (the net force zero, the casing still), and a small net command moves both faces."""
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    import compensation_tracks as ct
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    P = ct.prot(1.0, 1e-3, cap_up=4e5, cap_down=4e5)                             # never triggers; sets the caps the no-sled lanes need
    lanes = [pz.lane(no_sled=True, cancel=False, gap_dep=False, steady_loop_hz=0.0, prot=P, bias=20e3, label="bias at rest"),
             pz.lane(no_sled=True, cancel=False, gap_dep=False, steady_loop_hz=0.0, prot=P, bias=20e3, false_pull=(0.0, 0.2, 5e3), label="bias, small pull"),
             pz.lane(no_sled=True, cancel=False, gap_dep=False, steady_loop_hz=0.0, prot=P, false_pull=(0.0, 0.2, 5e3), label="no bias, small pull")]
    res, hist = pz.point_pulsed(lanes, 30000.0, tabs, t_settle=0.3)
    assert res[0]["max_exc_all"] < 1e-6                                        # the bias is balanced: nothing moves
    assert res[0]["pk_hu"] == pytest.approx(20e3, rel=1e-3) and res[0]["pk_hd"] == pytest.approx(20e3, rel=1e-3)
    assert res[1]["pk_hu"] == pytest.approx(22.5e3, rel=1e-2)                   # bias + c/2
    assert res[2]["pk_hu"] == pytest.approx(5e3, rel=1e-2)                      # no bias: the whole command on one face
    # the same impulse, faster from a working current: the biased lane reaches its force sooner, so it moves more in the first 50 ms
    j = int(np.argmax(hist["t"] >= 0.05))
    assert abs(hist["w"][1, j] - hist["wc"][1, j]) > abs(hist["w"][2, j] - hist["wc"][2, j])


def test_bias_beam_model_runs():
    import sled_passage_pulsed as pz
    import sled_passage_simulation as sp
    import compensation_tracks as ct
    tabs = dict(sp.limiter_tables(), F0=0.0, Q0=0.0)
    cases = [dict(label="false pull", no_sled=True, no_pulse=True, false_pull=(0.1, 20e-3, -50e3))]
    P = ct.prot(0.010, 1e-3, cap_up=1e5, cap_down=1e5)
    r0, _ = pz.beam_pulsed(30000.0, 100.0, [sp.truss(10.0, 0.1)], cases, tabs, tau=20e-3, L_win=4000.0, t_settle=0.3, prot=P)
    r1, _ = pz.beam_pulsed(30000.0, 100.0, [sp.truss(10.0, 0.1)], cases, tabs, tau=20e-3, L_win=4000.0, t_settle=0.3, prot=P, bias=10e3)
    assert 0 < r1[0]["max_exc"] < 0.1 and abs(r1[0]["max_exc"] - r0[0]["max_exc"]) < 0.02


def test_two_track_coupling_falls_with_the_spacing():
    g1, legs_s, legs_l, xl1 = sb.two_track_linear(0.1, NI_l=1e3, h=0.02)
    g2, _, _, xl2 = sb.two_track_linear(1.0, NI_l=1e3, h=0.02)
    M1 = abs(sb.fig8_flux(g1, legs_s)[0]); M2 = abs(sb.fig8_flux(g2, legs_s)[0])
    assert M2 < M1
    gs, legs_s, _, _ = sb.two_track_linear(0.5, NI_s=1e3, h=0.02)
    L_s = sb.fig8_flux(gs, legs_s)[0] / 1e3
    assert 2e-6 < L_s < 4e-6                                                   # the figure-eight: twice Part 1d's 1.445 uWb/m per ampere-turn


def test_plate_heat_falls_with_thinner_plates():
    h1, e1, _ = sb.plate_heat_pass(0.8, 30000.0, 50.0, 0.20, 1.4)
    h2, e2, _ = sb.plate_heat_pass(0.8, 30000.0, 50.0, 0.10, 1.4)
    assert e2 < e1 and h2 < h1
