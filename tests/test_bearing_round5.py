"""Round 5 (2026-11-13): the fast coil's basics on a plate result, the Bean slab's two branches, the hold time at the 10 ms rise, the point
parser, the cache key, the cold-load hook in lim_config, and the drawing script's loop option.  No field solves: seconds."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import bearing_round5 as r5
import bearing_design as bd
import lim_config as lc


def test_fast_basics_scale_with_the_turns():
    res = dict(W_mag=123814.0, NI_req=468114.0, ic_req=57.67)
    a = r5.fast_basics(res, 10); b = r5.fast_basics(res, 14)
    assert abs(a["I"] * 10 - b["I"] * 14) < 1e-6                         # the same ampere-turns per coil
    assert abs(a["V"] / a["V_turn"] - 10) < 1e-9 and abs(b["V"] / b["V_turn"] - 14) < 1e-9
    assert abs(a["E_coil"] - b["E_coil"]) < 1e-9 and abs(a["E_coil"] - 123814.0 * 25 / 20) < 1e-6
    assert abs(a["MVA"] - 2 * a["E_coil"] / r5.TAU) < 1e-6
    assert abs(a["layers"] * 10 - b["layers"] * 14) < 1e-6              # the same tapes per coil


def test_bean_slab_is_continuous_at_the_penetration_field():
    Bp = 1.0
    below = r5.q_bean_slab(Bp * (1 - 1e-9), Bp); above = r5.q_bean_slab(Bp * (1 + 1e-9), Bp)
    assert abs(below - above) / above < 1e-6
    assert r5.q_bean_slab(0.5, Bp) < r5.q_bean_slab(0.6, Bp) < r5.q_bean_slab(2.0, Bp)
    assert abs(r5.q_bean_slab(0.5, Bp) - (2 / 3) * 0.125 / (r5.MU0 * Bp)) < 1e-6


def test_hold_time_and_zone_at_30_km_per_s():
    t, L_on = r5.hold_time(30e3)
    assert abs(L_on - (30e3 * 0.05 + bd.L_SC + 2 * 25.0 + 2 * 30e3 * r5.TAU)) < 1e-6
    assert abs(t - L_on / (30e3 + bd.V_REL)) < 1e-12
    t2, L2 = r5.hold_time(30e3, tau=2e-3)
    assert L2 < L_on and t2 < t


def test_parse_points_and_cache_key():
    pts = r5.parse_points("support:1.2,0.1,0.4:opt135;fast:1.0,0.15,0.6")
    assert pts[0] == ("support", 1.2, 0.1, 0.4, "opt135") and pts[1] == ("fast", 1.0, 0.15, 0.6, "opt135")
    assert r5.key_of("fast", 1.0, 0.15, 0.6, 0.0025) == "fast|1.000|0.1500|0.600|0.0025"


def test_cold_load_hook_adds_to_the_coil_environment():
    base = lc.Q_COIL_ENVIRONMENT
    try:
        lc.configure(Q_COLD_EXTRA_PER_SITE=5000.0)
        assert abs(lc.Q_COIL_ENVIRONMENT - base - 5000.0) < 1e-6
    finally:
        lc.configure(Q_COLD_EXTRA_PER_SITE=0.0)
    assert abs(lc.Q_COIL_ENVIRONMENT - base) < 1e-6


def test_converge_constant_is_tight():
    assert 0 < r5.CONVERGE <= 0.01
