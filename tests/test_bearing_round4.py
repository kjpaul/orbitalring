"""Round 4 (2026-11-13): the fast coils' series turns and the supply invariant, the slow coils' allowance, the pull of a track against the
sled speed, the plate patch, and part_F's family constant.  No field solves: seconds."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import bearing_round4 as r4
import bearing_electrical as be
import bearing_design as bd
import bearing_optimize as bo
import final_bearing as fb


def test_series_turns_multiply_the_terminal_voltage():
    r_sup, r_lift = be.load_records()
    fc = r4.fast_coil_basics(r_lift)
    V1 = fc["Phi"] * r4.L_FAST / 2e-3
    assert 6.5e3 < V1 < 7.7e3                                      # one turn at 2 ms: about 7.1 kV
    for N in r4.SERIES_TURNS:
        I = fc["NI_coil"] / N
        assert abs(N * V1 * I - 2 * fc["E_coil"] / 2e-3) / (2 * fc["E_coil"] / 2e-3) < 1e-9   # the rating does not depend on N
    N_paul = round(fc["NI_coil"] / r4.I_SUPPLY_PAUL)
    assert N_paul == 14
    assert 95e3 < N_paul * V1 < 105e3                             # Paul's hand estimate: 100 kV at 2 ms


def test_supply_rating_per_metre_is_invariant():
    r_sup, r_lift = be.load_records()
    ref = 2 * r_lift["W_req"] / 2e-3
    for L in (12.5, 25.0, 50.0):
        for nf in (2, 5, 10, 20):
            E = r_lift["W_req"] * L / nf
            assert abs((2 * E / 2e-3) * (nf / L) - ref) / ref < 1e-12


def test_allowance_scales_with_the_inductance():
    V_d, t = 1.0, 10.0
    L_common, L_leak = 2135.0, 5.94
    assert abs((V_d * t / L_leak) / (V_d * t / L_common) - L_common / L_leak) < 1e-9
    assert V_d * 60.0 / L_common < 0.03                           # ten coils together: under 30 mA in a minute


def test_pull_of_a_track_against_speed():
    F1, k1 = r4.pull_track_of_u(1000.0); F30, k30 = r4.pull_track_of_u(30000.0); F_orb, _ = r4.pull_track_of_u(be.V_ORBIT)
    assert k1 == "load" and k30 == "lift"
    assert abs(F1 - 285e3 / 2) < 1.0 and abs(F30 - 393e3 / 2) < 1.0
    assert F_orb == 0.0


def test_plate_patch_restores_the_globals():
    W0, T0 = be.PLATE_W, be.PLATE_T
    with r4.plate_patch(1.2, 0.15):
        assert be.PLATE_W == 1.2 and bo.PLATE_T == 0.15
    assert be.PLATE_W == W0 and be.PLATE_T == T0 and bo.PLATE_W == W0


def test_part_F_family_constant_default():
    assert bd.N_FAM_F == 20 and fb.N_FAM == 20


def test_candidate_d_holds_the_hardware():
    hw, m_s, m_l = r4.hw_mass(r4.PLATE_CANDIDATES["d"])
    assert abs(hw - 10686.0) < 0.5
    hw_a = r4.hw_mass(r4.PLATE_CANDIDATES["a"])[0]
    assert abs(hw_a - 10371.6) < 1.0
