"""Reference check for Orbital Ring Engineering, Volume I: anyone with the repository and the packages of requirements.txt can
confirm that the code still prints the numbers in the book.

    pip install pytest
    python -m pytest tests            (about two minutes; the slow simulations are checked from your own run's logs)

Three kinds of test:
1. ring_altitude at the default design point (700 km, retrograde cable, cable hardware 10,686 kg/m), asserted to the
   precision the book quotes, then the prograde 700 km and retrograde 250 km cases through ring_altitude.design().
2. Every simulation that runs in well under a minute is run here, as run_all.py runs it, and its headline lines
   (tests/reference/<log>.txt, copied from the reference run of 2026-10-12) must appear in its output unchanged.
3. The slow simulations (minutes to half an hour each) are not run here.  Their headline lines are checked against the
   log of your own full run, <ORBITAL_RING_OUT>/<log>.txt (default runs/), made with  python run_all.py.  Those tests
   are skipped, with a message, when the log is not there.

The reference lines carry no run-time stamps or paths, so a difference means a number changed.
"""
import pytest

import ring_altitude as ring
from conftest import assert_lines_present, run_script, saved_log


# ---------------------------------------------------------------- 1. the design point

def test_default_design_point():
    """700 km, retrograde, hardware 10,686 kg/m: the values of Chapter 3 and Table 7.1."""
    assert ring.ALTITUDE_KM == 700.0 and ring.RETROGRADE and ring.M_HARDWARE_M == 10686.0
    assert ring.M_CABLE_STRUCTURAL == pytest.approx(115_005, abs=0.5)        # kg/m
    assert ring.M_CABLE_M == pytest.approx(125_691, abs=0.5)                 # kg/m, with hardware
    assert ring.M_RING_M == pytest.approx(137_691, abs=0.5)                  # kg/m, cable and casing
    assert ring.V_CABLE == pytest.approx(8_270.0, abs=0.05)                  # m/s
    assert ring.V_REL == pytest.approx(8_786.2, abs=0.05)                    # m/s, cable to casing
    assert ring.CABLE_AREA == pytest.approx(67.65, abs=0.005)                # m^2
    assert ring.CABLE_TENSION / 1e9 == pytest.approx(845.6, abs=0.05)        # GN
    assert ring.V_CABLE**2 / ring.R_ORBIT - ring.G_LOCAL == pytest.approx(1.7065, abs=0.00005)   # m/s^2, net outward
    assert ring.M_RING_M * ring.L_RING / 1e9 == pytest.approx(6_124, abs=0.5)                   # Mt
    assert ring.E_DEPLOY / 1e18 == pytest.approx(18.804, abs=0.0005)                            # EJ


RISING_M_STRUCT, RISING_V_CABLE, RISING_TENSION_GN, RISING_E_DEPLOY_EJ = 92883.313, 8330.820, 682.966, 19.1753   # the round-4 report, 2026-11-06


def test_rising_quadratic_design_point():
    """The design point since 2026-11-06: the quadratic with the build orbit (ORBITAL_RING_QUADRATIC=rising, the default
    outside these tests), 700 km, retrograde, hardware 10,686 kg/m; the values of the round-4 report."""
    d = ring.design(quadratic="rising")
    assert d["v_mom"] == pytest.approx(7_412.2, abs=0.05)                    # m/s, sqrt(GM r0) / R
    assert d["altitude_build_km"] == pytest.approx(527.36, abs=0.005)        # km, R / (1 + 2.5 %)
    assert d["m_struct"] == pytest.approx(RISING_M_STRUCT, abs=0.5)          # kg/m
    assert d["v_cable"] == pytest.approx(RISING_V_CABLE, abs=0.05)           # m/s
    assert d["tension"] / 1e9 == pytest.approx(RISING_TENSION_GN, abs=0.05)  # GN
    assert d["E_deploy"] / 1e18 == pytest.approx(RISING_E_DEPLOY_EJ, abs=0.0005)
    out = run_script(["ring_altitude.py"], env={"ORBITAL_RING_QUADRATIC": "rising"})
    assert f"Cable structural                       {RISING_M_STRUCT:>14,.3f} kg/m" in out


def test_prograde_700km():
    d = ring.design(altitude_km=700.0, retrograde=False)
    assert d["m_struct"] == pytest.approx(89_266, abs=0.5)
    assert d["v_cable"] == pytest.approx(8_343.3, abs=0.05)


def test_retrograde_250km():
    d = ring.design(altitude_km=250.0, retrograde=True)
    assert d["m_struct"] == pytest.approx(120_863, abs=0.5)
    assert d["v_cable"] == pytest.approx(8_506.3, abs=0.05)


def test_design_matches_module_constants():
    d = ring.design()
    assert d["m_struct"] == pytest.approx(ring.M_CABLE_STRUCTURAL, rel=1e-12)
    assert d["v_rel"] == pytest.approx(ring.V_REL, rel=1e-12)
    assert d["tension"] == pytest.approx(ring.CABLE_TENSION, rel=1e-12)


# ---------------------------------------------------------------- 2. the fast simulations, run here

FAST = [
    # log name, command after "python", extra environment
    ("ring_altitude_800", ["ring_altitude.py"], {"ORBITAL_RING_ALT_KM": "800"}),
    ("anchor_700_i0", ["anchor_line_analysis.py", "--inclination=0"], {}),
    ("anchor_700_i20", ["anchor_line_analysis.py", "--inclination=20"], {}),
    ("j2_retrograde_i20", ["j2_simulation.py", "--inclination=20"], {}),
    ("j2_prograde_i20", ["j2_simulation.py", "--inclination=20", "--cable=prograde"], {}),
    ("power_deployment", ["power_simulation.py"], {}),
    ("power_ops", ["power_simulation.py", "--demand=ops"], {}),
    ("power_deployment_ring_mounted_16MW", ["power_simulation.py", "--power=16"], {"ORBITAL_RING_ARRAY_OFFSET_KM": "0"}),
    ("levitation", ["levitation_simulation.py", "all"], {}),
    ("tether_stress", ["tether_stress.py"], {}),
    ("rendezvous", ["rendezvous_analysis.py"], {}),
    ("freefall", ["freefall_trajectory.py"], {}),
    ("atmospheric_loading", ["atmospheric_loading.py"], {}),
    ("climber_descent", ["climber_descent.py"], {}),
    ("lsm_chapter_metrics", ["tools/lsm_chapter_metrics.py"], {}),
    ("lsm_n_optimization", ["tools/lsm_n_optimization_demo.py"], {}),
    ("lsm_simulation", ["lsm_simulation.py"], {}),
    # round 13: the 130 m stator of the earlier Table 9.6 behind its switch must print what the script printed before
    ("lsm_simulation_table96", ["lsm_simulation.py"], {"ORBITAL_RING_LSM_STATOR": "table96"}),
    ("md_simulation", ["md_simulation.py", "--no-graphs"], {}),
    ("md_trade_study", ["md_trade_study.py"], {}),
    ("debris_simulation", ["debris_simulation.py"], {}),
    ("debris_orbits", ["debris_orbits.py"], {}),
    ("ch4_numbers", ["tools/ch4_numbers.py"], {}),
]


@pytest.mark.parametrize("name,args,env", FAST, ids=[f[0] for f in FAST])
def test_fast_simulation(name, args, env):
    assert_lines_present(run_script(args, env), name)


def test_lsm_ruled_stator_is_the_default():
    """Round 13: lsm_config's default is the stator ruled on 2026-11-17 (6 m, 17 turns, 2,088 A, 14,190 coils), and it
    gives the launch's thrust and power; the voltages are those of round 12's final stator."""
    import math
    import lsm_config as cfg
    import lsm_physics as phys
    assert cfg.STATOR == "ruled"
    assert (cfg.TAU_P, cfg.N_STATOR, cfg.N_HTS_LAYERS, cfg.TAPE_WIDTH_MM, cfg.ALPHA_SLED, cfg.V_CORNER) == (6.0, 17, 30, 3.0, 0.4, 5000.0)
    assert cfg.I_PEAK == pytest.approx(2087.6, abs=0.05)                                    # A
    assert cfg.A_ACTIVE_TOTAL == pytest.approx(39_000.0, abs=1e-6)                          # m^2
    assert cfg.N_LIM_SIDES * 3 * cfg.COILS_Q * cfg.N_UNITS == pytest.approx(14_190, abs=0.5)   # coils facing the sled
    K1 = phys.calc_sheet_current_full_pitch(cfg.N_STATOR, cfg.I_PEAK, cfg.COILS_Q, cfg.TAU_P, cfg.K_W)
    assert K1 == pytest.approx(49_323, abs=0.5)                                             # A/m
    F = phys.calc_thrust_current_sheet(K1, cfg.B_SLED_NOMINAL, cfg.DELTA_MAX, cfg.A_ACTIVE_TOTAL)
    assert F / 1e6 == pytest.approx(83.294, abs=0.0005)                                     # MN
    assert F * cfg.V_CORNER / 1e9 == pytest.approx(416.47, abs=0.005)                       # GW
    # the coil at 5 km/s and at 30 km/s: back-EMF, inductive voltage, terminal voltage (kV)
    for v, e_ref, x_ref, v_ref in ((5000.0, 32.47, 18.39, 28.20), (30000.0, 95.14, 37.66, 82.98)):
        s = cfg.V_CORNER / v
        E = phys.calc_EMF(cfg.N_STATOR, v, cfg.B_SLED_NOMINAL * s**cfg.ALPHA_SLED, cfg.W_COIL, cfg.K_W)
        X = 2 * math.pi * v / (2 * cfg.TAU_P) * cfg.L_COIL * cfg.I_PEAK * s**(1 - cfg.ALPHA_SLED)
        assert E / 1e3 == pytest.approx(e_ref, abs=0.005)
        assert X / 1e3 == pytest.approx(x_ref, abs=0.005)
        assert phys.calc_coil_voltage(E, X, cfg.DELTA_MAX) / 1e3 == pytest.approx(v_ref, abs=0.005)
    # the default formula of calc_EMF is the old one
    assert phys.calc_EMF(50, 5000.0, 0.10, 2.0) == pytest.approx(100_000.0, rel=1e-12)


# ---------------------------------------------------------------- 3. the slow simulations, from the reader's own run

SLOW = ["altitude_trade", "sled_passage", "sled_passage_limits", "sled_passage_pulsed", "sled_passage_spread",
        "compensation_fields", "compensation_tracks", "bearing_fem_check", "bearing_design", "bearing_mass",
        "final_bearing", "lim_simulation", "lim_simulation_16MW"]


@pytest.mark.parametrize("name", SLOW)
def test_slow_simulation_log(name):
    assert_lines_present(saved_log(name), name)
