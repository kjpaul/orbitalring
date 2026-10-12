"""Reference check for The Orbital Ring, Volume IV: anyone with the repository and the packages of requirements.txt can
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
    """700 km, retrograde, hardware 10,686 kg/m: the values of Chapter 3 and Table 8.1."""
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
    ("lsm_chapter_metrics", ["lsm_chapter_metrics.py"], {}),
    ("lsm_n_optimization", ["lsm_n_optimization_demo.py"], {}),
    ("lsm_simulation", ["lsm_simulation.py"], {}),
    ("md_simulation", ["md_simulation.py", "--no-graphs"], {}),
    ("md_trade_study", ["md_trade_study.py"], {}),
    ("debris_simulation", ["debris_simulation.py"], {}),
    ("debris_orbits", ["debris_orbits.py"], {}),
    ("ch4_numbers", ["ch4_numbers.py"], {}),
]


@pytest.mark.parametrize("name,args,env", FAST, ids=[f[0] for f in FAST])
def test_fast_simulation(name, args, env):
    assert_lines_present(run_script(args, env), name)


# ---------------------------------------------------------------- 3. the slow simulations, from the reader's own run

SLOW = ["altitude_trade", "sled_passage", "sled_passage_limits", "sled_passage_pulsed", "sled_passage_spread",
        "compensation_fields", "compensation_tracks", "bearing_fem_check", "bearing_design", "bearing_mass",
        "final_bearing", "lim_simulation", "lim_simulation_16MW"]


@pytest.mark.parametrize("name", SLOW)
def test_slow_simulation_log(name):
    assert_lines_present(saved_log(name), name)
