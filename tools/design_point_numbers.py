#!/usr/bin/env python3
"""
The hand-arithmetic rows of the design-point numbers sheet (docs/design_point_numbers_2026-11.md): every quoted number of the
book that follows from ring_altitude.py's values by a line of arithmetic and has no script of its own, printed October
(the fixed-radius quadratic) beside new (the rising quadratic), for both cable directions where the book quotes both.
Chapters 1, 3, 4, 12, 15, 16 and 17.  Nothing here is typed in: every figure is ring_altitude.design() under the two
quadratics and the constants of the book (Young's modulus 500 GPa, the 12,000 kg/m casing, the 1,200 kg/m array).
Book: Orbital Ring Engineering, Volume I: the re-keying of 2026-11 to the rising quadratic.
Run: python tools/design_point_numbers.py
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import math
import ring_altitude as R

E_CNT = R.E_CNT
M_ARRAY = 1200.0


def row(label, old, new, unit="", fmt="{:,.1f}"):
    f = lambda x: x if isinstance(x, str) else fmt.format(x)
    print(f"  {label:70} {f(old):>16} {f(new):>16} {unit}")


def numbers(retro):
    o = R.design(retrograde=retro, quadratic="fixed")
    n = R.design(retrograde=retro, quadratic="rising")
    return o, n


print("Design-point arithmetic, October (fixed-radius quadratic) beside new (rising quadratic); 700 km, 12,000 kg/m casing, 10,686 kg/m hardware")
for retro in (True, False):
    o, n = numbers(retro)
    name = "RETROGRADE" if retro else "PROGRADE"
    print(f"\n-- {name} cable --" + " " * 50 + "October              new")
    row("Chapter 1 Table 1.1, Chapter 3 Table 3.1: structural cable", o["m_struct"], n["m_struct"], "kg/m", "{:,.0f}")
    row("  cable with hardware", o["m_cable"], n["m_cable"], "kg/m", "{:,.0f}")
    row("  ring (cable, hardware, casing)", o["m_ring"], n["m_ring"], "kg/m", "{:,.0f}")
    row("  cable speed", o["v_cable"], n["v_cable"], "m/s")
    row("  speed past the casing", o["v_rel"], n["v_rel"], "m/s")
    row("  tension", o["tension"] / 1e9, n["tension"] / 1e9, "GN")
    row("  cross-section", o["area"], n["area"], "m^2", "{:.2f}")
    row("  side", o["side"], n["side"], "m", "{:.2f}")
    row("  cable mass, structural, around the ring", o["m_struct"] * o["L_ring"] / 1e9, n["m_struct"] * n["L_ring"] / 1e9, "Mt", "{:,.0f}")
    row("  cable mass with hardware, around the ring", o["m_cable"] * o["L_ring"] / 1e9, n["m_cable"] * n["L_ring"] / 1e9, "Mt", "{:,.0f}")
    row("  ring mass", o["M_ring"] / 1e9, n["M_ring"] / 1e9, "Mt", "{:,.0f}")
    row("  deployment energy (closed form)", o["E_deploy"] / 1e18, n["E_deploy"] / 1e18, "EJ", "{:.2f}")
    row("  mass ratio, structural cable per kg of casing (Table 3.1)", o["m_struct"] / R.M_LOAD_M, n["m_struct"] / R.M_LOAD_M, "", "{:.2f}")
    row("  mass ratio, cable with hardware per kg of casing (Ch 12, 16)", o["m_cable"] / R.M_LOAD_M, n["m_cable"] / R.M_LOAD_M, "", "{:.2f}")
    # Chapter 3 checks
    s_free = lambda d: d["m_cable"] / d["area"] * (d["v_cable"] ** 2 - d["v_orbit"] ** 2) / 1e9
    row("Chapter 3: hoop stress with no casing load, (m_c/A)(v_c^2 - v_orb^2)", s_free(o), s_free(n), "GPa", "{:.1f}")
    add0 = lambda d: d["m_cable"] * (d["v_cable"] ** 2 / d["r"] - d["g"]) / d["g_net"] - R.M_LOAD_M
    row("  added load at which the hoop stress reaches zero", add0(o), add0(n), "kg/m", "{:,.0f}")
    I = lambda d: d["side"] ** 4 / 12
    row("  second moment of area, side^4/12 (Eq 3.9)", I(o), I(n), "m^4", "{:.0f}")
    Lc = lambda d: math.sqrt(E_CNT * I(d) / d["tension"])
    row("  static bending length sqrt(EI/T) (Eq 3.10)", Lc(o), Lc(n), "m", "{:.1f}")
    row("  cable centrifugal acceleration v_c^2/r (Eq 4.3)", o["v_cable"] ** 2 / o["r"], n["v_cable"] ** 2 / n["r"], "m/s^2", "{:.3f}")
    row("  net outward acceleration (Eq 4.4)", o["v_cable"] ** 2 / o["r"] - o["g"], n["v_cable"] ** 2 / n["r"] - n["g"], "m/s^2", "{:.3f}")
    # Chapter 4, Table 4.2: reduced mass and collision energy per meter for three casing shares
    for share, lab in ((0.1, "a tenth of the casing"), (0.5, "half the casing"), (1.0, "the whole casing")):
        mu = lambda d: d["m_cable"] * share * R.M_LOAD_M / (d["m_cable"] + share * R.M_LOAD_M)
        row(f"Chapter 4 Table 4.2, {lab}: reduced mass", mu(o), mu(n), "kg/m", "{:,.0f}")
        row(f"  collision energy 1/2 mu v_rel^2", 0.5 * mu(o) * o["v_rel"] ** 2 / 1e9, 0.5 * mu(n) * n["v_rel"] ** 2 / 1e9, "GJ/m", "{:.0f}")
        row(f"  energy per kg of the casing share", 0.5 * mu(o) * o["v_rel"] ** 2 / (share * R.M_LOAD_M) / 1e6, 0.5 * mu(n) * n["v_rel"] ** 2 / (share * R.M_LOAD_M) / 1e6, "MJ/kg", "{:.2f}")
    row("Chapter 4: the cable's excess over circular speed at the break, v_c - v_orb", o["v_cable"] - o["v_orbit"], n["v_cable"] - n["v_orbit"], "m/s", "{:,.0f}")
    row("  ballistic coefficient of the full section broadside, rho side / C_D (C_D 2.2)", R.RHO_CNT * o["side"] / 2.2, R.RHO_CNT * n["side"] / 2.2, "kg/m^2", "{:,.0f}")
    # Chapter 12, Section 12.14: the array's share of the cable
    row("Chapter 12: cable claimed by the 1,200 kg/m array, at the with-hardware ratio", M_ARRAY * o["m_cable"] / R.M_LOAD_M, M_ARRAY * n["m_cable"] / R.M_LOAD_M, "kg/m", "{:,.0f}")
    row("  the same around the ring", M_ARRAY * o["m_cable"] / R.M_LOAD_M * o["L_ring"] / 1e9, M_ARRAY * n["m_cable"] / R.M_LOAD_M * n["L_ring"] / 1e9, "Mt", "{:,.0f}")
    row("  a tenth of the deployment energy", 0.1 * o["E_deploy"] / 1e18, 0.1 * n["E_deploy"] / 1e18, "EJ", "{:.2f}")
    # Chapter 15: the cable's speed past the casing (the drift-in particle speed)
    row("Chapter 15: the cable past the casing interior", o["v_rel"], n["v_rel"], "m/s", "{:,.0f}")
    # Chapter 17: the velocity-squared scaling quotes the Earth tension and cable
    row("Chapter 17: the Earth ring's tension quoted in the scaling", o["tension"] / 1e9, n["tension"] / 1e9, "GN")
    row("  the Earth ring's cable with hardware quoted in the scaling", o["m_cable"], n["m_cable"], "kg/m", "{:,.0f}")

print("\n-- the build orbit (both cables) --")
n = R.design(quadratic="rising")
print(f"  build orbit {n['altitude_build_km']:.2f} km, circular velocity {n['v_build']:,.1f} m/s, v_mom {n['v_mom']:,.1f} m/s, rise {(n['r'] - n['r_build']) / 1e3:.2f} km, strain {R.STRAIN_DESIGN * 100:.1f} %")
