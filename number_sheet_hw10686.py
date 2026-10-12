#!/usr/bin/env python3
"""
Builds the number sheet of the cable hardware change (4,076 to 10,686 kg/m): every quantity the book quotes, old value beside new, from two folders of logs.
Book: The Orbital Ring, Volume IV: The renumbering pass of 2026-10-12 (a tool, not a simulation).
Run: python number_sheet_hw10686.py [--old=FOLDER] [--new=FOLDER] [--out=FILE]
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)

Build the number sheet for the cable hardware change (4,076 -> 10,686 kg/m, 2026-10-12): every quantity the book quotes, old value
(the logs of the run before the change) beside new value (the logs of the run after it), with the log and line each comes from,
and whether it changed.  Writes a Markdown sheet and prints the largest changes.
Usage: python number_sheet_hw10686.py [--old=FOLDER] [--new=FOLDER] [--out=FILE] > log
  --old   folder with the old logs (default: <ORBITAL_RING_OUT>/old)
  --new   folder with the new logs (default: ORBITAL_RING_OUT, default runs)
  --out   the sheet to write (default: <ORBITAL_RING_OUT>/numbers_hw10686.md)
"""
import os
import re
import sys

import ring_altitude as ring

OLD = os.path.join(ring.OUT_DIR, "old")
NEW = ring.OUT_DIR
OUT = os.path.join(ring.OUT_DIR, "numbers_hw10686.md")
for _arg in sys.argv[1:]:
    if _arg in ("--help", "-h"):
        print(__doc__); sys.exit(0)
    elif _arg.startswith("--old="):
        OLD = _arg[6:]
    elif _arg.startswith("--new="):
        NEW = _arg[6:]
    elif _arg.startswith("--out="):
        OUT = _arg[6:]
    else:
        print(f"Unknown option {_arg}"); sys.exit(2)
_cache = {}


def lines(folder, name):
    p = os.path.join(folder, name)
    if p not in _cache:
        try:
            with open(p, encoding="utf-8", errors="replace") as f:
                _cache[p] = f.read().splitlines()
        except OSError:
            _cache[p] = None
    return _cache[p]


def grab(folder, name, pattern, nth=1, after=None, start=0):
    """(groups, line number) of the nth line matching pattern (regex), searching after the line that matches `after` if given."""
    ls = lines(folder, name)
    if ls is None:
        return None, None
    i0 = start
    if after:
        for i, l in enumerate(ls):
            if re.search(after, l):
                i0 = i + 1
                break
        else:
            return None, None
    k = 0
    for i in range(i0, len(ls)):
        m = re.search(pattern, ls[i])
        if m:
            k += 1
            if k == nth:
                return m.groups(), i + 1
    return None, None


def num(s):
    if s is None:
        return None
    try:
        return float(str(s).replace(",", "").replace(">", ""))
    except ValueError:
        return None


ROWS = []        # (section, quantity, unit, file, old_text, old_line, new_text, new_line, old_val, new_val)


def add(section, quantity, unit, name, pattern, group=1, nth=1, after=None, fmt=None, joiner=None):
    go, lo = grab(OLD, name, pattern, nth, after)
    gn, ln = grab(NEW, name, pattern, nth, after)

    def pick(g):
        if g is None:
            return None
        if joiner:
            return joiner.join(g)
        return g[group - 1]
    to, tn = pick(go), pick(gn)
    vo, vn = (num(to), num(tn)) if not joiner else (None, None)
    if fmt and vo is not None and vn is not None:
        to, tn = fmt.format(vo), fmt.format(vn)
    ROWS.append((section, quantity, unit, name, to, lo, tn, ln, vo, vn))


RING_KEYS = [("Radius", "km"), ("Circumference", "km"), ("g at ring", "m/s^2"), ("g net on casing", "m/s^2"), ("Circular velocity", "m/s"),
             ("Escape velocity", "m/s"), ("Orbital period", "min"), ("Ground-sync velocity", "m/s"), ("Delta-v load", "m/s"),
             ("Cable structural", "kg/m"), ("Cable with hardware", "kg/m"), ("Ring", "kg/m"), ("Cable velocity", "m/s"), ("Cable-casing relative", "m/s"),
             ("Cable area", "m^2"), ("Cable side", "m"), ("Cable tension", "GN"), ("Cable centrifugal accel.", "m/s^2"), ("Cable net outward accel.", "m/s^2"),
             ("Casing weight", "kN/m"), ("Cable mass ratio (with hw / casing)", ""), ("Cable mass, structural", "Mt"), ("Cable mass, with hardware", "Mt"),
             ("Casing mass", "Mt"), ("Ring mass", "Mt"), ("Deployment energy", "EJ"), ("LIM sites at 500 m", "")]


def ring_rows(section, name, label):
    for key, unit in RING_KEYS:
        add(section, f"{label}: {key}", unit, name, r"^\s*" + re.escape(key) + r"\s+([\d,.\-]+)")


# ---------------------------------------------------------------- 1. ring tables
ring_rows("1. Ring reference table", "ring_altitude_700.txt", "700 km retrograde")
ring_rows("1b. Ring reference table, 700 km prograde", "ring_altitude_700_prograde.txt", "700 km prograde")
ring_rows("1c. Ring reference table, 250 km", "ring_altitude_250.txt", "250 km")
ring_rows("1d. Ring reference table, 800 km", "ring_altitude_800.txt", "800 km")
# derived from the table: LIM supply frequency at the relative speed, pole pitch 100 m
for folder in (OLD, NEW):
    pass
add("1. Ring reference table", "LIM supply frequency at the relative speed, pole pitch 100 m (computed v_rel / 200 m)", "Hz", "ring_altitude_700.txt",
    r"^\s*Cable-casing relative\s+([\d,.]+)", fmt="{:.2f}")
# fix the frequency row: convert speed to Hz
for i, r in enumerate(ROWS):
    if r[1].startswith("LIM supply frequency"):
        vo, vn = r[8] / 200.0 if r[8] else None, r[9] / 200.0 if r[9] else None
        ROWS[i] = (r[0], r[1], r[2], r[3], f"{vo:.2f}", r[5], f"{vn:.2f}", r[7], vo, vn)

# ---------------------------------------------------------------- 2. derived values (hand numbers, equation logs, chapter 4 scripts)
S2 = "2. Chapters 2 to 4 derived values"
add(S2, "Relative radial acceleration at a break", "m/s^2", "numbers_700_hand.txt", r"relative radial acceleration at a break: .* = ([\d.]+) m/s2")
add(S2, "Time to cross the 100 mm gap", "s", "numbers_700_hand.txt", r"gap time .* = ([\d.]+) s")
add(S2, "Radial speed at contact", "m/s", "numbers_700_hand.txt", r"contact speed ([\d.]+) m/s")
add(S2, "Second moment of area", "m^4", "numbers_700_hand.txt", r"second moment .* = (\d+) m4")
add(S2, "L_char", "m", "numbers_700_hand.txt", r"L_char = .* = ([\d.]+) m")
add(S2, "Lateral kinetic energy at 30 deg", "J", "numbers_700_hand.txt", r"lateral KE at 30 deg: .* = ([\d.e+]+) J")
add(S2, "Net momentum per meter p", "MN s/m", "numbers_700_hand.txt", r"net momentum per metre p = .* = ([\d.]+) MN s/m")
for deg in (5, 10, 20, 30):
    add(S2, f"Table 2.1 lateral force at {deg} deg", "kN/m", "numbers_700_hand.txt", rf"Table 2.1, {deg} deg: F = 2 w p sin i = ([\d.]+) kN/m")
    add(S2, f"Table 2.1 lateral force at {deg} deg as a share of the vertical load", "%", "numbers_700_hand.txt", rf"Table 2.1, {deg} deg: .* = (\d+)% of vertical load")
for k, deg in enumerate((5, 10, 20, 30), start=1):
    add(S2, f"Lateral force from anchor_line_analysis at {deg} deg", "kN/m", "anchor_700_i30.txt", r"Lateral force/m \(kN/m\)\s+([\d.]+)", nth=k)
add(S2, "J2 nodal rate at 5 deg inclination (ring)", "deg/day", "j2_retrograde_i5.txt", r"Ring:\s+(-?[\d.]+) °/day")
add(S2, "Released cable, nominal: perigee x apogee", "km", "debris_orbits.txt", r"swept\s+0 kg/m, recoil\s+\+0: v = ([\d,]+) m/s, perigee (\d+) km, apogee ([\d,]+) km", joiner=" / ")
add(S2, "Slow half apogee (recoil -429)", "km", "debris_orbits.txt", r"swept\s+0 kg/m, recoil\s+-429: .* apogee ([\d,]+) km")
add(S2, "Fast half apogee (recoil +429)", "km", "debris_orbits.txt", r"swept\s+0 kg/m, recoil\s+\+429: .* apogee ([\d,]+) km")
add(S2, "Worst case (slow half plus radial recoil): perigee x apogee, lifetimes full / 0.8 m / 8 cm (Vallado)", "km, years", "altitude_trade.txt",
    r"700 \| slow half \+ radial recoil \|\s+(\d+) x\s+(\d+)\s+\|\s+([\d,>]+) / ([\d,]+) / ([\d,]+)", joiner=" / ")
add(S2, "Same on the solar-cycle average atmosphere", "km, years", "altitude_trade.txt",
    r"700 \| slow half \+ radial recoil \|\s+(\d+) x\s+(\d+)\s+\|\s+([\d,>]+) / ([\d,]+) / ([\d,]+)", nth=2, joiner=" / ")
add(S2, "Worst case with 10% of the casing swept up (Vallado)", "km, years", "altitude_trade.txt",
    r"700 \| same \+ 10% of casing swept\|\s+(\d+) x\s+(\d+)\s+\|\s+([\d,>]+) / ([\d,]+) / ([\d,]+)", joiner=" / ")
add(S2, "Worst case with 25% of the casing swept up (Vallado)", "km, years", "altitude_trade.txt",
    r"700 \| same \+ 25% of casing swept\|\s+(\d+) x\s+(\d+)\s+\|\s+([\d,>]+) / ([\d,]+) / ([\d,]+)", joiner=" / ")
add(S2, "Rod ballistic coefficient: side, B, m/A", "m, kg/m^2", "ch4_numbers.txt", r"side ([\d.]+) Brod ([\d.]+) m/A ([\d.]+)", joiner=" / ")
add(S2, "Chapter 4: 700 km cable total / casing / ring (Mt), structural (kg/m)", "", "ch4_numbers.txt",
    r"700 km: cable total Mt ([\d.]+) casing Mt ([\d.]+) ring Mt ([\d.]+) .* struct ([\d.]+)", joiner=" / ")
add(S2, "Table 4.1, casing 12,000: reduced mass, energy per m, per kg, ring total, Gt TNT", "", "ch4_extra_700.txt",
    r"casing\s+12000: reduced mass ([\d,]+) kg/m\s+E ([\d.]+) GJ/m\s+per kg ([\d.]+) MJ/kg\s+ring total ([\d.e+]+) J = ([\d.]+) Gt TNT = (\d+) x 50 Mt", joiner=" / ")
add(S2, "Chapter 4 prograde comparison: v_cable, v_rel, a_cent, a_net, side", "", "ch4_extra_700.txt",
    r"prograde: v_cable ([\d.]+) v_rel ([\d.]+) .* a_cent ([\d.]+) net ([\d.]+) rel ([\d.]+); side ([\d.]+)", joiner=" / ")
add(S2, "Casing kinetic energy per m, t TNT per m, total Gt TNT", "", "ch4_numbers.txt", r"casing KE per m MJ ([\d.]+) t TNT per m ([\d.]+) total Gt TNT ([\d.]+)", joiner=" / ")
# the equation log values (old log is the 2026-10-04 generator log; new is the rerun)
EQ_OLD = "equations_generator_log.txt"
for key, unit in (("a_cent", "m/s^2"), ("a_net", "m/s^2"), ("t_gap", "s"), ("v_radial", "m/s"), ("I_sec", "m^4"), ("L_char", "m"), ("M_ring", "kg"),
                  ("KE_lat", "J"), ("r_apo_km", "km"), ("h_apo_km", "km")):
    go, lo = grab(OLD, EQ_OLD, rf"^\s*{key}\s+([\d,.e+\-]+)")
    gn, ln = grab(NEW, "equations/generator_log_new.txt", rf"^\s*{key}\s+([\d,.e+\-]+)")
    ROWS.append(("2b. Equation generator values", f"Equation value {key}", unit, "equations_generator_log.txt / equations/generator_log_new.txt",
                 go[0] if go else None, lo, gn[0] if gn else None, ln, num(go[0]) if go else None, num(gn[0]) if gn else None))

# ---------------------------------------------------------------- 3. subsystems
S3 = "3. Subsystems"
add(S3, "LIM deployment time at 8 MW per site", "days", "lim_simulation.txt", r"Deployment time: ([\d.]+) days")
add(S3, "LIM deployment time at 8 MW per site", "years", "lim_simulation.txt", r"Deployment time: [\d.]+ days \(([\d.]+) years\)")
add(S3, "LIM deployment time at 16 MW per site", "days", "lim_simulation_16MW.txt", r"Deployment time: ([\d.]+) days")
add(S3, "Deployment energy delivered (8 MW)", "EJ", "lim_simulation.txt", r"Total energy delivered: ([\d.]+) EJ")
add(S3, "Energy per LIM site", "TJ", "lim_simulation.txt", r"Energy per site: ([\d.]+) TJ")
add(S3, "Cable mass used by lim_simulation", "kg/m", "lim_simulation.txt", r"M_CABLE_M=(\d+) kg/m")
add(S3, "Peak thrust per LIM", "N", "lim_simulation.txt", r"thrust_max\s+([\d.]+) N")
add(S3, "Peak cable temperature", "K", "lim_simulation.txt", r"Max cable temperature: ([\d.]+) K")
add(S3, "Deployment demand, sites x 8 MW", "GW", "power_deployment.txt", r"Total demand\s+([\d.]+) GW")
add(S3, "Array at 600 km, 53 m wide: total generation", "GW", "power_deployment.txt", r"Total generation\s+([\d.]+) GW")
add(S3, "Array orbit-averaged output", "W/m^2", "power_deployment.txt", r"Orbit-averaged output\s+([\d.]+) W/m\^2")
add(S3, "Array peak output (noon)", "W/m^2", "power_deployment.txt", r"Peak output \(noon\)\s+([\d.]+) W/m\^2")
add(S3, "Shadow half-angle at 600 km", "deg", "power_deployment.txt", r"Shadow half-angle\s+([\d.]+) deg")
add(S3, "HVDC peak power flow", "GW", "power_deployment.txt", r"Peak power flow\s+([\d.]+) GW")
add(S3, "HVDC conductor per meter", "kg/m", "power_deployment.txt", r"Conductor mass/m\s+([\d.]+) kg/m")
add(S3, "HVDC total conductor mass", "Mt", "power_deployment.txt", r"Total conductor mass\s+([\d.]+) Mt")
add(S3, "Ring-mounted panels, 8 MW: orbit-averaged output", "W/m^2", "power_deployment_ring_mounted_8MW.txt", r"Orbit-averaged output\s+([\d.]+) W/m\^2")
add(S3, "Ring-mounted panels, 8 MW: panel width needed", "m", "power_deployment_ring_mounted_8MW.txt", r"Width needed to match\s+([\d.]+) m")
add(S3, "Ring-mounted panels, 16 MW: panel width needed", "m", "power_deployment_ring_mounted_16MW.txt", r"Width needed to match\s+([\d.]+) m")
add(S3, "Operations demand", "GW", "power_ops.txt", r"Total demand\s+([\d.]+) GW")
add(S3, "Vertical anchor line: self-weight stress", "GPa", "anchor_700_i0.txt", r"Self-weight stress \(GPa\)\s+([\d.]+)")
add(S3, "Vertical anchor line: taper ratio", "", "anchor_700_i0.txt", r"Taper ratio\s+([\d.]+)")
add(S3, "Anchor lines of 10 cm per station", "", "anchor_700_i0.txt", r"N lines \(d=10cm\) per stn\s+(\d+)")
add(S3, "Anchor combined stress", "GPa", "anchor_700_i0.txt", r"Combined stress \(GPa\)\s+([\d.]+)")
add(S3, "Anchor lateral force per anchor at 30 deg", "GN", "anchor_700_i30.txt", r"F_lateral/anchor \(GN\)\s+([\d.]+)")
add(S3, "Cable transverse wave speed (anchor script, structural mass)", "m/s", "anchor_700_i0.txt", r"Cable transverse wave speed: ([\d,.]+) m/s")
add(S3, "Bearing: Halbach EML capacity (old levitation model)", "kN/m", "levitation.txt", r"Halbach EML:\s+([\d.]+) kN/m")
add(S3, "Bearing: cable velocity used by levitation_simulation", "m/s", "levitation.txt", r"Cable velocity: ([\d.]+) m/s")
add(S3, "Old limiter plate drag at the nominal gap", "N/m", "levitation.txt", r"Drag per meter: ([\d.e\-]+) N/m")
add(S3, "Old limiter plate drag power for the ring", "MW", "levitation.txt", r"Total ring drag power: .*\(([\d.]+) MW\)")
add(S3, "Mass driver, LSM cargo launch to 30 km/s: time", "min", "lsm_simulation.txt", r"Total time:\s+([\d.]+) minutes")
add(S3, "Mass driver, LSM: distance", "km", "lsm_chapter_metrics.txt", r"Distance traveled:\s+([\d.]+) km")
add(S3, "Mass driver, LSM: laps (hand script; its distances are hard-coded from the 2026-10-04 run)", "", "numbers_700_hand.txt", r"LSM launch [\d,.]+ km = ([\d.]+) laps")
for _label, _file, _pat in (("Mass driver, LSM: laps (distance from lsm_chapter_metrics over the circumference)", "lsm_chapter_metrics.txt", r"Distance traveled:\s+([\d.]+) km"),
                            ("Mass driver, LIM stages: laps (distance from md_simulation over the circumference)", "md_simulation.txt", r"Distance\s+(\d+) km")):
    _go, _lo = grab(OLD, _file, _pat); _gn, _ln = grab(NEW, _file, _pat)
    _vo = num(_go[0]) / 44473.246 if _go else None; _vn = num(_gn[0]) / 44473.246 if _gn else None
    ROWS.append((S3, _label, "", _file, f"{_vo:.2f}" if _vo else None, _lo, f"{_vn:.2f}" if _vn else None, _ln, _vo, _vn))
add(S3, "Mass driver, LSM: peak power", "GW", "lsm_simulation.txt", r"Peak power:\s+([\d.]+) GW")
add(S3, "Mass driver, LSM: peak thrust", "MN", "lsm_simulation.txt", r"Peak thrust:\s+([\d.]+) MN")
add(S3, "Peak radial load at 30 km/s (lsm_simulation, inertial frame)", "g", "lsm_simulation.txt", r"Peak g-load:\s+([\d.]+) g")
add(S3, "Mass driver, LIM stages to 30 km/s: time", "min", "md_simulation.txt", r"Total time\s+\d+ s\s+\(([\d.]+) min")
add(S3, "Mass driver, LIM stages: distance", "km", "md_simulation.txt", r"Distance\s+(\d+) km")
add(S3, "Mass driver, LIM stages: laps", "", "md_simulation.txt", r"Number of laps\s+([\d.]+)")
add(S3, "Elevator tether, GEO to 725 km: taper", "", "tether_stress.txt", r"Taper ratio A\(GEO\)/A\(bot\)\s+[\d.]+\s+([\d.]+)")
add(S3, "Elevator tether, GEO to ground: taper", "", "tether_stress.txt", r"Taper ratio A\(GEO\)/A\(bot\)\s+([\d.]+)")
add(S3, "Ring to elevator: lowest rendezvous dv and altitude", "m/s, km", "rendezvous.txt", r"Minimum rendezvous dv: ([\d.]+) m/s at altitude (\d+) km", joiner=" / ")
add(S3, "Free fall from 725 km: peak velocity and altitude", "m/s, km", "freefall.txt", r"Peak velocity:\s+([\d.]+) m/s\s+at\s+([\d.]+) km", joiner=" / ")
add(S3, "Climber descent: total time", "hours", "climber_descent.txt", r"Total time:\s+([\d.]+) hours")
add(S3, "Atmospheric loading: peak drag per meter", "N/m", "atmospheric_loading.txt", r"Peak drag per meter\s+([\d.e+]+) N/m")
add(S3, "Sled passage: neutral speed, zone, inertial (hand / run)", "km/s", "sled_passage.txt", r"Neutral speed, zone, inertial \(km/s\)\s+([\d.]+)\s+([\d.]+)", joiner=" / ")
add(S3, "Sled passage: peak radial load (g) hand / run / casing frame", "g", "sled_passage.txt", r"Peak radial load \(g\)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", joiner=" / ")
add(S3, "Sled passage: radial load at 30 km/s in the casing frame", "g", "sled_passage.txt", r"radial load is ([\d.]+) g")
add(S3, "Spread pull (rule 4): pull at 1 km/s, excursion", "kN/m, mm", "sled_passage_spread.txt", r"^\s+1\s+[\d.]+\s+[\d.]+\s+(\d+)\s+[\d.]+\s+\d+\s+[\d.]+\s+([\d.]+)", after=r"Rule 4: fixed T = 50 ms", joiner=" / ")
add(S3, "Spread pull (rule 4): pull at 30 km/s, excursion", "kN/m, mm", "sled_passage_spread.txt", r"^\s+30\s+[\d.]+\s+[\d.]+\s+(\d+)\s+[\d.]+\s+\d+\s+[\d.]+\s+([\d.]+)", after=r"Rule 4: fixed T = 50 ms", joiner=" / ")
add(S3, "Compensation tracks (2026-10-11 layout): pull at 30 km/s, movement", "kN/m, mm", "compensation_tracks.txt", r"^\s+30\s+\d+\s+lift\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+\d+\s+\d+\s+[\d,]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)", joiner=" / ")
add(S3, "Bearing design (iron, 2026-10-12): planned pass at 30 km/s, pull, movement", "kN/m, mm", "bearing_design.txt", r"^\s+30\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+$", after=r"F\.1  Sled passage", joiner=" / ")
add(S3, "Bearing design (iron): planned pass at 1 km/s, pull, movement", "kN/m, mm", "bearing_design.txt", r"^\s+1\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+$", after=r"F\.1  Sled passage", joiner=" / ")
add(S3, "Bearing mass: set 4 air, casing / hardware / cable / x book / tension", "", "bearing_mass.txt", r"^\s+4\s+1\.50\s+air\s+([\d,]+)\s+([\d,]+)\s+([\d,]+)\s+([\d.]+)\s+(\d+)", joiner=" / ")
add(S3, "Final bearing: caps lift / pair, load / pair, support one track", "kN/m", "final_bearing.txt", r"Caps: lift (\d+) kN/m per track \((\d+) pair\), load (\d+) \((\d+)\); support one track alone (\d+)", joiner=" / ")
add(S3, "Final bearing: smallest gap in any single failure", "mm", "final_bearing.txt", r"Smallest gap in any single failure: ([\d.]+) mm")
add(S3, "Final bearing: fastest launch with one lift track out", "km/s", "final_bearing.txt", r"fastest launch inside 20 mm with one lift track out: (\d+) km/s")
add(S3, "Final bearing: planned pass at 30 km/s, pull, movement", "kN/m, mm", "final_bearing.txt", r"^\s+30\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+$", after=r"F\.1  Sled passage", joiner=" / ")
add(S3, "Final bearing: planned pass at 1 km/s, pull, movement", "kN/m, mm", "final_bearing.txt", r"^\s+1\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+[\d.]+$", after=r"F\.1  Sled passage", joiner=" / ")
# ---------------------------------------------------------------- 4. array and panels (hand)
S4 = "4. The array at 600 km and the ring panels"
add(S4, "Sunlit share at equinox, array at 600 km", "%", "numbers_700_hand.txt", r"array sunlit share at equinox: .* = ([\d.]+)%")
add(S4, "Sunlit share at the ring", "%", "numbers_700_hand.txt", r"sunlit share at the ring: ([\d.]+)%")
add(S4, "Ring panel output after deployment, 8 MW / 16 MW", "GW", "numbers_700_hand.txt", r"LIM sites \d+; x 8 MW = ([\d.]+) GW; x 16 MW = ([\d.]+) GW", joiner=" / ")
add(S4, "Latitude reach of a 45 deg line and its length", "deg, km", "numbers_700_hand.txt", r"latitude reach of a 45 deg line: \+-([\d.]+) deg, length ([\d,.]+)", joiner=" / ")
add(S4, "Anchor transit at 0.1 g: time and peak speed", "s, m/s", "numbers_700_hand.txt", r"anchor transit at 0.1 g: .* = (\d+) s = [\d.]+ min, peak (\d+) m/s", joiner=" / ")

# ---------------------------------------------------------------- write
def fmt_change(vo, vn):
    if vo is None or vn is None:
        return "see text" if (vo is None and vn is None) else "n/a"
    if vo == vn:
        return "unchanged"
    if vo == 0:
        return "changed"
    return f"{(vn - vo) / abs(vo) * 100:+.2f} %"


sections = []
for r in ROWS:
    if r[0] not in sections:
        sections.append(r[0])
with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write("# The 700 km numbers with the cable hardware at 10,686 kg/m (2026-10-12)\n\n")
    f.write("Design point: ring at 700 km, retrograde cable, suspended array at 600 km; casing 12,000 kg/m; cable hardware 10,686 kg/m (LIM plates 1,872 kg/m plus four ")
    f.write("ferromagnetic plates 200 mm thick and 1.4 m wide, 8,814 kg/m; was 4,076 kg/m with two 100 mm plates). This sheet replaces the sheet of 2026-10-04 ")
    f.write("for the renumbering pass.\n\n")
    f.write("Old values come from the logs of the run before the change, new values from the same-named logs of the run after it ")
    f.write("(commands in that folder's `run_700.sh` and `run_700_lim.sh`, the extra scripts in the report). 'Log, line' gives the file and the line number in the new ")
    f.write("folder (old line in parentheses when it differs). Built by `Python/orbitalring/number_sheet_hw10686.py`.\n\n")
    for sec in sections:
        f.write(f"## {sec}\n\n")
        f.write("| Quantity | Unit | Old | New | Change | Log, line |\n|---|---|---:|---:|---:|---|\n")
        for (s, q, unit, name, to, lo, tn, ln, vo, vn) in ROWS:
            if s != sec:
                continue
            to_s = to if to is not None else "not found"
            tn_s = tn if tn is not None else "not found"
            lines_s = f"`{name}` {ln}" + (f" ({lo})" if lo and ln and lo != ln else "") if ln else f"`{name}` not found"
            ch = fmt_change(vo, vn) if (vo is not None or vn is not None) else ("unchanged" if to == tn and to is not None else ("changed" if to != tn and to is not None and tn is not None else "n/a"))
            f.write(f"| {q} | {unit} | {to_s} | {tn_s} | {ch} | {lines_s} |\n")
        f.write("\n")
    # changes summary
    numeric = [(abs((vn - vo) / vo) * 100, (vn - vo) / vo * 100, q, to, tn, unit, name, ln) for (s, q, unit, name, to, lo, tn, ln, vo, vn) in ROWS
               if vo not in (None, 0) and vn is not None and not q.startswith("800 km") and not q.startswith("250 km") and not q.startswith("700 km prograde")]
    numeric.sort(reverse=True)
    f.write("## 5. The largest changes (700 km retrograde and the subsystem logs)\n\n")
    f.write("| Quantity | Old | New | Change |\n|---|---:|---:|---:|\n")
    for a, pct, q, to, tn, unit, name, ln in numeric[:10]:
        f.write(f"| {q} ({unit}) | {to} | {tn} | {pct:+.2f} % |\n")
    f.write("\nEverything that moved by more than 5 percent:\n\n")
    big = [x for x in numeric if x[0] > 5.0]
    if big:
        for a, pct, q, to, tn, unit, name, ln in big:
            f.write(f"- {q} ({unit}): {to} to {tn}, {pct:+.1f} % (`{name}` line {ln})\n")
    else:
        f.write("- nothing moved by more than 5 percent\n")
    unchanged = sum(1 for x in numeric if abs(x[1]) < 1e-9)
    f.write(f"\n{len(numeric)} numeric quantities compared for the 700 km retrograde case and the subsystems; {unchanged} unchanged, {len(big)} moved by more than 5 percent.\n")
print(f"wrote {OUT}: {len(ROWS)} rows")
missing = [(s, q, name) for (s, q, unit, name, to, lo, tn, ln, vo, vn) in ROWS if tn is None]
print(f"{len(missing)} rows with no new value:")
for s, q, name in missing:
    print(f"  {q}  [{name}]")
print("\nTen largest changes:")
numeric = [(abs((vn - vo) / vo) * 100, (vn - vo) / vo * 100, q, to, tn) for (s, q, unit, name, to, lo, tn, ln, vo, vn) in ROWS
           if vo not in (None, 0) and vn is not None and not q.startswith("800 km") and not q.startswith("250 km") and not q.startswith("700 km prograde")]
numeric.sort(reverse=True)
for a, pct, q, to, tn in numeric[:10]:
    print(f"  {q}: {to} -> {tn} ({pct:+.2f} %)")
print("\nMoved by more than 5 percent:")
for a, pct, q, to, tn in numeric:
    if a > 5:
        print(f"  {q}: {to} -> {tn} ({pct:+.1f} %)")
