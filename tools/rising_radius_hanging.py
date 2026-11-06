#!/usr/bin/env python3
"""
What hangs from a ring that rises: the altitude of a point a fixed distance below the casing through the rising-radius
deployment (the suspended array of Chapter 12 and the spreader arms), from the daily profile lim_simulation.py --rising-radius
writes (rising_radius_daily.csv in the output folder).
Book: Orbital Ring Engineering, Volume I: Chapter 11 (Section 11.1, 150 km as the lower limit for what hangs during the
deployment), Chapter 12 (Section 12.11, the array's lowering and the arms held above the air).
Run: python tools/rising_radius_hanging.py [rising_radius_daily.csv]     (writes rising_radius_hanging.csv next to it)
Environment: ORBITAL_RING_OUT (output folder, default runs)

A point d below the casing is at altitude h(t) - d.  The depths are 20, 50 and 100 km (the array support lines are about
100 km, Section 12.4) and the spreader arm as Section 12.11 hangs it during the deployment: vertical from the casing with its
lower end held at 150 km, which from the 700 km ring is 550 km of arm (tools/anchor_lines_deployment.py pays out the same
550 km).  The full arm, 100 km to the array and 893 km below the kink (Section 12.3), is longer than the ring is high and is
never hung straight down.  For each depth the days on which the point would be below 150 km are listed, and for each day the
longest hang that stays above 150 km, h(t) - 150 km, which is the pay-out schedule the lines and arms would follow if they
are let down as far as the air allows from the first day.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import csv
import os
import sys

import numpy as np

import ring_altitude as ring

FLOOR_KM = 150.0                      # Section 11.1: the lines are let down to about 150 km early, no lower
DEPTHS_KM = [20.0, 50.0, 100.0]       # fixed distances below the casing
ARM_HANG_KM = ring.ALTITUDE / 1e3 - FLOOR_KM   # the spreader arm as hung during the deployment: casing to 150 km from the 700 km ring


def main():
    src = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(ring.out_dir(), "rising_radius_daily.csv")
    rows = list(csv.DictReader(open(src, newline="")))
    day = np.array([float(r["day"]) for r in rows]); h = np.array([float(r["altitude_km"]) for r in rows]); phase = np.array([int(r["phase"]) for r in rows])
    t2 = day[phase >= 2][0] if (phase >= 2).any() else None
    print(f"Rising-radius deployment profile: {os.path.basename(src)}, {len(rows)} rows, day 0 at {h[0]:.2f} km to day {day[-1]:.1f} at {h[-1]:.2f} km"
          + (f"; the casing at ground speed from day {t2:.1f} ({float(h[phase >= 2][0]):.2f} km)." if t2 is not None else "."))
    print(f"Lower limit for what hangs: {FLOOR_KM:.0f} km (Section 11.1).  Spreader arm as hung during the deployment: {ARM_HANG_KM:.0f} km below the casing (Section 12.11).")
    print()
    depths = DEPTHS_KM + [ARM_HANG_KM]
    names = [f"{d:.0f} km below" for d in DEPTHS_KM] + [f"arm, {ARM_HANG_KM:.0f} km"]
    print(f"{'day':>6} {'ring km':>8} " + " ".join(f"{n:>14}" for n in names) + f" {'max hang km':>12}")
    marks = [0, 1, 10, 30, 60, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700]
    for d_ in marks:
        if d_ > day[-1]:
            continue
        hk = float(np.interp(d_, day, h))
        print(f"{d_:6d} {hk:8.2f} " + " ".join(f"{hk - d:14.2f}" for d in depths) + f" {hk - FLOOR_KM:12.2f}")
    if day[-1] not in marks:
        hk = h[-1]
        print(f"{day[-1]:6.1f} {hk:8.2f} " + " ".join(f"{hk - d:14.2f}" for d in depths) + f" {hk - FLOOR_KM:12.2f}")
    print()
    for n, d in zip(names, depths):
        below = day[h - d < FLOOR_KM]
        if len(below) == 0:
            print(f"  {n:>16}: never below {FLOOR_KM:.0f} km (lowest {h.min() - d:.1f} km on day {day[np.argmin(h)]:.0f}).")
        else:
            above = day[h - d >= FLOOR_KM]
            first_ok = f"above {FLOOR_KM:.0f} km from day {above[0]:.0f} on" if len(above) else f"never above {FLOOR_KM:.0f} km during the run"
            ground = day[h - d < 0]
            print(f"  {n:>16}: below {FLOOR_KM:.0f} km on days {below[0]:.0f} to {below[-1]:.0f} ({len(below)} days; {first_ok})"
                  + (f"; below the ground on days {ground[0]:.0f} to {ground[-1]:.0f}" if len(ground) else "") + ".")
    print()
    for target in (400.0, 450.0, 500.0, 550.0, 600.0):
        ok = day[h - FLOOR_KM >= target]
        print(f"  a hang of {target:.0f} km stays above {FLOOR_KM:.0f} km from day {ok[0]:.0f} on (ring at {float(np.interp(ok[0], day, h)):.1f} km)." if len(ok) else
              f"  a hang of {target:.0f} km never stays above {FLOOR_KM:.0f} km during the run.")
    out = os.path.join(os.path.dirname(os.path.abspath(src)), "rising_radius_hanging.csv")
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["day", "ring_altitude_km", "phase"] + [f"alt_{d:.0f}km_below" for d in DEPTHS_KM] + [f"alt_arm_{ARM_HANG_KM:.0f}km_below", "max_hang_above_150km"])
        for k in range(len(day)):
            wr.writerow([f"{day[k]:.3f}", f"{h[k]:.3f}", phase[k]] + [f"{h[k] - d:.3f}" for d in depths] + [f"{h[k] - FLOOR_KM:.3f}"])
    print(f"\nWritten {os.path.basename(out)} ({len(day)} rows) next to the daily table.")


if __name__ == "__main__":
    main()
