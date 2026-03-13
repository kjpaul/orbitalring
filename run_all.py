#!/usr/bin/env python3
"""
Space Elevator Tether Dynamics — Master Script
================================================

Runs all five analysis sections and generates all plots and tables.

Usage:
    python run_all.py                  Run all sections, save all graphs
    python run_all.py --show           Display graphs interactively
    python run_all.py --save           Save graphs to files (default)
    python run_all.py --section=1      Run only section 1
    python run_all.py --section=1,3    Run sections 1 and 3
    python run_all.py -h / --help      Show this help

Sections:
    1  Tether stress and taper profiles     (tether_stress.py)
    2  Free-fall trajectory from 275 km     (freefall_trajectory.py)
    3  Climber descent from GEO             (climber_descent.py)
    4  Atmospheric loading on tether        (atmospheric_loading.py)
    5  Rendezvous velocity analysis         (rendezvous_analysis.py)

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import os
import subprocess
import time


# =============================================================================
# SECTION DEFINITIONS
# =============================================================================

SECTIONS = {
    1: ("Tether Stress & Taper Profiles",    "tether_stress.py"),
    2: ("Free-fall Trajectory (275 km)",      "freefall_trajectory.py"),
    3: ("Climber Descent from GEO",           "climber_descent.py"),
    4: ("Atmospheric Loading",                "atmospheric_loading.py"),
    5: ("Rendezvous Velocity Analysis",       "rendezvous_analysis.py"),
}


# =============================================================================
# MAIN
# =============================================================================

def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')

    # Parse arguments
    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__)
        return

    show_mode = "--show" if "--show" in sys.argv else "--save"
    sections_to_run = list(SECTIONS.keys())

    for arg in sys.argv[1:]:
        if arg.startswith("--section="):
            parts = arg.split("=")[1].split(",")
            sections_to_run = [int(p.strip()) for p in parts]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    python_exe = sys.executable

    streq = "=" * 70
    print(f"\n{streq}")
    print("SPACE ELEVATOR TETHER DYNAMICS — FULL ANALYSIS")
    print(f"{streq}")
    print(f"  Sections to run: {sections_to_run}")
    print(f"  Graph mode:      {show_mode}")
    print(f"{streq}\n")

    total_start = time.time()
    results = {}

    for sec_num in sections_to_run:
        if sec_num not in SECTIONS:
            print(f"  WARNING: Unknown section {sec_num}, skipping.")
            continue

        title, script = SECTIONS[sec_num]
        script_path = os.path.join(script_dir, script)

        if not os.path.exists(script_path):
            print(f"  ERROR: {script} not found, skipping section {sec_num}.")
            continue

        print(f"\n{'─' * 70}")
        print(f"  SECTION {sec_num}: {title}")
        print(f"  Running: {script}")
        print(f"{'─' * 70}\n")

        start = time.time()
        try:
            result = subprocess.run(
                [python_exe, script_path, "all", show_mode],
                cwd=script_dir,
                timeout=300,
            )
            elapsed = time.time() - start
            results[sec_num] = {
                "title": title,
                "returncode": result.returncode,
                "elapsed": elapsed,
            }
            if result.returncode != 0:
                print(f"\n  WARNING: {script} exited with code {result.returncode}")
        except subprocess.TimeoutExpired:
            elapsed = time.time() - start
            results[sec_num] = {
                "title": title,
                "returncode": -1,
                "elapsed": elapsed,
            }
            print(f"\n  ERROR: {script} timed out after 300s")
        except Exception as e:
            elapsed = time.time() - start
            results[sec_num] = {
                "title": title,
                "returncode": -2,
                "elapsed": elapsed,
            }
            print(f"\n  ERROR: {script} failed: {e}")

    total_elapsed = time.time() - total_start

    # Summary
    print(f"\n{streq}")
    print("EXECUTION SUMMARY")
    print(streq)
    for sec_num in sections_to_run:
        if sec_num in results:
            r = results[sec_num]
            status = "OK" if r["returncode"] == 0 else f"FAIL ({r['returncode']})"
            print(f"  Section {sec_num}: {r['title']:40s} {status:8s} {r['elapsed']:6.1f}s")
    print(f"{'─' * 70}")
    print(f"  {'Total':49s} {total_elapsed:6.1f}s")
    print(f"{streq}\n")


if __name__ == "__main__":
    main()
