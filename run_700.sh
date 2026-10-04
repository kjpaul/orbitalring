#!/bin/bash
# Reruns every subsystem script at the 2026-10-04 design point: ring 700 km retrograde, suspended array at 600 km.
# Run from Python/orbitalring. Logs go to runs/. The same commands reproduce runs_800km_2026-10-03 with ORBITAL_RING_ALT_KM=800 ORBITAL_RING_ARRAY_OFFSET_KM=50.
export MPLBACKEND=Agg
O=runs; mkdir -p $O
python3 ring_altitude.py > $O/ring_altitude_700.txt
ORBITAL_RING_CABLE=prograde python3 ring_altitude.py > $O/ring_altitude_700_prograde.txt
ORBITAL_RING_ALT_KM=250 python3 ring_altitude.py > $O/ring_altitude_250.txt
ORBITAL_RING_ALT_KM=800 python3 ring_altitude.py > $O/ring_altitude_800.txt
for i in 0 5 10 20 30; do python3 anchor_line_analysis.py --inclination=$i > $O/anchor_700_i$i.txt 2>&1; done
for i in 5 10 20 30; do python3 j2_simulation.py --inclination=$i > $O/j2_retrograde_i$i.txt 2>&1; python3 j2_simulation.py --inclination=$i --cable=prograde > $O/j2_prograde_i$i.txt 2>&1; done
python3 power_simulation.py > $O/power_deployment.txt 2>&1
python3 power_simulation.py --demand=ops > $O/power_ops.txt 2>&1
ORBITAL_RING_ARRAY_OFFSET_KM=0 python3 power_simulation.py > $O/power_deployment_ring_mounted_8MW.txt 2>&1
ORBITAL_RING_ARRAY_OFFSET_KM=0 python3 power_simulation.py --power=16 > $O/power_deployment_ring_mounted_16MW.txt 2>&1
python3 levitation_simulation.py all > $O/levitation.txt 2>&1
python3 tether_stress.py > $O/tether_stress.txt 2>&1
python3 rendezvous_analysis.py > $O/rendezvous.txt 2>&1
python3 freefall_trajectory.py > $O/freefall.txt 2>&1
python3 atmospheric_loading.py > $O/atmospheric_loading.txt 2>&1
python3 climber_descent.py > $O/climber_descent.txt 2>&1
python3 lsm_chapter_metrics.py > $O/lsm_chapter_metrics.txt 2>&1
python3 lsm_n_optimization_demo.py > $O/lsm_n_optimization.txt 2>&1
python3 lsm_simulation.py > $O/lsm_simulation.txt 2>&1
python3 md_simulation.py --no-graphs > $O/md_simulation.txt 2>&1
python3 md_trade_study.py > $O/md_trade_study.txt 2>&1
python3 debris_simulation.py > $O/debris_simulation.txt 2>&1
python3 debris_orbits.py > $O/debris_orbits.txt 2>&1
python3 ch4_numbers.py > $O/ch4_numbers.txt 2>&1
python3 altitude_trade.py > $O/altitude_trade.txt 2>&1
echo ALLDONE > $O/_done
