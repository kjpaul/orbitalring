#!/bin/bash
# The deployment of the ring with a 6,000 kg/m casing (Chapter 20, Table 20.2; 2026-10-28), at 8 MW per site and the default
# winding of lim_config.py (flat 3 mm, 208 turns, Brandt-Indenbom hysteresis).  The cable hardware is set through the plate
# thicknesses: D_LEV (four bearing plates 1.4 m wide, iron 7,870 kg/m3) and T_PLATE (two LIM plates 1.2 m wide, gamma-TiAl
# 3,900 kg/m3).  Phase 1: cases (a) and (b) and the LIM plate test; phase 2 (after casing_6000.py --bearing picks the plates):
# case (c) at 8 and 16 MW.  Usage: bash run_700_casing_6000.sh 1   |   bash run_700_casing_6000.sh 2 <kept|halved> <D_LEV m> <T_PLATE m>
export MPLBACKEND=Agg
LOG=../../Vol_III/runs_700km_hw10686_2026-10-12
run() {  # label, D_LEV m, T_PLATE m, power W
  mkdir -p runs/casing_6000/$1
  ORBITAL_RING_OUT=runs/casing_6000/$1 timeout 7200 python -c "import lim_config as c; c.configure(D_LEV=$2, T_PLATE=$3, MAX_SITE_POWER=$4); import runpy, sys; sys.argv=['lim_simulation.py', '--m_load=6000']; runpy.run_path('lim_simulation.py', run_name='__main__')" > $LOG/casing_6000_lim_$1.txt 2>&1
}
if [ "$1" == "1" ]; then
  run a_8MW            0.20 0.20 8e6 &     # (a) hardware unchanged, 10,686 kg/m
  run b_8MW            0.10 0.10 8e6 &     # (b) hardware halved with the casing, 5,343 kg/m
  run b_limplate200_8MW 0.10 0.20 8e6 &    # (b) with the LIM plate kept at 200 mm: what the plate thickness does to the thrust and the heating
  wait
  echo PHASE1_DONE
else
  run c_$2_8MW  $3 $4 8e6 &                 # (c) plates sized for the lighter casing: $2 = kept or halved (the standing load)
  run c_$2_16MW $3 $4 16e6 &
  wait
  echo PHASE2_DONE
fi
