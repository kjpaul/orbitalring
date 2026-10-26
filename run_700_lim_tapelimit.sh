#!/bin/bash
# The deployment on the measured tape limit (2026-10-26): 208 turns at 8 and 16 MW per site, and the turns trade at 100, 150 and
# 300 turns, four runs at a time (about 10 minutes each).  Logs to Vol_III/runs_700km_hw10686_2026-10-12/; each run gets its own
# output folder under runs/tapelimit/ for its CSV and graphs.
export MPLBACKEND=Agg
LOG=../../Vol_III/runs_700km_hw10686_2026-10-12
run() {  # name, turns, power (W)
  mkdir -p runs/tapelimit/$1
  ORBITAL_RING_OUT=runs/tapelimit/$1 timeout 5400 python3 -c "import lim_config as c; c.configure(TAPE_LIMIT_MODE='measured', N_TURNS=$2, MAX_SITE_POWER=$3); import runpy, sys; sys.argv=['lim_simulation.py']; runpy.run_path('lim_simulation.py', run_name='__main__')" > $LOG/$1.txt 2>&1
}
run lim_simulation_tapelimit_N100_8MW 100 8e6 &
run lim_simulation_tapelimit_N100_16MW 100 16e6 &
run lim_simulation_tapelimit_N150_8MW 150 8e6 &
run lim_simulation_tapelimit_N150_16MW 150 16e6 &
wait
run lim_simulation_tapelimit 208 8e6 &
run lim_simulation_tapelimit_16MW 208 16e6 &
run lim_simulation_tapelimit_N300_8MW 300 8e6 &
run lim_simulation_tapelimit_N300_16MW 300 16e6 &
wait
echo TAPELIMIT_DONE
