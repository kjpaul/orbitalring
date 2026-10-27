#!/bin/bash
# The deployment on the Brandt-Indenbom hysteresis model (2026-10-27) for the candidates the trade table picks out, at 8 and 16 MW per
# site, four runs at a time (about 10 minutes each).  Each candidate's winding constants are those tools/lim_hysteresis_model.py
# prints (the coil cap at the 50 percent rule and the worst tape's limit, field and angle at that cap).
export MPLBACKEND=Agg
LOG=../../Vol_III/runs_700km_hw10686_2026-10-12
run() {  # label, layout, width mm, layers, N, cap A, worst tape A, B T, theta deg, power W
  mkdir -p runs/brandt/$1_${10}
  ORBITAL_RING_OUT=runs/brandt/$1_${10} timeout 5400 python3 -c "import lim_config as c; c.configure(TAPE_LIMIT_MODE='measured', HYSTERESIS_MODEL='brandt', WINDING_LAYOUT_MEASURED='$2', HTS_TAPE_WIDTH_MM=$3, HTS_TAPE_LAYERS_MEASURED=$4, N_TURNS=$5, I_COIL_LIMIT_MEASURED=$6, I_TAPE_LIMIT_MEASURED=$7, B_TAPE_WORST_MEASURED=$8, THETA_TAPE_WORST_MEASURED=$9, MAX_SITE_POWER=${10}); import runpy, sys; sys.argv=['lim_simulation.py']; runpy.run_path('lim_simulation.py', run_name='__main__')" > $LOG/lim_simulation_brandt_$1_${11}.txt 2>&1
}
run flat3mm_N208    flat    3 2 208 183.1 183.1 0.072 7.6  8e6  8MW &
run flat6mm_N208    flat    6 1 208 202.8 405.5 0.043 0.0  8e6  8MW &
run flat3mm_N150    flat    3 2 150 184.3 184.3 0.070 7.9  8e6  8MW &
run pancake3mm_N208 pancake 3 3 208 219.6 146.4 0.186 70.9 8e6  8MW &
wait
run flat3mm_N208    flat    3 2 208 183.1 183.1 0.072 7.6  16e6 16MW &
run flat6mm_N208    flat    6 1 208 202.8 405.5 0.043 0.0  16e6 16MW &
run flat3mm_N150    flat    3 2 150 184.3 184.3 0.070 7.9  16e6 16MW &
run pancake3mm_N208 pancake 3 3 208 219.6 146.4 0.186 70.9 16e6 16MW &
wait
echo BRANDT_DONE
