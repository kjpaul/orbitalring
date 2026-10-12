#!/bin/bash
# LIM deployment at 8 MW and 16 MW per site (about 30 and 10 minutes).  Run from the repository root; logs go under
# ORBITAL_RING_OUT (default runs).  run_all.py runs the same two cases without bash.
export MPLBACKEND=Agg
export ORBITAL_RING_OUT=${ORBITAL_RING_OUT:-runs}
O=$ORBITAL_RING_OUT; mkdir -p $O
python3 lim_simulation.py > $O/lim_simulation.txt 2>&1
python3 -c "import lim_config; lim_config.MAX_SITE_POWER=16e6; import runpy,sys; sys.argv=['lim_simulation.py']; runpy.run_path('lim_simulation.py', run_name='__main__')" > $O/lim_simulation_16MW.txt 2>&1
echo LIMDONE > $O/_limdone
