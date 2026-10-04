#!/bin/bash
# LIM deployment at 8 MW and 16 MW per site (about 10 minutes each)
export MPLBACKEND=Agg; mkdir -p runs
python3 lim_simulation.py > runs/lim_simulation.txt 2>&1
python3 -c "import lim_config; lim_config.MAX_SITE_POWER=16e6; import runpy,sys; sys.argv=['lim_simulation.py']; runpy.run_path('lim_simulation.py', run_name='__main__')" > runs/lim_simulation_16MW.txt 2>&1
echo LIMDONE > runs/_limdone
