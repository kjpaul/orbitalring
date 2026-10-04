python3 - <<'P' > lowE_jobs.txt
import numpy as np
ps=np.loadtxt('pstar_water.txt'); 
for e in np.geomspace(12,2000,72)[:32]:
    R=np.exp(np.interp(np.log(e),np.log(ps[:,0]),np.log(ps[:,4])))
    print(f"{e:.2f} {1.25*R:.4f}")
P
while read e d; do DEPTH_CM=$d NBIN=250 python3 g4_depthdose.py QGSP_BIC_HP proton out_trap_lowE 6000:$e 2>&1 | grep DONE; done < lowE_jobs.txt
