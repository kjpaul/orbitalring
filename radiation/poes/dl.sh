base=https://www.ncei.noaa.gov/data/poes-metop-space-environment-monitor/access/l1b/v01r00/2026/metop03
for d in $(seq -w 1 31); do for mth in 08 09; do f=poes_m03_2026${mth}${d}_proc.nc; [ -s $f ] || curl -sS -m 180 -f -O $base/$f || rm -f $f; done; done
echo ALLDONE
