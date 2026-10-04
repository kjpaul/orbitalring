import numpy as np
from sweetspot_run import *
alt = 800
print("Under the cable at 800 km, clean arc, GCR dose equivalent mSv/yr (model), point on the centre line, room half-width 2 m, 1 m above floor")
print(" wall and floor g/cm2:      5     10     20     30     50")
for h in (0.5, 1.0, 2.0, 3.0):
    print(f"  h = {h} m:           " + " ".join(f"{gcr_point(alt, thickness(h, t_wall=tw, t_floor=tw)[0])*1e3:6.1f}" for tw in (5,10,20,30,50)))
print(" narrower room (half-width 1 m), h = 1 m, walls 20: ", round(gcr_point(alt, thickness(1.0, a=1.0, t_wall=20, t_floor=20)[0])*1e3, 1))
print(" no cable, uniform sphere: " + " ".join(f"{tw}: {gcr_point(alt, np.full(len(U), float(tw)))*1e3:.1f}" for tw in (5,10,20,30,50)))
# contributions at h=1, walls 20
t, hit = thickness(1.0, t_wall=20, t_floor=20); rc = stormer(alt); eb = earth_blocked(alt)
d = gcr_at(t, rc)*(~eb)*dOm*1e3
print(f" h=1 m, walls 20: through the cable {d[hit].sum():.1f}, through walls and floor {d[~hit].sum():.1f} mSv/yr; cable directions with less than 600 g/cm2 of carbon in the path: {100*((t<600)&hit).mean():.1f}% of the sky")
