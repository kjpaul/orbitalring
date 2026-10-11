"""Hand estimate (2026-10-12): an air-core coil as wide as the track, facing an iron plate.
Method of images, plate taken as unsaturated with infinite permeability. Not a simulation.
A coil at height g above the plate is attracted to its image at 2g. Each leg carries NI."""
import math
MU0 = 4e-7*math.pi
g, W = 0.100, 1.0
def pull(NI):
    own = MU0*NI**2/(2*math.pi*2*g)                 # leg to its own image, N/m
    d = math.hypot(W, 2*g)
    cross = MU0*NI**2/(2*math.pi*d)*(2*g/d)        # leg to the other leg's image, vertical part, repels
    return 2*(own-cross)
for F in (47.5e3, 95.0e3, 190e3, 349e3):
    NI = math.sqrt(F/pull(1.0))
    Bn_peak = MU0*NI/(2*math.pi*g)                  # largest normal field on the plate
    flux = 2*MU0*NI/(2*math.pi)*math.log(((W/2)**2+g**2)/g**2)   # Wb per m carried sideways between the legs
    Bself = MU0*NI/(2*math.pi*0.05)                 # field at the surface of a 50 mm radius bundle
    print(f"pull {F/1e3:6.1f} kN/m: {NI/1e3:5.0f} kA-turns per leg, peak normal field on plate {Bn_peak:.2f} T, "
          f"sideways flux {flux:.3f} Wb/m -> plate {flux/1.2*1e3:.0f} mm thick at 1.2 T, field at the winding about {Bself:.2f} T")
print("12 mm tape at 0.9 T, 70 K: Ic about 330 A (THEVA data), run at 165 A")
NI = math.sqrt(95.0e3/pull(1.0)); turns = NI/165
print(f"  95 kN/m: {turns:,.0f} turns per leg, {2*turns:,.0f} m of tape per m of track, {2*turns*0.08e-3*0.012*8890:.0f} kg/m of tape, no iron on the casing")
print("iron design of bearing_design_2026-10-12: 3,213 kg/m of laminations per 1 m track, 5,317 per 1.5 m track")
