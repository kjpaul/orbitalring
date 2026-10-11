"""Hand estimates for Chapter 8 design choices (2026-10-10). Not a simulation.
Copper against HTS tape for the compensation coils, and passive against active gap protection.
Inputs are from sled_passage_spread_2026-10-10.md, levitation.txt (700 km run), Table 6.2 and Chapter 18."""
import math
MU0 = 4e-7*math.pi
L_RING = 44.473e6          # m
G_NET, M_CASING = 7.918, 12000.0
W = M_CASING*G_NET         # N/m casing weight
GAP = 0.100
NI = 1.2*2*GAP/MU0         # ampere-turns per pole pair for 1.2 T across two gaps
PITCH = 0.5                # m, winding layout of the loss run [UNVERIFIED]
tracks = {"lift": 0.76, "load": 0.62}
print(f"NI = {NI:,.0f} A-turns, casing weight {W/1e3:.1f} kN/m")
print("\n-- copper --")
RHO_CU, D_CU, C_CU = 1.68e-8, 8960.0, 385.0
for j in (50e6, 10e6):
    tot = 0
    for n, w in tracks.items():
        Lt = 2*(w+PITCH); m = NI/j*Lt*D_CU; tot += m
        print(f" j {j/1e6:.0f} A/mm2 {n}: turn {Lt:.2f} m, {m:.0f} kg/m")
    print(f"   both tracks {tot:.0f} kg/m = {tot/M_CASING*100:.1f}% of casing, {tot*L_RING/1e9:.1f} Mt around the ring; "
          f"heating at full field {RHO_CU*j*j/(D_CU*C_CU):.2f} K/s")
print("\n-- HTS tape, 12 mm, 70 K --")
IC_1T = 316.0              # A, Table 6.2 at 1 T
FRAC = 0.5
I_op = IC_1T*FRAC
T_TAPE, W_TAPE, D_TAPE = 0.08e-3, 0.012, 8890.0   # Ch 6 thickness; Hastelloy density
COP = 0.0106               # Equation 18.1
tot_m = tot_len = 0
for n, w in tracks.items():
    Lt = 2*(w+PITCH); turns = NI/I_op; length = turns*Lt; m = length*T_TAPE*W_TAPE*D_TAPE
    tot_m += m; tot_len += length
    print(f" {n}: {turns:,.0f} tape-turns at {I_op:.0f} A, {length:,.0f} m of tape per m of track, {m:.1f} kg/m")
print(f"   both {tot_m:.0f} kg/m, {tot_m*L_RING/1e9:.2f} Mt, {tot_len*L_RING/1e9:,.0f} million km of tape")
print(f"   casing support tracks for scale: 2 tracks x 40 coils x 52 turns x 2 legs = {2*40*52*2:,} m of tape per m of ring")
print(" pulsed loss in the lift track, one pulse 0 -> B -> 0, strip in perpendicular field, fully penetrated:")
print("   Q = 1/2 x B_tape x w_tape x Ic per meter of tape (half the AC-cycle value B w Ic)")
Lt = 2*(0.76+PITCH); length = NI/I_op*Lt
cu50 = 20.9e3; cu10 = cu50/5
for Bt in (1.2, 0.3, 0.1):
    q = 0.5*Bt*W_TAPE*IC_1T*length
    print(f"   tape sees {Bt:.1f} T: {q/1e3:.2f} kJ per m of track at 70 K = {q/COP/1e3:,.0f} kJ electrical; "
          f"copper 50 A/mm2 {cu50/1e3:.1f} kJ, 10 A/mm2 {cu10/1e3:.1f} kJ; ratio {q/COP/cu50:.1f} and {q/COP/cu10:.1f}")
KE = 0.5*(5e6+525*10000)*30000**2
print(f"   sled kinetic energy {KE/1e9:,.0f} GJ; 5 full pulses per meter per launch (upper bound):")
for Bt in (1.2, 0.3, 0.1):
    q = 0.5*Bt*W_TAPE*IC_1T*length*5*L_RING
    print(f"   tape sees {Bt:.1f} T: {q/1e9:,.0f} GJ at 70 K, {q/COP/1e9:,.0f} GJ electrical, {q/COP/KE*100:.2f}% of kinetic energy")
print("\n-- passive Halbach limiter, 100 mm wavelength --")
k2 = 2*2*math.pi/0.1
F20_now = 2*26.4e3         # two top limiters, 1 m wide each, levitation.txt
for dmin in (0.020, 0.030, 0.040):
    Fd = k2*W*(GAP-dmin)                 # energy balance: F(dmin)/2k = W (gap - dmin)
    F20 = Fd*math.exp(k2*(dmin-0.020))
    print(f" turn the casing at {dmin*1e3:.0f} mm: F there {Fd/1e3:,.0f} kN/m, at 20 mm {F20/1e3:,.0f} kN/m, {F20/F20_now:,.0f} x the present pair")
print(f" speed after an 80 mm fall {math.sqrt(2*G_NET*0.08):.2f} m/s; 15 degree segments against 90: "
      f"field x {(math.sin(math.pi/24)/(math.pi/24))/(math.sin(math.pi/4)/(math.pi/4)):.2f}")
print("\n-- active: attractive pulsed track on the opposite face --")
F_load = 353e3
a = (F_load-W)/M_CASING
for trig in (0.005, 0.010, 0.020):
    t = math.sqrt(2*trig/G_NET); v = G_NET*t
    v2 = v+G_NET*0.002; d = trig+v*0.002+v2*v2/(2*a)
    print(f" trigger at {trig*1e3:.0f} mm ({t*1e3:.0f} ms after loss, {v:.2f} m/s), 2 ms rise, pull {F_load/1e3:.0f} kN/m "
          f"({F_load/W:.1f} x weight): stops after {d*1e3:.1f} mm, gap never below {100-d*1e3:.0f} mm")
