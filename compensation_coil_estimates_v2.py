"""Hand estimates, second pass (2026-10-11), after Paul's corrections. Not a simulation.
Changes from the first pass: the tape current is taken at the field the winding sees in an iron-core coil
(low, the flux is in the iron), copper is also rated for steady duty, 3 mm tape is added, and the
staggered 200 m coil layout is checked against the section-out failure case.
Tape data: Vol_III/bin THEVA AP 70 K field dependence (A per cm width x 1.2 for 12 mm)."""
import math
MU0 = 4e-7*math.pi
L_RING = 44.473e6
NI = 1.2*2*0.100/MU0
LT = 2*(0.76+0.5)          # m, lift track turn length, layout of the loss run [UNVERIFIED]
COP = 0.0106
print("Measured Ic of 12 mm AP tape at 70 K by field angle (A):")
print("  0.3 T: 487, 470, 529, 552 at 0, 30, 90, 120 deg;  1 T: 316, 300, 318, 323;  0.1 T at 0 deg: 684;  0 T: 892")
print("  -> at 70 K the angle changes Ic by under 10% in this data; the field strength at the winding is what matters.")
print(f"\nLift compensation track, NI = {NI:,.0f} A-turns per pole pair, turn length {LT:.2f} m")
print(" copper:")
cu = {}
for j in (50e6, 10e6, 5e6):
    m = NI/j*LT*8960; cu[j] = m
    print(f"   {j/1e6:4.0f} A/mm2: {m:6.0f} kg/m, heating at full field {1.68e-8*j*j/(8960*385):6.3f} K/s")
print(" HTS tape 0.08 mm thick, 8,890 kg/m3, run at half of Ic:")
for Bw, Ic12 in ((1.0, 316.0), (0.3, 487.0), (0.1, 684.0)):
    I = Ic12/2; length = NI/I*LT; m = length*0.08e-3*0.012*8890
    print(f"   winding field {Bw:.1f} T, Ic {Ic12:.0f} A, {I:.0f} A per tape ({I/0.96:.0f} A/mm2): {length:,.0f} m of tape per m, {m:5.1f} kg/m; "
          f"copper is {cu[50e6]/m:.0f}, {cu[10e6]/m:.0f}, {cu[5e6]/m:.0f} times heavier at 50, 10, 5 A/mm2")
print("\n Pulsed loss per meter of track at 30 km/s, upper estimate Q = 1/2 B_perp w Ic per meter of tape (fully penetrated strip):")
print("   it scales with tape width for the same total current, so 3 mm tape has a quarter of the 12 mm loss.")
for Bp, Ic12 in ((0.3, 487.0), (0.1, 684.0)):
    for w in (0.012, 0.003):
        Ic = Ic12*w/0.012; length = NI/(Ic/2)*LT
        q = 0.5*Bp*w*Ic*length
        print(f"   B_perp {Bp:.1f} T, {w*1e3:.0f} mm tape: {q:6.0f} J at 70 K = {q/COP/1e3:5.1f} kJ electrical "
              f"(copper: 20.9 kJ at 50 A/mm2, 4.2 at 10, 2.1 at 5)")
print("\nStaggered coils: 20 coils per track, each 200 m long and 5 cm wide, a new coil starting every 10 m.")
print(f"   one coil failed removes 1/20 of the pull over 200 m = {200/20:.0f} m of full pull; the run survives 50 m missing and fails at 200 m.")
print(f"   supplies per track: one per 10 m = {L_RING/10/1e6:.2f} million (50 m sections: {L_RING/50/1e6:.2f} million)")
print(f"   force footprint grows by 150 m on 1,770 m energized at 30 km/s = {150/1770*100:.0f}% more field energy; "
      f"coil transit at 30 km/s {200/30000*1e3:.1f} ms against a 50 ms pull")
print(f"   20 x 5 cm = 1.0 m of track against 0.76 m needed at 1.2 T: field needed {1.2*math.sqrt(0.76/1.0):.2f} T")
