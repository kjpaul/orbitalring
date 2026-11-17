#!/usr/bin/env python3
"""
Round 14: a truss spreader for the suspended solar array, in place of the ground spreader arms (Chapters 11 and 12 of Volume III).

Paul's design (17 November 2026): two truss legs joined by a hinge under the center line of the ring, a center line from the ring to the
hinge, a suspension line from the ring to each leg at two thirds of the way from the hinge to the tip, the retention line of each edge
fastened to a tip.  The frame is a strut in compression between the two edges.  Nothing goes to the ground.

What this script does:
  A. the statics of the open frame (flat, V or inverted V; what each line carries; which lines can go slack);
  B. the sizing of the frame as a three-chord truss of carbon-fibre tubes with CNT diagonals (and a stayed spine as a second row);
  C. Paul's question 2: panel mass x frame spacing;  question 1: how wide;  question 3: the 1,200 kg/m budget;
  D. opening the strip (quasi-static), and what holds the frame in roll and in yaw, in the simplest models;
  E. the ring without ground arms: the swing along the ring and the conductors.

The sheet, the retention line and the budget are round 9's (tools/round9_lines_array.py: v2, catenary, the constants), loaded from its
source by name as round 13 did, without running round 9's study.  Checked first against round 9's log.

Run from Python/orbitalring:  python tools/round14_truss_spreader.py
Log: runs/round14_logs/round14_truss_spreader.txt and .json
"""
import ast, io, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root
import ring_altitude as R, power_config as cfg, power_physics as phys

LOG_DIR = os.path.join(ROOT, "runs", "round14_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n")
OUT = {}

# ---------------------------------------------------------------------------------------------------------------------
# Round 9's sheet and budget, loaded by name
# ---------------------------------------------------------------------------------------------------------------------
R9 = os.path.join(HERE, "round9_lines_array.py")
NAMES = ["mu", "RE", "W_E", "RHO", "S_OP", "S_BR", "S75", "E_CNT", "h_r", "h_a", "L", "La", "r_top", "g_a", "g_r", "dP_arr", "kc", "S_SP"]
FUNCS = ["geff", "pot", "dia", "catenary", "parabola", "v2"]
def load_round9():
    tree = ast.parse(open(R9, encoding="utf-8").read())
    ns = dict(np=np, math=math, brentq=brentq, root=root, R=R, cfg=cfg, phys=phys, P=lambda *a: None)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in FUNCS:
            exec(compile(ast.Module(body=[node], type_ignores=[]), R9, "exec"), ns)
        elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id in NAMES for t in node.targets):
            exec(compile(ast.Module(body=[node], type_ignores=[]), R9, "exec"), ns)
    miss = [n for n in NAMES + FUNCS if n not in ns]
    if miss: raise RuntimeError("not found in round9_lines_array.py: %s" % miss)
    return ns
r9 = load_round9()
RHO = r9["RHO"]; S_OP = r9["S_OP"]; S75 = r9["S75"]; E_CNT = r9["E_CNT"]; g_a = r9["g_a"]; g_r = r9["g_r"]; L_RING = r9["L"]; La = r9["La"]
kc = r9["kc"]; h_r = r9["h_r"]; h_a = r9["h_a"]; v2 = r9["v2"]; dia = r9["dia"]; catenary = r9["catenary"]
LOG9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))
P9 = LOG9["power"]["catenary"]["rows"]                       # the book's rows of Table 12.1 (catenary sheet)
HS45, DR45 = catenary(45.0)                                  # half-span / half-sheet, drop / half-sheet
TW_PER_KM = P9["0.5"]["P_TW"]/P9["0.5"]["sheet_km"]          # TW for each km of sheet width
L_LINE = h_r - h_a                                           # the lines from the ring to the array, m

P("Round 14: a truss spreader for the suspended solar array.")
P("ring %.0f km, array %.0f km ; g at the array %.4f, at the ring %.4f m/s2 ; CNT %.0f kg/m3, %.2f GPa (cable, arms), %.2f GPa (steady lines), %.0f GPa modulus"
  % (h_r/1e3, h_a/1e3, g_a, g_r, RHO, S_OP/1e9, S75/1e9, E_CNT/1e9))
P("sheet (catenary, 45 deg at the ridge, flat at the edge): footprint %.4f of the sheet, edge %.4f of the sheet below the ridge ; line fraction kc %.4f ; %.3f TW per km of sheet"
  % (HS45, DR45/2, kc, TW_PER_KM))
d5 = v2(0.5)
ok0 = abs(d5["W"]/1e3/P9["0.5"]["sheet_km"] - 1) < 1e-9
P("CHECK 0. v2(0.5) here: sheet %.4f m, edge pull %.3f kN/m, sideways %.3f MN ; round 9's log: sheet %.4f m  -> %s"
  % (d5["W"], d5["H"]/1e3, d5["FR"]/1e6, P9["0.5"]["sheet_km"]*1e3, "reproduced" if ok0 else "NOT reproduced"))
OUT["check0"] = dict(sheet_m=d5["W"], H_kN_m=d5["H"]/1e3, FR_MN=d5["FR"]/1e6, reproduced=bool(ok0))

# ---------------------------------------------------------------------------------------------------------------------
# The compression material and the sizing rules (all stated; the assumed ones are marked)
# ---------------------------------------------------------------------------------------------------------------------
KSI = 6.894757e6
E_UD = 21.7e3*KSI          # Pa  0 deg compression modulus, IM7 UD, 77 F dry: 21.7 msi  (Hexcel, HexPly 8552 Product Data Sheet, US, copyright 2020, p. 3)
SC_UD = 245.0*KSI          # Pa  0 deg compression strength, IM7 UD, 77 F dry: 245 ksi  (same table)
RHO_RESIN = 0.0470*27679.9 # kg/m3  neat resin density 0.0470 lb/in3 (same data sheet, p. 1)
RHO_FIBRE = 1770.0         # kg/m3  IM7 fibre [from a search snippet of Hexcel's EU data sheet ; not opened, UNVERIFIED]
VF = 0.60                  # fibre volume [ASSUMED]
RHO_CF = VF*RHO_FIBRE + (1 - VF)*RHO_RESIN
AXIAL = 0.75               # share of the plies along the member ; modulus and strength taken as this share of the UD values [ASSUMED]
E_CF = AXIAL*E_UD; SC_CF = AXIAL*SC_UD
FOS = 2.0                  # on every buckling load and on the compressive strength [CHOSEN ; the book's yarn works at half its strength]
T_MIN = 1.0e-3             # m, smallest tube wall (8 plies of 0.125 mm) [CHOSEN]
D_MIN_CNT = 2.0e-3         # m, smallest CNT diagonal [CHOSEN]
JOINTS = 0.20              # joints, fittings, the hinge and its latch, as a share of the tube mass [ASSUMED]
BRACE = 0.02               # a brace is sized for 2 % of the force in the member it holds straight [ASSUMED, the common rule of thumb]
NU = 0.3
def gamma_shell(r_over_t):
    """Knockdown on the classical buckling stress of a thin cylinder in axial compression, NASA SP-8007 (1968):
    1 - 0.901 (1 - exp(-sqrt(r/t)/16)).  [Quoted from memory ; the document was not opened: UNVERIFIED.  It is for isotropic
    cylinders ; using it on a laminate is an approximation.]"""
    return 1 - 0.901*(1 - np.exp(-np.sqrt(r_over_t)/16.0))
P("\nCompression material: carbon fibre and epoxy, IM7 / 8552 unidirectional, room temperature, dry (Hexcel, HexPly 8552 Product Data Sheet, 2020):")
P("  0 deg compression modulus 21.7 msi = %.1f GPa ; 0 deg compression strength 245 ksi = %.0f MPa ; resin %.0f kg/m3" % (E_UD/1e9, SC_UD/1e6, RHO_RESIN))
P("  tube laminate: %.0f %% of the plies along the member [ASSUMED]: modulus %.1f GPa, strength %.0f MPa ; density %.0f kg/m3 at %.0f %% fibre [ASSUMED], fibre %.0f kg/m3 [UNVERIFIED]"
  % (AXIAL*100, E_CF/1e9, SC_CF/1e6, RHO_CF, VF*100, RHO_FIBRE))
P("  factor of safety %.1f on buckling and on strength ; smallest wall %.1f mm ; smallest CNT diagonal %.0f mm ; joints and hinge %.0f %% of the tubes ; braces sized for %.0f %% of the member force"
  % (FOS, T_MIN*1e3, D_MIN_CNT*1e3, JOINTS*100, BRACE*100))
OUT["material"] = dict(E_UD_GPa=E_UD/1e9, SC_UD_MPa=SC_UD/1e6, E_GPa=E_CF/1e9, SC_MPa=SC_CF/1e6, rho=RHO_CF, FOS=FOS, t_min_mm=T_MIN*1e3, joints=JOINTS, brace=BRACE, axial_share=AXIAL, Vf=VF)

_T = np.geomspace(T_MIN, 0.08, 90)
def tube(N, Lb, Amin=0.0, E=None, SC=None):
    """Lightest round tube for a compression N (N, already factored) over a pin-ended length Lb (m), with at least the area Amin.
    Checks: Euler N <= pi^3 E r^3 t / Lb^2 ; strength N <= 2 pi r t SC ; wall buckling N/(2 pi r t) <= gamma(r/t) E t / (r sqrt(3(1-nu^2))) ; t >= T_MIN.
    Returns r, t, area."""
    E = E_CF if E is None else E; SC = SC_CF if SC is None else SC
    t = _T
    r = np.maximum.reduce([(N*Lb*Lb/(math.pi**3*E*t))**(1/3.0), N/(2*math.pi*t*SC), Amin/(2*math.pi*t), 5*t])
    crit = gamma_shell(r/t)*E*t/(r*math.sqrt(3*(1 - NU*NU)))
    ok = N/(2*math.pi*r*t) <= crit
    A = np.where(ok, 2*math.pi*r*t, np.inf); i = int(np.argmin(A))
    if not np.isfinite(A[i]): return float("nan"), float("nan"), float("inf")
    return float(r[i]), float(t[i]), float(A[i])

def truss(Pc, Lf, E=None, SC=None, rho=None, detail=False, nsup=1):
    """Lightest three-chord truss strut (equilateral section of side d, bays of length lb) for a compression Pc (N) between pinned ends Lf apart.
    Chords and battens are tubes of the compression material ; each face has crossed CNT diagonals that only pull.
    Loads: Pc ; the frame's own weight at the array, carried by the center line and by the two suspension lines at two thirds of each leg,
    with the bending that leaves (w l^2 / 18 over a suspension line, l the leg), magnified by 1 / (1 - 1/FOS).
    nsup: suspension lines on each leg ; 1 is Paul's layout ; with more, evenly spread, the bending is taken as w l^2 / (18 nsup^2) [the same form, shorter spans].
    Checks: Euler buckling of the whole strut over Lf with the shear give of the diagonals ; Euler buckling of a chord between nodes ; chord strength ;
    chord wall buckling ; the battens as struts for the brace force ; the smallest wall."""
    E = E_CF if E is None else E; SC = SC_CF if SC is None else SC; rho = RHO_CF if rho is None else rho
    leg = Lf/2; best = None
    for d in np.geomspace(max(Lf/400, 0.5), Lf/4, 34):
        for lb in d*np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0]):
            nb = max(2, int(round(Lf/lb))); lb_ = Lf/nb
            ld = math.hypot(d, lb_); sph = d/ld; cph = lb_/ld
            m = 0.0; res = None
            for _ in range(6):
                w = m*g_a/Lf; M = w*leg*leg/18.0/nsup**2/(1 - 1/FOS); Vw = 0.417*w*leg/nsup/(1 - 1/FOS)
                Q = BRACE*Pc + Vw
                Ad = max(Q/(1.5*sph)/S_OP, math.pi*D_MIN_CNT**2/4, 5*FOS*Pc/(1.5*E_CNT*sph*sph*cph))   # strength, smallest size, shear stiffness 5 x FOS x Pc
                GA = 1.5*E_CNT*Ad*sph*sph*cph
                PE = 1.0/(1.0/(FOS*Pc) - 1.0/GA)                                                        # Euler load the chords have to give
                Aeul = 2*PE*Lf*Lf/(math.pi**2*E*d*d)                                                    # I = A d^2 / 2 for three chords
                Nc = FOS*(Pc/3 + M/(0.866*d))
                rc, tc, Ac = tube(Nc, lb_, Aeul, E, SC)
                Nb = BRACE*Nc + FOS*Q/(1.5*sph)*sph
                rb, tb, Ab = tube(Nb, d, 0.0, E, SC)
                m_ch = 3*Ac*Lf*rho; m_bat = 3*Ab*d*rho*(nb + 1); m_di = 6*Ad*ld*RHO*nb
                m_new = (m_ch + m_bat)*(1 + JOINTS) + m_di
                res = dict(d=d, lb=lb_, bays=nb, r_chord=rc, t_chord=tc, A_chord=Ac, r_batten=rb, t_batten=tb, d_diag_mm=2e3*math.sqrt(Ad/math.pi),
                           m_chords=m_ch, m_battens=m_bat, m_diag=m_di, m_joints=(m_ch + m_bat)*JOINTS, mass=m_new,
                           stress_chord_MPa=Nc/FOS/Ac/1e6, euler_governs=bool(Aeul >= Ac*0.999), M_weight_MNm=M/1e6, N_chord_MN=Nc/FOS/1e6)
                if not np.isfinite(m_new) or abs(m_new - m) < 1e-4*m_new: m = m_new; break
                m = m_new
            if res and np.isfinite(m) and (best is None or m < best["mass"]): best = res
    if best is None: best = dict(mass=float('inf'), d=float('nan'))
    return best if detail else best["mass"]

def stayed(Pc, Lf):
    """A second row, rougher: one spine tube held straight by three CNT stays that run from tip to tip over cross-arms of length a.
    The stays give the bending stiffness, E_CNT x A_stay x 1.5 a^2 x 0.8 (the 0.8 for the stays' taper to the tips [ASSUMED]).
    They have to stay taut, so each is pre-tensioned to 1.25 times the strain the spine has under its working load [ASSUMED rule],
    and that pre-tension adds to the spine's compression.  Cross-arms are tubes sized for the brace force.  Factor of safety as for the truss."""
    best = None
    for a in np.geomspace(Lf/200, Lf/6, 30):
        for n in (3, 5, 7, 11, 15, 23, 31, 47, 63):
            bay = Lf/(n + 1)
            As = FOS*Pc*Lf*Lf/(math.pi**2*E_CNT*1.5*a*a*0.8)
            for over in (1.0, 1.5, 2.0, 3.0, 5.0):
                T0 = 0.0
                for _ in range(40):
                    Nsp = FOS*(Pc + 3*T0); rs, ts, Asp = tube(Nsp, bay, over*Nsp/SC_CF)
                    if not np.isfinite(Asp): break
                    T0n = 1.25*E_CNT*As*(Pc + 3*T0)/(E_CF*Asp)
                    if abs(T0n - T0) < 1e-6*Pc: T0 = T0n; break
                    T0 = 0.5*(T0 + T0n)
                if not np.isfinite(Asp) or T0 > 20*Pc: continue
                ra, ta, Aa = tube(BRACE*Nsp, a)
                m = (Asp*Lf + 3*n*Aa*a)*RHO_CF*(1 + JOINTS) + 3*As*Lf*RHO*1.02
                if best is None or m < best["mass"]:
                    best = dict(mass=m, a=a, stations=n, r_spine=rs, t_spine=ts, A_spine=Asp, d_stay_mm=2e3*math.sqrt(As/math.pi), T0_MN=T0/1e6,
                                spine_load_MN=(Pc + 3*T0)/1e6, m_spine=Asp*Lf*RHO_CF, m_arms=3*n*Aa*a*RHO_CF, m_stays=3*As*Lf*RHO*1.02)
    return best

def frame_system(sA, W, S, **kw):
    """One frame for a sheet W wide with panels sA (kg/m2) and frames every S (m): the push it carries, its mass, its lines, the retention line."""
    q = sA*g_a; H = q*W/2; F = H*S; Lf = 2*HS45*(W/2)
    t = truss(F, Lf, detail=True, **kw); m = t["mass"]
    Wt = m*g_a
    Ts = 0.375*Wt; Tc = 0.25*Wt                               # flat frame with the hinge free: each suspension line 3/8, the center line 1/4
    m_lines = kc*m                                            # one thickness, 18.75 GPa at the top, as the array support lines (Equation 12.1)
    stay = 4*BRACE*F/S75*RHO*math.hypot(S, Lf/2)*1.0          # four yaw stays, tip to the next frames' hinges, for 2 % of the push [ASSUMED]
    TR = q*S*S/2; mR = 2*RHO*TR/S75                           # retention line with a bow of W/8, as v2
    casing = (m*(1 + kc) + stay)/S
    return dict(F=F, Lf=Lf, H=H, truss=t, mass_t=m/1e3, Ts_MN=Ts/1e6, Tc_MN=Tc/1e6, d_s_mm=dia(Ts/(1 - kc/(1 + kc)), S75), d_c_mm=dia(Tc/(1 - kc/(1 + kc)), S75),
                lines_t=m_lines/1e3, stays_t=stay/1e3, casing_kgm=casing, TR_MN=TR/1e6, ret_mm=dia(TR, S75), ret_kgm=mR, bow_m=W/8,
                lean_deg=math.degrees(math.atan((2/3.0)*(Lf/2)/L_LINE)))

# ---------------------------------------------------------------------------------------------------------------------
# PART A.  The statics of the open frame
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nPART A.  The statics of the open frame (one leg ; l its length ; F the push at the tip ; Wl the weight of the leg)\n" + "="*110)
P("  The lines from the ring are 100 km long and the frame is under a kilometer, so they are within a fraction of a degree of the vertical: they can only lift.")
P("  Moments about the hinge, with the hinge a height e below the line that joins the two tips (e > 0 a shallow V, e < 0 an inverted V), suspension line at 2l/3:")
P("      suspension line  Ts = 0.75 Wl - 1.5 F e / (l cos)      center line (its half for this leg)  Tc/2 = 0.25 Wl + 1.5 F e / (l cos)")
P("  Both lines pull only while  - Wl / 6  <  F e / l  <  Wl / 2 :  the hinge has to sit within a small distance of the line of thrust.")
fs = frame_system(0.5, d5["W"], 10e3); Wl = fs["mass_t"]*1e3*g_a/2; leg = fs["Lf"]/2
P("  Design case (0.5 kg/m2, sheet %.1f m, frames every 10 km): push %.2f MN, leg %.1f m, frame %.0f t, leg weight %.2f MN."
  % (d5["W"], fs["F"]/1e6, leg, fs["mass_t"], Wl/1e6))
e_up = Wl/6*leg/fs["F"]; e_dn = Wl/2*leg/fs["F"]
P("     the hinge may sit from %.2f m above the tips' line (center line slack) to %.2f m below it (suspension lines slack). Flat: suspension lines %.2f MN each, center line %.2f MN."
  % (e_up, e_dn, 0.75*Wl/1e6, 0.5*Wl/1e6))
P("     For comparison the sheet's ridge is %.0f m above the tips' line.  An inverted V with its hinge at the ridge is not possible with lines that only lift." % (DR45*d5["W"]/2))
k_line = E_CNT*(0.5*Wl/(1 - kc/(1 + kc))/S75)/L_LINE; k_need = fs["F"]/leg
P("  How stiff the lines are.  The center line, sized for its %.2f MN, gives %.1f N for each meter the hinge moves (E A / length, 100 km)." % (0.5*Wl/1e6, k_line))
P("     A hinged strut under %.2f MN needs about F / l = %.0f kN/m at the hinge to keep it from folding.  The line gives %.1e of that." % (fs["F"]/1e6, k_need/1e3, k_line/k_need))
A_need = 2*k_need*L_LINE/E_CNT
P("     A center line stiff enough (2 F / l) would be %.0f mm across and %.0f t, %.0f kg for each meter of ring at 10 km.  So the hinge has to latch when the frame is open,"
  % (2e3*math.sqrt(A_need/math.pi), RHO*A_need*L_LINE/1e3, RHO*A_need*L_LINE/10e3))
P("     and the open frame is sized as one strut from tip to tip, with no help from the lines against buckling.  The lines carry its weight and nothing else.")
OUT["statics"] = dict(push_MN=fs["F"]/1e6, leg_m=leg, frame_t=fs["mass_t"], hinge_above_max_m=e_up, hinge_below_max_m=e_dn, Ts_flat_MN=0.75*Wl/1e6, Tc_flat_MN=0.5*Wl/1e6,
                      k_center_line_N_m=k_line, k_needed_N_m=k_need, stiff_line_mm=2e3*math.sqrt(A_need/math.pi), stiff_line_t=RHO*A_need*L_LINE/1e3)

# ---------------------------------------------------------------------------------------------------------------------
# PART B.  The frame of the design case, piece by piece
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nPART B.  The frame of the design case (0.5 kg/m2, sheet %.1f m, footprint %.1f m, frames every 10 km, push %.2f MN)\n" % (d5["W"], fs["Lf"], fs["F"]/1e6) + "="*110)
def show_truss(t, nm):
    P("  %s: depth (side of the triangle) %.1f m, %d bays of %.1f m ; chords: tubes %.0f mm across, wall %.1f mm, %.0f MPa at work (%s governs) ; battens: tubes %.0f mm across, wall %.1f mm ; diagonals %.1f mm CNT"
      % (nm, t["d"], t["bays"], t["lb"], 2e3*t["r_chord"], 1e3*t["t_chord"], t["stress_chord_MPa"], "Euler buckling of the whole strut" if t["euler_governs"] else "the chord between nodes",
         2e3*t["r_batten"], 1e3*t["t_batten"], t["d_diag_mm"]))
    P("     mass: chords %.1f t + battens %.1f t + joints and hinge %.1f t + diagonals %.2f t = %.1f t ; bending from its own weight %.1f MN m ; chord force %.2f MN"
      % (t["m_chords"]/1e3, t["m_battens"]/1e3, t["m_joints"]/1e3, t["m_diag"]/1e3, t["mass"]/1e3, t["M_weight_MNm"], t["N_chord_MN"]))
show_truss(fs["truss"], "truss, IM7 / 8552")
t2 = truss(fs["F"], fs["Lf"], E=2*E_CF, detail=True); show_truss(t2, "WHAT IF the modulus were twice (%.0f GPa), same strength and density" % (2*E_CF/1e9))
t3 = truss(fs["F"], fs["Lf"], E=E_CNT, SC=S_OP*2, rho=RHO, detail=True); show_truss(t3, "WHAT IF [SPECULATIVE] the CNT yarn could be made into tubes that take compression (500 GPa, 25 GPa, 1,700 kg/m3)")
st = stayed(fs["F"], fs["Lf"])
P("  stayed spine (rough): cross-arms %.1f m at %d stations ; spine tube %.0f mm across, wall %.1f mm, carrying %.1f MN (the push and the stays' pre-tension of %.2f MN each) ; stays %.1f mm"
  % (st["a"], st["stations"], 2e3*st["r_spine"], 1e3*st["t_spine"], st["spine_load_MN"], st["T0_MN"], st["d_stay_mm"]))
P("     mass: spine %.1f t + cross-arms %.1f t + joints + stays %.1f t = %.1f t" % (st["m_spine"]/1e3, st["m_arms"]/1e3, st["m_stays"]/1e3, st["mass"]/1e3))
P("  Lines of the truss frame: suspension lines %.2f MN each, %.1f mm ; center line %.2f MN, %.1f mm ; all three %.1f t (%.1f %% of the frame, as Equation 12.1) ; they lean %.2f deg from the vertical"
  % (fs["Ts_MN"], fs["d_s_mm"], fs["Tc_MN"], fs["d_c_mm"], fs["lines_t"], kc*100, fs["lean_deg"]))
P("  On the casing: frame %.0f t + lines %.1f t + yaw stays %.2f t every 10 km = %.1f kg for each meter of ring.   The ground arms: 572 kg/m (book, straight), 1,195 kg/m (round 13, hanging)."
  % (fs["mass_t"], fs["lines_t"], fs["stays_t"], fs["casing_kgm"]))
OUT["design_frame"] = dict(truss=fs["truss"], what_if_modulus_x2=t2, what_if_cnt_tubes=t3, stayed=st, Ts_MN=fs["Ts_MN"], Tc_MN=fs["Tc_MN"], d_s_mm=fs["d_s_mm"], d_c_mm=fs["d_c_mm"],
                           lines_t=fs["lines_t"], stays_t=fs["stays_t"], casing_kgm=fs["casing_kgm"], lean_deg=fs["lean_deg"], push_MN=fs["F"]/1e6, Lf=fs["Lf"])

# ---------------------------------------------------------------------------------------------------------------------
# PART C.  Question 2: panel mass x spacing, at the book's widths ; then the budget ; then how wide
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nPART C2.  Panel mass and frame spacing, each at the book's width for that panel (Table 12.1)\n" + "="*110)
SP = (0.5e3, 1e3, 2e3, 5e3, 10e3)
OUT["spacing"] = {}
for sA in (1.0, 0.5, 0.2):
    W = P9[str(sA)]["sheet_km"]*1e3
    P("  %.1f kg/m2, sheet %.3f km (footprint %.3f km, bow of the retention line %.0f m), panels and support lines %.1f kg/m:" % (sA, W/1e3, 2*HS45*W/2e3, W/8, sA*W*La/L_RING*(1 + kc)))
    P("    spacing | push     | depth   | frame    | frame + lines on the casing | retention line          | retention lines | suspension, center line | frames + retention")
    OUT["spacing"][str(sA)] = {}
    for S in SP:
        f = frame_system(sA, W, S); tot = f["casing_kgm"] + f["ret_kgm"]
        P("    %4.1f km | %6.2f MN | %5.1f m | %6.1f t | %7.1f kg/m                | %7.1f MN, %5.1f mm      | %6.2f kg/m     | %.2f, %.2f MN          | %7.1f kg/m"
          % (S/1e3, f["F"]/1e6, f["truss"]["d"], f["mass_t"], f["casing_kgm"], f["TR_MN"], f["ret_mm"], f["ret_kgm"], f["Ts_MN"], f["Tc_MN"], tot))
        OUT["spacing"][str(sA)]["%g" % (S/1e3)] = dict(push_MN=f["F"]/1e6, depth_m=f["truss"]["d"], frame_t=f["mass_t"], casing_kgm=f["casing_kgm"], TR_MN=f["TR_MN"], ret_mm=f["ret_mm"],
                                                        ret_kgm=f["ret_kgm"], Ts_MN=f["Ts_MN"], Tc_MN=f["Tc_MN"], total_kgm=tot, lines_t=f["lines_t"],
                                                        chord_mm=2e3*f["truss"]["r_chord"], wall_mm=1e3*f["truss"]["t_chord"], bays=f["truss"]["bays"], bay_m=f["truss"]["lb"])
    Sg = np.geomspace(0.3e3, 30e3, 41); tt = [(lambda f: f["casing_kgm"] + f["ret_kgm"])(frame_system(sA, W, S)) for S in Sg]; i = int(np.argmin(tt))
    fb = frame_system(sA, W, Sg[i])
    P("    best spacing %.1f km: frames %.1f kg/m + retention lines %.1f kg/m = %.1f kg/m (frame %.0f t, push %.1f MN)   [within 10 %% of that from %.1f to %.1f km]"
      % (Sg[i]/1e3, fb["casing_kgm"], fb["ret_kgm"], tt[i], fb["mass_t"], fb["F"]/1e6, min(S for S, t in zip(Sg, tt) if t < 1.1*tt[i])/1e3, max(S for S, t in zip(Sg, tt) if t < 1.1*tt[i])/1e3))
    OUT["spacing"][str(sA)]["best"] = dict(S_km=Sg[i]/1e3, frames_kgm=fb["casing_kgm"], ret_kgm=fb["ret_kgm"], total_kgm=tt[i], frame_t=fb["mass_t"], push_MN=fb["F"]/1e6)

P("\n" + "="*110 + "\nPART C3.  The budget: 1,200 kg for each meter of ring held\n" + "="*110)
def budget(sA, S, M=1200.0, **kw):
    """The widest sheet for which panels and support lines + retention lines + frames, their lines and stays come to M."""
    def tot(W):
        f = frame_system(sA, W, S, **kw); return sA*W*La/L_RING*(1 + kc) + f["ret_kgm"] + f["casing_kgm"] - M
    if tot(50.0) > 0: return None
    W = brentq(tot, 50.0, M/(sA*La/L_RING*(1 + kc)), xtol=0.01); f = frame_system(sA, W, S, **kw)
    return dict(W=W, pan=sA*W*La/L_RING*(1 + kc), ret=f["ret_kgm"], frames=f["casing_kgm"], P_TW=TW_PER_KM*W/1e3, f=f)
OUT["budget"] = {}
BOOK = {"1.0": (0.54, 5.9), "0.5": (1.11, 12.2), "0.2": (2.83, 31.1)}
for sA in (1.0, 0.5, 0.2):
    P("  %.1f kg/m2:" % sA); rows = {}
    for S in (1e3, 2e3, 3e3, 5e3, 7e3, 10e3, 15e3, 25e3):
        b = budget(sA, S)
        if b is None: continue
        rows[S] = b
        P("    frames every %4.1f km: sheet %.3f km ; panels and support lines %.1f + retention %.1f + frames, lines, stays %.1f kg/m ; frame %.0f t, %.1f m deep, push %.1f MN ; power %.2f TW"
          % (S/1e3, b["W"]/1e3, b["pan"], b["ret"], b["frames"], b["f"]["mass_t"], b["f"]["truss"]["d"], b["f"]["F"]/1e6, b["P_TW"]))
    Sb = max(rows, key=lambda s: rows[s]["W"]); b = rows[Sb]
    # refine the best spacing
    Sg = np.geomspace(Sb/1.6, Sb*1.6, 9); bb = [budget(sA, s) for s in Sg]; j = int(np.argmax([x["W"] if x else 0 for x in bb])); b = bb[j]; Sb = Sg[j]
    P("    BEST: frames every %.1f km: sheet %.3f km, footprint %.3f km, %.2f TW   (book, straight ground arms: %.2f km, %.1f TW)   panels %.1f + retention %.1f + frames %.1f kg/m"
      % (Sb/1e3, b["W"]/1e3, 2*HS45*b["W"]/2e3, b["P_TW"], BOOK[str(sA)][0], BOOK[str(sA)][1], b["pan"], b["ret"], b["frames"]))
    f = b["f"]
    P("          frame %.0f t, %.1f m deep, %.0f m long, push %.1f MN ; retention line %.0f MN, %.0f mm, bow %.0f m ; suspension lines %.2f MN (%.1f mm), center line %.2f MN (%.1f mm) ; %.2f GW for each frame"
      % (f["mass_t"], f["truss"]["d"], f["Lf"], f["F"]/1e6, f["TR_MN"], f["ret_mm"], f["bow_m"], f["Ts_MN"], f["d_s_mm"], f["Tc_MN"], f["d_c_mm"], b["P_TW"]*1e3/(L_RING/Sb)))
    OUT["budget"][str(sA)] = dict(S_km=Sb/1e3, sheet_km=b["W"]/1e3, footprint_km=2*HS45*b["W"]/2e3, drop_m=DR45*b["W"]/2, P_TW=b["P_TW"], pan=b["pan"], ret=b["ret"], frames=b["frames"],
                                  frame_t=f["mass_t"], depth_m=f["truss"]["d"], Lf=f["Lf"], push_MN=f["F"]/1e6, TR_MN=f["TR_MN"], ret_mm=f["ret_mm"], bow_m=f["bow_m"],
                                  Ts_MN=f["Ts_MN"], Tc_MN=f["Tc_MN"], d_s_mm=f["d_s_mm"], d_c_mm=f["d_c_mm"], GW_per_frame=b["P_TW"]*1e3/(L_RING/Sb), n_frames=L_RING/Sb,
                                  book_sheet_km=BOOK[str(sA)][0], book_TW=BOOK[str(sA)][1], kgm_per_TW=1200/b["P_TW"], edge_pull_kN_m=f["H"]/1e3,
                                  rows={"%g" % (s/1e3): dict(sheet_km=r["W"]/1e3, P_TW=r["P_TW"], pan=r["pan"], ret=r["ret"], frames=r["frames"], frame_t=r["f"]["mass_t"]) for s, r in rows.items()})
    for ns in (3,):
        Sg = (5e3, 7e3, 10e3, 15e3, 20e3, 25e3, 35e3); bb = [budget(sA, s, nsup=ns) for s in Sg]; j = int(np.argmax([x["W"] if x else 0 for x in bb])); b3 = bb[j]
        P("    WITH %d suspension lines on each leg in place of one: frames every %.0f km: sheet %.3f km, %.2f TW ; panels %.1f + retention %.1f + frames %.1f kg/m ; frame %.0f t, %.1f m deep"
          % (ns, Sg[j]/1e3, b3["W"]/1e3, b3["P_TW"], b3["pan"], b3["ret"], b3["frames"], b3["f"]["mass_t"], b3["f"]["truss"]["d"]))
        OUT["budget"][str(sA)]["nsup%d" % ns] = dict(S_km=Sg[j]/1e3, sheet_km=b3["W"]/1e3, P_TW=b3["P_TW"], pan=b3["pan"], ret=b3["ret"], frames=b3["frames"], frame_t=b3["f"]["mass_t"], depth_m=b3["f"]["truss"]["d"])
P("  Round 13's hanging ground arm at 0.5 kg/m2 with 1,200 kg/m held: 0.72 km, 8.0 TW.")

P("\n" + "="*110 + "\nPART C1.  How wide: one frame and its lines against the width of the sheet\n" + "="*110)
OUT["width"] = {}
for sA in (1.0, 0.5, 0.2, 0.05):
    P("  %.2f kg/m2%s:" % (sA, "  (not a row of Table 12.1 ; a much lighter sheet)" if sA < 0.2 else ""))
    P("    sheet   | frame length | best spacing | push     | frame     | depth   | frames + lines | retention | panels + support lines | all      | frames' share | kg/m for each TW")
    OUT["width"][str(sA)] = {}
    for W in (1.11e3, 2e3, 3e3, 5e3, 7e3, 10e3):
        Sg = np.geomspace(0.3e3, 40e3, 29); ff = [frame_system(sA, W, S) for S in Sg]; tt = [f["casing_kgm"] + f["ret_kgm"] for f in ff]; i = int(np.argmin(tt)); f = ff[i]
        pan = sA*W*La/L_RING*(1 + kc); al = pan + tt[i]
        if not np.isfinite(al): P("    %5.2f km | %7.0f m    | with one suspension line on each leg no frame was found in the search (tubes up to 80 mm wall)" % (W/1e3, f["Lf"]))
        else: P("    %5.2f km | %7.0f m    | %5.1f km     | %6.1f MN | %7.0f t | %5.1f m | %7.1f kg/m   | %6.1f    | %7.1f kg/m           | %7.1f  | %4.1f %%        | %.1f"
          % (W/1e3, f["Lf"], Sg[i]/1e3, f["F"]/1e6, f["mass_t"], f["truss"]["d"], f["casing_kgm"], f["ret_kgm"], pan, al, 100*f["casing_kgm"]/al, al/(TW_PER_KM*W/1e3)))
        f3s = [frame_system(sA, W, S_, nsup=3) for S_ in Sg]; t3 = [x["casing_kgm"] + x["ret_kgm"] for x in f3s]; i3 = int(np.argmin(t3))
        P("              with 3 suspension lines on each leg: best spacing %.1f km, frame %.0f t, frames + lines %.1f kg/m, retention %.1f, all %.1f kg/m, %.1f kg/m for each TW"
          % (Sg[i3]/1e3, f3s[i3]["mass_t"], f3s[i3]["casing_kgm"], f3s[i3]["ret_kgm"], pan + t3[i3], (pan + t3[i3])/(TW_PER_KM*W/1e3)))
        OUT["width"][str(sA)]["%g" % (W/1e3)] = dict(nsup3=dict(S_km=Sg[i3]/1e3, frame_t=f3s[i3]["mass_t"], frames_kgm=f3s[i3]["casing_kgm"], ret_kgm=f3s[i3]["ret_kgm"], all_kgm=pan + t3[i3], kgm_per_TW=(pan + t3[i3])/(TW_PER_KM*W/1e3)),
                                                      Lf=f["Lf"], S_km=Sg[i]/1e3, push_MN=f["F"]/1e6, frame_t=f["mass_t"], depth_m=f["truss"]["d"], frames_kgm=f["casing_kgm"], ret_kgm=f["ret_kgm"],
                                                      pan_kgm=pan, all_kgm=al, share=f["casing_kgm"]/al, kgm_per_TW=al/(TW_PER_KM*W/1e3), lines_t=f["lines_t"])
P("  Panels and support lines alone are %.1f kg/m for each TW at 1.0 kg/m2, %.1f at 0.5, %.1f at 0.2 (the floor that no spreader can go under)."
  % tuple(sA*1e3*La/L_RING*(1 + kc)/TW_PER_KM for sA in (1.0, 0.5, 0.2)))

# ---------------------------------------------------------------------------------------------------------------------
# PART D.  Opening the strip, and what holds the frame in roll and in yaw
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nPART D.  Opening the strip (quasi-static), design case\n" + "="*110)
def sheet_at(y, z, s, q):
    """A hanging sheet of arc s and weight q per meter of arc, from the ridge at (0, 0) to its edge at (y, z).  Returns the pull on whatever holds
    the edge: H inward (N per m of ring) and V upward (N per m of ring ; negative when the sheet hangs below the edge), and the height of its center of mass."""
    c = math.hypot(y, z)
    if c >= s*(1 - 1e-9): return None
    k = math.sqrt(s*s - z*z)/y
    a = y/(2*brentq(lambda u: math.sinh(u)/u - k, 1e-9, 60.0))
    x0 = y/2 - a*math.atanh(z/s)
    ue = (y - x0)/a; u0 = -x0/a
    H = q*a; V = -H*math.sinh(ue)
    zc = a*a/s*((ue - u0)/2 + (math.sinh(2*ue) - math.sinh(2*u0))/4) - a*math.cosh(u0)          # mean of z over the arc, z measured from the ridge
    return H, V, zc
W = d5["W"]; s_half = W/2; q = 0.5*g_a; b = HS45*s_half; drop = DR45*s_half; S = 10e3; leg = b
chk = sheet_at(b, -drop, s_half, q)
P("  check of the sheet solver at the open position: pull %.4f kN/m (v2: %.4f), vertical %.2e kN/m (flat edge: 0), center of mass %.1f m below the ridge" % (chk[0]/1e3, d5["H"]/1e3, chk[1]/1e3, -chk[2]))
Wl = fs["mass_t"]*1e3*g_a/2
E_sheet = q*s_half*S*(chk[2] - (-s_half/2)); E_frame = Wl*(2*leg/3)*0.5                 # leg's center of mass drops l/6... see the text below
P("  Work to draw one half of the sheet out from hanging straight down, for one 10 km span: weight %.2f MN x rise of its center of mass %.1f m = %.2f GJ."
  % (q*s_half*S/1e6, chk[2] + s_half/2, E_sheet/1e9))
P("  If the hinge is lowered, the only thing that drives the legs apart is the frame's own weight coming down: a leg of %.2f MN whose center of mass falls l/6 = %.0f m gives %.3f GJ,"  % (Wl/1e6, leg/6, Wl*leg/6/1e9))
P("  %.1f %% of the work.  The rest has to come from a line that is hauled in under load." % (100*Wl*leg/6/E_sheet))
P("\n  The same, step by step.  The suspension points stay at the height of the open tips ; phi is the leg's angle from the vertical (0 folded, 90 open).")
for path, sgn in (("hinge above (Paul's: legs hang down, the hinge is lowered %.0f m)" % (2*leg/3), -1), ("hinge below (legs point up, the hinge is raised %.0f m)" % (2*leg/3), +1)):
    P("  Path: %s" % path)
    P("    phi    | tip: out, below the ridge | sheet pull at the tip: inward, downward | suspension line | center line")
    rows = []
    for ph in (10, 20, 30, 40, 50, 60, 70, 80, 85, 88, 89, 89.5, 90):
        p_ = math.radians(ph); y = leg*math.sin(p_); z = -drop + sgn*(leg/3)*math.cos(p_)
        sh = sheet_at(y, z, s_half, q)
        if sh is None: continue
        Hh, Vv, _ = sh; F = Hh*S; Vup = Vv*S; cot = math.cos(p_)/math.sin(p_)
        if sgn < 0: Ts = 1.5*F*cot + 0.75*Wl - 1.5*Vup
        else:       Ts = 0.75*Wl - 1.5*Vup - 1.5*F*cot
        Tc = 2*(Wl - Vup - Ts)
        rows.append(dict(phi=ph, y=y, z=z, F_MN=F/1e6, Vdown_MN=-Vup/1e6, Ts_MN=Ts/1e6, Tc_MN=Tc/1e6))
        P("    %5.1f  | %6.1f m, %6.1f m        | %7.2f MN, %7.2f MN                   | %8.2f MN%s | %8.2f MN%s"
          % (ph, y, -z, F/1e6, -Vup/1e6, Ts/1e6, "  SLACK" if Ts < 0 else "       ", Tc/1e6, "  SLACK (it would have to push)" if Tc < 0 else ""))
    OUT.setdefault("opening", {})["above" if sgn < 0 else "below"] = rows
ok_rows = [r for r in OUT["opening"]["below"] if r["Ts_MN"] > 0 and r["Tc_MN"] > 0]
Ts_max = max(r["Ts_MN"] for r in OUT["opening"]["below"]); Tc_max = max(r["Tc_MN"] for r in OUT["opening"]["below"])
m_open = (2*Ts_max + Tc_max)*1e6/(1 - kc/(1 + kc))/S75*RHO*L_LINE
P("  Result: on Paul's path the center line would have to push down on the hinge, with up to %.0f MN, from the start until the legs are within about 1 deg of flat.  A line cannot do that."
  % max(-r["Tc_MN"] for r in OUT["opening"]["above"]))
P("  On the other path (the legs point up from the hinge and the hinge is raised) both lines pull from phi = %.0f deg to open: the tips carry part of the sheet's weight, and hauling in the center line does the work."
  % min(r["phi"] for r in ok_rows))
P("  The largest pulls on that path: suspension line %.1f MN (at the start, when half of each half-sheet hangs from the tips), center line %.1f MN.  Open, they carry %.2f and %.2f MN."
  % (Ts_max, Tc_max, fs["Ts_MN"], fs["Tc_MN"]))
P("  Lines sized for the opening and not only for the open frame: %.1f mm and %.1f mm, %.0f t for the three in place of %.0f t, which is %.1f kg for each meter of ring at 10 km in place of %.1f."
  % (dia(Ts_max*1e6/(1 - kc/(1 + kc)), S75), dia(Tc_max*1e6/(1 - kc/(1 + kc)), S75), m_open/1e3, fs["lines_t"], m_open/10e3, fs["lines_t"]*1e3/10e3))
P("  The other way is to open the frame first with the tethers slack (then only its weight is on the lines) and draw the sheet out with the tethers afterwards, as Section 12.11 does with the retention lines.")
OUT["opening"]["below_Ts_max_MN"] = Ts_max; OUT["opening"]["below_Tc_max_MN"] = Tc_max; OUT["opening"]["lines_for_opening_t"] = m_open/1e3
OUT["opening"]["above_push_needed_MN"] = max(-r["Tc_MN"] for r in OUT["opening"]["above"])
OUT["opening"]["work_sheet_GJ"] = E_sheet/1e9; OUT["opening"]["work_frame_GJ"] = Wl*leg/6/1e9; OUT["opening"]["hinge_travel_m"] = 2*leg/3
OUT["opening"]["fold"] = dict(leg_m=leg, drop_m=drop, hinge_above_ridge_folded_m=2*leg/3 - drop, tip_below_ridge_folded_m=drop + leg/3, sheet_half_m=s_half)
P("  The fold: leg %.0f m ; the hinge travels 2l/3 = %.0f m against the suspension lines.  Open, the frame lies %.0f m below the ridge, at the height of the edges, under the sheet."
  % (leg, 2*leg/3, drop))
P("  Folded Paul's way the hinge is %.0f m above the ridge and the tips %.0f m below it (the folded sheet hangs %.0f m).  The center line is the array support line of that place, carried on down."
  % (2*leg/3 - drop, drop + leg/3, s_half))

P("\n" + "="*110 + "\nPART D2.  Roll and yaw of the open, latched frame, simplest models (design case)\n" + "="*110)
dz = 0.5
Hp, Vp, _ = sheet_at(b, -drop + dz, s_half, q); Hm, Vm, _ = sheet_at(b, -drop - dz, s_half, q)
kv = -(Vp - Vm)/(2*dz)                                    # N per m of ring per m: the sheet pulls a raised edge down
dy = 0.5
Hp2, _, _ = sheet_at(b + dy, -drop, s_half, q); Hm2, _, _ = sheet_at(b - dy, -drop, s_half, q); ky = (Hp2 - Hm2)/(2*dy)
P("  The sheet as a spring at its edge: %.2f N/m2 downward for each meter the edge is raised (H / %.0f m) ; %.2f N/m2 more pull for each meter the edge is drawn out (H / %.0f m)."
  % (kv, d5["H"]/kv, ky, d5["H"]/ky))
roll_restore = 2*kv*S*leg*leg; roll_upset = 2*fs["F"]*leg; k_s = E_CNT*(fs["Ts_MN"]*1e6/S75)/L_LINE; roll_lines = 2*k_s*(2*leg/3)**2
P("  Roll (one tip up, the other down): the push on a tilted strut upsets it with 2 F l = %.2e N m per radian ; the sheet restores %.2e, the suspension lines %.2e.  Ratio restoring / upsetting %.2f."
  % (roll_upset, roll_restore, roll_lines, (roll_restore + roll_lines)/roll_upset))
P("  Yaw (one tip east, the other west): the push upsets with F / l = %.0f N for each meter a tip moves along the ring ; the sheet, whose pull points at the ridge over a horizontal reach equal to the leg," % (fs["F"]/leg))
P("     restores with F / (half footprint) = %.0f N/m, the same.  The lines add %.0f N/m (tension / 100 km).  Net: neutral.  The frame is not held in yaw by the sheet and the lines." % (fs["F"]/b, fs["Ts_MN"]*1e6/L_LINE))
P("     The script adds four CNT stays from each frame's tips to the hinges of the next frames, sized for %.0f %% of the push [ASSUMED]: %.2f t for each frame." % (BRACE*100, fs["stays_t"]))
P("  Drift of the whole array along the ring: as in the book, the pendulum of the support lines, weight / 100 km = %.3f N for each meter of ring for each meter of drift.  The frames add their weight to that and change nothing else."
  % ((0.5*d5["W"]*La/L_RING*g_a)/L_LINE))
OUT["stability"] = dict(kv_N_m2=kv, ky_N_m2=ky, roll_ratio=(roll_restore + roll_lines)/roll_upset, roll_upset=roll_upset, roll_sheet=roll_restore, roll_lines=roll_lines,
                        yaw_upset_N_m=fs["F"]/leg, yaw_sheet_N_m=fs["F"]/b, yaw_lines_N_m=fs["Ts_MN"]*1e6/L_LINE, stays_t=fs["stays_t"])

# ---------------------------------------------------------------------------------------------------------------------
# PART E.  The ring without ground arms
# ---------------------------------------------------------------------------------------------------------------------
P("\n" + "="*110 + "\nPART E.  The ring with no ground arms for the array\n" + "="*110)
M_cas = 12000.0*L_RING
def swing(k): om = math.sqrt(k/M_cas); return 2*math.pi/om, 0.576/om, 0.281/om
for nm, k in (("book, with the array's arms (Table 11.3: 277,000 N/m)", 277000.0), ("no array arms: the lines of Table 11.3's first column (5,800 N/m, the book's figure, not derived again)", 5800.0)):
    T, a1, a2 = swing(k)
    P("  %-98s period %.0f s (%.1f h) ; swing after the design launch %.0f m, %.0f m with the sled braked ; reaction line at that swing %.1f kN of its 93.6 kN"
      % (nm, T, T/3600, a1, a2, 1.55*a1/1e3))
    OUT.setdefault("swing", {})[nm] = dict(k=k, period_s=T, swing_m=a1, swing_braked_m=a2, reaction_kN=1.55*a1/1e3)
a1mv = P9["0.5"]["anchor_1MV"]/P9["0.5"]["P_TW"]; a2mv = P9["0.5"]["anchor_2MV"]/P9["0.5"]["P_TW"]
P("  Conductors.  Anchor lines, all the power by wire: %.1f kg/m for each TW at 1 MV, %.1f at 2 MV (round 9's rows).  Ring lines, upper bound: 749 / 12.2 = %.1f kg/m for each TW at 10 MV, %.1f at 20 MV (the book's rows, in proportion)."
  % (a1mv, a2mv, 749/12.2, 187/12.2))
OUT["conductors"] = {}
for sA in (1.0, 0.5, 0.2):
    pw = OUT["budget"][str(sA)]["P_TW"]
    r = dict(P_TW=pw, anchor_1MV=a1mv*pw, anchor_2MV=a2mv*pw, ring_10MV=749/12.2*pw, ring_20MV=187/12.2*pw, md_pct=41.647/pw)
    OUT["conductors"][str(sA)] = r
    P("    %.1f kg/m2, %.1f TW: anchor lines %.0f kg/m at 1 MV, %.0f at 2 MV ; ring lines %.0f kg/m at 10 MV, %.0f at 20 MV ; the design launch is %.1f %% of the output"
      % (sA, pw, r["anchor_1MV"], r["anchor_2MV"], r["ring_10MV"], r["ring_20MV"], r["md_pct"]))
P("  With ground arms, 420 kg/m of structural arm was there to conduct in any case.  With none, every kilogram of anchor conductor is added mass, on lines that have no other job.")

with open(os.path.join(LOG_DIR, "round14_truss_spreader.txt"), "w", encoding="utf-8") as f: f.write(_buf.getvalue())
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    return o
with open(os.path.join(LOG_DIR, "round14_truss_spreader.json"), "w", encoding="utf-8") as f: json.dump(clean(OUT), f, indent=1)
print("\nlog written: runs/round14_logs/round14_truss_spreader.txt and .json")
