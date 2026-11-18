#!/usr/bin/env python3
"""
Round 15: the truss spreader as the array's design (Chapters 11 and 12 of Volume III).  Seven items of the brief of 17 November 2026:

  7. the shell-buckling formula, from NASA SP-8007 (August 1968) as opened and read, and what it changes;
  1. Paul's opening of the frame (legs hang from the hinge, the center line is paid out), quasi-static, lines at 2/3 and at 4/5,
     and the drawing out of the sheet by the tethers afterwards;
  2. one, two and three lines on each leg, the outermost at 2/3 and at 4/5: frame, widest sheet, power with 1,200 kg/m held;
  3. the array 2, 5, 10, 20 and 100 km below the ring: leaning lines as hangers and as stays ; a ring at 600 km;
  4. the array as the deployment's power: demand and supply, weight, drag and atomic oxygen (NRLMSIS), mass freed, radiation;
  5. the swing along the ring with no ground arms, seven and fourteen design launches a day, with the LIMs damping it.

The sheet and the budget are round 9's (tools/round9_lines_array.py), loaded by name as rounds 13 and 14 did.  The truss model is
round 14's, written out again here with the line layout, the buckling length and the wall-buckling rule as arguments ; it is checked
against round 14's log first.

Run from Python/orbitalring:  python tools/round15_truss_array.py
Log: runs/round15_logs/round15_truss_array.txt and .json
"""
import ast, csv, io, json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import numpy as np
from scipy.optimize import brentq, root, linprog
import ring_altitude as R, power_config as cfg, power_physics as phys

QUICK = "--quick" in sys.argv          # a short pass over every part, to test the script ; its log goes to a folder of its own
LOG_DIR = os.path.join(ROOT, "runs", "round15_logs", "quick") if QUICK else os.path.join(ROOT, "runs", "round15_logs"); os.makedirs(LOG_DIR, exist_ok=True)
_buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
OUT = {}
def head(s): P("\n" + "="*118 + "\n" + s + "\n" + "="*118)

# ---------------------------------------------------------------------------------------------------------------------
# Round 9's sheet and budget, loaded by name ; round 14's log for the checks
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
MU = r9["mu"]; RE = r9["RE"]; W_E = r9["W_E"]
RHO = r9["RHO"]; S_OP = r9["S_OP"]; S75 = r9["S75"]; E_CNT = r9["E_CNT"]
geff = r9["geff"]; pot = r9["pot"]; dia = r9["dia"]; catenary = r9["catenary"]; v2 = r9["v2"]
LOG9 = json.load(open(os.path.join(ROOT, "runs", "round9_logs", "round9_lines_array.json"), encoding="utf-8"))
P9 = LOG9["power"]["catenary"]["rows"]
LOG14 = json.load(open(os.path.join(ROOT, "runs", "round14_logs", "round14_truss_spreader.json"), encoding="utf-8"))
HS45, DR45 = catenary(45.0)
TW_PER_KM_600 = P9["0.5"]["P_TW"]/P9["0.5"]["sheet_km"]
BOOK = {"1.0": (0.54, 5.9), "0.5": (1.11, 12.2), "0.2": (2.83, 31.1)}

def lit(h):
    """Share of the circumference at altitude h that is outside the Earth's shadow (Sun in the plane of the ring)."""
    return 1 - math.degrees(math.asin(RE/(RE + h)))/180.0
class Env:
    """A ring at h_ring and an array at h_array.  For the book's 700 and 600 km it gives round 9's own numbers.
    The power for each km of sheet is round 9's figure at 600 km, scaled by the array's circumference only [the watts for each
    square meter are held ; the book has 293 W/m2 at 600 km and 291 at 700 km]."""
    def __init__(self, h_ring, h_array):
        self.h_r = h_ring; self.h_a = h_array; self.D = h_ring - h_array
        self.g_a = geff(h_array); self.g_r = geff(h_ring)
        self.Lr = 2*math.pi*(RE + h_ring); self.La = 2*math.pi*(RE + h_array)
        self.kc = RHO*self.D*self.g_a/(S75 - RHO*pot(h_array, h_ring))
        self.TW_PER_KM = TW_PER_KM_600*self.La/(2*math.pi*(RE + 600e3))
BOOKENV = Env(r9["h_r"], r9["h_a"])

P("Round 15: the truss spreader as the array's design.")
P("ring %.0f km, array %.0f km ; g at the array %.4f m/s2 ; kc %.4f (round 9: %.4f) ; %.3f TW for each km of sheet ; CNT %.0f kg/m3, %.2f GPa steady lines, %.0f GPa modulus"
  % (BOOKENV.h_r/1e3, BOOKENV.h_a/1e3, BOOKENV.g_a, BOOKENV.kc, r9["kc"], BOOKENV.TW_PER_KM, RHO, S75/1e9, E_CNT/1e9))
assert abs(BOOKENV.kc/r9["kc"] - 1) < 1e-9 and abs(BOOKENV.g_a/r9["g_a"] - 1) < 1e-12
d5 = v2(0.5)

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 7.  Material, sizing rules, and the wall-buckling formula
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 7.  The buckling formula of the tube wall, from the source")
KSI = 6.894757e6
E_UD = 21.7e3*KSI; SC_UD = 245.0*KSI                 # Hexcel HexPly 8552 Product Data Sheet (US, 2020), IM7 UD, 77 F dry: 0 deg compression modulus and strength
E2_UD = 1.70e3*KSI; NU12 = 0.316                     # same table: 90 deg tensile modulus 1.70 msi ; major Poisson's ratio, tension, 0.316
G12 = 5.0e9                                          # in-plane shear modulus of the ply: NOT on the data sheet [UNVERIFIED, a typical value for this class of ply]
RHO_RESIN = 0.0470*27679.9; RHO_FIBRE = 1770.0; VF = 0.60
RHO_CF = VF*RHO_FIBRE + (1 - VF)*RHO_RESIN
AXIAL = 0.75
E_CF = AXIAL*E_UD; SC_CF = AXIAL*SC_UD
FOS = 2.0; T_MIN = 1.0e-3; D_MIN_CNT = 2.0e-3; JOINTS = 0.20; BRACE = 0.02; NU = 0.3

def gamma_shell(r_over_t):
    """NASA SP-8007, 'Buckling of Thin-Walled Circular Cylinders' (revised August 1968), Section 4.2.1, page 5, Equations (5) and (6):
    gamma = 1 - 0.901 (1 - exp(-phi)),  phi = (1/16) sqrt(r/t)  for r/t < 1500.  Read from the scanned document (NTRS 19690013955)."""
    return 1 - 0.901*(1 - np.exp(-np.sqrt(r_over_t)/16.0))
def ply_Q(E1, E2, nu12, G):
    den = 1 - nu12*nu12*E2/E1
    return E1/den, E2/den, nu12*E2/den, G
def laminate(share0, other):
    """Wall stiffnesses of a laminate with a share of 0 deg plies and the rest at 90 deg ('90') or at +-45 deg ('45'), the plies mixed
    through the wall so that the bending stiffness is the membrane stiffness times t^2/12 [ASSUMED].  Returns Qxx, Qyy, Qxy, G (Pa)."""
    Q11, Q22, Q12, Q66 = ply_Q(E_UD, E2_UD, NU12, G12)
    if other == "90": o = (Q22, Q11, Q12, Q66)
    else:
        a = (Q11 + Q22 + 2*Q12 + 4*Q66)/4.0; o = (a, a, (Q11 + Q22 + 2*Q12 - 4*Q66)/4.0, (Q11 + Q22 - 2*Q12)/4.0)
    s = share0
    return s*Q11 + (1 - s)*o[0], s*Q22 + (1 - s)*o[1], s*Q12 + (1 - s)*o[2], s*Q66 + (1 - s)*o[3]
def K_wall(Qxx, Qyy, Qxy, G):
    """SP-8007 Equation (44), page 18, for a layered orthotropic cylinder with no coupling, for a wall of one material through its
    thickness: N l^2 / (pi^2 Dx) = m^2 f(beta) + gamma^2 l^4 g(beta) / (pi^4 m^2 Dx r^2).  Minimised over m (taken as continuous) the load
    is N = (2 gamma / r) sqrt(Dx f g), so the buckling stress is gamma K t / r with
        K = min over beta of sqrt(Qxx f g') / sqrt(3),
        f = 1 + 2 (Qxy + 2 G) / Qxx beta^2 + Qyy / Qxx beta^4,   g' = (Qxx Qyy - Qxy^2) / (Qxx + ((Qxx Qyy - Qxy^2)/G - 2 Qxy) beta^2 + Qyy beta^4).
    For an isotropic wall K = E / sqrt(3 (1 - nu^2)), which is Equation (4).
    NOTE: the twisting stiffness in f is taken as (Qxy + 2 G) t^3 / 12.  Equation (62) as printed gives twice that, and with it Equation (44)
    does not reduce to Equation (2) for an isotropic wall ; the value used here does [my reading of the document]."""
    best = None
    for b in np.concatenate(([0.0], np.geomspace(1e-3, 1e3, 2401))):
        f = 1 + 2*(Qxy + 2*G)/Qxx*b*b + Qyy/Qxx*b**4
        g = (Qxx*Qyy - Qxy*Qxy)/(Qxx + ((Qxx*Qyy - Qxy*Qxy)/G - 2*Qxy)*b*b + Qyy*b**4)
        v = math.sqrt(Qxx*f*g)/math.sqrt(3.0)
        if best is None or v < best[0]: best = (v, b)
    return best
K_ISO14 = E_CF/math.sqrt(3*(1 - NU*NU))               # round 14's rule: isotropic formula with the axial modulus
Eiso = 100e9; chk = K_wall(Eiso/(1 - NU*NU), Eiso/(1 - NU*NU), NU*Eiso/(1 - NU*NU), Eiso/(2*(1 + NU)))
P("  Source opened and read: NASA SP-8007, 'Buckling of Thin-Walled Circular Cylinders', NASA Space Vehicle Design Criteria (Structures),")
P("  revised August 1968, by Weingarten, Seide and Peterson ; NTRS document 19690013955 (a scan ; its pages were read as images).")
P("    Equation (4), p. 4:  sigma = gamma E / sqrt(3 (1 - mu^2)) x t / r  ( = 0.6 gamma E t / r for mu = 0.3 )")
P("    Equations (5), (6), p. 5:  gamma = 1 - 0.901 (1 - e^-phi),  phi = (1/16) sqrt(r/t),  for r/t < 1500.")
P("    The same page: the curve 'should be used with caution for cylinders with length-radius ratios greater than about 5'.")
P("    Equations (44) to (46), pp. 18-19, for layered orthotropic cylinders: the same gamma with phi = (1/29.8) sqrt(r / (Dx Dy / (Ex Ey))^(1/4)),")
P("    which for a wall of one material is sqrt(r/t) / %.2f, the same as Equation (6)." % (29.8/12**0.25))
P("  So the knockdown that round 14 quoted from memory is the document's, to the digit.  What round 14 had wrong is the stiffness that goes with it:")
P("  it used the isotropic Equation (4) with the tube's axial modulus.  A tube wall with most of its fibers along the member is far softer around")
P("  the hoop and in shear, and Equation (44) is the one that applies.")
P("  check of the Equation (44) routine on an isotropic wall (100 GPa, 0.3): K = %.3f GPa against E / sqrt(3 (1 - nu^2)) = %.3f GPa" % (chk[0]/1e9, Eiso/math.sqrt(3*(1 - NU*NU))/1e9))
P("  Ply (Hexcel HexPly 8552 data sheet, IM7, 77 F, dry): E1 %.1f GPa (0 deg compression), E2 %.2f GPa (90 deg tension, 1.70 msi), Poisson 0.316 ; G12 %.1f GPa [UNVERIFIED: not on the sheet]"
  % (E_UD/1e9, E2_UD/1e9, G12/1e9))
OUT["item7"] = dict(source="NASA SP-8007 (rev. August 1968), Eqs (4), (5), (6), (44)-(46); NTRS 19690013955", K_round14_GPa=K_ISO14/1e9, laminates={})
P("    wall                                            | Qxx     | Qyy     | G       | K (sigma = gamma K t / r) | against round 14's %.1f GPa" % (K_ISO14/1e9))
LAMS = {}
for nm, s0, oth in (("75 % at 0 deg, 25 % at 90 deg", 0.75, "90"), ("75 % at 0 deg, 25 % at +-45 deg", 0.75, "45"), ("60 % at 0 deg, 40 % at +-45 deg", 0.60, "45"), ("50 % at 0 deg, 50 % at +-45 deg", 0.50, "45")):
    Q = laminate(s0, oth); K, b = K_wall(*Q); LAMS[nm] = dict(K=K, share0=s0, Q=Q)
    P("    %-47s | %5.1f   | %5.1f   | %5.1f   | %5.1f GPa (beta %.2f)     | %.2f" % (nm, Q[0]/1e9, Q[1]/1e9, Q[3]/1e9, K/1e9, b, K/K_ISO14))
    OUT["item7"]["laminates"][nm] = dict(Qxx_GPa=Q[0]/1e9, Qyy_GPa=Q[1]/1e9, Qxy_GPa=Q[2]/1e9, G_GPa=Q[3]/1e9, K_GPa=K/1e9, ratio_to_round14=K/K_ISO14)
WALL_NAME = "75 % at 0 deg, 25 % at +-45 deg"
K_SRC = LAMS[WALL_NAME]["K"]
P("  Round 15 uses the wall '%s': axial modulus and strength as round 14 (%.1f GPa, %.0f MPa, the 0 deg share of the ply values), K = %.1f GPa." % (WALL_NAME, E_CF/1e9, SC_CF/1e6, K_SRC/1e9))

_T = np.geomspace(T_MIN, 0.08, 90)
def tube(N, Lb, Amin=0.0, E=None, SC=None, K=None):
    """Lightest round tube for a factored compression N over a pin-ended length Lb, with at least the area Amin.
    Euler N <= pi^3 E r^3 t / Lb^2 ; strength N <= 2 pi r t SC ; wall N / (2 pi r t) <= gamma(r/t) K t / r ; t >= T_MIN ; r >= 5 t."""
    E = E_CF if E is None else E; SC = SC_CF if SC is None else SC; K = K_SRC if K is None else K
    t = _T
    r = np.maximum.reduce([(N*Lb*Lb/(math.pi**3*E*t))**(1/3.0), N/(2*math.pi*t*SC), Amin/(2*math.pi*t), 5*t])
    ok = N/(2*math.pi*r*t) <= gamma_shell(r/t)*K*t/r
    A = np.where(ok, 2*math.pi*r*t, np.inf); i = int(np.argmin(A))
    if not np.isfinite(A[i]): return float("nan"), float("nan"), float("inf")
    return float(r[i]), float(t[i]), float(A[i])

def lp_bending(n, f, hinge_free=False):
    """Bending of one leg under its own weight (unit weight for each unit length, unit length), with n lines at f k / n from the hinge
    (k = 1..n) and the center line at the hinge.  The line pulls are chosen to make the largest bending moment least (a linear program):
    the lines are long and soft, so each carries what its length is set to make it carry.  hinge_free: no moment at the hinge (round 14).
    Returns cM (moment / (w l^2)), cV (shear / (w l)), the lines' shares of the leg's weight, the center line's share for this leg."""
    xk = np.array([f*k/n for k in range(1, n + 1)])
    xs = np.unique(np.concatenate((np.linspace(0, 1, 801), xk)))
    def Mrow(x): return np.maximum(xk - x, 0.0)
    A = []; b = []
    for x in xs:
        m0 = -(1 - x)**2/2.0; row = Mrow(x)
        A.append(np.concatenate((row, [-1.0]))); b.append(-m0)          #  m0 + row.R <= t
        A.append(np.concatenate((-row, [-1.0]))); b.append(m0)          # -m0 - row.R <= t
    A.append(np.concatenate((np.ones(n), [0.0]))); b.append(1.0)        # the center line cannot push: sum R <= 1
    Aeq = beq = None
    if hinge_free: Aeq = [np.concatenate((Mrow(0.0), [0.0]))]; beq = [0.5]
    res = linprog(np.concatenate((np.zeros(n), [1.0])), A_ub=np.array(A), b_ub=np.array(b), A_eq=Aeq, b_eq=beq, bounds=[(0, None)]*(n + 1), method="highs")
    Rk = res.x[:n]
    V = []
    for x in np.concatenate((np.linspace(0, 1, 801), xk - 1e-9, xk + 1e-9)):
        V.append(abs((1 - x) - Rk[xk > x].sum()))
    return dict(cM=float(res.x[n]), cV=float(max(V)), R=[float(v) for v in Rk], Rc=float(1 - Rk.sum()), x=[float(v) for v in xk], M0=float(-0.5 + (Rk*xk).sum()))

def truss(Pc, Ltot, Lb=None, cM=1/18.0, cV=0.417, g=None, K=None, E=None, detail=False):
    """Round 14's three-chord truss (equilateral section of side d, bays lb, tubes of the compression material, crossed CNT diagonals),
    for a compression Pc over a total length Ltot (two legs).  Lb: the length over which the whole strut can buckle (Ltot when the lines
    only hang it).  cM, cV: bending and shear from its own weight on its lines, as shares of w l^2 and w l (l the leg).  g: gravity at the array."""
    Lb = Ltot if Lb is None else Lb; E = E_CF if E is None else E; g = BOOKENV.g_a if g is None else g
    leg = Ltot/2; best = None
    for d in np.geomspace(max(Lb/400, 0.3), Lb/4, 34):
        for lb in d*np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0]):
            nb = max(2, int(round(Ltot/lb))); lb_ = Ltot/nb
            if lb_ > Lb: continue
            ld = math.hypot(d, lb_); sph = d/ld; cph = lb_/ld
            m = 0.0; res = None
            for _ in range(6):
                w = m*g/Ltot; M = cM*w*leg*leg/(1 - 1/FOS); Vw = cV*w*leg/(1 - 1/FOS)
                Q = BRACE*Pc + Vw
                Ad = max(Q/(1.5*sph)/S_OP, math.pi*D_MIN_CNT**2/4, 5*FOS*Pc/(1.5*E_CNT*sph*sph*cph))
                GA = 1.5*E_CNT*Ad*sph*sph*cph
                PE = 1.0/(1.0/(FOS*Pc) - 1.0/GA)
                Aeul = 2*PE*Lb*Lb/(math.pi**2*E*d*d)
                Nc = FOS*(Pc/3 + M/(0.866*d))
                rc, tc, Ac = tube(Nc, lb_, Aeul, E, None, K)
                Nb = BRACE*Nc + FOS*Q/(1.5*sph)*sph
                rb, tb, Ab = tube(Nb, d, 0.0, E, None, K)
                m_ch = 3*Ac*Ltot*RHO_CF; m_bat = 3*Ab*d*RHO_CF*(nb + 1); m_di = 6*Ad*ld*RHO*nb
                m_new = (m_ch + m_bat)*(1 + JOINTS) + m_di
                res = dict(d=d, lb=lb_, bays=nb, r_chord=rc, t_chord=tc, A_chord=Ac, r_batten=rb, t_batten=tb, d_diag_mm=2e3*math.sqrt(Ad/math.pi),
                           m_chords=m_ch, m_battens=m_bat, m_diag=m_di, m_joints=(m_ch + m_bat)*JOINTS, mass=m_new,
                           stress_chord_MPa=Nc/FOS/Ac/1e6, euler_governs=bool(Aeul >= Ac*0.999), M_weight_MNm=M/1e6, N_chord_MN=Nc/FOS/1e6)
                if not np.isfinite(m_new) or abs(m_new - m) < 1e-4*m_new: m = m_new; break
                m = m_new
            if res and np.isfinite(m) and (best is None or m < best["mass"]): best = res
    if best is None: best = dict(mass=float("inf"), d=float("nan"))
    return best if detail else best["mass"]

_LP = {}
def layout(n, f, hinge_free=False):
    k = (n, round(f, 4), hinge_free)
    if k not in _LP: _LP[k] = lp_bending(n, f, hinge_free)
    return _LP[k]
BETA_BRACE = 4.0     # a strut on many equal elastic supports is fully braced when each has a stiffness of 4 P / s (s the spacing) [from memory of Timoshenko and Gere ; no source opened: UNVERIFIED]
WINTER = 2.0         # and the stiffness is doubled so that the crookedness of the strut does not overload the braces [Winter's rule, from memory: UNVERIFIED]
def frame_system(sA, W, S, env=BOOKENV, n=1, f=2/3.0, hinge_free=False, braced=False, K=None):
    """One frame for a sheet W wide with panels sA and frames every S.
    braced = False: the lines only hang the frame (round 14) ; it can buckle over its whole length ; four yaw stays are added.
    braced = True: each line point and the hinge is held by a pair of lines splayed along the ring at 45 deg, stiff enough to brace the strut
    in both planes, so that it can buckle only between line points (or over twice the overhang past the last line)."""
    q = sA*env.g_a; H = q*W/2; F = H*S; leg = HS45*(W/2); Lf = 2*leg; Dz = env.D + DR45*(W/2)        # the frame hangs at the height of the edges
    lay = layout(n, f, hinge_free)
    yk = [x*leg for x in lay["x"]]
    Lb = max(f*leg/n, 2*(1 - f)*leg) if braced else Lf
    m = 0.0; t = None
    for _ in range(3 if (braced or env.D < 50e3) else 1):       # at 100 km the lean adds a fifth of a percent ; one pass
        Wl = m*env.g_a/2
        add = sum(Rk*Wl*y/Dz for Rk, y in zip(lay["R"], yk))                 # a leaning line pulls its point of the leg toward the center line
        t = truss(F + add, Lf, Lb, lay["cM"], lay["cV"], env.g_a, K, detail=True)
        if not np.isfinite(t["mass"]) or abs(t["mass"] - m) < 1e-3*t["mass"]: m = t["mass"]; break
        m = t["mass"]
    Wl = m*env.g_a/2; add = sum(Rk*Wl*y/Dz for Rk, y in zip(lay["R"], yk))
    Ts = [Rk*Wl*math.hypot(y, Dz)/Dz for Rk, y in zip(lay["R"], yk)]; Tc = 2*lay["Rc"]*Wl
    TR = q*S*S/2; mR = 2*RHO*TR/S75
    if not braced:
        m_lines = env.kc*m
        stay = 4*BRACE*F/S75*RHO*math.hypot(S, Lf/2)
        k_req = 0.0
    else:
        s = f*leg/n; k_req = WINTER*BETA_BRACE*FOS*(F + add)/s; m_lines = 0.0; stay = 0.0
        for Rv, y in [(lay["Rc"]*2*Wl, 0.0)] + [(Rk*Wl, y) for Rk, y in zip(lay["R"], yk)]*2:
            Lp = math.sqrt(y*y + 2*Dz*Dz); nz = Dz/Lp                         # a pair, each line splayed Dz along the ring
            A = max(k_req*Lp/(2*E_CNT*nz*nz), Rv/(2*nz)/S75*(1 + env.kc))
            m_lines += 2*RHO*A*Lp
    casing = (m + m_lines + stay)/S
    return dict(F=F, Pc=F + add, add=add, Lf=Lf, leg=leg, H=H, truss=t, mass_t=m/1e3, Ts_MN=[v/1e6 for v in Ts], Tc_MN=Tc/1e6, lines_t=m_lines/1e3, stays_t=stay/1e3,
                casing_kgm=casing, TR_MN=TR/1e6, ret_mm=dia(TR, S75), ret_kgm=mR, bow_m=W/8, Lb=Lb, k_req=k_req, Dz=Dz,
                lean_deg=[math.degrees(math.atan(y/Dz)) for y in yk])
def budget(sA, S, M=1200.0, env=BOOKENV, **kw):
    pw = sA*env.La/env.Lr*(1 + env.kc)
    def tot(W):
        fr = frame_system(sA, W, S, env, **kw); return pw*W + fr["ret_kgm"] + fr["casing_kgm"] - M
    if tot(50.0) > 0: return None
    W = brentq(tot, 50.0, M/pw, xtol=0.05); fr = frame_system(sA, W, S, env, **kw)
    return dict(W=W, pan=pw*W, ret=fr["ret_kgm"], frames=fr["casing_kgm"], P_TW=env.TW_PER_KM*W/1e3, f=fr)
def best_budget(sA, grid, **kw):
    bb = [(S, budget(sA, S, **kw)) for S in grid]; bb = [(S, b) for S, b in bb if b]
    if not bb: return None, None
    S, b = max(bb, key=lambda x: x[1]["W"])
    for S2 in (S/1.25, S*1.25):
        b2 = budget(sA, S2, **kw)
        if b2 and b2["W"] > b["W"]: S, b = S2, b2
    return S, b

# ---- checks against round 14, and what the sourced wall formula changes ----
W5 = d5["W"]
c14 = frame_system(0.5, W5, 10e3, hinge_free=True, K=K_ISO14)
r14 = LOG14["design_frame"]["truss"]["mass"]/1e3
P("\n  CHECK. Round 14's design frame with round 14's rules (hinge free of moment, isotropic wall formula): %.1f t here, %.1f t in round 14's log  -> %s"
  % (c14["mass_t"], r14, "reproduced" if abs(c14["mass_t"]/r14 - 1) < 2e-3 else "NOT reproduced"))
OUT["check14"] = dict(here_t=c14["mass_t"], round14_t=r14)
def show(t, nm):
    P("    %-58s %6.1f t ; %.1f m deep, %d bays of %.1f m ; chords %.0f mm x %.1f mm at %.0f MPa (%s) ; battens %.0f mm x %.1f mm"
      % (nm, t["mass"]/1e3, t["d"], t["bays"], t["lb"], 2e3*t["r_chord"], 1e3*t["t_chord"], t["stress_chord_MPa"], "whole strut" if t["euler_governs"] else "chord in its bay", 2e3*t["r_batten"], 1e3*t["t_batten"]))
P("  The design frame (0.5 kg/m2, sheet %.1f m, frames every 10 km, push %.2f MN) under each wall rule, lines at 2/3, hinge free of moment as round 14:" % (W5, c14["F"]/1e6))
show(c14["truss"], "round 14's rule (isotropic, axial modulus)")
OUT["item7"]["design_frame"] = {"round14_rule": c14["mass_t"]}
for nm in LAMS:
    fr = frame_system(0.5, W5, 10e3, hinge_free=True, K=LAMS[nm]["K"]); show(fr["truss"], "SP-8007 Eq. (44), " + nm)
    OUT["item7"]["design_frame"][nm] = dict(mass_t=fr["mass_t"], depth_m=fr["truss"]["d"], chord_mm=2e3*fr["truss"]["r_chord"], wall_mm=1e3*fr["truss"]["t_chord"], bays=fr["truss"]["bays"])
P("    (the 50 and 60 % walls are listed for their K only: the script keeps round 14's axial modulus for every row, and a wall with fewer")
P("     0 deg plies would have a lower one, which is not modelled.)")
base = frame_system(0.5, W5, 10e3, hinge_free=True)
c_r14 = LOG14["design_frame"]["casing_kgm"]
P("  With the wall of round 15 the design frame is %.0f t in place of %.0f t, and it puts %.1f kg/m on the casing in place of %.1f."
  % (base["mass_t"], c14["mass_t"], base["casing_kgm"], c_r14))
P("  The chords' bays are %.0f to %.0f times their radius, far past the 'about 5' to which SP-8007 says its curve was checked ; no better curve was looked for."
  % (base["truss"]["lb"]/base["truss"]["r_chord"], c14["truss"]["lb"]/c14["truss"]["r_chord"]))
OUT["item7"]["design_case"] = dict(frame_t=base["mass_t"], casing_kgm=base["casing_kgm"], round14_frame_t=c14["mass_t"], round14_casing_kgm=c_r14, bay_over_radius=base["truss"]["lb"]/base["truss"]["r_chord"])
P("  The budget row of 0.5 kg/m2 (1,200 kg/m held, one line at 2/3 on each leg, hinge free of moment), round 14's rule against the source's:")
for nm, K in (("round 14's rule", K_ISO14), ("SP-8007 Eq. (44)", None)):
    S, b = best_budget(0.5, (5e3, 7e3, 10e3, 14e3), hinge_free=True, K=K)
    P("    %-18s frames every %.1f km: sheet %.3f km, %.2f TW ; frame %.0f t, %.1f m deep" % (nm, S/1e3, b["W"]/1e3, b["P_TW"], b["f"]["mass_t"], b["f"]["truss"]["d"]))
    OUT["item7"]["budget_0.5_" + nm] = dict(S_km=S/1e3, sheet_km=b["W"]/1e3, P_TW=b["P_TW"], frame_t=b["f"]["mass_t"])

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 1.  Paul's opening
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 1.  Paul's opening: the legs hang from the hinge, the center line is paid out (quasi-static, design case)")
def sheet_at(y, z, s, q):
    """Round 14's hanging sheet: arc s, weight q for each meter of arc, from the ridge at (0, 0) to its edge at (y, z).
    Returns the pull on what holds the edge: H inward and V upward (N for each m of ring), and the height of its center of mass."""
    c = math.hypot(y, z)
    if c >= s*(1 - 1e-9) or y <= 0: return None
    k = math.sqrt(s*s - z*z)/y
    a = y/(2*brentq(lambda u: math.sinh(u)/u - k, 1e-9, 60.0))
    x0 = y/2 - a*math.atanh(z/s)
    ue = (y - x0)/a; u0 = -x0/a
    Hh = q*a; V = -Hh*math.sinh(ue)
    zc = a*a/s*((ue - u0)/2 + (math.sinh(2*ue) - math.sinh(2*u0))/4) - a*math.cosh(u0)
    return Hh, V, zc
s_half = W5/2; qs = 0.5*BOOKENV.g_a; b_half = HS45*s_half; drop = DR45*s_half; S10 = 10e3
OUT["item1"] = {}
P("  The reading: the tethers are slack and the sheet hangs folded from the ridge.  The two legs hang down from the hinge, side by side, with a gap between")
P("  their tips.  Each leg has one suspension line, fastened at a share f of its length from the hinge, of fixed length.  The center line holds the hinge")
P("  and is paid out.  phi is the leg's angle from the vertical: 0 folded, 90 deg open and flat.  Each leg weighs Wl, uniform along its length.")
P("  With the lines vertical, moments about the hinge give the same pulls at every angle:  suspension line Wl / (2 f),  center line 2 Wl (1 - 1 / (2 f)).")
P("  The center line pulls as long as f > 1/2, which is Paul's condition (more of the leg inside the line than outside).  The frame is then a mechanism")
P("  with one freedom, the length of the center line.  It does not snap: it opens as fast as the line is paid out and stops when the line stops.")
def opening(f, D, frame_t, leg, hinge_extra_t=0.0):
    Wl = frame_t*1e3*BOOKENV.g_a/2; Wh = hinge_extra_t*1e3*BOOKENV.g_a
    Zt = D + drop; Ls = math.hypot(f*leg, Zt)
    def zH(p): return Zt - math.sqrt(Ls*Ls - (f*leg*math.sin(p))**2) + f*leg*math.cos(p)
    def U(p): return 2*Wl*(zH(p) - leg/2*math.cos(p)) + Wh*zH(p)
    rows = []
    for ph in (2, 5, 10, 20, 30, 45, 60, 75, 85, 90):
        p_ = math.radians(ph); e = 1e-5
        Tc = (U(p_ + e) - U(p_ - e))/(zH(p_ + e) - zH(p_ - e))            # the line's pull is the energy the frame gives up for each meter the hinge comes down
        yS = f*leg*math.sin(p_); zS = zH(p_) - f*leg*math.cos(p_); nz = (Zt - zS)/Ls
        Tsl = (2*Wl + Wh - Tc)/(2*nz)
        rows.append(dict(phi=ph, hinge_above_open_m=zH(p_), tip_out_m=leg*math.sin(p_), tip_below_open_m=-(zH(p_) - leg*math.cos(p_)), yS_m=yS, lean_deg=math.degrees(math.asin(yS/Ls)),
                         Ts_MN=Tsl/1e6, Tc_MN=Tc/1e6, hinge_squeeze_kN=Tsl*yS/Ls/1e3))
    return dict(rows=rows, Wl_MN=Wl/1e6, travel_m=zH(math.radians(0.0)) - zH(math.pi/2), energy_GJ=(U(0.0) - U(math.pi/2))/1e9)
for f in (2/3.0, 4/5.0):
    fr = frame_system(0.5, W5, 10e3, n=1, f=f); leg = fr["leg"]; d_tr = fr["truss"]["d"]
    P("\n  LINES AT %.3f OF THE LEG.  Frame %.0f t (leg %.1f m, %.1f m deep, each leg %.2f MN) ; lines to the ring 100 km."
      % (f, fr["mass_t"], leg, d_tr, fr["mass_t"]*1e3*BOOKENV.g_a/2/1e6))
    o = opening(f, BOOKENV.D, fr["mass_t"], leg)
    P("    phi   | hinge above its open place | tip: out from the center, below the open frame | line point out | line's lean | suspension line | center line | legs press at the hinge")
    for r_ in o["rows"]:
        P("    %4.0f  | %7.1f m                  | %6.1f m, %6.1f m                              | %6.1f m       | %5.3f deg   | %6.3f MN       | %6.3f MN   | %6.2f kN"
          % (r_["phi"], r_["hinge_above_open_m"], r_["tip_out_m"], r_["tip_below_open_m"], r_["yS_m"], r_["lean_deg"], r_["Ts_MN"], r_["Tc_MN"], r_["hinge_squeeze_kN"]))
    P("    The hinge comes down %.1f m ; the frame's center of mass comes down %.1f m and gives up %.3f GJ, which the brake of the center line takes." % (o["travel_m"], (f - 0.5)*leg, o["energy_GJ"]))
    P("    It falls outward: the center line pulls at every angle (least %.3f MN), so the hinge follows the line down, and the tips can only go apart." % min(r_["Tc_MN"] for r_ in o["rows"]))
    # short lines: the lean of the suspension lines pulls the line points back toward the center
    for D in (20e3, 5e3, 2e3):
        o2 = opening(f, D, fr["mass_t"], leg); tcmin = min(r_["Tc_MN"] for r_ in o2["rows"])
        P("      with lines of %3.0f km: the center line's least pull is %.3f MN (at phi %d deg) %s" % (D/1e3, tcmin, min(o2["rows"], key=lambda r_: r_["Tc_MN"])["phi"], "" if tcmin > 0 else " -> it goes slack: the legs stop opening on their own there"))
        o.setdefault("short_lines", {})["%g" % (D/1e3)] = tcmin
    # the hinge and its stop
    Wl = o["Wl_MN"]*1e6; lever = 0.866*d_tr
    M_slack = Wl*(f - 0.5)*leg
    P("    The hinge.  For the legs to fold downward the hinge pins have to be on the lower side of the truss, and the upper chords then meet end to end when the")
    P("    frame is flat.  That butt joint is a stop: the frame cannot go past flat.  With the center line let go slack, each leg hangs from its own line and its")
    P("    inner part is heavier, so the weight alone holds the stop shut: moment %.0f MN m at the hinge, %.2f MN across the stop and %.2f MN of pull in the pins"
      % (M_slack/1e6, M_slack/lever/1e6, M_slack/lever/1e6))
    P("    (lever %.1f m).  The suspension lines then carry %.2f MN each and the center line nothing.  With the center line at its neutral pull (%.2f MN) the moment is zero."
      % (lever, Wl/1e6, 2*Wl*(1 - 1/(2*f))/1e6))
    P("    Once the sheet is drawn out the push of %.2f MN goes through the stop and the pins as through the rest of the chords (%.2f MN in each chord)." % (fr["F"]/1e6, fr["truss"]["N_chord_MN"]))
    P("    The joint opens again only if the line of thrust leaves the section upward, more than %.1f m above its centroid, which is %.0f MN m of moment." % (0.2887*d_tr, 0.2887*d_tr*fr["F"]/1e6))
    P("    So a latch does not have to carry the push.  It has to hold the stop shut before the sheet is drawn out, if the center line is kept above its")
    P("    neutral pull, and during the drawing out (below).  If the center line were let go instead of paid out, the %.3f GJ would arrive at the stop as a blow." % o["energy_GJ"])
    # drawing the sheet out with the tethers, from the open frame
    P("    Drawing the sheet out.  The retention lines are on the tips ; the tethers run from them to the edges of the hanging sheet and are hauled in.")
    P("    tether  | edge: from the center, below the ridge | tether's slope | push on the frame | weight on each tip | suspension line | center line | moment at the line point | chord force / sized")
    rows = []
    for tl in (585, 570, 540, 500, 450, 400, 300, 200, 100, 50, 20, 5, 0.5):
        def g_(al):
            sh = sheet_at(b_half - tl*math.cos(al), -drop - tl*math.sin(al), s_half, qs)
            if sh is None: return None
            return (-sh[1]) - sh[0]*math.tan(al)
        al_s = np.linspace(1e-4, math.pi/2 - 1e-4, 2000); vals = [g_(a_) for a_ in al_s]; al0 = None
        for i in range(len(al_s) - 1):
            if vals[i] is not None and vals[i + 1] is not None and vals[i]*vals[i + 1] <= 0: al0 = brentq(g_, al_s[i], al_s[i + 1]); break
        if al0 is None: continue
        sh = sheet_at(b_half - tl*math.cos(al0), -drop - tl*math.sin(al0), s_half, qs)
        Fp = sh[0]*S10; Pt = -sh[1]*S10
        Tsl = (Wl/2 + Pt)/f; Tcl = 2*(Wl + Pt - Tsl); Mh = 0.0
        if Tcl < 0: Tsl = Wl + Pt; Tcl = 0.0; Mh = Tsl*f*leg - Wl*leg/2 - Pt*leg          # negative: the tips are heavy and the stop wants to open
        Mtip = Pt*(1 - f)*leg + Wl*(1 - f)**2*leg/2
        Nch = Fp/3 + max(Mtip, abs(Mh))/lever
        rows.append(dict(tether_m=tl, edge_y=b_half - tl*math.cos(al0), edge_below_ridge=drop + tl*math.sin(al0), slope_deg=math.degrees(al0), push_MN=Fp/1e6, tip_load_MN=Pt/1e6,
                         Ts_MN=Tsl/1e6, Tc_MN=Tcl/1e6, M_hinge_MNm=Mh/1e6, M_linepoint_MNm=Mtip/1e6, chord_MN=Nch/1e6, chord_ratio=Nch/1e6/fr["truss"]["N_chord_MN"]))
        P("    %5.1f m | %6.1f m, %6.1f m                      | %5.1f deg      | %6.2f MN         | %6.2f MN          | %6.2f MN       | %5.2f MN    | %7.0f MN m%s | %5.2f MN, %.2f"
          % (tl, rows[-1]["edge_y"], rows[-1]["edge_below_ridge"], rows[-1]["slope_deg"], Fp/1e6, Pt/1e6, Tsl/1e6, Tcl/1e6, Mtip/1e6, " (hinge %+.0f)" % (Mh/1e6) if Mh else "             ", Nch/1e6, rows[-1]["chord_ratio"]))
    tsm = max(r_["Ts_MN"] for r_ in rows); ptm = max(r_["tip_load_MN"] for r_ in rows); crm = max(r_["chord_ratio"] for r_ in rows); mhm = min(r_["M_hinge_MNm"] for r_ in rows)
    m_draw = 2*tsm*1e6/(1 - BOOKENV.kc/(1 + BOOKENV.kc))/S75*RHO*BOOKENV.D
    P("    While it is drawn out, each tip carries up to %.2f MN of the sheet's weight (the open frame carries none).  That is %.1f times the weight of a leg." % (ptm, ptm*1e6/Wl))
    P("    The suspension line then carries up to %.2f MN against %.2f MN open, the center line is slack, the stop is pulled open with up to %.0f MN m (a tie of %.1f MN on the upper chord),"
      % (tsm, fr["Ts_MN"][0], -mhm, -mhm*1e6/lever/1e6))
    P("    and the chords see up to %.2f times the force they are sized for at work (the factor of safety is %.1f)." % (crm, FOS))
    P("    Suspension lines sized for the drawing out: %.1f mm, %.0f t the pair, %.1f kg for each meter of ring at 10 km (open frame alone: %.1f t for all three lines)."
      % (dia(tsm*1e6/(1 - BOOKENV.kc/(1 + BOOKENV.kc)), S75), m_draw/1e3, m_draw/10e3, fr["lines_t"]))
    o.update(frame_t=fr["mass_t"], leg_m=leg, depth_m=d_tr, drop_m=drop, sheet_half_m=s_half, push_MN=fr["F"]/1e6, Ts_open_MN=fr["Ts_MN"][0], Tc_open_MN=fr["Tc_MN"], lines_open_t=fr["lines_t"],
             stop_slack_MNm=M_slack/1e6, stop_force_MN=M_slack/lever/1e6, kern_m=0.2887*d_tr, draw=rows, draw_Ts_max_MN=tsm, draw_tip_max_MN=ptm, draw_chord_ratio=crm,
             draw_hinge_MNm=mhm, draw_tie_MN=-mhm/lever, draw_lines_t=m_draw/1e3, draw_lines_kgm=m_draw/10e3, N_chord_MN=fr["truss"]["N_chord_MN"])
    OUT["item1"]["%.3f" % f] = o
E_sheet = qs*s_half*S10*(sheet_at(b_half, -drop, s_half, qs)[2] + s_half/2)
P("\n  The work of drawing one half of the sheet out over one 10 km span is %.2f GJ (round 14), done by the tether winches: %.0f kJ for each meter of edge." % (E_sheet/1e9, E_sheet/S10/1e3))
P("  The weight on the tips comes from the geometry, not from the frame: the retention line is at the height of the open edge, the hanging edge starts %.0f m lower,"
  % (s_half - drop))
P("  so the tether slopes down and the retention line carries part of the sheet until the edge is up.  The same was found in round 14 for the other way of opening.")
P("  Ways around it that were NOT modelled: a line from the ring to each tip, used only for the drawing out ; hauling the tethers in from the ridge end first ;")
P("  drawing the sheet out before the casing has its full weight (item 4: the strip weighs a twenty-second of its final weight on the first day).")
OUT["item1"]["work_sheet_GJ"] = E_sheet/1e9

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 2.  More lines, wider frames
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 2.  One, two and three lines on each leg, the outermost at 2/3 and at 4/5 of the leg (lines 100 km long, hanging the frame only)")
P("  The lines sit at equal steps from the hinge out to the outermost.  The pulls are set to make the bending from the frame's own weight least ;")
P("  the hinge carries moment once it is shut (item 1).  cM is the largest moment as a share of w l^2, cV the largest shear as a share of w l.")
P("    lines | outermost | cM       | cV    | shares of a leg's weight on its lines (from the hinge out) | center line (both legs)")
OUT["item2"] = dict(layouts={}, rows={})
LAYS = [(1, 2/3.0), (1, 0.8), (2, 2/3.0), (2, 0.8), (3, 2/3.0), (3, 0.8)]
for n, f in LAYS:
    l_ = layout(n, f)
    P("    %d     | %.3f     | %.5f  | %.3f | %-58s | %.3f" % (n, f, l_["cM"], l_["cV"], ", ".join("%.3f" % v for v in l_["R"]), 2*l_["Rc"]))
    OUT["item2"]["layouts"]["%d_%.3f" % (n, f)] = l_
l0 = layout(1, 2/3.0, True)
P("    (round 14: one line at 2/3 with no moment at the hinge: cM %.5f = 1/18, cV %.3f)" % (l0["cM"], l0["cV"]))
xs = np.linspace(0.55, 0.95, 81); cm = [lp_bending(1, x)["cM"] for x in xs]; i = int(np.argmin(cm))
P("    One line does least bending at %.3f of the leg (cM %.5f) ; 2/3 gives %.5f and 4/5 gives %.5f." % (xs[i], cm[i], layout(1, 2/3.0)["cM"], layout(1, 0.8)["cM"]))
OUT["item2"]["best_single"] = dict(f=float(xs[i]), cM=float(cm[i]))
GRID = {1.0: (3e3, 5e3, 8e3), 0.5: (7e3, 11e3, 17e3), 0.2: (18e3, 28e3, 42e3)}
if QUICK: GRID = {1.0: (5e3,), 0.5: (11e3,), 0.2: (28e3,)}; LAYS = LAYS[:2]
for sA in (1.0, 0.5, 0.2):
    P("\n  %.1f kg/m2, 1,200 kg/m held   (book, ground arms: %.2f km, %.1f TW ; round 14, one line at 2/3: %.2f km, %.1f TW)" % (sA, BOOK[str(sA)][0], BOOK[str(sA)][1], LOG14["budget"][str(sA)]["sheet_km"], LOG14["budget"][str(sA)]["P_TW"]))
    P("    lines | outermost | frames every | sheet    | power    | frame: mass, length, depth  | push     | panels + retention + frames (kg/m) | pulls: lines (MN, hinge out), center")
    OUT["item2"]["rows"][str(sA)] = {}
    for n, f in LAYS:
        S, b = best_budget(sA, GRID[sA], n=n, f=f)
        if b is None: P("    %d     | %.3f     | no frame found" % (n, f)); continue
        fr = b["f"]
        P("    %d     | %.3f     | %5.1f km     | %.3f km | %5.2f TW | %6.0f t, %5.0f m, %5.1f m    | %5.1f MN | %.0f + %.0f + %.0f                  | %s ; %.2f"
          % (n, f, S/1e3, b["W"]/1e3, b["P_TW"], fr["mass_t"], fr["Lf"], fr["truss"]["d"], fr["F"]/1e6, b["pan"], b["ret"], b["frames"], ", ".join("%.2f" % v for v in fr["Ts_MN"]), fr["Tc_MN"]))
        OUT["item2"]["rows"][str(sA)]["%d_%.3f" % (n, f)] = dict(S_km=S/1e3, sheet_km=b["W"]/1e3, P_TW=b["P_TW"], frame_t=fr["mass_t"], Lf=fr["Lf"], depth_m=fr["truss"]["d"], push_MN=fr["F"]/1e6,
                                                                 pan=b["pan"], ret=b["ret"], frames=b["frames"], Ts_MN=fr["Ts_MN"], Tc_MN=fr["Tc_MN"])
P("\n  Reading: more lines cut the bending from the frame's own weight, which matters only where the frame is long.  They do nothing for the buckling of the strut")
P("  under the sheet's push, because a line 100 km long is too soft to hold a point of the frame in place (round 14).  That is item 3.")

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 3.  The array a few kilometers below the ring
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 3.  The array 2, 5, 10, 20 and 100 km below the ring: the lines lean, and short lines are stiff")
P("  What a leaning line does.  Every line runs from the center line of the ring down and outward to its point on the leg, so its pull on the leg has a part")
P("  that points back toward the center.  That ADDS to the compression of the frame inside the line point ; it never takes from it.  Leaning lines cannot")
P("  spread the sheet.  What a short line gives is stiffness: E A / length.  A point of the frame held by a stiff enough line cannot move, and the strut can")
P("  then buckle only between such points.  A single line holds its point up and down only.  To hold it along the ring as well (the frame buckles that way")
P("  just as easily, and round 14 found it has no hold in yaw), the script gives each point a pair of lines splayed along the ring, each reaching as far")
P("  along the ring as the array is below it.  Stiffness asked of each pair, in both directions: %.0f x %.0f x (factor of safety) x compression / spacing"
  % (WINTER, BETA_BRACE))
P("  [the bracing rule is from memory, UNVERIFIED ; the top ends are taken as fixed, which the force spreader cable on its hangers is not: ASSUMED].")
OUT["item3"] = dict(design={}, budget={})
P("\n  The design case first (0.5 kg/m2, sheet %.1f m, frames every 10 km, push %.2f MN):" % (W5, c14["F"]/1e6))
P("    below  | lines as hangers only (one at 2/3 on each leg, vertical plane)              | lines as stays (pairs, both planes): best layout")
P("    ring   | lean    | added compression | frame   | lines  | on casing  | lines, outermost | buckling length | pair stiffness | added compr. | frame   | lines     | on casing")
for D in (1e3, 2e3, 5e3, 10e3, 20e3, 100e3):
    env = Env(700e3, 700e3 - D)
    h = frame_system(0.5, W5, 10e3, env)
    best = None
    for n in (1, 2, 3, 4, 6, 8):
        for f in (2/3.0, 0.8, 1.0):
            x = frame_system(0.5, W5, 10e3, env, n=n, f=f, braced=True)
            if best is None or x["casing_kgm"] < best[2]["casing_kgm"]: best = (n, f, x)
    n, f, x = best
    P("    %3.0f km | %5.2f d | %6.3f MN (%.2f %%) | %5.0f t | %5.1f t | %5.1f kg/m | %d, at %.2f       | %6.1f m        | %7.0f kN/m   | %5.2f MN     | %5.0f t | %7.1f t | %6.1f kg/m"
      % (D/1e3, h["lean_deg"][0], h["add"]/1e6, 100*h["add"]/h["F"], h["mass_t"], h["lines_t"], h["casing_kgm"], n, f, x["Lb"], x["k_req"]/1e3, x["add"]/1e6, x["mass_t"], x["lines_t"], x["casing_kgm"]))
    OUT["item3"]["design"]["%g" % (D/1e3)] = dict(hang=dict(lean_deg=h["lean_deg"][0], add_MN=h["add"]/1e6, frame_t=h["mass_t"], lines_t=h["lines_t"], casing_kgm=h["casing_kgm"]),
                                                   stay=dict(n=n, f=f, Lb=x["Lb"], k_kN_m=x["k_req"]/1e3, add_MN=x["add"]/1e6, frame_t=x["mass_t"], lines_t=x["lines_t"], casing_kgm=x["casing_kgm"],
                                                             depth_m=x["truss"]["d"], chord_mm=2e3*x["truss"]["r_chord"]))
P("  The stiff lines' mass goes with the square of their length for the same stiffness, which is why they pay only when the array is close under the ring.")
P("\n  The widest sheet with 1,200 kg/m held (hangers: three lines, the outermost at 4/5 ; stays: the layout that is lightest for that panel and distance, outermost at 4/5 or at the tip):")
P("    panels    | below ring | lines as hangers: frames every, sheet, power, frame           | lines as stays: layout, frames every, sheet, power, frame, lines")
def stay_best(sA, env, W0):
    """The stay layout for a panel mass and a distance: the lightest on the casing at the width W0, then the widest sheet over the spacings."""
    best = None
    for n in (1, 2, 3, 4, 6, 8):
        for f in (0.8, 1.0):
            x = frame_system(sA, W0, GRID[sA][len(GRID[sA])//2], env, n=n, f=f, braced=True)
            if best is None or x["casing_kgm"] + x["ret_kgm"] < best[2]: best = (n, f, x["casing_kgm"] + x["ret_kgm"])
    n, f, _ = best; sb = None
    for S in GRID[sA]:
        b = budget(sA, S, env=env, n=n, f=f, braced=True)
        if b and (sb is None or b["W"] > sb[3]["W"]): sb = (n, f, S, b)
    return sb
def hang_best(sA, env):
    hb = None
    for S in GRID[sA]:
        b = budget(sA, S, env=env, n=3, f=0.8)
        if b and (hb is None or b["W"] > hb[2]["W"]): hb = (3, S, b)
    return hb
for sA in ((0.5,) if QUICK else (1.0, 0.5, 0.2)):
    OUT["item3"]["budget"][str(sA)] = {}
    for D in ((5e3, 100e3) if QUICK else (2e3, 5e3, 10e3, 20e3, 100e3)):
        env = Env(700e3, 700e3 - D)
        n_h, S_h, b_h = hang_best(sA, env); sb = stay_best(sA, env, LOG14["budget"][str(sA)]["sheet_km"]*1e3)
        row = dict(hang=dict(n=n_h, S_km=S_h/1e3, sheet_km=b_h["W"]/1e3, P_TW=b_h["P_TW"], frame_t=b_h["f"]["mass_t"], Lf=b_h["f"]["Lf"]))
        if sb is None: st = "no sheet fits: the stiff lines alone are over the 1,200 kg/m"
        else:
            n_s, f_s, S_s, b_s = sb
            st = ("%d at %.1f, %4.1f km: %.3f km, %5.2f TW, frame %5.0f t (%.0f m long, %.1f m deep, buckling length %.0f m), lines %.0f t"
                  % (n_s, f_s, S_s/1e3, b_s["W"]/1e3, b_s["P_TW"], b_s["f"]["mass_t"], b_s["f"]["Lf"], b_s["f"]["truss"]["d"], b_s["f"]["Lb"], b_s["f"]["lines_t"]))
            row["stay"] = dict(n=n_s, f=f_s, S_km=S_s/1e3, sheet_km=b_s["W"]/1e3, P_TW=b_s["P_TW"], frame_t=b_s["f"]["mass_t"], lines_t=b_s["f"]["lines_t"],
                               Lf=b_s["f"]["Lf"], depth_m=b_s["f"]["truss"]["d"], pan=b_s["pan"], ret=b_s["ret"], frames=b_s["frames"], Lb=b_s["f"]["Lb"])
        P("    %.1f kg/m2 | %3.0f km     | %d lines, %4.1f km: %.3f km, %5.2f TW, frame %5.0f t            | %s"
          % (sA, D/1e3, n_h, S_h/1e3, b_h["W"]/1e3, b_h["P_TW"], b_h["f"]["mass_t"], st))
        OUT["item3"]["budget"][str(sA)]["%g" % (D/1e3)] = row
    P("      floor: with no frames and no retention lines at all, 1,200 kg/m of panels and support lines is %.3f km at 100 km below and %.3f km at 2 km below."
      % (1200/(sA*BOOKENV.La/BOOKENV.Lr*(1 + BOOKENV.kc))/1e3, 1200/(sA*Env(700e3, 698e3).La/BOOKENV.Lr*(1 + Env(700e3, 698e3).kc))/1e3))
env6 = Env(600e3, 595e3)
P("\n  A ring at 600 km with the array 5 km below it (595 km), 0.5 kg/m2, 1,200 kg/m held.  g at the array %.4f m/s2 ; support lines %.2f %% of what they hold ; %.3f TW for each km of sheet."
  % (env6.g_a, env6.kc*100, env6.TW_PER_KM))
_n6, S6, b6 = hang_best(0.5, env6)
x6 = stay_best(0.5, env6, LOG14["budget"]["0.5"]["sheet_km"]*1e3)
P("    lines as hangers (three, the outermost at 4/5): frames every %.1f km, sheet %.3f km, %.2f TW, frame %.0f t" % (S6/1e3, b6["W"]/1e3, b6["P_TW"], b6["f"]["mass_t"]))
P("    lines as stays (%d on each leg, outermost at %.1f): frames every %.1f km, sheet %.3f km, footprint %.3f km, %.2f TW ; frame %.0f t, %.0f m long, %.1f m deep ; lines %.0f t ; push %.1f MN"
  % (x6[0], x6[1], x6[2]/1e3, x6[3]["W"]/1e3, HS45*x6[3]["W"]/1e3, x6[3]["P_TW"], x6[3]["f"]["mass_t"], x6[3]["f"]["Lf"], x6[3]["f"]["truss"]["d"], x6[3]["f"]["lines_t"], x6[3]["f"]["F"]/1e6))
P("    panels and support lines %.0f + retention %.0f + frames and lines %.0f kg/m.   [The cable of a ring at 600 km was not sized here ; only the array's 1,200 kg/m.]" % (x6[3]["pan"], x6[3]["ret"], x6[3]["frames"]))
OUT["item3"]["ring600"] = dict(g_a=env6.g_a, kc=env6.kc, TW_per_km=env6.TW_PER_KM, hang=dict(S_km=S6/1e3, sheet_km=b6["W"]/1e3, P_TW=b6["P_TW"], frame_t=b6["f"]["mass_t"]),
                               stay=dict(n=x6[0], f=x6[1], S_km=x6[2]/1e3, sheet_km=x6[3]["W"]/1e3, P_TW=x6[3]["P_TW"], frame_t=x6[3]["f"]["mass_t"], lines_t=x6[3]["f"]["lines_t"], Lf=x6[3]["f"]["Lf"],
                                         depth_m=x6[3]["f"]["truss"]["d"], pan=x6[3]["pan"], ret=x6[3]["ret"], frames=x6[3]["frames"], push_MN=x6[3]["f"]["F"]/1e6))
# the clearance
SUN = math.radians(0.533)
P("\n  The 20 km of Section 12.7.  It is not in the Python repository.  The only trace is in the old chapter script (Vol_III/renumber_scripts_2026-10-12/ch16_numbers.py),")
P("  which prints the support lines' mass for an array 20, 50 and 100 km below, and the shadow of the ring-mounted panels.  The full shadow of a strip of width b")
P("  ends b / (the Sun's angular width, 0.533 deg) below it: %.1f km for 55 m, %.1f km for 83 m, %.1f km for 110 m.  That is the '6 to 12 km' of Section 12.7,"
  % (55/SUN/1e3, 83/SUN/1e3, 110/SUN/1e3))
P("  and 20 km is a round figure beyond it [my reading ; the book gives no derivation].  With no panels on the ring what is left is the casing: a body %.0f m wide"
  % 12.0)
P("  [ASSUMED width] has a full shadow %.1f km long.  Below that it takes at most its own width out of the strip's light, %.1f %% of a 0.98 km footprint,"
  % (12/SUN/1e3, 100*12/(2*b_half)))
OUT["item3"]["clearance"] = dict(umbra_km={"55": 55/SUN/1e3, "83": 83/SUN/1e3, "110": 110/SUN/1e3, "casing_12m": 12/SUN/1e3}, days={})
ds = []
for D in (2e3, 5e3, 10e3, 20e3, 100e3):
    x = math.degrees(math.atan(b_half/D))/23.44; dd = 365.25*(2/math.pi)*math.asin(min(1.0, x)); ds.append("%.0f km: %.0f" % (D/1e3, dd)); OUT["item3"]["clearance"]["days"]["%g" % (D/1e3)] = dd
P("  on the days when the Sun stands close enough to the equator for the shadow to fall on the strip.  Days a year for the design strip: " + " ; ".join(ds) + ".")
P("  So with no panels on the ring the shadow does not ask for 20 km.  What does set a least distance was not worked out: the lines' own geometry, the radiators,")
P("  hardware falling from the ring, and the radiation, which at 2 to 20 km below a 700 km ring is nearly the ring's own (1,708 against 731 protons/cm2/s at the worst longitude).")

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 4.  The array as the deployment's power
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 4.  The array as the deployment's power, with no panels on the ring")
RUN = os.path.join(ROOT, "runs", "round9_design", "ch9_optimal", "L2_N165_P12MW")
rows = list(csv.DictReader(open(os.path.join(RUN, "rising_radius_daily.csv"))))
day = np.array([float(r_["day"]) for r_ in rows]); rr = np.array([float(r_["r_m"]) for r_ in rows]); vc = np.array([float(r_["v_casing_m_s"]) for r_ in rows])
wk = np.array([float(r_["casing_net_weight_kN_m"]) for r_ in rows])
txt = open(os.path.join(RUN, "rising_radius_lim_simulation.txt"), encoding="utf-8").read()
mon = [(float(m.group(1)), float(m.group(2)), float(m.group(3))) for m in re.finditer(r"^(\d+) mth\s+(-?[\d.]+) m/s.*?([\d.]+) MW\s*$", txt, re.M)]
fin = re.search(r"^final\s+(-?[\d.]+) m/s.*?([\d.]+) MW\s*$", txt, re.M)
md = [float(np.interp(-v, -vc, day)) for _, v, _ in mon] + [day[-1]]; mp = [p for _, _, p in mon] + [float(fin.group(2))]
N_SITES = 88946
P("  Run: runs/round9_design/ch9_optimal/L2_N165_P12MW (165 turns, 12 MW a site, the book's deployment): %d days, ring from %.1f to %.1f km." % (round(day[-1]), (rr[0] - RE)/1e3, (rr[-1] - RE)/1e3))
P("  The site power is the run's monthly table (13 rows and the last day), placed on the days by the casing's speed ; between them it is interpolated.")
import pymsis, warnings
warnings.filterwarnings("ignore")
ALT = np.arange(400.0, 721.0, 10.0)
def msis(F):
    acc = []
    for hr in range(0, 24, 3):
        o = pymsis.calculate(np.array([np.datetime64("2020-03-21T00:00") + np.timedelta64(hr, "h")]), 0., 0., ALT, F, F, [[15]*7])
        acc.append(np.squeeze(o))
    return np.mean(acc, axis=0)
phase = (np.arange(22) + 0.5)/22*np.pi
CYC = np.mean([msis(70 + 130*math.sin(p)**2) for p in phase], axis=0); SMAX = msis(230.0); SMIN = msis(70.0)
KB = 1.380649e-23
def air(h_km, tab=None):
    tab = CYC if tab is None else tab
    rho = math.exp(np.interp(h_km, ALT, np.log(tab[:, 0]))); nO = math.exp(np.interp(h_km, ALT, np.log(tab[:, 3]))); T = float(np.interp(h_km, ALT, tab[:, 10]))
    ntot = sum(math.exp(np.interp(h_km, ALT, np.log(np.maximum(tab[:, j], 1e-30)))) for j in (1, 2, 3, 4, 5, 6, 7))
    mbar = rho/ntot
    return rho, nO, T, math.sqrt(8*KB*T/(math.pi*mbar)), math.sqrt(8*KB*T/(math.pi*16*1.66054e-27))
P("  Air: NRLMSIS 2.1 through pymsis as altitude_trade.py runs it (equator, 21 March, eight local times, averaged over a solar cycle, F10.7 from 70 to 200).")
P("  The strip lies along its own motion, so the air only grazes it.  In free molecular flow a face along the stream is struck by rho c/4 for each square meter")
P("  (c the mean thermal speed), whatever the speed of the stream, and each molecule leaves its share of the stream's momentum (full accommodation ASSUMED):")
P("  drag on each face rho (c/4) v.  Faces that look forward (frame tubes, panel edges, wrinkles) are struck by n v instead.")
SH = W5                                           # sheet width
P("\n    day | ring   | strip  | casing speed | strip through | weight at  | of final | slowing leans | demand  | strip gives   | air at strip | drag on the strip      | O atoms on a face  | struck at | O on a forward")
P("        | km     | km     | (west +) m/s | the air m/s   | strip m/s2 |          | the lines by  | GW      | GW (x demand) | kg/m3        | N/m of ring, MN, GW    | /cm2/s (thermal)   | eV        | face /cm2/s")
tab4 = []; flu = 0.0; flu_ram = 0.0; imp = 0.0; e_drag = 0.0; prev = None
g_fin = None
for d_ in list(range(0, int(day[-1]) + 1)) + [float(day[-1])]:
    r_c = float(np.interp(d_, day, rr)); v = float(np.interp(d_, day, vc)); om = v/r_c
    r_s = r_c - 100e3; h_s = (r_s - RE)/1e3
    g_s = MU/r_s**2 - om*om*r_s
    v_air = (om + W_E)*r_s
    a_dec = -(float(np.interp(d_ + 0.5, day, vc)) - float(np.interp(d_ - 0.5, day, vc)))/86400.0 if 0.5 <= d_ <= day[-1] - 0.5 else float("nan")
    rho, nO, T, cb, cO = air(h_s)
    tau = rho*cb/4*abs(v_air); drag_m = 2*SH*tau; drag = drag_m*2*math.pi*r_s; pdrag = drag*abs(v_air)
    fth = nO*cO/4/1e4; fram = nO*abs(v_air)/1e4; eV = 0.5*16*1.66054e-27*v_air**2/1.602177e-19
    dem = float(np.interp(d_, md, mp))*N_SITES/1e3
    sup = TW_PER_KM_600*(SH/1e3)*1e3*(2*math.pi*r_s)/(2*math.pi*(RE + 600e3))*lit(h_s*1e3)/lit(600e3)
    cur = dict(day=d_, ring_km=(r_c - RE)/1e3, strip_km=h_s, v_casing=v, v_air=v_air, g_strip=g_s, lean_rad=a_dec/g_s, demand_GW=dem, supply_GW=sup, rho=rho, T=T, drag_N_m=drag_m, drag_MN=drag/1e6,
               drag_GW=pdrag/1e9, AO_thermal=fth, AO_ram=fram, eV=eV)
    if prev is not None:
        dt = (d_ - prev["day"])*86400.0
        flu += 0.5*(fth + prev["AO_thermal"])*dt; flu_ram += 0.5*(fram + prev["AO_ram"])*dt; imp += 0.5*(drag + prev["drag_MN"]*1e6)*dt; e_drag += 0.5*(pdrag + prev["drag_GW"]*1e9)*dt
    prev = cur; tab4.append(cur)
g_fin = tab4[-1]["g_strip"]
for c in tab4:
    if c["day"] in (0, 1, 7, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 350) or c is tab4[-1]:
        P("    %3.0f | %6.1f | %6.1f | %8.1f     | %8.1f      | %6.3f     | 1/%-5.1f  | %s | %7.1f | %7.0f (%4.1f) | %.2e     | %8.5f, %6.3f, %6.3f | %.2e           | %5.2f     | %.2e"
          % (c["day"], c["ring_km"], c["strip_km"], c["v_casing"], c["v_air"], c["g_strip"], g_fin/c["g_strip"], ("%.1e rad  " % c["lean_rad"]) if np.isfinite(c["lean_rad"]) else "n/a          ",
             c["demand_GW"], c["supply_GW"], c["supply_GW"]/c["demand_GW"], c["rho"], c["drag_N_m"], c["drag_MN"], c["drag_GW"], c["AO_thermal"], c["eV"], c["AO_ram"]))
OUT["item4"] = dict(table=[c for c in tab4 if c["day"] in (0, 1, 7, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 350) or c is tab4[-1]])
dmax = max(tab4, key=lambda c: c["demand_GW"]); ratio_min = min(c["supply_GW"]/c["demand_GW"] for c in tab4)
need_m = dmax["demand_GW"]/ (dmax["supply_GW"]/SH)
E_dem = sum(0.5*(a["demand_GW"] + b_["demand_GW"])*(b_["day"] - a["day"])*86400e9 for a, b_ in zip(tab4[:-1], tab4[1:]))
P("\n  (a) Power.  The deployment draws at most %.0f GW (from about day %d) and %.1f EJ in all.  A spread strip %.2f km wide gives %.0f to %.0f GW, never less than %.1f times the demand."
  % (dmax["demand_GW"], dmax["day"], E_dem/1e18, SH/1e3, min(c["supply_GW"] for c in tab4), max(c["supply_GW"] for c in tab4), ratio_min))
P("      %.0f m of the sheet's width would cover the peak (the ring-mounted panels of Chapter 8 are 83 m, lying flat).  The strip's output is scaled from 600 km by its"
  % need_m)
P("      circumference and by its share outside the Earth's shadow (%.1f %% at %.0f km, %.1f %% at 600 km) ; the light the Earth reflects was not scaled." % (100*lit(tab4[0]["strip_km"]*1e3), tab4[0]["strip_km"], 100*lit(600e3)))
P("      This holds only if the strip is spread.  Folded, as Section 12.11 has it, it shows the Sun almost nothing.  With truss frames it can be spread before the ring is anchored.")
P("      The power still has to go around the ring to the night side, on the HVDC lines of Chapter 8, and up the array support lines to the casing.")
P("  (b) Weight.  On the first day the strip, 100 km below the casing, weighs %.3f m/s2, 1/%.1f of its final %.3f m/s2.  Section 12.11 has 0.34 m/s2 and 'a twenty-fourth':"
  % (tab4[0]["g_strip"], g_fin/tab4[0]["g_strip"], g_fin))
P("      that is 3 GM / r^3 x 100 km at 700 km (%.3f).  The ring is built at %.0f km, where the same is %.3f, and the exact difference (not the first term only) is %.3f."
  % (3*MU/(RE + 700e3)**3*1e5, (rr[0] - RE)/1e3, 3*MU/rr[0]**3*1e5, tab4[0]["g_strip"]))
c30 = next(c for c in tab4 if c["day"] == 30); c90 = next(c for c in tab4 if c["day"] == 90); c180 = next(c for c in tab4 if c["day"] == 180)
P("      It never goes to zero or changes sign: %.2f m/s2 on day 30, %.2f on day 90, %.2f on day 180.  So the sheet hangs in its tent shape and pulls on its tethers from the"
  % (c30["g_strip"], c90["g_strip"], c180["g_strip"]))
P("      first day, and the frames are in compression from the first day.  The shape does not depend on the weight ; the forces go with it: push on a frame %.2f MN on day one,"
  % (c14["F"]/1e6*tab4[0]["g_strip"]/BOOKENV.g_a))
lean_max = np.nanmax([abs(c["lean_rad"]) for c in tab4])
P("      %.2f MN at the end.  The casing's slowing leans the 100 km lines along the ring by at most %.1e rad (%.0f m at the strip)." % (c14["F"]/1e6*g_fin/BOOKENV.g_a, lean_max, lean_max*1e5))
yr = 365.25*86400
r600 = air(600.0); f600 = r600[1]*r600[4]/4/1e4
P("  (c) Air.  The casing starts at orbital speed going west, and the air turns east with the Earth, so on the first day the strip goes through the air at %.0f m/s."
  % tab4[0]["v_air"])
P("      Drag on the spread strip (both faces of a %.2f km sheet): %.3f MN on the first day, %.3f MN on day 90, against the LIMs' %.0f MN ; over the deployment %.2e N s,"
  % (SH/1e3, tab4[0]["drag_MN"], c90["drag_MN"], 71.4, imp))
P("      which is %.2f %% of the casing's change of momentum (%.2e N s), and %.3f EJ of heat in all.  The drag HELPS the deployment: it slows the casing." % (100*imp/(12000*BOOKENV.Lr*(vc[0] - vc[-1])), 12000*BOOKENV.Lr*(vc[0] - vc[-1]), e_drag/1e18))
P("      Heating of a face on the first day: %.3f W/m2.  Nothing." % (tab4[0]["rho"]*air(tab4[0]["strip_km"])[3]/4*tab4[0]["v_air"]**2/2))
P("      Atomic oxygen on a face of the sheet over the %d days: %.2e atoms/cm2.  The same face at rest at 600 km takes %.2e in a year (the book: about 8e18)."
  % (round(day[-1]), flu, f600*yr))
P("      So in number the deployment costs the sheet %.1f years of its life at 600 km.  But for the first %d days the atoms strike at more than 1 eV (5.3 eV on the first day,"
  % (flu/(f600*yr), max(c["day"] for c in tab4 if c["eV"] > 1.0)))
P("      which is what a satellite takes), against about 0.1 eV at rest.  How much faster a grazing 5 eV atom erodes a coated film than a thermal one is NOT known to me.")
P("      A forward face (a frame tube, a panel's edge) takes %.2e atoms/cm2 over the deployment at those energies ; a satellite at 600 km takes about %.1e in a year."
  % (flu_ram, r600[1]*7560.0/1e4*yr))
P("      Solar maximum and minimum at the first day's %.0f km: density %.1e and %.1e kg/m3 against the cycle average %.1e." % (tab4[0]["strip_km"], air(tab4[0]["strip_km"], SMAX)[0], air(tab4[0]["strip_km"], SMIN)[0], tab4[0]["rho"]))
P("      A FOLDED strip escapes none of this: its sheets also lie along the motion.  Folding saves the Sun's side, not the air's.")
CAB = 103569.71/12000.0
P("  (d) Mass.  Chapter 8 gives the ring-mounted panels as 83 m wide at about 0.5 kg/m2 (film, substrate, wiring): %.1f kg for each meter of ring, %.2f Mt around the ring."
  % (83*0.5, 83*0.5*BOOKENV.Lr/1e9))
P("      It gives no mass for the ribs that carry them, so their structure cannot be counted from Chapter 8.  At %.2f kg of cable for each kg of casing (Section 12.14),"
  % CAB)
P("      the %.1f kg/m frees %.0f kg/m of cable, %.1f Mt.  That is %.2f %% of the 12,000 kg/m.  The HVDC lines around the ring (65 kg/m) stay."
  % (41.5, 41.5*CAB, 41.5*CAB*BOOKENV.Lr/1e9, 41.5/120))
FL = {"700": (1708.0, 196.0), "600": (731.0, 62.0)}
def loss(fl):
    """Share of the output left after a fluence of protons above 10 MeV, between the three points Chapter 5 takes from the cell's data sheet:
    0.975 at 1.2e10, 0.90 at 9e10, 0.80 at 4e11, straight lines against the logarithm of the fluence [my interpolation ; Chapter 5 uses only the two thresholds]."""
    xs_ = np.log([1.2e10, 9e10, 4e11]); ys = [0.975, 0.90, 0.80]
    if fl <= 1.2e10: return 1 - 0.025*fl/1.2e10
    return float(np.interp(math.log(fl), xs_, ys)) if fl <= 4e11 else float(0.80 - 0.10*(math.log(fl/4e11)/math.log(4e11/9e10)))
P("  (e) Radiation, with Chapter 5's figures (AP-8, field of 2026, solar minimum ; 10 % lost at 9e10 protons/cm2 above 10 MeV, 20 % at 4e11):")
P("      where                 | flux, worst longitude | 10 % / 20 % lost after  | flux, mean of the ring | 10 % / 20 % lost after | output lost after 30 years, mean of the ring [my interpolation]")
OUT["item4"]["radiation"] = {}
for nm, (fw, fm) in (("700 km, on the ring", FL["700"]), ("600 km, on the strip", FL["600"])):
    P("      %-21s | %5.0f /cm2/s          | %4.1f / %4.1f years       | %4.0f /cm2/s           | %4.0f / %4.0f years      | %.1f %%"
      % (nm, fw, 9e10/(fw*yr), 4e11/(fw*yr), fm, 9e10/(fm*yr), 4e11/(fm*yr), 100*(1 - loss(fm*yr*30))))
    OUT["item4"]["radiation"][nm] = dict(worst=fw, mean=fm, y10_worst=9e10/(fw*yr), y20_worst=4e11/(fw*yr), y10_mean=9e10/(fm*yr), y20_mean=4e11/(fm*yr), loss30_mean=1 - loss(fm*yr*30))
P("      The strip at 600 km takes %.2f times less at the worst longitude and %.2f times less on the mean than panels on the ring.  Chapter 5: 77 %% of the ring's circumference"
  % (FL["700"][0]/FL["600"][0], FL["700"][1]/FL["600"][1]))
P("      is untouched for 30 years at 700 km and over 85 % at 600 km.  A strip a few kilometers under the ring (item 3) gives this advantage up.")
OUT["item4"].update(demand_peak_GW=dmax["demand_GW"], demand_EJ=E_dem/1e18, supply_min_GW=min(c["supply_GW"] for c in tab4), supply_max_GW=max(c["supply_GW"] for c in tab4), ratio_min=ratio_min,
                    sheet_for_peak_m=need_m, g_day0=tab4[0]["g_strip"], g_final=g_fin, g_linear_700=3*MU/(RE + 700e3)**3*1e5, lean_max_rad=float(lean_max), drag_impulse_Ns=imp,
                    drag_share=imp/(12000*BOOKENV.Lr*(vc[0] - vc[-1])), drag_heat_EJ=e_drag/1e18, AO_fluence_face=flu, AO_rest600_per_year=f600*yr, AO_years_equiv=flu/(f600*yr),
                    AO_ram_fluence=flu_ram, days_above_1eV=max(c["day"] for c in tab4 if c["eV"] > 1.0), freed_kgm=41.5, freed_cable_kgm=41.5*CAB, push_day0_MN=c14["F"]/1e6*tab4[0]["g_strip"]/BOOKENV.g_a)

# ---------------------------------------------------------------------------------------------------------------------
# ITEM 5.  The swing with no arms
# ---------------------------------------------------------------------------------------------------------------------
head("ITEM 5.  The swing along the ring with no ground arms: seven and fourteen design launches a day, all one way")
M_CAS = 12000.0*BOOKENV.Lr; K_RING = 5800.0; K_LINE = 1.55; J_FULL = 0.576*M_CAS; J_NET = 0.281*M_CAS
F_LIM = 2*401.2*N_SITES; V_REL = 8846.96
om0 = math.sqrt(K_RING/M_CAS); C_CRIT = 2*math.sqrt(K_RING*M_CAS)
P("  Casing %.3e kg ; hold 5,800 N for each meter (Table 11.3, 'without the array') ; period %.1f h ; each reaction line %.2f N for each meter (round 14) and is sized for 93.6 kN."
  % (M_CAS, 2*math.pi/om0/3600, K_LINE))
P("  A design launch pushes the casing with %.2e N s ; %.2e N s is left once the sled is braked on the ring (0.576 and 0.281 m/s of the casing, the book's 0.58 and 0.28)." % (J_FULL, J_NET))
P("  The LIMs: Chapter 11 only says they can damp the swing by pushing the casing against the cable.  The push they have is taken from the deployment run's last day:")
P("  401 N from each of two LIMs at %d sites, %.1f MN in all, at 100 kV with the cable passing at %.0f m/s.  Used at full push that is %.0f GW." % (N_SITES, F_LIM/1e6, V_REL, F_LIM*V_REL/1e9))
P("  Two things are separate.  The MEAN push of the traffic is a steady force, launches x impulse / day.  Damping cannot remove it.  Either the lines carry it, and then")
P("  the casing sits displaced by force / 5,800 N/m, or the LIMs carry it against the cable, and then the cable's speed changes until something gives the momentum to the Earth.")
P("  The SWING is the motion about that place, and damping does remove it.")
def sim(N, J, c, hold=False, days=(8 if QUICK else 40), pulse=3690.0):
    """Casing on the spring, N launches a day, each a push J spread evenly over 'pulse' seconds.  LIM force - c v, limited to F_LIM ; with hold the LIMs also
    carry the mean push, as far as their limit allows.  Returns the largest displacement, the mean and the swing about it over the last 5 days."""
    dt = 20.0; T = 86400.0/N; n = int(days*86400/dt); x = 0.0; v = 0.0; xs_ = []; fl_ = []
    Fm = N*J/86400.0
    for i in range(n):
        t = i*dt; Fl = J/pulse if (t % T) < pulse else 0.0
        Fd = -c*v - (Fm if hold else 0.0)
        Fd = max(-F_LIM, min(F_LIM, Fd))
        a = (Fl + Fd - K_RING*x)/M_CAS
        v += a*dt; x += v*dt
        if t > (days - 5)*86400: xs_.append(x); fl_.append(abs(Fd))
    xs_ = np.array(xs_)
    return dict(x_max=float(np.abs(xs_).max()), x_mean=float(xs_.mean()), swing=float((xs_.max() - xs_.min())/2), F_lim_peak_MN=max(fl_)/1e6, F_lim_mean_MN=float(np.mean(fl_))/1e6)
P("  The launch is taken as a steady push over 3,690 s (the impulse over the 83.3 MN peak thrust) [ASSUMED shape].  Time step 20 s ; the last 5 of 40 days are read.")
P("\n    launches | sled     | mean push | LIMs                                   | casing sits at | swings about it by | farthest | reaction line there  | LIM push: peak, mean | LIM power at peak")
OUT["item5"] = dict(rows=[], F_LIM_MN=F_LIM/1e6, period_h=2*math.pi/om0/3600, c_crit=C_CRIT)
for N in (7, 14):
    for nm, J in (("braked", J_NET), ("not braked", J_FULL)):
        Fm = N*J/86400.0
        for cn, c, hold in (("none", 0.0, False), ("damping, critical (%.2e N s/m)" % C_CRIT, C_CRIT, False), ("damping and holding the mean push", C_CRIT, True)):
            r_ = sim(N, J, c, hold) if c else None
            if r_ is None:
                # undamped, from rest: the casing swings about the mean place with the mean displacement as amplitude, plus the ripple of each launch
                xm = Fm/K_RING; r_ = dict(x_max=2*xm + J/(M_CAS*om0), x_mean=xm, swing=xm + J/(M_CAS*om0), F_lim_peak_MN=0.0, F_lim_mean_MN=0.0); note = " (bound, from rest)"
            else: note = ""
            P("    %2d a day | %-8s | %5.1f MN  | %-38s | %7.0f m      | %7.0f m          | %6.0f m | %5.1f kN (%4.1f %%)%s | %5.1f, %5.1f MN       | %5.0f GW"
              % (N, nm, Fm/1e6, cn, r_["x_mean"], r_["swing"], r_["x_max"], K_LINE*r_["x_max"]/1e3, 100*K_LINE*r_["x_max"]/93.6e3, note, r_["F_lim_peak_MN"], r_["F_lim_mean_MN"], r_["F_lim_peak_MN"]*1e6*V_REL/1e9))
            OUT["item5"]["rows"].append(dict(N=N, sled=nm, mean_push_MN=Fm/1e6, lims=cn, **r_, line_kN=K_LINE*r_["x_max"]/1e3))
dv_cable = J_NET/(103569.71*BOOKENV.Lr)
P("\n  If the LIMs hold the mean push, the cable takes the momentum: %.4f m/s for each braked design launch, %.2f m/s a day at seven a day, %.0f m/s in a year, of its 8,331 m/s."
  % (dv_cable, 7*dv_cable, 7*dv_cable*365.25))
P("  (A launch to the east pushes the casing west ; the LIMs hold it by pushing the cable west, and the cable runs west, so it gains speed.)  The cable's speed, and with it the load on the bearing, drifts,")
P("  and the momentum has to go to the Earth in the end.  Only lines that lean along the ring do that.  With 5,800 N/m the mean push of seven braked launches a day")
P("  moves the casing %.1f km before the lines carry it ; a lean of %.1f km in 700 km is %.2f deg, and nothing in the lines' strength is near its limit."
  % (7*J_NET/86400/K_RING/1e3, 7*J_NET/86400/K_RING/1e3, math.degrees(7*J_NET/86400/K_RING/700e3)))
P("  What the swing is still not: a load on the mass driver (the stator moves with the casing) or a strain on the casing (it moves as a whole).")
OUT["item5"]["dv_cable_per_launch"] = dv_cable

with open(os.path.join(LOG_DIR, "round15_truss_array.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.ndarray): return [clean(v) for v in o.tolist()]
    return o
with open(os.path.join(LOG_DIR, "round15_truss_array.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
print("\nlog written: runs/round15_logs/round15_truss_array.txt and .json")
