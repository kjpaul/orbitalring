#!/usr/bin/env python3
"""
Round 26: two tables that do not agree with their neighbors.
Brief: Vol_III/Claude_Code_instructions_round26_2026-11-21.md.   Book: Orbital Ring Engineering, Volume III, Chapters 4 and 5.

1. Chapter 5.  Table 5.6 (1 m of water, worst longitude, 1.88 Sv/yr) against the last column of Table 5.11 (1 m of water all around,
   no cable, 0.72 Sv/yr at 56 W).  Both are rerun from their scripts (radiation/tables_trapped.py igrf ; radiation/sweetspot_trapped.py)
   and set against the files the chapter was keyed from.  The difference is then taken apart: the same spectrum and depth-dose response
   (radiation/fold_trapped.py), the shield's shape, and the direction the protons come from.
2. Chapter 4.  Table 4.11 with the side of the 250 km cable's structural section (ring_altitude.design(250)), by the method of the
   chapter's script (Vol_III/renumber_scripts_2026-10-12/ch04_numbers.py: a square that loses the depth from each side), and the masses
   behind Table 4.8's 250 km row (tools/ch4_numbers.py).

Run from Python/orbitalring:  python tools/round26_two_tables.py        (needs aep8, ppigrf, astropy ; about 2 minutes)
Log: runs/round26_logs/round26_two_tables.txt and .json ; the two reruns' own printing beside them.
"""
import io, json, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); RAD = os.path.join(ROOT, "radiation")
LOG_DIR = os.path.join(ROOT, "runs", "round26_logs"); os.makedirs(LOG_DIR, exist_ok=True)
T0 = time.time(); _buf = io.StringIO()
def P(s=""):
    print(s); _buf.write(s + "\n"); sys.stdout.flush()
OUT = {}

# =====================================================================================================================
P("Round 26.  Two tables that do not agree with their neighbors.")
P("\n" + "="*118 + "\n1. Chapter 5: Table 5.6 against the last column of Table 5.11\n" + "="*118)
def rerun(script, args, kept):
    o = subprocess.run([sys.executable, script] + args, cwd=RAD, capture_output=True, text=True, timeout=1800)
    txt = o.stdout.replace("\r\n", "\n"); old = open(os.path.join(RAD, kept), encoding="utf-8").read().replace("\r\n", "\n")
    open(os.path.join(LOG_DIR, "rerun_" + kept), "w", encoding="utf-8").write(txt)
    same = txt.strip() == old.strip()
    P("   python %s %s: %d lines ; the same as radiation/%s, character for character: %s" % (script, " ".join(args), len(txt.splitlines()), kept, same))
    if not same and o.stderr: P("   stderr: " + o.stderr[-400:])
    return txt, same
t56, same56 = rerun("tables_trapped.py", ["igrf"], "tables_trapped_igrf2026.txt")
t511, same511 = rerun("sweetspot_trapped.py", [], "sweetspot_trapped.txt")
OUT["reruns"] = dict(table_5_6_same=same56, table_5_11_same=same511)
if not (same56 and same511):
    P("   STOP: a rerun does not give the kept file."); sys.exit(1)
def block(txt, head):
    ls = txt.split("\n"); i = next(k for k, l in enumerate(ls) if l.startswith(head)); return ls[i:]
b56 = block(t56, "700 km, solar min, spectrum above 300 MeV: exponential")
r100 = next(l for l in b56 if l.strip().startswith("100 ")); P("   Table 5.6's row (700 km, IGRF 2026, spectrum continued above 300 MeV): " + " ".join(r100.split()))
b511 = block(t511, "Ring at 700 km"); r56 = next(l for l in b511 if l.strip().startswith("-56 ")); P("   Table 5.11's row at 56 W:   " + " ".join(r56.split()))

os.chdir(RAD); sys.path.insert(0, RAD); sys.path.insert(0, ROOT)
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from sweetspot import U, thickness, W as CABLE_W
from sweetspot_trapped import weights, d as SP, Eg, lons
from fold_trapped import dose, T
ALT = 700; N = len(U); ux, uy, uz = U.T
def Hcurve(k):
    D, H = dose(Eg, np.nan_to_num(SP["min_%d" % ALT][:, k]), "exp"); return H
def fold(k, t, w):
    """sweetspot_trapped.point_dose's sum, with the weights given: the sphere-center dose equivalent at the water path of each direction."""
    H = Hcurve(k); lt = np.log(np.maximum(H, 1e-12)); Ht = np.exp(np.interp(np.clip(t, T[0], None), T, lt, right=-30))
    sl = (lt[-1] - lt[-3])/(T[-1] - T[-3]); Ht = np.where(t > T[-1], np.exp(lt[-1] + sl*(t - T[-1])), Ht)
    return float((w*Ht).sum()), Ht
def room_tube(tw=100.0):
    """Table 5.11's last column, as sweetspot_trapped.py builds it: the room is a tube along the ring, 4 m wide, flat slabs of tw for the ceiling 1 m above
    the point, the floor 1 m below and the two side walls ; no end walls."""
    t, _ = thickness(1.0, t_wall=tw, t_floor=tw, cable=False)
    up = (uz > 0) & (np.abs(1.0*uy/np.where(uz > 0, uz, 1)) < 2.0)
    return np.where(up, tw/np.maximum(uz, 1e-9), t)
def room_box(half_len, tw=100.0):
    """The same room closed with two end walls half_len m east and west of the point: each direction leaves by the first face it meets."""
    big = 1e30
    dx = np.where(np.abs(ux) > 1e-12, half_len/np.maximum(np.abs(ux), 1e-12), big); dy = np.where(np.abs(uy) > 1e-12, 2.0/np.maximum(np.abs(uy), 1e-12), big); dz = np.where(np.abs(uz) > 1e-12, 1.0/np.maximum(np.abs(uz), 1e-12), big)
    first = np.argmin(np.stack([dx, dy, dz]), axis=0); comp = np.choose(first, [np.abs(ux), np.abs(uy), np.abs(uz)])
    return tw/np.maximum(comp, 1e-9)
iso = np.full(N, 1.0/N); t_tube = room_tube(); t_sph = np.full(N, 100.0)
P("\n   The two tables use the same AP-8 spectra (ap8_igrf2026_spectra.npz, solar minimum, continued above 300 MeV) and the same Geant4 depth-dose")
P("   response (fold_trapped.dose: the dose equivalent at the center of a water sphere against its radius).  They differ in the shield's shape and in the direction of the protons.")
P("   Table 5.6: a solid sphere, 100 g/cm2 of water in every direction, flux the same from every direction.")
P("   Table 5.11, last column: a room that is a tube along the ring, 4 m wide, the point 1 m under the ceiling and 1 m above the floor, flat slabs of 100 g/cm2 (ceiling, floor,")
P("   two side walls, no end walls), so a direction that crosses a slab on a slant goes through 100 / cos of water ; the protons come from the directions AP-8 gives at that place.")
P("\n   longitude | sphere, 100 everywhere (Table 5.6) | the tube room, flux from every direction alike | the tube room, AP-8's directions (Table 5.11) | factors: slant, direction, both | closed room, end walls 2 m away | 5 m | 10 m | 50 m")
dec = {}
for lon in (-60, -56, -50, -45, -40, -30):
    k = int(np.where(lons == lon)[0][0]); w, bh = weights(ALT, k)
    a, _ = fold(k, t_sph, iso); a2, _ = fold(k, t_sph, w); b, _ = fold(k, t_tube, iso); c, Ht = fold(k, t_tube, w)
    box = {hl: fold(k, room_box(hl), w)[0] for hl in (2.0, 5.0, 10.0, 50.0)}; box_iso = fold(k, room_box(2.0), iso)[0]
    assert abs(a - a2) < 1e-9*a
    # the water path, weighted by the flux
    order = np.argsort(t_tube); cw = np.cumsum(w[order]); med = float(t_tube[order][np.searchsorted(cw, 0.5)]); cwi = np.cumsum(iso[order]); medi = float(t_tube[order][np.searchsorted(cwi, 0.5)])
    sh = {x: float(w[t_tube > x].sum()) for x in (110, 150, 200, 300)}; shi = {x: float(iso[t_tube > x].sum()) for x in (110, 150, 200, 300)}
    dz_ = float((w*Ht)[t_tube <= 150].sum()/c)
    dec[lon] = dict(sphere=a, tube_iso=b, tube_dir=c, f_slant=a/b, f_dir=b/c, f_both=a/c, box=box, box2_iso=box_iso, median_path_dir=med, median_path_iso=medi, share_flux_path_over=sh, share_flux_path_over_iso=shi, share_dose_path_under_150=dz_,
                    field_enu=[float(v) for v in bh], below_horizontal=float(w[uz < 0].sum()))
    P("   %5d     | %.4f Sv/yr                       | %.4f Sv/yr                                    | %.4f Sv/yr                                  | %.2f x %.2f = %.2f              | %.4f                          | %.4f | %.4f | %.4f"
      % (lon, a, b, c, a/b, b/c, a/c, box[2.0], box[5.0], box[10.0], box[50.0]))
q = dec[-56]; q45 = dec[-45]
assert abs(q["tube_dir"] - 0.7210) < 5e-4
# Table 5.6's "worst longitude" is the largest over all longitudes at each thickness (tables_trapped.py: H[:, i].max()) ; where is it?
Hall = np.array([Hcurve(k) for k in range(len(lons))]); worst = {}
P("\n   Table 5.6's column 'worst longitude' is the largest value over all longitudes, taken again at each thickness (tables_trapped.py: H[:, i].max()).  Where that is:")
P("   shield g/cm2 | largest, Sv/yr | at longitude | at 56 W, the flux peak | largest / 56 W")
i56 = int(np.where(lons == -56)[0][0])
for sgc in (1, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300):
    i = int(np.where(T == sgc)[0][0]); kx = int(Hall[:, i].argmax()); worst[sgc] = dict(H=float(Hall[kx, i]), lon=float(lons[kx]), H56=float(Hall[i56, i]))
    P("   %6d       | %9.4f      | %5.0f        | %9.4f              | %.3f" % (sgc, Hall[kx, i], lons[kx], Hall[i56, i], Hall[kx, i]/Hall[i56, i]))
assert abs(worst[100]["H"] - 1.8786) < 5e-4
OUT["worst_by_thickness"] = {str(k): v for k, v in worst.items()}
P("   So Table 5.6's 1.8786 Sv/yr for 100 cm is at %.0f W, where the spectrum is harder, and not at 56 W.  At 56 W the same sphere gives %.4f Sv/yr." % (-worst[100]["lon"], q["sphere"]))
P("   The factor between the two printed numbers, 1.8786 / 0.7210 = %.2f, is: the longitude %.3f (45 W against 56 W in the sphere) x the slant through flat slabs %.3f x the protons' directions %.3f."
  % (1.8786/0.7210, worst[100]["H"]/q["sphere"], q["f_slant"], q["f_dir"]))
P("   At 45 W itself: the sphere %.4f Sv/yr, the tube room %.4f Sv/yr (not in Table 5.11, which prints 56 W and 40 W), a factor %.2f." % (q45["sphere"], q45["tube_dir"], q45["f_both"]))
P("\n   At 56 W.  The water path in the tube room: at least 100 g/cm2 in every direction ; half of the flux goes through more than %.0f g/cm2 if it came from every direction alike, and through more than %.0f g/cm2 on AP-8's directions." % (q["median_path_iso"], q["median_path_dir"]))
P("   Share of the flux that goes through more than 110, 150, 200, 300 g/cm2: every direction alike %s %% ; AP-8's directions %s %%."
  % (", ".join("%.0f" % (100*q["share_flux_path_over_iso"][x]) for x in (110, 150, 200, 300)), ", ".join("%.0f" % (100*q["share_flux_path_over"][x]) for x in (110, 150, 200, 300))))
P("   %.0f %% of the tube room's dose comes in through the directions with less than 150 g/cm2." % (100*q["share_dose_path_under_150"]))
P("   The field there (east, north, up): %s.  The protons move within a few degrees of the plane across the field, a nearly vertical plane close to the plane of the ring's own length." % np.round(q["field_enu"], 3))
P("   The directions shift where the flux meets the slabs, but the dose comes out within %.0f %% of the dose for a flux from all sides (the factor %.3f).  The slant is the whole of it." % (100*abs(q["f_dir"] - 1), q["f_dir"]))
P("   The closed room with end walls is not in the chapter.  It is run here to show how much of the 0.72 depends on the room having no ends: end walls 2 m away (a room 4 m long) give %.2f Sv/yr at 56 W, 50 m away %.2f, no end walls %.2f."
  % (q["box"][2.0], q["box"][50.0], q["tube_dir"]))
km = max(dec, key=lambda l: dec[l]["tube_dir"]); ks = max(dec, key=lambda l: dec[l]["sphere"])
tube_all = np.array([fold(k, t_tube, weights(ALT, k)[0])[0] if SP["min_%d" % ALT][0, k] > 0 and weights(ALT, k) is not None else 0.0 for k in range(len(lons))]); kt = int(tube_all.argmax())
P("   The tube room over every longitude: largest %.4f Sv/yr at %.0f W.  Table 5.11 prints 0.72 at 56 W and 0.80 at 40 W, and does not print that longitude." % (tube_all[kt], -lons[kt]))
OUT["tube_worst"] = dict(H=float(tube_all[kt]), lon=float(lons[kt]))
OUT["decomposition"] = {str(k): v for k, v in dec.items()}; OUT["cable_side_in_sweetspot"] = float(CABLE_W)
P("  (%.0f s so far)" % (time.time() - T0))

# =====================================================================================================================
os.chdir(ROOT)
P("\n" + "="*118 + "\n2. Chapter 4: Table 4.11 and the 250 km row of Table 4.8\n" + "="*118)
import ring_altitude as R
d2 = R.design(250, True); L = d2["L_ring"]; ms = d2["m_struct"]; mc = d2["m_cable"]; side = d2["side"]; side_hw = math.sqrt(mc/R.RHO_CNT)
M_s = ms*L/1e9; M_c = mc*L/1e9
P("   ring_altitude.design(250, retrograde), quadratic '%s', hardware %.0f kg/m: structural CNT %.1f kg/m, with the hardware %.1f kg/m ; ring %.1f km ; structural cable %.1f Mt, with the hardware %.1f Mt."
  % (R.QUADRATIC, R.M_HARDWARE_M, ms, mc, L/1e3, M_s, M_c))
P("   Side of the structural section sqrt(%.1f / %.0f) = %.3f m (design()['side']).  Side of a square with the hardware's mass counted as CNT: %.3f m, the book's 7.94." % (ms, R.RHO_CNT, side, side_hw))
P("   The chapter's script for this table is Vol_III/renumber_scripts_2026-10-12/ch04_numbers.py (its last lines ; output kept in Vol_III/runs_700km_hw10686_2026-10-12/ch04_numbers_extra.txt).")
P("   Its method: a square of side design(250)['side'] that loses the depth from each side, share = 1 - ((side - 2 depth) / side)^2.  No corners, no round section.  The side it takes is the structural one.")
def share(s, dep): return 1 - ((s - 2*dep)/s)**2
old_side = 8.43; kept = dict(d10=0.046876543394601455, d50=0.22313034655682518, d100=0.41812976707319327)
s_back = 2*0.5/(1 - math.sqrt(1 - kept["d50"]))
P("   Its kept output (October, the fixed quadratic): 10 cm %.4f, 50 cm %.4f, 100 cm %.4f.  Those are a square of side %.3f m, the structural section of that time (the chapter then printed 8.43 m)." % (kept["d10"], kept["d50"], kept["d100"], s_back))
P("   So the printed 5 %, 22 % and 40 % are the October shares (4.7, 22.3, 41.8 %) rounded, and they were not run again when the cable became lighter.")
rows = {}
P("\n   depth from each side | share ablated, side %.3f m (structural, the script's) | side %.3f m (with the hardware) | side %.2f m (October)" % (side, side_hw, old_side))
for nm, dep in (("10 cm", 0.1), ("50 cm", 0.5), ("100 cm", 1.0)):
    rows[nm] = dict(structural=share(side, dep), with_hw=share(side_hw, dep), october=share(s_back, dep))
    P("   %-20s | %.4f  (%.1f %%)                                    | %.4f (%.1f %%)                  | %.4f (%.1f %%)" % (nm, rows[nm]["structural"], 100*rows[nm]["structural"], rows[nm]["with_hw"], 100*rows[nm]["with_hw"], rows[nm]["october"], 100*rows[nm]["october"]))
P("\n   The soot row, by the script's own line: structural cable x the plausible share x the soot fraction.  The script takes 5 %, 15 % and 25 % for the plausible shares and 30 %, 50 % and 70 % for the soot fractions.")
soot = dict(low=M_s*0.05*0.3, central=M_s*0.15*0.5, high=M_s*0.25*0.7, central_10=M_s*0.10*0.5)
P("   low %.1f Mt ; central %.1f Mt (at 10 %% in place of 15 %%: %.1f Mt) ; high %.1f Mt.  The book prints 60, 300 and 720 Mt." % (soot["low"], soot["central"], soot["central_10"], soot["high"]))
P("   The 720 is not the script's.  It is October's printed 900 Mt (the script then gave %.0f Mt, rounded up) times 0.80 (final_pass_scripts_2026-10-31/ch04_sheet_250km_2026-11-07.py).  The script's line gives %.0f Mt." % (5033*0.25*0.7, soot["high"]))
P("   The plausible shares (5 %%, 10 to 15 %%, 25 %%) are the script's constants.  They are a judgment and are not computed from the ablated shares: the low one equals the 10 cm share (%.1f %%), the central is %.0f to %.0f %% of the 50 cm share, the high is %.0f %% of the 100 cm share."
  % (100*rows["10 cm"]["structural"], 100*0.10/rows["50 cm"]["structural"], 100*0.15/rows["50 cm"]["structural"], 100*0.25/rows["100 cm"]["structural"]))
P("   The two loading rows have no script.  The table's note gives the rule (deposition rate x residence time), and the printed rows are October's 13 to 23 and 44 to 76 Mt times 0.80.  By the rule, on the central %.1f Mt:" % soot["central"])
load = dict(century=(soot["central"]/100*3.5, soot["central"]/100*6.0), thirty=(soot["central"]/30*3.5, soot["central"]/30*6.0), century_300=(300/100*3.5, 300/100*6.0), thirty_300=(300/30*3.5, 300/30*6.0))
P("   over a century, residence 3.5 to 6 years: %.1f to %.1f Mt (on the rounded 300 Mt: %.1f to %.1f) ; within thirty years: %.1f to %.1f Mt (on 300 Mt: %.1f to %.1f).  The book prints 10 to 18 and 35 to 61."
  % (load["century"] + load["century_300"] + load["thirty"] + load["thirty_300"]))
P("   Of the 50 Mt threshold the loading rows are %.2f to %.2f (the chapter: 'about a fifth of the threshold to 1.2 times it')." % (load["century_300"][0]/50, load["thirty_300"][1]/50))
P("\n   The structural cable is %.1f Mt here.  The book prints 4,019 Mt.  With the hardware it is %.1f Mt (the book: 4,466), and %.1f less the hardware's %.1f Mt is %.1f.  The 4,019 is not reproduced ; the present code gives %.0f."
  % (M_s, M_c, M_c, R.M_HARDWARE_M*L/1e9, M_c - R.M_HARDWARE_M*L/1e9, M_s))
P("\n   Table 4.8's 250 km row.  tools/ch4_numbers.py (lines 59 to 62) takes the ejecta as 5 % of the cable with its hardware (m_cable_total) and the soot as 30 % of that:")
ej = 0.05*M_c
P("   0.05 x %.1f Mt = %.1f Mt of ejecta, x 0.30 = %.1f Mt.  On the structural CNT alone it would be 0.05 x %.1f x 0.30 = %.1f Mt.  At the shock-layer 0.3 %%: %.2f Mt with the hardware, %.2f Mt structural." % (M_c, ej, 0.3*ej, M_s, 0.3*0.05*M_s, 0.3*0.003*M_c, 0.3*0.003*M_s))
P("   So the two tables do stand on different masses, and it is what the two scripts do: the ejecta are a share of everything in the cable that is hit, the hardware with it ; the slow soot is counted on the carbon alone.")
P("   Neither script says why.  Whether the hardware's mass should make soot at the 30 % particulate rate is a choice for the author ; the scripts do not settle it.")
OUT["ch4"] = dict(m_struct=ms, m_cable=mc, L_ring=L, M_struct_Mt=M_s, M_cable_Mt=M_c, side_struct=side, side_with_hw=side_hw, october_side=s_back, october_kept=kept, shares=rows, soot=soot, loading=load,
                  table_4_8=dict(ejecta_Mt=ej, soot_with_hw=0.3*ej, soot_struct=0.3*0.05*M_s, shock_with_hw=0.3*0.003*M_c, shock_struct=0.3*0.003*M_s))

P("\n(%.0f s)" % (time.time() - T0))
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    try:
        import numpy as _np
        if isinstance(o, (_np.floating, _np.integer, _np.bool_)): return o.item()
    except Exception: pass
    return o
with open(os.path.join(LOG_DIR, "round26_two_tables.txt"), "w", encoding="utf-8") as fh: fh.write(_buf.getvalue())
with open(os.path.join(LOG_DIR, "round26_two_tables.json"), "w", encoding="utf-8") as fh: json.dump(clean(OUT), fh, indent=1)
print("log written: runs/round26_logs/round26_two_tables.txt and .json")
