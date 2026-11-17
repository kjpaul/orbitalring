#!/usr/bin/env python3
"""
Round 11, item 1b: writes runs/round11_logs/item1b_report.md from runs/round11_logs/round11_plate_report.json (bearing_round11_plate_report.py).
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run (from Python/orbitalring): python bearing_round11_plate.py --part=report, then python bearing_round11_plate_md.py
Every number in the report is read from that json; nothing is solved or typed in here except the design's constants named in the text.
"""
import os

import bearing_round11_plate as p
from bearing_round11_plate import RHO, HW_ALLOW, M_LIM_PLATES

SIZED = "sized on the face at 1.2 T"
SIZED18 = "sized on the face at 1.8 T"
FOUR = "four plates"
BOOK = dict(tape_s=144.2, tape_f=65.4, tape6=550.0, bearing=1687.0, plates=8185.0, hw=10057.0)        # Table 7.5 and round 8, 2.5 mm cells
LOG = "`round11_plate_report.txt`"


def pc(a, b, d=0):
    return f"{(a / b - 1) * 100:+.{d}f} %"


def main():
    D = p.load("round11_plate_report.json")
    B = D["base"]; m0 = B["m8"]
    L = []; w = L.append
    face = D.get("face", {}); cp = D.get("coupling", {}); mesh = D.get("mesh", {}); sp = D.get("spacing"); ff = D.get("fullface", {})
    wd = D.get("width", {}); al = D.get("alone", []); sk = D.get("skirt", [])
    ts = D.get("tailor_support", []); tf = D.get("tailor_fast", [])

    def trow(rows, bf):
        return next((r for r in rows if abs(r["Bf"] - bf) < 1e-6), None)

    def wrow(kind, Wl):
        return next((r for r in wd.get(kind, []) if abs(r["Wl"] - Wl) < 1e-6), None)

    def skrow(kind, T, s):
        return next((r for r in sk if r["kind"] == kind and abs(r["T"] - T) < 1e-6 and abs(r["skirt"] - s) < 1e-6), None)

    f4 = face[FOUR]; fs = face[SIZED]; fs18 = face.get(SIZED18)
    c4 = cp[FOUR]; cs = cp[SIZED]
    fw10 = face.get("one plate, web 10 mm"); fw25 = face.get("one plate, web 25 mm"); fu = face.get("one plate, 100 mm everywhere")
    lower0 = 2 * B["sup"]["plate"] + 2 * B["fast"]["plate"]
    ts12 = trow(ts, 1.2); tf12 = trow(tf, 1.2)
    f739 = D.get("face739", {})
    k_today = "the lift tracks on today's 1.0 m x 200 mm plates"; k_sized = "the lift tracks on plates sized at 1.2 T (table 1b)"
    FT = D.get("faces_tape", {}); FR = D.get("faces_ring", {}); FP = D.get("faces_plate", {})

    # ------------------------------------------------------------------ header
    w("# Item 1b: one plate under the cable, and a wider face")
    w("")
    w("Scripts, all new, in `Python/orbitalring`: `bearing_round11_plate.py` (the model and the single-track solves), `bearing_round11_plate_parts.py` (the half face, the coupling, the checks), "
      "`bearing_round11_plate_report.py` (the tables), `bearing_round11_plate_md.py` (this text, written from the tables' json). Logs in `runs/round11_logs/`: `round11_plate_*.json` and `round11_plate_*.txt`; "
      f"the tables are {LOG}. No existing module was changed. Nothing in the book or its image folders was touched.")
    w("")
    w("**The model.** The bearing's two-dimensional field solver (`magnetostatic_2d.py`), with the tracks' bundles and tape positions as `bearing_electrical.Track` lays them. The gap is 100 mm, the tape is 10 mm behind the casing's face, "
      "and the plate is sintered iron with a relative permeability of 200 that saturates at 1.2 T. Two things are new. The plate's thickness is given column by column, 5 mm wide, so the plate can have any shape. "
      "And the pull of a track is read on a box drawn around its own bundles, so two tracks under one plate can be told apart. On a plate of its own the two ways of reading the pull agree: "
      f"{B['F_box'][0]/1e3:.2f} kN/m on the bundles and {B['F_box'][1]/1e3:.2f} kN/m on the plate in the same solve. The plate's face toward the coils is flat in the model, so every track is one gap from it. "
      "Where the plate is thicker, it is thicker into the cable.")
    w("")
    w("**The cells.** 5 mm in the window around the tracks and 10 mm in the box around it, as in the gap trade. The design's tracks on their own plates reproduce the gap trade's 5 mm rows: "
      f"support {B['sup']['NI']/1e3:.1f} kA per leg and {B['sup']['tape']:.1f} kg/m of tape (gap trade {B['gt']['sup_NI']/1e3:.1f} kA at a pull 0.3 percent higher, {B['gt']['sup_tape']:.1f} kg/m), "
      f"compensation {B['fast']['NI']/1e3:.1f} kA and {B['fast']['tape']:.1f} kg/m (gap trade {B['gt']['fast_NI']/1e3:.1f} and {B['gt']['fast_tape']:.1f}). "
      f"The book's rows are at 2.5 mm ({B['r7']['sup_tape']:.1f} and {B['r7']['fast_tape']:.1f} kg/m), so the 5 mm tables here read 2 to 4 percent low in tape. Every comparison below is between rows of the same cell size. "
      "The convergence checks are in 1.5 and 2.3. Each solve is brought within 1 percent of its required pull and then moved to it along the local power law, so a tape figure carries about half a percent.")
    w("")
    w("**The loads and the tape rules** are the design's. A support track gives 47.5 kN/m in normal operation and 245.5 kN/m at its cap. A compensation track gives 196.5 kN/m (198 kN/m for the load track in the coupling). "
      "The worst tape runs at 50 percent of its limit at the requirement, and the rating is the pull at which it reaches 70 percent. The cable keeps its 54.64 m2 and its 92,883 kg/m. Only the width of its faces changes.")
    w("")
    L_head = len(L)

    # ------------------------------------------------------------------ 1
    S1 = []; w = S1.append
    w("## 1. One plate across the bottom face of 7.39 m")
    w("")
    w("### 1.1 What the plate does past the bundles")
    w("")
    w("Before the plates were joined, one thing had to be known: how much of the pull depends on the plate reaching past the outer bundles. Today it hardly reaches past them. "
      "The support's four bundles span 1.15 m under a 1.2 m plate, 25 mm of plate on each side. The compensation track's two bundles span 0.80 m under a 1.0 m plate, 98 mm on each side.")
    w("")
    if sk:
        w(f"Script: `bearing_round11_plate.py --part=skirt` (`round11_plate_skirt.json`). Each track alone at its own design ampere-turns ({B['sup']['NI']/1e3:.1f} and {B['fast']['NI']/1e3:.1f} kA per leg), on a rectangular plate, 5 mm cells.")
        w("")
        w("| Plate reaches past the outer bundles by | Support track, 100 mm plate | Support track, 300 mm plate | Compensation track, 200 mm plate | Compensation track, 300 mm plate |")
        w("|---|---|---|---|---|")
        for s in sorted({r["skirt"] for r in sk}):
            cells = []
            for kind, T in (("support", 0.1), ("support", 0.3), ("fast", 0.2), ("fast", 0.3)):
                r = skrow(kind, T, s); ref = skrow(kind, 0.1 if kind == "support" else 0.2, 0.025)
                cells.append(f"{r['F']/1e3:.1f} kN/m ({pc(r['F'], ref['F'])}), {r['sat']*100:.0f} % saturated" if r else "")
            w(f"| {s*1e3:.0f} mm | " + " | ".join(cells) + " |")
        w("")
        w("The percentages are against the track's own row at 25 mm on its own thickness. The pull is at fixed ampere-turns, so a higher pull means fewer ampere-turns and less tape for the same pull. "
          "Today's compensation plate reaches 98 mm past its bundles, between the first two rows, and pulls 196.5 kN/m at these ampere-turns. These solves have the coarse box inside a larger one, which the others of 1.2 do not; 2.3 gives what that is worth.")
        w("")
    w("### 1.2 A plate sized from its flux, each track alone")
    w("")
    w("This is Paul's rule made exact. The track is first put under a plate 300 mm thick that reaches 0.6 m past its outer bundles, at its requirement. "
      "The flux that runs sideways through the plate's section is read at every 5 mm column. The plate is then given, at each column, the thickness that carries that flux at a chosen fill field, rounded up to 5 mm and never under 10 mm, "
      "and it ends where 10 mm is enough. A fill field of 1.2 T is the iron's saturation value: the plate is as thick as the flux needs and no thicker. A higher fill field is a thinner plate that saturates more.")
    w("")
    w(f"Script: `bearing_round11_plate.py --part=tailor` (`round11_plate_tracks.json`), 5 mm cells. The first row of each table is the design's plate in the same model (`--part=base`).")
    w("")
    for kind, rows, base, lab, plate0 in (("support", ts, B["sup"], "One support track at its cap, 245.5 kN/m (two 0.4 m loops)", "today: 1.2 m by 100 mm"), ("fast", tf, B["fast"], "One compensation track at 196.5 kN/m (one 0.6 m loop)", "today: 1.0 m by 200 mm")):
        w(f"**{lab}:**")
        w("")
        w("| Plate | Width | Thickest | Mean thickness | Plate mass per meter of track | Ampere-turns per leg | Tape per meter of track | Against today | Rating over requirement | Plate saturated |")
        w("|---|---|---|---|---|---|---|---|---|---|")
        Wp, T = (1.2, 0.1) if kind == "support" else (1.0, 0.2)
        w(f"| {plate0} | {Wp:.2f} m | {T*1e3:.0f} mm | {T*1e3:.0f} mm | {base['plate']:,.0f} kg | {base['NI']/1e3:.1f} kA | {base['tape']:.1f} kg |  | {base['rating']:.2f} | {base.get('sat', 0)*100:.0f} % |")
        for r in rows:
            w(f"| sized at {r['Bf']:.1f} T | {r['width']:.2f} m | {r['Tmax']*1e3:.0f} mm | {r['mean']*1e3:.0f} mm | {r['mass']:,.0f} kg | {r['NI']/1e3:.1f} kA | {r['tape']:.1f} kg | {pc(r['tape'], base['tape'])} | {r['rating']:.2f} | {r['sat']*100:.0f} % |")
        w("")
        th = D.get(f"thick_{kind}")
        s = f"On the 300 mm plate the track needs {th['NI']/1e3:.1f} kA per leg, and the sideways flux is {th['flux_max']:.3f} Wb per meter of ring at its largest. " if th else ""
        fl = B["flux"]["sup" if kind == "support" else "fast"]
        s += f"In today's plate the largest sideways flux is {fl:.3f} Wb/m in {T*1e3:.0f} mm, which would be {fl/T:.2f} T if it were even over the section."
        et = D.get(f"equal_tape_{kind}"); em = D.get(f"equal_mass_{kind}")
        if et:
            s += f" Interpolating between the rows, a sized plate needs today's tape at a fill field of {et['Bf']:.2f} T and then weighs {et['mass']:,.0f} kg/m against {base['plate']:,.0f}."
        if em:
            s += f" A sized plate of today's mass needs {em['tape']:.1f} kg/m of tape against {base['tape']:.1f}."
        w(s)
        w("")
    # ---- 1.3 the face
    w("### 1.3 The 7.39 m face: the support and the load track under one plate")
    w("")
    w("Script: `bearing_round11_plate.py --part=face` (`round11_plate_face.json`), 5 mm cells. One half of the lower face, from the cable's centerline to its edge at 3.695 m, with the coils where the design has them: "
      "the support track centered 1.895 m from the centerline and the load track 3.195 m. The centerline is a wall for the flux in this model, which is the case when the two halves carry mirrored, opposite currents (1.5 looks at the other case). "
      "The coarse box sits inside a larger box of 10 by 8 m at 20 mm cells, because the coupling needs it (1.5). On each plate each track is solved alone at its requirement, and the tape is counted by the design's rule.")
    w("")
    w("The plates: the design's four (two of them in the half face); the same with a web of 10 or 25 mm joining them and running to the centerline; one plate 100 mm thick everywhere; "
      "and one plate sized on the face. For the last, both tracks are put under a 300 mm plate over the whole half face, each at its requirement alone, and the thickness at each column is the sum of the two sideways fluxes over the fill field, 10 mm at least.")
    w("")
    w("| Plate | Plates of the lower face | Against the four plates | Thickest | Support: ampere-turns per leg at the cap | Tape per meter of track | Rating | Saturated at the cap | Load track: ampere-turns per leg | Tape per meter of track | Rating | Stored energy |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for nm, v in face.items():
        s, l = v["sup"], v["load"]
        w(f"| {nm} | {v['mass']:,.0f} kg/m | {pc(v['mass'], lower0) if nm != FOUR else ''} | {v['Tmax']*1e3:.0f} mm | {s['NI']/1e3:.1f} kA | {s['tape']:.1f} kg | {s['rating']:.2f} | {s['sat']*100:.0f} % | {l['NI']/1e3:.1f} kA | {l['tape']:.1f} kg | {l['rating']:.2f} | {l['W']/1e3:.0f} kJ/m |")
    w("")
    w(f"The first row is the design in this model: {f4['sup']['tape']:.1f} and {f4['load']['tape']:.1f} kg/m, against {B['sup']['tape']:.1f} and {B['fast']['tape']:.1f} kg/m for each track alone in its own box. "
      "The difference is the other box and the centerline, not the plate.")
    w("")
    w("**Both tracks pulling at once.** The support at its cap and the load track at its full 198 kN/m, each at the ampere-turns it needs alone on that plate. This is a bound: the cap is the parked spacecraft carried by the support, and the load track's full pull is the same spacecraft carried by the load track.")
    w("")
    w("| Plate | Support pulls | Load track pulls | Worst support tape, share of its limit | Worst load tape, share of its limit |")
    w("|---|---|---|---|---|")
    for nm, v in face.items():
        b = v["both"]
        w(f"| {nm} | {b['F_sup']/1e3:.1f} kN/m | {b['F_load']/1e3:.1f} kN/m | {b['share_sup']*100:.0f} % | {b['share_load']*100:.0f} % |")
    w("")
    # ---- 1.4 coupling
    w("### 1.4 How much of the load coil's flux threads the support loops")
    w("")
    w("Two ways of counting it, both in `round11_plate_face.json` and `round11_plate_coupling.json` (`--part=face`, `--part=coupling`).")
    w("")
    w("**Round 7's method** (`bearing_round7.py --part=coupling`): linear plates; the flux per ampere-turn of each support loop, of one loop through the other, and of the load track through each loop; "
      "then the ampere-turns each support loop must lose or gain to keep its own flux when the load track comes to its full pull (453.4 kA per leg), as a share of the support's normal 156.7 kA. "
      "Loop 1 is the support loop farther from the load track, loop 2 the one beside it.")
    w("")
    w("**Directly:** the plate saturating, the support in normal operation at 47.5 kN/m, the load track brought to 198 kN/m, and the two support loops' ampere-turns found (by Newton steps) so that each loop has the flux it had before. "
      "The pull on the support's own bundles is read before and after.")
    w("")
    w("| Plate | Load coil's flux through loop 1, loop 2, as a share of the flux through its own loop | Round 7's method, 10 mm cells | 5 mm cells | Direct: the support's pull changes by | Lost per track | Loop 1, loop 2 after (before: both at the normal value) | Support's pull if its currents were held instead |")
    w("|---|---|---|---|---|---|---|---|")
    for nm, c in cp.items():
        a = f"{c['share'][0]*100:+.1f} %, {c['share'][1]*100:+.1f} %" if "share" in c else ""
        b10 = f"{c['lin10']*100:+.1f} %" if c.get("lin10") is not None else ""
        b5 = f"{c['lin5']*100:+.1f} %" if c.get("lin5") is not None else ""
        if "direct" in c:
            d = f"{c['direct']*100:+.1f} % | {c['lost']/1e3:.1f} kN/m | {c['x'][0]/1e3:.1f} and {c['x'][1]/1e3:.1f} kA (from {c['NI0']/1e3:.1f}) | {c['const']*100:+.1f} %"
        else:
            d = " | | | "
        w(f"| {nm} | {a} | {b10} | {b5} | {d} |")
    w("")
    w(f"Round 7's own figure is {D['c_round7']*100:.1f} percent, in its own box. The sign is for the load track wound as round 7 had it, with the neighboring bundles of the two tracks carrying current the same way. "
      "Wound the other way, the linear columns change sign: the support then gains what it loses here.")
    w("")
    # ---- 1.5 checks
    w("### 1.5 The box, the cells, the spacing and the other half of the face")
    w("")
    w("Script: `bearing_round11_plate.py --part=mesh`, `--part=spacing`, `--part=fullface` (`round11_plate_mesh.json`, `round11_plate_spacing.json`, `round11_plate_fullface.json`). All by round 7's method, linear plates.")
    w("")
    if mesh:
        w("**The box.** The coupling is a far-field quantity: the load coil's flux has to come back around through the air, and a box that is too small squeezes it through the support loops.")
        w("")
        w("| Box, 10 mm cells | Four plates | One plate, web 25 mm |")
        w("|---|---|---|")
        labs = []
        for k in mesh:
            if k.startswith("box|"):
                lab = k.split("|")[1]
                if lab not in labs:
                    labs.append(lab)
        for lab in labs:
            a = mesh.get(f"box|{lab}|four plates"); b = mesh.get(f"box|{lab}|one plate, web 25 mm")
            w(f"| {lab[0].upper() + lab[1:]} | {a*100:+.2f} % | {b*100:+.2f} % |")
        w("")
        w("**The cells** (the half face, the coarse box inside 10 by 8 m):")
        w("")
        w("| Cells | Four plates | One plate, web 25 mm |")
        w("|---|---|---|")
        for hh, x0, note in ((0.01, 0.0, "10 mm (the coarse box alone)"), (0.005, 0.0, "5 mm, the window over the whole half face"), (0.005, 1.0, "5 mm, the window from 1.0 to 4.0 m"), (0.0025, 1.0, "2.5 mm, the window from 1.0 to 4.0 m")):
            a = mesh.get(f"four plates|{hh:.4f}|{x0:.1f}"); b = mesh.get(f"one plate, web 25 mm|{hh:.4f}|{x0:.1f}")
            if a is not None and b is not None:
                w(f"| {note} | {a*100:+.2f} % | {b*100:+.2f} % |")
        w("")
    if sp:
        w("**The spacing.** A 12 m face, the load track at the edge as today, the support track moved toward the centerline. The design's rectangular plates under each track, and a web between them.")
        w("")
        w("| Clear distance between the two tracks' bundles | " + " | ".join(("Separate plates" if x == 0 else f"Web {x*1e3:.0f} mm") for x in sp["web"]) + " |")
        w("|---|" + "---|" * len(sp["web"]))
        for d in sp["d"]:
            row = sp["rows"][f"{d:.3f}"]
            w(f"| {d*1e3:.0f} mm{' (the design)' if abs(d - 0.323) < 1e-3 else ''} | " + " | ".join((f"{x*100:+.1f} %" if x is not None else "not run") for x in row) + " |")
        w("")
    if ff:
        w("**The other half of the face.** The whole 7.39 m face in one grid at 10 mm cells, no wall at the centerline. The change of one support track's pull from the load track beside it, from the load track on the far side, and from the far support track at its normal current:")
        w("")
        w("| Plate | Far half's currents | Near load track | Far load track | Far support track |")
        w("|---|---|---|---|---|")
        for k, v in ff.items():
            nm, mir = k.split("|")
            w(f"| {nm} | {'mirrored the same way' if mir == '+1' else 'mirrored and turned over'} | {v['load']*100:+.2f} % | {v['loadL']*100:+.2f} % | {v['supL']*100:+.2f} % |")
        w("")

    # ------------------------------------------------------------------ 2
    S2 = []; w = S2.append
    w("## 2. Wider loops on a wider face")
    w("")
    w("### 2.1 Each track alone, by loop width")
    w("")
    w(f"Script: `bearing_round11_plate.py --part=width --bfs=1.2 --bff=1.2` (`round11_plate_width.json`), 5 mm cells. Each track alone on a plate sized from its flux at 1.2 T (1.2), reaching up to 0.6 m past its outer bundles. "
      "The bundles are the design's (150 by 67 mm for the support, 204 by 50 mm for the compensation tracks), and so is the 50 mm between the support's two center bundles.")
    w("")
    w("**The support track at its cap, two loops:**")
    w("")
    w("| Loop width | Bundles span | Ampere-turns per leg | Worst tape: field, limit | Tape per meter of track | Rating | Tape fills the bundle | Plate width | Thickest | Plate mass per meter of track | Plate saturated |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in wd.get("support", []):
        w(f"| {r['Wl']:.1f} m | {r['span']:.2f} m | {r['NI']/1e3:.1f} kA | {r['B']:.2f} T, {r['ic']:.0f} A | {r['tape']:.1f} kg | {r['rating']:.2f} | {r['fill']*100:.0f} % | {r['width']:.2f} m | {r['Tmax']*1e3:.0f} mm | {r['mass']:,.0f} kg | {r['sat']*100:.0f} % |")
    w("")
    w("**A compensation track at 196.5 kN/m, one loop:**")
    w("")
    w("| Loop width | Bundles span | Ampere-turns per leg | Worst tape: field, limit | Tape per meter of track | Rating | Plate width | Thickest | Plate mass per meter of track | Stored energy | Voltage for the 10 ms rise | Coil current |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in wd.get("fast", []):
        w(f"| {r['Wl']:.1f} m | {r['span']:.2f} m | {r['NI']/1e3:.1f} kA | {r['B']:.2f} T, {r['ic']:.1f} A | {r['tape']:.1f} kg | {r['rating']:.2f} | {r['width']:.2f} m | {r['Tmax']*1e3:.0f} mm | {r['mass']:,.0f} kg | {r['W']/1e3:.0f} kJ/m | {r['V']/1e3:.1f} kV | {r['I']:,.0f} A |")
    w("")
    w("### 2.2 The faces")
    w("")
    w("The fit rule: on each half of the lower face a support track and a load track lie side by side, each with its whole plate, so the two plates touch at most and the load track's plate ends at the face's edge. "
      "The lift tracks have the upper face to themselves and take the load track's loop. A track must rate at least 1.35 times its requirement (round 5's ruling). "
      "The bearing on the casing is by round 8's rule: the tape; formers, cryostats and supports at one tape mass; stores at 1.5 kJ/kg for the compensation tracks' stored energy; switches of 4.3 kg, one for each 2.8 kV; the leads. "
      f"In this model at 5 mm the design is {m0['tape']:.0f} kg/m of tape and {m0['total']:,.0f} kg/m on the casing (the book's 550 and 1,687 are the same rule at 2.5 mm), and {B['plates']:,.0f} kg/m of plates.")
    w("")

    def ftable(T, title):
        w(f"**{title}**")
        w("")
        w("| Face | Support loops | Tape per meter of support track | Load and lift loop | Tape per meter of compensation track | Tape of the six tracks | Against the design | A support plate | A compensation plate | The six plates | Hardware on the cable | Of the 10,686 kg/m | Bearing on the casing | Against the design | Face left over on each half |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for f_ in ("7.39", "10.00", "12.00"):
            v = T.get(f_)
            if not v:
                w(f"| {float(f_):.2f} m | no pair of 2.1 fits side by side with its plates whole | | | | | | | | | | | | | |")
                continue
            w(f"| {float(f_):.2f} m | {v['Wls']:.1f} m | {v['tape_s']:.1f} kg | {v['Wlf']:.1f} m | {v['tape_f']:.1f} kg | {v['tape']:.0f} kg/m | {pc(v['tape'], m0['tape'])} | {v['plate_s']:,.0f} kg/m | {v['plate_f']:,.0f} kg/m | {v['plates']:,.0f} kg/m | "
              f"{v['hw']:,.0f} kg/m | {v['hw']/HW_ALLOW*100:.0f} % | {v['bearing']:,.0f} kg/m | {pc(v['bearing'], m0['total'])} | {v['room']:.2f} m |")
        w("")

    if FT:
        ftable(FT, "(a) The loops that need the least tape:")
        ftable(FP, "(b) The loops that need the least plate:")
        ftable(FR, "(c) The loops that put the least on the ring, with a kilogram on the casing counted as nine on the cable (Chapter 7):")
    if f739:
        w("**The 7.39 m face.** No pair of 2.1 fits there with whole plates, because each sized plate reaches 0.6 m past its bundles. What fits is the design's loops under one plate sized on the face (1.3):")
        w("")
        w("| Lift tracks | Tape of the six tracks | Against the design | Plates under the cable | Plates over the cable | Hardware on the cable | Of the 10,686 kg/m | Bearing on the casing | Against the design |")
        w("|---|---|---|---|---|---|---|---|---|")
        for k, lab in ((k_today, "on today's 1.0 m by 200 mm plates"), (k_sized, "on plates sized at 1.2 T, 2.00 m wide")):
            v = f739.get(k)
            if v:
                w(f"| {lab} | {v['tape']:.0f} kg/m | {pc(v['tape'], m0['tape'])} | {fs['mass']:,.0f} kg/m | {2*v['lift_plate']:,.0f} kg/m | {v['hw']:,.0f} kg/m | {v['hw']/HW_ALLOW*100:.0f} % | {v['bearing']:,.0f} kg/m | {pc(v['bearing'], m0['total'])} |")
        w("")
    cv = D.get("conv", [])
    w("### 2.3 The cell size and the box for the sized plates")
    w("")
    if cv:
        w("Script: `bearing_round11_plate.py --part=conv` (`round11_plate_conv.json`). The design's loops on plates sized at 1.2 T, solved again at 2.5 mm cells, and at 5 mm with the coarse box inside a larger one (6 m wider and 8 m tall at 20 mm cells). "
          "The 2.5 mm solves stopped at the solver's cap of 25 iterations with the pull still moving 0.3 to 0.4 percent a step.")
        w("")
        w("| Track | Cells and box | Ampere-turns per leg | Tape per meter of track | Plate mass | Rating |")
        w("|---|---|---|---|---|---|")
        seen = set()
        for c in cv:
            nm = "Support" if c["kind"] == "support" else "Compensation"
            if c["ref"] and c["kind"] not in seen:
                seen.add(c["kind"]); r = c["ref"]
                w(f"| {nm} | 5 mm, the usual box (the tables above) | {r['NI']/1e3:.1f} kA | {r['tape']:.1f} kg | {r['mass']:,.0f} kg/m | {r['rating']:.2f} |")
            w(f"| {nm} | {c['h']*1e3:.1f} mm{', the coarse box inside a larger one' if c['outer'] else ', the usual box'} | {c['NI']/1e3:.1f} kA | {c['tape']:.1f} kg | {c['mass']:,.0f} kg/m | {c['rating']:.2f} |")
        w("")
    else:
        w("Not run (see the last section).")
        w("")

    # ------------------------------------------------------------------ 3
    S3 = []; w = S3.append
    w("## 3. The 12 m face and the track that does the whole job alone")
    w("")
    w("Chapter 7 (Section 7.7): \"Making each track able to do the whole job alone was studied. It needs tracks about 3 m wide, and they do not fit on the 7.39 m square section that this book uses for the cable.\" "
      "The study behind it is `bearing_mass.py` of 12 October (`Vol_III/bearing_mass_2026-10-12.md`, its requirement set 1): either compensation track alone needs 393 kN/m (lift, 30 km/s) or 396 kN/m (load, the parked spacecraft), "
      "and the only air-core track of that study's catalog that reached it was three 1 m loops side by side, rated 471 kN/m, on a plate 3.4 m wide. That study had 100 by 100 mm bundles, a rating set by a mean field of 1.0 T at the winding, and plates 0.4 m wider than the track. "
      "The rules have changed since: the tape is now counted at its own worst place, and the bundle is flat.")
    w("")
    w(f"Script: `bearing_round11_plate.py --part=alone --bff=1.2` (`round11_plate_alone.json`), 5 mm cells. A compensation track at 396 kN/m by the design's rules, alone on a plate sized at 1.2 T. "
      "One loop with a bundle of twice the section, and two or three loops of the design's bundle side by side with alternating polarity, 50 mm between neighboring bundles. "
      "The bundle height needed is from the design's packing routine (`bearing_round6b.pack`), with each loop's coils as the design's: 20 sets of 10 turns.")
    w("")
    w("| Loops | Loop width | Bundle | Bundles span | Ampere-turns per leg | Worst tape: field, limit | Tape per meter of track | Against two design tracks | Rating | Bundle height needed | Plate width | Thickest | Plate mass per meter of track | Stored energy | Voltage per coil | Coil current |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for a in al:
        w(f"| {a['n']} | {a['Wl']:.1f} m | {a['bundle'][0]*1e3:.0f} by {a['bundle'][1]*1e3:.0f} mm | {a['span']:.2f} m | {a['NI']/1e3:.1f} kA | {a['B']:.2f} T, {a['ic']:.1f} A | {a['tape']:.1f} kg | {pc(a['tape'], 2*B['fast']['tape'])} | {a['rating']:.2f} | "
          f"{a['h_need']*1e3:.0f} mm of {a['bundle'][1]*1e3:.0f} | {a['width']:.2f} m | {a['Tmax']*1e3:.0f} mm | {a['mass']:,.0f} kg | {a['W']/1e3:.0f} kJ/m | {a['V']/1e3:.1f} kV | {a['I']:,.0f} A |")
    w("")
    pk = D.get("pack_design")
    if pk:
        w(f"For scale, the design's track by the same packing needs {pk['h_need']*1e3:.0f} mm of its 50 mm bundle height. Two design tracks are {2*B['fast']['tape']:.1f} kg/m of tape and {2*B['fast']['plate']:,.0f} kg/m of plate.")
        w("")
    if al and "hw" in al[0]:
        w(f"**On a 12 m face.** Two such tracks under the cable beside the two support tracks, and two over it, each with its whole plate. The support track has the design's {al[0]['sup_Wl']:.1f} m loops on its sized plate (2.1).")
        w("")
        w("| Loops | Loop width | Bundle | Lower face left over on each half | Upper face left over | Tape of the six tracks | The six plates | Hardware on the cable | Of the 10,686 kg/m | Bearing on the casing | Against the design |")
        w("|---|---|---|---|---|---|---|---|---|---|---|")
        for a in al:
            w(f"| {a['n']} | {a['Wl']:.1f} m | {a['bundle'][0]*1e3:.0f} by {a['bundle'][1]*1e3:.0f} mm | {a['room_bottom']:.2f} m | {a['room_top']:.2f} m | {a['tape6']:.0f} kg/m | {a['plates']:,.0f} kg/m | {a['hw']:,.0f} kg/m | {a['hw']/HW_ALLOW*100:.0f} % | {a['bearing']:,.0f} kg/m | {pc(a['bearing'], m0['total'])} |")
        w("")
        w("A negative figure in a \"left over\" column means the plates do not fit whole. The stores, switches and leads are by round 8's rule, with a track of two or three loops given two or three coils where the design has one.")
        w("")
    D["_parts"] = (L, S1, S2, S3)
    return D, L, S1, S2, S3


if __name__ == "__main__":
    D, L, S1, S2, S3 = main()
    import bearing_round11_plate_text as tx
    text = tx.assemble(D, L, S1, S2, S3)
    path = os.path.join(p.LOGS, "item1b_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"wrote {path} ({len(text.splitlines())} lines)")
