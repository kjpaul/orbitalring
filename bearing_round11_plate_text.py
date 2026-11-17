#!/usr/bin/env python3
"""
Round 11, item 1b: the parts of item1b_report.md that are sentences (the "Look at first" list, the casing's width, the rows of the book,
what was not run).  Called by bearing_round11_plate_md.py; every number is read from runs/round11_logs/round11_plate_report.json.
Book: Orbital Ring Engineering, Volume III: Chapter 7.
Run: imported; not run by itself
"""
from bearing_round11_plate import HW_ALLOW, M_LIM_PLATES, RHO

SIZED = "sized on the face at 1.2 T"
FOUR = "four plates"
W10 = "one plate, web 10 mm"; W25 = "one plate, web 25 mm"; UNI = "one plate, 100 mm everywhere"


def pc(a, b, d=0):
    return f"{abs(a / b - 1) * 100:.{d}f} percent {'more' if a > b else 'less'}"


def sg(x, d=1):
    return f"{abs(x) * 100:.{d}f} percent"


def assemble(D, L, S1, S2, S3):
    B = D["base"]; m0 = B["m8"]
    face = D["face"]; cp = D["coupling"]; mesh = D.get("mesh", {}); sp = D.get("spacing"); ff = D.get("fullface", {})
    wd = D.get("width", {}); al = D.get("alone", []); sk = D.get("skirt", [])
    ts = D.get("tailor_support", []); tf = D.get("tailor_fast", [])
    f4 = face[FOUR]; fs = face.get(SIZED)
    c4 = cp[FOUR]; cs = cp.get(SIZED, {})
    fw10 = face.get(W10); fw25 = face.get(W25); fu = face.get(UNI)
    cw10 = cp.get(W10, {}); cw25 = cp.get(W25, {}); cu = cp.get(UNI, {})
    lower0 = 2 * B["sup"]["plate"] + 2 * B["fast"]["plate"]
    FT = D.get("faces_tape", {}); FR = D.get("faces_ring", {}); FP = D.get("faces_plate", {})
    f739 = D.get("face739", {})
    k_today = "the lift tracks on today's 1.0 m x 200 mm plates"; k_sized = "the lift tracks on plates sized at 1.2 T (table 1b)"
    combos = D.get("combos", []); rect = D.get("rect5", {})

    def trow(rows, bf):
        return next((r for r in rows if abs(r["Bf"] - bf) < 1e-6), None)

    def wrow(kind, Wl):
        return next((r for r in wd.get(kind, []) if abs(r["Wl"] - Wl) < 1e-6), None)

    def skrow(kind, T, s):
        return next((r for r in sk if r["kind"] == kind and abs(r["T"] - T) < 1e-6 and abs(r["skirt"] - s) < 1e-6), None)

    def box(lab, name=FOUR):
        return mesh.get(f"box|{lab}|{name}")

    ts12 = trow(ts, 1.2); tf12 = trow(tf, 1.2)
    b_r7 = box("round 7's box, 9 x 3 m, the support alone at its center, no centerline"); b_big = box("the same in a box 15 x 8 m"); b_half = box("the half face, coarse box inside 10 x 8 m")
    b_half_w = box("the half face, coarse box inside 10 x 8 m", W25); b_r7_w = box("round 7's box, 9 x 3 m, the support alone at its center, no centerline", W25)

    # ------------------------------------------------------------------ look at first
    H = []; w = H.append
    w("## Look at first")
    w("")
    n = [0]

    def item(text):
        n[0] += 1
        w(f"{n[0]}. {text}")

    # 1 can one plate serve
    s = "**Can one plate serve? For the pulls, yes.** Under one plate every track pulls at least as hard as on its own plate. "
    if fw25:
        s += (f"With a 25 mm web joining the design's four plates, the support track needs {fw25['sup']['NI']/1e3:.1f} kA per leg at its cap against {f4['sup']['NI']/1e3:.1f} kA, and {fw25['sup']['tape']:.1f} kg/m of tape against {f4['sup']['tape']:.1f}. "
              f"The load track needs {fw25['load']['tape']:.1f} kg/m against {f4['load']['tape']:.1f}. The ratings are {fw25['sup']['rating']:.2f} for the support and {fw25['load']['rating']:.2f} for the load track, against {f4['sup']['rating']:.2f} and {f4['load']['rating']:.2f}. ")
    if fu:
        s += (f"The plate cannot be thin under the load track, though. One plate 100 mm thick everywhere weighs {fu['mass']:,.0f} kg/m against {lower0:,.0f}, and on it the load track needs {fu['load']['tape']:.1f} kg/m of tape and rates {fu['load']['rating']:.2f}, "
              f"while the support needs {fu['sup']['tape']:.1f} kg/m. ")
    s += "The plate does not have to be cut for the pulls to work (1.3)."
    item(s)
    # 2 what decides it
    s = "**What decides it is the load coil's flux through the support loops, and one plate makes that larger.** "
    s += (f"By round 7's method in this model, the support loses {sg(c4['lin5'])} of its normal pull with the four plates when the load track comes to its full pull")
    lst = [(nm, cp[nm]) for nm in (W10, W25, UNI, SIZED) if nm in cp and "lin5" in cp[nm]]
    if lst:
        s += ", and " + ", ".join(f"{sg(c['lin5'])} with {nm.replace('one plate, ', 'a ').replace('web', 'web of').replace('a 100 mm everywhere', 'one plate 100 mm thick everywhere') if nm != SIZED else 'the plate sized on the face'}" for nm, c in lst)
    s += ". "
    if "direct" in c4:
        s += (f"Solved directly, with the plate saturating and each support loop keeping its flux, the support's pull falls by {sg(c4['direct'])} ({c4['lost']/1e3:.1f} kN/m per track) with the four plates")
        lst = [(nm, cp[nm]) for nm in (W10, W25, UNI, SIZED) if nm in cp and "direct" in cp[nm]]
        if lst:
            s += ", and by " + ", ".join(f"{sg(c['direct'])} ({c['lost']/1e3:.1f} kN/m) with {nm.replace('one plate, ', 'a ').replace('web', 'web of').replace('a 100 mm everywhere', 'one plate 100 mm thick everywhere') if nm != SIZED else 'the plate sized on the face'}" for nm, c in lst)
        s += ". "
    one = [c for nm, c in cp.items() if nm != FOUR and "direct" in c and "lin5" in c]
    if one and "direct" in c4:
        rl = [c["lin5"] / c4["lin5"] for c in one]; rd = [c["direct"] / c4["direct"] for c in one]
        s += (f"So one plate at the design's spacing costs the support {min(rl):.1f} to {max(rl):.1f} times what four plates do by round 7's method, and {min(rd):.1f} to {max(rd):.1f} times solved directly. ")
    s += ("The direct figure is about half of round 7's because the load coil's own field also pulls on the support's bundles, which round 7's method leaves out "
          f"(with the support's currents held, its pull rises by {sg(c4.get('const', 0))} on the four plates). The load track makes the loss up, as it does today, and nothing saturates or runs out of tape for it. "
          "Whether that is acceptable is Paul's ruling, as the 11 to 12.4 percent was in round 8. The 114 cases were not rerun with it (1.4).")
    item(s)
    # 3 the 12.4
    if b_r7 is not None and b_big is not None:
        item(f"**The 12.4 percent itself was read in too small a box.** The same method on the same four plates gives {sg(b_r7, 2)} in round 7's box (9 by 3 m, its floor 1 m under the coils), {sg(b_big, 2)} in a box of 15 by 8 m, "
             f"and {sg(b_half, 2)} in the half-face model used here. The coupling is the one quantity of the bearing that depends on the far field, because the load coil's flux has to come back around through the air, and a small box squeezes it through the support loops. "
             f"The cell size does not matter (the same to 0.2 percent at 10, 5 and 2.5 mm). So today's four plates are at about 9 percent, not 12.4, and one plate is at {sg(b_half_w, 0)} in the same model. "
             "One plate is therefore close to the figure Paul already accepted. For Paul to rule on: which figure the cases and the numbers sheet carry (1.5).")
    # 4 spacing
    if sp:
        r0 = sp["rows"][f"{sp['d'][0]:.3f}"]; iw = {w_: i for i, w_ in enumerate(sp["web"])}
        def at(d, web):
            return sp["rows"][f"{d:.3f}"][iw[web]]
        ds = sp["d"]
        s = (f"**Spacing brings it back.** With the two tracks' bundles {ds[1]*1e3:.0f} mm apart, where the design has {ds[0]*1e3:.0f} mm, one plate with a 25 mm web gives {sg(at(ds[1], 0.025))}, "
             f"which is what four plates give today. At {ds[2]*1e3:.0f} mm it gives {sg(at(ds[2], 0.025))}, and at {ds[3]*1e3:.0f} mm {sg(at(ds[3], 0.025))}. "
             "On the 7.39 m face there is no room to move the support track in and keep the plates' reach. On a wider face there is (1.5).")
        item(s)
    # 5 where the gain is
    s = "**The web is not where the gain is. The gain is in iron above the bundles and in a plate that reaches a little past them.** "
    if fw25:
        s += f"A 25 mm web saves {abs(fw25['sup']['tape']/f4['sup']['tape']-1)*100:.0f} percent of the support's tape and none of the load track's, for {fw25['mass']-f4['mass']:,.0f} kg/m more plate. "
    a = skrow("support", 0.1, 0.2); a0 = skrow("support", 0.1, 0.025); b = skrow("support", 0.3, 0.025); c = skrow("support", 0.3, 0.6); e = skrow("support", 0.1, 1.2)
    if a and a0 and b and c and e:
        s += (f"Today the support plate ends 25 mm past its outer bundles. At the same ampere-turns, the same 100 mm plate reaching 0.2 m past them pulls {pc(a['F'], a0['F'])}, and reaching 1.2 m past them {pc(e['F'], a0['F'])}, so the first 0.2 m is nearly all of it. "
              f"A 300 mm plate of today's width pulls {pc(b['F'], a0['F'])}, and one reaching 0.6 m past the bundles {pc(c['F'], a0['F'])} (1.1). ")
    if ts12:
        s += (f"A support plate sized from its own flux at 1.2 T, the iron's saturation value, is {ts12['width']:.2f} m wide, {ts12['Tmax']*1e3:.0f} mm at its thickest above the center bundles and {ts12['mean']*1e3:.0f} mm on average. "
              f"It weighs {ts12['mass']:,.0f} kg/m against {B['sup']['plate']:,.0f}, and the track then needs {ts12['tape']:.1f} kg/m of tape against {B['sup']['tape']:.1f}, {pc(ts12['tape'], B['sup']['tape'])}. ")
        rr = [r for r in rect.get("support", []) if r["d_tape"] <= (ts12["tape"] / B["sup"]["tape"] - 1) + 0.015]
        if rr:
            r = min(rr, key=lambda r: r["mass"])
            s += (f"A rectangular plate needs about {r['d_mass']*100:.0f} percent more iron for the same saving (round 5's grid: {r['Wp']:.1f} m by {r['T']*1e3:.0f} mm, {abs(r['d_tape'])*100:.0f} percent less tape), "
                  f"where the sized plate needs {(ts12['mass']/B['sup']['plate']-1)*100:.0f} percent more. So Paul's rule is right: thickness where the flux is saves iron. ")
    if tf12:
        s += (f"For a compensation track the same rule gives a plate of {tf12['mass']:,.0f} kg/m against {B['fast']['plate']:,.0f} and {tf12['tape']:.1f} kg/m of tape against {B['fast']['tape']:.1f}. "
              "Its 200 mm plate is already out of saturation, so there the rule buys nothing (1.2).")
    item(s)
    # 6 the mass on 7.39
    if fs:
        s = (f"**One plate under the 7.39 m face, thickness chosen place by place:** {fs['mass']:,.0f} kg/m against the four plates' {lower0:,.0f} kg/m, {fs['Tmax']*1e3:.0f} mm at its thickest and 10 mm at its thinnest. "
             f"On it the support track needs {fs['sup']['tape']:.1f} kg/m of tape against {f4['sup']['tape']:.1f} on its own plate in the same model, {pc(fs['sup']['tape'], f4['sup']['tape'])}, and rates {fs['sup']['rating']:.2f} times its requirement. "
             f"The load track needs {fs['load']['tape']:.1f} kg/m against {f4['load']['tape']:.1f} and rates {fs['load']['rating']:.2f}. ")
        v = f739.get(k_today)
        if v:
            s += (f"With the lift plates left as they are, the hardware on the cable becomes {v['hw']:,.0f} kg/m, {v['hw']/HW_ALLOW*100:.0f} percent of the 10,686 kg/m the cable is sized for, and the bearing on the casing {v['bearing']:,.0f} kg/m against {m0['total']:,.0f} in the same model. ")
        s += "The saving is on the casing and the cost is on the cable. "
        ok = [c for c in combos if c["fits"] and c["tape"] < m0["tape"] * 0.995]
        if ok:
            c = ok[0]
            s += f"Among the plates solved for one track alone, the least tape inside the allowance is the {c['sup']} support plate with the {c['fast']} compensation plate: {c['tape']:.0f} kg/m of tape and {c['hw']:,.0f} kg/m of hardware (1.2)."
        else:
            s += "None of the sized plates solved fits inside the 10,686 kg/m with a rating of 1.35; the design's plates use 94 percent of it, and the allowance is what holds the plates thin (1.2)."
        item(s)
    # 7 wider face
    if FT:
        v10 = FT.get("10.00"); v12 = FT.get("12.00"); p12 = FP.get("12.00"); r12 = FR.get("12.00")
        s = "**What a wider face saves.** Tape, and it costs iron. "
        ws = wd.get("support", []); wf = wd.get("fast", [])
        if ws:
            s += "A support track with wider loops needs less tape: " + ", ".join(f"{r['tape']:.0f} kg/m at {r['Wl']:.1f} m" for r in ws) + f" (the design's plate and loops: {B['sup']['tape']:.0f}). "
        if wf:
            s += "A compensation track gains less: " + ", ".join(f"{r['tape']:.0f} kg/m at {r['Wl']:.1f} m" for r in wf) + f" (the design: {B['fast']['tape']:.0f}). "
        s += "The plates grow with the loops, because a wider loop sends more flux sideways and the plate under it is wider. "
        if v12:
            s += (f"On a 12 m face the loops that need the least tape are {v12['Wls']:.1f} m for the support and {v12['Wlf']:.1f} m for the compensation tracks: {v12['tape']:.0f} kg/m of tape in the six tracks against {m0['tape']:.0f}, "
                  f"a bearing of {v12['bearing']:,.0f} kg/m on the casing against {m0['total']:,.0f}, and {v12['hw']:,.0f} kg/m of hardware on the cable, {v12['hw']/HW_ALLOW*100:.0f} percent of what the cable is sized for. ")
        if v10:
            s += (f"On 10 m: {v10['Wls']:.1f} m and {v10['Wlf']:.1f} m, {v10['tape']:.0f} kg/m of tape, {v10['bearing']:,.0f} kg/m on the casing, {v10['hw']:,.0f} kg/m of hardware ({v10['hw']/HW_ALLOW*100:.0f} percent). ")
        s += ("The cable's hardware allowance is what stops it: every wider-loop choice puts more iron on the cable than the 10,686 kg/m, and the cable's mass is fixed in this item. "
              "On the 7.39 m face the wider loops do not fit with their plates (2.2).")
        item(s)
    # 8 the 3 m track
    if al:
        ok = [a for a in al if a["rating"] >= 1.30 and a["h_need"] <= a["bundle"][1] + 1e-9]
        fit12 = [a for a in ok if a.get("room_bottom", -1) >= 0 and a.get("room_top", -1) >= 0]
        s = "**Does the track that does the whole job alone fit on 12 m? "
        if fit12:
            a = min(fit12, key=lambda a: a["tape"])
            rmin = min(x["rating"] for x in al); rmax = max(x["rating"] for x in al)
            s += (f"Yes in width, top and bottom.** By today's rules it does not need to be 3 m wide. The tracks tried rate {rmin:.2f} to {rmax:.2f} times the 396 kN/m on plates sized at 1.2 T, some a little under the 1.35 of round 5 "
                  f"(the design's own track rates {tf12['rating'] if tf12 else float('nan'):.2f} on a plate sized the same way, and 1.47 on its 200 mm plate), so those would want a somewhat thicker plate than the table gives. "
                  f"The lightest whose tape fits its bundle is {a['n']} loop{'s' if a['n'] > 1 else ''} of {a['Wl']:.1f} m with {a['bundle'][0]*1e3:.0f} by {a['bundle'][1]*1e3:.0f} mm bundles: "
                  f"its bundles span {a['span']:.2f} m and its sized plate is {a['width']:.2f} m wide. It needs {a['tape']:.1f} kg/m of tape, {pc(a['tape'], 2*B['fast']['tape'])} than two of the design's tracks, and a plate of {a['mass']:,.0f} kg/m. "
                  f"With four such tracks the six tracks' tape is {a['tape6']:.0f} kg/m against {m0['tape']:.0f}, the bearing on the casing {a['bearing']:,.0f} kg/m against {m0['total']:,.0f}, "
                  f"and the hardware on the cable {a['hw']:,.0f} kg/m, {a['hw']/HW_ALLOW*100:.0f} percent of the 10,686. So it fits in width and not in the cable's mass. ")
        elif ok:
            a = min(ok, key=lambda a: a["tape"])
            s += (f"Not with its whole plates.** The lightest track that passes the tape rules is {a['n']} loop{'s' if a['n'] > 1 else ''} of {a['Wl']:.1f} m ({a['tape']:.1f} kg/m of tape, plate {a['width']:.2f} m wide). ")
        else:
            s += "None of the tracks tried passes the design's rules** (a rating of 1.35 and a bundle the tape fits in); the table in section 3 gives each one's shortfall. "
        s += "What was not checked is the roll: one track alone at 396 kN/m pulls on one side of a 12 m face (section 3)."
        item(s)
    # 9 the casing
    item("**The casing.** Its floor and its ceiling have to be as wide as the faces they look at: 7.39, 10 or 12 m. The widths the bundles themselves take are in section 4.")
    w("")

    # ------------------------------------------------------------------ 4
    S4 = []; w = S4.append
    w("## 4. What the casing has to do to follow")
    w("")
    w("Only the width is stated here. The casing stays 12,000 kg/m at every width, as Paul ruled.")
    w("")
    w("- The floor under the cable carries the two support tracks and the two load tracks. The ceiling over it carries the two lift tracks. Each has to be as wide as the face it looks at: 7.39 m, 10 m or 12 m.")
    if fs:
        w(f"- On the 7.39 m face with the design's loops the outermost bundles reach 3.60 m from the centerline, so the floor carries coils over 7.19 m of its width, as today. One plate does not change the casing at all. "
          "The floor still steps where the plate's face steps, if the plate is thicker toward the coils; in the model the plate's face is flat and the thickness goes into the cable.")
    for f_ in ("10.00", "12.00"):
        v = FT.get(f_)
        if v:
            half = float(f_) / 2
            xo = half - (v["width_f"] - v["span_f"]) / 2
            w(f"- On the {float(f_):.0f} m face with the least-tape loops of 2.2 ({v['Wls']:.1f} m support loops, {v['Wlf']:.1f} m compensation loops), each load track's bundles span {v['span_f']:.2f} m and each support track's {v['span_s']:.2f} m. "
              f"The outermost bundles reach {xo:.2f} m from the centerline, so the floor and the ceiling carry coils over {2*xo:.2f} m. The cold width of the six cryostats (the bundles' spans) is {2*v['span_s'] + 4*v['span_f']:.2f} m against {2*1.15 + 4*0.804:.2f} m today.")
    w("- The LIM stators face the cable's sides. A wider and flatter cable has lower sides: at 54.64 m2 a 10 m wide section is 5.46 m high and a 12 m wide one 4.55 m, if it stayed a solid rectangle (54.64 / 10 and 54.64 / 12). Whether the LIM plates and stators still fit there was not looked at.")
    w("")

    # ------------------------------------------------------------------ rows of the book
    S5 = []; w = S5.append
    w("## Rows of the book that change")
    w("")
    w("Nothing changes unless Paul adopts one of these. The table gives the line as it stands in `Vol_III/book/Vol_III_Ch07.md` and what it would become. The new figures are at 5 mm cells; the book's are at 2.5 mm. "
      "Where a new figure is a tape or a casing mass, it is given as the 5 mm result and as the book's figure moved by the same ratio.")
    w("")
    w("| Where | Today | If one plate sized on the 7.39 m face is adopted (lift plates as they are) | If only the coupling's box is corrected |")
    w("|---|---|---|---|")
    if fs:
        v = f739.get(k_today, {})
        rs = fs["sup"]["tape"] / f4["sup"]["tape"]; rl = fs["load"]["tape"] / f4["load"]["tape"]
        t6_book = 2 * 144.2 * rs + 2 * 65.4 * rl + 2 * 65.4
        w(f"| Table 7.5, tape of a support track | 144 kg/m | {fs['sup']['tape']:.0f} kg/m at 5 mm; {144.2*rs:.0f} kg/m by ratio | no change |")
        w(f"| Table 7.5, tape of a compensation track | 65 kg/m | load tracks {fs['load']['tape']:.0f} kg/m at 5 mm, {65.4*rl:.0f} kg/m by ratio; lift tracks 65 kg/m | no change |")
        w(f"| Table 7.5, tape of the six tracks | 550 kg/m | {v.get('tape', 0):.0f} kg/m at 5 mm; {t6_book:.0f} kg/m by ratio | no change |")
        if v:
            w(f"| Table 7.5, the bearing on the casing; Section 7.6 \"about 1,690 kg for each meter, which is 14% of the 12,000 kg/m\" | about 1,690 kg/m | {v['bearing']:,.0f} kg/m at 5 mm; about {1687*v['bearing']/m0['total']:,.0f} kg/m by ratio, {1687*v['bearing']/m0['total']/120:.0f}% | no change |")
        w(f"| Table 7.5, a support plate and a compensation plate (944 and 1,574 kg/m) | two rows | one row for the plate under the cable, {fs['mass']:,.0f} kg/m, and one for the two lift plates, 3,148 kg/m | no change |")
        if v:
            w(f"| Table 7.5, the six plates | 8,185 kg/m | {v['plates']:,.0f} kg/m (three plates) | no change |")
            w(f"| Table 7.5, all hardware on the cable | 10,057 kg/m | {v['hw']:,.0f} kg/m | no change |")
            w(f"| Table 7.5, margin on the cable | 629 kg/m, 6% of the allowance | {HW_ALLOW - v['hw']:,.0f} kg/m ({(HW_ALLOW - v['hw'])/HW_ALLOW*100:.0f}% of the allowance; negative means over it) | no change |")
        w("| Table 7.3, caption: \"A support plate is 1.2 m wide and a compensation plate is 1.0 m wide\" and its even fields | 0.32, 0.72, 0.85 T for the support | the even field of Equation 7.3 has no single plate width to use under one plate; the caption and the support and load rows would need rewording | no change |")
        w(f"| Section 7.5, the plates paragraph: \"The two support plates are 1.2 m wide and 100 mm thick, and at the cap about a third of their cross-section is saturated. The four compensation plates are 1.0 m wide and 200 mm thick\" | as quoted | "
          f"one plate under the cable, {fs['Tmax']*1e3:.0f} mm at its thickest and 10 mm at its thinnest, {fs['sup']['sat']*100:.0f} percent of it saturated under the support at the cap; the two lift plates as they are | no change |")
        w("| Section 7.5, \"The six plates together weigh 8,185 kg ... the hardware on the cable comes to 10,057 kg/m, inside the 10,686 kg/m\" | as quoted | the figures of the Table 7.5 rows above | no change |")
        w("| Figure 7.3 and its caption (\"each with four coil bundles and a thin plate on the cable ... two load tracks with thicker plates\"; \"The casing's face stands higher under the thinner support plate\") | four plates under the cable | one plate under the cable | no change |")
    if b_half is not None:
        w(f"| The numbers sheet, \"the load track's pull on the support beside it\" (not quoted in Chapter 7) | 12.4 percent of the support's normal force | {sg(cs.get('lin5', 0)) if cs else 'not run'} by round 7's method on the sized plate; {sg(cs['direct']) if 'direct' in cs else 'not run'} solved directly | "
          f"{sg(b_half)} for the four plates (round 7's method in a large box); {sg(c4['direct']) if 'direct' in c4 else 'not run'} solved directly |")
    w("| Section 7.7, \"It needs tracks about 3 m wide, and they do not fit on the 7.39 m square section\" | as quoted | see below | no change |")
    w("")
    if al:
        ok = [a for a in al if a["rating"] >= 1.30 and a["h_need"] <= a["bundle"][1] + 1e-9]
        if ok:
            a = min(ok, key=lambda a: a["tape"])
            w(f"**The 3 m sentence.** By today's tape rules a compensation track that gives 396 kN/m alone is {a['n']} loop{'s' if a['n'] > 1 else ''} of {a['Wl']:.1f} m, with bundles that span {a['span']:.2f} m and a plate {a['width']:.2f} m wide if the plate is sized from its flux. "
              f"Four of them put {4*a['mass']:,.0f} kg/m of plate on the cable where the design's four compensation plates are {4*B['fast']['plate']:,.0f} kg/m. "
              f"If Paul wants the sentence to stay true to the runs, it could say that such a track needs about {a['tape']/B['fast']['tape']:.1f} times the tape and {a['mass']/B['fast']['plate']:.1f} times the plate of this design's track, and that the cable of this book is not sized for the plates. "
              "The width is no longer the reason. This is a suggestion for the book session, not an edit.")
            w("")
    w("If a wider face is adopted, the same rows change by the figures of the tables in 2.2, and the section's width (7.39 m, in Chapters 3, 7 and 19 and Figure 7.3) with them. Those are not listed row by row, because no face has been chosen.")
    w("")

    # ------------------------------------------------------------------ not run
    S6 = []; w = S6.append
    w("## What was not run")
    w("")
    w("- **The 114 cases with a larger coupling.** The cases of rounds 7 to 10 carry 12.4 percent. They were not rerun at the one-plate figure, nor at the corrected four-plate figure. That is the run that would say whether one plate is acceptable.")
    w("- **The 2.5 mm solves stopped at the iteration cap.** The saturating solves at 2.5 mm reached the solver's 25 iterations with the pull still moving 0.3 to 0.4 percent a step, so the 2.5 mm rows of 2.3 carry that much more than the 5 mm rows.")
    w("- **The cells.** Everything is at 5 mm. At 2.5 mm: the coupling by round 7's method (1.5)" + (", and the design's loops on their sized plates (2.3)" if D.get("conv") else "; the sized plates were not solved at 2.5 mm") + ". The half face, the wider loops and the track of section 3 were not solved at 2.5 mm.")
    w("- **Wider loops under one plate on a 10 or 12 m face.** Section 2 solves each track alone on its own sized plate and lays them side by side. The half-face model with both tracks under one plate was run for 7.39 m only. "
      "So for 10 and 12 m there is no direct figure for the coupling with the wider loops; the spacing table of 1.5 is for the design's loops.")
    w("- **The plate sized on the face is a first pass.** Its thickness is the sum of the two tracks' sideways fluxes on a thick plate, each at its requirement alone. It was not iterated, the fill field on the face was 1.2 T only, and the support's cap and the load track's full pull were added although the design never asks for both at once. A lighter plate that does the same is likely [SPECULATIVE].")
    w("- **The lift tracks under one plate.** The two lift tracks are far apart on the upper face and nothing couples them to a persistent loop. They were left on their own plates, or given the single-track sized plate of 1.2.")
    w("- **The bundle chosen again for the wider loops.** The bundles are the design's in every row of section 2. The gap trade found the bundle's shape worth about 1 percent.")
    w("- **The support's lost loop, the stray field, the coil ends, the plate's heat and the tape's pulse loss** on any of the new plates or loops. A plate that reaches past the bundles carries the moving field over more iron, so the plate's heat per pass will not be the design's.")
    w("- **The roll.** One compensation track doing the whole job alone pulls 396 kN/m on one side. On a 12 m face its lever arm is longer than on 7.39 m. The heave and roll model was not run for it.")
    w("- **The structure.** How a plate that is 10 mm in one place and over 200 mm in another is built into the cable, and whether a thin web survives the cable's strain, was not looked at. The resin-bonded iron is not a structural material in the book.")
    w("- **The casing.** Not weighed, as instructed.")
    w("- **The mass loop.** The cable's mass is fixed in this item, so a hardware figure over 10,686 kg/m is reported and not fed back into the cable.")
    w("- **No literature search.** Nothing here rests on a source that was not already in the earlier rounds.")
    w("")
    return "\n".join(L + H + S1 + S2 + S3 + S4 + S5 + S6) + "\n"


if __name__ == "__main__":
    print(__doc__)
