"""
Draws Figure 7.4, the cross section of the suspension bearing (cable, the two casing support tracks below it, the two
lift compensation tracks above it, the four plates), as an SVG; widths across the cable are to scale, heights are not,
and the cable is drawn as a wide low rectangle so that the labels have room (Paul, 2026-10-24).
Book: Orbital Ring Engineering, Volume I: Chapter 7, Figure 7.4 (file name fig07_04_bearing_cross_section).  Redrawn 2026-10-24 to the
Part 3b coil design (coil_optimization_2026-10-12.md), text enlarged for print and EML renamed EMS on 2026-10-25.
Run: python tools/generate_fig_bearing_cross_section.py   (writes <output folder>/figures_ch07/fig07_04_bearing_cross_section.svg,
     .png and .pdf: cairosvg if it has a cairo library, else the PNG through a headless Edge or Chrome and the PDF through svglib)
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import ring_altitude as ring
# Widths across the cable are in proportion (1680 px for the cable's side, 8.22 m from ring_altitude); gaps, heights and the
# bundles are not to scale, and the cable's height is drawn far below its width.
SIDE = ring.CABLE_SIDE
V_CABLE = ring.V_CABLE
PXM = 1680/SIDE
# the Part 3b coil design: support tracks two 0.50 m loops with 150 x 67 mm bundles, the centre pair touching; lift tracks one 0.80 m loop with 200 x 50 mm bundles
SUP_OFF = (-0.575, -0.075, 0.075, 0.575)
LIFT_W = 0.80
SUP_B, LIFT_B = (150, 67), (200, 50)   # mm
CX0, CX1, CT, CB = 560, 2240, 760, 1480   # cable box: 1680 px wide, 720 px tall (not to scale)
WIDTH = 2800
def w(m): return round(m*PXM)
out = []
A = out.append
RED, RED2 = "#c0392b", "#922b21"
NAVY, NAVY2 = "#1f3a93", "#152a6b"
OLV, OLV2 = "#7d6608", "#5b4a06"
GRN, GRN2 = "#27ae60", "#1e8449"
ORG, ORG2 = "#f39c12", "#d68910"
# font sizes (px): raised about half from the 2026-10-24 figure so that the figure reads at 5 inches wide in print
F_TITLE, F_SUB, F_HEAD, F_TEXT, F_TAG, F_SMALL, F_CABLE = 72, 52, 50, 48, 46, 44, 72

def text(x, y, s, fs, fill="#333", anchor="middle", bold=False, extra=""):
    A(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{fs}" fill="{fill}"{" font-weight=\"bold\"" if bold else ""}{extra}>{s}</text>')
def box(x, y, ww, hh, fill, stroke):
    A(f'<rect x="{x}" y="{y}" width="{ww}" height="{hh}" rx="3" fill="{fill}" stroke="{stroke}" stroke-width="2.5"/>')
def bracket(x, y0, y1, col, txt, fs=F_TAG):
    A(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="{col}" stroke-width="2.5"/>')
    for y in (y0, y1):
        A(f'<line x1="{x-18}" y1="{y}" x2="{x+18}" y2="{y}" stroke="{col}" stroke-width="2.5"/>')
    if txt:
        text(x+28, (y0+y1)/2+fs*0.35, txt, fs, col, "start")
PXMM = 0.5   # bundle rectangles keep their real proportions at 0.5 px per mm (150 x 67 mm -> 75 x 34 px)
BS = round(SUP_B[1]*PXMM)
def bundles(xc, y_bottom, offsets, col, col2, size):
    bw, bh = size[0]*PXMM, size[1]*PXMM
    for o in offsets:
        A(f'<rect x="{xc + w(o) - bw/2}" y="{y_bottom - bh}" width="{bw}" height="{bh}" rx="4" fill="{col}" stroke="{col2}" stroke-width="2.5"/>')

GAP = 100
W_PLATE = 1.4
C = 1400
x_sup = [C - w(1.32), C + w(1.32)]
x_lift = [C - w(3.17), C + w(3.17)]

# ---------- header and casing ----------
body = []
B = body.append
text(C, 84, "Orbital Ring Bearing Cross Section (Not to Scale)", F_TITLE, "#000", bold=True)
text(C, 150, "Cable moves perpendicular to page. Casing is stationary relative to the ground.", F_SUB, "#555")
CAS_TOP = 200
text(C, 282, "CASING (stationary relative to the ground, 12,000 kg/m)", 56, "#333", bold=True)
text(C, 342, "Four air-core EMS tracks hold the casing to the cable across a 100 mm gap", F_TEXT, "#555")

# ---------- top: lift compensation tracks ----------
text(C, 430, "Lift Compensation Tracks (2, driven superconducting coils, attractive)", F_HEAD, NAVY, bold=True)
coil_b = CT - w(0.2) - GAP                  # the bundles' bottom face, one gap above the plate
for x, lab in zip(x_lift, ("Lift L", "Lift R")):
    bundles(x, coil_b, (-LIFT_W/2, LIFT_W/2), NAVY, NAVY2, LIFT_B)
    text(x, coil_b - LIFT_B[1]*PXMM - 18, lab, F_TAG, NAVY, bold=True)
bracket(x_lift[0] + w(W_PLATE)/2 + 44, coil_b, CT - w(0.2), NAVY, "100 mm gap")
A(f'<rect x="{CX0}" y="{CT}" width="{CX1-CX0}" height="{CB-CT}" rx="12" fill="#d4e6f1" stroke="#2471a3" stroke-width="4"/>')
for x in x_lift:
    box(x - w(W_PLATE)/2, CT - w(0.2), w(W_PLATE), w(0.2), OLV, OLV2)
text(C, CT + 96, "Ferromagnetic Plates (2)", F_TEXT, OLV, bold=True)
text(C, CT + 146, "(sintered iron in non-conductive resin, 200 mm thick)", F_SMALL, OLV)

# ---------- cable text ----------
yc = (CT + CB)/2
text(C, yc - 50, "CNT CABLE", F_CABLE, "#1a5276", bold=True)
text(C, yc + 18, f"{SIDE:.2f} m x {SIDE:.2f} m cross section", F_HEAD, "#1a5276")
text(C, yc + 76, f"Moving at {V_CABLE:,.0f} m/s (retrograde)", F_HEAD, "#1a5276")
text(C, yc + 134, "25 GPa CNT yarn, 12.5 GPa operating stress", F_HEAD, "#1a5276")
text(C, CB - 108, "Ferromagnetic Plates (2)", F_TEXT, OLV, bold=True)
text(C, CB - 58, "(sintered iron in non-conductive resin, 200 mm thick)", F_SMALL, OLV)

# ---------- bottom: casing support tracks ----------
coil_t = CB - 12 + w(0.2) + GAP             # the bundles' top face, one gap below the plate
for x, lab in zip(x_sup, ("Support L", "Support R")):
    box(x - w(W_PLATE)/2, CB - 12, w(W_PLATE), w(0.2), OLV, OLV2)
    bundles(x, coil_t + BS, SUP_OFF, RED, RED2, SUP_B)
    text(x, coil_t + BS + 56, lab, F_TAG, RED, bold=True)
bracket(C, CB - 12 + w(0.2), coil_t, RED, "100 mm gap")
y = coil_t + BS + 56
text(C, y + 100, "Casing Support Tracks (2, steady superconducting coils, attractive)", F_HEAD, RED, bold=True)
text(C, y + 160, "Each also carries a driven load compensation winding", F_HEAD, NAVY, bold=True)
text(C, y + 216, "Rectangles are coil bundles (not to scale). There is no iron on the casing. All tracks act vertically only.", F_SMALL, "#555")
yd = y + 290
A(f'<line x1="{CX0}" y1="{yd}" x2="{CX1}" y2="{yd}" stroke="#1a5276" stroke-width="3"/>')
for x in (CX0, CX1):
    A(f'<line x1="{x}" y1="{yd-16}" x2="{x}" y2="{yd+16}" stroke="#1a5276" stroke-width="3"/>')
text(C, yd + 62, f"{SIDE:.2f} m", F_HEAD, "#1a5276")
CASB = yd + 110

# ---------- LIM ----------
# The LIM reaction plates (1.2 m tall, 200 mm thick), the LIM stators and the EMS plates (1.4 m wide, 200 mm thick) share one
# scale, the figure's width scale PXM (204 px per metre), whatever the cable rectangle's own compression (Paul, 2026-10-25).
# The stator is drawn 1.0 m tall (the book says only that it is narrower than the reaction plate) and 0.3 m thick [UNVERIFIED].
T_PLATE = w(0.2); H_RPLATE = w(1.2); H_STATOR = w(1.0); T_STATOR = w(0.3); GAP_LIM = w(0.2)
PT, PB = round(yc - H_RPLATE/2), round(yc + H_RPLATE/2)        # the reaction plates, centred on the drawn cable
ST, SB = round(yc - H_STATOR/2), round(yc + H_STATOR/2)
xL_plate, xR_plate = CX0 - T_PLATE, CX1                           # reaction plates on the cable's faces
xL_stat, xR_stat = xL_plate - GAP_LIM - T_STATOR, xR_plate + T_PLATE + GAP_LIM
for x in (xL_plate, xR_plate):
    A(f'<rect x="{x}" y="{PT}" width="{T_PLATE}" height="{PB-PT}" rx="3" fill="{ORG}" stroke="{ORG2}" stroke-width="2.5"/>')
for x in (xL_stat, xR_stat):
    A(f'<rect x="{x}" y="{ST}" width="{T_STATOR}" height="{SB-ST}" rx="4" fill="{GRN}" stroke="{GRN2}" stroke-width="2.5"/>')
    text(x + T_STATOR/2, ST - 26, "LIM", 52, GRN2, bold=True)
# the 200 mm gap: a dimension line across it under the parts, the labels in a block at the lower left with colour swatches
yg = PB + 36
for xa, xb in ((xL_stat + T_STATOR, xL_plate), (xR_plate + T_PLATE, xR_stat)):
    A(f'<line x1="{xa}" y1="{yg}" x2="{xb}" y2="{yg}" stroke="{GRN}" stroke-width="3"/>')
    for x in (xa, xb):
        A(f'<line x1="{x}" y1="{yg-12}" x2="{x}" y2="{yg+12}" stroke="{GRN}" stroke-width="3"/>')
LX, LY0 = 230, CB + 76
lab = [("LIM stator", GRN2, (xL_stat + T_STATOR/2, SB)), ("200 mm gap", GRN, ((xL_stat + T_STATOR + xL_plate)/2, yg)),
       ("TiAl reaction plate, 1.2 m", ORG2, (xL_plate + T_PLATE/2, PB)), ("offset above center for lateral shear", GRN2, None)]
for i, (t, col, pt) in enumerate(lab):            # the parts are small at this scale: a colour swatch before each label instead of leader lines
    yy = LY0 + 54*i
    A(f'<rect x="{LX}" y="{yy-28}" width="34" height="28" fill="{col}"/>')
    text(LX + 48, yy, t, 40, col, anchor="start")
A(f'<line x1="100" y1="{yc-230}" x2="100" y2="{yc+230}" stroke="#666" stroke-width="4"/><polygon points="100,{yc+230} 86,{yc+200} 114,{yc+200}" fill="#666"/>')
text(100, yc + 300, "Toward", F_TAG, "#666"); text(100, yc + 352, "Earth", F_TAG, "#666")

# ---------- legend ----------
rows = [
 (RED, "Casing Support Tracks (2)", "Two loops of 12 mm HTS tape per track, 0.5 m wide each, four 150 x 67 mm bundles, in persistent mode. They carry the casing."),
 (NAVY, "Load Compensation Windings (2)", "Driven 3 mm HTS tape on the support track loops. They add lift under a slow or parked sled."),
 (NAVY, "Lift Compensation Tracks (2)", "One driven loop 0.8 m wide per track, two 200 x 50 mm bundles, on the casing ceiling. They pull the casing down as a fast sled passes."),
 (NAVY, "Gap Protection", "The driven windings also guard the gap. If it closes on one face, the windings on the other face pull back."),
 (GRN, "LIM Stators", "Deployment thrust and lateral centering. Narrower than the reaction plates (edge effect)."),
 (ORG, "LIM Reaction Plates", "γ-TiAl alloy, 200 mm thick, 1.2 m wide, one on each side. Offset above center with the LIM stators."),
 (OLV, "Ferromagnetic Plates (4)", "Sintered iron in non-conductive resin, 1.4 m wide and 200 mm thick. One for each track."),
]
F_LEG, LINE, ROW_GAP, MAXC = 44, 54, 26, 100     # legend font, line pitch, gap between rows, characters per line
def wrap(name, desc):
    words = (name + ": " + desc).split(); lines = []; cur = ""
    for wd in words:
        if cur and len(cur) + 1 + len(wd) > MAXC:
            lines.append(cur); cur = wd
        else:
            cur = (cur + " " + wd).strip()
    lines.append(cur); return lines
LY = CASB + 60
wrapped = [(col, name, wrap(name, desc)) for col, name, desc in rows]
leg_h = 110 + sum(LINE*len(ls) + ROW_GAP for _, _, ls in wrapped) + 10
H = LY + leg_h + 40
A(f'<rect x="140" y="{LY}" width="2520" height="{leg_h}" rx="10" fill="#f8f9f9" stroke="#bbb" stroke-width="2"/>')
text(C, LY + 76, "Component Summary", 56, "#000", bold=True)
yy = LY + 150
for col, name, lines in wrapped:
    A(f'<rect x="200" y="{yy-30}" width="46" height="32" fill="{col}"/>')
    first = lines[0][len(name)+1:]           # the name goes in bold, the rest of the first line after it
    A(f'<text x="272" y="{yy}" font-size="{F_LEG}" fill="#333"><tspan font-weight="bold" fill="{col}">{name}:</tspan>{first}</text>')
    for ln in lines[1:]:
        yy += LINE
        A(f'<text x="272" y="{yy}" font-size="{F_LEG}" fill="#333">{ln}</text>')
    yy += LINE + ROW_GAP
A('</svg>')

# the header, the casing box and the white background go first; they need the heights computed above
head = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {H}" font-family="Arial, sans-serif">',
        '<title>Orbital Ring Bearing Cross Section</title>',
        f'<rect width="{WIDTH}" height="{H}" fill="white"/>',
        f'<rect x="200" y="{CAS_TOP}" width="2400" height="{CASB-CAS_TOP}" rx="20" fill="#f5f5f5" stroke="#333" stroke-width="5"/>']
out = head + out
OUT = ring.out_dir("figures_ch07")   # under ORBITAL_RING_OUT (default runs)
SVG = _os.path.join(OUT, "fig07_04_bearing_cross_section.svg")
open(SVG, "w", encoding="utf-8").write("\n".join(out) + "\n")
PNG = _os.path.join(OUT, "fig07_04_bearing_cross_section.png"); PDF = _os.path.join(OUT, "fig07_04_bearing_cross_section.pdf")
try:
    import cairosvg
    cairosvg.svg2png(url=SVG, write_to=PNG, output_width=2*WIDTH)
    cairosvg.svg2pdf(url=SVG, write_to=PDF)
except (ImportError, OSError):
    # without cairo: the PNG through a headless Edge or Chrome (keeps the Arial text), else PyMuPDF (substitutes a serif font);
    # the PDF through svglib and reportlab (vector, Helvetica), else PyMuPDF
    import shutil, subprocess
    browser = next((p for p in (r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe", r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
                                shutil.which("msedge"), shutil.which("chrome"), shutil.which("google-chrome"), shutil.which("chromium")) if p and _os.path.exists(p)), None)
    done_png = False
    if browser:
        r = subprocess.run([browser, "--headless=new", "--disable-gpu", "--hide-scrollbars", f"--window-size={WIDTH},{H}", "--force-device-scale-factor=2",
                            f"--screenshot={_os.path.abspath(PNG)}", "file:///" + _os.path.abspath(SVG).replace("\\", "/")], capture_output=True)
        done_png = _os.path.exists(PNG) and _os.path.getsize(PNG) > 0
    try:
        from svglib.svglib import svg2rlg
        from reportlab.graphics import renderPDF
        renderPDF.drawToFile(svg2rlg(SVG), PDF); done_pdf = True
    except ImportError:
        done_pdf = False
    if not (done_png and done_pdf):
        try:
            import fitz                                   # PyMuPDF renders SVG directly
            doc = fitz.open(SVG); page = doc[0]
            if not done_png: page.get_pixmap(matrix=fitz.Matrix(2, 2)).save(PNG)
            if not done_pdf: open(PDF, "wb").write(doc.convert_to_pdf())
        except ImportError:
            print("no renderer found (cairosvg, a headless browser, svglib or PyMuPDF): the SVG was written, no PNG or PDF")
print(f"canvas {WIDTH} x {H}, casing box to {CASB}, legend {leg_h} tall; wrote {SVG}")
