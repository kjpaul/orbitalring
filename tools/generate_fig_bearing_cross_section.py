"""
Draws Figure 8.1, the cross section of the levitation bearing (cable, the two casing support tracks below it, the two
lift compensation tracks above it, the four plates), as an SVG; widths across the cable are to scale, heights are not.
Book: The Orbital Ring, Volume IV: Chapter 8, Figure 8.1.
Run: python tools/generate_fig_bearing_cross_section.py   (writes <output folder>/figures_ch08/fig10_bearing_cross_section.svg, and .png when cairosvg is installed)
Environment: ORBITAL_RING_OUT (output folder, default runs)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import ring_altitude as ring
# Figure 8.1 bearing cross section, 700 km design. Widths across the cable are in proportion
# (204 px per meter, cable 8.24 m); gaps and heights are not to scale.
PXM = 1680/8.24
CX0, CX1, CT, CB = 560, 2240, 640, 1940   # cable box
def w(m): return round(m*PXM)
out = []
A = out.append
RED, RED2 = "#c0392b", "#922b21"
NAVY, NAVY2 = "#1f3a93", "#152a6b"
PUR, PUR2 = "#8e44ad", "#6c3483"
OLV, OLV2 = "#7d6608", "#5b4a06"
TEAL, TEAL2 = "#0e6655", "#0b5345"
GRN, GRN2 = "#27ae60", "#1e8449"
ORG, ORG2 = "#f39c12", "#d68910"
H = 3200
A(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 2800 {H}" font-family="Arial, sans-serif">')
A('<title>Orbital Ring Bearing Cross Section</title>')
A(f'<rect width="2800" height="{H}" fill="white"/>')
A('<text x="1400" y="75" text-anchor="middle" font-size="58" font-weight="bold">Orbital Ring Bearing Cross Section (Not to Scale)</text>')
A('<text x="1400" y="130" text-anchor="middle" font-size="42" fill="#555">Cable moves perpendicular to page. Casing is stationary relative to the ground.</text>')
CASB = 2460
A(f'<rect x="200" y="185" width="2400" height="{CASB-185}" rx="20" fill="#f5f5f5" stroke="#333" stroke-width="5"/>')
A('<text x="1400" y="252" text-anchor="middle" font-size="46" font-weight="bold" fill="#333">CASING (stationary relative to the ground, 12,000 kg/m)</text>')
A(f'<text x="1400" y="304" text-anchor="middle" font-size="38" fill="#555">Four air-core EML tracks hold the casing to the cable across a 100 mm gap</text>')

def box(x, y, ww, hh, fill, stroke, label=None, fs=32):
    A(f'<rect x="{x}" y="{y}" width="{ww}" height="{hh}" rx="3" fill="{fill}" stroke="{stroke}" stroke-width="2.5"/>')
    if label:
        A(f'<text x="{x+ww/2}" y="{y+hh/2+fs*0.35}" text-anchor="middle" font-size="{fs}" fill="white" font-weight="bold">{label}</text>')
def arrow(x, y0, y1, col):
    d = 1 if y1 > y0 else -1
    A(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="{col}" stroke-width="4"/>')
    A(f'<polygon points="{x},{y1} {x-9},{y1-18*d} {x+9},{y1-18*d}" fill="{col}"/>')
def bracket(x, y0, y1, col, txt):
    A(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="{col}" stroke-width="2"/>')
    for y in (y0, y1):
        A(f'<line x1="{x-16}" y1="{y}" x2="{x+16}" y2="{y}" stroke="{col}" stroke-width="2"/>')
    A(f'<text x="{x+24}" y="{(y0+y1)/2+12}" font-size="34" fill="{col}">{txt}</text>')

GAP = 100
W_PLATE = 1.4
C = 1400
x_sup = [C - w(1.32), C + w(1.32)]
x_lift = [C - w(3.17), C + w(3.17)]
cable = f'<rect x="{CX0}" y="{CT}" width="{CX1-CX0}" height="{CB-CT}" rx="12" fill="#d4e6f1" stroke="#2471a3" stroke-width="4"/>'
BS = 34   # drawn size of a coil bundle
def bundles(xc, y, offsets, col, col2):
    for o in offsets:
        A(f'<rect x="{xc + w(o) - BS/2}" y="{y}" width="{BS}" height="{BS}" rx="4" fill="{col}" stroke="{col2}" stroke-width="2.5"/>')
def tag(xc, y, txt, col):
    A(f'<text x="{xc}" y="{y}" text-anchor="middle" font-size="32" fill="{col}" font-weight="bold">{txt}</text>')

# ---------- top ----------
A(f'<text x="{C}" y="392" text-anchor="middle" font-size="38" fill="{NAVY}" font-weight="bold">Lift Compensation Tracks (2, driven superconducting coils, attractive)</text>')
coil_b = CT - 40 - GAP
for x, lab in zip(x_lift, ("Lift L", "Lift R")):
    bundles(x, coil_b - BS, (-0.5, 0.5), NAVY, NAVY2)
    tag(x, coil_b - BS - 14, lab, NAVY)
bracket(x_lift[0] + w(W_PLATE)/2 + 40, coil_b, CT - 40, NAVY, "100 mm gap")
A(cable)
for x in x_lift:
    box(x - w(W_PLATE)/2, CT - 40, w(W_PLATE), 52, OLV, OLV2)
A(f'<text x="{C}" y="{CT+84}" text-anchor="middle" font-size="38" fill="{OLV}" font-weight="bold">Ferromagnetic Plates (2)</text>')
A(f'<text x="{C}" y="{CT+124}" text-anchor="middle" font-size="34" fill="{OLV}">(sintered iron in non-conductive resin, 200 mm thick)</text>')

# ---------- cable text ----------
A(f'<text x="{C}" y="1160" text-anchor="middle" font-size="62" font-weight="bold" fill="#1a5276">CNT CABLE</text>')
A(f'<text x="{C}" y="1235" text-anchor="middle" font-size="44" fill="#1a5276">8.24 m x 8.24 m cross section</text>')
A(f'<text x="{C}" y="1285" text-anchor="middle" font-size="44" fill="#1a5276">Moving at 8,309 m/s (retrograde)</text>')
A(f'<text x="{C}" y="1335" text-anchor="middle" font-size="44" fill="#1a5276">25 GPa CNT yarn, 12.5 GPa operating stress</text>')
A(f'<text x="{C}" y="{CB-140}" text-anchor="middle" font-size="38" fill="{OLV}" font-weight="bold">Ferromagnetic Plates (2)</text>')
A(f'<text x="{C}" y="{CB-100}" text-anchor="middle" font-size="34" fill="{OLV}">(sintered iron in non-conductive resin, 200 mm thick)</text>')

# ---------- bottom ----------
for x, lab in zip(x_sup, ("Support L", "Support R")):
    box(x - w(W_PLATE)/2, CB - 12, w(W_PLATE), 52, OLV, OLV2)
    bundles(x, CB + 40 + GAP, (-0.5, 0.0, 0.5), RED, RED2)
    tag(x, CB + 40 + GAP + BS + 40, lab, RED)
bracket(C, CB + 40, CB + 40 + GAP, RED, "")
A(f'<text x="{C}" y="{CB+40+GAP+44}" text-anchor="middle" font-size="32" fill="{RED}">100 mm gap</text>')
y = CB + 40 + GAP + BS + 40
A(f'<text x="{C}" y="{y+84}" text-anchor="middle" font-size="38" fill="{RED}" font-weight="bold">Casing Support Tracks (2, steady superconducting coils, attractive)</text>')
A(f'<text x="{C}" y="{y+134}" text-anchor="middle" font-size="38" fill="{NAVY}" font-weight="bold">Each also carries a driven load compensation winding</text>')
A(f'<text x="{C}" y="{y+180}" text-anchor="middle" font-size="34" fill="#555">Squares are coil bundles. There is no iron on the casing. All tracks act vertically only.</text>')
yd = y + 240
A(f'<line x1="{CX0}" y1="{yd}" x2="{CX1}" y2="{yd}" stroke="#1a5276" stroke-width="3"/>')
for x in (CX0, CX1):
    A(f'<line x1="{x}" y1="{yd-15}" x2="{x}" y2="{yd+15}" stroke="#1a5276" stroke-width="3"/>')
A(f'<text x="{C}" y="{yd+52}" text-anchor="middle" font-size="42" fill="#1a5276">8.24 m</text>')

# ---------- LIM ----------
for x in (CX0 - 36, CX1):
    A(f'<rect x="{x}" y="840" width="36" height="640" rx="3" fill="{ORG}" stroke="{ORG2}" stroke-width="2.5"/>')
for x in (248, 2492):
    A(f'<rect x="{x}" y="900" width="60" height="520" rx="4" fill="{GRN}" stroke="{GRN2}" stroke-width="2.5"/>')
    A(f'<text x="{x+30}" y="876" text-anchor="middle" font-size="42" fill="{GRN2}" font-weight="bold">LIM</text>')
A(f'<line x1="308" y1="1240" x2="524" y2="1240" stroke="{GRN}" stroke-width="2" stroke-dasharray="6,4"/>')
A(f'<text x="416" y="1130" text-anchor="middle" font-size="36" fill="{GRN}" font-weight="bold" transform="rotate(-90, 416, 1130)">200 mm gap</text>')
for yy in (1000, 1300):
    A(f'<line x1="248" y1="{yy}" x2="416" y2="{yy}" stroke="{GRN}" stroke-width="5"/><polygon points="416,{yy} 396,{yy-10} 396,{yy+10}" fill="{GRN}"/>')
    A(f'<line x1="2552" y1="{yy}" x2="2384" y2="{yy}" stroke="{GRN}" stroke-width="5"/><polygon points="2384,{yy} 2404,{yy-10} 2404,{yy+10}" fill="{GRN}"/>')
A(f'<text x="464" y="1160" text-anchor="middle" font-size="36" fill="{ORG2}" font-weight="bold" transform="rotate(-90, 464, 1160)">TiAl Reaction Plate (1.2 m)</text>')
A(f'<text x="2336" y="1160" text-anchor="middle" font-size="36" fill="{ORG2}" font-weight="bold" transform="rotate(90, 2336, 1160)">TiAl Reaction Plate (1.2 m)</text>')
for i, t in enumerate(("Offset above", "center for", "lateral shear")):
    A(f'<text x="406" y="{1540+42*i}" text-anchor="middle" font-size="36" fill="{GRN2}">{t}</text>')
A('<line x1="100" y1="1180" x2="100" y2="1640" stroke="#666" stroke-width="4"/><polygon points="100,1640 88,1612 112,1612" fill="#666"/>')
A('<text x="100" y="1700" text-anchor="middle" font-size="38" fill="#666">Toward</text><text x="100" y="1744" text-anchor="middle" font-size="38" fill="#666">Earth</text>')

# ---------- legend ----------
LY = CASB + 70
rows = [
 (RED, "Casing Support Tracks (2)", "Two loops of 12 mm HTS tape per track, 0.5 m wide each, in persistent mode. They carry the casing."),
 (NAVY, "Load Compensation Windings (2)", "Driven 3 mm HTS tape on the support track loops. They add lift under a slow or parked sled."),
 (NAVY, "Lift Compensation Tracks (2)", "One driven loop 1 m wide per track, on the casing ceiling. They pull the casing down as a fast sled passes."),
 (NAVY, "Gap Protection", "The driven windings also guard the gap. If it closes on one face, the windings on the other face pull back."),
 (GRN, "LIM Stators", "Deployment thrust and lateral centering. Narrower than the reaction plates (edge effect)."),
 (ORG, "LIM Reaction Plates", "γ-TiAl alloy, 200 mm thick, 1.2 m wide, one on each side. Offset above center with the LIM stators."),
 (OLV, "Ferromagnetic Plates (4)", "Sintered iron in non-conductive resin, 1.4 m wide and 200 mm thick. One for each track."),
]
A(f'<rect x="140" y="{LY}" width="2520" height="{110+70*len(rows)+30}" rx="10" fill="#f8f9f9" stroke="#bbb" stroke-width="2"/>')
A(f'<text x="1400" y="{LY+70}" text-anchor="middle" font-size="44" font-weight="bold">Component Summary</text>')
for i, (col, name, desc) in enumerate(rows):
    yy = LY + 126 + 70*i
    A(f'<rect x="200" y="{yy-24}" width="44" height="30" fill="{col}"/>')
    A(f'<text x="266" y="{yy}" font-size="31" fill="#333"><tspan font-weight="bold" fill="{col}">{name}:</tspan> {desc}</text>')
A('</svg>')
OUT = ring.out_dir("figures_ch08")   # under ORBITAL_RING_OUT (default runs)
SVG = _os.path.join(OUT, "fig10_bearing_cross_section.svg")
open(SVG, "w", encoding="utf-8").write("\n".join(out) + "\n")
try:
    import cairosvg
    cairosvg.svg2png(url=SVG, write_to=_os.path.join(OUT, "fig10_bearing_cross_section.png"), output_width=5600)
    cairosvg.svg2png(url=SVG, write_to=_os.path.join(OUT, "preview.png"), output_width=1500)
except ImportError:
    print("cairosvg is not installed: the SVG was written, no PNG")
print(LY+110+70*len(rows)+30, H)
