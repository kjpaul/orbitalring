#!/usr/bin/env python3
"""
Working drawings of the bearing's coils for round 5 (2026-11-13), Section 1 of Vol_III/Claude_Code_instructions_bearing_round5_2026-11-13.md.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 8 (Figure 8.4).  These are Paul's working drawings, not the book's figures.
Run: python tools/round5_drawings.py [--out=FOLDER] [--sup=W,T] [--fast=W,T]
     writes six drawings as SVG and PNG into Vol_III/bearing_drawings_2026-11-13/ (default), from the numbers the round 4 and round 5 runs use:
       01 the bearing in cross-section (to scale), 02 a fast track's bundle (to scale), 03 a fast coil in plan (length compressed),
       04 a support track's bundle (to scale), 05 a support coil in plan (length compressed), 06 the two circuit diagrams.
Environment: none.  The plates default to round 4's provisional 1.2 m x 150 mm and are redrawn with --sup/--fast when the optimum is known.

The words, as Paul ruled them: a tape; a layer (one tape in a soldered stack); a turn (a stack of soldered layers, once around the loop);
a coil (its series turns); the coils that overlap in a bundle (20 fast coils in 20 staggered families, 10 support coils staggered by 100 m).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Polygon
import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be
import bearing_optimize as bo
import final_bearing as fb

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "Vol_III", "bearing_drawings_2026-11-13")
SUP_PLATE = (1.2, 0.150)                  # m, round 4's provisional plates
FAST_PLATE = (1.2, 0.150)
SIDE = ring.CABLE_SIDE                    # 7.39 m
GAP = bd.GAP0                             # 0.100 m
W_SUP_LOOP = bo.BEST_SUPPORT["W"]         # 0.5 m per loop, two loops
W_FAST_LOOP = bo.BEST_LIFT["W"]           # 0.8 m
B_SUP = (0.150, 0.0667)                   # support bundle, w x h
B_FAST = (0.200, 0.050)                   # fast bundle
L_FAST, N_FAM = 25.0, 20
L_SUP, N_SUP = 1000.0, 10
TAPE_PITCH = bd.TAPE_T / be.F_CABLE       # 0.08 mm / 0.75 = 0.107 mm per layer in a soldered stack
FS = 12                                   # the smallest label, points

plt.rcParams.update({"font.size": FS, "font.family": "DejaVu Sans", "svg.fonttype": "none"})

BLUE, RED, GREY, IRON, COPPER, GREEN, ORANGE = "#2266aa", "#cc4400", "#666666", "#9a8a4a", "#b87333", "#228833", "#e08a00"


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=150 if ext == "png" else None, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.svg / .png")


def rect(ax, x, y, w, h, fc, ec="black", lw=1.0, **kw):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=lw, **kw))


def dim(ax, x1, y1, x2, y2, text, offset=0.0, fs=FS, color=GREY, side="above"):
    """A dimension line with arrow heads and its text."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="<->", color=color, lw=0.9, shrinkA=0, shrinkB=0))
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    if x1 == x2:
        ax.text(xm + offset, ym, text, fontsize=fs, color=color, ha="left" if offset >= 0 else "right", va="center")
    else:
        ax.text(xm, ym + offset, text, fontsize=fs, color=color, ha="center", va="bottom" if offset >= 0 else "top")


# ------------------------------------------------------------------------------------------------------------------------------
# 01  the bearing in cross-section, to scale
# ------------------------------------------------------------------------------------------------------------------------------
def track_positions(sup_plate, fast_plate):
    """Plate centres from the cable's centreline, bearing_design.part_F's rule: the load plates at the bottom face's outer edges, the support
    plates inboard with 0.2 m between plate edges, the lift plates at the top face's outer edges.  Plate widths as given (not loop + 0.4)."""
    half = SIDE / 2
    x_load = half - fast_plate[0] / 2
    x_sup = half - fast_plate[0] - 0.2 - sup_plate[0] / 2
    x_lift = half - fast_plate[0] / 2
    return x_sup, x_load, x_lift


def draw_cross_section(sup_plate=SUP_PLATE, fast_plate=FAST_PLATE, provisional=True):
    x_sup, x_load, x_lift = track_positions(sup_plate, fast_plate)
    fig, ax = plt.subplots(figsize=(13, 12))
    half = SIDE / 2
    # the cable
    rect(ax, -half, -half, SIDE, SIDE, "#dbe8f3", BLUE, 1.5)
    ax.text(0, 1.6, "CNT cable", fontsize=18, ha="center", color="#1b4a73", weight="bold")
    ax.text(0, 1.1, f"{SIDE:.2f} m square, {ring.M_CABLE_M:,.0f} kg/m with hardware", fontsize=FS + 1, ha="center", color="#1b4a73")
    ax.text(0, 0.6, "moving into the page at 8,847 m/s relative to the casing", fontsize=FS, ha="center", color="#1b4a73")
    # EMS plates on the cable's top and bottom faces
    tS, tF = sup_plate[1], fast_plate[1]
    for sx in (-1, 1):
        rect(ax, sx * x_sup - sup_plate[0] / 2, -half - tS, sup_plate[0], tS, IRON, "#5c5020")
        rect(ax, sx * x_load - fast_plate[0] / 2, -half - tF, fast_plate[0], tF, IRON, "#5c5020")
        rect(ax, sx * x_lift - fast_plate[0] / 2, half, fast_plate[0], tF, IRON, "#5c5020")
    # LIM reaction plates on the side faces (1.2 m tall, 200 mm thick, offset above centre as Figure 8.4 has them) and the stators outside them
    for sx in (-1, 1):
        x0 = sx * half if sx > 0 else -half - 0.2
        rect(ax, x0, 0.3, 0.2, 1.2, "#c8b273", "#5c5020")
        xs = half + 0.2 + 0.2 if sx > 0 else -half - 0.2 - 0.2 - 0.3
        rect(ax, xs, 0.4, 0.3, 1.0, "#cfe6d2", GREEN)
        ax.text(sx * (half + 0.1), 1.62, "LIM plate\n1.2 m, 200 mm", fontsize=FS, ha="center", va="bottom", color="#5c5020")
        ax.text(sx * (half + 0.55), 1.52, "LIM\nstator", fontsize=FS, ha="center", va="bottom", color=GREEN)
    # the bundles across the 100 mm gap: support (four 150 x 67, the centre pair touching) and load under the cable, lift over it
    def sup_bundles(xc, y_face, below=True):
        wb, hb = B_SUP
        xs = [xc - W_SUP_LOOP - wb / 2, xc - wb / 2, xc + wb / 2, xc + W_SUP_LOOP - wb / 2]
        for x in xs:
            y = (y_face - GAP - hb) if below else (y_face + GAP)
            rect(ax, x, y, wb, hb, "#cfe0f0", BLUE, 1.0)
    def fast_bundles(xc, y_face, below=True):
        wb, hb = B_FAST
        for x in (xc - W_FAST_LOOP / 2 - wb / 2, xc + W_FAST_LOOP / 2 - wb / 2):
            y = (y_face - GAP - hb) if below else (y_face + GAP)
            rect(ax, x, y, wb, hb, "#f6d5c3", RED, 1.0)
    y_bot = -half - tS; y_botF = -half - tF; y_top = half + tF
    for sx in (-1, 1):
        sup_bundles(sx * x_sup, y_bot, True)
        fast_bundles(sx * x_load, y_botF, True)
        fast_bundles(sx * x_lift, y_top, False)
    # the casing: floor under the lower bundles, ceiling over the upper, walls outside the LIM stators (schematic thickness)
    y_floor = y_bot - GAP - B_SUP[1] - 0.05
    y_ceil = y_top + GAP + B_FAST[1] + 0.05
    x_wall = half + 0.2 + 0.2 + 0.3 + 0.15
    rect(ax, -x_wall - 0.08, y_floor - 0.08, 2 * x_wall + 0.16, 0.08, "#e8e8e8", "#444")
    rect(ax, -x_wall - 0.08, y_ceil, 2 * x_wall + 0.16, 0.08, "#e8e8e8", "#444")
    rect(ax, -x_wall - 0.08, y_floor - 0.08, 0.08, y_ceil - y_floor + 0.16, "#e8e8e8", "#444")
    rect(ax, x_wall, y_floor - 0.08, 0.08, y_ceil - y_floor + 0.16, "#e8e8e8", "#444")
    prov = " (provisional)" if provisional else ""
    # labels: the lift tracks above the casing, the load and support tracks below it in two rows, each with a leader line to its bundles
    def leader(x0, y0, x1, y1, color):
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.8)
    for sx, lab in ((-1, "Lift track L"), (1, "Lift track R")):
        ax.text(sx * x_lift, y_ceil + 0.22, lab, fontsize=FS + 1, ha="center", va="bottom", color=RED, weight="bold")
        leader(sx * x_lift, y_ceil + 0.2, sx * x_lift, y_top + GAP + B_FAST[1] / 2, RED)
    for sx, lab in ((-1, "Load track L"), (1, "Load track R")):
        ax.text(sx * x_load, y_floor - 0.25, lab, fontsize=FS + 1, ha="center", va="top", color=RED, weight="bold")
        leader(sx * x_load, y_floor - 0.25, sx * x_load, y_botF - GAP - B_FAST[1] / 2, RED)
    for sx, lab in ((-1, "Support track L"), (1, "Support track R")):
        ax.text(sx * x_sup, y_floor - 0.62, lab, fontsize=FS + 1, ha="center", va="top", color=BLUE, weight="bold")
        leader(sx * x_sup, y_floor - 0.62, sx * x_sup, y_bot - GAP - B_SUP[1] / 2, BLUE)
    ax.text(0, y_ceil + 0.62, "casing, 12,000 kg/m (walls schematic; the bundles' cryostats and mounts not drawn).  Lift tracks: bundles on the casing's ceiling over the cable's top plates.",
            fontsize=FS, ha="center", va="bottom", color="#444")
    ax.text(0, y_floor - 1.5, "Support tracks: two 0.5 m loops, four 150 x 67 mm bundles of 12 mm tape, on the casing's floor under the cable's bottom plates;  "
                              "load tracks: one 0.8 m loop, two 200 x 50 mm bundles of 3 mm tape, beside them.", fontsize=FS - 1, ha="center", va="top", color="#444")
    ax.text(0, y_floor - 1.75, f"Plates on the cable: six, all {sup_plate[0]:.1f} m x {sup_plate[1]*1e3:.0f} mm{prov} here (support) and {fast_plate[0]:.1f} m x {fast_plate[1]*1e3:.0f} mm{prov} (lift and load); "
                              "0.2 m between a support plate and the load plate beside it; the load and lift plates at the cable's edges.", fontsize=FS - 1, ha="center", va="top", color="#5c5020")
    # dimensions
    dim(ax, -half, y_floor - 2.2, half, y_floor - 2.2, f"{SIDE:.2f} m", offset=-0.04)
    dim(ax, x_wall + 0.6, -half, x_wall + 0.6, half, f"{SIDE:.2f} m", offset=0.06)
    dim(ax, -x_sup, y_floor - 1.1, x_sup, y_floor - 1.1, f"{2*x_sup:.2f} m between the support tracks' centres", offset=-0.03)
    dim(ax, -x_load, y_ceil + 0.42, x_load, y_ceil + 0.42, f"{2*x_load:.2f} m between the lift tracks' centres (and the load tracks')", offset=0.02)
    xg = -x_lift - fast_plate[0] / 2 - 0.08
    dim(ax, xg, half + tF, xg, half + tF + GAP, "100 mm gap", offset=-0.06, fs=FS)
    # a detail of the right support and load tracks at 4 x, drawn inside the cable's outline
    k = 4.0; x_c = (x_sup + x_load) / 2; y_c = -half - tS - GAP / 2
    X0, Y0 = -1.4, -1.75                                                   # the detail's centre on the page
    def M(x, y):
        return X0 + k * (x - x_c), Y0 + k * (y - y_c)
    def rect_m(x, y, w, h, fc, ec, lw=1.0):
        X, Y = M(x, y); rect(ax, X, Y, k * w, k * h, fc, ec, lw)
    span = (x_sup - sup_plate[0] / 2 - 0.15, x_load + fast_plate[0] / 2 + 0.15)
    X1, Y1 = M(span[0], y_floor - 0.02); X2, Y2 = M(span[1], -half + 0.03)
    rect(ax, X1, Y1, X2 - X1, Y2 - Y1, "white", "#444", 1.0)
    rect_m(span[0], -half, span[1] - span[0], 0.03, "#dbe8f3", BLUE, 1.0)
    rect_m(x_sup - sup_plate[0] / 2, -half - tS, sup_plate[0], tS, IRON, "#5c5020")
    rect_m(x_load - fast_plate[0] / 2, -half - tF, fast_plate[0], tF, IRON, "#5c5020")
    for x in (x_sup - W_SUP_LOOP - B_SUP[0] / 2, x_sup - B_SUP[0] / 2, x_sup + B_SUP[0] / 2, x_sup + W_SUP_LOOP - B_SUP[0] / 2):
        rect_m(x, y_bot - GAP - B_SUP[1], B_SUP[0], B_SUP[1], "#cfe0f0", BLUE, 1.0)
    for x in (x_load - W_FAST_LOOP / 2 - B_FAST[0] / 2, x_load + W_FAST_LOOP / 2 - B_FAST[0] / 2):
        rect_m(x, y_botF - GAP - B_FAST[1], B_FAST[0], B_FAST[1], "#f6d5c3", RED, 1.0)
    rect_m(span[0], y_floor - 0.02, span[1] - span[0], 0.02, "#e8e8e8", "#444")
    ax.text((X1 + X2) / 2, Y2 + 0.08, "detail at 4 x: support track R and load track R (cable above, casing floor below)", fontsize=FS, ha="center", va="bottom", color="#444")
    Xa, Ya = M(x_sup - sup_plate[0] / 2, -half - tS); Xb, _ = M(x_sup + sup_plate[0] / 2, 0)
    dim(ax, Xa, Ya - 0.12, Xb, Ya - 0.12, f"support plate {sup_plate[0]:.1f} m x {sup_plate[1]*1e3:.0f} mm{prov}", offset=-0.04, fs=FS - 1)
    Xa, Ya = M(x_load - fast_plate[0] / 2, -half - tF); Xb, _ = M(x_load + fast_plate[0] / 2, 0)
    dim(ax, Xa, Ya - 0.12, Xb, Ya - 0.12, f"load plate {fast_plate[0]:.1f} m x {fast_plate[1]*1e3:.0f} mm{prov}", offset=-0.04, fs=FS - 1)
    Xa, Ya = M(x_sup + sup_plate[0] / 2, -half - tS / 2); Xb, _ = M(x_load - fast_plate[0] / 2, 0)
    dim(ax, Xa, Ya, Xb, Ya, "0.2 m", offset=0.05, fs=FS - 1)
    Xa, Ya = M(x_load + fast_plate[0] / 2 + 0.05, -half - tF); _, Yb = M(0, -half - tF - GAP)
    dim(ax, Xa, Ya, Xa, Yb, "100 mm gap", offset=0.05, fs=FS - 1)
    Xa, Ya = M(x_sup - W_SUP_LOOP, y_bot - GAP - B_SUP[1] - 0.04); Xb, _ = M(x_sup + W_SUP_LOOP, 0)
    dim(ax, Xa, Ya, Xb, Ya, "two 0.5 m loops: four 150 x 67 mm bundles, the centre pair touching", offset=-0.04, fs=FS - 1)
    Xa, Ya = M(x_load - W_FAST_LOOP / 2, y_botF - GAP - B_FAST[1] - 0.04); Xb, _ = M(x_load + W_FAST_LOOP / 2, 0)
    dim(ax, Xa, Ya, Xb, Ya, "0.8 m loop: two 200 x 50 mm bundles", offset=-0.04, fs=FS - 1)
    ax.plot([x_sup - sup_plate[0] / 2 - 0.15, X1], [y_floor, Y1], color="#444", lw=0.6, ls=":"); ax.plot([x_load + fast_plate[0] / 2 + 0.15, X2], [y_floor, Y1], color="#444", lw=0.6, ls=":")
    ax.set_xlim(-x_wall - 0.5, x_wall + 1.3); ax.set_ylim(y_floor - 2.5, y_ceil + 1.0)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Drawing 1: the bearing in cross-section, to scale (the casing walls and the LIM stators schematic; the detail box is the same geometry at 4 x).\n"
                 "Six tracks: two support and two load tracks under the cable, two lift tracks over it; each track's bundles face their own plate across the 100 mm gap.",
                 fontsize=FS + 1, loc="left")
    save(fig, "01_bearing_cross_section")


# ------------------------------------------------------------------------------------------------------------------------------
# 02  a fast track's bundle in cross-section, to scale
# ------------------------------------------------------------------------------------------------------------------------------
def fast_cells(N_turns, E_kVmm, tau=10e-3, n_fam=N_FAM, r_lift=None):
    r_sup, r_lift = be.load_records()
    Phi = 2 * r_lift["W_req"] / r_lift["NI_req"]; V_turn = Phi * L_FAST / tau; V_fam = N_turns * V_turn
    NI_coil = r_lift["NI_req"] / n_fam; tapes = NI_coil / (0.5 * r_lift["ic3"]); layers = tapes / N_turns
    A_fam = tapes * 0.003 * bd.TAPE_T; E_coil = r_lift["W_req"] * L_FAST / n_fam
    lay = be.bundle_layout(B_FAST[0], B_FAST[1], n_fam, A_fam, N_turns, V_turn, V_fam, E_kVmm, 10, 2)
    return dict(V_turn=V_turn, V_fam=V_fam, I=NI_coil / N_turns, tapes=tapes, layers=layers, A_fam=A_fam, lay=lay, E_coil=E_coil)


def draw_bundle_panel(ax, wb, hb, lay, n_x, n_y, N_turns, layers, w_tape, title, turn_color, h_draw=None):
    """One bundle to scale: the family cells, the walls, and inside one cell the series turns as stacks of soldered layers."""
    h_draw = h_draw or hb
    t_g, t_fam, t_turn, c_w, c_h = lay["t_g"], lay["t_fam"], lay["t_turn"], lay["c_w"], lay["c_h"]
    rect(ax, 0, 0, wb, h_draw, "#f3f3f3", "black", 1.5)                    # the bundle's wall to ground (the cryostat's cold wall)
    rect(ax, t_g, t_g, wb - 2 * t_g, h_draw - 2 * t_g, "#ffffff", GREY, 0.8)
    for j in range(n_y):
        for i in range(n_x):
            x = t_g + i * (c_w + t_fam); y = t_g + j * (c_h + t_fam)
            rect(ax, x, y, c_w, c_h, "#eef3f8", GREY, 0.6)
    # inside the first cell: the turns as stacks, w_tape wide and layers x pitch thick, laid in rows with the turn insulation between them
    x0 = t_g; y0 = t_g
    stack_t = layers * TAPE_PITCH
    # columns of stacks along the cell's width, rows along its height
    n_col = max(1, int((c_w - 2e-4) // (w_tape + t_turn)))
    n_row = int(math.ceil(N_turns / n_col))
    pitch_x = (c_w - 1e-3) / n_col; pitch_y = stack_t + t_turn + 0.3e-3
    k = 0
    for r in range(n_row):
        for c in range(n_col):
            if k >= N_turns:
                break
            x = x0 + 0.5e-3 + c * pitch_x; y = y0 + 0.5e-3 + r * pitch_y
            rect(ax, x, y, w_tape, stack_t, turn_color, "black", 0.5)
            k += 1
    ax.set_xlim(-0.004, wb + 0.004); ax.set_ylim(-0.004, h_draw + 0.004); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=FS, loc="left")
    return dict(n_col=n_col, n_row=n_row, stack_t=stack_t)


def draw_fast_bundle():
    c10 = fast_cells(10, 5.0); c14a = fast_cells(14, 10.0); c14b = fast_cells(14, 5.0)
    wb, hb = B_FAST
    fig, axs = plt.subplots(3, 1, figsize=(13, 13.5))
    l10 = c10["lay"]
    s = draw_bundle_panel(axs[0], wb, hb, l10, 10, 2, 10, c10["layers"], 0.003,
        f"Panel A: 200 x 50 mm bundle, 20 family cells (10 x 2), each the cross-section of one 25 m coil; 10 series turns of {c10['layers']:.0f} soldered layers of 3 mm tape.\n"
        f"At 10 ms a turn sees {c10['V_turn']/1e3:.1f} kV and the terminals {c10['V_fam']/1e3:.1f} kV: at 5 kV/mm the wall between cells and to ground is {l10['t_fam']*1e3:.2f} mm, the film between turns "
        f"{l10['t_turn']*1e3:.2f} mm; the cells {l10['c_w']*1e3:.1f} x {l10['c_h']*1e3:.1f} mm; height needed {l10['h_need']*1e3:.1f} mm of the 50: fits.", "#9fc5e8")
    l14 = c14a["lay"]
    draw_bundle_panel(axs[1], wb, hb, l14, 10, 2, 14, c14a["layers"], 0.003,
        f"Panel B: the same bundle with 14 series turns of {c14a['layers']:.0f} layers (Paul's 1,500 A case, {c14a['I']:,.0f} A per turn), terminals {c14a['V_fam']/1e3:.1f} kV at 10 ms, insulation at 10 kV/mm:\n"
        f"walls {l14['t_fam']*1e3:.2f} mm, turn film {l14['t_turn']*1e3:.2f} mm, cells {l14['c_w']*1e3:.1f} x {l14['c_h']*1e3:.1f} mm, height needed {l14['h_need']*1e3:.1f} mm of the 50: fits at 10 kV/mm.", "#f4b183")
    l14b = c14b["lay"]
    draw_bundle_panel(axs[2], wb, l14b["h_need"], l14b, 10, 2, 14, c14b["layers"], 0.003,
        f"Panel C: 14 turns at 20 kV with the insulation kept at 5 kV/mm: walls {l14b['t_fam']*1e3:.2f} mm, the bundle must be {l14b['h_need']*1e3:.0f} mm tall instead of 50 (drawn at that height).\n"
        f"In every panel the first cell (lower left) shows its turns as stacks of soldered layers, {s['stack_t']*1e3:.1f} mm thick for 70 layers at a 0.107 mm layer pitch (0.08 mm tape at a 0.75 fill); "
        f"the other cells are the other 19 overlapping coils.", "#f4b183", h_draw=l14b["h_need"])
    # a scale bar on the first panel
    axs[0].plot([0.150, 0.200], [-0.0025, -0.0025], color="black", lw=2); axs[0].text(0.175, -0.0038, "50 mm", ha="center", va="top", fontsize=FS)
    axs[0].text(0.0, -0.0038, "coolant: the 15 percent of each cell's area left open (not drawn separately); the wall to ground is the cryostat's cold wall", ha="left", va="top", fontsize=FS - 1, color=GREY)
    fig.suptitle("Drawing 2: a fast track's bundle (lift or load) in cross-section, to scale.  20 coils overlap in every bundle; each cell is one coil.", fontsize=FS + 2, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "02_fast_bundle_cross_section")
    return dict(c10=c10, c14a=c14a, c14b=c14b)


# ------------------------------------------------------------------------------------------------------------------------------
# 03  a fast coil in plan
# ------------------------------------------------------------------------------------------------------------------------------
def draw_fast_coil_plan():
    c10 = fast_cells(10, 5.0)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw=dict(height_ratios=(1.4, 1)))
    # the loop: 25 m by 0.8 m, drawn with the length compressed 5 x (x in m / 5)
    cx = 5.0; W = W_FAST_LOOP; L = L_FAST / cx
    for k in range(10):                                                   # the 10 series turns, drawn as nested rectangles
        d = 0.018 * k
        rect(ax, 0 + d, -W / 2 + d, L - 2 * d, W - 2 * d, "none", BLUE, 1.2)
    ax.text(L / 2, 0, f"one 25 m coil: a 0.8 m loop of 10 series turns\n(each turn a stack of {c10['layers']:.0f} soldered layers of 3 mm tape; drawn as 10 nested lines)",
            ha="center", va="center", fontsize=FS + 1, color=BLUE)
    # the two ends at the left, the leads out to the warm switch and store
    ax.plot([0.0, -0.6], [-W / 2 + 0.0, -0.9], color=COPPER, lw=2.5); ax.plot([0.0 + 0.162, -0.5], [-W / 2 + 0.162, -0.9], color=COPPER, lw=2.5)
    rect(ax, -1.5, -1.45, 1.1, 0.5, "#fff2cc", ORANGE, 1.2)
    ax.text(-0.95, -1.2, "warm: switches\nand store", ha="center", va="center", fontsize=FS, color="#7a4a00")
    ax.text(-0.55, -0.72, "the two ends of the coil;\npulse-sized leads out of the cold", ha="center", va="center", fontsize=FS - 1, color=COPPER)
    ax.text(L / 2, W / 2 + 0.12, "length compressed 5 x: the loop is 25 m by 0.8 m, drawn 5 m by 0.8 m", ha="center", fontsize=FS, color=GREY)
    dim(ax, 0, -W / 2 - 0.12, L, -W / 2 - 0.12, "25 m (compressed)", offset=-0.03)
    dim(ax, L + 0.12, -W / 2, L + 0.12, W / 2, "0.8 m", offset=0.05)
    ax.text(L / 2, W / 2 + 0.26, f"{c10['I']:,.0f} A per turn at full pull, {c10['V_turn']/1e3:.1f} kV per turn and {c10['V_fam']/1e3:.1f} kV at the terminals at 10 ms; {c10['E_coil']/1e6:.3f} MJ stored", ha="center", fontsize=FS, color=GREY)
    ax.set_xlim(-1.7, L + 0.9); ax.set_ylim(-1.6, 0.95); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Drawing 3: a fast coil (lift or load) in plan.  Not to scale along the track: the length is compressed 5 x; the loop's width is to that same scale.", fontsize=FS + 1, loc="left")
    # the 20 staggered coils of a track: each 25 m, the next one 1.25 m along
    cx2 = 5.0; L2 = L_FAST / cx2; s = L_FAST / N_FAM / cx2
    for k in range(N_FAM):
        y = 0.9 - k * 0.09
        ax2.plot([k * s, k * s + L2], [y, y], color=BLUE if k % 2 == 0 else "#88aacc", lw=3, solid_capstyle="butt")
    for k in range(N_FAM, N_FAM + 8):
        y = 0.9 - (k - N_FAM) * 0.09
        ax2.plot([k * s, k * s + L2], [y, y], color="#cccccc", lw=3, solid_capstyle="butt")
    ax2.text(L2 / 2, 1.05, "the 20 coils (families) of one track, each 25 m, staggered by 25 / 20 = 1.25 m; every point of the track lies under 20 coils", ha="left", fontsize=FS, color=BLUE)
    ax2.text(N_FAM * s + L2 / 2, 0.98 - 0 * 0.09, "the next 20 coils", ha="left", fontsize=FS, color="#888888")
    dim(ax2, 0, -1.0, L2, -1.0, "25 m (compressed 5 x)", offset=-0.04)
    dim(ax2, 0, -0.85, s, -0.85, "1.25 m", offset=-0.04)
    ax2.plot([L2 * 0.6, L2 * 0.6], [-0.95, 1.0], color=RED, lw=1.0, ls="--"); ax2.text(L2 * 0.6 + 0.05, -0.75, "one point of the track:\n20 coils over it", fontsize=FS - 1, color=RED, va="bottom")
    ax2.set_xlim(-0.2, (N_FAM + 8) * s + L2 + 0.2); ax2.set_ylim(-1.15, 1.2); ax2.axis("off")
    ax2.set_title("The stagger of the 20 coils of a track (the same 5 x compression; each line is one coil; the rows are only to separate them, the coils share one bundle)", fontsize=FS, loc="left")
    fig.tight_layout()
    save(fig, "03_fast_coil_plan")


# ------------------------------------------------------------------------------------------------------------------------------
# 04  a support track's bundle in cross-section, to scale
# ------------------------------------------------------------------------------------------------------------------------------
def support_cells(E_kVmm=5.0, V_charge=110.0, V_dump=1000.0):
    r_sup, r_lift = be.load_records()
    turns_leg = r_sup["turns12"] / 4                                     # Part 3b's turns12 counts the four bundles: 2,753 tapes through one bundle (an outer leg, or one of the centre pair)
    tapes_bundle = turns_leg                                             # 2,753 tapes of 12 mm per bundle, shared among the 10 coils
    tapes_coil = tapes_bundle / N_SUP                                    # 275 per coil: 27.5 turns of 10 layers
    A_coil = tapes_coil * 0.012 * bd.TAPE_T
    # cells: 10 coils as 5 x 2 cells; the wall between cells and to ground for the 1 kV dump, the film between turns for the 110 V charge
    lay = be.bundle_layout(B_SUP[0], B_SUP[1], N_SUP, A_coil, 28, V_charge, V_dump, E_kVmm, 5, 2)
    return dict(tapes_bundle=tapes_bundle, tapes_coil=tapes_coil, A_coil=A_coil, lay=lay, turns_coil=tapes_coil / 10)


def draw_support_bundle():
    sc = support_cells()
    lay = sc["lay"]; wb, hb = B_SUP
    fig, ax = plt.subplots(figsize=(13, 7.5))
    s = draw_bundle_panel(ax, wb, hb, lay, 5, 2, 28, 10, 0.012,
        f"Drawing 4: a support track's bundle, 150 x 67 mm, to scale.  10 cells (5 x 2), one per coil of the 10 staggered 1 km coils that overlap here; each cell holds that coil's\n"
        f"{sc['turns_coil']:.1f} turns through this leg (55 series turns around the two loops, 27 or 28 through each leg), each turn 10 soldered layers of 12 mm tape ({10*TAPE_PITCH*1e3:.2f} mm thick).\n"
        f"Walls between cells and to ground {lay['t_fam']*1e3:.2f} mm for the 1 kV of the dump at 5 kV/mm; the film between turns {lay['t_turn']*1e3:.3f} mm (110 V while charging, a 50 um film minimum);\n"
        f"cells {lay['c_w']*1e3:.1f} x {lay['c_h']*1e3:.1f} mm, height needed {lay['h_need']*1e3:.1f} of the 67 mm.  {sc['tapes_bundle']:,.0f} tapes through the bundle in all.", "#a9d18e")
    ax.plot([0.100, 0.150], [-0.0025, -0.0025], color="black", lw=2); ax.text(0.125, -0.0038, "50 mm", ha="center", va="top", fontsize=FS)
    ax.text(0.0, -0.0038, "the first cell (lower left) shows its 28 turns as 12 mm x 1.07 mm stacks; the other cells are the other nine coils; coolant is the 15 percent of each cell left open",
            ha="left", va="top", fontsize=FS - 1, color=GREY)
    fig.tight_layout()
    save(fig, "04_support_bundle_cross_section")
    return sc


# ------------------------------------------------------------------------------------------------------------------------------
# 05  a support coil in plan
# ------------------------------------------------------------------------------------------------------------------------------
def draw_support_coil_plan():
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(14, 10.5), gridspec_kw=dict(height_ratios=(1.3, 1)))
    cx = 200.0; L = L_SUP / cx                                           # 1 km drawn as 5 m: compressed 200 x
    W = W_SUP_LOOP
    # the figure-eight: loop 1 (left lobe in y) and loop 2 (right lobe), sharing the centre leg; drawn as two rectangles 0.5 m wide stacked in y
    for k in range(6):                                                    # a few of the 55 turns, nested
        d = 0.012 * k
        rect(ax, d, 0 + d, L - 2 * d, W - 2 * d, "none", BLUE, 1.1)
        rect(ax, d, -W + d, L - 2 * d, W - 2 * d, "none", BLUE, 1.1)
    ax.annotate("", xy=(L * 0.55, W * 0.3), xytext=(L * 0.45, W * 0.3), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=(L * 0.45, -W * 0.3), xytext=(L * 0.55, -W * 0.3), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=(L * 0.55, 0.0), xytext=(L * 0.45, 0.0), arrowprops=dict(arrowstyle="->", color=RED, lw=2.5))
    ax.text(L * 0.5, W * 0.62, "loop 1 (0.5 m wide)", ha="center", va="center", fontsize=FS, color=BLUE)
    ax.text(L * 0.5, -W * 0.62, "loop 2 (0.5 m wide), the opposite sense", ha="center", va="center", fontsize=FS, color=BLUE)
    ax.text(L * 0.5, -W - 0.62, "red arrow: the centre leg, both loops' turns in the same direction (two bundles side by side)", ha="center", va="top", fontsize=FS - 1, color=RED)
    ax.text(L * 0.5, W + 0.22, "one 1 km support coil: 55 series turns around the two loops in series (a figure-eight),\neach turn 10 soldered layers of 12 mm tape; six turns drawn",
            ha="center", va="bottom", fontsize=FS, color=BLUE)
    ax.text(L * 0.5, W + 0.72, "length compressed 200 x: the coil is 1,000 m by 1.0 m, drawn 5 m by 1.0 m", ha="center", va="bottom", fontsize=FS, color=GREY)
    dim(ax, 0, -W - 0.16, L, -W - 0.16, "1,000 m (compressed 200 x)", offset=-0.04)
    dim(ax, L + 0.15, -W, L + 0.15, W, "2 x 0.5 m", offset=0.05)
    # the coil's end: at the left it bends down (out of the plane, toward the casing floor) to its diode, strip, dump switch and resistor and lead contacts
    ax.text(-2.4, W + 0.22, "the coil's end, in the cold\n(bent down out of the loop's plane):", ha="center", va="bottom", fontsize=FS - 1, color=GREY)
    items = (("free-wheel diode D2", "#444"), ("strip and its heater", RED), ("dump switch and resistor", "#444"), ("lead contacts (to the warm supply)", COPPER))
    for i, (lab, col) in enumerate(items):
        y = W - 0.05 - i * 0.42
        rect(ax, -1.05, y - 0.3, 0.35, 0.3, "white", col, 1.2)
        ax.text(-1.15, y - 0.15, lab, ha="right", va="center", fontsize=FS - 1, color=col)
    ax.plot([0.0, -0.4, -0.4, -0.7], [W * 0.1, W * 0.1, -W - 0.2, -W - 0.2], color=COPPER, lw=2)
    ax.plot([0.012, -0.3, -0.3, -0.7], [-W * 0.1, -W * 0.1, -W - 0.1, -W - 0.1], color=COPPER, lw=2)
    ax.set_xlim(-4.3, L + 1.0); ax.set_ylim(-W - 1.0, W + 1.1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Drawing 5: a support coil in plan.  Not to scale along the track (200 x compressed); the loop widths are to the scale of the drawn length.", fontsize=FS + 1, loc="left")
    # the stagger of the ten coils in a bundle: 1 km each, 100 m apart
    s = (L_SUP / N_SUP) / cx; L2 = L
    for k in range(N_SUP):
        y = 0.9 - k * 0.17
        ax2.plot([k * s, k * s + L2], [y, y], color=BLUE if k % 2 == 0 else "#88aacc", lw=4, solid_capstyle="butt")
        ax2.plot([k * s, k * s], [y - 0.05, y + 0.05], color="black", lw=1)
    for k in range(N_SUP, N_SUP + 4):
        y = 0.9 - (k - N_SUP) * 0.17
        ax2.plot([k * s, k * s + L2], [y, y], color="#cccccc", lw=4, solid_capstyle="butt")
    ax2.text(0, 1.1, "the 10 coils of one bundle, each 1 km, staggered by 100 m: every point of the track is inside 10 coils; a coil's end (tick) lies in a 100 m cell where only its 9 neighbours run on",
             ha="left", fontsize=FS, color=BLUE)
    dim(ax2, 0, -0.85, L2, -0.85, f"1 km (compressed {cx:.0f} x)", offset=-0.04)
    dim(ax2, 0, -0.7, s, -0.7, "100 m", offset=-0.04)
    ax2.set_xlim(-0.2, (N_SUP + 4) * s + L2 + 0.2); ax2.set_ylim(-1.0, 1.25); ax2.axis("off")
    ax2.set_title("The stagger of the ten coils of a bundle (each line is one coil; the rows only separate them, the coils share the bundle)", fontsize=FS, loc="left")
    fig.tight_layout()
    save(fig, "05_support_coil_plan")


# ------------------------------------------------------------------------------------------------------------------------------
# 06  the two circuit diagrams
# ------------------------------------------------------------------------------------------------------------------------------
def wire(ax, pts, color="black", lw=1.6):
    xs, ys = zip(*pts); ax.plot(xs, ys, color=color, lw=lw, solid_joinstyle="miter")


def box(ax, x, y, w, h, label, fc="white", ec="black", fs=FS, color="black"):
    rect(ax, x, y, w, h, fc, ec, 1.4); ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fs, color=color)


def coil_symbol(ax, x, y, n=6, r=0.12, color=BLUE, vertical=True):
    t = np.linspace(0, math.pi, 30)
    for k in range(n):
        if vertical:
            ax.plot(x + r * np.sin(t), y + k * 2 * r + r - r * np.cos(t), color=color, lw=2)
        else:
            ax.plot(x + k * 2 * r + r - r * np.cos(t), y + r * np.sin(t), color=color, lw=2)


def diode(ax, x, y, up=True, color="black"):
    """A diode symbol on a vertical wire at x, centred on y; conducting upward when up."""
    s = 0.12
    tri = Polygon([(x - s, y - s if up else y + s), (x + s, y - s if up else y + s), (x, y + s if up else y - s)], closed=True, facecolor=color, edgecolor=color)
    ax.add_patch(tri); ax.plot([x - s, x + s], [y + s if up else y - s] * 2, color=color, lw=2)


def switch(ax, x, y, closed=False, color="black", label=None, horizontal=True):
    if horizontal:
        ax.plot([x, x + 0.1], [y, y], color=color, lw=1.6); ax.plot([x + 0.4, x + 0.5], [y, y], color=color, lw=1.6)
        ax.plot([x + 0.1, x + 0.4], [y, y if closed else y + 0.18], color=color, lw=1.8)
        ax.plot([x + 0.1], [y], "o", color=color, ms=4); ax.plot([x + 0.4], [y], "o", color=color, ms=4)
    else:
        ax.plot([x, x], [y, y + 0.1], color=color, lw=1.6); ax.plot([x, x], [y + 0.4, y + 0.5], color=color, lw=1.6)
        ax.plot([x, x if closed else x + 0.18], [y + 0.1, y + 0.4], color=color, lw=1.8)
    if label:
        ax.text(x + 0.25 if horizontal else x + 0.25, y + 0.3 if horizontal else y + 0.25, label, ha="center" if horizontal else "left", va="bottom" if horizontal else "center", fontsize=FS - 1, color=color)


def draw_circuits():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    # ---- the support coil ----
    ax.set_title("Drawing 6a: the support coil's circuit (one 1 km coil; each of the ten coils of a bundle has its own).\nWarm above the dashed line, 70 K below.  Not to scale.", fontsize=FS, loc="left")
    ax.axhline(3.2, color=GREY, ls="--", lw=1); ax.text(6.9, 3.28, "warm", fontsize=FS, color=GREY, ha="right"); ax.text(6.9, 3.0, "cold, 70 K", fontsize=FS, color=GREY, ha="right", va="top")
    box(ax, 0.3, 3.6, 2.0, 1.0, "DC supply\nfrom the 100 kV bus\n(current-regulated,\nabout 110 V at 1 %/s)", fs=FS - 1)
    box(ax, 2.9, 3.6, 1.6, 1.0, "lead contacts\n(withdrawn or\nopened between\ncharges)", fs=FS - 1, color=COPPER, ec=COPPER)
    wire(ax, [(2.3, 4.3), (2.9, 4.3)]); wire(ax, [(2.3, 3.9), (2.9, 3.9)])
    # leads down into the cold
    wire(ax, [(3.3, 3.6), (3.3, 2.4)], COPPER, 2.2); wire(ax, [(4.1, 3.6), (4.1, 2.4)], COPPER, 2.2)
    ax.text(3.7, 2.9, "a pair of copper leads", ha="center", fontsize=FS - 1, color=COPPER)
    # the cold node rails
    wire(ax, [(1.0, 2.4), (6.2, 2.4)]); wire(ax, [(1.0, 0.4), (6.2, 0.4)])
    # the coil (right)
    coil_symbol(ax, 6.2, 0.55, n=7, r=0.12, color=BLUE); wire(ax, [(6.2, 0.4), (6.2, 0.55)], BLUE); wire(ax, [(6.2, 2.23), (6.2, 2.4)], BLUE)
    ax.text(6.4, 1.4, "the coil:\n55 turns of\n10 layers,\n21 H (ten\ntogether),\n512 A normal,\n1,173 A cap", ha="left", va="center", fontsize=FS - 1, color=BLUE)
    # D2 free-wheel diode across the coil
    wire(ax, [(5.4, 0.4), (5.4, 2.4)]); diode(ax, 5.4, 1.4, up=True)
    ax.text(5.45, 0.5, "D2, free-\nwheel diode,\n1.0 V at 70 K", ha="left", va="bottom", fontsize=FS - 1)
    # the strip with heater across D2 (and so across the coil)
    wire(ax, [(4.6, 0.4), (4.6, 2.4)], RED, 1.6); rect(ax, 4.45, 1.1, 0.3, 0.6, "#ffd9cc", RED, 1.4)
    ax.text(4.6, 1.05, "strip: bare REBCO\ntape, held normal\nby its heater\nwhile charging", ha="center", va="top", fontsize=FS - 1, color=RED)
    rect(ax, 4.8, 1.15, 0.22, 0.5, "white", ORANGE, 1.2); ax.text(4.75, 1.75, "heater", ha="left", va="bottom", fontsize=FS - 1, color=ORANGE)
    # the dump branch: cold switch and 20 ohm resistor
    wire(ax, [(2.3, 0.4), (2.3, 2.4)]); switch(ax, 2.3, 1.55, closed=False, horizontal=False, label="dump\nswitch")
    rect(ax, 2.2, 0.7, 0.2, 0.6, "white", "black", 1.4); ax.text(2.45, 1.0, "dump\nresistor", ha="left", va="center", fontsize=FS - 1)
    ax.text(0.95, 1.4, "the dump branch:\nclosed only to\ndischarge a\nquenched coil;\nthe resistor sized\nfor 1 kV at the\ncoil's current\n(2 ohm at 512 A)", ha="left", va="center", fontsize=FS - 1)
    ax.text(0.3, 0.05, "Normal hold: leads open, strip superconducting (heater off), D2 carries nothing, dump open: the current circulates coil -> strip with no loss.\n"
                       "Charge: heater on (strip resistive), leads closed, supply drives the coil at 110 V; D2 blocks.  Close (Section 8): heater off, strip cools and takes the current.",
            fontsize=FS - 1, va="top", ha="left")
    ax.set_xlim(0, 8.3); ax.set_ylim(-0.9, 4.8); ax.axis("off")
    # ---- the fast coil ----
    ax2.set_title("Drawing 6b: a fast coil's circuit (one 25 m lift or load coil; 20 per 25 m of track).\nWarm above the dashed line, 70 K below.  Not to scale.", fontsize=FS, loc="left")
    ax2.axhline(2.6, color=GREY, ls="--", lw=1); ax2.text(7.9, 2.68, "warm", fontsize=FS, color=GREY, ha="right"); ax2.text(7.9, 2.4, "cold, 70 K", fontsize=FS, color=GREY, ha="right", va="top")
    box(ax2, 0.3, 3.0, 1.7, 1.2, "charger from\nthe 100 kV bus\n(refills the store's\nlosses)", fs=FS - 1)
    box(ax2, 2.4, 3.0, 1.8, 1.2, "store C\n0.150 MJ at full pull\n(capacitor bank, or\nthe pulse bus of 3.3)", fs=FS - 1, fc="#fff2cc", ec=ORANGE, color="#7a4a00")
    box(ax2, 4.6, 3.0, 2.0, 1.2, "H-bridge of pulse\nswitches: charges\nthe coil in 10 ms,\nholds, returns the\nenergy to C", fs=FS - 1, fc="#e8f1e8", ec=GREEN, color="#1e5c2a")
    wire(ax2, [(2.0, 3.8), (2.4, 3.8)]); wire(ax2, [(2.0, 3.4), (2.4, 3.4)]); wire(ax2, [(4.2, 3.8), (4.6, 3.8)]); wire(ax2, [(4.2, 3.4), (4.6, 3.4)])
    wire(ax2, [(5.2, 3.0), (5.2, 1.6)], COPPER, 2.4); wire(ax2, [(6.0, 3.0), (6.0, 1.6)], COPPER, 2.4)
    ax2.text(5.6, 2.1, "pulse-sized leads\n(2,107 A for 10 turns,\n1,505 A for 14)", ha="center", va="center", fontsize=FS - 1, color=COPPER)
    coil_symbol(ax2, 5.2, 1.6, n=4, r=0.1, color=RED, vertical=False)
    ax2.text(5.6, 1.15, "the coil: 10 (or 14) series turns, 14.3 (20) kV at the terminals at 10 ms", ha="center", va="top", fontsize=FS - 1, color=RED)
    ax2.text(0.3, 0.6, "The switch must block the terminal voltage (14.3 kV at 10 turns, 20 kV at 14), carry the turn current for the pass (2,107 or 1,505 A for up to 46 ms\n"
                       "at 30 km/s, 10 ms rise and fall), and return the coil's energy to the store: a four-quadrant bridge, not a thyristor crowbar (Section 4).\n"
                       "No persistent switch, no diode: the coil holds nothing between pulses.  A small local store keeps the gap control alive when the bus is lost (Section 3.4).",
             fontsize=FS - 1, va="top", ha="left")
    ax2.set_xlim(0, 8.3); ax2.set_ylim(-0.6, 4.6); ax2.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.02, wspace=0.04)
    save(fig, "06_circuits")


def main():
    global SUP_PLATE, FAST_PLATE, OUT
    provisional = True
    for a in sys.argv[1:]:
        if a.startswith("--out="): OUT = a[6:]
        elif a.startswith("--sup="): SUP_PLATE = tuple(float(v) for v in a[6:].split(",")); provisional = False
        elif a.startswith("--fast="): FAST_PLATE = tuple(float(v) for v in a[7:].split(",")); provisional = False
        elif a in ("-h", "--help"): print(__doc__); return
    print(f"Drawings into {OUT}: support plate {SUP_PLATE[0]:.2f} m x {SUP_PLATE[1]*1e3:.0f} mm, fast plate {FAST_PLATE[0]:.2f} m x {FAST_PLATE[1]*1e3:.0f} mm{' (provisional)' if provisional else ''}")
    draw_cross_section(SUP_PLATE, FAST_PLATE, provisional)
    fb_ = draw_fast_bundle()
    draw_fast_coil_plan()
    sc = draw_support_bundle()
    draw_support_coil_plan()
    draw_circuits()
    c10, c14a, c14b = fb_["c10"], fb_["c14a"], fb_["c14b"]
    print(f"  fast bundle: 10 turns x {c10['layers']:.0f} layers, {c10['I']:,.0f} A, {c10['V_fam']/1e3:.1f} kV: cells {c10['lay']['c_w']*1e3:.1f} x {c10['lay']['c_h']*1e3:.1f} mm, walls {c10['lay']['t_fam']*1e3:.2f} mm, turn film {c10['lay']['t_turn']*1e3:.2f} mm, height {c10['lay']['h_need']*1e3:.1f} mm")
    print(f"               14 turns x {c14a['layers']:.0f} layers, {c14a['I']:,.0f} A, {c14a['V_fam']/1e3:.1f} kV at 10 kV/mm: height {c14a['lay']['h_need']*1e3:.1f} mm; at 5 kV/mm {c14b['lay']['h_need']*1e3:.1f} mm")
    print(f"  support bundle: {sc['tapes_bundle']:,.0f} tapes, {sc['tapes_coil']:.0f} per coil ({sc['turns_coil']:.1f} turns of 10 layers through this leg); cells {sc['lay']['c_w']*1e3:.1f} x {sc['lay']['c_h']*1e3:.1f} mm, height {sc['lay']['h_need']*1e3:.1f} of 67 mm")


if __name__ == "__main__":
    main()
