#!/usr/bin/env python3
"""
The bearing's working drawings, round 6 part A (2026-11-14): the layout rules of the brief in one place (RULES below), the round 5 optimum drawn
("round 5 optimum, not yet accepted" on every sheet), the cells of the bundles laid out stack by stack.
Book: Orbital Ring Engineering, Volume III: Chapter 7 (the bearing), Chapter 8 (the circuits).
Run: python tools/round6_drawings.py [--out=FOLDER] [--only=1a,2,...]
     writes SVG and PNG into Vol_III/bearing_drawings_2026-11-14/ (default): 01a_bearing_cross_section, 01b_underside_right_detail, 01c_lift_track_detail,
     02_fast_bundle_cross_section, 03_fast_coil_plan, 04_support_bundle_cross_section, 05_support_coil_plan, 06a_support_coil_circuit, 06b_fast_coil_circuit.
Reads runs/bearing_round5_final.json (the optimum plates and their solves).  Environment: ORBITAL_RING_OUT (default runs).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import json
import math
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Arc
import numpy as np

import ring_altitude as ring
import bearing_design as bd
import bearing_electrical as be

# ------------------------------------------------------------------------------------------------------------------------------
# RULES: the brief's layout rules, in one place
# ------------------------------------------------------------------------------------------------------------------------------
RULES = dict(
    font="DejaVu Sans",                 # one sans-serif font
    fs_title=16, fs_label=12, fs_note=11,
    lw_outline=1.6, lw_dim=0.8,         # two line weights
    width_in=12.0,                      # every sheet 12 inches wide; the height follows the content
    margin=0.04,                        # the even margin around the content, as a fraction of the content's larger extent
    note_chars=150,                     # the notes wrap at this many characters (11 pt on a 12 in sheet)
    caption_chars=135,                  # a panel caption wraps at this many characters (12 pt bold)
    status="round 5 optimum, not yet accepted",
)
COL = dict(support="#163f73", fast="#e58a2a", plate="#5d4b1e", plate_fill="#eadfb8", lim="#3f7f32", lim_fill="#d3e7cc", copper="#9c4a1e",
           cable="#dde9f3", cable_edge="#3a5f85", grey="#666666", black="#111111", white="#ffffff", coolant="#f3f3f3")
HATCH = dict(plate="//", lim="xx")
LEGEND_ITEMS = dict(support=("support track (persistent 12 mm tape)", COL["support"], None), fast=("lift and load tracks (driven 3 mm tape)", COL["fast"], None),
                    plate=("iron plates on the cable", COL["plate_fill"], HATCH["plate"]), lim=("LIM reaction plate and stator", COL["lim_fill"], HATCH["lim"]),
                    copper=("copper leads", COL["copper"], None))
plt.rcParams["font.family"] = RULES["font"]
plt.rcParams["hatch.linewidth"] = 0.6

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "Vol_III", "bearing_drawings_2026-11-14")
FINAL = os.path.join(os.environ.get("ORBITAL_RING_OUT", "runs"), "bearing_round5_final.json")
SIDE = ring.CABLE_SIDE                   # 7.39 m
GAP = bd.GAP0                            # 0.100 m
B_SUP = (0.150, 0.0667)                  # support bundle w x h
B_FAST = (0.200, 0.050)                  # fast bundle w x h
L_FAST, N_FAM = 25.0, 20
L_SUP, N_SUP = 1000.0, 10
D_SPACING = 0.2                          # round 5's spacing between a support plate and the load plate beside it
TAPE_PITCH = bd.TAPE_T / be.F_CABLE      # 0.107 mm per soldered layer
TAU = 10e-3
F_NORMAL_SHARE = 47.5 / 245.5
L_COMMON_55 = 21.3
FILM_MIN = 50e-6                         # the thinnest film between turns [assumption, as bearing_electrical.bundle_layout]
V_REL_OPT = 8862.0                       # m/s past the casing with the optimum's hardware (the round 5 numbers sheet)
M_CABLE_OPT = 101879.0                   # kg/m with hardware, the same


# ------------------------------------------------------------------------------------------------------------------------------
# the design numbers from the round 5 optimum
# ------------------------------------------------------------------------------------------------------------------------------
def design():
    fin = json.load(open(FINAL, encoding="utf-8"))
    rs, rl = fin["results"]["support:opt135"], fin["results"]["fast:opt135"]
    d = dict(rs=rs, rl=rl, hw=fin["hw"])
    NI_coil = rl["NI_req"] / N_FAM; tapes = NI_coil / (0.5 * rl["ic_req"]); Phi = 2 * rl["W_mag"] / rl["NI_req"]
    d["fast"] = {}
    for N in (10, 14):
        V_turn = Phi * L_FAST / TAU
        d["fast"][N] = dict(N=N, I=NI_coil / N, layers=tapes / N, V_turn=V_turn, V=N * V_turn, E=rl["W_mag"] * L_FAST / N_FAM, t_stack=tapes / N * TAPE_PITCH)
    tapes_bundle = rs["turns_leg"]; tapes_coil = tapes_bundle / N_SUP; turns_leg = tapes_coil / 10
    I_cap = rs["NI_req"] / (N_SUP * turns_leg); I_norm = I_cap * math.sqrt(F_NORMAL_SHARE)
    d["sup"] = dict(tapes_bundle=tapes_bundle, tapes_coil=tapes_coil, turns_leg=turns_leg, turns=2 * turns_leg, I_cap=I_cap, I_norm=I_norm,
                    L=L_COMMON_55 * (2 * turns_leg / 55) ** 2, n1=math.ceil(turns_leg), n2=math.floor(turns_leg))
    d["x_load"] = SIDE / 2 - rl["Wp"] / 2; d["x_lift"] = d["x_load"]; d["x_sup"] = SIDE / 2 - rl["Wp"] - D_SPACING - rs["Wp"] / 2
    return d


def cell_layout(N_turns, layers, V, V_turn, E_kVmm, w_tape, bundle_w, n_cx, n_cy, films_min=FILM_MIN):
    """Whole stacks in a cell: n_x stacks (w_tape wide, layers x pitch thick, faces horizontal) across, rows up, the film t_turn between stacks,
    walls t_fam between cells and to ground.  Returns the cell and bundle sizes (m)."""
    E = E_kVmm * 1e6
    t_fam = max(V / E, films_min); t_turn = max(V_turn / E, films_min) if N_turns > 1 else 0.0
    t_s = layers * TAPE_PITCH
    w_c = (bundle_w - (n_cx + 1) * t_fam) / n_cx
    n_x = max(1, int((w_c + t_turn) // (w_tape + t_turn)))
    n_x = min(n_x, N_turns)
    n_y = math.ceil(N_turns / n_x)
    w_used = n_x * w_tape + (n_x - 1) * t_turn
    h_c = n_y * t_s + (n_y - 1) * t_turn
    h_need = n_cy * h_c + (n_cy + 1) * t_fam
    return dict(t_fam=t_fam, t_turn=t_turn, t_s=t_s, w_c=w_c, n_x=n_x, n_y=n_y, w_used=w_used, h_c=h_c, h_need=h_need, w_tape=w_tape, N=N_turns, layers=layers)


# ------------------------------------------------------------------------------------------------------------------------------
# the sheet: title, panels with equal aspect, dimensions, labels with leaders, legend, notes, scale bar
# ------------------------------------------------------------------------------------------------------------------------------
def measure(fig, ax):
    """The union of the axes' artists' extents in data units (texts included), after a draw."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    xs, ys = [], []
    for a in list(ax.patches) + list(ax.lines) + list(ax.texts):
        try:
            bb = a.get_window_extent(r)
        except Exception:
            continue
        if bb.width == 0 and bb.height == 0:
            continue
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        xs += [x0, x1]; ys += [y0, y1]
    return min(xs), max(xs), min(ys), max(ys)


class Panel:
    def __init__(self, draw_fn, caption=None):
        self.draw_fn = draw_fn; self.caption = caption


class Sheet:
    def __init__(self, title, name, scale_note, legend=(), notes=()):
        self.title = title; self.name = name; self.scale_note = scale_note; self.legend = legend; self.notes = list(notes); self.panels = []
        self.dpu = 1.0

    def panel(self, draw_fn, caption=None):
        self.panels.append(Panel(draw_fn, caption)); return self

    def text_w(self, s, fs):
        return len(s) * fs / 72.0 * 0.62 * self.dpu

    def text_h(self, fs, lines=1):
        return lines * fs / 72.0 * 1.3 * self.dpu

    def dim(self, ax, p1, p2, offset, text, fs=None, color=None):
        """An engineering dimension between p1 and p2 (horizontal or vertical): extension lines, a dimension line with arrowheads at the signed offset,
        the figure centered beside the line when it fits along it, else beside the line at its middle (horizontal: to the right of the far end)."""
        fs = fs or RULES["fs_label"]; color = color or COL["black"]; lw = RULES["lw_dim"]
        (x1, y1), (x2, y2) = p1, p2
        ext = 0.03 * abs(offset) if offset else 0
        th = self.text_h(fs)
        if abs(y2 - y1) < 1e-12:
            yd = y1 + offset
            for x in (x1, x2):
                ax.plot([x, x], [y1 + (0.12 * offset if offset else 0), yd + math.copysign(ext, offset)], color=color, lw=lw)
            ax.annotate("", xy=(x2, yd), xytext=(x1, yd), arrowprops=dict(arrowstyle="<|-|>", color=color, lw=lw, mutation_scale=9, shrinkA=0, shrinkB=0))
            if self.text_w(text, fs) < 0.9 * abs(x2 - x1):
                ax.text((x1 + x2) / 2, yd + (0.25 * th if offset >= 0 else -0.25 * th), text, fontsize=fs, color=color, ha="center", va="bottom" if offset >= 0 else "top")
            else:
                ax.text(max(x1, x2) + 0.35 * th, yd, text, fontsize=fs, color=color, ha="left", va="center")
        else:
            xd = x1 + offset
            for y in (y1, y2):
                ax.plot([x1 + (0.12 * offset if offset else 0), xd + math.copysign(ext, offset)], [y, y], color=color, lw=lw)
            ax.annotate("", xy=(xd, y2), xytext=(xd, y1), arrowprops=dict(arrowstyle="<|-|>", color=color, lw=lw, mutation_scale=9, shrinkA=0, shrinkB=0))
            if self.text_w(text, fs) < 0.9 * abs(y2 - y1):
                ax.text(xd + (0.25 * th if offset >= 0 else -0.25 * th), (y1 + y2) / 2, text, fontsize=fs, color=color, ha="left" if offset >= 0 else "right", va="center", rotation=90)
            else:
                ax.text(xd + (0.35 * th if offset >= 0 else -0.35 * th), (y1 + y2) / 2, text, fontsize=fs, color=color, ha="left" if offset >= 0 else "right", va="center")

    def label(self, ax, target, anchor, text, color=None, ha="left", va="center", fs=None, leader=True):
        fs = fs or RULES["fs_label"]; color = color or COL["black"]
        ax.text(anchor[0], anchor[1], text, fontsize=fs, color=color, ha=ha, va=va)
        if leader:
            ax.plot([target[0], anchor[0]], [target[1], anchor[1]], color=color, lw=RULES["lw_dim"])
            ax.plot([target[0]], [target[1]], marker="o", ms=2.5, color=color)

    def rect(self, ax, x, y, w, h, fc, ec, lw=None, hatch=None, z=2):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, lw=lw if lw is not None else RULES["lw_outline"], hatch=hatch, zorder=z))

    def scale_bar(self, ax, x, y, length, text):
        th = self.text_h(RULES["fs_label"])
        ax.plot([x, x + length], [y, y], color=COL["black"], lw=2.2, solid_capstyle="butt")
        for xx in (x, x + length):
            ax.plot([xx, xx], [y - 0.12 * th, y + 0.12 * th], color=COL["black"], lw=1.2)
        ax.text(x + length / 2, y - 0.35 * th, text, fontsize=RULES["fs_label"], ha="center", va="top")

    def render(self, out):
        W = RULES["width_in"]; fs_t, fs_l, fs_n = RULES["fs_title"], RULES["fs_label"], RULES["fs_note"]
        title_lines = textwrap.wrap(self.title + ".  " + RULES["status"] + ".", 95)
        note_lines = []
        for n in [self.scale_note] + self.notes:
            note_lines += textwrap.wrap(n, RULES["note_chars"]) + [""]
        if note_lines and note_lines[-1] == "":
            note_lines.pop()
        h_title = 0.32 * len(title_lines) + 0.25
        h_notes = 0.19 * len(note_lines) + 0.25
        leg_rows = []
        if self.legend:
            row = []; x = 0.03
            for key in self.legend:
                w = 0.02 + len(LEGEND_ITEMS[key][0]) * fs_n / 72.0 * 0.6 / W + 0.025
                if row and x + w > 0.98:
                    leg_rows.append(row); row = []; x = 0.03
                row.append(key); x += w
            leg_rows.append(row)
        h_legend = 0.3 * len(leg_rows) + (0.08 if leg_rows else 0.0)
        caps = [textwrap.wrap(p.caption, RULES["caption_chars"]) if p.caption else [] for p in self.panels]
        # pass 1: measure each panel on a scratch figure at a provisional scale
        bboxes = []
        for p in self.panels:
            fig = plt.figure(figsize=(W, W)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_aspect("equal"); ax.axis("off")
            self.dpu = 1.0; p.draw_fn(ax); x0, x1, y0, y1 = measure(fig, ax)
            self.dpu = (x1 - x0) / W
            ax.cla(); ax.set_aspect("equal"); ax.axis("off"); p.draw_fn(ax); x0, x1, y0, y1 = measure(fig, ax)
            plt.close(fig)
            bboxes.append([x0, x1, y0, y1])
        # the real figure, rebuilt until every panel's content (texts at their final size) lies inside its box
        for it in range(4):
            heights = []
            for (x0, x1, y0, y1), cl in zip(bboxes, caps):
                m = RULES["margin"] * max(x1 - x0, y1 - y0)
                heights.append((W * (y1 - y0 + 2 * m) / (x1 - x0 + 2 * m), 0.22 * len(cl) + (0.1 if cl else 0.0), m))
            H = h_title + sum(h + c for h, c, _ in heights) + h_legend + h_notes
            fig = plt.figure(figsize=(W, H))
            fig.text(0.03, 1 - 0.18 / H, "\n".join(title_lines), fontsize=fs_t, ha="left", va="top", linespacing=1.25)
            y_top = 1 - h_title / H
            axes = []; grow = False
            for p, (x0, x1, y0, y1), (h_ax, h_cap, m), cl in zip(self.panels, bboxes, heights, caps):
                if cl:
                    fig.text(0.03, y_top - 0.05 / H, "\n".join(cl), fontsize=fs_l, ha="left", va="top", weight="bold", linespacing=1.25)
                    y_top -= h_cap / H
                ax = fig.add_axes([0, y_top - h_ax / H, 1, h_ax / H]); ax.set_aspect("equal"); ax.axis("off")
                self.dpu = (x1 - x0 + 2 * m) / W
                p.draw_fn(ax); ax.set_xlim(x0 - m, x1 + m); ax.set_ylim(y0 - m, y1 + m)
                axes.append(ax)
                y_top -= h_ax / H
            # measure at the final scale
            for k, ax in enumerate(axes):
                mx0, mx1, my0, my1 = measure(fig, ax)
                x0, x1, y0, y1 = bboxes[k]
                tol = 0.002 * max(x1 - x0, y1 - y0)
                nb = [min(x0, mx0), max(x1, mx1), min(y0, my0), max(y1, my1)]
                if any(abs(a - b) > tol for a, b in zip(nb, bboxes[k])):
                    bboxes[k] = nb; grow = True
            if not grow or it == 3:
                break
            plt.close(fig)
        for row in leg_rows:
            x = 0.03
            for key in row:
                lab, c, hatch = LEGEND_ITEMS[key]
                fig.patches.append(Rectangle((x, y_top - 0.2 / H), 0.016, 0.13 / H, transform=fig.transFigure, facecolor=c,
                                             edgecolor=COL["plate"] if key == "plate" else (COL["lim"] if key == "lim" else c), hatch=hatch, lw=1.0))
                fig.text(x + 0.02, y_top - 0.135 / H, lab, fontsize=fs_n, ha="left", va="center")
                x += 0.02 + len(lab) * fs_n / 72.0 * 0.6 / W + 0.025
            y_top -= 0.3 / H
        if leg_rows:
            y_top -= 0.08 / H
        fig.text(0.03, y_top - 0.1 / H, "\n".join(note_lines), fontsize=fs_n, ha="left", va="top", linespacing=1.3)
        os.makedirs(out, exist_ok=True)
        for ext in ("svg", "png"):
            fig.savefig(os.path.join(out, f"{self.name}.{ext}"), dpi=150, facecolor="white")
        plt.close(fig)
        print(f"  wrote {self.name}.svg / .png ({W:.0f} x {H:.1f} in)")


# ------------------------------------------------------------------------------------------------------------------------------
# 1a  the bearing in cross-section
# ------------------------------------------------------------------------------------------------------------------------------
def sheet_1a(D):
    rs, rl = D["rs"], D["rl"]; half = SIDE / 2; tS, tF = rs["T"], rl["T"]
    x_sup, x_load, x_lift = D["x_sup"], D["x_load"], D["x_lift"]
    y_floor = -half - tS - GAP - B_SUP[1] - 0.25; y_ceil = half + tF + GAP + B_FAST[1] + 0.25
    x_wall = half + 1.25

    def draw(ax):
        s = sh
        s.rect(ax, -half, -half, SIDE, SIDE, COL["cable"], COL["cable_edge"])
        ax.text(0, 0.6, "CNT cable", fontsize=RULES["fs_title"], ha="center", va="center", color=COL["cable_edge"], weight="bold")
        ax.text(0, 0.0, f"{SIDE:.2f} m square, {M_CABLE_OPT:,.0f} kg/m with hardware", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["cable_edge"])
        ax.text(0, -0.5, f"into the page at {V_REL_OPT:,.0f} m/s past the casing", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["cable_edge"])
        for sx in (-1, 1):
            s.rect(ax, sx * x_sup - rs["Wp"] / 2, -half - tS, rs["Wp"], tS, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
            s.rect(ax, sx * x_load - rl["Wp"] / 2, -half - tF, rl["Wp"], tF, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
            s.rect(ax, sx * x_lift - rl["Wp"] / 2, half, rl["Wp"], tF, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
            xp = half if sx > 0 else -half - 0.2
            s.rect(ax, xp, 0.3, 0.2, 1.2, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
            xs = half + 0.4 if sx > 0 else -half - 0.7
            s.rect(ax, xs, 0.4, 0.3, 1.0, COL["lim_fill"], COL["lim"], hatch=HATCH["lim"])
            for x in (sx * x_sup - rs["Wl"] - B_SUP[0] / 2, sx * x_sup - B_SUP[0] / 2, sx * x_sup + B_SUP[0] / 2, sx * x_sup + rs["Wl"] - B_SUP[0] / 2):
                s.rect(ax, x, -half - tS - GAP - B_SUP[1], B_SUP[0], B_SUP[1], COL["support"], COL["support"], lw=0.8)
            for x in (sx * x_load - rl["Wl"] / 2 - B_FAST[0] / 2, sx * x_load + rl["Wl"] / 2 - B_FAST[0] / 2):
                s.rect(ax, x, -half - tF - GAP - B_FAST[1], B_FAST[0], B_FAST[1], COL["fast"], COL["fast"], lw=0.8)
                s.rect(ax, x, half + tF + GAP, B_FAST[0], B_FAST[1], COL["fast"], COL["fast"], lw=0.8)
        s.rect(ax, -x_wall, y_floor - 0.08, 2 * x_wall, 0.08, COL["white"], COL["grey"], lw=1.0)
        s.rect(ax, -x_wall, y_ceil, 2 * x_wall, 0.08, COL["white"], COL["grey"], lw=1.0)
        s.rect(ax, -x_wall - 0.08, y_floor - 0.08, 0.08, y_ceil - y_floor + 0.16, COL["white"], COL["grey"], lw=1.0)
        s.rect(ax, x_wall, y_floor - 0.08, 0.08, y_ceil - y_floor + 0.16, COL["white"], COL["grey"], lw=1.0)
        ax.text(0, y_ceil - 0.17, "casing, 12,000 kg/m (schematic)", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["grey"])
        xa = x_wall + 0.45
        s.label(ax, (x_lift, half + tF + GAP + B_FAST[1] / 2), (xa, half + 1.0), "lift track (one of two)", COL["fast"])
        s.label(ax, (half + 0.1, 1.5), (xa, 2.5), "LIM reaction plate, 1.2 m by 200 mm", COL["plate"])
        s.label(ax, (half + 0.55, 1.4), (xa, 1.75), "LIM stator (schematic)", COL["lim"])
        s.label(ax, (x_load, -half - tF - GAP - B_FAST[1] / 2), (xa, -half - 1.0), "load track (one of two)", COL["fast"])
        s.label(ax, (x_sup + rs["Wl"] - B_SUP[0] / 2, -half - tS - GAP - B_SUP[1] / 2), (xa, -half - 1.9), "support track (one of two)", COL["support"])
        s.dim(ax, (-half, y_floor - 0.08), (half, y_floor - 0.08), -1.0, f"{SIDE:.2f} m")
        s.dim(ax, (x_wall + 0.08, -half), (x_wall + 0.08, half), 0.3, f"{SIDE:.2f} m")
        s.dim(ax, (-x_sup, y_floor - 0.08), (x_sup, y_floor - 0.08), -0.5, f"{2*x_sup:.2f} m, support tracks")
        s.dim(ax, (-x_lift, y_ceil + 0.08), (x_lift, y_ceil + 0.08), 0.5, f"{2*x_lift:.2f} m, lift and load tracks")
        s.scale_bar(ax, -x_wall, y_floor - 1.9, 1.0, "1 m")
    sh = Sheet("Drawing 1a: the bearing in cross-section, the whole", "01a_bearing_cross_section",
               "To scale (1 m bar): the cable, the plates, the bundles and the gap; the casing's walls, floor and ceiling, the LIM stator and the labels are schematic.",
               legend=("support", "fast", "plate", "lim"),
               notes=[f"Six tracks. Two support tracks on the casing floor under the cable's underside, each two {rs['Wl']:.1f} m loops of 12 mm tape in four 150 by 67 mm bundles facing a "
                      f"{rs['Wp']:.1f} m by {rs['T']*1e3:.0f} mm plate; two load tracks beside them at the underside's edges and two lift tracks on the ceiling over the top face's edges, "
                      f"each one {rl['Wl']:.1f} m loop of 3 mm tape in two 200 by 50 mm bundles facing a {rl['Wp']:.1f} m by {rl['T']*1e3:.0f} mm plate. {D_SPACING:.1f} m between a support plate's edge and "
                      f"the load plate's beside it (round 5; round 6a's Section 3 asks what it should be). The 100 mm gap is between a plate and its bundles. Drawings 1b and 1c enlarge the right-hand underside and one lift track.",
                      f"The cable is drawn {SIDE:.2f} m square (the round 5 geometry; with the optimum's hardware {ring.M_CABLE_M:,.0f} becomes {M_CABLE_OPT:,.0f} kg/m and the side 7.41 m). Cable hardware "
                      f"{D['hw']:,.0f} kg/m (the LIM plates 1,872, two support plates {rs['plate_m']:,.0f}, four lift and load plates {rl['plate_m']:,.0f}); tape {2*rs['tape']+4*rl['tape']:.0f} kg per meter of casing. "
                      "Cryostats, mounts and the LIM stator's real shape are not drawn."])
    sh.panel(draw)
    return sh


# ------------------------------------------------------------------------------------------------------------------------------
# 1b  the right-hand underside, enlarged; 1c  one lift track, enlarged
# ------------------------------------------------------------------------------------------------------------------------------
def sheet_1b(D):
    rs, rl = D["rs"], D["rl"]; half = SIDE / 2; tS, tF = rs["T"], rl["T"]
    x_sup, x_load = D["x_sup"], D["x_load"]
    y_cab = -half; y_floor = y_cab - tS - GAP - B_SUP[1] - 0.06
    xl0 = x_sup - rs["Wp"] / 2 - 0.2; xl1 = half + 0.05

    def draw(ax):
        s = sh
        s.rect(ax, xl0, y_cab, xl1 - xl0, 0.1, COL["cable"], COL["cable_edge"])
        ax.text((x_sup + x_load) / 2, y_cab + 0.05, "the cable's underside", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["cable_edge"])
        s.rect(ax, x_sup - rs["Wp"] / 2, y_cab - tS, rs["Wp"], tS, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
        s.rect(ax, x_load - rl["Wp"] / 2, y_cab - tF, rl["Wp"], tF, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
        xs_sup = [x_sup - rs["Wl"] - B_SUP[0] / 2, x_sup - B_SUP[0] / 2, x_sup + B_SUP[0] / 2, x_sup + rs["Wl"] - B_SUP[0] / 2]
        yb_s = y_cab - tS - GAP - B_SUP[1]
        for x in xs_sup:
            s.rect(ax, x, yb_s, B_SUP[0], B_SUP[1], COL["support"], COL["support"], lw=0.8)
        xs_load = [x_load - rl["Wl"] / 2 - B_FAST[0] / 2, x_load + rl["Wl"] / 2 - B_FAST[0] / 2]
        yb_l = y_cab - tF - GAP - B_FAST[1]
        for x in xs_load:
            s.rect(ax, x, yb_l, B_FAST[0], B_FAST[1], COL["fast"], COL["fast"], lw=0.8)
        s.rect(ax, xl0, y_floor - 0.04, xl1 - xl0, 0.04, COL["white"], COL["grey"], lw=1.0)
        # above the cable band: the plate widths, the spacing above them, the cable edge above that
        s.dim(ax, (x_sup - rs["Wp"] / 2, y_cab + 0.1), (x_sup + rs["Wp"] / 2, y_cab + 0.1), 0.2, f"support plate {rs['Wp']*1e3:.0f}")
        s.dim(ax, (x_load - rl["Wp"] / 2, y_cab + 0.1), (x_load + rl["Wp"] / 2, y_cab + 0.1), 0.2, f"load plate {rl['Wp']*1e3:.0f}")
        s.dim(ax, (x_sup + rs["Wp"] / 2, y_cab + 0.1), (x_load - rl["Wp"] / 2, y_cab + 0.1), 0.48, f"{D_SPACING*1e3:.0f} between the plates")
        s.label(ax, (half, y_cab + 0.1), (half - 0.05, y_cab + 0.68), "cable edge: the load plate's outer edge is flush with it", COL["cable_edge"], ha="right")
        # the thicknesses beside the plates, the gap between the load bundles
        s.dim(ax, (x_sup - rs["Wp"] / 2, y_cab - tS), (x_sup - rs["Wp"] / 2, y_cab), -0.08, f"{tS*1e3:.0f}")
        s.dim(ax, (half, y_cab - tF), (half, y_cab), 0.08, f"{tF*1e3:.0f}")
        s.dim(ax, (x_load, y_cab - tF - GAP), (x_load, y_cab - tF), 0.0, f"gap {GAP*1e3:.0f}")
        # the bundles: width under the first of each kind, height beside it
        s.dim(ax, (xs_sup[0], yb_s), (xs_sup[0] + B_SUP[0], yb_s), -0.05, f"{B_SUP[0]*1e3:.0f}")
        s.dim(ax, (xs_sup[0], yb_s), (xs_sup[0], yb_s + B_SUP[1]), -0.06, f"{B_SUP[1]*1e3:.0f}")
        s.dim(ax, (xs_load[1], yb_l), (xs_load[1] + B_FAST[0], yb_l), -0.05, f"{B_FAST[0]*1e3:.0f}")
        s.dim(ax, (xs_load[1] + B_FAST[0], yb_l), (xs_load[1] + B_FAST[0], yb_l + B_FAST[1]), 0.06, f"{B_FAST[1]*1e3:.0f}")
        # under the floor: the loops (leg center to leg center), then the labels with leaders, in rows
        s.dim(ax, (xs_sup[0] + B_SUP[0] / 2, y_floor - 0.04), (xs_sup[1] + B_SUP[0] / 2, y_floor - 0.04), -0.14, f"loop {rs['Wl']*1e3:.0f}")
        s.dim(ax, (xs_sup[2] + B_SUP[0] / 2, y_floor - 0.04), (xs_sup[3] + B_SUP[0] / 2, y_floor - 0.04), -0.14, f"loop {rs['Wl']*1e3:.0f}")
        s.dim(ax, (xs_load[0] + B_FAST[0] / 2, y_floor - 0.04), (xs_load[1] + B_FAST[0] / 2, y_floor - 0.04), -0.14, f"loop {rl['Wl']*1e3:.0f}")
        s.label(ax, (xs_sup[1] + B_SUP[0] / 2, yb_s), (xs_sup[1] - 0.3, y_floor - 0.3), "support bundles (four, the center pair touching), 150 by 67 mm", COL["support"], ha="left")
        s.label(ax, (xs_load[0] + B_FAST[0] / 2, yb_l), (x_load - 0.2, y_floor - 0.42), "load bundles (two), 200 by 50 mm", COL["fast"], ha="left")
        s.label(ax, (xl1, y_floor - 0.02), (xl1 - 0.1, y_floor - 0.54), "casing floor (schematic; cryostats and mounts not drawn)", COL["grey"], ha="right")
        s.scale_bar(ax, xl0, y_floor - 0.62, 0.2, "200 mm")
    sh = Sheet("Drawing 1b: the right-hand underside, enlarged (one support track and one load track)", "01b_underside_right_detail",
               "To scale (200 mm bar): the plates, the gap, the bundles and the loop widths; the casing floor and the labels are schematic. Dimensions in millimeters; the loop widths are between leg centers.",
               legend=("support", "fast", "plate"),
               notes=[f"The support track: two {rs['Wl']:.1f} m loops in series as a figure-eight, their four legs the four bundles (the center pair, one leg of each loop, touching). "
                      f"The load track: one {rl['Wl']:.1f} m loop, its two legs the two bundles. Each bundle is 100 mm below its plate. The {D_SPACING*1e3:.0f} mm between the "
                      f"support plate's edge and the load plate's is round 5's; the plates and their spacings take {2*rs['Wp']+2*rl['Wp']+2*D_SPACING:.2f} m of the {SIDE:.2f} m underside and {SIDE-2*rs['Wp']-2*rl['Wp']-2*D_SPACING:.2f} m is left between the support plates."])
    sh.panel(draw)
    return sh


def sheet_1c(D):
    rs, rl = D["rs"], D["rl"]; half = SIDE / 2; tF = rl["T"]
    x_lift = D["x_lift"]; y_top = half; y_ceil = y_top + tF + GAP + B_FAST[1] + 0.06
    xl0 = x_lift - rl["Wp"] / 2 - 0.25; xl1 = half + 0.05

    def draw(ax):
        s = sh
        s.rect(ax, xl0, y_top - 0.1, xl1 - xl0, 0.1, COL["cable"], COL["cable_edge"])
        ax.text(x_lift, y_top - 0.05, "the cable's top face", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["cable_edge"])
        s.rect(ax, x_lift - rl["Wp"] / 2, y_top, rl["Wp"], tF, COL["plate_fill"], COL["plate"], hatch=HATCH["plate"])
        xs = [x_lift - rl["Wl"] / 2 - B_FAST[0] / 2, x_lift + rl["Wl"] / 2 - B_FAST[0] / 2]
        yb = y_top + tF + GAP
        for x in xs:
            s.rect(ax, x, yb, B_FAST[0], B_FAST[1], COL["fast"], COL["fast"], lw=0.8)
        s.rect(ax, xl0, y_ceil, xl1 - xl0, 0.04, COL["white"], COL["grey"], lw=1.0)
        s.dim(ax, (xs[0] + B_FAST[0] / 2, y_ceil + 0.04), (xs[1] + B_FAST[0] / 2, y_ceil + 0.04), 0.14, f"loop {rl['Wl']*1e3:.0f}")
        s.dim(ax, (x_lift - rl["Wp"] / 2, y_ceil + 0.04), (x_lift + rl["Wp"] / 2, y_ceil + 0.04), 0.4, f"lift plate {rl['Wp']*1e3:.0f}")
        s.dim(ax, (x_lift - rl["Wp"] / 2, y_top), (x_lift - rl["Wp"] / 2, y_top + tF), -0.08, f"{tF*1e3:.0f}")
        s.dim(ax, (xs[0] - 0.04, y_top + tF), (xs[0] - 0.04, yb), 0.0, f"gap {GAP*1e3:.0f}")
        s.dim(ax, (xs[1], yb), (xs[1] + B_FAST[0], yb), -0.04, f"{B_FAST[0]*1e3:.0f}")
        s.dim(ax, (xs[0] + B_FAST[0], yb), (xs[0] + B_FAST[0], yb + B_FAST[1]), 0.06, f"{B_FAST[1]*1e3:.0f}")
        s.label(ax, (xs[1] + B_FAST[0], yb + B_FAST[1] / 2), (half + 0.3, yb + B_FAST[1] / 2), "lift bundles (two), 200 by 50 mm, on the casing ceiling", COL["fast"], ha="left")
        s.label(ax, (xl1, y_ceil + 0.02), (half + 0.3, y_ceil + 0.02), "casing ceiling (schematic; cryostats and mounts not drawn)", COL["grey"], ha="left")
        s.label(ax, (half, y_top - 0.05), (half - 0.05, y_top - 0.25), "cable edge: the plate's outer edge is flush with it", COL["cable_edge"], ha="right")
        s.scale_bar(ax, xl0, y_top - 0.42, 0.2, "200 mm")
    sh = Sheet("Drawing 1c: one lift track, enlarged", "01c_lift_track_detail",
               "To scale (200 mm bar); the casing ceiling and the labels are schematic. Dimensions in millimeters; the loop width is between leg centers.",
               legend=("fast", "plate"),
               notes=[f"The lift track is the load track turned over: the same {rl['Wp']:.1f} m by {rl['T']*1e3:.0f} mm plate, on the top face at the cable's edge, the same {rl['Wl']:.1f} m loop in two 200 by 50 mm bundles, "
                      f"on the casing ceiling 100 mm above the plate. It pulls the casing down toward the cable; the load track below pulls it up. Each holds 196.5 kN/m (lift) or 198 kN/m (load) "
                      f"at full pull with the worst tape at half its limit, 1.46 times that at the rating."])
    sh.panel(draw)
    return sh


# ------------------------------------------------------------------------------------------------------------------------------
# 2  the fast bundle, stack by stack;  4  the support bundle
# ------------------------------------------------------------------------------------------------------------------------------
def draw_bundle(s, ax, lay, bundle, n_cx, n_cy, color, h_draw=None):
    """A bundle of n_cx by n_cy cells, every cell with its whole stacks (faces horizontal: w_tape wide, t_s thick), the walls and films; to scale, in mm."""
    wb, hb = bundle[0] * 1e3, (h_draw if h_draw else bundle[1]) * 1e3
    t_fam, t_turn, t_s, w_c, w_used, h_c = (lay[k] * 1e3 for k in ("t_fam", "t_turn", "t_s", "w_c", "w_used", "h_c"))
    n_x, n_y = lay["n_x"], lay["n_y"]
    w_tape = lay["w_tape"] * 1e3; N = lay["N"]
    s.rect(ax, 0, 0, wb, hb, COL["coolant"], COL["black"])
    for j in range(n_cy):
        for i in range(n_cx):
            x0 = t_fam + i * (w_c + t_fam); y0 = t_fam + j * (h_c + t_fam)
            s.rect(ax, x0, y0, w_c, h_c, COL["white"], COL["grey"], lw=0.6)
            k = 0
            for r in range(n_y):
                for c in range(n_x):
                    if k >= N:
                        break
                    s.rect(ax, x0 + c * (w_tape + t_turn), y0 + r * (t_s + t_turn), w_tape, t_s, color, color, lw=0.3, z=3)
                    k += 1
    s.dim(ax, (0, 0), (wb, 0), -4.0, f"{wb:.0f}")
    s.dim(ax, (wb, 0), (wb, hb), 3.0, f"{hb:.1f}" if abs(hb - round(hb)) > 0.05 else f"{hb:.0f}")
    i = n_cx - 1; j = n_cy - 1
    x0 = t_fam + i * (w_c + t_fam); y0 = t_fam + j * (h_c + t_fam)
    s.dim(ax, (x0, y0 + h_c), (x0 + w_c, y0 + h_c), 2.0 + (hb - y0 - h_c), f"cell {w_c:.1f}")
    s.dim(ax, (x0 + w_c, y0), (x0 + w_c, y0 + h_c), (wb - x0 - w_c) + 13, f"{h_c:.1f}")
    s.label(ax, (t_fam + w_tape / 2, t_fam + t_s / 2), (-4, -13), f"one turn: a stack of {lay['layers']:.0f} soldered layers of {w_tape:.0f} mm tape, {w_tape:.0f} by {t_s:.1f} mm; film between turns {t_turn:.2f} mm", color, ha="left")
    s.label(ax, (t_fam / 2, hb / 2), (-4, hb + 7), f"walls between cells and to ground {t_fam:.2f} mm", COL["black"], ha="left")


def sheet_2(D):
    rl = D["rl"]; f10, f14 = D["fast"][10], D["fast"][14]
    A = cell_layout(10, f10["layers"], f10["V"], f10["V_turn"], 5.0, 0.003, B_FAST[0], 10, 2)
    B = cell_layout(14, f14["layers"], f14["V"], f14["V_turn"], 10.0, 0.003, B_FAST[0], 10, 2)
    C = cell_layout(14, f14["layers"], f14["V"], f14["V_turn"], 5.0, 0.003, B_FAST[0], 10, 2)
    D["lay2"] = dict(A=A, B=B, C=C)

    def panel(lay):
        def draw(ax):
            draw_bundle(sh, ax, lay, B_FAST, 10, 2, COL["fast"], h_draw=max(B_FAST[1], lay["h_need"]))
            sh.scale_bar(ax, 0, -22, 50, "50 mm")
        return draw

    def cap(lay, name, N, E, f):
        fit = "fits" if lay["h_need"] <= B_FAST[1] + 1e-9 else f"does not fit: the bundle must be {lay['h_need']*1e3:.1f} mm tall, drawn so"
        return (f"Panel {name}: {N} series turns of {lay['layers']:.0f} layers, {f['I']:,.0f} A per turn, {f['V']/1e3:.1f} kV at the terminals at 10 ms, insulation at {E:.0f} kV/mm: "
                f"{lay['n_x']} stacks across and {lay['n_y']} up per cell, cell {lay['w_c']*1e3:.1f} by {lay['h_c']*1e3:.1f} mm, {lay['h_need']*1e3:.1f} mm of the 50 needed ({fit}).")
    sh = Sheet("Drawing 2: a fast track's bundle (lift or load) in cross-section, 200 by 50 mm, 20 cells (10 by 2), one cell per overlapping 25 m coil", "02_fast_bundle_cross_section",
               "To scale in every panel (50 mm bar); the dimensions in millimeters. The stacks are laid out one by one with their films and walls: the cell and bundle heights are what the stacks need, not an area fill.",
               legend=("fast",),
               notes=[f"Every cell is one 25 m coil's cross-section: its series turns, each turn a stack of soldered 3 mm tape layers laid face up (the tape's width across, the layers piled up; faces vertical "
                      f"does not fit a 10 by 2 cell grid in 200 mm at 10 turns of {f10['layers']:.0f} layers). The film between turns and the walls between cells and to ground are the terminal and turn voltages "
                      f"at 10 ms over the stated stress, at least 50 um. The space in a cell beside the stacks and under the 50 mm is the coolant's. Ampere-turns per coil {rl['NI_req']/N_FAM/1e3:.1f} kA, "
                      f"{rl['NI_req']/N_FAM/(0.5*rl['ic_req']):.0f} tapes per coil, stored energy {f10['E']/1e6:.3f} MJ per coil at full pull.",
                      "Panel A is the design case; B and C keep 14 turns in view (Paul has not ruled it out). The report's Section 1 gives the numbers the stack-by-stack layout moved."])
    sh.panel(panel(A), cap(A, "A", 10, 5.0, f10)); sh.panel(panel(B), cap(B, "B", 14, 10.0, f14)); sh.panel(panel(C), cap(C, "C", 14, 5.0, f14))
    return sh


def sheet_4(D):
    rs, S = D["rs"], D["sup"]
    n_turn = S["n1"]
    lay = cell_layout(n_turn, 10, 1000.0, 110.0, 5.0, 0.012, B_SUP[0], 5, 2)
    D["lay4"] = lay

    def draw(ax):
        draw_bundle(sh, ax, lay, B_SUP, 5, 2, COL["support"], h_draw=max(B_SUP[1], lay["h_need"]))
        sh.scale_bar(ax, 0, -22, 50, "50 mm")
    fit = "fits" if lay["h_need"] <= B_SUP[1] + 1e-9 else f"does not fit: the bundle must be {lay['h_need']*1e3:.1f} mm tall"
    sh = Sheet("Drawing 4: a support track's bundle in cross-section, 150 by 67 mm, 10 cells (5 by 2), one cell per overlapping 1 km coil", "04_support_bundle_cross_section",
               "To scale (50 mm bar); the dimensions in millimeters. The stacks laid out one by one with their films and walls.",
               legend=("support",),
               notes=[f"Each cell holds one coil's turns through this leg: {n_turn} stacks (the coil's {S['turns']:.0f} series turns go {S['n1']} times around one loop and {S['n2']} times around the other, "
                      f"Drawing 5), each a stack of 10 soldered layers of 12 mm tape, 12 by {lay['t_s']*1e3:.2f} mm, laid face up two across and {lay['n_y']} up: cell {lay['w_c']*1e3:.1f} by {lay['h_c']*1e3:.1f} mm, "
                      f"{lay['h_need']*1e3:.1f} mm of the 67 needed ({fit}). Walls {lay['t_fam']*1e3:.2f} mm for the 1 kV of the dump at 5 kV/mm; the film between turns {lay['t_turn']*1e3:.3f} mm "
                      f"(the 50 um minimum; 110 V across the whole coil while charging). {S['tapes_bundle']:,.0f} tapes through the bundle in the solve, {S['tapes_coil']:.0f} per coil (40.7 turns): "
                      f"the drawn 41 and 40 whole turns per leg carry 405 per bundle on average, 0.5 percent under the solve.",
                      f"The coil's current: {S['I_norm']:,.0f} A in normal operation (the support at 47.5 kN/m per track), {S['I_cap']:,.0f} A at the cap, the worst tape at half its limit."])
    sh.panel(draw)
    return sh


# ------------------------------------------------------------------------------------------------------------------------------
# 3  the fast coil in plan;  5  the support coil in plan
# ------------------------------------------------------------------------------------------------------------------------------
def sheet_3(D):
    rl, f10 = D["rl"], D["fast"][10]
    cx = 5.0; W = rl["Wl"]; L = L_FAST / cx

    def draw_loop(ax):
        s = sh
        for k in range(10):
            d = 0.012 * k
            s.rect(ax, d, -W / 2 + d, L - 2 * d, W - 2 * d, "none", COL["fast"], lw=1.1)
        ax.plot([0.0, -0.45], [-W / 2, -W / 2 - 0.45], color=COL["copper"], lw=2.4); ax.plot([0.12, -0.33], [-W / 2 + 0.11, -W / 2 - 0.45], color=COL["copper"], lw=2.4)
        s.rect(ax, -1.25, -W / 2 - 0.95, 1.0, 0.5, COL["white"], COL["copper"], lw=1.0)
        ax.text(-0.75, -W / 2 - 0.7, "warm: the switches\nand the store", fontsize=RULES["fs_label"], ha="center", va="center", color=COL["copper"])
        s.label(ax, (-0.2, -W / 2 - 0.2), (0.2, -W / 2 - 0.6), f"the coil's two ends and its pair of leads ({f10['I']:,.0f} A at full pull), at one corner", COL["copper"])
        s.dim(ax, (0, -W / 2 - 1.1), (L, -W / 2 - 1.1), -0.08, "25 m (drawn compressed 5 times)")
        s.dim(ax, (L + 0.15, -W / 2), (L + 0.15, W / 2), 0.1, f"{W:.1f} m")
        ax.text(L / 2, W / 2 + 0.15, f"one 25 m coil: a {W:.1f} m loop of 10 series turns, {f10['V']/1e3:.1f} kV at the terminals at 10 ms, {f10['E']/1e6:.3f} MJ stored",
                fontsize=RULES["fs_label"], ha="center", va="bottom", color=COL["fast"])
        s.scale_bar(ax, 0, -W / 2 - 1.55, 1.0, "5 m along the track (1 m drawn)")

    def draw_stagger(ax):
        s = sh
        st = L_FAST / N_FAM / cx
        for k in range(N_FAM):
            y = -k * 0.07
            ax.plot([k * st, k * st + L], [y, y], color=COL["fast"] if k % 2 == 0 else "#f2c08a", lw=3, solid_capstyle="butt")
            ax.plot([k * st], [y], marker="s", ms=4, color=COL["copper"])
        for k in range(N_FAM, N_FAM + 6):
            y = -(k - N_FAM) * 0.07
            ax.plot([k * st, k * st + L], [y, y], color="#cccccc", lw=3, solid_capstyle="butt")
            ax.plot([k * st], [y], marker="s", ms=4, color=COL["copper"])
        y_axis = -(N_FAM - 1) * 0.07 - 0.25
        x_end = (N_FAM + 6) * st + L
        ax.plot([0, x_end], [y_axis, y_axis], color=COL["black"], lw=1.0)
        for k in range(N_FAM + 6):
            ax.plot([k * st, k * st], [y_axis - 0.03, y_axis + 0.03], color=COL["copper"], lw=1.5)
        ax.text(x_end + 0.1, y_axis, "the track: a lead pair\nevery 1.25 m (copper ticks)", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["copper"])
        s.dim(ax, (0, y_axis - 0.12), (L, y_axis - 0.12), -0.1, "25 m")
        s.dim(ax, (0, y_axis - 0.5), (st, y_axis - 0.5), -0.1, "1.25 m")
        ax.plot([L * 0.6, L * 0.6], [y_axis, 0.05], color=COL["black"], lw=0.8, ls="--")
        ax.text(L * 0.6 - 0.05, 0.08, "any point of the track lies under 20 coils", fontsize=RULES["fs_label"], ha="right", va="bottom")
        ax.text(x_end + 0.1, 0.0, "the next coils (gray)", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["grey"])
    sh = Sheet("Drawing 3: a fast coil (lift or load) in plan, and the 20 coils' stagger along the track", "03_fast_coil_plan",
               "Not to scale along the track: the length is drawn compressed 5 times; the loop's width is to that same compressed scale (bar: 5 m of track = 1 m drawn; the 0.6 m width is drawn 0.12 m). The 10 turns are drawn as nested lines; the real turns are 3 mm stacks side by side in a 200 by 50 mm bundle (Drawing 2).",
               legend=("fast", "copper"),
               notes=["The 20 coils (families) of a track are each 25 m long and start 1.25 m apart, so every point of the track is under 20 coils and the coils' ends, with their lead pairs, "
                      "come every 1.25 m along the track (one coil's two ends are at the same place, its loop closing on itself). The squares mark each coil's end; the ticks on the track line are the lead pairs.",
                      f"Each coil's leads go out of the cold to its own switches and store (Drawing 6b); {f10['I']:,.0f} A per turn at full pull, {f10['V']/1e3:.1f} kV at the terminals with the 10 ms rise."])
    sh.panel(draw_loop, "The coil"); sh.panel(draw_stagger, "The stagger of a track's 20 coils (each line one coil; the rows only separate them)")
    return sh


def sheet_5(D):
    rs, S = D["rs"], D["sup"]
    cx = 200.0; WX = 5.0                                       # the length compressed 200 times; the widths drawn 5 times that compressed scale
    W = rs["Wl"] * WX; b = B_SUP[0] * WX; L = L_SUP / cx
    y_c1 = b / 2; y_o1 = y_c1 + W; y_c2 = -b / 2; y_o2 = y_c2 - W

    def draw_coil(ax):
        s = sh
        n_draw = 3; sp = 0.13
        for k in range(n_draw):
            d = sp * k
            s.rect(ax, d, y_c1 + d, L - 2 * d, W - 2 * d, "none", COL["support"], lw=1.4)
            s.rect(ax, d, y_o2 + d, L - 2 * d, W - 2 * d, "none", COL["support"], lw=1.4)
        for k in range(n_draw):
            d = sp * k
            ax.annotate("", xy=(L * 0.55, y_o1 - d), xytext=(L * 0.45, y_o1 - d), arrowprops=dict(arrowstyle="-|>", color=COL["support"], lw=1.2, mutation_scale=12))
            ax.annotate("", xy=(L * 0.45, y_c1 + d), xytext=(L * 0.55, y_c1 + d), arrowprops=dict(arrowstyle="-|>", color=COL["support"], lw=1.2, mutation_scale=12))
            ax.annotate("", xy=(L * 0.45, y_c2 - d), xytext=(L * 0.55, y_c2 - d), arrowprops=dict(arrowstyle="-|>", color=COL["support"], lw=1.2, mutation_scale=12))
            ax.annotate("", xy=(L * 0.55, y_o2 + d), xytext=(L * 0.45, y_o2 + d), arrowprops=dict(arrowstyle="-|>", color=COL["support"], lw=1.2, mutation_scale=12))
        xc = L - sp * (n_draw - 1)
        ax.plot([xc, xc], [y_c1 + sp * (n_draw - 1), y_c2 - sp * (n_draw - 1)], color=COL["support"], lw=1.4, ls=":")
        s.label(ax, (xc, 0), (L + 0.9, 0.0), f"the crossover (dotted): after {S['n1']} turns around loop 1\nthe winding goes on to loop 2 for {S['n2']}", COL["support"])
        ax.text(L * 0.5, y_o1 + 0.12, f"loop 1: {S['n1']} turns, clockwise as drawn (3 of them drawn)", fontsize=RULES["fs_label"], ha="center", va="bottom", color=COL["support"])
        ax.text(L * 0.5, y_o2 - 0.12, f"loop 2: {S['n2']} turns, counterclockwise (3 drawn): the center legs carry both loops' current the same way", fontsize=RULES["fs_label"], ha="center", va="top", color=COL["support"])
        s.dim(ax, (L + 0.25, y_c1), (L + 0.25, y_o1), 0.1, f"{rs['Wl']:.1f} m")
        s.dim(ax, (L + 0.25, y_o2), (L + 0.25, y_c2), 0.1, f"{rs['Wl']:.1f} m")
        s.dim(ax, (-0.25, y_c2), (-0.25, y_c1), -0.1, "150 mm, the center pair")
        s.dim(ax, (0, y_o2 - 0.75), (L, y_o2 - 0.75), -0.1, "1,000 m (drawn compressed 200 times)")
        xh = -2.4
        items = (("free-wheel diode D2", COL["black"]), ("strip and its heater", COL["copper"]), ("dump switch and resistor", COL["black"]), ("lead contacts (to the warm supply)", COL["copper"]))
        for i, (lab, c) in enumerate(items):
            y = y_o1 - 0.2 - i * 0.5
            s.rect(ax, xh, y - 0.15, 0.3, 0.3, COL["white"], c, lw=1.0)
            ax.text(xh - 0.1, y, lab, fontsize=RULES["fs_label"], ha="right", va="center", color=c)
        ax.plot([0, -0.6, -0.6, xh + 0.3], [y_c1, y_c1, y_o1 - 0.2 - 0.5, y_o1 - 0.2 - 0.5], color=COL["copper"], lw=1.6)
        ax.plot([0, -0.45, -0.45, xh + 0.3], [y_o1, y_o1, y_o1 - 0.2 - 1.5, y_o1 - 0.2 - 1.5], color=COL["copper"], lw=1.6)
        ax.text(xh - 0.1, y_o1 + 0.35, "the coil's end, in the cold:", fontsize=RULES["fs_label"], ha="right", va="center", color=COL["grey"])
        s.scale_bar(ax, 0, y_o2 - 1.3, 0.5, "100 m along the track (0.5 m drawn)")

    def draw_stagger(ax):
        s = sh
        st = L_SUP / N_SUP / cx
        for k in range(N_SUP):
            y = -k * 0.12
            ax.plot([k * st, k * st + L], [y, y], color=COL["support"] if k % 2 == 0 else "#7f9fc4", lw=3.5, solid_capstyle="butt")
            ax.plot([k * st, k * st], [y - 0.05, y + 0.05], color=COL["black"], lw=1.2)
        for k in range(N_SUP, N_SUP + 4):
            y = -(k - N_SUP) * 0.12
            ax.plot([k * st, k * st + L], [y, y], color="#cccccc", lw=3.5, solid_capstyle="butt")
        y_d = -(N_SUP - 1) * 0.12 - 0.3
        s.dim(ax, (0, y_d), (L, y_d), -0.12, "1 km")
        s.dim(ax, (0, y_d - 0.45), (st, y_d - 0.45), -0.12, "100 m")
        x_end = (N_SUP + 3) * st + L
        ax.text(x_end + 0.1, -0.2, "the next coils (gray)", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["grey"])
        ax.text(0, 0.15, "a bundle's ten coils, each 1 km, starting 100 m apart: every point of the track is inside 10 coils;\na coil's end (tick) lies where only its 9 neighbors run on", fontsize=RULES["fs_label"], ha="left", va="bottom", color=COL["support"])
    sh = Sheet("Drawing 5: a support coil in plan (a figure-eight of two loops in series), and the ten coils' stagger", "05_support_coil_plan",
               "Not to scale: along the track the length is drawn compressed 200 times (bar: 100 m of track = 0.5 m drawn); across the track the loop widths and the bundle pair are drawn 5 times that compressed scale, to be seen. Six of the turns are drawn as nested lines.",
               legend=("support", "copper"),
               notes=[f"One coil is {S['turns']:.0f} series turns of 10 soldered layers of 12 mm tape: {S['n1']} turns around loop 1, then the crossover, then {S['n2']} around loop 2 in the opposite sense, so that loop 1's "
                      f"inner leg and loop 2's inner leg (the center pair of bundles) carry the current the same way and the two outer legs the other way. Through loop 1's two legs run {S['n1']} of this coil's "
                      f"turns, through loop 2's {S['n2']}; the next coil of the bundle puts its {S['n1']} on loop 2, so that each of the four bundles carries 405 turns of the ten coils, the solve's 40.7 per "
                      f"coil (4,071 per bundle) within 0.5 percent.",
                      f"Current {S['I_norm']:,.0f} A in normal operation, {S['I_cap']:,.0f} A at the cap; about {S['L']:.0f} H with the ten coils of a bundle charged together (Stage B's 21.3 H at 55 turns scaled as the turns squared). "
                      "The coil's end carries its free-wheel diode, its persistent strip with heater, its dump branch and its lead contacts (Drawing 6a)."])
    sh.panel(draw_coil, "The coil"); sh.panel(draw_stagger, "The stagger of a bundle's ten coils (each line one coil)")
    return sh


# ------------------------------------------------------------------------------------------------------------------------------
# 6a, 6b  the circuits, standard symbols
# ------------------------------------------------------------------------------------------------------------------------------
def wire(ax, pts, color=COL["black"], lw=1.4):
    xs, ys = zip(*pts); ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round")


def node(ax, x, y, color=COL["black"]):
    ax.plot([x], [y], marker="o", ms=4, color=color)


def resistor(ax, x, y, vertical=True, size=0.6, color=COL["black"]):
    if vertical:
        ax.add_patch(Rectangle((x - 0.1, y - size / 2), 0.2, size, facecolor="white", edgecolor=color, lw=1.4, zorder=3))
    else:
        ax.add_patch(Rectangle((x - size / 2, y - 0.1), size, 0.2, facecolor="white", edgecolor=color, lw=1.4, zorder=3))


def diode(ax, x, y, up=True, color=COL["black"], size=0.26):
    h = size
    tri = [(x - h * 0.75, y - h / 2), (x + h * 0.75, y - h / 2), (x, y + h / 2)] if up else [(x - h * 0.75, y + h / 2), (x + h * 0.75, y + h / 2), (x, y - h / 2)]
    ax.add_patch(Polygon(tri, closed=True, facecolor=color, edgecolor=color, zorder=3))
    yb = y + h / 2 if up else y - h / 2
    ax.plot([x - h * 0.75, x + h * 0.75], [yb, yb], color=color, lw=1.8, zorder=3)


def capacitor(ax, x, y, color=COL["black"], size=0.5):
    ax.plot([x - size / 2, x + size / 2], [y + 0.07, y + 0.07], color=color, lw=2.2); ax.plot([x - size / 2, x + size / 2], [y - 0.07, y - 0.07], color=color, lw=2.2)


def switch(ax, x, y, vertical=True, color=COL["black"], size=0.5):
    if vertical:
        ax.plot([x], [y - size / 2], marker="o", ms=4, color=color); ax.plot([x], [y + size / 2], marker="o", ms=4, color=color)
        ax.plot([x, x + 0.22], [y - size / 2, y + size / 2 - 0.05], color=color, lw=1.6)
    else:
        ax.plot([x - size / 2], [y], marker="o", ms=4, color=color); ax.plot([x + size / 2], [y], marker="o", ms=4, color=color)
        ax.plot([x - size / 2, x + size / 2 - 0.05], [y, y + 0.22], color=color, lw=1.6)


def igbt(ax, x, y, color=COL["black"], size=0.6):
    h = size
    ax.plot([x - 0.12, x - 0.12], [y - h * 0.35, y + h * 0.35], color=color, lw=2.4)
    ax.plot([x - 0.26, x - 0.26], [y - h * 0.3, y + h * 0.3], color=color, lw=1.6)
    ax.plot([x - 0.5, x - 0.26], [y, y], color=color, lw=1.4)
    ax.plot([x - 0.12, x], [y + h * 0.25, y + h * 0.5], color=color, lw=1.6); ax.plot([x, x], [y + h * 0.5, y + h / 2], color=color, lw=1.4)
    ax.plot([x - 0.12, x], [y - h * 0.25, y - h * 0.5], color=color, lw=1.6); ax.plot([x, x], [y - h * 0.5, y - h / 2], color=color, lw=1.4)
    ax.annotate("", xy=(x - 0.01, y - h * 0.48), xytext=(x - 0.11, y - h * 0.27), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2, mutation_scale=9))


def hop(ax, x, y, color=COL["copper"], r=0.12, lw=2.4):
    """A wire crossing another without a connection: a half circle over it."""
    ax.add_patch(Arc((x, y), 2 * r, 2 * r, theta1=90, theta2=270, edgecolor=color, lw=lw, zorder=4))


def sheet_6a(D):
    S = D["sup"]

    def draw(ax):
        s = sh
        yw = 6.2
        s.rect(ax, 0.2, yw, 4.5, 1.3, COL["white"], COL["black"], lw=1.4)
        ax.text(2.45, yw + 0.65, "DC supply from the 100 kV bus\ncurrent-regulated, about 110 V\nfor 1 percent per second", fontsize=RULES["fs_label"], ha="center", va="center")
        wire(ax, [(4.7, yw + 0.95), (5.4, yw + 0.95)]); wire(ax, [(4.7, yw + 0.35), (5.4, yw + 0.35)])
        switch(ax, 5.65, yw + 0.95, vertical=False); switch(ax, 5.65, yw + 0.35, vertical=False)
        ax.text(5.65, yw + 1.5, "lead contacts, open between charges", fontsize=RULES["fs_label"], ha="center", va="bottom", color=COL["copper"])
        y1, y0 = 4.9, 1.0
        wire(ax, [(5.9, yw + 0.95), (6.3, yw + 0.95), (6.3, y1)], COL["copper"], 2.4)
        wire(ax, [(5.9, yw + 0.35), (7.1, yw + 0.35), (7.1, y1 + 0.12)], COL["copper"], 2.4); hop(ax, 7.1, y1); wire(ax, [(7.1, y1 - 0.12), (7.1, y0)], COL["copper"], 2.4)
        ax.axhline(5.6, color=COL["grey"], ls="--", lw=1.0); ax.text(14.1, 5.68, "warm", fontsize=RULES["fs_label"], color=COL["grey"], ha="right", va="bottom"); ax.text(14.1, 5.5, "cold, 70 K", fontsize=RULES["fs_label"], color=COL["grey"], ha="right", va="top")
        ax.text(7.4, 5.3, "a pair of copper leads", fontsize=RULES["fs_label"], color=COL["copper"], ha="left", va="center")
        wire(ax, [(1.0, y1), (12.4, y1)]); wire(ax, [(1.0, y0), (12.4, y0)])
        for x in (1.0, 3.5, 8.5, 11.4):
            node(ax, x, y1); node(ax, x, y0)
        node(ax, 6.3, y1, COL["copper"]); node(ax, 7.1, y0, COL["copper"])
        wire(ax, [(1.0, y0), (1.0, y1)]); switch(ax, 1.0, 3.9, vertical=True); resistor(ax, 1.0, 2.1, vertical=True, size=0.8)
        ax.text(0.7, 3.9, "dump switch", fontsize=RULES["fs_label"], ha="right", va="center"); ax.text(0.7, 1.9, f"dump resistor\n{1000/S['I_norm']:.1f} ohm: 1 kV\nat {S['I_norm']:,.0f} A", fontsize=RULES["fs_label"], ha="right", va="center")
        wire(ax, [(3.5, y0), (3.5, y1)]); s.rect(ax, 3.3, 2.4, 0.4, 1.4, "#f6d5c3", COL["copper"], lw=1.4)
        resistor(ax, 4.3, 3.1, vertical=True, size=0.8, color=COL["copper"]); wire(ax, [(4.3, 3.5), (4.3, 3.9), (4.8, 3.9)], COL["copper"], 1.2); wire(ax, [(4.3, 2.7), (4.3, 2.3), (4.8, 2.3)], COL["copper"], 1.2)
        ax.text(4.9, 3.1, "heater and\nits supply", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["copper"])
        ax.text(3.5, 2.1, "the strip: bare REBCO tape,\nheld normal while charging", fontsize=RULES["fs_label"], ha="center", va="top", color=COL["copper"])
        wire(ax, [(8.5, y0), (8.5, y1)]); diode(ax, 8.5, 3.0, up=True)
        ax.text(8.85, 3.0, "D2, free-wheel\ndiode, 1.0 V at 70 K", fontsize=RULES["fs_label"], ha="left", va="center")
        wire(ax, [(11.4, y0), (11.4, 2.0)]); wire(ax, [(11.4, 4.0), (11.4, y1)])
        t = np.linspace(0, math.pi, 20)
        for k in range(5):
            ax.plot(11.4 + 0.2 * np.sin(t), 2.2 + k * 0.4 - 0.2 * np.cos(t), color=COL["support"], lw=2.0)
        ax.text(11.8, 3.0, f"the coil:\n{S['turns']:.0f} turns of 10 layers\nabout {S['L']:.0f} H (ten coils\ntogether)\n{S['I_norm']:,.0f} A normal\n{S['I_cap']:,.0f} A at the cap", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["support"])
    sh = Sheet("Drawing 6a: a support coil's circuit (one 1 km coil; each of a bundle's ten has its own)", "06a_support_coil_circuit",
               "Not to scale: a circuit diagram with standard symbols. Warm above the dashed line, 70 K below.",
               legend=("support", "copper"),
               notes=["Normal hold: the lead contacts open, the strip superconducting (heater off), D2 carrying nothing, the dump switch open; the current circulates coil to strip with no loss. "
                      "Charge: heater on (the strip resistive, 21 ohm at Stage B's rule), the contacts closed, the supply drives the coil at about 110 V; D2 blocks. "
                      "Close (round 5, Section 8): order (a) removes the supply, the coil free-wheels through D2 while the strip cools and takes the current; order (b) holds the current until the strip is "
                      "closed, then ramps the supply down, in current regulation throughout. The dump branch is closed only to discharge a quenched coil: 1 kV across the coil at its current.",
                      "The second lead crosses the upper rail without touching it (the half circle) and lands on the lower rail. Dump switch and lead contacts: mechanical or semiconductor, not chosen here."])
    sh.panel(draw)
    return sh


def sheet_6b(D):
    f10, f14 = D["fast"][10], D["fast"][14]

    def draw(ax):
        s = sh
        s.rect(ax, 0.2, 5.0, 2.4, 1.4, COL["white"], COL["black"], lw=1.4)
        ax.text(1.4, 5.7, "charger from the\n100 kV bus: refills\nthe store's losses", fontsize=RULES["fs_label"], ha="center", va="center")
        xp = 3.8; yp, ym = 6.2, 4.2
        wire(ax, [(2.6, 6.0), (xp, 6.0), (xp, yp)]); wire(ax, [(2.6, 4.4), (xp, 4.4), (xp, ym)])
        capacitor(ax, xp, 5.2); ax.text(xp + 0.35, 5.2, f"store C: {f10['E']/1e6:.3f} MJ at full pull\n(film capacitors, warm)", fontsize=RULES["fs_label"], ha="left", va="center", color="#7a4a00")
        xa, xb = 7.6, 9.9; yA = 5.2
        wire(ax, [(xp, yp), (xb + 0.6, yp)]); wire(ax, [(xp, ym), (xb + 0.6, ym)]); node(ax, xp, yp); node(ax, xp, ym)
        ax.text(xb + 0.75, yp, "+", fontsize=RULES["fs_label"] + 2, ha="left", va="center"); ax.text(xb + 0.75, ym, "−", fontsize=RULES["fs_label"] + 2, ha="left", va="center")
        wire(ax, [(xa, yp), (xa, yA)]); igbt(ax, xa, 5.7); ax.text(xa - 0.6, 5.7, "S1", fontsize=RULES["fs_label"], ha="right", va="center")
        wire(ax, [(xb, yA), (xb, ym)]); igbt(ax, xb, 4.7); ax.text(xb - 0.6, 4.7, "S2", fontsize=RULES["fs_label"], ha="right", va="center")
        wire(ax, [(xb, yA), (xb, yp)]); diode(ax, xb, 5.7, up=True); ax.text(xb + 0.3, 5.7, "D1", fontsize=RULES["fs_label"], ha="left", va="center")
        wire(ax, [(xa, ym), (xa, yA)]); diode(ax, xa, 4.7, up=True); ax.text(xa + 0.3, 4.7, "D2", fontsize=RULES["fs_label"], ha="left", va="center")
        node(ax, xa, yA); node(ax, xb, yA)
        ax.axhline(3.3, color=COL["grey"], ls="--", lw=1.0); ax.text(0.2, 3.38, "warm", fontsize=RULES["fs_label"], color=COL["grey"], ha="left", va="bottom"); ax.text(0.2, 3.2, "cold, 70 K", fontsize=RULES["fs_label"], color=COL["grey"], ha="left", va="top")
        wire(ax, [(xa, yA), (xa, 2.2)], COL["copper"], 2.4); wire(ax, [(xb, yA), (xb, 2.2)], COL["copper"], 2.4)
        ax.text(xb + 0.3, 2.75, f"pulse-sized leads: {f10['I']:,.0f} A at full pull\n(10 turns), {f14['I']:,.0f} A (14)", fontsize=RULES["fs_label"], ha="left", va="center", color=COL["copper"])
        t = np.linspace(0, math.pi, 20)
        for k in range(5):
            ax.plot(xa + 0.23 + k * 0.46 - 0.23 * np.cos(t), 2.2 + 0.23 * np.sin(t), color=COL["fast"], lw=2.0)
        ax.text((xa + xb) / 2, 1.8, f"the coil: 10 (or 14) series turns,\n{f10['V']/1e3:.1f} ({f14['V']/1e3:.1f}) kV at the terminals at 10 ms", fontsize=RULES["fs_label"], ha="center", va="top", color=COL["fast"])
    sh = Sheet("Drawing 6b: a fast coil's circuit (one 25 m lift or load coil; 20 per 25 m of track): the asymmetric two-quadrant bridge and the capacitor store", "06b_fast_coil_circuit",
               "Not to scale: a circuit diagram with standard symbols (IGBT, diode, capacitor). Warm above the dashed line, 70 K below.",
               legend=("fast", "copper"),
               notes=["Charge: S1 and S2 on, +V across the coil, C discharges into it in 10 ms. Hold: one switch on, the other off, the current free-wheels through that switch and a diode at zero voltage. "
                      "Return: both switches off, the current flows through D2 and D1 back into C, −V across the coil, 10 ms. One current direction in the coil, both voltage signs across it (round 5, Section 4).",
                      f"Each switch position is 5 StakPak press-pack IGBTs in series at 10 turns ({f10['V']/1e3:.1f} kV), 7 at 14 ({f14['V']/1e3:.1f} kV); the diodes the same count. "
                      "The store is a film capacitor bank at every coil (round 5's way (a)); it also carries the gap control's small corrections and the ride-through when the bus is lost. "
                      "No persistent switch and no free-wheel path outside the bridge: the coil holds nothing between pulses."])
    sh.panel(draw)
    return sh


# ------------------------------------------------------------------------------------------------------------------------------
def main():
    global OUT
    only = None
    for a in sys.argv[1:]:
        if a.startswith("--out="): OUT = a[6:]
        elif a.startswith("--only="): only = set(a[7:].split(","))
        elif a in ("-h", "--help"): print(__doc__); return
    D = design()
    print(f"Drawings into {OUT}: support {D['rs']['Wp']:.2f} m x {D['rs']['T']*1e3:.0f} mm with {D['rs']['Wl']:.1f} m loops, fast {D['rl']['Wp']:.2f} m x {D['rl']['T']*1e3:.0f} mm with a {D['rl']['Wl']:.1f} m loop ({RULES['status']})")
    sheets = {"1a": sheet_1a, "1b": sheet_1b, "1c": sheet_1c, "2": sheet_2, "3": sheet_3, "4": sheet_4, "5": sheet_5, "6a": sheet_6a, "6b": sheet_6b}
    for key, fn in sheets.items():
        if only and key not in only:
            continue
        fn(D).render(OUT)
    f10, f14 = D["fast"][10], D["fast"][14]
    print(f"  fast coil: 10 turns x {f10['layers']:.0f} layers, {f10['I']:,.0f} A, {f10['V']/1e3:.1f} kV; 14 turns x {f14['layers']:.0f} layers, {f14['I']:,.0f} A, {f14['V']/1e3:.1f} kV; {f10['E']/1e6:.3f} MJ")
    if "lay2" in D:
        for name, lay in D["lay2"].items():
            print(f"  bundle 2{name}: walls {lay['t_fam']*1e3:.2f} mm, film {lay['t_turn']*1e3:.3f} mm, stack {lay['w_tape']*1e3:.0f} x {lay['t_s']*1e3:.2f} mm, {lay['n_x']} x {lay['n_y']} per cell, cell {lay['w_c']*1e3:.2f} x {lay['h_c']*1e3:.2f} mm, height needed {lay['h_need']*1e3:.2f} of 50 mm")
    if "lay4" in D:
        lay = D["lay4"]
        print(f"  support bundle: walls {lay['t_fam']*1e3:.2f} mm, film {lay['t_turn']*1e3:.3f} mm, stack 12 x {lay['t_s']*1e3:.2f} mm, {lay['n_x']} x {lay['n_y']} per cell ({lay['N']} stacks), cell {lay['w_c']*1e3:.2f} x {lay['h_c']*1e3:.2f} mm, height needed {lay['h_need']*1e3:.2f} of 67 mm")
    S = D["sup"]
    print(f"  support coil: {S['turns']:.1f} series turns ({S['n1']} + {S['n2']}), {S['I_norm']:,.0f} A normal, {S['I_cap']:,.0f} A cap, about {S['L']:.0f} H")


if __name__ == "__main__":
    main()
