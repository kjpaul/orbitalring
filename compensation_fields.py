#!/usr/bin/env python3
"""
Solves the field cases of the compensation and support tracks (linear and saturating iron) and writes the cache data/compensation_fields.json that the bearing scripts read.
Book: Orbital Ring Engineering, Volume I: Chapter 7 (library, and the air-core case that a flat array gives 0.2 kN/m).
Run: python compensation_fields.py (recomputes the cache, about 4 minutes)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

Part B of the 2026-10-11 task: can the compensation coils make the field?

Cross-section of one 1 m track facing its 1.4 m plate across a 100 mm gap,
solved with magnetostatic_2d.py.  Cases:

  A1  air core, 20 flat racetrack coils 5 cm wide side by side, all one polarity
      (the inner legs cancel: it is one 1 m loop)
  A2  air core, 20 coils on edge (legs stacked 100 mm apart): a 1 m current sheet
      with its return 100 mm behind
  A3  air core, 20 flat coils with alternating polarity (pitch 50 mm)
  A4  air core, 1 m current sheet with its return far away (the Chapter 7 picture)
  I1  iron: laminated yoke with alternating poles, pitch 250 mm (4 poles of 200 mm)
  I2  iron: pitch 500 mm (2 poles of 400 mm)
  I3  iron: pitch 125 mm (8 poles of 100 mm)
  I4  iron: one pole pair at the track edges, 200 mm poles (C-core)

For each: field at the plate (mean of |B_y| over the pole faces and over the
whole track width), the pull per meter (Maxwell stress across the gap), the
field inside the plate (saturation check against 1.2 T), the field at the
winding (components across the track, B_x, and normal to the plate, B_y), all
per 100 kA-turns, then scaled to 0.74 T (both tracks of a pair sharing the
437 kN/m pull at 30 km/s) and 1.2 T.  The same for the casing support track
(40 coils of 25 mm, 52 turns at 400 A).

The field is linear in the current, so one solve per geometry gives every
current.  Iron: relative permeability 1,000 for the laminated yoke and 200 for
the sintered plate (levitation_config.py), linear; saturation is checked, not
modeled.  Cell size 5 mm.  [UNVERIFIED] the yoke and pole dimensions are a
design sketch, not a Chapter 7 value.

Reference: "Orbital Ring, Mass Driver and Space Elevator Engineering" by Paul G de Jong
"""

import json
import math
import os
import sys
import time

import numpy as np

import levitation_config as lc
import magnetostatic_2d as ms

MU0 = ms.MU0
GAP = 0.100
PLATE_T = lc.FERRO_PLATE_THICKNESS      # 0.1 m
PLATE_W = 1.4
PLATE_MU = lc.FERRO_MU_R                # 200
YOKE_MU = 1000.0
YOKE_B_SAT = 1.8                        # T, laminated steel yoke and poles [UNVERIFIED]
PLATE_B_SAT = lc.FERRO_B_SAT            # 1.2 T
TRACK_W = 1.0
NI_REF = 1e5                            # ampere-turns per track used in the solve (results scale linearly)
H = 0.005


def base_grid():
    g = ms.Grid(3.0, 1.6, H)                                   # x from -1.5 to 1.5, y from -0.8 to 0.8
    g.add_iron(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T, PLATE_MU, B_sat=PLATE_B_SAT)   # plate above the gap, coil face at y = 0
    return g


def evaluate(g, pole_faces, winding_boxes, label, cuts=None):
    """Fields from a solved grid.  pole_faces: list of (x_lo, x_hi) at y just below the plate; winding_boxes: list of
    (x_lo, x_hi, y_lo, y_hi); cuts: x positions where the flux along the plate is read (default: midway between poles)."""
    y_gap = GAP - 2 * H
    xs, By, Bx = g.line_By(y_gap, -TRACK_W / 2, TRACK_W / 2)
    B_track_mean = float(np.abs(By).mean())
    pole_sel = np.zeros_like(xs, dtype=bool)
    for (a, b) in pole_faces:
        pole_sel |= (xs >= a) & (xs <= b)
    B_pole_mean = float(np.abs(By[pole_sel]).mean()) if pole_sel.any() else B_track_mean
    F = g.force_on_box(-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, y_gap, GAP + PLATE_T + 2 * H)       # closed surface around the plate
    # flux carried along the plate: mean |B_x| through the plate thickness at the cuts between pole faces, the largest
    # of them; this is what saturates the plate
    if cuts is None:
        cuts = [0.5 * (pole_faces[k][1] + pole_faces[k + 1][0]) for k in range(len(pole_faces) - 1)] or [0.0]
    j_lo = int(round((GAP - g.y0) / g.h)) + 1; j_hi = int(round((GAP + PLATE_T - g.y0) / g.h)) - 1
    B_plate_int = 0.0
    for xc in cuts:
        i = int(round((xc - g.x0) / g.h))
        B_plate_int = max(B_plate_int, float(np.abs(g.Bx[i, j_lo:j_hi]).mean()))
    # stored magnetic energy per meter along the track, J/m: 1/2 integral of A J (all of it, gap, slots, leakage)
    W_mag = float(0.5 * np.sum(g.A * g.J) * g.h * g.h)
    wf = [g.region_field(*wb) for wb in winding_boxes]
    Bx_w = max(w["Bx_max"] for w in wf); By_w = max(w["By_max"] for w in wf)
    Bx_wm = float(np.mean([w["Bx_mean"] for w in wf])); By_wm = float(np.mean([w["By_mean"] for w in wf]))
    return dict(label=label, B_track=B_track_mean, B_pole=B_pole_mean, F=F, B_plate=B_plate_int, W_mag=W_mag, NI=NI_REF, n_groups=1, NI_total=NI_REF,
                Bx_w=Bx_w, By_w=By_w, Bx_wm=Bx_wm, By_wm=By_wm, iters=g.iters, t=g.solve_time, pole_frac=float(pole_sel.mean()))


def case_air_flat_same(NI=NI_REF):
    g = base_grid()
    n = 20; wc = TRACK_W / n
    for k in range(n):
        x0 = -TRACK_W / 2 + k * wc
        g.add_current(x0, x0 + 0.010, -0.03, 0.0, NI / n)          # go leg, 10 mm wide conductor at the coil's edge
        g.add_current(x0 + wc - 0.010, x0 + wc, -0.03, 0.0, -NI / n)
    g.solve()
    return evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.03, 0.0)], "A1 air, 20 flat coils, one polarity")


def case_air_edge(NI=NI_REF):
    g = base_grid()
    n = 20; wc = TRACK_W / n
    for k in range(n):
        x0 = -TRACK_W / 2 + k * wc
        g.add_current(x0, x0 + wc, -0.02, 0.0, NI / n)             # leg nearest the plate
        g.add_current(x0, x0 + wc, -0.12, -0.10, -NI / n)          # return leg 100 mm behind
    g.solve()
    return evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.02, 0.0), (-TRACK_W / 2, TRACK_W / 2, -0.12, -0.10)],
                    "A2 air, 20 coils on edge, return 100 mm behind")


def case_air_flat_alt(NI=NI_REF):
    g = base_grid()
    n = 20; wc = TRACK_W / n
    for k in range(n):
        x0 = -TRACK_W / 2 + k * wc
        s = 1 if k % 2 == 0 else -1
        g.add_current(x0, x0 + 0.010, -0.03, 0.0, s * NI / n)
        g.add_current(x0 + wc - 0.010, x0 + wc, -0.03, 0.0, -s * NI / n)
    g.solve()
    return evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.03, 0.0)], "A3 air, 20 flat coils, alternating")


def case_air_sheet_far(NI=NI_REF):
    g = ms.Grid(4.0, 3.0, 0.01)
    g.add_iron(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T, PLATE_MU)
    g.add_current(-TRACK_W / 2, TRACK_W / 2, -0.02, 0.0, NI)
    g.add_current(-TRACK_W / 2, TRACK_W / 2, -1.22, -1.20, -NI)   # return 1.2 m behind
    g.solve()
    return evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.02, 0.0)], "A4 air, 1 m sheet, return 1.2 m behind")


def case_iron_poles(pitch, pole_w, label, NI=NI_REF, slot_depth=0.30, yoke_t=0.20, plate_t=None, end_half=False, nonlinear=False):
    """Alternating poles across the track with the windings in the slots between them and a yoke below.  NI is the
    ampere-turns of the coil wound on each pole, so a slot between two poles carries 2 NI (same direction, since
    neighboring poles alternate) and the slot outside an end pole carries NI.
    end_half=False: n = W / pitch full poles of width pole_w, the first centered at pitch / 2 from the track edge.
    end_half=True: half-width poles at both track edges and n - 1 full poles between them, every slot the same width;
    the end poles then send all their flux inward and every cut of the plate carries half a pole's flux."""
    global PLATE_T
    plate_t_saved = PLATE_T
    if plate_t is not None:
        PLATE_T = plate_t
    g = base_grid()
    n_pole = int(round(TRACK_W / pitch))
    slot_w = pitch - pole_w
    y_top = 0.0
    faces = []
    if end_half:
        # half, slot, full, slot, ..., full, slot, half across W; outer slots (NI) just outside the track edges
        x = -TRACK_W / 2
        faces.append((x, x + pole_w / 2)); x += pole_w / 2 + slot_w
        for k in range(n_pole - 1):
            faces.append((x, x + pole_w)); x += pole_w + slot_w
        faces.append((x, x + pole_w / 2))
        slot_edges = [(-TRACK_W / 2 - slot_w, -TRACK_W / 2)] + [(faces[k][1], faces[k + 1][0]) for k in range(len(faces) - 1)] + [(TRACK_W / 2, TRACK_W / 2 + slot_w)]
        n_groups = len(faces)                      # every pole, half or full, has its own coil of NI
    else:
        centers = [-TRACK_W / 2 + (k + 0.5) * pitch for k in range(n_pole)]
        faces = [(xc - pole_w / 2, xc + pole_w / 2) for xc in centers]
        slot_edges = [(-TRACK_W / 2 - 0.05, faces[0][0])] + [(faces[k][1], faces[k + 1][0]) for k in range(n_pole - 1)] + [(faces[-1][1], TRACK_W / 2 + 0.05)]
        n_groups = n_pole
    x_min = min(slot_edges[0][0], faces[0][0]) - 0.05; x_max = max(slot_edges[-1][1], faces[-1][1]) + 0.05
    g.add_iron(x_min, x_max, -slot_depth - yoke_t, -slot_depth, YOKE_MU, B_sat=YOKE_B_SAT)   # yoke
    for (a, b) in faces:
        g.add_iron(a, b, -slot_depth, y_top, YOKE_MU, B_sat=YOKE_B_SAT)
    slots = []
    n_slots = len(slot_edges)
    for k, (x_lo, x_hi) in enumerate(slot_edges):
        if x_hi - x_lo < H:
            continue
        s = 1 if k % 2 == 0 else -1
        amp = NI if k in (0, n_slots - 1) else 2 * NI
        g.add_current(x_lo, x_hi, -slot_depth, y_top, s * amp)
        slots.append((x_lo, x_hi, -slot_depth, y_top))
    meta = dict(n_pole=n_pole, n_groups=n_groups, pole_w=pole_w, pitch=pitch, slot_w=slot_w, plate_t=PLATE_T, end_half=end_half)
    if nonlinear == "grid":
        PLATE_T = plate_t_saved
        return g, faces, slots, meta
    if nonlinear:
        g.solve_nonlinear(force_box=(-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, GAP - 2 * H, GAP + PLATE_T + 2 * H))
    else:
        g.solve()
    r = finish_case(g, faces, slots, label, NI, meta, bool(nonlinear))
    PLATE_T = plate_t_saved
    return r


def finish_case(g, faces, slots, label, NI, meta, nonlinear):
    r = evaluate(g, faces, slots, label)
    r["NI"] = NI; r["nonlinear"] = nonlinear; r["h"] = g.h
    if nonlinear:
        r["outer"] = g.outer_iters; r["outer_change"] = g.outer_change; r["outer_dF"] = g.outer_dF
        m_pl = g.box(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + meta["plate_t"])
        r["B_plate_max"] = float(np.sqrt(g.Bx[m_pl]**2 + g.By[m_pl]**2).max())
        r["mu_plate_min"] = float(g.mu_r[m_pl].min())
        # how much of the pole face and of the plate is saturated (mu_r under a tenth of its start value)
        r["plate_sat_frac"] = float(np.mean(g.mu_r[m_pl] < 0.1 * PLATE_MU))
    r.update(meta)
    r["NI_total"] = meta["n_groups"] * NI
    return r


def case_iron_ccore(NI=NI_REF, pole_w=0.20, slot_depth=0.30, yoke_t=0.20):
    """One pole pair at the track edges with the winding in the window between them (one coil around the yoke,
    return below the yoke outside the circuit)."""
    g = base_grid()
    g.add_iron(-TRACK_W / 2, TRACK_W / 2, -slot_depth - yoke_t, -slot_depth, YOKE_MU)
    g.add_iron(-TRACK_W / 2, -TRACK_W / 2 + pole_w, -slot_depth, 0.0, YOKE_MU)
    g.add_iron(TRACK_W / 2 - pole_w, TRACK_W / 2, -slot_depth, 0.0, YOKE_MU)
    g.add_current(-TRACK_W / 2 + pole_w + 0.02, TRACK_W / 2 - pole_w - 0.02, -slot_depth, 0.0, NI)
    g.add_current(-TRACK_W / 2 + pole_w + 0.02, TRACK_W / 2 - pole_w - 0.02, -slot_depth - yoke_t - 0.3, -slot_depth - yoke_t - 0.05, -NI)
    g.solve()
    return evaluate(g, [(-TRACK_W / 2, -TRACK_W / 2 + pole_w), (TRACK_W / 2 - pole_w, TRACK_W / 2)],
                    [(-TRACK_W / 2 + pole_w + 0.02, TRACK_W / 2 - pole_w - 0.02, -slot_depth, 0.0)], "I4 iron C-core, 200 mm poles at the edges")


# ---------------------------------------------------------------- casing support tracks (40 coils of 25 mm, 52 turns x 400 A)
NI_SUPPORT = lc.EML_COILS_PER_METER * lc.EML_TURNS_PER_COIL * lc.EML_I_OPERATING     # 832,000 ampere-turns per 1 m track


def case_support_flat(NI=NI_SUPPORT):
    """S1: 40 flat racetrack coils 25 mm wide side by side, one polarity (go and return legs 25 mm apart)."""
    g = base_grid()
    n = 40; wc = TRACK_W / n
    for k in range(n):
        x0 = -TRACK_W / 2 + k * wc
        g.add_current(x0, x0 + 0.010, -0.025, 0.0, NI / n)
        g.add_current(x0 + wc - 0.010, x0 + wc, -0.025, 0.0, -NI / n)
    g.solve()
    r = evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.025, 0.0)], "S1 support, 40 flat coils, one polarity")
    r["NI"] = NI
    return r


def case_support_edge(NI=NI_SUPPORT, back=0.025):
    """S2: 40 coils on edge, the return leg `back` m behind the plate-side leg: a 1 m sheet with a close return."""
    g = base_grid()
    n = 40; wc = TRACK_W / n
    for k in range(n):
        x0 = -TRACK_W / 2 + k * wc
        g.add_current(x0, x0 + wc, -0.025, 0.0, NI / n)
        g.add_current(x0, x0 + wc, -back - 0.050, -back - 0.025, -NI / n)
    g.solve()
    r = evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.025, 0.0), (-TRACK_W / 2, TRACK_W / 2, -back - 0.050, -back - 0.025)],
                 f"S2 support, 40 coils on edge, return {back*1e3:.0f} mm behind")
    r["NI"] = NI
    return r


def case_support_sheet_far(NI=NI_SUPPORT):
    """S3: the Chapter 7 picture, a 1 m current sheet with its return far away (mu0 K / 2 for an infinite sheet)."""
    g = ms.Grid(4.0, 3.0, 0.01)
    g.add_iron(-PLATE_W / 2, PLATE_W / 2, GAP, GAP + PLATE_T, PLATE_MU)
    g.add_current(-TRACK_W / 2, TRACK_W / 2, -0.025, 0.0, NI)
    g.add_current(-TRACK_W / 2, TRACK_W / 2, -1.225, -1.20, -NI)
    g.solve()
    r = evaluate(g, [(-TRACK_W / 2, TRACK_W / 2)], [(-TRACK_W / 2, TRACK_W / 2, -0.025, 0.0)], "S3 support, 1 m sheet, return 1.2 m behind")
    r["NI"] = NI
    return r


def case_support_iron(NI=NI_SUPPORT):
    """S4: the support ampere-turns on iron poles, 4 poles of 200 mm (the I1 geometry), NI per pole = NI_track / 4."""
    r = case_iron_poles(0.25, 0.20, "S4 support, iron, 3 full + 2 half poles of 200 mm", NI=NI / 5, end_half=True)
    r["NI"] = NI
    return r


# ---------------------------------------------------------------- general track builder (2026-10-12 task)
def build_track(W=1.0, n_full=1, pole_w=0.30, slot_w=0.20, slot_depth=0.30, yoke_t=0.20, plate_t=0.10, plate_over=0.20,
                plate_split=False, split_gap=0.05, gap=GAP, h=0.01, NI=NI_REF, domain_x=None, domain_y=None,
                plate_mu=PLATE_MU, plate_bsat=PLATE_B_SAT, yoke_mu=YOKE_MU, yoke_bsat=YOKE_B_SAT):
    """A track of width W facing a plate across `gap`: n_full full poles of width pole_w and a half pole at each edge,
    slots of width slot_w between them (the geometry is consistent when W = (n_full + 1) (pole_w + slot_w)); the
    windings fill the slots to slot_depth; a yoke of yoke_t below.  Each pole is wound with NI, so an inner slot carries
    2 NI and the outer slots NI.  The plate is plate_over wider than the track on each side, plate_t thick; with
    plate_split it is cut into one strip per pole pair (the split gaps sit over the full poles' centers for an even number
    of pole pairs, else over every second slot), split_gap wide.  Returns (grid, pole faces, slot boxes, meta)."""
    pitch = pole_w + slot_w
    W_eff = (n_full + 1) * pitch
    plate_w = W_eff + 2 * plate_over
    Lx = domain_x or max(3.0, 2.0 * plate_w + 1.0)
    Ly = domain_y or max(1.6, 2 * (slot_depth + yoke_t + gap + plate_t) + 0.6)
    g = ms.Grid(Lx, Ly, h)
    # plate
    y_lo, y_hi = gap, gap + plate_t
    strips = []
    if plate_split and n_full >= 1:
        # strips between the slot centers: one strip per pole pair means n_pairs = (n_full + 1) / 2 strips
        n_poles = n_full + 2
        n_strips = max(1, (n_full + 1) // 2)
        edges = np.linspace(-plate_w / 2, plate_w / 2, n_strips + 1)
        for k in range(n_strips):
            a = edges[k] + (split_gap / 2 if k > 0 else 0.0); b = edges[k + 1] - (split_gap / 2 if k < n_strips - 1 else 0.0)
            g.add_iron(a, b, y_lo, y_hi, plate_mu, B_sat=plate_bsat); strips.append((a, b))
    else:
        g.add_iron(-plate_w / 2, plate_w / 2, y_lo, y_hi, plate_mu, B_sat=plate_bsat); strips.append((-plate_w / 2, plate_w / 2))
    # poles and slots
    faces = []
    x = -W_eff / 2
    faces.append((x, x + pole_w / 2)); x += pole_w / 2 + slot_w
    for k in range(n_full):
        faces.append((x, x + pole_w)); x += pole_w + slot_w
    faces.append((x, x + pole_w / 2))
    slot_edges = [(-W_eff / 2 - slot_w, -W_eff / 2)] + [(faces[k][1], faces[k + 1][0]) for k in range(len(faces) - 1)] + [(W_eff / 2, W_eff / 2 + slot_w)]
    x_min = slot_edges[0][0] - 0.05; x_max = slot_edges[-1][1] + 0.05
    g.add_iron(x_min, x_max, -slot_depth - yoke_t, -slot_depth, yoke_mu, B_sat=yoke_bsat)
    for (a, b) in faces:
        g.add_iron(a, b, -slot_depth, 0.0, yoke_mu, B_sat=yoke_bsat)
    slots = []
    n_slots = len(slot_edges)
    for k, (x_lo, x_hi) in enumerate(slot_edges):
        if x_hi - x_lo < h:
            continue
        sgn = 1 if k % 2 == 0 else -1
        amp = NI if k in (0, n_slots - 1) else 2 * NI
        g.add_current(x_lo, x_hi, -slot_depth, 0.0, sgn * amp)
        slots.append((x_lo, x_hi, -slot_depth, 0.0))
    n_groups = len(faces)
    iron_casing = ((x_max - x_min) * yoke_t + sum(b - a for a, b in faces) * slot_depth)        # m^2 of yoke and poles per m
    iron_cable = sum(b - a for a, b in strips) * plate_t                                        # m^2 of plate per m
    meta = dict(W=W_eff, n_full=n_full, n_groups=n_groups, pole_w=pole_w, slot_w=slot_w, slot_depth=slot_depth, yoke_t=yoke_t,
                plate_t=plate_t, plate_w=plate_w, plate_split=plate_split, strips=strips, gap=gap, h=h, NI_total=n_groups * NI,
                iron_casing_m2=iron_casing, iron_cable_m2=iron_cable, winding_m2=sum((b - a) * slot_depth for a, b, _, _ in slots))
    return g, faces, slots, meta


def evaluate_track(g, faces, slots, meta, label=""):
    """Fields and force for a built track: like evaluate() but using the geometry's own gap, plate and width."""
    gap, plate_t, plate_w, W = meta["gap"], meta["plate_t"], meta["plate_w"], meta["W"]
    h = g.h
    y_gap = gap - 2 * h
    xs, By, Bx = g.line_By(y_gap, -W / 2, W / 2)
    B_track = float(np.abs(By).mean())
    pole_sel = np.zeros_like(xs, dtype=bool)
    for (a, b) in faces:
        pole_sel |= (xs >= a) & (xs <= b)
    B_pole = float(np.abs(By[pole_sel]).mean()) if pole_sel.any() else B_track
    F = g.force_on_box(-plate_w / 2 - 0.1, plate_w / 2 + 0.1, y_gap, gap + plate_t + 2 * h)
    cuts = [0.5 * (faces[k][1] + faces[k + 1][0]) for k in range(len(faces) - 1)]
    j_lo = int(round((gap - g.y0) / h)) + 1; j_hi = int(round((gap + plate_t - g.y0) / h)) - 1
    B_plate = 0.0
    for xc in cuts:
        i = int(round((xc - g.x0) / h))
        B_plate = max(B_plate, float(np.abs(g.Bx[i, j_lo:j_hi]).mean()))
    W_mag = float(0.5 * np.sum(g.A * g.J) * h * h)
    wf = [g.region_field(*wb) for wb in slots]
    m_pl = np.zeros_like(g.mu_r, dtype=bool)
    for (a, b) in meta["strips"]:
        m_pl |= g.box(a, b, gap, gap + plate_t)
    Bm = np.sqrt(g.Bx[m_pl]**2 + g.By[m_pl]**2)
    return dict(label=label, B_track=B_track, B_pole=B_pole, F=F, B_plate=B_plate, W_mag=W_mag,
                Bx_w=max(w["Bx_max"] for w in wf), By_w=max(w["By_max"] for w in wf),
                Bx_wm=float(np.mean([w["Bx_mean"] for w in wf])), By_wm=float(np.mean([w["By_mean"] for w in wf])),
                B_plate_max=float(Bm.max()) if Bm.size else 0.0, plate_sat_frac=float(np.mean(g.mu_r[m_pl] < 0.1 * PLATE_MU)) if Bm.size else 0.0,
                pole_frac=float(pole_sel.mean()), **meta)


def track_curve(levels, verbose=False, **geom):
    """Saturating force curve for a built track: the current stepped up through `levels` (ampere-turns per pole) with the
    permeability carried from level to level.  Returns the list of evaluate_track() dicts plus a linear solve at NI_REF."""
    g, faces, slots, meta = build_track(NI=NI_REF, **geom)
    J_ref = g.J.copy()
    g.solve()
    lin = evaluate_track(g, faces, slots, meta, "linear"); lin["NI"] = NI_REF; lin["nonlinear"] = False
    out = [lin]
    mu_lin = g.mu_r.copy()
    for NI in levels:
        g.J = J_ref * (NI / NI_REF)
        box = (-meta["plate_w"] / 2 - 0.1, meta["plate_w"] / 2 + 0.1, meta["gap"] - 2 * g.h, meta["gap"] + meta["plate_t"] + 2 * g.h)
        g.solve_nonlinear(force_box=box)
        r = evaluate_track(g, faces, slots, meta, f"{NI/1e3:.0f} kA"); r["NI"] = NI; r["nonlinear"] = True
        r["outer"] = g.outer_iters; r["outer_dF"] = g.outer_dF; r["F_lin"] = lin["F"] * (NI / NI_REF) ** 2
        out.append(r)
        if verbose:
            print(f"      {NI/1e3:5.0f} kA: F {r['F']/1e3:6.1f} kN/m (linear {r['F_lin']/1e3:6.0f}), B_track {r['B_track']:.3f}, plate max {r['B_plate_max']:.2f} T, "
                  f"sat {r['plate_sat_frac']*100:3.0f} %, winding {r['Bx_wm']:.2f}/{r['By_wm']:.2f} T, W {r['W_mag']/1e3:6.1f} kJ/m ({r['outer']} outer, dF {r['outer_dF']:.4f})", flush=True)
    return out


def run_support(verbose=True):
    cases = []
    for fn in (case_support_flat, case_support_edge, case_support_sheet_far, case_support_iron):
        cases.append(fn())
        if verbose:
            print(f"    solved {cases[-1]['label']} ({cases[-1]['iters']} iterations, {cases[-1]['t']:.0f} s)")
    return cases


GEOMS = {"H1": (0.25, 0.20, dict(end_half=True)), "H3": (0.3333, 0.20, dict(end_half=True)), "H4": (0.3333, 0.15, dict(end_half=True)),
         "H6": (0.50, 0.25, dict(end_half=True)), "H7": (0.50, 0.30, dict(end_half=True)),
         "H3t": (0.3333, 0.20, dict(end_half=True, plate_t=0.20)), "H6t": (0.50, 0.25, dict(end_half=True, plate_t=0.20)),
         "H7t": (0.50, 0.30, dict(end_half=True, plate_t=0.20))}
NL_LEVELS = (25e3, 50e3, 75e3, 100e3, 150e3, 200e3, 300e3)     # ampere-turns per pole for the saturating solves (continuation)
NL_H = 0.01                                           # cell size for the saturating solves (10 mm: 10 cells through the plate)


def run_nonlinear(verbose=True):
    """Saturating-iron force curves for the leading geometries on the 10 mm grid: one grid per geometry, the current stepped
    up through NL_LEVELS with the permeability distribution and the potential carried from level to level (continuation),
    plus the linear solve on the same grid for the grid-size comparison, plus the support track at its installed ampere-turns."""
    global H, PLATE_T
    H_saved = H
    H = NL_H
    out = []
    t0 = time.time()
    try:
        for name, (pitch, pw, kw) in GEOMS.items():
            lin = case_iron_poles(pitch, pw, f"{name} linear 10 mm grid", NI=NI_REF, **kw)
            out.append(lin)
            g, faces, slots, meta = case_iron_poles(pitch, pw, name, NI=NI_REF, nonlinear="grid", **kw)
            J_ref = g.J.copy()
            pt_saved = PLATE_T; PLATE_T = meta["plate_t"]
            for NI in NL_LEVELS:
                g.J = J_ref * (NI / NI_REF)
                g.solve_nonlinear(force_box=(-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, GAP - 2 * H, GAP + PLATE_T + 2 * H))
                r = finish_case(g, faces, slots, f"{name} saturating, {NI/1e3:.0f} kA per pole", NI, meta, True)
                r["geom"] = name
                out.append(r)
                if verbose:
                    print(f"    {r['label']}: F {r['F']/1e3:.0f} kN/m (linear {lin['F']*(NI/NI_REF)**2/1e3:.0f}), B_track {r['B_track']:.3f}, plate max {r['B_plate_max']:.2f} T, "
                          f"plate saturated {r['plate_sat_frac']*100:.0f} %, winding Bx/By mean {r['Bx_wm']:.2f}/{r['By_wm']:.2f} T ({r['outer']} outer, mu change {r['outer_change']:.4f}, "
                          f"F change {r['outer_dF']:.4f}, {g.solve_time:.0f} s)", flush=True)
            PLATE_T = pt_saved
        g, faces, slots, meta = case_iron_poles(0.25, 0.20, "S4", NI=NI_REF, nonlinear="grid", end_half=True)
        J_ref = g.J.copy()
        for NI_track in (NI_SUPPORT / 2, NI_SUPPORT):
            NI = NI_track / 5
            g.J = J_ref * (NI / NI_REF)
            g.solve_nonlinear(force_box=(-PLATE_W / 2 - 0.1, PLATE_W / 2 + 0.1, GAP - 2 * H, GAP + PLATE_T + 2 * H))
            r = finish_case(g, faces, slots, f"S4 support saturating, {NI_track/1e3:.0f} kA per track", NI, meta, True)
            r["geom"] = "S4"; r["NI_track"] = NI_track
            out.append(r)
            if verbose:
                print(f"    {r['label']}: F {r['F']/1e3:.0f} kN/m, B_track {r['B_track']:.3f} T, plate max {r['B_plate_max']:.2f} T ({r['outer']} outer, mu change {r['outer_change']:.4f}, "
                      f"F change {r['outer_dF']:.4f}, {g.solve_time:.0f} s)", flush=True)
    finally:
        H = H_saved
    if verbose:
        print(f"    saturating solves done ({time.time()-t0:.0f} s)")
    return out


def run_all(verbose=True):
    cases = []
    t0 = time.time()
    for fn in (case_air_flat_same, case_air_edge, case_air_flat_alt, case_air_sheet_far):
        cases.append(fn())
        if verbose:
            print(f"    solved {cases[-1]['label']} ({cases[-1]['iters']} iterations, {cases[-1]['t']:.0f} s)")
    for pitch, pw, lab, kw in ((0.25, 0.20, "I1 iron, 4 full poles of 200 mm, slots 50 mm", {}),
                               (0.50, 0.40, "I2 iron, 2 full poles of 400 mm, slots 100 mm", {}),
                               (0.125, 0.10, "I3 iron, 8 full poles of 100 mm, slots 25 mm", {}),
                               (0.25, 0.20, "H1 iron, 3 full + 2 half poles of 200 mm, slots 50 mm", dict(end_half=True)),
                               (0.25, 0.125, "H2 iron, 3 full + 2 half poles of 125 mm, slots 125 mm", dict(end_half=True)),
                               (0.3333, 0.20, "H3 iron, 2 full + 2 half poles of 200 mm, slots 133 mm", dict(end_half=True)),
                               (0.3333, 0.15, "H4 iron, 2 full + 2 half poles of 150 mm, slots 183 mm", dict(end_half=True)),
                               (0.3333, 0.25, "H5 iron, 2 full + 2 half poles of 250 mm, slots 83 mm", dict(end_half=True)),
                               (0.50, 0.25, "H6 iron, 1 full + 2 half poles of 250 mm, slots 250 mm", dict(end_half=True)),
                               (0.50, 0.30, "H7 iron, 1 full + 2 half poles of 300 mm, slots 200 mm", dict(end_half=True)),
                               (0.50, 0.20, "H8 iron, 1 full + 2 half poles of 200 mm, slots 300 mm", dict(end_half=True)),
                               (0.3333, 0.20, "H3t iron, as H3 with a 200 mm plate", dict(end_half=True, plate_t=0.20)),
                               (0.50, 0.25, "H6t iron, as H6 with a 200 mm plate", dict(end_half=True, plate_t=0.20)),
                               (0.50, 0.30, "H7t iron, as H7 with a 200 mm plate", dict(end_half=True, plate_t=0.20))):
        cases.append(case_iron_poles(pitch, pw, lab, **kw))
        if verbose:
            print(f"    solved {cases[-1]['label']} ({cases[-1]['iters']} iterations, {cases[-1]['t']:.0f} s)")
    cases.append(case_iron_ccore())
    if verbose:
        print(f"    solved {cases[-1]['label']} ({cases[-1]['iters']} iterations, {cases[-1]['t']:.0f} s); all {time.time()-t0:.0f} s")
    return cases


# The solved cases, tracked in data/ so the scripts that read them need no 4 minute solve; running this file recomputes it.
JSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "compensation_fields.json")


def _without_times(res):
    """The cases without their solve-time field, so the tracked cache does not change from one recompute to the next."""
    return {k: [{kk: vv for kk, vv in c.items() if kk != "t"} for c in v] for k, v in res.items()}


def load_or_run(verbose=True):
    """The solved cases (compensation and support), from the JSON of the last run of this file if it exists."""
    if os.path.exists(JSON_PATH):
        with open(JSON_PATH, encoding="utf-8") as f:
            return json.load(f)
    comp = run_all(verbose); sup = run_support(verbose); nl = run_nonlinear(verbose)
    res = dict(compensation=comp, support=sup, nonlinear=nl)
    os.makedirs(os.path.dirname(JSON_PATH), exist_ok=True)
    with open(JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(_without_times(res), f, indent=1, default=float)
    return res


def print_cases(cases):
    for c in cases:
        print(f"{c['label']:52} NI {c['NI']/1e3:7.1f} kA: B_pole {c['B_pole']:.4f} T  B_track {c['B_track']:.4f} T  F {c['F']/1e3:7.1f} kN/m  "
              f"B_plate {c['B_plate']:.3f} T  W_mag {c['W_mag']:9.1f} J/m  winding Bx {c['Bx_w']:.3f}/{c['Bx_wm']:.3f}  By {c['By_w']:.3f}/{c['By_wm']:.3f} T (max/mean)")


if __name__ == "__main__":
    if "--help" in sys.argv[1:] or "-h" in sys.argv[1:]:
        print(__doc__); sys.exit(0)
    if "--add-nonlinear" in sys.argv and os.path.exists(JSON_PATH):
        # keep the saved linear cases, add (or redo) the saturating solves
        with open(JSON_PATH, encoding="utf-8") as f:
            res = json.load(f)
        res["nonlinear"] = run_nonlinear()
        with open(JSON_PATH, "w", encoding="utf-8") as f:
            json.dump(_without_times(res), f, indent=1, default=float)
    else:
        if os.path.exists(JSON_PATH):
            os.remove(JSON_PATH)
        res = load_or_run()
    print()
    print("Compensation tracks, per 100 kA-turns per pole (per track for the air cases); linear in the current:")
    print_cases(res["compensation"])
    print()
    print("Casing support tracks at their design ampere-turns (40 coils x 52 turns x 400 A = 832 kA per 1 m track):")
    print_cases(res["support"])
    print()
    print("Saturating iron (plate 1.2 T, yoke 1.8 T), 10 mm grid; the first row of each geometry is the linear solve on the same grid for the grid check:")
    for c in res["nonlinear"]:
        extra = (f"  plate max {c['B_plate_max']:.2f} T, saturated {c['plate_sat_frac']*100:.0f} %, outer {c['outer']}, mu change {c['outer_change']:.4f}, F change {c['outer_dF']:.4f}"
                 if c.get("nonlinear") else "")
        print(f"{c['label']:50} NI {c['NI']/1e3:6.1f} kA: B_pole {c['B_pole']:.3f} T  B_track {c['B_track']:.3f} T  F {c['F']/1e3:6.1f} kN/m  W_mag {c['W_mag']/1e3:7.1f} kJ/m  "
              f"winding Bx {c['Bx_w']:.2f}/{c['Bx_wm']:.2f}  By {c['By_w']:.2f}/{c['By_wm']:.2f} T{extra}")
    print()
    print(f"Saved {JSON_PATH}")
