#!/usr/bin/env python3
"""
First check of the tape current in the deployment LIM coil against the THEVA 70 K limit for three winding layouts and 1 to 12 soldered layers of 3 mm tape (a 2-D solver check is a later task).
Book: Orbital Ring Engineering, Volume I: Chapter 9 (Table 9.2, the tape current rows; Section 9.5).
Run: python tools/lim_tape_current_check.py   (reads the THEVA file in data/ through bearing_design.py)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE, ORBITAL_RING_HW_KGM (the ring, through ring_altitude.py)

First check of the tape current in the deployment LIM coil (Chapter 9, Table 9.2 and Section 9.5).
Coil: 208 turns at 163 A peak (lim_config.py).  The code has no winding layout, so two limiting layouts are tried:
  A  one pancake: the turns stacked face to face, bundle 3 mm wide (the tape width) and 208 turn pitches high
  B  one flat sheet: the turns laid edge to edge, 208 tape widths wide and one conductor thick
  C  eight pancakes side by side, 26 turns each
Each turn is N layers of 3 mm tape soldered together (0.08 mm per layer) plus 0.05 mm of turn insulation.
Field: 2-D, straight conductor, uniform current density in each turn, return side of the coil 2 m away ignored (a few mT).
Tape limit: THEVA 70 K curve, lower of the 0 and 90 degree curves at the peak field in the winding, as bearing_design.py does.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))   # the repository root, for the imports below
import os, sys, math
import numpy as np
import bearing_design as bd
MU0 = 4e-7 * math.pi
theva = bd.read_theva()
N_TURNS, I_COIL, W_TAPE, T_LAYER, T_INS = 208, 163.0, 3e-3, 0.08e-3, 0.05e-3

def peak_field(rects, I_each, n_sub=(6, 2)):
    """rects: list of (x0, y0, w, h) conductor rectangles, each carrying I_each. Returns peak |B| at the conductor centres and corners (T)."""
    src = []; pts = []
    for (x0, y0, w, h) in rects:
        nx = max(2, int(round(w / 0.5e-3)) if w > h else n_sub[1]); ny = max(2, int(round(h / 0.5e-3)) if h >= w else n_sub[1])
        xs = x0 + (np.arange(nx) + 0.5) * w / nx; ys = y0 + (np.arange(ny) + 0.5) * h / ny
        X, Y = np.meshgrid(xs, ys); src.append(np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, I_each / X.size)]))
        for px in (x0, x0 + w / 2, x0 + w):
            for py in (y0, y0 + h / 2, y0 + h):
                pts.append((px, py))
    src = np.vstack(src); pts = np.array(pts)
    Bmax = 0.0
    for i in range(0, len(pts), 400):
        P = pts[i:i + 400]
        dx = P[:, None, 0] - src[None, :, 0]; dy = P[:, None, 1] - src[None, :, 1]
        r2 = dx * dx + dy * dy + (0.03e-3) ** 2
        bx = (-MU0 / (2 * math.pi) * src[None, :, 2] * dy / r2).sum(1); by = (MU0 / (2 * math.pi) * src[None, :, 2] * dx / r2).sum(1)
        Bmax = max(Bmax, float(np.sqrt(bx * bx + by * by).max()))
    return Bmax

def layout(kind, n_layers):
    t = n_layers * T_LAYER; pitch = t + T_INS
    if kind == "A":   # tape width along x, stack along y
        return [(0.0, k * pitch, W_TAPE, t) for k in range(N_TURNS)], "one pancake, %.0f mm high" % (N_TURNS * pitch * 1e3)
    if kind == "B":   # tapes edge to edge along x
        return [(k * (W_TAPE + 0.2e-3), 0.0, W_TAPE, t) for k in range(N_TURNS)], "flat sheet, %.0f mm wide" % (N_TURNS * (W_TAPE + 0.2e-3) * 1e3)
    r = []
    for p in range(8):
        r += [(p * (W_TAPE + 1e-3), k * pitch, W_TAPE, t) for k in range(26)]
    return r, "8 pancakes of 26 turns, %.0f mm wide, %.1f mm high" % (8 * (W_TAPE + 1e-3) * 1e3, 26 * pitch * 1e3)

print("Tape at 70 K, 3 mm wide: %.0f A with no field (THEVA), 200 A in lim_config.py" % bd.ic_tape(theva, 0.0, W_TAPE))
print("Coil: %d turns x %.0f A = %.0f ampere-turns\n" % (N_TURNS, I_COIL, N_TURNS * I_COIL))
for kind in "ABC":
    print("Layout %s" % kind)
    print("  layers  layout                                   peak field   limit per tape   current per tape   share")
    done50 = done80 = None
    for n in range(1, 13):
        rects, desc = layout(kind, n)
        B = peak_field(rects, I_COIL)
        ic = bd.ic_tape(theva, B, W_TAPE); i_t = I_COIL / n; share = i_t / ic
        print("  %4d    %-40s %6.3f T     %6.0f A         %6.1f A        %5.0f %%" % (n, desc, B, ic, i_t, share * 100))
        if done80 is None and share <= 0.80: done80 = n
        if done50 is None and share <= 0.50: done50 = n; break
    print("  -> layers for 80%%: %s ; for 50%%: %s ; tape per LIM coil turn-length multiplies by the layer count\n" % (done80, done50))
