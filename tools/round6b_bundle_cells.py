#!/usr/bin/env python3
"""Round 6b, Section 2.2: the tape loss with each cell of the 200 x 50 mm bundle lying the better way (run from Python/orbitalring after bearing_round6b.py --part=bundle)."""
import json, os, sys, math
sys.path.insert(0, os.getcwd())
import numpy as np
import bearing_round6b as b, bearing_round5 as r5, bearing_electrical as be
done = json.load(open("runs/round6b_logs/round6b_bundle.json"))
theva, _ = be.load_theva()
out = {}
for (Wp, T, Wl) in ((1.0, 0.15, 0.6), (1.0, 0.20, 0.6)):
    res = done[f"{Wp:.2f}|{T:.3f}|{Wl:.2f}|0.200|0.050|0.0050"]
    fields, s, tr = b.leg_fields(Wp, T, Wl, res["NI_req"], 0.005)
    li = dict(b.launch_profile(res["W_mag"], be.load_records()[1]["F_req"]), F_req=be.load_records()[1]["F_req"])
    fh = res["faces"]["h"]
    for (ncx, ncy) in ((10, 2), (4, 5)):
        mixed = []; n_v = 0; pattern = []
        for (x, sgn), (xs, ys), (Bx, By) in zip(tr.legs, tr.lattices, fields):
            ix = np.clip(((xs - (x - 0.1)) / 0.2 * ncx).astype(int), 0, ncx - 1)
            iy = np.clip(((ys - ys.min()) / (ys.max() - ys.min() + 1e-12) * ncy).astype(int), 0, ncy - 1)
            Bx2, By2 = Bx.copy(), By.copy(); pat = [["h"] * ncx for _ in range(ncy)]
            for cx in range(ncx):
                for cy in range(ncy):
                    m = (ix == cx) & (iy == cy)
                    if (np.abs(Bx[m])**2).mean() < (np.abs(By[m])**2).mean():      # less field through the face if the faces are vertical
                        Bx2[m], By2[m] = By[m], Bx[m]; n_v += 1; pat[cy][cx] = "v"
            mixed.append((Bx2, By2)); pattern.append(pat)
        ic, B, th = b.worst_ic(mixed, "h", theva)
        tapes = res["NI_req"] / (0.5 * ic); layers = (res["NI_req"] / 20) / (0.5 * ic) / 10
        Ep, Es = b.tape_loss_launch(mixed, "h", tapes, layers, li)
        out[f"{Wp}x{T}|{ncx}x{ncy}"] = dict(n_v=n_v, n_cells=2 * ncx * ncy, ic=ic, B=B, th=th, Ep=Ep, Es=Es, tape=2 * tapes * 0.003 * 8e-5 * 0 + 2 * tapes * 0.003 * __import__("bearing_design").TAPE_T * __import__("bearing_design").TAPE_RHO, pattern=pattern)
        print(f"plate {Wp} x {T}: cells {ncx} x {ncy}: {n_v} of {2*ncx*ncy} cells faces vertical; worst tape {B:.2f} T at {th:.0f} deg, {ic:.0f} A; loss {(Ep+Es)/1e9:.0f} GJ ({Ep/1e9:.0f} + {Es/1e9:.0f}); all horizontal {(fh['Ep']+fh['Es'])/1e9:.0f}")
        for leg, pat in zip(tr.legs, pattern):
            print(f"   leg at {leg[0]*1e3:+.0f} mm: top to bottom " + " / ".join("".join(r) for r in reversed(pat)))
b.save("round6b_bundle_cells.json", out)
