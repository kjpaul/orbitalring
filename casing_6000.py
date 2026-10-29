#!/usr/bin/env python3
"""
The ring with a 6,000 kg/m casing, for Chapter 16's Table 16.2 (2026-10-28): the Chapter 3 quadratic for three hardware cases,
the levitation plates and the bearing tape sized for the halved casing weight with the field solver of bearing_optimize.py,
the LIM reaction plate judged from the deployment runs, the anchor attachments of Chapter 11 and the standing cryogenic load
of Chapter 14, and Table 16.2 assembled from the deployment logs of run_700_casing_6000.sh.
Book: Orbital Ring Engineering, Volume I: Chapter 16 (Section 16.6, Table 16.2); Chapters 3, 7, 9, 11 and 14 for the rules.
Run: python casing_6000.py --bearing            (the plate and tape sweep, about fifteen minutes; writes runs/casing_6000/bearing_sweep.json)
     python casing_6000.py --table <log folder>  (the ring cases, the LIM plate, the rest, and Table 16.2 from the deployment logs)
Environment: ORBITAL_RING_ALT_KM, ORBITAL_RING_CABLE (the ring, through ring_altitude.py); ORBITAL_RING_OUT (output folder, default runs)
"""
import os, sys, re, json, math, time
import ring_altitude as ring

OUT = ring.out_dir("casing_6000")
SWEEP = os.path.join(OUT, "bearing_sweep.json")
M_CASING_NEW = 6000.0
M_CASING_DESIGN = 12000.0
# the cable hardware as ring_altitude.py and lim_config.py count it: four bearing plates 1.4 m wide of sintered iron and two LIM
# reaction plates 1.2 m wide of gamma-TiAl, the thickness of each the variable
RHO_FE, W_LEV, N_LEV = 7870.0, 1.4, 4
RHO_TIAL, W_LIM, N_LIM = 3900.0, 1.2, 2
T_LEV_DESIGN, T_LIM_DESIGN = 0.20, 0.20
THICKNESSES = (0.20, 0.175, 0.15, 0.125, 0.10, 0.075, 0.05)
SAT_LIMIT = 0.12                      # bearing_optimize Part 2: the plate may not be saturated beyond 12 percent at the requirement
RATE_MARGIN = 1.13                    # bearing_optimize: the rating at the 70 percent share must be at least 1.13 x the requirement
SHARE = 0.7
CRYO_FACTOR = 1.0                     # formers, cryostat and supports at one tape mass (bearing_mass.py) [UNVERIFIED]
# Chapter 11's attachments: a 30 m attachment may add a tenth of the casing weight at any one place; the trolley line pulls
# 2.24 MN at the ring (four branches at 0.56 MN) and the spreader arm 23 MN (Chapter 11, Section 11.4)
L_ATTACH, F_TROLLEY, F_ARM = 30.0, 2.24e6, 23e6
# Chapter 14: 5 W per square metre of levitation coil facing the cable, 4 m of coil width, 500 m per site
Q_COIL_PER_M2, LEV_COIL_WIDTH, SITE_SPACING = 5.0, 4.0, 500.0


def hardware(t_lev, t_lim):
    return N_LEV * RHO_FE * t_lev * W_LEV + N_LIM * RHO_TIAL * t_lim * W_LIM


def fmt(x, nd=0):
    return f"{x:,.{nd}f}"


def ring_row(label, m_load, m_hw):
    d = ring.design(m_load=m_load, m_hw=m_hw)
    print(f"  {label:<58} {m_hw:>9,.0f} {d['m_struct']:>10,.0f} {d['m_cable']:>10,.0f} {d['v_cable']:>8,.0f} {d['v_rel']:>8,.0f} {ring.V_GROUND_SYNC:>8,.0f} "
          f"{d['tension']/1e9:>8.1f} {d['area']:>7.2f} {d['side']:>6.2f} {d['E_deploy']/1e18:>7.2f} {d['M_ring']/1e12:>7.3f}")
    return d


# ======================================================================================================================
# the bearing sweep: the support track at the halved casing weight (standing load kept or halved) and the lift track,
# against the plate thickness, on bearing_optimize's evaluator at 5 mm cells
# ======================================================================================================================
def bearing_sweep():
    import bearing_optimize as bo
    import bearing_design as bd
    import numpy as np
    theva4 = bd.read_theva(); theva77 = bo.read_theva_file(bo.THEVA_77)
    rob = bo.read_robinson_angles(); theva = dict(theva4); theva["full"] = rob      # as Part 3b: the full Robinson angle data
    h = 0.005
    g = ring.G_NET
    cache = json.load(open(SWEEP)) if os.path.exists(SWEEP) else {}

    def save():
        json.dump(cache, open(SWEEP, "w"), indent=1)

    def record(kind, label, t, F_req, F_load, geom):
        key = f"{kind}|{label}|{t:.3f}"
        if key in cache:
            return cache[key]
        bo.PLATE_T = t; bo.F_LOAD = F_load
        t0 = time.time()
        ev = bo.Evaluator(geom, h, theva, theva77)
        r = ev.evaluate_rerate(F_req=F_req, kind=kind)
        rated = {s: bo.rated_tape(r, s) for s in (0.6, 0.7, 0.8)}
        out = dict(kind=kind, label=label, t=t, F_req=F_req, F_load=F_load, NI_req=r["NI_req"], sat=r["sat"], B_plate=r["B_plate"],
                   sat13=r["sat13"], sat16=r["sat16"], rated60=rated[0.6], rated70=rated[0.7], rated80=rated[0.8], F_1T=r["F_rated"],
                   seconds=time.time() - t0, solves=r["solves"], track_width=r["track_width"])
        if kind == "support":
            out.update(ic12=r["ic12"], I12=r["I12"], turns12=r["turns12"], tape12=r["tape12"], NI_load=r["NI_load"], ic3=r["ic3"], I3=r["I3"], tape3=r["tape3"],
                       NI_446=r["NI_446"], F_446=r["F_446"], sat_446=r["sat_446"], ic3_446=r["ic3_446"], I3_446=r["I3_446"], tape3_446=r["tape3_446"],
                       F_comb=r["F_comb"], sat_comb=r["sat_comb"], tape3_comb=r["tape3_comb"], W_req=r["W_req"], E_coil=r["E_coil"])
        else:
            out.update(ic3=r["ic3"], I3=r["I3"], turns3=r["turns3"], tape3=r["tape3"], W_req=r["W_req"], E_coil=r["E_coil"])
        out = {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in out.items()}
        cache[key] = out; save()
        return out

    def ok(o):
        return o["sat"] <= SAT_LIMIT and o["rated70"] >= RATE_MARGIN * o["F_req"]

    print("=" * 150)
    print("THE LEVITATION PLATES AND THE BEARING TAPE AT A 6,000 kg/m CASING (bearing_optimize's evaluator, 5 mm cells, the Part 3b rules)")
    print("=" * 150)
    print(f"  Support track: the best geometry of the optimization ({bo.geom_key(bo.BEST_SUPPORT)}), requirement = half the casing weight plus the")
    print(f"  track's share of the standing load (kept at {bo.fb.F_LOAD_TRACK/1e3:.0f} kN/m, or halved); lift track: {bo.geom_key(bo.BEST_LIFT)} at its {bo.LIFT_REQ/1e3:.1f} kN/m")
    print(f"  (the 393 kN/m pull of the design launch is the sled's, not the casing's, so it does not change).  Plate 1.4 m wide, iron saturating at")
    print(f"  {bo.cf.PLATE_B_SAT} T (mu {bo.cf.PLATE_MU:.0f}), gap {bo.GAP*1e3:.0f} mm, cryostat wall {bo.WALL_DESIGN*1e3:.0f} mm.  Rules: plate saturation at the requirement <= {SAT_LIMIT:.0%}; the rating at the")
    print(f"  {SHARE:.0%} share (the pull at which the worst 12 mm tape reaches {SHARE:.0%} of its limit) >= {RATE_MARGIN:.2f} x the requirement; the winding built at 50 percent.")
    print(f"  Casing weight {M_CASING_DESIGN:,.0f} kg/m x {g:.3f} = {M_CASING_DESIGN*g/1e3:.1f} kN/m (design), {M_CASING_NEW:,.0f} kg/m x {g:.3f} = {M_CASING_NEW*g/1e3:.1f} kN/m; half of each per track.")
    print()
    sup_cases = [("design case, 12,000 kg/m, standing load 198 kN/m", M_CASING_DESIGN * g / 2 + 198e3, 198e3, (0.20,)),
                 ("6,000 kg/m, standing load kept at 198 kN/m", M_CASING_NEW * g / 2 + 198e3, 198e3, THICKNESSES),
                 ("6,000 kg/m, standing load halved to 99 kN/m", M_CASING_NEW * g / 2 + 99e3, 99e3, THICKNESSES)]
    hdr = f"  {'plate mm':>8} {'req kN/m':>8} {'NI kA':>6} {'sat %':>6} {'B_pl T':>6} {'rated 60/70/80 kN/m':>22} {'x req':>5} {'Ic12 A':>6} {'I50 A':>5} {'tape12':>7} {'load tape own/446/both':>24} {'sat 446/both %':>14} {'s':>4} {'rule':>5}"
    for label, F_req, F_load, ts in sup_cases:
        print(f"  Support track, {label}: requirement {F_req/1e3:.1f} kN/m, design state {(F_req+F_load)/1e3:.0f} kN/m together")
        print(hdr)
        for t in ts:
            o = record("support", label, t, F_req, F_load, bo.BEST_SUPPORT)
            print(f"  {t*1e3:>8.0f} {o['F_req']/1e3:>8.1f} {o['NI_req']/1e3:>6.0f} {o['sat']*100:>6.1f} {o['B_plate']:>6.2f} {o['rated60']/1e3:>7.1f}/{o['rated70']/1e3:>6.1f}/{o['rated80']/1e3:>6.1f} "
                  f"{o['rated70']/o['F_req']:>5.2f} {o['ic12']:>6.0f} {o['I12']:>5.0f} {o['tape12']:>7.1f} {o['tape3']:>7.1f}/{o['tape3_446']:>7.1f}/{o['tape3_comb']:>7.1f} "
                  f"{o['sat_446']*100:>6.1f}/{o['sat_comb']*100:>6.1f} {o['seconds']:>4.0f} {'ok' if ok(o) else 'FAIL':>5}", flush=True)
        print()
    print(f"  Lift track at {bo.LIFT_REQ/1e3:.1f} kN/m (unchanged):")
    print(f"  {'plate mm':>8} {'req kN/m':>8} {'NI kA':>6} {'sat %':>6} {'B_pl T':>6} {'rated 60/70/80 kN/m':>22} {'x req':>5} {'Ic3 A':>6} {'I50 A':>5} {'tape3':>7} {'sat 1.3/1.6 x %':>16} {'s':>4} {'rule':>5}")
    for t in THICKNESSES:
        o = record("lift", "lift", t, bo.LIFT_REQ, 198e3, bo.BEST_LIFT)
        print(f"  {t*1e3:>8.0f} {o['F_req']/1e3:>8.1f} {o['NI_req']/1e3:>6.0f} {o['sat']*100:>6.1f} {o['B_plate']:>6.2f} {o['rated60']/1e3:>7.1f}/{o['rated70']/1e3:>6.1f}/{o['rated80']/1e3:>6.1f} "
              f"{o['rated70']/o['F_req']:>5.2f} {o['ic3']:>6.0f} {o['I3']:>5.0f} {o['tape3']:>7.1f} {o['sat13']*100:>7.1f}/{o['sat16']*100:>7.1f} {o['seconds']:>4.0f} {'ok' if ok(o) else 'FAIL':>5}", flush=True)
    print()
    choose(cache, verbose=True)


def choose(cache, verbose=False):
    """The thinnest plate that meets both rules, for each requirement; the bearing tape and casing side that go with it."""
    out = {}
    for label in ("design case, 12,000 kg/m, standing load 198 kN/m", "6,000 kg/m, standing load kept at 198 kN/m", "6,000 kg/m, standing load halved to 99 kN/m"):
        sup = sorted([o for o in cache.values() if o["kind"] == "support" and o["label"] == label], key=lambda o: o["t"])
        okk = [o for o in sup if o["sat"] <= SAT_LIMIT and o["rated70"] >= RATE_MARGIN * o["F_req"]]
        lift = sorted([o for o in cache.values() if o["kind"] == "lift"], key=lambda o: o["t"])
        okl = [o for o in lift if o["sat"] <= SAT_LIMIT and o["rated70"] >= RATE_MARGIN * o["F_req"]]
        if not okk or not okl:
            continue
        s = okk[0]; l = okl[0]
        tape = 2 * (s["tape12"] + s["tape3_446"]) + 2 * l["tape3"]
        out[label] = dict(t_sup=s["t"], t_lift=l["t"], tape12=s["tape12"], tape3_446=s["tape3_446"], tape3_lift=l["tape3"], tape=tape, casing_side=tape * (1 + CRYO_FACTOR),
                          sat=s["sat"], rated70=s["rated70"], F_req=s["F_req"], m_lev=N_LEV * RHO_FE * W_LEV * (0.5 * s["t"] + 0.5 * l["t"]), m_sup_plates=2 * RHO_FE * W_LEV * s["t"], m_lift_plates=2 * RHO_FE * W_LEV * l["t"])
        if verbose:
            print(f"  {label}: thinnest support plate meeting both rules {s['t']*1e3:.0f} mm (saturation {s['sat']*100:.1f} %, rated {s['rated70']/1e3:.0f} kN/m = {s['rated70']/s['F_req']:.2f} x), "
                  f"lift plate {l['t']*1e3:.0f} mm (saturation {l['sat']*100:.2f} %, rated {l['rated70']/1e3:.0f} kN/m = {l['rated70']/l['F_req']:.2f} x);")
            print(f"     tape per metre of casing: support {s['tape12']:.1f} + load {s['tape3_446']:.1f} per track, lift {l['tape3']:.1f} per track, two of each = {tape:.1f} kg/m; casing side {tape*(1+CRYO_FACTOR):,.0f} kg/m;")
            print(f"     levitation plates on the cable: two support plates {2*RHO_FE*W_LEV*s['t']:,.0f} kg/m + two lift plates {2*RHO_FE*W_LEV*l['t']:,.0f} kg/m = {out[label]['m_lev']:,.0f} kg/m "
                  f"(four 200 mm plates: {N_LEV*RHO_FE*W_LEV*0.2:,.0f} kg/m)")
    return out


# ======================================================================================================================
# the deployment logs
# ======================================================================================================================
def parse_log(path):
    if not os.path.exists(path):
        return None
    txt = open(path, encoding="utf-8", errors="replace").read()
    def grab(pat, conv=float):
        m = re.search(pat, txt)
        return conv(m.group(1)) if m else None
    return dict(path=path,
                days=grab(r"Deployment time:\s*([\d.]+) days"), energy_EJ=grab(r"Total energy delivered:\s*([\d.]+) EJ"),
                m_struct=grab(r"M_CABLE_STRUCT\s+kg/m\s+([\d.]+)"), m_hw=grab(r"M_HARDWARE\s+kg/m\s+([\d.]+)"), m_cable=grab(r"M_CABLE_TOTAL\s+kg/m\s+([\d.]+)"),
                m_load=grab(r"M_LOAD\s+kg/m\s+([\d.]+)"), power_MW=grab(r"MAX_SITE_POWER\s+MW\s+([\d.]+)"),
                thrust_max_N=grab(r"thrust_max\s+([\d.]+) N"), thrust_when=grab(r"thrust_max\s+[\d.]+ N\s+([\d.]+ \w+)", str),
                p_eddy_max_kW=grab(r"p_eddy_max\s+([\d.]+) kW"), p_heat_max_kW=grab(r"p_heat_max\s+([\d.]+) kW"), current_max_A=grab(r"current_max\s+([\d.]+) A"),
                cable_T_max=grab(r"cable_temp_max\s+([\d.]+) K"), wall_T_max=grab(r"wall_temp_max\s+([\d.]+) K"),
                q_cold_max_kW=grab(r"q_cryo_cold_max\s+([\d.]+) kW"), q_cold_when=grab(r"q_cryo_cold_max\s+[\d.]+ kW\s+([\d.]+ \w+)", str),
                p_cryo_max_MW=grab(r"p_cryo_max\s+([\d.]+) MW"), p_thrust_max_MW=grab(r"p_thrust_max\s+([\d.]+) MW"), rad_width=grab(r"Max cryo radiator width:\s*([\d.]+) m"),
                t_plate_mm=grab(r"T_PLATE:\s*([\d.]+) mm"), v_slip_max=grab(r"v_slip_max\s+([\d.]+) m/s"))


def table(logdir):
    print("=" * 150)
    print("1. THE RING WITH A 6,000 kg/m CASING: THE CHAPTER 3 QUADRATIC AT 700 km, RETROGRADE, 12.5 GPa")
    print("=" * 150)
    print("  Structural cable mass from ring_altitude.cable_mass_structural (momentum balance of cable, hardware and casing with the cable at 12.5 GPa);")
    print("  cable velocity in the inertial frame, relative speed over the casing, casing speed (ground-synchronous); tension and section of the")
    print("  structural CNT at 1,700 kg/m3; deployment energy 1/2 L [m_cable (v_c^2 - v_o^2) + m_casing (v_g^2 - v_o^2)]; ring mass = (cable + hardware + casing) x L.")
    print(f"  {'case':<58} {'hw kg/m':>9} {'struct':>10} {'cable':>10} {'v_c m/s':>8} {'v_rel':>8} {'v_g':>8} {'T GN':>8} {'A m2':>7} {'side':>6} {'E EJ':>7} {'ring Tt':>7}")
    sweep = json.load(open(SWEEP)) if os.path.exists(SWEEP) else {}
    pick = choose(sweep)
    hw_c = {}
    lim_t_c = T_LIM_DESIGN      # the LIM plate: from the deployment runs below (the plate test); 200 mm unless the run shows otherwise
    d0 = ring_row("12,000 kg/m design case (check: 115,005 / 125,691 kg/m, 8,270 m/s, 845.6 GN, 67.65 m2, 18.8 EJ)", M_CASING_DESIGN, hardware(0.20, 0.20))
    da = ring_row("(a) 6,000 kg/m, hardware unchanged (four 200 mm bearing plates, two 200 mm LIM plates)", M_CASING_NEW, hardware(0.20, 0.20))
    db = ring_row("(b) 6,000 kg/m, hardware halved with the casing (100 mm plates throughout)", M_CASING_NEW, hardware(0.10, 0.10))
    dc = {}
    for label, short in (("6,000 kg/m, standing load kept at 198 kN/m", "kept"), ("6,000 kg/m, standing load halved to 99 kN/m", "halved")):
        if label in pick:
            p = pick[label]
            hw = 2 * RHO_FE * W_LEV * p["t_sup"] + 2 * RHO_FE * W_LEV * p["t_lift"] + N_LIM * RHO_TIAL * lim_t_c * W_LIM
            hw_c[short] = hw
            dc[short] = ring_row(f"(c) 6,000 kg/m, plates sized: support {p['t_sup']*1e3:.0f} mm, lift {p['t_lift']*1e3:.0f} mm, LIM {lim_t_c*1e3:.0f} mm, standing load {short}", M_CASING_NEW, hw)
    print(f"  (hardware per metre: four bearing plates 1.4 m wide of iron at {RHO_FE:,.0f} kg/m3, two LIM plates 1.2 m wide of gamma-TiAl at {RHO_TIAL:,.0f} kg/m3; 200 mm throughout = {hardware(0.2, 0.2):,.0f} kg/m)")
    print()
    print("  What scales and what does not.  In the quadratic the hardware is part of the rotor: it moves at the cable's speed, so it supplies momentum and")
    print("  centrifugal lift like the structural cable does, and the structural cable is set by the casing almost alone.  Halving the casing halves it whether the")
    print(f"  plates shrink or not: (a) keeps {da['m_struct']/d0['m_struct']*100:.1f} % of the design case's structural cable, (b) {db['m_struct']/d0['m_struct']*100:.1f} %" +
          (f", (c) {dc['kept']['m_struct']/d0['m_struct']*100:.1f} % with the standing load kept and {dc['halved']['m_struct']/d0['m_struct']*100:.1f} % halved." if dc else "."))
    print("  Heavier hardware lowers the cable velocity a little (more rotor mass carries the casing's momentum at less speed) and adds its own mass to the ring;")
    print("  the deployment energy, which is the casing's momentum deficit, hardly moves.  The hardware's cost is its own mass in the ring total and the LIM plate's")
    print("  thickness in the deployment (section 3), not structural cable.  [QUESTION FOR AUTHOR: Section 16.6 says the cable falls by much less than half when the")
    print("  hardware stays at 10,686 kg/m; the Chapter 3 model does not say that.]")
    print()

    print("=" * 150)
    print("2. THE LEVITATION PLATES AND THE BEARING TAPE (from the sweep; python casing_6000.py --bearing for the full tables)")
    print("=" * 150)
    for label in ("design case, 12,000 kg/m, standing load 198 kN/m", "6,000 kg/m, standing load kept at 198 kN/m", "6,000 kg/m, standing load halved to 99 kN/m"):
        rows = sorted([o for o in sweep.values() if o["kind"] == "support" and o["label"] == label], key=lambda o: o["t"])
        if rows:
            print(f"  {label}: requirement {rows[0]['F_req']/1e3:.1f} kN/m per support track")
            for o in rows:
                print(f"    plate {o['t']*1e3:>3.0f} mm: saturation {o['sat']*100:4.1f} % at the requirement ({o['sat_446']*100:4.1f} % in the design state), rated {o['rated70']/1e3:5.1f} kN/m = {o['rated70']/o['F_req']:.2f} x, "
                      f"tape support {o['tape12']:5.1f} + load {o['tape3_446']:5.1f} kg/m per track  {'ok' if (o['sat'] <= SAT_LIMIT and o['rated70'] >= RATE_MARGIN*o['F_req']) else 'FAIL'}")
    rows = sorted([o for o in sweep.values() if o["kind"] == "lift"], key=lambda o: o["t"])
    if rows:
        print(f"  lift track at {rows[0]['F_req']/1e3:.1f} kN/m:")
        for o in rows:
            print(f"    plate {o['t']*1e3:>3.0f} mm: saturation {o['sat']*100:4.2f} %, rated {o['rated70']/1e3:5.1f} kN/m = {o['rated70']/o['F_req']:.2f} x, tape {o['tape3']:5.1f} kg/m per track  "
                  f"{'ok' if (o['sat'] <= SAT_LIMIT and o['rated70'] >= RATE_MARGIN*o['F_req']) else 'FAIL'}")
    print()
    choose(sweep, verbose=True)
    print()

    print("=" * 150)
    print("3. THE DEPLOYMENT (lim_simulation.py, --m_load=6000, the plates through lim_config.configure; run_700_casing_6000.sh)")
    print("=" * 150)
    logs = {}
    for key, name in (("design 8 MW", "lim_simulation_brandt_flat3mm_N208_8MW.txt"), ("design 16 MW", "lim_simulation_brandt_flat3mm_N208_16MW.txt"),
                      ("a 8 MW", "casing_6000_lim_a_8MW.txt"), ("b 8 MW", "casing_6000_lim_b_8MW.txt"), ("b, LIM plate 200 mm, 8 MW", "casing_6000_lim_b_limplate200_8MW.txt"),
                      ("c kept 8 MW", "casing_6000_lim_c_kept_8MW.txt"), ("c kept 16 MW", "casing_6000_lim_c_kept_16MW.txt"),
                      ("c halved 8 MW", "casing_6000_lim_c_halved_8MW.txt"), ("c halved 16 MW", "casing_6000_lim_c_halved_16MW.txt")):
        logs[key] = parse_log(os.path.join(logdir, name))
    print(f"  {'run':<28} {'hw kg/m':>8} {'cable':>8} {'LIM mm':>6} {'MW':>4} {'days':>6} {'EJ':>6} {'thrust N':>8} {'when':>12} {'eddy kW':>8} {'I A':>6} {'cable K':>7} {'wall K':>6} {'cold kW':>7} {'when':>12} {'cryo MW':>7}")
    for key, L in logs.items():
        if L is None or L["days"] is None:
            print(f"  {key:<28} (log not found or run not finished)"); continue
        print(f"  {key:<28} {L['m_hw']:>8,.0f} {L['m_cable']:>8,.0f} {L['t_plate_mm']:>6.0f} {L['power_MW']:>4.0f} {L['days']:>6.1f} {L['energy_EJ']:>6.2f} {L['thrust_max_N']:>8,.0f} {L['thrust_when']:>12} "
              f"{L['p_eddy_max_kW']:>8.1f} {L['current_max_A']:>6.1f} {L['cable_T_max']:>7.1f} {L['wall_T_max']:>6.1f} {L['q_cold_max_kW']:>7.1f} {L['q_cold_when']:>12} {L['p_cryo_max_MW']:>7.2f}")
    print("  (thrust per LIM; eddy heating per LIM on the first day; cold-side load per site, the two LIMs together; the 12,000 kg/m rows are the reference runs of 2026-10-27)")
    print()
    print("  The LIM reaction plate.  In the thrust model (lim_physics.calc_thrust_model1) the plate is an eddy loop of resistance R = rho tau_p / (d w) with d the")
    print("  smaller of the plate thickness and the skin depth, and inductance L = (mu0 tau_p / pi) ln(tau_p / w).  gamma-TiAl at 75e-8 ohm m, 100 m pitch, 2 m coil:")
    rho, tau, w, d = 75e-8, 100.0, 2.0, 0.2
    Lind = 4e-7 * math.pi * tau / math.pi * math.log(tau / w)
    for dd in (0.2, 0.1):
        R = rho * tau / (dd * w); f_match = R / (2 * math.pi * Lind)
        delta1 = math.sqrt(rho / (math.pi * 4e-7 * math.pi * 0.96)); delta_min = math.sqrt(rho / (math.pi * 4e-7 * math.pi * 5 / 200))
        print(f"    plate {dd*1e3:.0f} mm: R {R*1e6:.0f} micro-ohm, L {Lind*1e6:.0f} micro-H, R = omega L at {f_match:.2f} Hz (slip {f_match*2*tau:.0f} m/s); skin depth {delta1:.2f} m at the first day's 0.96 Hz and {delta_min:.1f} m at the 5 m/s minimum slip,")
    print("    so the whole plate conducts at every slip the controller uses and the thickness enters only through R.  The thermal model carries no plate heat capacity")
    print("    (the plate's temperature rise per pass is 0.01 K in the reference run); the first-day heating is the thrust times the slip the controller chooses, and the")
    print("    cable temperature follows from that heat over a fixed 7.5 m2/m of radiating surface.  The plate's width (1.2 m) does not enter the thrust model, whose eddy")
    print("    loop is the coil's 2 m; nothing in the repository sizes the width, so whether it can fall is not decided here [UNVERIFIED].  The (b) rows above, with the")
    print("    LIM plate at 100 and 200 mm on nearly the same cable, show what halving the thickness does to the deployment time, the peak thrust and the heating.")
    print()

    print("=" * 150)
    print("4. THE REST: ANCHOR ATTACHMENTS, BEARING TAPE, STANDING CRYOGENIC LOAD")
    print("=" * 150)
    for m in (M_CASING_DESIGN, M_CASING_NEW):
        w = m * ring.G_NET; f_att = 0.1 * w * L_ATTACH
        print(f"  Casing {m:,.0f} kg/m weighs {w/1e3:.1f} kN/m; a tenth over a {L_ATTACH:.0f} m attachment is {f_att/1e6:.3f} MN; trolley line {F_TROLLEY/1e6:.2f} MN at the ring -> {F_TROLLEY/f_att:.1f}, so {math.ceil(F_TROLLEY/f_att)} "
              f"attachments (branches) per line; spreader arm {F_ARM/1e6:.0f} MN -> {F_ARM/f_att:.0f} attachments (Chapter 11: 8 or more and about 80 at 12,000 kg/m).")
    print(f"  Standing cryogenic load: Chapter 14 takes {Q_COIL_PER_M2:.0f} W/m2 of levitation coil facing the cable over {LEV_COIL_WIDTH:.0f} m of coil width, {Q_COIL_PER_M2*LEV_COIL_WIDTH*SITE_SPACING/1e3:.0f} kW per {SITE_SPACING:.0f} m site.")
    print(f"  The sized bearing keeps the optimization's track geometries (two 0.5 m loops per support track, one 0.8 m loop per lift track), so the coil area and the 10 kW")
    print(f"  do not change; only the ampere-turns and the tape fall.  The eddy heat of the LIMs onto the coils and the compressor power are in the deployment rows above.")
    print()

    print("=" * 150)
    print("5. TABLE 20.2 (the book's form; 12,000 kg/m design case, then 6,000 kg/m cases (a), (b), (c))")
    print("=" * 150)
    pk = pick.get("6,000 kg/m, standing load kept at 198 kN/m"); ph = pick.get("6,000 kg/m, standing load halved to 99 kN/m"); p0 = pick.get("design case, 12,000 kg/m, standing load 198 kN/m")
    L0, L016, La, Lb = logs["design 8 MW"], logs["design 16 MW"], logs["a 8 MW"], logs["b 8 MW"]
    Lk, Lk16, Lh, Lh16 = logs["c kept 8 MW"], logs["c kept 16 MW"], logs["c halved 8 MW"], logs["c halved 16 MW"]
    def days(L): return f"{L['days']:.0f} days" if L and L["days"] else "not run"
    def thrust(L): return f"{L['thrust_max_N']/1e3:.2f} kN per LIM; {L['p_eddy_max_kW']:.0f} kW per LIM on day one" if L and L["days"] else "not run"
    def coldK(L): return f"{L['cable_T_max']:.0f} K cable; {L['q_cold_max_kW']:.1f} kW cold side" if L and L["days"] else "not run"
    def att(m): return f"{math.ceil(F_TROLLEY/(0.1*m*ring.G_NET*L_ATTACH))} per trolley line; {F_ARM/(0.1*m*ring.G_NET*L_ATTACH):.0f} per spreader arm"
    dcc = dc.get("kept"); dch = dc.get("halved")
    def rd(d, key, f): return f.format(d[key]) if d else "-"
    rows = [
        ("Casing", "12,000 kg/m", "6,000 kg/m", "6,000 kg/m", "6,000 kg/m", "6,000 kg/m"),
        ("Standing load allowance per support track", "198 kN/m", "198 kN/m", "198 kN/m", "198 kN/m", "99 kN/m"),
        ("Cable hardware", f"{hardware(0.2,0.2):,.0f} kg/m", f"{hardware(0.2,0.2):,.0f} kg/m, unchanged", f"{hardware(0.1,0.1):,.0f} kg/m, halved",
         f"{hw_c.get('kept', float('nan')):,.0f} kg/m", f"{hw_c.get('halved', float('nan')):,.0f} kg/m"),
        ("Structural cable mass", f"{d0['m_struct']:,.0f} kg/m", f"{da['m_struct']:,.0f} kg/m", f"{db['m_struct']:,.0f} kg/m", rd(dcc, "m_struct", "{:,.0f} kg/m"), rd(dch, "m_struct", "{:,.0f} kg/m")),
        ("Cable with hardware", f"{d0['m_cable']:,.0f} kg/m", f"{da['m_cable']:,.0f} kg/m", f"{db['m_cable']:,.0f} kg/m", rd(dcc, "m_cable", "{:,.0f} kg/m"), rd(dch, "m_cable", "{:,.0f} kg/m")),
        ("Cable velocity", f"{d0['v_cable']:,.0f} m/s", f"{da['v_cable']:,.0f} m/s", f"{db['v_cable']:,.0f} m/s", rd(dcc, "v_cable", "{:,.0f} m/s"), rd(dch, "v_cable", "{:,.0f} m/s")),
        ("Tension", f"{d0['tension']/1e9:.1f} GN", f"{da['tension']/1e9:.1f} GN", f"{db['tension']/1e9:.1f} GN", f"{dcc['tension']/1e9:.1f} GN" if dcc else "-", f"{dch['tension']/1e9:.1f} GN" if dch else "-"),
        ("Cross-section", f"{d0['area']:.2f} m2", f"{da['area']:.2f} m2", f"{db['area']:.2f} m2", rd(dcc, "area", "{:.2f} m2"), rd(dch, "area", "{:.2f} m2")),
        ("Deployment energy", f"{d0['E_deploy']/1e18:.2f} EJ", f"{da['E_deploy']/1e18:.2f} EJ", f"{db['E_deploy']/1e18:.2f} EJ", f"{dcc['E_deploy']/1e18:.2f} EJ" if dcc else "-", f"{dch['E_deploy']/1e18:.2f} EJ" if dch else "-"),
        ("Deployment time at 8 MW per site", days(L0), days(La), days(Lb), days(Lk), days(Lh)),
        ("Deployment time at 16 MW", days(L016), "not run", "not run", days(Lk16), days(Lh16)),
        ("LIM peak thrust; first-day heating", thrust(L0), thrust(La), thrust(Lb), thrust(Lk), thrust(Lh)),
        ("Peak cable temperature; cold-side load", coldK(L0), coldK(La), coldK(Lb), coldK(Lk), coldK(Lh)),
        ("Levitation plates", f"four 200 mm, {N_LEV*RHO_FE*W_LEV*0.2:,.0f} kg/m", f"four 200 mm, {N_LEV*RHO_FE*W_LEV*0.2:,.0f} kg/m", f"four 100 mm, {N_LEV*RHO_FE*W_LEV*0.1:,.0f} kg/m",
         f"support {pk['t_sup']*1e3:.0f} mm, lift {pk['t_lift']*1e3:.0f} mm, {pk['m_lev']:,.0f} kg/m" if pk else "-", f"support {ph['t_sup']*1e3:.0f} mm, lift {ph['t_lift']*1e3:.0f} mm, {ph['m_lev']:,.0f} kg/m" if ph else "-"),
        ("LIM reaction plates", f"two 200 mm, {N_LIM*RHO_TIAL*W_LIM*0.2:,.0f} kg/m", f"two 200 mm, {N_LIM*RHO_TIAL*W_LIM*0.2:,.0f} kg/m", f"two 100 mm, {N_LIM*RHO_TIAL*W_LIM*0.1:,.0f} kg/m",
         f"two {lim_t_c*1e3:.0f} mm, {N_LIM*RHO_TIAL*W_LIM*lim_t_c:,.0f} kg/m", f"two {lim_t_c*1e3:.0f} mm, {N_LIM*RHO_TIAL*W_LIM*lim_t_c:,.0f} kg/m"),
        ("Bearing tape, two support and two lift tracks", f"{p0['tape']:.0f} kg/m" if p0 else "-", f"{pk['tape']:.0f} kg/m" if pk else "-", f"{pk['tape']:.0f} kg/m on 200 mm plates; the 100 mm plates fail the saturation rule" if pk else "-",
         f"{pk['tape']:.0f} kg/m" if pk else "-", f"{ph['tape']:.0f} kg/m" if ph else "-"),
        ("Bearing on the casing, tape doubled for formers and cryostats", f"{p0['casing_side']:,.0f} kg/m" if p0 else "-", f"{pk['casing_side']:,.0f} kg/m" if pk else "-", f"{pk['casing_side']:,.0f} kg/m" if pk else "-",
         f"{pk['casing_side']:,.0f} kg/m" if pk else "-", f"{ph['casing_side']:,.0f} kg/m" if ph else "-"),
        ("Anchor attachments", att(M_CASING_DESIGN), att(M_CASING_NEW), att(M_CASING_NEW), att(M_CASING_NEW), att(M_CASING_NEW)),
        ("Standing cryogenic load per site", f"{Q_COIL_PER_M2*LEV_COIL_WIDTH*SITE_SPACING/1e3:.0f} kW", "10 kW", "10 kW", "10 kW", "10 kW"),
    ]
    print(f"  {'Row':<46} | {'12,000 kg/m design case':<40} | {'(a) hardware unchanged':<40} | {'(b) hardware halved':<40} | {'(c1) sized, standing load kept':<40} | (c2) sized, standing load halved")
    for r in rows:
        print(f"  {r[0]:<46} | {r[1]:<40} | {r[2]:<40} | {r[3]:<40} | {r[4]:<40} | {r[5]}")
    print()

if __name__ == "__main__":
    if "--bearing" in sys.argv:
        bearing_sweep()
    elif "--table" in sys.argv:
        table(sys.argv[sys.argv.index("--table") + 1] if len(sys.argv) > sys.argv.index("--table") + 1 else "../../Vol_III/runs_700km_hw10686_2026-10-12")
    else:
        print(__doc__)
