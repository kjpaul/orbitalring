#!/usr/bin/env python3
"""
Mass Driver Launch Simulation — 4-stage LIM, Model 1 (eddy current)

Usage:
    python md_simulation.py [options] [graphs...]

Options:
    --model=1       (only model 1 supported)
    --v_launch=N    Target velocity (m/s)
    --accel=N       Max acceleration (g)
    --mass=N        Spacecraft mass (tonnes)
    --quick         Larger time steps
    --no-graphs     Skip graphs
    --save-csv      Save time-series to CSV

Graph keywords (or 'all'):
    velocity accel thrust power eddy slip temp centrifugal energy
    b_field frequency current voltage skin_depth ring_force occupant_g
    combined stage_detail

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""
import sys, os, math, datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

import md_config as cfg
import md_physics as phys

# ═════════════════════════════════════════════════════════════════════
# DATA COLLECTION
# ═════════════════════════════════════════════════════════════════════

class SimData:
    FIELDS = [
        "time", "velocity", "position", "acceleration",
        "thrust", "P_thrust", "P_eddy", "P_total",
        "slip_ratio", "slip_velocity", "T_plate",
        "F_centrifugal", "F_ring_net", "ring_net_g", "occupant_g",
        "KE", "B_field", "f_supply", "f_slip",
        "current", "voltage", "skin_depth",
        "stage_index", "stage_name", "E_eddy_cumulative",
    ]
    def __init__(self):
        for f in self.FIELDS:
            setattr(self, f, [])
    def record(self, **kw):
        for k, v in kw.items():
            getattr(self, k).append(v)

# ═════════════════════════════════════════════════════════════════════
# CONTROLLER
# ═════════════════════════════════════════════════════════════════════

def run_controller(v_sled, T_plate, stage, n_active, m_total):
    """Determine I and v_slip for this timestep.

    Strategy:
    1. Start at I_target (80% I_c)
    2. Set v_slip to achieve ~target thrust, floor by f_slip_min
    3. Iterate: if thrust > target, reduce v_slip; if < target, increase
    4. Throttle I if V_pfc > V_max
    5. Throttle I if accel > max_g
    6. Thermal protection: reduce v_slip if plate near T_max
    """
    mat = cfg.get_material()
    I = stage.I_target
    F_target = m_total * cfg.MAX_ACCEL_G * cfg.G_ACCEL

    # Slip velocity: need f_slip high enough for EM coupling
    # f_slip_min ~ 1 Hz → v_slip_min = 2 × τ_p × 1
    v_slip_floor = 2.0 * stage.tau_p * 1.0  # Hz

    # Initial estimate based on velocity regime
    if v_sled < 10.0:
        v_slip = max(v_slip_floor, 100.0)
    elif v_sled < 500.0:
        v_slip = max(v_slip_floor, 0.10 * v_sled)
    elif v_sled < 3000.0:
        v_slip = max(v_slip_floor, 0.05 * v_sled)
    else:
        v_slip = max(v_slip_floor, 0.02 * v_sled)

    # Calculate initial thrust
    r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Iterate to find v_slip that gives ~F_target
    for _ in range(8):
        if r["thrust"] < 0.3 * F_target:
            v_slip *= 1.5
        elif r["thrust"] > 1.5 * F_target:
            v_slip *= 0.7
        else:
            break
        v_slip = max(v_slip, v_slip_floor)
        r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Voltage throttle (PFC voltage)
    if r["V_pfc"] > cfg.VOLTS_MAX and r["V_pfc"] > 0:
        ratio = cfg.VOLTS_MAX / r["V_pfc"]
        I *= ratio * 0.95
        I = max(I, 10.0)
        r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Acceleration throttle: thrust ∝ I² → I_new = I × √(F_target/F)
    if r["thrust"] > 0 and m_total > 0:
        a = r["thrust"] / m_total
        a_max = cfg.MAX_ACCEL_G * cfg.G_ACCEL
        if a > a_max * 1.05:
            ratio = math.sqrt(a_max * m_total / r["thrust"])
            I *= ratio * 0.98
            I = max(I, 10.0)
            r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Thermal protection
    T_margin = (mat["T_max"] - T_plate) / (mat["T_max"] - cfg.T_AMBIENT)
    if T_margin < 0.10:
        v_slip *= max(T_margin / 0.10, 0.01)
        v_slip = max(v_slip, 1.0)
        r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Final clamp
    I = min(I, stage.I_peak)
    r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    return {"I": I, "v_slip": v_slip, "r": r}

# ═════════════════════════════════════════════════════════════════════
# MAIN SIMULATION
# ═════════════════════════════════════════════════════════════════════

def run_simulation(v_launch=None, max_accel_g=None, m_spacecraft=None,
                   quick=False):
    if v_launch is not None:
        cfg.V_LAUNCH = v_launch
    if max_accel_g is not None:
        cfg.MAX_ACCEL_G = max_accel_g
    if m_spacecraft is not None:
        cfg.M_SPACECRAFT = m_spacecraft

    dt = cfg.DT_QUICK if quick else cfg.DT
    cfg.print_config()

    d = cfg.calc_derived()
    m_total = d["total_mass"]
    mat = d["mat"]
    n_active = phys.count_active_segments()

    v_sled = cfg.V_INITIAL
    x_sled = 0.0
    T_plate = cfg.T_AMBIENT
    E_eddy_cum = 0.0
    t = 0.0

    data = SimData()
    prev_stage_idx = -1
    handoffs = []

    print(f"\nSimulation: dt={dt}s, model=1 (eddy), target={cfg.V_LAUNCH/1000:.1f} km/s, "
          f"mass={m_total/1000:.0f}t, max_a={cfg.MAX_ACCEL_G}g")
    hdr = (f"{'t(s)':>8} {'v(km/s)':>8} {'stage':>5} {'F(MN)':>7} "
           f"{'a(g)':>6} {'T(K)':>7} {'P(GW)':>7} {'I(A)':>7} "
           f"{'Vpfc(kV)':>9} {'v_slip':>7} {'f_sup':>6}")
    print(hdr)
    print("-" * len(hdr))

    step = 0
    max_steps = int(2e7)
    report_interval = max(1, int(30.0 / dt))

    while v_sled < cfg.V_LAUNCH and step < max_steps:
        # Stage selection
        v_slip_est = max(1.0, 0.02 * max(v_sled, 50.0))
        stage_idx = phys.select_active_stage(v_sled, v_slip_est)
        stage = cfg.LIM_STAGES[stage_idx]

        if stage_idx != prev_stage_idx:
            if prev_stage_idx >= 0:
                pn = cfg.LIM_STAGES[prev_stage_idx].name
                handoffs.append(dict(time=t, frm=pn, to=stage.name,
                    velocity=v_sled, T_plate=T_plate, position=x_sled))
                print(f"  >>> HANDOFF {pn}→{stage.name} at v={v_sled/1000:.2f} km/s, "
                      f"T={T_plate:.1f} K, t={t:.0f}s")
            prev_stage_idx = stage_idx

        # Controller
        ctrl = run_controller(v_sled, T_plate, stage, n_active, m_total)
        r = ctrl["r"]
        I = ctrl["I"]
        v_slip = ctrl["v_slip"]

        thrust = r["thrust"]
        P_eddy = r["P_eddy"]
        B = r["B_peak"]
        f_sup = r["f_supply"]
        f_sl = r["f_slip"]
        delta = r["delta"]
        V_pfc = r["V_pfc"]

        # Dynamics
        accel = thrust / m_total if m_total > 0 else 0.0
        v_sled_new = v_sled + accel * dt
        x_sled += v_sled * dt + 0.5 * accel * dt**2

        # Thermal
        T_plate = phys.plate_temp_rise(P_eddy, dt, T_plate)
        T_plate = phys.radiation_cooling(T_plate, dt)
        E_eddy_cum += P_eddy * dt

        # Forces
        cent = phys.centrifugal_force(v_sled, m_total)
        occ_g = phys.occupant_g(thrust, m_total)
        P_thrust = thrust * v_sled
        P_total = P_thrust + P_eddy
        s_ratio = v_slip / (v_sled + v_slip) if (v_sled + v_slip) > 0 else 1.0

        # Record
        data.record(
            time=t, velocity=v_sled, position=x_sled,
            acceleration=accel, thrust=thrust,
            P_thrust=P_thrust, P_eddy=P_eddy, P_total=P_total,
            slip_ratio=s_ratio, slip_velocity=v_slip, T_plate=T_plate,
            F_centrifugal=cent["F_centrifugal"], F_ring_net=cent["F_net"],
            ring_net_g=cent["net_g"], occupant_g=occ_g,
            KE=0.5 * m_total * v_sled**2,
            B_field=B, f_supply=f_sup, f_slip=f_sl,
            current=I, voltage=V_pfc, skin_depth=delta,
            stage_index=stage_idx, stage_name=stage.name,
            E_eddy_cumulative=E_eddy_cum,
        )

        v_sled = v_sled_new
        t += dt
        step += 1

        if step % report_interval == 0:
            print(f"{t:>8.0f} {v_sled/1000:>8.2f} {stage.name:>5} "
                  f"{thrust/1e6:>7.2f} {accel/cfg.G_ACCEL:>6.2f} "
                  f"{T_plate:>7.1f} {P_total/1e9:>7.2f} {I:>7.0f} "
                  f"{V_pfc/1000:>9.1f} {v_slip:>7.1f} {f_sup:>6.1f}")

    # Final report
    sep = "=" * 76
    print(f"\n{sep}")
    print("SIMULATION COMPLETE")
    print(sep)
    print(f"  Final velocity       {v_sled/1000:>10.3f} km/s")
    print(f"  Total time           {t:>10.0f} s  ({t/60:.1f} min, {t/3600:.2f} hr)")
    print(f"  Distance             {x_sled/1000:>10.0f} km  ({x_sled/cfg.L_RING*100:.2f}% of ring)")
    print(f"  Plate temp (final)   {T_plate:>10.1f} K  (limit {mat['T_max']} K)")
    print(f"  Total eddy heat      {E_eddy_cum:>12.3e} J  ({E_eddy_cum/3.6e9:.1f} MWh)")
    print(f"  Kinetic energy       {0.5*m_total*v_sled**2:>12.3e} J  ({0.5*m_total*v_sled**2/3.6e12:.1f} GWh)")
    print(f"  Peak thrust          {max(data.thrust)/1e6:>10.2f} MN")
    print(f"  Peak power           {max(data.P_total)/1e9:>10.2f} GW")
    print(f"  Peak plate temp      {max(data.T_plate):>10.1f} K")
    print(f"  Peak occupant g      {max(data.occupant_g):>10.3f} g")
    if handoffs:
        print(f"\n  HANDOFFS:")
        for h in handoffs:
            print(f"    {h['frm']}→{h['to']} at t={h['time']:.0f}s, "
                  f"v={h['velocity']/1000:.2f} km/s, T={h['T_plate']:.1f} K")
    print(sep)
    return data

# ═════════════════════════════════════════════════════════════════════
# PLOTTING
# ═════════════════════════════════════════════════════════════════════

def _t_min(data):
    return [t/60 for t in data.time]

def _handoff_lines(ax, data):
    prev = data.stage_index[0] if data.stage_index else 0
    for i, s in enumerate(data.stage_index):
        if s != prev:
            ax.axvline(data.time[i]/60, color='red', ls='--', alpha=0.4, lw=0.8)
            prev = s

def _plot1(data, ydata, ylabel, title, fname, color='#1976D2', yscale='linear'):
    if not HAS_MPL: return
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(_t_min(data), ydata, color=color, lw=1.2)
    _handoff_lines(ax, data)
    ax.set_xlabel("Time (min)"); ax.set_ylabel(ylabel)
    ax.set_title(title); ax.set_yscale(yscale); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(cfg.GRAPH_DIR, fname), dpi=200)
    plt.close(fig); print(f"  Saved {fname}")

def plot_combined(data):
    if not HAS_MPL: return
    tm = _t_min(data)
    fig, axes = plt.subplots(3, 3, figsize=(18, 14))
    fig.suptitle("Mass Driver Launch — 4-Stage LIM (Model 1)", fontsize=15, y=0.98)
    plots = [
        (axes[0,0], [v/1000 for v in data.velocity], "Velocity (km/s)", '#1976D2'),
        (axes[0,1], [a/cfg.G_ACCEL for a in data.acceleration], "Accel (g)", '#E65100'),
        (axes[0,2], [F/1e6 for F in data.thrust], "Thrust (MN)", '#2E7D32'),
        (axes[1,0], [P/1e9 for P in data.P_total], "Total Power (GW)", '#7B1FA2'),
        (axes[1,1], [P/1e6 for P in data.P_eddy], "Eddy Losses (MW)", '#C62828'),
        (axes[1,2], data.T_plate, "Plate Temp (K)", '#FF6F00'),
        (axes[2,0], data.current, "Current (A)", '#1565C0'),
        (axes[2,1], [V/1000 for V in data.voltage], "V_pfc (kV)", '#B71C1C'),
        (axes[2,2], data.slip_velocity, "v_slip (m/s)", '#00838F'),
    ]
    for ax, yd, yl, c in plots:
        ax.plot(tm, yd, color=c, lw=1.0)
        _handoff_lines(ax, data)
        ax.set_ylabel(yl, fontsize=9); ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)
    for ax in axes[2]:
        ax.set_xlabel("Time (min)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(cfg.GRAPH_DIR, "md_combined.png"), dpi=200)
    plt.close(fig); print("  Saved md_combined.png")

def plot_stage_detail(data):
    if not HAS_MPL: return
    snames = ["S1", "S2", "S3", "S4"]
    scolors = ['#2196F3', '#FF9800', '#4CAF50', '#9C27B0']
    sdata = {}
    for nm in snames:
        idx = [i for i, s in enumerate(data.stage_name) if s == nm]
        if idx: sdata[nm] = idx
    if not sdata: return

    ncol = len(sdata)
    fig, axes = plt.subplots(4, ncol, figsize=(6*ncol, 16))
    if ncol == 1: axes = axes.reshape(-1, 1)
    fig.suptitle("Mass Driver — Stage Detail", fontsize=14, y=0.99)

    for col, (nm, idx) in enumerate(sdata.items()):
        tm = [data.time[i]/60 for i in idx]
        si = snames.index(nm)
        c = scolors[si]
        tp = cfg.LIM_STAGES[si].tau_p

        axes[0,col].plot(tm, [data.velocity[i]/1000 for i in idx], color=c, lw=1.2)
        axes[0,col].set_ylabel("Velocity (km/s)")
        axes[0,col].set_title(f"Stage {nm} (τ_p = {tp} m)")

        ax1 = axes[1,col]; ax2 = ax1.twinx()
        ax1.plot(tm, [data.current[i] for i in idx], color='blue', lw=1)
        ax2.plot(tm, [data.voltage[i]/1000 for i in idx], color='red', lw=1)
        ax1.set_ylabel("Current (A)", color='blue')
        ax2.set_ylabel("V_pfc (kV)", color='red')

        ax3 = axes[2,col]; ax4 = ax3.twinx()
        ax3.plot(tm, [data.thrust[i]/1e6 for i in idx], color='green', lw=1)
        ax4.plot(tm, [data.P_total[i]/1e9 for i in idx], color='purple', lw=1)
        ax3.set_ylabel("Thrust (MN)", color='green')
        ax4.set_ylabel("Power (GW)", color='purple')

        ax5 = axes[3,col]; ax6 = ax5.twinx()
        ax5.plot(tm, [data.T_plate[i] for i in idx], color='orange', lw=1)
        ax6.plot(tm, [data.slip_velocity[i] for i in idx], color='cyan', lw=1)
        ax5.set_ylabel("T_plate (K)", color='orange')
        ax6.set_ylabel("v_slip (m/s)", color='cyan')
        ax5.set_xlabel("Time (min)")

    for row in axes:
        for ax in row:
            ax.grid(True, alpha=0.2); ax.tick_params(labelsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(cfg.GRAPH_DIR, "md_stage_detail.png"), dpi=200)
    plt.close(fig); print("  Saved md_stage_detail.png")

GRAPHS = {
    "velocity": lambda d: _plot1(d, [v/1000 for v in d.velocity], "km/s", "Sled Velocity", "md_velocity.png"),
    "accel": lambda d: _plot1(d, [a/cfg.G_ACCEL for a in d.acceleration], "g", "Acceleration", "md_accel.png", '#E65100'),
    "thrust": lambda d: _plot1(d, [F/1e6 for F in d.thrust], "MN", "Thrust", "md_thrust.png", '#2E7D32'),
    "power": lambda d: _plot1(d, [P/1e9 for P in d.P_total], "GW", "Total Power", "md_power.png", '#7B1FA2'),
    "eddy": lambda d: _plot1(d, [P/1e6 for P in d.P_eddy], "MW", "Eddy Losses", "md_eddy.png", '#C62828'),
    "slip": lambda d: _plot1(d, d.slip_velocity, "m/s", "Slip Velocity", "md_slip.png", '#00838F'),
    "temp": lambda d: _plot1(d, d.T_plate, "K", "Plate Temperature", "md_temp.png", '#FF6F00'),
    "centrifugal": lambda d: _plot1(d, [F/1e6 for F in d.F_centrifugal], "MN", "Centrifugal Force", "md_centrifugal.png", '#795548'),
    "ring_force": lambda d: _plot1(d, [F/1e6 for F in d.F_ring_net], "MN", "Net Ring Force", "md_ring_force.png"),
    "occupant_g": lambda d: _plot1(d, d.occupant_g, "g", "Occupant G-Force", "md_occupant_g.png", '#E65100'),
    "energy": lambda d: _plot1(d, [E/1e12 for E in d.KE], "TJ", "Kinetic Energy", "md_energy.png", '#0D47A1'),
    "b_field": lambda d: _plot1(d, d.B_field, "T", "B Field at Plate", "md_b_field.png", '#4527A0'),
    "frequency": lambda d: _plot1(d, d.f_supply, "Hz", "Supply Frequency", "md_frequency.png", '#00695C'),
    "current": lambda d: _plot1(d, d.current, "A", "Phase Current", "md_current.png", '#1565C0'),
    "voltage": lambda d: _plot1(d, [V/1000 for V in d.voltage], "kV", "V_pfc", "md_voltage.png", '#B71C1C'),
    "skin_depth": lambda d: _plot1(d, [sd*1000 for sd in d.skin_depth], "mm", "Skin Depth", "md_skin_depth.png", '#4E342E'),
    "combined": lambda d: plot_combined(d),
    "stage_detail": lambda d: plot_stage_detail(d),
}

def save_csv(data, fname="md_launch_data.csv"):
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    path = os.path.join(cfg.GRAPH_DIR, fname)
    with open(path, 'w') as f:
        f.write(",".join(SimData.FIELDS) + "\n")
        for i in range(len(data.time)):
            row = [getattr(data, fld)[i] for fld in SimData.FIELDS]
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  Saved {path}")

# ═════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════

def main():
    v_launch = max_accel = mass = None
    quick = no_graphs = do_csv = False
    graph_keys = []

    for arg in sys.argv[1:]:
        if arg.startswith("--v_launch="): v_launch = float(arg.split("=")[1])
        elif arg.startswith("--accel="): max_accel = float(arg.split("=")[1])
        elif arg.startswith("--mass="): mass = float(arg.split("=")[1]) * 1000
        elif arg.startswith("--model="): pass  # only model 1
        elif arg == "--quick": quick = True
        elif arg == "--no-graphs": no_graphs = True
        elif arg == "--save-csv": do_csv = True
        elif arg == "--help": print(__doc__); sys.exit(0)
        else: graph_keys.append(arg)

    if not graph_keys and not no_graphs:
        graph_keys = ["combined", "stage_detail"]

    data = run_simulation(v_launch=v_launch, max_accel_g=max_accel,
                          m_spacecraft=mass, quick=quick)

    if not no_graphs:
        print("\nGenerating graphs...")
        keys = list(GRAPHS.keys()) if "all" in graph_keys else graph_keys
        for k in keys:
            if k in GRAPHS:
                try: GRAPHS[k](data)
                except Exception as e: print(f"  [error] {k}: {e}")

    if do_csv:
        save_csv(data)

if __name__ == "__main__":
    main()
