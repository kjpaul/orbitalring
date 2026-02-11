#!/usr/bin/env python3
"""
Mass Driver Launch Simulation — 4-stage LIM, Model 1 (eddy current)

Usage:
    python md_simulation.py [options] [graphs...]

Options:
    --v_launch=N    Target velocity (m/s)
    --accel=N       Max acceleration (g)
    --mass=N        Spacecraft mass (tonnes)
    --quick         Larger time steps (1s vs 0.1s)
    --no-graphs     Skip graph generation
    --save-csv      Save time-series to CSV

Graph keywords (or 'all'):
    velocity accel thrust power eddy slip temp occupant_g
    b_field frequency current voltage skin_depth ring_force
    combined stage_detail

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""
import sys, os, math

# Ensure UTF-8 output on Windows (γ-TiAl, τ_p, etc.)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

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
        "F_centrifugal", "F_ring_net", "ring_net_g",
        "g_thrust", "g_radial", "g_total",
        "KE", "B_field", "f_supply", "f_slip",
        "current", "voltage", "V_coil", "skin_depth",
        "stage_index", "stage_name", "E_eddy_cumulative",
        "P_hts_ac",
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
    """Determine I and v_slip for this timestep."""
    mat = cfg.get_material()
    I = stage.I_target
    F_target = m_total * cfg.MAX_ACCEL_G * cfg.G_ACCEL

    v_slip_floor = 2.0 * stage.tau_p * 1.0

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

    # Voltage throttle — supply voltage
    if r["V_supply"] > cfg.VOLTS_MAX and r["V_supply"] > 0:
        ratio = cfg.VOLTS_MAX / r["V_supply"]
        I *= ratio * 0.95
        I = max(I, 10.0)
        r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Voltage throttle — coil insulation (V_coil = ωLI must be < 100 kV)
    if r["V_coil"] > cfg.VOLTS_MAX and r["V_coil"] > 0:
        # V_coil ∝ I (at fixed frequency), so scale I directly
        ratio = cfg.VOLTS_MAX / r["V_coil"]
        I *= ratio * 0.95
        I = max(I, 10.0)
        r = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate, n_active)

    # Acceleration throttle
    a_max = cfg.MAX_ACCEL_G * cfg.G_ACCEL
    for _ in range(8):
        if r["thrust"] <= 0 or m_total <= 0:
            break
        a = r["thrust"] / m_total
        if a <= a_max * 1.02:
            break
        # Strategy: reduce both I and v_slip
        overshoot = a / a_max
        if overshoot > 2.0:
            # Large overshoot — cut v_slip aggressively
            v_slip *= 0.5
        elif overshoot > 1.2:
            v_slip *= 0.8
        else:
            v_slip *= 0.95
        # Also reduce I if above minimum
        I *= min(0.95, math.sqrt(a_max * m_total / r["thrust"]))
        I = max(I, 10.0)
        v_slip = max(v_slip, v_slip_floor)
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
    hdr = (f"{'t(s)':>8} {'v(km/s)':>8} {'stg':>4} {'F(MN)':>7} "
           f"{'a(g)':>7} {'g_rad':>7} {'g_tot':>7} {'T(K)':>7} "
           f"{'P(GW)':>6} {'I(A)':>6} {'Vsup':>5} {'Vcoil':>5} {'vslip':>6}")
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
        V_supply = r["V_supply"]
        V_coil = r["V_coil"]

        # Dynamics
        accel = thrust / m_total if m_total > 0 else 0.0
        v_sled_new = v_sled + accel * dt
        x_sled += v_sled * dt + 0.5 * accel * dt**2

        # Thermal
        T_plate = phys.plate_temp_rise(P_eddy, dt, T_plate)
        T_plate = phys.radiation_cooling(T_plate, dt)
        E_eddy_cum += P_eddy * dt

        # Forces — occupant feels thrust + centrifugal - gravity
        occ = phys.occupant_forces(thrust, v_sled, m_total)
        cent = phys.centrifugal_force(v_sled, m_total)

        P_thrust = thrust * v_sled
        P_total = P_thrust + P_eddy
        s_ratio = v_slip / (v_sled + v_slip) if (v_sled + v_slip) > 0 else 1.0

        # HTS AC losses (from calc_thrust)
        P_hts = r.get("P_hts_coil", 0.0)

        # Record
        data.record(
            time=t, velocity=v_sled, position=x_sled,
            acceleration=accel, thrust=thrust,
            P_thrust=P_thrust, P_eddy=P_eddy, P_total=P_total,
            slip_ratio=s_ratio, slip_velocity=v_slip, T_plate=T_plate,
            F_centrifugal=cent["F_centrifugal"], F_ring_net=cent["F_net"],
            ring_net_g=cent["net_g"],
            g_thrust=occ["g_thrust"], g_radial=occ["g_radial"],
            g_total=occ["g_total"],
            KE=0.5 * m_total * v_sled**2,
            B_field=B, f_supply=f_sup, f_slip=f_sl,
            current=I, voltage=V_supply, V_coil=V_coil, skin_depth=delta,
            stage_index=stage_idx, stage_name=stage.name,
            E_eddy_cumulative=E_eddy_cum, P_hts_ac=P_hts,
        )

        v_sled = v_sled_new
        t += dt
        step += 1

        if step % report_interval == 0:
            print(f"{t:>8.0f} {v_sled/1000:>8.2f} {stage.name:>4} "
                  f"{thrust/1e6:>7.2f} {occ['g_thrust']:>7.3f} "
                  f"{occ['g_radial']:>+7.3f} {occ['g_total']:>7.3f} "
                  f"{T_plate:>7.1f} {P_total/1e9:>6.1f} {I:>6.0f} "
                  f"{V_supply/1000:>5.1f} {V_coil/1000:>5.1f} {v_slip:>6.1f}")

    # Final report
    sep = "=" * 76
    occ_final = phys.occupant_forces(0, v_sled, m_total)  # coast after launch
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
    print(f"\n  OCCUPANT G-FORCE:")
    print(f"    At start:   g_thrust={data.g_thrust[0]:>+6.3f}, g_radial={data.g_radial[0]:>+6.3f}, g_total={data.g_total[0]:.3f}")
    if data.g_total:
        i_peak = max(range(len(data.g_total)), key=lambda i: data.g_total[i])
        print(f"    Peak total: g_thrust={data.g_thrust[i_peak]:>+6.3f}, g_radial={data.g_radial[i_peak]:>+6.3f}, g_total={data.g_total[i_peak]:.3f}")
    print(f"    At launch:  g_thrust={data.g_thrust[-1]:>+6.3f}, g_radial={data.g_radial[-1]:>+6.3f}, g_total={data.g_total[-1]:.3f}")
    print(f"    After release (coast): g_radial={occ_final['g_radial']:>+6.3f}")

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

def _param_subtitle():
    """Build parameter subtitle string from current config."""
    s0 = cfg.LIM_STAGES[0]
    tp_min = min(s.tau_p for s in cfg.LIM_STAGES)
    tp_max = max(s.tau_p for s in cfg.LIM_STAGES)
    tp_str = f"{tp_min:.0f}\u2013{tp_max:.0f}" if tp_min != tp_max else f"{tp_min:.0f}"
    return (f"(\u03c4\u209a: {tp_str} m, N: {s0.n_turns}, "
            f"{s0.hts_layers}\u00d7{s0.hts_width_mm:.0f} mm HTS, "
            f"m: {cfg.M_SPACECRAFT/1000:,.0f} t, "
            f"v: {cfg.V_LAUNCH/1000:.1f} km/s)")

def _plot1(data, ydata, ylabel, title, fname, color='#1976D2'):
    if not HAS_MPL: return
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(_t_min(data), ydata, color=color, lw=1.2)
    _handoff_lines(ax, data)
    ax.set_xlabel("Time (min)"); ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\n{_param_subtitle()}", fontsize=12)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(cfg.GRAPH_DIR, fname), dpi=200)
    plt.close(fig); print(f"  Saved {fname}")

def plot_combined(data):
    if not HAS_MPL: return
    tm = _t_min(data)
    fig, axes = plt.subplots(3, 3, figsize=(18, 14))
    fig.suptitle(f"Mass Driver Launch \u2014 4-Stage LIM (Model 1 Eddy Current)\n{_param_subtitle()}",
                 fontsize=13, y=0.99)

    # Disable scientific notation offset on all axes
    for row in axes:
        for ax in row:
            ax.ticklabel_format(useOffset=False)

    # Row 0: velocity, thrust accel, thrust force
    axes[0,0].plot(tm, [v/1000 for v in data.velocity], color='#1976D2', lw=1.0)
    _handoff_lines(axes[0,0], data); axes[0,0].set_ylabel("Velocity (km/s)")

    axes[0,1].plot(tm, data.g_thrust, color='#E65100', lw=1.0)
    _handoff_lines(axes[0,1], data); axes[0,1].set_ylabel("Thrust (g)")
    axes[0,1].set_ylim(0, max(data.g_thrust) * 1.3)

    axes[0,2].plot(tm, [F/1e6 for F in data.thrust], color='#2E7D32', lw=1.0)
    _handoff_lines(axes[0,2], data); axes[0,2].set_ylabel("Thrust (MN)")
    axes[0,2].set_ylim(0, max(data.thrust)/1e6 * 1.3)

    # Row 1: KE delivery rate (P_thrust), LIM eddy losses, plate temp
    axes[1,0].plot(tm, [P/1e9 for P in data.P_thrust], color='#7B1FA2', lw=1.0, label='KE rate')
    axes[1,0].plot(tm, [P/1e9 for P in data.P_eddy], color='#C62828', lw=1.0, label='Eddy loss')
    _handoff_lines(axes[1,0], data)
    axes[1,0].set_ylabel("Power (GW)"); axes[1,0].legend(fontsize=8)

    axes[1,1].plot(tm, [P/1e6 for P in data.P_eddy], color='#C62828', lw=1.0)
    _handoff_lines(axes[1,1], data); axes[1,1].set_ylabel("Eddy Losses (MW)")

    axes[1,2].plot(tm, data.T_plate, color='#FF6F00', lw=1.0)
    _handoff_lines(axes[1,2], data); axes[1,2].set_ylabel("Plate Temp (K)")
    axes[1,2].set_ylim(200, max(data.T_plate) * 1.1)

    # Row 2: radial g, total g, V_coil
    axes[2,0].plot(tm, data.g_radial, color='#00695C', lw=1.0)
    axes[2,0].axhline(0, color='gray', ls='-', alpha=0.3, lw=0.5)
    _handoff_lines(axes[2,0], data); axes[2,0].set_ylabel("Radial g (+=outward)")

    axes[2,1].plot(tm, data.g_total, color='#D32F2F', lw=1.0)
    _handoff_lines(axes[2,1], data); axes[2,1].set_ylabel("Total g-load")

    axes[2,2].plot(tm, [V/1000 for V in data.V_coil], color='#B71C1C', lw=1.0)
    axes[2,2].axhline(100, color='red', ls='--', alpha=0.5, lw=1, label='100 kV limit')
    _handoff_lines(axes[2,2], data)
    axes[2,2].set_ylabel("V_coil (kV)"); axes[2,2].legend(fontsize=8)

    for ax in axes.flat:
        ax.grid(True, alpha=0.3); ax.tick_params(labelsize=8)
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
    fig.suptitle(f"Mass Driver \u2014 Stage Detail\n{_param_subtitle()}", fontsize=12, y=0.99)

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
        ax2.set_ylabel("V_supply (kV)", color='red')

        ax3 = axes[2,col]; ax4 = ax3.twinx()
        ax3.plot(tm, [data.thrust[i]/1e6 for i in idx], color='green', lw=1)
        ax4.plot(tm, [data.P_total[i]/1e9 for i in idx], color='purple', lw=1)
        ax3.set_ylabel("Thrust (MN)", color='green')
        ax4.set_ylabel("Power (GW)", color='purple')

        ax5 = axes[3,col]; ax6 = ax5.twinx()
        ax5.plot(tm, [data.T_plate[i] for i in idx], color='orange', lw=1)
        ax6.plot(tm, [data.g_total[i] for i in idx], color='#D32F2F', lw=1)
        ax5.set_ylabel("T_plate (K)", color='orange')
        ax6.set_ylabel("Total g-load", color='#D32F2F')
        ax5.set_xlabel("Time (min)")

    for row in axes:
        for ax in row:
            ax.grid(True, alpha=0.2); ax.tick_params(labelsize=8)
            ax.ticklabel_format(useOffset=False)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(cfg.GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(cfg.GRAPH_DIR, "md_stage_detail.png"), dpi=200)
    plt.close(fig); print("  Saved md_stage_detail.png")

GRAPHS = {
    "velocity": lambda d: _plot1(d, [v/1000 for v in d.velocity], "km/s", "Sled Velocity", "md_velocity.png"),
    "accel": lambda d: _plot1(d, d.g_thrust, "g", "Thrust Acceleration", "md_accel.png", '#E65100'),
    "thrust": lambda d: _plot1(d, [F/1e6 for F in d.thrust], "MN", "Thrust", "md_thrust.png", '#2E7D32'),
    "power": lambda d: _plot1(d, [P/1e9 for P in d.P_total], "GW", "Total Power", "md_power.png", '#7B1FA2'),
    "eddy": lambda d: _plot1(d, [P/1e6 for P in d.P_eddy], "MW", "Eddy Losses", "md_eddy.png", '#C62828'),
    "slip": lambda d: _plot1(d, d.slip_velocity, "m/s", "Slip Velocity", "md_slip.png", '#00838F'),
    "temp": lambda d: _plot1(d, d.T_plate, "K", "Plate Temperature", "md_temp.png", '#FF6F00'),
    "ring_force": lambda d: _plot1(d, [F/1e6 for F in d.F_ring_net], "MN", "Net Ring Force", "md_ring_force.png"),
    "occupant_g": lambda d: _plot1(d, d.g_total, "g", "Total Occupant G-Load", "md_occupant_g.png", '#D32F2F'),
    "energy": lambda d: _plot1(d, [E/1e12 for E in d.KE], "TJ", "Kinetic Energy", "md_energy.png", '#0D47A1'),
    "b_field": lambda d: _plot1(d, d.B_field, "T", "B Field at Plate", "md_b_field.png", '#4527A0'),
    "frequency": lambda d: _plot1(d, d.f_supply, "Hz", "Supply Frequency", "md_frequency.png", '#00695C'),
    "current": lambda d: _plot1(d, d.current, "A", "Phase Current", "md_current.png", '#1565C0'),
    "voltage": lambda d: _plot1(d, [V/1000 for V in d.voltage], "kV", "V_supply", "md_voltage.png", '#B71C1C'),
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
