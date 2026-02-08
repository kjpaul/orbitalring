#!/usr/bin/env python3
"""
Mass Driver Launch Simulation - Main simulation loop and plotting

Simulates electromagnetic launch of a spacecraft from the orbital ring
mass driver. Three sequential LIM stages accelerate a sled from rest
to the target launch velocity (e.g. 15 km/s for Mars transfer).

Usage:
    python md_simulation.py [options] [graphs...]

Options:
    --model=N       Thrust model: 1=eddy (default), 2=goodness, 3=slip×pressure
    --v_launch=N    Target launch velocity (m/s), default from config
    --accel=N       Max acceleration (g), default from config
    --mass=N        Spacecraft mass (tonnes), default from config
    --quick         Larger time steps for fast testing
    --no-graphs     Skip all graph generation
    --save-csv      Save time-series data to CSV

Graph keywords (or 'all'):
    velocity  accel  thrust  power  eddy  slip  temp  centrifugal
    energy  b_field  frequency  current  voltage  skin_depth  goodness
    ring_force  occupant_g  combined

Reference: "Orbital Ring Engineering" by Paul G de Jong
"""

import sys
import os
import math
import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False

import md_config as cfg
import md_physics as phys


GRAPH_DIR = "graphs"


# =============================================================================
# DATA COLLECTION
# =============================================================================

class SimData:
    """Collects time-series data during the simulation."""

    def __init__(self):
        self.time = []
        self.velocity = []
        self.position = []
        self.acceleration = []
        self.thrust = []
        self.P_thrust = []          # Thrust power (F × v)
        self.P_eddy = []            # Eddy current loss in plate
        self.P_total = []           # Total electrical power input
        self.slip_ratio = []
        self.slip_velocity = []
        self.T_plate = []
        self.F_centrifugal = []
        self.F_ring_net = []        # Net force on ring from sled
        self.ring_net_g = []
        self.occupant_g = []
        self.KE = []
        self.B_field = []
        self.f_supply = []
        self.f_slip = []
        self.current = []
        self.voltage = []
        self.skin_depth = []
        self.goodness = []
        self.stage_index = []
        self.stage_name = []
        self.E_eddy_cumulative = []  # Cumulative heat deposited in plate (J)

    def record(self, **kwargs):
        for key, val in kwargs.items():
            getattr(self, key).append(val)


# =============================================================================
# CONTROLLER
# =============================================================================

def run_controller(v_sled, T_plate, stage, n_active, m_total, model, dt):
    """Main control logic for one timestep.

    Determines slip velocity and current for the active stage, enforcing:
    - Voltage limit (100 kV)
    - Current limit (HTS tape)
    - Acceleration limit (0.5g human-rated)
    - Thermal limit (plate T_max)

    Strategy:
    - Start with target current (80% of I_c)
    - Calculate slip velocity to give desired thrust
    - Throttle current if voltage or accel limits exceeded
    - If plate nearing T_max, reduce slip to cut eddy losses

    Returns:
        dict with: I_peak, v_slip, thrust_result (from calc_thrust),
                   stage, n_active
    """
    mat = cfg.get_material()

    # Start with target current
    I = stage.I_target

    # Target thrust for desired acceleration
    F_target = m_total * cfg.MAX_ACCEL_G * cfg.G_ACCEL

    # Slip velocity strategy:
    # We want enough thrust for ~0.5g but with minimal eddy heating.
    # The slip must produce enough f_slip for eddy current coupling.
    # f_slip = v_slip / (2 × τ_p)
    # For useful coupling in γ-TiAl (high ρ), we need f_slip > ~1 Hz.
    # Minimum v_slip = 2 × τ_p × f_slip_min

    f_slip_min = 1.0  # Hz — minimum for useful coupling
    v_slip_floor = 2.0 * stage.tau_p * f_slip_min

    if v_sled < 10.0:
        # Near-stationary: use generous slip for startup
        v_slip = max(v_slip_floor, 50.0)
    elif v_sled < 500.0:
        v_slip = max(v_slip_floor, 0.05 * v_sled)  # 5% slip
    elif v_sled < 2000.0:
        v_slip = max(v_slip_floor, 0.02 * v_sled)  # 2% slip
    else:
        v_slip = max(v_slip_floor, 0.005 * v_sled)  # 0.5% slip at high speed

    # Calculate thrust with initial parameters
    result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                              n_active, model)

    # --- Thrust feedback: adjust slip if thrust too low ---
    # If thrust is less than 50% of target, increase slip to get more coupling
    iterations = 0
    while result["thrust"] < 0.5 * F_target and iterations < 5:
        v_slip *= 1.5  # Increase slip by 50%
        result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                                  n_active, model)
        iterations += 1

    # --- Voltage throttle (using PFC voltage, scales with v_slip) ---
    if result["V_pfc"] > cfg.VOLTS_MAX and result["V_pfc"] > 0:
        ratio = cfg.VOLTS_MAX / result["V_pfc"]
        I = I * ratio * 0.95
        I = max(I, 10.0)  # Floor to avoid zero current
        result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                                  n_active, model)

    # --- Acceleration throttle ---
    if result["thrust"] > 0 and m_total > 0:
        a = result["thrust"] / m_total
        a_max = cfg.MAX_ACCEL_G * cfg.G_ACCEL
        if a > a_max * 1.05:  # 5% tolerance
            # Thrust ∝ I², so I_new = I × sqrt(F_target / F_actual)
            ratio = math.sqrt(a_max * m_total / result["thrust"])
            I = I * ratio * 0.98
            I = max(I, 10.0)
            result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                                      n_active, model)

    # --- Thermal protection ---
    # If plate is within 10% of T_max, start reducing slip
    T_margin = (mat["T_max"] - T_plate) / (mat["T_max"] - cfg.T_AMBIENT)
    if T_margin < 0.10:
        # Reduce slip velocity proportionally
        v_slip = v_slip * (T_margin / 0.10)
        v_slip = max(v_slip, 0.1)
        result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                                  n_active, model)

    # --- Ensure current doesn't exceed HTS limit ---
    I = min(I, stage.I_peak)

    # Final thrust calculation with all throttling applied
    result = phys.calc_thrust(stage, I, v_slip, v_sled, T_plate,
                              n_active, model)

    return {
        "I_peak": I,
        "v_slip": v_slip,
        "result": result,
    }


# =============================================================================
# MAIN SIMULATION LOOP
# =============================================================================

def run_simulation(model=None, v_launch=None, max_accel_g=None,
                   m_spacecraft=None, quick=False):
    """Run the full mass driver launch simulation.

    Args:
        model: Thrust model override (1, 2, or 3)
        v_launch: Target launch velocity override (m/s)
        max_accel_g: Max acceleration override (g)
        m_spacecraft: Spacecraft mass override (kg)
        quick: Use larger time steps

    Returns:
        SimData object with all time-series data
    """
    # Apply overrides
    if model is not None:
        cfg.THRUST_MODEL = model
    if v_launch is not None:
        cfg.V_LAUNCH = v_launch
    if max_accel_g is not None:
        cfg.MAX_ACCEL_G = max_accel_g
    if m_spacecraft is not None:
        cfg.M_SPACECRAFT = m_spacecraft

    dt = cfg.DT_QUICK if quick else cfg.DT
    _model = cfg.THRUST_MODEL

    # Print config
    cfg.print_config()

    # Derived parameters
    d = cfg.calc_derived()
    m_total = d["total_mass"]
    mat = d["mat"]

    # Initial state
    v_sled = cfg.V_INITIAL
    x_sled = 0.0
    T_plate = cfg.T_AMBIENT
    E_eddy_cum = 0.0
    t = 0.0

    data = SimData()

    # Stage tracking
    prev_stage_idx = -1
    handoff_data = []  # Record handoff events

    print(f"\nStarting simulation (model={_model}, dt={dt}s)...")
    print(f"Target: {cfg.V_LAUNCH/1000:.1f} km/s, Mass: {m_total/1000:.0f} t, "
          f"Max accel: {cfg.MAX_ACCEL_G}g")
    print("-" * 72)

    step = 0
    max_steps = int(1e8)  # Safety limit

    while v_sled < cfg.V_LAUNCH and step < max_steps:
        # --- Stage selection ---
        # Estimate slip velocity for stage selection
        v_slip_est = max(1.0, 0.005 * max(v_sled, 100.0))
        stage_idx = phys.select_active_stage(v_sled, v_slip_est)
        stage = cfg.LIM_STAGES[stage_idx]

        # Count active segments
        n_active = phys.count_active_segments(stage_idx)

        # Report handoffs
        if stage_idx != prev_stage_idx:
            if prev_stage_idx >= 0:
                prev_name = cfg.LIM_STAGES[prev_stage_idx].name
                handoff_data.append({
                    "time": t,
                    "from": prev_name,
                    "to": stage.name,
                    "velocity": v_sled,
                    "T_plate": T_plate,
                    "position": x_sled,
                })
                print(f"  HANDOFF at t={t:.1f}s: {prev_name} → {stage.name} "
                      f"@ v={v_sled:.1f} m/s ({v_sled/1000:.2f} km/s), "
                      f"T_plate={T_plate:.1f} K, x={x_sled/1000:.1f} km")
            else:
                print(f"  START on {stage.name} stage")
            prev_stage_idx = stage_idx

        # --- Controller ---
        ctrl = run_controller(v_sled, T_plate, stage, n_active,
                              m_total, _model, dt)

        result = ctrl["result"]
        I = ctrl["I_peak"]
        v_slip = ctrl["v_slip"]

        thrust = result["thrust"]
        P_eddy = result["P_eddy"]
        B = result["B_peak"]
        f_sup = result["f_supply"]
        f_sl = result["f_slip"]
        delta = result["delta"]
        V_back = result["V_back"]
        V_pfc = result.get("V_pfc", V_back)
        G = result["goodness"]

        # --- Dynamics ---
        accel = thrust / m_total if m_total > 0 else 0.0
        v_sled_new = v_sled + accel * dt
        x_sled += v_sled * dt + 0.5 * accel * dt**2

        # --- Thermal ---
        T_plate = phys.plate_temp_rise(P_eddy, dt, T_plate)
        T_plate = phys.radiation_cooling(T_plate, dt)
        E_eddy_cum += P_eddy * dt

        # --- Centrifugal ---
        cent = phys.centrifugal_force(v_sled, m_total)
        occ_g = phys.occupant_g_force(thrust, m_total)

        # --- Power ---
        P_thrust = thrust * v_sled
        P_total = P_thrust + P_eddy

        # --- Record data ---
        slip_ratio = v_slip / (v_sled + v_slip) if (v_sled + v_slip) > 0 else 1.0

        data.record(
            time=t,
            velocity=v_sled,
            position=x_sled,
            acceleration=accel,
            thrust=thrust,
            P_thrust=P_thrust,
            P_eddy=P_eddy,
            P_total=P_total,
            slip_ratio=slip_ratio,
            slip_velocity=v_slip,
            T_plate=T_plate,
            F_centrifugal=cent["F_centrifugal"],
            F_ring_net=cent["F_net"],
            ring_net_g=cent["net_g"],
            occupant_g=occ_g,
            KE=0.5 * m_total * v_sled**2,
            B_field=B,
            f_supply=f_sup,
            f_slip=f_sl,
            current=I,
            voltage=V_pfc,
            skin_depth=delta,
            goodness=G,
            stage_index=stage_idx,
            stage_name=stage.name,
            E_eddy_cumulative=E_eddy_cum,
        )

        # --- Update state ---
        v_sled = v_sled_new
        t += dt
        step += 1

        # Progress report every 60 seconds of sim time
        if step % max(1, int(60.0 / dt)) == 0:
            print(f"  t={t:>8.1f}s  v={v_sled/1000:>7.2f} km/s  "
                  f"stage={stage.name}  F={thrust/1e6:>7.2f} MN  "
                  f"T={T_plate:>7.1f} K  a={accel/cfg.G_ACCEL:>5.2f}g  "
                  f"P={P_total/1e9:>7.2f} GW")

    # --- Final report ---
    print(f"\n{'=' * 72}")
    print("SIMULATION COMPLETE")
    print(f"{'=' * 72}")
    print(f"  Final velocity         {v_sled/1000:>10.3f} km/s")
    print(f"  Total time             {t:>10.1f} s  ({t/60:.1f} min)")
    print(f"  Distance               {x_sled/1000:>10.0f} km  "
          f"({x_sled/cfg.L_RING*100:.2f}% of ring)")
    print(f"  Final plate temp       {T_plate:>10.1f} K  "
          f"(limit: {mat['T_max']} K)")
    print(f"  Total eddy heat        {E_eddy_cum:>12.3e} J  "
          f"({E_eddy_cum/3.6e9:.1f} MWh)")
    print(f"  Final kinetic energy   {0.5*m_total*v_sled**2:>12.3e} J  "
          f"({0.5*m_total*v_sled**2/3.6e12:.1f} GWh)")
    print(f"  Peak thrust            {max(data.thrust)/1e6:>10.2f} MN")
    print(f"  Peak power             {max(data.P_total)/1e9:>10.2f} GW")
    print(f"  Peak plate temp        {max(data.T_plate):>10.1f} K")
    print(f"  Peak occupant g        {max(data.occupant_g):>10.3f} g")
    print(f"  Steps                  {step:>10}")

    if handoff_data:
        print(f"\n  HANDOFF SUMMARY:")
        for h in handoff_data:
            print(f"    {h['from']} → {h['to']} at t={h['time']:.1f}s, "
                  f"v={h['velocity']/1000:.2f} km/s, T={h['T_plate']:.1f} K")

    return data


# =============================================================================
# PLOTTING
# =============================================================================

def make_time_array(data):
    """Convert time to minutes for plotting."""
    return [t / 60.0 for t in data.time]


def stage_colors(data):
    """Return color array based on active stage for scatter-style plots."""
    color_map = {0: '#2196F3', 1: '#FF9800', 2: '#4CAF50'}  # Blue, Orange, Green
    return [color_map.get(s, 'gray') for s in data.stage_index]


def add_handoff_lines(ax, data):
    """Add vertical lines at stage handoff points."""
    prev = data.stage_index[0] if data.stage_index else 0
    for i, s in enumerate(data.stage_index):
        if s != prev:
            t_min = data.time[i] / 60.0
            ax.axvline(t_min, color='red', linestyle='--', alpha=0.4, linewidth=0.8)
            prev = s


def plot_single(data, x_data, y_data, xlabel, ylabel, title, filename,
                yscale='linear', color=None):
    """Generic single-panel plot."""
    if not HAS_MPL:
        print(f"  [skip] {filename} — matplotlib not available")
        return
    fig, ax = plt.subplots(figsize=(12, 6))
    t_min = make_time_array(data)

    if color == 'stage':
        colors = stage_colors(data)
        ax.scatter(t_min, y_data, c=colors, s=1, alpha=0.7)
    else:
        ax.plot(t_min, y_data, color=color or '#1976D2', linewidth=1.2)

    add_handoff_lines(ax, data)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.set_yscale(yscale)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    os.makedirs(GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(GRAPH_DIR, filename), dpi=200)
    plt.close(fig)
    print(f"  Saved {filename}")


def plot_combined(data):
    """Combined 3×3 overview plot showing all stages."""
    if not HAS_MPL:
        return
    t_min = make_time_array(data)

    fig, axes = plt.subplots(3, 3, figsize=(18, 14))
    fig.suptitle("Mass Driver Launch — All Stages", fontsize=15, y=0.98)

    plots = [
        (axes[0, 0], [v/1000 for v in data.velocity], "Velocity (km/s)", '#1976D2'),
        (axes[0, 1], [a/cfg.G_ACCEL for a in data.acceleration], "Acceleration (g)", '#E65100'),
        (axes[0, 2], [F/1e6 for F in data.thrust], "Thrust (MN)", '#2E7D32'),
        (axes[1, 0], [P/1e9 for P in data.P_total], "Total Power (GW)", '#7B1FA2'),
        (axes[1, 1], [P/1e6 for P in data.P_eddy], "Eddy Losses (MW)", '#C62828'),
        (axes[1, 2], data.T_plate, "Plate Temp (K)", '#FF6F00'),
        (axes[2, 0], [s*100 for s in data.slip_ratio], "Slip Ratio (%)", '#00838F'),
        (axes[2, 1], data.B_field, "B Field (T)", '#4527A0'),
        (axes[2, 2], data.ring_net_g, "Ring Net Loading (g)", '#33691E'),
    ]

    for ax, ydata, ylabel, color in plots:
        ax.plot(t_min, ydata, color=color, linewidth=1.0)
        add_handoff_lines(ax, data)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)

    for ax in axes[2]:
        ax.set_xlabel("Time (min)", fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(GRAPH_DIR, "md_combined.png"), dpi=200)
    plt.close(fig)
    print("  Saved md_combined.png")


def plot_stage_detail(data):
    """Per-stage detail plots — 3 columns (LV, MV, HV) showing key parameters."""
    if not HAS_MPL:
        return

    stage_names = ["LV", "MV", "HV"]
    stage_colors_map = ['#2196F3', '#FF9800', '#4CAF50']

    # Split data by stage
    stage_data = {}
    for name in stage_names:
        indices = [i for i, s in enumerate(data.stage_name) if s == name]
        if indices:
            stage_data[name] = indices

    if not stage_data:
        return

    fig, axes = plt.subplots(4, len(stage_data), figsize=(6*len(stage_data), 16))
    if len(stage_data) == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle("Mass Driver — Stage Detail", fontsize=14, y=0.99)

    for col, (name, indices) in enumerate(stage_data.items()):
        t_min = [data.time[i]/60 for i in indices]
        color = stage_colors_map[["LV", "MV", "HV"].index(name)]

        # Velocity
        axes[0, col].plot(t_min, [data.velocity[i]/1000 for i in indices],
                         color=color, lw=1.2)
        axes[0, col].set_ylabel("Velocity (km/s)")
        axes[0, col].set_title(f"Stage {name} (τ_p = {cfg.LIM_STAGES[['LV','MV','HV'].index(name)].tau_p} m)")

        # Current & Voltage
        ax_I = axes[1, col]
        ax_V = ax_I.twinx()
        ax_I.plot(t_min, [data.current[i] for i in indices],
                 color='blue', lw=1.0, label='Current (A)')
        ax_V.plot(t_min, [data.voltage[i]/1000 for i in indices],
                 color='red', lw=1.0, label='Voltage (kV)')
        ax_I.set_ylabel("Current (A)", color='blue')
        ax_V.set_ylabel("Voltage (kV)", color='red')

        # Thrust & Power
        ax_F = axes[2, col]
        ax_P = ax_F.twinx()
        ax_F.plot(t_min, [data.thrust[i]/1e6 for i in indices],
                 color='green', lw=1.0, label='Thrust (MN)')
        ax_P.plot(t_min, [data.P_total[i]/1e9 for i in indices],
                 color='purple', lw=1.0, label='Power (GW)')
        ax_F.set_ylabel("Thrust (MN)", color='green')
        ax_P.set_ylabel("Power (GW)", color='purple')

        # Temperature & Slip
        ax_T = axes[3, col]
        ax_s = ax_T.twinx()
        ax_T.plot(t_min, [data.T_plate[i] for i in indices],
                 color='orange', lw=1.0)
        ax_s.plot(t_min, [data.slip_velocity[i] for i in indices],
                 color='cyan', lw=1.0)
        ax_T.set_ylabel("T_plate (K)", color='orange')
        ax_s.set_ylabel("v_slip (m/s)", color='cyan')
        ax_T.set_xlabel("Time (min)")

    for row in axes:
        for ax in row:
            ax.grid(True, alpha=0.2)
            ax.tick_params(labelsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(GRAPH_DIR, exist_ok=True)
    fig.savefig(os.path.join(GRAPH_DIR, "md_stage_detail.png"), dpi=200)
    plt.close(fig)
    print("  Saved md_stage_detail.png")


GRAPH_REGISTRY = {
    "velocity": lambda d: plot_single(
        d, None, [v/1000 for v in d.velocity],
        "Time (min)", "Velocity (km/s)",
        "Sled Velocity", "md_velocity.png"),
    "accel": lambda d: plot_single(
        d, None, [a/cfg.G_ACCEL for a in d.acceleration],
        "Time (min)", "Acceleration (g)",
        "Sled Acceleration", "md_accel.png", color='#E65100'),
    "thrust": lambda d: plot_single(
        d, None, [F/1e6 for F in d.thrust],
        "Time (min)", "Thrust (MN)",
        "Total Thrust", "md_thrust.png", color='#2E7D32'),
    "power": lambda d: plot_single(
        d, None, [P/1e9 for P in d.P_total],
        "Time (min)", "Power (GW)",
        "Total Electrical Power", "md_power.png", color='#7B1FA2'),
    "eddy": lambda d: plot_single(
        d, None, [P/1e6 for P in d.P_eddy],
        "Time (min)", "Eddy Losses (MW)",
        "Eddy Current Losses in Plate", "md_eddy.png", color='#C62828'),
    "slip": lambda d: plot_single(
        d, None, [s*100 for s in d.slip_ratio],
        "Time (min)", "Slip Ratio (%)",
        "Slip Ratio", "md_slip.png", color='#00838F'),
    "temp": lambda d: plot_single(
        d, None, d.T_plate,
        "Time (min)", "Temperature (K)",
        "Reaction Plate Temperature", "md_temp.png", color='#FF6F00'),
    "centrifugal": lambda d: plot_single(
        d, None, [F/1e6 for F in d.F_centrifugal],
        "Time (min)", "Centrifugal Force (MN)",
        "Centrifugal Force on Ring", "md_centrifugal.png", color='#795548'),
    "ring_force": lambda d: plot_single(
        d, None, [F/1e6 for F in d.F_ring_net],
        "Time (min)", "Net Ring Force (MN)",
        "Net Force on Orbital Ring (+ = outward)", "md_ring_force.png"),
    "occupant_g": lambda d: plot_single(
        d, None, d.occupant_g,
        "Time (min)", "Occupant g-force",
        "Occupant G-Force (tangential)", "md_occupant_g.png", color='#E65100'),
    "energy": lambda d: plot_single(
        d, None, [E/1e12 for E in d.KE],
        "Time (min)", "Kinetic Energy (TJ)",
        "Sled Kinetic Energy", "md_energy.png", color='#0D47A1'),
    "b_field": lambda d: plot_single(
        d, None, d.B_field,
        "Time (min)", "B Field (T)",
        "Peak Magnetic Field at Plate", "md_b_field.png", color='#4527A0'),
    "frequency": lambda d: plot_single(
        d, None, d.f_supply,
        "Time (min)", "Frequency (Hz)",
        "Supply Frequency", "md_frequency.png", color='#00695C'),
    "current": lambda d: plot_single(
        d, None, d.current,
        "Time (min)", "Current (A)",
        "Phase Current", "md_current.png", color='#1565C0'),
    "voltage": lambda d: plot_single(
        d, None, [V/1000 for V in d.voltage],
        "Time (min)", "Voltage (kV)",
        "Back-EMF Voltage", "md_voltage.png", color='#B71C1C'),
    "skin_depth": lambda d: plot_single(
        d, None, [sd*1000 for sd in d.skin_depth],
        "Time (min)", "Skin Depth (mm)",
        "Skin Depth", "md_skin_depth.png", color='#4E342E'),
    "goodness": lambda d: plot_single(
        d, None, d.goodness,
        "Time (min)", "Goodness Factor G",
        "Goodness Factor", "md_goodness.png", color='#37474F'),
    "combined": lambda d: plot_combined(d),
    "stage_detail": lambda d: plot_stage_detail(d),
}


def generate_graphs(data, graph_keys):
    """Generate requested graphs."""
    if not HAS_MPL:
        print("\nMatplotlib not available — skipping graphs.")
        return

    print(f"\nGenerating graphs...")
    if "all" in graph_keys:
        graph_keys = list(GRAPH_REGISTRY.keys())

    for key in graph_keys:
        if key in GRAPH_REGISTRY:
            try:
                GRAPH_REGISTRY[key](data)
            except Exception as e:
                print(f"  [error] {key}: {e}")
        else:
            print(f"  [unknown] {key}")


# =============================================================================
# CSV EXPORT
# =============================================================================

def save_csv(data, filename="md_launch_data.csv"):
    """Save time-series data to CSV."""
    path = os.path.join(GRAPH_DIR, filename)
    os.makedirs(GRAPH_DIR, exist_ok=True)

    headers = [
        "time_s", "velocity_ms", "position_m", "acceleration_ms2",
        "thrust_N", "P_thrust_W", "P_eddy_W", "P_total_W",
        "slip_ratio", "slip_velocity_ms", "T_plate_K",
        "F_centrifugal_N", "F_ring_net_N", "ring_net_g", "occupant_g",
        "KE_J", "B_field_T", "f_supply_Hz", "f_slip_Hz",
        "current_A", "voltage_V", "skin_depth_m", "goodness",
        "stage_index", "stage_name", "E_eddy_cum_J"
    ]

    with open(path, 'w') as f:
        f.write(",".join(headers) + "\n")
        for i in range(len(data.time)):
            row = [
                data.time[i], data.velocity[i], data.position[i],
                data.acceleration[i], data.thrust[i], data.P_thrust[i],
                data.P_eddy[i], data.P_total[i], data.slip_ratio[i],
                data.slip_velocity[i], data.T_plate[i],
                data.F_centrifugal[i], data.F_ring_net[i],
                data.ring_net_g[i], data.occupant_g[i],
                data.KE[i], data.B_field[i], data.f_supply[i],
                data.f_slip[i], data.current[i], data.voltage[i],
                data.skin_depth[i], data.goodness[i],
                data.stage_index[i], data.stage_name[i],
                data.E_eddy_cumulative[i],
            ]
            f.write(",".join(str(x) for x in row) + "\n")

    print(f"  Saved {path}")


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    model = None
    v_launch = None
    max_accel = None
    mass = None
    quick = False
    no_graphs = False
    save_csv_flag = False
    graph_keys = []

    for arg in sys.argv[1:]:
        if arg.startswith("--model="):
            model = int(arg.split("=")[1])
        elif arg.startswith("--v_launch="):
            v_launch = float(arg.split("=")[1])
        elif arg.startswith("--accel="):
            max_accel = float(arg.split("=")[1])
        elif arg.startswith("--mass="):
            mass = float(arg.split("=")[1]) * 1000  # tonnes → kg
        elif arg == "--quick":
            quick = True
        elif arg == "--no-graphs":
            no_graphs = True
        elif arg == "--save-csv":
            save_csv_flag = True
        elif arg == "--help":
            print(__doc__)
            sys.exit(0)
        else:
            graph_keys.append(arg)

    # Default: show combined + stage detail
    if not graph_keys and not no_graphs:
        graph_keys = ["combined", "stage_detail"]

    return {
        "model": model,
        "v_launch": v_launch,
        "max_accel_g": max_accel,
        "m_spacecraft": mass,
        "quick": quick,
        "no_graphs": no_graphs,
        "save_csv": save_csv_flag,
        "graph_keys": graph_keys,
    }


def main():
    args = parse_args()

    data = run_simulation(
        model=args["model"],
        v_launch=args["v_launch"],
        max_accel_g=args["max_accel_g"],
        m_spacecraft=args["m_spacecraft"],
        quick=args["quick"],
    )

    if not args["no_graphs"]:
        generate_graphs(data, args["graph_keys"])

    if args["save_csv"]:
        save_csv(data)


if __name__ == "__main__":
    main()
