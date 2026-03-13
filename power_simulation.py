#!/usr/bin/env python3
"""
Power Generation and HVDC Distribution Simulation

Computes solar power generation around the orbital ring and sizes the
circumferential HVDC transmission system for power distribution from
dayside to nightside.

Reference: "Orbital Ring Engineering" by Paul G de Jong

Usage:
    python power_simulation.py [options] [graph_keywords...]

Options:
    --width=N       Panel width in metres (default: 53)
    --voltage=N     HVDC pole-to-pole voltage in MV (default: 10)
    --loss=N        Loss budget as fraction (default: 0.05)
    --diameter=N    Cable diameter in metres (default: 0.10)
    --demand=MODE   Demand mode: "deployment" or "ops" (default: deployment)
    --power=N       LIM power per site in MW (default: 8, use 16 for high-power)
    --save          Save graphs to files (default)
    --show          Display graphs interactively
    --outdir=PATH   Graph output directory (default: ./graphs_power)
    --dpi=N         Graph resolution (default: 300)
    --help, -h      Show this help

Graph keywords:
    ring_profile    Generation, demand, net, and power flow vs angle
    panel_sweep     Parametric panel width sweep (generation, flow, conductor)
    voltage_sweep   Parametric HVDC voltage sweep (conductor, mass, current)
    all             Show all graphs
"""

import sys
import os
import math

import matplotlib.pyplot as plt

import power_config as cfg
import power_physics as phys


def format_power(watts):
    """Format power value with appropriate SI prefix."""
    if abs(watts) >= 1e12:
        return f"{watts/1e12:.2f} TW"
    elif abs(watts) >= 1e9:
        return f"{watts/1e9:.2f} GW"
    elif abs(watts) >= 1e6:
        return f"{watts/1e6:.2f} MW"
    elif abs(watts) >= 1e3:
        return f"{watts/1e3:.2f} kW"
    else:
        return f"{watts:.2f} W"


def format_mass(kg):
    """Format mass value with appropriate prefix."""
    if abs(kg) >= 1e9:
        return f"{kg/1e9:.2f} Mt"
    elif abs(kg) >= 1e6:
        return f"{kg/1e6:.2f} kt"
    elif abs(kg) >= 1e3:
        return f"{kg/1e3:.2f} t"
    else:
        return f"{kg:.2f} kg"


def run_analysis(panel_width, demand_mode, v_hvdc, loss_budget, lim_power=None):
    """Run complete power analysis for a single configuration.

    Returns dict with all results.
    """
    # Solar generation profile
    profile = phys.compute_ring_power_profile(panel_width)

    # Demand
    if lim_power is None:
        lim_power = cfg.LIM_POWER_PER_SITE

    if demand_mode == "ops":
        p_demand_total = cfg.P_DEMAND_OPS
        demand_label = "Post-deployment ops"
    else:
        p_demand_total = cfg.LIM_SITES * lim_power
        demand_label = f"Deployment (LIM {lim_power/1e6:.0f} MW)"

    demand_per_m = p_demand_total / cfg.L_RING

    # Power flow
    dphi = 2 * math.pi / len(profile['phi'])
    flow = phys.compute_power_flow(profile['p_gen_per_m'], demand_per_m, dphi)

    # Conductor sizing
    conductor = phys.size_conductor(flow['p_flow_peak'], v_hvdc, loss_budget)

    # Panel width needed to match demand
    if profile['p_avg_per_m2'] > 0:
        width_needed = p_demand_total / (profile['p_avg_per_m2'] * cfg.L_RING)
    else:
        width_needed = float('inf')

    return {
        'profile': profile,
        'flow': flow,
        'conductor': conductor,
        'demand_per_m': demand_per_m,
        'demand_label': demand_label,
        'p_demand_total': p_demand_total,
        'width_needed': width_needed,
        'panel_width': panel_width,
    }


def print_results(results):
    """Print formatted results to console."""
    streq = "=" * 70

    profile = results['profile']
    flow = results['flow']
    cond = results['conductor']

    # --- Solar Generation ---
    print(f"\n{streq}")
    print("SOLAR POWER GENERATION")
    print(streq)
    print(f"  {'Panel width':24} {results['panel_width']:>12.1f} m")
    print(f"  {'Cell efficiency':24} {cfg.CELL_EFFICIENCY*100:>12.1f} %")
    print(f"  {'Ring circumference':24} {cfg.L_RING/1e3:>12.0f} km")
    print(f"  {'Shadow half-angle':24} {math.degrees(cfg.SHADOW_HALF_ANGLE):>12.1f} deg")
    print()
    print(f"  {'Peak output (noon)':24} {profile['p_peak_per_m2']:>12.1f} W/m^2")
    print(f"  {'Orbit-averaged output':24} {profile['p_avg_per_m2']:>12.1f} W/m^2")
    print(f"  {'Reference (Chapter 7)':24} {cfg.P_AVG_REFERENCE:>12.1f} W/m^2")

    ratio = profile['p_avg_per_m2'] / cfg.P_AVG_REFERENCE
    if abs(ratio - 1.0) > 0.05:
        print(f"  *** MISMATCH: computed average is {ratio:.3f}x reference ***")
        print(f"  *** First-principles gives {profile['p_avg_per_m2']:.1f} W/m^2, "
              f"reference is {cfg.P_AVG_REFERENCE:.1f} W/m^2 ***")
    else:
        print(f"  {'Match ratio':24} {ratio:>12.3f}x (OK)")

    print()
    print(f"  {'Total generation':24} {format_power(profile['p_total']):>12}")
    print(f"  {'Generation per metre':24} {profile['p_total']/cfg.L_RING:>12.0f} W/m")

    # --- Power Demand ---
    print(f"\n{streq}")
    print("POWER DEMAND")
    print(streq)
    print(f"  {'Mode':24} {results['demand_label']:>12}")
    print(f"  {'LIM sites':24} {cfg.LIM_SITES:>12,}")
    if results['demand_label'].startswith("Post"):
        print(f"  {'Power per site':24} {cfg.OPS_POWER_PER_SITE/1e3:>12.1f} kW")
    else:
        site_power = results['p_demand_total'] / cfg.LIM_SITES
        print(f"  {'Power per site':24} {site_power/1e6:>12.1f} MW")
    print(f"  {'Total demand':24} {format_power(results['p_demand_total']):>12}")
    print(f"  {'Demand per metre':24} {results['demand_per_m']:>12.0f} W/m")

    surplus = profile['p_total'] - results['p_demand_total']
    print()
    print(f"  {'Generation - Demand':24} {format_power(surplus):>12}")
    if surplus >= 0:
        print(f"  {'Surplus ratio':24} {profile['p_total']/results['p_demand_total']:>12.2f}x")
    else:
        print(f"  {'Width needed to match':24} {results['width_needed']:>12.1f} m")

    # --- HVDC Transmission ---
    print(f"\n{streq}")
    print("HVDC CIRCUMFERENTIAL TRANSMISSION")
    print(streq)
    print(f"  {'HVDC voltage (p-to-p)':24} {cond['v_hvdc']/1e6:>12.1f} MV")
    print(f"  {'Per-pole voltage':24} {cond['v_pole']/1e6:>12.1f} MV")
    print(f"  {'Loss budget':24} {cond['loss_budget']*100:>12.1f} %")
    print(f"  {'Avg transmission dist':24} {cond['l_trans']/1e3:>12.0f} km")
    print()
    print(f"  {'Peak power flow':24} {format_power(flow['p_flow_peak']):>12}")
    print(f"  {'Peak current per pole':24} {cond['i_peak']:>12,.0f} A")
    print(f"  {'Current density':24} {cond['j_peak']/1e6:>12.2f} MA/m^2")
    print()
    print(f"  {'Conductor cross-sect':24} {cond['a_total']*1e4:>12.2f} cm^2")
    print(f"  {'Per-pole cross-sect':24} {cond['a_per_pole']*1e4:>12.2f} cm^2")
    print(f"  {'Cable diameter':24} {cfg.CABLE_DIAMETER*100:>12.1f} cm")
    print(f"  {'Cables per pole':24} {cond['n_cables_per_pole']:>12,}")
    print(f"  {'Total cables':24} {cond['n_cables_total']:>12,}")
    print()
    print(f"  {'Conductor mass/m':24} {cond['m_per_m']:>12.2f} kg/m")
    print(f"  {'Total conductor mass':24} {format_mass(cond['m_total']):>12}")
    print()
    print(f"  {'Peak loss':24} {format_power(cond['p_loss_peak']):>12}")
    print(f"  {'Actual loss fraction':24} {cond['loss_fraction']*100:>12.2f} %")

    # --- Key Configuration Summary ---
    print(f"\n{streq}")
    print("KEY CONFIGURATIONS CROSS-CHECK")
    print(streq)

    configs = [
        ("Deployment 8 MW (53 m)", 53.0, 666e9),
        ("Deployment 16 MW (107 m)", 107.0, 1332e9),
        ("Expanded (400 m)", 400.0, None),
    ]

    print(f"  {'Configuration':26} {'Width':>8} {'Generation':>14} {'Target':>14}")
    print(f"  {'-'*26} {'-'*8} {'-'*14} {'-'*14}")
    for label, w, target in configs:
        p = phys.compute_ring_power_profile(w)
        tgt = format_power(target) if target else "---"
        print(f"  {label:26} {w:>7.0f}m {format_power(p['p_total']):>14} {tgt:>14}")

    print(f"{streq}\n")


def plot_ring_profile(results, show_graphs):
    """Plot ring profile: generation, demand, net, and power flow vs angle."""
    profile = results['profile']
    flow = results['flow']

    phi_deg = [math.degrees(p) for p in profile['phi']]
    gen_kw = [g / 1e3 for g in profile['p_gen_per_m']]
    demand_kw = [results['demand_per_m'] / 1e3] * len(phi_deg)
    net_kw = [n / 1e3 for n in flow['net_per_m']]
    flow_gw = [f / 1e9 for f in flow['p_flow']]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(cfg.GRAPH_WIDTH_INCHES,
                                                    cfg.GRAPH_HEIGHT_INCHES * 2))

    # Top: generation, demand, net
    ax1.plot(phi_deg, gen_kw, 'gold', linewidth=1.5, label='Generation')
    ax1.plot(phi_deg, demand_kw, 'r--', linewidth=1.0, label='Demand')
    ax1.fill_between(phi_deg, net_kw, 0, where=[n > 0 for n in net_kw],
                     color='green', alpha=0.3, label='Surplus')
    ax1.fill_between(phi_deg, net_kw, 0, where=[n <= 0 for n in net_kw],
                     color='red', alpha=0.3, label='Deficit')
    ax1.set_ylabel("Power per metre (kW/m)", fontsize=14)
    ax1.set_title(f"Ring Power Profile  |  Panel width = {results['panel_width']:.0f} m  |  "
                  f"{results['demand_label']}", fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-180, 180)

    # Bottom: cumulative power flow
    ax2.plot(phi_deg, flow_gw, 'blue', linewidth=1.5)
    ax2.axhline(0, color='gray', linewidth=0.5)
    ax2.set_xlabel("Angle from subsolar point (degrees)", fontsize=14)
    ax2.set_ylabel("Circumferential power flow (GW)", fontsize=14)
    ax2.set_title(f"HVDC Transmission Load  |  Peak = {format_power(flow['p_flow_peak'])}",
                  fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-180, 180)

    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"01-ring_profile.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


def plot_panel_sweep(demand_mode, v_hvdc, loss_budget, lim_power, show_graphs):
    """Parametric sweep: panel width vs generation, flow, conductor sizing."""
    widths = []
    gen_gw = []
    flow_gw = []
    a_cm2 = []
    mass_kg_m = []
    loss_pct = []

    for i in range(cfg.PANEL_WIDTH_STEPS + 1):
        w = cfg.PANEL_WIDTH_MIN + (cfg.PANEL_WIDTH_MAX - cfg.PANEL_WIDTH_MIN) * i / cfg.PANEL_WIDTH_STEPS
        profile = phys.compute_ring_power_profile(w)

        if demand_mode == "ops":
            p_demand = cfg.P_DEMAND_OPS
        else:
            p_demand = cfg.LIM_SITES * lim_power

        demand_per_m = p_demand / cfg.L_RING
        dphi = 2 * math.pi / len(profile['phi'])
        flow_result = phys.compute_power_flow(profile['p_gen_per_m'], demand_per_m, dphi)
        cond = phys.size_conductor(flow_result['p_flow_peak'], v_hvdc, loss_budget)

        widths.append(w)
        gen_gw.append(profile['p_total'] / 1e9)
        flow_gw.append(flow_result['p_flow_peak'] / 1e9)
        a_cm2.append(cond['a_total'] * 1e4)
        mass_kg_m.append(cond['m_per_m'])
        loss_pct.append(cond['loss_fraction'] * 100)

    fig, axes = plt.subplots(2, 2, figsize=(cfg.GRAPH_WIDTH_INCHES * 1.2,
                                             cfg.GRAPH_HEIGHT_INCHES * 2))

    ax = axes[0, 0]
    ax.plot(widths, gen_gw, 'gold', linewidth=2)
    ax.set_ylabel("Total generation (GW)", fontsize=12)
    ax.set_title("Generation vs Panel Width", fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(widths, flow_gw, 'blue', linewidth=2)
    ax.set_ylabel("Peak HVDC flow (GW)", fontsize=12)
    ax.set_title("Peak Transmission Load", fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(widths, a_cm2, 'darkgreen', linewidth=2)
    ax.set_ylabel("Conductor area (cm^2)", fontsize=12)
    ax.set_title("Required Conductor Cross-Section", fontsize=13)
    ax.set_xlabel("Panel width (m)", fontsize=12)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(widths, mass_kg_m, 'brown', linewidth=2)
    ax.set_ylabel("Conductor mass (kg/m)", fontsize=12)
    ax.set_title("Conductor Mass per Metre", fontsize=13)
    ax.set_xlabel("Panel width (m)", fontsize=12)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"Panel Width Parametric Sweep  |  V_HVDC = {v_hvdc/1e6:.0f} MV  |  "
                 f"Loss budget = {loss_budget*100:.0f}%", fontsize=14, y=1.02)
    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"02-panel_sweep.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


def plot_voltage_sweep(panel_width, demand_mode, loss_budget, lim_power, show_graphs):
    """Parametric sweep: HVDC voltage vs conductor sizing."""
    # Compute power flow once for the fixed panel width
    profile = phys.compute_ring_power_profile(panel_width)

    if demand_mode == "ops":
        p_demand = cfg.P_DEMAND_OPS
    else:
        p_demand = cfg.LIM_SITES * lim_power

    demand_per_m = p_demand / cfg.L_RING
    dphi = 2 * math.pi / len(profile['phi'])
    flow_result = phys.compute_power_flow(profile['p_gen_per_m'], demand_per_m, dphi)

    voltages_mv = []
    a_cm2 = []
    mass_kg_m = []
    current_ka = []
    loss_pct = []

    for i in range(cfg.V_HVDC_STEPS + 1):
        v = cfg.V_HVDC_MIN + (cfg.V_HVDC_MAX - cfg.V_HVDC_MIN) * i / cfg.V_HVDC_STEPS
        cond = phys.size_conductor(flow_result['p_flow_peak'], v, loss_budget)

        voltages_mv.append(v / 1e6)
        a_cm2.append(cond['a_total'] * 1e4)
        mass_kg_m.append(cond['m_per_m'])
        current_ka.append(cond['i_peak'] / 1e3)
        loss_pct.append(cond['loss_fraction'] * 100)

    fig, axes = plt.subplots(2, 2, figsize=(cfg.GRAPH_WIDTH_INCHES * 1.2,
                                             cfg.GRAPH_HEIGHT_INCHES * 2))

    ax = axes[0, 0]
    ax.plot(voltages_mv, a_cm2, 'darkgreen', linewidth=2)
    ax.set_ylabel("Conductor area (cm^2)", fontsize=12)
    ax.set_title("Conductor Cross-Section vs Voltage", fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(voltages_mv, mass_kg_m, 'brown', linewidth=2)
    ax.set_ylabel("Conductor mass (kg/m)", fontsize=12)
    ax.set_title("Conductor Mass per Metre", fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(voltages_mv, current_ka, 'red', linewidth=2)
    ax.set_ylabel("Peak current per pole (kA)", fontsize=12)
    ax.set_title("Peak Current vs Voltage", fontsize=13)
    ax.set_xlabel("HVDC voltage (MV, pole-to-pole)", fontsize=12)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(voltages_mv, loss_pct, 'purple', linewidth=2)
    ax.set_ylabel("Loss fraction (%)", fontsize=12)
    ax.set_title("Transmission Loss vs Voltage", fontsize=13)
    ax.set_xlabel("HVDC voltage (MV, pole-to-pole)", fontsize=12)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"HVDC Voltage Parametric Sweep  |  Panel = {panel_width:.0f} m  |  "
                 f"Peak flow = {format_power(flow_result['p_flow_peak'])}", fontsize=14, y=1.02)
    plt.tight_layout()

    if cfg.SAVE_GRAPHS:
        filename = f"03-voltage_sweep.{cfg.GRAPH_FORMAT}"
        filepath = os.path.join(cfg.GRAPH_OUTPUT_DIR, filename)
        plt.savefig(filepath, dpi=cfg.GRAPH_DPI, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {filename}")
        plt.close(fig)
    else:
        plt.show()


def main():
    """Main entry point."""
    # Defaults
    panel_width = cfg.PANEL_WIDTH
    v_hvdc = cfg.V_HVDC
    loss_budget = cfg.LOSS_BUDGET
    demand_mode = "deployment"
    lim_power = cfg.LIM_POWER_PER_SITE
    show_graphs = []

    # Parse CLI arguments
    if len(sys.argv) > 1:
        if "--help" in sys.argv or "-h" in sys.argv:
            print(__doc__)
            return

        for arg in sys.argv[1:]:
            if arg.startswith("--width="):
                panel_width = float(arg.split("=")[1])
            elif arg.startswith("--voltage="):
                v_hvdc = float(arg.split("=")[1]) * 1e6  # Input in MV
            elif arg.startswith("--loss="):
                loss_budget = float(arg.split("=")[1])
            elif arg.startswith("--diameter="):
                cfg.CABLE_DIAMETER = float(arg.split("=")[1])
            elif arg.startswith("--demand="):
                demand_mode = arg.split("=")[1].lower()
            elif arg.startswith("--power="):
                lim_power = float(arg.split("=")[1]) * 1e6  # Input in MW
            elif arg == "--save":
                cfg.SAVE_GRAPHS = True
            elif arg == "--show":
                cfg.SAVE_GRAPHS = False
            elif arg.startswith("--outdir="):
                cfg.GRAPH_OUTPUT_DIR = arg.split("=")[1]
            elif arg.startswith("--dpi="):
                cfg.GRAPH_DPI = int(arg.split("=")[1])
            else:
                show_graphs.append(arg)

    # Create output directory
    if cfg.SAVE_GRAPHS:
        os.makedirs(cfg.GRAPH_OUTPUT_DIR, exist_ok=True)

    # Run analysis
    results = run_analysis(panel_width, demand_mode, v_hvdc, loss_budget, lim_power)

    # Print results
    print_results(results)

    # Graphs
    if show_graphs:
        print("Generating graphs...")
        if "ring_profile" in show_graphs or "all" in show_graphs:
            plot_ring_profile(results, show_graphs)
        if "panel_sweep" in show_graphs or "all" in show_graphs:
            plot_panel_sweep(demand_mode, v_hvdc, loss_budget, lim_power, show_graphs)
        if "voltage_sweep" in show_graphs or "all" in show_graphs:
            plot_voltage_sweep(panel_width, demand_mode, loss_budget, lim_power, show_graphs)
        print("Done.")


if __name__ == "__main__":
    main()
