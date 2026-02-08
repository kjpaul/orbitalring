#!/usr/bin/env python3
"""
Mass Driver Parametric Trade Study

Sweeps sled length (1–5 km) and plate width (0.5–2.0 m) at fixed 200mm
thickness, computing for each combination:
  - Plate mass, sled mass, payload capacity
  - Number of simultaneous LIM interactions
  - Electromagnetic coupling (goodness factor, skin depth)
  - Required slip velocity and thermal budget
  - Whether the design closes (thrust meets requirement within thermal limit)
  - Payload fraction

Reference mission: Mars transfer, 15 km/s at 0.5g

Usage: python md_trade_study.py
"""

import math
import os
import sys

# ── Physical constants ──────────────────────────────────────────────────────
MU0 = 4 * math.pi * 1e-7
G0 = 9.81

# ── Orbital ring geometry ───────────────────────────────────────────────────
R_EARTH_EQ = 6_378_137.0
R_ORBIT = R_EARTH_EQ + 250_000
V_ORBIT = 7754.866
L_RING = 41_645_813.012

# ── γ-TiAl properties ──────────────────────────────────────────────────────
RHO_293K = 75e-8          # resistivity at 293 K (Ω·m)
DENSITY = 3900            # kg/m³
C_P = 570                 # J/(kg·K)
ALPHA_R = 0.0012          # temperature coefficient of resistivity
EM = 0.60                 # emissivity

# ── Fixed design parameters ─────────────────────────────────────────────────
PLATE_THICKNESS = 0.200   # m (fixed at 200 mm)
N_PLATES = 2              # one per stator face
RIB_FACTOR = 1.15         # 15% mass adder for reinforcing ribs
B_PLATE = 1.5             # T — traveling wave field at plate surface
GAP = 0.10                # m — air gap

# Mission
V_TARGET = 15_000.0       # m/s — Mars transfer
ACCEL = 4.905             # m/s² — 0.5g

# Temperature limits
T_INITIAL = 293.0         # K
T_MAX = 800.0             # K
DELTA_T_MAX = T_MAX - T_INITIAL  # 507 K

# LIM stage layout (repeating unit)
# Stage 3: τ_p=200m, 3 poles → 600m active + 134m overlap = 734m
# Stage 2: τ_p=30m,  3 poles →  90m active +  20m overlap = 110m
# Stage 1: τ_p=3m,   5 poles →  15m active +   2m overlap =  17m
S3_TAU_P = 200.0;  S3_L_ACTIVE = 600.0;  S3_L_SEG = 734.0
S2_TAU_P = 30.0;   S2_L_ACTIVE =  90.0;  S2_L_SEG = 110.0
S1_TAU_P = 3.0;    S1_L_ACTIVE =  15.0;  S1_L_SEG =  17.0

L_UNIT = S3_L_SEG + S2_L_SEG + S1_L_SEG  # 861 m
DUTY_S3 = S3_L_SEG / L_UNIT
DUTY_S2 = S2_L_SEG / L_UNIT
DUTY_S1 = S1_L_SEG / L_UNIT

# Sled structure mass scales with length
STRUCTURE_MASS_PER_M = 20.0   # kg/m — trusses, rails, clamps, fairings

# Minimum spacecraft mass (what we're trying to maximize)
SPACECRAFT_MASS_MIN = 100_000  # 100 t minimum useful payload


# ── Physics functions ────────────────────────────────────────────────────────

def resistivity(T):
    return RHO_293K * (1 + ALPHA_R * (T - 293))

def skin_depth(tau_p, v_slip, T):
    rho = resistivity(T)
    f_slip = max(v_slip, 0.1) / (2 * tau_p)
    return math.sqrt(rho / (math.pi * MU0 * max(f_slip, 0.01)))

def goodness_factor(tau_p, v_slip, T, t_plate):
    if v_slip <= 0:
        return 0.0
    rho = resistivity(T)
    sigma = 1.0 / rho
    f_slip = v_slip / (2 * tau_p)
    omega = 2 * math.pi * f_slip
    delta = skin_depth(tau_p, v_slip, T)
    d_eff = min(t_plate, delta)
    return (omega * MU0 * sigma * d_eff * tau_p) / math.pi

def slip_efficiency(s, G):
    if s <= 0 or G <= 0:
        return 0.0
    sG = s * G
    return (2 * sG) / (1 + sG**2)

def thrust_one_segment(tau_p, l_active, w_coil, v_sled, v_slip, T, t_plate, n_plates):
    """Thrust from one LIM segment (one pole-pitch group)."""
    if v_slip <= 0 or v_sled <= 0:
        return 0.0
    v_sync = v_sled + v_slip
    s = v_slip / v_sync
    G = goodness_factor(tau_p, v_slip, T, t_plate)
    eta = slip_efficiency(s, G)
    a_lim = l_active * w_coil * n_plates
    F_max = (B_PLATE**2) * a_lim / (2 * MU0)
    return F_max * eta

def total_thrust_at_speed(v_sled, v_slip, T, w_coil, plate_length):
    """
    Total thrust from all stages, accounting for how many repeating
    units the sled plate overlaps simultaneously.
    
    At high speed (v > 3000 m/s), only Stage 3 is active.
    At Mars transfer speed, Stage 3 dominates.
    """
    plate_units = plate_length / L_UNIT  # how many repeating units sled spans
    
    # Stage 3 contribution (dominant at high speed)
    F_s3_seg = thrust_one_segment(S3_TAU_P, S3_L_ACTIVE, w_coil,
                                   v_sled, v_slip, T, PLATE_THICKNESS, N_PLATES)
    n_s3_active = plate_units * DUTY_S3
    F_s3 = F_s3_seg * n_s3_active
    
    # Stage 2 (active 300-3000 m/s, contributes at transition)
    F_s2 = 0.0
    if v_sled < 4000:
        F_s2_seg = thrust_one_segment(S2_TAU_P, S2_L_ACTIVE, w_coil,
                                       v_sled, v_slip, T, PLATE_THICKNESS, N_PLATES)
        n_s2_active = plate_units * DUTY_S2
        F_s2 = F_s2_seg * n_s2_active
    
    # Stage 1 (active 0-300 m/s)
    F_s1 = 0.0
    if v_sled < 400:
        F_s1_seg = thrust_one_segment(S1_TAU_P, S1_L_ACTIVE, w_coil,
                                       v_sled, v_slip, T, PLATE_THICKNESS, N_PLATES)
        n_s1_active = plate_units * DUTY_S1
        F_s1 = F_s1_seg * n_s1_active
    
    return F_s3 + F_s2 + F_s1, F_s3, F_s2, F_s1

def find_slip_for_thrust(v_sled, T, w_coil, plate_length, F_required):
    """Binary search for v_slip that produces required thrust at given speed."""
    v_lo = 0.1
    v_hi = min(500.0, v_sled * 0.05)  # cap at 5% slip ratio
    v_hi = max(v_hi, 1.0)
    
    # Check if max slip is enough
    F_max, _, _, _ = total_thrust_at_speed(v_sled, v_hi, T, w_coil, plate_length)
    if F_max < F_required:
        return v_hi, F_max, False  # can't meet thrust requirement
    
    for _ in range(80):
        v_mid = (v_lo + v_hi) / 2
        F_mid, _, _, _ = total_thrust_at_speed(v_sled, v_mid, T, w_coil, plate_length)
        if F_mid < F_required:
            v_lo = v_mid
        else:
            v_hi = v_mid
        if abs(v_hi - v_lo) < 0.001:
            break
    
    v_slip = (v_lo + v_hi) / 2
    F_actual, _, _, _ = total_thrust_at_speed(v_sled, v_slip, T, w_coil, plate_length)
    return v_slip, F_actual, True


# ── Simplified launch thermal simulation ─────────────────────────────────────

def simulate_launch_thermal(sled_mass, w_coil, plate_length, plate_mass):
    """
    Run a simplified launch simulation to check if thermal budget closes.
    
    Steps through velocity from 1 m/s to V_TARGET at constant acceleration,
    tracking cumulative heat deposited in plate from slip losses.
    
    Returns (success, T_final, E_heat, max_slip_pct, details_dict)
    """
    F_required = sled_mass * ACCEL
    E_kinetic = 0.5 * sled_mass * V_TARGET**2
    thermal_capacity = plate_mass * C_P * DELTA_T_MAX
    
    # Quick analytical check first
    slip_budget = thermal_capacity / E_kinetic if E_kinetic > 0 else 1.0
    
    # Step through representative velocity points
    # Use coarse steps for speed — this is a trade study, not final sim
    T_plate = T_INITIAL
    E_heat = 0.0
    max_slip_pct = 0.0
    thrust_ok = True
    
    # Velocity steps: finer at low speed, coarser at high speed
    v_points = []
    v = 1.0
    while v < V_TARGET:
        v_points.append(v)
        if v < 100:
            v += 10
        elif v < 1000:
            v += 50
        elif v < 5000:
            v += 200
        else:
            v += 500
    v_points.append(V_TARGET)
    
    for i in range(len(v_points) - 1):
        v = v_points[i]
        v_next = v_points[i + 1]
        dv = v_next - v
        dt_seg = dv / ACCEL  # time for this velocity segment
        v_mid = (v + v_next) / 2
        
        # Find required slip at midpoint velocity
        v_slip, F_actual, ok = find_slip_for_thrust(v_mid, T_plate, w_coil,
                                                      plate_length, F_required)
        if not ok:
            thrust_ok = False
        
        slip_pct = v_slip / v_mid * 100 if v_mid > 0 else 0
        max_slip_pct = max(max_slip_pct, slip_pct)
        
        # Heat deposited in this segment
        P_slip = F_actual * v_slip
        dE = P_slip * dt_seg
        E_heat += dE
        
        # Temperature rise (adiabatic — conservative)
        if plate_mass > 0:
            T_plate += dE / (plate_mass * C_P)
        
        if T_plate >= T_MAX:
            return False, T_plate, E_heat, max_slip_pct, {
                'v_fail': v_mid, 'slip_budget': slip_budget,
                'F_required': F_required, 'E_kinetic': E_kinetic,
                'thermal_capacity': thermal_capacity, 'thrust_ok': thrust_ok,
            }
    
    return thrust_ok, T_plate, E_heat, max_slip_pct, {
        'slip_budget': slip_budget,
        'F_required': F_required, 'E_kinetic': E_kinetic,
        'thermal_capacity': thermal_capacity, 'thrust_ok': thrust_ok,
        'eff_slip': E_heat / E_kinetic if E_kinetic > 0 else 0,
    }


# ── Trade study sweep ────────────────────────────────────────────────────────

def run_trade_study():
    print("=" * 90)
    print("MASS DRIVER PARAMETRIC TRADE STUDY")
    print(f"Mission: {V_TARGET/1000:.0f} km/s Mars transfer at {ACCEL/G0:.1f}g")
    print(f"Plate: 200mm γ-TiAl, B={B_PLATE:.1f} T, 2 plates per sled")
    print(f"Ring: {L_RING/1000:,.0f} km, repeating unit {L_UNIT:.0f} m")
    print("=" * 90)
    
    # Sweep ranges
    sled_lengths = [1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]  # m
    plate_widths = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]  # m
    
    results = []
    
    print(f"\n{'L_sled':>7s} {'W_plate':>7s} │ {'m_plate':>8s} {'m_struct':>8s} "
          f"{'m_sled':>8s} {'F_req':>7s} │ {'Th_cap':>7s} {'E_kin':>7s} "
          f"{'s_bud':>6s} │ {'#LIM':>5s} │ {'Closes':>6s} {'T_fin':>6s} "
          f"{'m_pay':>8s} {'η_pay':>6s}")
    print(f"{'(m)':>7s} {'(m)':>7s} │ {'(t)':>8s} {'(t)':>8s} "
          f"{'(t)':>8s} {'(MN)':>7s} │ {'(GJ)':>7s} {'(TJ)':>7s} "
          f"{'(%)':>6s} │ {'S3':>5s} │ {'':>6s} {'(K)':>6s} "
          f"{'(t)':>8s} {'(%)':>6s}")
    print("─" * 120)
    
    for L in sled_lengths:
        for W in plate_widths:
            # ── Masses ──
            plate_mass_single = DENSITY * L * W * PLATE_THICKNESS * RIB_FACTOR
            plate_mass = plate_mass_single * N_PLATES
            struct_mass = STRUCTURE_MASS_PER_M * L
            sled_empty = plate_mass + struct_mass
            
            # ── How many LIM units does the sled span? ──
            plate_units = L / L_UNIT
            n_s3_active = plate_units * DUTY_S3
            
            # ── Thermal capacity ──
            thermal_cap = plate_mass * C_P * DELTA_T_MAX
            
            # ── What's the maximum spacecraft we can carry? ──
            # We need to find the max spacecraft mass where the design closes.
            # Start with a target and iterate.
            
            # First check: can we even close with minimum payload?
            sled_mass_min = sled_empty + SPACECRAFT_MASS_MIN
            F_req_min = sled_mass_min * ACCEL
            E_kin_min = 0.5 * sled_mass_min * V_TARGET**2
            slip_budget_min = thermal_cap / E_kin_min if E_kin_min > 0 else 0
            
            # Quick EM check at cruise speed (15 km/s)
            # What thrust can we get at reasonable slip?
            v_slip_test, F_test, can_thrust = find_slip_for_thrust(
                V_TARGET * 0.9, T_INITIAL, W, L, F_req_min)
            
            if not can_thrust:
                # Can't even thrust for minimum payload
                results.append({
                    'L': L, 'W': W, 'plate_mass': plate_mass,
                    'struct_mass': struct_mass, 'sled_empty': sled_empty,
                    'closes': False, 'reason': 'thrust',
                    'n_s3': n_s3_active, 'thermal_cap': thermal_cap,
                    'spacecraft_mass': 0, 'payload_frac': 0,
                })
                print(f"{L:>7.0f} {W:>7.2f} │ {plate_mass/1000:>8.0f} {struct_mass/1000:>8.0f} "
                      f"{sled_empty/1000:>8.0f} {F_req_min/1e6:>7.1f} │ {thermal_cap/1e9:>7.0f} "
                      f"{E_kin_min/1e12:>7.1f} {slip_budget_min*100:>6.3f} │ {n_s3_active:>5.1f} │ "
                      f"{'NO-F':>6s} {'---':>6s} {'---':>8s} {'---':>6s}")
                continue
            
            # Binary search for maximum spacecraft mass
            m_lo = SPACECRAFT_MASS_MIN
            m_hi = sled_empty * 10  # generous upper bound
            m_hi = min(m_hi, 50_000_000)  # 50,000 t cap
            
            best_spacecraft = 0
            best_T_final = T_INITIAL
            best_details = {}
            
            for _ in range(30):
                m_mid = (m_lo + m_hi) / 2
                sled_mass = sled_empty + m_mid
                
                success, T_final, E_heat, max_slip, details = simulate_launch_thermal(
                    sled_mass, W, L, plate_mass)
                
                if success and T_final < T_MAX:
                    best_spacecraft = m_mid
                    best_T_final = T_final
                    best_details = details
                    m_lo = m_mid
                else:
                    m_hi = m_mid
                
                if abs(m_hi - m_lo) < 1000:  # 1 tonne precision
                    break
            
            sled_total = sled_empty + best_spacecraft
            F_req = sled_total * ACCEL
            E_kin = 0.5 * sled_total * V_TARGET**2
            slip_budget = thermal_cap / E_kin if E_kin > 0 else 0
            payload_frac = best_spacecraft / sled_total if sled_total > 0 else 0
            
            closes = best_spacecraft >= SPACECRAFT_MASS_MIN
            
            results.append({
                'L': L, 'W': W, 'plate_mass': plate_mass,
                'struct_mass': struct_mass, 'sled_empty': sled_empty,
                'spacecraft_mass': best_spacecraft,
                'sled_total': sled_total, 'F_req': F_req,
                'thermal_cap': thermal_cap, 'E_kin': E_kin,
                'slip_budget': slip_budget,
                'n_s3': n_s3_active,
                'closes': closes, 'T_final': best_T_final,
                'payload_frac': payload_frac,
                'details': best_details,
            })
            
            status = "YES" if closes else "NO-T"
            print(f"{L:>7.0f} {W:>7.2f} │ {plate_mass/1000:>8.0f} {struct_mass/1000:>8.0f} "
                  f"{sled_total/1000:>8.0f} {F_req/1e6:>7.1f} │ {thermal_cap/1e9:>7.0f} "
                  f"{E_kin/1e12:>7.1f} {slip_budget*100:>6.3f} │ {n_s3_active:>5.1f} │ "
                  f"{status:>6s} {best_T_final:>6.0f} "
                  f"{best_spacecraft/1000:>8.0f} {payload_frac*100:>6.1f}")
        
        print()  # blank line between sled lengths
    
    # ── Summary: best designs ────────────────────────────────────────────────
    valid = [r for r in results if r['closes']]
    
    if valid:
        print("\n" + "=" * 90)
        print("TOP DESIGNS BY PAYLOAD FRACTION")
        print("=" * 90)
        valid.sort(key=lambda r: r['payload_frac'], reverse=True)
        
        for i, r in enumerate(valid[:15]):
            print(f"  {i+1:>2d}. L={r['L']:>5.0f}m  W={r['W']:.2f}m  "
                  f"plate={r['plate_mass']/1000:>6.0f}t  "
                  f"payload={r['spacecraft_mass']/1000:>7.0f}t  "
                  f"total={r['sled_total']/1000:>7.0f}t  "
                  f"η={r['payload_frac']*100:.1f}%  "
                  f"T={r['T_final']:.0f}K  "
                  f"#S3={r['n_s3']:.1f}")
        
        print("\n" + "=" * 90)
        print("TOP DESIGNS BY ABSOLUTE PAYLOAD MASS")
        print("=" * 90)
        valid.sort(key=lambda r: r['spacecraft_mass'], reverse=True)
        
        for i, r in enumerate(valid[:15]):
            print(f"  {i+1:>2d}. L={r['L']:>5.0f}m  W={r['W']:.2f}m  "
                  f"plate={r['plate_mass']/1000:>6.0f}t  "
                  f"payload={r['spacecraft_mass']/1000:>7.0f}t  "
                  f"total={r['sled_total']/1000:>7.0f}t  "
                  f"η={r['payload_frac']*100:.1f}%  "
                  f"T={r['T_final']:.0f}K  "
                  f"#S3={r['n_s3']:.1f}")
    else:
        print("\n*** NO VALID DESIGNS FOUND ***")
        print("All combinations either fail thrust or thermal limits.")
    
    # ── Detailed analysis of a few key points ────────────────────────────────
    print("\n" + "=" * 90)
    print("EM COUPLING ANALYSIS AT 15 km/s")
    print("(Stage 3, τ_p = 200m, T = 293K)")
    print("=" * 90)
    
    for W in [0.5, 1.0, 1.5, 2.0]:
        for L in [1200, 2500, 5000]:
            plate_units = L / L_UNIT
            n_s3 = plate_units * DUTY_S3
            a_lim = S3_L_ACTIVE * W * N_PLATES
            F_max = B_PLATE**2 * a_lim / (2 * MU0)
            
            # At various slip velocities
            print(f"\n  W={W:.1f}m, L={L}m ({n_s3:.1f} S3 segments):")
            print(f"    F_max (one seg) = {F_max/1e6:.2f} MN")
            for vs in [1, 5, 10, 20, 50, 100]:
                delta = skin_depth(S3_TAU_P, vs, 293)
                G = goodness_factor(S3_TAU_P, vs, 293, PLATE_THICKNESS)
                s = vs / (15000 + vs)
                eta = slip_efficiency(s, G)
                F_seg = F_max * eta
                F_total = F_seg * n_s3
                print(f"    v_slip={vs:>3d} m/s: δ={delta*1000:>6.0f}mm  "
                      f"t/δ={PLATE_THICKNESS/delta:.3f}  G={G:>6.1f}  "
                      f"sG={s*G:.3f}  η={eta*100:>5.1f}%  "
                      f"F_seg={F_seg/1e6:.3f}  F_tot={F_total/1e6:.2f} MN")
    
    return results


# ── Plotting ─────────────────────────────────────────────────────────────────

def make_plots(results):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not available for plotting")
        return
    
    os.makedirs('./graphs', exist_ok=True)
    
    valid = [r for r in results if r['closes']]
    if not valid:
        print("No valid designs to plot")
        return
    
    # Extract unique sweep values
    lengths = sorted(set(r['L'] for r in results))
    widths = sorted(set(r['W'] for r in results))
    
    # ── Plot 1: Payload fraction heatmap ─────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    pf_grid = np.zeros((len(widths), len(lengths)))
    for r in results:
        i = widths.index(r['W'])
        j = lengths.index(r['L'])
        pf_grid[i, j] = r['payload_frac'] * 100 if r['closes'] else 0
    
    im = ax.imshow(pf_grid, aspect='auto', origin='lower', cmap='RdYlGn',
                    extent=[lengths[0]/1000, lengths[-1]/1000,
                            widths[0], widths[-1]],
                    vmin=0, vmax=max(pf_grid.flatten()) * 1.1)
    
    # Annotate cells
    for r in results:
        j = lengths.index(r['L'])
        i = widths.index(r['W'])
        x = r['L'] / 1000
        y = r['W']
        if r['closes']:
            ax.text(x, y, f"{r['payload_frac']*100:.0f}%\n{r['spacecraft_mass']/1000:.0f}t",
                    ha='center', va='center', fontsize=7, fontweight='bold')
        else:
            ax.text(x, y, "✗", ha='center', va='center', fontsize=12, color='red')
    
    ax.set_xlabel('Sled Length (km)', fontsize=12)
    ax.set_ylabel('Plate Width (m)', fontsize=12)
    ax.set_title('Mass Driver Trade Study: Payload Fraction & Mass\n'
                 f'200mm γ-TiAl, B={B_PLATE}T, {V_TARGET/1000:.0f} km/s @ {ACCEL/G0:.1f}g',
                 fontsize=13)
    fig.colorbar(im, ax=ax, label='Payload Fraction (%)')
    fig.tight_layout()
    fig.savefig('./graphs/trade_payload_heatmap.png', dpi=300)
    plt.close(fig)
    
    # ── Plot 2: Payload mass vs sled length, lines for each width ────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(widths)))
    
    for k, W in enumerate(widths):
        subset = [r for r in results if r['W'] == W]
        Ls = [r['L']/1000 for r in subset]
        pay = [r['spacecraft_mass']/1000 if r['closes'] else 0 for r in subset]
        pf = [r['payload_frac']*100 if r['closes'] else 0 for r in subset]
        
        ax1.plot(Ls, pay, 'o-', color=colors[k], label=f'W={W:.2f}m', lw=1.5, ms=5)
        ax2.plot(Ls, pf, 's-', color=colors[k], label=f'W={W:.2f}m', lw=1.5, ms=5)
    
    ax1.set_ylabel('Max Spacecraft Mass (tonnes)', fontsize=12)
    ax1.set_title('Payload Capacity vs Sled Length', fontsize=13)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Sled Length (km)', fontsize=12)
    ax2.set_ylabel('Payload Fraction (%)', fontsize=12)
    ax2.set_title('Payload Fraction vs Sled Length', fontsize=13)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    fig.tight_layout()
    fig.savefig('./graphs/trade_payload_vs_length.png', dpi=300)
    plt.close(fig)
    
    # ── Plot 3: Plate mass vs payload mass (dead weight ratio) ───────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for k, W in enumerate(widths):
        subset = [r for r in results if r['W'] == W and r['closes']]
        if not subset:
            continue
        pm = [r['plate_mass']/1000 for r in subset]
        pay = [r['spacecraft_mass']/1000 for r in subset]
        ax.plot(pm, pay, 'o-', color=colors[k], label=f'W={W:.2f}m', lw=1.5, ms=5)
        for r in subset:
            ax.annotate(f"{r['L']/1000:.1f}km", 
                       (r['plate_mass']/1000, r['spacecraft_mass']/1000),
                       fontsize=6, alpha=0.7)
    
    # 1:1 line
    mx = max(r['plate_mass']/1000 for r in valid)
    ax.plot([0, mx*1.2], [0, mx*1.2], 'k--', alpha=0.3, label='1:1')
    ax.set_xlabel('Plate Mass (tonnes)', fontsize=12)
    ax.set_ylabel('Max Spacecraft Mass (tonnes)', fontsize=12)
    ax.set_title('Payload vs Plate Dead Weight', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('./graphs/trade_payload_vs_plate.png', dpi=300)
    plt.close(fig)
    
    # ── Plot 4: Thermal margin ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for k, W in enumerate(widths):
        subset = [r for r in results if r['W'] == W and r['closes']]
        if not subset:
            continue
        Ls = [r['L']/1000 for r in subset]
        margin = [T_MAX - r['T_final'] for r in subset]
        ax.plot(Ls, margin, 'o-', color=colors[k], label=f'W={W:.2f}m', lw=1.5, ms=5)
    
    ax.axhline(0, color='r', ls='--', alpha=0.5, label='Thermal limit')
    ax.set_xlabel('Sled Length (km)', fontsize=12)
    ax.set_ylabel('Thermal Margin (K)', fontsize=12)
    ax.set_title('Thermal Margin at Maximum Payload', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('./graphs/trade_thermal_margin.png', dpi=300)
    plt.close(fig)
    
    print("Plots saved to ./graphs/")


if __name__ == "__main__":
    results = run_trade_study()
    make_plots(results)
