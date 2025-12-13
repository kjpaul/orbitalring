#!/usr/bin/env python3
"""
Orbital Ring Deployment Simulation - Version 2

Key fixes from v1:
- Renamed t_plate → d_plate (thickness) to avoid confusion with T_plate (temperature)
- Fixed cryo power calculation (was dividing by COP, giving 20× inflation)
- Added debug output for diagnosing controller issues
- Added sanity checks for unusual physics values
"""
from __future__ import annotations
import argparse
import sys
import os
import math
import json
from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, List, Any, Tuple

try:
    import yaml
    HAVE_YAML = True
except Exception:
    HAVE_YAML = False

try:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
MU0 = 4 * math.pi * 1e-7          # Permeability of free space (H/m)
STEFAN_BOLTZMANN = 5.67e-8        # Stefan-Boltzmann constant (W/(m²·K⁴))
SEC_PER_MONTH = 30.0 * 86400.0    # Seconds per 30-day month

# =============================================================================
# ORBITAL PARAMETERS (250 km altitude)
# =============================================================================
R_EARTH = 6_378_137.0
ALT = 250_000.0
R_ORBIT = R_EARTH + ALT
G_LOCAL = 9.073
L_ORBIT = 2 * math.pi * R_ORBIT
V_ORBIT = 7_754.866
V_CASING_FINAL = 483.331
GRAPH_DIR = "graphs"

# =============================================================================
# CONFIGURATION DATACLASSES
# =============================================================================
@dataclass
class LIM:
    """Linear Induction Motor configuration parameters."""
    # Coil geometry
    n_turns: int = 100
    tau_p: float = 50.0             # Pole pitch (m)
    w_coil: float = 2.0             # Coil width (m)
    gap: float = 0.05               # Air gap (m)
    pitch_count: int = 3
    spacing: float = 250.0          # LIM site spacing (m)
    
    # Electrical
    i_peak_target: float = 650.0
    i_peak_min: float = 1.0
    ic: float = 800.0
    volts_max_user: float = 200000.0
    
    # Insulation
    d_kapton_mm: float = 1.0
    e_allow_kv_per_mm: float = 10.0
    k_fill: float = 0.449
    w_tape: float = 0.012
    
    # Slip control
    v_slip_min: float = 1.0
    v_slip_max: float = 200.0
    slip_min: float = 0.01
    v_rel_min: float = 1.0
    
    # Reaction plate geometry (d = thickness, not temperature!)
    d_plate: float = 0.080          # Plate THICKNESS (m) - renamed from t_plate
    heat_sink_l: float = 30.0
    casing_outer_w: float = 10.0
    
    # Reaction plate temperature
    T_plate_min: float = 100.0      # Minimum plate TEMPERATURE (K) - renamed from t_plate_min
    
    # Aluminum properties
    rho_alu_e_20C: float = 2.65e-8  # Resistivity at 20°C (Ω·m)
    alpha_alu_e: float = 0.00366    # Temperature coefficient (1/K)
    
    # Thermal
    em_alu: float = 0.85
    em_heat_sink: float = 0.90
    
    # External heat
    q_sun_1au: float = 1361.0
    q_earth_day: float = 650.0
    q_shielding: float = 0.005
    
    # Cryogenic - NOTE: This is a multiplier, not COP!
    # P_cryo = cryo_eff * heat_load (like original code)
    # For COP interpretation, set cryo_cop and use that instead
    cryo_eff: float = 0.05          # Legacy: multiplier for cryo overhead
    
    # Efficiencies
    inv_eff: float = 0.90
    lim_eff: float = 0.95
    
    # Power limit
    max_site_power: float = 8.0e6
    
    # Controller
    di_max_per_s: float = 1.0
    dvslip_max_per_s: float = 2.0
    
    # Hysteresis
    alpha_angle_deg: float = 20.0
    
    # Legacy
    thrust_coeff: float = 1.0


@dataclass
class Params:
    """Simulation parameters."""
    m_cable_per_m: float = 99_198.0
    m_casing_per_m: float = 12_000.0
    cable_area_m2: float = 56.9
    sigma_break: float = 25e9
    
    r_orbit: float = R_ORBIT
    l_orbit: float = L_ORBIT
    g_local: float = G_LOCAL
    v_orbit: float = V_ORBIT
    v_casing_final: float = V_CASING_FINAL
    
    lim: LIM = field(default_factory=LIM)
    
    dt: float = 1.0
    max_time: float = 2.0e8
    tol_v: float = 1e-3
    
    verbose: bool = False
    explain: bool = False
    debug: bool = False             # New: extra debug output
    log: str = "coarse"
    log_interval: float = 600.0
    plot_skip_s: float = 600.0
    plot: List[str] = field(default_factory=list)
    save_plots: bool = False
    outdir: str = "."
    graph_dir: str = GRAPH_DIR


# =============================================================================
# STATE AND OUTPUT DATACLASSES
# =============================================================================
@dataclass
class LIMState:
    """State variables that evolve each timestep."""
    i_peak: float
    v_slip: float
    T_plate: float      # Temperature (K), not thickness!


@dataclass
class LIMOutputs:
    """Computed outputs for current state."""
    F_site: float
    P_thrust: float
    P_eddy: float
    P_hyst: float
    P_cryo: float
    P_site: float
    V_site: float
    f_supply: float
    f_slip: float
    B_plate: float
    rho_alu: float
    delta_skin: float
    G: float
    eta_slip: float


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def ring_mass(per_m: float, p: Params) -> float:
    return per_m * p.l_orbit

def lim_sites(p: Params) -> int:
    return int(round(p.l_orbit / p.lim.spacing))

def l_active(p: Params) -> float:
    return p.lim.tau_p * p.lim.pitch_count

def a_coil(p: Params) -> float:
    return p.lim.tau_p * p.lim.w_coil

def kapton_safe_v(p: Params) -> float:
    e_allow = p.lim.e_allow_kv_per_mm * 1e6
    d_ins = p.lim.d_kapton_mm * 1e-3
    return e_allow * d_ins

def volts_max_allowed(p: Params) -> float:
    kapton_v = max(0, p.lim.n_turns - 1) * kapton_safe_v(p)
    return min(p.lim.volts_max_user, kapton_v if kapton_v > 0 else p.lim.volts_max_user)


# =============================================================================
# PHYSICS FUNCTIONS
# =============================================================================
def rho_alu_at_temp(T_plate: float, p: Params) -> float:
    """Aluminum resistivity at temperature T_plate (Ω·m)."""
    T_20C = 293.0
    rho = p.lim.rho_alu_e_20C * (1.0 + p.lim.alpha_alu_e * (T_plate - T_20C))
    return max(rho, 1.0e-10)


def b_plate_peak(i_peak: float, p: Params) -> float:
    """Magnetic field at reaction plate (T)."""
    W = p.lim.w_coil
    g = p.lim.gap
    return (2 * MU0 * p.lim.n_turns * i_peak / (math.pi * W)) * math.atan(W / (2 * g))


def b_coil_peak(i_peak: float, p: Params) -> float:
    """Peak field in coil for hysteresis calc (T)."""
    return MU0 * p.lim.n_turns * i_peak / (math.pi * (p.lim.w_coil + p.lim.tau_p))


def f_from_v(v: float, p: Params) -> float:
    """Frequency from velocity: f = v / (2τp)."""
    return max(v, 1e-9) / (2.0 * p.lim.tau_p)


def skin_depth(f_slip: float, rho_alu: float) -> float:
    """Skin depth in aluminum (m)."""
    if f_slip <= 0:
        return 1.0
    return math.sqrt(rho_alu / (math.pi * MU0 * f_slip))


def slip_ratio(v_slip: float, v_rel: float, p: Params) -> float:
    """Slip ratio s = f_slip / f_supply."""
    f_slip = f_from_v(v_slip, p)
    f_supply = f_from_v(v_rel + v_slip, p)
    return f_slip / max(f_supply, 1e-9)


def goodness_factor(f_slip: float, rho_alu: float, delta_eff: float, p: Params) -> float:
    """LIM goodness factor G."""
    if f_slip <= 0:
        return 1.0
    omega_slip = 2.0 * math.pi * f_slip
    sigma_alu = 1.0 / rho_alu
    return (omega_slip * MU0 * sigma_alu * delta_eff * p.lim.tau_p) / math.pi


def reaction_plate_temperature(P_eddy: float, p: Params) -> float:
    """Equilibrium plate temperature from radiation (K)."""
    A_active = 2.0 * l_active(p) * p.lim.w_coil
    A_heatsink = 2.0 * p.lim.heat_sink_l * p.lim.w_coil
    A_total = A_active + A_heatsink
    
    if P_eddy <= 0 or A_total <= 0:
        return p.lim.T_plate_min
    
    T_eq = (P_eddy / (p.lim.em_alu * STEFAN_BOLTZMANN * A_total)) ** 0.25
    return max(T_eq, p.lim.T_plate_min)


# =============================================================================
# MAIN LIM CALCULATION
# =============================================================================
def compute_lim_outputs(state: LIMState, v_rel: float, p: Params) -> LIMOutputs:
    """Compute all LIM outputs from current state."""
    
    # Temperature-dependent resistivity
    rho_alu = rho_alu_at_temp(state.T_plate, p)
    
    # Frequencies
    f_slip = f_from_v(state.v_slip, p)
    f_supply = f_from_v(v_rel + state.v_slip, p)
    
    # Skin depth - use d_plate (thickness), not T_plate (temperature)!
    delta_raw = skin_depth(f_slip, rho_alu)
    delta_eff = min(p.lim.d_plate, delta_raw)
    
    # Magnetic field
    B = b_plate_peak(state.i_peak, p)
    
    # Goodness factor and efficiency
    G = goodness_factor(f_slip, rho_alu, delta_eff, p)
    s = slip_ratio(state.v_slip, v_rel, p)
    
    # Slip efficiency: η = 2sG / (1 + s²G²)
    denom = 1.0 + s * s * G * G
    eta_slip = (2.0 * s * G) / denom if denom > 0 else 0.0
    
    # Clamp eta_slip to physical range [0, 1]
    eta_slip = max(0.0, min(1.0, eta_slip))
    
    # Thrust
    A_active = p.lim.w_coil * l_active(p)
    F_max = (B * B / (2.0 * MU0)) * A_active
    F_site = F_max * eta_slip * p.lim.thrust_coeff
    
    # Thrust power
    v_rel_eff = max(v_rel, 0.1)
    P_thrust = F_site * v_rel_eff
    
    # Eddy losses (in plate)
    V_eddy = p.lim.w_coil * l_active(p) * delta_eff
    P_eddy = (math.pi**2 / (6.0 * rho_alu)) * (B * B) * (delta_eff**2) * (f_slip**2) * V_eddy
    
    # Hysteresis losses (in HTS)
    Bc = b_coil_peak(state.i_peak, p)
    q = Bc * p.lim.ic * p.lim.w_tape * math.sin(math.radians(p.lim.alpha_angle_deg))
    L_HTS_coil = 2.0 * (p.lim.w_coil + p.lim.tau_p) * p.lim.n_turns
    P_hyst = q * L_HTS_coil * f_supply * 3 * p.lim.pitch_count
    
    # External heat absorption
    Q_abs_m = (p.lim.q_sun_1au + p.lim.q_earth_day) * p.lim.casing_outer_w * p.lim.q_shielding
    Q_abs_site = Q_abs_m * p.lim.spacing
    
    # Cryo power - FIXED: use multiplier like original code, not division by COP
    # P_cryo represents the overhead of running the cryogenic system
    Q_cryo_load = P_hyst + Q_abs_site  # Heat that cryo must remove (eddy radiates away)
    P_cryo = p.lim.cryo_eff * Q_cryo_load  # Small fraction, not 20× multiplication!
    
    # Total site power
    P_site = (P_thrust + P_eddy + P_hyst + P_cryo) / (p.lim.inv_eff * p.lim.lim_eff)
    
    # Voltage estimate
    omega = 2.0 * math.pi * f_supply
    A_link = a_coil(p)
    psi_peak = p.lim.n_turns * B * A_link * p.lim.k_fill
    V_phase_peak = omega * psi_peak
    V_site = V_phase_peak / math.sqrt(2.0)
    
    return LIMOutputs(
        F_site=F_site, P_thrust=P_thrust, P_eddy=P_eddy, P_hyst=P_hyst,
        P_cryo=P_cryo, P_site=P_site, V_site=V_site,
        f_supply=f_supply, f_slip=f_slip, B_plate=B,
        rho_alu=rho_alu, delta_skin=delta_eff, G=G, eta_slip=eta_slip
    )


# =============================================================================
# CONTROLLER
# =============================================================================
def find_optimal_slip(state: LIMState, v_rel: float, p: Params, P_target: float) -> float:
    """Find slip velocity that gives approximately P_target power.
    
    Uses coarse-to-fine search since P(v_slip) is not monotonic.
    """
    # Coarse search
    best_slip = state.v_slip
    best_error = float('inf')
    
    # First pass: coarse grid
    for v_slip_try in [1, 2, 5, 10, 20, 50, 100, 150, 200]:
        if v_slip_try < p.lim.v_slip_min or v_slip_try > p.lim.v_slip_max:
            continue
        test_state = LIMState(state.i_peak, float(v_slip_try), state.T_plate)
        test_out = compute_lim_outputs(test_state, v_rel, p)
        error = abs(test_out.P_site - P_target)
        if error < best_error:
            best_error = error
            best_slip = float(v_slip_try)
    
    # Second pass: fine search around best
    for offset in [-20, -10, -5, 5, 10, 20]:
        v_slip_try = best_slip + offset
        if v_slip_try < p.lim.v_slip_min or v_slip_try > p.lim.v_slip_max:
            continue
        test_state = LIMState(state.i_peak, float(v_slip_try), state.T_plate)
        test_out = compute_lim_outputs(test_state, v_rel, p)
        error = abs(test_out.P_site - P_target)
        if error < best_error:
            best_error = error
            best_slip = float(v_slip_try)
    
    return best_slip


def update_controller(state: LIMState, outputs: LIMOutputs, v_rel: float, p: Params, debug: bool = False) -> LIMState:
    """Controller that finds optimal slip for target power."""
    V_cap = volts_max_allowed(p)
    P_target = p.lim.max_site_power * 0.95  # Target 95% of limit
    
    power_ratio = outputs.P_site / p.lim.max_site_power if p.lim.max_site_power > 0 else 0
    voltage_ratio = outputs.V_site / V_cap if V_cap > 0 else 0
    limit_ratio = max(power_ratio, voltage_ratio)
    
    if debug:
        which_limit = "POWER" if power_ratio >= voltage_ratio else "VOLTAGE"
        print(f"  DEBUG: limit_ratio={limit_ratio:.3f} ({which_limit}) "
              f"P={outputs.P_site/1e6:.2f}MW V={outputs.V_site/1e3:.1f}kV "
              f"I={state.i_peak:.0f}A slip={state.v_slip:.1f}m/s "
              f"η={outputs.eta_slip:.4f} G={outputs.G:.1f} F={outputs.F_site:.0f}N")
    
    # Current control: ramp toward target, back off if over limit
    if limit_ratio > 1.05:
        # Over limit - reduce current proportionally
        scale = 0.95 / limit_ratio
        want_i = max(p.lim.i_peak_min, state.i_peak * scale)
    elif limit_ratio > 0.95:
        # At limit - hold current
        want_i = state.i_peak
    else:
        # Under limit - ramp current toward target
        want_i = p.lim.i_peak_target
    
    # Apply current ramp rate limit
    di = want_i - state.i_peak
    new_i = state.i_peak + math.copysign(min(abs(di), p.lim.di_max_per_s * p.dt), di)
    new_i = min(max(new_i, p.lim.i_peak_min), p.lim.i_peak_target)
    
    # Slip control: find optimal slip for current power target
    # But only adjust slip slowly to avoid oscillations
    if limit_ratio > 1.05:
        # Over limit - find slip that gives target power
        optimal_slip = find_optimal_slip(state, v_rel, p, P_target)
        want_vslip = optimal_slip
    elif limit_ratio < 0.5:
        # Well under limit - ramp toward max slip
        want_vslip = p.lim.v_slip_max
    else:
        # Near limit - hold slip steady
        want_vslip = state.v_slip
    
    # Apply slip ramp rate limit
    dv = want_vslip - state.v_slip
    new_vslip = state.v_slip + math.copysign(min(abs(dv), p.lim.dvslip_max_per_s * p.dt), dv)
    new_vslip = min(max(new_vslip, p.lim.v_slip_min), p.lim.v_slip_max)
    
    if debug:
        print(f"         -> want_i={want_i:.0f} want_slip={want_vslip:.0f} "
              f"-> new_i={new_i:.0f} new_slip={new_vslip:.1f}")
    
    return LIMState(i_peak=new_i, v_slip=new_vslip, T_plate=state.T_plate)


# =============================================================================
# YAML LOADING
# =============================================================================
def _coerce_like(cur: Any, val: Any) -> Any:
    if val is None:
        return cur
    if isinstance(cur, list):
        if isinstance(val, list):
            return list(val)
        if isinstance(val, str):
            return [s.strip() for s in val.split(",") if s.strip()]
        return [val]
    if isinstance(cur, bool):
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return bool(val)
        if isinstance(val, str):
            return val.strip().lower() in ("true", "yes", "on", "1")
        return bool(val)
    if isinstance(cur, float):
        if isinstance(val, (int, float)):
            return float(val)
        if isinstance(val, str):
            return float(val.strip())
    if isinstance(cur, int):
        if isinstance(val, (int, float)):
            return int(val)
        if isinstance(val, str):
            return int(float(val.strip()))
    return val


def load_yaml(path: str, base: Params) -> Params:
    if not HAVE_YAML:
        raise RuntimeError("PyYAML not installed.")
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    
    # Handle legacy parameter names
    if 'lim' in data and isinstance(data['lim'], dict):
        lim_data = data['lim']
        # Map old names to new names
        if 't_plate' in lim_data and 'd_plate' not in lim_data:
            lim_data['d_plate'] = lim_data.pop('t_plate')
        if 't_plate_min' in lim_data and 'T_plate_min' not in lim_data:
            lim_data['T_plate_min'] = lim_data.pop('t_plate_min')
    
    for k, v in data.items():
        if hasattr(base, k):
            cur = getattr(base, k)
            if isinstance(cur, LIM) and isinstance(v, dict):
                lim = cur
                for kk, vv in v.items():
                    if hasattr(lim, kk):
                        setattr(lim, kk, _coerce_like(getattr(lim, kk), vv))
                setattr(base, k, lim)
            else:
                setattr(base, k, _coerce_like(cur, v))
        else:
            print(f"Warning: unknown config key '{k}' ignored.", file=sys.stderr)
    return base


# =============================================================================
# OUTPUT
# =============================================================================
def banner(p: Params):
    print("=" * 79)
    print("Orbital Ring Deployment Simulation v2")
    print("=" * 79)
    print(f"v_orbit={p.v_orbit:.3f} m/s  v_casing_final={p.v_casing_final:.3f} m/s")
    Mc, Ms = ring_mass(p.m_cable_per_m, p), ring_mass(p.m_casing_per_m, p)
    print(f"Mass: cable={Mc:.3e} kg  casing={Ms:.3e} kg")
    print(f"LIM: τp={p.lim.tau_p:.0f}m × w={p.lim.w_coil:.1f}m × {p.lim.pitch_count} = {l_active(p):.0f}m")
    print(f"     N={p.lim.n_turns} gap={p.lim.gap*100:.1f}cm d_plate={p.lim.d_plate*1000:.0f}mm")
    print(f"     Sites={lim_sites(p)} spacing={p.lim.spacing:.0f}m")
    print(f"     I_target={p.lim.i_peak_target:.0f}A v_slip=[{p.lim.v_slip_min:.0f},{p.lim.v_slip_max:.0f}]m/s")
    print(f"     V_max={volts_max_allowed(p)/1000:.1f}kV P_max={p.lim.max_site_power/1e6:.1f}MW")
    print()


def explain_text():
    print("Physics Model (v2 - fixed):")
    print("  B_plate = (2μ₀NI)/(πW) × arctan(W/(2g))")
    print("  δ = √(ρ_alu(T) / (πμ₀f_slip))")
    print("  G = (ω_slip × μ₀ × σ_alu × δ_eff × τp) / π")
    print("  η_slip = 2sG / (1 + s²G²)")
    print("  F = (B²/2μ₀) × A_active × η_slip")
    print("  P_cryo = cryo_eff × Q_load  (NOT Q_load/cryo_eff)")
    print()


def next_log_target(t: float) -> float:
    if t < 60: return 60
    if t < 3600: return 3600
    if t < 86400: return 86400
    # Log every week up to 2 months
    week = 604800.0
    if t < 8 * week:
        return (math.floor(t / week) + 1) * week
    month = SEC_PER_MONTH
    return (math.floor(t / month) + 1) * month


# =============================================================================
# PLOTTING
# =============================================================================
def downsample(x: List[float], y: List[float], max_pts: int = 5000):
    n = len(x)
    if n <= max_pts:
        return x, y
    step = max(1, n // max_pts)
    return x[::step], y[::step]


def annotate_last(ax, x, y, unit="", fmt=".0f", dx=-40, dy=-40, fontsize=11):
    if not x or not y:
        return
    val = f"{y[-1]:{fmt}} {unit}".strip()
    ax.annotate(val, xy=(x[-1], y[-1]), xytext=(dx, dy), textcoords="offset points",
                fontsize=fontsize, ha="left", va="bottom",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0))


def maybe_plot(p: Params, history: Dict[str, List[float]]):
    if not p.plot or not HAVE_MPL:
        return
    
    graph_path = f"{p.outdir.rstrip('/')}/{p.graph_dir}"
    os.makedirs(graph_path, exist_ok=True)
    
    ts = history['t']
    i0 = 0
    while i0 < len(ts) and ts[i0] < p.plot_skip_s:
        i0 += 1
    
    t_months = [t / SEC_PER_MONTH for t in ts[i0:]]
    
    def get_data(key, scale=1.0):
        if key not in history:
            return [], []
        data = [v * scale for v in history[key][i0:]]
        n = min(len(t_months), len(data))
        return downsample(t_months[:n], data[:n])
    
    def save_or_show(fig, name):
        if p.save_plots:
            plt.savefig(f"{graph_path}/{name}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()
    
    # Velocities
    if "velocities" in p.plot or "all" in p.plot:
        tm, vc = get_data('v_casing')
        _, vb = get_data('v_cable')
        fig, ax = plt.subplots()
        ax.plot(tm, vc, label="v_casing")
        ax.plot(tm, vb, label="v_cable")
        annotate_last(ax, tm, vc, "m/s")
        annotate_last(ax, tm, vb, "m/s", dy=10)
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Velocity [m/s]")
        ax.legend()
        ax.set_title("Velocities")
        save_or_show(fig, "velocities")
    
    # Power
    if "power" in p.plot or "all" in p.plot:
        tm, pwr = get_data('P_site', 1e-6)
        fig, ax = plt.subplots()
        ax.plot(tm, pwr)
        ax.axhline(p.lim.max_site_power/1e6, color='r', ls='--', label='Limit')
        annotate_last(ax, tm, pwr, "MW", fmt=".2f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Power [MW]")
        ax.legend()
        ax.set_title("Site Power")
        save_or_show(fig, "power")
    
    # Voltage
    if "voltage" in p.plot or "all" in p.plot:
        tm, v = get_data('V_site', 1e-3)
        fig, ax = plt.subplots()
        ax.plot(tm, v)
        ax.axhline(volts_max_allowed(p)/1e3, color='r', ls='--', label='Limit')
        annotate_last(ax, tm, v, "kV", fmt=".1f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Voltage [kV]")
        ax.legend()
        ax.set_title("Phase Voltage")
        save_or_show(fig, "voltage")
    
    # Current and slip
    if "current_slip" in p.plot or "all" in p.plot:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        tm, curr = get_data('i_peak')
        _, vslip = get_data('v_slip')
        ax1.plot(tm, curr, 'b-', label='I_peak')
        ax2.plot(tm, vslip, 'r-', label='v_slip')
        ax1.set_xlabel("Time [months]")
        ax1.set_ylabel("Current [A]", color='b')
        ax2.set_ylabel("Slip velocity [m/s]", color='r')
        ax1.set_title("Controller Output")
        save_or_show(fig, "current_slip")
    
    # Plate temperature
    if "plate_temp" in p.plot or "all" in p.plot:
        tm, T = get_data('T_plate')
        fig, ax = plt.subplots()
        ax.plot(tm, T, color='orange')
        annotate_last(ax, tm, T, "K", fmt=".0f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Temperature [K]")
        ax.set_title("Reaction Plate Temperature")
        save_or_show(fig, "plate_temp")
    
    # Losses
    if "losses" in p.plot or "all" in p.plot:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        tm, p_eddy = get_data('P_eddy', 1e-3)
        _, p_hyst = get_data('P_hyst')
        ax1.plot(tm, p_eddy, 'r-', label='Eddy')
        ax2.plot(tm, p_hyst, 'b-', label='Hyst')
        ax1.set_xlabel("Time [months]")
        ax1.set_ylabel("Eddy [kW]", color='r')
        ax2.set_ylabel("Hysteresis [W]", color='b')
        ax1.set_title("Losses")
        save_or_show(fig, "losses")
    
    # Relative speed
    if "relative_speed" in p.plot or "all" in p.plot:
        tm, vrel = get_data('v_rel')
        fig, ax = plt.subplots()
        ax.plot(tm, vrel)
        annotate_last(ax, tm, vrel, "m/s", fmt=".0f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("v_rel [m/s]")
        ax.set_title("Relative Velocity")
        save_or_show(fig, "relative_speed")
    
    # Stress
    if "stress" in p.plot or "all" in p.plot:
        tm, tension = get_data('tension')
        stress = [T / p.cable_area_m2 / 1e9 for T in tension]
        tm, stress = downsample(tm, stress)
        fig, ax = plt.subplots()
        ax.plot(tm, stress)
        annotate_last(ax, tm, stress, "GPa", fmt=".2f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Stress [GPa]")
        ax.set_title("Cable Stress")
        save_or_show(fig, "stress")
    
    # Efficiency
    if "efficiency" in p.plot or "all" in p.plot:
        tm, eta = get_data('eta_slip')
        eta_pct = [e * 100 for e in eta]
        tm, eta_pct = downsample(tm, eta_pct)
        fig, ax = plt.subplots()
        ax.plot(tm, eta_pct)
        annotate_last(ax, tm, eta_pct, "%", fmt=".1f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Slip Efficiency [%]")
        ax.set_title("LIM Slip Efficiency")
        save_or_show(fig, "efficiency")
    
    # Goodness factor
    if "goodness" in p.plot or "all" in p.plot:
        tm, G = get_data('G')
        fig, ax = plt.subplots()
        ax.plot(tm, G)
        annotate_last(ax, tm, G, "", fmt=".1f")
        ax.set_xlabel("Time [months]")
        ax.set_ylabel("Goodness Factor G")
        ax.set_title("LIM Goodness Factor")
        save_or_show(fig, "goodness")


# =============================================================================
# MAIN SIMULATION
# =============================================================================
def simulate(p: Params) -> Dict:
    if p.verbose:
        banner(p)
        if p.explain:
            explain_text()
    
    Mcasing = ring_mass(p.m_casing_per_m, p)
    Mcable = ring_mass(p.m_cable_per_m, p)
    Mtotal = Mcasing + Mcable
    sites = lim_sites(p)
    
    v_casing = p.v_orbit
    v_cable = p.v_orbit
    
    state = LIMState(
        i_peak=p.lim.i_peak_min,
        v_slip=p.lim.v_slip_min,
        T_plate=p.lim.T_plate_min
    )
    
    history = {k: [] for k in [
        't', 'v_casing', 'v_cable', 'tension', 'net_per_m',
        'P_site', 'V_site', 'i_peak', 'v_slip', 'v_rel', 'a_out',
        'E_site', 'E_total', 'E_site_ke', 'E_total_ke',
        'P_eddy', 'P_hyst', 'T_plate', 'f_supply', 'f_slip',
        'rho_alu', 'delta_skin', 'G', 'eta_slip', 'F_site'
    ]}
    
    E_site = E_total = E_site_ke = E_total_ke = 0.0
    
    min_T = (float("inf"), None)
    max_T = (-float("inf"), None)
    
    next_log = next_log_target(0.0) if p.verbose and p.log == "coarse" else float("inf")
    
    t = 0.0
    steps = 0
    broke = False
    reason = ""
    
    while t < p.max_time and (v_casing - p.v_casing_final) > p.tol_v:
        v_rel = max(v_cable - v_casing, 0.0)
        
        # Compute outputs
        outputs = compute_lim_outputs(state, v_rel, p)
        
        # Apply thrust
        F_tot = outputs.F_site * sites * 2.0
        v_casing += (-F_tot / Mcasing) * p.dt
        v_cable += (+F_tot / Mcable) * p.dt
        
        # Orbital mechanics
        a_out_cable = (v_cable ** 2) / p.r_orbit
        a_out_casing = (v_casing ** 2) / p.r_orbit
        mpm = p.m_cable_per_m + p.m_casing_per_m
        a_out_weighted = (p.m_cable_per_m * a_out_cable + p.m_casing_per_m * a_out_casing) / mpm
        
        net_per_m = (
            p.m_cable_per_m * (a_out_cable - p.g_local) +
            p.m_casing_per_m * (a_out_casing - p.g_local)
        )
        
        T = max(0.0, net_per_m * p.r_orbit)
        sigma = T / p.cable_area_m2
        
        if sigma > p.sigma_break:
            broke = True
            reason = f"Cable failure: {sigma:.3e} > {p.sigma_break:.3e} Pa"
            if p.verbose:
                print(reason)
            break
        
        # Energy
        E_site += outputs.P_site * 2.0 * p.dt
        E_total += outputs.P_site * sites * 2.0 * p.dt
        E_site_ke += outputs.P_thrust * 2.0 * p.dt
        E_total_ke += outputs.P_thrust * sites * 2.0 * p.dt
        
        # Update controller
        state = update_controller(state, outputs, v_rel, p, debug=(p.debug and steps % 1000 == 0))
        
        # Update temperature
        new_outputs = compute_lim_outputs(state, v_rel, p)
        state = LIMState(state.i_peak, state.v_slip, 
                        reaction_plate_temperature(new_outputs.P_eddy, p))
        
        # Track extremes
        if T < min_T[0]:
            min_T = (T, t)
        if T > max_T[0]:
            max_T = (T, t)
        
        t += p.dt
        steps += 1
        
        # Record history
        if steps % 10 == 0:
            history['t'].append(t)
            history['v_casing'].append(v_casing)
            history['v_cable'].append(v_cable)
            history['tension'].append(T)
            history['net_per_m'].append(net_per_m)
            history['P_site'].append(outputs.P_site)
            history['V_site'].append(outputs.V_site)
            history['i_peak'].append(state.i_peak)
            history['v_slip'].append(state.v_slip)
            history['v_rel'].append(v_rel)
            history['a_out'].append(a_out_weighted)
            history['E_site'].append(E_site)
            history['E_total'].append(E_total)
            history['E_site_ke'].append(E_site_ke)
            history['E_total_ke'].append(E_total_ke)
            history['P_eddy'].append(outputs.P_eddy)
            history['P_hyst'].append(outputs.P_hyst)
            history['T_plate'].append(state.T_plate)
            history['f_supply'].append(outputs.f_supply)
            history['f_slip'].append(outputs.f_slip)
            history['rho_alu'].append(outputs.rho_alu)
            history['delta_skin'].append(outputs.delta_skin)
            history['G'].append(outputs.G)
            history['eta_slip'].append(outputs.eta_slip)
            history['F_site'].append(outputs.F_site)
        
        # Logging
        if p.verbose and t >= next_log:
            print(f"t={t:12.0f}s | {t/SEC_PER_MONTH:5.2f}mo | "
                  f"v_cas={v_casing:8.2f} v_cab={v_cable:8.2f} | "
                  f"I={state.i_peak:5.0f}A slip={state.v_slip:5.0f}m/s | "
                  f"P={outputs.P_site/1e6:5.2f}MW η={outputs.eta_slip:.3f} G={outputs.G:.0f}")
            # Always show controller state at log points
            V_cap = volts_max_allowed(p)
            power_ratio = outputs.P_site / p.lim.max_site_power
            voltage_ratio = outputs.V_site / V_cap
            print(f"         P_ratio={power_ratio:.3f} V_ratio={voltage_ratio:.3f} "
                  f"V_site={outputs.V_site/1e3:.1f}kV V_cap={V_cap/1e3:.0f}kV "
                  f"F={outputs.F_site:.0f}N")
            next_log = next_log_target(next_log)
    
    if p.verbose:
        print(f"Complete. {'FAILED: ' + reason if broke else ''}")
    
    out = {
        "results": {
            "months_elapsed": t / SEC_PER_MONTH,
            "v_casing_end": v_casing,
            "v_cable_end": v_cable,
            "failed": broke,
            "stress_max_GPa": max_T[0] / p.cable_area_m2 / 1e9,
            "total_energy_EJ": E_total / 1e18,
            "total_ke_EJ": E_total_ke / 1e18,
        }
    }
    print(json.dumps(out, indent=2))
    
    maybe_plot(p, history)
    return out


# =============================================================================
# CLI
# =============================================================================
def parse_args(argv=None) -> Params:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=str)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--explain", action="store_true")
    ap.add_argument("--plot", type=str)
    ap.add_argument("--save-plots", action="store_true")
    ap.add_argument("--outdir", type=str, default=".")
    ap.add_argument("--plot-skip-s", type=float)
    
    args = ap.parse_args(argv)
    p = Params()
    p.verbose = args.verbose
    p.debug = args.debug
    p.explain = args.explain
    p.save_plots = args.save_plots
    p.outdir = args.outdir
    
    if args.plot:
        p.plot = [s.strip() for s in args.plot.split(",")]
    if args.plot_skip_s is not None:
        p.plot_skip_s = args.plot_skip_s
    if args.config:
        p = load_yaml(args.config, p)
    
    return p


def main(argv=None) -> int:
    p = parse_args(argv)
    simulate(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
