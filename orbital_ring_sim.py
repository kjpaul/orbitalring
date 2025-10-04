#!/usr/bin/env python3
"""
orbital_ring_sim.py — Open-source orbital ring deployment simulator (CLI)

Purpose
-------
This tool models the *tangential momentum exchange* between an outer casing and an inner cable
during an orbital-ring deployment phase. Linear induction motors (LIMs) couple the two,
decelerating the casing while accelerating the cable. We track velocities over time and stop
when the casing reaches ground-synchronous speed.

This code is designed for clarity, reproducibility, and learning. It prints what it is doing,
explains the governing equations (with --explain), and can optionally plot selected outputs
(--plot velocities[,tension,weight]) to visualize the evolution.

Key Physical Setup (defaults reflect a circular orbit ~250 km above the equator)
--------------------------------------------------------------------------------
- Orbital altitude (above mean equator):       h = 250,000 m
- Orbital radius:                              r = R_EARTH + h  (R_Earth ~ 6,378,137 m -> r ~ 6,628,137 m)
- Circular-orbit tangential speed:             v_orbit = 7,754.866 m/s    (given)
- Final casing speed (ground synchronous):     v_casing_final = 483.331 m/s (given at equator, same altitude frame)
- Local g at altitude:                         g = 9.073 m/s^2 (given; for future tension/weight modeling)
- Ring circumference (circular):               L = 2π r  (~ 41,645,813 m for r ≈ 6,628,137 m)

Two fixed parameters in this script (by default) are v_orbit and v_casing_final.
Both are set by orbit altitude (assuming a circular orbit at h = 250 km).

What this simulator *does*
--------------------------
- Runs a time-dependent loop that integrates the coupled ODE system for casing and cable velocities:
  d(v_casing)/dt =  F_casing / M_casing
  d(v_cable)/dt  =  F_cable  / M_cable
  with F_cable = -F_casing  (Newton's third law; LIM force is internal to the ring system).
- You can choose a thrust model:
  * constant_force: apply a fixed tangential force magnitude |F| until target casing speed is reached
  * constant_power: apply a fixed power P; the instantaneous force is F = P / v_rel, where v_rel is the
                    relative speed of the actuator. In this simple model we use casing speed for decel
                    and cable speed for accel sides respectively.
- Stops when v_casing <= v_casing_final (to within tolerance).

What this simulator *does not yet do*
-------------------------------------
- Full structural tension bookkeeping and net weight (radial force balance) over time.
  Placeholders are provided to compute first-order estimates; see TODO notes in code.
- Spatially varying properties around the ring, losses, thermal/EM detail, LIM slip, etc.

Equations (concise)
-------------------
1) Momentum exchange (internal LIM force only; no external torque):
   F_cable = -F_casing  and  d/dt(M_casing v_casing + M_cable v_cable) = 0  (≈ conserved)
2) Kinematics update with time step Δt:
   v_casing(t+Δt) = v_casing(t) + (F_casing / M_casing) Δt
   v_cable(t+Δt)  = v_cable(t)  + (F_cable  / M_cable) Δt
3) Constant-power case (instantaneous):
   F = P / v, using v_casing to decelerate casing and v_cable to accelerate cable, capped by F_max if given.
4) Angular-momentum check (end state; neglecting external torques):
   M_tot v_orbit = M_casing v_casing_final + M_cable v_cable_final  =>  v_cable_final = (M_tot v_orbit - M_casing v_casing_final)/M_cable

Usage
-----
$ python orbital_ring_sim.py --help

Examples
--------
Run with defaults, verbose step-by-step prints and a small constant force:
$ python orbital_ring_sim.py --model constant_force --force 1.0e6 --dt 0.5 --verbose --explain

Run constant power with a plot of velocities:
$ python orbital_ring_sim.py --model constant_power --power 5.0e9 --dt 0.25 --plot velocities

Load parameters from YAML:
$ python orbital_ring_sim.py --config example_config.yaml --plot velocities

Outputs
-------
- A summary JSON (printed) of start/end states
- Optional plots saved as PNG files if --save-plots is provided, otherwise shown interactively

License
-------
MIT (see README.md).

"""

from __future__ import annotations

import argparse
import math
import sys
import json
#import time

from dataclasses import dataclass, asdict
from typing import List, Optional, Dict, Any #, Tuple

# Optional YAML support (only if PyYAML is installed). We fail gracefully if absent.
try:
    import yaml  # type: ignore
    _HAVE_YAML = True
except Exception:
    _HAVE_YAML = False

# Optional plotting
try:
    import matplotlib.pyplot as plt  # noqa: F401
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False


# ------------------------------
# Constants and default settings
# ------------------------------

R_EARTH = 6_378_137.0         # m
ALTITUDE = 250_000.0          # m
R_ORBIT = R_EARTH + ALTITUDE  # m
G_LOCAL = 9.073               # m/s^2 (given for 250 km altitude)
L_ORBIT = 2.0 * math.pi * R_ORBIT

V_ORBIT = 7_754.866           # m/s  (initial casing & cable speed)
V_CASING_FINAL = 483.331      # m/s  (ground-synchronous speed at equator frame)

# Defaults for mass distribution (can be overridden via CLI or YAML)
M_CABLE_PER_M = 6_700.0       # kg/m (inner cable)
M_CASING_PER_M = 12_000.0     # kg/m (outer casing)

# Derived total masses (uniform around the ring)
def ring_mass(per_meter: float) -> float:
    return per_meter * L_ORBIT


@dataclass
class Params:
    # Geometry / environment
    r_orbit: float = R_ORBIT
    l_orbit: float = L_ORBIT
    g_local: float = G_LOCAL

    # Kinematic targets (fixed by altitude)
    v_orbit: float = V_ORBIT
    v_casing_final: float = V_CASING_FINAL

    # Mass distribution
    m_cable_per_m: float = M_CABLE_PER_M
    m_casing_per_m: float = M_CASING_PER_M

    # Simulation controls
    model: str = "constant_force"   # or "constant_power"
    force: float = 1.0e6            # N  (used if model == "constant_force")
    power: float = 5.0e9            # W  (used if model == "constant_power")
    fmax: Optional[float] = None    # N, optional cap in constant_power model
    dt: float = 0.5                 # s
    max_time: float = 1.0e8         # s (failsafe)
    tol_v: float = 1e-3             # m/s tolerance for stopping condition
    verbose: bool = False           # print per-step logs
    explain: bool = False           # print equations/explanations at start

    # NEW: logging controls
    log: str = "coarse"             # "coarse", "interval", or "none"
    log_interval: float = 60.0      # seconds, used only if log == "interval"

    # Plotting options
    plot: List[str] = None          # e.g., ["velocities", "tension", "weight"]
    save_plots: bool = False        # if True, saves PNGs instead of interactive show
    outdir: str = "."               # where to save plots

    # Tension/materials (future expansion; partial stubs provided)
    cable_area_m2: float = 15.0     # m^2 CNT effective area (aggregate)
    cable_sigma_break: float = 25e9 # Pa, breaking strength


@dataclass
class State:
    t: float
    v_casing: float
    v_cable: float
    # Placeholders for future: tension and net radial loads
    tension_total: Optional[float] = None
    net_weight_sign: Optional[int] = None


@dataclass
class LimParams:
    tau_p: float = 50.0          # pole pitch [m]
    w_coil: float = 1.0          # coil width [m]
    gap: float = 0.2             # coil-to-plate gap [m]
    n_turns: int = 8
    lim_spacing: float = 500.0   # spacing between LIMs [m]
    max_site_power: float = 8.0e6  # per-site cap [W]
    inv_eff: float = 0.90
    lim_eff: float = 0.95

# Add to Params:
lim: LimParams = LimParams()
i_target: float = 650.0
i_min: float = 10.0
v_slip_max: float = 460.0
v_slip_min: float = 10.0
volts_max: float = 2000.0

# LIM helpers (trimmed, mirroring your formulas)
def _lim_sites(params: Params) -> int:
    return int(round(params.l_orbit / params.lim.lim_spacing))

MU0 = 4 * math.pi * 1e-7

def _b_plate_peak(i_peak: float, p: Params) -> float:
    return MU0 * p.lim.n_turns * i_peak / p.lim.gap

def _f_slip_from_v_slip(v_slip: float, p: Params) -> float:
    return v_slip / (2.0 * p.lim.tau_p)

def _f_supply_from_v_wave(v_wave: float, p: Params) -> float:
    return v_wave / (2.0 * p.lim.tau_p)

def _plate_eddy_volume(p: Params) -> float:
    L_active = p.lim.tau_p * 3  # pitch_count=3
    return p.lim.w_coil * L_active * min(0.040,  # T_PLATE
                                         math.sqrt((2*2.86e-8)/(2*math.pi*MU0*1.0)))  # crude eff thickness

def _set_r_r(p: Params, f_slip: float) -> float:
    L_p = math.pi/2 * p.lim.tau_p
    t_eff = min(0.040, math.sqrt((2*2.86e-8)/(2*math.pi*MU0*max(f_slip, 0.01))))
    return (2.86e-8 * L_p) / (t_eff * p.lim.tau_p)

def _f_thrust(v_slip: float, v_rel: float, i_peak: float, p: Params) -> float:
    b_plate = _b_plate_peak(i_peak, p)
    f_slip = _f_slip_from_v_slip(v_slip, p)
    f_supply = _f_supply_from_v_wave(v_slip + v_rel, p)
    slip_ratio = f_slip / max(f_supply, 1e-6)
    rr = _set_r_r(p, f_slip)
    v_plate = _plate_eddy_volume(p)
    return ((b_plate**2) * (slip_ratio/(1+slip_ratio**2)) / rr) * v_plate


def print_banner(params: Params) -> None:
    print("=" * 79)
    print("Orbital Ring Deployment Simulator (tangential momentum exchange via LIMs)")
    print("=" * 79)
    print("Fixed-by-altitude kinematics:")
    print(f"  - v_orbit (initial):           {params.v_orbit:.3f} m/s")
    print(f"  - v_casing_final (target):     {params.v_casing_final:.3f} m/s")
    print(f"  - r_orbit:                     {params.r_orbit:.0f} m")
    print(f"  - l_orbit:                     {params.l_orbit:.0f} m")
    print(f"  - g_local:                     {params.g_local:.3f} m/s^2")
    print()
    print("Mass model (uniform around ring):")
    print(f"  - m_cable_per_m:               {params.m_cable_per_m:.1f} kg/m")
    print(f"  - m_casing_per_m:              {params.m_casing_per_m:.1f} kg/m")
    M_cable = ring_mass(params.m_cable_per_m)
    M_casing = ring_mass(params.m_casing_per_m)
    print(f"  - M_cable (total):             {M_cable:.3e} kg")
    print(f"  - M_casing (total):            {M_casing:.3e} kg")
    print(f"  - M_total:                     {(M_cable+M_casing):.3e} kg")
    print()
    print("Simulation controls:")
    print(f"  - model:                       {params.model}")
    if params.model == "constant_force":
        print(f"  - force:                       {params.force:.3e} N")
    else:
        print(f"  - power:                       {params.power:.3e} W")
        if params.fmax is not None:
            print(f"  - fmax (cap):                  {params.fmax:.3e} N")
    print(f"  - dt:                          {params.dt:.3f} s")
    print(f"  - max_time:                    {params.max_time:.3e} s")
    print()


def print_explanations() -> None:
    lines = [
        "We integrate the coupled first-order ODEs for casing and cable tangential velocities under an",
        "internal equal-and-opposite linear-induction force:",
        "",
        "    F_cable = -F_casing  (Newton's third law)",
        "",
        "    dv_casing/dt = F_casing / M_casing",
        "    dv_cable/dt  = F_cable  / M_cable",
        "",
        "For the constant-force model, F_casing is a negative constant (deceleration), and F_cable is the",
        "positive of that constant (acceleration).",
        "",
        "For the constant-power model, we use instantaneous forces:",
        "",
        "    F_casing = - P / max(v_casing, v_eps)",
        "    F_cable  = + P / max(v_cable,  v_eps)",
        "",
        "with an optional cap |F| <= fmax. P is an input power that the LIMs transfer between the two",
        "subsystems through internal reaction forces. This is a simplified 'work-limited' scheme.",
        "",
        "Stopping condition:",
        "",
        "    v_casing <= v_casing_final + tol_v",
        "",
        "Angular-momentum consistency check (end state; neglecting external torques):",
        "",
        "    M_total * v_orbit = M_casing * v_casing_final + M_cable * v_cable_final",
        "",
        "so that",
        "",
        "    v_cable_final = (M_total * v_orbit - M_casing * v_casing_final) / M_cable",
        "",
        "Future extension (tension & radial balance — placeholders in code):",
        "",
        "- Cable tension ~ ∫ (m_per_m * v_cable^2 / r) dl  (first-order, centripetal requirement)",
        "- Net 'weight' sign: compare outward centrifugal term m_per_m * v^2 / r against inward gravity",
        "  and casing/coupler loads; track sign changes over time to flag dangerous regimes."
    ]
    print("\n".join(lines))
    print()

def _downsample(x, y, max_points=5000):
    """Return roughly-uniformly sampled series with at most max_points."""
    n = len(x)
    if n <= max_points:
        return x, y
    step = max(1, n // max_points)
    return x[::step], y[::step]

def maybe_plot(params: Params, ts: List[float], vc: List[float], vb: List[float],
               tensions: Optional[List[float]] = None,
               weight_sign: Optional[List[int]] = None) -> None:
    if not params.plot:
        return
    if not _HAVE_MPL:
        print("matplotlib is not available; skipping plots.", file=sys.stderr)
        return

    # Convert seconds -> 30-day “months”
    SEC_PER_MONTH = 30.0 * 86400.0
    t_months = [t / SEC_PER_MONTH for t in ts]

    # Downsample to speed up plotting
    t_m_ds, vc_ds = _downsample(t_months, vc, max_points=5000)
    _,       vb_ds = _downsample(t_months, vb, max_points=5000)

    if "velocities" in params.plot:
        import matplotlib.pyplot as plt  # noqa: F811
        from matplotlib.ticker import MaxNLocator

        plt.figure()
        plt.plot(t_m_ds, vc_ds, label="v_casing [m/s]")
        plt.plot(t_m_ds, vb_ds, label="v_cable [m/s]")
        plt.xlabel("time [30-day months]")
        plt.ylabel("speed [m/s]")
        plt.title("Casing and Cable Speeds vs Time")
        ax = plt.gca()
        # Force integer-like monthly ticks (e.g., 0, 1, 2, …)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins='auto', prune=None))
        plt.legend(loc="upper right")  # avoid slow 'best'

        if params.save_plots:
            fname = f"{params.outdir.rstrip('/')}/velocities.png"
            plt.savefig(fname, dpi=150, bbox_inches="tight")
            print(f"Saved plot: {fname}")
            plt.close()
        else:
            plt.show()

    if "tension" in params.plot and tensions is not None:
        import matplotlib.pyplot as plt  # noqa: F811
        t_m_ds, tens_ds = _downsample(t_months, tensions, max_points=5000)
        plt.figure()
        plt.plot(t_m_ds, tens_ds, label="approx. total cable tension [N]")
        plt.xlabel("time [30-day months]")
        plt.ylabel("tension [N]")
        plt.title("Approximate Cable Tension vs Time (first-order)")
        plt.legend(loc="upper right")
        if params.save_plots:
            fname = f"{params.outdir.rstrip('/')}/tension.png"
            plt.savefig(fname, dpi=150, bbox_inches="tight")
            print(f"Saved plot: {fname}")
            plt.close()
        else:
            plt.show()

    if "weight" in params.plot and weight_sign is not None:
        import matplotlib.pyplot as plt  # noqa: F811
        t_m_ds, ws_ds = _downsample(t_months, weight_sign, max_points=5000)
        plt.figure()
        plt.plot(t_m_ds, ws_ds, label="net radial sign (+out, -in)")
        plt.xlabel("time [30-day months]")
        plt.ylabel("sign")
        plt.title("Net Radial Sign vs Time (stub)")
        plt.legend(loc="upper right")
        if params.save_plots:
            fname = f"{params.outdir.rstrip('/')}/weight_sign.png"
            plt.savefig(fname, dpi=150, bbox_inches="tight")
            print(f"Saved plot: {fname}")
            plt.close()
        else:
            plt.show()


def approximate_tension(params: Params, v_cable: float) -> float:
    """
    Placeholder: integrate centripetal load of the cable around the ring.
    For a uniform ring, total inward radial force requirement ~ (M_cable * v_cable^2) / r.
    Tension is the internal force carrying this load; for a closed loop, max section tension is of order:
        T ~ (M_cable * v_cable^2) / (2π)  (very rough, TODO: refine with proper ring mechanics)
    We return a conservative proxy using total inward requirement M_cable*v^2/r as 'effective' tension.
    """
    M_cable = ring_mass(params.m_cable_per_m)
    return M_cable * v_cable * v_cable / params.r_orbit


def approximate_weight_sign(params: Params, v_cable: float, v_casing: float) -> int:
    """
    Stub: sign of net radial load (+1 outward, -1 inward), comparing centrifugal to gravity at altitude.
    This is schematic and should be expanded with proper mass-weighting and coupler mechanics.
    """
    mpm_total = params.m_cable_per_m + params.m_casing_per_m
    centripetal_per_m = mpm_total * (0.5*(v_cable + v_casing))**2 / params.r_orbit
    gravity_per_m = mpm_total * params.g_local
    return +1 if centripetal_per_m > gravity_per_m else -1


def simulate(params: Params) -> Dict:
    if params.verbose:
        print_banner(params)
        if params.explain:
            print_explanations()

    # Initial conditions
    v_casing = params.v_orbit
    v_cable = params.v_orbit

    M_casing = ring_mass(params.m_casing_per_m)
    M_cable  = ring_mass(params.m_cable_per_m)
    M_total  = M_casing + M_cable

    # For logging
    ts: List[float] = [0.0]
    v_casing_hist: List[float] = [v_casing]
    v_cable_hist:  List[float] = [v_cable]
    tensions: List[float] = []
    weight_signs: List[int] = []

    t = 0.0
    steps = 0
    v_eps = 1e-3

    if params.verbose:
        print("Starting time loop... (Ctrl+C to stop)")
    
    # Progress logging scheduler
    def _next_coarse_target(current_t: float) -> float:
        # 1 min, 1 hr, 1 day, 1 week, then every 30 days
        if current_t < 60.0:
            return 60.0
        if current_t < 3600.0:
            return 3600.0
        if current_t < 86400.0:
            return 86400.0
        if current_t < 604800.0:
            return 604800.0
        month = 30.0 * 86400.0
        k = math.floor(current_t / month) + 1
        return k * month

    if params.verbose and params.log != "none":
        if params.log == "interval":
            next_log_t = params.log_interval
        else:  # coarse
            next_log_t = _next_coarse_target(0.0)
    elif params.model == "lim_simple":
        # very simple controller like yours
        sites = _lim_sites(params)
        i_peak = params.i_target * 0.7
        v_slip = params.v_slip_min
        # per-step:
        v_rel = max(v_cable - v_casing, 0.0)
        thrust_one = _f_thrust(v_slip, v_rel, i_peak, params)
        thrust_total = thrust_one * sites * 2

        # crude site power estimate and ramp towards cap
        p_thrust_one = thrust_one * max(v_rel, 1.0)
        p_heat_one = 0.2 * p_thrust_one              # placeholder
        p_site = (2 * (p_thrust_one + p_heat_one)) * (2 - params.lim.lim_eff)
        p_site *= (2 - params.lim.inv_eff)

        if p_site < params.lim.max_site_power * 0.9:
            # ramp up
            i_peak = min(params.i_target, i_peak + 0.5)
            v_slip = min(params.v_slip_max, v_slip + 1.0)
        elif p_site > params.lim.max_site_power:
            # ramp down
            i_peak = max(params.i_min, i_peak * 0.98)
            v_slip = max(params.v_slip_min, v_slip * 0.98)

        # update velocities with total thrust
        v_casing += (-thrust_total / ring_mass(params.m_casing_per_m)) * params.dt
        v_cable  += (+thrust_total / ring_mass(params.m_cable_per_m))  * params.dt

    else:
        next_log_t = float("inf")


    try:
        while t < params.max_time and (v_casing - params.v_casing_final) > params.tol_v:
            # Determine instantaneous forces
            if params.model == "constant_force":
                F_casing = -abs(params.force)
                F_cable  = +abs(params.force)
            elif params.model == "constant_power":
                F_mag_casing = abs(params.power) / max(v_casing, v_eps)
                F_mag_cable  = abs(params.power) / max(v_cable,  v_eps)
                if params.fmax is not None:
                    F_mag_casing = min(F_mag_casing, params.fmax)
                    F_mag_cable  = min(F_mag_cable,  params.fmax)
                F_casing = -F_mag_casing
                F_cable  = +F_mag_cable
            else:
                raise ValueError(f"Unknown model: {params.model}")

            # Update velocities
            v_casing += (F_casing / M_casing) * params.dt
            v_cable  += (F_cable  / M_cable)  * params.dt

            # Advance time
            t += params.dt
            steps += 1

            # Store
            ts.append(t)
            v_casing_hist.append(v_casing)
            v_cable_hist.append(v_cable)

            # Progress log (lightweight)
            if params.verbose and params.log != "none":
                if t >= next_log_t:
                    print(f"t={t:12.2f} s | v_casing={v_casing:9.3f} m/s | v_cable={v_cable:10.3f} m/s")
                    if params.log == "interval":
                        next_log_t += params.log_interval
                    else:
                        next_log_t = _next_coarse_target(next_log_t)

            # First-order "tension" and "weight sign" (placeholders)
            T_eff = approximate_tension(params, v_cable)
            W_sig = approximate_weight_sign(params, v_cable, v_casing)
            tensions.append(T_eff)
            weight_signs.append(W_sig)

        if params.verbose:
            print("Time loop complete.")
    except KeyboardInterrupt:
        print("\nSimulation interrupted by user.", file=sys.stderr)

    # Angular momentum check for v_cable_final (closed system, no external torque)
    v_cable_final_AM = (M_total * params.v_orbit - M_casing * params.v_casing_final) / M_cable

    summary = {
        "params": asdict(params),
        "results": {
            "time_elapsed_s": t,
            "steps": steps,
            "v_casing_end": v_casing,
            "v_cable_end": v_cable,
            "v_cable_end_angular_momentum_check": v_cable_final_AM,
            "conservation_delta": (M_total*params.v_orbit - (M_casing*v_casing + M_cable*v_cable))
        }
    }

    # Plot if requested
    maybe_plot(params, ts, v_casing_hist, v_cable_hist,
               tensions=tensions if tensions else None,
               weight_sign=weight_signs if weight_signs else None)

    # Print JSON summary
    print(json.dumps(summary, indent=2))

    return summary


def _coerce_like(current: Any, value: Any) -> Any:
    """
    Coerce YAML-loaded 'value' to the same kind of type as 'current'.
    Handles float/int/bool, Optional[float]/Optional[int], lists (including comma strings),
    and leaves other types as-is.
    """
    if value is None:
        return None

    # If current is None, try to infer from the YAML value
    if current is None:
        # try numeric strings like "5e9"
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("none", "null"):
                return None
            # bools
            if low in ("true", "false", "yes", "no", "on", "off", "1", "0"):
                return low in ("true", "yes", "on", "1")
            # numbers
            try:
                if any(ch in low for ch in ".e"):
                    return float(low)
                return int(low)
            except Exception:
                # maybe a comma list
                if "," in value:
                    return [s.strip() for s in value.split(",") if s.strip()]
                return value
        # lists pass through
        if isinstance(value, list):
            return value
        return value

    # If we know the current attribute type, mimic it
    t = type(current)

    # Lists: allow list OR comma-separated string
    if isinstance(current, list):
        if isinstance(value, list):
            return list(value)
        if isinstance(value, str):
            return [s.strip() for s in value.split(",") if s.strip()]
        # single item -> wrap
        return [value]

    # Bools
    if isinstance(current, bool):
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        if isinstance(value, str):
            return value.strip().lower() in ("true", "yes", "on", "1")
        return bool(value)

    # Floats (and allow string scientific notation)
    if isinstance(current, float):
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("none", "null"):
                return None
            try:
                return float(low)
            except Exception:
                raise ValueError(f"Expected float-like value, got: {value!r}")
        raise ValueError(f"Expected float-like value, got: {type(value).__name__}")

    # Ints
    if isinstance(current, int):
        if isinstance(value, (int, float)):
            return int(value)
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("none", "null"):
                return None
            # allow scientific notation by casting through float
            try:
                return int(float(low))
            except Exception:
                raise ValueError(f"Expected int-like value, got: {value!r}")
        raise ValueError(f"Expected int-like value, got: {type(value).__name__}")

    # Fallback: just return it
    return value

def load_params_from_yaml(path: str, base: Params) -> Params:
    if not _HAVE_YAML:
        raise RuntimeError("PyYAML is not installed; cannot load YAML configs.")
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    for k, v in data.items():
        if hasattr(base, k):
            current = getattr(base, k)
            try:
                coerced = _coerce_like(current, v)
            except Exception as e:
                print(f"Warning: could not coerce '{k}' value {v!r}: {e}", file=sys.stderr)
                coerced = v  # last resort
            setattr(base, k, coerced)
        else:
            print(f"Warning: unknown config key '{k}' ignored.", file=sys.stderr)
    return base

def parse_args(argv: Optional[List[str]] = None) -> Params:
    p = argparse.ArgumentParser(
        description="Orbital ring deployment simulator (casing-cable momentum exchange via LIMs)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument("--config", type=str, default=None,
                   help="YAML config file with parameters")

    # Mass model
    p.add_argument("--m-cable-per-m", type=float, default=M_CABLE_PER_M,
                   help="Cable mass per meter [kg/m]")
    p.add_argument("--m-casing-per-m", type=float, default=M_CASING_PER_M,
                   help="Casing mass per meter [kg/m]")

    # Simulation controls
    p.add_argument("--model", type=str, choices=["constant_force", "constant_power"],
                   default="constant_force", help="Thrust/power scheduling model")
    p.add_argument("--force", type=float, default=1.0e6,
                   help="Force magnitude [N] if model=constant_force")
    p.add_argument("--power", type=float, default=5.0e9,
                   help="Power [W] if model=constant_power")
    p.add_argument("--fmax", type=float, default=None,
                   help="Optional cap on instantaneous |F| in constant_power model [N]")
    p.add_argument("--dt", type=float, default=0.5, help="Time step [s]")
    p.add_argument("--max-time", type=float, default=1.0e7, help="Max simulated time [s]")
    p.add_argument("--tol-v", type=float, default=1e-3, help="Stop when v_casing <= v_final + tol-v")

    # Output / verbosity
    p.add_argument("--verbose", action="store_true", help="Print per-second status lines")
    p.add_argument("--explain", action="store_true", help="Print governing equations and notes")
    p.add_argument("--log", type=str, choices=["coarse", "interval", "none"],
               default="coarse", help="Progress logging mode")
    p.add_argument("--log-interval", type=float, default=60.0,
               help="If --log=interval, print every N seconds")

    # Plotting
    p.add_argument("--plot", type=str, default=None,
                   help="Comma-separated list of plots to show: velocities,tension,weight")
    p.add_argument("--save-plots", action="store_true", help="Save plots to PNG files instead of showing")
    p.add_argument("--outdir", type=str, default=".", help="Directory for saved plots")

    # Materials (for future tension checks)
    p.add_argument("--cable-area-m2", type=float, default=15.0, help="Effective CNT area [m^2]")
    p.add_argument("--cable-sigma-break", type=float, default=25e9, help="Breaking strength [Pa]")

    args = p.parse_args(argv)

    params = Params()
    params.m_cable_per_m = args.m_cable_per_m
    params.m_casing_per_m = args.m_casing_per_m
    params.model = args.model
    params.force = args.force
    params.power = args.power
    params.fmax = args.fmax
    params.dt = args.dt
    params.max_time = args.max_time
    params.tol_v = args.tol_v
    params.verbose = args.verbose
    params.explain = args.explain
    params.plot = [s.strip() for s in args.plot.split(",")] if args.plot else []
    params.save_plots = args.save_plots
    params.outdir = args.outdir
    params.cable_area_m2 = args.cable_area_m2
    params.cable_sigma_break = args.cable_sigma_break
    params.log = args.log
    params.log_interval = args.log_interval

    # Load YAML overrides last
    if args.config:
        params = load_params_from_yaml(args.config, params)

    return params


def main(argv: Optional[List[str]] = None) -> int:
    params = parse_args(argv)
    simulate(params)
    return 0


if __name__ == "__main__":
    sys.exit(main())
