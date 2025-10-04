#!/usr/bin/env python3
from __future__ import annotations
import argparse, sys, math, json
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any
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

MU0 = 4*math.pi*1e-7
SEC_PER_MONTH = 30.0 * 86400.0
R_EARTH = 6_378_137.0
ALT = 250_000.0
R_ORBIT = R_EARTH + ALT
G_LOCAL = 9.073
L_ORBIT = 2*math.pi*R_ORBIT
V_ORBIT = 7_754.866
V_CASING_FINAL = 483.331

@dataclass
class LIM:
    tau_p: float = 50.0
    w_coil: float = 1.0
    gap: float = 0.20
    pitch_count: int = 3
    spacing: float = 500.0
    casing_outer_w: float = 10.0
    n_turns: int = 8
    d_kapton_mm: float = 0.01
    k_fill: float = 0.1667
    ic: float = 800.0
    w_tape: float = 0.012
    alpha_angle_deg: float = 20.0
    volts_max_user: float = 2000.0
    v_slip_min: float = 10.0
    v_slip_max: float = 460.0
    i_peak_target: float = 650.0
    i_peak_min: float = 10.0
    slip_min: float = 0.01
    v_rel_min: float = 10.0
    inv_eff: float = 0.90
    lim_eff: float = 0.95
    max_site_power: float = 8.0e6
    t_plate: float = 0.040
    rho_alu_e: float = 2.86e-8
    em_alu: float = 0.85
    em_heat_sink: float = 0.90
    cryo_eff: float = 0.05
    heat_sink_l: float = 30.0
    q_sun_1au: float = 1361.0
    q_earth_day: float = 650.0
    q_shielding: float = 0.005
    di_max_per_s: float = 1.0
    dvslip_max_per_s: float = 2.0
    thrust_coeff: float = 1.0

@dataclass
class Params:
    m_cable_per_m: float = 26_300.0
    m_casing_per_m: float = 5_000.0
    r_orbit: float = R_ORBIT
    l_orbit: float = L_ORBIT
    g_local: float = G_LOCAL
    v_orbit: float = V_ORBIT
    v_casing_final: float = V_CASING_FINAL
    cable_area_m2: float = 15.0
    sigma_break: float = 25e9
    lim: LIM = field(default_factory=LIM)
    dt: float = 1.0
    max_time: float = 2.0e8
    tol_v: float = 1e-3
    verbose: bool = False
    explain: bool = False
    log: str = "coarse"
    log_interval: float = 600.0
    plot_skip_s: float = 600.0
    plot: List[str] = field(default_factory=list)
    save_plots: bool = False
    outdir: str = "."

def ring_mass(per_m: float, p: Params) -> float:
    return per_m * p.l_orbit

def _coerce_like(cur: Any, val: Any) -> Any:
    if val is None:
        return cur
    if isinstance(cur, list):
        if isinstance(val, list): return list(val)
        if isinstance(val, str): return [s.strip() for s in val.split(",") if s.strip()]
        return [val]
    if isinstance(cur, bool):
        if isinstance(val, bool): return val
        if isinstance(val, (int, float)): return bool(val)
        if isinstance(val, str): return val.strip().lower() in ("true","yes","on","1")
        return bool(val)
    if isinstance(cur, float):
        if isinstance(val, (int,float)): return float(val)
        if isinstance(val, str): return float(val.strip())
    if isinstance(cur, int):
        if isinstance(val, (int,float)): return int(val)
        if isinstance(val, str): return int(float(val.strip()))
    return val

def load_yaml(path: str, base: Params) -> Params:
    if not HAVE_YAML: raise RuntimeError("PyYAML not installed.")
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
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

def banner(p: Params):
    print("="*79)
    print("Orbital Ring Deployment — per-site LIM physics + coil limits + power caps")
    print("="*79)
    print(f"Altitude-fixed speeds: v_orbit={p.v_orbit:.3f} m/s  v_casing_final={p.v_casing_final:.3f} m/s")
    print(f"r_orbit={p.r_orbit:.0f} m  L={p.l_orbit:.0f} m  g={p.g_local:.3f} m/s^2")
    print("Mass model (kg/m): cable={:.1f}  casing={:.1f}".format(p.m_cable_per_m, p.m_casing_per_m))
    Mc, Ms = ring_mass(p.m_cable_per_m, p), ring_mass(p.m_casing_per_m, p)
    print("Totals (kg): cable={:.3e} casing={:.3e} total={:.3e}".format(Mc, Ms, Mc+Ms))
    print("LIM: tau_p={:.2f} m  w={:.2f} m  gap={:.3f} m  turns={}  spacing={} m  pitches={}".format(
        p.lim.tau_p, p.lim.w_coil, p.lim.gap, p.lim.n_turns, int(p.lim.spacing), p.lim.pitch_count))
    print("LIM: volts_max_user={:.0f} V  i_target={:.1f} A  v_slip=[{:.1f},{:.1f}] m/s  P_site_max={:.2e} W".format(
        p.lim.volts_max_user, p.lim.i_peak_target, p.lim.v_slip_min, p.lim.v_slip_max, p.lim.max_site_power))
    print()

def explain_text():
    print("Force & power skeleton (your structure):")
    print("  B_peak = μ0 * N_turns * I_peak / gap")
    print("  f_slip = v_slip/(2τp), f_sup = (v_rel+v_slip)/(2τp), s = f_slip/f_sup")
    print("  Eddy skin depth δ = sqrt(2 ρ_Al / (2π μ0 f_slip)) ; t_eff = min(t_plate, δ)")
    print("  V_eddy = w_coil * (pitch_count τp) * t_eff")
    print("  r_r ~ (ρ_Al * L_p)/(t_eff τp), L_p = (π/2) τp")
    print("  F_site ≈ thrust_coeff * (B^2) * (s/(1+s^2)) * V_eddy / r_r")
    print("  P_site includes thrust + eddy + hysteresis + cryo, / (η_inv η_lim)")
    print("  V_site (per phase) ≈ (2π f_sup L_coil I_peak)/√3, L_coil ≈ N^2 μ0 (A_coil/gap) k_fill")
    print()
    print("Net load and tension (outward/up is positive):")
    print("  a_out,cable  = v_cable^2 / r ; a_out,casing = v_casing^2 / r")
    print("  w' [N/m] = m'_cable (a_out,cable - g) + m'_casing (a_out,casing - g)")
    print("  Hoop tension: T = max(0, w' * r), stress σ = T/A_cable")
    print()

# Helpers
def lim_sites(p: Params) -> int: return int(round(p.l_orbit / p.lim.spacing))
def l_active(p: Params) -> float: return p.lim.tau_p * p.lim.pitch_count
def a_coil(p: Params) -> float: return p.lim.tau_p * p.lim.w_coil
def l_p(p: Params) -> float: return math.pi/2.0 * p.lim.tau_p
def kapton_safe_v(p: Params) -> float: return 1e5 * p.lim.d_kapton_mm
def volts_max_allowed(p: Params) -> float:
    kapton_v = max(0, p.lim.n_turns-1) * kapton_safe_v(p)
    return min(p.lim.volts_max_user, kapton_v if kapton_v > 0 else p.lim.volts_max_user)
def b_plate_peak(i_peak: float, p: Params) -> float: return MU0 * p.lim.n_turns * i_peak / p.lim.gap
def b_coil_peak(i_peak: float, p: Params) -> float:
    return MU0 * p.lim.n_turns * i_peak / (math.pi * (p.lim.w_coil + p.lim.tau_p))
def f_from_v(v: float, p: Params) -> float: return max(v, 1e-9) / (2.0 * p.lim.tau_p)
def skin_depth(f_slip: float, p: Params) -> float:
    return math.sqrt((2.0 * p.lim.rho_alu_e) / (2.0 * math.pi * MU0 * max(f_slip, 1e-9)))
def eff_plate_depth(f_slip: float, p: Params) -> float: return min(p.lim.t_plate, skin_depth(f_slip, p))
def plate_eddy_volume(f_slip: float, p: Params) -> float: return p.lim.w_coil * l_active(p) * eff_plate_depth(f_slip, p)
def set_r_r(f_slip: float, p: Params) -> float: return (p.lim.rho_alu_e * l_p(p)) / (eff_plate_depth(f_slip, p) * p.lim.tau_p)
def slip_ratio(v_slip: float, v_rel: float, p: Params) -> float:
    return f_from_v(v_slip, p) / max(f_from_v(v_rel + v_slip, p), 1e-9)

def site_thrust(i_peak: float, v_rel: float, v_slip: float, p: Params) -> float:
    B = b_plate_peak(i_peak, p); s = slip_ratio(v_slip, v_rel, p)
    rr = set_r_r(f_from_v(v_slip, p), p); V_eddy = plate_eddy_volume(f_from_v(v_slip, p), p)
    return p.lim.thrust_coeff * (B*B) * (s/(1.0 + s*s)) * (V_eddy / max(rr, 1e-12))

def site_voltage_estimate(i_peak: float, v_rel: float, v_slip: float, p: Params) -> float:
    f_sup = f_from_v(v_rel + v_slip, p)
    L_coil = (p.lim.n_turns**2) * MU0 * (a_coil(p) / p.lim.gap) * p.lim.k_fill
    return (2.0 * math.pi * f_sup * L_coil * i_peak) / math.sqrt(3.0)

def site_power_estimate(i_peak: float, v_rel: float, v_slip: float, p: Params) -> float:
    F = site_thrust(i_peak, v_rel, v_slip, p)
    v_rel_eff = max(v_rel, 0.1)
    P_thrust = F * v_rel_eff
    B = b_plate_peak(i_peak, p)
    f_sl = f_from_v(v_slip, p)
    t_eff = eff_plate_depth(f_sl, p)
    V_eddy = plate_eddy_volume(f_sl, p)
    P_eddy = (math.pi**2 / (6.0 * p.lim.rho_alu_e)) * (B*B) * (t_eff**2) * (f_sl**2) * V_eddy
    Bc = b_coil_peak(i_peak, p)
    q = Bc * p.lim.ic * p.lim.w_tape * math.sin(math.radians(p.lim.alpha_angle_deg))
    L_HTS_COIL = 2.0 * (p.lim.w_coil + p.lim.tau_p) * p.lim.n_turns
    P_hyst = q * L_HTS_COIL * f_from_v(v_rel + v_slip, p) * 3 * p.lim.pitch_count
    Q_abs_m = (p.lim.q_sun_1au + p.lim.q_earth_day) * p.lim.casing_outer_w * p.lim.q_shielding
    Q_abs_site = Q_abs_m * p.lim.spacing
    P_cryo = p.lim.cryo_eff * (P_eddy + P_hyst + Q_abs_site)
    return (P_thrust + P_eddy + P_hyst + P_cryo) / (p.lim.inv_eff * p.lim.lim_eff)

def downsample(x: List[float], y: List[float], max_points: int = 5000):
    n = len(x)
    if n <= max_points: return x, y
    step = max(1, n // max_points)
    return x[::step], y[::step]

def _align_tail(x, y):
    n = min(len(x), len(y))
    return x[-n:], y[-n:]

def _prep_xy_for_plot(x_months, y_series, max_points=5000):
    xx, yy = _align_tail(x_months, y_series)
    return downsample(xx, yy, max_points=max_points)

def annotate_last(ax, x, y):
    if not x or not y: return
    ax.annotate(f"{y[-1]:.3g}", xy=(x[-1], y[-1]), xytext=(5,5), textcoords="offset points", fontsize=9)

def maybe_plot(p: Params, ts: List[float], vc: List[float], vb: List[float],
               tension: List[float], net_per_m_hist: List[float],
               extras: Dict[str, List[float]] = None) -> None:
    if not p.plot or not HAVE_MPL: return
    extras = extras or {}
    i0 = 0
    while i0 < len(ts) and ts[i0] < p.plot_skip_s:
        i0 += 1
    ts = ts[i0:]; vc = vc[i0:]; vb = vb[i0:]
    tension = tension[i0:] if tension else tension
    net_per_m_hist = net_per_m_hist[i0:] if net_per_m_hist else net_per_m_hist
    t_months = [t/SEC_PER_MONTH for t in ts]

    if "velocities" in p.plot and t_months:
        tm, vc_ds = downsample(t_months, vc)
        _,  vb_ds = downsample(t_months, vb)
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, vc_ds, label="v_casing [m/s]")
        ax.plot(tm, vb_ds, label="v_cable [m/s]")
        annotate_last(ax, tm, vc_ds); annotate_last(ax, tm, vb_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("speed [m/s]")
        ax.set_title("Casing and Cable Speeds vs Time")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/velocities.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "stress" in p.plot and t_months and tension:
        stress = [T/p.cable_area_m2/1e9 for T in tension]
        tm, st_ds = downsample(t_months, stress)
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, st_ds, label="cable stress [GPa]")
        annotate_last(ax, tm, st_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("GPa")
        ax.set_title("Cable Stress vs Time")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/stress.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "net_weight_sign" in p.plot and t_months and net_per_m_hist:
        sign_series = [ -1 if v < 0 else (1 if v > 0 else 0) for v in net_per_m_hist ]
        tm, s_ds = downsample(t_months, sign_series)
        any_neg = any(v < 0 for v in s_ds)
        color = "red" if any_neg else "green"
        fig = plt.figure(); ax = fig.gca()
        ax.step(tm, s_ds, where="post", label="sign (−1,0,+1)", color=color)
        annotate_last(ax, tm, s_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("sign")
        ax.set_yticks([-1, 0, 1]); ax.set_ylim(-1.2, 1.2)
        ax.set_title("Net Weight Sign (out + / in −)")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/net_weight_sign.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "max_load_per_m" in p.plot and t_months and net_per_m_hist:
        load_t_per_m = [max(0.0, v)/p.g_local/1000.0 for v in net_per_m_hist]
        tm, y_ds = downsample(t_months, load_t_per_m)
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, y_ds, label="Maximum Load per Meter [tonne/m]")
        annotate_last(ax, tm, y_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("tonne/m")
        ax.set_title("Maximum Load Capacity per Meter")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/max_load_per_m.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "accel_excess" in p.plot and t_months and extras.get("a_out"):
        a_excess = [v - p.g_local for v in extras["a_out"]]
        tm, y_ds = _prep_xy_for_plot(t_months, a_excess)
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, y_ds, label="a_out - g [m/s²]")
        annotate_last(ax, tm, y_ds)
        ax.axhline(0.0, linestyle="--", linewidth=1, label="zero net")
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("m/s²")
        ax.set_title("Net Radial Acceleration (Outward Positive)")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/accel_excess.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "power" in p.plot and extras.get("power"):
        tm, y_ds = _prep_xy_for_plot(t_months, extras["power"])
        fig = plt.figure(); ax = fig.gca()
        y_mw = [v/1e6 for v in y_ds]
        ax.plot(tm, y_mw, label="P_site per LIM [MW]"); annotate_last(ax, tm, y_mw)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("MW")
        ax.set_title("Per-site Input Power"); ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/power.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "voltage" in p.plot and extras.get("voltage"):
        tm, y_ds = _prep_xy_for_plot(t_months, extras["voltage"])
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, y_ds, label="V_site per phase [V]"); annotate_last(ax, tm, y_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("V")
        ax.set_title("Coil Inductive Voltage"); ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/voltage.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "current_slip" in p.plot and (extras.get("current") or extras.get("v_slip")):
        fig = plt.figure(); ax = fig.gca()
        if extras.get("current"):
            tm_i, i_ds = _prep_xy_for_plot(t_months, extras["current"])
            ax.plot(tm_i, i_ds, label="I_peak [A]"); annotate_last(ax, tm_i, i_ds)
        if extras.get("v_slip"):
            tm_vs, vs_ds = _prep_xy_for_plot(t_months, extras["v_slip"])
            ax.plot(tm_vs, vs_ds, label="v_slip [m/s]"); annotate_last(ax, tm_vs, vs_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("A or m/s")
        ax.set_title("Controller: Current & Slip"); ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/current_slip.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

    if "relative_speed" in p.plot and extras.get("v_rel"):
        tm, y_ds = _prep_xy_for_plot(t_months, extras["v_rel"])
        fig = plt.figure(); ax = fig.gca()
        ax.plot(tm, y_ds, label="v_rel = v_cable - v_casing [m/s]"); annotate_last(ax, tm, y_ds)
        ax.set_xlabel("time [30-day months]"); ax.set_ylabel("m/s")
        ax.set_title("Relative Speed"); ax.legend(loc="upper left")
        if p.save_plots: plt.savefig(f"{p.outdir.rstrip('/')}/relative_speed.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        else: plt.show()

def next_log_target(t: float) -> float:
    if t < 60.0: return 1.0
    if t < 3600.0: return 3600.0
    if t < 86400.0: return 86400.0
    if t < 604800.0: return 604800.0
    month = 30.0 * 86400.0
    k = math.floor(t / month) + 1
    return k * month

def parse_args(argv=None) -> Params:
    ap = argparse.ArgumentParser(description="Orbital ring deployment with per-site LIM physics (v2 plots)")
    ap.add_argument("--config", type=str, default=None, help="YAML config file")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--explain", action="store_true")
    ap.add_argument("--plot", type=str, default=None, help="Comma list e.g. velocities,stress,max_load_per_m,accel_excess,net_weight_sign")
    ap.add_argument("--save-plots", action="store_true")
    ap.add_argument("--outdir", type=str, default=".")
    ap.add_argument("--plot-skip-s", type=float, default=None, help="Seconds to skip at start of plots")
    args = ap.parse_args(argv)
    p = Params()
    p.verbose = args.verbose; p.explain = args.explain
    p.plot = [s.strip() for s in args.plot.split(",")] if args.plot else []
    p.save_plots = args.save_plots; p.outdir = args.outdir
    if args.plot_skip_s is not None: p.plot_skip_s = args.plot_skip_s
    if args.config: p = load_yaml(args.config, p)
    return p

def simulate(p: Params) -> Dict:
    if p.verbose:
        banner(p)
        if p.explain: explain_text()

    v_casing = p.v_orbit
    v_cable  = p.v_orbit

    Mcasing = ring_mass(p.m_casing_per_m, p)
    Mcable  = ring_mass(p.m_cable_per_m, p)
    Mtotal  = Mcasing + Mcable
    sites   = lim_sites(p)

    i_peak = p.lim.i_peak_min
    v_slip = p.lim.v_slip_min

    next_log = (p.log_interval if p.verbose and p.log=="interval"
                else (next_log_target(0.0) if p.verbose and p.log=="coarse" else float("inf")))

    ts=[0.0]; vch=[v_casing]; vbh=[v_cable]
    Thist=[0.0]; NetPerMhist=[- (p.m_cable_per_m + p.m_casing_per_m)*p.g_local]
    Phist=[]; Vhist=[]; Ihist=[]; Vsliphist=[]; Vrelhist=[]; Aouthist=[]

    min_T = (float("inf"), None); max_T = (-float("inf"), None)
    min_net = (float("inf"), None); max_net = (-float("inf"), None)

    t=0.0; steps=0; broke=False; reason=""
    while t < p.max_time and (v_casing - p.v_casing_final) > p.tol_v:
        v_rel = max(v_cable - v_casing, 0.0)

        P_site = site_power_estimate(i_peak, v_rel, v_slip, p)
        V_site = site_voltage_estimate(i_peak, v_rel, v_slip, p)
        V_cap  = volts_max_allowed(p)

        want_i = p.lim.i_peak_target
        want_vslip = p.lim.v_slip_max
        if (P_site > p.lim.max_site_power) or (V_site > V_cap):
            want_i = max(p.lim.i_peak_min, i_peak * 0.98)
            want_vslip = max(p.lim.v_slip_min, v_slip * 0.98)
        di = want_i - i_peak; dv = want_vslip - v_slip
        i_peak += math.copysign(min(abs(di), p.lim.di_max_per_s * p.dt), di)
        v_slip += math.copysign(min(abs(dv), p.lim.dvslip_max_per_s * p.dt), dv)
        i_peak = min(max(i_peak, p.lim.i_peak_min), p.lim.i_peak_target)
        v_slip = min(max(v_slip, p.lim.v_slip_min), p.lim.v_slip_max)

        F_tot = site_thrust(i_peak, v_rel, v_slip, p) * sites * 2.0

        v_casing += (-F_tot / Mcasing) * p.dt
        v_cable  += (+F_tot / Mcable)  * p.dt

        a_out_cable  = (v_cable  * v_cable)  / p.r_orbit
        a_out_casing = (v_casing * v_casing) / p.r_orbit

        net_per_m = (
            p.m_cable_per_m  * (a_out_cable  - p.g_local) +
            p.m_casing_per_m * (a_out_casing - p.g_local)
        )

        T = max(0.0, net_per_m * p.r_orbit)
        sigma = T / p.cable_area_m2
        if sigma > p.sigma_break:
            broke = True; reason = f"Cable failure: stress {sigma:.3e} > {p.sigma_break:.3e} Pa"
            if p.verbose: print(reason); break

        mpm_total = p.m_cable_per_m + p.m_casing_per_m
        a_out_mass_weighted = (p.m_cable_per_m*a_out_cable + p.m_casing_per_m*a_out_casing)/mpm_total

        if T < min_T[0]: min_T = (T, t)
        if T > max_T[0]: max_T = (T, t)
        if net_per_m < min_net[0]: min_net = (net_per_m, t)
        if net_per_m > max_net[0]: max_net = (net_per_m, t)

        t += p.dt; steps += 1

        if steps % 10 == 0:
            ts.append(t); vch.append(v_casing); vbh.append(v_cable)
            Thist.append(T); NetPerMhist.append(net_per_m)
            Phist.append(P_site); Vhist.append(V_site); Ihist.append(i_peak)
            Vsliphist.append(v_slip); Vrelhist.append(v_rel); Aouthist.append(a_out_mass_weighted)

        if p.verbose and t >= next_log:
            print(f"t={t:12.2f} s | months={t/SEC_PER_MONTH:7.2f} | "
                  f"v_casing={v_casing:9.3f} | v_cable={v_cable:10.3f} | "
                  f"Ip={i_peak:6.1f} A | v_slip={v_slip:7.1f} m/s | "
                  f"P_site~{P_site/1e6:6.1f} MW | V_site~{V_site:7.1f} V | V_cap~{V_cap:7.1f} V")
            next_log = (next_log + p.log_interval) if p.log=="interval" else next_log_target(next_log)

    if p.verbose:
        print("Time loop complete." + (" (BROKE)" if broke else ""))

    v_cable_AM = ((Mtotal * p.v_orbit) - (Mcasing * max(v_casing, p.v_casing_final))) / Mcable

    out = {
        "params": asdict(p),
        "results": {
            "time_elapsed_s": t,
            "months_elapsed": t/SEC_PER_MONTH,
            "steps": steps,
            "v_casing_end": v_casing,
            "v_cable_end": v_cable,
            "v_cable_end_angular_momentum_check": v_cable_AM,
            "failed": broke,
            "reason": reason,
            "tension_min_N": min_T[0], "tension_min_at_s": min_T[1],
            "tension_max_N": max_T[0], "tension_max_at_s": max_T[1],
            "net_weight_per_m_min_Npm": min_net[0], "net_weight_min_at_s": min_net[1],
            "net_weight_per_m_max_Npm": max_net[0], "net_weight_max_at_s": max_net[1]
        }
    }
    print(json.dumps(out, indent=2))

    maybe_plot(
        p, ts, vch, vbh, Thist, NetPerMhist,
        extras={
            "power": Phist, "voltage": Vhist, "current": Ihist, "v_slip": Vsliphist,
            "v_rel": Vrelhist, "a_out": Aouthist
        }
    )
    return out

def main(argv=None) -> int:
    p = parse_args(argv); simulate(p); return 0

if __name__ == "__main__":
    sys.exit(main())
