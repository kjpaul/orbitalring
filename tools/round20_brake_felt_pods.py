#!/usr/bin/env python3
"""
Round 20: the braking run of lsm_simulation.py, the felt load taken at the speed in space, and the return launch of empty pods.
Brief: Vol_III/Claude_Code_instructions_round20_2026-11-18.md.

1. lsm_simulation.run_lsm_simulation(brake=True) is run for Chapter 10's design launch.  The braking of the empty sled is now a run
   of the simulation (run_lsm_braking): the LSM as a generator on the launch's limits, from the release with no pause.  It is checked
   here against the braking that round 19 built by hand from Chapter 10's words.
2. The simulation now takes the centrifugal term of the felt load at the speed in space.  The crewed launch at the 3 g limit is run
   again.
3. The empty mass of the 100-seat pod is taken from elevator_system2.size_pod2 (the 78 t pod of Chapter 13) less the script's 130 kg
   for each seat.  One pod and ten pods are launched on the 5,250 t sled to the Hohmann transfer of Table 13.4.
4. The casing's swing (round 18's hold, round 19's integrator, copied here unchanged) on the simulated launch and the simulated
   braking, and one table against round 19.

Run from Python/orbitalring:  python tools/round20_brake_felt_pods.py
Log: runs/round20_logs/round20_brake_felt_pods.txt and .json ; the simulation's own printing for each run: runs/round20_logs/lsm_<case>.txt
"""
import contextlib, io, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.environ.pop("ORBITAL_RING_LSM_STATOR", None)                      # the ruled stator, the default
import numpy as np
_argv = sys.argv; sys.argv = [sys.argv[0]]
import lsm_config as cfg
import lsm_simulation as sim
import ring_altitude as ring
import elevator_system as e1
import elevator_system2 as e2
sys.argv = _argv
assert cfg.STATOR == "ruled"

LOG_DIR = os.path.join(ROOT, "runs", "round20_logs"); os.makedirs(LOG_DIR, exist_ok=True)
L18 = json.load(open(os.path.join(ROOT, "runs", "round18_logs", "round18_swing_hold.json"), encoding="utf-8"))
L19 = json.load(open(os.path.join(ROOT, "runs", "round19_logs", "round19_swing_simulated.json"), encoding="utf-8"))
buf = io.StringIO()
def P(s=""):
    print(s); buf.write(s + "\n"); sys.stdout.flush()
OUT = {}; T0 = time.time()
DEFAULTS = dict(M_SPACECRAFT=cfg.M_SPACECRAFT, V_LAUNCH=cfg.V_LAUNCH, G_LOAD_MAX=cfg.G_LOAD_MAX, A_MAX_G=cfg.A_MAX_G, L_SLED=cfg.L_SLED, B_SLED_NOMINAL=cfg.B_SLED_NOMINAL, P_HVDC_MAX=cfg.P_HVDC_MAX,
                A_BRAKE_MAX_G=cfg.A_BRAKE_MAX_G)

def lsm_run(tag, **over):
    """One run of lsm_simulation.run_lsm_simulation(brake=True) at full time step.  'over' changes only what the command line of lsm_simulation.py can change."""
    for k, v in DEFAULTS.items(): setattr(cfg, k, v)
    for k, v in over.items():
        assert k in DEFAULTS, k; setattr(cfg, k, v)
    cfg.calc_derived()
    cap = io.StringIO()
    with contextlib.redirect_stdout(cap):
        ok = sim.run_lsm_simulation(quick_mode=False, brake=True)
    with open(os.path.join(LOG_DIR, "lsm_%s.txt" % tag), "w", encoding="utf-8") as fh: fh.write(cap.getvalue())
    t = np.array(sim.data_t); F = np.array(sim.data_F); v = np.array(sim.data_v); Pm = np.array(sim.data_P_mech); g = np.array(sim.data_a_occupant_g); ar = np.array(sim.data_a_radial_net)
    lr = dict(sim.launch_result); br = dict(sim.brake_result); dt = cfg.DT
    return dict(tag=tag, ok=bool(ok), t=t, F=F, v=v, g=g, ar=ar, dt=dt, dur=lr["t"], v_end=lr["v"], x_end=lr["x"], v_space=lr["v_space"], a_out_g=lr["a_out_g"], g_peak=lr["g_peak"],
                g_limited=lr["g_load_limited"], g_engage_v=lr["g_load_engage_v"], E_heat=lr["E_heat"], KE=lr["KE"],
                m_total=cfg.M_TOTAL, m_sled=cfg.M_SLED_HARDWARE, F_peak=float(F.max()), F_end=float(F[-1]), P_peak=float(Pm.max()), impulse=float(F.sum()*dt),
                FB=np.array(sim.brake_F), vB=np.array(sim.brake_v), brake=br, settings={k: getattr(cfg, k) for k in DEFAULTS})

def strip(d):
    return {k: v for k, v in d.items() if not isinstance(v, np.ndarray)}

# ---------------------------------------------------------------------------------------------------------------------
# the casing: round 18's hold, round 19's integrator, unchanged
M_CAS = L18["const"]["M_CAS"]; K = L18["hold"]["k_ring"]; k_r = L18["hold"]["k_r_hang"]; F_LIM = L18["const"]["F_LIM_MN"]*1e6; F_LINE = L18["const"]["F_LINE_kN"]*1e3; Lr = L18["const"]["Lr"]
OM = math.sqrt(K/M_CAS); C_CRIT = 2*math.sqrt(K*M_CAS)

def rebin(F, dt_in, dt_out):
    n = int(round(dt_out/dt_in)); m = int(math.ceil(len(F)/n)); G = np.zeros(m*n); G[:len(F)] = F
    return G.reshape(m, n).mean(axis=1)

def swing_once(Fc, dt, c, k=None, days=12):
    k = K if k is None else k
    n = int(days*86400/dt); x = v = 0.0; xs = np.empty(n); fm = 0.0; nF = len(Fc)
    for i in range(n):
        Fl = Fc[i] if i < nF else 0.0; Fd = max(-F_LIM, min(F_LIM, -c*v)); a = (Fl + Fd - k*x)/M_CAS; v += a*dt; x += v*dt
        xs[i] = x; fm = max(fm, abs(Fd))
    ax = np.abs(xs); xm = float(ax.max()); late = np.where(ax > 0.01*xm)[0]
    return dict(x_max=xm, peak_h=float(np.argmax(ax))*dt/3600.0, die_h=(late[-1] + 1)*dt/3600.0, F_lim_MN=fm/1e6)

def undamped(Fc, dt, k=None):
    k = K if k is None else k; om = math.sqrt(k/M_CAS)
    x = v = 0.0; xm = 0.0
    for Fl in Fc:
        a = (Fl - k*x)/M_CAS; v += a*dt; x += v*dt; xm = max(xm, abs(x))
    amp = math.hypot(x, v/om)
    return dict(amp_after=amp, x_max_during=xm, x_max=max(xm, amp))

def swing_daily(Fc, dtF, N, c, hold, k=None, days=40, J_exact=None, dt=20.0):
    k = K if k is None else k
    T = 86400.0/N; nF = len(Fc); kmax = int(nF*dtF/T) + 1
    Fm = (N*J_exact/86400.0) if J_exact is not None else float(Fc.sum()*dtF)/T
    n = int(days*86400/dt); x = v = 0.0; xs_ = []; fl_ = []
    for i in range(n):
        t = i*dt; ph = t % T; Fl = 0.0
        for m in range(kmax):
            j = int((ph + m*T)/dtF)
            if j < nF: Fl += Fc[j]
        Fd = -c*v - (Fm if hold else 0.0); Fd = max(-F_LIM, min(F_LIM, Fd)); a = (Fl + Fd - k*x)/M_CAS; v += a*dt; x += v*dt
        if t > (days - 5)*86400: xs_.append(x); fl_.append(abs(Fd))
    xs_ = np.array(xs_)
    return dict(mean_push_MN=Fm/1e6, x_mean=float(xs_.mean()), swing=float((xs_.max() - xs_.min())/2), x_max=float(np.abs(xs_).max()), F_lim_peak_MN=max(fl_)/1e6, F_lim_mean_MN=float(np.mean(fl_))/1e6, overlap=nF*dtF/T)

def derived(o):
    o["line_kN"] = k_r*o["x_max"]/1e3; o["line_share"] = k_r*o["x_max"]/F_LINE; o["lean_deg"] = math.degrees(math.atan(o["x_max"]/700e3)); o["rot_deg"] = 360*o["x_max"]/Lr; return o

def casing_rows(Fl, Fb, daily=True):
    """Fl: thrust of the launch, Fb: braking force, both in 0.1 s steps and both from the simulation.  x is counted positive the way the launch pushes the casing."""
    res = {}
    for nm, F01 in (("not braked", Fl), ("braked", np.concatenate((Fl, -Fb)))):
        F10 = rebin(F01, 0.1, 10.0); r = derived(swing_once(F10, 10.0, C_CRIT)); u = undamped(F10, 10.0)
        r.update(undamped_after=u["amp_after"], undamped_max=u["x_max"], push_s=len(F01)*0.1, impulse=float(F01.sum()*0.1)); res[nm] = r
        P("   %-10s: push lasts %.0f s, net impulse %.4e N s.  LIMs damping: farthest %.0f m, %.2f h after the launch begins ; within 1 %% of that after %.1f h ; LIM push at most %.2f MN of %.1f MN ;"
          % (nm, r["push_s"], r["impulse"], r["x_max"], r["peak_h"], r["die_h"], r["F_lim_MN"], F_LIM/1e6))
        P("   %s  reaction line %.3f kN (%.2f %% of %.1f kN) ; lean %.4f deg ; casing turned %.5f deg.  No damping: farthest %.0f m, left swinging %.0f m." % (" "*10, r["line_kN"], 100*r["line_share"], F_LINE/1e3, r["lean_deg"], r["rot_deg"], r["undamped_max"], r["undamped_after"]))
    if daily:
        res["daily"] = []
        P("     launches | sled       | mean push | under way at once | LIMs                              | casing sits at | swings about it | farthest | reaction line there | LIM push: peak, mean")
        for N in (7, 14):
            for nm, F01 in (("braked", np.concatenate((Fl, -Fb))), ("not braked", Fl)):
                F10 = rebin(F01, 0.1, 10.0)
                for cn, hold in (("damping, critical", False), ("damping and holding the mean push", True)):
                    r = derived(swing_daily(F10, 10.0, N, C_CRIT, hold, J_exact=float(F01.sum()*0.1))); r.update(N=N, sled=nm, lims=cn); res["daily"].append(r)
                    P("     %2d a day | %-10s | %6.3f MN | %4.2f              | %-33s | %7.0f m      | %6.0f m        | %6.0f m | %5.2f kN (%4.1f %%)    | %5.1f, %5.1f MN"
                      % (N, nm, r["mean_push_MN"], r["overlap"], cn, r["x_mean"], r["swing"], r["x_max"], r["line_kN"], 100*r["line_share"], r["F_lim_peak_MN"], r["F_lim_mean_MN"]))
    return res

# ---------------------------------------------------------------------------------------------------------------------
P("Round 20.  The braking run, the felt load in space, and the return launch of empty pods.")
P("\n1. The design launch and the braking of its sled, lsm_simulation.run_lsm_simulation(brake=True)  ('python lsm_simulation.py --brake')")
D = lsm_run("design_brake"); laps = D["x_end"]/cfg.L_RING; B = D["brake"]
P("   Launch: %.1f minutes ; %.1f km, %.2f laps ; peak thrust %.2f MN ; peak power %.2f GW ; impulse %.4e N s.  The book: 190 minutes, five laps, 83.3 MN, 416 GW." % (D["dur"]/60, D["x_end"]/1e3, laps, D["F_peak"]/1e6, D["P_peak"]/1e9, D["impulse"]))
if not (abs(D["dur"]/60 - 190) < 0.5 and abs(laps - 5.0) < 0.05 and abs(D["F_peak"]/1e6 - 83.3) < 0.05 and abs(D["P_peak"]/1e9 - 416) < 0.5):
    P("   STOP: the simulated launch is not the book's."); sys.exit(1)
o19 = L19["design"]
assert abs(D["dur"] - o19["dur"]) < 1e-6 and abs(D["impulse"] - o19["impulse"]) < 1.0 and abs(D["x_end"] - o19["x_end"]) < 1e-3
P("   Reproduced, and the thrust is round 19's to the last digit (the felt-load change does not touch the design launch's thrust).")
P("   Braking (simulated): the %.0f t sled from %.1f m/s ; %.1f s, %.1f minutes ; %.1f km, %.2f laps ; force %.2f MN at the release, rising as 1 / speed to %.2f MN when the sled is down to %.0f km/s"
  % (B["m_sled"]/1e3, B["v0"], B["duration_s"], B["duration_s"]/60, B["distance_m"]/1e3, B["distance_m"]/cfg.L_RING, B["F_first"]/1e6, B["F_peak"]/1e6, cfg.V_CORNER/1e3))
P("     (%.1f minutes after the release), then held there for the last %.1f minutes.  Peak power %.2f GW ; peak coil voltage %.1f kV (limit %.0f kV) ; peak deceleration %.3f g ; peak load on the frame %.3f g (limit %.1f g)."
  % (B["t_corner_s"]/60, (B["duration_s"] - B["t_corner_s"])/60, B["P_peak"]/1e9, B["V_peak"]/1e3, cfg.V_COIL_LIMIT/1e3, B["decel_peak_g"], B["g_peak"], cfg.G_FRAME_MAX))
P("     Limits met: voltage %s, frame load %s, HVDC %s." % (B["voltage_limited"], B["frame_limited"], B["hvdc_limited"]))
mins = (0, 15, 30, 45, 60, 75, 90)
P("     Braking force at %s minutes after the release: %s MN ; at the end %.1f MN." % (", ".join(str(m) for m in mins), ", ".join("%.1f" % (D["FB"][int(m*600)]/1e6) for m in mins), D["FB"][-2]/1e6))
P("     Energy the sled gives up %.1f TJ (the book: 2,360 TJ) ; heat in the stator's tape %.3f TJ, %.3f of the %.3f TJ the launch left there (the book: about half as much again) ;"
  % (B["E_mech"]/1e12, B["E_heat"]/1e12, B["E_heat"]/D["E_heat"], D["E_heat"]/1e12))
P("     returned at the coils %.1f TJ.  Inverter losses are not in lsm_simulation.py, for the launch or the braking, so the book's 'a few percent in the inverters' is not checked here." % (B["E_returned"]/1e12))
P("     Impulse given back %.4e N s (sled mass x release speed %.4e) ; net impulse left on the casing %.4e N s." % (B["impulse"], B["m_sled"]*B["v0"], D["impulse"] - B["impulse"]))
b19 = L19["brake"]
P("     Round 19's built braking: %.1f s, %.4e N s.  The simulated braking: %.1f s, %.4e N s." % (b19["dur_s"], b19["impulse"], B["duration_s"], B["impulse"]))
P("     Launch and braking together keep the track for %.1f minutes (%.2f h) ; one sled at a time, that is %.2f launches a day." % ((D["dur"] + B["duration_s"])/60, (D["dur"] + B["duration_s"])/3600, 86400/(D["dur"] + B["duration_s"])))
OUT["design"] = strip(D); OUT["design"].update(laps=laps, brake_F_MN_every_15_min=[float(D["FB"][int(m*600)]/1e6) for m in mins])
# the same braking with the launch's 1 g cap on the acceleration along the track (not the ruling ; for the record)
D1 = lsm_run("design_brake_cap1g", A_BRAKE_MAX_G=1.0); B1 = D1["brake"]
P("   For the record, the same braking with the deceleration capped at 1 g (the launch's default cap along the track ; not the ruling): %.1f minutes, peak force %.2f MN." % (B1["duration_s"]/60, B1["F_peak"]/1e6))
OUT["design_brake_cap1g"] = strip(D1)

# ---------------------------------------------------------------------------------------------------------------------
P("\n2. The felt load, at the speed in space")
R = cfg.R_ORBIT; v_rot = cfg.V_CASING_SPACE
P("   lsm_simulation.py now takes the centrifugal term at the speed over the casing plus the casing's own %.1f m/s (ring_altitude.V_GROUND_SYNC)." % v_rot)
P("   Design launch: release %.1f m/s over the casing, %.1f m/s in space ; net outward there %.3f g ; the simulation's peak load (vector sum with the %.3f g along the track) %.3f g.  Before the change it printed %.3f g."
  % (D["v_end"], D["v_space"], D["a_out_g"], D["F_end"]/D["m_total"]/cfg.G_0, D["g_peak"], L19["design"]["g_peak"]))
P("   The book's 12.6 g for the design launch stands ; it was already worked out at 30.5 km/s in space.")
v3_space = math.sqrt((3*cfg.G_0 + cfg.G_250)*R); v3_cas = v3_space - v_rot
P("   Crewed limit: a net outward acceleration of 3 g is reached at %.1f m/s in space, %.1f m/s over the casing.  lsm_simulation.py --v_launch=%.1f --g-limit=3 --brake" % (v3_space, v3_cas, v3_cas))
Cw = lsm_run("crewed", V_LAUNCH=v3_cas, G_LOAD_MAX=3.0); c19 = L19["crewed"]
i_thr = np.where(Cw["g"] >= 3.0*0.999999)[0]; t_thr = (len(Cw["F"]) - i_thr[0])*Cw["dt"] if len(i_thr) else 0.0
P("     %.0f t in all ; %.2f minutes (round 19: %.2f) ; %.1f km, %.2f laps ; peak thrust %.2f MN ; impulse %.4e N s (round 19: %.4e) ; release %.1f m/s over the casing, %.1f m/s in space ;"
  % (Cw["m_total"]/1e3, Cw["dur"]/60, c19["dur"]/60, Cw["x_end"]/1e3, Cw["x_end"]/cfg.L_RING, Cw["F_peak"]/1e6, Cw["impulse"], c19["impulse"], Cw["v_end"], Cw["v_space"]))
P("     net outward at release %.3f g ; peak load %.3f g (round 19's run read %.3f g)." % (Cw["a_out_g"], Cw["g_peak"], c19["g_peak"]))
P("     The simulation's limit acts on the vector sum of the outward load and the acceleration along the track.  So the 3 g limit now binds: from %.1f m/s over the casing the thrust is throttled, for the last %.1f s,"
  % (Cw["g_engage_v"] if Cw["g_engage_v"] else float("nan"), t_thr))
P("     from %.2f MN to %.2f MN at the release.  In round 19 it did not bind and the thrust at the release was %.2f MN." % (Cw["F"][i_thr[0] - 1]/1e6 if len(i_thr) and i_thr[0] > 0 else float("nan"), Cw["F_end"]/1e6, c19["F_end_MN"]))
v_esc = math.sqrt(2*e1.GM/R)
P("     Hyperbolic excess at %.1f m/s in space, escape speed %.1f m/s: %.2f km/s (Table 10.1: 12.3 km/s)." % (Cw["v_space"], v_esc, math.sqrt(Cw["v_space"]**2 - v_esc**2)/1e3))
Bc = Cw["brake"]
P("     Braking of its sled (simulated): %.1f minutes ; impulse given back %.4e N s ; %.1f TJ." % (Bc["duration_s"]/60, Bc["impulse"], Bc["E_mech"]/1e12))
OUT["crewed"] = strip(Cw); OUT["crewed"].update(v3_space=v3_space, v3_casing=v3_cas, throttled_s=t_thr, v_escape=v_esc)

# ---------------------------------------------------------------------------------------------------------------------
P("\n3. The return launch of empty pods to the platform")
rates = e1.load_rates(); T_sh = rates["T"]
d4 = e1.descent(4000.0, 1000.0); _p, _e, tot = e1.dose_per_descent(rates, d4); x_sh = e1.shield_for(T_sh, tot, 2.0)
rows = e2.arrangements(); mag = dict(min([r for r in rows if r["kind"] == "facing pair"], key=lambda r: r["kg_per_N"]))
g = e2.size_pod2(100, x_sh, mag); m_pod = g["m"]; m_pax = g["m_pax"]; m_empty = m_pod - m_pax
P("   The 78 t pod of Chapter 13 is elevator_system2.size_pod2(100 seats, the 2 mSv shield at 4,000 m/s, the facing-pair magnets) ; elevator_system.py by itself stops at 3,000 m/s.")
P("   It counts: seats and passengers %d x %.0f kg = %.2f t ; shield %.2f t ; hull and skin %.2f t ; magnets and carriage %.2f t ; loop, pumps, wheels, avionics %.2f t ; in all %.2f t."
  % (g["N"], e1.M_SEAT, m_pax/1e3, g["m_sh"]/1e3, g["m_hull"]/1e3, g["m_car"]/1e3, e1.M_FIXED/1e3, m_pod/1e3))
P("   Empty pod = %.2f t less the %.2f t = %.2f t.  The script's 130 kg is one figure for the passenger, the seat and the supplies, so the seats go out with the passengers here ; the script has no separate seat mass."
  % (m_pod/1e3, m_pax/1e3, m_empty/1e3))
P("   Left out: any transfer craft, engine or propellant for the %.2f km/s burn at the platform (Table 13.4), and any cradle that holds the pods on the sled.  Neither is in a script." % 1.43)
r_p = R; r_a = e1.r_geo; v_h = math.sqrt(e1.GM*(2/r_p - 2/(r_p + r_a))); v_h_cas = v_h - v_rot
P("   Hohmann transfer from the ring to the platform (%.0f km radius): %.1f m/s in space, %.1f m/s over the casing (Table 13.4: 9.82 and 9.30 km/s)." % (r_a/1e3, v_h, v_h_cas))
assert abs(v_h/1e3 - 9.82) < 0.005 and abs(v_h_cas/1e3 - 9.30) < 0.005
P("   The sled: lsm_simulation.py has one sled, 525 kg/m, with its length an option.  The 10 km, 5,250 t sled is used.  The thrust is the limit: --accel=2 only lifts the default 1 g cap along the track,")
P("   as Chapter 10 has it for a light load ('the same 83.3 MN on the 5,250 tonne sled alone is 1.6 g').  The runs with the 1 g cap kept are below them, for the record.")
PODS = {}
for n, tag in ((1, "pod_1"), (10, "pod_10")):
    for cap, sfx in ((2.0, ""), (1.0, "_cap1g")):
        c = lsm_run(tag + sfx, M_SPACECRAFT=n*m_empty, V_LAUNCH=v_h_cas, A_MAX_G=cap, A_BRAKE_MAX_G=(None if cap > 1.5 else 1.0)); b = c["brake"]; PODS[tag + sfx] = c
        P("   %s%2d empty pod%s, %.1f t on the sled, %.1f t in all%s: %.2f minutes ; %.0f km (%.3f laps) ; peak thrust %.2f MN (%.3f g along the track) ; thrust at release %.2f MN ; peak power %.2f GW ; impulse %.4e N s ;"
          % ("" if not sfx else "  [1 g cap] ", n, "" if n == 1 else "s", n*m_empty/1e3, c["m_total"]/1e3, "" if not sfx else "", c["dur"]/60, c["x_end"]/1e3, c["x_end"]/cfg.L_RING, c["F_peak"]/1e6, c["F_peak"]/c["m_total"]/cfg.G_0, c["F_end"]/1e6, c["P_peak"]/1e9, c["impulse"]))
        P("      release %.1f m/s over the casing, %.1f m/s in space ; net outward there %.3f g ; peak felt load %.3f g (at the start, the %.3f g along the track with the %.3f g of weight on the casing) ;"
          % (c["v_end"], c["v_space"], c["a_out_g"], c["g_peak"], c["F"][0]/c["m_total"]/cfg.G_0, -c["ar"][0]/cfg.G_0))
        P("      braking of the sled: %.2f minutes, peak %.2f MN, %.3f g, impulse given back %.4e N s, %.1f TJ." % (b["duration_s"]/60, b["F_peak"]/1e6, b["decel_peak_g"], b["impulse"], b["E_mech"]/1e12))
        OUT[tag + sfx] = strip(c)
OUT["pod"] = dict(m_pod=m_pod, m_pax=m_pax, m_empty=m_empty, m_shield=g["m_sh"], m_hull=g["m_hull"], m_car=g["m_car"], m_fixed=e1.M_FIXED, shield_g_cm2=x_sh, v_hohmann=v_h, v_hohmann_casing=v_h_cas, r_geo=r_a)

# ---------------------------------------------------------------------------------------------------------------------
P("\n4. The casing's swing.  Hold %.1f N/m ; casing %.4e kg ; period %.0f s (%.2f h) ; LIMs damping at the critical rate %.3e N s/m, limited to %.1f MN.  Launch and braking both from the simulation." % (K, M_CAS, 2*math.pi/OM, 2*math.pi/OM/3600, C_CRIT, F_LIM/1e6))
P("   The design launch:")
OUT["design"]["casing"] = casing_rows(D["F"], D["FB"])
for nm, key, c in (("The crewed limit", "crewed", Cw), ("Ten empty pods", "pod_10", PODS["pod_10"]), ("One empty pod", "pod_1", PODS["pod_1"]), ("Ten empty pods, 1 g cap", "pod_10_cap1g", PODS["pod_10_cap1g"]), ("One empty pod, 1 g cap", "pod_1_cap1g", PODS["pod_1_cap1g"])):
    P("   %s:" % nm); OUT[key]["casing"] = casing_rows(c["F"], c["FB"], daily=False)

P("\n   Round 19 (braking built by hand) against round 20 (braking simulated), the printed numbers of Section 11.6")
c19 = L19["design"]["casing"]; c20 = OUT["design"]["casing"]; cmp_ = []
def row(what, a, b, fmt):
    cmp_.append(dict(what=what, round19=a, round20=b, p19=fmt % a, p20=fmt % b)); P("   %-66s | %-12s | %-12s | %s" % (what, fmt % a, fmt % b, "moves" if fmt % a != fmt % b else "same"))
P("   %-66s | %-12s | %-12s |" % ("", "round 19", "round 20"))
row("Swing with no damping, not braked (km)", c19["not braked"]["undamped_after"]/1e3, c20["not braked"]["undamped_after"]/1e3, "%.1f")
row("Swing with no damping, braked (km)", c19["braked"]["undamped_after"]/1e3, c20["braked"]["undamped_after"]/1e3, "%.1f")
for nm in ("not braked", "braked"):
    a = c19[nm]; b = c20[nm]
    row("Farthest, LIMs damping, %s (km)" % nm, a["x_max"]/1e3, b["x_max"]/1e3, "%.1f"); row("  when, after the launch begins (h)", a["peak_h"], b["peak_h"], "%.1f"); row("  within 1 % after (h)", a["die_h"], b["die_h"], "%.0f")
    row("  LIMs' peak push (MN)", a["F_lim_MN"], b["F_lim_MN"], "%.0f"); row("  reaction line (kN)", a["line_kN"], b["line_kN"], "%.1f"); row("  lean of the lines (deg)", a["lean_deg"], b["lean_deg"], "%.2f"); row("  rotation of the casing (deg)", a["rot_deg"], b["rot_deg"], "%.3f")
for N in (7, 14):
    for nm in ("braked", "not braked"):
        a = [q for q in c19["daily"] if q["N"] == N and q["sled"] == nm and q["lims"].startswith("damping, critical")][0]; b = [q for q in c20["daily"] if q["N"] == N and q["sled"] == nm and q["lims"].startswith("damping, critical")][0]
        row("%d a day, %s: mean push (MN)" % (N, nm), a["mean_push_MN"], b["mean_push_MN"], "%.1f"); row("  casing sits at (km)", a["x_mean"]/1e3, b["x_mean"]/1e3, "%.1f"); row("  swings about it (m)", a["swing"], b["swing"], "%.0f")
        row("  reaction line at the farthest (kN)", a["line_kN"], b["line_kN"], "%.1f"); row("  LIMs' peak push (MN)", a["F_lim_peak_MN"], b["F_lim_peak_MN"], "%.1f")
a = L19["crewed"]["casing"]; b = OUT["crewed"]["casing"]
row("Crewed limit: farthest, not braked (m)", a["not braked"]["x_max"], b["not braked"]["x_max"], "%.0f"); row("Crewed limit: farthest, braked (m)", a["braked"]["x_max"], b["braked"]["x_max"], "%.0f")
OUT["compare"] = cmp_
OUT["const"] = dict(K=K, M_CAS=M_CAS, period_s=2*math.pi/OM, c_crit=C_CRIT, F_LIM_MN=F_LIM/1e6, F_LINE_kN=F_LINE/1e3, k_r=k_r, v_casing_space=v_rot)
P("\n   run time %.0f s" % (time.time() - T0))
with open(os.path.join(LOG_DIR, "round20_brake_felt_pods.txt"), "w", encoding="utf-8") as fh: fh.write(buf.getvalue())
with open(os.path.join(LOG_DIR, "round20_brake_felt_pods.json"), "w", encoding="utf-8") as fh: json.dump(OUT, fh, indent=1, default=float)
print("ok")
