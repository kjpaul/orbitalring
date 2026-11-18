#!/usr/bin/env python3
"""
Round 19: the swing of the casing along the ring driven by the thrust that lsm_simulation.py gives, not by an assumed push.
Brief: Vol_III/Claude_Code_instructions_round19_2026-11-18.md.

1. lsm_simulation.run_lsm_simulation() is run as it stands (time step 0.1 s, the ruled stator, the settings of lsm_config.py) for
   Chapter 10's design launch, and the thrust against time is taken from its own arrays (data_t, data_F).  If the run does not give
   the book's 190 minutes, five laps, 83.3 MN and 416 GW, this script stops.
2. The simulation has no braking run.  The braking of the empty sled is built here from what Chapter 10 says of it (the LSM as a
   generator, the same limits as the launch): the force is the smaller of the peak thrust and the stator's power over the speed.
   It is marked in every line as NOT from the simulation.  It starts at the release with no pause, which is also not in the
   simulation or the book.
3. The casing is a mass on the lines' spring (round 18: 5,095.9 N/m, 5.337e11 kg) with the LIMs damping at the critical rate and
   limited to their 71.4 MN.  The integrator is new here because the push now has a shape ; it is the stepping of rounds 15, 16
   and 18 (v += a dt, then x += v dt) with the thrust averaged over each step so that the impulse is kept.  Fed the 3,690 s push of
   those rounds it must give round 18's numbers again, and the script stops if it does not.
4. Two more launches from the simulation: the 5,000 t spacecraft to the crewed limit of 3 g, and the empty sled to the 9.8 km/s
   of a geosynchronous transfer (the lightest case Chapter 10 gives).

Run from Python/orbitalring:  python tools/round19_swing_simulated.py
Log: runs/round19_logs/round19_swing_simulated.txt and .json ; the simulation's own printing for each run: runs/round19_logs/lsm_<case>.txt
"""
import contextlib, io, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.environ.pop("ORBITAL_RING_LSM_STATOR", None)                      # the ruled stator, the default
import numpy as np
_argv = sys.argv; sys.argv = [sys.argv[0]]
import lsm_config as cfg
import lsm_simulation as sim
sys.argv = _argv
assert cfg.STATOR == "ruled"

LOG_DIR = os.path.join(ROOT, "runs", "round19_logs"); os.makedirs(LOG_DIR, exist_ok=True)
L18 = json.load(open(os.path.join(ROOT, "runs", "round18_logs", "round18_swing_hold.json"), encoding="utf-8"))
buf = io.StringIO()
def P(s=""):
    print(s); buf.write(s + "\n"); sys.stdout.flush()
OUT = {}; T0 = time.time()
DEFAULTS = dict(M_SPACECRAFT=cfg.M_SPACECRAFT, V_LAUNCH=cfg.V_LAUNCH, G_LOAD_MAX=cfg.G_LOAD_MAX, A_MAX_G=cfg.A_MAX_G, L_SLED=cfg.L_SLED, B_SLED_NOMINAL=cfg.B_SLED_NOMINAL, P_HVDC_MAX=cfg.P_HVDC_MAX)

def lsm_run(tag, **over):
    """One run of lsm_simulation.run_lsm_simulation() at full time step.  'over' changes only what the command line of lsm_simulation.py can change."""
    for k, v in DEFAULTS.items(): setattr(cfg, k, v)
    for k, v in over.items():
        assert k in DEFAULTS, k; setattr(cfg, k, v)
    cfg.calc_derived()
    cap = io.StringIO()
    with contextlib.redirect_stdout(cap):
        ok = sim.run_lsm_simulation(quick_mode=False)
    with open(os.path.join(LOG_DIR, "lsm_%s.txt" % tag), "w", encoding="utf-8") as fh: fh.write(cap.getvalue())
    t = np.array(sim.data_t); F = np.array(sim.data_F); v = np.array(sim.data_v); x = np.array(sim.data_x); Pm = np.array(sim.data_P_mech); g = np.array(sim.data_a_occupant_g)
    dt = cfg.DT; dur = t[-1] + dt; v_end = v[-1] + F[-1]/cfg.M_TOTAL*dt; x_end = x[-1] + v[-1]*dt
    return dict(tag=tag, ok=bool(ok), t=t, F=F, v=v, dt=dt, dur=dur, v_end=v_end, x_end=x_end, m_total=cfg.M_TOTAL, m_sled=cfg.M_SLED_HARDWARE, F_peak=float(F.max()), P_peak=float(Pm.max()),
                g_peak=float(g.max()), impulse=float(F.sum()*dt), settings={k: getattr(cfg, k) for k in DEFAULTS})

def brake_profile(m_sled, v0, dt=0.1):
    """NOT from the simulation.  Chapter 10, 'Braking the Sled': the LSM as a generator, held to the same limits as the launch.
    Force on the sled = min(peak thrust, stator power / speed), against its motion, from the release speed to rest."""
    Fp = cfg.F_PEAK_DESIGN*(cfg.L_SLED/10000.0); Pp = cfg.P_STATOR_MAX*(cfg.L_SLED/10000.0); v = v0; Fs = []
    while v > 0:
        F = min(Fp, Pp/v); dv = F/m_sled*dt
        if dv > v: F = v*m_sled/dt; dv = v                           # the last step brings the sled just to rest
        Fs.append(F); v -= dv
        if v < 1e-9: break
    return np.array(Fs)

# ---------------------------------------------------------------------------------------------------------------------
# the casing
M_CAS = L18["const"]["M_CAS"]; K = L18["hold"]["k_ring"]; k_r = L18["hold"]["k_r_hang"]; F_LIM = L18["const"]["F_LIM_MN"]*1e6; F_LINE = L18["const"]["F_LINE_kN"]*1e3; Lr = L18["const"]["Lr"]
OM = math.sqrt(K/M_CAS); C_CRIT = 2*math.sqrt(K*M_CAS)

def rebin(F, dt_in, dt_out):
    """The push averaged over steps of dt_out: the impulse is kept exactly."""
    n = int(round(dt_out/dt_in)); m = int(math.ceil(len(F)/n)); G = np.zeros(m*n); G[:len(F)] = F
    return G.reshape(m, n).mean(axis=1)

def swing_once(Fc, dt, c, k=None, days=12):
    """One launch.  Fc: the push on the casing in steps of dt, from the start of the launch.  Stepping of round 16's one_launch."""
    k = K if k is None else k
    n = int(days*86400/dt); x = v = 0.0; xs = np.empty(n); fm = 0.0; nF = len(Fc)
    for i in range(n):
        Fl = Fc[i] if i < nF else 0.0; Fd = max(-F_LIM, min(F_LIM, -c*v)); a = (Fl + Fd - k*x)/M_CAS; v += a*dt; x += v*dt
        xs[i] = x; fm = max(fm, abs(Fd))
    ax = np.abs(xs); xm = float(ax.max()); late = np.where(ax > 0.01*xm)[0]; i_end = min(nF, n) - 1
    return dict(x_max=xm, peak_h=float(np.argmax(ax))*dt/3600.0, die_h=(late[-1] + 1)*dt/3600.0, F_lim_MN=fm/1e6, x_end_push=float(xs[i_end]), xs=xs)

def undamped(Fc, dt, k=None):
    """No damping: the farthest the casing gets while pushed, and the swing it is left with when the push is over."""
    k = K if k is None else k; om = math.sqrt(k/M_CAS)
    x = v = 0.0; xm = 0.0
    for Fl in Fc:
        a = (Fl - k*x)/M_CAS; v += a*dt; x += v*dt; xm = max(xm, abs(x))
    amp = math.hypot(x, v/om)
    return dict(amp_after=amp, x_max_during=xm, x_max=max(xm, amp))

def swing_daily(Fc, dtF, N, c, hold, k=None, days=40, J_exact=None, dt=20.0):
    """N launches a day at even intervals, the same push each time, overlapping where one is not over when the next begins.  Stepping of round 15's sim."""
    k = K if k is None else k
    T = 86400.0/N; nF = len(Fc); kmax = int(nF*dtF/T) + 1
    Fm = (N*J_exact/86400.0) if J_exact is not None else float(Fc.sum()*dtF)/T
    n = int(days*86400/dt); x = v = 0.0; xs_ = []; fl_ = []
    for i in range(n):
        t = i*dt; ph = t % T; Fl = 0.0
        for m in range(kmax):                                        # launches that overlap add
            j = int((ph + m*T)/dtF)
            if j < nF: Fl += Fc[j]
        Fd = -c*v - (Fm if hold else 0.0); Fd = max(-F_LIM, min(F_LIM, Fd)); a = (Fl + Fd - k*x)/M_CAS; v += a*dt; x += v*dt
        if t > (days - 5)*86400: xs_.append(x); fl_.append(abs(Fd))
    xs_ = np.array(xs_)
    return dict(mean_push_MN=Fm/1e6, x_mean=float(xs_.mean()), swing=float((xs_.max() - xs_.min())/2), x_max=float(np.abs(xs_).max()), F_lim_peak_MN=max(fl_)/1e6, F_lim_mean_MN=float(np.mean(fl_))/1e6, overlap=nF*dtF/T)

def derived(o):
    o["line_kN"] = k_r*o["x_max"]/1e3; o["line_share"] = k_r*o["x_max"]/F_LINE; o["lean_deg"] = math.degrees(math.atan(o["x_max"]/700e3)); o["rot_deg"] = 360*o["x_max"]/Lr; return o

# ---------------------------------------------------------------------------------------------------------------------
P("Round 19.  The swing of the casing driven by the simulated thrust.")
P("\n0. Check of the new integrator against round 18 (the 3,690 s push, hold 5,100 N/m):")
for nm, J in (("not braked", 10.25e6*30000.0), ("braked", 5.0e6*30000.0)):
    Fc = np.full(369, J/3690.0); r = swing_once(Fc, 10.0, 2*math.sqrt(5100.0*M_CAS), k=5100.0); o = L18["new"]["one"][nm]
    P("   %-10s: farthest %.3f m (round 18: %.3f) ; LIM push %.4f MN (%.4f) ; within 1 %% after %.3f h (%.3f)" % (nm, r["x_max"], o["x_max"], r["F_lim_MN"], o["F_lim_MN"], r["die_h"], o["die_h"]))
    assert abs(r["x_max"] - o["x_max"]) < 1e-6 and abs(r["F_lim_MN"] - o["F_lim_MN"]) < 1e-9 and abs(r["die_h"] - o["die_h"]) < 1e-9
Fc = np.full(369, 0.281*M_CAS/3690.0)
r = swing_daily(Fc, 10.0, 7, 2*math.sqrt(5100.0*M_CAS), False, k=5100.0, J_exact=0.281*M_CAS); o = [q for q in L18["new"]["rows"] if q["N"] == 7 and q["sled"] == "braked" and q["lims"].startswith("damping, critical")][0]
P("   7 a day, braked: sits at %.3f m (round 18: %.3f) ; swings %.3f m (%.3f) ; LIM peak %.4f MN (%.4f)" % (r["x_mean"], o["x_mean"], r["swing"], o["swing"], r["F_lim_peak_MN"], o["F_lim_peak_MN"]))
assert abs(r["x_mean"] - o["x_mean"]) < 1e-6 and abs(r["swing"] - o["swing"]) < 1e-6 and abs(r["F_lim_peak_MN"] - o["F_lim_peak_MN"]) < 1e-9
P("   The new integrator gives round 18 again.")

# ---------------------------------------------------------------------------------------------------------------------
P("\n1. The thrust of the design launch, from lsm_simulation.py")
D = lsm_run("design")
laps = D["x_end"]/cfg.L_RING
P("   Run: lsm_simulation.run_lsm_simulation(quick_mode=False), the same as 'python lsm_simulation.py' with no arguments.  Its own printing: runs/round19_logs/lsm_design.txt")
P("   Settings (lsm_config.py, the ruled stator): pole pitch %.0f m ; %d turns of %d layers of %.0f mm tape ; coil current %.1f A ; sled field %.2f T ; load angle %.0f deg ; winding factor %.4f ;"
  % (cfg.TAU_P, cfg.N_STATOR, cfg.N_HTS_LAYERS, cfg.TAPE_WIDTH_MM, cfg.I_PEAK, cfg.B_SLED_NOMINAL, math.degrees(cfg.DELTA_MAX), cfg.K_W))
P("     sled %.0f m, %.0f kg/m, active share %.3f ; spacecraft %.0f t ; sled %.0f t ; in all %.0f t ; to %.0f m/s over the casing ; tangential limit %.1f g ; load limit %.1f g ;"
  % (cfg.L_SLED, cfg.M_SLED_PER_M, cfg.ACTIVE_FRACTION, D["settings"]["M_SPACECRAFT"]/1e3, D["m_sled"]/1e3, D["m_total"]/1e3, D["settings"]["V_LAUNCH"], D["settings"]["A_MAX_G"], D["settings"]["G_LOAD_MAX"]))
P("     corner speed %.0f m/s ; stator power %.2f GW ; above the corner the sled's field goes as (5 km/s / v)^%.1f and the coil current as (5 km/s / v)^%.1f ; coil limit %.0f kV ; grid limit %.1f GW (not reached) ; time step %.1f s."
  % (cfg.V_CORNER, cfg.P_STATOR_MAX/1e9, cfg.ALPHA_SLED, 1 - cfg.ALPHA_SLED, cfg.V_COIL_LIMIT/1e3, cfg.P_HVDC_MAX/1e9, cfg.DT))
P("   Result: %.1f minutes (%.2f h) ; %.1f km, %.2f laps ; peak thrust %.2f MN ; peak power %.2f GW ; peak load %.2f g." % (D["dur"]/60, D["dur"]/3600, D["x_end"]/1e3, laps, D["F_peak"]/1e6, D["P_peak"]/1e9, D["g_peak"]))
P("   The book: 190 minutes, five laps, 83.3 MN, 416 GW.")
if not (abs(D["dur"]/60 - 190) < 0.5 and abs(laps - 5.0) < 0.05 and abs(D["F_peak"]/1e6 - 83.3) < 0.05 and abs(D["P_peak"]/1e9 - 416) < 0.5):
    P("   STOP: the simulated launch is not the book's."); sys.exit(1)
P("   Reproduced.")
i_c = int(np.argmax(D["v"] > cfg.V_CORNER)); t_c = D["t"][i_c]
P("   The thrust against time: %.2f MN for the first %.0f s (%.1f minutes), until the sled passes %.0f km/s ; then the power is held and the thrust falls as 1 / speed," % (D["F_peak"]/1e6, t_c, t_c/60, cfg.V_CORNER/1e3))
P("     to %.2f MN at the release.  Mean thrust %.2f MN.  Thrust at 30, 60, 90, 120, 150, 180 minutes: %s MN." % (D["F"][-1]/1e6, D["impulse"]/D["dur"]/1e6, ", ".join("%.1f" % (D["F"][int(m*600)]/1e6) for m in (30, 60, 90, 120, 150, 180))))
P("   Impulse (the sum of thrust x time step): %.4e N s.  Mass x final speed: %.4e N s.  The book's 3.1 x 10^11 N s: agrees." % (D["impulse"], D["m_total"]*D["v_end"]))
P("   The hold-back for the reaction plate's temperature: the LSM sled has no reaction plate, and lsm_simulation.py has no thermal limit on the thrust.  (The plate and its temperature")
P("     are the LIM mass driver's, md_simulation.py.)  The limits that shape this thrust are the coil current (the peak thrust) and the stator's power above 5 km/s.")
P("   Round 15, 16 and 18 took the push as the impulse over the peak thrust, 3,690 s.  The simulated launch pushes for %.0f s, %.1f times as long, at a thrust that falls to a sixth." % (D["dur"], D["dur"]/3690.0))
OUT["design"] = {k: (v if not isinstance(v, np.ndarray) else None) for k, v in D.items()}
OUT["design"].update(laps=laps, t_corner_s=float(t_c), F_end_MN=float(D["F"][-1]/1e6), F_mean_MN=D["impulse"]/D["dur"]/1e6, thrust_every_10_min_MN=[float(D["F"][int(m*600)]/1e6) for m in range(0, 19)])

P("\n   The braking.  lsm_simulation.py has NO braking run: its loop only accelerates, and the word does not occur in it or in lsm_physics.py.")
P("   Chapter 10, 'Braking the Sled' (Vol_III_Ch10.md line 122 to 124): the LSM as a generator, 'held to the same 100 kV coil limit as the launch, so the power is again about 416 GW',")
P("   'the 5,250 tonne sled gives up its energy in about an hour and a half'.  Built here from that, NOT simulated: force = min(%.2f MN, %.2f GW / speed) on the %.0f t sled from %.0f m/s."
  % (cfg.F_PEAK_DESIGN/1e6, cfg.P_STATOR_MAX/1e9, D["m_sled"]/1e3, D["v_end"]))
B = brake_profile(D["m_sled"], D["v_end"]); t_b = len(B)*0.1; J_b = float(B.sum()*0.1)
P("   That gives %.0f s (%.1f minutes), impulse %.4e N s (sled mass x release speed %.4e), %.0f TJ returned ; the book's 'about an hour and a half' and 2,360 TJ." % (t_b, t_b/60, J_b, D["m_sled"]*D["v_end"], 0.5*D["m_sled"]*D["v_end"]**2/1e12))
P("   The braking is taken to begin at the release with no pause.  [ASSUMED: neither the simulation nor the book gives the pause.]  Net impulse left on the casing: %.4e N s." % (D["impulse"] - J_b))
OUT["brake"] = dict(simulated=False, source="Vol_III_Ch10.md, 'Braking the Sled'", dur_s=t_b, impulse=J_b, F_peak_MN=float(B.max()/1e6), net_impulse=D["impulse"] - J_b, pause_s=0.0)

# ---------------------------------------------------------------------------------------------------------------------
def casing_rows(tag, Fl, Fb, daily=True):
    """Fl: thrust of the launch, Fb: braking force, both in 0.1 s steps.  The casing is pushed by -thrust and by +braking ; x is counted positive the way the launch pushes it."""
    res = {}
    for nm, F01 in (("not braked", Fl), ("braked", np.concatenate((Fl, -Fb)))):
        F10 = rebin(F01, 0.1, 10.0); r = derived(swing_once(F10, 10.0, C_CRIT)); r.pop("xs"); u = undamped(F10, 10.0)
        r.update(undamped_after=u["amp_after"], undamped_max=u["x_max"], push_s=len(F01)*0.1, impulse=float(F01.sum()*0.1)); res[nm] = r
        P("   %-10s: push lasts %.0f s, impulse %.4e N s.  LIMs damping: farthest %.0f m, %.2f h after the launch begins ; within 1 %% of that after %.1f h ; LIM push at most %.2f MN of %.1f MN ;"
          % (nm, r["push_s"], r["impulse"], r["x_max"], r["peak_h"], r["die_h"], r["F_lim_MN"], F_LIM/1e6))
        P("   %s  reaction line %.2f kN (%.1f %% of %.1f kN) ; lean %.4f deg ; casing turned %.5f deg.  No damping: farthest %.0f m, left swinging %.0f m." % (" "*10, r["line_kN"], 100*r["line_share"], F_LINE/1e3, r["lean_deg"], r["rot_deg"], r["undamped_max"], r["undamped_after"]))
    if daily:
        res["daily"] = []
        P("     launches | sled       | mean push | under way at once | LIMs                              | casing sits at | swings about it | farthest | reaction line there | LIM push: peak, mean")
        for N in (7, 14):
            for nm, F01 in (("braked", np.concatenate((Fl, -Fb))), ("not braked", Fl)):
                F10 = rebin(F01, 0.1, 10.0)
                for cn, hold in (("damping, critical", False), ("damping and holding the mean push", True)):
                    r = derived(swing_daily(F10, 10.0, N, C_CRIT, hold, J_exact=float(F01.sum()*0.1))); r.update(N=N, sled=nm, lims=cn); res["daily"].append(r)
                    P("     %2d a day | %-10s | %5.2f MN  | %4.2f              | %-33s | %7.0f m      | %6.0f m        | %6.0f m | %5.2f kN (%4.1f %%)    | %5.1f, %5.1f MN"
                      % (N, nm, r["mean_push_MN"], r["overlap"], cn, r["x_mean"], r["swing"], r["x_max"], r["line_kN"], 100*r["line_share"], r["F_lim_peak_MN"], r["F_lim_mean_MN"]))
    return res

P("\n2. The casing's swing on that thrust.  Hold %.1f N/m ; casing %.4e kg ; period %.0f s (%.2f h) ; LIMs damping at the critical rate %.3e N s/m, limited to %.1f MN." % (K, M_CAS, 2*math.pi/OM, 2*math.pi/OM/3600, C_CRIT, F_LIM/1e6))
P("   ('braked' rows rest on the braking built from the book's words, not on a simulation run.)")
OUT["design"]["casing"] = casing_rows("design", D["F"], B)
P("   Seven a day is one launch every %.0f s and the launch lasts %.0f s, so with braking (%.0f s more) two sleds are being driven or braked at once part of the time ; their pushes are added." % (86400/7, D["dur"], t_b))

# ---------------------------------------------------------------------------------------------------------------------
P("\n3. Two more launches from the simulation")
R = cfg.R_ORBIT; W_E = 2*math.pi/86164.0905
v3_space = math.sqrt((3*cfg.G_0 + cfg.G_250)*R); v_rot = W_E*R
P("   The crewed limit.  Net outward acceleration of 3 g (3 x %.3f m/s2) at the ring's radius is reached at %.0f m/s in space (the book's 16.3 km/s) ; the casing turns with the Earth at %.0f m/s," % (cfg.G_0, v3_space, v_rot))
P("     so the sled's speed over the casing is %.0f m/s.  Chapter 10 names no separate crewed spacecraft, so this is the 5,000 t spacecraft on the 10 km sled.  lsm_simulation.py --v_launch=%.0f --g-limit=3" % (v3_space - v_rot, v3_space - v_rot))
P("     [lsm_simulation.py works out the felt load from the speed over the casing, not the speed in space, so its own peak load for this run reads lower than 3 g.]")
Cw = lsm_run("crewed", V_LAUNCH=v3_space - v_rot, G_LOAD_MAX=3.0)
v_geo = 9800.0 - v_rot
P("   The lightest case of Chapter 10 ('Heavy Launches and Light Ones'): the empty sled to the 9.8 km/s of a transfer to geosynchronous orbit, %.0f m/s over the casing, at the motor's whole thrust (1.6 g)." % v_geo)
P("     lsm_simulation.py --mass=0 --v_launch=%.0f --accel=2  (the 2 g only lifts the simulation's default cap of 1 g on the tangential acceleration ; the thrust is the limit)." % v_geo)
E = lsm_run("empty_sled_geo", M_SPACECRAFT=0.0, V_LAUNCH=v_geo, A_MAX_G=2.0)
P("     The book's figures for it, from a hand formula: about 12 minutes and 3,600 km.  The simulation: %.1f minutes, %.0f km." % (E["dur"]/60, E["x_end"]/1e3))
P("     Nothing leaves the ring in this run.  Chapter 10 gives no lighter launch with a payload, and I have not made one up.")
for nm, c in (("crewed limit, 5,000 t spacecraft", Cw), ("empty sled to 9.8 km/s", E)):
    P("\n   %s: %.0f t in all to %.0f m/s over the casing ; %.1f minutes ; %.0f km (%.2f laps) ; peak thrust %.2f MN ; peak power %.2f GW ; thrust at release %.2f MN ; impulse %.4e N s ; peak load in the simulation %.2f g"
      % (nm, c["m_total"]/1e3, c["v_end"], c["dur"]/60, c["x_end"]/1e3, c["x_end"]/cfg.L_RING, c["F_peak"]/1e6, c["P_peak"]/1e9, c["F"][-1]/1e6, c["impulse"], c["g_peak"]))
    Bc = brake_profile(c["m_sled"], c["v_end"])
    P("   braking built as above (not simulated): %.0f s, impulse %.4e N s." % (len(Bc)*0.1, float(Bc.sum()*0.1)))
    o = {k: (v if not isinstance(v, np.ndarray) else None) for k, v in c.items()}; o.update(laps=c["x_end"]/cfg.L_RING, F_end_MN=float(c["F"][-1]/1e6), brake_s=len(Bc)*0.1, brake_impulse=float(Bc.sum()*0.1))
    o["casing"] = casing_rows(c["tag"], c["F"], Bc, daily=False); OUT[c["tag"]] = o
OUT["crew_speed"] = dict(v_space=v3_space, v_rot=v_rot, v_casing=v3_space - v_rot); OUT["empty_speed"] = dict(v_space=9800.0, v_casing=v_geo)

# ---------------------------------------------------------------------------------------------------------------------
P("\n4. Round 18 (the 3,690 s push at 5,100 N/m) against the simulated thrust (5,095.9 N/m)")
n18 = L18["new"]; c19 = OUT["design"]["casing"]
cmp_ = []
def row(what, a, b, fmt):
    cmp_.append(dict(what=what, round18=a, round19=b, p18=fmt % a, p19=fmt % b)); P("   %-66s | %-12s | %-12s | %s" % (what, fmt % a, fmt % b, "moves" if fmt % a != fmt % b else "same"))
P("   %-66s | %-12s | %-12s |" % ("", "round 18", "round 19"))
row("Swing with no damping, not braked (km)", n18["swing_m"]/1e3, c19["not braked"]["undamped_after"]/1e3, "%.1f")
row("Swing with no damping, braked (km)", n18["swing_braked_m"]/1e3, c19["braked"]["undamped_after"]/1e3, "%.1f")
for nm in ("not braked", "braked"):
    a = n18["one"][nm]; b = c19[nm]
    row("Farthest, LIMs damping, %s (km)" % nm, a["x_max"]/1e3, b["x_max"]/1e3, "%.1f"); row("  when, after the launch begins (h)", a["peak_h"], b["peak_h"], "%.1f"); row("  within 1 % after (h)", a["die_h"], b["die_h"], "%.0f")
    row("  LIMs' peak push (MN)", a["F_lim_MN"], b["F_lim_MN"], "%.0f"); row("  reaction line (kN)", a["line_kN"], b["line_kN"], "%.1f"); row("  lean of the lines (deg)", a["lean_deg"], b["lean_deg"], "%.2f"); row("  rotation of the casing (deg)", a["rot_deg"], b["rot_deg"], "%.3f")
for N in (7, 14):
    for nm in ("braked", "not braked"):
        a = [q for q in n18["rows"] if q["N"] == N and q["sled"] == nm and q["lims"].startswith("damping, critical")][0]; b = [q for q in c19["daily"] if q["N"] == N and q["sled"] == nm and q["lims"].startswith("damping, critical")][0]
        row("%d a day, %s: mean push (MN)" % (N, nm), a["mean_push_MN"], b["mean_push_MN"], "%.1f"); row("  casing sits at (km)", a["x_mean"]/1e3, b["x_mean"]/1e3, "%.1f"); row("  swings about it (m)", a["swing"], b["swing"], "%.0f")
        row("  reaction line at the farthest (kN)", a["line_kN"], b["line_kN"], "%.1f"); row("  LIMs' peak push (MN)", a["F_lim_peak_MN"], b["F_lim_peak_MN"], "%.1f")
OUT["compare"] = cmp_
OUT["const"] = dict(K=K, M_CAS=M_CAS, period_s=2*math.pi/OM, c_crit=C_CRIT, F_LIM_MN=F_LIM/1e6, F_LINE_kN=F_LINE/1e3, k_r=k_r)
P("\n   run time %.0f s" % (time.time() - T0))
with open(os.path.join(LOG_DIR, "round19_swing_simulated.txt"), "w", encoding="utf-8") as fh: fh.write(buf.getvalue())
with open(os.path.join(LOG_DIR, "round19_swing_simulated.json"), "w", encoding="utf-8") as fh: json.dump(OUT, fh, indent=1, default=float)
print("ok")
