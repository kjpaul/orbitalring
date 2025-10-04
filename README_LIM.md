# Orbital Ring LIM Physics Simulator

Run:
```bash
python orbital_ring_lim.py --config example_lim_config.yaml --plot velocities,tension,weight --verbose --explain
```

- Per-site physics from coil geometry (turns, gap, pitch), slip, and power/voltage caps
- Controller starts cold and ramps when limits allow
- Plots use 30-day months, downsampled; skip initial warmup (plot_skip_s)

To change physics, edit `example_lim_config.yaml` under `lim:` (no hard-coded constants).
