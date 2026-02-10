# LSM Mass Driver Simulation

This directory contains a complete LSM (Linear Synchronous Motor) mass driver simulation for launching spacecraft from an orbital ring at 250 km altitude.

## Files Created

1. **lsm_config.py** - Configuration parameters and constants
2. **lsm_physics.py** - Physics calculations (fields, thrust, EMF, dynamics)
3. **lsm_simulation.py** - Main simulation loop, controller, and plotting

## Key Differences from LIM

| Feature | LIM | LSM |
|---------|-----|-----|
| Operation | Asynchronous (slip required) | Synchronous (no slip) |
| Sled | Passive γ-TiAl reaction plate | Superconducting DC field coils |
| Losses | Eddy current heating in plate | No eddy losses |
| Thermal limit | Plate temperature (1200 K max) | None (SC coils stay cold) |
| Control variable | Slip velocity | Load angle δ |
| Primary constraint | Thermal limit | Voltage limit (back-EMF) |
| Pole pitch | Multiple stages (50-200 m) | Single pitch (default 20 m) |
| Stator | Air-core HTS | Iron-core (default) or air-core |

## Quick Start

### Basic Run (with adjustable B_sled to reach 15 km/s)
```bash
python lsm_simulation.py --adjustable-B combined
```

### Fast Test Run
```bash
python lsm_simulation.py --quick --adjustable-B combined
```

### Fixed B_sled (voltage-limited to ~8 km/s)
```bash
python lsm_simulation.py combined
```

## Command-Line Options

### Configuration Options
- `--v_launch=N` - Target velocity in m/s (default: 15000)
- `--accel=N` - Maximum acceleration in g (default: 0.5)
- `--mass=N` - Spacecraft mass in tonnes (default: 500)
- `--sled-length=N` - Sled length in m (default: 5000)

### Stator Options
- `--iron-core` - Use iron-core stator (default)
- `--air-core` - Use air-core stator
- `--N_stator=N` - Stator turns per coil (default: 10)
- `--tau_p=N` - Pole pitch in m (default: 20)

### Sled Field Options
- `--B_sled=N` - Sled DC field in T (default: 0.10)
- `--adjustable-B` - Allow controller to reduce B_sled at high speed (REQUIRED to reach 15 km/s)

### Output Options
- `--quick` - Use 1-second timesteps (faster, less accurate)
- `--save-csv` - Export time-series data to CSV

### Graph Keywords
- `combined` - 9-panel combined plot
- `all` - All individual plots
- `velocity`, `accel`, `thrust`, `power`, `delta`, `frequency`, `voltage`, `current`, `bfield`, `gload` - Individual plots

## Example Results (Default Configuration with --adjustable-B)

```
Final velocity:       15,000.3 m/s
Target velocity:      15,000.0 m/s
Achievement:          100.0%
Total time:           51.0 minutes (0.85 hours)
Distance traveled:    22,943.0 km
Fraction of ring:     55.1%
Final KE:             129.380 TJ

Peak thrust:          5.64 MN
Peak power:           84.58 GW
Peak frequency:       375.0 Hz
Peak voltage:         100.0 kV
Peak current:         131.1 A
Peak load angle:      60.0°
Peak g-load:          2.585 g
```

## Physics Notes

### Thrust Equation
```
F = (B_stator × B_sled × sin(δ) / μ₀) × A_active
```

Where:
- `B_stator` = Stator AC field (T)
- `B_sled` = Sled DC field (T)
- `δ` = Load angle (electrical phase between fields)
- `A_active` = Total active area (m²)

### Iron-Core vs Air-Core Stator

**Iron-core** (default):
```
B_stator = min(B_sat, μ₀ × N × I / (2 × g))
```
- Much stronger field (up to 2.0 T saturation limit)
- Iron teeth concentrate the flux
- More efficient for LSM

**Air-core**:
```
B_stator = (2/π) × (μ₀ × N × I / w_coil) × arctan(w_coil / (2 × g))
```
- Weaker field (typically < 0.5 T)
- No saturation limit
- Simpler geometry

### Voltage Limit

The back-EMF induced in each stator coil is:
```
EMF = N × 2π × v × B_sled × w_coil
```

**Key insight**: This is independent of pole pitch τ_p!

At v = 15,000 m/s with N = 10, B_sled = 0.1 T, w_coil = 2.0 m:
```
EMF = 10 × 2π × 15000 × 0.1 × 2.0 = 188,500 V
```

This exceeds the 100 kV insulation limit! The controller must reduce B_sled to:
```
B_sled_max = V_limit / (N × 2π × v × w_coil)
            = 100,000 / (10 × 2π × 15000 × 2.0)
            = 0.053 T
```

### Why Adjustable B_sled is Critical

Without adjustable B_sled, the voltage limit is hit at:
```
v_max = V_limit / (N × 2π × B_sled × w_coil)
      = 100,000 / (10 × 2π × 0.10 × 2.0)
      ≈ 7,958 m/s
```

This is only 53% of the target 15 km/s!

With adjustable B_sled, the controller reduces the sled field at high velocities to stay within the voltage limit, allowing the full 15 km/s to be achieved.

## Controller Logic

The LSM controller adjusts two variables:

1. **Stator current `I`** (0 to I_peak = 2160 A)
   - Ramps up if thrust is limited (δ ≥ δ_max)
   - Ramps down if operating well below target load angle

2. **Load angle `δ`** (0 to δ_max = 60°)
   - Calculated to achieve desired thrust
   - Limited to δ_max for stability margin
   - Pull-out occurs at δ = 90°

3. **Sled field `B_sled`** (only if adjustable)
   - Reduced at high velocity to meet voltage limit
   - Starts at B_sled_nominal = 0.10 T
   - Reduces to ~0.053 T at 15 km/s

## Comparison with LIM

| Metric | LIM (Gamma-TiAl) | LSM (Iron-Core, Adjustable B) |
|--------|------------------|-------------------------------|
| Launch time | 103.3 minutes | 51.0 minutes |
| Distance | 53,605 km | 22,943 km |
| Peak thrust | 24.5 MN | 5.64 MN |
| Peak power | 177.5 GW | 84.58 GW |
| Peak frequency | ~1 Hz (S4) | 375 Hz |
| Peak voltage | 100 kV | 100 kV |
| Final plate temp | 1,164 K | N/A (no heating) |
| Eddy heating | 5.14 TJ | 0 (no eddy losses) |
| Primary constraint | Thermal limit | Voltage limit |
| Sled mass | Heavier (dense metal plate) | Lighter (SC coils + cryo) |

## Files and Directories

- `lsm_config.py` - Configuration module
- `lsm_physics.py` - Physics module
- `lsm_simulation.py` - Main simulation
- `graphs_lsm/` - Output directory for plots
- `output/lsm_data.csv` - Time-series data (if --save-csv used)

## Future Enhancements

Potential improvements:
1. Add inverter and electrical system modeling
2. Add cryogenic system heat loads (SC coils)
3. Multi-section sled (partial activation)
4. Optimized pole pitch sweep
5. Current waveform harmonics
6. Iron core saturation effects
7. Magnetic bearing / levitation forces
8. Start/stop transients
9. Emergency braking scenarios
10. Track layout optimization

## References

- "Orbital Ring Engineering" by Paul G de Jong
- LIM simulation: `lim_simulation.py`, `lim_config.py`, `lim_physics.py`
