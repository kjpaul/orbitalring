# Input files the scripts read

Everything in this folder is tracked by git. Nothing here is written by a normal run; the three solver caches are rewritten only by the recompute commands named below.

## HTS tape data

`THEVA_Pro-Line_Advanced_Pinning_2G_HTS_70K_Field_Dependence.csv`

Source: THEVA Dünnschichttechnik GmbH, measured data for the product Pro-Line Advanced Pinning 2G HTS tape (coated conductor), critical current against applied magnetic field at 70 K. Used with permission of THEVA.

Columns, one row per measurement point:

| Column | Meaning |
|---|---|
| Temperature (K) | sample temperature, about 70 K |
| Applied field (T) | magnitude of the applied magnetic field |
| Applied field angle (°) | angle between the field and the tape, 0, 30, 90 and 120 degrees (0 is the field perpendicular to the tape face, the convention Chapter 6 uses) |
| Critical current per cm width (A/cm) | critical current of the tape per centimetre of tape width |

Read by `bearing_design.py` (`read_theva`), `bearing_mass.py`, `final_bearing.py` and `compensation_tracks.py`. The scripts multiply the per-centimetre value by the tape width (3 mm or 12 mm) and take the lower of the 0 and 90 degree curves at the field in the winding. `bearing_optimize.py` reads the angle too.

`THEVA_Pro-Line_Advanced_Pinning_2G_HTS_77.5K_Field_Dependence.csv`: the same measurement at 77.5 K, same columns and source, used with permission of THEVA. Read by `bearing_optimize.py` for the warming a tape can take before its limit falls to the operating current (Part 2 of the coil optimization).

`THEVA_ProLine_TPL4420c_Robinson_2019_30K_50K_60K_70K.csv`: a different sample, added 2026-11-17 for the temperature dependence only. It is the Robinson Research Institute's published measurement of a THEVA Pro-Line 2G HTS wire of 2019 (4 mm wide, designation TPL4420c, not the Advanced Pinning grade): Wimbush, S. and Strickland, N. (2019), "Critical current characterisation of THEVA Pro-Line 2G HTS superconducting wire", figshare, dataset, https://doi.org/10.6084/m9.figshare.3759327.v2, licence CC BY 4.0. The file holds the rows at 30, 50, 60 and 70 K of the dataset's `nnK angle dependences.xlsx` files (set temperature, set field 0 to 7 T, set angle 0 to 240 degrees, Ic per cm of width). Read by `tools/round11_stator_model.py`, which uses only the ratio Ic(T) / Ic(70 K) at each field and angle and applies it to the Advanced Pinning tape's 70 K curve. The full dataset has every 2.5 K from 20 to 90 K.

## The mass driver's stator

`lsm_stator_ruled.json`: the ruled stator's tape current, coil inductance and loss table (hysteresis, soldered layers, and the mean square field for the eddy loss, per meter of track of one side and per cycle, at 24 speeds). Written by `python tools/round13_stator_table.py` from the tape-by-tape model of `tools/round12_stator_model.py`; read by `lsm_config.py` when the stator is the default (`ORBITAL_RING_LSM_STATOR=ruled`).

## Solver caches

Results of the two-dimensional magnetostatic solver (`magnetostatic_2d.py`) that several scripts read instead of solving again. Each one can be recomputed from scratch; the recompute overwrites the file here.

| File | Written by | Read by | Recompute | Time |
|---|---|---|---|---|
| `compensation_fields.json` | `compensation_fields.py` | `compensation_tracks.py`, `bearing_design.py`, `bearing_mass.py`, `final_bearing.py` | `python compensation_fields.py` (every run of that file recomputes it) | about 4 minutes |
| `bearing_sweep.json` | `bearing_design.py` Part B | `bearing_design.py`, `bearing_mass.py`, `final_bearing.py` | `python bearing_design.py --resweep` | about 42 minutes |
| `bearing_air.json` | `bearing_mass.py` Part B | `bearing_mass.py`, `final_bearing.py` | `python bearing_mass.py --resolve` | about 24 minutes |

## Design-history logs

Logs of earlier runs that two scripts quote. They are kept so that the scripts print what the book's logs print; a current run does not regenerate them.

| File | What it is | Read by |
|---|---|---|
| `levitation_8620_2026-10-04.txt` | `levitation_simulation.py` log of 2026-10-04 with the cable speed relative to the casing hard coded at 8,620 m/s (the 250 km value) | `compensation_tracks.py` Part E, the comparison of the two logs |
| `levitation_2026-10-11.txt` | the same script after the fix (speed from `ring_altitude.py`), run on 2026-10-11 with the cable hardware at 4,076 kg/m; the output folder in its five path lines was shortened | `compensation_tracks.py` Part E |
| `compensation_tracks_2026-10-11.txt` | `compensation_tracks.py` log of 2026-10-11 (cable hardware 4,076 kg/m) | `bearing_design.py` Part F reads the Part A pull at each speed from it |
