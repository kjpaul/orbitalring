# Radiation calculations for Chapter 5 (2026-10-04)

Run everything from this folder. Needs: `pip install aep8 ppigrf astropy numpy scipy netCDF4 geant4-pybind`.
The first import of geant4_pybind downloads about 2 GB of Geant4 data.

| Script | What it does | Output |
|---|---|---|
| `ap8_spectra.py` | AP-8 integral proton spectra on the equator in AP-8's own 1960s field | `ap8_spectra.npz`, `ap8_spectra.txt` |
| `igrf_sector.py` | McIlwain L and B/B0 by field-line tracing in IGRF, checked against the AP-8 package field at 1960; field direction | `igrf_sector.txt`, `igrf_LB.npz` |
| `igrf_altscan.py` | AP-8 MIN and MAX evaluated with the IGRF 2026 field, 500 to 1,000 km | `igrf_altscan_clean.txt`, `ap8_igrf2026_spectra.npz` |
| `poes_equator.py` | MetOp-C MEPED proton flux on the equator at 823 km, August and September 2026 (files from NOAA NCEI, see `poes/dl.sh`) | `poes_equator.txt` |
| `igrf_823.py` | AP-8 with the IGRF 2026 field at 823 km next to the MetOp-C measurement | `igrf_823.txt` |
| `g4_depthdose.py` | Geant4 depth-dose runs in a water slab (protons 12 MeV to 1.1 GeV with QGSP_BIC_HP; cosmic ray H, He, O, Si, Fe with FTFP_BERT) | `out_trap/`, `out_trap_lowE/`, `out_gcr/` |
| `response.py` | PSTAR tables, quality factor, loader for the Geant4 runs | |
| `fold_trapped.py`, `tables_trapped.py` | trapped-proton dose behind a water sphere | `tables_trapped_igrf2026.txt`, `tables_trapped_epoch.txt` |
| `pitch.py` | pitch-angle distribution from AP-8 | |
| `sweetspot.py`, `sweetspot_trapped.py` | dose under the cable by direction, trapped protons | `sweetspot_trapped.txt` |
| `gcr.py`, `sweetspot_run.py`, `sweetspot_scan.py` | cosmic ray dose behind shielding and under the cable, with the ISS and atmosphere checks | `gcr_per_sr.txt`, `sweetspot_gcr.txt` |
| `solar_cal.py`, `ch5_tables.py` | solar cell thresholds, electronics dose, panel life | `solar_cal.txt`, `ch5_tables.txt` |
| `pstar_*.txt`, `estar_al.txt` | NIST PSTAR and ESTAR tables as downloaded | |

Geant4 commands used: `g4_depthdose.py QGSP_BIC_HP proton out_trap <events:MeV ...>` (see `jobs_p.txt`), `lowE.sh`, `q1.sh`, `q2.sh`.
