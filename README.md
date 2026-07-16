# UCMR5 PFAS Occurrence Data

Cleaning pipeline for EPA's [Fifth Unregulated Contaminant Monitoring Rule (UCMR5)](https://www.epa.gov/dwucmr/fifth-unregulated-contaminant-monitoring-rule) occurrence data, focused on per- and polyfluoroalkyl substances (PFAS) detected in U.S. public water systems.

## Data Sources

The raw data files are too large for GitHub. Download them and place them in `raw_data/`:

| File | Source | Description |
|------|--------|-------------|
| `ucmr5_all.csv` | [EPA UCMR5 Results](https://www.epa.gov/dwucmr/occurrence-data-unregulated-contaminant-monitoring-rule) | ~1.86M sample-level occurrence results |
| `UCMR5_ZIPCodes.txt` | EPA UCMR5 Results | ZIP codes served by each public water system |
| `UCMR5_AddtlDataElem.txt` | EPA UCMR5 Results | Self-reported survey data (prior PFAS occurrence, nearby sources, treatment) |
| `SDWA_PUB_WATER_SYSTEMS.csv` | [EPA SDWIS](https://www.epa.gov/enviro/sdwis-model) | Population served, water source type, and service connections |

## Pipeline

`scripts/clean_ucmr5.py` processes the raw data into a single analysis-ready CSV:

1. Keeps key columns and imputes below-detection-limit values with MRL/sqrt(2)
2. Aggregates samples to one row per water system (PWSID) + contaminant
3. Joins 3-digit ZIP codes and explodes multi-ZIP systems into separate rows
4. Resolves numeric EPA region codes to state abbreviations
5. Drops unresolvable tribal systems, missing ZIPs, incomplete contaminant coverage, and U.S. territories
6. Joins self-reported survey data (prior PFAS detections, known nearby sources, treatment status)
7. Joins SDWIS data (population served, water source type, service connections)
8. Aggregates to one row per ZIP3 + contaminant, using population served as weights for concentration and detection rate metrics

**Output:** `output/ucmr5_analysis_zip3.csv` — one row per ZIP3 + contaminant combination, with population-weighted statistics.

## Setup

```bash
python -m venv my_env
source my_env/bin/activate
pip install -r requirements.txt
```

## Usage

```bash
python scripts/clean_ucmr5.py
python scripts/inspect_output.py   # preview the output
```
