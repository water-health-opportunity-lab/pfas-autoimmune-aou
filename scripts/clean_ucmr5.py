# cleanv5.py — UCMR5 PFAS occurrence data cleaning pipeline
#
# Input files (from raw_data/):
#   - ucmr5_all.csv: raw UCMR5 occurrence results (~1.86M sample rows)
#   - UCMR5_ZIPCodes.txt: zip codes served by each public water system
#   - UCMR5_AddtlDataElem.txt: self-reported survey data including prior PFAS
#     occurrence, known nearby PFAS sources, source types, and treatment info
#   - SDWA_PUB_WATER_SYSTEMS.csv: EPA SDWIS data with population served,
#     water source type, and service connections per public water system
#
# Processing steps:
#   1. Keep only key columns from raw data
#   2. Impute below-detection-limit values with MRL/sqrt(2)
#   3. Aggregate multiple samples into one row per PWSID + Contaminant
#   4. Join Zip3 codes and explode multi-zip PWSIDs into separate rows
#   5. Resolve numeric EPA region codes to state abbreviations via zip lookup
#   6. Drop unresolvable tribal PWSIDs, missing zips, incomplete contaminant
#      coverage, and US territories
#   7. Join self-reported survey data from UCMR5_AddtlDataElem.txt:
#      - PriorPFAS: whether PFAS was previously detected (Yes/DK/No)
#      - KnownPFASSources: whether nearby PFAS sources are known (Yes/DK/No)
#      - Binary flags for each reported source type (airport, military base,
#        manufacturing, farms, landfills, fire training, etc.)
#      - HasPFASTreatment: whether the system treats for PFAS
#   8. Join SDWIS population and water system data:
#      - PopulationServed: number of people served by the water system
#      - WaterSourceType: GW (groundwater) or SW (surface water)
#      - ServiceConnections: number of service connections
#   9. Save to output/ucmr5_analysis_zip3.csv
#
# Output: one row per PWSID + Contaminant + Zip3 combination

import os
import pandas as pd
import numpy as np

# Resolve paths relative to this script's location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(SCRIPT_DIR, "..", "raw_data")
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "output")

# ── Load raw occurrence data ──
df = pd.read_csv(
    os.path.join(RAW_DIR, "ucmr5_all.csv"),
    dtype={"PWSID": str, "State": str, "FacilityID": str},
    low_memory=False,
)

# Keep only the important columns
df = df[[
    "PWSID",
    "CollectionDate",
    "Contaminant",
    "MRL",
    "Units",
    "AnalyticalResultsSign",
    "AnalyticalResultValue",
    "State",
]]

# Impute below-MRL samples with MRL / sqrt(2)
below_mrl = df["AnalyticalResultsSign"] == "<"
df.loc[below_mrl, "AnalyticalResultValue"] = df.loc[below_mrl, "MRL"] / np.sqrt(2)

# Create detection flag: 1 = detected (exceeds MRL), 0 = below MRL
df["Detected"] = (~below_mrl).astype(int)

# Aggregate to one row per PWSID + Contaminant
agg = df.groupby(["PWSID", "Contaminant"]).agg(
    MeanConcentration=("AnalyticalResultValue", "mean"),
    MaxConcentration=("AnalyticalResultValue", "max"),
    MRL=("MRL", "first"),
    Units=("Units", "first"),
    DetectionRate=("Detected", "mean"),
    NumSamples=("Detected", "count"),
    NumDetections=("Detected", "sum"),
    State=("State", "first"),
).reset_index()

# ── Join Zip3 codes ──
zips = pd.read_csv(os.path.join(RAW_DIR, "UCMR5_ZIPCodes.txt"), sep="\t", dtype=str)
zips["ZIP3"] = zips["ZIPCODE"].str[:3]
zips_grouped = zips.groupby("PWSID")["ZIP3"].apply(
    lambda x: ",".join(sorted(x.unique()))
).reset_index()
zips_grouped.rename(columns={"ZIP3": "Zip3"}, inplace=True)

agg = agg.merge(zips_grouped, on="PWSID", how="left")

# ── Resolve numeric state codes ──
is_numeric_state = agg["State"].str.match(r"^\d+$")
has_alpha_and_zip = (~is_numeric_state) & agg["Zip3"].notna()

prefix_records = []
for _, row in agg.loc[has_alpha_and_zip, ["State", "Zip3"]].drop_duplicates().iterrows():
    for z in row["Zip3"].split(","):
        prefix_records.append({"prefix": z, "State": row["State"]})
prefix_map = (
    pd.DataFrame(prefix_records)
    .groupby("prefix")["State"]
    .agg(lambda x: x.mode()[0])
    .to_dict()
)

needs_fix = is_numeric_state & agg["Zip3"].notna()
agg.loc[needs_fix, "State"] = agg.loc[needs_fix, "Zip3"].apply(
    lambda z: prefix_map.get(z.split(",")[0], None)
)

# ── Drop rows ──
agg = agg[~agg["State"].str.match(r"^\d+$", na=False)]
agg = agg[agg["Zip3"].notna()]
contam_counts = agg.groupby("PWSID")["Contaminant"].transform("count")
agg = agg[contam_counts == 29]
territories = {"AS", "GU", "MP", "VI", "PR", "NN"}
agg = agg[~agg["State"].isin(territories)]

# Replace raw Units with symbol
agg["Units"] = "µg/L"

# Explode comma-separated Zip3 into separate rows
agg["Zip3"] = agg["Zip3"].str.split(",")
agg = agg.explode("Zip3").reset_index(drop=True)

# ── Join self-reported survey data from UCMR5_AddtlDataElem.txt ──
addtl = pd.read_csv(
    os.path.join(RAW_DIR, "UCMR5_AddtlDataElem.txt"), sep="\t", dtype=str,
)

# Helper: for each PWSID, pick the most informative response across sample events
# Priority: Yes > DK > No (if they ever said Yes, that matters most)
RESPONSE_PRIORITY = {"Yes": 2, "DK": 1, "No": 0}

def best_response(series):
    return max(series, key=lambda x: RESPONSE_PRIORITY.get(x, -1))

# PriorPFAS: has the system previously detected PFAS?
pfas_occ = addtl[addtl["AdditionalDataElement"] == "PFASOccurrence"]
pfas_occ_agg = pfas_occ.groupby("PWSID")["Response"].agg(best_response).reset_index()
pfas_occ_agg.rename(columns={"Response": "PriorPFAS"}, inplace=True)

# KnownPFASSources: does the system know of nearby PFAS sources?
pfas_src = addtl[addtl["AdditionalDataElement"] == "PotentialPFASSources"]
pfas_src_agg = pfas_src.groupby("PWSID")["Response"].agg(best_response).reset_index()
pfas_src_agg.rename(columns={"Response": "KnownPFASSources"}, inplace=True)

# Source type detail: create binary columns for each reported source type
SOURCE_CODES = {
    "AO": "Src_Airport",
    "MB": "Src_MilitaryBase",
    "MF": "Src_Manufacturing",
    "FT": "Src_FireTraining",
    "CW": "Src_ChemicalWaste",
    "WM": "Src_WasteMgmt",
    "FM": "Src_Farm",
    "FP": "Src_FireStation",
    "CC": "Src_ChemicalCompany",
    "MM": "Src_Mining",
    "PC": "Src_PlatingChrome",
    "SC": "Src_Semiconductor",
    "TA": "Src_TextileApparel",
    "HW": "Src_HazardousWaste",
    "PR": "Src_PaperRecycling",
    "OG": "Src_OilGas",
    "PP": "Src_PaperPulp",
    "CE": "Src_ChemicalElectronics",
    "FF": "Src_FoodProcessing",
    "UW": "Src_UpstreamWWTP",
    "CT": "Src_CoatingTanning",
    "PS": "Src_PFASSpill",
    "UT": "Src_UpstreamDischarge",
    "OT": "Src_Other",
}

detail = addtl[addtl["AdditionalDataElement"] == "PotentialPFASSourcesDetail"]
# Get unique source codes per PWSID
detail_pivot = detail.groupby("PWSID")["Response"].apply(set).reset_index()

# Create binary columns
for code, col_name in SOURCE_CODES.items():
    detail_pivot[col_name] = detail_pivot["Response"].apply(lambda s: int(code in s))
detail_pivot = detail_pivot.drop(columns=["Response"])

# HasPFASTreatment: is the system treating for PFAS? (anything other than NMT)
pfas_treat = addtl[addtl["AdditionalDataElement"] == "PFASTreatment"]
treat_agg = pfas_treat.groupby("PWSID")["Response"].apply(
    lambda x: int(any(r != "NMT" for r in x))
).reset_index()
treat_agg.rename(columns={"Response": "HasPFASTreatment"}, inplace=True)

# Merge all survey data onto main dataset
agg = agg.merge(pfas_occ_agg, on="PWSID", how="left")
agg = agg.merge(pfas_src_agg, on="PWSID", how="left")
agg = agg.merge(detail_pivot, on="PWSID", how="left")
agg = agg.merge(treat_agg, on="PWSID", how="left")

# Fill NaN for source detail columns (PWSIDs not in survey get 0)
src_cols = list(SOURCE_CODES.values())
agg[src_cols] = agg[src_cols].fillna(0).astype(int)
agg["HasPFASTreatment"] = agg["HasPFASTreatment"].fillna(0).astype(int)

# ── Join SDWIS population and water system data ──
sdwa = pd.read_csv(
    os.path.join(RAW_DIR, "SDWA_PUB_WATER_SYSTEMS.csv"),
    dtype=str,
    usecols=["PWSID", "POPULATION_SERVED_COUNT", "GW_SW_CODE", "SERVICE_CONNECTIONS_COUNT"],
)
sdwa.rename(columns={
    "POPULATION_SERVED_COUNT": "PopulationServed",
    "GW_SW_CODE": "WaterSourceType",
    "SERVICE_CONNECTIONS_COUNT": "ServiceConnections",
}, inplace=True)
sdwa["PopulationServed"] = pd.to_numeric(sdwa["PopulationServed"], errors="coerce")
sdwa["ServiceConnections"] = pd.to_numeric(sdwa["ServiceConnections"], errors="coerce")

agg = agg.merge(sdwa, on="PWSID", how="left")

# ── Save ──
agg.to_csv(os.path.join(OUT_DIR, "ucmr5_analysis_zip3.csv"), index=False)

print("Output shape:", agg.shape)
print("\nColumns:")
print(agg.dtypes.to_string())
print("\nFirst 5 rows (core columns):")
print(agg[["PWSID", "Contaminant", "MeanConcentration", "State", "Zip3",
           "PopulationServed", "WaterSourceType", "ServiceConnections"]].head(5).to_string())
print("\nMissing values:")
missing = agg.isnull().sum()
missing = missing[missing > 0]
if len(missing) > 0:
    print(missing.to_string())
else:
    print("  None")
print(f"\nSDWIS coverage:")
print(f"  PopulationServed filled: {agg['PopulationServed'].notna().sum():,} / {len(agg):,}")
print(f"  WaterSourceType filled: {agg['WaterSourceType'].notna().sum():,} / {len(agg):,}")
print(f"\nPopulationServed stats:")
print(agg["PopulationServed"].describe().to_string())
print(f"\nWaterSourceType breakdown:")
print(agg["WaterSourceType"].value_counts().to_string())
