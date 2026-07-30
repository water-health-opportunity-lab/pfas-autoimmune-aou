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
#   10. Aggregate to one row per Zip3 + Contaminant, using population-served
#       as weights for concentration and detection rate metrics
#
# Output: one row per Zip3 + Contaminant combination

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
df = df[
    [
        "PWSID",
        "CollectionDate",
        "Contaminant",
        "MRL",
        "Units",
        "AnalyticalResultsSign",
        "AnalyticalResultValue",
        "State",
    ]
]

# Impute below-MRL samples with MRL / sqrt(2)
below_mrl = df["AnalyticalResultsSign"] == "<"
df.loc[below_mrl, "AnalyticalResultValue"] = df.loc[below_mrl, "MRL"] / np.sqrt(2)

# Create detection flag: 1 = detected (exceeds MRL), 0 = below MRL
df["Detected"] = (~below_mrl).astype(int)

# Aggregate to one row per PWSID + Contaminant
agg = (
    df.groupby(["PWSID", "Contaminant"])
    .agg(
        MeanConcentration=("AnalyticalResultValue", "mean"),
        MaxConcentration=("AnalyticalResultValue", "max"),
        MRL=("MRL", "first"),
        Units=("Units", "first"),
        DetectionRate=("Detected", "mean"),
        NumSamples=("Detected", "count"),
        NumDetections=("Detected", "sum"),
        State=("State", "first"),
    )
    .reset_index()
)

# ── Join Zip3 codes ──
zips = pd.read_csv(os.path.join(RAW_DIR, "UCMR5_ZIPCodes.txt"), sep="\t", dtype=str)
zips["ZIP3"] = zips["ZIPCODE"].str[:3]
zips_unique = zips[["PWSID", "ZIP3"]].drop_duplicates().rename(columns={"ZIP3": "Zip3"})

agg = agg.merge(zips_unique, on="PWSID", how="left")

# Drop rows without Zip3 (PWSIDs not in the zip code file)
# drops 5192 out of 355067
agg = agg[agg["Zip3"].notna()]


# ── Resolve numeric state codes ──
is_numeric_state = agg["State"].str.match(r"^\d+$")

prefix_map = (
    agg.loc[~is_numeric_state, ["Zip3", "State"]]
    .drop_duplicates()
    .groupby("Zip3")["State"]
    .agg(lambda x: x.mode()[0])
    .to_dict()
)

agg.loc[is_numeric_state, "State"] = agg.loc[is_numeric_state, "Zip3"].map(prefix_map)

# ── Drop rows ──

contam_counts = agg.groupby(["PWSID", "Zip3"])["Contaminant"].transform("count")


# drops 95310 out of 349872 rows(a lot)
# agg = agg[contam_counts == 29]

print("\n── Detections above MRL per contaminant ──")
detections_per_contam = (
    agg.groupby("Contaminant")["NumDetections"].sum().sort_values(ascending=False)
)
print(detections_per_contam.to_string())

# remove the contaminants that are above the mrl <1000 times, 8 contaminants are left out of 29 (21 removed)
low_detection_contaminants = {
    "PFDA",
    "PFUnA",
    "8:2 FTS",
    "PFHpS",
    "NFDHA",
    "ADONA",
    "PFDoA",
    "4:2 FTS",
    "PFMPA",
    "NMeFOSAA",
    "PFMBA",
    "NEtFOSAA",
    "9Cl-PF3ONS",
    "PFTrDA",
    "11Cl-PF3OUdS",
    "PFTA",
    "PFEESA",
    "PFPeS",
    "HFPO-DA",
    "PFNA",
    "6:2 FTS",
}
agg = agg[~agg["Contaminant"].isin(low_detection_contaminants)]

territories = {"AS", "GU", "MP", "VI", "PR", "NN"}
agg = agg[~agg["State"].isin(territories)]

# ── Join self-reported survey data from UCMR5_AddtlDataElem.txt ──
addtl = pd.read_csv(
    os.path.join(RAW_DIR, "UCMR5_AddtlDataElem.txt"),
    sep="\t",
    dtype=str,
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
treat_agg = (
    pfas_treat.groupby("PWSID")["Response"]
    .apply(lambda x: int(any(r != "NMT" for r in x)))
    .reset_index()
)
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


agg.to_csv(os.path.join(OUT_DIR, "ucmr5_pwsid_contaminant.csv"), index=False)
# start of second half, there are 251053 rows of pwsid & contaminant combinations


# ── Join SDWIS population and water system data ──
sdwa = pd.read_csv(
    os.path.join(RAW_DIR, "SDWA_PUB_WATER_SYSTEMS.csv"),
    dtype=str,
    usecols=[
        "PWSID",
        "POPULATION_SERVED_COUNT",
        "GW_SW_CODE",
        "SERVICE_CONNECTIONS_COUNT",
    ],
)
sdwa.rename(
    columns={
        "POPULATION_SERVED_COUNT": "PopulationServed",
        "GW_SW_CODE": "WaterSourceType",
        "SERVICE_CONNECTIONS_COUNT": "ServiceConnections",
    },
    inplace=True,
)
sdwa["PopulationServed"] = pd.to_numeric(sdwa["PopulationServed"], errors="coerce")
sdwa["ServiceConnections"] = pd.to_numeric(sdwa["ServiceConnections"], errors="coerce")

agg = agg.merge(sdwa, on="PWSID", how="left")

# Drop rows without population data (can't weight them)
# drops 232 rows out of 251053
agg = agg[agg["PopulationServed"].notna() & (agg["PopulationServed"] > 0)]


# ── Aggregate to one row per Zip3 + Contaminant (population-weighted) ──
def pop_weighted_mean(group, col):
    weights = group["PopulationServed"]
    return np.average(group[col], weights=weights)


zip3_agg = (
    agg.groupby(["Zip3", "Contaminant"])
    .apply(
        lambda g: pd.Series(
            {
                "MeanConcentration": pop_weighted_mean(g, "MeanConcentration"),
                "MaxConcentration": pop_weighted_mean(g, "MaxConcentration"),
                "MRL": g["MRL"].iloc[0],
                "Units": g["Units"].iloc[0],
                "DetectionRate": pop_weighted_mean(g, "DetectionRate"),
                "NumSamples": g["NumSamples"].sum(),
                "NumDetections": g["NumDetections"].sum(),
                "NumSystems": g["PWSID"].nunique(),
                "State": g["State"].mode().iloc[0],
                "PopulationServed": g["PopulationServed"].sum(),
                "ServiceConnections": g["ServiceConnections"].sum(),
                "WaterSourceType": g["WaterSourceType"].mode().iloc[0]
                if g["WaterSourceType"].notna().any()
                else np.nan,
                "PriorPFAS": "Yes"
                if (g["PriorPFAS"] == "Yes").any()
                else (
                    "DK"
                    if (g["PriorPFAS"] == "DK").any()
                    else ("No" if g["PriorPFAS"].notna().any() else np.nan)
                ),
                "KnownPFASSources": "Yes"
                if (g["KnownPFASSources"] == "Yes").any()
                else (
                    "DK"
                    if (g["KnownPFASSources"] == "DK").any()
                    else ("No" if g["KnownPFASSources"].notna().any() else np.nan)
                ),
                "HasPFASTreatment": int((g["HasPFASTreatment"] == 1).any()),
                **{col: pop_weighted_mean(g, col) for col in src_cols},
            }
        ),
        include_groups=False,
    )
    .reset_index()
)
# ── Save ──
zip3_agg.to_csv(os.path.join(OUT_DIR, "ucmr5_zip3_contaminant.csv"), index=False)

print("Output shape:", zip3_agg.shape)
print("\nColumns:")
print(zip3_agg.dtypes.to_string())
print("\nFirst 5 rows (core columns):")
print(
    zip3_agg[
        [
            "Zip3",
            "Contaminant",
            "MeanConcentration",
            "MaxConcentration",
            "State",
            "NumSystems",
            "PopulationServed",
        ]
    ]
    .head(5)
    .to_string()
)
print("\nMissing values:")
missing = zip3_agg.isnull().sum()
missing = missing[missing > 0]
if len(missing) > 0:
    print(missing.to_string())
else:
    print("  None")
print(f"\nUnique Zip3 values: {zip3_agg['Zip3'].nunique()}")
print(f"Unique contaminants: {zip3_agg['Contaminant'].nunique()}")
print(f"\nPopulationServed stats:")
print(zip3_agg["PopulationServed"].describe().to_string())
