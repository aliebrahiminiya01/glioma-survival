"""
Step 1 - Download TCGA LGG+GBM clinical/molecular data and build one clean table.

Source: cBioPortal study `lgggbm_tcga_pub` (Ceccarelli et al., Cell 2016),
1122 diffuse glioma patients (grade 2-4) with overall survival + IDH, 1p/19q,
MGMT and other molecular markers.

Output: data/processed/glioma.csv  - one row per patient
    time   : overall survival in months
    event  : 1 = died, 0 = alive at last follow-up (right-censored)
    + clinical / molecular features

Usage:
    python src/prepare_data.py
"""
import argparse
import json
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://www.cbioportal.org/api/studies/{study}/clinical-data?clinicalDataType={kind}&projection=SUMMARY"

# cBioPortal attribute -> our column name
PATIENT_COLS = {
    "OS_MONTHS": "time",
    "OS_STATUS": "os_status",
    "AGE": "age",
    "SEX": "sex",
    "KARNOFSKY_PERFORMANCE_SCORE": "kps",
    "HISTOLOGICAL_DIAGNOSIS": "histology",
}
SAMPLE_COLS = {
    "GRADE": "grade",
    "IDH_CODEL_SUBTYPE": "subtype",
    "IDH_STATUS": "idh",
    "MGMT_PROMOTER_STATUS": "mgmt",
    "CHR_7_GAIN_CHR_10_LOSS": "chr7_10",
    "ATRX_STATUS": "atrx",
    "TERT_PROMOTER_STATUS": "tert",
}


def download(study: str, kind: str, out: Path) -> None:
    """Save the raw JSON once; later runs reuse the local copy."""
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    url = API.format(study=study, kind=kind)
    print(f"downloading {kind} data ...")
    urllib.request.urlretrieve(url, out)


def to_wide(path: Path, cols: dict) -> pd.DataFrame:
    """cBioPortal returns long format (patient, attribute, value) -> pivot to one row per patient."""
    long = pd.DataFrame(json.load(open(path)))
    long = long[long["clinicalAttributeId"].isin(cols)]
    wide = long.pivot_table(index="patientId", columns="clinicalAttributeId",
                            values="value", aggfunc="first")
    return wide.rename(columns=cols)


def clean(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["kps"] = pd.to_numeric(df["kps"], errors="coerce")
    # "1:DECEASED" / "0:LIVING" -> 1 / 0
    df["event"] = df.pop("os_status").str.startswith("1").astype("float")

    # Survival models need both time and event. Patients without them are useless here.
    # time == 0 means no follow-up at all, so drop those too.
    df = df.dropna(subset=["time", "event"])
    df = df[df["time"] > 0]
    df["event"] = df["event"].astype(int)

    df["grade"] = df["grade"].str.replace("G", "").astype("Int64")
    df["chr7_10"] = df["chr7_10"].map({"Gain chr 7 & loss chr 10": "yes", "No combined CNA": "no"})
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", default="lgggbm_tcga_pub")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--out", default="data/processed/glioma.csv")
    args = parser.parse_args()

    raw = Path(args.raw_dir)
    download(args.study, "PATIENT", raw / "patient.json")
    download(args.study, "SAMPLE", raw / "sample.json")

    # every patient has exactly one sample in this study, so joining on patientId is safe
    df = to_wide(raw / "patient.json", PATIENT_COLS).join(to_wide(raw / "sample.json", SAMPLE_COLS))
    df = clean(df)
    df.columns.name = None

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index_label="patient_id")

    print(f"saved {len(df)} patients -> {out}")
    print(f"events (deaths): {df['event'].sum()}  |  censored: {(df['event'] == 0).sum()}")
    print("missing values per column:")
    print(df.isna().sum().to_string())


if __name__ == "__main__":
    main()
