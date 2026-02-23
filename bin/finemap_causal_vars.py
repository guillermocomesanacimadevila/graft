#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

base = Path("../outputs/susie")
# pheno1_id
# pheno2_ic
# coords -> need a way in nf to pull this from process

def overlap(pheno1_id: str, pheno2_id: str, coords: str):
    pheno1_path = base / pheno1_id / f"locus_{coords}" / f"cs_{pheno1_id}.tsv"
    pheno2_path = base / pheno2_id / f"locus_{coords}" / f"cs_{pheno2_id}.tsv"
    pheno1_df = pd.read_csv(pheno1_path, sep="\t")
    pheno2_df = pd.read_csv(pheno2_path, sep="\t")
    merged = pheno1_df.merge(
        pheno2_df,
        how="inner",
        on="SNP",
        suffixes=(f"_{pheno1_id}", f"_{pheno2_id}")
    )
    out_dir = base / "overlap"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pheno1_id}_{pheno2_id}_{coords}.tsv"
    merged.to_csv(out_path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno1_id", required=True)
    ap.add_argument("--pheno2_id", required=True)
    ap.add_argument("--coords", required=True)
    args = ap.parse_args()
    overlap(args.pheno1_id, args.pheno2_id, args.coords)

if __name__ == "__main__":
    main()