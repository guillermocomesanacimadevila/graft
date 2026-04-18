#!/usr/bin/env python3
import argparse
import pandas as pd 
from pathlib import Path

def align_alleles(pheno1_id, pheno2_id, sumstats1, sumstats2, out_dir):
    p1 = Path(sumstats1)
    p2 = Path(sumstats2)
    out_base = Path(out_dir)
    out1 = out_base / "QC" / pheno1_id / f"{pheno1_id}.ldsc_ready_neff.tsv"
    out2 = out_base / "QC" / pheno2_id / f"{pheno2_id}.ldsc_ready_neff.tsv"
    df1 = pd.read_csv(p1, sep="\t", dtype=str, low_memory=False)
    df2 = pd.read_csv(p2, sep="\t", dtype=str, low_memory=False)
    print(f"{pheno1_id}: {len(df1):,}")
    print(f"{pheno2_id}: {len(df2):,}")
    df1["SNP"] = df1["SNP"].str.strip()
    df2["SNP"] = df2["SNP"].str.strip()
    df1["A1"] = df1["A1"].str.strip().str.upper()
    df1["A2"] = df1["A2"].str.strip().str.upper()
    df2["A1"] = df2["A1"].str.strip().str.upper()
    df2["A2"] = df2["A2"].str.strip().str.upper()
    df1["BETA"] = pd.to_numeric(df1["BETA"], errors="coerce")
    df2["BETA"] = pd.to_numeric(df2["BETA"], errors="coerce")
    df1["FRQ"] = pd.to_numeric(df1["FRQ"], errors="coerce")
    df2["FRQ"] = pd.to_numeric(df2["FRQ"], errors="coerce")
    df1 = df1.drop_duplicates("SNP").dropna(subset=["SNP","A1","A2","BETA","FRQ"])
    df2 = df2.drop_duplicates("SNP").dropna(subset=["SNP","A1","A2","BETA","FRQ"])
    ref = df1[["SNP","A1","A2"]].rename(columns={"A1":"A1_ref","A2":"A2_ref"})
    m = ref.merge(df2, on="SNP", how="inner")
    direct = (m["A1"] == m["A1_ref"]) & (m["A2"] == m["A2_ref"])
    swapped = (m["A1"] == m["A2_ref"]) & (m["A2"] == m["A1_ref"])
    m = m[direct | swapped]
    flip = (m["A1"] == m["A2_ref"]) & (m["A2"] == m["A1_ref"])
    m.loc[flip, ["A1","A2"]] = m.loc[flip, ["A2","A1"]].values
    m.loc[flip, "BETA"] = -m.loc[flip, "BETA"]
    m.loc[flip, "FRQ"] = 1 - m.loc[flip, "FRQ"]
    keep = m["SNP"].drop_duplicates()
    order = pd.DataFrame({"SNP": keep})
    df1 = order.merge(df1[df1["SNP"].isin(keep)], on="SNP", how="left")
    df2 = order.merge(m[df2.columns], on="SNP", how="left")
    print(f"aligned SNPs: {len(df1):,}")
    out1.parent.mkdir(parents=True, exist_ok=True)
    out2.parent.mkdir(parents=True, exist_ok=True)
    df1.to_csv(out1, sep="\t", index=False)
    df2.to_csv(out2, sep="\t", index=False)
    print(f"overwritten: {out1}")
    print(f"overwritten: {out2}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno1_id", required=True)
    p.add_argument("--pheno2_id", required=True)
    p.add_argument("--pheno1_gwas", required=True)
    p.add_argument("--pheno2_gwas", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    align_alleles(
        pheno1_id=args.pheno1_id,
        pheno2_id=args.pheno2_id,
        sumstats1=args.pheno1_gwas,
        sumstats2=args.pheno2_gwas,
        out_dir=args.out_dir,
    )

if __name__ == "__main__":
    main()