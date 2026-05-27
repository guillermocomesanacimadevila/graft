#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import os
import argparse

def read_coord_table(path):
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    pos_col = "POS" if "POS" in header else "BP"
    df = pd.read_csv(path, sep="\t", usecols=["SNP", "CHR", pos_col])
    if pos_col == "POS":
        df.rename(columns={"POS": "BP"}, inplace=True)
    return df

def read_full_table(path):
    df = pd.read_csv(path, sep="\t")
    if "POS" in df.columns:
        df.rename(columns={"POS": "BP"}, inplace=True)
    return df

def merge_overlapping_loci(lead_snps_df, pheno1_df, window):
    rows = []
    for _, r in lead_snps_df.iterrows():
        snp = r["lead_snp"]
        lead1 = pheno1_df.loc[pheno1_df["SNP"] == snp]
        if lead1.empty:
            continue
        lead1 = lead1.iloc[0]
        chr_ = int(lead1["CHR"])
        bp = int(lead1["BP"])
        start = max(0, bp - window)
        end = bp + window
        row = r.to_dict()
        row["CHR"] = chr_
        row["BP"] = bp
        row["START"] = start
        row["END"] = end
        row["lead_snps_merged"] = str(snp)
        rows.append(row)
    if len(rows) == 0:
        return pd.DataFrame()
    loci = pd.DataFrame(rows).sort_values(["CHR", "START", "END"]).reset_index(drop=True)
    merged = []
    for _, row in loci.iterrows():
        row = row.to_dict()
        if not merged:
            merged.append(row)
            continue
        last = merged[-1]
        if row["CHR"] == last["CHR"] and row["START"] <= last["END"]:
            last["START"] = min(last["START"], row["START"])
            last["END"] = max(last["END"], row["END"])
            last["lead_snps_merged"] = str(last["lead_snps_merged"]) + "," + str(row["lead_snp"])
        else:
            merged.append(row)
    merged_df = pd.DataFrame(merged)
    merged_df["locus"] = merged_df.apply(
        lambda x: f"locus_chr{int(x['CHR'])}_{int(x['START'])}_{int(x['END'])}",
        axis=1
    )
    return merged_df

def define_loci(pheno1: str,
                pheno2: str,
                pheno1_prefix: str,
                pheno2_prefix: str,
                clump_path: str,
                out_dir: str,
                window: int):

    pheno1_df = read_coord_table(pheno1)
    pheno2_df = read_coord_table(pheno2)
    pheno1_full = read_full_table(pheno1)
    pheno2_full = read_full_table(pheno2)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    clump_path = Path(clump_path)
    rows = []
    for dir in os.listdir(clump_path):
        full_path = clump_path / dir
        if os.path.isdir(full_path) and dir.startswith("locus_"):
            locus_coords = "_".join(dir.split("_")[-3:])
            lead_file = full_path / "lead_snps.tsv"
            lead_df = pd.read_csv(lead_file, sep="\t")
            lead_snp = lead_df["lead_snp"].iloc[0]
            rows.append({
                "trait1": pheno1_prefix,
                "trait2": pheno2_prefix,
                "pair": f"{pheno1_prefix}_{pheno2_prefix}",
                "locus": dir,
                "coords": locus_coords,
                "lead_snp": lead_snp
            })
    lead_snps_df = pd.DataFrame(rows, columns=["trait1", "trait2", "pair", "locus", "coords", "lead_snp"])
    lead_snps_df = merge_overlapping_loci(lead_snps_df, pheno1_df, window)
    if lead_snps_df.empty:
        lead_snps_df.to_csv(out / "lead_snps.tsv", sep="\t", index=False)
        return lead_snps_df
    lead_snps_df = lead_snps_df.sort_values(["CHR", "START", "END"]).reset_index(drop=True)
    lead_snps_df.to_csv(out / "lead_snps.tsv", sep="\t", index=False)
    for _, r in lead_snps_df.iterrows():
        chr_ = int(r["CHR"])
        start = int(r["START"])
        end = int(r["END"])
        pheno1_loc = pheno1_df.loc[(pheno1_df["CHR"] == chr_) & (pheno1_df["BP"] >= start) & (pheno1_df["BP"] <= end)].copy()
        pheno2_loc = pheno2_df.loc[(pheno2_df["CHR"] == chr_) & (pheno2_df["BP"] >= start) & (pheno2_df["BP"] <= end)].copy()
        common = pheno1_loc[["SNP"]].merge(pheno2_loc[["SNP"]], on="SNP", how="inner").drop_duplicates()
        common_snps = set(common["SNP"].astype(str))
        pheno1_common = pheno1_full.loc[pheno1_full["SNP"].astype(str).isin(common_snps)].copy()
        pheno2_common = pheno2_full.loc[pheno2_full["SNP"].astype(str).isin(common_snps)].copy()
        pheno1_common = pheno1_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        pheno2_common = pheno2_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        locus_dir = out / r["locus"]
        locus_dir.mkdir(parents=True, exist_ok=True)
        pheno1_common.to_csv(locus_dir / f"gwas_{pheno1_prefix}.ldgwas.tsv", sep="\t", index=False)
        pheno2_common.to_csv(locus_dir / f"gwas_{pheno2_prefix}.ldgwas.tsv", sep="\t", index=False)
        common.sort_values("SNP").to_csv(locus_dir / "common_snps.tsv", sep="\t", index=False)
    return lead_snps_df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno1", required=True)
    ap.add_argument("--pheno2", required=True)
    ap.add_argument("--pheno1_prefix", required=True)
    ap.add_argument("--pheno2_prefix", required=True)
    ap.add_argument("--clump_path", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--window", type=int, default=500000)
    # Need to change this param to 250kb probs
    args = ap.parse_args()
    define_loci(args.pheno1,
                args.pheno2,
                args.pheno1_prefix,
                args.pheno2_prefix,
                args.clump_path,
                args.out_dir,
                args.window)

if __name__ == "__main__":
    main()