#!/usr/bin/env python3
import argparse
import pandas as pd
import re
import os 
from pathlib import Path

# for locus in dir 
# grab its coordinates 
# window param - from the ends +/- X kb
# Grab those genes and append onto a per locus list 
# save list as txt (sep="\t")
# results/TargetGenes/LAVA/pheno1_pheno2/coords/txt
# resutls/TargetGenes/defined_loci/pheno1_pheno2/coods/txt

def load_gencode_pc_gtf(gtf_path: str):
    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=["seqname","source","feature","start","end","score","strand","frame","attributes"],
        low_memory=False
    )
    gtf = gtf[gtf["feature"].eq("gene")].copy()
    attrs = gtf["attributes"].astype(str).str.strip().str.split(";")
    parsed = []
    for row in attrs:
        d = {}
        for part in row:
            part = part.strip()
            if part and " " in part:
                k, v = part.split(" ", 1)
                d[k] = v.strip().strip('"')
        parsed.append(d)
    gtf["gene_type"] = [x.get("gene_type", x.get("gene_biotype")) for x in parsed]
    gtf["gene_name"] = [x.get("gene_name") for x in parsed]
    gtf = gtf[gtf["gene_type"].eq("protein_coding")].copy()
    gtf["gene_name"] = gtf["gene_name"].astype(str).str.strip()
    gtf = gtf[gtf["gene_name"].ne("") & gtf["gene_name"].str.lower().ne("nan")].copy()
    gtf["gene_name"] = gtf["gene_name"].replace({"FAM63B": "MINDY1", "FAM63A": "MINDY2"})
    gtf["seqname"] = gtf["seqname"].astype(str).str.replace("^chr", "", regex=True)
    gtf["start"] = pd.to_numeric(gtf["start"], errors="coerce")
    gtf["end"] = pd.to_numeric(gtf["end"], errors="coerce")
    gtf = gtf[["seqname","start","end","gene_name"]].dropna().sort_values(["seqname","start","end"]).reset_index(drop=True)
    return gtf

def map_genes(pheno1: str, pheno2: str, locus_dir: str, window: int, out_dir: str, gencode_file: str):
    out_dir = Path(out_dir)
    locus_dir = Path(locus_dir)
    # dir = locus_dir / f"{pheno1}_{pheno2}"
    dir = Path(locus_dir)
    gencode_ref = load_gencode_pc_gtf(gencode_file)

    loci = []
    files = os.listdir(dir)

    if len(files) == 0:
        print("Yowza! No files round here")

    if "lead_snps.tsv" in files:
        source = "conjFDR"
    else:
        source = "LAVA"

    for file in files:
        if file.startswith("locus_chr"):
            loci.append(file)

    loci = [i.replace("locus_chr", "").replace(".txt", "").replace(".tsv", "").strip() for i in loci]
    mapped = []

    for locus in loci:
        locus_clean = locus.replace("chr", "").strip()
        if locus_clean == ".DS_Store":
            continue

        bits = locus_clean.split("_")

        if len(bits) != 3:
            print(f"Skipping weird locus name: {locus}")
            continue

        chrom = bits[0]
        start = int(bits[1])
        end = int(bits[2])
        window_bp = window * 1000
        lower = max(0, start - window_bp)
        upper = end + window_bp

        genes = gencode_ref[
            (gencode_ref["seqname"].astype(str).eq(str(chrom))) &
            (gencode_ref["start"].le(upper)) &
            (gencode_ref["end"].ge(lower))
        ].copy()

        gene_names = sorted(genes["gene_name"].dropna().unique())
        mapped.append({
            "source": source,
            "locus": locus,
            "chr": chrom,
            "start": start,
            "end": end,
            "window_kb": window,
            "window_start": lower,
            "window_end": upper,
            "genes": ",".join(gene_names),
            "n_genes": len(gene_names)
        })

    expected_cols = [
        "source", "locus", "chr", "start", "end",
        "window_kb", "window_start", "window_end",
        "genes", "n_genes"
    ]

    mapped = pd.DataFrame(mapped)

    if mapped.empty:
        mapped = pd.DataFrame(columns=expected_cols)
    else:
        mapped = mapped[expected_cols]

    target_dir = out_dir / "TargetGenes" / f"{pheno1}_{pheno2}" / source
    target_dir.mkdir(parents=True, exist_ok=True)
    out_file = target_dir / f"{pheno1}_{pheno2}_{source}_target_genes_{window}kb.tsv"
    mapped.to_csv(out_file, sep="\t", index=False)
    return mapped

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1", required=True)
    parser.add_argument("--pheno2", required=True)
    parser.add_argument("--locus_dir", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--gencode_file", required=True)
    args = parser.parse_args()
    map_genes(
        pheno1=args.pheno1,
        pheno2=args.pheno2,
        locus_dir=args.locus_dir,
        window=args.window,
        out_dir=args.out_dir,
        gencode_file=args.gencode_file
    )

if __name__ == "__main__":
    main()