#!/usr/bin/env python3
import argparse
import polars as pl
from pathlib import Path

# Instructions
# Look at LAVA result file
# FDR significant hits for trait pair X
# Grab coordinates
# Go to aligned GWAS for each trait and grab all SNPs within those coords
# Save ready WITH COLOC FORMATTING
# Save to results/LAVA/loci/trait1_trait2/locus_coords/...

def map_lava_loci(trait1: str, trait2: str, lava_file: str, out_dir: str, qc_dir: str):
    # AD_ldsc_ready_neff.tsv
    # ${params.outdir}/LAVA/${meta.trait1}_${meta.trait2}
    qc_dir = Path(qc_dir)
    out_dir = Path(out_dir)
    lava_path = Path(lava_file)
    trait1_file = qc_dir / trait1 / f"{trait1}_ldsc_ready_neff.tsv"
    trait2_file = qc_dir / trait2 / f"{trait2}_ldsc_ready_neff.tsv"

    if not trait1_file.exists():
        raise FileNotFoundError(f"Missing aligned GWAS for {trait1}: {trait1_file}")

    if not trait2_file.exists():
        raise FileNotFoundError(f"Missing aligned GWAS for {trait2}: {trait2_file}")

    if not lava_path.exists():
        raise FileNotFoundError(f"Yowza! LAVA result file not found: {lava_path}. Run LAVA!")

    trait1_sumstats = pl.read_csv(trait1_file, separator="\t")
    trait2_sumstats = pl.read_csv(trait2_file, separator="\t")
    df = pl.read_csv(lava_path, separator="\t")

    # here
    df2 = df.filter(pl.col("q_fdr") < 0.05)  # remember to keep this at < .05
    for row in df2.iter_rows(named=True):
        chrom = row["chr"]
        start = row["start"]
        end = row["stop"]
        locus_name = f"locus_chr{chrom}_{start}_{end}"
        locus_dir = Path(out_dir) / "LAVA" / "loci" / f"{trait1}_{trait2}" / locus_name
        locus_dir.mkdir(parents=True, exist_ok=True)

        trait1_locus = trait1_sumstats.filter(
            (pl.col("CHR") == chrom) &
            (pl.col("POS") >= start) &
            (pl.col("POS") <= end)
        )

        trait2_locus = trait2_sumstats.filter(
            (pl.col("CHR") == chrom) &
            (pl.col("POS") >= start) &
            (pl.col("POS") <= end)
        )

        trait1_locus.write_csv(locus_dir / f"gwas_{trait1}.ldgwas.tsv", separator="\t")
        trait2_locus.write_csv(locus_dir / f"gwas_{trait2}.ldgwas.tsv", separator="\t")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trait1", required=True)
    p.add_argument("--trait2", required=True)
    p.add_argument("--lava_file", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--qc_dir", required=True)
    args = p.parse_args()
    map_lava_loci(
        trait1=args.trait1,
        trait2=args.trait2,
        lava_file=args.lava_file,
        out_dir=args.out_dir,
        qc_dir=args.qc_dir,
    )

if __name__ == "__main__":
    main()