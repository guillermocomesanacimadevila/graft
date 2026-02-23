#!/usr/bin/env bash
set -euo pipefail

PHENO1_ID="${1:-}"
PHENO2_ID="${2:-}"
[ -n "${PHENO1_ID}" ] && [ -n "${PHENO2_ID}" ] || { echo "usage: bin/gen_ld_matrix.sh PHENO1 PHENO2"; exit 1; }
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOCI_DIR="${BASE_DIR}/outputs/defined_loci/${PHENO1_ID}_${PHENO2_ID}"

for LOCUS_DIR in "${LOCI_DIR}"/locus_chr*_*_*; do
  [ -d "${LOCUS_DIR}" ] || continue
  LOCUS_NAME="$(basename "${LOCUS_DIR}")"
  CHR_VAR="$(echo "${LOCUS_NAME}" | cut -d'_' -f2 | sed 's/^chr//')"
  cd "${LOCUS_DIR}"
  awk -F'\t' 'NR>1{print $1}' "gwas_${PHENO1_ID}.ldgwas.tsv" | sort -u > "${PHENO1_ID}.snps"
  awk -F'\t' 'NR>1{print $1}' "gwas_${PHENO2_ID}.ldgwas.tsv" | sort -u > "${PHENO2_ID}.snps"
  comm -12 "${PHENO1_ID}.snps" "${PHENO2_ID}.snps" > "gwas.snps.intersection"
  mkdir -p ld_gwas
  plink \
    --bfile "${BASE_DIR}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR_VAR}" \
    --extract "gwas.snps.intersection" \
    --maf 0.01 \
    --geno 0.02 \
    --make-bed \
    --out "ld_gwas/locus_tmp"
  cut -f2 "ld_gwas/locus_tmp.bim" > "ld_gwas/snps_in_ld_order.txt"

  plink \
    --bfile "ld_gwas/locus_tmp" \
    --r square gz \
    --out "ld_gwas/locus"

  python - <<PY
import pandas as pd
from pathlib import Path

order = pd.read_csv("ld_gwas/snps_in_ld_order.txt", header=None)[0].astype(str).tolist()
def reorder(f):
    f = Path(f)
    df = pd.read_csv(f, sep="\\t", dtype={"SNP": str})
    df = df.set_index("SNP")
    keep = [x for x in order if x in df.index]
    out = f.with_name(f.name.replace(".ldgwas.tsv", ".ldorder.tsv"))
    df.loc[keep].reset_index().to_csv(out, sep="\\t", index=False)

reorder("gwas_${PHENO1_ID}.ldgwas.tsv")
reorder("gwas_${PHENO2_ID}.ldgwas.tsv")
PY

  cd "${BASE_DIR}"
done
