#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ALIGN_GWAS {

  tag "${meta.trait1}_vs_${meta.trait2}"

  input:
  tuple val(meta), path(gwas1), path(gwas2)
  file(align_script)

  output:
  tuple val(meta), path("${meta.trait1}_ldsc_ready_neff.tsv"), path("${meta.trait2}_ldsc_ready_neff.tsv"), emit: aligned_pair

  script:
  """
  set -euo pipefail

  PYBIN="${params.pybin}"
  if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
    PYBIN="python3"
  fi

  "\$PYBIN" "${align_script}" \
    --pheno1_id "${meta.trait1}" \
    --pheno2_id "${meta.trait2}" \
    --pheno1_gwas "${gwas1}" \
    --pheno2_gwas "${gwas2}" \
    --out_dir "."

  cp "QC/${meta.trait1}/${meta.trait1}.ldsc_ready_neff.tsv" "${meta.trait1}_ldsc_ready_neff.tsv"
  cp "QC/${meta.trait2}/${meta.trait2}.ldsc_ready_neff.tsv" "${meta.trait2}_ldsc_ready_neff.tsv"
  """
}