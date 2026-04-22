#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LDSC {

  tag "${meta.trait1}_${meta.trait2}"

  input:
  tuple val(meta), path(trait1_sumstats), path(trait2_sumstats)
  path ldsc_r
  path hm3_snplist
  path ld_chr_dir
  path wld_dir

  output:
  tuple val(meta), path("*"), emit: ldsc_outdir

  script:
  """
  set -euo pipefail

  RBIN="${params.rbin}"
  if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
    RBIN="Rscript"
  fi

  "\$RBIN" "${ldsc_r}" \
    "${trait1_sumstats}" \
    "${trait2_sumstats}" \
    "${meta.trait1}" \
    "${meta.trait2}" \
    "${meta.cases1}" \
    "${meta.controls1}" \
    "${meta.cases2}" \
    "${meta.controls2}" \
    "${meta.pop_prev1}" \
    "${meta.pop_prev2}" \
    "." \
    "${hm3_snplist}" \
    "${ld_chr_dir}" \
    "${wld_dir}" \
    "${params.ldsc_maf_filter}" \
    "${params.ldsc_info_filter}" \
    "${params.ldsc_mean_chisq_min}" \
    "${params.ldsc_intercept_min}" \
    "${params.ldsc_intercept_max}" \
    "${params.ldsc_h2_z_min}"
  """
}