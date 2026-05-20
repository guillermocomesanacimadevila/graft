#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process HDL {

  tag "${meta.trait1}_${meta.trait2}_hdl"

  input:
  tuple val(meta), path(trait1_sumstats), path(trait2_sumstats)
  file(hdl_r)
  path(ld_path)

  output:
  tuple val(meta), path("${meta.trait1}.hdl.tsv"), emit: hdl_trait1
  tuple val(meta), path("${meta.trait2}.hdl.tsv"), emit: hdl_trait2
  tuple val(meta), path("${meta.trait1}_${meta.trait2}.global_rg.tsv"), emit: hdl_rg
  tuple val(meta), path("${meta.trait1}_${meta.trait2}.global_fit.rds"), emit: hdl_fit

  script:
  """
  set -euo pipefail

  RBIN="${params.rbin}"
  if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
    RBIN="Rscript"
  fi

  "\$RBIN" "${hdl_r}" \
    "${trait1_sumstats}" \
    "${trait2_sumstats}" \
    "${meta.trait1}" \
    "${meta.trait2}" \
    "${meta.beta1}" \
    "${meta.se1}" \
    "${meta.beta2}" \
    "${meta.se2}" \
    "${ld_path}" \
    "." \
    "${params.hdl_eigen_cut}" \
    "${params.hdl_N0}" \
    "${params.hdl_lim}" \
    "${params.hdl_jackknife}"
  """
}