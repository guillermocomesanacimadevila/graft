#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ADD_NEFF {

  tag "${meta.id}"

  input:
  tuple val(meta), path(ldsc)
  file(neff_script)

  output:
  tuple val(meta), path("${meta.id}_ldsc_ready_neff.tsv"), emit: ldsc_neff

  script:
  """
  set -euo pipefail

  PYBIN="${params.pybin}"
  if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
    PYBIN="python3"
  fi

  CASES="${meta.cases}"
  CONTROLS="${meta.controls}"

  if [ -z "\$CASES" ] || [ "\$CASES" = "null" ]; then
    CASES="0"
  fi

  if [ -z "\$CONTROLS" ] || [ "\$CONTROLS" = "null" ]; then
    CONTROLS="0"
  fi

  # If not case-control (or missing values), just pass through
  if [ "\$CASES" = "0" ] || [ "\$CONTROLS" = "0" ]; then
    cp "${ldsc}" "${meta.id}_ldsc_ready_neff.tsv"
    exit 0
  fi

  "\$PYBIN" "${neff_script}" \
    --in "${ldsc}" \
    --pheno_id "${meta.id}" \
    --outdir "." \
    --cases "\$CASES" \
    --controls "\$CONTROLS"

  cp "${meta.id}/post-qc/${meta.id}.ldsc_ready_neff.tsv" \
     "${meta.id}_ldsc_ready_neff.tsv"
  """
}