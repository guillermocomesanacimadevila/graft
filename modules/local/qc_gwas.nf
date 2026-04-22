#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process QC_GWAS {
  
  tag "${meta.id}"

  input:
  tuple val(meta), path(gwas)
  file(qc_script)

  output:
  tuple val(meta), path("${meta.id}_ldsc_ready.tsv"), emit: ldsc_ready

  script:
  """
  set -euo pipefail

  PYBIN="${params.pybin}"
  if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
    PYBIN="python3"
  fi

  MAF_MIN="${params.maf_min}"
  if [ -z "\$MAF_MIN" ] || [ "\$MAF_MIN" = "null" ]; then
    MAF_MIN="0.01"
  fi

  INFO_MIN="${params.info_min}"
  if [ -z "\$INFO_MIN" ] || [ "\$INFO_MIN" = "null" ]; then
    INFO_MIN="0.90"
  fi

  OUTDIR="${params.outdir}"
  if [ -z "\$OUTDIR" ] || [ "\$OUTDIR" = "null" ]; then
    OUTDIR="results"
  fi

  EXTRA=""

  if [ -n "${meta.eaf_col}" ] && [ "${meta.eaf_col}" != "null" ]; then
    EXTRA="\$EXTRA --eaf_col ${meta.eaf_col}"
  fi

  if [ -n "${meta.n_col}" ] && [ "${meta.n_col}" != "null" ]; then
    EXTRA="\$EXTRA --n_col ${meta.n_col}"
  fi

  if [ -n "${meta.info_col}" ] && [ "${meta.info_col}" != "null" ]; then
    EXTRA="\$EXTRA --info_col ${meta.info_col}"
  fi

  if [ -n "${meta.freq_case_col}" ] && [ "${meta.freq_case_col}" != "null" ]; then
    EXTRA="\$EXTRA --freq_case_col ${meta.freq_case_col}"
  fi

  if [ -n "${meta.freq_ctrl_col}" ] && [ "${meta.freq_ctrl_col}" != "null" ]; then
    EXTRA="\$EXTRA --freq_ctrl_col ${meta.freq_ctrl_col}"
  fi

  if [ -n "${meta.n_case_col}" ] && [ "${meta.n_case_col}" != "null" ]; then
    EXTRA="\$EXTRA --n_case_col ${meta.n_case_col}"
  fi

  if [ -n "${meta.n_ctrl_col}" ] && [ "${meta.n_ctrl_col}" != "null" ]; then
    EXTRA="\$EXTRA --n_ctrl_col ${meta.n_ctrl_col}"
  fi

  if [ "${meta.require_info}" = "true" ]; then
    EXTRA="\$EXTRA --require_info"
  fi

  if [ "${meta.exclude_mhc}" = "true" ]; then
    EXTRA="\$EXTRA --exclude_mhc"
  fi

  if [ "${meta.exclude_apoe}" = "true" ]; then
    EXTRA="\$EXTRA --exclude_apoe --apoe_chr ${meta.apoe_chr} --apoe_start ${meta.apoe_start} --apoe_end ${meta.apoe_end}"
  fi

  if [ "${meta.drop_palindromes}" = "true" ]; then
    EXTRA="\$EXTRA --drop_palindromes"
  fi

  if [ "${meta.keep_snps_only}" = "true" ]; then
    EXTRA="\$EXTRA --keep_snps_only"
  fi

  "\$PYBIN" "${qc_script}" \
    --in "${gwas}" \
    --pheno_id "${meta.id}" \
    --outdir "." \
    --sep '${meta.sep}' \
    --snp_col "${meta.snp_col}" \
    --chr_col "${meta.chr_col}" \
    --pos_col "${meta.pos_col}" \
    --a1_col "${meta.a1_col}" \
    --a2_col "${meta.a2_col}" \
    --beta_col "${meta.beta_col}" \
    --se_col "${meta.se_col}" \
    --p_col "${meta.p_col}" \
    --maf_min "\$MAF_MIN" \
    --info_min "\$INFO_MIN" \
    \$EXTRA

  cp "${meta.id}/post-qc/${meta.id}_ldsc_ready.tsv" "${meta.id}_ldsc_ready.tsv"
  """
}