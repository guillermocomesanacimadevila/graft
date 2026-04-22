#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils

include { PREP_SUSIE_INPUT } from '../../../subworkflows/local/susie/main'
include { SUSIE_OVERLAP_MAP } from '../../../modules/local/susie'

workflow STAGE1_SUSIE {

  main:
  ParamUtils.requireParam(params.pairs, "pairs")

  ref_chr_files = Channel
    .fromPath("${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.*.{bed,bim,fam}")
    .collect()

  ref_chr_base = "1000G.EUR.QC"
  gen_ld_bash  = file("${projectDir}/bin/gen_ld_matrix.sh")
  susie_r      = file("${projectDir}/bin/susie.R")
  overlap_py   = file("${projectDir}/bin/finemap_causal_vars.py")
  map_gene_py  = file("${projectDir}/bin/map_credible_snps_to_closest_gene.py")
  gtf          = file("${projectDir}/ref/GENCODE/gencode.v37lift37.annotation.gtf")

  if( !gen_ld_bash.exists() ) error "Missing gen_ld_matrix.sh: ${gen_ld_bash}"
  if( !susie_r.exists() ) error "Missing susie.R: ${susie_r}"
  if( !overlap_py.exists() ) error "Missing finemap_causal_vars.py: ${overlap_py}"
  if( !map_gene_py.exists() ) error "Missing map_credible_snps_to_closest_gene.py: ${map_gene_py}"
  if( !gtf.exists() ) error "Missing GTF: ${gtf}"

  ch_pairs = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: "\t")

  ch_susie_in = PREP_SUSIE_INPUT(
    ch_pairs,
    gen_ld_bash,
    ref_chr_files,
    ref_chr_base
  ).susie_input

  SUSIE_OVERLAP_MAP(
    ch_susie_in,
    susie_r,
    overlap_py,
    map_gene_py,
    gtf,
    params.susie_L
  )
}