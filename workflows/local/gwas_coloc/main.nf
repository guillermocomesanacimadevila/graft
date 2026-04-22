#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils

include { PREP_GWAS_COLOC_INPUT } from '../../../subworkflows/local/gwas_coloc/main'
include { GWAS_COLOCALISATION }   from '../../../modules/local/gwas_coloc'

workflow STAGE1_GWAS_COLOC {

  main:
  ParamUtils.requireParam(params.pairs, "pairs")

  gwas_coloc_r = file("${projectDir}/bin/coloc.R")

  if( !gwas_coloc_r.exists() ) {
    error "Missing coloc script: ${gwas_coloc_r}"
  }

  ch_pairs = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: '\t')

  ch_coloc_in = PREP_GWAS_COLOC_INPUT(ch_pairs).coloc_input

  GWAS_COLOCALISATION(
    ch_coloc_in,
    gwas_coloc_r
  )
}