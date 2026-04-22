#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils
import neurobridge.GWASSamplesheet
import neurobridge.PairSamplesheet

include { PREP_LAVA_INPUT }         from '../../../subworkflows/local/lava'
include { LAVA; LAVA_HITS_TO_LOCI } from '../../../modules/local/lava'

workflow STAGE1_LAVA {

  main:
  ParamUtils.requireParam(params.input, "input")
  ParamUtils.requireParam(params.pairs, "pairs")

  ch_in = Channel
    .fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
      def meta = GWASSamplesheet.buildMeta(row)
      def qc_tsv = file("${params.outdir}/QC/${meta.id}/${meta.id}_ldsc_ready_neff.tsv")

      if( !qc_tsv.exists() ) {
        error "Missing QC/NEFF file for ${meta.id}: ${qc_tsv}"
      }

      tuple(meta, qc_tsv)
    }

  ch_pairs = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
      def meta = PairSamplesheet.buildMeta(row)

      def overlap_csv = file("${params.outdir}/LDSC/${meta.trait1}_${meta.trait2}/overlap_corr_for_LAVA_${meta.trait1}_${meta.trait2}.csv")
      if( !overlap_csv.exists() ) {
        error "Missing LDSC overlap file for ${meta.trait1}_${meta.trait2}: ${overlap_csv}"
      }

      tuple(meta.trait1, meta.trait2, meta)
    }

  lava_data_prep = file("${projectDir}/bin/prep_data.py")
  lava_r         = file("${projectDir}/bin/lava_pair.R")
  lava_ref_dir   = file("${projectDir}/ref/lava_ref/lava_ref")
  loci_file      = file("${projectDir}/ref/lava_ref/lava_blocks.coords.loci")
  lava_hits_py   = file("${projectDir}/bin/lava_hits_to_loci.py")

  if( !lava_ref_dir.exists() ) {
    error "Missing LAVA reference directory: ${lava_ref_dir}"
  }

  if( !loci_file.exists() ) {
    error "Missing LAVA loci file: ${loci_file}"
  }

  if( !lava_hits_py.exists() ) {
    error "Missing LAVA loci-mapping script: ${lava_hits_py}"
  }

  ch_lava_in = PREP_LAVA_INPUT(ch_in, ch_pairs, lava_data_prep).lava_input

  ch_lava_out = LAVA(ch_lava_in, lava_r, lava_ref_dir, loci_file).lava_out

  LAVA_HITS_TO_LOCI(ch_lava_out, lava_hits_py)
}