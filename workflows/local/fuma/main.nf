#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREP_FUMA_INPUT } from '../../../subworkflows/local/fuma/main'

workflow STAGE1_FUMA {

  main:
  neurobridge.ParamUtils.requireParam(params.pairs, "pairs")

  file(params.pairs)
    .splitCsv(header: true, sep: "\t")
    .each { row ->
      def trait1      = row.trait1.toString().trim()
      def trait2      = row.trait2.toString().trim()
      def pair        = "${trait1}_${trait2}"
      def defined_dir = file("${params.outdir}/Defined_loci/${pair}")
      def lava_dir    = file("${params.outdir}/LAVA/loci/${pair}")

      if( !defined_dir.exists() && !lava_dir.exists() ) {
        error "Missing loci for ${pair}: neither ${defined_dir} nor ${lava_dir} exists ; run STAGE1_DEF_LOCI and/or STAGE1_LAVA first!"
      }
    }

  println "At least one loci source found for all pairs!"

  ch_pairs_raw = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: "\t")

  PREP_FUMA_INPUT(ch_pairs_raw)
}