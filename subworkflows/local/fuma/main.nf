#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FUMA_PREP }      from '../../../modules/local/fuma_prep'
include { FUMA_POST_DIRS } from '../../../modules/local/fuma_post_dirs'

workflow PREP_FUMA_INPUT {

  take:
  ch_pairs_raw

  main:

  ch_defined = ch_pairs_raw
    .map { row ->
      def trait1   = row.trait1.toString().trim()
      def trait2   = row.trait2.toString().trim()
      def pair     = "${trait1}_${trait2}"
      def loci_dir = file("${params.outdir}/Defined_loci/${pair}")

      if( loci_dir.exists() ) {
        def meta = [
          trait1 : trait1,
          trait2 : trait2,
          pair   : pair,
          traits : [trait1, trait2],
          source : "defined_loci"
        ]
        tuple(meta, loci_dir)
      } else {
        null
      }
    }
    .filter { it != null }

  ch_lava = ch_pairs_raw
    .map { row ->
      def trait1   = row.trait1.toString().trim()
      def trait2   = row.trait2.toString().trim()
      def pair     = "${trait1}_${trait2}"
      def loci_dir = file("${params.outdir}/LAVA/loci/${pair}")

      if( loci_dir.exists() ) {
        def meta = [
          trait1 : trait1,
          trait2 : trait2,
          pair   : pair,
          traits : [trait1, trait2],
          source : "lava"
        ]
        tuple(meta, loci_dir)
      } else {
        null
      }
    }
    .filter { it != null }

  ch_pairs = ch_defined.mix(ch_lava)

  ch_fuma_ready = FUMA_PREP(ch_pairs).fuma_ready

  ch_post = ch_pairs.flatMap { meta, loci_dir ->
    meta.traits.collect { trait ->
      tuple(
        [
          pair   : meta.pair,
          trait  : trait,
          source : meta.source
        ],
        loci_dir
      )
    }
  }

  ch_post_dirs = FUMA_POST_DIRS(ch_post).post_dirs

  emit:
  fuma_ready = ch_fuma_ready
  post_dirs  = ch_post_dirs
}