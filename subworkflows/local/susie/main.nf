#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GEN_LD_MATRIX }     from '../../../modules/local/ld_matrix'
include { SUSIE_OVERLAP_MAP } from '../../../modules/local/susie'

workflow PREP_SUSIE_INPUT {

  take:
  ch_pairs
  gen_ld_bash
  ref_chr_files
  ref_chr_base

  main:

  ch_defined = ch_pairs.map { row ->
    def trait1 = row.trait1.toString().trim()
    def trait2 = row.trait2.toString().trim()
    def loci_dir = file("${params.outdir}/Defined_loci/${trait1}_${trait2}")
    
    if( loci_dir.exists() ) {
      def meta = [
        trait1   : trait1,
        trait2   : trait2,
        pair     : "${trait1}_${trait2}",
        source   : "defined_loci",
        cases1   : row.cases1,
        controls1: row.controls1,
        cases2   : row.cases2,
        controls2: row.controls2
      ]
      tuple(meta, loci_dir)
    } else {
      null
    }
  }.filter { it != null }

  ch_lava = ch_pairs.map { row ->
    def trait1 = row.trait1.toString().trim()
    def trait2 = row.trait2.toString().trim()
    def loci_dir = file("${params.outdir}/LAVA/loci/${trait1}_${trait2}")

    if( loci_dir.exists() ) {
      def meta = [
        trait1   : trait1,
        trait2   : trait2,
        pair     : "${trait1}_${trait2}",
        source   : "lava",
        cases1   : row.cases1,
        controls1: row.controls1,
        cases2   : row.cases2,
        controls2: row.controls2
      ]
      tuple(meta, loci_dir)
    } else {
      null
    }
  }.filter { it != null }

  ch_loci = ch_defined.mix(ch_lava)

  ch_ld_ready = GEN_LD_MATRIX(
    ch_loci,
    gen_ld_bash,
    ref_chr_files,
    ref_chr_base
  ).ld_ready

  ch_susie_in = ch_ld_ready
    .flatMap { meta, loci_dir ->
      def dir = loci_dir.toFile()

      def n1 = (
        meta.cases1.toString().toInteger() +
        meta.controls1.toString().toInteger()
      ).toString()

      def n2 = (
        meta.cases2.toString().toInteger() +
        meta.controls2.toString().toInteger()
      ).toString()

      def locus_dirs = []
      dir.eachDirRecurse { d ->
        if( d.name.startsWith("locus_chr") ) {
          locus_dirs << d
        }
      }

      locus_dirs.collect { locus_dir ->
        def locus  = locus_dir.name
        def coords = locus.replaceFirst(/^locus_/, "")
        def gwas1  = file("${locus_dir}/gwas_${meta.trait1}.ldorder.tsv")
        def gwas2  = file("${locus_dir}/gwas_${meta.trait2}.ldorder.tsv")
        def ld     = file("${locus_dir}/ld_gwas/locus.ld.gz")

        if( !gwas1.exists() || !gwas2.exists() || !ld.exists() ) {
          return null
        }

        tuple(
          [
            pair   : meta.pair,
            trait1 : meta.trait1,
            trait2 : meta.trait2,
            source : meta.source,
            locus  : locus,
            coords : coords,
            n1     : n1,
            n2     : n2
          ],
          gwas1,
          gwas2,
          ld
        )
      }.findAll { it != null }
    }

  emit:
  susie_input = ch_susie_in
}