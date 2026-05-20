#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
> GWAS coloc for:
> conjFDR loci (pleiotropic loci)
> LAVA loci (correlation-based loci)
*/

workflow PREP_GWAS_COLOC_INPUT {

  take:
  ch_pairs

  main:

  ch_defined = ch_pairs.map { row ->
    def trait1 = row.trait1.toString().trim()
    def trait2 = row.trait2.toString().trim()
    def loci_dir = file("${params.outdir}/Defined_loci/${trait1}_${trait2}")

    if( loci_dir.exists() ) {
      def meta = [
        trait1: trait1,
        trait2: trait2,
        type1 : row.type1 ? row.type1.toString().trim() : "cc",
        type2 : row.type2 ? row.type2.toString().trim() : "cc",
        s1    : row.s1    ? row.s1.toString().trim()    : "0.5",
        s2    : row.s2    ? row.s2.toString().trim()    : "0.5",
        source: "defined_loci"
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
        trait1: trait1,
        trait2: trait2,
        type1 : row.type1 ? row.type1.toString().trim() : "cc",
        type2 : row.type2 ? row.type2.toString().trim() : "cc",
        s1    : row.s1    ? row.s1.toString().trim()    : "0.5",
        s2    : row.s2    ? row.s2.toString().trim()    : "0.5",
        source: "lava"
      ]
      tuple(meta, loci_dir)
    } else {
      null
    }
  }.filter { it != null }

  ch_coloc_in = ch_defined.mix(ch_lava)

  emit:
  coloc_input = ch_coloc_in
}
