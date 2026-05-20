#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow PREP_MAGMA_INPUT {

  take:
  ch_pairs

  main:

  ch_traits = ch_pairs
    .flatMap { row ->
      [ row.trait1.toString().trim(), row.trait2.toString().trim() ]
    }
    .unique()
    .map { trait ->
      def meta = [ id: trait ]
      tuple(trait, meta)
    }

  ch_ldsc = Channel
    .fromPath("${params.outdir}/LDSC/*/*.sumstats.gz", checkIfExists: true)
    .map { f ->
      def trait = f.getName().replaceFirst(/\.sumstats\.gz$/, "")
      tuple(trait, f)
    }

  ch_magma_in = ch_traits
    .join(ch_ldsc, by: 0)
    .map { trait, meta, sumstats ->
      tuple(meta, sumstats)
    }

  emit:
  magma_input = ch_magma_in
}