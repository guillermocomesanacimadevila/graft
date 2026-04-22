#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LD_CLUMP }    from '../../../modules/local/ld_clump'
include { DEFINE_LOCI } from '../../../modules/local/def_loci'

workflow PREP_DEFINED_LOCI_INPUT {

  take:
  ch_pairs
  clump_py
  define_loci_py
  ref_bed
  ref_bim
  ref_fam

  main:

  ch_clumped = LD_CLUMP(
    ch_pairs,
    clump_py,
    ref_bed,
    ref_bim,
    ref_fam
  ).loci

  ch_pairs_key = ch_pairs
    .map { meta, gwas1, gwas2, hits ->
      tuple("${meta.trait1}_${meta.trait2}", meta, gwas1, gwas2)
    }

  ch_clumped_key = ch_clumped
    .map { meta, clump_dir ->
      tuple("${meta.trait1}_${meta.trait2}", meta, clump_dir)
    }

  ch_define_loci_in = ch_pairs_key
    .join(ch_clumped_key, by: 0)
    .map { pair_id, meta_a, gwas1, gwas2, meta_b, clump_dir ->
      tuple(meta_a, gwas1, gwas2, clump_dir)
    }

  ch_defined = DEFINE_LOCI(
    ch_define_loci_in,
    define_loci_py
  ).loci

  emit:
  clumped_loci = ch_clumped
  defined_loci = ch_defined
}
