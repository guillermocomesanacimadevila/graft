#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { ALIGN_GWAS } from '../../../modules/local/align'

workflow PREP_ALIGNED_LDSC_INPUT {

  take:
  ch_sumstats_neff   
  ch_pairs         
  align_script

  main:

  ch_sum_by_id = ch_sumstats_neff
    .map { meta, f ->
      tuple(meta.id, meta, f)
    }

  ch_with_t1 = ch_pairs
    .join(ch_sum_by_id, by: 0)
    .map { trait1, trait2, pair_meta, meta1, f1 ->
      tuple(trait2, pair_meta, f1)
    }

  ch_align_in = ch_with_t1
    .join(ch_sum_by_id, by: 0)
    .map { trait2, pair_meta, f1, meta2, f2 ->
      tuple(pair_meta, f1, f2)
    }

  ch_aligned = ALIGN_GWAS(ch_align_in, align_script).aligned_pair

  emit:
  aligned_ldsc_input = ch_aligned
}