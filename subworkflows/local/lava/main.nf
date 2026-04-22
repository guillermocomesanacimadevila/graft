#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LAVA_GWAS_PREP } from '../../../modules/local/lava_prep'
include { MAKE_INFO_FILE } from '../../../modules/local/lava_prep'

workflow PREP_LAVA_INPUT {

  take:
  ch_qc_inputs      // tuple(meta, path(qc_neff_tsv))
  ch_pairs          // tuple(trait1, trait2, pair_meta)
  lava_data_prep

  main:

  ch_lava = LAVA_GWAS_PREP(ch_qc_inputs, lava_data_prep).lava_tsv

  ch_lava_key = ch_lava.map { meta, tsv ->
    tuple(meta.id, meta, tsv)
  }

  ch_pairs_with_t1 = ch_pairs
    .join(ch_lava_key, by: 0)
    .map { trait1, trait2, pair_meta, meta1, t1_tsv ->
      tuple(trait2, pair_meta, t1_tsv)
    }

  ch_info_in = ch_pairs_with_t1
    .join(ch_lava_key, by: 0)
    .map { trait2, pair_meta, t1_tsv, meta2, t2_tsv ->
      tuple(pair_meta, t1_tsv, t2_tsv)
    }

  ch_info = MAKE_INFO_FILE(ch_info_in).info_tsv

  ch_overlap = ch_pairs.map { trait1, trait2, meta ->
    def overlap_csv = file("${params.outdir}/LDSC/${trait1}_${trait2}/overlap_corr_for_LAVA_${trait1}_${trait2}.csv")
    tuple(meta, overlap_csv)
  }

  ch_lava_in = ch_info_in
    .join(ch_info, by: 0)
    .map { meta, t1_tsv, t2_tsv, info_tsv ->
      tuple(meta, t1_tsv, t2_tsv, info_tsv)
    }
    .join(ch_overlap, by: 0)
    .map { meta, t1_tsv, t2_tsv, info_tsv, overlap_csv ->
      tuple(meta, t1_tsv, t2_tsv, info_tsv, overlap_csv)
    }

  emit:
  lava_input = ch_lava_in
}