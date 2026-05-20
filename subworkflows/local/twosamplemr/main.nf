#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }  from '../../../modules/local/qc_gwas'
include { ADD_NEFF } from '../../../modules/local/add_neff'

workflow TWOSAMPLEMR_PREP {

    take:
    ch_in
    ch_pairs
    qc_script
    neff_script

    main:
    ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
    ch_neff = ADD_NEFF(ch_qc, neff_script)

    /*
      ch_neff should be:
      tuple val(meta), path(neff_file)
    */

    ch_trait_files = ch_neff.map { meta, neff_file ->
        tuple(meta.id, neff_file)
    }

    ch_exp = ch_pairs
        .map { meta -> tuple(meta.trait1, meta) }
        .join(ch_trait_files, by: 0)
        .map { trait1, pair_meta, exposure_file ->
            tuple(pair_meta.trait2, pair_meta, exposure_file)
        }

    ch_mr = ch_exp
        .join(ch_trait_files, by: 0)
        .map { trait2, pair_meta, exposure_file, outcome_file ->
            tuple(pair_meta, exposure_file, outcome_file)
        }

    emit:
    mr_input = ch_mr
    ldsc_ready = ch_qc
    ldsc_ready_neff = ch_neff
}