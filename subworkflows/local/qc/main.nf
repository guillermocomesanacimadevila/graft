#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }    from '../../../modules/local/qc_gwas'
include { ADD_NEFF }   from '../../../modules/local/add_neff'

workflow QC_PIPELINE {

    take:
    ch_in
    qc_script
    neff_script

    main:
    ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
    ch_neff = ADD_NEFF(ch_qc, neff_script)

    emit:
    ldsc_ready = ch_qc
    ldsc_ready_neff = ch_neff
}