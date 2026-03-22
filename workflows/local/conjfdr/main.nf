#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils
import neurobridge.GWASSamplesheet
import neurobridge.PairSamplesheet

include { PREP_CONJFDR_INPUT }         from '../../../subworkflows/local/conjfdr/main'
include { CONJFDR_DATA_PREP; CONJFDR } from '../../../modules/local/conjFDR'

workflow STAGE1_CONJFDR {

    main:
    ParamUtils.requireParam(params.input, "input")
    ParamUtils.requireParam(params.pairs, "pairs")

    ch_in = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = GWASSamplesheet.buildMeta(row)
            tuple(meta, file(row.gwas))
        }

    ch_pairs_cfdr = Channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = PairSamplesheet.buildMeta(row)
            tuple(meta.trait1, meta.trait2, meta)
        }

    /*
    refs / scripts
    */
    qc_script      = file("${projectDir}/bin/qc_gwas.py")
    neff_script    = file("${projectDir}/bin/compute_neff.py")
    cfdr_prep      = file("${projectDir}/bin/conjFDR_prep.py")
    cfdr_r         = file("${projectDir}/bin/conjFDR.R")
    conjfdr_refdir = file("${projectDir}/ref/conjFDR")
    ref_bfile      = "${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC"

    /*
    QC + NEFF reuse / prep
    */
    ch_sum = PREP_CONJFDR_INPUT(
        ch_in,
        qc_script,
        neff_script
    ).cfdr_input

    /*
    pairwise join
    */

    ch_conjfdr_with_t1 = ch_pairs_cfdr
        .join(ch_sum, by: 0)
        .map { trait1, trait2, meta, f1 ->
            tuple(trait2, meta, f1)
        }

    ch_conjfdr_in = ch_conjfdr_with_t1
        .join(ch_sum, by: 0)
        .map { trait2, meta, f1, f2 ->
            tuple(meta, f1, f2)
        }

    /*
    conjFDR prep
    */

    ch_harmonised = CONJFDR_DATA_PREP(
        ch_conjfdr_in,
        cfdr_prep
    ).harmonised_tsv

    /*
    conjFDR run
    */
    
    CONJFDR(
        ch_harmonised,
        cfdr_r,
        conjfdr_refdir,
        ref_bfile
    )
}