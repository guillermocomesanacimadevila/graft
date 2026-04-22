#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils
import neurobridge.GWASSamplesheet
import neurobridge.PairSamplesheet

include { PREP_HDL_INPUT } from '../../../subworkflows/local/hdl/main'
include { HDL }            from '../../../modules/local/hdl' 

workflow STAGE1_HDL {

    main:
    ParamUtils.requireParam(params.input, "input")
    ParamUtils.requireParam(params.pairs, "pairs")

    ch_in = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = GWASSamplesheet.buildMeta(row)
            tuple(meta, file(row.gwas))
        }

    ch_pairs_hdl = Channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = PairSamplesheet.buildMeta(row)
            tuple(meta.trait1, meta.trait2, meta)
        }

    /*
    refs and scr
    */

    qc_script   = file("${projectDir}/bin/qc_gwas.py")
    neff_script = file("${projectDir}/bin/compute_neff.py")
    hdl_script  = file("${workflow.launchDir}/bin/hdl.R")
    ld_path     = file("${workflow.launchDir}/ref/HDL_ref/UKB_imputed_SVD_eigen99_extraction")
    
    // run checks for pre-HDL QC
    ch_sum  = PREP_HDL_INPUT(
        ch_in,
        qc_script,
        neff_script
    ).hdl_input

    ch_hdl_with_t1 = ch_pairs_hdl
        .join(ch_sum, by: 0)
        .map { trait1, trait2, meta, f1 ->
            def m = meta + [ beta1: "BETA", se1: "SE" ]
            tuple(trait2, m, f1)
        }
    
    ch_hdl_in = ch_hdl_with_t1
        .join(ch_sum, by: 0)
        .map { trait2, meta, f1, f2 ->
            def m = meta + [ beta2: "BETA", se2: "SE" ]
            tuple(m, f1, f2)
        }

    HDL(
        ch_hdl_in,
        hdl_script,
        ld_path
    )
}