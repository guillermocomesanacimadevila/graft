#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TWOSAMPLEMR_PREP } from '../../../subworkflows/local/twosamplemr/main'
include { TWOSAMPLEMR }      from '../../../modules/local/twosamplemr'

workflow STAGE1_TWOSAMPLEMR {

    main:
    neurobridge.ParamUtils.requireParam(params.input, "input")
    neurobridge.ParamUtils.requireParam(params.pairs, "pairs")
    neurobridge.ParamUtils.requireParam(params.outdir, "outdir")

    mr_script   = file("${projectDir}/bin/mr-pipeline.R")
    qc_script   = file("${projectDir}/bin/qc_gwas.py")
    neff_script = file("${projectDir}/bin/compute_neff.py")

    if( !mr_script.exists() ) {
        error "MR script not found: ${mr_script}"
    }
    if( !qc_script.exists() ) {
        error "QC script not found: ${qc_script}"
    }
    if( !neff_script.exists() ) {
        error "Neff script not found: ${neff_script}"
    }

    ch_in = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def meta = neurobridge.GWASSamplesheet.buildMeta(row)
            tuple(meta, file(row.gwas))
        }

    ch_pairs = Channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def meta = neurobridge.PairSamplesheet.buildMeta(row)
            meta.id = "${meta.trait1}_${meta.trait2}"
            meta
        }

    ch_mr_in = TWOSAMPLEMR_PREP(ch_in, ch_pairs, qc_script, neff_script).mr_input

    TWOSAMPLEMR(ch_mr_in, mr_script)
}