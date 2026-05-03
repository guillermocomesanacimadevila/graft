#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SMR_BULK } from '../../../subworkflows/local/smr_bulk/main'
include { SMR_SC }   from '../../../subworkflows/local/smr_sc/main'

workflow STAGE1_SMR {

    take:
    qtl_manifest_ch

    main:
    neurobridge.ParamUtils.requireParam(params.input, "input")
    neurobridge.ParamUtils.requireParam(params.gencode_file, "gencode_file")
    neurobridge.ParamUtils.requireParam(params.smr_bfile, "smr_bfile")

    prep_qtl_script       = file("${projectDir}/bin/preprocess_qtl_for_smr.py")
    format_gwas_script    = file("${projectDir}/bin/format_gwas.py")
    grab_bulk_hits_script = file("${projectDir}/bin/grab_smr_hits_bulk.py")
    grab_sc_hits_script   = file("${projectDir}/bin/grab_smr_hits_sc.py")
    gtf_file              = file(params.gencode_file)

    qtl_manifest = qtl_manifest_ch.map { meta_qtl, manifest_file ->
        manifest_file
    }

    ch_gwas = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = neurobridge.GWASSamplesheet.buildMeta(row)
            tuple(meta, file(row.gwas))
        }

    ch_gwas_bulk = ch_gwas
    ch_gwas_sc   = ch_gwas

    SMR_BULK(
        ch_gwas_bulk,
        qtl_manifest,
        prep_qtl_script,
        format_gwas_script,
        grab_bulk_hits_script,
        gtf_file
    )

    SMR_SC(
        ch_gwas_sc,
        qtl_manifest,
        prep_qtl_script,
        format_gwas_script,
        grab_sc_hits_script,
        gtf_file
    )

    emit:
    bulk_results = SMR_BULK.out.smr_results
    bulk_hits    = SMR_BULK.out.smr_hits
    sc_results   = SMR_SC.out.smr_results
    sc_hits      = SMR_SC.out.smr_hits
}