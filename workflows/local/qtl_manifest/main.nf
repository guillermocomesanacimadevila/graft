#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils

include { QTL_MANIFEST_PREP } from '../../../subworkflows/local/qtl_manifest/main'

workflow STAGE1_QTL_MANIFEST {

    main:
    ParamUtils.requireParam(params.qtls, "qtls")
    ParamUtils.requireParam(params.qtl_ids, "qtl_ids")

    if (!file(params.qtls).exists()) {
        error "QTL samplesheet not found: ${params.qtls}"
    }

    manifest_py = file("${projectDir}/bin/prepare_qtl_manifest.py")

    if (!manifest_py.exists()) {
        error "QTL manifest script not found: ${manifest_py}"
    }

    QTL_MANIFEST_PREP()

    emit:
    manifest = QTL_MANIFEST_PREP.out.manifest
}