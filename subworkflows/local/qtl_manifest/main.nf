#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { GEN_QTL_MANIFEST } from '../../../modules/local/qtl_manifest.nf'

/*
meta_qtl = [
    qtl_ids        : params.qtl_ids,
    all_cell_types : params.all_cell_types,
    cell_type_list : params.cell_type_list
]

ch_qtl_manifest_input = Channel.of(
    tuple(
        meta_qtl,
        file(params.qtls),
        file("${projectDir}/bin/prepare_qtl_manifest.py"),
        "${params.outdir}/QTLManifest"
    )
)

GEN_QTL_MANIFEST(ch_qtl_manifest_input)

*/

workflow QTL_MANIFEST_PREP {

    main:

    meta_qtl = [
        qtl_ids        : params.qtl_ids,
        all_cell_types : params.all_cell_types,
        cell_type_list : params.cell_type_list
    ]

    qtls_file = file(params.qtls)
    manifest_py = file("${projectDir}/bin/prepare_qtl_manifest.py")

    ch_manifest_input = Channel.of(
        tuple(meta_qtl, qtls_file, manifest_py)
    )

    GEN_QTL_MANIFEST(ch_manifest_input)

    emit:
    manifest = GEN_QTL_MANIFEST.out.manifest
}