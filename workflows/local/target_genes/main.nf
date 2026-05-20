#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TARGET_GENES_PREP } from '../../../subworkflows/local/target_genes/main'

workflow STAGE1_TARGET_GENES {

    main:
    neurobridge.ParamUtils.requireParam(params.pairs, "pairs")
    neurobridge.ParamUtils.requireParam(params.gencode_file, "gencode_file")
    neurobridge.ParamUtils.requireParam(params.outdir, "outdir")

    map_script = file("${projectDir}/bin/map_target_genes.py")
    gencode_file = file(params.gencode_file)

    if (!map_script.exists()) {
        error "Target gene mapping script not found: ${map_script}"
    }

    if (!gencode_file.exists()) {
        error "GENCODE file not found: ${gencode_file}"
    }

    ch_pairs = Channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = neurobridge.PairSamplesheet.buildMeta(row)
            meta.id = "${meta.trait1}_${meta.trait2}"
            meta
        }

    TARGET_GENES_PREP(ch_pairs, map_script, gencode_file)
}