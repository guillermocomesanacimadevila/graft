#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GEN_LD_MATRIX {

    tag "${meta.pair}_${meta.source}_gen_ld_matrix"

    input:
    tuple val(meta), path(loci_dir)
    path gen_ld_bash
    path ref_chr_files, stageAs: 'refpanel/*'
    val  ref_chr_base

    output:
    tuple val(meta), path("${meta.pair}_${meta.source}"), emit: ld_ready

    script:
    """
    set -euo pipefail

    BASHBIN="bash"
    OUTDIR="${meta.pair}_${meta.source}"

    "\$BASHBIN" "${gen_ld_bash}" \\
        "${meta.trait1}" \\
        "${meta.trait2}" \\
        "${loci_dir}" \\
        "refpanel/${ref_chr_base}" \\
        "\$OUTDIR"
    """
}