#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FUMA_POST_DIRS {

    tag "${meta.pair}_${meta.trait}_fuma_post"

    input:
    tuple val(meta), path(loci_dir)

    output:
    tuple val(meta), path("${meta.trait}"), emit: post_dirs

    script:
    """
    set -euo pipefail

    mkdir -p "${meta.trait}"

    for d in ${loci_dir}/locus_*; do
        [ -d "\$d" ] || continue
        locus_dir_name="\$(basename "\$d")"
        mkdir -p "${meta.trait}/\${locus_dir_name}"
    done
    """
}