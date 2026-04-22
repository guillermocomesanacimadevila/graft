#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DEFINE_LOCI {

    tag "${meta.trait1}_${meta.trait2}_define_loci"

    input:
    tuple val(meta), path(gwas1), path(gwas2), path(clump_dir)
    path define_loci_py

    output:
    tuple val(meta), path("${meta.trait1}_${meta.trait2}"), emit: loci

    script:
    """
    set -euo pipefail

    PYBIN="${params.pybin}"
    if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
        PYBIN="python3"
    fi

    OUTDIR="${meta.trait1}_${meta.trait2}"
    mkdir -p "\$OUTDIR"

    "\$PYBIN" "${define_loci_py}" \\
        --pheno1 "${gwas1}" \\
        --pheno2 "${gwas2}" \\
        --pheno1_prefix "${meta.trait1}" \\
        --pheno2_prefix "${meta.trait2}" \\
        --clump_path "${clump_dir}" \\
        --out_dir "\$OUTDIR" \\
        --window "${params.window}"
    """
}