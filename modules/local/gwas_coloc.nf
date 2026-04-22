#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GWAS_COLOCALISATION {

    tag "${meta.trait1}_${meta.trait2}_${meta.source}_gwas_coloc"

    input:
    tuple val(meta), path(loci_dir)
    path gwas_coloc_r

    output:
    tuple val(meta), path("coloc_out/${meta.trait1}_${meta.trait2}_${meta.source}_coloc.tsv"), emit: coloc
    tuple val(meta), path("coloc_out"), emit: coloc_dir

    script:
    """
    set -euo pipefail

    RBIN="${params.rbin}"
    if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
        RBIN="Rscript"
    fi

    mkdir -p coloc_out

    "\$RBIN" "${gwas_coloc_r}" \\
        "${meta.trait1}_${meta.trait2}" \\
        "${meta.trait1}" \\
        "${meta.trait2}" \\
        "${loci_dir}" \\
        "coloc_out" \\
        "${meta.type1}" \\
        "${meta.type2}" \\
        "${meta.s1}" \\
        "${meta.s2}"

    if [ -f "coloc_out/${meta.trait1}_${meta.trait2}_coloc.tsv" ]; then
        mv "coloc_out/${meta.trait1}_${meta.trait2}_coloc.tsv" \\
           "coloc_out/${meta.trait1}_${meta.trait2}_${meta.source}_coloc.tsv"
    fi
    """
}