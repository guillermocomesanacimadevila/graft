#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MAP_TARGET_GENES {
    tag "${meta.trait1}_${meta.trait2}_${source}_target_genes"

    input:
    tuple val(meta), path(locus_dir), val(source)
    path map_script
    path gencode_file

    output:
    tuple val(meta), val(source), path("${meta.trait1}_${meta.trait2}_${source}_target_genes_${params.mapping_range}kb.tsv"), emit: target_genes

    script:
    """
    set -euo pipefail

    python ${map_script} \\
        --pheno1 ${meta.trait1} \\
        --pheno2 ${meta.trait2} \\
        --locus_dir ${locus_dir} \\
        --window ${params.mapping_range} \\
        --out_dir . \\
        --gencode_file ${gencode_file}

    mv TargetGenes/${meta.trait1}_${meta.trait2}/${source}/${meta.trait1}_${meta.trait2}_${source}_target_genes_${params.mapping_range}kb.tsv .
    """
}