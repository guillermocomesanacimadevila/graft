#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TWOSAMPLEMR {
    tag "${meta.trait1}_${meta.trait2}_mr"

    input:
    tuple val(meta), path(exposure, stageAs: 'exposure.tsv'), path(outcome, stageAs: 'outcome.tsv')
    path mr_script

    output:
    tuple val(meta), path("forward_expToOut"), path("reverse_outToExp"), emit: mr_dir

    script:
    """
    set -euo pipefail

    if [ -z "${params.opengwas_jwt}" ] || [ "${params.opengwas_jwt}" = "null" ]; then
        echo "ERROR: --opengwas_jwt not provided"
        exit 1
    fi

    export OPENGWAS_JWT='${params.opengwas_jwt}'

    Rscript ${mr_script} \\
        --exposure ${exposure} \\
        --outcome ${outcome} \\
        --exposure_label ${meta.trait1} \\
        --outcome_label ${meta.trait2} \\
        --outdir . \\
        --pval_threshold ${params.mr_pval_threshold} \\
        --f_threshold ${params.f_stat_thresh} \\
        --clump_r2 ${params.mr_clump_r2} \\
        --clump_kb ${params.mr_clump_kb} \\
        --clump_p1 ${params.clump_p1} \\
        --clump_p2 ${params.clump_p2} \\
        --phenoscanner_p ${params.phenoscanner_p} \\
        --mrpresso_nb_distribution ${params.mrpresso_nb_distribution} \\
        --mrpresso_signif_threshold ${params.mrpresso_signif_threshold}
    """
}

// phenoscanner_p

/*
if [ -z "${params.opengwas_jwt}" ] || [ "${params.opengwas_jwt}" = "null" ]; then
        echo "ERROR: --opengwas_jwt not provided"
        exit 1
fi

    export OPENGWAS_JWT='${params.opengwas_jwt}'
*/
