#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process H2 {

    tag "${meta.id}_SumHer_h2"

    input:
    tuple val(meta), path(sumstats)
    path tagfile
    path ldak_bin

    output:
    tuple val(meta), path("${meta.id}.ldak.summaries"), emit: sums
    path "${meta.id}_h2.*", emit: h2
    path "${meta.id}_sumher.log", emit: log

    script:
    """
    set -euo pipefail

    chmod +x "${ldak_bin}" || true

    awk -v OFS='\\t' 'NR==1{print "Predictor","A1","A2","Z","N"; next} {print \$1,\$2,\$3,\$6/\$7,\$5}' \\
        ${sumstats} > ${meta.id}.ldak.summaries

    "./${ldak_bin}" \\
        --sum-hers ${meta.id}_h2 \\
        --summary ${meta.id}.ldak.summaries \\
        --tagfile "${tagfile}" \\
        --cutoff ${params.thresh} \\
        --check-sums NO \\
        > ${meta.id}_sumher.log 2>&1
    """
}