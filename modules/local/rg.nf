#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RG {

    tag "${meta.trait1}-${meta.trait2}_SumHer_rg"

    input:
    tuple val(meta), path(sum1), path(sum2)
    path tagfile
    path ldak_bin
    path calc_p

    output:
    path "${meta.trait1}-${meta.trait2}.*", emit: rg_all
    path "res/*", emit: rg_res

    script:
    """
    set -euo pipefail
    mkdir -p res

    chmod +x "${ldak_bin}" || true

    "./${ldak_bin}" \\
      --sum-cors "${meta.trait1}-${meta.trait2}" \\
      --summary "${sum1}" \\
      --summary2 "${sum2}" \\
      --tagfile "${tagfile}" \\
      --cutoff ${params.thresh} \\
      --check-sums NO \\
      > "${meta.trait1}-${meta.trait2}.log" 2>&1

    ${params.pybin} "${calc_p}" \\
        --pheno1_prefix ${meta.trait1} \\
        --pheno2_prefix ${meta.trait2} \\
        --rg_column_name Value \\
        --se_column_name SE \\
        --results "${meta.trait1}-${meta.trait2}.cors" \\
        --out_dir res
    """
}