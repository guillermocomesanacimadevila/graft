#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CONJFDR_DATA_PREP {

    tag "${meta.trait1}_${meta.trait2}_prep4_conjfdr"
    publishDir "${params.outdir}/conjFDR/ready/${meta.trait1}_${meta.trait2}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(gwas1), path(gwas2)
    path conjfdr_prep

    output:
    tuple val(meta), path("${meta.trait1}_${meta.trait2}.harmonised_${meta.trait1}_${meta.trait2}.tsv"), emit: harmonised_tsv

    script:
    """
    set -euo pipefail

    PYBIN="${params.pybin}"
    if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
        PYBIN="python3"
    fi

    "\$PYBIN" "${conjfdr_prep}" \\
        --trait1 "${gwas1}" \\
        --trait2 "${gwas2}" \\
        --prefix1 "${meta.trait1}" \\
        --prefix2 "${meta.trait2}" \\
        --out_prefix "${meta.trait1}_${meta.trait2}" \\
        --out_dir .
    """
}

process CONJFDR {

    tag "${meta.trait1}_${meta.trait2}_conjfdr"
    publishDir "${params.outdir}/conjFDR/res/${meta.trait1}_${meta.trait2}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(harmonised_tsv)
    path conjfdr_r
    path conjfdr_refdir

    output:
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_cfdr_results.tsv"), emit: cfdr_results
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_shared_hits.tsv"), optional: true, emit: shared_hits
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_shared_leads.tsv"), optional: true, emit: shared_leads
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_condFDR_${meta.trait1}_given_${meta.trait2}_0.01.tsv"), optional: true, emit: cond_trait1
    tuple val(meta), path("${meta.trait1}_${meta.trait2}_condFDR_${meta.trait2}_given_${meta.trait1}_0.01.tsv"), optional: true, emit: cond_trait2

    script:
    """
    set -euo pipefail

    RBIN="${params.rbin}"
    if [ -z "\$RBIN" ] || [ "\$RBIN" = "null" ]; then
        RBIN="Rscript"
    fi

    export PLINK="plink"
    export REF_BFILE="${params.conjfdr_ref_bfile}"

    "\$RBIN" "${conjfdr_r}" \\
        "${harmonised_tsv}" \\
        "${meta.trait1}_${meta.trait2}" \\
        "${conjfdr_refdir}" \\
        .
    """
}