#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LD_CLUMP {

    tag "${meta.trait1}_${meta.trait2}_clump"
    publishDir "${params.outdir}/LD_clumping", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(gwas1), path(gwas2), path(hits)
    path clump_py
    path ref_bed
    path ref_bim
    path ref_fam

    output:
    tuple val(meta), path("${meta.trait1}-${meta.trait2}"), emit: loci

    script:
    """
    set -euo pipefail

    PYBIN="${params.pybin}"
    if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
        PYBIN="python3"
    fi

    REF_PREFIX=\$(basename "${ref_bed}" .bed)

    "\$PYBIN" "${clump_py}" \\
        --mode pairwise \\
        --hits "${hits}" \\
        --out-dir . \\
        --pheno1 "${meta.trait1}" \\
        --pheno2 "${meta.trait2}" \\
        --plink "${params.plink}" \\
        --ref-bfile "\$REF_PREFIX" \\
        --score-col "${params.clump_score_col}" \\
        --clump-r2 "${params.clump_r2}" \\
        --clump-kb "${params.clump_kb}" \\
        --merge-kb "${params.merge_kb}" \\
        --ld-r2 "${params.ld_r2}" \\
        --ld-kb "${params.ld_kb}" \\
        --sumstats "${gwas1}" "${gwas2}" \\
        --sumstats-names "${meta.trait1}" "${meta.trait2}"
    """
}