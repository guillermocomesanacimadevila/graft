#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LAVA_GWAS_PREP {
    
    tag "${meta.id}_prep4_lava"

    input:
    tuple val(meta), path(gwas)
    path lava_data_prep

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: lava_tsv

    script:
    """
    set -euo pipefail

    PYBIN="${params.pybin}"
    if [ -z "\$PYBIN" ] || [ "\$PYBIN" = "null" ]; then
        PYBIN="python3"
    fi

    "\$PYBIN" "${lava_data_prep}" \\
        --input "${gwas}" \\
        --output "${meta.id}.tsv"
    """
}

process MAKE_INFO_FILE {

    tag "${meta.trait1}_${meta.trait2}_info_prep"
    publishDir "${params.outdir}/LAVA/info", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(t1_tsv), path(t2_tsv)

    output:
    tuple val(meta), path("${meta.trait1}_${meta.trait2}.input.info.txt"), emit: info_tsv

    script:
    """
    set -euo pipefail
    cat << EOF > "${meta.trait1}_${meta.trait2}.input.info.txt"
phenotype\tcases\tcontrols\tprevalence\tfilename
${meta.trait1}\t${meta.cases1}\t${meta.controls1}\t${meta.pop_prev1}\t${t1_tsv}
${meta.trait2}\t${meta.cases2}\t${meta.controls2}\t${meta.pop_prev2}\t${t2_tsv}
EOF
    """
}