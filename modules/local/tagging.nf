#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TAGGING {

    tag "tagging_ldak"

    input:
    path ldak_bin
    path plink_dir

    output:
    tuple path("${params.tag_prefix}.tagging"), path("${params.tag_prefix}.taglist"), emit: tagfiles

    script:
    """
    set -euo pipefail
    OUT="${params.tag_prefix}"

    chmod +x "${ldak_bin}" || true

    REF="${plink_dir}/1000G.EUR.QC"

    for chr in {1..22}; do
        "./${ldak_bin}" \\
          --calc-tagging \${OUT}.chr\${chr} \\
          --bfile "\${REF}.\${chr}" \\
          --power ${params.power} \\
          --chr \${chr}
    done

    ls \${OUT}.chr*.tagging > \${OUT}.taglist
    "./${ldak_bin}" --join-tagging \${OUT} --taglist \${OUT}.taglist
    """
}