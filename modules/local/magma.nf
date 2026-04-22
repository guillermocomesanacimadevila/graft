#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MAGMA {

    tag "${meta.id}_magma"

    input:
    tuple val(meta), path(sumstats)
    path magma_bash
    path g100eur_bed
    path g100eur_bim
    path g100eur_fam
    path g100eur_snploc
    path grch37_gene_loc
    path magma_bin

    output:
    path("${meta.id}_magma.genes.out"), emit: genes_out
    path("${meta.id}_magma.genes.mapped.tsv"), emit: mapped
    path("${meta.id}_for_magma.txt"), emit: snp_pvals
    path("${meta.id}_01magma.done"), emit: done

    script:
    def bfile_prefix = "g1000_eur"

    """
    set -euo pipefail

    chmod +x ${magma_bash}
    chmod +x ${magma_bin}

    MAGMA_BIN="./${magma_bin.getName()}" \\
    bash ${magma_bash} \\
      . \\
      "${bfile_prefix}" \\
      "${g100eur_snploc}" \\
      "${grch37_gene_loc}" \\
      "${meta.id}" \\
      "${sumstats}"

    touch "${meta.id}_01magma.done"
    """
}