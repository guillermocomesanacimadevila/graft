#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
workflows/neurobridge/main_magma.nf
*/

include { MAGMA } from '../../modules/local/magma'

workflow {

    if( !params.pairs )
        error "Missing --pairs (trait pairs tsv)"

    // if LDSC

    ch_traits = Channel
        .fromPath(params.pairs)
        .splitCsv(header: true, sep: "\t")
        .flatMap { row ->
            [ row.trait1.toString().trim(), row.trait2.toString().trim() ]
        }
        .unique()
        .map { trait ->
            def meta = [:]
            meta.id  = trait
            tuple(trait, meta)
        }

    ch_ldsc = Channel
        .fromPath("${workflow.launchDir}/results/ldsc/*/*_ldsc/*.sumstats.gz", checkIfExists: true)
        .map { f ->
            def trait = f.getName().replaceFirst(/\.sumstats\.gz$/, "")
            tuple(trait, f)
        }

    ch_in = ch_traits
        .join(ch_ldsc)
        .map { trait, meta, sumstats ->
            tuple(meta, sumstats)
        }

    // default params & scr
    magma_bash = file("${workflow.launchDir}/bin/magma_genome_wide.sh")

    /* reference
    path g100eur_bed
    path g100eur_bim
    path g100eur_fam
    path g100eur_snploc
    path grch37_gene_loc
    path magma_bin
    */

    g100eur_bed     = file("${workflow.launchDir}/ref/magma/g1000_eur.bed")
    g100eur_bim     = file("${workflow.launchDir}/ref/magma/g1000_eur.bim")
    g100eur_fam     = file("${workflow.launchDir}/ref/magma/g1000_eur.fam")
    g100eur_snploc  = file("${workflow.launchDir}/ref/magma/g1000_eur.snp.loc")
    grch37_gene_loc = file("${workflow.launchDir}/ref/magma/NCBI37.3.gene.loc")
    magma_bin       = file("${workflow.launchDir}/ref/magma/software/linux/magma")

    MAGMA(
        ch_in,
        magma_bash,
        g100eur_bed,
        g100eur_bim,
        g100eur_fam,
        g100eur_snploc,
        grch37_gene_loc,
        magma_bin
    )
}