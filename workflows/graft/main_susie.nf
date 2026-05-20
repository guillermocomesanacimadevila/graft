#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GEN_LD_MATRIX }     from '../../modules/local/ld_matrix'
include { SUSIE_OVERLAP_MAP } from '../../modules/local/susie'

workflow {

    // defaults
    ref_chr_files = Channel
    .fromPath("${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.*.{bed,bim,fam}")

    ref_chr_base = "1000G.EUR.QC"
    gen_ld_bash = file("${workflow.launchDir}/bin/gen_ld_matrix.sh")
    susie_r = file("${workflow.launchDir}/bin/susie.R")
    overlap_py = file("${workflow.launchDir}/bin/finemap_causal_vars.py")
    map_gene_py = file("${workflow.launchDir}/bin/map_credible_snps_to_closest_gene.py")
    gtf = file("${workflow.launchDir}/ref/GENCODE/gencode.v37lift37.annotation.gtf")

    if ( !params.pairs) {
        error "Missing --pairs (pairs tsv)"
    }

    file(params.pairs)
        .splitCsv(header: true, sep: "\t")
        .each { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def loci_dir = file("${params.outdir}/defined_loci/${trait1}_${trait2}/defined_loci")

            if( !loci_dir.exists() ) {
                error "Missing defined_loci for ${trait1}_${trait2}: ${loci_dir} ; run main_clump.nf first!"
            }
        }

    println "All defined_loci dirs found!"

    ch_pairs = Channel
        .fromPath(params.pairs)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()

            def meta = [
                trait1 : trait1,
                trait2 : trait2,
                pair : "${trait1}_${trait2}",
                cases1 : row.cases1,
                controls1 : row.controls1,
                cases2 : row.cases2,
                controls2 : row.controls2
            ]

            def loci_dir = file("${params.outdir}/defined_loci/${trait1}_${trait2}/defined_loci")
            tuple(meta, loci_dir)
        }

    ch_ld_ready = GEN_LD_MATRIX(
        ch_pairs,
        gen_ld_bash,
        ref_chr_files.collect(),
        ref_chr_base
    ).ld_ready

    ch_susie_in = ch_ld_ready
        .flatMap { meta, loci_dir ->
            def dir = loci_dir.toFile()
            def n1 = (
                meta.cases1.toString().toInteger() +
                meta.controls1.toString().toInteger()
            ).toString()

            def n2 = (
                meta.cases2.toString().toInteger() +
                meta.controls2.toString().toInteger()
            ).toString()

            def locus_dirs = []
            dir.eachDirRecurse { d ->
                if( d.name.startsWith("locus_chr") ) {
                    locus_dirs << d
                }
            }

            locus_dirs.collect { locus_dir ->
                def locus  = locus_dir.name
                def coords = locus.replaceFirst(/^locus_/, "")
                def gwas1 = file("${locus_dir}/gwas_${meta.trait1}.ldorder.tsv")
                def gwas2 = file("${locus_dir}/gwas_${meta.trait2}.ldorder.tsv")
                def ld = file("${locus_dir}/ld_gwas/locus.ld.gz")

                if( !gwas1.exists() || !gwas2.exists() || !ld.exists() ) {
                    return null
                }

                tuple(
                    [
                        pair   : meta.pair,
                        trait1 : meta.trait1,
                        trait2 : meta.trait2,
                        locus  : locus,
                        coords : coords,
                        n1     : n1,
                        n2     : n2
                    ],
                    gwas1,
                    gwas2,
                    ld
                )
            }.findAll { it != null }
        }

    SUSIE_OVERLAP_MAP(
        ch_susie_in,
        susie_r,
        overlap_py,
        map_gene_py,
        gtf,
        params.susie_L
    )

    
}