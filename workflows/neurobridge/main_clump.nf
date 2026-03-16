#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* 
> Need to check whether conjFDR outputs/ are there, if not print(run "main_conjfdr.nf first!")
*/

include { LD_CLUMP }    from '../../modules/local/ld_clump'
include { DEFINE_LOCI } from '../../modules/local/def_loci'

workflow {

    if ( !params.pairs ) {
        error "Missing --pairs (pairs tsv)"
    }

    file(params.pairs)
        .splitCsv(header:true, sep: "\t")
        .each { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def hits = file("${params.outdir}/conjFDR/${trait1}_${trait2}/${trait1}_${trait2}_shared_hits.tsv")

            if( !hits.exists() || hits.size() == 0 ) {
                error "Missing conjFDR hits for ${trait1}_${trait2}: ${hits} ; run main_conjfdr.nf first!"
            }
        }

    println "All conjFDR shared_hits files found!"

    // scr & ref
    clump_py = file("${workflow.launchDir}/bin/clump.py")
    define_loci_py = file("${workflow.launchDir}/bin/define_loci.py")
    ref_bed = file("${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bed")
    ref_bim = file("${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bim")
    ref_fam = file("${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.fam")

    ch_pairs = Channel
        .fromPath(params.pairs)
        .splitCsv(header:true, sep: "\t")
        .map { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def meta = [
                trait1: trait1,
                trait2: trait2
            ]
            def hits = file("${params.outdir}/conjFDR/${trait1}_${trait2}/${trait1}_${trait2}_shared_hits.tsv")
            def gwas1 = file("${params.outdir}/qc/${trait1}/${trait1}_ldsc_ready_neff.tsv")
            def gwas2 = file("${params.outdir}/qc/${trait2}/${trait2}_ldsc_ready_neff.tsv")

            tuple(meta, gwas1, gwas2, hits)
        }

    ch_clumped = LD_CLUMP(
        ch_pairs,
        clump_py,
        ref_bed,
        ref_bim,
        ref_fam
    ).loci

    ch_gwas_for_loci = Channel
        .fromPath(params.pairs)
        .splitCsv(header:true, sep: "\t")
        .map { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def meta = [
                trait1: trait1,
                trait2: trait2
            ]
            def pair_id = "${trait1}_${trait2}"
            def gwas1 = file("${params.outdir}/qc/${trait1}/${trait1}_ldsc_ready_neff.tsv")
            def gwas2 = file("${params.outdir}/qc/${trait2}/${trait2}_ldsc_ready_neff.tsv")

            tuple(pair_id, meta, gwas1, gwas2)
        }

    ch_clumped_key = ch_clumped
        .map { meta, clump_dir ->
            tuple("${meta.trait1}_${meta.trait2}", meta, clump_dir)
        }

    ch_define_loci_in = ch_gwas_for_loci
        .join(ch_clumped_key, by: 0)
        .map { pair_id, meta_a, gwas1, gwas2, meta_b, clump_dir ->
            tuple(meta_a, gwas1, gwas2, clump_dir)
        }

    DEFINE_LOCI(
        ch_define_loci_in,
        define_loci_py
    )
}