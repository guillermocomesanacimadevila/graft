#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
=> Modules to import
> GWAS_COLOC
> if ( !defined_loci )
>    error
*/

include { GWAS_COLOCALISATION } from '../../modules/local/gwas_coloc'

workflow {

    if ( !params.pairs ) {
        error "Missing --pairs (pairs tsv)"
    }

    file(params.pairs)
        .splitCsv(header:true, sep:"\t")
        .each { row ->
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def loci_dir = file("${params.outdir}/defined_loci/${trait1}_${trait2}/defined_loci")

            if( !loci_dir.exists() ) {
                error "Missing defined loci for ${trait1}_${trait2}: ${loci_dir} ; run main_clump.nf first!"
            }
        }

        println "All defined_loci dirs found!"

        // scr
        gwas_coloc_r = file("${workflow.launchDir}/bin/coloc.R")

        // need to add cols to pairs csv 
        ch_loci = Channel
            .fromPath(params.pairs)
            .splitCsv(header:true, sep:"\t")
            .map { row -> 
            def trait1 = row.trait1.toString().trim()
            def trait2 = row.trait2.toString().trim()
            def loci_dir = file("${params.outdir}/defined_loci/${trait1}_${trait2}/defined_loci")
            def meta = [
                trait1: trait1,
                trait2: trait2,
                type1 : row.type1 ? row.type1.toString().trim() : "cc",
                type2 : row.type2 ? row.type2.toString().trim() : "cc",
                s1 : row.s1 ? row.s1.toString().trim() : "0.5",
                s2 : row.s2 ? row.s2.toString().trim() : "0.5"
                ]
                tuple(meta, loci_dir)
            }
        
        GWAS_COLOCALISATION(
        ch_loci,
        gwas_coloc_r
    )
}