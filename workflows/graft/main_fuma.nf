#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FUMA_PREP }      from '../../modules/local/fuma_prep'
include { FUMA_POST_DIRS } from '../../modules/local/fuma_post_dirs'

/*
Checks => 
if defined_loci exists and is not empty
*/

workflow {

    if ( !params.pairs )
        error "Missing --pairs (Configure assets/ldsc_pairs.csv)"

    file(params.pairs)
        .splitCsv(header: true, sep: "\t")
        .each { row -> 
            def trait1   = row.trait1.toString().trim()
            def trait2   = row.trait2.toString().trim()
            def loci_dir = file("${params.outdir}/defined_loci/${trait1}_${trait2}/defined_loci")
            
            // check wether loci == mapped
            if ( !loci_dir.exists() ) {
                error "Missing defined_loci for ${trait1}_${trait2}: ${loci_dir} ; run main_clump.nf first!"
            } 
        }

    println "All defined_loci dirs found!"

    ch_pairs = Channel
        .fromPath(params.pairs)
        .splitCsv(header: true, sep: "\t")
        .map { row -> 
            def trait1   = row.trait1.toString().trim()
            def trait2   = row.trait2.toString().trim()
            def pair     = "${trait1}_${trait2}"
            def loci_dir = file("${params.outdir}/defined_loci/${pair}/defined_loci")
            def meta = [
                trait1 : trait1,
                trait2 : trait2,
                pair : pair,
                traits : [trait1, trait2]
            ]

            tuple(meta, loci_dir)
        }
        
        // parse channel onto FUMA_PREP module/
        FUMA_PREP(ch_pairs)

    /*
    > ASSEMBLE POST-FUMA DIRS
    */

    ch_post = ch_pairs.flatMap { meta, loci_dir ->
        meta.traits.collect { trait ->
            tuple(
                [
                    pair  : meta.pair,
                    trait : trait
                ],
                loci_dir
            )
        }
    }

    FUMA_POST_DIRS(ch_post)
}