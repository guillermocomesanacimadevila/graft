#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MAP_TARGET_GENES } from '../../../modules/local/targets'

workflow TARGET_GENES_PREP {

    take:
    ch_pairs
    map_script
    gencode_file

    main:

    ch_lava_loci = ch_pairs.map { meta ->
        def locus_dir = file("${params.outdir}/LAVA/loci/${meta.trait1}_${meta.trait2}")
        tuple(meta, locus_dir, "LAVA")
    }

    ch_defined_loci = ch_pairs.map { meta ->
        def locus_dir = file("${params.outdir}/Defined_loci/${meta.trait1}_${meta.trait2}")
        tuple(meta, locus_dir, "conjFDR")
    }
    ch_loci = ch_lava_loci.mix(ch_defined_loci)

    MAP_TARGET_GENES(ch_loci, map_script, gencode_file)

    emit:
    target_genes = MAP_TARGET_GENES.out.target_genes
}