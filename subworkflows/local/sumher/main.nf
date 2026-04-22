#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TAGGING } from '../../../modules/local/tagging'
include { H2 }      from '../../../modules/local/h2'
include { RG }      from '../../../modules/local/rg'

workflow SUMHER_RUN {

    take:
    ch_sumstats
    ch_pairs
    plink_dir
    calc_p
    ldak_bin

    main:
    /*
    TAGGING
    */
    def tag_path = "${params.outdir}/SumHer/tagging/${params.tag_prefix}.tagging"
    def taglist_path = "${params.outdir}/SumHer/tagging/${params.tag_prefix}.taglist"

    def tag_file = new File(tag_path)
    def taglist_file = new File(taglist_path)

    if( !params.do_tagging && tag_file.exists() && tag_file.length() > 0 && taglist_file.exists() && taglist_file.length() > 0 ) {
        println "[SKIP] tagging -> using existing ${tag_path}"
        ch_tag = Channel.of( tuple(file(tag_path), file(taglist_path)) )
    }
    else {
        println "[RUN ] tagging -> running LDAK tagging"
        ch_tag = TAGGING(
            ldak_bin,
            plink_dir
        ).tagfiles
    }

    ch_tagfile = ch_tag.map { tagfile, taglist -> tagfile }

    /*
    H2 per trait
    */
    ch_h2_in = ch_sumstats
        .map { meta, sumstats ->
            tuple(meta, sumstats)
        }

    ch_h2_sums = H2(
        ch_h2_in,
        ch_tagfile,
        ldak_bin
    ).sums

    /*
    RG per pair
    */
    ch_h2_key = ch_h2_sums
        .map { meta, sums ->
            tuple(meta.id, sums)
        }

    ch_rg_with_t1 = ch_pairs
        .join(ch_h2_key, by: 0)
        .map { trait1, trait2, meta, sum1 ->
            tuple(trait2, meta, sum1)
        }

    ch_rg_in = ch_rg_with_t1
        .join(ch_h2_key, by: 0)
        .map { trait2, meta, sum1, sum2 ->
            tuple(meta, sum1, sum2)
        }

    RG(
        ch_rg_in,
        ch_tagfile,
        ldak_bin,
        calc_p
    )

    emit:
    tagging = ch_tag
    h2      = ch_h2_sums
}