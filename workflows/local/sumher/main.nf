#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils
import neurobridge.GWASSamplesheet
import neurobridge.PairSamplesheet

include { PREP_SUMHER_INPUT } from '../../../subworkflows/local/sumher/main_prep'
include { SUMHER_RUN }        from '../../../subworkflows/local/sumher/main'

workflow STAGE1_SUMHER {

    main:
    ParamUtils.requireParam(params.input, "input")
    ParamUtils.requireParam(params.pairs, "pairs")

    ch_in = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = GWASSamplesheet.buildMeta(row)
            tuple(meta, file(row.gwas))
        }

    ch_pairs = Channel
        .fromPath(params.pairs, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = PairSamplesheet.buildMeta(row)
            tuple(meta.trait1, meta.trait2, meta)
        }

    qc_script   = file("${projectDir}/bin/qc_gwas.py")
    neff_script = file("${projectDir}/bin/compute_neff.py")
    calc_p      = file("${projectDir}/bin/calc_p.py")

    plink_dir   = file("${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink")
    ldak_bin    = file("${projectDir}/ref/SumHer/LDAK/${params.ldak_os}")

    ch_sumstats = PREP_SUMHER_INPUT(
        ch_in,
        qc_script,
        neff_script
    ).sumher_input

    SUMHER_RUN(
        ch_sumstats,
        ch_pairs,
        plink_dir,
        calc_p,
        ldak_bin
    )
}