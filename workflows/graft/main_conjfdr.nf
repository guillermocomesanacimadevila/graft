#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }                    from '../../modules/local/qc_gwas'
include { ADD_NEFF }                   from '../../modules/local/add_neff'
include { CONJFDR_DATA_PREP; CONJFDR } from "../../modules/local/conjFDR"

workflow {

    // scripts ~~~~~~~~
    qc_script = file("${workflow.launchDir}/bin/qc_gwas.py")
    neff_script = file("${workflow.launchDir}/bin/compute_neff.py")
    conjfdr_prep = file("${workflow.launchDir}/bin/conjFDR_prep.py")
    conjfdr_r = file("${workflow.launchDir}/bin/conjFDR.R")

    // refs ~~~~~~~~
    conjfdr_refdir = file("${workflow.launchDir}/ref/conjFDR")
    ref_bfile = "${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC"

    if( !params.input )
        error "Missing --input (samplesheet tsv)"

    if( !params.pairs )
        error "Missing --pairs (pairs tsv)"

    ch_in = Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:'\t')
        .map { row ->

            def meta              = [:]
            meta.id               = row.id.toString().trim()
            meta.sep              = row.sep ? row.sep.replace('\\t','\t') : '\t'
            meta.snp_col          = row.snp_col
            meta.chr_col          = row.chr_col
            meta.pos_col          = row.pos_col
            meta.a1_col           = row.a1_col
            meta.a2_col           = row.a2_col
            meta.beta_col         = row.beta_col
            meta.se_col           = row.se_col
            meta.p_col            = row.p_col
            meta.eaf_col          = row.eaf_col
            meta.n_col            = row.n_col
            meta.info_col         = row.info_col
            meta.freq_case_col    = row.freq_case_col
            meta.freq_ctrl_col    = row.freq_ctrl_col
            meta.n_case_col       = row.n_case_col
            meta.n_ctrl_col       = row.n_ctrl_col
            meta.require_info     = (row.require_info ?: "false").toString()
            meta.exclude_mhc      = (row.exclude_mhc ?: "false").toString()
            meta.exclude_apoe     = (row.exclude_apoe ?: "false").toString()
            meta.drop_palindromes = (row.drop_palindromes ?: "false").toString()
            meta.keep_snps_only   = (row.keep_snps_only ?: "false").toString()
            meta.apoe_chr         = (row.apoe_chr ?: "19").toString()
            meta.apoe_start       = (row.apoe_start ?: "44000000").toString()
            meta.apoe_end         = (row.apoe_end ?: "46500000").toString()
            meta.cases            = (row.cases ?: "").toString().trim()
            meta.controls         = (row.controls ?: "").toString().trim()

            tuple(meta, file(row.gwas))
        }

    // QC + NEFF
    ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
    ch_neff = ADD_NEFF(ch_qc, neff_script).ldsc_neff
    ch_neff_key = ch_neff.map { meta, f -> tuple(meta.id, f) }

    ch_pairs = Channel
        .fromPath(params.pairs)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def meta = [
                trait1    : row.trait1.toString().trim(),
                trait2    : row.trait2.toString().trim(),
                cases1    : row.cases1,
                controls1 : row.controls1,
                cases2    : row.cases2,
                controls2 : row.controls2,
                pop_prev1 : row.pop_prev1,
                pop_prev2 : row.pop_prev2
            ]
            tuple(meta.trait1, meta.trait2, meta)
        }

    ch_conjfdr_with_t1 = ch_pairs
        .join(ch_neff_key, by: 0)
        .map { trait1, trait2, meta, f1 ->
            tuple(trait2, meta, f1)
        }

    ch_conjfdr_in = ch_conjfdr_with_t1
        .join(ch_neff_key, by: 0)
        .map { trait2, meta, f1, f2 ->
            tuple(meta, f1, f2)
        }

    ch_harmonised = CONJFDR_DATA_PREP(
        ch_conjfdr_in,
        conjfdr_prep
    ).harmonised_tsv

    CONJFDR(
        ch_harmonised,
        conjfdr_r,
        conjfdr_refdir,
        ref_bfile
    )
}

// nextflow run workflows/neurobridge/main_conjfdr.nf \
//   -profile docker \
//   -c conf/local/nextflow.config \
//   --input assets/gwas.tsv \
//   --pairs assets/ldsc_pairs.tsv \
//   --outdir results \
//   -resume