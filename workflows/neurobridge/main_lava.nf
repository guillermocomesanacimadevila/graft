#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }        from '../../modules/local/qc_gwas'
include { ADD_NEFF }       from '../../modules/local/add_neff'
include { LAVA_GWAS_PREP } from '../../modules/local/lava_prep'
include { MAKE_INFO_FILE } from '../../modules/local/lava_prep'
include { LAVA }           from '../../modules/local/lava'

workflow {

  if( !params.input )
    error "Missing --input (samplesheet tsv)"

  if( !params.pairs )
    error "Missing --pairs (pairs tsv)"

  ch_in = Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t')
    .map { row ->

      def meta                = [:]
      meta.id                 = row.id.toString().trim()
      meta.sep                = row.sep ? row.sep.replace('\\t','\t') : '\t'
      meta.snp_col            = row.snp_col
      meta.chr_col            = row.chr_col
      meta.pos_col            = row.pos_col
      meta.a1_col             = row.a1_col
      meta.a2_col             = row.a2_col
      meta.beta_col           = row.beta_col
      meta.se_col             = row.se_col
      meta.p_col              = row.p_col
      meta.eaf_col            = row.eaf_col
      meta.n_col              = row.n_col
      meta.info_col           = row.info_col
      meta.freq_case_col      = row.freq_case_col
      meta.freq_ctrl_col      = row.freq_ctrl_col
      meta.n_case_col         = row.n_case_col
      meta.n_ctrl_col         = row.n_ctrl_col
      meta.require_info       = (row.require_info ?: "false").toString()
      meta.exclude_mhc        = (row.exclude_mhc ?: "false").toString()
      meta.exclude_apoe       = (row.exclude_apoe ?: "false").toString()
      meta.drop_palindromes   = (row.drop_palindromes ?: "false").toString()
      meta.keep_snps_only     = (row.keep_snps_only ?: "false").toString()
      meta.apoe_chr           = (row.apoe_chr ?: "19").toString()
      meta.apoe_start         = (row.apoe_start ?: "44000000").toString()
      meta.apoe_end           = (row.apoe_end ?: "46500000").toString()
      meta.cases              = (row.cases ?: "").toString().trim()
      meta.controls           = (row.controls ?: "").toString().trim()

      tuple(meta, file(row.gwas))
    }

  ch_pairs = Channel
    .fromPath(params.pairs)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
      def meta = [
        trait1: row.trait1.toString().trim(),
        trait2: row.trait2.toString().trim(),
        cases1: row.cases1,
        controls1: row.controls1,
        cases2: row.cases2,
        controls2: row.controls2,
        pop_prev1: row.pop_prev1,
        pop_prev2: row.pop_prev2
      ]
      tuple(meta.trait1, meta.trait2, meta)
    }

  qc_script = file("${workflow.launchDir}/bin/qc_gwas.py")
  neff_script = file("${workflow.launchDir}/bin/compute_neff.py")
  lava_data_prep = file("${workflow.launchDir}/bin/prep_data.py")
  lava_r = file("${workflow.launchDir}/bin/lava_pair.R")

  // Reference data
  lava_ref_dir = file("${workflow.launchDir}/ref/lava/lava_ref")
  lava_ref_prefix = "lava-ukb-v1.1"
  loci_file = file("${workflow.launchDir}/ref/lava/hdll_blocks.coords.loci")

  ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
  ch_neff = ADD_NEFF(ch_qc, neff_script).ldsc_neff

  ch_lava = LAVA_GWAS_PREP(
    ch_neff,
    lava_data_prep
  ).lava_tsv

  ch_lava_key = ch_lava.map { meta, tsv -> tuple(meta.id, tsv) }

  ch_pairs_with_t1 = ch_pairs
    .join(ch_lava_key, by: 0)
    .map { trait1, trait2, meta, t1_tsv ->
      tuple(trait2, meta, t1_tsv)
    }

  ch_info_in = ch_pairs_with_t1
    .join(ch_lava_key, by: 0)
    .map { trait2, meta, t1_tsv, t2_tsv ->
      tuple(meta, t1_tsv, t2_tsv)
    }

  ch_info = MAKE_INFO_FILE(ch_info_in).info_tsv

  // need to declare sample overlap matrix
  // overlap_corr_for_LAVA_SCZ_LON.csv
  ch_overlap = ch_pairs.map { trait1, trait2, meta ->
    def so_matrix = file("${workflow.launchDir}/results/ldsc/${trait1}_${trait2}/overlap_corr_for_LAVA_${trait1}_${trait2}.csv")
    tuple(meta, so_matrix)
  }

  ch_lava_in = ch_info_in
    .join(ch_info, by: 0)
    .map { meta, t1_tsv, t2_tsv, info_tsv ->
      tuple(meta, t1_tsv, t2_tsv, info_tsv)
    }
    .join(ch_overlap, by: 0)
    .map { meta, t1_tsv, t2_tsv, info_tsv, overlap_csv ->
      tuple(meta, t1_tsv, t2_tsv, info_tsv, overlap_csv)
    }

  LAVA(
  ch_lava_in,
  lava_r,
  lava_ref_dir,
  loci_file
  )  
}

// nextflow run workflows/neurobridge/main_lava.nf \
//   -profile docker \
//   -c conf/local/nextflow.config \
//   --input assets/gwas.tsv \
//   --pairs assets/ldsc_pairs.tsv \
//   --outdir results \
//   -resume