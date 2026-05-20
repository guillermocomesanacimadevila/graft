#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// nextflow run workflows/neurobridge/main_qc.nf -profile docker -c conf/local/nextflow.config --input assets/gwas.tsv --outdir results

include { QC_GWAS }  from '../../modules/local/qc_gwas'
include { ADD_NEFF } from '../../modules/local/add_neff'

workflow {
  if( !params.input )
    error "Missing --input (samplesheet tsv)"

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
      meta.cases              = row.cases
      meta.controls           = row.controls
      tuple(meta, file(row.gwas))
    }

  qc_script = file("${workflow.launchDir}/bin/qc_gwas.py")
  neff_script = file("${workflow.launchDir}/bin/compute_neff.py")

  ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
  ADD_NEFF(ch_qc, neff_script)
}