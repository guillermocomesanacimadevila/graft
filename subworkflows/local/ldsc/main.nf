#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }  from '../../../modules/local/qc_gwas'
include { ADD_NEFF } from '../../../modules/local/add_neff'

workflow PREP_LDSC_INPUT {

  take:
  ch_in
  qc_script
  neff_script

  main:
  ch_in
    .map { meta, gwas ->

      def out_path = "${params.outdir}/QC/${meta.id}/${meta.id}.ldsc_ready_neff.tsv"
      def out_file = new File(out_path)
      def exists   = out_file.exists() && out_file.length() > 0

      if( exists ) {
        println "[SKIP] ${meta.id} -> using existing ${out_path}"
        [ 'done', tuple(meta, file(out_path)) ]
      } else {
        println "[RUN ] ${meta.id} -> running QC_GWAS + ADD_NEFF"
        [ 'missing', tuple(meta, gwas) ]
      }
    }
    .set { ch_checked }

  ch_checked
    .filter { tag, val -> tag == 'done' }
    .map    { tag, val -> val }
    .set { ch_existing }

  ch_checked
    .filter { tag, val -> tag == 'missing' }
    .map    { tag, val -> val }
    .set { ch_missing }

  ch_qc   = QC_GWAS(ch_missing, qc_script).ldsc_ready
  ch_neff = ADD_NEFF(ch_qc, neff_script).ldsc_neff

  ch_sum = ch_existing.mix(ch_neff)

  emit:
  ldsc_input = ch_sum
}