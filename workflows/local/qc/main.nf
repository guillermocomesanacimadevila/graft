#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_GWAS }  from '../../../modules/local/qc_gwas'
include { ADD_NEFF } from '../../../modules/local/add_neff'

workflow STAGE1_QC {

  neurobridge.ParamUtils.requireParam(params.input, "input")

  ch_in = Channel
    .fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
      def meta = neurobridge.GWASSamplesheet.buildMeta(row)
      tuple(meta, file(row.gwas))
    }

  qc_script   = file("${projectDir}/bin/qc_gwas.py")
  neff_script = file("${projectDir}/bin/compute_neff.py")

  ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
  ADD_NEFF(ch_qc, neff_script)
}