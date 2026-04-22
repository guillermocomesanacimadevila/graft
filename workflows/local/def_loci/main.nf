#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import neurobridge.ParamUtils
import neurobridge.PairSamplesheet

include { PREP_DEFINED_LOCI_INPUT } from '../../../subworkflows/local/def_loci/main'

workflow STAGE1_DEF_LOCI {

  main:
  ParamUtils.requireParam(params.pairs, "pairs")

  ch_pairs = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
      def meta = PairSamplesheet.buildMeta(row)

      def hits = file("${params.outdir}/ConjFDR/${meta.trait1}_${meta.trait2}/${meta.trait1}_${meta.trait2}_shared_hits.tsv")
      if( !hits.exists() || hits.size() == 0 ) {
        error "Missing conjFDR hits for ${meta.trait1}_${meta.trait2}: ${hits} ; run STAGE1_CONJFDR first!"
      }

      def gwas1 = file("${params.outdir}/QC/${meta.trait1}/${meta.trait1}_ldsc_ready_neff.tsv")
      def gwas2 = file("${params.outdir}/QC/${meta.trait2}/${meta.trait2}_ldsc_ready_neff.tsv")

      if( !gwas1.exists() ) {
        error "Missing QC/NEFF file for ${meta.trait1}: ${gwas1}"
      }

      if( !gwas2.exists() ) {
        error "Missing QC/NEFF file for ${meta.trait2}: ${gwas2}"
      }

      tuple(meta, gwas1, gwas2, hits)
    }

  clump_py       = file("${projectDir}/bin/clump.py")
  define_loci_py = file("${projectDir}/bin/define_loci.py")

  ref_bed = file("${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bed")
  ref_bim = file("${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bim")
  ref_fam = file("${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.fam")

  if( !clump_py.exists() ) {
    error "Missing LD clumping script: ${clump_py}"
  }

  if( !define_loci_py.exists() ) {
    error "Missing define_loci script: ${define_loci_py}"
  }

  if( !ref_bed.exists() || !ref_bim.exists() || !ref_fam.exists() ) {
    error "Missing PLINK reference files under ${projectDir}/ref/ldsc/1000G_EUR_Phase3_plink"
  }

  PREP_DEFINED_LOCI_INPUT(
    ch_pairs,
    clump_py,
    define_loci_py,
    ref_bed,
    ref_bim,
    ref_fam
  )
}