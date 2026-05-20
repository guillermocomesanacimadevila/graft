#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREP_MAGMA_INPUT } from '../../../subworkflows/local/magma/main'
include { MAGMA }            from '../../../modules/local/magma'

workflow STAGE1_MAGMA {

  main:
  neurobridge.ParamUtils.requireParam(params.pairs, "pairs")

  ch_pairs = Channel
    .fromPath(params.pairs, checkIfExists: true)
    .splitCsv(header: true, sep: "\t")

  magma_bash      = file("${projectDir}/bin/magma_genome_wide.sh")
  g100eur_bed     = file("${projectDir}/ref/magma/g1000_eur.bed")
  g100eur_bim     = file("${projectDir}/ref/magma/g1000_eur.bim")
  g100eur_fam     = file("${projectDir}/ref/magma/g1000_eur.fam")
  g100eur_snploc  = file("${projectDir}/ref/magma/g1000_eur.snp.loc")
  grch37_gene_loc = file("${projectDir}/ref/magma/NCBI37.3.gene.loc")
  magma_bin       = file("${projectDir}/ref/magma/software/linux/magma")

  if( !magma_bash.exists() ) {
    error "Missing MAGMA wrapper script: ${magma_bash}"
  }

  if( !g100eur_bed.exists() || !g100eur_bim.exists() || !g100eur_fam.exists() ) {
    error "Missing MAGMA 1000G EUR reference BED/BIM/FAM files under ${projectDir}/ref/magma"
  }

  if( !g100eur_snploc.exists() ) {
    error "Missing MAGMA SNP location file: ${g100eur_snploc}"
  }

  if( !grch37_gene_loc.exists() ) {
    error "Missing MAGMA gene location file: ${grch37_gene_loc}"
  }

  if( !magma_bin.exists() ) {
    error "Missing MAGMA binary: ${magma_bin}"
  }

  ch_magma_in = PREP_MAGMA_INPUT(ch_pairs).magma_input

  MAGMA(
    ch_magma_in,
    magma_bash,
    g100eur_bed,
    g100eur_bim,
    g100eur_fam,
    g100eur_snploc,
    grch37_gene_loc,
    magma_bin
  )
}