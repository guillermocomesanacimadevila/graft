#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/graft
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started on October 2025
    Escott-Price Lab; UK Dementia Research Institute
    Dev: Guillermo Comesaña Cimadevila
    GitHub: https://github.com/guillermocomesanacimadevila/graft

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { STAGE1_QC }          from './workflows/local/qc/main'
include { STAGE1_LDSC }        from './workflows/local/ldsc/main'
include { STAGE1_HDL }         from './workflows/local/hdl/main'
include { STAGE1_SUMHER }      from './workflows/local/sumher/main'
include { STAGE1_TWOSAMPLEMR } from './workflows/local/twosamplemr/main'
include { STAGE1_MAGMA }       from './workflows/local/magma/main'
include { STAGE1_LAVA }        from './workflows/local/lava/main'
include { STAGE1_CONJFDR }     from './workflows/local/conjfdr/main'
include { STAGE1_DEF_LOCI }    from './workflows/local/def_loci/main'
include { STAGE1_GWAS_COLOC }  from './workflows/local/gwas_coloc/main'
include { STAGE1_SUSIE }       from './workflows/local/susie/main'
include { STAGE1_FUMA }        from './workflows/local/fuma/main'
// include { TWOSAMPLEMR } from
// include { STAGE1_FUMA } from './workflows/local/fuma/main'
// include { STAGE1_SMR_HEIDI } from './workflows/local/fuma/main'
// include { STAGE1_POST_SMR } from './workflows/local/fuma/main'
// include { STAGE1_QTL_COLOC } from './workflows/local/fuma/main'
// include { STAGE1_SMR_SUSIE_LD } from './workflows/local/fuma/main'
// MAYBE INCLUDE MR? (Between SumHer and LAVA?)

// GRAFT: GWAS Relatedness, Architecture & Functional Trait mapping
// cross-trait localised pipeline - bridges genetic and biological layers
workflow {
    STAGE1_QC()
    STAGE1_LDSC()
    STAGE1_HDL()
    STAGE1_SUMHER()
    STAGE1_TWOSAMPLEMR()
    STAGE1_MAGMA()
    STAGE1_LAVA()
    STAGE1_CONJFDR()
    STAGE1_DEF_LOCI()
    STAGE1_GWAS_COLOC()
    STAGE1_SUSIE()
    STAGE1_FUMA()   
    // SUSIE -
    // FUMA - 
    // MAGMA - 
    // Bulk - SMR + HEIDI
    // Bulk - Coloc
    // sc - SMR + HEIDI
    // sc - Coloc
}
