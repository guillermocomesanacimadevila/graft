#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/graft
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started on October 2025
    Escott-Price Lab; UK Dementia Research Institute
    Dev: Guillermo Comesaña Cimadevila
    GitHub: https://github.com/guillermocomesanacimadevila/neurobridge

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

include { STAGE1_QC }      from './workflows/local/qc/main'
include { STAGE1_ALIGN }   from './workflows/local/align/main'
include { STAGE1_LDSC }    from './workflows/local/ldsc/main'
include { STAGE1_HDL }     from './workflows/local/hdl/main'
include { STAGE1_SUMHER }  from './workflows/local/sumher/main'
// include { STAGE1_LAVA } from './workflows/local/lava/main'
include { STAGE1_CONJFDR } from './workflows/local/conjfdr/main'
// include { STAGE1_DEF_LOCI } from './workflows/def_loci/lava/main'
// include { STAGE1_COLOCALIZATION } from './workflows/local/coloc/main'
// include { STAGE1_SuSiE } from './workflows/local/susie/main'
// include { STAGE1_FUMA } from './workflows/local/fuma/main'


// GRAFT: GWAS Relatedness, Architecture & Functional Trait mapping
workflow {
    STAGE1_QC()
    STAGE1_ALIGN()
    //STAGE1_LDSC()
    // STAGE1_HDL()
    //STAGE1_SUMHER()
    // STAGE1_LAVA()
    // STAGE1_MIXER()
    // STAGE1_CONJFDR()
    // CLUMP - 
    // COLOC -
    // SUSIE -
    // FUMA - 
    // MAGMA - 
    // Bulk - SMR + HEIDI
    // Bulk - Coloc
    // sc - SMR + HEIDI
    // sc - Coloc
}
