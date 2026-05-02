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

include { STAGE1_QC }           from './workflows/local/qc/main'
include { STAGE1_LDSC }         from './workflows/local/ldsc/main'
include { STAGE1_HDL }          from './workflows/local/hdl/main'
include { STAGE1_SUMHER }       from './workflows/local/sumher/main'
include { STAGE1_TWOSAMPLEMR }  from './workflows/local/twosamplemr/main'
include { STAGE1_MAGMA }        from './workflows/local/magma/main'
include { STAGE1_LAVA }         from './workflows/local/lava/main'
include { STAGE1_CONJFDR }      from './workflows/local/conjfdr/main'
include { STAGE1_DEF_LOCI }     from './workflows/local/def_loci/main'
include { STAGE1_GWAS_COLOC }   from './workflows/local/gwas_coloc/main'
include { STAGE1_SUSIE }        from './workflows/local/susie/main'
include { STAGE1_FUMA }         from './workflows/local/fuma/main'
include { STAGE1_TARGET_GENES } from './workflows/local/target_genes/main'
include { STAGE1_QTL_MANIFEST } from './workflows/local/qtl_manifest/main'
include { STAGE1_SMR }          from './workflows/local/smr/main'


workflow {

    log.info "---- Welcome to nf-core/graft!----"

    if (params.run_qc) {
        STAGE1_QC()
    }

    if (params.run_ldsc) {
        STAGE1_LDSC()
    }

    if (params.run_hdl) {
        STAGE1_HDL()
    }

    if (params.run_sumher) {
        STAGE1_SUMHER()
    }

    if (params.run_twosamplemr) {
        STAGE1_TWOSAMPLEMR()
    }

    if (params.run_magma) {
        STAGE1_MAGMA()
    }

    if (params.run_lava) {
        STAGE1_LAVA()
    }

    if (params.run_conjfdr) {
        STAGE1_CONJFDR()
    }

    if (params.run_def_loci) {
        STAGE1_DEF_LOCI()
    }

    if (params.run_gwas_coloc) {
        STAGE1_GWAS_COLOC()
    }

    if (params.run_susie) {
        STAGE1_SUSIE()
    }

    if (params.run_fuma) {
        STAGE1_FUMA()
    }

    if (params.run_target_genes) {
        STAGE1_TARGET_GENES()
    }

    if (params.run_qtl_manifest) {
        STAGE1_QTL_MANIFEST()
    }

    if (params.run_smr) {
        STAGE1_SMR()
    }
}
