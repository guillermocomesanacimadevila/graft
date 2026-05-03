#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
PULL STUFF FROM THE MANIFEST 
WHERE TYPE == BULK 
AND RUN IT FOR ALL WHERE TYPE == BULK
*/

process PREP_BULK_QTL_FOR_SMR {
    tag "${meta.id}_${meta.qtl_modality}_prep"

    input:
    tuple val(meta), path(qtl_file)
    path prep_script
    path gtf_file

    output:
    tuple val(meta), path("${meta.id}*"), emit: smr_qtl

    script:
    """
    set -euo pipefail

    python ${prep_script} \\
        --input ${qtl_file} \\
        --gtf ${gtf_file} \\
        --outdir . \\
        --qtl-type bulk \\
        --dataset-id ${meta.id} \\
        --cell-type bulk \\
        --gene-col ${meta.gene_col} \\
        --snp-col ${meta.snp_col} \\
        --chr-col ${meta.chr_col} \\
        --bp-col ${meta.bp_col} \\
        --a1-col ${meta.ea_col} \\
        --a2-col ${meta.oa_col} \\
        --freq-col ${meta.eaf_col} \\
        --beta-col ${meta.beta_col} \\
        --se-col ${meta.se_col} \\
        --p-col ${meta.pval_col} \\
        --default-freq ${params.default_qtl_freq} \\
        --input-build ${meta.build} \\
        --output-build hg19 \\
        --run-smr \\
        --smr-bin smr
    """
}

process PREP_GWAS_FOR_SMR {
    tag "${meta.id}_gwas_smr"

    input:
    tuple val(meta), path(gwas)
    path format_script

    output:
    tuple val(meta), path("${meta.id}.smr.ma"), emit: smr_gwas

    script:
    """
    set -euo pipefail

    python ${format_script} ${meta.id} \\
        --infile ${gwas} \\
        --out_dir . \\
        --snp_col ${meta.snp_col} \\
        --a1_col ${meta.a1_col} \\
        --a2_col ${meta.a2_col} \\
        --freq_col ${meta.eaf_col} \\
        --beta_col ${meta.beta_col} \\
        --se_col ${meta.se_col} \\
        --p_col ${meta.p_col} \\
        --n_col ${meta.n_col}
    """
}

/*

process PREP_GWAS_FOR_SMR {
    tag "${meta.id}_gwas_smr"

    input:
    tuple val(meta), path(gwas)
    path format_script

    output:
    tuple val(meta), path("${meta.id}.smr.ma"), emit: smr_gwas

    script:
    """
    set -euo pipefail

    python ${format_script} ${meta.id} \\
        --infile ${gwas} \\
        --out_dir . \\
        --snp_col SNP \\
        --a1_col A1 \\
        --a2_col A2 \\
        --freq_col FRQ \\
        --beta_col BETA \\
        --se_col SE \\
        --p_col P \\
        --n_col N \\
        --sep "	"
    """
}

*/

process BULK_SMR_HEIDI {
    tag "${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.qtl_modality}"

    input:
    tuple val(gwas_meta), path(gwas_smr), val(qtl_meta), path(qtl_smr_files)

    output:
    tuple val(gwas_meta), val(qtl_meta), path("${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.qtl_modality}.smr"), emit: smr_res

    script:
    """
    set -euo pipefail

    smr \\
        --bfile ${params.smr_bfile} \\
        --gwas-summary ${gwas_smr} \\
        --beqtl-summary ${qtl_meta.id} \\
        --maf ${params.maf_min} \\
        --peqtl-smr ${params.peqtl_smr} \\
        --peqtl-heidi ${params.peqtl_heidi} \\
        --thread-num ${params.thread_num} \\
        --out ${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.qtl_modality}
    """
}

process GRAB_BULK_SMR_HITS {
    tag "${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.qtl_modality}_hits"

    input:
    tuple val(gwas_meta), val(qtl_meta), path(smr_file)
    path grab_script

    output:
    tuple val(gwas_meta), val(qtl_meta), path("SMR"), emit: hits

    script:
    """
    set -euo pipefail

    python ${grab_script} \\
        --qtl_dataset ${qtl_meta.id} \\
        --qtl_type ${qtl_meta.qtl_modality} \\
        --pheno_id ${gwas_meta.id} \\
        --smr_res ${smr_file} \\
        --heidi_thresh ${params.heidi_thresh} \\
        --fdr_thresh ${params.smr_fdr_thresh} \\
        --out_dir .
    """
}

/*
smr \
  --bfile /Users/c24102394/Desktop/neurobridge/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL \
  --gwas-summary /Users/c24102394/Desktop/neurobridge/data/AD_rep/SMR/AD_bellenguez_smr_ready.dedup.ma \
  --beqtl-summary /Users/c24102394/Desktop/neurobridge/ref/sc_eQTLs/fujita2/SMR_ready/Oli \
  --maf 0.01 \
  --peqtl-smr 5e-8 \
  --peqtl-heidi 1e-5 \
  --thread-num 8 \
  --out /Users/c24102394/Desktop/neurobridge/data/AD_rep/SMR/res/AD_bellenguez_OLI_dedup
*/