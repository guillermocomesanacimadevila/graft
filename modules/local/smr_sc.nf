#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PREP_SC_QTL_FOR_SMR {
    tag "${meta.id}_${meta.cell_type}_${meta.qtl_modality}_prep"

    input:
    tuple val(meta), path(qtl_file)
    path prep_script
    path gtf_file

    output:
    tuple val(meta), path("${meta.id}_${meta.cell_type}*"), emit: smr_qtl

    script:
    """
    set -euo pipefail

    python ${prep_script} \\
        --input ${qtl_file} \\
        --gtf ${gtf_file} \\
        --outdir . \\
        --qtl-type sc \\
        --dataset-id ${meta.id} \\
        --cell-type ${meta.cell_type} \\
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

process SC_SMR_HEIDI {
    tag "${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.cell_type}_${qtl_meta.qtl_modality}"

    input:
    tuple val(gwas_meta), path(gwas_smr), val(qtl_meta), path(qtl_smr_files)

    output:
    tuple val(gwas_meta), val(qtl_meta), path("${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.cell_type}_${qtl_meta.qtl_modality}.smr"), emit: smr_res

    script:
    """
    set -euo pipefail

    smr \\
        --bfile ${params.smr_bfile} \\
        --gwas-summary ${gwas_smr} \\
        --beqtl-summary ${qtl_meta.id}_${qtl_meta.cell_type} \\
        --maf ${params.maf_min} \\
        --peqtl-smr ${params.peqtl_smr} \\
        --peqtl-heidi ${params.peqtl_heidi} \\
        --thread-num ${params.thread_num} \\
        --out ${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.cell_type}_${qtl_meta.qtl_modality}
    """
}

process GRAB_SC_SMR_HITS {
    tag "${gwas_meta.id}_${qtl_meta.id}_${qtl_meta.cell_type}_${qtl_meta.qtl_modality}_hits"

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
        --cell_type ${qtl_meta.cell_type} \\
        --pheno_id ${gwas_meta.id} \\
        --smr_res ${smr_file} \\
        --heidi_thresh ${params.heidi_thresh} \\
        --fdr_thresh ${params.smr_fdr_thresh} \\
        --out_dir .
    """
}


/*
parser = argparse.ArgumentParser()
    parser.add_argument("--qtl_dataset", required=True)
    parser.add_argument("--qtl_type", required=True)
    parser.add_argument("--cell_type", required=True)
    parser.add_argument("--pheno_id", required=True)
    parser.add_argument("--smr_res", required=True)
    parser.add_argument("--heidi_thresh", type=float, required=True)
    parser.add_argument("--fdr_thresh", type=float, required=True)
    parser.add_argument("--out_dir", required=True)
*/