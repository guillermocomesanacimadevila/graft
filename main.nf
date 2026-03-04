#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/neurobridge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started on October 2025
    Escott-Price Lab; UK Dementia Research Institute
    GitHub: https://github.com/guillermocomesanacimadevila/neurobridge

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

def GREEN = '\033[38;5;40m'   
def DARK = '\033[38;5;236m'  
def BOLD = '\033[1m'
def RESET = '\033[0m'

println """
${GREEN}${BOLD}
 ███╗   ██╗███████╗       ██████╗ ██████╗ ██████╗ ███████╗
 ████╗  ██║██╔════╝      ██╔════╝██╔═══██╗██╔══██╗██╔════╝
 ██╔██╗ ██║█████╗  █████╗██║     ██║   ██║██████╔╝█████╗
 ██║╚██╗██║██╔══╝  ╚════╝██║     ██║   ██║██╔══██╗██╔══╝
 ██║ ╚████║██║            ╚██████╗╚██████╔╝██║  ██║███████╗
 ╚═╝  ╚═══╝╚═╝             ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝
${RESET}

${GREEN}${BOLD}
███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ██████╗ ██████╗ ██╗██████╗  ██████╗ ███████╗
████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔══██╗██╔══██╗██║██╔══██╗██╔════╝ ██╔════╝
██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██████╔╝██████╔╝██║██║  ██║██║  ███╗█████╗
██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██╔══██╗██╔══██╗██║██║  ██║██║   ██║██╔══╝
██║ ╚████║███████╗╚██████╔╝██   ██╔╝╚█████╔╝███████║██║  ██║██║██████╔╝╚██████╔╝███████╗
╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═════╝  ╚═════╝ ╚══════╝
${RESET}

    Escott-Price Lab | UKDRI Cardiff
    v${workflow.manifest.version ?: 'dev'}

------------------------------------------------------------
"""

include { QC_GWAS }        from './modules/qc_gwas'
include { ADD_NEFF }       from './modules/add_neff'
include { LDSC_PAIRWISE }  from './subworkflows/ldsc_pairwise'
include { HDL_L_PAIRS }    from './subworkflows/hdl_pairs'
include { SUMHER_RUN }     from './subworkflows/sumher_run'
include { LAVA_GWAS_PREP } from './modules/lava_prep'
include { LAVA_RUN }       from './subworkflows/lava_pairwise'

// nextflow run main.nf -profile local  -c conf/local/nextflow.config --input assets/gwas.tsv --pairs assets/ldsc_pairs.tsv --outdir results
// nextflow run main.nf -profile docker -c conf/local/nextflow.config --input assets/gwas.tsv --pairs assets/ldsc_pairs.tsv --outdir results

/*
pre-run check;
-> if workflow.launchDir != nerurobridge -> cd neurobridge/
*/


workflow {

  if( !params.input )
    error "Missing --input (samplesheet tsv)"

  if( !params.pairs )
    error "Missing --pairs (pairs tsv)"

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
      meta.cases              = (row.cases ?: "").toString().trim()
      meta.controls           = (row.controls ?: "").toString().trim()

      tuple(meta, file(row.gwas))
    }

  // scr
  qc_script = file("${workflow.launchDir}/bin/qc_gwas.py")
  neff_script = file("${workflow.launchDir}/bin/compute_neff.py")
  ldsc_r = file("${workflow.launchDir}/bin/ldsc.R")
  hdl_l_r = file("${workflow.launchDir}/bin/hdl_l.R")
  calc_p = file("${workflow.launchDir}/bin/calc_p.py")
  lava_data_prep = file("${workflow.launchDir}/bin/prep_data.py")
  lava_r = file("${workflow.launchDir}/bin/lava_pair.R")

  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~ REFERENCE FILES ~~~~~~~~ //
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 

  // refs (LDSC)
  hm3_snplist = file("${workflow.launchDir}/ref/ldsc/w_hm3.snplist")
  ld_chr_dir = file("${workflow.launchDir}/ref/ldsc/eur_w_ld_chr")
  wld_dir = file("${workflow.launchDir}/ref/ldsc/weights_hm3_no_hla")

  // refs (HDL-L)
  ld_path = file("${workflow.launchDir}/ref/HDL-L_ref/LD.path")
  bim_path = file("${workflow.launchDir}/ref/HDL-L_ref/bimfile")

  // ref (SumHer)
  plink_dir = file("${workflow.launchDir}/ref/ldsc/1000G_EUR_Phase3_plink")
  ldak_bin = file("${workflow.launchDir}/ref/SumHer/LDAK/${params.ldak_os}")

  // ref (LAVA)
  lava_ref_dir = file("${workflow.launchDir}/ref/lava/lava_ref")
  lava_ref_prefix = "lava-ukb-v1.1"
  loci_file = file("${workflow.launchDir}/ref/lava/hdll_blocks.coords.loci")

  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~  ~~~~~~~~~~~~~~ ~~~~~~~~ //
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 
  // ~~~~~~~~ ~~~~~~~~~~~~~~~ ~~~~~~~~ // 

  // QC + NEFF
  ch_qc = QC_GWAS(ch_in, qc_script).ldsc_ready
  ch_neff = ADD_NEFF(ch_qc, neff_script).ldsc_neff

  ch_sum = ch_neff.map { meta, f -> tuple(meta.id, f) }
  ch_sumstats = ch_neff.map { meta, f -> tuple(meta, f) }

  ch_pairs = Channel
    .fromPath(params.pairs)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
      def meta = [
        trait1: row.trait1.toString().trim(),
        trait2: row.trait2.toString().trim(),
        cases1: row.cases1,
        controls1: row.controls1,
        cases2: row.cases2,
        controls2: row.controls2,
        pop_prev1: row.pop_prev1,
        pop_prev2: row.pop_prev2
      ]
      tuple(meta.trait1, meta.trait2, meta)
    }

  // join -> (meta, sum1, sum2)
  ch_with_t1 = ch_pairs
    .join(ch_sum, by: 0)
    .map { trait1, trait2, meta, f1 -> tuple(trait2, meta, f1) }

  ch_ldsc_in = ch_with_t1
    .join(ch_sum, by: 0)
    .map { trait2, meta, f1, f2 -> tuple(meta, f1, f2) }

  LDSC_PAIRWISE(
    ch_ldsc_in,
    ldsc_r,
    hm3_snplist,
    ld_chr_dir,
    wld_dir
  )

  ch_neff_key = ch_neff.map { meta, f -> tuple(meta.id, f) }

  ch_pairs_hdl = Channel
    .fromPath(params.pairs)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
      def meta = [
        trait1: row.trait1.toString().trim(),
        trait2: row.trait2.toString().trim()
      ]
      tuple(meta.trait1, meta.trait2, meta)
    }

  ch_hdl_with_t1 = ch_pairs_hdl
    .join(ch_neff_key, by: 0)
    .map { trait1, trait2, meta, f1 ->
      def m = meta + [ beta1: "BETA", se1: "SE" ]
      tuple(trait2, m, f1)
    }

  ch_hdl_in = ch_hdl_with_t1
    .join(ch_neff_key, by: 0)
    .map { trait2, meta, f1, f2 ->
      def m = meta + [ beta2: "BETA", se2: "SE" ]
      tuple(m, f1, f2)
    }

  HDL_L_PAIRS(
    ch_hdl_in,
    hdl_l_r,
    ld_path,
    bim_path
  )

  /*
  SumHer (LDAK)
  */

  ch__sumher_pairs = Channel
    .fromPath(params.pairs)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
      def meta = [
        trait1: row.trait1.toString().trim(),
        trait2: row.trait2.toString().trim()
      ]
      tuple(meta.trait1, meta.trait2, meta)
    }

  SUMHER_RUN(
    ch_sumstats,
    ch__sumher_pairs,
    plink_dir,
    calc_p,
    ldak_bin
  )

  /*
  LAVA
  */

  // ch_qc
  // ch_neff
  // TO Do´s
  // -> ADAPT from lava_main.nf

  
}
