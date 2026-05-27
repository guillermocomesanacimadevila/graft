#!/usr/bin/env Rscript
# Run global HDL rg for a trait pair
# Writes outputs directly into the out_dir provided by the user
# Does NOT create a nested trait-pair directory

suppressPackageStartupMessages({
  library(data.table)
  library(HDL)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 14) {
  stop(paste(
    "Usage:",
    "hdl_global.R <gwas1_ldsc.tsv> <gwas2_ldsc.tsv> <trait1> <trait2> <beta1_col> <se1_col> <beta2_col> <se2_col> <LD.path> <out_dir> <eigen_cut> <N0> <lim> <jackknife>",
    sep = "\n"
  ))
}

gwas1_in  <- args[1]
gwas2_in  <- args[2]
trait1    <- args[3]
trait2    <- args[4]
beta1_col <- args[5]
se1_col   <- args[6]
beta2_col <- args[7]
se2_col   <- args[8]
LD.path   <- args[9]
out_dir   <- args[10]
eigen_cut <- as.numeric(args[11])
N0        <- as.numeric(args[12])
lim       <- as.numeric(args[13])
jackknife <- as.logical(args[14])

pair <- paste0(trait1, "_", trait2)

stopifnot(file.exists(gwas1_in))
stopifnot(file.exists(gwas2_in))
stopifnot(dir.exists(out_dir))
stopifnot(dir.exists(LD.path))

compute_z <- function(df, beta_col, se_col) {
  df <- as.data.table(df)
  df[, Z := get(beta_col) / get(se_col)]
  df
}

to_hdl <- function(df, trait, beta_col, se_col) {
  df <- as.data.table(df)

  need <- c("SNP", "A1", "A2", beta_col, se_col, "N")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) {
    stop(paste(trait, "missing columns:", paste(miss, collapse = ", ")))
  }

  df <- df[!is.na(SNP) & !is.na(A1) & !is.na(A2)]
  df <- df[!is.na(get(beta_col)) & !is.na(get(se_col)) & !is.na(N)]

  df[, SNP := as.character(SNP)]
  df[, A1 := toupper(as.character(A1))]
  df[, A2 := toupper(as.character(A2))]
  df[, N := as.numeric(N)]
  df[, (beta_col) := as.numeric(get(beta_col))]
  df[, (se_col) := as.numeric(get(se_col))]

  df <- df[
    is.finite(N) &
    is.finite(get(beta_col)) &
    is.finite(get(se_col)) &
    get(se_col) > 0
  ]

  df <- compute_z(df, beta_col, se_col)
  df <- df[is.finite(Z)]
  df <- df[!duplicated(SNP)]

  df[, .(SNP, A1, A2, N, Z)]
}

gwas1_raw <- fread(gwas1_in, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
gwas2_raw <- fread(gwas2_in, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)

gwas1 <- to_hdl(gwas1_raw, trait1, beta1_col, se1_col)
gwas2 <- to_hdl(gwas2_raw, trait2, beta2_col, se2_col)

fwrite(gwas1, file.path(out_dir, paste0(trait1, ".hdl.tsv")), sep = "\t")
fwrite(gwas2, file.path(out_dir, paste0(trait2, ".hdl.tsv")), sep = "\t")

global_fit <- HDL.rg(
  gwas1.df = gwas1,
  gwas2.df = gwas2,
  LD.path = LD.path,
  N0 = N0,
  eigen.cut = eigen_cut,
  lim = lim,
  jackknife = jackknife
)

rg_se <- if (!is.null(global_fit$rg.se)) global_fit$rg.se else NA_real_
rg_z  <- if (is.finite(global_fit$rg) && is.finite(rg_se) && rg_se > 0) global_fit$rg / rg_se else NA_real_
rg_p  <- if (!is.null(global_fit$P)) global_fit$P else if (is.finite(rg_z)) 2 * pnorm(-abs(rg_z)) else NA_real_

global_tab <- data.table(
  trait1 = trait1,
  trait2 = trait2,
  rg = global_fit$rg,
  SE = rg_se,
  z = rg_z,
  p = rg_p,
  eigen_cut = eigen_cut,
  N0 = N0,
  lim = lim,
  jackknife = jackknife
)

fwrite(
  global_tab,
  file.path(out_dir, paste0(pair, ".global_rg.tsv")),
  sep = "\t"
)
saveRDS(
  global_tab,
  file.path(out_dir, paste0(pair, ".global_rg.rds"))
)

if (!is.null(global_fit$estimates.df)) {
  est <- as.data.frame(global_fit$estimates.df)
  est$parameter <- rownames(est)
  est <- est[, c("parameter", setdiff(colnames(est), "parameter")), drop = FALSE]
  est <- as.data.table(est)

  fwrite(
    est,
    file.path(out_dir, paste0(pair, ".global_estimates.tsv")),
    sep = "\t"
  )
  saveRDS(
    est,
    file.path(out_dir, paste0(pair, ".global_estimates.rds"))
  )
}

if (!is.null(global_fit$jackknife.df)) {
  jk <- as.data.table(global_fit$jackknife.df)
  fwrite(
    jk,
    file.path(out_dir, paste0(pair, ".jackknife.tsv")),
    sep = "\t"
  )
  saveRDS(
    jk,
    file.path(out_dir, paste0(pair, ".jackknife.rds"))
  )
}

saveRDS(
  global_fit,
  file.path(out_dir, paste0(pair, ".global_fit.rds"))
)

print(global_tab)