#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(remotes)
  library(data.table)
  library(mvsusieR)
})

# install.packages("remotes")
# remotes::install_github("stephenslab/mvsusieR")

ad_file <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/gwas_AD.ldorder.tsv"
scz_file <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/gwas_SCZ.ldorder.tsv"
ld_file <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/ld_gwas/locus.ld.gz"
ld_snps_file <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/ld_gwas/snps_in_ld_order.txt"

out_rds <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/mvsusie_ad_scz.rds"
out_pip <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/mvsusie_ad_scz.pip.tsv"
out_cs <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/mvsusie_ad_scz.cs.tsv"
out_alpha <- "/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174/mvsusie_ad_scz.alpha.tsv"

ad <- fread(ad_file)
scz <- fread(scz_file)
ld_snps <- fread(ld_snps_file, header = FALSE)
setnames(ld_snps, "V1", "SNP")

ad <- ad[, .(SNP, A1_AD = A1, A2_AD = A2, N_AD = N, BETA_AD = BETA, SE_AD = SE, P_AD = P, CHR_AD = CHR, BP_AD = BP)]
scz <- scz[, .(SNP, A1_SCZ = A1, A2_SCZ = A2, N_SCZ = N, BETA_SCZ = BETA, SE_SCZ = SE, P_SCZ = P, CHR_SCZ = CHR, BP_SCZ = BP)]

dt <- merge(ld_snps, ad, by = "SNP", all.x = TRUE, sort = FALSE)
dt <- merge(dt, scz, by = "SNP", all.x = TRUE, sort = FALSE)

if (any(is.na(dt$BETA_AD)) || any(is.na(dt$BETA_SCZ))) stop("Missing SNPs after merging with LD SNP order")

same <- dt$A1_AD == dt$A1_SCZ & dt$A2_AD == dt$A2_SCZ
flip <- dt$A1_AD == dt$A2_SCZ & dt$A2_AD == dt$A1_SCZ

if (any(!(same | flip))) stop("Found SNPs with alleles that do not match between traits")

dt[flip, `:=`(BETA_SCZ = -BETA_SCZ, A1_SCZ = A1_AD, A2_SCZ = A2_AD)]

dt[, Z_AD := BETA_AD / SE_AD]
dt[, Z_SCZ := BETA_SCZ / SE_SCZ]

Z <- as.matrix(dt[, .(Z_AD, Z_SCZ)])
R <- as.matrix(fread(ld_file, header = FALSE))

if (nrow(Z) != nrow(R)) stop("Z matrix row count does not match LD matrix dimension")

N_use <- floor(min(median(dt$N_AD, na.rm = TRUE), median(dt$N_SCZ, na.rm = TRUE)))

prior <- create_mixture_prior(R = 2)

fit <- mvsusie_rss(
  Z = Z,
  R = R,
  N = N_use,
  L = 10,
  prior_variance = prior
)

saveRDS(fit, out_rds)

pip_dt <- data.table(
  SNP = dt$SNP,
  CHR = dt$CHR_AD,
  BP = dt$BP_AD,
  A1 = dt$A1_AD,
  A2 = dt$A2_AD,
  PIP = fit$pip
)

fwrite(pip_dt, out_pip, sep = "\t")

alpha_dt <- as.data.table(t(fit$alpha))
setnames(alpha_dt, paste0("alpha_L", seq_len(ncol(alpha_dt))))
alpha_dt <- cbind(data.table(SNP = dt$SNP), alpha_dt)
fwrite(alpha_dt, out_alpha, sep = "\t")

if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
  cs_dt <- rbindlist(lapply(seq_along(fit$sets$cs), function(i) {
    idx <- fit$sets$cs[[i]]
    data.table(
      CS = paste0("CS", i),
      SNP = dt$SNP[idx],
      CHR = dt$CHR_AD[idx],
      BP = dt$BP_AD[idx],
      A1 = dt$A1_AD[idx],
      A2 = dt$A2_AD[idx],
      PIP = fit$pip[idx],
      alpha = apply(fit$alpha[, idx, drop = FALSE], 2, max)
    )
  }))
} else {
  cs_dt <- data.table(
    CS = character(),
    SNP = character(),
    CHR = integer(),
    BP = integer(),
    A1 = character(),
    A2 = character(),
    PIP = numeric(),
    alpha = numeric()
  )
}

fwrite(cs_dt, out_cs, sep = "\t")

print(fit$sets)
print(pip_dt[order(-PIP)][1:20, .(SNP, CHR, BP, A1, A2, PIP)])