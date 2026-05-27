#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(coloc) })

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript coloc_prior_files.R <prefix> <file1> <file2> <out_dir> <label1> <label2> <type1> <type2> [s1] [s2]")
}

prefix  <- args[1]
file1   <- args[2]
file2   <- args[3]
out_dir <- args[4]
label1  <- args[5]
label2  <- args[6]
type1   <- args[7]
type2   <- args[8]

s1 <- if (length(args) >= 9) as.numeric(args[9]) else 0.5
s2 <- if (length(args) >= 10) as.numeric(args[10]) else 0.5

p1 <- 1e-4
p2 <- 1e-4
p12_grid <- c(5e-6, 1e-5, 5e-5, 1e-4)

if (!file.exists(file1)) stop("File not found: ", file1)
if (!file.exists(file2)) stop("File not found: ", file2)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
out_file <- file.path(out_dir, paste0(prefix, "_coloc_p12_sensitivity.tsv"))

getN <- function(mm, suff) {
  cand <- c(paste0("N", suff), paste0("N_eff", suff), paste0("Neff", suff), paste0("NEFF", suff))
  cand <- unique(cand)
  for (nm in cand) if (nm %in% names(mm)) return(mm[[nm]])
  rep(NA_real_, nrow(mm))
}

empty_result <- function(n_snps) {
  do.call(rbind, lapply(p12_grid, function(p12) {
    data.frame(
      nsnps=NA,
      PP.H0.abf=NA,
      PP.H1.abf=NA,
      PP.H2.abf=NA,
      PP.H3.abf=NA,
      PP.H4.abf=NA,
      p1=p1,
      p2=p2,
      p12=p12,
      file1=basename(file1),
      file2=basename(file2),
      label1=label1,
      label2=label2,
      lead_snp=NA,
      n_snps=n_snps,
      stringsAsFactors=FALSE
    )
  }))
}

d1 <- read.table(file1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
d2 <- read.table(file2, sep="\t", header=TRUE, stringsAsFactors=FALSE)

m <- merge(d1, d2, by="SNP", suffixes=c(".1", ".2"))

if (nrow(m) < 10) {
  final <- empty_result(nrow(m))
  write.table(final, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote:", out_file, "\n")
  quit(save="no", status=0)
}

same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2

m$BETA.2[flip] <- -m$BETA.2[flip]
m <- m[same | flip, ]

if (nrow(m) < 10) {
  final <- empty_result(nrow(m))
  write.table(final, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Wrote:", out_file, "\n")
  quit(save="no", status=0)
}

lead_snp <- m$SNP[which.min(m$P.1)]

N1 <- unique(na.omit(getN(m, ".1")))
N2 <- unique(na.omit(getN(m, ".2")))
N1 <- if (length(N1)) as.numeric(N1[1]) else NA_real_
N2 <- if (length(N2)) as.numeric(N2[1]) else NA_real_

ds1 <- list(snp=m$SNP, beta=m$BETA.1, varbeta=m$SE.1^2, N=N1, type=type1)
ds2 <- list(snp=m$SNP, beta=m$BETA.2, varbeta=m$SE.2^2, N=N2, type=type2)

if (type1 == "cc") ds1$s <- s1
if (type2 == "cc") ds2$s <- s2
if (type1 == "quant") ds1$sdY <- 1
if (type2 == "quant") ds2$sdY <- 1

results <- vector("list", length(p12_grid))

for (j in seq_along(p12_grid)) {
  p12_val <- p12_grid[j]
  
  co <- coloc.abf(ds1, ds2, p1=p1, p2=p2, p12=p12_val)
  s <- as.data.frame(t(co$summary), stringsAsFactors=FALSE)
  
  s$p1 <- p1
  s$p2 <- p2
  s$p12 <- p12_val
  s$file1 <- basename(file1)
  s$file2 <- basename(file2)
  s$label1 <- label1
  s$label2 <- label2
  s$lead_snp <- lead_snp
  s$n_snps <- nrow(m)
  
  results[[j]] <- s
}

final <- do.call(rbind, results)
final <- final[order(final$p12), ]

write.table(final, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote:", out_file, "\n")