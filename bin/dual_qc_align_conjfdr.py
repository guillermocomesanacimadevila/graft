#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl

bases = ["A", "C", "G", "T"]

def count(label: str, df: pl.DataFrame) -> None:
    print(f"{label}: {df.height:,}")

def read_table(path: Path, sep: str, has_header: bool = True) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator=sep,
        has_header=has_header,
        infer_schema_length=5000,
        ignore_errors=True,
        null_values=["NA", "NaN", "nan", "", "null", "NULL"],
    )

def as_upper(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(
        pl.col(col).cast(pl.Utf8, strict=False).str.strip_chars().str.to_uppercase().alias(col)
    )

def to_int(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(pl.col(col).cast(pl.Int64, strict=False).alias(col))

def to_float(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(pl.col(col).cast(pl.Float64, strict=False).alias(col))

def exclude_region(df: pl.DataFrame, chr_col: str, pos_col: str, chrom: int, start: int, end: int) -> pl.DataFrame:
    if chr_col not in df.columns or pos_col not in df.columns:
        return df
    bad = (pl.col(chr_col) == chrom) & (pl.col(pos_col) >= start) & (pl.col(pos_col) <= end)
    return df.filter(~bad)

def keep_snps_only(df: pl.DataFrame, a1_col: str, a2_col: str) -> pl.DataFrame:
    if a1_col not in df.columns or a2_col not in df.columns:
        return df
    a1 = pl.col(a1_col)
    a2 = pl.col(a2_col)
    ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
    ok_bases = a1.is_in(bases) & a2.is_in(bases)
    no_gap = ~a1.str.contains("-") & ~a2.str.contains("-")
    return df.filter(ok_len & ok_bases & no_gap)

def drop_palindromes(df: pl.DataFrame, a1_col: str, a2_col: str) -> pl.DataFrame:
    if a1_col not in df.columns or a2_col not in df.columns:
        return df
    a1 = pl.col(a1_col)
    a2 = pl.col(a2_col)
    pal = (
        ((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A")) |
        ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C"))
    )
    return df.filter(~pal)

def add_frq_maf(
    df: pl.DataFrame,
    frq_col: str,
    eaf_col: str | None,
    fcas_col: str | None,
    fcon_col: str | None,
    ncas_col: str | None,
    ncon_col: str | None,
) -> pl.DataFrame:
    if frq_col in df.columns:
        pass
    elif eaf_col and (eaf_col in df.columns):
        df = df.with_columns(pl.col(eaf_col).cast(pl.Float64, strict=False).alias(frq_col))
    elif fcas_col and fcon_col and (fcas_col in df.columns) and (fcon_col in df.columns):
        fcas = pl.col(fcas_col)
        fcon = pl.col(fcon_col)
        if ncas_col and ncon_col and (ncas_col in df.columns) and (ncon_col in df.columns):
            ncas = pl.col(ncas_col)
            ncon = pl.col(ncon_col)
            w = ncas + ncon
            frq = pl.when(w > 0).then((fcas * ncas + fcon * ncon) / w).otherwise((fcas + fcon) / 2.0)
        else:
            frq = (fcas + fcon) / 2.0
        df = df.with_columns(frq.alias(frq_col))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias(frq_col))
    df = df.with_columns(
        pl.min_horizontal(pl.col(frq_col), (1.0 - pl.col(frq_col))).alias("MAF")
    )
    return df

def add_n(
    df: pl.DataFrame,
    n_out_col: str,
    n_col: str | None,
    ncas_col: str | None,
    ncon_col: str | None,
) -> pl.DataFrame:
    if n_out_col in df.columns:
        return df
    if n_col and (n_col in df.columns):
        return df.with_columns(pl.col(n_col).cast(pl.Float64, strict=False).alias(n_out_col))
    if ncas_col and ncon_col and (ncas_col in df.columns) and (ncon_col in df.columns):
        return df.with_columns(
            (pl.col(ncas_col).cast(pl.Float64, strict=False) + pl.col(ncon_col).cast(pl.Float64, strict=False)).alias(n_out_col)
        )
    return df.with_columns(pl.lit(None, dtype=pl.Float64).alias(n_out_col))

def add_info(df: pl.DataFrame, info_out_col: str, info_col: str | None, require_info: bool) -> pl.DataFrame:
    if info_out_col in df.columns:
        return df
    if info_col and (info_col in df.columns):
        return df.with_columns(pl.col(info_col).cast(pl.Float64, strict=False).alias(info_out_col))
    if require_info:
        raise SystemExit("ERROR: INFO required but --info_col not provided or not found")
    return df.with_columns(pl.lit(None, dtype=pl.Float64).alias(info_out_col))

def filter_maf_info(df: pl.DataFrame, maf_min: float, info_min: float, info_col: str, require_info: bool) -> pl.DataFrame:
    if "MAF" in df.columns:
        df = df.filter(pl.col("MAF") >= maf_min)
    if require_info:
        df = df.filter(pl.col(info_col).is_not_null() & (pl.col(info_col) >= info_min))
    else:
        if info_col in df.columns:
            df = df.filter(pl.col(info_col).is_null() | (pl.col(info_col) >= info_min))
    return df

def drop_missing_required(df: pl.DataFrame, required_cols: list[str]) -> pl.DataFrame:
    have = [c for c in required_cols if c in df.columns]
    return df.drop_nulls(subset=have)

def comp_expr(colname: str) -> pl.Expr:
    return (
        pl.when(pl.col(colname) == "A").then(pl.lit("T"))
        .when(pl.col(colname) == "T").then(pl.lit("A"))
        .when(pl.col(colname) == "C").then(pl.lit("G"))
        .when(pl.col(colname) == "G").then(pl.lit("C"))
        .otherwise(pl.lit(None))
    )

def qc_gwas(
    df: pl.DataFrame,
    label: str,
    snp_col: str,
    chr_col: str,
    pos_col: str,
    a1_col: str,
    a2_col: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    eaf_col: str | None,
    freq_case_col: str | None,
    freq_ctrl_col: str | None,
    n_case_col: str | None,
    n_ctrl_col: str | None,
    n_col: str | None,
    info_col: str | None,
    require_info: bool,
    maf_min: float,
    info_min: float,
    exclude_mhc_flag: bool,
    exclude_apoe_flag: bool,
    apoe_chr: int,
    apoe_start: int,
    apoe_end: int,
    drop_palindromes_flag: bool,
    keep_snps_only_flag: bool,
) -> tuple[pl.DataFrame, str, str, str]:
    count(f"{label}_loaded", df)
    df = as_upper(df, a1_col)
    df = as_upper(df, a2_col)
    df = to_int(df, chr_col)
    df = to_int(df, pos_col)
    df = to_float(df, beta_col)
    df = to_float(df, se_col)
    df = to_float(df, p_col)
    if n_col:
        df = to_float(df, n_col)
    if info_col:
        df = to_float(df, info_col)
    if eaf_col:
        df = to_float(df, eaf_col)
    if freq_case_col:
        df = to_float(df, freq_case_col)
    if freq_ctrl_col:
        df = to_float(df, freq_ctrl_col)
    if n_case_col:
        df = to_float(df, n_case_col)
    if n_ctrl_col:
        df = to_float(df, n_ctrl_col)
    count(f"{label}_after_types", df)
    if exclude_mhc_flag:
        df = exclude_region(df, chr_col, pos_col, 6, 25_000_000, 34_000_000)
        count(f"{label}_after_exclude_MHC", df)
    if exclude_apoe_flag:
        df = exclude_region(df, chr_col, pos_col, apoe_chr, apoe_start, apoe_end)
        count(f"{label}_after_exclude_APOE", df)
    if keep_snps_only_flag:
        df = keep_snps_only(df, a1_col, a2_col)
        count(f"{label}_after_remove_indels", df)
    if drop_palindromes_flag:
        df = drop_palindromes(df, a1_col, a2_col)
        count(f"{label}_after_remove_palindromes", df)
    frq_col = f"__frq_{label}__"
    n_out_col = f"__n_{label}__"
    info_out_col = f"__info_{label}__"
    df = add_frq_maf(
        df,
        frq_col=frq_col,
        eaf_col=eaf_col,
        fcas_col=freq_case_col,
        fcon_col=freq_ctrl_col,
        ncas_col=n_case_col,
        ncon_col=n_ctrl_col,
    )
    df = add_n(df, n_out_col=n_out_col, n_col=n_col, ncas_col=n_case_col, ncon_col=n_ctrl_col)
    df = add_info(df, info_out_col=info_out_col, info_col=info_col, require_info=require_info)
    df = filter_maf_info(df, maf_min=maf_min, info_min=info_min, info_col=info_out_col, require_info=require_info)
    count(f"{label}_after_maf_info", df)
    required = [
        snp_col, a1_col, a2_col, frq_col,
        beta_col, se_col, p_col, chr_col, pos_col
    ]
    df = drop_missing_required(df, required)
    count(f"{label}_after_dropna_required", df)
    df = df.unique(subset=[snp_col], keep="first")
    count(f"{label}_after_dedup_snp", df)
    return df, frq_col, n_out_col, info_out_col

def harmonise_two_gwas(
    df1: pl.DataFrame,
    df2: pl.DataFrame,
    snp1: str,
    a11: str,
    a21: str,
    beta1: str,
    snp2: str,
    a12: str,
    a22: str,
    beta2: str,
) -> pl.DataFrame:
    x = df1.rename({
        snp1: "SNP",
        a11: "A1_1",
        a21: "A2_1",
        beta1: "BETA_1",
    })
    y = df2.rename({
        snp2: "SNP",
        a12: "A1_2",
        a22: "A2_2",
        beta2: "BETA_2",
    })

    merged = x.join(y, on="SNP", how="inner")
    count("after_intersection", merged)

    merged = merged.with_columns(
        comp_expr("A1_2").alias("__A1_2C"),
        comp_expr("A2_2").alias("__A2_2C"),
    )

    same = (pl.col("A1_1") == pl.col("A1_2")) & (pl.col("A2_1") == pl.col("A2_2"))
    swap = (pl.col("A1_1") == pl.col("A2_2")) & (pl.col("A2_1") == pl.col("A1_2"))
    same_comp = (pl.col("A1_1") == pl.col("__A1_2C")) & (pl.col("A2_1") == pl.col("__A2_2C"))
    swap_comp = (pl.col("A1_1") == pl.col("__A2_2C")) & (pl.col("A2_1") == pl.col("__A1_2C"))

    merged = merged.with_columns(
        pl.when(same | same_comp).then(pl.lit("same"))
        .when(swap | swap_comp).then(pl.lit("swap"))
        .otherwise(pl.lit(None))
        .alias("__align_status")
    )

    merged = merged.filter(pl.col("__align_status").is_not_null())
    count("after_alignment_filter", merged)

    merged = merged.with_columns(
        pl.when(pl.col("__align_status") == "swap")
        .then(-pl.col("BETA_2"))
        .otherwise(pl.col("BETA_2"))
        .alias("BETA_2"),
        pl.col("A1_1").alias("A1"),
        pl.col("A2_1").alias("A2"),
    )

    merged = merged.drop(["__A1_2C", "__A2_2C", "__align_status"])
    return merged

def write_aligned_output(
    merged: pl.DataFrame,
    out: Path,
    source: str,
    snp_col: str,
    a1_col: str,
    a2_col: str,
    frq_col: str,
    n_col_out: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    chr_col: str,
    pos_col: str,
    info_col_out: str,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)

    src_suffix = "_1" if source == "1" else "_2"

    def col_or_null(name: str) -> pl.Expr:
        if name in merged.columns:
            return pl.col(name)
        return pl.lit(None)

    beta_name = "BETA_1" if source == "1" else "BETA_2"

    ldsc = merged.select(
        col_or_null("SNP").alias(snp_col),
        col_or_null("A1").alias(a1_col),
        col_or_null("A2").alias(a2_col),
        col_or_null(frq_col + src_suffix).alias("FRQ") if (frq_col + src_suffix) in merged.columns else col_or_null(frq_col).alias("FRQ"),
        col_or_null(n_col_out + src_suffix).alias("N") if (n_col_out + src_suffix) in merged.columns else col_or_null(n_col_out).alias("N"),
        col_or_null(beta_name).alias(beta_col),
        col_or_null(se_col + src_suffix).alias("SE") if (se_col + src_suffix) in merged.columns else col_or_null(se_col).alias("SE"),
        col_or_null(p_col + src_suffix).alias("P") if (p_col + src_suffix) in merged.columns else col_or_null(p_col).alias("P"),
        col_or_null(chr_col + src_suffix).alias("CHR") if (chr_col + src_suffix) in merged.columns else col_or_null(chr_col).alias("CHR"),
        col_or_null(pos_col + src_suffix).alias("POS") if (pos_col + src_suffix) in merged.columns else col_or_null(pos_col).alias("POS"),
        col_or_null(info_col_out + src_suffix).alias("INFO") if (info_col_out + src_suffix) in merged.columns else col_or_null(info_col_out).alias("INFO"),
    )
    ldsc.write_csv(out, separator="\t")

def main():
    ap = argparse.ArgumentParser(description="Master dual-GWAS QC + allele harmonisation for conjFDR (polars)")
    ap.add_argument("--in1", required=True)
    ap.add_argument("--in2", required=True)
    ap.add_argument("--pheno1", required=True)
    ap.add_argument("--pheno2", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sep1", default="\t")
    ap.add_argument("--sep2", default="\t")

    ap.add_argument("--snp_col1", required=True)
    ap.add_argument("--chr_col1", required=True)
    ap.add_argument("--pos_col1", required=True)
    ap.add_argument("--a1_col1", required=True)
    ap.add_argument("--a2_col1", required=True)
    ap.add_argument("--beta_col1", required=True)
    ap.add_argument("--se_col1", required=True)
    ap.add_argument("--p_col1", required=True)
    ap.add_argument("--eaf_col1", default=None)
    ap.add_argument("--freq_case_col1", default=None)
    ap.add_argument("--freq_ctrl_col1", default=None)
    ap.add_argument("--n_case_col1", default=None)
    ap.add_argument("--n_ctrl_col1", default=None)
    ap.add_argument("--n_col1", default=None)
    ap.add_argument("--info_col1", default=None)
    ap.add_argument("--require_info1", action="store_true")

    ap.add_argument("--snp_col2", required=True)
    ap.add_argument("--chr_col2", required=True)
    ap.add_argument("--pos_col2", required=True)
    ap.add_argument("--a1_col2", required=True)
    ap.add_argument("--a2_col2", required=True)
    ap.add_argument("--beta_col2", required=True)
    ap.add_argument("--se_col2", required=True)
    ap.add_argument("--p_col2", required=True)
    ap.add_argument("--eaf_col2", default=None)
    ap.add_argument("--freq_case_col2", default=None)
    ap.add_argument("--freq_ctrl_col2", default=None)
    ap.add_argument("--n_case_col2", default=None)
    ap.add_argument("--n_ctrl_col2", default=None)
    ap.add_argument("--n_col2", default=None)
    ap.add_argument("--info_col2", default=None)
    ap.add_argument("--require_info2", action="store_true")

    ap.add_argument("--maf_min", type=float, default=0.01)
    ap.add_argument("--info_min", type=float, default=0.90)
    ap.add_argument("--exclude_mhc", action="store_true")
    ap.add_argument("--exclude_apoe", action="store_true")
    ap.add_argument("--apoe_chr", type=int, default=19)
    ap.add_argument("--apoe_start", type=int, default=44_000_000)
    ap.add_argument("--apoe_end", type=int, default=46_500_000)
    ap.add_argument("--drop_palindromes", action="store_true")
    ap.add_argument("--keep_snps_only", action="store_true")

    a = ap.parse_args()

    df1 = read_table(Path(a.in1), a.sep1)
    df2 = read_table(Path(a.in2), a.sep2)

    df1, frq1, n1, info1 = qc_gwas(
        df=df1,
        label=a.pheno1,
        snp_col=a.snp_col1,
        chr_col=a.chr_col1,
        pos_col=a.pos_col1,
        a1_col=a.a1_col1,
        a2_col=a.a2_col1,
        beta_col=a.beta_col1,
        se_col=a.se_col1,
        p_col=a.p_col1,
        eaf_col=a.eaf_col1,
        freq_case_col=a.freq_case_col1,
        freq_ctrl_col=a.freq_ctrl_col1,
        n_case_col=a.n_case_col1,
        n_ctrl_col=a.n_ctrl_col1,
        n_col=a.n_col1,
        info_col=a.info_col1,
        require_info=a.require_info1,
        maf_min=a.maf_min,
        info_min=a.info_min,
        exclude_mhc_flag=a.exclude_mhc,
        exclude_apoe_flag=a.exclude_apoe,
        apoe_chr=a.apoe_chr,
        apoe_start=a.apoe_start,
        apoe_end=a.apoe_end,
        drop_palindromes_flag=a.drop_palindromes,
        keep_snps_only_flag=a.keep_snps_only,
    )

    df2, frq2, n2, info2 = qc_gwas(
        df=df2,
        label=a.pheno2,
        snp_col=a.snp_col2,
        chr_col=a.chr_col2,
        pos_col=a.pos_col2,
        a1_col=a.a1_col2,
        a2_col=a.a2_col2,
        beta_col=a.beta_col2,
        se_col=a.se_col2,
        p_col=a.p_col2,
        eaf_col=a.eaf_col2,
        freq_case_col=a.freq_case_col2,
        freq_ctrl_col=a.freq_ctrl_col2,
        n_case_col=a.n_case_col2,
        n_ctrl_col=a.n_ctrl_col2,
        n_col=a.n_col2,
        info_col=a.info_col2,
        require_info=a.require_info2,
        maf_min=a.maf_min,
        info_min=a.info_min,
        exclude_mhc_flag=a.exclude_mhc,
        exclude_apoe_flag=a.exclude_apoe,
        apoe_chr=a.apoe_chr,
        apoe_start=a.apoe_start,
        apoe_end=a.apoe_end,
        drop_palindromes_flag=a.drop_palindromes,
        keep_snps_only_flag=a.keep_snps_only,
    )

    rename1 = {
        a.se_col1: a.se_col1 + "_1",
        a.p_col1: a.p_col1 + "_1",
        a.chr_col1: a.chr_col1 + "_1",
        a.pos_col1: a.pos_col1 + "_1",
        frq1: frq1 + "_1",
        n1: n1 + "_1",
        info1: info1 + "_1",
    }
    rename2 = {
        a.se_col2: a.se_col2 + "_2",
        a.p_col2: a.p_col2 + "_2",
        a.chr_col2: a.chr_col2 + "_2",
        a.pos_col2: a.pos_col2 + "_2",
        frq2: frq2 + "_2",
        n2: n2 + "_2",
        info2: info2 + "_2",
    }

    df1 = df1.rename({k: v for k, v in rename1.items() if k in df1.columns})
    df2 = df2.rename({k: v for k, v in rename2.items() if k in df2.columns})

    merged = harmonise_two_gwas(
        df1=df1,
        df2=df2,
        snp1=a.snp_col1,
        a11=a.a1_col1,
        a21=a.a2_col1,
        beta1=a.beta_col1,
        snp2=a.snp_col2,
        a12=a.a1_col2,
        a22=a.a2_col2,
        beta2=a.beta_col2,
    )

    out1 = Path(a.outdir) / a.pheno1 / "post-qc" / f"{a.pheno1}_aligned.tsv"
    out2 = Path(a.outdir) / a.pheno2 / "post-qc" / f"{a.pheno2}_aligned.tsv"

    write_aligned_output(
        merged=merged,
        out=out1,
        source="1",
        snp_col=a.snp_col1,
        a1_col=a.a1_col1,
        a2_col=a.a2_col1,
        frq_col=frq1,
        n_col_out=n1,
        beta_col=a.beta_col1,
        se_col=a.se_col1,
        p_col=a.p_col1,
        chr_col=a.chr_col1,
        pos_col=a.pos_col1,
        info_col_out=info1,
    )

    write_aligned_output(
        merged=merged,
        out=out2,
        source="2",
        snp_col=a.snp_col2,
        a1_col=a.a1_col2,
        a2_col=a.a2_col2,
        frq_col=frq2,
        n_col_out=n2,
        beta_col=a.beta_col2,
        se_col=a.se_col2,
        p_col=a.p_col2,
        chr_col=a.chr_col2,
        pos_col=a.pos_col2,
        info_col_out=info2,
    )

    print("wrote:", out1)
    print("wrote:", out2)

if __name__ == "__main__":
    main()
