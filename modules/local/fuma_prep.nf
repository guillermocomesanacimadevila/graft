#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FUMA_PREP {

    tag "${meta.pair}_fuma_prep"

    input:
    tuple val(meta), path(loci_dir)

    output:
    tuple val(meta), path("${meta.pair}"), emit: fuma_ready

    script:
    def traits_str = meta.traits.join(' ')
    """
    set -euo pipefail

    mkdir -p "${meta.pair}/loci" "${meta.pair}/meta"
    echo -e "pair\\ttrait\\tlocus_dir\\tfile_gz" > "${meta.pair}/meta/manifest.tsv"

    for t in ${traits_str}; do
        mkdir -p "${meta.pair}/loci/\$t"
        for d in ${loci_dir}/locus_*; do
            [ -d "\$d" ] || continue
            in_ldgwas="\$d/gwas_\${t}.ldgwas.tsv"
            in_ldorder="\$d/gwas_\${t}.ldorder.tsv"
            in_tsv="\$d/gwas_\${t}.tsv"

            if [ -f "\$in_ldgwas" ]; then
                in="\$in_ldgwas"
            elif [ -f "\$in_ldorder" ]; then
                in="\$in_ldorder"
            elif [ -f "\$in_tsv" ]; then
                in="\$in_tsv"
            else
                continue
            fi

            locus_dir_name="\$(basename "\$d")"
            out="${meta.pair}/loci/\$t/\${locus_dir_name}_gwas_\${t}.tsv.gz"

            python3 - "\$in" "\$out" <<'PY'
import sys
import gzip
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t")

for col in ("N", "CHR", "POS", "BP"):
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        if col in ("N", "CHR", "POS", "BP"):
            df[col] = df[col].astype("Int64")

if "POS" not in df.columns and "BP" in df.columns:
    df["POS"] = df["BP"]

with gzip.open(outfile, "wt") as out:
    df.to_csv(out, sep="\t", index=False)
PY

            echo -e "${meta.pair}\\t\$t\\t\${locus_dir_name}\\t\$out" >> "${meta.pair}/meta/manifest.tsv"
        done
    done

    [ -f ${loci_dir}/lead_snps.tsv ] && cp ${loci_dir}/lead_snps.tsv "${meta.pair}/meta/"
    [ -f ${loci_dir}/locus_coords.tsv ] && cp ${loci_dir}/locus_coords.tsv "${meta.pair}/meta/" || true
    """
}