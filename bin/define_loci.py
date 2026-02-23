from pathlib import Path
import pandas as pd

out_dir = Path("/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ")
lead_snp = "rs383902"
lead_snp2 = "rs11039131"
lead_snp3 = "rs11777131"
window = 500_000

ad_gwas = pd.read_csv("/Users/c24102394/Desktop/neurobridge/data/AD/post-qc/AD.ldsc_ready_neff.tsv", sep="\t")
scz_gwas = pd.read_csv("/Users/c24102394/Desktop/neurobridge/data/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv", sep="\t")
lon_gwas = pd.read_csv("/Users/c24102394/Desktop/neurobridge/data/LON/post-qc/LON.ldsc_ready_neff.tsv", sep="\t")

ad_gwas.rename(columns={"POS": "BP"}, inplace=True)
scz_gwas.rename(columns={"POS": "BP"}, inplace=True)
lon_gwas.rename(columns={"POS": "BP"}, inplace=True)

out_dir.mkdir(parents=True, exist_ok=True)

lead_snps = [
    (lead_snp,  ["AD", "SCZ"]),
    (lead_snp2, ["AD", "SCZ"]),
    (lead_snp3, ["AD", "SCZ", "LON"]),
]

for snp, traits in lead_snps:
    ad_lead = ad_gwas.loc[ad_gwas["SNP"] == snp].iloc[0]
    chr_ = int(ad_lead["CHR"])
    bp = int(ad_lead["BP"])
    start = max(0, bp - window)
    end = bp + window

    ad_loc = ad_gwas.loc[(ad_gwas["CHR"] == chr_) & (ad_gwas["BP"] >= start) & (ad_gwas["BP"] <= end)].copy()
    scz_loc = scz_gwas.loc[(scz_gwas["CHR"] == chr_) & (scz_gwas["BP"] >= start) & (scz_gwas["BP"] <= end)].copy()

    if traits == ["AD", "SCZ"]:
        common = ad_loc.merge(scz_loc[["SNP"]], on="SNP", how="inner")
        ad_common = ad_loc.loc[ad_loc["SNP"].isin(common["SNP"])].copy()
        scz_common = scz_loc.loc[scz_loc["SNP"].isin(common["SNP"])].copy()

        ad_common = ad_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        scz_common = scz_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)

        locus_dir = out_dir / f"locus_chr{chr_}_{start}_{end}"
        locus_dir.mkdir(parents=True, exist_ok=True)

        ad_common.to_csv(locus_dir / "gwas_AD.ldgwas.tsv", sep="\t", index=False)
        scz_common.to_csv(locus_dir / "gwas_SCZ.ldgwas.tsv", sep="\t", index=False)
        pd.DataFrame({"SNP": common["SNP"]}).sort_values("SNP").to_csv(locus_dir / "common_snps.tsv", sep="\t", index=False)

    else:
        lon_loc = lon_gwas.loc[(lon_gwas["CHR"] == chr_) & (lon_gwas["BP"] >= start) & (lon_gwas["BP"] <= end)].copy()

        common = ad_loc.merge(scz_loc[["SNP"]], on="SNP", how="inner").merge(lon_loc[["SNP"]], on="SNP", how="inner")
        ad_common = ad_loc.loc[ad_loc["SNP"].isin(common["SNP"])].copy()
        scz_common = scz_loc.loc[scz_loc["SNP"].isin(common["SNP"])].copy()
        lon_common = lon_loc.loc[lon_loc["SNP"].isin(common["SNP"])].copy()

        ad_common = ad_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        scz_common = scz_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        lon_common = lon_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)

        locus_dir = out_dir / f"locus_chr{chr_}_{start}_{end}"
        locus_dir.mkdir(parents=True, exist_ok=True)

        ad_common.to_csv(locus_dir / "gwas_AD.ldgwas.tsv", sep="\t", index=False)
        scz_common.to_csv(locus_dir / "gwas_SCZ.ldgwas.tsv", sep="\t", index=False)
        lon_common.to_csv(locus_dir / "gwas_LON.ldgwas.tsv", sep="\t", index=False)
        pd.DataFrame({"SNP": common["SNP"]}).sort_values("SNP").to_csv(locus_dir / "common_snps.tsv", sep="\t", index=False)
