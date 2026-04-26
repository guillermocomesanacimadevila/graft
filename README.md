<p align="center">
  <img src="https://github.com/user-attachments/assets/24c5c653-8a83-4a47-9604-6427ab210c65"
       alt="nf-core-neurobridge_logo_light"
       width="1600">
</p>

---

## Introduction

**GRAFT** is a two-stage bioinformatics pipeline which takes two **genome-wide association study (GWAS) summary statistics** corresponding to two distinct traits (typically trait pairs where a pleiotropic relationship may be hypothesised). The pipeline covers analyses from estimating global heritability to functional approaches using **bulk and single-cell xQTL datasets**. In **Stage I**, neurobridge performs GWAS QC and computes **SNP heritability**, polygenicity and discoverability, **global genetic correlation** (across multiple methods), and local genetic correlation. It then identifies loci showing evidence of pleiotropic association (via **conjFDR**) and performs **Bayesian colocalisation** and **fine-mapping** to prioritise candidate causal variants/loci. Following this stage, loci can be annotated using **FUMA** (performed manually) to map variants to genes and obtain basic functional annotation. In **Stage II**, neurobridge focuses on functional interpretation of the prioritised loci by integrating bulk and single-cell xQTL datasets and performing **gene prioritisation (via SMR)** and consequent GWAS-to-xQTL colocalisation of priorisited genes. In essence, Stage II adds multi-layer biological evidence to the genetic findings from **Stage I**.

---

## Authors

**Guillermo Comesaña Cimadevila¹²³** · Alfie Thain⁴ · Marie-Joe Dib⁵ · Philip Taylor¹² · Jeremy Hall²³ · Dervis Salih¹⁶ · Nicholas Bray² · Emily Simmonds¹ · Valentina Escott-Price¹²†  

 ¹ UK Dementia Research Institute, UK  
 ² MRC Centre for Neuropsychiatric Genetics & Genomics, Cardiff University, Cardiff, UK  
 ³ MeOmics Technology Ltd, Cardiff, UK  
 ⁴ Division of Psychiatry, University College London, London, UK  
 ⁵ Nascent Studio, London, UK  
 ⁶ Department of Neuroscience, Physiology and Pharmacology, University College London, London, UK  

---

[![Python](https://img.shields.io/badge/Python-3.12%2B-3776AB)]() 
[![R](https://img.shields.io/badge/R-4.4%2B-276DC3)]() 
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.18986935-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.18986935)

---

## Pipeline

<img width="1247" height="1209" alt="graft" src="https://github.com/user-attachments/assets/e76e3f19-76f1-4dcd-a734-9c09820f0f64" />

---

## Steps

### Download reference data

> Download it **[here](https://doi.org/10.5281/zenodo.18986935)**.

---

### Adjust parameters

> Configure [`assets/params.stage1.yaml`] to your desired parameter set (**default parameters provided**).

---

### Build Docker image

```bash
cd neurobridge/
```

```bash
docker build --no-cache -t neurobridge:1 -f env/Dockerfile .
```

---

### Run neurobridge! (with Docker)

```bash
nextflow run . \
  -profile docker \
  -c conf/local/nextflow.config \
  -params-file assets/params.stage1.yaml 
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

---

## Tools

- **[LDSC](https://github.com/bulik/ldsc)** – Linkage Disequilibrium Score Regression  
- **[SumHer](https://github.com/zietzm/sumher_rs)** – Summary-based Heritability Estimation  
- **[HDL / HDL-L](https://github.com/zhenin/HDL)** – High-Definition Likelihood (global and local)  
- **[MiXeR](https://github.com/precimed/mixer)** – Mixture of Regressions  
- **[LAVA](https://github.com/josefin-werme/LAVA)** – Local analysis of genetic covariance  
- **[condFDR / conjFDR](https://github.com/alexploner/cfdr.pleio)** – Conditional / conjunctional false discovery rate  
- **[PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)** – Whole-genome association analysis toolkit  
- **[MAGMA](https://ctg.cncr.nl/software/magma)** – Gene-level and pathway enrichment analysis  
- **[COLOC](https://cran.r-project.org/web/packages/coloc/)** – Bayesian colocalisation analysis  
- **[SuSiE](https://github.com/stephenslab/susieR)** – Bayesian fine-mapping and credible set estimation  
- **[FUMA](https://fuma.ctglab.nl/)** – Functional mapping and annotation of GWAS loci  
- **[SMR + HEIDI](https://yanglab.westlake.edu.cn/software/smr/)** – Summary-data-based Mendelian randomization for GWAS–QTL integration  

---

## Get Started!

1. Install [`Nextflow`](https://www.nextflow.io/) `(>=25.10.0)`

2. Install [`Docker`](https://www.docker.com/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/)

3. Clone repo and run!
   
```bash
git clone https://github.com/guillermocomesanacimadevila/neurobridge.git
```

```bash
cd neurobridge/
```

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

```bash
nextflow run main.nf \
  -profile docker \
  -c conf/local/nextflow.config \
  --input assets/gwas.tsv \
  --pairs assets/ldsc_pairs.tsv \
  --outdir results
```

---

** Note ** 

> You can also run each method within [`workflows/neurobridge`], individually like so =>

```bash
nextflow run workflows/neurobridge/main_<input_method>.nf \
  -profile docker \
  -c conf/local/nextflow.config \
  -params-file assets/params.stage1.yaml 
```

---

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/mrcope for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

