# WGS/WES Germline Secondary Analysis Pipeline

A clinical-grade Nextflow DSL2 pipeline for whole-genome and whole-exome sequencing germline variant analysis, developed for the Department of Genomic Medicine and Neurology, National Cheng Kung University Hospital.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.x-brightgreen)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Container: Apptainer](https://img.shields.io/badge/container-Apptainer-blue)](https://apptainer.org/)

---

## Overview

This pipeline performs secondary analysis of germline variants from short-read sequencing data (Illumina). It supports both WGS and WES modes, and is optimized for GPU-accelerated computing using NVIDIA Clara Parabricks.

A single entry point (`main.nf`) is used, selected with `--seq_type WGS|WES`.
Behaviour is controlled by flags — the table shows their **current defaults**:

| Flag | Effect | Default | License |
|------|--------|---------|---------|
| `--run_phasing` | WhatsHap phasing + compound (MNV) merging, per caller, before the ensemble | **on** | — |
| `--run_gcnv` | GATK germline gCNV (WES only; requires a PON) | **on** | — |
| `--run_roh` | ROH via bcftools roh | **on** | MIT/GPL ✅ |
| `--run_manta` | Manta SV calling | off | PolyForm Strict 1.0.0 ⚠️ |
| `--run_expansionhunter` | ExpansionHunter STR | off | PolyForm Strict 1.0.0 ⚠️ |
| `--run_automap` | ROH via AutoMap | off | none published ⚠️ |

> **Commercial-safe by default.** The default path uses only commercially-usable tools
> (DeepVariant, HaplotypeCaller, Delly, CNVkit, gCNV, GangSTR, bcftools roh, mtDNA Mutect2).
> The three non-commercial tools (Manta / ExpansionHunter / AutoMap) stay **off** unless
> explicitly enabled. `--run_gcnv false` is an escape hatch for WES when the PON is not
> built yet.

---

## Pipeline Flowcharts

### `main.nf` — Clinical Pipeline

`main.nf` is composed of one sub-workflow per stage (`modules/*.nf`); the workflow
body is pure composition. Per-sample lanes below run in parallel off the same BAM.

```
FASTQ (R1, R2 — multi-lane per sample supported)
   │
   ▼  FASTP ── adapter trim + quality filter ─────────────────────► 01_preprocessing
   │
   ▼  PARABRICKS_FQ2BAM (GPU) ── BWA-MEM + BQSR, lanes merged ─────► 02_alignment
   │
   ├── ALIGNMENT_QC ──────────────────────────────────────────────► 03_alignment_qc
   │     samtools stats · mosdepth (WGS: autosome BED / WES: capture BED)
   │     PLOIDY_CHECK  ── sex + per-contig ploidy QC (mosdepth-based)
   │
   ├── CALL_SNV ──────────────────────────────────────────────────► 04_snv_indel
   │     DEEPVARIANT (GPU) · HAPLOTYPECALLER (GPU) · VQSR (WGS only)
   │     [--run_phasing] WhatsHap phase + combine_phased.py  (per caller)
   │        └─► BCFTOOLS_ENSEMBLE  (DV+HC union, +fixploidy) ─► *.ensemble.fixed.vcf.gz
   │
   ├── CALL_CNV_SV ───────────────────────────────────────────────► 05_cnv_sv
   │     DELLY (SV, PASS-only) · CNVKIT (WGS) · GCNV (WES, --run_gcnv) · [MANTA]
   │
   ├── CALL_STR ──────────────────────────────────────────────────► 06_repeat
   │     GangSTR (24-contig scatter → merge) · [ExpansionHunter]
   │
   ├── CALL_MITO ─────────────────────────────────────────────────► 07_mitochondria
   │     Mutect2 2-pass (normal + shifted) → liftover → merge → filter
   │
   └── CALL_ROH ──────────────────────────────────────────────────► 08_roh
         bcftools roh (--run_roh) · [AutoMap --run_automap]
                                          │
                                          ▼
   BCFTOOLS_STATS (DV + ensemble) ─► 09_postprocessing   ·   MULTIQC ─► pipeline_info/
```

### Callers by stage (sub-workflow)

| Stage | Default (always runs) | Optional (opt-in flag) |
|-------|-----------------------|------------------------|
| SNV/indel · `CALL_SNV` | DeepVariant + HaplotypeCaller → ensemble; VQSR (WGS only) | — |
| CNV/SV · `CALL_CNV_SV` | Delly (SV, BSD-3) + CNVkit **[WGS]** / gCNV **[WES]** | Manta — PolyForm Strict ⚠️ (`--run_manta`) |
| STR · `CALL_STR` | GangSTR v2.5.0 — GPL | ExpansionHunter — PolyForm Strict ⚠️ (`--run_expansionhunter`) |
| mtDNA · `CALL_MITO` | GATK Mutect2 (2-pass, shifted origin) | — |
| ROH · `CALL_ROH` | bcftools roh — MIT/GPL | AutoMap — no license ⚠️ (`--run_automap`) |

---

## Hardware Requirements

| Environment | CPU | GPU | RAM | Role |
|-------------|-----|-----|-----|------|
| Local (dev) | R9 9950X 16c | RTX PRO 6000 96GB | 128GB | Development & testing |
| DGM Server | Xeon w7-3565X 32c | RTX 2000 Ada 16GB | 125GB | Clinical deployment |
| DGX-2 | Xeon Platinum 8168 48c | V100 × 6 (GPU 10–15) | 1.5TB | Batch processing |

---

## Installation

### Step 1 — Install dependencies

**Option A: Online (DGM Server)**

```bash
# Install Apptainer
sudo apt update && sudo apt install -y apptainer

# Install Miniforge + Nextflow
sudo mkdir -p /opt/NGS2ndAnalysis
sudo chown -R $USER /opt/NGS2ndAnalysis
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p /opt/NGS2ndAnalysis/miniforge
source /opt/NGS2ndAnalysis/miniforge/bin/activate
mamba create -n NGS2ndAnalysis openjdk=17 nextflow procps-ng -y
```

**Option B: Offline (DGX-2)**

```bash
# Java 17: https://adoptium.net/temurin/releases/?os=linux&arch=x64&package=jdk&version=17
sudo mkdir -p /opt/java
sudo tar -xzf OpenJDK17U-jdk_x64_linux_hotspot_17.0.17_10.tar.gz -C /opt/java

# Apptainer .deb + dependencies: https://github.com/apptainer/apptainer/releases
sudo dpkg -i uidmap*.deb fakeroot*.deb libfakeroot*.deb apptainer_*.deb

# Nextflow: https://github.com/nextflow-io/nextflow/releases
sudo mkdir -p /opt/nextflow
sudo cp nextflow nextflow-*-dist /opt/nextflow/
sudo chmod +x /opt/nextflow/nextflow
sudo mv /opt/nextflow/nextflow-25.10.2-dist /opt/nextflow/nextflow-all.jar
```

### Step 2 — Build Apptainer containers

> ⚠️ Build containers **before** downloading reference data (mitochondria reference setup requires containers).

```bash
# Parabricks (requires NGC login: https://ngc.nvidia.com)
apptainer registry login -u '$oauthtoken' docker://nvcr.io
apptainer build parabricks_4.7.0-1.sif docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
apptainer build parabricks_4.4.0.sif docker://nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1
# Note: V100 (compute 7.0) only supports up to Parabricks 4.4.0

# GATK (requires docker.io login)
apptainer registry login -u <username> docker://docker.io
apptainer build gatk_4.6.2.0.sif docker://broadinstitute/gatk:4.6.2.0

# Preprocessing & QC
apptainer build fastp_1.3.0.sif      docker://quay.io/biocontainers/fastp:1.3.0--h43da1c4_0
apptainer build samtools_1.23.1.sif  docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0
apptainer build mosdepth_0.3.13.sif  docker://quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0
apptainer build bcftools_1.23.1.sif  docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0
apptainer build multiqc_1.33.sif     docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0

# CNV/SV (clinical pipeline)
apptainer build cnvkit_0.9.12.sif    docker://quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1
apptainer build delly_1.7.3.sif      docker://quay.io/biocontainers/delly:1.7.3--hd6466ae_0

# STR (clinical pipeline)
apptainer build gangstr_2.5.0.sif    docker://quay.io/biocontainers/gangstr:2.5.0--h7337834_10

# Mitochondria alignment
apptainer build bwa_0.7.19.sif       docker://quay.io/biocontainers/bwa:0.7.19--h577a1d6_1

# Phasing (default on: --run_phasing) + Python QC scripts
apptainer build whatshap_2.8.sif     docker://quay.io/biocontainers/whatshap:2.8--py39h2de1943_0
# combine_phased.py (compound merge) and ploidy_check.py (sex/ploidy QC) need only
# Python 3 stdlib + bcftools; the config names this image tertiary_python_1.0.0.sif.
# Reuse the tertiary repo's container, or build any Python-3-plus-bcftools image:
cat > tertiary_python.def << 'EOF'
Bootstrap: docker
From: python:3.11-slim
%post
    apt-get update && apt-get install -y bcftools && rm -rf /var/lib/apt/lists/*
EOF
apptainer build tertiary_python_1.0.0.sif tertiary_python.def

# ROH
# AutoMap is not on bioconda — build from .def file
cat > automap.def << 'EOF'
Bootstrap: docker
From: rocker/r-base:4.4.2
%post
    apt-get update && apt-get install -y wget bcftools bedtools perl git procps bc
    git clone https://github.com/mquinodo/AutoMap.git /opt/AutoMap
    chmod +x /opt/AutoMap/AutoMap_v1.3.sh
%runscript
    exec bash "$@"
EOF
apptainer build automap_1.3.sif automap.def

# Research pipeline only (non-commercial)
apptainer build manta_1.6.0.sif         docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6
apptainer build expansionhunter_5.0.0.sif docker://quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5
```

### Step 3 — Download reference data

```bash
REF_DIR="/path/to/reference/hg38"
mkdir -p ${REF_DIR}
cd ${REF_DIR}

# --- Broad Institute hg38 bundle ---
BASE_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/"
for FILE in \
    Homo_sapiens_assembly38.fasta \
    Homo_sapiens_assembly38.fasta.fai \
    Homo_sapiens_assembly38.dict \
    Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
    Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
    Homo_sapiens_assembly38.known_indels.vcf.gz \
    Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
    1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
    wgs_calling_regions.hg38.interval_list \
    hapmap_3.3.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz.tbi \
    1000G_omni2.5.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz.tbi \
    Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi \
    Homo_sapiens_assembly38.fasta.64.alt \
    Homo_sapiens_assembly38.fasta.64.amb \
    Homo_sapiens_assembly38.fasta.64.ann \
    Homo_sapiens_assembly38.fasta.64.bwt \
    Homo_sapiens_assembly38.fasta.64.pac \
    Homo_sapiens_assembly38.fasta.64.sa
do
    wget -c "${BASE_URL}${FILE}"
done

# BWA symlinks (required by Parabricks)
ln -sf Homo_sapiens_assembly38.fasta.64.bwt Homo_sapiens_assembly38.fasta.bwt
ln -sf Homo_sapiens_assembly38.fasta.64.pac Homo_sapiens_assembly38.fasta.pac
ln -sf Homo_sapiens_assembly38.fasta.64.ann Homo_sapiens_assembly38.fasta.ann
ln -sf Homo_sapiens_assembly38.fasta.64.amb Homo_sapiens_assembly38.fasta.amb
ln -sf Homo_sapiens_assembly38.fasta.64.sa  Homo_sapiens_assembly38.fasta.sa

# --- WES capture kit BED (Illumina Exome Panel v1.2) ---
wget "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/trusight/hg38/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"

# --- gCNV contig ploidy priors ---
wget "https://storage.googleapis.com/gatk-sv-resources-public/gcnv-exome/contig_ploidy_prior_hg38.tsv"

# --- Delly exclude list (telomere/centromere regions) ---
wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv

# --- GangSTR TR reference ---
wget https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz
gunzip hg38_ver17.bed.gz
# WES-intersected version (run after downloading capture BED)
bedtools intersect \
    -a hg38_ver17.bed \
    -b Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
    > gangstr_hg38_ver17_WES.bed
mv hg38_ver17.bed gangstr_hg38_ver17.bed

# --- Autosome BED for WGS depth QC (chr1-22 only) ---
# Used by mosdepth --by in WGS mode so the mean-coverage QC metric is not skewed by
# chrM (high copy), sex chromosomes, or decoy/unplaced contigs.
awk 'BEGIN{OFS="\t"} /^chr([1-9]|1[0-9]|2[0-2])\t/{print $1, 0, $2}' \
    Homo_sapiens_assembly38.fasta.fai \
    > hg38_autosome_primary.bed

# --- Sex/ploidy map for bcftools +fixploidy (single source of truth) ---
# Format: CHROM FROM TO SEX PLOIDY, with GRCh38 PAR coordinates. Versioned template
# lives in the repo at assets/sex_ploidy_GRCh38.txt; copy it into the reference dir.
cp /path/to/pipeline_code/assets/sex_ploidy_GRCh38.txt sex_ploidy_GRCh38.txt

# --- CNV blacklist (PAR + centromere + telomere) ---
# Download gap table from UCSC Table Browser:
# https://genome.ucsc.edu/cgi-bin/hgTables (track: gap, table: gap)
# Save as raw_blacklist.bed, then:
cat > raw_blacklist.bed << 'EOF'
chrX    10000   2781479
chrY    10000   2781479
chrX    155701383   156030895
chrY    56887903    57217415
# Paste UCSC gap table rows below:
EOF
# (add gap table rows, then run:)
awk '{print $1"\t"$2"\t"$3}' raw_blacklist.bed \
    | bedtools sort -i - | bedtools merge -i - \
    > hg38_clinical_blacklist.bed
grep -E "^chr([0-9]+|X|Y|M)\s" hg38_clinical_blacklist.bed \
    > hg38_clinical_blacklist.main.bed
rm raw_blacklist.bed

# --- Mappability (for CNV PON correction) ---
# Download from: https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz
wget https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz
zcat k100.umap.bed.gz | bedtools sort -i - | bedtools merge -i - \
    | bgzip > hg38_k100_umap_merged.bed.gz
tabix -p bed hg38_k100_umap_merged.bed.gz
apptainer exec /path/to/gatk_4.6.2.0.sif gatk IndexFeatureFile \
    -I hg38_k100_umap_merged.bed.gz
rm k100.umap.bed.gz

# --- Segmental duplication (for CNV PON correction) ---
# Download from UCSC Table Browser:
# https://genome.ucsc.edu/cgi-bin/hgTables (track: Segmental Dups, table: genomicSuperDups)
# Save as seg_dup.bed, then:
bedtools sort -i seg_dup.bed | bedtools merge -i - | bgzip > hg38_seg_dup.bed.gz
tabix -p bed hg38_seg_dup.bed.gz
rm seg_dup.bed
```

### Step 4 — Build mitochondria reference

> Several pitfalls documented from experience:
> - GATK 4.6+ renamed `ShiftFastaForMitochondria` to `ShiftFasta`
> - `ShiftFasta` only accepts single-chromosome FASTA (not whole genome)
> - The blacklist BED is stored in Git LFS — use `media.githubusercontent.com`, not `raw.githubusercontent.com`
> - `VariantFiltration --mask` requires a GATK `IndexFeatureFile` index

```bash
mkdir -p ${REF_DIR}/chrM

# Blacklist (Git LFS — must use media.githubusercontent.com)
wget -O ${REF_DIR}/chrM/blacklist_sites.hg38.chrM.bed \
    "https://media.githubusercontent.com/media/broadinstitute/gatk/master/src/test/resources/large/mitochondria_references/blacklist_sites.hg38.chrM.bed"

apptainer exec --bind /path/to/ref gatk_4.6.2.0.sif \
    gatk IndexFeatureFile -I ${REF_DIR}/chrM/blacklist_sites.hg38.chrM.bed

# Extract chrM FASTA
apptainer exec --bind /path/to/ref samtools_1.23.1.sif \
    samtools faidx ${REF_DIR}/Homo_sapiens_assembly38.fasta chrM \
    > ${REF_DIR}/chrM/chrM_only.fasta

apptainer exec --bind /path/to/ref samtools_1.23.1.sif \
    samtools faidx ${REF_DIR}/chrM/chrM_only.fasta
apptainer exec --bind /path/to/ref gatk_4.6.2.0.sif \
    gatk CreateSequenceDictionary -R ${REF_DIR}/chrM/chrM_only.fasta

# Shifted reference (for GATK mito pipeline)
apptainer exec --bind /path/to/ref gatk_4.6.2.0.sif \
    gatk ShiftFasta \
    -R ${REF_DIR}/chrM/chrM_only.fasta \
    -O ${REF_DIR}/chrM/chrM_shifted.fasta \
    --shift-back-output ${REF_DIR}/chrM/chrM_shift_back.chain

# BWA index for both references
apptainer exec --bind /path/to/ref bwa_0.7.19.sif \
    bwa index ${REF_DIR}/chrM/chrM_only.fasta
apptainer exec --bind /path/to/ref bwa_0.7.19.sif \
    bwa index ${REF_DIR}/chrM/chrM_shifted.fasta
```

Expected files in `${REF_DIR}/chrM/`:
- `chrM_only.fasta` + `.fai` + `.dict` + BWA index files
- `chrM_shifted.fasta` + `.fai` + `.dict` + BWA index files
- `chrM_shift_back.chain`
- `blacklist_sites.hg38.chrM.bed` + `.idx`

---

## Quick Start

### 1. Load environment

```bash
source /path/to/pipeline_code/NGS2ndAnalysis_env.sh
```

### 2. Prepare samplesheet

```csv
sample,fastq_1,fastq_2,sex,lane
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,female,L001
SAMPLE002,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,male,L001
```

For multi-lane samples, list each lane as a separate row with the same `sample` value and different `lane` values. `fq2bam` will merge lanes automatically with correct read group IDs.

### 3. Run

```bash
# WES + gCNV (requires PON)
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv /path/to/samplesheet.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir /path/to/output \
    -resume

# WGS
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv /path/to/samplesheet.csv \
    --seq_type WGS \
    --out_dir /path/to/output \
    -resume

# Single-sample accelerated mode (6 GPUs on DGX-2, ~40 min alignment)
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx_single \
    --input_csv /path/to/samplesheet.csv \
    --seq_type WGS \
    --out_dir /path/to/output \
    -resume
```
**Flag notes**

- Phasing (`--run_phasing`), ROH via bcftools (`--run_roh`), and — for WES — gCNV
  (`--run_gcnv`) are **on by default**; you do not need to pass them. `--run_gcnv true`
  above is explicit for clarity.
- Add `--run_manta --run_expansionhunter --run_automap` only to enable the three
  non-commercial research tools (all off by default).
- Escape hatches: `--run_phasing false` skips WhatsHap/compound merging;
  `--run_gcnv false` skips WES gCNV when the PON is not built yet.

---

## Output Structure

Per-sample outputs go under `{out_dir}/{SAMPLE_ID}/NN_stage/`; run-level reports go to
`{out_dir}/pipeline_info/`. Files tagged **[WGS]** / **[WES]** appear only in that mode;
**[opt]** only when the corresponding flag is on.

```
{out_dir}/
├── {SAMPLE_ID}/
│   ├── 01_preprocessing/     fastp trimmed reads + QC
│   ├── 02_alignment/         analysis-ready BAM + BQSR table + dup metrics
│   ├── 03_alignment_qc/      samtools stats, mosdepth depth, sex/ploidy QC
│   ├── 04_snv_indel/         DeepVariant, HaplotypeCaller, (VQSR), ensemble
│   ├── 05_cnv_sv/            Delly SV + CNVkit [WGS] / gCNV [WES]
│   ├── 06_repeat/            GangSTR STR
│   ├── 07_mitochondria/      chrM variants (Mutect2)
│   ├── 08_roh/               runs of homozygosity
│   └── 09_postprocessing/    per-VCF bcftools stats
└── pipeline_info/            MultiQC report + Nextflow execution reports
```

### 01_preprocessing — fastp
| File | Description |
|------|-------------|
| `<id>_{1,2}.fastp.fastq.gz` | adapter-trimmed, quality-filtered reads (one pair per lane if multi-lane) |
| `<id>.fastp.json` | machine-readable QC (consumed by MultiQC) |
| `<id>.fastp.html` | human-readable QC report |

Filtering: adapter auto-detect, `--cut_front/--cut_tail` mean Q20, min length 50, qualified Q15.

### 02_alignment — Parabricks fq2bam
| File | Description |
|------|-------------|
| `<id>.aligned.sorted.bam` (+ `.bai`) | BWA-MEM aligned, duplicate-marked, BQSR-applied; all lanes merged |
| `<id>.recal.txt` | BQSR recalibration table |
| `<id>.duplicate_metrics.txt` | duplicate metrics (MultiQC) |
| `qc_metrics_dir/` | Parabricks built-in QC (insert size, coverage) |

### 03_alignment_qc — samtools + mosdepth + ploidy
| File | Description |
|------|-------------|
| `<id>.stats` | `samtools stats` (mapping/error rate, insert size) |
| `<id>.mosdepth.summary.txt` | per-contig depth — columns `chrom length bases mean min max` |
| `<id>.mosdepth.global.dist.txt` | cumulative coverage distribution (MultiQC) |
| `<id>.regions.bed.gz` (+ `.csi`) | per-region mean depth (`--by` = autosome BED [WGS] / capture BED [WES]) |
| `<id>.thresholds.bed.gz` (+ `.csi`) | bases covered ≥ 1,10,15,20,30,50,100× per region |
| `<id>.ploidy.vcf.gz` | sex/ploidy QC, DRAGEN-style (see below) |
| `<id>.ploidy_qc.txt` | human-readable sex-check + per-contig NDC/RATIO + warnings |

**Sex/ploidy QC** (`ploidy_check.py`, mosdepth-based — used for **both** WGS and WES; for
WES, gCNV also computes its own internal contig-ploidy, which is *not* published here).
`ploidy.vcf.gz` has one record per **primary contig only** (chr1-22, X, Y, M):

| Field | Meaning |
|-------|---------|
| `FORMAT/DC` | mean depth of coverage |
| `FORMAT/NDC` | depth normalized to the **expected** ploidy for the estimated karyotype (~1.0 = as-expected; `.` for chrM, which is high-copy and not a ploidy unit) |
| `FORMAT/RATIO` | raw depth ÷ autosomal median (~1.0 diploid, ~0.5 male hemizygous chrX/Y) |
| `FILTER` | `PASS`, or `SUSPECT` when NDC deviates → possible aneuploidy |
| `##estimatedSexKaryotype` / `##referenceSexKaryotype` (header) | data-inferred vs samplesheet-declared karyotype (keys aligned with DRAGEN) |

Warn-only: a sex mismatch or aneuploidy prints a WARNING but never changes calls or fails the run.

### 04_snv_indel — SNV / indel
| File | Description |
|------|-------------|
| `<id>.deepvariant.vcf.gz` (+ `.tbi`) | DeepVariant calls |
| `<id>.haplotypecaller.vcf.gz` (+ `.tbi`) | HaplotypeCaller calls (GT/AD preserved; feeds ROH) |
| `<id>.vqsr_snp.vcf.gz`, `<id>.vqsr_indel.vcf.gz` **[WGS]** | VQSR-filtered HC (WES skips VQSR) |
| `<id>.snp.recal`, `<id>.snp.tranches` **[WGS]** | VQSR model + tranches |
| `<id>.ensemble.fixed.vcf.gz` (+ `.tbi`) | **main output** — DV + HC merged |

**`ensemble.fixed.vcf.gz`** has **two sample columns**, `<id>_DV` and `<id>_HC` (provenance =
which caller populated each genotype). It is biallelic-split, sex-ploidy corrected
(`bcftools +fixploidy`), and — with `--run_phasing` — phased with adjacent *cis* variants
merged into single MNVs. FORMAT `GT:GQ:DP:AD:VAF:PL:PS` (`AD`/`VAF` depth preserved on merged
records; `PS` = phase-set, `|` = phased genotype). This is the file tertiary analysis reads.

### 05_cnv_sv — copy number & structural variants
| File | Mode | Description |
|------|------|-------------|
| `<id>.delly.vcf.gz` (+ `.tbi`) | WGS + WES | Delly SV, **PASS-only**. INFO `SVTYPE` (DEL/DUP/INV/BND/INS), `END`, `SVLEN`; FORMAT includes `RDCN` (read-depth copy number) |
| `<id>.call.cns` | **[WGS]** | CNVkit absolute-CN segments — cols `chromosome start end gene log2 baf ci_hi ci_lo cn cn1 cn2 depth probes weight` (`cn` = total integer copy number, `cn1`/`cn2` = allele-specific CN, `baf` = B-allele frequency from the DeepVariant VCF; sex-aware; the file tertiary consumes) |
| `<id>.aligned.sorted.cns` / `.cnr` | **[WGS]** | CNVkit segmented / per-bin log2 ratios |
| `<id>-scatter.pdf` / `-diagram.pdf` | **[WGS]** | CNVkit plots |
| `<id>.gcnv.vcf.gz` (+ `.tbi`) | **[WES, --run_gcnv]** | GATK gCNV segments — FORMAT `GT:CN:NP:QA:QS:QSE:QSS` (`CN` = copy number, `NP` = # bins, `QS` = quality score) |
| `<id>.denoisedCR.tsv` | **[WES, --run_gcnv]** | gCNV denoised copy-ratio matrix |
| `manta_results/results/variants/diploidSV.vcf.gz` | **[opt --run_manta]** | Manta SV calls |

> Depth-CNV caller differs by mode: **WGS → CNVkit**, **WES → gCNV**; Delly (SV) runs in both.

### 06_repeat — short tandem repeats
| File | Description |
|------|-------------|
| `<id>.str.vcf` | GangSTR genotypes (24-contig scatter → merge). FORMAT `GT:DP:Q:REPCN:REPCI:…` (`REPCN` = repeat copy number per allele, `REPCI` = confidence interval) |
| `<id>.expansionhunter.{vcf,json}`, `<id>.expansionhunter_realigned.bam` | **[opt --run_expansionhunter]** ExpansionHunter genotypes + evidence BAMlet (named `.expansionhunter.*` so it never clobbers the GangSTR `.str.vcf`) |

### 07_mitochondria — chrM (Mutect2)
| File | Description |
|------|-------------|
| `<id>.mito.vcf.gz` (+ `.tbi`) | chrM variants: 2-pass (normal + shifted origin) → liftover → merge → FilterMutectCalls (`--mitochondria-mode`) + blacklist mask. FORMAT `GT:AD:AF:DP:F1R2:F2R1:FAD:SB` (`AF` = heteroplasmy fraction, `AD` = allele depths, `SB` = strand bias). `FILTER=PASS` = confident; `weak_evidence` / `strand_bias` / `base_qual` / `blacklisted_site` = filtered |

### 08_roh — runs of homozygosity
| File | Description |
|------|-------------|
| `<id>.roh.txt` | bcftools roh (`-O r` region format). Cols: `RG  Sample  Chromosome  Start  End  Length(bp)  #markers  Quality` |
| `<id>.HomRegions.tsv` / `.pdf` | **[opt --run_automap]** AutoMap ROH table + plot |

Both consume the HaplotypeCaller raw VCF (needs GT + AD; not VQSR/DeepVariant).

### 09_postprocessing — per-VCF stats
| File | Description |
|------|-------------|
| `<id>.deepvariant.vcf.stats` | `bcftools stats` on DeepVariant (counts, Ti/Tv, indel distribution) — MultiQC |
| `<id>.ensemble.fixed.vcf.stats` | `bcftools stats` on the ensemble VCF — MultiQC |

### pipeline_info/ (run-level)
| File | Description |
|------|-------------|
| `multiqc_report.html` (+ `multiqc_report_data/`) | aggregated QC across fastp, samtools, mosdepth, bcftools stats |
| Nextflow `report` / `timeline` / `trace` | run provenance (when enabled in the config) |

---

## Profiles

| Profile | Target | GPU strategy |
|---------|--------|--------------|
| `local` | Development machine | Single GPU, no lock |
| `dgm` | DGM Server (RTX 2000 Ada 16GB) | Single GPU, no lock |
| `dgx` | DGX-2 batch mode | 1 GPU per sample, up to 6 parallel |
| `dgx_single` | DGX-2 single-sample | 6 GPUs for 1 sample (~40 min alignment) |

---

## Building a gCNV Panel of Normals (PON)

WES gCNV calling requires a PON built from ≥30 normal samples (≥100 recommended). Build once; reuse for all subsequent samples.

```bash
# Prepare PON samplesheet (one row per sample, no multi-lane needed)
# Exclude: cancer samples, samples with known large CNVs, duplicates
# sample,fastq_1,fastq_2,sex

source /path/to/pipeline_code/NGS2ndAnalysis_env.sh
cd /raid/DGM/work

nextflow -c /path/to/nextflow_pon.config \
    run /path/to/main_pon.nf \
    --input_csv /path/to/pon_samplesheet.csv \
    --pon_out_dir /path/to/reference/hg38/gcnv_pon \
    -resume
```

After completion, copy `filtered.interval_list` to the PON directory:

```bash
INTERVAL=$(find /raid/DGM/work -name "filtered.interval_list" | head -1)
cp ${INTERVAL} /path/to/reference/hg38/gcnv_pon/filtered.interval_list
```

---

## Validation with NA12878

Download test data:

```bash
# WES (NA12878, HiSeq Exome)
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# WGS (NA12878, ERR194147 — use HTTPS, not ftp://)
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# GIAB benchmark VCF
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

Expected QC values:

| Metric | WES | WGS |
|--------|-----|-----|
| DeepVariant SNPs | ~270,000 | ~4,000,000–5,000,000 |
| HaplotypeCaller SNPs | ~240,000 | ~3,800,000–4,500,000 |
| Ti/Tv ratio | — | 2.0–2.1 |
| Mito PASS variants | 35–100 | 35–100 |
| Mapping rate | >99% | >98% |
| Mean depth | >100x (clinical) | >30x |

Structural / QC invariants to check on every run:

| Check | Expected |
|-------|----------|
| `ensemble.fixed.vcf.gz` sample columns | `<id>_DV` **and** `<id>_HC` |
| Ensemble FORMAT (with `--run_phasing`) | `GT:GQ:DP:AD:VAF:PL:PS` — `AD`/`VAF` non-empty, some genotypes phased (`\|`) |
| Delly `FILTER` | `PASS` only |
| CNV caller | WGS → CNVkit `call.cns`; WES → `gcnv.vcf.gz` |
| Sex-check (`ploidy_qc.txt`) | `estimated == declared` karyotype → `sex_check: OK` |

> Validated (2026-07) on NA12878 (WES) and an internal WGS sample after the sub-workflow
> refactor: all nine output stages populate; the ensemble preserves `AD`/`VAF` and phasing;
> the WGS/WES CNV split is correct; and sex-check returns the right karyotype (WES `XX`, WGS
> `XY`). Note: mosdepth reports a `*_region=0` row for contigs outside the `--by` BED, so
> `ploidy_check.py` uses the whole-contig mean when a contig has no on-target region — without
> this, WGS males were mis-called `X0?`.

---

## License and Third-party Tools

This pipeline is released under the [GNU General Public License v3](LICENSE) (GPL v3). You are free to use, modify, and distribute this pipeline, including for commercial purposes, provided that any derivative works are also released under GPL v3.

> ⚠️ **Optional non-commercial tools (all default off):** `--run_manta` and `--run_expansionhunter` enable Manta / ExpansionHunter, licensed under [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0/) (prohibits commercial use); `--run_automap` enables AutoMap, which publishes **no license** (all rights reserved). For ROH prefer `--run_roh` (bcftools roh, MIT/GPL, commercial-safe). Leave the non-commercial flags off for commercial/clinical service.

| Tool | Version | License |
|------|---------|---------|
| [Nextflow](https://github.com/nextflow-io/nextflow) | ≥ 23.x | Apache 2.0 |
| [Apptainer](https://github.com/apptainer/apptainer) | ≥ 1.x | BSD 3-Clause |
| [NVIDIA Clara Parabricks](https://docs.nvidia.com/clara/parabricks/4.3.0/Documentation/EULA.html) | 4.4.0 / 4.7.0 | [NVIDIA AI Product Agreement](https://www.nvidia.com/en-us/data-center/products/nvidia-ai-enterprise/eula/) — free, commercial use permitted |
| [GATK](https://github.com/broadinstitute/gatk) | 4.6.2.0 | Apache 2.0 |
| [fastp](https://github.com/OpenGene/fastp) | 1.3.0 | MIT |
| [SAMtools](https://github.com/samtools/samtools) | 1.23.1 | MIT |
| [BCFtools](https://github.com/samtools/bcftools) | 1.23.1 | MIT |
| [Mosdepth](https://github.com/brentp/mosdepth) | 0.3.13 | MIT |
| [Delly](https://github.com/dellytools/delly) | 1.7.3 | BSD ✅ |
| [CNVkit](https://github.com/etal/cnvkit) | 0.9.12 | Apache 2.0 |
| [GangSTR](https://github.com/gymreklab/GangSTR) | 2.5.0 | GPL v3 ✅ |
| [BWA](https://github.com/lh3/bwa) | 0.7.19 | GPL v3 |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 1.33 | GPL v3 |
| [AutoMap](https://github.com/mquinodo/AutoMap) | 1.3 | ⚠️ none published (all rights reserved) — opt-in `--run_automap` |
| [Manta](https://github.com/Illumina/manta) | 1.6.0 | PolyForm Strict 1.0.0 ⚠️ |
| [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | 5.0.0 | PolyForm Strict 1.0.0 ⚠️ |

Users are responsible for compliance with each tool's license terms.

---

## Citation

If you use this pipeline in your research, please cite the relevant tools listed above. Reference data is sourced from the [Broad Institute GCS bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references) subject to [Broad Institute data use terms](https://software.broadinstitute.org/gatk/download/bundle).
