# WGS/WES Germline Secondary Analysis Pipeline

A clinical-grade Nextflow DSL2 pipeline for whole-genome and whole-exome sequencing germline variant analysis, developed for the Department of Genomic Medicine, National Cheng Kung University Hospital.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.x-brightgreen)](https://www.nextflow.io/)
[![Container: Apptainer](https://img.shields.io/badge/container-Apptainer-blue)](https://apptainer.org/)

---

## Overview

This pipeline performs secondary analysis of germline variants from short-read sequencing data (Illumina). It supports both WGS and WES modes, and is optimized for GPU-accelerated computing using NVIDIA Clara Parabricks.

Two pipeline entry points are provided:

| Entry point | SV caller | STR caller | License | Use case |
|-------------|-----------|------------|---------|----------|
| `main.nf` | **Delly** (BSD) | **GangSTR** (GPL v3) | ✅ GPL v3 compatible | Clinical / commercial use |
| `main_research.nf` | **Manta** (PolyForm Strict) | **ExpansionHunter** (PolyForm Strict) | ⚠️ Non-commercial only | Research use |

> **Why two pipelines?** Manta and ExpansionHunter are licensed under [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0/), which restricts commercial use. If your institution charges for sequencing services, use `main.nf` with Delly and GangSTR.

---

## Pipeline Flowcharts

### `main.nf` — Clinical Pipeline

```
FASTQ (R1, R2)
       │
       ▼
┌─────────────────────────────────────────────────────────┐
│  Step 0 · Preprocessing                                 │
│  FASTP (adapter trimming, QC)                           │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Step 1 · Alignment (GPU)                               │
│  Parabricks fq2bam (BWA-MEM2 + BQSR)                    │
└───────────────────────┬─────────────────────────────────┘
                        │
          ┌─────────────┴─────────────┐
          ▼                           ▼
┌──────────────────┐       ┌──────────────────────────────┐
│  Step 2 · QC     │       │  Step 3 · Parallel Calling   │
│  SAMtools stats  │       │                              │
│  Mosdepth        │       │  Lane 1 ── DeepVariant (GPU) │
└──────────────────┘       │  Lane 2a── HaplotypeCaller   │
                           │           (GPU)              │
                           │  Lane 2b── GATK VQSR         │
                           │           (WGS only)         │
                           │  Lane 3a── CNVkit            │
                           │  Lane 3b── DELLY ◀ clinical  │
                           │  Lane 3c── gCNV              │
                           │           (WES + PON only)   │
                           │  Lane 4 ── GANGSTR ◀ clinical│
                           │  Lane 5 ── GATK Mutect2      │
                           │           (mitochondria)     │
                           └──────────────┬───────────────┘
                                          │
                                          ▼
                           ┌──────────────────────────────┐
                           │  Step 4 · Post-processing    │
                           │  BCFtools Ensemble           │
                           │  (DV + HC/VQSR merge)        │
                           │  AutoMap ROH                 │
                           │  BCFtools Stats              │
                           │  MultiQC                     │
                           └──────────────────────────────┘
```

### Key Differences: `main.nf` vs `main_research.nf`

| Step | `main.nf` (clinical) | `main_research.nf` (research) |
|------|----------------------|-------------------------------|
| **SV calling** (Lane 3b) | Delly v1.7.3 — BSD | Manta v1.6.0 — PolyForm Strict ⚠️ |
| **STR calling** (Lane 4) | GangSTR v2.5.0 — GPL v3 | ExpansionHunter v5.0.0 — PolyForm Strict ⚠️ |
| **SV output** | `{sample}.delly.vcf.gz` | `manta_results/results/variants/diploidSV.vcf.gz` |
| **STR output** | `{sample}.str.vcf` + `.samplestats.tab` | `{sample}.str.vcf` + `.str.json` |
| **Config** | `nextflow_main.config` | `nextflow_main_research.config` |
| **Commercial use** | ✅ | ❌ |

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
awk 'BEGIN{OFS="\t"} /^chr([1-9]|1[0-9]|2[0-2])\t/{print $1, 0, $2}' \
    Homo_sapiens_assembly38.fasta.fai \
    > hg38_autosome_primary.bed

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

---

## Output Structure

```
{out_dir}/{SAMPLE_ID}/
├── 01_preprocessing/         FASTP QC reports
├── 02_alignment/             BAM + recalibration table
├── 03_alignment_qc/          SAMtools stats + Mosdepth
├── 04_snv_indel/             DeepVariant + HaplotypeCaller + Ensemble VCF
├── 05_cnv_sv/                CNVkit + Delly (or Manta) + gCNV
├── 06_repeat/                GangSTR (or ExpansionHunter) STR
├── 07_mitochondria/          mtDNA variants (Mutect2)
├── 08_roh/                   ROH regions (AutoMap)
└── pipeline_info/            Execution reports
```

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
---

## License and Third-party Tools

This pipeline is released under the [GNU General Public License v3](LICENSE) (GPL v3). You are free to use, modify, and distribute this pipeline, including for commercial purposes, provided that any derivative works are also released under GPL v3.

> ⚠️ **Research pipeline warning:** `main_research.nf` uses **Manta** and **ExpansionHunter**, both licensed under [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0/) which **prohibits commercial use**. If your institution charges for sequencing services, use `main.nf` (with Delly and GangSTR) instead, or obtain a separate commercial license from Illumina.

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
| [AutoMap](https://github.com/mquinodo/AutoMap) | 1.3 | MIT |
| [Manta](https://github.com/Illumina/manta) | 1.6.0 | PolyForm Strict 1.0.0 ⚠️ |
| [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | 5.0.0 | PolyForm Strict 1.0.0 ⚠️ |

Users are responsible for compliance with each tool's license terms.

---

## Citation

If you use this pipeline in your research, please cite the relevant tools listed above. Reference data is sourced from the [Broad Institute GCS bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references) subject to [Broad Institute data use terms](https://software.broadinstitute.org/gatk/download/bundle).
