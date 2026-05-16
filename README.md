# WGS/WES Germline Secondary Analysis Pipeline

A clinical-grade Nextflow DSL2 pipeline for whole-genome and whole-exome sequencing germline variant analysis, developed for the Department of Genomic Medicine and Neurology, National Cheng Kung University Hospital.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.x-brightgreen)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Container: Apptainer](https://img.shields.io/badge/container-Apptainer-blue)](https://apptainer.org/)

---

## Overview

This pipeline performs secondary analysis of germline variants from short-read sequencing data (Illumina). It supports both WGS and WES modes, and is optimized for GPU-accelerated computing using NVIDIA Clara Parabricks.

Two pipeline entry points are provided:

| Entry point | SV caller | STR caller | License | Use case |
|-------------|-----------|------------|---------|----------|
| `main.nf` | **Delly** (BSD) | **GangSTR** (GPL v3) | ✅ MIT-compatible | Clinical / commercial use |
| `main_research.nf` | **Manta** (PolyForm Strict) | **ExpansionHunter** (PolyForm Strict) | ⚠️ Non-commercial only | Research use |

> **Why two pipelines?** Manta and ExpansionHunter are licensed under [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0/), which restricts commercial use. If your institution charges for sequencing services, use `main.nf` with Delly and GangSTR.

---

## Pipeline Flowcharts

### `main.nf` — Clinical Pipeline (MIT-compatible)

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
┌──────────────────┐       ┌────────────────────────────┐
│  Step 2 · QC     │       │  Step 3 · Parallel         │
│  SAMtools stats  │       │  Variant Calling           │
│  Mosdepth        │       │                            │
└──────────────────┘       │  Lane 1 ─ DeepVariant(GPU) │
                           │  Lane 2a─ HaplotypeCaller  │
                           │         (GPU)              │
                           │  Lane 2b─ GATK VQSR        │
                           │         (WGS only)         │
                           │  Lane 3a─ CNVkit           │
                           │  Lane 3b─ DELLY ◀ new      │
                           │  Lane 3c─ gCNV             │
                           │         (WES + PON only)   │
                           │  Lane 4 ─ GANGSTR ◀ new    │
                           │  Lane 5 ─ GATK Mutect2     │
                           │         (mitochondria)     │
                           └──────────────┬─────────────┘
                                          │
                                          ▼
                           ┌────────────────────────────┐
                           │  Step 4 · Post-processing  │
                           │  BCFtools Ensemble         │
                           │  (DV + HC/VQSR merge)      │
                           │  AutoMap ROH               │
                           │  BCFtools Stats            │
                           │  MultiQC                   │
                           └────────────────────────────┘
```

### `main_research.nf` — Research Pipeline (non-commercial only)

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
┌──────────────────┐       ┌────────────────────────────┐
│  Step 2 · QC     │       │  Step 3 · Parallel         │
│  SAMtools stats  │       │  Variant Calling           │
│  Mosdepth        │       │                            │
└──────────────────┘       │  Lane 1 ─ DeepVariant(GPU) │
                           │  Lane 2a─ HaplotypeCaller  │
                           │         (GPU)              │
                           │  Lane 2b─ GATK VQSR        │
                           │         (WGS only)         │
                           │  Lane 3a─ CNVkit           │
                           │  Lane 3b─ MANTA ⚠️         │
                           │  Lane 3c─ gCNV             │
                           │         (WES + PON only)   │
                           │  Lane 4 ─ EXPANSIONHUNTER  │
                           │          ⚠️                │
                           │  Lane 5 ─ GATK Mutect2     │
                           │         (mitochondria)     │
                           └──────────────┬─────────────┘
                                          │
                                          ▼
                           ┌────────────────────────────┐
                           │  Step 4 · Post-processing  │
                           │  BCFtools Ensemble         │
                           │  (DV + HC/VQSR merge)      │
                           │  AutoMap ROH               │
                           │  BCFtools Stats            │
                           │  MultiQC                   │
                           └────────────────────────────┘

⚠️ Manta and ExpansionHunter are licensed under PolyForm Strict
   License 1.0.0 — non-commercial use only.
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

> **Note on STR output format:** GangSTR's `REPCN` field contains repeat copy numbers as a comma-separated string (e.g. `"13,13"`), while ExpansionHunter returns a tuple. Downstream tertiary analysis tools must handle this difference accordingly.

---

## Hardware Requirements

| Environment | CPU | GPU | RAM | Role |
|-------------|-----|-----|-----|------|
| Local (dev) | R9 9950X 16c | RTX PRO 6000 96GB | 128GB | Development & testing |
| DGM Server | Xeon w7-3565X 32c | RTX 2000 Ada 16GB | 125GB | Clinical deployment (single sample) |
| DGX-2 | Xeon Platinum 8168 48c | V100 × 6 (GPU 10–15) | 1.5TB | Batch processing |

---

## Quick Start

### 1. Load environment

```bash
source /path/to/pipeline_code/NGS2ndAnalysis_env.sh
```

### 2. Prepare samplesheet

```csv
sample,fastq_1,fastq_2,sex
SAMPLE001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,female
SAMPLE002,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,male
```

For multi-lane samples, add a `lane` column (`L001`, `L002`, etc.).

### 3. Run (clinical pipeline)

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

# Single-sample accelerated mode (6 GPUs on DGX-2)
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx_single \
    --input_csv /path/to/samplesheet.csv \
    --seq_type WGS \
    --out_dir /path/to/output \
    -resume
```

### 3a. Run (research pipeline)

```bash
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main_research.nf \
    -profile dgx \
    --input_csv /path/to/samplesheet.csv \
    --seq_type WES \
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
| `dgm` | DGM Server (RTX 2000 Ada) | Single GPU, no lock |
| `dgx` | DGX-2 batch mode | 1 GPU per sample, up to 6 parallel |
| `dgx_single` | DGX-2 single-sample | 6 GPUs for 1 sample (~40 min alignment) |

---

## Additional Reference Files Required

Beyond the standard Broad Institute hg38 bundle, the following files are needed:

```bash
# Delly exclude list (telomere/centromere regions)
wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv

# GangSTR TR reference (hg38)
wget https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz
gunzip hg38_ver17.bed.gz

# GangSTR WES-intersected TR reference
bedtools intersect \
    -a gangstr_hg38_ver17.bed \
    -b Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
    > gangstr_hg38_ver17_WES.bed
```

---

## Container List

### Clinical pipeline (`main.nf`)

```
parabricks_4.7.0-1.sif    (local)
parabricks_4.4.0.sif      (DGM / DGX-2)
gatk_4.6.2.0.sif
fastp_1.3.0.sif
samtools_1.23.1.sif
mosdepth_0.3.13.sif
bcftools_1.23.1.sif
multiqc_1.33.sif
delly_1.7.3.sif            ← Delly (replaces Manta)
cnvkit_0.9.12.sif
gangstr_2.5.0.sif          ← GangSTR (replaces ExpansionHunter)
bwa_0.7.19.sif
automap_1.3.sif
```

### Research pipeline (`main_research.nf`) — additional containers

```
manta_1.6.0.sif            ⚠️ PolyForm Strict — non-commercial only
expansionhunter_5.0.0.sif  ⚠️ PolyForm Strict — non-commercial only
```

---

## License and Third-party Tools

This pipeline is released under the [GNU General Public License v3](LICENSE) (GPL v3).
You are free to use, modify, and distribute this pipeline, including for commercial purposes,
provided that any derivative works are also released under GPL v3.

> ⚠️ **Research pipeline warning:** `main_research.nf` uses **Manta** and **ExpansionHunter**,
> both licensed under [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0/)
> which **prohibits commercial use**. If your institution charges for sequencing services,
> use `main.nf` (with Delly and GangSTR) instead, or obtain a separate commercial license from Illumina.

The pipeline orchestrates the following third-party tools, each subject to their own license terms:

| Tool | Version | License |
|------|---------|---------|
| [Nextflow](https://github.com/nextflow-io/nextflow) | ≥ 23.x | Apache 2.0 |
| [Apptainer](https://github.com/apptainer/apptainer) | ≥ 1.x | BSD 3-Clause |
| [NVIDIA Clara Parabricks](https://docs.nvidia.com/clara/parabricks/4.3.0/Documentation/EULA.html) | 4.4.0 / 4.7.0 | [NVIDIA AI Product Agreement](https://www.nvidia.com/en-us/data-center/products/nvidia-ai-enterprise/eula/) (free to use) |
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

If you use this pipeline in your research, please cite the relevant tools listed above. Reference data is sourced from the [Broad Institute GCS bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references) and is subject to [Broad Institute data use terms](https://software.broadinstitute.org/gatk/download/bundle).

---

*For detailed setup instructions, see the environment setup guide in the repository.*
