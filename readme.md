# NGS WGS/WES Germline Analysis Pipeline - Environment Setup

```
# 本機開發環境
pylin1991@192.168.61.49

# DGM Server（部署目標）
n101569@192.168.84.91

# DGX-2（部署目標）
n101569@10.11.33.75
```

---
# Part 1：本機開發環境建立

## 1-1. 資料夾結構

```
# Reference（HDD，主要存放）
/data/pylin1991/GenomicReference/hg38/

# Reference（SSD，自動 rsync cache，pipeline 實際讀取）
/scratch/pylin1991/GenomicReference_Cache/hg38/

# 容器存放
/data/pylin1991/nf-containers/

# Pipeline 程式碼
/data/pylin1991/nf-containers/NGSSecondary/1_0_0/

# Nextflow 工作目錄
/scratch/pylin1991/nextflow_workspace/
├── home/        # NXF_HOME
├── work/        # NXF_WORK（中間檔）
├── temp/        # NXF_TEMP
├── apptainer_tmp/
└── apptainer_cache/

# 測試資料
/scratch/pylin1991/Pipeline_test/NA12878/
```

## 1-2. Conda nextflow 環境 (原本的nextflow環境)

本機使用 miniforge3，nextflow 環境位於 `/home/pylin1991/miniforge3/envs/nextflow/`。

activate 時自動執行以下腳本（已設定於 conda activate.d）：
- 設定所有 NXF_* 和 APPTAINER_* 環境變數
- 執行 `rsync -a --update` 將 HDD reference 同步到 SSD cache

```bash
conda activate nextflow
# 確認環境變數
env | grep -E "NXF|APPTAINER|JAVA"
```

## 1-3. 建立 Apptainer 容器

> ⚠️ 容器需在下載 Reference 之前建立，因為 1-6 的 Mitochondria reference 步驟需要用到容器。

```bash
mkdir -p /data/pylin1991/nf-containers
cd /data/pylin1991/nf-containers

# Parabricks 4.7.0（需要 NGC 登入）
apptainer registry login -u '$oauthtoken' docker://nvcr.io
apptainer build parabricks_4.7.0-1.sif \
    docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1

# DGX2 V100只能用到4.4.0
apptainer build parabricks_4.4.0.sif \
    docker://nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1

# GATK
# broadinstitute/gatk 在 docker.io，需要先登入才能拉取
apptainer registry login -u n101569 docker://docker.io
apptainer build gatk_4.6.2.0.sif \
    docker://broadinstitute/gatk:4.6.2.0

# 前處理與 QC
apptainer build fastp_1.3.0.sif \
    docker://quay.io/biocontainers/fastp:1.3.0--h43da1c4_0
apptainer build samtools_1.23.1.sif \
    docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0
apptainer build mosdepth_0.3.13.sif \
    docker://quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0

# 後處理
apptainer build bcftools_1.23.1.sif \
    docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0
apptainer build multiqc_1.33.sif \
    docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0

# Whatshap
apptainer build /data/pylin1991/nf-containers/whatshap_2.8.sif \
  docker://quay.io/biocontainers/whatshap:2.8--py39h2de1943_0
# 用途：對 NCKUH ensemble VCF 做 read-backed phasing 補 PS，供三級正確處理 compound
#       （相鄰 cis del+ins，如 SUZ12 delAAAinsTT）。以 --run_phasing true 開啟（預設關），
#       僅 NCKUH 路徑需要（DRAGEN VCF 自帶 PS）。此容器只含 whatshap；切檔/合併用既有
#       bcftools 容器，phase 用此容器，依 contig 平行（見 modules/phasing.nf）。

# Lane 3: SV/CNV（Parabricks 4.0+ 已移除這兩個工具）
apptainer build manta_1.6.0.sif \
    docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6
apptainer build cnvkit_0.9.12.sif \
    docker://quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1

# Delly
apptainer build delly_1.7.3.sif \
    docker://quay.io/biocontainers/delly:1.7.3--hd6466ae_0

# Lane 4: STR
apptainer build expansionhunter_5.0.0.sif \
    docker://quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5

# GangSTR
apptainer build gangstr_2.5.0.sif \
    docker://quay.io/biocontainers/gangstr:2.5.0--h7337834_10

# Lane 5: Mitochondria（bwa 單獨容器，chrM alignment 專用）
apptainer build bwa_0.7.19.sif \
    docker://quay.io/biocontainers/bwa:0.7.19--h577a1d6_1



# Post-processing: AutoMap ROH
# AutoMap 不在 bioconda，需自行建立容器
# 注意：apptainer 不能直接讀 Dockerfile，需改寫成 .def 格式
# 依賴：BCFtools、BEDTools、Perl、R
mkdir -p /tmp/automap_docker
cat > /tmp/automap_docker/automap.def << 'EOF'
Bootstrap: docker
From: rocker/r-base:4.4.2

%environment
    export AUTOMAP_HOME=/opt/AutoMap

%post
    apt-get update && apt-get install -y \
        wget bcftools bedtools perl git procps bc \
        && rm -rf /var/lib/apt/lists/*

    git clone https://github.com/mquinodo/AutoMap.git /opt/AutoMap
    chmod +x /opt/AutoMap/AutoMap_v1.3.sh

%runscript
    exec bash "$@"
EOF

APPTAINER_BIND="" apptainer build \
    /data/pylin1991/nf-containers/automap_1.3.sif \
    /tmp/automap_docker/automap.def

apptainer build /data/pylin1991/nf-containers/automap_1.3.sif /tmp/automap_docker/automap.def

# 確認容器清單
ls -lh /data/pylin1991/nf-containers/*.sif
```

完成後容器清單應為：

```
parabricks_4.7.0-1.sif
parabricks_4.4.0.sif
gatk_4.6.2.0.sif
fastp_1.3.0.sif
samtools_1.23.1.sif
mosdepth_0.3.13.sif
bcftools_1.23.1.sif
multiqc_1.33.sif
manta_1.6.0.sif
cnvkit_0.9.12.sif
expansionhunter_5.0.0.sif
bwa_0.7.19.sif
automap_1.3.sif
whatshap_2.8.sif
```

## 1-4. 下載 hg38 Reference

```bash
REF_DIR="/data/pylin1991/GenomicReference/hg38"
mkdir -p ${REF_DIR}
cd ${REF_DIR}

cat > download_refs.sh << 'EOF'
#!/bin/bash
DEST_DIR="/data/pylin1991/GenomicReference/hg38"
BASE_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/"

FILES=(
    "Homo_sapiens_assembly38.fasta"
    "Homo_sapiens_assembly38.fasta.fai"
    "Homo_sapiens_assembly38.dict"
    "Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
    "Homo_sapiens_assembly38.known_indels.vcf.gz"
    "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
    "wgs_calling_regions.hg38.interval_list"
    "hapmap_3.3.hg38.vcf.gz"
    "hapmap_3.3.hg38.vcf.gz.tbi"
    "1000G_omni2.5.hg38.vcf.gz"
    "1000G_omni2.5.hg38.vcf.gz.tbi"
    "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
    "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
    "Homo_sapiens_assembly38.fasta.64.alt"
    "Homo_sapiens_assembly38.fasta.64.amb"
    "Homo_sapiens_assembly38.fasta.64.ann"
    "Homo_sapiens_assembly38.fasta.64.bwt"
    "Homo_sapiens_assembly38.fasta.64.pac"
    "Homo_sapiens_assembly38.fasta.64.sa"
)

mkdir -p "$DEST_DIR"
cd "$DEST_DIR"
for file in "${FILES[@]}"; do
    wget -c --show-progress "${BASE_URL}${file}"
done
ls -lh
EOF

bash download_refs.sh
```

## 1-5. 建立 BWA Symbolic Links（Parabricks 需要）
```bash
cd /data/pylin1991/GenomicReference/hg38/
ln -sf Homo_sapiens_assembly38.fasta.64.bwt Homo_sapiens_assembly38.fasta.bwt
ln -sf Homo_sapiens_assembly38.fasta.64.pac Homo_sapiens_assembly38.fasta.pac
ln -sf Homo_sapiens_assembly38.fasta.64.ann Homo_sapiens_assembly38.fasta.ann
ln -sf Homo_sapiens_assembly38.fasta.64.amb Homo_sapiens_assembly38.fasta.amb
ln -sf Homo_sapiens_assembly38.fasta.64.sa  Homo_sapiens_assembly38.fasta.sa

# 同步到 SSD cache（之後 conda activate 會自動同步）
rsync -a --update /data/pylin1991/GenomicReference/ \
    /scratch/pylin1991/GenomicReference_Cache/
```

## 1-6. 下載 Mitochondria Reference（Lane 5 必備）

> ⚠️ 踩坑紀錄：
> 1. GATK 4.6+ 已將 `ShiftFastaForMitochondria` 改名為 `ShiftFasta`
> 2. `ShiftFasta` 只能接受單一染色體 fasta，不能傳入全基因組 fasta
> 3. `ShiftFasta` 執行前需要 `.fai` 和 `.dict`，必須先建好
> 4. `chrM_numt_regions.bed` 已無法從 Broad GCS 下載
> 5. GATK 容器和 Samtools 容器都沒有內建 bwa，需要獨立的 bwa 容器
> 6. 所有 apptainer exec 都需要加 `--bind /data`
> 7. blacklist BED 在 GitHub 用 Git LFS 儲存，必須用 `media.githubusercontent.com` 下載
> 8. VariantFiltration --mask 需要 GATK IndexFeatureFile 建立的 index，舊格式不相容

```bash
REF_DIR="/data/pylin1991/GenomicReference/hg38"
mkdir -p ${REF_DIR}/chrM

# Step 1: 下載 blacklist（BED 格式）
# 必須用 media.githubusercontent.com，raw.githubusercontent.com 只會下載 LFS 指標
wget -O ${REF_DIR}/chrM/blacklist_sites.hg38.chrM.bed \
    "https://media.githubusercontent.com/media/broadinstitute/gatk/master/src/test/resources/large/mitochondria_references/blacklist_sites.hg38.chrM.bed"

# 確認檔案大小合理（應為 ~132 bytes，包含 6 個已知 artifact 位點）
ls -lh ${REF_DIR}/chrM/blacklist_sites.hg38.chrM.bed

# 建立 BED index（VariantFiltration --mask 需要）
apptainer exec --bind /data /data/pylin1991/nf-containers/gatk_4.6.2.0.sif \
    gatk IndexFeatureFile \
    -I ${REF_DIR}/chrM/blacklist_sites.hg38.chrM.bed

# Step 2: 抽出 chrM only fasta
apptainer exec --bind /data /data/pylin1991/nf-containers/samtools_1.23.1.sif \
    samtools faidx ${REF_DIR}/Homo_sapiens_assembly38.fasta chrM \
    > ${REF_DIR}/chrM/chrM_only.fasta

# Step 3: 建立 fai 和 dict（ShiftFasta 的前置需求）
apptainer exec --bind /data /data/pylin1991/nf-containers/samtools_1.23.1.sif \
    samtools faidx ${REF_DIR}/chrM/chrM_only.fasta

apptainer exec --bind /data /data/pylin1991/nf-containers/gatk_4.6.2.0.sif \
    gatk CreateSequenceDictionary \
    -R ${REF_DIR}/chrM/chrM_only.fasta

# Step 4: 產生 shifted reference
apptainer exec --bind /data /data/pylin1991/nf-containers/gatk_4.6.2.0.sif \
    gatk ShiftFasta \
    -R ${REF_DIR}/chrM/chrM_only.fasta \
    -O ${REF_DIR}/chrM/chrM_shifted.fasta \
    --shift-back-output ${REF_DIR}/chrM/chrM_shift_back.chain

# Step 5: 建立 BWA index（正常版和 shifted 版各一份）
apptainer exec --bind /data /data/pylin1991/nf-containers/bwa_0.7.19.sif \
    bwa index ${REF_DIR}/chrM/chrM_only.fasta

apptainer exec --bind /data /data/pylin1991/nf-containers/bwa_0.7.19.sif \
    bwa index ${REF_DIR}/chrM/chrM_shifted.fasta

# Step 6: 同步到 SSD cache
rsync -a --update ${REF_DIR}/chrM/ \
    /scratch/pylin1991/GenomicReference_Cache/hg38/chrM/
```

產生後 `${REF_DIR}/chrM/` 目錄應包含：
- `chrM_only.fasta` + `.fai` + `.dict` + BWA index
- `chrM_shifted.fasta` + `.fai` + `.dict` + BWA index
- `chrM_shift_back.chain`
- `blacklist_sites.hg38.chrM.bed` + `.idx`

## 1-7. 下載其他 Reference 檔案

```bash
# ExpansionHunter variant catalog（STR 位點定義，Lane 4 必備）
wget -P /data/pylin1991/GenomicReference/hg38/ \
    https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/hg38/variant_catalog.json

# WES Capture Kit Target BED
# Illumina Exome Panel v1.2 (CEX)，對應 Cat. No. 15050026
# CNVkit WES hybrid 模式必須提供此檔案
wget -P /data/pylin1991/GenomicReference/hg38/ \
    "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/trusight/hg38/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"

# WGS primary chromosome bed
# 從 fasta.fai 產生 chr1-22 的完整 BED
awk 'BEGIN{OFS="\t"} /^chr([1-9]|1[0-9]|2[0-2])\t/{print $1, 0, $2}' \
    /data/pylin1991/GenomicReference/hg38/Homo_sapiens_assembly38.fasta.fai \
    > /data/pylin1991/GenomicReference/hg38/hg38_autosome_primary.bed

# Contig ploidy priors（gCNV 必備，GATK 官方公用檔案）
# 定義各染色體正常 copy number 的先驗機率，與樣本無關，直接使用官方版本
wget -P /data/pylin1991/GenomicReference/hg38/ \
    "https://storage.googleapis.com/gatk-sv-resources-public/gcnv-exome/contig_ploidy_prior_hg38.tsv"

# Delly需要exclude的地方
wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv \
    -O /data/pylin1991/GenomicReference/hg38/human.hg38.excl.tsv
```

```bash
# 建立blacklist，用來在cnv校正black list: PAR, centromere, telomere，或是之後再加入本實驗室常常CNV會false positive的地方，這些地方不call CNV
nano raw_blacklist.bed
```
```
chrX	10000	2781479
chrY	10000	2781479
chrX	155701383	156030895
chrY	56887903	57217415
# 下面再貼上UCSC gap table https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3894100131_aVK6G3X9e1lp4sDcciK0bMVZJKMY&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_gap&hgta_ctDesc=table+browser+query+on+gap&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED
```
```bash
awk '{print $1"\t"$2"\t"$3}' raw_blacklist.bed | bedtools sort -i - | bedtools merge -i - > hg38_clinical_blacklist.bed

# 只保留標準染色體的 blacklist
grep -E "^chr([0-9]+|X|Y|M)\s" \
    hg38_clinical_blacklist.bed \
    > hg38_clinical_blacklist.main.bed

rm raw_blacklist.bed
```
```bash
# 用來在cnv pon校正mappability
wget https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz
zcat k100.umap.bed.gz | bedtools sort -i - | bedtools merge -i - | bgzip > hg38_k100_umap_merged.bed.gz
tabix -p bed hg38_k100_umap_merged.bed.gz

apptainer exec /data/pylin1991/nf-containers/gatk_4.6.2.0.sif \
    gatk IndexFeatureFile \
    -I hg38_k100_umap_merged.bed.gz

ls hg38_k100_umap_merged.bed.gz*

rm k100.umap.bed.gz

# 用來在cnv pon校正segment duplication
nano seg_dup.bed
#貼上 UCSC seg dup table https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3894176447_46FlTUv6KJ3waAqmabK2Ex6DHlw9&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_genomicSuperDups&hgta_ctDesc=table+browser+query+on+genomicSuperDups&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED

bedtools sort -i seg_dup.bed | bedtools merge -i - | bgzip > hg38_seg_dup.bed.gz
tabix -p bed hg38_seg_dup.bed.gz
rm seg_dup.bed
```
```bash
# 同步到 SSD cache
rsync -a --update /data/pylin1991/GenomicReference/ \
    /scratch/pylin1991/GenomicReference_Cache/
```

## 1-8. Pipeline 程式碼結構

```
/data/pylin1991/nf-containers/NGSSecondary/1_0_0/
├── main.nf
├── nextflow_main.config
├── main_research.nf
├── nextflow_main_research.config
├── run_pipeline.sh
└── modules/
    ├── preprocessing.nf     (FASTP)
    ├── alignment.nf         (PARABRICKS_FQ2BAM)
    ├── alignment_qc.nf      (SAMTOOLS_STATS, MOSDEPTH)
    ├── variant_calling.nf   (Lane 1, 2a, 2b, 2c)
    ├── cnv_sv.nf            (Lane 3: CNVKIT_BATCH, MANTA_GERMLINE, gCNV)
    ├── repeat.nf            (Lane 4: EXPANSIONHUNTER)
    ├── mitochondria.nf      (Lane 5: MITO_*)
    └── postprocessing.nf    (BCFTOOLS_ENSEMBLE, BCFTOOLS_ROH, BCFTOOLS_STATS, MULTIQC)
```

```bash
mkdir -p /data/pylin1991/nf-containers/NGSSecondary/1_0_0/modules
cd /data/pylin1991/nf-containers/NGSSecondary/1_0_0
```
- main.nf
- nextflow_main.config
- modules/preprocessing.nf
- modules/alignment.nf
- modules/alignment_qc.nf
- modules/variant_calling.nf
- modules/cnv_sv.nf
- modules/repeat.nf
- modules/mitochondria.nf
- modules/postprocessing.nf
- modules/roh.nf

確認沒有crlf換行
```bash
find . -type f \( -name "*.nf" -o -name "*.config" -o -name "*.sh" \) \
    -exec sed -i 's/\r//' {} +
```

## 1-9. 測試資料

```bash
mkdir -p /scratch/pylin1991/Pipeline_test/NA12878
cd /scratch/pylin1991/Pipeline_test/NA12878

# WES HiSeq（NA12878，female）
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# WGS（NA12878，ERR194147）
# 注意：EBI 只支援 HTTPS，ftp:// 協定無法連線
# SRR622457 品質較差（Q20=Q30, 大量 N reads），不建議用於 pipeline 驗證
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# NIST 黃金標準 VCF（驗證用）
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

Samplesheet 格式（`sample,fastq_1,fastq_2,sex,lane`）

```bash
# samplesheetWES.csv
cat > /scratch/pylin1991/Pipeline_test/NA12878/samplesheetWES.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WES,/scratch/pylin1991/Pipeline_test/NA12878/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NA12878/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz,female
EOF

# samplesheetWGS.csv
cat > /scratch/pylin1991/Pipeline_test/NA12878/samplesheetWGS.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WGS,/scratch/pylin1991/Pipeline_test/NA12878/ERR194147_1.fastq.gz,/scratch/pylin1991/Pipeline_test/NA12878/ERR194147_2.fastq.gz,female
EOF

# samplesheetNCKUH.csv
cat > /scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/samplesheetVAL55.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L001_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L001_R2_001.fastq.gz,male,L001
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L002_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L002_R2_001.fastq.gz,male,L002
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L003_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L003_R2_001.fastq.gz,male,L003
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L004_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L004_R2_001.fastq.gz,male,L004
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L005_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L005_R2_001.fastq.gz,male,L005
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L006_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L006_R2_001.fastq.gz,male,L006
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L007_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L007_R2_001.fastq.gz,male,L007
VAL55,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L008_R1_001.fastq.gz,/scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/VAL-55_S47_L008_R2_001.fastq.gz,male,L008
EOF
```

## 1-10. 測試跑

```bash
tmux new -s pipeline_test

conda activate nextflow

cd /data/pylin1991/nf-containers/NGSSecondary/1_0_0

# 語法檢查
nextflow -c nextflow_main.config  run main.nf -profile local --input_csv /dev/null -preview

# 實際執行（WES）
nextflow -c nextflow_main.config \
    run main.nf \
    -profile local \
    --input_csv /scratch/pylin1991/Pipeline_test/NA12878/samplesheetWES.csv \
    --seq_type WES \
    --run_gcnv true \
    --run_manta --run_expansionhunter --run_automap \
    --out_dir /scratch/pylin1991/Pipeline_test/NA12878_WES_PON \
    -resume

# 實際執行（WGS）
nextflow -c nextflow_main.config \
    run main.nf \
    -profile local \
    --input_csv /scratch/pylin1991/Pipeline_test/NA12878/samplesheetWGS.csv \
    --seq_type WGS \
    --out_dir /scratch/pylin1991/Pipeline_test/NA12878_WGS \
    -resume

nextflow -c nextflow_main.config \
    run main.nf \
    -profile local \
    --input_csv /scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55/samplesheetVAL55.csv \
    --run_manta --run_expansionhunter --run_automap \
    --run_phasing true \
    --seq_type WGS \
    --out_dir /scratch/pylin1991/Pipeline_test/NCKUH_WGS_VAL55 \
    -resume
```

## 1-11. 輸出結果驗證

```bash
# WES 驗證
SAMPLE="NA12878_WES"
OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878_WES_PON/${SAMPLE}"
BCFTOOLS="apptainer exec /data/pylin1991/nf-containers/bcftools_1.23.1.sif bcftools"

# WGS 驗證時改為：
# SAMPLE="NA12878_WGS"
# OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878/${SAMPLE}"

# in DGX2 BCFTOOLS="apptainer exec /datalake_Intermediate/pipeline/nextflow_containers/bcftools_1.23.1.sif bcftools"
# in DGM
```

### Variant Count

```bash
echo "=== DeepVariant ==="
$BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | grep "^SN"

# WGS：VQSR 後；WES：直接 HaplotypeCaller 輸出（無 VQSR）
if [ -f "${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz" ]; then
    echo "=== HaplotypeCaller (post-VQSR, WGS) ==="
    $BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.vqsr_indel.vcf.gz | grep "^SN"
else
    echo "=== HaplotypeCaller (WES, no VQSR) ==="
    $BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.haplotypecaller.vcf.gz | grep "^SN"
fi

echo "=== Ensemble ==="
$BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.ensemble.fixed.vcf.gz | grep "^SN"

echo "=== Mitochondria PASS ==="
$BCFTOOLS view -f PASS ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz | grep -v "^#" | wc -l

echo "=== Mito FILTER breakdown ==="
$BCFTOOLS view ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz | grep -v "^#" \
    | awk '{print $7}' | sort | uniq -c | sort -rn

```

**預期值（NA12878，高品質資料）：**

| 工具 | SNPs | Indels | 備註 |
|------|------|--------|------|
| DeepVariant WES | ~270,000 | ~39,000 | 單 lane 資料 |
| HaplotypeCaller WES | ~240,000 | ~33,000 | |
| DeepVariant WGS | ~4,000,000–5,000,000 | ~700,000–900,000 | |
| HaplotypeCaller WGS (VQSR) | ~3,800,000–4,500,000 | ~700,000–900,000 | |
| Mitochondria PASS | 35–100 | — | |

> ⚠️ 若 WGS 資料來自 SRR622457，DeepVariant SNPs 可能高達 930 萬，Ti/Tv 可能偏低至 1.73，此為該資料品質問題，非 pipeline 錯誤。建議改用 ERR194147 進行 WGS 驗證。

### Alignment QC

```bash
echo "=== Mosdepth Summary ==="
cat ${OUTDIR}/03_alignment_qc/${SAMPLE}.mosdepth.summary.txt

echo "=== Mapping Rate ==="
grep -E "^SN.*(raw total sequences|reads mapped:)" \
    ${OUTDIR}/03_alignment_qc/${SAMPLE}.stats

echo "=== Error Rate & Read Length ==="
grep -E "^SN.*(error rate|average length)" \
    ${OUTDIR}/03_alignment_qc/${SAMPLE}.stats
```

**預期值：**

| 指標 | WES 預期 | WGS 預期 |
|------|----------|----------|
| Mapping rate | >99% | >98% |
| Mean depth (target) | >100x（臨床）/ 24x（單 lane 測試）| >30x |
| Error rate | <0.3% | <0.3% |

### Ti/Tv Ratio（WGS 特有）

```bash
if [ -f "${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz" ]; then
    echo "=== VQSR Ti/Tv ==="
    $BCFTOOLS stats ${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz | grep "^TSTV"
else
    echo "=== WES 模式，跳過 Ti/Tv（無 VQSR）==="
fi
```

**預期值：**

| 範圍 | 評估 |
|------|------|
| 2.0–2.1 | ✅ 正常 WGS |
| <1.9 | ⚠️ 假陽性偏多，確認資料品質 |
| >2.2 | ⚠️ 可能僅計算 coding region |

### Fastp QC 解析

```bash
cat ${OUTDIR}/01_preprocessing/${SAMPLE}.fastp.json \
    | python3 -c "
import json, sys
d = json.load(sys.stdin)
bf = d['summary']['before_filtering']
fr = d['filtering_result']
print(f'Q20 (before): {bf[\"q20_rate\"]:.3f}')
print(f'Q30 (before): {bf[\"q30_rate\"]:.3f}')
print(f'GC  (before): {bf[\"gc_content\"]:.3f}')
print(f'Passed reads: {fr[\"passed_filter_reads\"]}')
print(f'Low quality:  {fr[\"low_quality_reads\"]}')
print(f'Too short:    {fr[\"too_short_reads\"]}')
print(f'Too many N:   {fr[\"too_many_N_reads\"]}')
"
```

**注意：fastp 1.0+ after filtering 的 Q20/Q30 會顯示為 1.0**，屬正常現象，請以 before filtering 數值評估資料品質。

**預期值：**

| 指標 | 良好資料 |
|------|----------|
| Q20 (before) | >95% |
| Q30 (before) | >90% |
| Passed rate | >95% |

---

### CNV / SV / STR / ROH 驗證

```bash
SAMPLE="NA12878_WES"
OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878_WES_PON/${SAMPLE}"
BCFTOOLS="apptainer exec /data/pylin1991/nf-containers/bcftools_1.23.1.sif bcftools"

echo "=== CNVkit CN 分布（第 7 欄為絕對 CN）==="
grep -v "^chromosome" ${OUTDIR}/05_cnv_sv/${SAMPLE}.call.cns \
    | awk '{print $7}' | sort | uniq -c | sort -rn | head -10

echo "=== Delly PASS SV ==="
$BCFTOOLS view -f PASS \
    ${OUTDIR}/05_cnv_sv/${SAMPLE}.delly.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Delly PASS SV:"

echo "=== Delly SV type 分布 ==="
$BCFTOOLS view -f PASS \
    ${OUTDIR}/05_cnv_sv/${SAMPLE}.delly.vcf.gz \
    | grep -v "^#" | grep -oP 'SVTYPE=\K\w+' | sort | uniq -c | sort -rn

echo "=== STR (GangSTR) ==="
grep -v "^#" ${OUTDIR}/06_repeat/${SAMPLE}.str.vcf | wc -l | xargs echo "STR loci genotyped:"

echo "=== Mitochondria ==="
$BCFTOOLS view -f PASS ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Mito PASS:"
$BCFTOOLS view ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | awk '{print $7}' | sort | uniq -c | sort -rn | head -5

echo "=== ROH ==="
cat ${OUTDIR}/08_roh/${SAMPLE}.HomRegions.tsv
```

```bash
SAMPLE="NA12878_WES"
OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878_WES_PON/${SAMPLE}"
BCFTOOLS="apptainer exec /data/pylin1991/nf-containers/bcftools_1.23.1.sif bcftools"

# echo "=== gCNV ==="
# $BCFTOOLS view ${OUTDIR}/05_cnv_sv/${SAMPLE}.gcnv.vcf.gz \
#     | grep -v "^#" | wc -l | xargs echo "gCNV total:"
# $BCFTOOLS view -f PASS ${OUTDIR}/05_cnv_sv/${SAMPLE}.gcnv.vcf.gz \
#     | grep -v "^#" | wc -l | xargs echo "gCNV PASS:"

echo "=== CNVkit CN 分布（第 7 欄為絕對 CN）==="
grep -v "^chromosome" ${OUTDIR}/05_cnv_sv/${SAMPLE}.call.cns \
    | awk '{print $7}' | sort | uniq -c | sort -rn | head -10

echo "=== Manta PASS SV ==="
$BCFTOOLS view -f PASS \
    ${OUTDIR}/05_cnv_sv/manta_results/results/variants/diploidSV.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Manta PASS SV:"

echo "=== STR ==="
grep -v "^#" ${OUTDIR}/06_repeat/${SAMPLE}.str.vcf | wc -l | xargs echo "STR loci:"

echo "=== Mitochondria ==="
$BCFTOOLS view -f PASS ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Mito PASS:"
$BCFTOOLS view ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | awk '{print $7}' | sort | uniq -c | sort -rn | head -5

echo "=== ROH ==="
cat ${OUTDIR}/08_roh/${SAMPLE}.HomRegions.tsv
```

**WES NA12878 單 lane（24x）預期值：**

| 項目 | 預期 | 備註 |
|------|------|------|
| gCNV PASS | 0 | 正常樣本 + 低深度，臨床樣本深度足夠時才有 PASS |
| CNVkit CN=2 | 主要 | 低深度時 CN 估計不穩定，臨床樣本再驗證 |
| Manta PASS SV | ~61 | |
| STR loci | ~38 | WES capture 外的 loci 無法偵測 |
| Mito PASS | 35-100 | |
| ROH total | <100 Mb | NA12878 非近親，不應有大片 ROH |


---
# Part 2：移植到 DGX-2

> 前提：local已完整跑通。DGX-2 完全離線，所有檔案需從本機傳入。

## 2-1. DGX-2 環境需求

| 項目 | 規格 |
|------|------|
| CPU | Xeon Platinum 8168（48 cores）|
| GPU | V100 × 16（分配 GPU 10-15，各 32GB VRAM）|
| RAM | 1.5TB |
| OS | Ubuntu |
| 帳號 | n101569@10.11.33.75 |

## 2-2. 建立資料夾結構

```bash
ssh n101569@10.11.33.75

mkdir -p /datalake_Intermediate/pipeline/reference/hg38
mkdir -p /datalake_Intermediate/pipeline/pipeline_code
mkdir -p /datalake_Intermediate/pipeline/nextflow_containers
mkdir -p /datalake_Intermediate/pipeline/nextflow_home
mkdir -p /datalake_Intermediate/pipeline/nextflow_output
mkdir -p /datalake_Intermediate/pipeline/install
mkdir -p /raid/DGM/reference
mkdir -p /raid/DGM/work
mkdir -p /raid/DGM/nextflow_temp
mkdir -p /raid/DGM/apptainer_temp
mkdir -p /raid/DGM/pytensor_cache
```

## 2-3. 安裝 JAVA（離線）

```bash
# 從本機下載並傳送
# https://adoptium.net/zh-CN/temurin/releases?os=linux&arch=x64&package=jdk&version=17

scp OpenJDK17U-jdk_x64_linux_hotspot_17.0.17_10.tar.gz \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/install/

# 在 DGX-2 執行
sudo mkdir -p /opt/java
cd /datalake_Intermediate/pipeline/install
sudo tar -xzf OpenJDK17U-jdk_x64_linux_hotspot_17.0.17_10.tar.gz -C /opt/java
ls /opt/java
```

## 2-4. 安裝 Apptainer（離線）

```bash
# 下載以下套件後傳送到 DGX-2：
# - apptainer .deb：https://github.com/apptainer/apptainer/releases
# - uidmap .deb：https://launchpad.net/ubuntu/jammy/+package/uidmap
# - libfakeroot .deb：https://launchpad.net/ubuntu/jammy/amd64/libfakeroot/1.28-1ubuntu1
# - fakeroot .deb：https://launchpad.net/ubuntu/jammy/amd64/fakeroot/1.28-1ubuntu1

scp *.deb n101569@10.11.33.75:/datalake_Intermediate/pipeline/install/

# 在 DGX-2 執行
cd /datalake_Intermediate/pipeline/install
sudo dpkg -i uidmap*.deb fakeroot*.deb libfakeroot*.deb
sudo dpkg -i apptainer_*.deb
apptainer --version
```

## 2-5. 安裝 Nextflow（離線）

```bash
# 下載 nextflow 和 nextflow-dist jar：
# https://github.com/nextflow-io/nextflow/releases (25.10.2)

scp nextflow nextflow-*-dist \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/install/

# 在 DGX-2 執行
sudo mkdir -p /opt/nextflow
cd /datalake_Intermediate/pipeline/install
sudo cp nextflow nextflow-*-dist /opt/nextflow/
sudo chmod +x /opt/nextflow/nextflow
sudo mv /opt/nextflow/nextflow-25.10.2-dist /opt/nextflow/nextflow-all.jar
```

## 2-6. 傳送容器（從本機）

```bash
scp /data/pylin1991/nf-containers/*.sif \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/nextflow_containers/
```

## 2-7. 傳送 Reference（從本機）

```bash
# scp -r /data/pylin1991/GenomicReference/hg38/* \
#     n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/
rsync -avz --progress \
    /data/pylin1991/GenomicReference/hg38/ \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/
```

## 2-8. 傳送 Pipeline 程式碼（從本機）

```bash
scp -r /data/pylin1991/nf-containers/NGSSecondary/1_0_0/* \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/pipeline_code/
```

## 2-10. 傳送測試資料並執行測試

```bash
# 在本機執行：傳送 NA12878 測試 FASTQ 到 DGM
scp /scratch/pylin1991/Pipeline_test/NA12878/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz \
    /scratch/pylin1991/Pipeline_test/NA12878/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz \
    /scratch/pylin1991/Pipeline_test/NA12878/ERR194147_1.fastq.gz \
    /scratch/pylin1991/Pipeline_test/NA12878/ERR194147_2.fastq.gz \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/nextflow_output/NA12878/
```

```bash
# 在 DGM 執行：建立輸出資料夾和 samplesheet
mkdir -p /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES
mkdir -p /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS

cat > /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES/samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WES,/datalake_Intermediate/pipeline/nextflow_output/NA12878/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz,/datalake_Intermediate/pipeline/nextflow_output/NA12878/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz,female
EOF

cat > /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS/samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WGS,/datalake_Intermediate/pipeline/nextflow_output/NA12878/ERR194147_1.fastq.gz,/datalake_Intermediate/pipeline/nextflow_output/NA12878/ERR194147_2.fastq.gz,female
EOF

# 執行測試（WES）
tmux 

source /datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh

cd /raid/DGM/work

nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf -profile dgx --input_csv /dev/null -preview


nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES/samplesheet.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES \
    -resume

nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgx_single \
    --input_csv /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES/samplesheet.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output/NA12878_WES \
    -resume

# 執行測試（WGS）
nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS/samplesheet.csv \
    --seq_type WGS \
    --run_gcnv false \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS \
    -resume

nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgx_single \
    --input_csv /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS/samplesheet.csv \
    --seq_type WGS \
    --run_gcnv false \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output/NA12878_WGS \
    -resume
```


---
# Part 3：建立 gCNV Panel of Normals（PON）

> 前提：DGX-2 環境已完整建立（Part 3 完成）。
> PON 建立只需要做一次，之後所有新樣本都用同一個 PON 跑 case mode。

## 3-1. 準備 PON Samplesheet

請依照以下格式準備 CSV 檔案：

```csv
sample,fastq_1,fastq_2,sex
SAMPLE001,/path/to/SAMPLE001_R1.fastq.gz,/path/to/SAMPLE001_R2.fastq.gz,female
SAMPLE002,/path/to/SAMPLE002_R1.fastq.gz,/path/to/SAMPLE002_R2.fastq.gz,unknown
```
```bash
# 取最近的150男150女
python3 subsample_pon.py /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/pon_samplesheet.csv -n 150 --by-sex --systematic -o /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/pon_150.csv
```
**注意事項：**
- `sample` 欄位每個樣本名稱必須唯一
- `sex` 填 `male`、`female` 或 `unknown`；不確定可填 `unknown`，模型會從 chrX/chrY depth 自動推斷，但建議盡量填正確性別以提高 chrX/chrY CNV 的準確度
- 路徑必須是 DGX-2 上的完整絕對路徑
- 檔案格式必須是 `.fastq.gz`
- **排除**：癌症樣本、已確診大片段 CNV 的樣本、重複樣本只留一個

## 3-2 gpu lock避免nextflow分配錯誤

- gpu_lock.sh
- gpu_unlock.sh 
- NGS2ndAnalysis_env.sh

## 3-3. 建立 PON Pipeline 程式碼

- main_pon.nf
- nextflow_pon.config
- modules/pon.nf

## 3-4. 建立輸出目錄

```bash
mkdir -p /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon
mkdir -p /raid/DGM/pon_work
mkdir -p /raid/DGM/pon_temp
```

## 3-5. 執行 PON 建立

```bash
tmux new -s gcnv_pon

cd /raid/DGM/pon_work

source /datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh

#nextflow -c /datalake_Intermediate/pipeline/pipeline_code/nextflow_pon.config \
#    run /datalake_Intermediate/pipeline/pipeline_code/main_pon.nf \
#    --input_csv /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/pon_samplesheet.csv \
#    --pon_out_dir /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon \
#    -work-dir /raid/DGM/pon_work \
#    -resume

nextflow -c /datalake_Intermediate/pipeline/pipeline_code/nextflow_pon.config \
    run /datalake_Intermediate/pipeline/pipeline_code/main_pon.nf \
    --input_csv /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/pon_150.csv \
    --pon_out_dir /datalake_Intermediate/pipeline/reference/hg38/pon_output \
    -work-dir /raid/DGM/pon_work \
    -resume
```

## 3-6. 確認 PON 輸出

### Step 1：在 DGX-2 確認 PON 輸出完整

```bash
PON_DIR="/datalake_Intermediate/pipeline/reference/hg38/gcnv_pon"
PON_DIR="/datalake_Intermediate/pipeline/reference/hg38/pon_output"

echo "=== 目錄結構 ==="
du -sh ${PON_DIR}/*/

echo "=== ploidy_model ==="
ls -lh ${PON_DIR}/gcnv_model/ploidy_model/

echo "=== gcnv_model shard ==="
ls -lh ${PON_DIR}/gcnv_model/shards/
ls -lh ${PON_DIR}/gcnv_model/shards/gcnv_model_shard_0/
# shard 已改為各自 index（gcnv_model_shard_0..38，內層 cohort_0-model..cohort_38-model）；
# 舊版曾全部撞名成 gcnv_model_shard_scattered（見踩坑 #33）。全部 *-model 應為 39：
find ${PON_DIR}/gcnv_model/shards -maxdepth 2 -type d -name '*-model' | wc -l

echo "=== cnvkit_reference ==="
ls -lh ${PON_DIR}/cnvkit_reference/

echo "=== filtered.interval_list ==="
# 從 pon_work 找到並複製到正式路徑
find /raid/DGM/pon_work -name "filtered.interval_list" | head -1
```

### Step 2: 把 filtered.interval_list 複製到正式路徑
```bash
# 在 DGX-2 執行
INTERVAL_PATH=$(find /raid/DGM/pon_work -name "filtered.interval_list" | head -1)
echo "找到：${INTERVAL_PATH}"

cp ${INTERVAL_PATH} \
    /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/filtered.interval_list

# 確認
ls -lh /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/filtered.interval_list
wc -l /datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/filtered.interval_list
# 應該約 196,286 行
```

## 3-7. 從 DGX-2 拉回 gCNV PON（PON 建立完成後執行）

```bash
# 在本機執行
mkdir -p /data/pylin1991/GenomicReference/hg38/gcnv_pon

rsync -avz --progress \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/ \
    /data/pylin1991/GenomicReference/hg38/gcnv_pon/

# gcnv_model（含 ploidy_model 和 shards）
rsync -avz --progress \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/gcnv_model/ \
    /data/pylin1991/GenomicReference/hg38/gcnv_pon/gcnv_model/

# cnvkit_reference
rsync -avz --progress \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/cnvkit_reference/ \
    /data/pylin1991/GenomicReference/hg38/gcnv_pon/cnvkit_reference/

# filtered.interval_list
scp n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/gcnv_pon/filtered.interval_list \
    /data/pylin1991/GenomicReference/hg38/gcnv_pon/filtered.interval_list

# 同步到 SSD cache
rsync -a --update \
    /data/pylin1991/GenomicReference/hg38/gcnv_pon/ \
    /scratch/pylin1991/GenomicReference_Cache/hg38/gcnv_pon/

# 確認整體大小（應約 2.8GB）
du -sh /scratch/pylin1991/GenomicReference_Cache/hg38/gcnv_pon/

# 確認 model shard 路徑（pipeline 會用 *-model glob）
ls /scratch/pylin1991/GenomicReference_Cache/hg38/gcnv_pon/gcnv_model/shards/gcnv_model_shard_0/
# 應該看到 cohort_0-model 和 cohort_0-tracking（每個 shard 各自 index：gcnv_model_shard_0..38）
# 全部 *-model 數量應為 39（= scatter 分片數）：
# find .../gcnv_pon/gcnv_model/shards -maxdepth 2 -type d -name '*-model' | wc -l

# 確認 filtered.interval_list
wc -l /scratch/pylin1991/GenomicReference_Cache/hg38/gcnv_pon/filtered.interval_list
```

## 3-8 測試pon使用

```bash
# 實際執行（WES + gCNV，需先完成 Part 4 建立 PON 並拉回本機）

tmux

conda activate nextflow
mkdir /scratch/pylin1991/Pipeline_test/NA12878_WES

cd /scratch/pylin1991/nextflow_workspace/work

nextflow -c /data/pylin1991/nf-containers/NGSSecondary/1_0_0/nextflow_main.config \
    run /data/pylin1991/nf-containers/NGSSecondary/1_0_0/main.nf \
    -profile local \
    --input_csv /scratch/pylin1991/Pipeline_test/NA12878/samplesheetWES.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir /scratch/pylin1991/Pipeline_test/NA12878_WES \
    -resume
```

## 3-9. 清理中間檔

```bash
# PON 建立完成確認無誤後，清理 work 目錄（節省空間）
rm -rf /raid/DGM/pon_work
```

## 3-10. 輸出結果驗證

```bash
# WES 驗證
SAMPLE="NA12878_WES"
OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878_WES_PON/${SAMPLE}"
BCFTOOLS="apptainer exec /data/pylin1991/nf-containers/bcftools_1.23.1.sif bcftools"

# WGS 驗證時改為：
# SAMPLE="NA12878_WGS"
# OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878/${SAMPLE}"
```

### Variant Count

```bash
echo "=== DeepVariant ==="
$BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | grep "^SN"

# WGS：VQSR 後；WES：直接 HaplotypeCaller 輸出（無 VQSR）
if [ -f "${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz" ]; then
    echo "=== HaplotypeCaller (post-VQSR, WGS) ==="
    $BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.vqsr_indel.vcf.gz | grep "^SN"
else
    echo "=== HaplotypeCaller (WES, no VQSR) ==="
    $BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.haplotypecaller.vcf.gz | grep "^SN"
fi

echo "=== Ensemble ==="
$BCFTOOLS stats ${OUTDIR}/04_snv_indel/${SAMPLE}.ensemble.fixed.vcf.gz | grep "^SN"

echo "=== Mitochondria PASS ==="
$BCFTOOLS view -f PASS ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz | grep -v "^#" | wc -l

echo "=== Mito FILTER breakdown ==="
$BCFTOOLS view ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz | grep -v "^#" \
    | awk '{print $7}' | sort | uniq -c | sort -rn

```

**預期值（NA12878，高品質資料）：**

| 工具 | SNPs | Indels | 備註 |
|------|------|--------|------|
| DeepVariant WES | ~270,000 | ~39,000 | 單 lane 資料 |
| HaplotypeCaller WES | ~240,000 | ~33,000 | |
| DeepVariant WGS | ~4,000,000–5,000,000 | ~700,000–900,000 | |
| HaplotypeCaller WGS (VQSR) | ~3,800,000–4,500,000 | ~700,000–900,000 | |
| Mitochondria PASS | 35–100 | — | |

> ⚠️ 若 WGS 資料來自 SRR622457，DeepVariant SNPs 可能高達 930 萬，Ti/Tv 可能偏低至 1.73，此為該資料品質問題，非 pipeline 錯誤。建議改用 ERR194147 進行 WGS 驗證。

### Alignment QC

```bash
echo "=== Mosdepth Summary ==="
cat ${OUTDIR}/03_alignment_qc/${SAMPLE}.mosdepth.summary.txt

echo "=== Mapping Rate ==="
grep -E "^SN.*(raw total sequences|reads mapped:)" \
    ${OUTDIR}/03_alignment_qc/${SAMPLE}.stats

echo "=== Error Rate & Read Length ==="
grep -E "^SN.*(error rate|average length)" \
    ${OUTDIR}/03_alignment_qc/${SAMPLE}.stats
```

**預期值：**

| 指標 | WES 預期 | WGS 預期 |
|------|----------|----------|
| Mapping rate | >99% | >98% |
| Mean depth (target) | >100x（臨床）/ 24x（單 lane 測試）| >30x |
| Error rate | <0.3% | <0.3% |

### Ti/Tv Ratio（WGS 特有）

```bash
if [ -f "${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz" ]; then
    echo "=== VQSR Ti/Tv ==="
    $BCFTOOLS stats ${OUTDIR}/04_variant_calling/${SAMPLE}.vqsr_indel.vcf.gz | grep "^TSTV"
else
    echo "=== WES 模式，跳過 Ti/Tv（無 VQSR）==="
fi
```

**預期值：**

| 範圍 | 評估 |
|------|------|
| 2.0–2.1 | ✅ 正常 WGS |
| <1.9 | ⚠️ 假陽性偏多，確認資料品質 |
| >2.2 | ⚠️ 可能僅計算 coding region |

### Fastp QC 解析

```bash
cat ${OUTDIR}/01_preprocessing/${SAMPLE}.fastp.json \
    | python3 -c "
import json, sys
d = json.load(sys.stdin)
bf = d['summary']['before_filtering']
fr = d['filtering_result']
print(f'Q20 (before): {bf[\"q20_rate\"]:.3f}')
print(f'Q30 (before): {bf[\"q30_rate\"]:.3f}')
print(f'GC  (before): {bf[\"gc_content\"]:.3f}')
print(f'Passed reads: {fr[\"passed_filter_reads\"]}')
print(f'Low quality:  {fr[\"low_quality_reads\"]}')
print(f'Too short:    {fr[\"too_short_reads\"]}')
print(f'Too many N:   {fr[\"too_many_N_reads\"]}')
"
```

**注意：fastp 1.0+ after filtering 的 Q20/Q30 會顯示為 1.0**，屬正常現象，請以 before filtering 數值評估資料品質。

**預期值：**

| 指標 | 良好資料 |
|------|----------|
| Q20 (before) | >95% |
| Q30 (before) | >90% |
| Passed rate | >95% |

---

### CNV / SV / STR / ROH 驗證

```bash
SAMPLE="NA12878_WES"
OUTDIR="/scratch/pylin1991/Pipeline_test/NA12878_WES_PON2/${SAMPLE}"
BCFTOOLS="apptainer exec /data/pylin1991/nf-containers/bcftools_1.23.1.sif bcftools"

echo "=== gCNV ==="
$BCFTOOLS view ${OUTDIR}/05_cnv_sv/${SAMPLE}.gcnv.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "gCNV total:"
$BCFTOOLS view -f PASS ${OUTDIR}/05_cnv_sv/${SAMPLE}.gcnv.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "gCNV PASS:"

echo "=== CNVkit CN 分布（第 7 欄為絕對 CN）==="
grep -v "^chromosome" ${OUTDIR}/05_cnv_sv/${SAMPLE}.call.cns \
    | awk '{print $7}' | sort | uniq -c | sort -rn | head -10

echo "=== Manta PASS SV ==="
$BCFTOOLS view -f PASS \
    ${OUTDIR}/05_cnv_sv/manta_results/results/variants/diploidSV.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Manta PASS SV:"

echo "=== STR ==="
grep -v "^#" ${OUTDIR}/06_repeat/${SAMPLE}.str.vcf | wc -l | xargs echo "STR loci:"

echo "=== Mitochondria ==="
$BCFTOOLS view -f PASS ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | wc -l | xargs echo "Mito PASS:"
$BCFTOOLS view ${OUTDIR}/07_mitochondria/${SAMPLE}.mito.vcf.gz \
    | grep -v "^#" | awk '{print $7}' | sort | uniq -c | sort -rn | head -5

echo "=== ROH ==="
cat ${OUTDIR}/08_roh/${SAMPLE}.HomRegions.tsv
```

**WES NA12878 單 lane（24x）預期值：**

| 項目 | 預期 | 備註 |
|------|------|------|
| gCNV PASS | 0 | 正常樣本 + 低深度，臨床樣本深度足夠時才有 PASS |
| CNVkit CN=2 | 主要 | 低深度時 CN 估計不穩定，臨床樣本再驗證 |
| Manta PASS SV | ~61 | |
| STR loci | ~38 | WES capture 外的 loci 無法偵測 |
| Mito PASS | 35-100 | |
| ROH total | <100 Mb | NA12878 非近親，不應有大片 ROH |


---
# Part 4：移植到 DGM Server

> 前提：本機 pipeline 已在測試中完全跑通。

## 4-1. DGM Server 環境需求

| 項目 | 規格 |
|------|------|
| CPU | Xeon w7-3565X（32 cores）|
| GPU | RTX 2000 Ada（16GB VRAM）|
| RAM | 125GB |
| OS | Ubuntu |
| 帳號 | n101569@192.168.84.91 |

## 4-2. 建立資料夾結構

```bash
ssh n101569@192.168.84.91

mkdir -p /home/pipeline/reference/hg38
mkdir -p /home/pipeline/pipeline_code
mkdir -p /home/pipeline/nextflow_containers
mkdir -p /home/pipeline/nextflow_home
mkdir -p /home/pipeline/nextflow_output
mkdir -p /home/pipeline/pipeline_info
mkdir -p /home/pipeline/nextflow_temp
mkdir -p /home/pipeline/apptainer_temp
mkdir -p /home/pipeline/work
mkdir -p /home/pipeline/pytensor_cache


# 1. 把擁有者改為 n101569，群組改為 dgm_nckuh
sudo chown -R n101569:dgm_nckuh /home/pipeline

# 2. 設定權限：
# 擁有者 rwx，群組 rwx，其他人 r-x
# setgid (2) 讓新建立的檔案自動繼承 dgm_nckuh 群組
sudo chmod -R 2775 /home/pipeline

# 3. 確認結果
ls -la /home/ | grep pipeline
ls -la /home/pipeline/

sudo usermod -aG dgm_nckuh <新使用者帳號>
```

## 4-3. 安裝 Apptainer

```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
apptainer --version
```

## 4-4. 安裝 Miniforge 與 Nextflow

```bash
sudo mkdir -p /opt/NGS2ndAnalysis
sudo chown -R n101569:n101569 /opt/NGS2ndAnalysis
sudo chmod -R 755 /opt/NGS2ndAnalysis

mkdir ~/Download && cd ~/Download
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p /opt/NGS2ndAnalysis/miniforge
chmod -R o+rx /opt/NGS2ndAnalysis/miniforge

source /opt/NGS2ndAnalysis/miniforge/bin/activate

conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda config --remove channels defaults
conda config --show channels

mamba create -n NGS2ndAnalysis openjdk=17 nextflow procps-ng -y
conda init bash
```

## 4-5. 傳送容器（從本機）

```bash
scp /data/pylin1991/nf-containers/*.sif \
    n101569@192.168.84.91:/home/pipeline/nextflow_containers/
```

## 4-6. 傳送 Reference（從本機）

```bash
# scp -r /data/pylin1991/GenomicReference/hg38/* \
#     n101569@192.168.84.91:/home/pipeline/reference/hg38/

rsync -avz --progress \
    /data/pylin1991/GenomicReference/hg38/ \
    n101569@192.168.84.91:/home/pipeline/reference/hg38/
```

## 4-7. 傳送 Pipeline 程式碼（從本機）

```bash
scp -r /data/pylin1991/nf-containers/NGSSecondary/1_0_0/* \
    n101569@192.168.84.91:/home/pipeline/pipeline_code/
```

## 4-8. 傳送測試資料並執行測試

```bash
# 在 DGM 執行：建立輸出資料夾和 samplesheet
ssh n101569@192.168.84.91

mkdir -p /home/pipeline/nextflow_output/NA12878_WES
mkdir -p /home/pipeline/nextflow_output/NA12878_WGS

cat > /home/pipeline/nextflow_output/NA12878_WES/samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WES,/home/datalake_Intermediate/pipeline/nextflow_output/NA12878/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz,/home/datalake_Intermediate/pipeline/nextflow_output/NA12878/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz,female
EOF

cat > /home/pipeline/nextflow_output/NA12878_WGS/samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,sex,lane
NA12878_WGS,/home/datalake_Intermediate/pipeline/nextflow_output/NA12878/ERR194147_1.fastq.gz,/home/datalake_Intermediate/pipeline/nextflow_output/NA12878/ERR194147_2.fastq.gz,female
EOF

```
```bash
# WES
tmux 

source /home/datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh

cd /home/pipeline/work

nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgm \
    --input_csv /home/pipeline/nextflow_output/NA12878_WES/samplesheet.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir /home/pipeline/nextflow_output/NA12878_WES \
    -resume

# 執行測試（WGS）
nextflow -c ${PIPELINE_CONFIG} \
    run ${PIPELINE_CODE}/main.nf \
    -profile dgm \
    --input_csv /home/pipeline/nextflow_output/NA12878_WGS/samplesheet.csv \
    --seq_type WGS \
    --run_gcnv false \
    --out_dir /home/pipeline/nextflow_output/NA12878_WGS \
    -resume
```

---

# Appendix：Pipeline 計畫

## 分析流程

```
FASTQ (R1, R2)
    ↓
[Step 1] Preprocessing
    ├── FASTP（adapter removal, quality filter）
    ├── Parabricks fq2bam（GPU alignment + BQSR）
    └── Samtools stats + Mosdepth（alignment QC）
         ↓
[Step 2] Parallel Variant Calling（五路並進）
    ├── Lane 1: Parabricks DeepVariant（GPU）→ deepvariant.vcf.gz
    ├── Lane 2a: Parabricks HaplotypeCaller（GPU）→ haplotypecaller.vcf.gz
    │   └── Lane 2b: GATK VQSR（WGS only）→ vqsr_snp/indel.vcf.gz
    ├── Lane 3a: CNVkit（WGS/WES）→ CNV cns/cnr
    ├── Lane 3b: Manta（WGS/WES）→ SV VCF
    ├── Lane 3c: gCNV（WES only，需 PON）→ gcnv.vcf.gz
    ├── Lane 4: ExpansionHunter（CPU）→ str.vcf
    └── Lane 5: GATK Mutect2 mito mode（CPU）→ mito.vcf.gz
         ↓
[Step 3] Post-processing
    ├── BCFtools Ensemble（合併 DV + HC/VQSR）
    ├── BCFtools ROH（隱性遺傳診斷輔助）
    ├── BCFtools Stats（VCF QC）
    └── MultiQC（整合報告）
```

## Variant Classification

CNV、SV 和 Mitochondria 的 variant classification 留給三級分析：
- CNV/SV：需對照 OMIM/ClinVar/DGV 資料庫
- Mitochondria：需對照 MITOMAP 資料庫
- 二級分析目標是產生乾淨可信的 VCF

## 模式切換

```bash
# WGS（預設）
--seq_type WGS

# WES
--seq_type WES

# WES + gCNV（需要 PON）
--seq_type WES --run_gcnv true \
    --gcnv_pon_dir /path/to/pon.hdf5 \
    --gcnv_model_dir /path/to/gcnv_model \
    --gcnv_ploidy_model_dir /path/to/ploidy_model
```

---

# 踩坑紀錄彙整

1. **Parabricks 4.0+ 已移除 Manta 和 CNVkit**，需用獨立容器
2. **CNNScoreVariants 在 GATK 4.6.1.0 移除**，WES 直接用 DeepVariant + HaplotypeCaller Ensemble
3. **GATK ShiftFasta**（原 ShiftFastaForMitochondria）只能接受單一染色體 fasta
4. **chrM_numt_regions.bed** 已無法從 Broad GCS 下載
5. **--median-autosomal-coverage** 在新版 GATK 已移除
6. **--blacklisted-sites** 在 GATK 4.6 已移除，改用 VariantFiltration --mask
7. **blacklist BED 在 GitHub 用 Git LFS 儲存**，必須用 `media.githubusercontent.com` 下載
8. **VariantFiltration --mask 需要 GATK IndexFeatureFile 建立的 index**，舊格式不相容
9. **Nextflow process 只能指定一個容器**，chrM alignment 需拆成三個 process
10. **BWA index 檔案需明確宣告在 input 裡**，Nextflow 不會自動帶入同目錄的 index 檔
11. **WES mosdepth** 需傳入 capture BED（`--by`），否則 `_region` 統計等同整條染色體
12. **EBI FTP 只支援 HTTPS**，`ftp://` 協定無法連線
13. **SRR622457 資料品質問題**：Q20=Q30、大量 too_short reads、Ti/Tv 僅 1.73，建議改用 ERR194147
14. **fastp 1.0+ after filtering Q20/Q30 = 1.0**：屬正常現象，請以 before filtering 數值評估
15. **Broad GCS bucket 已更換**：舊網址 `genomics-public-data` 已停用，請改用 `gcp-public-data--broad-references`
16. **docker.io 需要登入**才能用 apptainer 拉取 broadinstitute/gatk，quay.io 不需要
17. Parabricks 4.7.0 不支援 V100（compute 7.0），DGX-2 需使用 4.4.0
18. AutoMap 執行時需要對 Resources 目錄有寫入權限，需先 cp -r /opt/AutoMap ./AutoMap_local
19. AutoMap 需要未壓縮的 VCF（先用 bcftools view 解壓）
20. CUDA_VISIBLE_DEVICES 在 Apptainer 容器內不生效，需在 singularity runOptions 加 --env CUDA_VISIBLE_DEVICES=...
21. Nextflow local executor 的 process_gpu maxForks 是 per-process 限制，不同 process 間不互相等待，需用 lock file 或 channel dependency 控制 GPU 使用順序
22. bcftools fixploidy plugin 需設定 BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools
23. DGX-2 執行 nextflow 時需在沒有 nextflow.config 的目錄下執行，或將 nextflow.config 改名
24. **gCNV case mode 的 --model 需指向 *-model 子目錄**，不是 shard 根目錄（例如應指向 `gcnv_model_shard_0/cohort_0-model`，而非 `gcnv_model_shard_0`）
25. **BCFTOOLS_STATS 被呼叫兩次（DV + Ensemble）時輸出檔名會撞名**，需從 VCF 檔名自動產生 stats 檔名（`vcf.name.replace('.vcf.gz', '.vcf.stats')`）
26. **COLLECT_GATK_COUNTS input 需宣告 fasta_fai 和 fasta_dict**，否則 GATK 找不到 .fai index
27. **FILTER_INTERVALS 的 -L 參數需用 preprocessed.interval_list**，不能用 annotated.tsv（GATK 不認識 .tsv 格式作為 interval）
28. **pon.nf 的 CNVKIT_REFERENCE input 需宣告 fasta_fai**，CNVkit 計算 GC content 時需要
29. **PLOIDY_COHORT 在容器內需要寫入 ~/.pytensor/compiledir**，需在 singularity runOptions 加 `--env PYTENSOR_FLAGS=compiledir=/raid/DGM/pytensor_cache` 並預先建立該目錄
30. **PON samplesheet 不能有重複 sample ID**，多 lane 樣本需只保留一個（PON 不需要合併 lane）；重複樣本會造成 FILTER_INTERVALS 的 input file name collision
31. **gCNV model 的 ch_model_shards glob 應為 `*-model`（單層）**，不是 `*/*-model`（雙層）
32. PYTENSOR_FLAGS 在 DGX-2 的所有 gCNV 相關 process 都需要（PON 和 case mode），原因是 DGX-2 的 singularity runOptions 沒有 bind /home，導致容器內 /home/n101569 是唯讀的。本機因為有 bind /home 所以不受影響。
33. **GCNV_COHORT scatter shard 撞名**：IntervalListTools 切出的每個 shard 檔名都叫 `scattered.interval_list`；若用 `interval_shard.baseName` 命名輸出，39 個 shard 會全部輸出到 `gcnv_model_shard_scattered` / `cohort_scattered-model` 互相覆蓋 → 模型只剩 1 個 shard。解法：channel 帶入 index（`tuple val(idx), path(interval_shard)`），用 idx 命名 → `gcnv_model_shard_0..38`、`cohort_0-model..cohort_38-model`。case mode 的 `*-model` glob 靠這些唯一名稱才收得齊 39 個。
34. **FilterMutectCalls 在 GATK 4.6 沒有 `--autosomal-coverage`**（舊版才有，用於 polymorphic NuMT filter；亦見第 5 點的 median 版）。誤加會報 `autosomal-coverage is not a recognized option`，讓每個 case 的 MITO_FILTER 掛掉。Broad 現行 mito WDL 也不用它 → mito 過濾只靠 `--mitochondria-mode` + VariantFiltration blacklist mask。
35. **WhatsHap phasing（`--run_phasing`，NCKUH 專用）**：ensemble 是雙樣本(_DV/_HC)，whatshap 需 `--ignore-read-groups --sample <id>_HC`（phase HaplotypeCaller 欄，它 local assembly 最會把 compound 拆成相鄰兩筆）。biocontainer 只含 whatshap，故切 contig/合併/索引用 bcftools 容器、phase 用 whatshap 容器（per-contig scatter）。非破壞性（只加 PS），DRAGEN 自帶 PS 不走這條。
36. **PON 就位建議用 `install_pon.sh`**（verify → 備份舊版 → mv 新版就位 → rollback）；PON samplesheet 用 `subsample_pon.py`（依 run 日期取最近 + 男女均衡 + 依 sample 去重，建議 ~100–150 個同 assay 樣本）；PON 建置監控用 `monitor_pon.sh`。
