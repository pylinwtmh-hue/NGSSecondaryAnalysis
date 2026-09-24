# NGS WGS/WES Germline Analysis Pipeline — 開發筆記 / Environment Setup

> 這份是**內部開發紀錄**（環境、部署、踩雷記錄、驗證過程）。
> 對外的使用說明在 `README.md`。
>
> 檔名原本是小寫的 `readme.md`，但 Windows 檔案系統不分大小寫，
> 跟 `README.md` 撞名會讓 `git pull` / `git clone` 出錯，
> 2026-08 改名為 `DEVELOPMENT_NOTES.md`。

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
#       （相鄰 cis del+ins，如 SUZ12 delAAAinsTT）。由 --run_phasing 控制（DGX 驗證後已改為預設開；--run_phasing false 可關），
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
    --run_phasing true \
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

> ⚠️ **DeepVariant 的數字包含 RefCall**：`bcftools stats` 不看 GT，DV 否決的候選
> （`FILTER=RefCall`、GT `./.` / `0/0`）也算在 SNPs / indels 裡。要看 DV 真正 call 了幾個：
> ```bash
> $BCFTOOLS view -H -i 'GT="alt"' ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | wc -l
> $BCFTOOLS view -H -f RefCall     ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | wc -l
> ```
> **2026-09 起 `ensemble.fixed` 不再包含 DV 的 RefCall**（`BCFTOOLS_ENSEMBLE` 在 merge 前
> 只留 `GT="alt"`，見「踩雷記錄 → SUZ12」），所以 Ensemble 的紀錄數會比之前少 ——
> **這是預期的，不是 regression**。

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

> ⚠️ **DeepVariant 的數字包含 RefCall**：`bcftools stats` 不看 GT，DV 否決的候選
> （`FILTER=RefCall`、GT `./.` / `0/0`）也算在 SNPs / indels 裡。要看 DV 真正 call 了幾個：
> ```bash
> $BCFTOOLS view -H -i 'GT="alt"' ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | wc -l
> $BCFTOOLS view -H -f RefCall     ${OUTDIR}/04_snv_indel/${SAMPLE}.deepvariant.vcf.gz | wc -l
> ```
> **2026-09 起 `ensemble.fixed` 不再包含 DV 的 RefCall**（`BCFTOOLS_ENSEMBLE` 在 merge 前
> 只留 `GT="alt"`，見「踩雷記錄 → SUZ12」），所以 Ensemble 的紀錄數會比之前少 ——
> **這是預期的，不是 regression**。

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
35. **WhatsHap phasing（`--run_phasing`，NCKUH 專用）** — ⚠️ **架構已改為「ensemble merge 前、各 caller 各自 phase + combine」，見 #40；以下為舊版 post-merge 做法，保留作歷史脈絡**：ensemble 是雙樣本(_DV/_HC)，whatshap 需 `--ignore-read-groups --sample <id>_HC`（phase HaplotypeCaller 欄，它 local assembly 最會把 compound 拆成相鄰兩筆）。biocontainer 只含 whatshap，故切 contig/合併/索引用 bcftools 容器、phase 用 whatshap 容器（per-contig scatter）。非破壞性（只加 PS），DRAGEN 自帶 PS 不走這條。
36. **PON 就位建議用 `install_pon.sh`**（verify → 備份舊版 → mv 新版就位 → rollback）；PON samplesheet 用 `subsample_pon.py`（依 run 日期取最近 + 男女均衡 + 依 sample 去重，建議 ~100–150 個同 assay 樣本）；PON 建置監控用 `monitor_pon.sh`。
37. **WhatsHap 對混合倍體染色體會報 `PloidyError: Inconsistent ploidy (2 and 1)`** — ⚠️ **已被 #40 取代：phasing 移到 `+fixploidy` 之前，原始 caller VCF 全基因體皆 diploid → 不再有此錯、也不需 sex-aware 切分。以下為歷史脈絡**：ensemble 經二級 `bcftools +fixploidy` 後，男性 chrX 為 PAR diploid + 非PAR haploid（混合倍體），whatshap 要求單一染色體倍體一致 → 解析 chrX 時崩潰（跑完所有體染色體後才爆，浪費數小時）。解法（sex-aware，見 `buildPhaseShards()`）：只把 **diploid 區段**送 whatshap、haploid 段 passthrough —— 體染色體全 phase；女性/unknown chrX 全長 phase；男性 chrX 只 phase PAR1(1-2781479)+PAR2(155701383-156030895)、非PAR/chrY passthrough；chrM passthrough。倍體切法與 `postprocessing.nf` 的 `hg38_ploidy.txt` 一致，**兩邊要改需一起改**。分片數因性別而異（男 26 / 女 25），故 groupTuple 不指定 size。
38. **在 whatshap 容器裡做 `python3 -c` 空檔判斷會靜默失敗**：原本 WHATSHAP_PHASE 用
    `if python3 -c "...有無變異..."; then whatshap; else cp; fi` 想跳過空 contig。但 whatshap
    biocontainer 內 `python3` 環境不一定可用/該 one-liner 可能出錯 → `if` 為假 → 每個分片都走
    `cp` passthrough → **輸出完全沒 phase（連 autosome 都沒 PS）卻不報錯**。症狀：`ensemble.phased.vcf.gz`
    所有 het 的 PS 都是 `.`、header 沒有 `##FORMAT=<ID=PS>`。解法：拿掉空檔判斷，diploid 分片直接
    跑 whatshap（phase 分片本就有變異；真空檔就讓它 fail loud）。診斷：`bcftools view -h ... | grep -iE 'whatshap|ID=PS'` 若空 = whatshap 從沒執行。
39. **WhatsHap phasing 會漏掉「非主要 contig」的變異（非破壞性破功）**：`buildPhaseShards()` 舊版只切
    `chr1-22 / X / Y / M`，其餘 `*_alt / *_random / chrUn_* / *_decoy / HLA-*` 上的變異被 `bcftools view -r`
    整批丟掉 —— 實測一個 WGS 樣本掉了 **128,827 個（≈2.2%）**。症狀：`bcftools view -H` 比對
    `ensemble.fixed`（phasing 前）與 `ensemble.phased`（phasing 後）的變異數不相等（後者較少）。解法：加一個
    「其餘所有 contig」的 passthrough 分片，用 `bcftools -t` 的**補集語法** `^chr1,…,chrM`（`-r` 不支援 `^`）；
    `WHATSHAP_SUBSET` 用 POSIX `case` 判斷 —— `^` 開頭走 `-t`（補集）、其餘走 `-r`（索引較快）。contig 層級
    的補集不會有邊界重疊，這些非主要 contig 原樣保留、不 phase（不加 PS）。驗收：變異數應相等（註：
    此「數量相等」針對 contig 涵蓋；開啟 compound 合成後最終數量會**略降**，見 #40）。
40. **compound 合成架構（現行；取代 #35/#37 的 post-merge 做法）**：把被拆開的相鄰 cis 變異（如 SUZ12
    `c.2168_2170delAAAinsTT`）在進三級 VEP 前合成單一 canonical MNV，讓 HGVS p. 正確
    （`p.Glu723_Thr724delinsAla`）。
    - **為何在 ensemble merge「之前」、各 caller 各自做**：`bcftools merge`(DV+HC) 會把兩 caller 對同一
      compound 的不同表示法變成 **multiallelic**，而 whatshap **跳過 multiallelic** → compound 拿不到
      phase、也無法合。故必須在「還是單一 caller、還是 biallelic」時 phase+combine。
    - **二級（NCKUH）**：`main.nf` 對 DV/HC 原始單樣本 VCF 各自跑 `PHASE_COMBINE`（`modules/phasing.nf`）
      = whatshap phase（單樣本用 `--ignore-read-groups`，per-contig scatter）→ `scripts/combine_phased.py`，
      **再**進 `BCFTOOLS_ENSEMBLE`。產出的 `ensemble.fixed` 直接帶 phase(PS/`|`) + 已合成 compound；三級
      `prepare_vcf` 照舊讀 `*.ensemble.fixed.vcf.gz`（名稱不變）。只 publish 這個 fixed，中間檔留 `work/`。
    - **三級（DRAGEN）**：`prepare_vcf_dragen.nf` 的 `COMBINE_DRAGEN` 在 norm/tag 前跑同一支
      `combine_phased.py`，用 DRAGEN **原生 PS**（不需 whatshap）。`params.combine_phased`（預設 true）開關。
      2026-09 起只拿 PASS（與全部 chrM）進 combine：合成紀錄沿用 anchor 的 FILTER、三級只收 PASS，舊版會讓
      PASS 變異被同叢較寬的 non-PASS 紀錄一起丟掉（VAL-10 發現，見三級 DEVELOPMENT_NOTES）。
    - **phasing 在 `+fixploidy` 之前** → 原始 VCF 皆 diploid → 無 PloidyError → 不需 sex-aware 切分
      （#37 作廢），只依 contig 分片（主要 contig phase、其餘 passthrough，#39 的 catch-all 仍在）；倍體由
      後面 `+fixploidy` 校正。
    - **`combine_phased.py`（二級/三級同一支，只用 Python 標準庫，兩 repo 各存一份需同步）**：依「參考足跡」
      重疊或 cis gap≤`combine_max_gap`(預設 2，對齊 DRAGEN) 叢集，做局部單體重建成 MNV；只合 het cis / hom，
      trans 不重疊不合。合前先 ref-free trim 去 padding（否則 DV 的 `GAAA>GAA`+`A>T` 會假性重疊、漏掉 SNV；
      trim 後正確重建 `GAAA>GAT`）。實測 SUZ12：HC→`GAAA>GTT`、DV→`GAAA>GAT`，兩 caller 重建後常一致而塌成
      單筆（COMBINED tag 記合併筆數）。單元測試：`python3 scripts/test_combine_phased.py`（13 例）。
    - **注意**：combine 會**降低**變異數（compound 多筆→一筆），故 ensemble 變異數比未 phase 時略少（設計如此，
      非漏變異）。驗證看特定位點（SUZ12）+ combine 的 stderr `in/out/merged_clusters`，不要用「總數相等」。
41. **whatshap 把 DV 的 `FORMAT/AD` header 重新宣告成 `Number=.` → 與 HC 的 `Number=R` 不一致 → norm/merge 把 AD 弄壞**：
    原始 DeepVariant 與 HaplotypeCaller 的 AD **都是** `Number=R`；但 **whatshap 輸出 phased VCF 時，把
    DeepVariant 的 AD 重新宣告成 `Number=.`（Description 變 "Observed allele depths"），HaplotypeCaller 那邊
    維持 `Number=R`**（實測 NA12878 的 `*.DV.phased.combined` = `Number=.`、`*.HC.phased.combined` = `Number=R`）。
    於是兩個 `*.phased.combined.vcf.gz`（merge 的輸入）AD header 不一致，`BCFTOOLS_ENSEMBLE` 合併時 bcftools
    警告 `Trying to combine "AD" tag definitions of different lengths`。只要有一邊非 `Number=R`，bcftools 就
    無法把 AD 依 allele 正確處理，造成兩種當機：
    - **拆 multiallelic 時非 R 的 AD 不被 re-size** → biallelic 卻帶多個 AD 值：NA12878 `chr1:111241360`
      「2 alleles 卻 3 個 AD」→ `bcftools merge` 直接 `Incorrect number of FORMAT/AD values ... cannot merge`(255)。
    - **直接合併時 AD 沒依新 ALT union 補齊** → VAL-55 `chr1:83829`（3 ALT）HC 的 AD 只 3 值卻需 4 →
      三級 `bcftools norm -m -any` 報 `wrong number of fields in FMT/AD, expected 8, found 6`(255)，三級全掛。
    開 phasing 會**加重**（combine 產生的 MNV/`1|2` 逼出更多 multiallelic 合併），但**根因是 AD header
    Number 不一致**，raw DV+HC 也會中。解法（`modules/postprocessing.nf`，順序很重要）：
    1. **先把兩 caller 的 AD header 強制成 `Number=R`**（`bcftools view -h | sed 's/…ID=AD,Number=[^,]*,/…Number=R,/'`
       → `bcftools reheader -h`）。AD 本就是 per-allele，強制 R 正確、非 hack。
    2. 再各自 `bcftools norm -m -any` 拆 biallelic（此時 AD 會被正確 re-size）。
    3. 再 `bcftools merge --merge all`（標準 biallelic→multiallelic 聯集，Number=R/A/G 正確；缺的 allele 補 `.`）。
    4. 發布前 **preflight** `bcftools norm -m -any … -Ou -o /dev/null`，壞掉就 **fail loud**、不把壞檔丟三級。
    **不要**用 `norm --force`（那是丟棄壞 tag，會靜默掉 AD/VAF）。若之後 `PL`(Number=G) 等其他 tag 也報同類錯，
    比照把該 tag 的 header 補成正確 Number。
42. **`combine_phased.py` 合成紀錄只輸出 `GT:PS` → `AD`/`DP`/`VAF` 全丟成 `.`（三級 DRAGEN「AD 消失」bug）**：
    早期版本把重建出的 MNV 只寫 `GT:PS`，深度/等位分數欄位（`AD`/`DP`/`VAF`/`GQ`/`PL`）一律不帶。三級
    `add_dragen_tag.py` 用 `variant.format("AD")` 依名字讀 → 讀到空 → `AD_DRAGEN=.`；再經 `bcftools norm -m -any`
    拆開後兩筆 biallelic 都 `AD=. DP=. VAF=.`（同事回報 VAL-58 `chr17:80260571`；實測 VAL-10 共 **145,428** 筆
    phased 紀錄 AD 不見）。**根因不是 `bcftools norm`、也不是 `add_dragen_tag`，是 combine 設計把深度丟了**。
    解法（`scripts/combine_phased.py`，二級三級同一支）：
    - 合成的 biallelic MNV **繼承 anchor（叢集內足跡最寬的 biallelic diploid 顆）整組 FORMAT**，只覆寫 `GT`（重建的
      phased GT）與 `PS`，`QUAL`/`FILTER` 也沿用 anchor。anchor 是 biallelic diploid，其 `AD`(R)/`PL`(G) 元素數
      與合成後 biallelic 一致 → 長度正確、不會再壞。符合 bug report §4：保留**原始 locus** 的 VAF/AD，不用 2 元 AD
      重算成誤導的 `1.0`。
    - **ploidy-aware 重建**：diploid 叢集重建兩條單體（phased `0|1`/`1|1`…）；**全 haploid 非 chrM** 叢集
      （男性 non-PAR `chrX`/`chrY`）重建單套 → hemizygous `GT=1`（不帶 PS），一樣繼承 anchor 的 AD。
    - 四種情況**不重建、原封通過**（來源紀錄的 AD 原樣保住，交下游 `norm -m -any` 拆）：(a) 重建成 2 個 ALT
      （`1|2`，含原生 multiallelic `1/2`，正是 reporter 的 `CCGGCGG→CCGG,C AD=0,28,20`）；(b) 找不到 biallelic
      anchor；(c) 混 ploidy（haploid+diploid 同叢）；(d) `chrM` haploid 叢集（多拷貝異質性，不宜當單一分子合）。
    - stderr 由 `merged_clusters`（含 `haploid=`）加報 `passthrough_clusters`。回歸測試：`test_combine_phased.py`
      （15 例：合成保留 AD/DP/AF、`1|2` 退回、haploid 合成 hemizygous、混 ploidy 退回、chrM 退回、孤立逐字通過）。
    - **NCKUH combine 跑在 `+fixploidy` 之前 = 一律 diploid**，故 haploid 路實務上只在三級 DRAGEN 觸發；男性性染色體
      compound 在二級是「diploid 合 → 之後 fixploidy 轉 haploid」，本來就有處理,三級這條是補上它原本缺的。
43. **`combine_phased.py` 讓「沒有 ALT」的紀錄參與叢集 → DV 否決的候選被挑成 anchor（2026-09）**：
    叢集規則「足跡重疊必合」與 anchor「足跡最寬的 biallelic」都不看 GT；DeepVariant 否決的較寬候選
    （`FILTER=RefCall`、`./.`，常是缺失）蓋住真的 call（如 SNV）時，合成紀錄帶著被否決候選的
    `QUAL`/`FILTER=RefCall`/`GQ`/`DP`/`AD`/`VAF`/`PL`，POS 被撐到候選的起點（和 HC 對不上 → merge 合不起來 →
    三級 norm 後同一變異拆成 CALLERS=DV + CALLERS=HC 兩列），未 phase 的 het 被寫成 `0|1`＋假 PS，
    還可能把兩顆不相干的 call 串成一個 MNV。合成紀錄 GT 有 ALT，ensemble 的 DV `GT="alt"` 擋不掉。
    解法：`process()` 讀檔時 GT 沒有 ALT（`./.`、`0/0`、單套 `0`/`.`）的紀錄**不進叢集、原行輸出**。
    stderr 行尾加 `nocall_passthrough=`（DV ≈ RefCall 數，HC 應為 0）。回歸測試：`test_combine_phased.py`
    18 例（新增 `test_nocall_wider_candidate_not_anchor`、`test_nocall_does_not_bridge`、
    `test_haploid_nocall_passthrough`；三個在舊版都會失敗）。三級同一支（md5 一致）。實例見下方踩雷記錄。
44. **`combine_phased.py` 合成結果沒有最小化 → 同一變異在 DV、HC 寫法不同、merge 合不起來（2026-09）**：
    叢集範圍從第一個成分的 POS 起算，成分若是 indel 就帶著它的前導鹼基；DV、HC 拆成分的方式不同 →
    前導鹼基長度不同（例：HC `chr2:130206583 ACTT>AACC` vs DV `130206584 CTT>ACC`）→ POS 不同 →
    三級 norm 修剪後才一樣，報告裡拆成 DV 一列 + HC 一列（VAL55 修完第 43 條後還剩 820 個）。
    解法：`flush()` 在輸出前對重建結果跑 `trim_alleles()`（右修剪→左修剪，至少各留 1 個鹼基，所以
    純缺失 `GAAA>G` 的 anchor 會保留）。SUZ12 因此寫成 `31998951 AAA>TT`（原本 `31998950 GAAA>GTT`，
    與三級 norm 後相同）。回歸測試 `test_merged_output_minimised`（舊版會失敗）；
    `test_merged_keeps_format` 的預期改為最小化後的寫法。共 19 例。
45. **男性單倍體區的 het 被 `+fixploidy` 依 allele 順序截斷（2026-09）**：`+fixploidy` 只留第一個 allele，
    `0/1`、`0|1` → `0`（消失）、`1|0` → `1`（hemizygous）；開 phasing 後結果取決於 whatshap 任意定的方向。
    男性 het 多半是比對假象，但也可能是 47,XXY 或體細胞嵌合（X-linked 顯性、男性通常致死疾病的男性病人；
    PCDH19 只有嵌合男性發病），所以不能藏、也不能看方向。解法：`BCFTOOLS_ENSEMBLE` 對**男性**在
    `+fixploidy` 前 pipe 過 `scripts/haploid_het.awk`（staged input；讀同一份 ploidy 檔）：
    chrX ploidy=1 區間的 het `0/k` → `k/k` 並加 `INFO/HAPLOID_HET=<DV,HC>`（三級輸出 `HAPLOID_HET` 欄，
    AD/VAF 不動）；chrY 的 het → `./.`（不進報告）；PAR、體染色體、chrM、女性不動。只用 POSIX awk
    （bcftools 容器是 busybox 基底，沒有 Python）；`set -o pipefail` 的 subshell 裡跑，失敗會中止。
    測試：`test_haploid_het.py`（mawk 與 busybox awk 都過；`AWK="busybox awk"` 可切換）。
    重跑二級 `-resume` 會從 `COMBINE_PHASED`（第 44 條）與 `BCFTOOLS_ENSEMBLE`（新 input）往後重跑。
    **VAL55 重跑確認**：`chrX_het_to_alt_DV=2243 chrX_het_to_alt_HC=2699 chrY_het_to_missing_DV=2395
    chrY_het_to_missing_HC=10066`；帶標記的 4,326 筆全在 chrX；三級 NONE 10,159 筆全在 chrY；第 44 條的
    「拆兩列」820 → 2（剩下是重複序列 indel 左對齊差異；要消除需在 merge 前 `norm -f`，未做）。
    **chrM 不處理**（實驗室決定）：ensemble 裡 chrM 的 het 仍會被截斷，粒線體以 `07_mitochondria` 為準。
46. **`combine_phased.py` 把 phase 未知的重疊 het 當成同一條單體合成（2026-09）**：「足跡重疊一律合」＋
    未 phase 的 GT 在 `reconstruct()` 裡依位置都落在同一條單體 → 等於假設 cis。缺失範圍內的 SNV、包在大缺失裡
    的小缺失被吃掉（報告裡消失）；兩個重疊缺失被併成更長的缺失、錨在缺失中間的插入被截斷（寫出沒人 call 的
    allele）；還給假的 `0|1`+PS。phased 也會：同一條單體上互相矛盾的 call、同一鹼基的 SNV+插入（舊版字串排序先
    套插入，SNV 被吃）。**VAL-10（DRAGEN 女性 WGS）實測**：104,277 個會出報告的合成中 9,299 個是 phase 未知；
    **7,927 個 PASS allele 消失**、1,435 叢寫出沒人 call 的 allele。修法：(1) ≥2 顆 het 不在同一 phase set →
    不合（`_phase_unknown()`；phased 無 PS = 同一隱含 set）；(2) `build_hap()` 重疊只准落在沒改變的錨定鹼基，
    否則 `RebuildConflict` → 整叢原封通過（`*`/symbolic ALT 同）；(3) 同一 POS 先套 SNV/MNV 再套 indel
    （`_edit_order()`）；(4) 判斷抽成 `plan_cluster()`，stderr 加 `phase_unknown=` `overlap_conflict=`。
    unphased 但乾淨的 del+ins（1,096 叢）也不合：沒有 PS 就無法確定 cis；SUZ12 有 whatshap PS 照合。
    測試 19 → 26 個（新 7 個舊版全失敗）；40 組隨機資料 2,707 個合成沒有 allele 被吃掉或截斷。
    二級下次重跑生效（staged input，`-resume` 從 `COMBINE_PHASED` 往後）：看 stderr 兩個新數字與三級拆兩列數。
    同批三級修正（見三級 DEVELOPMENT_NOTES）：`COMBINE_DRAGEN` 只拿 PASS；`ADD_DRAGEN_TAG` 丟掉拆多等位後
    樣本沒帶的 allele（GT `0/0`，如 CYP21A2 `C>G,A 2/2` 拆出的 `C>G 0/0`）。

---

## 踩雷記錄

### ⚠️ SUZ12：ensemble 把 DV 否決的候選併進 HC 的 call（2026-09）

**症狀**：與 DRAGEN 比對 SUZ12 `c.2168_2170delinsTT` 時，三級報告多出一個錯誤的
`c.2170del`（`chr17:31998950 GA>G`，AD 10,14、VAF 0.583），正確的 delinsTT 反而顯示
`AD 10,0 / VAF 0`。DRAGEN 只有正確的那一個。

**ensemble 在這個位點實際長這樣**：

```
chr17 31998950 GAAA GAA,GTT RefCall COMBINED=2
  DV: ./.  DP 24  AD 10,14,.  VAF 0.583,.     ← ALT1 只有 DV 有資料（DV 否決的候選）
  HC: 0|2  DP 5   AD 3,.,2    PS 31998950     ← ALT2 只有 HC 有資料（HC 合成的 delinsTT）
chr17 31998952 A T RefCall   DV ./.（AD 11,14）  HC ./.
chr17 31998953 A T RefCall   DV ./.（AD 10,15）  HC ./.
```

- DV 看到 delinsTT 的三個片段（950 缺一個 A、952 A>T、953 A>T），**三個都判 RefCall**
- HC 的 combine_phased **正確**合成 `GAAA>GTT`（`0|1`、`COMBINED=2`）—— **phasing 沒有問題**
- `merge --merge all` 把同一 POS 的 DV `GAAA>GAA` 與 HC `GAAA>GTT` 併成一筆多等位，
  FILTER 取了 DV 的 `RefCall`
- 三級 norm 拆開後，DV 否決的 allele 變成「兩邊都沒 call」的獨立紀錄，又碰上三級
  `determine_callers()` 把「兩邊都沒 call」標成 HC 的 bug → 四筆全進 ACMG 表

**二級的修正**：`BCFTOOLS_ENSEMBLE` 的 DV arm 在 `norm -m -any` 之後加
`bcftools view -i 'GT="alt"'`，DV 否決的候選根本不進 merge。三級另外修了
`determine_callers()`（→ `NONE`）與 `get_ad()`（缺值保留 `.`），兩邊是彼此獨立的防線。
詳細推導、測試與重現方法見三級的 `DEVELOPMENT_NOTES.md`「SUZ12 幽靈變異」。

**實作細節（都實測過）**：
- **先 norm 再 filter**：DV 的多等位 `0/2` 拆開後是 `0/0`（丟）+ `0/1`（留）。
- **不用 pipe**：shell 是 `bash -ue`、沒有 `pipefail`，`norm | view` 在 norm 中途失敗時
  可能以 0 結束並寫出截斷檔 → 拆成 `norm -o norm_dv.bcf` 與 `view` 兩步。
  （`modules/cnv_sv.nf` 的 `bcftools view -f PASS | bcftools sort` 也有同樣的暴露，本次未動。）
- **`GT="alt"`** 丟掉 `./. 0/0 0|0 ./0` 以及半缺失 `./1 1/.`，與三級 `is_called()` 一致。
- **`--merge both`（預設）其實不會併這兩筆** —— delins 不是單純 indel，只有 `all` 會併。
  本 pipeline 用 `all`；加上 DV 過濾後，這種「否決候選併進真 call」已不會再發生。

**影響**：`ensemble.fixed` 紀錄數下降（不再含 DV RefCall），上面 Variant Count 已註明。
RefCall 仍保留在已發布的 `<id>.deepvariant.vcf.gz`；CNVkit 的 b-allele 讀那份，不受影響。

### ⚠️ VAL55 重跑：ensemble 仍有 28,050 筆 RefCall → combine_phased 的 anchor bug（2026-09）

**症狀**：加上 DV `GT="alt"` 過濾後重跑 VAL55（WGS，男性），`ensemble.fixed` 仍有 **28,050** 筆
`FILTER=RefCall`（header 確認 `GT="alt"` 那一步有跑）。**全部帶 `INFO/COMBINED`** → 都是 combine_phased 合成的。

**實際紀錄（節錄）**：

```
chr1 744865 CG>CA  2.8 RefCall COMBINED=2   DV 1|1 GQ 3 DP 12 AD 6,6 VAF 0.5 PL 0,2,4   HC ./.
chr1 602156 CA>CG,GG   RefCall;VQSRTrancheSNP99.90to100.00 COMBINED=2
     DV 1|1 GQ 6 AD 8,5,. VAF 0.38 PL 0,11,4        HC 2|0 AD 6,.,3
```

- 兩筆 DV 的 GT 都是 `1|1`（hom），但繼承來的 AD/VAF/PL 是 het 樣或「最可能是 REF」（`PL 0,…`）
  → 這些數字屬於被否決的候選，不是這個 call。
- 三級報告有 **23,023** 個變異被拆成「DV 一列 + HC 一列」（上限估計）。

**原因與修法**：見上方第 43 條。

**修正後重跑（`-resume`，從 `COMBINE_PHASED` 往後）**：ensemble 的 RefCall 28,050 → **0**；
combine stderr `nocall_passthrough` DV **865,130**（≈ DV 否決的候選數，約佔 DV 紀錄 16%）、HC **0**；
三級「拆兩列」23,023 → **820**；ADD_CALLERS_TAG DV+HC 89.9% / DV only 2.6% / HC only 7.4%。
剩下 820 個：抽查 5 例都是兩個 caller 對同一變異的寫法不同：combine 合成時叢集範圍從第一個成分的 POS 開始，成分若是 indel 就帶著它的前導鹼基，而合成結果沒有再最小化；DV、HC 拆成分的方式不同 → 前導鹼基長度不同 → POS 不同，merge 合不起來，三級 norm 修剪後才一樣。兩邊各自最小化後 5 例完全相同。→ **已修**：合成結果輸出前先最小化（上方第 44 條）。

**另外兩件同時確認的事**：
- **SUZ12 在這次重跑已不是當初出錯的情境**：同一個缺失成分 DP/AD/VAF 完全相同（24 / 10,14 / 0.583），
  但 DV 的 GQ 9→15、PL `0,8,34`→`15,0,24`，由 RefCall 翻成 het；HC 的 anchor 也由 DP 5 變 DP 25。
  原因是兩次用的 **Parabricks 版本不同**。這個位點 DV 本來就在邊緣（GQ 都很低），版本一換就翻。
  → 評鑑用的樣本應固定同一個 Parabricks 版本重跑，並把版本記進驗證紀錄。
- **男性 chrX/chrY：phasing + `+fixploidy` 的交互作用（2026-09 已處理，見上方第 45 條）**：
  phasing 把 chrX/chrY 也送 whatshap（原始 VCF 一律 diploid），之後 `+fixploidy` 把男性 non-PAR 的 GT
  **截成第一個 allele**（bcftools `plugins/fixploidy.c`：`0/1`、`0|1`→`0`，`1|0`→`1`）。沒開 phasing 時 het
  一律變 REF；開了之後 phase 方向決定結果 → 一部分男性 chrX「het」（多半是誤比對）在報告裡顯示成
  **hemizygous**。VAL55 實測（非合併、GT=1、有 PS，即原本 phased het 被截成 ALT）：DV 234、HC 403 筆；
  被截成 REF 的（GT=0、有 PS）DV 292、HC 875 筆。修正後重跑的三級 NONE 8,656 筆**全部**在 chrX（2,607）
  與 chrY（6,049）＝被截成 REF、報告裡看不到的男性 het。
  ⚠️ 不能一律藏掉：46,XY 男性遺傳來的變異（不論顯性或隱性）都是 hemizygous，不受影響；但男性出現的 het
  除了比對假象，也可能是 47,XXY 或**體細胞嵌合** —— X-linked dominant、男性通常致死的疾病（IKBKG、MECP2、
  CDKL5、PORCN、OFD1…）存活的男性病人常是嵌合；PCDH19 則是 hemizygous 男性通常不發病、嵌合男性才發病。
  → 決定（2026-09）：chrX 非 PAR 的 het 一律保留成 ALT 並加 `HAPLOID_HET` 標記給人工複核；chrY 的 het 不進報告。

## 未來進步方向（Roadmap；尚未實作，備忘）

討論過、但還沒落地的中長期方向。總目標：**只要有 FASTQ，性別 / ploidy / karyotype 都應該能自動推定**，
且**像 DRAGEN 那樣把「各有強項的多個 caller」依區域組合起來**，而不是單一 caller 全域打天下。

### 1. 從 FASTQ 全自動判定性別與 ploidy（免手動填 samplesheet 的 sex）
現況：二級 `sex` 由 samplesheet 手動填、三級沒有 sex 欄。目標從對齊後 BAM 直接推：
- **性別**：mosdepth 的 X:autosome 覆蓋比、Y 覆蓋、X 上雜合率 → XX / XY / 其他（mosdepth 已在 `03_alignment_qc` 跑，資料現成）。
- **karyotype / ploidy**：見第 3 點。推得的值回頭餵 `+fixploidy`、gCNV contig-ploidy、與第 2 點的 caller ploidy。
- DRAGEN 已內建 Ploidy Estimator（輸出 `*.ploidy.vcf.gz`），NCKUH 端要自建。

### 2. Ploidy-aware variant calling（目前兩個 caller 都一律 diploid，只靠事後 `+fixploidy`）
現況：`modules/variant_calling.nf` 的 DeepVariant / HaplotypeCaller **都沒有**傳 ploidy / 性別旗標，一律
diploid call，最後才 `bcftools +fixploidy` 修 GT。**事後修 GT 無法回復**用錯 ploidy 的誤判（男性 non-PAR
chrX/chrY 會出現假 het，尤其 DeepVariant）。方向：
- **HaplotypeCaller**：男性 non-PAR chrX/chrY 用 `-ploidy 1`（PAR 維持 2），interval 分片再合。
- **DeepVariant**：`--haploid-contigs chrX,chrY` + `--par-regions-bed`（Parabricks 4.x 應支援，需確認版本旗標）。
- 已知三倍體染色體：該染色體 `-ploidy 3` 重跑（GATK 可；DV 只 diploid/haploid；DRAGEN 存疑，非開源）。
- ⚠️ **與 `--run_phasing` 的交互**：目前 phasing/combine 跑在 `+fixploidy` 之前、假設「一律 diploid」；改成
  sex-aware calling 後 raw VCF 會有 haploid 區（combine 的 haploid 路已實作、會被用到），phasing 分片也要重設計。
  屬較大的一次改動，需獨立設計；先做第 1 點（可靠判性別）才有基礎。

### 3. Aneuploidy 自動偵測 + 提示（讓人知道要不要手動 `-ploidy N` 重跑）
- **read-depth**：per-chromosome 覆蓋 z-score（chr21 ≈ 1.5×→ 三體）；沿用 mosdepth / gCNV / CNVkit。
- **SNV BAF**：het SNV 的 B-allele fraction 正常一條帶 ~0.5，三體會分裂成 ~0.33 / 0.67 兩條帶（用到現在保住的 AD/VAF）。
- 兩訊號一致 → QC/報告標「chrN 疑似非整倍體，本染色體 SNV 基因型可能不準，考慮 `-ploidy N` 重跑」。
- 短期先做「sex 防呆」把性染色體 aneuploidy（XXY 等）先擋起來（見下）。

### 4. 多 caller 組合，善用各自強項（像 DRAGEN）
目前 NCKUH = DeepVariant + HaplotypeCaller ensemble（全域）。方向是納入「在特定區域更強」的 caller，
依**區域**取捨而非全域平均：
- **dark / 高同源 / 低複雜度區**：加對這些區域較強的 caller（或 graph-based / long-read-aware 方法）。
- 建立 region-aware 的信心度與合併策略（confidence by region），而非單純多數決。

### 近期落地步驟（近 → 遠）
- [x] **sex 防呆 + per-chromosome ploidy 提示**：`PLOIDY_CHECK`（`modules/alignment_qc.nf` +
      `scripts/ploidy_check.py`）從 mosdepth summary 推 sex/ploidy，與 samplesheet sex 比對，不符或疑似
      非整倍體 → **warn-only** + 出 `03_alignment_qc/*.ploidy.vcf.gz` ＋ `*.ploidy_qc.txt`。
      **NDC 對估計核型正規化（與 DRAGEN 一致：男 chrX ≈ 1.0）**，另存 `RATIO`＝相對體染色體的原始比
      （男 chrX ≈ 0.5，性別/劑量證據）；VCF header 用 `##estimatedSexKaryotype`/`##referenceSexKaryotype`
      對齊 DRAGEN。（het BAF 交叉驗證仍待做；WES 覆蓋較吵，建議搭 gCNV。）
- [x] **strand-bias 警示欄**：三級 `parse_vep_csq.py` 的 `STRAND_BIAS` 欄（FS/SOR；DRAGEN + NCKUH-HC 有，
      DeepVariant-only → `.` 人工複核）。
- [x] **三級接 DRAGEN `*.ploidy.vcf.gz`**：`PLOIDY_REPORT_DRAGEN` + `scripts/parse_dragen_ploidy.py`
      讀 DRAGEN 原生 ploidy.vcf（`##estimatedSexKaryotype` + 每 contig `NDC`）→ 產「與二級同一套」的
      `00_prepare/*.ploidy_qc.txt`（性別、sex_check、NDC、aneuploidy），warn-only。兩條路 NDC 語意已統一。
- [ ] （長期）第 2、4 點：sex-aware / ploidy-aware calling 與多 caller region-aware 組合。
