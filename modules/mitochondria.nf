/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Mitochondria Module
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin (林伯昱)
 * 
 *  * This program is free software: you can redistribute it and/or modify
 *  * it under the terms of the GNU General Public License as published by
 *  * the Free Software Foundation, either version 3 of the License, or
 *  * (at your option) any later version.
 *  *
 *  * This program is distributed in the hope that it will be useful,
 *  * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  * GNU General Public License for more details.
 *  *
 *  * You should have received a copy of the GNU General Public License
 *  * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  *
 *  * THIRD-PARTY TOOLS NOTICE:
 *  * This pipeline orchestrates third-party tools subject to their own licenses.
 *  * Users of main_research.nf must comply with:
 *  *   - Manta (Illumina): PolyForm Strict License 1.0.0 (non-commercial only)
 *  *   - ExpansionHunter (Illumina): PolyForm Strict License 1.0.0 (non-commercial only)
 *  * See README.md and LICENSE for details.
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results. This software shall not be held liable for any
 * direct, indirect, or consequential damages arising from
 * its use.
 * =========================================================
 */

// =========================================================
// Lane 5: Mitochondria variant calling
//
// 流程說明：
//   1. 從全基因組 BAM 抽取 chrM reads + unmapped mates
//   2. 分兩路重新 alignment：
//      - 正常版（chrM_only）：處理大部分變異
//      - Shifted 版（chrM_shifted）：處理環狀 DNA 邊界區域
//   3. 各跑一次 Mutect2 mitochondria mode
//   4. Shifted 結果用 LiftoverVcf 轉回正常座標
//   5. MergeVcfs 合併兩份結果
//   6. FilterMutectCalls 最終過濾
//
// 為什麼需要 Shifted Reference？
//   粒線體是環狀 DNA，hg38 的 chrM 從任意位置切開成線性，
//   導致邊界區域（開頭與結尾）的 reads 無法正確比對。
//   Shifted reference 把切點移到中間，讓兩個版本互補。
//
// 容器分工（每個 process 只能用一個容器）：
//   MITO_EXTRACT_READS    → samtools
//   MITO_BAM2FASTQ_*      → samtools
//   MITO_BWA_*            → bwa
//   MITO_SORT_INDEX       → samtools（alias 成 NORMAL / SHIFTED）
//   MITO_MUTECT2_*        → gatk
//   MITO_LIFTOVER         → gatk
//   MITO_MERGE            → gatk
//   MITO_FILTER           → gatk
// =========================================================

// Step 1: 抽取 chrM 相關 reads
// 實測：RAM 28MB，CPU ~100%（單核）
process MITO_EXTRACT_READS {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   alignment_bundle - 全基因組 [BAM, BAI, recal.txt]
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)

    // OUTPUT:
    //   reads - 抽取出的 chrM + unmapped mates BAM（後續兩路 alignment 共用此輸入）
    output:
    tuple val(meta), path("${meta.id}.chrM_candidates.bam"),
                     path("${meta.id}.chrM_candidates.bam.bai"), emit: reads

    script:
    // samtools view -b ${bam} chrM -o chrM.bam
    // # 抓取 unmapped (未對齊, -f 4) 且不是 secondary/supplementary (-F 264) 的 reads
    // samtools view -b -f 4 -F 264 ${bam} -o unmapped.bam
    // samtools merge -f merged.bam chrM.bam unmapped.bam
    // samtools sort -@ ${task.cpus} -o ${meta.id}.chrM_candidates.bam merged.bam
    // samtools index ${meta.id}.chrM_candidates.bam
    """
    samtools view -b -@ ${task.cpus} ${bam} chrM -o ${meta.id}.chrM_candidates.bam
    samtools index -@ ${task.cpus} ${meta.id}.chrM_candidates.bam
    """
}

// Step 2a-1: BAM 轉 FASTQ（正常版）
// 實測：RAM 1.1GB，CPU 623%（~6 cores）
process MITO_BAM2FASTQ_NORMAL {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   bam/bai - chrM candidate BAM
    input:
    tuple val(meta), path(bam), path(bai)

    // OUTPUT:
    //   reads - FASTQ pair，準備重新比對到 chrM_only reference
    output:
    tuple val(meta), path("${meta.id}.chrM_normal_R1.fastq.gz"),
                     path("${meta.id}.chrM_normal_R2.fastq.gz"), emit: reads

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} -o qname_sorted.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${meta.id}.chrM_normal_R1.fastq.gz \
        -2 ${meta.id}.chrM_normal_R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        -n qname_sorted.bam
    """
}

// Step 2a-2: BWA alignment 到正常版 chrM
// 實測：RAM 856MB，CPU 722%（~7 cores）
// 注意：BWA index 檔案需明確宣告在 input 裡，Nextflow 不會自動帶入
process MITO_BWA_NORMAL {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   r1/r2                - chrM FASTQ pair
    //   chrM_only_fasta/...  - 正常版 chrM reference + BWA index 全套檔案
    input:
    tuple val(meta), path(r1), path(r2)
    path chrM_only_fasta
    path chrM_only_fai
    path chrM_only_dict
    path chrM_only_amb
    path chrM_only_ann
    path chrM_only_bwt
    path chrM_only_pac
    path chrM_only_sa

    // OUTPUT:
    //   bam - unsorted BAM（待 MITO_SORT_INDEX_NORMAL 排序）
    output:
    tuple val(meta), path("${meta.id}.chrM_normal.unsorted.bam"), emit: bam

    script:
    """
    bwa mem \
        -M \
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1" \
        -t ${task.cpus} \
        ${chrM_only_fasta} \
        ${r1} ${r2} \
        -o ${meta.id}.chrM_normal.unsorted.bam
    """
}

// Step 2b-1: BAM 轉 FASTQ（Shifted 版）
// 實測：RAM 1.1GB，CPU 576%（~6 cores）
process MITO_BAM2FASTQ_SHIFTED {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   bam/bai - chrM candidate BAM（與 NORMAL 路共用同一輸入）
    input:
    tuple val(meta), path(bam), path(bai)

    // OUTPUT:
    //   reads - FASTQ pair，準備重新比對到 chrM_shifted reference
    output:
    tuple val(meta), path("${meta.id}.chrM_shifted_R1.fastq.gz"),
                     path("${meta.id}.chrM_shifted_R2.fastq.gz"), emit: reads

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} -o qname_sorted.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${meta.id}.chrM_shifted_R1.fastq.gz \
        -2 ${meta.id}.chrM_shifted_R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        -n qname_sorted.bam
    """
}

// Step 2b-2: BWA alignment 到 Shifted 版 chrM
// 實測：RAM 856MB，CPU 725%（~7 cores）
process MITO_BWA_SHIFTED {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   r1/r2                   - chrM FASTQ pair
    //   chrM_shifted_fasta/...  - shifted 版 chrM reference + BWA index 全套檔案
    input:
    tuple val(meta), path(r1), path(r2)
    path chrM_shifted_fasta
    path chrM_shifted_fai
    path chrM_shifted_dict
    path chrM_shifted_amb
    path chrM_shifted_ann
    path chrM_shifted_bwt
    path chrM_shifted_pac
    path chrM_shifted_sa

    // OUTPUT:
    //   bam - unsorted BAM（待 MITO_SORT_INDEX_SHIFTED 排序）
    output:
    tuple val(meta), path("${meta.id}.chrM_shifted.unsorted.bam"), emit: bam

    script:
    """
    bwa mem \
        -M \
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1" \
        -t ${task.cpus} \
        ${chrM_shifted_fasta} \
        ${r1} ${r2} \
        -o ${meta.id}.chrM_shifted.unsorted.bam
    """
}

// Step 2c: Sort + MarkDuplicates + Index（純 GATK/Picard 實作）
// 同時 alias 成 MITO_SORT_INDEX_NORMAL 和 MITO_SORT_INDEX_SHIFTED
// 使用 GATK 容器內建的 Picard 工具，全程不需要 samtools：
//   SortSam        → 座標排序
//   MarkDuplicates → 標記 PCR duplicate（高深度粒線體 reads 的必要步驟）
//   BuildBamIndex  → 產生 .bam.bai 索引（輸出格式符合後續 Mutect2 輸入要求）

process MITO_SORT_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.*.sorted.markdup.bam"),
                     path("${meta.id}.*.sorted.markdup.bai"), emit: bam

    script:
    def prefix   = bam.name.replace('.unsorted.bam', '')
    def avail_mem = task.memory ? (task.memory.toGiga() - 2) : 6
    """
    # 1. 座標排序
    # --TMP_DIR . 讓暫存檔寫到 Nextflow work 目錄（/scratch 或 /raid），
    # 避免寫到系統 /tmp 導致空間不足
    gatk --java-options "-Xmx${avail_mem}g" SortSam \
        -I ${bam} \
        -O ${prefix}.sorted.bam \
        --SORT_ORDER coordinate \
        --TMP_DIR .

    # 2. 標記重複讀段
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \
        -I ${prefix}.sorted.bam \
        -O ${prefix}.sorted.markdup.bam \
        -M ${prefix}.markdup.metrics.txt \
        --ASSUME_SORT_ORDER coordinate \
        --TMP_DIR .

    # 3. 建立 BAM index
    # BuildBamIndex 產生 .bai（與 BAM 同名，副檔名為 .bai 而非 .bam.bai）
    gatk --java-options "-Xmx2g" BuildBamIndex \
        -I ${prefix}.sorted.markdup.bam \
        --TMP_DIR .

    # 清理暫存
    rm ${prefix}.sorted.bam
    """
}

// Step 3a: Mutect2 - 正常版 chrM
// 實測：RAM 4.1GB，CPU 400%（~4 cores）
// 記憶體由 config withName 個別設定為 6GB
process MITO_MUTECT2_NORMAL {
    tag "$meta.id"
    label 'process_medium'

    // INPUT:
    //   bam/bai        - 正常版 chrM sorted BAM
    //   fasta/fai/dict - chrM_only reference + index
    //   blacklist      - 已知 artifact 位點 BED
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path blacklist

    // OUTPUT:
    //   vcf - [VCF, TBI, stats]（stats 供 MergeMutectStats 使用）
    output:
    tuple val(meta), path("${meta.id}.mito_normal.vcf.gz"),
                     path("${meta.id}.mito_normal.vcf.gz.tbi"),
                     path("${meta.id}.mito_normal.vcf.gz.stats"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \
        -R ${fasta} \
        -I ${bam} \
        -O ${meta.id}.mito_normal.vcf.gz \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
    """
}

// Step 3b: Mutect2 - Shifted 版 chrM
// 實測：RAM 4.4GB，CPU 412%（~4 cores）
process MITO_MUTECT2_SHIFTED {
    tag "$meta.id"
    label 'process_medium'

    // INPUT:
    //   bam/bai        - shifted 版 chrM sorted BAM
    //   fasta/fai/dict - chrM_shifted reference + index
    //   blacklist      - 已知 artifact 位點 BED
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path blacklist

    // OUTPUT:
    //   vcf - [VCF, TBI, stats]（shifted 座標，需 LiftoverVcf 轉換）
    output:
    tuple val(meta), path("${meta.id}.mito_shifted.vcf.gz"),
                     path("${meta.id}.mito_shifted.vcf.gz.tbi"),
                     path("${meta.id}.mito_shifted.vcf.gz.stats"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \
        -R ${fasta} \
        -I ${bam} \
        -O ${meta.id}.mito_shifted.vcf.gz \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
    """
}

// Step 4: LiftoverVcf
// 實測：RAM 343MB，CPU 238%（~2-3 cores）
process MITO_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   shifted_vcf/tbi/stats - Mutect2 shifted 版輸出
    //   chrM_only_fasta/...   - 目標座標系（正常版 chrM reference）
    //   shift_back_chain      - shifted → normal 的座標轉換 chain 檔
    input:
    tuple val(meta), path(shifted_vcf), path(shifted_tbi), path(shifted_stats)
    path chrM_only_fasta
    path chrM_only_fai
    path chrM_only_dict
    path shift_back_chain

    // OUTPUT:
    //   vcf - [lifted VCF, TBI, stats]（已轉回正常 chrM 座標）
    output:
    tuple val(meta), path("${meta.id}.mito_shifted_lifted.vcf.gz"),
                     path("${meta.id}.mito_shifted_lifted.vcf.gz.tbi"),
                     path(shifted_stats), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" LiftoverVcf \
        -I ${shifted_vcf} \
        -O ${meta.id}.mito_shifted_lifted.vcf.gz \
        --CHAIN ${shift_back_chain} \
        --REJECT rejected.vcf.gz \
        -R ${chrM_only_fasta}

    gatk IndexFeatureFile -I ${meta.id}.mito_shifted_lifted.vcf.gz
    """
}

// Step 5: MergeVcfs + MergeMutectStats
// 實測：RAM 320MB，CPU 244%（~2-3 cores）
process MITO_MERGE {
    tag "$meta.id"
    label 'process_low'

    // INPUT:
    //   normal_vcf    - [normal VCF, TBI, stats]
    //   lifted_vcf    - [lifted shifted VCF, TBI, stats]（已轉回正常座標）
    //   chrM_only_dict - sequence dictionary（MergeVcfs 需要）
    input:
    // main.nf 用 .join(by: 0) 確保同一樣本配對後，合併成單一 tuple 傳入
    // 避免多樣本非同步完成時跨樣本錯配的 bug
    tuple val(meta),
        path(normal_vcf), path(normal_tbi), path(normal_stats),
        path(lifted_vcf), path(lifted_tbi), path(lifted_stats)
    path chrM_only_dict

    // OUTPUT:
    //   vcf - [merged VCF, TBI, merged stats]（兩路合併後的完整 chrM variants）
    output:
    tuple val(meta), path("${meta.id}.mito_merged.vcf.gz"),
                     path("${meta.id}.mito_merged.vcf.gz.tbi"),
                     path("${meta.id}.mito_merged.stats"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" MergeVcfs \
        -I ${normal_vcf} \
        -I ${lifted_vcf} \
        -O ${meta.id}.mito_merged.vcf.gz \
        -D ${chrM_only_dict}

    gatk MergeMutectStats \
        -stats ${normal_stats} \
        -stats ${lifted_stats} \
        -O ${meta.id}.mito_merged.stats
    """
}

// Step 6: FilterMutectCalls + VariantFiltration（blacklist mask）
// 實測：RAM 358MB，CPU 292%（~3 cores）
// 注意：--blacklisted-sites 在 GATK 4.6 已移除，改用 VariantFiltration --mask
process MITO_FILTER {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/07_mitochondria", mode: 'copy'

    // INPUT:
    //   merged_vcf/tbi/stats - MITO_MERGE 輸出
    //   fasta/fai/dict       - chrM_only reference + index
    //   blacklist            - 已知 artifact 位點 BED
    //   blacklist_idx        - blacklist BED 的 GATK index（需用 IndexFeatureFile 建立）
    //   mosdepth_summary     - MOSDEPTH 的 summary.txt reading depth 4.6.2 不再需要了
    input:
    tuple val(meta), path(merged_vcf), path(merged_tbi), path(merged_stats), path(mosdepth_summary)
    path fasta
    path fasta_fai
    path fasta_dict
    path blacklist
    path blacklist_idx

    // OUTPUT:
    //   vcf - 最終過濾後的 chrM variant VCF（PASS = 可信變異）
    output:
    tuple val(meta), path("${meta.id}.mito.vcf.gz"),
                     path("${meta.id}.mito.vcf.gz.tbi"), emit: vcf

    script:
    def avail_mem = task.memory ? (task.memory.toGiga() - 1) : 4
    """
    # 從 mosdepth summary 取 autosome (chr1-22) 平均深度，作為 NuMT filter 的核基因體覆蓋度。
    # GATK mito best-practice 用中位數；mosdepth summary 提供平均值，為合理近似。
    # 取不到（例如部分 WES summary）則省略 --autosomal-coverage，退回原本行為。
    AUTO_COV=\$(awk -F'\\t' '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {sum+=\$4; n++} END{ if(n>0) printf "%d", (sum/n)+0.5 }' ${mosdepth_summary})
    AUTOCOV_ARG=""
    if [ -n "\$AUTO_COV" ] && [ "\$AUTO_COV" -gt 0 ] 2>/dev/null; then
        AUTOCOV_ARG="--autosomal-coverage \$AUTO_COV"
        echo "[MITO_FILTER] ${meta.id} autosomal-coverage=\$AUTO_COV" >&2
    else
        echo "[MITO_FILTER] ${meta.id} 無 autosomal coverage，略過 --autosomal-coverage" >&2
    fi

    # 執行過濾
    gatk --java-options "-Xmx${avail_mem}g" FilterMutectCalls \
        -R ${fasta} \
        -V ${merged_vcf} \
        -O ${meta.id}.mito_filtered_tmp.vcf.gz \
        --stats ${merged_stats} \
        \$AUTOCOV_ARG \
        --mitochondria-mode

    # 套用黑名單過濾 (NuMTs 同源熱區)
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \
        -R ${fasta} \
        -V ${meta.id}.mito_filtered_tmp.vcf.gz \
        -O ${meta.id}.mito.vcf.gz \
        --mask ${blacklist} \
        --mask-name "blacklisted_site"
        
    # 清理暫存
    rm ${meta.id}.mito_filtered_tmp.vcf.gz*
    """
}