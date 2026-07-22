/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Alignment QC Module
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

process SAMTOOLS_STATS {
    tag "$meta.id"
    // 實測：RAM 9.5MB，CPU ~100%（單核）
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/03_alignment_qc", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   （recal.txt 在此不使用，但需維持 tuple 結構）
    input:
    tuple val(meta), path(bam), path(bai), path(recal)

    // OUTPUT:
    //   stats - samtools stats 文字報告（mapping rate, error rate 等）（供 MultiQC 使用）
    output:
    path "*.stats", emit: stats

    script:
    // 保留 1 個 thread 給主程式，剩下的給解壓縮
    def threads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    samtools stats -@ ${threads} ${bam} > ${meta.id}.stats
    """
}


process MOSDEPTH {
    tag "$meta.id"
    // 實測：RAM 1.9GB，CPU 239%（~2-3 cores）
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/03_alignment_qc", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   targets          - WES: capture BED；WGS: 'NO_FILE'（pipeline 自動判斷）
    input:
    tuple val(meta), path(bam), path(bai), path(recal)
    path targets
    path autosome_bed   // WGS 深度 QC 用的 autosome BED（chr1-22），由 main.nf 傳入
                        // WES 傳 NO_FILE（WES 用 capture BED 就已限制範圍）

    // OUTPUT:
    //   global_dist    - 全基因組深度分布（供 MultiQC 使用）
    //   summary        - 各染色體深度摘要（WES 的 _region 行 = target region 真實深度）
    //   修改輸出：移除了 per-base，新增了 thresholds 和 regions
    output:
    path "*.global.dist.txt",     emit: global_dist
    tuple val(meta), path("*.summary.txt"), emit: summary
    path "*.thresholds.bed.gz*",  emit: thresholds   // 包含 .gz 和 .csi 索引
    path "*.regions.bed.gz*",     emit: regions      // 這是由 --by 產生的區塊深度

    script:
    // WES：用 capture BED 做 --by，summary 的 _region 行才是真實 target 深度
    // WGS：用 autosome BED（chr1-22）做 --by，排除 chrM/chrX/chrY/unplaced
    //      這是國際標準做法，避免 chrM 高拷貝數和性別染色體深度差異污染 QC 指標
    def by_cmd = (targets.name != 'NO_FILE')         ? "--by ${targets}"         : "--by ${autosome_bed}"
    """
    mosdepth \
        --threads ${task.cpus} \
        ${by_cmd} \
        --no-per-base \
        --thresholds 1,10,15,20,30,50,100 \
        --flag 1796 \
        ${meta.id} \
        ${bam}
    """
}


// ──────────────────────────────────────────────────────────────
// PLOIDY_CHECK：從 mosdepth summary 推 sex/ploidy（sex 防呆 + aneuploidy 提示）。
//   warn-only：只印警示 + 出 QC 檔，不改 ploidy、不讓 pipeline 失敗。
//   輸出（對齊 DRAGEN *.ploidy.vcf.gz 風格）：
//     - <id>.ploidy.vcf.gz ：每條 contig 一列 FORMAT=DC:NDC，header 帶 estimated/declared 核型
//     - <id>.ploidy_qc.txt ：人可讀摘要 + WARNINGS
//   ploidy_check.py 以 staged path input 傳入（content-hash → 改 script 後 -resume 正確重跑）。
// ──────────────────────────────────────────────────────────────
process PLOIDY_CHECK {

    tag "${meta.id}"

    publishDir "${params.out_dir}/${meta.id}/03_alignment_qc", mode: 'copy'

    input:
    tuple val(meta), path(summary)
    path ploidy_py

    output:
    tuple val(meta),
          path("${meta.id}.ploidy.vcf.gz"),
          path("${meta.id}.ploidy_qc.txt"), emit: ploidy

    script:
    """
    python3 ${ploidy_py} \
        --summary ${summary} \
        --sample ${meta.id} \
        --declared-sex ${meta.sex} \
        --seq-type ${params.seq_type} \
        --out-vcf ${meta.id}.ploidy.vcf \
        --out-qc ${meta.id}.ploidy_qc.txt
    # 對齊 DRAGEN 的 *.ploidy.vcf.gz（用 bcftools 壓縮；tertiary_python 容器含 bcftools）
    bcftools view ${meta.id}.ploidy.vcf -Oz -o ${meta.id}.ploidy.vcf.gz
    rm -f ${meta.id}.ploidy.vcf
    """
}