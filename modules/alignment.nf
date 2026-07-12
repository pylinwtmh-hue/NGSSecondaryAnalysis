/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Alignment Module
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

process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    // 實測：RAM 114GB（30x WGS），CPU 2106%（~21 cores）
    // DGX：--low-memory 模式，VRAM 壓到 16GB 以下，單卡跑一個 process，maxForks=6
    // 記憶體由 config withName 個別設定
    label 'process_high'
    label 'process_gpu'

    publishDir "${params.out_dir}/${meta.id}/02_alignment", mode: 'copy'

    // INPUT:
    //   reads         - fastp 過濾後的 FASTQ pair [R1, R2]
    //   fasta         - hg38 reference genome
    //   known_sites   - BQSR 用的已知變異位點
    //                   四份：dbSNP + Mills indels + GATK known indels + 1000G SNPs
    input:
    tuple val(meta), path(reads)
    path fasta
    path fasta_fai
    path fasta_dict
    path bwa_amb
    path bwa_ann
    path bwa_bwt
    path bwa_pac
    path bwa_sa
    path known_sites
    path known_sites_tbi
    path known_sites_2
    path known_sites_2_tbi
    path known_sites_3
    path known_sites_3_tbi
    path known_sites_4
    path known_sites_4_tbi

    // OUTPUT:
    //   alignment_bundle - [BAM, BAI, BQSR recalibration table]（後續所有步驟的主要輸入）
    //   qc_metrics       - duplicate metrics（供 MultiQC 使用）
    //   parabricks_qc    - Parabricks 內建 QC metrics 目錄
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), path("*.recal.txt"), emit: alignment_bundle
    path "*.duplicate_metrics.txt",                                     emit: qc_metrics
    path "qc_metrics_dir",                                              emit: parabricks_qc

    script:
    // 1. 動態判斷是否為 WES，如果是，則加入 interval-file 參數
    // 註：${params.wes_targets} 是透過 config 定義的絕對路徑，且已被 Singularity bind
    def interval_arg  = params.seq_type == "WES" ? "--interval-file ${params.wes_targets}" : ""

    // 2. 組合標準的 GATK Read Group 字串
    def rg_string     = "\"@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:${meta.id}\""

    // 3. 接收環境變數中的低記憶體設定等參數
    def custom_args   = params.fq2bam_args ?: ""

    // 4. reads 現在可能是 [R1a, R2a, R1b, R2b, ...]（多 lane）
    // meta.lanes 是 lane list（如 ['L001', 'L002']），由 main.nf 在 groupTuple 時保留
    // 每兩個 reads 一組配對對應的 lane，產生 unique RG ID
    def reads_list = reads instanceof List ? reads : [reads]
    def lanes_list = meta.lanes instanceof List ? meta.lanes : [meta.lanes ?: 'L001']
    def in_fq_args = reads_list.collate(2).withIndex().collect { pair, idx ->
        def r1      = pair[0]
        def r2      = pair[1]
        def lane    = lanes_list[idx] ?: "L${String.format('%03d', idx + 1)}"
        def rg_id   = "${meta.id}.${lane}"
        def rg_str  = "\"@RG\\tID:${rg_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:${rg_id}\""
        "--in-fq ${r1} ${r2} ${rg_str}"
    }.join(" \\\n        ")

    // GPU lock 開關：只有 DGX profile 會設定 use_gpu_lock = true
    def use_lock      = params.use_gpu_lock ?: false
    def lock_script   = "/datalake_Intermediate/pipeline/pipeline_code/gpu_lock.sh"
    def unlock_script = "/datalake_Intermediate/pipeline/pipeline_code/gpu_unlock.sh"

    // 要使用的 GPU 數量：dgx_single 設 6，dgx 設 1（預設）
    def num_gpus      = params.fq2bam_num_gpus ?: 1

    // GPU lock 區塊：use_lock=true 時才插入，否則給空字串
    // 注意：lock/unlock 腳本在容器外的 host 上執行，不受 Singularity 影響
    def lock_block = use_lock ? """
    # -------------------------------------------------------
    # GPU Lock：搶 ${num_gpus} 張空閒的 V100
    # -------------------------------------------------------
    eval \$(bash ${lock_script} ${num_gpus})
    echo "[fq2bam] ${meta.id} 取得 GPU \${MY_GPUS}"
    trap "bash ${unlock_script} \${MY_GPUS}; echo '[fq2bam] ${meta.id} 釋放 GPU \${MY_GPUS}'" EXIT
    """ : """
    # -------------------------------------------------------
    # 非 DGX 環境：使用 config 指定的預設 GPU
    # -------------------------------------------------------
    echo "[fq2bam] ${meta.id} 使用預設 GPU（非 DGX 環境）"
    """

    """
    mkdir -p qc_metrics_dir

    ${lock_block}

    pbrun fq2bam \
        --ref ${fasta} \
        ${in_fq_args} \
        --out-bam ${meta.id}.aligned.sorted.bam \
        --out-recal-file ${meta.id}.recal.txt \
        --out-duplicate-metrics ${meta.id}.duplicate_metrics.txt \
        --out-qc-metrics-dir qc_metrics_dir \
        --knownSites ${known_sites} \
        --knownSites ${known_sites_2} \
        --knownSites ${known_sites_3} \
        --knownSites ${known_sites_4} \
        ${interval_arg} \
        --bwa-options="-Y -K 100000000" \
        --fix-mate \
        --optical-duplicate-pixel-distance 2500 \
        --num-gpus ${num_gpus} \
        ${custom_args}

    # trap 會在此自動觸發，釋放 GPU lock
    """
    // --optical-duplicate-pixel-distance 2500：NovaSeq/NextSeq 適用
    // HiSeq2000 請改為 100
}
