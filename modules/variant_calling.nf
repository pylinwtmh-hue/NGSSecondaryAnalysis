/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Variant Calling Module
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

// Lane 1: DeepVariant（GPU）
// 深度學習 SNP/Indel caller，對 WGS/WES 均有高 sensitivity
// 實測：RAM 27.6GB，CPU 3111%（~31 cores）
// 記憶體由 config withName 個別設定為 32GB
process PARABRICKS_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   fasta/fai/dict   - hg38 reference + index
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict

    // OUTPUT:
    //   vcf - DeepVariant genotyped VCF（壓縮 + tabix index）
    output:
    // 輸出標準 VCF，以便後續與 HaplotypeCaller 進行 Ensemble Calling (bcftools isec)
    tuple val(meta), path("*.deepvariant.vcf"), emit: vcf

    script:
    // 動態判斷 WES 專用參數
    def interval_arg  = params.seq_type == "WES" ? "--interval-file ${params.wes_targets}" : ""
    def model_arg     = params.seq_type == "WES" ? "--use-wes-model" : ""

    def use_lock      = params.use_gpu_lock ?: false
    def lock_script   = "/datalake_Intermediate/pipeline/pipeline_code/gpu_lock.sh"
    def unlock_script = "/datalake_Intermediate/pipeline/pipeline_code/gpu_unlock.sh"
    def num_gpus      = params.dv_num_gpus ?: 1

    def lock_block = use_lock ? """
    eval \$(bash ${lock_script} ${num_gpus})
    echo "[deepvariant] ${meta.id} 取得 GPU \${MY_GPUS}"
    trap "bash ${unlock_script} \${MY_GPUS}; echo \'[deepvariant] ${meta.id} 釋放 GPU \${MY_GPUS}\'" EXIT
    """ : """
    echo "[deepvariant] ${meta.id} 使用預設 GPU（非 DGX 環境）"
    """

    """
    echo "Waiting 30s for GPU memory cleanup..."
    sleep 30

    ${lock_block}

    pbrun deepvariant \
        --ref ${fasta} \
        --in-bam ${bam} \
        --out-variants ${meta.id}.deepvariant.vcf \
        ${interval_arg} \
        ${model_arg} \
        --num-gpus ${num_gpus}

    # trap 會在此自動觸發，釋放 GPU lock
    """
}

// Lane 2a: HaplotypeCaller（GPU）
// GATK HaplotypeCaller 的 GPU 加速版
// 實測：RAM 20.7GB，CPU 1435%（~14 cores）
// 記憶體由 config withName 個別設定為 24GB
process PARABRICKS_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]（recal.txt 用於 BQSR on-the-fly）
    //   fasta/fai/dict   - hg38 reference + index
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict

    // OUTPUT:
    //   vcf - HaplotypeCaller genotyped VCF（壓縮 + tabix index）
    output:
    tuple val(meta), path("*.haplotypecaller.vcf"), emit: vcf

    script:
    // 動態判斷：如果是 WES，自動加入 Capture Kit BED 檔限制運算範圍
    def interval_arg  = params.seq_type == "WES" ? "--interval-file ${params.wes_targets}" : ""
    // 抓取 config 裡的 htvc_args，若無定義則留空
    def custom_args   = params.htvc_args ?: ""

    def use_lock      = params.use_gpu_lock ?: false
    def lock_script   = "/datalake_Intermediate/pipeline/pipeline_code/gpu_lock.sh"
    def unlock_script = "/datalake_Intermediate/pipeline/pipeline_code/gpu_unlock.sh"
    def num_gpus      = params.htvc_num_gpus ?: 1

    def lock_block = use_lock ? """
    eval \$(bash ${lock_script} ${num_gpus})
    echo "[haplotypecaller] ${meta.id} 取得 GPU \${MY_GPUS}"
    trap "bash ${unlock_script} \${MY_GPUS}; echo \'[haplotypecaller] ${meta.id} 釋放 GPU \${MY_GPUS}\'" EXIT
    """ : """
    echo "[haplotypecaller] ${meta.id} 使用預設 GPU（非 DGX 環境）"
    """

    """
    echo "Waiting 30s for GPU memory cleanup..."
    sleep 30

    ${lock_block}

    pbrun haplotypecaller \
        --ref ${fasta} \
        --in-bam ${bam} \
        --in-recal-file ${recal_table} \
        --out-variants ${meta.id}.haplotypecaller.vcf \
        -G StandardAnnotation \
        -G StandardHCAnnotation \
        ${interval_arg} \
        --num-gpus ${num_gpus} \
        ${custom_args}

    # trap 會在此自動觸發，釋放 GPU lock
    """
}

// Lane 2a（備用）: HaplotypeCaller CPU
// 無 GPU 時的替代方案，結果與 GPU 版相同
process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   fasta/fai/dict   - hg38 reference + index
    //   dbsnp            - dbSNP vcf + index (用於標註 rsID)
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict
    path dbsnp
    path dbsnp_idx

    // OUTPUT:
    //   vcf - HaplotypeCaller genotyped VCF（壓縮 + tabix index）
    output:
    tuple val(meta), path("*.haplotypecaller.vcf.gz"), path("*.haplotypecaller.vcf.gz.tbi"), emit: vcf

    script:
    // 1. 動態判斷：如果是 WES，加入 -L 與 Capture Kit BED 檔
    def interval_arg = params.seq_type == "WES" ? "-L ${params.wes_targets}" : ""
    
    // 2. 記憶體安全閥：總記憶體扣除 4GB，保留給作業系統與 C++ PairHMM 使用
    def avail_mem = task.memory ? (task.memory.toGiga() - 4) : 12

    """
# GATK4 已移除 HaplotypeCaller 的 on-the-fly --BQSR（broadinstitute/gatk#6041）。
    # 改為先 ApplyBQSR 產生 recalibrated BAM，再進 HaplotypeCaller，
    # 使 CPU 備援結果與 GPU 主線（Parabricks --in-recal-file）一致。
    gatk --java-options "-Xmx${avail_mem}g" ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${meta.id}.recal.bam

    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
        -R ${fasta} \\
        -I ${meta.id}.recal.bam \\
        -O ${meta.id}.haplotypecaller.vcf.gz \\
        -D ${dbsnp} \\
        -G StandardAnnotation \\
        -G StandardHCAnnotation \\
        ${interval_arg} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --max-reads-per-alignment-start 0

    rm -f ${meta.id}.recal.bam ${meta.id}.recal.bai
    """
}

// Lane 2b: VQSR SNP（WGS only）
// 機器學習方法對 SNP 進行品質分層過濾
// WES 因變異數不足，跳過 VQSR，直接進 Ensemble
// 實測：RAM 9.4GB，CPU 123%（~1-2 cores）
// 記憶體由 config withName 個別設定為 12GB
process GATK_VQSR_SNP {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    // INPUT:
    //   vcf            - HaplotypeCaller 輸出的原始 VCF
    //   fasta/fai/dict - hg38 reference + index
    //   resources      - SNP 機器學習訓練所需的資料庫 (宣告於此以確保 Singularity 正確掛載)
    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fasta_fai
    path fasta_dict
    path hapmap
    path hapmap_idx
    path omni
    path omni_idx
    path known_snps
    path known_snps_idx
    path dbsnp
    path dbsnp_idx

    // OUTPUT:
    //   vcf      - SNP VQSR 過濾後的 VCF
    //   recal    - SNP recalibration model
    //   tranches - SNP sensitivity tranches
    output:
    tuple val(meta), path("*.vqsr_snp.vcf.gz"), path("*.vqsr_snp.vcf.gz.tbi"), emit: vcf
    path "*.snp.recal",    emit: recal
    path "*.snp.tranches", emit: tranches

    script:
    def prefix = "${meta.id}"
    
    // 記憶體安全閥：保留 2GB 給系統與 Java 虛擬機外的開銷
    def avail_mem = task.memory ? (task.memory.toGiga() - 2) : 14

    """
    # 步驟一：建立 SNP 高斯混合模型 (VariantRecalibrator)
    gatk --java-options "-Xmx${avail_mem}g" VariantRecalibrator \
        -R ${fasta} \
        -V ${vcf} \
        -O ${prefix}.snp.recal \
        --tranches-file ${prefix}.snp.tranches \
        --trust-all-polymorphic \
        -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
        -mode SNP \
        --max-gaussians 8 \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${known_snps} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp}

    # 步驟二：套用模型並執行過濾 (ApplyVQSR)
    gatk --java-options "-Xmx${avail_mem}g" ApplyVQSR \\
        -R ${fasta} \\
        -V ${vcf} \\
        -O ${prefix}.vqsr_snp.vcf.gz \\
        --recal-file ${prefix}.snp.recal \\
        --tranches-file ${prefix}.snp.tranches \\
        --truth-sensitivity-filter-level 99.5 \\
        -mode SNP
    """
}

// Lane 2b: VQSR INDEL（WGS only）
// SNP VQSR 完成後，對同一 VCF 再做 INDEL VQSR
// 實測：RAM 6.5GB，CPU 132%（~1-2 cores）
// 記憶體由 config withName 個別設定為 10GB
process GATK_VQSR_INDEL {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    // INPUT:
    //   vcf            - VQSR_SNP 輸出的 VCF（已完成 SNP 過濾）
    //   fasta/fai/dict - hg38 reference + index
    //   resources      - 機器學習訓練所需的資料庫 (必須宣告在這裡才能掛載進 Container)
    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fasta_fai
    path fasta_dict
    path mills
    path mills_idx
    path axiom
    path axiom_idx
    path dbsnp
    path dbsnp_idx

    // OUTPUT:
    //   vcf - SNP + INDEL VQSR 雙重過濾後的最終 VCF
    output:
    tuple val(meta), path("*.vqsr_indel.vcf.gz"), path("*.vqsr_indel.vcf.gz.tbi"), emit: vcf

    script:
    def prefix = "${meta.id}"
    
    // 記憶體安全閥：保留 2GB 給系統與 Java 虛擬機外的開銷
    def avail_mem = task.memory ? (task.memory.toGiga() - 2) : 10

    """
    # 步驟一：建立 INDEL 高斯混合模型 (VariantRecalibrator)
    gatk --java-options "-Xmx${avail_mem}g" VariantRecalibrator \
        -R ${fasta} \
        -V ${vcf} \
        -O ${prefix}.indel.recal \
        --tranches-file ${prefix}.indel.tranches \
        --trust-all-polymorphic \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode INDEL \
        --max-gaussians 4 \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
        --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${axiom} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp}

    # 步驟二：套用模型並執行過濾 (ApplyVQSR)
    gatk --java-options "-Xmx${avail_mem}g" ApplyVQSR \
        -R ${fasta} \
        -V ${vcf} \
        -O ${prefix}.vqsr_indel.vcf.gz \
        --recal-file ${prefix}.indel.recal \
        --tranches-file ${prefix}.indel.tranches \
        --truth-sensitivity-filter-level 99.0 \
        -mode INDEL
    """
}