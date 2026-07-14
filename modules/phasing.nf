/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Phasing Module
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin (林伯昱)
 * Licensed under the GNU General Public License v3.0
 *
 * DISCLAIMER: This pipeline is provided "as is" without warranty of any
 * kind. Users are solely responsible for validating and interpreting all
 * results. See LICENSE for details.
 * =========================================================
 * modules/phasing.nf
 * ==================
 * 用 WhatsHap 對 NCKUH ensemble VCF 做 read-backed phasing，補上 PS phase set，
 * 讓三級能判斷「相鄰變異是否 in cis」，正確合併/註解 compound（如相鄰 del+ins）。
 *
 *   - 只用於 NCKUH（DV+HC ensemble）路徑；DRAGEN 自帶 PS，不走這裡。
 *   - 非破壞性：只加 PS + 把可 phase 的 GT 由 '/' 改為 '|'，不改變任何 variant 內容。
 *   - 由 params.run_phasing 開關，預設 false（在 DGX 驗證後再開）。
 *   - ensemble 是雙樣本(_DV/_HC)；WhatsHap 需單一樣本，故 --sample <id>_HC
 *     （HaplotypeCaller 有 local assembly，最可能把 compound 拆成相鄰兩筆）。
 *   - 依 contig 平行(scatter)以控制 WGS 執行時間；phase block 受 read 連通性限制，
 *     本來就是 local，per-contig 與全基因體結果等價（phase set 不跨 contig）。
 *
 * 容器：whatshap + bcftools（見 nextflow_main.config）。授權皆 MIT/GPL，可商用。
 */

// ─────────────────────────────────────────────────────────────
// 每個 contig 各自 phase（scatter）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_PHASE {
    tag "${meta.id}:${contig}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai), val(contig)
    path fasta
    path fasta_fai

    output:
    tuple val(meta),
          path("${meta.id}.${contig}.phased.vcf.gz"),
          path("${meta.id}.${contig}.phased.vcf.gz.tbi")

    script:
    def hc_sample = "${meta.id}_HC"
    """
    # 取出此 contig 的變異（WhatsHap 只需讀該區間的 reads）
    bcftools view -r ${contig} ${vcf} -Oz -o sub.vcf.gz
    tabix -p vcf sub.vcf.gz

    if [ \$(bcftools view -H sub.vcf.gz | head -c1 | wc -c) -eq 0 ]; then
        # 此 contig 沒有變異：直接輸出（避免 WhatsHap 對空檔報錯）
        cp sub.vcf.gz ${meta.id}.${contig}.phased.vcf.gz
    else
        # 雙樣本 ensemble：--ignore-read-groups 時必須指定單一 sample（選 _HC）
        # phase 失敗（該 contig 無可用 read 等）則保留未 phase 版本，不中斷整條 pipeline
        whatshap phase \\
            --reference ${fasta} \\
            --ignore-read-groups \\
            --sample ${hc_sample} \\
            -o ${meta.id}.${contig}.phased.vcf.gz \\
            sub.vcf.gz ${bam} 2> ${meta.id}.${contig}.whatshap.log \\
        || cp sub.vcf.gz ${meta.id}.${contig}.phased.vcf.gz
    fi
    tabix -f -p vcf ${meta.id}.${contig}.phased.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 收集各 contig 的 phased VCF，concat + sort 成單一檔
// ─────────────────────────────────────────────────────────────
process WHATSHAP_CONCAT {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta),
          path("${meta.id}.ensemble.phased.vcf.gz"),
          path("${meta.id}.ensemble.phased.vcf.gz.tbi"), emit: vcf

    script:
    """
    # -a 允許任意順序 / 重疊；再 sort 保證輸出座標有序
    bcftools concat -a ${vcfs} -Oz -o concat.vcf.gz
    bcftools sort concat.vcf.gz -Oz -o ${meta.id}.ensemble.phased.vcf.gz
    tabix -p vcf ${meta.id}.ensemble.phased.vcf.gz
    rm -f concat.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 組合 workflow：吃 ensemble VCF + alignment BAM，輸出 phased ensemble VCF
// ─────────────────────────────────────────────────────────────
workflow PHASE_ENSEMBLE {
    take:
    ensemble_ch   // tuple(meta, vcf, tbi)
    bam_ch        // tuple(meta, bam, bai, recal)  ← alignment_bundle
    fasta
    fasta_fai

    main:
    // 每個樣本 × 每個 contig 展開
    ch_contigs = Channel.fromList(params.phasing_contigs)
    ch_in = ensemble_ch
        .join(bam_ch, by: 0)
        .map { meta, vcf, tbi, bam, bai, recal -> [meta, vcf, tbi, bam, bai] }
        .combine(ch_contigs)   // → (meta, vcf, tbi, bam, bai, contig)

    WHATSHAP_PHASE(ch_in, fasta, fasta_fai)

    // 依 meta 收齊所有 contig（size = contig 數）後 concat
    ch_grouped = WHATSHAP_PHASE.out
        .groupTuple(by: 0, size: params.phasing_contigs.size())
    WHATSHAP_CONCAT(ch_grouped)

    emit:
    vcf = WHATSHAP_CONCAT.out.vcf   // tuple(meta, phased_vcf, tbi)
}
