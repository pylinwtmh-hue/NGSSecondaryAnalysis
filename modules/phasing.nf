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
 *   - 拆成兩個 process、各用單一工具容器（符合本 pipeline 慣例）：
 *       WHATSHAP_PHASE → whatshap 容器（可直接由 biocontainer 轉 sif，無需自建）
 *       WHATSHAP_INDEX → 既有 bcftools 容器（tabix 建索引 + publish）
 */

// ─────────────────────────────────────────────────────────────
// WhatsHap phasing（whatshap 容器；只需 whatshap）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_PHASE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.ensemble.phased.vcf.gz")

    script:
    """
    # ensemble 為雙樣本(_DV/_HC)：--ignore-read-groups 時需指定單一 sample（選 _HC）。
    # whatshap 會輸出 bgzip 壓縮 VCF（副檔名 .gz）；index 交給下一個 process（bcftools 容器）。
    # --reference 開啟 re-alignment，對 indel phasing 較準（需 ${fasta}.fai）。
    whatshap phase \\
        --reference ${fasta} \\
        --ignore-read-groups \\
        --sample ${meta.id}_HC \\
        -o ${meta.id}.ensemble.phased.vcf.gz \\
        ${vcf} ${bam}
    """
}

// ─────────────────────────────────────────────────────────────
// 建 tabix 索引 + publish（既有 bcftools 容器）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_INDEX {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcf)

    output:
    tuple val(meta),
          path("${meta.id}.ensemble.phased.vcf.gz"),
          path("${meta.id}.ensemble.phased.vcf.gz.tbi"), emit: vcf

    script:
    """
    tabix -p vcf ${meta.id}.ensemble.phased.vcf.gz
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
    ch_in = ensemble_ch
        .join(bam_ch, by: 0)
        .map { meta, vcf, tbi, bam, bai, recal -> [meta, vcf, tbi, bam, bai] }

    WHATSHAP_PHASE(ch_in, fasta, fasta_fai)
    WHATSHAP_INDEX(WHATSHAP_PHASE.out)

    emit:
    vcf = WHATSHAP_INDEX.out.vcf   // tuple(meta, phased_vcf, tbi)
}
