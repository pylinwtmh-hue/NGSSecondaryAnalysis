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
 *   - 三個 process、各用單一工具容器（符合本 pipeline 慣例）：
 *       WHATSHAP_SUBSET → bcftools（切出單一 contig）
 *       WHATSHAP_PHASE  → whatshap（可直接由 biocontainer 轉 sif，無需自建）
 *       WHATSHAP_CONCAT → bcftools（合併各 contig + 建索引 + publish）
 */

// ─────────────────────────────────────────────────────────────
// 切出單一 contig（bcftools 容器）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_SUBSET {
    tag "${meta.id}:${contig}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi), val(contig)

    output:
    tuple val(meta), val(contig),
          path("${meta.id}.${contig}.sub.vcf.gz"),
          path("${meta.id}.${contig}.sub.vcf.gz.tbi")

    script:
    """
    # header 含全部 contig，故即使該 contig 無變異也只是輸出空 VCF（exit 0）。
    bcftools view -r ${contig} ${vcf} -Oz -o ${meta.id}.${contig}.sub.vcf.gz
    tabix -p vcf ${meta.id}.${contig}.sub.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 單一 contig phasing（whatshap 容器；只需 whatshap）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_PHASE {
    tag "${meta.id}:${contig}"
    label 'process_medium'

    input:
    tuple val(meta), val(contig), path(sub_vcf), path(sub_tbi), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.${contig}.phased.vcf.gz")

    script:
    """
    # 空 contig（header 有但無變異，如 chrM / 女性 chrY）直接輸出，避免 whatshap 對空檔報錯。
    # 用容器內建的 python（whatshap 依賴）判斷有無非表頭列，不需 bcftools。
    if python3 -c "import gzip,sys; sys.exit(0 if any(not l.startswith('#') for l in gzip.open('${sub_vcf}','rt')) else 1)"; then
        # ensemble 為雙樣本(_DV/_HC)：--ignore-read-groups 時需指定單一 sample（選 _HC）。
        # --reference 開啟 re-alignment，對 indel phasing 較準（需 ${fasta}.fai）。
        whatshap phase \\
            --reference ${fasta} \\
            --ignore-read-groups \\
            --sample ${meta.id}_HC \\
            -o ${meta.id}.${contig}.phased.vcf.gz \\
            ${sub_vcf} ${bam}
    else
        cp ${sub_vcf} ${meta.id}.${contig}.phased.vcf.gz
    fi
    """
}

// ─────────────────────────────────────────────────────────────
// 合併各 contig + 建索引 + publish（bcftools 容器）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_CONCAT {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    input:
    tuple val(meta), path(phased_vcfs)

    output:
    tuple val(meta),
          path("${meta.id}.ensemble.phased.vcf.gz"),
          path("${meta.id}.ensemble.phased.vcf.gz.tbi"), emit: vcf

    script:
    """
    # whatshap 輸出的分片未建索引；concat -a 需要索引，先各自建 tbi。
    for f in ${phased_vcfs}; do tabix -f -p vcf \$f; done
    # -a 允許任意順序/重疊；再 sort 保證輸出座標有序。
    bcftools concat -a ${phased_vcfs} -Oz -o concat.vcf.gz
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
    ch_contigs = Channel.fromList(params.phasing_contigs)

    // 每個樣本 × 每個 contig 切檔
    WHATSHAP_SUBSET(ensemble_ch.combine(ch_contigs))

    // 併回該樣本的 BAM 後 phasing
    ch_phase_in = WHATSHAP_SUBSET.out
        .combine(bam_ch, by: 0)
        .map { meta, contig, sv, st, bam, bai, recal -> [meta, contig, sv, st, bam, bai] }
    WHATSHAP_PHASE(ch_phase_in, fasta, fasta_fai)

    // 依 meta 收齊所有 contig（size = contig 數，因每個 contig 都會產出一份）後合併
    ch_concat_in = WHATSHAP_PHASE.out
        .groupTuple(by: 0, size: params.phasing_contigs.size())
    WHATSHAP_CONCAT(ch_concat_in)

    emit:
    vcf = WHATSHAP_CONCAT.out.vcf   // tuple(meta, phased_vcf, tbi)
}
