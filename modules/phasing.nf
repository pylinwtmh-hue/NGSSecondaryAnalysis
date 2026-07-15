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
 *
 * ── 依性別的倍體切分（sex-aware，與二級 +fixploidy 一致）─────────────────────
 *   whatshap 要求「單一染色體倍體一致」，否則報 PloidyError(2 and 1)。二級
 *   `bcftools +fixploidy` 依性別設定 chrX/chrY/chrM 倍體，所以 phasing 必須只餵
 *   「diploid 區段」給 whatshap，其餘 passthrough（保留變異、不加 PS）：
 *     * 體染色體 chr1-22            → diploid → phase
 *     * 男性 chrX PAR1(1-2781479) + PAR2(155701383-156030895) → diploid → phase
 *     * 男性 chrX 非PAR、chrY       → haploid → passthrough
 *     * 女性/unknown chrX（全長）    → diploid → phase（fixploidy 對 unknown 視為 F）
 *     * chrM                        → haploid → passthrough
 *   座標與 postprocessing.nf 的 hg38_ploidy.txt 一致；要改請兩邊一起改。
 *
 *   分片各用單一工具容器（符合本 pipeline 慣例）：
 *     WHATSHAP_SUBSET → bcftools（切出該分片的區段，可含多段以逗號分隔）
 *     WHATSHAP_PHASE  → whatshap（type=phase 才 phase，type=pass 直接 cp）
 *     WHATSHAP_CONCAT → bcftools（合併所有分片 + 索引 + publish）
 */

// ─────────────────────────────────────────────────────────────
// 依性別產生 phasing 分片：每筆 [region, type, idx]
//   type: 'phase'（diploid，送 whatshap）/ 'pass'（haploid，直接 passthrough）
//   男性判定與 BCFTOOLS_ENSEMBLE 一致（M / MALE），其餘（F / unknown / 空）視為女性。
// ─────────────────────────────────────────────────────────────
def buildPhaseShards(sex) {
    def s = (sex ?: '').toString().toUpperCase()
    def male = (s == 'M' || s == 'MALE')
    def sh = []
    (1..22).each { sh << ["chr${it}".toString(), 'phase'] }
    if (male) {
        // PAR1 + PAR2 為 diploid → phase；非PAR + 尾段 + chrY 為 haploid → passthrough
        sh << ['chrX:1-2781479,chrX:155701383-156030895'.toString(), 'phase']
        sh << ['chrX:2781480-155701382,chrX:156030896-156040895'.toString(), 'pass']
        sh << ['chrY', 'pass']
    } else {
        sh << ['chrX', 'phase']
        sh << ['chrY', 'pass']
    }
    sh << ['chrM', 'pass']
    return sh.withIndex().collect { e, i -> [e[0], e[1], i] }
}

// ─────────────────────────────────────────────────────────────
// 切出該分片的區段（bcftools 容器）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_SUBSET {
    tag "${meta.id}:${region}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi), val(region), val(type), val(idx)

    output:
    tuple val(meta), val(idx), val(type),
          path("${meta.id}.s${idx}.sub.vcf.gz"),
          path("${meta.id}.s${idx}.sub.vcf.gz.tbi")

    script:
    """
    # region 可含多段（逗號分隔，如男性 chrX PAR1,PAR2）；bcftools -r 支援。
    bcftools view -r ${region} ${vcf} -Oz -o ${meta.id}.s${idx}.sub.vcf.gz
    tabix -p vcf ${meta.id}.s${idx}.sub.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 分片 phasing（whatshap 容器；type=phase 才 phase）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_PHASE {
    tag "${meta.id}:s${idx}(${type})"
    label 'process_medium'

    input:
    tuple val(meta), val(idx), val(type), path(sub_vcf), path(sub_tbi), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.s${idx}.phased.vcf.gz")

    script:
    if (type == 'phase')
        """
        # 有變異才 phase，空分片直接輸出（避免 whatshap 對空檔報錯）。
        # 用容器內建的 python（whatshap 依賴）判斷有無非表頭列，不需 bcftools。
        if python3 -c "import gzip,sys; sys.exit(0 if any(not l.startswith('#') for l in gzip.open('${sub_vcf}','rt')) else 1)"; then
            # --ignore-read-groups 時需指定單一 sample（選 _HC）；--reference 開 re-alignment 對 indel 較準。
            whatshap phase \\
                --reference ${fasta} \\
                --ignore-read-groups \\
                --sample ${meta.id}_HC \\
                -o ${meta.id}.s${idx}.phased.vcf.gz \\
                ${sub_vcf} ${bam}
        else
            cp ${sub_vcf} ${meta.id}.s${idx}.phased.vcf.gz
        fi
        """
    else
        """
        # haploid 區段（男性非PAR chrX / chrY / chrM）：passthrough，不 phase、不加 PS。
        cp ${sub_vcf} ${meta.id}.s${idx}.phased.vcf.gz
        """
}

// ─────────────────────────────────────────────────────────────
// 合併所有分片 + 建索引 + publish（bcftools 容器）
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
    # 分片未建索引；concat -a 需要索引，先各自建 tbi。
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
    // 依每個樣本的性別展開分片（涵蓋全基因體：diploid→phase、haploid→pass）
    ch_shards = ensemble_ch.flatMap { meta, vcf, tbi ->
        buildPhaseShards(meta.sex).collect { r -> tuple(meta, vcf, tbi, r[0], r[1], r[2]) }
    }
    WHATSHAP_SUBSET(ch_shards)

    // 併回該樣本的 BAM 後 phasing
    ch_phase_in = WHATSHAP_SUBSET.out
        .combine(bam_ch, by: 0)
        .map { meta, idx, type, sv, st, bam, bai, recal -> tuple(meta, idx, type, sv, st, bam, bai) }
    WHATSHAP_PHASE(ch_phase_in, fasta, fasta_fai)

    // 依 meta 收齊所有分片後合併（分片數因性別而異，故不指定 size）
    WHATSHAP_CONCAT(WHATSHAP_PHASE.out.groupTuple(by: 0))

    emit:
    vcf = WHATSHAP_CONCAT.out.vcf   // tuple(meta, phased_vcf, tbi)
}
