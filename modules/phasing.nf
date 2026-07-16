/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Phasing + Combine Module
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
 * DISCLAIMER: Provided "as is" without warranty. Users are solely responsible
 * for validating and interpreting all results. See LICENSE.
 * =========================================================
 * modules/phasing.nf
 * ==================
 * NCKUH 專用：對「各 caller 原始單樣本 VCF」做 read-backed phasing（WhatsHap）
 * 後，用 combine_phased.py 把相鄰/重疊的 cis 變異合成 canonical MNV —— 全部在
 * BCFTOOLS_ENSEMBLE 合併「之前」完成（見討論）。
 *
 * 為什麼在 merge 之前、且各 caller 各自做？
 *   - ensemble 是 DV/HC 聯集，兩 caller 對同一 compound 的表示法常不同，merge 會
 *     產生 multiallelic，whatshap 跳過 multiallelic → compound 拿不到 phase。
 *   - 在「還是單一 caller、還是 biallelic」時 phase+combine，SUZ12 這類 del+ins 才
 *     phase 得到、也才合得成單筆（HC→GAAA>GTT、DV→GAAA>GAT）。merge 後各自 PS 仍保留。
 *
 * 為什麼不再需要 sex-aware 倍體切分？
 *   - phasing 在 fixploidy「之前」，原始 caller VCF 全基因體皆 diploid（無混合倍體）
 *     → 不會有 PloidyError → 不需依性別切 chrX PAR。倍體由後面 BCFTOOLS_ENSEMBLE 的
 *     +fixploidy 校正。
 *
 * 分片（純為平行度 + 完整性，無性別邏輯）：
 *   - 主要 contig chr1-22/X/Y/M → 'phase'（送 whatshap）
 *   - 其餘所有 contig（alt/decoy/random/Un/HLA）→ 'pass'（原樣保留、不 phase，非破壞性）
 *
 * 各分片用單一工具容器（符合本 pipeline 慣例）：
 *   WHATSHAP_SUBSET → bcftools（切出該分片；補集用 -t ^）
 *   WHATSHAP_PHASE  → whatshap（type=phase 才 phase；單樣本用 --ignore-read-groups）
 *   WHATSHAP_CONCAT → bcftools（合併分片 + 索引）
 *   COMBINE_PHASED  → tertiary_python（python3 跑 combine_phased.py + bcftools 排序/索引）
 *
 * 由 params.run_phasing 開關（預設 false）；DRAGEN 路徑自帶 PS，其 combine 在三級做。
 */

// ─────────────────────────────────────────────────────────────
// phasing 分片：主要 contig → phase；其餘 contig（補集）→ pass
//   回傳每筆 [region, type, idx]
// ─────────────────────────────────────────────────────────────
def phaseRegions() {
    def sh = []
    (1..22).each { sh << ["chr${it}".toString(), 'phase'] }
    sh << ['chrX', 'phase']
    sh << ['chrY', 'phase']
    sh << ['chrM', 'phase']
    // 其餘所有 contig 一律 passthrough（非破壞性；-r 不支援 ^ 補集，SUBSET 用 -t）
    def primary = ((1..22).collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']).join(',')
    sh << ["^${primary}".toString(), 'pass']
    return sh.withIndex().collect { e, i -> [e[0], e[1], i] }
}

// ─────────────────────────────────────────────────────────────
// 切出該分片區段（bcftools）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_SUBSET {
    tag "${meta.id}:${caller}:${region}"
    label 'process_low'

    input:
    tuple val(meta), val(caller), path(vcf), path(tbi), val(region), val(type), val(idx)

    output:
    tuple val(meta), val(caller), val(idx), val(type),
          path("${meta.id}.${caller}.s${idx}.sub.vcf.gz"),
          path("${meta.id}.${caller}.s${idx}.sub.vcf.gz.tbi")

    script:
    """
    # 一般 contig 用 -r（走索引）；補集（^ 開頭，收其餘 contig）用 -t（-r 不支援 ^）。
    case "${region}" in
      ^*) bcftools view -t "${region}" ${vcf} -Oz -o ${meta.id}.${caller}.s${idx}.sub.vcf.gz ;;
      *)  bcftools view -r "${region}" ${vcf} -Oz -o ${meta.id}.${caller}.s${idx}.sub.vcf.gz ;;
    esac
    tabix -p vcf ${meta.id}.${caller}.s${idx}.sub.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 分片 phasing（whatshap；type=phase 才 phase）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_PHASE {
    tag "${meta.id}:${caller}:s${idx}(${type})"
    label 'process_medium'

    input:
    tuple val(meta), val(caller), val(idx), val(type), path(sub_vcf), path(sub_tbi), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val(caller), path("${meta.id}.${caller}.s${idx}.phased.vcf.gz")

    script:
    if (type == 'phase')
        """
        # 單一樣本 VCF：--ignore-read-groups 就把全部 reads 指到該唯一樣本，不必 --sample。
        # --reference 開 re-alignment 對 indel phasing 較準。空 contig 不會報錯（原樣輸出）。
        whatshap phase \\
            --reference ${fasta} \\
            --ignore-read-groups \\
            -o ${meta.id}.${caller}.s${idx}.phased.vcf.gz \\
            ${sub_vcf} ${bam}
        """
    else
        """
        # 非主要 contig：passthrough，不 phase、不加 PS。
        cp ${sub_vcf} ${meta.id}.${caller}.s${idx}.phased.vcf.gz
        """
}

// ─────────────────────────────────────────────────────────────
// 合併分片（bcftools）
// ─────────────────────────────────────────────────────────────
process WHATSHAP_CONCAT {
    tag "${meta.id}:${caller}"
    label 'process_low'

    input:
    tuple val(meta), val(caller), path(phased_vcfs)

    output:
    tuple val(meta), val(caller),
          path("${meta.id}.${caller}.phased.vcf.gz"),
          path("${meta.id}.${caller}.phased.vcf.gz.tbi")

    script:
    """
    for f in ${phased_vcfs}; do tabix -f -p vcf \$f; done
    bcftools concat -a ${phased_vcfs} -Oz -o concat.vcf.gz
    bcftools sort concat.vcf.gz -Oz -o ${meta.id}.${caller}.phased.vcf.gz
    tabix -p vcf ${meta.id}.${caller}.phased.vcf.gz
    rm -f concat.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────
// 合成相鄰/重疊 cis 變異為 MNV（combine_phased.py；tertiary_python 容器）
// ─────────────────────────────────────────────────────────────
process COMBINE_PHASED {
    tag "${meta.id}:${caller}"
    label 'process_low'

    input:
    tuple val(meta), val(caller), path(phased_vcf), path(phased_tbi)
    path fasta
    path fasta_fai
    path combine_py

    output:
    tuple val(meta), val(caller),
          path("${meta.id}.${caller}.phased.combined.vcf.gz"),
          path("${meta.id}.${caller}.phased.combined.vcf.gz.tbi")

    script:
    """
    # combine_phased.py 只用 Python 標準庫，讀 bgzip VCF、自帶 faidx（讀 \${fasta}.fai）。
    python3 ${combine_py} \\
        --in ${phased_vcf} \\
        --out ${meta.id}.${caller}.combined.vcf \\
        --fasta ${fasta} \\
        --max-gap ${params.combine_max_gap}
    # 用 bcftools 排序 + 索引（-Oz 自帶 bgzip、index -t 免 tabix），故 tertiary_python 即可。
    bcftools sort ${meta.id}.${caller}.combined.vcf \\
        -Oz -o ${meta.id}.${caller}.phased.combined.vcf.gz
    bcftools index -t ${meta.id}.${caller}.phased.combined.vcf.gz
    rm -f ${meta.id}.${caller}.combined.vcf
    """
}

// ─────────────────────────────────────────────────────────────
// 組合 workflow：吃「各 caller 單樣本 VCF」+ BAM，輸出 phased+combined 單樣本 VCF
// ─────────────────────────────────────────────────────────────
workflow PHASE_COMBINE {
    take:
    caller_ch     // tuple(meta, caller, vcf, tbi)  caller ∈ {'DV','HC'}
    bam_ch        // tuple(meta, bam, bai, recal)  ← alignment_bundle
    fasta
    fasta_fai
    combine_py    // file: scripts/combine_phased.py（staged，免綁定顧慮）

    main:
    // 每個 (sample, caller) 展開成分片
    ch_shards = caller_ch.flatMap { meta, caller, vcf, tbi ->
        phaseRegions().collect { r -> tuple(meta, caller, vcf, tbi, r[0], r[1], r[2]) }
    }
    WHATSHAP_SUBSET(ch_shards)

    // 併回該樣本 BAM（by meta）後 phasing
    ch_phase_in = WHATSHAP_SUBSET.out
        .combine(bam_ch, by: 0)
        .map { meta, caller, idx, type, sv, st, bam, bai, recal ->
               tuple(meta, caller, idx, type, sv, st, bam, bai) }
    WHATSHAP_PHASE(ch_phase_in, fasta, fasta_fai)

    // 依 (meta, caller) 收齊分片後合併
    WHATSHAP_CONCAT(WHATSHAP_PHASE.out.groupTuple(by: [0, 1]))

    // 合成 MNV
    COMBINE_PHASED(WHATSHAP_CONCAT.out, fasta, fasta_fai, combine_py)

    emit:
    vcf = COMBINE_PHASED.out   // tuple(meta, caller, combined_vcf, tbi)
}
