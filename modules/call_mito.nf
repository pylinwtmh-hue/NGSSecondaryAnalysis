/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
 * modules/call_mito.nf
 * =========================================================
 * 粒線體變異呼叫 sub-workflow（Lane 5）。把 modules/mitochondria.nf 的 12 個 MITO_* process
 * 組合起來：extract → (normal / shifted 兩路 bwa+sort) → Mutect2 → liftover → merge → filter。
 *
 * 對外只吃 bam_ch；chrM 參考檔一律從 params 於內部建立（自成一體）。emit 過濾後的 mito VCF
 * （MITO_FILTER 自己 publishDir 到 07_mitochondria，主流程通常不需再接）。
 *
 * 註：MITO_SORT_MARKDUP 需用兩次（normal / shifted），靠 include-time 別名區分 —— 這也是為何
 *     本 sub-workflow 獨立成一個 module 檔（Nextflow 無法對「同檔內定義」的 process 取別名），
 *     而不是寫進 mitochondria.nf 本身。
 */

include { MITO_EXTRACT_READS;
          MITO_BAM2FASTQ_NORMAL;
          MITO_BAM2FASTQ_SHIFTED;
          MITO_BWA_NORMAL;
          MITO_BWA_SHIFTED;
          MITO_SORT_MARKDUP as MITO_SORT_INDEX_NORMAL;
          MITO_SORT_MARKDUP as MITO_SORT_INDEX_SHIFTED;
          MITO_MUTECT2_NORMAL;
          MITO_MUTECT2_SHIFTED;
          MITO_LIFTOVER;
          MITO_MERGE;
          MITO_FILTER }                     from './mitochondria'

workflow CALL_MITO {
    take:
    bam_ch      // tuple(meta, bam, bai, ...)（與原 Lane 5 的 ch_bam 相同）

    main:
    ch_chrM_only_fasta    = file(params.chrM_only_fasta)
    ch_chrM_only_fai      = file("${params.chrM_only_fasta}.fai")
    ch_chrM_only_dict     = file(params.chrM_only_fasta.replace('.fasta', '.dict'))
    ch_chrM_shifted_fasta = file(params.chrM_shifted_fasta)
    ch_chrM_shifted_fai   = file("${params.chrM_shifted_fasta}.fai")
    ch_chrM_shifted_dict  = file(params.chrM_shifted_fasta.replace('.fasta', '.dict'))
    ch_shift_back         = file(params.chrM_shift_back)
    ch_chrM_blacklist     = file(params.chrM_blacklist)

    MITO_EXTRACT_READS(bam_ch)

    MITO_BAM2FASTQ_NORMAL(MITO_EXTRACT_READS.out.reads)
    MITO_BWA_NORMAL(
        MITO_BAM2FASTQ_NORMAL.out.reads,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        file("${params.chrM_only_fasta}.amb"),
        file("${params.chrM_only_fasta}.ann"),
        file("${params.chrM_only_fasta}.bwt"),
        file("${params.chrM_only_fasta}.pac"),
        file("${params.chrM_only_fasta}.sa")
    )
    MITO_SORT_INDEX_NORMAL(MITO_BWA_NORMAL.out.bam)

    MITO_BAM2FASTQ_SHIFTED(MITO_EXTRACT_READS.out.reads)
    MITO_BWA_SHIFTED(
        MITO_BAM2FASTQ_SHIFTED.out.reads,
        ch_chrM_shifted_fasta, ch_chrM_shifted_fai, ch_chrM_shifted_dict,
        file("${params.chrM_shifted_fasta}.amb"),
        file("${params.chrM_shifted_fasta}.ann"),
        file("${params.chrM_shifted_fasta}.bwt"),
        file("${params.chrM_shifted_fasta}.pac"),
        file("${params.chrM_shifted_fasta}.sa")
    )
    MITO_SORT_INDEX_SHIFTED(MITO_BWA_SHIFTED.out.bam)

    MITO_MUTECT2_NORMAL(
        MITO_SORT_INDEX_NORMAL.out.bam,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        ch_chrM_blacklist
    )
    MITO_MUTECT2_SHIFTED(
        MITO_SORT_INDEX_SHIFTED.out.bam,
        ch_chrM_shifted_fasta, ch_chrM_shifted_fai, ch_chrM_shifted_dict,
        ch_chrM_blacklist
    )

    MITO_LIFTOVER(
        MITO_MUTECT2_SHIFTED.out.vcf,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        ch_shift_back
    )
    // join() 確保同一樣本的 normal VCF 與 lifted VCF 配對（多樣本非同步不跨樣本錯配）
    ch_mito_merge_input = MITO_MUTECT2_NORMAL.out.vcf
        .join(MITO_LIFTOVER.out.vcf, by: 0)
        .map { meta, normal_vcf, normal_tbi, normal_stats,
                      lifted_vcf, lifted_tbi, lifted_stats ->
            [meta, normal_vcf, normal_tbi, normal_stats,
                   lifted_vcf, lifted_tbi, lifted_stats]
        }
    MITO_MERGE(ch_mito_merge_input, ch_chrM_only_dict)

    // MITO_FILTER：GATK 4.6 的 FilterMutectCalls 已無 --autosomal-coverage，
    // NuMT 過濾靠 --mitochondria-mode + blacklist mask（不需 mosdepth summary）。
    MITO_FILTER(
        MITO_MERGE.out.vcf,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        ch_chrM_blacklist,
        file("${params.chrM_blacklist}.idx")
    )

    emit:
    vcf = MITO_FILTER.out.vcf
}
