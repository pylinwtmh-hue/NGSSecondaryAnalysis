/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
 * modules/call_snv.nf
 * =========================================================
 * SNV/indel calling ＋ ensemble sub-workflow（原 Lane 1 / 2 ＋ Step 5 的 phasing/ensemble）。
 *
 * 多 process 的 caller 各自為第二層 sub-workflow：
 *   DEEPVARIANT     = PARABRICKS_DEEPVARIANT → BGZIP_VCF_DV
 *   HAPLOTYPECALLER = PARABRICKS_HAPLOTYPECALLER → BGZIP_VCF_HC
 *   VQSR（WGS only） = GATK_VQSR_SNP → GATK_VQSR_INDEL
 * 之後 PHASE_COMBINE（現有 sub-workflow，--run_phasing 時，各 caller 先 phase+combine）
 *   → BCFTOOLS_ENSEMBLE（reheader/norm/merge/+fixploidy）→ ensemble.fixed。
 *
 * 需跨 module 且對 BGZIP_VCF / BCFTOOLS_STATS 取兩次別名，故獨立成一個 module 檔
 * （Nextflow 無法對同檔內定義的 process 取別名）。
 *
 * emit：
 *   dv_vcf         → CALL_CNV_SV（CNVkit 的 SNP b-allele 輸入）
 *   hc_vcf_raw     → CALL_ROH（保留 GT/AD 的 HC 原始 VCF）
 *   ensemble       → 發布 04_snv_indel + 三級讀取
 *   dv_stats / ensemble_stats → MultiQC
 */

include { PARABRICKS_DEEPVARIANT;
          PARABRICKS_HAPLOTYPECALLER;
          GATK_VQSR_SNP;
          GATK_VQSR_INDEL }                 from './variant_calling'
include { BGZIP_VCF as BGZIP_VCF_DV;
          BGZIP_VCF as BGZIP_VCF_HC;
          BCFTOOLS_ENSEMBLE;
          BCFTOOLS_STATS;
          BCFTOOLS_STATS as BCFTOOLS_STATS_ENSEMBLE } from './postprocessing'
include { PHASE_COMBINE }                   from './phasing'


// ── DeepVariant（GPU）→ bgzip ─────────────────────────────────
workflow DEEPVARIANT {
    take:
    bam_ch
    fasta
    fasta_fai
    fasta_dict

    main:
    PARABRICKS_DEEPVARIANT(bam_ch, fasta, fasta_fai, fasta_dict)
    BGZIP_VCF_DV(PARABRICKS_DEEPVARIANT.out.vcf)

    emit:
    vcf = BGZIP_VCF_DV.out.vcf
}

// ── HaplotypeCaller（GPU）→ bgzip。bam_ch 應已 join 過 DV 輸出（同張 GPU 序列化）──
workflow HAPLOTYPECALLER {
    take:
    bam_ch
    fasta
    fasta_fai
    fasta_dict

    main:
    PARABRICKS_HAPLOTYPECALLER(bam_ch, fasta, fasta_fai, fasta_dict)
    BGZIP_VCF_HC(PARABRICKS_HAPLOTYPECALLER.out.vcf)

    emit:
    vcf = BGZIP_VCF_HC.out.vcf
}

// ── VQSR（WGS）：SNP → INDEL；resource 檔內部依 params 建立 ──────
workflow VQSR {
    take:
    hc_vcf
    fasta
    fasta_fai
    fasta_dict

    main:
    GATK_VQSR_SNP(
        hc_vcf,
        fasta, fasta_fai, fasta_dict,
        file(params.hapmap),     file("${params.hapmap}.tbi"),
        file(params.omni),       file("${params.omni}.tbi"),
        file(params.known_snps), file("${params.known_snps}.tbi"),
        file(params.dbsnp),      file("${params.dbsnp}.tbi")
    )
    GATK_VQSR_INDEL(
        GATK_VQSR_SNP.out.vcf,
        fasta, fasta_fai, fasta_dict,
        file(params.known_indels), file("${params.known_indels}.tbi"),
        file(params.axiom),        file("${params.axiom}.tbi"),
        file(params.dbsnp),        file("${params.dbsnp}.tbi")
    )

    emit:
    vcf = GATK_VQSR_INDEL.out.vcf
}


workflow CALL_SNV {
    take:
    bam_ch      // tuple(meta, bam, bai, recal)

    main:
    ch_fasta      = file(params.fasta)
    ch_fasta_fai  = file("${params.fasta}.fai")
    ch_fasta_dict = file(params.fasta.replace('.fasta', '.dict'))

    // Lane 1: DeepVariant
    DEEPVARIANT(bam_ch, ch_fasta, ch_fasta_fai, ch_fasta_dict)
    ch_dv_vcf = DEEPVARIANT.out.vcf

    // Lane 2a: HaplotypeCaller（等 DV 發射後才觸發，避免同張 GPU 同時佔用）
    ch_bam_after_dv = bam_ch
        .join(ch_dv_vcf.map { meta, vcf, tbi -> [meta, meta.id] })
        .map { meta, bam, bai, recal, dummy -> [meta, bam, bai, recal] }
    HAPLOTYPECALLER(ch_bam_after_dv, ch_fasta, ch_fasta_fai, ch_fasta_dict)
    ch_hc_vcf_raw = HAPLOTYPECALLER.out.vcf

    // Lane 2b: VQSR（WGS）；WES 直接用 raw
    if (params.seq_type == "WGS") {
        VQSR(ch_hc_vcf_raw, ch_fasta, ch_fasta_fai, ch_fasta_dict)
        ch_filtered_hc_vcf = VQSR.out.vcf
    } else {
        ch_filtered_hc_vcf = ch_hc_vcf_raw
    }

    // Step 5:（選用 --run_phasing）各 caller 先 phase+combine 再進 ensemble；否則直接 join。
    //   在 merge「之前」、各 caller 還單樣本 biallelic 時做，compound（SUZ12）才 phase/合得成。
    if (params.run_phasing) {
        ch_combine_py = file("${projectDir}/scripts/combine_phased.py")
        ch_callers = ch_dv_vcf.map { meta, vcf, tbi -> tuple(meta, 'DV', vcf, tbi) }
            .mix( ch_filtered_hc_vcf.map { meta, vcf, tbi -> tuple(meta, 'HC', vcf, tbi) } )
        PHASE_COMBINE(ch_callers, bam_ch, ch_fasta, ch_fasta_fai, ch_combine_py)

        ch_dv_ready = PHASE_COMBINE.out.vcf
            .filter { it[1] == 'DV' }.map { meta, c, vcf, tbi -> tuple(meta, vcf, tbi) }
        ch_hc_ready = PHASE_COMBINE.out.vcf
            .filter { it[1] == 'HC' }.map { meta, c, vcf, tbi -> tuple(meta, vcf, tbi) }
        ch_ensemble_input = ch_dv_ready.join(ch_hc_ready, by: 0)
    } else {
        // join() 確保同一樣本的 DV/HC 配對（多樣本非同步完成不會跨樣本錯配）
        ch_ensemble_input = ch_dv_vcf.join(ch_filtered_hc_vcf, by: 0)
    }

    // 性別感知倍體定義（單一真相來源；+fixploidy 用）
    ch_sex_ploidy = file(params.sex_ploidy_file)
    BCFTOOLS_ENSEMBLE(ch_ensemble_input, ch_sex_ploidy)

    // 統計（→ MultiQC）
    BCFTOOLS_STATS(ch_dv_vcf)
    BCFTOOLS_STATS_ENSEMBLE(BCFTOOLS_ENSEMBLE.out.vcf)

    emit:
    dv_vcf         = ch_dv_vcf
    hc_vcf_raw     = ch_hc_vcf_raw
    ensemble       = BCFTOOLS_ENSEMBLE.out.vcf
    dv_stats       = BCFTOOLS_STATS.out.stats
    ensemble_stats = BCFTOOLS_STATS_ENSEMBLE.out.stats
}
