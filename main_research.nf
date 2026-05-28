#!/usr/bin/env nextflow

/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
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

nextflow.enable.dsl = 2

include { FASTP }                           from './modules/preprocessing'
include { PARABRICKS_FQ2BAM }               from './modules/alignment'
include { SAMTOOLS_STATS; MOSDEPTH }        from './modules/alignment_qc'
include { PARABRICKS_DEEPVARIANT;
          PARABRICKS_HAPLOTYPECALLER;
          GATK_HAPLOTYPECALLER;
          GATK_VQSR_SNP;
          GATK_VQSR_INDEL }                 from './modules/variant_calling'
include { CNVKIT_BATCH;
          MANTA_GERMLINE;
          GATK_COLLECT_READ_COUNTS;
          GATK_PLOIDY_CASE;
          GATK_GERMLINE_CNV_CASE;
          GATK_POSTPROCESS_CNV }            from './modules/cnv_sv'
include { EXPANSIONHUNTER }                 from './modules/repeat'
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
          MITO_FILTER }                     from './modules/mitochondria'
include { BGZIP_VCF as BGZIP_VCF_DV;
          BGZIP_VCF as BGZIP_VCF_HC;
          BCFTOOLS_ENSEMBLE;
          BCFTOOLS_STATS;
          BCFTOOLS_STATS as BCFTOOLS_STATS_ENSEMBLE;
          MULTIQC }                         from './modules/postprocessing'
include { AUTOMAP }                         from './modules/roh'

if (!params.input_csv) {
    error "錯誤：請提供 --input_csv 參數"
}

workflow {

    // =========================================================
    // (A) 讀取 Samplesheet
    // 格式：sample,fastq_1,fastq_2,sex
    // sex 填 male / female / unknown
    // =========================================================
    ch_input = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def lane  = row.lane ?: 'L001'  // 沒有 lane 欄位就預設 L001
            def meta  = [
                id:        row.lane ? "${row.sample}_${lane}" : row.sample,
                sample_id: row.sample,               // 真實 sample ID（FQ2BAM 用）
                sex:       row.sex ?: 'unknown',
                lane:      lane
            ]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // =========================================================
    // (B) Step 1: Preprocessing
    // =========================================================
    FASTP(ch_input)

    // =========================================================
    // (C) Step 2: Alignment + BQSR
    // =========================================================
    ch_fasta      = file(params.fasta)
    ch_fasta_fai  = file("${params.fasta}.fai")
    ch_fasta_dict = file(params.fasta.replace('.fasta', '.dict'))

    ch_fq2bam_input = FASTP.out.reads
        .map { meta, reads -> 
            // 把 meta.id 換回 sample_id，保留 sex
            def new_meta = [id: meta.sample_id, sex: meta.sex]
            [new_meta, reads]
        }
        .groupTuple()   // 同一 sample_id 的所有 lane reads 合併成 list
        .map { meta, reads_list ->
            [meta, reads_list.flatten()]  // [[R1a,R2a],[R1b,R2b]] → [R1a,R2a,R1b,R2b]
        }

    PARABRICKS_FQ2BAM(
        ch_fq2bam_input,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        file("${params.fasta}.amb"),
        file("${params.fasta}.ann"),
        file("${params.fasta}.bwt"),
        file("${params.fasta}.pac"),
        file("${params.fasta}.sa"),
        file(params.dbsnp),        file("${params.dbsnp}.tbi"),
        file(params.known_indels), file("${params.known_indels}.tbi"),
        file(params.known_indels2),file("${params.known_indels2}.tbi"),
        file(params.known_snps),   file("${params.known_snps}.tbi")
    )
    ch_bam = PARABRICKS_FQ2BAM.out.alignment_bundle

    // =========================================================
    // (D) Step 3: Alignment QC
    // =========================================================
    ch_mosdepth_targets = (params.seq_type == "WES") ?
        file(params.wes_targets) : file("NO_FILE")

    SAMTOOLS_STATS(ch_bam)
    MOSDEPTH(ch_bam, ch_mosdepth_targets)

    // =========================================================
    // (E) Step 4: Parallel Variant Calling
    // =========================================================

    // Lane 1: DeepVariant (GPU)
    PARABRICKS_DEEPVARIANT(ch_bam, ch_fasta, ch_fasta_fai, ch_fasta_dict)
    BGZIP_VCF_DV(PARABRICKS_DEEPVARIANT.out.vcf)
    ch_dv_vcf = BGZIP_VCF_DV.out.vcf

    // Lane 2a: HaplotypeCaller (GPU)
    // 依賴 DeepVariant 完成後才啟動，確保同一張 GPU 不被同時佔用
    // join 讓同一個 sample 的 HaplotypeCaller 等 DeepVariant 的 channel 發射後才觸發
    ch_bam_after_dv = ch_bam
        .join(BGZIP_VCF_DV.out.vcf.map { meta, vcf, tbi -> [meta, meta.id] })
        .map { meta, bam, bai, recal, dummy -> [meta, bam, bai, recal] }

    PARABRICKS_HAPLOTYPECALLER(ch_bam_after_dv, ch_fasta, ch_fasta_fai, ch_fasta_dict)
    BGZIP_VCF_HC(PARABRICKS_HAPLOTYPECALLER.out.vcf)
    ch_hc_vcf_raw = BGZIP_VCF_HC.out.vcf

    // Lane 2b: VQSR (WGS only)
    if (params.seq_type == "WGS") {
        GATK_VQSR_SNP(
            ch_hc_vcf_raw, 
            ch_fasta, ch_fasta_fai, ch_fasta_dict,
            file(params.hapmap),     file("${params.hapmap}.tbi"),
            file(params.omni),       file("${params.omni}.tbi"),
            file(params.known_snps), file("${params.known_snps}.tbi"),
            file(params.dbsnp),      file("${params.dbsnp}.tbi")
        )
        GATK_VQSR_INDEL(
            GATK_VQSR_SNP.out.vcf,
            ch_fasta, ch_fasta_fai, ch_fasta_dict,
            file(params.known_indels), file("${params.known_indels}.tbi"),
            file(params.axiom),        file("${params.axiom}.tbi"),
            file(params.dbsnp),        file("${params.dbsnp}.tbi")
        )
        ch_filtered_hc_vcf = GATK_VQSR_INDEL.out.vcf
    } else {
        ch_filtered_hc_vcf = ch_hc_vcf_raw  // ← WES 直接用壓縮後的
    }
    
    // Lane 3a: CNVkit (WGS only)
    if (params.seq_type == "WGS") {
        ch_cnvkit_pon = params.cnvkit_pon ? file(params.cnvkit_pon) : file("NO_FILE")
        CNVKIT_BATCH(
            ch_bam,
            ch_dv_vcf,
            ch_fasta,
            ch_cnvkit_pon
        )
    }

    // Lane 3b: Manta SV calling
    MANTA_GERMLINE(ch_bam, ch_fasta, ch_fasta_fai)

    // Lane 3c: gCNV (WES only，需 --run_gcnv true 且已有 PON)
    if (params.seq_type == "WES" && params.run_gcnv) {
        ch_gcnv_intervals    = file(params.gcnv_pon_dir)
        ch_ploidy_model      = file(params.gcnv_ploidy_model_dir)
        ch_model_shards_list = Channel.fromPath("${params.gcnv_model_dir}/**/*-model", type: 'dir').collect()
        ch_model_shards_flat = Channel.fromPath("${params.gcnv_model_dir}/**/*-model", type: 'dir')
        
        GATK_COLLECT_READ_COUNTS(
            ch_bam, ch_fasta, ch_fasta_fai, ch_fasta_dict,
            ch_gcnv_intervals
        )
        GATK_PLOIDY_CASE(
            GATK_COLLECT_READ_COUNTS.out.counts,
            ch_ploidy_model
        )
        ch_gcnv_caller_in = GATK_COLLECT_READ_COUNTS.out.counts
            .join(GATK_PLOIDY_CASE.out.ploidy_calls)
            .combine(ch_model_shards_flat)

        GATK_GERMLINE_CNV_CASE(ch_gcnv_caller_in)

        ch_postprocess_in = GATK_GERMLINE_CNV_CASE.out.call_shard
            .groupTuple()
            .join(GATK_PLOIDY_CASE.out.ploidy_calls)

        GATK_POSTPROCESS_CNV(
            ch_postprocess_in,
            ch_model_shards_list,
            ch_fasta_dict
        )
    } else if (params.seq_type == "WGS" && params.run_gcnv) {
    log.warn "WGS 模式不支援 gCNV，忽略 --run_gcnv 參數"
    }

    // Lane 4: STR (ExpansionHunter)
    EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fasta_fai, file(params.str_catalog))

    // =========================================================
    // Lane 5: Mitochondria variant calling
    // =========================================================
    ch_chrM_only_fasta    = file(params.chrM_only_fasta)
    ch_chrM_only_fai      = file("${params.chrM_only_fasta}.fai")
    ch_chrM_only_dict     = file(params.chrM_only_fasta.replace('.fasta', '.dict'))
    ch_chrM_shifted_fasta = file(params.chrM_shifted_fasta)
    ch_chrM_shifted_fai   = file("${params.chrM_shifted_fasta}.fai")
    ch_chrM_shifted_dict  = file(params.chrM_shifted_fasta.replace('.fasta', '.dict'))
    ch_shift_back         = file(params.chrM_shift_back)
    ch_chrM_blacklist     = file(params.chrM_blacklist)

    MITO_EXTRACT_READS(ch_bam)

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
    MITO_MERGE(
        MITO_MUTECT2_NORMAL.out.vcf,
        MITO_LIFTOVER.out.vcf,
        ch_chrM_only_dict
    )
    MITO_FILTER(
        MITO_MERGE.out.vcf,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        ch_chrM_blacklist,
        file("${params.chrM_blacklist}.idx")
    )

    // =========================================================
    // (F) Step 5: Post-processing
    // =========================================================
    // join() 確保同一樣本的 DV 和 HC VCF 配對
    ch_ensemble_input = ch_dv_vcf
        .join(ch_filtered_hc_vcf, by: 0)
        .map { meta, dv_vcf, dv_tbi, hc_vcf, hc_tbi ->
            [meta, dv_vcf, dv_tbi, hc_vcf, hc_tbi]
        }

    BCFTOOLS_ENSEMBLE(ch_ensemble_input)

    // AutoMap ROH：用 HaplotypeCaller VCF（保留 AD 欄位，VQSR 後可能移除）
    AUTOMAP(ch_hc_vcf_raw)

    BCFTOOLS_STATS(ch_dv_vcf)
    BCFTOOLS_STATS_ENSEMBLE(BCFTOOLS_ENSEMBLE.out.vcf)

    ch_multiqc = Channel.empty()
        .mix(FASTP.out.json)
        .mix(PARABRICKS_FQ2BAM.out.qc_metrics)
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(MOSDEPTH.out.global_dist)
        .mix(MOSDEPTH.out.summary)
        .mix(BCFTOOLS_STATS.out.stats)
        .mix(BCFTOOLS_STATS_ENSEMBLE.out.stats)
        .collect()

    MULTIQC(ch_multiqc)
}