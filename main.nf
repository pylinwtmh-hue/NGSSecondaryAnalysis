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
 *  * Optional callers are controlled by flags (all default OFF):
 *  *   - Manta (Illumina, PolyForm Strict 1.0.0)            --run_manta
 *  *   - ExpansionHunter (Illumina, PolyForm Strict 1.0.0)  --run_expansionhunter
 *  *   - ROH via bcftools roh (MIT/GPL, commercial OK)      --run_roh
 *  *   - ROH via AutoMap (no license published)             --run_automap
 *  * Defaults reproduce the evaluated outputs (SNV/indel, CNV/SV, STR, mito).
 *  * ROH is off by default and is not part of the evaluation.
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
          DELLY_GERMLINE;
          BCFTOOLS_CONVERT_DELLY;
          MANTA_GERMLINE;
          GATK_COLLECT_READ_COUNTS;
          GATK_PLOIDY_CASE;
          GATK_GERMLINE_CNV_CASE;
          GATK_POSTPROCESS_CNV }            from './modules/cnv_sv'
include { GANGSTR_CHROM;
          GANGSTR_MERGE;
          EXPANSIONHUNTER }                from './modules/repeat'
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
include { AUTOMAP; BCFTOOLS_ROH }           from './modules/roh'

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
            // 把 meta.id 換回 sample_id，保留 sex 和 lane
            // lane 資訊保留在 lanes list，供 alignment.nf 組 RG 字串用
            def new_meta = [id: meta.sample_id, sex: meta.sex]
            [new_meta, meta.lane, reads]
        }
        .groupTuple()   // 同一 sample_id 的所有 lane reads 合併成 list
        .map { meta, lanes_list, reads_list ->
            // 按 lane 排序，確保 lanes 和 reads 順序一致
            def sorted = [lanes_list, reads_list].transpose()
                .sort { a, b -> a[0] <=> b[0] }
            def sorted_lanes = sorted.collect { it[0] }
            def sorted_reads = sorted.collect { it[1] }.flatten()
            def new_meta = meta + [lanes: sorted_lanes] // 把 lanes 存進 meta，reads 展平
            [new_meta, sorted_reads]   // [[R1a,R2a],[R1b,R2b]] → [R1a,R2a,R1b,R2b]
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

    // WGS 深度 QC 只看 autosome primary contig（chr1-22）
    // 排除 chrM（高拷貝數）、chrX/chrY（受性別影響）、unplaced contig
    ch_autosome_bed = (params.seq_type == "WGS") ?
        file(params.autosome_bed) : file("NO_FILE")

    SAMTOOLS_STATS(ch_bam)
    MOSDEPTH(ch_bam, ch_mosdepth_targets, ch_autosome_bed)
    // (meta, summary.txt)：mito NuMT filter 與 MultiQC 共用同一輸出
    ch_mosdepth_summary = MOSDEPTH.out.summary

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

    // Lane 3b: Delly SV calling（替代 Manta，BSD license）
    ch_delly_excl = params.delly_excl ? file(params.delly_excl) : file("NO_FILE")
    DELLY_GERMLINE(ch_bam, ch_fasta, ch_fasta_fai, ch_delly_excl)
    BCFTOOLS_CONVERT_DELLY(DELLY_GERMLINE.out.bcf)

    // Lane 3b（選用）: Manta SV calling（--run_manta，預設關閉；非商用授權）
    if (params.run_manta) {
        MANTA_GERMLINE(ch_bam, ch_fasta, ch_fasta_fai)
    }

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

    // Lane 4: STR (GangSTR，替代 ExpansionHunter，GPL v3 license)
    // WGS：按染色體平行化（24 個 process），大幅縮短執行時間
    // WES：也平行化，但 loci 較少，效果有限
    def gangstr_regions = params.seq_type == "WES"
        ? file(params.gangstr_regions_wes)
        : file(params.gangstr_regions_wgs)

    // 展開 24 個染色體，每個樣本 × 每條染色體 = 一個 GANGSTR_CHROM process
    def chroms = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
    ch_bam_chrom = ch_bam.combine(Channel.from(chroms))

    GANGSTR_CHROM(ch_bam_chrom, ch_fasta, ch_fasta_fai, gangstr_regions)

    // 按樣本收集 24 個 VCF，按染色體順序排序後傳入 GANGSTR_MERGE
    ch_gangstr_vcfs = GANGSTR_CHROM.out.vcf
        .map { meta, chrom, vcf -> [meta, chrom, vcf] }
        .groupTuple(by: 0)
        .map { meta, chroms_list, vcfs ->
            // 按染色體順序排序
            def order = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
            def sorted_vcfs = [chroms_list, vcfs].transpose()
                .sort { a, b -> order.indexOf(a[0]) <=> order.indexOf(b[0]) }
                .collect { it[1] }
            [meta, sorted_vcfs]
        }

    GANGSTR_MERGE(ch_gangstr_vcfs)

    // Lane 4（選用）: ExpansionHunter（--run_expansionhunter，預設關閉；非商用授權）
    if (params.run_expansionhunter) {
        EXPANSIONHUNTER(ch_bam, ch_fasta, ch_fasta_fai, file(params.str_catalog))
    }

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
    // join() 確保同一個樣本的 normal VCF 和 lifted VCF 配對
    // 沒有 join 的話，Nextflow 按 queue 順序配對，多樣本非同步完成時會跨樣本錯配
    ch_mito_merge_input = MITO_MUTECT2_NORMAL.out.vcf
        .join(MITO_LIFTOVER.out.vcf, by: 0)
        .map { meta, normal_vcf, normal_tbi, normal_stats,
                      lifted_vcf, lifted_tbi, lifted_stats ->
            [meta, normal_vcf, normal_tbi, normal_stats,
                   lifted_vcf, lifted_tbi, lifted_stats]
        }
    MITO_MERGE(
        ch_mito_merge_input,
        ch_chrM_only_dict
        )
    // MITO_FILTER 需要 autosomal coverage（NuMT filter），join mosdepth summary（依樣本配對）
    ch_mito_filter_input = MITO_MERGE.out.vcf.join(ch_mosdepth_summary, by: 0)
    MITO_FILTER(
        ch_mito_filter_input,
        ch_chrM_only_fasta, ch_chrM_only_fai, ch_chrM_only_dict,
        ch_chrM_blacklist,
        file("${params.chrM_blacklist}.idx")
    )

    // =========================================================
    // (F) Step 5: Post-processing
    // =========================================================
    // join() 確保同一樣本的 DV 和 HC VCF 配對
    // 沒有 join 的話，DV 和 HC 執行時間不同，多樣本時會跨樣本錯配
    ch_ensemble_input = ch_dv_vcf
        .join(ch_filtered_hc_vcf, by: 0)
        .map { meta, dv_vcf, dv_tbi, hc_vcf, hc_tbi ->
            [meta, dv_vcf, dv_tbi, hc_vcf, hc_tbi]
        }

    BCFTOOLS_ENSEMBLE(ch_ensemble_input)

    // ROH（選用，皆預設關閉；ROH 不納入評鑑）：用 HaplotypeCaller VCF（保留 GT/AD）。
    //   --run_roh     → bcftools roh（MIT/GPL，可商用）
    //   --run_automap → AutoMap（無公開授權，僅非商用/研究）
    if (params.run_roh) {
        BCFTOOLS_ROH(ch_hc_vcf_raw)
    }
    if (params.run_automap) {
        AUTOMAP(ch_hc_vcf_raw)
    }

    BCFTOOLS_STATS(ch_dv_vcf)
    BCFTOOLS_STATS_ENSEMBLE(BCFTOOLS_ENSEMBLE.out.vcf)

    ch_multiqc = Channel.empty()
        .mix(FASTP.out.json)
        .mix(PARABRICKS_FQ2BAM.out.qc_metrics)
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(MOSDEPTH.out.global_dist)
        .mix(ch_mosdepth_summary.map { meta, f -> f })
        .mix(BCFTOOLS_STATS.out.stats)
        .mix(BCFTOOLS_STATS_ENSEMBLE.out.stats)
        .collect()

    MULTIQC(ch_multiqc)
}