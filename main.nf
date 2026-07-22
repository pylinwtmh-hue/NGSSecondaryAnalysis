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
include { ALIGNMENT_QC }                    from './modules/alignment_qc'
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
include { CALL_STR }                        from './modules/repeat'
include { CALL_MITO }                       from './modules/call_mito'
include { BGZIP_VCF as BGZIP_VCF_DV;
          BGZIP_VCF as BGZIP_VCF_HC;
          BCFTOOLS_ENSEMBLE;
          BCFTOOLS_STATS;
          BCFTOOLS_STATS as BCFTOOLS_STATS_ENSEMBLE;
          MULTIQC }                         from './modules/postprocessing'
include { AUTOMAP; BCFTOOLS_ROH }           from './modules/roh'
include { PHASE_COMBINE }                   from './modules/phasing'

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
    // (D) Step 3: Alignment QC（sub-workflow：SAMTOOLS_STATS + MOSDEPTH + PLOIDY_CHECK）
    // =========================================================
    ALIGNMENT_QC(ch_bam)
    // (meta, summary.txt)：MultiQC 用；PLOIDY_CHECK 在 sub-workflow 內部已消費同一輸出
    ch_mosdepth_summary = ALIGNMENT_QC.out.summary

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

    // Lane 4: STR（sub-workflow：GangSTR 依染色體平行化 → 合併；選用 ExpansionHunter）
    CALL_STR(ch_bam)

    // =========================================================
    // Lane 5: Mitochondria variant calling（sub-workflow）
    // =========================================================
    CALL_MITO(ch_bam)

    // =========================================================
    // (F) Step 5: Post-processing
    // =========================================================
    // (選用，--run_phasing，僅 NCKUH) 各 caller 先 phase + combine，再進 ensemble：
    //   在 merge「之前」、各 caller 還是單樣本 biallelic 時做，SUZ12 這類 del+ins 才
    //   phase 得到、也才合得成單筆 MNV（merge 後 multiallelic 會讓 whatshap 跳過）。
    //   非破壞性：產出的 ensemble.fixed 直接帶 phase(PS/|) + 已合成的 compound MNV，
    //   三級 prepare_vcf 照舊讀 ensemble.fixed 即可。DRAGEN 自帶 PS，其 combine 在三級做。
    //   預設關閉；在 DGX 驗證後以 --run_phasing true 開啟。
    if (params.run_phasing) {
        ch_combine_py = file("${projectDir}/scripts/combine_phased.py")
        ch_callers = ch_dv_vcf.map { meta, vcf, tbi -> tuple(meta, 'DV', vcf, tbi) }
            .mix( ch_filtered_hc_vcf.map { meta, vcf, tbi -> tuple(meta, 'HC', vcf, tbi) } )
        PHASE_COMBINE(ch_callers, ch_bam, ch_fasta, ch_fasta_fai, ch_combine_py)

        ch_dv_ready = PHASE_COMBINE.out.vcf
            .filter { it[1] == 'DV' }.map { meta, c, vcf, tbi -> tuple(meta, vcf, tbi) }
        ch_hc_ready = PHASE_COMBINE.out.vcf
            .filter { it[1] == 'HC' }.map { meta, c, vcf, tbi -> tuple(meta, vcf, tbi) }
        ch_ensemble_input = ch_dv_ready.join(ch_hc_ready, by: 0)
    } else {
        // join() 確保同一樣本的 DV/HC 配對（多樣本非同步完成不會跨樣本錯配）
        ch_ensemble_input = ch_dv_vcf.join(ch_filtered_hc_vcf, by: 0)
    }

    // 性別感知倍體定義（單一真相來源；+fixploidy 用，未來 ploidy-aware calling 也由此推導）
    //   不用 checkIfExists：讓 -preview 不需先放好參考檔；缺檔時仍會在實跑 staging 階段 fail-loud。
    ch_sex_ploidy = file(params.sex_ploidy_file)
    BCFTOOLS_ENSEMBLE(ch_ensemble_input, ch_sex_ploidy)

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
        .mix(ALIGNMENT_QC.out.stats)
        .mix(ALIGNMENT_QC.out.global_dist)
        .mix(ch_mosdepth_summary.map { meta, f -> f })
        .mix(BCFTOOLS_STATS.out.stats)
        .mix(BCFTOOLS_STATS_ENSEMBLE.out.stats)
        .collect()

    MULTIQC(ch_multiqc)
}