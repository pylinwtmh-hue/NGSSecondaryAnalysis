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
include { CALL_SNV }                        from './modules/call_snv'
include { CALL_CNV_SV }                     from './modules/cnv_sv'
include { CALL_STR }                        from './modules/repeat'
include { CALL_MITO }                       from './modules/call_mito'
include { MULTIQC }                         from './modules/postprocessing'
include { CALL_ROH }                        from './modules/roh'

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
    // (E) Step 4-5: SNV/indel calling + ensemble（sub-workflow）
    //   內含 DeepVariant / HaplotypeCaller / VQSR（各為第二層 sub-workflow）、
    //   PHASE_COMBINE（--run_phasing）、BCFTOOLS_ENSEMBLE、DV/ensemble stats。
    // =========================================================
    CALL_SNV(ch_bam)

    // Lane 3: CNV / SV（CNVkit(WGS) + Delly + Manta(選) + gCNV(WES)）；CNVkit 用 DeepVariant VCF
    CALL_CNV_SV(ch_bam, CALL_SNV.out.dv_vcf)

    // Lane 4: STR（GangSTR 依染色體平行化 → 合併；選用 ExpansionHunter）
    CALL_STR(ch_bam)

    // Lane 5: Mitochondria variant calling
    CALL_MITO(ch_bam)

    // Lane 6: ROH（sub-workflow：bcftools roh + AutoMap，皆選用、預設關閉；ROH 不納入評鑑）
    //   用 HaplotypeCaller raw VCF（保留 GT/AD）；flag 收在 CALL_ROH 內部。
    CALL_ROH(CALL_SNV.out.hc_vcf_raw)

    // =========================================================
    // (F) MultiQC 匯總
    // =========================================================
    ch_multiqc = Channel.empty()
        .mix(FASTP.out.json)
        .mix(PARABRICKS_FQ2BAM.out.qc_metrics)
        .mix(ALIGNMENT_QC.out.stats)
        .mix(ALIGNMENT_QC.out.global_dist)
        .mix(ch_mosdepth_summary.map { meta, f -> f })
        .mix(CALL_SNV.out.dv_stats)
        .mix(CALL_SNV.out.ensemble_stats)
        .collect()

    MULTIQC(ch_multiqc)
}