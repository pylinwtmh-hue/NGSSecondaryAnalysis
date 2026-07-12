/*
 * =========================================================
 * gCNV & CNVkit Panel of Normals (PON) Modules
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin
 * Licensed under the MIT License
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

// =========================================================
// GATK: Preprocess & Annotate Intervals (動態支援 WES/WGS)
// =========================================================
process PREP_GATK_INTERVALS {
    label 'process_low'

    input:
    path fasta
    path fasta_fai
    path fasta_dict
    path targets           // WES 放 capture BED；WGS 可以放 primary chr BED
    path blacklist_bed     // 黑名單 BED 檔 (Centromere, PAR 等)
    path mappability_bed
    path mappability_bed_tbi
    path segdup_bed
    path segdup_bed_tbi

    output:
    path "preprocessed.interval_list", emit: preprocessed
    path "annotated.tsv",              emit: annotated

    script:
    // 動態判斷：如果是 WES 就 padding 250 且不切 bin；WGS 則切 1000bp 且不 padding
    def bin_length = params.seq_type == "WES" ? 0 : 1000
    def padding    = params.seq_type == "WES" ? 250 : 0
    
    // 如果有提供 blacklist，就加上 -XL 參數
    def exclude_cmd = blacklist_bed.name != 'NO_FILE' ? "-XL ${blacklist_bed}" : ""

    """
    gatk PreprocessIntervals \
        -R ${fasta} \
        -L ${targets} \
        ${exclude_cmd} \
        --bin-length ${bin_length} \
        --padding ${padding} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O preprocessed.interval_list

    gatk AnnotateIntervals \
        -R ${fasta} \
        -L preprocessed.interval_list \
        --mappability-track ${mappability_bed} \
        --segmental-duplication-track ${segdup_bed} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O annotated.tsv
    """
}

// =========================================================
// CNVkit: Prepare Target/Antitarget BEDs
// =========================================================
process PREP_CNVKIT_BEDS {
    label 'process_low'

    input:
    path fasta
    path fasta_fai
    path wes_targets

    output:
    path "targets.bed",     emit: target_bed
    path "antitargets.bed", emit: antitarget_bed

    script:
    """
    cnvkit.py target ${wes_targets} --split -o targets.bed
    cnvkit.py access ${fasta} -o access.bed
    cnvkit.py antitarget targets.bed -g access.bed -o antitargets.bed
    """
}

// =========================================================
// 收集深度 (GATK 專用)
// =========================================================
process COLLECT_GATK_COUNTS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(recal)
    path fasta
    path fasta_fai 
    path fasta_dict
    path gatk_intervals

    output:
    path "${meta.id}.counts.hdf5", emit: gatk_counts

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-1}g" CollectReadCounts \
        -I ${bam} \
        -R ${fasta} \
        -L ${gatk_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${meta.id}.counts.hdf5
    """
}

// =========================================================
// 收集深度 (CNVkit 專用)
// =========================================================
process COLLECT_CNVKIT_COV {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(recal)
    path cnvkit_target_bed
    path cnvkit_antitarget_bed

    output:
    path "${meta.id}.targetcoverage.cnn", emit: cnvkit_t_cov
    path "${meta.id}.antitargetcoverage.cnn", emit: cnvkit_a_cov

    script:
    """
    cnvkit.py coverage ${bam} ${cnvkit_target_bed} -o ${meta.id}.targetcoverage.cnn
    cnvkit.py coverage ${bam} ${cnvkit_antitarget_bed} -o ${meta.id}.antitargetcoverage.cnn
    """
}

// =========================================================
// CNVkit: 建立 Pooled Reference (.cnn)
// =========================================================
process CNVKIT_REFERENCE {
    label 'process_high'
    publishDir "${params.pon_out_dir}/cnvkit_reference", mode: 'copy'

    input:
    path fasta
    path fasta_fai
    path t_covs
    path a_covs

    output:
    path "cnvkit_pooled_reference.cnn", emit: reference

    script:
    """
    cnvkit.py reference *.targetcoverage.cnn *.antitargetcoverage.cnn \
        -f ${fasta} \
        -o cnvkit_pooled_reference.cnn
    """
}

// =========================================================
// GATK gCNV: Filter Intervals (結合 Mappability 與 SegDup 過濾)
// =========================================================
process FILTER_INTERVALS {
    label 'process_medium'
    publishDir "${params.pon_out_dir}", mode: 'copy'

    input:
    path counts
    path preprocessed_intervals
    path annotated_intervals

    output:
    path "filtered.interval_list", emit: intervals

    script:
    def counts_args = counts.collect { "-I ${it}" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-1}g" FilterIntervals \
        -L ${preprocessed_intervals} \
        --annotated-intervals ${annotated_intervals} \
        ${counts_args} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --minimum-gc-content 0.1 \
        --maximum-gc-content 0.9 \
        --minimum-mappability 0.9 \
        --maximum-mappability 1.0 \
        --minimum-segmental-duplication-content 0.0 \
        --maximum-segmental-duplication-content 0.5 \
        --low-count-filter-count-threshold 5 \
        --low-count-filter-percentage-of-samples 90.0 \
        --extreme-count-filter-minimum-percentile 1.0 \
        --extreme-count-filter-maximum-percentile 99.0 \
        --extreme-count-filter-percentage-of-samples 90.0 \
        -O filtered.interval_list
    """
}

// =========================================================
// GATK gCNV: Ploidy Cohort Model (自動偵測為 Cohort 模式)
// =========================================================
process PLOIDY_COHORT {
    label 'process_high'
    publishDir "${params.pon_out_dir}/gcnv_model", mode: 'copy'

    input:
    path counts
    path filtered_intervals
    path contig_ploidy_priors

    output:
    path "ploidy_model", emit: ploidy_model
    path "ploidy_calls", emit: ploidy_calls

    script:
    def counts_args = counts.collect { "-I ${it}" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-2}g" DetermineGermlineContigPloidy \
        -L ${filtered_intervals} \
        ${counts_args} \
        --contig-ploidy-priors ${contig_ploidy_priors} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ploidy_calls \
        --output-prefix cohort \
        --verbosity INFO
        
    mv ploidy_calls/cohort-model ploidy_model
    """
}

// =========================================================
// GATK gCNV: Scatter Intervals
// =========================================================
process SCATTER_INTERVALS {
    label 'process_low'

    input:
    path filtered_intervals

    output:
    path "scatter_dir/*/*.interval_list", emit: scattered_lists

    script:
    """
    gatk IntervalListTools \
        --INPUT ${filtered_intervals} \
        --SUBDIVISION_MODE INTERVAL_COUNT \
        --SCATTER_CONTENT 5000 \
        --OUTPUT scatter_dir
    """
}

// =========================================================
// GATK gCNV: Cohort Model (散佈執行)
// 超參數由 nextflow_pon.config 提供（貼齊 Broad germline CNV WDL 預設）：
//   gcnv_p_alt / gcnv_cnv_coherence / gcnv_class_coherence / gcnv_p_active
// ⚠️ 改動後必須重跑 main_pon.nf 重建 model；敏感度取捨見 CLAUDE.md。
// =========================================================
process GCNV_COHORT {
    label 'process_high'
    publishDir "${params.pon_out_dir}/gcnv_model/shards", mode: 'copy'

    input:
    path interval_shard
    path counts
    path annotated_intervals
    path ploidy_calls

    output:
    path "gcnv_model_shard_*", emit: model_shard
    path "gcnv_calls_shard_*", emit: call_shard

    script:
    def shard_name = interval_shard.baseName
    def counts_args = counts.collect { "-I ${it}" }.join(" \\\n        ")
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g" GermlineCNVCaller \
        --run-mode COHORT \
        -L ${interval_shard} \
        ${counts_args} \
        --contig-ploidy-calls ${ploidy_calls}/cohort-calls \
        --annotated-intervals ${annotated_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --cnv-coherence-length ${params.gcnv_cnv_coherence} \
        --class-coherence-length ${params.gcnv_class_coherence} \
        --p-alt ${params.gcnv_p_alt} \
        --p-active ${params.gcnv_p_active} \
        -O gcnv_model_shard_${shard_name} \
        --output-prefix cohort_${shard_name} \
        --verbosity INFO
        
    mv gcnv_model_shard_${shard_name}/cohort_${shard_name}-calls gcnv_calls_shard_${shard_name}
    """
}