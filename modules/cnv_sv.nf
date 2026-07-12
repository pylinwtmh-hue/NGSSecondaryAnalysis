/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - CNV SV Module
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

// Lane 3a: CNVkit（WGS/WES）
// Copy number variant calling
// WGS: method=wgs，不需要 targets BED
// WES: method=hybrid，需要 capture BED 限定分析區域
// 注意：WES 建議優先使用 gCNV（更標準），CNVkit 作為補充對照
process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.out_dir}/${meta.id}/05_cnv_sv", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   meta_vcf         - 引入我們做好的 SNV VCF
    //   fasta            - hg38 reference（用於計算 GC content 等 bias 校正）
    //   cnvkit_pon       - 雖然傳入 WES 的 .cnn，但我們靠腳本決定要不要用它
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    tuple val(meta_vcf), path(vcf), path(tbi)  
    path fasta
    path cnvkit_pon

    // OUTPUT:
    //   cns   - segmented copy number calls（每個 segment 的 CN 估計）
    //   call.cns - 帶有 CN=1,2,3 絕對拷貝數的結果
    //   cnr   - per-bin copy number ratios（未 segment 的原始 CN 比值）
    //   plots - diagram + scatter plots（PDF 格式）
    output:
    tuple val(meta), path("*.aligned.sorted.cns"), emit: cns
    tuple val(meta), path("*.call.cns"),            emit: call_cns
    tuple val(meta), path("*.aligned.sorted.cnr"),  emit: cnr
    path "*.pdf",                        emit: plots

    script:
    def prefix = "${meta.id}"
    def method = params.seq_type == "WGS" ? "wgs" : "hybrid"
    
    def ref_cmd = ""
    if (params.seq_type == "WES") {
        if (cnvkit_pon.name != 'NO_FILE') {
            // 情境 1：WES 且有 PON -> 最完美狀態，全靠 PON
            ref_cmd = "-r ${cnvkit_pon}"
        } else {
            // 情境 2：WES 但缺 PON -> 退回 Flat Reference，且「必須」補上 targets 參數
            ref_cmd = "--normal --fasta ${fasta} --targets ${params.wes_targets}"
        }
    } else {
        // 情境 3：WGS -> 不需要 targets，直接用 Flat Reference
        ref_cmd = "--normal --fasta ${fasta}"
    }

    def sex_arg = ""
    if (meta.sex && (meta.sex == 'M' || meta.sex.toString().toLowerCase() == 'male')) {
        sex_arg = "-y"
    }

    """
    # -----------------------------------------------------------------
    # 1. 核心 Batch (自動切換 WES 模型 vs WGS Flat 模式)
    # -----------------------------------------------------------------
    cnvkit.py batch ${bam} \
        ${ref_cmd} \
        --method ${method} \
        --output-dir . \
        --segment-method hmm-germline \
        -p ${task.cpus}

    # -----------------------------------------------------------------
    # 2. 絕對拷貝數轉換 (Call)
    # -----------------------------------------------------------------
    cnvkit.py call ${prefix}.aligned.sorted.cns \
        -v ${vcf} \
        -m clonal \
        ${sex_arg} \
        -o ${prefix}.call.cns

    # -----------------------------------------------------------------
    # 3. 繪製臨床報告
    # -----------------------------------------------------------------
    cnvkit.py scatter ${prefix}.aligned.sorted.cnr -s ${prefix}.call.cns -v ${vcf} -o ${prefix}-scatter.pdf
    cnvkit.py diagram ${prefix}.aligned.sorted.cnr -s ${prefix}.call.cns ${sex_arg} -o ${prefix}-diagram.pdf
    """
    // 已移除 --filter cn：評鑑診斷確認它對 germline 過濾太激進、會漏掉真實 CNV。
    // 現在輸出所有 segment（含 CN=2），交由下游（三級 AnnotSV / 臨床審閱）判讀。
}

// Lane 3b: Delly（WGS/WES, MANTA商用政策改變，加入Delly）
// Delly 建議加入 exclude list（telomere/centromere 等高重複區域）
//   下載：https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
// 輸出：BCF 格式，接著由 BCFTOOLS_CONVERT_DELLY 轉成 VCF.gz

process DELLY_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path delly_excl   // exclude list TSV（可選，傳空檔案則不套用）

    output:
    tuple val(meta), path("${meta.id}.delly.bcf"), emit: bcf

    script:
    def excl_arg = delly_excl.name != 'NO_FILE' ? "--exclude ${delly_excl}" : ""

    """
    delly call \
        --genome ${fasta} \
        --outfile ${meta.id}.delly.bcf \
        ${excl_arg} \
        ${bam}
    """
}

// BCFTOOLS_CONVERT_DELLY - 將 Delly BCF 轉成 VCF.gz
process BCFTOOLS_CONVERT_DELLY {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/05_cnv_sv", mode: 'copy'

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta),
        path("${meta.id}.delly.vcf.gz"),
        path("${meta.id}.delly.vcf.gz.tbi"), emit: vcf

    script:
    """
    # 只保留 FILTER=PASS（delly call 已依 PE>=3、MAPQ>=20 標記 PASS/LowQual；
    # 單樣本不能用 delly filter -f germline，故在轉檔時直接過濾 PASS，砍掉 LowQual 噪音）
    bcftools view -f PASS ${bcf} | bcftools sort -O z -o ${meta.id}.delly.vcf.gz -
    bcftools index --tbi ${meta.id}.delly.vcf.gz
    """
}

// Lane 3b: Manta（WGS/WES）
// Structural variant（SV）calling：deletion, insertion, inversion, translocation 等
// WES 加 --exome 提高特異性；WGS 不加
process MANTA_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.out_dir}/${meta.id}/05_cnv_sv", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   fasta/fai        - hg38 reference + fai index（dict 不需要）
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai

    // OUTPUT:
    //   vcf          - diploid SV calls（主要結果）
    //   small_indels - 候選小 indels（可傳給 HaplotypeCaller 做 realignment）
    //   evidence_bams - 支持 SV 斷點的迷你 BAM 檔與索引 (供 IGV 視覺化)
    output:
    tuple val(meta),
        path("manta_results/results/variants/diploidSV.vcf.gz"),
        path("manta_results/results/variants/diploidSV.vcf.gz.tbi"), emit: vcf
    tuple val(meta),
        path("manta_results/results/variants/candidateSmallIndels.vcf.gz"),
        path("manta_results/results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: small_indels
    tuple val(meta),
        path("manta_results/results/evidence/*.bam"),
        path("manta_results/results/evidence/*.bam.bai"), optional: true, emit: evidence_bams

    script:
    def exome_flag = params.seq_type == "WES" ? "--exome" : ""
    """
    # bgzip/tabix 在 manta 容器裡位於 /usr/local/libexec/，不在預設 PATH
    export PATH="/usr/local/libexec:\$PATH"

    # -------------------------------------------------------------
    # 1. 建立 hg38 標準染色體 BED (避開 decoy 碎片導致的效能災難)
    # -------------------------------------------------------------
    cat <<EOF > hg38_main_chroms.bed
chr1\t0\t248956422
chr2\t0\t242193529
chr3\t0\t198295559
chr4\t0\t190214555
chr5\t0\t181538259
chr6\t0\t170805979
chr7\t0\t159345973
chr8\t0\t145138636
chr9\t0\t138394717
chr10\t0\t133797422
chr11\t0\t135086622
chr12\t0\t133275309
chr13\t0\t114364328
chr14\t0\t107043718
chr15\t0\t101991189
chr16\t0\t90338345
chr17\t0\t83257441
chr18\t0\t80373285
chr19\t0\t58617616
chr20\t0\t64444167
chr21\t0\t46709983
chr22\t0\t50818468
chrX\t0\t156040895
chrY\t0\t57227415
chrM\t0\t16569
EOF

    # 壓縮並建立索引 (Manta 要求 --callRegions 必須是 bgzip 壓縮且有 tabix 索引的 BED)
    bgzip hg38_main_chroms.bed
    tabix -p bed hg38_main_chroms.bed.gz

    # -------------------------------------------------------------
    # 2. 執行 Configuration
    # -------------------------------------------------------------
    configManta.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        --runDir manta_results \
        --callRegions hg38_main_chroms.bed.gz \
        --outputContig \
        --generateEvidenceBam \
        ${exome_flag}

    # -------------------------------------------------------------
    # 3. 執行 Workflow
    # -------------------------------------------------------------
    manta_results/runWorkflow.py -m local -j ${task.cpus}
    """
}

// Lane 3c: gCNV - CollectReadCounts（WES only，需 --run_gcnv true）
// 計算每個 capture region bin 的 read depth，作為 gCNV 的輸入
// 同時用於 PON 建立和 case 分析
process GATK_COLLECT_READ_COUNTS {
    tag "$meta.id"
    label 'process_medium'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   fasta/fai/dict   - hg38 reference + index
    //   intervals        - WES capture BED 轉換的 interval_list
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict
    path intervals // 這裡必須餵入 PON 產生的「已經扣掉黑名單」的 preprocessed.interval_list

    // OUTPUT:
    //   counts - per-bin read counts HDF5 檔（供 gCNV case mode 使用）
    output:
    tuple val(meta), path("${meta.id}.counts.hdf5"), emit: counts

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-1}g" CollectReadCounts \
        -I ${bam} \
        -R ${fasta} \
        -L ${intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${meta.id}.counts.hdf5
    """
}

// Lane 3c: gCNV - DetermineGermlineContigPloidy case mode（WES only）
// 根據 PON 建立時的 ploidy model，推斷當前樣本各染色體的 ploidy
// 注意：sex 若未提供（填 unknown 或空白），模型會從 chrX/chrY depth 自動推斷
//       但建議提供正確性別，避免 chrX/chrY CNV 誤判
process GATK_PLOIDY_CASE {
    tag "$meta.id"
    label 'process_medium'

    // INPUT:
    //   counts           - CollectReadCounts 輸出的 HDF5
    //   ploidy_model_dir - PON 建立時產生的 ploidy model 目錄
    input:
    tuple val(meta), path(counts)
    path ploidy_model_dir 

    // OUTPUT:
    //   ploidy_calls - 各染色體 ploidy 推斷結果目錄
    output:
    tuple val(meta), path("${meta.id}_ploidy_calls"), emit: ploidy_calls

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-1}g" DetermineGermlineContigPloidy \
        --model ${ploidy_model_dir} \
        -I ${counts} \
        -O ${meta.id}_ploidy_calls \
        --output-prefix ${meta.id}
    """
}

// Lane 3c: gCNV - GermlineCNVCaller case mode (切塊平行運算，WES only）
// 使用 PON model 對當前樣本進行 CNV calling
process GATK_GERMLINE_CNV_CASE {
    tag "${meta.id} - ${model_shard.baseName}"
    label 'process_high'

    // INPUT:
    //   counts         - CollectReadCounts 輸出的 HDF5
    //   ploidy_calls   - DetermineGermlineContigPloidy case mode 輸出
    //   model_shard    - 來自 PON 的單一 gcnv_model_shard
    input:
    tuple val(meta), path(counts), path(ploidy_calls), path(model_shard)

    // OUTPUT:
    //   gcnv_calls - GermlineCNVCaller case mode 輸出目錄
    output:
    tuple val(meta), path("gcnv_calls_${model_shard.baseName}"), emit: call_shard

    script:
    def shard_name = model_shard.baseName
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g" GermlineCNVCaller \
        --run-mode CASE \
        -I ${counts} \
        --contig-ploidy-calls ${ploidy_calls}/${meta.id}-calls \
        --model ${model_shard} \
        -O gcnv_calls_${shard_name} \
        --output-prefix ${meta.id}_${shard_name}
    """
}

// Lane 3c: gCNV - PostprocessGermlineCNVCalls（WES only）
// 將 GermlineCNVCaller 的輸出轉換為標準 VCF 格式
// =========================================================
// Lane 3c: gCNV - PostprocessGermlineCNVCalls (縫合產出 VCF)
// =========================================================
process GATK_POSTPROCESS_CNV {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.out_dir}/${meta.id}/05_cnv_sv", mode: 'copy'

    // INPUT:
    //   call_shards  - GermlineCNVCaller 該病人的所有 call_shard
    //   ploidy_calls - DetermineGermlineContigPloidy case mode 輸出
    //   model_shards - 收集 PON 的所有 model_shard
    //   fasta_dict   - hg38 sequence dictionary
    input:
    tuple val(meta), path(call_shards), path(ploidy_calls)
    path model_shards
    path fasta_dict

    // OUTPUT:
    //   vcf - gCNV calls VCF（含 genotype 和 copy number 資訊）
    //   denoisedCR.tsv - 降噪矩陣輸出
    output:
    tuple val(meta), path("${meta.id}.gcnv.vcf.gz"), path("${meta.id}.gcnv.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.denoisedCR.tsv"), emit: denoised_cr 

    script:
    // 將多個資料夾轉化為多個參數，例如：--calls-shard-path shard1 --calls-shard-path shard2 ...
    def calls_args  = call_shards.collect { "--calls-shard-path ${it}/${meta.id}_${it.baseName.replace('gcnv_calls_', '')}-calls" }.join(" ")
    def models_args = model_shards.collect { "--model-shard-path ${it}" }.join(" ")
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-1}g" PostprocessGermlineCNVCalls \
        ${calls_args} \
        --contig-ploidy-calls ${ploidy_calls}/${meta.id}-calls \
        ${models_args} \
        --allosomal-contig chrX \
        --allosomal-contig chrY \
        --output-denoised-copy-ratios ${meta.id}.denoisedCR.tsv \
        --output-genotyped-intervals ${meta.id}.gcnv_intervals.vcf.gz \
        --output-genotyped-segments ${meta.id}.gcnv.vcf.gz \
        --sequence-dictionary ${fasta_dict}

    gatk IndexFeatureFile -I ${meta.id}.gcnv.vcf.gz
    """
}