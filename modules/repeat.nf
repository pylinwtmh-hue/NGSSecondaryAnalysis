/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Repeat Expansion Module
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

// Lane 4: GangSTR (ExpansionHunter liscence改變，改用GangSTR)
// 需要額外的 TR reference BED：https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz
//   gunzip hg38_ver17.bed.gz
//   mv hg38_ver17.bed ${ref_dir}/gangstr_hg38.bed

// ============================================================
// GANGSTR_CHROM - 單一染色體 STR genotyping（平行化用）
// 由 main.nf 展開 24 個染色體，各自平行執行
// ============================================================
process GANGSTR_CHROM {
    tag "${meta.id} ${chrom}"
    label 'process_medium'

    // 中間檔案不 publish，由 GANGSTR_MERGE 合併後再 publish

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table), val(chrom)
    path fasta
    path fasta_fai
    path str_regions

    output:
    tuple val(meta), val(chrom),
        path("${meta.id}.${chrom}.str.vcf"),         emit: vcf
    path "${meta.id}.${chrom}.str.samplestats.tab",  emit: stats

    script:
    def sex_str = "F"
    if (meta.sex && meta.sex != 'unknown') {
        def s = meta.sex.toString().toLowerCase()
        if (s == 'm' || s == 'male') { sex_str = 'M' }
    }
    def exome_args = params.seq_type == "WES" ? "--nonuniform" : ""

    """
    GangSTR \
        --bam ${bam} \
        --ref ${fasta} \
        --regions ${str_regions} \
        --out ${meta.id}.${chrom}.str \
        --bam-samps ${meta.id} \
        --samp-sex ${sex_str} \
        --chrom ${chrom} \
        ${exome_args}

    # GangSTR 在該染色體無 loci 時可能不輸出檔案，建立空檔避免 Nextflow 報錯
    [ -f ${meta.id}.${chrom}.str.vcf ] || touch ${meta.id}.${chrom}.str.vcf
    [ -f ${meta.id}.${chrom}.str.samplestats.tab ] || touch ${meta.id}.${chrom}.str.samplestats.tab
    """
}

// ============================================================
// GANGSTR_MERGE - 合併 24 個染色體的 GangSTR VCF
// 用 bcftools concat 按染色體順序合併，輸出最終結果
// ============================================================
process GANGSTR_MERGE {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/06_repeat", mode: 'copy'

    input:
    tuple val(meta), path(vcfs)   // 24 個 per-chrom VCF，已按染色體順序排列

    output:
    tuple val(meta), path("${meta.id}.str.vcf"), emit: vcf

    script:
    """
    # 過濾掉空的 VCF（某些染色體可能無 loci）
    VALID_VCFS=""
    for vcf in ${vcfs}; do
        if [ -s "\${vcf}" ] && grep -q "^#" "\${vcf}" 2>/dev/null; then
            VALID_VCFS="\${VALID_VCFS} \${vcf}"
        fi
    done

    if [ -z "\${VALID_VCFS}" ]; then
        echo "ERROR: No valid GangSTR VCF files to merge" >&2
        exit 1
    fi

    # bcftools concat 按輸入順序合併（不排序，保持染色體順序）
    bcftools concat \${VALID_VCFS} -o ${meta.id}.str.vcf
    """
}

// Lane 4: ExpansionHunter（CPU）
// Short tandem repeat（STR）genotyping
// 使用 Illumina 官方完整 variant catalog（hg38），涵蓋所有已知致病 STR 位點
// 包含但不限於：HD (HTT), SCA1-17, DM1/2 (DMPK/CNBP), FRDA (FXN),
//               C9orf72 ALS/FTD, FCMTE (MARCHF6), DRPLA, SBMA (AR),
//               CANVAS (RFC1), Fragile X (FMR1) 等
// 注意：WES 資料 STR 位點通常不在 capture region 內，結果可信度較低
// 實測：RAM 369MB，CPU ~145%（~1-2 cores）
process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/06_repeat", mode: 'copy'

    // INPUT:
    //   alignment_bundle - [BAM, BAI, recal.txt]
    //   fasta            - hg38 reference
    //   fasta_fai        - 宣告 index，Nextflow 才會把它掛載進來
    //   str_catalog      - variant catalog JSON（定義所有 STR 位點的座標與結構）
    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path str_catalog

    // OUTPUT:
    //   vcf  - STR genotype VCF（每個 locus 的 repeat 數估計）
    //   json - STR genotype JSON（詳細資訊，含 read support 圖）
    //   bam  - 可以用IGV視覺化證據檔
    output:
    tuple val(meta), path("*.str.vcf"), emit: vcf
    path "*.str.json",                  emit: json
    path "*.str_realigned.bam",         emit: bamlet

    script:
    // 嚴格的性別字串轉換 (M/F -> male/female)
    def sex_str = "female" // 保守做法，預設 female
    if (meta.sex && meta.sex != 'unknown') {
        def s = meta.sex.toString().toLowerCase()
        if (s == 'm' || s == 'male') { sex_str = 'male' }
    }
    
    """
    ExpansionHunter \
        --reads ${bam} \
        --reference ${fasta} \
        --variant-catalog ${str_catalog} \
        --output-prefix ${meta.id}.str \
        --sex ${sex_str} \
        --threads ${task.cpus}
    """
}


// ──────────────────────────────────────────────────────────────
// CALL_STR sub-workflow（Lane 4）：GangSTR（依染色體平行化 → 合併）＋ 選用 ExpansionHunter。
//   對外只吃 bam_ch；fasta / gangstr_regions / str_catalog 於內部依 params 建立。
//   GangSTR / ExpansionHunter 各自 publishDir 到 06_repeat，主流程不需再接。
// ──────────────────────────────────────────────────────────────
workflow CALL_STR {
    take:
    bam_ch      // tuple(meta, bam, bai, ...)

    main:
    ch_fasta     = file(params.fasta)
    ch_fasta_fai = file("${params.fasta}.fai")
    def gangstr_regions = params.seq_type == "WES"
        ? file(params.gangstr_regions_wes)
        : file(params.gangstr_regions_wgs)

    // 展開 24 條染色體：每個樣本 × 每條染色體 = 一個 GANGSTR_CHROM
    def chroms = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
    ch_bam_chrom = bam_ch.combine(Channel.from(chroms))
    GANGSTR_CHROM(ch_bam_chrom, ch_fasta, ch_fasta_fai, gangstr_regions)

    // 按樣本收集 24 個 VCF，依染色體順序排序後併入 GANGSTR_MERGE
    ch_gangstr_vcfs = GANGSTR_CHROM.out.vcf
        .map { meta, chrom, vcf -> [meta, chrom, vcf] }
        .groupTuple(by: 0)
        .map { meta, chroms_list, vcfs ->
            def order = (1..22).collect { "chr${it}" } + ["chrX", "chrY"]
            def sorted_vcfs = [chroms_list, vcfs].transpose()
                .sort { a, b -> order.indexOf(a[0]) <=> order.indexOf(b[0]) }
                .collect { it[1] }
            [meta, sorted_vcfs]
        }
    GANGSTR_MERGE(ch_gangstr_vcfs)

    // 選用：ExpansionHunter（--run_expansionhunter，預設關閉；非商用授權）
    if (params.run_expansionhunter) {
        EXPANSIONHUNTER(bam_ch, ch_fasta, ch_fasta_fai, file(params.str_catalog))
    }

    emit:
    vcf = GANGSTR_MERGE.out.vcf
}
