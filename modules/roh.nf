/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - ROH Module
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

// AutoMap ROH（隱性遺傳診斷輔助）
// 比 bcftools roh 更準確，專為臨床遺傳診斷設計
// 輸出 HomRegions.tsv（ROH 清單）和 HomRegions.pdf（圖形報告）
// 必須用 HaplotypeCaller VCF（含 GT + AD 欄位），不能用 VQSR 後或 DeepVariant VCF
process AUTOMAP {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/08_roh", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    path "*.HomRegions.*", emit: roh_results

    script:
    """
    # AutoMap 需要對 Resources 目錄有寫入權限（解壓縮暫存檔）
    # 容器內 /opt/AutoMap 是唯讀，先複製到 work 目錄
    cp -r /opt/AutoMap ./AutoMap_local

    # AutoMap 不支援 gzip 壓縮的 VCF，需先解壓縮
    bcftools view ${vcf} -O v -o input.vcf

    bash ./AutoMap_local/AutoMap_v1.3.sh \
        --vcf input.vcf \
        --out . \
        --genome hg38 \
        --id ${meta.id} \
        --chrX \
        --minsize 1.0  

    # AutoMap 輸出在 ./<id>/ 子目錄，移到當前目錄去掉多餘的一層
    mv ${meta.id}/* .
    rm -rf ${meta.id}
    """
    // --chrX 加入 X 染色體 ROH 分析（對近親結婚診斷很重要）
    // --minsize# 最小 ROH 長度 1Mb（過濾掉短片段噪音）
}

// BCFtools ROH（備用，不在主流程中）
process BCFTOOLS_ROH {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/08_roh", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    path "*.roh.txt", emit: roh

    script:
    """
    bcftools roh \
        --AF-dflt 0.4 \
        -O r \
        -o ${meta.id}.roh.txt \
        ${vcf}
    """
}


// ──────────────────────────────────────────────────────────────
// CALL_ROH sub-workflow（Lane 6）：獨立的 ROH 車道，兩支互不相依、皆為選用（預設關閉，
//   ROH 不納入評鑑）。輸入用 HaplotypeCaller raw VCF（保留 GT/AD；VQSR 後或 DeepVariant
//   VCF 不可用）。flag 收在 sub-workflow 內部，主流程無條件呼叫即可（與 CALL_STR gate
//   ExpansionHunter、CALL_CNV_SV gate Manta/gCNV 的作法一致）。
//     --run_roh     → bcftools roh（MIT/GPL，可商用）
//     --run_automap → AutoMap（無授權，非商用／研究）
//   兩者各自 publishDir 到 08_roh，為終端輸出，故無 emit。
// ──────────────────────────────────────────────────────────────
workflow CALL_ROH {
    take:
    hc_vcf_ch      // tuple(meta, vcf, tbi)：HaplotypeCaller raw VCF

    main:
    if (params.run_roh) {
        BCFTOOLS_ROH(hc_vcf_ch)
    }
    if (params.run_automap) {
        AUTOMAP(hc_vcf_ch)
    }
}