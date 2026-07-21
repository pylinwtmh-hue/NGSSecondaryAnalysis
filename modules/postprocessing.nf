/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - PostProcessing Module
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

// Compress the vcf (since parabrick 4.4.0 cannot output vcf.gz )
process BGZIP_VCF {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    input:
    tuple val(meta), path(vcf)  // 未壓縮的 .vcf

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf

    script:
    """
    bgzip -@ ${task.cpus} ${vcf}
    tabix -p vcf ${vcf}.gz
    """
}

// BCFtools Ensemble（合併 DeepVariant + HaplotypeCaller）
// 標註每個 variant 的來源（SOURCE tag）
// WGS：DV + VQSR；WES：DV + HaplotypeCaller（直接）
// 實測：RAM 66MB，CPU 514%（~5 cores）
process BCFTOOLS_ENSEMBLE {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.out_dir}/${meta.id}/04_snv_indel", mode: 'copy'

    // INPUT:
    //   dv_vcf - DeepVariant VCF（Lane 1 輸出）
    //   hc_vcf - HaplotypeCaller VCF（WGS: VQSR 後；WES: 直接輸出）
    input:
    // main.nf 用 .join(by: 0) 確保同一樣本配對後合併成單一 tuple 傳入
    // 避免多樣本非同步完成時 DV/HC 跨樣本錯配的 bug
    tuple val(meta),
        path(dv_vcf), path(dv_tbi),
        path(hc_vcf), path(hc_tbi)
        
    // OUTPUT:
    //   vcf - 兩個 caller 合併且校正過 ploidy後的 ensemble VCF（含 SOURCE INFO tag，2 samples）
    output:
    tuple val(meta), path("*.ensemble.fixed.vcf.gz"), path("*.ensemble.fixed.vcf.gz.tbi"), emit: vcf
    
    script:
    def prefix = "${meta.id}"

    // 嚴格的性別字串轉換 (轉為 bcftools 規定的 M / F)
    def sex = "F" // 預設女性最安全 (二倍體)
    if (meta.sex && meta.sex != 'unknown') {
        def s = meta.sex.toString().toUpperCase()
        if (s == 'M' || s == 'MALE') { 
            sex = 'M' 
        }
    }
    """
    # fixploidy.so 存在於容器的 /usr/local/libexec/bcftools/，
    # 但 BCFTOOLS_PLUGINS 環境變數預設未設定，需手動指定
    export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools
    # -------------------------------------------------------------
    # 1. 修改 Sample ID (加上 _DV 和 _HC 後綴)
    # -------------------------------------------------------------
    echo "${prefix} ${prefix}_DV" > rename_dv.txt
    echo "${prefix} ${prefix}_HC" > rename_hc.txt

    bcftools reheader -s rename_dv.txt ${dv_vcf} -o rn_dv.vcf.gz
    bcftools reheader -s rename_hc.txt ${hc_vcf} -o rn_hc.vcf.gz

    # -------------------------------------------------------------
    # 2. 統一 FORMAT/AD header 為 Number=R → 各自拆 biallelic → 聯集合併 (Union)
    # -------------------------------------------------------------
    # 根因：DeepVariant 與 HaplotypeCaller 對 FORMAT/AD 的 header Number 定義「不一致」
    #   （bcftools 警告 "combine AD tag definitions of different lengths"）。只要有一邊不是
    #   Number=R，bcftools norm/merge 就無法把 AD 依 allele 正確拆分/重排：
    #     - 拆 multiallelic 時，非 R 的 AD 不會被 re-size → biallelic 卻帶多個 AD 值
    #       （NA12878 chr1:111241360：2 alleles 卻 3 個 AD → merge 失敗）；
    #     - 直接合併時 AD 沒依新 ALT union 補齊 → 三級 norm 報 "wrong number of fields"
    #       （VAL-55 chr1:83829）。
    #   解法：合併前先把兩邊 header 的 AD 強制成 Number=R（AD 本就是 per-allele），之後
    #   norm -m -any 才會正確 re-size、merge 也不再衝突。sed 對 ID=AD 那行不論原本
    #   Number 是 . / 數字 / R 一律改 R（已是 R 則無副作用）。
    # 一併把 PL 補成 Number=G（同理，PL 本就是 per-genotype；DV/HC 若對 PL 也定義不一致，
    # 會在 AD 修好後換 PL 報同類錯。已是 G / 無 PL 行則無副作用）。
    bcftools view -h rn_dv.vcf.gz \\
        | sed 's/##FORMAT=<ID=AD,Number=[^,]*,/##FORMAT=<ID=AD,Number=R,/' \\
        | sed 's/##FORMAT=<ID=PL,Number=[^,]*,/##FORMAT=<ID=PL,Number=G,/' > hdr_dv.txt
    bcftools reheader -h hdr_dv.txt rn_dv.vcf.gz -o fx_dv.vcf.gz
    bcftools view -h rn_hc.vcf.gz \\
        | sed 's/##FORMAT=<ID=AD,Number=[^,]*,/##FORMAT=<ID=AD,Number=R,/' \\
        | sed 's/##FORMAT=<ID=PL,Number=[^,]*,/##FORMAT=<ID=PL,Number=G,/' > hdr_hc.txt
    bcftools reheader -h hdr_hc.txt rn_hc.vcf.gz -o fx_hc.vcf.gz

    # 各自拆成 biallelic（AD 已是 Number=R，會被正確 re-size），再走 bcftools 標準的
    # biallelic→multiallelic 聯集路徑（Number=R/A/G 正確處理；某 caller 缺的 allele 補 '.'）。
    # phasing 開啟時，combine_phased.py 產生的 MNV / 1|2 記錄也在此一併拆開。
    bcftools norm -m -any fx_dv.vcf.gz -O z -o temp_dv.vcf.gz
    bcftools index --tbi temp_dv.vcf.gz
    bcftools norm -m -any fx_hc.vcf.gz -O z -o temp_hc.vcf.gz
    bcftools index --tbi temp_hc.vcf.gz

    bcftools merge \\
        --merge all \\
        -O z -o ${prefix}.ensemble.raw.vcf.gz \\
        temp_dv.vcf.gz temp_hc.vcf.gz

    bcftools index --tbi ${prefix}.ensemble.raw.vcf.gz

    # -------------------------------------------------------------
    # 3. 準備 fixploidy 所需的設定檔 (適應雙樣本)
    # -------------------------------------------------------------
    # A. 建立病人的性別檔 (必須把 DV 和 HC 兩欄都指定性別)
    echo "${prefix}_DV ${sex}" > sample_sex.txt
    echo "${prefix}_HC ${sex}" >> sample_sex.txt

    # B. 建立 hg38 的倍體定義檔 (Ploidy Map)
    # 定義男性 (M) 的 chrX 非 PAR 區和 chrY 為單倍體 (1)，其餘皆為二倍體 (2)
    # (1-10000 是 N，所以從 1 開始寫也沒差，結尾精準對齊你的 2781479)
    cat <<EOF > hg38_ploidy.txt
chrX 1 2781479 M 2
chrX 2781480 155701382 M 1
chrX 155701383 156030895 M 2
chrX 156030896 156040895 M 1
chrY 1 57227415 M 1
chrX 1 156040895 F 2
chrM 1 16569 * 1
* * * * 2
EOF

    # -------------------------------------------------------------
    # 4. 執行 bcftools +fixploidy 進行優雅校正
    # -------------------------------------------------------------
    bcftools +fixploidy ${prefix}.ensemble.raw.vcf.gz \\
        -O z -o ${prefix}.ensemble.fixed.vcf.gz \\
        -- -s sample_sex.txt -p hg38_ploidy.txt

    bcftools index --tbi ${prefix}.ensemble.fixed.vcf.gz

    # -------------------------------------------------------------
    # 5. 發布前 preflight：確認 ensemble 可在「不用 --force」下通過 norm -m
    #    （Number=R/A/G 欄位數正確）。壞掉就讓二級 fail loud，不把壞檔丟給三級（見回報 §8）。
    # -------------------------------------------------------------
    bcftools norm -m -any ${prefix}.ensemble.fixed.vcf.gz -O u -o /dev/null

    # -------------------------------------------------------------
    # 6. 清理所有暫存檔
    # -------------------------------------------------------------
    rm -f rename_dv.txt rename_hc.txt rn_dv.vcf.gz* rn_hc.vcf.gz* hdr_dv.txt hdr_hc.txt fx_dv.vcf.gz* fx_hc.vcf.gz* temp_dv.vcf.gz* temp_hc.vcf.gz* sample_sex.txt hg38_ploidy.txt ${prefix}.ensemble.raw.vcf.gz*
    """
    // # -------------------------------------------------------------
    // # 方案 B：嚴格取交集 (Intersection) -> 產出 1 個 Sample 欄位的 VCF
    // # -n =2 代表只要兩個 Caller 都有的位點
    // # -w 1  代表遇到交集時，保留檔案 1 (也就是 DeepVariant) 的紀錄
    // # -------------------------------------------------------------
    // bcftools isec -p isec_dir -n =2 -w 1 -O z ${dv_vcf} ${hc_vcf}
    // mv isec_dir/0000.vcf.gz ${prefix}.ensemble.vcf.gz
    // bcftools index --tbi ${prefix}.ensemble.vcf.gz
}

// BCFtools Stats（VCF QC）
// 計算 VCF 統計量（variant count, Ti/Tv ratio, indel size distribution 等）
// 實測：RAM 29MB，CPU 182%（~2 cores）
process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    // 把獨立的 stats 報告也保留下來
    publishDir "${params.out_dir}/${meta.id}/09_postprocessing", mode: 'copy'

    // INPUT:
    //   vcf - DeepVariant VCF（用於最終 QC 統計）
    input:
    tuple val(meta), path(vcf), path(tbi)

    // OUTPUT:
    //   stats - bcftools stats 文字報告（供 MultiQC 使用）
    output:
    path "*.vcf.stats", emit: stats

    script:
    // 運用 process_low 配給的 CPU 來加速讀取
    def threads = task.cpus > 1 ? task.cpus - 1 : 1
    // 從 VCF 檔名自動產生 stats 檔名，避免兩個 caller 的 stats 撞名
    def stats_name = vcf.name.replace('.vcf.gz', '.vcf.stats')
    """
    bcftools stats --threads ${threads} ${vcf} > ${stats_name}
    """
}

// MultiQC
// 整合所有 QC 報告（fastp, samtools stats, mosdepth, bcftools stats）為單一 HTML
// 實測：RAM 136MB，CPU 13%（單核，I/O bound）
process MULTIQC {
    label 'process_low'

    publishDir "${params.out_dir}/pipeline_info", mode: 'copy'

    // INPUT:
    //   multiqc_files - 所有 QC 報告的集合（fastp JSON, samtools stats,
    //                   mosdepth dist/summary, bcftools stats）
    input:
    path multiqc_files

    // OUTPUT:
    //   report - MultiQC HTML 報告
    //   data   - MultiQC 原始數據目錄
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}
