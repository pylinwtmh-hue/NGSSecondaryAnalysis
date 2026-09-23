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
    // 性別感知倍體定義（單一真相來源；params.sex_ploidy_file，staged）
    path ploidy_file
    // 男性單倍體區 het 的前處理（scripts/haploid_het.awk，staged → 內容改了 -resume 會重跑）
    path haploid_awk

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
    #
    # ⚠️ DV arm 拆完後只留「DV 真的有 ALT call」的紀錄（-i 'GT="alt"'），再進 merge。
    #   DeepVariant 的 VCF 會保留它考慮過但否決的候選（FILTER=RefCall，GT ./. 或 0/0）。
    #   若讓它們進 merge --merge all，只要 HC 在同一 POS 有不同 allele，兩者就會被併成
    #   一筆多等位紀錄 —— ALT1 是 DV 否決的候選、ALT2 是 HC 真正的 call —— 而且合併後的
    #   FILTER 會變成 DV 的 RefCall。三級 norm 再拆開時，被否決的那個 allele 會變成一筆
    #   「兩邊都沒 call」的獨立紀錄。
    #   實例（SUZ12 chr17:31998950）：DV 把 delinsTT 的三個片段（GAAA>GAA、952 A>T、
    #   953 A>T）全判 RefCall；HC 經 combine_phased 正確合成 GAAA>GTT（0|1）。merge 後
    #   變成 GAAA  GAA,GTT、FILTER=RefCall，三級拆開後報告多出錯誤的 c.2170del，且正確的
    #   delinsTT 帶著 DV 的 AD「10,0」（DV 從未評估過這個 allele）。
    #   先 norm 再 filter 的順序是必要的：DV 的多等位 0/2 拆開後是 0/0（被否決的 ALT）+
    #   0/1（真的 call），這樣才能只丟掉被否決的那個 allele。
    #   GT="alt" 的實測語義（bcftools）：保留 0/1、0|1、1/1、1|1、1/0；丟掉 ./.、0/0、0|0、
    #   ./0，以及半缺失的 ./1、1/. —— 與三級 add_callers_tag.is_called()（任一 allele 缺失
    #   即視為沒 call）一致。DV 本身不產生半缺失 GT。
    #   RefCall 仍完整保留在已發布的 <id>.deepvariant.vcf.gz（BGZIP_VCF_DV），可供稽核；
    #   CNVkit 的 b-allele 輸入讀的也是那份原始 DV VCF，不受影響。
    #   HC 不需要同樣處理：HC 的 VCF 模式只輸出有 ALT 的位點。
    #   這裡擋不到 combine_phased.py 合成的紀錄（它的 GT 是重建的、帶 ALT）。舊版 combine
    #   會把被否決的較寬候選挑成 anchor，合成紀錄因此帶著 FILTER=RefCall 與它的 AD/VAF 留下來
    #   （VAL55：28,050 筆）；已在 combine_phased.py 修正（沒有 ALT 的紀錄不參與叢集）。
    #   不用 pipe：本 pipeline 的 shell 是 bash -ue（沒有 pipefail），norm 若中途失敗，
    #   接在後面的 view 仍可能以 0 結束並寫出截斷的檔案。拆成兩步，各自被 -e 檢查。
    bcftools norm -m -any fx_dv.vcf.gz -O u -o norm_dv.bcf
    bcftools view -i 'GT="alt"' norm_dv.bcf -O z -o temp_dv.vcf.gz
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

    # B. 倍體定義檔 (Ploidy Map) 改由 config 的 params.sex_ploidy_file 傳入（單一真相來源，
    #    格式 CHROM FROM TO SEX PLOIDY，GRCh38 PAR 座標）。以 staged input 進來，見上方 input。

    # -------------------------------------------------------------
    # 3.5 男性單倍體區的 het（只對男性；在 +fixploidy 之前）
    # -------------------------------------------------------------
    # +fixploidy 把男性 chrX 非 PAR / chrY 的 GT 改成單套時「只留第一個 allele」：
    #   0/1、0|1 → 0（變 REF，報告看不到）；1|0 → 1（hemizygous）。開了 phasing 之後，
    #   同樣是 het，結果取決於 whatshap 任意定的 phase 方向（VAL55：chrX 2,607 筆 + chrY 6,049 筆
    #   被截成 REF 而消失，另有數百筆 phase 過的 het 被顯示成 hemizygous）。
    # 男性的 het 多半是比對假象，但也可能是 47,XXY 或體細胞嵌合 —— X-linked 顯性、男性通常致死的
    #   疾病（IKBKG、MECP2、CDKL5、PORCN、OFD1…）存活的男性病人常是嵌合；PCDH19 則是嵌合男性才
    #   發病。所以不能一律藏掉，也不能讓 phase 方向決定。
    # haploid_het.awk（讀同一份 ploidy 檔）：
    #   chrX 非 PAR：het 0/k → k/k，INFO 加 HAPLOID_HET=<DV,HC> → 截斷後一律保留 ALT，
    #     三級報告的 HAPLOID_HET 欄提示人工複核；AD/VAF 原封不動，判讀者可看比例。
    #   chrY：het → ./.（不進報告；chrY 的 het 幾乎都是比對假象）。chrM 不處理。
    # 用 subshell + pipefail：本 pipeline 的 shell 是 bash -ue（沒有 pipefail），直接接 pipe 的話
    #   上游失敗會被吞掉；subshell 讓 pipefail 只作用在這一段。
    FIX_IN=${prefix}.ensemble.raw.vcf.gz
    if [ "${sex}" = "M" ]; then
        command -v awk >/dev/null || { echo "[BCFTOOLS_ENSEMBLE] awk not found in container" >&2; exit 1; }
        ( set -o pipefail
          bcftools view ${prefix}.ensemble.raw.vcf.gz \\
            | awk -v SEX=M -f ${haploid_awk} ${ploidy_file} - \\
            | bcftools view -O z -o ${prefix}.ensemble.hh.vcf.gz - )
        FIX_IN=${prefix}.ensemble.hh.vcf.gz
    fi

    # -------------------------------------------------------------
    # 4. 執行 bcftools +fixploidy 進行優雅校正
    # -------------------------------------------------------------
    bcftools +fixploidy \$FIX_IN \\
        -O z -o ${prefix}.ensemble.fixed.vcf.gz \\
        -- -s sample_sex.txt -p ${ploidy_file}

    bcftools index --tbi ${prefix}.ensemble.fixed.vcf.gz

    # -------------------------------------------------------------
    # 5. 發布前 preflight：確認 ensemble 可在「不用 --force」下通過 norm -m
    #    （Number=R/A/G 欄位數正確）。壞掉就讓二級 fail loud，不把壞檔丟給三級（見回報 §8）。
    # -------------------------------------------------------------
    bcftools norm -m -any ${prefix}.ensemble.fixed.vcf.gz -O u -o /dev/null

    # -------------------------------------------------------------
    # 6. 清理所有暫存檔
    # -------------------------------------------------------------
    rm -f rename_dv.txt rename_hc.txt rn_dv.vcf.gz* rn_hc.vcf.gz* hdr_dv.txt hdr_hc.txt fx_dv.vcf.gz* fx_hc.vcf.gz* norm_dv.bcf temp_dv.vcf.gz* temp_hc.vcf.gz* sample_sex.txt ${prefix}.ensemble.raw.vcf.gz* ${prefix}.ensemble.hh.vcf.gz
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
