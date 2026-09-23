#!/usr/bin/awk -f
# =========================================================
# WGS/WES Germline Analysis Pipeline - haploid_het.awk
# =========================================================
# Author   : Po-Yu Lin (林伯昱)
# Institute: Department of Neurology and
#            Department of Genomic Medicine,
#            National Cheng Kung University Hospital
# Contact  : p88124019@gs.ncku.edu.tw
#
# Copyright (c) 2026, Po-Yu Lin (林伯昱)
# Licensed under the GNU General Public License v3.0
#
# DISCLAIMER: Provided "as is" without warranty. Users are solely responsible
# for validating and interpreting all results.
# =========================================================
# scripts/haploid_het.awk
# =======================
# 男性單倍體區的 het：在 bcftools +fixploidy 之前先處理（BCFTOOLS_ENSEMBLE 呼叫）。
#
# 為什麼要做
# ----------
# DeepVariant / HaplotypeCaller 一律用 diploid call，男性 chrX 非 PAR 與 chrY 也可能叫出 het。
# +fixploidy 把 GT 改成單套時「只留第一個 allele」（bcftools plugins/fixploidy.c）：
#   0/1、0|1 → 0（變 REF，報告看不到）；1|0 → 1（hemizygous）。
# 開了 phasing 之後，同樣是 het，結果取決於 whatshap 任意定的 phase 方向。
# 男性的 het 多半是比對假象（X/Y 高度相似區、片段重複、偽基因），但也可能是 47,XXY 或
# 體細胞嵌合 —— X-linked 顯性、男性通常致死的疾病（IKBKG、MECP2、CDKL5、PORCN、OFD1…）存活的
# 男性病人常是嵌合；PCDH19 則是 hemizygous 男性通常不發病、嵌合男性才發病。
# 所以既不能一律藏掉，也不能讓 phase 方向決定。
#
# 做法（只在 -v SEX=M、且 ploidy 檔中該區間 ploidy=1 的位置；chrM 不處理）
# --------------------------------------------------------------------
#   chrY：het → ./.（不進報告；chrY 的 het 幾乎都是比對假象）
#   其他（chrX 非 PAR）：het 0/k → k/k（保留 / 或 |），INFO 加 HAPLOID_HET=<caller,...>
#     → +fixploidy 截斷後一律是 k，call 會進報告；AD/VAF 等其他 FORMAT 欄位原封不動，
#       三級報告的 HAPLOID_HET 欄提示需人工複核（嵌合／XXY／比對假象）。
#     兩個不同 ALT 的 het（j/k，j、k 都 > 0）：只加標記，GT 不動（實務上不會出現，
#     因為 ensemble 的每個 caller 在 merge 前都已拆成 biallelic）。
#   半缺失（./1）、hom（1/1）、REF（0/0）、缺失（./.）：不動。
#
# 用法
# ----
#   awk -v SEX=M -f haploid_het.awk <ploidy.txt> <in.vcf | ->  > out.vcf
#     ploidy.txt：CHROM FROM TO SEX PLOIDY（與 +fixploidy -p 同一份，空白或 tab 分隔）
#     sample 名稱結尾的 _DV / _HC 當作 caller 名稱（BCFTOOLS_ENSEMBLE 的命名方式）
#   stderr：各 caller 在 chrX 轉成 ALT 的筆數、chrY 設成 missing 的筆數。
#   只用 POSIX awk 語法（容器裡是 busybox awk；mawk / gawk 也測過）。
# =========================================================

BEGIN { FS = "\t"; OFS = "\t"; nreg = 0 }

# ---- 第一個檔：ploidy 定義 -------------------------------------------------
FILENAME == ARGV[1] {
    n = split($0, a, /[ \t]+/)
    if (n >= 5 && a[1] !~ /^#/ && a[4] == SEX && a[5] == "1" \
        && a[1] != "chrM" && a[1] != "MT" && a[1] != "M") {
        nreg++
        rchr[nreg] = a[1]; rlo[nreg] = a[2] + 0; rhi[nreg] = a[3] + 0
    }
    next
}

# ---- VCF header -------------------------------------------------------------
/^##/ { print; next }
/^#CHROM/ {
    print "##INFO=<ID=HAPLOID_HET,Number=.,Type=String,Description=\"Male haploid region (chrX non-PAR): caller(s) whose original GT was heterozygous. GT was set to the ALT allele before +fixploidy so the call is kept; possible mosaicism, XXY or mis-mapping - review manually\">"
    print
    for (i = 10; i <= NF; i++) { c = $i; sub(/^.*_/, "", c); caller[i] = c }
    next
}

# ---- 資料列 -----------------------------------------------------------------
{
    chr = $1; pos = $2 + 0
    inreg = 0
    for (r = 1; r <= nreg; r++)
        if (chr == rchr[r] && pos >= rlo[r] && pos <= rhi[r]) { inreg = 1; break }
    if (!inreg || $9 !~ /^GT(:|$)/) { print; next }     # GT 依 VCF 規範必須是第一個 FORMAT 欄

    isY = (chr == "chrY" || chr == "Y")
    flag = ""
    for (i = 10; i <= NF; i++) {
        s = $i
        p = index(s, ":")
        if (p > 0) { gt = substr(s, 1, p - 1); rest = substr(s, p) }
        else       { gt = s; rest = "" }
        sep = (index(gt, "|") > 0) ? "|" : "/"
        na = split(gt, al, /[\/|]/)
        if (na != 2 || al[1] == "." || al[2] == "." || al[1] == al[2]) continue   # 不是 het

        if (isY) {                                   # chrY：het → missing
            $i = "./." rest
            ny[caller[i]]++
            continue
        }
        if (al[1] == "0")      k = al[2]
        else if (al[2] == "0") k = al[1]
        else                   k = ""                # j/k：只標記
        if (k != "") { $i = k sep k rest; nx[caller[i]]++ }
        else         { nxm[caller[i]]++ }
        flag = (flag == "") ? caller[i] : flag "," caller[i]
    }
    if (flag != "") $8 = (($8 == "." || $8 == "") ? "" : $8 ";") "HAPLOID_HET=" flag
    print
}

END {
    msg = "[haploid_het] SEX=" SEX " haploid_regions=" nreg
    for (c in nx)  msg = msg " chrX_het_to_alt_" c "=" nx[c]
    for (c in nxm) msg = msg " chrX_multi_alt_het_" c "=" nxm[c]
    for (c in ny)  msg = msg " chrY_het_to_missing_" c "=" ny[c]
    print msg > "/dev/stderr"
}
