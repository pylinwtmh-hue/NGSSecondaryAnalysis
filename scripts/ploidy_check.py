#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=========================================================
WGS/WES Germline Analysis Pipeline - ploidy_check.py
=========================================================
Author   : Po-Yu Lin (林伯昱)
Institute: Department of Neurology and
           Department of Genomic Medicine,
           National Cheng Kung University Hospital
Contact  : p88124019@gs.ncku.edu.tw

Copyright (c) 2026, Po-Yu Lin
Licensed under the GNU General Public License v3.0

DISCLAIMER: Provided "as is" without warranty. Users are solely responsible
for validating and interpreting all results.
=========================================================
scripts/ploidy_check.py
=======================

從 mosdepth 的 `*.summary.txt` 推「性別 + 每條染色體的相對 ploidy」，做兩件事：
  1. sex 防呆：把「資料推得的性別」和 samplesheet 宣告的 sex 比對，不符 → WARN。
  2. aneuploidy 提示：每條體染色體的 normalized coverage（相對體染色體中位數），
     偏離 1.0 太多（如三體 ≈1.5、單體 ≈0.5）→ WARN「考慮手動 -ploidy N 重跑」。

**warn-only**：只印警示 + 寫進 QC 檔，永不改 ploidy、永不讓 pipeline 失敗（exit 0）。
最保守：把判斷權留給人。

輸出（對齊 DRAGEN 的 *.ploidy.vcf 風格，交由 Nextflow bgzip 成 .gz）：
  - <sample>.ploidy.vcf ：每條 contig 一列，FORMAT=DC:NDC（DC=mean depth、NDC=normalized），
    header 帶 ##estimatedSexKaryotype / ##declaredSexKaryotype；疑似非整倍體的 contig
    FILTER=SUSPECT。
  - <sample>.ploidy_qc.txt：人可讀摘要 + WARNINGS 區塊。

相依：只用 Python 標準庫。

mosdepth summary 欄位：chrom  length  bases  mean  min  max
  - WGS：用整條染色體 mean；--by autosome_bed 產生的 *_region 列（autosome 專屬）優先，
    X/Y 無 region 列 → 回退整條 mean（WGS 全基因體覆蓋，兩者接近）。
  - WES：整條 mean 會被 off-target 稀釋，故優先用 *_region（on-target）mean；X/Y 若無
    target 則不可靠 → 於 QC 檔註記（WES 另有 gCNV 的 ploidy 可交叉驗證）。
"""

import argparse
import statistics
import sys

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
# 分類門檻（保守、warn-only）
X_MALE_MAX = 0.75     # X normalized < 此 → 傾向單套（男）
Y_PRESENT_MIN = 0.15  # Y normalized ≥ 此 → 有 Y
ANEUPLOIDY_LO = 0.75  # 體染色體 normalized < 此 或
ANEUPLOIDY_HI = 1.25  #                    > 此 → 疑似非整倍體


def parse_mosdepth_summary(path):
    """回傳 {chrom: mean}。每條 contig 優先取 *_region（on-target），否則整條 mean。"""
    whole, region, length = {}, {}, {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        ci = idx.get("chrom", 0)
        li = idx.get("length", 1)
        mi = idx.get("mean", 3)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= mi:
                continue
            name = p[ci]
            if name in ("total", "total_region"):
                continue
            try:
                mean = float(p[mi])
            except ValueError:
                continue
            if name.endswith("_region"):
                base = name[:-len("_region")]
                region[base] = mean
            else:
                whole[name] = mean
                try:
                    length[name] = int(p[li])
                except (ValueError, IndexError):
                    length[name] = 0
    means = {}
    for c in set(list(whole) + list(region)):
        means[c] = region[c] if c in region else whole.get(c)
    return means, length


def normalize_declared_sex(raw):
    """samplesheet 的 sex（male/female/unknown 或 M/F/XY/XX/1/2）→ XX / XY / unknown。"""
    s = (raw or "").strip().lower()
    if s in ("male", "m", "xy", "1"):
        return "XY"
    if s in ("female", "f", "xx", "2"):
        return "XX"
    return "unknown"


def infer_karyotype(means, baseline):
    """回傳 (karyotype_str, x_ratio, y_ratio)。"""
    x = means.get("chrX")
    y = means.get("chrY")
    x_ratio = (x / baseline) if (x is not None and baseline) else None
    y_ratio = (y / baseline) if (y is not None and baseline) else None
    if x_ratio is None:
        return "unknown", x_ratio, y_ratio
    has_y = y_ratio is not None and y_ratio >= Y_PRESENT_MIN
    x_single = x_ratio < X_MALE_MAX
    if x_single and has_y:
        kary = "XY"
    elif not x_single and not has_y:
        kary = "XX"
    elif not x_single and has_y:
        kary = "XXY?"     # X 兩套 + 有 Y（Klinefelter-like）
    elif x_single and not has_y:
        kary = "X0?"      # 單 X、無 Y（Turner-like）
    else:
        kary = "ambiguous"
    return kary, x_ratio, y_ratio


def aneuploidy_flags(means, baseline):
    """回傳 [(chrom, ndc), ...]，體染色體 normalized coverage 偏離 1.0 太多者。"""
    out = []
    if not baseline:
        return out
    for c in AUTOSOMES:
        m = means.get(c)
        if m is None:
            continue
        ndc = m / baseline
        if ndc < ANEUPLOIDY_LO or ndc > ANEUPLOIDY_HI:
            out.append((c, ndc))
    return out


def analyze(means, length, declared_sex, seq_type="WGS"):
    """回傳 dict 結果（含 warnings）。純函式、好測試。"""
    auto_means = [means[c] for c in AUTOSOMES if means.get(c) is not None]
    baseline = statistics.median(auto_means) if auto_means else 0.0
    est_kary, x_ratio, y_ratio = infer_karyotype(means, baseline)
    declared_kary = normalize_declared_sex(declared_sex)
    flags = aneuploidy_flags(means, baseline)

    warnings = []
    if declared_kary != "unknown" and est_kary not in ("unknown", "ambiguous") \
            and est_kary != declared_kary:
        warnings.append("SEX MISMATCH：samplesheet 宣告 %s，資料推得 %s"
                        "（可能 sample swap 或性染色體 aneuploidy，請人工確認）"
                        % (declared_kary, est_kary))
    if est_kary in ("XXY?", "X0?", "ambiguous"):
        warnings.append("性染色體核型異常/不明：推得 %s（X_ratio=%s, Y_ratio=%s）"
                        % (est_kary, _fmt(x_ratio), _fmt(y_ratio)))
    for c, ndc in flags:
        warnings.append("%s normalized coverage=%.2f → 疑似非整倍體，"
                        "本染色體 SNV 基因型可能不準，考慮手動 -ploidy N 重跑" % (c, ndc))
    if seq_type == "WES":
        warnings.append("注意：WES 以 target 覆蓋估計，性別/ploidy 判斷較粗；"
                        "aneuploidy 建議以 gCNV 交叉驗證")

    return {
        "baseline": baseline,
        "estimated_karyotype": est_kary,
        "declared_karyotype": declared_kary,
        "x_ratio": x_ratio,
        "y_ratio": y_ratio,
        "aneuploidy": flags,
        "warnings": warnings,
        "means": means,
        "length": length,
    }


def _fmt(x):
    return "%.2f" % x if isinstance(x, float) else "NA"


# ─────────────────────────────────────────────────────────────
# 輸出
# ─────────────────────────────────────────────────────────────
def write_vcf(path, sample, res, seq_type):
    baseline = res["baseline"]
    means, length = res["means"], res["length"]
    flagged = {c for c, _ in res["aneuploidy"]}
    # 依標準染色體排序，其餘照名稱
    order = AUTOSOMES + ["chrX", "chrY", "chrM"]
    chroms = [c for c in order if c in means] + \
             sorted(c for c in means if c not in order)
    with open(path, "w") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write("##source=NCKUH_PLOIDY_MOSDEPTH\n")
        w.write("##seqType=%s\n" % seq_type)
        w.write("##estimatedSexKaryotype=%s\n" % res["estimated_karyotype"])
        w.write("##declaredSexKaryotype=%s\n" % res["declared_karyotype"])
        w.write('##FILTER=<ID=SUSPECT,Description="Normalized coverage suggests '
                'possible aneuploidy for this contig">\n')
        w.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Contig end">\n')
        w.write('##FORMAT=<ID=DC,Number=1,Type=Float,Description="Mean depth of coverage '
                '(mosdepth)">\n')
        w.write('##FORMAT=<ID=NDC,Number=1,Type=Float,Description="Normalized depth of '
                'coverage relative to autosomal median (~1.0 diploid, ~0.5 haploid)">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        for c in chroms:
            m = means[c]
            ndc = (m / baseline) if baseline else 0.0
            filt = "SUSPECT" if c in flagged else "PASS"
            end = length.get(c, 0)
            w.write("%s\t1\t.\tN\t.\t.\t%s\tEND=%d\tDC:NDC\t%.4f:%.4f\n"
                    % (c, filt, end, m, ndc))


def write_qc(path, sample, res, seq_type):
    with open(path, "w") as w:
        w.write("# Ploidy QC — %s (%s)\n" % (sample, seq_type))
        w.write("declared_sex_karyotype : %s\n" % res["declared_karyotype"])
        w.write("estimated_sex_karyotype: %s\n" % res["estimated_karyotype"])
        sex_ok = (res["declared_karyotype"] == "unknown"
                  or res["estimated_karyotype"] in ("unknown", "ambiguous")
                  or res["declared_karyotype"] == res["estimated_karyotype"])
        w.write("sex_check              : %s\n" % ("OK" if sex_ok else "MISMATCH"))
        w.write("autosomal_baseline_mean: %.4f\n" % res["baseline"])
        w.write("chrX_ratio             : %s\n" % _fmt(res["x_ratio"]))
        w.write("chrY_ratio             : %s\n" % _fmt(res["y_ratio"]))
        w.write("\n--- per-chromosome normalized coverage (NDC) ---\n")
        baseline = res["baseline"]
        order = AUTOSOMES + ["chrX", "chrY", "chrM"]
        for c in [x for x in order if x in res["means"]]:
            ndc = (res["means"][c] / baseline) if baseline else 0.0
            w.write("%-6s %.3f\n" % (c, ndc))
        w.write("\n--- WARNINGS ---\n")
        if res["warnings"]:
            for msg in res["warnings"]:
                w.write("WARN: %s\n" % msg)
        else:
            w.write("(none)\n")


def main():
    ap = argparse.ArgumentParser(description="Infer sex/ploidy from mosdepth summary; warn-only.")
    ap.add_argument("--summary", required=True, help="mosdepth *.summary.txt")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--declared-sex", default="unknown",
                    help="samplesheet sex (male/female/unknown/…)")
    ap.add_argument("--seq-type", default="WGS", help="WGS or WES")
    ap.add_argument("--out-vcf", required=True, help="output ploidy VCF (uncompressed)")
    ap.add_argument("--out-qc", required=True, help="output human-readable QC txt")
    a = ap.parse_args()

    means, length = parse_mosdepth_summary(a.summary)
    res = analyze(means, length, a.declared_sex, a.seq_type)
    write_vcf(a.out_vcf, a.sample, res, a.seq_type)
    write_qc(a.out_qc, a.sample, res, a.seq_type)

    # warn-only：印到 stderr（→ nextflow log），永不 fail
    for msg in res["warnings"]:
        sys.stderr.write("[ploidy_check] WARN: %s\n" % msg)
    sys.stderr.write("[ploidy_check] %s declared=%s estimated=%s baseline=%.2f\n"
                     % (a.sample, res["declared_karyotype"],
                        res["estimated_karyotype"], res["baseline"]))


if __name__ == "__main__":
    main()
