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

從 mosdepth 的 `*.summary.txt` 推「性別 + 每條染色體 ploidy」，做兩件事：
  1. sex 防呆：把資料推得的性別與 samplesheet 宣告的 sex 比對，不符 → WARN。
  2. aneuploidy 提示：某染色體 normalized coverage（NDC）偏離 1.0 太多 → WARN「考慮手動
     -ploidy N 重跑」。

**warn-only**：只印警示 + 寫 QC 檔，永不改 ploidy、永不讓 pipeline 失敗（exit 0）。

輸出（**與三級 DRAGEN parse_dragen_ploidy.py 統一**：同 VCF header key、同 NDC 語意）：
  - <sample>.ploidy.vcf ：每 contig 一列 FORMAT=DC:NDC:RATIO；header 帶
    ##estimatedSexKaryotype（資料推得）/ ##referenceSexKaryotype（samplesheet 宣告）；
    疑似非整倍體的 contig FILTER=SUSPECT。交由 Nextflow 壓成 .gz。
  - <sample>.ploidy_qc.txt：人可讀摘要 + WARNINGS。

⚠️ NDC 語意（對齊 DRAGEN）：NDC = 觀測 ÷「**估計核型下的期望**」，已對性別正規化：
    - autosome 期望 = 1.0；chrX/chrY 期望依估計核型（男 chrX 期望 0.5）。
    - 正常樣本每 contig（含 chrX/chrY）NDC ≈ 1.0；偏離 1.0 = 非預期（如三體 ≈1.5）。
  另存 RATIO = 相對體染色體中位數的**原始**覆蓋比（男 chrX ≈ 0.5），保留性別/劑量證據。
  （aneuploidy 判定用 NDC 偏離 1.0；性染色體非整倍體另由核型判定 XXY?/X0? 反映。）

mosdepth summary 欄位：chrom length bases mean min max。每 contig 優先取 *_region
（on-target）mean 且該值 >0，否則回退整條 mean。WES 用捕獲 BED，chrX 有 region → 用 region；
WGS 用 autosome BED，chrX/chrY/chrM 的 *_region=0 代表「無區間」而非深度 0 → 回退 whole
（否則男生性染色體會被誤判成 X0?/MISMATCH）。輸出只列主要 contig（chr1-22,X,Y,M）。
相依：只用 Python 標準庫。
"""

import argparse
import statistics
import sys

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
MITO = {"chrM", "chrMT", "MT", "M"}
# 分類門檻（保守、warn-only）
X_MALE_MAX = 0.75     # X 原始比 < 此 → 傾向單套（男）
Y_PRESENT_MIN = 0.15  # Y 原始比 ≥ 此 → 有 Y
ANEUPLOIDY_LO = 0.75  # NDC < 此 或
ANEUPLOIDY_HI = 1.25  #     > 此 → 疑似非整倍體


def parse_mosdepth_summary(path):
    """回傳 ({chrom: mean}, {chrom: length})。每 contig 優先取 *_region，否則整條 mean。"""
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
                region[name[:-len("_region")]] = mean
            else:
                whole[name] = mean
                try:
                    length[name] = int(p[li])
                except (ValueError, IndexError):
                    length[name] = 0
    means = {}
    for c in set(list(whole) + list(region)):
        # 優先取 on-target 的 *_region mean；但 mosdepth 會對「不在 --by BED 裡的 contig」
        # 也吐一行 *_region=0（例如 WGS 用 autosome BED 時的 chrX/chrY/chrM），那代表「此
        # contig 沒有區間」而非「深度 0」。若照抄 0 會把性染色體判成單套缺失（男生→X0?）。
        # 因此 *_region 僅在 >0 時採用，否則回退整條 contig 的 whole mean。
        r = region.get(c)
        means[c] = r if (r is not None and r > 0) else whole.get(c, r)
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
    """從**原始**覆蓋比推核型。回傳 (karyotype, x_ratio, y_ratio)。"""
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


def _xy_counts(kary):
    """從核型字串數 X / Y 拷貝數；unknown/ambiguous → (None, None)（不正規化性染色體）。"""
    k = kary.replace("?", "").upper()
    if k in ("UNKNOWN", "AMBIGUOUS", ""):
        return None, None
    return k.count("X"), k.count("Y")


def expected_ratio(chrom, kary):
    """該 contig 在估計核型下的期望原始比（相對體染色體 diploid=1.0）。
    autosome=1.0；chrX=nX*0.5；chrY=nY*0.5。回 None 代表不正規化（未知核型）。"""
    if chrom in AUTOSOMES:
        return 1.0
    nx, ny = _xy_counts(kary)
    if nx is None:
        return 1.0                      # 未知核型 → 不對性染色體正規化
    if chrom == "chrX":
        return nx * 0.5
    if chrom == "chrY":
        return ny * 0.5
    return 1.0


def analyze(means, length, declared_sex, seq_type="WGS"):
    """回傳 dict 結果（含 per-contig ratio/ndc 與 warnings）。純函式、好測試。"""
    auto_means = [means[c] for c in AUTOSOMES if means.get(c) is not None]
    baseline = statistics.median(auto_means) if auto_means else 0.0
    est_kary, x_ratio, y_ratio = infer_karyotype(means, baseline)
    declared_kary = normalize_declared_sex(declared_sex)

    ratio, ndc = {}, {}
    for c, m in means.items():
        r = (m / baseline) if baseline else 0.0
        ratio[c] = r
        if c in MITO:
            ndc[c] = None                # chrM 多拷貝，NDC 無意義
            continue
        exp = expected_ratio(c, est_kary)
        ndc[c] = (r / exp) if (exp and exp > 0) else None   # 期望 0（如女 chrY）→ NA

    # aneuploidy：NDC 偏離 1.0（非 mito、有 NDC 值）。性染色體正常會 ≈1.0（已正規化）。
    flags = [(c, ndc[c]) for c in (AUTOSOMES + ["chrX", "chrY"])
             if ndc.get(c) is not None and (ndc[c] < ANEUPLOIDY_LO or ndc[c] > ANEUPLOIDY_HI)]

    warnings = []
    if declared_kary != "unknown" and est_kary not in ("unknown", "ambiguous") \
            and est_kary != declared_kary:
        warnings.append("SEX MISMATCH：samplesheet 宣告 %s，資料推得 %s"
                        "（可能 sample swap 或性染色體 aneuploidy，請人工確認）"
                        % (declared_kary, est_kary))
    if est_kary in ("XXY?", "X0?", "ambiguous"):
        warnings.append("性染色體核型異常/不明：推得 %s（X_ratio=%s, Y_ratio=%s）"
                        % (est_kary, _fmt(x_ratio), _fmt(y_ratio)))
    for c, v in flags:
        warnings.append("%s NDC=%.2f → 疑似非整倍體，本染色體 SNV 基因型可能不準，"
                        "考慮手動 -ploidy N 重跑" % (c, v))
    if seq_type == "WES":
        warnings.append("注意：WES 以 target 覆蓋估計，性別/ploidy 判斷較粗；"
                        "aneuploidy 建議以 gCNV 交叉驗證")

    return {
        "baseline": baseline,
        "estimated_karyotype": est_kary,
        "declared_karyotype": declared_kary,
        "x_ratio": x_ratio,
        "y_ratio": y_ratio,
        "ratio": ratio,
        "ndc": ndc,
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
def _ordered_chroms(keys):
    # 只輸出主要 contig（chr1-22, X, Y, M）。alt/decoy/unplaced 對 sex/ploidy 判斷無意義、
    # 只是雜訊，一律不列（也與 DRAGEN ploidy.vcf 只含主要 contig 對齊）。analyze() 的
    # baseline / 核型 / aneuploidy 本就只用主要 contig，這裡純粹是「輸出過濾」，不影響判定。
    order = AUTOSOMES + ["chrX", "chrY", "chrM"]
    return [c for c in order if c in keys]


def write_vcf(path, sample, res, seq_type):
    means, length = res["means"], res["length"]
    flagged = {c for c, _ in res["aneuploidy"]}
    with open(path, "w") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write("##source=NCKUH_PLOIDY_MOSDEPTH\n")
        w.write("##seqType=%s\n" % seq_type)
        w.write("##estimatedSexKaryotype=%s\n" % res["estimated_karyotype"])
        w.write("##referenceSexKaryotype=%s\n" % res["declared_karyotype"])
        w.write('##FILTER=<ID=SUSPECT,Description="Normalized coverage suggests possible '
                'aneuploidy for this contig">\n')
        w.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Contig end">\n')
        w.write('##FORMAT=<ID=DC,Number=1,Type=Float,Description="Mean depth of coverage '
                '(mosdepth)">\n')
        w.write('##FORMAT=<ID=NDC,Number=1,Type=Float,Description="Normalized depth of coverage '
                'relative to expected ploidy for the estimated karyotype (~1.0 = as-expected). '
                'NA when expected copy number is 0">\n')
        w.write('##FORMAT=<ID=RATIO,Number=1,Type=Float,Description="Raw depth ratio relative to '
                'autosomal median (~1.0 diploid, ~0.5 haploid); sex/dosage evidence">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        for c in _ordered_chroms(means.keys()):
            dc = means[c]
            r = res["ratio"].get(c)
            v = res["ndc"].get(c)
            filt = "SUSPECT" if c in flagged else "PASS"
            end = length.get(c, 0)
            ndc_s = ("%.4f" % v) if v is not None else "."
            ratio_s = ("%.4f" % r) if r is not None else "."
            w.write("%s\t1\t.\tN\t.\t.\t%s\tEND=%d\tDC:NDC:RATIO\t%.4f:%s:%s\n"
                    % (c, filt, end, dc, ndc_s, ratio_s))


def write_qc(path, sample, res, seq_type):
    with open(path, "w") as w:
        w.write("# Ploidy QC — %s (%s)\n" % (sample, seq_type))
        w.write("declared_sex_karyotype : %s\n" % res["declared_karyotype"])
        w.write("estimated_sex_karyotype: %s\n" % res["estimated_karyotype"])
        sex_ok = (res["declared_karyotype"] == "unknown"
                  or res["estimated_karyotype"] in ("unknown", "ambiguous")
                  or res["declared_karyotype"] == res["estimated_karyotype"])
        w.write("sex_check              : %s\n" % ("OK" if sex_ok else "MISMATCH"))
        w.write("source                 : NCKUH mosdepth (NDC normalized to expected karyotype; "
                "~1.0 = as-expected. RATIO = raw ratio vs autosomes)\n")
        w.write("autosomal_baseline_mean: %.4f\n" % res["baseline"])
        w.write("chrX_ratio             : %s\n" % _fmt(res["x_ratio"]))
        w.write("chrY_ratio             : %s\n" % _fmt(res["y_ratio"]))
        w.write("\n--- per-chromosome NDC / RATIO ---\n")
        for c in _ordered_chroms(res["means"].keys()):
            v = res["ndc"].get(c)
            r = res["ratio"].get(c)
            w.write("%-6s NDC=%-6s RATIO=%s\n"
                    % (c, ("%.3f" % v) if v is not None else "NA",
                       ("%.3f" % r) if r is not None else "NA"))
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

    for msg in res["warnings"]:
        sys.stderr.write("[ploidy_check] WARN: %s\n" % msg)
    sys.stderr.write("[ploidy_check] %s declared=%s estimated=%s baseline=%.2f\n"
                     % (a.sample, res["declared_karyotype"],
                        res["estimated_karyotype"], res["baseline"]))


if __name__ == "__main__":
    main()
