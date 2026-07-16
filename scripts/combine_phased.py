#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=========================================================
WGS/WES Germline Analysis Pipeline - combine_phased.py
=========================================================
Author   : Po-Yu Lin (林伯昱)
Institute: Department of Neurology and
           Department of Genomic Medicine,
           National Cheng Kung University Hospital
Contact  : p88124019@gs.ncku.edu.tw

Copyright (c) 2026, Po-Yu Lin (林伯昱)
Licensed under the GNU General Public License v3.0

DISCLAIMER: Provided "as is" without warranty. Users are solely responsible
for validating and interpreting all results.
=========================================================
scripts/combine_phased.py
=========================

把「同一單體上、重疊或鄰近」的變異，用**局部單體重建**合成單一 canonical MNV。
供二級 NCKUH（各 caller phase 後、進 BCFTOOLS_ENSEMBLE 前）與三級 DRAGEN 共用。
單一樣本 VCF；相依零套件（只用 Python 標準庫 + 自帶 faidx 隨機存取）。

語義（對齊 DRAGEN --vc-combine-phased-variants-distance，並依臨床討論擴充）
----------------------------------------------------------------------
以「參考足跡」而非 POS 距離判斷是否相連：每顆變異影響 ref 區間 [POS, POS+len(REF)-1]。
兩顆變異會被歸到同一叢集（cluster）若：
  (1) 足跡**重疊**（物理上必須共列）；或
  (2) 足跡空隙 gap ≤ --max-gap 且**在同一單體上（cis）**。
其中「cis」定義（保守）：任一顆為 homozygous（在兩條單體上都在），或兩顆皆
phased、同一個 PS、且非參考等位都落在同一條單體。trans（同 PS 但落在不同單體）
且**不重疊** → 不合（本來就是兩顆）。無法判定 phase 的（未 phase 的 het）→ 只在
**重疊**時才合，否則原封通過。

叢集重建：對 hapA / hapB 兩條單體，各自沿參考游標套用該單體上的變異，重建序列；
再依 {ref, hapA_seq, hapB_seq} 去重得到 REF / ALT，並重新給 GT（phased）。因此：
  * het 相鄰 cis（含 SUZ12 del+ins 重疊）→ 一條單體改、另一條為 ref → 0|1 的 MNV
  * hom + hom            → 兩條單體都改成同序列 → 1|1 的 MNV
  * hom + het 重疊       → 兩條單體不同 → 1|2（揭露「其實不是 hom」）
  * 孤立變異（叢集只含 1 顆）→ 原行輸出，完全不動

重疊套用規則（處理如 SUZ12 的 GAAA>G 與 A>ATT 在同一單體重疊）：沿參考游標套用，
遇到 POS < 游標（重疊已消耗的 ref）時，只補上該變異 ALT 中「尚未輸出」的尾段
（= 插入的部分）。SUZ12：hapB = "G"(delAAA) 之後補 "TT" → GAAA>GTT，符合
c.2168_2170delAAAinsTT。

輸出與相依
----------
* 輸入：單樣本 VCF(.gz)；輸出：未壓縮 VCF（交由 Nextflow bgzip+tabix）。
* 合成出的紀錄 FORMAT 只保留 GT:PS（PS = 叢集第一顆的 PS，或 span 起點），並在
  INFO 加 COMBINED=<n>；原始（未合）紀錄原封輸出。合成紀錄的深度/品質欄位（AD/DP…）
  故意不帶（對 MNV 語義不明），下游若需要再議。
* header 會補上 ##INFO=<ID=COMBINED> 與 ##FORMAT=<ID=PS>（若原本沒有）。
"""

import argparse
import gzip
import sys
from dataclasses import dataclass, field
from typing import Callable, List, Optional


# ─────────────────────────────────────────────────────────────
# FASTA 隨機存取（自帶 faidx，免 pysam/samtools）
# ─────────────────────────────────────────────────────────────
class Faidx:
    def __init__(self, fasta: str):
        self.idx = {}
        with open(fasta + ".fai") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                # name, length, offset, linebases, linewidth
                self.idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
        self.fh = open(fasta, "rb")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive [start, end]。"""
        length, offset, linebases, linewidth = self.idx[chrom]
        if start < 1:
            start = 1
        if end > length:
            end = length
        if end < start:
            return ""

        def byte_of(pos: int) -> int:
            z = pos - 1
            return offset + (z // linebases) * linewidth + (z % linebases)

        b0 = byte_of(start)
        b1 = byte_of(end) + 1
        self.fh.seek(b0)
        raw = self.fh.read(b1 - b0)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()


# ─────────────────────────────────────────────────────────────
# 變異資料結構 + GT 解析
# ─────────────────────────────────────────────────────────────
@dataclass
class Var:
    chrom: str
    pos: int
    ref: str
    alts: List[str]
    line: str                       # 原始 VCF 行（passthrough 用）
    alleles: List[int] = field(default_factory=list)   # GT 索引，如 [0,1]；0=ref、-1=缺失
    phased: bool = False
    ps: Optional[str] = None

    @property
    def end(self) -> int:
        return self.pos + len(self.ref) - 1

    def allele_on(self, hap: int) -> int:
        """該單體(hap=0/1)的等位索引；不足或缺失回 0（視為 ref）。"""
        if hap < len(self.alleles) and self.alleles[hap] > 0:
            return self.alleles[hap]
        return 0

    def is_hom_alt(self) -> bool:
        return (len(self.alleles) == 2 and self.alleles[0] > 0
                and self.alleles[0] == self.alleles[1])

    def het_hap(self) -> Optional[int]:
        """若為單一單體帶非參考（het），回該 hap（0/1）；否則 None。"""
        haps = [h for h in (0, 1) if self.allele_on(h) > 0]
        return haps[0] if len(haps) == 1 else None


def parse_gt(fmt_keys: List[str], sample_vals: List[str]):
    d = dict(zip(fmt_keys, sample_vals))
    gt = d.get("GT", ".")
    phased = "|" in gt
    parts = gt.replace("|", "/").split("/")
    alleles = [int(p) if p.isdigit() else -1 for p in parts]
    ps = d.get("PS")
    if ps in (None, "", "."):
        ps = None
    return alleles, phased, ps


# ─────────────────────────────────────────────────────────────
# 叢集判定
# ─────────────────────────────────────────────────────────────
def _linkable_cis(v: Var, cluster: List[Var]) -> bool:
    """v 與 cluster 是否同一單體(cis)可合（保守）。"""
    if v.is_hom_alt() or any(m.is_hom_alt() for m in cluster):
        return True
    # 兩邊都必須是 phased、同一 PS、同一 hap
    if not v.phased or v.ps is None:
        return False
    vhap = v.het_hap()
    if vhap is None:
        return False
    chaps, cps = set(), set()
    for m in cluster:
        if not m.phased or m.ps is None:
            return False
        cps.add(m.ps)
        mh = m.het_hap()
        if mh is None:
            return False
        chaps.add(mh)
    return cps == {v.ps} and chaps == {vhap}


def cluster_vars(variants: List[Var], max_gap: int) -> List[List[Var]]:
    """依 POS 排序後分叢；重疊必合、cis 且 gap≤max_gap 合，其餘斷開。"""
    clusters: List[List[Var]] = []
    cur: List[Var] = []
    cur_end = 0
    for v in sorted(variants, key=lambda x: (x.pos, x.end)):
        if not cur:
            cur, cur_end = [v], v.end
            continue
        overlap = v.pos <= cur_end
        near = v.pos <= cur_end + max_gap + 1
        if overlap or (near and _linkable_cis(v, cur)):
            cur.append(v)
            cur_end = max(cur_end, v.end)
        else:
            clusters.append(cur)
            cur, cur_end = [v], v.end
    if cur:
        clusters.append(cur)
    return clusters


# ─────────────────────────────────────────────────────────────
# 單體重建
# ─────────────────────────────────────────────────────────────
def build_hap(span_start: int, ref_seq: str, edits: List[tuple]) -> str:
    """
    edits: list of (pos, ref, alt)（此單體上的變異，已排序、可能重疊）。
    沿參考游標重建；遇重疊只補 ALT 尚未輸出的尾段。
    """
    out = []
    cursor = span_start                       # 下一個要輸出的 ref 位置（1-based）
    span_end = span_start + len(ref_seq) - 1
    for pos, r, a in sorted(edits):
        if pos >= cursor:
            out.append(ref_seq[cursor - span_start: pos - span_start])   # 中間未變 ref
            out.append(a)
            cursor = pos + len(r)
        else:
            consumed = cursor - pos           # 該變異 ref 已被前一顆消耗的長度
            if consumed < len(a):
                out.append(a[consumed:])       # 只補尾段（多半是插入的部分）
            cursor = max(cursor, pos + len(r))
    if cursor <= span_end:
        out.append(ref_seq[cursor - span_start:])
    return "".join(out)


def reconstruct(cluster: List[Var], fetch: Callable[[str, int, int], str]):
    """回傳 (pos, ref, alt_list, gt_list) 或 None（無非參考、無法合）。"""
    chrom = cluster[0].chrom
    span_start = min(v.pos for v in cluster)
    span_end = max(v.end for v in cluster)
    ref_seq = fetch(chrom, span_start, span_end)
    if not ref_seq:
        return None

    hap_seqs = []
    for h in (0, 1):
        edits = [(v.pos, v.ref, v.alts[v.allele_on(h) - 1])
                 for v in cluster if v.allele_on(h) > 0]
        hap_seqs.append(build_hap(span_start, ref_seq, edits))

    if hap_seqs[0] == ref_seq and hap_seqs[1] == ref_seq:
        return None                            # 兩條都 ref，不該合

    alleles = [ref_seq]
    for s in hap_seqs:
        if s != ref_seq and s not in alleles:
            alleles.append(s)
    gt = [alleles.index(hap_seqs[0]), alleles.index(hap_seqs[1])]
    return span_start, ref_seq, alleles[1:], gt


# ─────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────
def _open(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def process(in_vcf: str, out_vcf: str, fetch: Callable, max_gap: int,
            sample_col: int = 0) -> dict:
    """主流程；回傳統計。fetch 可注入（測試用）。"""
    stats = {"clusters_merged": 0, "records_in": 0, "records_out": 0}
    header, chrom_vars, order_chrom = [], {}, []
    fmt_extra = ['##INFO=<ID=COMBINED,Number=1,Type=Integer,'
                 'Description="Number of source records combined into this MNV '
                 'by combine_phased.py">',
                 '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">']

    def flush(w):
        """把累積的每染色體變異分叢、重建、輸出（維持座標順序）。"""
        for chrom in order_chrom:
            vs = chrom_vars[chrom]
            clusters = cluster_vars(vs, max_gap)
            recs = []
            for cl in clusters:
                if len(cl) == 1:
                    recs.append((cl[0].pos, cl[0].line))
                    continue
                res = reconstruct(cl, fetch)
                if res is None:
                    for v in cl:               # 無法合 → 原行輸出
                        recs.append((v.pos, v.line))
                    continue
                pos, ref, alt_list, gt = res
                gtstr = "|".join(str(g) for g in gt)
                ps = next((v.ps for v in cl if v.ps), str(pos))
                info = "COMBINED=%d" % len(cl)
                recs.append((pos, "\t".join([
                    chrom, str(pos), ".", ref, ",".join(alt_list),
                    ".", "PASS", info, "GT:PS", "%s:%s" % (gtstr, ps)])))
                stats["clusters_merged"] += 1
            for _, line in sorted(recs, key=lambda x: x[0]):
                w.write(line if line.endswith("\n") else line + "\n")
                stats["records_out"] += 1

    with _open(in_vcf) as fin, open(out_vcf, "wt") as w:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    for h in fmt_extra:        # 在 #CHROM 前補 header 定義
                        w.write(h + "\n")
                    w.write(line)
                else:
                    w.write(line)
                continue
            stats["records_in"] += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            fmt_keys = f[8].split(":") if len(f) > 8 else []
            sample_vals = f[9 + sample_col].split(":") if len(f) > 9 + sample_col else []
            alleles, phased, ps = parse_gt(fmt_keys, sample_vals)
            v = Var(chrom, pos, ref, alt.split(","), line,
                    alleles=alleles, phased=phased, ps=ps)
            if chrom not in chrom_vars:
                chrom_vars[chrom] = []
                order_chrom.append(chrom)
            chrom_vars[chrom].append(v)
        flush(w)
    return stats


def main():
    ap = argparse.ArgumentParser(description="Combine phased/overlapping variants into MNVs.")
    ap.add_argument("--in", dest="inp", required=True, help="input single-sample VCF(.gz)")
    ap.add_argument("--out", required=True, help="output VCF (uncompressed)")
    ap.add_argument("--fasta", required=True, help="reference FASTA (needs .fai)")
    ap.add_argument("--max-gap", type=int, default=1,
                    help="max untouched-ref gap (bp) between cis variants to combine "
                         "(overlaps always combine). Default 1.")
    ap.add_argument("--sample-index", type=int, default=0, help="0-based sample column")
    a = ap.parse_args()
    fa = Faidx(a.fasta)
    st = process(a.inp, a.out, fa.fetch, a.max_gap, a.sample_index)
    sys.stderr.write("[combine_phased] in=%d out=%d merged_clusters=%d\n"
                     % (st["records_in"], st["records_out"], st["clusters_merged"]))


if __name__ == "__main__":
    main()
