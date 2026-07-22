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
  * haploid（男性 non-PAR chrX/chrY；整叢皆單套）→ 只重建單一單體 → GT=1 的 hemizygous
    MNV（不帶 PS）；chrM 不走此路（多拷貝異質性，見下）
  * 孤立變異（叢集只含 1 顆）→ 原行輸出，完全不動

重疊套用規則（處理如 SUZ12 的 GAAA>G 與 A>ATT 在同一單體重疊）：沿參考游標套用，
遇到 POS < 游標（重疊已消耗的 ref）時，只補上該變異 ALT 中「尚未輸出」的尾段
（= 插入的部分）。SUZ12：hapB = "G"(delAAA) 之後補 "TT" → GAAA>GTT，符合
c.2168_2170delAAAinsTT。

輸出與相依
----------
* 輸入：單樣本 VCF(.gz)；輸出：未壓縮 VCF（交由 Nextflow bgzip+tabix）。
* 合成出的 biallelic 紀錄會「繼承一顆代表變異（anchor＝叢集內參考足跡最寬、ploidy 與合成
  結果一致的 biallelic 顆）」的整組 FORMAT：AD/DP/GQ/VAF/PL… 原封保留，當作該 compound 的
  讀取支持與等位分數（符合 bug report §4：保留原始 locus 的 VAF/AD，不用 2 元 AD 重算成誤導
  的 1.0），只覆寫 GT（diploid → 重建的 phased GT；haploid → GT=1，不帶 PS）與 PS，並在
  INFO 標 COMBINED=<n>。QUAL/FILTER 沿用 anchor。
* 下列情形「不重建、原封通過」（保留各來源紀錄原本的 AD/DP/VAF，交由下游 bcftools norm
  拆分），以免捏造深度：
    (a) 叢集重建後有 2 個 ALT（1|2，如原生 multiallelic 1/2）；
    (b) 叢集內找不到可當 anchor 的 biallelic 紀錄；
    (c) 叢集內「混 ploidy」（haploid 與 diploid 同叢）；
    (d) chrM 的 haploid 叢集（多拷貝異質性，不宜當單一分子合）。
* 原始（未合、孤立）紀錄一律原封輸出（保留所有 FORMAT）。
* header 會補上 ##INFO=<ID=COMBINED> 與 ##FORMAT=<ID=PS>（若原本沒有）。

⚠️ 早期版本合成紀錄只輸出 GT:PS，會把 AD/DP/VAF 丟成 '.'（三級 DRAGEN AD 消失 bug、
   145k+ 筆受影響）。現改為「anchor 繼承 + 多等位/混ploidy/chrM 退回原封通過」，都保住深度。
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


def trim_alleles(pos: int, ref: str, alts: List[str]):
    """
    Ref-free 最小化（右修剪→左修剪），只處理 biallelic。移除 caller 表示法裡的
    「padding」參考鹼基，避免非真正重疊被誤判成重疊（如 DV 的 GAAA>GAA →修剪為
    GA>G，就不會假性蓋到 31998953 的 SNV）。symbolic/spanning ALT 不動。
    """
    if len(alts) != 1:
        return pos, ref, alts
    alt = alts[0]
    if not ref or not alt or alt.startswith("<") or alt == "*" or "]" in alt or "[" in alt:
        return pos, ref, alts
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, [alt]


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
    # haploid（單一拷貝，如男性 non-PAR chrX/chrY）：同一條上的變異一律 cis。
    if len(v.alleles) == 1 and all(len(m.alleles) == 1 for m in cluster):
        return True
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


def reconstruct_haploid(cluster: List[Var], fetch: Callable[[str, int, int], str]):
    """單套（haploid）重建：只有一條單體，把叢集內所有帶 alt 的變異沿參考游標套上，
    重建那唯一一條序列。回傳 (pos, ref, [alt], [1])（hemizygous）或 None（無 alt/無法合）。
    永遠只會有 1 個 ALT（單套 → 不可能 1|2）。"""
    chrom = cluster[0].chrom
    span_start = min(v.pos for v in cluster)
    span_end = max(v.end for v in cluster)
    ref_seq = fetch(chrom, span_start, span_end)
    if not ref_seq:
        return None
    edits = []
    for v in cluster:
        a = v.alleles[0] if v.alleles else 0
        if 0 < a <= len(v.alts):
            edits.append((v.pos, v.ref, v.alts[a - 1]))
    if not edits:
        return None
    hap = build_hap(span_start, ref_seq, edits)
    if hap == ref_seq:
        return None
    return span_start, ref_seq, [hap], [1]


# ─────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────
_MITO_CONTIGS = {"chrM", "chrMT", "MT", "M"}
def _open(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


# ─────────────────────────────────────────────────────────────
# 合成紀錄的 FORMAT 繼承（保住 AD/DP/VAF…；只覆寫 GT、PS）
# ─────────────────────────────────────────────────────────────
def _split_sample(line: str, sample_col: int):
    """把一行 VCF 拆成 (fields, FORMAT keys, {key: value})（指定 sample 欄）。"""
    f = line.rstrip("\n").split("\t")
    keys = f[8].split(":") if len(f) > 8 else []
    vals = f[9 + sample_col].split(":") if len(f) > 9 + sample_col else []
    return f, keys, dict(zip(keys, vals))


def _fmt_anchor(cluster: List[Var], sample_col: int) -> Optional[Var]:
    """挑 FORMAT 捐贈者：叢集內「biallelic（len(alts)==1）且有 FORMAT」中足跡最寬
    （tie → 最左）的一顆。找不到回 None（呼叫端會退回原封通過）。挑最寬是因為
    compound 的主事件（如 SUZ12 的 GAAA>G 缺失）通常足跡最寬，其 AD/DP 最能代表
    整個 compound 的讀取支持。只在「同 ploidy」的叢集上呼叫，故 anchor 的 ploidy 必與
    合成結果一致 → AD(Number=R)/PL(Number=G) 長度天生對得上（diploid 或 haploid 皆然）。"""
    cands = [v for v in cluster
             if len(v.alts) == 1 and len(_split_sample(v.line, sample_col)[1]) > 0]
    if not cands:
        return None
    cands.sort(key=lambda v: (-len(v.ref), v.pos))
    return cands[0]


def _render_merged(chrom: str, pos: int, ref: str, alt: str, gtstr: str,
                   ps: Optional[str], n_combined: int, anchor: Var, sample_col: int) -> str:
    """組出合成後的 biallelic 紀錄：沿用 anchor 的 QUAL/FILTER/FORMAT，只覆寫 GT、PS，
    INFO 設 COMBINED=<n>。因 anchor 與合成結果同 ploidy 且 biallelic，其 AD(Number=R)/
    PL(Number=G) 元素數與合成後的 biallelic 一致，直接繼承即為正確長度。ps 為 None
    （haploid hemizygous）時不加、也不覆寫 PS。"""
    f, keys, d = _split_sample(anchor.line, sample_col)
    qual = f[5] if len(f) > 5 else "."
    filt = f[6] if len(f) > 6 else "."
    if "GT" not in keys:
        keys = ["GT"] + keys
    if ps is not None and "PS" not in keys:
        keys = keys + ["PS"]
    vals = []
    for k in keys:
        if k == "GT":
            vals.append(gtstr)
        elif k == "PS":
            vals.append(ps if ps is not None else d.get("PS", "."))
        else:
            vals.append(d.get(k, "."))
    return "\t".join([chrom, str(pos), ".", ref, alt, qual, filt,
                      "COMBINED=%d" % n_combined, ":".join(keys), ":".join(vals)])


def process(in_vcf: str, out_vcf: str, fetch: Callable, max_gap: int,
            sample_col: int = 0) -> dict:
    """主流程；回傳統計。fetch 可注入（測試用）。"""
    stats = {"clusters_merged": 0, "clusters_haploid": 0, "clusters_fallback": 0,
             "records_in": 0, "records_out": 0}
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
            is_mito = chrom in _MITO_CONTIGS
            recs = []
            for cl in clusters:
                if len(cl) == 1:                       # 孤立顆 → 原行輸出，完全不動
                    recs.append((cl[0].pos, cl[0].line))
                    continue
                ploidies = {len(v.alleles) for v in cl}
                if ploidies == {2}:                    # 全 diploid → 雙單體重建
                    res = reconstruct(cl, fetch)
                    multi = res is not None and len(res[2]) > 1
                elif ploidies == {1} and not is_mito:  # 全 haploid（非 chrM）→ 單套重建
                    res, multi = reconstruct_haploid(cl, fetch), False
                else:                                  # 混 ploidy、chrM haploid、其他 → 不合
                    res, multi = None, False
                anchor = _fmt_anchor(cl, sample_col)
                # 無法重建 / 重建成多 ALT（1|2）/ 無 biallelic anchor → 原封通過（保留 AD）。
                if res is None or multi or anchor is None:
                    for v in cl:
                        recs.append((v.pos, v.line))
                    stats["clusters_fallback"] += 1
                    continue
                pos, ref, alt_list, gt = res
                if len(gt) > 1:                        # diploid：phased GT + PS
                    gtstr = "|".join(str(g) for g in gt)
                    ps = next((v.ps for v in cl if v.ps), str(pos))
                else:                                  # haploid hemizygous：GT=1，不帶 PS
                    gtstr, ps = str(gt[0]), None
                    stats["clusters_haploid"] += 1
                recs.append((pos, _render_merged(
                    chrom, pos, ref, alt_list[0], gtstr, ps, len(cl),
                    anchor, sample_col)))
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
            # 最小化後再參與叢集/重建（移除 padding 假性重疊）；passthrough 仍用原始 line。
            tpos, tref, talts = trim_alleles(pos, ref, alt.split(","))
            v = Var(chrom, tpos, tref, talts, line,
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
    sys.stderr.write(
        "[combine_phased] in=%d out=%d merged_clusters=%d (haploid=%d) "
        "passthrough_clusters=%d\n"
        % (st["records_in"], st["records_out"], st["clusters_merged"],
           st["clusters_haploid"], st["clusters_fallback"]))


if __name__ == "__main__":
    main()
