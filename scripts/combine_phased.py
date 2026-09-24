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
且**不重疊** → 不合（本來就是兩顆）。重疊的一定歸同一叢，但**只有 phase 已知才合成**：
叢集內 ≥2 顆 het 且不全是「同一個 phase set 的 phased」（未 phase、或 PS 不同）→ 原封通過
（2026-09，見下方 ⚠️）。phased 但沒有 PS 欄的，依 VCF 規格屬於同一個隱含的 phase set。

叢集重建：對 hapA / hapB 兩條單體，各自沿參考游標套用該單體上的變異，重建序列；
再依 {ref, hapA_seq, hapB_seq} 去重得到 REF / ALT，並重新給 GT（phased）。因此：
  * het 相鄰 cis（含 SUZ12 del+ins 重疊）→ 一條單體改、另一條為 ref → 0|1 的 MNV
  * hom + hom            → 兩條單體都改成同序列 → 1|1 的 MNV
  * hom + het 重疊       → 兩條單體不同 → 1|2（揭露「其實不是 hom」）
  * haploid（男性 non-PAR chrX/chrY；整叢皆單套）→ 只重建單一單體 → GT=1 的 hemizygous
    MNV（不帶 PS）；chrM 不走此路（多拷貝異質性，見下）
  * 孤立變異（叢集只含 1 顆）→ 原行輸出，完全不動
  * 該 sample 的 GT 沒有 ALT（./.、0/0、單套 0 或 .；典型是 DeepVariant 否決的候選
    FILTER=RefCall）→ **不參與叢集**，原行輸出（二級 BCFTOOLS_ENSEMBLE 之後會把 DV 的丟掉）

重疊套用規則（處理如 SUZ12 的 GAAA>G 與 A>ATT 在同一單體重疊）：沿參考游標套用，
遇到 POS < 游標（這顆的前 k 個 ref 鹼基已被前面的變異處理過）時，只有當這 k 個鹼基在這顆裡
**沒有改變**（ALT 前 k 個 == REF 前 k 個，即 VCF 的前導／錨定鹼基）才把其餘部分接上去。
SUZ12：hapB = "G"(delAAA) 之後，A>ATT 的 k=1、錨定鹼基 A 不變 → 補上 "TT" → GAAA>GTT，符合
c.2168_2170delAAAinsTT（輸出時再最小化成 31998951 AAA>TT，見下）。其他重疊（缺失範圍內的 SNV、
包在大缺失裡的小缺失、錨在缺失中間的插入…）代表兩顆在同一條單體上互相矛盾 → 整叢原封通過，
不吃掉、不截斷任何一個 allele（2026-09）。同一位置時先套改變錨定鹼基的 SNV/MNV、再套 indel；
完全相同的重複紀錄只算一次。

輸出與相依
----------
* 輸入：單樣本 VCF(.gz)；輸出：未壓縮 VCF（交由 Nextflow bgzip+tabix）。
* 合成出的 REF/ALT 在輸出前再**最小化**（trim_alleles：去掉 REF/ALT 共同的前後鹼基，至少各留 1 個，
  所以純 indel 的前導鹼基會保留）。叢集範圍從第一個成分的 POS 起算，成分若是 indel 就帶著它的
  前導鹼基；不修剪的話，DV、HC 對同一個 compound 拆成分的方式不同 → 前導鹼基長度不同 → POS 不同，
  ensemble merge 合不起來（2026-09，VAL55：820 個變異在報告裡拆成 DV 一列 + HC 一列）。
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
    (d) chrM 的 haploid 叢集（多拷貝異質性，不宜當單一分子合）；
    (e) phase 未知：≥2 顆 het 不在同一個 phase set（stderr 的 phase_unknown=）；
    (f) 重疊互相矛盾：某顆的 allele 會被吃掉或截斷（stderr 的 overlap_conflict=）。
* 原始（未合、孤立）紀錄一律原封輸出（保留所有 FORMAT）。
* header 會補上 ##INFO=<ID=COMBINED> 與 ##FORMAT=<ID=PS>（若原本沒有）。

⚠️ 早期版本合成紀錄只輸出 GT:PS，會把 AD/DP/VAF 丟成 '.'（三級 DRAGEN AD 消失 bug、
   145k+ 筆受影響）。現改為「anchor 繼承 + 多等位/混ploidy/chrM 退回原封通過」，都保住深度。

⚠️ 沒有 ALT 的紀錄曾經會參與叢集（2026-09 修正）。DeepVariant 的 VCF 保留它否決的候選
   （FILTER=RefCall）；舊版讓它們照「足跡重疊必合」進叢集，又因 anchor 只挑「足跡最寬」、
   不看有沒有被 call，於是被否決的較寬候選（常是缺失）蓋住真的 call（如 SNV）時：
     - 合成紀錄沿用被否決候選的 QUAL / FILTER=RefCall / GQ / DP / AD / VAF / PL；
     - 從被否決候選的 POS 起、補上參考鹼基重寫 → 與 HC 同一變異的 POS 對不上，merge 合不起來，
       三級 norm 後又變成同一個變異 → 報告裡拆成 CALLERS=DV（深度是被否決候選的）+ CALLERS=HC 兩列；
     - 未 phase 的 het 也被寫成 0|1 並給假的 PS；被否決候選還會把兩顆不相干的 call 串成同一叢。
   合成紀錄的 GT 有 ALT，所以二級 BCFTOOLS_ENSEMBLE 的 DV `GT="alt"` 過濾擋不掉。實例 VAL55：ensemble 有
   28,050 筆 FILTER=RefCall 全部帶 COMBINED；三級 23,023 個變異被拆成 DV 一列 + HC 一列。

⚠️ 重疊但 phase 未知的 het 曾被當成在同一條單體上合成（2026-09 修正，上方 (e)(f)）。未 phase 的
   GT 在 reconstruct() 裡依位置都落在同一條單體；於是缺失範圍內的 SNV、包在大缺失裡的小缺失被吃掉
   （PASS 變異從報告消失），兩個重疊的缺失被併成更長的缺失、錨在缺失中間的插入被截斷（寫出兩個 caller
   都沒 call 的 allele），還給了假的 0|1 + PS。實例 VAL-10（DRAGEN 女性 WGS，只拿 PASS）：104,277 個合成中
   9,299 個是 phase 未知的 het；7,927 個 PASS allele 在合成時消失（unphased 6,778 叢、phased 88 叢），
   1,435 叢寫出沒人 call 的 allele。phased 也會發生：同一條單體上互相矛盾的 call（如 hom 缺失裡的 het 缺失）、
   同一位置的 SNV + 插入（舊版依字串排序先套插入，SNV 被吃掉）。
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
class RebuildConflict(Exception):
    """這一叢無法在不吃掉、不截斷任何 allele 的情況下重建（重疊互相矛盾，或 ALT 是 * / symbolic）。
    呼叫端（plan_cluster）接到後整叢原封通過。"""


def _edit_order(e):
    """build_hap 的套用順序：依 POS；同一位置先套「改變錨定鹼基」的 SNV/MNV，再套保留錨定鹼基的
    indel。舊版用字串排序，同一鹼基的 SNV + 插入（A>G、A>AT）會先套插入，SNV 被吃掉。"""
    pos, r, a = e
    return (pos, a[:1] == r[:1], r, a)


def build_hap(span_start: int, ref_seq: str, edits: List[tuple]) -> str:
    """
    edits: list of (pos, ref, alt)（此單體上的變異，可能重疊）。沿參考游標重建。
    重疊時，這顆的前 k 個 ref 鹼基已被前面的變異處理過（k = 游標 - POS）。只有這 k 個鹼基在這顆裡
    沒有改變（k ≤ len(ref)、k ≤ len(alt)、alt[:k] == ref[:k]，即 VCF 的前導／錨定鹼基）才接上其餘部分：
      SUZ12 的 A>ATT 接在 GAAA>G 之後：k=1、錨定 A 不變 → 補上 TT；
      同一 POS 的 SNV + 缺失：缺失的錨定鹼基已被 SNV 改寫，其餘照刪。
    其他重疊（缺失範圍內的 SNV、包在大缺失裡的小缺失、兩個互相重疊的缺失、錨在缺失中間的插入）
    raise RebuildConflict —— 舊版會吃掉或截斷其中一個 allele。完全相同的重複 edit 只算一次。
    """
    out = []
    cursor = span_start                       # 下一個要輸出的 ref 位置（1-based）
    span_end = span_start + len(ref_seq) - 1
    for pos, r, a in sorted(set(edits), key=_edit_order):
        if a == "*" or a.startswith("<"):     # spanning deletion / symbolic：不是鹼基序列
            raise RebuildConflict("symbolic ALT %s at %d" % (a, pos))
        if pos >= cursor:
            out.append(ref_seq[cursor - span_start: pos - span_start])   # 中間未變 ref
            out.append(a)
        else:
            k = cursor - pos                  # 這顆的 ref 已被前面處理過的長度
            if k > len(r) or k > len(a) or a[:k] != r[:k]:
                raise RebuildConflict("overlap at %d: %s>%s" % (pos, r, a))
            out.append(a[k:])                 # 其餘部分（多半是插入的序列）
        cursor = pos + len(r)                 # 重疊時 k ≤ len(r)，游標不會倒退
    if cursor <= span_end:
        out.append(ref_seq[cursor - span_start:])
    return "".join(out)


def reconstruct(cluster: List[Var], fetch: Callable[[str, int, int], str]):
    """回傳 (pos, ref, alt_list, gt_list) 或 None（無非參考、無法合）；
    重疊互相矛盾時 build_hap 會 raise RebuildConflict。"""
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


# ─────────────────────────────────────────────────────────────
# 合不合成的決定（process() 與診斷腳本共用同一份邏輯）
# ─────────────────────────────────────────────────────────────
def _phase_unknown(cluster: List[Var]) -> bool:
    """叢集內有 ≥2 顆雙套 het，但它們不全是「同一個 phase set 的 phased」→ 相對 phase 未知。
    這種叢集只可能因足跡重疊而形成（不重疊的 het 要 phased 同 PS 才會連進來）；未 phase 的 GT 在
    reconstruct() 裡依位置都落在同一條單體，等於憑空假設 cis。phased 但沒有 PS 欄的，依 VCF 規格
    屬於同一個隱含的 phase set。"""
    hets = [v for v in cluster if len(v.alleles) == 2 and v.het_hap() is not None]
    if len(hets) < 2:
        return False
    if not all(v.phased for v in hets):
        return True
    return len({v.ps or "." for v in hets}) != 1


def plan_cluster(cluster: List[Var], fetch: Callable, is_mito: bool, sample_col: int = 0):
    """決定一個多顆的叢集要不要合成。回傳 (res, anchor, reason)：
    reason == "merged" 時 res = (pos, ref, alt_list, gt) 為合成結果；其餘 reason 代表原封通過：
      phase_unknown     ≥2 顆 het 的相對 phase 未知（上方 (e)）
      overlap_conflict  重疊互相矛盾，或含 * / symbolic ALT（(f)）
      multi             重建後有 2 個 ALT，如 1|2（(a)）
      no_anchor         找不到 biallelic、有 FORMAT 的 anchor（(b)）
      mixed_or_mito     混 ploidy，或 chrM 的 haploid 叢集（(c)(d)）
      no_alt            重建結果兩條都是 ref（正常不會發生）"""
    ploidies = {len(v.alleles) for v in cluster}
    anchor = _fmt_anchor(cluster, sample_col)
    try:
        if ploidies == {2}:                            # 全 diploid → 雙單體重建
            if _phase_unknown(cluster):
                return None, anchor, "phase_unknown"
            res = reconstruct(cluster, fetch)
            if res is not None and len(res[2]) > 1:
                return None, anchor, "multi"
        elif ploidies == {1} and not is_mito:          # 全 haploid（非 chrM）→ 單套重建
            res = reconstruct_haploid(cluster, fetch)
        else:                                          # 混 ploidy、chrM haploid、其他 → 不合
            return None, anchor, "mixed_or_mito"
    except RebuildConflict:
        return None, anchor, "overlap_conflict"
    if res is None:
        return None, anchor, "no_alt"
    if anchor is None:
        return None, anchor, "no_anchor"
    return res, anchor, "merged"


def process(in_vcf: str, out_vcf: str, fetch: Callable, max_gap: int,
            sample_col: int = 0) -> dict:
    """主流程；回傳統計。fetch 可注入（測試用）。"""
    stats = {"clusters_merged": 0, "clusters_haploid": 0, "clusters_fallback": 0,
             "clusters_phase_unknown": 0, "clusters_overlap_conflict": 0,
             "records_in": 0, "records_out": 0, "records_nocall": 0}
    header, chrom_vars, order_chrom = [], {}, []
    chrom_nocall = {}                  # chrom -> [(pos, line)]：沒有 ALT、不參與叢集的原行
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
            recs = list(chrom_nocall.get(chrom, []))   # 沒有 ALT 的紀錄：原行輸出
            for cl in clusters:
                if len(cl) == 1:                       # 孤立顆 → 原行輸出，完全不動
                    recs.append((cl[0].pos, cl[0].line))
                    continue
                res, anchor, reason = plan_cluster(cl, fetch, is_mito, sample_col)
                # phase 未知 / 重疊矛盾 / 多 ALT（1|2）/ 無 anchor / 混 ploidy → 原封通過（保留 AD）。
                if reason != "merged":
                    for v in cl:
                        recs.append((v.pos, v.line))
                    stats["clusters_fallback"] += 1
                    if reason in ("phase_unknown", "overlap_conflict"):
                        stats["clusters_" + reason] += 1
                    continue
                pos, ref, alt_list, gt = res
                # 輸出前最小化（見檔頭「輸出與相依」）：DV/HC 的前導鹼基長度不同時，才會對得上同一個 POS。
                pos, ref, alt_list = trim_alleles(pos, ref, alt_list)
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
            if chrom not in chrom_vars:
                chrom_vars[chrom] = []
                order_chrom.append(chrom)
            # 這個 sample 沒有 ALT（./.、0/0、單套 0 或 .，如 DV 的 RefCall）→ 不參與叢集，
            # 原行輸出。否則被否決的候選會當 anchor（FORMAT/QUAL/FILTER 被換成它的）、
            # 把合成紀錄的 POS/REF 撐寬、或把兩顆不相干的 call 串成同一叢（見檔頭 ⚠️）。
            if not any(a > 0 for a in alleles):
                chrom_nocall.setdefault(chrom, []).append((pos, line))
                stats["records_nocall"] += 1
                continue
            # 最小化後再參與叢集/重建（移除 padding 假性重疊）；passthrough 仍用原始 line。
            tpos, tref, talts = trim_alleles(pos, ref, alt.split(","))
            v = Var(chrom, tpos, tref, talts, line,
                    alleles=alleles, phased=phased, ps=ps)
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
        "passthrough_clusters=%d nocall_passthrough=%d phase_unknown=%d overlap_conflict=%d\n"
        % (st["records_in"], st["records_out"], st["clusters_merged"],
           st["clusters_haploid"], st["clusters_fallback"], st["records_nocall"],
           st["clusters_phase_unknown"], st["clusters_overlap_conflict"]))


if __name__ == "__main__":
    main()
