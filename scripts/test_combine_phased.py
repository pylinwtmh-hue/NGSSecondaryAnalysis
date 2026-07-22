#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for combine_phased.py — dependency-free (stdlib only).
Run:  python3 scripts/test_combine_phased.py

Covers the cases discussed for the NCKUH compound-merging design:
  - footprint-based clustering (overlap beats POS distance)
  - het cis del+ins (SUZ12 HC representation) -> single MNV
  - two cis SNVs within gap; kept separate at max_gap=0
  - trans non-overlapping -> not merged
  - hom+hom -> 1|1 MNV
  - overlapping opposite-hap -> 1|2 co-representation
  - isolated -> passthrough

NOTE (known limitation, see module docstring): reconstruction of *overlapping*
edits on the SAME haplotype (padded/complex caller splits) is not universally
reliable from the VCF alone — e.g. DV's SUZ12 split (GAAA>GAA + A>T) is NOT
correctly rebuilt here. That case is handled by the chosen haplotype engine,
not this stdlib reconstructor; these tests pin only the well-defined behaviour.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combine_phased as C  # noqa: E402


def mkfetch(windows):
    def fetch(chrom, s, e):
        start, seq = windows[chrom]
        return seq[s - start: e - start + 1]
    return fetch


def _run_process(records, fetch, max_gap=2, sample_col=0):
    """把 records（完整 VCF 資料行）寫成暫存 VCF，跑 process()，回 (非#行, stats)。"""
    hdr = ["##fileformat=VCFv4.2",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
           '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">',
           '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
    d = tempfile.mkdtemp()
    ip, op = os.path.join(d, "in.vcf"), os.path.join(d, "out.vcf")
    with open(ip, "w") as fh:
        fh.write("\n".join(hdr + records) + "\n")
    st = C.process(ip, op, fetch, max_gap, sample_col)
    with open(op) as fh:
        out = [ln.rstrip("\n") for ln in fh if not ln.startswith("#")]
    return out, st


def test_suz12_hc():
    v1 = C.Var("chr17", 31998950, "GAAA", ["G"], "L1", alleles=[0, 1], phased=True, ps="P")
    v2 = C.Var("chr17", 31998953, "A", ["ATT"], "L2", alleles=[0, 1], phased=True, ps="P")
    fetch = mkfetch({"chr17": (31998950, "GAAA")})
    cls = C.cluster_vars([v1, v2], 1)
    assert len(cls) == 1, "SUZ12 HC records overlap at 31998953 -> must cluster"
    assert C.reconstruct(cls[0], fetch) == (31998950, "GAAA", ["GTT"], [0, 1])
    print("PASS test_suz12_hc -> GAAA>GTT 0|1")


def test_cis_two_snv_gap():
    v1 = C.Var("chr1", 100, "C", ["T"], "L", alleles=[0, 1], phased=True, ps="P")
    v2 = C.Var("chr1", 102, "G", ["A"], "L", alleles=[0, 1], phased=True, ps="P")
    fetch = mkfetch({"chr1": (100, "CAG")})
    assert C.reconstruct(C.cluster_vars([v1, v2], 1)[0], fetch) == (100, "CAG", ["TAA"], [0, 1])
    assert len(C.cluster_vars([v1, v2], 0)) == 2
    print("PASS test_cis_two_snv_gap -> CAG>TAA (max_gap>=1); separate at max_gap=0")


def test_trans_not_merged():
    v1 = C.Var("chr1", 100, "C", ["T"], "L", alleles=[0, 1], phased=True, ps="P")
    v2 = C.Var("chr1", 102, "G", ["A"], "L", alleles=[1, 0], phased=True, ps="P")
    assert len(C.cluster_vars([v1, v2], 5)) == 2
    print("PASS test_trans_not_merged")


def test_hom_hom():
    v1 = C.Var("chr2", 200, "A", ["G"], "L", alleles=[1, 1], phased=False, ps=None)
    v2 = C.Var("chr2", 201, "C", ["T"], "L", alleles=[1, 1], phased=False, ps=None)
    fetch = mkfetch({"chr2": (200, "AC")})
    assert C.reconstruct(C.cluster_vars([v1, v2], 1)[0], fetch) == (200, "AC", ["GT"], [1, 1])
    print("PASS test_hom_hom -> AC>GT 1|1")


def test_overlap_opposite_hap_1_2():
    v1 = C.Var("chr4", 400, "A", ["AT"], "L", alleles=[1, 0], phased=True, ps="P")
    v2 = C.Var("chr4", 400, "A", ["AG"], "L", alleles=[0, 1], phased=True, ps="P")
    fetch = mkfetch({"chr4": (400, "A")})
    assert C.reconstruct(C.cluster_vars([v1, v2], 0)[0], fetch) == (400, "A", ["AT", "AG"], [1, 2])
    print("PASS test_overlap_opposite_hap_1_2 -> A>AT,AG 1|2")


def test_isolated_passthrough():
    v1 = C.Var("chr1", 100, "C", ["T"], "L", alleles=[0, 1], phased=True, ps="P")
    v2 = C.Var("chr1", 500, "G", ["A"], "L", alleles=[0, 1], phased=True, ps="Q")
    cls = C.cluster_vars([v1, v2], 2)
    assert len(cls) == 2 and all(len(c) == 1 for c in cls)
    print("PASS test_isolated_passthrough")


def test_trim_alleles():
    # padding removal: GAAA>GAA -> GA>G (pos unchanged); SNV stays
    assert C.trim_alleles(31998950, "GAAA", ["GAA"]) == (31998950, "GA", ["G"])
    assert C.trim_alleles(31998953, "A", ["T"]) == (31998953, "A", ["T"])
    print("PASS test_trim_alleles")


def test_suz12_dv_after_trim():
    # DV raw: GAAA>GAA (del one A) + A>T (SNV); WITHOUT trim they falsely overlap
    # at 31998953 and the SNV is dropped. Trim first -> correct GAAA>GAT.
    fetch = mkfetch({"chr17": (31998950, "GAAA")})
    raw = [(31998950, "GAAA", ["GAA"]), (31998953, "A", ["T"])]
    vs = []
    for p, r, a in raw:
        tp, tr, ta = C.trim_alleles(p, r, a)
        vs.append(C.Var("chr17", tp, tr, ta, "L", alleles=[0, 1], phased=True, ps="P"))
    assert C.reconstruct(C.cluster_vars(vs, 2)[0], fetch) == (31998950, "GAAA", ["GAT"], [0, 1])
    print("PASS test_suz12_dv_after_trim -> GAAA>GAT (SNV preserved)")


def test_footprint_beats_pos_distance():
    v1 = C.Var("chr5", 1000, "ACGT", ["A"], "L", alleles=[0, 1], phased=True, ps="P")
    v2 = C.Var("chr5", 1003, "T", ["TA"], "L", alleles=[0, 1], phased=True, ps="P")
    assert len(C.cluster_vars([v1, v2], 0)) == 1, "overlap must merge regardless of POS distance"
    print("PASS test_footprint_beats_pos_distance")


def test_merged_keeps_format():
    # SUZ12-like het cis del+ins -> 0|1 MNV must INHERIT depth from the widest
    # (deletion) anchor, not drop it. Regression test for the DRAGEN AD-loss bug.
    fetch = mkfetch({"chr17": (31998950, "GAAA")})
    recs = [
        "chr17\t31998950\t.\tGAAA\tG\t60\tPASS\t.\tGT:AD:DP:AF\t0|1:30,12:42:0.29",
        "chr17\t31998953\t.\tA\tATT\t55\tPASS\t.\tGT:AD:DP:AF\t0|1:31,11:42:0.26",
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 1 and st["clusters_fallback"] == 0
    assert len(out) == 1, out
    f = out[0].split("\t")
    assert (f[3], f[4]) == ("GAAA", "GTT"), (f[3], f[4])
    assert f[5] == "60" and f[6] == "PASS", "QUAL/FILTER from anchor"
    assert "COMBINED=2" in f[7]
    d = dict(zip(f[8].split(":"), f[9].split(":")))
    assert d["GT"] == "0|1"
    assert d["AD"] == "30,12", d          # inherited from widest (deletion) anchor
    assert d["DP"] == "42" and d["AF"] == "0.29"
    print("PASS test_merged_keeps_format -> AD/DP/AF preserved on combined MNV")


def test_triallelic_passthrough_keeps_ad():
    # opposite-hap overlap reconstructs to 1|2 -> DO NOT fabricate 3-allele AD;
    # pass both source records through untouched (their AD survives for norm to split).
    fetch = mkfetch({"chr4": (400, "A")})
    recs = [
        "chr4\t400\t.\tA\tAT\t50\tPASS\t.\tGT:AD:DP\t1|0:10,5:15",
        "chr4\t400\t.\tA\tAG\t50\tPASS\t.\tGT:AD:DP\t0|1:9,6:15",
    ]
    out, st = _run_process(recs, fetch, max_gap=0)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 1
    assert len(out) == 2 and all("COMBINED" not in ln for ln in out)
    ads = sorted(ln.split("\t")[9].split(":")[1] for ln in out)
    assert ads == ["10,5", "9,6"] or ads == ["10,5", "9,6"][::-1], ads
    print("PASS test_triallelic_passthrough_keeps_ad -> 1|2 kept as 2 records w/ AD")


def test_haploid_passthrough_keeps_ad():
    # haploid (chrX/Y/M) components can't be diploid-reconstructed -> pass through.
    fetch = mkfetch({"chrX": (1000, "CG")})
    recs = [
        "chrX\t1000\t.\tCG\tC\t50\tPASS\t.\tGT:AD:DP\t1:8:8",   # haploid del, ftpt 1000-1001
        "chrX\t1001\t.\tG\tT\t50\tPASS\t.\tGT:AD:DP\t1:7:7",    # overlaps at 1001
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 1
    assert len(out) == 2 and all("COMBINED" not in ln for ln in out)
    print("PASS test_haploid_passthrough_keeps_ad -> haploid kept as-is")


def test_isolated_still_untouched():
    # a lone record must be byte-identical on the way out (no COMBINED, full FORMAT).
    fetch = mkfetch({"chr1": (100, "C")})
    recs = ["chr1\t100\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t0/1:20,18:38"]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 0
    assert out == recs, out
    print("PASS test_isolated_still_untouched -> lone record passthrough verbatim")


if __name__ == "__main__":
    test_suz12_hc()
    test_cis_two_snv_gap()
    test_trans_not_merged()
    test_hom_hom()
    test_overlap_opposite_hap_1_2()
    test_isolated_passthrough()
    test_trim_alleles()
    test_suz12_dv_after_trim()
    test_footprint_beats_pos_distance()
    test_merged_keeps_format()
    test_triallelic_passthrough_keeps_ad()
    test_haploid_passthrough_keeps_ad()
    test_isolated_still_untouched()
    print("\nALL TESTS PASSED")
