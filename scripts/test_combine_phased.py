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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combine_phased as C  # noqa: E402


def mkfetch(windows):
    def fetch(chrom, s, e):
        start, seq = windows[chrom]
        return seq[s - start: e - start + 1]
    return fetch


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
    print("\nALL TESTS PASSED")
