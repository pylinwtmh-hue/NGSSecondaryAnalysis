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
  - records with no ALT in the sample's GT (DeepVariant RefCall ./., 0/0, haploid 0/.)
    never join a cluster: not an anchor, no padding, no bridging
  - combined REF/ALT are written minimised, so two callers' representations of one
    event get the same POS (a pure deletion keeps its anchor base)
  - overlapping het calls whose relative phase is unknown (unphased, or different PS) are
    not merged; a merge that would swallow or truncate an allele (SNV inside a deletion,
    deletion inside a deletion, insertion anchored inside a deletion, * allele) is not
    written either -> the source records pass through untouched (VAL-10: 7,927 alleles lost)
  - an SNV and an indel at the same base are applied SNV first (no SNV swallowed)

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
    # reconstructed GAAA>GTT, written minimised (shared leading G trimmed)
    assert (f[1], f[3], f[4]) == ("31998951", "AAA", "TT"), (f[1], f[3], f[4])
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


def test_haploid_cis_merge():
    # male non-PAR chrX: two adjacent hemizygous SNVs on the single copy -> one haploid MNV.
    # GT must be hemizygous "1" (not 1|1, not phased) with NO PS, inheriting an anchor's AD.
    fetch = mkfetch({"chrX": (1000, "CAG")})
    recs = [
        "chrX\t1000\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t1:8:8",
        "chrX\t1002\t.\tG\tA\t50\tPASS\t.\tGT:AD:DP\t1:9:9",
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 1 and st["clusters_haploid"] == 1 and st["clusters_fallback"] == 0
    assert len(out) == 1, out
    f = out[0].split("\t")
    assert (f[3], f[4]) == ("CAG", "TAA"), (f[3], f[4])
    d = dict(zip(f[8].split(":"), f[9].split(":")))
    assert d["GT"] == "1", d                       # hemizygous single allele
    assert "PS" not in d, "no phase set on a single-copy call"
    assert d["AD"] == "8" and d["DP"] == "8"       # inherited from leftmost (widest-tie) anchor
    assert "COMBINED=2" in f[7]
    print("PASS test_haploid_cis_merge -> CAG>TAA hemizygous GT=1, AD kept")


def test_mixed_ploidy_passthrough():
    # haploid + diploid in one overlapping cluster -> don't mix models, passthrough (keep AD).
    fetch = mkfetch({"chrX": (1000, "CG")})
    recs = [
        "chrX\t1000\t.\tCG\tC\t50\tPASS\t.\tGT:AD:DP\t1:8:8",      # haploid, ftpt 1000-1001
        "chrX\t1001\t.\tG\tT\t50\tPASS\t.\tGT:AD:DP\t0/1:5,6:11",  # diploid, overlaps at 1001
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 1
    assert len(out) == 2 and all("COMBINED" not in ln for ln in out)
    print("PASS test_mixed_ploidy_passthrough -> mixed ploidy kept as-is")


def test_haploid_mito_passthrough():
    # chrM haploid adjacent variants -> NOT combined (multi-copy heteroplasmy), passthrough.
    fetch = mkfetch({"chrM": (300, "CAG")})
    recs = [
        "chrM\t300\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t1:80:100",
        "chrM\t302\t.\tG\tA\t50\tPASS\t.\tGT:AD:DP\t1:60:100",
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 1
    assert len(out) == 2 and all("COMBINED" not in ln for ln in out)
    print("PASS test_haploid_mito_passthrough -> chrM haploid kept as-is")


def test_isolated_still_untouched():
    # a lone record must be byte-identical on the way out (no COMBINED, full FORMAT).
    fetch = mkfetch({"chr1": (100, "C")})
    recs = ["chr1\t100\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t0/1:20,18:38"]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 0
    assert out == recs, out
    print("PASS test_isolated_still_untouched -> lone record passthrough verbatim")


def test_nocall_wider_candidate_not_anchor():
    # DeepVariant keeps candidates it REJECTED (FILTER=RefCall, GT ./.). A wider rejected
    # deletion overlapping a real SNV call used to join the cluster and become the anchor:
    # the SNV came out as CGTA>CGGA with the deletion's QUAL/FILTER/AD/VAF and a fake 0|1+PS,
    # at the deletion's POS (so it no longer lined up with HC's call of the same SNV).
    # Regression test (VAL55: 28,050 such RefCall-anchored records).
    fetch = mkfetch({"chrT": (10, "CGTA")})
    recs = [
        "chrT\t10\t.\tCGTA\tC\t0.8\tRefCall\t.\tGT:AD:DP\t./.:20,3:23",   # rejected, widest
        "chrT\t12\t.\tT\tG\t40.1\tPASS\t.\tGT:AD:DP\t0/1:11,12:23",       # real call
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["records_nocall"] == 1, st
    assert all("COMBINED" not in ln for ln in out), out
    assert sorted(out) == sorted(recs), out        # both records verbatim, SNV keeps its own AD
    print("PASS test_nocall_wider_candidate_not_anchor -> rejected candidate passes through")


def test_nocall_does_not_bridge():
    # A rejected candidate spanning two real cis SNVs must not glue them into one MNV:
    # the SNVs are 4 bp apart (> max_gap=2), so on their own they stay separate.
    fetch = mkfetch({"chr1": (100, "CAGGTA")})
    recs = [
        "chr1\t100\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:100",
        "chr1\t100\t.\tCAGGTA\tC\t1.2\tRefCall\t.\tGT:AD:DP:PS\t0/0:17,2:19:.",
        "chr1\t105\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:9,10:19:100",
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["records_nocall"] == 1, st
    assert all("COMBINED" not in ln for ln in out), out
    assert sorted(out) == sorted(recs), out
    print("PASS test_nocall_does_not_bridge -> no MNV built across a rejected candidate")


def test_merged_output_minimised():
    # One caller splits an event into an insertion (with its anchor base) + an SNV, the other
    # reports it as one minimal record. Before: combine wrote the cluster span 100 AC>ATG, which
    # never lined up with the other caller's 101 C>TG at merge (VAL55: 820 split variants).
    fetch = mkfetch({"chrT": (100, "AC")})
    split_caller = [
        "chrT\t100\t.\tA\tAT\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:9,8:17:100",
        "chrT\t101\t.\tC\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:9,8:17:100",
    ]
    minimal_caller = ["chrT\t101\t.\tC\tTG\t60\tPASS\t.\tGT:AD:DP\t0/1:10,9:19"]
    out_a, st_a = _run_process(split_caller, fetch, max_gap=2)
    out_b, _ = _run_process(minimal_caller, fetch, max_gap=2)
    assert st_a["clusters_merged"] == 1 and len(out_a) == 1, out_a
    a, b = out_a[0].split("\t"), out_b[0].split("\t")
    assert (a[1], a[3], a[4]) == ("101", "C", "TG"), a[:5]
    assert (a[1], a[3], a[4]) == (b[1], b[3], b[4]), (a[:5], b[:5])

    # a pure-deletion result keeps its VCF anchor base (never trimmed to an empty allele)
    fetch = mkfetch({"chr17": (950, "GAAA")})
    dels = [
        "chr17\t950\t.\tGA\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:950",
        "chr17\t952\t.\tAA\tA\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:950",
    ]
    out, st = _run_process(dels, fetch, max_gap=2)
    assert st["clusters_merged"] == 1 and len(out) == 1, out
    f = out[0].split("\t")
    assert (f[1], f[3], f[4]) == ("950", "GAA", "G"), f[:5]
    print("PASS test_merged_output_minimised -> 100 AC>ATG written as 101 C>TG; deletion keeps anchor")


def test_haploid_nocall_passthrough():
    # haploid no-calls ("0" / ".") never join a cluster either; the real hemizygous
    # call next to them stays a lone, untouched record.
    fetch = mkfetch({"chrX": (1000, "CAG")})
    recs = [
        "chrX\t1000\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t1:0,8:8",
        "chrX\t1001\t.\tA\tG\t3\tRefCall\t.\tGT:AD:DP\t0:7,1:8",
        "chrX\t1002\t.\tG\tA\t2\tRefCall\t.\tGT:AD:DP\t.:6,2:8",
    ]
    out, st = _run_process(recs, fetch, max_gap=2)
    assert st["clusters_merged"] == 0 and st["records_nocall"] == 2, st
    assert sorted(out) == sorted(recs), out
    print("PASS test_haploid_nocall_passthrough -> haploid 0 / . not clustered")


def _passthrough(recs, fetch, reason, max_gap=2):
    """recs 必須原封輸出、不合成，且統計歸在 reason（phase_unknown / overlap_conflict）。"""
    out, st = _run_process(recs, fetch, max_gap=max_gap)
    assert st["clusters_merged"] == 0 and st["clusters_fallback"] == 1, st
    assert st["clusters_" + reason] == 1, st
    assert sorted(out) == sorted(recs), out
    return st


def test_unphased_overlap_not_merged():
    # VAL-10 chr4:115927671: two overlapping het deletions without PS. Old combine put both on
    # one haplotype and wrote a 5-bp deletion (CTGTTT>C 0|1 + fake PS) that neither call made.
    fetch = mkfetch({"chr4": (100, "CTGTTTA")})
    _passthrough(["chr4\t100\t.\tCTGT\tC\t40\tPASS\t.\tGT:AD:DP\t0/1:20,8:28",
                  "chr4\t102\t.\tGTTT\tG\t40\tPASS\t.\tGT:AD:DP\t0/1:20,9:29"],
                 fetch, "phase_unknown")
    # same pair phased but in two different phase sets: relative phase still unknown
    _passthrough(["chr4\t100\t.\tCTGT\tC\t40\tPASS\t.\tGT:AD:DP:PS\t0|1:20,8:28:100",
                  "chr4\t102\t.\tGTTT\tG\t40\tPASS\t.\tGT:AD:DP:PS\t0|1:20,9:29:102"],
                 fetch, "phase_unknown")
    print("PASS test_unphased_overlap_not_merged -> no 5-bp deletion invented")


def test_unphased_delins_not_merged():
    # The SUZ12 del+ins shape, but with no phase: could be cis (delins) or trans (del on one
    # copy, ins on the other). Not merged; with a common PS it still is (test_merged_keeps_format).
    fetch = mkfetch({"chr17": (31998950, "GAAA")})
    _passthrough(["chr17\t31998950\t.\tGAAA\tG\t60\tPASS\t.\tGT:AD:DP\t0/1:30,12:42",
                  "chr17\t31998953\t.\tA\tATT\t55\tPASS\t.\tGT:AD:DP\t0/1:31,11:42"],
                 fetch, "phase_unknown")
    print("PASS test_unphased_delins_not_merged -> unphased del+ins kept as 2 records")


def test_snv_inside_deletion_not_swallowed():
    # An SNV inside a deletion cannot sit on the same copy. Old combine wrote only the deletion
    # and the SNV vanished from the report (VAL-10 chr1:1746439, chr1:2981045 ...).
    fetch = mkfetch({"chr1": (500, "TTCAT")})
    _passthrough(["chr1\t500\t.\tTTCAT\tT\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:500",
                  "chr1\t502\t.\tC\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:500"],
                 fetch, "overlap_conflict")
    print("PASS test_snv_inside_deletion_not_swallowed")


def test_deletion_inside_hom_deletion_not_swallowed():
    # VAL-10 chr1:94814: hom 11-bp deletion + het 1-bp deletion inside it (contradictory calls).
    # Old combine returned the hom deletion alone; the het call disappeared.
    fetch = mkfetch({"chr1": (94814, "CTTTTCTTTTCT")})
    _passthrough(["chr1\t94814\t.\tCTTTTCTTTTCT\tC\t50\tPASS\t.\tGT:AD:DP\t1/1:0,20:20",
                  "chr1\t94824\t.\tCT\tC\t30\tPASS\t.\tGT:AD:DP\t0/1:9,8:17"],
                 fetch, "overlap_conflict")
    print("PASS test_deletion_inside_hom_deletion_not_swallowed")


def test_insertion_inside_deletion_not_truncated():
    # An insertion anchored in the middle of a deletion: old combine appended alt[k:] and
    # dropped (here all of) the inserted sequence.
    fetch = mkfetch({"chr1": (1300, "CGTTTC")})
    _passthrough(["chr1\t1300\t.\tCGTTTC\tC\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:1300",
                  "chr1\t1302\t.\tT\tTCCC\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:1300"],
                 fetch, "overlap_conflict")
    print("PASS test_insertion_inside_deletion_not_truncated")


def test_star_allele_not_rebuilt():
    # A '*' (spanning-deletion) allele is not sequence; never write it into a haplotype.
    fetch = mkfetch({"chr1": (900, "ACG")})
    _passthrough(["chr1\t900\t.\tACG\tA\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:900",
                  "chr1\t901\t.\tC\tT,*\t50\tPASS\t.\tGT:AD:DP:PS\t1|2:1,9,9:19:900"],
                 fetch, "overlap_conflict")
    print("PASS test_star_allele_not_rebuilt")


def test_snv_and_indel_at_same_base():
    # Same copy, same base: SNV A>G + insertion A>AT. Old combine sorted by string ("AT" < "G"),
    # applied the insertion first and swallowed the SNV (A>AT). SNV first -> A>GT.
    fetch = mkfetch({"chr1": (700, "A")})
    out, st = _run_process(["chr1\t700\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:700",
                            "chr1\t700\t.\tA\tAT\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:700"],
                           fetch, max_gap=2)
    assert st["clusters_merged"] == 1 and len(out) == 1, (st, out)
    f = out[0].split("\t")
    assert (f[1], f[3], f[4]) == ("700", "A", "GT"), f[:5]
    # SNV + deletion at the same base: the deletion's anchor is rewritten, the rest deleted
    fetch = mkfetch({"chr1": (800, "AC")})
    out, st = _run_process(["chr1\t800\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:800",
                            "chr1\t800\t.\tAC\tA\t50\tPASS\t.\tGT:AD:DP:PS\t0|1:10,9:19:800"],
                           fetch, max_gap=2)
    assert st["clusters_merged"] == 1 and len(out) == 1, (st, out)
    f = out[0].split("\t")
    assert (f[1], f[3], f[4]) == ("800", "AC", "G"), f[:5]
    print("PASS test_snv_and_indel_at_same_base -> A>GT, AC>G")


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
    test_haploid_cis_merge()
    test_mixed_ploidy_passthrough()
    test_haploid_mito_passthrough()
    test_isolated_still_untouched()
    test_nocall_wider_candidate_not_anchor()
    test_nocall_does_not_bridge()
    test_haploid_nocall_passthrough()
    test_merged_output_minimised()
    test_unphased_overlap_not_merged()
    test_unphased_delins_not_merged()
    test_snv_inside_deletion_not_swallowed()
    test_deletion_inside_hom_deletion_not_swallowed()
    test_insertion_inside_deletion_not_truncated()
    test_star_allele_not_rebuilt()
    test_snv_and_indel_at_same_base()
    print("\nALL TESTS PASSED")
