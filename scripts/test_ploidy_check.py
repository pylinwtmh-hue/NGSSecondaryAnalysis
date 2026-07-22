#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for ploidy_check.py — dependency-free (stdlib only).
Run:  python3 scripts/test_ploidy_check.py

NDC is normalized to the estimated karyotype (male chrX -> ~1.0, matching DRAGEN);
RATIO keeps the raw autosome ratio (male chrX -> ~0.5) as sex/dosage evidence.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ploidy_check as P  # noqa: E402


def _auto(mean):
    return {c: mean for c in P.AUTOSOMES}


def _approx(a, b, tol=0.03):
    return a is not None and abs(a - b) <= tol


def test_normalize_declared_sex():
    for raw in ("male", "M", "XY", "1"):
        assert P.normalize_declared_sex(raw) == "XY", raw
    for raw in ("female", "f", "xx", "2"):
        assert P.normalize_declared_sex(raw) == "XX", raw
    for raw in ("unknown", "", None, "weird"):
        assert P.normalize_declared_sex(raw) == "unknown", raw
    print("PASS test_normalize_declared_sex")


def test_male_ndc_normalized_ratio_raw():
    # THE point of option B: male chrX NDC ~1.0 (matches DRAGEN), RATIO ~0.5 (raw evidence)
    m = _auto(42.8); m["chrX"] = 21.4; m["chrY"] = 21.0
    res = P.analyze(m, {}, "male", "WGS")
    assert res["estimated_karyotype"] == "XY"
    assert _approx(res["ndc"]["chrX"], 1.0), res["ndc"]["chrX"]
    assert _approx(res["ndc"]["chrY"], 1.0, 0.05), res["ndc"]["chrY"]
    assert _approx(res["ratio"]["chrX"], 0.5), res["ratio"]["chrX"]
    assert res["aneuploidy"] == [] and res["warnings"] == []
    print("PASS test_male_ndc_normalized_ratio_raw -> chrX NDC~1.0, RATIO~0.5")


def test_female_xx():
    m = _auto(42.8); m["chrX"] = 42.8; m["chrY"] = 0.4
    res = P.analyze(m, {}, "female", "WGS")
    assert res["estimated_karyotype"] == "XX"
    assert _approx(res["ndc"]["chrX"], 1.0)
    assert res["ndc"]["chrY"] is None, "female chrY expected 0 -> NDC NA (not flagged)"
    assert res["aneuploidy"] == [] and res["warnings"] == []
    print("PASS test_female_xx -> chrX NDC~1.0, chrY NDC NA")


def test_sex_mismatch():
    m = _auto(42.8); m["chrX"] = 21.0; m["chrY"] = 20.5
    res = P.analyze(m, {}, "female", "WGS")     # data male, declared female
    assert res["estimated_karyotype"] == "XY"
    assert any("SEX MISMATCH" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_sex_mismatch -> flagged")


def test_trisomy21():
    m = _auto(42.8); m["chr21"] = 42.8 * 1.5; m["chrX"] = 42.8; m["chrY"] = 0.3
    res = P.analyze(m, {}, "female", "WGS")
    assert any(c == "chr21" for c, _ in res["aneuploidy"]), res["aneuploidy"]
    assert _approx(res["ndc"]["chr21"], 1.5, 0.05)
    print("PASS test_trisomy21 -> chr21 NDC~1.5 flagged")


def test_xxy_klinefelter():
    # 2 copies X + Y present -> XXY?; NDC self-normalizes to ~1.0, karyotype warning fires
    m = _auto(42.8); m["chrX"] = 42.8; m["chrY"] = 21.0
    res = P.analyze(m, {}, "male", "WGS")
    assert res["estimated_karyotype"] == "XXY?"
    assert _approx(res["ndc"]["chrX"], 1.0), res["ndc"]["chrX"]
    assert res["aneuploidy"] == [], "sex-chrom aneuploidy shows as karyotype, not NDC flag"
    assert any("XXY" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_xxy_klinefelter -> XXY? flagged via karyotype")


def test_normal_no_false_aneuploidy():
    m = {c: 42.8 * (1.0 + (0.05 if i % 2 else -0.05)) for i, c in enumerate(P.AUTOSOMES)}
    m["chrX"] = 42.8; m["chrY"] = 0.3
    res = P.analyze(m, {}, "female", "WGS")
    assert res["aneuploidy"] == [], res["aneuploidy"]
    print("PASS test_normal_no_false_aneuploidy")


def test_header_keys_match_dragen():
    # VCF header must use the SAME keys as DRAGEN (##estimatedSexKaryotype/##referenceSexKaryotype)
    m = _auto(42.8); m["chrX"] = 21.4; m["chrY"] = 21.0
    res = P.analyze(m, {"chr1": 100, "chrX": 200}, "male", "WGS")
    d = tempfile.mkdtemp(); vcf = os.path.join(d, "s.ploidy.vcf")
    P.write_vcf(vcf, "S1", res, "WGS")
    txt = open(vcf).read()
    assert "##estimatedSexKaryotype=XY" in txt
    assert "##referenceSexKaryotype=XY" in txt
    assert "DC:NDC:RATIO" in txt, "FORMAT must expose DC, NDC and raw RATIO"
    # chrX line: NDC ~1.0, RATIO ~0.5
    xline = [ln for ln in txt.splitlines() if ln.startswith("chrX")][0].split("\t")
    dc, ndc, ratio = xline[9].split(":")
    assert _approx(float(ndc), 1.0) and _approx(float(ratio), 0.5), xline[9]
    print("PASS test_header_keys_match_dragen -> unified header + DC:NDC:RATIO")


def test_parse_prefers_region():
    txt = ("chrom\tlength\tbases\tmean\tmin\tmax\n"
           "chr1\t1000\t5000\t5.0\t0\t99\n"
           "chr1_region\t900\t8100\t9.0\t0\t99\n"
           "chrX\t500\t1000\t2.0\t0\t99\n"
           "total\t0\t0\t0\t0\t0\n")
    d = tempfile.mkdtemp(); p = os.path.join(d, "s.summary.txt")
    open(p, "w").write(txt)
    means, length = P.parse_mosdepth_summary(p)
    assert means["chr1"] == 9.0 and means["chrX"] == 2.0 and "total" not in means
    assert length["chr1"] == 1000
    print("PASS test_parse_prefers_region")


if __name__ == "__main__":
    test_normalize_declared_sex()
    test_male_ndc_normalized_ratio_raw()
    test_female_xx()
    test_sex_mismatch()
    test_trisomy21()
    test_xxy_klinefelter()
    test_normal_no_false_aneuploidy()
    test_header_keys_match_dragen()
    test_parse_prefers_region()
    print("\nALL TESTS PASSED")
