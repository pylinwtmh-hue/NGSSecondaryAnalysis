#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for ploidy_check.py — dependency-free (stdlib only).
Run:  python3 scripts/test_ploidy_check.py
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ploidy_check as P  # noqa: E402


def _auto(mean):
    return {c: mean for c in P.AUTOSOMES}


def test_normalize_declared_sex():
    for raw in ("male", "M", "XY", "1"):
        assert P.normalize_declared_sex(raw) == "XY", raw
    for raw in ("female", "f", "xx", "2"):
        assert P.normalize_declared_sex(raw) == "XX", raw
    for raw in ("unknown", "", None, "weird"):
        assert P.normalize_declared_sex(raw) == "unknown", raw
    print("PASS test_normalize_declared_sex")


def test_male_xy():
    m = _auto(42.8); m["chrX"] = 21.4; m["chrY"] = 21.0
    res = P.analyze(m, {}, "male", "WGS")
    assert res["estimated_karyotype"] == "XY", res["estimated_karyotype"]
    assert res["warnings"] == [], res["warnings"]        # declared male matches
    print("PASS test_male_xy -> XY, no warnings")


def test_female_xx():
    m = _auto(42.8); m["chrX"] = 42.8; m["chrY"] = 0.4
    res = P.analyze(m, {}, "female", "WGS")
    assert res["estimated_karyotype"] == "XX", res["estimated_karyotype"]
    assert res["warnings"] == [], res["warnings"]
    print("PASS test_female_xx -> XX, no warnings")


def test_sex_mismatch():
    # data looks male, samplesheet says female -> WARN
    m = _auto(42.8); m["chrX"] = 21.0; m["chrY"] = 20.5
    res = P.analyze(m, {}, "female", "WGS")
    assert res["estimated_karyotype"] == "XY"
    assert any("SEX MISMATCH" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_sex_mismatch -> flagged")


def test_trisomy21():
    m = _auto(42.8); m["chr21"] = 42.8 * 1.5; m["chrX"] = 42.8; m["chrY"] = 0.3
    res = P.analyze(m, {}, "female", "WGS")
    assert any(c == "chr21" for c, _ in res["aneuploidy"]), res["aneuploidy"]
    assert any("chr21" in w and "非整倍體" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_trisomy21 -> chr21 flagged")


def test_xxy_klinefelter():
    m = _auto(42.8); m["chrX"] = 42.8; m["chrY"] = 21.0   # 2 copies X + Y present
    res = P.analyze(m, {}, "male", "WGS")
    assert res["estimated_karyotype"] == "XXY?", res["estimated_karyotype"]
    assert any("XXY" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_xxy_klinefelter -> XXY? flagged")


def test_normal_no_false_aneuploidy():
    # ±5% noise on autosomes must NOT trip the aneuploidy flag
    m = {c: 42.8 * (1.0 + (0.05 if i % 2 else -0.05)) for i, c in enumerate(P.AUTOSOMES)}
    m["chrX"] = 42.8; m["chrY"] = 0.3
    res = P.analyze(m, {}, "female", "WGS")
    assert res["aneuploidy"] == [], res["aneuploidy"]
    print("PASS test_normal_no_false_aneuploidy")


def test_parse_prefers_region():
    txt = ("chrom\tlength\tbases\tmean\tmin\tmax\n"
           "chr1\t1000\t5000\t5.0\t0\t99\n"
           "chr1_region\t900\t8100\t9.0\t0\t99\n"     # region preferred over whole
           "chrX\t500\t1000\t2.0\t0\t99\n"            # no region -> fall back to whole
           "total\t0\t0\t0\t0\t0\n")
    d = tempfile.mkdtemp()
    p = os.path.join(d, "s.summary.txt")
    open(p, "w").write(txt)
    means, length = P.parse_mosdepth_summary(p)
    assert means["chr1"] == 9.0, means           # region wins
    assert means["chrX"] == 2.0, means           # whole-chrom fallback
    assert "total" not in means
    assert length["chr1"] == 1000
    print("PASS test_parse_prefers_region")


if __name__ == "__main__":
    test_normalize_declared_sex()
    test_male_xy()
    test_female_xx()
    test_sex_mismatch()
    test_trisomy21()
    test_xxy_klinefelter()
    test_normal_no_false_aneuploidy()
    test_parse_prefers_region()
    print("\nALL TESTS PASSED")
