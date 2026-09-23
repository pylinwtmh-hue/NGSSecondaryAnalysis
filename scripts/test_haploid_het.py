#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for haploid_het.awk — stdlib only (needs an `awk` on PATH).
Run:  python3 scripts/test_haploid_het.py
      AWK="busybox awk" python3 scripts/test_haploid_het.py   # the container's awk

Male haploid regions come from the real assets/sex_ploidy_GRCh38.txt (the same file
+fixploidy reads), so a PAR coordinate change there is picked up here too.

Pins:
  - chrX non-PAR het 0/k -> k/k (separator kept), INFO/HAPLOID_HET=<callers>
  - two different ALTs (j/k): flagged, GT untouched
  - chrY het -> ./.
  - PAR1 / PAR2, autosomes, chrM, hom, REF, missing, half-missing: untouched
  - other FORMAT fields and existing INFO kept; header gains the INFO definition
  - SEX=F: no data line changes
"""
import os
import shlex
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
AWK_SCRIPT = os.path.join(HERE, "haploid_het.awk")
PLOIDY = os.path.join(HERE, "..", "assets", "sex_ploidy_GRCh38.txt")
AWK = shlex.split(os.environ.get("AWK", "awk"))

HEADER = [
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVAL55_DV\tVAL55_HC",
]

# (chrom, pos, alt, info, dv, hc)
RECORDS = [
    ("chrX", 1000000,   "T",   ".",          "0/1:9,8:.",       "./.:.:."),      # PAR1
    ("chrX", 5000000,   "T",   ".",          "0|1:9,8:4999000", "1|0:7,6:4999000"),
    ("chrX", 5000100,   "T",   "COMBINED=2", "0/1:9,8:.",       "1/1:0,15:."),
    ("chrX", 5000200,   "T,G", ".",          "0/1:9,8,.:.",     "0/2:7,.,6:."),  # fused
    ("chrX", 5000300,   "T,G", ".",          "1/2:1,8,7:.",     "./.:.:."),
    ("chrX", 5000400,   "T",   ".",          "./.:.:.",         "0/0:12,0:."),
    ("chrX", 5000500,   "T",   ".",          "./1:.,8:.",       "./.:.:."),
    ("chrX", 155800000, "T",   ".",          "./.:.:.",         "0/1:9,8:."),    # PAR2
    ("chrX", 156035000, "T",   ".",          "./.:.:.",         "0/1:9,8:."),    # ploidy-1 tail
    ("chrY", 3000000,   "T",   ".",          "0/1:9,8:.",       "./.:.:."),
    ("chrY", 3000100,   "T",   ".",          "1/1:0,9:.",       "./.:.:."),
    ("chr1", 100,       "T",   ".",          "0/1:9,8:.",       "0/1:9,8:."),
    ("chrM", 100,       "T",   ".",          "0/1:90,80:.",     "./.:.:."),
]


def _line(r):
    chrom, pos, alt, info, dv, hc = r
    return "\t".join([chrom, str(pos), ".", "A", alt, "50", "PASS", info, "GT:AD:PS", dv, hc])


def _run(sex):
    d = tempfile.mkdtemp()
    vcf = os.path.join(d, "in.vcf")
    with open(vcf, "w") as fh:
        fh.write("\n".join(HEADER + [_line(r) for r in RECORDS]) + "\n")
    p = subprocess.run(AWK + ["-v", "SEX=" + sex, "-f", AWK_SCRIPT, PLOIDY, vcf],
                       capture_output=True, text=True, check=True)
    lines = p.stdout.rstrip("\n").split("\n")
    head = [ln for ln in lines if ln.startswith("#")]
    body = {(f[0], int(f[1])): f for f in (ln.split("\t") for ln in lines if not ln.startswith("#"))}
    return head, body, p.stderr


def _gts(f):
    return f[9], f[10]


def test_male_chrx_nonpar():
    head, body, err = _run("M")
    f = body[("chrX", 5000000)]
    assert _gts(f) == ("1|1:9,8:4999000", "1|1:7,6:4999000"), _gts(f)
    assert f[7] == "HAPLOID_HET=DV,HC", f[7]
    f = body[("chrX", 5000100)]
    assert _gts(f) == ("1/1:9,8:.", "1/1:0,15:."), _gts(f)
    assert f[7] == "COMBINED=2;HAPLOID_HET=DV", f[7]          # existing INFO kept
    f = body[("chrX", 5000200)]
    assert _gts(f) == ("1/1:9,8,.:.", "2/2:7,.,6:."), _gts(f)   # each caller keeps its own ALT
    assert f[7] == "HAPLOID_HET=DV,HC", f[7]
    f = body[("chrX", 5000300)]
    assert _gts(f) == ("1/2:1,8,7:.", "./.:.:."), _gts(f)
    assert f[7] == "HAPLOID_HET=DV", f[7]
    f = body[("chrX", 156035000)]
    assert _gts(f)[1] == "1/1:9,8:." and f[7] == "HAPLOID_HET=HC", f
    assert "chrX_het_to_alt_DV=3" in err and "chrX_het_to_alt_HC=3" in err, err
    assert "chrX_multi_alt_het_DV=1" in err, err
    print("PASS test_male_chrx_nonpar -> het to ALT + HAPLOID_HET, AD/PS kept")


def test_male_untouched():
    head, body, _ = _run("M")
    for key in [("chrX", 1000000), ("chrX", 5000400), ("chrX", 5000500),
                ("chrX", 155800000), ("chrY", 3000100), ("chr1", 100), ("chrM", 100)]:
        orig = [r for r in RECORDS if (r[0], r[1]) == key][0]
        assert "\t".join(body[key]) == _line(orig), (key, body[key])
    print("PASS test_male_untouched -> PAR1/PAR2, missing, half-missing, hom, autosome, chrM unchanged")


def test_male_chry_het_to_missing():
    head, body, err = _run("M")
    f = body[("chrY", 3000000)]
    assert _gts(f) == ("./.:9,8:.", "./.:.:."), _gts(f)
    assert f[7] == ".", f[7]                                      # hidden, not flagged
    assert "chrY_het_to_missing_DV=1" in err, err
    print("PASS test_male_chry_het_to_missing -> chrY het becomes ./.")


def test_header():
    head, _, _ = _run("M")
    assert head[-1].startswith("#CHROM"), head[-1]
    assert head[-2].startswith("##INFO=<ID=HAPLOID_HET,Number=.,Type=String,"), head[-2]
    assert head[:-2] == HEADER[:-1], head
    print("PASS test_header -> INFO/HAPLOID_HET declared right before #CHROM")


def test_female_no_change():
    _, body, err = _run("F")
    for r in RECORDS:
        assert "\t".join(body[(r[0], r[1])]) == _line(r), r
    assert "haploid_regions=0" in err, err
    print("PASS test_female_no_change -> SEX=F leaves every record as is")


if __name__ == "__main__":
    test_male_chrx_nonpar()
    test_male_untouched()
    test_male_chry_het_to_missing()
    test_header()
    test_female_no_change()
    print("\nALL TESTS PASSED")
