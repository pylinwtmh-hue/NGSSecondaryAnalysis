# CLAUDE.md — NGS Secondary Analysis Pipeline

Guidance for Claude Code and developers working in this repository.

## Overview

WGS/WES **germline secondary analysis**: FASTQ → BAM → SNV/indel, CNV/SV, STR, mtDNA, ROH.
GPU-accelerated with NVIDIA Clara Parabricks; orchestrated with Nextflow (DSL2) + Apptainer.

- `main.nf` — the single clinical pipeline (the old `main_research.nf` is retired; optional
  research tools are now flags on `main.nf`).
- `main_pon.nf` + `modules/pon.nf` + `nextflow_pon.config` — builds the gCNV + CNVkit
  **Panel of Normals (PON)**. Run once; case runs reuse the model.
- `modules/*.nf` — one file per stage (preprocessing, alignment, variant_calling, cnv_sv,
  repeat, mitochondria, alignment_qc, postprocessing, roh).
- `nextflow_main.config` — profiles `local` / `dgm` / `dgx` / `dgx_single` + params.

### Optional callers (flags)

| Flag | Tool | License | Default |
|------|------|---------|---------|
| `--run_roh` | bcftools roh | MIT/GPL ✅ commercial | **on** (lab decision) |
| `--run_automap` | AutoMap | none published ⚠️ | off |
| `--run_manta` | Manta | PolyForm Strict ⚠️ non-commercial | off |
| `--run_expansionhunter` | ExpansionHunter | PolyForm Strict ⚠️ non-commercial | off |

The default clinical path uses only commercially-usable tools (Delly BSD-3, GangSTR GPL,
bcftools roh MIT/GPL, GATK/Parabricks/fastp/samtools/bcftools/mosdepth/CNVkit).

## Run

```bash
# case (per-sample). samplesheet columns: sample,fastq_1,fastq_2,sex[,lane]
nextflow -c nextflow_main.config run main.nf -profile local \
    --input_csv samplesheet.csv --seq_type WES --run_gcnv true --out_dir <out>
# syntax check: append  --input_csv /dev/null -preview

# Panel of Normals (run once; rebuild after changing any gCNV hyperparameter — see below)
nextflow -c nextflow_pon.config run main_pon.nf \
    --input_csv pon_samplesheet.csv --pon_out_dir <pon_dir>
```

Output tree: `01_preprocessing 02_alignment 03_alignment_qc 04_snv_indel 05_cnv_sv
06_repeat 07_mitochondria 08_roh`.

---

## Evaluation feedback (2026-07) and this round of changes

Two issues were raised by the department after the first evaluation:

1. **Delly emitted far too many SVs.** `delly call` publishes every SV (PASS + LowQual).
   Fix: publish only `FILTER=PASS` (delly call already sets PASS from PE≥3 & MAPQ≥20).
   `delly filter -f germline` is **not usable on a single sample** (it needs ≥10 samples to
   compare depth ratios), so PASS-filtering is the correct single-sample approach. The
   clinically-relevant events were confirmed present in the PASS set.

2. **CNV under-called vs other platforms.** Diagnosis: the missing events are absent from
   CNVkit's raw `.cns` (not dropped by `--filter cn`) but **are** caught by Delly PASS →
   they are SV-type events that depth-based CNVkit inherently misses; Delly covers them.
   Decision: stop ad-hoc tuning, align the CNV/gCNV settings to **Broad's published
   germline-CNV defaults**, and **rebuild the PON**.

### Broad alignment applied (secondary)

| Area | file | before | after | rationale |
|------|------|--------|-------|-----------|
| Delly output | `modules/cnv_sv.nf` | all calls | `FILTER=PASS` only | cut LowQual noise |
| CNVkit call | `modules/cnv_sv.nf` | `--filter cn` | (removed) | over-aggressive for germline |
| VQSR SNP | `modules/variant_calling.nf` | no `DP` | `-an DP` | Broad WGS SNP recommendation |
| CPU HaplotypeCaller | `modules/variant_calling.nf` | `--BQSR` (removed in GATK4) | `ApplyBQSR` step | gatk#6041 (dead fallback path) |
| mtDNA filter | `modules/mitochondria.nf` | (no change) | `--mitochondria-mode` + blacklist mask | `--autosomal-coverage` was **tried but reverted**: it was removed from GATK 4.6 `FilterMutectCalls` (errors out) and Broad's current mito WDL doesn't use it. NuMT filtering = mitochondria-mode + blacklist mask |
| Alignment | `modules/alignment.nf` | `-Y` | `-Y -K 100000000` | thread-deterministic bwa |
| gCNV hyperparams | `nextflow_pon.config` + `modules/pon.nf` | p-alt 1e-3, coherence 1000 | **p-alt 5e-4, coherence 10000, p-active 1e-1** | Broad germline-CNV WDL defaults |

gCNV hyperparameters are now **parameterised in `nextflow_pon.config`** (`gcnv_p_alt`,
`gcnv_cnv_coherence`, `gcnv_class_coherence`, `gcnv_p_active`) so the PON can be re-tuned
without editing code.

### ⚠️ gCNV sensitivity trade-off — READ before rebuilding the PON

Broad's WDL defaults (`p-alt 5e-4`, `coherence 10000`) are **LESS sensitive** than the
pipeline's prior ad-hoc high-sensitivity values (`1e-3` / `1000`). The evaluation issue was
UNDER-calling, so Broad-aligned gCNV will call **fewer**, not more, CNVs — but it is
validated, and the specific missing events are covered by **Delly (PASS)**, not gCNV.

To restore higher sensitivity, set in `nextflow_pon.config`:
`gcnv_p_alt = "1e-3"`, `gcnv_cnv_coherence = "10000.0"→"1000.0"`, `gcnv_class_coherence = "1000.0"`.
**Decide this before rebuilding the PON.** (Source: Broad `cnv_germline_cohort_workflow.wdl`.)

### Rebuilding the PON (required)

gCNV hyperparameters are baked into the cohort model; case mode only reads it. After
changing them you **must** re-run `main_pon.nf` to regenerate `gcnv_model/` and
`cnvkit_reference/`. CNVkit `--filter cn` removal and the Delly/VQSR/mito/alignment fixes
take effect on the **next case run** and do NOT need a PON rebuild.

---

## Phasing + compound merging (`--run_phasing`) — 2026-07

**Goal:** merge caller-split adjacent/overlapping *cis* variants (e.g. SUZ12
`c.2168_2170delAAAinsTT`) into one canonical MNV so tertiary VEP reports the correct
combined `p.` (`p.Glu723_Thr724delinsAla`), matching outside labs. Default **off**.

**NCKUH — per caller, BEFORE the ensemble merge.** `main.nf` runs `PHASE_COMBINE`
(`modules/phasing.nf`) on each raw single-sample caller VCF (DV, HC): `whatshap phase`
(single-sample `--ignore-read-groups`, per-contig scatter) → `scripts/combine_phased.py`
→ then `BCFTOOLS_ENSEMBLE`. So `ensemble.fixed` already carries phase (PS/`|`) + combined
compounds; tertiary `prepare_vcf` reads `*.ensemble.fixed.vcf.gz` unchanged.
- *Why before merge:* `bcftools merge` turns DV/HC's differing compound representations
  into multiallelic, and whatshap **skips multiallelic** → the compound never gets phased.
  Must phase+combine while still single-caller & biallelic.
- *Why no sex-aware ploidy sharding (supersedes the old approach):* phasing runs on the
  raw **pre-`+fixploidy`** VCFs (uniformly diploid) → no `PloidyError` → plain per-contig
  scatter (primary contigs phased, everything else passthrough). `+fixploidy` still runs
  in `BCFTOOLS_ENSEMBLE`.

**DRAGEN — in tertiary.** `NGSTertiaryAnalysis/modules/prepare_vcf_dragen.nf`'s
`COMBINE_DRAGEN` runs the same `combine_phased.py` using DRAGEN's **native PS** (no
whatshap), gated by `params.combine_phased` (default true). Does NOT touch NCKUH's
`prepare_vcf`.

**`combine_phased.py`** (stdlib-only; **duplicated byte-identical in the tertiary repo —
keep in sync**; md5 must match): clusters variants by reference footprint (overlap, or
*cis* gap ≤ `combine_max_gap`, default 2), local-haplotype-reconstructs each cluster into
an MNV; het-cis / hom only, non-overlapping trans left alone; ref-free trim first to drop
caller "padding". Tests: `python3 scripts/test_combine_phased.py`.

### ⚠️ Combined records must keep depth (`AD`/`DP`/`VAF`) — do NOT emit `GT:PS` only

A combined MNV **inherits the full FORMAT of an anchor** (= the cluster's widest biallelic
record, e.g. the deletion in a del+ins compound), overwriting only `GT` (the reconstructed
GT) and `PS`; `QUAL`/`FILTER` also come from the anchor. The anchor shares the cluster's
ploidy, so its `AD`(Number=R)/`PL`(Number=G) element counts already match the biallelic MNV,
so inheritance is length-correct. Per the DRAGEN bug report §4 this keeps the **original
locus** `VAF`/`AD` rather than recomputing a misleading `1.0` from a 2-element `AD`.

Reconstruction is **ploidy-aware**: **diploid** clusters rebuild two haplotypes (phased
`0|1`/`1|1`/…); **all-haploid non-mito** clusters (male non-PAR `chrX`/`chrY`) rebuild the
single copy → hemizygous `GT=1` (no `PS`), still inheriting the anchor's AD. Four cases
**do not reconstruct** — they pass the source records through untouched so their `AD`
survives (downstream `bcftools norm -m -any` splits them): (a) reconstruction yields 2 ALTs
(`1|2`, incl. native multiallelic `1/2`); (b) no biallelic anchor in the cluster; (c) mixed
ploidy (haploid + diploid in one cluster); (d) `chrM` haploid clusters (multi-copy
heteroplasmy — not safe to treat as one molecule). The stderr line reports `merged_clusters`
(with `haploid=`) and `passthrough_clusters`. (NCKUH combine runs pre-`+fixploidy` =
uniformly diploid, so the haploid path is in practice DRAGEN-only.)

> The earlier version emitted only `GT:PS` on combined records, dropping `AD`/`DP`/`VAF` to
> `.` — this silently killed depth on 145k+ phased records (DRAGEN tertiary "AD 消失" bug,
> VAL-58 `chr17:80260571` / confirmed on VAL-10). Fixed by anchor inheritance + passthrough.

### ⚠️ Ensemble `FORMAT/AD` header reconcile (required, or tertiary dies)

`whatshap` re-declares **DeepVariant's** `AD` header as `Number=.` in the phased VCF (HC
stays `Number=R`). That mismatch makes `bcftools norm`/`merge` mishandle `AD` → malformed
`AD` that crashes the merge (`cannot merge`) or tertiary (`wrong number of fields in
FMT/AD`). So `BCFTOOLS_ENSEMBLE`, per caller before merging: force `AD`→`Number=R` and
`PL`→`Number=G` (sed the `##FORMAT` line + `bcftools reheader -h`) → `bcftools norm -m
-any` to biallelic → `bcftools merge --merge all`. A **pre-publish preflight**
(`bcftools norm -m -any … -Ou -o /dev/null`) makes secondary **fail loud** if any
`Number=A/R/G` field is still malformed (protects against, e.g., `VAF` too). Never use
`norm --force` (drops the tag → silently loses AD/VAF).

**Side effects to know:** `ensemble.fixed` is now **biallelic-split** at former
multiallelic sites (benign — tertiary's `norm -m -any` becomes a no-op). Combining
**lowers** the variant count (compound multi-records → one MNV); validate by specific
sites (SUZ12) + `combine_phased.py` stderr, **not** by total count.

**Validation status (2026-07):** secondary confirmed (VAL55 SUZ12 → `GAAA>GTT`; NA12878
`chr1:111241360` AD well-formed; preflight passes). Combined-record depth-preservation fix
confirmed by unit+integration tests (`test_combine_phased.py`, 13 cases) and a CLI smoke run
(SUZ12 compound keeps `AD=30,12`; reporter's `1/2 AD=0,28,20` passes through intact). Pending:
tertiary NCKUH end-to-end `-resume` (`ADD_CALLERS_TAG`); a real DRAGEN sample re-run to confirm
`AD_DRAGEN` now populates; broader multi-sample validation before clinical use.

---

## Conventions & gotchas

- **Commercial licensing is a hard constraint.** Every default-path tool must be free for
  commercial use. Non-commercial tools (Manta, ExpansionHunter, AutoMap) stay behind opt-in
  flags, default OFF. See README license table.
- **`--optical-duplicate-pixel-distance 2500`** in `alignment.nf` targets NovaSeq/NextSeq
  patterned flowcells (use 100 for HiSeq2000).
- **Per-base quality is already binned by the instrument** (NovaSeq X ~{2,9,24,40}, NextSeq
  2000 ~{2,12,26,34}). Losslessly re-encoding those bins is safe; the pipeline thresholds
  (fastp Q15/Q20, min-base-quality 10) all sit in the Q9–Q24 gap. Don't merge the top bins.
- **DeepVariant reads the BAM's original base qualities** (no recal file); HaplotypeCaller
  applies the recal table on-the-fly. Keep that in mind for any base-quality change.
- CRLF: keep files LF only (`sed -i 's/\r//'` after editing on Windows).

## Verify a change

No GPU/containers here means static review only; validate on a real sample. Syntax check
with `-preview`. For CNV, compare `.cns` (pre-filter) vs `.call.cns` (post-filter) segment
counts to tell whether a miss is a filter issue (fixable now) or a coverage/PON issue
(needs PON rebuild).
