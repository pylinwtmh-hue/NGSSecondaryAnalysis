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
combined `p.` (`p.Glu723_Thr724delinsAla`), matching outside labs. Default **on**
(`params.run_phasing = true` since DGX validation; `--run_phasing false` skips it — the config
comment still says "預設 OFF", which is stale).

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

### ⚠️ Records with no ALT never join a cluster (2026-09)

`process()` writes records whose sample GT has no ALT allele (`./.`, `0/0`, haploid `0`/`.` —
typically DeepVariant's rejected candidates, `FILTER=RefCall`) straight through; they are
not clustered. Before this, "overlap always merges" + "anchor = widest biallelic record"
ignored call status, so a wider rejected candidate (usually a deletion) covering a real call
(e.g. an SNV) became the anchor: the combined record carried the rejected candidate's
`QUAL`/`FILTER=RefCall`/`GQ`/`DP`/`AD`/`VAF`/`PL`, was re-anchored at its `POS` (so it no longer
lined up with HC's call and `merge` could not join them — tertiary `norm` then produced the same
variant twice, `CALLERS=DV` with the wrong depth + `CALLERS=HC`), got a fake `0|1`+`PS`, and could
bridge two unrelated calls into one MNV. Its GT has an ALT, so the ensemble's DV `GT="alt"` filter
kept it. VAL55: all 28,050 `FILTER=RefCall` records left in `ensemble.fixed` carried `COMBINED`;
23,023 variants were split into a DV row + an HC row in the tertiary report. Tests:
`test_nocall_wider_candidate_not_anchor`, `test_nocall_does_not_bridge`,
`test_haploid_nocall_passthrough`. stderr now ends with `nocall_passthrough=` (≈ RefCall count
for DV, 0 for HC). `COMBINE_PHASED` stages the script as an input, so `-resume` re-runs from
there. Check after a run: `bcftools view -H -f RefCall <id>.ensemble.fixed.vcf.gz | wc -l` → 0.
**Verified on the VAL55 re-run:** RefCall 28,050 → 0; `nocall_passthrough` DV 865,130 / HC 0;
split variants 23,023 → 820; tertiary DV+HC / DV only / HC only = 89.9% / 2.6% / 7.4%. The 820
left (5 checked) are combine outputs that keep leading bases shared by REF/ALT: the span starts at
the first component's POS, which for an indel includes its anchor base, and the result is not
minimised. DV and HC split one event differently, so the padding and POS differ and `merge`
cannot join them; trimming both sides makes all 5 identical. **Fixed:** `flush()` now runs
`trim_alleles()` on the reconstructed REF/ALT before rendering (a pure deletion keeps its anchor
base); SUZ12 is therefore written `31998951 AAA>TT`, not `31998950 GAAA>GTT`. Test:
`test_merged_output_minimised`.

### ⚠️ Male haploid-region hets are handled before `+fixploidy` (2026-09)

`+fixploidy` keeps only the **first** allele when it makes a GT haploid (`plugins/fixploidy.c`):
`0/1`, `0|1` → `0` (REF, vanishes from the report), `1|0` → `1` (hemizygous). With phasing the
outcome depended on whatshap's arbitrary orientation (VAL55: chrX 2,607 + chrY 6,049 hets
truncated to REF = the whole tertiary NONE count; hundreds more shown as hemizygous). A male het is
usually mis-mapping, but can be 47,XXY or **somatic mosaicism** — how males are affected by
X-linked dominant, male-lethal disorders (IKBKG, MECP2, CDKL5, PORCN, OFD1; PCDH19 affects mosaic
males only) — so it must be neither hidden nor phase-dependent. For **male** samples only,
`BCFTOOLS_ENSEMBLE` pipes the merged VCF through `scripts/haploid_het.awk` (staged input; reads the
same ploidy file) before `+fixploidy`: chrX ploidy-1 regions: het `0/k` → `k/k` plus
`INFO/HAPLOID_HET=<DV,HC>` (tertiary shows it as the `HAPLOID_HET` review column; AD/VAF untouched);
chrY: het → `./.` (not reported); PAR, autosomes, chrM, female/unknown: untouched. POSIX awk only
(the bcftools container is busybox-based; tested with mawk and busybox awk:
`AWK="busybox awk" python3 scripts/test_haploid_het.py`). The pipe runs in a `set -o pipefail`
subshell so a failure is not swallowed by `bash -ue`. stderr: `[haploid_het] … chrX_het_to_alt_*
chrY_het_to_missing_*`. Germline callers still miss low-level mosaics; this only rescues calls that
were made as het. Verified on VAL55: 4,326 tagged records, all on chrX; tertiary NONE 10,159,
all on chrY; split DV/HC rows 820 → 2 (left-alignment differences of repeat indels — adding
`norm -f` before the merge would remove them; not done). **chrM is deliberately untouched**
(lab decision): its hets are still truncated by `+fixploidy`; mtDNA is read from `07_mitochondria`.

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

### ⚠️ DV non-calls are dropped before the merge (2026-09)

After the DV arm's `norm -m -any`, `BCFTOOLS_ENSEMBLE` keeps only `GT="alt"` records
(two steps via `norm_dv.bcf`, not a pipe — the shell is `bash -ue` without `pipefail`).
DeepVariant's VCF keeps candidates it **rejected** (`FILTER=RefCall`, GT `./.`/`0/0`). If
they reach `merge --merge all` and HC has a *different* allele at the same POS, the two are
fused into one multiallelic record (ALT1 = DV's rejected candidate, ALT2 = HC's real call)
whose FILTER becomes DV's `RefCall`; tertiary `norm` then splits the rejected allele back
out as a record **neither caller called**. Real case — SUZ12 `chr17:31998950`: DV rejected
all three delinsTT components; HC's combine produced the correct `GAAA>GTT 0|1`; the merge
produced `GAAA GAA,GTT`, and the report gained a bogus `c.2170del` while the true delinsTT
carried DV's meaningless `AD 10,0`. Filter *after* `norm` so a DV `0/2` splits into `0/0`
(dropped) + `0/1` (kept). Measured `GT="alt"`: keeps `0/1 0|1 1/1 1|1 1/0`; drops `./. 0/0
0|0 ./0` and half-missing `./1 1/.` — matching tertiary `add_callers_tag.is_called()`. DV's
RefCalls stay in the published raw `<id>.deepvariant.vcf.gz` (auditable; CNVkit reads that
file, not the ensemble). HC needs no equivalent (VCF mode emits ALT sites only). Tertiary
separately fixed `determine_callers()` (no-call → `NONE`, was `HC`) and `get_ad()` (keeps
missing as `.`, was `0`) — the ensemble change and those fixes are independent safeguards.

**Side effects to know:** `ensemble.fixed` is now **biallelic-split** at former
multiallelic sites (benign — tertiary's `norm -m -any` becomes a no-op), except where the
two callers carry different alleles at one POS (still fused by `--merge all`, still split in
tertiary). It also **no longer contains DV RefCall records**, so its record count drops. Combining
**lowers** the variant count (compound multi-records → one MNV); validate by specific
sites (SUZ12) + `combine_phased.py` stderr, **not** by total count.

**Validation status (2026-07):** secondary confirmed (VAL55 SUZ12 → `GAAA>GTT`; NA12878
`chr1:111241360` AD well-formed; preflight passes). Combined-record depth-preservation fix
confirmed by unit+integration tests (`test_combine_phased.py`, now 19 cases) and a CLI smoke run
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
