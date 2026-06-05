# RiceG2G V1 Usage Guide

This guide documents the user-facing V1 commands, common workflows, and expected outputs.

## Help Commands

Python entrypoints:

```powershell
.\.tools\python\python\python.exe .\V1\bin\riceg2g.py --help
.\.tools\python\python\python.exe .\V1\expand\Select_FarmCPU_Peak.py --help
.\.tools\python\python\python.exe .\V1\bin\annotate_recommendations.py --help
```

PowerShell launcher:

```powershell
Get-Help .\V1\bin\run_v1.ps1 -Detailed
Get-Help .\V1\bin\run_v1.ps1 -Examples
```

## Which Command To Use

- Use `V1/bin/run_v1.ps1` on Windows for normal runs.
- Use `V1/bin/riceg2g.py` when you want direct Python control.
- Use `V1/expand/Select_FarmCPU_Peak.py` when you only want to filter FarmCPU significant sites and select peaks.
- Use `V1/bin/annotate_recommendations.py` when Step2 outputs already exist and you only want annotated tables plus Top10 recommendations.

## FarmCPU Input Discovery

V1 checks these FarmCPU inputs:

- `data_flowering/<trait>.FarmCPU_signals.csv`
- `data/farmCPU/<trait>.FarmCPU_signals.csv`
- `data/farmCPU/<trait>.FarmCPU/<trait>.FarmCPU_signals.csv`
- `data/farmCPU/<trait>.FarmCPU.csv`
- `data/farmCPU/<trait>.FarmCPU/<trait>.FarmCPU.csv`
- `data_flowering/<trait>.FarmCPU.csv`
- `data_*/<trait>.FarmCPU.csv`

For example, `data_TGW/all_TGW_quality.FarmCPU.csv` is found automatically when `--trait all_TGW_quality` is used.

## FarmCPU Significant-Site Filtering

Full FarmCPU CSV inputs are filtered by p-value before peak merging. The default is:

```text
p <= 1e-5
```

The filter output is:

```text
<trait>.FarmCPU.significant.tsv
```

Then nearby significant sites are merged into peaks with `--farmcpu-window-bp` or `--window-bp`, default `2000000` bp:

```text
<trait>.FarmCPU.peak
```

Filtering-only example:

```powershell
.\.tools\python\python\python.exe .\V1\expand\Select_FarmCPU_Peak.py all_TGW_quality --p-threshold 1e-5
```

Use a stricter threshold:

```powershell
.\.tools\python\python\python.exe .\V1\expand\Select_FarmCPU_Peak.py all_TGW_quality --p-threshold 1e-6
```

## Candidate Modes

`--candidate-mode` controls how candidate genes are collected around FarmCPU peaks.

- `auto`: use WinQTLcart if available, otherwise use the window fallback.
- `winqtl`: require WinQTLcart peaks and LOD files.
- `window`: use a local window around each FarmCPU peak; default is 200kb.
- `both`: write parallel WinQTLcart-backed and 200kb-window outputs.

Use `window` for FarmCPU-only datasets such as TGW when no matching WinQTLcart inputs exist.

## HNHZ Example

Run both WinQTLcart-backed and 200kb-window candidate sets:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 `
  -Trait HNHZ `
  -Keyword Heading_date `
  -Tissue leaves_flowering_stage `
  -WinQtlTrait HZ_Heading_date `
  -CandidateMode both `
  -ForceFarmCpu `
  -ForceRap
```

Expected main outputs:

```text
HNHZ.FarmCPU.significant.tsv
HNHZ.FarmCPU.peak
RAP_HNHZ_winqtl.FarmCPUpeak_info
RAP_Step2_HNHZ_winqtl.FarmCPUpeak_info
RAP_Step2_HNHZ_winqtl.annotated.tsv
RAP_HNHZ_window200kb.FarmCPUpeak_info
RAP_Step2_HNHZ_window200kb.FarmCPUpeak_info
RAP_Step2_HNHZ_window200kb.annotated.tsv
HNHZ_top10_candidate_genes.tsv
```

## TGW Example

Run TGW from `data_TGW/all_TGW_quality.FarmCPU.csv` using 200kb-window fallback:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 `
  -Trait all_TGW_quality `
  -Keyword TGW `
  -Tissue developing_seeds `
  -CandidateMode window `
  -SkipMlm `
  -SkipWinQtl `
  -ForceFarmCpu `
  -ForceRap `
  -FarmCpuPThreshold 1e-5
```

Expected main outputs:

```text
all_TGW_quality.FarmCPU.significant.tsv
all_TGW_quality.FarmCPU.peak
RAP_all_TGW_quality.FarmCPUpeak_info
RAP_Step2_all_TGW_quality.FarmCPUpeak_info
RAP_Step2_all_TGW_quality.annotated.tsv
all_TGW_quality_top10_candidate_genes.tsv
```

## Annotation And Top10 Recommendation

Annotation uses:

```text
basic_data/all.locus_brief_info.7.0.with_keyword.tsv
basic_data/RAP-MSU_2021-11-11.txt
```

Annotated outputs add:

```text
Annotation_locus
Symbol
Keyword
Locus_annotation
Locus_is_TE
Locus_is_expressed
Locus_is_representative
Recommendation_score
Recommendation_reason
```

`Recommendation_score` starts from the Step2 `Score` and adds transparent bonuses for:

- support from both WinQTLcart and 200kb-window modes
- Symbol and Keyword annotation evidence
- strong SIFT, INDEL, or SV evidence
- genotype-match evidence

The reason string is written in `Recommendation_reason`.

## Common Issues

- If there is no WinQTLcart input for a trait, use `--candidate-mode window` or `-CandidateMode window`.
- If a full FarmCPU CSV produces too many or too few sites, change `--farmcpu-p-threshold` or `-FarmCpuPThreshold`.
- If existing intermediate files should be reused, omit `--force-farmcpu` and `--force-rap`.
- If existing intermediate files should be overwritten, pass `--force-farmcpu --force-rap`.
- If Step2 scores look weak for a new trait keyword, check whether the keyword is represented in `basic_data/key_word.txt` or in the Step2 fallback rules.
