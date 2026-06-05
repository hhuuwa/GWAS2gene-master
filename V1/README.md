# RiceG2G Python Toolkit V1

This is the packaged **V1** workflow for the RiceG2G project.

V1 is organized as a lighter, easier-to-run layout for demonstration, testing, and code sharing. It focuses on the core candidate-prioritization flow and includes small example inputs plus generated demo outputs.

## What V1 includes

- Unified CLI wrapper: `V1/bin/riceg2g.py`
- Annotation and top10 recommender: `V1/bin/annotate_recommendations.py`
- Visualization builder: `V1/bin/make_demo_assets.py`
- Windows easy launcher: `V1/bin/run_v1.ps1`
- One-click demo runner: `V1/demo/run_demo.ps1`
- Demo and visualization runner: `V1/demo/run_demo_with_viz.ps1`
- Frozen V1 scripts:
  - `V1/expand/Select_MLM_Peak.py`
  - `V1/expand/Select_Winqtlcart_Peak.py`
  - `V1/expand/WinQTLcart_lod_geno.py`
  - `V1/expand/WinQTLcart_lod_geno_rap.py`
  - `V1/code/select_FarmCPUpeak_info_RAP.py`
  - `V1/code/RAP_Stpe2_select_FarmCPUpeak_info.py`

## Directory guide

- `V1/bin`: entrypoints and Windows launchers
- `V1/code`: Step2 candidate-prioritization scripts
- `V1/expand`: helper scripts for peak and QTL processing
- `V1/demo`: runnable examples and bundled result snapshots
- `V1/docs`: supplementary usage notes
- `V1/docs/USAGE.md`: command help, workflows, examples, and troubleshooting
- `V1/source`: keyword annotation source used in summaries
- `V1/test_data`: lightweight files for quick validation

## Pipeline note

V1 is configured to **skip `SelectFarmCPUPeakByCor.pl`**, so the packaged path does not require `Parent_unmiss.tped`.
For newer rMVP outputs, V1 now prefers `data_flowering/<trait>.FarmCPU_signals.csv` as the FarmCPU peak source. If that file is absent, it falls back to full `*.FarmCPU.csv` input, including files placed under `data_*` directories such as `data_TGW`, and first extracts significant sites with the default threshold `p <= 1e-5`.

## Quick demo (Windows PowerShell)

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo.ps1
```

Expected output file:

- `RAP_Step2_HZ_Awn_length.FarmCPUpeak_info`
- `RAP_Step2_HZ_Awn_length.annotated.tsv`
- `HZ_Awn_length_top10_candidate_genes.tsv`

## Demo with visualization

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo_with_viz.ps1
```

Generated assets:

- `V1/demo/results/top_candidates.tsv`
- `V1/demo/results/qtl_summary.tsv`
- `V1/demo/results/top_candidates.svg`
- `V1/demo/results/qtl_max_score.svg`
- `V1/demo/results/report.html`
- `V1/demo/results/score_components_heatmap.svg`

The demo summary tables include:

- `Symbol`
- `Keyword`

These annotations are merged from:

- `V1/source/rice_new_keyword.txt`

## Inputs expected by V1

V1 works from prepared locus-level input files derived from GWAS and related evidence sources. The packaged demo already includes a small example, so you can inspect the workflow without reconstructing the full research environment first.

The most relevant example files in this repository are:

- `V1/test_data/RAP_HZ_Awn_length_demo.FarmCPUpeak_info`
- `V1/test_data/RAP_Step2_HZ_Awn_length_demo.FarmCPUpeak_info`
- `V1/source/rice_new_keyword.txt`

## General run (Windows PowerShell)

For detailed command help:

```powershell
Get-Help .\V1\bin\run_v1.ps1 -Detailed
.\.tools\python\python\python.exe .\V1\bin\riceg2g.py --help
.\.tools\python\python\python.exe .\V1\expand\Select_FarmCPU_Peak.py --help
.\.tools\python\python\python.exe .\V1\bin\annotate_recommendations.py --help
```

See `V1/docs/USAGE.md` for the full usage guide.

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait HZ_Awn_length -Keyword Awn_length -Tissue young_panicle -Step2Trait HZ_Awn_length
```

HNHZ rMVP signals run with both WinQTLcart-backed and 200kb-window candidate sets:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait HNHZ -Keyword Heading_date -Tissue leaves_flowering_stage -WinQtlTrait HZ_Heading_date -CandidateMode both -ForceFarmCpu -ForceRap
```

This also writes:

- `RAP_Step2_HNHZ_winqtl.annotated.tsv`
- `RAP_Step2_HNHZ_window200kb.annotated.tsv`
- `HNHZ_top10_candidate_genes.tsv`

TGW full FarmCPU CSV filtering-only example:

```powershell
.\.tools\python\python\python.exe .\V1\expand\Select_FarmCPU_Peak.py all_TGW_quality --p-threshold 1e-5
```

The FarmCPU filtering step writes `<trait>.FarmCPU.significant.tsv` before selecting `<trait>.FarmCPU.peak`.

## Direct Python run

```bash
python V1/bin/riceg2g.py --trait HZ_Awn_length --keyword Awn_length --tissue young_panicle --step2-trait HZ_Awn_length
```

HNHZ direct run with both candidate modes:

```bash
python V1/bin/riceg2g.py --trait HNHZ --keyword Heading_date --tissue leaves_flowering_stage --winqtl-trait HZ_Heading_date --candidate-mode both --force-farmcpu --force-rap
```

Annotation and top10 recommendation run automatically by default. Add `--skip-annotation` if you only want the raw Step2 outputs.

## Lightweight test dataset

See `V1/test_data/README.md` for a small subset input and an expected Step2 output file for quick validation.

## Workflow notes

- Demo result files under `V1/demo/results` are example outputs for checking report format and pipeline shape.
- `V1/test_data` is intended for quick inspection and smoke testing without loading full-size files.
- The packaged V1 workflow is designed for easier reuse than the original full research workspace.
