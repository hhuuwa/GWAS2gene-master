# RiceG2G Python Toolkit V1

This is the packaged **V1** bioinformatics tool version for your project.

## What V1 includes

- Unified CLI wrapper: `V1/bin/riceg2g.py`
- Visualization builder: `V1/bin/make_demo_assets.py`
- Windows easy launcher: `V1/bin/run_v1.ps1`
- One-click demo runner: `V1/demo/run_demo.ps1`
- Demo + visualization runner: `V1/demo/run_demo_with_viz.ps1`
- Frozen V1 scripts:
  - `V1/expand/Select_MLM_Peak.py`
  - `V1/expand/Select_Winqtlcart_Peak.py`
  - `V1/expand/WinQTLcart_lod_geno.py`
  - `V1/expand/WinQTLcart_lod_geno_rap.py`
  - `V1/code/select_FarmCPUpeak_info_RAP.py`
  - `V1/code/RAP_Stpe2_select_FarmCPUpeak_info.py`

## Pipeline note

V1 is configured to **skip `SelectFarmCPUPeakByCor.pl`** (no `Parent_unmiss.tped` required).

## Quick demo (Windows PowerShell)

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo.ps1
```

Expected output file:

- `RAP_Step2_HZ_Awn_length.FarmCPUpeak_info`

## Demo with visualization (recommended)

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo_with_viz.ps1
```

Generated report:

- `V1/demo/results/report.html`
- `V1/demo/results/top_candidates.tsv`
- `V1/demo/results/qtl_summary.tsv`
- `V1/demo/results/score_components_heatmap.svg`

The demo summary tables now include:

- `Symbol`
- `Keyword`

These annotations are merged from:

- `V1/source/rice_new_keyword.txt`

## Demo with visualization (recommended)

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo_with_viz.ps1
```

Generated assets:

- `V1/demo/results/top_candidates.tsv`
- `V1/demo/results/qtl_summary.tsv`
- `V1/demo/results/top_candidates.svg`
- `V1/demo/results/qtl_max_score.svg`
- `V1/demo/results/report.html`

## Test data

See `V1/test_data` for a lightweight subset input and corresponding expected Step2 output.

## General run (Windows PowerShell)

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait HZ_Awn_length -Keyword Awn_length -Tissue young_panicle -Step2Trait HZ_Awn_length
```

## Direct Python run

```bash
python V1/bin/riceg2g.py --trait HZ_Awn_length --keyword Awn_length --tissue young_panicle --step2-trait HZ_Awn_length
```

## Lightweight test dataset

See `V1/test_data/README.md` for a small subset input and an expected Step2 output file for quick validation.

## Lightweight test data

Use `V1/test_data` to quickly inspect example input/output formats without loading the full-size files.
