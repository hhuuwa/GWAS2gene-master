# RiceG2G

RiceG2G is a rice GWAS-to-gene prioritization workflow for ranking likely causal genes and candidate variants at trait-associated loci.

This repository currently publishes the packaged **V1** workflow from the project. The upload is focused on runnable code, lightweight example inputs, and demo outputs rather than the full original research workspace.

## Overview

RiceG2G integrates multiple evidence types around a trait-associated locus, including:

- GWAS association signals
- gene annotation
- expression patterns
- homolog functional evidence
- variant effect interpretation
- linkage mapping support across families

The goal is to produce a ranked candidate list that helps narrow a GWAS signal toward likely causal genes and variants.

## What is included in this repository

This GitHub repository is intentionally lightweight and centered on **V1**.

Included:

- packaged V1 command-line entrypoints
- helper scripts for V1 preprocessing and scoring
- demo runners for Windows PowerShell
- lightweight test data
- example output tables, plots, and HTML report

Not included:

- the full large raw datasets from the original project workspace
- bulky local tool caches and installers
- the complete historical research directory structure

## Repository structure

Key paths:

- `V1/bin/riceg2g.py`: main Python entrypoint
- `V1/bin/run_v1.ps1`: Windows launcher for a standard V1 run
- `V1/demo/run_demo.ps1`: quick demo execution
- `V1/demo/run_demo_with_viz.ps1`: demo execution plus visualization assets
- `V1/docs/INPUT_OUTPUT.md`: input and output notes
- `V1/test_data`: small example files for validation

## Quick start

Clone the repository:

```bash
git clone git@github.com:hhuuwa/GWAS2gene-master.git
cd GWAS2gene-master
```

Run the lightweight demo on Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo.ps1
```

Run the demo and generate the report assets:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\demo\run_demo_with_viz.ps1
```

Expected demo outputs include:

- `V1/demo/results/report.html`
- `V1/demo/results/top_candidates.tsv`
- `V1/demo/results/qtl_summary.tsv`
- `V1/demo/results/top_candidates.svg`
- `V1/demo/results/qtl_max_score.svg`

## Typical V1 run

For a direct PowerShell run:

```powershell
powershell -ExecutionPolicy Bypass -File .\V1\bin\run_v1.ps1 -Trait HZ_Awn_length -Keyword Awn_length -Tissue young_panicle -Step2Trait HZ_Awn_length
```

For a direct Python run:

```bash
python V1/bin/riceg2g.py --trait HZ_Awn_length --keyword Awn_length --tissue young_panicle --step2-trait HZ_Awn_length
```

## Included example data

This repository includes a lightweight example based on `HZ_Awn_length`, so the workflow can be inspected without downloading the full original datasets first.

See:

- `V1/test_data/README.md`
- `V1/test_data/RAP_HZ_Awn_length_demo.FarmCPUpeak_info`
- `V1/test_data/RAP_Step2_HZ_Awn_length_demo.FarmCPUpeak_info`

## Full data background

The original RiceG2G project depends on larger reference and variant resources that are not bundled in this GitHub snapshot.

Historical external data references from the original project:

- `RiceG2G_Basicdata1`: <https://figshare.com/articles/dataset/RiceG2G_Basicdata1/21115783>
- `RiceG2G_Basicdata2`: <https://figshare.com/articles/dataset/RiceG2G_Basicdata2/21117007>
- `RiceG2G_Basicdata3`: <https://figshare.com/articles/dataset/RiceG2G_Basicdata3/21117010>

## Notes

- This repository upload is centered on the maintainable V1 package.
- V1 is configured to skip `SelectFarmCPUPeakByCor.pl` in the packaged demo flow.
- Some original Perl-based steps were converted or wrapped for the packaged V1 workflow.

## Citation

If you use RiceG2G in research, please cite the related publication or project source from your lab records.
