#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def run_help(script: Path) -> str:
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=REPO_ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def test_main_help_includes_workflow_examples_and_outputs() -> None:
    help_text = run_help(REPO_ROOT / "V1" / "bin" / "riceg2g.py")

    assert "Examples:" in help_text
    assert "TGW 200kb-window fallback" in help_text
    assert "Outputs:" in help_text
    assert "RAP_Step2_<trait>.annotated.tsv" in help_text
    assert "--farmcpu-p-threshold" in help_text
    assert "--candidate-mode {auto,winqtl,window,both}" in help_text


def test_farmcpu_help_describes_significant_filtering() -> None:
    help_text = run_help(REPO_ROOT / "V1" / "expand" / "Select_FarmCPU_Peak.py")

    assert "Extract significant FarmCPU sites" in help_text
    assert "data_TGW/all_TGW_quality.FarmCPU.csv" in help_text
    assert "<trait>.FarmCPU.significant.tsv" in help_text
    assert "--p-threshold" in help_text


def test_annotation_help_describes_top10_scoring() -> None:
    help_text = run_help(REPO_ROOT / "V1" / "bin" / "annotate_recommendations.py")

    assert "Top10 candidate genes" in help_text
    assert "Recommendation_score" in help_text
    assert "all.locus_brief_info.7.0.with_keyword.tsv" in help_text
    assert "Examples:" in help_text


if __name__ == "__main__":
    test_main_help_includes_workflow_examples_and_outputs()
    test_farmcpu_help_describes_significant_filtering()
    test_annotation_help_describes_top10_scoring()
