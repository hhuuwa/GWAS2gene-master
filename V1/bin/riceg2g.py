#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import subprocess
import sys
from textwrap import dedent


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def run(cmd, cwd):
    print("[RUN]", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def has_mlm_input(root: str, trait: str) -> bool:
    return any(
        Path(root, path).exists()
        for path in [
            f"data/farmCPU/{trait}.MLM_signals.csv",
            f"data/farmCPU/{trait}.MLM/{trait}.MLM_signals.csv",
            f"data/farmCPU/{trait}.MLM.csv",
            f"data/farmCPU/{trait}.MLM/{trait}.MLM.csv",
            f"data_flowering/{trait}.MLM_signals.csv",
            f"data_flowering/{trait}.MLM.csv",
        ]
    )


def has_winqtl_input(root: str, winqtl_trait: str) -> bool:
    area = winqtl_trait.split("_", 1)[0]
    return Path(root, "data", "WinQTLcart", f"{area}_Nam1-C.qrt").exists()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=dedent(
            """\
            RiceG2G V1 end-to-end wrapper.

            This command runs the practical V1 workflow:
              1. extract FarmCPU peaks from rMVP signals or full FarmCPU CSV input
              2. optionally generate MLM and WinQTLcart evidence
              3. collect candidate genes with WinQTLcart intervals or a local window
              4. score candidates with association, annotation, expression, SIFT, INDEL, SV, and genotype-match evidence
              5. annotate results with all.locus_brief_info.7.0.with_keyword.tsv and write Top10 candidate genes
            """
        ),
        formatter_class=HelpFormatter,
        epilog=dedent(
            """\
            Inputs:
              FarmCPU signals: data_flowering/<trait>.FarmCPU_signals.csv
              FarmCPU full CSV: data_flowering/<trait>.FarmCPU.csv, data/farmCPU/<trait>.FarmCPU.csv, or data_*/<trait>.FarmCPU.csv
              WinQTLcart: data/WinQTLcart/<AREA>_Nam1-C.qrt ... Nam15-C.qrt
              Annotation: basic_data/all.locus_brief_info.7.0.with_keyword.tsv

            Outputs:
              <trait>.FarmCPU.significant.tsv
              <trait>.FarmCPU.peak
              RAP_<trait>.FarmCPUpeak_info
              RAP_Step2_<trait>.FarmCPUpeak_info
              RAP_Step2_<trait>.annotated.tsv
              <trait>_top10_candidate_genes.tsv

            Examples:
              HNHZ with WinQTLcart and 200kb-window outputs:
                python V1/bin/riceg2g.py --trait HNHZ --keyword Heading_date --tissue leaves_flowering_stage --winqtl-trait HZ_Heading_date --candidate-mode both --force-farmcpu --force-rap

              TGW 200kb-window fallback from data_TGW/all_TGW_quality.FarmCPU.csv:
                python V1/bin/riceg2g.py --trait all_TGW_quality --keyword TGW --tissue developing_seeds --candidate-mode window --skip-mlm --skip-winqtl --force-farmcpu --force-rap --farmcpu-p-threshold 1e-5

              Reuse existing FarmCPU/Step1 files and only refresh Step2 + annotation:
                python V1/bin/riceg2g.py --trait all_TGW_quality --keyword TGW --tissue developing_seeds --candidate-mode window --skip-farmcpu --skip-mlm --skip-winqtl
            """
        ),
    )
    parser.add_argument("--trait", required=True, help="Trait id and file prefix, e.g. HNHZ or all_TGW_quality.")
    parser.add_argument("--keyword", required=True, help="Step2 trait keyword, e.g. Heading_date, Awn_length, TGW.")
    parser.add_argument("--tissue", required=True, help="Main tissue used by Step1 compatibility, e.g. young_panicle or developing_seeds.")
    parser.add_argument(
        "--step2-trait",
        default=None,
        help="Step2 input/output prefix. Usually leave empty; use only when Step2 files use a different prefix.",
    )
    parser.add_argument(
        "--skip-mlm",
        action="store_true",
        help="Skip MLM peak generation when MLM data is absent or not needed.",
    )
    parser.add_argument(
        "--skip-winqtl",
        action="store_true",
        help="Skip WinQTLcart peak + LOD generation. Use with --candidate-mode window for FarmCPU-only data.",
    )
    parser.add_argument(
        "--winqtl-trait",
        default=None,
        help="Trait key for WinQTLcart inputs/outputs, e.g. HZ_Heading_date when --trait is HNHZ.",
    )
    parser.add_argument(
        "--candidate-mode",
        choices=["auto", "winqtl", "window", "both"],
        default="auto",
        help="Candidate interval mode. auto uses WinQTLcart when available; window uses the --candidate-window-bp fallback; both writes parallel winqtl/window200kb results.",
    )
    parser.add_argument(
        "--skip-farmcpu",
        action="store_true",
        help="Skip FarmCPU peak extraction and reuse an existing <trait>.FarmCPU.peak.",
    )
    parser.add_argument(
        "--force-farmcpu",
        action="store_true",
        help="Regenerate <trait>.FarmCPU.significant.tsv and <trait>.FarmCPU.peak even when they already exist.",
    )
    parser.add_argument(
        "--force-rap",
        action="store_true",
        help="Regenerate RAP_<trait>.FarmCPUpeak_info candidate tables even when they already exist.",
    )
    parser.add_argument(
        "--farmcpu-window-bp",
        type=int,
        default=2_000_000,
        help="Window used to merge nearby significant FarmCPU sites into one peak.",
    )
    parser.add_argument(
        "--farmcpu-p-threshold",
        type=float,
        default=1e-5,
        help="P-value threshold used to extract significant sites from full FarmCPU CSV input. rMVP *_signals.csv files are already filtered and are not thresholded again.",
    )
    parser.add_argument(
        "--candidate-window-bp",
        type=int,
        default=200_000,
        help="FarmCPU-only candidate gene window around each selected peak.",
    )
    parser.add_argument(
        "--demo-viz",
        action="store_true",
        help="Also build demo visualization assets after the pipeline run.",
    )
    parser.add_argument(
        "--skip-annotation",
        action="store_true",
        help="Skip locus keyword annotation and Top10 candidate-gene recommendation.",
    )
    parser.add_argument(
        "--annotation-output-prefix",
        default=None,
        help="Prefix for <prefix>_top10_candidate_genes.tsv.",
    )
    args = parser.parse_args()

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    v1_dir = os.path.join(root, "V1")
    py = sys.executable
    step2_trait = args.step2_trait or args.trait
    winqtl_trait = args.winqtl_trait or args.trait

    farmcpu_peak_path = Path(root, f"{args.trait}.FarmCPU.peak")
    farmcpu_regenerated = False
    if not args.skip_farmcpu and (args.force_farmcpu or not farmcpu_peak_path.exists()):
        farmcpu_cmd = [
            py,
            os.path.join(v1_dir, "expand", "Select_FarmCPU_Peak.py"),
            args.trait,
            "--window-bp",
            str(args.farmcpu_window_bp),
            "--p-threshold",
            str(args.farmcpu_p_threshold),
        ]
        run(farmcpu_cmd, cwd=root)
        farmcpu_regenerated = True

    if not args.skip_mlm and has_mlm_input(root, args.trait):
        run([py, os.path.join(v1_dir, "expand", "Select_MLM_Peak.py"), args.trait], cwd=root)
    elif not args.skip_mlm:
        print(f"[INFO] No MLM input found for {args.trait}; skipping MLM peak generation.")

    if not args.skip_winqtl and has_winqtl_input(root, winqtl_trait):
        run([py, os.path.join(v1_dir, "expand", "Select_Winqtlcart_Peak.py"), winqtl_trait], cwd=root)
        run([py, os.path.join(v1_dir, "expand", "WinQTLcart_lod_geno.py"), winqtl_trait], cwd=root)
        run([py, os.path.join(v1_dir, "expand", "WinQTLcart_lod_geno_rap.py"), winqtl_trait], cwd=root)
    elif not args.skip_winqtl:
        print(f"[INFO] No WinQTLcart input found for {winqtl_trait}; 200kb window mode remains available.")

    rap_cmd = [
        py,
        os.path.join(v1_dir, "code", "select_FarmCPUpeak_info_RAP.py"),
        args.trait,
        args.tissue,
        "--candidate-mode",
        args.candidate_mode,
        "--winqtl-trait",
        winqtl_trait,
        "--candidate-window-bp",
        str(args.candidate_window_bp),
    ]
    if args.force_rap or farmcpu_regenerated:
        rap_cmd.append("--force")
    run(rap_cmd, cwd=root)
    step2_outputs = []
    if args.candidate_mode == "both":
        for suffix in ["winqtl", "window200kb"]:
            trait_key = f"{args.trait}_{suffix}"
            run(
                [py, os.path.join(v1_dir, "code", "RAP_Stpe2_select_FarmCPUpeak_info.py"), trait_key, args.keyword],
                cwd=root,
            )
            step2_outputs.append(f"RAP_Step2_{trait_key}.FarmCPUpeak_info")
        print("[DONE] Outputs:", ", ".join(step2_outputs))
    else:
        run(
            [py, os.path.join(v1_dir, "code", "RAP_Stpe2_select_FarmCPUpeak_info.py"), step2_trait, args.keyword],
            cwd=root,
        )
        step2_outputs.append(f"RAP_Step2_{step2_trait}.FarmCPUpeak_info")
        print("[DONE] Output:", ", ".join(step2_outputs))

    if not args.skip_annotation:
        annotation_prefix = args.annotation_output_prefix or args.trait
        annotation_cmd = [
            py,
            os.path.join(v1_dir, "bin", "annotate_recommendations.py"),
            "--inputs",
            *step2_outputs,
            "--keyword-table",
            os.path.join(root, "basic_data", "all.locus_brief_info.7.0.with_keyword.tsv"),
            "--rap-msu-map",
            os.path.join(root, "basic_data", "RAP-MSU_2021-11-11.txt"),
            "--output-dir",
            root,
            "--output-prefix",
            annotation_prefix,
        ]
        run(annotation_cmd, cwd=root)
        annotated_outputs = [
            f"{Path(path).name.removesuffix('.FarmCPUpeak_info')}.annotated.tsv" for path in step2_outputs
        ]
        print("[DONE] Annotated outputs:", ", ".join(annotated_outputs))
        print("[DONE] Top10 recommendations:", f"{annotation_prefix}_top10_candidate_genes.tsv")
    if args.demo_viz:
        run(
            [
                py,
                os.path.join(v1_dir, "bin", "make_demo_assets.py"),
                "--input",
                step2_outputs[0],
                "--output-dir",
                os.path.join("V1", "demo", "results"),
            ],
            cwd=root,
        )


if __name__ == "__main__":
    main()
