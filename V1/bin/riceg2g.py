#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys


def run(cmd, cwd):
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="RiceG2G Python Toolkit V1 (bioinformatics-friendly wrapper)."
    )
    parser.add_argument("--trait", required=True, help="Trait id, e.g. HZ_Awn_length")
    parser.add_argument("--keyword", required=True, help="Trait keyword for Step2, e.g. Awn_length")
    parser.add_argument("--tissue", required=True, help="Main tissue, e.g. young_panicle")
    parser.add_argument(
        "--step2-trait",
        default=None,
        help="Step2 file prefix trait, default same as --trait",
    )
    parser.add_argument(
        "--skip-mlm",
        action="store_true",
        help="Skip MLM peak generation.",
    )
    parser.add_argument(
        "--skip-winqtl",
        action="store_true",
        help="Skip WinQTLcart peak + LOD generation.",
    )
    parser.add_argument(
        "--demo-viz",
        action="store_true",
        help="Also build demo visualization assets after pipeline run.",
    )
    args = parser.parse_args()

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    v1_dir = os.path.join(root, "V1")
    py = sys.executable
    step2_trait = args.step2_trait or args.trait

    if not args.skip_mlm:
        run([py, os.path.join(v1_dir, "expand", "Select_MLM_Peak.py"), args.trait], cwd=root)

    if not args.skip_winqtl:
        run([py, os.path.join(v1_dir, "expand", "Select_Winqtlcart_Peak.py"), args.trait], cwd=root)
        run([py, os.path.join(v1_dir, "expand", "WinQTLcart_lod_geno.py"), args.trait], cwd=root)
        run([py, os.path.join(v1_dir, "expand", "WinQTLcart_lod_geno_rap.py"), args.trait], cwd=root)

    run(
        [py, os.path.join(v1_dir, "code", "select_FarmCPUpeak_info_RAP.py"), args.trait, args.tissue],
        cwd=root,
    )
    run(
        [py, os.path.join(v1_dir, "code", "RAP_Stpe2_select_FarmCPUpeak_info.py"), step2_trait, args.keyword],
        cwd=root,
    )

    print("[DONE] Output:", f"RAP_Step2_{step2_trait}.FarmCPUpeak_info")
    if args.demo_viz:
        run(
            [
                py,
                os.path.join(v1_dir, "bin", "make_demo_assets.py"),
                "--input",
                f"RAP_Step2_{step2_trait}.FarmCPUpeak_info",
                "--output-dir",
                os.path.join("V1", "demo", "results"),
            ],
            cwd=root,
        )


if __name__ == "__main__":
    main()
