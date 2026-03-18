#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys


def run(cmd):
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Python pipeline for RiceG2G (skips SelectFarmCPUPeakByCor.pl)."
    )
    parser.add_argument("--trait", required=True, help="e.g. HZ_Awn_length")
    parser.add_argument("--keyword", required=True, help="e.g. Awn_length")
    parser.add_argument("--tissue", required=True, help="e.g. young_panicle")
    parser.add_argument(
        "--step2-trait",
        default=None,
        help="Trait key used by Step2 input/output names (default: same as --keyword).",
    )
    parser.add_argument(
        "--skip-winqtl",
        action="store_true",
        help="Skip WinQTLcart scripts when outputs already exist.",
    )
    args = parser.parse_args()

    py = sys.executable
    trait = args.trait
    step2_trait = args.step2_trait or args.keyword

    if not os.path.exists(f"{trait}.MLM.peak"):
        run([py, "expand/Select_MLM_Peak.py", trait])

    if not args.skip_winqtl:
        if not os.path.exists(f"{trait}.winqtlcart.peak"):
            run([py, "expand/Select_Winqtlcart_Peak.py", trait])
        if not os.path.exists(f"{trait}_NAM_LOD_geno.info"):
            run([py, "expand/WinQTLcart_lod_geno.py", trait])
        if not os.path.exists(f"{trait}_NAM_LOD_geno_rap.info"):
            run([py, "expand/WinQTLcart_lod_geno_rap.py", trait])

    run([py, "code/select_FarmCPUpeak_info_RAP.py", trait, args.tissue])
    run([py, "code/RAP_Stpe2_select_FarmCPUpeak_info.py", step2_trait, args.keyword])

    print("[DONE] Output:", f"RAP_Step2_{step2_trait}.FarmCPUpeak_info")


if __name__ == "__main__":
    main()
