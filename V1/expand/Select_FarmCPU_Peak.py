#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
from textwrap import dedent

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gwas_common import DEFAULT_SIGNIFICANT_P_THRESHOLD, find_gwas_input, read_gwas_records, select_peak_records


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def write_significant_records(output_path: str, records) -> None:
    with open(output_path, "w", encoding="utf-8", newline="") as out:
        out.write("SNP\tChromosome\tPosition\tPvalue\t-log10(Pvalue)\tREF\tALT\tEffect\tSE\n")
        for record in records:
            out.write(
                f"{record.snp}\t{record.chromosome}\t{record.position}\t{record.pvalue:.12g}\t"
                f"{record.logp:.6f}\t{record.ref}\t{record.alt}\t{record.effect}\t{record.se}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=dedent(
            """\
            Extract significant FarmCPU sites and select peak markers.

            The script supports two input styles:
              1. rMVP signal files, e.g. data_flowering/<trait>.FarmCPU_signals.csv
                 These are treated as already significant.
              2. full FarmCPU CSV files, e.g. data_TGW/all_TGW_quality.FarmCPU.csv
                 These are filtered with --p-threshold before peak merging.
            """
        ),
        formatter_class=HelpFormatter,
        epilog=dedent(
            """\
            Outputs:
              <trait>.FarmCPU.significant.tsv   all retained significant sites
              <trait>.FarmCPU.peak              one best marker per merged peak window

            Examples:
              Filter TGW full FarmCPU CSV at the default p <= 1e-5:
                python V1/expand/Select_FarmCPU_Peak.py all_TGW_quality

              Use a stricter threshold:
                python V1/expand/Select_FarmCPU_Peak.py all_TGW_quality --p-threshold 1e-6

              Use an explicit input path:
                python V1/expand/Select_FarmCPU_Peak.py all_TGW_quality --input data_TGW/all_TGW_quality.FarmCPU.csv
            """
        ),
    )
    parser.add_argument("trait", help="Trait id/file prefix, e.g. HNHZ or all_TGW_quality.")
    parser.add_argument(
        "--window-bp",
        type=int,
        default=2_000_000,
        help="Merge significant sites within this distance and keep the lowest p-value as the peak marker.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=DEFAULT_SIGNIFICANT_P_THRESHOLD,
        help="P-value threshold for full FarmCPU CSV input. Use None-like behavior only by passing an rMVP *_signals.csv input.",
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Optional explicit FarmCPU CSV path when automatic discovery is not enough.",
    )
    args = parser.parse_args()

    if args.input:
        input_path = Path(args.input)
        is_signals = input_path.name.endswith("_signals.csv")
    else:
        input_path, is_signals = find_gwas_input(args.trait, "FarmCPU")
    threshold = None if is_signals else args.p_threshold
    records = read_gwas_records(input_path, args.trait, "FarmCPU", p_threshold=threshold)
    peaks = select_peak_records(records, args.window_bp)

    significant_output_path = f"{args.trait}.FarmCPU.significant.tsv"
    write_significant_records(significant_output_path, records)

    output_path = f"{args.trait}.FarmCPU.peak"
    with open(output_path, "w", encoding="utf-8", newline="") as out:
        out.write("Chromosome\tPosition\t-log10(Pvalue)\tParent_geno\t1_genotype\t2_genotype\tEffect\tSE\n")
        for peak in peaks:
            out.write(
                f"{peak.chromosome}\t{peak.position}\t{peak.logp:.2f}\t*\t"
                f"{peak.ref}\t{peak.alt}\t{peak.effect}\t{peak.se}\n"
            )
    print(f"[INFO] FarmCPU input: {input_path}")
    if threshold is not None:
        print(f"[INFO] FarmCPU significant threshold: p <= {threshold:g}")
    print(f"[INFO] FarmCPU significant sites: {len(records)} -> {significant_output_path}")
    print(f"[INFO] FarmCPU peaks: {len(peaks)} -> {output_path}")


if __name__ == "__main__":
    main()
