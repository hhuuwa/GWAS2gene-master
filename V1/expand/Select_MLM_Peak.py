#!/usr/bin/env python3
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gwas_common import DEFAULT_P_THRESHOLD, find_gwas_input, read_gwas_records, select_peak_records


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python expand/Select_MLM_Peak.py <TRAIT_NAME>")

    trait = sys.argv[1]
    input_path, _ = find_gwas_input(trait, "MLM")
    records = read_gwas_records(input_path, trait, "MLM", p_threshold=DEFAULT_P_THRESHOLD)
    peaks = select_peak_records(records, window_bp=1_000_000)

    output_path = f"{trait}.MLM.peak"
    with open(output_path, "w", encoding="utf-8", newline="") as out:
        for peak in peaks:
            out.write(f"{peak.snp}\t{peak.chromosome}\t{peak.position}\t{peak.logp}\n")
    print(f"[INFO] MLM input: {input_path}")
    print(f"[INFO] MLM peaks: {len(peaks)} -> {output_path}")


if __name__ == "__main__":
    main()
