#!/usr/bin/env python3
import csv
import math
import os
import sys
from collections import defaultdict


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python expand/Select_MLM_Peak.py <TRAIT_NAME>")

    trait = sys.argv[1]
    input_path = f"./data/farmCPU/{trait}.MLM.csv"
    if not os.path.exists(input_path):
        alt = f"./data/farmCPU/{trait}.MLM/{trait}.MLM.csv"
        if os.path.exists(alt):
            input_path = alt
    output_path = f"{trait}.MLM.peak"

    snps = []
    chromosome = {}
    position = {}
    pvalue = {}

    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) < 8:
                continue
            pv_raw = row[7]
            if pv_raw == "NA":
                continue
            pv = float(pv_raw)
            if pv <= 1.67e-8:
                snp = row[0].replace('"', "")
                snps.append(snp)
                chromosome[snp] = int(float(row[1]))
                position[snp] = int(float(row[2]))
                pvalue[snp] = -math.log10(pv)

    bins = defaultdict(list)
    if snps:
        chr_now = 0
        bin_id = 0
        pos = 0
        prev_snp = ""
        for snp in snps:
            if chromosome[snp] == chr_now:
                if position[snp] - pos <= 1_000_000:
                    bins[bin_id].append(prev_snp)
                    pos = position[snp]
                    prev_snp = snp
                else:
                    bins[bin_id].append(prev_snp)
                    pos = position[snp]
                    prev_snp = snp
                    bin_id += 1
            else:
                bins[bin_id].append(prev_snp)
                chr_now = chromosome[snp]
                pos = position[snp]
                prev_snp = snp
                bin_id += 1

        if chromosome[snps[-1]] == chr_now:
            bins[bin_id].append(snps[-1])
        else:
            bins[bin_id].append(snps[-1])

    with open(output_path, "w", encoding="utf-8", newline="") as out:
        for i in range(1, max(bins.keys(), default=0) + 1):
            items = bins[i]
            if not items:
                continue
            if len(items) == 1:
                key = items[0]
            else:
                best = max(pvalue[j] for j in items)
                key = next(j for j in items if pvalue[j] == best)
            out.write(f"{key}\t{chromosome[key]}\t{position[key]}\t{pvalue[key]}\n")


if __name__ == "__main__":
    main()
