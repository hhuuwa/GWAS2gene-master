#!/usr/bin/env python3
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _winqtlcart_common import (
    fill_forward_lod,
    parse_bin_map,
    parse_qrt_lod,
    parse_trait_and_area,
)


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python expand/Select_Winqtlcart_Peak.py <AREA_TRAIT>")

    trait, area = parse_trait_and_area(sys.argv[1])
    out_path = f"{sys.argv[1]}.winqtlcart.peak"

    with open(out_path, "w", encoding="utf-8", newline="") as out:
        for nam in range(1, 16):
            qrt = f"./data/WinQTLcart/{area}_Nam{nam}-C.qrt"
            map_file = f"./data/WinQTLcart/N{nam}.bin.list.map"

            lod_raw = parse_qrt_lod(qrt, trait)
            position_star, position_end, marker_num = parse_bin_map(map_file)
            lod_ff = fill_forward_lod(lod_raw, marker_num)

            snp_chr = []
            snp_win = []
            for chr_id in range(1, 13):
                for i in range(1, marker_num.get(chr_id, 0) + 1):
                    if lod_ff[chr_id][i] >= 3.5:
                        snp_chr.append(chr_id)
                        snp_win.append(i)

            bins = defaultdict(list)
            if snp_chr:
                chr_now = 0
                bin_id = 0
                pos = 0.0
                prev = 0
                for idx in range(len(snp_chr)):
                    cur_chr = snp_chr[idx]
                    cur_pos = position_star[cur_chr].get(snp_win[idx], 0.0)
                    if cur_chr == chr_now:
                        if cur_pos - pos <= 10:
                            bins[bin_id].append(idx)
                            pos = cur_pos
                            prev = idx
                        else:
                            pos = cur_pos
                            prev = idx
                            bin_id += 1
                    else:
                        bins[bin_id].append(prev)
                        chr_now = cur_chr
                        pos = cur_pos
                        prev = idx
                        bin_id += 1
                bins[bin_id].append(len(snp_chr) - 1)

                for i in range(1, bin_id + 1):
                    cur = bins[i]
                    if not cur:
                        continue
                    if len(cur) == 1:
                        key = cur[0]
                    else:
                        best = max(lod_ff[snp_chr[j]][snp_win[j]] for j in cur)
                        key = next(j for j in cur if lod_ff[snp_chr[j]][snp_win[j]] == best)
                    c = snp_chr[key]
                    w = snp_win[key]
                    pos_s = position_star[c][w] / 10
                    pos_e = position_end[c][w] / 10
                    out.write(f"NAM{nam}\t{c}\t{pos_s}\t{pos_e}\t{lod_ff[c][w]}\n")


if __name__ == "__main__":
    main()
