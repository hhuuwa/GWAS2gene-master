#!/usr/bin/env python3
import os
import re
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _winqtlcart_common import (
    build_lod_windows,
    fill_forward_lod,
    parse_bin_map,
    parse_qrt_lod,
    parse_trait_and_area,
)


def parse_genes_rap(gff_path: str) -> List[Tuple[str, int, int, int, int]]:
    out = []
    pat = re.compile(r"chr(\d+)\tirgsp1_locus\tgene\t(\d+)\t(\d+)\t")
    id_pat = re.compile(r"ID=(\S+);Name=(\S+);")
    with open(gff_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = pat.search(line)
            if not m:
                continue
            chr_id = int(m.group(1))
            start = int(m.group(2))
            end = int(m.group(3))
            im = id_pat.search(line)
            if not im:
                continue
            gid = im.group(1)
            mid = (start + end) // 2
            out.append((gid, chr_id, start, end, mid))
    return out


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python expand/WinQTLcart_lod_geno_rap.py <AREA_TRAIT>")

    trait, area = parse_trait_and_area(sys.argv[1])
    lod_win_by_nam: Dict[int, Dict[int, Dict[int, float]]] = {}

    for nam in range(1, 16):
        qrt = f"./data/WinQTLcart/{area}_Nam{nam}-C.qrt"
        map_file = f"./data/WinQTLcart/N{nam}.bin.list.map"

        lod_raw = parse_qrt_lod(qrt, trait)
        position_star, position_end, marker_num = parse_bin_map(map_file)
        lod_ff = fill_forward_lod(lod_raw, marker_num)
        lod_win_by_nam[nam] = build_lod_windows(lod_ff, position_star, position_end, marker_num)

    genes = parse_genes_rap("./basic_data/Rice_IRGSP-1.0.gff3")
    out_path = f"{sys.argv[1]}_NAM_LOD_geno_rap.info"
    with open(out_path, "w", encoding="utf-8", newline="") as out:
        out.write(
            "Geno_id\tChr\tGeno_star\tGeno_end\tNAM1\tNAM2\tNAM3\tNAM4\tNAM5\tNAM6\tNAM7\tNAM8\tNAM9\tNAM10\tNAM11\tNAM12\tNAM13\tNAM14\tNAM15\n"
        )
        for gid, chr_id, start, end, mid in genes:
            bin_pos = int(mid / 100000)
            vals = []
            for nam in range(1, 16):
                v = lod_win_by_nam[nam].get(chr_id, {}).get(bin_pos, "NA")
                vals.append(str(v))
            out.write(f"{gid}\t{chr_id}\t{start}\t{end}\t" + "\t".join(vals) + "\n")


if __name__ == "__main__":
    main()
