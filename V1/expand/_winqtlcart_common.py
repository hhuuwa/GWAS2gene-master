#!/usr/bin/env python3
import re
from typing import Dict, List, Tuple


def _safe_float(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        if text.startswith("1.#"):
            return 1.0
        if text.startswith("-1.#"):
            return -1.0
        return 0.0


def parse_trait_and_area(trait_arg: str) -> Tuple[str, str]:
    m = re.match(r"^(SH|HN|HZ)_(\S+)$", trait_arg)
    if not m:
        raise ValueError(f"Invalid trait argument: {trait_arg}")
    return m.group(2), m.group(1)


def parse_qrt_lod(file_path: str, trait: str) -> Dict[int, Dict[int, float]]:
    lod: Dict[int, Dict[int, float]] = {}
    hit_trait = False
    hit_count = 0
    trait_id = None
    trait_pat = re.compile(r"-trait\s+(\d+)\s+Analyzed trait \[(\S+)\]")
    row_pat = re.compile(r"^(\s+)?(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+")
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = trait_pat.search(line)
            if m:
                cur_id = int(m.group(1))
                cur_name = m.group(2)
                if cur_name == trait:
                    hit_trait = True
                    hit_count += 1
                    trait_id = cur_id
                if trait_id is not None and hit_count > 1:
                    break
                continue
            if not hit_trait:
                continue
            rm = row_pat.match(line)
            if not rm:
                continue
            chr_id = int(rm.group(2))
            marker = int(rm.group(3))
            score = _safe_float(rm.group(5)) * 2.5 / 11.5
            lod.setdefault(chr_id, {})
            if marker not in lod[chr_id] or score > lod[chr_id][marker]:
                lod[chr_id][marker] = score
    return lod


def parse_bin_map(file_path: str) -> Tuple[Dict[int, Dict[int, float]], Dict[int, Dict[int, float]], Dict[int, int]]:
    position: Dict[int, Dict[int, float]] = {}
    position_star: Dict[int, Dict[int, float]] = {}
    position_end: Dict[int, Dict[int, float]] = {}
    marker_num: Dict[int, int] = {}

    chr_now = 0
    num = 0
    pos = 0.0
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = re.match(r"chr(\d+)\t(\S+)", line)
            if not m:
                continue
            chr_id = int(m.group(1))
            cur_pos = float(m.group(2))
            if chr_id == chr_now:
                num += 1
                position.setdefault(chr_id, {})[num] = cur_pos
                position_star.setdefault(chr_id, {})[num] = pos
                position_end.setdefault(chr_id, {})[num] = cur_pos * 2 - pos
                pos = position_end[chr_id][num]
            else:
                marker_num[chr_now] = num
                chr_now += 1
                num = 1
                pos = 0.0
                position.setdefault(chr_now, {})[num] = cur_pos
                position_star.setdefault(chr_now, {})[num] = pos
                position_end.setdefault(chr_now, {})[num] = cur_pos * 2
                pos = position_end[chr_now][num]
    marker_num[chr_now] = num
    return position_star, position_end, marker_num


def fill_forward_lod(lod_raw: Dict[int, Dict[int, float]], marker_num: Dict[int, int]) -> Dict[int, Dict[int, float]]:
    out: Dict[int, Dict[int, float]] = {}
    for chr_id in range(1, 13):
        out[chr_id] = {}
        last = 0.0
        for i in range(1, marker_num.get(chr_id, 0) + 1):
            if i in lod_raw.get(chr_id, {}):
                last = lod_raw[chr_id][i]
            out[chr_id][i] = last
    return out


def build_lod_windows(
    lod_ff: Dict[int, Dict[int, float]],
    position_star: Dict[int, Dict[int, float]],
    position_end: Dict[int, Dict[int, float]],
    marker_num: Dict[int, int],
) -> Dict[int, Dict[int, float]]:
    lod_win: Dict[int, Dict[int, float]] = {}
    for chr_id in range(1, 13):
        lod_win[chr_id] = {}
        for i in range(1, marker_num.get(chr_id, 0) + 1):
            value = lod_ff[chr_id][i]
            s = int(position_star[chr_id][i])
            e = int(position_end[chr_id][i])
            for j in range(s, e + 1):
                lod_win[chr_id][j] = value
    return lod_win
