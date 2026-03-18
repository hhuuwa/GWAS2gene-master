#!/usr/bin/env python3
import math
import re
import sys
from statistics import fmean


def to_num(v: str, default: float = 0.0) -> float:
    try:
        return float(v)
    except Exception:
        return default


def load_keywords(trait_keyword: str):
    key_words = []
    key_words_masked = []
    with open("./basic_data/key_word.txt", "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            m = re.match(r"(\S+):(.*);(.*)?", line)
            if not m:
                continue
            trait = m.group(1)
            key = m.group(2)
            key2 = m.group(3) or ""
            if re.search(trait, trait_keyword, flags=re.IGNORECASE):
                key_words = [k for k in key.split(",") if k]
                key_words_masked = [k for k in key2.split(",") if k]
                break
    return key_words, key_words_masked


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: python code/RAP_Stpe2_select_FarmCPUpeak_info.py <trait_short_name> <trait_keyword>"
        )

    trait_short = sys.argv[1]
    trait_keyword = sys.argv[2]

    key_words, key_words_masked = load_keywords(trait_keyword)

    tissue = {
        "root_seeding_stage": "",
        "shoot_seeding_stage": "Leaf_length,Leaf_width,Plant_height,Culm_length,Panicle_enclosure,Panicle_length",
        "leaf_sheath_seeding_stage": "",
        "pistils": "Recombination_times_per_sample",
        "anthers": "Recombination_times_per_sample",
        "Leaves_seeding_stage": "Leaf_length,Leaf_width",
        "leaves_tillering_stage": "",
        "leaves_flowering_stage": "Leaf_angle,Leaf_length,Leaf_width,Heading_data,Plant_height,Culm_length,Panicle_enclosure,Protein,Seed_length,Seed_width",
        "young_panicle": "Hull_color,Awn_length,Panicle_length,Yield,Seed_length,Seed_width",
        "panicle_filling_stage": "Panicle_length,Panicle_enclosure",
        "tiller_buds": "Tiller_number,Mutiple_panicles_per_tiller,Yield",
        "embryos": "",
        "developing_seeds": "Protein",
        "seeds_after_ageing": "",
        "seed_Germinating_stage": "48H,60H,72H",
        "shoot_apical_meristem": "Leaf_length,Leaf_width,Plant_height,Culm_length,Panicle_enclosure,Panicle_length",
        "lamina_join": "Leaf_angle",
    }

    ranseq_name = [
        "root_seeding_stage",
        "shoot_seeding_stage",
        "leaf_sheath_seeding_stage",
        "pistils",
        "anthers",
        "Leaves_seeding_stage",
        "leaves_tillering_stage",
        "leaves_flowering_stage",
        "young_panicle",
        "panicle_flowering_stage",
        "panicle_filling_stage",
        "tiller_buds",
        "embryos",
        "developing_seeds",
        "seeds_after_ageing",
        "seed_Germinating_stage",
        "shoot_apical_meristem",
        "lamina_joint",
    ]

    main_rna_cols = []
    for idx, name in enumerate(ranseq_name):
        if re.search(trait_keyword, tissue.get(name, "")):
            main_rna_cols.append(idx + 17)

    src = f"RAP_{trait_short}.FarmCPUpeak_info"
    farm_scores = []
    win_scores = []
    with open(src, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 5:
                continue
            if row[2] not in farm_scores:
                farm_scores.append(row[2])
            for seg in row[4].split(";"):
                m = re.match(r"NAM\d+\:chr\d+:\d+\-\d+\:(\S+)", seg)
                if m and m.group(1) not in win_scores:
                    win_scores.append(m.group(1))

    sort_farm = sorted(float(x) for x in farm_scores) if farm_scores else [0.0]
    sort_win = sorted(float(x) for x in win_scores) if win_scores else [0.0]

    if len(sort_farm) >= 5:
        rank_farm = 2.0
        group_farm = 5
    else:
        rank_farm = 10.0 / max(len(sort_farm), 1)
        group_farm = max(len(sort_farm), 1)

    if len(sort_win) >= 5:
        rank_win = 2.0
        group_win = 5
    else:
        rank_win = 10.0 / max(len(sort_win), 1)
        group_win = max(len(sort_win), 1)

    out_path = f"RAP_Step2_{trait_short}.FarmCPUpeak_info"
    with open(src, "r", encoding="utf-8", errors="ignore") as fin, open(
        out_path, "w", encoding="utf-8", newline=""
    ) as out:
        out.write(
            "QTL\tFarmCPU_pos\tFarmCPU_-logP\tMLM_INFO(Chr-pos-(-log10))\tGeno\tGeno_pos\tSIFT_Score\tINDEL_Score\tSV_Score\tScore_peak\tScore_TE\tScore_annotate\tScore_expression\tScore_match\tScore\n"
        )
        for i, line in enumerate(fin):
            if i == 0:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 17:
                continue

            p = to_num(row[2], 0.0)
            idx = 0
            for j, v in enumerate(sort_farm):
                if p == v:
                    idx = j
                    break
            div = max(int(len(sort_farm) / max(group_farm, 1)), 1)
            group_seq = int(idx / div)
            score1_1 = (group_seq + 1) * rank_farm if group_seq <= group_farm else 10.0

            score1_2 = 0.0
            m = re.search(r"Chr\d+\-\d+\-(\S+)", row[3] if len(row) > 3 else "")
            if m:
                v = to_num(m.group(1), 0.0)
                if v < 10:
                    score1_2 = 0.0
                elif v < 20:
                    score1_2 = 5.0
                else:
                    score1_2 = 10.0

            s_win = []
            for seg in (row[4] if len(row) > 4 else "").split(";"):
                wm = re.match(r"NAM\d+\:chr\d+:\d+\-\d+\:(\S+)", seg)
                if not wm:
                    continue
                wv = to_num(wm.group(1), 0.0)
                widx = 0
                for j, sv in enumerate(sort_win):
                    if wv == sv:
                        widx = j
                        break
                wdiv = max(int(len(sort_win) / max(group_win, 1)), 1)
                wseq = int(widx / wdiv)
                s_win.append((wseq + 1) * rank_win if wseq <= group_win else 10.0)
            score1_3 = fmean(s_win) if s_win else 0.0
            score1 = score1_1 + score1_2 + score1_3

            trans = 0 if (len(row) > 8 and row[8] in ("transposon", "retrotransposon")) else 1
            score2 = 5.0 if (len(row) > 10 and row[10] == "Important") else 0.0

            num_at = 0
            num_rice = 0
            at_info = ""
            if len(row) > 11:
                at_parts = row[11].split("=>")
                p1 = at_parts[1] if len(at_parts) > 1 else ""
                p3 = at_parts[3] if len(at_parts) > 3 else ""
                p4 = at_parts[4] if len(at_parts) > 4 else ""
                at_info = f"{p1}{p3}{p4}"
            rice_ann = row[12] if len(row) > 12 else ""
            for k in key_words:
                if re.search(k, at_info):
                    num_at += 1
                if re.search(k, rice_ann):
                    num_rice += 1
            for k in key_words_masked:
                if re.search(k, at_info):
                    num_at -= 1
                if re.search(k, rice_ann):
                    num_rice -= 1

            if num_at >= 2 and num_rice >= 2:
                score3 = 25.0
            elif num_at >= 3 and num_rice == 1:
                score3 = 23.0
            elif num_rice >= 3 and num_at == 1:
                score3 = 23.0
            elif num_at >= 3 and num_rice == 0:
                score3 = 18.0
            elif num_rice >= 3 and num_at == 0:
                score3 = 18.0
            elif num_at >= 2 and num_rice == 1:
                score3 = 22.0
            elif num_rice >= 2 and num_at == 1:
                score3 = 22.0
            elif num_at >= 2 and num_rice == 0:
                score3 = 17.0
            elif num_rice >= 2 and num_at == 0:
                score3 = 17.0
            elif num_at == 1 and num_rice >= 1:
                score3 = 20.0
            elif num_at == 1 and num_rice >= 0:
                score3 = 10.0
            elif num_at == 0 and num_rice >= 1:
                score3 = 10.0
            else:
                score3 = 0.0

            selected = [0.0]
            for col in main_rna_cols:
                if col < len(row):
                    selected.append(to_num(row[col], 0.0))
            main_col = main_rna_cols[0] if main_rna_cols else 17
            max_val = max(selected) if selected else 0.0
            for col in main_rna_cols:
                if col < len(row) and to_num(row[col], 0.0) == max_val:
                    main_col = col
                    break

            main_expr = to_num(row[main_col], 0.0) if main_col < len(row) else 0.0
            score4 = 10.0 if main_expr >= 0.05 else 0.0
            others = []
            for ci in range(17, len(row)):
                if ci == main_col:
                    continue
                others.append(to_num(row[ci], 0.0))
            others_sorted = sorted(others, reverse=True)
            num = sum(1 for v in others_sorted if (v / 10.0) > main_expr)
            if num >= 5:
                score4 = 5.0
            else:
                trim = others_sorted[2:] if len(others_sorted) > 2 else []
                average = fmean(trim) if trim else 0.0
                if main_expr > 0.05:
                    if average == 0:
                        average = 0.01
                    ratio = main_expr / average
                    if ratio >= 10:
                        score4 = 22.0
                    elif ratio >= 5:
                        score4 = 20.0
                    elif ratio >= 2:
                        score4 = 18.0
                    elif ratio >= 1:
                        score4 = 15.0
                    else:
                        score4 = 10.0

            if trans == 0:
                score = 0.0
            else:
                var = []
                if len(row) > 13 and row[13] != "*":
                    var.append(to_num(row[13], 0.0))
                if len(row) > 14 and row[14] != "*":
                    var.append(to_num(row[14], 0.0))
                if len(row) > 15 and row[15] != "*":
                    var.append(to_num(row[15], 0.0))
                max_var = max(var) if var else 0.0
                score = score1 + score2 + score3 + score4 + max_var
                if len(row) > 16 and row[16] != "*":
                    score += to_num(row[16], 0.0)

            out.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[5]}\t{row[6]}\t{row[13] if len(row) > 13 else '*'}\t{row[14] if len(row) > 14 else '*'}\t{row[15] if len(row) > 15 else '*'}\t{score1}\t{score2}\t{score3}\t{score4}\t{row[16] if len(row) > 16 else '*'}\t{score}\n"
            )


if __name__ == "__main__":
    main()
