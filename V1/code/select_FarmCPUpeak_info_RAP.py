#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


RNA_TISSUES = [
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


OUTPUT_HEADER = (
    "QTL\tFarmCPU_pos\tFarmCPU_-logP\tMLM_INFO(Chr-pos-(-log10)\tWinqtlcart_peak_info\t"
    "Geno\tGeno_pos\tLOD\ttransposon\tHouse_keeping\tImportance\tAT\tRAPDB_annotation\t"
    "SIFT_Score\tINDEL_Score\tSV_Score\tMatch_Score\t"
    + "\t".join(RNA_TISSUES)
)


@dataclass(frozen=True)
class Peak:
    chromosome: int
    position: int
    logp: str


@dataclass(frozen=True)
class Gene:
    gene_id: str
    chromosome: int
    start: int
    end: int
    expression: list[str]


@dataclass(frozen=True)
class WinQtlOverlap:
    text: str
    nam_ids: list[int]


@dataclass(frozen=True)
class Candidate:
    qtl_idx: int
    peak: Peak
    gene: Gene
    mlm_info: str
    winqtl: WinQtlOverlap


def first_existing(paths: list[str]) -> Path | None:
    for path in paths:
        p = Path(path)
        if p.is_file():
            return p
    return None


def find_resource(name: str) -> Path | None:
    return first_existing(
        [
            f"./basic_data/{name}",
            f"./basic_data/{name}/{name}",
            f"./21117007/{name}",
            f"./21115783/{name}",
            f"./{name}",
            f"./{name}/{name}",
        ]
    )


def parse_chr(value: str) -> int:
    match = re.search(r"(\d+)", value)
    if not match:
        raise ValueError(value)
    return int(match.group(1))


def gene_base(gene_id: str) -> str:
    return re.sub(r"\.\d+$", "", gene_id.strip())


def score_text(value: float | int | None) -> str:
    if value is None:
        return "*"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:g}"


def parse_score(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_farmcpu_peaks(trait: str) -> list[Peak]:
    path = Path(f"{trait}.FarmCPU.peak")
    if not path.exists():
        raise FileNotFoundError(
            f"Cannot find {path}. Run V1/expand/Select_FarmCPU_Peak.py first."
        )

    peaks: list[Peak] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 3:
                continue
            try:
                peaks.append(Peak(parse_chr(row[0]), int(float(row[1])), row[2]))
            except ValueError:
                continue
    return peaks


def normalize_expression(values: list[str]) -> list[str]:
    expr = values[: len(RNA_TISSUES)]
    if len(expr) < len(RNA_TISSUES):
        expr.extend(["*"] * (len(RNA_TISSUES) - len(expr)))
    return expr


def load_expression_genes() -> dict[int, list[Gene]]:
    path = Path("./basic_data/combineRAPvsMSU_RNAseq.txt")
    genes_by_chr: dict[int, list[Gene]] = {}
    if not path.exists():
        return genes_by_chr

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 4:
                continue
            try:
                gene = Gene(
                    gene_id=row[0],
                    chromosome=parse_chr(row[1]),
                    start=int(float(row[2])),
                    end=int(float(row[3])),
                    expression=normalize_expression(row[4:]),
                )
            except ValueError:
                continue
            genes_by_chr.setdefault(gene.chromosome, []).append(gene)

    for genes in genes_by_chr.values():
        genes.sort(key=lambda g: (g.start, g.end, g.gene_id))
    return genes_by_chr


def parse_gff_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def load_gff_annotations() -> dict[str, str]:
    annotations: dict[str, str] = {}
    for path in [
        Path("./basic_data/Rice_IRGSP-1.0.gff3"),
        Path("./basic_data/Rice_MSUv7.gff3"),
        Path("./basic_data/Rice_MSUv7/Rice_MSUv7.gff3"),
    ]:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                row = line.rstrip("\n").split("\t")
                if len(row) < 9 or row[2] != "gene":
                    continue
                attrs = parse_gff_attrs(row[8])
                gene_id = attrs.get("ID")
                if gene_id:
                    annotations[gene_id] = attrs.get("Note", "*") or "*"
    return annotations


def load_locus_info() -> tuple[dict[str, str], dict[str, str]]:
    transposons: dict[str, str] = {}
    annotations: dict[str, str] = {}
    path = Path("./basic_data/all.locus_brief_info.7.0")
    if not path.exists():
        return transposons, annotations

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            row = line.rstrip("\n").split("\t")
            if i == 0 or len(row) < 10:
                continue
            locus = row[1]
            if row[6] == "Y":
                transposons[locus] = f"{parse_chr(row[0])}transposon"
            annotations[locus] = row[9] or "*"
    return transposons, annotations


def load_gene_set(path: str) -> set[str]:
    p = Path(path)
    if not p.exists():
        return set()
    genes: set[str] = set()
    with p.open("r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            gene = line.strip().split("\t")[0]
            if gene:
                genes.add(gene)
                genes.add(gene.split(".")[0])
    return genes


def load_arabidopsis_info() -> dict[str, str]:
    at: dict[str, str] = {}
    paths = [
        "./basic_data/Basic_Information_Arabidopsis_to_rice.txt",
        "./basic_data/Basic_Information_Arabidopsis_to_rice_rap.txt",
        "./21115783/Basic_Information_Arabidopsis_to_rice.txt",
        "./21115783/Basic_Information_Arabidopsis_to_rice_rap.txt",
    ]
    for path in paths:
        p = Path(path)
        if not p.exists():
            continue
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if not row:
                    continue
                padded = row + [""] * 7
                fields = [padded[1], padded[2], padded[3], padded[5], padded[6]]
                at[padded[0]] = "=>".join(field if field else "*" for field in fields)
    return at


def load_mlm_peaks(trait: str) -> list[Peak]:
    path = Path(f"{trait}.MLM.peak")
    if not path.exists():
        return []
    peaks: list[Peak] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            row = line.rstrip("\n").split()
            if len(row) < 4:
                continue
            try:
                peaks.append(Peak(parse_chr(row[1]), int(float(row[2])), row[3]))
            except ValueError:
                continue
    return peaks


def nearby_mlm_info(mlm_peaks: list[Peak], peak: Peak) -> str:
    parts = [
        f"Chr{m.chromosome}-{m.position}-{m.logp}"
        for m in mlm_peaks
        if m.chromosome == peak.chromosome and abs(m.position - peak.position) <= 300_000
    ]
    return "".join(parts) if parts else "*"


@dataclass(frozen=True)
class WinPeak:
    nam: str
    nam_id: int
    chromosome: int
    start: int
    end: int
    lod: float


def load_winqtl_peaks(winqtl_trait: str) -> list[WinPeak]:
    path = Path(f"{winqtl_trait}.winqtlcart.peak")
    if not path.exists():
        raise FileNotFoundError(
            f"Cannot find {path}. Run V1/expand/Select_Winqtlcart_Peak.py {winqtl_trait} first."
        )
    peaks: list[WinPeak] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            row = line.rstrip("\n").split()
            if len(row) < 5:
                continue
            nam_match = re.search(r"NAM(\d+)", row[0])
            if not nam_match:
                continue
            try:
                peaks.append(
                    WinPeak(
                        nam=row[0],
                        nam_id=int(nam_match.group(1)),
                        chromosome=parse_chr(row[1]),
                        start=int(float(row[2]) * 1_000_000),
                        end=int(float(row[3]) * 1_000_000),
                        lod=float(row[4]),
                    )
                )
            except ValueError:
                continue
    return peaks


def overlapping_winqtl_peaks(win_peaks: list[WinPeak], peak: Peak) -> list[WinPeak]:
    return [
        w
        for w in win_peaks
        if w.chromosome == peak.chromosome
        and (abs(w.start - peak.position) <= 500_000 or abs(w.end - peak.position) <= 500_000)
    ]


def winqtl_overlap_text(win_peaks: list[WinPeak]) -> str:
    if not win_peaks:
        return "*"
    return "".join(
        f"{w.nam}:chr{w.chromosome}:{w.start}-{w.end}:{w.lod:g};" for w in win_peaks
    )


def candidate_region_from_winqtl(win_peaks: list[WinPeak], chr_length: int) -> tuple[int, int]:
    if len(win_peaks) == 1:
        win = win_peaks[0]
        return max(0, win.start - 200_000), min(chr_length, win.end + 200_000)

    min_lod = min(w.lod for w in win_peaks)
    max_lod = max(w.lod for w in win_peaks)
    if min_lod >= 15:
        selected = win_peaks
    else:
        selected = [w for w in win_peaks if max_lod and (w.lod / max_lod) >= 0.5]
    if not selected:
        selected = win_peaks

    starts = [w.start for w in selected]
    ends = [w.end for w in selected]
    start = duplicated_or_extreme(starts, min)
    end = duplicated_or_extreme(ends, max)
    return max(0, start - 200_000), min(chr_length, end + 200_000)


def duplicated_or_extreme(values: list[int], fallback):
    counts: dict[int, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    duplicates = [value for value, count in counts.items() if count > 1]
    if duplicates:
        return min(duplicates) if fallback is min else max(duplicates)
    return fallback(values)


def load_chr_lengths() -> dict[int, int]:
    return {
        1: 43270923,
        2: 35937250,
        3: 36413819,
        4: 35502694,
        5: 29958434,
        6: 31248787,
        7: 29697621,
        8: 28443022,
        9: 23012720,
        10: 23207287,
        11: 29021106,
        12: 27531856,
    }


def load_lod_rows(winqtl_trait: str) -> dict[str, tuple[str, list[str]]]:
    out: dict[str, tuple[str, list[str]]] = {}
    for path in [Path(f"{winqtl_trait}_NAM_LOD_geno.info"), Path(f"{winqtl_trait}_NAM_LOD_geno_rap.info")]:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                row = line.rstrip("\n").split()
                if len(row) < 19:
                    continue
                gene = gene_base(row[0])
                pos = f"{row[1]}:{row[2]}-{row[3]}"
                out[gene] = (pos, row[4:19])
    return out


def nearby_genes(genes: list[Gene], peak: Peak, window_bp: int) -> list[Gene]:
    start = max(0, peak.position - window_bp)
    end = peak.position + window_bp
    return [gene for gene in genes if gene.end >= start and gene.start <= end]


def genes_in_region(genes: list[Gene], start: int, end: int) -> list[Gene]:
    return [gene for gene in genes if gene.end >= start and gene.start <= end]


def load_sift_by_position() -> tuple[dict[tuple[int, int], float], set[tuple[int, int]]]:
    path = find_resource("SIFTtolerantScore_MSUvsRAP.txt")
    scores: dict[tuple[int, int], float] = {}
    stops: set[tuple[int, int]] = set()
    if not path:
        return scores, stops
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            row = line.split()
            if len(row) < 8:
                continue
            try:
                key = (parse_chr(row[0]), int(float(row[1])))
                scores[key] = float(row[5])
                if "STOP" in row[7] or "START" in row[7]:
                    stops.add(key)
            except ValueError:
                continue
    return scores, stops


def sift_impact_score(raw_sift: float, is_stop: bool = False) -> int:
    value = 1 - raw_sift
    if value < 0.9:
        score = 1
    elif value < 0.95:
        score = 5
    elif value < 0.98:
        score = 15
    elif value < 0.99:
        score = 18
    else:
        score = 20
    if is_stop and value >= 1:
        score = 25
    return score


def match_score_from_sift(raw_sift: float, is_stop: bool = False) -> int:
    value = 1 - raw_sift
    if value < 0.9:
        return 1
    if value < 0.95:
        return 5
    if value < 0.98:
        return 15
    if value < 0.99:
        return 18
    if value < 1:
        return 20
    return 22 if is_stop else 20


def as_bit(value: str) -> int:
    try:
        return 1 if int(float(value)) != 0 else 0
    except ValueError:
        return 1 if value not in ("0", "0|0", "*", "") else 0


def coding_snp_mis_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 7:
        return out
    p1 = row[5]
    ref = row[2]
    for idx in range(6, len(row)):
        if idx == 12:
            out.append(0 if ref == p1 else 1)
        out.append(0 if row[idx] == p1 else 1)
    return out[:15]


def coding_snp_match_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 6:
        return out
    ref = row[2]
    for idx in range(5, len(row)):
        if idx == 12:
            out.append(0)
        out.append(0 if row[idx] == ref else 1)
    return out[:16]


def coding_indel_mis_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 8:
        return out
    base = row[6]
    for idx in range(7, len(row) - 1):
        if idx == 13:
            out.append(0 if base == "0" else 1)
        out.append(0 if row[idx] == base else 1)
    return out[:15]


def coding_indel_match_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 7:
        return out
    for idx in range(6, len(row) - 1):
        if idx == 13:
            out.append(0)
        out.append(as_bit(row[idx]))
    return out[:16]


def coding_sv_mis_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 7:
        return out
    base = row[5]
    for idx in range(6, len(row) - 1):
        if idx == 12:
            out.append(0 if base == "0|0" else 1)
        out.append(0 if row[idx] == base else 1)
    return out[:15]


def coding_sv_match_vector(row: list[str]) -> list[int]:
    out: list[int] = []
    if len(row) < 6:
        return out
    for idx in range(5, len(row) - 1):
        if idx == 12:
            out.append(0)
        out.append(0 if row[idx] == "0|0" else 1)
    return out[:16]


def noncoding_mis(row: list[str], nam_ids: list[int]) -> bool:
    for nam_id in nam_ids:
        idx = nam_id + 3
        if idx < len(row) and as_bit(row[idx]) == 1:
            return True
    return False


def noncoding_match_vector(row: list[str]) -> list[int]:
    if len(row) <= 4:
        return []
    pivot = row[11] if len(row) > 11 else "0"
    body = [as_bit(v) for v in row[4:]]
    if pivot == "0":
        return ([0] + body)[:16]
    return ([1] + [0 if v else 1 for v in body])[:16]


def variant_gene_entries(text: str) -> list[tuple[str, float]]:
    entries: list[tuple[str, float]] = []
    for part in text.split(";"):
        if not part:
            continue
        m = re.match(r"([^:]+):(\S+)", part)
        if not m:
            continue
        score = parse_score(m.group(2))
        if score is not None:
            entries.append((gene_base(m.group(1)), score))
    return entries


def has_mis(genotypes: list[int], nam_ids: list[int]) -> bool:
    for nam_id in nam_ids:
        idx = nam_id - 1
        if 0 <= idx < len(genotypes) and genotypes[idx] == 1:
            return True
    return False


def find_candidate_indices_at_pos(candidates_by_chr: dict[int, list[tuple[int, int, int]]], chr_id: int, pos: int) -> list[int]:
    hits = []
    for start, end, idx in candidates_by_chr.get(chr_id, []):
        if start <= pos <= end:
            hits.append(idx)
        elif start > pos:
            break
    return hits


def score_match(variants: list[tuple[float, list[int]]], nam_ids: list[int]) -> str:
    if not variants:
        return "*"
    max_score = max(score for score, _ in variants)
    nam_set = sorted(set(nam_ids))
    if max_score >= 20:
        true_geno: set[int] = set()
        for score, vector in variants:
            if score >= 22:
                true_geno.update(idx + 1 for idx, value in enumerate(vector[:16]) if value == 1)
        base = [1 if i in true_geno else 0 for i in range(1, 17)]
        diff = [idx for idx in range(1, 16) if base[idx] != base[0]]
        return "5" if diff == nam_set else "*"
    if 10 <= max_score <= 20:
        best = [(score, vector) for score, vector in variants if score == max_score]
        vector = best[-1][1]
        true_geno = {idx + 1 for idx, value in enumerate(vector[:16]) if value == 1}
        base = [1 if i in true_geno else 0 for i in range(1, 17)]
        diff = [idx for idx in range(1, 16) if base[idx] != base[0]]
        if diff == nam_set:
            return "5" if len(best) == 1 else "abnormal"
    return "*"


def compute_variant_scores(candidates: list[Candidate]) -> tuple[dict[int, dict[str, str]], list[str]]:
    warnings: list[str] = []
    scores: dict[int, dict[str, str]] = {
        idx: {"sift": "*", "indel": "*", "sv": "*", "match": "*"} for idx in range(len(candidates))
    }
    sift_raw_values: dict[int, list[float]] = {idx: [] for idx in range(len(candidates))}
    sift_stop: dict[int, bool] = {idx: False for idx in range(len(candidates))}
    noncoding_snp_scores: dict[int, list[float]] = {idx: [] for idx in range(len(candidates))}
    indel_scores: dict[int, list[float]] = {idx: [] for idx in range(len(candidates))}
    sv_scores: dict[int, list[float]] = {idx: [] for idx in range(len(candidates))}
    match_variants: dict[int, list[tuple[float, list[int]]]] = {idx: [] for idx in range(len(candidates))}

    candidates_by_gene: dict[str, list[int]] = {}
    candidates_by_chr: dict[int, list[tuple[int, int, int]]] = {}
    for idx, candidate in enumerate(candidates):
        candidates_by_gene.setdefault(gene_base(candidate.gene.gene_id), []).append(idx)
        candidates_by_chr.setdefault(candidate.gene.chromosome, []).append((candidate.gene.start, candidate.gene.end, idx))
    for rows in candidates_by_chr.values():
        rows.sort()

    sift_by_pos, stop_by_pos = load_sift_by_position()

    for chr_id in range(1, 13):
        snp_path = find_resource(f"Chr{chr_id}_combineSNPeffect_MSUvsRAP.txt")
        if snp_path:
            with snp_path.open("r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    row = line.split()
                    if len(row) < 7 or row[0].startswith("#"):
                        continue
                    try:
                        pos = int(float(row[1]))
                    except ValueError:
                        continue
                    key = (chr_id, pos)
                    hit_indices = find_candidate_indices_at_pos(candidates_by_chr, chr_id, pos)
                    if not hit_indices:
                        continue
                    mis_vector = coding_snp_mis_vector(row)
                    match_vector = coding_snp_match_vector(row)
                    raw_sift = sift_by_pos.get(key)
                    for idx in hit_indices:
                        if has_mis(mis_vector, candidates[idx].winqtl.nam_ids):
                            if raw_sift is not None:
                                sift_raw_values[idx].append(raw_sift)
                                if key in stop_by_pos:
                                    sift_stop[idx] = True
                        if raw_sift is not None and match_vector:
                            match_variants[idx].append((match_score_from_sift(raw_sift, key in stop_by_pos), match_vector))
        else:
            warnings.append(f"Missing coding SNP resource for Chr{chr_id}")

        indel_path = find_resource(f"Chr{chr_id}_combineINDELtolerantScore_MSUvsRAP.txt")
        if indel_path:
            with indel_path.open("r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    row = line.split()
                    if len(row) < 8 or row[0].startswith("#"):
                        continue
                    gene = gene_base(row[5])
                    score = parse_score(row[-1])
                    if score is None:
                        continue
                    mis_vector = coding_indel_mis_vector(row)
                    match_vector = coding_indel_match_vector(row)
                    for idx in candidates_by_gene.get(gene, []):
                        if has_mis(mis_vector, candidates[idx].winqtl.nam_ids):
                            indel_scores[idx].append(score)
                        if match_vector:
                            match_variants[idx].append((score, match_vector))
        else:
            warnings.append(f"Missing coding INDEL resource for Chr{chr_id}")

        sv_path = find_resource(f"NAM_Chr{chr_id}combineRAPvsMSU_codingSV.txt")
        if sv_path:
            with sv_path.open("r", encoding="utf-8", errors="ignore") as f:
                reader = csv.reader(f, delimiter="\t")
                next(reader, None)
                for row in reader:
                    if len(row) < 7 or row == ["*"]:
                        continue
                    score = parse_score(row[-1])
                    if score is None:
                        continue
                    mis_vector = coding_sv_mis_vector(row)
                    match_vector = coding_sv_match_vector(row)
                    gene_field = row[0]
                    for gene, idxs in candidates_by_gene.items():
                        if gene not in gene_field:
                            continue
                        for idx in idxs:
                            if has_mis(mis_vector, candidates[idx].winqtl.nam_ids):
                                sv_scores[idx].append(score)
                            if match_vector:
                                match_variants[idx].append((score, match_vector))
        else:
            warnings.append(f"Missing coding SV resource for Chr{chr_id}")

    nocoding_snp_path = find_resource("combineRAPvsMSU_nocodingSNP.txt")
    if nocoding_snp_path:
        with nocoding_snp_path.open("r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if len(row) < 5:
                    continue
                gene = gene_base(row[2])
                score = parse_score(row[3])
                if score is None:
                    continue
                vector = noncoding_match_vector(row)
                for idx in candidates_by_gene.get(gene, []):
                    if noncoding_mis(row, candidates[idx].winqtl.nam_ids):
                        noncoding_snp_scores[idx].append(score)
                    if vector:
                        match_variants[idx].append((score, vector))
    else:
        warnings.append("Missing noncoding SNP resource")

    for resource_name, target in [
        ("combineRAPvsMSU_nocodingINDEL.txt", indel_scores),
        ("combineRAPvsMSU_nocodingSV.txt", sv_scores),
    ]:
        path = find_resource(resource_name)
        if not path:
            warnings.append(f"Missing {resource_name}")
            continue
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                if len(row) < 5:
                    continue
                vector = noncoding_match_vector(row)
                for gene, score in variant_gene_entries(row[3]):
                    for idx in candidates_by_gene.get(gene, []):
                        if noncoding_mis(row, candidates[idx].winqtl.nam_ids):
                            target[idx].append(score)
                        if vector:
                            match_variants[idx].append((score, vector))

    for idx in range(len(candidates)):
        sift_score = None
        if sift_raw_values[idx]:
            sift_score = sift_impact_score(min(sift_raw_values[idx]), sift_stop[idx])
        if noncoding_snp_scores[idx]:
            sift_score = max([sift_score or 0, max(noncoding_snp_scores[idx])])
        scores[idx]["sift"] = score_text(sift_score)
        scores[idx]["indel"] = score_text(max(indel_scores[idx]) if indel_scores[idx] else None)
        scores[idx]["sv"] = score_text(max(sv_scores[idx]) if sv_scores[idx] else None)
        scores[idx]["match"] = score_match(match_variants[idx], candidates[idx].winqtl.nam_ids)
    return scores, warnings


def build_winqtl_candidates(
    peaks: list[Peak],
    genes_by_chr: dict[int, list[Gene]],
    mlm_peaks: list[Peak],
    win_peaks: list[WinPeak],
    lod_rows: dict[str, tuple[str, list[str]]],
    chr_lengths: dict[int, int],
) -> list[Candidate]:
    candidates: list[Candidate] = []
    qtl_idx = 0
    for peak in peaks:
        overlaps = overlapping_winqtl_peaks(win_peaks, peak)
        if not overlaps:
            continue
        qtl_idx += 1
        region_start, region_end = candidate_region_from_winqtl(
            overlaps, chr_lengths.get(peak.chromosome, peak.position + 200_000)
        )
        winqtl = WinQtlOverlap(
            text=winqtl_overlap_text(overlaps),
            nam_ids=sorted({w.nam_id for w in overlaps}),
        )
        mlm_info = nearby_mlm_info(mlm_peaks, peak)
        for gene in genes_in_region(genes_by_chr.get(peak.chromosome, []), region_start, region_end):
            if gene_base(gene.gene_id) not in lod_rows:
                continue
            candidates.append(Candidate(qtl_idx, peak, gene, mlm_info, winqtl))
    return candidates


def build_window_candidates(
    peaks: list[Peak],
    genes_by_chr: dict[int, list[Gene]],
    mlm_peaks: list[Peak],
    window_bp: int,
) -> list[Candidate]:
    candidates: list[Candidate] = []
    fallback_nam_ids = list(range(1, 16))
    for qtl_idx, peak in enumerate(peaks, start=1):
        mlm_info = nearby_mlm_info(mlm_peaks, peak)
        winqtl = WinQtlOverlap(text="*", nam_ids=fallback_nam_ids)
        for gene in nearby_genes(genes_by_chr.get(peak.chromosome, []), peak, window_bp):
            candidates.append(Candidate(qtl_idx, peak, gene, mlm_info, winqtl))
    return candidates


def write_candidates(
    out_file: Path,
    candidates: list[Candidate],
    variant_scores: dict[int, dict[str, str]],
    lod_rows: dict[str, tuple[str, list[str]]],
    transposons: dict[str, str],
    locus_annotations: dict[str, str],
    gff_annotations: dict[str, str],
    house_keeping: set[str],
    important: set[str],
    at_info: dict[str, str],
) -> int:
    with out_file.open("w", encoding="utf-8", newline="") as out:
        out.write(OUTPUT_HEADER + "\n")
        for idx, candidate in enumerate(candidates):
            gene = candidate.gene
            base = gene_base(gene.gene_id)
            annotation = locus_annotations.get(gene.gene_id) or gff_annotations.get(gene.gene_id, "*")
            lod_pos, lod_values = lod_rows.get(base, (f"{gene.chromosome}:{gene.start}-{gene.end}", ["*"] * 15))
            scores = variant_scores[idx]
            row = [
                f"QTL{candidate.qtl_idx}",
                f"Chr{candidate.peak.chromosome}-{candidate.peak.position}",
                candidate.peak.logp,
                candidate.mlm_info,
                candidate.winqtl.text,
                gene.gene_id,
                lod_pos,
                " ".join(lod_values),
                transposons.get(gene.gene_id, "*"),
                "House_keeping" if base in house_keeping or gene.gene_id in house_keeping else "*",
                "Important" if base in important or gene.gene_id in important else "*",
                at_info.get(gene.gene_id) or at_info.get(base, "*"),
                annotation,
                scores["sift"],
                scores["indel"],
                scores["sv"],
                scores["match"],
                *gene.expression,
            ]
            out.write("\t".join(row) + "\n")
    return len(candidates)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create RAP FarmCPU peak info from FarmCPU peaks with WinQTLcart or 200kb-window candidates."
    )
    parser.add_argument("trait", help="Trait id, e.g. HNHZ")
    parser.add_argument("tissue", help="Main tissue name kept for V1 compatibility")
    parser.add_argument(
        "--candidate-mode",
        choices=["auto", "winqtl", "window", "both"],
        default="auto",
        help="Candidate interval mode. auto uses WinQTLcart if available, otherwise 200kb window.",
    )
    parser.add_argument(
        "--winqtl-trait",
        default=None,
        help="Trait key used for WinQTLcart outputs, e.g. HZ_Heading_date. Default: same as trait.",
    )
    parser.add_argument(
        "--candidate-window-bp",
        type=int,
        default=200_000,
        help="Window around each FarmCPU peak used when no WinQTL interval is available.",
    )
    parser.add_argument("--force", action="store_true", help="Regenerate even if RAP_<trait>.FarmCPUpeak_info exists.")
    args = parser.parse_args()

    peaks = load_farmcpu_peaks(args.trait)
    genes_by_chr = load_expression_genes()
    if not genes_by_chr:
        raise SystemExit(
            "Cannot build RAP table because basic_data/combineRAPvsMSU_RNAseq.txt is missing."
        )

    transposons, locus_annotations = load_locus_info()
    gff_annotations = load_gff_annotations()
    house_keeping = load_gene_set("./basic_data/house-keeping_gene.list")
    important = load_gene_set("./basic_data/repressed_geno.list")
    at_info = load_arabidopsis_info()
    mlm_peaks = load_mlm_peaks(args.trait)
    chr_lengths = load_chr_lengths()
    winqtl_trait = args.winqtl_trait or args.trait

    can_use_winqtl = Path(f"{winqtl_trait}.winqtlcart.peak").exists()
    if args.candidate_mode in ("winqtl", "both") and not can_use_winqtl:
        raise SystemExit(
            f"WinQTLcart mode requested but {winqtl_trait}.winqtlcart.peak is missing. "
            f"Run Select_Winqtlcart_Peak.py {winqtl_trait} first."
        )

    modes: list[tuple[str, str]]
    if args.candidate_mode == "both":
        modes = [("winqtl", f"{args.trait}_winqtl"), ("window", f"{args.trait}_window200kb")]
    elif args.candidate_mode == "winqtl":
        modes = [("winqtl", args.trait)]
    elif args.candidate_mode == "window":
        modes = [("window", args.trait)]
    elif can_use_winqtl:
        modes = [("winqtl", args.trait)]
    else:
        modes = [("window", args.trait)]

    for mode, output_trait in modes:
        out_file = Path(f"RAP_{output_trait}.FarmCPUpeak_info")
        if out_file.exists() and not args.force:
            print(f"[INFO] Reusing existing output: {out_file}")
            continue

        lod_rows = load_lod_rows(winqtl_trait) if mode == "winqtl" else {}
        if mode == "winqtl":
            win_peaks = load_winqtl_peaks(winqtl_trait)
            candidates = build_winqtl_candidates(peaks, genes_by_chr, mlm_peaks, win_peaks, lod_rows, chr_lengths)
        else:
            candidates = build_window_candidates(peaks, genes_by_chr, mlm_peaks, args.candidate_window_bp)

        variant_scores, warnings = compute_variant_scores(candidates)
        row_count = write_candidates(
            out_file,
            candidates,
            variant_scores,
            lod_rows,
            transposons,
            locus_annotations,
            gff_annotations,
            house_keeping,
            important,
            at_info,
        )
        for warning in sorted(set(warnings)):
            print(f"[WARN] {warning}", file=sys.stderr)
        print(f"[INFO] RAP candidates ({mode}): {row_count} -> {out_file}")


if __name__ == "__main__":
    main()
