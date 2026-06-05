#!/usr/bin/env python3
import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from textwrap import dedent


ANNOTATION_COLUMNS = [
    "Source_File",
    "Mode",
    "Annotation_locus",
    "Symbol",
    "Keyword",
    "Locus_annotation",
    "Locus_is_TE",
    "Locus_is_expressed",
    "Locus_is_representative",
    "Recommendation_score",
    "Recommendation_reason",
]


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


@dataclass
class LocusAnnotation:
    locus: str
    symbol: str = ""
    keyword: str = ""
    annotation: str = ""
    is_te: str = ""
    is_expressed: str = ""
    is_representative: str = ""
    symbols: set[str] = field(default_factory=set)
    keywords: set[str] = field(default_factory=set)


@dataclass
class AnnotatedRow:
    source_file: str
    mode: str
    row: dict[str, str]
    gene: str
    locus: str
    annotation: LocusAnnotation
    base_score: float
    recommendation_score: float = 0.0
    recommendation_reason: str = ""


def gene_base(gene_id: str) -> str:
    return re.sub(r"\.\d+$", "", (gene_id or "").strip())


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        if value is None or value == "*" or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def mode_from_path(path: Path) -> str:
    name = path.name.lower()
    if "winqtl" in name:
        return "winqtl"
    if "window200kb" in name:
        return "window200kb"
    if "window" in name:
        return "window"
    return path.stem


def output_stem(path: Path) -> str:
    name = path.name
    for suffix in [".FarmCPUpeak_info", ".tsv", ".txt", ".csv"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def merge_text(values: set[str]) -> str:
    clean = sorted(v for v in values if v)
    return ";".join(clean)


def load_keyword_annotations(path: Path) -> dict[str, LocusAnnotation]:
    annotations: dict[str, LocusAnnotation] = {}
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            locus = gene_base(row.get("locus", ""))
            if not locus:
                continue
            cur = annotations.setdefault(locus, LocusAnnotation(locus=locus))
            symbol = (row.get("Symbol") or "").strip()
            keyword = (row.get("Keyword") or "").strip()
            if symbol:
                cur.symbols.add(symbol)
            if keyword:
                cur.keywords.add(keyword)

            is_rep = (row.get("is_representative") or "").strip()
            should_replace = not cur.annotation or is_rep == "Y" or (not cur.symbol and symbol)
            if should_replace:
                cur.symbol = symbol or cur.symbol
                cur.keyword = keyword or cur.keyword
                cur.annotation = (row.get("annotation") or "").strip()
                cur.is_te = (row.get("is_TE") or "").strip()
                cur.is_expressed = (row.get("is_expressed") or "").strip()
                cur.is_representative = is_rep

    for cur in annotations.values():
        cur.symbol = merge_text(cur.symbols) or cur.symbol
        cur.keyword = merge_text(cur.keywords) or cur.keyword
    return annotations


def load_rap_msu_map(path: Path | None) -> dict[str, str]:
    if not path or not path.exists():
        return {}
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[1] == "None":
                continue
            rap = gene_base(parts[0])
            locs = [gene_base(p) for p in parts[1].split(",") if p]
            if locs:
                mapping[rap] = locs[0]
    return mapping


def resolve_locus(gene: str, annotations: dict[str, LocusAnnotation], rap_msu: dict[str, str]) -> str:
    base = gene_base(gene)
    if base in annotations:
        return base
    mapped = rap_msu.get(base)
    if mapped and mapped in annotations:
        return mapped
    return base


def variant_max(row: dict[str, str]) -> float:
    return max(
        [
            to_float(row.get("SIFT_Score")),
            to_float(row.get("INDEL_Score")),
            to_float(row.get("SV_Score")),
        ]
    )


def build_reason(row: dict[str, str], ann: LocusAnnotation, modes: set[str], rec_score: float) -> str:
    parts = [f"pipeline_score={to_float(row.get('Score')):g}"]
    farm = to_float(row.get("FarmCPU_-logP"))
    if farm:
        parts.append(f"farm_logP={farm:g}")
    vmax = variant_max(row)
    if vmax:
        parts.append(f"max_variant={vmax:g}")
    if len(modes) > 1:
        parts.append("supported_by_both_modes")
    if ann.symbol:
        parts.append(f"symbol={ann.symbol}")
    if ann.keyword:
        parts.append(f"keyword={ann.keyword}")
    match = row.get("Score_match") or row.get("Match_Score")
    if match and match != "*":
        parts.append(f"match={match}")
    parts.append(f"recommendation_score={rec_score:g}")
    return "; ".join(parts)


def row_recommendation_score(row: dict[str, str], ann: LocusAnnotation, modes: set[str]) -> float:
    score = to_float(row.get("Score"))
    vmax = variant_max(row)
    bonus = 0.0
    if len(modes) > 1:
        bonus += 5.0
    if ann.keyword:
        bonus += 5.0
    if ann.symbol:
        bonus += 2.0
    if vmax >= 20:
        bonus += 5.0
    elif vmax >= 10:
        bonus += 3.0
    match = row.get("Score_match") or row.get("Match_Score")
    if match and match != "*":
        bonus += 3.0
    return score + bonus


def read_result(path: Path, annotations: dict[str, LocusAnnotation], rap_msu: dict[str, str]) -> list[AnnotatedRow]:
    rows: list[AnnotatedRow] = []
    mode = mode_from_path(path)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("Geno", "")
            locus = resolve_locus(gene, annotations, rap_msu)
            ann = annotations.get(locus, LocusAnnotation(locus=locus))
            rows.append(
                AnnotatedRow(
                    source_file=path.name,
                    mode=mode,
                    row=row,
                    gene=gene_base(gene),
                    locus=locus,
                    annotation=ann,
                    base_score=to_float(row.get("Score")),
                )
            )
    return rows


def write_annotated(inputs: list[Path], rows_by_file: dict[str, list[AnnotatedRow]], output_dir: Path) -> None:
    for input_path in inputs:
        rows = rows_by_file[input_path.name]
        if not rows:
            continue
        original_fields = list(rows[0].row.keys())
        out_path = output_dir / f"{output_stem(input_path)}.annotated.tsv"
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t", fieldnames=original_fields + ANNOTATION_COLUMNS)
            writer.writeheader()
            for item in rows:
                ann = item.annotation
                out = dict(item.row)
                out.update(
                    {
                        "Source_File": item.source_file,
                        "Mode": item.mode,
                        "Annotation_locus": item.locus,
                        "Symbol": ann.symbol,
                        "Keyword": ann.keyword,
                        "Locus_annotation": ann.annotation,
                        "Locus_is_TE": ann.is_te,
                        "Locus_is_expressed": ann.is_expressed,
                        "Locus_is_representative": ann.is_representative,
                        "Recommendation_score": f"{item.recommendation_score:g}",
                        "Recommendation_reason": item.recommendation_reason,
                    }
                )
                writer.writerow(out)


def write_top10(all_rows: list[AnnotatedRow], output_dir: Path, prefix: str) -> Path:
    grouped: dict[str, list[AnnotatedRow]] = {}
    for item in all_rows:
        grouped.setdefault(item.locus, []).append(item)

    best_rows: list[AnnotatedRow] = []
    for items in grouped.values():
        best_rows.append(
            max(
                items,
                key=lambda item: (
                    item.recommendation_score,
                    item.base_score,
                    variant_max(item.row),
                    to_float(item.row.get("FarmCPU_-logP")),
                ),
            )
        )

    best_rows.sort(
        key=lambda item: (
            item.recommendation_score,
            item.base_score,
            variant_max(item.row),
            to_float(item.row.get("FarmCPU_-logP")),
        ),
        reverse=True,
    )
    out_path = output_dir / f"{prefix}_top10_candidate_genes.tsv"
    fields = [
        "Rank",
        "Gene",
        "Annotation_locus",
        "Modes",
        "Best_QTL",
        "Best_FarmCPU_pos",
        "Best_FarmCPU_-logP",
        "Best_pipeline_score",
        "Recommendation_score",
        "SIFT_Score",
        "INDEL_Score",
        "SV_Score",
        "Score_peak",
        "Score_annotate",
        "Score_expression",
        "Score_match",
        "Symbol",
        "Keyword",
        "Locus_annotation",
        "Recommendation_reason",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for rank, item in enumerate(best_rows[:10], start=1):
            group = grouped[item.locus]
            modes = ";".join(sorted({g.mode for g in group}))
            writer.writerow(
                {
                    "Rank": rank,
                    "Gene": item.gene,
                    "Annotation_locus": item.locus,
                    "Modes": modes,
                    "Best_QTL": item.row.get("QTL", ""),
                    "Best_FarmCPU_pos": item.row.get("FarmCPU_pos", ""),
                    "Best_FarmCPU_-logP": item.row.get("FarmCPU_-logP", ""),
                    "Best_pipeline_score": item.row.get("Score", ""),
                    "Recommendation_score": f"{item.recommendation_score:g}",
                    "SIFT_Score": item.row.get("SIFT_Score", ""),
                    "INDEL_Score": item.row.get("INDEL_Score", ""),
                    "SV_Score": item.row.get("SV_Score", ""),
                    "Score_peak": item.row.get("Score_peak", ""),
                    "Score_annotate": item.row.get("Score_annotate", ""),
                    "Score_expression": item.row.get("Score_expression", ""),
                    "Score_match": item.row.get("Score_match", ""),
                    "Symbol": item.annotation.symbol,
                    "Keyword": item.annotation.keyword,
                    "Locus_annotation": item.annotation.annotation,
                    "Recommendation_reason": item.recommendation_reason,
                }
            )
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=dedent(
            """\
            Annotate RiceG2G candidate tables and rank Top10 candidate genes.

            This script reads one or more RAP_Step2_*.FarmCPUpeak_info files, adds
            Symbol/Keyword/locus annotations from all.locus_brief_info.7.0.with_keyword.tsv,
            computes Recommendation_score, and writes a compact Top10 candidate-gene table.
            """
        ),
        formatter_class=HelpFormatter,
        epilog=dedent(
            """\
            Added row-level columns:
              Annotation_locus, Symbol, Keyword, Locus_annotation
              Locus_is_TE, Locus_is_expressed, Locus_is_representative
              Recommendation_score, Recommendation_reason

            Top10 candidate genes:
              Recommendation_score starts from the final Step2 Score and adds transparent bonuses
              for multi-mode support, Symbol/Keyword evidence, strong SIFT/INDEL/SV evidence,
              and genotype-match evidence. The reasoning is written in Recommendation_reason.

            Examples:
              Annotate one TGW window result:
                python V1/bin/annotate_recommendations.py --inputs RAP_Step2_all_TGW_quality.FarmCPUpeak_info --output-prefix all_TGW_quality

              Annotate WinQTLcart and 200kb-window results together:
                python V1/bin/annotate_recommendations.py --inputs RAP_Step2_HNHZ_winqtl.FarmCPUpeak_info RAP_Step2_HNHZ_window200kb.FarmCPUpeak_info --output-prefix HNHZ
            """
        ),
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Step2 or peak_info result files to annotate.")
    parser.add_argument(
        "--keyword-table",
        default="./basic_data/all.locus_brief_info.7.0.with_keyword.tsv",
        help="Keyword-enriched locus annotation TSV, usually basic_data/all.locus_brief_info.7.0.with_keyword.tsv.",
    )
    parser.add_argument(
        "--rap-msu-map",
        default="./basic_data/RAP-MSU_2021-11-11.txt",
        help="Optional RAP-to-MSU mapping file for Os/LOC gene id bridging.",
    )
    parser.add_argument("--output-dir", default=".", help="Directory for annotated TSV and Top10 outputs.")
    parser.add_argument("--output-prefix", default="RiceG2G", help="Prefix for <prefix>_top10_candidate_genes.tsv.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    annotations = load_keyword_annotations(Path(args.keyword_table))
    rap_msu = load_rap_msu_map(Path(args.rap_msu_map) if args.rap_msu_map else None)

    inputs = [Path(p) for p in args.inputs]
    rows_by_file: dict[str, list[AnnotatedRow]] = {}
    all_rows: list[AnnotatedRow] = []
    for path in inputs:
        rows = read_result(path, annotations, rap_msu)
        rows_by_file[path.name] = rows
        all_rows.extend(rows)

    modes_by_locus: dict[str, set[str]] = {}
    for item in all_rows:
        modes_by_locus.setdefault(item.locus, set()).add(item.mode)

    for item in all_rows:
        modes = modes_by_locus.get(item.locus, {item.mode})
        item.recommendation_score = row_recommendation_score(item.row, item.annotation, modes)
        item.recommendation_reason = build_reason(item.row, item.annotation, modes, item.recommendation_score)

    write_annotated(inputs, rows_by_file, output_dir)
    top10 = write_top10(all_rows, output_dir, args.output_prefix)
    print(f"[INFO] Annotated files: {len(inputs)}")
    print(f"[INFO] Top10 recommendations: {top10}")


if __name__ == "__main__":
    main()
