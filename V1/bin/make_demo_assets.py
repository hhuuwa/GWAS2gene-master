#!/usr/bin/env python3
import argparse
import csv
import html
import os
from collections import defaultdict


COMPONENT_FIELDS = [
    "SIFT_Score",
    "INDEL_Score",
    "SV_Score",
    "Score_peak",
    "Score_TE",
    "Score_annotate",
    "Score_expression",
    "Score_match",
    "Score",
]


def to_num(v: str) -> float:
    try:
        return float(v)
    except Exception:
        return 0.0


def read_rows(path: str):
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row["_score"] = to_num(row.get("Score", "0"))
            rows.append(row)
    return rows


def read_gene_keywords(path: str):
    info = {}
    if not os.path.exists(path):
        return info
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("MSU") or "").strip()
            if not gene:
                continue
            info[gene] = {
                "Symbol": (row.get("Symbol") or "").strip(),
                "Keyword": (row.get("Keyword") or "").strip(),
            }
    return info


def attach_gene_keywords(rows, gene_info):
    for row in rows:
        gene = row.get("Geno", "").strip()
        info = gene_info.get(gene, {})
        row["Symbol"] = info.get("Symbol", "")
        row["Keyword"] = info.get("Keyword", "")


def write_tsv(path: str, rows, fields):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def make_bar_svg(labels, values, title, width=1180, height=520, bar_color="#216869"):
    margin_left = 360
    margin_top = 60
    margin_bottom = 40
    bar_gap = 8
    n = max(len(values), 1)
    plot_h = height - margin_top - margin_bottom
    bar_h = max(10, int((plot_h - bar_gap * max(0, n - 1)) / n))
    plot_w = width - margin_left - 60
    vmax = max(values) if values else 1.0
    if vmax <= 0:
        vmax = 1.0

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        "<style>"
        "text{font-family:Segoe UI,Arial,sans-serif;font-size:12px;fill:#14213d}"
        ".title{font-size:18px;font-weight:700}"
        ".axis{fill:#6b7280}"
        ".grid{stroke:#e5e7eb;stroke-width:1}"
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fffdf8"/>',
        f'<text x="22" y="30" class="title">{html.escape(title)}</text>',
    ]

    for t in [0.25, 0.5, 0.75, 1.0]:
        x = margin_left + plot_w * t
        out.append(
            f'<line class="grid" x1="{x:.1f}" y1="{margin_top}" x2="{x:.1f}" y2="{height - margin_bottom}"/>'
        )
        out.append(f'<text class="axis" x="{x-12:.1f}" y="{height-12}">{vmax * t:.1f}</text>')

    for idx, (label, value) in enumerate(zip(labels, values)):
        y = margin_top + idx * (bar_h + bar_gap)
        bar_w = (value / vmax) * plot_w
        out.append(
            f'<rect x="{margin_left}" y="{y}" width="{bar_w:.1f}" height="{bar_h}" rx="3" fill="{bar_color}"/>'
        )
        out.append(f'<text x="12" y="{y + bar_h - 2}">{html.escape(label)}</text>')
        out.append(f'<text x="{margin_left + bar_w + 8:.1f}" y="{y + bar_h - 2}">{value:.2f}</text>')

    out.append("</svg>")
    return "\n".join(out)


def color_for_value(value: float, vmax: float) -> str:
    if vmax <= 0:
        vmax = 1.0
    ratio = max(0.0, min(1.0, value / vmax))
    # warm sandstone -> deep teal
    r = int(245 - 142 * ratio)
    g = int(238 - 85 * ratio)
    b = int(220 - 110 * ratio)
    return f"rgb({r},{g},{b})"


def make_heatmap_svg(rows, title):
    labels = []
    matrix = []
    for row in rows:
        symbol = row.get("Symbol", "")
        gene = row.get("Geno", "")
        label = f"{row.get('QTL', '?')} | {gene}"
        if symbol:
            label = f"{label} ({symbol})"
        labels.append(label)
        matrix.append([to_num(row.get(field, "0")) for field in COMPONENT_FIELDS])

    vmax = max([max(r) for r in matrix], default=1.0)
    cell_w = 108
    cell_h = 28
    left = 400
    top = 80
    width = left + cell_w * len(COMPONENT_FIELDS) + 30
    height = top + cell_h * max(len(rows), 1) + 50

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        "<style>"
        "text{font-family:Segoe UI,Arial,sans-serif;fill:#14213d}"
        ".title{font-size:18px;font-weight:700}"
        ".head{font-size:12px;font-weight:700}"
        ".label{font-size:12px}"
        ".value{font-size:11px}"
        ".border{stroke:#ffffff;stroke-width:1}"
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fffdf8"/>',
        f'<text x="22" y="30" class="title">{html.escape(title)}</text>',
    ]

    for col, field in enumerate(COMPONENT_FIELDS):
        x = left + col * cell_w + cell_w / 2
        out.append(f'<text x="{x:.1f}" y="58" text-anchor="middle" class="head">{html.escape(field)}</text>')

    for row_idx, (label, values) in enumerate(zip(labels, matrix)):
        y = top + row_idx * cell_h
        out.append(f'<text x="12" y="{y + 19}" class="label">{html.escape(label)}</text>')
        for col_idx, value in enumerate(values):
            x = left + col_idx * cell_w
            fill = color_for_value(value, vmax)
            out.append(
                f'<rect x="{x}" y="{y}" width="{cell_w-2}" height="{cell_h-2}" fill="{fill}" class="border" rx="3"/>'
            )
            text_color = "#0f172a" if value < vmax * 0.55 else "#ffffff"
            out.append(
                f'<text x="{x + cell_w/2 - 1:.1f}" y="{y + 18}" text-anchor="middle" class="value" fill="{text_color}">{value:.2f}</text>'
            )

    out.append("</svg>")
    return "\n".join(out)


def build_qtl_summary(rows, gene_info):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.get("QTL", "")].append(row)

    summary = []
    for qtl in sorted(grouped.keys()):
        cur = grouped[qtl]
        best = max(cur, key=lambda r: r["_score"])
        gene = best.get("Geno", "")
        info = gene_info.get(gene, {})
        summary.append(
            {
                "QTL": qtl,
                "Candidate_Gene_Count": len(cur),
                "Best_Gene": gene,
                "Best_Gene_Symbol": info.get("Symbol", ""),
                "Best_Gene_Keyword": info.get("Keyword", ""),
                "Max_Score": f"{best['_score']:.3f}",
            }
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build demo tables and rich visualizations for RiceG2G V1.")
    parser.add_argument("--input", required=True, help="Step2 result file")
    parser.add_argument("--output-dir", required=True, help="Directory to write demo assets")
    parser.add_argument("--topn", type=int, default=20, help="Top-N genes to summarize")
    parser.add_argument(
        "--gene-keywords",
        default=os.path.join("V1", "source", "rice_new_keyword.txt"),
        help="MSU/Symbol/Keyword annotation file",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    rows = read_rows(args.input)
    if not rows:
        raise SystemExit(f"No rows found in {args.input}")

    gene_info = read_gene_keywords(args.gene_keywords)
    attach_gene_keywords(rows, gene_info)

    rows_sorted = sorted(rows, key=lambda r: r["_score"], reverse=True)
    top_rows = rows_sorted[: args.topn]

    top_fields = [
        "QTL",
        "FarmCPU_pos",
        "Geno",
        "Symbol",
        "Keyword",
        "Geno_pos",
        "SIFT_Score",
        "INDEL_Score",
        "SV_Score",
        "Score_peak",
        "Score_TE",
        "Score_annotate",
        "Score_expression",
        "Score_match",
        "Score",
    ]
    write_tsv(os.path.join(args.output_dir, "top_candidates.tsv"), top_rows, top_fields)

    qtl_rows = build_qtl_summary(rows, gene_info)
    qtl_fields = [
        "QTL",
        "Candidate_Gene_Count",
        "Best_Gene",
        "Best_Gene_Symbol",
        "Best_Gene_Keyword",
        "Max_Score",
    ]
    write_tsv(os.path.join(args.output_dir, "qtl_summary.tsv"), qtl_rows, qtl_fields)

    score_labels = []
    score_values = []
    for row in top_rows:
        gene = row.get("Geno", "")
        symbol = row.get("Symbol", "")
        label = f"{row.get('QTL','?')} | {gene}"
        if symbol:
            label += f" ({symbol})"
        score_labels.append(label)
        score_values.append(row["_score"])

    top_svg = make_bar_svg(score_labels, score_values, f"Top {len(top_rows)} Candidate Genes by Total Score")
    with open(os.path.join(args.output_dir, "top_candidates.svg"), "w", encoding="utf-8") as f:
        f.write(top_svg)

    heatmap_svg = make_heatmap_svg(top_rows, "Top Candidate Genes: Multi-Score Heatmap")
    with open(os.path.join(args.output_dir, "score_components_heatmap.svg"), "w", encoding="utf-8") as f:
        f.write(heatmap_svg)

    qtl_svg = make_bar_svg(
        [row["QTL"] for row in qtl_rows],
        [to_num(row["Max_Score"]) for row in qtl_rows],
        "QTL Max Score Overview",
        width=1020,
        height=460,
        bar_color="#b56576",
    )
    with open(os.path.join(args.output_dir, "qtl_max_score.svg"), "w", encoding="utf-8") as f:
        f.write(qtl_svg)

    report_path = os.path.join(args.output_dir, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(
            "<!doctype html><html><head><meta charset='utf-8'><title>RiceG2G V1 Demo Report</title>"
            "<style>"
            "body{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#14213d;background:#fffdf8}"
            "h1{margin:0 0 8px 0} h2{margin:0 0 12px 0}"
            "p{color:#5c677d}"
            "code{background:#f3efe5;padding:2px 6px;border-radius:4px}"
            ".card{border:1px solid #e9dcc9;border-radius:12px;padding:16px;margin:16px 0;background:#fffaf0}"
            "a{color:#216869;text-decoration:none}"
            "a:hover{text-decoration:underline}"
            "</style></head><body>"
        )
        f.write("<h1>RiceG2G V1 Demo Report</h1>")
        f.write(f"<p>Input result file: <code>{html.escape(args.input)}</code></p>")
        f.write("<div class='card'><h2>Top Candidate Genes by Total Score</h2>")
        f.write("<p>This panel ranks candidate genes using the integrated RiceG2G final score.</p>")
        f.write("<img src='top_candidates.svg' style='max-width:100%;height:auto;'></div>")
        f.write("<div class='card'><h2>Multi-Score Explanation Panel</h2>")
        f.write("<p>Each row is a top candidate gene. Each column shows one evidence component used in prioritization.</p>")
        f.write("<img src='score_components_heatmap.svg' style='max-width:100%;height:auto;'></div>")
        f.write("<div class='card'><h2>QTL Overview</h2>")
        f.write("<p>The strongest candidate score detected in each QTL region.</p>")
        f.write("<img src='qtl_max_score.svg' style='max-width:100%;height:auto;'></div>")
        f.write("<div class='card'><h2>Files</h2><ul>")
        f.write("<li><a href='top_candidates.tsv'>top_candidates.tsv</a></li>")
        f.write("<li><a href='qtl_summary.tsv'>qtl_summary.tsv</a></li>")
        f.write("<li><a href='top_candidates.svg'>top_candidates.svg</a></li>")
        f.write("<li><a href='score_components_heatmap.svg'>score_components_heatmap.svg</a></li>")
        f.write("<li><a href='qtl_max_score.svg'>qtl_max_score.svg</a></li>")
        f.write("</ul></div>")
        f.write("</body></html>")

    print("[DONE] Demo assets generated at:", os.path.abspath(args.output_dir))
    print("[DONE] Report:", os.path.abspath(report_path))


if __name__ == "__main__":
    main()
