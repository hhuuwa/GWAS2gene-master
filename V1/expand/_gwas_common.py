#!/usr/bin/env python3
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_P_THRESHOLD = 1.67e-8
DEFAULT_SIGNIFICANT_P_THRESHOLD = 1e-5


@dataclass(frozen=True)
class GwasRecord:
    snp: str
    chromosome: int
    position: int
    pvalue: float
    logp: float
    ref: str = "*"
    alt: str = "*"
    effect: str = "*"
    se: str = "*"


def clean_header(value: str) -> str:
    return value.strip().strip('"')


def normalized(value: str) -> str:
    return clean_header(value).strip().lower()


def parse_chr(value: str) -> int:
    text = clean_header(value)
    match = re.search(r"(\d+)", text)
    if not match:
        raise ValueError(f"Cannot parse chromosome from {value!r}")
    return int(match.group(1))


def parse_int_position(value: str) -> int:
    return int(float(clean_header(value)))


def parse_pvalue(value: str) -> float | None:
    text = clean_header(value)
    if not text or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def log10_score(pvalue: float) -> float:
    if pvalue <= 0:
        return 309.0
    value = -math.log10(pvalue)
    if math.isinf(value):
        return 309.0
    return value


def find_named_column(header: list[str], names: Iterable[str]) -> int | None:
    wanted = {name.lower() for name in names}
    for idx, name in enumerate(header):
        if normalized(name) in wanted:
            return idx
    return None


def find_pvalue_column(header: list[str], trait: str, method: str) -> int:
    method_norm = method.lower()
    exact = {
        f"{trait}.{method}".lower(),
        f"{trait}_{method}".lower(),
        method_norm,
        "p",
        "pvalue",
        "p.value",
        "p-value",
    }
    for idx, name in enumerate(header):
        norm = normalized(name)
        if norm in exact or norm.endswith(f".{method_norm}"):
            return idx

    metadata = {
        "snp",
        "marker",
        "chrom",
        "chr",
        "chromosome",
        "pos",
        "position",
        "ref",
        "alt",
        "maf",
        "effect",
        "se",
    }
    for idx in range(len(header) - 1, -1, -1):
        if normalized(header[idx]) not in metadata:
            return idx

    raise ValueError(f"Cannot identify p-value column in header: {header}")


def find_gwas_input(trait: str, method: str) -> tuple[Path, bool]:
    method = method.strip()
    signals = method.lower() == "farmcpu"
    candidates: list[tuple[str, bool]] = []
    if signals:
        candidates.extend(
            [
                (f"./data_flowering/{trait}.{method}_signals.csv", True),
                (f"./data/farmCPU/{trait}.{method}_signals.csv", True),
                (f"./data/farmCPU/{trait}.{method}/{trait}.{method}_signals.csv", True),
            ]
        )
    candidates.extend(
        [
            (f"./data/farmCPU/{trait}.{method}.csv", False),
            (f"./data/farmCPU/{trait}.{method}/{trait}.{method}.csv", False),
            (f"./data_flowering/{trait}.{method}.csv", False),
        ]
    )

    exact_seen = {path for path, _ in candidates}
    for path in sorted(Path(".").glob(f"data_*/{trait}.{method}_signals.csv")):
        text = str(path)
        if text not in exact_seen:
            candidates.append((text, True))
            exact_seen.add(text)
    for path in sorted(Path(".").glob(f"data_*/{trait}.{method}.csv")):
        text = str(path)
        if text not in exact_seen:
            candidates.append((text, False))
            exact_seen.add(text)

    trait_dir = Path(f"data_{trait}")
    if trait_dir.exists():
        signal_matches = sorted(trait_dir.glob(f"*.{method}_signals.csv"))
        full_matches = sorted(trait_dir.glob(f"*.{method}.csv"))
        if len(signal_matches) == 1:
            candidates.append((str(signal_matches[0]), True))
        if len(full_matches) == 1:
            candidates.append((str(full_matches[0]), False))

    for path, is_signals in candidates:
        if os.path.exists(path):
            return Path(path), is_signals
    checked = "\n  - ".join(path for path, _ in candidates)
    raise FileNotFoundError(f"Cannot find {method} input for {trait}. Checked:\n  - {checked}")


def read_gwas_records(
    path: Path,
    trait: str,
    method: str,
    p_threshold: float | None = DEFAULT_P_THRESHOLD,
) -> list[GwasRecord]:
    records: list[GwasRecord] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        raw_header = next(reader, None)
        if not raw_header:
            return records
        header = [clean_header(col) for col in raw_header]
        snp_col = find_named_column(header, ["SNP", "Marker"])
        chr_col = find_named_column(header, ["CHROM", "CHR", "Chromosome"])
        pos_col = find_named_column(header, ["POS", "Position"])
        p_col = find_pvalue_column(header, trait, method)
        ref_col = find_named_column(header, ["REF"])
        alt_col = find_named_column(header, ["ALT"])
        effect_col = find_named_column(header, ["Effect"])
        se_col = find_named_column(header, ["SE"])
        if snp_col is None or chr_col is None or pos_col is None:
            raise ValueError(f"Missing SNP/CHROM/POS columns in {path}")

        for row in reader:
            if len(row) <= max(snp_col, chr_col, pos_col, p_col):
                continue
            pvalue = parse_pvalue(row[p_col])
            if pvalue is None:
                continue
            if p_threshold is not None and pvalue > p_threshold:
                continue
            try:
                chromosome = parse_chr(row[chr_col])
                position = parse_int_position(row[pos_col])
            except ValueError:
                continue
            records.append(
                GwasRecord(
                    snp=clean_header(row[snp_col]),
                    chromosome=chromosome,
                    position=position,
                    pvalue=pvalue,
                    logp=log10_score(pvalue),
                    ref=clean_header(row[ref_col]) if ref_col is not None and ref_col < len(row) else "*",
                    alt=clean_header(row[alt_col]) if alt_col is not None and alt_col < len(row) else "*",
                    effect=clean_header(row[effect_col]) if effect_col is not None and effect_col < len(row) else "*",
                    se=clean_header(row[se_col]) if se_col is not None and se_col < len(row) else "*",
                )
            )
    return records


def select_peak_records(records: list[GwasRecord], window_bp: int) -> list[GwasRecord]:
    if not records:
        return []

    selected: list[GwasRecord] = []
    current: list[GwasRecord] = []
    last_chr: int | None = None
    last_pos: int | None = None

    for record in sorted(records, key=lambda r: (r.chromosome, r.position, r.pvalue)):
        if (
            current
            and last_chr == record.chromosome
            and last_pos is not None
            and record.position - last_pos <= window_bp
        ):
            current.append(record)
        else:
            if current:
                selected.append(min(current, key=lambda r: r.pvalue))
            current = [record]
        last_chr = record.chromosome
        last_pos = record.position

    if current:
        selected.append(min(current, key=lambda r: r.pvalue))
    return selected
