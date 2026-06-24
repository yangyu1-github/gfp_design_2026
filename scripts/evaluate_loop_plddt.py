#!/usr/bin/env python3
"""Rank Run05 loop-redesign candidates by Boltz pLDDT in designed loops.

Run05 changed only three loop modules, so whole-protein pLDDT and TM-score are
secondary audits. This evaluator reads Boltz per-residue pLDDT arrays and ranks
candidates by mean pLDDT over the designed loop windows.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from statistics import mean
from typing import Iterable

import numpy as np


DEFAULT_WINDOWS = "129-147,188-198,209-216"


def parse_windows(spec: str) -> list[tuple[int, int]]:
    windows: list[tuple[int, int]] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" not in chunk:
            pos = int(chunk)
            windows.append((pos, pos))
            continue
        start_s, end_s = chunk.split("-", 1)
        start, end = int(start_s), int(end_s)
        if start < 1 or end < start:
            raise ValueError(f"Invalid loop window: {chunk}")
        windows.append((start, end))
    if not windows:
        raise ValueError("At least one loop window is required")
    return windows


def window_positions(windows: Iterable[tuple[int, int]]) -> list[int]:
    positions: list[int] = []
    for start, end in windows:
        positions.extend(range(start, end + 1))
    return positions


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def sanitize(header: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", header)[:180]


def sequence_for_prediction(name: str, seqs: dict[str, str]) -> str:
    by_sanitized = {sanitize(header): seq for header, seq in seqs.items()}
    exact = by_sanitized.get(name)
    if exact is not None:
        return exact
    for header, seq in seqs.items():
        sanitized = sanitize(header)
        if name.startswith(sanitized[:80]) or sanitized.startswith(name[:80]):
            return seq
    return ""


def find_prediction_dirs(boltz_out: Path) -> list[Path]:
    pred_dirs: list[Path] = []
    for pred_root in boltz_out.rglob("predictions"):
        pred_dirs.extend([p for p in pred_root.iterdir() if p.is_dir()])
    return sorted(set(pred_dirs))


def load_confidence(pred_dir: Path) -> dict:
    conf_files = sorted(pred_dir.glob("confidence_*_model_0.json"))
    if not conf_files:
        return {}
    return json.loads(conf_files[0].read_text(encoding="utf-8"))


def load_plddt(pred_dir: Path) -> np.ndarray | None:
    plddt_files = sorted(pred_dir.glob("plddt_*_model_0.npz"))
    if not plddt_files:
        return None
    data = np.load(plddt_files[0])
    return np.asarray(data["plddt"], dtype=float)


def model_path(pred_dir: Path) -> str:
    files = sorted(pred_dir.glob("*_model_0.pdb")) or sorted(pred_dir.glob("*_model_0.cif"))
    return str(files[0]) if files else ""


def window_values(plddt: np.ndarray, start: int, end: int) -> list[float]:
    values: list[float] = []
    for pos in range(start, end + 1):
        idx = pos - 1
        if 0 <= idx < len(plddt):
            values.append(float(plddt[idx]))
    return values


def analyze(args: argparse.Namespace) -> list[dict]:
    windows = parse_windows(args.loop_windows)
    loop_positions = window_positions(windows)
    seqs = read_fasta(args.fasta) if args.fasta else {}

    rows: list[dict] = []
    for pred_dir in find_prediction_dirs(args.boltz_out):
        plddt = load_plddt(pred_dir)
        if plddt is None:
            continue

        all_loop_values = [
            float(plddt[pos - 1])
            for pos in loop_positions
            if 0 <= pos - 1 < len(plddt)
        ]
        if not all_loop_values:
            continue

        conf = load_confidence(pred_dir)
        row = {
            "name": pred_dir.name,
            "loop_plddt_mean": mean(all_loop_values),
            "loop_plddt_min": min(all_loop_values),
            "whole_plddt_mean": float(plddt.mean()),
            "confidence": conf.get("confidence_score", 0),
            "ptm": conf.get("ptm", 0),
            "model": model_path(pred_dir),
            "sequence": sequence_for_prediction(pred_dir.name, seqs),
        }
        for start, end in windows:
            vals = window_values(plddt, start, end)
            label = f"loop_{start}_{end}"
            row[f"{label}_mean"] = mean(vals) if vals else 0.0
            row[f"{label}_min"] = min(vals) if vals else 0.0
        rows.append(row)

    rows.sort(
        key=lambda r: (
            r["loop_plddt_mean"],
            r["loop_plddt_min"],
            r["whole_plddt_mean"],
            r["confidence"],
        ),
        reverse=True,
    )
    return rows


def write_csv(rows: list[dict], out_csv: Path, windows: list[tuple[int, int]]) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank",
        "name",
        "loop_plddt_mean",
        "loop_plddt_min",
        "whole_plddt_mean",
        "confidence",
        "ptm",
    ]
    for start, end in windows:
        fields.extend([f"loop_{start}_{end}_mean", f"loop_{start}_{end}_min"])
    fields.extend(["model", "sequence"])

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, 1):
            out = {"rank": rank, **row}
            writer.writerow(out)


def write_fasta(rows: list[dict], out_fasta: Path, top_n: int) -> None:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as handle:
        for row in rows[:top_n]:
            seq = row.get("sequence", "")
            if not seq:
                continue
            header = (
                f"{row['name']} loop_plddt={row['loop_plddt_mean']:.4f} "
                f"loop_min={row['loop_plddt_min']:.4f} "
                f"whole_plddt={row['whole_plddt_mean']:.4f} "
                f"conf={float(row['confidence']):.4f}"
            )
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boltz-out", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--top-fasta", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=48)
    parser.add_argument("--loop-windows", default=DEFAULT_WINDOWS)
    args = parser.parse_args()

    windows = parse_windows(args.loop_windows)
    rows = analyze(args)
    write_csv(rows, args.out_csv, windows)
    write_fasta(rows, args.top_fasta, args.top_n)

    print(f"Predictions analyzed: {len(rows)}")
    print(f"Loop windows: {args.loop_windows}")
    print(f"Wrote loop pLDDT CSV: {args.out_csv}")
    print(f"Wrote top loop pLDDT FASTA: {args.top_fasta}")
    print("Top 20 by designed-loop pLDDT:")
    for row in rows[:20]:
        print(
            f"{row['name']:<82} "
            f"loop_pLDDT={row['loop_plddt_mean']:.3f} "
            f"loop_min={row['loop_plddt_min']:.3f} "
            f"whole_pLDDT={row['whole_plddt_mean']:.3f} "
            f"conf={float(row['confidence']):.3f}"
        )


if __name__ == "__main__":
    main()
