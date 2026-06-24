#!/usr/bin/env python3
"""Utilities for Run09 loop-seeded barrel PSSM refinement."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

import numpy as np
from Bio import BiopythonWarning
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser


AA20 = set("ACDEFGHIKLMNPQRSTVWY")
LOOP_WINDOWS = [(129, 147), (188, 198), (209, 216)]
CHROMOPHORE_POSITIONS = {65, 66, 67}
BETA_RANGES = [
    (12, 22),
    (25, 36),
    (40, 48),
    (92, 100),
    (103, 115),
    (118, 128),
    (148, 155),
    (159, 171),
    (174, 187),
    (199, 208),
    (217, 227),
]
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def slug(text: str, limit: int = 150) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)[:limit].strip("_")


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def beta_positions() -> set[int]:
    return {pos for start, end in BETA_RANGES for pos in range(start, end + 1)}


def loop_positions() -> set[int]:
    return {pos for start, end in LOOP_WINDOWS for pos in range(start, end + 1)}


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def read_fasta_dict(path: Path) -> dict[str, str]:
    return {header.split()[0]: seq for header, seq in read_fasta_records(path)}


def read_named_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        elif not line.startswith("#"):
            chunks.append("".join(c for c in line.upper() if c in AA20))
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{wrap(seq)}\n")


def load_exclusion_sequences(path: Path) -> set[str]:
    seqs: set[str] = set()
    if not path.exists():
        return seqs
    text = path.read_text(encoding="utf-8", errors="ignore")
    for match in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]{50,}", text):
        seqs.add(match.group(0))
    return seqs


def global_align_maps(ref: str, query: str) -> tuple[dict[int, int], dict[int, int]]:
    """Needleman-Wunsch mapping, 1-based positions, excluding gaps."""
    n = len(ref)
    m = len(query)
    gap = -10
    match = 2
    mismatch = -1
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
        trace[i][0] = "U"
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
        trace[0][j] = "L"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + (match if ref[i - 1] == query[j - 1] else mismatch)
            up = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap
            best = max(diag, up, left)
            dp[i][j] = best
            trace[i][j] = "D" if best == diag else ("U" if best == up else "L")
    i, j = n, m
    ref_to_query: dict[int, int] = {}
    query_to_ref: dict[int, int] = {}
    while i > 0 or j > 0:
        t = trace[i][j]
        if t == "D":
            ref_to_query[i] = j
            query_to_ref[j] = i
            i -= 1
            j -= 1
        elif t == "U":
            i -= 1
        elif t == "L":
            j -= 1
        else:
            raise RuntimeError("Alignment traceback failed")
    return ref_to_query, query_to_ref


def coordinate_maps(ref: str, query: str) -> tuple[dict[int, int], dict[int, int]]:
    """Return coordinate maps, using direct maps for fixed-length designs."""
    if len(ref) == len(query):
        direct = {i: i for i in range(1, len(ref) + 1)}
        return direct, direct.copy()
    return global_align_maps(ref, query)


def parse_pdb_residues(path: Path) -> tuple[str, list[tuple[str, int, str]]]:
    seq: list[str] = []
    residues: list[tuple[str, int, str]] = []
    seen: set[tuple[str, int, str]] = set()
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if resname not in THREE_TO_ONE:
            continue
        chain = line[21].strip() or " "
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        key = (chain, resseq, icode)
        if key in seen:
            continue
        seen.add(key)
        residues.append(key)
        seq.append(THREE_TO_ONE[resname])
    return "".join(seq), residues


def trim_pdb_by_keys(source: Path, dest: Path, keep_keys: set[tuple[str, int, str]]) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with source.open(encoding="utf-8", errors="ignore") as inp, dest.open("w", encoding="utf-8") as out:
        for line in inp:
            if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
                chain = line[21].strip() or " "
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                icode = line[26].strip()
                if (chain, resseq, icode) in keep_keys:
                    out.write(line)
            elif line.startswith(("MODEL", "ENDMDL", "REMARK")):
                out.write(line)
        out.write("TER\nEND\n")


def position_list_csv(value: str) -> list[int]:
    return [int(x) for x in value.split(";") if x]


def mutation_positions(ref: str, seq: str) -> list[int]:
    return [i for i, (a, b) in enumerate(zip(ref, seq), start=1) if a != b]


def longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def longest_dipeptide_repeat(seq: str) -> int:
    best = 0
    for i in range(max(0, len(seq) - 3)):
        di = seq[i : i + 2]
        units = 1
        j = i + 2
        while j + 1 < len(seq) and seq[j : j + 2] == di:
            units += 1
            j += 2
        best = max(best, units)
    return best


def low_complexity_fraction(seq: str, window: int = 12, entropy_thresh: float = 2.0) -> float:
    if len(seq) < window:
        return 0.0
    lows = 0
    total = 0
    for i in range(len(seq) - window + 1):
        w = seq[i : i + window]
        counts = Counter(w)
        entropy = 0.0
        for count in counts.values():
            p = count / len(w)
            entropy -= p * math.log(p, 2)
        total += 1
        if entropy < entropy_thresh:
            lows += 1
    return lows / total if total else 0.0


def net_charge(seq: str) -> int:
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def parse_float(pattern: str, text: str) -> float | None:
    match = re.search(pattern, text)
    return float(match.group(1)) if match else None


def parse_temperature(header: str) -> str:
    match = re.search(r"T=([0-9.]+)", header)
    if match:
        return f"{float(match.group(1)):.2f}"
    match = re.search(r"_T([0-9p]+)__seq", header)
    if match:
        return match.group(1).replace("p", ".")
    return "unknown"


def parse_backbone(header: str) -> str:
    token = header.split()[0]
    match = re.match(r"(.+)_T[0-9p]+__seq\d+$", token)
    if match:
        return match.group(1)
    return token.split("__seq")[0]


def jaccard(a: set[int], b: set[int]) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b) if (a | b) else 0.0


def prepare_inputs(args: argparse.Namespace) -> None:
    args.input_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir = args.input_dir / "trimmed_backbones"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    args.pssm_dir.mkdir(parents=True, exist_ok=True)

    sf_seq = read_named_fasta(args.ref_fasta)["sfGFP"]
    with args.loop_csv.open(newline="", encoding="utf-8") as handle:
        loop_rows = list(csv.DictReader(handle))[: args.top_loop_n]
    if len(loop_rows) != args.top_loop_n:
        raise SystemExit(f"Expected {args.top_loop_n} Run05 loop donors, found {len(loop_rows)}")

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest_by_name = {row["name"]: row for row in manifest_rows}

    metadata: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str]] = []
    mapping_rows: list[dict[str, str | int]] = []
    required_positions = sorted(CHROMOPHORE_POSITIONS | loop_positions() | {args.final_length})

    for loop_row in loop_rows:
        donor_name = loop_row["name"]
        donor_rank = int(loop_row["rank"])
        donor_seq = loop_row["sequence"].strip().upper()
        if donor_name not in manifest_by_name:
            raise SystemExit(f"Run05 donor not found in pLDDT manifest: {donor_name}")
        manifest_row = manifest_by_name[donor_name]
        manifest_seq = manifest_row["sequence"].strip().upper()
        if donor_seq != manifest_seq:
            raise SystemExit(f"Loop CSV sequence and manifest sequence differ for {donor_name}")

        sf_to_donor, donor_to_sf = global_align_maps(sf_seq, donor_seq)
        for pos in required_positions:
            donor_idx = sf_to_donor.get(pos)
            if donor_idx is None:
                raise SystemExit(f"sfGFP position {pos} does not map to donor {donor_name}")
        if "".join(donor_seq[sf_to_donor[p] - 1] for p in [65, 66, 67]) != "TYG":
            raise SystemExit(f"Donor {donor_name} does not map exact TYG65-67")

        source_pdb = Path(manifest_row["pdb_file"])
        if not source_pdb.is_absolute():
            source_pdb = args.root_dir / source_pdb
        if not source_pdb.exists():
            raise SystemExit(f"Source PDB missing for {donor_name}: {source_pdb}")

        pdb_seq, pdb_residues = parse_pdb_residues(source_pdb)
        donor_to_pdb, pdb_to_donor = coordinate_maps(donor_seq, pdb_seq)
        donor_idx_to_pdb_key: dict[int, tuple[str, int, str]] = {}
        for donor_idx, pdb_idx in donor_to_pdb.items():
            donor_idx_to_pdb_key[donor_idx] = pdb_residues[pdb_idx - 1]

        keep_donor_indices: list[int] = []
        for donor_idx in range(1, len(donor_seq) + 1):
            sf_pos = donor_to_sf.get(donor_idx)
            if sf_pos is not None and sf_pos <= args.final_length:
                keep_donor_indices.append(donor_idx)
                continue
            if sf_pos is None:
                prev_sf = max((donor_to_sf.get(i, 0) or 0) for i in range(1, donor_idx))
                next_sf_values = [donor_to_sf.get(i) for i in range(donor_idx + 1, len(donor_seq) + 1) if donor_to_sf.get(i)]
                next_sf = min(next_sf_values) if next_sf_values else None
                if prev_sf <= args.final_length and (next_sf is None or next_sf <= args.final_length + 1):
                    keep_donor_indices.append(donor_idx)

        keep_keys: set[tuple[str, int, str]] = set()
        sf_to_pdb_key: dict[int, tuple[str, int, str]] = {}
        donor_idx_to_trimmed_idx: dict[int, int] = {}
        for trimmed_idx, donor_idx in enumerate(keep_donor_indices, start=1):
            donor_idx_to_trimmed_idx[donor_idx] = trimmed_idx
            key = donor_idx_to_pdb_key.get(donor_idx)
            if key is None:
                raise SystemExit(f"Donor index {donor_idx} does not map to PDB in {donor_name}")
            keep_keys.add(key)
            sf_pos = donor_to_sf.get(donor_idx)
            if sf_pos is not None:
                sf_to_pdb_key[sf_pos] = key

        for pos in required_positions:
            if pos not in sf_to_pdb_key:
                raise SystemExit(f"sfGFP position {pos} does not map to a PDB residue in {donor_name}")

        trimmed_seq = "".join(donor_seq[donor_idx - 1] for donor_idx in keep_donor_indices)
        trimmed_index_to_sf = [donor_to_sf.get(donor_idx) for donor_idx in keep_donor_indices]
        tyg_indices = [donor_idx_to_trimmed_idx[sf_to_donor[pos]] for pos in (65, 66, 67)]
        if "".join(trimmed_seq[i - 1] for i in tyg_indices) != "TYG":
            raise SystemExit(f"Trimmed donor does not preserve TYG65-67: {donor_name}")

        backbone_id = f"run09_L{donor_rank:02d}_{slug(donor_name, 45)}"
        trimmed_pdb = trimmed_dir / f"{backbone_id}_trim{args.final_length}.pdb"
        trim_pdb_by_keys(source_pdb, trimmed_pdb, keep_keys)
        trimmed_pdb_seq, trimmed_pdb_residues = parse_pdb_residues(trimmed_pdb)
        if trimmed_pdb_seq != trimmed_seq:
            raise SystemExit(f"Trimmed PDB sequence does not match trimmed donor sequence for {donor_name}")
        if len(trimmed_pdb_residues) != len(trimmed_seq):
            raise SystemExit(f"Trimmed PDB residue count does not match trimmed sequence: {donor_name}")

        missing_sf_positions = [pos for pos in range(1, args.final_length + 1) if pos not in sf_to_donor]
        inherited_mutations = []
        for trimmed_idx, sf_pos in enumerate(trimmed_index_to_sf, start=1):
            if sf_pos is not None and sf_pos <= args.final_length and trimmed_seq[trimmed_idx - 1] != sf_seq[sf_pos - 1]:
                inherited_mutations.append(sf_pos)
        inherited_source_differences = sorted(set(inherited_mutations + missing_sf_positions))

        loop_strings: dict[str, str] = {}
        loop_index_fields: dict[str, str] = {}
        for start, end in LOOP_WINDOWS:
            donor_start = sf_to_donor[start]
            donor_end = sf_to_donor[end]
            trim_start = donor_idx_to_trimmed_idx[donor_start]
            trim_end = donor_idx_to_trimmed_idx[donor_end]
            loop_strings[f"loop_{start}_{end}"] = trimmed_seq[trim_start - 1 : trim_end]
            loop_index_fields[f"loop_{start}_{end}_seq_start"] = str(trim_start)
            loop_index_fields[f"loop_{start}_{end}_seq_end"] = str(trim_end)
        metadata.append(
            {
                "backbone_id": backbone_id,
                "loop_rank": str(donor_rank),
                "donor_name": donor_name,
                "source_pdb": str(source_pdb),
                "trimmed_pdb": str(trimmed_pdb),
                "length": str(len(trimmed_seq)),
                "trimmed_tail_removed": "".join(
                    donor_seq[i - 1] for i in range(1, len(donor_seq) + 1) if i not in set(keep_donor_indices)
                ),
                "trim_boundary_sf_position": str(args.final_length),
                "trim_boundary_donor_index": str(sf_to_donor[args.final_length]),
                "missing_sf_positions_1_to_trim": ";".join(map(str, missing_sf_positions)),
                "inherited_source_differences": ";".join(map(str, inherited_source_differences)),
                "index_to_sf": ";".join("" if sf_pos is None else str(sf_pos) for sf_pos in trimmed_index_to_sf),
                "sequence": trimmed_seq,
                **loop_strings,
                **loop_index_fields,
            }
        )
        fasta_records.append((backbone_id, trimmed_seq))
        for sf_pos in required_positions:
            donor_idx = sf_to_donor[sf_pos]
            pdb_key = sf_to_pdb_key[sf_pos]
            mapping_rows.append(
                {
                    "backbone_id": backbone_id,
                    "sfGFP_position": sf_pos,
                    "donor_sequence_index": donor_idx,
                    "pdb_chain": pdb_key[0].strip() or " ",
                    "pdb_resseq": pdb_key[1],
                    "pdb_icode": pdb_key[2],
                    "sfGFP_residue": sf_seq[sf_pos - 1],
                    "donor_residue": donor_seq[donor_idx - 1],
                }
            )

    fieldnames = [
        "backbone_id",
        "loop_rank",
        "donor_name",
        "source_pdb",
        "trimmed_pdb",
        "length",
        "trimmed_tail_removed",
        "trim_boundary_sf_position",
        "trim_boundary_donor_index",
        "missing_sf_positions_1_to_trim",
        "inherited_source_differences",
        "index_to_sf",
        "loop_129_147",
        "loop_129_147_seq_start",
        "loop_129_147_seq_end",
        "loop_188_198",
        "loop_188_198_seq_start",
        "loop_188_198_seq_end",
        "loop_209_216",
        "loop_209_216_seq_start",
        "loop_209_216_seq_end",
        "redesigned_residues_txt",
        "bias_json",
        "omit_json",
        "sequence",
    ]
    with (args.input_dir / "mapping_qc.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "backbone_id",
                "sfGFP_position",
                "donor_sequence_index",
                "pdb_chain",
                "pdb_resseq",
                "pdb_icode",
                "sfGFP_residue",
                "donor_residue",
            ],
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    selected_positions = beta_positions() - loop_positions() - CHROMOPHORE_POSITIONS
    selected_positions = {p for p in selected_positions if p <= args.final_length}
    with args.conservative_positions.open(newline="", encoding="utf-8") as handle:
        conservative_rows = list(csv.DictReader(handle))
    base_pssm_rows = [
        row
        for row in conservative_rows
        if int(row["position"]) in selected_positions
        and int(row["position"]) <= args.final_length
        and row["pdb_token"].strip()
    ]
    if len(base_pssm_rows) != args.expected_editable_count:
        raise SystemExit(
            f"Expected {args.expected_editable_count} barrel PSSM positions, found {len(base_pssm_rows)}"
        )
    if any(int(row["position"]) in loop_positions() or int(row["position"]) in CHROMOPHORE_POSITIONS for row in base_pssm_rows):
        raise SystemExit("Editable PSSM positions overlap fixed loops or chromophore")

    shutil.copy2(args.source_matrix, args.pssm_dir / "2b3p_sequence_variance_matrix.csv")
    shutil.copy2(args.source_position_tiers, args.pssm_dir / "2b3p_position_tiers.csv")
    source_bias = json.loads(args.source_bias.read_text(encoding="utf-8"))
    source_omit = json.loads(args.source_omit.read_text(encoding="utf-8"))
    pssm_out_rows: list[dict[str, str]] = []
    common_tokens: list[str] | None = None
    for row in metadata:
        backbone_id = row["backbone_id"]
        donor_seq = row["sequence"]
        index_to_sf = [int(x) if x else None for x in row["index_to_sf"].split(";")]
        sf_to_trimmed = {sf_pos: idx for idx, sf_pos in enumerate(index_to_sf, start=1) if sf_pos is not None}
        pdb_seq, pdb_residues = parse_pdb_residues(Path(row["trimmed_pdb"]))
        if pdb_seq != donor_seq:
            raise SystemExit(f"Prepared PDB sequence changed unexpectedly for {backbone_id}")
        mapped_tokens: list[str] = []
        mapped_bias: dict[str, dict[str, float]] = {}
        mapped_omit: dict[str, str] = {}
        for base in base_pssm_rows:
            sf_pos = int(base["position"])
            trimmed_idx = sf_to_trimmed.get(sf_pos)
            if trimmed_idx is None:
                raise SystemExit(f"Editable sfGFP position {sf_pos} is missing in {backbone_id}")
            chain, resseq, icode = pdb_residues[trimmed_idx - 1]
            target_token = f"{chain.strip() or 'A'}{resseq}{icode}"
            source_token = base["pdb_token"]
            mapped_tokens.append(target_token)
            mapped_bias[target_token] = source_bias[source_token]
            if source_token in source_omit:
                mapped_omit[target_token] = source_omit[source_token]
            pssm_out_rows.append(
                {
                    **base,
                    "backbone_id": backbone_id,
                    "sfGFP_position": str(sf_pos),
                    "source_pdb_token": source_token,
                    "target_sequence_index": str(trimmed_idx),
                    "target_pdb_token": target_token,
                }
            )
        if common_tokens is None:
            common_tokens = mapped_tokens
        elif common_tokens != mapped_tokens:
            print(f"WARNING: mapped editable tokens differ for {backbone_id}; per-backbone files will be used")
        residue_file = args.pssm_dir / f"{backbone_id}_redesigned_residues_barrel_pssm.txt"
        bias_file = args.pssm_dir / f"{backbone_id}_ligandmpnn_bias_AA_per_residue_barrel_pssm.json"
        omit_file = args.pssm_dir / f"{backbone_id}_ligandmpnn_omit_AA_per_residue_barrel_pssm.json"
        residue_file.write_text(" ".join(mapped_tokens) + "\n", encoding="utf-8")
        bias_file.write_text(json.dumps(mapped_bias, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        omit_file.write_text(json.dumps(mapped_omit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        row["redesigned_residues_txt"] = str(residue_file)
        row["bias_json"] = str(bias_file)
        row["omit_json"] = str(omit_file)

    meta_tsv = args.input_dir / "trimmed_backbones.tsv"
    with meta_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata)
    write_fasta(args.input_dir / "trimmed_backbones.fasta", fasta_records)
    with (args.input_dir / "input_backbones.txt").open("w", encoding="utf-8") as handle:
        for row in metadata:
            handle.write(f"{row['trimmed_pdb']}\n")

    with (args.pssm_dir / "barrel_pssm_positions.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "backbone_id",
                "sfGFP_position",
                "source_pdb_token",
                "target_sequence_index",
                "target_pdb_token",
                *conservative_rows[0].keys(),
            ],
        )
        writer.writeheader()
        writer.writerows(pssm_out_rows)
    tokens = common_tokens or []
    (args.pssm_dir / "redesigned_residues_barrel_pssm.txt").write_text(" ".join(tokens) + "\n", encoding="utf-8")
    (args.pssm_dir / "ligandmpnn_bias_AA_per_residue_barrel_pssm.json").write_text(
        (args.pssm_dir / f"{metadata[0]['backbone_id']}_ligandmpnn_bias_AA_per_residue_barrel_pssm.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (args.pssm_dir / "ligandmpnn_omit_AA_per_residue_barrel_pssm.json").write_text(
        (args.pssm_dir / f"{metadata[0]['backbone_id']}_ligandmpnn_omit_AA_per_residue_barrel_pssm.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    print(f"Prepared trimmed backbones: {len(metadata)}")
    print(f"Trim boundary sfGFP position: {args.final_length}")
    print(f"Trimmed lengths: {', '.join(row['length'] for row in metadata)}")
    print(f"Barrel PSSM editable residues: {len(tokens)}")
    print(f"Wrote {meta_tsv}")


def merge_mpnn(args: argparse.Namespace) -> None:
    fa_files = sorted(args.mpnn_out.rglob("*.fa"))
    total = 0
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fasta.open("w", encoding="utf-8") as out:
        for fa in fa_files:
            parts = fa.relative_to(args.mpnn_out).parts
            backbone = parts[0] if len(parts) > 0 else fa.stem
            temp_tag = parts[1] if len(parts) > 1 else "Tunknown"
            for i, (header, seq) in enumerate(read_fasta_records(fa), start=1):
                out.write(f">{backbone}_{temp_tag}__seq{i} {header}\n{wrap(seq)}\n")
                total += 1
    print(f"Merged {total} sequences from {len(fa_files)} files into {args.out_fasta}")
    if total == 0:
        raise SystemExit("No LigandMPNN FASTA records merged")


def load_trimmed_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["backbone_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def prefilter(args: argparse.Namespace) -> None:
    records = read_fasta_records(args.raw_fasta)
    metadata = load_trimmed_metadata(args.metadata_tsv)
    sf_seq = read_named_fasta(args.ref_fasta)["sfGFP"]
    exclusion = load_exclusion_sequences(args.exclusion_list)
    editable_by_backbone: dict[str, set[int]] = defaultdict(set)
    editable_sf_by_backbone: dict[str, set[int]] = defaultdict(set)
    with args.barrel_positions.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            editable_by_backbone[row["backbone_id"]].add(int(row["target_sequence_index"]))
            editable_sf_by_backbone[row["backbone_id"]].add(int(row["sfGFP_position"]))

    posterior: dict[tuple[int, str], float] = {}
    with args.matrix_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            posterior[(int(row["position"]), row["candidate_residue"])] = float(row["posterior_delta"])

    prepared: dict[str, dict[str, object]] = {}
    for backbone, row in metadata.items():
        index_to_sf = [int(x) if x else None for x in row["index_to_sf"].split(";")]
        sf_to_index = {sf_pos: idx for idx, sf_pos in enumerate(index_to_sf, start=1) if sf_pos is not None}
        missing_sf = set(position_list_csv(row["missing_sf_positions_1_to_trim"]))
        inherited_source_differences = set(position_list_csv(row["inherited_source_differences"]))
        loop_ranges = {
            f"loop_{start}_{end}": (
                int(row[f"loop_{start}_{end}_seq_start"]),
                int(row[f"loop_{start}_{end}_seq_end"]),
                row[f"loop_{start}_{end}"],
            )
            for start, end in LOOP_WINDOWS
        }
        tyg_indices = [sf_to_index[pos] for pos in (65, 66, 67)]
        prepared[backbone] = {
            "row": row,
            "index_to_sf": index_to_sf,
            "sf_to_index": sf_to_index,
            "missing_sf": missing_sf,
            "inherited_source_differences": inherited_source_differences,
            "loop_ranges": loop_ranges,
            "tyg_indices": tyg_indices,
            "editable_indices": editable_by_backbone.get(backbone, set()),
            "editable_sf_positions": editable_sf_by_backbone.get(backbone, set()),
        }

    seen: set[str] = set()
    rows: list[dict[str, object]] = []
    processed: list[tuple[str, str]] = []
    reason_counts: Counter[str] = Counter()
    for header, seq in records:
        seq = seq.strip().upper().replace(":", "")
        processed.append((header, seq))
        backbone = parse_backbone(header)
        source = prepared.get(backbone)
        reasons: list[str] = []
        if source is None:
            reasons.append("unknown_backbone")
            source_seq = ""
            expected_len = args.final_length
            editable_indices: set[int] = set()
            editable_sf_positions: set[int] = set()
            index_to_sf: list[int | None] = []
            source_differences: set[int] = set()
        else:
            source_row = source["row"]  # type: ignore[index]
            source_seq = str(source_row["sequence"])
            expected_len = int(source_row["length"])
            editable_indices = set(source["editable_indices"])  # type: ignore[arg-type]
            editable_sf_positions = set(source["editable_sf_positions"])  # type: ignore[arg-type]
            index_to_sf = list(source["index_to_sf"])  # type: ignore[arg-type]
            source_differences = set(source["inherited_source_differences"])  # type: ignore[arg-type]

        if len(seq) != expected_len:
            reasons.append("length")
        if not seq.startswith("M"):
            reasons.append("start_m")
        if set(seq) - AA20:
            reasons.append("nonstandard_aa")
        if source is not None and len(seq) == expected_len:
            tyg_indices = list(source["tyg_indices"])  # type: ignore[arg-type]
            if "".join(seq[i - 1] for i in tyg_indices) != "TYG":
                reasons.append("tyg65_67")
        elif len(seq) >= 67 and seq[64:67] != "TYG":
            reasons.append("tyg65_67")
        if seq in exclusion:
            reasons.append("exclusion_hit")
        if seq in seen:
            reasons.append("duplicate")
        seen.add(seq)

        changed_indices: set[int] = set()
        illegal_indices: list[int] = []
        mut_sf_positions: set[int] = set(source_differences)
        if source is not None and len(seq) == expected_len:
            loop_ranges = source["loop_ranges"]  # type: ignore[assignment]
            for loop_key, loop_data in loop_ranges.items():  # type: ignore[union-attr]
                loop_start, loop_end, loop_seq = loop_data
                if seq[loop_start - 1 : loop_end] != loop_seq:
                    reasons.append(f"fixed_{loop_key}_changed")
            changed_indices = {
                i for i, (source_aa, candidate_aa) in enumerate(zip(source_seq, seq), start=1) if source_aa != candidate_aa
            }
            illegal_indices = sorted(changed_indices - editable_indices)
            if illegal_indices:
                reasons.append("mutation_outside_allowed")
            for idx, candidate_aa in enumerate(seq, start=1):
                sf_pos = index_to_sf[idx - 1]
                if sf_pos is None or sf_pos > args.final_length:
                    continue
                if candidate_aa != sf_seq[sf_pos - 1]:
                    mut_sf_positions.add(sf_pos)
                else:
                    mut_sf_positions.discard(sf_pos)
        else:
            changed_indices = set()
            illegal_indices = []
            mut_sf_positions = set()

        lm_conf = parse_float(r"overall_confidence=([0-9.]+)", header)
        ligand_conf = parse_float(r"ligand_confidence=([0-9.]+)", header)
        seq_rec = parse_float(r"seq_rec=([0-9.]+)", header)
        if lm_conf is not None and lm_conf < args.lm_conf_min:
            reasons.append("mpnn_confidence")

        homo = longest_homopolymer(seq)
        dipep = longest_dipeptide_repeat(seq)
        lowc = low_complexity_fraction(seq)
        charge = net_charge(seq)
        cys_count = seq.count("C")
        if homo > args.homopolymer_max:
            reasons.append("homopolymer")
        if dipep > args.dipeptide_repeat_max:
            reasons.append("dipeptide_repeat")
        if lowc > args.low_complexity_max:
            reasons.append("low_complexity")
        if charge < args.net_charge_min or charge > args.net_charge_max:
            reasons.append("net_charge")
        if cys_count > args.cysteine_max:
            reasons.append("cysteine_count")

        editable_mutations = sorted(mut_sf_positions & editable_sf_positions)
        new_barrel_mutation_indices = sorted(changed_indices & editable_indices)
        new_barrel_pairs = [
            (idx, int(index_to_sf[idx - 1]))
            for idx in new_barrel_mutation_indices
            if index_to_sf and index_to_sf[idx - 1] is not None
        ]
        new_barrel_sf_positions = [sf_pos for _, sf_pos in new_barrel_pairs]
        deltas = [posterior.get((sf_pos, seq[idx - 1]), 0.0) for idx, sf_pos in new_barrel_pairs]
        pssm_sum = sum(deltas)
        pssm_mean = pssm_sum / len(deltas) if deltas else 0.0
        pssm_min = min(deltas) if deltas else 0.0

        lm = lm_conf if lm_conf is not None else 0.0
        ligand = ligand_conf if ligand_conf is not None else 0.0
        rec = seq_rec if seq_rec is not None else 0.0
        mutation_penalty = abs(len(mut_sf_positions) - args.mutation_target) / max(1, args.mutation_target)
        barrel_penalty = len(new_barrel_mutation_indices) / max(1, len(editable_indices))
        composite = 4.0 * lm + 1.5 * ligand + pssm_mean + 0.25 * rec - 0.5 * mutation_penalty - 0.4 * barrel_penalty

        passed = not reasons
        for reason in reasons or ["pass"]:
            reason_counts[reason] += 1
        rows.append(
            {
                "header": header,
                "backbone": backbone,
                "temperature": parse_temperature(header),
                "length": len(seq),
                "lm_confidence": lm_conf if lm_conf is not None else "",
                "ligand_confidence": ligand_conf if ligand_conf is not None else "",
                "seq_rec": seq_rec if seq_rec is not None else "",
                "mutation_count": len(mut_sf_positions),
                "editable_mutation_count": len(editable_mutations),
                "new_barrel_mutation_count": len(new_barrel_mutation_indices),
                "mutation_positions": ";".join(map(str, sorted(mut_sf_positions))),
                "editable_mutation_positions": ";".join(map(str, editable_mutations)),
                "new_barrel_sf_positions": ";".join(map(str, sorted(new_barrel_sf_positions))),
                "illegal_mutation_positions": ";".join(
                    map(
                        str,
                        [
                            index_to_sf[idx - 1] if index_to_sf and index_to_sf[idx - 1] is not None else f"idx{idx}"
                            for idx in illegal_indices
                        ],
                    )
                ),
                "illegal_target_indices": ";".join(map(str, illegal_indices)),
                "pssm_sum": pssm_sum,
                "pssm_mean": pssm_mean,
                "pssm_min": pssm_min,
                "homopolymer_max": homo,
                "dipeptide_repeat_max": dipep,
                "low_complexity_fraction": lowc,
                "net_charge": charge,
                "cysteine_count": cys_count,
                "composite_score": composite,
                "pass": passed,
                "reason": "|".join(reasons) if reasons else "pass",
                "sequence": seq,
                "loop_129_147_seq_start": source["loop_ranges"]["loop_129_147"][0] if source is not None else "",
                "loop_129_147_seq_end": source["loop_ranges"]["loop_129_147"][1] if source is not None else "",
                "loop_188_198_seq_start": source["loop_ranges"]["loop_188_198"][0] if source is not None else "",
                "loop_188_198_seq_end": source["loop_ranges"]["loop_188_198"][1] if source is not None else "",
                "loop_209_216_seq_start": source["loop_ranges"]["loop_209_216"][0] if source is not None else "",
                "loop_209_216_seq_end": source["loop_ranges"]["loop_209_216"][1] if source is not None else "",
                "_mut_set": set(mut_sf_positions),
            }
        )

    write_fasta(args.processed_fasta, processed)
    args.prefilter_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "header",
        "backbone",
        "temperature",
        "length",
        "lm_confidence",
        "ligand_confidence",
        "seq_rec",
        "mutation_count",
        "editable_mutation_count",
        "new_barrel_mutation_count",
        "mutation_positions",
        "editable_mutation_positions",
        "new_barrel_sf_positions",
        "illegal_mutation_positions",
        "illegal_target_indices",
        "pssm_sum",
        "pssm_mean",
        "pssm_min",
        "homopolymer_max",
        "dipeptide_repeat_max",
        "low_complexity_fraction",
        "net_charge",
        "cysteine_count",
        "composite_score",
        "pass",
        "reason",
        "sequence",
        "loop_129_147_seq_start",
        "loop_129_147_seq_end",
        "loop_188_198_seq_start",
        "loop_188_198_seq_end",
        "loop_209_216_seq_start",
        "loop_209_216_seq_end",
    ]
    with args.prefilter_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})

    passing = [row for row in rows if row["pass"]]
    passing.sort(key=lambda row: float(row["composite_score"]), reverse=True)
    by_lane: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in passing:
        by_lane[(str(row["backbone"]), str(row["temperature"]))].append(row)

    selected: list[dict[str, object]] = []
    selected_ids: set[int] = set()
    selected_sets: list[set[int]] = []

    def try_add(row: dict[str, object], require_diversity: bool = True) -> bool:
        ident = id(row)
        if ident in selected_ids:
            return False
        mut_set = row["_mut_set"]  # type: ignore[index]
        if require_diversity and any(jaccard(mut_set, old) > args.diversity_max_jaccard for old in selected_sets):
            return False
        selected.append(row)
        selected_ids.add(ident)
        selected_sets.append(mut_set)
        return True

    for lane in sorted(by_lane):
        kept = 0
        for row in by_lane[lane]:
            if kept >= args.boltz_per_lane:
                break
            if try_add(row, True):
                kept += 1
        if kept < args.boltz_per_lane:
            for row in by_lane[lane]:
                if kept >= args.boltz_per_lane:
                    break
                if try_add(row, False):
                    kept += 1
    for row in passing:
        if len(selected) >= args.boltz_target_total:
            break
        try_add(row, True)
    for row in passing:
        if len(selected) >= args.boltz_target_total:
            break
        try_add(row, False)
    selected = selected[: args.boltz_target_total]

    boltz_records: list[tuple[str, str]] = []
    manifest_rows: list[dict[str, object]] = []
    for i, row in enumerate(selected, start=1):
        temp_tag = str(row["temperature"]).replace(".", "p")
        donor_rank = metadata.get(str(row["backbone"]), {}).get("loop_rank", "X")
        name = (
            f"run09_candidate_{i:04d}_L{donor_rank}_T{temp_tag}"
            f"_score{float(row['composite_score']):.4f}_pssm{float(row['pssm_mean']):.4f}"
            f"_lm{float(row['lm_confidence'] or 0.0):.4f}_mut{int(row['mutation_count'])}"
        )
        boltz_records.append((name, str(row["sequence"])))
        manifest_rows.append({"candidate_id": name, **{key: row[key] for key in fieldnames}})
    write_fasta(args.boltz_fasta, boltz_records)
    with args.boltz_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", *fieldnames])
        writer.writeheader()
        writer.writerows(manifest_rows)

    by_temp_raw = Counter(parse_temperature(header) for header, _ in records)
    by_lane_pass = Counter((str(r["backbone"]), str(r["temperature"])) for r in passing)
    by_lane_selected = Counter((str(r["backbone"]), str(r["temperature"])) for r in selected)
    with args.summary_md.open("w", encoding="utf-8") as handle:
        handle.write("# Run09 Sequence Prefilter Summary\n\n")
        handle.write(f"- Raw records: {len(records)}\n")
        handle.write(f"- Processed records: {len(processed)}\n")
        handle.write(f"- Passing records: {len(passing)}\n")
        handle.write(f"- Boltz selected: {len(selected)}\n")
        handle.write(f"- Boltz target: {args.boltz_target_total}\n")
        handle.write(f"- Boltz per lane target: {args.boltz_per_lane}\n")
        editable_counts = sorted({len(v) for v in editable_by_backbone.values()})
        handle.write(f"- Editable positions per backbone: {', '.join(map(str, editable_counts))}\n\n")
        handle.write("## Raw By Temperature\n\n")
        for key in sorted(by_temp_raw):
            handle.write(f"- {key}: {by_temp_raw[key]}\n")
        handle.write("\n## Pass By Backbone And Temperature\n\n")
        for key in sorted(by_lane_pass):
            handle.write(f"- {key[0]} / {key[1]}: {by_lane_pass[key]}\n")
        handle.write("\n## Selected By Backbone And Temperature\n\n")
        for key in sorted(by_lane_selected):
            handle.write(f"- {key[0]} / {key[1]}: {by_lane_selected[key]}\n")
        handle.write("\n## Reason Counts\n\n")
        for key, value in reason_counts.most_common():
            handle.write(f"- {key}: {value}\n")

    print(f"Raw records: {len(records)}")
    print(f"Passing prefilter: {len(passing)}")
    print(f"Boltz selected: {len(selected)}")
    print(f"Wrote {args.prefilter_csv}")
    print(f"Wrote {args.boltz_fasta}")
    if len(selected) == 0:
        raise SystemExit("No sequences selected for Boltz")
    if len(selected) < args.boltz_target_total:
        print(f"WARNING: selected {len(selected)} below target {args.boltz_target_total}")


def make_yamls(args: argparse.Namespace) -> None:
    args.yaml_dir.mkdir(parents=True, exist_ok=True)
    args.chunk_dir.mkdir(parents=True, exist_ok=True)
    args.smoke_dir.mkdir(parents=True, exist_ok=True)
    for old in args.yaml_dir.glob("*.yaml"):
        old.unlink()
    for old in args.chunk_dir.glob("chunk_*"):
        if old.is_dir():
            shutil.rmtree(old)
    for old in args.smoke_dir.glob("*.yaml"):
        old.unlink()

    yaml_paths: list[Path] = []
    for header, seq in read_fasta_records(args.fasta):
        name = slug(header.split()[0], 180)
        path = args.yaml_dir / f"{name}.yaml"
        path.write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {seq}\n"
            "      msa: empty\n",
            encoding="utf-8",
        )
        yaml_paths.append(path)

    for i, path in enumerate(yaml_paths):
        chunk = args.chunk_dir / f"chunk_{i // args.chunk_size + 1:04d}"
        chunk.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, chunk / path.name)
        if i < args.smoke_n:
            shutil.copy2(path, args.smoke_dir / path.name)
    print(f"Generated YAMLs: {len(yaml_paths)}")
    print(f"Generated chunks: {math.ceil(len(yaml_paths) / args.chunk_size) if yaml_paths else 0}")
    print(f"Generated smoke YAMLs: {min(args.smoke_n, len(yaml_paths))}")
    if not yaml_paths:
        raise SystemExit("No YAMLs generated")


def run_tmalign(model: Path, ref_pdb: Path, tmalign_bin: Path) -> tuple[float, float, float, float]:
    result = subprocess.run([str(tmalign_bin), str(model), str(ref_pdb)], capture_output=True, text=True)
    tm1 = tm2 = rmsd = aligned = 0.0
    for line in result.stdout.splitlines():
        if "TM-score=" in line and "normalized by length of Chain_1" in line:
            tm1 = float(line.split("TM-score=")[1].split()[0])
        elif "TM-score=" in line and "normalized by length of Chain_2" in line:
            tm2 = float(line.split("TM-score=")[1].split()[0])
        elif "Aligned length=" in line:
            aligned = float(line.split("Aligned length=")[1].split(",")[0].strip())
            rmsd = float(line.split("RMSD=")[1].split(",")[0].strip())
    return tm1, tm2, rmsd, aligned


def plddt_window_values(plddt: np.ndarray, start: int, end: int) -> list[float]:
    return [float(plddt[pos - 1]) for pos in range(start, end + 1) if 0 <= pos - 1 < len(plddt)]


def analyze(args: argparse.Namespace) -> None:
    seqs = read_fasta_dict(args.fasta)
    manifest: dict[str, dict[str, str]] = {}
    with args.candidate_manifest.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            manifest[row["candidate_id"]] = row

    pred_dirs: list[Path] = []
    for pred_root in args.boltz_out.rglob("predictions"):
        pred_dirs.extend([p for p in pred_root.iterdir() if p.is_dir()])
    pred_dirs = sorted(set(pred_dirs))

    rows: list[dict[str, object]] = []
    for pred_dir in pred_dirs:
        model_files = sorted(pred_dir.glob("*_model_0.pdb")) or sorted(pred_dir.glob("*_model_0.cif"))
        if not model_files:
            continue
        conf_files = sorted(pred_dir.glob("confidence_*_model_0.json"))
        conf = json.loads(conf_files[0].read_text(encoding="utf-8")) if conf_files else {}
        plddt_files = sorted(pred_dir.glob("plddt_*_model_0.npz"))
        plddt_arr = np.asarray([], dtype=float)
        if plddt_files:
            plddt_arr = np.asarray(np.load(plddt_files[0])["plddt"], dtype=float)
        tm1, tm2, rmsd, aligned = run_tmalign(model_files[0], args.ref_pdb, args.tmalign_bin)
        source = manifest.get(pred_dir.name, {})
        whole_plddt = float(plddt_arr.mean()) if len(plddt_arr) else float(conf.get("complex_plddt", 0.0))
        mapped_loop_windows: list[tuple[int, int, int, int]] = []
        for start, end in LOOP_WINDOWS:
            seq_start = int(source.get(f"loop_{start}_{end}_seq_start") or start)
            seq_end = int(source.get(f"loop_{start}_{end}_seq_end") or end)
            mapped_loop_windows.append((start, end, seq_start, seq_end))
        loop_vals = [
            v
            for _, _, seq_start, seq_end in mapped_loop_windows
            for v in plddt_window_values(plddt_arr, seq_start, seq_end)
        ]
        fold_gate = (
            whole_plddt >= args.plddt_gate
            and tm2 >= args.tm2_gate
            and aligned >= args.aligned_gate
            and rmsd <= args.rmsd_gate
        )
        row: dict[str, object] = {
            "name": pred_dir.name,
            "confidence": conf.get("confidence_score", 0.0),
            "ptm": conf.get("ptm", 0.0),
            "plddt": whole_plddt,
            "tm1_audit": tm1,
            "tm2_audit": tm2,
            "rmsd_audit": rmsd,
            "aligned_audit": aligned,
            "fold_gate": fold_gate,
            "loop_plddt_mean": mean(loop_vals) if loop_vals else 0.0,
            "loop_plddt_min": min(loop_vals) if loop_vals else 0.0,
            "backbone": source.get("backbone", ""),
            "temperature": source.get("temperature", ""),
            "mutation_count": source.get("mutation_count", ""),
            "editable_mutation_count": source.get("editable_mutation_count", ""),
            "new_barrel_mutation_count": source.get("new_barrel_mutation_count", ""),
            "pssm_mean": source.get("pssm_mean", ""),
            "model": str(model_files[0]),
            "sequence": seqs.get(pred_dir.name, source.get("sequence", "")),
        }
        for start, end, seq_start, seq_end in mapped_loop_windows:
            vals = plddt_window_values(plddt_arr, seq_start, seq_end)
            row[f"loop_{start}_{end}_mean"] = mean(vals) if vals else 0.0
            row[f"loop_{start}_{end}_min"] = min(vals) if vals else 0.0
            row[f"loop_{start}_{end}_seq_start"] = seq_start
            row[f"loop_{start}_{end}_seq_end"] = seq_end
        rows.append(row)

    rows.sort(key=lambda r: (float(r["plddt"]), float(r["confidence"]), float(r["ptm"]), float(r["tm2_audit"])), reverse=True)
    fieldnames = [
        "name",
        "confidence",
        "ptm",
        "plddt",
        "tm1_audit",
        "tm2_audit",
        "rmsd_audit",
        "aligned_audit",
        "fold_gate",
        "loop_plddt_mean",
        "loop_plddt_min",
        "loop_129_147_mean",
        "loop_129_147_min",
        "loop_188_198_mean",
        "loop_188_198_min",
        "loop_209_216_mean",
        "loop_209_216_min",
        "loop_129_147_seq_start",
        "loop_129_147_seq_end",
        "loop_188_198_seq_start",
        "loop_188_198_seq_end",
        "loop_209_216_seq_start",
        "loop_209_216_seq_end",
        "backbone",
        "temperature",
        "mutation_count",
        "editable_mutation_count",
        "new_barrel_mutation_count",
        "pssm_mean",
        "model",
        "sequence",
    ]
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    def write_ranked_fasta(path: Path, selected: list[dict[str, object]]) -> int:
        with path.open("w", encoding="utf-8") as handle:
            written = 0
            for row in selected:
                seq = str(row["sequence"])
                if not seq:
                    continue
                handle.write(
                    f">{row['name']} plddt={float(row['plddt']):.4f} "
                    f"conf={float(row['confidence']):.4f} tm2={float(row['tm2_audit']):.4f} "
                    f"fold_gate={row['fold_gate']}\n{wrap(seq)}\n"
                )
                written += 1
        return written

    fold_rows = [row for row in rows if row["fold_gate"]]
    top_written = write_ranked_fasta(args.top_fasta, rows[: args.top_n])
    fold_written = write_ranked_fasta(args.fold_fasta, fold_rows[: args.top_n])
    print(f"Predictions analyzed: {len(rows)}")
    print(f"pLDDT >= {args.plddt_gate}: {sum(float(r['plddt']) >= args.plddt_gate for r in rows)}")
    print(f"pLDDT >= {args.plddt_preferred}: {sum(float(r['plddt']) >= args.plddt_preferred for r in rows)}")
    print(f"TM2 >= {args.tm2_gate}: {sum(float(r['tm2_audit']) >= args.tm2_gate for r in rows)}")
    print(f"Fold gate pass: {len(fold_rows)}")
    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.top_fasta} records={top_written}")
    print(f"Wrote {args.fold_fasta} records={fold_written}")
    if not rows:
        raise SystemExit("No predictions analyzed")


def package(args: argparse.Namespace) -> None:
    warnings.simplefilter("ignore", BiopythonWarning)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out_dir / "pdbs"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for old in pdb_dir.glob("rank*_*.pdb"):
        old.unlink()

    with args.analysis_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))[: args.top_n]
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    packaged: list[dict[str, str]] = []
    for rank, row in enumerate(rows, 1):
        source = Path(row["model"])
        pdb_path = pdb_dir / f"rank{rank:03d}_{slug(row['name'], 120)}.pdb"
        if source.suffix.lower() == ".pdb":
            shutil.copy2(source, pdb_path)
        else:
            structure = parser.get_structure(pdb_path.stem, str(source))
            io.set_structure(structure)
            io.save(str(pdb_path))
        packaged.append({**row, "rank": str(rank), "pdb_file": str(pdb_path)})

    fieldnames = [
        "rank",
        "name",
        "plddt",
        "confidence",
        "ptm",
        "tm2_audit",
        "rmsd_audit",
        "aligned_audit",
        "fold_gate",
        "loop_plddt_mean",
        "backbone",
        "temperature",
        "mutation_count",
        "new_barrel_mutation_count",
        "model",
        "pdb_file",
        "sequence",
    ]
    with (args.out_dir / "manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in packaged:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    with (args.out_dir / "top_ranked_by_plddt.fasta").open("w", encoding="utf-8") as handle:
        for row in packaged:
            handle.write(
                f">rank{int(row['rank']):03d}|{row['name']}|pdb={Path(row['pdb_file']).name}|"
                f"plddt={float(row['plddt']):.4f}|conf={float(row['confidence']):.4f}|"
                f"tm2={float(row['tm2_audit']):.4f}|fold_gate={row['fold_gate']}\n"
                f"{wrap(row['sequence'])}\n"
            )
    with (args.out_dir / "load_top_by_plddt.pml").open("w", encoding="utf-8") as handle:
        handle.write("reinitialize\n")
        for row in packaged:
            name = Path(row["pdb_file"]).name
            obj = Path(name).stem.replace(".", "_").replace("-", "_")
            handle.write(f"load pdbs/{name}, {obj}\n")
        handle.write("hide everything\nshow cartoon\nspectrum count, rainbow\nzoom\n")
    print(f"Packaged {len(packaged)} structures into {args.out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare-inputs")
    p.add_argument("--root-dir", type=Path, required=True)
    p.add_argument("--loop-csv", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path, required=True)
    p.add_argument("--conservative-positions", type=Path, required=True)
    p.add_argument("--source-matrix", type=Path, required=True)
    p.add_argument("--source-position-tiers", type=Path, required=True)
    p.add_argument("--source-bias", type=Path, required=True)
    p.add_argument("--source-omit", type=Path, required=True)
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--pssm-dir", type=Path, required=True)
    p.add_argument("--top-loop-n", type=int, default=2)
    p.add_argument("--final-length", type=int, default=230)
    p.add_argument("--expected-editable-count", type=int, default=52)
    p.set_defaults(func=prepare_inputs)

    p = sub.add_parser("merge-mpnn")
    p.add_argument("--mpnn-out", type=Path, required=True)
    p.add_argument("--out-fasta", type=Path, required=True)
    p.set_defaults(func=merge_mpnn)

    p = sub.add_parser("prefilter")
    p.add_argument("--raw-fasta", type=Path, required=True)
    p.add_argument("--processed-fasta", type=Path, required=True)
    p.add_argument("--prefilter-csv", type=Path, required=True)
    p.add_argument("--summary-md", type=Path, required=True)
    p.add_argument("--boltz-fasta", type=Path, required=True)
    p.add_argument("--boltz-manifest", type=Path, required=True)
    p.add_argument("--metadata-tsv", type=Path, required=True)
    p.add_argument("--barrel-positions", type=Path, required=True)
    p.add_argument("--matrix-csv", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path, required=True)
    p.add_argument("--exclusion-list", type=Path, required=True)
    p.add_argument("--final-length", type=int, default=230)
    p.add_argument("--lm-conf-min", type=float, default=0.43)
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--dipeptide-repeat-max", type=int, default=3)
    p.add_argument("--low-complexity-max", type=float, default=0.05)
    p.add_argument("--net-charge-min", type=int, default=-20)
    p.add_argument("--net-charge-max", type=int, default=20)
    p.add_argument("--cysteine-max", type=int, default=3)
    p.add_argument("--boltz-target-total", type=int, default=2000)
    p.add_argument("--boltz-per-lane", type=int, default=200)
    p.add_argument("--mutation-target", type=int, default=40)
    p.add_argument("--diversity-max-jaccard", type=float, default=0.90)
    p.set_defaults(func=prefilter)

    p = sub.add_parser("make-yamls")
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--yaml-dir", type=Path, required=True)
    p.add_argument("--chunk-dir", type=Path, required=True)
    p.add_argument("--smoke-dir", type=Path, required=True)
    p.add_argument("--chunk-size", type=int, default=50)
    p.add_argument("--smoke-n", type=int, default=10)
    p.set_defaults(func=make_yamls)

    p = sub.add_parser("analyze")
    p.add_argument("--boltz-out", type=Path, required=True)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--candidate-manifest", type=Path, required=True)
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--top-fasta", type=Path, required=True)
    p.add_argument("--fold-fasta", type=Path, required=True)
    p.add_argument("--ref-pdb", type=Path, required=True)
    p.add_argument("--tmalign-bin", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=64)
    p.add_argument("--plddt-gate", type=float, default=0.85)
    p.add_argument("--plddt-preferred", type=float, default=0.90)
    p.add_argument("--tm2-gate", type=float, default=0.75)
    p.add_argument("--aligned-gate", type=float, default=210)
    p.add_argument("--rmsd-gate", type=float, default=3.5)
    p.set_defaults(func=analyze)

    p = sub.add_parser("package")
    p.add_argument("--analysis-csv", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=64)
    p.set_defaults(func=package)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
