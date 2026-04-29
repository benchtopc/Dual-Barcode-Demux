#!/usr/bin/env python3
# bin/demux_scm.py

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

M13_FWD = "GTAAACGACGGCCAGT"
M13_REV = "AGCGGATAACAATTTCACACAGGA"
SPACER = "GAT"
RC_SPACER = "ATC"


@dataclass(frozen=True)
class BarcodeCall:
    barcode_id: str
    barcode_seq: str
    source: str
    total_score: int
    linker_mismatches: int
    spacer_mismatches: int
    observed_barcode: str
    position: int


@dataclass(frozen=True)
class Assignment:
    read_id: str
    status: str
    orientation: str
    row_barcode_id: str
    row_barcode_seq: str
    row_source: str
    row_score: str
    column_barcode_id: str
    column_barcode_seq: str
    column_source: str
    column_score: str
    well_id: str
    sample_id: str
    alias: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Demultiplex reads with the SCM workflow implementing the published 96-well combinatorial row/column barcode design from FASTQ."
    )
    parser.add_argument("--fastq", required=True, help="Input FASTQ or FASTQ.GZ.")
    parser.add_argument("--plate-map", required=True, help="Plate map CSV.")
    parser.add_argument("--row-barcodes", required=True, help="Row barcode CSV.")
    parser.add_argument("--column-barcodes", required=True, help="Column barcode CSV.")
    parser.add_argument("--front-window", type=int, default=200, help="5' search window size.")
    parser.add_argument("--rear-window", type=int, default=200, help="3' search window size.")
    parser.add_argument(
        "--max-linker-mismatches",
        type=int,
        default=2,
        help="Maximum mismatches allowed in the linker match.",
    )
    parser.add_argument(
        "--max-total-score",
        type=int,
        default=4,
        help="Maximum combined mismatch score across spacer, linker, and barcode.",
    )
    parser.add_argument(
        "--min-margin",
        type=int,
        default=2,
        help="Minimum distance between the best and second-best barcode score.",
    )
    parser.add_argument(
        "--keep-partial",
        action="store_true",
        help="Track row-only and column-only reads in the report.",
    )
    parser.add_argument("--outdir", default=".", help="Output directory.")
    return parser.parse_args()


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def hamming(a: str, b: str) -> int:
    return sum(left != right for left, right in zip(a, b))


def sanitize_alias(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned or "sample"


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def iter_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with open_text_auto(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            if not quality:
                raise ValueError("Malformed FASTQ: truncated record encountered.")
            yield header.rstrip(), sequence.upper(), plus, quality


def read_barcode_table(path: Path, id_field: str, seq_field: str) -> Dict[str, str]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {id_field, seq_field}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        rows: Dict[str, str] = {}
        for row in reader:
            barcode_id = row[id_field].strip()
            barcode_seq = row[seq_field].strip().upper()
            if len(barcode_seq) != 8:
                raise ValueError(f"{barcode_id} in {path} does not have an 8 nt barcode.")
            if barcode_id in rows:
                raise ValueError(f"Duplicate barcode identifier in {path}: {barcode_id}")
            rows[barcode_id] = barcode_seq
    return rows


def read_plate_map(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "well_id",
            "row_id",
            "column_id",
            "row_barcode_id",
            "column_barcode_id",
            "sample_id",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

        plate_map: Dict[Tuple[str, str], Dict[str, str]] = {}
        aliases: set[str] = set()
        for row in reader:
            row_id = row["row_barcode_id"].strip()
            column_id = row["column_barcode_id"].strip()
            key = (row_id, column_id)
            if key in plate_map:
                raise ValueError(f"Duplicate barcode combination in plate map: {key}")
            well_id = row["well_id"].strip()
            sample_id = row["sample_id"].strip()
            alias = sanitize_alias(row.get("alias") or sample_id or well_id)
            if alias in aliases:
                raise ValueError(f"Alias must be unique. Duplicate alias found: {alias}")
            aliases.add(alias)
            plate_map[key] = {
                "plate_id": (row.get("plate_id") or "").strip(),
                "well_id": well_id,
                "row_id": row["row_id"].strip(),
                "column_id": row["column_id"].strip(),
                "row_barcode_id": row_id,
                "column_barcode_id": column_id,
                "sample_id": sample_id,
                "alias": alias,
                "notes": (row.get("notes") or "").strip(),
            }
    return plate_map


def candidate_calls_front(
    window: str,
    linker: str,
    barcodes: Dict[str, str],
    max_linker_mismatches: int,
) -> List[BarcodeCall]:
    calls: List[BarcodeCall] = []
    linker_len = len(linker)
    for position in range(len(window) - linker_len + 1):
        linker_slice = window[position : position + linker_len]
        linker_mismatches = hamming(linker_slice, linker)
        if linker_mismatches > max_linker_mismatches or position < 11:
            continue
        observed_spacer = window[position - 11 : position - 8]
        observed_barcode = window[position - 8 : position]
        spacer_mismatches = hamming(observed_spacer, SPACER)
        for barcode_id, barcode_seq in barcodes.items():
            total_score = (
                linker_mismatches
                + spacer_mismatches
                + hamming(observed_barcode, barcode_seq)
            )
            calls.append(
                BarcodeCall(
                    barcode_id=barcode_id,
                    barcode_seq=barcode_seq,
                    source="front",
                    total_score=total_score,
                    linker_mismatches=linker_mismatches,
                    spacer_mismatches=spacer_mismatches,
                    observed_barcode=observed_barcode,
                    position=position,
                )
            )
    return calls


def candidate_calls_rear(
    window: str,
    linker: str,
    barcodes: Dict[str, str],
    max_linker_mismatches: int,
) -> List[BarcodeCall]:
    calls: List[BarcodeCall] = []
    linker_rc = reverse_complement(linker)
    linker_len = len(linker_rc)
    max_start = len(window) - linker_len - 11
    if max_start < 0:
        return calls
    for position in range(max_start + 1):
        linker_slice = window[position : position + linker_len]
        linker_mismatches = hamming(linker_slice, linker_rc)
        if linker_mismatches > max_linker_mismatches:
            continue
        barcode_start = position + linker_len
        barcode_end = barcode_start + 8
        spacer_end = barcode_end + 3
        observed_barcode_rc = window[barcode_start:barcode_end]
        observed_barcode = reverse_complement(observed_barcode_rc)
        observed_spacer = window[barcode_end:spacer_end]
        spacer_mismatches = hamming(observed_spacer, RC_SPACER)
        for barcode_id, barcode_seq in barcodes.items():
            total_score = (
                linker_mismatches
                + spacer_mismatches
                + hamming(observed_barcode, barcode_seq)
            )
            calls.append(
                BarcodeCall(
                    barcode_id=barcode_id,
                    barcode_seq=barcode_seq,
                    source="rear",
                    total_score=total_score,
                    linker_mismatches=linker_mismatches,
                    spacer_mismatches=spacer_mismatches,
                    observed_barcode=observed_barcode,
                    position=position,
                )
            )
    return calls


def choose_barcode_call(
    front_window: str,
    rear_window: str,
    linker: str,
    barcodes: Dict[str, str],
    max_linker_mismatches: int,
    max_total_score: int,
    min_margin: int,
) -> Optional[BarcodeCall]:
    candidates = candidate_calls_front(
        window=front_window,
        linker=linker,
        barcodes=barcodes,
        max_linker_mismatches=max_linker_mismatches,
    )
    candidates.extend(
        candidate_calls_rear(
            window=rear_window,
            linker=linker,
            barcodes=barcodes,
            max_linker_mismatches=max_linker_mismatches,
        )
    )
    if not candidates:
        return None

    candidates.sort(
        key=lambda call: (
            call.total_score,
            call.barcode_id,
            call.source,
            call.position,
        )
    )
    best = candidates[0]
    second = next(
        (candidate for candidate in candidates if candidate.barcode_id != best.barcode_id),
        None,
    )
    if best.total_score > max_total_score:
        return None
    if second and second.total_score - best.total_score < min_margin:
        return None
    return best


def determine_orientation(
    row_call: Optional[BarcodeCall],
    column_call: Optional[BarcodeCall],
) -> str:
    if row_call and column_call:
        if row_call.source == "front" and column_call.source == "rear":
            return "forward_strand"
        if row_call.source == "rear" and column_call.source == "front":
            return "reverse_strand"
        return "mixed"
    return "unknown"


def write_validated_plate_map(
    plate_map: Dict[Tuple[str, str], Dict[str, str]],
    out_path: Path,
) -> None:
    rows = list(plate_map.values())
    fieldnames = [
        "plate_id",
        "well_id",
        "row_id",
        "column_id",
        "row_barcode_id",
        "column_barcode_id",
        "sample_id",
        "alias",
        "notes",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: item["well_id"]):
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def render_report(
    report_path: Path,
    fastq_name: str,
    total_reads: int,
    status_counts: Counter,
    sample_counts: Counter,
    partial_counts: Counter,
    assignments_preview: List[Assignment],
) -> None:
    assignment_rate = 0.0
    if total_reads:
        assignment_rate = 100.0 * status_counts.get("assigned", 0) / total_reads

    sample_rows = "\n".join(
        f"<tr><td>{html.escape(alias)}</td><td>{count}</td></tr>"
        for alias, count in sorted(sample_counts.items())
    ) or "<tr><td colspan='2'>No assigned samples</td></tr>"

    status_rows = "\n".join(
        f"<tr><td>{html.escape(status)}</td><td>{count}</td></tr>"
        for status, count in sorted(status_counts.items())
    ) or "<tr><td colspan='2'>No reads</td></tr>"

    partial_rows = "\n".join(
        f"<tr><td>{html.escape(status)}</td><td>{count}</td></tr>"
        for status, count in sorted(partial_counts.items())
    ) or "<tr><td colspan='2'>No partial assignments recorded</td></tr>"

    preview_rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(item.read_id)}</td>"
        f"<td>{html.escape(item.status)}</td>"
        f"<td>{html.escape(item.orientation)}</td>"
        f"<td>{html.escape(item.row_barcode_id)}</td>"
        f"<td>{html.escape(item.column_barcode_id)}</td>"
        f"<td>{html.escape(item.well_id)}</td>"
        f"<td>{html.escape(item.alias)}</td>"
        "</tr>"
        for item in assignments_preview[:25]
    ) or "<tr><td colspan='7'>No assignments to preview</td></tr>"

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SCM demultiplex summary</title>
<style>
body {{
  font-family: Arial, sans-serif;
  margin: 2rem;
  color: #222;
}}
h1, h2 {{
  margin-bottom: 0.5rem;
}}
.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 1rem;
}}
.card {{
  border: 1px solid #ddd;
  border-radius: 12px;
  padding: 1rem;
}}
.metric {{
  font-size: 2rem;
  font-weight: 700;
}}
table {{
  border-collapse: collapse;
  width: 100%;
}}
th, td {{
  border: 1px solid #ddd;
  padding: 0.4rem 0.6rem;
  text-align: left;
}}
th {{
  background: #f5f5f5;
}}
code {{
  background: #f5f5f5;
  padding: 0.1rem 0.3rem;
  border-radius: 4px;
}}
</style>
</head>
<body>
<h1>SCM 96-well demultiplex summary</h1>
<p><strong>Input:</strong> <code>{html.escape(fastq_name)}</code></p>

<div class="grid">
  <div class="card">
    <h2>Total reads</h2>
    <div class="metric">{total_reads}</div>
  </div>
  <div class="card">
    <h2>Assigned reads</h2>
    <div class="metric">{status_counts.get("assigned", 0)}</div>
  </div>
  <div class="card">
    <h2>Assignment rate</h2>
    <div class="metric">{assignment_rate:.2f}%</div>
  </div>
</div>

<div class="grid">
  <div class="card">
    <h2>Status counts</h2>
    <table>
      <thead><tr><th>Status</th><th>Reads</th></tr></thead>
      <tbody>{status_rows}</tbody>
    </table>
  </div>
  <div class="card">
    <h2>Assigned reads per sample</h2>
    <table>
      <thead><tr><th>Alias</th><th>Reads</th></tr></thead>
      <tbody>{sample_rows}</tbody>
    </table>
  </div>
  <div class="card">
    <h2>Partial assignment counts</h2>
    <table>
      <thead><tr><th>Status</th><th>Reads</th></tr></thead>
      <tbody>{partial_rows}</tbody>
    </table>
  </div>
</div>

<div class="card">
  <h2>Assignment preview</h2>
  <table>
    <thead>
      <tr>
        <th>Read ID</th>
        <th>Status</th>
        <th>Orientation</th>
        <th>Row barcode</th>
        <th>Column barcode</th>
        <th>Well</th>
        <th>Alias</th>
      </tr>
    </thead>
    <tbody>{preview_rows}</tbody>
  </table>
</div>
</body>
</html>
"""
    report_path.write_text(document, encoding="utf-8")


def main() -> None:
    args = parse_args()

    fastq_path = Path(args.fastq)
    plate_map_path = Path(args.plate_map)
    row_barcodes_path = Path(args.row_barcodes)
    column_barcodes_path = Path(args.column_barcodes)
    outdir = Path(args.outdir)

    demux_dir = outdir / "demultiplexed_fastq"
    tables_dir = outdir / "tables"
    reports_dir = outdir / "reports"
    metadata_dir = outdir / "metadata"

    for directory in (demux_dir, tables_dir, reports_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)

    row_barcodes = read_barcode_table(
        path=row_barcodes_path,
        id_field="row_barcode_id",
        seq_field="barcode_seq",
    )
    column_barcodes = read_barcode_table(
        path=column_barcodes_path,
        id_field="column_barcode_id",
        seq_field="barcode_seq",
    )
    plate_map = read_plate_map(plate_map_path)

    missing_rows = sorted({row["row_barcode_id"] for row in plate_map.values()} - set(row_barcodes))
    missing_columns = sorted(
        {row["column_barcode_id"] for row in plate_map.values()} - set(column_barcodes)
    )
    if missing_rows:
        raise ValueError(f"Plate map references unknown row barcode IDs: {', '.join(missing_rows)}")
    if missing_columns:
        raise ValueError(
            f"Plate map references unknown column barcode IDs: {', '.join(missing_columns)}"
        )

    write_validated_plate_map(
        plate_map=plate_map,
        out_path=metadata_dir / "plate_map.validated.csv",
    )

    output_handles: Dict[str, gzip.GzipFile] = {}
    status_counts: Counter = Counter()
    partial_counts: Counter = Counter()
    sample_counts: Counter = Counter()
    assignment_preview: List[Assignment] = []
    total_reads = 0

    assignments_path = tables_dir / "read_assignments.tsv.gz"
    with gzip.open(assignments_path, "wt") as assignments_handle:
        assignment_header = [
            "read_id",
            "status",
            "orientation",
            "row_barcode_id",
            "row_barcode_seq",
            "row_source",
            "row_score",
            "column_barcode_id",
            "column_barcode_seq",
            "column_source",
            "column_score",
            "well_id",
            "sample_id",
            "alias",
        ]
        assignments_handle.write("\t".join(assignment_header) + "\n")

        for header, sequence, plus, quality in iter_fastq(fastq_path):
            total_reads += 1
            read_id = header.split()[0].removeprefix("@")
            front_window = sequence[: args.front_window]
            rear_window = sequence[-args.rear_window :]

            row_call = choose_barcode_call(
                front_window=front_window,
                rear_window=rear_window,
                linker=M13_FWD,
                barcodes=row_barcodes,
                max_linker_mismatches=args.max_linker_mismatches,
                max_total_score=args.max_total_score,
                min_margin=args.min_margin,
            )
            column_call = choose_barcode_call(
                front_window=front_window,
                rear_window=rear_window,
                linker=M13_REV,
                barcodes=column_barcodes,
                max_linker_mismatches=args.max_linker_mismatches,
                max_total_score=args.max_total_score,
                min_margin=args.min_margin,
            )

            orientation = determine_orientation(row_call=row_call, column_call=column_call)

            if row_call and column_call:
                mapping = plate_map.get((row_call.barcode_id, column_call.barcode_id))
                if mapping:
                    status = "assigned"
                    alias = mapping["alias"]
                    well_id = mapping["well_id"]
                    sample_id = mapping["sample_id"]
                else:
                    status = "barcode_pair_not_in_plate_map"
                    alias = "unclassified"
                    well_id = ""
                    sample_id = ""
            elif row_call and not column_call:
                status = "row_only"
                alias = "unclassified"
                well_id = ""
                sample_id = ""
                if args.keep_partial:
                    partial_counts[status] += 1
            elif column_call and not row_call:
                status = "column_only"
                alias = "unclassified"
                well_id = ""
                sample_id = ""
                if args.keep_partial:
                    partial_counts[status] += 1
            else:
                status = "unclassified"
                alias = "unclassified"
                well_id = ""
                sample_id = ""

            status_counts[status] += 1
            if status == "assigned":
                sample_counts[alias] += 1

            assignment = Assignment(
                read_id=read_id,
                status=status,
                orientation=orientation,
                row_barcode_id=row_call.barcode_id if row_call else "",
                row_barcode_seq=row_call.barcode_seq if row_call else "",
                row_source=row_call.source if row_call else "",
                row_score=str(row_call.total_score) if row_call else "",
                column_barcode_id=column_call.barcode_id if column_call else "",
                column_barcode_seq=column_call.barcode_seq if column_call else "",
                column_source=column_call.source if column_call else "",
                column_score=str(column_call.total_score) if column_call else "",
                well_id=well_id,
                sample_id=sample_id,
                alias=alias,
            )
            if len(assignment_preview) < 250:
                assignment_preview.append(assignment)
            assignments_handle.write(
                "\t".join(
                    [
                        assignment.read_id,
                        assignment.status,
                        assignment.orientation,
                        assignment.row_barcode_id,
                        assignment.row_barcode_seq,
                        assignment.row_source,
                        assignment.row_score,
                        assignment.column_barcode_id,
                        assignment.column_barcode_seq,
                        assignment.column_source,
                        assignment.column_score,
                        assignment.well_id,
                        assignment.sample_id,
                        assignment.alias,
                    ]
                )
                + "\n"
            )

            if alias not in output_handles:
                output_dir = demux_dir / alias
                output_dir.mkdir(parents=True, exist_ok=True)
                output_handles[alias] = gzip.open(output_dir / "reads.fastq.gz", "wt")
            output_handles[alias].write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

    for handle in output_handles.values():
        handle.close()

    with (tables_dir / "status_counts.tsv").open("w") as handle:
        handle.write("status\treads\n")
        for status, count in sorted(status_counts.items()):
            handle.write(f"{status}\t{count}\n")

    with (tables_dir / "sample_counts.tsv").open("w") as handle:
        handle.write("alias\treads\n")
        for alias, count in sorted(sample_counts.items()):
            handle.write(f"{alias}\t{count}\n")

    with (tables_dir / "run_metadata.json").open("w") as handle:
        json.dump(
            {
                "fastq": str(fastq_path),
                "plate_map": str(plate_map_path),
                "row_barcodes": str(row_barcodes_path),
                "column_barcodes": str(column_barcodes_path),
                "front_window": args.front_window,
                "rear_window": args.rear_window,
                "max_linker_mismatches": args.max_linker_mismatches,
                "max_total_score": args.max_total_score,
                "min_margin": args.min_margin,
                "keep_partial": args.keep_partial,
                "total_reads": total_reads,
                "status_counts": dict(status_counts),
                "sample_counts": dict(sample_counts),
                "partial_counts": dict(partial_counts),
            },
            handle,
            indent=2,
        )

    render_report(
        report_path=reports_dir / "demux_summary.html",
        fastq_name=fastq_path.name,
        total_reads=total_reads,
        status_counts=status_counts,
        sample_counts=sample_counts,
        partial_counts=partial_counts,
        assignments_preview=assignment_preview,
    )


if __name__ == "__main__":
    main()
