#!/usr/bin/env python3
# bin/align_samples_scm.py

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align demultiplexed per-sample FASTQ files with minimap2 and summarize coverage/gaps."
    )
    parser.add_argument("--demux-dir", required=True, help="Directory containing per-sample FASTQ folders.")
    parser.add_argument("--reference", required=True, help="Reference FASTA.")
    parser.add_argument("--aligner-preset", default="map-ont", help="minimap2 preset.")
    parser.add_argument("--threads", type=int, default=4, help="Threads for minimap2/samtools.")
    parser.add_argument("--min-mapq", type=int, default=10, help="Minimum MAPQ retained in filtered BAM.")
    parser.add_argument("--min-depth-ko", type=int, default=10, help="Depth threshold for gap summaries.")
    parser.add_argument("--keep-unfiltered-bam", action="store_true", help="Keep unfiltered BAMs.")
    parser.add_argument("--outdir", default=".", help="Output directory.")
    return parser.parse_args()


def run_command(command: List[str] | str, *, shell: bool = False) -> None:
    printable = command if isinstance(command, str) else " ".join(command)
    print(f"[scm-align] running: {printable}", flush=True)
    if shell:
        subprocess.run(["bash", "-o", "pipefail", "-lc", command], check=True)
    else:
        subprocess.run(command, check=True)


def ensure_tool(name: str) -> None:
    if not shutil.which(name):
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def find_fastq_files(demux_dir: Path) -> List[Tuple[str, Path]]:
    pairs: List[Tuple[str, Path]] = []
    for sample_dir in sorted(demux_dir.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name == "unclassified":
            continue
        candidates = sorted(sample_dir.glob("*.fastq.gz"))
        candidates += sorted(sample_dir.glob("*.fq.gz"))
        candidates += sorted(sample_dir.glob("*.fastq"))
        candidates += sorted(sample_dir.glob("*.fq"))
        for candidate in candidates:
            if candidate.is_file():
                pairs.append((sample_dir.name, candidate))
                break
    return pairs


def count_fastq_reads(path: Path) -> int:
    with open_text_auto(path) as handle:
        return sum(1 for _ in handle) // 4


def parse_reference_lengths(reference: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    current: Optional[str] = None
    current_len = 0
    with reference.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = current_len
                current = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
    if current is not None:
        lengths[current] = current_len
    if not lengths:
        raise ValueError(f"No reference contigs found in {reference}")
    return lengths


def parse_depth_file(path: Path) -> Dict[Tuple[str, int], int]:
    depth: Dict[Tuple[str, int], int] = {}
    with path.open("r") as handle:
        for line in handle:
            if not line.strip():
                continue
            contig, pos, value = line.rstrip("\n").split("\t")[:3]
            depth[(contig, int(pos))] = int(value)
    return depth


def write_gap_summary(
    depth_path: Path,
    reference_lengths: Dict[str, int],
    min_depth: int,
    output_path: Path,
) -> None:
    depth = parse_depth_file(depth_path)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", "gap_start", "gap_end", "gap_length", "max_depth_in_gap", "min_depth_threshold"])

        for contig, length in reference_lengths.items():
            gap_start: Optional[int] = None
            gap_depths: List[int] = []
            for pos in range(1, length + 1):
                value = depth.get((contig, pos), 0)
                if value < min_depth:
                    if gap_start is None:
                        gap_start = pos
                        gap_depths = []
                    gap_depths.append(value)
                elif gap_start is not None:
                    writer.writerow([contig, gap_start, pos - 1, pos - gap_start, max(gap_depths), min_depth])
                    gap_start = None
                    gap_depths = []
            if gap_start is not None:
                writer.writerow([contig, gap_start, length, length - gap_start + 1, max(gap_depths), min_depth])


def parse_flagstat_mapped(path: Path) -> Tuple[str, str]:
    total = "0"
    mapped = "0"
    with path.open("r") as handle:
        for line in handle:
            if " in total " in line:
                total = line.split()[0]
            elif " mapped (" in line and "primary" not in line:
                mapped = line.split()[0]
                break
    return total, mapped


def render_alignment_summary(report_path: Path, rows: List[Dict[str, str]]) -> None:
    body_rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(row['sample'])}</td>"
        f"<td>{html.escape(row['fastq_reads'])}</td>"
        f"<td>{html.escape(row['mapped_reads'])}</td>"
        f"<td>{html.escape(row['bam'])}</td>"
        "</tr>"
        for row in rows
    )
    report_path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SCM alignment summary</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; }}
th {{ background: #f5f5f5; }}
</style>
</head>
<body>
<h1>SCM alignment summary</h1>
<table>
<thead><tr><th>Sample</th><th>FASTQ reads</th><th>Mapped reads</th><th>BAM</th></tr></thead>
<tbody>{body_rows or "<tr><td colspan='4'>No samples found</td></tr>"}</tbody>
</table>
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()

    for tool in ("minimap2", "samtools", "bash"):
        ensure_tool(tool)

    demux_dir = Path(args.demux_dir)
    reference = Path(args.reference)
    outdir = Path(args.outdir)

    alignments_dir = outdir / "alignments"
    gaps_dir = outdir / "gaps"
    reports_dir = outdir / "analysis_reports"
    for directory in (alignments_dir, gaps_dir, reports_dir):
        directory.mkdir(parents=True, exist_ok=True)

    reference_lengths = parse_reference_lengths(reference)
    run_command(["samtools", "faidx", str(reference)])

    sample_fastqs = find_fastq_files(demux_dir)
    if not sample_fastqs:
        raise RuntimeError(f"No per-sample FASTQ files found under {demux_dir}")

    manifest_rows: List[Dict[str, str]] = []
    for sample, fastq in sample_fastqs:
        sample_align_dir = alignments_dir / sample
        sample_gap_dir = gaps_dir / sample
        sample_align_dir.mkdir(parents=True, exist_ok=True)
        sample_gap_dir.mkdir(parents=True, exist_ok=True)

        unfiltered_bam = sample_align_dir / f"{sample}.unfiltered.bam"
        filtered_bam = sample_align_dir / f"{sample}.bam"
        flagstat_path = sample_align_dir / f"{sample}.flagstat.txt"
        idxstats_path = sample_align_dir / f"{sample}.idxstats.tsv"
        depth_path = sample_align_dir / f"{sample}.depth.tsv"
        gap_path = sample_gap_dir / f"{sample}.gaps.tsv"

        align_cmd = (
            f"minimap2 -ax {args.aligner_preset} -t {args.threads} {reference} {fastq} "
            f"| samtools sort -@ {args.threads} -o {unfiltered_bam} -"
        )
        run_command(align_cmd, shell=True)

        filter_cmd = (
            f"samtools view -@ {args.threads} -b -q {args.min_mapq} -F 2308 {unfiltered_bam} "
            f"| samtools sort -@ {args.threads} -o {filtered_bam} -"
        )
        run_command(filter_cmd, shell=True)

        run_command(["samtools", "index", str(filtered_bam)])

        with flagstat_path.open("w") as handle:
            subprocess.run(["samtools", "flagstat", str(filtered_bam)], check=True, stdout=handle)
        with idxstats_path.open("w") as handle:
            subprocess.run(["samtools", "idxstats", str(filtered_bam)], check=True, stdout=handle)
        with depth_path.open("w") as handle:
            subprocess.run(["samtools", "depth", "-aa", str(filtered_bam)], check=True, stdout=handle)

        write_gap_summary(depth_path, reference_lengths, args.min_depth_ko, gap_path)

        if not args.keep_unfiltered_bam:
            unfiltered_bam.unlink(missing_ok=True)

        flagstat_total, mapped_reads = parse_flagstat_mapped(flagstat_path)
        manifest_rows.append(
            {
                "sample": sample,
                "fastq": str(fastq),
                "fastq_reads": str(count_fastq_reads(fastq)),
                "bam": f"{sample}/{sample}.bam",
                "bai": f"{sample}/{sample}.bam.bai",
                "flagstat": f"{sample}/{sample}.flagstat.txt",
                "idxstats": f"{sample}/{sample}.idxstats.tsv",
                "depth": f"{sample}/{sample}.depth.tsv",
                "gaps": f"../gaps/{sample}/{sample}.gaps.tsv",
                "flagstat_total": flagstat_total,
                "mapped_reads": mapped_reads,
            }
        )

    manifest_path = alignments_dir / "sample_bams.tsv"
    fieldnames = [
        "sample",
        "fastq",
        "fastq_reads",
        "bam",
        "bai",
        "flagstat",
        "idxstats",
        "depth",
        "gaps",
        "flagstat_total",
        "mapped_reads",
    ]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    with (reports_dir / "alignment_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    render_alignment_summary(reports_dir / "alignment_summary.html", manifest_rows)

    with (reports_dir / "alignment_manifest.json").open("w") as handle:
        json.dump(
            {
                "reference": str(reference),
                "aligner": "minimap2",
                "aligner_preset": args.aligner_preset,
                "min_mapq": args.min_mapq,
                "samples": manifest_rows,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
