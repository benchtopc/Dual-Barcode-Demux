#!/usr/bin/env python3
# bin/call_variants_bcftools_scm.py

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterator, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Call per-sample small variants with bcftools.")
    parser.add_argument("--alignments-dir", required=True, help="Directory containing sample_bams.tsv and BAM folders.")
    parser.add_argument("--reference", required=True, help="Reference FASTA.")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-baseq", type=int, default=7)
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


def run_command(command: List[str] | str, *, shell: bool = False) -> None:
    printable = command if isinstance(command, str) else " ".join(command)
    print(f"[scm-bcftools] running: {printable}", flush=True)
    if shell:
        subprocess.run(["bash", "-o", "pipefail", "-lc", command], check=True)
    else:
        subprocess.run(command, check=True)


def ensure_tool(name: str) -> None:
    if not shutil.which(name):
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def safe_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def first_float(value: str) -> float:
    try:
        return float(value.split(",")[0])
    except (TypeError, ValueError, IndexError):
        return 0.0


def summarize_vcf(sample: str, vcf_path: Path, output_path: Path) -> None:
    with open_vcf(vcf_path) as handle, output_path.open("w") as out:
        out.write("sample\tcontig\tpos\tref\talt\tqual\tfilter\tgenotype\tdepth\tad_ref\tad_alt\tallele_fraction\n")
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            contig, pos, _vid, ref, alts, qual, filt, info, fmt, sample_field = fields[:10]
            alt = alts.split(",")[0]
            fmt_map = dict(zip(fmt.split(":"), sample_field.split(":")))
            gt = fmt_map.get("GT", "./.")
            dp = safe_int(fmt_map.get("DP", "0"))
            ad_ref = 0
            ad_alt = 0
            if "AD" in fmt_map:
                ad_values = [safe_int(item) for item in fmt_map["AD"].split(",")]
                ad_ref = ad_values[0] if ad_values else 0
                ad_alt = ad_values[1] if len(ad_values) > 1 else 0
            af = first_float(fmt_map.get("AF", "0"))
            if not af and dp:
                af = ad_alt / dp if dp else 0.0
            out.write(
                f"{sample}\t{contig}\t{pos}\t{ref}\t{alt}\t{qual}\t{filt}\t{gt}\t{dp}\t{ad_ref}\t{ad_alt}\t{af:.4f}\n"
            )


def main() -> None:
    args = parse_args()
    for tool in ("bcftools", "samtools", "bash"):
        ensure_tool(tool)

    alignments_dir = Path(args.alignments_dir)
    reference = Path(args.reference)
    outdir = Path(args.outdir)
    variants_dir = outdir / "variants"
    variants_dir.mkdir(parents=True, exist_ok=True)

    run_command(["samtools", "faidx", str(reference)])

    rows = read_manifest(alignments_dir / "sample_bams.tsv")
    vcf_rows: List[Dict[str, str]] = []

    for row in rows:
        sample = row["sample"]
        bam = alignments_dir / row["bam"]
        sample_variant_dir = variants_dir / sample
        sample_variant_dir.mkdir(parents=True, exist_ok=True)

        vcf_gz = sample_variant_dir / f"{sample}.vcf.gz"
        call_cmd = (
            f"bcftools mpileup -Ou -f {reference} -Q {args.min_baseq} -q {args.min_mapq} {bam} "
            f"| bcftools call -mv -Oz -o {vcf_gz}"
        )
        run_command(call_cmd, shell=True)
        run_command(["bcftools", "index", "-t", str(vcf_gz)])

        summary_path = sample_variant_dir / f"{sample}.variants.tsv"
        summarize_vcf(sample, vcf_gz, summary_path)

        vcf_rows.append(
            {
                "sample": sample,
                "caller": "bcftools",
                "vcf": f"{sample}/{sample}.vcf.gz",
                "tbi": f"{sample}/{sample}.vcf.gz.tbi",
                "summary": f"{sample}/{sample}.variants.tsv",
            }
        )

    with (variants_dir / "sample_vcfs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["sample", "caller", "vcf", "tbi", "summary"])
        writer.writeheader()
        writer.writerows(vcf_rows)


if __name__ == "__main__":
    main()
