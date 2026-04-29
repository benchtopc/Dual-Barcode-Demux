#!/usr/bin/env python3
# bin/run_clair3_batch_scm.py

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Clair3 per sample from an alignment manifest.")
    parser.add_argument("--alignments-dir", required=True, help="Directory containing sample_bams.tsv and BAM folders.")
    parser.add_argument("--reference", required=True, help="Reference FASTA.")
    parser.add_argument("--model-path", required=True, help="Clair3 model directory visible inside the container.")
    parser.add_argument("--platform", default="ont", choices=["ont", "hifi", "ilmn"], help="Clair3 platform.")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--use-gpu", action="store_true", help="Pass --use_gpu to Clair3.")
    parser.add_argument("--device", default="", help="Clair3 GPU device string, e.g. cuda:0 or cuda:0,1.")
    parser.add_argument("--bed-fn", default="", help="Optional BED file for region-limited calling.")
    parser.add_argument("--ctg-name", default="", help="Optional Clair3 ctg_name value.")
    parser.add_argument("--enable-dwell-time", action="store_true", help="Enable mv-tag/dwell-time model support.")
    parser.add_argument("--include-all-ctgs", action="store_true", help="Pass --include_all_ctgs.")
    parser.add_argument("--haploid-mode", choices=["none", "haploid_precise", "haploid_sensitive"], default="none")
    parser.add_argument("--no-phasing-for-fa", action="store_true")
    parser.add_argument("--pileup-only", action="store_true")
    parser.add_argument("--print-ref-calls", action="store_true")
    parser.add_argument("--gvcf", action="store_true")
    parser.add_argument("--qual", type=int, default=2)
    parser.add_argument("--chunk-size", type=int, default=0)
    parser.add_argument("--remove-intermediate-dir", action="store_true")
    parser.add_argument("--extra-args", default="", help="Advanced Clair3 args, quoted as one string.")
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


def run_command(command: List[str]) -> None:
    print(f"[scm-clair3] running: {' '.join(command)}", flush=True)
    subprocess.run(command, check=True)


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_clair3_runner() -> str:
    for candidate in ("/opt/bin/run_clair3.sh", "/opt/bin/run_clair3.py", "run_clair3.sh", "run_clair3.py"):
        resolved = shutil.which(candidate) if not candidate.startswith("/") else candidate
        if resolved and Path(resolved).exists():
            return resolved
    raise RuntimeError("Could not find run_clair3.sh or run_clair3.py in the Clair3 container.")


def ensure_reference_index(reference: Path) -> None:
    if Path(f"{reference}.fai").exists():
        return
    samtools = shutil.which("samtools")
    if not samtools:
        raise RuntimeError("Reference FASTA is not indexed and samtools is not available in the Clair3 container.")
    run_command([samtools, "faidx", str(reference)])


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()

    alignments_dir = Path(args.alignments_dir)
    reference = Path(args.reference)
    model_path = Path(args.model_path)
    outdir = Path(args.outdir)
    variants_dir = outdir / "variants"
    clair3_runs_dir = outdir / "clair3_runs"
    variants_dir.mkdir(parents=True, exist_ok=True)
    clair3_runs_dir.mkdir(parents=True, exist_ok=True)

    if not model_path.exists():
        raise RuntimeError(
            f"Clair3 model path does not exist inside the container: {model_path}. "
            "Use a bundled model path such as /opt/models/r1041_e82_400bps_sup_v500 or stage a local model directory."
        )

    bed_fn = Path(args.bed_fn) if args.bed_fn else None
    if bed_fn and (not bed_fn.exists() or bed_fn.stat().st_size == 0):
        bed_fn = None

    runner = find_clair3_runner()
    ensure_reference_index(reference)

    manifest_rows = read_manifest(alignments_dir / "sample_bams.tsv")
    vcf_rows: List[Dict[str, str]] = []
    run_records: List[Dict[str, str]] = []

    for row in manifest_rows:
        sample = row["sample"]
        bam = alignments_dir / row["bam"]
        if not bam.exists():
            raise RuntimeError(f"BAM listed in manifest does not exist: {bam}")

        sample_run_dir = clair3_runs_dir / sample
        sample_variant_dir = variants_dir / sample
        sample_run_dir.mkdir(parents=True, exist_ok=True)
        sample_variant_dir.mkdir(parents=True, exist_ok=True)

        command = [
            runner,
            f"--bam_fn={bam}",
            f"--ref_fn={reference}",
            f"--threads={args.threads}",
            f"--platform={args.platform}",
            f"--model_path={model_path}",
            f"--output={sample_run_dir}",
            f"--sample_name={sample}",
            f"--qual={args.qual}",
        ]

        if args.use_gpu:
            command.append("--use_gpu")
            if args.device:
                command.append(f"--device={args.device}")
        if bed_fn:
            command.append(f"--bed_fn={bed_fn}")
        if args.ctg_name:
            command.append(f"--ctg_name={args.ctg_name}")
        if args.enable_dwell_time:
            command.append("--enable_dwell_time")
        if args.include_all_ctgs:
            command.append("--include_all_ctgs")
        if args.haploid_mode != "none":
            command.append(f"--{args.haploid_mode}")
        if args.no_phasing_for_fa:
            command.append("--no_phasing_for_fa")
        if args.pileup_only:
            command.append("--pileup_only")
        if args.print_ref_calls:
            command.append("--print_ref_calls")
        if args.gvcf:
            command.append("--gvcf")
        if args.chunk_size > 0:
            command.append(f"--chunk_size={args.chunk_size}")
        if args.remove_intermediate_dir:
            command.append("--remove_intermediate_dir")
        if args.extra_args.strip():
            command.extend(shlex.split(args.extra_args))

        run_command(command)

        final_vcf = sample_run_dir / "merge_output.vcf.gz"
        if not final_vcf.exists():
            raise RuntimeError(f"Clair3 completed but final VCF was not found: {final_vcf}")

        sample_vcf = sample_variant_dir / f"{sample}.vcf.gz"
        shutil.copy2(final_vcf, sample_vcf)

        source_tbi = sample_run_dir / "merge_output.vcf.gz.tbi"
        sample_tbi = sample_variant_dir / f"{sample}.vcf.gz.tbi"
        if source_tbi.exists():
            shutil.copy2(source_tbi, sample_tbi)
        elif shutil.which("tabix"):
            run_command(["tabix", "-p", "vcf", str(sample_vcf)])

        copy_if_exists(sample_run_dir / "pileup.vcf.gz", sample_variant_dir / f"{sample}.clair3.pileup.vcf.gz")
        copy_if_exists(sample_run_dir / "full_alignment.vcf.gz", sample_variant_dir / f"{sample}.clair3.full_alignment.vcf.gz")
        copy_if_exists(sample_run_dir / "run_clair3.log", sample_variant_dir / f"{sample}.clair3.log")

        vcf_rows.append(
            {
                "sample": sample,
                "caller": "clair3_gpu" if args.use_gpu else "clair3",
                "vcf": f"{sample}/{sample}.vcf.gz",
                "tbi": f"{sample}/{sample}.vcf.gz.tbi" if sample_tbi.exists() else "",
                "summary": "",
            }
        )
        run_records.append(
            {
                "sample": sample,
                "bam": str(bam),
                "vcf": str(sample_vcf),
                "clair3_output": str(sample_run_dir),
                "model_path": str(model_path),
                "platform": args.platform,
                "use_gpu": str(args.use_gpu).lower(),
                "device": args.device,
            }
        )

    with (variants_dir / "sample_vcfs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["sample", "caller", "vcf", "tbi", "summary"])
        writer.writeheader()
        writer.writerows(vcf_rows)

    with (clair3_runs_dir / "clair3_manifest.json").open("w") as handle:
        json.dump(
            {
                "model_path": str(model_path),
                "platform": args.platform,
                "use_gpu": args.use_gpu,
                "device": args.device,
                "include_all_ctgs": args.include_all_ctgs,
                "enable_dwell_time": args.enable_dwell_time,
                "samples": run_records,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
