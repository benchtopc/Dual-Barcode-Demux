#!/usr/bin/env python3
# bin/analyze_ko_igv_scm.py

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class Variant:
    contig: str
    pos: int
    ref: str
    alt: str
    qual: str
    filt: str
    gt: str
    dp: int
    ad_ref: int
    ad_alt: int
    af: float
    info: str


@dataclass(frozen=True)
class KoTarget:
    target_id: str
    contig: str
    start: int
    end: int
    wt_allele: str
    ko_allele: str
    ko_type: str
    expected_effect: str
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify KO targets and generate IGV outputs from BAM/VCF outputs.")
    parser.add_argument("--alignments-dir", required=True)
    parser.add_argument("--variants-dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--ko-targets", required=True)
    parser.add_argument("--min-depth-ko", type=int, default=10)
    parser.add_argument("--het-min-af", type=float, default=0.25)
    parser.add_argument("--hom-alt-min-af", type=float, default=0.80)
    parser.add_argument("--igv-locus", default="all")
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


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


def parse_info(info: str) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            data[key] = value
    return data


def variant_end(variant: Variant) -> int:
    if variant.alt.startswith("<") and "END=" in variant.info:
        return safe_int(parse_info(variant.info).get("END", str(variant.pos)))
    return variant.pos + max(len(variant.ref), 1) - 1


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_ko_targets(path: Path) -> List[KoTarget]:
    targets: List[KoTarget] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_id", "contig", "start", "end", "wt_allele", "ko_allele"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            ko_type = (row.get("ko_type") or "allele").strip().lower()
            targets.append(
                KoTarget(
                    target_id=row["target_id"].strip(),
                    contig=row["contig"].strip(),
                    start=int(row["start"]),
                    end=int(row["end"]),
                    wt_allele=row["wt_allele"].strip().upper(),
                    ko_allele=row["ko_allele"].strip().upper(),
                    ko_type=ko_type,
                    expected_effect=(row.get("expected_effect") or "").strip(),
                    notes=(row.get("notes") or "").strip(),
                )
            )
    return targets


def parse_depth_file(path: Path) -> Dict[Tuple[str, int], int]:
    depth: Dict[Tuple[str, int], int] = {}
    with path.open("r") as handle:
        for line in handle:
            if not line.strip():
                continue
            contig, pos, value = line.rstrip("\n").split("\t")[:3]
            depth[(contig, int(pos))] = int(value)
    return depth


def parse_sample_vcf(vcf_path: Path) -> List[Variant]:
    variants: List[Variant] = []
    if not vcf_path.exists():
        return variants

    with open_vcf(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            contig, pos_s, _vid, ref, alt_s, qual, filt, info, fmt, sample = fields[:10]
            alt = alt_s.split(",")[0]
            fmt_keys = fmt.split(":")
            sample_values = sample.split(":")
            fmt_map = dict(zip(fmt_keys, sample_values))
            gt = fmt_map.get("GT", "./.")
            info_map = parse_info(info)

            dp = safe_int(fmt_map.get("DP", info_map.get("DP", "0")))
            ad_ref = 0
            ad_alt = 0
            if "AD" in fmt_map:
                ad_values = [safe_int(item) for item in fmt_map["AD"].split(",")]
                ad_ref = ad_values[0] if ad_values else 0
                ad_alt = ad_values[1] if len(ad_values) > 1 else 0

            af = first_float(fmt_map.get("AF", "0"))
            if af == 0.0:
                af = first_float(info_map.get("AF", "0"))
            if af == 0.0 and dp:
                af = ad_alt / dp if dp else 0.0
            if af == 0.0:
                normalized_gt = gt.replace("|", "/")
                if normalized_gt == "1/1":
                    af = 1.0
                elif normalized_gt in {"0/1", "1/0"}:
                    af = 0.5

            variants.append(
                Variant(
                    contig=contig,
                    pos=int(pos_s),
                    ref=ref.upper(),
                    alt=alt.upper(),
                    qual=qual,
                    filt=filt,
                    gt=gt,
                    dp=dp,
                    ad_ref=ad_ref,
                    ad_alt=ad_alt,
                    af=af,
                    info=info,
                )
            )
    return variants


def select_variants_for_target(target: KoTarget, variants: List[Variant]) -> List[Variant]:
    selected: List[Variant] = []
    for variant in variants:
        if variant.contig != target.contig:
            continue
        if intervals_overlap(target.start, target.end, variant.pos, variant_end(variant)):
            selected.append(variant)
    return selected


def variant_matches_ko(target: KoTarget, variant: Variant) -> bool:
    ko_allele = target.ko_allele.upper()
    ko_type = target.ko_type.lower()

    if ko_type == "deletion" or ko_allele in {"<DEL>", "DEL", "-"}:
        if variant.alt == "<DEL>":
            return intervals_overlap(target.start, target.end, variant.pos, variant_end(variant))
        return len(variant.ref) > len(variant.alt) and intervals_overlap(target.start, target.end, variant.pos, variant_end(variant))

    if ko_type == "insertion" or ko_allele in {"<INS>", "INS"}:
        if variant.alt == "<INS>":
            return intervals_overlap(target.start, target.end, variant.pos, variant_end(variant))
        return len(variant.alt) > len(variant.ref) and (ko_allele in {"<INS>", "INS"} or ko_allele in variant.alt)

    if ko_type == "any_variant":
        return intervals_overlap(target.start, target.end, variant.pos, variant_end(variant))

    return variant.alt == ko_allele or variant.alt.endswith(ko_allele)


def classify_ko(
    target: KoTarget,
    candidate_variants: List[Variant],
    depth: int,
    min_depth: int,
    het_min_af: float,
    hom_alt_min_af: float,
) -> Tuple[str, str, str, float, str, Optional[Variant]]:
    if depth < min_depth:
        return "no_call_low_depth", "low", str(depth), 0.0, "Depth below threshold", None

    if not candidate_variants:
        return "WT", "high", str(depth), 0.0, "No variant overlaps target interval", None

    ko_variants = [variant for variant in candidate_variants if variant_matches_ko(target, variant)]
    if not ko_variants:
        observed = ";".join(f"{v.pos}:{v.ref}>{v.alt}" for v in candidate_variants[:5])
        return "variant_non_ko", "medium", str(depth), 0.0, f"Variant(s) present but not expected KO allele: {observed}", candidate_variants[0]

    best = max(ko_variants, key=lambda item: item.af)
    genotype = best.gt.replace("|", "/")
    call_depth = str(best.dp or depth)

    if genotype == "1/1" or best.af >= hom_alt_min_af:
        return "homozygous_KO", "high", call_depth, best.af, "KO allele dominant", best
    if genotype in {"0/1", "1/0"} or best.af >= het_min_af:
        return "heterozygous_KO", "medium", call_depth, best.af, "KO allele present below homozygous threshold", best

    return "WT", "medium", call_depth, best.af, "KO allele below heterozygous threshold", best


def summarize_vcf(sample: str, vcf_path: Path, output_path: Path) -> int:
    variants = parse_sample_vcf(vcf_path)
    with output_path.open("w") as out:
        out.write("sample\tcontig\tpos\tend\tref\talt\tqual\tfilter\tgenotype\tdepth\tad_ref\tad_alt\tallele_fraction\n")
        for variant in variants:
            out.write(
                f"{sample}\t{variant.contig}\t{variant.pos}\t{variant_end(variant)}\t{variant.ref}\t{variant.alt}\t"
                f"{variant.qual}\t{variant.filt}\t{variant.gt}\t{variant.dp}\t{variant.ad_ref}\t{variant.ad_alt}\t{variant.af:.4f}\n"
            )
    return len(variants)


def write_fai(reference: Path, fai_path: Path) -> None:
    rows: List[Tuple[str, int, int, int, int]] = []
    with reference.open("rb") as handle:
        name: Optional[str] = None
        seq_len = 0
        seq_offset = 0
        line_bases = 0
        line_width = 0
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                if name is not None:
                    rows.append((name, seq_len, seq_offset, line_bases, line_width))
                name = line[1:].decode().strip().split()[0]
                seq_len = 0
                seq_offset = handle.tell()
                line_bases = 0
                line_width = 0
            else:
                stripped = line.rstrip(b"\r\n")
                if line_bases == 0 and stripped:
                    line_bases = len(stripped)
                    line_width = len(line)
                seq_len += len(stripped)
        if name is not None:
            rows.append((name, seq_len, seq_offset, line_bases, line_width))

    with fai_path.open("w") as out:
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")


def write_ko_calls(
    sample: str,
    targets: List[KoTarget],
    depth_map: Dict[Tuple[str, int], int],
    variants: List[Variant],
    min_depth: int,
    het_min_af: float,
    hom_alt_min_af: float,
    output_path: Path,
) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "target_id",
                "contig",
                "start",
                "end",
                "wt_allele",
                "ko_allele",
                "ko_type",
                "classification",
                "confidence",
                "depth",
                "allele_fraction",
                "genotype",
                "variant_pos",
                "ref",
                "alt",
                "qual",
                "expected_effect",
                "notes",
                "interpretation",
            ]
        )
        for target in targets:
            interval_depths = [depth_map.get((target.contig, pos), 0) for pos in range(target.start, target.end + 1)]
            depth = min(interval_depths) if interval_depths else 0
            candidates = select_variants_for_target(target, variants)
            classification, confidence, call_depth, af, interpretation, variant = classify_ko(
                target=target,
                candidate_variants=candidates,
                depth=depth,
                min_depth=min_depth,
                het_min_af=het_min_af,
                hom_alt_min_af=hom_alt_min_af,
            )
            counts[classification] += 1
            writer.writerow(
                [
                    sample,
                    target.target_id,
                    target.contig,
                    target.start,
                    target.end,
                    target.wt_allele,
                    target.ko_allele,
                    target.ko_type,
                    classification,
                    confidence,
                    call_depth,
                    f"{af:.4f}",
                    variant.gt if variant else "0/0",
                    variant.pos if variant else "",
                    variant.ref if variant else target.wt_allele,
                    variant.alt if variant else ".",
                    variant.qual if variant else ".",
                    target.expected_effect,
                    target.notes,
                    interpretation,
                ]
            )
    return counts


def write_igv_session(sample: str, reference: Path, bam: str, vcf: str, output_path: Path) -> None:
    output_path.write_text(
        f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{html.escape(str(reference))}" hasGeneTrack="false" hasSequenceTrack="true" version="8">
  <Resources>
    <Resource path="{html.escape(str(reference))}"/>
    <Resource path="{html.escape(bam)}"/>
    <Resource path="{html.escape(vcf)}"/>
  </Resources>
  <Panel name="DataPanel">
    <Track attributeKey="{html.escape(sample)} BAM" id="{html.escape(bam)}" name="{html.escape(sample)} alignment"/>
    <Track attributeKey="{html.escape(sample)} VCF" id="{html.escape(vcf)}" name="{html.escape(sample)} variants"/>
  </Panel>
</Session>
""",
        encoding="utf-8",
    )


def write_igv_html(samples: List[Dict[str, str]], reference_fasta: str, reference_fai: str, locus: str, output_path: Path) -> None:
    tracks = []
    for item in samples:
        tracks.append(
            {
                "name": f"{item['sample']} alignment",
                "type": "alignment",
                "format": "bam",
                "url": item["bam"],
                "indexURL": item["bai"],
            }
        )
        if item.get("vcf"):
            track = {
                "name": f"{item['sample']} variants",
                "type": "variant",
                "format": "vcf",
                "url": item["vcf"],
            }
            if item.get("tbi"):
                track["indexURL"] = item["tbi"]
            tracks.append(track)

    config = {
        "reference": {
            "id": "custom_reference",
            "name": "Custom reference",
            "fastaURL": reference_fasta,
            "indexURL": reference_fai,
        },
        "locus": locus,
        "tracks": tracks,
    }

    output_path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>SCM IGV alignment viewer</title>
  <script src="https://cdn.jsdelivr.net/npm/igv@2.15.12/dist/igv.min.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 1.5rem; color: #222; }}
    .note {{ border: 1px solid #ddd; border-radius: 8px; padding: 1rem; background: #fafafa; }}
    #igv-div {{ height: 760px; border: 1px solid #ddd; margin-top: 1rem; }}
    code {{ background: #f1f1f1; padding: 0.1rem 0.3rem; border-radius: 4px; }}
  </style>
</head>
<body>
  <h1>SCM per-sample alignment viewer</h1>
  <div class="note">
    <p>This IGV.js report loads the copied reference FASTA plus workflow BAM/BAI and VCF/TBI files.</p>
    <p>If your browser blocks local loading, run <code>python3 -m http.server 8000</code> inside the output directory and open <code>http://localhost:8000/igv/igv_viewer.html</code>.</p>
  </div>
  <div id="igv-div"></div>
  <script>
    const config = {json.dumps(config, indent=2)};
    igv.createBrowser(document.getElementById("igv-div"), config);
  </script>
</body>
</html>
""",
        encoding="utf-8",
    )


def render_summary(report_path: Path, rows: List[Dict[str, str]], variant_caller: str) -> None:
    body_rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(row['sample'])}</td>"
        f"<td>{html.escape(row['fastq_reads'])}</td>"
        f"<td>{html.escape(row['mapped_reads'])}</td>"
        f"<td>{html.escape(row['variant_count'])}</td>"
        f"<td>{html.escape(row['ko_status_summary'])}</td>"
        "</tr>"
        for row in rows
    )
    report_path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SCM KO and IGV summary</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; }}
th {{ background: #f5f5f5; }}
.card {{ border: 1px solid #ddd; border-radius: 10px; padding: 1rem; margin-bottom: 1rem; }}
</style>
</head>
<body>
<h1>SCM variant, gap, and KO summary</h1>
<div class="card">
<p><strong>Variant caller:</strong> {html.escape(variant_caller)}</p>
<p>Open <code>igv/igv_viewer.html</code> for alignment and variant review.</p>
</div>
<table>
<thead><tr><th>Sample</th><th>FASTQ reads</th><th>Mapped reads</th><th>Variants</th><th>KO calls</th></tr></thead>
<tbody>{body_rows or "<tr><td colspan='5'>No samples found</td></tr>"}</tbody>
</table>
</body>
</html>
""",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()

    alignments_dir = Path(args.alignments_dir)
    variants_dir = Path(args.variants_dir)
    reference = Path(args.reference)
    ko_targets_path = Path(args.ko_targets)
    outdir = Path(args.outdir)

    ko_dir = outdir / "ko_calls"
    igv_dir = outdir / "igv"
    igv_sessions_dir = igv_dir / "sessions"
    igv_ref_dir = igv_dir / "reference"
    reports_dir = outdir / "analysis_reports"
    for directory in (ko_dir, igv_dir, igv_sessions_dir, igv_ref_dir, reports_dir):
        directory.mkdir(parents=True, exist_ok=True)

    alignment_rows = read_manifest(alignments_dir / "sample_bams.tsv")
    vcf_rows = {row["sample"]: row for row in read_manifest(variants_dir / "sample_vcfs.tsv")}
    targets = read_ko_targets(ko_targets_path)

    reference_copy = igv_ref_dir / reference.name
    shutil.copy2(reference, reference_copy)
    fai_copy = Path(f"{reference_copy}.fai")
    source_fai = Path(f"{reference}.fai")
    if source_fai.exists():
        shutil.copy2(source_fai, fai_copy)
    else:
        write_fai(reference_copy, fai_copy)

    summary_rows: List[Dict[str, str]] = []
    igv_samples: List[Dict[str, str]] = []
    variant_caller = "unknown"

    merged_ko_path = ko_dir / "ko_calls.summary.tsv"
    with merged_ko_path.open("w", newline="") as merged_handle:
        merged_writer: Optional[csv.writer] = None

        for row in alignment_rows:
            sample = row["sample"]
            sample_ko_dir = ko_dir / sample
            sample_ko_dir.mkdir(parents=True, exist_ok=True)

            vcf_row = vcf_rows.get(sample)
            if not vcf_row:
                raise RuntimeError(f"No VCF listed for sample {sample}")
            variant_caller = vcf_row.get("caller", variant_caller)
            vcf_path = variants_dir / vcf_row["vcf"]
            variant_summary_path = variants_dir / sample / f"{sample}.variants.tsv"
            variant_count = summarize_vcf(sample, vcf_path, variant_summary_path)

            depth_path = alignments_dir / row["depth"]
            depth_map = parse_depth_file(depth_path)
            variants = parse_sample_vcf(vcf_path)

            ko_path = sample_ko_dir / f"{sample}.ko_calls.tsv"
            counts = write_ko_calls(
                sample=sample,
                targets=targets,
                depth_map=depth_map,
                variants=variants,
                min_depth=args.min_depth_ko,
                het_min_af=args.het_min_af,
                hom_alt_min_af=args.hom_alt_min_af,
                output_path=ko_path,
            )

            with ko_path.open("r") as sample_ko_handle:
                reader = csv.reader(sample_ko_handle, delimiter="\t")
                header = next(reader)
                if merged_writer is None:
                    merged_writer = csv.writer(merged_handle, delimiter="\t")
                    merged_writer.writerow(header)
                for ko_row in reader:
                    merged_writer.writerow(ko_row)

            ko_summary = ", ".join(f"{key}:{value}" for key, value in sorted(counts.items())) or "none"
            summary_rows.append(
                {
                    "sample": sample,
                    "fastq_reads": row.get("fastq_reads", "0"),
                    "mapped_reads": row.get("mapped_reads", "0"),
                    "variant_count": str(variant_count),
                    "ko_status_summary": ko_summary,
                }
            )

            bam_rel = f"../alignments/{row['bam']}"
            bai_rel = f"../alignments/{row['bai']}"
            vcf_rel = f"../variants/{vcf_row['vcf']}"
            tbi_rel = f"../variants/{vcf_row['tbi']}" if vcf_row.get("tbi") else ""
            write_igv_session(sample, reference_copy, bam_rel, vcf_rel, igv_sessions_dir / f"{sample}.igv_session.xml")
            igv_samples.append({"sample": sample, "bam": bam_rel, "bai": bai_rel, "vcf": vcf_rel, "tbi": tbi_rel})

    write_igv_html(
        samples=igv_samples,
        reference_fasta=f"reference/{reference_copy.name}",
        reference_fai=f"reference/{reference_copy.name}.fai",
        locus=args.igv_locus,
        output_path=igv_dir / "igv_viewer.html",
    )
    render_summary(reports_dir / "ko_igv_summary.html", summary_rows, variant_caller)

    with (reports_dir / "analysis_manifest.json").open("w") as handle:
        json.dump(
            {
                "reference": str(reference),
                "ko_targets": str(ko_targets_path),
                "variant_caller": variant_caller,
                "min_depth_ko": args.min_depth_ko,
                "het_min_af": args.het_min_af,
                "hom_alt_min_af": args.hom_alt_min_af,
                "samples": summary_rows,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
