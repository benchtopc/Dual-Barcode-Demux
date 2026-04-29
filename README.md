# Dual-Barcode-Demux

**Maintainer:** SCM  
**License:** Apache-2.0  
**Workflow type:** EPI2ME-compatible Nextflow workflow  
**Purpose:** Demultiplex 96-well combinatorial row/column nanopore barcodes, align per-sample reads to custom references, call variants with either `bcftools` or Clair3/Clair3 GPU, classify KO status, and generate IGV viewing reports.

> Research and educational use only. Not for clinical diagnosis, patient care, or regulated use without independent validation.

## Repository layout required by EPI2ME

EPI2ME GitHub import expects these files directly at the repository root:

```text
README.md
main.nf
nextflow.config
nextflow_schema.json
```

Do not upload this project inside a nested folder. If EPI2ME says `README.md not found`, the workflow files are probably one directory too deep or the repository root is incomplete.

## Features

- 96-well combinatorial dual-barcode demultiplexing from one untrimmed FASTQ/FASTQ.GZ
- Per-sample FASTQ output
- `minimap2` per-sample long-read alignment with default ONT preset `map-ont`
- Sorted/indexed BAM output with `samtools` flagstat, idxstats, depth, and gap tables
- Variant calling options:
  - `bcftools` CPU fallback
  - `clair3` CPU deep-learning caller
  - `clair3_gpu` NVIDIA GPU deep-learning caller
- Full Clair3 model-path support:
  - bundled container model path, for example `/opt/models/r1041_e82_400bps_sup_v500`
  - optional staged local model directory with `--clair3_model_dir`
  - platform selection: `ont`, `hifi`, `ilmn`
  - `--include_all_ctgs`, `--bed_fn`, `--ctg_name`, haploid modes, dwell-time model support, gVCF/ref-call options, and advanced extra args
- KO target classification:
  - `homozygous_KO`
  - `heterozygous_KO`
  - `WT`
  - `variant_non_ko`
  - `no_call_low_depth`
- KO target modes:
  - `allele`
  - `deletion`
  - `insertion`
  - `any_variant`
- IGV outputs:
  - browser-based `igv/igv_viewer.html`
  - per-sample IGV Desktop XML sessions

## Required demultiplexing inputs

| Parameter | Description |
|---|---|
| `--fastq` | One untrimmed FASTQ or FASTQ.GZ file |
| `--plate_map` | CSV mapping wells to row barcode IDs, column barcode IDs, and sample IDs |

Default barcode tables:

```text
data/scm_row_barcodes.csv
data/scm_column_barcodes.csv
```

Required plate-map columns:

```text
well_id,row_id,column_id,row_barcode_id,column_barcode_id,sample_id,alias
```

Editable template:

```text
templates/plate_map.template.csv
```

## Alignment and KO inputs

Enable the analysis stage:

```bash
--run_alignment true
```

Then provide:

| Parameter | Description |
|---|---|
| `--reference` | Custom reference FASTA |
| `--ko_targets` | TSV of target sites and expected KO alleles |
| `--variant_caller` | `bcftools`, `clair3`, or `clair3_gpu` |

KO target template:

```text
templates/ko_targets.template.tsv
```

Required KO columns:

```text
target_id,contig,start,end,wt_allele,ko_allele
```

Optional KO columns:

```text
ko_type,expected_effect,notes
```

Coordinates are 1-based reference coordinates.

## EPI2ME import

Use GitHub import, not a locally built `.2me` file:

```text
Launch / Workflows -> Import workflow -> GitHub URL
```

Paste your repository URL without `.git`:

```text
https://github.com/benchtopc/Dual-Barcode-Demux
```

## Demux-only run

```bash
nextflow run . \
  -profile docker \
  --fastq /path/to/reads.fastq.gz \
  --plate_map /path/to/plate_map.csv \
  --out_dir output
```

## Demux + alignment + bcftools

```bash
nextflow run . \
  -profile docker \
  --fastq /path/to/reads.fastq.gz \
  --plate_map /path/to/plate_map.csv \
  --run_alignment true \
  --reference /path/to/reference.fa \
  --ko_targets /path/to/ko_targets.tsv \
  --variant_caller bcftools \
  --out_dir output
```

## Demux + alignment + Clair3 CPU

Use a bundled model path from the Clair3 container:

```bash
nextflow run . \
  -profile docker \
  --fastq /path/to/reads.fastq.gz \
  --plate_map /path/to/plate_map.csv \
  --run_alignment true \
  --reference /path/to/reference.fa \
  --ko_targets /path/to/ko_targets.tsv \
  --variant_caller clair3 \
  --clair3_platform ont \
  --clair3_model_path /opt/models/r1041_e82_400bps_sup_v500 \
  --out_dir output
```

Use a local model directory instead:

```bash
nextflow run . \
  -profile docker \
  --fastq /path/to/reads.fastq.gz \
  --plate_map /path/to/plate_map.csv \
  --run_alignment true \
  --reference /path/to/reference.fa \
  --ko_targets /path/to/ko_targets.tsv \
  --variant_caller clair3 \
  --clair3_model_dir /path/to/clair3_model_directory \
  --out_dir output
```

The local model directory should contain the Clair3 v2 PyTorch model files expected by `run_clair3.sh`, including model files such as `pileup.pt` and `full_alignment.pt`.

## Demux + alignment + Clair3 GPU

Requires Docker, NVIDIA driver, NVIDIA Container Toolkit, and GPU passthrough into Docker/WSL.

```bash
nextflow run . \
  -profile docker,gpu \
  --fastq /path/to/reads.fastq.gz \
  --plate_map /path/to/plate_map.csv \
  --run_alignment true \
  --reference /path/to/reference.fa \
  --ko_targets /path/to/ko_targets.tsv \
  --variant_caller clair3_gpu \
  --clair3_platform ont \
  --clair3_model_path /opt/models/r1041_e82_400bps_sup_v500 \
  --clair3_device cuda:0 \
  --out_dir output
```

Multiple GPUs:

```bash
nextflow run . \
  -profile docker,gpu \
  --gpu_devices '"device=0,1"' \
  --clair3_device cuda:0,1 \
  --variant_caller clair3_gpu \
  --fastq reads.fastq.gz \
  --plate_map plate_map.csv \
  --run_alignment true \
  --reference reference.fa \
  --ko_targets ko_targets.tsv
```

## Clair3 accuracy controls

Recommended for custom/non-human references:

```bash
--clair3_include_all_ctgs true
```

Call selected contigs:

```bash
--clair3_ctg_name "chr1,chr2"
```

Call selected BED regions:

```bash
--clair3_bed /path/to/regions.bed
```

Amplicon/high-depth custom target options to consider after validation:

```bash
--clair3_no_phasing_for_fa true
--clair3_extra_args "--var_pct_full=1 --ref_pct_full=1"
```

Signal-aware ONT model support:

```bash
--clair3_enable_dwell_time true
```

This requires BAMs containing Dorado `mv` tags and a compatible signal-aware Clair3 model. FASTQ-only workflows normally do not preserve `mv` tags.

Haploid organisms:

```bash
--clair3_haploid_mode haploid_precise
```

or:

```bash
--clair3_haploid_mode haploid_sensitive
```

## Demo run

```bash
nextflow run . -profile docker,demo
```

The demo uses the portable `bcftools` path so CI can run without a GPU or Clair3 model download.

## Main outputs

```text
demultiplexed_fastq/
tables/
reports/
metadata/
alignments/
gaps/
variants/
clair3_runs/
ko_calls/
igv/
analysis_reports/
```

Key files:

```text
tables/read_assignments.tsv.gz
reports/demux_summary.html
alignments/sample_bams.tsv
alignments/<sample>/<sample>.bam
alignments/<sample>/<sample>.bam.bai
gaps/<sample>/<sample>.gaps.tsv
variants/sample_vcfs.tsv
variants/<sample>/<sample>.vcf.gz
variants/<sample>/<sample>.variants.tsv
ko_calls/ko_calls.summary.tsv
igv/igv_viewer.html
igv/sessions/<sample>.igv_session.xml
analysis_reports/alignment_summary.html
analysis_reports/ko_igv_summary.html
```

## IGV viewing

Open:

```text
igv/igv_viewer.html
```

If your browser blocks local BAM/VCF loading, serve the output directory:

```bash
cd output
python3 -m http.server 8000
```

Then open:

```text
http://localhost:8000/igv/igv_viewer.html
```

You can also open the XML sessions in IGV Desktop:

```text
igv/sessions/<sample>.igv_session.xml
```

## License and disclaimer

This workflow is licensed under Apache-2.0 and includes a separate `DISCLAIMER.md`.

This workflow is provided for research and educational use only. It is not intended for clinical diagnosis, medical decision-making, patient care, or regulated laboratory use without independent validation.

The software is provided "AS IS", without warranties or conditions of any kind. To the maximum extent permitted by applicable law, SCM and contributors are not liable for claims, damages, data loss, analysis errors, barcode misclassification, variant-calling errors, KO misclassification, regulatory issues, or other liability arising from use, modification, distribution, or inability to use this workflow.


## Replace broken GitHub upload

If EPI2ME reports `README.md not found` or `Error loading workflow from file`, your repo was likely uploaded with broken/nested files.

From the unzipped clean workflow folder, run:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\replace_github_repo.ps1
```

Then import this exact URL in EPI2ME:

```text
https://github.com/benchtopc/Dual-Barcode-Demux
```

