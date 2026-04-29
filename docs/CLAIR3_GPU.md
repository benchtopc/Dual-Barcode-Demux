# Clair3 GPU configuration

This workflow supports three variant-calling modes:

```text
bcftools
clair3
clair3_gpu
```

Use `bcftools` for the portable CPU fallback and CI/demo runs. Use `clair3` or `clair3_gpu` for higher-accuracy long-read small-variant calling after selecting the right model for your chemistry/basecaller.

## Bundled container model path

The default Clair3 v2 GPU/CPU containers include models under:

```text
/opt/models/
```

Default:

```text
--clair3_model_path /opt/models/r1041_e82_400bps_sup_v500
```

## Local model directory

Use this when you downloaded or trained a model:

```bash
--clair3_model_dir /absolute/path/to/model_directory
```

The directory is staged by Nextflow and passed to Clair3 as `--model_path`.

## GPU run

```bash
nextflow run . \
  -profile docker,gpu \
  --variant_caller clair3_gpu \
  --enable_gpu true \
  --gpu_devices all \
  --clair3_device cuda:0 \
  --fastq reads.fastq.gz \
  --plate_map plate_map.csv \
  --run_alignment true \
  --reference reference.fa \
  --ko_targets ko_targets.tsv
```

For multiple GPUs:

```bash
--gpu_devices '"device=0,1"' --clair3_device cuda:0,1
```

## Recommended custom-reference settings

For bacterial, plasmid, amplicon, custom target, or other non-human references:

```bash
--clair3_include_all_ctgs true
```

For a specific region:

```bash
--clair3_bed regions.bed
```

For a specific contig list:

```bash
--clair3_ctg_name "contig1,contig2"
```

## Signal-aware models

Use only when the BAM contains Dorado `mv` tags and the model is compatible:

```bash
--clair3_enable_dwell_time true
```

FASTQ-only workflows normally cannot retain move-table signal tags.
