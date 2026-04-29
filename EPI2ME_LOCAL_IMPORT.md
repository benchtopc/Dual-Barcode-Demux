# Local EPI2ME Import Notes

This folder is a complete Nextflow workflow application for **Dual-Barcode-Demux**.

## Preferred import path

EPI2ME Desktop's documented third-party workflow path is GitHub import:

```text
Launch -> Import workflow -> GitHub URL
```

Use:

```text
https://github.com/benchtopc/Dual-Barcode-Demux
```

## Local file import path

If your EPI2ME Desktop version shows:

```text
Launch -> Import workflow -> Import a 2ME file
```

try the file:

```text
Dual-Barcode-Demux-local-import-v3.2me
```

This package is a compressed TAR-format source package with `main.nf`,
`nextflow.config`, and `nextflow_schema.json` at the archive root.

## If local import says "Failed to unpack"

That means the installed EPI2ME version is enforcing ONT's internal 2ME
package structure. In that case, use GitHub import. The public documentation
does not currently provide a community 2ME package-builder specification.

## Command-line test

From this folder:

```bash
nextflow run . -profile standard,demo --out_dir test_output
```

Or run the Python demultiplexer directly:

```bash
python3 bin/demux_scm.py   --fastq tests/data/demo.fastq.gz   --plate-map templates/plate_map.template.csv   --row-barcodes data/scm_row_barcodes.csv   --column-barcodes data/scm_column_barcodes.csv   --outdir test_output
```
