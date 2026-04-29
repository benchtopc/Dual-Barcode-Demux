# install_or_test_local.ps1
# Runs a local smoke test for Dual-Barcode-Demux.
# This does not replace EPI2ME Desktop import; it confirms the workflow files are usable.

$ErrorActionPreference = "Stop"

Write-Host "Testing Dual-Barcode-Demux local workflow..."

if (-not (Test-Path "main.nf")) {
    throw "Run this script from the workflow root folder where main.nf exists."
}

python --version
python .\bin\demux_scm.py `
  --fastq .\tests\data\demo.fastq.gz `
  --plate-map .\templates\plate_map.template.csv `
  --row-barcodes .\data\scm_row_barcodes.csv `
  --column-barcodes .\data\scm_column_barcodes.csv `
  --outdir .\local_test_output

Write-Host "Local smoke test complete. Check local_test_output\tables\status_counts.tsv"
