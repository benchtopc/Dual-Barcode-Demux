@echo off
REM Runs a local smoke test for Dual-Barcode-Demux.
python --version
python bin\demux_scm.py ^
  --fastq tests\data\demo.fastq.gz ^
  --plate-map templates\plate_map.template.csv ^
  --row-barcodes data\scm_row_barcodes.csv ^
  --column-barcodes data\scm_column_barcodes.csv ^
  --outdir local_test_output
echo Local smoke test complete. Check local_test_output\tables\status_counts.tsv
pause
