# Test data

The demo dataset contains two synthetic barcoded reads:

- `demo_forward_A01`
- `demo_reverse_H12`

The demo profile runs demultiplexing, reference alignment, gap/variant summaries, KO classification, and IGV report generation.

Run:

```bash
nextflow run . -profile docker,demo
```
