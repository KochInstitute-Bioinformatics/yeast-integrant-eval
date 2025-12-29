# Quick Start: Configuring Assembly Depth Thresholds

## Overview

The workflow includes a **Flye Preflight** step that estimates sequencing coverage before running full assemblies. This prevents wasting computational resources on samples with insufficient or excessive coverage.

Assembly depth thresholds determine which samples proceed to assembly:

- **`min_assembly_depth`**: Minimum coverage required (default: 20x)
- **`max_assembly_depth`**: Maximum coverage allowed (default: 100x)

## TL;DR

Add these two lines to your `nextflow.config` under the `params` section:

```groovy
params {
    min_assembly_depth = 20    // Minimum coverage for assembly (default: 20)
    max_assembly_depth = 100   // Maximum coverage for assembly (default: 100)
}
```

## Three Ways to Configure Thresholds

### 1. Edit Configuration File (Permanent)

Edit `nextflow.config` to set new defaults:

```groovy
params {
    min_assembly_depth = 10
    max_assembly_depth = 200
}
```

Run normally:

```bash
nextflow run main.nf --samples batch.csv --outdir results
```

### 2. Command Line Override (One-time)

Keep your config as-is, override when running:

```bash
nextflow run main.nf \
    --samples batch.csv \
    --min_assembly_depth 10 \
    --max_assembly_depth 200 \
    --outdir results
```

### 3. Update Execution Script (Convenient for Repeated Runs)

Add parameters to your workflow execution script:

```bash
#!/bin/bash
# run_workflow.sh

nextflow run main.nf \
    --samples batch.csv \
    --min_assembly_depth 10 \
    --max_assembly_depth 200 \
    --outdir results
```

Then run:

```bash
bash run_workflow.sh
```

## Common Use Cases

### Include Low-Coverage Samples

For exploratory analysis or when working with limited sequencing data:

```bash
--min_assembly_depth 5 --max_assembly_depth 300
```

### Handle High-Coverage Samples

For deeply sequenced samples or organisms with smaller genomes:

```bash
--min_assembly_depth 10 --max_assembly_depth 500
```

### No Upper Limit

Remove upper threshold to assemble all samples regardless of coverage:

```bash
--min_assembly_depth 10 --max_assembly_depth 999999
```

### Only High-Quality Assemblies

Strict thresholds for production-quality assemblies:

```bash
--min_assembly_depth 50 --max_assembly_depth 150
```

### Broad Range (Recommended for Initial Exploration)

Accept a wide range of coverages to see what your data looks like:

```bash
--min_assembly_depth 10 --max_assembly_depth 300
```

## Understanding the Output

After running the Flye Preflight step, check these files to see which samples were selected or filtered:

### Assembly Candidates (Passed Filtering)

**File:** `results/flye_preflight/assembly_candidates.csv`

Contains samples that passed depth thresholds and will be assembled:

```csv
sample_name,fastq_path,estimated_depth,status
sample_001_40k_Plus,/path/to/filtered.fastq.gz,75.3,candidate
sample_002_50k_Plus_ds0.5,/path/to/downsampled.fastq.gz,45.2,candidate
```

### Assembly Filtered (Did Not Pass)

**File:** `results/flye_preflight/assembly_filtered.csv`

Contains samples filtered out with reasons:

```csv
sample_name,fastq_path,estimated_depth,status,reason
sample_003_30k_Plus,/path/to/filtered.fastq.gz,15.2,filtered,Below minimum depth (15.2 < 20)
sample_004_40k_Plus,/path/to/filtered.fastq.gz,250.8,filtered,Above maximum depth (250.8 > 100)
```

## Choosing Appropriate Thresholds

### Consider Your Organism

| Organism Type | Genome Size | Recommended Min | Recommended Max |
| --------------- | ------------- | ----------------- | ----------------- |
| Yeast (S. cerevisiae) | ~12 Mb | 20-30x | 100-150x |
| Yeast (P. pastoris) | ~9-10 Mb | 20-30x | 100-150x |
| Small plasmids | <100 kb | 50-100x | 500-1000x |

### Consider Your Goals

**Exploratory Analysis:**

- Lower minimum (10-15x) to include more samples
- Higher maximum (200-300x) to see full range

**Production Assemblies:**

- Higher minimum (30-50x) for quality
- Moderate maximum (100-150x) to avoid over-coverage issues

**Coverage Optimization Study:**

- Very low minimum (5x) to include all samples
- No maximum (999999) to see full spectrum

### Flye Coverage Recommendations

From Flye documentation:

- **Minimum:** 10-20x coverage typically required
- **Optimal:** 30-100x coverage for best results
- **Maximum:** >200x may cause memory issues or assembly artifacts

## Troubleshooting

### All Samples Filtered Out

**Symptom:** `assembly_candidates.csv` is empty

**Causes:**

1. Thresholds too strict for your data
2. Read filtering removed too many reads
3. Downsampling rates too aggressive

**Solutions:**

- Lower `min_assembly_depth` (try 10 or 15)
- Increase `max_assembly_depth` (try 200 or 300)
- Review NanoPlot reports to check actual coverage
- Reduce read length filtering thresholds
- Increase downsample rates (use 0.75 or 1.0)

### No Samples Filtered Out

**Symptom:** All samples in `assembly_candidates.csv`, none filtered

**Causes:**

1. Thresholds too permissive
2. All samples happen to be in range

**Solutions:**

- If intentional, no action needed
- If trying to filter, adjust thresholds more strictly
- Review `estimated_depth` column to see actual coverage distribution

### Inconsistent Assembly Quality

**Symptom:** Some assemblies are fragmented or poor quality

**Causes:**

1. Including samples near minimum threshold
2. Coverage estimation slightly off

**Solutions:**

- Increase `min_assembly_depth` by 5-10x
- Review individual sample NanoPlot reports
- Consider higher read length thresholds

## Advanced: Per-Sample Depth Requirements

While the workflow doesn't support per-sample depth thresholds directly, you can achieve this by:

### Option 1: Multiple Workflow Runs

Run workflow separately for different sample groups:

```bash
# High-quality samples
nextflow run main.nf \
    --samples high_quality_samples.csv \
    --min_assembly_depth 30 \
    --max_assembly_depth 150

# Exploratory samples
nextflow run main.nf \
    --samples exploratory_samples.csv \
    --min_assembly_depth 10 \
    --max_assembly_depth 300
```

### Option 2: Pre-filter Sample CSV

Run preflight once, examine results, then create filtered CSV with only desired samples:

```bash
# First run to get coverage estimates
nextflow run main.nf --samples all_samples.csv --outdir initial_run

# Examine results/flye_preflight/assembly_candidates.csv
# Create new CSV with selected samples

# Re-run with selected samples only
nextflow run main.nf --samples selected_samples.csv --outdir final_run
```

## Related Parameters

These parameters also affect which assemblies are created:

### Read Filtering

- `--min_quality`: Minimum read quality score (default: 10)
- `size_ranges` in CSV: Read length thresholds to test

### Downsampling

- `downsample_rates` in CSV: Coverage fractions to test
- `--replicates`: Number of assembly replicates per condition

See the main [README.md](README.md) for complete parameter documentation.

## Additional Resources

- **Main Documentation:** [README.md](README.md)
- **Assembly Evaluation:** [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md)
- **Flye Documentation:** [https://github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)
