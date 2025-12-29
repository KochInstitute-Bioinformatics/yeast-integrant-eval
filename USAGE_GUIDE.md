# Usage Guide

Comprehensive guide for running the yeast-integrant-eval workflow with detailed examples and best practices.

## Table of Contents

- [Quick Reference](#quick-reference)
- [Single Sample Processing](#single-sample-processing)
- [Batch Processing](#batch-processing)
- [Configuring Parameters](#configuring-parameters)
- [Working with Transgenes](#working-with-transgenes)
- [Advanced Workflows](#advanced-workflows)
- [Best Practices](#best-practices)
- [Interpreting Results](#interpreting-results)

## Quick Reference

### Minimal Commands

```bash
# Single sample with defaults
nextflow run main.nf --input_fastq sample.fastq.gz --name sample_001

# Batch processing
nextflow run main.nf --samples batch.csv --outdir results

# With assembly evaluation
nextflow run main.nf --samples batch.csv --run_assembly_evaluation true --reference_genome ref.fa
```

## Single Sample Processing

### Basic Single Sample Run

Process one FASTQ file with default parameters:

```bash
nextflow run main.nf \
  --input_fastq /path/to/sample.fastq.gz \
  --name sample_001 \
  --outdir results
```

**What this does:**

- Filters reads by quality (Q≥10) and length (40kb+, 50kb+)
- Creates downsampled datasets at 25% and 50%
- Runs 1 assembly replicate per condition
- Detects default transgene (A-vector_herceptin_pEY345)
- Generates QC reports and summary statistics

### Specify Custom Transgene

```bash
nextflow run main.nf \
  --input_fastq /path/to/sample.fastq.gz \
  --name sample_001 \
  --transgene superfolder_GFP \
  --outdir results
```

**Note:** Transgene must exist in `transgenes/transgene_library.csv`

### Single Sample with Multiple Replicates

Generate multiple assembly replicates for statistical analysis:

```bash
nextflow run main.nf \
  --input_fastq /path/to/sample.fastq.gz \
  --name sample_001 \
  --replicates 5 \
  --outdir results
```

**What this does:**

- Creates 5 independent assemblies at 25% downsample rate
- Creates 5 independent assemblies at 50% downsample rate
- Enables calculation of confidence intervals

### Single Sample with Custom Genome Size

For non-yeast organisms or unusual genome sizes:

```bash
nextflow run main.nf \
  --input_fastq /path/to/sample.fastq.gz \
  --name sample_001 \
  --genome_size 15m \
  --outdir results
```

**Common genome sizes:**

- S. cerevisiae: `12m`
- P. pastoris: `9m` to `10m`
- Larger yeast: `15m` to `20m`

## Batch Processing

### Basic Batch Processing

Process multiple samples with a CSV configuration file:

```bash
nextflow run main.nf \
  --samples my_samples.csv \
  --outdir results
```

### CSV File Format

Create a semicolon-separated CSV file:

```csv
name;fastq;transgene;size_ranges;downsample_rates
sample_001;/data/sample_001.fastq.gz;A-vector_herceptin_pEY345;[{"min":40000,"name":"40k_Plus"},{"min":50000,"name":"50k_Plus"}];[0.25,0.5]
sample_002;/data/sample_002.fastq.gz;superfolder_GFP;[{"min":30000,"name":"30k_Plus"},{"min":40000,"name":"40k_Plus"}];[0.5,0.75]
sample_003;/data/sample_003.fastq.gz;E1_boNT;[{"min":20000,"name":"20k_Plus"},{"min":30000,"name":"30k_Plus"},{"min":40000,"name":"40k_Plus"}];[0.1,0.25,0.5]
```

**Key points:**

- Use **semicolons** (`;`) as separators, not commas
- Transgene must exist in transgene library
- `size_ranges` and `downsample_rates` are optional (will use defaults if omitted)

### Example 1: Uniform Parameters for All Samples

All samples use same settings (40kb+ and 50kb+ filtering, 25% and 50% downsampling):

```csv
name;fastq;transgene
sample_001;/data/sample_001.fastq.gz;A-vector_herceptin_pEY345
sample_002;/data/sample_002.fastq.gz;A-vector_herceptin_pEY345
sample_003;/data/sample_003.fastq.gz;A-vector_herceptin_pEY345
```

### Example 2: Per-Sample Custom Parameters

Each sample has custom read length thresholds and downsample rates:

```csv
name;fastq;transgene;size_ranges;downsample_rates
low_coverage;/data/low_cov.fastq.gz;transgene_A;[{"min":20000,"name":"20k_Plus"},{"min":30000,"name":"30k_Plus"}];[0.5,0.75,1.0]
high_coverage;/data/high_cov.fastq.gz;transgene_A;[{"min":50000,"name":"50k_Plus"},{"min":60000,"name":"60k_Plus"}];[0.1,0.25]
standard;/data/standard.fastq.gz;transgene_B;[{"min":40000,"name":"40k_Plus"}];[0.5]
```

**Use cases:**

- **low_coverage**: Includes shorter reads and more aggressive downsampling
- **high_coverage**: Uses only long reads and less downsampling
- **standard**: Single condition for quick testing

### Example 3: Mixed Transgenes

Different samples with different transgenes:

```csv
name;fastq;transgene;size_ranges;downsample_rates
GFP_sample_01;/data/GFP_01.fastq.gz;superfolder_GFP;[{"min":40000,"name":"40k_Plus"}];[0.5]
GFP_sample_02;/data/GFP_02.fastq.gz;superfolder_GFP;[{"min":40000,"name":"40k_Plus"}];[0.5]
herceptin_01;/data/herceptin_01.fastq.gz;A-vector_herceptin_pEY345;[{"min":40000,"name":"40k_Plus"}];[0.5]
boNT_sample_01;/data/boNT_01.fastq.gz;E1_boNT;[{"min":30000,"name":"30k_Plus"}];[0.25,0.5]
```

## Configuring Parameters

### Read Length Thresholds (size_ranges)

Control which read length cutoffs to test:

```json
[
  {"min":20000,"name":"20k_Plus"},
  {"min":30000,"name":"30k_Plus"},
  {"min":40000,"name":"40k_Plus"},
  {"min":50000,"name":"50k_Plus"}
]
```

**Choosing thresholds:**

- **Lower thresholds (20-30kb)**: Retain more data, lower assembly quality
- **Medium thresholds (40kb)**: Balanced (recommended for most cases)
- **Higher thresholds (50-60kb)**: Fewer reads, higher assembly quality

**When to use each:**

| Read Length | Typical Use Case | Coverage Needed |
| ------------- | ------------------ | ----------------- |
| 20kb+ | Low-coverage samples, exploratory | 30-50x |
| 30kb+ | Standard yeast assembly | 40-80x |
| 40kb+ | High-quality assembly (recommended) | 50-100x |
| 50kb+ | Ultra-high quality | 80-150x |

### Downsample Rates

Control what fractions of reads to use for assembly:

```json
[0.1, 0.25, 0.5, 0.75, 1.0]
```

**Values:**

- `0.1` = 10% of reads
- `0.25` = 25% of reads
- `0.5` = 50% of reads
- `0.75` = 75% of reads
- `1.0` = 100% of reads (no downsampling)

**Strategy recommendations:**

**Coverage Optimization Study:**

```json
[0.1, 0.25, 0.5, 0.75, 1.0]
```

Test wide range to find minimum required coverage.

**Production Assemblies:**

```json
[0.75, 1.0]
```

Focus on high-coverage assemblies only.

**Rapid Screening:**

```json
[0.25]
```

Quick single-point assembly for initial assessment.

**Statistical Robustness:**

```json
[0.5, 0.75]
```

Combine with `--replicates 10` for confidence intervals.

### Assembly Depth Thresholds

Control which samples proceed to assembly based on estimated coverage:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --min_assembly_depth 20 \
  --max_assembly_depth 100 \
  --outdir results
```

See [QUICK_START_DEPTH_CONFIG.md](QUICK_START_DEPTH_CONFIG.md) for detailed guidance.

## Working with Transgenes

### Using Existing Transgenes

List available transgenes:

```bash
cat transgenes/transgene_library.csv
```

Example output:

```csv
transgene_name,fasta_file
A-vector_herceptin_pEY345,A-vector_herceptin_pEY345.fasta
superfolder_GFP,superfolder_GFP.fasta
E1_boNT,E1_boNT.fasta
```

Use in your samples CSV:

```csv
name;fastq;transgene
sample_001;/data/sample.fastq.gz;superfolder_GFP
```

### Adding New Transgenes

#### Step 1: Create FASTA File

Create a FASTA file with your transgene sequence:

```bash
cat > transgenes/my_new_transgene.fasta << 'EOF'
>my_new_transgene description
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
EOF
```

#### Step 2: Add to Library CSV

Add entry to `transgenes/transgene_library.csv`:

```bash
echo "my_new_transgene,my_new_transgene.fasta" >> transgenes/transgene_library.csv
```

#### Step 3: Use in Workflow

Reference in your samples CSV:

```csv
name;fastq;transgene
sample_001;/data/sample.fastq.gz;my_new_transgene
```

### Multiple Transgenes Per Sample

To detect multiple transgenes in a single sample, run the workflow separately for each transgene:

```csv
# First CSV: sample_GFP.csv
name;fastq;transgene
sample_001;/data/sample_001.fastq.gz;superfolder_GFP

# Second CSV: sample_mCherry.csv
name;fastq;transgene
sample_001;/data/sample_001.fastq.gz;mCherry
```

Then run:

```bash
nextflow run main.nf --samples sample_GFP.csv --outdir results_GFP
nextflow run main.nf --samples sample_mCherry.csv --outdir results_mCherry
```

## Advanced Workflows

### Replicate-Based Statistical Analysis

Generate confidence intervals for assembly metrics:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --replicates 10 \
  --outdir results_statistical
```

**What this does:**

- Creates 10 independent assemblies per condition
- Each uses different random seed for downsampling
- Enables calculation of mean, standard deviation, confidence intervals

**Analysis example:**

After workflow completes, analyze results:

```bash
# Extract assembly statistics for all replicates
cat results_statistical/summary/assembly_stats_summary.csv

# Calculate mean and std dev for sample_001 at 50% downsample
grep "sample_001.*ds0.5" results_statistical/summary/assembly_stats_summary.csv | \
  awk -F',' '{sum+=$5; sumsq+=$5*$5} END {print "Mean:", sum/NR, "StdDev:", sqrt(sumsq/NR - (sum/NR)^2)}'
```

### Assembly Evaluation with Reference

Enable comprehensive assembly validation:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --run_assembly_evaluation true \
  --reference_genome /path/to/WT_reference.fa \
  --transcripts_fasta /path/to/transcripts.fa \
  --outdir results_evaluated
```

**Outputs include:**

- Reference-aligned assemblies
- Oriented contigs (corrected for strand)
- Final merged chromosomes
- Read mapping (BAM files for IGV)
- Transgene locations (BED files)
- Transcript mapping (BED files)

See [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md) for complete documentation.

### Combined Statistical + Evaluation Workflow

For maximum confidence:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --replicates 5 \
  --run_assembly_evaluation true \
  --reference_genome /path/to/WT_reference.fa \
  --min_assembly_depth 30 \
  --max_assembly_depth 150 \
  --outdir results_comprehensive
```

**Resource considerations:**

- High CPU/memory usage
- Long runtime (multiply by replicate count)
- Large storage requirements

### Parameter Sweep for Optimization

Test multiple parameter combinations systematically:

**Create CSV with graduated parameters:**

```csv
name;fastq;transgene;size_ranges;downsample_rates
sample_sweep_20k;/data/sample.fastq.gz;transgene;[{"min":20000,"name":"20k_Plus"}];[0.25,0.5,0.75,1.0]
sample_sweep_30k;/data/sample.fastq.gz;transgene;[{"min":30000,"name":"30k_Plus"}];[0.25,0.5,0.75,1.0]
sample_sweep_40k;/data/sample.fastq.gz;transgene;[{"min":40000,"name":"40k_Plus"}];[0.25,0.5,0.75,1.0]
sample_sweep_50k;/data/sample.fastq.gz;transgene;[{"min":50000,"name":"50k_Plus"}];[0.25,0.5,0.75,1.0]
```

This creates 16 assemblies (4 read lengths × 4 coverage levels) from one sample.

## Best Practices

### Workflow Design

✅ **DO:**

- Start with 1-2 samples to validate parameters
- Use default parameters first, then optimize
- Include positive controls (samples with known results)
- Document your parameter choices
- Archive raw FASTQ files separately

❌ **DON'T:**

- Run hundreds of samples without testing
- Use extremely aggressive filtering without justification
- Ignore quality control reports
- Delete intermediate files before validation

### Resource Management

**For typical yeast genome (10-15 Mb):**

- Allocate 32+ GB RAM for Flye
- Use 8-24 CPUs for assembly
- Budget 100-150 GB storage per sample

**For batch processing:**

- Nextflow automatically parallelizes
- Limit concurrent jobs with executor settings
- Monitor with `execution_report.html`

### Quality Control Checkpoints

1. **After Filtering** - Review NanoPlot reports
   - Check read length distributions
   - Verify quality scores
   - Assess total data retained

2. **After Preflight** - Review assembly candidates
   - Check `assembly_candidates.csv`
   - Verify coverage estimates reasonable
   - Adjust thresholds if needed

3. **After Assembly** - Review Flye outputs
   - Check assembly sizes (expect ~12 Mb for yeast)
   - Verify contig counts (fewer is better)
   - Look for circular contigs

4. **After BLAST** - Review transgene detection
   - Check copy number estimates
   - Verify coverage of transgene
   - Look for multiple integration sites

### Organizing Large Projects

**Directory structure:**

```text
project/
├── raw_data/
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   └── ...
├── configs/
│   ├── batch_exploratory.csv
│   ├── batch_production.csv
│   └── custom.config
├── results/
│   ├── run_20231201_exploratory/
│   ├── run_20231215_production/
│   └── run_20240110_replicates/
└── analysis/
    ├── notebooks/
    └── scripts/
```

**Sample CSV organization:**

```bash
# Exploratory phase
batch_exploratory.csv      # All samples, permissive parameters

# Production phase  
batch_production.csv       # Selected samples, optimized parameters

# Statistical analysis
batch_statistical.csv      # Final samples with replicates
```

## Interpreting Results

### Key Output Files

#### Assembly Statistics Summary

**File:** `results/summary/assembly_stats_summary.csv`

**Key columns:**

- `sample_name`: Sample identifier with conditions
- `total_length`: Total assembly size (expect ~12 Mb for yeast)
- `num_contigs`: Number of contigs (fewer = better)
- `N50`: Assembly contiguity metric (higher = better)
- `circular_contigs`: Number of circular contigs (chromosomes/plasmids)

**Good assembly indicators:**

- Total length: 10-15 Mb (for yeast)
- Num contigs: <50 (ideally <20)
- N50: >500 kb (ideally >1 Mb)
- Circular contigs: ≥1

#### BLAST Results Summary

**File:** `results/summary/blast_results_summary.csv`

**Key columns:**

- `sample_name`: Sample identifier
- `transgene`: Transgene detected
- `copy_number_estimate`: Estimated transgene copies
- `total_coverage`: Total bases aligned
- `num_hits`: Number of BLAST alignments

**Interpreting copy numbers:**

- 1.0: Single integration
- 2.0: Two copies (tandem or separate)
- 0.5: Partial integration or contamination
- >3.0: Multiple integrations or amplification

#### Quality Control Reports

**Files:** `results/nanoplot/*/NanoPlot-report.html`

**Key metrics:**

- Mean read length: Should match your filtering threshold
- Mean read quality: Should be >10 (typically 12-15)
- Total bases: Indicates sequencing depth
- Read length distribution: Should show clear cutoff

### Comparing Conditions

**Compare read length thresholds:**

```bash
# Extract assembly sizes for different size filters
grep "sample_001" results/summary/assembly_stats_summary.csv | \
  awk -F',' '{print $1, $3}' | \
  sort
```

**Compare downsample rates:**

```bash
# Extract N50 values for different coverage levels
grep "sample_001" results/summary/assembly_stats_summary.csv | \
  awk -F',' '{print $1, $6}' | \
  grep "ds" | \
  sort
```

### Red Flags

⚠️ **Assembly Issues:**

- Total length >>15 Mb or <<10 Mb
- >100 contigs
- N50 <100 kb
- No circular contigs

⚠️ **Transgene Detection Issues:**

- Copy number <0.3 (possible detection failure)
- Copy number >10 (possible artifact)
- No BLAST hits (transgene not present or wrong sequence)

⚠️ **QC Issues:**

- Mean quality <8
- Very low read counts after filtering
- Bimodal read length distribution

## Troubleshooting Common Issues

See main [README.md](README.md#troubleshooting) for detailed troubleshooting guide.

## Additional Resources

- **Main Documentation:** [README.md](README.md)
- **Output Structure:** [OUTPUT_STRUCTURE.md](OUTPUT_STRUCTURE.md)
- **Assembly Evaluation:** [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md)
- **Depth Configuration:** [QUICK_START_DEPTH_CONFIG.md](QUICK_START_DEPTH_CONFIG.md)
- **Nextflow Documentation:** [https://www.nextflow.io/docs/latest/](https://www.nextflow.io/docs/latest/)
- **Flye Documentation:** [https://github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)
