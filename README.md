# Yeast Integrant Evaluation

A Nextflow workflow for assembling and evaluating integrant-containing yeast genomes using Oxford Nanopore Technologies (ONT) long-read sequencing data. This workflow performs comprehensive assembly quality assessment, transgene detection, and copy number evaluation.

## Overview

This workflow processes ONT long-read sequencing data through a pipeline that:

1. **Quality Filters Reads** - Removes low-quality reads (Q<10) using Chopper
2. **Applies Read Length Thresholds** - Tests multiple length cutoffs to identify optimal parameters
3. **Downsamples Dnput Data** - Ensure that depth of coverage is within usable range
4. **Creates Multiple Assembly Replicates** - Generates assemblies from downsampled read sets at various coverage levels
5. **Performs Assembly QC** - Comprehensive quality control with NanoPlot and Flye statistics
6. **Detects Transgenes** - BLAST-based identification and copy number estimation
7. **Evaluates Assemblies** - Reference-based validation and visualization

## Key Features

### **Flexible Read Length Testing**

Configure custom read length thresholds per sample to identify optimal parameters for your dataset. Default: 40kb+ and 50kb+ cutoffs.

### **Multi-Coverage Assembly Strategy**

Generate assemblies from multiple downsampled datasets to:

- Assess assembly quality vs. sequencing depth
- Identify minimum coverage requirements
- Evaluate assembly consistency across replicates

### **Transgene Detection & Quantification**

Automated BLAST-based transgene detection with:

- Copy number estimation
- Integration site identification
- Multi-transgene support via transgene library

### **Replicate-Based Robustness Testing**

Create multiple assembly replicates (default: 1 per downsample rate, configurable up to N replicates) to:

- Assess assembly reproducibility
- Calculate confidence intervals
- Identify spurious assembly artifacts

### **Assembly Evaluation Module**

Comprehensive assembly validation including:

- Reference genome alignment
- Contig orientation correction
- Chromosome coverage analysis
- Read mapping and visualization
- Transcript mapping (optional)

## Requirements

### Software Dependencies

- **Nextflow** (≥22.10.0) - [Installation guide](https://www.nextflow.io/docs/latest/getstarted.html)
- **Singularity** or **Docker** (for containerized execution)
  - Singularity recommended for HPC environments
  - All tools run in pre-built containers (no manual installation needed)

### Computational Resources

**Recommended minimum:**

- **CPU**: 8+ cores (24+ cores for optimal Flye performance)
- **Memory**: 32+ GB RAM
- **Storage**: ~10-15x input file size for intermediate files

**Example resource needs** for a 10 GB FASTQ file:

- Storage: ~100-150 GB free space
- Runtime: 4-12 hours (depends on coverage and replicate count)

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/KochInstitute-Bioinformatics/yeast-integrant-eval.git
cd yeast-integrant-eval
```

### 2. Verify Nextflow Installation

```bash
nextflow -version
# Should show version 22.10.0 or higher
```

### 3. Test with Example Data (Optional)

```bash
# Quick test with provided example configuration
nextflow run main.nf --samples samples_example.csv --outdir test_results
```

## Quick Start

### Single Sample Analysis

Process a single FASTQ file with default parameters:

```bash
nextflow run main.nf \
  --input_fastq sample.fastq.gz \
  --name sample_001 \
  --outdir results
```

### Batch Processing (Recommended)

Process multiple samples using a CSV configuration file:

```bash
nextflow run main.nf \
  --samples batch_samples.csv \
  --outdir results
```

## Input Configuration

### CSV File Format

The workflow uses a **semicolon-separated** CSV file with the following columns:

```csv
name;fastq;transgene;size_ranges;downsample_rates
sample_001;/path/to/sample_001.fastq.gz;transgene_name;[{"min":40000,"name":"40k_Plus"},{"min":50000,"name":"50k_Plus"}];[0.25,0.5]
sample_002;/path/to/sample_002.fastq.gz;transgene_name;[{"min":30000,"name":"30k_Plus"},{"min":40000,"name":"40k_Plus"}];[0.5,0.75]
```

### Column Descriptions

| Column | Description | Example |
| -------- | ------------- | --------- |
| `name` | Unique sample identifier | `sample_001` |
| `fastq` | Full path to FASTQ file (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`) | `/data/sample_001.fastq.gz` |
| `transgene` | Transgene name (must exist in `transgenes/transgene_library.csv`) | `A-vector_herceptin_pEY345` |
| `size_ranges` | JSON array of read length thresholds to test | `[{"min":40000,"name":"40k_Plus"}]` |
| `downsample_rates` | JSON array of downsample fractions | `[0.25,0.5,0.75]` |

### Configuring Read Length Thresholds

The `size_ranges` parameter specifies which read length cutoffs to test. Each entry creates a separate filtered dataset:

```json
[
  {"min":30000,"name":"30k_Plus"},
  {"min":40000,"name":"40k_Plus"},
  {"min":50000,"name":"50k_Plus"}
]
```

- **`min`**: Minimum read length in base pairs
- **`name`**: Descriptive name (used in output filenames)

**Default if not specified:** `40k_Plus` and `50k_Plus`

### Configuring Downsample Rates

The `downsample_rates` parameter specifies what fractions of the filtered reads to use for assembly:

```json
[0.25, 0.5, 0.75]
```

This will create assemblies using:

- 25% of filtered reads
- 50% of filtered reads  
- 75% of filtered reads

**Default if not specified:** `[0.25, 0.5]`

### Complete Example CSV

See `samples_example.csv` for a working example:

```csv
name;fastq;transgene;size_ranges;downsample_rates
S-1077-1;/net/bmc-lab3/data/bcc/projects/cmelen-Love/072423_ONT/230718Lov/data/all_S-1077-1.fastq;superfolder_GFP;[{"min":30000,"name":"30k_Plus"},{"min":40000,"name":"40k_Plus"},{"min":50000,"name":"50k_Plus"}];[0.25,0.5]
S-1030_3;/net/bmc-lab3/data/bcc/projects/cmelen-Love/091123_ONT/230828LovA/all_S-1030_3.fastq;E1_boNT;[{"min":20000,"name":"20k_Plus"},{"min":30000,"name":"30k_Plus"},{"min":40000,"name":"40k_Plus"},{"min":50000,"name":"50k_Plus"}];[0.1,0.25,0.5]
```

## Parameters Reference

### Required Parameters (Single Sample Mode)

| Parameter | Description | Example |
| ----------- | ------------- | --------- |
| `--input_fastq` | Path to input FASTQ file | `sample.fastq.gz` |
| `--name` | Sample name | `sample_001` |

### Required Parameters (Batch Mode)

| Parameter | Description | Example |
| ----------- | ------------- | --------- |
| `--samples` | Path to CSV configuration file | `batch_samples.csv` |

### Optional Global Parameters

| Parameter | Default | Description |
| ----------- | --------- | ------------- |
| `--outdir` | `results` | Output directory |
| `--genome_size` | `10m` | Estimated genome size for Flye (e.g., `10m`, `12m`, `15m`) |
| `--min_quality` | `10` | Minimum read quality score for filtering |
| `--replicates` | `1` | Number of assembly replicates per downsample rate |
| `--default_transgene` | `A-vector_herceptin_pEY345` | Default transgene if not specified in CSV |

### Assembly Evaluation Parameters (Advanced)

| Parameter | Default | Description |
| ----------- | --------- | ------------- |
| `--run_assembly_evaluation` | `false` | Enable assembly evaluation module |
| `--reference_genome` | `null` | Path to wild-type reference genome (FASTA) |
| `--reference_chromosomes` | `null` | Optional: Chromosome-only reference for coverage |
| `--transcripts_fasta` | `null` | Optional: Transcripts for mapping (FASTA) |
| `--min_assembly_depth` | `20` | Minimum depth threshold for assembly candidates |
| `--max_assembly_depth` | `100` | Maximum depth threshold for assembly candidates |

## Output Structure

```text
results/
├── pipeline_info/                  # Execution reports and logs
│   ├── execution_timeline.html     # Visual timeline of process execution
│   ├── execution_report.html       # Detailed resource usage report
│   ├── execution_trace.txt         # Machine-readable execution trace
│   └── pipeline_dag.svg            # Pipeline DAG visualization
│
├── chopper/                        # Quality and size-filtered reads
│   └── <sample>_<size_range>/
│       └── filtered_reads.fastq.gz
│
├── downsample/                     # Downsampled read sets
│   └── <sample>_<size>_ds<rate>_rep<N>/
│       └── downsampled_reads.fastq.gz
│
├── nanoplot/                       # Quality control reports
│   ├── <sample>_original/          # Original FASTQ QC
│   │   └── NanoPlot-report.html
│   └── <sample>_<size>/            # Filtered FASTQ QC
│       └── NanoPlot-report.html
│
├── flye_preflight/                 # Assembly candidate evaluation
│   ├── assembly_candidates.csv     # Samples passing depth thresholds
│   └── assembly_filtered.csv       # Samples filtered out (with reasons)
│
├── flye/                           # Assembly outputs
│   └── <sample>_<condition>/
│       ├── assembly.fasta          # Final assembly
│       ├── assembly_info.txt       # Contig statistics
│       └── flye.log                # Assembly log
│
├── blast/                          # Transgene detection results
│   └── <sample>_<transgene>/
│       ├── blast_results.txt       # Raw BLAST output
│       └── blast_summary.csv       # Parsed copy number estimates
│
├── summary/                        # Aggregated results
│   ├── assembly_stats_summary.csv  # All assembly statistics
│   └── blast_results_summary.csv   # All transgene detection results
│
└── assembly_evaluation/            # Optional evaluation module outputs
    └── <sample>/
        ├── alignments/             # Reference alignments
        ├── repaired/               # Oriented contigs
        ├── final/                  # Final merged assembly
        ├── read_mapping/           # BAM files for visualization
        ├── transgene_blast/        # Transgene locations (BED format)
        └── transcript_mapping/     # Optional transcript locations
```

### Key Output Files

#### Assembly Results

- **`flye/<sample>/assembly.fasta`** - Final assembled genome
- **`flye/<sample>/assembly_info.txt`** - Contig lengths, coverage, circularity
- **`summary/assembly_stats_summary.csv`** - Compiled statistics for all assemblies

#### Quality Control

- **`nanoplot/<sample>/NanoPlot-report.html`** - Interactive QC report with read length/quality distributions
- **`flye_preflight/assembly_candidates.csv`** - Which samples proceed to assembly
- **`flye_preflight/assembly_filtered.csv`** - Samples filtered out (low/high coverage)

#### Transgene Analysis

- **`blast/<sample>/blast_results.txt`** - Raw BLAST alignments
- **`blast/<sample>/blast_summary.csv`** - Copy number estimates and statistics
- **`summary/blast_results_summary.csv`** - All transgene results compiled

## Transgene Library

The workflow includes a transgene library for automated detection. To add new transgenes:

### 1. Add FASTA File

Place your transgene sequence in `transgenes/`:

```bash
cp my_transgene.fasta transgenes/
```

### 2. Update Library CSV

Add an entry to `transgenes/transgene_library.csv`:

```csv
transgene_name,fasta_file
my_transgene,my_transgene.fasta
```

### 3. Reference in Samples CSV

Use the transgene name in your samples CSV:

```csv
name;fastq;transgene;size_ranges;downsample_rates
sample_001;/path/to/data.fastq.gz;my_transgene;[{"min":40000,"name":"40k_Plus"}];[0.5]
```

## Advanced Usage

### Running with Multiple Replicates

Generate statistical confidence by creating multiple assembly replicates:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --replicates 10 \
  --outdir results_replicate_analysis
```

This creates 10 independent assemblies per downsample rate, enabling:

- Calculation of assembly metric confidence intervals
- Identification of spurious assembly features
- Robust transgene copy number estimation

### Enabling Assembly Evaluation

For comprehensive assembly validation against a reference:

```bash
nextflow run main.nf \
  --samples batch.csv \
  --run_assembly_evaluation true \
  --reference_genome /path/to/WT_reference.fa \
  --transcripts_fasta /path/to/transcripts.fa \
  --outdir results_evaluated
```

See [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md) for detailed documentation.

### Custom Resource Configuration

Edit `nextflow.config` to adjust resources for your system:

```groovy
process {
    withName: 'FLYE' {
        cpus = 24          // Increase for faster assembly
        memory = '64 GB'   // Increase for large genomes
        time = '48.h'      // Adjust based on dataset size
    }
}
```

### HPC Execution (SLURM Example)

Create a custom configuration file `my_cluster.config`:

```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    clusterOptions = '--account=my_account'
}
```

Run with:

```bash
nextflow run main.nf \
  --samples batch.csv \
  -c my_cluster.config
```

## Workflow Details

### Phase 1: Quality Filtering (CHOPPER)

- Filters reads by quality score (default: Q≥10)
- Applies per-sample read length thresholds
- Generates separate filtered datasets for each size cutoff

### Phase 2: Downsampling

- Creates random subsets of filtered reads at specified fractions
- Generates N replicates per downsample rate (configurable)
- Uses reproducible random seeds for each replicate

### Phase 3: Assembly Candidate Evaluation (Flye Preflight)

- Runs fast Flye preflight check to estimate coverage
- Filters assemblies by coverage depth thresholds
- Prevents resource waste on low-quality or over-sequenced samples

### Phase 4: Assembly (Flye)

- De novo assembly of candidate datasets using Flye
- Optimized for long-read yeast genome assembly
- Automatic circularization detection

### Phase 5: Transgene Detection (BLAST)

- BLAST search against transgene library
- Copy number estimation from alignment coverage
- Multiple hit detection for complex integrations

### Phase 6: Quality Control & Summarization

- NanoPlot reports for read quality assessment
- Assembly statistics compilation
- Integrated summary tables for downstream analysis

### Phase 7: Assembly Evaluation (Optional)

See [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md) for details on:

- Reference-based contig orientation
- Chromosome coverage analysis
- Read mapping and visualization
- Transcript mapping

## Troubleshooting

### Common Issues

#### Assembly Failures

**Symptom:** Flye process fails or produces no output

**Solutions:**

- Check read quality and quantity with NanoPlot reports
- Verify `--genome_size` parameter matches your organism (yeast ≈ 10-15m)
- Increase memory allocation in `nextflow.config`
- Review `flye_preflight/assembly_filtered.csv` for coverage issues

#### Missing Transgene Files

**Symptom:** Error about missing transgene FASTA

**Solutions:**

- Verify transgene name in CSV matches entry in `transgenes/transgene_library.csv`
- Check that FASTA file exists in `transgenes/` directory
- Ensure no typos in transgene names (case-sensitive)

#### Low Assembly Quality

**Symptom:** Fragmented assemblies or low BLAST scores

**Solutions:**

- Increase minimum read length threshold (try 50kb+ cutoff)
- Use higher downsample rates (0.75-1.0) for better coverage
- Check NanoPlot reports for read quality issues
- Verify sufficient sequencing depth (aim for 50-100x)

#### Resource Errors (Out of Memory)

**Symptom:** Process killed due to memory limits

**Solutions:**

- Increase memory allocation in `nextflow.config` for specific processes
- Reduce `--genome_size` if over-estimated
- Process fewer samples simultaneously (reduce parallelization)

### Getting Help

- **Documentation:**

  - [USAGE_GUIDE.md](USAGE_GUIDE.md) - Complete parameter reference and workflow details
  - [OUTPUT_STRUCTURE.md](OUTPUT_STRUCTURE.md) - Understanding all output files
  - [ASSEMBLY_EVALUATION.md](ASSEMBLY_EVALUATION.md) - Reference-based validation
  - [QUICK_START_DEPTH_CONFIG.md](QUICK_START_DEPTH_CONFIG.md) - Depth filtering configuration
- **Issues:** Open an issue on [GitHub Issues](https://github.com/KochInstitute-Bioinformatics/yeast-integrant-eval/issues)
- **Contact:** Charlie Whittaker ([charliew@mit.edu](mailto:charliew@mit.edu))

### Display Help Message

```bash
nextflow run main.nf --help
```

## Performance Tips

### Optimizing Runtime

1. **Use appropriate read length cutoffs** - Higher thresholds (50kb+) reduce data volume and assembly time
2. **Limit downsample rates** - Focus on 1-2 key coverage levels instead of many
3. **Reduce replicates** - Start with `--replicates 1` for initial exploration
4. **Increase CPU allocation** - Flye scales well with more cores (up to 24 cpus)

### Managing Storage

1. **Use compressed FASTQ** - Input files as `.fastq.gz` save space
2. **Clean intermediate files** - Delete `downsample/` and `chopper/` after successful runs
3. **Archive old results** - Compress result directories after analysis

### Batch Processing Strategy

For large cohorts:

1. **Test with subset** - Run 2-3 samples first to validate parameters
2. **Use conservative thresholds** - Start with default `min_assembly_depth=20`
3. **Parallelize by compute** - Let Nextflow handle parallelization automatically
4. **Monitor resources** - Use `execution_report.html` to identify bottlenecks

## Citation

If you use this workflow in your research, please cite:

- **Nextflow:** Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4):316-319.
- **Flye:** Kolmogorov, M., et al. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37:540-546.
- **Minimap2:** Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18):3094-3100.
- **NanoPlot:** De Coster, W., et al. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15):2666-2669.
- **BLAST:** Camacho, C., et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10:421.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Support

For questions, issues, or feature requests:

- **GitHub Issues:** [yeast-integrant-eval/issues](https://github.com/KochInstitute-Bioinformatics/yeast-integrant-eval/issues)
- **Email:** Charlie Whittaker - [charliew@mit.edu](mailto:charliew@mit.edu)

## Acknowledgments

Developed at the Koch Institute for Integrative Cancer Research, MIT.
