# KOCheck

A reproducible Nextflow pipeline for validating genetic knockouts. KOCheck analyzes Illumina sequencing reads to confirm gene deletions, detect correct resistance marker insertion, and identify ectopic integrations across the genome. Designed primarily for microbiological and fungal genetic workflows.

## Features

- **Quality Control**: Automatic adapter trimming and quality filtering with fastp
- **Deletion Detection**: Validates gene deletion by comparing target gene coverage to chromosomal mean
- **Marker Validation**: Confirms presence and copy number of resistance marker
- **Ectopic Integration Detection**: Identifies off-target marker insertions using junction support analysis
- **Coverage Visualization**: Generates coverage plots for the target gene region
- **Comprehensive Reporting**: CSV outputs for deletion status and marker presence
- **Summary Reports**: Aggregated CSV and interactive HTML report with all sample results and embedded coverage plots

## Requirements

- Nextflow (>= 22.10.0)
- Docker or Conda (for containerized tool execution)
- Java 11 or higher (for Nextflow)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/meadm/KOCheck
cd KOCheck
```

2. Ensure Nextflow is installed:
```bash
# Install Nextflow if needed
curl -s https://get.nextflow.io | bash
```

The pipeline uses Docker containers and Conda environments as specified in `conf/modules.config`. Make sure Docker is running or Conda is available.

## Usage

### Basic Usage

```bash
nextflow run main.nf \
  --reads "path/to/*_{R1,R2}.fq.gz" \
  --reference path/to/reference.fasta \
  --gene_bed path/to/target_gene.bed \
  --marker_fasta path/to/marker.fasta
```

### Example with Test Data

```bash
nextflow run main.nf \
  --reads "assets/testdata/*_{R1,R2}.fq.gz" \
  --reference assets/testdata/ref.fasta \
  --gene_bed assets/testdata/target_gene.bed \
  --marker_fasta assets/testdata/kanR.fasta \
  --marker_contig kanR
```

## Parameters

### Required Parameters

- `--reads`: Path pattern to paired-end FASTQ files (e.g., `"data/*_{R1,R2}.fq.gz"`)
- `--reference`: Path to reference genome FASTA file
- `--gene_bed`: Path to BED file with target gene coordinates (must contain exactly one line: `chrom start end`)
- `--marker_fasta`: Path to FASTA file containing the resistance marker sequence

### Optional Parameters

- `--marker_contig`: Name of the marker contig in the combined reference (default: `kanR`)
- `--outdir`: Output directory (default: `results`)
- `--delete_ratio`: Coverage ratio threshold for deletion classification (default: `0.1`)
- `--intact_ratio`: Coverage ratio threshold for intact gene classification (default: `0.5`)
- `--flank`: Flanking region size in bp for coverage plots and junction analysis (default: `500`)
- `--min_mapq`: Minimum mapping quality for junction support analysis (default: `20`)

## Workflow

The pipeline consists of the following steps:

1. **FASTP**: Quality control and adapter trimming of paired-end reads
2. **APPEND_MARKER**: Combines reference genome with marker sequence for marker detection
3. **BWA_MEM2**: Alignment of trimmed reads to the combined reference (alignment, sorting, indexing)
4. **MOSDEPTH**: Coverage calculation across the genome
5. **DELETION_CHECK**: Classifies target gene as `deleted`, `intact`, or `ambiguous` based on coverage ratios
6. **COVERAGE_PLOT**: Generates visualization of coverage across the target gene region
7. **MARKER_CHECK**: Validates marker presence, estimates copy number, and detects ectopic integrations
8. **AGGREGATE_RESULTS**: Combines all sample results into a single CSV summary and generates an interactive HTML report

## Output Files

Results are organized in the `results/` directory (or custom `--outdir`):

- `qc/`: Fastp quality control reports (HTML and JSON)
- `mapping/`: Aligned BAM files and indices
- `coverage/`: Mosdepth coverage summaries and per-base coverage files
- `deletion_check/`: CSV files with deletion status for each sample
- `marker_check/`: CSV files with marker validation results including:
  - Marker presence (true/false)
  - Deletion status
  - Marker copy number estimate (0, 1, or >1)
  - Junction support score
  - Final classification status
- `plots/`: Coverage plots for each sample showing the target gene region
- `summary/`: Aggregated results including:
  - `kocheck_summary.csv`: Combined CSV with all sample results
  - `kocheck_report.html`: Interactive HTML report with summary statistics, sample table, and embedded coverage plots

### Overall Status Classifications

The marker check module classifies samples into the following categories:

- `PASS`: Successful knockout - gene deleted, marker present in one copy, and strong junction support (>10)
- `FAIL`: Failed knockout - gene is still present (intact)
- `REVIEW`: Requires manual review - includes cases with:
  - Gene deleted but marker in multiple copies (including tandem insertions and/or ectopic integrations)
  - Gene deleted with single-copy marker but weak/moderate junction support
  - Other ambiguous situations where classification is unclear

## Configuration

Process-specific resource requirements and container configurations are defined in `conf/modules.config`. The pipeline supports both Docker and Conda environments.

## License

See LICENSE file for details.
