# KOCheck Wrapper Scripts

This directory contains user-friendly wrapper scripts for the KOCheck pipeline, designed for biologists without command-line experience.

## Prerequisites

Before using either wrapper script, you need to have these installed:

### Required
- **Nextflow** - The workflow manager
  - Install: https://www.nextflow.io/
  - Quick install: `curl -s https://get.nextflow.io | bash`

### Required (choose one)
- **Docker** - For containerized tool execution (recommended)
  - Install: https://www.docker.com/get-started
  - Must be running when you execute the pipeline
  
  **OR**
  
- **Conda** - Alternative to Docker for managing software environments
  - Install: https://docs.conda.io/en/latest/miniconda.html

**The wrapper scripts will check for these dependencies and tell you what's missing before running.**

---

## Two Options Available

### Option 1: Graphical User Interface (GUI) - Easiest for Beginners

**File:** `run_kocheck_gui.py`

The GUI provides a point-and-click interface with file browsers for selecting inputs.

#### Requirements
- Python 3.6 or higher
- Tkinter (usually included with Python)

#### Usage

```bash
python3 run_kocheck_gui.py
```

A window will open with input fields and browse buttons. Fill in the required fields, adjust optional parameters if needed, and click "Run Pipeline".

**Features:**
- File browser dialogs for easy file selection
- Clear labeling of required vs. optional parameters
- Real-time output display
- Helpful tooltips with default values
- Clear success/error messages

---

### Option 2: Command-Line Script - For Those Comfortable with Terminal

**File:** `run_kocheck.sh`

A bash script wrapper that simplifies command-line usage with helpful defaults and validation.

#### Usage

```bash
bash run_kocheck.sh -r "data/*_{R1,R2}.fq.gz" \
                    -f reference.fasta \
                    -b target_gene.bed \
                    -m marker.fasta
```

#### Required Options
- `-r, --reads PATTERN` - Path pattern to paired-end FASTQ files
- `-f, --reference FILE` - Reference genome FASTA file
- `-b, --bed FILE` - Target gene BED file
- `-m, --marker FILE` - Marker sequence FASTA file

#### Optional Options
- `-c, --marker-contig NAME` - Marker contig name (default: kanR)
- `-o, --outdir DIR` - Output directory (default: results)
- `--delete-ratio NUM` - Deletion threshold (default: 0.1)
- `--intact-ratio NUM` - Intact threshold (default: 0.5)
- `--flank SIZE` - Flanking region in bp (default: 500)
- `--min-mapq QUAL` - Minimum mapping quality (default: 20)
- `-h, --help` - Show help message

#### Examples

Test with provided test data:
```bash
bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
                    -f assets/testdata/ref.fasta \
                    -b assets/testdata/target_gene.bed \
                    -m assets/testdata/kanR.fasta
```

Run with custom data:
```bash
bash run_kocheck.sh -r "my_data/*_{R1,R2}.fq.gz" \
                    -f my_data/reference.fasta \
                    -b my_data/target.bed \
                    -m my_data/resistance_marker.fasta \
                    -o my_results \
                    --delete-ratio 0.15 \
                    --intact-ratio 0.6
```

---

## Making Scripts Executable

### GUI Script
```bash
chmod +x run_kocheck_gui.py
```

Then you can run it as:
```bash
./run_kocheck_gui.py
```

### Shell Script
```bash
chmod +x run_kocheck.sh
```

Then you can run it as:
```bash
./run_kocheck.sh [OPTIONS]
```

---

## Understanding the Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| FASTQ Files Pattern | Glob pattern to your paired-end reads | `data/*_{R1,R2}.fq.gz` |
| Reference Genome | FASTA file of reference sequence | `reference.fasta` |
| Target Gene BED | Single-line BED file with target coordinates | `target_gene.bed` |
| Marker Sequence | FASTA file with resistance marker | `kanR.fasta` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Marker Contig Name | `kanR` | Name of the marker in combined reference |
| Output Directory | `results` | Where to save results |
| Deletion Ratio | `0.1` | Coverage threshold for detecting deletions |
| Intact Ratio | `0.5` | Coverage threshold for intact genes |
| Flank Region | `500` | bp on each side of target for visualization |
| Min Mapping Quality | `20` | Quality threshold for reads |

---

## Understanding Results

The pipeline generates several output files in the results directory:

- `qc/` - Quality control reports from read trimming
- `mapping/` - Aligned BAM files and indices
- `coverage/` - Coverage analysis files
- `deletion_check/` - CSV with deletion status for each sample
- `marker_check/` - CSV with marker validation results
- `plots/` - Coverage plots for visual inspection
- `summary/` - Final summary report (HTML and CSV)

The most important file is typically the summary report which shows:
- **WT**: Wild-type sample (no marker, gene intact)
- **OK**: Correct knockout (marker present, gene deleted)
- **WRONG_SITE**: Marker inserted at wrong location
- **ECTOPIC**: Marker has integrated elsewhere in genome
- **AMBIGUOUS**: Unclear classification

---

## Troubleshooting

### GUI won't open or shows "Missing Dependencies"
**If you see this error, you need to install the missing tools:**

1. **Nextflow not installed:**
   ```bash
   curl -s https://get.nextflow.io | bash
   mv nextflow ~/bin/  # Or add to PATH
   ```

2. **Docker not installed:**
   - Visit https://www.docker.com/get-started and follow platform-specific instructions
   - On Mac/Linux, verify it's running: `docker ps`

3. **Conda not installed (as alternative to Docker):**
   - Download Miniconda: https://docs.conda.io/en/latest/miniconda.html
   - Verify installation: `conda --version`

### GUI won't open but dependencies are installed
- Make sure Python 3 is installed: `python3 --version`
- Check that tkinter is available: `python3 -m tkinter`
- Try running from terminal to see error messages: `python3 run_kocheck_gui.py`

### Shell script shows "command not found"
- Verify Nextflow is in your PATH: `which nextflow`
- Verify Docker or Conda is available: `which docker` or `which conda`
- Make sure the script is executable: `chmod +x run_kocheck.sh`

### Pipeline fails to start
The wrappers will check dependencies before running, but if you see errors:
- Ensure Nextflow is installed: `nextflow -version`
- Make sure Docker is running (if using Docker): `docker ps`
- Verify Conda is available (if using Conda): `conda --version`
- Check that all input files exist and are readable

### Docker/Conda related errors during pipeline
- **Docker:** Make sure Docker Desktop is running (on Mac/Windows) or Docker daemon is active (Linux)
- **Conda:** The pipeline will automatically download and set up required environments
- Check available disk space (environments take several GB)

### Results look wrong
- Verify your BED file contains exactly one line with format: `chrom start end`
- Check that FASTQ files match the pattern you provided
- Ensure reference genome and marker sequence are in FASTA format
- Check the quality control reports in `results/qc/` for data quality issues

For more detailed help, see the main README.md or run the pipeline with `--help` if using direct Nextflow commands.
