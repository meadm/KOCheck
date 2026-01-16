# KOCheck - Quick Start Guide

## What is KOCheck?

KOCheck is a reproducible Nextflow pipeline for validating genetic knockouts. It analyzes Illumina sequencing reads to confirm gene deletions, detect correct resistance marker insertion, and identify ectopic integrations.

## Quick Installation (Choose One)

### 🐳 Option 1: Docker (Recommended - Easiest)
```bash
# Clone and build
git clone https://github.com/meadm/KOCheck
cd KOCheck
docker build -t kocheck:latest .

# Run with test data
docker run --rm -v $(pwd):/home/kocheck/project kocheck:latest \
  bash run_kocheck.sh \
  -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta \
  -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta
```

### 🔧 Option 2: Conda (Alternative)
```bash
# Clone and setup
git clone https://github.com/meadm/KOCheck
cd KOCheck
conda env create -f environment.yml
conda activate kocheck

# Run with test data
bash run_kocheck.sh \
  -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta \
  -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta
```

### 🤖 Option 3: Automated Setup
```bash
bash setup.sh
```

## Running the Pipeline

### Method 1: Graphical Interface (GUI)
```bash
python3 run_kocheck_gui.py
```

### Method 2: Command Line
```bash
bash run_kocheck.sh \
  --reads "path/to/*_{R1,R2}.fq.gz" \
  --reference path/to/reference.fasta \
  --bed path/to/target_gene.bed \
  --marker path/to/marker.fasta
```

See `WRAPPER_USAGE.md` for full documentation.

## Documentation

- **INSTALLATION.md** - Detailed setup instructions for all methods
- **WRAPPER_USAGE.md** - GUI and shell script usage guide
- **README.md** - Full pipeline documentation
- **WRAPPER_USAGE.md** - Complete parameter reference

## Features

✅ Quality control and adapter trimming
✅ Gene deletion detection
✅ Resistance marker validation
✅ Ectopic integration detection
✅ Coverage visualization
✅ Comprehensive reporting (CSV + HTML)

## Results

The pipeline generates:
- **CSV files** with deletion and marker status for each sample
- **Coverage plots** showing target gene coverage
- **HTML report** with interactive results and embedded plots
- **Quality control reports** from sequence trimming

## Quick Test

All methods above include test data (`assets/testdata/`) so you can verify your setup works:

```bash
# Docker
docker run --rm -v $(pwd):/home/kocheck/project kocheck:latest \
  bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta

# Conda
conda activate kocheck
bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta
```

## Requirements

**Choose one:**
- Docker (v20.10+) - https://www.docker.com/get-started
- Conda - https://docs.conda.io/en/latest/miniconda.html

**Nextflow is installed automatically with Docker/Conda**

## Getting Help

- See detailed installation instructions: `INSTALLATION.md`
- See wrapper script documentation: `WRAPPER_USAGE.md`
- Report issues: https://github.com/meadm/KOCheck/issues

## License

See LICENSE file

## Citation

If you use KOCheck, please cite:
[Citation information to be added]
