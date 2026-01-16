# Installation Guide for KOCheck

This document describes multiple ways to set up KOCheck depending on your preferences and environment.

## Quick Summary

| Method | Best For | Difficulty |
|--------|----------|-----------|
| **Docker** | Most users, any OS | Easy |
| **Docker Compose** | Simplified Docker setup | Easy |
| **Conda** | Linux/Mac users, no Docker | Medium |
| **Manual** | Development/customization | Hard |

---

## Method 1: Docker (Recommended for Most Users)

### Prerequisites
- Docker installed: https://www.docker.com/get-started

### Steps

1. **Clone or download KOCheck**
   ```bash
   git clone https://github.com/meadm/KOCheck
   cd KOCheck
   ```

2. **Build the Docker image** (first time only)
   ```bash
   docker build -t kocheck:latest .
   ```

3. **Run the pipeline**
   ```bash
   docker run --rm -v $(pwd):/home/kocheck/project kocheck:latest \
     bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
     -f assets/testdata/ref.fasta \
     -b assets/testdata/target_gene.bed \
     -m assets/testdata/kanR.fasta
   ```

### Using with Your Own Data

```bash
docker run --rm -v /path/to/your/data:/data kocheck:latest \
  bash run_kocheck.sh -r "/data/*_{R1,R2}.fq.gz" \
  -f /data/reference.fasta \
  -b /data/target.bed \
  -m /data/marker.fasta \
  -o /data/results
```

---

## Method 2: Docker Compose (Easiest Setup)

### Prerequisites
- Docker and Docker Compose installed

### Steps

1. **Clone KOCheck**
   ```bash
   git clone https://github.com/meadm/KOCheck
   cd KOCheck
   ```

2. **Build and enter container**
   ```bash
   docker-compose run kocheck /bin/bash
   ```

3. **Inside container, run the pipeline**
   ```bash
   bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
     -f assets/testdata/ref.fasta \
     -b assets/testdata/target_gene.bed \
     -m assets/testdata/kanR.fasta
   ```

4. **Exit container**
   ```bash
   exit
   ```

### Using the GUI with Docker Compose (Mac/Linux with X11)

```bash
# On Mac, first install XQuartz and allow connections:
# XQuartz > Preferences > Security > Allow connections from network clients

# Start the GUI
docker-compose -f docker-compose.gui.yml run kocheck-gui
```

---

## Method 3: Conda (Linux/Mac Alternative to Docker)

### Prerequisites
- Conda installed: https://docs.conda.io/en/latest/miniconda.html

### Steps

1. **Clone KOCheck**
   ```bash
   git clone https://github.com/meadm/KOCheck
   cd KOCheck
   ```

2. **Create the Conda environment**
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment**
   ```bash
   conda activate kocheck
   ```

4. **Run the pipeline**
   ```bash
   bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
     -f assets/testdata/ref.fasta \
     -b assets/testdata/target_gene.bed \
     -m assets/testdata/kanR.fasta
   ```

5. **Deactivate when done**
   ```bash
   conda deactivate
   ```

### Using the GUI with Conda

```bash
conda activate kocheck
python3 run_kocheck_gui.py
```

---

## Method 4: Manual Installation (For Development)

### Prerequisites
- Nextflow >= 22.10.0: https://www.nextflow.io/
- Java 11+
- Docker or Conda (for tools)
- Bioinformatics tools: samtools, bwa, mosdepth, fastp

### Steps

1. **Install Nextflow**
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```

2. **Install dependencies** (choose one)
   
   **Option A: Using Conda**
   ```bash
   conda create -n kocheck-dev \
     -c bioconda -c conda-forge \
     fastp bwa samtools mosdepth seqkit bedtools
   conda activate kocheck-dev
   ```
   
   **Option B: Using your package manager (Ubuntu/Debian)**
   ```bash
   sudo apt-get install samtools bwa mosdepth fastp seqkit bedtools
   ```

3. **Clone KOCheck**
   ```bash
   git clone https://github.com/meadm/KOCheck
   cd KOCheck
   ```

4. **Run the pipeline**
   ```bash
   bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
     -f assets/testdata/ref.fasta \
     -b assets/testdata/target_gene.bed \
     -m assets/testdata/kanR.fasta
   ```

---

## Testing Your Installation

### Quick Test (< 5 minutes)

```bash
# Docker
docker run --rm -v $(pwd):/home/kocheck/project kocheck:latest \
  bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta \
  -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta

# Conda
conda activate kocheck
bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta \
  -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta

# Direct (if manually installed)
bash run_kocheck.sh -r "assets/testdata/*_{R1,R2}.fq.gz" \
  -f assets/testdata/ref.fasta \
  -b assets/testdata/target_gene.bed \
  -m assets/testdata/kanR.fasta
```

Expected output: Results in `results/` directory with CSV files and HTML report.

---

## Troubleshooting

### Docker Issues

**"Docker daemon is not running"**
- On Mac: Open Docker Desktop app
- On Windows: Open Docker Desktop app
- On Linux: `sudo systemctl start docker`

**"Cannot find image kocheck:latest"**
- Build the image: `docker build -t kocheck:latest .`

**Permissions error on Linux**
- Add your user to docker group: `sudo usermod -aG docker $USER`
- Then log out and back in

### Conda Issues

**"conda: command not found"**
- Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
- Make sure conda is in your PATH

**Environment fails to create**
- Update conda: `conda update conda`
- Try: `conda env create -f environment.yml --force`

**Tools not found after activating environment**
- Make sure you're in the right environment: `conda info --envs`
- Reactivate: `conda deactivate && conda activate kocheck`

### General Issues

**"Nextflow not found"**
- Install Nextflow: `curl -s https://get.nextflow.io | bash`
- Add to PATH: `export PATH=$PATH:$(pwd)`

**Pipeline fails during execution**
- Check dependencies are installed: `nextflow -version`
- Ensure input files exist: `ls -la assets/testdata/`
- Check disk space: `df -h`

---

## Getting Help

For more information, see:
- Main README: `README.md`
- Wrapper usage: `WRAPPER_USAGE.md`
- GitHub issues: https://github.com/meadm/KOCheck/issues

---

## Advanced: Building Custom Docker Image

To customize the Docker image (e.g., add more tools), edit `Dockerfile` and rebuild:

```bash
docker build -t kocheck:custom -f Dockerfile .
```

To push to Docker Hub:

```bash
docker tag kocheck:latest your-username/kocheck:latest
docker push your-username/kocheck:latest
```

Others can then pull it with:

```bash
docker pull your-username/kocheck:latest
```
