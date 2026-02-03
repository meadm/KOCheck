# KOCheck Pipeline Workflow

## Data Flow Diagram

```mermaid
graph TD
    A["📁 FASTQ Files<br/>(Paired-end reads)"] --> B["✂️ FASTP<br/>(QC & Trimming)"]
    B --> C["🗺️ BWA-MEM2<br/>(Alignment)"]
    C --> D["🔍 MOSDEPTH<br/>(Coverage Analysis)"]
    D --> E["🔍 DELETION_CHECK<br/>(Gene Coverage)"]
    E --> F["🔍 MARKER_CHECK<br/>(Marker & Junctions Analysis)"]
    F --> G["📈 COVERAGE_PLOT<br/>(Visualization)"]
    G --> H["📋 AGGREGATE_RESULTS<br/>(Summary Report)"]
    H --> H2["🎯 Final Classification<br/>(PASS / REVIEW / FAIL)"]
    
    B -.-> B1["📊 QC Reports<br/>(*.fastp.html/.json)"]
    C -.-> C1["📦 BAM Files<br/>(*.bam, *.bam.bai)"]
    D -.-> D1["📚 Coverage Files<br/>(*.mosdepth.summary.txt)"]
    E -.-> E1["📚 deletion_status.csv<br/>(Coverage Ratio)"]
    F -.-> F1["📚 marker_status.csv<br/>(Marker Present?<br/>Junction Support)"]
    G -.-> G1["📈 Coverage Plots<br/>(*.coverage.png)"]
    H -.-> H1["📚 summary.csv<br/>🖼️ summary.html"]
    
    style A fill:#1e88e5,color:#fff
    style H2 fill:#00c853,color:#fff
    style B1 fill:#ff6f00,color:#fff
    style C1 fill:#ff6f00,color:#fff
    style D1 fill:#ff6f00,color:#fff
    style E1 fill:#ff6f00,color:#fff
    style F1 fill:#ff6f00,color:#fff
    style G1 fill:#ff6f00,color:#fff
    style H1 fill:#ff6f00,color:#fff
```

## Processing Steps

### 1. **FASTP - Quality Control** 
- Removes adapters from reads
- Trims low-quality bases
- Generates QC reports
- Output: Cleaned FASTQ files

### 2. **BWA-MEM2 - Alignment**
- Maps reads to reference genome
- Creates BAM files with alignment coordinates
- Outputs sorted, indexed BAM files
- Output: `*.bam`, `*.bam.bai`

### 3. **MOSDEPTH - Coverage Analysis**
- Calculates per-base coverage across genome
- Generates coverage summary statistics
- Creates BED files for coverage regions
- Output: Coverage summary files

### 4. **DELETION_CHECK - Gene Status**
- Compares coverage in target gene vs. flanking regions
- Calculates coverage ratio
- **Classification:**
  - **"deleted"**: coverage_ratio < delete_ratio (default 0.1)
  - **"intact"**: coverage_ratio > intact_ratio (default 0.5)
  - **"ambiguous"**: between thresholds

### 5. **MARKER_CHECK - Marker Validation**
- Checks if marker sequence is present in reads
- Analyzes junction support (proper marker integration)
- Detects ectopic insertions (marker elsewhere in genome)
- **Classification:**
  - **"present"**: marker detected with junction support
  - **"absent"**: no marker detected
  - **"ectopic"**: marker present but wrong location

### 6. **COVERAGE_PLOT - Visualization**
- Generates PNG plots of coverage across target region
- Shows target gene, flanking regions, coverage profile
- Output: `*.png` coverage plots

### 7. **AGGREGATE_RESULTS - Final Report**
- Combines results from all samples
- Creates summary CSV and HTML reports
- Final classification for each sample
- Output: `summary.csv`, `summary.html`

## Output Files

```
results/
├── qc/                          # Quality control reports
│   ├── *.fastp.html
│   └── *.fastp.json
├── mapping/                     # Aligned reads
│   ├── *.bam
│   └── *.bam.bai
├── coverage/                    # Coverage analysis
│   ├── *.mosdepth.summary.txt
│   ├── *.per-base.bed.gz
│   └── *.regions.bed.gz
├── deletion_check/              # Deletion classification
│   └── deletion_status.csv
├── marker_check/                # Marker validation
│   └── marker_status.csv
├── plots/                       # Coverage visualizations
│   └── *.coverage.png
└── summary/                     # Final reports
    ├── summary.csv              # ⭐ Main results table
    └── summary.html             # ⭐ Interactive HTML report
```

## Sample Classification Logic

See `CLASSIFICATION_LOGIC.md` for detailed classification rules.
