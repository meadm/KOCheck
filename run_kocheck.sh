#!/bin/bash
#
# KOCheck Pipeline Wrapper Script
# A simple command-line interface for running the KOCheck Nextflow pipeline
#

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
MARKER_CONTIG="kanR"
OUTPUT_DIR="results"
DELETE_RATIO="0.1"
INTACT_RATIO="0.5"
FLANK="500"
MIN_MAPQ="20"

# Function to print usage
usage() {
    cat << EOF
${BLUE}KOCheck Pipeline - Knockout Validation${NC}

Usage: $0 [OPTIONS]

${YELLOW}Required Options:${NC}
  -r, --reads PATTERN         Path pattern to paired-end FASTQ files
                              Example: 'data/*_{R1,R2}.fq.gz'
  
  -f, --reference FILE        Path to reference genome FASTA file
  
  -b, --bed FILE              Path to target gene BED file
                              (must contain exactly one line: chrom start end)
  
  -m, --marker FILE           Path to marker sequence FASTA file

${YELLOW}Optional Options:${NC}
  -c, --marker-contig NAME    Name of marker contig (default: kanR)
  
  -o, --outdir DIR            Output directory (default: results)
  
  --delete-ratio NUM          Deletion ratio threshold (default: 0.1)
  
  --intact-ratio NUM          Intact ratio threshold (default: 0.5)
  
  --flank SIZE                Flanking region in bp (default: 500)
  
  --min-mapq QUAL             Minimum mapping quality (default: 20)
  
  -h, --help                  Show this help message

${YELLOW}Examples:${NC}
  # Run with test data
  $0 -r "assets/testdata/*_{R1,R2}.fq.gz" \\
     -f assets/testdata/ref.fasta \\
     -b assets/testdata/target_gene.bed \\
     -m assets/testdata/kanR.fasta

  # Run with custom data and output directory
  $0 -r "data/*_{R1,R2}.fq.gz" \\
     -f reference.fasta \\
     -b target_gene.bed \\
     -m marker.fasta \\
     -o my_results

EOF
    exit 1
}

# Function to print messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Function to validate file exists
validate_file() {
    local file=$1
    local name=$2
    
    if [ ! -f "$file" ]; then
        print_error "$name not found: $file"
        exit 1
    fi
}

# Function to check if a command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        return 1
    fi
    return 0
}

# Parse command-line arguments
if [ $# -eq 0 ]; then
    usage
fi

READS=""
REFERENCE=""
GENE_BED=""
MARKER_FASTA=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--reads)
            READS="$2"
            shift 2
            ;;
        -f|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -b|--bed)
            GENE_BED="$2"
            shift 2
            ;;
        -m|--marker)
            MARKER_FASTA="$2"
            shift 2
            ;;
        -c|--marker-contig)
            MARKER_CONTIG="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --delete-ratio)
            DELETE_RATIO="$2"
            shift 2
            ;;
        --intact-ratio)
            INTACT_RATIO="$2"
            shift 2
            ;;
        --flank)
            FLANK="$2"
            shift 2
            ;;
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required parameters
if [ -z "$READS" ]; then
    print_error "FASTQ files pattern is required (use -r or --reads)"
    exit 1
fi

if [ -z "$REFERENCE" ]; then
    print_error "Reference genome file is required (use -f or --reference)"
    exit 1
fi

if [ -z "$GENE_BED" ]; then
    print_error "Target gene BED file is required (use -b or --bed)"
    exit 1
fi

if [ -z "$MARKER_FASTA" ]; then
    print_error "Marker sequence file is required (use -m or --marker)"
    exit 1
fi

# Validate files exist
validate_file "$REFERENCE" "Reference genome"
validate_file "$GENE_BED" "Gene BED file"
validate_file "$MARKER_FASTA" "Marker FASTA file"

# Print configuration
echo ""
print_info "KOCheck Pipeline Configuration"
echo ""
echo "  FASTQ Pattern:        $READS"
echo "  Reference:            $REFERENCE"
echo "  Target Gene BED:      $GENE_BED"
echo "  Marker FASTA:         $MARKER_FASTA"
echo "  Marker Contig:        $MARKER_CONTIG"
echo "  Output Directory:     $OUTPUT_DIR"
echo "  Delete Ratio:         $DELETE_RATIO"
echo "  Intact Ratio:         $INTACT_RATIO"
echo "  Flank Region:         $FLANK bp"
echo "  Min Mapping Quality:  $MIN_MAPQ"
echo ""

# Ask for confirmation
read -p "$(echo -e ${YELLOW}Continue with these settings? [y/N]${NC} )" -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_warning "Pipeline cancelled by user"
    exit 0
fi

echo ""
print_info "Starting KOCheck pipeline..."
echo ""

# Check for required dependencies
print_info "Checking for required dependencies..."
echo ""

MISSING_DEPS=0

if ! check_command "nextflow"; then
    print_error "Nextflow is not installed"
    echo "  Please install Nextflow from: https://www.nextflow.io/"
    echo "  Quick install: curl -s https://get.nextflow.io | bash"
    MISSING_DEPS=1
else
    echo "  ✓ Nextflow found: $(nextflow -v)"
fi

if ! check_command "docker"; then
    print_warning "Docker is not available"
    echo "  The pipeline uses Docker containers for tools."
    echo "  Install Docker from: https://www.docker.com/get-started"
    
    # Check if Conda is available as fallback
    if ! check_command "conda"; then
        print_error "Neither Docker nor Conda is installed!"
        echo "  Please install one of the following:"
        echo "  1. Docker: https://www.docker.com/get-started"
        echo "  2. Conda: https://docs.conda.io/en/latest/miniconda.html"
        MISSING_DEPS=1
    else
        echo "  ✓ Conda found: $(conda --version)"
        echo "  (Pipeline will use Conda environments instead of Docker)"
    fi
else
    echo "  ✓ Docker is running"
fi

# If Docker wasn't found but Conda was, that's okay - print info
if check_command "docker" > /dev/null 2>&1; then
    echo "  ✓ Docker is available"
fi

echo ""

if [ $MISSING_DEPS -eq 1 ]; then
    print_error "Missing required dependencies. Please install them and try again."
    echo ""
    exit 1
fi

print_success "All dependencies are available!"
echo ""

# Run the Nextflow pipeline
nextflow run main.nf \
    --reads "$READS" \
    --reference "$REFERENCE" \
    --gene_bed "$GENE_BED" \
    --marker_fasta "$MARKER_FASTA" \
    --marker_contig "$MARKER_CONTIG" \
    --outdir "$OUTPUT_DIR" \
    --delete_ratio "$DELETE_RATIO" \
    --intact_ratio "$INTACT_RATIO" \
    --flank "$FLANK" \
    --min_mapq "$MIN_MAPQ"

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    print_success "Pipeline completed successfully!"
    echo ""
    print_info "Results are available in: $OUTPUT_DIR"
    echo ""
else
    echo ""
    print_error "Pipeline failed. Check the output above for details."
    exit 1
fi
