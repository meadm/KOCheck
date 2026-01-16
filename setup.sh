#!/bin/bash
#
# KOCheck Setup Script
# Automatically sets up KOCheck with your choice of Docker or Conda
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_header() {
    echo -e "\n${BLUE}═══════════════════════════════════════${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}═══════════════════════════════════════${NC}\n"
}

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

print_header "KOCheck Setup Wizard"

# Check what's available
HAS_DOCKER=0
HAS_CONDA=0
HAS_NEXTFLOW=0

if command -v docker &> /dev/null; then
    HAS_DOCKER=1
    DOCKER_VERSION=$(docker --version)
    echo -e "${GREEN}✓${NC} Docker found: $DOCKER_VERSION"
else
    echo -e "${RED}✗${NC} Docker not found"
fi

if command -v conda &> /dev/null; then
    HAS_CONDA=1
    CONDA_VERSION=$(conda --version)
    echo -e "${GREEN}✓${NC} Conda found: $CONDA_VERSION"
else
    echo -e "${RED}✗${NC} Conda not found"
fi

if command -v nextflow &> /dev/null; then
    HAS_NEXTFLOW=1
    NXF_VERSION=$(nextflow -version | head -1)
    echo -e "${GREEN}✓${NC} Nextflow found: $NXF_VERSION"
else
    echo -e "${RED}✗${NC} Nextflow not found"
fi

echo ""

# Determine setup method
if [ $HAS_DOCKER -eq 0 ] && [ $HAS_CONDA -eq 0 ]; then
    print_error "Neither Docker nor Conda found. Please install one:"
    echo "  Docker: https://www.docker.com/get-started"
    echo "  Conda:  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

print_info "Choose setup method:"
echo ""

if [ $HAS_DOCKER -eq 1 ]; then
    echo "  1) Docker (recommended - works on all platforms)"
    OPTION_DOCKER=1
else
    OPTION_DOCKER=0
fi

if [ $HAS_CONDA -eq 1 ]; then
    echo "  2) Conda (faster, runs on your system)"
    OPTION_CONDA=2
else
    OPTION_CONDA=0
fi

echo ""
read -p "Enter your choice (1 or 2): " CHOICE

case $CHOICE in
    1)
        if [ $HAS_DOCKER -eq 0 ]; then
            print_error "Docker not available"
            exit 1
        fi
        
        print_header "Setting up with Docker"
        
        print_info "Building Docker image..."
        docker build -t kocheck:latest .
        
        print_success "Docker image built successfully!"
        echo ""
        print_info "You can now run the pipeline with:"
        echo ""
        echo "  docker run --rm -v \$(pwd):/home/kocheck/project kocheck:latest \\"
        echo "    bash run_kocheck.sh -r \"assets/testdata/*_{R1,R2}.fq.gz\" \\"
        echo "    -f assets/testdata/ref.fasta \\"
        echo "    -b assets/testdata/target_gene.bed \\"
        echo "    -m assets/testdata/kanR.fasta"
        echo ""
        print_info "Or use docker-compose:"
        echo "  docker-compose run kocheck bash run_kocheck.sh [OPTIONS]"
        echo ""
        ;;
    
    2)
        if [ $HAS_CONDA -eq 0 ]; then
            print_error "Conda not available"
            exit 1
        fi
        
        print_header "Setting up with Conda"
        
        print_info "Creating Conda environment..."
        conda env create -f environment.yml
        
        print_success "Conda environment created!"
        echo ""
        print_info "To use it, run:"
        echo "  conda activate kocheck"
        echo ""
        print_info "Then you can run the pipeline:"
        echo "  bash run_kocheck.sh -r \"assets/testdata/*_{R1,R2}.fq.gz\" \\"
        echo "    -f assets/testdata/ref.fasta \\"
        echo "    -b assets/testdata/target_gene.bed \\"
        echo "    -m assets/testdata/kanR.fasta"
        echo ""
        print_info "When done, deactivate with:"
        echo "  conda deactivate"
        echo ""
        ;;
    
    *)
        print_error "Invalid choice"
        exit 1
        ;;
esac

print_header "Setup Complete!"
print_info "For more information, see INSTALLATION.md or WRAPPER_USAGE.md"
echo ""
