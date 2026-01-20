# Build stage - compile and prepare everything
FROM eclipse-temurin:17-jdk-jammy as builder

WORKDIR /build

# Install base system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    git \
    bzip2 \
    build-essential \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (ARM-safe)
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MINIFORGE=Miniforge3-Linux-x86_64.sh ; \
    else \
        MINIFORGE=Miniforge3-Linux-aarch64.sh ; \
    fi && \
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/${MINIFORGE} -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

# Add conda/mamba to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Update conda and install mamba
RUN conda install -y mamba -n base -c conda-forge && conda clean -afy

# Create single GUI environment with all KOCheck dependencies
RUN mamba create -y -n kocheck_env \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    fastp \
    samtools \
    bwa-mem2 \
    mosdepth \
    python=3.10 \
    pandas \
    matplotlib \
    && mamba clean -afy

# Activate environment by default
ENV CONDA_DEFAULT_ENV=kocheck_env
ENV PATH="/opt/conda/envs/kocheck_env/bin:${PATH}"

# Clean up conda and pip caches
RUN mamba clean --all -y

# Install Nextflow in builder stage
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Runtime stage - minimal image
FROM eclipse-temurin:17-jre-jammy

WORKDIR /home/kocheck

# Copy Nextflow from builder
COPY --from=builder /usr/local/bin/nextflow /usr/local/bin/nextflow

# Copy other installed packages and tools from builder
COPY --from=builder /opt/conda /opt/conda

# Activate conda environment
ENV CONDA_DEFAULT_ENV=kocheck_env
ENV PATH="/opt/conda/envs/kocheck_env/bin:/opt/conda/bin:${PATH}"

# Create a non-root user
RUN useradd -m -s /bin/bash kocheck

# Copy the KOCheck repository
COPY . /home/kocheck/

# Set permissions
RUN chown -R kocheck:kocheck /home/kocheck

# Switch to non-root user
USER kocheck

# Set environment variables
ENV PATH="/home/kocheck:${PATH}"
ENV NXF_HOME="/home/kocheck/.nextflow"

# Expose port for any web services
EXPOSE 8080

# Default command
CMD ["/bin/bash"]
