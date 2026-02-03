# Build stage - compile and prepare everything
FROM eclipse-temurin:17-jdk-jammy AS builder

WORKDIR /build

# Install only build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    bzip2 \
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

# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Update conda and install mamba
RUN conda install -y mamba -n base -c conda-forge && conda clean -afy

# Create kocheck environment with pinned versions
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
    && mamba clean --all -y && \
    rm -rf /opt/conda/pkgs /opt/conda/share/doc /opt/conda/share/man

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Runtime stage - slim base
FROM eclipse-temurin:17-jre-jammy

WORKDIR /home/kocheck

# Install only runtime dependencies (not build tools)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy Nextflow from builder
COPY --from=builder /usr/local/bin/nextflow /usr/local/bin/nextflow

# Copy ONLY the kocheck_env (not all of /opt/conda)
COPY --from=builder /opt/conda/envs/kocheck_env /opt/conda/envs/kocheck_env

# Minimal conda setup - symlink binaries
RUN mkdir -p /opt/conda/bin && \
    ln -s /opt/conda/envs/kocheck_env/bin/* /opt/conda/bin/ 2>/dev/null || true

# Activate environment
ENV CONDA_DEFAULT_ENV=kocheck_env
ENV PATH="/opt/conda/envs/kocheck_env/bin:/opt/conda/bin:${PATH}"

# Create non-root user
RUN useradd -m -s /bin/bash kocheck

# Copy KOCheck repository
COPY --chown=kocheck:kocheck . /home/kocheck/

# Create Nextflow directories and set proper permissions
# Need to do this as root before switching to kocheck user
RUN mkdir -p /home/kocheck/.nextflow /home/kocheck/work && \
    chown -R kocheck:kocheck /home/kocheck && \
    chmod -R u+w /home/kocheck

# Switch to non-root user
USER kocheck

# Environment variables
ENV PATH="/home/kocheck:${PATH}"
ENV NXF_HOME="/home/kocheck/.nextflow"

EXPOSE 8080

CMD ["/bin/bash"]
