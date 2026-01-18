# Build stage - compile and prepare everything
FROM eclipse-temurin:17-jdk-jammy as builder

WORKDIR /build

# Install only build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow in builder stage
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Runtime stage - minimal image
FROM eclipse-temurin:17-jre-jammy

WORKDIR /home/kocheck

# Install only runtime dependencies (jre instead of jdk, --no-install-recommends to skip extras)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    samtools \
    bwa \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/* /var/tmp/*

# Copy Nextflow from builder
COPY --from=builder /usr/local/bin/nextflow /usr/local/bin/nextflow

# Install Python wrapper dependencies with minimal cache
RUN pip3 install --no-cache-dir --no-deps \
    pillow \
    matplotlib

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
