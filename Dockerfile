FROM ubuntu:22.04

# Set working directory
WORKDIR /home/kocheck

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    openjdk-11-jre-headless \
    python3 \
    python3-pip \
    samtools \
    bwa \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Install Python wrapper dependencies
RUN pip3 install --no-cache-dir \
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
