process COVERAGE_PLOT {
    tag "$sample_id"

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(gene_bed), path(reference)

    output:
    path "${sample_id}.gene_coverage.png"

    script:
    """
    set -euo pipefail

    # Read gene coordinates (BED format: 0-based, half-open)
    read chrom start end < ${gene_bed}

    # Make samtools index for total genome size calculation
    samtools faidx ${reference}

    # Pull chromosome length
    chrom_len=\$(awk -v c="\$chrom" '\$1==c {print \$2}' ${reference}.fai)

    # Convert BED coordinates to 1-based for samtools
    # BED start is 0-based, so add 1 for 1-based
    # BED end is exclusive, so it's already the correct 1-based end
    samtools_start=\$(( start + 1 ))
    samtools_end=\${end}

    # Define flanking region (bp) - work in 1-based coordinates
    FLANK=${params.flank}
    raw_start=\$(( samtools_start - FLANK ))
    raw_end=\$(( samtools_end + FLANK ))

    # Make sure the flanks aren't outside the chromosome
    plot_start=\$(( raw_start < 1 ? 1 : raw_start ))
    plot_end=\$(( raw_end > chrom_len ? chrom_len : raw_end ))

    # Extract per-base coverage for the plotting window
    samtools depth -a -r "\$chrom:\${plot_start}-\${plot_end}" ${bam} >  ${sample_id}_coverage.tsv  

    # Generate coverage plot
    # Note: samtools depth outputs 1-based coordinates
    python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt

sample_id = "${sample_id}"
chrom = "\${chrom}"
# Use 1-based coordinates for plotting (matching samtools depth output)
start = int(\${samtools_start})
end = int(\${samtools_end})
coverage_file = "${sample_id}_coverage.tsv"

# Load coverage data
cols = ["chrom", "start", "cov"]
df = pd.read_csv(coverage_file, sep="\t", header=None, names=cols)

plt.figure(figsize=(10, 3))
plt.plot(df["start"], df["cov"], linewidth=1)

# Highlight gene region (using 1-based coordinates)
plt.axvspan(start, end, alpha=0.1)

plt.xlabel("Genomic position")
plt.ylabel("Coverage")
plt.title(f"{sample_id} – {chrom}:{start}-{end}")

plt.tight_layout()
plt.savefig(f"{sample_id}.gene_coverage.png", dpi=200)
EOF
    """
}

