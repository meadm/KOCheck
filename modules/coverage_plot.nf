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

    # Read gene coordinates
    read chrom start end < ${gene_bed}

    # Make samtools index for total genome size calculation
    samtools faidx ${reference}

    # Pull chromosome length
    chrom_len=\$(awk -v c="\$chrom" '\$1==c {print \$2}' ${reference}.fai)

    # Define flanking region (bp)
    FLANK=${params.flank}
    raw_start=\$(( start - FLANK ))
    raw_end=\$(( end + FLANK ))

    # Make sure the flanks aren't outside the chromosome
    plot_start=\$(( raw_start < 0 ? 0 : raw_start ))
    plot_end=\$(( raw_end > chrom_len ? chrom_len : raw_end ))

    # Extract per-base coverage for the plotting window
    samtools depth -a -r "\$chrom:\$plot_start-\$plot_end" ${bam} >  ${sample_id}_coverage.tsv  

    # Generate coverage plot
    python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt

sample_id = "${sample_id}"
chrom = "\${chrom}"
start = int(\${start})
end = int(\${end})
coverage_file = "${sample_id}_coverage.tsv"

# Load coverage data
cols = ["chrom", "start", "cov"]
df = pd.read_csv(coverage_file, sep="\t", header=None, names=cols)

plt.figure(figsize=(10, 3))
plt.plot(df["start"], df["cov"], linewidth=1)

# Highlight gene region
plt.axvspan(start, end, alpha=0.1)

plt.xlabel("Genomic position")
plt.ylabel("Coverage")
plt.title(f"{sample_id} – {chrom}:{start}-{end}")

plt.tight_layout()
plt.savefig(f"{sample_id}.gene_coverage.png", dpi=200)
EOF
    """
}

