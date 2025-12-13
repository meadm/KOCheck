process DELETION_CHECK {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(summary_txt), path(gene_bed)

    output:
    tuple val(sample_id), path("${sample_id}.deletion_check.txt")

    script:
    """
    # Extract gene coordinates (BED file must have 3 columns: chrom start end)
    chrom=\$(awk '{print \$1}' ${gene_bed})
    start=\$(awk '{print \$2}' ${gene_bed})
    end=\$(awk '{print \$3}' ${gene_bed})

    # Calculate gene length
    gene_len=\$(( end - start ))

    # Extract gene coverage summary from mosdepth summary file
    # Filter for exact chromosome match (not the _region line)
    chrom_mean=\$(grep -P "^${chrom}\t" ${summary_txt} | awk '{print \$4}')

    # Extract gene coverage using mosdepth per-region stats (quick & dirty):
    # mosdepth always creates <prefix>.regions.bed.gz for per-region coverage.
    # We'll approximate gene region coverage using samtools depth.
    mean_gene_cov=\$(samtools depth -r ${chrom}:\$start-\$end ${bam} | \
                     awk '{sum+=\$3} END { if (NR>0) print sum/NR; else print 0 }')

    # Decide deletion using a simple threshold:
    #   If gene coverage < 10% of chromosome mean → deleted
    status="intact"
    threshold=\$(echo "\$chrom_mean * 0.1" | bc)

    below=\$(echo "\$mean_gene_cov < \$threshold" | bc)
    if [ "\$below" -eq 1 ]; then
        status="deleted"
    fi

    # Write output
    echo -e "sample\\tstatus\\tchrom_mean_cov\\tgene_mean_cov" > ${sample_id}.deletion_check.txt
    echo -e "${sample_id}\\t\$status\\t\$chrom_mean\\t\$mean_gene_cov" >> ${sample_id}.deletion_check.txt
    """
}

