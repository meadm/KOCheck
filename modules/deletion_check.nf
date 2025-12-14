process DELETION_CHECK {

    tag "${sample_id}"

    publishDir "${params.outdir}/deletion_check", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(summary_txt), path(gene_bed), path(bai)

    output:
    path "${sample_id}.deletion_check.csv"

    script:
    """
    set -euo pipefail

    ############################
    # Validate BED file
    ############################

    # Must contain exactly one line
    if [ \$(wc -l < ${gene_bed}) -ne 1 ]; then
        echo "ERROR: gene BED file must contain exactly one line" >&2
        exit 1
    fi

    read chrom start end < ${gene_bed}

    # Start must be less than end
    if [ "\$start" -ge "\$end" ]; then
        echo "ERROR: BED coordinates invalid (start >= end): \$start \$end" >&2
        exit 1
    fi

    ############################
    # Get mean coverage from the mosdepth summary of the chromosome defined in the chrom variable
    ############################

    chrom_mean=\$(awk -v chr="\$chrom" '\$1 == chr {print \$4}' ${summary_txt})

    if [ -z "\$chrom_mean" ]; then
        echo "ERROR: Chromosome \$chrom not found in mosdepth summary" >&2
        exit 1
    fi

    ############################
    # Compute mean gene coverage
    ############################

    mean_gene_cov=\$(samtools depth -a -r "\$chrom:\$start-\$end" ${bam} | \
        awk '{sum+=\$3} END { if (NR>0) print sum/NR; else print 0 }')
    
    ############################
    # Classify deletion status
    ############################

    ratio=\$(awk -v g="\$mean_gene_cov" -v c="\$chrom_mean" \
        'BEGIN { if (c>0) print g/c; else print 0 }')

    status="ambiguous"

    is_deleted=\$(awk -v r="\$ratio" -v t="${params.delete_ratio}" \
        'BEGIN { print (r < t) }')

    is_intact=\$(awk -v r="\$ratio" -v t="${params.intact_ratio}" \
        'BEGIN { print (r > t) }')

    if [ "\$is_deleted" -eq 1 ]; then
        status="deleted"
    elif [ "\$is_intact" -eq 1 ]; then
        status="intact"
    fi

    ############################
    # Write CSV output
    ############################

    echo "sample_id,deletion_status" > ${sample_id}.deletion_check.csv
    echo "${sample_id},\${status}" >> ${sample_id}.deletion_check.csv
    """
}

