process MOSDEPTH {
    tag "$sample_id"

    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.per-base.bed.gz"),
          path("${sample_id}.regions.bed.gz"),
          path("${sample_id}.mosdepth.summary.txt")

    script:
    """
    mosdepth \
        --threads ${task.cpus} \
        --by 1000 \
        ${sample_id} \
        ${bam}
    """
}


