process BWA_MEM2 {

    tag "$sample_id"

    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2), path(reference)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    # Index reference if necessary
    if [ ! -f "${reference}.bwt.2bit.64" ]; then
        bwa-mem2 index ${reference}
    fi

    # Align
    bwa-mem2 mem -t ${task.cpus} ${reference} ${r1} ${r2} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """
}

