// modules/fastp.nf
process FASTP {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
        path "${sample_id}.fastp.html", emit: html
        path "${sample_id}.fastp.json", emit: json

    publishDir "${params.outdir}/qc", mode: 'copy'    

    script:
    """
    fastp \
      -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_R1.trim.fastq.gz -O ${sample_id}_R2.trim.fastq.gz \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --trim_poly_g \
      --trim_poly_x \
      --unqualified_percent_limit 40 \
      --thread ${task.cpus} \
      --html ${sample_id}.fastp.html \
      --json ${sample_id}.fastp.json
    """
}

