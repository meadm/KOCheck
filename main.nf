nextflow.enable.dsl=2

workflow {

    // Minimal input channel for one test sample
    reads = Channel.fromFilePairs(params.reads)
    reference = file('assets/testdata/ref.fasta')
    gene_bed = file('assets/testdata/target_gene.bed')

    // For now, just echo inputs and outputs to test workflow
    reads.view { sample_id, files ->
        "Sample: $sample_id\nReads: ${files.join(', ')}"
    }

    println "Reference genome: $reference"
    println "Gene BED: $gene_bed"

    // Placeholder step: just copy input reads to output folder
    DUMMY_COPY(reads)
    
    DUMMY_COPY.out[0].view { f -> "Generated dummy output R1: $f" }
    DUMMY_COPY.out[1].view { f -> "Generated dummy output R2: $f" }
}

// Placeholder step: just copy input reads to output folder
process DUMMY_COPY {
    input:
        tuple val(sample_id), path(reads)
    output:
        path("${sample_id}_processed_R1.fq.gz")
        path("${sample_id}_processed_R2.fq.gz")

    script:
        """
        cp ${reads[0]} ${sample_id}_processed_R1.fq.gz
        cp ${reads[1]} ${sample_id}_processed_R2.fq.gz
        """
}

