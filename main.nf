nextflow.enable.dsl=2

// Import modules
include { FASTP } from './modules/fastp.nf'

// Workflow params
workflow {

    // Minimal input channel for one test sample
    reads = Channel.fromFilePairs(params.reads)
    reference = file('assets/testdata/ref.fasta')
    gene_bed = file('assets/testdata/target_gene.bed')

    // Echo initial inputs
    reads.view { sample_id, files ->
        "Sample: $sample_id\nReads: ${files.join(', ')}"
    }

    println "Reference genome: $reference"
    println "Gene BED: $gene_bed"
    
    // Run FASTP module

trimmed_reads = FASTP(reads)

}
