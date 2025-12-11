nextflow.enable.dsl=2

// Import modules
include { FASTP } from './modules/fastp.nf'
include { BWA_MEM2 } from './modules/bwa_mem2.nf'

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
    FASTP(reads)
    
    // Run BWA-MEM2 alignment (aligns, sorts, and indexes in one step)
    // BWA_MEM2 expects: (sample_id, r1, r2, reference) in a single tuple
    // Access the trimmed_reads output specifically, then map to add reference
    trimmed_reads_with_ref = FASTP.out.trimmed_reads.map { sample_id, r1, r2 -> 
        [sample_id, r1, r2, reference] 
    }
    aligned_bam = BWA_MEM2(trimmed_reads_with_ref)

}
