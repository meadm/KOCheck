nextflow.enable.dsl=2

// Import modules
include { FASTP } from './modules/fastp.nf'
include { BWA_MEM2 } from './modules/bwa_mem2.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { DELETION_CHECK } from './modules/deletion_check.nf'
include { COVERAGE_PLOT } from './modules/coverage_plot.nf'

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
    
    // Run MOSDEPTH module
    coverage = MOSDEPTH(aligned_bam)

    // Run DELETION_CHECK module
    // DELETION_CHECK needs: (sample_id, bam, summary_txt, gene_bed)
    // Combine BAM from BWA_MEM2, summary from MOSDEPTH, and gene_bed
    deletion_check_input = aligned_bam.join(coverage).map { sample_id, bam, bai, per_base, regions, summary_txt ->
        [sample_id, bam, summary_txt, gene_bed, bai]
    }
    deletion_results = DELETION_CHECK(deletion_check_input)

    //Run COVERAGE_PLOT module
    //COVERAGE_PLOT needs: (sample_id, bam, bai, gene_bed)
    //Combine BAM from BWA_MEM2 and the gene bed
    coverage_plot_input = aligned_bam.map { sample_id, bam, bai ->
        tuple(sample_id, bam, bai, file(params.gene_bed), file(params.reference))
    }

    COVERAGE_PLOT(coverage_plot_input)

}
