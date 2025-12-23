nextflow.enable.dsl=2

// Import modules
include { FASTP } from './modules/fastp.nf'
include { APPEND_MARKER } from './modules/append_marker.nf'
include { BWA_MEM2 } from './modules/bwa_mem2.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { DELETION_CHECK } from './modules/deletion_check.nf'
include { COVERAGE_PLOT } from './modules/coverage_plot.nf'
include { MARKER_CHECK } from './modules/marker_check.nf'

// Workflow 
workflow {
    
    reads = Channel.fromFilePairs(params.reads)
    reference = file(params.reference)
    gene_bed = file(params.gene_bed)
    marker_fasta = file(params.marker_fasta)
    
    // Run FASTP module
    FASTP(reads)
    
    //Combine Reference with marker sequence for later marker detection
    APPEND_MARKER(reference, marker_fasta)

    // Run BWA-MEM2 alignment (aligns, sorts, and indexes in one step)
    // BWA_MEM2 expects: (sample_id, r1, r2, reference) in a single tuple
    // Access the trimmed_reads output specifically, then map to add reference
    trimmed_reads_with_ref = FASTP.out.trimmed_reads
        .combine(APPEND_MARKER.out.combined_ref)
    
    aligned_bam = BWA_MEM2(trimmed_reads_with_ref)
    
    // Run MOSDEPTH module
    coverage = MOSDEPTH(aligned_bam)

    // Run DELETION_CHECK module
    // DELETION_CHECK needs: (sample_id, bam, summary_txt, gene_bed)
    // Combine BAM from BWA_MEM2, summary from MOSDEPTH, and gene_bed
    deletion_check_input = aligned_bam.join(coverage).map { sample_id, bam, bai, per_base, regions, summary_txt ->
        [sample_id, bam, summary_txt, gene_bed, bai]
    }

    DELETION_CHECK(deletion_check_input)

    //Run COVERAGE_PLOT module
    //COVERAGE_PLOT needs: (sample_id, bam, bai, gene_bed)
    //Combine BAM from BWA_MEM2 and the gene bed
    coverage_plot_input = aligned_bam.map { sample_id, bam, bai ->
        tuple(sample_id, bam, bai, gene_bed, reference)
    }

    COVERAGE_PLOT(coverage_plot_input)

    //Make input variable for MARKER_CHECK module
    // Use join() instead of combine() to match samples 1-to-1
    marker_check_input = deletion_check_input.join(DELETION_CHECK.out.deletion_csv)

    MARKER_CHECK(marker_check_input)
}
