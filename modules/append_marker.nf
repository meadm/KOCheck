process APPEND_MARKER {

    tag "append_marker"

    publishDir "${params.outdir}/append_marker", mode: 'copy'

    input:
    path reference_fasta
    path marker_fasta

    output:
    path "reference_with_marker.fasta", emit: combined_ref
    path "reference_with_marker.fasta.fai", emit: combined_fai

    script:
    """
    set -euo pipefail

    # Concatenate reference and marker into a single FASTA
    cat ${reference_fasta} ${marker_fasta} > reference_with_marker.fasta

    # Index the combined FASTA for downstream tools
    samtools faidx reference_with_marker.fasta
    """
}

