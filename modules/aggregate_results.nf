process AGGREGATE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path deletion_csvs

    output:
    path "kocheck_deletions.csv"

    script:
    """
    echo "sample_id,deletion_status" > kocheck_deletions.csv
    tail -q -n +2 ${deletion_csvs} >> kocheck_deletions.csv
    """
}

