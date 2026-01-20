process MARKER_CHECK {

    tag "${sample_id}"

    label 'samtools'

    publishDir "${params.outdir}/marker_check", mode: 'copy'

    input:
    tuple val(sample_id),
          path(bam),
          path(summary_txt),
          path(gene_bed),
          path(bai),
          path(deletion_csv) 

    output:
    path "${sample_id}.marker_check.csv"

    script:
    """
    set -euo pipefail

    MARKER_CONTIG="${params.marker_contig}"
    FLANK=${params.flank}
    MIN_MAPQ=${params.min_mapq}
    DELETION_CHECK=\$(awk -F',' 'NR==2 {print \$2}' "${deletion_csv}")
    
    # --------------------------------------------------
    # Read target gene coordinates (BED format: 0-based, half-open)
    # --------------------------------------------------
    read chrom start end < ${gene_bed}

    # Convert BED coordinates to 1-based for samtools
    # BED start is 0-based, so add 1 for 1-based
    # BED end is exclusive, so it's already the correct 1-based end
    samtools_start=\$(( start + 1 ))
    samtools_end=\${end}

    # Define window with flanking regions (work in 1-based coordinates)
    window_start=\$(( samtools_start - FLANK ))
    window_end=\$(( samtools_end + FLANK ))

    # Clamp start to >=1
    if [ \$window_start -lt 1 ]; then
        window_start=1
    fi

    # --------------------------------------------------
    # Genome-wide mean coverage (exclude marker + region rows)
    # --------------------------------------------------
    genome_mean=\$(awk '
        \$1 !~ /_region/ && \$1 != "'"\$MARKER_CONTIG"'" {
            sum += \$4; n++
        }
        END {
            if (n>0) print sum/n;
            else print 0
        }
    ' ${summary_txt})

    # --------------------------------------------------
    # Marker mean coverage
    # --------------------------------------------------
    marker_mean=\$(awk -v m="\$MARKER_CONTIG" '
        \$1 == m { print \$4 }
    ' ${summary_txt})

    marker_mean=\${marker_mean:-0}

    # --------------------------------------------------
    # Marker presence
    # --------------------------------------------------
    marker_present=\$(awk -v m=\$marker_mean 'BEGIN {
        print (m > 0.1) ? "true" : "false"
    }')

    # --------------------------------------------------
    # Marker copy number estimate
    # --------------------------------------------------
    marker_ratio=\$(awk -v m=\$marker_mean -v g=\$genome_mean 'BEGIN {
        if (g>0) print m/g;
        else print 0
    }')

    marker_copy=\$(awk -v r=\$marker_ratio 'BEGIN {
        if (r < 0.3) print "0";
        else if (r < 1.5) print "1";
        else print ">1";
    }')

    # --------------------------------------------------
    # Junction support score (paired-end mates)
    # --------------------------------------------------
    junction_support=\$(samtools view ${bam} "\$MARKER_CONTIG" | awk -v chr="\$chrom" -v start=\$window_start -v end=\$window_end -v minq=\$MIN_MAPQ '\$5 >= minq && \$7 == chr && \$8 >= start && \$8 <= end { count++ } END { print count+0 }')

    # --------------------------------------------------
    # Final classification logic
    # --------------------------------------------------
    # PASS: gene deleted, marker present in one copy, strong junction support (>10)
    # FAIL: gene still present (intact)
    # REVIEW: all other situations
    status=\$(awk -v present="\$marker_present" -v del="\$DELETION_CHECK" -v copy="\$marker_copy" -v support=\$junction_support '
        BEGIN {
            if (del == "intact") {
                print "FAIL"
            } else if (present == "true" && del == "deleted" && copy == "1" && support > 10) {
                print "PASS"
            } else {
                print "REVIEW"
            }
        }
    ')

    # --------------------------------------------------
    # Output CSV
    # --------------------------------------------------
    echo "sample_id,status,deletion_status,marker_present,marker_ratio,marker_copy,junction_support" \\
        > ${sample_id}.marker_check.csv

    echo "${sample_id},\${status},\${DELETION_CHECK},\${marker_present},\${marker_ratio},\${marker_copy},\${junction_support}" \\
        >> ${sample_id}.marker_check.csv
    """
}

