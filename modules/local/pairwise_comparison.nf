process PAIRWISE_COMPARISON {
    tag "$pair_id"
    publishDir "${params.output_dir}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(subject_genome), path(query_genome)
    val common_region_ident
    val unique_region_ident
    val min_unique_length
    val max_unique_length
    val min_common_length
    val max_common_length

    output:
    tuple val(pair_id), path("${pair_id}_*.fa"), emit: result_fastas
    tuple val(pair_id), path("${pair_id}.db"), emit: database_dir

    script:
    """
    echo "Starting pairwise comparison"
   
    bash "${baseDir}/bin/pairwise_comparison.sh" \
        -s "${subject_genome}" \
        -q "${query_genome}" \
        -Ci "${common_region_ident}" \
        -Ui "${unique_region_ident}" \
        -u ${min_unique_length},${max_unique_length} \
        -c ${min_common_length},${max_common_length}

    echo "Pairwise comparison completed"

    # Rename output files to include pair_id
    for file in *_Filtered_Unique_sequences.fa *_Common_regions_sequences.fa; do
        mv "\$file" "${pair_id}_\$file"
    done

    # Rename the database directory
    mv *_vs_*.db "${pair_id}.db"
    """
}
