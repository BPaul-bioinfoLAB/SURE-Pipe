process ALIGNMENT_FILTER_MODULE {
    tag "Alignment filtering process"
    debug params.debug_mode

    input:
    path core_genome_alignment_with_intraspecies
    path unaligned_core_bed_file
    val unique_qcov
    val unique_ident
    val shared_qcov
    val shared_ident
    val unique_qcov_op
    val unique_ident_op
    val shared_qcov_op
    val shared_ident_op
    val shared_qcov_mode
    val shared_ident_mode

    output:
    path "core_vs_neighbour_shared_regions.csv", emit: core_vs_neighbour_shared_regions
    path "core_vs_neighbour_unique_regions.tsv", emit: core_vs_neighbour_unique_regions
    path "core_unaligned_coordinates.bed", emit: core_unaligned_coordinates, optional: true
    path "versions.yml", emit: versions

    script:
    """
    echo "Debug: Starting ALIGNMENT_FILTER_MODULE"
    alignment_filter_module.py --input '${core_genome_alignment_with_intraspecies}' --output_dir . --ibed '${unaligned_core_bed_file}' --unique_qcov '${unique_qcov}' --unique_ident '${unique_ident}' --shared_qcov '${shared_qcov}' --shared_ident '${shared_ident}' --unique_qcov_op '${unique_qcov_op}' --unique_ident_op '${unique_ident_op}' --shared_qcov_op '${shared_qcov_op}' --shared_ident_op '${shared_ident_op}' --shared_qcov_mode '$shared_qcov_mode' --shared_ident_mode '$shared_ident_mode'

    echo "Debug: Alignment filtering process completed"
    ls -l

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
