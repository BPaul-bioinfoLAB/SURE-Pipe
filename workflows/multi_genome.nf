include { GENOME_ANALYSIS } from '../subworkflows/local/genome_analysis'

workflow MULTI_GENOME {
    take:
    channel_reference_genome
    channel_target_dir
    channel_neighbour_dir
    output_dir
    core_genome_ident
    unique_qcov
    unique_ident
    shared_qcov
    shared_ident
    unique_qcov_op
    unique_ident_op
    shared_qcov_op
    shared_ident_op
    shared_qcov_mode
    shared_ident_mode
    min_seq_length
    max_seq_length
    accessory_unique
    shared_region
    modules_to_run

    main:
    GENOME_ANALYSIS(
        channel_reference_genome,
        channel_target_dir.collect(),
        channel_neighbour_dir.collect(),
        output_dir,
        core_genome_ident,
        unique_qcov,
        unique_ident,
        shared_qcov,
        shared_ident,
        unique_qcov_op,
        unique_ident_op,
        shared_qcov_op,
        shared_ident_op,
        shared_qcov_mode,
    	shared_ident_mode,
        min_seq_length,
        max_seq_length,
        accessory_unique,
        shared_region,
        modules_to_run
    )

    emit:
    unique_regions = GENOME_ANALYSIS.out.unique_regions
    shared_regions = GENOME_ANALYSIS.out.shared_regions
}
