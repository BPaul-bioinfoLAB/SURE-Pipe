include { INTRA_GENOME_MODULE } from '../../modules/local/intra_genome_module'
include { INTER_GENOME_MODULE } from '../../modules/local/inter_genome_module'
include { ALIGNMENT_FILTER_MODULE } from '../../modules/local/alignment_filter_module'
include { UNIQUE_AND_SHARED_GENOME_MODULE } from '../../modules/local/unique_and_shared_genome_module'

workflow GENOME_ANALYSIS {
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
    def resolved_modules = []

    if (modules_to_run.contains('all')) {
        resolved_modules = ['all']
        log.info "Mode: ALL - Running complete pipeline (INTRA → INTER → FILTER → ANALYSIS)"
    } else if (modules_to_run.contains('intra')) {
        resolved_modules = ['intra']
        log.info "Mode: INTRA - Running intra-genome analysis only"
    } else if (modules_to_run.contains('inter')) {
        resolved_modules = ['inter', 'filter', 'analysis']
        log.info "Mode: INTER - Running from inter to analysis (INTER → FILTER → ANALYSIS)"
    } else {
        error "Invalid module specified: '${modules_to_run.join(', ')}'. Valid options: 'all', 'intra', or 'inter'"
    }

    // Default empty channels
    core_genome_file = Channel.empty()
    core_genome_alignment_with_intraspecies = Channel.empty()
    unaligned_core_bed_file = Channel.empty()
    core_vs_neighbour_unique_regions = Channel.empty()
    core_unaligned_coordinates = Channel.empty()
    core_vs_neighbour_shared_regions = Channel.empty()
    unique_regions = Channel.empty()
    shared_regions = Channel.empty()
    accessory_genome_file = Channel.empty()
    intra_unique_file = Channel.empty()

    // -------------------------
    // Module 1: Intra-genome
    // -------------------------
    if (resolved_modules.contains('all') || resolved_modules.contains('intra')) {
        INTRA_GENOME_MODULE(
            channel_reference_genome, 
            channel_target_dir,  
            core_genome_ident, 
            accessory_unique
        )
        core_genome_file = INTRA_GENOME_MODULE.out.core_genome
    } else {
        // Skip INTRA, use reference genome directly
        core_genome_file = channel_reference_genome
        log.info "✓ Using reference genome directly (INTRA skipped)"
    }

    // -------------------------
    // Module 2: Inter-genome
    // -------------------------
    if (resolved_modules.contains('all') || resolved_modules.contains('inter')) {
        INTER_GENOME_MODULE(core_genome_file, channel_neighbour_dir)
        core_genome_alignment_with_intraspecies = INTER_GENOME_MODULE.out.core_genome_alignment_with_intraspecies
        unaligned_core_bed_file = INTER_GENOME_MODULE.out.unaligned_core_bed_file
    }

    // -------------------------
    // Module 3: Filter BLAST
    // -------------------------
    if (resolved_modules.contains('all') || resolved_modules.contains('filter')) {
        ALIGNMENT_FILTER_MODULE(
            core_genome_alignment_with_intraspecies,
            unaligned_core_bed_file,
            unique_qcov,
            unique_ident,
            shared_qcov,
            shared_ident,
            unique_qcov_op,
            unique_ident_op,
            shared_qcov_op,
            shared_ident_op,
            shared_qcov_mode,
            shared_ident_mode
        )
        core_vs_neighbour_unique_regions = ALIGNMENT_FILTER_MODULE.out.core_vs_neighbour_unique_regions
        core_unaligned_coordinates = ALIGNMENT_FILTER_MODULE.out.core_unaligned_coordinates
        core_vs_neighbour_shared_regions = ALIGNMENT_FILTER_MODULE.out.core_vs_neighbour_shared_regions
    }

    // -------------------------
    // Module 4: Final analysis
    // -------------------------
    if (resolved_modules.contains('all') || resolved_modules.contains('analysis')) {
        UNIQUE_AND_SHARED_GENOME_MODULE(
            core_genome_file,
            channel_neighbour_dir,
            core_vs_neighbour_unique_regions,
            core_unaligned_coordinates,
            core_vs_neighbour_shared_regions,
            min_seq_length,
            max_seq_length,
            shared_region
        )
        unique_regions = UNIQUE_AND_SHARED_GENOME_MODULE.out.unique_regions
        shared_regions = UNIQUE_AND_SHARED_GENOME_MODULE.out.shared_regions
    }

    // -------------------------
    // Final process to organize results
    // -------------------------
    // Collect all available output files
    all_files = Channel.empty()
        .mix(core_genome_file)
        .mix(accessory_genome_file.ifEmpty([]))
        .mix(intra_unique_file.ifEmpty([]))
        .mix(unique_regions.ifEmpty([]))
        .mix(shared_regions.ifEmpty([]))
        .collect()

    PROCESS_FINAL_RESULTS(
        all_files,
        output_dir
    )

    emit:
    unique_regions = unique_regions
    shared_regions = shared_regions
}

// -------------------------
// Organize final outputs
// -------------------------
process PROCESS_FINAL_RESULTS {
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    path(files)
    val(output_dir)

    output:
    path("intra-genome-analysis/*"), optional: true
    path("inter-genome-analysis/*"), optional: true

    script:
    """
    mkdir -p intra-genome-analysis
    mkdir -p inter-genome-analysis

    # Debug: List all input files
    echo "Input files received:"
    for file in ${files}; do
        echo "  - \$file"
    done

    # Process each file based on its name/type
    for file in ${files}; do
        echo "Processing: \$file"
        
        # Intra-genome files
        if [[ "\$file" == *"core_genome"* ]]; then
            echo "  -> Copying to intra-genome-analysis/"
            cp "\$file" intra-genome-analysis/
        elif [[ "\$file" == *"accessory_genome"* ]]; then
            echo "  -> Copying to intra-genome-analysis/"
            cp "\$file" intra-genome-analysis/
        elif [[ "\$file" == *"intra_unique_"* ]]; then
            echo "  -> Copying to intra-genome-analysis/"
            cp "\$file" intra-genome-analysis/
            
        # Inter-genome files - FIXED PATTERNS
        elif [[ "\$file" == *"unique_genomic_regions"* ]]; then
            echo "  -> Copying to inter-genome-analysis/"
            cp "\$file" inter-genome-analysis/
        elif [[ "\$file" == *"shared_genomic_regions"* ]]; then
            echo "  -> Copying to inter-genome-analysis/"
            cp "\$file" inter-genome-analysis/
        else
            echo "  -> No matching pattern for \$file"
        fi
    done
    
    # Show what was copied
    echo "Final results:"
    echo "Intra-genome analysis:"
    ls -la intra-genome-analysis/ || echo "  (empty)"
    echo "Inter-genome analysis:"
    ls -la inter-genome-analysis/ || echo "  (empty)"
    """
}
