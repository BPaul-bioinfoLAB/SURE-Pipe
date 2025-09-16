process UNIQUE_AND_SHARED_GENOME_MODULE {
    tag "unique and shared genome module"
    debug params.debug_mode
    
    cpus params.threads ?: 4   // default 4 if not provided

    input:
    path core_genome_file
    path neighbour_genomes_dir
    path core_vs_neighbour_unique_regions
    path(core_unaligned_coordinates, stageAs: 'core_unaligned_coordinates.bed')
    path core_vs_neighbour_shared_regions
    val min_seq_length
    val max_seq_length
    val shared_region
    
    output:
    path "unique_genomic_regions.fasta", emit: unique_regions, optional: true
    path "shared_genomic_regions.fasta", emit: shared_regions, optional: true
    path "versions.yml", emit: versions

    script:
    def neighbour_group_list = neighbour_genomes_dir.collect { "'${it.name}'" }.join(' ')
    def unaligned_core_genome_arg = core_unaligned_coordinates.name != 'core_unaligned_coordinates.bed' ? "--final_valid_queries ${core_unaligned_coordinates}" : ''
    """
    set -e
    
    bash "${baseDir}/bin/unique_and_shared_genome_module.sh" \
        --core_genome_fasta "${core_genome_file}" \
        --neighbour_group ${neighbour_group_list} \
        --core_unique_out "${core_vs_neighbour_unique_regions}" \
        --core_shared_genome "${core_vs_neighbour_shared_regions}" \
        --t "$task.cpus" \
        --output_dir . \
        --min_seq_length ${min_seq_length} \
        --max_seq_length ${max_seq_length} \
        --shared_region ${shared_region} \
        ${unaligned_core_genome_arg}

    echo "Debug: UNIQUE_AND_SHARED_GENOME_MODULE completed"
    ls -l 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unique and shared genome module: \$(bash "${baseDir}/bin/unique_and_shared_genome_module.sh" --version 2>&1 | sed 's/^.*version //; s/Using.*\$//' || echo "unknown")
        samtools: \$(samtools --version | grep -oP 'samtools \\K[0-9.]+')
        bedtools: \$(bedtools --version | grep -oP 'bedtools v\\K[0-9.]+')
    END_VERSIONS
    """
}
