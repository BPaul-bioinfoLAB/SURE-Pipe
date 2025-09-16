process INTER_GENOME_MODULE {
    tag "Inter genome analysis"
    debug params.debug_mode
    cpus params.threads ?: 4   // default 4 if not provided

    input:
    path core_genome_file
    path neighbour_dir

    output:
    path "core_genome_alignment_with_intraspecies.tsv", emit: core_genome_alignment_with_intraspecies
    path "unaligned_core_bed_file.bed", emit: unaligned_core_bed_file
    path "versions.yml", emit: versions

    script:
    """
    echo "Debug: Starting INTER_GENOME_MODULE"
  
    inter_genome_module.sh \\
        --query ${core_genome_file} \\
        --pool ${neighbour_dir} \\
        --output . \\
        --threads '$task.cpus'
        

    echo "Debug: INTER_GENOME_MODULE completed"
    ls -l 

    echo "Debug: Generating versions.yml"
    VERSION=\$(inter_genome_module.sh --version 2>&1 || echo "unknown")
    echo "Debug: Version output: \$VERSION"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        inter_genome_module: \$VERSION
    END_VERSIONS

    echo "Debug: versions.yml content:"
    cat versions.yml
    """
}
