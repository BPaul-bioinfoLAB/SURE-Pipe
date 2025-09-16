#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process INTRA_GENOME_MODULE {
    tag "Intra genome comparison"
    debug params.debug_mode
    
    cpus params.threads ?: 4   // default 4 if not provided
    
    input:
    path reference_genome
    path target_dir
    val core_genome_ident
    val accessory_unique
    
    output:
    path "core_genome.fasta", emit: core_genome
    path "accessory_genome.fasta", emit: accessory_genome, optional: true
    path "intra_unique_*", emit: intra_unique_genome, optional: true
    path "versions.yml", emit: versions

    script:
    def target_genomes_list = target_dir.collect { "'$it'" }.join(' ')
    """
    echo "Debug: Starting INTRA_GENOME_MODULE"
    
    intra_genome_module.sh \\
        -m '${reference_genome}' \\
        -n ${target_genomes_list} \\
        -o '.' \\
        -Ci '${core_genome_ident}' \\
        -t $task.cpus \\
        -A '${accessory_unique}'
    echo "Debug: INTRA_GENOME_MODULE completed"
    ls -l
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        intra_genome_module: \$(intra_genome_module.sh --version 2>&1 | sed 's/^.*version //; s/Using.*\$//' || echo "unknown")
    END_VERSIONS
    """
}
