#!/usr/bin/env nextflow
params.debug_mode = false // default

nextflow.enable.dsl = 2

include { MULTI_GENOME } from './workflows/multi_genome'
include { GENOME_ANALYSIS } from './subworkflows/local/genome_analysis'
include { PAIRWISE_COMPARISON } from './modules/local/pairwise_comparison'


// Define parameters with defaults
params.reference_genome = null // optional now
params.target_dir       = "data/target"
params.neighbour_dir    = "data/non_target"
params.output_dir = "results"
params.threads = Runtime.runtime.availableProcessors()
params.mode = "pairwise"  // Changed default to pairwise
params.subject_genome = null
params.query_genome = null
params.common_region_ident = 98
params.unique_region_ident = 85
params.min_unique_length = 50
params.max_unique_length = 0
params.min_common_length = 50
params.max_common_length = 0
params.genome_pairs = null  // Add this line
params.modules = 'all'  // Add this line
params.core_genome_ident = 98
params.unique_qcov = 0
params.unique_ident = 85
params.shared_region = 'yes'
params.shared_qcov = 90
params.shared_ident = 95
params.unique_qcov_op = "<"
params.unique_ident_op = "<"
params.shared_qcov_op = ">"
params.shared_ident_op = "<"
params.min_seq_length = 50
params.max_seq_length = 10000000
params.shared_qcov_mode = 'inside'
params.shared_ident_mode = 'inside'
params.accessory_unique = 'no'
params.help = false
params.version = false

// ------------------------------
// Banner & Help
// ------------------------------
def printBanner() {
    def banner = """
============================================================================================================================================
.oOOOo.  O       o `OooOOo.  o.OOoOoo                       OooOOo.                 
o     o  o       O  o     `o  O                             O     `O               
O.       O       o  O      O  o                             o      O o              
 `OOoo.  o       o  o     .O  ooOO                          O     .o                
      `O o       O  OOooOO'   O             ooooooooo       oOooOO'  O  .oOo. .oOo. 
       o O       O  o    o    o                             o        o  O   o OooO' 
O.    .O `o     Oo  O     O   O                             O        O  o   O O     
 `oooO'   `OoooO'O  O      o ooOooOoO                       o'       o' oOoO' `OoO' 
                                                                        O           
                                                                        o'          
============================================================================================================================================
                🧬 SURE-Pipe: Shared and Unique Region Extraction Pipeline - V.1.0.0 🧬
                
                🛠️  Mode: ${params.mode}
                ⚙️  Modules: ${params.modules ?: 'all'}
                📁 Output: ${params.output_dir}
                🧵 Threads: ${params.threads}
                🆔 Run ID: ${workflow.runName}
                🕐 Started: ${new Date().format('yyyy-MM-dd HH:mm:ss')}
                
============================================================================================================================================
"""
    log.info banner
}

// Help message
def helpMessage() {
    log.info"""
    ============================================================================================================================================
                    SURE-Pipe: Shared and Unique Region Extraction Pipeline - V.1.0.0
    ============================================================================================================================================
    Developed by Infant thomas S; Dr.Bobby paul 
    (For a commercial usage please contact : infantbiotech@gmail.com; bobby.paul@manipal.edu), Manipal School of Life Sciences, MAHE-Manipal,INDIA
    Usage:
    nextflow run main.nf [options]

    Required Options:
      --mode                    	Analysis mode: 'group' or 'pairwise' (default: $params.mode)
      --modules    			Execution mode (required) (default: $params.modules)
					Options:
					     - 'all'   : Run complete pipeline (INTRA → INTER → FILTER → FINAL ANALYSIS)
					     - 'intra' : Run only intra-genome analysis 
					     - 'inter' : Run from inter-genome to final analysis (INTER → FILTER → FINAL ANALYSIS); Run when you have only one genome in target group.
      --help                    	Display this help message
      --version                    	Show pipeline version
      --threads                 	assign the processer (defaul=$params.threads = Runtime.runtime.availableProcessors())   
    
    
    Input/Output:
      --reference_genome		Choose from target_dir [ex. data/intra_other_genomes/target_genome_1.fasta or target_genome_1.fasta] (default: $params.reference_genome)
      --target_dir          	Path to the intra genomes directory and it should contain minimum two genomes including reference (default: $params.target_dir)
      --neighbour_dir  		Path to the Neighbour_genomes directory (default: $params.neighbour_dir)
      --output_dir			Path to the output directory (default: $params.output_dir)
    
    Filtering Options:
      --min_seq_length      	Minimum length for unique & shared markers (default: $params.min_seq_length)
      --max_seq_length      	Maximum length for unique & shared markers (default: $params.max_seq_length)
      --unique_qcov         	Unique regions's maximum percentage of query coverage in inter comparsion (default: $params.unique_qcov)
      --unique_ident        	Unique regions's maximum percentage of identity in inter comparsion (default: $params.unique_ident)
      --shared_region       	Extract shared genomic regions between target & neighbour genomes: 'yes' or 'no' (default: $params.shared_region)
      --shared_qcov         	shared regions's percentage of query coverage in inter comparsion (default: $params.shared_qcov)
      --shared_ident        	shared regions's percentage of identity in inter comparsion (default: $params.shared_ident)
      --core_genome_ident		core_genome regions's percentage of identity in intra comparison (default: $params.core_genome_ident)
      --unique_qcov_op		choices=ops.keys(), help="Operator for unique region query coverage . (Choices: >, <, >=, <=; default:$params.unique_qcov_op)
      --unique_ident_op		choices=ops.keys(), help="Operator for unique region identity . (Choices: >, <, >=, <=; default:$params.unique_ident_op)
      --shared_qcov_op		choices=ops.keys(), help="Operator for shared region query coverage (will be inverted internally). (Choices: >, <, >=, <=, ==, !=; default:$params.shared_qcov_op)
      --shared_ident_op		choices=ops.keys(), help="Operator for shared region identity (will be inverted internally). (Choices: >, <, >=, <=, ==, !=; default:$params.shared_ident_op)
      --shared_qcov_mode    	filtering option for shared regions qcov is in interval (ex. 80,70; 'inside' interval or 'outside' interval)
      --shared_ident_mode   	filtering option for shared regions ident is in interval (ex. 80,70; 'inside' interval or 'outside' interval)	
      --accessory_unique    	Extract accessory & unique regions for genomes in target_group: 'yes' or 'no' (default: $params.accessory_unique)
       
    
    Pairwise Mode Options:
      --subject_genome		Path to the subject genome file for pairwise mode
      --query_genome		Path to the query genome file for pairwise mode
      --common_region_ident 		Common regions's minimum percentage of identity in pair wise mode (default: $params.common_region_ident)
      --unique_region_ident 		Unique regions's maximum percentage of identity in pair wise mode (default: $params.unique_region_ident)
      --min_unique_length		Minimum length for unique sequences (default: $params.min_unique_length) for pairwise mode
      --max_unique_length		Maximum length for unique sequences (default: $params.max_unique_length) for pairwise mode
      --min_common_length		Minimum length for common sequences (default: $params.min_common_length) for pairwise mode
      --max_common_length		Maximum length for common sequences (default: $params.max_common_length) for pairwise mode

    Nextflow options:
      -resume                   	Resume a previous run
    """.stripIndent()
}

// =====================
// Handle Help & Version
// =====================

if (params.version) {
    log.info "SURE-Pipe version 1.0.0"
    System.exit(0)
}

if (params.help || 
    (params.mode == "group" && (!params.reference_genome || !params.target_dir || !params.neighbour_dir)) ||
    (params.mode == "pairwise" && (!params.subject_genome || !params.query_genome))) {
    helpMessage()
    System.exit(0)
}

// Validate input parameters
def error_message = ""

if (params.mode == "group") {
    def missing = []
    if (!params.target_dir) missing << "--target_dir"
    if (!params.neighbour_dir) missing << "--neighbour_dir"

    if (missing) {
        error_message += "Missing required parameters for group mode: ${missing.join(', ')}\n"
    }

    if (params.target_dir && !file(params.target_dir).exists()) {
        error_message += "Target genome directory not found: ${params.target_dir}\n"
    }
    if (params.neighbour_dir && !file(params.neighbour_dir).exists()) {
        error_message += "Neighbour genome directory not found: ${params.neighbour_dir}\n"
    }

    if (!params.reference_genome) {
        log.warn "No reference genome provided. Proceeding with auto-detection."
    }

} else if (params.mode == "pairwise") {
    if (!params.subject_genome || !params.query_genome) {
        error_message += "Both --subject_genome and --query_genome must be specified for pairwise mode\n"
    }
} else {
    error_message += "Invalid mode: ${params.mode}. Use 'group' or 'pairwise'\n"
}

// Exit if validation failed
if (error_message) {
    log.error "Parameter validation failed:\n${error_message}"
    exit 1
}

// Main workflow
workflow {
    printBanner()
    log_info "Starting pipeline in ${params.mode} mode"

    if (params.mode == "group") {

        // --- Supported genome extensions ---
        def exts = ['fa','fasta','fna','ffn','frn']

        // --- Collect all target genome files that exist ---
        def target_files = []
        exts.each { ext ->
            def matches = new File(params.target_dir).listFiles({ f -> f.name.endsWith(".${ext}") } as FileFilter)
            if (matches) target_files.addAll(matches)
        }

        if (!target_files) error "No genome files found in target_dir: ${params.target_dir}"

        // --- Reference genome selection ---
        def channel_reference_genome
        def channel_target_dir

	if (params.reference_genome) {
	    // Check if user gave full path or just file name
	    def ref_file = new File(params.reference_genome)	
	    
	    if (!ref_file.exists()) {
		// Try to find it in the target_dir
		ref_file = new File(params.target_dir, params.reference_genome)
		if (!ref_file.exists()) {
		    error "Reference genome not found: ${params.reference_genome} (not in target_dir either)"
		}
	    }
	    
	    // Collect other files in target_dir excluding reference
	    def other_files = target_files.findAll { it.absolutePath != ref_file.absolutePath }

	    // Wrap as Nextflow channels
	    channel_reference_genome = Channel.of(ref_file.path).map { file(it) }
	    channel_target_dir       = Channel.from(other_files*.path).map { file(it) }
	} else {
	    // Auto-select smallest genome using seqkit
	    def existing_files = target_files*.path
	    def cmd = (["seqkit","stats","-T"] + existing_files).join(" ")
	    def output_text = cmd.execute().text.readLines().drop(1) // skip header

	    def file_sizes = output_text.collect { line ->
		def fields = line.split("\\s+")
		[file: new File(fields[0]), size: fields[4].toLong()]
	    }

	    def min_size = file_sizes*.size.min()
	    def smallest_files = file_sizes.findAll { it.size == min_size }*.file
	    def ref_file = smallest_files[0]  // pick first if multiple identical
	    def other_files = target_files.findAll { it != ref_file }

	    channel_reference_genome = Channel.of(ref_file.path).map { file(it) }
	    channel_target_dir       = Channel.from(other_files.collect { file(it) })
	}

        // --- Neighbour genomes (optional) ---
        def neighbour_files = []
        if (params.neighbour_dir) {
            exts.each { ext ->
                def matches = new File(params.neighbour_dir).listFiles({ f -> f.name.endsWith(".${ext}") } as FileFilter)
                if (matches) neighbour_files.addAll(matches)
            }
        }
        channel_neighbour_dir = neighbour_files ? Channel.from(neighbour_files.collect { file(it) }) : Channel.empty()

        // --- Modules to run ---
        def modules_to_run = params.modules?.tokenize(',') ?: []

        MULTI_GENOME(
            channel_reference_genome,
            channel_target_dir,
            channel_neighbour_dir,
            params.output_dir,
            params.core_genome_ident,
            params.unique_qcov,
            params.unique_ident,
            params.shared_qcov,
            params.shared_ident,
            params.unique_qcov_op,
            params.unique_ident_op,
            params.shared_qcov_op,
            params.shared_ident_op,
            params.shared_qcov_mode,
            params.shared_ident_mode,
            params.min_seq_length,
            params.max_seq_length,
            params.accessory_unique,    
            params.shared_region,
            modules_to_run
        )

    } else if (params.mode == "pairwise") {

        def subject_file = file(params.subject_genome)
        def query_file   = file(params.query_genome)

        if (!subject_file.exists() || !query_file.exists()) {
            error "One or both genome files do not exist."
        }

        genome_pair = Channel.of(tuple("direct_comparison", subject_file, query_file))

        PAIRWISE_COMPARISON(
            genome_pair,
            params.common_region_ident,
            params.unique_region_ident,
            params.min_unique_length,
            params.max_unique_length,
            params.min_common_length,
            params.max_common_length
        )
    }
}


// Custom logging function
def log_info(message) {
    log.info(message)
    // Ensure output directory exists
    new File(params.output_dir).mkdirs()
    file("${params.output_dir}/pipeline_execution.log").append("${new Date()} - $message\n")
}

// Workflow completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"
    log.info "Execution duration: $workflow.duration"
}
