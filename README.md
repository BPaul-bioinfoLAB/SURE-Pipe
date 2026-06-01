[![Release](https://img.shields.io/badge/release-v1.1-blue.svg)](https://github.com/BPaul-bioinfoLAB/SURE-Pipe/releases/tag/v.1.1)
[![Downloads](https://img.shields.io/github/downloads/BPaul-bioinfoLAB/SURE-Pipe/v.1.1/total?logo=github)](https://github.com/BPaul-bioinfoLAB/SURE-Pipe/releases/tag/v.1.1)
	
	"""
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
# SURE-Pipe
SURE-Pipe is a Nextflow pipeline for genome comparison and analysis.

SURE-Pipe:(Shared and Unique Region Extraction Pipeline) is a Nextflow-based pipeline for genome comparison and analysis. It integrates BLAST, BEDTools, Python scripts, and Nextflow workflows to identify unique and shared genomic regions across multiple genomes.

	**Developed by:** Infant Thomas S & Dr. Bobby Paul  
	**Contact for commercial use:** infantbiotech@gmail.com, bobby.paul@manipal.edu  
	**Affiliation:** Manipal School of Life Sciences, MAHE-Manipal, INDIA

## 🔹 Features

- Identify **unique and shared genomic regions**
- Perform **intra- and inter-genome comparisons**
- Automated **Nextflow workflows** with modular scripts
- Portable **conda environment** for reproducible results
- Parallel processing support for faster execution.

## ⚡ Requirements

	- Linux or macOS
	- [Git](https://git-scm.com/)
	- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda/Anaconda)
	- Optional: [Mamba](https://mamba.readthedocs.io/) for faster environment setup

## Installation
## Clone from git-hub
	git clone https://github.com/BPaul-bioinfoLAB/SURE-Pipe.git
	cd SURE-Pipe
## Setup the environment	
	bash setup_surepipe.sh	
	chmod +x SURE-Pipe
⚠️ The setup script will install Python, pandas, BEDTools, BLAST, Biopython, SeqKit, parallel, samtools, and graphviz. Nextflow is recommended to install globally rather than inside Conda.

Install Nextflow if not already present.
Create or update the Conda environment SURE-Pipe with all dependencies.
Make all scripts executable.
Create a global symlink for SURE-Pipe.


## Check all the dependancies were properly intsalled
	Ex.
	python --version
	nextflow --version
	bedtools --version
	blastn -version

## You can add the wrapper to your $PATH to run it from anywhere (it will done by setup, if not do separately):
   ##Global execution
	sudo ln -s /path/to/SURE-Pipe/SURE-Pipe /usr/local/bin/SURE-Pipe 
## Example run_cmd
	SURE-Pipe -h
	SURE-Pipe -V
	SURE-Pipe run --mode group (options) ; See help
	SURE-Pipe run --mode group --reference_genome file_in_target_dir  --target_dir path_to_file --neigbour_dir path_to_file
	SURE-Pipe run --mode pairwise --subject_genome file_in_target_dir  --query_genome path_to_file

## 📌 Genome Architecture Compatibility

SURE-Pipe is designed to flexibly handle diverse genome architectures and input formats:

  🔹 Plasmids: Fully supported when included in the FASTA file. In our Bacillus case study, plasmids were excluded intentionally due to their frequent involvement in horizontal gene transfer (HGT), which   can confound species-level marker discovery. This reflects a biological choice, not a technical limitation.

  🔹 Multipartite genomes (e.g., multi-chromosome species): Supported without restriction. SURE-Pipe treats each chromosome or contig as an independent sequence within the same genome, and all region-based comparisons proceed normally.

  🔹 Draft / fragmented assemblies: Draft genomes in multiFASTA format are supported. However, highly fragmented assemblies may lead to shorter or interrupted shared/unique regions due to contig boundaries—an inherent limitation of alignment-based analyses rather than a constraint of SURE-Pipe.

  🔹 Input formats: Accepts common FASTA extensions (.fa, .fasta, .fna, .ffn, .frn). Each genome should be supplied as a separate FASTA file (single- or multi-FASTA), regardless of chromosome number or plasmid content. Should keep the filenames & folder names without space, (). instead use underscore(_)
  🔹 Verification: Always check the shared & unique regions by doing BLAST against curated genome databases before proceeding downstream analysis

## Parameter Selection Guidelines
## Core Genome Identity (--core_genome_ident)

🔹This parameter should be selected based on the Average Nucleotide Identity (ANI) within the target genome cluster.

## Recommended Strategy

🔹Use the lowest ANI value observed among all target genomes as the baseline threshold.

## Why?

🔹Using the lower ANI boundary captures a larger and more representative conserved core genome, improving downstream discovery of:

robust unique regions,
conserved shared regions,
and diagnostic-quality markers.

Overly strict thresholds may fragment the core genome and exclude biologically valid orthologous regions.

## Unique Region Identity (--unique_ident)
	Recommended Default --unique_ident 85

🔹An 85% identity cutoff provides the best balance between: 
			specificity,
			marker length,
			and diagnostic robustness.
🔹Lower thresholds generate overly short hyper-variable regions, while thresholds above 85% increase the risk of cross-reactivity with non-target genomes. A minimum ~15% sequence divergence (<85% identity) helps reduce false-positive primer binding.

## Shared Region Identity (`--shared_ident`)

🔹This parameter controls homologous region detection between target and neighboring genome groups.

## SURE-Pipe supports two directional modes:

| Mode | Meaning       | Biological Interpretation                |
| ---- | ------------- | ---------------------------------------- |
| `gt` | `> threshold` | Detect highly similar homologous regions |
| `lt` | `< threshold` | Detect divergent homologous regions      |

### Recommended Usage

## Detect highly conserved orthologous blocks shared between genome groups.
	--shared_ident 95 --shared_operator gt
## Detect divergent homologous regions shared between genome groups.
	--shared_ident 95 --shared_operator lt

The `gt` mode enriches for highly conserved genomic structures while filtering fragmented low-identity background alignments.
The `lt` mode captures more divergent homologous components and enables flexible exclusion of highly conserved non-target regions.

## Biological Design Principle

🔹SURE-Pipe identifies conserved core regions within the target species and compares them against neighboring genomes. Shared or conserved regions above the selected thresholds are removed, enriching for biologically stable discriminatory markers.

🔹To reduce artifacts from horizontal gene transfer, plasmids, and accessory genome variation, the workflow requires multiple genomes per target species. Only regions consistently conserved across independent genomes are retained.

## Project Directory structure
	SURE-Pipe
	├── bin
	│   ├── alignment_filter_module.py
	│   ├── inter_genome_module.sh
	│   ├── intra_genome_module.sh
	│   ├── pairwise_comparison.sh
	│   └── unique_and_shared_genome_module.sh
	├── completion
	│   ├── surepipe.bash
	├── data
	│   ├── non_target
	│   ├── pair
	│   └── target
	│   │   ├── non_target_genome_1.fasta
	│   │   ├── non_target_genome_2.fasta
	│   │   ├── non_target_genome_3.fasta
	│   │   ├── non_target_genome_4.fasta
	│   │   └── non_target_genome_5.fasta
	│   ├── pair
	│   │   ├── GCA_029857465.2_Bacillus_licheniformis_strain_DSM_26543.fa
	│   │   └── GCA_034478925.1_Bacillus_licheniformis_strain_ATCC_14580.fa
	│   └── target
	│       ├── target_genome_1.fasta
	│       ├── target_genome_2.fasta
	│       ├── target_genome_3.fasta
	│       └── target_genome_4.fasta
	├── env
	│   └── environment.yml
	├── LICENSE
	├── main.nf
	├── modules
	│   └── local
	│       ├── alignment_filter_module.nf
	│       ├── inter_genome_module.nf
	│       ├── intra_genome_module.nf
	│       ├── pairwise_comparison.nf
	│       └── unique_and_shared_genome_module.nf
	├── nextflow.config
	├── README.md
	├── results
	│   └── pipeline_trace.txt
	├── setup_surepipe.sh
	├── subworkflows
	│   └── local
	│       └── genome_analysis.nf
	├── SURE-Pipe
	├── work
	└── workflows
	    └── multi_genome.nf


## activate conda environment
	conda activate SURE-Pipe

## Help message
    ============================================================================================================================================
                    SURE-Pipe: Shared and Unique Region Extraction Pipeline - V.1.0.0
    ============================================================================================================================================
    Developed by Infant thomas S; Dr.Bobby paul 
    (For a commercial usage please contact : infantbiotech@gmail.com; bobby.paul@manipal.edu), Manipal School of Life Sciences, MAHE-Manipal,INDIA
    Usage:

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
	params.shared_qcov = 90
	params.shared_ident = 95
	params.unique_qcov_op = "<"
	params.unique_ident_op = "<"
	params.shared_qcov_op = ">"
	params.shared_ident_op = "<"
	params.min_seq_length = 0
	params.max_seq_length = 1000000000
	params.shared_qcov_mode = 'inside'
	params.shared_ident_mode = 'inside'
	params.accessory_unique = 'no'
	params.help = false
	params.version = false

    SURE-Pipe run [options]

    Required Options:
      --mode                Analysis mode: 'group' or 'pairwise' (default: $params.mode)
      --modules    			Execution mode (required) (default: $params.modules)
							Options:
							    - 'all'   : Run complete pipeline (INTRA → INTER → FILTER → FINAL ANALYSIS)
							    - 'intra' : Run only intra-genome analysis 
							    - 'inter' : Run from inter-genome to final analysis (INTER → FILTER → FINAL ANALYSIS); Run when you have only one genome in target group.
      --help                    Display this help message
      --version                 Show pipeline version
      --threads                 assign the processer (defaul=$params.threads = Runtime.runtime.availableProcessors())   
    
    
    Input/Output:
      --reference_genome		Choose from target_dir [ex. data/intra_other_genomes/target_genome_1.fasta or target_genome_1.fasta] (default: $params.reference_genome)
      --target_dir          	Path to the intra genomes directory and it should contain minimum two genomes including reference (default: $params.target_dir)
      --neighbour_dir  		    Path to the Neighbour_genomes directory (default: $params.neighbour_dir)
      --output_dir			    Path to the output directory (default: $params.output_dir)
    
    Filtering Options:
      --min_seq_length      	Minimum length for unique & shared markers (default: $params.min_seq_length)
      --max_seq_length      	Maximum length for unique & shared markers (default: $params.max_seq_length)
      --unique_qcov         	Unique regions's maximum percentage of query coverage in inter comparsion (default: $params.unique_qcov)
      --unique_ident        	Unique regions's maximum percentage of identity in inter comparsion (default: $params.unique_ident)
      --shared_region       	Extract shared genomic regions between target & neighbour genomes: 'yes' or 'no' (default: $params.shared_region) 
      --shared_qcov         	shared regions's percentage of query coverage in inter comparsion (default: $params.shared_qcov)
      --shared_ident        	shared regions's percentage of identity in inter comparsion (default: $params.shared_ident)
      --core_genome_ident		core_genome regions's percentage of identity in intra comparison (default: $params.core_genome_ident)
      --unique_qcov_op		    choices=ops.keys(), help="Operator for unique region query coverage . (Choices: >, <, >=, <=; default:$params.unique_qcov_op)
      --unique_ident_op		    choices=ops.keys(), help="Operator for unique region identity . (Choices: >, <, >=, <=; default:$params.unique_ident_op)
      --shared_qcov_op		    choices=ops.keys(), help="Operator for shared region query coverage . (Choices: >, <, >=, <=, ==, !=; default:$params.shared_qcov_op)
      --shared_ident_op		    choices=ops.keys(), help="Operator for shared region identity . (Choices: >, <, >=, <=, ==, !=; default:$params.shared_ident_op)	
      --accessory_unique    	Extract accessory & unique regions for genomes in target_group: 'yes' or 'no' (default: $params.accessory_unique)
    
    Pairwise Mode Options:
      --subject_genome		    Path to the subject genome file for pairwise mode
      --query_genome		    Path to the query genome file for pairwise mode
      --common_region_ident 	Common regions's minimum percentage of identity in pair wise mode (default: $params.common_region_ident)
      --unique_region_ident 	Unique regions's maximum percentage of identity in pair wise mode (default: $params.unique_region_ident)
      --min_unique_length		Minimum length for unique sequences (default: $params.min_unique_length) for pairwise mode
      --max_unique_length		Maximum length for unique sequences (default: $params.max_unique_length) for pairwise mode
      --min_common_length		Minimum length for common sequences (default: $params.min_common_length) for pairwise mode
      --max_common_length		Maximum length for common sequences (default: $params.max_common_length) for pairwise mode

    Nextflow options:
      -resume                   Resume a previous run
