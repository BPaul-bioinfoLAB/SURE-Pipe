#!/bin/bash
set -e

echo "Script started"
echo "Arguments: $@"

# Define input files and output directory
core_genome_file=""
neighbour_genomes_dir=""
OUTPUT_DIR=""
THREADS=""
POOL_FILES=()

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --query)
            core_genome_file="$2"
            shift 2
            ;;
        --pool)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
                if [[ -d "$1" ]]; then
                    neighbour_genomes_dir="$1"
                else
                    POOL_FILES+=("$1")
                fi
                shift
            done
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Step 1: Check required arguments
if [[ -z "$core_genome_file" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --query <core_genome_file> --pool <genome_files_directory|genome_files...> --output <output_dir> [--threads <num_threads>]"
    exit 1
fi

# Step 2: Make output directory
mkdir -p "$OUTPUT_DIR"

# Step 3: Collect genome files
if [[ -n "$neighbour_genomes_dir" ]]; then
    echo "🔍 Pool is a directory: $neighbour_genomes_dir"
    genome_files=$(find "$neighbour_genomes_dir" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.ffn" -o -name "*.frn" \))
else
    echo "🔍 Pool is a list of files: ${POOL_FILES[@]}"
    genome_files="${POOL_FILES[@]}"
fi

if [[ -z "$genome_files" ]]; then
    echo "❌ Error: No genome files found!"
    exit 1
fi

echo "📊 Found genome files:"
echo "$genome_files"

# Step 4: Run BLAST
run_blast_single() {
    local subject="$1"
    local core_genome_file="$2"
    local outdir="$3"

    local subj_basename=$(basename "$subject")
    local out_file="$outdir/blast_${subj_basename}.tsv"

    blastn -query "$core_genome_file" -subject "$subject" -out "$out_file" \
        -outfmt "7 qseqid sseqid stitle qcovs pident qstart qend sstart send length evalue bitscore" \
        -task megablast -subject_besthit

    echo "✅ BLAST completed for $subj_basename"
}

export -f run_blast_single
export core_genome_file
export OUTPUT_DIR

echo "$genome_files" | tr ' ' '\n' | parallel -j "$THREADS" run_blast_single {} "$core_genome_file" "$OUTPUT_DIR"

echo "🔹 All BLAST runs completed."

# Step 5: Merge BLAST outputs
merged_output="$OUTPUT_DIR/merged_blast.tsv"
awk 'FNR>0{print}' "$OUTPUT_DIR"/blast_*.tsv > "$merged_output"

awk '/^# Query:/ {split($3,a, "[:-]"); chr=a[1]; start=a[2]-1; end=a[3]; name=$3} /^# 0 hits found/ {print chr "\t" start "\t" end "\t" name}' "$merged_output" > "$OUTPUT_DIR/unaligned_core_bed_file.bed"

grep -v '^#' "$merged_output" > "$OUTPUT_DIR/core_genome_alignment_with_intraspecies.tsv"

echo "✅ BLAST processing complete. Output saved as '$OUTPUT_DIR/core_genome_alignment_with_intraspecies.tsv'."
