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


# -----------------------------
# STEP 4: Build BLAST database
# -----------------------------
echo "🧬 Creating combined BLAST FASTA..."
cat $genome_files > "$OUTPUT_DIR/inter_genome_db.fasta"

echo "🧬 Building BLAST database..."
makeblastdb \
    -in "$OUTPUT_DIR/inter_genome_db.fasta" \
    -dbtype nucl \
    -out "$OUTPUT_DIR/inter_genome_db"


# -----------------------------
# STEP 5: Run BLAST
# -----------------------------
echo "🚀 Running BLAST against combined database..."

blastn \
    -query "$core_genome_file" \
    -db "$OUTPUT_DIR/inter_genome_db" \
    -out "$OUTPUT_DIR/merged_blast.tsv" \
    -outfmt "7 qseqid sseqid stitle qcovs pident qstart qend sstart send length evalue bitscore" \
    -task megablast \
    -num_threads "$THREADS"

# -----------------------------
# STEP 6: Extract unaligned core regions
# -----------------------------
merged_output="$OUTPUT_DIR/merged_blast.tsv"

awk '/^# Query:/ {split($3,a, "[:-]"); chr=a[1]; start=a[2]-1; end=a[3]; name=$3} /^# 0 hits found/ {print chr "\t" start "\t" end "\t" name}' "$merged_output" > "$OUTPUT_DIR/unaligned_core_bed_file.bed"

grep -v '^#' "$merged_output" > "$OUTPUT_DIR/core_genome_alignment_with_intraspecies.tsv"

echo "✅ BLAST processing complete. Output saved as '$OUTPUT_DIR/core_genome_alignment_with_intraspecies.tsv'."
