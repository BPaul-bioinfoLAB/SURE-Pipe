#!/bin/bash
set -e
echo "Script started"
# Define input files and output directory
target_ref=""
target_group=()
OUTPUT_DIR=""
core_genome_ident=""
accessory_unique=""
THREADS=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m)
            target_ref="$2"
            shift 2
            ;;
        -n)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
                target_group+=("$1")
                shift
            done
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -Ci)
            core_genome_ident="$2"
            shift 2
            ;;
        -A)
            accessory_unique="$2"
            shift 2
            ;;
        -t)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Ensure output directory is not empty
if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Output directory not specified"
    exit 1
fi

# Normalize accessory option (convert to lowercase)
accessory_unique=$(echo "$accessory_unique" | tr '[:upper:]' '[:lower:]')
echo "Accessory unique: $accessory_unique"

# Remove any trailing slashes from OUTPUT_DIR
OUTPUT_DIR="${OUTPUT_DIR%/}"

echo "Final OUTPUT_DIR: $OUTPUT_DIR"

# Ensure output directory exists and is writable
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create output directory: $OUTPUT_DIR"; exit 1; }

echo "Output directory created: $OUTPUT_DIR"

# Ensure required subdirectories exist
for dir in "intra_master_bed_FAI" "intra_first_blast" "intra_2nd_blast" "intra_first_blast/genome_bed" "intra_2nd_blast/genome_bed"; do
    mkdir -p "${OUTPUT_DIR}/${dir}" || { echo "Failed to create directory: ${OUTPUT_DIR}/${dir}"; exit 1; }
done

echo "Required directories were created"

# Step 1: Generate FAI and BED files for all genomes
process_genome() {
  local genome="$1"
  echo "Processing genome: $genome"
  samtools faidx "$genome" &&
  mv "${genome}.fai" "$OUTPUT_DIR/intra_master_bed_FAI/" &&
  awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' \
    "$OUTPUT_DIR/intra_master_bed_FAI/$(basename "$genome").fai" \
    > "$OUTPUT_DIR/intra_master_bed_FAI/$(basename "$genome").bed"
}
export -f process_genome
export OUTPUT_DIR

parallel -j "$THREADS" process_genome ::: "$target_ref" "${target_group[@]}"
echo "FAI & BED file conversion were completed"

# Step 2: Function to run pairwise BLAST and generate BED
run_pairwise_blast() {
    local query="$1"
    local subject="$2"
    local outdir="$3"   # now this will have a proper value
    local query_basename=$(basename "$query")
    local subject_basename=$(basename "$subject")

    mkdir -p "$outdir/intra_first_blast" "$outdir/intra_first_blast/genome_bed"

    local blast_out="$outdir/intra_first_blast/${query_basename}_vs_${subject_basename}.csv"
    local filtered_out="$outdir/intra_first_blast/${query_basename}_vs_${subject_basename}_filtered.csv"

    blastn -query "$query" -subject "$subject" \
        -out "$blast_out" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs qcovhsp"

    awk -v cg_ident="$core_genome_ident" '($4 > 300 && $3 > cg_ident )' "$blast_out" > "$filtered_out"

    awk -v subj="$subject_basename" 'BEGIN {OFS="\t"} {print $1, ($5<=$6?$5:$6), ($5>=$6?$5:$6), subj}' "$filtered_out" \
        | sort -k1,1 -k2,2n > "$outdir/intra_first_blast/genome_bed/${subject_basename}.bed"
}


export -f run_pairwise_blast
export OUTPUT_DIR

# Parallel execution: all queries x all subjects
parallel -j "$THREADS" run_pairwise_blast "$target_ref" {} "$OUTPUT_DIR" ::: "${target_group[@]}"

echo "All pairwise BLAST & BED files generated successfully."

# Step 3: Filter common and accessory regions
num_genomes=$(ls "$OUTPUT_DIR/intra_first_blast/genome_bed"/*.bed | wc -l)

if [ "$num_genomes" -lt 2 ]; then
    echo "Only one genome BED file found ($num_genomes)."
    echo "Skipping Accessory region analysis..."
    bedtools getfasta -fi "$target_ref" -bed "$OUTPUT_DIR/intra_first_blast/genome_bed"/*.bed -fo "$OUTPUT_DIR/core_genome.fasta"
    if [[ "$accessory_unique" == "yes" ]]; then
        echo "Accessory genome analysis enabled..."
        # Extract common regions from target_ref
        subject_BED_basename=$(basename "$OUTPUT_DIR"/intra_first_blast/genome_bed/*.bed .bed)
        awk 'BEGIN {OFS="\t"} {print $1, ($5<=$6?$5:$6), ($5>=$6?$5:$6)}' "$OUTPUT_DIR"/intra_first_blast/*_filtered.csv | sort -k1,1 -k2,2n  > "$OUTPUT_DIR/intra_2nd_blast/genome_bed/${subject_BED_basename}.bed"
       
        for genome in "${target_group[@]}"; do
            bedtools subtract -a "$OUTPUT_DIR/intra_master_bed_FAI/$(basename "$genome").bed" -b "$OUTPUT_DIR/intra_2nd_blast/genome_bed/${subject_BED_basename}.bed" > "$OUTPUT_DIR/intra_2nd_blast/Unique_$(basename "$genome").bed"
            bedtools getfasta -fi "$genome" -bed "$OUTPUT_DIR/intra_2nd_blast/Unique_$(basename "$genome").bed" -fo "$OUTPUT_DIR/unique_$(basename "$genome")"
        done
        genome_basename=$(basename "$target_ref")
        bedtools subtract -a "$OUTPUT_DIR/intra_master_bed_FAI/${genome_basename}.bed" -b "$OUTPUT_DIR/intra_first_blast/genome_bed/${subject_BED_basename}.bed" > "$OUTPUT_DIR/intra_2nd_blast/Unique_${genome_basename}.bed"
        bedtools getfasta -fi "$target_ref" -bed "$OUTPUT_DIR/intra_2nd_blast/Unique_${genome_basename}.bed" -fo "$OUTPUT_DIR/unique_${genome_basename}"
        echo "Accessory/Unique genome extraction completed."
    else
        echo "Skipping accessory/unique analysis (accessory_unique=$accessory_unique)."
    fi
else
    echo "$num_genomes genome BED files found."
    echo "Proceeding with Accessory region analysis..."
    # Step 3: Perform bedtools multiintersect
    bedtools multiinter -i "$OUTPUT_DIR/intra_first_blast/genome_bed"/*.bed > "$OUTPUT_DIR/intra_first_blast/common_regions_raw.bed"
    echo "The aligned regions intersection was completed"
    awk -v num_genomes="$num_genomes" -v common_regions="$OUTPUT_DIR/intra_first_blast/common_regions_unmerged.bed" -v accessory_genome="$OUTPUT_DIR/intra_first_blast/accessory_genome_unmerged.bed" 'BEGIN {OFS="\t"} {if ($4 == num_genomes) print > common_regions; else print > accessory_genome}' "$OUTPUT_DIR/intra_first_blast/common_regions_raw.bed"
    awk 'OFS="\t" {print $1, $2, $3}' "$OUTPUT_DIR/intra_first_blast/common_regions_unmerged.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTPUT_DIR/intra_first_blast/common_regions.bed"
    bedtools getfasta -fi "$target_ref" -bed "$OUTPUT_DIR/intra_first_blast/common_regions.bed" -fo "$OUTPUT_DIR/core_genome.fasta"
    if [[ "$accessory_unique" == "yes" ]]; then
       echo "Accessory genome analysis enabled..."
       # Step 4: Extract fasta sequences
       awk 'OFS="\t" {print $1, $2, $3}' "$OUTPUT_DIR/intra_first_blast/accessory_genome_unmerged.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTPUT_DIR/intra_first_blast/accessory_genome.bed"
       bedtools getfasta -fi "$target_ref" -bed "$OUTPUT_DIR/intra_first_blast/accessory_genome.bed" -fo "$OUTPUT_DIR/accessory_genome.fasta"
       echo "The common and accessory  genomics regions fasta sequneces parsed"
       # Step 5: Concatenate common & accessory_genome.fasta
       cat "$OUTPUT_DIR"/core_genome.fasta "$OUTPUT_DIR"/accessory_genome.fasta > "$OUTPUT_DIR/intra_first_blast/concatenated_common_accessory.fasta"
       # Step 6: Function for pairwise BLAST
       run_pairwise_blast_2() {
           local query="$1"
           local subject="$2"
           local outdir="$3"
           local query_basename
           local subject_basename

           query_basename=$(basename "$query")
           subject_basename=$(basename "$subject")

           mkdir -p "$outdir/intra_2nd_blast/genome_bed"

           local blast_out="$outdir/intra_2nd_blast/${query_basename}_vs_${subject_basename}.csv"
           local filtered_out="$outdir/intra_2nd_blast/${query_basename}_vs_${subject_basename}_filtered.csv"

           # Run BLAST directly without a database
           blastn -query "$query" -subject "$subject" -out "$blast_out" -outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs qcovhsp"

           # Filter hits:  >98% identity
           awk '($3 > 98)' "$blast_out" > "$filtered_out"

           # Convert to BED format
           awk -v subj="$subject_basename" 'BEGIN {OFS="\t"} {print $2, ($7 < $8 ? $7 : $8), ($7 > $8 ? $7 : $8)}' "$filtered_out" | sort -k1,1 -k2,2n | bedtools merge > "$outdir/intra_2nd_blast/genome_bed/${subject_basename}.bed"
        }

        export -f run_pairwise_blast_2
        export OUTPUT_DIR

        # Step 7: Run pairwise BLAST in parallel
        query_file="$OUTPUT_DIR/intra_first_blast/concatenated_common_accessory.fasta"
        parallel -j "$THREADS" run_pairwise_blast_2 "$query_file" {} "$OUTPUT_DIR" ::: "$target_ref" "${target_group[@]}"


        echo "All pairwise BLAST & BED files generated successfully."

        # Step 15: Function to subtract unique regions
        subtract_unique_regions() {
            local genome="$1"
            local outdir="$2"
            local genome_basename=$(basename "$genome")

            bedtools subtract -a "$outdir/intra_master_bed_FAI/${genome_basename}.bed" -b "$outdir/intra_2nd_blast/genome_bed/${genome_basename}.bed" > "$outdir/intra_2nd_blast/Unique_${genome_basename}.bed"
        }

        export -f subtract_unique_regions
        export OUTPUT_DIR

        # Run unique region subtraction in parallel
        parallel -j "$THREADS" subtract_unique_regions {} "$OUTPUT_DIR" ::: "$target_ref" "${target_group[@]}"

        echo "Unique regions BED files created in parallel."
        extract_unique_fasta() {
            local genome="$1"
            local bed_file="$OUTPUT_DIR/intra_2nd_blast/Unique_$(basename "$genome").bed"
            local output_file="$OUTPUT_DIR/intra_unique_$(basename "$genome")"
            echo "Extracting unique sequences from $genome → $output_file"
            bedtools getfasta -fi "$genome" -bed "$bed_file" -fo "$output_file"
            echo "Completed processing $genome"
        }
        export -f extract_unique_fasta
        export OUTPUT_DIR
        parallel -j "$THREADS" extract_unique_fasta {} ::: "$target_ref" "${target_group[@]}"

        echo "Unique regions sequnces generated for all genomes."
    else
        echo "Accessory genome analysis skipped (accessory_unique=no)."
        echo "Pipeline completed with common regions only."
        exit 0
    fi
fi

