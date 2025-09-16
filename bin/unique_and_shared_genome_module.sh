#!/bin/bash
set -euo pipefail
trap 'echo "ERROR on line $LINENO"' ERR
set -e
set -x
echo "Script started"
echo "Arguments: $@"
# Initialize variables
core_genome_fasta=""
neighbour_group=()
core_unique_out=""
OUTPUT_DIR=""
core_shared_genome=""
unaligned_core_genome_bed=""
MIN_SEQ_LENGTH=100
MAX_SEQ_LENGTH=1000000
MIN_COUNT=0
shared_region=yes
THREADS=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --core_genome_fasta)
            core_genome_fasta="$2"
            shift 2
            ;;
        --neighbour_group)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
                neighbour_group+=("$1")
                shift
            done
            ;;
        --core_unique_out)
            core_unique_out="$2"
            shift 2
            ;;
        --core_shared_genome)
            core_shared_genome="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --shared_region)
            shared_region="$2"
            shift 2
            ;;
        --unaligned_core_genome_bed)
            unaligned_core_genome_bed="$2"
            shift 2
            ;;
        --min_seq_length)
            MIN_SEQ_LENGTH="$2"
            shift 2
            ;;
        --max_seq_length)
            MAX_SEQ_LENGTH="$2"
            shift 2
            ;;
        --t)
            THREADS="$2"
            shift 2
            ;; 
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$core_genome_fasta" ||  ${#neighbour_group[@]} -eq 0 || -z "$core_unique_out" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --core_genome_fasta <file> --neighbour_group <dir> --core_unique_out <file> --output_dir <dir>"
    exit 1
fi

echo "core_genome_fasta: $core_genome_fasta"
echo "neighbour_group: ${neighbour_group[@]}"
echo "core_unique_out: $core_unique_out"
echo "OUTPUT_DIR: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Ensure required directories exist
mkdir -p "$OUTPUT_DIR/inter_master_bed_FAI" "$OUTPUT_DIR/inter_first_blast" "$OUTPUT_DIR/inter_first_blast/genome_bed"

echo "Required directories were created"

# Step 1: Generate FAI and BED files for all genomes
process_genome() {
    local genome="$1"
    echo "Processing genome: $core_genome_fasta"
    samtools faidx "$genome" || { echo "samtools faidx failed for $genome"; return 1; }
    mv "${genome}.fai" "$OUTPUT_DIR/inter_master_bed_FAI/" || { echo "Failed to move FAI file for $genome"; return 1; }
    awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$genome").fai" > "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$genome").bed" || { echo "Failed to create BED file for $genome"; return 1; }
    echo "Completed: $genome"
}

export -f process_genome
export OUTPUT_DIR

parallel -j "$THREADS" process_genome ::: "$core_genome_fasta" "${neighbour_group[@]}"
echo "FAI & BED files generated for all genomes."

# Step 2: Filter BLAST results (retain only columns 1, 2, 4, 5, 6, 7)
awk '($4 > 10) {print $1, $2, $3, $4, $5, $6, $7, $8}' OFS="\t" "$core_unique_out" > "$OUTPUT_DIR/inter_first_blast/filtered_blast_results.csv"

# Step 3: Create BED files for filtered results
awk '{print $2}' "$OUTPUT_DIR/inter_first_blast/filtered_blast_results.csv" | sort | uniq | \
parallel -j "$THREADS" '
  subject={} 
  awk -v subj="$subject" "BEGIN {OFS=\"\t\"} (\$2==subj){print \$1, (\$5<\$6?\$5:\$6), (\$5>\$6?\$5:\$6)}" "'"$OUTPUT_DIR"'/inter_first_blast/filtered_blast_results.csv" | sort -k1,1 -k2,2n | bedtools merge -i - > "'"$OUTPUT_DIR"'/inter_first_blast/${subject}.bed"'

echo "✅ Chromosomal BED files created for all subjects in parallel."


# Step 4: Concatenate BED files based on FAI files and chromosome IDs
concat_bed_for_genome() {
    local genome="$1"
    local genome_basename
    genome_basename=$(basename "$genome")
    local fai_file="$OUTPUT_DIR/inter_master_bed_FAI/${genome_basename}.fai"
    local output_bed="$OUTPUT_DIR/inter_first_blast/genome_bed/${genome_basename}.bed"
    # Skip if FAI not found
    [[ -f "$fai_file" ]] || { echo "Warning: FAI not found for $genome_basename"; return; }

    # Concatenate BED files for chromosomes
    while read -r chr_id _; do
        bed_file="$OUTPUT_DIR/inter_first_blast/${chr_id}.bed"
        [[ -f "$bed_file" ]] && cat "$bed_file"
    done < "$fai_file" | bedtools sort -i - | bedtools merge -i - > "$output_bed"

    echo "✅ Concatenated and sorted BED for $genome_basename"
}
export -f concat_bed_for_genome
export OUTPUT_DIR
parallel -j "$THREADS" concat_bed_for_genome ::: "${neighbour_group[@]}"

# Step 3: Filter common and accessory regions
shopt -s nullglob  # makes *.bed expand to empty if no match
bed_files=("$OUTPUT_DIR/intra_first_blast/genome_bed/*.bed")
num_genomes=${#bed_files[@]}
if [ "$num_genomes" -lt 2 ]; then
    echo "Only one genome BED file found ($num_genomes)."
    echo "Skipping multiintersect analysis..."
    # For the single neighbour genome
    bedtools subtract -a "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$core_genome_fasta").bed" -b $OUTPUT_DIR/inter_first_blast/genome_bed/*.bed > "$OUTPUT_DIR/inter_master_bed_FAI/Trimmed_unique_regions_$(basename "$core_genome_fasta").bed"
    bedtools getfasta -fi "$core_genome_fasta" -bed "$OUTPUT_DIR/inter_master_bed_FAI/Trimmed_unique_regions_$(basename "$core_genome_fasta").bed" -fo "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" -fullHeader
    # Check if the BED file from the third module exists
    if [[ -f "$unaligned_core_genome_bed" ]]; then
        echo "Third module BED file found: $unaligned_core_genome_bed"
    
        # Extract fasta sequences using bedtools
        bedtools getfasta -fi "$core_genome_fasta" -bed "$unaligned_core_genome_bed" -fo "$OUTPUT_DIR/inter_first_blast/sequences_queries.fasta" -fullHeader
        
        echo "Filtered sequences saved to: sequences_queries.fasta"
    else
        echo "No valid BED file found from the third module. Skipping sequence extraction."
    fi
    # Initialize concatenated FASTA file
    FINAL_FASTA="$OUTPUT_DIR/inter_first_blast/concatenated_sequences.fasta"
    > "$FINAL_FASTA"  # Create an empty file
    # Concatenate extracted FASTA files if they exist
    if [[ -f "$OUTPUT_DIR/inter_first_blast/sequences_queries.fasta" ]]; then
        cat "$OUTPUT_DIR/inter_first_blast/sequences_queries.fasta" >> "$FINAL_FASTA"
    fi

    if [[ -f "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" ]]; then
        cat "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" >> "$FINAL_FASTA"
    fi
    # Filter sequences longer than 200 bp
    awk -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" '/^>/ {if (seq && length(seq) >= min_len && length(seq) <= max_len) print header " size=" length(seq) " bp\n" seq; header=$0; seq=""} /^[^>]/ {seq = seq $0} END {if (seq && length(seq) >= min_len && length(seq) <= max_len) print header " size=" length(seq) " bp\n" seq}' "$FINAL_FASTA" > "$OUTPUT_DIR/unique_genomic_regions.fasta"
else
    echo "Multiple genome BED file found ($num_genomes)."
    echo "Processing multiintersect analysis..."
    # Step 5: Perform bedtools multiintersect
    ls "$OUTPUT_DIR"/inter_conserved_blast/genome_bed/*.bed > bedlist.txt
    bedtools multiinter -i $(cat bedlist.txt)/*.bed > "$OUTPUT_DIR/inter_first_blast/core_genome_raw.bed"
    echo "The aligned regions intersection was completed"

    # Extract only the first three columns
    awk '{print $1, $2, $3}' OFS="\t" "$OUTPUT_DIR/inter_first_blast/core_genome_raw.bed" | sort -k1,1 -k2,2n |bedtools merge -i -> "$OUTPUT_DIR/inter_first_blast/core_genome.bed"

    # Step 6: Subtract unique regions
    export OUTPUT_DIR
    # For the single reference genome
    bedtools subtract -a "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$core_genome_fasta").bed" -b "$OUTPUT_DIR/inter_first_blast/core_genome.bed" > "$OUTPUT_DIR/inter_master_bed_FAI/Trimmed_unique_regions_$(basename "$core_genome_fasta").bed"

    # Check if the BED file from the third module exists
    if [[ -f "$unaligned_core_genome_bed" ]]; then
        echo "Third module BED file found: $unaligned_core_genome_bed"
    
        # Extract fasta sequences using bedtools
        bedtools getfasta -fi "$core_genome_fasta" -bed "$unaligned_core_genome_bed" -fo "$OUTPUT_DIR/inter_first_blast/extracted_sequences_queries.fasta" -fullHeader
        
        echo "Filtered sequences saved to: extracted_sequences_queries.fasta"
    else
        echo "No valid BED file found from the third module. Skipping sequence extraction."
    fi

    # Initialize concatenated FASTA file
    FINAL_FASTA="$OUTPUT_DIR/inter_first_blast/concatenated_sequences.fasta"
    > "$FINAL_FASTA"  # Create an empty file

    # Step 7: Extract fasta sequences
    export OUTPUT_DIR
    # Extract fasta sequences with bedtools
    bedtools getfasta -fi "$core_genome_fasta" -bed "$OUTPUT_DIR/inter_master_bed_FAI/Trimmed_unique_regions_core_genome.fasta.bed" -fo "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" -fullHeader

    # Concatenate extracted FASTA files if they exist
    if [[ -f "$OUTPUT_DIR/inter_first_blast/extracted_sequences_queries.fasta" ]]; then
        cat "$OUTPUT_DIR/inter_first_blast/extracted_sequences_queries.fasta" >> "$FINAL_FASTA"
    fi

    if [[ -f "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" ]]; then
        cat "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" >> "$FINAL_FASTA"
    fi

    # Filter sequences longer than 200 bp
    awk -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" '/^>/ {if (seq && length(seq) >= min_len && length(seq) <= max_len) print header " size=" length(seq) " bp\n" seq; header=$0; seq=""} /^[^>]/ {seq = seq $0} END {if (seq && length(seq) >= min_len && length(seq) <= max_len) print header " size=" length(seq) " bp\n" seq}' $FINAL_FASTA > "$OUTPUT_DIR/unique_genomic_regions.fasta"
fi

# Step for conserved: Check if core_shared_genome exists and is non-empty
# Step for conserved: process only if shared_region=yes
if [[ "$shared_region" == "yes" ]]; then
    if [[ -s "$core_shared_genome" && $(wc -l < "$core_shared_genome") -gt 0 ]]; then
        echo "Filtered BLAST conserved file found: $core_shared_genome"

        mkdir -p "$OUTPUT_DIR/inter_conserved_blast/genome_bed"
        cp "$core_shared_genome" "$OUTPUT_DIR/inter_conserved_blast/filtered_blast_conserved_results.csv"
        echo "Filtered BLAST conserved results saved to: filtered_blast_conserved_results.csv"

        # Step 2 : create and merge BED files per subject in parallel
        awk -F'\t' '{print $2}' "$OUTPUT_DIR/inter_conserved_blast/filtered_blast_conserved_results.csv" | sort -u | \
        parallel -j "$THREADS" '
            subject={};
            awk -v subj="$subject" "BEGIN{OFS=\"\t\"} (\$2==subj){print \$1, (\$5<\$6?\$5:\$6), (\$5>\$6?\$5:\$6)}" \
            "$OUTPUT_DIR/inter_conserved_blast/filtered_blast_conserved_results.csv" | \
            sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTPUT_DIR/inter_conserved_blast/${subject}.bed"
        '

        echo "The conserved regions chromosomal BED files were created"

	process_genome() {
	    genome="$1"
	    OUTPUT_DIR="$2"

	    genome_basename=$(basename "$genome")
	    fai_file="$OUTPUT_DIR/inter_master_bed_FAI/${genome_basename}.fai"
	    tmp_concat="$OUTPUT_DIR/inter_conserved_blast/${genome_basename}_concat.bed"
	    tmp_sorted="$OUTPUT_DIR/inter_conserved_blast/${genome_basename}_sorted.bed"
	    output_bed="$OUTPUT_DIR/inter_conserved_blast/genome_bed/${genome_basename}_conserved.bed"

	    if [[ -f "$fai_file" ]]; then
		# Step 1: Concatenate
		> "$tmp_concat"
		while read -r chrom _; do
		    bed_file="$OUTPUT_DIR/inter_conserved_blast/${chrom}.bed"
		    [[ -f "$bed_file" ]] && cat "$bed_file" >> "$tmp_concat"
		done < "$fai_file"

		# Step 2: Sort
		bedtools sort -i "$tmp_concat" > "$tmp_sorted"

		# Step 3: Merge
		bedtools merge -i "$tmp_sorted" > "$output_bed"

		echo "Conserved BED file created for $genome_basename: $output_bed"
	    else
		echo "Warning: FAI file not found for $genome_basename, skipping genome."
	    fi
	}

	export -f process_genome

	# Run in parallel (adjust -j for number of CPUs)
	parallel -j 8 process_genome {} "$OUTPUT_DIR" ::: "${neighbour_group[@]}"

        
        # Count the number of chromosome BED files for this genome
        shopt -s nullglob
        bed_files=("$OUTPUT_DIR"/inter_conserved_blast/genome_bed/*_conserved.bed)
        num_bed_files=${#bed_files[@]}


        if [[ "$num_bed_files" -lt 2 ]]; then
            echo "Only one BED file found ($num_bed_files). Running the single-file logic..."
            # Find the single conserved BED file
            single_bed=$(ls "$OUTPUT_DIR/inter_conserved_blast/genome_bed/"*_conserved.bed | head -n 1)
            genome_basename=$(basename "$single_bed" "_conserved.bed")
            
            # Run bedtools getfasta and prepend genome_basename to FASTA headers
            bedtools getfasta -fi "$core_genome_fasta" -bed "$single_bed" -fo - -name+ | sed "s/^>/>${genome_basename}|/" > "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta"
            awk -v prefix="$genome_basename" -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" -v min_count="$MIN_COUNT" '/^>/ {if (seq && length(seq) >= min_len && length(seq) <= max_len && split(header,a,";") >= min_count){c=split(header,a,";");sub(/^>/,"",header);printf ">%s|%s count=%d size=%d bp\n%s\n",prefix,header,c,length(seq),seq} header=$0;seq="";next} {seq=seq $0} END{if(seq && length(seq) >= min_len && length(seq) <= max_len && split(header,a,";") >= min_count){c=split(header,a,";");sub(/^>/,"",header);printf ">%s|%s count=%d size=%d bp\n%s\n",prefix,header,c,length(seq),seq}}' "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta" > "$OUTPUT_DIR/shared_genomic_regions.fasta"

        else
            echo "Multiple BED files found ($num_bed_files). Proceeding with further processing..."
            # Step 3: Perform bedtools multiintersect
            ls "$OUTPUT_DIR"/inter_conserved_blast/genome_bed/*.bed > bedlist.txt
            bedtools multiinter -i $(cat bedlist.txt)/*.bed -header > "$OUTPUT_DIR/inter_first_blast/core_genome_raw.bed"
            bedtools multiinter -i "$OUTPUT_DIR/inter_conserved_blast/genome_bed"/*.bed -header > "$OUTPUT_DIR/inter_conserved_blast/conserved_regions_raw.bed"
            echo "The conserved regions intersection was completed"

            CONSRV_RAW="$OUTPUT_DIR/inter_conserved_blast/conserved_regions_raw.bed"
            FINAL="$OUTPUT_DIR/inter_conserved_blast/max_overlap_conserved_regions.bed"

            # Build IDX_KV_FILE only here
            IDX_KV_FILE="${PWD}/idx_kv.txt"
            # ——— PARSE HEADER FOR INDEX→FILENAME MAP —————————————————————————
            read -r hdr_line < "$CONSRV_RAW"
            IFS=$'\t' read -r -a hdr_cols <<< "$hdr_line"
            list_col=-1
            for i in "${!hdr_cols[@]}"; do
                [[ "${hdr_cols[i]}" == "list" ]] && { list_col=$i; break; }
            done
            if (( list_col < 0 )); then
                echo "ERROR: 'list' column not found in $CONSRV_RAW" >&2
                exit 1
            fi
            declare -A IDX2NAME
            idx=1
            for ((i = list_col + 1; i < ${#hdr_cols[@]}; i++)); do
                fname="${hdr_cols[i]##*/}"
                name="${fname%%.fa_conserved.bed}"
                IDX2NAME[$idx]="$name"
                echo "Mapped IDX $idx → $name"
                ((idx++))
            done
            # Write IDX_KV_FILE for AWK
            > "$IDX_KV_FILE"
            for k in "${!IDX2NAME[@]}"; do
                echo "$k:${IDX2NAME[$k]}" >> "$IDX_KV_FILE"
            done
            # Clear FINAL file
            > "$FINAL"
            # Export variables for GNU parallel
            export CONSRV_RAW FINAL IDX_KV_FILE
            # Extract unique chromosome IDs
            chroms=$(tail -n +2 "$CONSRV_RAW" | cut -f1 | sort -u)
   

            # Define function to process one chromosome
            process_chr() {
                chr="$1"
                tmp_file=$(mktemp)
                awk -v c="$chr" '$1==c' "$CONSRV_RAW" > "$tmp_file"

                awk -v OFS="\t" -v kv_file="$IDX_KV_FILE" '
                    BEGIN {
                        while ((getline line < kv_file) > 0) {
                             split(line, kv, ":")
                             idx2name[kv[1]] = kv[2]
                        }
                    }
                    {
                        split($5, ids, ",")
                        out = ""
                        for (j = 1; j <= length(ids); j++) {
                            nm = (ids[j] in idx2name) ? idx2name[ids[j]] : "MISSING_ID_" ids[j]
                            out = (out == "" ? nm : out ";" nm)
                        }
                        print $1, $2, $3, out
                    }
                ' "$tmp_file" > "$tmp_file.out"
                 echo "$tmp_file.out"

                rm -f "$tmp_file"
            }
            export -f process_chr

            # Run in parallel and collect list of temp files
	    tmp_files=$(echo "$chroms" | parallel -j "$THREADS" process_chr {})
            
            # Merge serially to final BED
	    cat $tmp_files | sort -k1,1 -k2,2n > "$FINAL"
	    
	    # Clean up temp files
	    rm -f $tmp_files

            echo "Done. Max-overlap conserved regions → $FINAL"

            # Step 5: Extract fasta sequences
            for genome in "$core_genome_fasta"; do
                bedtools getfasta -fi "$genome" -bed "$FINAL" -fo "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta" -name+
            done

            # Step 6: Append size to fasta headers
            awk -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" -v min_count="$MIN_COUNT" '/^>/ {if (seq && length(seq) >= min_len && length(seq) <= max_len && split(header,a,";") >= min_count){c=split(header,a,";"); printf "%s count=%d size=%d bp\t%s\t%d\t%d\n", header, c, length(seq), seq, c, length(seq);} header=$0; seq=""; next} {seq=seq $0} END {if (seq && length(seq) >= min_len && length(seq) <= max_len && split(header,a,";") >= min_count){c=split(header,a,";"); printf "%s count=%d size=%d bp\t%s\t%d\t%d\n", header, c, length(seq), seq, c, length(seq);}}' "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta" | sort -k3,3nr | awk -F'\t' '{print $1"\n"$2}' > "$OUTPUT_DIR/shared_genomic_regions.fasta"
        fi 
    else
        echo "WARNING: '$core_shared_genome' is not a valid .csv or is empty. Skipping conserved regions analysis."
    fi

else
    echo "Shared region analysis skipped (shared_region=$shared_region)."
fi
  echo "The unique genomics regions fasta sequneces parsed"
# Final message
  echo "Pipeline completed successfully. Output saved in $OUTPUT_DIR."
  exit 0
