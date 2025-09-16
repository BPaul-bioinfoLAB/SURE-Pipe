	#!/bin/bash

	# Default values for sequence length thresholds
	MIN_UNIQUE_LENGTH=50
	MAX_UNIQUE_LENGTH=0  # Default max 0 indicates no upper limit
	MIN_COMMON_LENGTH=50
	MAX_COMMON_LENGTH=0  # Default max 0 indicates no upper limit
	QUERY_GENOME=""
	SUB_GENOME=""
	common_region_ident=98
	unique_region_ident=85

	# Parse arguments for custom min/max lengths
	while [[ $# -gt 0 ]]; do
	  case $1 in
	    -q) 
	      QUERY_GENOME="$2"
	      shift 2
	      ;;
	    -s) 
	      SUB_GENOME="$2"
	      shift 2
	      ;;	      
	    -Ci) 
	      common_region_ident="$2"
	      shift 2
	      ;;
	    -Ui) 
	      unique_region_ident="$2"
	      shift 2
	      ;;	      
	    -u) 
	      IFS=',' read -r MIN_UNIQUE_LENGTH MAX_UNIQUE_LENGTH <<< "$2"
	      # Remove leading/trailing whitespace from values
	      MIN_UNIQUE_LENGTH=$(echo "$MIN_UNIQUE_LENGTH" | xargs)
	      MAX_UNIQUE_LENGTH=$(echo "$MAX_UNIQUE_LENGTH" | xargs)
	      echo "Parsed unique lengths: min = $MIN_UNIQUE_LENGTH, max = $MAX_UNIQUE_LENGTH"
	      shift 2
	      ;;
	    -c)
	      IFS=',' read -r MIN_COMMON_LENGTH MAX_COMMON_LENGTH <<< "$2"
	      # Remove leading/trailing whitespace from values
	      MIN_COMMON_LENGTH=$(echo "$MIN_COMMON_LENGTH" | xargs)
	      MAX_COMMON_LENGTH=$(echo "$MAX_COMMON_LENGTH" | xargs)
	      echo "Parsed COMMON lengths: min = $MIN_COMMON_LENGTH, max = $MAX_COMMON_LENGTH"
	      shift 2
	      ;;
	    -h)
	      echo "Usage: $0 -q <query_genome.fa> -s <subject_genome.fa> [-Ci ident] [-Ui ident] [-u min,max] [-c min,max]"
      	      echo "-q : Query genome fasta"
              echo "-s : Subject genome fasta"
              echo "-Ci: Common region identity cutoff"
              echo "-Ui: Unique region identity cutoff"
              echo "-u : Unique sequence length min,max (default 50,0)"
              echo "-c : Common sequence length min,max (default 50,0)"
              exit 0
	      ;;
	    *)
	      echo "Invalid option:$1"
	      exit 1
	      ;;
	  esac
	done

	if [[ -z "$QUERY_GENOME" || -z "$SUB_GENOME" ]]; then
            echo "Error: Both -q <query.fa> and -s <subject.fa> must be provided."
            exit 1
        fi

	# Debug output to check remaining arguments
	echo "Subject genome: $SUB_GENOME"
	echo "Query genome: $QUERY_GENOME"

	# Extract base names for output files
	SUB_BASE="${SUB_GENOME%.fa}"
	QUERY_BASE="${QUERY_GENOME%.fa}"

	# Function to convert FASTA to BED format using shell
	fasta_to_bed() {
	    local fasta_file=$1
	    local bed_file=$2
	    echo "Indexing $fasta_file to generate .fai file..."

	    # Generate .fai index file
	    samtools faidx "$fasta_file"
	    if [[ $? -ne 0 ]]; then
	        echo "Error: Failed to create .fai file for $fasta_file"
	        exit 1
	    fi

	    # Convert .fai file to BED format
	    echo "Converting .fai file to BED format for $fasta_file..."
	    awk -v OFS='\t' '{print $1, 0, $2}' "${fasta_file}.fai" > "$bed_file"

	    echo "BED file created: $bed_file"
	}

	# Get file names for output BED files based on input FASTA file names
	BED1="${SUB_BASE}.bed"
	BED2="${QUERY_BASE}.bed"

	# Step 1: Convert first genome to BED using .fai file
	fasta_to_bed "$SUB_GENOME" "$BED1"

	# Step 2: Convert second genome to BED using .fai file
	fasta_to_bed "$QUERY_GENOME" "$BED2"

	echo "FASTA to BED conversion complete for both genomes."

	# Check if the reference and query genome files exist
	if [[ ! -f "$SUB_GENOME" || ! -f "$QUERY_GENOME" ]]; then	
	    echo "Reference or query genome file does not exist."
	    exit 1
	fi

	# Step 4: Perform BLAST search with the second genome
	blastn -query "$QUERY_GENOME" -subject "${SUB_GENOME}" -out "$SUB_GENOME"_vs_"$QUERY_GENOME"_blast_results.txt -outfmt '6 qseqid sseqid pident length qstart qend sstart send qcovs qcovhsp'

	# Step 5: Filter results for >= 98% identity
	awk -v Ci="$common_region_ident" '$3 >= Ci' "$SUB_GENOME"_vs_"$QUERY_GENOME"_blast_results.txt > "$SUB_GENOME"_vs_"$QUERY_GENOME"_common_regions.txt
	# Step 5a: Filter results for >= 85% identity
	awk -v Ui="$unique_region_ident" '$3 >= Ui' "$SUB_GENOME"_vs_"$QUERY_GENOME"_blast_results.txt > "$SUB_GENOME"_vs_"$QUERY_GENOME"_unique_regions.txt

	# Step 6: Extract aligned regions from the filtered results for the query, including strand
	awk '{
	    start = ($5 < $6 ? $5 : $6);
	    end = ($5 > $6 ? $5 : $6);
	    print $1 "\t" start "\t" end
	}' "$SUB_GENOME"_vs_"$QUERY_GENOME"_unique_regions.txt | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -> "${QUERY_BASE}_unique_regions.bed"

	# Step 6a: Extract aligned regions from the filtered results for the query, including strand
	awk '{
	    start = ($5 < $6 ? $5 : $6);
	    end = ($5 > $6 ? $5 : $6);
	    print $1 "\t" start "\t" end
	}' "$SUB_GENOME"_vs_"$QUERY_GENOME"_common_regions.txt | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -> "${QUERY_BASE}_Common_regions.bed"

	# Step 7: Extract aligned regions from the filtered results for the subject, including strand
	awk '{
	    start = ($7 < $8 ? $7 : $8);
	    end = ($7 > $8 ? $7 : $8);
	    print $2 "\t" start "\t" end
	}' "$SUB_GENOME"_vs_"$QUERY_GENOME"_unique_regions.txt | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -> "${SUB_BASE}_unique_regions.bed"

	# Step 7a: Extract aligned regions from the filtered results for the subject, including strand
	awk '{
	    start = ($7 < $8 ? $7 : $8);
	    end = ($7 > $8 ? $7 : $8);
	    print $2 "\t" start "\t" end
	}' "$SUB_GENOME"_vs_"$QUERY_GENOME"_common_regions.txt | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i -> "${SUB_BASE}_Common_regions.bed"

	# Step 8: Subtract aligned regions from the query and save to a new file
	bedtools subtract -a "$BED2" -b "${QUERY_BASE}_unique_regions.bed" > "${QUERY_BASE}_Unique.bed"

	# Step 9: Subtract aligned regions from the reference and save to a new file
	bedtools subtract -a "$BED1" -b "${SUB_BASE}_unique_regions.bed" > "${SUB_BASE}_Unique.bed"

	echo "Subtraction complete. Results saved to ${QUERY_BASE}_Unique.bed and ${SUB_BASE}_Unique.bed."

	# Step 10: Get FASTA sequences for Unique.bed from the query genome and filter for user-specified length
	bedtools getfasta -fi "$QUERY_GENOME" -bed "${QUERY_BASE}_Unique.bed" -fo "${QUERY_BASE}_Unique_sequences.fa" -fullHeader -s
	awk -v min_len="$MIN_UNIQUE_LENGTH" -v max_len="$MAX_UNIQUE_LENGTH" 'BEGIN { RS=">"; FS="\n"; seq_num = 0 } NR > 1 { header = $1; seq = ""; for (i = 2; i <= NF; i++) seq = seq $i; if (length(seq) >= min_len && (max_len == 0 || length(seq) <= max_len)) { seq_num++; printf ">%d:%s (%d bp)\n%s\n", seq_num,header, length(seq), seq } }' "${QUERY_BASE}_Unique_sequences.fa" > "${QUERY_BASE}_Filtered_Unique_sequences.fa"

	# Step 11: Get FASTA sequences for Unique.bed from the reference genome and filter for user-specified length
	bedtools getfasta -fi "$SUB_GENOME" -bed "${SUB_BASE}_Unique.bed" -fo "${SUB_BASE}_Unique_sequences.fa" -fullHeader -s
	awk -v min_len="$MIN_UNIQUE_LENGTH" -v max_len="$MAX_UNIQUE_LENGTH" 'BEGIN { RS=">"; FS="\n"; seq_num = 0 } NR > 1 { header = $1; seq = ""; for (i = 2; i <= NF; i++) seq = seq $i; if (length(seq) >= min_len && (max_len == 0 || length(seq) <= max_len)) { seq_num++; printf ">%d:%s (%d bp)\n%s\n", seq_num,header, length(seq), seq } }' "${SUB_BASE}_Unique_sequences.fa" > "${SUB_BASE}_Filtered_Unique_sequences.fa"
		
	# Step 12: Get FASTA sequences for the query genome from Common_regions.bed (query genome)
	bedtools getfasta -fi "$QUERY_GENOME" -bed "${QUERY_BASE}_Common_regions.bed" -fo "${QUERY_BASE}_Common_regions_sequences_all.fa"
	awk '/^>/{if(!seen[$0]++){header[++s_no]=$0; seq[s_no]=""}next} {seq[s_no]=seq[s_no] $0} END {for(i=1;i<=s_no;i++) printf ">%d:%s (%d bp)\n%s\n", i, substr(header[i], 2), length(seq[i]), seq[i]}' "${QUERY_BASE}_Common_regions_sequences_all.fa" > "${QUERY_BASE}_Common_regions_sequences.fa"


	# Step 13: Get FASTA sequences for Common_regions.bed (reference genome)
	bedtools getfasta -fi "$SUB_GENOME" -bed "${SUB_BASE}_Common_regions.bed" -fo "${SUB_BASE}_Common_regions_sequences_all.fa"  
	awk '/^>/{if(!seen[$0]++){header[++s_no]=$0; seq[s_no]=""}next} {seq[s_no]=seq[s_no] $0} END {for(i=1;i<=s_no;i++) printf ">%d:%s (%d bp)\n%s\n", i, substr(header[i], 2), length(seq[i]), seq[i]}' "${SUB_BASE}_Common_regions_sequences_all.fa" > "${SUB_BASE}_Common_regions_sequences.fa"
	 

	# Step 14: Create folder database.db to store intermediate files
	mkdir -p "$SUB_GENOME"_vs_"$QUERY_GENOME".db

	# Step 15: Move all relevant files to the database directory
	mv "$SUB_GENOME"_vs_"$QUERY_GENOME"_blast_results.txt "$SUB_GENOME"_vs_"$QUERY_GENOME"_common_regions.txt "$SUB_GENOME"_vs_"$QUERY_GENOME"_unique_regions.txt \
	   "${SUB_BASE}.bed" "${QUERY_BASE}.bed" "${SUB_BASE}_Unique.bed" "${QUERY_BASE}_Unique.bed" \
	   "${SUB_BASE}_Common_regions.bed" "${QUERY_BASE}_Common_regions.bed" \
	   "${QUERY_BASE}_unique_regions.bed" "${SUB_BASE}_unique_regions.bed"\
	   "${QUERY_BASE}_Unique_sequences.fa" "${SUB_BASE}_Unique_sequences.fa" \
	   "${SUB_BASE}_Common_regions_sequences_all.fa" "${QUERY_BASE}_Common_regions_sequences_all.fa" \
	   "${SUB_GENOME}.fai" "${QUERY_GENOME}.fai" \
	   "$SUB_GENOME"_vs_"$QUERY_GENOME".db/

	# Notify user about the successful operations
	echo "FASTA sequences extracted for unique regions & common regions."
	echo "Results saved to ${QUERY_BASE}_Filtered_Unique_sequences.fa, ${SUB_BASE}_Filtered_Unique_sequences.fa"

