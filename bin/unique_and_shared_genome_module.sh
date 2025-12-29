#!/bin/bash
set -euo pipefail
trap 'echo "❌ ERROR on line $LINENO"; exit 1' ERR

############################################
# Argument parsing
############################################
core_genome_fasta=""
neighbour_group=()
core_unique_out=""
core_shared_genome=""
unaligned_core_genome_bed=""
OUTPUT_DIR=""
THREADS=1
MIN_SEQ_LENGTH=0
MAX_SEQ_LENGTH=1000000000
MIN_COUNT=0
shared_region=yes

while [[ $# -gt 0 ]]; do
    case "$1" in
        --core_genome_fasta) core_genome_fasta="$2"; shift 2;;
        --neighbour_group)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                neighbour_group+=("$1"); shift
            done;;
        --core_unique_out) core_unique_out="$2"; shift 2;;
        --core_shared_genome) core_shared_genome="$2"; shift 2;;
        --unaligned_core_genome_bed) unaligned_core_genome_bed="$2"; shift 2;;
        --output_dir) OUTPUT_DIR="$2"; shift 2;;
        --t) THREADS="$2"; shift 2;;
        --min_seq_length) MIN_SEQ_LENGTH="$2"; shift 2;;
        --max_seq_length) MAX_SEQ_LENGTH="$2"; shift 2;;
        --min_count) MIN_COUNT="$2"; shift 2;;
        --shared_region) shared_region="$2"; shift 2;;
        *) echo "Unknown argument $1"; exit 1;;
    esac
done

mkdir -p "$OUTPUT_DIR"/{inter_master_bed_FAI,inter_first_blast/genome_bed,inter_conserved_blast/genome_bed}

############################################
# STEP 1: FAI + genome BED
############################################
process_genome() {
    g="$1"
    samtools faidx "$g"
    mv "$g.fai" "$OUTPUT_DIR/inter_master_bed_FAI/"
    awk 'BEGIN{OFS="\t"}{print $1,0,$2}' \
        "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$g").fai" \
        > "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$g").bed"
}
export -f process_genome
export OUTPUT_DIR

parallel -j "$THREADS" process_genome ::: "$core_genome_fasta" "${neighbour_group[@]}"

############################################
# STEP 2: FILTER BLAST (UNIQUE)
############################################
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4,$5,$6,$7,$8,$9}' \
    "$core_unique_out" \
    > "$OUTPUT_DIR/inter_first_blast/filtered_blast_results.tsv"

############################################
# STEP 3: SUBJECT BEDs
############################################
sort -k2,2 -k1,1 -k5,5n "$OUTPUT_DIR/inter_first_blast/filtered_blast_results.tsv" |
awk 'BEGIN{FS=OFS="\t"}
{
    s=($5<$6?$5:$6)-1; if(s<0)s=0
    e=($5>$6?$5:$6)
    print $1,s,e >> "'"$OUTPUT_DIR"'/inter_first_blast/"$2".bed"
}'
############################################
# STEP 4: CONCAT PER GENOME
############################################
concat_bed() {
    g="$1"
    out="$OUTPUT_DIR/inter_first_blast/genome_bed/$(basename "$g").bed"
    fai="$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$g").fai"
    > "$out"
    while read -r c _; do
        [[ -f "$OUTPUT_DIR/inter_first_blast/$c.bed" ]] && \
            cat "$OUTPUT_DIR/inter_first_blast/$c.bed"
    done < "$fai" | bedtools sort -i - | bedtools merge -i - > "$out"
}
export -f concat_bed
export OUTPUT_DIR

parallel -j "$THREADS" concat_bed ::: "${neighbour_group[@]}"

############################################
# STEP 5: UNIQUE REGIONS (num_genomes logic)
############################################
bed_files=("$OUTPUT_DIR/inter_first_blast/genome_bed/"*.bed)
num_genomes=${#bed_files[@]}

if (( num_genomes < 2 )); then
    bedtools subtract \
        -a "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$core_genome_fasta").bed" \
        -b "${bed_files[0]}" \
        > "$OUTPUT_DIR/unique_regions.bed"
else
    bedtools multiinter -i "${bed_files[@]}" |
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' |
    sort -k1,1 -k2,2n | bedtools merge -i - \
        > "$OUTPUT_DIR/core_genome.bed"

    bedtools subtract \
        -a "$OUTPUT_DIR/inter_master_bed_FAI/$(basename "$core_genome_fasta").bed" \
        -b "$OUTPUT_DIR/core_genome.bed" \
        > "$OUTPUT_DIR/unique_regions.bed"
fi

bedtools getfasta \
    -fi "$core_genome_fasta" \
    -bed "$OUTPUT_DIR/unique_regions.bed" \
    -fo "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" \
    -fullHeader

############################################
# STEP 6: UNIQUE FASTA FINAL (RESTORED)
############################################
awk -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" '
BEGIN{OFS=""}
/^>/{if(seq && length(seq)>=min_len && length(seq)<=max_len){
    abs_s=ps+ss; abs_e=ps+se
    print ">"chr":"abs_s"-"abs_e" size="length(seq)" bp\n"seq}
 split($0,a,"[:>-]"); chr=a[2]; ps=a[3]; ss=a[5]; se=a[6]; seq=""; next}
{seq=seq$0}
END{if(seq && length(seq)>=min_len && length(seq)<=max_len){
    abs_s=ps+ss; abs_e=ps+se
    print ">"chr":"abs_s"-"abs_e" size="length(seq)" bp\n"seq}}
' "$OUTPUT_DIR/inter_first_blast/all_extracted_regions.fasta" \
> "$OUTPUT_DIR/unique_genomic_regions.fasta"

############################################
# STEP 7: SHARED / CONSERVED REGIONS
############################################
if [[ "$shared_region" == "yes" && -s "$core_shared_genome" ]]; then
    echo "Processing shared / conserved regions..."

    mkdir -p "$OUTPUT_DIR/inter_conserved_blast/genome_bed"

    # --------------------------------------------------
    # 7.1 Prepare conserved BLAST BEDs
    # --------------------------------------------------
    cp "$core_shared_genome" \
       "$OUTPUT_DIR/inter_conserved_blast/filtered_blast_conserved_results.tsv"

    awk 'BEGIN{FS=OFS="\t"}
    {
        s=($5<$6?$5:$6)-1; if(s<0)s=0
        e=($5>$6?$5:$6)
        print $1,s,e >> "'"$OUTPUT_DIR"'/inter_conserved_blast/"$2".bed"
    }' "$OUTPUT_DIR/inter_conserved_blast/filtered_blast_conserved_results.tsv"

    # --------------------------------------------------
    # 7.2 Concatenate conserved BEDs per genome
    # --------------------------------------------------
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
	parallel -j "$THREADS" process_genome {} "$OUTPUT_DIR" ::: "${neighbour_group[@]}"

    # --------------------------------------------------
    # 7.3 Branch by number of conserved genomes
    # --------------------------------------------------
    bed_files=("$OUTPUT_DIR/inter_conserved_blast/genome_bed/"*_conserved.bed)
    num_bed_files=${#bed_files[@]}

    if (( num_bed_files < 2 )); then
        echo "Only one conserved genome found ($num_bed_files)"

        bedtools getfasta \
            -fi "$core_genome_fasta" \
            -bed "${bed_files[0]}" \
            -name+ \
            > "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta"

    else
        echo "Multiple conserved genomes found ($num_bed_files). Running multiintersect..."

        # --------------------------------------------------
        # 7.4 multiintersect
        # --------------------------------------------------
        bedtools multiinter \
            -i "$OUTPUT_DIR"/inter_conserved_blast/genome_bed/*.bed \
            -header \
            > "$OUTPUT_DIR/inter_conserved_blast/conserved_regions_raw.bed"

        CONSRV_RAW="$OUTPUT_DIR/inter_conserved_blast/conserved_regions_raw.bed"
        FINAL="$OUTPUT_DIR/inter_conserved_blast/max_overlap_conserved_regions.bed"
        IDX_KV_FILE="$OUTPUT_DIR/inter_conserved_blast/idx_kv.txt"

        # --------------------------------------------------
        # 7.5 Build IDX → genome-name map
        # --------------------------------------------------
        read -r hdr < "$CONSRV_RAW"
        IFS=$'\t' read -r -a cols <<< "$hdr"

        list_col=-1
        for i in "${!cols[@]}"; do
            [[ "${cols[i]}" == "list" ]] && { list_col=$i; break; }
        done

        if (( list_col < 0 )); then
            echo "ERROR: 'list' column not found in conserved_regions_raw.bed" >&2
            exit 1
        fi

        > "$IDX_KV_FILE"
        idx=1
        for ((i=list_col+1; i<${#cols[@]}; i++)); do
            nm=$(basename "${cols[i]}" | sed 's/_conserved\.bed//')
            echo "$idx:$nm" >> "$IDX_KV_FILE"
            ((idx++))
        done

        # --------------------------------------------------
        # 7.6 Resolve max-overlap conserved regions
        # --------------------------------------------------
        tail -n +2 "$CONSRV_RAW" \
        | awk -v kv="$IDX_KV_FILE" -v OFS="\t" '
            BEGIN{
                while((getline<kv)>0){
                    split($0,a,":"); id[a[1]]=a[2]
                }
            }
            {
                split($5,ids,",")
                out=""
                for(i in ids){
                    out=(out?out";":"") id[ids[i]]
                }
                print $1,$2,$3,out
            }
        ' \
        | sort -k1,1 -k2,2n \
        > "$FINAL"

        echo "Max-overlap conserved regions written to $FINAL"

        # --------------------------------------------------
        # 7.7 Extract conserved FASTA
        # --------------------------------------------------
        bedtools getfasta \
            -fi "$core_genome_fasta" \
            -bed "$FINAL" \
            -name+ \
            > "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta"
    fi

    # --------------------------------------------------
    # 7.8 Final shared FASTA formatting
    # --------------------------------------------------
    awk -v min_len="$MIN_SEQ_LENGTH" -v max_len="$MAX_SEQ_LENGTH" -v min_count="$MIN_COUNT" '/^>/{if(seq){split(h,p,"::"); sh=p[1]; tg=p[2]; c=split(sh,s,";"); if(length(seq)>=min_len && length(seq)<=max_len && c>=min_count){split(tg,t,"[:\\-]"); chr=t[1]; ps=t[2]; ss=t[4]; se=t[5]; as=ps+ss; ae=ps+se; sz=length(seq); printf "%d\t%d\t>%s:%d-%d | size=%d bp | count=%d | shared: %s\t%s\n",c,sz,chr,as,ae,sz,c,sh,seq}} h=substr($0,2); seq=""; next}{seq=seq$0}END{if(seq){split(h,p,"::"); sh=p[1]; tg=p[2]; c=split(sh,s,";"); if(length(seq)>=min_len && length(seq)<=max_len && c>=min_count){split(tg,t,"[:\\-]"); chr=t[1]; ps=t[2]; ss=t[4]; se=t[5]; as=ps+ss; ae=ps+se; sz=length(seq); printf "%d\t%d\t>%s:%d-%d | size=%d bp | count=%d | shared: %s\t%s\n",c,sz,chr,as,ae,sz,c,sh,seq}}}' "$OUTPUT_DIR/inter_conserved_blast/all_extracted_conserved_regions.fasta" | sort -t $'\t' -k1,1nr -k2,2nr | awk -F'\t' '{print $3"\n"$4}' > "$OUTPUT_DIR/shared_genomic_regions.fasta"


else
    echo "Shared region analysis skipped (shared_region=$shared_region or file empty)"
fi

echo "Pipeline completed successfully. Output saved in $OUTPUT_DIR"
exit 0
